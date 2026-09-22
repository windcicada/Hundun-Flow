# HUNDUN-FLOW

English | [简体中文](README.zh-CN.md)

HUNDUN-FLOW is an open-source simulation code for turbulent reacting flows and dilute sprays. It provides a computational foundation for turning physical models into coupled numerical methods and developing research software. Its C++17 implementation combines Message Passing Interface (MPI) parallelism with explicit module interfaces and coding-agent guidance. The Apache License 2.0 supports research, commercial use, and proprietary extensions under its terms.

The solver uses finite-volume discretization on Cartesian grids, local-reconstruction immersed boundary methods (IBM), and pressure-based coupling. Large-eddy simulation (LES), finite-rate chemistry, partially stirred reactor (PaSR) closure, Eulerian stochastic fields (ESF), and Lagrangian parcels share thermodynamic services and conservation accounting. This organization supports AI-assisted model implementation while keeping equations, state ownership, and numerical checks explicit.

The name comes from the description of Dijiang in the “Western Mountains” chapter of *Shan Hai Jing* (*Classic of Mountains and Seas*). Its pouch-shaped, fiery-red appearance provides an image for a bounded combustion chamber containing interacting flow and chemical processes. *Hundun* also denotes an undifferentiated whole, from which modeling seeks to establish physical understanding.

[Version 1.1.0](https://github.com/windcicada/Hundun-Flow/releases/tag/v1.1.0) covers fixed Cartesian grids, static immersed boundaries, and dilute parcels. Release checks cover selected flow, transport, chemistry, dynamic micromixing, gas–liquid exchange, and restart combinations. Long-duration combustor calculations and statistical assessment remain separate validation tasks; the recorded scope is given in [release acceptance](docs/accept.md).

## Software architecture

The implementation has a flat source layout with modules grouped by responsibility. Public headers reside in `versions/v0.4/include/hundun`, and implementations reside in `versions/v0.4/src`. The directory name retains the implementation lineage; the root `VERSION` file defines the product version.

| Layer / source prefix | Responsibility | Interface boundary |
| --- | --- | --- |
| Application: `app_*` | Parse case input, compile execution plans, expose the command-line interface, and record build identity. | Case preparation resolves model choices and resource requirements before time integration. |
| State and orchestration: `core_*` | Own fields, workspaces, time histories, the product driver, and coupled step acceptance. | Trial states become accepted states after the configured checks succeed across MPI ranks. |
| Physical models: `models_*` | Implement combustion closures, ESF, micromixing, chemistry adapters, and parcel kernels. | Value inputs, borrowed views, and providers return candidate changes and reports. |
| Material services: `physics_*` | Provide thermodynamics, transport properties, LES closures, and source-contribution interfaces. | State revisions and material identities determine which properties and sources are valid. |
| Discrete equations: `solver_*` | Assemble momentum, pressure, enthalpy, and scalar equations; apply linear solvers and preconditioners. | Operators consume mesh metrics, boundary states, and consistent face fluxes. |
| Mesh and boundaries: `mesh_*`, `bc_*` | Build Cartesian meshes, query geometry, prepare IBM reconstruction, and impose physical boundary conditions. | Geometry and reconstruction plans supply the data used by the discrete equations. |
| Communication: `parallel_*` | Manage CPU resources, halo exchange, and distributed donor communication. | Prepared communication plans define ownership, buffers, and completion requirements. |
| Persistence and output: `io_*` | Write visualization, checkpoints, monitoring, and run evidence. | Output refers to accepted states; restart validates state, assets, and method history. |

## Numerical methods

### Mesh, transport, and pressure coupling

The mesh is fixed during a run and may use uniform or stretched Cartesian spacing. STL surfaces or imported Cartesian markers define static immersed boundaries. Local reconstruction supplies boundary-adjacent states; adaptive reconstruction order refers to the interpolation order on the fixed mesh.

| Component | Method and role |
| --- | --- |
| Spatial discretization | Cell-centered finite volumes with Rhie–Chow face-flux construction couple pressure and velocity. IBM reconstruction supplies wall-adjacent values. |
| Time integration | The default CN/BE scheme combines Crank–Nicolson (CN) momentum integration with backward Euler (BE) reactive-scalar updates. Ordinary single-fluid enthalpy uses CN; reacting and spray configurations use BE for enthalpy and composition. |
| Pressure coupling | CN/BE uses `outer_corrected`. Explicit `backward_euler` configurations use PISO or SIMPLE scheduling. Coupling checks evaluate the original discrete equations. |
| LES | Smagorinsky, WALE, and Vreman closures provide subgrid-scale (SGS) viscosity. The supplied reacting-flow examples use Vreman. |
| Species and enthalpy | Independent species fields reconstruct the full composition. Total thermochemical enthalpy includes species formation enthalpies; material queries recover temperature and transport properties. |
| Linear solution | Iterative solvers use configured iteration budgets and report independently evaluated residuals. Reusable transport coefficients and preconditioners reduce repeated assembly work. |

### Combustion and micromixing

The turbulence–chemistry closure and chemistry representation are separate choices. Closure identifiers are `finite_rate_mean`, `pasr_algebraic_v1`, and `esf_tpdf`. The `direct_cantera` representation supplies the reacting backend; `analytic_isomer` provides an analytic chemistry option for controlled tests.

| Component | Method and coupling |
| --- | --- |
| Mean-state finite rate | Integrate chemistry at the mean thermochemical state and register its species source. This path supplies the mean-state baseline. |
| Algebraic PaSR | Combine named chemical and mixing timescales to obtain a reaction fraction. The same fraction weights all reacting-species contributions and their associated heat-release accounting. Unavailable timescales produce an explicit status. |
| Transported PDF / ESF | Represent the transported probability density function (TPDF) with an ensemble of stochastic composition fields. Each spatial MPI rank holds all local fields. Production checks cover 2, 4, 8, and 16 fields, subject to species and memory limits. |
| Micromixing | Interaction by Exchange with the Mean (IEM) supplies the baseline mixing operator. The Turbulence–Chemistry Recursive (TCR) model recovers mixing-state parameters from reaction-rate information and supplies history-dependent mixing controls. |
| Dynamic TCR | `cdphyso_dynamic_v1` and `dyn711_v1` retain distinct rate statistics and update schedules. `experimental` applies accepted controls; `shadow` records dynamic history while retaining baseline IEM mixing. |
| Chemistry services | Shared thermodynamic and kinetic interfaces support state recovery, rate queries, and interval integration. Prepared batch workspaces reuse backend resources. Mechanism identity includes phase, species ordering, and configured rate regularization. |

PaSR uses the scalar-dissipation mixing timescale

$$
\tau_{\mathrm{mix}}=\frac{C_Z\Delta^2}{2(D+D_t)},
\qquad D_t=\frac{\nu_t}{Sc_t},
\qquad \kappa_{\mathrm{PaSR}}=\frac{\tau_{\mathrm{chem}}}{\tau_{\mathrm{chem}}+\tau_{\mathrm{mix}}}.
$$

Here, $C_Z$ is an explicit model coefficient, $\Delta$ is the filter width, $D$ is molecular diffusivity, $\nu_t$ is SGS kinematic viscosity, and $Sc_t$ is the turbulent Schmidt number. The chemical-timescale definition belongs to the selected model contract. The PaSR reaction fraction $\kappa_{\mathrm{PaSR}}$ and the TCR mixing-state parameter $\kappa$ have separate definitions and roles.

For the current CN/BE ESF schedule, transport, paired Wiener increments, implicit IEM mixing, and chemistry precede the flow corrections. Outer corrections reuse the resulting chemistry state. The physical ensemble mean carries composition and enthalpy; an independent `field0` supplies the reduced-noise state used for pressure coupling. TCR statistics follow the selected model's schedule, and accepted histories participate in the common step decision. [ESF configuration](docs/esf.md) and [TCR definitions](docs/tcr.md) specify these contracts.

### Lagrangian sprays

A parcel represents a weighted population of single-component droplets. Stable IDs identify injection, breakup, migration, and restart records. Gas sampling and source deposition connect parcel motion to the Eulerian fields.

| Process | Method and accounting |
| --- | --- |
| Injection and trajectories | Injection candidates retain deterministic identities. Adaptive parcel substeps resolve cell crossings, wall impacts, exits, and evaporation endpoints. |
| Momentum and transfer | Schiller–Naumann drag uses local gas and liquid properties. Heat and mass transfer follow the selected evaporation formulation. Wall rebound uses the static-geometry intersection and surface normal. |
| Evaporation | `abramzon_sirignano` and `thick_exchange` select Abramzon–Sirignano and THICK_EX formulations. Liquid properties, vapor-species mapping, and absolute liquid enthalpy share explicit material references. |
| Breakup | `tab` selects the Taylor analogy breakup (TAB) model; `stochastic_sgs` selects the SGS-exposure model. Parent replacement, child IDs, and model histories enter candidate state. |
| Two-way exchange | Parcel mass, momentum, enthalpy, and kinetic-energy changes define gas-side exchange. Each physical exchange is counted once; every ESF receives the same gas source. |
| Migration and persistence | Endpoint ownership determines MPI migration. Parcel state and histories are accepted or rejected together with the gas state and retained in native restart records. |

Liquid-property conventions and evaporation details are given in [liquid properties](docs/liquid.md).

## Design rationale

| Design choice | Numerical or software basis |
| --- | --- |
| Shared conservation accounting | Consistent species, element, momentum, and energy definitions make chemical and interphase sources traceable. Raw budget defects and discrete-equation residuals remain separately visible. |
| Common step acceptance and rollback | Flow, stochastic fields, model clocks, and parcels represent the same accepted time endpoint. A rejected attempt retains the previous accepted state for retry. |
| Explicit module contracts | State revisions, units, material identities, borrowed views, and workspace capacities define each call. This supports isolated model tests and bounded code changes. |
| Deterministic stochastic identities | Counter-based random numbers, stable parcel IDs, and accepted-step clocks preserve retry identities and support cross-partition restart checks. |
| Prepared execution | HUNDUN-owned workspaces and communication buffers are prepared before use. Shared transport assembly, preconditioner reuse, and batched chemistry limit repeated work. Third-party resource behavior remains backend-specific. |
| Agent-assisted development | Named interfaces, case checks, contribution rules, and run evidence give coding agents the same verifiable numerical contracts used in manual development. |

## Build and run

### Packaged runtime

The [v1.1.0 release](https://github.com/windcicada/Hundun-Flow/releases/tag/v1.1.0) provides source, Linux x86-64 runtime, and example archives, together with acceptance notes, dependency licenses, and `SHA256SUMS`. The runtime contains its own MPI and chemistry dependencies. From the extracted runtime directory:

```sh
./run --version
./run check case
./mpirun -n 2 ./run run case --restart restart --output next --steps 1
```

### Source build

The documented release profile uses CMake 3.21+, Ninja, Clang 15, lld, libstdc++ ABI1, MPI 3+, and FP64 with ThinLTO. `HUNDUN_CANTERA_ROOT` points to the project's pinned, hash-verified Cantera 3.2.0 SDK. Its Linux x86-64 build targets glibc 2.35 or later. Select the intended `clang` and `clang++` through `PATH` before configuring.

Run from the repository root:

```sh
HUNDUN_CANTERA_ROOT=/path/to/cantera cmake --preset release
cmake --build --preset release -j 2
export PATH="$PWD/b3/versions/v0.4:$PATH"
hundun --version
```

### New calculation and restart

Generate a starter case, inspect its compiled plan, and run a short calculation:

```sh
hundun init-case --output case
mpiexec -n 4 hundun check case --dry-plan
mpiexec -n 4 hundun run case --output run --steps 10 \
  --output-interval 0 --restart-interval 10 --diagnostics-interval 1
```

`hundun check` reports the time scheme, pressure coupling, model selection, and resource plan. The default Courant–Friedrichs–Lewy (CFL) target is 0.30, with a holding interval of 0.25–0.35; explicit case settings take precedence. Case checks verify configuration and resource admission; numerical validation uses the observables and tolerances defined for each case.

Continue from the saved state, using a new output directory:

```sh
mpiexec -n 2 hundun run case --restart run/Restart --output next \
  --steps 10 --output-interval 0 --restart-interval 10 --diagnostics-interval 1
```

Native restart retains fields, time histories, face fluxes, stochastic fields, model clocks, and parcels. It supports the documented repartitioning path when the MPI rank count changes. The case and physical assets retain their identities. Detailed instructions are in [installation, running, and restart](docs/run.md).

| Run output | Purpose |
| --- | --- |
| `monitor.jsonl` | Time-step controls, equation iterations, and module timings. |
| `diagnostics.jsonl` | Mass, species, element, and energy accounting. |
| `evidence.jsonl` | Program, input, and accepted-state identities. |
| `Restart/` | Verified checkpoint generations for continuation. |

## Model examples

| Example | Configuration |
| --- | --- |
| [GTMC gas model](examples/g/README.md) | Methane/JL4 chemistry, four ESF, Vreman, and dynamic TCR. |
| [624CF gas–liquid model](examples/cf/README.md) | Four-step kerosene chemistry, two ESF, Vreman, dynamic TCR, static IBM, THICK_EX evaporation, SGS breakup, and fixed thermodynamic pressure. |

Both directories contain small native model demonstrations and fresh-run/restart commands. Full combustor geometries and long-duration reference calculations have their own inputs and acceptance records.

## GTMC combustion-field demonstration

The following images show a native GTMC restart calculation at accepted step 32016 and $t=21.815583$ ms on a $153\times325\times151$ mesh. The snapshot uses methane/JL4 chemistry, four-field ESF/TPDF, IEM, Vreman, and static IBM. It illustrates an instantaneous field; statistical assessment and release acceptance have separate evidence records. [Image provenance](docs/hot.md) records the source state and rendering checks.

The three-dimensional view combines translucent combustor geometry, 1400 K and 1800 K temperature isosurfaces, and instantaneous streamlines colored by speed.

![GTMC geometry, temperature isosurfaces, and speed-colored streamlines](docs/images/gtmc-3d.png)

The central plane at $z=-0.075$ mm shows speed and temperature recovered from the physical mean enthalpy and composition. Gray regions mark IBM solids. The plane retains the original mesh resolution; the three-dimensional view samples every second cell along each axis.

![GTMC central-plane speed and temperature](docs/images/gtmc-mid.png)

## Documentation and contributions

| Topic | Reference |
| --- | --- |
| Installation, case checks, execution, and restart | [Run guide](docs/run.md) |
| Release contents and verification scope | [Release notes](docs/release.md), [acceptance record](docs/accept.md) |
| Chemistry and material definitions | [Reaction-rate regularization](docs/chem.md), [liquid properties](docs/liquid.md) |
| Stochastic fields and micromixing | [ESF](docs/esf.md), [dynamic TCR](docs/tcr.md) |
| Development and contribution provenance | [Repository guidance](AGENTS.md), [contributing](CONTRIBUTING.md), [Developer Certificate of Origin](DCO.md) |

The linked model and development documents are maintained in the repository, primarily in Chinese. Contributions require documented provenance and a DCO sign-off (`git commit -s`).

## License

HUNDUN-FLOW is distributed under the [Apache License 2.0](LICENSE). Derivative programs may use open-source or proprietary licensing, subject to the applicable license conditions. Third-party dependencies retain their own licenses; see [THIRD_PARTY.md](THIRD_PARTY.md).
