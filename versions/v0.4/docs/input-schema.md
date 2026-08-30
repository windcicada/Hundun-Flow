<!-- SPDX-License-Identifier: Apache-2.0 -->

# HUNDUN-FLOW v0.4 case-input schema

Status: Task 20 production controls and wire-v9 schema, 2026-08-21

## Authority and directory layout

`case.json` is the only configuration authority: it selects every input and no
referenced file may add or override JSON structure. Referenced `.d` files carry
typed table data and STL carries geometry bytes. Every referenced file must be a
unique, regular direct child of the case root. Absolute paths, parent traversal,
nested paths, incorrect suffixes, missing files, symbolic links, and hard-link
aliases are rejected.

Rank 0 alone opens the case-root descriptor, `case.json`, `.d`, and STL files. It
reads and parses the thermophysical `.d` file exactly once, normalizes the closed
schema, hashes the complete typed model and referenced bytes, and broadcasts a
bounded wire-v9 model. The runtime model retains compact typed thermophysical
data, never the source text. Other ranks never inspect the filesystem or parse
JSON or `.d` text. Failed compilation leaves the caller's output model unchanged
on every rank.

## Complete closed schema

Every object is closed. Unknown, duplicate, missing, or incorrectly typed keys
are errors. The root has exactly these required keys:

`schema_version`, `units`, `mesh`, `flow`, `solver`, `thermophysics`,
`transported_scalars`, `boundaries`, `schemes`, and `time`.

`turbulence` is the only optional key anywhere in the root schema. Omitting it
selects `vreman_wall_function`. JSON `schema_version` remains `1`; wire version
`9` is an internal typed-model format.

`CaseCompiler` closes and types the JSON, checks local enum, finite-number, and
range constraints, canonicalizes the mesh, and publishes the typed model. It
does not decide compatibility between different faces or between a boundary
kind and its parameters. Those cross-field and cross-face rules belong to the
later `BoundaryCompiler` plan-compilation stage.

This is a complete valid open, low-Mach case. All six boundary faces and every
required boundary, scheme, and time-control field are explicit.

```json
{
  "schema_version": 1,
  "units": "SI",
  "mesh": {
    "kind": "uniform",
    "domain": {
      "lower": [0.0, 0.0, 0.0],
      "upper": [1.0, 0.5, 0.5]
    },
    "exact_cells": [8, 4, 4],
    "base_spacing": null,
    "minimum_spacing": [0.125, 0.125, 0.125],
    "max_growth_ratio": 1.0,
    "focus_regions": [],
    "limits": {
      "max_global_cells": 4096,
      "max_memory_bytes_per_rank": 67108864
    },
    "data_files": [],
    "immersed_boundary": null
  },
  "flow": {
    "model": "single_phase_low_mach_compressible",
    "pressure_reference": "boundary_absolute",
    "reacting": false
  },
  "solver": {
    "coupling": "PISO",
    "pressure_correctors": 2,
    "pressure_linear": {
      "absolute_tolerance": 1e-13,
      "relative_tolerance": 1e-13,
      "maximum_iterations": 400,
      "true_residual_interval": 4,
      "krylov_restart": 12
    },
    "terminal_tolerances": {
      "eos": 1e-10,
      "continuity": 1e-10,
      "closed_mass": 1e-10,
      "gauge": 1e-10
    }
  },
  "turbulence": {
    "model": "none"
  },
  "thermophysics": {
    "data_file": "thermophysics.d"
  },
  "transported_scalars": [],
  "boundaries": {
    "x_min": {
      "flow_kind": "velocity_inlet",
      "thermal_kind": "none",
      "velocity": [0.1, 0.0, 0.0],
      "direction": [1.0, 0.0, 0.0],
      "backflow_velocity": [0.0, 0.0, 0.0],
      "mass_flow_rate": 0.0,
      "pressure": 0.0,
      "temperature": 300.0,
      "total_pressure": 0.0,
      "total_temperature": 0.0,
      "backflow_temperature": 0.0,
      "heat_flux": 0.0,
      "relaxation": 0.0,
      "mach_limit": 0.95,
      "allow_backflow": false,
      "scalars": []
    },
    "x_max": {
      "flow_kind": "pressure_outlet",
      "thermal_kind": "none",
      "velocity": [0.0, 0.0, 0.0],
      "direction": [1.0, 0.0, 0.0],
      "backflow_velocity": [0.0, 0.0, 0.0],
      "mass_flow_rate": 0.0,
      "pressure": 101325.0,
      "temperature": 0.0,
      "total_pressure": 0.0,
      "total_temperature": 0.0,
      "backflow_temperature": 300.0,
      "heat_flux": 0.0,
      "relaxation": 0.0,
      "mach_limit": 0.95,
      "allow_backflow": false,
      "scalars": []
    },
    "y_min": {
      "flow_kind": "symmetry",
      "thermal_kind": "none",
      "velocity": [0.0, 0.0, 0.0],
      "direction": [0.0, -1.0, 0.0],
      "backflow_velocity": [0.0, 0.0, 0.0],
      "mass_flow_rate": 0.0,
      "pressure": 0.0,
      "temperature": 0.0,
      "total_pressure": 0.0,
      "total_temperature": 0.0,
      "backflow_temperature": 0.0,
      "heat_flux": 0.0,
      "relaxation": 0.0,
      "mach_limit": 0.95,
      "allow_backflow": false,
      "scalars": []
    },
    "y_max": {
      "flow_kind": "symmetry",
      "thermal_kind": "none",
      "velocity": [0.0, 0.0, 0.0],
      "direction": [0.0, 1.0, 0.0],
      "backflow_velocity": [0.0, 0.0, 0.0],
      "mass_flow_rate": 0.0,
      "pressure": 0.0,
      "temperature": 0.0,
      "total_pressure": 0.0,
      "total_temperature": 0.0,
      "backflow_temperature": 0.0,
      "heat_flux": 0.0,
      "relaxation": 0.0,
      "mach_limit": 0.95,
      "allow_backflow": false,
      "scalars": []
    },
    "z_min": {
      "flow_kind": "symmetry",
      "thermal_kind": "none",
      "velocity": [0.0, 0.0, 0.0],
      "direction": [0.0, 0.0, -1.0],
      "backflow_velocity": [0.0, 0.0, 0.0],
      "mass_flow_rate": 0.0,
      "pressure": 0.0,
      "temperature": 0.0,
      "total_pressure": 0.0,
      "total_temperature": 0.0,
      "backflow_temperature": 0.0,
      "heat_flux": 0.0,
      "relaxation": 0.0,
      "mach_limit": 0.95,
      "allow_backflow": false,
      "scalars": []
    },
    "z_max": {
      "flow_kind": "symmetry",
      "thermal_kind": "none",
      "velocity": [0.0, 0.0, 0.0],
      "direction": [0.0, 0.0, 1.0],
      "backflow_velocity": [0.0, 0.0, 0.0],
      "mass_flow_rate": 0.0,
      "pressure": 0.0,
      "temperature": 0.0,
      "total_pressure": 0.0,
      "total_temperature": 0.0,
      "backflow_temperature": 0.0,
      "heat_flux": 0.0,
      "relaxation": 0.0,
      "mach_limit": 0.95,
      "allow_backflow": false,
      "scalars": []
    }
  },
  "schemes": {
    "momentum": "limited_central2",
    "enthalpy": "limited_central2",
    "species": "tvd2",
    "passive_scalar": "tvd2",
    "diffusion": "central2",
    "limiter": 1.0
  },
  "time": {
    "control": "adaptive_flow",
    "scheme": "variable_bdf2",
    "initial_dt": 0.0001,
    "minimum_dt": 1e-10,
    "maximum_dt": 0.01,
    "convective_cfl": 0.8,
    "viscous_cfl": 0.5,
    "thermal_cfl": 0.5,
    "species_cfl": 0.5,
    "acoustic_cfl": 0.8,
    "maximum_growth": 1.25,
    "retry_factor": 0.5,
    "maximum_retries": 8,
    "minimum_bdf_ratio": 0.2,
    "maximum_bdf_ratio": 5.0
  }
}
```

## Mesh

The `mesh` object has exactly `kind`, `domain`, `exact_cells`,
`base_spacing`, `minimum_spacing`, `max_growth_ratio`, `focus_regions`,
`limits`, `data_files`, and `immersed_boundary`. `domain` has exactly `lower` and
`upper`; `limits` has exactly `max_global_cells` and
`max_memory_bytes_per_rank`; each focus-region object has exactly `lower`,
`upper`, and `target_spacing`.

All mesh numbers must be finite IEEE-754 values. Each domain upper coordinate
must be greater than its lower coordinate. Every spacing is positive,
`max_growth_ratio >= 1`, every cell count is a positive signed-32-bit
representable integer, and both limits are positive. When `exact_cells` is
present, its checked product cannot exceed `max_global_cells`.

### Uniform Cartesian mesh

For `kind: "uniform"`:

- `exact_cells` is a required three-integer vector;
- `base_spacing` is `null`;
- `focus_regions` is empty;
- `max_growth_ratio` is exactly `1.0`; and
- each `minimum_spacing` component is no greater than the corresponding
  domain span divided by the exact cell count.

### Tensor-stretched Cartesian mesh

For `kind: "tensor_stretched"`, `base_spacing` is a required positive vector
and each component is at least `minimum_spacing`. `exact_cells` may be `null`
or a three-integer vector; when present it fixes the global counts while base
spacing and focus regions define the target distribution. Target spacing must
be component-wise within `[minimum_spacing, base_spacing]`.

A focus box must overlap the domain with positive extent in all three axes.
Partly external boxes are clipped to the domain; wholly external or merely
touching boxes are rejected. After clipping, regions are sorted
lexicographically by lower bound, upper bound, and target spacing, then exact
duplicates are removed. Negative zero is normalized to positive zero.

## Flow, solver, and turbulence

These objects are closed and have the following exact fields and values:

| Object | Exact fields | Accepted values |
| --- | --- | --- |
| `flow` | `model`, `pressure_reference`, `reacting` | `model` is `single_phase_low_mach_compressible`; `pressure_reference` is `boundary_absolute` or `closed_mass`; `reacting` is `false` |
| `solver` | `coupling`, `pressure_correctors`, `pressure_linear`, `terminal_tolerances` | `coupling` is `PISO`; `pressure_correctors` is unsigned integer `2`; the two nested control objects are described below |
| `turbulence` | `model` | Optional object; `model` is `vreman_wall_function`, `wale`, or `none` |

`pressure_linear` is closed and has the exact fields
`absolute_tolerance`, `relative_tolerance`, `maximum_iterations`,
`true_residual_interval`, and `krylov_restart`. Both tolerances are finite and
strictly between zero and one. `maximum_iterations` is in `[1, 1000000]`, the
true-residual interval is in `[1, maximum_iterations]`, and the fixed-workspace
FGMRES restart is in `[2, 64]`. Residual convergence is always decided from a
recomputed true residual at the configured interval; the recursive residual is
not an acceptance authority.

`terminal_tolerances` is closed and has the exact finite positive fields `eos`,
`continuity`, `closed_mass`, and `gauge`, each less than one. The pressure
relative tolerance cannot be looser than either the continuity or gauge gate.
These values are compiled into the immutable product plan and its semantic
fingerprint, so changing the numerical work budget changes product identity.

For migration only, the compiler also accepts the former two-field `solver`
object and assigns the strict values shown in the example. Newly generated,
validated, and performance-candidate cases must use the explicit four-field
form; there is no runtime environment or command-line tolerance override.

`pressure_closure` is not a schema key. `CaseCompiler` accepts either listed
`pressure_reference` enum without deriving pressure authority from the faces.
During boundary-plan compilation, `boundary_absolute` must have at least one
actual absolute-pressure authority among `pressure_outlet`, `nscbc_outlet`,
`static_state_inlet`, and `total_state_inlet`; `closed_mass` must have none. A
subsonic `nscbc_inlet` preserves the outgoing acoustic characteristic and is
therefore not, by itself, an absolute-pressure authority. The later
thermodynamics plan supplies the global mass closure for a valid closed-mass
boundary plan.

## Thermophysics

`thermophysics` is required and closed, with exactly one key:

```json
"thermophysics": {"data_file": "thermophysics.d"}
```

`data_file` is a unique direct-root filename ending in `.d`. Rank 0 reads that
file once and parses the strict `HUNDUN_THERMOPHYSICS_V1` token format. Whitespace
separates tokens and `#` begins a line comment. A minimal one-species file is:

```text
HUNDUN_THERMOPHYSICS_V1
temperature_bounds 200 3000
temperature_inversion 1e-12 32
closed_mass_newton 1e-12 24 0.25
species_count 1
species air
molecular_weight 28.96546
temperature_switch 1000
nasa7_low 3.5 0 0 0 0 0 0
nasa7_high 3.5 0 0 0 0 0 0
transport_sutherland 1.716e-5 273.15 110.4 0.71
end_species
end
```

The three global records respectively define the inclusive temperature interval,
bounded enthalpy-to-temperature inversion tolerance/iteration count, and closed
mass pressure-correction tolerance/iteration count/maximum relative step. The
file declares 1 to 65 uniquely named species. Each species supplies molecular
weight, an interior NASA-7 switch temperature, low/high seven-coefficient
records, and exactly one transport record:

- `transport_constant <mu> <k>`; or
- `transport_sutherland <mu_ref> <T_ref> <S> <Pr>`.

All numbers are finite and the physical scales and tolerances are positive.
Iteration counts are bounded by 256. Validation proves `cp/R > 1` throughout
both coefficient intervals and requires continuous `cp` at the switch. An
input `h` jump within the strict acceptance tolerance is removed during cold
compilation by replacing only the high-interval NASA integration constant;
larger jumps are rejected. The canonical typed model is then shared by the
thermodynamics and transport plans. The same validator guards root
serialization and wire-v9 decoding, so a
structurally or scientifically invalid typed thermophysical payload is never
published. Only the compact typed records and source filename survive in
`ValidatedModel`; source text is discarded after compilation.

## Transported-scalar catalog

`transported_scalars` is a required closed array and may be empty. It is the
single canonical declaration-order catalog for inert species and passive
scalars. A catalog entry is a closed object with exactly `stable_name`,
`role`, `molecular_schmidt`, and `turbulent_schmidt`:

```json
"transported_scalars": [
  {"stable_name": "O2", "role": "species",
   "molecular_schmidt": 0.7, "turbulent_schmidt": 0.9},
  {"stable_name": "mixture_fraction", "role": "passive_scalar",
   "molecular_schmidt": 0.7, "turbulent_schmidt": 0.9}
]
```

`role` is exactly `species` or `passive_scalar`. A `stable_name` is 1 to 255
bytes and contains only letters, digits, `_`, or `-`. Names are globally
unique within the catalog; `U`, `pi`, and `h` are reserved state-field names
and cannot be declared. At most 64 transported scalars are accepted. The
two Schmidt numbers are finite, dimensionless, and strictly positive. The
runtime forms molecular and SGS diffusivity from their corresponding dynamic
viscosity divided by these compiled constants. The catalog, including
declaration order, role, and both closures, is retained in the typed model,
wire payload, and semantic fingerprint.

Boundary scalar arrays are not a second field catalog. Boundary-plan
compilation requires every nonperiodic physical face to contain exactly one
closure for every catalog entry. Periodic faces contain no scalar entries and
use the halo topology. The compiler maps catalog names and roles to stable
field IDs and static species/passive-scalar dispatch before hot-loop execution.

## Boundaries

`boundaries` has exactly six keys: `x_min`, `x_max`, `y_min`, `y_max`,
`z_min`, and `z_max`. Each value is a closed boundary-face object with exactly
these fields:

`flow_kind`, `thermal_kind`, `velocity`, `direction`, `backflow_velocity`,
`mass_flow_rate`, `pressure`, `temperature`, `total_pressure`,
`total_temperature`, `backflow_temperature`, `heat_flux`, `relaxation`,
`mach_limit`, `allow_backflow`, and `scalars`.

The three velocity/direction fields are finite three-number vectors. The nine
scalar parameters from `mass_flow_rate` through `mach_limit` are finite JSON
numbers, `allow_backflow` is a Boolean, `relaxation >= 0`, and
`0 < mach_limit < 1`.

`heat_flux` is the outward heat flux in W/m2: positive removes energy from
the fluid.  At a `heat_flux_wall`, the runtime converts it to the temperature
normal gradient using the local conductivity before applying the boundary
stencil; it is not interpreted as an enthalpy gradient.

Both `flow_kind` and `thermal_kind` use the closed `BoundaryKind` enum:

- `none`, `velocity_inlet`, `mass_flow_inlet`, `static_state_inlet`,
  `total_state_inlet`, `pressure_outlet`, `nscbc_inlet`, `nscbc_outlet`;
- `no_slip_wall`, `moving_wall`, `slip`, `symmetry`, or `periodic`; and
- `adiabatic_wall`, `isothermal_wall`, or `heat_flux_wall`.

The case compiler checks only that each spelling is in this enum and that the
face's numeric fields satisfy the local finite/range constraints above. In
particular, it does not reject `flow_kind: "none"`, a nonthermal
`thermal_kind`, or a kind paired with inactive/insufficient parameters.

`BoundaryCompiler` subsequently applies boundary-plan compatibility. At that
stage, `flow_kind: "none"` is invalid; `thermal_kind` must be `none`,
`adiabatic_wall`, `isothermal_wall`, or `heat_flux_wall`; and a non-`none`
thermal condition requires a wall-like flow kind. An isothermal wall requires
`temperature > 0`. A velocity inlet also requires `temperature > 0` and
compiles a face-cell-local enthalpy target resolver, so the thermodynamics
service can map its declared temperature and species composition to `h`. A
pressure or NSCBC outlet requires `pressure > 0`; a
mass-flow inlet requires `mass_flow_rate > 0`; a static-state or NSCBC inlet
requires positive pressure and temperature; and a total-state inlet requires
positive total pressure and total temperature.

Also during `BoundaryCompiler`, opposite periodic faces must be paired. A
periodic face requires `thermal_kind: "none"` and an empty `scalars` array.
`no_slip_wall` requires an exactly zero velocity; nonzero wall velocity uses
`moving_wall`. Outlet backflow may be enabled only with a positive declared
backflow temperature and an inward-pointing backflow velocity. Every
transported scalar on such a face must also declare a `dirichlet` backflow
closure. Independent species values lie in `[0, 1]` and sum to at most one;
passive-scalar backflow values need only be finite.

`velocity_inlet`, `mass_flow_inlet`, `static_state_inlet`,
`total_state_inlet`, and `nscbc_inlet` own the entering composition. On each
such face every catalog entry whose role is `species` must use `dirichlet`;
each independent species fraction lies in `[0, 1]`, and their sum is at most
one. On walls and ordinary outlets a species may instead use a flux,
convective, or zero-gradient closure. There, only species that actually use
`dirichlet` participate in the same bounds and sum check; inactive numeric
values on non-Dirichlet closures are not mixed into a composition. Passive
scalars remain finite but otherwise unbounded.

Each entry in `scalars` is a closed object with exactly `stable_name`, `kind`,
`value`, `backflow_kind`, and `backflow_value`. `stable_name` is a nonempty,
per-face-unique identifier of at most 255 bytes containing only letters,
digits, `_`, or `-`. Both kinds use `dirichlet`, `normal_flux`,
`zero_gradient`, or `convective`; both values are finite. Outside a
backflow-enabled outlet, `backflow_kind: "zero_gradient"` and
`backflow_value: 0` are the explicit inactive defaults. A face has at most 64
scalar entries. Empty scalar arrays remain explicit when the case declares no
species or passive-scalar boundary data.

For `normal_flux`, `value` is the outward diffusive flux of the conserved
scalar. Positive values leave the fluid. The boundary resolver converts it to
`normal_gradient = -value / effective_diffusivity` after the local
diffusivity is available.

## Schemes

`schemes` is a closed object with exactly `momentum`, `enthalpy`, `species`,
`passive_scalar`, `diffusion`, and `limiter`.

The four convection fields accept `central2`, `limited_central2`, or `tvd2`.
`diffusion` accepts only `central2`. `limiter` is finite and lies in `[0, 1]`.
The selected convection schemes determine the compiled ghost width; no scheme
name is inspected in the hot boundary or solver path.

## Time control

`time` is a closed object with every `TimeControlSpec` field, exactly:

`control`, `scheme`, `initial_dt`, `minimum_dt`, `maximum_dt`,
`convective_cfl`, `viscous_cfl`, `thermal_cfl`, `species_cfl`, `acoustic_cfl`,
`maximum_growth`, `retry_factor`, `maximum_retries`, `minimum_bdf_ratio`, and
`maximum_bdf_ratio`.

`control` is `fixed`, `adaptive_flow`, or `adaptive_acoustic`; `scheme` is
`backward_euler` or `variable_bdf2`. All real-valued fields are finite. The
following bounds apply:

- `0 < minimum_dt <= initial_dt <= maximum_dt`;
- all five CFL values are positive;
- `maximum_growth >= 1`;
- `0 < retry_factor < 1`;
- `maximum_retries` is a positive unsigned 32-bit integer; and
- `0 < minimum_bdf_ratio <= maximum_bdf_ratio`.

`adaptive_flow` applies convective, viscous, thermal, and species hard limits
but deliberately has no acoustic hard limit. `adaptive_acoustic` also applies
the acoustic limit. `fixed` proposes `initial_dt`. Growth, retry, and BDF-ratio
bounds remain explicit for every control kind.

## Remaining accepted values and exclusions

| Key | Accepted value |
| --- | --- |
| `units` | `SI` |
| `thermophysics.data_file` | unique direct-root name ending in `.d` |
| `mesh.data_files[]` | unique direct-root names ending in `.d` |
| `mesh.immersed_boundary` | `null` or exactly `{ "stl_file": "body.stl", "fluid_side": "outside" | "inside" }`; the STL remains a unique direct-root `.stl` file |

The STL may extend beyond the Cartesian domain. Cold product compilation
preserves exact parsed vertices for topology, clips surface quadrature to the
open Cartesian domain, and discards domain-exterior closure facets. An
immersed surface may intersect an external face only when that face is
`symmetry` or `slip`; the compiler then emits an explicit one-sided,
full-quadratic boundary-intersection stencil. Intersections with inlet,
outlet, wall, or NSCBC faces are rejected rather than assigned an implicit
corner model.

The closed case input rejects body-fitted, multiblock, AMR/nonmatching
refinement, constant-density flow, SIMPLE/PIMPLE, reacting physics, other
pressure-corrector counts, and unlisted enum values. Boundary-plan compilation
separately rejects incompatible face combinations and unsupported
shock/supersonic settings; those are not cross-field decisions made by
`CaseCompiler`.

## Cold-path bounds and fingerprint

| Resource | Limit |
| --- | ---: |
| `case.json` | 1 MiB |
| JSON nesting depth | 32 |
| focus regions before normalization | 1024 |
| scalar entries per face | 64 |
| scalar stable name | 255 bytes |
| referenced files | 256 |
| direct-root filename | 255 bytes |
| each referenced file | 64 MiB |
| thermophysical `.d` file | 4 MiB |
| typed MPI wire payload | 256 KiB |

The deterministic nonzero fingerprint hashes a versioned semantic contract
followed by the complete canonical typed model: mesh geometry and normalized
focus regions; limits; turbulence and pressure-reference kinds; all six faces
and all scalar entries in face order; all schemes; every time-control field;
the validated typed thermophysical model; normalized reference kinds/names; and
exact referenced file bytes plus byte counts. Floating-point values are hashed
by their canonical IEEE-754 binary64 bit patterns. JSON whitespace, object-key
order, normalization-equivalent focus input, and `-0.0` therefore do not affect
identity, while any typed setting or referenced-byte change does.

`FieldRegistry::freeze()` publishes one immutable field schema after all
capabilities register. No case-input compilation freezes this registry early;
the production bundle compiler owns the single freeze point.
