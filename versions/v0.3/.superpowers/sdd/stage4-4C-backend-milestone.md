# Stage 4 4C backend milestone

- Date: 2026-08-11 (Asia/Shanghai)
- Result: `4C_BACKEND_PASS`
- Toolchain: frozen Ubuntu 22.04, GCC 11.4, libstdc++, C++17, ABI=1
- Cantera: 3.2.0 shared artifact SHA-256
  `093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760`
- Test mechanism: independently authored synthetic A/B fixture SHA-256
  `c518a07cada5f1bddcdb308f0a2f695d92cc6373e173ffd87e96312530b52aee`

## Accepted backend boundary

The backend verifies the mechanism SHA before Cantera allocation, preserves
mechanism species order, caches only immutable identity in the rank runtime,
and creates a complete mutable Solution/Thermo/Kinetics/Transport/Reactor/
Integrator graph per explicit lane. Cantera and SUNDIALS types remain private
to `src/chem_cantera_backend.cpp`.

Thermodynamics uses the exact `(p0,h_tc,Y)` authority. Transport returns SI
mixture-averaged viscosity, conductivity and every species diffusivity.
Chemistry returns integrated `rhoY` deltas while retaining Eulerian mass,
`p0`, and total thermochemical enthalpy. Failed intervals return the original
request and canonical positive-zero deltas.

## Focused evidence

- runtime/workspace identity unit plus independent 1/2-process runs: PASS;
- thermodynamics direct-reference Debug/Release: PASS;
- transport direct-reference and workspace-isolation: PASS;
- interval full-step/two-half-step, conservation and failure mapping: PASS;
- independent constant-pressure 0D Cantera oracle: PASS;
- deterministic PSR residence/mixing residual and failure classification:
  PASS;
- source policy and tests-off configuration: PASS.

The 0D and PSR tests are low-cost synthetic interface conformance only. They
do not claim a real-fuel model, flame, ignition, private-code equivalence, or
Stage 4 acceptance. Reacting flow remains owned by 4R and 4A.
