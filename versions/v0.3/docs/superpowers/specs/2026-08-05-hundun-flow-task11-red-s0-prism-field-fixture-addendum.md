# Task 11 RED-S0 Prism Field-Fixture Addendum

Status: **FROZEN — THREE-FILE TEST-FIXTURE IMPLEMENTATION AUTHORIZED**

Date: 2026-08-05

## 1. Scope and authority

This addendum resolves one candidate-independent fixture-quality blocker in
Task 11 RED-S0. It does not change any product algorithm, formal error norm,
scientific threshold, body geometry, solver option or Stage 3 capability.

The controlling RED-S0 design and packet remain authoritative except where
this addendum replaces the third real fixture from the unscaled approved
oblique prism with the same geometry carrying the exact MMS factor multiplier
defined below. All historical failed evidence remains valid and retained.

Independent science/requirements and plan/test-quality reviews both passed
before this status was changed to frozen. Implementation remains limited to
the three files and evidence gates below; this status does not accept Task 3
or RED-S0.

## 2. Natural RED and rejected geometry alternative

After the frozen minimum arithmetic/classifier was implemented, every pure and
synthetic observation passed. Sphere, finite cylinder, translated sphere and
inside-sphere cavity passed all four real observations. The approved centred
oblique prism alone failed the full wrapper.

Its 48/96 references are finite and stable. Here and below,
`viscous_force_scale = kDynamicViscosityPaS*kReferenceVelocityMPerS*
surface_area_m2/kReferenceLengthM`. The failed values are approximately:

```text
F_net                 = 2.0596547622679542e-18 N
minimum A_abs         = 1.7229953960060171e-16 N
4096*epsilon*viscous_force_scale
                       = 2.7575879357754987e-15 N
max_abs(F48-F96)      = 3.7902301305540802e-31 N
abs(T_cert48-T_cert96)= 1.3016204936146695e-29 N
```

Thus the failure is non-degeneracy, not quadrature instability. A one-shot
candidate precommitted before measurement reused the existing translated
sphere displacement `{+0.013,-0.009,+0.007}`. It remained equally degenerate.
No translation or geometry parameter scan followed.

## 3. Mathematical field normalization

For the prism, the current manufactured scalar uses

```text
f = product_i(local_i^2 - h_i^2)
psi = envelope * asymmetry * f^2
u = fixed divergence-free linear map of grad(psi).
```

At each face, the extra transverse quadratic factors make `grad(f)` and hence
the viscous traction roundoff-scale for this small body. Translation cannot
change those local factors.

Normalize the extra planar quadratic factor by the permutation-symmetric mean
pair area, against the already fixed reference length:

```text
A_char = (h_x*h_y + h_y*h_z + h_z*h_x) / 3
m = kReferenceLengthM^2 / A_char
  = 3*kReferenceLengthM^2
    / (h_x*h_y + h_y*h_z + h_z*h_x).
```

For `{h_x,h_y,h_z}={0.14,0.11,0.09} m`:

```text
A_char = 0.012633333333333335 m^2
m      = 79.155672823218993
m^2    = 6265.62054009649.
```

The formula was frozen before its analytic prediction was calculated. It is
not selected from a product result, error order or parameter sweep.

Both the prism factor and its three gradient components receive `m` before
`psi_gradient()` uses them. The factor is a `Dual3`; multiplying it scales its
stored value, gradient and Hessian, not only its scalar value. Therefore
`psi`, velocity, velocity gradient,
viscous traction, signed viscous force, absolute-component integrals and
viscous RMS certificate scale exactly by `m^2`; pressure is unchanged. The
wall velocity remains zero, divergence-free construction is unchanged and the
implicit-factor sign/zero set is unchanged.

## 4. Exact test-support contract

Extend test-only `BodySpec` with:

```cpp
double prism_mms_factor_multiplier{1.0};
```

Expose:

```cpp
BodySpec force_certified_oblique_prism();
```

The helper returns `approved_body(BodyKind::oblique_rectangular_prism)` and
sets only `prism_mms_factor_multiplier` to the exact expression in Section 3.
The centre, axes and half lengths remain byte-for-byte the approved values.

Only the prism branch of `factor_gradient()` consumes the coefficient. Its
returned `Dual3` factor, including its automatic-differentiation derivatives,
and each returned gradient-component `Dual3`, likewise including its stored
derivatives, are multiplied once as whole objects by the same coefficient.
Neither path may scale only `.value`. Generic approved-prism MMS/surface
oracles keep the default multiplier 1.

The existing focused selector's third body and `formal_case()`'s `prism_*`
branch must both call `force_certified_oblique_prism()`. Append one named
observation after the existing frozen vector:

```text
prism_fixture_route_exact
```

It compares the complete BodySpec returned by the formal prism route with the
focused third body, including the multiplier. This check constructs a case but
does not run a manufactured solve. The third diagnostic label becomes
`force_certified_oblique_prism`; all earlier observation names and their order
remain unchanged.

## 5. Pre-implementation analytic prediction

Applying the exact `m^2` transformation to the existing viscous primitives
and recomputing `T_abs` from the scaled components gives:

```text
48:
  F_viscous = { 5.7958318674513133e-15,
                -1.2905015183971273e-14,
                -9.2190883808026967e-15 } N
  A_abs     = {1.7618301502019749e-12,
               1.4062747410223192e-12,
               1.0795635343707129e-12} N
  T_abs     = 4.2476684255950067e-12 N
  T_cert    = 6.9452126242715438e-12 N

96:
  F_viscous = { 5.7958318674491266e-15,
                -1.2905015183973647e-14,
                -9.2190883808041056e-15 } N
  A_abs     = {1.7618301502019864e-12,
               1.4062747410223060e-12,
               1.0795635343706987e-12} N
  T_abs     = 4.2476684255949914e-12 N
  T_cert    = 6.9452126242714630e-12 N
```

At 96 points:

```text
F_net                       = 1.2905015183973647e-14 N
4096*epsilon*viscous_force_scale
                             = 2.7575879357754987e-15 N
1024*sqrt(epsilon)*T_cert   = 1.0597553442797032e-16 N
omitted dimensionless error = 4.2562714986720775e-12
reversed dimensionless error= 8.5125429973441550e-12
D_net                       = 0.00303814090248016
D_cert                      = 0.0018581166455400426
max_abs(F48-F96)            = 2.3744713247152455e-27 N
abs(T_cert48-T_cert96)      = 8.0779356694631609e-26 N
classifier accepts          = true
```

The prediction calls the existing shared classifier. It is not a substitute
for the post-implementation analytic path. The first prediction harness
incorrectly scaled the old L1 scalar independently and correctly failed the
exact composition predicate; that harness RED is retained. The corrected
prediction scales the three primitives then recomputes fixed `(x+y)+z`.

## 6. Allowed files and forbidden changes

Only these tracked source files may change:

```text
tests/support/stage3_test_contracts.hpp
tests/support/stage3_mms.cpp
tests/numerical/test_laminar_ibm_order.cpp
```

`CMakeLists.txt`, product sources/headers, formal error rows, thresholds,
solver controls, body geometry and all other fixtures are forbidden. The
inherited dirty candidate and empty index must be preserved; no worker stages,
commits, stashes, resets, restores or cleans.

## 7. TDD, acceptance and mutations

1. Reproduce the current focused RED and prism-only false observations.
2. Add only the field/helper/factor scaling and shared focused/formal routing.
3. Build Debug `-j2`; run the focused selector. Every existing observation and
   `prism_fixture_route_exact` must be true.
4. Run the existing oracle plus focused selector in Debug and Release with the
   frozen clang runtime path.
5. Report actual 48/96 full analytic values from the implemented path. In
   addition to the focused log, run a serial `/tmp` evidence diagnostic that
   calls only the implemented `analytic_force_reference()` at 48/96 and the
   shared classifier, prints all primitives/floors/ratios/stability at
   `setprecision(17)`, and records its exact source, binary, argv, exit and log
   SHA-256 values. This is evidence formatting, not a second algorithm. Do not
   change the formula if values differ by ordinary evaluation roundoff from
   the prediction; stop if any frozen predicate fails.
6. In separate disposable exact-candidate worktrees, require successful build
   and focused failure for each mutation:

   ```text
   F-M0 helper multiplier forced to 1
     -> force_certified_oblique_prism_wrapper_accepted
   F-M1 omit multiplier from returned Dual3 factor/derivatives only
     -> force_certified_oblique_prism_wrapper_accepted
   F-M2 omit multiplier from returned gradient only
     -> force_certified_oblique_prism_wrapper_accepted
   F-M3 omit the numerator factor 3 from the exact coefficient
     -> force_certified_oblique_prism_wrapper_accepted
   F-M4 focused third body reverted to approved prism
     -> force_certified_oblique_prism_wrapper_accepted
   F-M5 formal prism branch reverted to approved prism
     -> prism_fixture_route_exact
   ```

   At a wall where `f.value==0`, differentiating the existing `2*f*g` term
   gives the traction-bearing `2*df*g` term. The returned Dual3 factor supplies
   `df` and the separately returned gradient supplies `g`; each multiplier is
   therefore independently required. Omitting either leaves only `m` rather
   than `m^2` traction scaling. The predicted corresponding
   `F_net` is about `1.63e-16 N`, below the unchanged
   `2.76e-15 N` floor.

7. The main agent reviews all deltas and mutation evidence. Only then is Task 3
   GREEN and the original RED-S0 Task 4 mutation matrix authorized.

No product fast/screen/formal run is allowed during this fixture subcluster.

## 8. Evidence invalidation and downstream consequence

Old prism formal numerical logs cannot be promoted after this field fixture is
adopted because the exact velocity/source amplitude changes. Geometry,
classification and topology fingerprints are expected to remain unchanged;
the formal prism error/force rows must be rerun through the normal later
fast/screen/acceptance chain. Other bodies and generic approved-prism oracles
remain valid.

## 9. Pre-edit review closure

The final reviewed-content candidate had SHA-256
`385b1bd5b9a565b654f48eb0988d2d6b8daa6169c0801d3723c12dc2b0757ad3`;
the matching implementation packet candidate had SHA-256
`35cec010353428e63328cb226864189cb1121cfd6fc794a26a3bc92d432a7180`.

Independent science/requirements review passed with no Critical, Important or
Minor defect after verifying the whole-Dual3 `m^2` derivation, invariants,
precommit-before-prediction chronology and fixed-L1 harness RED. Independent
plan/test-quality review passed with no Critical or Important defect after
verifying allowed hashes/files, direct field-by-field routing, high-precision
implemented-path evidence, F-M0--F-M5 sensitivity and NUL-safe worktree
invariance. Its implementation-time Minor checks are: never use `memcmp` for
BodySpec padding, compare every named field directly, and change only the third
body's diagnostic prefix while preserving the existing observation order.

Post-edit evidence, Debug/Release GREEN and mutations remain mandatory. Their
absence at freeze time is expected and is not Task 3 acceptance.
