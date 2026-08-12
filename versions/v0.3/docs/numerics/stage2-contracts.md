# Stage 2 numerical and measurement contracts

This document freezes vocabulary and acceptance values for Stage 2. It does
not claim that the Stage 2 mesh, linear solvers, flow solver, Checkpoint v2,
or execution backends have been implemented.

## Face orientation and mass flux

For an internal face, the face area vector `S_f` points from owner to
neighbour. For a physical boundary face, it points outward from the owner.
Positive `mdot_f` therefore means owner-to-neighbour flow on an internal face
and outward flow on a boundary face. A shared internal face is stored exactly
once. The neighbour consumes the exact negative of that stored flux.

## Linear residual

The linear residual is `r=b-Ax`. The success norm is the project-owned FP64
L2 norm, and the stopping contract is

```text
||r||_2 <= max(atol, rtol*||b||_2).
```

A recursively updated residual may be reported as a diagnostic, but it cannot
establish success. Successful completion requires an independent FP64
recomputation of `b-Ax` and the same stopping test.

## Finite-volume equation residual

For cell `i`, the raw residual `R_i` is the signed sum of the discrete time,
face-flux, diffusion, and source contributions for the equation. Its unit is
the conserved quantity per time. Let `a_i` be the sum of the absolute values
of those same contributions. Define

```text
R2 = sqrt(sum_i R_i^2)
A2 = sqrt(sum_i a_i^2)
normalized_L2 = R2/A2.
```

If `R2` and `A2` are both exactly zero, `normalized_L2` is zero. If only `A2`
is zero, it is positive infinity. Global sums and norms use FP64 and include
every owned cell exactly once.

## Global conservation defect

For a globally conserved quantity `Q`, the signed defect over a step is

```text
D = Q_(n+1) - Q_n
    + dt*(sum boundary outward flux - sum volume source).
```

The relative-defect denominator is the maximum of `|Q_n|`, `|Q_(n+1)|`, the
time-integrated absolute boundary flux, the time-integrated absolute source,
and `numeric_limits<double>::min()`. An exact zero defect with all physical
scales exactly zero reports a zero relative defect.

## Approved analytic mapping

For `xi,eta,zeta in [0,1]`, the mapping is exactly

```text
X = x0 + Lx[xi + ax sin(2*pi*xi) sin(pi*eta) sin(pi*zeta)]
Y = y0 + Ly[eta + ay sin(pi*xi) sin(2*pi*eta) sin(pi*zeta)]
Z = z0 + Lz[zeta + az sin(pi*xi) sin(pi*eta) sin(2*pi*zeta)]
```

Each amplitude satisfies `|a|<=0.02`. The amplitude bound alone never
substitutes for runtime checks of positive Jacobian, positive cell volume,
shared-face area-vector reciprocity, and per-cell closure.

## Frozen acceptance cases and thresholds

| Acceptance item | Fixed case and threshold |
|---|---|
| Geometry | Warped amplitude `[0.02,-0.015,0.01]`; closure `<=256*eps*sum(|S|)`; shared face area vectors reciprocal; every Jacobian and volume positive |
| Linear | 63x63 SPD and nonsymmetric manufactured systems; `atol=1e-12`, `rtol=1e-10`; reported residual versus independent recompute differs by `<=64*eps*max(1,||r||)` |
| Poisson | Periodic sine solution, `16^3/32^3/64^3`; L2 convergence order `>=1.8`; zero mean `<=1e-12` |
| Taylor--Green | Periodic `[0,2*pi]^2`, 4 cells in z; spatial and temporal convergence order `>=1.8`; relative total-mass error `<=5e-12` |
| Checkerboard | Initial parity-pressure amplitude 1; after two correctors parity amplitude `<=1e-8`, normalized continuity L2 `<=1e-10`, corrector count exactly 2 |
| Density wave | `rho=1+0.2 sin(2*pi*(x-t))`, grids `32/64/128`; L1 order `>=1.8`; mass error `<=5e-12`; density remains positive |
| Variable-density vortex | Approved steady analytic field; L2 order `>=1.8`; 1/2/4-rank maximum field difference `<=5e-12*max(1,||q||_infinity)` |
| Ideal gas | Relative errors in `h/(cp*T)` and `rho*R*T/p0` `<=1e-12`; closed-domain mass error `<=5e-12`; nonpositive state rejected collectively |
| Final residuals | Normalized continuity L2 `<=1e-10`; normalized momentum/enthalpy/scalar L2 `<=1e-9`; global conservation relative error `<=5e-11` |
| Restart v2 | With the same partition, continuous and resumed fields, history, next `dt`, order, and final flux are bitwise identical; every corrupt or mismatched case is explicitly rejected |
| MPI | Every applicable case covers 1/2/4 ranks; no hang; failure category and lowest failing rank agree |
| Performance | Exact counters equal independent formulas; every raw timing is positive and finite; incompatible metadata produces only `incomparable` |

The approved steady variable-density vortex field is

```text
psi = sin(x) sin(y)
rho = 1 + 0.1 psi
u = sin(x) cos(y)
v = -cos(x) sin(y)
w = 0
pi = 0.
```

## Performance measurement and scaling

For each repetition, a rank sample is the measured-phase elapsed time divided
by the common measured-step count. The repetition value is the maximum rank
step time. `T_p` is the median of the repetition maxima; for an even number of
repetitions the median is the arithmetic mean of the two middle values. Every
raw rank sample is retained.

Strong scaling uses the same global problem:

```text
S_p = T_1/T_p
E_p = T_1/(p*T_p).
```

Weak scaling uses the same per-rank problem:

```text
W_p = T_1/T_p.
```

The manual matrix is ranks `1/2/4`, strong global `64^3`, weak per-rank
`32^3`, 5 warmup steps, 20 measured steps, and 5 repetitions. Portable CI
hard-gates only exact counters. Wall time, RSS, bandwidth, and throughput are
recorded but have no portable hard threshold.

Comparisons require equal hardware and node identity, MPI identity, compiler
identity, complete compiler and link options, rank placement and affinity,
resolved-case or approved scaling-case fingerprint, numerical tolerance
contract, measurement method and policy, and execution backend. Rank count
and process grid may differ only for a declared strong- or weak-scaling
comparison. Strong scaling keeps the global problem fixed; weak scaling keeps
the per-rank problem fixed. Any required mismatch makes the result
`incomparable`; it is never interpreted as a performance-regression pass or
failure.
