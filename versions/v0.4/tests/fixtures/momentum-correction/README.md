# Momentum correction regression

Synthetic 8 mm box, 32³ cells, with a centred 2 mm cubic solid, periodic transverse boundaries, mass-flow inlet and pressure outlet. The STL is a generated cube, not an application geometry. A non-reacting 300 K mixture supplies temperature-dependent transport; no chemistry mechanism is loaded.

The inlet drives a finite velocity from a zero-velocity initialization. With all three independent outer reference gates at `1e-12`, the total-velocity linear solve can accept a residual relative to its large transient RHS while the terminal momentum equation still fails. Solving the velocity correction must let the actual product driver accept the step without relaxing those gates.

This is a numerical regression fixture, not a validated combustion case or performance benchmark.
