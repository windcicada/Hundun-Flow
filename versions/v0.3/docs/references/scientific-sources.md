# Scientific and runtime sources

The Stage 4--6 Linux CPU v1 source catalog is maintained in
`docs/references/2026-08-09-hundun-flow-stage4-6-reference-catalog.md`. That
catalog supersedes this file for reacting-flow, Cantera packaging, ESF/TPDF,
TCR and spray decisions. The Stage 0--3 entries below remain the historical
authority for capabilities already accepted under those stages.

HUNDUN-FLOW uses published equations and standards, not external
implementation source. The approved references and the decisions they support
are:

1. R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems* (2002),
   DOI [10.1017/CBO9780511791253](https://doi.org/10.1017/CBO9780511791253):
   conservative finite-volume formulation.
2. B. van Leer, “Towards the Ultimate Conservative Difference Scheme. V. A
   Second-Order Sequel to Godunov's Method” (1979), DOI
   [10.1016/0021-9991(79)90145-1](https://doi.org/10.1016/0021-9991(79)90145-1):
   MUSCL reconstruction.
3. B. van Leer, “Towards the Ultimate Conservative Difference Scheme. IV. A
   New Approach to Numerical Convection” (1977), DOI
   [10.1016/0021-9991(77)90095-X](https://doi.org/10.1016/0021-9991(77)90095-X):
   monotonized slope limiting.
4. S. Gottlieb and C.-W. Shu, “Total Variation Diminishing Runge-Kutta
   Schemes” (1998), DOI
   [10.1090/S0025-5718-98-00913-2](https://doi.org/10.1090/S0025-5718-98-00913-2):
   SSPRK2 time integration.
5. MPI Forum, *MPI: A Message-Passing Interface Standard, Version 3.1*
   (2015), [standard report](https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report.pdf):
   communicator, Cartesian topology, and nonblocking communication contracts.
6. C. M. Rhie and W. L. Chow, “Numerical Study of the Turbulent Flow Past an
   Airfoil with Trailing Edge Separation” (1983), DOI
   [10.2514/3.8284](https://doi.org/10.2514/3.8284): input for future
   collocated pressure-velocity coupling, not a Stage 1 capability.
7. R. I. Issa, “Solution of the Implicitly Discretised Fluid Flow Equations by
   Operator-Splitting” (1986), DOI
   [10.1016/0021-9991(86)90099-9](https://doi.org/10.1016/0021-9991(86)90099-9):
   input for future PISO pressure correction, not a Stage 1 capability.
8. P. D. Thomas and C. K. Lombard (1979), DOI
   [10.2514/3.61273](https://doi.org/10.2514/3.61273): metric identities and
   geometric-conservation-law checks for independent curvilinear-grid work.
9. J. B. Bell and D. L. Marcus (1992), DOI
   [10.1016/0021-9991(92)90011-M](https://doi.org/10.1016/0021-9991(92)90011-M):
   variable-density projection-method background.
10. S. K. Choi (1999), DOI
    [10.1080/104077899274679](https://doi.org/10.1080/104077899274679):
    time-step-consistent momentum-interpolation background.
11. S. Zhang, X. Zhao, and S. Bayyuk (2014), DOI
    [10.1016/j.jcp.2013.11.006](https://doi.org/10.1016/j.jcp.2013.11.006):
    transient momentum interpolation and checkerboard-suppression background.

These citations support independent derivation and tests. They do not
authorize copying an implementation or treating software source as an
algorithm reference.
