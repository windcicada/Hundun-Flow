# PaSR stiff-source depletion

The timescale fraction alone does not make an instantaneous chemical source
stable over a finite flow step. A 128-rank continuation reproduced rejection
10217 in the species target solve. Its worst row had accepted mass fraction
3.4145790e-7, density 0.12398668 kg/m3, source -3.2794336 kg/m3/s and attempted
duration 2.7364084e-8 s. The source consumed 2.11967 times the available species
inventory. The bounded search approached zero while the unconstrained row
required -3.81412e-7. More linear iterations could not admit a physical solution.

CN/BE PaSR now integrates the reaction zone once per active cell per attempted
flow interval, before fluid coupling, and freezes the common-kappa interval
increment. The stored accepted density supplies its conservative source.
The same kappa scales all species; no clipping, source floor or convergence
relaxation is introduced. This follows the reaction-zone recombination in
COAST `c4c1f788c0ca7359653a5acd621dd355bb5c7db3` `SRC.Coast/fieldpdf.F90`,
while retaining Hundun's explicit timescale definitions and equation audits.
It is not a full COAST PDF-method equivalence claim.

The failure was minimized to the production driver with uniform A -> B,
rate 2/s, a fixed 1 s step and nonzero mixing time. The old implementation
rejected its first step with 10217. MPI-1 and MPI-2 now complete three steps;
the maximum composition error against the independent convex-exponential
oracle is below 1.33e-13. Original energy and independent composition/element
budgets pass. Tests assert one interval call per active cell per attempt, reject
a missing external interval provider, and verify complete rollback after
partial query/interval failures. Zero-reaction states and MPI-2 WALE/retry
checks also pass. Existing mean-interval and BE source paths are retained.

The source-method fingerprint changes: old instantaneous-source native
checkpoints are not silently restored as this method. Reinitialization or an
explicit physical-state import is required. As with the preceding change,
validation binaries rebuild the modified translation units over the V1.1.0
archive. The production CLI's 27 model-selection/restart checks and real-Cantera
pure-species initialization/restart also pass.

The reinitialized real-geometry 128-rank run accepted its first two steps without
retry, with nonzero reaction. Maximum relative physical energy error was
1.60591e-9 and element error 7.89719e-12. The first advance took 23.5181 s;
interval integration has a cost compared with the unstable instantaneous source.
Larger-step throughput and long-test acceptance remain pending; the small stiff
test establishes the depletion fix, not those goals.
