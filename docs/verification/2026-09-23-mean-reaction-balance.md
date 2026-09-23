# Independent mean-reaction task distribution

Spatially concentrated finite-rate chemistry left most MPI ranks waiting: on a
128-rank reacting workload, an interval took about 6.9 s on the slowest rank,
3.0 s on average, and less than 0.05 s on the fastest rank.

The built-in Cantera mean-interval path now distributes independent PH/Y requests
cyclically and returns every response to its original cell. The fluid partition,
provider inputs, integrator controls, response checks, equation assembly order,
and final acceptance criteria are unchanged. Caller-owned providers retain their
local call order. Buffers are reserved during freeze and included in reported
reaction workspace bytes. Input preparation and remote evaluation failures are
agreed collectively before results are consumed.

Validation:

- MPI-1/2/4 routing tests cover empty ranks, unequal counts, concentrated work,
  exact result order, inconsistent record shapes, remote failures, and reuse of
  reserved storage.
- Full product reaction tests pass on MPI-1/2, including failed-provider rollback.
- Same checkpoint, mesh, time step, 128 ranks and core binding: a clean baseline
  advance took 172.9301703 s; the production candidate took 116.0512267 s
  (1.49x speedup, 32.9% less time). Both used 16 outer iterations with no retry.
- Native restart comparison checked 340,408,760 values. All current/history
  fields, source histories, face fluxes and pressure-mass references were
  bitwise identical; metadata matched. Physical mass, element and energy
  budgets passed the existing thresholds.

Build scope: changed runtime translation units were rebuilt and linked over the
existing V1.1.0 archive; this was not a fresh full-repository CMake build. The
focused test was also registered in CTest. The measured speedup concerns a
single-node, heterogeneous chemistry workload; communication costs and gains on
other machines require their own measurement. Long-run reacting acceptance is
still pending. Case inputs, geometry and results remain outside the repository.
