# Arena budget without optional payloads

The product compiler previously checked arena plus model bytes only when a
spray, ESF, derived-output, ICCG or passive-scalar payload was present. A CN/BE
case without these payloads could therefore allocate its state arena despite a
budget below that arena's size. The check now runs whenever arena compilation
succeeds, before allocating state layers. Existing checked byte arithmetic and
collective failure handling remain unchanged.

The regression uses CN/BE, an analytic gas, an adiabatic wall, a pressure outlet,
periodic lateral faces and no optional models. On 16^3 cells, even seven owned
primitive components exceed 32 KiB at each tested decomposition. It checks that
the low-budget candidate fails without publishing a plan, that an existing
independent plan remains unchanged, and that sufficient-budget compilation
still succeeds. CTest registers `v04_product_arena_budget_mpi_{1,2,4}`.

Validation with the focused release-archive overlay toolchain:

- Before the fix, the final regression returned failure because low-budget
  compilation incorrectly succeeded (`status=0/0`, sufficient-budget retry=1).
- After the fix, MPI 1/2/4 all passed: `allocation_failure/10204`, retry=1.
- An independent real-NASA hydrogen PaSR/IBM box, MPI2, three accepted steps,
  produced byte-identical per-rank exported current/previous fields, rate
  histories, face fluxes, metadata and step diagnostics relative to the old
  core at the same sufficient budget.

This is a bounded memory-admission fix, not a full RSS cap or a claim that every
transient allocation is budgeted. It makes no performance or long-run CFD
acceptance claim. The case data and private predictor experiments remain local.
