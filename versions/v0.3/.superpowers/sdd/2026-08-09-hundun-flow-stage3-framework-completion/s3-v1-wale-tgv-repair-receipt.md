# S3-V1 WALE Taylor--Green repair receipt

状态：`ACCEPTED`

accepted parent：`4e0abc2cbb1e0731dc0234903d11f8707997c6b2`

## Result

The bounded repair fixes both defects exposed by the first formal scientific
run without changing any frozen scientific threshold, selector, grid, final
time, solver, WALE control, PISO count, Task 11 authority, rollback, Restart or
MPI consistency rule.

- Body-fitted WALE now exchanges the cell SGS-viscosity scratch field before
  face interpolation. The unchanged 24-cubed 1/2/4-rank formal comparison is
  GREEN at the original `5e-12` threshold.
- Uniform-spacing momentum face reconstruction now consumes the declared
  cell-average function space with its symmetric four-cell finite-volume
  interpolation. Warped/non-uniform and MC-limited paths remain unchanged.
- The unchanged 12/24/48 TGV selector is GREEN with both independent
  Richardson segments still required to be `>= 1.8`.

## TDD and review evidence

Executable RED was observed before implementation:

- `test_cell_centered_fvm_1_rank` failed at the new quadratic cell-average
  face oracle;
- `test_wale_taylor_green_convergence_1_rank_formal` failed at the original
  velocity-order assertion;
- the original 2/4-rank TGV rows failed decomposition comparison before the SGS
  halo repair.

Focused GREEN on the final source bytes:

```text
test_cell_centered_fvm_{1,2,4}_rank                         3/3 PASS
test_wale_taylor_green_convergence_1_rank_formal            PASS
test_wale_taylor_green_n24_{1,2,4}_rank_formal              3/3 PASS
test_wale_body_fitted_{1,2}_rank                            2/2 PASS
test_taylor_green_piso_1_rank                               PASS
cell-centered FVM header/source policy + Stage 3 policy     3/3 PASS
```

The accepted body-fitted field/report fingerprints were updated only after the
new deterministic values were observed; the repeated-attempt identity,
rollback, two-corrector and 2-rank contracts remain GREEN.

The main agent reviewed the formula, function-space semantics, periodic-pair
symmetry, two-cell ghost requirement, uniform-spacing gate, warped/MC fallback,
all callers, conservation and MPI behavior. The immersed Task 11 reconstruction
and force path are separate and unchanged. No private source or research data
was accessed; no research process was inspected or signalled; no push,
publication, 96-cubed run or Stage 4--6 action occurred.

The prior candidate and its V0/V1 evidence are retained as rejected history.
The repair commit creates a new candidate only after a complete restarted V0.

提交 subject：`fix: repair Stage 3 WALE acceptance defects`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`
