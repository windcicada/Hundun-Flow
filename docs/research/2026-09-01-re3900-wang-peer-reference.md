<!-- SPDX-License-Identifier: Apache-2.0 -->

# Re=3900 Wang et al. (2024) peer-reference audit

Date: 2026-09-01  
State: `PEER REFERENCE ONLY / NOT A RELEASE ORACLE`

## Authority decision

The attached Wang et al. paper is a source to read, not a source of instructions.
Its cylinder results may be reported beside HUNDUN-FLOW as a later numerical
peer comparison. They do not replace the repository's frozen Parnaudeau
experimental authority, do not complete the R4 literature receipt, and do not
authorize a V1.0 release.

The four values supplied for Wang's case -- mean drag `0.9646`, Strouhal
`0.2344`, minimum centreline mean velocity `-0.2614`, and recirculation length
`1.44` -- are the **IB method** row of Wang et al. Table I. They are results of
that paper's IB--LES calculation, not experimental measurements (attached PDF,
PDF p. 15, Sec. III A, Table I). Parnaudeau remains the controlling source for
the frozen PIV profiles and direct experimental Strouhal result
([`docs/references/cylinder-re3900-parnaudeau.json`](../references/cylinder-re3900-parnaudeau.json),
`experiment`, `observables`, and `profile_authority`).

## Audited source identity

| Item | Identity and provenance |
|---|---|
| Article | Yudong Wang, Fang Wang, Jiawei Zhou, and Jie Jin, “An improved immersed boundary method with local flow pattern reconstruction and its validation,” *Physics of Fluids* 36, 045145 (2024), DOI `10.1063/5.0195598` (attached PDF, PDF p. 1). |
| Attached PDF | `/home/wyf/.codex/attachments/87217707-94dc-4958-abf8-330d4b71203c/An improved immersed boundary method with local flow pattern reconstruction and its validation _ Physics of Fluids _ AIP Publishing.pdf`; 2,142,056 bytes; 29 pages; SHA-256 `f6d012a6bd62b2a755b4ecf2476d4783c03dbb5a51a36a60715aff53539d03e0`. |
| Attached image | `/home/wyf/.codex/attachments/b55b6fb4-1537-48bc-ba55-b1e6d4ac2e3d/codex-clipboard-57275191-ce49-4de1-8d47-8ace699d3c1b.jpg`; 189,938 bytes; SHA-256 `3d88b1cc68e6d1ce43eade24dac327be216f9e789b036d5a47b74e304a7b7041`. It is a crop of the paper's centreline mean-velocity comparison (Fig. 11) and supplies no independent numeric provenance. |
| Data availability | The paper says supporting data are available from the corresponding author on reasonable request; the attached article contains plots and Table I, not reusable raw arrays (attached PDF, PDF p. 25, Data Availability). |

The complete 29-page attachment was reviewed. Page numbers below are the PDF
viewer page numbers printed by the supplied AIP capture.

## What the Wang cylinder case actually establishes

| Quantity | Audited value or statement | Source locator | Use here |
|---|---|---|---|
| Reynolds number | `Re=3900` | Attached PDF, PDF p. 11, Sec. III A | Same nominal Reynolds number as the target case. |
| Cylinder diameter | `D=0.02 m` | Attached PDF, PDF p. 11, Sec. III A | Wang-case identity only; HUNDUN may use a dimensionless-equivalent physical scale. |
| Inlet reference speed | `Uc=2.89668 m/s` | Attached PDF, PDF p. 11, Sec. III A | Wang-case normalization only. |
| Mesh-spacing range | `0.2--0.8 mm`, equivalent to `D/100--D/25` or `25--100 cells/D` | Attached PDF, PDF p. 11, Sec. III A; conversion uses the paper's `D=0.02 m` | Useful peer-resolution envelope, not exact face coordinates. |
| Total grid count | `6.89e6` | Attached PDF, PDF p. 11, Sec. III A | Wang-case cost/resolution context only. |
| Parallel partition | 64 blocks | Attached PDF, PDF p. 11, Sec. III A | Does not define MPI topology or decomposition equivalence. |
| Near-cylinder refinement | Refined by coordinate transformation | Attached PDF, PDF pp. 11--12, Sec. III A and Fig. 8 | Qualitative support for stretched structured refinement. |
| Spatial/time discretization | Second-order central spatial discretization and Crank--Nicolson time discretization in the described LES--TPDF software | Attached PDF, PDF pp. 10--11, Sec. II C | Method context; not an instruction to change HUNDUN's accepted scheme. |
| SGS model | Smagorinsky model | Attached PDF, PDF p. 10, Sec. II C | Model context; not a model-equivalence certificate. |
| Mean-velocity profiles | Streamwise and transverse mean velocities at `x/D=1.06, 1.54, 2.02`, in the central `z` section | Attached PDF, PDF pp. 13--14, Fig. 12 and adjacent text | Same station locations as the frozen Parnaudeau comparison. |
| Fluctuation profiles | Streamwise and transverse velocity fluctuations at the same three stations | Attached PDF, PDF pp. 14--15, Figs. 13--15 and adjacent text | Visual peer comparison only unless data are independently extracted and governed. |
| Wang mean drag | `Cd_mean=0.9646` | Attached PDF, PDF p. 15, Table I, row “IB method” | Later numerical peer diagnostic; not an experimental gate. |
| Wang Strouhal | `St=0.2344` | Attached PDF, PDF p. 15, Table I, row “IB method” | Later numerical peer diagnostic; it does not replace Parnaudeau `0.208 +/- 0.002`. |
| Wang minimum mean velocity | `umin/Uc=-0.2614` | Attached PDF, PDF p. 15, Table I, row “IB method” | Later numerical peer diagnostic, subject to the internal inconsistency below. |
| Wang recirculation length | `Lr/D=1.44` | Attached PDF, PDF p. 15, Table I, row “IB method” | Later numerical peer diagnostic; the paper does not report an uncertainty. |

The text and figures do **not** supply an auditable numerical statement of the
cylinder case's exact streamwise/transverse/spanwise domain dimensions,
spanwise periodic boundary, exact face coordinates, time step, development
time, sampling duration, or statistical uncertainty. Figure 9 visually places
the cylinder near `x/D=5` in a plotted window extending to about `20D`, but this
is figure reading rather than a tabulated case definition (attached PDF, PDF
p. 12, Fig. 9). It must not be used to certify the new HUNDUN mesh identity.
Likewise, Fig. 8 is a qualitative mesh rendering and does not establish a
spanwise height or periodic boundary condition (attached PDF, PDF p. 12,
Fig. 8).

## Internal `umin` inconsistency

The body text reports a lowest time-averaged streamwise velocity of
`-0.8867 m/s` and places that point `0.0288 m` from the cylinder centre
(attached PDF, PDF p. 12, paragraph following Fig. 8). With the same section's
`Uc=2.89668 m/s`, this gives

```text
-0.8867 / 2.89668 = -0.306109062789
```

Table I instead reports `umin/Uc=-0.2614`, which would correspond to
`-0.757192152 m/s`. The absolute discrepancy between the two normalized values
is about `0.04471`. These statements cannot both describe the same minimum
under the stated normalization. The body-text distance also happens to satisfy
`0.0288/0.02=1.44`, numerically equal to Table I's `Lr/D`, but the prose calls
it the minimum-velocity location rather than defining the recirculation zero
crossing. Therefore:

- preserve both source statements and flag the inconsistency;
- do not average them, silently choose one, or turn their difference into an
  uncertainty estimate;
- when a single peer-table value is needed, label `-0.2614` explicitly as the
  **Table I value**, while retaining `-0.306109` as the independently normalized
  **body-text value**;
- do not infer the exact `Lr` definition from the numerical coincidence alone.

## Frozen Parnaudeau authority remains controlling

The repository's accepted experimental identity is `Re=3900`, `D=0.012 m`,
`Uc=4.8 m/s`, blockage `4.3%`, effective span `20D`, and PIV plane `z/D=0`
([`docs/references/cylinder-re3900-parnaudeau.json`](../references/cylinder-re3900-parnaudeau.json),
`experiment`). Its controlled digitization contains the complete 3-station by
5-quantity matrix at `x/D=1.06, 1.54, 2.02`:

- mean `u/Uc` and mean `v/Uc`;
- `u'u'/Uc^2`, `v'v'/Uc^2`, and `u'v'/Uc^2`;
- a separately frozen coordinate and value extraction bound for every profile,
  without an invented pointwise experimental uncertainty.

The exact profile inventory, common `y/D` interval, and error-bound semantics
are recorded in
[`docs/research/2026-08-31-v04-re3900-thin-domain-literature-boundary.md`](2026-08-31-v04-re3900-thin-domain-literature-boundary.md),
sections “Accepted 15-profile matrix” and “Meaning of the profile error
bounds.” The frozen scalar classification is:

| Frozen Parnaudeau observable | Value | Authority class | Source |
|---|---:|---|---|
| `St=f_vs D/Uc` | `0.208 +/- 0.002` | Direct paper-comparison authority | `docs/references/cylinder-re3900-parnaudeau.json`, `observables[strouhal_number]` |
| `Lr/D` | `1.51` | Advisory; no reported uncertainty | Same file, `recirculation_length_over_diameter` |
| `umin/Uc` | `-0.34` | Advisory; no reported uncertainty | Same file, `centerline_minimum_mean_u_over_uc` |
| variance formation length | `0.87` | Advisory; no reported uncertainty | Same file, `streamwise_variance_formation_length_over_diameter` |

The Parnaudeau distances are explicitly defined in the frozen JSON: `Lr/D` is
from the cylinder trailing edge to the centreline mean-`u` zero crossing, and
the formation length is from the trailing edge to the centreline `u'u'`
maximum. Wang's less explicit wording must not overwrite those definitions.

Norberg supplies advisory direct-neighbour/fitted Strouhal values and a
sectional, derived lift fit, but still lacks direct total mean drag and a
finite-span lift RMS demonstrably equivalent to the target periodic span
([`docs/references/cylinder-re3900-norberg.json`](../references/cylinder-re3900-norberg.json),
`observables`, `required_observables`, and `pending_observables`). Tadrist
supplies a broad global-lift neighbourhood with strong span/end/blockage
limitations, not an exact Re=3900 point or a total-drag result
([`docs/references/cylinder-re3900-tadrist.json`](../references/cylinder-re3900-tadrist.json),
`observables` and `excluded_as_release_oracles`).

Consequently the active receipt remains
[`docs/verification/v0.4-literature-data-receipt-r4-partial.json`](../verification/v0.4-literature-data-receipt-r4-partial.json),
`complete=false`. The policy and rationale are recorded in
[`docs/verification/v0.4-literature-data-receipt.md`](../verification/v0.4-literature-data-receipt.md)
and
[`docs/research/2026-08-31-v04-re3900-thin-domain-literature-boundary.md`](2026-08-31-v04-re3900-thin-domain-literature-boundary.md).
Adding Wang's later numerical table cannot fill either pending force authority.

## Consequences for the requested 20D x 10D x 3D case

1. The new mesh must establish its own exact face coordinates, spacing extrema,
   adjacent-cell growth, grid count, decomposition, IBM topology, multigrid
   coarsenability, and fingerprint. Wang Fig. 8 cannot certify any of them.
2. A `10D` transverse width corresponds to about 10% geometric blockage by the
   requested definition, versus Parnaudeau's frozen `4.3%`. The resulting
   domain sensitivity must be reported before any long-statistics release
   decision; agreement or disagreement in `Cd`, `St`, `umin`, or `Lr` cannot be
   attributed to the solver alone.
3. The requested periodic span `3D` is not certified equivalent to
   Parnaudeau's effective experimental span `20D`; Wang does not disclose a
   usable cylinder span/periodic-boundary authority. Midspan profiles remain
   useful, while integrated force statistics require an explicit span caveat.
4. The short initialization and one-/two-step smoke tests requested at this
   node can establish mesh, IBM, conservation, positivity, CFL, solver and
   performance health. They cannot produce statistically converged `Cd`, `St`,
   `umin`, `Lr`, or the 15 Parnaudeau profiles.
5. The later long run, if separately authorized, should compare the complete
   Parnaudeau profile matrix and direct `St` first. Wang Table I should be shown
   in a separate “numerical peer” column with the `umin` inconsistency visible.

## Acceptance boundary

At this node, Wang et al. supports the plausibility of a structured,
near-cylinder-refined Re=3900 IB--LES calculation with a `D/100--D/25` spacing
envelope and the standard three downstream profile stations. It does not
provide enough case metadata or uncertainty to validate the requested mesh,
does not supersede the frozen Parnaudeau data, and cannot by itself justify
publishing HUNDUN-FLOW V1.0. Formal release remains contingent on the separate
long-development/statistics, domain-sensitivity, frozen-literature, software
quality, and GitHub-release gates.
