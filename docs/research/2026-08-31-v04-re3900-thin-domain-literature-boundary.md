<!-- SPDX-License-Identifier: Apache-2.0 -->

# V0.4 Re=3900 thin-domain comparison: frozen literature boundary

Date: 2026-08-31
State: `SCOPED BENCHMARK AUTHORITY / NOT A V0.4 RELEASE GATE`

This note fixes the literature and definition boundary for the quick/medium-length
HUNDUN--COAST comparison in a periodic thin cylinder domain. It introduces no new
literature value. It uses only the frozen repository references, their extraction
rule, and the R4 partial receipt. In particular, it does not turn an unavailable
total drag or a non-equivalent lift statistic into an acceptance target.

## Controlling decision

| Result class | Frozen quantities | Permitted use in this benchmark |
|---|---|---|
| Paper-comparison authority | The 15 Parnaudeau PIV profiles at three stations; Parnaudeau HWA `St=0.208 +/- 0.002` | Scoped paper comparison, with the exact normalization and extraction/experimental error status below |
| Paper advisory | Parnaudeau `Lr/D=1.51`, `Umin/Uc=-0.34`, formation length `0.87`; Norberg bracketing/direct and formula-derived Strouhal values; Norberg sectional lift fit; Tadrist global-lift neighborhood | Report separately; no hard pass/fail band may be invented where the source reports no uncertainty or exact target point |
| Solver-to-solver diagnostic | `Cd`, periodic-`pi D` integrated `Cl_rms`, force histories, instantaneous/mean fields and spectra not covered by a frozen paper observable | Compare HUNDUN with COAST under identical definitions; never label the comparison a paper validation |
| Non-blocking unavailable authority | Direct total `Cd` at exactly Re=3900; lift RMS demonstrably equivalent to a periodic `pi D` active span | Must not block this scoped comparison and must not be fabricated from pressure drag, sectional lift, a broad prose envelope, or a later CFD table |

The machine boundary remains
[`docs/verification/v0.4-literature-data-receipt-r4-partial.json`](../verification/v0.4-literature-data-receipt-r4-partial.json):
`complete=false`, with only the Norberg reference incomplete. Validation without
`--require-complete` succeeds; validation with `--require-complete` rejects. That
is the correct state for this explicitly scoped benchmark and does not authorize
the V0.4 `literature` or `final` release gates.

## Parnaudeau experiment and coordinate identity

The frozen experiment fields in
[`docs/references/cylinder-re3900-parnaudeau.json`](../references/cylinder-re3900-parnaudeau.json)
are `Re=3900`, `D=0.012 m`, `Uc=4.8 m/s`, blockage `4.3%`, effective
span `20D`, and PIV plane `z/D=0` (`experiment.*`). The primary-source audit
records the coordinate origin at the cylinder centre, `x` along the wind-tunnel
axis, `z` along the cylinder axis, and right-handed `y`; it also fixes
`Re=Uc D/nu` and `St=f_vs D/Uc`
([`docs/research/2026-08-21-v04-cylinder3900-backstep-primary-data.md`](2026-08-21-v04-cylinder3900-backstep-primary-data.md),
section 2.1).

For the comparison:

- the profile coordinate is `y/D`, at the paper stations `x/D=1.06, 1.54,
  2.02`; no nearest-station substitution or streamwise interpolation is an
  authority;
- first moments are normalized by `Uc`; second moments are normalized by
  `Uc^2`;
- a plot may say `U/Uinf` only when the run manifest explicitly binds
  `Uinf=Uc`. Otherwise the paper label remains `Uc` and the solver output is
  converted using that declared reference;
- an axis-orientation conversion, if needed, must be declared once before
  comparison and applied consistently to `y`, mean `v`, and `u'v'`. A sign flip
  chosen after seeing the curves is inadmissible;
- the thin computational span `pi D` is not geometrically equivalent to the
  experiment's effective `20D` span. The midspan profiles remain admissible
  observables, but a systematic mismatch may be classified as a thin-span
  limitation rather than hidden or used to claim full-span validation.

The exact source is the institution-hosted 14-page primary PDF, SHA-256
`b8a775e5a5078e19fc47d9c5f47e95b81b4a68d4449e4444459088ed9befcdd4`.
The accepted target is the filled-circle present-PIV series; article scope fixes
those statistical profiles as PIV case 2. The source and series identity are
frozen in
[`docs/verification/v0.4-parnaudeau-source-layout-freeze.md`](../verification/v0.4-parnaudeau-source-layout-freeze.md),
sections "Source decision" and "Target-series identity".

## Accepted 15-profile matrix

Every array below is embedded under `profile_authority.profiles` in the
Parnaudeau JSON with status `controlled_digitization_complete`. The first data
column is `y/D`; the second is the quantity shown. `Ey` and `Eq` are absolute,
profile-level controlled-digitization bounds for the coordinate and value,
respectively. They are not PIV experimental confidence intervals.

| PDF page / figure | `x/D` | Frozen quantity | Points | `Ey` | `Eq` |
|---|---:|---|---:|---:|---:|
| 11 / Fig. 11 | 1.06 | mean `u/Uc` | 130 | 0.012578726241965817 | 0.0125810536008402 |
| 11 / Fig. 11 | 1.54 | mean `u/Uc` | 101 | 0.012578726241965817 | 0.0125810536008402 |
| 11 / Fig. 11 | 2.02 | mean `u/Uc` | 158 | 0.012578726241965817 | 0.0125810536008402 |
| 11 / Fig. 12 | 1.06 | mean `v/Uc` | 158 | 0.01257872802513601 | 0.009501425666658002 |
| 11 / Fig. 12 | 1.54 | mean `v/Uc` | 157 | 0.01257872802513601 | 0.009501425666658002 |
| 11 / Fig. 12 | 2.02 | mean `v/Uc` | 158 | 0.01257872802513601 | 0.009501425666658002 |
| 12 / Fig. 13 | 1.06 | `u'u'/Uc^2` | 128 | 0.012580626629830642 | 0.00285489287736669 |
| 12 / Fig. 13 | 1.54 | `u'u'/Uc^2` | 91 | 0.012580626629830642 | 0.00285489287736669 |
| 12 / Fig. 13 | 2.02 | `u'u'/Uc^2` | 118 | 0.012580626629830642 | 0.00285489287736669 |
| 12 / Fig. 14 | 1.06 | `v'v'/Uc^2` | 158 | 0.012585738784388262 | 0.0038028904236447098 |
| 12 / Fig. 14 | 1.54 | `v'v'/Uc^2` | 139 | 0.012585738784388262 | 0.0038028869582144302 |
| 12 / Fig. 14 | 2.02 | `v'v'/Uc^2` | 97 | 0.012540618983676834 | 0.0038028904236447098 |
| 12 / Fig. 15 | 1.06 | `u'v'/Uc^2` | 152 | 0.01258062662983153 | 0.0022131179871697233 |
| 12 / Fig. 15 | 1.54 | `u'v'/Uc^2` | 111 | 0.01258062662983153 | 0.0022131159523322203 |
| 12 / Fig. 15 | 2.02 | `u'v'/Uc^2` | 102 | 0.01258062662983153 | 0.0022131179871697233 |

The source locators are primary PDF pages 11--12, Figs. 11--15, as stored in
each profile's `source_locator`. The preregistered common transverse interval is

```text
-0.9117483947994405 <= y/D <= 0.9391120663474339
```

(`profile_authority.common_y_over_d_interval`). All profile points should be
shown in the plots. A scalar aggregate across all 15 profiles may use the
common interval because it was frozen before the solver results; it must say so
and must not selectively trim individual panels.

### Meaning of the profile error bounds

The exact-vector procedure is frozen in
[`docs/verification/v0.4-parnaudeau-controlled-digitization-rule.md`](../verification/v0.4-parnaudeau-controlled-digitization-rule.md),
especially "Panel calibration and sampling" and "Exact-vector selection and
extraction-error bound". For axis `z` and cleared extraction `k`, it uses

```text
B_z,k = r_z,k + |a_z,k| (q_k + w_k/2)
E_z   = max(B_z,1, B_z,2)
```

where `r` is the actual-tick affine-fit residual, `q` is half the converted
coordinate's decimal quantization step, and `w/2` is half the printed marker
width. The composition is conservative L1. Both conversions have exactly equal
canonical marker inventories; no smoothing, densification, symmetry fill, or
station interpolation is present.

Consequently:

- paper markers should carry horizontal `+/- Ey` and vertical `+/- Eq` bars or
  an equivalent rectangle/band;
- solver interpolation error and temporal/block sampling confidence intervals
  remain separate from `Ey/Eq` and should be reported alongside them;
- the paper reports no pointwise experimental uncertainty separately for these
  PIV curves (`experimental_uncertainty_status=not_reported_separately`);
- HWA first-/second-moment percentages and the paper's finite-duration LES
  estimates must not be reassigned to the PIV points;
- no single percentage tolerance is authorized. Report per-profile errors with
  the paper extraction bounds and the solver sampling uncertainty, and expose
  the full residual curve rather than reducing the decision to a favorable
  subrange.

The extractor enforces the exact 3-by-5 matrix, ordered `y/D`, at least eight
points per profile, common coverage, locators, method fields, and separate
coordinate/value bounds
([`tools/v04_literature_extract.py`](../../tools/v04_literature_extract.py),
`validate_parnaudeau_profiles`). The authoritative R3 replay and independent
audit are recorded in
[`docs/verification/2026-08-22-v04-parnaudeau-r2-independent-audit-resolution.md`](../verification/2026-08-22-v04-parnaudeau-r2-independent-audit-resolution.md).

## Scalar and force authority classification

### Parnaudeau

The following fields are in `observables` in the Parnaudeau JSON:

- `strouhal_number = 0.208`, definition `f_vs D/Uc`, with reported absolute
  uncertainty `0.002`: admissible as the benchmark's paper Strouhal comparison;
- `recirculation_length_over_diameter = 1.51`, Table II PIV;
- `centerline_minimum_mean_u_over_uc = -0.34`, Table II PIV;
- `streamwise_variance_formation_length_over_diameter = 0.87`, Table II PIV.

The last three have `uncertainty.status=not_reported`; they may be reported as
paper advisory diagnostics but do not acquire an invented pass band. The
recirculation and formation distances are measured from the cylinder trailing
edge, while the profile-station coordinate origin remains the cylinder centre.

### Norberg

[`docs/references/cylinder-re3900-norberg.json`](../references/cylinder-re3900-norberg.json)
freezes:

- direct `St=0.2108` at `Re=3704.6` and `St=0.2104` at `Re=4211.1` from the
  1987 appendix table; neither row reports measurement uncertainty, and the
  `-3 dB` bandwidth is explicitly not uncertainty;
- `St=0.208884` at Re=3900, derived from the Norberg 2003 Table 4 fit; its fit
  interval/uncertainty is not reported;
- sectional `Cl_rms=0.083046` at Re=3900, also derived from the Table 4 fit.

The Strouhal rows are useful corroboration but are advisory rather than a second
hard paper gate. The lift value is sectional, pressure-derived in the audited
Norberg lineage, and is not a finite-span total coefficient. It cannot be
compared as though it were periodic-`pi D` integrated `Cl_rms`. The controlling
primary-source audit is
[`docs/research/2026-08-22-v04-norberg2003-finite-span-lift-primary-audit.md`](2026-08-22-v04-norberg2003-finite-span-lift-primary-audit.md):
the Norberg present lift row begins at `Re=10^4`, and no signed spanwise
covariance or `pi D` active-span measurement closes the mapping.

Norberg's `required_observables` therefore remain
`total_mean_drag_coefficient_re3900` and
`finite_span_lift_rms_coefficient_re3900`, both represented exactly in
`pending_observables`. The exclusions in `excluded_as_release_oracles` are
controlling: the later-paper `Cd about 0.98 +/- 0.05`, `Cl_rms about 0.10 +/-
0.05`, and the 1987 `C_Dp=0.98` at Re=3000 are not admissible replacements.

### Tadrist

[`docs/references/cylinder-re3900-tadrist.json`](../references/cylinder-re3900-tadrist.json)
and the page-by-page
[`docs/research/2026-08-22-v04-tadrist1990-primary-pdf-inventory.md`](2026-08-22-v04-tadrist1990-primary-pdf-inventory.md)
support only this scoped interpretation:

- printed pp. 2176--2179 establish a direct two-end-sensor global fluctuating-
  lift measurement over Re=3000--30000;
- printed p. 2179 and Fig. 6 give only the approximate prose envelope
  `C_z=0.02--0.08` over Re=3000--7000, not a labelled Re=3900 point;
- printed p. 2179 reports a 7% full-band versus midheight-band difference; this
  is a spectrum-band sensitivity, not experimental uncertainty;
- printed p. 2180 adopts `lambda_c=4D` from other studies to derive local lift;
  it is not a Tadrist measurement and does not establish `pi D` equivalence;
- the apparatus has `L/D=12.5--30`, blockage `5--12%`, upstream `Tu=5%`, and
  end restraints. It supplies neither an exact active-span mapping to periodic
  `pi D` nor a Tadrist total mean-drag result.

Thus the broad global-lift range is a visibly labelled target-neighborhood
advisory observation only. It is not a pointwise uncertainty band and cannot
accept or reject periodic-`pi D` `Cl_rms`. HUNDUN--COAST `Cl_rms` and `Cd`
remain solver-to-solver diagnostics for this benchmark.

## Comparison and verdict boundary

The post-processing may issue one of the task's three scoped verdicts only
after this ordering is respected:

1. Establish complete, numerically healthy HUNDUN and COAST sampling windows.
2. Compare all 15 Parnaudeau profiles with their frozen coordinate/value
   digitization bounds and each solver's independently estimated sampling
   uncertainty. Preserve per-station and per-quantity results.
3. Compare shedding Strouhal with Parnaudeau's direct `0.208 +/- 0.002`; list
   Norberg results separately as advisory corroboration.
4. Report `Cd`, periodic-`pi D` `Cl_rms`, field differences, and any quantities
   lacking a frozen paper uncertainty as solver-to-solver diagnostics or paper
   advisory observations, as classified above.
5. If a profile mismatch persists, classify evidence in the order requested by
   the benchmark: sampling insufficiency, thin-span bias, boundary/IBM
   difference, or solver difference. Do not widen the literature band or choose
   a favorable time subwindow.

Neither a missing direct total `Cd` nor a missing periodic-`pi D`-equivalent
`Cl_rms` is a rejection condition for this thin-domain comparison. Conversely,
agreement of the two solvers on those quantities cannot be promoted to paper
validation or V0.4 release acceptance.

## Frozen hashes and reproducibility

The R4 receipt binds these exact repository artifacts:

| Artifact | SHA-256 |
|---|---|
| `docs/references/cylinder-re3900-parnaudeau.json` | `8a90682cf6d5d503538dabf47f8097a7e18f962eeeb5e53ce39858bcbb6361b9` |
| `docs/references/cylinder-re3900-norberg.json` | `4c8db4c241c72b383d68828d391e30b9f721128e4d1b0f83758b63df3ce6696f` |
| `docs/references/cylinder-re3900-tadrist.json` | `6788b727e72558dd87fd19b41e63e0030cc570ad772dc0a5e4cc5b2468e1705b` |
| `tools/v04_literature_extract.py` | `7336f5dc1c697ea514b9b8a3011a14e3ad2a2465c9957bb10e7b96dfd714ead6` |
| R4 partial receipt | `1dfdeb16b92f4a29a5dc2a0c4daa146499e8decb1af5df197b652203e36b2e3f` |

Recheck the scoped authority with:

```sh
python3 tools/v04_literature_extract.py validate \
  docs/references/cylinder-re3900-parnaudeau.json \
  docs/references/cylinder-re3900-norberg.json \
  docs/references/cylinder-re3900-tadrist.json
python3 tools/v04_literature_extract.py receipt-validate \
  docs/verification/v0.4-literature-data-receipt-r4-partial.json
```

Do not add `--require-complete` for this scoped comparison; that mode is the
separate release-authority check and correctly fails while the two Norberg
force rows remain pending.
