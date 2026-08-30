# Tadrist et al. (1990) primary-PDF inventory

审计日期：2026-08-22（Asia/Shanghai）。状态：`RESEARCH_ONLY / PRIMARY-PDF
INVENTORY`。

本文件只登记用户授权提供的完整 primary PDF 中可逐页核验的 source facts。它不作
`piD` scientific-work equivalence、不选择生产数值、不设 tolerance、不更新
literature receipt 或 candidate ledger，也不作 `ACCEPT/REJECT` 判断。

## Source identity and byte inventory

论文：H. Tadrist, R. Martin, L. Tadrist, and P. Seguin, “Experimental investigation
of fluctuating forces exerted on a cylindrical tube (Reynolds numbers from 3000 to
30000),” *Physics of Fluids A: Fluid Dynamics* **2**, 2176--2182 (1990), DOI
[`10.1063/1.857804`](https://doi.org/10.1063/1.857804)。PDF 第 1 页是 AIP
front/copyright page；正文是 PDF 第 2--8 页，对应印刷页 2176--2182。

用户授权的 PDF 原始 bytes 位于仓库外：

~~~text
/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/
  tadrist-1990-user-authorized-e62f0c6b/tadrist1990.pdf
~~~

| field | observed value |
|---|---|
| bytes | 543,517 |
| SHA-256 | `e62f0c6bfc2b2bea7e89de6ffb075eee0740bcc1ee393254e44c2bb84ab4e1d3` |
| PDF version / producer | 1.4 / iText 2.1.7 by 1T3XT |
| pages | 8 |
| page size | 577.44 × 810.96 pt |
| encrypted | no |
| structure check | `qpdf --check`: no syntax or stream encoding errors reported |
| rendered review | 8 page JPEGs plus contact sheet already sealed beside the PDF; their hashes are in `SHA256SUMS.txt` |

The PDF contains the AIP copyright/download footer. That footer is provenance metadata,
not an experimental result.

## Page-by-page source inventory

The locators below use both the PDF page number and the printed journal page. All eight
pages were extracted with `pdftotext -layout` and rendered; the page-1--8 contact sheet
was reviewed, with high-resolution checks of the equations/apparatus on PDF pages 2--5.

### PDF p. 1 (no printed article page)

* AIP article landing/copyright page. It gives the title, four authors, journal citation,
  DOI, and related-article links.
* It contains no experimental apparatus, coefficient definition, target-Re data point,
  calibration result, uncertainty, or mean-drag result.

### PDF p. 2 (printed p. 2176)

* The abstract states that the experiment measures a **global lift coefficient** over
  Reynolds numbers 3000--30,000, with upstream velocity and turbulence measured by
  laser-Doppler anemometry. It says aspect ratio, blockage, and end effects are
  discussed, and that global lift varies substantially over Re 3000--10,000.
* For a cylindrical section at axial position `z`, the paper defines the per-unit-length
  pressure-force components (Eqs. 1a--1b) as

  ~~~text
  F_d(z,t) = (D/2) ∫_0^(2π) p(z,θ,t) cos(θ) dθ
  F_l(z,t) = (D/2) ∫_0^(2π) p(z,θ,t) sin(θ) dθ
  ~~~

  Here `p(z,θ,t)` is pressure; the text explicitly calls these forces per unit
  length. The local coefficients (Eqs. 2a--2b) are
  `C_d = F_d / (0.5 ρ U² D)` and
  `C_l = F_l / (0.5 ρ U² D)`.
* The paper introduces force intercorrelation spectra for two axial positions (z_1,z_2),
  and describes the fluctuating loads as random steady-state space/time functions.
  No numerical uncertainty or mean-drag result is given on this page.

### PDF p. 3 (printed p. 2177)

* Eq. 5 defines force variance from the force power spectrum. Eqs. 6--7 use an exponential
  axial coherence model and define a correlation length `λ_c`. Eqs. 10--12 relate the
  global lift coefficient to the local coefficient through a double integral over the
  cylinder length `L`:

  ~~~text
  C_z² = C_l² / L² ∫_0^L ∫_0^L exp(-|z₁-z₂|/λ_c) dz₁ dz₂
  ~~~

  followed by the paper's generalized-correlation-length expression. The paper then
  states that, allowing for three-dimensional shedding, the global lift coefficient was
  determined by **direct measurement of the fluctuating forces exerted on the cylinder**.
* The suction wind tunnel operated at atmospheric pressure. The rectangular test section
  is 300 mm high × 200 mm wide × 400 mm deep. LDA found a flat upstream profile in a
  section perpendicular to the main flow.
* For upstream velocities above 4 m/s, turbulence is defined as
  `Tu = σ_u/U` and is reported as velocity-independent and equal to 5%.
* The force sensor is a ceramic piezoelectric detector between two undeformable bodies;
  the mechanical media is stiffer than the ceramic media and the ceramic is preloaded.
  Sensors were dynamically calibrated with a Brüel & Kjær load amplifier, vibrating pot,
  force sensor, and accelerometer (Fig. 2). The next page continues the calibration
  principle. No calibration residual, uncertainty, repeatability, or error bar is reported.

### PDF p. 4 (printed p. 2178)

* The calibration force was applied at the middle of a tube with known mechanical
  properties; vibration frequency and amplitude were monitored with an accelerometer,
  and applied force was inferred from acceleration and tube mechanical properties.
* Measurement procedure: a rigid cylinder restrained at **both ends** was placed vertically
  at the test-section center. It had **two force sensors, one at each end**; the assembly
  was supported by an arm isolated from the tunnel. The insertion opening was blocked with
  a flexible foam plug.
* Four rigid Pyrex tubes had bending natural frequencies above 200 Hz. The paper lists,
  in order, aspect ratios `L/D = 30, 18.75, 16, 12.5` and cross-sectional blockage ratios
  `D/l_u = 5%, 8%, 9%, 12%`, respectively. The paper does not print the physical length
  `L` in millimetres, exposed/effective span, end-plate dimensions, or a separate
  end-condition drawing.
* A Hewlett Packard 3562A analyzer acquired and processed the piezoelectric-balance output.
  The authors averaged about 20 spectra over a 500-Hz full-scale range with 0.625-Hz
  resolution. The vortex-shedding relation is `f_s = S U / D` (Eq. 14).
* The two sensor signals were reported as statistically symmetric about the cylinder
  midplane; the interspectrum modulus was unitary and phase null over the analysis band.
  Fig. 3 shows force power spectra for a 24-mm cylinder at (U=7.2) and 20 m/s; Fig. 4
  plots Strouhal number versus Re. Neither figure is a Re=3900 force-coefficient table.

### PDF p. 5 (printed p. 2179)

* Section III.3 is titled “Standard deviation of the lift coefficient.” Eq. 15 defines
  the fluctuating-force variance as

  ~~~text
  σ_F² = ∫_(-∞)^(+∞) S_F(f) df
  ~~~

  `S_F(f)` is stated to be the fluctuating-force power spectral density exerted over
  the **entire length of the cylinder**. For the digitized signal, the integral is a sum
  over the frequency band containing the aerodynamic force; the entire characteristic
  vortex-shedding band was included. The paper reports about 7% difference between using
  the full band and the midheight band. This 7% is only a band-choice sensitivity in the
  force-spectrum integration; the paper does not identify it as measurement uncertainty
  and it must not be used as an error bar. The paper does not provide a sampling time,
  confidence interval, repeatability, or uncertainty budget.
* Measurements used four cylinders with diameters 10, 16, 18, and 24 mm over Re 3000--30,000.
  The text says that over Re 3000--7000 the lift coefficient rises sharply from about 0.02
  to 0.08; over Re 8000--10,000 it rises slightly from about 0.08 to 0.2; above Re 10,000
  it depends on tube diameter. The 24-mm and 16-mm results differ by roughly 50% in
  10,000 < Re < 15,000, attributed to aspect ratio, blockage, end effects, and turbulence.
* Fig. 5 is a 24-mm-cylinder two-sensor intercorrelation example at (U=13.5) m/s; it
  reports unit modulus and zero phase across the band.
* Fig. 6 is the global-lift plot with vertical axis `C_z`, logarithmic Re axis from
  `10³` to `10⁵`, and marker legend `D=10, 16, 18, 24` mm. Its caption says
  `D=10, 14, 18, 24` mm, which conflicts with both the apparatus text and the legend's
  `D=16` label. This is a source-internal caption typo and is retained here verbatim as
  an audit fact.
* There is no table and no marker labelled exactly Re=3900. The only source-readable
  target-neighborhood value is the prose envelope `C_z` (the paper's plotted lift
  coefficient), approximately 0.02--0.08 over Re 3000--7000. A point-specific value at
  Re≈3900 is not reproducibly printed: the graph has no point labels, its markers are
  crowded, and the source gives no tabulated coordinates. No digitized or interpolated
  value is recorded in this inventory.

### PDF p. 6 (printed p. 2180)

* Fig. 7 compares the force spectrum for a 16-mm cylinder at (U=10) m/s with a 32-mm
  tunnel-wall clearance opening unsealed versus sealed. The caption says the unsealed
  opening inhibits vortex shedding and the sealed opening restores the shedding band.
  The text states that flexible foam plugs were used to minimize this end effect.
* The text says global coefficient depends on aspect ratio and that, for a given correlation
  length, the generalized correlation length changes with `L/D`. It uses `λ_c = 4D`,
  attributed to Surry and Sonneville, to **deduce local** lift-coefficient variation from
  the global measurement (Fig. 8). The 4D value is an adopted correlation-length model
  input, not a correlation length measured by this paper. That local curve is not an
  additional direct sensor measurement.
* The authors state that blockage and wall effects were small relative to aspect-ratio
  effects for the 10-, 16-, and 18-mm tubes, while blockage above 10% changes longitudinal
  correlation length; high blockage can also alter shedding. The 24-mm tube is the listed
  12% blockage case, but the paper does not give a numerical blockage correction.
* End-effect discussion cites other studies: 15--20D at Re<300 and approximately 6D in a
  cited Re=45,000 case with `L/D = 12.5`. These are not a measured end-plate/span value for
  the present experiment. The present apparatus still has no printed end-plate geometry,
  exposed span, or effective-force-length definition.
* Fig. 8's caption also prints (D=10,14,18,24) mm while the apparatus and legend use
  (D=16) mm. The plot is explicitly “local” and is deduced from global measurements.

### PDF p. 7 (printed p. 2181)

* The discussion compares the present results with other authors. The cited Bishop--Hassan
  setup is described as a free-surface water stream, (L/D=5), blockage 23%, and Re 4000--
  10,000; those are **other-study** facts, not Tadrist apparatus facts.
* For the present experiment, the authors describe a lift-coefficient increase over
  3000<Re<10,000 and a slight plateau near Re≈4000. The conclusion calls Re 4000--10,000
  a transition zone and says the global lift coefficient was determined with the developed
  experimental unit. No exact Re=3900 coordinate, uncertainty, or table is supplied.
* The page discusses drag only by comparison with literature (including Whitaker); it does
  not report a Tadrist mean drag measurement or a Tadrist mean-drag coefficient. The paper's
  earlier (C_d) equation is a local pressure-force definition, not a reported mean-drag
  result.

### PDF p. 8 (printed p. 2182)

* Acknowledgments and the 26-item reference list only. It adds no apparatus dimensions,
  force coefficients, Re=3900 value, calibration uncertainty, or mean drag.
* The references identify Bishop and Hassan, Gerrard, Norberg, Surry, Sonneville, and other
  works used in the discussion. Those references are not substituted for the Tadrist
  primary measurements in this inventory.

## Requested-field matrix

| requested field | source fact in this PDF | explicit gap retained |
|---|---|---|
| experimental apparatus | Atmospheric-pressure suction wind tunnel; 300×200×400 mm rectangular section; LDA upstream velocity/profile/turbulence; vertical rigid cylinder at center | No full dimensional drawing of the force assembly or end restraints |
| cylinder diameter | 10, 16, 18, 24 mm in apparatus/results text and Fig. 6 legend | Fig. 6 and Fig. 8 captions print 14 mm in place of 16 mm |
| physical length / aspect ratio | (L/D=30,18.75,16,12.5), respectively, for four tubes | No printed mm lengths; exact diameter-to-aspect pairing is only given by the listed order |
| effective/active span | The force PSD is described over the entire cylinder length `L` | No exposed span, effective span, force-collecting length, or `πD` normalization statement |
| ends / end plates | Cylinder restrained at both ends; insertion opening sealed by flexible foam plug; foam used to reduce end effect | No end-plate dimensions, dummy-end geometry, leakage area, or end correction |
| blockage | Cross-sectional `D/l_u = 5%, 8%, 9%, 12%`, listed in the same order as aspect ratios | No correction equation or uncertainty for blockage; wall/test-section interpretation is qualitative |
| Re coverage | Overall Re 3000--30,000; text bands 3000--7000, 8000--10,000, and >10,000 | No exact source coordinate at Re=3900 and no per-run Re table |
| global lift sensor | Two force sensors, one at each end; piezoelectric balance; arm isolated from tunnel; direct fluctuating-force measurement | No explicit algebra for combining the two signals into the plotted `C_z` |
| coefficient definition | Local pressure-integral `C_d, C_l` in Eqs. 1--2; global relation `C_z²/C_l²` in Eq. 10 | The paper does not explicitly print a modern `C_L,rms` normalization line for Fig. 6 |
| RMS / standard deviation | Section III.3 and Eq. 15 derive force variance from the entire-band PSD; ~7% full-band/midheight-band band-choice sensitivity | The 7% is not a reported measurement uncertainty; no sampling window, degrees of freedom, confidence interval, or error bar |
| calibration | Dynamic calibration with Brüel & Kjær load amplifier, vibrating pot, force sensor, accelerometer; force inferred from acceleration and tube mechanics | No calibration curve, residual, drift, sensitivity, or numerical uncertainty |
| target Re≈3900 plot/table value | Fig. 6 has no labelled/table point; prose gives only `C_z ≈ 0.02--0.08` over Re 3000--7000 and a near-4000 plateau | No reproducible exact `C_z(3900)`; no digitization/interpolation made |
| mean drag | `C_d` is introduced as a local pressure-force definition; discussion cites other drag work | No Tadrist mean-drag measurement, mean-drag coefficient, or drag data table/figure |

## Controlled interpretation boundary

The PDF supports the narrow source statement “a two-sensor piezoelectric apparatus directly
measured fluctuating lift over the full tube length, with results plotted as a global lift
coefficient over Re 3000--30,000.” It does not, by itself, close the following fields:

1. a reproducible finite-span/effective-span definition suitable for a `πD` target;
2. a point value at Re=3900;
3. a modern `C_L,rms` normalization and sampling/uncertainty protocol; or
4. a Tadrist mean-drag result.

The prose value range and the plotted `C_z` are retained as source facts only. No inference
from the paper's correlation model, no conversion to a different force normalization, and no
comparison to HUNDUN output is performed here.

## Reproducibility, complete diff, and validation

The byte and page artifacts were checked with the following commands:

~~~sh
pdf=/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/tadrist-1990-user-authorized-e62f0c6b/tadrist1990.pdf
sha256sum "$pdf"
pdfinfo "$pdf"
qpdf --check "$pdf"
pdftotext -layout "$pdf" /tmp/tadrist1990-layout.txt
render_dir=$(mktemp -d /tmp/tadrist1990-pages.XXXXXX)
pdftoppm -png -r 220 "$pdf" "$render_dir/page"
sha256sum -c /home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/tadrist-1990-user-authorized-e62f0c6b/SHA256SUMS.txt
~~~

Observed controls:

* SHA-256 exactly matched `e62f0c6bfc2b2bea7e89de6ffb075eee0740bcc1ee393254e44c2bb84ab4e1d3`.
* `pdfinfo` reported 8 pages, 543,517 bytes, PDF 1.4, unencrypted.
* `qpdf --check` reported no syntax or stream-encoding errors.
* `pdftotext -layout` produced 538 lines; its SHA-256 was
  `e7dc806ccbf8beffdace3ef76868a9ce3a06496fb14cba37483e1003257cbf38`, matching the
  pre-sealed `tadrist1990-layout.txt` beside the PDF.
* `sha256sum -c` returned success for the PDF, all eight rendered page files, the contact
  sheet, and the extracted layout text in the artifact directory's `SHA256SUMS.txt`.

The complete repository diff for this subtask is exactly one new file, this file, and can be
reproduced without changing any other path:

~~~sh
git diff --no-index -- /dev/null \
  docs/research/2026-08-22-v04-tadrist1990-primary-pdf-inventory.md
git diff --check -- \
  docs/research/2026-08-22-v04-tadrist1990-primary-pdf-inventory.md
git status --short --untracked-files=all
~~~

No receipt, ledger, reference, solver, COAST, Stage5, or other repository path was modified;
no commit was created by this worker.
