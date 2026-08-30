# Re=3900 primary-authority gap audit

Date: 2026-08-21

## Decision

The literature authority is still incomplete. No checked primary source makes
the Parnaudeau et al. (2008) Fig. 11--15 PIV arrays or a direct Re=3900 total
force record available as a small machine-readable authority. The 2026 INRAE
PIV deposit is primary and public, but its published file inventory does not
match the catalogue description closely enough to substitute it without a
separate reconstruction and experiment-mapping receipt.

Accordingly, `docs/verification/v0.4-literature-data-receipt.json` must remain
`complete=false`. This audit did not inspect HUNDUN literature outputs and did
not start long solver statistics.

## Parnaudeau et al. (2008)

Primary article: P. Parnaudeau, J. Carlier, D. Heitz and E. Lamballais,
*Experimental and numerical studies of the flow over a circular cylinder at
Reynolds number 3900*, Physics of Fluids 20, 085101 (2008),
[DOI 10.1063/1.2957018](https://doi.org/10.1063/1.2957018).

- Figures 11--15 publish the three stations `x/D=1.06,1.54,2.02` as curves,
  not numeric tables.
- Footnote 31 says that the experimental and numerical statistics are
  available by contacting the authors. This is evidence that the article
  itself does not provide a public array attachment.
- The author-uploaded full text on
  [ResearchGate](https://www.researchgate.net/publication/234130406_Experimental_and_numerical_studies_of_the_flow_over_a_circular_cylinder_at_Reynolds_number_3900)
  is an admissible source for controlled figure digitization, but digitization
  would still have to bind the exact PDF, crop/calibration data, separated PIV
  curve picks, symmetry checks and extraction error.
- No author/DOI/hosting-page attachment exposing the Fig. 11--15 arrays was
  found in this audit. The lawful next routes are author-provided arrays or a
  controlled digitization receipt; absence of an array is not permission to
  invent points or reuse LES curves.

## 2026 INRAE PIV deposit

Primary dataset: P. Georgeault and D. Heitz,
*Non-time-resolved PIV dataset of flow over a circular cylinder at Reynolds
number 3900*, Recherche Data Gouv V1,
[DOI 10.57745/DHJXM6](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/DHJXM6),
Etalab Open License 2.0.

The catalogue says the deposit contains 10,000 instantaneous fields plus mean,
Reynolds-stress, skewness and flatness fields. The public Dataverse V1 metadata
API was queried directly at:

```text
https://entrepot.recherche.data.gouv.fr/api/datasets/:persistentId/?persistentId=doi:10.57745/DHJXM6
```

The returned inventory instead contains 998 public files, all named
`Serie_10001.txt` through `Serie_10999.txt` with `Serie_10902.txt` absent. Their
aggregate size is 13,176,725,066 bytes. No separately named mean or
Reynolds-stress file, directory label, category or file description is present.
A range request to the first file, Dataverse data-file id `720393`, begins:

```text
#DaVis 10.x 2C vector field 4 545 740 "x [mm]";"y [mm]";"Vx [m/s]";"Vy [m/s]"
```

Thus the exposed files are instantaneous 2D2C fields, not ready-made profile
arrays. The deposit can become a primary authority only after all of the
following are frozen:

1. a complete API inventory and content-hash receipt;
2. an explanation of the 998-file/10,000-field and missing-statistics mismatch;
3. cylinder origin, diameter, reference velocity and experimental-batch
   reconciliation with the 2008 case;
4. a deterministic mean/Reynolds-stress reconstruction over the actual
   published ensemble;
5. interpolation at the registered stations with sampling/interpolation
   uncertainty and an explicit statement that these are 2026 dataset profiles,
   not relabelled 2008 Fig. 11--15 arrays.

Until those steps are complete, this dataset is useful primary evidence but is
not a drop-in replacement for the pending profile authority.

## Force authority

The checked primary Norberg material does not supply a machine-readable direct
total-force record exactly at Re=3900:

- Norberg (1987) tabulates Strouhal points around Re=3900, while its directly
  tabulated pressure-drag cases are Re=3000 and Re=8000; pressure drag is not
  total drag.
- Norberg (2003) provides a sectional fluctuating-lift correlation. Its value
  evaluated at Re=3900 is a fitted, local sectional quantity, not a direct
  finite-span total-force measurement.
- Values subsequently summarized by CFD papers are secondary transcriptions
  unless traced to a primary table/array or controlled primary-figure
  digitization.

No primary machine-readable `Cd_mean`/finite-span `Cl_rms` authority exactly at
Re=3900 was found. Such a force oracle therefore remains optional only if the
final physical-accuracy policy is explicitly based on the direct Strouhal,
recirculation and registered velocity/stress profiles; it must not be presented
as a completed direct-force comparison.

## Commands used

Read-only source checks used the Dataverse metadata API, a byte-range request
to data-file id `720393`, the DOI/author full-text pages, and the existing
repository references. No external message was sent to an author, no dataset
bundle was downloaded, and no solver or pairing run was launched.
