<!-- SPDX-License-Identifier: Apache-2.0 -->

# Third-party software

## yyjson

- Version: 0.12.0
- Repository: https://github.com/ibireme/yyjson
- Tag: 0.12.0
- Annotated tag object: `7871d321ff4cd8068c1f777c97975dc2fb640ab3`
- Peeled commit: `8b4a38dc994a110abaec8a400615567bd996105f`
- License: MIT
- Vendored files: `third_party/yyjson/yyjson.c`, `third_party/yyjson/yyjson.h`
- License text: `LICENSES/yyjson-MIT.txt`

## Design references (not vendored or linked)

The following projects were inspected only for public mathematics, data-layout concepts,
resource-lifecycle ideas, or replacement capability. They are not copied into, linked by, or
required at runtime by HUNDUN-FLOW. Exact evidence paths and adoption decisions are frozen in
`docs/references/2026-08-13-hundun-v04-public-lifecycle-survey.md` and its TSV ledger.

### OpenFOAM-dev

- Repository: https://github.com/OpenFOAM/OpenFOAM-dev.git
- Commit: `b9da51ab0673423aa2af6a45a72a3fbec9c66f9f`
- Upstream license: GPL-3.0-or-later (`COPYING` at the fixed commit)
- Reference scope: public PISO mathematics and lifecycle of `rAU`, `rAtU`, `HbyA`,
  `phiHbyA`, pressure-equation flux, and final velocity update.
- Boundary: GPL-reference-only. HUNDUN independently derives and implements its equations,
  data structures, schedule, and tests; no OpenFOAM source or translation is included.

### AMReX

- Repository: https://github.com/AMReX-Codes/amrex.git
- Commit: `59d066aab774bc388cc6ed944f7beaf645607ed3`
- Upstream license: BSD 3-Clause terms (`LICENSE` at the fixed commit)
- Reference scope: box-local layout, cached halo/fill-patch metadata, begin/finish exchange,
  and EB factory/resource lifetime.
- Boundary: independently reimplemented HUNDUN ideas only. No AMReX source, API, or runtime
  type is included; AMR fill-patch/reflux/average-down are outside v0.4.

### IncFlo

- Repository: https://github.com/AMReX-Fluids/incflo.git
- Commit: `7307d8725c2a538f09cafbeacbfeb63e0fb11d22`
- Upstream license: BSD-3-Clause (`LICENSE` at the fixed commit)
- Reference scope: projection, EB, and operator/projector reconstruction boundaries.
- Boundary: independently reimplemented HUNDUN ideas only. No IncFlo source, API, or AMR
  regrid implementation is included.

### AMReX-Hydro

- Repository: https://github.com/AMReX-Fluids/amrex-hydro.git
- Commit: `e49df248aabd2cc11865eb5be734a2f5f2f65ee5`
- Upstream license: BSD 3-Clause terms (`LICENSE` at the fixed commit)
- Reference scope: projector, multigrid operator, coefficient refresh, solver, and workspace
  reuse.
- Boundary: independently reimplemented behind HUNDUN-owned linear interfaces. No upstream
  source, object model, or runtime dependency is included.

### COAST functional reference

- Read-only location: `/home/wyf/code_dev/Coast_software`
- Observation date: 2026-08-13
- Reference scope: replacement capabilities and observed lifecycle only, including local
  absolute-pressure EOS, pressure density derivative, sound speed/Mach, NSCBC, SIMPLE-style
  transient iterations, Vreman/wall-function paths, flat `.d` inputs, STL scan, Restart,
  Visit, screen, and ICCG.
- Boundary: user-specified read-only reference; no general redistribution or license grant is
  presumed. Task 1 copies no COAST source. Old COAST Fortran may not be copied or translated,
  and no other COAST C++ reuse is authorized. Its SIMPLE loop and lack of periodic product
  capability are reference limitations, not HUNDUN v0.4 algorithm requirements.
- Sole exception: for Task 5, the user authorizes porting or adapting the `imb_mesh_y.cpp`
  scanning mathematics/method without provenance boilerplate. This narrow authorization does
  not extend to any other COAST C++ or to old COAST Fortran.

### Cantera example-data scientific reference

- Locator: Cantera 3.1.0 points `data/example_data` to the independent
  `Cantera/cantera-example-data` gitlink
  `1a5d27e508a38b1791543e9fded80ffd5c5b8d75`; the referenced record is
  `n-hexane-NUIG-2015.yaml` at that commit.
- Reference scope: the record's O2 and N2 NASA7 values were used as an
  independently documented scientific source for the fixed 0.21/0.79
  mole-fraction pseudo-air coefficients.
- Licensing boundary: that pinned example-data repository does not provide a
  license grant for the record, and the parent Cantera repository's license is
  not asserted to cover it. HUNDUN neither redistributes the YAML/mechanism nor
  links Cantera at runtime; it records only the directly mole-weighted O2/N2
  numerical scientific values and makes no upstream-license claim for them.

### Published transport-model mathematics

- P. Yoon and G. Thodos, "Viscosity of Nonpolar Gaseous Mixtures at Normal
  Pressures," AIChE Journal 16(2), 300--304 (1970), DOI
  `10.1002/aic.690160225`.
- C. R. Wilke, "A Viscosity Equation for Gas Mixtures," Journal of Chemical
  Physics 18 (1950), DOI `10.1063/1.1747673`.
- A. W. Vreman, "An eddy-viscosity subgrid-scale model for turbulent shear
  flow: Algebraic theory and applications," Physics of Fluids 16 (2004), DOI
  `10.1063/1.1785131`.
- Reference scope: public mathematical definitions for the Yoon--Thodos
  low-pressure species viscosity, Wilke gas-mixture viscosity weighting, and
  the Vreman SGS invariant. The fixed O2/N2 inputs use `Tc` in K, `Pc` in atm,
  and molecular weight in g/mol; the correlation returns micropoise before the
  exact `1e-7 Pa s` conversion. The retained O2/N2 inputs
  `Tc=154.4/126.2 K` and `Pc=49.7/33.5 atm` are the deliberately rounded
  comparison-contract values also tabulated in R. B. Bird, W. E. Stewart, and
  E. N. Lightfoot, *Transport Phenomena*, 2nd ed. (Wiley, 2002), Appendix E,
  Table E.1. NIST Chemistry WebBook SRD 69 was used only as an independent
  reasonableness check, not as the source of those exact rounded values.
  HUNDUN owns the implementation and tests; no publication or third-party
  source text is copied.
