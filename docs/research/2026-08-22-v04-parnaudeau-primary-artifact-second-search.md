# Parnaudeau et al. (2008) primary-artifact second search

**DOI:** [10.1063/1.2957018](https://doi.org/10.1063/1.2957018)
**Paper:** P. Parnaudeau, J. Carlier, D. Heitz, and E. Lamballais, “Experimental and numerical studies of the flow over a circular cylinder at Reynolds number 3900,” *Physics of Fluids* 20(8), 085101 (2008).
**Search date:** 2026-08-22 (UTC)
**Scope:** second-round search of public publisher, author, institutional, and repository primary artifacts for Figs. 11–15. No login, paywall, anti-bot, or author contact was used. No HUNDUN/COAST/Stage 5 material was read. The only repository file written by this worker is this note.

## Result

One public institutional/author-hosted copy of the paper was acquired from the IRISA Fluminance site. It is a complete 14-page publisher-layout PDF and contains Figs. 11–15, but it is **not established as the AIP publisher VOR**. The source page states that copyright remains with the authors/rightsholders and that the work must not be reposted without permission; the PDF therefore remains only in `/tmp` and is not committed.

The AIP VOR and the exact Crossref-linked AIP PDF were tested but returned Cloudflare HTTP 403 challenge responses. The HAL/INRAE record exposes an archival PDF pointer, but that pointer redirects to the INRA CAS login. ResearchGate was HTTP 403/login restricted. No public author/institution numeric arrays or independent figure files were found in the checked pages. These findings are provenance evidence only: this worker does not decide `DIGITIZATION_SOURCE_ADMISSIBLE`, receipt completeness, tolerance, `full20`, or any final gate.

## Source and access matrix

| Source | Public pointer tested | Result and provenance classification |
|---|---|---|
| Crossref | [Crossref DOI record](https://api.crossref.org/works/10.1063%2F1.2957018) | Primary metadata confirms AIP Publishing, DOI, title/authors, and the exact AIP PDF pointer `https://pubs.aip.org/aip/pof/article-pdf/doi/10.1063/1.2957018/15878439/085101_1_online.pdf`. Metadata pointer only. |
| AIP | [DOI redirect](https://doi.org/10.1063/1.2957018), [article page](https://pubs.aip.org/aip/pof/article/20/8/085101/256405/Experimental-and-numerical-studies-of-the-flow), [PDF URL](https://pubs.aip.org/aip/pof/article-pdf/doi/10.1063/1.2957018/15878439/085101_1_online.pdf) | DOI returned HTTP 302 to AIP, then the article returned HTTP 403 with `cf-mitigated: challenge`; the exact PDF URL returned HTTP 403, `Content-Type: text/html`, not a PDF. No challenge bypass was attempted. VOR not acquired. |
| IRISA/Fluminance | [Dominique Heitz publication list](https://www.irisa.fr/fluminance/team/Dominique.Heitz.html), [publication record](https://www.irisa.fr/fluminance/publi/one_page_per_publi/Parnaudeau08a-eng.html), [institutional PDF](https://www.irisa.fr/fluminance/team/Carlier/publications/ParnaudeauCarlierHeitzLamballais_2008_POF.pdf), [linked PDF](https://www.irisa.fr/fluminance/publi/papers/2008_parnaudeau_etal.pdf) | Both public PDF URLs returned HTTP 200 and the same 1,892,798-byte payload. The page is an author/institution static mirror, not a publisher-VOR claim. This is the acquired candidate artifact. |
| HAL original | [HAL record](https://hal.science/hal-00383669v1), [HAL API record](https://api.archives-ouvertes.fr/search/?q=halId_s:hal-00383669&wt=json&fl=*) | API metadata identifies the DOI and `uri_s=https://hal.science/hal-00383669v1`, but reports `openAccess_bool=false`; no public file field was acquired. Metadata pointer only. |
| HAL/INRAE duplicate | [HAL/INRAE JSON](https://hal.inrae.fr/hal-02590739/json) | Metadata reports `openAccess_bool=false` and contains `seeAlso` pointer `https://archives-publications.inrae.fr/pub00024478.pdf`. This is a metadata pointer, not an acquired file. |
| INRAE archive | [archive PDF pointer](https://archives-publications.inrae.fr/pub00024478.pdf) | HTTP 302 to `https://cas.inra.fr/cas/login?...`; final response was an HTML CAS login page (HTTP 200, 8,298 bytes), not a PDF. Authentication was not attempted. |
| ResearchGate | [public record](https://www.researchgate.net/publication/234130406_Experimental_and_numerical_studies_of_the_flow_over_a_circular_cylinder_at_Reynolds_number_3900) | HTTP 403 Cloudflare response states access is restricted and asks for ResearchGate login/signup. No login or bypass was attempted; no artifact acquired. |
| Author/institution pages | [Pprime Eric Lamballais page](https://pprime.fr/lamballais-eric/), IRISA publication index and Heitz page above | Citation/author metadata only beyond the IRISA PDF. No public arrays or independent Fig. 11–15 files were found. |

## Acquired artifact and provenance receipt

Acquisition command (public static GET, no credentials):

```sh
curl -fsSL --retry 2 --connect-timeout 20 \
  -D /tmp/parnaudeau-irisa-institutional.headers \
  'https://www.irisa.fr/fluminance/team/Carlier/publications/ParnaudeauCarlierHeitzLamballais_2008_POF.pdf' \
  -o /tmp/parnaudeau-irisa-institutional.pdf
```

Recorded response (from `/tmp/parnaudeau-irisa-institutional.headers`):

- UTC `Date`: `Sat, 22 Aug 2026 05:30:24 GMT`
- HTTP `200 OK`; `Content-Type: application/pdf`; `Content-Length: 1892798`
- `Last-Modified: Wed, 03 Mar 2010 16:17:22 GMT`
- `ETag: "1ce1be-480e7d12c8880"`
- `Server: Apache`; `Accept-Ranges: bytes`
- first bytes: `%PDF-1.4\r\n%\xe2\xe3\xcf\xd3\r\n`
- `file`: `PDF document, version 1.4`; MIME: `application/pdf`

The second IRISA URL returned HTTP 200 with the same byte count and PDF hash (`Last-Modified: Fri, 27 Nov 2009 11:28:11 GMT`, `ETag: "1ce1be-47958963548c0"`). The PDF first page identifies the DOI, article, four authors, and AIP copyright notice. `pdfinfo` reports:

```text
Creator:        XPP
Producer:       Acrobat Distiller 6.0.1 (Windows)
CreationDate:   Fri Aug  1 00:35:08 2008 CST
ModDate:        Fri Aug  1 00:35:08 2008 CST
Pages:          14
Encrypted:      no
Page size:      612 x 792 pts (letter)
File size:      1892798 bytes
PDF version:    1.4
```

Local evidence and SHA256 (all under `/tmp`; no PDF is in the repository):

| Local file | Meaning | SHA256 |
|---|---|---|
| `/tmp/parnaudeau-irisa-institutional.pdf` | acquired PDF, 1,892,798 bytes | `b8a775e5a5078e19fc47d9c5f47e95b81b4a68d4449e4444459088ed9befcdd4` |
| `/tmp/parnaudeau-irisa-institutional.headers` | acquisition headers, 461 bytes | `c386e13dc9a0c98a27d365c04805c509327b3105360d731a1f78e3b2338934f5` |
| `/tmp/parnaudeau-irisa-institutional.txt` | `pdftotext -layout` extraction | `67d011221b5916ffafe4757c3071e70ca0513ee93929fd8303ee1bc26dba0110` |
| `/tmp/parnaudeau-irisa-heitz.html` | author/institution page | `ac98493680be2f9e154c7f795edb9d4cfaf8d16d2a206bb2d54e813d354216e6` |
| `/tmp/parnaudeau-irisa-publication-page.html` | publication record and copyright notice | `2ea975c3ccb194868131c7ade8cc63c48fe0b2be027b4772cf75ed38ce07fc8f` |
| `/tmp/parnaudeau-irisa-index-publications.html` | IRISA publication index | `16acc52f43703a76898ae7117e8826d6d05bac4fb300d9facb80cb3b5dbce6ca` |
| `/tmp/parnaudeau-pprime-lamballais.html` | Pprime author page | `d415812aaadbc26a17f1c9c84b68e2a08b46311b24460a001ffbebaa302eafaa` |

The IRISA publication page links the same PDF and explicitly presents a copyright notice retaining rights with the authors/rightsholders. This note links to the public source but does not redistribute the PDF.

## Fig. 11–15 location and controlled-inspection facts

Text extraction and page scan commands:

```sh
pdftotext -layout /tmp/parnaudeau-irisa-institutional.pdf \
  /tmp/parnaudeau-irisa-institutional.txt
for p in $(seq 1 14); do
  pdftotext -f "$p" -l "$p" -layout /tmp/parnaudeau-irisa-institutional.pdf - |
    rg -q 'FIG\. 1[1-5]|Figure 1[1-5]' && echo "page=$p"
done
```

The scan returned `page=11` and `page=12`. The extracted captions and station labels establish:

| PDF page | Figures | Observable and stations |
|---|---|---|
| 11 | Fig. 11 and Fig. 12 | mean streamwise velocity and mean normal velocity; each has `x/D = 1.06`, `1.54`, `2.02` profiles versus `y/D` |
| 12 | Fig. 13 and Fig. 15 | streamwise fluctuation variance and velocity-fluctuation covariance; each has the same three `x/D` stations versus `y/D` |
| 12 | Fig. 14 | normal fluctuation variance; the same three `x/D` stations versus `y/D` |

All five captions refer to the Fig. 9 caption for the series legend. The extracted Fig. 9 legend identifies the solid line as present LES, filled circles as “present PIV experiments (cases 1 and 2),” open circles as present HWA, and other symbols/line styles as prior experiments or simulations. The body text then states that the statistical profiles presented thereafter come from **case 2**, with a small `0.025D` field/window versus `0.045D` for case 1. Therefore:

- case 2 is identified textually as the source of the statistical profiles in Figs. 11–15;
- the plotted filled-circle legend marker is shared by “present PIV experiments (cases 1 and 2)” and is **not a case-2-only marker**; the marker glyph alone cannot separate cases 1 and 2;
- at the rendered inspection scale, every one of the three profiles in each Fig. 11–15 panel has visibly more than eight plotted marker locations for the present PIV series (and the comparison `x` series is also multi-point). This is a visual sufficiency observation, not a recovered numeric array or an exact point count;
- each composite profile panel has a shared `y/D` horizontal coordinate with labeled ticks from approximately `-2` to `2` and a labeled ordinate with substantially more than three ticks. The three profiles are vertically offset within one panel; they are not three separately boxed axes;
- `pdfimages -f 11 -l 12 -list /tmp/parnaudeau-irisa-institutional.pdf` reports no embedded raster image objects for pages 11–12. The figures and text are vector PDF content with embedded fonts, so the source has no finite native pixel resolution;
- for a reproducible raster inspection only, `pdftoppm -f 11 -l 12 -r 300 -png ...` produced 300-dpi, `2550 x 3301` PNGs. These are derived review images, not source artifacts.

Derived inspection image hashes:

| Local file | Resolution | SHA256 |
|---|---:|---|
| `/tmp/parnaudeau-irisa-page-11.png` | 2550 x 3301 (300 dpi) | `02a49ba3850da0f18a8750b7ccfc25ddfc17bfad4d784c26b16babcfbd7f200c` |
| `/tmp/parnaudeau-irisa-page-12.png` | 2550 x 3301 (300 dpi) | `ad55df6203d0969d9d6305d372834c4bc7486aaa203f57486690312b6fc13a79` |

## Numeric-array and independent-figure search

The checked IRISA author/institution pages, IRISA publication index, and Pprime author page were scanned for direct links ending in `.dat`, `.txt`, `.csv`, `.zip`, `.tar`, `.h5`, `.hdf5`, or `.mat`; each of the four page scans returned no match. The only relevant IRISA publication link was the PDF and a `.ref` citation file. No public author/institution arrays, point tables, or independent Fig. 11–15 raster/vector files were located in this round. This is a bounded search result, not a claim that no copy can exist elsewhere.

## Exact blockers and negative evidence

### AIP VOR

Crossref returned the exact AIP PDF URL, but the recorded response was HTTP 403, `Content-Type: text/html; charset=UTF-8`, with `cf-mitigated: challenge`. The DOI redirect was HTTP 302 followed by the same AIP Cloudflare 403. The attempted AIP body was not a PDF:

- `/tmp/parnaudeau-aip-publisher-pdf.body`: HTML, SHA256 `95927d2a1cb44498607a2c7f59a0a43525b050ba042351b3c3dbd0889f311942`
- `/tmp/parnaudeau-aip-publisher-pdf.headers`: SHA256 `fb83bb2e1f9ed551fbc73d5542964f10bd392e55399317d1fcad1b19fe92881d`

### HAL/INRAE archive

The HAL/INRAE `seeAlso` pointer redirects to CAS login. The recorded final response is HTML, not a PDF; no authentication or bypass was attempted.

### ResearchGate

The public record returned HTTP 403 and an “Access restricted”/unusual-activity page requesting login or signup. No login or bypass was attempted.

### Arrays and figure assets

No numeric arrays, point tables, or independent figure files were exposed by the checked author/institution pages. The acquired PDF contains curves and marker glyphs only; it does not itself provide pointwise numeric tables or PIV uncertainty arrays.

## Reproduction commands and log hashes

The following commands were used for the public metadata/artifact checks (URLs are also recorded in the matrix above):

```sh
curl -sSIL --max-time 30 --retry 1 'https://doi.org/10.1063/1.2957018'
curl -sSIL --max-time 30 --retry 1 'https://pubs.aip.org/pof/article/20/8/085101/256405/Experimental-and-numerical-studies-of-the-flow'
curl -sSIL --max-time 30 --retry 1 'https://pubs.aip.org/aip/pof/article-pdf/doi/10.1063/1.2957018/15878439/085101_1_online.pdf'
curl -fsSL 'https://api.crossref.org/works/10.1063%2F1.2957018'
curl -fsSL 'https://api.archives-ouvertes.fr/search/?q=halId_s:hal-00383669&wt=json&fl=*'
curl -fsSL 'https://hal.inrae.fr/hal-02590739/json'
curl -sSIL --max-time 30 --retry 1 'https://archives-publications.inrae.fr/pub00024478.pdf'
curl -sSIL --max-time 30 --retry 1 'https://www.researchgate.net/publication/234130406_Experimental_and_numerical_studies_of_the_flow_over_a_circular_cylinder_at_Reynolds_number_3900'
curl -fsSL 'https://www.irisa.fr/fluminance/team/Dominique.Heitz.html'
curl -fsSL 'https://www.irisa.fr/fluminance/publi/one_page_per_publi/Parnaudeau08a-eng.html'
curl -fsSL 'https://pprime.fr/lamballais-eric/'
```

The existing complete response logs and SHA256 values are:

| Log | SHA256 |
|---|---|
| `/tmp/parnaudeau-primary-second-search-aip-blocker.log` | `2efaebcdd02d3455763e7a13b8f24a35143c8be3485f54d16773c6636bdcb227` |
| `/tmp/parnaudeau-primary-second-search-author-institution-pages.log` | `63b50895180a7b7e34a41b0d09d1a11d1d615b59b66dac5e990a741d1de89d95` |
| `/tmp/parnaudeau-primary-second-search-crossref-hal.log` | `1258251f4fd342c0d3da4b6d4be173143b21abfe55689f3c9194053c04046be1` |
| `/tmp/parnaudeau-primary-second-search-doi-redirect.log` | `c4d441417e2ce1ebee8cdb313420be338af3c547104eb325af27d13055d361ba` |
| `/tmp/parnaudeau-primary-second-search-hal.log` | `21ffa09ab05b8ffd9ba49592c1b5972e2b2410d233cee5022de3f627182d8afe` |
| `/tmp/parnaudeau-primary-second-search-heads.log` | `4d62051b3665f4fd710592ea9bd754c738f62f112f52d7bb4dda93f45b3ecfea` |
| `/tmp/parnaudeau-primary-second-search-inrae-archive-head.log` | `0c232923b29ed3364c5ab79c3d481661186b9b167016c68d675d31c79a638f16` |
| `/tmp/parnaudeau-primary-second-search-irisa-linked-corrected-head.log` | `62b91624e40b647225b37d9e859f2797292aaa92a22c14b7c7e29cfa77a834d1` |
| `/tmp/parnaudeau-primary-second-search-irisa-linked-head.log` | `baa3993e82a11c95e25674d92ed5c79ea71abf43f277771afb676258b1797d87` |
| `/tmp/parnaudeau-primary-second-search-irisa-pages.log` | `0d12c6959a0829d1693494b6e0ac43930929906838eec055cb1aff47caa2d878` |
| `/tmp/parnaudeau-primary-second-search-irisa.log` | `ed00d7d38e43a8abb0c0753bf0383f187cd88d09380dcc0cce70a0fc4df40624` |
| `/tmp/parnaudeau-primary-second-search-publisher-hal-json.log` | `73851f70c57174f22ee67bbfa6615878d594d88b84cc1dc67484f7cff0f13d7a` |
| `/tmp/parnaudeau-primary-second-search-researchgate.log` | `2b49492745852f48aa123cae2fd1e9866cc4bb31d0c6c5b3f25683905ae331f5` |
| `/tmp/parnaudeau-pdftoppm.log` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` (empty successful rasterization log) |

The AnySearch runtime used for discovery was `/usr/bin/python3` with CLI `/home/wyf/.codex/skills/anysearch/scripts/anysearch_cli.py` (CLI SHA256 `912c8ea187c5cbdb7966a67ab1450523171ce99902de269f1a52a4520d9584e2`; Python 3.6.9). Discovery results were treated as URL pointers and then verified by direct public GET/HEAD responses above; no search-result page was treated as an acquired artifact.

## Worker handoff conclusion

The second search produced one publicly accessible, provenance-checkable institutional/author-hosted paper PDF with Figs. 11–15 and vector source content. It does not produce author arrays or a case-2-specific marker encoding: case 2 is stated in the article text, while the legend combines cases 1 and 2 under one filled-circle symbol. AIP VOR access, HAL/INRAE archive access, and ResearchGate access remain blocked by the exact responses recorded above. The parent agent must independently review whether the IRISA copy meets the project’s admissibility and controlled-digitization rules; this worker makes no such determination.
