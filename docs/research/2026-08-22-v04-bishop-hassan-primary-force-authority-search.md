# Bishop--Hassan (1964) Re=3900 force-authority search

检索日期：2026-08-22（UTC）。

本笔记是 `RESEARCH_ONLY` 的 primary-source legwork。它只记录 Bishop & Hassan
的公开身份元数据、合法获取尝试和目前可审计的事实；不修改 literature receipt、
candidate、ledger、plan 或 solver，不读取 HUNDUN long-stat output，也不作
`ACCEPT/REJECT`、等价性或容差决定。

## 结论摘要

目标文献是 R. E. D. Bishop and A. Y. Hassan, *The lift and drag forces on a
circular cylinder in a flowing fluid*, Proceedings of the Royal Society of London
A 277(1368), 32--50 (1964), DOI
[10.1098/rspa.1964.0004](https://doi.org/10.1098/rspa.1964.0004)。Crossref 的一手
出版元数据和摘要确认：作者描述了直接测量 fluctuating lift、fluctuating drag
以及 steady mean drag 的 apparatus；圆柱轴线横向于水道流向；stationary-cylinder
结果覆盖 Re=3600--11000。Re=3900 因而落在文献声明的区间内。

但是，本次没有合法取得可读的完整 primary paper/scan。Royal Society 的 DOI
landing page、DOI PDF 和 Crossref 提供的 syndication PDF URL 均在本环境返回
Cloudflare 403 challenge；OpenAlex、Unpaywall 和 Semantic Scholar 都把该工作
标为 closed、没有 repository PDF。没有完整纸本文本就不能逐页核对 force balance、
active span、coefficient normalization、end plates、free-surface/blockage
correction、Re=3900 的具体点、统计窗口或 uncertainty。因此目前只能确认
“直接测力 apparatus + stationary Re range”的方法级事实，不能把它登记为
`Cd_mean` 或 `Cl_rms` 的数值 authority。

这也意味着下列问题保持未决，而本报告不替 root 作等价性判断：

| 问题 | 当前 primary evidence | 状态 |
|---|---|---|
| (a) direct total mean drag `Cd` | 摘要只确认直接测量 steady mean drag **force**；没有系数定义、归一化、数值表/图或 measured force segment 的可读 primary 页面 | 未核实 |
| (b) direct fluctuating lift RMS coefficient | 摘要确认直接测量 fluctuating lift **force**；没有 RMS/statistical definition、coefficient normalization 或数值数组 | 未核实 |
| measured active span | 未取得 apparatus/force-balance 页面；active segment、dummy ends、end plates 和 spanwise support 未核实 | 未核实 |
| boundary/free surface | 摘要只说 water channel；水深、submergence、free-surface criterion/correction 未核实 | 未核实 |
| blockage | channel width/depth、blockage ratio 和 correction 未核实 | 未核实 |
| Re=3900 value | 只确认声明的 range 3600--11000；没有 Re=3900 具体点或可控图表 digitization | 未核实 |
| uncertainty/repeatability | 未取得 measurement-accuracy/error section | 未核实 |
| HUNDUN periodic span `πD` equivalence | active span、boundary/end treatment、normalization、spanwise correlation 都缺 primary evidence；本笔记不作等价性结论 | 未评定 |

## 允许使用的 primary metadata

### Publisher/DOI identity

Crossref 的 DOI record（该 record 的 `resource.primary.URL` 指向 Royal Society
article）给出以下 identity：

* title: *The lift and drag forces on a circular cylinder in a flowing fluid*;
* authors: Richard Evelyn Donohue Bishop and A. Y. Hassan，均标为 Department of
  Mechanical Engineering, University College London；
* published: 1964-01-07；volume 277, issue 1368, pages 32--50；
* publisher: The Royal Society；
* publisher PDF links:
  `https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.1964.0004` and
  `https://royalsocietypublishing.org/rspa/article-pdf/277/1368/32/54549/rspa.1964.0004.pdf`;
* Crossref abstract states direct fluctuating lift/drag and steady mean drag
  measurement on a cylinder normal to the water-channel flow, with stationary
  results for Re 3600--11000.

Sources:

* [Crossref DOI record](https://api.crossref.org/works/10.1098/rspa.1964.0004)
  (publisher-deposited metadata and abstract; not a paper scan).
* [Royal Society article landing URL](https://royalsocietypublishing.org/rspa/article/277/1368/32/11602/The-lift-and-drag-forces-on-a-circular-cylinder-in)
  (publisher primary landing page; inaccessible from this retrieval environment).
* [DOI resolver](https://doi.org/10.1098/rspa.1964.0004).

The [McMaster institutional scholarly-work record](https://experts.mcmaster.ca/scholarly-works/3519678)
repeats the same article identity and abstract and exposes only a DOI link; it does
not expose a primary file. It is used here only as an institutional metadata
cross-check, not as a numerical source.

### What the abstract does and does not establish

The abstract is enough to establish these narrow facts:

1. the experiment used a water channel and a cylinder whose central axis was
   perpendicular to the flow direction;
2. the apparatus was intended to measure fluctuating lift and drag directly and
   steady mean drag directly;
3. the reported stationary-cylinder Reynolds-number interval includes Re=3900.

It does **not** establish any of the following: whether the reported mean force is
the whole exposed cylinder or an active finite segment; whether dummy ends/end
plates were used; how the force balance separates fluid force from support/inertia;
whether a mean drag coefficient is tabulated; whether fluctuating lift is reported
as RMS, amplitude, spectral component, or another statistic; the reference area and
velocity used for normalization; sampling duration and filtering; uncertainty;
blockage/free-surface corrections; or the spanwise correlation needed to compare a
finite measured segment with HUNDUN's periodic span πD.

## Retrieval and access evidence

No target artifact directory was created because no primary PDF/scan was obtained.
Before searching, the requested directory was checked and did not exist:

```sh
find /home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature \
  -maxdepth 1 -type d -name 'bishop-hassan-1964-primary-search*' -print
# no output
```

The following records were captured under `/tmp/bishop-search-20260822/` during
the search. They are retrieval logs only and are not repository artifacts.

| Probe | HTTP/content result | bytes | body SHA-256 |
|---|---|---:|---|
| Crossref DOI JSON | 200, `application/json` | 3962 | `e63905e1ddf27afa1ec9bf541acb7951ede1c71ed3123f481b91a1c449caf78b` |
| OpenAlex DOI record | 200, JSON; `is_oa=false`, `any_repository_has_fulltext=false` | 13080 | `330fd40c95ff5938cd1b9b4817c12c68c1ea8db527931ac147449220595dfd67` |
| Unpaywall DOI record | 200, JSON; `is_oa=false`, `oa_locations=[]`, `has_repository_copy=false` | 720 | `6154e5508f21232a56b6fdfc0fba0a432dc8c523ccf0a46df62cd52dd535ea32` |
| Semantic Scholar DOI record | 200, JSON; `openAccessPdf.status=CLOSED`, empty PDF URL | 508 | `81360d1ec5e2f821a1d18bcb18625d95c06667c6b2dddb041058107aa0b36384` |
| McMaster institutional record | 200, HTML; metadata/abstract only | 85683 | `0a52f49d0d5e685f4a6082212b909459798a306d80a90be14d882088eb14e5f0` |
| Royal Society landing URL | 403, `cf-mitigated: challenge`, HTML | 5852 | `a0ba5ea22a690c824e9982e27a21f865f8f9c97d5f039365c87b137003861fdb` |
| Royal Society DOI PDF URL | 403, `cf-mitigated: challenge`, HTML | 5635 | `2762010a6ba79c43ace48ddd52cfcc578bd003e5e593578f3240aa6f91bd8c24` |
| Royal Society syndication PDF URL | 403, `cf-mitigated: challenge`, HTML | 5725 | `0089a8ae1947c64d63686eb36b1b8a24167ccc8fa992ce5c5f81041c5f5e721f` |

The Crossref record supplies `54549` as the current publisher syndication asset
identifier; trying the older guessed `54536` path was also a 403 and did not return
a PDF. A `r.jina.ai` text proxy returned the Royal Society Cloudflare challenge
text, not article content, so it was not retained as an artifact.

Internet Archive advanced-search probes for the exact title, author identifiers,
and DOI timed out with no bytes (`HTTP=000`, `curl` exit 28). This is a network
access failure, not evidence that the archive contains no copy. Searches of UCL
Discovery/EThOS and general DOI/title indexes found metadata or citations but no
author/institution primary scan. No external author request was sent.

## Reproducible search commands

The core retrieval probes can be reproduced as follows. `-D` captures response
headers; `-o` captures the response body; the output must be placed in a new
temporary directory and must never overwrite a candidate or receipt artifact.

```sh
set -o pipefail
mkdir -p /tmp/bishop-search-20260822

curl -L --max-time 30 \
  -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D /tmp/bishop-search-20260822/crossref.headers \
  -o /tmp/bishop-search-20260822/crossref.body \
  'https://api.crossref.org/works/10.1098/rspa.1964.0004'

curl -L --max-time 30 \
  -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D /tmp/bishop-search-20260822/royal-landing.headers \
  -o /tmp/bishop-search-20260822/royal-landing.body \
  'https://royalsocietypublishing.org/rspa/article/277/1368/32/11602/The-lift-and-drag-forces-on-a-circular-cylinder-in'

curl -L --max-time 30 \
  -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D /tmp/bishop-search-20260822/royal-pdf.headers \
  -o /tmp/bishop-search-20260822/royal-pdf.body \
  'https://royalsocietypublishing.org/rspa/article-pdf/277/1368/32/54549/rspa.1964.0004.pdf'

curl -L --max-time 30 \
  -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D /tmp/bishop-search-20260822/openalex.headers \
  -o /tmp/bishop-search-20260822/openalex.body \
  'https://api.openalex.org/works/https://doi.org/10.1098/rspa.1964.0004'

curl -L --max-time 30 \
  -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D /tmp/bishop-search-20260822/unpaywall.headers \
  -o /tmp/bishop-search-20260822/unpaywall.body \
  'https://api.unpaywall.org/v2/10.1098/rspa.1964.0004?email=research@example.org'

curl -L --max-time 30 \
  -A 'mailto:research@example.org' \
  -D /tmp/bishop-search-20260822/semanticscholar.headers \
  -o /tmp/bishop-search-20260822/semanticscholar.body \
  'https://api.semanticscholar.org/graph/v1/paper/DOI:10.1098/rspa.1964.0004?fields=paperId,title,year,authors,openAccessPdf,url,externalIds'

sha256sum /tmp/bishop-search-20260822/*.body
```

## Evidence matrix for a future page-level audit

The matrix below is intentionally a gap register, not an inferred reconstruction.
Each row needs the primary article page/figure/table and a stable artifact hash before
it can be used as a numerical authority.

| Required field | Needed primary evidence | Current result |
|---|---|---|
| stationary condition | explicit fixed-cylinder setup and run protocol | abstract confirms “stationary cylinder” results, setup details unavailable |
| direct force observable | force-balance diagram, transducer arrangement, sign convention, separation of support/inertial terms | abstract confirms direct fluctuating lift/drag and steady mean drag measurement only |
| `Cd_mean` | equation, dynamic-pressure/reference-area/span normalization, actual Re=3900 value/table/figure | missing |
| `Cl_rms` | time-series/statistical definition, window/filter, RMS vs amplitude, normalization and actual value | missing |
| active span | cylinder length, active measuring segment, dummy ends, end plates, force-transducer attachment | missing |
| spanwise physics | evidence for 2-D suppression, end effects, correlation/phase over active segment | missing |
| boundaries | channel width/depth, free surface/submergence, inlet/outlet and side-wall condition | missing |
| blockage | ratio definition and raw/corrected values | missing |
| Re mapping | `U`, `D`, viscosity/temperature, exact sampled Re values and interpolation rule | only broad 3600--11000 range is known |
| uncertainty | calibration, sensitivity, repeatability, scatter, digitization/reading error | missing |
| HUNDUN comparison | a declared mapping from measured finite span to periodic span πD, including normalization and spanwise correlation | not assessed here |

## Controlling gaps and lawful next retrieval

The controlling gap is access to the primary pages, not a lack of evidence that the
authors intended direct force measurements. A future continuation should obtain a
licensed library scan, publisher-accessible PDF, or author/institution deposit and
first freeze its bytes, MIME type, PDF metadata, page count and SHA-256 in a new
external directory such as:

```text
/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/
  bishop-hassan-1964-primary-search-20260822/
```

The directory did not exist at the start of this search and must not be created by
copying any secondary PDF. Once a genuine primary artifact is obtained, the page audit
order should be: apparatus/force balance; active and dummy span; end plates and
boundary corrections; coefficient equations; Re=3900 table/figure; time-series/RMS
definition; calibration and uncertainty. Only those page-level facts can determine
whether the measured work is suitable for a separate HUNDUN periodic-span comparison.

## Validation log

* The requested new report path did not exist before this write.
* No existing repository file was edited.
* No primary PDF/scan was obtained or copied.
* No solver, COAST case, long-stat output, receipt, candidate, ledger or Stage 5
  content was read or changed.
* The external retrieval response bodies and headers listed above are temporary
  search evidence; their identities are recorded so a later agent can reproduce the
  access result without treating the 403 HTML as a paper artifact.
