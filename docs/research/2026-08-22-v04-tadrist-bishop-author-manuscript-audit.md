# Re=3900 force-source author/repository manuscript audit

检索日期：2026-08-22（UTC）。本报告是 RESEARCH_ONLY 的一次性 primary-source access
audit，只回答两个 DOI 是否有作者或机构 repository 中合法公开、可读的完整
manuscript/PDF：

* R. E. D. Bishop and A. Y. Hassan, The lift and drag forces on a circular cylinder in
  a flowing fluid, Proceedings of the Royal Society of London A 277(1368), 32--50
  (1964), DOI https://doi.org/10.1098/rspa.1964.0004。
* H. Tadrist, R. Martin, L. Tadrist and P. Seguin, Experimental investigation of
  fluctuating forces exerted on a cylindrical tube (Reynolds numbers from 3000 to
  30 000), Physics of Fluids A 2(12), 2176--2182 (1990), DOI
  https://doi.org/10.1063/1.857804。

不读取 HUNDUN 输出，不提取或决定任何生产数值，不修改 literature receipt、candidate
ledger、reference、solver、COAST 或 Stage 5，也不绕过付费墙、登录、DRM 或 JavaScript
challenge。普通搜索结果和二手论文未被当作 authority，也未写入数值结论。

## 结论

本轮没有取得任一目标的完整 primary paper、作者 manuscript 或可读 PDF。两个目标都不
应因此升级为 literature authority，且没有创建 benchmark literature artifact directory。

* Bishop--Hassan：Royal Society DOI landing、DOI PDF 和 publisher syndication PDF
  均返回 Cloudflare 403 challenge；JSTOR 原刊稳定条目的 PDF 入口返回 3038-byte
  Client Challenge HTML。UCL Discovery 的首页、DOI/title 搜索和 OAI-PMH 入口均
  返回 Cloudflare 403 challenge。
* Tadrist et al.：AIP DOI landing 和 Crossref 提供的 PDF URL 均返回 Cloudflare
  403 challenge；旧 Scitation PDF/EPDF URL 返回 200 HTML 壳而非 PDF。HAL 的 DOI
  和精确题名查询均为 numFound=0；AMU HAL 页面、AMU institution repository 和
  IUSTI 历史站点分别受到 bot challenge、TLS 连接失败或 gateway 失败。
* OpenAlex、Unpaywall、Semantic Scholar 和 OpenAIRE 的机器记录未给出可读
  author/repository file URL；DataCite 对两个 DOI 均返回 404。

控制性 blocker 仍是合法取得一份完整原始页面。后续若由用户或机构提供有授权的 PDF，
应先把原始 bytes、MIME、页数、pdfinfo、qpdf --check 和 SHA-256 封存，再做
apparatus/span/normalization/uncertainty 审计。本轮没有 PDF，所以没有执行 pdfinfo
或 qpdf，也没有伪造其结果。

## DOI identity and registry records

Crossref 是两篇期刊记录的出版方沉积元数据来源；它确认题名、作者、卷期页码及
publisher PDF URL，但不是论文正文。

| record | URL | Date / type | elapsed | bytes | body SHA-256 | result |
|---|---|---|---:|---:|---|---|
| Bishop Crossref | https://api.crossref.org/works/10.1098/rspa.1964.0004 | 09:31:53Z / JSON | 1.091436 s | 3962 | e63905e1ddf27afa1ec9bf541acb7951ede1c71ed3123f481b91a1c449caf78b | identity + publisher URLs |
| Tadrist Crossref | https://api.crossref.org/works/10.1063/1.857804 | 09:29:35Z / JSON | 2.113915 s | 7422 | 8996ee96c75d0520ecccae924a144685d26f02cafa8198fe12b70279b678f22e | identity + AIP PDF URL |

Crossref first-party PDF links:

~~~text
https://pubs.aip.org/pof/article-pdf/2/12/2176/12512312/2176_1_online.pdf
https://royalsocietypublishing.org/rspa/article-pdf/277/1368/32/54549/rspa.1964.0004.pdf
~~~

Independent registry/aggregator checks were used only as availability corroboration:

| record | status/type, Date | elapsed | bytes | body SHA-256 | machine result |
|---|---|---:|---:|---|---|
| Bishop OpenAlex | 200 JSON, 09:31:57Z | 0.946813 s | 13080 | 2b37ec40621240859d52d87e123356f8c05355e0040bfe196f96204ab0070bd4 | is_oa=false; no repository full text/PDF URL |
| Tadrist OpenAlex | 200 JSON, 09:29:38Z | 1.054791 s | 14417 | 63bf2b0b1856c60408e0f5a368e0c473ccd1d612bd5207f73527815baf5e8e31 | is_oa=false; no repository full text/PDF URL |
| Bishop Unpaywall | 200 JSON, 09:31:58Z | 1.266219 s | 720 | 6154e5508f21232a56b6fdfc0fba0a432dc8c523ccf0a46df62cd52dd535ea32 | oa_status=closed; no locations/repository copy |
| Tadrist Unpaywall | 200 JSON, 09:29:40Z | 1.614107 s | 2024 | 286c653b20576c71ca52abed9225dff8081b30711842bd14ec2ce6187fcfe327 | oa_status=closed; no locations/repository copy |
| Bishop Semantic Scholar | 200 JSON, 09:32:00Z | 1.408123 s | 508 | 81360d1ec5e2f821a1d18bcb18625d95c06667c6b2dddb041058107aa0b36384 | openAccessPdf CLOSED; empty URL |
| Tadrist Semantic Scholar | 200 JSON, 09:29:41Z | 1.017944 s | 648 | 894af5a1fd2d01f8e2be804ea7ec34d304ed0ed06be1f3f9555e6604747a676b | openAccessPdf CLOSED; empty URL |
| Bishop DataCite | 404 JSON, 09:39:54Z | 1.772658 s | 87 | 30b443bb0e590a597609b71b83a2424d39b6f050b18d234c6543cdd1f8791b6b | no DataCite record |
| Tadrist DataCite | 404 JSON, 09:39:56Z | 2.706931 s | 87 | 30b443bb0e590a597609b71b83a2424d39b6f050b18d234c6543cdd1f8791b6b | no DataCite record |
| Bishop OpenAIRE | 200 XML, 09:39:58Z | 1.976291 s | 13390 | 0c9dc3aa8ca047a1e5a3ce1395d416d8a6287db305b9f0e31e979b91c4a24448 | metadata only, CLOSED, DOI URL |
| Tadrist OpenAIRE | 200 XML, 09:40:00Z | 1.950211 s | 14803 | fd6f7c729277e2df1aeed2f9e8674887eb5cd02ace49375e364c2110789c4a57 | metadata only, UNKNOWN, DOI URL |

这些机器记录不能证明私有馆藏扫描不存在，只证明审计时被查询记录没有宣传 repository
file。

## Primary and first-party paths

### Bishop--Hassan: publisher and journal archive

| path | Date / status / headers | elapsed | bytes | response-body SHA-256 | result |
|---|---|---:|---:|---:|---|
| DOI resolver https://doi.org/10.1098/rspa.1964.0004 | 09:31:55Z / 403 HTML; cf-mitigated: challenge | 2.835018 s | 5852 | 0482edb8eb6cbcfed9c53b7711bd6aefbb8e5a404c5bda95769c02b66e1dcaf7 | Royal Society redirect; no paper |
| Royal Society landing https://royalsocietypublishing.org/rspa/article/277/1368/32/11602/The-lift-and-drag-forces-on-a-circular-cylinder-in | 09:36:24Z / 403 HTML; cf-mitigated: challenge | 0.973547 s | 5852 | 9fff55e3f5237355bcb571fe00ce812ed58b8c6addd29fc43d7af7be7c4a4b60 | no landing page |
| Royal Society DOI PDF https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.1964.0004 | 09:36:25Z / 403 HTML; cf-mitigated: challenge | 1.115656 s | 5635 | a6240db191832a5bb53b2c0fc225ecd174fced6cb3098aa617ecc78f05742481 | no PDF bytes |
| Royal Society syndication https://royalsocietypublishing.org/rspa/article-pdf/277/1368/32/54549/rspa.1964.0004.pdf | 09:36:25Z / 403 HTML; cf-mitigated: challenge | 0.636290 s | 5725 | 1c55111e7cf724c59ffe2cf4ce90c9615a0235ddfc8685867eb5737da5c98da4 | no PDF bytes |
| JSTOR item https://www.jstor.org/stable/2414647, PDF https://www.jstor.org/stable/pdf/2414647.pdf?acceptTC=true | 09:36:26Z / 200 HTML; title Client Challenge | 0.894662 s | 3038 | 32ed63159c77e21ee19ca1b9aa3213ccf0218eb59539560b132a8e68ef0e18ea | challenge HTML, not PDF |

Royal Society responses contained server: cloudflare and cf-mitigated: challenge. JSTOR
contained content-type: text/html and no PDF magic header. No challenge-solving or
access-control bypass was attempted.

### Bishop--Hassan: author/institution routes

UCL Discovery is the relevant institutional repository for the UCL affiliation. Direct probes:

| path | Date / status / headers | elapsed | bytes | response-body SHA-256 | result |
|---|---|---:|---:|---:|---|
| UCL Discovery home https://discovery.ucl.ac.uk/ | 09:32:32Z / 403 HTML; Cloudflare challenge | 1.518395 s | 5495 | 83d54515f721628b24ce7a13b5785348a5972696c2641ab9346e605ba786c352 | inaccessible |
| UCL DOI search https://discovery.ucl.ac.uk/cgi/search/simple?q=10.1098%2Frspa.1964.0004 | 09:32:28Z / 403 HTML; Cloudflare challenge | 1.579578 s | 5728 | c8fed93393417cf77c8f98ab21404f52eaed54f3d77d6bcba3695bd90dda957b | no result page |
| UCL title search https://discovery.ucl.ac.uk/cgi/search/simple?q=The%20lift%20and%20drag%20forces%20on%20a%20circular%20cylinder | 09:32:29Z / 403 HTML; Cloudflare challenge | 0.568268 s | 5839 | 5870020e50a5782ad9422592060120e20405d2a419a7e3b41866a0fbaac808bb | no result page |
| UCL advanced search https://discovery.ucl.ac.uk/cgi/search/archive/advanced?title=The%20lift%20and%20drag%20forces%20on%20a%20circular%20cylinder | 09:32:29Z / 403 HTML; Cloudflare challenge | 0.717316 s | 5903 | 0ea376c948f66fbd59ab1598ba1ff22e0568b2bf2371e601b6d8706eab12933b | no result page |
| UCL OAI https://discovery.ucl.ac.uk/cgi/oai2?verb=ListRecords&metadataPrefix=oai_dc | 09:32:30Z / 403 HTML; Cloudflare challenge | 0.784020 s | 5730 | 549b70a1d8dd75d4fd56c629ebd3ef212b3484b68803597111bd250aa47e31ca | no OAI records |
| UCL Profiles DOI search https://profiles.ucl.ac.uk/search?query=10.1098%2Frspa.1964.0004 | 09:34:05Z / 200 HTML shell | 1.867498 s | 3327 | 326bb43cb363c72f150d3824ce6dcb5456603c9d0acd51539f0927e0bf5cfe87 | no author/PDF record in shell |

UCL Profiles is an official staff-profile service, not a repository. Its DOI/title requests
returned the same shell and no matching record; this is not proof that UCL has no internal
record.

### Tadrist et al.: publisher routes

| path | Date / status / headers | elapsed | bytes | response-body SHA-256 | result |
|---|---|---:|---:|---:|---|
| DOI resolver https://doi.org/10.1063/1.857804 | 09:29:38Z / 403 HTML; cf-mitigated: challenge | 1.984697 s | 5783 | 5dcae4f76cbf0d932597d517826850465df450a38d551a28e9156cd02102f3be | AIP redirect; no paper |
| AIP landing https://pubs.aip.org/pof/article/2/12/2176/401766/Experimental-investigation-of-fluctuating-forces | 09:36:15Z / 403 HTML; cf-mitigated: challenge | 0.612152 s | 5783 | cfb77863dbf7f1a2511cd5fae7763c4e1c81a36e90404f80060103d8f5c5440e | no landing page |
| AIP PDF https://pubs.aip.org/pof/article-pdf/2/12/2176/12512312/2176_1_online.pdf | 09:36:16Z / 403 HTML; cf-mitigated: challenge | 0.806710 s | 5687 | 63c76ab95409497ea5fe1cfd4177c1a400535b106ef89624f19022d334412cd3 | no PDF bytes |
| AIP PDF with download=true | 09:36:17Z / 403 HTML; cf-mitigated: challenge | 0.688189 s | 5744 | e2603226b75f6303e0de194bd6f18953f4cd5039f4b39341e4b6f59d4d2ca9c6 | no PDF bytes |
| Legacy Scitation https://aip.scitation.org/doi/pdf/10.1063/1.857804 | request batch / 200 HTML | 4.318306 s | 1434 | 7a0a0663d95e2f3a358e9a7076ca7472ea4a1c6acd26058a74de46d7bc4ccb20 | HTML shell, not PDF |
| Legacy Scitation https://aip.scitation.org/doi/epdf/10.1063/1.857804 | request batch / 200 HTML | 1.125789 s | 1426 | bd7fda6839106e7cd9e9c546e92c38b37362ec09305d580cff55207ade9df95c | HTML shell, not PDF |

AIP 403 bodies contained server: cloudflare and cf-mitigated: challenge. Scitation bodies
were ASCII HTML, not PDF bytes. No alternate host, token or access-control workaround was
attempted.

### Tadrist et al.: author/institution routes

HAL is the relevant national/institutional repository for the French affiliation:

| path | Date / status / headers | elapsed | bytes | body SHA-256 | result |
|---|---|---:|---:|---:|---|
| HAL DOI query https://api.archives-ouvertes.fr/search/?q=doiId_s:%2210.1063%2F1.857804%22&fl=*&wt=json | 09:30:39Z / 200 JSON | 1.301847 s | 116 | 316b596864c92e18c8f6688107ad49ecfc3fe632ebdc85efe82c65ed28ba8064 | numFound=0 |
| HAL exact-title query https://api.archives-ouvertes.fr/search/?q=title_t:%22Experimental%20investigation%20of%20fluctuating%20forces%22&fl=*&wt=json | 09:30:41Z / 200 JSON | 1.386444 s | 116 | 316b596864c92e18c8f6688107ad49ecfc3fe632ebdc85efe82c65ed28ba8064 | numFound=0 |
| HAL author query authFullName_t:Tadrist, rows=1000 | 09:31:03Z / 200 JSON | 2.847890 s | 51192 | c0769495bad0c48899ab7700b52d7a91bfc11e92cb4b608e606b8cee0301a7dc | 163 hits; no target DOI/title |
| AMU HAL API exact DOI | 09:35:55Z / 200 JSON | 1.896840 s | 116 | 316b596864c92e18c8f6688107ad49ecfc3fe632ebdc85efe82c65ed28ba8064 | numFound=0 |
| AMU HAL web https://amu.hal.science/search/index/?q=10.1063%2F1.857804 | request batch / 200 HTML bot challenge | 2.018926 s | 12550 | temporary body in log dir | no result |
| AMU institutional repository https://entrepot.recherche.univ-amu.fr/search/index/?q=10.1063%2F1.857804 | Date absent / HTTP 000 | 2.007983 s | 0 | not applicable | TLS SSL_ERROR_SYSCALL; no inference |
| IUSTI current site https://iusti.cnrs.fr/ | Date absent / HTTP 000 | 30.001650 s | 0 | not applicable | timeout; no inference |
| IUSTI historical URL http://iusti.univ-provence.fr/document.php?pagendx=11201&project=iusti | Date absent / 502 | 2.060553 s | 0 | e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855 | gateway failure; no document |

The exact HAL zero-result records are stronger than the challenged HTML page, but do not rule
out an unindexed private author copy. The broad author query did not license any other paper's
file for this audit.

## Reproducible commands and logs

All response bodies, headers, curl logs and stderr are retained outside the repository in:

~~~text
/tmp/tadrist-bishop-author-manuscript-audit-20260822/
~~~

Representative commands:

~~~sh
probe_dir=/tmp/tadrist-bishop-author-manuscript-audit-20260822
mkdir -p "$probe_dir"

curl -L --max-time 30 -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D "$probe_dir/bishop-crossref.headers" \
  -o "$probe_dir/bishop-crossref.body" \
  'https://api.crossref.org/works/10.1098/rspa.1964.0004'

curl -L --max-time 30 -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D "$probe_dir/tadrist-crossref.headers" \
  -o "$probe_dir/tadrist-crossref.body" \
  'https://api.crossref.org/works/10.1063/1.857804'

curl -L --max-time 30 -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D "$probe_dir/bishop-royal-pdf.headers" \
  -o "$probe_dir/bishop-royal-pdf.body" \
  'https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.1964.0004'

curl -L --max-time 30 -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D "$probe_dir/tadrist-aip-pdf.headers" \
  -o "$probe_dir/tadrist-aip-pdf.body" \
  'https://pubs.aip.org/pof/article-pdf/2/12/2176/12512312/2176_1_online.pdf'

curl -L --max-time 30 -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D "$probe_dir/ucl-oai.headers" \
  -o "$probe_dir/ucl-oai.body" \
  'https://discovery.ucl.ac.uk/cgi/oai2?verb=ListRecords&metadataPrefix=oai_dc'

curl -L --max-time 30 -A 'Mozilla/5.0 (research retrieval; contact: open-source)' \
  -D "$probe_dir/hal-doi.headers" \
  -o "$probe_dir/hal-doi.body" \
  'https://api.archives-ouvertes.fr/search/?q=doiId_s:%2210.1063%2F1.857804%22&fl=*&wt=json'

sha256sum "$probe_dir"/*.body
file "$probe_dir"/{bishop-royal-pdf,tadrist-aip-pdf,bishop-jstor-pdf}.body
~~~

完整 per-request log 是临时目录中的 paired *.headers、*.curl-log 和 *.stderr 文件。
本报告中的 body hashes 是 challenge、metadata 或 error response 的 hashes，均不是
source-paper artifact hash。

## Artifact, diff and validation status

No full bytes were obtained, so this audit intentionally created no directory below
/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/. Consequently:

* no PDF was copied, rendered, digitized or assigned a candidate hash;
* pdfinfo and qpdf --check are not applicable and were not claimed as passed;
* no force value, span, normalization, blockage, boundary or uncertainty fact was extracted;
* no existing repository file was edited.

Validation after writing this report:

~~~sh
git diff --check -- docs/research/2026-08-22-v04-tadrist-bishop-author-manuscript-audit.md
git status --short --untracked-files=all
git diff --stat -- docs/research/2026-08-22-v04-tadrist-bishop-author-manuscript-audit.md
~~~

Expected diff scope is one new file at the requested path. The negative result is an access
finding, not a statement that a licensed library or private author copy cannot exist.

