# Bishop--Hassan (1964) archival-primary follow-up

检索日期：2026-08-22（UTC）。本报告是 `RESEARCH_ONLY` 的 archival-primary
检索记录，只处理 R. E. D. Bishop and A. Y. Hassan, *The lift and drag forces on a
circular cylinder in a flowing fluid*, *Proceedings of the Royal Society of London
A* **277**(1368), 32--50 (1964), DOI
[`10.1098/rspa.1964.0004`](https://doi.org/10.1098/rspa.1964.0004)。不读取
HUNDUN long-stat 输出，不修改 references、receipt、ledger、solver、COAST 或
Stage 5，不做数字化、容差、等价性或 `ACCEPT/REJECT` 决策。

## 结论

这次检索发现了一个此前未记录的原刊档案入口：JSTOR 的原始期刊条目
[`2414647`](https://www.jstor.org/stable/2414647)，题名与目标论文一致。它证明
JSTOR 保存了该原刊条目的数字档案记录，但本检索环境对原刊 PDF/全文入口返回
JSTOR client challenge；得到的是 3038-byte HTML，不是 PDF。没有绕过该挑战、登录
或付费墙，也没有把 challenge 页面当作 primary artifact。

Google Books 的合订本目录也给出一个明确的 1964 年卷记录：
[`l4kkAAAAMAAJ`](https://books.google.com/books/feeds/volumes/l4kkAAAAMAAJ)。该
记录标为 646 pages、Harrison and Son，但同时明确
`openAccess=disabled`、`viewability=view_no_pages`、`not_embeddable`；它是公开
目录元数据，不是可读扫描。按目标题名、作者和 277(1368) 查询时，没有发现另一份
公开可读原刊或作者/机构稿件。

因此本轮没有取得合法可读的完整 primary scan/manuscript，没有创建
`/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/` 下的
Bishop--Hassan artifact directory，也没有任何可封印的 PDF hash、`pdfinfo`、页级
locators 或数值。`Cd_mean`、`Cl_rms`、归一化、active span、端部/支撑、边界、
blockage/correction 和 uncertainty 仍不可从该文 primary pages 核验；该文不能
因本次检索而成为数值 authority。

## 新发现的原始目录与受限入口

### JSTOR 原刊记录

[`https://www.jstor.org/stable/2414647`](https://www.jstor.org/stable/2414647) 的
页面标题为 **The Lift and Drag Forces on a Circular Cylinder in a Flowing Fluid**。
这是目标论文的原始期刊档案稳定标识（JSTOR item `2414647`），不是后续论文或二手
表格。JSTOR 的 XML 入口
[`https://www.jstor.org/doi/xml/10.2307/2414647`](https://www.jstor.org/doi/xml/10.2307/2414647)
重定向回同一稳定条目。

原始 PDF 入口按 JSTOR 页面给出的稳定 ID 试探：

```text
https://www.jstor.org/stable/pdf/2414647.pdf?acceptTC=true
```

响应为 `HTTP/2 200`、`Content-Type: text/html; charset=utf-8`，页面标题为
`Client Challenge`，并要求 JavaScript；不是 PDF。没有继续尝试破解 challenge、
模拟浏览器验证或使用绕过服务。

### Google Books 合订本目录

Google Books legacy feed 查询
[`Proceedings Royal Society London 277 no 1368`](https://books.google.com/books/feeds/volumes?q=Proceedings%20Royal%20Society%20London%20277%20no%201368)
返回了 volume ID `l4kkAAAAMAAJ`。其公开 machine record
[`l4kkAAAAMAAJ`](https://books.google.com/books/feeds/volumes/l4kkAAAAMAAJ) 的
关键字段为：

```text
title: Proceedings of the Royal Society of London
subtitle: Containing papers of a mathematical and physical character. Series A
creator: Royal Society (Great Britain)
date: 1964
publisher: Harrison and Son
format: 646 pages; book
embeddability: not_embeddable
openAccess: disabled
viewability: view_no_pages
```

按题名/作者的 legacy feed
[`Bishop Hassan circular cylinder`](https://books.google.com/books/feeds/volumes?q=Bishop%20Hassan%20circular%20cylinder)
没有返回目标论文或可读原刊卷；结果是其他书籍/后续资料，未被本报告用于任何数值
判断。Google Books API 另一次查询只返回 quota `429`，所以 legacy feed 的公开
记录是本项的可复现证据。

## 访问证据

以下 body 是临时检索日志，保留在 `/tmp/bishop-hassan-archival-20260822/`，不属于
候选 artifact，也没有复制到仓库。时间取响应 `Date` header；`sha256` 是返回 body
的 hash。对于 challenge/错误 body，hash 只用于证明访问结果，绝不作为论文内容
hash。

| 入口 | UTC 时间 | 响应/类型 | bytes | 返回 body SHA-256 | 判定 |
|---|---|---|---:|---|---|
| JSTOR PDF `stable/pdf/2414647.pdf?acceptTC=true` | 2026-08-22 08:37:12 | 200, `text/html`, Client Challenge | 3038 | `32ed63159c77e21ee19ca1b9aa3213ccf0218eb59539560b132a8e68ef0e18ea` | 无 PDF |
| JSTOR XML `doi/xml/10.2307/2414647` | 2026-08-22 08:38:57 | 302 到 stable 条目，再 200 challenge HTML | 3038 | `32ed63159c77e21ee19ca1b9aa3213ccf0218eb59539560b132a8e68ef0e18ea` | 无 XML/全文 |
| Royal Society landing | 2026-08-22 08:42:03 | 403, `cf-mitigated: challenge` | 5831 | `c3522a74ae536410374777806a0327b98dc96e4494953cc8af78b4e30d1a00a1` | 无页面 |
| Royal Society DOI PDF | 2026-08-22 08:42:04 | 403, `cf-mitigated: challenge` | 5614 | `8e8e8e69a613962624a7d63631a53031c93650d2fddc1140d2d65645fa51dd02` | 无 PDF |
| Royal Society syndication PDF (`54549`) | 2026-08-22 08:42:05 | 403, `cf-mitigated: challenge` | 5704 | `f48203afdbf9d0d855377ebfdc15f09eaaed25ef9cfce9980d6d249361ac7758` | 无 PDF |
| Royal Society `doi/epdf` | 2026-08-22 08:42:07 | 403, `cf-mitigated: challenge` | 5617 | `e773ae6c77257e49a911f873ceab3bfe88c8274caa2bdf590c7d215399216ba6` | 无 PDF |
| Royal Society old HTTP landing | 2026-08-22 08:42:06 | 502, empty body | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | 无页面 |
| Google Books volume `l4kkAAAAMAAJ` | 2026-08-22 08:41:04 | 200, `application/atom+xml` | 2102 | `17c886575fb34a99c58edb430abdf8870e0bbff458004838369c02e3c9cd5dfb` | 目录；`view_no_pages` |
| Google Books 277(1368) query | 2026-08-22 08:44:17 | 200, `application/atom+xml` | 24824 | `57947699d15bb507e3cd6fb0c470d7c815e21537e03877ce9440edd0314e880c` | 只有合订本目录，无扫描 |
| HathiTrust title query | 2026-08-22 08:39:19 | 403, `cf-mitigated: challenge` | 5827 | `36f335506891e4ed859036ffd581f279165c35767d7d491860d75e9923af8910` | 无可读记录 |
| Internet Archive DOI advanced search | 2026-08-22 08:45（约） | `curl` exit 28, HTTP 000, 8.001 s；仅收到 proxy connection line | 0 | 未生成 body | 网络超时；不推断 archive 中不存在 |
| WorldCat exact-title query | 2026-08-22 08:39:28 | 200, `text/html` | 262745 | `f16a1969567da4f4280767e3c80057da40492183ca0d1295c94db43d617cfefd` | 页面/嵌入数据未给出目标数字 artifact |
| UCL Discovery author query | 2026-08-22 08:42:51 | 403, `cf-mitigated: challenge` | 6056 | `c0e04d8cd2444d3b811a690551a000a2e90e79eb040b9491bbe48740fbde7f04` | 无 manuscript |

Royal Society 的四个 403 body 均含 `cf-mitigated: challenge`；HathiTrust 和 UCL
也同样是服务端 challenge。Internet Archive 的超时是检索环境的网络失败，不是
“无结果”的证据。WorldCat 的页面查询只用于馆藏/目录线索；它没有提供可下载的
原刊页面，不能充当 primary artifact。

## 检索命令

以下命令可在不覆盖候选/receipt 的临时目录中复现核心检查。它们只保存响应，不
绕过登录、付费墙、DRM 或 JavaScript challenge。

```sh
tmp_dir=/tmp/bishop-hassan-archival-20260822
mkdir -p "$tmp_dir"

curl -L --max-time 20 \
  -A 'Mozilla/5.0 (archival research; contact: open-source)' \
  -D "$tmp_dir/jstor.headers" \
  -o "$tmp_dir/jstor.body" \
  'https://www.jstor.org/stable/pdf/2414647.pdf?acceptTC=true'

curl -L --max-time 20 \
  -A 'Mozilla/5.0 (archival research; contact: open-source)' \
  -D "$tmp_dir/royal-pdf.headers" \
  -o "$tmp_dir/royal-pdf.body" \
  'https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.1964.0004'

curl -L --max-time 15 \
  -A 'Mozilla/5.0 (archival research; contact: open-source)' \
  -D "$tmp_dir/google-volume.headers" \
  -o "$tmp_dir/google-volume.xml" \
  'https://books.google.com/books/feeds/volumes/l4kkAAAAMAAJ'

curl -L --max-time 8 \
  -A 'Mozilla/5.0 (archival research; contact: open-source)' \
  -D "$tmp_dir/ia.headers" \
  -o "$tmp_dir/ia.body" \
  -w 'curl_exit=%{http_code} time_total=%{time_total}\n' \
  'https://archive.org/advancedsearch.php?q=%2210.1098%2Frspa.1964.0004%22&fl%5B%5D=identifier&fl%5B%5D=title&rows=20&output=json'

sha256sum "$tmp_dir"/{jstor.body,royal-pdf.body,google-volume.xml}
```

## Page-level audit status

Because no complete paper artifact was lawfully obtained, the following rows remain
unresolved and intentionally contain no inferred values:

| Required primary fact | Status |
|---|---|
| Re≈3900 direct total mean drag value | 未核验 |
| Re≈3900 fluctuating total lift RMS value | 未核验 |
| force-balance observable and inertia/support separation | 未核验 |
| coefficient normalization, reference velocity, area and span | 未核验 |
| active measured span, dummy ends, end plates and support arrangement | 未核验 |
| channel/free-surface/boundary condition | 未核验 |
| blockage ratio and correction | 未核验 |
| calibration, sampling/statistical window and uncertainty | 未核验 |
| mapping to HUNDUN periodic span `πD` | 未评定 |

The previous primary-source search already established only the narrow Crossref/publisher
abstract facts (directly measured fluctuating lift/drag and steady mean drag; stationary
range 3600--11000). The archival follow-up adds the JSTOR stable record and the
Google Books volume record, but adds no page-level force evidence. A future lawful
continuation needs an institutionally licensed scan, an openly deposited author
manuscript, or a publisher session that makes the original pages readable; then freeze
the source bytes, MIME type, page count, `pdfinfo`, SHA-256 and page locators before any
scientific audit.

## Validation

* Only this new Markdown report was added; no existing file was edited.
* No Bishop--Hassan PDF/scan/manuscript was copied into the benchmark directory.
* No numeric value was digitized or accepted.
* `git diff --check` passes for this report.
