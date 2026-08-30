# Parnaudeau et al. (2008) Fig. 11--15：一手 authority 检索与冻结结论

检索日期：2026-08-21（UTC）。本笔记只处理 `10.1063/1.2957018` 的五类近尾迹
profile、可用于受控数字化的原始图件，以及 `Re≈3900` 的直接 total force authority。
没有读取 HUNDUN 数值或性能输出，没有使用后续 CFD 的汇总数值，也没有联系作者。ResearchGate
页面只用于确认作者上传 provenance 和可发现路径；作者上传全文原则上可以成为受控数字化源，
但本轮没有取得其 bytes/hash，因此不把“页面可发现”提升为已取得的 source artifact。

## 冻结判定

| 对象 | 判定 | 判定理由 |
|---|---|---|
| Fig. 11--15 的 author arrays（`x/D=1.06, 1.54, 2.02`） | `INSUFFICIENT` | 没有找到由作者、机构或出版社托管且同时给出下载 URL、版本/bytes/hash、字段、归一化和 station mapping 的 arrays。论文脚注只说 statistics 可通过联系作者获得。 |
| 可供受控数字化的 primary PDF/figure source | `INSUFFICIENT` | Crossref 给出 AIP VOR PDF 的稳定 URL，但本次请求得到 403 HTML 而非 PDF；两个 HAL 记录明确无文件；ResearchGate 明确标注作者上传的全文原则上可作 primary digitization source，但本次未取得 PDF bytes/hash，故当前仍不能晋级。 |
| `Re≈3900` 直接测量的 total `Cd` | `INSUFFICIENT` | Parnaudeau DOI/HAL metadata 只描述 PIV/HWA/LES 的 wake statistics；Norberg 的可核验机构记录明确是 pressure coefficients/pressure forces，未证明 `pressure+viscous` 的 total force、精确 Re=3900 和对应原始页/文件。 |
| `Re≈3900` 有限跨距 total `Cl_rms` | `INSUFFICIENT` | Lund 的 Norberg 2003 一手机构记录明确为 nominal infinitely-long/nonconfined cylinder 的 **sectional** r.m.s. lift；这不能替代 finite-span total lift。没有找到带跨距、端部条件、直接力测量和精确 Re 的一手记录。 |

Norberg 1987 报告、Norberg/Sundén 1987 论文、Norberg 1993 pressure-force chapter 和 Norberg
2003 lift review 都可作背景或后续人工核验的 `SUPPLEMENT_ONLY` 候选，不能解除以上三项
`INSUFFICIENT`。

## 判据和术语

本轮沿用如下冻结边界：

* `ARRAY_ADMISSIBLE` 需要托管主体是作者/机构/出版社，且可复核 URL、version/date、bytes、
  SHA-256、字段定义、归一化和三个 station 的对应关系。论文里的曲线、截图、OCR 或后续
  论文复绘都不满足这个条件。
* `DIGITIZATION_SOURCE_ADMISSIBLE` 需要实际取得稳定、provenance 可核验的 publisher VOR 或
  author-uploaded primary PDF/figure source，并确认 exact bytes/hash、printed page、
  figure/legend、坐标轴、曲线来源和分辨率；不能仅凭 DOI metadata、下载按钮或 preview 晋级。
* 本笔记中的 `direct total Cd` 指同一有限圆柱参考面积上的 pressure **加** viscous/shear
  contribution 的总阻力（并明确 Re、跨距/端板和测量方式）。pressure drag、base-pressure
  coefficient、sectional/per-unit-span coefficient、公式派生值和后续 CFD 表格均不等价。
* `finite-span total Cl_rms` 指整个暴露有限跨距上的时间变化升力积分；sectional lift、
  vanishing-segment lift 或从局部压力/速度公式推导的值不等价。

## 2008 论文的 first-party/author/institution/publisher 路径

### DOI/Crossref（只作 DOI 与 publisher 文件指针）

* DOI landing：[`https://doi.org/10.1063/1.2957018`](https://doi.org/10.1063/1.2957018)。
* Crossref REST record：[`https://api.crossref.org/works/10.1063%2F1.2957018`](https://api.crossref.org/works/10.1063%2F1.2957018)。
  2026-08-21 UTC 返回 HTTP 200、`application/json`，12,312 bytes，SHA-256
  `8b62696e861c99825b94d28d95aecdb80261a62900ec75b8a1682292ea79a634`。record 给出：AIP
  Publishing、Physics of Fluids 20(8)、article number `085101`、published 2008-08-01、
  DOI，以及 `relation: {}`。
* Crossref 的 `resource.primary.URL` 为
  [`https://pubs.aip.org/pof/article/20/8/085101/256405/Experimental-and-numerical-studies-of-the-flow`](https://pubs.aip.org/pof/article/20/8/085101/256405/Experimental-and-numerical-studies-of-the-flow)。
  同一 record 的 `link` 明确标为 `content-type=application/pdf`、`content-version=vor`：
  [`https://pubs.aip.org/aip/pof/article-pdf/doi/10.1063/1.2957018/15878439/085101_1_online.pdf`](https://pubs.aip.org/aip/pof/article-pdf/doi/10.1063/1.2957018/15878439/085101_1_online.pdf)。
  这是稳定的 publisher file pointer，不是本次已取得的 PDF；Crossref 没有 supplementary
  relation 或可下载数组字段。

### AIP publisher retrieval

在本机以 `Mozilla/5.0` 和 `Accept: application/pdf` 访问 DOI redirect、AIP landing 和上述
VOR URL：

| URL | 观察结果 |
|---|---|
| publisher landing | HTTP 403，`text/html; charset=UTF-8`，5,698 bytes；`file` 判定 HTML |
| VOR PDF URL | HTTP 403，`text/html; charset=UTF-8`，5,655 bytes；`file` 判定 HTML，不是 PDF |
| 旧 `https://aip.scitation.org/doi/pdf/10.1063/1.2957018` | HTTP 200 但返回 JavaScript/challenge HTML（约 1.4--1.5 kB），不是 PDF |

本轮不会把 challenge HTML 的 hash 当成论文 hash；没有取得合法的 PDF bytes，因此不能确认
Fig. 11--15 的原始 page raster/vector resolution，也不能建立 `DIGITIZATION_SOURCE_ADMISSIBLE`。

### HAL/INRAE institution records

两个机构记录指向同一个 DOI，但均是 metadata-only：

1. 原始 HAL：[`https://hal.science/hal-00383669v1`](https://hal.science/hal-00383669v1)，
   API：[`https://api.archives-ouvertes.fr/search/?q=halId_s:hal-00383669&wt=json&fl=*`](https://api.archives-ouvertes.fr/search/?q=halId_s:hal-00383669&wt=json&fl=*)。
   record 为 `hal-00383669`, version 1, DOI `10.1063/1.2957018`, `openAccess_bool=false`。
2. INRAE 镜像：[`https://hal.inrae.fr/hal-02590739v1`](https://hal.inrae.fr/hal-02590739v1)，
   API：[`https://api.archives-ouvertes.fr/search/?q=halId_s:hal-02590739&wt=json&fl=*`](https://api.archives-ouvertes.fr/search/?q=halId_s:hal-02590739&wt=json&fl=*)。
   record 为 `hal-02590739`, version 1, DOI 同上，`openAccess_bool=false`。
3. 用只返回文件字段的 API 查询核对：
   `fl=uri_s,openAccess_bool,fileMain_s,files_s` 对两个 record 都返回一个 document、
   `openAccess_bool:false`，只有 `uri_s`，没有 `fileMain_s` 或 `files_s`。原始 JSON 各为
   210 bytes；`hal-00383669` SHA-256
   `1618e07fae947d9a6748dd8376ac825e8befa50d3031a6e1efbef77faa672caf`，
   `hal-02590739` SHA-256
   `d24889456b95caeeba407416ce3f8547ac9dcd6b1dffb531195f5d5745e538e6`。
4. Data HAL record：[`https://data.hal.science/document/hal-02590739v1`](https://data.hal.science/document/hal-02590739v1)，
   以及其官方 RDF/JSON exports：
   [`https://data.archives-ouvertes.fr/document/hal-02590739v1.json`](https://data.archives-ouvertes.fr/document/hal-02590739v1.json)、
   [`https://data.archives-ouvertes.fr/document/hal-02590739v1.rdf`](https://data.archives-ouvertes.fr/document/hal-02590739v1.rdf)。
   metadata 记录 DOI、`hal-02590739v1`、14 pages、PIV/HWA/LES subjects 和 abstract，但没有
   file relation；JSON 15,639 bytes，SHA-256
   `d4be5d9887a9d0b2e6becf2d201cf8751c214ab1e3064e6a02f261ad5ed0e145`；RDF 6,613 bytes，
   SHA-256 `dde32edc50f37df0822fb5a8c91ce903315f3832e43dd71118f7ead54d5e2e33`。
5. HAL record 的 external `sameAs` 还指向 [`https://irsteadoc.irstea.fr/cemoa/PUB00024478`](https://irsteadoc.irstea.fr/cemoa/PUB00024478)。
   该旧 IRSTEA endpoint 的 HTTP 路径现在 302 到 HAL，直接 HTTPS 路径证书失效；其他 web
   fetch 还出现 502。三条路径都没有返回 attachment bytes，不能将它作为另一份可下载
   primary PDF；HAL `sameAs` 关系本身也不提供 arrays。

HAL 页面在可视化页面中显示 `Fichier non déposé`，与 API 没有 file fields 一致。用 `/document`
路径的自动请求遇到 HAL anti-bot HTML；该 HTML 不是论文文件，也没有被用作证据。

## ResearchGate author upload：原则上可用，当前未取得

发现页：[`https://www.researchgate.net/publication/234130406_Experimental_and_numerical_studies_of_the_flow_over_a_circular_cylinder_at_Reynolds_number_3900`](https://www.researchgate.net/publication/234130406_Experimental_and_numerical_studies_of_the_flow_over_a_circular_cylinder_at_Reynolds_number_3900)。
页面明确列出 `Public Full-texts 2`，并标注两个文件由 Dominique Heitz 上传为 `Author content`：

* 显示名 `2008_parnaudeau_etal_PoF.pdf`，attachment link id/path
  `0912f50f729d29a832000000`；
* 显示名 `0912f50a3828de19d1000000.pdf`，attachment link id/path
  `0912f50a3828de19d1000000`。

页面解析出的两个 download href（平台可能动态重写）分别是：

* [`.../links/0912f50f729d29a832000000/Experimental-and-numerical-studies-of-the-flow-over-a-circular-cylinder-at-Reynolds-number-3900.pdf`](https://www.researchgate.net/publication/profile/Dominique-Heitz/publication/234130406_Experimental_and_numerical_studies_of_the_flow_over_a_circular_cylinder_at_Reynolds_number_3900/links/0912f50f729d29a832000000/Experimental-and-numerical-studies-of-the-flow-over-a-circular-cylinder-at-Reynolds-number-3900.pdf)
* [`.../links/0912f50a3828de19d1000000/Experimental-and-numerical-studies-of-the-flow-over-a-circular-cylinder-at-Reynolds-number-3900.pdf`](https://www.researchgate.net/publication/profile/Dominique-Heitz/publication/234130406_Experimental_and_numerical_studies_of_the_flow_over_a_circular_cylinder_at_Reynolds_number_3900/links/0912f50a3828de19d1000000/Experimental-and-numerical-studies-of-the-flow-over-a-circular-cylinder-at-Reynolds-number-3900.pdf)

本机对这两个 href 得到 Cloudflare/anti-bot HTTP 403 HTML（约 20,178 bytes；web fetch 也出现
404），没有 PDF magic bytes、文件长度或 SHA-256。现有 authority-gap policy 允许经作者 provenance
核验的 author-uploaded full text 作为 controlled digitization source；所以这里的否决对象是“本轮
尚未取得 artifact”，不是 ResearchGate 这个托管域名本身。以后若取得文件，仍须由主 agent 保存
exact bytes/hash、核验作者上传 provenance、页码/legend/分辨率和可复核下载记录，才能重新评估；
本轮不把它列为 `DIGITIZATION_SOURCE_ADMISSIBLE`，也不会把受版权保护的全文提交进仓库。

## Fig. 11--15 的可核验 locator（只来自 discovery copy）

ResearchGate 页面可显示 article text，但按上节规则不把它作为冻结 primary PDF。为了让后续主
agent 在取得 AIP/正式 author PDF 后能做逐项复核，记录其 locator，不把曲线数值抄入任何 receipt：

| figure | quantity（discovery text） | stations/axes | printed-page locator in that copy |
|---|---|---|---|
| Fig. 11 | mean streamwise velocity `⟨u⟩/Uc` | `x/D=1.06, 1.54, 2.02`; vertical `y/D` | `085101-11` |
| Fig. 12 | mean normal velocity `⟨v⟩/Uc` | 同上；vertical `y/D` | `085101-11` |
| Fig. 13 | streamwise fluctuation variance `⟨u′u′⟩/Uc²` | 同上；vertical `y/D` | `085101-11` |
| Fig. 14 | normal fluctuation variance `⟨v′v′⟩/Uc²` | 同上；vertical `y/D` | `085101-11` |
| Fig. 15 | covariance `⟨u′v′⟩/Uc²` | 同上；vertical `y/D` | `085101-11` |

同一 discovery text 的 Fig. 9 caption 将曲线来源区分为 present LES、present PIV cases 1/2、
present HWA 和若干既有实验/模拟；因此未来数字化时不能把 panel 中所有 marker 当成 PIV。
当前不冻结“哪个 PIV case 属于哪一条曲线”，必须从取得的 primary PDF 的完整 legend 和
caption 逐字复核。Data HAL metadata 的“14 pages”只确认书目信息，不确认这些 figure pages。

脚注 31 的 discovery locator 是 printed page `085101-14`，短释义为“本文的 experimental and
numerical statistics 可通过联系作者获得”。这说明论文没有在该页面随附公开 arrays；也不构成
本轮联系作者的授权。

### 若主 agent 后续取得合规 primary PDF：受控数字化门槛

以下只是待主 agent 决定的执行规则，不是本轮通过证据：

1. 保存原始下载 URL、HTTP headers、UTC 时间、文件 bytes、SHA-256，并用 `file`/PDF parser
   确认真实 PDF；记录 printed page、figure number、panel layout、全量 legend/caption 和
   `pdfimages -list` 或等价的 raster/vector 分辨率。
2. 每个 panel 独立用至少两个 x/y tick 标定 pixel→data 仿射变换，记录轴刻度与拟合残差；
   对垂直平移的三个 station 逐 panel 解偏移，不能借用另一 panel 的变换。仅采样明确标为
   target PIV/experimental 的曲线，保留原始 pixel 坐标、曲线厚度和遮挡记录。
3. station mapping 只接受图中明确的 `1.06/1.54/2.02`。若未来文件只有别的 station，是否
   插值、允许的 `Δ(x/D)` 和插值误差必须由主 agent 另行书面决定；本笔记不授权最近邻替代。
4. 误差至少拆为轴标定/刻度读数、pixel quantization、线宽/抗锯齿、扫描几何畸变（若有）和
   重复数字化差异，按主 agent 选定的 RSS 或保守上界合成；不能把论文 abstract 中约 10%
   的统计不确定度直接冒充图像数字化误差。若曲线来源/legend 或坐标原点方向无法确认，
   结果必须退回 `INSUFFICIENT`。

## Re≈3900 force authority 分层

### 一手/机构记录

| source | official locator and primary fact | category boundary | result |
|---|---|---|---|
| Parnaudeau et al. 2008 | [AIP DOI/Crossref record](https://doi.org/10.1063/1.2957018)；Crossref abstract 说研究 Re=3900 的 near-wake turbulence statistics/power spectra，实验为 HWA/PIV、数值为 LES；HAL/Data HAL 记录同一 DOI、PIV/HWA/LES、14 pages，但无 force file/field | metadata 没有 direct total force；不能从 wake profile 公式派生 total Cd 或 finite-span Cl 作为“直接测量” | `INSUFFICIENT` |
| Norberg 1987 report | [Lund/Chalmers institutional record](https://portal.research.lu.se/en/publications/effect-of-reynolds-number-and-a-low-intensity-freestream-turbulen/)；report `Publikation 87/2`, 54 pages, Re about 50--2×10⁵, Tu=1.4%, abstract 明确 mean/RMS **pressure coefficients**；official PDF pointer [Norberg_Publikation_87_2.pdf](https://portal.research.lu.se/files/137949583/Norberg_Publikation_87_2.pdf) 当前 403 | pressure coefficient/pressure-force category；公开 metadata 没有 pressure+viscous total decomposition、exact Re=3900 row 或 finite-span integration；PDF 未取得 | `SUPPLEMENT_ONLY` / target `INSUFFICIENT` |
| Norberg & Sundén 1987 | [Lund institutional record](https://portal.research.lu.se/sv/publications/turbulence-and-reynolds-number-e%EF%AC%80ects-on-the-%EF%AC%82ow-and-fluid-forces-/)，实验 Re `2×10⁴--3×10⁵`，端板圆柱，pressure microphones/hot film，描述的是 pressure forces | Re 不在 3900 附近；且 metadata 未给 total finite-span coefficient | `SUPPLEMENT_ONLY` |
| Norberg 1993 | [Lund institutional record](https://portal.research.lu.se/en/publications/pressure-forces-on-a-circular-cylinder-in-cross-%EF%AC%82ow) 与 [Springer chapter DOI](https://doi.org/10.1007/978-3-662-00414-2_60)，标题和记录明确是 **Pressure forces**，pp. 275--278；[Springer PDF endpoint](https://link.springer.com/content/pdf/10.1007/978-3-662-00414-2_60.pdf) 返回 challenge HTML（3,038 bytes） | pressure-only chapter；没有从机构/出版社可取得的 total pressure+viscous finite-span table | `SUPPLEMENT_ONLY` |
| Norberg 2003 | [Lund institutional record](https://portal.research.lu.se/en/publications/fluctuating-lift-on-a-circular-cylinder-review-and-new-measuremen/) 与 [Elsevier DOI](https://doi.org/10.1016/S0889-9746(02)00099-3)，官方摘要/record 明确 nominal infinitely-long/nonconfined cylinder、sectional r.m.s. lift，Re≈47--2×10⁵ | sectional/per-unit-span lift 不是 finite-span total `Cl_rms`; 不能替代目标 authority | `SUPPLEMENT_ONLY` |

### 三种 force 结果

* **Direct total `Cd`：`INSUFFICIENT`。** 1987 report 是最接近的 primary lead，且其 Re range
  包含 3900，但当前可验证的 institution metadata 只保证 pressure-coefficient material，
  不能证明目标是包含 viscous/shear 的 total drag，也没有可锁定的 exact-Re primary table
  bytes/page。没有采用后续 CFD 文献中常见的 `≈0.98` 等汇总值。
* **Finite-span total `Cl_rms`：`INSUFFICIENT`。** Norberg 2003 的一手机构记录明确把结果称为
  sectional r.m.s. lift；Parnaudeau 2008 记录是 wake PIV/HWA/LES statistics，不是已公开的
  finite-span force record。没有采用任何后续 CFD comparison table 或从 sectional coefficient
  乘跨距的公式猜值。
* **Pressure/sectional/formula-derived 值：`SUPPLEMENT_ONLY`。** 它们可帮助主 agent 设计
  未来 force-equivalence 查询，但不满足本任务的 direct total/finite-span gate。

这项 force 缺口不自行改变 literature receipt 的 completion schema：当前 machine receipt 只把
Parnaudeau profiles 列为 `incomplete_references`。若最终 physical-accuracy policy 明确只以直接
Strouhal、recirculation 和已注册 velocity/stress profiles 判定，direct total-force oracle 可以保持
未采用；但任何最终文件都必须明确写成“未做 direct-force comparison”，不能暗示已完成该门。

## 查询范围、命令和日志摘要

查询覆盖：

* DOI、文章题名、`10.1063/1.2957018` + `supplementary/data/array/profile`；AIP publisher
  landing/PDF 与 Crossref `link/relation`；
* `hal-00383669v1`、`hal-02590739v1`、Data HAL JSON/RDF/Turtle/N3、作者 HAL CV pages；
* Lund/Chalmers institution records 和 DOI：Norberg 1987 report、Norberg/Sundén 1987、
  Norberg 1993 pressure-force chapter、Norberg 2003 lift review；
* author/institution/PIV profile terms、exact station strings `1.06`, `1.54`, `2.02`；
* ResearchGate author-upload page 及两个 attachment link ids（仅记录，不晋级）。

关键只读命令（所有缓存均在 `/tmp`，没有把大文件写入仓库）：

```sh
curl -L -A 'Mozilla/5.0' -o /tmp/audit-crossref.json \
  'https://api.crossref.org/works/10.1063%2F1.2957018'
sha256sum /tmp/audit-crossref.json

curl -L -o /tmp/audit-hal1-files.json \
  'https://api.archives-ouvertes.fr/search/?q=halId_s:hal-00383669&wt=json&fl=uri_s,openAccess_bool,fileMain_s,files_s'
curl -L -o /tmp/audit-hal2-files.json \
  'https://api.archives-ouvertes.fr/search/?q=halId_s:hal-02590739&wt=json&fl=uri_s,openAccess_bool,fileMain_s,files_s'

curl -L -A 'Mozilla/5.0' -o /tmp/audit-aip.pdf \
  'https://pubs.aip.org/aip/pof/article-pdf/doi/10.1063/1.2957018/15878439/085101_1_online.pdf'
file /tmp/audit-aip.pdf
sha256sum /tmp/audit-aip.pdf

curl -L -A 'Mozilla/5.0' -o /tmp/audit-norberg-report.pdf \
  'https://portal.research.lu.se/files/137949583/Norberg_Publikation_87_2.pdf'
file /tmp/audit-norberg-report.pdf
```

本轮代表性输出：AIP VOR `HTTP 403`, `text/html`, 5,655 bytes；Norberg report pointer
`HTTP 403`, `text/html`, 5,608 bytes；两者均经 `file` 判定为 HTML。挑战页 hash 仅用于
审计日志，不是可复用的 scientific artifact。所有 primary metadata URL、HTTP/attachment
事实、bytes/hash 均已在上文列出。

## 未解决点与下一步

1. 需要从 AIP、作者/机构静态路径，或 provenance 可核验的 ResearchGate 作者上传取得真实
   VOR/author PDF，并在主 agent 规则下封存 bytes/hash、页码、legend、坐标轴和分辨率；
   ResearchGate 的页面记录和受阻下载按钮本身不能代替该 artifact。
2. 若要把 Fig. 11--15 补成 arrays，需要原作者/机构合法发布的 profile file 或主 agent
   认可的 primary-PDF digitization；当前脚注“contacting the authors”不是公开 receipt。
3. 需要一手的 exact-Re≈3900、`pressure+viscous` total Cd 和有限跨距 total `Cl_rms` 记录，
   明确参考面积、端板/跨距、测量 vs 公式派生；当前三类 force authority 均未满足。
4. 主 agent 仍需决定 station interpolation tolerance、digitization coordinate/error budget
   组合规则，以及是否允许任何第三方平台托管的 author upload 进入后续候选。此笔记不做
   candidate freeze、scientific-work equivalence 或最终 gate 判定。
