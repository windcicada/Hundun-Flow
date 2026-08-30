# Norberg (2003) finite-span lift primary-source audit

审计日期：2026-08-22（Asia/Shanghai）。本文只记录可复核的 primary-source 证据和
控制性缺口；不选择 HUNDUN 的最终数值、容差、scientific-work equivalence 或
ACCEPT/REJECT。全文检索和下载 artifact 保存在仓库外，避免把不可复核的网页缓存当作
论文原件。

## 结论先行

| 冻结问题 | 当前可审计结论 |
| --- | --- |
| 作者直接测得的是 sectional 还是 finite-span total lift RMS？ | 1987 作者机构 OA thesis 的 Table 2/key 逐字区分 `sect.`（`ell/D -> 0`）与 `total`（whole cylinder），并把 Norberg 1987 present row 标为 `sect., C'_p, FST`；第 12 页说明 present RMS lift 用 RMS pressure 和对置周向点的 cross-correlation 估计。因此 Norberg 自己的 lift evidence 是 sectional pressure-derived，而不是 finite-span total-force RMS。 |
| Re=3900 是否有作者直接值？ | Thesis 第 9--10 页的 Fig.1 是包含 present markers 的 compiled mean-drag plot，Re 轴覆盖下游/低 Re 区间，但没有 Re=3900 行或精确坐标；第 9 页明确 present mean drag 是 mean-pressure integration 的 pressure drag、sectional force，并与 Wieselberger 的 whole-length total drag 区分。因此不能把图中约低 Re 的 present marker 当作 Re=3900 的 direct total mean drag。lift Table 2 的 Norberg present row 仅为 `Re·10^-4=1--30`（`10^4--3×10^5`），没有 Re=3900 lift row。 |
| 能否映射到 HUNDUN periodic `pi D`？ | 不能直接映射。Thesis 给出的 Norberg 1987 lift row 是 `L/D=9,12`、`B=4,11%`、`ell/D=sect.`（`ell/D -> 0`），并说明实验圆柱加端板；没有 `ell/D=pi` whole-cylinder lift、可重建的 spanwise `R_LL(s)` 或 πD active-span measurement。 |
| 完整 primary artifact | 1987 thesis 已从 Lund 官方 OA record 合法取得并逐页核查；2003 article 的完整 PDF 仍未取得（Elsevier API 最小 metadata、ScienceDirect 403、Lund 仅 DOI/portal）。两类 artifact 的路径、headers、时间和 SHA 均在外部目录和本节记录。 |

上述结论的证据等级是“可审计的否定/缺口”，不是对论文不存在某个数值的绝对断言：若取得全文，需逐页复核 Table/Fig/公式后再升级结论。

## 文献身份和 primary 入口

- C. Norberg, “Fluctuating lift on a circular cylinder: review and new measurements,”
  *Journal of Fluids and Structures* **17**(1), 57–96 (2003), DOI
  [10.1016/S0889-9746(02)00099-3](https://doi.org/10.1016/S0889-9746(02)00099-3)。
- 出版方记录：[ScienceDirect article page](https://www.sciencedirect.com/science/article/pii/S0889974602000993)。
  该页是 primary publisher 入口，但本次请求被 Cloudflare 403 拦截，不能作为本地全文 artifact。
- Elsevier 官方 API：[DOI endpoint](https://api.elsevier.com/content/article/doi/10.1016/S0889-9746(02)00099-3)。
  2026-08-22 的匿名 GET 返回 200 和最小 XML metadata；响应明确警告 unauthorized request，且
  `openaccess=0/openaccessArticle=false/openArchiveArticle=false`。
- 作者机构记录：[Lund University Publications](https://lup.lub.lu.se/search/publication/48a4572b-eeaf-4749-9823-a53e5ffe0813)、
  [Lund Research Portal](https://portal.research.lu.se/en/publications/fluctuating-lift-on-a-circular-cylinder-review-and-new-measuremen/)
  和 [Lund publication page](https://www.lu.se/publikation/48a4572b-eeaf-4749-9823-a53e5ffe0813)。三者给出作者、期刊、页码、摘要和 DOI；
  “Links/Final published version”只指向 DOI，没有 author manuscript 或下载文件。

Lund 作者机构摘要可以直接核验以下范围性事实：研究对象是 stationary circular cylinder in cross-flow，
名义情形是 infinitely long、nonconfined、smooth oncoming flow；覆盖约 `Re=47` 到 `2e5`；在
`Re≈0.3e3` 给出约 `30D` 的 spanwise correlation-length 指示，并说 `Re≈1.6e3--20e3`
间 sectional RMS lift coefficient 约增加一个数量级，转变约从 `Re≈5e3` 开始、`Re≈8e3` 完成。
这些是摘要结论，不是 Re=3900 的直接数值。

## 证据矩阵

证据级别：`P0` 为官方 metadata/机构记录，`P1` 为出版方文章页被索引的正文摘录（本次未能保存完整
页面/PDF），`G` 为控制性缺口。所有 `P1` 结论都不得替代逐页全文核验。

| 问题 | 当前 primary 证据 | 状态/限制 |
| --- | --- | --- |
| 研究范围 | Lund 摘要（P0）明确为 stationary circular cylinder、nominally infinitely long/nonconfined、smooth flow，`Re≈47--2e5`。 | 可用作范围边界；不能把名义无限长解释为实验装置真的无限长。 |
| RMS 定义和测量段 | ScienceDirect 正文索引摘录（P1）区分 `ell_c` 段上的 fluctuating lift coefficient：`ell_c/d -> 0` 是 sectional limit，`ell_c=ell` 是 exposed cylinder 全段 total fluctuation。 | 精确符号、分母、RMS/峰值约定需 PDF 逐页核验。 |
| Norberg present measurement | 1987 OA thesis Table 2/key（printed pp.13--14）直接把 Norberg present row 标为 `sect., C'_p, FST`；printed p.12 说明 RMS pressure at `alpha=90°` 与 opposite circumferential-point cross-correlation 的估计方法。2003 出版方摘录给出同一 sectional/near-sectional pressure-method 语境。 | 这支持“Norberg present lift result 为 sectional pressure-derived”，不能排除 2003 review 汇总其他作者的 total-force 行；不能把 review table 混为 Norberg direct total。 |
| Table 1 | 出版方索引摘录（P1）称 Table 1 汇总 `Re<3e5`、`Tu<=2%` 下 sectional 和 total lift fluctuations。 | 表格的作者/实验来源行、`Re=3900` 是否有点、点对应哪种 `ell_c`，全部是 G。 |
| Re=3900 直接值 | Lund 2003 摘要和不可访问的 ScienceDirect PDF 没有 exact `Re=3900` row。1987 OA thesis printed p.9--10（PDF pp.17--18）有 lower-Re compiled mean-drag plot 和 present pressure-derived sectional markers，但无 Re=3900 坐标/表值；其 printed p.13--14 lift row 仅 `Re·10^-4=1--30`。 | 仍不能得到 direct total mean drag 或 total lift RMS at Re=3900；不得采用二手表格或后续 CFD 回引值。 |
| Spanwise correlation | 出版方正文/索引摘录（P1）使用 lift cross-correlation `R_LL`，并区分 one-sided correlation length `Lambda` 与 centroid `sigma`；完全相关时给出 `Lambda=ell`、`sigma=ell/2` 的边界解释。 | 该摘录没有提供 Re=3900 的完整 `R_LL(s)` 数据或精确方程排版；G。 |
| 相关长度的 Re 位置 | 机构摘要（P0）给出约 `30D` 的低 Re 指示；出版方摘录（P1）显示 `Re≈5.1e3` 附近 correlation peak、`Re≈8.0e3` 的 exponential 例及约 `5.3D` 的相关长度。 | 这些不是 `Re=3900` 的可替代数据，也不能外推到 `pi D`。 |
| 仪器/来流 | 出版方实验段摘录（P1）称在 Chalmers L2 low-turbulence wind tunnel；free-stream turbulence `<0.06%`、acoustic noise `<0.6%` dynamic pressure；圆柱水平安装在 vertical supporting plates 之间，工作段约 2.9 m。 | 精确圆柱直径、active span/aspect ratio、板尺寸、测点布置、blockage 比和 tunnel correction 数值是 G。 |
| 壁面/压力修正 | 出版方方法摘录（P1）说明 pressure method neglects fluctuating wall friction，并给出 pressure-RMS 与 total-RMS 的 Re 相关近似比例 `r=1-1.2/sqrt(Re)`（示例 `r≈0.92` at `Re=200`）。 | 公式的精确上下文、是否用于 Norberg present data、Re=3900 的适用性和任何 blockage/end correction 均需全文；不得自行套用。 |
| 不确定度 | 官方 metadata/机构页没有逐测量 uncertainty、重复性、置信区间或 calibration budget。 | G。 |
| 归一化 | 正文摘录表明 `C'_L` 按 dynamic pressure、diameter 和 measured span segment `ell_c` 归一化。 | 精确分母、使用 `rho U_inf^2/2`、sectional limit 的单位长度处理和 RMS 与 amplitude 的区别需 PDF；当前不作数值比较。 |

## 关于 `pi D` 的可推导性（不是 Norberg 的等价结论）

对均匀 span 的 sectional lift covariance，公开的二阶统计恒等式可写成（这里只作为检查
所需数据的数学说明，不声称这是已从全文逐字恢复的 Norberg 方程）：

```text
Var(F_ell) = 2 * integral_0^ell (ell - s) C_LL(s) ds
C'_{L,ell}^2 = (2 / ell^2) * integral_0^ell (ell - s) R_LL(s) * C'_{L,sec}^2 ds
```

其中第二式假定均匀 sectional variance、`R_LL(s)=C_LL(s)/C_LL(0)`，且 `C'_{L,ell}`
按全 span `ell` 归一化。要把它评价在 `ell=pi D`，至少需要：

1. 同一 `Re` 和同一 normalization 下的 `C'_{L,sec}`；
2. `0<=s<=pi D` 的 signed `R_LL(s)` 或足以重建它的 `Lambda/sigma` 数据；
3. active span 和端部/支撑条件，确认实验的 `ell` 与周期计算的 `pi D` 可比；
4. blockage、wall-friction、pressure-to-total 和不确定度处理。

本审计只确认 Norberg 提供了这种“从 sectional correlation 到 finite-span RMS”的方法语境，
没有确认上述四项在 `Re=3900, ell/D=pi` 均有 primary 数值。因此 `pi D` 映射当前是
**未封印的派生路径**，不是直接实验结果，也不构成 COAST/HUNDUN equivalence 判断。

## Artifact、headers、时间和 hash

外部 artifact 根目录：

```text
/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/norberg2003-20260822-audit/
```

时间记录在 `audit-start-time.txt`（`2026-08-22T16:12:43+08:00`），每个文件的 SHA-256
也收录在同目录 `SHA256SUMS.txt`。表中的“paper pages”仅是官方 metadata 的页码范围
`57--96`；由于没有 PDF，**PDF page count = unknown**。

| 文件 | URL/用途 | HTTP / MIME | bytes | SHA-256 | 结论 |
| --- | --- | --- | ---: | --- | --- |
| `publisher-api.xml` | Elsevier DOI API | 200 / `text/xml;charset=UTF-8` | 1830 | `4ccc897d9bbe90d0160015355c42055aeb371b3f0d3d09568c78d62ece53fb7e` | 官方最小 metadata；不是全文 |
| `publisher-api.headers` | 同上 headers | 200 | 1058 | `f7810e229052270d997371d9ae84e2280abd8804a40431fd9d9a25838fe0c3c2` | 含时间、`openaccess=0` 警告等响应上下文 |
| `publisher-api-full.response` | API `?view=FULL` | 401 / `text/xml;charset=UTF-8` | 133 | `a199c5505f5d84ca6206475b3a116d509aa4528f154f0d3cdbc57669ec400bb5` | 未授权，非全文 |
| `publisher-api-pdf.response` | API `?httpAccept=application/pdf` | 406 / `text/xml` | 162 | `8896659c8776e3601cba7215e75d52d13a0aa4cde5150bbf00660945a33f39b7` | 不接受，非 PDF |
| `lund-record.html` | Lund institutional record | 200 / `text/html; charset=utf-8` | 37756 | `6dcd23df70e85b1243abd4c28f7ec84ff4f709f636e676bb4f4e557104bc8ddd` | metadata/abstract；无文件 |
| `lund-portal.html` | Lund Research Portal | 200 / `text/html;charset=UTF-8` | 48293 | `d5d4f57a3eca74df0c4056b49ef168f3cbdf58260299d74891efd761d18f5eac` | metadata/abstract；无文件 |
| `lund-page.html` | Lund public page | 200 / `text/html; charset=UTF-8` | 112776 | `342d236eb4d47c9661e5a5bfb06b6f001ad3e2f6e20e31c08cc177b956919e1a` | 明列 final published version 为 DOI；无文件 |
| `sciencedirect-pdfft.response` | publisher PDF endpoint `/pdfft` | 403 / `text/html; charset=UTF-8` | 832805 | `87adc5c87ea8952443fc8789c60eb01f76e22947d4af83c0879f32fb0d646b73` | Cloudflare challenge/HTML；明确不是 PDF |
| `sciencedirect-pdfft.headers` | 同上 headers | 403 | 936 | `87dcfc0cd453ef3f4aeb1c343eb7387828cf0e0b2a1e05dfc07f8b6e466358a1` | 访问证据 |

曾作为检索导航发现一个非 publisher/author/institution 的 mirror URL：

```text
https://www.electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20Fluids%20and%20Structures/2003%20Volume%2017/1/57-96.pdf
```

它未被当作 primary artifact；2026-08-22 的 HTTPS `curl` 失败，exit code 35
（`OpenSSL SSL_connect: SSL_ERROR_SYSCALL`）。`mirror-attempt.log`、`mirror.headers` 和
`mirror.stderr` 只保存失败证据，**没有 mirror PDF bytes、MIME、page count 或 paper SHA**。

## 可复核命令

以下是本次实际使用的下载/访问命令（artifact 根目录如上；所有写入均在仓库外）：

```bash
AUDIT_DIR=/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/norberg2003-20260822-audit
mkdir -p "$AUDIT_DIR"
date -Is > "$AUDIT_DIR/audit-start-time.txt"

curl -L --http1.1 --silent --show-error --dump-header "$AUDIT_DIR/publisher-api.headers" \
  --output "$AUDIT_DIR/publisher-api.xml" \
  'https://api.elsevier.com/content/article/doi/10.1016/S0889-9746(02)00099-3'
curl -L --http1.1 --silent --show-error --dump-header "$AUDIT_DIR/publisher-api-full.headers" \
  --output "$AUDIT_DIR/publisher-api-full.response" \
  'https://api.elsevier.com/content/article/doi/10.1016/S0889-9746(02)00099-3?view=FULL'
curl -L --http1.1 --silent --show-error --dump-header "$AUDIT_DIR/publisher-api-pdf.headers" \
  --output "$AUDIT_DIR/publisher-api-pdf.response" \
  'https://api.elsevier.com/content/article/doi/10.1016/S0889-9746(02)00099-3?httpAccept=application/pdf'

curl -L --http1.1 --silent --show-error --dump-header "$AUDIT_DIR/lund-record.headers" \
  --output "$AUDIT_DIR/lund-record.html" \
  'https://lup.lub.lu.se/search/publication/48a4572b-eeaf-4749-9823-a53e5ffe0813'
curl -L --http1.1 --silent --show-error --dump-header "$AUDIT_DIR/lund-portal.headers" \
  --output "$AUDIT_DIR/lund-portal.html" \
  'https://portal.research.lu.se/en/publications/fluctuating-lift-on-a-circular-cylinder-review-and-new-measuremen/'
curl -L --http1.1 --silent --show-error --dump-header "$AUDIT_DIR/lund-page.headers" \
  --output "$AUDIT_DIR/lund-page.html" \
  'https://www.lu.se/publikation/48a4572b-eeaf-4749-9823-a53e5ffe0813'

curl -L --http1.1 --silent --show-error --dump-header "$AUDIT_DIR/sciencedirect-pdfft.headers" \
  --output "$AUDIT_DIR/sciencedirect-pdfft.response" \
  'https://www.sciencedirect.com/science/article/pii/S0889974602000993/pdfft'

MIRROR_URL='https://www.electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20Fluids%20and%20Structures/2003%20Volume%2017/1/57-96.pdf'
curl -L --http1.1 --connect-timeout 20 --max-time 60 --silent --show-error \
  --dump-header "$AUDIT_DIR/mirror.headers" --output /dev/null "$MIRROR_URL"
# 实际返回：exit 35，SSL_ERROR_SYSCALL；stderr 记录于 mirror.stderr。

find "$AUDIT_DIR" -maxdepth 1 -type f ! -name SHA256SUMS.txt -print0 \
  | sort -z | xargs -0 sha256sum > /tmp/norberg2003-SHA256SUMS.txt
mv /tmp/norberg2003-SHA256SUMS.txt "$AUDIT_DIR/SHA256SUMS.txt"
```

## Controlling gaps / hand-off boundary

在以下证据补齐前，Norberg receipt 不能标记为 complete，也不能用它封印 HUNDUN 的
`Re=3900, span=pi D` 统计 gate：

1. 合法取得并 hash 完整 primary PDF 或作者机构 manuscript；
2. 逐页核对 §2/§3 的 `C'_L`、`R_LL`、`Lambda`、`sigma` 精确公式、Table 1、相关图和实验 apparatus；
3. 证明 `Re=3900` 的候选数值是 Norberg 直接 measured 值，而不是 review 中转引的其他作者数据；
4. 获得同一 normalization 下的 sectional RMS、`R_LL(s)`/correlation-length 数据、active span 和 `ell/D`；
5. 核验 blockage、端板/支撑、wall-friction/pressure correction、force calibration、采样/统计长度和 uncertainty。

这些是 evidence gaps，不是本 worker 对 HUNDUN/COAST 等价规则或最终容差的判断。

## Thesis follow-up：Norberg (1987) 官方 OA PDF

### 来源、artifact 和逐页核查状态

本 follow-up 使用 Lund 官方记录 [Norberg 1987 thesis record](https://lup.lub.lu.se/search/publication/4a9470ca-10cf-4bb8-9f80-9ea6c000642d)
及其明确标为 Open Access、Publisher's PDF 的
[Norberg_PhD_1987.pdf](https://lup.lub.lu.se/search/files/137948976/Norberg_PhD_1987.pdf)。记录给出的正式 thesis
页数为 28；下载 PDF 包含 8 页 Lund/版权/封面前置页，因此 `pdfinfo` 总页数为 36，打印的 thesis 正文页为 1--28。

artifact 外部目录仍为：

```text
/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/norberg2003-20260822-audit/
```

下载时间、headers、原始 PDF、`pdfinfo`、`qpdf --check`、逐页核查记录和独立 hash manifest 均在该目录。关键结果：

| artifact | bytes / MIME | SHA-256 | 备注 |
| --- | ---: | --- | --- |
| `norberg1987-phd-thesis.pdf` | 7,853,503 / `application/pdf` | `beb7b87242ac3d3693e6f9dd1990b38d289e8b87cd4892ed011c9718c7d37495` | `pdfinfo`: PDF 1.5, 36 pages, A4, unencrypted |
| `norberg1987-thesis.headers` | 262 bytes; 200 response, `application/pdf`, chunked | `d2261ffb2daf67b697f58c3ade01f30f43cfce45ede2a2f1e61844a4a2028eb6` | `Content-Disposition: Norberg-PhD-1987.pdf`; Date `Sat, 22 Aug 2026 08:22:51 GMT` |
| `norberg1987-thesis.pdfinfo` | 332 bytes | `0f4ef2216f66d5dbd9014126d6e731c0c8abfbfd2df3c6de55e04920f72bde70` | records page count/size/version |
| `norberg1987-thesis.qpdf-check` | 292 bytes | `65f0a0fdae91d70f6125f404cd89e213d73a4c832497a3eedb6fccff2e34b0a7` | no syntax/stream encoding errors reported |
| `norberg1987-record.html` | 32,518 / `text/html; charset=utf-8` | `bf81c20f3fc1c92d430c364a170cb621a4c75e23a413f69db7b3b194b5791ac5` | record confirms 28 pages, Open Access PDF and ISBN |
| `norberg1987-record.headers` | 362 bytes; 200 response, `text/html; charset=utf-8` | `9818a057e6090dab8914ee3d84c80fc2a16783a716d266cfeb7d0fa0513a52a7` | Date `Sat, 22 Aug 2026 08:29:43 GMT` |
| `norberg1987-page-review.txt` | page-by-page locator log | `be30f7dedbc7670920be164716bfbf9b5187d029cc9bd31fdfecc71948e544cf` | PDF pages 1--36 individually reviewed |
| `norberg1987-thesis.SHA256SUMS.txt` | hash manifest | — | hashes above plus request logs and qpdf result |

`pdftotext` 只能读出 Lund wrapper，正文页是扫描图；因此每一页均以 100 dpi raster 逐页视觉核查，
并在 `thesis-pages-100dpi/page-01.jpg` ... `page-36.jpg` 留有可复核渲染。核心论文页又以 200 dpi
复核了 Table 2/key、Fig. 1、Fig. 4、pressure-distribution 图和相关文字。PDF 服务器请求时间为
`2026-08-22T16:22:49+08:00`--`16:23:03+08:00`，record 请求时间为
`2026-08-22T16:29:41+08:00`--`16:29:44+08:00`。

本 follow-up 的完整增量 diff 另保存为
`/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/norberg2003-20260822-audit/norberg2003-audit-followup.diff`；其当前 SHA-256 记录在
`norberg1987-thesis.SHA256SUMS.txt`。

本 follow-up 的下载与完整性核验命令为：

```bash
AUDIT_DIR=/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/norberg2003-20260822-audit
THESIS_URL='https://lup.lub.lu.se/search/files/137948976/Norberg_PhD_1987.pdf'
RECORD_URL='https://lup.lub.lu.se/search/publication/4a9470ca-10cf-4bb8-9f80-9ea6c000642d'

curl -L --http1.1 --silent --show-error \
  --dump-header "$AUDIT_DIR/norberg1987-thesis.headers" \
  --output "$AUDIT_DIR/norberg1987-phd-thesis.pdf" "$THESIS_URL"
curl -L --http1.1 --silent --show-error \
  --dump-header "$AUDIT_DIR/norberg1987-record.headers" \
  --output "$AUDIT_DIR/norberg1987-record.html" "$RECORD_URL"
pdfinfo "$AUDIT_DIR/norberg1987-phd-thesis.pdf" \
  > "$AUDIT_DIR/norberg1987-thesis.pdfinfo"
qpdf --check "$AUDIT_DIR/norberg1987-phd-thesis.pdf" \
  > "$AUDIT_DIR/norberg1987-thesis.qpdf-check" 2>&1
pdftoppm -r 100 -jpeg -jpegopt quality=85 -f 1 -l 36 \
  "$AUDIT_DIR/norberg1987-phd-thesis.pdf" \
  "$AUDIT_DIR/thesis-pages-100dpi/page"
```

### Thesis 对冻结问题的直接回答

| 冻结问题 | thesis primary 证据和结论 |
| --- | --- |
| Re≈3900 direct total mean drag？ | **没有 direct total mean drag。** printed p.9（PDF p.17）说 present mean drag coefficients 来自 mean pressure distributions 的 integration，并明确 present data 是 “only pressure drag” 的 sectional forces；它把 Wieselberger 的 whole-length total drag（pressure + friction）作为对照。printed p.10（PDF p.18）Fig.1 是 compiled mean-drag plot，曲线图例包含 Wieselberger/Tritton/Roshko/Schewe，另叠加 `Present Tu=0.1%` 和 `Present Tu=1.4%` markers；这些 present markers 不能升级为 whole-cylinder total-force measurement。图轴覆盖低 Re，但没有 Re=3900 的行、坐标标注或可复核表值。 |
| Re≈3900 finite-span/total lift RMS？ | **没有。** printed p.13--14（PDF p.21--22）Table 2/key 将 `sect.` 定义为 `ell/D -> 0`、`total` 定义为 whole-cylinder force；Norberg 1987 present row 是 `Re·10^-4=1--30`、`L/D=9,12`、`B=4,11%`、`ell/D=sect.`、remark `C'_p, FST`。因此 present lift RMS 既不是 total，也不覆盖 `Re=0.39·10^4`。Fig.4（PDF p.22）present curves likewise begin around `Re~10^4`, not 3900. |
| 测量方法 | printed p.7（PDF p.15）说明 P1/P2 用 pinhole microphones 测 wall pressure，并以两点 wall-pressure 的 phase/coherence 估计 sectional forces；printed p.12（PDF p.20）说明 present RMS lift coefficients 用 `RMS pressure coefficient` at circumferential `alpha=90°` 与 opposite points `alpha=-90°/90°` 的 cross-correlation 估计。printed p.14（PDF p.22）明确 `C_L` 是 force measurement、`C_LP` 是 pressure integration、`C'_P` 是 RMS-pressure lift。没有 whole-span force balance 记录。 |
| 归一化 | thesis summary/key 只定义 `C_L`/`C_LP`/`C'_P` 的测量来源和 RMS 标记，没有给出完整 `rho U_inf^2/2` 分母、时间 RMS/幅值转换或有限段归一化式；所引用 P3/P4 原论文不嵌入此 thesis。故不能从 thesis 单独封印与 HUNDUN coefficient 的数值归一化相同。 |
| active span / aspect ratio | Table 2 present P3 给 `L/D=9,12`（aspect ratio），但 active force-sensing `ell/D` 明确为 `sect.`，即 `ell/D -> 0`；这不是 `ell/D=pi` 的 total span。 |
| 端板和边界 | printed p.6（PDF p.14）说本工作圆柱尽可能刚性、光滑，并加 end plates 以改善 end conditions。该摘要没有给端板尺寸、端部间隙或端部修正。 |
| tunnel / blockage | printed p.7（PDF p.15）P1/P2 tunnel cross-section `0.5 x 0.4 m^2`，`D=41 mm`，reference `Tu<0.1%`，测试 `Tu=1.3--4.1%`、`Lambda/D=0.1--0.5`；printed p.8（PDF p.16）P3/P4 tunnel cross-section `1.25 x 1.80 m^2`，P3 `D=41,120 mm`、`Re~2e4--3e5`。Table 2 的 P3 `B=4,11%`。printed p.9/p.18/p.19（PDF p.17/p.26/p.27）说明 mean coefficients 按 Allen--Vincenti blockage method correction；printed p.12 footnote（PDF p.20）说明较大 `D=120 mm, B=11%` 系列作 blockage correction 以接到 `D=41 mm, B=4%` 系列。精确 correction equation/uncertainty 未给出。 |
| uncertainty | printed p.9（PDF p.17）说图中曲线为压低 reported-data scatter 而平滑；printed p.12（PDF p.20）称 present RMS results 尤其在 critical regime 应视为 tentative；没有 error bars、重复性、confidence interval 或 calibration budget。 |
| spanwise correlation / `pi D` | printed p.8（PDF p.16）只说 axial correlations 在 `Re~5e3`（有 turbulence 时约 `4e3`）附近有 local maximum，随后随 Re 快速降低；没有给 `R_LL(s)` 方程、数据表或可积到 `pi D` 的曲线。p.12 的 opposite-point cross-correlation 是圆周压力点相关，不是 spanwise `R_LL(s)`。因此 thesis 不能提供 `ell=pi D` 的 total RMS 映射。 |

### 逐页 locators（PDF 页 → printed thesis 页）

下面列出本次逐页检查的定位，不把扫描页中未出现的内容推断为缺失实验数据；完整逐页日志和 raster 文件位于外部目录。

| PDF page | printed page / 内容定位 |
| ---: | --- |
| 1 | Lund publisher wrapper、rights、title/year/author；非 thesis data |
| 2 | Lund download-date sheet；无 data |
| 3 | Chalmers title cover；无 data |
| 4 | blank separator |
| 5 | dissertation title/submission/opponent page；无 data |
| 6 | ISBN/printer page；无 data |
| 7 | internal title page；无 data |
| 8 | blank separator |
| 9 | contents；无 data |
| 10 | printed p.2 abstract：experimental scope、mean/RMS pressure、near-wake、transition around `Re~5e3` (`~4e3` with turbulence) |
| 11 | printed p.3 abstract continuation：subcritical transition/pressure-force summary；无 total-force result |
| 12 | printed p.4 dissertation：P1--P4 list |
| 13 | printed p.5 introduction：Re definition and flow regimes |
| 14 | printed p.6 Table 1：aspect ratio/blockage/end conditions；rigid smooth cylinders with end plates |
| 15 | printed p.7 Papers 1/2：`0.5x0.4 m^2` tunnel, `D=41 mm`, microphone sectional-force method |
| 16 | printed p.8 Papers 3/4：`1.25x1.80 m^2` tunnel, P3 diameters/range, P4 down to `Re~50`, axial-correlation trend |
| 17 | printed p.9 summary：mean drag from pressure integration; blockage correction; present sectional pressure drag vs Wieselberger whole-length total |
| 18 | printed p.10 Fig.1/2：compiled mean drag/base pressure; present markers and background legend; no exact Re=3900 total row |
| 19 | printed p.11 Fig.3：St versus Re and transition discussion |
| 20 | printed p.12 lift method/limitations: RMS `C'_P` formula inputs, tentative results, blockage footnote |
| 21 | printed p.13 Table 2: present 1987 P3 `Re·10^-4=1--30`, `L/D=9,12`, `B=4,11%`, `ell/D=sect.` |
| 22 | printed p.14 Table 2 key/Fig.4: `sect.=ell/D->0`, `total=whole cylinder`, `C_L` force, `C_LP` pressure integration, `C'_P` RMS pressure; present curve starts near `Re~1e4` |
| 23 | printed p.15 Fig.5 transition frequency and near-wake lengths |
| 24 | printed p.16 Fig.6 diffusion-length scaling |
| 25 | printed p.17 Fig.7 vortex-formation/transition lengths; scatter/need for further measurements |
| 26 | printed p.18 mean/RMS pressure-distribution discussion; mean blockage correction |
| 27 | printed p.19 Fig.8 mean/RMS pressure distributions; mean coefficients corrected for blockage |
| 28 | printed p.20 Fig.9 RMS pressure vs turbulence intensity/scale; sectional lift relation |
| 29 | printed p.21 Fig.10 critical-regime pressure distributions |
| 30 | printed p.22 conditional sampling: pinhole wall pressure and conditional pressure/shear |
| 31 | printed p.23 Fig.11 conditional pressure/force-vector sketch; no span-integrated force |
| 32 | printed p.24 acknowledgements/references 1--8 |
| 33 | printed p.25 references 9--22 |
| 34 | printed p.26 references 23--35 |
| 35 | printed p.27 references 36--49 |
| 36 | printed p.28 references 50--53 |

### Follow-up conclusion and remaining control gaps

This OA thesis tightens, rather than overturns, the 2003 audit:

1. Norberg 1987 **does provide author-primary pressure-derived sectional mean-drag and sectional RMS-lift evidence**, but the thesis explicitly distinguishes it from whole-cylinder total force. The lower-Re Fig.1 compiled plot therefore cannot be used as a direct total mean-drag receipt at Re=3900.
2. The author-primary lift summary has no Re=3900 row: the present row begins at `10^4` and is `sect., C'_p`, with `L/D=9,12`, blockage `4/11%`, and end plates. It cannot supply finite-span/total RMS lift at `ell=pi D`.
3. The thesis gives no exact coefficient normalization equation, spanwise covariance function, active `ell/D=pi` measurement, end-plate correction, uncertainty budget, or enough signed axial correlation data for a `pi D` derivation.

Accordingly the thesis is a complete primary artifact for the method/definition audit, but it does **not** seal the Re=3900 finite-span total-lift or total-drag evidence needed by the parent gate. This follow-up makes no HUNDUN tolerance, COAST equivalence, candidate, or ACCEPT/REJECT decision.
