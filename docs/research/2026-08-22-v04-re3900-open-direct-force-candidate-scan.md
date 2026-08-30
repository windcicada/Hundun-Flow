# v0.4 Re=3900 open direct-force candidate scan

审计日期：2026-08-22（Asia/Shanghai）。这是 RESEARCH_ONLY 的 broader primary-source
scan，不是候选冻结、数值选择、容差选择或最终 ACCEPT/REJECT。目标是寻找 stationary
circular-cylinder cross-flow 在 Re 3000--5000（或明确近邻 3900）处的直接总平均阻力
和/或直接有限跨度总升力 RMS。只采用 publisher、作者/机构 repository、官方政府
记录或官方 OA bytes；不把后续 CFD 汇总、review 数值、pressure-only、base-pressure 或
sectional-only 结果当作目标证据。

## 结论

本轮没有取得一个同时满足“完整 primary apparatus + 目标 Re 的直接总力点 + 可审计跨度、
边界、归一化和不确定度”的候选。因此不能形成生产 reference，也不能支撑 candidate
freeze。最值得合法取得完整原文的两个目标是：

1. Bishop--Hassan (1964)：publisher/Crossref 一手摘要明确写出直接测量 fluctuating
   lift、fluctuating drag 和 steady mean drag，stationary-cylinder Re=3600--11000；
2. Tadrist et al. (1990)：AIP publisher 一手摘要明确写出 global lift coefficient、专门
   的 sensor/signal-conditioning measurement unit，Re=3000--30000，并说明 aspect ratio、
   blockage、ends effects。

Wieselsberger/NACA TN 84 是本轮唯一有官方完整 PDF 且方法级直接总阻力证据的来源，但只
有混合直径的双对数图，没有印刷的 Re=3900 数值表；任何图读数都必须另作受控
DERIVED 工作，不能在本扫描中完成。Fukuoka et al. (2016) 的官方 OA 完整 PDF 反而
提供了一个有用的负证据：PIV 覆盖 Re=4100，但三轴 load-cell 的力测量只在
Re=28000--317000，不能把 PIV Re 当成目标总力。

## 冻结字段和状态语义

每个候选都按以下字段记录：source identity/provenance；total force 与 pressure/sectional
区分；Re target；mean/RMS normalization；active span/end/BC；blockage/correction；
uncertainty；向周期跨度 piD 映射的可能性；ACQUIRE、EXCLUDE 或 FOLLOW-UP。
ACQUIRE 只表示值得取得合法完整 primary bytes，不表示已通过 scientific-work
equivalence；FOLLOW-UP 表示已有足够方法级线索但仍有关键审计缺口；EXCLUDE 表示
本扫描边界下有明确的 Re 或 observable 排除理由。

## 候选矩阵（8 个）

| ID / source identity | Provenance and observable | Re target | Mean / RMS normalization | Active span, end, BC | Blockage / correction / uncertainty | piD mapping possibility | Status |
|---|---|---|---|---|---|---|---|
| C1 Bishop & Hassan (1964), Proc. R. Soc. A 277, 32--50, DOI [10.1098/rspa.1964.0004](https://doi.org/10.1098/rspa.1964.0004) | Royal Society/Crossref 一手摘要：apparatus 直接测 fluctuating lift、fluctuating drag、steady mean drag；圆柱轴线横向于水道流向。不是 pressure-only；但完整 apparatus 未取得。 | 覆盖；摘要给 stationary cylinder Re 3600--11000，包含 3900。 | 只确认 steady mean drag force 和 fluctuating force；系数定义、dynamic-pressure/reference-area/span、RMS 算法和 Re=3900 点均未取得。 | 水道 stationary cylinder 已确认；active span、dummy ends、end plates、support/BC 未核实。 | 摘要未给 blockage、correction 或 uncertainty；不能从二手图表补齐。 | 若完整原文给出 active span 和 force-balance reference area，原则上可检查总力到 piD；当前不可映射。 | ACQUIRE |
| C2 Tadrist, Martin, Tadrist & Seguin (1990), Phys. Fluids A 2, 2176--2182, DOI [10.1063/1.857804](https://doi.org/10.1063/1.857804) | AIP publisher 一手页面明确写 global lift coefficient；unsteady lift 为毫牛量级，专门设计 sensor 与 signal-conditioning circuit；并测上游速度/湍流。摘要没有把 observable 说成 pressure integral，属于直接总力强候选，但仍待完整原文核验。 | 覆盖；AIP 摘要给 Re 3000--30000，3900 在内；全局 lift 在 Re 3000--10000 有显著变化。 | 只确认 global lift coefficient；没有 coefficient denominator、active force area、mean subtraction 或 RMS/statistical window。 | 摘要只说讨论 tube aspect ratio 和 ends effects；span、端板/自由端、支撑和 BC 未给。 | 讨论 blockage，但 blockage 定义、修正和误差预算未给。 | 取得直径、有效 span、传感器受力段和 reference area 后可评估；当前不能映射。 | ACQUIRE |
| C3 C. Wieselsberger (1922), NACA TN 84, [NASA NTRS 19930080855](https://ntrs.nasa.gov/citations/19930080855) | NASA 官方翻译扫描。PDF pp. 3, 6--7 定义 D=cSρV²/2；小圆柱用悬挂摆偏转、较粗圆柱经导线接天平测 drag；是总机械阻力，不是 Cp/base-pressure。 | 宽范围覆盖；正文给 R=Vd/ν=4.2--800000，图上有目标 decade，但没有印刷 Re=3900 表点。 | 系数 c 的总阻力定义明确；机械读数的 time-mean protocol、采样窗口和重复性未报告。 | 作者目标为 infinitely long cylinders；长径比最不利约 380，较粗圆柱夹于与流向平行的板间并以迷宫填塞端部；不是 HUNDUN finite-span piD 试件。 | 试验段尺寸、阻塞比、端部/支撑修正、湍流、粗糙度和 uncertainty 未给；不可从图读数推断。 | 只能作为长/二维总 drag 背景；与有限 piD 的映射目前不支持。 | FOLLOW-UP |
| C4 Richter & Naudascher (1976), JFM 78, 561--576, DOI [10.1017/S0022112076002607](https://doi.org/10.1017/S0022112076002607) | Cambridge publisher 一手摘要：long rigidly supported cylinder in narrow rectangular duct；报告 mean-drag coefficient 与 drag/lift fluctuation RMS。摘要未披露 force balance 是否直接承载整段力，故暂不升级为已核验 direct method。 | 目标未知；publisher 只说 wide Re range around critical value，未列 exact endpoints；不使用 secondary numeric citation。 | mean-drag coefficient、RMS drag/lift 明确存在；reference area、free-stream choice、RMS window 未给。 | long rigid support、duct symmetry 已知；有效 span、端板/漏流和端部 BC 未知。 | 多种 blockage percentages，且文章目标包括 blockage correction；具体 ratio、修正方程和 uncertainty 未取得。 | 完整原文若给 force segment/span，可评估；当前 piD 映射不可审计。 | ACQUIRE |
| C5 Benitz et al. (2016), Computers & Fluids 136, 247--259, DOI [10.1016/j.compfluid.2016.06.013](https://doi.org/10.1016/j.compfluid.2016.06.013)，University of Southampton record [ePrints 507242](https://eprints.soton.ac.uk/507242/) | 机构一手 record 摘要说明 surface-piercing finite cylinders 的 numerical and experimental loads（drag、lift、frequency）；但摘要没有 force sensor/pressure integration 细节，不能认证为 direct total force。 | 近邻但不覆盖；摘要明确 Re=2900、Fr=0.65，低于本轮下界且没有 Re=3900。 | 只说 drag/lift loads 和 drag coefficient comparison；归一化、RMS、参考面积均未给。 | free-surface piercing、free lower end、AR 约 1--19；几何与 BC 明显不同于 fully submerged periodic span。 | blockage、free-surface correction、端部修正、uncertainty 未给；record 无完整附件。 | 近邻只能作为边界/自由表面敏感性线索；不能映射到 piD 总力。 | FOLLOW-UP（近邻） |
| C6 Fukuoka, Hirabayashi & Suzuki (2016), J. Mar. Sci. Technol. 21, 145--153, DOI [10.1007/s00773-015-0338-x](https://doi.org/10.1007/s00773-015-0338-x) | Springer 官方 OA 完整 PDF。PIV 与 force experiment 分开：三轴 load cell 直接测 semi-submerged cylinder 的 total Fx,Fy，但 force Re 不在目标。 | 明确排除：PIV 表 1（PDF p.2 / printed p.146）有 Re=4100--12600；force Table 2（PDF p.5 / printed p.149）为 Re=28000--317000。 | PDF p.4 / printed p.148 给 CD=F̄x/[AR(½ρD²U²)]，CDrms/CLrms 为绕时间均值的 RMS，同一 frontal-area normalization；这部分是 total-force 定义。 | semi-submerged finite cylinder，underwater AR=0.25--1.5，自由表面和 free/end-cell；不是 fully submerged periodic span。 | PIV channel 2 m × 0.3 m × 0.3 m；force towing tank 85 m × 3.5 m × 2.4 m；未给 blockage/correction 或 force uncertainty。 | 即使归一化完整，free-surface/free-end force 且 Re out-of-range，不能映射为 piD target。 | EXCLUDE |
| C7 Gerrard (1961), JFM 11, 244--256, DOI [10.1017/S0022112061000494](https://doi.org/10.1017/S0022112061000494) | Cambridge publisher 一手摘要明确：oscillating lift/drag 由 cylinder-surface fluctuating pressure measurements 确定；因此是 pressure-derived/sectional-style observable，不是本轮要求的 direct total-force authority。 | Re 从 4×10³ 到略高于 10⁵，3900 还低于其声明起点；即使近邻也不满足 direct method。 | 摘要有 RMS lift coefficient，但没有 direct-force normalization；pressure integration/section span 未给。 | active span/end/BC 未给。 | blockage、correction、uncertainty 未给。 | pressure-only source 不作为 piD total-force mapping 输入。 | EXCLUDE |
| C8 Szepessy & Bearman (1992), JFM 234, 191--217, DOI [10.1017/S0022112092000752](https://doi.org/10.1017/S0022112092000752) | Cambridge publisher 一手摘要明确测 fluctuating forces、shedding frequency、spanwise correlation，并改变 moveable end plates；是有限 span force 候选，但不是目标 Re。 | 明确排除；publisher 给 8×10³<Re<1.4×10⁵，高于 3900。 | force observable 和 sectional fluctuating lift trend 明确；完整 coefficient/RMS normalization 未从摘要取得。 | AR=0.25--12；moveable end plates；finite-span/end effects 是研究变量。 | blockage/correction/uncertainty 未给。 | 几何字段较适合后续 piD 比较，但 Re 越界，不能进入目标统计。 | EXCLUDE |

## Primary evidence 与审计缺口

### C1 — Bishop--Hassan

Royal Society landing page 和 Crossref DOI record 是本候选的 provenance。Crossref 的
primary abstract 明确：apparatus 直接测 fluctuating lift/drag 和 steady mean drag；
圆柱轴线垂直于水道流向；stationary-cylinder results cover Re 3600--11000。这足以把
Re=3900 记为“摘要声明范围内”，并把 source 置于 ACQUIRE，但不提供 force balance
图、受力段、端板、CD 定义、RMS 算法、Re=3900 表/图点或 uncertainty。

本环境没有合法取得可读完整论文：Royal Society landing、DOI PDF 和 syndication PDF
均返回访问挑战；没有保存不完整 HTML/PDF 到 literature 目录。因此此项不能填入
数值 receipt，也不能替代完整 primary apparatus。

### C2 — Tadrist et al.

AIP 的 primary article page [Experimental investigation of fluctuating forces exerted on a
cylindrical tube](https://pubs.aip.org/aip/pof/article/2/12/2176/401766/Experimental-investigation-of-fluctuating-forces)
在摘要中明确给出：Re=3000--30000；global lift coefficient；毫牛量级 unsteady lift
requires a dedicated sensor and signal-conditioning circuit；上游 laser-Doppler 测速度场
与湍流；并讨论 aspect ratio、blockage、ends effects。该摘要是本轮最直接的“目标 Re +
直接传感器 + global force coefficient”组合。

但 AIP 页面标注 Available to Purchase，当前未取得合法完整 primary bytes。摘要没有
给出 mean/RMS 的 exact definition、受力 span、reference area、端部装置、blockage
correction、传感器校准或 uncertainty。它应进入后续合法 acquire，而不是当前
scientific-work equivalence。

### C3 — Wieselsberger / NACA TN 84

这是本轮已取得的官方完整 bytes。PDF pp. 3, 6--7 定义总阻力并描述摆/天平的机械测量；
PDF p.16 的 Figure 1 给混合直径双对数 c–Vd/ν 图。正文 broad range 覆盖目标 decade，
但没有逐点表或印刷 Re=3900 数值。该 source 是 DIRECT-METHOD、不是
DIRECT-TARGET(Re=3900)。

本扫描不读取/数字化 Figure 1，不选择图点、不做插值、不设 tolerance。长/无限 cylinder
与 HUNDUN finite periodic span 的差异、blockage、支撑、粗糙度和 uncertainty 都保留
为缺口。

### C4 — Richter--Naudascher

Cambridge publisher primary abstract 直接列出 long rigidly supported cylinder、narrow
rectangular duct、various blockage percentages，以及 mean drag coefficient 和 RMS drag/
lift fluctuations。由于 abstract 没有 exact Re endpoints，也没有 force-balance layout，
本扫描不采用镜像 PDF 或 later paper 的二手数值，不声称它覆盖 3900；只列为需取得
完整原文的 ACQUIRE。

### C5 — Benitz et al.

University of Southampton institutional record 是作者/机构 primary metadata。摘要明确
surface-piercing finite cylinders、free lower ends、AR 约 1--19、Re=2900、Fr=0.65，以及
experimental loads。Re=2900 是近邻负证据，不是 Re=3900；摘要未说明载荷传感器、总力
还是压力/sectional reconstruction，也没有完整附件。因此只作 FOLLOW-UP 近邻，不作
目标数据。

### C6 — Fukuoka et al.

官方 Springer OA PDF 的分离实验是本轮最清楚的排除证据：

- PIV experiment：PDF p.2 / printed p.146，D=26/38/60 mm，Re=4100--12600；这是流场，
  不是 force observable；
- force experiment：PDF p.3 / printed p.147 描述 towing tank 与 three-axis load cell；
  PDF p.5 / printed p.149 的 Table 2 给 force Re=28000--317000；
- PDF p.4 / printed p.148 给 CD、CDrms、CLrms 的 total-force equations。

所以“Re=4100 出现”不能被提升成“Re=4100 direct force”。该 source 明确 EXCLUDE，
但保留为实验分支/observable 分离的审计例子。

### C7 — Gerrard

Cambridge publisher primary abstract 的排除理由是方法而非数值：lift/drag oscillation
由 cylinder-surface fluctuating-pressure measurements 确定；Re 起点为 4×10³。
它既不是 direct force-balance authority，3900 也低于声明起点，故 EXCLUDE。

### C8 — Szepessy--Bearman

Cambridge publisher primary abstract 明确 force measurements、moveable end plates、AR
0.25--12 和 Re 8×10³ 到 1.4×10⁵。它对 finite-span/end-effect 语义有参考价值，但不
覆盖 Re=3900，故 EXCLUDE；不从全文镜像提取数值。

## 已取得完整 primary bytes

仅列本扫描实际依赖且完整、公开合法取得的 primary PDF；其他候选仅保存 publisher/机构
landing URL，没有把 paywall/challenge 响应当 artifact。

### Fukuoka et al. 2016

目录：

~~~
/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/fukuoka-hirabayashi-suzuki-2016-open-20260822/
~~~

| Field | Value |
|---|---|
| URL | https://link.springer.com/content/pdf/10.1007/s00773-015-0338-x.pdf |
| acquired | final response 2026-08-22T09:05:53Z (17:05:53 +08:00), headers in response.headers |
| HTTP / MIME | final HTTP/2 200, application/pdf |
| headers | content-length: 2019725; x-goog-stored-content-length: 2019725; last-modified: Tue, 29 Aug 2017 16:24:17 GMT; ETag "838197a1f0ad36cc4bc216d6792059ac" |
| bytes / SHA-256 | 2,019,725; 560e38492ab89349fa997f209cd5fdd34bc79f24740ceaaeea52c4a06a1933cc |
| pdfinfo | PDF 1.6, 9 pages, 595.276 × 790.866 pt, unencrypted; pdfinfo.txt |
| structure / visual review | `qpdf --check` reports malformed linearization hint-table lengths; `pdftoppm` nevertheless rendered all 9 pages.  The root agent visually reviewed the 3x3 contact sheet and the registered pp. 2--5 locators; render hashes are in `render-review.SHA256SUMS.txt`.  This supports readable-content review, not a claim that the container is structurally clean. |
| locators | PIV Re table PDF p.2/printed 146; force setup p.3/printed 147; force equations p.4/printed 148; force Re Table 2 p.5/printed 149 |

### Wieselsberger / NACA TN 84

该完整官方 artifact 已由同一 lane 的 primary audit 取得，未覆盖或重写：

~~~
/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/wieselsberger-naca-tn84-audit-20260822/NACA-TN-84.pdf
~~~

| Field | Value |
|---|---|
| URL | https://ntrs.nasa.gov/api/citations/19930080855/downloads/19930080855.pdf |
| acquired | 2026-08-22T07:52:48Z--07:52:55Z (headers、时间文件和 manifest 在目录内) |
| HTTP / MIME | final HTTP/2 200, application/pdf |
| bytes / SHA-256 | 12,440,785; 92561c3402f473ded94440ef93229e0d75c6b197e1fcab28206a26b84e35814c |
| pdfinfo | PDF 1.6, 16 pages, 550 × 729 pt, unencrypted; pdfinfo.txt |
| locators | total-drag definition pp.3, 6--7; broad Re range p.6; Figure 1 p.16 |

## Reproducibility commands

以下命令只重取/核验公开 primary metadata 或已取得 bytes；不启动 solver、benchmark、
digitization、long statistics 或 COAST/Stage5 工作。

~~~
# Metadata-only primary checks.
curl -fsSL --max-time 30 'https://api.crossref.org/works/10.1098/rspa.1964.0004' -o /tmp/bishop-crossref.json
curl -fsSL --max-time 30 'https://api.crossref.org/works/10.1063/1.857804' -o /tmp/tadrist-crossref.json

# Verify the two complete PDFs without changing them.
f=/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/fukuoka-hirabayashi-suzuki-2016-open-20260822/fukuoka-2016.pdf
sha256sum "$f"; pdfinfo "$f"; file --mime-type "$f"

n=/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/wieselsberger-naca-tn84-audit-20260822/NACA-TN-84.pdf
sha256sum "$n"; pdfinfo "$n"; file --mime-type "$n"

# Re-acquire Fukuoka only into a new, unique audit directory if bytes are needed.
audit_dir=/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/fukuoka-hirabayashi-suzuki-2016-open-20260822-recheck
mkdir -p "$audit_dir"
curl -L --fail --retry 2 --connect-timeout 20 --max-time 120 -sS \
  -D "$audit_dir/response.headers" \
  -o "$audit_dir/fukuoka-2016.pdf" \
  'https://link.springer.com/content/pdf/10.1007/s00773-015-0338-x.pdf'
sha256sum "$audit_dir/fukuoka-2016.pdf"; pdfinfo "$audit_dir/fukuoka-2016.pdf"
~~~

## Unresolved gaps / handoff boundary

- No candidate in this scan closes the direct Re=3900 total-force gate.
- The legal next input is a complete Bishop--Hassan and/or Tadrist primary paper, with apparatus,
  force segment, coefficient equations, target-Re records, and uncertainty retained as source
  facts. Do not fill missing fields from later reviews.
- Richter--Naudascher needs its complete publisher-author artifact before target Re or direct
  balance status can be classified.
- NACA TN 84 may be eligible for a separately governed controlled figure extraction, but this
  scan intentionally performs no digitization and selects no target value.
- No tolerance, piD equivalence rule, candidate freeze, pressure-performance route, solver,
  COAST equivalence decision, receipt/ledger update, or Stage 5 action was performed.

## Diff boundary

The only repository path intended to change in this lane is this file. External PDF/text/header
files listed above are outside Git and are not solver/reference/receipt/ledger/COAST/Stage5
changes.
