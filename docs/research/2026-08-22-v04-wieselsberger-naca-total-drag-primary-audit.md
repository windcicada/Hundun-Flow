# NACA TN 84 Wieselsberger total-drag primary audit

审计日期：2026-08-22（UTC）

审计范围：只核对 C. Wieselsberger, *New Data on the Laws of Fluid Resistance*, NACA Technical Note No. 84 (1922) 是否能为 "total_mean_drag_coefficient_re3900" 提供公开、可审计的 primary authority。本文不选取 Re=3900 的数值，不设计容差，不决定 ACCEPT，也不修改 solver、receipt、ledger、candidate 或 Stage 5 内容。

## 结论摘要

NASA NTRS 保存的是一份公开可下载的 NACA 翻译扫描件。原始方法层面具备直接总阻力测量证据：文中把阻力 \(D\) 定义为

\[
D = c S \rho V^2/2
\]

并说明小圆柱由悬挂摆的偏转测得阻力，较粗圆柱通过导线连接到天平测得阻力（PDF pp. 3, 6--7；打印页码 "-2-", "-5-", "-6-"）。这不是压力孔积分或仅有 \(C_p\) 的结果，因此可标为 "DIRECT-METHOD"。

但是，"DIRECT-TARGET(Re=3900)" 尚未成立。报告正文给出实验范围 \(R=Vd/\nu=4.2\)--\(800{,}000\)，并在 Figure 1（扫描 PDF p. 16）以双对数坐标绘制所有结果；图上没有逐点表、没有 Re=3900 的印刷数值，也没有单点的不确定度。因而 Re=3900 的任何数值都必须从扫描图控制数字化/插值取得，只能先标为 "DERIVED"。该 source 足以开启受控数字化工作项，不足以在当前状态封印 direct target。

更重要的可比性缺口是：主实验意在“infinitely long cylinders”，长径比最不利约 380；较粗圆柱夹在与气流平行的平板之间并用迷宫式填塞防漏，但没有给出试验段尺寸、阻塞比、端部修正、来流湍流、表面粗糙度、支撑/导线修正或测量不确定度（PDF pp. 5--7, 13）。这些字段必须在后续 candidate metadata 中保留为缺失，不能用图表数字化掩盖。

## 1. Primary artifact、公开性和完整性

### 1.1 官方记录

- 官方记录：[NASA NTRS citation 19930080855](https://ntrs.nasa.gov/citations/19930080855)
- 官方元数据 API：[NASA NTRS citation API](https://ntrs.nasa.gov/api/citations/19930080855)
- 官方 PDF：[19930080855.pdf](https://ntrs.nasa.gov/api/citations/19930080855/downloads/19930080855.pdf)
- 记录字段：title "New Data on the Laws of Fluid Resistance"；author "Wieselsberger, C"，Aerodynamic Institute, Gottingen；publication date "1922-03-01"；report number "NACA-TN-84"；distribution "PUBLIC"；copyright determination "GOV_PUBLIC_USE_PERMITTED"；"downloadsAvailable=true"。
- PDF 扉页（PDF p. 2）明确写 "TECHNICAL NOTE No. 84"、题名和 "By C. Wieselsberger"。结尾（PDF p. 13，打印页 "-12-"）写 "Aerodynamic Institute, Göttingen, April, 1921" 及 "Translated by National Advisory Committee for Aeronautics"。这里保留 NTRS 的 1922 publication date 与原文机构日期，不将两者混为同一日期。

### 1.2 本次独立取得的字节 artifact

本次从 NASA 官方 PDF URL 重新下载，没有登录、绕过访问控制或使用第三方镜像。为避免覆盖已有工作，保存到新的外部目录：

~~~text
/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/wieselsberger-naca-tn84-audit-20260822/
~~~

| 项目 | 记录 |
|---|---|
| PDF | "NACA-TN-84.pdf" |
| 获取开始/结束 | "2026-08-22T07:52:48Z" / "2026-08-22T07:52:55Z" |
| HTTP | final response "HTTP/2 200" |
| MIME | "application/pdf" |
| Content-Length | "12440785" bytes |
| Last-Modified | "Thu, 06 Aug 2020 20:18:04 GMT" |
| ETag | "fa01da400b196f10d64b2a7b934321ab" |
| SHA-256 | "92561c3402f473ded94440ef93229e0d75c6b197e1fcab28206a26b84e35814c" |
| PDF pages | 16 |
| page size / rotation | "550 x 729 pt", rotation "0" |
| PDF version / producer | PDF 1.6; Adobe Acrobat 9.54 Paper Capture Plug-in |
| encryption | none |

外部目录还保存了 "response.headers"、"ntrs-citation.json"、"ntrs-fulltext.txt"、"pdfinfo.txt"、"file-mime.txt" 和 180 dpi 全 16 页渲染及其 "SHA256SUMS"。NASA 的 fulltext endpoint 只有 18 字节、不能提供可依赖的 OCR；正文判断以原始扫描页逐页目视核对为准。

## 2. 逐页 primary evidence

下表的 PDF page 是下载文件的 1-based page；扫描页内的 "-2-" 等是原报告打印页码。

| Locator | source observation | 对 HUNDUN target 的含义 |
|---|---|---|
| PDF p. 2, title page | NACA TN 84；题名和作者；正文说明流体粘性使阻力需由实验确定 | primary bibliographic identity；不是数值表 |
| PDF p. 3, printed "-2-" | 对横向来流圆柱给出 \(D=cS\rho V^2/2\)；\(D\) 是 drag，\(S\) 为投影面积（或选定 characteristic area），\(c\) 为 dimensionless coefficient of drag | 明确定义总阻力系数语义；现代 "Cd" 对 "c" 的名称映射仍应在 candidate metadata 说明 |
| PDF p. 4, printed "-3-" | \(d\) 为线性 characteristic dimension，\(\nu=\mu/\rho\)，\(R=Vd/\nu\)；对圆柱后文采用直径为 \(d\) | Re 归一化可审计 |
| PDF p. 5, printed "-4-" | 圆柱轴与气流正交；所有测量为 "uniplanar flow"，流线形态在与轴平行的平面中相同；作者因此说所得系数是 infinitely long cylinders 的系数 | 不是 HUNDUN 的周期跨度 \(\pi D\) finite-span 试件；只能作为长/二维总阻力背景 |
| PDF p. 5--6, printed "-4-"/"-5-" | 速度约 1.2--36 m/s；9 个圆柱直径 0.05--300 mm | 说明通过尺寸/速度改变 Re；未给出 Re=3900 单点表 |
| PDF p. 6, printed "-5-" | 这些实验覆盖 \(R=4.2\)--\(800{,}000\)；直径是 characteristic length；直径不超过 8 mm 的小圆柱垂直悬挂在长线上，底部带重物形成摆，以气流下的偏转确定 drag | direct mechanical force method；目标 Re 落在总体范围内 |
| PDF p. 6--7, printed "-5-"/"-6-" | 圆柱长度远大于直径，最不利约 380 倍；较粗圆柱置于气流内、与气流平行的两块刚性平板之间，端部与平板之间有 labyrinth packing；导线接到 balance 测 drag | span/end boundary 有定性说明，但没有平板尺寸、阻塞比或 correction |
| PDF p. 7, printed "-6-" | Figure 1 将全部实验的 \(c\) 对 \(R=Vd/\nu\) 作 logarithmically divided coordinates；实验值延伸到约 4.2 | 目标 Re 附近有图表覆盖；没有 direct tabulated Re=3900 value |
| PDF p. 8, printed "-7-" | 约 \(R=2000\) 有明显向下偏离；\(R=15{,}000\)--\(180{,}000\) 时 \(c\) 约 1.2；约 \(R=300{,}000\) 发生快速下降 | 3900 位于曲线的低 Re 区间；不能据此替代精确 target extraction |
| PDF p. 9, printed "-8-" | 对某一 30 cm 圆柱，速度 15--30 m/s 时绝对 drag 约从 4 kg 降到 3.5 kg；并讨论不同直径曲线覆盖 | 进一步表明测的是绝对力；不是 Re=3900 表格 |
| PDF p. 13, printed "-12-" | 讨论人为粗糙度改变分离点和阻力系数；Figure 7 是另一种后部锥削圆柱的 smooth/rough 对比 | 主 Figure 1 圆柱的粗糙度没有定量报告；surface condition 是 candidate gap |
| PDF p. 16, Figure 1 | 图例列出直径 0.05、0.1、0.3、1.0、3.0、7.9、42.0、80.0、300 mm；横轴是 \(Vd/\nu\) 的 log scale，纵轴为 \(c\) 的 log scale；扫描内容需旋转后阅读 | 具有可复核的图例和双对数轴，可开启受控 digitization；图上点/符号拥挤，不能直接读作表值 |

## 3. Direct total drag、Re=3900 和“mean”的拆分

### 3.1 "DIRECT-METHOD": 方法级通过

这份 primary source 没有以压力孔的周向积分来构造 "C_Dp"。它直接通过：

1. 小圆柱悬挂系统的摆偏转；
2. 较粗圆柱经导线连接的 balance；

测量圆柱受到的 drag，再按 \(D=cS\rho V^2/2\) 定义系数。只要把 "drag" 理解为沿来流方向的总流体力，该量包含压力和表面摩擦的合力，方法级可标为 "DIRECT-METHOD"。这与只报告 "C_p"、base pressure 或 pressure drag 的来源不同。

但原文没有现代实验报告所称的时间序列、采样窗口或 "mean" 算法。机械摆/天平的静态偏转可合理理解为稳态或时间平均力的读数，却不能在没有额外证据时写成已审计的 "mean" protocol。因此应记录：

~~~text
direct total-drag force: supported
modern time-mean protocol: not reported
~~~

### 3.2 "DIRECT-TARGET(Re=3900)": 未成立

以下事实只证明目标在 source 的总体范围内，不证明有 Re=3900 的直接点：

- PDF p. 6 写实验 \(R\) 范围从 4.2 到 800,000；
- PDF p. 7 写 Figure 1 的实验值延伸至约 4.2；
- Figure 1 只是一张多直径、多点的双对数扫描图，没有逐点表或每个 marker 的 "R", "c", "d" 数字列表。

因此本审计不输出 Re=3900 的 "c"/"Cd" 数值，也不把图中目测值写入 direct target。若后续从 Figure 1 取得某个目标点或邻近点，它必须标为 "DERIVED-FROM-FIGURE"，并另行说明是否做了同一直径的 interpolation；不能把图读数升级为 "DIRECT-TARGET"。

## 4. Geometry / boundary / blockage / correction / uncertainty matrix

| 字段 | primary evidence | 当前状态 |
|---|---|---|
| body orientation | cylinder axis normal to stream（PDF p. 5, printed "-4-"） | "SUPPORTED" |
| diameter and span intent | 0.05--300 mm；长径比最不利约 380；作者明确称 infinitely long cylinders（PDF pp. 5--6） | "SUPPORTED, but not HUNDUN finite-span equivalent" |
| flow idealization | "uniplanar flow"；同一平面族内流线相同（PDF p. 5） | "SUPPORTED, 2-D/long-cylinder context" |
| small-body support | long wire + lower weight, pendulum deflection (PDF p. 6) | "SUPPORTED" |
| large-body support | wires to a balance; two rigid walls parallel to flow; labyrinth packing at ends (PDF pp. 6--7) | "SUPPORTED qualitatively" |
| fluid / kinematic viscosity | footnote gives air \(\nu=0.145\ {\rm cm^2/s}\) at 760 mm pressure and 15°C (PDF p. 6) | "PARTIAL"; density/temperature/velocity calibration per point absent |
| blockage / test-section dimensions | no width, height, open-area ratio, wall spacing or blockage correction found in 16-page artifact | "MISSING" |
| end/support interference correction | qualitative assertion that edge nonuniformity is confined to a few diameters; no measured correction or uncertainty | "MISSING" |
| surface roughness | roughness is discussed for a separate tapered-cylinder illustration (PDF p. 13); main Figure 1 specimen roughness is not quantified | "MISSING" |
| inlet turbulence / freestream quality | no turbulence intensity or boundary-layer specification found | "MISSING" |
| force calibration / tare / wire drag | balance/pendulum principle is stated; calibration, tare, support drag and repeatability are not reported | "MISSING" |
| uncertainty | no error bars, standard deviations, confidence interval or uncertainty budget found | "MISSING" |
| compressibility | air and velocities are stated, but no Mach number or compressibility correction is reported | "MISSING"; likely low-speed, not an accepted equivalence claim |

The phrase in PDF p. 6 that edge effects “could not materially affect the main flow” is the author’s qualitative argument based on the long cylinder; it is not a numeric blockage or end correction. It must not be converted into a HUNDUN equivalence tolerance.

## 5. Figure 1 controlled-digitization readiness

### 5.1 What is recoverable

Figure 1 (PDF p. 16) has a visible logarithmic grid, a legend keyed to the nine diameters, and labels \(R=Vd/\nu\) and \(c\). It is therefore technically possible to prepare a controlled digitization artifact. The source scan is the official 16-page PDF above, not a screenshot or a secondary plot.

The chart page is embedded in a portrait PDF page but the printed graph is rotated in the scan. A future extractor must record the exact rendering command, image hash, rotation, crop and the source PDF hash. The following are controls, not a selected data value:

- calibrate both axes from the actual printed grid/tick centers in Figure 1;
- use logarithmic transforms for both \(R\) and \(c\);
- bind each marker to the diameter legend rather than treating the connected visual curve as a data table;
- preserve raw marker coordinates and the calibration record;
- do not synthesize evenly spaced tick arrays from panel bounds;
- predeclare how a point at exactly 3900 is obtained (an identifiable marker, a bracketed same-diameter interpolation, or no value);
- preserve the distinction between direct force method and digitized figure output;
- independently rerun the extractor from the immutable PDF hash and audit negative mutations (changed tick/marker/common-mode inputs must fail closed).

No value, interpolation result, or uncertainty/tolerance is produced here. The conclusion is only:

~~~text
chart scale: controlled-digitization feasible in principle
target value: not extracted
target status if extracted: DERIVED, never DIRECT-TARGET by figure alone
~~~

### 5.2 Why the chart does not close the target by itself

The plotted family combines several diameters and symbols, and the target decade is visually dense. The paper does not provide a point ledger associating every plotted marker with exact \(R\), \(c\), diameter, velocity and force. A smooth line through markers is not a primary table. Even a reproducible pixel coordinate would remain a figure-derived value and would not supply the missing uncertainty, blockage, surface, turbulence or support-correction fields.

## 6. Gap matrix for "total_mean_drag_coefficient_re3900"

| Gate field | Status from this audit | Permitted use |
|---|---|---|
| public primary artifact | "PASS" | may be cited and hash-pinned |
| direct total-force method | "PASS" at method level | may support a "DIRECT-METHOD" authority row |
| total-drag normalization | "PASS" via \(D=cS\rho V^2/2\) | may map "c" to modern "Cd" with notation note |
| \(R=Vd/\nu\) and broad target range | "PASS" | proves target lies in reported range |
| exact Re=3900 direct point | "MISSING" | cannot fill direct target field |
| figure-derived target possibility | "PENDING CONTROLLED DIGITIZATION" | may produce a "DERIVED" candidate only |
| finite-span/periodic-span equivalence | "MISSING" | source is infinitely-long/long-cylinder context |
| blockage and end corrections | "MISSING" | must remain explicit unknowns |
| surface/turbulence/support metadata | "MISSING" | cannot infer from figure |
| uncertainty / tolerance | "MISSING" | no tolerance may be invented from this source |

## 7. Follow-up preconditions and non-decision

This source is suitable for the root agent to decide whether to pre-freeze a controlled Figure 1 digitization rule. Before any derived value is used in a receipt, the root-owned rule should bind:

1. the official PDF URL and SHA-256 above;
2. the exact PDF page and printed Figure 1;
3. deterministic renderer/version, rotation and crop;
4. actual log-axis tick/grid coordinates and marker inventory;
5. the target-point/interpolation decision made before extracting the point;
6. raw and calibrated outputs plus an independent audit.

Those controls would make the extraction reproducible; they would not repair the missing source uncertainty, blockage, finite-span equivalence or correction data. The present audit therefore hands off the following frozen-status wording:

~~~text
Wieselsberger/NACA TN 84:
  primary artifact: admissible
  direct total-drag method: supported
  total_mean_drag_coefficient_re3900 direct target: not established
  controlled Figure 1 digitization: eligible to be pre-frozen
  candidate value/tolerance/ACCEPT: intentionally not decided here
~~~

## 8. Reproducibility commands and local evidence

The commands used for the independent artifact check were:

~~~text
url='https://ntrs.nasa.gov/api/citations/19930080855/downloads/19930080855.pdf'
curl -L --fail --max-time 90 -sS -D response.headers -o NACA-TN-84.pdf "$url"
sha256sum NACA-TN-84.pdf
pdfinfo NACA-TN-84.pdf
file --mime-type NACA-TN-84.pdf
curl -L --fail --max-time 30 -sS \
  -D ntrs-api.headers -o ntrs-citation.json \
  'https://ntrs.nasa.gov/api/citations/19930080855'
pdftoppm -jpeg -r 180 -f 1 -l 16 NACA-TN-84.pdf render-180dpi/page
sha256sum render-180dpi/*.jpg
~~~

The downloaded PDF, headers, metadata, fulltext response, deterministic page render and render hashes are outside the repository at the external directory recorded in §1.2. No solver, benchmark, long statistic, candidate, receipt or Stage 5 process was run or modified for this audit.
