# Tritton (1959) total-drag primary-source audit

审计日期：2026-08-22（UTC）。

审计范围：只核查 D. J. Tritton, *Experiments on the flow past a circular
cylinder at low Reynolds numbers*, *Journal of Fluid Mechanics* **6**(4),
547--567 (1959), DOI `10.1017/S0022112059000829`，是否能为
`total_mean_drag_coefficient_re3900` 提供公开、可审计的 primary authority。
本文不选取 Re=3900 数值、不做图表数字化、不选择 HUNDUN 容差、不决定
ACCEPT/REJECT，也不修改 solver、plans/spec、references、receipt、ledger、
candidate、COAST 或 Stage 5。

## 结论摘要

| 冻结问题 | 当前可审计结论 |
| --- | --- |
| 是否为一手文献 | 是。Cambridge University Press 的文章记录确认作者、题名、期刊卷期、页码和 DOI。 |
| 是否直接测 drag | 仅在方法摘要层面成立：publisher abstract 明确说用 quartz-fibre bending 测量 circular-cylinder drag。完整 force equation 和校准细节不在公开文章页。 |
| 是否有 Re=3900 direct target | 否。publisher abstract 明确给出的测量范围是 `Re=0.5--100`，没有 Re=3900 的表、图点或数值。 |
| 是否能把值标为 derived | 本审计没有提取或插值任何值；不能用该来源制造 Re=3900 derived target。 |
| 是否提供 total/mean 语义 | 否。公开摘要没有 total-force normalization、时间平均协议或测量窗口定义。 |
| 是否支持 HUNDUN finite-span/periodic-span 等价 | 否。公开文章页没有 active span、端部/支撑、blockage 或 correction 证据。 |

冻结结论：Tritton 1959 可以作为“低 Re 圆柱直接拖曳测量方法存在、且测量范围远低于
Re=3900”的历史 primary evidence；不能填充或封印
`total_mean_drag_coefficient_re3900`，也不能解除 literature receipt 的 force gap。

## 1. Primary artifact 与 provenance

### 1.1 Publisher record

- [Cambridge Core article record](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/experiments-on-the-flow-past-a-circular-cylinder-at-low-reynolds-numbers/0386A4A3A98750248AEB772532863BB6)
  是 Cambridge University Press 的一手文章记录。
- 记录的作者为 D. J. Tritton，affiliation 为 Cavendish Laboratory, Cambridge；
  题名为 *Experiments on the flow past a circular cylinder at low Reynolds numbers*。
- 记录确认 *Journal of Fluid Mechanics*, Volume 6, Issue 4, November 1959,
  pp. 547--567，DOI `10.1017/S0022112059000829`。
- publisher abstract 的关键范围声明是：Part I 测量 circular-cylinder drag，
  方法为观察 quartz fibres 的弯曲，Reynolds-number range 为 `0.5--100`；
  Part II 讨论约 Re=90 的 vortex-street transition。
- 该文章记录列出一条同作者、同年、明确引用的 *Variation of young's modulus
  of fused quartz fibre with diameter*（*Philosophical Magazine* 4(42), 780）。
  该条目是 fibre calibration 背景，不是 Re=3900 direct drag artifact，故不作为
  本目标的 force authority。

### 1.2 本次取得的字节 artifact

为了保留访问状态且不覆盖其他 lane 的内容，artifact 放在新的仓库外目录：

```text
/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/
  tritton-1959-s0022112059000829-20260822/
```

| 文件 | URL/用途 | 获取时间（UTC） | MIME / bytes | SHA-256 |
| --- | --- | --- | ---: | --- |
| `cambridge-landing.html` | [publisher article record](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/S0022112059000829) 的最终文章页 | `2026-08-22T08:15:14Z` | `text/html`, 819346 | `5e45e5c0655c03da25fc0c15728a72ce17a500475f8ea23504084874b4b1a9d2` |
| `cambridge-article-access-response.html` | publisher 的 citation-PDF endpoint 返回的访问响应；文件类型核验为 HTML，不是 PDF | `2026-08-22T08:15:28Z` | `text/html`, 819346 | `1ef25b4ba91b4338fb514fa07a192d4b9398e9ea6f0f6463d621c910341c969a` |

完整正文 PDF 未能从 publisher 公开 endpoint 合法取得：该 endpoint 返回 Cambridge
文章访问页，`file` 核验不是 `application/pdf`。本审计没有把第三方镜像、后续 CFD
汇总、搜索摘要或二手表格当作 numeric authority，也没有保存或使用任何 Re=3900
数值。

## 2. 证据矩阵

证据等级：`P0` 为 publisher article record/abstract；`G` 为缺失或未能由公开
publisher artifact 审计的字段。

| 冻结字段 | Primary evidence | 状态 / permitted use |
| --- | --- | --- |
| artifact provenance | Cambridge Core article record，DOI、作者、期刊、卷期和页码可核验 | `P0 PASS`；可作文献身份与范围证据 |
| direct force method | abstract 明确写 drag measurements are made by observing bending of quartz fibres | `P0 PARTIAL`；支持“直接拖曳测量方法”存在，不足以审计 force equation、tare 或 calibration |
| Re definition / range | abstract 给出 measurement range `Re=0.5--100`；公开页未给逐点 Re definition equation | `P0 RANGE-PASS`；明确排除 Re=3900 direct target |
| Re=3900 direct/derived status | abstract range 不覆盖 3900；公开记录无 Re=3900 row、marker 或 printed value | `MISSING / REJECT`；本审计不派生任何值 |
| total-drag semantics | 摘要只写 drag on circular cylinders，没有 total-force decomposition 或 per-length convention | `MISSING`；不能升级为 HUNDUN total mean force authority |
| `C_D` normalization | 公开文章页未提供 coefficient equation、reference area、dynamic-pressure convention 或 diameter/span denominator | `MISSING`；不得套用现代标准归一化 |
| geometry / orientation | 公开摘要仅说 circular cylinders；未给 axis orientation、diameter、length/span、supports、surface finish 或 boundary conditions | `MISSING`；不能建立 finite-span/periodic-span equivalence |
| active span / end conditions | 未见 active span、端板、支撑、端部修正或无限长近似的可审计 apparatus 数据 | `MISSING`；不能映射 HUNDUN `span=pi D` |
| blockage / wall correction | 未见试验段尺寸、blockage ratio、wall correction 或 wake-interference correction | `MISSING` |
| inlet / fluid state | 公开摘要未给流体状态、来流湍流、Mach/compressibility 或 velocity calibration | `MISSING` |
| uncertainty / repeatability | 公开摘要未给误差、置信区间、重复性、采样时间或 mean-window protocol | `MISSING`；不能从它推导 HUNDUN tolerance |
| permitted use | 方法历史、范围否定、force-authority gap record | 允许作为 negative/qualifying literature evidence；不允许填 receipt 数值字段 |

## 3. Re=3900 与“直接 total mean drag”的拆分

### 3.1 方法级证据不等于 target-level authority

quartz-fibre bending 是力学读数，因而 publisher abstract 至少支持
`DIRECT-DRAG-METHOD` 的方法级标签。但摘要没有说明：

1. 读数对应有限圆柱的总力还是单位长度力；
2. 参考面积和 `C_D` 分母如何定义；
3. 是静态/时间平均读数，还是某个明确采样窗口的 mean；
4. 支撑、端部、阻塞、壁面和校准如何处理。

因此即使方法名称是 direct drag，也不能写成已经满足
`DIRECT-TOTAL-MEAN-FORCE`。

### 3.2 目标 Re 明确不在该 primary 的报告范围

publisher abstract 给出的 Re 范围上限是 100。Re=3900 超出该范围约两个数量级；
公开 article record 没有 Re=3900 的表项、图点、原始力或系数。因此：

```text
Tritton 1959 direct Re=3900 value: ABSENT
Tritton 1959 figure-derived Re=3900 value: NOT EXTRACTED
Tritton 1959 permitted role: historical method/range evidence only
```

不能把“同一圆柱问题”或后续文献对 Tritton 的引用当作 Re=3900 authority，也不能
用该 paper 的低 Re curve 外推到目标。

## 4. 同系列引用的边界

Cambridge 记录列出 Tritton 同年关于 fused-quartz-fibre Young's modulus 的
*Philosophical Magazine* 条目。它与 quartz-fibre calibration 有关，但不是圆柱
drag-vs-Re 数据集；即使取得全文，也只能支持仪器材料参数，不能补足 Re=3900
total mean drag、span、blockage 或 uncertainty。本审计不把它加入 force receipt，
不以它替代 Tritton 1959 正文，也不从中做任何数字化。

## 5. Gate hand-off

```text
Tritton 1959:
  publisher primary record: admissible
  direct quartz-fibre drag method: supported at abstract level only
  Re=3900 direct target: absent (publisher range Re=0.5--100)
  total/mean normalization: not auditable from public artifact
  finite-span / blockage / correction / uncertainty: not auditable
  controlled digitization: not started and not authorized by this audit
  receipt use: negative/qualifying evidence only
```

该结论不声称全文中不存在超出摘要的任何图表；它只基于公开 publisher artifact
能直接审计的范围，并在正文未合法取得时保持 fail-closed。若后续取得合法、完整的
publisher 或作者/机构版本，必须逐页复核其 apparatus、equation、table/figure 和
uncertainty；在此之前不能改变上述状态。

## 6. Reproducibility commands

以下命令只写入仓库外的上述 audit directory；不构建、不测试、不运行 MPI：

```bash
AUDIT_DIR=/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/literature/tritton-1959-s0022112059000829-20260822
mkdir -p "$AUDIT_DIR"

curl -L --fail --retry 2 --connect-timeout 20 --max-time 120 \
  -A 'Mozilla/5.0' \
  -o "$AUDIT_DIR/cambridge-landing.html" \
  'https://www.cambridge.org/core/services/aop-cambridge-core/content/view/S0022112059000829'

curl -L --fail --retry 2 --connect-timeout 20 --max-time 180 \
  -A 'Mozilla/5.0' \
  -o "$AUDIT_DIR/cambridge-article-access-response.html" \
  'https://www.cambridge.org/core/services/aop-cambridge-core/content/view/0386A4A3A98750248AEB772532863BB6/S0022112059000829a.pdf/div-class-title-experiments-on-the-flow-past-a-circular-cylinder-at-low-reynolds-numbers-div.pdf'

file "$AUDIT_DIR"/*.html
stat --format='%n|bytes=%s|mtime=%y' "$AUDIT_DIR"/*.html
sha256sum "$AUDIT_DIR"/*.html
```

本 audit 仅新增本文档；外部 artifact 的内容和 hash 如 §1.2 所列。完整文件 diff、
`git diff --check` 和最终 review 由 parent agent 执行；本 worker 不提交、不修改其他
文件。
