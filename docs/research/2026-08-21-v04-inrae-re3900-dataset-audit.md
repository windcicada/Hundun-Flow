<!-- SPDX-License-Identifier: Apache-2.0 -->

# 2026 INRAE Re=3900 PIV 数据集审计

日期：2026-08-21

范围：只审计 Recherche Data Gouv/Dataverse DOI `10.57745/DHJXM6` 的官方
V1 metadata、文件清单和抽取到 `/tmp` 的数据文件；Parnaudeau 等 (2008) 原论文
仅用于实验身份、坐标和 Fig. 11--15 归一化的对照。不读取 HUNDUN 输出，不把本数据
集改名为 Parnaudeau 2008 数组。

## 判定摘要

这份数据集是公开的一手数据源，文件许可和下载权限清楚，因而就“可合法取得并作为
独立 2026 实验补充”而言是 **ADMISSIBLE**。但是，官方 metadata 没有给出圆柱直径
`D`、参考速度 `Uc`、坐标原点/方向、实验装置型号或视场与 2008 批次的对应关系；V1
文件清单也没有描述中宣称的均值、Reynolds stress、偏度和峰度文件。当前证据不足以
把它重构成或替代 Parnaudeau Fig. 11--15 的三站原数组，结论为 **INSUFFICIENT**。

在冻结上述缺口后，可以把由 2026 瞬时场重建的三个站剖面作为另一个实验数据源，并
明确标注“2026 INRAE 数据集”，不能标成 2008 PIV 曲线；这种用途是
**SUPPLEMENT_ONLY**，不完成当前 2008 Fig. 11--15 literature authority。

逐项结论：

| 审计项 | 结论 | 证据和限制 |
|---|---|---|
| DOI、版本、日期、发布方 | `ADMISSIBLE` | 官方 API 返回 dataset `720391`、V1.0、`RELEASED`，发布时间 `2026-02-19T18:24:44Z`。 |
| 公开许可和文件权限 | `ADMISSIBLE` | API 的 `license.name=etalab 2.0`，`fileAccessRequest=false`，998 个文件均 `restricted=false`。 |
| 2026 数据确为 2D2C PIV 瞬时场 | `ADMISSIBLE` | 三个完整下载文件的 DaVis header 都是四列 `x [mm]`,`y [mm]`,`Vx [m/s]`,`Vy [m/s]`。 |
| 2026 与 2008 是否同一实验批次/装置 | `INSUFFICIENT` | 2026 只说“INRAE wind tunnel”；2008 说 Rennes Regional Center of Cemagref，并给出 D=12 mm、端板和仪器细节。官方 2026 metadata 没有 venue、装置型号、圆柱/端板或批次字段；虽有作者 Dominique Heitz 重合，也不能据此证明相同。 |
| `D`、`Uc`、坐标原点及正方向 | `INSUFFICIENT` | 文件只有带单位的绝对 mm 坐标，没有 `D`、`Uc`、`x/D`、`y/D`、圆柱中心、来流方向或 `z`/平面声明。 |
| 采样平面、有效视场、样本数 | `INSUFFICIENT` | 描述声称 10,000 instantaneous fields，但 V1 清单只有 998 个文件且缺 `Serie_10902.txt`；无 plane/FOV/采样时间或缺失文件说明。 |
| 均值、Reynolds stress 的归一化 | `INSUFFICIENT` | V1 实际清单没有统计文件；瞬时文件的速度单位是 m/s，不能从文件本身证明 `Uc` 或 `Uc^2` 归一化。 |
| 三站 `x/D=1.06,1.54,2.02` 覆盖 | `INSUFFICIENT` | 原始坐标在物理范围上可能覆盖条件站点，但缺少 `x0,y0,D` 映射，不能合法声明三站覆盖。 |
| 作为独立 2026 补充数据 | `SUPPLEMENT_ONLY` | 可在实验映射、缺失文件、有效向量掩码、统计公式和插值误差规则冻结后重建；不得替换 2008 数组。 |

## 一手来源和定位

### 2026 INRAE/Recherche Data Gouv

* [数据集落地页](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/DHJXM6)
  是官方展示页；其 JSON-LD distribution 同样列出 998 个 `text/plain` 文件。
* [Dataverse V1 metadata API](https://entrepot.recherche.data.gouv.fr/api/datasets/:persistentId/?persistentId=doi:10.57745/DHJXM6)
  是本审计的权威清单来源。读取 `data.latestVersion`：
  `id=283288`、`datasetId=720391`、`versionNumber=1`、
  `versionMinorNumber=0`、`versionState=RELEASED`、
  `createTime=2026-02-19T18:03:37Z`、`releaseTime=2026-02-19T18:24:44Z`、
  `lastUpdateTime=2026-02-19T18:24:44Z`、`publicationDate=2026-02-19`。
* metadata 的标题是 *Non-time-resolved PIV dataset of flow over a circular
  cylinder at Reynolds number 3900*，作者为 `GEORGEAULT, Philippe` 和
  `HEITZ, Dominique`，联系人为 Dominique Heitz/INRAE。`dsDescription` 声称：
  INRAE 风洞中的 planar 2D2C PIV、10,000 instantaneous velocity fields、mean
  velocity/Reynolds stress、Skewness/Flatness，并称其补充和改进 Parnaudeau et al.
  (2008) 数据库。这个描述是作者的主张，不等于文件清单已实现这些类别。
* metadata 的许可对象为
  [`etalab 2.0`](https://spdx.org/licenses/etalab-2.0.html)，文件访问请求为
  `false`。这是合法公开使用判断的直接依据；不推断实验身份或统计含义。

### 2008 对照论文

* [Parnaudeau, Carlier, Heitz & Lamballais (2008), DOI
  10.1063/1.2957018](https://doi.org/10.1063/1.2957018) 是文章的 DOI 身份。
* [作者/机构 HAL 记录 hal-00383669v1](https://hal.science/hal-00383669v1) 确认
  题名、作者、期刊、2008 年和 DOI；HAL 记录当前标为“Fichier non déposé”。
* [作者上传的全文副本](https://www.researchgate.net/publication/234130406_Experimental_and_numerical_studies_of_the_flow_over_a_circular_cylinder_at_Reynolds_number_3900)
  用于按页核对原文的实验段落和 Fig. 11--15；其内容是该论文全文而非 CFD 二手汇总。
  关键位置为原文 pp. 085101-2--3（实验设置）、p. 085101-11（Fig. 11--15）和
  p. 085101-14 footnote 31。

## 2026 文件清单和数据文件审计

### 完整 V1 清单

API 的 `data.latestVersion.files` 有 998 条记录：

* labels 的整数范围为 `10001..10999`，唯一缺口是 `10902`，即缺
  `Serie_10902.txt`；不是连续的 10,000 个 snapshot 文件。
* 998 个文件全为公开 `text/plain`，唯一 content type 为 `text/plain`，总字节数
  `13,176,725,066`（约 13.18 GB）。每条记录含 data-file id、文件名、每文件 DOI
  persistent ID、`filesize`、官方 `MD5` checksum 和访问 URL。
* 没有文件夹、tag、description 或统计类别字段；清单中没有独立的 mean、stress、
  skewness 或 flatness 文件。不能以 dataset description 代替文件证据。

V1 raw API response 保存到 `/tmp/dhjxm6-metadata.json` 的字节数为 `499,910`，
SHA-256 为
`02f0a612e515788073e3380be99f21df73e18a01634fbd81afaecaa9c2f55aed`。
按文件名排序、逐行拼接如下八列：
`label|dataFile.id|persistentId|pidURL|filesize|checksum.type|checksum.value|downloadURL`，
得到 187,624 bytes 的完整可重建清单，SHA-256 为
`1ccd37fc99e36c944c2f737e6c80648a7956e7937a50d50f61e7957df058ea26`。
这两个 digest 绑定的是官方 API 响应和清单，而不是把 13 GB 数据复制进仓库。

三条代表性记录（也实际完整下载核对）如下：

| 文件 | 官方 data-file URL | persistent ID | bytes | 官方 MD5 |
|---|---|---|---:|---|
| `Serie_10001.txt` | [`/api/access/datafile/720393`](https://entrepot.recherche.data.gouv.fr/api/access/datafile/720393) | [`doi:10.57745/LPNPL5`](https://doi.org/10.57745/LPNPL5) | 13,198,574 | `ad8e65c5f4e028763db96a957629d644` |
| `Serie_10500.txt` | [`/api/access/datafile/720862`](https://entrepot.recherche.data.gouv.fr/api/access/datafile/720862) | [`doi:10.57745/TE2GZ5`](https://doi.org/10.57745/TE2GZ5) | 13,221,466 | `ca6fba7773f3555aafcdd0757df992f8` |
| `Serie_10999.txt` | [`/api/access/datafile/721086`](https://entrepot.recherche.data.gouv.fr/api/access/datafile/721086) | [`doi:10.57745/0KEKUM`](https://doi.org/10.57745/0KEKUM) | 13,216,581 | `34cfb2ff14efdf2cd049a2ec4fc766c5` |

三次下载的本地 MD5 分别与上述官方 checksum 完全相同。其余 995 条记录的
URL/persistent ID/bytes/MD5 可由下面的只写 `/tmp` 命令逐字重建，不需要猜测或手工
录入：

```bash
API='https://entrepot.recherche.data.gouv.fr/api/datasets/:persistentId/?persistentId=doi:10.57745/DHJXM6'
curl -fsSL --retry 2 "$API" -o /tmp/dhjxm6-metadata.json
python3 - <<'PY'
import json
v=json.load(open('/tmp/dhjxm6-metadata.json'))['data']['latestVersion']
rows=[]
for f in sorted(v['files'], key=lambda x:x['label']):
    d=f['dataFile']
    rows.append('|'.join((f['label'], str(d['id']), d['persistentId'], d['pidURL'],
                          str(d['filesize']), d['checksum']['type'], d['checksum']['value'],
                          f'https://entrepot.recherche.data.gouv.fr/api/access/datafile/{d["id"]}')))
open('/tmp/dhjxm6-inventory.psv','w').write('\n'.join(rows)+'\n')
PY
sha256sum /tmp/dhjxm6-metadata.json /tmp/dhjxm6-inventory.psv
```

### 实际数据文件格式

命令 `curl -fL URL -o /tmp/Serie_*.txt` 完整取得上表三文件后，三者均有相同首行：

```text
#DaVis 10.x 2C vector field 4 545 740 "x [mm]";"y [mm]";"Vx [m/s]";"Vy [m/s]"
```

`Serie_10001.txt` 的核对结果为 403,300 个数据行（加 header 共 403,301 行）、
740 个 x 点 × 545 个 y 点；坐标按文件顺序为 x 增、y 减，步长均约 `0.1833 mm`：

* x：`-60.1328 .. 75.3276 mm`，最近零点 `-0.0096616 mm`；
* y：`60.1803 .. -39.5362 mm`，最近零点 `0.0571616 mm`；
* 四列均可按分号解析为有限浮点数；三文件没有额外列、均值、应力或站点标签。

文件边缘可见成片的 `0;-0` 数值，但 API/文件 header 没有声明这是圆柱遮罩、无效向量
还是有效的零速度；统计重建不能把它们自动当作有效样本。也没有官方缺失值掩码或
异常向量替换规则。原始物理窗口的坐标差约 `135.4604 mm × 99.7165 mm`，但有效
视场必须由后续实验映射和掩码审计确定，不能仅以矩形边界作为 FOV 结论。

## 2008 实验身份和 Fig. 11--15 的可比量

按论文实验设置（pp. 085101-2--3）：

* 装置是 Rennes Regional Center of Cemagref 的风洞；测试段横截面宽 `H=28 cm`、
  长 `100 cm`，自由来流湍流强度小于 `0.2%`。
* 圆柱长 `L=280 mm`、直径 `D=12 mm`，带薄矩形端板；端板间距 `240 mm`、
  `L/D=20`，blockage `4.3%`，圆柱在测试段入口后 `3.5D`。
* 轴定义是风洞纵向轴 x、圆柱轴 z，法向轴按右手系，原点在圆柱中心。PIV 主测量
  为 `z=0` 平面、圆柱后方的 2D2C；两台相机视场分别 `3.6D×2.9D` 和
  `1.9D×1.6D`，统计剖面主要来自小视场 case 2。
* `Uc=4.8 m/s` 对应 `Re=3900`；每台相机 5,000 image pairs，帧间隔 `25 μs`。
  论文还说明了相关窗口、50% overlap、错误向量 median 检查/局部均值替换，误差
  小于 `0.3 pixel`、错误向量少于 `0.1%`。这些是 2008 的处理规则，不是 2026
  文件的 metadata。

Fig. 11--15 的三站明确为 `x/D=1.06, 1.54, 2.02`（论文 pp. 085101-11--12）：

| 图 | 量 | 原文归一化 |
|---|---|---|
| Fig. 11 | mean streamwise velocity | `⟨u⟩/Uc`，横坐标 `y/D` |
| Fig. 12 | mean normal velocity | `⟨v⟩/Uc`，横坐标 `y/D` |
| Fig. 13 | streamwise fluctuation variance | `⟨u'u'⟩/Uc²`，横坐标 `y/D` |
| Fig. 14 | normal fluctuation variance | `⟨v'v'⟩/Uc²`，横坐标 `y/D` |
| Fig. 15 | fluctuation covariance | `⟨u'v'⟩/Uc²`，横坐标 `y/D` |

论文 footnote 31 明确说实验和数值统计可通过联系作者获得；文章本身没有公开的
数值数组附件。因此，2026 数据的作者描述“completes and improves”不能反向证明
它就是 2008 Fig. 11--15 的原始数组。

## 为什么当前不能认定同批次或三站已覆盖

2026 文件的 740×545、约 135×100 mm 原始网格与 2008 论文中最终 `160×128` PIV
处理网格和两种小视场描述不同；这可能是不同采集/导出层级，也可能是重新导出的
原始场，官方材料没有说明。更关键的是：

1. 2026 metadata 没有 `D=12 mm`、`Uc=4.8 m/s`、空气黏度/温度、圆柱中心、端板、
   `z=0` 或流向正负；只能看到物理单位 mm 和 m/s。
2. 若**假设** `x0=y0=0` 且沿用 2008 `D=12 mm`，三站会对应 `x=12.72,18.48,
   24.24 mm`；若 D 或原点不同，站点完全不同。这个条件换算不是数据集证据，
   不能进入冻结 receipt。
3. 998 个实际文件并没有统计类别；描述中的“10,000 + mean/stress/skewness/
   flatness”与公开 V1 清单冲突。缺 `Serie_10902.txt` 的影响也不能靠重复某个
   snapshot 或把三文件样本数当成 10,000 来补齐。
4. 瞬时文件的 `Vx/Vy` 是 m/s，未定义来流参考值、零点、有效向量 mask、时间权重、
   样本是否独立或是否来自同一相机/批次，所以不能合法宣称与 2008 的
   `⟨u⟩/Uc`、`⟨u'u'⟩/Uc²` 等量已经一致。

## 若以后接纳为补充数据，machine receipt 必须冻结的内容

以下是可执行的证据要求，不是本审计替主 agent 决定坐标或插值误差规则：

1. 保存上述 API raw JSON SHA-256 和完整 998-row inventory digest；由数据发布方解释
   `10,000` 描述、998 文件、缺失 `Serie_10902.txt` 以及统计文件缺失之间的关系。
2. 由官方实验资料或作者确认 `x0,y0`、x 正方向、`D`、`Uc`、运动黏度/温度、圆柱
   端板/阻塞和 `z=0` 平面，并明确 2026 批次是否为 2008 同一装置/几何。没有这组
   绑定，三站结论保持 `INSUFFICIENT`。
3. 明确有效向量与 `0;-0` 的 mask、异常向量处理、实际样本数及时间/相机权重。重建
   统计的候选定义应至少记录：
   `ū=mean(Vx)`, `v̄=mean(Vy)`,
   `u'u'=mean((Vx-ū)^2)`, `v'v'=mean((Vy-v̄)^2)`,
   `u'v'=mean((Vx-ū)(Vy-v̄))`，并由主 agent 冻结总体 `N` 还是无偏 `N-1` 分母。
4. 站点映射应先用冻结的 `(x0,y0,D)` 计算 `x/D`、`y/D`，再由主 agent 选择一维剖面
   插值/网格积分方法，记录原始相邻点、插值残差和空间/采样误差；不能把未证明的
   坐标零点或统一 `±10%` 误差写入 receipt。
5. 最终字段必须带 provenance：`dataset DOI/version`、每个输入文件的 PID/MD5/bytes、
   统计公式、mask、站点坐标、插值方法和误差；结果名称保持“2026 INRAE PIV”，
   不得覆盖或重命名 Parnaudeau 2008 曲线。

## 可复核命令和日志摘要

本审计使用的只读命令为：

```bash
curl -fsSL --retry 2 \
  'https://entrepot.recherche.data.gouv.fr/api/datasets/:persistentId/?persistentId=doi:10.57745/DHJXM6' \
  -o /tmp/dhjxm6-metadata.json
curl -fL --retry 2 https://entrepot.recherche.data.gouv.fr/api/access/datafile/720393 \
  -o /tmp/dhjxm6-Serie_10001.txt
curl -fL --retry 2 https://entrepot.recherche.data.gouv.fr/api/access/datafile/720862 \
  -o /tmp/dhjxm6-Serie_10500.txt
curl -fL --retry 2 https://entrepot.recherche.data.gouv.fr/api/access/datafile/721086 \
  -o /tmp/dhjxm6-Serie_10999.txt
md5sum /tmp/dhjxm6-Serie_10001.txt /tmp/dhjxm6-Serie_10500.txt /tmp/dhjxm6-Serie_10999.txt
```

日志摘要：三次 MD5 分别为 `ad8e65c5f4e028763db96a957629d644`、
`ca6fba7773f3555aafcdd0757df992f8`、`34cfb2ff14efdf2cd049a2ec4fc766c5`，与官方
API 完全一致；三个 header/网格结构一致，解析没有坏行。没有下载其余约 13.14 GB，
因为当前数据已经判定不足以替代 2008 数组，官方 per-file MD5 已由 metadata 绑定。
缓存只在 `/tmp`，没有把数据文件写入仓库。

未解决点保持为：实验批次/装置同一性、`x0,y0,D,Uc`、平面/有效 FOV、缺失文件和
统计文件解释、有效向量 mask、统计分母、三站插值和误差规则。它们解决前，本审计
不授权 literature receipt complete，也不授权长统计或最终 gate。
