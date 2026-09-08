# 产品耦合接线进度

用户要求以实际产品接线验收为合并完成条件，然后才能续算圆柱绕流。
`f4f04ec` / `de0e77a` 仅归集源码与验证单相兼容，不代表完整耦合验收。
圆柱服务 `hundun-re3900-integrated-long-20260908.service` 保持 SIGSTOP；
暂停时最后完整健康记录为 10577，续算源仍为独立封存的 10510 checkpoint。

## 已接通的产品路径

| 路径 | 产品行为 |
| --- | --- |
| 显式源项准入 | `ContributionAdmission` 绑定模型/气相身份；编译计划和实际方程消费同时校验 capability 与 source identity。默认无准入，旧惰性计划指纹保留。 |
| finite-rate mean | 通过完整 `(p,h,Y)` 查询得到净组分质量生成率，按名称映射到 Ns−1 产品字段；核查机理、组分、元素、分子量、焓参考、修订号、输出缓冲和原生 EOS/热容一致性。源项进入原生 BDF/EX2 非对流速率历史。 |
| PaSR | 使用同一已核查净反应率；`reactant_depletion_l1_v1`、单元体积立方根滤宽、气相质量分数加权分子扩散率和实际 `mu_eff-mu` 接入既有混合时间及反应分数函数。显式配置 C_Z 与 Sc_t。 |
| 压力/能量 | 组分改变后由现有两次 PISO / `rho*h-p` 联立闭合响应；总热化学焓不重复加入形成焓放热。 |
| 初始与失败状态 | 均匀反应初值也构造首步化学速率，修复只在 IBM 初值重建速率导致首步漏反应的问题。单 rank provider 失败或输出合同损坏使所有 rank 拒绝；已接受物理字段、两层历史和通量不变。 |
| 输入与恢复 | 严格 reaction JSON、跨 rank wire、case/模型指纹与方法历史签名；现有 Restart V3 保存反应速率历史。载入值精确一致，带放热续算按求解器无量纲终止容差比较；惰性同分异构试验连续/恢复值精确相同。 |
| 实际后端 | `HUNDUN_ENABLE_REACTING_CANTERA=ON` 将已验证 Cantera 3.2.0 包接入产品；关闭时 direct_cantera 明确拒绝。analytic_isomer 是具名合成验收表示，不是燃料机理。 |

当前 finite-rate / PaSR 使用原生显式 BDF/EX2 反应率历史。
Cantera 区间积分器在后续 ESF 接线中消费；此节点不宣称已接通 ESF 两个化学半区间。

## 当前验证

- Clang/libc++ focused 25/25 PASS：产品冻结、原生方程、压力能量重试、历史恢复、
  反应准入与 finite-rate/PaSR 的 1/2/4-rank 产品测试。
- GCC11/Cantera 产品矩阵 12/12 PASS：外部 provider、具名解析表示、PaSR、真实 Cantera，
  各 1/2/4 ranks；包含带放热的密度与 `rho*h-p` 守恒、实际磁盘 Restart 及续算比较。
- 两个真实 CLI 两步计算完成：Clang 解析表示与 GCC11 Cantera；两组 runtime evidence
  validator 均返回 0。Cantera 机理为原创 A→B 合成夹具，其 SHA 随 case.json 冻结。
- GCC 编译确认并修复了既有代码依赖间接包含 `<algorithm>` / `<climits>` 的问题。
- 日志和运行目录在自有工作树的 `build-coupling-evidence/`、`build-integrated-cantera/`；
  未覆盖封存的圆柱数据，也未运行新的圆柱长算。

## 正在接入的质量交换入口

- 原生连续性与压力装配增加了目标时间步质量源，单位 kg/(m³ s)，冷计划绑定模型身份，
  C1/C2 和终止审计使用同一不可变修订。关闭质量交换时沿用空入口和既有指纹。
- 终止审计的 V2 单元公式显式减去质量源；V1 ABI 仍按零源调用。
  压力线性收敛检查和失败见证采用同一源项与归一化尺度。
- 非均匀网格、全局单元编号构造的空间变化源通过 1/2/4 ranks 连续性测试；
  错误模型身份与非有限末单元源在输出发布前被拒绝。
- 当前 focused 20/20 PASS，覆盖压力权限/数值变更、连续性、原生产品冻结、
  压力能量重试和 finite-rate/PaSR 产品测试。
- 预测器另设 `PredictorRateHistory.current`，使当前步质量、组分、焓交换只进入一次，
  不经 EX2 历史外推。预测密度与配对质量通量审计均含同一质量源。
  当前步耦合源导致状态不可接受时，整个尝试失败，不进入只缩放气相源的低阶端点。
  当前 focused 扩展到 22/22 PASS，新增 BDF2 精确源项代数、过量消耗回退、
  非有限源与过期时间身份拒绝，并覆盖原有 IDP 路径。
- 这些是实际方程的质量源消费入口；产品喷雾尚未提供该源，尚未形成完整两相接线。

## ESF 持续场与过滤化学源接线（周期域增量验收）

- ESF/TCR 配置进入严格 JSON、模型指纹和跨 rank wire。字段数仅 2/4，seed 用 uint64
  保存；测试使用超过 2^53 的 seed，确认没有浮点化。TCR reactants/progress weights/
  初始分支设置可序列化；这不等同于启用 validated 科学反馈。
- ESF 化学报告补充两次半区间的等权组分质量密度增量，单位 kg/m³，可直接向 MeanState
  输出过滤化学源；失败不发布这个借用结果。新增元素/总质量/形成焓收支检查通过，
  ESF 与 P8 相关测试 6/6 PASS。
- `reacting-esf` 已从公共 ProductCompiler / ProductDriver 进入原生状态层：四个全组分加
  总焓随机场与共同输运缓存同时提交和回退，accepted/previous 写入 Restart V3；Philox
  地址由 seed 与 accepted step 重建，没有保存浮点 RNG 游标。
- 当前步使用原生已提交质量通量、配置的对流重构及中心扩散/梯度算子，经过 IEM 后执行
  两个连续化学半区间。过滤组分密度增量进入预测器 current 源，不写入旧 EX2 化学速率历史。
  原生压力能量闭合后统一平移场均值到 MeanState，越界时整步拒绝，不逐场裁剪。
- 区间积分当前明确要求 backward_euler；BDF2、非周期物理边界、IBM 和 TCR 非 off
  仍明确拒绝。未将这个增量作为完整 Stage 5/6 合并验收。
- Clang focused 26/26 PASS；GCC11/Cantera 产品矩阵 15/15 PASS。ESF 1/2/4 ranks
  对独立解析两半区间的密度加权源、等权均值、方差衰减作检查；真实磁盘重启后继续计算。
  rank 0 第三次 chemistry 调用失败后所有已接受场/历史/两层通量精确保持，清除故障后
  与不中断路径等价；非有限重启随机场在安装任何状态前被拒绝。
- 共同扩散后续增量：MeanState 的总焓/组分扩散与随机场统一到 Gamma=lambda_m/cp+mu_t/Sc_t。
  明确冻结 unity-Lewis 总焓方程，原温度导热方程保持默认；压力焓线性算子使用同一
  Gamma*grad(delta h)，非均匀 cp 的独立有限差分与热源历史/能量残差测试通过。
  输运系数由实际组成、温度重新求本征分子项，不从已加过湍流项的缓存重复累加。
- 周期域正弦扰动经过实际公共 Restart/Driver，跨 rank 的 Halo、配置对流算子、
  共同扩散、IEM、两段化学和均值重整后的方差符合独立离散 Fourier 预期，1/2/4 ranks
  均通过。Clang 扩展 focused 27/27 PASS；GCC/Cantera 产品矩阵 15/15 PASS。
  尚需湍流随机项、物理/IBM 边界与 TCR 持续历史接线。
- 接线参考 HUNDUN governance 只读提交 `8ffdf2b` 的
  `docs/numerics/stage5-esf-tcr-equations.md`、`stage5-coast-esf-semantic-audit.md`，
  以及 ESF transport/element consistency 的公开算法实现。保留当前供体 P8 的
  source→transport/IEM/TCR→两次连续 chemistry 半区间顺序，不移植旧 driver/历史。

## 模型整数历史的原生重启入口

- Restart V4 支持模型定义的定宽单元字节记录，模型身份与宽度纳入完整性和读取合同。
  记录随单元重分区搬运，不借用 FP64 保存 uint64 分支/折点计数。
- 原有分块完整性、读取峰值内存预算和失败时输出保持合同覆盖新增记录；V1–V3 保持原格式。
- 公共 RestartWriter/Reader 四 rank 测试通过：V4 同分区、1→4、4→1 恢复超过 2^53 的
  计数逐字节一致；单 rank 身份/宽度不符和不足内存统一拒绝，已有图像不变。
  格式化后 `v04_io_restart_mpi_4` 1/1 PASS（28.40 s）。
- 此入口尚不等于 TCR 产品历史已经接线；变长 parcel 记录也不在本增量中。

## TCR 周期域产品历史与混合反馈

- TCR 现在从已接受 MeanState 与全部随机场的实际 PH 查询构造 eta/R，使用显式初始
  分支及已接受分支延续；experimental 的 kappa 实际进入 IEM，shadow 只记录观察历史。
  该显式统计时间位置写入方法签名。详细合同见 [产品 TCR 历史](../numerics/product-tcr-history.md)。
- 原生 Restart V4 保存每单元 120 字节的 TCR 历史；分支/折点/修订整数保留整数形式。
  恢复校验补上固定统计映射身份及初始分支与配置的匹配，再统一安装气相和模型历史。
- 历史与编码预分配两套存储。只有原生气相时间步接受后一起交换；在第三个化学半区间
  或全部化学已完成后的最后重整查询中注入 rank 0 失败，所有已接受字段/通量/模型记录
  精确不变。重试后的 TCR 浮点量按原压力求解容差比较，整数仍逐字节比较。
- A→B 独立解析测试 eta=1/4、R=1、signed root=-1/2、kappa=1/3，验证实际方差衰减与
  密度加权气相源。weak denominator 下 shadow 保留历史继续，experimental 整步拒绝。
  产品恢复后的 input revision `9007199254740993` 正确加一，没有 FP64 整数损失。
- 最终 Clang focused 26/26 PASS，包含原生产品/压力能量重试、case、Restart 与反应模式；
  GCC11/Cantera 产品矩阵 21/21 PASS。TCR 和 ESF 都覆盖 1/2/4 ranks。
- 本增量仍为周期域 BE、无 IBM；显式折点证据及折点后续接尚未接入产品。
  validated 不因产品测试通过而获得科学授权。完整 Stage 5/6 合并与圆柱恢复仍未执行。

## 真实 Cantera ESF 与 CLI 重启验收

- 新增真实 Cantera ESF 原创合成机理 fixture，经真实 PH 两段化学、共同输运、空间扰动与
  磁盘恢复验证。原 `rtol=1e-10` 的化学积分使气相组分解析误差为 1.27e-10；该解析验证
  fixture 显式改用 rtol=1e-13、atol=1e-17，原 2e-12 断言保持不变。1/2/4 ranks 通过。
- CLI 首次实测发现 Restart 冻结服务容量遗漏 ESF accepted/previous 元组、输运缓存和
  TCR 字节记录；现已按最大 patch 包络计入，保留每项乘加溢出检查。未绕开写出容量限制。
- 第二个 CLI 阻断为运行证据仍只接受 Restart V1–V3。C++ 证据生产者和独立 Python
  validator 现在接受 V4，验证非零模型记录身份/宽度和原方法签名/manifest 完整性。
- 三种真实 CLI：ESF、TCR experimental、Cantera ESF，均新算两步、写 checkpoint、
  以 1→4 ranks 续算两步，再将证据与起算 manifest 独立核对。
- 最终 Clang focused 29/29 PASS；GCC11/Cantera focused 28/28 PASS，包括原模式回归、
  新真实 Cantera ESF 3 项、CLI 3 项和 evidence workflow 自检。
- 仍未作为完整合并发布。喷雾 parcel 的变长持久状态和原生质量/动量/焓耦合正在接线。

## 构建身份排除生成目录

- 双工具链并行复核暴露：GLOB_RECURSE 顶层 CMakeLists.txt 模式实际递归匹配整个仓库，
  包含另一套构建目录的 CMakeScratch。摘要读取期间文件消失会使配置失败，也会让生成
  文件污染源码内容身份。
- 输入收集现在只递归明确的 v0.4 src/include、cmake 和 yyjson 目录，顶层/版本
  CMakeLists.txt 显式追加。新建/删除嵌套构建临时文件不再改变摘要，真实产品源变化仍改变摘要。
- 目标身份测试在宿主 Clang 和 GCC rootfs 通过；两套真实产品重新配置/构建均成功。
  GCC 测试 fixture 初次因 PATH 缺少已有 /tmp/make 失败，加入该已有工具后复核通过。

## 喷雾变长持久状态的原生入口

- Restart V5 为每个单元提供 u32 长度表与拼接模型字节，不固定每个单元的 parcel 数量。
  V1–V4 保持原格式；V4 TCR 专用恢复拒绝意外的变长记录。
- 读入先统计目标分区长度、确认精确载荷预算，再恢复字段/字节；保留旧图像、模型长度/
  临时偏移、最大源块等均计入峰值上限。详情见 [V5 合同](../numerics/restart-variable-cell-records.md)。
- 1→4、4→1、4→4 分区恢复，混合/全空集合、单 rank 长度和身份不符、整数及预算边界通过。
  Clang 最终 focused 13/13 PASS（40.38 s），GCC Restart+原 ESF/TCR/Cantera CLI 4/4 PASS；
  另一个 GCC 构建身份 fixture 由工具 PATH 修正后单独 1/1 PASS。
- 此增量只提供原生持久化入口，尚未将 parcel/TAB/注射器与气相源消费及共同事务发布。
  V5 运行证据与产品安装会随实际喷雾配置接入；当前 CLI 验证仍是 ESF V3/TCR V4。

## 当前步交换与被动标量范围

- 独立测试复现：原范围构造遗漏 current 源和新增气相质量的稀释，会拒绝合法的有符号
  被动标量。耦合端点现在按 `(rho*q+dt*S_q)/(rho+dt*S_m)` 构造；关闭交换时保持原运算顺序。
- 零或负端点气相库存统一拒绝，已经发布的范围不变；补充输入字段/数组与时间身份检查。
- Clang focused 7/7 PASS，覆盖当前步源的代数、IDP、压力能量重试及 ESF/TCR CLI 续算。
  日志 `build-coupling-evidence/passive-current-{red,green,final}.log`。

## 喷雾冻结配置入口

- `SpraySpec` 保存液体根目录资产/FNV-1a64、整数 seed、parcel/segment 容量、子步控制、
  TAB 开关与点/锥注射器物理参数。液滴单体质量后续由资产密度和配置直径确定。
- 严格 JSON 与 wire 的独立 64 标志位支持该可选节，旧无喷雾 wire 字节不变。资产在根目录
  以有界文件读取并核对精确字节身份；模型指纹包含资产和全部喷雾参数。
- 独立 JSON/wire 测试覆盖 >2^53 整数、未知字段、重复 ID、非法时间控制/方向、路径越界和
  资产字节改变。Clang focused 24/24 PASS，包括各反应模式 MPI 和 ESF/TCR CLI 回归。
- 当前 ProductCompiler 仍明确拒绝 spray 配置，待生命周期和所有气相消费者共同接通后解除。
  此配置增量不代表已运行两相产品，更不代表完整合并完成。

## 完整合并前仍须完成

| 待接线 | 必须验收的产品合同 |
| --- | --- |
| ESF / TCR | 实际空间输运/梯度/Halo、持续随机场、两个化学半区间、TCR 分支/折点历史、共同接受/拒绝、完整重启。科学证据不足的 validated 模式保持明确拒绝。 |
| 喷雾 P4–P8 | 实际气相采样与几何/IBM 事件、液体资产身份、parcel 生命周期/TAB/注射器、跨 rank 迁移与源路由、质量源进入连续性和压力闭合、动量/总焓/动能收支。 |
| 共同事务和运行切换 | 全部模型状态与气相一次提交及回退；完整产品 Restart 和真实 CLI 验收后，再重建冻结程序并从封存 10510 做短续算验收，最后恢复圆柱长算。 |

本节点是接线中的可验证增量，完整合并尚未完成。
