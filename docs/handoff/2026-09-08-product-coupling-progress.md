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

## 原生气相采样、周期模板与角点数据

- `ProductParcelGas` 适配器读取原生 MeanState 的 p/h/U/Ns−1 字段，按完整组分索引映射
  构造 PH 样本；实际网格面坐标用于 parcel 单元/owner 定位。过期修订、无效组成或超出
  可用 Halo 的采样明确不可用，不部分覆盖调用者的完整组分缓冲。
- 共同三线性模板新增可选周期轴，使用真实端部单元中心间距跨周期边界插值，保留全局
  单元排序、非负权重与总和一。物理域退出保持独立判断，关闭周期时原模板行为保留。
- 4-rank 实测发现主 Halo 只填六个面，二维分区的角点仍为 NaN。新增显式目标入口复用
  `RemoteDonorExchangePlan` 补齐边/角点，未改主 Halo 或既有 quadratic 入口的身份路径。
  新入口冷阶段核对所有 rank 的字段顺序/组分宽度、域尺寸和周期属性，防止错序消费。
- 17×11×7 非均匀网格的真实 Halo+角点收集+原生场采样在 1/2/4 ranks 通过；热阶段
  通信/采样无 C++ 动态分配。原始失败与定位证据见 `spray-native-gas-mpi-diagnostic.log`。
- FilmEnvironmentBridge 可重绑目标修订并复用既有液体/气膜工作区，实际 A–S 区间在
  第二个已接受步继续推进；旧修订仍拒绝。无需逐步重建和分配气膜工作区。
- 最终 Clang 11/11 PASS，GCC11/Cantera 9/9 PASS，包括原 IBM donor 计划、压力权限、
  原生喷雾采样/迁移与 ESF/TCR CLI。日志 `spray-sampling-final-{clang,gcc}.log`。
- 这些是 native product 所需的采样/通信适配器；ProductCompiler 的 spray 拒绝仍保留，
  尚未共同提交 parcel/注射器，也尚未发布喷雾气相源。完整合并和圆柱恢复仍未执行。

## 喷雾与 TCR 的组合持久状态

- 原生 `ProductSprayHistory` 将 parcel 完整整数 ID、TAB、破碎序号和注射器余量/序号按
  所在全局单元编码为 Restart V5；可同时携带每单元 120 字节 TCR 分支历史。
- 注射器恢复与 TCR 字节恢复均先暂存，统一 preflight 后才发布；迟到的 TCR 身份错误
  可撤回全部暂存状态。热阶段编码、恢复与提交不分配，初始容量计入所有自有缓冲。
- 实际 RestartWriter/Reader 的 1→4、4→1 恢复逐字节核对组合状态，覆盖大于 2^53 的
  ID/序号。原 ESF/TCR CLI 路径同步复核。Clang 13/13、GCC11/Cantera 14/14 PASS，
  日志 `spray-history-final-{clang,gcc}.log`。
- 父运行时仍须核对恢复位置/液体身份、全局 ID，并将这些暂存状态接入气相唯一事务。
  当前仍未启用 spray 产品配置，完整合并和圆柱恢复尚未执行。

## 按实际单元 owner 路由喷雾交换

- 新增有界 `OwnerExchangeRoutingPlan`：先把同一 parcel ID/段序号的完整模板送往确定的
  auditor 核对归一性、重复单元和一致预算，再仅发送到实际单元 owner。使用 Alltoallv，
  不收集全局单元表或复制全部交换段；17 个整数/IEEE 位模式载荷保留完整整数身份。
- 发现并修正跨 owner 结算问题：owner 收到的是部分模板，不能再次要求本地权重和为 1。
  只有路由计划生成的修订/代际绑定视图可进入分布式批处理；普通批处理仍要求完整模板。
  任意 prepare/discard 或成功重配都使旧视图失效，失败冷重配保留旧计划。
- 已验证质量/动量、完整组分和在累计后计算的气相动能/总焓改变量；wall 账本独立，
  不投到气相。接收/auditor 容量不足、重复跨 rank 段、坏权重、非有限预算、过期修订
  均集体拒绝；空重试可继续。自有缓冲精确预算边界和热阶段零分配通过。
- 批处理按已排序单元/交换段线性遍历，避免对每个本地单元扫描全部交换段。
  Clang 12/12、GCC11/Cantera 13/13 PASS，包含原 portable composition 和真实 CLI 续算。
  日志 `spray-owner-{red,green,final-clang,final-gcc}.log`。
- 运行时轨迹推进与气相场写入尚待接入；本增量不解除 spray 产品拒绝或圆柱暂停。

## 实际网格与 IBM 事件查询

- `ProductParcelGeometry` 使用真实非均匀网格面、周期属性与已有静态 IBM surface 查询，
  为自适应轨迹段产生换单元、壁面和物理出口事件。周期轨迹保留展开坐标，采样/最终
  定位按同一周期映射；单单元周期轴不制造无意义的 owner 改变。
- 原事件积分器只支持零时长破碎，实际面上反向运动因此失败。现支持零时长拓扑变更、
  出口和壁面冲量，不创建虚假的交换时段；重复事件仍受事件总容量限制。
- 验证覆盖周期跨界、同一时刻两个面换单元、精确面起点、独立出口/壁面账本、IBM 与
  单元面重合时的优先级、非均匀坐标和过期修订。Clang 7/7、GCC11/Cantera 8/8 PASS，
  包含原 portable composition 与 CLI 回归，日志 `spray-events-native-final-{clang,gcc}.log`。
- 这是轨迹段的几何适配器，运行时还须协调子步/跨 rank 迁移和域外预测采样；它本身
  不启动 parcel 推进或写入气相源。完整合并与圆柱恢复仍未执行。

## 原生 parcel 推进与交换结算流程

- `ProductParcelAdvance` 借用已绑定的原生气相场，串起实际气膜/A–S、TAB 与子滴剩余时间
  推进、注射器暂存、分段插值投递、按 owner 迁移、完整出生 ID 审计及交换结算。
- 同步推进按全局速度与最小网格宽度限制子区间；不可用的局部预测可有界缩短区间重试。
  每个子区间结束迁移，避免逐 parcel MPI 查询，也无需将尝试内时间写进永久 parcel ABI。
  汇总全部交换段后才计算气相动能修正。壁面、出口、外部注射和 TAB 能量库分别记账。
- 1/2/4 ranks 的实际原生场 + 液体资产 + A–S/TAB 测试覆盖蒸发、跨 owner/周期迁移、
  破碎子滴续进、质量与 H_tc+K 收支、整数序号；热阶段无 C++ 动态分配。
  单 rank 后期 PH 查询失败撤回所有候选并回退注射器，重试与控制结果一致。
- Clang 12/12、GCC11/Cantera 13/13 PASS，包含原采样、portable composition 与真实
  ESF/TCR CLI 回归；日志 `spray-advance-final-{clang,gcc}.log`。GCC 暴露的测试 `<climits>`
  直接依赖已补齐；首次 Make 自动重配后新增目标需重新调用 build，未复用失败构建结果。
- 当前推进测试使用周期原生场；实际非周期域外预测采样仍需产品边界约束。此流程尚未
  绑定 ProductDriver 气相源字段或共同提交点，喷雾产品入口仍关闭，圆柱仍暂停。

## 原生平均气相与喷雾接线（2026-09-09）

- 新增产品所有者 `ProductSpray`：实际液体资产、按原生 owner 配置的注射器、
  全采样 fringe 交换、A–S/TAB 推进、迁移/源路由和原生 V5 历史已接入 ProductDriver。
  当前准入为周期、BE、平均气相反应；ESF+喷雾、非周期和 IBM 仍明确拒绝。
- Sm、N−1 蒸气组分、动量和总焓源进入实际预测器与方程装配；Sm 同时进入连续性、
  C1/C2 和最终审计，封闭体总气相质量目标仅随共同事务提交。
  总焓源使用当前 stage 2，进入候选能量残差，避免被存入 EX2 历史重复计算。
  独立 128 W/m³ × 单元体积残差 oracle 与完整装配逐位相同。
- 原生共同回退测试在气相后期 PH 查询注入失败，已接受字段/通量/质量目标及 V5
  记录逐字节不变；重试气相按 2e−12 求解精度与控制比较，喷雾记录逐字节一致。
  压力暖启动在拒绝后清除，因此新求解不要求气相浮点结果逐位相同。
- 1/2/4 ranks 实际产品三步算例验证总气液质量误差 <1e−10 kg、动量误差
  <1e−9 kg·m/s、总能量误差 <1e−5 J；观察能量误差约 6.47e−7 J。
  能量 oracle 为积分 rho*h − p + rho*|U|²/2 加液滴 H+K，注入账单独计入。
- 真正的 RestartWriter/Reader/initialize_restart 检查 V5 完整恢复；有限正直径但
  与液体密度/质量不符的记录在原生恢复预检拒绝，再恢复原文件成功。
  CLI 新算两步后 1→4 ranks 续算两步，证据文件及 manifest 对照校验通过。
- 修正产品初始化入口遗留的 V4 上限，以及 C++/Python 证据读取器的 V5 上限。
  增补源字段的计算图读依赖，并在 arena 分配前累计喷雾所有者和 arena 的容量。
  全产品其他冷服务/模型的统一 resident budget 仍需后续审计。
- Clang 25/25、GCC11/Cantera 32/32 相关原生回归 PASS，日志
  `build-coupling-evidence/spray-product-final-clang.log` 与
  `build-integrated-cantera/spray-product-tests.log`。
  这只是平均气相接线节点；P8 的源后随机场/输运/TCR 顺序仍须完成。主线未推进，
  PID 165884 的圆柱服务继续 SIGSTOP，后续仍从封存 10510 进行新版本验收。

## P8 原生共同更新（2026-09-09）

- `ProductSpray` 现已接受周期 BE 的 ESF/TCR 组合。先把同一气液交换按旧/新气相
  质量沉积到每个随机场，再使用独立原生工作场和完整 donor fringe 交换。
  更新后的 PH 状态生成共享 Gamma/Dt，随后执行空间输运、TCR/IEM、两次化学半区间。
  零交换单元保留原字段值，避免无意义的乘除舍入。
- `PredictorTransportState` 让气相输运读取源后 h/Y/被动标量，同时保留时间导数中
  原 accepted rho*q；当前 Sm/组分/总焓源仍只加入一次。接口仅接受相同源身份的 BE，
  校验时间、存储、Halo、有限性和输出别名。独立线性梯度 oracle 经 RED→GREEN。
- 有符号被动标量按新库存稀释并使用自己的分子/湍流 Schmidt 数，气液总质量、
  动量、总能量与无外源标量积分均有原生 1/2/4 ranks 检查。
  TCR 首步 eta 在恰好八个实际沉积单元升高，证明读取了源后的统计量。
- 修复化学第二半区间 rank-local 失败后，其他进程提前进入 parcel-history collective
  的错位：进入共同历史阶段前先归一失败状态。MPI 故障测试不再超时，已接受气相、
  随机场、TCR、粒子、注射器和通量全部撤回；原值恢复仍逐字节核对。
- 重启后继续求解采用已有原生反应测试的无量纲尺度比较。TCR 的六个浮点历史量按
  数值容差比较，整数/粒子记录保持精确；避免把约 1 ULP 的 R 变化当作整数字节错误。
  两种故障均核对重试结果，并保持 2e−12 数值容差。
- 源身份/时间的一致性合并进已有 high-state 集体通信，保持快路径一次通信的证据
  契约。更新方法历史签名、原生 arena + ESF/TCR/喷雾 owner 内存界限及采样 Halo 预算。
  采样 U/p 使用独立 donor fringe，不冒充主 face-Halo 包内的字段。
- 实际 CLI 覆盖平均喷雾、ESF/TCR 喷雾、负值被动标量和真实 Cantera ESF/TCR 喷雾：
  两步新算、V5 检查点、1→4 ranks 再续两步及证据/manifest 校验。
  Clang 40/40、GCC11/Cantera 48/48 PASS，日志 `p8-clang-final.log` 和
  `build-integrated-cantera/p8-final-tests.log`。GCC 日志位于该构建目录；Clang 日志在
  `build-coupling-evidence/`。最初扩大回归时漏建 ghost-authority 测试目标，补建后通过。
- 完整合并仍未完成：非周期/IBM 的粒子采样边界延续和实际产品绑定仍未准入；
  ESF 非周期/IBM 也仍拒绝。未推进主线，未恢复圆柱。

## 原生粒子物理边界（2026-09-09）

- 确认并复现 A–S 试探终点先采样、后裁剪事件的失败原因。气相采样新增明确启用的
  一单元宽有界延续，使用既有非负单侧插值；域外点仍不能成为沉积位置或迁移目标。
  精确落在上壁面的已反弹粒子归属相邻内部单元，出口粒子在定位前移除。
- ProductSpray 从实际成对周期边界和静止壁面配置事件；完整 donor fringe 仅交换
  有效域内单元，保留原生物理边界 ghost。移动壁面仍拒绝，现有事件账仅支持静止壁面。
- 实际 A–S 在 1/2/4 ranks 验证跨单元、壁面反弹、出口移除、精确面起点、质量/动量/
  总能量与外部账闭合，且热阶段无动态分配。壁面 ProductDriver 检查实际反弹方向、
  总质量/能量、后期 PH 失败共同撤回以及真实 V5 恢复；壁面外力存在时不使用封闭
  气液系统的零外力 x 动量 oracle，该方向冲量由 A–S 事件账单独验证。
- 壁面与贯流出口 CLI 均完成两步新算、1→4 ranks 续两步、证据与 manifest 校验。
  首个静止气体出口配置触发现有 1713 回流方向保护：局部源造成单元速度与面通量
  方向不一致。无喷雾对照成功；有明确入口/初始 0.3 m/s 的贯流配置成功，保护未修改。
- Clang 32/32、GCC11/Cantera 34/34 PASS；新增精确面起点两种路径各 1/2/4 ranks，
  两个工具链分别再通过 6/6。日志 `build-coupling-evidence/spray-boundary-final-clang.log`、
  `spray-boundary-exact-clang.log`，以及 `build-integrated-cantera/spray-boundary-final-tests.log`、
  `spray-boundary-exact-tests.log`。临时嵌套失败码诊断已撤回，测试保留候选拒绝原因输出。
- IBM 仍需流体侧公共插值/沉积、静态几何实际绑定；ESF 非周期/IBM 仍未准入。
  主线未推进，圆柱 PID 165884 仍为 Ts；不把本节点称为完整合并。

## 原生静态 IBM 喷雾（2026-09-09）

- ProductSpray 在原生 STL/EBTopology 编译后绑定实际表面与流体侧。单独的冷阶段
  donor 交换封存二进制流体掩码，包括 MPI 周期角点；插值和沉积共用同一组可见流体
  donor 及归一权重，固体状态不参与 PH 查询。试探采样允许一单元宽表面投影延续，
  域内定位和实际沉积仍拒绝固体位置。相关内存计入冷阶段容量预检。
- 16³ 静态立方体实际产品在 1/2/4 ranks 完成三步蒸发/反弹、气液质量和总能量检查、
  后期 PH 故障共同撤回、V5 恢复与重试；CLI 完成两步新算和 1→4 ranks 再续两步。
  第二步观察质量误差约 −7.1e−15 kg、能量误差约 −1.55e−9 J。
- 实际第三步复现了微小子步的大背景库存相减消去误差。质量、H、动量差改为先求
  状态差再乘 multiplicity 的长双精度表达式；独立微小库存 oracle 经 RED→GREEN。
  自适应 H/动量误差判据加入由端点质量、U、T、h 的 ULP 和液体 cp 推导的舍入误差
  下限，配置的截断误差容差不变；错误热预算及虚假动能增量仍被拒绝。
- 几何事件时间容差曾允许已反弹终点进入固体约 6.6e−13 m。事件现携带实际接触点，
  仅在位置绝对容差、时间定位误差及 ULP 界限内修正位置；超过界限则缩步再积分，
  不改变质量、速度、热力学状态或交换账。壁面与 IBM 精确接触用例经 RED→GREEN。
- Clang 72/72、GCC11/Cantera 74/74 PASS，覆盖 portable 事件/组合、实际平均及
  ESF/TCR 喷雾、物理边界、IBM 与 CLI 重启。日志位于
  `build-coupling-evidence/spray-ibm-final-clang.log` 和
  `build-coupling-evidence/spray-ibm-final-gcc-tests.log`。临时诊断均已撤回。
- ESF 非周期/IBM 的边界闭合与算子修正仍须完成；主线未推进，圆柱仍暂停。

## 完整合并前仍须完成

| 待接线 | 必须验收的产品合同 |
| --- | --- |
| ESF / TCR | 实际空间输运/梯度/Halo、持续随机场、两个化学半区间、TCR 分支/折点历史、共同接受/拒绝、完整重启。科学证据不足的 validated 模式保持明确拒绝。 |
| ESF 边界 | 平均喷雾的静态 IBM 已验收；仍须接入随机场的物理边界闭合、IBM 零通量扩散和固体单元处理。 |
| 共同事务和运行切换 | 全部模型状态与气相一次提交及回退；完整产品 Restart 和真实 CLI 验收后，再重建冻结程序并从封存 10510 做短续算验收，最后恢复圆柱长算。 |

本节点是接线中的可验证增量，完整合并尚未完成。
