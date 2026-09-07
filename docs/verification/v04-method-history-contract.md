# v0.4 方法历史与普通应用入口

历史签名定义于 `core_product_freeze.cpp::method_history_signature()`，
只描述已保存状态、速率、面通量和时间历史的数值含义；物性、网格、边界参数
仍由 product plan / schema / geometry 绑定。完整 Git SHA 不参与历史兼容判断。

## 签名组件清单

| 组件 | 当前版本 | 约束对象 |
|---|---|---|
| history / bdf2-ex2 / rho-h-p / scalar-split | v1 | 历史格式含义、时间权重、状态量定义 |
| accepted-ibm-thermal-zero-normal | v3 | link 热重建使用流体侧物性，不借用固体占位导热系数 |
| thermal-inverse-representable | v1 | NASA 相邻可表示温度终止及浮点误差界 |
| stationary-ibm-placeholder | v1 | 无固体传热模型的占位状态不独立推进，非对流速率与固体内部面通量为零 |
| momentum-rates | v1 | 动量非对流速率 |
| simple-fresh-flux / c1-joint-target | v2 | SIMPLE 新通量、C1 联立目标 |
| open-periodic-flux | v3 | 开口及周期面通量 |
| periodic-metrics | v2 | 周期接缝度量 |
| momentum-afc-arithmetic | v4 | 动量通量校正算术 |
| conditional-boundary | v2 | 条件性边界分支 |
| scalar-paired-mass-remap / composition-picard / mass-roundoff-closure | v1，仅有标量 | 预测库存配对、组成迭代、连续性精度 |
| passive-envelope / composition-inner-accuracy | v1，仅有标量 | 被动标量范围、组成内迭代精度 |
| physical-donor | v2，仅有标量 | 物理边界 donor 权威 |
| ibm-scalar-impermeable-flux | v1，仅有标量 | 二值控制体切面零标量通量 |

改变任一组件的离散含义时更新相应组件，而不是替换完整签名为当前 Git SHA。
容量修复、无效对象保护、初场入口和不改变有效状态算术的观测修改不更改组件。
新增方法必须补精确续算拒绝、显式恢复、恢复后再精确续算的链式回归。

## 已知迁移表

| 来源 | 默认策略 | 显式策略及限制 |
|---|---|---|
| V1，确实没有时间/速率历史 | 既有缺历史 BE 恢复 | 不伪装为完整 V2/V3；按源格式处理 |
| 完整 V2，没有方法签名 | 拒绝当作同方法历史 | `rebuild_method_history` 重算当前方法速率，保留源质量目标，先 BE 后 BDF2 |
| 完整 V3，签名一致，plan/schema/geometry 一致 | 精确续算，继承有效 BDF2 历史 | 主动恢复仍是独立显式选择 |
| 完整 V3，方法签名不同，但物理 plan/schema/geometry 一致 | 拒绝静默复用旧速率 | 显式恢复重建速率，不能覆盖回源旧速率 |
| 已登记 AFC v3 旧方法 plan | 拒绝 | 只允许 `compatible_method_plan` 所列旧 plan，strict storage；不能任意豁免物理配置 |
| 已登记 MG bundle ghost v1 旧存储 | strict 默认拒绝 | 独立 `mg_bundle_ghost_v1` 存储迁移；不得与方法恢复同时开启 |
| 未知 plan/schema/geometry 或损坏文件 | 拒绝 | 方法恢复不跳过文件长度、校验和、字段、几何及物性完整性检查 |

## 普通应用与 runner

普通 `ApplicationRunOptions.restart_history_policy` 默认 `require_compatible`；
普通 `hundun run` 新增 `--restart-method-recovery`，必须同时给 `--restart`。
重复恢复开关、无来源恢复以及存储/方法迁移混用均拒绝。
Evidence 记录来源格式、来源/目标签名及恢复策略。专用 runner 的统计新 epoch
规则不变：显式方法恢复不能继承旧样本，一个 BE 步不代表流场已充分发展。

普通应用允许 `ApplicationRunOptions.initial_state` 指定初场（t=0，标量借用值在同步
run 调用内存活、按冻结目录排列）；与 restart 互斥。缺省只接受一致的边界提示值，
冲突返回 initialize/10505，不能由最后一个边界决定温度/速度/压力参考。
CLI 使用 `--initial-state p,T,Ux,Uy,Uz[,q...]`（Pa、K、m/s，逗号分隔）指定
均匀初场，q 按冻结标量目录排列。标量数量、组分范围和 EOS 仍由 ProductDriver
验证；非有限值、缺少分量、重复选项、与 restart 混用均拒绝。暂不支持非均匀初场文件。
ProductDriver 允许有限有符号被动标量；只有组分施加质量分数范围及组成/EOS 检查。

## 尚未自动化的时间尺度

普通应用已用已接受流场收紧对流限制。`LocalTimeLimits` 其余分量仍由调用者按所用
离散算子提供，并非自动声学/黏性/热/组分稳定性估计。没有把所有扩散都当显式处理。
自动非对流尺度是独立方法设计工作，需分别定义显式稳定性、预测器可接受性和隐式精度
约束，之后用相应算子回归验收；不能靠填入一个看似物理的通用常数完成此项。
当前固定 dt 的 Re3900 不使用这些候选尺度，也未更改固定 dt。
