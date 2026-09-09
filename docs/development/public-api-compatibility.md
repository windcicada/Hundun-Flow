# 公共 API 兼容性

公共 C/C++ 头文件位于 `versions/v0.4/include/hundun`。当前命名空间 `hundun::v04` 和 `v04_*` 名称为兼容性保留，不是退休源码。以下变化按兼容性变化管理：

- 删除或重命名公共类型、函数、枚举值和宏；
- 改变参数、返回值、异常语义、单位或符号；
- 改变公共结构体布局或 C ABI；
- 改变 JSON key、默认值和约束；
- 改变 Restart、diagnostics 或插件 schema。

私有头文件位于 `versions/v0.4/src`，其中 `_detail.hpp` 不构成兼容承诺。即便如此，私有重构也不得悄悄改变数值结果、MPI collective 顺序或持久化格式。

新增可选字段通常比改变既有字段安全。必须发生不兼容变化时，应提升对应 schema 版本，保留明确的错误信息，并在版本演进文档中记录迁移方法。

## v0.4 thermophysical physical-ghost certificate 迁移

GTMC G26 修复为所有输运模型交换 live/candidate 的 MPI/周期 k 与 k/cp，
通用混合物的物理边界仍使用原 EOS 闭合，COAST 的有效输运物理零梯度不变。
这是离散算子语义修复，不增加公共 API 或 Restart wire 字段；method-history
signature 增加 `generic-thermal-neighbor-material-v1`。旧 V3 的速率历史不得
静默当作新算子的精确 BDF 历史。本次开发验收从原始均值场重新导入 V1，
明确采用 BE 恢复；旧 V3 的其他迁移必须显式重建方法历史。

GTMC 开发修复增加显式 `BoundaryThermophysicalClosureKind::physical_inlet_face`
模式；默认 `ghost_state` 保留旧数值契约。新模式仅对完整固定组分、定焓、
压力外推的物理入口启用：p/h/Y 保持离散镜像，热物性由合法物理面状态计算。
输出边界 material 槽存放面值，T 槽仍为 `2*T_face-T_owner` 离散镜像。
`EquationPlanSpec::physical_inlet_material` 和 Cartesian kernel 的显式编译参数
必须同时启用，才能将这些系数作为面值消费；旧独立 kernel 默认为单元外延系数。
PISO / candidate 消费者检查两者契约一致。证书绑定新增模式和重建面值所依赖的
owner/source 数值，不能沿用只绑定 ghost 的旧解释。

`refresh_inlet_material` 用于 halo/零梯度及湍流更新后的入口系数重发布，不颁发
跨阶段密度证书；空 output view 不写入。有效黏度中的 SGS 增量从 owner 外推，
分子部分重新计算。内部物理状态、入口目标和 EOS 可行域不作截断或放宽。
以上新字段追加在原有 aggregate 后，旧源码仍可编译；需完整重编译，不可混用
新旧二进制布局。产品 method-history signature 已变更，不能把旧 V3 速率历史
作为此算子的精确 BDF 历史；原始均值转移需重新生成当前 plan 的 V1 恢复点。

`BoundaryGhostFieldAuthority` 保留原有 `field`、`revision`、`storage`、
`revision_domain` 四字段的顺序和语义；精确 `base`、`replica` 身份只追加在
其后。旧四字段 aggregate 初始化和五参数
`BoundaryThermophysicalFaceClosure::close` 继续执行原数值闭合，但不会发布可供
PISO 重用的证书。需要跨阶段消费时，生产者应使用
`make_boundary_ghost_field_authority` 构造精确 authority，并调用带 context 和
certificate 的重载。

PISO intermediate/state/terminal、pressure-correction 和 pressure-energy flux
certificate 均在既有字段之后追加 collective/rank-local lineage。手工 aggregate
初始化仍可编译，但缺少 lineage 时 `valid()` 会 fail closed；调用方应传播生产者
返回的完整 certificate，不应自行拼装 token。这些 certificate 不是 Restart、网络
或磁盘格式，因此本次变化不改变持久化 schema；若外部代码把公共结构体原始字节
当作 wire ABI，必须改为按字段序列化并显式版本化。

`PressureEnergyGlobalizationIterationReport` 在既有字段之后追加
`jacobian_scope_valid` 和 `jacobian_scope`。旧 aggregate 初始化继续成立，新增字段
默认为无效；消费者必须先检查 `jacobian_scope_valid`。该观测只存在于进程内报告，
不改变 runtime evidence、Restart 或网络格式，但混用新旧头文件的二进制需要完整
重编译。

## GTMC conservative transport integration

The GTMC merge retains native reaction/ESF/TCR/spray model and wire fields.
The imported-marker/patch envelope precedes the existing flagged case wire.
`EquationPlanSpec::conservative_total_energy` defaults to false; the product
currently enables it for inert, spray-free cases. The h-primary residual uses
rho(h+K)-p, with explicitly typed temporal K*rho derivatives and a declared
spatial quasi-Newton scope. It does not certify a complete convective-density
continuity Jacobian. Species enthalpy diffusion preserves the ESF Gamma-grad-h
route by avoiding an additional species correction there.

New material and historical-temperature halos, mixture/face authorities and
report members require rebuilding every client. Method-history fingerprints
bind the new arithmetic; create development V1 starts with the matching importer
or use an explicitly supported recovery policy. The original GTMC V1 source
was retained, and the mainline importer regenerated a matching V1 start with
zero field/flux readback difference. No restart manifest was edited.
