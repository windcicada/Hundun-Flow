# 公共 API 兼容性

公共 C/C++ 头文件位于 `include/hundun`。以下变化按兼容性变化管理：

- 删除或重命名公共类型、函数、枚举值和宏；
- 改变参数、返回值、异常语义、单位或符号；
- 改变公共结构体布局或 C ABI；
- 改变 JSON key、默认值和约束；
- 改变 Restart、diagnostics 或插件 schema。

私有头文件位于 `src` 并以 `_detail.hpp` 结尾，不构成兼容承诺。即便如此，私有重构也不得悄悄改变数值结果、MPI collective 顺序或持久化格式。

新增可选字段通常比改变既有字段安全。必须发生不兼容变化时，应提升对应 schema 版本，保留明确的错误信息，并在版本演进文档中记录迁移方法。

## v0.4 thermophysical physical-ghost certificate 迁移

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
