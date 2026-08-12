# Task 19A 接受回执

结论：`ACCEPT`

基线：`f037fd7`

## 已交付

- 同一 `hundun` 可执行程序保持 schema v1/v2 行为，并新增 schema v3
  validate、canonical print-resolved 和 constant static LFP-GCIBM 正常运行。
- driver 在时间循环外一次性构造 STL query、immersed domain、ghost plan、
  wall quadrature 和 halo；每个 trial 只调用 Task 11 已接受的
  `FixedStepImmersedFlow`。
- committed step 强制检查恰好两次 PISO corrector、finite continuity/
  residual 和四字段 force；recoverable failure 使用同一 state 的事务回退后
  缩小 `dt` 重试。
- 19A 对 WALE、非 constant density、scalars、Checkpoint v3 restore 和
  performance mode 明确拒绝，不提前实现后续任务。

## 验证证据

- `git diff --check`：通过。
- focused Debug build：通过。
- v3 loader/header、1/2-rank broadcast、Stage 1/2 dispatch、19A 1/2-rank
  smoke、Task 11 signed-force/PISO/transaction：11/11 通过。
- 最终格式/错误文本修正后的 Stage 1/2 dispatch 与 19A 1/2-rank：4/4
  通过。
- Clang 15.0.6 + libc++ Release、`HUNDUN_BUILD_TESTS=OFF`、`-j32`：通过。
- 既有 `ib_ghost_stencil_plan.cpp` misleading-indentation warning 未变化，
  不属于 19A。

## 审查结论

- requirements：满足 19A constant IBM MVP；Checkpoint、完整 diagnostics、
  WALE 和三密度仍按权威计划后续实现。
- 数学/物理：没有新增离散算子；pressure、final flux、force 和 PISO 均消费
  Task 11 单一权威。
- MPI/事务：v3 config collective 广播复用既有 fault-sensitive 实现；新
  1/2-rank 产品 smoke 通过；Task 11 rollback/collective 证据保持有效。
- API/ABI/schema/Restart：只新增私有 app driver；公共头、Checkpoint 和字段
  identity 未变。
- 版权：只复用 HUNDUN-FLOW 自有 public API 和 accepted tests 的构造语义；
  未复制外部或私有源码，未新增依赖。
- caller impact：Stage 1/2 dispatch 回归通过；codegraphf 对未提交新符号未能
  建立节点，调用闭包由 `rg` 精确核对，提交后重新索引。
