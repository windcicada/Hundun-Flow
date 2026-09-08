# Re3900 稳定版本长测交接

本轮完成长测准备与启动，不是 35000 步稳定性、统计充分发展或 COAST 替代验收。
两相仍暂挂；网格、原生变物性、固定 dt、残差要求、refinement 容量和 checkpoint
持久化合同均保持。没有把失败的性能实验带入生产。

## 候选与回归

生产源码来自 `7c03a54a9f1061d1de0e0922ed38a3bcc8f2412f`，干净 Release 的
31/31 回归与 128-rank 单轮 9500→9510 验收见
[恢复原因观测报告](2026-09-08-product-fgmres-recovery-observation.md)。
冻结 runner SHA-256 为
`345cca4802286a7ae1b2bd39b7c25afceeedb8a13fc60b5fdf71f80d92f0f67a`。
manifest 为 `c352d1c4ff748b5ed33d3d068115d8de54f3b8e61e5efe70cfb3ffb10a3b2e21`。

[后续列保留实验](2026-09-08-fgmres-salvage-trial.md)被局部反例拒绝：合法的小缩放
非对称算例原路径 44 次收敛，候选到 80 次仍拒绝。实验源码已撤回，保留公开接口
保护回归；撤回后 Release 与 ASan/UBSan 各 1/2/4 ranks 均通过。
`d706a41f37801553de29bc52a5ce85209a10aeee` 的生产 src/include、runner 与 observer
和 `7c03a54` 无 diff；该提交新增测试与证据，不冒充重新构建了生产程序。
GitHub main 已正常快进到 `d706a41` 并通过远端引用回查；未覆盖未知并发提交。

## 交接过程和实际证据

基准根目录：
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52`。
证据根目录为其下 `long-handoff-20260908`，紧凑记录归档于
[data/2026-09-08-long-handoff](data/2026-09-08-long-handoff/)。
脚本固定绑定该冻结目录，不宣称可在任意目录直接重放。

| 阶段 | 实际结果 |
|---|---|
| 旧 checkpoint 保存 | 原 9000/9500 两个 generation、current、统计与 accumulator 另存 `source-before-resume`；全部文件 SHA-256 比对一致，不移动或手改原件 |
| 内存进度交接 | 2026-09-08 13:59:09 +08:00 恢复旧 128 ranks；14:05:21 检出完整 10000 checkpoint 后重新 SIGSTOP。9951–10000 共 50 步全部 BDF2、attempt=1、无 retry |
| 物理检查 | 50 步最大连续性 `2.657905180139039e-7`、能量 `9.992124619488932e-7`，均在原门槛内；区域数量正确，固体占位极值未漂移。不是新增全载荷固体比对 |
| 持久点 | `generation-10000-171280289222885`，t=`0.13808912271981852 s`；128 rank 文件、manifest、complete、statistics、accumulator 和 128 Visit 文件齐全。另存 `source-10000` 并逐文件比对、同步；manifest 完整性通过公开验证器 |
| 旧作业结束 | 保存并验证交接点后 administratively stop；旧 unit inactive、MainPID=0，所有旧 MPI 进程退出。不把停止的日志前缀说成流关闭验收或 COMPLETED |
| 新长测启动 | 14:08:32 +08:00，`hundun-re3900-observed-long-20260908.service`，新目录 `long-observed-35000-20260908`。128 ranks，精确续算 10000→35000；无其他 CFD/MPI 或编译测试竞争 |
| 恢复合同 | V3 source/target signature 同为 `12213963202598979269`，policy=`require_compatible`，无 method recovery；源 rank 内容由真实 Reader 校验后进入推进。启动后源 10000 文件 SHA-256 未变 |
| 统计 | epoch=7000，development=10000，sampling_start=17001，旧样本数=0；未重新清空统计。元数据中的 reset_reason=method_recovery 是继承的 epoch 来源，不代表本次恢复又执行了 reset |
| 启动健康前缀 | 10001–10016 共 16 个完整记录全部 BDF2、无恢复/重试；最大连续性 `2.0864669393733558e-7`、能量 `9.576842172929002e-7`，固体区域极值不变 |
| Evidence / 观测 | 启动 Evidence 对源 10000 manifest 的真实 CLI validator 返回 0。V6 observer 校验 17 步、102 loops、128 ranks，`complete=false`，唯一 issue 为预期步 10018 尚不完整 |

运行参数沿用 `456×256×52` 网格，`D=0.02 m`、`Uc=2.89668 m/s`、
dt=`1.3808912271980336e-5 s`，目标 35000；checkpoint/Visit 间隔均为 500。
启用同验收 pilot 的 performance、MG cost、FGMRES recovery 观测，cell trace 关闭。
互斥锁覆盖计算与 maintenance 槽；未同时启动两个 128-rank 作业。

## 证据边界与下一检查点

- 原 9951 步计时包含 SIGSTOP 间隔，不能用于吞吐比较；本轮不声称新加速。
- 启动前缀从运行中的文件只读复制，各流可能停在不同步。健康记录完整到 10016，
  观测完整到 10017；不裁改 RUN.meta 中的 35000 目标来伪造 complete。
- 原作业的正常 writer 在 10000 发布时更新 current 并按两代合同清理旧代。
  旧收据仍描述当时的源路径/hash，不重写；9000/9500 原内容保存在另存副本。
- 新运行首次 checkpoint/Visit 应在 10500；在此之前可恢复点是只读 10000 源。
  后续核对新 checkpoint/统计附件、RSS 趋势、残差/工作量与输出完整性。
- 本轮未运行到 35000，未完成实验统计对比、COAST 等窗口性能验收或燃烧接入。
  待长测可用证据形成后再选一个有依据的优化；不恢复已撤回的实验或放宽物理门槛。
