# Re=3900 三维圆柱绕流 Stage 4 融合验收报告

```text
schema=hundun.cylinder_re3900.stage4_acceptance.v1
decision=ACCEPT
accepted_product_head=9c18c7ad4a084dd3d9a5e481d2d2a3af52b1990c
accepted_product_tree=cc5b7be1302a42b3930d20642edb92db2ca661df
stage4_governance_base=033a685c90c1f9c674e93a4b82db10db4c381abe
stage4_tested_code=6407cd7c591ce088db7f1dd7e296d77acd18da1c
stage4_to_product_diff_sha256=1fe5674c8f08d08594d578cf3ff9f4d8a5273b15857bbf9eeb51afb23667f1ca
governance_receipt_commit=commit_containing_this_report
long_time_statistics_status=DEFERRED
```

## 1. 结论和能力边界

接受上述 product candidate。它已完成 Stage 4 治理基线与周期圆柱 IBM
开发线的融合，并通过低成本功能、`1/2/4-rank` MPI、力诊断、
`48³/64-rank` fast 和冻结 `480×480×48/64-rank` 一步启动门。

本结论证明静止周期 IBM、WALE、两次 PISO、精确 Schur 压力通路、MPI
分解和四字段力报告可以在冻结工程网格上共同完成一个 committed step。
它不证明时均阻力、脉动升力、Strouhal 数、尾流剖面、湍流统计或网格
无关性已经收敛。这些长时统计不属于本功能门。

## 2. exact-candidate 身份

- 代码 HEAD：`9c18c7ad4a084dd3d9a5e481d2d2a3af52b1990c`；
- tree：`cc5b7be1302a42b3930d20642edb92db2ca661df`；
- Stage 4 governance seal 到 candidate 的 binary diff SHA-256：
  `1fe5674c8f08d08594d578cf3ff9f4d8a5273b15857bbf9eeb51afb23667f1ca`；
- tests-on Release `hundun` SHA-256：
  `13eb07cc40e056bd6151b7cbf19cb7b12febaf5a1cda00fbd97663b94c3e06c0`；
- tests-off Release `hundun` SHA-256：
  `70e5324b8cd864a55177c09011f964d566a8ddf8622c280f657a311f7b86643b`；
- full-grid case SHA-256：
  `8f3066e5d58e530d0d22ec5c9033360743cd33ff0785cef36a13a93d70704a4a`；
- STL SHA-256：
  `bd264c586543de4ec330f53cb2d3d9dfba550823db69451ac3176c603d248f46`；
- 验收前后 `worktree-status` 均为空，SHA-256 均为
  `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`。

数值测试完成后没有修改上述 candidate。本报告及执行计划收口只属于后续
治理提交；该治理提交不冒充经过数值验收的产品 HEAD。

## 3. 低成本和 MPI 证据

| 证据簇 | 结果 | 证据 |
| --- | --- | --- |
| Release focused | `10/10 PASS`：FGMRES、immersed PISO、transaction 的 `1/2/4-rank` 及 public header | 外部日志 `ctest-release-focused.txt` SHA-256 `b207844d3bfb45f4096e09d611d78514cda127c67e9fab08402c199be6612f20` |
| Debug focused | `10/10 PASS`，同上矩阵 | `build/cylinder-debug/Testing/Temporary/LastTest.log` SHA-256 `1b54ad88e7e8cc1e212362439fb414d20865e7d893eb16b8f28c5828076ed563` |
| UBSan focused | `4/4 PASS`：FGMRES、PISO、transaction、header | `build/cylinder-ubsan/Testing/Temporary/LastTest.log` SHA-256 `31f40bf3939d0cc50d244ad42fbb4c833616b3b8d760339527dd45f205aef72c` |
| ASan focused | 关闭已知 OpenMPI/OpenPAL/PMIx 退出泄漏检测后 FGMRES 与 PISO `2/2 PASS`；transaction 在先前开启泄漏检测的执行中已到达功能终点，失败只含 MPI 运行时退出泄漏，没有 HUNDUN 的 OOB/UAF/double-free 栈 | `build/cylinder-asan/Testing/Temporary/LastTest.log` SHA-256 `cfc36fb8000f97eef681caf108ceae2518ec06d698335bb90435660008963b08` |
| exact 48³ tests-on | `PASS`，33.20 s | `runs/hundun-fast48-9c18c7a-tW3DA8`，`run.log` SHA-256 `d5bd428850bfef3b92e6f0acff853ed3815713dfc6f10ca99f53927fb9e3a40f` |
| exact 48³ tests-off | `PASS`，32.08 s | `runs/hundun-fast48-tests-off-9c18c7a-GffkqO`，`run.log` SHA-256 同上 |

tests-on 与 tests-off 的 `48³` 产品输出逐字节相同，证明移除测试访问宏
不改变数值路径。分解矩阵覆盖稀疏 active halo、精确 Schur 齐次/仿射
权威、两次 corrector、rollback 和 collective failure。

## 4. 冻结全网格启动门

证据目录：
`/home/wyf/code_dev/.benchmarks/cylinder-re3900-stage4/runs/hundun-fullgrid-tests-off-9c18c7a-s3QpZa`。

- 配置：`480×480×48`，64 ranks，`4×4×4`，`dt=0.006 s`，WALE，两次 PISO；
- 执行：`attempts=1`，`correctors=2`，`exit_status=0`，`terminated=0`；
- 时间：26 min 51.50 s；
- continuity：`7.4994793093741969e-15`；
- pressure residual：`1.7799043776016363e-11`；
- `force_operator=(1272.1640409893714,-0.042943683757230941,-0.11973548955862613)`；
- `force_budget_reaction=-force_operator`，`force_surface_traction` 和
  `force_consistency` 全部有限；
- WALE identity 完整；初始均匀速度场的第一步
  `nu_t_min=nu_t_max=0`，`zero_count=172800`；
- 资源采样的 64-rank 聚合 RSS 峰值 `139338260 KiB`，即
  `132.883320 GiB`；最低可用内存 `119943020 KiB`，watchdog 未触发；
- `manifest.txt` SHA-256：
  `e0421cf2d50ddf3df60d95d4e1d93bffd2dbdbbe45a6613eec1a0612201c8b62`；
- `run.log` SHA-256：
  `a37d38e3c36c752acee6a4222d9994b16aa2752131a2e668f59afbfa9a0b5fc5`；
- `time.txt` SHA-256：
  `d09e1bfa5cdf867c3ae8f08209753d9204a2d9969bdbbeb8492a179aa6dcb58c`；
- `resource-samples.tsv` SHA-256：
  `928ac7f2a4edafd29a97d25dfd32f467862a4756e06f6094c24f7c113bbff97e`。

按 `rho=U=D=1`、`Lz=pi D` 和运算符力计算，该瞬时启动步为：

```text
Cd = 809.8847821888761
Cl = -0.027338798178153764
```

`Cd` 的大幅瞬态值对应静止初场的脉冲启动；本门只要求符号、量纲和
输出链完整且数值有限，不将它与文献的稳态或时均阻力直接比较。

## 5. 内存故障区分

同一 candidate 的第一次全网格使用 tests-on 开发构建，在约 17 分钟
后被资源 watchdog 停止。根因是 `HUNDUN_FLOW_ENABLE_TEST_ACCESS=1`
使产品可执行程序包含测试专用的全活动单元宽 `MPI_Allgatherv` oracle，
并在每个 rank 复制约 1100 万单元的宽记录。该 oracle 不在
`HUNDUN_BUILD_TESTS=OFF` 的正式产品构建中。

正式 tests-off 构建不含该宏，且它的 `48³` 输出与 tests-on 构建逐字节
一致。全网格 tests-off 作业的聚合 RSS 在旧失败时间窗口内保持稳定，
最后正常退出。因此旧 OOM 属于误用 tests-on 开发二进制的测试 oracle
开销，不是正式 product path 的数值内存增长。

## 6. 完整 diff 审查

- **数学和物理：** 受控 RED 证明精确 immersed Schur 是非对称算子；
  外层使用 restarted right-preconditioned FGMRES，紧凑压力算子只作不精确
  inverse，没有替代精确权威。齐次响应为 `F(p,0)-F(0,0)`，仿射 RHS
  保留 `F(0,g)`。壁面 `mu_eff` 在冻结 donor 正值范围内有界。
- **离散和一致性：** 周期分割 STL 维持单一 parent-surface 语义；active
  pressure/momentum 通路使用稀疏 peer halo，不在产品 runtime 做全局
  active-value gather。压力、最终通量和力保留同一 authority，PISO 仍恰好两次。
- **API/ABI/schema：** 新 FGMRES 公共头通过 standalone contract；已有
  失败阶段 enum 的 `0..4` 序号冻结，只追加
  `pressure_compact_preconditioner=5`。没有改动 Restart/schema/field identity。
- **ownership/lifetime/rollback：** restarted FGMRES 持久拥有 Buffer
  workspace；相同 layout 重复求解不再分配 Buffer，layout 替换失败保留旧
  workspace。内层失败的 phase/component/reason/report 可确定转换，retry
  与 rollback 语义不变。
- **MPI：** 稀疏 layout 检查 duplicate、locally-owned ghost、
  missing/inactive ID 和 collective disagreement；`1/2/4-rank` 覆盖运行和失败路径。
- **caller/build impact：** 压力 solver 只在 immersed-flow driver 选用
  FGMRES，SPD 动量路径仍用 CG。正式基准使用 tests-off 产品二进制，无
  Python 运行时依赖。
- **版权和范围：** 变更代码中未发现私有 COAST/BOFFIN 路径或上游项目
  标识；未访问私有源码和研究数据，未引入 vendor solver 或新运行时依赖。

## 7. DCO、进程和最终状态

`033a685c..9c18c7a` 共 12 个提交，author/committer 均为
`WANG YUDONG <wangyudong@buaa.edu.cn>`，每个提交均已有同一
`Signed-off-by`。本次只验证已有 DCO，没有伪造或为旧历史自动添加 sign-off。

验收结束时，全网格主 unit 和 watchdog unit 均为 `inactive/dead`，退出码
均为 0；没有遗留属于本基准的 `hundun`/`mpiexec` 进程。验收期间未启动
第二个高内存作业，未干扰其他研究进程。

**最终结论：`ACCEPT`。** 分支可用于后续由用户决定的本地融合或长时统计；
默认保留当前分支和工作树，不自动 push、发布或清理。
