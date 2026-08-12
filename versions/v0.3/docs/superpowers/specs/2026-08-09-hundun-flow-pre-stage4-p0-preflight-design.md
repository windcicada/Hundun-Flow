# HUNDUN-FLOW Pre-Stage-4 P0 并行预研设计

日期：2026-08-09

状态：用户已完成书面审阅并批准本设计，授权编写实施计划后连续执行 P0，不再设置
例行人工确认门。该授权不改变产品、Stage 3、私有来源、破坏性操作、push 或发布边界。

规划父节点：governance commit
`0a1c22b46d0e6dda218b75caed8e7f8080beccdd`。本规格补充而不取代
`2026-08-09-hundun-flow-stage4-6-linux-cpu-v1-architecture-design.md`。Stage 4 产品实现的
入口门保持不变：只有 Stage 3 正式接受且用户明确启动 Stage 4 后，才能执行 `4F-0`
及任何产品集成。本文只允许一个与活跃 Stage 3 隔离的非产品预研通道。

## 1. 目标

Stage 3 由另一 agent 推进期间，提前消除 Stage 4 最可能造成关键路径返工的外部风险：

- Cantera 及传递依赖的版本、来源、许可证、ABI 和构建可复现性；
- Ubuntu 22.04、GCC 11、libstdc++、C++17、ABI=1 下的共享库产物；
- 无 Python 用户环境中的 C++ 链接、加载、线程和小规模 MPI 生命周期；
- 相对 RPATH、移动安装前缀和完整 bundle 升级；
- 两套性质明显不同的代理燃料数据和机理的来源独立性；
- Stage 4--6 可复用的公开方程测试向量；
- Stage 3 接受后 `4F-0` 所需的只读 intake 清单。

预研的成功含义是“风险已经被独立实验定位并形成可复现候选输入”，不是 Stage 4 已实现、
已链接或已接受。预计它可以把 Stage 4 前端的包构建、ABI 和数据许可探索移出正式关键路径。

## 2. 不可突破的边界

### 2.1 禁止修改的范围

P0 不得修改、暂存、提交或清理：

- `/home/wyf/code_dev/hundun-flow` 产品仓库；
- 任一 Stage 3 worktree，包括 framework、infrastructure 和旧 Task 11 证据树；
- Stage 3 的产品源码、测试、CMake、schema、Checkpoint、diagnostics、driver 或 registry；
- 私有 COAST、COAST-2、BOFFIN、研究算例、研究数据和研究进程；
- 宿主机网络、firewall、routing、NetworkManager 或 Codex 服务。

P0 不得以活跃 Stage 3 worktree 的 dirty/untracked 内容推断未来接口。最终接口只来自正式
接受的 Stage 3 HEAD、tree、receipt 和安装结果。

### 2.2 禁止提前实现的能力

P0 不创建 HUNDUN 产品头文件、`src/` 实现、中央 CMake target、schema v4、Checkpoint v4、
diagnostics provider、ChemistryBackend、flow coupling、driver 或 capability claim。它也不运行
反应流、ESF/TCR、喷雾或 48^3 长算例。

### 2.3 版权边界

- Cantera 和每个传递依赖保持原 upstream identity、许可证、版权和免责声明；
- 机理文件和液体数据是独立版权对象，不能借用 Cantera 许可证；
- 没有明确再分发许可的数据只登记为 user-supplied candidate，不进入安装包；
- 不从 Cantera 零散复制函数到 HUNDUN 源码；
- 本阶段不访问 COAST。ESF/TCR 和 `EXEC/Fuels` 的私有 oracle gate 仍留在 Stage 5/6，
  并要求用户届时重新确认 exact realpath 和版本。

## 3. 隔离拓扑

### 3.1 Git

书面规格和后续预研计划位于 governance 独立分支。用户批准书面规格后，从该规划提交创建：

```text
branch:   coast/stage4-p0-preflight
worktree: /home/wyf/code_dev/.worktrees/hundun-flow-stage4-p0-preflight
```

该 worktree 只允许治理文档、预研 manifest、独立测试协议和 receipt。P0 worker 不提交、
不添加 DCO；主 agent 完成完整 diff、来源和许可证审查后创建签署提交。

### 3.2 外部生成目录

下载、解压、构建、安装、临时源码和二进制不进入 Git，统一放在：

```text
/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/
  inputs/
  source/
  build/
  install/
  spikes/
  logs/
  manifests/
```

每个子目录用内容哈希和 profile 标识，不覆盖旧产物。任何删除、替换或清理另行核对精确
目标；默认保留全部候选产物和失败日志。正式 Stage 4 不直接信任目录名，只信任 manifest
中的文件哈希。

### 3.3 资源

P0 只使用低/中负载资源：

- 任意时刻最多一个 third-party 构建；
- 默认 `-j16`，确认不与 HUNDUN 高负载作业竞争时最多 `-j32`；
- 独立 C++ spike 最多 2 MPI ranks，不运行数值网格；
- 构建进程使用较低调度优先级，并保留完整日志；
- 启动前只检查属于 HUNDUN 的构建/MPI 作业，不检查或干扰研究进程。

## 4. 预研工作包

### P0-0：冻结预研身份与证据格式

记录规划 HEAD/tree、worktree 状态、宿主 profile、编译器、MPI、CMake、CPU ISA policy、
外部根目录和允许文件。定义统一 manifest、日志命名和结果枚举：

```text
PREFLIGHT_PASS
PREFLIGHT_PARTIAL
PREFLIGHT_REJECT
```

任何结果都必须带 `stage4_product_accepted=false` 和 `product_changes=none`。

### P0-1：Cantera provenance 与依赖 BOM

从官方公开来源核验 Stage 4 规格冻结的 Cantera release/tag/commit/source SHA。逐项登记实际
构建所消费的 SUNDIALS、yaml-cpp、fmt、Eigen、BLAS/LAPACK 等依赖：

- 官方 URL、版本、source archive SHA-256；
- copyright、SPDX/许可证全文和再分发条件；
- 本地 patch 及独立 diff；
- shared/static、编译器、标准库、exceptions、RTTI、ABI、ISA；
- 是否进入最终 bundle，以及不进入的理由。

候选版本在 archive 和许可证核验前不得写成 accepted。P0-1 完成前不得开始二进制构建。

### P0-2：外部 Release artifact prototype

在冻结的 Ubuntu 22.04/GCC 11 profile 中，以本地 hash-verified 输入构建 Cantera Release
共享库和必要传递依赖。上游若需要 Python/SCons，只允许存在于 maintainer builder 环境；
普通 HUNDUN 用户侧不消费该工具链。

产物 manifest 至少记录：

- 所有输入 archive、patch、工具和环境镜像的 identity；
- 完整 configure/build/install 命令和退出状态；
- C++ 标准、`_GLIBCXX_USE_CXX11_ABI`、exceptions、RTTI 和 ISA flags；
- 每个 header/data/shared library 的 SHA-256；
- `readelf`、`ldd`、RPATH/RUNPATH 和未解析符号结果；
- 日志 SHA-256、wall time、峰值 RSS 和 package size。

“network-independent consumer”只通过本地输入、命令 trace 和 fetch-string policy 验证，
绝不切断宿主机网络。

### P0-3：独立 C++、线程、MPI 和 relocation spike

在外部 `spikes/` 编写最小、非产品 C++ 程序，通过 Cantera 公共 C++ API 完成：

- 加载一个具有独立许可与 hash 的小机制；
- 计算 thermo、mixture-averaged transport 和一个短 0D chemistry interval；
- 每 rank 建立独立 runtime，每线程建立完整独立 mutable workspace；
- 1/2-rank 输出固定顺序的数值摘要和对象身份检查；
- 普通 GCC Release 与 ABI 兼容的 HUNDUN-style Debug consumer 链接；
- 移动 install prefix 后重新运行；
- 检查没有 Python、Conda、绝对 build-tree RPATH、ABI=0、`_GLIBCXX_DEBUG` 或
  Clang/libc++ 混链。

spike 不链接 HUNDUN，不证明 Stage 4 flow coupling，只验证 third-party 边界。

### P0-4：公开数据、测试向量与许可候选

建立两条互不依赖的代理燃料候选链：n-dodecane 类航空煤油代理和 iso-octane 类汽油代理。
每条链分别登记气相机理、热化学数据、输运数据和纯液体性质来源。低成本 C++ spike 只证明
解析、species identity、基本物性和非硬编码接口，不宣称真实航空煤油或汽油验证。

同时依据公开论文和标准方程生成候选测试向量：

- Stage 4：理想气体混合物、总热化学焓、0D/PSR、Strang 时间对称性和失败原子性；
- Stage 5：Philox published vectors、N=2/4 antithetic Wiener、IEM 精确指数混合、
  simplex/element/mean 一致性；
- Stage 6：parcel/gas 动量与总热化学能符号、D^2-law、drag、heat/mass transfer 和
  source deposition 守恒。

Stage 5 TCR 数值向量只依据两篇已登记论文形成独立推导候选；在 COAST gate 前不得称为
COAST-equivalent 或 validated。

### P0-5：Stage 3 到 Stage 4 intake dry-run

只编写检查清单和只读命令，不冻结当前 dirty Stage 3 接口。清单覆盖：

- accepted code/product HEAD、tree、version 和 governance receipt；
- installed public headers、targets、symbols、CMake exports 和 driver commands；
- schema v1--v3、Checkpoint v1--v3、diagnostic/provider ID；
- field/final-flux/transaction/rollback/MPI authority；
- Stage 3 后台作业与工作树状态；
- `4F-0` 必须重新确认的 planned paths 和 collision list。

真正的 intake 只能在 Stage 3 接受后针对 exact accepted HEAD 执行。

## 5. 顺序和并行

默认执行图为：

```text
P0-0 -> P0-1 -> P0-2 -> P0-3
             \-> P0-4
P0-0 ---------------------> P0-5
```

P0-2 和公开来源的 P0-4 在 P0-1 冻结输入后可以并行：前者是单个低优先级构建，后者是
只读来源/数学审查。P0-3 等待 artifact。P0-5 可与两者并行，但主 agent 保持 central
ownership。不得同时运行两个 third-party 构建或任何大型 MPI/数值任务。

边界清晰的 archive/hash/license 收集、独立 build log 监看和机械 manifest 校验可委派给
默认 worker。总体依赖、ABI 判断、数学推导、许可证结论、完整 diff 和是否接受 artifact
由主 agent 负责。

## 6. 证据复用与失效

| P0 证据 | Stage 4 可直接复用的条件 | 必须重做的部分 |
|---|---|---|
| source/license/BOM | URL、revision、archive、patch、license 未变 | 任一输入或许可变化即重审 |
| Release artifact | 输入、builder profile、compiler、ABI、ISA、命令和二进制 hash 全同 | artifact producer 的正式治理接入 |
| C++/thread/MPI spike | 同一 artifact、机制、workspace policy 和 executable hash | HUNDUN adapter、CMake target 和 flow integration |
| relocation evidence | 完整 bundle 和 RPATH manifest 全同 | 正式 HUNDUN install/CPack smoke |
| public math vectors | 方程、单位、ordering 和 contract 未变 | 对应 mutation-sensitive product RED |
| intake checklist | 命令模板可复用 | 必须针对 accepted Stage 3 HEAD 重新执行 |

治理文档变化本身不使二进制证据失效；产品 tree 或接口变化不能把 standalone spike 变成
产品证据。P0 的所有 receipt 都必须列出“可复用”和“未证明”两栏。

## 7. Stage 4 合流门

Stage 3 正式接受且用户明确启动 Stage 4 后：

1. 从 accepted Stage 3 governance/code head 创建正式 Stage 4 worktree；
2. 执行 `4F-0`，不得从 P0 猜测 accepted Stage 3 identity；
3. 对 P0 provenance、artifact、机制和测试向量逐哈希复核；
4. 通过正式 mutation RED 实现 `4P-1..4P-4`；
5. 只有正式 HUNDUN CMake、安装和 C++ adapter gate 通过后，P0 artifact 才能升级为
   Stage 4 accepted input；
6. P0 分支不作为 accepted Stage 3 的祖先，不把实验源码或构建树 cherry-pick 进产品；
   只移植经审查的治理记录和后来按正式文件边界实现的内容。

若 P0 证明方案 A 无法满足无 Python 用户环境、ABI 或可重定位包要求，立即以
`PREFLIGHT_REJECT` 停止，不擅自切换到方案 B/C。主 agent提交设计修订并取得用户批准后
才能改变打包路线。

## 8. 完成标准

P0 完成必须同时满足：

- Stage 3 和 product worktree 未被修改、清理或干扰；
- Cantera/依赖 provenance 和许可证不存在未解释缺口；
- 至少一个 frozen Linux CPU artifact candidate 可由本地输入重建；
- standalone C++、线程、1/2-rank、RPATH 和 relocation smoke 有 hash-bound 证据；
- normal-consumer 环境不需要 Python/Conda/Cantera 系统安装；
- 双代理燃料候选明确区分“可再分发”和“必须由用户提供”；
- 公开数学向量只声明候选 oracle，不越权声明产品科学接受；
- intake 清单已准备，但没有把未接受 Stage 3 当作基线；
- governance worktree clean，提交带合法 DCO；
- 没有遗留 P0 构建、MPI 或后台进程。

完成 P0 后停在边界。若 Stage 3 尚未接受，只报告预研结果并等待；不得自动进入 `4F-0`
或 Stage 4 产品实现。
