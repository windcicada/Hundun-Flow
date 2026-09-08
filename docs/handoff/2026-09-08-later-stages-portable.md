# Stage 4–6 统一分支与主线移植交接

统一入口为 `codex/stage4-6-coast-replacement`，后续 Stage 4/5/6 开发与主线移植
均从该分支推进，不再分别合并旧 Stage 分支。目标为约定业务范围内替代 COAST；
本次完成代码归集，不宣称完整 Stage 接受或 COAST 替代接受。
主线工作树只读，不修改其 driver、准入、字段、schema、Restart 或验收台账。

## 身份与边界

| 对象 | 本轮开始时记录 |
| --- | --- |
| 本任务 | 复用自己的 `/home/wyf/code_dev/.worktrees/hundun-flow-stage6-two-phase`；由 `codex/stage6-two-phase` 新建 `codex/stage4-6-coast-replacement`，旧 ref 保留 |
| HEAD / tree / parent | `216ac889bc3f600cdcb7af9467e2dfb5beea37de` / `d1e926f84f1f19d5e698b3b860c1901634c57554` / `86542bb96678ec844ae5ac95d8f6391993da239e` |
| dirty / untracked | 统一前两份 CMake 与 public_headers 测试已修改；combustion、parcel、mechanics、transfer、source、breakup、migration、transaction 新文件及本交接未提交；全部保留并收拢 |
| 主线只读快照 | 建立统一分支时 `codex/v04-restart-receipt-observability@461d7631aa251969ddfc4e639bcadb74c640070c`；tree `e1e46205ca609ad32e35875bb30097c1b787b008`；parent `2ea65b60e8b45264232037e0bf2430d84e884186`。主线仍在独立推进，不将此快照当最终接入 head |
| 主线范围 | 2026-09-08 台账仍将两相**产品接入**暂挂；本任务按用户新指示继续独立模块开发，不改变主线范围 |
| README | 独立文档工作树的改动暂停，未混入本分支，未提交或推送 |

V04-2 是旧任务阶段名；本轮以实际主线快照和当前接口为依据，不沿用旧 accepted 标签。

## 兼容性清单

| 分类 | 供体/行为 | 处理 |
| --- | --- | --- |
| Stage 4 显式移植 | governance `6407cd7c591ce088db7f1dd7e296d77acd18da1c` 的 chem_composition、chem_reports、backend/service/workspace 与 chem_cantera_backend | 合并纯声明到 `v04_chemistry.hpp`；`v04_cantera.hpp` + `models_cantera.cpp` 使用独立配置值，不引入 cfg_resolved_case_v4 或旧 schema |
| 直接移植 | Stage 5 finite-rate/PaSR、命名时间尺度与 candidate/report；Stage 6 值类型、SN/RM、简化蒸发、TAB 振子 | 从 `codex/stage5-stage6-portable@a20b85b` 显式取文件；不合并治理历史 |
| 模块内重写 | SoA 容器、injector、真实 A–S、共享 stencil、轨迹、IBM 事件、候选迁移 | 新建 `models_spray_*` 局部模块；仅值输入、候选输出和聚焦测试 |
| accepted-head 适配 | Stage 4 backend 到 Stage 5 representation/rate-query 的映射、ESF 字段与 TCR mapper、气膜物性查询 | 独立 Cantera backend 已归集，但尚未连接 Stage 5/6 或主线；接入前核对 thermo/transport/源项区间 |
| 产品重写 | driver、ContributionRegistry、压力/焓/密度顺序、field/halo、common-source、schema/Restart/diagnostics | 保留给主线接受后的适配任务；局部 finish 测试不代表这些面已接入 |

旧供体提交：Stage 5 `25287ece`；Stage 6 `a3cdd925`。本分支已有 DCO 提交
`216ac88 feat(v0.4): add portable spray kernels`。供体简化蒸发为 Spalding/d²，
与真实 A–S 分开命名；其旧 `cp*T` 能量报告不能直接用作产品总热化学焓源。
SplitMix parcel RNG 与 Stage 5 历史供体中的 Philox 身份分别记录，不互相冒充。

统一代码 DCO 提交：`e8939c8 feat(v0.4): consolidate Stage 4-6 portable development modules`。
本文件的后续 DCO 提交只封存交接记录，不改变已测试的模型源码。

Stage 4 seal `033a685c`、Stage 5 framework `41b2aac9` / seal `02b57cce`、
field follow-up `coast/stage5-field-validation@8ffdf2b` 是只读供体身份，
不是本分支祖先或新验收。历史 ESF 输运/IEM/Philox/TCR 与 field executor 尚未整体移植；
它们和旧反应流 driver 的重写继续在本统一分支完成，不能将归集当成功能齐备。

## 当前实施包台账（2026-09-08）

当前任务基线：`codex/stage4-6-coast-replacement@3d281086182865493a2daca63f4686e79147acf7`；
tree `8006da6602526ca976c28d440e71b78e20bc609c`，parent `e8939c863d2deea1c1047a5742fc9cd358a3116f`。
开工时 tracked/untracked 均为空。主线与 pressure-enthalpy 工作树只读，不自动合并或推送。

共同接口冻结：`v04_portable.hpp` 的 Revision、GasIdentity、GasQueryProvider、
ExchangeDelta；已有 composition/closure 指纹含义不变。数组均为完整 Ns；显式
SHA/phase/组分顺序/元素矩阵/焓参考映射在准备阶段校验。查询支持 p,h,Y 与 p,T,Y，
parcel 不接触 Cantera；HUNDUN 容量预分配，第三方内部资源行为单列。
候选失败不可用；accepted step/revision/算法版本确定重试；快照只含已接受普通值。
所有交换为区间积分总量、新减旧，multiplicity 与沉积权重各一次。

协调者独占公开头、中央 CMake、共同源合同、本台账和归集提交。
worker 独占其局部 cpp/detail.hpp 与 focused tests，在相同冻结提交的独立临时
worktree 工作；不得写产品主线/旧供体，不提交、不改中央文件。接口变更先由协调者处理。
物性/气膜已从 transfer 拆入 `models_spray_properties.cpp`，原 transfer selector 1/1 PASS。

| 包 | 依赖 | 状态 | 实现与验收出口 | 移植后消费接口 / 持久状态 / 产品验收 |
| --- | --- | --- | --- | --- |
| P1 化学 | 冻结接口 | 待实施 | direct+解析 provider、显式映射、区间/净率/焓、批量容量和失败 | Ns−1/物性、源准入；机制身份；真实反应时序 |
| P2 ESF | P1 真实组合 | 待实施 | Philox/Wiener、N=2/4 输运、exact IEM、通量修正 | 通量/梯度/边界；accepted fields/RNG；真实空间输运 |
| P3 TCR | P2 | 待实施 | 统计映射、根可接受性/历史、折点/拒绝、证据准入 | cell/revision；branch history；真实统计反馈 |
| P4 物性/气膜 | P1 | 待实施 | 焓积分、包读取/身份、两套合成包、完整气膜 | 气相采样；资产身份；用户真实燃料验收 |
| P5 轨迹/事件 | P4/P6 | 待实施 | 自适应蒸发、重采样、统一事件/出口 H+K | 几何/采样/沉积；parcel accepted clocks；真实边界 |
| P6 TAB | 冻结接口 | 实现/本地验证完成 | 粒径预测、成对速度扰动、表面/变形/体动能、子候选与 ID | 共同父子替换；TAB/ordinal；产品破碎预算 |
| P7 生命周期 | P5/P6 组合 | 独立实现通过，事件组合待 P5 | 定位预检、容量、普通值快照、1/2/4 rank | 分区/共同提交/Restart；injector/parcel/RNG；恢复一致性 |
| P8 组合 | P1–P7 | 汇总核通过，组合待依赖 | cell/ID/segment 有序汇总、非线性 K、公共源、共同候选 | rho*h-p/压力功/PISO/执行图；全参与者状态；联立验收 |

接口冻结 DCO：`5149e92`。P6 来源：
[O'Rourke–Amsden 原始报告](https://digital.library.unt.edu/ark:/67531/metadc1106123/m2/1/high_res_d/6118786.pdf)
PDF p7，式14、19–22；独立实现 K=10/3、Cb=1/2、Ck=8、y=1 的代表等径模型，
偶数2..16子 parcel 成对横向速度闭合预算；不兼容阈值明确拒绝。
P6 RED 为新增函数链接缺失，GREEN 为 `^v04_models_spray_breakup$` 1/1 PASS。
P8 汇总核 RED 为新接口链接缺失；`^v04_models_exchange_batch$` 1/1 PASS：
两滴→单元动能解析预算、顺序/重试、重复项、revision、容量、权重及出口分账。
该核固定容量，无场写入，借用候选在下一次调用（即使失败）后失效。
P6/汇总核 DCO：`bbd2505`。
P7 新 seam 链接 RED 后实现；根代理复核
`^v04_models_spray_(parcel|migration_mpi_[124])$` **4/4 PASS**。
checked migration 在发送/接收方均验证坐标、cell、owner与跨rank revision；
旧 prepare 仅保留低层范围检查。快照不含 trial，恢复是 owning 候选、冷准备可分配，
实际容器/Restart 发布仍由后续产品事务承担。

OpenFOAM 测试问题参考固定为
`43eea1d4b6ae2fdf67a638cf0f452bc2ece123c4`（只读 ls-remote 核实）。
仅参考 [chemistry 测试组织](https://github.com/OpenFOAM/OpenFOAM-dev/tree/43eea1d4b6ae2fdf67a638cf0f452bc2ece123c4/test/chemistry)
和 [Lagrangian 问题分类](https://github.com/OpenFOAM/OpenFOAM-dev/tree/43eea1d4b6ae2fdf67a638cf0f452bc2ece123c4/test/Lagrangian)；
不复制其实现、字典或脚本，不运行上游测试。

所有 P→I 失败由产品适配映射为整次尝试拒绝；不得部分源发布。
V1/V2 仅在包全部完成后运行：合成小算例，最多 12³/4 ranks，每配置一轮，
显式超时与预登记绝对+相对容差。真实燃料、field/Halo、schema/Checkpoint/Restart/
diagnostics/CLI、产品共同提交与 COAST 替代验收均只登记，不实施。
以下是基线历史证据，不代表当前 P1–P8 已完成。

## 基线切片与退出条件（历史）

| 切片 | 当前工作 | 本地退出条件 / 后续适配 |
| --- | --- | --- |
| Stage 4 backend | 组分/元素身份、总焓值与 interval/service 接口；独立 Cantera runtime/lane/thermo/transport/integrate 实现 | 纯值测试与一个现有 backend conformance selector 已适配；真实 interval/物理场数值验收未重跑 |
| Stage 5 fast path | 已显式移植并验证；新增分配失败状态 | 有界 kappa、0/1 精确退化、质量/元素/能量、失败无部分结果、PaSR 2 与 N=4 TPDF 2N 调用 |
| S0 容器与注入 | SoA、TAB 续算态、稳定 ID、注入余量、trial/preflight/commit/rollback | 值快照与 retry 一致性、容量及热路径分配检查；产品持久化另行适配 |
| S0 力学 | stencil、受限轨迹、静态 IBM 反弹、出口收支 | 解析极限/有限性/权重与最早事件；field/halo 采样接线另行适配 |
| S0 候选迁移与共同完成 | 17-lane 无填充打包、全局 ID 审核、固定字段与 parcel 共同 accept/reject 测试 | 小型 1/2/4-rank；无 driver 调用、无 chemistry/inert_source 准入变更 |
| S1 传递与蒸发 | 液体性质包、one-third film、A–S、预测–校正积分；已验证 | 相变/传热与机械能分开报告；陈旧表面温度样本明确拒绝 |
| S1 气相源代数 | 新增 `closed_isobaric_exchange_v1` 候选核 | 在固定热力学压力、无外力/壁面/出口功的交换子问题中封闭总 H+K；热通量与实际液滴焓差超过显式预算则拒绝 |
| S2 common-source | 已完成 N=2/4 纯值映射 | 同一物理源用于每个随机场；均值/方差、物种与焓预算、任一场失败整组撤回；无 parcel 循环或生产字段写入 |
| S2 子滴候选 | 给定统一粒径的守恒拆分，稳定子 ID、质量/动量/体动能/热化学焓预算 | 不预测 TAB 子滴粒径；表面能与变形能缺额单列，不能视为完整能量闭合的 breakup 模型 |
| S2 产品组合 | 尚未形成产品组合 | 完整 TAB、持久化与真实 CLI 仍需完成；不继承历史 seal |

## 验证约定

按已确认的纯值/候选接口执行 RED→GREEN，不测试内部私有方法。
仅运行本次模型的 focused selectors、公共头与小型 MPI 合同，禁止 full ctest、
Stage 3/4 数值回归、正式火焰/喷雾计算、COAST 对比及性能长测。Halo P0 不在本任务测试范围。

Clang/libc++ 构建目录：`/tmp/hundun-flow-stage6-two-phase-build-libcxx`。
系统 GCC 7 缺少既有源码所需的 `<filesystem>`，不因此更改产品源码。

```sh
cmake -S versions/v0.4 -B /tmp/hundun-flow-stage6-two-phase-build-libcxx \
  -DHUNDUN_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_COMPILER=/home/wyf/.local/bin/clang \
  -DCMAKE_CXX_COMPILER=/home/wyf/.local/bin/clang++ \
  -DCMAKE_CXX_FLAGS=-stdlib=libc++ -DCMAKE_EXE_LINKER_FLAGS=-stdlib=libc++
git diff --check
```

本轮 RED 证据：移植测试先于头/实现时缺声明或符号；Stage 5 分配失败注入使旧
`noexcept` 路径以 exit 90 终止，修复后返回 `workspace_failure`；陈旧 A–S 气膜表面温度
测试在修复前输出明确 FAIL。S1 总 H+K 代数与 S2 common-source 新接口均先 RED 后 GREEN。
容器、力学和 transfer 子切片另有各自严格编译与 ASan/UBSan 快测；本轮根代理也独立运行了
common-source 与气相源代数的 ASan/UBSan，均 exit 0。

```sh
cmake --build /tmp/hundun-flow-stage6-two-phase-build-libcxx --target \
  v04_models_chemistry_test v04_models_spray_breakup_test \
  v04_models_combustion_test v04_models_combustion_allocation_test \
  v04_models_combustion_common_source_test v04_models_spray_test \
  v04_models_spray_parcel_test v04_models_spray_mechanics_test \
  v04_models_spray_source_test v04_models_spray_transfer_test \
  v04_models_spray_migration_test v04_models_spray_transaction_test \
  v04_public_headers_test -j4
LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu \
  ctest --test-dir /tmp/hundun-flow-stage6-two-phase-build-libcxx \
  -R '^v04_(models_chemistry|models_combustion(_allocation|_common_source)?|models_spray(_breakup|_parcel|_mechanics|_source|_transfer|_migration_mpi_[124]|_transaction_mpi_[124])?|public_headers)$' \
  --output-on-failure
```

统一组合上述 **17/17 PASS**。子滴测试第一次组合运行暴露 `pow(d,3)` 与连乘的
浮点逐位相等误断言；改为 1e-14 相对/1e-24 kg 绝对容差后通过，没有改变模型公式。

Cantera 目标 `hundun_v04_portable_cantera` 默认关闭，不链接产品 core，不扩展准入。
从供体保留 `HundunPortableCantera.cmake` 中的依赖包 hash、license、ABI 与符号链接检查，
以及 `third_party/cantera/` 元数据和 license；未复制二进制、机理或第三方源码。
在兼容 Ubuntu 22.04/GCC 11/libstdc++ 环境下：

```sh
cmake -S versions/v0.4 -B <backend-build> -DHUNDUN_BUILD_TESTS=ON \
  -DHUNDUN_BUILD_PORTABLE_CANTERA=ON \
  -DHUNDUN_CANTERA_PACKAGE_ROOT=<verified-package> \
  -DHUNDUN_PORTABLE_CANTERA_TEST_MECHANISM=<synthetic-mechanism-v2.yaml>
cmake --build <backend-build> --target v04_models_cantera_backend_test -j2
ctest --test-dir <backend-build> -R '^v04_models_cantera_backend$' --output-on-failure
```

本机宿主 glibc 2.31 不兼容固定包；复用 P0 的只读 Jammy rootfs，以 bwrap 隔离，
仅统一源码只读挂载和 `/tmp/hundun-unified-cantera.XGRtBh` 构建输出可写，
无网络、无安装、无重建依赖包。使用独立挂载的 CMake/make 与 rootfs gcc-11/g++-11。
最终日志：该目录 `Testing/Temporary/LastTest.log`，selector **1/1 PASS**。
验证仅覆盖 backend 身份、工作区隔离/容量/生命周期与无 schema 的 solver controls 拒绝，
不代表 interval/thermo/transport 数值重新接受。无 schema 控制校验先 RED
（非法 tolerance 被准入）再 GREEN（工厂在领取 lane 前拒绝且不消耗 lane）。

## 后续单分支接入顺序

主线仍只读；合并权由主线任务在 accepted head 上执行，不自动 merge/rebase/push。
共同产品祖先为 `86542bb96678ec844ae5ac95d8f6391993da239e`，
本分支是其线性后继，无 governance merge，旧分支无需再合入。

1. Stage 4：backend 身份/焓基准/源区间适配 → 化学准入 → 真实 mean reacting driver。
2. Stage 5：finite-rate/PaSR → 完整 ESF/IEM/TCR → 共同回退与持久化。
3. Stage 6：主线恢复两相接入范围后，再接 parcel/gas/ESF 公共源与 Restart。
4. 冻结 COAST 业务案例和物理/性能/I/O 验收；不得把本分支提交或单测作为替代接受。

合并前先核对 CMake 与公开头，当前主线会继续演进，不承诺未来无冲突。
本次没有改 flow/product driver、ContributionRegistry、field/halo、schema、
checkpoint/Restart、diagnostics 或 Halo P0，没有启动 full ctest/长算/COAST 比较。

## 移植时不可省略的检查

- Stage 5 是完整 Ns 值接口，主线是 Ns−1 视图；native-air 的参考平移与生成焓基准须显式转换。
  当前 vector 版本报告分配失败，但不是批量零热分配的产品工作区。
- A–S 使用 [PeleMP 公开方程](https://amrex-combustion.github.io/PeleMP/Equations.html)
  中的气膜/Stefan/Nu0/Sh0 形式，来源为 [Abramzon–Sirignano](https://doi.org/10.1016/0017-9310(89)90043-4)。
  Ranz–Marshall 独立核保留，二者基准传递关联式不冒充同一实现。
- A–S 常系数事件步与有限子步积分可能产生 `thermal_exchange_state_residual_j`；
  必须显式检验、细化或拒绝，不用 `gas_h=-parcel_h` 隐去误差。气相源核要求先按单元累计
  实际 Δm/ΔP/ΔH/ΔK，再做非线性气相动能修正；产品的 rho*h-p/压力功约定仍需适配。
- stencil 在首/末单元中心之外显式夹持；IBM 反弹是分段事件处理，碰撞后未重新查询加速度。
  出口能量需液体焓服务；目前出口 ledger 仅质量/动量。两者均不是完整联立轨迹 driver。
- 17-lane 迁移验证 global-cell 范围与当前 Cartesian owner，不代替位置到网格的几何定位。
  模型局部共同完成 helper 只验证已准备参与者的共同提交，不证明化学准入、源项物理或 PISO 时序。
- 容器快照包含 parcel/TAB 值，injector 余量与 ordinal 可查询；尚无产品持久化编码、机制包加载、
  field/halo 采样、driver/schema/Restart/diagnostics 接线及真实双 surrogate CLI 验收。
