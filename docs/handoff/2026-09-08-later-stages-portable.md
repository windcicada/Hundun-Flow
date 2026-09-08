# HUNDUN-FLOW 可移植模块实施包台账

交付入口：`codex/stage4-6-coast-replacement`。本任务只完成移植前模块，
不修改主线，不自动合并或推送，不宣称产品燃烧/喷雾验收或完整替代 COAST。
本文件取代原分 Stage 的后续开发台账；旧记录仍可从基线 Git 历史读取。

## 基线、归属与提交

| 项目 | 记录 |
| --- | --- |
| 自有 worktree | `/home/wyf/code_dev/.worktrees/hundun-flow-stage6-two-phase` |
| 开工 HEAD | `3d281086182865493a2daca63f4686e79147acf7` |
| tree / parent | `8006da6602526ca976c28d440e71b78e20bc609c` / `e8939c863d2deea1c1047a5742fc9cd358a3116f` |
| 开工 dirty / untracked | 均为空；没有覆盖既有用户修改 |
| 产品共同祖先 | `86542bb96678ec844ae5ac95d8f6391993da239e` |
| 主线观察 | 统一分支建立时为 `codex/v04-restart-receipt-observability@461d7631aa251969ddfc4e639bcadb74c640070c`；只是历史快照，不作为新的接入 head |
| 独立队列 | A `/tmp/hundun-portable-queue-a`：P1→P4→气膜桥；B `/tmp/hundun-portable-queue-b`：P2→P3及验证输入；C `/tmp/hundun-portable-queue-c`：P6→P7→P5。均从同一冻结提交分支，只显式归集局部文件 |
| DCO 节点 | `5149e92` 接口冻结/物性拆分；`bbd2505` TAB/单元汇总；`b15aefa` 生命周期/定位快照。后续归集节点见下方验证记录及本文件 Git 历史 |

协调者独占公开头、中央 CMake、模块间合同和本文件；worker 不提交主线，
不合并 governance 历史、旧 driver/schema/checkpoint，也不碰 Halo 优化。

## 共同合同

- TCI closure：finite-rate mean、PaSR、ESF/TCR；chemistry representation：
  direct Cantera 与解析验收 provider。FGM 只保留表示边界，无表模型实现。
- `GasIdentity` 逐字段绑定机理 SHA、phase、组分顺序、元素矩阵、分子量、焓参考及
  原有两个指纹；不改变旧指纹含义。查询与推进 provider 均声明完整身份；
  backend/query 配对及跨 rank 身份不一致明确拒绝。
- 普通值/借用视图/provider 输入，candidate/report 输出。完整 Ns 组成支持
  `(p,h,Y)`、`(p,T,Y)`；parcel 不直接调用 Cantera。材料到 vapor species 的映射
  同时绑定液体内容身份与气相身份，不只传裸索引。
- 交换采用区间积分总量：kg、kg·m/s、J，新状态减旧状态；体积量/瞬时率另作显式转换。
  multiplicity 已含在 parcel 交换中，沉积权重只用一次；同一气源作用于每个 ESF，不乘除 N。
- 热路径 HUNDUN 缓冲预分配、容量不足拒绝。旧 vector backend 桥和普通值快照属于冷准备，
  可分配；Cantera/CVODE/solverStats、MPI 内部资源行为不宣称零分配。
- revision、accepted step、算法版本、区间及 RNG 地址固定重试。任何模块失败，
  整组候选不可用；借用候选在下一次 prepare（即使失败）后失效。无隐式换模型、
  裁剪归一化、补造时间尺度或部分发布。
- 快照只含已接受普通值：fields、TCR 历史、parcel/TAB/ordinal、injector 余量/ordinal、
  ESF/parcel RNG 和 clocks。完整验证后生成 owning 恢复候选；无产品磁盘格式或 Restart section。

## 实施包与移植记录

下列代码已归集并通过模块 focused 检查；V1/V2 尚待执行。所有 P→I 的产品适配
均将失败映射为本次尝试整体拒绝，而非部分 source admission。

| 包 / 依赖 | 移植前实现 | 移植后消费接口、持久状态与验收 |
| --- | --- | --- |
| P1 / 冻结接口 | 独立 Cantera 与解析 backend；representation/rate-query；净质量生成率、组分焓、thermo/transport；有界查询/推进批量；完整身份、异常与 pool 生命周期；finite-rate/PaSR/ensemble | 主线物性及 Ns−1 视图、单位体积率/积分转换、机理身份持久化、化学源准入与反应流时序 |
| P2 / P1 | Philox/Wiener counter；N=2/4 持续候选；外部确定性通量/梯度、零湍流随机退化、通量修正、exact IEM、均值/元素/方差；真实解析/Cantera 桥连接 | 最终面通量 kg/(m²·s)、梯度、边界和 accepted field storage；RNG 与 cell/revision；真实空间输运验收 |
| P3 / P2 | 具名 TCR 根/统计映射、accepted/trial 历史、折点证据、off/shadow 无反馈；实验反馈未解拒绝；validated 缺科学证据拒绝 | cell/revision、跨步分支/折点历史、真实统计反馈与持久化；独立科学证据的产品准入 |
| P4 / P1 | 严格 SI 资产加载/范围/来源/FNV 内容身份；绝对液相焓积分、潜热与气相焓一致性；两套合成包、不同相关式和组分顺序；one-third full-Y film 与中立采样桥 | 主线气相 p/h/Y/速度采样；冻结真实燃料资产及材料映射；真实燃料独立验收 |
| P5 / P4/P6 | 固定一次 predictor/corrector 外层；整步/两半步误差及有界细分；碰撞后重采样；跨 cell、完全蒸发、出口、壁面、破碎统一事件；实际 Δm/ΔP/ΔH/ΔK；动态 TAB 与子滴剩余时间 | 几何/IBM/边界查询、真实气相采样与沉积；accepted parcel clocks；真实边界与步长收敛 |
| P6 / 冻结接口 | 公开 TAB 方程、代表等径子 parcel、成对横向速度；表面/形变/体动能分账；完整 parent ID+step+ordinal；保留独立给定粒径模型 | 共同事务父删子建、所有权；TAB/ordinal/ID 持久化；真实破碎预算 |
| P7 / P5/P6 | 坐标→cell/owner 双端预检、容量与全局 ID；18-lane v2 普通值候选迁移；accepted lifecycle snapshot 与恢复 | 真实分区/MPI 调度、共同提交、产品 Restart；injector/TAB/RNG/clock 全状态；连续/恢复一致性 |
| P8 / P1–P7 | cell/ID/segment 确定性汇总；先求实际交换和，再算单元非线性 ΔK；完整身份/区间/global ID 审计（含消失父滴/全部子滴）；候选源路由、共同接受/拒绝和普通值快照；壁面/出口/TAB 能量单列 | `rho*h-p`/压力功、源准入、两次 PISO 与执行图；真实 field/Halo/stencil；所有产品参与者一致提交与回退 |

P8 是规定气相状态上的普通数组组合器，顺序明确为
source → transport/IEM/TCR → 两个连续 chemistry 半区间；**不是 Strang 分裂或第二套流动 driver**。
cell inventory density 与 EOS 查询 density 分开报告，不能据此宣称已实现联立压力/密度闭合。
源路由是最多 4 ranks、固定全局容量的参考 Allgatherv 普通值通道，每个 cell 唯一 owner；
不复制全局气相库存，不是 field Halo 或可扩展生产通信实现。当前 P8 每轨迹段沉积到一个
cell；真实跨分区加权 stencil 的归一性/所有权留给产品适配，不静默接受缺失权重。

PaSR 保留 `finite_rate_mean_shadow` 基线，`kappa=0/1` 精确退化；
`tau_mix=C_Z*Delta²/[2(D+nu_t/Sc_t)]`，显式报告 C_Z，不构造 RANS epsilon。
κ 同时缩放物种反应增量和形成焓放热报告；总热化学焓不重复加放热。
TPDF 每场两半区间：PaSR 2 次、N=4 ensemble 8 次 chemistry 调用。

## 来源与科学边界

| 来源 | 固定身份 / 用途 |
| --- | --- |
| HUNDUN governance 只读供体 | Stage 4 `6407cd7c591ce088db7f1dd7e296d77acd18da1c`；ESF/TCR `8ffdf2b5673374fb14639fc4dce09a1b586ee5db`。只移植纯算法/测试行为，无历史合并 |
| TAB 原始方程 | [O'Rourke–Amsden 原始报告](https://digital.library.unt.edu/ark:/67531/metadc1106123/m2/1/high_res_d/6118786.pdf)，PDF p7 式14、19–22，K=10/3、Cb=1/2、Ck=8、y=1；代表粒径/等径子 parcel 是本实现明确约定 |
| A–S | [Abramzon–Sirignano 1989](https://doi.org/10.1016/0017-9310(89)90043-4)；解析低 Re/恒定膜极限与独立参考，不以 OpenFOAM 输出作真值 |
| TCR | 作者理论文档 `Coast_software/docs/TCR_Model/Thesis_Part2.md` 式2.80、2.88、2.90、2.91；只读方程，不复制 COAST 实现。统计映射版本 `ideal_gas_reactant_mole_fraction_v1`；不继承旧 oracle 指纹作为 validated 证据 |
| OpenFOAM 仅问题组织参考 | 固定 `43eea1d4b6ae2fdf67a638cf0f452bc2ece123c4` 的 [chemistry](https://github.com/OpenFOAM/OpenFOAM-dev/tree/43eea1d4b6ae2fdf67a638cf0f452bc2ece123c4/test/chemistry) / [Lagrangian](https://github.com/OpenFOAM/OpenFOAM-dev/tree/43eea1d4b6ae2fdf67a638cf0f452bc2ece123c4/test/Lagrangian)，未复制实现、字典或脚本，未运行上游测试 |

合成 alpha：k=2/s、cp=1000 J/(kg·K)、A/B 顺序、常 cp/Antoine；
beta：k=3/s、cp=1200、B/A 顺序、多项式 cp/Clausius。资产位于
`versions/v0.4/tests/fixtures/synthetic-liquid-{alpha,beta}.asset`，
原始字节 FNV1a64 分别为 `6004043157121730787`、`668675689539421851`。
机理 SHA 是合成方程规范的 SHA，不是假称真实燃料机理。真实 Cantera conformance 使用
独立合成机理 `synthetic-mechanism-v2.yaml`，SHA
`c518a07cada5f1bddcdb308f0a2f695d92cc6373e173ffd87e96312530b52aee`。

## 验证与执行入口

RED 先于实现：新增 seam 缺声明/符号；失配 backend/query、跨 rank 区间/机理、
损坏 RNG snapshot、不同 next revision、无反应时负 EOS density、伪造 ΔK、
陈旧采样位置及丢失 TAB 预算均出现预期失败，再实现 GREEN。解析参考在
`tests/validation/portable_v1_manifest.json` / `portable_v2_manifest.json`
中预登记公式、常数和绝对+相对容差；不从被测实现生成参考。

当前增量 focused 32 项及真实 backend selector 已通过；最终精确 HEAD 重验记录待补。
V1/V2 只在 P1–P8 完成后执行；空间夹具 4³、最多4 ranks，每配置一轮，
单项5/10/20秒超时，超时记失败，不扩容重跑。无 full ctest、既有流动回归、
正式火焰/喷雾长算、COAST 比较或 OpenFOAM Allrun。

```sh
cmake -S versions/v0.4 -B /tmp/hundun-flow-stage6-two-phase-build-libcxx \
  -DHUNDUN_BUILD_TESTS=ON -DHUNDUN_BUILD_PORTABLE_CANTERA=OFF -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_COMPILER=/home/wyf/.local/bin/clang \
  -DCMAKE_CXX_COMPILER=/home/wyf/.local/bin/clang++ \
  -DCMAKE_CXX_FLAGS=-stdlib=libc++ -DCMAKE_EXE_LINKER_FLAGS=-stdlib=libc++
# 只构建对应 v04_models_*_test / v04_public_headers_test 与 v04_portable_v1/v2 目标。
LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu \
ctest --test-dir /tmp/hundun-flow-stage6-two-phase-build-libcxx --output-on-failure \
  -R '^v04_(models_(chemistry(_adapter)?|combustion(_allocation|_common_source)?|spray(_breakup|_parcel|_mechanics|_source|_transfer|_properties|_film_bridge|_events|_migration_mpi_[124]|_transaction_mpi_[124])?|esf(_reaction|_backend)?|tcr|exchange_batch|exchange_routing_mpi_[124]|portable_composition_mpi_[124])|public_headers)$'
# 完成门后，分别执行：
# ctest ... -L '^portable_v1$' --output-on-failure
# ctest ... -L '^portable_v2$' --output-on-failure
git diff --check
```

Cantera 使用既有只读 Jammy/GCC11 兼容 rootfs 与已验证3.2.0依赖包，不下载/重建依赖。
`HUNDUN_BUILD_PORTABLE_CANTERA=ON` 单独构建纯模型静态库，不链接产品 driver/core。
本机源映射到 `/mnt`、依赖包到 `/tmp/cantera`、合成机理到 `/tmp/mechanism`，
构建目录 `/tmp/hundun-unified-cantera.XGRtBh` 映射到 `/tmp/build`：

```sh
cmake -S /mnt/versions/v0.4 -B /tmp/build -DCMAKE_MAKE_PROGRAM=/tmp/make \
  -DCMAKE_C_COMPILER=/usr/bin/gcc-11 -DCMAKE_CXX_COMPILER=/usr/bin/g++-11 \
  -DHUNDUN_BUILD_TESTS=ON -DHUNDUN_BUILD_PORTABLE_CANTERA=ON \
  -DHUNDUN_CANTERA_PACKAGE_ROOT=/tmp/cantera \
  -DHUNDUN_PORTABLE_CANTERA_TEST_MECHANISM=/tmp/mechanism/synthetic-mechanism-v2.yaml \
  -DCMAKE_BUILD_TYPE=Debug
cmake --build /tmp/build --target v04_models_cantera_backend_test -j2
ctest --test-dir /tmp/build -R '^v04_models_cantera_backend$' --output-on-failure
```

## 只登记、不实施的接入项

主线 source admission/ContributionRegistry、最终面质量通量/场/Halo/边界、
化学—输运—喷雾—压力/焓/密度顺序、全部参与者共同提交与回退、
schema/Checkpoint/Restart/diagnostics/CLI，以及恢复两相接入范围后的真实耦合算例和
COAST 替代验收。FGM、稠密/多组分液滴、液膜、碰并、移动 IBM、AMR、Halo 优化不在本次范围。
