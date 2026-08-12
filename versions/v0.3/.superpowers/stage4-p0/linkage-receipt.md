# Stage 4 P0 Cantera C++ 边界验收记录

- 记录时间：`2026-08-09T23:58:09+08:00`
- 结果：`PREFLIGHT_PASS`
- 结果范围：`standalone_cantera_cxx_thread_mpi_abi_and_relocation_boundary_only`
- Stage 4 产品已接受：`false`
- 产品修改：`none`
- 治理父 HEAD：`8923015c3e489d45d66fbff9cc44bc4967b553a9`
- P0-2 artifact manifest 原始 SHA-256：`efd2fbbc9f497b7b0f7212104497591464a5e0ae76d84527528e9586f8338296`
- 追加 P0-3 证据后的 artifact manifest SHA-256：`62962ba6217c52d0b12fa6b27c3b5da0e7d55c5443ea95222614f5f84706521b`
- 外部证据根：`/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0`

artifact manifest 原有顶层 `relocation_status=pending_p0_3` 和 `review.next_gate=P0-3` 是 P0-2
冻结快照，按 append-only 规则不回写；新增的 `p0_3_linkage.result=PREFLIGHT_PASS` 与
`p0_3_linkage.relocation.status=pass` 是后续权威证据。

本记录只接受冻结 Linux CPU profile 上的 Cantera C++ 消费边界。它证明独立 C++
consumer 能链接 P0-2 Release 共享库，在每 rank、每 OpenMP lane 使用互不别名的完整
mutable workspace，完成 thermo、mixture-averaged transport 和短时 0D chemistry，并能随
完整 bundle 移动。它不是 Stage 4 产品接受，也不证明 HUNDUN CMake 集成、
`ChemistryBackend`、反应输运、Checkpoint v4、diagnostics v4 或任何真实燃料科学结论。

## 工具链与独立性

MPI spike 使用经过审计的 Ubuntu 22.04 rootfs、GCC 11.4.0、libstdc++ 和 Open MPI 4.1.2。
OpenMPI wrapper 在隔离 rootfs 中显式选择 `OMPI_CXX=g++-11`；未修改 rootfs。MPI 工具链
审计脚本 SHA-256 为
`ceed60e9272141e019fe8c309b8eefdb31659c5206bdcdc0e1232abf253de0ad`，执行日志 SHA-256 为
`7d8a55b70b6af075614e96fc1f251538392e6ff0c484616e57b1b29ad5eca900`，退出 0。

standalone spike 不包含 HUNDUN 头文件或产品源码。每个 MPI rank 只保存 immutable 的机制
路径、phase 名和机制 SHA；每个 OpenMP lane 独立创建 `Solution`、thermo、kinetics、
transport、Reactor 和 ReactorNet。Reactor 显式消费本 lane 已独立拥有的 Solution graph，
不再做隐藏的第二次 clone。MPI 只在主线程调用，初始化要求 `MPI_THREAD_FUNNELED`。

合成机制 v2 SHA-256 为
`c518a07cada5f1bddcdb308f0a2f695d92cc6373e173ffd87e96312530b52aee`。它由 HUNDUN-FLOW
P0 独立编写，使用人工 A/B 状态、独立选取的常数 NASA7 系数和人工可逆反应，只用于接口
测试；不代表真实组分、燃料、点火或化学机理。其 provenance SHA-256 为
`59b6d08311e82115b1f3d125ba745e10056f74e8ad9f12d7274c9a0c79fdecaa`。

## 构建与数值摘要

| 项目 | 证据 |
|---|---|
| Spike source SHA-256 | `d1f7890d6846643d2adc6aebf104b811b541225ed023b239ce33f68eeaadde10` |
| 可复现构建脚本 SHA-256 | `0e5779c2e8a48a19eea66c8ddd744aa7aea109463c9d8805d2d6bd2cbc5ac466` |
| Release binary SHA-256 | `9083224db30c2f750e9365fcae20e1bace39a2bd92ed09c7d2077d9061d5385a` |
| 普通 Debug binary SHA-256 | `a5fadbb62929e4f400b3d935b915d5d8949b249cec6f4f8bfa4825769583bd0c` |
| 构建日志 / 状态 | `bd23e6ee0d2a66ab3c179fd7f944cab2ad30d6904d5413f3f916902438a4d830` / exit 0 |
| 构建 wall / wrapper RSS | 4.69 s / 3,300 KiB |

固定输出 schema 为 `hundun.stage4_p0.cantera_spike.v1`。候选摘要为：

```text
temperature_K=1.00000000000000000e+03
pressure_Pa=1.01325000000000000e+05
density_kg_m3=1.22840891457018497e-02
cp_mass_J_kg_K=2.10336107899709968e+04
viscosity_Pa_s=1.57224613101764864e-05
net_production_rate_hash=6c9b671f67e0611d
reactor_final_state_hash=d82f7aa6bc186c8e
workspace_alias_count=0
```

两个 16-hex hash 是按固定顺序对 binary64 bit pattern 做 FNV-1a-64 得到的确定性接口摘要，
不是密码学证据，也不作为真实机理科学判据。

## 线程、MPI、Debug 和 mutation

低成本矩阵覆盖：

```text
Release: 1 thread / 1 rank, repeat 2
Release: 2 threads / 1 rank, repeat 2
Release: 1 thread / 2 ranks, repeat 2
Release: 2 threads / 2 ranks, repeat 2
Debug:   1 thread / 1 rank
```

九次运行全部退出 0、stderr 为空；固定 lane 顺序 hash 重复一致，相同输入的 rank 摘要一致，
所有跨 lane workspace alias 计数为零。矩阵 runner SHA-256 为
`d58adcde9ed22884fbd8988616f29a709d1ada7b6e9ace4efc769392d1d28c5d`，总日志 SHA-256 为
`7c798626be22a3aa8dfac49d64dc1b035effee4f2bbb05cd28063650e6285492`，证据索引 SHA-256 为
`df6d6875bd35a8eaf0923e072e64666ee8c0521c083da713adc66679fce52241`。

普通 `-O0 -g3`、ABI=1 consumer 成功链接同一 Release Cantera artifact。编译命令中不存在
`_GLIBCXX_DEBUG`、Clang 或 libc++。ABI=0 mutation 在执行前被 policy 以 exit 65 拒绝；
policy SHA-256 为
`fc81e27f4cb785ca20af90bb0ae61d9aa9ae45f08459c2c51c133675f4a1e905`。

固定输出 validator v2 SHA-256 为
`457ae8171eb91ff0a4b5fab6f14c62594daa79103d3363b3804814e8a905da34`。
它对缺字段、字段换序、workspace alias、非有限压力、非法 hash 和错误 schema 六个 mutation
全部 RED，unexpected pass 为零。mutation summary SHA-256 为
`e462f8c3311867ec61d56cd7ead553655616eda335e2b3b2c28a5e1339147977`。

MPI 2-rank/2-thread 的 `execve` trace SHA-256 为
`8ade95d07cf07233c17929e6d221234207cbd11eca9bf479bb883fe278a195ce`。
执行链中没有 Python、Conda、libpython 或 site-packages。

## 完整 bundle 重定位

P0-2 artifact 先完整复制到：

```text
install/relocated/cantera-3.2.0-gcc11-release-v4-p0-3-v1
```

复制前后的 653 个 artifact regular file manifest 完全相同，SHA-256 均为
`5fd187220f23d11af133f43f4c1208f9744943bc793e89bc7543171363fa3004`；两个 symlink 也完全
相同。随后只在移动副本中加入 standalone consumer 和已声明来源的 synthetic fixture。
最终副本包含 656 个 regular files、2 个 symlinks，`du -sb` 为 21,691,508 bytes。

移动后的 consumer SHA-256 为
`0b409c6124a435250e5851f8f92f80ae68a6e2966643fbe906bb35fe0aa93c58`，RUNPATH 只有
`$ORIGIN/../lib`。Jammy `ldd` 解析到移动副本的 `bin/../lib/libcantera_shared.so.3`；1-rank
和 2-rank 运行均退出 0，数值摘要与原始 consumer 一致。

重定位 validation v2 SHA-256 为
`7155067c27fb19dbeff83e33aeb8c1126782e610918ea5fb36c1428770976e9a`，执行日志 SHA-256 为
`835bc7f469fc380f561938207e3df8d6e660728ec94d2a2e4d3d93e226af4f4d`，退出 0。动态元数据、
数值输出和 binary exec 之后的 file trace 均不引用原 artifact/build/source 根；trace 只读取
移动后的共享库与 synthetic mechanism。因此 P0-2 记录的编译期 prefix 字符串未被运行时消费。

## 保留的 RED 与修订理由

历史失败未删除：

1. OpenMPI wrapper 找不到未版本化 `g++`，改为显式 `OMPI_CXX=g++-11`；
2. `-Werror` 误作用于 upstream headers，改用 system include 并关闭旧 MPI C++ bindings；
3. synthetic mechanism v1 的 `constant-cp` 组合在 `newSolution` 阶段得到零温度，v2 改为独立
   构造的 NASA7 常数系数；
4. ABI policy validator v1 没有给缺失 ABI=1 统一返回 policy code，v2 修正；
5. relocation validator v1 未接受动态加载器等价的 `bin/../lib` 路径，v2 按实际解析验证；
6. output validator v1 未拒绝字符串 `nan`，mutation suite 暴露后由 v2 增加数值词法门。

追加 manifest 的严格校验确认 P0-2 对象未改变，并重新计算 30 个 P0-3 外部证据 hash。
validator v2 SHA-256 为
`6bdb1b27564dc1a30449720533384d6eaf3e1ed4f4abc3377746ad3caef4f87b`，日志 SHA-256 为
`e7380e3156c3a18f23dd6001eb486cc54a6d653dcb308b6bcdb650d6e340bb19`，退出 0。

## 明确不证明的能力

- 未修改或链接 HUNDUN 产品，未执行 Stage 4 `4F-0`；
- 未建立 HUNDUN CMake/CPack/RPATH 集成或 `ChemistryBackend`；
- 未证明 reacting-flow coupling、operator splitting、retry/rollback、Checkpoint 或 diagnostics；
- 未证明真实航空煤油、汽油、点火延迟、火焰速度、PSR 或 COAST 相似性；
- 未访问 COAST、BOFFIN、私有研究数据或研究进程；
- 未 push、发布或改变宿主网络。
