# HUNDUN-FLOW Stage 4 P0 最终预研封印

- 记录时间：`2026-08-10T00:43:29+08:00`
- 结果：`PREFLIGHT_PASS`
- 结果范围：`pre_stage4_non_product_public_input_artifact_boundary_and_oracle_preflight_only`
- Stage 4 产品已接受：`false`
- 产品修改：`none`
- `4F-0`：`not_run`
- Stage 3 intake：`deferred_until_formal_stage3_acceptance_and_explicit_stage4_start`

`PREFLIGHT_PASS` 表示 P0 计划要求的来源锁、Linux CPU artifact candidate、独立 C++/线程/MPI/
relocation 边界、公开向量、双代理燃料接口候选和 Stage 3→4 intake 模板均已形成可复核证据。
它不是 Stage 4 产品接受，也不把双燃料候选链写成 bundle-ready。后者仍为
`PREFLIGHT_PARTIAL`：机制需要用户合法提供，PRF transport pairing 因重复记录被拒绝，
property packs 只是维护者环境生成的外部候选。

## 1. 冻结候选

| 字段 | 值 |
|---|---|
| P0 plan HEAD | `813670efc2a0ce6adb1a033fb98b7582b05a0fce` |
| sealed evidence candidate HEAD | `401a5ff285aa1a6283d68a241155bbfab405b9e4` |
| candidate parent | `10d20177f8be2d93fb654698ebeae6fd53296387` |
| candidate tree | `3e0819270bc8542a7f9b72b2e0c3103e496e167f` |
| P0 design SHA-256 | `6dd94a9f9aeb26361af0b47b58866cae9e4f1f6fcc69d73e06295a4d827352b8` |
| P0 plan SHA-256 | `a2c7967f141ba356d8bac12c84d9bc16d4fb8755dbc9b9885c4d8e92fe69b325` |
| `plan..candidate` diff SHA-256 | `44607f38e23c7c429f9cb097215e26b56f84fe4682e8f0d89ecdcda3a6ce385e` |
| changed tracked paths | `13` governance/policy/reference paths；`0` product source/build/test paths |
| diff size | `5160 insertions / 30 deletions` |
| candidate worktree before seal files | clean |

候选冻结后没有修改 P0-1--P0-6 消费的产品或测试路径。此 receipt 和两个 manifest 的
`p0_final` 对象是结果治理层；它们不改变经过复验的 candidate HEAD。

## 2. DCO 与提交边界

`plan..candidate` 共 20 个提交。逐提交核对 author、committer 和唯一现有
`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`，结果 `20/20 PASS`：

| Commit | Subject | DCO |
|---|---|---|
| `95630e9` | docs: activate Stage 4 P0 preflight | PASS |
| `d7471d2` | docs: clarify Stage 4 P0 isolation | PASS |
| `21b5a95` | docs: order Stage 4 P0 revalidation | PASS |
| `7fb02ed` | docs: lock Stage 4 P0 third-party inputs | PASS |
| `dd85bf1` | docs: harden Stage 4 P0 input plan | PASS |
| `09b0097` | docs: harden Stage 4 P0 input locks | PASS |
| `fe2aed5` | docs: prepare Stage 3 to Stage 4 intake | PASS |
| `312af28` | docs: align P0 artifact build with Cantera 3.2 | PASS |
| `d9e1579` | docs: split P0 fuel provenance from parse evidence | PASS |
| `2f280e2` | docs: lock Stage 4 P0 Boost builder input | PASS |
| `a19e4c8` | docs: clarify Stage 4 P0 probe artifacts | PASS |
| `67de677` | docs: isolate Stage 4 P0 Boost headers | PASS |
| `3ed7f79` | docs: audit Stage 4 P0 fuel data candidates | PASS |
| `14f7b52` | docs: freeze Stage 4-6 P0 oracle vectors | PASS |
| `0a788a9` | docs: lock Stage 4 P0 Doxygen builder inputs | PASS |
| `3139044` | docs: lock Stage 4 P0 sourcegen builder input | PASS |
| `3ef96dd` | docs: lock Stage 4 P0 YAML builder inputs | PASS |
| `8923015` | build: record Stage 4 P0 Cantera artifact candidate | PASS |
| `10d2017` | test: prove Stage 4 P0 Cantera boundary | PASS |
| `401a5ff` | docs: validate Stage 4 P0 fuel data candidates | PASS |

没有替第三方、worker 或旧历史伪造 sign-off。包含本 receipt 的 seal commit 在提交后由 Git 外
post-commit receipt 记录自身 HEAD、parent、tree 和 DCO，避免 tracked file 自引用。

## 3. P0-0--P0-6 完成矩阵

| Task | 结论 | 权威证据 |
|---|---|---|
| P0-0 identity/isolation | PASS | baseline receipt；唯一外部写根与 governance worktree 已冻结 |
| P0-1 provenance/BOM | PASS | input manifest、严格 validator、source/license/consumed-file records |
| P0-2 Release artifact | PASS | Ubuntu 22.04/GCC 11/ABI=1 artifact；653 files、2 symlinks；无机制/Python runtime payload |
| P0-3 C++ boundary | PASS | Release/Debug、1/2 threads、1/2 ranks、workspace isolation、ABI mutation、relative RPATH、relocation |
| P0-4 dual surrogate | COMPLETE；chain PARTIAL | 100/553 Reitz、65/363 Hybrid、89/480 iso-octane parse；两类液体 typed-pack range/mutation checks |
| P0-5 public vectors | PASS | 23 structured vectors、20 closed-form reviews、43 C++ checks、三工具链一致、39/39 document mutations |
| P0-6 intake template | PASS | `template_only`、29 deferred placeholders、14/14 omission mutations；未执行 accepted-state command |

P0-5 的两个 property identity vectors 保持绑定 P0-4 初始 source-provenance commit `3ed7f79`；
它们验证来源指纹，不是 P0-4 后续 output pack。P0-4 的 pack SHA 与范围验证单独登记，二者
不得互相冒充。

## 4. P0-7 低成本复验

最终 v2 runner SHA-256：
`e725ee14ef7d9fe3cb0425f1d2bb1c28bddff0ac9c7a7a01e2b0f950f08c0e83`。
它在 frozen candidate 上重新执行：

- input/artifact manifest validators 和 P0-3 30/30 append evidence hashes；
- artifact/relocated file trees、`readelf`、`ldd -r` 和相对 RUNPATH；
- Release 1/2-thread × 1/2-rank、ordinary Debug、relocated 1/2-rank；
- 三份 public mechanism C++ parse 和双 property pack generic C++ consumption；
- GCC 7、Clang 15、Jammy GCC 11 oracle outputs 与文档 validator；
- intake base validator 与 14 类 omission mutations。

结果：exit `0`。关键 seal：

| Evidence | SHA-256 |
|---|---|
| v2 overall log | `84d9a33ee831200f2038c12a55b1e38d850985b73d84191fe2967bea9bef537a` |
| v2 status | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| v2 logs manifest | `f83622eddb63a901d5ade48a66a8e7c2d933a880ca675742c6c5ec7289023cb0` |
| v2 primary manifest | `2af6bb03a03cf1a315cf47805373eff928836141e1cf39f1631a796c25f663f0` |
| v2 seal-input manifest | `cf827bdc6ec9067e72b8804e40ac4c591c532769b4ab9515c6578dfe90853cbe` |

v2 primary manifest 绑定 frozen candidate `401a5ff`。seal 阶段追加 `p0_final` 后，当前磁盘上的
input/artifact manifests 有意变化；因此最终校验从 `401a5ff` 的 Git blobs 核对这两个旧摘要，
并对 primary manifest 的其余路径直接执行磁盘 SHA 校验。只有这两个 result-only 路径允许按
candidate blob 复核，不能把它扩展为一般“旧证据可忽略”的规则。

### 4.1 保留的 final-validator RED

首次 P0-7 runner 在 input-manifest validator 立即 RED：维护者 Python import 曾在 frozen
wheel-extract 树生成清单外 `cpython-310.pyc`。原 wheel、源码和所有 manifest 内文件哈希
未改变。修复脚本只接受 `*/__pycache__/*.cpython-310.pyc`，把 657 个 cache files 移入
`logs/p0-7-retained-python-cache-contamination-v1/files/`，没有删除证据；随后逐树要求文件集
与原 manifest 完全相同。恢复脚本 exit `0`，其 retained-evidence manifest SHA-256 为
`3fad3b5a32a33a797175c2a095e50de035d63bb79aa2bf24abb26c30eebe9ab9`。

第一次 runner 的 exit `1`、空 stdout、validator error 和 timing 均保留，不计为 GREEN。
这次 RED 说明未来 maintainer builder 必须对 frozen inputs 使用只读 mount 或
`PYTHONDONTWRITEBYTECODE=1`；它不改变“普通 HUNDUN 不依赖 Python”的产品边界。

## 5. Consolidated review

主 agent 对 `813670e..401a5ff` 完成一次合并审查，覆盖 requirements、质量、来源/许可证、
host/product/Stage 3 隔离、ABI、能力声明和任务范围，不再以不同名称重复扫描同一 diff。

结论：

- 数学：公开向量的单位、符号、operation order、canonical `+0.0`、retry clock 和 tolerances
  与已登记公开方程一致；TCR 仍明确为 pre-COAST candidate status；
- ABI/package：artifact 绑定 Jammy/GCC 11/libstdc++/C++17/ABI=1/x86-64，普通 Debug 消费
  同一 Release ABI；没有 Clang/libc++ 混链；
- runtime：standalone/relocated consumer 没有 Python/Conda/libpython/site-packages；
- provenance：Cantera 与传递依赖保留 upstream identity、revision、license、notice、实际
  consumed-file closure；机理许可不从 Cantera BSD 推断；
- scope：13 个 tracked paths 全为 governance/policy/reference，未修改产品 CMake、schema、
  Checkpoint、diagnostics、driver、源码或测试；
- claims：没有 Stage 4/6、真实燃料、COAST equivalence、科学 acceptance 或全平台声明；
- retained failures：构建、validator、relocation、NaN 和 Python-cache RED 均保留并解释。

未发现阻断 P0 目标的 requirements、quality、license、ABI 或边界问题。双燃料完整链的缺口
是显式候选限制，不是被隐藏的 PASS 前提。

## 6. 可复用与不证明

| 可在正式 Stage 4 intake 后按哈希复核复用 | P0 不证明，必须由正式任务建立 |
|---|---|
| 官方来源、revision、archive/wheel、license/notice hashes | HUNDUN CMake/CPack/RPATH integration |
| Ubuntu 22.04/GCC 11/ABI=1 builder profile | `ChemistryBackend`、workspace API 和失败映射 |
| Cantera Release artifact candidate 和 relocation protocol | reacting species/enthalpy transport、Strang/PISO/IBM/WALE coupling |
| standalone C++ thread/MPI/lifecycle fixture与输出 schema | Checkpoint v4、Restart、diagnostics v4、driver/capability acceptance |
| public mathematical vectors、comparison classes、mutations | Stage 4/5/6 scientific acceptance 或性能 |
| mechanism/property candidate identities 与外部 parse/range evidence | real kerosene/gasoline fidelity、ignition/flame/spray validity |
| Stage 3→4 read-only intake command template | COAST equivalence 或 private oracle acceptance |

## 7. Repository 与进程边界

封印前只读快照：

- product `/home/wyf/code_dev/hundun-flow`：`main`，HEAD
  `ae3d08bbb220d1d3b28ec070d1cba9c33fb85877`，tree
  `833b633939765365a07fc8d49e34802399954399`，clean；
- governance main：仍只有用户原有的 untracked
  `docs/plans/2026-07-21-hundun-flow-stage7-neural-warm-start.md`；
- Stage 3 worktrees：旧 stage3 `52 modified / 258 untracked`，framework
  `2 modified / 0 untracked`，infrastructure `0/0`；P0 未 stage、clean 或修改它们；
- P0 exact-path process：none；P0 systemd user service：none；
- 未查看、枚举、停止或干扰无关研究进程；未访问 COAST/BOFFIN/private research tree；
- 未 push、发布、cherry-pick、执行 `4F-0` 或进入 Stage 4 产品实现。

最终 input manifest SHA-256：
`c96666be5c392925b582faab90652e9229c0f1126b966449f4c6394f66b91ac1`；最终 artifact
manifest SHA-256：`aa6eb2a1af5a3dc73a55b92f41999e88747ceca43b0b6762d39d484f45f9d452`。

P0 在此边界停止。正式 Stage 3 接受且用户明确启动 Stage 4 后，`4F-0` 必须从 accepted
governance/code/product identities 重新执行 intake；不得把 P0 branch 当作 Stage 3 ancestor
或直接进入产品 history。
