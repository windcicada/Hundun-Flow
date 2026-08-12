# HUNDUN-FLOW 仓库拆分切换前封存

## 结论

状态：`READY_FOR_DIRECTORY_SWITCH`

本记录封存仓库拆分、产品投影和最小迁移验收。它不改变 Task 11 的数值行为、阈值、selector 或科学接受结论，也不进入 Stage 3 后续功能实现。

封存时间：`2026-08-09T01:34:08+08:00`

## 基线与提交

- Task 11 accepted product HEAD：`66080e324089599711fdb26082af9b330bfdb5ce`
- Task 11 accepted tree：`ab071a61f00eba9ec973beb0fe600066a33ef74f`
- Task 11 acceptance：`CORE_ACCEPT`
- Task 11 scientific status：`ACCEPTED_FOR_CURRENT_REQUIREMENTS`
- 签署扁平化提交：`15001fbbf2acba3c61adae54f829bff01266ba74`
- 扁平化 tree：`0f3a38b058ad0c59f54a877b7a6037e3d67619b9`
- 迁移基线提交：`b06590a86d40369d2b178547de46c69cbb4f78ea`
- 术语收口提交：`56b716314e661cf4c7917642cc584bc5e6c80ebe`
- governance pre-switch tree：`b7c946d91c5fe303a51b01f9774d3fced959aa6f`
- product root commit：`ae3d08bbb220d1d3b28ec070d1cba9c33fb85877`
- product tree：`833b633939765365a07fc8d49e34802399954399`

`15001fb` 只给未接受的布局提交补上获授权的 DCO。修改前后 tree 均为 `0f3a38b058ad0c59f54a877b7a6037e3d67619b9`；Task 11 及更早历史未改写。上述 governance 三个迁移提交和 product root commit 均含：

```text
Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>
```

## 产品投影

- product candidate：`/home/wyf/code_dev/hundun-flow-product.candidate`
- governance candidate：`/home/wyf/code_dev/hundun-flow-governance.candidate`
- product tracked files：245
- product public documents：40
- product public headers：72
- product Git roots：1
- product Git commits：1
- product remotes：0
- product tests、治理记录和历史验收日志：未投影
- product `CMakePresets.json`：全部默认 `HUNDUN_BUILD_TESTS=OFF`，不含 test preset

投影 manifest：`.superpowers/product-projection-manifest-2026-08-09.tsv`

- 数据行：245
- 与 governance 相同 blob：244
- 受控差异：1，仅 `CMakePresets.json` 的 product tests-off preset
- manifest SHA-256：`bd736f22a78f68ac857353d87ac1ea4e3ee4dc6cd06cb1d19edfa30f442fbf64`

稳定兼容标识 allowlist：`.superpowers/product-stable-identity-allowlist-2026-08-09.tsv`

- allowlist SHA-256：`e5387a5cbd60746a203b3a9c50748681bde7ffcfc61ba1269cc15d1031ce7da0`
- 保留内容仅为 field ownership、numerical fingerprint、已登记 field identity 和 legacy VTK title。

产品 CMake guard 已移除治理专用敏感词表和具体内部目录名。对应 token mutation 由不进入产品投影的 `tests/cmake/HundunPrivateProvenanceGuard.cmake` 承担。

## 审查结论

`Task11..56b7163` 的迁移差异已由主 agent 完整审查。产品数值文件的变化限于：

- 公开类型和私有调用方由开发阶段名改为领域名；
- 用户可见错误和 MPI operation label 改为领域名；
- 测试、诊断目标和私有符号去除任务编号；
- 版本、安装规则、产品文档和 tests-free 投影边界；
- 测试专用 provenance 规则迁出产品 guard。

没有修改公式、阈值、PISO corrector 数量、配置 key、Restart/diagnostics schema、field identity、单位、符号、MPI collective 顺序或默认数值路径。公开文档明确区分 Task 11 IBM 数值核心与尚未完成的 schema 3 driver、Checkpoint v3 和组合 diagnostics。

## 最小迁移验收

工具链：Clang 15.0.6、libc++、Open MPI 3.1、C++17；构建并行度 `-j32`。

- governance Release/tests-off 干净构建：通过。
- 安装、`HUNDUN-FLOW 0.1.0`、公共头 standalone 编译、`nm`/`ldd`：通过。
- 最小模板 `--validate`、`--print-resolved`、单 rank 一步运行：通过。
- include-authority fixture 与 mutation fixture：2/2 通过。
- 受影响的低成本 unit/header/policy：7/7 通过。
- provenance/source-policy：39/39 通过。
- schema v3 broadcast 1-rank 与 failure：2/2 通过。
- 受领域重命名影响的 IBM 测试目标：只做编译闭包，全部通过；未执行数值矩阵。
- product candidate 独立 Release/tests-off 构建、安装、版本、校验和单步运行：通过。
- 48³、96³、warped/prism、sanitizer 和大型 MPI：未运行。

关键日志 SHA-256：

| 证据 | SHA-256 |
| --- | --- |
| governance configure | `6cd8a788e62d249308b2204269ab58c5dc9f109fb9bf7b561d78965a3d48609f` |
| governance build | `5f535f513a19b2f1d88c7695d9174dfe3d9864f5943ecdd5b9db65f667514713` |
| governance install | `23ed7803025800979493b6e74b5fd78010868bf5edbaf25d0842917d2cace4b6` |
| focused 7-test CTest | `0d7fc8c3dc3feaa9680c2901f7dbe8b0129322990302d5a5c07a261a88f942ac` |
| provenance 39-test CTest | `6a44c6fac50e64f4b0b79e2c35e598678463a2fb131fab363f5547c9f56b0490` |
| v3 broadcast CTest | `5dbc47fc7cdb1104b3ce02571793000d1efd81b87ca386286c3e14221db7be41` |
| product configure | `62914d24dfc40243611fdb26082af9b60e306362503286b64db15845003445` |
| product build | `6bd3bae09799ab8881504f689f40f132322090a210417e30a64d61b1f691525c` |
| product install | `85ab7615aa4a32bc3717f3cecbbd56fb97efb3141d49f04cf5697ae8d4c1f96a` |
| product one-step run | `5f14de00d38425c4d2c3990afe0d670eaf92c1cabe74ff39d1c1e0cc00960450` |

二进制 SHA-256：

- governance candidate binary：`f006eead5fcef3aed9473a125942bff204659a5fb7aaa449dc3eddeabdeb3815`
- product candidate binary：`c3a4d8cbbaf8db1fe26d9865b84adfb308414581fc255554f3ecc579d54a3957`

二者由不同源码根构建，因此调试路径会使二进制哈希不同；产品投影 manifest 证明所有编译输入 blob 相同，唯一投影差异是未参与命令行配置的 preset 文件。

## 扫描和边界

产品首次提交前和提交后均完成 tracked-text 扫描；提交后又扫描完整单提交历史。结果：

- 私有绝对路径：0
- 私有来源 token：0
- 内部治理目录名：0
- `Task N`/`taskN` 产品术语：0
- 未列入 allowlist 的 `Stage 1/2/3` 产品术语：0
- product remote：0
- NOTICE：与法律冻结文本完全一致
- 第三方 yyjson 来源和 MIT 许可证：已保留

未访问私有参考源码或研究数据，未检查、停止或干扰研究进程，未 push 或发布。

## 切换边界

记录生成时，governance 工作树除本 manifest、本 seal 和旧审计的 superseding notice 外无修改；product candidate 干净；未发现 HUNDUN-FLOW build、CTest、MPI 或数值进程。下一步只允许：

1. 签署本治理记录；
2. 将 mixed repository 迁到最终 governance 路径；
3. 将 product candidate 迁到 `/home/wyf/code_dev/hundun-flow`；
4. 修复并验证所有 linked-worktree `.git` 绝对指针；
5. 从本治理基线创建 `coast/stage3-framework-completion` 和新的 Stage 3 工作树。

旧 `/home/wyf/code_dev/.worktrees/hundun-flow-stage3` 及其脏修改和证据必须原样保留，不再作为新开发基线。
