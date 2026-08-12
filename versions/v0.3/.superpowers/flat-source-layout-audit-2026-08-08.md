# HUNDUN-FLOW 扁平源码布局候选审计

> 本文件保留 2026-08-08 的批准前审计状态。后续签署提交、产品投影和切换前结论由 `.superpowers/repository-split-pre-switch-seal-2026-08-09.md` 接续；本历史记录不回写为新的最终状态。

## 结论

状态：`PAUSED_FOR_APPROVAL`

此前阻止 tests-free 产品投影的产品代码到 `tests/support/` 依赖已经修复。当前候选满足扁平目录、公共/私有头文件边界、产品目标包含目录权限和 tests-off 产品构建要求。Tasks 1--7 均已完成并通过独立审查，Task 8 的证据刷新已经完成，等待主 agent 复审和用户批准。

本轮没有提交、添加 DCO sign-off、切换目录、发布或 push，也没有运行数值、MPI、sanitizer 或其他高负载矩阵。

此前的阻塞审计原文已逐字保存在
`.superpowers/flat-source-layout-audit-pre-private-boundary-2026-08-08.md`。
当前文件是边界修复后的最终审计；两份文件分别记录修复前的
`PAUSED_FOR_APPROVAL` 停止点和本轮修复后的审计结论，不能互相覆盖。

## 冻结基线

- acceptance profile：`task11-core-functional-v1`
- Task 11 result：`CORE_ACCEPT`
- accepted product HEAD：`66080e324089599711fdb26082af9b330bfdb5ce`
- accepted product tree：`ab071a61f00eba9ec973beb0fe600066a33ef74f`
- migration worktree：`/home/wyf/code_dev/.worktrees/hundun-flow-flat-layout`
- branch：`governance/flat-source-layout`
- 当前 HEAD 仍是 accepted product HEAD
- governance layout commit：`PENDING_APPROVAL`
- product initial commit：`PENDING_APPROVAL`

工作树仍包含完整的未提交迁移候选。没有清理、覆盖或丢弃任何已有修改。

## 已完成的布局

- 产品公共头文件位于扁平的 `include/hundun/`。
- 产品实现和私有头文件位于扁平的 `src/`；私有头文件使用登记前缀和 `_detail.hpp`。
- 测试适配器、故障注入、typed snapshot 和测试夹具位于 `tests/support/`，不进入产品投影。
- 根 `CMakeLists.txt` 只保留项目级设置和子目录入口。
- 产品目标在 `src/CMakeLists.txt` 注册，测试与治理目标在 `tests/CMakeLists.txt` 注册。
- 产品默认 `HUNDUN_BUILD_TESTS=OFF`。
- 产品目标公开包含目录只有 `include/`，私有实现目录只有 `src/`；`hundun` 另有生成头目录。
- 仓库根目录和 `tests/` 不再是产品目标的包含目录。
- 中文规范 `docs/development/naming-and-style.md` 已补充产品代码与测试支撑的单向依赖规则。

## 原 blocker 的处理结果

旧审计记录的八个无条件产品到测试头文件依赖已经拆开：始终属于产品侧的声明或中性桥接放入 `src/*_detail.hpp` 或产品实现，测试专用枚举、快照和 typed adapter 留在 `tests/support/`。随后又按同一边界处理了 diagnostics、flow、FVM 和 linear 模块中的剩余依赖。

修复没有把测试头复制进产品，也没有给产品目标加入仓库根目录或 `tests/` 包含权限。tests-off 路径不依赖测试类型；tests-on 的 raw hook 或中性记录由测试适配器转换为原有测试接口。

## accepted-blob 对照

`.superpowers/flat-source-layout-move-map.tsv` 的 208 对移动逐一与 accepted HEAD blob 比较。比较同时保留字节级结果和删除预处理 `#include` 行后的归一化结果：

- 20 对为字节相同移动；
- 133 对只改 include 行，其中 2 对是 tests/support 文件；
- 55 对含非 include 行变化：31 个产品文件和 24 个测试适配器。

31 个产品文件没有继续标为 `include_path_only`。逐文件审查结果如下：

| 修复簇 | 产品文件 | 审查结论 |
| --- | --- | --- |
| Task 2 | `src/app_case_config_broadcast.cpp`, `src/bc_basic_boundary.cpp`, `src/exec_execution.cpp`, `src/fvm_cell_centered.cpp`, `src/fvm_immersed_reconstruction.cpp`, `src/rt_field_storage.cpp`, `src/rt_halo_exchange.cpp`, `src/rt_mpi_context.cpp` | 测试专用类型改为宏控 raw hook、中性标量或产品类型；默认产品分支、MPI 顺序和错误语义不变 |
| Task 4A/4B | `src/diag_checkpoint_v2.cpp`, `src/diag_ideal_gas_closure.cpp`, `src/diag_material_density_piso.cpp`, `src/diag_material_density_transport.cpp`, `src/diag_structured.cpp`, `src/diag_time_control.cpp` | typed diagnostic fault/work 类型留在 tests/support；产品侧使用中性编码和 raw adapter；记录、schema 和 collective 语义不变 |
| Task 5A/5B/5C | `include/hundun/flow_constant_density_piso.hpp`, `include/hundun/flow_immersed.hpp`, `src/flow_adaptive_time_control.cpp`, `src/flow_adaptive_time_control_detail.hpp`, `src/flow_checkpoint_v2.cpp`, `src/flow_constant_density_piso.cpp`, `src/flow_density_closure_detail.hpp`, `src/flow_state.cpp`, `src/flow_ideal_gas_closure.cpp`, `src/flow_ideal_gas_piso.cpp`, `src/flow_material_density_piso.cpp`, `src/flow_material_density_transport.cpp`, `src/flow_momentum_predictor.cpp`, `src/flow_immersed.cpp` | 测试 friend/authority 入口改接产品侧 detail/raw authority；PISO、rollback、checkpoint、pressure/flux/force authority 和公开运行接口不变 |
| Task 6A/6B | `src/fvm_immersed_boundary_authority_detail.hpp`, `src/fvm_immersed_operator.cpp`, `src/lin_preconditioners.cpp` | FVM authority 和 Jacobi 状态观察改为产品侧中性记录/raw seam；符号、所有权、缓存发布与失败回滚不变 |

另有两个没有 accepted 对应 blob 的新产品私有头：

- `src/flow_immersed_access_detail.hpp`
- `src/lin_preconditioners_detail.hpp`

它们只承载上述私有边界提取所需的中性声明，manifest 分类为 `new_private_detail`。

24 个移动后的测试适配器分类为 `private_boundary_test_adapter`。它们负责恢复原有 typed test API，不进入产品投影。已有同路径测试文件中，139 个只改 include 行；19 个 CMake/source-policy 或路径合同文件更新了布局路径；`tests/unit/test_basic_boundary_header_contract.cpp` 更新了被检查的产品路径；`tests/unit/test_preconditioners.cpp` 新增 mutation-sensitive raw seam 验证，单独分类为 `test_validation_extension_not_exported`。

上述 31 个产品文件分别由 Tasks 2、4A、4B、5A、5B、5C、6A、6B 的实现报告和独立 review 覆盖。没有发现数值算法、阈值、PISO corrector 数量、配置键、Restart/diagnostic schema、单位、符号、MPI collective 顺序或默认运行路径变化。

## migration manifest

manifest：`.superpowers/flat-source-layout-manifest.tsv`

数据行：378。`.superpowers/` 下的计划、历史报告和审计证据明确不属于产品迁移 manifest，未混入产品文件清单。

当前分类：

| 数量 | change_kind |
| ---: | --- |
| 139 | `test_include_path_update_not_exported` |
| 131 | `include_path_only` |
| 31 | `private_boundary_extraction` |
| 24 | `private_boundary_test_adapter` |
| 20 | `identical_move` |
| 19 | `governance_path_update_not_exported` |
| 5 | `governance_fixture_not_exported` |
| 2 | `new_private_detail` |
| 2 | `include_path_only_not_exported` |
| 1 | `test_validation_extension_not_exported` |
| 1 | `governance_registration_and_authority` |
| 1 | `documentation_and_boundary_policy` |
| 1 | `cmake_registration_include_authority_and_test_definitions` |
| 1 | `cmake_registration` |

每个非 `-` old hash 都从 accepted HEAD blob 重新计算；每个 new hash 都从当前候选文件重新计算。逐行校验结果为 0 个 mismatch，没有重复 new path，208 个 move-map 目标全部有记录，新增的两个私有头和四个治理 fixture 也全部纳入。

manifest SHA-256：`d5bb79e246d832b2e30642e4102ef62929871ace47769a145b2628aa8bac805e`

## Tasks 1--7 证据

| Task | 结果 | 报告 SHA-256 |
| --- | --- | --- |
| 1 | 边界 RED 与 tests-free projection fixture 建立，review PASS | `e5aaf9d7dd473c3143f443854eefec9dc826c5031caf17a848eab79e3ac9e59b` |
| 2 | 八个硬依赖拆分完成，review PASS | `de6e707e2a5375e54415ca5b33e196e33c36ddce7ddd25e06c8869dacfdcd9e8` |
| 3 | macro-hidden tests 路径 mutation-sensitive，review PASS | `8850235364d7f2329c4f5bc81e0db374e30f67f9006830a7a18db87f28c5ae37` |
| 4A | checkpoint/structured/material-PISO diagnostics，review PASS | `8e09837773cbf368b29df6fabc676aecbaca07a645ca559ff0c4e8c8b3238b58` |
| 4B | ideal-gas/material transport/time diagnostics，review PASS | `d20a130c466617f0d92b4bb8f48fac41955fc8227a8febd26029f59b5faae198` |
| 5A | flow state/checkpoint，review PASS | `ad535bcf0e174f6cd86e6de248e982f9fd5a238bc66363562ce7f4c59a2761f8` |
| 5B | adaptive/PISO/density/momentum，review PASS | `777f474a254a1c71d07d9c65a9468512412051092f78c68dd238e8aed66faf86` |
| 5C | immersed flow authority adapter，review PASS | `7475e193b53c9e6c700979ee1d06146a20ccffab54e83b31fdeae47e2c380a38` |
| 6A | FVM immersed authority bridge，review PASS | `abed70bc56388b179ad6fd388b76d7f018c6b3bc57a36c35e7077d043f059c11` |
| 6B | Jacobi/preconditioner boundary，mutation 后 review PASS | `3e380a3140ebdd4594d95c0c9dfb066f40752b593f429df43bd8e83593a1ac14` |
| 7 | configured target graph + 五类 CMake authority mutation，独立复审 PASS | `f4a2f20944b5a93881c60649385bf62984712411d4dd3749fe0cb83af9bcf1f4` |

旧版 blocker 审计已保留为
`.superpowers/flat-source-layout-audit-pre-private-boundary-2026-08-08.md`；
本文件只记录边界修复后的最终审计，不再把当前修订内容当作旧版记录。

## 本轮低成本验证

本轮直接运行：

1. `source_layout_fixture.cmake`：exit 0。
2. `source_layout_mutation_fixture.cmake`：exit 0。
3. `cmake_include_authority_fixture`：PASS，0.01 s。
4. `cmake_include_authority_mutation_fixture`：PASS，1.39 s。
5. `rg -n 'tests/|tests/support' src --glob '*.{cpp,hpp,h}'`：exit 1，无匹配。

CMake authority mutation 覆盖 direct root、direct tests、directory-scope root、target property tests 和 included-module/transitive root 五类注入。正常 tests-on configure 还会读取真实 target graph 的 `INCLUDE_DIRECTORIES`、`INTERFACE_INCLUDE_DIRECTORIES` 和可达 link target，静态 CMake 扫描不再是唯一证据。

## 复用的构建与产品投影证据

没有重跑 product projection 或产品构建。

Task 7 已建 tests-free projection：

- 根：`/tmp/hundun-task7-tests-on.YskPEg/product-projection-fixture`
- `CMAKE_BUILD_TYPE=Release`
- `HUNDUN_BUILD_TESTS=OFF`
- 不含 `tests/`、`AGENTS.md`、`CONTRIBUTING.md`、`.github/`、`.superpowers/`
- 编译命令中的项目包含目录只有 projection 的 `include/`、`src/`、`src/generated` 和 `third_party/yyjson`
- `hundun` binary SHA-256：`405079116d648afa5d8681ad616df4b8bab60420a576e6adc0af3dc2c4b2ca8d`

主 agent 已对 projection 的 `CMakeLists.txt`、`src/`、`include/`、`cmake/`、`third_party/`、`NOTICE` 和 `LICENSE` 与当前候选做 `diff -qr`，无输出。随后 Task 8 只修改了不参与编译的中文规范和治理证据，不使该产品构建失效。

Task 7 tests-off Debug build：

- 根：`/tmp/hundun-task7-tests-off.tdMitD`
- `HUNDUN_BUILD_TESTS=OFF`
- `hundun` binary SHA-256：`9f299b76af48e5b1ef8965965cf7aefef78db50e32f6fb387024c1fa51d9eae1`

Tasks 2--7 报告中的 clean tests-on/off、focused unit/header/linkage 和低成本 rollback/collective 证据按相同文件 hash 复用。没有因 manifest、审计或中文规范变化重跑这些测试。

## 文档治理

`docs/development/naming-and-style.md` SHA-256：`d62f98905f91bd95530b3fc822002094725ac1dd068136e715852d8cfa1d0423`

新增段落明确：

- 依赖只能从 tests 指向产品；
- 产品文件不得直接或通过宏引用 `tests/`；
- 产品侧中性 raw/detail 声明与 tests/support typed adapter 分开；
- tests-off 默认路径不依赖测试状态；
- 产品目标不得取得仓库根目录或 tests 包含权限。

文本按 `humanizer-zh` 和 `shuorenhua` 的 docs/minimal 规则回读。路径、宏、术语和技术条件保留原样，没有新增未核实的能力声明。

## 进程和资源状态

审计开始和结束时均只按本 worktree、Task 7 临时构建根以及明确的 HUNDUN 构建/MPI命令过滤进程。没有发现正在运行的 HUNDUN 构建、CTest、MPI 或数值模拟。没有检查、停止或干扰其他研究进程。

## 限制和批准边界

- 候选仍未提交，最终 commit、tree、diff SHA 和 DCO 只能在用户批准后冻结。
- 现有 accepted-tree 编译 warning 未在本布局任务中处理。
- 本轮没有重跑数值、MPI、sanitizer 或长矩阵；布局和私有边界修改使用 Tasks 2--7 的同文件构建/低成本测试证据。
- migration manifest 不包含 `.superpowers/` 历史和审计文件；这些文件作为治理证据保留在当前工作树。
- 未执行仓库目录切换、产品仓库初始化、发布或 push。

下一步必须先由主 agent 完整复审 Task 8 的治理 diff，再向用户报告。未经用户批准，不创建迁移 commit，也不切换目录或进入后续仓库拆分动作。
