# 当前源码整理与发布范围

整理基线：`e484cbbf302f971c92b1d32592ffe86079d1a6dc`。
推送前查询的远端 main：`0127616f3716b7d7d6fcb4a00831129c2d6631b6`，是该基线的祖先。

## 移除与保留

- 移除 `versions/v0.3` 的 760 个历史实现文件，以及根目录未参与当前构建的
  `src/`、`include/` 共 207 个文件。
- 移除三个仅服务旧源码的 CMake helper 和四份退休 Stage 2/3 能力/合同页面。
  已接受历史保留在 Git 中，不改写历史提交。
- 当前 `versions/v0.4` 的源码、头文件、测试、构建清单与 runner 完整保留。
  `v04_*`、`hundun::v04`、方法/文件格式和脚本名称是当前兼容合同，不做机械改名。
- 保留当前依赖的 `third_party/yyjson`、许可证、研究引用、历史 receipt、
  验证报告及原始附件；不为了消除版本文字而改写签名或历史哈希。
- 根构建只提供当前实现，旧的显式 `HUNDUN_SOURCE_VERSION=v0.4` 调用仍可用；
  退休或未知值明确拒绝。presets 不再呈现多源码选择。
- README 和使用文档聚焦产品 1.0.0，去除旧 CLI/schema、未接入的多物理耦合步骤
  及过时域/网格的当前能力声明。最小示例与既有 8³ 回归输入逐字节一致。

## 本轮实际检查

1. `cmake -P cmake/tests/current_source_layout_test.cmake` 通过：执行真实 selector，
   在测试中截获 `add_subdirectory`，不配置编译器、不启动 MPI；检查当前入口、
   v0.3/v0.2/v9 拒绝、退休目录移除、当前依赖存在及示例/夹具 SHA256 一致。
2. `cmake --list-presets` 成功读取 debug/release/asan/ubsan。
3. 新改文档的 69 个相对链接存在；当前实现的 248 个双引号 include 均可解析
   （生成头 `app_identity_build.hpp` 由原 CMake 规则生成）；JSON 解析通过。
4. `git diff --check HEAD` 通过。相对整理基线，`versions/v0.4`、`tools`、
   `third_party`、研究/验证证据目录无修改。
5. 冻结目录 `method-frozen-module-reviewed-20260907/FROZEN.sha256` 全部核对通过；
   长测服务仍为 active/running，MainPID=113610，没有暂停、替换或重启动。

检查中曾发现最小示例与夹具的物性文件首行注释不同，统一后哈希对照通过，
未改变物性数据。旧示例 JSON 的 `stl_file` 和缺失求解器字段已替换为当前测试输入。

本轮是仓库/构建入口/文档清理，没有重新编译或运行数值/MPI 测试，以免与
128 核冻结长测争抢资源。不能把脚本级路由检查说成完整配置、链接或科学验收；
先前数值回归仍只对应[原验收提交](../verification/2026-09-07-exclusive-module-acceptance.md)。
构建输入与 Git 身份因本次清理而改变，未来重新构建应生成新的 build manifest，
不得借用冻结程序的二进制身份；当前长测继续使用其原始已验证冻结程序。
