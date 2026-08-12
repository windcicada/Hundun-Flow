# S3-V2 product projector repair packet

状态：`ACCEPTED`

发现于：S3-V2 Step 2 preflight

rejected candidate：`7534c614f6fb81c0d03c596f871d1b71533c2ce6`

## Defect

`tests/cmake/stage3_product_projection.cmake` 只验证既有 manifest 存在，未实现已激活
v2 计划要求的产品投影、产品专用 `VERSION=0.2.0`、最终 manifest 或 fail-closed
校验。`test_stage3_product_projection_contract` 只搜索脚本文本中的 manifest 名称，因此
未能杀死空实现。

旧候选的 V0/V1 manifests 保留为 rejected history；修复提交产生新 exact candidate 后，
从 S3-V0 Step 1 重新冻结和验收。

## Allowed files

- Modify: `tests/cmake/stage3_product_projection.cmake`
- Modify: `tests/cmake/stage3_product_projection_contract.cmake`
- Create: `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/s3-v2-product-projector-repair-packet.md`
- Create: `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/s3-v2-product-projector-repair-receipt.md`

不得修改产品数值代码、Task 11 authority、阈值、PISO 次数、force sign、rollback、
Restart 或 MPI 一致性要求。

## Required behavior

1. `HUNDUN_BASE_MANIFEST` 和 `HUNDUN_OUTPUT_MANIFEST` 相对路径只相对
   `HUNDUN_GOVERNANCE_ROOT` 解析。
2. 投影前要求 governance/product 均为 clean Git worktree，product 为零 remote。
3. base manifest 精确认证当前产品 baseline；未知 tracked product path、缺失路径、blob
   不匹配、重复/非法 path、symlink 均 fail closed，且校验失败前不修改产品。
4. 产品路径只允许既有 base manifest 路径，以及治理树中受控的 public/product roots；
   `tests`、`.superpowers`、`.github`、cases、设计/交接/计划/参考资料不得进入产品。
5. `CMakePresets.json` 保留 product tests-off 版本；`VERSION` 仅在产品写为 `0.2.0`；
   其余目标文件与 governance candidate 相同。
6. 输出四列 final manifest：`path/product_blob/governance_blob/relation`，包含所有且仅包含
   投影后的 tracked product paths。

## TDD and verification

- RED：真实临时 Git fixture 调用占位 projector；必须因未生成产品变更/输出 manifest
  失败，而不是因 fixture 语法或配置失败。
- GREEN：同一 fixture 验证 copy/preserve/override/output manifest；mutation 验证未知产品
  path 和 baseline blob 错误均 nonzero，且不发生部分写入。
- 运行 focused contract、governance group、`git diff --check`、DCO/source-policy；签署一个
  repair commit 后重启 S3-V0。
