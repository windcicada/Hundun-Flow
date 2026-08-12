# S3-V2 product policy leak repair packet

状态：`ACCEPTED`

发现于：S3-V2 product pre-commit private-token scan

rejected candidate：`fe9065f8559e1367e8e112505bdd565f108d217f`

## Defect

正式 projector 将 `docs/development/source-policy.md` 归入通用 public development
目录。该文件记录 private legal comparison baseline 与内部来源边界，是治理文档；产品
pre-commit scan 正确拒绝它。其余 272 条 projected path 没有命中 private token。

本次未提交 product candidate 和 final manifest 均为 rejected attempt。修复产生新 exact
candidate 后，从 S3-V0 Step 1 重启；旧 V0/V1 evidence 保留为 rejected history。

## Allowed files

- Modify: `tests/cmake/stage3_product_projection.cmake`
- Modify: `tests/cmake/stage3_product_projection_contract.cmake`
- Create: `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/s3-v2-product-policy-repair-packet.md`
- Create: `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/s3-v2-product-policy-repair-receipt.md`

## Required behavior

- 保留既有 public development path policy；
- 精确排除 `docs/development/source-policy.md`；
- contract 在 governance fixture 中创建该文件，并断言 product 与 final manifest 都不含它；
- tests、`.superpowers`、private/token、baseline blob 与 fail-before-write 规则不变；
- 不修改产品数值代码、Task 11 authority、阈值、两次 PISO、force sign、rollback、
  Restart 或 MPI 语义。

## TDD and acceptance

- RED：当前 projector 将 fixture source-policy 写入产品，contract 精确失败；
- GREEN：同一 contract 通过，既有 copy/preserve/override 与三类 fail-closed mutation
  仍通过；
- 清理 rejected product worktree/output manifest，签署 repair commit，重启 V0。
