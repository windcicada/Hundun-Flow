# S3-V2 product policy leak repair receipt

状态：`ACCEPTED`

accepted parent：`fe9065f8559e1367e8e112505bdd565f108d217f`

## Defect and rejected attempt

正式 projector 在 273 条输出中包含 `docs/development/source-policy.md`。产品
pre-commit private-token scan 发现只有该路径命中内部 comparison baseline 名称；文档
内容明确属于治理来源政策，不是产品用户文档。因此未提交的 product worktree、branch 与
final manifest 全部作为 rejected attempt 清理，产品根 `main` 保持
`ae3d08bbb220d1d3b28ec070d1cba9c33fb85877`、clean、zero remotes。

同一次 scan 的 273 个 blob mismatch 是检查器把 Git blob SHA-1 当成文件 SHA-256，已由
`git hash-object` 样本确认；LICENSE 首字符断言也是检查器错误。这两类误报没有转化为
产品或 projector 修改。

## RED / GREEN

contract 的 governance fixture 新增 `docs/development/source-policy.md`。修复前真实
projector 运行后观察到：

```text
governance source-policy document leaked into product
0% tests passed, 1 tests failed out of 1
```

GREEN 只在 `product_path_allowed()` 精确排除该路径；其他 public development 文档仍按
原 policy 投影。contract 同时断言产品 filesystem 与 final manifest 均不含该文件，且
既有 copy/preserve/override、unexpected tracked path、baseline blob mismatch、非法 path
和 fail-before-write mutation 全部继续通过。

## Fresh verification

```text
ctest --test-dir build/stage3-final-debug --output-on-failure \
  -R '^(test_stage3_product_projection_contract|test_stage3_acceptance_contract|test_stage3_source_policy|provenance_clean)$'
4/4 PASS

git diff --check
exit 0
```

SHA-256：projector
`d0306fc405893aadf1ea350fd6991541cf7dbf7f20bec738d96558a740e4b479`；
contract
`af69b089cb4aad7cda39e161b62d62e4b65d55e8a08ebdea063cedcc185eba18`；
packet
`8baffc05be39230c22e9184e4794061ca63d9f00df3642cb26a2e0e6df22cc95`。

实际修改只有 packet/receipt、projector 与 contract。没有修改产品数值代码、public
API/ABI/schema、Task 11 authority、阈值、两次 PISO、force sign、rollback、Restart 或
MPI 语义。未访问私有源码或研究数据，未干扰研究进程，未运行 96³，未 push/publish。

修复提交后从 S3-V0 Step 1 重启 exact-candidate 验收。

提交 subject：`fix: exclude governance policy from product`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`
