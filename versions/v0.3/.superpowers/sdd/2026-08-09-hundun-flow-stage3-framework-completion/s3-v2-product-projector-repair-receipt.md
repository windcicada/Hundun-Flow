# S3-V2 product projector repair receipt

状态：`ACCEPTED`

accepted parent：`7534c614f6fb81c0d03c596f871d1b71533c2ce6`

## Defect and disposition

S3-V2 Step 2 preflight 发现 `stage3_product_projection.cmake` 只有 manifest-name
占位检查。它会对权威投影命令返回 0，却不复制产品文件、不保留 tests-off preset、不写
产品专用 `VERSION=0.2.0`，也不生成最终 manifest。旧候选的 V0 31/31 与 V1 26/26
虽已通过并完成 57/57 manifest 身份审计，但不能证明 V2；这些证据保留为 rejected
history，不绑定修复后的 exact candidate。

## RED / GREEN

contract 改为在临时、独立、clean、zero-remote Git fixtures 上运行真实 projector。
占位实现的 RED 为：

```text
product projection omitted .../product/include/hundun/new.hpp
0% tests passed, 1 tests failed out of 1
```

GREEN 验证以下可观察行为：governance 内容复制、新 public header 投影、产品
`CMakePresets.json` 原样保留、产品 `VERSION` 精确为 `0.2.0`、最终四列 manifest
生成，且 product 中没有 tests 或 `.superpowers`。三个 mutation 分别要求精确错误原因：

- `unexpected tracked product path`；
- `product baseline blob mismatch`；
- `illegal product projection path`。

每个失败 mutation 都在写入前终止，并检查产品 sentinel 未变化。

## Implementation boundary

projector 只接受四个显式参数；两个 manifest 路径只相对 governance root 解析。它先验证
两个 Git worktree、clean 状态、product zero remotes、base manifest schema/关系/重复与
非法路径、完整 product tracked set、baseline blobs、governance target set 和 symlink，
然后才执行：

- base manifest 路径与明确 public/product roots 投影；
- `CMakePresets.json` 使用 `product_tests_off_preset`；
- `VERSION` 使用 `product_only_override` 并写入 `0.2.0`；
- 其他路径使用 `identical`；
- tests、`.superpowers`、`.github`、cases、设计/交接/计划/参考资料及 private/token
  路径排除；
- 生成所有且仅有目标产品路径的四列 final manifest。

没有修改产品数值代码、public API/ABI/schema、Task 11 pressure/operator/final-flux/force
authority、阈值、两次 PISO、force sign、rollback、Restart 或 MPI 语义。

## Fresh verification

```text
ctest --test-dir build/stage3-final-debug --output-on-failure \
  -R '^(test_stage3_product_projection_contract|test_stage3_acceptance_contract|test_stage3_source_policy|provenance_clean)$'
4/4 PASS

git diff --check
exit 0
```

SHA-256：projector
`605603d4c1fb80943c814b4d637a87f9d771f39d1a8804d5f27e29348c325244`；
contract
`9fc422865e3cdfdf20ad9e6a54e2e8e35f1116e3951c2e778a9bf6b428bde640`；
packet
`8ce7ae76655b7ca7de0334441aae4ede96367637d0b37d3259e1041b16654bab`。

修复提交后从 S3-V0 Step 1 重启 exact-candidate 验收。未访问私有源码或研究数据，未干扰
研究进程，未运行 96³，未 push/publish。

提交 subject：`fix: implement Stage 3 product projection`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`
