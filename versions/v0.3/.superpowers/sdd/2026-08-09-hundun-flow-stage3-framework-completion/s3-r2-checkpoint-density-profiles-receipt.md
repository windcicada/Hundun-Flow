# S3-R2 Checkpoint density profiles 验收回执

状态：`ACCEPTED`

accepted at：`2026-08-09T21:33:03+08:00`

accepted parent：`d118b4702b2e20fd28794b844a521b565a755cca`

accepted commit：`8f869aafd299617a36876c3380768b22fa8af2ce`

## Profile、协议与 transaction

Checkpoint v3 additive profiles 4--9 已完成：material/ideal-gas 与 IBM/WALE 的六种合法
组合均严格校验 required/forbidden object、section count 与 byte count。旧 profile 1 byte
fixture 保持不变；profiles 1--3 不追加 ideal section 字段。profiles 4--9 追加固定 30-byte
ideal-gas authority section：mode、p0、canonical optional target mass 与 non-wrapping revision。

read 在单一 collective-ready 边界前准备 FlowState replacement、IBM pressure authority 与
optional ideal-gas closure replacement；随后按 FlowState→IBM→closure 的 allocation-free、
MPI-free、`noexcept` publication 顺序发布。closed-domain restore 重新验证 history/committed
EOS 与 target mass，open-domain 固定 p0/no target；snapshot preparation fault、authenticated
EOS mutation、CRC、wrong profile、inactive negative zero 均失败中性。driver 依据
density×IBM×WALE 生成精确 profile 1--9 module view，旧 Restart deferral 已关闭。

fast continuation 使用 8-cubed、2 accepted steps、step 1 restart；每次 attempt 仍恰两次
PISO。正式 selector 固定为 `formal 12`、4 accepted steps、step 2 restart，并精确注册：

- `checkpoint-continuation-n12-r1`；
- `checkpoint-continuation-n12-r2`；
- `checkpoint-continuation-n12-r4`。

三条 formal row 仅执行 `ctest -N` 列出；未在 R2 运行。

## 验证证据

- R2 focused gate：8/8 PASS，real 128.67 s；
- R1 WALE regression：2/2 PASS，real 30.57 s；
- Clang 15/libc++ tests-off `hundun`：PASS；
- `test_checkpoint_v3_density_profiles` SHA-256：
  `40c4159d6a127c5ee3bd7c37f167ff4c6acee15d20b7073c1d685e5a8277509e`。

未运行 24/48/96-cubed、正式 12-cubed continuation 或正式科学矩阵；未访问私有源码、
研究数据或研究进程，未 push/publish。Task 11 阈值、两次 PISO、force sign、rollback、
Restart 与 MPI 一致性 authority 未修改。

## 权威与 DCO

激活 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

提交 subject：`feat: complete Checkpoint v3 module profiles`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`
