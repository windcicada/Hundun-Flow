# S3-S1 最终科学 selector 冻结回执

状态：`ACCEPTED`

accepted parent：`05dc59884ccb44c0352c132f71f1a651e98845c8`

accepted at：`2026-08-09T18:42:36+08:00`

## 实现边界

本 task 只修改测试支持、科学 producer 与 `stage3_science_registration.cmake`，没有修改
product flow、Task 11 force authority、PISO 次数、阈值、rollback、Restart 或 MPI 一致性。
新增 canonical scientific row schema：每个数值字段均携带 availability bit；unavailable
字段不输出数值占位。row 同时绑定 cells、ranks、process grid、steps、dt、final time、WALE
identity 与 pass/fail status，parser 拒绝缺项、重复项、未知项、非法网格和不一致 final time。

WALE TGV 通过既有 body-fitted product facade，使用光滑周期 cell-average 初值、常物理黏度、
固定 WALE controls 和固定物理终止时间。12/24/48 的 dt 随 h 缩放，首步使用既有权威
backward-Euler startup，后续使用 BDF2；每一步仍精确两次 pressure corrector，且从 product
test-access 证明 WALE evaluation count 精确为 1、identity 与 accepted report 一致。均匀
spanwise Galilean 速度只为严格二维零 z 动量方程提供非退化归一化尺度，不改变二维
Taylor--Green 解或压力。

单 rank convergence selector 同时执行：

- 对解析 cell-average velocity 与 gauge-normalized pressure 计算 12/24/48 三层误差及两个
  observed-order segment；
- 对 24→12 和 48→24 数值场执行真正 cell-average restriction，并要求相邻自差阶次沿用
  `>=1.8`；
- 将同一次 accepted `WaleSummary` 转为全局体积加权 `nu_t` L2，以 WALE 的 delta-squared
  连续极限作为误差；不做第二次 WALE evaluate、替代 PDE 求解或过滤；
- 24-cubed 2/4-rank selector 在同一 executable 内用 `MPI_COMM_SELF` 跑相同 product reference，
  复用既有 Taylor--Green `5e-12` 字段/结果分解阈值。

channel、constant IBM、material IBM 与 ideal IBM producer 只压缩已认证的 accepted report。
`nu_t_l2` 统一为全局体积加权 L2；IBM force 使用 Task 11
`consistency.total_N` 范数；material closure 使用 authenticated transport normalized residual；
ideal closure 使用 authenticated rho/rho-h remap、enthalpy-temperature 与 EOS 误差。正式
material/ideal selector 只走 IBM+WALE，且在既有 rollback/retry/closure mutation 全部通过后
才输出一行 pass row。

## selector 与正式 inventory

全部正式 executable 只接受冻结的 `formal ...` argv，未知组合返回 2。静态 CTest inventory
共 21 行：

- TGV convergence 12/24/48 单 rank；TGV 24-cubed 1/2/4；
- body-fitted channel 48-cubed 单 rank；
- constant IBM+WALE 48-cubed 单 rank、24-cubed 1/2/4；
- material 与 ideal IBM+WALE 各自 12/24-cubed 1/2/4。

24-cubed 行全部设置 exact `PROCESSORS`、M lock、7200 s；48-cubed行设置 H lock、43200 s；
12-cubed variable-density 行设置 M lock、1800 s。未注册 warped、prism 或 96-cubed。
`ctest -N -V -R '_formal$'` 只读确认全部 argv；正式行在 S3-S1 开发期均未执行。

## TDD、RED 与 mutation 证据

初始 RED 先于实现建立：scientific row contract 在 `validate_scientific_row` 失败；TGV 12-cubed
smoke target 可编译但 selector 返回 2。产品 smoke 首次进入数值路径后，以
`final_momentum_residual` 拒绝严格零 z 分量；诊断显示 x/y 约 `5e-17`、z 为
`2.06714e-05`。加入物理等价的均匀 spanwise 速度后，同一 12-cubed product smoke GREEN，
没有改产品门禁或阈值。

单元 mutation 覆盖 mismatched final time、point-sample restriction、一个 segment 被平均
slope 隐藏、zero/epsilon error、wrong ranks/process grid、stale density revision、第二次
WALE evaluation、zero identity、缺失 availability bit，以及 partial optional `nu_t` layout。
每项均由同一个 executable contract 拒绝。

## GREEN 与主 agent review

最终 `build/debug` 六个 S1 target 全部编译链接成功。按计划只执行：

```text
test_stage3_scientific_row_contract       PASS 0.01 s
test_wale_taylor_green_12_smoke_1_rank    PASS 0.81 s
2/2 PASS, 0 failures, real 0.83 s
```

主 agent 完整审查 schema/parser、cell-average restriction、gauge normalization、体积加权
`nu_t`、两段 observed order、24-cubed decomposition、selector exclusivity、row 唯一性、
timeout/lock/PROCESSORS 和全部 S1 diff。`git diff --check`、staged diff check、Clang 15 新文件
全格式检查、既有文件 changed-lines 格式检查与编译器 warning policy均通过。

未运行正式 12-cubed variable-density、24/48-cubed、96-cubed、Release、ASan、UBSan 或
long matrix；没有私有源码、研究数据、研究进程、push、publish 或 Stage 4--6 行为。

## 证据身份

- implementation/test staged diff SHA-256（治理文件加入前）：
  `0b99f26b457aaebe3f381b4b7ae409e53aa68518fafc7eb750b1f9eff96ceda3`；
- Debug scientific-row contract：
  `f6abda3e2ce6e3a85c27792fd9d823d030234adaf5e7757e416dd3b3fda878bf`；
- Debug TGV producer：
  `67128adf0043094293d783b3b8541ebbcc860e48f613e08d2c8ad6b8e008ec88`；
- Debug body/channel、constant、material、ideal producer：
  `c9bac845e2778acbfe66ad752ee67ac07c4836284a631ee3dcacd4b25b97dce6`、
  `4478b81a472c79bdc92790a40fecb4e80fd8c8c98c5577fab6467b9c00a14c53`、
  `376564177d1f78f956f2d214add50e9ced4ef282f7b715ec2cc7b3739b416a6a`、
  `cb1019068c9e466e01f0898e8636c6828eeb0783f3adf9b18244088eb11be967`；
- final `LastTest.log`：
  `5f6e5416f0bb07aab62adf0f9a8e1cf8db5e05737bed526d2a235e6aac76415a`。

激活候选 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

提交 subject：`test: freeze Stage 3 scientific selectors`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`
