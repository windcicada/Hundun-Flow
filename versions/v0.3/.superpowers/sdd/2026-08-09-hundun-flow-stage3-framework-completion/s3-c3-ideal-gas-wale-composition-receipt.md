# S3-C3 ideal-gas body-fitted/IBM+WALE 与 Gate 5 验收回执

状态：`ACCEPTED`

accepted parent：`da9106a7d8ab54e513fc81b951f0d8336061895d`

accepted at：`2026-08-09T18:10:19+08:00`

## 实现边界

本 task 完成 schema-v3 profile 8（body-fitted/WALE/ideal-gas）与 profile 9
（LFP-GCIBM/WALE/ideal-gas），并关闭三种 density model 的 combined Gate 5。body-fitted
ideal-gas 只调用 C2 冻结的 `DensityClosureBridge::attempt_with_optional_wale`，把 D2 既有
`DensityClosureHooks` 传入同一个 `FixedStepMaterialDensityFlow` core；IBM ideal-gas 继续调用
D2 的 `ImmersedIdealGasDensityAdapter` 与既有 immersed PISO。没有新增 public attempt、第二
density field、第二 closure bridge、第二 pressure solver 或第三 corrector。

body-fitted closure collaborator 在 facade 构造期做 topology/geometry/boundary/MPI/registry/
field/active-layout 匹配；因此 `flow_ideal_gas_closure.cpp` 的白名单外小改动只实现既有 friend
adapter 的 body-fitted match，不增加 public API 或数值行为。driver 对 body-fitted ideal
构造 public all-active closure，对 IBM ideal 构造 active-domain closure；缺失 physics/WALE/
IBM setup 仍在 time step 前拒绝。

constant body facade 现在只在成功 commit 后转发真实 WALE evaluation count 与 coefficient
identity。material/ideal body 继续只在成功 commit 后发布 WALE summary；IBM 同时只在成功
commit 后发布 force。所有失败保持 force/WALE 未发布、FlowState 与 closure authority 逐位
回滚。

## exact conditional trace

`test_immersed_combined_retry` 对 profiles 2/3/5/6/8/9 逐项记录由真实认证权威支持的 trace：

```text
begin
density-predict-if-variable
closure-predict-if-ideal
lagged-gradient
wale-evaluate-once
momentum-predict
pressure-wall-authority-1-if-immersed
pressure-corrector-1
provisional-density-transport-if-variable
closure-provisional-if-ideal
pressure-wall-authority-2-if-immersed
pressure-corrector-2
final-density-transport-if-variable
closure-final-if-ideal
post-closure-assessment-if-ideal
final-residual
force-if-immersed
prepare-flow-state
prepare-density-transport-if-variable
prepare-closure-if-ideal
prepare-pressure-authority-if-immersed
collective-ready
publish-noexcept
```

每个 `*-if-*` event 只在对应 authenticated material/closure/pressure/force authority 存在时
记录；未激活模块必须 absent。证据同时要求 corrector count 精确为 2、ideal closure
evaluation count 精确为 3、WALE count 精确为 1、IBM pressure apply schedule 为 `1234`、
material provenance 为 `final_corrected`、commit metadata 已发布且 attempt inactive。

测试内定向 mutation 覆盖 third corrector、第二次/final WALE refresh、跳过 provisional
closure、provisional flux 冒充 final、collective-ready 后 fallible action，任何一项都会使
exact trace 不相等。`ExactPredictorResponse::Work` 与 `WaleSummary` 增加 no-throw copy
assignment 静态断言；collective-ready 后只执行已 prepare storage 的 `noexcept` publish、
平凡标量 authority 复制和 no-throw swap。

## TDD、RED 与 source mutation 证据

两个新 target 首先编译链接成功。有效初始 behavior RED：
`test_ideal_gas_wale_composition_1_rank` 2.67 s 后 exit 8，输出
`STAGE3_C3_UNSUPPORTED model=2 immersed=0 reason=1 rank=0`；失败来自 profile 8 的稳定
`invalid_input`，不是 missing symbol 或测试夹具。

每条 legal path 都注入 rank-local failure，并断言 lowest failing rank、history/committed/
trial/metadata bitwise rollback、同输入 retry 的 committed fields 与 WALE identity、失败
force/report isolation 及最终 commit。profile 9 另在 `before_commit` 注入 rank-local failure，
要求 authenticated non-retryable rejection、state/closure bitwise rollback 和无 force/WALE，
证明失败点位于 collective-ready 之前。

主 agent 用 `apply_patch` 施加并立即恢复两项可编译产品 mutation：

- 跳过 IBM provisional closure：1-rank CTest 10.62 s 后 exit 8，精确报
  `ideal-gas closure report is unavailable`；
- final transport 错用 provisional provenance：1-rank CTest 11.49 s 后 exit 8，精确报
  `material report bridge could not authenticate`。

最终源码已恢复，`finalize_from_corrector_two_flux` 明确获取
`MaterialFluxProvenance::final_corrected`，provisional closure 调用存在且成功报告的
evaluation count 为 3。C2 的一次性 WALE token/identity mutation 证据继续覆盖绕过标准
counter 的 final refresh。

## GREEN / Gate 5

最终 Debug Gate 5：`6/6 PASS`，0 failures，real 181.24 s：

- ideal-gas WALE composition 1/2-rank：29.75 / 17.38 s；
- six-profile combined retry 1/2-rank：69.01 / 40.45 s；
- schema-v3 fast acceptance 1/2-rank：15.04 / 9.60 s。

fast acceptance 最大 12-cubed、1 step；新增两个 8-cubed case 分别要求 profile 8 的
exactly-two correctors、WALE identity、无伪造 force，以及 profile 9 的 exactly-two
correctors、完整 force 与 WALE identity。原 profile 7 ideal IBM/no-WALE 仍在同一 gate 通过。

C1/C2/D2/header focused regression：`10/10 PASS`，0 failures，real 156.67 s；包括
constant IBM+WALE 1/2-rank、material body/IBM+WALE 1/2-rank、ideal IBM direct/transaction
1/2-rank 与两项 header contract。P0 registration/mutation contract `2/2 PASS`，real 0.66 s。
两个 C3 executable 对未知 selector 均精确返回 2 且无 stdout/stderr。

未运行正式 24/48³、96³、1/2/4 full matrix、Release、ASan、UBSan 或 long wake
statistics。环境无 `clang-format`；`bash -n`、`git diff --check`、staged diff check 与编译器
warning policy均通过。

## 调用图与单一 authority

CodeGraphF `0.9.8`；大函数 caller under-report 以 exact `rg` 交叉验证：

- `attempt_immersed` 仍恰有两处 `correct_active_pressure`、一处
  `prepare_wale_authority`、一处 provisional closure 与一处 final closure；
- body profile 8 只有一个 `attempt_with_optional_wale` caller，复用既有 material core 的
  predictor/provisional/final transport 与两次 PISO；
- IBM finalizer 只从 corrector #2 的 `final_corrected` flux 获取 authority；
- WALE 在 lagged accepted velocity 与 predictor/closure 后的 `rho_attempt` 上每 attempt
  求值一次，同一 authority 贯穿 momentum、wall viscosity、force 和 final checks；
- commit collective pass 后没有 MPI、allocation 或可失败动作。

从 C1 parent 到 C3 candidate 的完整 combined diff 已由主 agent 复审。没有改动 Task 11
force sign、scientific thresholds、exactly-two-PISO、rollback、Restart 或 MPI 最低失败 rank
要求。

## 证据身份

- implementation/test staged diff SHA-256（治理文件加入前）：
  `b1d47fe79ce03b14e872c0d68352212f04ca5c2ecf9620ffda0db004c9500bb5`；
- Debug `test_ideal_gas_wale_composition` SHA-256：
  `e3d62740cab077674e5a84d55d5b6ea71aa1257e538a49a641b2b1bd13d0d804`；
- Debug `test_immersed_combined_retry` SHA-256：
  `5d65a631e34adc44a44c6bc8a8321f3a9bf04856c16f2dceaf074219f64416b2`；
- Debug `hundun` SHA-256：
  `bf627a44e4952522ef8b70a3bc7a01dc5f7d9d42e1e79c889feddff395789374`；
- final focused regression `LastTest.log` SHA-256：
  `a1ae09fd22ea8619485264eb068c843dc5b83734f1d1dd25fa294a3b52518376`；
- CMake `3.31.12`；CodeGraphF `0.9.8`。

激活候选 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

## 主 agent review

- profile 8 进入既有 material PISO + closure hooks；profile 9 进入既有 immersed PISO；
  两者没有第二套 pressure/final-flux/closure authority；
- active-domain ideal closure 与 body all-active closure 在 facade 构造前区分，cp/R/p0 在
  attempt 时 bitwise 匹配；
- final WALE summary、force、material/ideal reports 与 pressure authority 都只属于成功
  attempt；失败不泄漏 persistent attempt report；
- R1 worker 结果仍在 `coast/stage3-infrastructure-lane` 隔离 dirty worktree，未集成；O1
  未启动；
- 无私有源码、研究数据、研究进程、push、publish 或 Stage 4–6 实现。

提交 subject：`feat: complete immersed density and WALE composition`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`
