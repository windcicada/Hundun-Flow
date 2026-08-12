# S3-A1 application combinations 验收回执

状态：`ACCEPTED`

accepted at：`2026-08-09T21:59:06+08:00`

accepted parent：`e2b8d1686758f47b23072d3b47e97e1b16e7db15`

accepted commit：`3c1b8dc30445be74c81b1f1b6cc09b97788fd221`

## Same-executable matrix 与 construction order

单一 `hundun` executable 现已覆盖 Checkpoint presence 1--9 的精确合法矩阵：
constant/material/ideal-gas × static IBM/body-fitted × optional WALE。driver 保持冻结构造顺序：
config broadcast 后构造 mesh/boundary，随后只在 presence 要求时构造 IBM static plans、
ideal-gas closure 和 WALE，最后构造 `FixedStepImmersedFlow`；Checkpoint restore 位于所有
optional numerical object 完成后、retry loop 前。

`stage3_dispatch_inventory_complete` 以实际 `CheckpointV3WriteModules` object view 为唯一
应用层完整性证据，并在 restore、CASE stdout、diagnostic output 和 step 前 collective
收敛。profiles 1--9 精确接受；复制 profile 6 后清空 IBM surface、任一 incomplete setup、
unknown presence 10 均被 pre-step 拒绝。不存在 runtime-generated fallback row。

诊断 session 仅在未来存在 scheduled step 时构造；IBM static summary 也只在该条件下构造。
accepted commit 后按 presence inventory 发布 kinds 18--22：IBM 仅 18--21，body-fitted WALE
仅 22，IBM+WALE 为 18--22。static/LFP/wall-force 使用 O2 collective provider，WALE 原样
复用 O1 summary；failed/retried attempt 不发布记录，Checkpoint 仍只在 commit 后发布。

## Main-owned O2 integration repair

真实 material/ideal driver 接线发现 O2 `diag_immersed_module.cpp` 的 base report extractor
硬编码 `std::get<StepAttemptReport>`，导致非 constant accepted variant 被错误归类为
wall-force unavailable。主 agent仅将该 read-only extractor 扩展为三种既有 authenticated
variant：constant 取自身、material 取 `flow()`、ideal 取 `flow().flow()`。没有修改 report、
force sample、seal、数值计算、阈值或 provider schema。A1 9-profile acceptance 覆盖该修复。

## RED、GREEN 与 product projection

- construction RED：`test_stage3_dispatch_inventory` 首次完整编译后精确因产品 validator
  未定义而链接失败；实现后 legal 1--9 与 row-6/ID-10/incomplete mutations 全部通过；
- 首次 diagnostics RED：static provider 拒绝非 construction frame；接线改为冻结的 step 0/
  time 0 static frame；
- 第二次 diagnostics RED：material profile 命中上述 O2 variant extractor 缺陷；scope repair
  后 constant/material/ideal profiles 全部通过；
- final combined gate：10/10 PASS，real 205.13 s，包含 P0 registration/mutation、dispatch、
  O2 provider 1/2-rank、9-profile flow 1/2-rank 与 restart 1/2-rank；
- 最终 include dependency 修复后，受影响 provider/dispatch gate：4/4 PASS，real 16.49 s；
- 9-profile flow：1-rank 25.28 s、2-rank 16.06 s；
- 9-profile continuation：1-rank 91.13 s、2-rank 56.05 s；每个 profile 的 continuous 与
  split-restart step-2 Checkpoint 逐文件相等，diagnostic provider 顺序与 state fingerprint
  相等；attempt provenance revision 保留各自真实身份，不被伪造为 byte-identical；
- Clang 15/libc++ Release tests-off `hundun`：PASS；`nm -C` 无 `TestAccess` 或
  `ENABLE_TEST_ACCESS`；`ldd` 仅系统 C/C++、MPI、pthread/dl/hwloc/numa 等系统依赖，无
  Python、vendor solver 或 GPU runtime。

内容 SHA-256：

- `test_stage3_dispatch_inventory.cpp`：
  `d3b11ec190d4f092fd7003871b178521ec5a4fecd36dced54c4f30f43d0e564e`；
- `stage3_flow_models_fast.sh`：
  `600c762efe9783a1aeac640cde9bdf7f43a195377709fa66960c7ab2c2d9ae25`；
- `stage3_restart_fast.sh`：
  `f32f598e3c6944d90ed9f460676a006a7e27be6746dd806ff3587d85f9814979`；
- tests-off `hundun`：
  `5345d80945780a5794aa98b0370fad89c523886e10646cf6374a48ce7832921d`。

未运行 24/48/96-cubed 或任何正式矩阵；未访问私有源码、研究数据或研究进程，未
push/publish。Task 11 科学 authority、阈值、两次 PISO、force sign、rollback、Restart
与 MPI 一致性要求均未修改。

## 权威与 DCO

激活 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

提交 subject：`feat: complete Stage 3 application combinations`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`
