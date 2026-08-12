# S3-R1 Checkpoint constant profiles 验收回执

状态：`ACCEPTED`

accepted at：`2026-08-09T19:21:04+08:00`

## 委派与集成身份

- worker worktree：`/home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure`；
- worker baseline HEAD：`3735742b65d8372d54d6ac517eac932281f37fe2`；
- signed handoff commit：`c2a3b71edf1fe5267e9a66642e191d4453bacb04`；
- main integration parent：`cee4c60d46816e3b4772ca485c2fbcd63efaa6cb`；
- main integration commit：`efcdd25c80639816f14f25f1fa933adca1ee93a9`；
- handoff staged diff SHA-256：
  `f2f0c77c296c45993cb3455f4f4edf84822fe9cc948b324fccbaf3be0a714e0d`。

worker 严格只修改 R1 packet 的七个路径，不暂存、不提交。主 agent 在隔离 worktree
完成完整 byte protocol、调用方、事务、MPI、失败中性与版权独立性审查，独立重跑
GREEN 后创建 DCO handoff，再按 R1→O1 顺序集成到 main。

## 协议与实现边界

Checkpoint v3 presence 追加并固定为 profile 2 `constant_body_fitted_wale` 与 profile 3
`constant_static_ibm_wale`，旧 profile 1 的 legacy encoder branch、manifest schema、endian
和旧 public overload 均保留。冻结 profile-1 fixture 仍精确为 558 bytes，CRC-64/ECMA-182
为 `8909770059348032994`。

profile 2 要求 canonical absent IBM authority；profile 3 要求恰一个 IBM identity section；
两者都要求恰一个 WALE section。WALE section 为 47 bytes，只保存三个 controls、既有
17A 规范化 resolved numerical-config CRC、transient schema/version 与
`nu_t/mu_sgs/mu_eff` 字段身份。不保存 transient 数组、summary、attempt identity、地址或
第二次 WALE evaluation。

主 agent 审查发现原 39-byte 候选没有给 profile 2 绑定 resolved numerical config：相同
layout/WALE controls 但不同合法黏度会被错误恢复。worker 先加入可执行 mutation，证明
profile 2 以 2× 动态黏度读取时错误 restored；随后原样抽取 17A normalization CRC helper，
将同一 CRC 加入 profile 2/3 WALE identity 并在 read publish 前精确比较。WALE controls 的
验证同时与 `WaleModel` inclusive ranges 对齐：`Cw [1e-6,1]`、`Prt/Sct [0.1,10]`。

read 仍先完成 marker/manifest、partition、identity、rank CRC、payload 和 collective
readiness，再准备 replacement 并 publish。错误配置和损坏 manifest 都要求
`rollback=passed`、FlowState bitwise 不变、last WALE identity 不变；随后原配置仍可成功
读取并在下一 accepted attempt 重算相同 WALE summary。write publish-last 顺序未改。

## TDD 与 GREEN 证据

- 初始 codec RED：header seam 编译，profile-2 manifest roundtrip 在
  `decoded.has_value()` 失败；
- 主审修复 RED：1-rank profile 2 在
  `!incompatible.restored()` 失败，证明不同合法黏度被错误接受；
- 修复后 worker focused：4/4 PASS；registration/layout：5/5 PASS；
- 主 agent 隔离 worktree 独立 focused：4/4 PASS，real 30.12 s；
- 集成 main 后 R1/O1 组合 gate：8/8 PASS，real 30.61 s；其中 R1 1-rank
  18.83 s、2-rank 11.70 s、codec/header 均 0.01 s；
- 集成 main registration/layout：5/5 PASS，real 0.13 s；
- 集成 main Clang 15/libc++ Release tests-off `hundun`：PASS。

Debug binary SHA-256：

- `test_checkpoint_v3_codec`：
  `89e4d636b1611fa4d8516ce14edf3fec77deeb9e1b56a838727aa646271c5693`；
- `test_checkpoint_v3_wale`：
  `b941e8892b8bd6a4c20e47b6f60f475f58d72ac56b57ae09e53fe7be8fcb0755`；
- `test_checkpoint_v3_header_contract`：
  `885963ee8b6304aa5b5ab2225247828985221df58f26e123af72731782cac4dd`。

未运行 24/48/96-cubed、sanitizer、rank-changing Restart 或正式矩阵；未访问私有源码、
研究数据或研究进程，未 push/publish。Task 11 阈值、两次 PISO、force sign、rollback、
Restart 与 MPI 一致性 authority 均未修改。

## 权威与 DCO

激活 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

提交 subject：`feat: extend Checkpoint v3 for constant WALE`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`
