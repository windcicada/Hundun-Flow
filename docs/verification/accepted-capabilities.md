# 已接受能力

本页区分“已有数值产品证据”和“V1.0 发布接受”。V04-2 已接受 pressure--enthalpy numerical product closure；V1.0 只有在 [`HUNDUN_V1_RE3900_MEDIUM_RELEASE_POLICY_V1`](v1.0-re3900-medium-release-policy.json) 的全部门通过后才可发布。

已有 focused、MPI 和真实产品路径证据支持：

- C++17、MPI 3、CPU tests-off 产品构建；
- tensor-stretched Cartesian 网格、持久 halo 与原生多重网格；
- 单相低马赫理想气体热力学和守恒方程；
- 两次 PISO 与 same-target continuity--energy pressure--enthalpy refinement；
- `p/h/rho/T/U`、EOS 和最终面质量通量的同步事务发布；
- EOS、continuity、energy、closed mass、gauge 五个独立终端门；
- common-face owner AFC 的 MPI partition-face 守恒/唯一 alpha；
- provisional/committed CFL typed certificate，fixed fatal 与 adaptive retry/rollback；
- 静止封闭 STL IBM 的严格二次及可审计 adaptive-order stencil；
- WALE、Vreman wall-function、压力/黏性分解的最终表面力；
- exact-history rank-changing Restart、Visit、Evidence V6 immutable candidate identity；
- 1/2/4-rank decomposition、MMS/temporal、rollback、I/O 和 CLI focused 覆盖。

这些证据不自动证明长程统计收敛、所有网格/时间步独立、任意燃烧室几何鲁棒或旧 v0.4 规格中的 full20/420-cycle release gate。V1.0 的最终状态以外部封存的 exact candidate receipt、Re=3900 medium evidence 和 GitHub tag 为准。
