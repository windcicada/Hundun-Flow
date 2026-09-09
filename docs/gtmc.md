# GTMC 主线开发记录

输入来自 `/home/wyf/code_dev/gtmc-hundun-handoff-20260909`。本地工作区为
`/home/wyf/code_dev/flow`，算例为 `/home/wyf/code_dev/cases/g`，证据为
`/home/wyf/code_dev/check`。算例沿用原网格、NASA 表、入口流量和数值门槛。

## Standards

G28 的 12 个提交具备 DCO，明确规范违规 0 项。G29 提升条件 1 项：
删除 `DEBUG-G29` 临时探针。设计建议 1 项：合并多处组分扩散系数算法。
后者属于维护性建议，独立修复的接纳以方程与验证结果为依据。

## Spec

审查发现 3 项缺口：原网格组分猜测收敛、非局部连续性 Jacobian、临时
探针清理。完整交接材料保存了原始失败证据。原网格单步、连续 short100、
medium3000 和参考焓偏移对照需由相同主线版本重新运行。

审查共计：Standards 1 项提升条件，Spec 3 项缺口；各轴的主要问题分别为
临时探针清理和原网格求解／Jacobian 完整性。

## 融合范围

基线 `f8b70b7` 与 G28 `0eb23ce` 的共同祖先为 `86542bb`。
融合保留 main 的反应、ESF、TCR、喷雾、源项及 Restart V4/V5 接线。
入口扩展在原 case wire 外包一层 version 19，内部继续携带反应／喷雾标志；
历史签名同时包含入口、热材料和原有模型组件。

G29 独立修复首先移植 NASA 反演 Newton polish，以及分配异常的集体状态
传播。NASA 反例在 G28 融合库上返回 dT=0、drho=0；隔离替换反演实现后
正负四个扰动全部通过。生产历史增加 `thermal-inverse-newton-polish-v1`。
组分猜测、AA3、总能量、面重建及其他 G29 组继续按依赖和反例审查接入。

初轮 250 项回归中，变物性标量 MPI 1/2/4 复现交接中的失败。
应用测试还保留“源码树包含任何 case.json 即失败”的旧假设，与现有反应／
喷雾夹具冲突；改为断言应用运行前后源码算例集合相同。
分配审计中的主动中止记为中止，后续独立复跑结果另记。

NASA 修复与融合后的 63 项检查全部通过，覆盖原失败的变物性标量 MPI
1/2/4、应用初始化／重启、case wire、原生反应／PaSR／ESF／TCR／喷雾。
本机 OpenMPI 的 unsigned MIN/MAX 独立合约在 4 ranks 通过。
初轮完整结果为 245/250；三项标量及应用测试已经复跑通过，一项主动中止的
immersed MPI2 分配扫描仍需完整复跑。
