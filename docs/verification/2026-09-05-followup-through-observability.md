# 重启、receipt与性能观测后续

起点：`0b0eff11c8365e570725d6220d031c84a6aef4e0`，
分支`codex/v04-restart-receipt-observability`。用户要求推进到剩余回归和性能观测，
然后给出优化方案；不实施内核优化、不推送、不长测、不修改原始checkpoint或历史凭据。
沿用上一轮Release/独立ASan+UBSan工具链。已有AGENTS.md和.codegraphf不纳入提交。

## 1. 重启兼容

- 默认严格读取原69d8eee小算例：10307，restart_load，0次推进。
- 最初只重建schema/product尾部的兼容实现能通过合成测试，但原始文件仍失败。
  独立文件重放发现BoundaryCompiler的semantic_hash也包含registry fingerprint；
  现复用同一哈希实现重建旧边界身份，不分配第二套边界/求解工作区。
- 兼容必须显式选择mg-bundle-ghost-v1；只接受完整格式2历史，默认行为不变。
- 公共合成测试验证旧身份拒绝/显式接受、错误身份拒绝、全部场与速率历史、面通量、
  时间控制、压力参考保持一致，以及新checkpoint使用当前身份并通过严格加载。
- 1/2/4 ranks：3/3通过。实际原二进制checkpoint在2 ranks从步1续到步2，
  随后新checkpoint默认严格路径4 ranks续到步3。
- thin runner的完成标记和accumulator仍核对来源身份，输出标明来源和迁移，
  新文件使用当前身份。尚需目标Re3900短诊断验证，不宣称已经迁移目标文件。

日志前缀`/tmp/hundun-followup-`；最终汇总时打包哈希与原始结果。

## 2. receipt

用户确认通过receipt-validate CLI验证。新增13项CLI检查，先确认旧程序缺失
历史模式，再实现显式`--historical`、根目录映射及内容寻址附件库。
损坏、缺失、重复映射、不完整却要求complete均拒绝；默认严格绑定当前工具不变。

- R4原文件不改写，SHA-256仍为
  `1dfdeb16b92f4a29a5dc2a0c4daa146499e8decb1af5df197b652203e36b2e3f`。
- 两份历史工具从59985db精确恢复为`.snapshot`附件，其哈希与R4一致。
  历史附件不执行、不补水印；当前生产工具水印保留。
- 原R4的14项附件通过历史校验，当前校验器与原提取器身份分别报告。
- 新生成[R5凭据](v0.4-literature-data-receipt-r5-partial.json)，绑定17项附件，
  包括父R4、历史提取工具、当前工具；默认严格校验通过。
  这次没有重新提取实验数据，也没有改变参考数值或完整性判定。
- R4/R5均为`complete=false`；R5的`--require-complete`按预期退出2。
  Norberg未完成项没有因路径修复而消失。
- [复现命令及历史快照说明](literature-r4-artifacts/README.md)。
  原始PDF/数字化附件仍是外部依赖；不将缺少附件的环境伪报为通过。

## 3. 剩余回归和观测

已补充并验证：

- `initialize_restart`真实new分配失败扫描：Cartesian/IBM、1/2/4 ranks，6/6。
  每次失败后在同一driver上重试，核对历史并检查跟踪的C++/通信资源归零。
- 新增`RestartReader::load`扫描时，原程序rank 0首次分配失败触发SIGABRT。
  根因包括read_current_name的noexcept、广播前本地分配未汇总、以及读入/重分区
  构造的外层catch。现按本地阶段汇总错误，保留原输出对象，并区分分配与文件错误。
  完整扫描Cartesian/IBM、1/2/4 ranks，6/6；2-rank示例扫描108/100个分配点。
- 真实温暖非均匀重启夹具：外推alpha=2被拒绝，alpha=1梯级候选被接受，
  整步一次完成；1/2/4 ranks均覆盖。与原候选计数及应用计时检查合计7/7。
- POSIX边界故障注入确认open/write/fsync的单次EINTR原先会终止输出。
  现重试EINTR、完成短写，但不重试close；ENOSPC等首个错误不会被close覆盖。
  Visit/Screen/Monitor/Evidence提供固定容量的errno、操作、路径、失败rank报告，
  即使仅rank 0请求报告也保持同一MPI流程。实际目录冲突在ApplicationService中
  报告为已接受步之后的Visit失败。相关新增及原分配回归24/24。
- ApplicationRunReport/CLI区分完整本地运行时间与原advance计时；计时包含初始化、
  输出、资源、Evidence和清理。数值候选细分8个不重叠的包含式阶段，Schur求解分为
  准备、Krylov、关闭/焓恢复。时间只在阶段边界采样，没有逐单元计时或额外候选MPI。
- thin runner可选`--observe-performance`保存各rank的整步及子阶段数据。
  其计时不包含写这份观测记录本身的gather/write，不能据此声称完全无I/O。
- 可选PMPI动态附件只在advance区间统计9种阻塞collective及调用点，不增加collective。
  独立2-rank检查确认区间内1次Allreduce、2次Bcast、1次Barrier被准确计数，区间外忽略。

待完成：最终Release/ASan+UBSan验证；目标Re3900从743起的单次3步诊断（含最终输出/
checkpoint，不是长测）；把观测原始文件与符号归类结果打包，再写优化优先级。
尚不宣称RestartWriter所有POSIX错误上下文/分配路径均已完成审计。
