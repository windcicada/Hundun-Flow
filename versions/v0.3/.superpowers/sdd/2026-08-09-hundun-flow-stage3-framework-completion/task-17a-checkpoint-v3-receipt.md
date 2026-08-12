# Task 17A：Checkpoint v3 IBM-only 实施回执

基线：`17b8434413a1dcc1c88070351f9669c6f6f7ed83`

状态：`TASK_GATE_ACCEPTED`

## 已实现边界

- constant-density、fixed-time、static LFP-GCIBM、LES none；
- per-rank accepted FlowState payload、统一 manifest、`COMPLETED` 最后发布；
- history/committed FlowState、accepted metadata、fixed control state；
- history/committed pressure authority，pending/force/cache 不持久化；
- exact inventory、CRC、presence、partition、geometry/plan/config identity；
- 1/2-rank collective write/read 和确定性失败收敛；
- state 与 facade 均准备完成后才发布，failed read 保持两者不变；
- `hundun` driver 的到期写入、相同分区读取和 BDF2 continuation。

## RED 与 mutation

- public header 缺失；
- manifest codec 缺失；
- `read_checkpoint_v3` 未定义；
- 1-rank accepted-state/authority continuation；
- 2-rank write/read collective；
- manifest CRC corruption rollback；
- viscosity identity mutation；
- 无本地 wall link 的合法 rank authority；
- 非法或不匹配的 BE/BDF2 metadata order；
- driver 未写 checkpoint、driver continuation identity 误含 I/O 编排字段。

## focused 证据

- `test_checkpoint_v3_header_contract`：PASS；
- `test_checkpoint_v3_codec`：PASS；
- `test_checkpoint_v3_1_rank`：PASS；
- `test_checkpoint_v3_2_rank`：PASS；
- `test_checkpoint_v3_collective_preflight_2_rank`：PASS；
- `test_checkpoint_v3_path_agreement_2_rank`：PASS；
- `test_checkpoint_v3_invalid_state_1_rank`：PASS；
- `test_checkpoint_v3_active_attempt_1_rank`：PASS；
- `test_task19a_immersed_flow_dispatch_1_rank`：PASS；
- `test_task19a_immersed_flow_dispatch_2_rank`：PASS。

clang/libc++ Release tests-off：PASS。

独立只读审查由 `agent_type="luna_worker"` 完成；原始 turn context 已核验为
`model=gpt-5.6-luna`、`reasoning_effort=max`。复审结论为无尚存
Critical/Important 阻断项。完整 diff、科学边界与最终 task 接受由主 agent
负责。

未运行 sanitizer、48³或完整矩阵；4-rank continuation 按计划进入 MVP gate。

## 版权独立性

仅复用本项目 Apache-2.0 的 Checkpoint v2 codec、CRC、文件发布和 collective
原语，以及 Task 11 自有 pressure-authority 数据结构。未复制外部项目源码，
未访问私有参考源码或研究数据，未新增运行时依赖。
