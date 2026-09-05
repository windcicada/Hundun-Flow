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

待实施。历史R4文件不改写；两份历史工具脚本可从59985db精确恢复，
当前脚本相对它们仅新增水印。新校验区分历史工具身份与当前校验器身份，
使用显式路径映射；历史complete=false不因工程修复变成true。

## 3. 剩余回归和观测

待实施：重启全分配失败、外推拒绝后梯级接受、I/O错误信息，
完整步墙钟/模块耗时及通信归类；随后提出可归因优化方案。
