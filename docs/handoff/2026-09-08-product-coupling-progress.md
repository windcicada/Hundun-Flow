# 产品耦合接线进度

用户要求以实际产品接线验收为合并完成条件，然后才能续算圆柱绕流。
`f4f04ec` / `de0e77a` 仅归集源码与验证单相兼容，不代表完整耦合验收。
圆柱服务 `hundun-re3900-integrated-long-20260908.service` 保持 SIGSTOP；
暂停时最后完整健康记录为 10577，续算源仍为独立封存的 10510 checkpoint。

## 已接通的产品路径

| 路径 | 产品行为 |
| --- | --- |
| 显式源项准入 | `ContributionAdmission` 绑定模型/气相身份；编译计划和实际方程消费同时校验 capability 与 source identity。默认无准入，旧惰性计划指纹保留。 |
| finite-rate mean | 通过完整 `(p,h,Y)` 查询得到净组分质量生成率，按名称映射到 Ns−1 产品字段；核查机理、组分、元素、分子量、焓参考、修订号、输出缓冲和原生 EOS/热容一致性。源项进入原生 BDF/EX2 非对流速率历史。 |
| PaSR | 使用同一已核查净反应率；`reactant_depletion_l1_v1`、单元体积立方根滤宽、气相质量分数加权分子扩散率和实际 `mu_eff-mu` 接入既有混合时间及反应分数函数。显式配置 C_Z 与 Sc_t。 |
| 压力/能量 | 组分改变后由现有两次 PISO / `rho*h-p` 联立闭合响应；总热化学焓不重复加入形成焓放热。 |
| 初始与失败状态 | 均匀反应初值也构造首步化学速率，修复只在 IBM 初值重建速率导致首步漏反应的问题。单 rank provider 失败或输出合同损坏使所有 rank 拒绝；已接受物理字段、两层历史和通量不变。 |
| 输入与恢复 | 严格 reaction JSON、跨 rank wire、case/模型指纹与方法历史签名；现有 Restart V3 保存反应速率历史。载入值精确一致，带放热续算按求解器无量纲终止容差比较；惰性同分异构试验连续/恢复值精确相同。 |
| 实际后端 | `HUNDUN_ENABLE_REACTING_CANTERA=ON` 将已验证 Cantera 3.2.0 包接入产品；关闭时 direct_cantera 明确拒绝。analytic_isomer 是具名合成验收表示，不是燃料机理。 |

当前 finite-rate / PaSR 使用原生显式 BDF/EX2 反应率历史。
Cantera 区间积分器在后续 ESF 接线中消费；此节点不宣称已接通 ESF 两个化学半区间。

## 当前验证

- Clang/libc++ focused 25/25 PASS：产品冻结、原生方程、压力能量重试、历史恢复、
  反应准入与 finite-rate/PaSR 的 1/2/4-rank 产品测试。
- GCC11/Cantera 产品矩阵 12/12 PASS：外部 provider、具名解析表示、PaSR、真实 Cantera，
  各 1/2/4 ranks；包含带放热的密度与 `rho*h-p` 守恒、实际磁盘 Restart 及续算比较。
- 两个真实 CLI 两步计算完成：Clang 解析表示与 GCC11 Cantera；两组 runtime evidence
  validator 均返回 0。Cantera 机理为原创 A→B 合成夹具，其 SHA 随 case.json 冻结。
- GCC 编译确认并修复了既有代码依赖间接包含 `<algorithm>` / `<climits>` 的问题。
- 日志和运行目录在自有工作树的 `build-coupling-evidence/`、`build-integrated-cantera/`；
  未覆盖封存的圆柱数据，也未运行新的圆柱长算。

## 完整合并前仍须完成

| 待接线 | 必须验收的产品合同 |
| --- | --- |
| ESF / TCR | 实际空间输运/梯度/Halo、持续随机场、两个化学半区间、TCR 分支/折点历史、共同接受/拒绝、完整重启。科学证据不足的 validated 模式保持明确拒绝。 |
| 喷雾 P4–P8 | 实际气相采样与几何/IBM 事件、液体资产身份、parcel 生命周期/TAB/注射器、跨 rank 迁移与源路由、质量源进入连续性和压力闭合、动量/总焓/动能收支。 |
| 共同事务和运行切换 | 全部模型状态与气相一次提交及回退；完整产品 Restart 和真实 CLI 验收后，再重建冻结程序并从封存 10510 做短续算验收，最后恢复圆柱长算。 |

本节点是接线中的可验证增量，完整合并尚未完成。
