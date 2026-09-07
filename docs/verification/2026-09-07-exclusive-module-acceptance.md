# 分模块修复：独占验收窗口

## 长测暂停边界

用户授权暂停长测、完成修改并验收后再恢复。2026-09-07 09:29 +0800
对 `hundun-re3900-scalar-repaired-20260906.service` 执行正常 service stop。
随后确认 `inactive/dead`、MainPID=0、无 mpirun、无冻结 runner 进程，独占锁可取得。
没有使用 SIGSTOP 后叠加另一组 CFD rank。

暂停前最后完整 health 行：step=5936，time=0.081969703246474651 s，
attempts=1、retry=0、BDF2，continuity=9.6084679016128186e-8、
energy=3.9564146455100481e-7、committed CFL(out)=0.34664011716605792。
这只是最后一条完整记录，不推断被中断的当前步已提交。

最新持久化恢复点为 `generation-5500-64323211746098`，
time=0.075949017495889531 s。5500 之后的非持久化推进不能恢复。
源 run 为 `long-scalar-repaired-35000-20260906`；源 checkpoint、
Visit、统计与日志保持原位只读，未来续算使用新 run 目录。

候选起点 `03087ed`，本地分支 `codex/v04-restart-receipt-observability`。
第一次构建沿用 `build-review-fixes`：Release、clang/libc++、
`-O3 -march=znver3 -mno-fma -ffp-contract=off`，不启用 fast-math，
构建最多 8 jobs，测试串行；尚未声明候选程序验收通过。

## 验收记录

远端只读核对：main 仍为 `0127616f3716b7d7d6fcb4a00831129c2d6631b6`。
本地起点完整 SHA 为 `03087ed158327e020542ebbba4ac26e6716851c9`。
旧报告 `2026-09-07-module-review.md` 保留为当时的待验收记录；本文件补充后续结果，
不将旧报告里的“未运行”倒写成当时已通过。

证据目录：
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/module-audit-20260907`。

### 标量容量、范围与回退

`product_scalar_contract_experiment.cpp` 扩展两组公共 ProductDriver 回归：

- capacity：4/5/12/60/61 个被动标量；允许的目录内再混排一个 EOS 组分；
  PISO/SIMPLE × FGMRES(restart=2)/BiCGStab；1/2/4 ranks。
  对负值和大于 1 的被动标量分别使用 -0.2、1.2；这不是组分质量分数。
- capacity-ranges：4/5/12 个被动标量的非均匀闭域输运，EOS 组分位于目录首、中、尾；
  两种耦合与两种 Krylov。独立计算 rho(q,h,p)、库存和范围；两步推进后检查常量、
  上下界和库存。另用单次压力迭代触发真实失败，逐项对比已接受全载荷没有变化。

此前 `2*N_passive` 的构建期归约容量修复已通过上述 1/2/4-rank 测试。
同时实际复现 fresh initialize 对合法有符号被动标量返回 5/10201；
`ProductDriver::initialize` 改为全部标量检查有限性，仅组分检查 [0,1] 及平衡组成。
这与源项/重启中的既有被动标量合同一致，未改变默认零标量 Re3900 数值路径。
最大目录受 IoServicePlan 的 64 字段上限约束，3 个主字段 + 61 个标量实测通过，
没有用增加 restart、减少合法标量数或热路径扩容绕过归约容量。

### 应用初场与方法历史

`ApplicationRunOptions` 显式接入 `restart_history_policy` 和可选 `initial_state`。
默认仍拒绝未知方法历史；初场与 restart 互斥，fresh evidence 从 t=0 开始。
`app_driver.cpp` 传递恢复策略至源预期与 ProductDriver，记录恢复首步 BE、后续 BDF2。
没有修改完整 image 的 `backward_euler_recovery`，也没有重定义闭域质量目标。

实际 RED→GREEN 包括：公共选项恢复策略、显式温度 350 K 的 checkpoint 焓、
歧义边界初场 initialize/10505 拒绝、有符号被动标量初始化。
CLI 新增 `--restart-method-recovery` 和 `--initial-state p,T,Ux,Uy,Uz[,q...]`；
解析非有限值、缺分量、重复选项、初场/重启混用及迁移策略混用均拒绝。
原 continuity witness 夹具的 400 K 初场现在显式给出；边界、压力最大迭代 1 和
预期 stage=44 的求解失败目标保持不变，不通过恢复 last-wins 来使测试通过。

方法签名组件及已知迁移表见 `v04-method-history-contract.md`。
非均匀闭域回归包含合法的“目标质量与库存差 1e-13”及毒化旧速率：
同方法精确续算 → 显式恢复 → 再次精确续算，逐阶段核对质量目标、速率和 BE/BDF2。
runner 的 1/2/4-rank 非零旧样本链已通过，恢复建立新 epoch，重新精确续算继承新样本，
旧 checkpoint/统计/附件哈希不变。没有把一个恢复步解释为充分发展。

### 唯一分配与峰值观测

新增六组 cartesian/immersed × 1/2/4 ranks 的标量内存 profile，包含五个被动标量和
一个混排 EOS 组分。逐阶段捕获 C++ 实际唯一分配，并分离 arena 与面通量的对齐分配；
覆盖 compile/create、冷/热 advance、Visit、restart write/read、owned image 存活期间
销毁/重建/恢复和再推进。所有六组释放后的 C++ 活跃字节/对象均为零，MPI 资源检查通过。

| 标量夹具 | ranks | 最大单 rank C++ 请求字节峰值 | 峰值阶段 |
|---|---:|---:|---|
| Cartesian 8³ | 1 | 11,786,222 | restart_read |
| Cartesian 8³ | 2 | 9,608,906 | restart_read |
| Cartesian 8³ | 4 | 8,244,226 | restart_read |
| IBM 16³ | 1 | 42,205,926 | restart_read |
| IBM 16³ | 2 | 28,199,290 | restart_read |
| IBM 16³ | 4 | 19,602,618 | restart_read |

证据 `scalar-memory-and-initial-red.log` 保存完整 profile；其中末尾初场测试的 RED 不隐去。
`scalar-allocation-sweeps.log` 记录 create/initialize_restart/read_restart 的 1/2/4 ranks
单 rank 分配失败逐位置扫描，共九项通过，包含 rank 0 与最后一个非零 rank。

这是测试夹具中的分配预算核对，不是新增生产硬内存上限。输入 model 的预先分配、
MPI/libc 内部 malloc、栈和 allocator 余量未纳入这些 C++ 请求字节。
`owned_payload_bytes` 仍是 remap 容量报告，不能称为整个进程 RSS 上限。
完整运行期硬预算仍需明确各库归属和跨阶段同时活跃集合，不能将此表外推到 Re3900。

### 扩展回归状态

首轮 `-L 'v04_focused|scalar-conservation' -j1`：247 项中 246 通过，耗时 1118.78 s。
唯一失败为新初场歧义拒绝挡住旧 CLI continuity witness，根因和修复如上。
该首轮不包含后来新增的 CLI 初场参数与 IBM 标量 MMS，不能替代其最终构建验收。
日志 `release-focused.log` 保留失败。后续 `expanded-contracts.log` 为 59 项、302.20 s，
其中 58 项通过；唯一失败为新增圆柱 MMS 夹具 sin(2π) 与 sin(0) 不精确相同，
形成未焊接的表面接缝，surface/2 正确拒绝。修正夹具复用同一顶点，不改表面校验。
随后 MMS 与 observer 两项均通过，17.44 s，见 `ibm-mms-observer.log`。
ordinary CLI 历史/初场 1/2/4 ranks、runner 统计链、15 项新增标量分配失败逐位置扫描、
六项内存 profile、27 项标量合同和非均匀闭域方法历史链均在这组日志中通过。

ASan+UBSan 检查构建：9/9 通过，304.65 s，覆盖 capacity/capacity-ranges 1/2/4 ranks、
IBM 公共接口、ApplicationService 初场/恢复及 ProductDriver 方法历史链。
日志 `sanitize-contracts.log`；`halt_on_error=1`，因 MPI 运行库进程退出分配使用
`detect_leaks=0`，因此不称为 LeakSanitizer 泄漏检查通过。正常销毁分配计数另有上述证据。

### IBM 标量空间精度：守恒通过不等于曲面二阶

`solver_ibm_scalar_mms_test.cpp` 使用真实公共 `cartesian_diffusion` 与
`correct_impermeable_scalar_diffusion` 组成 `(I-0.02 L)q=f`，独立解析体积平均 RHS，
测试侧 CG 解离散系统并核对真残差；不是更换生产压力求解器。
域为 [-1.5,1.5]³，N=24/48/96；内部流体分别为对齐方盒、旋转方盒、圆柱。
圆柱 STL 随网格细化到 4N 个周向分段，几何离散误差仍包含在测量中。

| 形状 | L2(N=24) | L2(N=48) | L2(N=96) | 两次观测阶数 |
|---|---:|---:|---:|---|
| 对齐平壁 | 0.0156581 | 0.00398562 | 0.00100084 | 1.974 / 1.994 |
| 旋转斜壁 | 0.0538830 | 0.0296574 | 0.0135142 | 0.861 / 1.134 |
| 圆柱壁 | 0.00549450 | 0.00174678 | 0.000961918 | 1.653 / 0.861 |

独立六邻居配对行与公共算子差最多 5.69e-13；封闭域净通量绝对值最多 4.94e-15；
常量误差为 0；线性真相对残差小于 9.59e-13。只有对齐平壁设置已知二阶门槛。
斜壁约一阶、圆柱未呈稳定二阶，是本次实测的方法限制，不能拿守恒通过覆盖它。
不撤销二值控制体零通量修复；若要求光滑曲面的二阶解，需要另行设计一致的几何/通量
离散并同时守恒、有界、EOS 与时间阶数验收，不能仅替换单个 ghost 插值。
当前零标量 Re3900 不使用这一标量扩散路径；其热传导 IBM 路径没有在这里改写。

新增界面支持域回归还检查 `interface_cells()` 无重复、覆盖每个 link 的 fluid owner；
从零数组出发的压力功与热扩散修正，在该集合外均为零。用于后续局部性能实验的前置证据。

### 配对残差与组成迭代核查

`mass_pairing_residual` 本身仍是诊断字段，未擅自新增经验接受门槛。
其上游约束位于 `pressure_energy_components_converged`：有标量时连续性目标收紧至
不高于 16 ε，配合 capture 的预测密度/通量与同一 a0，构成质量配对行和。
本次多标量/混排测试未得到违反库存或常量保持的反例；不能将“没有该单独门槛”写成
已发生守恒错误。`kRecouple` 在同一个 dt 上回退组成 iterate，与真正 dt retry 分开；
中间组成没有提交。缓存复用仍需证明 BDF、源项、limiter、边界和速率权威同时不变。

## 本轮唯一性能实验

选择候选能量残差的 IBM 修正累加，不改变压力算法。两个从零 scratch 出发的修正
只在 `interface_cells()` 中非零；旧实现各遍历一次完整局部域，实验改为遍历已有的
唯一界面单元列表。保留 scratch 全域清零、两个修正的原始算术顺序、solid 最终清零、
物性/边界/halo、候选审计及所有 consensus。没有增加缓存/全场数组，也没有删除 donor
计算或改为 link 累加而重复计入同一 owner。

代码为 `core_product_freeze.cpp::execute_attempt` 内的 `assemble_candidate_energy`。
67 项相关 Release 回归通过，137.75 s，见 `sparse-candidate-regressions.log`。
三个 1/2/4-rank 产品综合 sanitizer 回归通过；IBM 标量检查首次触及测试框架 180 s
时限（没有 sanitizer 报错），将检查构建的标量测试时限改为 600 s 后，212.27 s
完整通过四种 passive/EOS × constant/nonuniform 情况，见 `sanitize-scalar-ibm-complete.log`。
不把首次 timeout 隐藏或解释为生产数值失败；Release 测试时限和全部数值阈值不变。

独占物理窗口 5501–5510，128 ranks，两程序各一轮；冻结脚本、程序、build manifest、
case 和 spec 位于证据目录的 `preopt-frozen` / `postopt-frozen`。
同一 source generation-5500-64323211746098，case/spec 的 SHA256 完全一致。
新 V3 observer 对 preopt 窗口确认 complete=true、128 ranks、10 steps、61 loops；
10 次 C1 spatial + 51 次 C2 diagonal，无 C2 spatial，不用这组数据决定后备切换。
前测 `advance` rank 均值为 8.399503 s/步、候选为 1.973020 s/步、
候选 residual 为 0.542895 s/步；进程整段墙钟为 98.94 s。
后测 advance=8.416037 s/步、进程墙钟=99.59 s，无总耗时收益，因此已撤回实验源代码，
不进入默认生产路径。128 份 checkpoint 载荷、128 份 Visit 文件和守恒 CSV 逐字节一致，
所有 solver 非耗时列一致。详细分 refinement 表、最终动量分区数据、冻结哈希和实验 diff
见 `experiments/2026-09-07-sparse-ibm-candidate.md`。正确性修复保留。

### 明确的剩余边界

非对流自动时间尺度仍不是现有应用能力。需要按具体离散项区分显式稳定性、预测器
可接受性与隐式精度；当前仅记录调用方责任，不以通用扩散 CFL 缩小固定 dt。
标量 remap 系数缓存、reach=1、Picard 预测复用没有混入本轮，缺少有标量大规模成本依据。
64-loop 诊断溢出仍拒绝完整归因；高组成 sweep/retry 需要独立诊断窗口或有界分段设计。
生产内存硬预算、非对流时间尺度和分段观测不能标记为本轮已实现功能。

## 最终冻结与长测恢复

最终生产源提交 `3024bfd9ab4f54abbf858272d9037733ba303f1d`；
与 `e649489` 相比，src/include/runner 无新增差异，性能实验已经撤回。
最后构建的 27 项针对性验收全部通过，167.05 s，见 `final-acceptance.log`。
完整测试输出与先前失败记录均保留，不把最后 27 项说成仓库全部测试。

冻结目录 `method-frozen-module-reviewed-20260907` 位于上述 trial 根目录。
程序 SHA256 为 `45b311d391d8f724965671bf47ba347ab38c76738c1540448dd611306359645e`；
build manifest 的 core/target_source_clean 均为 true。
`FINAL_ACCEPTED.json` 明确本次是零标量 Re3900 恢复门槛，不是七模块全部能力或科学验收。
`FROZEN.sha256` 绑定程序、build manifest、case、spec、启动脚本和该门槛记录。

2026-09-07 11:16:54 +0800 已启动
`hundun-re3900-module-reviewed-20260907.service`，ActiveState=active、SubState=running。
新 run 为 `long-module-reviewed-35000-20260907`。
从已验证的 preopt-pilot generation-5510-72559809913240 精确恢复，剩余 29490 步至 35000；
不从未持久化的 5936 步假装恢复，不使用被撤回的实验程序。
首次核查有 128 个实际 runner 进程，独占锁不能被第二个作业取得，无并行测试/编译。

RUN.meta 核对：restart_method_recovery=0、restart_requires_recovery=0、
statistics_epoch_start_step=1000、sampling_start_step=11001。
继承元数据中的 statistics_reset_reason=method_recovery 指原 epoch 的来源，不是此次重置。
已确认运行接受到 step=5516、time=0.076169960092241279 s：BDF2、attempt=1、retry=0，
continuity=1.4727164461916333e-7、energy=6.2309360502208999e-7、EOS=0、
committed CFL(out)=0.29118850842390626。这里只记录启动核查，不代替长测最终结论。
Visit/checkpoint 周期仍为 500，下一次为绝对步 6000；现有 step5510 的 Visit 已做字节对照。
