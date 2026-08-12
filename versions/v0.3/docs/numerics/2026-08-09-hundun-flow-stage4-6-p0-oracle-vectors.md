# HUNDUN-FLOW Stage 4--6 P0 公开数学向量

日期：2026-08-09
状态：`preflight_candidate`
schema：`hundun.stage4_6_p0.oracle_vectors.v1`

本文只冻结公开方程的低成本输入输出向量。它不包含产品算法、不调用 COAST，也不证明
Stage 4--6 已实现或已通过科学验收。后续产品测试必须从相同公开方程独立实现 RED，不能把
本文的生成器作为运行时 oracle 或替代产品数值路径。

## 1. 向量合同

本文共有 23 个独立 oracle-vector 块。每块都以 `<!-- oracle-vector:ID -->` 开始，并逐行
包含且只包含以下九个字段：`vector_id`、`public_equation_and_citation`、`SI_input`、
`operation_order`、`expected_decimal`、`expected_binary64_hex_when_exact`、
`comparison_class`、`mutation_and_expected_failure`、`capability_not_proved`。块外的
说明不属于向量字段。

比较规则：

- `bitwise`：所有整数、状态和明确列出的 binary64 hex 必须逐位相同；失败报告中的零必须是
  `+0.0`，即 `0x0p+0` 且 `signbit=false`。
- `ulp_bounded(N)`：按本表 `operation_order` 计算，结果与列出的 binary64 值相差不超过
  `N` ULP。只有 `exp`、`sqrt`、`log1p` 等 libm 路径使用该类。
- `relative_absolute(r,a)`：要求 `abs(actual-expected) <= a + r*abs(expected)`。本表中的
  守恒 roundoff 默认 `r=0, a=1e-15`；它只是小型 algebra oracle 的舍入预算，不是产品科学阈值。

公开来源：

- 理想气体、组分焓和反应增量合同：Cantera 3.2
  [C++ 文档](https://cantera.org/3.2/cxx/)与 NASA/TP-2002-211556
  [NTRS 记录](https://ntrs.nasa.gov/citations/20020085330)；
- 对称 `C(dt/2)-T(dt)-C(dt/2)`：Strang
  [10.1137/0705041](https://doi.org/10.1137/0705041)，低马赫反应流时序参见
  Day--Bell [10.1088/1364-7830/4/4/309](https://doi.org/10.1088/1364-7830/4/4/309)；
- ESF/IEM：Valiño [10.1023/A:1009968902446](https://doi.org/10.1023/A:1009968902446)
  和低雷诺一致形式 [10.1007/s10494-015-9687-0](https://doi.org/10.1007/s10494-015-9687-0)；
- Philox：Salmon 等 [10.1145/2063384.2063405](https://doi.org/10.1145/2063384.2063405)及
  [Random123 官方仓库](https://github.com/DEShawResearch/random123)；
- 液滴蒸发、传热和阻力：Abramzon--Sirignano [10.1016/0017-9310(89)90043-4](https://doi.org/10.1016/0017-9310(89)90043-4)、
  Ranz--Marshall 与 Schiller--Naumann 的公开极限式。

## 2. Stage 4：组分、热化学和 source transaction

### 2.1 守恒状态和状态方程

<!-- oracle-vector:S4_SPECIES_SUM_001 -->
vector_id: S4_SPECIES_SUM_001
public_equation_and_citation: 守恒恒等式 sum_i(Y_i)=1；Cantera 3.2 C++ 文档与 NASA/TP-2002-211556
SI_input: Y=[0.1,0.2,0.7]（无量纲）
operation_order: 按索引 0→2 累加
expected_decimal: 1.00000000000000000e+00
expected_binary64_hex_when_exact: 0x1p+0
comparison_class: bitwise
mutation_and_expected_failure: P0 数值 mutation 省略最后一个组分，和变为 0.3，必须失败；字段数/Restart identity 属于延期产品结构 RED。
capability_not_proved: 未证明 simplex projection。

<!-- oracle-vector:S4_ELEMENT_DELTA_001 -->
vector_id: S4_ELEMENT_DELTA_001
public_equation_and_citation: 质量反应增量与元素守恒；manufactured integer molar masses for 2H2+O2->2H2O，来源为 Cantera 3.2 与 NASA/TP-2002-211556
SI_input: 质量反应增量 [-0.004,-0.032,+0.036] kg；元素权重按 H=4/36、O=32/36
operation_order: 先求质量和，再按 4/36、32/36 求元素差
expected_decimal: mass=-6.93889390390722838e-18 kg; H=-8.67361737988403547e-19 kg; O=-6.93889390390722838e-18 kg
expected_binary64_hex_when_exact: not_applicable_roundoff
comparison_class: relative_absolute(0,1e-16)
mutation_and_expected_failure: P0 数值 mutation 把 H 权重 4/36 改为 2/36，得到非零元素差，必须失败。
capability_not_proved: 未证明任意机理元素矩阵或 clipping/projection。

<!-- oracle-vector:S4_IDEAL_GAS_RHO_001 -->
vector_id: S4_IDEAL_GAS_RHO_001
public_equation_and_citation: 理想气体混合物密度 rho=(p0*W)/(R*T)；Cantera 3.2 C++ 文档
SI_input: p0=101325 Pa; W=0.02897 kg/mol; T=1200 K; R=8.31446261815324 J/(mol K)
operation_order: 严格计算 (p0*W)/(R*T)
expected_decimal: 2.94204747479317574e-01 kg/m3
expected_binary64_hex_when_exact: 0x1.2d44026301a56p-2
comparison_class: bitwise
mutation_and_expected_failure: P0 数值 mutation 把 W 放入分母，结果必须不同。
capability_not_proved: 未证明机械压力隔离、单位 schema 或 Cantera EOS 反演。

<!-- oracle-vector:S4_HTC_001 -->
vector_id: S4_HTC_001
public_equation_and_citation: 总热化学焓 h=sum_i(Y_i*(h_f,i+cp_i*(T-Tref)))；Cantera 3.2 与 NASA/TP-2002-211556
SI_input: Y=[0.25,0.75]; h_f=[-1e6,2e5] J/kg; cp=[1000,1200] J/(kg K); Tref=300 K; T=800 K
operation_order: 逐组分先算 h_f+cp*(T-Tref)，再乘 Y，按索引累加
expected_decimal: 4.75000000000000000e+05 J/kg
expected_binary64_hex_when_exact: 0x1.cfdep+18
comparison_class: bitwise
mutation_and_expected_failure: sensible-only mutation 得 5.75000000000000000e+05 J/kg，必须失败。
capability_not_proved: 未证明温变 cp 或 NASA 多项式。

### 2.2 积分源、失败发布和对称耦合

<!-- oracle-vector:S4_INTEGRATED_RATE_001 -->
vector_id: S4_INTEGRATED_RATE_001
public_equation_and_citation: 线性速率积分 delta=int r dt；Cantera source transaction 公开 algebra
SI_input: r(t)=r0+slope*t; r0=0.01 kg/(m3 s); slope=0.04 kg/(m3 s2); dt=0.5 s；A/B 取负/正
operation_order: delta=r0*dt+0.5*slope*dt*dt；先得到 delta，再写 A=-delta、B=+delta，并检查 A+B canonical zero
expected_decimal: A=-1.00000000000000002e-02 kg/m3; B=+1.00000000000000002e-02 kg/m3; A+B=+0.0
expected_binary64_hex_when_exact: A=-0x1.47ae147ae147bp-7; B=+0x1.47ae147ae147bp-7; sum=0x0p+0
comparison_class: bitwise
mutation_and_expected_failure: endpoint-rate mutation 得 1.49999999999999994e-02，或破坏 A+B canonical zero，必须失败。
capability_not_proved: 未证明 stiff chemistry 精度。

<!-- oracle-vector:S4_FAILED_ZERO_001 -->
vector_id: S4_FAILED_ZERO_001
public_equation_and_citation: 失败 chemistry interval 的 publish-last 原子性与 canonical zero 合同
SI_input: 任意失败 chemistry interval；accepted_time_before=1.25 s；rejected_dt=0.1 s
operation_order: standalone contract 先判定 rejected，再令全部 delta 为 canonical +0.0，且 accepted time 保持不变；真实 publish-last 由产品 RED 验证
expected_decimal: 每个数值成员 +0.0；accepted_time_after=1.25 s
expected_binary64_hex_when_exact: delta=0x0p+0, signbit=false; accepted_time_after=0x1.4p+0
comparison_class: bitwise
mutation_and_expected_failure: P0 数值 mutation 使用 1e-12、-0.0 或把 rejected interval 的时间推进到 1.35 s，值/状态必须失败。
capability_not_proved: 真实 partial publication、MPI collective rollback 与 scheduler 事务仍延期到产品结构 RED。

<!-- oracle-vector:S4_COMMUTING_CTC_001 -->
vector_id: S4_COMMUTING_CTC_001
public_equation_and_citation: 对称 Strang C(dt/2)-T(dt)-C(dt/2)，来源 Strang DOI 10.1137/0705041 与 Day--Bell DOI 10.1088/1364-7830/4/4/309
SI_input: q0=0.75; C:q'=-0.4 s^-1 q; T:q'=0.1 s^-1 q; dt=0.2 s
operation_order: 依次乘 exp(-0.4*dt/2)、exp(0.1*dt)、exp(-0.4*dt/2)；reverse 使用相反两个 half chemistry 与一个 full transport
expected_decimal: forward=7.06323400188186401e-01; analytic=7.06323400188186512e-01; reverse=7.49999999999999778e-01
expected_binary64_hex_when_exact: forward=0x1.69a338806a624p-1; analytic=not_applicable; reverse=not_applicable
comparison_class: ulp_bounded(2)
mutation_and_expected_failure: P0 数值 mutation 把两个 chemistry half 都推进 full dt，必须失败。
capability_not_proved: 只验证 half-duration/count 和可逆性；可交换标量不能识别 C/T 排序或一般 stage-time 错误，也不证明 reacting-flow 整体二阶。

`S4_COMMUTING_CTC_001` 的 reverse 依次使用相反的两个 half chemistry 和一个 full transport。它是
时序/步长 RED，不允许据此把 HUNDUN 的 BDF2 momentum map 称为经典 Strang 子流。

## 3. Stage 5：Philox、平衡 Wiener、IEM 和一致性

### 3.1 Philox 与平衡增量

Random123 固定为官方 tag `v1.14.0`，commit
`726a093cd9a73f3ec3c8d7a70ff10ed8efec8d13`。官方 `tests/kat_vectors` 的 SHA-256 为
`aab5ebabf40003f63d6d87b24cbd2c8a02652e00cf8bad64226fd50586929183`。

<!-- oracle-vector:S5_PHILOX4X32_10_ZERO_001 -->
vector_id: S5_PHILOX4X32_10_ZERO_001
public_equation_and_citation: Philox4x32-10 published KAT；Salmon 等 DOI 10.1145/2063384.2063405 与 Random123 官方仓库
SI_input: counter=[0,0,0,0]; key=[0,0]; rounds=10
operation_order: 每轮按论文规定的两个 32x32→64 乘法、xor 和 Weyl key update
expected_decimal: not_applicable_integer_word_contract; output words=6627e8d5 e169c58d bc57ac4c 9b00dbd8
expected_binary64_hex_when_exact: not_applicable_integer_word_contract
comparison_class: bitwise
mutation_and_expected_failure: 与下一非零 KAT 共同杀 9 rounds、漏 Weyl key bump 和反转 counter word order；三种 mutation 均必须不同。
capability_not_proved: 不证明序列化 byte endianness、统计质量或 Gaussian transform。

<!-- oracle-vector:S5_PHILOX4X32_10_NONZERO_001 -->
vector_id: S5_PHILOX4X32_10_NONZERO_001
public_equation_and_citation: Philox4x32-10 official KAT line 28；Salmon 等 DOI 10.1145/2063384.2063405 与 Random123 官方仓库
SI_input: counter=[243f6a88,85a308d3,13198a2e,03707344]; key=[a4093822,299f31d0]; rounds=10
operation_order: 10 rounds；保持 counter word order、乘法高低字、xor 与 Weyl key bump 的论文顺序
expected_decimal: not_applicable_integer_word_contract; output words=d16cfe09 94fdcceb 5001e420 24126ea1
expected_binary64_hex_when_exact: not_applicable_integer_word_contract
comparison_class: bitwise
mutation_and_expected_failure: P0 生成器实际执行 9-round、no-key-bump 和 reversed-counter-word mutations，三者均必须不同。
capability_not_proved: serialized byte encoding 延期；不证明并行地址策略。

<!-- oracle-vector:S5_WIENER_N2_001 -->
vector_id: S5_WIENER_N2_001
public_equation_and_citation: 平衡 antithetic Wiener increment；ESF/IEM 低雷诺一致形式 DOI 10.1007/s10494-015-9687-0
SI_input: dt=0.25 s；dichotomic signs=[+1,-1]；按 field 顺序
operation_order: 先算 sqrt(dt)=0.5 s^0.5，再逐 field 乘 signs
expected_decimal: [+0.5,-0.5] s^0.5; sum canonical +0.0
expected_binary64_hex_when_exact: [+0x1p-1,-0x1p-1]; sum=0x0p+0
comparison_class: bitwise
mutation_and_expected_failure: P0 数值 mutation 用 dt 代替 sqrt(dt)，幅值变 0.25，必须失败。
capability_not_proved: key/address/retry 是延期产品结构 RED。

<!-- oracle-vector:S5_WIENER_N4_001 -->
vector_id: S5_WIENER_N4_001
public_equation_and_citation: N=4 antithetic Wiener increments；Valiño DOI 10.1023/A:1009968902446 与低雷诺一致形式 DOI 10.1007/s10494-015-9687-0
SI_input: dt=0.25 s；两个 antithetic pairs；冻结 signs=[+1,-1,-1,+1]
operation_order: 先算共享 sqrt(dt)=0.5，再按 field/pair 顺序生成四个符号增量并求 pair mean
expected_decimal: [+0.5,-0.5,-0.5,+0.5] s^0.5; sum canonical +0.0
expected_binary64_hex_when_exact: [0x1p-1,-0x1p-1,-0x1p-1,0x1p-1]; sum=0x0p+0
comparison_class: bitwise
mutation_and_expected_failure: 同一 dt mutation 与 pair-sum 检查必须失败；same-sign pair 也必须失败。
capability_not_proved: mutable cursor、retry redraw、rank/cell/species key 污染和样本数充分性延期到产品结构 RED。

<!-- oracle-vector:S5_LAMINAR_NOISE_001 -->
vector_id: S5_LAMINAR_NOISE_001
public_equation_and_citation: stochastic coefficient sqrt(2*mu_t/(rho*Sc_t))；ESF/IEM DOI 10.1023/A:1009968902446
SI_input: mu=1.8e-5 Pa s; mu_t=+0.0 Pa s; rho=1.2 kg/m3; Sc_t=0.7
operation_order: 只以 sqrt(2*mu_t/(rho*Sc_t)) 形成随机系数，不把 molecular mu 加入根号
expected_decimal: +0.0 m/s^0.5
expected_binary64_hex_when_exact: 0x0p+0; signbit=false
comparison_class: bitwise
mutation_and_expected_failure: P0 数值 mutation 把 molecular mu 加入根号，得到 6.54653670707977184e-03 m/s^0.5，必须失败。
capability_not_proved: 未证明有限 mu_t 系数或完整 SPDE 离散。

### 3.2 exact IEM 与守恒扰动

对本 P0 冻结的 IEM 合同，低雷诺一致形式中的连续混合方程为

\[
\frac{d\phi_i}{dt}=-\frac{\phi_i-\bar\phi}{2\tau}.
\]

在单次 `dt` mixing interval 上精确积分，affine factor 是

\[
f=\exp\!\left(-\frac{\Delta t}{2\tau}\right),\qquad
\phi_i^{n+1}=\bar\phi+f(\phi_i^n-\bar\phi).
\]

<!-- oracle-vector:S5_IEM_N4_001 -->
vector_id: S5_IEM_N4_001
public_equation_and_citation: exact IEM affine relaxation dphi_i/dt=-(phi_i-mean)/(2*tau)；Valiño DOI 10.1023/A:1009968902446
SI_input: phi=[0.1,0.3,0.7,0.9]; mean=0.5; dt=0.4 s; tau=0.2 s
operation_order: 在完整 dt interval 上先算一个共享 factor=exp(-dt/(2*tau))，再按 field 顺序 affine update
expected_decimal: factor=3.67879441171442334e-01; final=[0.352848223531423089,0.426424111765711544,0.573575888234288511,0.647151776468576911]; mean=0.5
expected_binary64_hex_when_exact: not_applicable_transcendental
comparison_class: ulp_bounded(2)
mutation_and_expected_failure: P0 数值 mutation 删除模型内生的 1/2，factor 变 1.35335283236612702e-01，必须失败；共享 factor 被 field-specific 替换也必须失败。
capability_not_proved: field-specific mean 与 TCR feedback 延期到产品结构 RED。

<!-- oracle-vector:S5_ELEMENT_NULL_PAIR_001 -->
vector_id: S5_ELEMENT_NULL_PAIR_001
public_equation_and_citation: simplex/element/mean-preserving null pair；Cantera composition identity 与公开 element matrix algebra
SI_input: mean Y=[0.2,0.3,0.1,0.4]; delta=[0,0.02,0.02,-0.04]; fields=mean +/- delta; element row E0=[1,0,0.5,0.25]
operation_order: 分别形成 mean+delta 与 mean-delta，再逐组分求 sum、pair mean 和 E0*Y
expected_decimal: 两 field 的 sum(Y)=1；E0*Y=0.35（roundoff <=1e-15）；pair mean 逐组分恢复原 mean
expected_binary64_hex_when_exact: weight_sum=0x1p+0; deposited=0x1.0624dd2f1a9fcp-7
comparison_class: relative_absolute(0,1e-15)
mutation_and_expected_failure: P0 数值 mutation 将第二场也取 mean+delta，pair mean 不再等于 mean，必须失败。
capability_not_proved: clipping、逐场归一化、coupled projection 和 only-field-0 修复延期到产品结构 RED。

### 3.3 TCR 候选边界

TCR 只引用用户已确认的两篇论文：[10.1016/j.proci.2026.106128](https://doi.org/10.1016/j.proci.2026.106128)
和 [10.1016/j.cja.2026.104123](https://doi.org/10.1016/j.cja.2026.104123)。P0 不从摘要补写生产公式，也不访问
COAST。以下块是未执行的状态合同，不是数值向量。

<!-- oracle-vector:S5_TCR_KAPPA1_CANDIDATE_001 -->
vector_id: S5_TCR_KAPPA1_CANDIDATE_001
public_equation_and_citation: TCR kappa=1 退化极限候选，只引用 DOI 10.1016/j.proci.2026.106128 与 DOI 10.1016/j.cja.2026.104123
SI_input: mode=shadow; kappa=1；其余输入与 S5_IEM_N4_001 相同；只允许走共享 mixer
operation_order: 仅冻结 candidate_before_COAST_oracle 状态，不执行 TCR 数值更新或 COAST 对照
expected_decimal: not_applicable_contract
expected_binary64_hex_when_exact: not_applicable_contract
comparison_class: bitwise
mutation_and_expected_failure: 未执行 TCR mutation；任何把候选状态当作已复算输出等价的检查必须失败。
capability_not_proved: 明确不声称 kappa=1 数值等价已复算；独立 TCR/IEM 更新器、shadow feedback、algebra、root/history、feedback timing 和 COAST-equivalence 全部延期。

完整 algebra、九值 history、root/status 和 feedback timing 必须在 Stage 5 的公开全文推导与用户确认的进程外 COAST oracle gate 后冻结。

## 4. Stage 6：parcel/gas 交换与低成本极限

<!-- oracle-vector:S6_PAIRWISE_EXCHANGE_001 -->
vector_id: S6_PAIRWISE_EXCHANGE_001
public_equation_and_citation: parcel/gas pairwise conservation of mass, momentum and total thermochemical energy；公开液滴交换极限式
SI_input: parcel mass=-0.00125 kg; momentum=[-0.02,+0.04,-0.01] kg m/s; energy=-125 J；gas 为严格相反值
operation_order: 逐分量相加 parcel 与 gas；同时检查 energy direction 为 parcel loss/gas gain
expected_decimal: mass=+0.0; momentum=[+0.0,+0.0,+0.0] kg m/s; total thermochemical energy=+0.0 J
expected_binary64_hex_when_exact: 每个成员=0x0p+0; signbit=false
comparison_class: bitwise
mutation_and_expected_failure: P0 数值 mutation 把 gas mass 改成同号，或反转 energy direction，闭合差/方向必须失败。
capability_not_proved: 未证明实际相关式；遗漏 vapor momentum/enthalpy、double deposition 与 MPI rollback 延期到产品结构 RED。

<!-- oracle-vector:S6_D2_001 -->
vector_id: S6_D2_001
public_equation_and_citation: Abramzon--Sirignano D-squared evaporation law，DOI 10.1016/0017-9310(89)90043-4
SI_input: d0=1e-4 m; K=1e-8 m2/s; dt=0.25 s
operation_order: 严格先算 d2=d0*d0-K*dt，验证正值后 sqrt(d2)
expected_decimal: d2=7.50000000000000098e-09 m2; d=8.66025403784438641e-05 m
expected_binary64_hex_when_exact: d2=0x1.01b2b29a4692cp-27
comparison_class: ulp_bounded(2)
mutation_and_expected_failure: P0 数值 mutation 把蒸发减号改为加号，d2=1.24999999999999994e-08 m2，必须失败。
capability_not_proved: 负 d2、消失事件和跨事件步长延期到 Stage 6 产品 RED。

<!-- oracle-vector:S6_DRAG_RELAX_001 -->
vector_id: S6_DRAG_RELAX_001
public_equation_and_citation: constant-gas Stokes drag relaxation，Schiller--Naumann 公开低 Re 极限
SI_input: gas u=10 m/s; parcel v0=2 m/s; dt=0.1 s; tau=0.2 s
operation_order: v=u+(v0-u)*exp(-dt/tau)
expected_decimal: 5.14775472229893261 m/s
expected_binary64_hex_when_exact: 0x1.4974d039069f6p+2
comparison_class: ulp_bounded(2)
mutation_and_expected_failure: P0 数值 mutation 用 explicit Euler，得 6 m/s，必须失败。
capability_not_proved: 未证明 Schiller--Naumann finite-Re integration 或产品 drag sign wiring。

<!-- oracle-vector:S6_TRANSFER_LIMIT_001 -->
vector_id: S6_TRANSFER_LIMIT_001
public_equation_and_citation: dimensionless transfer limits；Abramzon--Sirignano DOI 10.1016/0017-9310(89)90043-4 与 Ranz--Marshall 静止球极限
SI_input: B_M=0；Re=0；球形 parcel
operation_order: 先算 log1p(B_M)，再取静止球 Nu=2、Sh=2
expected_decimal: log1p(B_M)=+0.0; Nu=2; Sh=2
expected_binary64_hex_when_exact: log1p(B_M)=0x0p+0; Nu=2; Sh=2
comparison_class: bitwise
mutation_and_expected_failure: P0 数值 mutation 使用 log(B_M) 产生非有限值，或把静止球极限常数设为 0，必须失败。
capability_not_proved: 不把 Nu=Sh=2 误称为零传热/传质；未证明 film-property interpolation。

<!-- oracle-vector:S6_DEPOSITION_001 -->
vector_id: S6_DEPOSITION_001
public_equation_and_citation: normalized source deposition weights；公开 parcel sampling/deposition conservation identity
SI_input: dimensionless weights=[0.125,0.375,0.25,0.25]; source=0.008 kg/s
operation_order: 按索引形成 sum(w_i*source)，先求 weight sum 再发布 owner contribution
expected_decimal: weight sum=1.0; deposited=8.00000000000000017e-03 kg/s
expected_binary64_hex_when_exact: not_applicable_roundoff
comparison_class: bitwise
mutation_and_expected_failure: P0 数值 mutation 重复登记 owner 0 contribution，结果必须不再等于 source。
capability_not_proved: sampling/deposition stencil identity、漏 cell 和 MPI halo 延期到产品结构 RED。

<!-- oracle-vector:S6_REJECT_ZERO_001 -->
vector_id: S6_REJECT_ZERO_001
public_equation_and_citation: rejected-trial publish contract：失败候选的 gas、parcel 与 injector residual 必须全为 canonical zero
SI_input: gas mass/energy/3-momentum、parcel-state delta 和 injector residual 共七个发布成员
operation_order: standalone contract array 逐成员初始化为 +0.0，并令 rejected_trial_is_published=false；真实 publish-last 顺序由产品 RED 验证
expected_decimal: 七个发布成员均为 canonical +0.0; rejected_trial_is_published=false
expected_binary64_hex_when_exact: 七个成员均为 0x0p+0, signbit=false; publish-state=false
comparison_class: bitwise
mutation_and_expected_failure: P0 数值 mutation 把 gas energy 改为 1e-12，完整数组比较必须失败；rejected publish-state token mutation 也必须失败。
capability_not_proved: 真实 publish-last、parcel deletion ordering、scheduler 状态和 collective rollback 延期到产品结构 RED。

### 4.1 双液体身份行

P0 要求至少两个性质明显不同的纯液体身份，但不允许把未获许可的生产系数写入 Git。本节的
两个 identity 合同只冻结来源、内容和消费状态。消费规则是：

1. n-dodecane/kerosene-like 与 iso-octane/gasoline-like 必须拥有不同的
   `source_revision + content_sha256 + phase + unit_schema` fingerprint；
2. dispatch 只使用 fingerprint，名称只是 label；
3. 缺少资产级再分发许可时状态只能是 `user_supplied_only`，数值系数不得进入本文；
4. 两个候选的 exact identity、许可状态和解析条件以
   `docs/references/2026-08-09-hundun-flow-stage4-p0-fuel-data-candidates.md` 为权威；
5. 在 P0-3 C++ parse smoke 前，不得把任何候选写成 `validated`。

<!-- oracle-vector:S6_N_DODECANE_PROPERTY_IDENTITY_001 -->
vector_id: S6_N_DODECANE_PROPERTY_IDENTITY_001
public_equation_and_citation: n-Dodecane pure-liquid property identity contract；P0-4 source candidate and public provenance record
SI_input: P0-4 commit=3ed7f79ee4348d9402d0b475ea455ff1724777b1; Git blob=a63d749da585d0a187ad2b89108e33bc874fc524; document_sha256=305dad119bbb3e0a75cb73071e7516f95c319102bf9fcb1da92706b14ab61d77
operation_order: canonical identity string=source_revision=CoolProp:v8.0.0:ae81610e7d23efc57f9d051c8e70a4d66e87537f;content_sha256=abf6ce7d6ef8a2a492e1ee3f4434d4bf7ebe2087e205602b5a351e2786b4ed1c;phase=pure_liquid_property_source:n-Dodecane;unit_schema=source_native_to_SI_pending_p0_4_step4
expected_decimal: not_applicable_contract; status=bundled_candidate; output_pack=unresolved
expected_binary64_hex_when_exact: not_applicable_contract; fingerprint_sha256=e20546afaa8f0366c4d4cc8ca562b393f0829229f29ea6ccbb0f26c21859047d
comparison_class: bitwise
mutation_and_expected_failure: Changing P0-4 commit/blob/document SHA, any canonical identity component, or dispatching by fuel name instead of fingerprint must fail.
capability_not_proved: 不证明 property coefficient 数值正确性、parse smoke、许可升级或最终 output pack；状态仅为 bundled_candidate。

<!-- oracle-vector:S6_ISO_OCTANE_PROPERTY_IDENTITY_001 -->
vector_id: S6_ISO_OCTANE_PROPERTY_IDENTITY_001
public_equation_and_citation: iso-octane (2,2,4-trimethylpentane) pure-liquid property identity contract；P0-4 generator/table candidate and public provenance record
SI_input: P0-4 commit=3ed7f79ee4348d9402d0b475ea455ff1724777b1; Git blob=a63d749da585d0a187ad2b89108e33bc874fc524; document_sha256=305dad119bbb3e0a75cb73071e7516f95c319102bf9fcb1da92706b14ab61d77
operation_order: canonical identity string=source_revision=FuelLib:v3.0.0:4e26169ccbc456f5f1e6b4f97d6832fa94192f91;content_sha256=sdist:92633a64168878e2b3694c67fc6b66a670cebfd8e954be9e1de102c3cbe3ee02+gcmTable:a5fa6fb586f7cbcf4282d826dc1fd09f367919437d002f9c2b0d6eeeb7d43e8a+fuel.py:cbeb754b4a79f774bf782dd497b7a30e4212f19ac09c27bb648cfe9a581aee4a;phase=pure_liquid_property_generator_candidate:2,2,4-trimethylpentane;unit_schema=generated_output_SI_pending_p0_4_step4
expected_decimal: not_applicable_contract; status=generator/table bundled_candidate; final_output_pack=unresolved
expected_binary64_hex_when_exact: not_applicable_contract; fingerprint_sha256=b36b8dfca48a245fa7c9034db2cefd31dce2750febeaa6f0f517da6579ffa646
comparison_class: bitwise
mutation_and_expected_failure: Changing P0-4 commit/blob/document SHA, any canonical identity component, or dispatching by fuel name instead of fingerprint must fail.
capability_not_proved: 不证明 generator/table 系数数值正确性、parse/property smoke、许可升级或最终 output pack；属性求值和 final output pack 均未完成。

## 5. 独立复算与 mutation 证据

### 5.1 P0 已执行数值 mutation 与延期结构 RED

v3 生成器对下列 mutation 实际计算了错误候选，并验证其结果或状态与 expected 不同；新增的
generator checks 也计入 43 个具名检查：

- Stage 4 composition/thermo：omit dependent species、wrong element weight、reciprocal molecular-weight placement、sensible-only enthalpy；
- Stage 4 source/time：endpoint rate、nonzero failed delta、两个 full chemistry intervals 代替两个 half intervals、integrated A+B canonical zero、rejected interval 错误推进 accepted time；
- Stage 5 RNG：Philox 9 rounds、missing Weyl key bump、reversed counter-word order、Wiener `dt` 代替 `sqrt(dt)`；
- Stage 5 closure：molecular viscosity enters noise、IEM missing intrinsic `1/2`、same-sign field pair、pair mean preservation；
- Stage 6 exchange/models：same-sign gas/parcel mass、energy direction inversion、D-squared wrong sign、explicit-Euler drag、`log(0)`、wrong static `Nu`；
- Stage 6 publication：duplicated owner deposition、one nonzero member in rejected publication array、rejected publish-state token。

真实 scheduler、publish-last 与 MPI collective rollback 仍延期到相应产品 mutation-sensitive RED。以下项目需要
registry、transaction、scheduler、stencil、MPI 或已确认 COAST oracle 上下文，因此明确延期，本文件不把它们记为已覆盖：

- schema/Restart field identity、真实 publish-last、parcel deletion ordering、collective rollback；
- RNG address 中的 rank/cell/species 污染、mutable cursor、retry redraw 和 serialized byte encoding；retry clock 仅由本 v3 生成器 token check 覆盖；
- clipping、逐场归一化、coupled element projection、TCR root/history/feedback；
- sampling/deposition stencil identity、MPI halo owner、负 d-squared 消失事件；
- C/T 非交换排序、一般 stage-time 错误和 reacting-flow 整体二阶。

### 5.2 手算/闭式独立复算

第二种方法不复用 C++ 生成器控制流：20 个数值向量逐项以手算代数、守恒恒等式或已发表
known-answer test 复核；TCR 与两个燃料 identity 明确排除在数值复算之外，只做状态/hash 合同。
Git 外记录 `closed-form-review-v3.md` 的 SHA-256 为
`51a055b9aaac16566519ebb92905c766ec56c0841ea87cf97b0191390e26ffa4`，覆盖清单为 20/20。

### 5.3 三工具链独立复算（v3 evidence）

Git 外独立 C++17 计算器：

```text
source:
/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/manifests/oracle-vectors-v1/oracle_vectors.cpp
source_sha256:
4396080e4c81b1f55bdc6954af602ab71eb5045dfeb9354ba0d3f9bc9786ad86

host_compiler: /usr/bin/g++ (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0
compiler_binary_sha256:
4f4a7fae7ebe488fd2ae29a2dc5b870eafc278cc1ccd2ca9e8ebea251875c56d
binary_sha256:
6c4bc8536d07b90b3c94bf564aaa307899a6c6c6a39844d009db89b517709e84
compile_log_sha256:
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
output_sha256:
ae39da1fc0cef680b8426784159c525112d8552d394a21ce04dc64d1ac5003c0
```

编译命令固定为：

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
  /home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/manifests/oracle-vectors-v1/oracle_vectors.cpp \
  -o /home/wyf/code_dev/.hundun-flow-preflight/stage4-p0/manifests/oracle-vectors-v1/oracle_vectors_v3_gcc7
```

`/home/wyf/.local/bin/clang++`（Clang 15.0.6）对同一 source SHA 独立编译：

```text
compiler_binary_sha256:
388be41dc565a891ced9e78da2e89a249ca9b9a26f71a3c912e8ba89585be89c
binary_sha256:
8ef15146a8b9af02bff05d3a9969266c8f079774b5e97108f25be69f8e79a555
compile_log_sha256:
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
output_sha256:
ae39da1fc0cef680b8426784159c525112d8552d394a21ce04dc64d1ac5003c0
```

目标交叉检查在 P0-2 的只读 Jammy rootfs 中完成，使用 rootfs 内的
`/usr/bin/g++-11 (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0`：

```text
compiler_binary_sha256:
2360901d864cf10bfd6296e261cb2c14053552a80377761ab07146ec9ec9a2c0
binary_sha256:
54ebb623c0a8790469254e6dd532ac4120671f3624bb159089e8c688b39692b3
compile_log_sha256:
e2018c8113b56c2954ea21138c9832061ea14c1c5b6b03831747145e03624581
output_sha256:
ae39da1fc0cef680b8426784159c525112d8552d394a21ce04dc64d1ac5003c0
```

GCC 7、Clang 15、Jammy GCC 11 的 43 项具名检查全部 `PASS`，`check.failure_count=0`，且
三工具链 output SHA 均为 `ae39da1fc0cef680b8426784159c525112d8552d394a21ce04dc64d1ac5003c0`。
新增 `integrated A+B canonical zero`、`pair mean`、`energy direction`、`rejected publish-state token`
和 `retry clock` checks 均包含在 43 项中；没有放宽任何 ULP 或 absolute budget。Philox、整数、
canonical zero、exact binary fractions 和本表列出的 libm 结果均保持不变。该交叉检查只接受本文件的
公开数学向量，不接受 P0-2 Cantera artifact 或任何 Stage 4--6 产品能力。

旧 source/binary/output 证据仍按原文件名保留在 Git 外，但不得混用于当前 v3 文档；一次使用不存在命令名
`clang++-15` 的失败调用也保留在 Git 外，不计入 Clang GREEN。有效 Clang 证据绑定其绝对 compiler 路径和 SHA。

### 5.4 v2 validator/hash 证据的历史边界

旧 v2 validator、hash 和 mutation 日志证据保留为历史记录，不删除其存在说明，也不把它们混入当前
schema。当前 v3 的 document、validator、runner、三工具链 output 与 mutation summary SHA 统一记录在
Git 外 `oracle-vectors-v1/document-freeze-v3.receipt`。该 receipt 在正文冻结后生成，正文不反向嵌入
receipt SHA，从而避免 document/validator/receipt 循环自引用。

## 6. 明确未证明的能力

本文不证明：Cantera 状态反演或 stiff integration、真实机理版权/正确性、PISO/IBM/WALE
耦合、reacting-flow 二阶、TCR feedback、COAST 相似性、Abramzon--Sirignano 产品实现、
parcel migration、Checkpoint/Restart、MPI decomposition、性能或任何 Stage 4--6 最终能力。
它只把公开数学的低成本 RED 起点前移到 P0；真实 scheduler、publish-last 和 MPI rollback 仍须由产品结构 RED 验证。
