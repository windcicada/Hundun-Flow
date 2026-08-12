# HUNDUN-FLOW Stage 4 P0 双代理燃料公开数据候选

日期：2026-08-09

状态：`P0-4_COMPLETE / PREFLIGHT_PARTIAL`

本文件先在 commit `3ed7f79` 冻结 Task P0-4 的 public-provenance/legal subcluster，随后在
已接受的 P0-3 边界上完成 Step 4 的外部机制解析和双液体属性接口验证。`P0-4_COMPLETE`
表示本任务要求的低成本接口检查已经闭合；`PREFLIGHT_PARTIAL` 表示两条公开候选链仍没有
同时满足再分发许可、气相输运完整性和产品接受门，不能解释为 bundle-ready 或科学接受。

本文只登记公开、官方或原始论文来源的候选身份、许可和低成本外部检查。解析结果和候选
property pack 全部位于 P0 外部生成树，不是 HUNDUN 产品文件，也不接受 Stage 4/6 产品能力，
不证明真实航空煤油、汽油、着火、火焰、喷雾燃烧或 COAST 相似性。公开检索没有访问 COAST、COAST-2、
BOFFIN、研究源码、算例、数据或可执行文件。

## 1. 状态和版权层次

本文使用以下保守状态：

| 状态 | 含义 |
|---|---|
| `bundled_candidate` | 公开资产有固定身份和明确再分发许可；仍须通过 consumed-file、数值、格式和 notice gate，不能解释为已经打包或接受 |
| `user_supplied_only` | HUNDUN 不再分发；用户只能在其合法获得并按精确 SHA-256 提供后使用 |
| `rejected` | 因许可、身份、完整性、成本或当前接口范围不合格，不进入候选链 |
| `unresolved` | 公开证据不足以决定或缺少完整资产；不得消费 |

必须分别判断三类对象：

1. Cantera/FuelLib/FDS 等软件及其源码许可证；
2. kinetics、thermo 和 gas-transport 机理数据自身的版权与许可；
3. 纯液体实验数据、拟合系数、生成表和最终 property pack 的版权与许可。

Cantera 的 BSD-3-Clause 不自动覆盖机理文件。CoolProp、FuelLib 或 FDS 的许可也不能自动
覆盖原始论文正文、图表、第三方数据库或另一个项目中的数值资产。

## 2. 冻结候选登记表

`sha256` 是从表内官方 HTTPS 资产读取后计算的完整 64 位摘要，除另有说明外不是上游签名
checksum。下表保留 commit `3ed7f79` 的初始检索状态，便于审计当时的选择依据；其中
`pending_later_external_*` 和 property `unresolved` 字样由第 5.1--5.4 节的后续 hash-bound
验证结果取代，但资产的再分发状态没有因此放宽。

| surrogate_family | asset_role | public_name | official_url | release_or_revision | sha256 | citation_or_doi | copyright_holder | license_or_terms_url | redistribution_status | Cantera_parse_status | intended_low_cost_check | excluded_claims |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `n_dodecane_kerosene_surrogate` | `gas_mechanism_kinetics_thermo` | Cantera `nDodecane_Reitz.yaml` | [raw file](https://raw.githubusercontent.com/Cantera/cantera/v3.2.0/data/nDodecane_Reitz.yaml); [release](https://github.com/Cantera/cantera/releases/tag/v3.2.0) | Cantera `v3.2.0`; peeled commit `4a8358eb80cfeb50474386b5f9ec0b3a83519889`; file conversion date `2019-12-11` | `3d3b59ed91dec0d0bcbac2fa2ef2cba13fbd565bfff8f266ca847ed6aa92f7f1` | Wang, Ra, Jia, Reitz, DOI [10.1016/j.fuel.2014.07.028](https://doi.org/10.1016/j.fuel.2014.07.028) | mechanism authors are named; asset-specific copyright holder is not established by the file | [Cantera code license](https://raw.githubusercontent.com/Cantera/cantera/v3.2.0/License.txt); [data disclaimer](https://raw.githubusercontent.com/Cantera/cantera/v3.2.0/data/README.md) | `user_supplied_only` | `pending_later_external_direct_yaml`; phases advertised as `nDodecane_IG` and `nDodecane_RK`; no `transport:` block found in the frozen file | future external hash-locked load; verify phase, exact parsed counts, units, thermo ranges and transport absence | no Cantera-BSD inheritance; no mixture-averaged transport claim; no ignition/soot/kerosene validation |
| `n_dodecane_kerosene_surrogate` | `alternate_gas_mechanism_kinetics` | LLNL `NC12H26_Hybrid_2019-10-17_mech.txt` | [mechanism](https://combustion.llnl.gov/sites/combustion/files/NC12H26_Hybrid_2019-10-17_mech.txt); [official page](https://combustion.llnl.gov/mechanisms/alkanes/n-dodecane) | file revision `2019-10-17`; `LLNL-MI-794264`; contact Simon Lapointe | `62b74c20946cec713e02e347ab03603e70fb20af18d8cd525c0ef7f1dba72d89` | Lapointe et al., DOI [10.1016/j.proci.2018.06.139](https://doi.org/10.1016/j.proci.2018.06.139) | LLNL/LLNS and named mechanism authors; asset-specific grant is absent | [LLNL copyright status](https://www.llnl.gov/disclaimer) | `user_supplied_only` | `pending_later_external_ck2yaml_kinetics`; requires the separately registered thermo asset and has no frozen transport asset | future external conversion; record exact parsed counts, command/output SHA and missing transport status | no redistribution; no count inferred from lexical scans; no validated Cantera transport or kerosene claim |
| `n_dodecane_kerosene_surrogate` | `alternate_gas_mechanism_thermo` | LLNL `NC12H26_Hybrid_2019-10-17_therm.txt` | [thermo](https://combustion.llnl.gov/sites/combustion/files/2020-10/NC12H26_Hybrid_2019-10-17_therm.txt); [official page](https://combustion.llnl.gov/mechanisms/alkanes/n-dodecane) | file revision `2019-10-17`; paired publication identity `LLNL-MI-794264` | `8a9a4324b87ee023d3971d5eb7f91daf49daa1b650015897503fb45fc38765c8` | Lapointe et al., DOI [10.1016/j.proci.2018.06.139](https://doi.org/10.1016/j.proci.2018.06.139) | LLNL/LLNS and named mechanism authors; asset-specific grant is absent | [LLNL copyright status](https://www.llnl.gov/disclaimer) | `user_supplied_only` | `pending_later_external_ck2yaml_thermo`; usable only with separately hash-locked kinetics | future external conversion; verify species coverage, thermo ranges, command/input ordering and output SHA | no redistribution; no standalone kinetics, transport or kerosene validation claim |
| `n_dodecane_kerosene_surrogate` | `pure_liquid_property_source` | CoolProp `n-Dodecane.json` | [raw JSON](https://raw.githubusercontent.com/CoolProp/CoolProp/v8.0.0/dev/fluids/n-Dodecane.json); [official fluid page](https://coolprop.org/fluid_properties/fluids/n-Dodecane.html) | CoolProp `v8.0.0`; peeled commit `ae81610e7d23efc57f9d051c8e70a4d66e87537f`; release `2026-06-28` | `abf6ce7d6ef8a2a492e1ee3f4434d4bf7ebe2087e205602b5a351e2786b4ed1c` | Lemmon--Huber DOI [10.1021/ef0341062](https://doi.org/10.1021/ef0341062); Huber--Laesecke--Perkins DOI [10.1021/ef034109e](https://doi.org/10.1021/ef034109e); Mulero et al. DOI [10.1063/1.4768782](https://doi.org/10.1063/1.4768782) | Ian H. Bell and CoolProp developers | [MIT license](https://raw.githubusercontent.com/CoolProp/CoolProp/v8.0.0/LICENSE) | `bundled_candidate` | `not_applicable_property_asset` | later independent evaluator: two in-range points, critical/triple boundaries and one out-of-range failure; preserve JSON/license identity | candidate only; no CoolProp runtime/Python dependency; no copied paper text; no accepted HUNDUN correlation pack yet |
| `n_dodecane_kerosene_surrogate` | `licensed_property_oracle` | NIST REFPROP 10 n-dodecane fluid | [official product](https://www.nist.gov/srd/refprop); [commercial terms](https://shop.nist.gov/ccrz__ProductDetails?cclcl=en_US&sku=SRD23DV10) | REFPROP `10.0`, 2018; DOI `10.18434/T4/1502528` | not computed: paid asset was not downloaded | DOI [10.18434/T4/1502528](https://doi.org/10.18434/T4/1502528); fluid database authors E. W. Lemmon, I. H. Bell, M. L. Huber and M. O. McLinden | U.S. Department of Commerce holds the SRD 23 copyright; NIST SRDP is the licensor; database authors are credited separately | [NIST SRD 23 site license](https://www.nist.gov/document/nist-srd-refprop-23-site-license-agreement); [product page](https://www.nist.gov/srd/refprop) | `user_supplied_only` | `not_applicable_property_asset` | user-licensed value-level cross-check only | no HUNDUN redistribution, fixture, source copy or required runtime dependency |
| `n_dodecane_kerosene_surrogate` | `property_compilation_reference` | NIST Chemistry WebBook SRD 69, CAS 112-40-3 | [official species page](https://webbook.nist.gov/cgi/cbook.cgi?ID=C112403&Mask=4) | continuously updated SRD 69 page; no immutable release supplied | unavailable: no immutable downloadable asset identity | cited primary references are listed on the species page | U.S. Secretary of Commerce on behalf of the U.S.A. | [page copyright/notes](https://webbook.nist.gov/cgi/cbook.cgi?ID=C112403&Mask=4) | `rejected` | `not_run_rejected_asset` | manual cross-check only | no bundled correlation; incomplete eight-property pack; no frozen revision/hash |
| `iso_octane_gasoline_surrogate` | `gas_mechanism_kinetics_thermo` | Osecky et al. 89-species reduced iso-octane mechanism | [ACS Figshare dataset](https://acs.figshare.com/articles/dataset/Investigation_of_Iso-octane_Ignition_and_Validation_of_a_Multizone_Modeling_Method_in_an_Ignition_Quality_Tester/3914538); [official file](https://ndownloader.figshare.com/files/6128106); [metadata API](https://api.figshare.com/v2/articles/3914538) | dataset DOI `10.1021/acs.energyfuels.6b01406.s001`, version `1`; posted `2016-08-18`; file header `LLNL-MI-421507`, review/release `2009-12-16` | `17ed9e579fccafff882c52134457b46d40fa0d953b2a783b0ee0a0216e220355` | Osecky et al., DOI [10.1021/acs.energyfuels.6b01406](https://doi.org/10.1021/acs.energyfuels.6b01406) | Osecky et al.; ACS repository; dataset metadata does not separately name a copyright holder | [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) as recorded by the official metadata API | `user_supplied_only` | `pending_later_external_ck2yaml`; one CHEMKIN text contains `ELEMENTS`, `SPECIES`, `THERMO ALL` and `REACTIONS`; no `TRANSPORT` section | future external conversion; verify 89 parsed species, reactions, units, thermo range and conversion output SHA | non-commercial license is not a general product redistribution grant; no gas transport, gasoline or ignition validation claim |
| `iso_octane_gasoline_surrogate` | `gas_transport_candidate` | LLNL `prf_tran_dat_v1b.txt` | [official transport file](https://combustion.llnl.gov/sites/combustion/files/prf_tran_dat_v1b.txt); [iso-octane v3 page](https://combustion.llnl.gov/mechanisms/alkanes/iso-octane-version-3) | transport `version 1b`; review/release `2004-05-19`; updated `2009-12-01` | `db1bbaa374153de3806e3e1dc8bb5325ddf7d2e3d981ced4218201678b629629` | Curran et al., DOI [10.1016/S0010-2180(01)00373-X](https://doi.org/10.1016/S0010-2180(01)00373-X) | H. J. Curran, W. J. Pitz, C. K. Westbrook and LLNL/LLNS; asset-specific holder/grant not established | [LLNL copyright status](https://www.llnl.gov/disclaimer) | `user_supplied_only` | `pending_later_external_ck2yaml_transport_coverage`; compatibility with every 89-species name is unproved | future external species-by-species coverage audit and conversion with the exact 89-species input | no standalone redistribution; no assumed pairing merely because both assets descend from LLNL PRF work |
| `iso_octane_gasoline_surrogate` | `rejected_detailed_mechanism_kinetics` | LLNL Iso-Octane Version 3 `ic8_ver3_mech.txt` | [kinetics](https://combustion.llnl.gov/sites/combustion/files/ic8_ver3_mech.txt); [official page](https://combustion.llnl.gov/mechanisms/alkanes/iso-octane-version-3) | Version 3; review/release `2009-12-16` | `2fc57d65d8dd7b4e9793f09639f41aca114aa878624f4b93b726ba2f39a68d5f` | Mehl et al. 2009; Curran et al., DOI [10.1016/S0010-2180(01)00373-X](https://doi.org/10.1016/S0010-2180(01)00373-X) | LLNL/LLNS and named authors; asset-specific grant absent | [LLNL copyright status](https://www.llnl.gov/disclaimer) | `rejected` | `not_run_rejected_asset` | none; this detailed kinetics file is not the low-cost interface candidate | no redistribution; detailed scale is outside the intended low-cost acceptance role |
| `iso_octane_gasoline_surrogate` | `rejected_detailed_mechanism_thermo` | LLNL Iso-Octane Version 3 `prf_v3_therm_dat.txt` | [thermo](https://combustion.llnl.gov/sites/combustion/files/prf_v3_therm_dat.txt); [official page](https://combustion.llnl.gov/mechanisms/alkanes/iso-octane-version-3) | Version 3; paired mechanism review/release `2009-12-16` | `acb98b7574d9568714b928444ebea8aae86e794168d1d597297dc2c00ff5ae4e` | Mehl et al. 2009; Curran et al., DOI [10.1016/S0010-2180(01)00373-X](https://doi.org/10.1016/S0010-2180(01)00373-X) | LLNL/LLNS and named authors; asset-specific grant absent | [LLNL copyright status](https://www.llnl.gov/disclaimer) | `rejected` | `not_run_rejected_asset` | none; this detailed thermo file is not part of the low-cost candidate | no redistribution; no standalone kinetics/transport capability or low-cost acceptance claim |
| `iso_octane_gasoline_surrogate` | `property_generator_and_tables` | FuelLib 3.0.0 GCM source/table candidate | [official NLR page](https://www.nlr.gov/computational-science/fuellib-jet-fuel-library); [PyPI 3.0.0](https://pypi.org/project/fuellib/3.0.0/); [source archive](https://files.pythonhosted.org/packages/7d/54/328cd9c99f842ec37d6887a6d55bb236eed59d14752050db24d355df19d2/fuellib-3.0.0.tar.gz) | FuelLib `3.0.0`; lightweight tag/commit `4e26169ccbc456f5f1e6b4f97d6832fa94192f91`; released `2026-06-25` | source archive `92633a64168878e2b3694c67fc6b66a670cebfd8e954be9e1de102c3cbe3ee02`; `gcmTable.csv` `a5fa6fb586f7cbcf4282d826dc1fd09f367919437d002f9c2b0d6eeeb7d43e8a`; `fuel.py` `cbeb754b4a79f774bf782dd497b7a30e4212f19ac09c27bb648cfe9a581aee4a` | software record DOI [10.11578/dc.20250317.1](https://doi.org/10.11578/dc.20250317.1); Constantinou--Gani and Govindaraju--Ihme references are listed upstream | Alliance for Energy Innovation, LLC | [BSD-3-Clause license](https://github.com/NatLabRockies/FuelLib/blob/v3.0.0/LICENSE) | `bundled_candidate` for the generator/table asset; final iso-octane pack remains `unresolved` | `not_applicable_property_asset` | maintainer-only generation after exact 2,2,4-trimethylpentane group decomposition is frozen; evaluate MW, Tb, Tc, density, viscosity, conductivity, cp, latent heat, surface tension and vapor pressure | no Python runtime; no claim that built-in `C08-Isoparaffin` is iso-octane—it is 2-methylheptane; no final pack/hash accepted yet |
| `iso_octane_gasoline_surrogate` | `partial_property_source` | NIST FDS 6.11.1 `Source/prop.f90` `ISOOCTANE` entry | [release](https://github.com/firemodels/fds/releases/tag/FDS-6.11.1); [frozen source](https://raw.githubusercontent.com/firemodels/fds/ff928dbf68a2623ced3a007be75391c6dc6c9f90/Source/prop.f90) | FDS `6.11.1`; peeled commit `ff928dbf68a2623ced3a007be75391c6dc6c9f90`; released `2026-07-10` | full `prop.f90` `ed9fb9deef75350d4c78fcead52611e8cb1931c34aba8942ea9959d715ab1b40` | FDS manuals and references embedded in the source | NIST-developed software is not subject to U.S. copyright; the notice does not establish provenance/holders for embedded property values | [NIST software notice](https://raw.githubusercontent.com/firemodels/fds/ff928dbf68a2623ced3a007be75391c6dc6c9f90/LICENSE.md) | `rejected`; the software notice alone is not an asset-level property-data grant | `not_run_rejected_asset` | none | source-internal identity/value conflicts and an incomplete property set make it unsafe even as a manual numerical oracle; no Fortran extraction |
| `iso_octane_gasoline_surrogate` | `licensed_property_oracle` | NIST REFPROP 10 iso-octane fluid | [official product](https://www.nist.gov/srd/refprop); [commercial terms](https://shop.nist.gov/ccrz__ProductDetails?cclcl=en_US&sku=SRD23DV10) | REFPROP `10.0`, 2018; DOI `10.18434/T4/1502528` | not computed: paid asset was not downloaded | DOI [10.18434/T4/1502528](https://doi.org/10.18434/T4/1502528); database authors are credited by the product | U.S. Department of Commerce holds the SRD 23 copyright; NIST SRDP is the licensor; database authors are credited separately | [NIST SRD 23 site license](https://www.nist.gov/document/nist-srd-refprop-23-site-license-agreement); [product page](https://www.nist.gov/srd/refprop) | `user_supplied_only` | `not_applicable_property_asset` | user-licensed value-level cross-check only | no redistribution, fixture, source copy or required runtime dependency |
| `iso_octane_gasoline_surrogate` | `property_compilation_reference` | NIST Chemistry WebBook SRD 69, CAS 540-84-1 | [official species page](https://webbook.nist.gov/cgi/cbook.cgi?ID=C540841&Mask=FFF) | continuously updated SRD 69 page; no immutable release supplied | unavailable: no immutable downloadable asset identity | primary references are listed on the species page | U.S. Secretary of Commerce on behalf of the U.S.A. | [page copyright/notes](https://webbook.nist.gov/cgi/cbook.cgi?ID=C540841&Mask=FFF) | `rejected` | `not_run_rejected_asset` | manual cross-check of identity, selected phase-change data and Antoine ranges only | compilation says all rights reserved; incomplete eight-property pack; no frozen revision/hash |

可冻结的许可证据也单独记录摘要，避免把“网页上看到许可证名称”误当成已固定的证据：

- Cantera v3.2.0 `License.txt` SHA-256 为
  `e92980b9712ce20e73898a97b0116889e84e07f548d6be8591e87dcad79c41bb`；同版本
  `data/README.md` SHA-256 为
  `067c221d01413542325424c7b9c6708674e67c4b5a346c58dfdea891ea81754f`。后者将示例数据
  定位为 illustrative material，并没有补出 Reitz 机理自身的再分发许可。
- CoolProp v8.0.0 `LICENSE` SHA-256 为
  `9bf835333ef602af4cb19338b9f9d43671e174fa029b00280e9bdba6ea4719b2`（MIT）。
- FuelLib v3.0.0 `LICENSE` SHA-256 为
  `2f912a6984c8115be8ba0b4535ae7cd50617743253148c4e2d4a07a7bddfb6fc`（BSD-3-Clause）。
- FDS commit `ff928dbf68a2623ced3a007be75391c6dc6c9f90` 的 `LICENSE.md` SHA-256 为
  `38c542304b97afc4171a9b67866499eaf222509cab45c095ea69bf88d57755b7`（NIST notice）。
- ACS Figshare version 1 的官方 metadata API 明确记录 CC BY-NC 4.0；许可 deed、LLNL
  disclaimer、NIST 商店条款和 WebBook 页面是动态网页，不伪造“冻结网页 SHA”。缺少资产级
  明确许可时，状态保持 `user_supplied_only` 或 `rejected`。
- REFPROP 的官方 SRD 23 site license 说明 U.S. Department of Commerce 取得数据库版权
  保护、NIST SRDP 是许可方，并限制复制与分发；数据库作者是技术署名，不替代法律权利主体。
- FDS notice 允许保留 notice 地使用、复制和分发 NIST-developed software，并说明 NIST
  雇员开发的软件不受美国版权保护；它没有单独证明 `prop.f90` 内嵌 property values 的来源、
  第三方权利或可作为 HUNDUN 数据资产再分发的范围。

## 3. 格式、范围和预期消费

> 历史说明：本节第 3.1、3.2 节保留初始 legal subcluster 在执行 Step 4 前记录的
> `unresolved` 字段状态。当前验证结论以第 5.1--5.4 节为准；保留原文是为了不回写或抹去
> P0-4 的决策历史。

### 3.1 n-dodecane 链

- `nDodecane_Reitz.yaml` 是已经转换为 Cantera YAML 的 reduced n-dodecane/PAH 机理；文件
  自述为 100 species、432 reactions，并给出 ideal-gas `nDodecane_IG` 和
  Redlich--Kwong `nDodecane_RK` 两个 phase。冻结文件中没有 species transport block，
  所以“YAML 可读”和“Stage 4 mixture-averaged transport 可用”是两个不同结论。
- LLNL Hybrid 是 CHEMKIN kinetics 与独立 thermo 文件；官方页面给出 650--1400 K、
  1--100 bar、`phi=0.25--3` 的模型范围。它是低成本备选，不补齐 transport 许可或格式。
- CoolProp JSON 覆盖 molecular weight、critical/triple state、EOS、饱和液密度、粘度、
  导热、液体热容、汽化潜热、表面张力和饱和蒸气压所需的模型/ancillary。各相关项的原生
  范围并不相同；下面分别记录，不能把 EOS 上限直接写成每个 liquid property 的接受范围。

| n-dodecane property | candidate source/key | intended pack SI unit / conversion | native or published validity evidence | P0-4 field status |
|---|---|---|---|---|
| critical temperature | CoolProp JSON `STATES.critical.T` and `EOS[0].STATES.reducing.T` | K | identity constant `658.1 K`; not a temperature interval | `unresolved`: source asset is `bundled_candidate`, but no pack field has passed Step 4 |
| normal boiling point | [NIST WebBook CAS 112-40-3](https://webbook.nist.gov/cgi/cbook.cgi?ID=C112403&Mask=4), `Tboil = 489 +/- 2 K`; CoolProp JSON has no standalone boiling-point input | K | dynamic compilation page, not an immutable correlation asset; a later EOS solve at `101325 Pa` would be a derived check | `unresolved`; WebBook asset remains `rejected` for bundling |
| molecular weight | CoolProp JSON `EOS[0].molar_mass` | kg/mol | constant `0.17033484 kg/mol`; no temperature interval | `unresolved`: source field frozen, pack acceptance deferred |
| saturated-liquid density | CoolProp JSON `ANCILLARIES.rhoL` | raw mol/m^3; convert with the frozen molar mass to kg/m^3 | JSON `Tmin=263.6 K`, `Tmax=658.0999999999985 K` | `unresolved`: evaluator, conversion and endpoint policy are Step 4 work |
| dynamic viscosity | CoolProp JSON `TRANSPORT.viscosity`; Huber--Laesecke--Perkins correlation | Pa s | publication range `263.59--800 K`, up to `200 MPa`; the same JSON EOS has `T_max=700 K`, `p_max=200 MPa`, so the common evaluator intersection is not yet frozen | `unresolved`: state convention and accepted intersection are not frozen |
| thermal conductivity | CoolProp JSON `TRANSPORT.conductivity`; Huber--Laesecke--Perkins correlation | W/(m K) | publication range `263.59--800 K`, up to `200 MPa`; the same JSON EOS has `T_max=700 K`, `p_max=200 MPa`, so the common evaluator intersection is not yet frozen | `unresolved`: state convention and accepted intersection are not frozen |
| liquid heat capacity | derivative of CoolProp JSON `EOS[0]` at a later explicitly frozen liquid state | J/(kg K) | JSON EOS declares `Ttriple=263.6 K`, `T_max=700 K`, `p_max=200 MPa`; that does not itself choose saturated versus fixed-pressure liquid | `unresolved`: evaluator/state and numerical evidence are deferred |
| latent heat of vaporization | CoolProp JSON `ANCILLARIES.hLV` | raw J/mol; divide by frozen molar mass for J/kg | JSON `Tmin=263.6 K`, `Tmax=658.0 K` | `unresolved`: evaluator, conversion and upper-end failure policy are deferred |
| surface tension | CoolProp JSON `ANCILLARIES.surface_tension`, Mulero et al. | N/m | entry supplies `Tc=658.1 K` but no native `Tmin`; a conservative lower bound cannot be attributed to this entry without a later HUNDUN gate | `unresolved`: lower bound and out-of-range behavior are not frozen |
| saturation vapor pressure | CoolProp JSON `ANCILLARIES.pS` | Pa | JSON `Tmin=263.6 K`, `Tmax=658.0999999999985 K` | `unresolved`: evaluator and endpoint failure policy are deferred |

这里的数值只是冻结源字段或上游声明范围，不是属性求值结果。尤其是 `263.6 K <= T <
658.1 K` 只能作为后续 consumed-field gate 的候选交集；P0-4 Step 4 尚未把它接受为
HUNDUN native range。

预期消费方式是：机制始终由用户在外部 hash-verified 路径提供；CoolProp JSON 只有在
notice、consumed-field 和数值 gate 通过后才可生成/成为 HUNDUN 的 typed property pack
候选。普通用户运行时不引入 CoolProp 或 Python。

### 3.2 iso-octane 链

- ACS Figshare 文件是一个 CHEMKIN 文本，含 kinetics 和 thermo，没有 transport；出版者
  说明它是从 874-species 详细机理降至 89 species 的接口成本候选。CC BY-NC 4.0 不适合
  HUNDUN 的不受非商业条款限制的产品分发，因此必须外置。
- LLNL PRF transport 文件是独立资产。只有未来逐 species 名称验证 89-species 文件全部被
  覆盖后，才可把两者作为一次 `ck2yaml` 输入组合；共同学术谱系不等于格式兼容或共同许可。
- FuelLib 3.0.0 的 BSD 资产提供 GCM tables，以及 MW、Tb、Tc、density、Dutt viscosity、
  liquid cp、Lee--Kesler/Ambrose--Walton vapor pressure、latent heat、Brock--Bird/Pitzer
  surface tension 和 Latini conductivity 的计算路径。它的 built-in `C08-Isoparaffin`
  明确代表 2-methylheptane，不是 2,2,4-trimethylpentane；因此本轮只接受 generator/table
  作为 `bundled_candidate`，不虚构 iso-octane output pack。
- FuelLib 源证据精确到 sdist 内 `fuellib/data/gcmTableData/gcmTable.csv` 与
  `fuellib/fuel.py`。v3.0.0 的这些方法没有统一的 native temperature-domain enforcement 或
  “out of range” 状态契约；所以文献关联名不能替代逐属性接受范围。

| iso-octane property | FuelLib v3.0.0 source/method | intended pack SI unit | native validity behavior in frozen code | P0-4 field status |
|---|---|---|---|---|
| critical temperature | GCM table plus exact group decomposition loaded by `fuel.__init__` | K | generated constant; no temperature domain, and the 2,2,4-trimethylpentane decomposition is not frozen | `unresolved` |
| normal boiling point | GCM table plus exact group decomposition loaded by `fuel.__init__` | K | generated constant; no temperature domain, and no output value/hash exists | `unresolved` |
| molecular weight | GCM table plus exact group decomposition loaded by `fuel.__init__` | kg/mol | generated constant; no temperature domain, and no output value/hash exists | `unresolved` |
| density | `fuel.density` through `fuel.molar_liquid_vol` | kg/m^3 | no general range rejection; `molar_liquid_vol` changes branch above `Tc` instead of establishing a pack failure boundary | `unresolved` |
| viscosity | `fuel.viscosity_kinematic` and `fuel.viscosity_dynamic` (Dutt path) | m^2/s and Pa s; pack would consume Pa s | no explicit temperature-range enforcement | `unresolved` |
| thermal conductivity | `fuel.thermal_conductivity` (Latini path) | W/(m K) | no explicit temperature-range enforcement; reduced-temperature powers do not provide a HUNDUN error contract | `unresolved` |
| liquid heat capacity | `fuel.Cp` and `fuel.Cl` | `Cp`: J/(mol K); `Cl`: J/(kg K) | polynomial evaluation has no explicit temperature-range enforcement | `unresolved` |
| latent heat of vaporization | `fuel.latent_heat_vaporization` | J/kg | no explicit range rejection; the frozen method returns zero above `Tc` rather than a HUNDUN out-of-range error | `unresolved` |
| surface tension | `fuel.surface_tension` (Brock--Bird or Pitzer) | N/m | no explicit range rejection; reduced-temperature powers do not provide a HUNDUN error contract | `unresolved` |
| saturation vapor pressure | `fuel.psat` (Lee--Kesler or Ambrose--Walton) | Pa | no explicit temperature-range enforcement | `unresolved` |

因此 FuelLib 的 exact generator/table files 可以保持 `bundled_candidate`，但每个 iso-octane
output field 都保持 `unresolved`。后续若接受 pack，必须由 HUNDUN 自己冻结并执行输入范围与
越界失败规则，不能把 FuelLib 当前的外推、分支、零值或非有限结果当成原生验证。

- FDS 的固定 `ISOOCTANE` block 不能用作可选数值 oracle：它标记 `FORMULA='C8H18'`，却以
  `CP_UNIT=8314.472/112.21264` 换算，并给出 `T_BOIL=398.44 K`；NIST WebBook 对同一 CAS
  记录的 molecular weight 为 `114.2285`、boiling point 为 `372.4 +/- 0.2 K`。同一
  `prop.f90` 的 `N-OCTANE` block 又复用 `398.44 K` 及 density/conductivity/viscosity
  常数。加上缺失完整 Tc、surface-tension 和 vapor-pressure package，这构成身份/数值一致性
  冲突；本轮不提取、不采用、也不拿其中任何值作 manual cross-check。

预期消费方式是：89-species/transport 均在外部由用户提供；FuelLib 只允许在 maintainer
环境中生成一个具有独立输入、命令、输出 SHA 和 BSD notice 的候选 pack。HUNDUN 正常
configure/build/install/runtime 不依赖 FuelLib、Python、NIST WebBook 或 REFPROP。

## 4. 选择、拒绝和未解决矩阵

| Family | Provisional low-cost mechanism | Gas transport | Pure-liquid properties | End-to-end result |
|---|---|---|---|---|
| n-dodecane / kerosene surrogate | primary: Cantera Reitz YAML `user_supplied_only`; fallback: LLNL Hybrid `user_supplied_only` | `unresolved`: no legally and technically complete frozen transport candidate | CoolProp v8.0.0 JSON source `bundled_candidate`; every output pack field remains `unresolved` | `PREFLIGHT_PARTIAL`; no bundle-ready chain |
| iso-octane / gasoline surrogate | ACS Figshare 89-species CHEMKIN+thermo `user_supplied_only` | LLNL PRF transport `user_supplied_only`, exact coverage unresolved | FuelLib generator/table `bundled_candidate`; exact iso-octane output pack `unresolved` | `PREFLIGHT_PARTIAL`; no bundle-ready chain |

Rejected alternatives include the LLNL detailed Iso-Octane Version 3 for the low-cost role, NIST
WebBook compilations as immutable redistributable packs, and the internally inconsistent/incomplete FDS
iso-octane entry as either a standalone pack or a numerical oracle. GRI-Mech is not a dual-fuel candidate
and supplies no evidence for either surrogate.

The important P0 result is that no public chain is presently safe to redistribute end-to-end. A legally
redistributable property component does not repair a non-redistributable mechanism or unresolved gas
transport asset.

## 5. P0-3 依赖和 P0-4 Step 4 后续条件

> 本节第 1--7 项是 commit `3ed7f79` 冻结的执行前条件，现已按第 5.1--5.4 节逐项处置。
> 它们不再是待办清单；未解决的 gas transport、再分发许可和产品接入仍按原状态保留。

P0-3 must not consume any fuel mechanism marked `user_supplied_only`. Its standalone Cantera linkage
smoke must instead use a HUNDUN-authored synthetic tiny mechanism with a separately recorded explicit
license, source SHA-256, phase and expected counts. This preserves P0-3 progress without copying one of
the fuel mechanisms into Git, the artifact or the install tree.

The fuel-specific checks below are Task P0-4 Step 4 work. Under the live plan they require an accepted
P0-3 first, remain pending for a later legal external-asset gate, and were not run in this public-source
subcluster. Therefore the whole Task P0-4 remains incomplete and requires a separate later validation
commit; these checks are not prerequisites for the synthetic P0-3 smoke:

1. Re-fetch each selected external asset by its exact official URL and require the SHA-256 above before
   any conversion or parse.
2. For Reitz YAML, require the requested phase, exact parsed species/reaction/element identities,
   thermo ranges, duplicate handling and an explicit `transport_model=none` result unless a separately
   licensed complete transport asset is added.
3. For LLNL Hybrid and the ACS 89-species candidate, record exact `ck2yaml` command/version, input
   order, diagnostics, output YAML SHA, phase, parsed counts and all unit conversions. A conversion
   warning or missing thermo/transport row is a rejection, not a silent default.
4. Before pairing ACS 89-species kinetics with LLNL PRF transport, compare every canonical species name
   and reject missing, duplicate or alias-only matches. Record any later alias mapping as a reviewed
   external manifest, not product hard-coding.
5. For CoolProp JSON, freeze a minimal consumed-field set and independently evaluate at least two
   in-range temperatures plus one out-of-range failure. Require finite positive density, viscosity,
   conductivity and cp; positive latent heat/surface tension below Tc; and explicit failure rather than
   clamp/extrapolation outside the accepted range.
6. For FuelLib, first freeze a reviewed exact 2,2,4-trimethylpentane group decomposition. Generate a
   candidate pack in a maintainer-only environment, record FuelLib/archive/table/code/input/output
   hashes, then run the same two-in-range/one-out-of-range checks. Until then, the final iso-octane
   property pack remains `unresolved`.
7. Cross-check a small number of values against NIST WebBook or a user-licensed REFPROP installation at
   value level only. Do not copy their database pages/files into fixtures or make them runtime
   dependencies.

Any later COAST `EXEC/Fuels` comparison is outside this P0-4 public audit. It still requires the user to
confirm the exact current realpath and version before any read; no candidate path may be guessed from
this document.

### 5.1 外部机制解析结果

全部资产重新从登记的官方 HTTPS URL 获取并按初始 SHA-256 拒绝漂移。CHEMKIN 转换使用
Cantera 3.2.0 的 `ck2yaml.py`、Python 3.10.12、NumPy 1.26.4、ruamel.yaml 0.17.21，命令、
输入顺序、构建器 wheel 和输出 YAML 均在外部 manifest 中绑定。转换没有 warning 或
error；由于最终验证由独立 C++ `Solution` 解析完成，转换命令使用 `--no-validate`，不把
Python 验证器当作产品或接受依赖。

| candidate | phase | C++ parsed species / reactions / elements | transport result | disposition |
|---|---|---:|---|---|
| Cantera Reitz n-dodecane | `nDodecane_IG` | `100 / 553 / 4` | `none` | parser interface passed; external user-supplied mechanism only |
| LLNL n-dodecane Hybrid | `gas` | `65 / 363 / 6` | `none` | strict conversion and parser interface passed; user-supplied only |
| ACS 89-species iso-octane | `gas` | `89 / 480 / 6` | `none` | strict conversion and parser interface passed; user-supplied only |

Reitz 文件说明文字中的 `432 reactions` 与冻结 YAML 实际列表不一致；词法计数和 Cantera C++
解析均得到 `553`，所以这里以可执行 phase 的解析结果为权威，不修改文件、不强行匹配说明文字。
三份输出在 `800 K / 101325 Pa` 得到有限、正的 density 和 cp。错误 species count mutation
被拒绝。ACS 89-species 列表在 LLNL PRF transport 中没有 missing name，但存在三个命中候选
各有重复记录：`NEOC5H10OOH`、`NEOC5H11`、`NEOC5H11O2`。因此 transport pairing 明确
`rejected`；没有使用 permissive 模式、alias 或任意选第一条记录。

关键证据：

- refetch script SHA-256：`7c1fd34951c6371f072d63f2bd84e853ac8056b73c1cc7c22bf49e509562214f`；
- conversion script SHA-256：`cd2662038daad8e430b7554375da571b02f6205eb010240b1cb93cfc249d4c5e`；
- C++ inspector SHA-256：`50b6618ebd9a9e024f58539649547bc7dc97b0c5e22075bf96a3bb0997621d04`；
- conversion evidence manifest SHA-256：`85e8ba9dc41acb38995723f81e5df7e33f359404c9bc5e93e83706dfdae621be`；
- parse evidence manifest SHA-256：`79c6471b414225ebdf1f18945ba740698cef461854d24a191c4257ab25550447`。

### 5.2 双液体候选 pack 和范围契约

维护者环境固定 CoolProp 8.0.0、FuelLib 3.0.0、NumPy 1.26.4、pandas 2.2.3、SciPy 1.15.3
及其纯 Python 传递 wheel。它们只生成外部候选值；HUNDUN 正常 configure、build、test、
install 和 runtime 不消费 Python、FuelLib 或 CoolProp。每个 builder wheel 的精确 URL、版本、
SHA-256 和 wheel 内 principal license 已冻结；license evidence manifest SHA-256 为
`4ae36a3033ec01fe8c200f55e009045d6a458bc9ac7b3cf0c9abf5349bad5129`。这些 builder
dependencies 不进入安装包，也不转化为 HUNDUN 的运行时传递依赖。

CoolProp n-dodecane 候选以饱和液体状态求值，候选范围为
`263.6 K <= T < 658.10002685820268 K`，在 `300 K` 和 `450 K` 生成 density、dynamic
viscosity、thermal conductivity、mass cp、latent heat、surface tension 和 saturation
pressure，`700 K` 由范围 wrapper 明确拒绝。候选 pack SHA-256 为
`15f6612239f852148e23e4120f951f57e12fda6733a40df8c2f2f983b4cc8c5f`。

iso-octane 候选使用 PubChem CID 10907 的公开 connectivity
`CC(C)CC(C)(C)C`，独立冻结一阶 groups `CH3=5, CH2=1, CH=1, C=1` 和二阶 groups
`(CH3)2CH=1, (CH3)3C=1`；它重构 `C8H18`，没有把 FuelLib 内置的 2-methylheptane
当作 iso-octane。decomposition CSV SHA-256 为
`c514a906f94926a7320a3db6316e4976967035b3df05edcfa04bd904b4e0dac9`。FuelLib 候选范围为
`298 K <= T < 540.33354268219875 K`，在 `300 K` 和 `350 K` 生成同一七字段，
`Tc + 1 K` 由范围 wrapper 明确拒绝。候选 pack SHA-256 为
`d8018ce92b0a6384a636ad7a15d8c89486a0b6f25fcff63bf799bc67285a2238`。

这两个 pack 是 future typed-pack 输入候选，不是最终 HUNDUN property correlation，也不表示
两套上游方法在全范围内准确。它们的作用只是在两个性质明显不同的纯液体上冻结相同字段、
单位和失败协议。

### 5.3 通用 C++ 消费与 mutation RED

standalone C++17 检查器只按 schema、SHA-256 identity、范围和 SI 数值字段消费 pack，不按
`n-dodecane` 或 `iso-octane` 名称分派。它要求两点都位于半开范围内、所有七个属性有限且
为正、一个范围外温度的状态恰为 `rejected`，并要求两个 family 和 fingerprint 不同。

GREEN 二进制 SHA-256 为
`80929478f44734c12cbc0c2859050e84d87941d8f4ce8532d1df9d5b3582bd83`。以下 mutations 均以
exit `65` 被拒绝：

- density 改为 `nan`；
- 范围外结果从 `rejected` 改为 `accepted`；
- 两个 pack 使用同一 fingerprint；
- 把一个有效点的温度移到接受范围外。

第一次 GREEN runner 已正确生成 pack 并通过通用检查，但其 mutation 日志匹配器错误地等待
下游“non-positive”消息；`nan` 实际在更早的 finite parser 被拒绝。该次 status `1` 和日志被
保留，v2 只修正测试判据并通过。v2 primary manifest SHA-256 为
`376b73b96fbe5e0a40206488853dfcce1b40c8b1fb3a0b0ed5913946e936e9a6`，evidence manifest
SHA-256 为 `e6f39829354b7b15ba5062fea532547a3e1e4e7a64fb3124a02b47e88ba9241a`。

### 5.4 完成处置

P0-4 Step 4 的 parser/property interface gate 已完成，且两种燃料身份没有产品硬编码。整个
候选链仍为 `PREFLIGHT_PARTIAL`：三份机制继续 `user_supplied_only`，PRF transport pairing
因重复记录被拒绝，两个 property pack 只是外部候选。NIST WebBook/REFPROP 没有复制进
fixture，也没有成为依赖。动态 WebBook 页面只在临时目录读取，源页面没有保留；派生 summary
以 molecular weight `1%`、normal boiling point `2%` 的预研数量级容差比较两个 pack，均通过。
该容差不是产品科学阈值，也不覆盖其他七字段。summary SHA-256 为
`15f14afdd74ad35dd9256ee6a2f746e362ea78bd99a0e8b9977931fba492295e`，evidence manifest
SHA-256 为 `a8b1131de75e536d8330f93195e230e15ad38ffd368486940f31b2b26a081288`；它不提升任何资产的
再分发状态。任何正式 Stage 4/6 消费都必须在 Stage 3 接受后的 intake 中重新绑定
accepted product HEAD、schema、notice、consumed fields 和产品测试。

## 6. Capability limitations

The two families can eventually test only that the parser/property interface is not hard-coded to one
fuel identity. Even after the pending low-cost checks pass, they do not establish:

- real aviation kerosene or real gasoline fidelity;
- ignition-delay, flame-speed, soot or emissions accuracy;
- spray atomization, evaporation or reacting-spray validation;
- COAST equivalence or similarity;
- HUNDUN Stage 4, Stage 6 or v1 scientific acceptance.

No Vblowoff, Flame D, 48^3/96^3 case, long chemistry run or coupling implementation belongs to this
candidate audit.
