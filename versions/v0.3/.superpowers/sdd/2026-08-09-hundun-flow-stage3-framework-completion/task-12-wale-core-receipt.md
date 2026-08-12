# Task 12：backend-neutral WALE core 实施回执

基线：`da61bf945b2ac313dea84fb9fdfb7f42abb68a8b`

状态：`TASK_GATE_ACCEPTED`

## 已实现边界

- 新增独立 `hundun_les` 静态库以及公开 `les_wale.hpp`；
- `WaleModel` 在构造期冻结 owned-first active layout、背景单元体积和配置指纹；
- `evaluate` 只读 lagged velocity gradient 与本次 attempt density，生成由对象拥有的
  `nu_t` 和 `mu_sgs`，不写 committed flow state；
- WALE 使用 `g g` 的对称无迹部分；以最大梯度分量作齐次归一化，避免中间不变量
  overflow/underflow，同时保持速度尺度一次齐次和长度尺度二次齐次；
- 零梯度和双零不变量返回 bitwise `+0.0`，没有 epsilon、clip 或经验滤波；
- 所有 active index 在进入 unchecked kernel view 前验证 ghost 覆盖，stale checked
  view、非正密度、非有限梯度和非 host backend 均被拒绝；
- coefficient identity 由模型、step/dt/order、committed/history、lagged gradient 和
  density fingerprint 派生，不伪装成 committed revision；
- 测试专用数学入口只在 `HUNDUN_BUILD_TESTS=ON` 编译，不进入 tests-off 产品库。

## RED 与 mutation 覆盖

- `g g` 变成 `g g^T`、遗漏对称或无迹投影；
- 分母指数、速度尺度或 filter-width 平方错误；
- 零场返回负零或固定 epsilon；
- `mu_sgs` 使用 stale density；
- tensor 旋转不变性和 `y^3` 近壁律被破坏；
- active order、metadata 或输入 fingerprint mutation 不改变 identity；
- stale view、非有限输入或 device context 被接受。

## focused 证据

- `test_wale_header_contract`：PASS；
- `test_wale` tensor oracle、旋转、尺度、exact-zero、`y^3` 与非有限 RED：PASS；
- `test_wale_mpi_1_rank`：PASS；
- `test_wale_mpi_2_rank`：PASS；
- Debug 测试二进制 SHA-256：
  - `test_wale`：`52cc95b8e440bce1adb7815fbc7bd9d9f12422e3b874cee99690798b1ca00ade`；
  - `test_wale_header_contract`：`07e03021425e82d7293193008793721c5728c4827c707063d610418edd1b62c7`；
  - `test_wale_mpi`：`df912ccbbde872d8671d3fa03dfdf25ae8ea7551ddc5521962d47414650ebddf`；
- tests-off `libhundun_les.a` SHA-256：
  `119733f4c7b24625ef3bc60f9a371617803da67c13813f1dadf52d70e8fd24a6`；
- tests-off `nm -C` 不含 `wale_kinematic_viscosity_for_test`：PASS；
- tests-off install 包含公开头和 `libhundun_les.a`：PASS；
- `git diff --check`：PASS。

## 延期与调用方

Task 12 不接入 flow composition，也不构造 lagged velocity gradient。body-fitted
variable-viscosity、`none/wale` driver、一次求值多消费者、TGV 12/24 screen 和
focused sanitizer 按批准计划进入 Task 13+19B 与 WALE milestone。没有运行 48³、
96³或大型 MPI 数值矩阵。

## 审查与版权

主 agent 完成 WALE 数学、量纲、public API、ownership/lifetime、checked-view、
allocation、tests-off/install、调用方影响和完整 task diff 审查。当前 Codex model
catalog 警告禁止重启前设置子代理 model/reasoning override，因此没有伪造 Luna
review 证据。

公式来自已批准规格和公开论文的数学定义；实现、命名、控制流和独立 oracle 均为
HUNDUN-FLOW 原创代码，未复制外部或私有源码，未新增运行时依赖。
