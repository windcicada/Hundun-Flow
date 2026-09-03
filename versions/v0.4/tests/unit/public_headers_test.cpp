// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_app.hpp"
#include "hundun/v04_boundary.hpp"
#include "hundun/v04_case.hpp"
#include "hundun/v04_execution.hpp"
#include "hundun/v04_field.hpp"
#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"
#include "hundun/v04_initialization.hpp"
#include "hundun/v04_io.hpp"
#include "hundun/v04_linear.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_parallel.hpp"
#include "hundun/v04_physics.hpp"
#include "hundun/v04_product.hpp"
#include "hundun/v04_status.hpp"
#include "hundun/v04_types.hpp"

#include <type_traits>

static_assert(
    std::is_nothrow_move_constructible<hundun::v04::CompiledCasePlan>::value,
    "public product plan must remain cheaply movable");
static_assert(
    std::is_nothrow_move_constructible<hundun::v04::ProductDriver>::value,
    "public driver must remain cheaply movable");

int main() { return 0; }
