// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/lin_vector_ops.hpp"

#include <type_traits>

static_assert(std::is_final_v<hundun::linear::VectorOps>);
static_assert(std::is_aggregate_v<hundun::linear::DotProductPair>);

int main() { return 0; }
