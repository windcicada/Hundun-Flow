// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <array>

namespace hundun::les::test {

double wale_kinematic_viscosity_for_test(const std::array<double, 9> &gradient,
                                         double coefficient,
                                         double filter_width_m);

} // namespace hundun::les::test
