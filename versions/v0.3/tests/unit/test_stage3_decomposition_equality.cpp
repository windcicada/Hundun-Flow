// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/stage3_decomposition_equality.hpp"
#include "tests/support/test_main.hpp"

namespace {

void run() {
  HUNDUN_CHECK(hundun::test::stage3::
                   decomposition_equality_oracle_is_mutation_sensitive());
}

} // namespace

int main() { return hundun::test::run(run); }
