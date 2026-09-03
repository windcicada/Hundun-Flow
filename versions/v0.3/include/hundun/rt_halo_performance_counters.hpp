// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstdint>

namespace hundun::runtime {

struct HaloPerformanceCounters final {
  std::uint64_t completed_exchanges{};
  std::uint64_t begin_calls{};
  std::uint64_t wait_calls{};
  std::uint64_t send_payload_bytes{};
  std::uint64_t receive_payload_bytes{};
  std::uint64_t pack_bytes{};
  std::uint64_t unpack_bytes{};
  std::uint64_t send_messages{};
  std::uint64_t receive_messages{};
  double completed_wait_seconds{};
};

}  // namespace hundun::runtime
