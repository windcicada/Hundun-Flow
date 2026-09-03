// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_status.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <type_traits>

namespace {

using hundun::v04::Status;
using hundun::v04::StatusCode;
using hundun::v04::status_message;

struct StatusExpectation {
  StatusCode code;
  std::string_view message;
};

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

}  // namespace

int main() {
  static_assert(std::is_same_v<hundun::v04::FieldId, std::uint16_t>);
  static_assert(std::is_same_v<hundun::v04::StageId, std::uint16_t>);
  static_assert(std::is_same_v<hundun::v04::RevisionToken, std::uint64_t>);
  static_assert(std::is_same_v<hundun::v04::PlanFingerprint, std::uint64_t>);

  const hundun::v04::Int3 integer_point{1, 2, 3};
  const hundun::v04::Real3 real_point{1.0, 2.0, 3.0};
  int values[] = {4, 5};
  const hundun::v04::Span<int> span{values, 2};

  bool passed = true;
  passed &= expect(integer_point.x == 1 && integer_point.y == 2 &&
                       integer_point.z == 3,
                   "Int3 preserves its three components");
  passed &= expect(real_point.x == 1.0 && real_point.y == 2.0 &&
                       real_point.z == 3.0,
                   "Real3 preserves its three components");
  passed &= expect(span.data == values && span.size == 2,
                   "Span preserves pointer and size");
  passed &= expect(static_cast<bool>(Status{}),
                   "a default Status represents success");
  passed &= expect(!static_cast<bool>(Status{StatusCode::invalid_plan, 9}),
                   "a non-ok Status represents failure");

  constexpr std::array expectations{
      StatusExpectation{StatusCode::ok, "ok"},
      StatusExpectation{StatusCode::invalid_case, "invalid case"},
      StatusExpectation{StatusCode::invalid_plan, "invalid plan"},
      StatusExpectation{StatusCode::allocation_failure, "allocation failure"},
      StatusExpectation{StatusCode::mpi_failure, "MPI failure"},
      StatusExpectation{StatusCode::numerical_failure, "numerical failure"},
      StatusExpectation{StatusCode::rejected_step, "rejected step"},
      StatusExpectation{StatusCode::io_failure, "I/O failure"},
  };

  for (const auto& expectation : expectations) {
    passed &= expect(status_message(Status{expectation.code, 37}) ==
                         expectation.message,
                     expectation.message);
  }
  passed &= expect(status_message(Status{static_cast<StatusCode>(65535), 0}) ==
                       "unknown status",
                   "an unknown status code has a safe message");

  return passed ? 0 : 1;
}
