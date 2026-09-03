// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "src/app_dispatch_order_detail.hpp"
#include "src/app_passive_scalar_driver_detail.hpp"

#include "tests/support/test_main.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace {

using EstablishRoot = std::filesystem::path (*)(
    const hundun::runtime::MpiContext&, const std::filesystem::path&);
using RunPassiveScalar = int (*)(
    const hundun::application::CliOptions&, hundun::runtime::MpiContext&,
    const hundun::config::CaseConfig&, const std::filesystem::path&);

static_assert(std::is_same_v<
              decltype(&hundun::application::establish_authoritative_case_root),
              EstablishRoot>);
static_assert(std::is_same_v<
              decltype(static_cast<RunPassiveScalar>(
                  &hundun::application::run_passive_scalar_case)),
              RunPassiveScalar>);

class PhaseFailure final : public std::runtime_error {
 public:
  explicit PhaseFailure(const std::string& message)
      : std::runtime_error(message) {}
};

enum class FailingPhase { none, root, config, rank };

struct DispatchProbeResult {
  std::vector<std::string> phases;
  std::string failure;
  int result{};
};

DispatchProbeResult run_probe(FailingPhase failing_phase) {
  DispatchProbeResult probe;
  try {
    probe.result = hundun::application::detail::dispatch_in_root_config_rank_order(
        [&]() -> std::filesystem::path {
          probe.phases.emplace_back("root");
          if (failing_phase == FailingPhase::root) {
            throw PhaseFailure("root-first-failure");
          }
          return "/tmp/case";
        },
        [&]() -> hundun::config::ResolvedCase {
          probe.phases.emplace_back("config");
          if (failing_phase == FailingPhase::config) {
            throw PhaseFailure("config-first-failure");
          }
          return hundun::config::CaseConfig{};
        },
        [&](const hundun::config::ResolvedCase& resolved,
            const std::filesystem::path&) -> int {
          probe.phases.emplace_back("rank");
          HUNDUN_CHECK(
              std::holds_alternative<hundun::config::CaseConfig>(resolved));
          if (failing_phase == FailingPhase::rank) {
            throw PhaseFailure("rank-first-failure");
          }
          return 17;
        });
  } catch (const PhaseFailure& error) {
    probe.failure = error.what();
  }
  return probe;
}

void check_phase_order_and_first_failure() {
  const DispatchProbeResult success = run_probe(FailingPhase::none);
  HUNDUN_CHECK(success.phases ==
               (std::vector<std::string>{"root", "config", "rank"}));
  HUNDUN_CHECK(success.failure.empty());
  HUNDUN_CHECK(success.result == 17);

  const DispatchProbeResult root_failure = run_probe(FailingPhase::root);
  HUNDUN_CHECK(root_failure.phases ==
               (std::vector<std::string>{"root"}));
  HUNDUN_CHECK(root_failure.failure == "root-first-failure");

  const DispatchProbeResult config_failure = run_probe(FailingPhase::config);
  HUNDUN_CHECK(config_failure.phases ==
               (std::vector<std::string>{"root", "config"}));
  HUNDUN_CHECK(config_failure.failure == "config-first-failure");

  const DispatchProbeResult rank_failure = run_probe(FailingPhase::rank);
  HUNDUN_CHECK(rank_failure.phases ==
               (std::vector<std::string>{"root", "config", "rank"}));
  HUNDUN_CHECK(rank_failure.failure == "rank-first-failure");
}

}  // namespace

int main() {
  return hundun::test::run([] { check_phase_order_and_first_failure(); });
}
