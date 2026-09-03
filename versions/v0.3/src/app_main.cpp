// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "app_cli_options_detail.hpp"
#include "app_dispatch_order_detail.hpp"
#include "app_flow_driver_detail.hpp"
#include "app_immersed_flow_driver_detail.hpp"
#include "app_passive_scalar_driver_detail.hpp"
#include "app_reacting_flow_driver_detail.hpp"
#include "app_resolved_case_v3_broadcast_detail.hpp"
#include "app_resolved_case_v4_broadcast_detail.hpp"
#include "app_version.hpp"

#include "hundun/cfg_resolved_case_v3_loader.hpp"
#include "hundun/cfg_resolved_case_v4_loader.hpp"
#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "rt_mpi_error_detail.hpp"

#include <mpi.h>

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <locale>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

namespace {

using hundun::runtime::Error;
using hundun::runtime::MpiContext;

enum class RunMode : int { normal = 0, validate = 1, print_resolved = 2 };

std::string exception_message_or_fallback(const std::exception& error,
                                          std::string_view fallback) {
  const char* message = error.what();
  return message != nullptr && *message != '\0' ? std::string(message)
                                                : std::string(fallback);
}

void require_collective_success(const MpiContext& context, bool local_ok,
                                std::string_view local_message) {
  const auto status =
      hundun::runtime::collective_status(context, local_ok, local_message);
  if (!status.ok) {
    throw Error(status.message);
  }
}

RunMode run_mode(const hundun::application::CliOptions& options) noexcept {
  if (options.validate_only) {
    return RunMode::validate;
  }
  if (options.print_resolved) {
    return RunMode::print_resolved;
  }
  return RunMode::normal;
}

void require_run_mode_agreement(const MpiContext& context, RunMode mode) {
  const int local_mode = static_cast<int>(mode);
  int minimum_mode = local_mode;
  int maximum_mode = local_mode;
  hundun::runtime::detail::check_mpi(MPI_Allreduce(&local_mode, &minimum_mode,
                                                   1, MPI_INT, MPI_MIN,
                                                   context.comm()),
                                     "MPI_Allreduce run-mode minimum");
  hundun::runtime::detail::check_mpi(MPI_Allreduce(&local_mode, &maximum_mode,
                                                   1, MPI_INT, MPI_MAX,
                                                   context.comm()),
                                     "MPI_Allreduce run-mode maximum");
  if (minimum_mode != maximum_mode) {
    throw Error("MPI run mode differs across communicator ranks");
  }
}

using DispatchCase =
    std::variant<hundun::config::ResolvedCaseV3,
                 hundun::config::ResolvedReactingCaseV4>;

DispatchCase load_and_broadcast_dispatch_case(
    const MpiContext& context, const std::filesystem::path& path) {
  std::optional<hundun::config::ResolvedCaseV3> legacy;
  std::optional<hundun::config::ResolvedReactingCaseV4> reacting;
  bool local_ok = true;
  std::string local_message;
  int kind = 0;
  if (context.rank() == 0) {
    try {
      reacting = hundun::config::load_resolved_reacting_case_v4(path);
      kind = 4;
    } catch (...) {
      try {
        legacy = hundun::config::load_resolved_case_v3(path);
        kind = 3;
      } catch (const std::exception& error) {
        local_ok = false;
        local_message = exception_message_or_fallback(
            error, "case configuration failed");
      } catch (...) {
        local_ok = false;
        local_message = "case configuration failed";
      }
    }
  }
  require_collective_success(context, local_ok, local_message);
  hundun::runtime::detail::check_mpi(
      MPI_Bcast(&kind, 1, MPI_INT, 0, context.comm()),
      "MPI_Bcast dispatch schema kind");
  if (kind == 4) {
    return hundun::config::broadcast_resolved_reacting_case_v4(
        context.comm(), 0,
        context.rank() == 0 ? &reacting.value() : nullptr);
  }
  return hundun::config::broadcast_resolved_case_v3(
      context.comm(), 0, context.rank() == 0 ? &legacy.value() : nullptr);
}

template <class Operation>
void write_root_output(const MpiContext& context,
                       std::string_view failure_message,
                       Operation&& operation) {
  bool local_ok = true;
  std::string local_message;
  if (context.rank() == 0) {
    try {
      std::forward<Operation>(operation)();
      std::cout.flush();
      if (!std::cout) {
        throw Error(std::string(failure_message));
      }
    } catch (const std::exception& error) {
      local_ok = false;
      local_message = exception_message_or_fallback(error, failure_message);
    } catch (...) {
      local_ok = false;
      local_message = failure_message;
    }
  }
  require_collective_success(context, local_ok, local_message);
}

int run_case(const hundun::application::CliOptions& options,
             MpiContext& context) {
  const RunMode mode = run_mode(options);
  require_run_mode_agreement(context, mode);
  return hundun::application::detail::dispatch_in_root_config_rank_order(
      [&] {
        return hundun::application::establish_authoritative_case_root(
            context, options.case_path);
      },
      [&] {
        return load_and_broadcast_dispatch_case(context, options.case_path);
      },
      [&](const DispatchCase& dispatch_case,
          const std::filesystem::path& authoritative_case_root) -> int {
        if (const auto* reacting =
                std::get_if<hundun::config::ResolvedReactingCaseV4>(
                    &dispatch_case)) {
          if (mode == RunMode::validate) {
            static_cast<void>(
                hundun::application::plan_reacting_flow_case(*reacting));
            write_root_output(context, "unable to write validation output",
                              [] { std::cout << "VALID\n"; });
            context.barrier();
            return EXIT_SUCCESS;
          }
          if (mode == RunMode::print_resolved) {
            write_root_output(
                context, "unable to write resolved case configuration",
                [&] {
                  std::cout
                      << hundun::config::to_resolved_reacting_json_v4(*reacting)
                      << '\n';
                });
            context.barrier();
            return EXIT_SUCCESS;
          }
          return hundun::application::run_reacting_flow_case(
              *reacting, context.comm());
        }
        const auto& resolved =
            std::get<hundun::config::ResolvedCaseV3>(dispatch_case);
        if (const auto* passive_scalar =
                std::get_if<hundun::config::CaseConfig>(&resolved)) {
          return hundun::application::run_passive_scalar_case(
              options, context, *passive_scalar, authoritative_case_root);
        }

        if (mode == RunMode::validate) {
          write_root_output(context, "unable to write validation output",
                            [] { std::cout << "VALID\n"; });
          context.barrier();
          return EXIT_SUCCESS;
        }
        if (mode == RunMode::print_resolved) {
          write_root_output(
              context, "unable to write resolved case configuration",
              [&resolved] {
                std::cout << hundun::config::to_resolved_json_v3(resolved)
                          << '\n';
              });
          context.barrier();
          return EXIT_SUCCESS;
        }
        if (const auto* flow =
                std::get_if<hundun::config::FlowCaseConfig>(&resolved)) {
          return hundun::application::run_flow_case(options, context, *flow,
                                                    authoritative_case_root);
        }
        return hundun::application::run_immersed_flow_case(
            options, context,
            std::get<hundun::config::ImmersedFlowCaseConfig>(resolved),
            authoritative_case_root);
      });
}

int run_with_mpi(int argc, char** argv,
                 const hundun::application::CliOptions& options) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  MpiContext context = MpiContext::duplicate(MPI_COMM_WORLD);
  try {
    return run_case(options, context);
  } catch (const std::exception& error) {
    if (context.rank() == 0) {
      std::cerr << exception_message_or_fallback(error, "HUNDUN-FLOW failed")
                << '\n';
    }
    return EXIT_FAILURE;
  } catch (...) {
    if (context.rank() == 0) {
      std::cerr << "HUNDUN-FLOW failed\n";
    }
    return EXIT_FAILURE;
  }
}

}  // namespace

int main(int argc, char** argv) {
  std::cout.imbue(std::locale::classic());
  std::cerr.imbue(std::locale::classic());
  try {
    const auto options = hundun::application::parse_cli(argc, argv);
    if (options.show_version) {
      std::cout << hundun::application::kHundunVersion << '\n';
      std::cout.flush();
      if (!std::cout) {
        throw Error("unable to write version output");
      }
      return EXIT_SUCCESS;
    }
    return run_with_mpi(argc, argv, options);
  } catch (const std::exception& error) {
    std::cerr << exception_message_or_fallback(error, "HUNDUN-FLOW failed")
              << '\n';
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "HUNDUN-FLOW failed\n";
    return EXIT_FAILURE;
  }
}
