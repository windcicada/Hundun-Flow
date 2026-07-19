// SPDX-License-Identifier: Apache-2.0

#include "applications/hundun/case_config_broadcast.hpp"
#include "applications/hundun/cli_options.hpp"
#include "applications/hundun/detail/dispatch_order.hpp"
#include "applications/hundun/stage1_driver.hpp"

#include "hundun/config/resolved_case_loader.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "runtime/src/mpi_error.hpp"

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
  hundun::runtime::detail::check_mpi(
      MPI_Allreduce(&local_mode, &minimum_mode, 1, MPI_INT, MPI_MIN,
                    context.comm()),
      "MPI_Allreduce run-mode minimum");
  hundun::runtime::detail::check_mpi(
      MPI_Allreduce(&local_mode, &maximum_mode, 1, MPI_INT, MPI_MAX,
                    context.comm()),
      "MPI_Allreduce run-mode maximum");
  if (minimum_mode != maximum_mode) {
    throw Error("MPI run mode differs across communicator ranks");
  }
}

hundun::config::ResolvedCase load_and_broadcast_resolved_case(
    const MpiContext& context, const std::filesystem::path& path) {
  std::optional<hundun::config::ResolvedCase> root_case;
  bool local_ok = true;
  std::string local_message;
  if (context.rank() == 0) {
    try {
      root_case = hundun::config::load_resolved_case(path);
    } catch (const std::exception& error) {
      local_ok = false;
      local_message =
          exception_message_or_fallback(error, "case configuration failed");
    } catch (...) {
      local_ok = false;
      local_message = "case configuration failed";
    }
  }
  require_collective_success(context, local_ok, local_message);
  return hundun::application::broadcast_resolved_case(
      context.comm(), 0, context.rank() == 0 ? &root_case.value() : nullptr);
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
        return load_and_broadcast_resolved_case(context, options.case_path);
      },
      [&](const hundun::config::ResolvedCase& resolved,
          const std::filesystem::path& authoritative_case_root) -> int {
        if (const auto* stage1 =
                std::get_if<hundun::config::CaseConfig>(&resolved)) {
          return hundun::application::run_stage1_case(
              options, context, *stage1, authoritative_case_root);
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
                std::cout << hundun::config::to_resolved_json(resolved)
                          << '\n';
              });
          context.barrier();
          return EXIT_SUCCESS;
        }
        throw Error(
            "Stage 2 variable-density flow driver is not implemented before Task 24");
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
