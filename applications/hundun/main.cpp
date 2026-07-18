// SPDX-License-Identifier: Apache-2.0

#include "applications/hundun/case_config_broadcast.hpp"
#include "applications/hundun/cli_options.hpp"

#include "hundun/config/case_config_loader.hpp"
#include "hundun/mesh/uniform_structured_mesh.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_descriptor.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/restart_binary.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "hundun/runtime/vtk_legacy.hpp"
#include "hundun/solver/passive_scalar.hpp"

#include <mpi.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <locale>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

namespace {

constexpr std::string_view kVersion = "HUNDUN-FLOW 0.0.0-stage1";
constexpr double kPi = 3.141592653589793238462643383279502884;

using hundun::config::CaseConfig;
using hundun::runtime::Error;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::MpiContext;
using hundun::runtime::OutputPolicy;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;

std::string exception_message_or_fallback(const std::exception &error,
                                          std::string_view fallback) {
  const char *message = error.what();
  return message != nullptr && *message != '\0' ? std::string(message)
                                                : std::string(fallback);
}

void require_collective_success(const MpiContext &context, bool local_ok,
                                std::string_view local_message) {
  const auto status =
      hundun::runtime::collective_status(context, local_ok, local_message);
  if (!status.ok) {
    throw Error(status.message);
  }
}

std::filesystem::path case_root(const std::filesystem::path &case_path) {
  try {
    return std::filesystem::absolute(case_path)
        .lexically_normal()
        .parent_path();
  } catch (const std::filesystem::filesystem_error &error) {
    throw Error("unable to resolve case directory: " +
                std::string(error.what()));
  }
}

std::string step_directory_name(std::int64_t step) {
  if (step < 0) {
    throw Error("Restart step must be nonnegative");
  }
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << "step" << std::setfill('0') << std::setw(8) << step;
  if (!output.good()) {
    throw Error("unable to format Restart step directory");
  }
  return output.str();
}

std::pair<FieldId, FieldId> declare_stage1_fields(FieldRegistry &registry) {
  const FieldId scalar = registry.declare_field(
      FieldDescriptor{"scalar", "1", "passive_scalar_solver",
                      FunctionSpace::cell_average, ScalarType::float64, 1U, 2,
                      true, RestartPolicy::persistent, OutputPolicy::selected});
  const FieldId stage = registry.declare_field(
      FieldDescriptor{"stage", "1", "passive_scalar_solver",
                      FunctionSpace::cell_average, ScalarType::float64, 1U, 2,
                      false, RestartPolicy::transient, OutputPolicy::never});
  registry.freeze();
  return {scalar, stage};
}

void initialize_sine(const CaseConfig &config,
                     const hundun::mesh::UniformStructuredMesh &mesh,
                     FieldStorage &storage, FieldId scalar, FieldId stage) {
  auto scalar_view = storage.view<double>(scalar);
  auto stage_view = storage.view<double>(stage);
  const auto extent = mesh.local_extent();
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        const double x = mesh.cell_center({i, j, k}).x;
        const double phase =
            2.0 * kPi * (x - config.mesh.origin_m.x) / config.mesh.length_m.x;
        scalar_view(i, j, k, 0) = 1.0 + 0.2 * std::sin(phase);
        stage_view(i, j, k, 0) = 0.0;
      }
    }
  }
}

void write_vtk_collectively(const MpiContext &context,
                            const std::filesystem::path &output_directory,
                            std::int64_t step,
                            const hundun::mesh::UniformStructuredMesh &mesh,
                            const FieldRegistry &registry,
                            const FieldStorage &storage, FieldId scalar) {
  bool local_ok = true;
  std::string local_message;
  try {
    hundun::runtime::write_vtk_rank(output_directory, step, context.rank(),
                                    mesh, registry, storage, scalar);
  } catch (const std::exception &error) {
    local_ok = false;
    local_message = exception_message_or_fallback(error, "VTK output failed");
  } catch (...) {
    local_ok = false;
    local_message = "VTK output failed";
  }
  require_collective_success(context, local_ok, local_message);
}

CaseConfig load_and_broadcast_case(const MpiContext &context,
                                   const std::filesystem::path &path) {
  std::optional<CaseConfig> root_config;
  bool local_ok = true;
  std::string local_message;
  if (context.rank() == 0) {
    try {
      root_config = hundun::config::load_case_config(path);
    } catch (const std::exception &error) {
      local_ok = false;
      local_message =
          exception_message_or_fallback(error, "case configuration failed");
    } catch (...) {
      local_ok = false;
      local_message = "case configuration failed";
    }
  }
  require_collective_success(context, local_ok, local_message);
  return hundun::application::broadcast_case_config(
      context.comm(), 0, context.rank() == 0 ? &root_config.value() : nullptr);
}

int run_case(const hundun::application::CliOptions &options,
             MpiContext &context) {
  const CaseConfig config = load_and_broadcast_case(context, options.case_path);
  hundun::runtime::require_expected_ranks(context, config.expected_ranks);

  if (options.validate_only) {
    if (context.rank() == 0) {
      std::cout << "VALID\n";
    }
    context.barrier();
    return EXIT_SUCCESS;
  }
  if (options.print_resolved) {
    bool local_ok = true;
    std::string local_message;
    if (context.rank() == 0) {
      try {
        const std::string resolved = hundun::config::to_resolved_json(config);
        std::cout << resolved << '\n';
        std::cout.flush();
        if (!std::cout) {
          local_ok = false;
          local_message = "unable to write resolved case configuration";
        }
      } catch (const std::exception &error) {
        local_ok = false;
        local_message = exception_message_or_fallback(
            error, "unable to serialize resolved case configuration");
      } catch (...) {
        local_ok = false;
        local_message = "unable to serialize resolved case configuration";
      }
    }
    require_collective_success(context, local_ok, local_message);
    context.barrier();
    return EXIT_SUCCESS;
  }

  const auto decomposition = hundun::runtime::StructuredDecomposition::create(
      context, config.mesh.cells, config.mesh.periodic,
      hundun::runtime::DecompositionOptions{config.process_grid});
  const hundun::mesh::UniformStructuredMesh mesh(
      config.mesh.cells, config.mesh.origin_m, config.mesh.length_m,
      decomposition);
  FieldRegistry registry;
  const auto fields = declare_stage1_fields(registry);
  const FieldId scalar = fields.first;
  const FieldId stage = fields.second;
  FieldStorage storage(registry, decomposition.local_extent());
  auto halo = hundun::runtime::HaloExchange::create(
      decomposition, hundun::runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 2));
  hundun::solver::PassiveScalarSolver solver(
      context, decomposition, mesh, halo, config.transport.velocity_m_per_s,
      config.transport.diffusivity_m2_per_s);

  const std::filesystem::path root = case_root(options.case_path);
  const std::filesystem::path output_directory =
      (root / config.output.directory).lexically_normal();
  const std::filesystem::path restart_write_directory =
      (root / config.restart.write_directory).lexically_normal();

  std::int64_t current_step = 0;
  double current_time_s = 0.0;
  if (config.restart.read) {
    const std::filesystem::path read_directory =
        (root / config.restart.read_directory.value()).lexically_normal();
    const auto metadata = hundun::runtime::read_restart_checkpoint(
        context, decomposition, registry, storage, read_directory);
    current_step = metadata.step;
    current_time_s = metadata.time_s;
  } else {
    initialize_sine(config, mesh, storage, scalar, stage);
  }
  if (current_step > static_cast<std::int64_t>(config.time.steps)) {
    throw Error("Restart step exceeds the configured absolute target step");
  }

  const double initial_mass =
      hundun::solver::global_mass(context, mesh, storage, scalar);
  if (!std::isfinite(initial_mass) || initial_mass == 0.0) {
    throw Error("initial passive-scalar mass must be finite and nonzero");
  }

  if (context.rank() == 0) {
    std::cout << kVersion << '\n'
              << "CASE name=" << config.case_name << " ranks=" << context.size()
              << " cells=" << config.mesh.cells.x << 'x' << config.mesh.cells.y
              << 'x' << config.mesh.cells.z << '\n';
  }

  for (std::int64_t step = current_step + 1;
       step <= static_cast<std::int64_t>(config.time.steps); ++step) {
    solver.advance_ssprk2(storage, scalar, stage, config.time.dt_s);
    current_step = step;
    current_time_s += config.time.dt_s;

    if (step % static_cast<std::int64_t>(config.output.write_interval) == 0) {
      const double mass =
          hundun::solver::global_mass(context, mesh, storage, scalar);
      const double relative_mass_error =
          std::abs(mass - initial_mass) / std::abs(initial_mass);
      if (!std::isfinite(mass) || !std::isfinite(relative_mass_error)) {
        throw Error("passive-scalar mass diagnostic is not finite");
      }
      write_vtk_collectively(context, output_directory, step, mesh, registry,
                             storage, scalar);
      if (context.rank() == 0) {
        std::cout << std::setprecision(17) << "STEP " << step
                  << " time_s=" << current_time_s << " mass=" << mass
                  << " relative_mass_error=" << relative_mass_error << '\n';
      }
    }

    if (step % static_cast<std::int64_t>(config.output.restart_interval) == 0) {
      hundun::runtime::write_restart_checkpoint(
          context, decomposition, registry, storage,
          restart_write_directory / step_directory_name(step), step,
          current_time_s);
    }
  }

  context.barrier();
  if (context.rank() == 0) {
    std::cout << std::setprecision(17) << "FINISHED step=" << current_step
              << " time_s=" << current_time_s << '\n';
  }
  return EXIT_SUCCESS;
}

int run_with_mpi(int argc, char **argv,
                 const hundun::application::CliOptions &options) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  MpiContext context = MpiContext::duplicate(MPI_COMM_WORLD);
  try {
    return run_case(options, context);
  } catch (const std::exception &error) {
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

} // namespace

int main(int argc, char **argv) {
  std::cout.imbue(std::locale::classic());
  std::cerr.imbue(std::locale::classic());
  try {
    const auto options = hundun::application::parse_cli(argc, argv);
    if (options.show_version) {
      std::cout << kVersion << '\n';
      return EXIT_SUCCESS;
    }
    return run_with_mpi(argc, argv, options);
  } catch (const std::exception &error) {
    std::cerr << exception_message_or_fallback(error, "HUNDUN-FLOW failed")
              << '\n';
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "HUNDUN-FLOW failed\n";
    return EXIT_FAILURE;
  }
}
