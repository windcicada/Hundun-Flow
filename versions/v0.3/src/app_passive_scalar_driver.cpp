// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "app_case_config_broadcast_detail.hpp"
#include "app_cli_options_detail.hpp"
#include "app_passive_scalar_driver_detail.hpp"
#include "app_version.hpp"

#include "hundun/cfg_case_config_loader.hpp"
#include "hundun/mesh_uniform_structured.hpp"
#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_restart_binary.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "hundun/rt_vtk_legacy.hpp"
#include "hundun/flow_passive_scalar.hpp"
#include "rt_mpi_error_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

namespace {

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

enum class RunMode : int { normal = 0, validate = 1, print_resolved = 2 };

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

RunMode run_mode(const hundun::application::CliOptions &options) noexcept {
  if (options.validate_only) {
    return RunMode::validate;
  }
  if (options.print_resolved) {
    return RunMode::print_resolved;
  }
  return RunMode::normal;
}

void broadcast_bytes(const MpiContext &context, char *bytes,
                     std::uint64_t byte_count) {
  constexpr std::uint64_t chunk_limit =
      static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  std::uint64_t offset = 0;
  while (offset < byte_count) {
    const std::uint64_t remaining = byte_count - offset;
    const int count = static_cast<int>(std::min(remaining, chunk_limit));
    hundun::runtime::detail::check_mpi(
        MPI_Bcast(bytes + static_cast<std::size_t>(offset), count, MPI_BYTE, 0,
                  context.comm()),
        "MPI_Bcast authoritative case directory bytes");
    offset += static_cast<std::uint64_t>(count);
  }
}

hundun::runtime::HaloExchange create_halo_exchange(
    const MpiContext &context,
    const hundun::runtime::StructuredDecomposition &decomposition,
    int ghost_width) {
  constexpr std::string_view fallback = "unable to create halo exchange plan";
  std::optional<hundun::runtime::ExchangePlan> local_plan;
  bool local_ok = true;
  std::string local_message;
  try {
    local_plan.emplace(hundun::runtime::ExchangePlan::create(
        decomposition, decomposition.local_extent(), ghost_width));
  } catch (const std::exception &error) {
    local_ok = false;
    local_message = exception_message_or_fallback(error, fallback);
  } catch (...) {
    local_ok = false;
    local_message = fallback;
  }
  require_collective_success(context, local_ok, local_message);
  return hundun::runtime::HaloExchange::create(decomposition,
                                               std::move(local_plan).value());
}

template <class WriteOperation>
void write_root_output(const MpiContext &context,
                       std::string_view write_failure_message,
                       std::string_view exception_fallback,
                       WriteOperation &&write_operation) {
  bool local_ok = true;
  std::string local_message;
  if (context.rank() == 0) {
    try {
      std::forward<WriteOperation>(write_operation)();
      std::cout.flush();
      if (!std::cout) {
        local_ok = false;
        local_message = write_failure_message;
      }
    } catch (const std::exception &error) {
      local_ok = false;
      local_message = exception_message_or_fallback(error, exception_fallback);
    } catch (...) {
      local_ok = false;
      local_message = exception_fallback;
    }
  }
  require_collective_success(context, local_ok, local_message);
}

std::filesystem::path
resolve_case_root(const std::filesystem::path &case_path) {
  try {
    return std::filesystem::absolute(case_path)
        .lexically_normal()
        .parent_path();
  } catch (const std::filesystem::filesystem_error &error) {
    throw Error("unable to resolve case directory: " +
                std::string(error.what()));
  }
}

std::filesystem::path
broadcast_case_root(const MpiContext &context,
                    const std::filesystem::path &local_case_path) {
  constexpr std::string_view resolve_fallback =
      "unable to resolve case directory";
  std::string root_text;
  bool local_ok = true;
  std::string local_message;
  if (context.rank() == 0) {
    try {
      root_text = resolve_case_root(local_case_path).generic_string();
      if (root_text.size() > std::numeric_limits<std::uint64_t>::max()) {
        local_ok = false;
        local_message = "case directory path exceeds the uint64 wire domain";
      }
    } catch (const std::exception &error) {
      local_ok = false;
      local_message = exception_message_or_fallback(error, resolve_fallback);
    } catch (...) {
      local_ok = false;
      local_message = resolve_fallback;
    }
  }
  require_collective_success(context, local_ok, local_message);

  std::uint64_t root_length = 0;
  if (context.rank() == 0) {
    root_length = static_cast<std::uint64_t>(root_text.size());
  }
  hundun::runtime::detail::check_mpi(
      MPI_Bcast(&root_length, 1, MPI_UINT64_T, 0, context.comm()),
      "MPI_Bcast authoritative case directory length");

  local_ok = root_length <= static_cast<std::uint64_t>(
                                std::numeric_limits<std::size_t>::max());
  if (local_ok) {
    try {
      root_text.resize(static_cast<std::size_t>(root_length));
    } catch (...) {
      local_ok = false;
    }
  }
  require_collective_success(context, local_ok,
                             "unable to allocate authoritative case directory");
  broadcast_bytes(context, root_text.data(), root_length);

  std::optional<std::filesystem::path> root;
  local_ok = true;
  try {
    root.emplace(root_text);
  } catch (...) {
    local_ok = false;
  }
  require_collective_success(
      context, local_ok, "unable to construct authoritative case directory");
  return std::move(root).value();
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

std::pair<FieldId, FieldId> declare_passive_scalar_fields(FieldRegistry &registry) {
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

void require_finite_owned_passive_scalar_state(const MpiContext &context,
                                               const FieldRegistry &registry,
                                               const FieldStorage &storage,
                                               FieldId scalar) {
  bool local_ok = true;
  std::string_view local_message;
  try {
    const auto &descriptor = registry.descriptor(scalar);
    if (descriptor.restart != RestartPolicy::persistent ||
        descriptor.scalar_type != ScalarType::float64 ||
        descriptor.components >
            static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
      throw Error("passive-scalar field descriptor is not inspectable");
    }
    const auto values = storage.view<double>(scalar);
    const auto extent = storage.interior_extent();
    const int components = static_cast<int>(descriptor.components);
    for (int k = 0; k < extent.z && local_ok; ++k) {
      for (int j = 0; j < extent.y && local_ok; ++j) {
        for (int i = 0; i < extent.x && local_ok; ++i) {
          for (int component = 0; component < components; ++component) {
            if (!std::isfinite(values(i, j, k, component))) {
              local_ok = false;
              local_message = "passive-scalar state is not finite";
              break;
            }
          }
        }
      }
    }
  } catch (...) {
    local_ok = false;
    local_message = "unable to inspect passive-scalar state";
  }
  require_collective_success(context, local_ok, local_message);
}

int run_passive_scalar_case_impl(const hundun::application::CliOptions &options,
                         MpiContext &context, const CaseConfig &config,
                         const std::filesystem::path &root) {
  const RunMode mode = run_mode(options);
  hundun::runtime::require_expected_ranks(context, config.expected_ranks);

  if (mode == RunMode::validate) {
    write_root_output(context, "unable to write validation output",
                      "unable to write validation output",
                      [] { std::cout << "VALID\n"; });
    context.barrier();
    return EXIT_SUCCESS;
  }
  if (mode == RunMode::print_resolved) {
    write_root_output(
        context, "unable to write resolved case configuration",
        "unable to serialize resolved case configuration", [&config] {
          const std::string resolved = hundun::config::to_resolved_json(config);
          std::cout << resolved << '\n';
        });
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
  const auto fields = declare_passive_scalar_fields(registry);
  const FieldId scalar = fields.first;
  const FieldId stage = fields.second;
  FieldStorage storage(registry, decomposition.local_extent());
  auto halo = create_halo_exchange(context, decomposition, 2);
  hundun::solver::PassiveScalarSolver solver(
      context, decomposition, mesh, halo, config.transport.velocity_m_per_s,
      config.transport.diffusivity_m2_per_s);

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

  write_root_output(context, "unable to write passive-scalar output",
                    "unable to write passive-scalar output", [&config, &context] {
                      std::cout << hundun::application::kHundunVersion << '\n'
                                << "CASE name=" << config.case_name
                                << " ranks=" << context.size()
                                << " cells=" << config.mesh.cells.x << 'x'
                                << config.mesh.cells.y << 'x'
                                << config.mesh.cells.z << '\n';
                    });

  for (std::int64_t step = current_step + 1;
       step <= static_cast<std::int64_t>(config.time.steps); ++step) {
    const double next_time_s = current_time_s + config.time.dt_s;
    require_collective_success(context, std::isfinite(next_time_s),
                               "physical time is not finite");
    solver.advance_ssprk2(storage, scalar, stage, config.time.dt_s);
    require_finite_owned_passive_scalar_state(context, registry, storage,
                                              scalar);
    current_step = step;
    current_time_s = next_time_s;

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
      write_root_output(context, "unable to write passive-scalar output",
                        "unable to write passive-scalar output",
                        [step, current_time_s, mass, relative_mass_error] {
                          std::cout
                              << std::setprecision(17) << "STEP " << step
                              << " time_s=" << current_time_s
                              << " mass=" << mass
                              << " relative_mass_error=" << relative_mass_error
                              << '\n';
                        });
    }

    if (step % static_cast<std::int64_t>(config.output.restart_interval) == 0) {
      hundun::runtime::write_restart_checkpoint(
          context, decomposition, registry, storage,
          restart_write_directory / step_directory_name(step), step,
          current_time_s);
    }
  }

  context.barrier();
  write_root_output(context, "unable to write passive-scalar output",
                    "unable to write passive-scalar output",
                    [current_step, current_time_s] {
                      std::cout << std::setprecision(17)
                                << "FINISHED step=" << current_step
                                << " time_s=" << current_time_s << '\n';
                    });
  return EXIT_SUCCESS;
}

} // namespace

namespace hundun::application {

std::filesystem::path establish_authoritative_case_root(
    const runtime::MpiContext &context,
    const std::filesystem::path &local_case_path) {
  return broadcast_case_root(context, local_case_path);
}

int run_passive_scalar_case(const CliOptions &options, runtime::MpiContext &context,
                    const config::CaseConfig &config,
                    const std::filesystem::path &authoritative_case_root) {
  return run_passive_scalar_case_impl(options, context, config,
                              authoritative_case_root);
}

} // namespace hundun::application
