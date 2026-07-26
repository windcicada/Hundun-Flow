// SPDX-License-Identifier: Apache-2.0

#include "applications/hundun/stage2_driver.hpp"

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/diagnostics/checkpoint_v2_diagnostics.hpp"
#include "hundun/diagnostics/diagnostic_session.hpp"
#include "hundun/diagnostics/ideal_gas_closure_diagnostics.hpp"
#include "hundun/diagnostics/material_density_piso_diagnostics.hpp"
#include "hundun/diagnostics/mesh_diagnostic_v2.hpp"
#include "hundun/diagnostics/stage2_module_diagnostics.hpp"
#include "hundun/diagnostics/time_control_diagnostics.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/adaptive_time_control.hpp"
#include "hundun/flow/checkpoint_v2.hpp"
#include "hundun/flow/constant_density_piso.hpp"
#include "hundun/flow/ideal_gas_closure.hpp"
#include "hundun/flow/ideal_gas_piso.hpp"
#include "hundun/flow/material_density_piso.hpp"
#include "hundun/linear/bicgstab.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/structured_decomposition.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace hundun::application {
namespace {

using config::DensityModel;
using diagnostics::DiagnosticBatch;
using diagnostics::DiagnosticBatchSink;
using diagnostics::DiagnosticLevel;
using diagnostics::DiagnosticRequest;
using diagnostics::DiagnosticScope;
using flow::FlowFieldIds;
using flow::FlowLayerValues;
using runtime::Error;
using runtime::FieldDescriptor;
using runtime::FieldId;
using runtime::FieldRegistry;
using runtime::FunctionSpace;
using runtime::MpiContext;

constexpr std::string_view kStage2Banner = "HUNDUN-FLOW 0.0.0-stage2";
constexpr std::string_view kFieldOwner = "stage2-flow-driver";

std::string density_model_name(DensityModel model) {
  switch (model) {
  case DensityModel::constant:
    return "constant";
  case DensityModel::material:
    return "material";
  case DensityModel::ideal_gas:
    return "ideal_gas";
  }
  throw Error("invalid Stage 2 density model");
}

std::string_view checkpoint_phase(flow::CheckpointV2Phase phase) {
  switch (phase) {
  case flow::CheckpointV2Phase::none:
    return "none";
  case flow::CheckpointV2Phase::preflight:
    return "preflight";
  case flow::CheckpointV2Phase::transaction_entry:
    return "transaction-entry";
  case flow::CheckpointV2Phase::rank_payload:
    return "rank-payload";
  case flow::CheckpointV2Phase::rank_temporary_file:
    return "rank-temporary-file";
  case flow::CheckpointV2Phase::rank_publish:
    return "rank-publish";
  case flow::CheckpointV2Phase::manifest:
    return "manifest";
  case flow::CheckpointV2Phase::completed_marker:
    return "completed-marker";
  case flow::CheckpointV2Phase::marker_read:
    return "marker-read";
  case flow::CheckpointV2Phase::manifest_read:
    return "manifest-read";
  case flow::CheckpointV2Phase::rank_read:
    return "rank-read";
  case flow::CheckpointV2Phase::restore_prepare:
    return "restore-prepare";
  case flow::CheckpointV2Phase::restore_publish:
    return "restore-publish";
  }
  throw Error("invalid Checkpoint v2 phase");
}

FieldDescriptor cell_field(std::string name, std::string unit,
                           std::uint32_t components, bool conservative) {
  return {std::move(name),
          std::move(unit),
          std::string(kFieldOwner),
          FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          2,
          conservative,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

FieldDescriptor face_field(std::string name, std::string unit,
                           std::uint32_t components) {
  return {std::move(name),
          std::move(unit),
          std::string(kFieldOwner),
          FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

struct DeclaredFields final {
  FieldRegistry registry;
  FlowFieldIds ids;
  std::vector<FieldId> transported;
};

DeclaredFields declare_fields(const config::FlowCaseConfig& config) {
  DeclaredFields result;
  result.ids.density =
      result.registry.declare_field(cell_field("rho", "kg/m3", 1U, true));
  result.ids.velocity =
      result.registry.declare_field(cell_field("u", "m/s", 3U, false));
  result.ids.mechanical_pressure =
      result.registry.declare_field(cell_field("pi", "Pa", 1U, false));
  result.ids.face_velocity =
      result.registry.declare_field(face_field("uf", "m/s", 3U));
  result.ids.face_mass_flux =
      finite_volume::declare_face_mass_flux(result.registry);

  if (config.density_model == DensityModel::constant) {
    result.transported.push_back(
        result.registry.declare_field(cell_field("h", "J/kg", 1U, false)));
    for (const auto& scalar : config.scalars)
      result.transported.push_back(result.registry.declare_field(
          cell_field(scalar.name, "1", 1U, false)));
  } else {
    result.transported.push_back(result.registry.declare_field(
        cell_field("rho_h", "J/m3", 1U, true)));
    for (const auto& scalar : config.scalars)
      result.transported.push_back(result.registry.declare_field(
          cell_field("rho_" + scalar.name, "kg/m3", 1U, true)));
  }
  result.ids.transported_cell_fields = result.transported;
  result.registry.freeze();
  return result;
}

std::array<bool, 3>
periodicity(const config::FlowCaseConfig& config) noexcept {
  std::array<bool, 3> result{};
  for (std::size_t axis = 0; axis < 3U; ++axis)
    result[axis] =
        config.boundaries[2U * axis].type == config::BoundaryType::periodic;
  return result;
}

runtime::Int3 local_extent(
    const runtime::StructuredDecomposition& decomposition) noexcept {
  const auto box = decomposition.owned_box();
  return {box.end.x - box.begin.x, box.end.y - box.begin.y,
          box.end.z - box.begin.z};
}

const boundary::ResolvedInletState* inlet_state(
    const boundary::BoundaryRegistry& boundaries) {
  const auto patch = boundaries.velocity_inlet_patch_id();
  if (!patch)
    return nullptr;
  const auto& state = boundaries.patch(*patch).inlet_state();
  if (!state)
    throw Error("velocity inlet has no resolved inlet state");
  return &*state;
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

FlowLayerValues initial_values(const config::FlowCaseConfig& config,
                               const mesh::MeshTopology& topology,
                               const mesh::MeshGeometry& geometry,
                               const boundary::BoundaryRegistry& boundaries) {
  const std::size_t cells = topology.owned_cell_count();
  const std::size_t faces = topology.local_face_count();
  const auto* inlet = inlet_state(boundaries);

  double density = config.physics.rho_ref_kg_per_m3;
  runtime::Real3 velocity{};
  double enthalpy = 0.0;
  std::vector<double> scalar_values(config.scalars.size(), 0.0);
  if (inlet != nullptr) {
    density = inlet->density_kg_per_m3;
    velocity = inlet->velocity_m_per_s;
    enthalpy = inlet->enthalpy_J_per_kg;
    scalar_values = inlet->scalar_values;
  } else if (config.density_model == DensityModel::ideal_gas) {
    const double pressure =
        config.physics.thermodynamic_pressure_pa.value();
    const double gas_constant =
        config.physics.gas_constant_J_per_kg_K.value();
    const double cp = config.physics.cp_J_per_kg_K.value();
    const double temperature =
        pressure / (gas_constant * config.physics.rho_ref_kg_per_m3);
    enthalpy = cp * temperature;
    if (!std::isfinite(enthalpy) || enthalpy <= 0.0)
      throw Error("invalid derived ideal-gas initial enthalpy");
  }
  if (!std::isfinite(density) || density <= 0.0)
    throw Error("invalid Stage 2 initial density");

  FlowLayerValues result;
  result.density.assign(cells, density);
  result.velocity.resize(cells * 3U);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    result.velocity[cell * 3U] = velocity.x;
    result.velocity[cell * 3U + 1U] = velocity.y;
    result.velocity[cell * 3U + 2U] = velocity.z;
  }
  result.mechanical_pressure.assign(cells, 0.0);
  result.face_velocity.resize(faces * 3U);
  result.face_mass_flux.resize(faces);
  for (std::size_t face = 0; face < faces; ++face) {
    result.face_velocity[face * 3U] = velocity.x;
    result.face_velocity[face * 3U + 1U] = velocity.y;
    result.face_velocity[face * 3U + 2U] = velocity.z;
    result.face_mass_flux[face] =
        density * dot(velocity, geometry.face_area_vector_m2(
                                    face, mesh::FaceSide::owner));
  }

  result.transported_cell_fields.reserve(config.scalars.size() + 1U);
  const bool conservative = config.density_model != DensityModel::constant;
  result.transported_cell_fields.push_back(
      std::vector<double>(cells, conservative ? density * enthalpy : enthalpy));
  for (std::size_t scalar = 0; scalar < config.scalars.size(); ++scalar)
    result.transported_cell_fields.push_back(std::vector<double>(
        cells, conservative ? density * scalar_values[scalar]
                            : scalar_values[scalar]));
  return result;
}

std::vector<flow::ConstantDensityTransportSpec>
constant_transport_specs(const DeclaredFields& fields,
                         const config::FlowCaseConfig& config) {
  std::vector<flow::ConstantDensityTransportSpec> result;
  result.push_back({fields.transported.front(),
                    finite_volume::FiniteVolumeQuantity::enthalpy(), 0.0});
  for (std::size_t scalar = 0; scalar < config.scalars.size(); ++scalar) {
    result.push_back(
        {fields.transported[scalar + 1U],
         finite_volume::FiniteVolumeQuantity::scalar(
             static_cast<std::uint32_t>(scalar)),
         config.physics.rho_ref_kg_per_m3 *
             config.scalars[scalar].diffusivity_m2_per_s});
  }
  return result;
}

flow::MaterialDensityTransportSpec
material_transport_spec(const DeclaredFields& fields,
                        const config::FlowCaseConfig& config) {
  flow::MaterialDensityTransportSpec result;
  result.enthalpy_density = fields.transported.front();
  result.enthalpy_diffusivity_kg_per_m_s = 0.0;
  for (std::size_t scalar = 0; scalar < config.scalars.size(); ++scalar) {
    result.scalar_densities.push_back(fields.transported[scalar + 1U]);
    result.scalar_diffusivities_kg_per_m_s.push_back(
        config.physics.rho_ref_kg_per_m3 *
        config.scalars[scalar].diffusivity_m2_per_s);
  }
  return result;
}

std::filesystem::path beneath(const std::filesystem::path& root,
                              const std::filesystem::path& relative) {
  return (root / relative).lexically_normal();
}

std::filesystem::path checkpoint_step_path(
    const std::filesystem::path& root, std::uint64_t step) {
  std::ostringstream name;
  name.imbue(std::locale::classic());
  name << "step" << std::setw(20) << std::setfill('0') << step;
  return root / name.str();
}

void ensure_directory(const MpiContext& mpi,
                      const std::filesystem::path& directory) {
  bool local_ok = true;
  std::string message;
  if (mpi.rank() == 0) {
    try {
      std::filesystem::create_directories(directory);
    } catch (const std::exception& error) {
      local_ok = false;
      message = error.what();
    }
  }
  const auto status = runtime::collective_status(mpi, local_ok, message);
  if (!status.ok)
    throw Error(status.message);
}

void require_collective_success(const MpiContext& mpi, bool local_ok,
                                std::string_view message) {
  const auto status = runtime::collective_status(mpi, local_ok, message);
  if (!status.ok)
    throw Error(status.message);
}

template <class Source>
void append_summary(const Source& source, std::uint64_t step, double time_s,
                    int rank, std::string_view phase, DiagnosticBatch& batch) {
  const auto descriptor = diagnostics::describe_diagnostics(source);
  DiagnosticRequest request{DiagnosticLevel::summary,
                            DiagnosticScope::local,
                            {rank, step, time_s, phase},
                            {},
                            0U};
  DiagnosticBatchSink sink(descriptor, request, batch);
  diagnostics::collect_diagnostics(source, request, sink);
}

diagnostics::FieldLayoutDiagnosticSource field_layout_source(
    const DeclaredFields& fields,
    const runtime::StructuredDecomposition& decomposition,
    const mesh::MeshTopology& topology) {
  diagnostics::FieldLayoutDiagnosticSource source{
      &fields.registry,
      {local_extent(decomposition), topology.local_face_count()},
      {{"density", fields.ids.density},
       {"velocity", fields.ids.velocity},
       {"mechanical-pressure", fields.ids.mechanical_pressure},
       {"face-velocity", fields.ids.face_velocity},
       {"face-mass-flux", fields.ids.face_mass_flux}}};
  for (std::size_t field = 0;
       field < fields.ids.transported_cell_fields.size(); ++field)
    source.roles.push_back(
        {"transported." + std::to_string(field),
         fields.ids.transported_cell_fields[field]});
  return source;
}

void append_static_prefix_diagnostics(
    const MpiContext& mpi,
    const runtime::StructuredDecomposition& decomposition,
    const diagnostics::FieldLayoutDiagnosticSource& field_layout,
    const runtime::ExchangePlan& exchange_plan,
    const mesh::MeshTopology& topology, const mesh::MeshGeometry& geometry,
    const execution::ExecutionContext& execution,
    std::uint64_t step, double time_s,
    DiagnosticBatch& batch) {
  append_summary(mpi, step, time_s, mpi.rank(), "accepted-step", batch);
  append_summary(decomposition, step, time_s, mpi.rank(), "accepted-step",
                 batch);
  append_summary(field_layout, step, time_s, mpi.rank(), "accepted-step",
                 batch);
  append_summary(exchange_plan, step, time_s, mpi.rank(), "accepted-step",
                 batch);
  append_summary(topology, step, time_s, mpi.rank(), "accepted-step", batch);
  append_summary(geometry, step, time_s, mpi.rank(), "accepted-step", batch);
  append_summary(execution, step, time_s, mpi.rank(), "accepted-step", batch);
}

void append_static_suffix_diagnostics(
    const MpiContext& mpi, const boundary::BoundaryRegistry& boundaries,
    const diagnostics::SharedFluxDiagnosticSource& shared_flux,
    std::uint64_t step, double time_s, DiagnosticBatch& batch) {
  append_summary(shared_flux, step, time_s, mpi.rank(), "accepted-step",
                 batch);
  append_summary(boundaries, step, time_s, mpi.rank(), "accepted-step", batch);
}

const flow::StepAttemptReport& base_flow_report(
    const flow::DensityStepAttemptReport& report) {
  if (const auto* constant = std::get_if<flow::StepAttemptReport>(&report))
    return *constant;
  if (const auto* material =
          std::get_if<flow::MaterialDensityStepAttemptReport>(&report))
    return material->flow();
  return std::get<flow::IdealGasStepAttemptReport>(report).flow().flow();
}

void append_linear_solve_diagnostics(const flow::StepAttemptReport& report,
                                     std::uint64_t step, double time_s,
                                     int rank, DiagnosticBatch& batch) {
  constexpr std::array<std::string_view, 3> momentum_ids{
      "momentum-x", "momentum-y", "momentum-z"};
  for (std::size_t component = 0; component < momentum_ids.size();
       ++component)
    append_summary(
        diagnostics::LinearSolveDiagnosticSource{
            momentum_ids[component], &report.momentum.components[component]},
        step, time_s, rank, "accepted-step", batch);
  constexpr std::array<std::string_view, 2> pressure_ids{
      "pressure-corrector-1", "pressure-corrector-2"};
  for (std::size_t corrector = 0; corrector < pressure_ids.size();
       ++corrector)
    append_summary(
        diagnostics::LinearSolveDiagnosticSource{
            pressure_ids[corrector], &report.pressure[corrector]},
        step, time_s, rank, "accepted-step", batch);
}

template <class AdvanceReport>
void append_time_diagnostics(const flow::Bdf2RetryController& controller,
                             const flow::FlowState& state,
                             const AdvanceReport& advance, int rank,
                             DiagnosticBatch& batch) {
  auto source = controller.diagnostic_source(state, advance);
  append_summary(source, state.metadata().step, state.metadata().time_s, rank,
                 "time-control.advance-result", batch);
}

template <class Flow>
std::optional<flow::IdealGasClosureState>
closure_state_for_checkpoint(const Flow&) {
  return std::nullopt;
}

template <>
std::optional<flow::IdealGasClosureState>
closure_state_for_checkpoint(const flow::FixedStepIdealGasFlow& value) {
  return value.closure_state();
}

void print_start(const config::FlowCaseConfig& config,
                 const MpiContext& mpi) {
  bool local_ok = true;
  if (mpi.rank() == 0) {
    std::cout << kStage2Banner << '\n'
              << "CASE name=" << config.case_name << " ranks=" << mpi.size()
              << " cells=" << config.mesh.cells.x << 'x'
              << config.mesh.cells.y << 'x' << config.mesh.cells.z
              << " density_model="
              << density_model_name(config.density_model) << '\n';
    std::cout.flush();
    local_ok = static_cast<bool>(std::cout);
  }
  require_collective_success(mpi, local_ok,
                             "unable to write Stage 2 output");
}

void print_step(const flow::FlowState& state, std::size_t attempts,
                const MpiContext& mpi) {
  bool local_ok = true;
  if (mpi.rank() == 0) {
    const auto metadata = state.metadata();
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "STEP " << metadata.step << " time_s=" << metadata.time_s
              << " dt_s=" << metadata.dt_s << " attempts=" << attempts
              << '\n';
    std::cout.flush();
    local_ok = static_cast<bool>(std::cout);
  }
  require_collective_success(mpi, local_ok,
                             "unable to write Stage 2 output");
}

void print_finished(const flow::FlowState& state, const MpiContext& mpi) {
  bool local_ok = true;
  if (mpi.rank() == 0) {
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "FINISHED step=" << state.metadata().step
              << " time_s=" << state.metadata().time_s << '\n';
    std::cout.flush();
    local_ok = static_cast<bool>(std::cout);
  }
  require_collective_success(mpi, local_ok,
                             "unable to write Stage 2 output");
}

bool has_scheduled_diagnostics(std::uint64_t current_step,
                               std::uint64_t target_step,
                               std::uint64_t interval,
                               bool include_current) noexcept {
  if (include_current && current_step != 0U &&
      current_step % interval == 0U)
    return true;
  if (current_step >= target_step)
    return false;
  const auto first =
      ((current_step / interval) + 1U) * interval;
  return first <= target_step;
}

template <class Flow>
void run_loop(
    const config::FlowCaseConfig& config, const MpiContext& mpi,
    const runtime::StructuredDecomposition& decomposition,
    const mesh::MeshTopology& topology, const mesh::MeshGeometry& geometry,
    const boundary::BoundaryRegistry& boundaries,
    const DeclaredFields& fields, const runtime::ExchangePlan& exchange_plan,
    execution::ExecutionContext& execution, flow::FlowState& state, Flow& flow,
    flow::Bdf2RetryController& controller,
    diagnostics::DiagnosticSession* diagnostic_session,
    const std::filesystem::path& checkpoint_root) {
  const linear::SolveControl solve_control{};
  while (state.metadata().step <
         static_cast<std::uint64_t>(config.time.steps)) {
    flow::TimeAdvanceReport advance = [&] {
      if constexpr (std::is_same_v<Flow,
                                   flow::FixedStepConstantDensityFlow>) {
        return controller.advance(
            state, flow, config.physics.rho_ref_kg_per_m3,
            config.physics.dynamic_viscosity_pa_s, solve_control,
            solve_control);
      } else {
        return controller.advance(state, flow,
                                  config.physics.dynamic_viscosity_pa_s,
                                  solve_control, solve_control);
      }
    }();
    if (advance.disposition() != flow::TimeAdvanceDisposition::committed)
      throw Error("Stage 2 time advance failed");

    std::optional<flow::CheckpointV2Report> checkpoint_report;
    if (state.metadata().step %
            static_cast<std::uint64_t>(config.restart.write_interval) ==
        0U) {
      checkpoint_report = flow::write_checkpoint_v2(
          mpi, decomposition, topology, geometry, boundaries, config, state,
          controller.state(), closure_state_for_checkpoint(flow),
          checkpoint_step_path(checkpoint_root, state.metadata().step));
      if (checkpoint_report->disposition() !=
          flow::CheckpointV2Disposition::completed)
        throw Error(
            "Stage 2 checkpoint write failed: reason=" +
            std::to_string(static_cast<int>(checkpoint_report->reason())) +
            " phase=" +
            std::to_string(static_cast<int>(checkpoint_report->phase())) +
            " rank=" +
            std::to_string(checkpoint_report->lowest_failing_rank()));
    }

    if (diagnostic_session != nullptr &&
        diagnostic_session->due(state.metadata().step)) {
      DiagnosticBatch batch;
      bool batch_ok = true;
      std::string batch_message;
      try {
        const auto& final_attempt = advance.final_attempt();
        const auto& report = base_flow_report(final_attempt);
        const auto reason = advance.reason();
        diagnostics::FlowDriverDiagnosticSource driver_source{
            config.density_model,
            state.metadata().step,
            state.metadata().time_s,
            advance.attempt_count(),
            advance.disposition(),
            reason,
            advance.lowest_failing_rank()};
        auto field_layout =
            field_layout_source(fields, decomposition, topology);
        append_static_prefix_diagnostics(
            mpi, decomposition, field_layout, exchange_plan, topology,
            geometry, execution, state.metadata().step,
            state.metadata().time_s, batch);
        append_linear_solve_diagnostics(
            report, state.metadata().step, state.metadata().time_s,
            mpi.rank(), batch);
        append_static_suffix_diagnostics(
            mpi, boundaries,
            {fields.ids.face_mass_flux, topology.local_face_count(), true},
            state.metadata().step, state.metadata().time_s, batch);
        if constexpr (std::is_same_v<
                          Flow, flow::FixedStepConstantDensityFlow>) {
          append_summary(
              diagnostics::ConstantDensityPisoDiagnosticSource{&report},
              state.metadata().step, state.metadata().time_s, mpi.rank(),
              "accepted-step", batch);
          append_summary(driver_source, state.metadata().step,
                         state.metadata().time_s, mpi.rank(), "accepted-step",
                         batch);
        } else if constexpr (std::is_same_v<
                                 Flow, flow::FixedStepMaterialDensityFlow>) {
          append_summary(driver_source, state.metadata().step,
                         state.metadata().time_s, mpi.rank(), "accepted-step",
                         batch);
          const auto& material_report =
              std::get<flow::MaterialDensityStepAttemptReport>(final_attempt);
          auto source = flow.diagnostic_source(state, material_report);
          append_summary(source, state.metadata().step,
                         state.metadata().time_s, mpi.rank(),
                         "material-density.attempt-result", batch);
        } else {
          const auto& ideal_report =
              std::get<flow::IdealGasStepAttemptReport>(final_attempt);
          append_summary(driver_source, state.metadata().step,
                         state.metadata().time_s, mpi.rank(), "accepted-step",
                         batch);
          {
            auto source =
                flow.flow_diagnostic_source(state, ideal_report);
            append_summary(source, state.metadata().step,
                           state.metadata().time_s, mpi.rank(),
                           "material-density.attempt-result", batch);
          }
          {
            auto source =
                flow.closure_diagnostic_source(state, ideal_report);
            append_summary(source, state.metadata().step,
                           state.metadata().time_s, mpi.rank(),
                           "ideal-gas-closure.attempt-result", batch);
          }
        }
        append_time_diagnostics(controller, state, advance, mpi.rank(), batch);
        if (checkpoint_report) {
          auto source =
              flow::checkpoint_v2_diagnostic_source(*checkpoint_report);
          append_summary(
              source, state.metadata().step, state.metadata().time_s,
              mpi.rank(), checkpoint_phase(checkpoint_report->phase()),
              batch);
        }
      } catch (const std::exception& error) {
        batch_ok = false;
        batch_message = error.what();
      }
      require_collective_success(
          mpi, batch_ok,
          batch_message.empty() ? "unable to construct diagnostic batch"
                                : batch_message);
      diagnostic_session->publish(mpi, state.metadata().step, batch);
    }
    print_step(state, advance.attempt_count(), mpi);
  }
}

}  // namespace

int run_stage2_case(const CliOptions&, MpiContext& mpi,
                    const config::FlowCaseConfig& config,
                    const std::filesystem::path& authoritative_case_root) {
  require_collective_success(
      mpi,
      !config.resources.expected_ranks ||
          *config.resources.expected_ranks == mpi.size(),
      "configured MPI rank count differs from runtime");
  require_collective_success(
      mpi, !config.performance.enabled,
      "Stage 2 performance artifacts are unavailable before Task 25");

  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, config.mesh.cells, periodicity(config),
      runtime::DecompositionOptions{config.resources.process_grid});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry =
      config.mesh.mapping == config::MeshMapping::uniform_box
          ? mesh::MeshGeometry(
                topology,
                mesh::UniformBoxMapping(config.mesh.origin_m,
                                        config.mesh.length_m))
          : mesh::MeshGeometry(
                topology,
                mesh::AnalyticWarpedBoxMapping(
                    config.mesh.origin_m, config.mesh.length_m,
                    config.mesh.warp_amplitude.value()));
  auto boundaries = boundary::BoundaryRegistry::create(config, topology);
  auto fields = declare_fields(config);
  auto state = flow::FlowState::create(
      fields.registry,
      {local_extent(decomposition), topology.local_face_count()}, fields.ids,
      {0U, 0.0, config.time.initial_dt_s, 0.0,
       flow::MomentumTimeOrder::backward_euler});
  const auto seed = initial_values(config, topology, geometry, boundaries);
  state.seed_accepted_layers(seed, seed);

  std::optional<flow::CheckpointV2ReadResult> restored;
  if (config.restart.read) {
    restored = flow::read_checkpoint_v2(
        mpi, decomposition, topology, geometry, boundaries, config, state,
        beneath(authoritative_case_root, *config.restart.read_directory));
    if (!restored->restored())
      throw Error("Stage 2 checkpoint read failed");
    if (state.metadata().step > static_cast<std::uint64_t>(config.time.steps))
      throw Error("restored step exceeds configured target step");
  }

  auto exchange_plan = runtime::ExchangePlan::create(
      decomposition, local_extent(decomposition), 2);
  auto halo =
      runtime::HaloExchange::create(decomposition, exchange_plan);
  execution::CpuReferenceContext execution;
  linear::BiCGStabSolver momentum_solver(execution, mpi);
  linear::ConjugateGradientSolver pressure_cg(execution, mpi);
  linear::BiCGStabSolver pressure_bicgstab(execution, mpi);
  const linear::LinearSolver& pressure_solver =
      config.mesh.mapping == config::MeshMapping::uniform_box
          ? static_cast<const linear::LinearSolver&>(pressure_cg)
          : static_cast<const linear::LinearSolver&>(pressure_bicgstab);
  linear::JacobiPreconditioner momentum_x(execution), momentum_y(execution),
      momentum_z(execution), pressure_preconditioner(execution);

  auto controller = restored
                        ? flow::Bdf2RetryController::restore(
                              config.time, config.density_model, topology,
                              geometry, mpi, state,
                              restored->time_control_state())
                        : flow::Bdf2RetryController::create(
                              config.time, config.density_model, topology,
                              geometry, mpi, state);
  const auto diagnostics_directory =
      beneath(authoritative_case_root, config.diagnostics.directory);
  std::optional<diagnostics::DiagnosticSession> diagnostic_session;
  const auto target_step =
      static_cast<std::uint64_t>(config.time.steps);
  if (config.diagnostics.write_mesh ||
      has_scheduled_diagnostics(
          state.metadata().step, target_step,
          static_cast<std::uint64_t>(
              config.diagnostics.write_interval),
          restored.has_value()))
    diagnostic_session.emplace(
        diagnostics_directory, config.diagnostics.write_interval,
        mpi.rank());
  const auto checkpoint_root =
      beneath(authoritative_case_root, config.restart.write_directory);
  ensure_directory(mpi, checkpoint_root);

  if (config.diagnostics.write_mesh) {
    std::optional<diagnostics::MeshDiagnosticV2> mesh_record;
    bool mesh_ok = true;
    std::string mesh_message;
    try {
      mesh_record = diagnostics::make_mesh_diagnostic_v2(
          mpi.rank(), mpi.size(), topology, geometry);
    } catch (const std::exception& error) {
      mesh_ok = false;
      mesh_message = error.what();
    }
    require_collective_success(
        mpi, mesh_ok,
        mesh_message.empty() ? "unable to construct mesh diagnostic"
                             : mesh_message);
    diagnostics::write_mesh_diagnostic_v2(
        mpi, diagnostics_directory, *mesh_record);
  }

  print_start(config, mpi);
  if (restored && diagnostic_session &&
      diagnostic_session->due(state.metadata().step)) {
    DiagnosticBatch batch;
    bool batch_ok = true;
    std::string batch_message;
    try {
      const diagnostics::FlowDriverDiagnosticSource driver_source{
          config.density_model,
          state.metadata().step,
          state.metadata().time_s,
          0U,
          flow::TimeAdvanceDisposition::committed,
          flow::StepFailureReason::none,
          -1};
      auto field_layout =
          field_layout_source(fields, decomposition, topology);
      append_static_prefix_diagnostics(
          mpi, decomposition, field_layout, exchange_plan, topology,
          geometry, execution, state.metadata().step,
          state.metadata().time_s, batch);
      append_static_suffix_diagnostics(
          mpi, boundaries,
          {fields.ids.face_mass_flux, topology.local_face_count(), true},
          state.metadata().step, state.metadata().time_s, batch);
      append_summary(driver_source, state.metadata().step,
                     state.metadata().time_s, mpi.rank(), "accepted-step",
                     batch);
      auto source =
          flow::checkpoint_v2_diagnostic_source(restored->report());
      append_summary(source, state.metadata().step,
                     state.metadata().time_s, mpi.rank(),
                     checkpoint_phase(restored->report().phase()), batch);
    } catch (const std::exception& error) {
      batch_ok = false;
      batch_message = error.what();
    }
    require_collective_success(
        mpi, batch_ok,
        batch_message.empty() ? "unable to construct diagnostic batch"
                              : batch_message);
    diagnostic_session->publish(mpi, state.metadata().step, batch);
  }
  if (config.density_model == DensityModel::constant) {
    auto facade = flow::FixedStepConstantDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&momentum_x, &momentum_y, &momentum_z},
        pressure_solver, pressure_preconditioner,
        constant_transport_specs(fields, config));
    run_loop(config, mpi, decomposition, topology, geometry, boundaries,
             fields, exchange_plan, execution, state, facade, controller,
             diagnostic_session ? &*diagnostic_session : nullptr,
             checkpoint_root);
  } else if (config.density_model == DensityModel::material) {
    auto facade = flow::FixedStepMaterialDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&momentum_x, &momentum_y, &momentum_z},
        pressure_solver, pressure_preconditioner, fields.registry, fields.ids,
        material_transport_spec(fields, config));
    run_loop(config, mpi, decomposition, topology, geometry, boundaries,
             fields, exchange_plan, execution, state, facade, controller,
             diagnostic_session ? &*diagnostic_session : nullptr,
             checkpoint_root);
  } else {
    const flow::IdealGasClosureSpec spec{
        fields.transported.front(),
        config.physics.cp_J_per_kg_K.value(),
        config.physics.gas_constant_J_per_kg_K.value(),
        config.physics.thermodynamic_pressure_pa.value()};
    auto closure =
        restored ? flow::IdealGasClosure::restore(
                       topology, geometry, boundaries, mpi, fields.registry,
                       fields.ids, state, spec,
                       restored->ideal_gas_closure_state())
                 : flow::IdealGasClosure::create(
                       topology, geometry, boundaries, mpi, fields.registry,
                       fields.ids, state, spec);
    auto facade = flow::FixedStepIdealGasFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&momentum_x, &momentum_y, &momentum_z},
        pressure_solver, pressure_preconditioner, fields.registry, fields.ids,
        material_transport_spec(fields, config), std::move(closure));
    run_loop(config, mpi, decomposition, topology, geometry, boundaries,
             fields, exchange_plan, execution, state, facade, controller,
             diagnostic_session ? &*diagnostic_session : nullptr,
             checkpoint_root);
  }
  print_finished(state, mpi);
  mpi.barrier();
  return EXIT_SUCCESS;
}

}  // namespace hundun::application
