// SPDX-License-Identifier: Apache-2.0

#include "applications/hundun/stage2_driver.hpp"
#include "applications/hundun/stage2_performance_build.hpp"
#include "applications/hundun/stage2_performance_crc64.hpp"
#include "applications/hundun/stage2_performance_diagnostic.hpp"
#include "applications/hundun/stage2_performance_fingerprint.hpp"
#include "applications/hundun/stage2_performance_timing.hpp"

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/diagnostics/checkpoint_v2_diagnostics.hpp"
#include "hundun/diagnostics/diagnostic_session.hpp"
#include "hundun/diagnostics/ideal_gas_closure_diagnostics.hpp"
#include "hundun/diagnostics/material_density_piso_diagnostics.hpp"
#include "hundun/diagnostics/mesh_diagnostic_v2.hpp"
#include "hundun/diagnostics/performance_artifact.hpp"
#include "hundun/diagnostics/performance_correctness.hpp"
#include "hundun/diagnostics/stage2_module_diagnostics.hpp"
#include "hundun/diagnostics/time_control_diagnostics.hpp"
#include "hundun/config/resolved_case_loader.hpp"
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

#include <algorithm>
#include <array>
#include <bitset>
#include <climits>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
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
  if (mpi.rank() == 0) {
    try {
      std::filesystem::create_directories(directory);
    } catch (const std::exception&) {
      local_ok = false;
    } catch (...) {
      local_ok = false;
    }
  }
  const auto status = runtime::collective_status(
      mpi, local_ok, "unable to create collective output directory");
  if (!status.ok)
    throw Error(status.message);
}

void require_collective_success(const MpiContext& mpi, bool local_ok,
                                std::string_view message) {
  const auto status = runtime::collective_status(mpi, local_ok, message);
  if (!status.ok)
    throw Error(status.message);
}

template <class Function>
void collective_transaction(const MpiContext& mpi, bool execute,
                            std::string_view description,
                            Function&& function) {
  bool local_ok = true;
  if (execute) {
    try {
      std::forward<Function>(function)();
    } catch (const std::exception&) {
      local_ok = false;
    } catch (...) {
      local_ok = false;
    }
  }
  require_collective_success(mpi, local_ok, description);
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
flow::TimeAdvanceReport advance_one(
    const config::FlowCaseConfig& config, flow::FlowState& state, Flow& flow,
    flow::Bdf2RetryController& controller) {
  const linear::SolveControl solve_control{};
  flow::TimeAdvanceReport advance = [&] {
    if constexpr (std::is_same_v<Flow, flow::FixedStepConstantDensityFlow>) {
      return controller.advance(
          state, flow, config.physics.rho_ref_kg_per_m3,
          config.physics.dynamic_viscosity_pa_s, solve_control, solve_control);
    } else {
      return controller.advance(state, flow,
                                config.physics.dynamic_viscosity_pa_s,
                                solve_control, solve_control);
    }
  }();
  if (advance.disposition() != flow::TimeAdvanceDisposition::committed) {
    throw Error("Stage 2 time advance failed");
  }
  return advance;
}

template <class Flow>
std::uint64_t execute_step_diagnostics(
    const config::FlowCaseConfig& config, const MpiContext& mpi,
    const runtime::StructuredDecomposition& decomposition,
    const mesh::MeshTopology& topology, const mesh::MeshGeometry& geometry,
    const boundary::BoundaryRegistry& boundaries,
    const DeclaredFields& fields, const runtime::ExchangePlan& exchange_plan,
    execution::ExecutionContext& execution, const flow::FlowState& state,
    Flow& flow, const flow::Bdf2RetryController& controller,
    const flow::TimeAdvanceReport& advance,
    const std::optional<flow::CheckpointV2Report>& checkpoint_report,
    diagnostics::DiagnosticSession* diagnostic_session) {
  return detail::execute_stage2_diagnostic_branch(
             diagnostic_session != nullptr,
             diagnostic_session != nullptr
                 ? static_cast<std::uint64_t>(
                       diagnostic_session->write_interval())
                 : 1U,
             state.metadata().step, [&]() -> std::uint64_t {
    DiagnosticBatch batch;
    bool batch_ok = true;
    try {
      const auto& final_attempt = advance.final_attempt();
      const auto& report = base_flow_report(final_attempt);
      diagnostics::FlowDriverDiagnosticSource driver_source{
          config.density_model,
          state.metadata().step,
          state.metadata().time_s,
          advance.attempt_count(),
          advance.disposition(),
          advance.reason(),
          advance.lowest_failing_rank()};
      auto field_layout = field_layout_source(fields, decomposition, topology);
      append_static_prefix_diagnostics(
          mpi, decomposition, field_layout, exchange_plan, topology, geometry,
          execution, state.metadata().step, state.metadata().time_s, batch);
      append_linear_solve_diagnostics(report, state.metadata().step,
                                      state.metadata().time_s, mpi.rank(),
                                      batch);
      append_static_suffix_diagnostics(
          mpi, boundaries,
          {fields.ids.face_mass_flux, topology.local_face_count(), true},
          state.metadata().step, state.metadata().time_s, batch);
      if constexpr (std::is_same_v<Flow,
                                   flow::FixedStepConstantDensityFlow>) {
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
        append_summary(source, state.metadata().step, state.metadata().time_s,
                       mpi.rank(), "material-density.attempt-result", batch);
      } else {
        const auto& ideal_report =
            std::get<flow::IdealGasStepAttemptReport>(final_attempt);
        append_summary(driver_source, state.metadata().step,
                       state.metadata().time_s, mpi.rank(), "accepted-step",
                       batch);
        {
          auto source = flow.flow_diagnostic_source(state, ideal_report);
          append_summary(source, state.metadata().step,
                         state.metadata().time_s, mpi.rank(),
                         "material-density.attempt-result", batch);
        }
        {
          auto source = flow.closure_diagnostic_source(state, ideal_report);
          append_summary(source, state.metadata().step,
                         state.metadata().time_s, mpi.rank(),
                         "ideal-gas-closure.attempt-result", batch);
        }
      }
      append_time_diagnostics(controller, state, advance, mpi.rank(), batch);
      if (checkpoint_report) {
        auto source =
            flow::checkpoint_v2_diagnostic_source(*checkpoint_report);
        append_summary(source, state.metadata().step, state.metadata().time_s,
                       mpi.rank(), checkpoint_phase(checkpoint_report->phase()),
                       batch);
      }
    } catch (const std::exception&) {
      batch_ok = false;
    } catch (...) {
      batch_ok = false;
    }
    require_collective_success(
        mpi, batch_ok, "unable to construct diagnostic batch");
    std::uint64_t logical_bytes = 0U;
    collective_transaction(
        mpi, true, "unable to serialize diagnostic batch", [&] {
          logical_bytes = static_cast<std::uint64_t>(
              batch.canonical_json_lines().size());
        });
    diagnostic_session->publish(mpi, state.metadata().step, batch);
    return logical_bytes;
  }).logical_bytes;
}

std::uint64_t fp64_bits(double value) noexcept;
#ifdef HUNDUN_APPLICATION_ENABLE_TEST_ACCESS
std::string local_state_fingerprint_bytes(
    const flow::FlowState& state,
    const flow::Bdf2RetryController& controller,
    const std::optional<flow::IdealGasClosureState>& closure);

struct Stage2DiagnosticObserverSnapshot final {
  std::array<std::uint64_t, 6> allocation{};
  std::array<std::uint64_t, 10> runtime_halo{};
  std::array<std::uint64_t, 10> pressure_halo{};
  std::array<std::uint64_t, 3> fp64{};
  std::array<std::array<std::uint64_t, 9>, 5> solves{};
  std::string state;
  std::vector<std::string> files;
};

std::array<std::uint64_t, 10> observer_halo(
    const runtime::HaloPerformanceCounters& counters) {
  return {counters.completed_exchanges,
          counters.begin_calls,
          counters.wait_calls,
          counters.send_payload_bytes,
          counters.receive_payload_bytes,
          counters.pack_bytes,
          counters.unpack_bytes,
          counters.send_messages,
          counters.receive_messages,
          fp64_bits(counters.completed_wait_seconds)};
}

template <class Flow>
runtime::HaloPerformanceCounters observer_pressure_halo(Flow& flow) {
  if constexpr (std::is_same_v<Flow, flow::FixedStepIdealGasFlow>) {
    return {};
  } else {
    return flow.pressure_halo_performance_counters();
  }
}

std::vector<std::string> observer_file_set(
    const std::filesystem::path& root) {
  std::vector<std::string> result;
  if (!std::filesystem::exists(root))
    return result;
  if (!std::filesystem::is_directory(root)) {
    result.push_back(".:" +
                     std::to_string(std::filesystem::file_size(root)));
    return result;
  }
  for (const auto& entry :
       std::filesystem::recursive_directory_iterator(root)) {
    const auto relative =
        std::filesystem::relative(entry.path(), root).generic_string();
    if (entry.is_regular_file()) {
      result.push_back(relative + ":" +
                       std::to_string(entry.file_size()));
    } else if (entry.is_directory()) {
      result.push_back(relative + "/");
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

std::array<std::array<std::uint64_t, 9>, 5> observer_solve_reports(
    const flow::TimeAdvanceReport& advance) {
  const auto& report = base_flow_report(advance.final_attempt());
  const std::array<const linear::SolveReport*, 5> solves{
      &report.momentum.components[0], &report.momentum.components[1],
      &report.momentum.components[2], &report.pressure[0],
      &report.pressure[1]};
  std::array<std::array<std::uint64_t, 9>, 5> result{};
  for (std::size_t index = 0; index < solves.size(); ++index) {
    const auto& solve = *solves[index];
    result[index] = {
        static_cast<std::uint64_t>(solve.reason),
        solve.iterations,
        solve.matvec_count,
        solve.preconditioner_apply_count,
        solve.global_reduction_count,
        fp64_bits(solve.initial_residual),
        fp64_bits(solve.recursive_residual),
        fp64_bits(solve.final_residual),
        static_cast<std::uint64_t>(
            static_cast<std::int64_t>(solve.lowest_failing_rank))};
  }
  return result;
}

template <class Flow>
Stage2DiagnosticObserverSnapshot capture_diagnostic_observer_snapshot(
    const MpiContext& mpi, const runtime::HaloExchange& halo,
    const flow::FlowState& state, Flow& flow,
    const flow::Bdf2RetryController& controller,
    const flow::TimeAdvanceReport& advance,
    const std::filesystem::path& diagnostic_directory) {
  const auto allocation = execution::allocation_counters();
  const auto fp64 = mpi.fp64_reduction_counters();
  return {{allocation.allocation_events,
           allocation.allocated_bytes,
           allocation.deallocation_events,
           allocation.deallocated_bytes,
           allocation.live_bytes,
          allocation.peak_live_bytes},
          observer_halo(halo.performance_counters()),
          observer_halo(observer_pressure_halo(flow)),
          {fp64.collective_calls, fp64.reduced_scalars,
           fp64.logical_payload_bytes},
          observer_solve_reports(advance),
          local_state_fingerprint_bytes(
              state, controller, closure_state_for_checkpoint(flow)),
          observer_file_set(diagnostic_directory)};
}

struct Stage2DiagnosticObserverComparison final {
  bool allocation{};
  bool runtime_halo{};
  bool pressure_halo{};
  bool fp64{};
  bool solves{};
  bool state{};
  bool files{};
};

Stage2DiagnosticObserverComparison compare_diagnostic_observer_snapshots(
    const Stage2DiagnosticObserverSnapshot& before,
    const Stage2DiagnosticObserverSnapshot& after) {
  return {before.allocation == after.allocation,
          before.runtime_halo == after.runtime_halo,
          before.pressure_halo == after.pressure_halo,
          before.fp64 == after.fp64,
          before.solves == after.solves,
          before.state == after.state,
          before.files == after.files};
}

void write_diagnostic_observer_record(
    const char* base, const MpiContext& mpi, std::string_view mode,
    const Stage2DiagnosticObserverComparison& comparison,
    std::uint64_t logical_bytes) {
  std::ostringstream suffix;
  suffix << ".rank-" << std::setw(6) << std::setfill('0') << mpi.rank()
         << ".jsonl";
  const std::filesystem::path path = std::string(base) + suffix.str();
  if (!path.parent_path().empty())
    std::filesystem::create_directories(path.parent_path());
  std::ofstream stream(path, std::ios::binary | std::ios::app);
  stream.imbue(std::locale::classic());
  stream << "{\"mode\":\"" << mode
         << "\",\"allocation_equal\":"
         << (comparison.allocation ? "true" : "false")
         << ",\"runtime_halo_equal\":"
         << (comparison.runtime_halo ? "true" : "false")
         << ",\"pressure_halo_equal\":"
         << (comparison.pressure_halo ? "true" : "false")
         << ",\"fp64_equal\":" << (comparison.fp64 ? "true" : "false")
         << ",\"solves_equal\":"
         << (comparison.solves ? "true" : "false")
         << ",\"state_equal\":"
         << (comparison.state ? "true" : "false")
         << ",\"files_equal\":"
         << (comparison.files ? "true" : "false")
         << ",\"logical_bytes\":" << logical_bytes << "}\n";
  stream.flush();
  if (!stream)
    throw Error("unable to write Task 25 diagnostic observer record");
}
#endif

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
    runtime::HaloExchange& halo,
    const std::filesystem::path& checkpoint_root,
    const std::filesystem::path& diagnostic_directory) {
#ifndef HUNDUN_APPLICATION_ENABLE_TEST_ACCESS
  static_cast<void>(halo);
  static_cast<void>(diagnostic_directory);
#endif
  while (state.metadata().step <
         static_cast<std::uint64_t>(config.time.steps)) {
    flow::TimeAdvanceReport advance =
        advance_one(config, state, flow, controller);

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

    std::uint64_t diagnostic_logical_bytes = 0U;
#ifdef HUNDUN_APPLICATION_ENABLE_TEST_ACCESS
    const char* observer = std::getenv("HUNDUN_TASK25_DIAGNOSTIC_OBSERVER");
    const bool observe = observer != nullptr && observer[0] != '\0';
    std::optional<Stage2DiagnosticObserverSnapshot> observer_before;
    if (observe) {
      observer_before = capture_diagnostic_observer_snapshot(
          mpi, halo, state, flow, controller, advance,
          diagnostic_directory);
    }
    try {
#endif
      diagnostic_logical_bytes = execute_step_diagnostics(
          config, mpi, decomposition, topology, geometry, boundaries, fields,
          exchange_plan, execution, state, flow, controller, advance,
          checkpoint_report, diagnostic_session);
#ifdef HUNDUN_APPLICATION_ENABLE_TEST_ACCESS
    } catch (...) {
      if (observe) {
        const auto after = capture_diagnostic_observer_snapshot(
            mpi, halo, state, flow, controller, advance,
            diagnostic_directory);
        const auto comparison =
            compare_diagnostic_observer_snapshots(*observer_before, after);
        write_diagnostic_observer_record(
            observer, mpi, "sink_failure", comparison, 0U);
        if (!comparison.allocation || !comparison.runtime_halo ||
            !comparison.pressure_halo || !comparison.fp64 ||
            !comparison.state || !comparison.solves || !comparison.files ||
            diagnostic_logical_bytes != 0U)
          throw Error("Task 25 diagnostic sink failure changed flow state");
      }
      throw;
    }
    if (observe) {
      const auto after = capture_diagnostic_observer_snapshot(
          mpi, halo, state, flow, controller, advance,
          diagnostic_directory);
      const auto comparison =
          compare_diagnostic_observer_snapshots(*observer_before, after);
      const bool disabled = diagnostic_session == nullptr;
      const bool due =
          !disabled &&
          state.metadata().step %
                  static_cast<std::uint64_t>(
                      diagnostic_session->write_interval()) ==
              0U;
      write_diagnostic_observer_record(
          observer, mpi, disabled ? "disabled" : (due ? "due" : "not_due"),
          comparison, diagnostic_logical_bytes);
      if ((!due && (!comparison.allocation || !comparison.runtime_halo ||
                    !comparison.pressure_halo || !comparison.fp64 ||
                    !comparison.solves || !comparison.state ||
                    !comparison.files || diagnostic_logical_bytes != 0U)) ||
          (due && (!comparison.allocation ||
                   !comparison.runtime_halo ||
                   !comparison.pressure_halo || !comparison.fp64 ||
                   !comparison.solves || !comparison.state ||
                   diagnostic_logical_bytes == 0U))) {
        throw Error("Task 25 diagnostic observer contract failed");
      }
    }
#endif
    static_cast<void>(diagnostic_logical_bytes);
    print_step(state, advance.attempt_count(), mpi);
  }
}

std::uint64_t fp64_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

std::uint64_t checked_add(std::uint64_t left, std::uint64_t right,
                          std::string_view description) {
  if (right > std::numeric_limits<std::uint64_t>::max() - left)
    throw Error(std::string(description) + " would overflow");
  return left + right;
}

std::uint64_t checked_multiply(std::uint64_t left, std::uint64_t right,
                               std::string_view description) {
  if (left != 0U &&
      right > std::numeric_limits<std::uint64_t>::max() / left)
    throw Error(std::string(description) + " would overflow");
  return left * right;
}

void validate_halo_counter_consistency(
    const MpiContext& mpi,
    const runtime::HaloPerformanceCounters& value,
    std::string_view family) {
  static_cast<void>(family);
  const bool consistent =
      value.completed_exchanges == value.begin_calls &&
      value.completed_exchanges == value.wait_calls &&
      value.pack_bytes == value.send_payload_bytes &&
      value.unpack_bytes == value.receive_payload_bytes &&
      std::isfinite(value.completed_wait_seconds) &&
      value.completed_wait_seconds >= 0.0;
  require_collective_success(
      mpi, consistent,
      "performance Halo counters are internally inconsistent");
}

bool performance_failure_injection(std::string_view selected) {
#ifdef HUNDUN_APPLICATION_ENABLE_TEST_ACCESS
  const char* value = std::getenv("HUNDUN_TASK25_PERFORMANCE_FAILURE");
  return value != nullptr && selected == value;
#else
  static_cast<void>(selected);
  return false;
#endif
}

template <class Map, std::size_t Size>
void validate_canonical_counter_keys(
    const Map& values, const std::array<std::string_view, Size>& expected) {
  if (values.size() != Size)
    throw Error("canonical performance counter key set differs");
  for (const auto key : expected) {
    if (values.find(std::string(key)) == values.end())
      throw Error("canonical performance counter key set differs");
  }
}

void validate_canonical_counter_maps(
    const diagnostics::ExactCounterMaps& counters) {
  validate_canonical_counter_keys(
      counters.allocated_bytes,
      std::array<std::string_view, 2>{"execution.allocated",
                                      "execution.peak-live"});
  validate_canonical_counter_keys(
      counters.halo_payload_bytes,
      std::array<std::string_view, 4>{"pack", "receive", "send", "unpack"});
  validate_canonical_counter_keys(
      counters.halo_messages,
      std::array<std::string_view, 2>{"receive", "send"});
  validate_canonical_counter_keys(
      counters.collectives,
      std::array<std::string_view, 3>{"checkpoint", "fp64-reduction",
                                      "linear-reduction"});
  validate_canonical_counter_keys(
      counters.collective_logical_payload_bytes,
      std::array<std::string_view, 1>{"fp64-reduction"});
  validate_canonical_counter_keys(
      counters.matvec,
      std::array<std::string_view, 2>{"momentum", "pressure"});
  validate_canonical_counter_keys(
      counters.preconditioner_applications,
      std::array<std::string_view, 2>{"momentum", "pressure"});
  validate_canonical_counter_keys(
      counters.logical_io_bytes,
      std::array<std::string_view, 2>{"checkpoint", "diagnostics"});
}

void append_u64_fingerprint(std::string& output, std::uint64_t value) {
  for (int shift = 56; shift >= 0; shift -= 8)
    output.push_back(static_cast<char>(
        (value >> static_cast<unsigned>(shift)) & UINT64_C(0xff)));
}

#ifdef HUNDUN_APPLICATION_ENABLE_TEST_ACCESS
std::string local_state_fingerprint_bytes(
    const flow::FlowState& state,
    const flow::Bdf2RetryController& controller,
    const std::optional<flow::IdealGasClosureState>& closure) {
  std::string encoding("hundun-performance-state-fp-v1");
  detail::append_performance_flow_layer(
      encoding, state.snapshot(flow::FlowLayer::history));
  detail::append_performance_flow_layer(
      encoding, state.snapshot(flow::FlowLayer::committed));
  detail::append_performance_flow_layer(
      encoding, state.snapshot(flow::FlowLayer::trial));
  const auto metadata = state.metadata();
  append_u64_fingerprint(encoding, metadata.step);
  append_u64_fingerprint(encoding, fp64_bits(metadata.time_s));
  append_u64_fingerprint(encoding, fp64_bits(metadata.dt_s));
  append_u64_fingerprint(encoding, fp64_bits(metadata.previous_dt_s));
  append_u64_fingerprint(
      encoding, static_cast<std::uint64_t>(metadata.order));
  const auto control = controller.state();
  append_u64_fingerprint(encoding, control.schema_version);
  append_u64_fingerprint(encoding, control.accepted_step);
  append_u64_fingerprint(encoding, fp64_bits(control.proposed_next_dt_s));
  append_u64_fingerprint(encoding, fp64_bits(control.last_accepted_dt_s));
  append_u64_fingerprint(
      encoding, static_cast<std::uint64_t>(control.last_accepted_order));
  append_u64_fingerprint(encoding, control.history_ready ? 1U : 0U);
  append_u64_fingerprint(
      encoding,
      control.last_all_linear_solves_within_half_limit ? 1U : 0U);
  append_u64_fingerprint(
      encoding, fp64_bits(control.last_convective_rate_per_s));
  append_u64_fingerprint(
      encoding, fp64_bits(control.last_diffusive_rate_per_s));
  append_u64_fingerprint(
      encoding, control.last_stability_metrics_available ? 1U : 0U);
  append_u64_fingerprint(encoding, control.last_retry_count);
  append_u64_fingerprint(encoding, control.revision);
  append_u64_fingerprint(encoding, control.state_seal);
  detail::append_performance_closure_state(encoding, closure);
  return encoding;
}
#endif

struct CanonicalStateRecord final {
  std::uint64_t layer{};
  std::uint64_t entity{};
  std::uint64_t field{};
  std::uint64_t nested{};
  std::uint64_t global_id{};
  std::uint64_t component{};
  std::uint64_t value_bits{};
};
static_assert(std::is_trivially_copyable_v<CanonicalStateRecord>);
static_assert(sizeof(CanonicalStateRecord) == 7U * sizeof(std::uint64_t));

bool canonical_state_record_less(const CanonicalStateRecord& left,
                                 const CanonicalStateRecord& right) noexcept {
  return std::tie(left.layer, left.entity, left.field, left.nested,
                  left.global_id, left.component) <
         std::tie(right.layer, right.entity, right.field, right.nested,
                  right.global_id, right.component);
}

std::vector<CanonicalStateRecord> local_canonical_state_records(
    const mesh::MeshTopology& topology, const flow::FlowState& state) {
  std::vector<CanonicalStateRecord> result;
  constexpr std::array<flow::FlowLayer, 3> layers{
      flow::FlowLayer::history, flow::FlowLayer::committed,
      flow::FlowLayer::trial};
  for (std::size_t layer_index = 0; layer_index < layers.size();
       ++layer_index) {
    const auto values = state.snapshot(layers[layer_index]);
    const auto append_cell_field =
        [&](std::uint64_t field, std::uint64_t nested,
            const std::vector<double>& source, std::size_t components) {
          if (source.size() != topology.owned_cell_count() * components)
            throw Error("performance cell field size is not canonical");
          for (mesh::LocalCellId cell = 0;
               cell < topology.owned_cell_count(); ++cell) {
            const auto global = topology.global_cell_id(cell);
            for (std::size_t component = 0; component < components;
                 ++component)
              result.push_back(
                  {static_cast<std::uint64_t>(layer_index), 0U, field, nested,
                   global, static_cast<std::uint64_t>(component),
                   fp64_bits(source[cell * components + component])});
          }
        };
    append_cell_field(0U, 0U, values.density, 1U);
    append_cell_field(1U, 0U, values.velocity, 3U);
    append_cell_field(2U, 0U, values.mechanical_pressure, 1U);
    for (std::size_t nested = 0;
         nested < values.transported_cell_fields.size(); ++nested)
      append_cell_field(3U, static_cast<std::uint64_t>(nested),
                        values.transported_cell_fields[nested], 1U);

    if (values.face_velocity.size() != topology.local_face_count() * 3U ||
        values.face_mass_flux.size() != topology.local_face_count())
      throw Error("performance face field size is not canonical");
    for (mesh::LocalFaceId face = 0; face < topology.local_face_count();
         ++face) {
      if (topology.face_ownership(face) !=
          mesh::EntityOwnership::owned)
        continue;
      const auto global = topology.global_face_id(face);
      for (std::size_t component = 0; component < 3U; ++component)
        result.push_back(
            {static_cast<std::uint64_t>(layer_index), 1U, 0U, 0U, global,
             static_cast<std::uint64_t>(component),
             fp64_bits(values.face_velocity[face * 3U + component])});
      result.push_back(
          {static_cast<std::uint64_t>(layer_index), 1U, 1U, 0U, global, 0U,
           fp64_bits(values.face_mass_flux[face])});
    }
  }
  return result;
}

std::size_t append_scalar_shape_fingerprint(
    std::string& encoding, const mesh::MeshTopology& topology,
    const flow::FlowState& state) {
  constexpr std::array<flow::FlowLayer, 3> layers{
      flow::FlowLayer::history, flow::FlowLayer::committed,
      flow::FlowLayer::trial};
  std::size_t words = 0U;
  for (const auto layer : layers) {
    const auto values = state.snapshot(layer);
    append_u64_fingerprint(
        encoding,
        static_cast<std::uint64_t>(
            values.transported_cell_fields.size()));
    ++words;
    for (const auto& nested : values.transported_cell_fields) {
      if (nested.size() != topology.owned_cell_count())
        throw Error(
            "performance transported scalar shape is not canonical");
      append_u64_fingerprint(encoding, topology.global_cell_count());
      ++words;
    }
  }
  return words;
}

std::size_t append_metadata_fingerprint(
    std::string& encoding, const flow::FlowState& state,
    const flow::Bdf2RetryController& controller,
    const std::optional<flow::IdealGasClosureState>& closure,
    const mesh::MeshTopology& topology) {
  const auto metadata = state.metadata();
  append_u64_fingerprint(encoding, metadata.step);
  append_u64_fingerprint(encoding, fp64_bits(metadata.time_s));
  append_u64_fingerprint(encoding, fp64_bits(metadata.dt_s));
  append_u64_fingerprint(encoding, fp64_bits(metadata.previous_dt_s));
  append_u64_fingerprint(encoding,
                         static_cast<std::uint64_t>(metadata.order));
  const auto control = controller.state();
  append_u64_fingerprint(encoding, control.schema_version);
  append_u64_fingerprint(encoding, control.accepted_step);
  append_u64_fingerprint(encoding, fp64_bits(control.proposed_next_dt_s));
  append_u64_fingerprint(encoding, fp64_bits(control.last_accepted_dt_s));
  append_u64_fingerprint(
      encoding, static_cast<std::uint64_t>(control.last_accepted_order));
  append_u64_fingerprint(encoding, control.history_ready ? 1U : 0U);
  append_u64_fingerprint(
      encoding,
      control.last_all_linear_solves_within_half_limit ? 1U : 0U);
  append_u64_fingerprint(
      encoding, fp64_bits(control.last_convective_rate_per_s));
  append_u64_fingerprint(
      encoding, fp64_bits(control.last_diffusive_rate_per_s));
  append_u64_fingerprint(
      encoding, control.last_stability_metrics_available ? 1U : 0U);
  append_u64_fingerprint(encoding, control.last_retry_count);
  append_u64_fingerprint(encoding, control.revision);
  append_u64_fingerprint(encoding, control.state_seal);
  detail::append_performance_closure_state(encoding, closure);
  const std::size_t closure_words =
      closure.has_value()
          ? (closure->target_mass_kg.has_value() ? 6U : 5U)
          : 1U;
  return 18U + closure_words +
         append_scalar_shape_fingerprint(encoding, topology, state);
}

std::string metadata_fingerprint_bytes(
    const flow::FlowState& state,
    const flow::Bdf2RetryController& controller,
    const std::optional<flow::IdealGasClosureState>& closure,
    const mesh::MeshTopology& topology) {
  std::string result;
  const auto words = append_metadata_fingerprint(
      result, state, controller, closure, topology);
  if (result.size() != words * sizeof(std::uint64_t))
    throw Error("performance metadata fingerprint size is malformed");
  return result;
}

std::string collective_state_fingerprint(
    const MpiContext& mpi, const mesh::MeshTopology& topology,
    const flow::FlowState& state,
    const flow::Bdf2RetryController& controller,
    const std::optional<flow::IdealGasClosureState>& closure) {
  std::vector<CanonicalStateRecord> local;
  bool local_ok = true;
  try {
    local = local_canonical_state_records(topology, state);
  } catch (const std::exception&) {
    local_ok = false;
  } catch (...) {
    local_ok = false;
  }
  require_collective_success(
      mpi, local_ok, "unable to encode local performance state");

  std::string local_metadata;
  local_ok = true;
  try {
    local_metadata =
        metadata_fingerprint_bytes(state, controller, closure, topology);
    if (mpi.rank() == 1 &&
        performance_failure_injection("nonroot_scalar_shape"))
      append_u64_fingerprint(local_metadata, 0U);
    if (mpi.rank() == 1 &&
        performance_failure_injection("nonroot_metadata"))
      local_metadata.back() =
          static_cast<char>(local_metadata.back() ^ static_cast<char>(1));
  } catch (const std::exception&) {
    local_ok = false;
  } catch (...) {
    local_ok = false;
  }
  require_collective_success(
      mpi, local_ok, "unable to encode performance metadata fingerprint");
  std::uint64_t metadata_size =
      static_cast<std::uint64_t>(local_metadata.size());
  std::uint64_t minimum_metadata_size = metadata_size;
  std::uint64_t maximum_metadata_size = metadata_size;
  require_collective_success(
      mpi,
      MPI_Allreduce(MPI_IN_PLACE, &minimum_metadata_size, 1, MPI_UINT64_T,
                    MPI_MIN, mpi.comm()) == MPI_SUCCESS,
      "unable to validate minimum performance metadata fingerprint size");
  require_collective_success(
      mpi,
      MPI_Allreduce(MPI_IN_PLACE, &maximum_metadata_size, 1, MPI_UINT64_T,
                    MPI_MAX, mpi.comm()) == MPI_SUCCESS,
      "unable to validate maximum performance metadata fingerprint size");
  require_collective_success(
      mpi, minimum_metadata_size == maximum_metadata_size,
      "performance scalar collection shape differs between ranks");
  require_collective_success(
      mpi, metadata_size <= static_cast<std::uint64_t>(INT_MAX),
      "performance metadata fingerprint exceeds MPI count range");
  std::vector<char> gathered_metadata;
  collective_transaction(
      mpi, mpi.rank() == 0,
      "unable to allocate gathered performance metadata", [&] {
        gathered_metadata.resize(
            local_metadata.size() * static_cast<std::size_t>(mpi.size()));
      });
  require_collective_success(
      mpi,
      MPI_Gather(local_metadata.data(),
                 static_cast<int>(local_metadata.size()), MPI_CHAR,
                 gathered_metadata.data(),
                 static_cast<int>(local_metadata.size()), MPI_CHAR, 0,
                 mpi.comm()) == MPI_SUCCESS,
      "unable to gather performance metadata");
  bool metadata_equal = true;
  collective_transaction(
      mpi, mpi.rank() == 0,
      "unable to compare gathered performance metadata", [&] {
        for (int rank = 1; rank < mpi.size(); ++rank) {
          const auto offset =
              static_cast<std::size_t>(rank) * local_metadata.size();
          metadata_equal =
              metadata_equal &&
              std::equal(
                  local_metadata.begin(), local_metadata.end(),
                  gathered_metadata.begin() +
                      static_cast<std::ptrdiff_t>(offset));
        }
      });
  require_collective_success(
      mpi, mpi.rank() != 0 || metadata_equal,
      "performance metadata/controller differs between ranks");

  std::uint64_t local_bytes_size = 0U;
  collective_transaction(
      mpi, true, "unable to derive local performance fingerprint byte count",
      [&] {
        local_bytes_size = checked_multiply(
            static_cast<std::uint64_t>(local.size()),
            static_cast<std::uint64_t>(sizeof(CanonicalStateRecord)),
            "local performance fingerprint byte count");
      });
  require_collective_success(
      mpi, local_bytes_size <= static_cast<std::uint64_t>(INT_MAX),
      "local performance fingerprint exceeds MPI count range");
  const int local_bytes = static_cast<int>(local_bytes_size);
  std::vector<int> byte_counts;
  local_ok = true;
  if (mpi.rank() == 0) {
    try {
      byte_counts.resize(static_cast<std::size_t>(mpi.size()));
    } catch (const std::exception&) {
      local_ok = false;
    } catch (...) {
      local_ok = false;
    }
  }
  require_collective_success(
      mpi, local_ok, "unable to allocate performance fingerprint counts");
  require_collective_success(
      mpi,
      MPI_Gather(&local_bytes, 1, MPI_INT, byte_counts.data(), 1, MPI_INT, 0,
                 mpi.comm()) == MPI_SUCCESS,
      "unable to gather performance fingerprint counts");

  std::vector<int> byte_offsets;
  std::vector<CanonicalStateRecord> gathered;
  local_ok = true;
  if (mpi.rank() == 0) {
    try {
      byte_offsets.resize(static_cast<std::size_t>(mpi.size()));
      std::uint64_t total_bytes = 0U;
      for (int rank = 0; rank < mpi.size(); ++rank) {
        const auto count = byte_counts[static_cast<std::size_t>(rank)];
        if (count < 0 ||
            total_bytes > static_cast<std::uint64_t>(INT_MAX) -
                              static_cast<std::uint64_t>(count))
          throw Error("performance fingerprint gather exceeds MPI range");
        byte_offsets[static_cast<std::size_t>(rank)] =
            static_cast<int>(total_bytes);
        total_bytes += static_cast<std::uint64_t>(count);
      }
      if (total_bytes % sizeof(CanonicalStateRecord) != 0U)
        throw Error("performance fingerprint record bytes are malformed");
      gathered.resize(
          static_cast<std::size_t>(total_bytes / sizeof(CanonicalStateRecord)));
    } catch (const std::exception&) {
      local_ok = false;
    } catch (...) {
      local_ok = false;
    }
  }
  require_collective_success(
      mpi, local_ok,
      "unable to allocate gathered performance fingerprint");
  require_collective_success(
      mpi,
      MPI_Gatherv(local.data(), local_bytes, MPI_BYTE, gathered.data(),
                  byte_counts.data(), byte_offsets.data(), MPI_BYTE, 0,
                  mpi.comm()) == MPI_SUCCESS,
      "unable to gather performance fingerprint records");

  std::string fingerprint;
  local_ok = true;
  if (mpi.rank() == 0) {
    try {
      std::sort(gathered.begin(), gathered.end(),
                canonical_state_record_less);
      for (std::size_t index = 1; index < gathered.size(); ++index) {
        const auto& left = gathered[index - 1U];
        const auto& right = gathered[index];
        if (!canonical_state_record_less(left, right))
          throw Error("performance fingerprint contains duplicate entity");
      }
      std::string encoding("hundun-performance-state-fp-v1");
      for (const auto& record : gathered) {
        append_u64_fingerprint(encoding, record.layer);
        append_u64_fingerprint(encoding, record.entity);
        append_u64_fingerprint(encoding, record.field);
        append_u64_fingerprint(encoding, record.nested);
        append_u64_fingerprint(encoding, record.global_id);
        append_u64_fingerprint(encoding, record.component);
        append_u64_fingerprint(encoding, record.value_bits);
      }
      encoding.append(local_metadata);
      fingerprint = detail::tagged_performance_crc64(
          "hundun-performance-state-fp-v1", encoding);
    } catch (const std::exception&) {
      local_ok = false;
    } catch (...) {
      local_ok = false;
    }
  }
  require_collective_success(
      mpi, local_ok, "unable to finalize performance fingerprint");

  std::uint64_t length = static_cast<std::uint64_t>(fingerprint.size());
  require_collective_success(
      mpi, MPI_Bcast(&length, 1, MPI_UINT64_T, 0, mpi.comm()) == MPI_SUCCESS,
      "unable to broadcast performance fingerprint length");
  require_collective_success(
      mpi, length <= static_cast<std::uint64_t>(INT_MAX),
      "performance fingerprint exceeds MPI count range");
  local_ok = true;
  try {
    fingerprint.resize(static_cast<std::size_t>(length));
  } catch (const std::exception&) {
    local_ok = false;
  } catch (...) {
    local_ok = false;
  }
  require_collective_success(
      mpi, local_ok, "unable to allocate performance fingerprint result");
  require_collective_success(
      mpi,
      MPI_Bcast(fingerprint.data(), static_cast<int>(length), MPI_CHAR, 0,
                mpi.comm()) == MPI_SUCCESS,
      "unable to broadcast performance fingerprint");
  return fingerprint;
}

std::string termination_token(linear::SolveTerminationReason reason) {
  switch (reason) {
  case linear::SolveTerminationReason::converged:
    return "converged";
  case linear::SolveTerminationReason::zero_right_hand_side:
    return "zero_right_hand_side";
  case linear::SolveTerminationReason::maximum_iterations:
    return "maximum_iterations";
  case linear::SolveTerminationReason::numerical_breakdown:
    return "numerical_breakdown";
  case linear::SolveTerminationReason::non_finite_value:
    return "non_finite_value";
  case linear::SolveTerminationReason::invalid_control:
    return "invalid_control";
  case linear::SolveTerminationReason::collective_failure:
    return "collective_failure";
  }
  throw Error("invalid linear solve termination reason");
}

std::array<diagnostics::PerformanceWorkRecord, 5>
performance_work_records(const flow::TimeAdvanceReport& advance,
                         std::uint64_t repetition, char phase,
                         std::uint64_t relative_step) {
  const auto& report = base_flow_report(advance.final_attempt());
  std::array<const linear::SolveReport*, 5> solves{
      &report.momentum.components[0], &report.momentum.components[1],
      &report.momentum.components[2], &report.pressure[0],
      &report.pressure[1]};
  std::array<diagnostics::PerformanceWorkRecord, 5> result;
  for (std::size_t slot = 0; slot < solves.size(); ++slot) {
    const auto& solve = *solves[slot];
    if ((solve.reason != linear::SolveTerminationReason::converged &&
         solve.reason != linear::SolveTerminationReason::zero_right_hand_side) ||
        !std::isfinite(solve.initial_residual) ||
        !std::isfinite(solve.recursive_residual) ||
        !std::isfinite(solve.final_residual)) {
      throw Error("performance solve report failed its residual contract");
    }
    result[slot] = {
        repetition,
        phase,
        relative_step,
        static_cast<std::uint64_t>(slot),
        termination_token(solve.reason),
        solve.iterations,
        solve.matvec_count,
        solve.preconditioner_apply_count,
        solve.global_reduction_count,
        fp64_bits(solve.initial_residual),
        fp64_bits(solve.recursive_residual),
        fp64_bits(solve.final_residual)};
  }
  return result;
}

std::uint64_t pack_work_termination(std::string_view value) {
  constexpr std::array<std::string_view, 7> tokens{
      "converged",           "zero_right_hand_side",
      "maximum_iterations",  "numerical_breakdown",
      "non_finite_value",    "invalid_control",
      "collective_failure"};
  const auto found = std::find(tokens.begin(), tokens.end(), value);
  if (found == tokens.end())
    throw Error("performance work termination is not canonical");
  return static_cast<std::uint64_t>(
      std::distance(tokens.begin(), found));
}

std::array<std::uint64_t, 12> pack_work(
    const diagnostics::PerformanceWorkRecord& work) {
  if (work.phase != 'W' && work.phase != 'M')
    throw Error("performance work phase is not canonical");
  return {work.repetition,
          static_cast<std::uint64_t>(work.phase == 'W' ? 0U : 1U),
          work.relative_step,
          work.slot,
          pack_work_termination(work.termination),
          work.iterations,
          work.matvec,
          work.preconditioner,
          work.reduction,
          work.initial_residual_bits,
          work.recursive_residual_bits,
          work.independent_final_residual_bits};
}

void validate_collective_work(
    const MpiContext& mpi,
    const std::vector<diagnostics::PerformanceWorkRecord>& local) {
  std::uint64_t minimum_count =
      static_cast<std::uint64_t>(local.size());
  std::uint64_t maximum_count = minimum_count;
  require_collective_success(
      mpi,
      MPI_Allreduce(MPI_IN_PLACE, &minimum_count, 1, MPI_UINT64_T, MPI_MIN,
                    mpi.comm()) == MPI_SUCCESS,
      "unable to validate minimum performance work count");
  require_collective_success(
      mpi,
      MPI_Allreduce(MPI_IN_PLACE, &maximum_count, 1, MPI_UINT64_T, MPI_MAX,
                    mpi.comm()) == MPI_SUCCESS,
      "unable to validate maximum performance work count");
  require_collective_success(
      mpi, minimum_count == maximum_count,
      "performance work report counts differ between ranks");
  for (const auto& work : local) {
    const auto packed = pack_work(work);
    std::vector<std::uint64_t> gathered;
    collective_transaction(
        mpi, mpi.rank() == 0,
        "unable to allocate gathered performance work reports", [&] {
          gathered.resize(packed.size() *
                          static_cast<std::size_t>(mpi.size()));
        });
    require_collective_success(
        mpi,
        MPI_Gather(packed.data(), static_cast<int>(packed.size()),
                   MPI_UINT64_T, gathered.data(),
                   static_cast<int>(packed.size()), MPI_UINT64_T, 0,
                   mpi.comm()) == MPI_SUCCESS,
        "unable to gather performance work reports");
    bool equal = true;
    collective_transaction(
        mpi, mpi.rank() == 0,
        "unable to validate gathered performance work reports", [&] {
          for (int rank = 1; rank < mpi.size(); ++rank) {
            const auto offset =
                static_cast<std::size_t>(rank) * packed.size();
            equal = equal &&
                    std::equal(packed.begin(), packed.end(),
                               gathered.begin() +
                                   static_cast<std::ptrdiff_t>(offset));
          }
        });
    require_collective_success(
        mpi, mpi.rank() != 0 || equal,
        "performance work reports differ between ranks");
  }
}

runtime::HaloPerformanceCounters halo_delta(
    runtime::HaloPerformanceCounters before,
    runtime::HaloPerformanceCounters after) {
  const auto subtract = [](std::uint64_t left, std::uint64_t right) {
    if (left < right)
      throw Error("Halo performance counter regressed");
    return left - right;
  };
  runtime::HaloPerformanceCounters result{
      subtract(after.completed_exchanges, before.completed_exchanges),
      subtract(after.begin_calls, before.begin_calls),
      subtract(after.wait_calls, before.wait_calls),
      subtract(after.send_payload_bytes, before.send_payload_bytes),
      subtract(after.receive_payload_bytes, before.receive_payload_bytes),
      subtract(after.pack_bytes, before.pack_bytes),
      subtract(after.unpack_bytes, before.unpack_bytes),
      subtract(after.send_messages, before.send_messages),
      subtract(after.receive_messages, before.receive_messages),
      after.completed_wait_seconds - before.completed_wait_seconds};
  if (!std::isfinite(result.completed_wait_seconds) ||
      result.completed_wait_seconds < 0.0)
    throw Error("Halo performance wait counter regressed");
  return result;
}

runtime::Fp64ReductionCounters fp64_delta(
    runtime::Fp64ReductionCounters before,
    runtime::Fp64ReductionCounters after) {
  if (after.collective_calls < before.collective_calls ||
      after.reduced_scalars < before.reduced_scalars ||
      after.logical_payload_bytes < before.logical_payload_bytes)
    throw Error("FP64 reduction counter regressed");
  return {after.collective_calls - before.collective_calls,
          after.reduced_scalars - before.reduced_scalars,
          after.logical_payload_bytes - before.logical_payload_bytes};
}

std::uint64_t sum_rank_u64(const MpiContext& mpi, std::uint64_t local,
                           std::string_view description) {
  std::vector<std::uint64_t> gathered;
  collective_transaction(
      mpi, mpi.rank() == 0,
      "unable to allocate gathered performance counter", [&] {
        gathered.resize(static_cast<std::size_t>(mpi.size()));
      });
  require_collective_success(
      mpi,
      MPI_Gather(&local, 1, MPI_UINT64_T, gathered.data(), 1, MPI_UINT64_T, 0,
                 mpi.comm()) == MPI_SUCCESS,
      "unable to gather performance counter");
  std::uint64_t sum = 0U;
  collective_transaction(mpi, mpi.rank() == 0, description, [&] {
    for (const auto value : gathered)
      sum = checked_add(sum, value, description);
  });
  return sum;
}

std::uint64_t common_rank_u64(const MpiContext& mpi, std::uint64_t local,
                              std::string_view description) {
  std::vector<std::uint64_t> gathered;
  collective_transaction(
      mpi, mpi.rank() == 0,
      "unable to allocate gathered communicator counter", [&] {
        gathered.resize(static_cast<std::size_t>(mpi.size()));
      });
  require_collective_success(
      mpi,
      MPI_Gather(&local, 1, MPI_UINT64_T, gathered.data(), 1, MPI_UINT64_T, 0,
                 mpi.comm()) == MPI_SUCCESS,
      "unable to gather performance communicator counter");
  bool equal = true;
  collective_transaction(mpi, mpi.rank() == 0, description, [&] {
    equal = std::all_of(gathered.begin(), gathered.end(),
                        [local](std::uint64_t value) {
                          return value == local;
                        });
  });
  require_collective_success(mpi, mpi.rank() != 0 || equal, description);
  return mpi.rank() == 0 ? local : 0U;
}

struct PerformanceTotals final {
  std::uint64_t assembly_live{};
  runtime::HaloPerformanceCounters halo;
  runtime::Fp64ReductionCounters fp64;
  std::uint64_t checkpoint_collectives{};
  std::uint64_t checkpoint_logical_bytes{};
  std::uint64_t checkpoint_actual_bytes{};
  std::uint64_t diagnostic_logical_bytes{};
  std::uint64_t diagnostic_actual_bytes{};
  double checkpoint_io_seconds{};
  double diagnostic_io_seconds{};
};

struct HaloEvidenceRecord final {
  std::uint64_t repetition{};
  std::uint64_t relative_rank{};
  runtime::HaloPerformanceCounters runtime;
  runtime::HaloPerformanceCounters pressure;
};

std::array<std::uint64_t, 20> pack_halo_evidence(
    const runtime::HaloPerformanceCounters& runtime,
    const runtime::HaloPerformanceCounters& pressure) {
  const auto pack = [](std::array<std::uint64_t, 20>& result,
                       std::size_t offset,
                       const runtime::HaloPerformanceCounters& value) {
    result[offset + 0U] = value.completed_exchanges;
    result[offset + 1U] = value.begin_calls;
    result[offset + 2U] = value.wait_calls;
    result[offset + 3U] = value.send_payload_bytes;
    result[offset + 4U] = value.receive_payload_bytes;
    result[offset + 5U] = value.pack_bytes;
    result[offset + 6U] = value.unpack_bytes;
    result[offset + 7U] = value.send_messages;
    result[offset + 8U] = value.receive_messages;
    result[offset + 9U] = fp64_bits(value.completed_wait_seconds);
  };
  std::array<std::uint64_t, 20> result{};
  pack(result, 0U, runtime);
  pack(result, 10U, pressure);
  return result;
}

runtime::HaloPerformanceCounters unpack_halo_evidence(
    const std::uint64_t* values) {
  double wait_seconds = 0.0;
  std::memcpy(&wait_seconds, values + 9U, sizeof(wait_seconds));
  return {values[0U], values[1U], values[2U], values[3U], values[4U],
          values[5U], values[6U], values[7U], values[8U], wait_seconds};
}

void gather_halo_evidence(
    const MpiContext& mpi, std::uint64_t repetition,
    const runtime::HaloPerformanceCounters& runtime,
    const runtime::HaloPerformanceCounters& pressure,
    std::vector<HaloEvidenceRecord>& records) {
  const auto packed = pack_halo_evidence(runtime, pressure);
  std::vector<std::uint64_t> gathered;
  collective_transaction(
      mpi, mpi.rank() == 0,
      "unable to allocate gathered Halo evidence", [&] {
        gathered.resize(packed.size() *
                        static_cast<std::size_t>(mpi.size()));
      });
  require_collective_success(
      mpi,
      MPI_Gather(packed.data(), static_cast<int>(packed.size()),
                 MPI_UINT64_T, gathered.data(),
                 static_cast<int>(packed.size()), MPI_UINT64_T, 0,
                 mpi.comm()) == MPI_SUCCESS,
      "unable to gather Halo evidence");
  collective_transaction(
      mpi, mpi.rank() == 0, "unable to retain gathered Halo evidence", [&] {
        for (int rank = 0; rank < mpi.size(); ++rank) {
          const auto offset =
              static_cast<std::size_t>(rank) * packed.size();
          records.push_back(
              {repetition, static_cast<std::uint64_t>(rank),
               unpack_halo_evidence(gathered.data() + offset),
               unpack_halo_evidence(gathered.data() + offset + 10U)});
        }
      });
}

double maximum_rank_seconds(const MpiContext& mpi, double local,
                            std::string_view description) {
  require_collective_success(
      mpi, std::isfinite(local) && local >= 0.0,
      "performance I/O duration contains an invalid rank sample");
  double maximum = 0.0;
  require_collective_success(
      mpi,
      MPI_Reduce(&local, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0,
                 mpi.comm()) == MPI_SUCCESS,
      description);
  return mpi.rank() == 0 ? maximum : 0.0;
}

void append_json_double(std::ostringstream& output, double value) {
  if (!std::isfinite(value) || value < 0.0)
    throw Error("performance evidence FP64 value is invalid");
  output << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
}

runtime::HaloPerformanceCounters combine_halo_evidence(
    const runtime::HaloPerformanceCounters& left,
    const runtime::HaloPerformanceCounters& right) {
  runtime::HaloPerformanceCounters result;
  result.completed_exchanges = checked_add(
      left.completed_exchanges, right.completed_exchanges,
      "evidence completed exchanges");
  result.begin_calls =
      checked_add(left.begin_calls, right.begin_calls,
                  "evidence begin calls");
  result.wait_calls =
      checked_add(left.wait_calls, right.wait_calls,
                  "evidence wait calls");
  result.send_payload_bytes = checked_add(
      left.send_payload_bytes, right.send_payload_bytes,
      "evidence send payload");
  result.receive_payload_bytes = checked_add(
      left.receive_payload_bytes, right.receive_payload_bytes,
      "evidence receive payload");
  result.pack_bytes = checked_add(
      left.pack_bytes, right.pack_bytes, "evidence pack bytes");
  result.unpack_bytes = checked_add(
      left.unpack_bytes, right.unpack_bytes, "evidence unpack bytes");
  result.send_messages = checked_add(
      left.send_messages, right.send_messages, "evidence send messages");
  result.receive_messages = checked_add(
      left.receive_messages, right.receive_messages,
      "evidence receive messages");
  result.completed_wait_seconds =
      left.completed_wait_seconds + right.completed_wait_seconds;
  if (!std::isfinite(result.completed_wait_seconds) ||
      result.completed_wait_seconds < 0.0)
    throw Error("evidence wait seconds are invalid");
  return result;
}

void append_halo_evidence_json(
    std::ostringstream& output,
    const runtime::HaloPerformanceCounters& value) {
  output << "{\"completed_exchanges\":" << value.completed_exchanges
         << ",\"begin_calls\":" << value.begin_calls
         << ",\"wait_calls\":" << value.wait_calls
         << ",\"send_payload_bytes\":" << value.send_payload_bytes
         << ",\"receive_payload_bytes\":" << value.receive_payload_bytes
         << ",\"pack_bytes\":" << value.pack_bytes
         << ",\"unpack_bytes\":" << value.unpack_bytes
         << ",\"send_messages\":" << value.send_messages
         << ",\"receive_messages\":" << value.receive_messages
         << ",\"successful_wait_seconds\":";
  append_json_double(output, value.completed_wait_seconds);
  output << '}';
}

void append_value_metric(std::ostringstream& output,
                         std::string_view unit, std::string_view status,
                         double value) {
  output << "{\"unit\":\"" << unit << "\",\"status\":\"" << status
         << "\",\"value\":";
  if (status == "not_applicable") {
    output << "null";
  } else {
    append_json_double(output, value);
  }
  output << '}';
}

std::string serialize_performance_evidence_input(
    const std::vector<HaloEvidenceRecord>& records,
    const PerformanceTotals& totals, int repetitions, int ranks) {
  if (records.size() !=
      static_cast<std::size_t>(repetitions) *
          static_cast<std::size_t>(ranks))
    throw Error("performance Halo evidence coverage is incomplete");
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << "{\"schema_version\":1,\"entries\":[";
  for (std::size_t index = 0; index < records.size(); ++index) {
    const auto& record = records[index];
    if (index != 0U)
      output << ',';
    output << "{\"repetition\":" << record.repetition
           << ",\"relative_rank\":" << record.relative_rank
           << ",\"runtime\":";
    append_halo_evidence_json(output, record.runtime);
    output << ",\"pressure\":";
    append_halo_evidence_json(output, record.pressure);
    output << ",\"combined\":";
    append_halo_evidence_json(
        output, combine_halo_evidence(record.runtime, record.pressure));
    output << '}';
  }
  output << "],\"halo_summaries\":[";
  bool first_summary = true;
  constexpr std::array<std::string_view, 3> families{
      "runtime", "pressure", "combined"};
  for (int repetition = 0; repetition < repetitions; ++repetition) {
    for (std::size_t family = 0; family < families.size(); ++family) {
      std::uint64_t receive_payload = 0U;
      double critical_wait = 0.0;
      for (const auto& record : records) {
        if (record.repetition !=
            static_cast<std::uint64_t>(repetition))
          continue;
        const auto value =
            family == 0U
                ? record.runtime
                : (family == 1U
                       ? record.pressure
                       : combine_halo_evidence(
                             record.runtime, record.pressure));
        receive_payload = checked_add(
            receive_payload, value.receive_payload_bytes,
            "evidence summary receive payload");
        critical_wait =
            std::max(critical_wait, value.completed_wait_seconds);
      }
      const bool applicable = receive_payload != 0U;
      if (applicable &&
          (!std::isfinite(critical_wait) || critical_wait <= 0.0))
        throw Error(
            "nonlocal performance evidence requires positive critical wait");
      if (!first_summary)
        output << ',';
      first_summary = false;
      output << "{\"repetition\":" << repetition
             << ",\"family\":\"" << families[family]
             << "\",\"critical_wait_seconds\":";
      append_value_metric(output, "s",
                          applicable ? "available" : "not_applicable",
                          critical_wait);
      output << ",\"effective_bandwidth_bytes_per_second\":";
      append_value_metric(
          output, "byte/s",
          applicable ? "available" : "not_applicable",
          applicable ? static_cast<double>(receive_payload) / critical_wait
                     : 0.0);
      output << '}';
    }
  }
  output << "],\"io\":{\"checkpoint\":{\"logical_bytes\":"
         << totals.checkpoint_logical_bytes
         << ",\"actual_bytes\":" << totals.checkpoint_actual_bytes
         << ",\"duration_seconds\":";
  const bool checkpoint_applicable =
      totals.checkpoint_actual_bytes != 0U ||
      totals.checkpoint_logical_bytes != 0U;
  if (checkpoint_applicable && totals.checkpoint_io_seconds <= 0.0)
    throw Error("checkpoint evidence requires positive I/O duration");
  append_value_metric(
      output, "s",
      checkpoint_applicable ? "available" : "not_applicable",
      totals.checkpoint_io_seconds);
  output << ",\"throughput_bytes_per_second\":";
  append_value_metric(
      output, "byte/s",
      checkpoint_applicable ? "available" : "not_applicable",
      checkpoint_applicable
          ? static_cast<double>(totals.checkpoint_actual_bytes) /
                totals.checkpoint_io_seconds
          : 0.0);
  output << "},\"diagnostics\":{\"logical_bytes\":"
         << totals.diagnostic_logical_bytes
         << ",\"actual_bytes\":" << totals.diagnostic_actual_bytes
         << ",\"duration_seconds\":";
  const bool diagnostic_applicable =
      totals.diagnostic_actual_bytes != 0U ||
      totals.diagnostic_logical_bytes != 0U;
  if (diagnostic_applicable && totals.diagnostic_io_seconds <= 0.0)
    throw Error("diagnostic evidence requires positive I/O duration");
  append_value_metric(
      output, "s",
      diagnostic_applicable ? "available" : "not_applicable",
      totals.diagnostic_io_seconds);
  output << ",\"throughput_bytes_per_second\":";
  append_value_metric(
      output, "byte/s",
      diagnostic_applicable ? "available" : "not_applicable",
      diagnostic_applicable
          ? static_cast<double>(totals.diagnostic_actual_bytes) /
                totals.diagnostic_io_seconds
          : 0.0);
  output << "}}}";
  return output.str();
}

void add_halo_totals(PerformanceTotals& totals,
                     runtime::HaloPerformanceCounters value) {
  totals.halo.completed_exchanges =
      checked_add(totals.halo.completed_exchanges,
                  value.completed_exchanges, "Halo completed exchanges");
  totals.halo.begin_calls =
      checked_add(totals.halo.begin_calls, value.begin_calls,
                  "Halo begin calls");
  totals.halo.wait_calls =
      checked_add(totals.halo.wait_calls, value.wait_calls,
                  "Halo wait calls");
  totals.halo.send_payload_bytes =
      checked_add(totals.halo.send_payload_bytes, value.send_payload_bytes,
                  "Halo send payload");
  totals.halo.receive_payload_bytes =
      checked_add(totals.halo.receive_payload_bytes,
                  value.receive_payload_bytes, "Halo receive payload");
  totals.halo.pack_bytes =
      checked_add(totals.halo.pack_bytes, value.pack_bytes, "Halo pack bytes");
  totals.halo.unpack_bytes = checked_add(
      totals.halo.unpack_bytes, value.unpack_bytes, "Halo unpack bytes");
  totals.halo.send_messages =
      checked_add(totals.halo.send_messages, value.send_messages,
                  "Halo send messages");
  totals.halo.receive_messages =
      checked_add(totals.halo.receive_messages, value.receive_messages,
                  "Halo receive messages");
  totals.halo.completed_wait_seconds += value.completed_wait_seconds;
  if (!std::isfinite(totals.halo.completed_wait_seconds))
    throw Error("Halo completed wait seconds would overflow");
}

void atomic_publish_performance(const MpiContext& mpi,
                                const std::filesystem::path& case_root,
                                const std::filesystem::path& relative_directory,
                                std::string_view json,
                                std::string_view evidence_input) {
  std::filesystem::path directory;
  bool path_ok = true;
  try {
    if (mpi.rank() == 1 &&
        performance_failure_injection("nonroot_path"))
      throw Error("injected performance path derivation failure");
    directory = beneath(case_root, relative_directory);
  } catch (...) {
    path_ok = false;
  }
  require_collective_success(
      mpi, path_ok, "unable to derive performance artifact path");

  bool local_ok = true;
  std::filesystem::path artifact_temporary;
  std::filesystem::path evidence_temporary;
  std::filesystem::path evidence_final;
  std::filesystem::path evidence_backup;
  bool evidence_backup_created = false;
  bool evidence_committed = false;
  if (mpi.rank() == 0) {
    try {
      artifact_temporary = directory / "performance.v1.json.tmp";
      const auto artifact_final = directory / "performance.v1.json";
      evidence_temporary =
          directory / "performance-evidence-input.v1.json.tmp";
      evidence_final =
          directory / "performance-evidence-input.v1.json";
      evidence_backup =
          directory / "performance-evidence-input.v1.json.backup";
      std::filesystem::create_directories(directory);
      if (performance_failure_injection("stage"))
        throw Error("injected performance artifact stage failure");
      {
        std::ofstream stream(artifact_temporary,
                             std::ios::binary | std::ios::trunc);
        stream.write(json.data(), static_cast<std::streamsize>(json.size()));
        stream.flush();
        if (!stream)
          throw Error("unable to stage performance artifact");
      }
      if (performance_failure_injection("evidence_stage"))
        throw Error("injected performance evidence stage failure");
      {
        std::ofstream stream(evidence_temporary,
                             std::ios::binary | std::ios::trunc);
        stream.write(evidence_input.data(),
                     static_cast<std::streamsize>(evidence_input.size()));
        stream.flush();
        if (!stream)
          throw Error("unable to stage performance evidence input");
      }
      if (performance_failure_injection("evidence_rename"))
        throw Error("injected performance evidence rename failure");
      std::error_code ignored;
      std::filesystem::remove(evidence_backup, ignored);
      if (std::filesystem::exists(evidence_final)) {
        std::filesystem::rename(evidence_final, evidence_backup);
        evidence_backup_created = true;
      }
      std::filesystem::rename(evidence_temporary, evidence_final);
      evidence_committed = true;
      if (performance_failure_injection("rename"))
        throw Error("injected performance artifact rename failure");
      std::filesystem::rename(artifact_temporary, artifact_final);
      if (evidence_backup_created) {
        std::error_code cleanup_error;
        std::filesystem::remove(evidence_backup, cleanup_error);
      }
    } catch (...) {
      local_ok = false;
      std::error_code ignored;
      std::filesystem::remove(artifact_temporary, ignored);
      std::filesystem::remove(evidence_temporary, ignored);
      if (evidence_committed)
        std::filesystem::remove(evidence_final, ignored);
      if (evidence_backup_created)
        std::filesystem::rename(evidence_backup, evidence_final, ignored);
    }
  }
  require_collective_success(
      mpi, local_ok, "unable to publish performance artifact");
}

std::string mpi_identity() {
  std::array<char, MPI_MAX_LIBRARY_VERSION_STRING> buffer{};
  int length = 0;
  if (MPI_Get_library_version(buffer.data(), &length) != MPI_SUCCESS ||
      length <= 0)
    return "unavailable";
  std::string result(buffer.data(), static_cast<std::size_t>(length));
  for (char& character : result)
    if (character == '\n' || character == '\r' || character == '\t')
      character = ' ';
  while (!result.empty() &&
         (result.back() == ' ' || result.back() == '\0'))
    result.pop_back();
  return result.empty() ? "unavailable" : result;
}

void broadcast_root_string(const MpiContext& mpi, std::string& value,
                           std::string_view description) {
  std::uint64_t length = static_cast<std::uint64_t>(value.size());
  require_collective_success(
      mpi, MPI_Bcast(&length, 1, MPI_UINT64_T, 0, mpi.comm()) == MPI_SUCCESS,
      "performance metadata string length broadcast failed");
  require_collective_success(
      mpi, length <= static_cast<std::uint64_t>(INT_MAX),
      "performance metadata string exceeds MPI count range");
  collective_transaction(mpi, mpi.rank() != 0, description, [&] {
    value.resize(static_cast<std::size_t>(length));
  });
  require_collective_success(
      mpi,
      MPI_Bcast(value.data(), static_cast<int>(length), MPI_CHAR, 0,
                mpi.comm()) == MPI_SUCCESS,
      "performance metadata string broadcast failed");
}

struct PerformancePlacement final {
  std::string node_identity;
  std::string rank_placement;
};

PerformancePlacement performance_placement(const MpiContext& mpi) {
  std::array<char, MPI_MAX_PROCESSOR_NAME> local_name{};
  int local_name_length = 0;
  const bool processor_name_ok =
      MPI_Get_processor_name(local_name.data(), &local_name_length) ==
          MPI_SUCCESS &&
      local_name_length > 0 &&
      local_name_length < static_cast<int>(local_name.size());
  if (!processor_name_ok) {
    local_name.fill('\0');
    constexpr std::string_view unavailable = "unavailable";
    std::copy(unavailable.begin(), unavailable.end(), local_name.begin());
  } else {
    local_name[static_cast<std::size_t>(local_name_length)] = '\0';
  }
  require_collective_success(
      mpi, true, "unable to materialize local processor name");
  std::vector<char> gathered;
  collective_transaction(
      mpi, mpi.rank() == 0,
      "unable to allocate gathered processor names", [&] {
        gathered.resize(local_name.size() *
                        static_cast<std::size_t>(mpi.size()));
      });
  require_collective_success(
      mpi,
      MPI_Gather(local_name.data(), static_cast<int>(local_name.size()),
                 MPI_CHAR, gathered.data(),
                 static_cast<int>(local_name.size()), MPI_CHAR, 0,
                 mpi.comm()) == MPI_SUCCESS,
      "unable to gather processor names");
  PerformancePlacement result;
  collective_transaction(
      mpi, mpi.rank() == 0,
      "unable to construct performance placement metadata", [&] {
        std::vector<std::string> nodes;
        nodes.reserve(static_cast<std::size_t>(mpi.size()));
        std::ostringstream placement;
        for (int rank = 0; rank < mpi.size(); ++rank) {
          const char* begin =
              gathered.data() +
              static_cast<std::size_t>(rank) * local_name.size();
          const auto length =
              std::find(begin, begin + local_name.size(), '\0') - begin;
          std::string name(begin, static_cast<std::size_t>(length));
          if (name.empty())
            name = "unavailable";
          if (rank != 0)
            placement << ',';
          placement << rank << ':' << name;
          if (name != "unavailable")
            nodes.push_back(std::move(name));
        }
        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        if (nodes.empty()) {
          result.node_identity = "unavailable";
        } else {
          std::ostringstream node_set;
          for (std::size_t index = 0; index < nodes.size(); ++index) {
            if (index != 0U)
              node_set << ',';
            node_set << nodes[index];
          }
          result.node_identity = node_set.str();
        }
        result.rank_placement = placement.str();
      });
  broadcast_root_string(mpi, result.node_identity,
                        "performance node identity");
  broadcast_root_string(mpi, result.rank_placement,
                        "performance rank placement");
  return result;
}

diagnostics::ArtifactMetadata performance_metadata(
    const config::FlowCaseConfig& config,
    const runtime::StructuredDecomposition& decomposition,
    const MpiContext& mpi) {
  const auto grid = decomposition.process_grid();
  const auto local = decomposition.local_extent();
  const auto placement = performance_placement(mpi);
  std::optional<diagnostics::ArtifactMetadata> result;
  collective_transaction(
      mpi, true, "unable to construct performance metadata", [&] {
    diagnostics::CompatibilityMetadata compatibility;
    compatibility.hardware_identity = "unavailable";
    compatibility.node_identity = placement.node_identity;
    compatibility.mpi_identity = mpi_identity();
#if defined(__clang__)
    compatibility.compiler_identity = "clang";
#elif defined(__GNUC__)
    compatibility.compiler_identity = "gcc";
#else
    compatibility.compiler_identity = "unavailable";
#endif
    compatibility.compiler_version = __VERSION__;
    compatibility.compiler_flags =
        std::string(detail::performance_compiler_flags);
    compatibility.link_flags = std::string(detail::performance_link_flags);
    compatibility.build_type = std::string(detail::performance_build_type);
    compatibility.cpu_affinity = "unavailable";
    compatibility.rank_placement = placement.rank_placement;
    compatibility.problem_fingerprint = detail::tagged_performance_crc64(
        "hundun-case-crc64-v1",
        config::to_resolved_json(config::ResolvedCase{config}));
    compatibility.numerical_tolerance_contract = "stage2-default-v1";
    compatibility.measurement_method = "mpi-wtime-v1";
    compatibility.warmup_steps =
        static_cast<std::uint64_t>(config.performance.warmup_steps);
    compatibility.measured_steps =
        static_cast<std::uint64_t>(config.performance.measured_steps);
    compatibility.repetitions = config.performance.repetitions;
    compatibility.execution_backend = "cpu_reference";
    compatibility.ranks = mpi.size();
    compatibility.threads = 1;
    compatibility.process_grid = {grid.x, grid.y, grid.z};
    compatibility.global_owned_cell_extents = {
        static_cast<std::uint64_t>(config.mesh.cells.x),
        static_cast<std::uint64_t>(config.mesh.cells.y),
        static_cast<std::uint64_t>(config.mesh.cells.z)};
    compatibility.per_rank_owned_cell_extents = {
        static_cast<std::uint64_t>(local.x),
        static_cast<std::uint64_t>(local.y),
        static_cast<std::uint64_t>(local.z)};
    result = diagnostics::ArtifactMetadata{
        std::string(detail::performance_source_commit),
        detail::performance_source_clean,
        std::string(detail::performance_source_dirty_summary),
        std::move(compatibility)};
  });
  return std::move(*result);
}

bool same_work(const diagnostics::PerformanceWorkRecord& left,
               const diagnostics::PerformanceWorkRecord& right) {
  return left.phase == right.phase &&
         left.relative_step == right.relative_step &&
         left.slot == right.slot &&
         left.termination == right.termination &&
         left.iterations == right.iterations && left.matvec == right.matvec &&
         left.preconditioner == right.preconditioner &&
         left.reduction == right.reduction &&
         left.initial_residual_bits == right.initial_residual_bits &&
         left.recursive_residual_bits == right.recursive_residual_bits &&
         left.independent_final_residual_bits ==
             right.independent_final_residual_bits;
}

int run_performance_case(
    MpiContext& mpi, const config::FlowCaseConfig& config,
    const std::filesystem::path& authoritative_case_root) {
  const auto expected_steps = checked_add(
      static_cast<std::uint64_t>(config.performance.warmup_steps),
      static_cast<std::uint64_t>(config.performance.measured_steps),
      "performance step count");
  require_collective_success(
      mpi, !config.restart.read,
      "performance mode does not accept restart input");
  require_collective_success(
      mpi,
      expected_steps == static_cast<std::uint64_t>(config.time.steps),
      "performance time.steps must equal warmup_steps + measured_steps");
  require_collective_success(
      mpi,
      config.density_model == DensityModel::constant &&
          config.mesh.mapping == config::MeshMapping::uniform_box,
      "performance mode requires the constant-density uniform-box path");
  require_collective_success(
      mpi, !config.diagnostics.write_mesh,
      "performance mode excludes MeshDiag from the measurement fixture");
  const auto entry_allocation = execution::allocation_counters();
  require_collective_success(
      mpi, entry_allocation.allocation_events == 0U &&
               entry_allocation.allocated_bytes == 0U &&
               entry_allocation.deallocation_events == 0U &&
               entry_allocation.deallocated_bytes == 0U &&
               entry_allocation.live_bytes == 0U &&
               entry_allocation.peak_live_bytes == 0U,
      "performance mode requires zero entry allocation counters");

  diagnostics::PerformanceCorrectnessRecord correctness;
  correctness.passed = true;
  correctness.repetitions =
      static_cast<std::uint64_t>(config.performance.repetitions);
  std::vector<diagnostics::RawSample> raw_samples;
  std::vector<HaloEvidenceRecord> halo_evidence;
  std::vector<diagnostics::PerformanceWorkRecord> reference_work;
  std::string reference_state;
  PerformanceTotals totals;
  std::optional<diagnostics::ArtifactMetadata> metadata;
  std::optional<execution::AllocationCounters> reference_local_assembly;

  for (int repetition = 0; repetition < config.performance.repetitions;
       ++repetition) {
    const auto repetition_entry = execution::allocation_counters();
    require_collective_success(
        mpi, repetition_entry.live_bytes == 0U,
        "performance repetition began with live Buffers");
    std::vector<diagnostics::PerformanceWorkRecord> local_work;
    double elapsed = 0.0;
    {
      auto decomposition = runtime::StructuredDecomposition::create(
          mpi, config.mesh.cells, periodicity(config),
          runtime::DecompositionOptions{config.resources.process_grid});
      mesh::MeshTopology topology(decomposition);
      mesh::MeshGeometry geometry(
          topology,
          mesh::UniformBoxMapping(config.mesh.origin_m, config.mesh.length_m));
      auto boundaries = boundary::BoundaryRegistry::create(config, topology);
      auto fields = declare_fields(config);
      auto state = flow::FlowState::create(
          fields.registry,
          {local_extent(decomposition), topology.local_face_count()},
          fields.ids,
          {0U, 0.0, config.time.initial_dt_s, 0.0,
           flow::MomentumTimeOrder::backward_euler});
      const auto seed = initial_values(config, topology, geometry, boundaries);
      state.seed_accepted_layers(seed, seed);
      auto exchange_plan = runtime::ExchangePlan::create(
          decomposition, local_extent(decomposition), 2);
      auto halo = runtime::HaloExchange::create(decomposition, exchange_plan);
      execution::CpuReferenceContext execution;
      linear::BiCGStabSolver momentum_solver(execution, mpi);
      linear::ConjugateGradientSolver pressure_solver(execution, mpi);
      linear::JacobiPreconditioner momentum_x(execution),
          momentum_y(execution), momentum_z(execution),
          pressure_preconditioner(execution);
      auto controller = flow::Bdf2RetryController::create(
          config.time, config.density_model, topology, geometry, mpi, state);
      std::filesystem::path diagnostics_directory;
      collective_transaction(
          mpi, true, "unable to derive performance diagnostics path", [&] {
            if (mpi.rank() == 1 &&
                performance_failure_injection(
                    "nonroot_diagnostic_path"))
              throw Error("injected diagnostics path failure");
            diagnostics_directory =
                beneath(authoritative_case_root,
                        config.diagnostics.directory);
          });
      std::optional<diagnostics::DiagnosticSession> diagnostic_session;
      if (has_scheduled_diagnostics(
              state.metadata().step, expected_steps,
              static_cast<std::uint64_t>(
                  config.diagnostics.write_interval),
              false))
        diagnostic_session.emplace(
            diagnostics_directory, config.diagnostics.write_interval,
            mpi.rank());
      std::filesystem::path checkpoint_root;
      collective_transaction(
          mpi, true, "unable to derive performance checkpoint path", [&] {
            if (mpi.rank() == 1 &&
                performance_failure_injection(
                    "nonroot_checkpoint_path"))
              throw Error("injected checkpoint path failure");
            checkpoint_root =
                beneath(authoritative_case_root,
                        config.restart.write_directory);
          });
      ensure_directory(mpi, checkpoint_root);
      auto facade = flow::FixedStepConstantDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&momentum_x, &momentum_y, &momentum_z},
          pressure_solver, pressure_preconditioner,
          constant_transport_specs(fields, config));
      if (!metadata) {
        auto candidate_metadata =
            performance_metadata(config, decomposition, mpi);
        collective_transaction(
            mpi, true, "unable to retain performance metadata", [&] {
              metadata = std::move(candidate_metadata);
            });
      }

      const auto assembly = execution::allocation_counters();
      const bool assembly_monotonic =
          assembly.allocation_events >= repetition_entry.allocation_events &&
          assembly.allocated_bytes >= repetition_entry.allocated_bytes &&
          assembly.deallocation_events >=
              repetition_entry.deallocation_events &&
          assembly.deallocated_bytes >=
              repetition_entry.deallocated_bytes;
      require_collective_success(
          mpi, assembly_monotonic,
          "performance assembly allocation counters regressed");
      execution::AllocationCounters assembly_delta;
      if (assembly_monotonic) {
        assembly_delta.allocation_events =
            assembly.allocation_events -
            repetition_entry.allocation_events;
        assembly_delta.allocated_bytes =
            assembly.allocated_bytes - repetition_entry.allocated_bytes;
        assembly_delta.deallocation_events =
            assembly.deallocation_events -
            repetition_entry.deallocation_events;
        assembly_delta.deallocated_bytes =
            assembly.deallocated_bytes -
            repetition_entry.deallocated_bytes;
        assembly_delta.live_bytes = assembly.live_bytes;
      }
      const bool assembly_footprint_exact =
          assembly_delta.allocation_events >=
              assembly_delta.deallocation_events &&
          assembly_delta.allocated_bytes >=
              assembly_delta.deallocated_bytes &&
          assembly_delta.allocated_bytes -
                  assembly_delta.deallocated_bytes ==
              assembly_delta.live_bytes;
      require_collective_success(
          mpi, assembly_footprint_exact,
          "performance assembly allocation footprint is not exact");
      if (!reference_local_assembly) {
        reference_local_assembly = assembly_delta;
      } else {
        require_collective_success(
            mpi,
            assembly_delta.allocation_events ==
                    reference_local_assembly->allocation_events &&
                assembly_delta.allocated_bytes ==
                    reference_local_assembly->allocated_bytes &&
                assembly_delta.deallocation_events ==
                    reference_local_assembly->deallocation_events &&
                assembly_delta.deallocated_bytes ==
                    reference_local_assembly->deallocated_bytes &&
                assembly_delta.live_bytes ==
                    reference_local_assembly->live_bytes,
            "performance assembly allocation footprint differs between repetitions");
      }
      const auto global_assembly =
          sum_rank_u64(mpi, assembly.live_bytes,
                       "performance assembly allocation footprint");
      collective_transaction(
          mpi, mpi.rank() == 0,
          "performance assembly allocation footprint", [&] {
            totals.assembly_live =
                checked_add(totals.assembly_live, global_assembly,
                            "performance assembly allocation footprint");
          });

      for (int step = 0; step < config.performance.warmup_steps; ++step) {
        auto advance = advance_one(config, state, facade, controller);
        require_collective_success(
            mpi, advance.attempt_count() == 1U,
            "performance warmup step retried");
        collective_transaction(
            mpi, true, "unable to retain performance warmup work", [&] {
              const auto work = performance_work_records(
                  advance, static_cast<std::uint64_t>(repetition), 'W',
                  static_cast<std::uint64_t>(step));
              local_work.insert(local_work.end(), work.begin(), work.end());
            });
      }

      const auto runtime_halo_before = halo.performance_counters();
      const auto pressure_halo_before =
          facade.pressure_halo_performance_counters();
      const auto fp64_before = mpi.fp64_reduction_counters();
      std::uint64_t checkpoint_collectives = 0U;
      std::uint64_t checkpoint_logical_bytes = 0U;
      std::uint64_t checkpoint_actual_bytes = 0U;
      std::uint64_t diagnostic_logical_bytes = 0U;
      double checkpoint_io_seconds = 0.0;
      double diagnostic_io_seconds = 0.0;

      mpi.barrier();
      for (int step = 0; step < config.performance.measured_steps; ++step) {
        std::optional<flow::TimeAdvanceReport> advance;
        elapsed += detail::measure_performance_region(
            [] { return MPI_Wtime(); }, [&] {
              advance.emplace(
                  advance_one(config, state, facade, controller));
            });
        require_collective_success(
            mpi, advance->attempt_count() == 1U,
            "performance measured step retried");
        collective_transaction(
            mpi, true, "unable to retain performance measured work", [&] {
              const auto work = performance_work_records(
                  *advance, static_cast<std::uint64_t>(repetition), 'M',
                  static_cast<std::uint64_t>(step));
              local_work.insert(local_work.end(), work.begin(), work.end());
            });

        std::optional<flow::CheckpointV2Report> checkpoint_report;
        if (state.metadata().step %
                static_cast<std::uint64_t>(
                    config.restart.write_interval) ==
            0U) {
          const double checkpoint_elapsed =
              detail::measure_performance_region(
                  [] { return MPI_Wtime(); }, [&] {
                    checkpoint_report = flow::write_checkpoint_v2(
                        mpi, decomposition, topology, geometry, boundaries,
                        config, state, controller.state(), std::nullopt,
                        checkpoint_step_path(checkpoint_root,
                                             state.metadata().step));
                  });
          elapsed += checkpoint_elapsed;
          checkpoint_io_seconds += checkpoint_elapsed;
          require_collective_success(
              mpi,
              checkpoint_report->disposition() ==
                  flow::CheckpointV2Disposition::completed,
              "performance checkpoint write failed");
          collective_transaction(
              mpi, true, "unable to retain performance checkpoint counters",
              [&] {
                checkpoint_collectives = checked_add(
                    checkpoint_collectives,
                    checkpoint_report->collective_count(),
                    "performance checkpoint collectives");
                checkpoint_logical_bytes = checked_add(
                    checkpoint_logical_bytes,
                    checkpoint_report->local_logical_bytes(),
                    "performance checkpoint logical bytes");
                checkpoint_actual_bytes = checked_add(
                    checkpoint_actual_bytes,
                    checkpoint_report->global_actual_bytes(),
                    "performance checkpoint actual bytes");
              });
        }
        const bool diagnostic_due =
            diagnostic_session.has_value() &&
            diagnostic_session->due(state.metadata().step);
        std::uint64_t step_diagnostic_bytes = 0U;
        const auto execute_diagnostics = [&] {
          step_diagnostic_bytes = execute_step_diagnostics(
              config, mpi, decomposition, topology, geometry, boundaries,
              fields, exchange_plan, execution, state, facade, controller,
              *advance, checkpoint_report,
              diagnostic_session ? &*diagnostic_session : nullptr);
        };
        if (diagnostic_due) {
          const double diagnostic_elapsed =
              detail::measure_performance_region(
                  [] { return MPI_Wtime(); }, execute_diagnostics);
          elapsed += diagnostic_elapsed;
          diagnostic_io_seconds += diagnostic_elapsed;
        } else {
          execute_diagnostics();
        }
        collective_transaction(
            mpi, true, "unable to retain performance diagnostic counters",
            [&] {
              diagnostic_logical_bytes = checked_add(
                  diagnostic_logical_bytes, step_diagnostic_bytes,
                  "performance diagnostic logical bytes");
            });
      }
      mpi.barrier();
      require_collective_success(
          mpi, std::isfinite(elapsed) && elapsed > 0.0,
          "performance elapsed time must be finite and positive");

      const auto runtime_halo_after = halo.performance_counters();
      const auto pressure_halo_after =
          facade.pressure_halo_performance_counters();
      const auto fp64_after = mpi.fp64_reduction_counters();
      runtime::HaloPerformanceCounters runtime_delta;
      runtime::HaloPerformanceCounters pressure_delta;
      runtime::Fp64ReductionCounters reduction_delta;
      collective_transaction(
          mpi, true, "unable to compute performance counter deltas", [&] {
            runtime_delta =
                halo_delta(runtime_halo_before, runtime_halo_after);
            pressure_delta =
                halo_delta(pressure_halo_before, pressure_halo_after);
            reduction_delta = fp64_delta(fp64_before, fp64_after);
          });

      validate_halo_counter_consistency(mpi, runtime_delta, "runtime");
      validate_halo_counter_consistency(mpi, pressure_delta, "pressure");
      gather_halo_evidence(
          mpi, static_cast<std::uint64_t>(repetition), runtime_delta,
          pressure_delta, halo_evidence);
      if (mpi.rank() == 1 &&
          performance_failure_injection("nonroot_work_key") &&
          !local_work.empty())
        ++local_work.front().relative_step;
      if (mpi.rank() == 1 &&
          performance_failure_injection("nonroot_work_count") &&
          !local_work.empty())
        local_work.pop_back();
      validate_collective_work(mpi, local_work);
      const std::string state_fingerprint =
          collective_state_fingerprint(
              mpi, topology, state, controller, std::nullopt);
      bool repetition_matches = true;
      if (repetition == 0) {
        collective_transaction(
            mpi, true, "unable to retain reference performance evidence",
            [&] {
              reference_work = local_work;
              reference_state = state_fingerprint;
            });
      } else {
        collective_transaction(
            mpi, true, "unable to compare performance repetitions", [&] {
              repetition_matches =
                  reference_state == state_fingerprint &&
                  reference_work.size() == local_work.size();
              if (repetition_matches) {
                for (std::size_t index = 0; index < local_work.size();
                     ++index)
                  repetition_matches =
                      repetition_matches &&
                      same_work(reference_work[index], local_work[index]);
              }
            });
      }
      require_collective_success(
          mpi, repetition_matches,
          "performance repetitions produced different work or state");
      collective_transaction(
          mpi, mpi.rank() == 0,
          "unable to retain performance correctness evidence", [&] {
            correctness.states.emplace_back(
                static_cast<std::uint64_t>(repetition), state_fingerprint);
            correctness.work.insert(correctness.work.end(), local_work.begin(),
                                    local_work.end());
          });

      std::vector<double> gathered_elapsed;
      collective_transaction(
          mpi, mpi.rank() == 0,
          "unable to allocate gathered performance elapsed samples", [&] {
            gathered_elapsed.resize(static_cast<std::size_t>(mpi.size()));
          });
      require_collective_success(
          mpi,
          MPI_Gather(&elapsed, 1, MPI_DOUBLE, gathered_elapsed.data(), 1,
                     MPI_DOUBLE, 0, mpi.comm()) == MPI_SUCCESS,
          "unable to gather performance elapsed samples");
      collective_transaction(
          mpi, mpi.rank() == 0,
          "unable to retain gathered performance elapsed samples", [&] {
            for (int rank = 0; rank < mpi.size(); ++rank)
              raw_samples.push_back(
                  {repetition, rank,
                   gathered_elapsed[static_cast<std::size_t>(rank)],
                   static_cast<std::uint64_t>(
                       config.performance.measured_steps)});
          });

      const auto add_rank_halo = [&](runtime::HaloPerformanceCounters value) {
        runtime::HaloPerformanceCounters global;
        global.completed_exchanges =
            sum_rank_u64(mpi, value.completed_exchanges,
                         "global Halo completed exchanges");
        global.begin_calls =
            sum_rank_u64(mpi, value.begin_calls, "global Halo begin calls");
        global.wait_calls =
            sum_rank_u64(mpi, value.wait_calls, "global Halo wait calls");
        global.send_payload_bytes =
            sum_rank_u64(mpi, value.send_payload_bytes,
                         "global Halo send payload");
        global.receive_payload_bytes =
            sum_rank_u64(mpi, value.receive_payload_bytes,
                         "global Halo receive payload");
        global.pack_bytes =
            sum_rank_u64(mpi, value.pack_bytes, "global Halo pack bytes");
        global.unpack_bytes =
            sum_rank_u64(mpi, value.unpack_bytes, "global Halo unpack bytes");
        global.send_messages =
            sum_rank_u64(mpi, value.send_messages,
                         "global Halo send messages");
        global.receive_messages =
            sum_rank_u64(mpi, value.receive_messages,
                         "global Halo receive messages");
        collective_transaction(
            mpi, mpi.rank() == 0,
            "unable to accumulate global Halo counters",
            [&] { add_halo_totals(totals, global); });
      };
      add_rank_halo(runtime_delta);
      add_rank_halo(pressure_delta);

      const auto common_fp64_calls = common_rank_u64(
          mpi, reduction_delta.collective_calls,
          "FP64 reduction calls differ between ranks");
      const auto common_fp64_scalars = common_rank_u64(
          mpi, reduction_delta.reduced_scalars,
          "FP64 reduction scalar counts differ between ranks");
      const auto common_fp64_bytes = common_rank_u64(
          mpi, reduction_delta.logical_payload_bytes,
          "FP64 reduction logical bytes differ between ranks");
      const auto common_checkpoint = common_rank_u64(
          mpi, checkpoint_collectives,
          "checkpoint collective counts differ between ranks");
      const auto common_checkpoint_actual = common_rank_u64(
          mpi, checkpoint_actual_bytes,
          "checkpoint actual bytes differ between ranks");
      const auto global_checkpoint_bytes = sum_rank_u64(
          mpi, checkpoint_logical_bytes,
          "checkpoint logical bytes across ranks");
      const auto global_diagnostic_bytes = sum_rank_u64(
          mpi, diagnostic_logical_bytes,
          "diagnostic logical bytes across ranks");
      if (mpi.rank() == 1 &&
          performance_failure_injection("nonroot_nonfinite_io"))
        checkpoint_io_seconds =
            std::numeric_limits<double>::quiet_NaN();
      if (mpi.rank() == 1 &&
          performance_failure_injection("nonroot_negative_io"))
        diagnostic_io_seconds = -1.0;
      const auto checkpoint_duration = maximum_rank_seconds(
          mpi, checkpoint_io_seconds,
          "unable to gather checkpoint I/O duration");
      const auto diagnostic_duration = maximum_rank_seconds(
          mpi, diagnostic_io_seconds,
          "unable to gather diagnostic I/O duration");
      collective_transaction(
          mpi, mpi.rank() == 0,
          "unable to accumulate performance totals", [&] {
            totals.fp64.collective_calls =
                checked_add(totals.fp64.collective_calls, common_fp64_calls,
                            "FP64 collective calls");
            totals.fp64.reduced_scalars =
                checked_add(totals.fp64.reduced_scalars, common_fp64_scalars,
                            "FP64 reduced scalars");
            totals.fp64.logical_payload_bytes = checked_add(
                totals.fp64.logical_payload_bytes, common_fp64_bytes,
                "FP64 logical payload");
            totals.checkpoint_collectives = checked_add(
                totals.checkpoint_collectives, common_checkpoint,
                "checkpoint collective count");
            totals.checkpoint_logical_bytes = checked_add(
                totals.checkpoint_logical_bytes, global_checkpoint_bytes,
                "checkpoint logical bytes");
            totals.checkpoint_actual_bytes = checked_add(
                totals.checkpoint_actual_bytes, common_checkpoint_actual,
                "checkpoint actual bytes");
            totals.diagnostic_logical_bytes = checked_add(
                totals.diagnostic_logical_bytes, global_diagnostic_bytes,
                "diagnostic logical bytes");
            totals.diagnostic_actual_bytes = checked_add(
                totals.diagnostic_actual_bytes, global_diagnostic_bytes,
                "diagnostic actual bytes");
            totals.checkpoint_io_seconds += checkpoint_duration;
            totals.diagnostic_io_seconds += diagnostic_duration;
            if (!std::isfinite(totals.checkpoint_io_seconds) ||
                !std::isfinite(totals.diagnostic_io_seconds))
              throw Error("performance I/O duration is invalid");
          });
    }
    const auto retired = execution::allocation_counters();
    const bool fully_retired =
        retired.allocation_events >= repetition_entry.allocation_events &&
        retired.allocated_bytes >= repetition_entry.allocated_bytes &&
        retired.deallocation_events >= repetition_entry.deallocation_events &&
        retired.deallocated_bytes >= repetition_entry.deallocated_bytes &&
        retired.allocation_events - repetition_entry.allocation_events ==
            retired.deallocation_events -
                repetition_entry.deallocation_events &&
        retired.allocated_bytes - repetition_entry.allocated_bytes ==
            retired.deallocated_bytes - repetition_entry.deallocated_bytes &&
        retired.live_bytes == 0U;
    require_collective_success(
        mpi, fully_retired,
        "performance repetition did not exactly retire every Buffer");
  }

  const auto final_allocation = execution::allocation_counters();
  const auto global_peak =
      sum_rank_u64(mpi, final_allocation.peak_live_bytes,
                   "performance run-level peak allocation");
  std::string json;
  collective_transaction(
      mpi, mpi.rank() == 0,
      "unable to serialize performance artifact", [&] {
    if (performance_failure_injection("root_overflow"))
      static_cast<void>(checked_add(
          std::numeric_limits<std::uint64_t>::max(), 1U,
          "injected root performance postprocessing"));
    if (performance_failure_injection("missing_work") &&
        !correctness.work.empty())
      correctness.work.pop_back();
    if (performance_failure_injection("malformed_work") &&
        !correctness.work.empty())
      correctness.work.front().phase = 'X';
    const auto cells_xy = checked_multiply(
        static_cast<std::uint64_t>(config.mesh.cells.x),
        static_cast<std::uint64_t>(config.mesh.cells.y),
        "global owned-cell count");
    const auto global_cells = checked_multiply(
        cells_xy, static_cast<std::uint64_t>(config.mesh.cells.z),
        "global owned-cell count");
    const auto repetitions =
        static_cast<std::uint64_t>(config.performance.repetitions);
    correctness.allocation_bytes_per_owned_cell =
        static_cast<double>(totals.assembly_live) /
        static_cast<double>(checked_multiply(
            repetitions, global_cells,
            "allocation normalization denominator"));
    correctness.peak_allocation_bytes_per_owned_cell =
        static_cast<double>(global_peak) / static_cast<double>(global_cells);
    if (!std::isfinite(correctness.allocation_bytes_per_owned_cell) ||
        !std::isfinite(
            correctness.peak_allocation_bytes_per_owned_cell))
      throw Error("performance allocation ratio is non-finite");

    std::uint64_t linear_reductions = 0U;
    std::uint64_t momentum_matvec = 0U;
    std::uint64_t pressure_matvec = 0U;
    std::uint64_t momentum_preconditioner = 0U;
    std::uint64_t pressure_preconditioner = 0U;
    for (const auto& work : correctness.work) {
      if (work.phase != 'M')
        continue;
      linear_reductions =
          checked_add(linear_reductions, work.reduction,
                      "linear reduction count");
      if (work.slot < 3U) {
        momentum_matvec =
            checked_add(momentum_matvec, work.matvec, "momentum matvec");
        momentum_preconditioner = checked_add(
            momentum_preconditioner, work.preconditioner,
            "momentum preconditioner");
      } else {
        pressure_matvec =
            checked_add(pressure_matvec, work.matvec, "pressure matvec");
        pressure_preconditioner = checked_add(
            pressure_preconditioner, work.preconditioner,
            "pressure preconditioner");
      }
    }

    diagnostics::Artifact artifact;
    artifact.metadata = std::move(*metadata);
    artifact.correctness = {
        true, diagnostics::serialize_performance_correctness(correctness)};
    artifact.aggregation = diagnostics::aggregate_samples(
        std::move(raw_samples), config.performance.repetitions, mpi.size());
    artifact.comparison = diagnostics::compare_artifact_metadata(
        artifact.metadata, artifact.metadata,
        diagnostics::ComparisonMode::identical);
    artifact.counters.allocated_bytes = {
        {"execution.allocated", totals.assembly_live},
        {"execution.peak-live", global_peak}};
    artifact.counters.halo_payload_bytes = {
        {"pack", totals.halo.pack_bytes},
        {"receive", totals.halo.receive_payload_bytes},
        {"send", totals.halo.send_payload_bytes},
        {"unpack", totals.halo.unpack_bytes}};
    artifact.counters.halo_messages = {
        {"receive", totals.halo.receive_messages},
        {"send", totals.halo.send_messages}};
    artifact.counters.collectives = {
        {"checkpoint", totals.checkpoint_collectives},
        {"fp64-reduction", totals.fp64.collective_calls},
        {"linear-reduction", linear_reductions}};
    artifact.counters.collective_logical_payload_bytes = {
        {"fp64-reduction", totals.fp64.logical_payload_bytes}};
    artifact.counters.matvec = {{"momentum", momentum_matvec},
                                {"pressure", pressure_matvec}};
    artifact.counters.preconditioner_applications = {
        {"momentum", momentum_preconditioner},
        {"pressure", pressure_preconditioner}};
    artifact.counters.logical_io_bytes = {
        {"checkpoint", totals.checkpoint_logical_bytes},
        {"diagnostics", totals.diagnostic_logical_bytes}};
    if (performance_failure_injection("missing_counter_key"))
      artifact.counters.allocated_bytes.erase("execution.allocated");
    if (performance_failure_injection("extra_counter_key"))
      artifact.counters.allocated_bytes.emplace("unapproved", 0U);
    validate_canonical_counter_maps(artifact.counters);
    json = diagnostics::to_json(artifact);
  });
  std::string evidence_input;
  collective_transaction(
      mpi, mpi.rank() == 0,
      "unable to serialize performance evidence input", [&] {
        evidence_input = serialize_performance_evidence_input(
            halo_evidence, totals, config.performance.repetitions,
            mpi.size());
      });
  atomic_publish_performance(
      mpi, authoritative_case_root, config.performance.directory, json,
      evidence_input);
  mpi.barrier();
  return EXIT_SUCCESS;
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
  if (config.performance.enabled)
    return run_performance_case(mpi, config, authoritative_case_root);

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
    try {
      mesh_record = diagnostics::make_mesh_diagnostic_v2(
          mpi.rank(), mpi.size(), topology, geometry);
    } catch (const std::exception&) {
      mesh_ok = false;
    } catch (...) {
      mesh_ok = false;
    }
    require_collective_success(
        mpi, mesh_ok, "unable to construct mesh diagnostic");
    diagnostics::write_mesh_diagnostic_v2(
        mpi, diagnostics_directory, *mesh_record);
  }

  print_start(config, mpi);
  if (restored && diagnostic_session &&
      diagnostic_session->due(state.metadata().step)) {
    DiagnosticBatch batch;
    bool batch_ok = true;
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
    } catch (const std::exception&) {
      batch_ok = false;
    } catch (...) {
      batch_ok = false;
    }
    require_collective_success(
        mpi, batch_ok, "unable to construct diagnostic batch");
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
             halo, checkpoint_root, diagnostics_directory);
  } else if (config.density_model == DensityModel::material) {
    auto facade = flow::FixedStepMaterialDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&momentum_x, &momentum_y, &momentum_z},
        pressure_solver, pressure_preconditioner, fields.registry, fields.ids,
        material_transport_spec(fields, config));
    run_loop(config, mpi, decomposition, topology, geometry, boundaries,
             fields, exchange_plan, execution, state, facade, controller,
             diagnostic_session ? &*diagnostic_session : nullptr,
             halo, checkpoint_root, diagnostics_directory);
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
             halo, checkpoint_root, diagnostics_directory);
  }
  print_finished(state, mpi);
  mpi.barrier();
  return EXIT_SUCCESS;
}

}  // namespace hundun::application
