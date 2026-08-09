// SPDX-License-Identifier: Apache-2.0

#include "app_immersed_flow_driver_detail.hpp"

#include "flow_density_closure_detail.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/diag_immersed_module.hpp"
#include "hundun/diag_immersed_static.hpp"
#include "hundun/diag_les_wale.hpp"
#include "hundun/diag_session.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/flow_checkpoint_v3.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_surface_query.hpp"
#include "hundun/lin_bicgstab.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/les_wale.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_structured_decomposition.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
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
#include <vector>

namespace hundun::application {

bool stage3_dispatch_inventory_complete(
    const flow::CheckpointV3WriteModules& modules) noexcept {
  bool expects_ibm = false;
  bool expects_wale = false;
  bool expects_ideal_gas = false;
  switch (modules.presence) {
  case flow::CheckpointV3Presence::constant_static_ibm:
    expects_ibm = true;
    break;
  case flow::CheckpointV3Presence::constant_body_fitted_wale:
    expects_wale = true;
    break;
  case flow::CheckpointV3Presence::constant_static_ibm_wale:
    expects_ibm = true;
    expects_wale = true;
    break;
  case flow::CheckpointV3Presence::material_static_ibm:
    expects_ibm = true;
    break;
  case flow::CheckpointV3Presence::material_body_fitted_wale:
    expects_wale = true;
    break;
  case flow::CheckpointV3Presence::material_static_ibm_wale:
    expects_ibm = true;
    expects_wale = true;
    break;
  case flow::CheckpointV3Presence::ideal_gas_static_ibm:
    expects_ibm = true;
    expects_ideal_gas = true;
    break;
  case flow::CheckpointV3Presence::ideal_gas_body_fitted_wale:
    expects_wale = true;
    expects_ideal_gas = true;
    break;
  case flow::CheckpointV3Presence::ideal_gas_static_ibm_wale:
    expects_ibm = true;
    expects_wale = true;
    expects_ideal_gas = true;
    break;
  default:
    return false;
  }

  const bool has_all_ibm =
      modules.surface != nullptr && modules.query != nullptr &&
      modules.domain != nullptr && modules.ghost_plan != nullptr &&
      modules.wall_plan != nullptr && modules.transform != nullptr &&
      modules.flow != nullptr;
  const bool has_no_ibm =
      modules.surface == nullptr && modules.query == nullptr &&
      modules.domain == nullptr && modules.ghost_plan == nullptr &&
      modules.wall_plan == nullptr && modules.transform == nullptr &&
      modules.flow == nullptr;
  return (expects_ibm ? has_all_ibm : has_no_ibm) &&
         (modules.wale != nullptr) == expects_wale &&
         (modules.ideal_gas != nullptr) == expects_ideal_gas;
}

namespace {

using runtime::Error;

constexpr std::string_view kFieldOwner = "immersed-flow-driver";

std::string_view density_model_name(config::DensityModel model) {
  switch (model) {
  case config::DensityModel::constant:
    return "constant";
  case config::DensityModel::material:
    return "material";
  case config::DensityModel::ideal_gas:
    return "ideal_gas";
  }
  throw Error("invalid immersed-flow density model");
}

runtime::FieldDescriptor cell_field(std::string name, std::string unit,
                                    std::uint32_t components, int ghost_width,
                                    bool conservative = false) {
  return {std::move(name),
          std::move(unit),
          std::string(kFieldOwner),
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost_width,
          conservative,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_field(std::string name, std::string unit,
                                    std::uint32_t components) {
  return {std::move(name),
          std::move(unit),
          std::string(kFieldOwner),
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

std::array<bool, 3> periodicity(const config::FlowCaseConfig& config) noexcept {
  std::array<bool, 3> result{};
  for (std::size_t axis = 0; axis < result.size(); ++axis) {
    result[axis] =
        config.boundaries[2U * axis].type == config::BoundaryType::periodic;
  }
  return result;
}

runtime::Int3
local_extent(const runtime::StructuredDecomposition& decomposition) noexcept {
  const auto box = decomposition.owned_box();
  return {box.end.x - box.begin.x, box.end.y - box.begin.y,
          box.end.z - box.begin.z};
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

void ensure_directory(const runtime::MpiContext& mpi,
                      const std::filesystem::path& directory) {
  bool local_ok = true;
  if (mpi.rank() == 0) {
    try {
      std::filesystem::create_directories(directory);
      local_ok = std::filesystem::is_directory(directory);
    } catch (...) {
      local_ok = false;
    }
  }
  const auto status = runtime::collective_status(
      mpi, local_ok, "unable to create immersed-flow checkpoint root");
  if (!status.ok) {
    throw Error("unable to create immersed-flow checkpoint root");
  }
}

void require_supported(const config::ImmersedFlowCaseConfig& config) {
  const auto& flow = config.common_flow;
  const bool immersed =
      config.immersed_boundary.model ==
      config::ImmersedBoundaryModel::local_flow_pattern_ghost_cell;
  const bool wale = config.les.model == config::LesModel::wale;
  if (immersed &&
      (!config.immersed_boundary.geometry || !config.immersed_boundary.wall))
    throw Error("immersed-flow driver requires complete static LFP-GCIBM");
  if (!immersed && !wale)
    throw Error("immersed-flow driver requires IBM or WALE");
  if (wale && !config.les.wale)
    throw Error("immersed-flow driver requires complete WALE controls");
  if (flow.density_model == config::DensityModel::constant &&
      !flow.scalars.empty()) {
    throw Error("immersed-flow driver does not yet transport scalars");
  }
  if (flow.time.mode != config::TimeMode::fixed) {
    throw Error("immersed-flow Checkpoint v3 currently requires fixed time");
  }
  if (flow.performance.enabled) {
    throw Error("immersed-flow driver does not yet expose performance "
                "mode");
  }
}

flow::CheckpointV3Presence checkpoint_presence(config::DensityModel density,
                                               bool immersed, bool wale) {
  switch (density) {
  case config::DensityModel::constant:
    if (immersed)
      return wale ? flow::CheckpointV3Presence::constant_static_ibm_wale
                  : flow::CheckpointV3Presence::constant_static_ibm;
    return flow::CheckpointV3Presence::constant_body_fitted_wale;
  case config::DensityModel::material:
    if (immersed)
      return wale ? flow::CheckpointV3Presence::material_static_ibm_wale
                  : flow::CheckpointV3Presence::material_static_ibm;
    return flow::CheckpointV3Presence::material_body_fitted_wale;
  case config::DensityModel::ideal_gas:
    if (immersed)
      return wale ? flow::CheckpointV3Presence::ideal_gas_static_ibm_wale
                  : flow::CheckpointV3Presence::ideal_gas_static_ibm;
    return flow::CheckpointV3Presence::ideal_gas_body_fitted_wale;
  }
  throw Error("invalid Checkpoint v3 density profile");
}

const boundary::ResolvedInletState*
inlet_state(const boundary::BoundaryRegistry& boundaries) {
  const auto patch = boundaries.velocity_inlet_patch_id();
  if (!patch) {
    return nullptr;
  }
  const auto& state = boundaries.patch(*patch).inlet_state();
  if (!state) {
    throw Error("velocity inlet has no resolved inlet state");
  }
  return &*state;
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

struct DeclaredFields final {
  runtime::FieldRegistry registry;
  flow::FlowFieldIds ids;
};

DeclaredFields declare_fields(const config::FlowCaseConfig& config,
                              int ghost_width) {
  DeclaredFields result;
  const bool variable =
      config.density_model != config::DensityModel::constant;
  result.ids.density = result.registry.declare_field(
      cell_field("rho", "kg/m3", 1U, ghost_width, variable));
  result.ids.velocity =
      result.registry.declare_field(cell_field("u", "m/s", 3U, ghost_width));
  result.ids.mechanical_pressure =
      result.registry.declare_field(cell_field("pi", "Pa", 1U, ghost_width));
  result.ids.face_velocity =
      result.registry.declare_field(face_field("uf", "m/s", 3U));
  result.ids.face_mass_flux =
      finite_volume::declare_face_mass_flux(result.registry);
  if (variable) {
    result.ids.transported_cell_fields.push_back(
        result.registry.declare_field(
            cell_field("rho_h", "J/m3", 1U, ghost_width, true)));
    for (const auto& scalar : config.scalars)
      result.ids.transported_cell_fields.push_back(
          result.registry.declare_field(cell_field(
              "rho_" + scalar.name, "kg/m3", 1U, ghost_width, true)));
  }
  result.registry.freeze();
  return result;
}

flow::MaterialDensityTransportSpec material_transport_spec(
    const DeclaredFields& fields, const config::FlowCaseConfig& config) {
  if (fields.ids.transported_cell_fields.size() !=
      1U + config.scalars.size())
    throw Error("immersed material transport field layout is invalid");
  flow::MaterialDensityTransportSpec result;
  result.enthalpy_density = fields.ids.transported_cell_fields.front();
  result.enthalpy_diffusivity_kg_per_m_s = 0.0;
  for (std::size_t scalar = 0U; scalar < config.scalars.size(); ++scalar) {
    result.scalar_densities.push_back(
        fields.ids.transported_cell_fields[scalar + 1U]);
    result.scalar_diffusivities_kg_per_m_s.push_back(
        config.physics.rho_ref_kg_per_m3 *
        config.scalars[scalar].diffusivity_m2_per_s);
  }
  return result;
}

flow::FlowLayerValues
initial_values(const config::FlowCaseConfig& config,
               const mesh::MeshTopology& topology,
               const mesh::MeshGeometry& geometry,
               const boundary::BoundaryRegistry& boundaries,
               const immersed::ImmersedDomain* domain) {
  const auto* inlet = inlet_state(boundaries);
  const double density = inlet != nullptr ? inlet->density_kg_per_m3
                                          : config.physics.rho_ref_kg_per_m3;
  const runtime::Real3 velocity =
      inlet != nullptr ? inlet->velocity_m_per_s : runtime::Real3{};
  double enthalpy = inlet != nullptr ? inlet->enthalpy_J_per_kg : 0.0;
  if (inlet == nullptr &&
      config.density_model == config::DensityModel::ideal_gas) {
    const double temperature =
        config.physics.thermodynamic_pressure_pa.value() /
        (config.physics.gas_constant_J_per_kg_K.value() * density);
    enthalpy = config.physics.cp_J_per_kg_K.value() * temperature;
  }
  const std::vector<double> scalar_values =
      inlet != nullptr ? inlet->scalar_values
                       : std::vector<double>(config.scalars.size(), 0.0);
  if (!std::isfinite(density) || density <= 0.0) {
    throw Error("invalid immersed-flow initial density");
  }
  if (!std::isfinite(enthalpy) ||
      scalar_values.size() != config.scalars.size() ||
      std::any_of(scalar_values.begin(), scalar_values.end(),
                  [](double value) { return !std::isfinite(value); }))
    throw Error("invalid immersed-flow initial transported state");

  flow::FlowLayerValues result;
  result.density.assign(topology.owned_cell_count(), 0.0);
  result.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  result.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  const bool variable =
      config.density_model != config::DensityModel::constant;
  if (variable) {
    result.transported_cell_fields.resize(1U + config.scalars.size());
    for (auto& transported : result.transported_cell_fields)
      transported.assign(topology.owned_cell_count(), 0.0);
  }
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain != nullptr &&
        domain->region(cell) != immersed::CellRegion::fluid) {
      continue;
    }
    result.density[cell] = density;
    result.velocity[cell * 3U] = velocity.x;
    result.velocity[cell * 3U + 1U] = velocity.y;
    result.velocity[cell * 3U + 2U] = velocity.z;
    if (variable) {
      result.transported_cell_fields.front()[cell] = density * enthalpy;
      for (std::size_t scalar = 0U; scalar < config.scalars.size(); ++scalar)
        result.transported_cell_fields[scalar + 1U][cell] =
            density * scalar_values[scalar];
    }
  }

  result.face_velocity.resize(topology.local_face_count() * 3U);
  result.face_mass_flux.resize(topology.local_face_count());
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto neighbour = topology.neighbour(face);
    const bool owner_active =
        domain == nullptr ||
        domain->region(topology.owner(face)) == immersed::CellRegion::fluid;
    const bool neighbour_active =
        !neighbour.has_value() || domain == nullptr ||
        domain->region(*neighbour) == immersed::CellRegion::fluid;
    if (!owner_active || !neighbour_active) {
      result.face_velocity[face * 3U] = 0.0;
      result.face_velocity[face * 3U + 1U] = 0.0;
      result.face_velocity[face * 3U + 2U] = 0.0;
      result.face_mass_flux[face] = 0.0;
      continue;
    }
    result.face_velocity[face * 3U] = velocity.x;
    result.face_velocity[face * 3U + 1U] = velocity.y;
    result.face_velocity[face * 3U + 2U] = velocity.z;
    result.face_mass_flux[face] =
        density * dot(velocity, geometry.face_area_vector_m2(
                                    face, mesh::FaceSide::owner));
  }
  return result;
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

bool finite(const immersed::ForceComponents& value) noexcept {
  return finite(value.pressure_N) && finite(value.viscous_N) &&
         finite(value.total_N);
}

void require_committed_report(
    const flow::StepAttemptReport& base,
    const flow::ImmersedFlowStepAttemptReport& report, bool immersed,
    bool wale,
    const flow::MaterialDensityStepAttemptReport* material = nullptr,
    const flow::IdealGasStepAttemptReport* ideal_gas = nullptr) {
  if (base.pressure_corrector_count != 2U) {
    throw Error("immersed-flow committed step did not use two PISO correctors");
  }
  if (!std::isfinite(base.final_continuity_normalized_l2) ||
      !std::isfinite(base.final_pressure_residual_l2) ||
      !std::all_of(base.final_momentum_normalized_l2.begin(),
                   base.final_momentum_normalized_l2.end(),
                   [](double value) { return std::isfinite(value); })) {
    throw Error("immersed-flow committed diagnostics are non-finite");
  }
  if (immersed &&
      (!report.force || !finite(report.force->operator_force) ||
       !finite(report.force->budget_reaction) ||
       !finite(report.force->surface_traction) ||
       !finite(report.force->consistency))) {
    throw Error("immersed-flow committed force report is incomplete");
  }
  if (!immersed && report.force)
    throw Error("body-fitted WALE unexpectedly reported immersed force");
  if (wale &&
      (!report.wale || report.wale->identity.value == 0U ||
       !std::isfinite(report.wale->minimum_nu_t_m2_per_s) ||
       !std::isfinite(report.wale->maximum_nu_t_m2_per_s) ||
       !std::isfinite(report.wale->l2_nu_t_m2_per_s)))
    throw Error("immersed-flow committed WALE report is incomplete");
  if (material != nullptr &&
      (!material->material_report_available() ||
       material->material_failure_reason() !=
           flow::MaterialTransportFailureReason::none ||
       material->flux_provenance() !=
           flow::MaterialFluxProvenance::final_corrected ||
       !material->mass_conservation_available() ||
       !material->material_report().density_residual_available() ||
       !material->material_report().minimum_density_available()))
    throw Error("immersed-flow committed material report is incomplete");
  if (ideal_gas != nullptr &&
      (!ideal_gas->closure_report_available() ||
       ideal_gas->closure_report().disposition() !=
           flow::IdealGasClosureDisposition::closed ||
       ideal_gas->closure_report().reason() !=
           flow::IdealGasClosureFailureReason::none ||
       ideal_gas->closure_report().stage() !=
           flow::IdealGasClosureStage::final ||
       ideal_gas->closure_report().evaluation_count() != 3U ||
       !ideal_gas->closure_report().final_metrics_available()))
    throw Error("immersed-flow committed ideal-gas report is incomplete");
}

void print_vector(std::ostream& output, runtime::Real3 value) {
  output << value.x << ',' << value.y << ',' << value.z;
}

void print_step(const flow::FlowState& state, std::uint32_t attempts,
                const flow::StepAttemptReport& base,
                const flow::ImmersedFlowStepAttemptReport& report,
                const runtime::MpiContext& mpi) {
  if (mpi.rank() != 0) {
    return;
  }
  std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "STEP step=" << state.metadata().step
            << " time_s=" << state.metadata().time_s
            << " dt_s=" << state.metadata().dt_s << " attempts=" << attempts
            << " correctors=" << base.pressure_corrector_count
            << " continuity=" << base.final_continuity_normalized_l2
            << " pressure_residual=" << base.final_pressure_residual_l2;
  if (report.force) {
    std::cout << " force_operator=";
    print_vector(std::cout, report.force->operator_force.total_N);
    std::cout << " force_budget_reaction=";
    print_vector(std::cout, report.force->budget_reaction.total_N);
    std::cout << " force_surface_traction=";
    print_vector(std::cout, report.force->surface_traction.total_N);
    std::cout << " force_consistency=";
    print_vector(std::cout, report.force->consistency.total_N);
  }
  if (report.wale) {
    std::cout << " wale_identity=" << report.wale->identity.value
              << " wale_nu_t_min="
              << report.wale->minimum_nu_t_m2_per_s
              << " wale_nu_t_max="
              << report.wale->maximum_nu_t_m2_per_s
              << " wale_zero_count=" << report.wale->exact_zero_count;
  }
  std::cout << '\n';
}

double next_attempt_dt(const config::FlowTimeConfig& time, double current,
                       const flow::StepAttemptReport& report) {
  double candidate = current * time.retry_factor;
  if (std::isfinite(report.suggested_dt_s) && report.suggested_dt_s > 0.0) {
    candidate = std::min(candidate, report.suggested_dt_s);
  }
  return std::max(time.min_dt_s, candidate);
}

double next_accepted_dt(const config::FlowTimeConfig& time, double current,
                        const flow::StepAttemptReport& report) {
  if (time.mode == config::TimeMode::fixed) {
    return time.initial_dt_s;
  }
  double candidate = current * time.growth_factor;
  if (std::isfinite(report.suggested_dt_s) && report.suggested_dt_s > 0.0) {
    candidate = report.suggested_dt_s;
  }
  return std::max(time.min_dt_s, std::min(time.max_dt_s, candidate));
}

bool has_future_diagnostic_step(std::uint64_t current, std::uint64_t target,
                                int write_interval) noexcept {
  if (current >= target || write_interval < 1)
    return false;
  const auto interval = static_cast<std::uint64_t>(write_interval);
  return current / interval < target / interval;
}

diagnostics::DiagnosticRequest diagnostic_request(
    const runtime::MpiContext& mpi, const flow::FlowState& state,
    diagnostics::DiagnosticScope scope, std::string_view phase) {
  diagnostics::DiagnosticRequest result;
  result.level = diagnostics::DiagnosticLevel::summary;
  result.scope = scope;
  result.frame = {mpi.rank(), state.metadata().step,
                  state.metadata().time_s, phase};
  return result;
}

void append_stage3_diagnostics(
    const runtime::MpiContext& mpi, const flow::FlowState& state,
    flow::CheckpointV3Presence presence,
    const diagnostics::ImmersedStaticDiagnosticSummary* static_summary,
    const flow::ImmersedFlowDiagnosticSource* immersed_source,
    const flow::ImmersedFlowStepAttemptReport& report,
    diagnostics::DiagnosticBatch& batch) {
  using diagnostics::DiagnosticModuleKind;
  const auto inventory = diagnostics::stage3_added_provider_inventory(presence);
  for (const auto kind : inventory) {
    if (kind == DiagnosticModuleKind::immersed_surface ||
        kind == DiagnosticModuleKind::ghost_stencil ||
        kind == DiagnosticModuleKind::local_flow_pattern) {
      if (static_summary == nullptr || immersed_source == nullptr)
        throw Error("Stage 3 static diagnostic source is unavailable");
      const auto phase =
          kind == DiagnosticModuleKind::immersed_surface
              ? std::string_view("immersed-static.surface")
              : kind == DiagnosticModuleKind::ghost_stencil
                    ? std::string_view("immersed-static.ghost-stencil")
                    : std::string_view("immersed-static.local-flow-pattern");
      const auto descriptor =
          diagnostics::describe_diagnostics(*static_summary, kind);
      auto request = diagnostic_request(
          mpi, state, diagnostics::DiagnosticScope::collective, phase);
      request.frame.step = 0U;
      request.frame.time_s = 0.0;
      diagnostics::DiagnosticBatchSink sink(descriptor, request, batch);
      diagnostics::collect_diagnostics(*static_summary, kind, mpi, request,
                                       sink);
      continue;
    }
    if (kind == DiagnosticModuleKind::wall_force) {
      if (immersed_source == nullptr)
        throw Error("Stage 3 wall-force diagnostic source is unavailable");
      const auto descriptor =
          diagnostics::describe_diagnostics(*immersed_source, kind);
      auto request = diagnostic_request(
          mpi, state, diagnostics::DiagnosticScope::collective,
          "immersed-flow.wall-force");
      diagnostics::DiagnosticBatchSink sink(descriptor, request, batch);
      diagnostics::collect_diagnostics(*immersed_source, kind, mpi, request,
                                       sink);
      continue;
    }
    if (kind == DiagnosticModuleKind::les) {
      if (!report.wale)
        throw Error("Stage 3 WALE diagnostic source is unavailable");
      const auto descriptor = diagnostics::describe_diagnostics(*report.wale);
      auto request = diagnostic_request(
          mpi, state, diagnostics::DiagnosticScope::local, "les.wale");
      diagnostics::DiagnosticBatchSink sink(descriptor, request, batch);
      diagnostics::collect_diagnostics(*report.wale, request, sink);
      continue;
    }
    throw Error("Stage 3 diagnostic inventory contains an invalid kind");
  }
}

} // namespace

int run_immersed_flow_case(
    const CliOptions&, runtime::MpiContext& mpi,
    const config::ImmersedFlowCaseConfig& immersed_config,
    const std::filesystem::path& authoritative_case_root) {
  require_supported(immersed_config);
  const auto& config = immersed_config.common_flow;
  if (config.resources.expected_ranks &&
      *config.resources.expected_ranks != mpi.size()) {
    throw Error("configured MPI rank count differs from runtime");
  }

  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, config.mesh.cells, periodicity(config),
      runtime::DecompositionOptions{config.resources.process_grid});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry =
      config.mesh.mapping == config::MeshMapping::uniform_box
          ? mesh::MeshGeometry(topology,
                               mesh::UniformBoxMapping(config.mesh.origin_m,
                                                       config.mesh.length_m))
          : mesh::MeshGeometry(topology,
                               mesh::AnalyticWarpedBoxMapping(
                                   config.mesh.origin_m, config.mesh.length_m,
                                   config.mesh.warp_amplitude.value()));
  auto boundaries = boundary::BoundaryRegistry::create(config, topology);
  const bool has_immersed =
      immersed_config.immersed_boundary.model ==
      config::ImmersedBoundaryModel::local_flow_pattern_ghost_cell;
  const bool has_wale = immersed_config.les.model == config::LesModel::wale;
  std::optional<immersed::ImmersedSurface> surface;
  std::optional<immersed::SurfaceQuery> query;
  std::optional<immersed::ImmersedDomain> domain;
  std::optional<immersed::GhostStencilPlan> ghost_plan;
  std::optional<immersed::WallQuadraturePlan> wall_plan;
  std::optional<immersed::LocalFlowPatternTransform> transform;
  int halo_width = 2;
  if (has_immersed) {
    const auto& geometry_config =
        *immersed_config.immersed_boundary.geometry;
    const auto surface_path =
        (authoritative_case_root / geometry_config.file).lexically_normal();
    surface.emplace(immersed::ImmersedSurface::load_collective(
        surface_path, geometry_config.length_scale_to_m, mpi, 0));
    query.emplace(immersed::SurfaceQuery::create(*surface));
    domain.emplace(immersed::ImmersedDomain::create(
        *surface, *query, geometry_config.fluid_side, topology, geometry,
        boundaries, mpi));
    ghost_plan.emplace(immersed::GhostStencilPlan::create(
        *surface, *query, *domain, topology, geometry, decomposition, mpi));
    wall_plan.emplace(immersed::WallQuadraturePlan::create(
        *surface, *query, *domain, topology, geometry, mpi));
    transform.emplace();
    halo_width = static_cast<int>(std::max(
        ghost_plan->maximum_halo_reach(), wall_plan->maximum_halo_reach()));
  }

  auto fields = declare_fields(config, halo_width);
  auto state = flow::FlowState::create(
      fields.registry,
      {local_extent(decomposition), topology.local_face_count()}, fields.ids,
      {0U, 0.0, config.time.initial_dt_s, 0.0,
       flow::MomentumTimeOrder::backward_euler});
  const auto seed =
      initial_values(config, topology, geometry, boundaries,
                     domain ? &*domain : nullptr);
  state.seed_accepted_layers(seed, seed);

  std::optional<flow::IdealGasClosure> ideal_gas_closure;
  if (config.density_model == config::DensityModel::ideal_gas) {
    const flow::IdealGasClosureSpec closure_spec{
        fields.ids.transported_cell_fields.front(),
        config.physics.cp_J_per_kg_K.value(),
        config.physics.gas_constant_J_per_kg_K.value(),
        config.physics.thermodynamic_pressure_pa.value()};
    if (domain.has_value()) {
      ideal_gas_closure.emplace(
          flow::detail::DensityClosureAdapter::create_immersed(
              topology, geometry, boundaries, *domain, mpi, fields.registry,
              fields.ids, state, closure_spec));
    } else {
      ideal_gas_closure.emplace(flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, fields.registry, fields.ids,
          state, closure_spec));
    }
  }

  auto exchange_plan = runtime::ExchangePlan::create(
      decomposition, local_extent(decomposition), halo_width);
  auto halo = runtime::HaloExchange::create(decomposition, exchange_plan);
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner momentum_x(execution), momentum_y(execution),
      momentum_z(execution), pressure_preconditioner(execution);
  std::optional<les::WaleModel> wale;
  if (has_wale) {
    std::vector<mesh::GlobalCellId> active;
    std::size_t owned_active_count = topology.owned_cell_count();
    if (domain.has_value()) {
      owned_active_count = domain->active_cells().owned_active_count();
      active = domain->active_cells().ordered_global_ids();
    } else {
      active.reserve(topology.local_cell_count());
      for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count();
           ++cell)
        active.push_back(topology.global_cell_id(cell));
    }
    const auto& control = *immersed_config.les.wale;
    wale.emplace(les::WaleModel::create(
        {control.coefficient, control.turbulent_prandtl,
         control.turbulent_schmidt},
        topology, geometry, owned_active_count, active, execution));
  }
  auto facade = [&] {
    if (config.density_model != config::DensityModel::constant) {
      const flow::ImmersedFlowDensitySetup density_setup{
          config.density_model, &fields.registry, fields.ids,
          material_transport_spec(fields, config),
          ideal_gas_closure ? &*ideal_gas_closure : nullptr};
      return flow::FixedStepImmersedFlow::create(
          decomposition, topology, geometry, boundaries,
          domain ? &*domain : nullptr, ghost_plan ? &*ghost_plan : nullptr,
          wall_plan ? &*wall_plan : nullptr,
          transform ? &*transform : nullptr, wale ? &*wale : nullptr,
          density_setup, mpi, execution, halo, momentum_solver,
          {&momentum_x, &momentum_y, &momentum_z}, pressure_solver,
          pressure_preconditioner);
    }
    return flow::FixedStepImmersedFlow::create(
        decomposition, topology, geometry, boundaries,
        domain ? &*domain : nullptr, ghost_plan ? &*ghost_plan : nullptr,
        wall_plan ? &*wall_plan : nullptr, transform ? &*transform : nullptr,
        wale ? &*wale : nullptr, mpi, execution, halo, momentum_solver,
        {&momentum_x, &momentum_y, &momentum_z}, pressure_solver,
        pressure_preconditioner);
  }();
  const bool ideal_gas =
      config.density_model == config::DensityModel::ideal_gas;
  const flow::ImmersedFlowPhysics physics{
      config.density_model,
      config.physics.rho_ref_kg_per_m3,
      config.physics.dynamic_viscosity_pa_s,
      ideal_gas ? config.physics.cp_J_per_kg_K : std::nullopt,
      ideal_gas ? config.physics.gas_constant_J_per_kg_K : std::nullopt,
      ideal_gas ? config.physics.thermodynamic_pressure_pa : std::nullopt};

  const auto presence =
      checkpoint_presence(config.density_model, has_immersed, has_wale);
  flow::CheckpointV3WriteModules write_checkpoint_modules;
  write_checkpoint_modules.presence = presence;
  write_checkpoint_modules.wale = wale ? &*wale : nullptr;
  write_checkpoint_modules.ideal_gas =
      ideal_gas_closure ? &*ideal_gas_closure : nullptr;
  flow::CheckpointV3ReadModules read_checkpoint_modules;
  read_checkpoint_modules.presence = presence;
  read_checkpoint_modules.wale = wale ? &*wale : nullptr;
  read_checkpoint_modules.ideal_gas =
      ideal_gas_closure ? &*ideal_gas_closure : nullptr;
  if (has_immersed) {
    write_checkpoint_modules.surface = &*surface;
    write_checkpoint_modules.query = &*query;
    write_checkpoint_modules.domain = &*domain;
    write_checkpoint_modules.ghost_plan = &*ghost_plan;
    write_checkpoint_modules.wall_plan = &*wall_plan;
    write_checkpoint_modules.transform = &*transform;
    write_checkpoint_modules.flow = &facade;
    read_checkpoint_modules.surface = &*surface;
    read_checkpoint_modules.query = &*query;
    read_checkpoint_modules.domain = &*domain;
    read_checkpoint_modules.ghost_plan = &*ghost_plan;
    read_checkpoint_modules.wall_plan = &*wall_plan;
    read_checkpoint_modules.transform = &*transform;
    read_checkpoint_modules.flow = &facade;
  }

  const auto dispatch_status = runtime::collective_status(
      mpi, stage3_dispatch_inventory_complete(write_checkpoint_modules),
      "Stage 3 dispatch inventory is incomplete");
  if (!dispatch_status.ok)
    throw Error(dispatch_status.message);

  double proposed_dt = config.time.initial_dt_s;
  std::uint32_t last_retry_count{};
  if (config.restart.read) {
    const auto restored = flow::read_checkpoint_v3(
        mpi, decomposition, topology, geometry, boundaries, immersed_config,
        read_checkpoint_modules, state,
        beneath(authoritative_case_root, *config.restart.read_directory));
    if (!restored.restored()) {
      throw Error("immersed-flow checkpoint read failed: reason=" +
                  std::to_string(static_cast<int>(restored.report().reason())) +
                  " phase=" +
                  std::to_string(static_cast<int>(restored.report().phase())) +
                  " rank=" +
                  std::to_string(restored.report().lowest_failing_rank()));
    }
    proposed_dt = restored.control_state().proposed_next_dt_s;
    last_retry_count = restored.control_state().last_retry_count;
    if (state.metadata().step >
        static_cast<std::uint64_t>(config.time.steps)) {
      throw Error("restored immersed-flow step exceeds configured target");
    }
  }

  const auto target_step = static_cast<std::uint64_t>(config.time.steps);
  std::optional<diagnostics::DiagnosticSession> diagnostic_session;
  std::optional<diagnostics::ImmersedStaticDiagnosticSummary> static_summary;
  if (has_future_diagnostic_step(state.metadata().step, target_step,
                                 config.diagnostics.write_interval)) {
    diagnostic_session.emplace(
        beneath(authoritative_case_root, config.diagnostics.directory),
        config.diagnostics.write_interval, mpi.rank());
    if (has_immersed) {
      static_summary.emplace(diagnostics::summarize_immersed_static(
          *surface, *query, *domain, *ghost_plan, *wall_plan, *transform));
    }
  }
  const auto checkpoint_root =
      beneath(authoritative_case_root, config.restart.write_directory);
  const bool writes_checkpoint =
      config.restart.write_interval > 0 &&
      config.time.steps >= config.restart.write_interval;
  if (writes_checkpoint)
    ensure_directory(mpi, checkpoint_root);

  if (mpi.rank() == 0) {
    std::cout << "CASE name=" << config.case_name << " ranks=" << mpi.size()
              << " cells=" << config.mesh.cells.x << 'x' << config.mesh.cells.y
              << 'x' << config.mesh.cells.z
              << " immersed_boundary="
              << (has_immersed ? "local_flow_pattern_ghost_cell" : "none")
              << " density_model=" << density_model_name(config.density_model)
              << " les="
              << (has_wale ? "wale" : "none") << '\n';
  }

  while (state.metadata().step <
         static_cast<std::uint64_t>(config.time.steps)) {
    double attempt_dt = proposed_dt;
    bool committed = false;
    for (int retry = 0; retry <= config.time.max_retries; ++retry) {
      const auto before = state.metadata();
      const auto order = before.step == 0U
                             ? flow::MomentumTimeOrder::backward_euler
                             : flow::MomentumTimeOrder::bdf2;
      const auto stencil = flow::make_momentum_time_stencil(
          order, attempt_dt, before.step == 0U ? 0.0 : before.dt_s);
      auto report = facade.attempt(state, physics, stencil, {}, {});
      const flow::MaterialDensityStepAttemptReport* material = nullptr;
      const flow::IdealGasStepAttemptReport* ideal = nullptr;
      const flow::StepAttemptReport* base = nullptr;
      if (const auto* constant =
              std::get_if<flow::StepAttemptReport>(&report.base)) {
        base = constant;
      } else if (const auto* candidate =
                     std::get_if<flow::MaterialDensityStepAttemptReport>(
                         &report.base)) {
        material = candidate;
        base = &candidate->flow();
      } else if (const auto* ideal_candidate =
                     std::get_if<flow::IdealGasStepAttemptReport>(
                         &report.base)) {
        ideal = ideal_candidate;
        material = &ideal_candidate->flow();
        base = &ideal_candidate->flow().flow();
      } else {
        throw Error("immersed-flow driver received an unexpected report type");
      }
      if (base->disposition == flow::StepAttemptDisposition::committed) {
        require_committed_report(*base, report, has_immersed, has_wale,
                                 material, ideal);
        if (diagnostic_session &&
            diagnostic_session->due(state.metadata().step)) {
          diagnostics::DiagnosticBatch batch;
          std::optional<flow::ImmersedFlowDiagnosticSource> source;
          std::optional<diagnostics::ImmersedStaticDiagnosticSummary>
              attempt_static_summary;
          if (has_immersed) {
            source.emplace(facade.diagnostic_source(state, report));
            attempt_static_summary.emplace(
                diagnostics::with_local_flow_pattern_snapshot(
                    *static_summary, *source));
          }
          append_stage3_diagnostics(
              mpi, state, presence,
              attempt_static_summary ? &*attempt_static_summary : nullptr,
              source ? &*source : nullptr, report, batch);
          diagnostic_session->publish(mpi, state.metadata().step, batch);
        }
        print_step(state, static_cast<std::uint32_t>(retry + 1), *base, report,
                   mpi);
        proposed_dt = next_accepted_dt(config.time, attempt_dt, *base);
        last_retry_count = static_cast<std::uint32_t>(retry);
        committed = true;
        break;
      }
      if (base->disposition !=
          flow::StepAttemptDisposition::recoverable_failure) {
        throw Error("immersed-flow step failed without retry: reason=" +
                    std::to_string(static_cast<int>(base->reason)) +
                    " rank=" + std::to_string(base->lowest_failing_rank));
      }
      if (retry == config.time.max_retries) {
        throw Error("immersed-flow step exhausted retry limit");
      }
      const double reduced = next_attempt_dt(config.time, attempt_dt, *base);
      if (!(reduced < attempt_dt)) {
        throw Error("immersed-flow step reached minimum dt");
      }
      attempt_dt = reduced;
    }
    if (!committed) {
      throw Error("immersed-flow step did not commit");
    }
    if (writes_checkpoint &&
        state.metadata().step %
                static_cast<std::uint64_t>(config.restart.write_interval) ==
            0U) {
      const auto checkpoint = flow::write_checkpoint_v3(
          mpi, decomposition, topology, geometry, boundaries, immersed_config,
          write_checkpoint_modules, state, {proposed_dt, last_retry_count},
          checkpoint_step_path(checkpoint_root, state.metadata().step));
      if (checkpoint.disposition() !=
          flow::CheckpointV3Disposition::completed) {
        throw Error("immersed-flow checkpoint write failed: reason=" +
                    std::to_string(static_cast<int>(checkpoint.reason())) +
                    " phase=" +
                    std::to_string(static_cast<int>(checkpoint.phase())) +
                    " rank=" +
                    std::to_string(checkpoint.lowest_failing_rank()));
      }
    }
  }

  if (mpi.rank() == 0) {
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "FINISHED step=" << state.metadata().step
              << " time_s=" << state.metadata().time_s << '\n';
  }
  mpi.barrier();
  return EXIT_SUCCESS;
}

} // namespace hundun::application
