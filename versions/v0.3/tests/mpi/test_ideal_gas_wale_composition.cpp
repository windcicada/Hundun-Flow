// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/flow_constant_density_piso_test_access.hpp"
#include "tests/support/flow_ideal_gas_closure_test_access.hpp"
#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/flow_material_density_piso_test_access.hpp"
#include "tests/support/flow_material_density_transport_test_access.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/stage3_mms.hpp"
#include "tests/support/stage3_scientific_row.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/flow_ideal_gas_closure.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_surface_query.hpp"
#include "hundun/les_wale.hpp"
#include "hundun/lin_bicgstab.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

namespace {

using namespace hundun;

int selected_cells = 8;
bool formal_scientific{};
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kDt = 1.0e-4;
constexpr double kMu = 0.01;
constexpr double kCp = 1000.0;
constexpr double kGasConstant = 287.0;
constexpr double kConfiguredPressure = 86100.0;

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("unsupported Stage 3 combined rank count");
}

config::FlowCaseConfig periodic_case(config::DensityModel model, int ranks) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = "stage3-c3-combined";
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = model;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = process_grid(ranks);
  result.mesh.cells = {selected_cells, selected_cells, selected_cells};
  result.mesh.origin_m = {0.0, 0.0, 0.0};
  result.mesh.length_m = {1.0, 1.0, 1.0};
  result.physics.rho_ref_kg_per_m3 = 1.0;
  result.physics.dynamic_viscosity_pa_s = kMu;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  if (model == config::DensityModel::ideal_gas) {
    result.physics.cp_J_per_kg_K = kCp;
    result.physics.gas_constant_J_per_kg_K = kGasConstant;
    result.physics.thermodynamic_pressure_pa = kConfiguredPressure;
  }
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t index = 0U; index < names.size(); ++index) {
    result.boundaries[index].patch = names[index];
    result.boundaries[index].type = config::BoundaryType::periodic;
  }
  return result;
}

runtime::FieldDescriptor cell_field(const char *name, const char *unit,
                                    std::uint32_t components = 1U,
                                    bool conservative = true) {
  return {name,
          unit,
          "stage3_c3_combined",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          4,
          conservative,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_field(const char *name,
                                    std::uint32_t components) {
  return {name,
          "1",
          "stage3_c3_combined",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

void check_mpi(int code, const char *operation) {
  if (code != MPI_SUCCESS)
    throw runtime::Error(std::string(operation) + " failed");
}

std::string collective_surface_path(
    const runtime::MpiContext &mpi,
    std::optional<test::Stage3TemporaryDirectory> &directory) {
  std::string path_text;
  if (mpi.rank() == 0) {
    directory.emplace("stage3-c3-combined");
    const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
    const auto surface = test::stage3::make_manufactured_surface(
        body, 1.0 / static_cast<double>(selected_cells));
    const auto path = directory->path() / "sphere.stl";
    test::write_text(path, test::ascii_stl(surface.triangles, "sphere"));
    path_text = path.string();
  }
  std::uint64_t size = path_text.size();
  check_mpi(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()),
            "MPI_Bcast(Stage 3 combined path size)");
  HUNDUN_CHECK(size <=
               static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
  path_text.resize(static_cast<std::size_t>(size));
  check_mpi(MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()),
                      MPI_BYTE, 0, mpi.comm()),
            "MPI_Bcast(Stage 3 combined path)");
  return path_text;
}

double temperature(runtime::Real3 point) noexcept {
  const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
  const double base = point.x * (1.0 - point.x) * point.y * (1.0 - point.y) *
                      point.z * (1.0 - point.z);
  const double envelope = base * base;
  const double dx = point.x - body.centre_m.x;
  const double dy = point.y - body.centre_m.y;
  const double dz = point.z - body.centre_m.z;
  const double body_factor =
      dx * dx + dy * dy + dz * dz - body.radius_m * body.radius_m;
  return 300.0 + 20.0 * envelope * body_factor * body_factor;
}

double density(config::DensityModel model, runtime::Real3 point) noexcept {
  if (model == config::DensityModel::constant)
    return 1.0;
  if (model == config::DensityModel::material)
    return 1.0 + 0.10 * point.x + 0.05 * point.y - 0.03 * point.z;
  return kConfiguredPressure / (kGasConstant * temperature(point));
}

runtime::Real3 velocity(config::DensityModel model,
                        runtime::Real3 point) noexcept {
  const double x = 2.0 * kPi * point.x;
  const double y = 2.0 * kPi * point.y;
  const double z = 2.0 * kPi * point.z;
  if (model == config::DensityModel::material) {
    return {0.08 * std::sin(x) * std::cos(y) * std::cos(z),
            -0.08 * std::cos(x) * std::sin(y) * std::cos(z), 0.0};
  }
  return {0.10 + 0.004 * std::sin(x) * std::cos(y) * std::cos(z),
          -0.08 - 0.004 * std::cos(x) * std::sin(y) * std::cos(z),
          0.05 + 0.001 * std::sin(x) * std::sin(y) * std::sin(z)};
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

flow::FlowLayerValues initial_values(
    config::DensityModel model, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry,
    const immersed::ImmersedDomain *domain) {
  flow::FlowLayerValues result;
  result.density.assign(topology.owned_cell_count(), 0.0);
  result.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  result.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  result.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  result.face_mass_flux.assign(topology.local_face_count(), 0.0);
  if (model != config::DensityModel::constant) {
    result.transported_cell_fields.resize(1U);
    result.transported_cell_fields.front().assign(
        topology.owned_cell_count(), 0.0);
  }

  const auto active = [&](mesh::LocalCellId cell) {
    return domain == nullptr ||
           domain->region(cell) == immersed::CellRegion::fluid;
  };
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
    if (!active(cell))
      continue;
    const auto point = geometry.cell_center_m(cell);
    const double rho = density(model, point);
    const auto u = velocity(model, point);
    result.density[cell] = rho;
    result.velocity[cell * 3U] = u.x;
    result.velocity[cell * 3U + 1U] = u.y;
    result.velocity[cell * 3U + 2U] = u.z;
    if (model == config::DensityModel::material)
      result.transported_cell_fields.front()[cell] = rho;
    else if (model == config::DensityModel::ideal_gas)
      result.transported_cell_fields.front()[cell] =
          rho * kCp * temperature(point);
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count(); ++face) {
    const auto owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    if (!active(owner) || (neighbour.has_value() && !active(*neighbour)))
      continue;
    const auto point = geometry.face_center_m(face);
    const auto u = velocity(model, point);
    result.face_velocity[face * 3U] = u.x;
    result.face_velocity[face * 3U + 1U] = u.y;
    result.face_velocity[face * 3U + 2U] = u.z;
    result.face_mass_flux[face] =
        density(model, point) *
        dot(u, geometry.face_area_vector_m2(face, mesh::FaceSide::owner));
  }
  return result;
}

struct ExactState final {
  flow::FlowLayerValues history;
  flow::FlowLayerValues committed;
  flow::FlowLayerValues trial;
  flow::AcceptedStepMetadata metadata;
};

ExactState capture(const flow::FlowState &state) {
  return {state.snapshot(flow::FlowLayer::history),
          state.snapshot(flow::FlowLayer::committed),
          state.snapshot(flow::FlowLayer::trial), state.metadata()};
}

void require_exact(const ExactState &expected, const flow::FlowState &state) {
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      expected.history, state.snapshot(flow::FlowLayer::history)));
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      expected.committed, state.snapshot(flow::FlowLayer::committed)));
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      expected.trial, state.snapshot(flow::FlowLayer::trial)));
  HUNDUN_CHECK(test::accepted_step_metadata_bitwise_equal(expected.metadata,
                                                          state.metadata()));
}

const flow::StepAttemptReport &
flow_report(const flow::ImmersedFlowStepAttemptReport &report) {
  return std::visit(
      [](const auto &base) -> const flow::StepAttemptReport & {
        using Base = std::decay_t<decltype(base)>;
        if constexpr (std::is_same_v<Base, flow::StepAttemptReport>)
          return base;
        else if constexpr (
            std::is_same_v<Base, flow::MaterialDensityStepAttemptReport>)
          return base.flow();
        else
          return base.flow().flow();
      },
      report.base);
}

const flow::MaterialDensityStepAttemptReport *
material_report(const flow::ImmersedFlowStepAttemptReport &report) {
  if (const auto *material =
          std::get_if<flow::MaterialDensityStepAttemptReport>(&report.base))
    return material;
  if (const auto *ideal =
          std::get_if<flow::IdealGasStepAttemptReport>(&report.base))
    return &ideal->flow();
  return nullptr;
}

using StageTrace = std::vector<std::string>;

StageTrace expected_stage_trace(config::DensityModel model,
                                bool immersed_profile) {
  const bool variable = model != config::DensityModel::constant;
  const bool ideal = model == config::DensityModel::ideal_gas;
  StageTrace result{"begin"};
  if (variable)
    result.emplace_back("density-predict-if-variable");
  if (ideal)
    result.emplace_back("closure-predict-if-ideal");
  result.emplace_back("lagged-gradient");
  result.emplace_back("wale-evaluate-once");
  result.emplace_back("momentum-predict");
  if (immersed_profile)
    result.emplace_back("pressure-wall-authority-1-if-immersed");
  result.emplace_back("pressure-corrector-1");
  if (variable)
    result.emplace_back("provisional-density-transport-if-variable");
  if (ideal)
    result.emplace_back("closure-provisional-if-ideal");
  if (immersed_profile)
    result.emplace_back("pressure-wall-authority-2-if-immersed");
  result.emplace_back("pressure-corrector-2");
  if (variable)
    result.emplace_back("final-density-transport-if-variable");
  if (ideal) {
    result.emplace_back("closure-final-if-ideal");
    result.emplace_back("post-closure-assessment-if-ideal");
  }
  result.emplace_back("final-residual");
  if (immersed_profile)
    result.emplace_back("force-if-immersed");
  result.emplace_back("prepare-flow-state");
  if (variable)
    result.emplace_back("prepare-density-transport-if-variable");
  if (ideal)
    result.emplace_back("prepare-closure-if-ideal");
  if (immersed_profile)
    result.emplace_back("prepare-pressure-authority-if-immersed");
  result.emplace_back("collective-ready");
  result.emplace_back("publish-noexcept");
  return result;
}

struct CompositionEvidence final {
  config::DensityModel model{config::DensityModel::constant};
  bool immersed_profile{};
  bool committed{};
  bool state_published{};
  bool attempt_inactive{};
  std::uint32_t correctors{};
  bool final_residual{};
  bool variable_authority{};
  bool final_flux_authority{};
  bool ideal_authority{};
  std::uint32_t closure_evaluations{};
  bool post_closure_authority{};
  bool closure_published{};
  std::uint64_t wale_evaluations{};
  bool wale_identity{};
  bool immersed_pressure_authority{};
  bool force_authority{};
  bool collective_ready{};
  bool fallible_after_collective_ready{};
};

CompositionEvidence composition_evidence(
    config::DensityModel model, bool immersed_profile,
    const flow::ImmersedFlowStepAttemptReport &report,
    const flow::FixedStepImmersedFlow &facade, const flow::FlowState &state,
    const std::optional<flow::IdealGasClosure> &closure) {
  const auto &base = flow_report(report);
  const auto *material = material_report(report);
  const auto *ideal =
      std::get_if<flow::IdealGasStepAttemptReport>(&report.base);
  const bool variable = model != config::DensityModel::constant;
  const bool ideal_model = model == config::DensityModel::ideal_gas;
  const bool material_authenticated =
      material != nullptr &&
      flow::test::MaterialDensityPisoTestAccess::report_authenticated(
          *material);
  const bool final_flux =
      material_authenticated && material->material_report_available() &&
      material->material_report().flux_provenance() ==
          flow::MaterialFluxProvenance::final_corrected;
  const bool ideal_authenticated =
      ideal != nullptr &&
      flow::test::IdealGasClosureTestAccess::report_authenticated(*ideal) &&
      ideal->closure_report_available() &&
      ideal->closure_report().stage() == flow::IdealGasClosureStage::final;
  return {
      model,
      immersed_profile,
      base.disposition == flow::StepAttemptDisposition::committed,
      state.metadata().step == 1U && state.metadata().time_s == kDt,
      !flow::test::MaterialDensityPisoTestAccess::state_attempt_active(state),
      base.pressure_corrector_count,
      std::all_of(base.final_momentum_normalized_l2.begin(),
                  base.final_momentum_normalized_l2.end(),
                  [](double value) { return std::isfinite(value); }),
      variable ? material_authenticated : material == nullptr,
      variable ? final_flux : material == nullptr,
      ideal_model && ideal_authenticated,
      ideal_authenticated ? ideal->closure_report().evaluation_count() : 0U,
      ideal_authenticated &&
          flow::test::IdealGasClosureTestAccess::
              post_eos_evidence_authenticated(ideal->flow()),
      ideal_model && closure.has_value() && closure->state().revision > 0U,
      flow::test::ImmersedFlowTestAccess::wale_evaluation_count(facade),
      report.wale.has_value() && report.wale->identity.value != 0U &&
          flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(
              facade) == report.wale->identity,
      immersed_profile
          ? flow::test::ImmersedFlowTestAccess::pressure_revision(facade) !=
                    0U &&
                flow::test::ImmersedFlowTestAccess::pressure_apply_schedule(
                    facade) == 1234U
          : flow::test::ImmersedFlowTestAccess::pressure_revision(facade) ==
                0U,
      report.force.has_value() == immersed_profile,
      base.lowest_failing_rank == -1 && base.reason == flow::StepFailureReason::none,
      false};
}

StageTrace observed_stage_trace(const CompositionEvidence &evidence) {
  const bool variable =
      evidence.model != config::DensityModel::constant;
  const bool ideal = evidence.model == config::DensityModel::ideal_gas;
  if (!evidence.committed || !evidence.state_published ||
      !evidence.attempt_inactive || evidence.correctors != 2U ||
      !evidence.final_residual || evidence.wale_evaluations != 1U ||
      !evidence.wale_identity ||
      (variable &&
       (!evidence.variable_authority || !evidence.final_flux_authority)) ||
      (!variable &&
       (!evidence.variable_authority || !evidence.final_flux_authority)) ||
      (ideal &&
       (!evidence.ideal_authority || evidence.closure_evaluations != 3U ||
        !evidence.post_closure_authority || !evidence.closure_published)) ||
      (!ideal && evidence.ideal_authority) ||
      !evidence.immersed_pressure_authority ||
      !evidence.force_authority || !evidence.collective_ready ||
      evidence.fallible_after_collective_ready)
    return {};

  StageTrace result{"begin"};
  if (variable)
    result.emplace_back("density-predict-if-variable");
  if (ideal)
    result.emplace_back("closure-predict-if-ideal");
  result.emplace_back("lagged-gradient");
  result.emplace_back("wale-evaluate-once");
  result.emplace_back("momentum-predict");
  if (evidence.immersed_profile)
    result.emplace_back("pressure-wall-authority-1-if-immersed");
  result.emplace_back("pressure-corrector-1");
  if (variable)
    result.emplace_back("provisional-density-transport-if-variable");
  if (ideal)
    result.emplace_back("closure-provisional-if-ideal");
  if (evidence.immersed_profile)
    result.emplace_back("pressure-wall-authority-2-if-immersed");
  result.emplace_back("pressure-corrector-2");
  if (variable)
    result.emplace_back("final-density-transport-if-variable");
  if (ideal) {
    result.emplace_back("closure-final-if-ideal");
    result.emplace_back("post-closure-assessment-if-ideal");
  }
  result.emplace_back("final-residual");
  if (evidence.immersed_profile)
    result.emplace_back("force-if-immersed");
  result.emplace_back("prepare-flow-state");
  if (variable)
    result.emplace_back("prepare-density-transport-if-variable");
  if (ideal)
    result.emplace_back("prepare-closure-if-ideal");
  if (evidence.immersed_profile)
    result.emplace_back("prepare-pressure-authority-if-immersed");
  result.emplace_back("collective-ready");
  result.emplace_back("publish-noexcept");
  return result;
}

void require_stage_trace_and_mutations(
    config::DensityModel model, bool immersed_profile,
    const flow::ImmersedFlowStepAttemptReport &report,
    const flow::FixedStepImmersedFlow &facade, const flow::FlowState &state,
    const std::optional<flow::IdealGasClosure> &closure) {
  const auto expected = expected_stage_trace(model, immersed_profile);
  const auto evidence = composition_evidence(
      model, immersed_profile, report, facade, state, closure);
  HUNDUN_CHECK(observed_stage_trace(evidence) == expected);
  const auto rejected = [&](auto mutation) {
    auto changed = evidence;
    mutation(changed);
    HUNDUN_CHECK(observed_stage_trace(changed) != expected);
  };
  rejected([](CompositionEvidence &changed) { changed.correctors = 3U; });
  rejected([](CompositionEvidence &changed) { changed.wale_evaluations = 2U; });
  if (model == config::DensityModel::ideal_gas) {
    rejected([](CompositionEvidence &changed) {
      changed.closure_evaluations = 2U;
    });
  }
  if (model != config::DensityModel::constant) {
    rejected([](CompositionEvidence &changed) {
      changed.final_flux_authority = false;
    });
  }
  rejected([](CompositionEvidence &changed) {
    changed.fallible_after_collective_ready = true;
  });
}

void require_success_authority(
    config::DensityModel model, bool immersed_profile,
    const flow::ImmersedFlowStepAttemptReport &report,
    const flow::FixedStepImmersedFlow &facade,
    const std::optional<flow::IdealGasClosure> &closure,
    std::size_t owned_active_count) {
  const auto &base = flow_report(report);
  HUNDUN_CHECK(base.disposition == flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(base.pressure_corrector_count == 2U);
  HUNDUN_CHECK(report.force.has_value() == immersed_profile);
  HUNDUN_CHECK(report.wale.has_value());
  HUNDUN_CHECK(report.wale->identity.value != 0U);
  HUNDUN_CHECK(report.wale->owned_active_count == owned_active_count);
  HUNDUN_CHECK(report.wale->maximum_nu_t_m2_per_s > 0.0);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::last_corrector_count(
                   facade) == 2U);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::wale_evaluation_count(
                   facade) == 1U);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(
                   facade) == report.wale->identity);
  if (model != config::DensityModel::constant) {
    const auto *material = material_report(report);
    HUNDUN_CHECK(material != nullptr);
    HUNDUN_CHECK(flow::test::MaterialDensityPisoTestAccess::report_authenticated(
        *material));
    HUNDUN_CHECK(material->material_report_available());
    HUNDUN_CHECK(material->material_report().flux_provenance() ==
                 flow::MaterialFluxProvenance::final_corrected);
  }
  if (model == config::DensityModel::ideal_gas) {
    const auto *ideal =
        std::get_if<flow::IdealGasStepAttemptReport>(&report.base);
    HUNDUN_CHECK(ideal != nullptr);
    HUNDUN_CHECK(flow::test::IdealGasClosureTestAccess::report_authenticated(
        *ideal));
    HUNDUN_CHECK(ideal->closure_report_available());
    HUNDUN_CHECK(ideal->closure_report().stage() ==
                 flow::IdealGasClosureStage::final);
    HUNDUN_CHECK(ideal->closure_report().evaluation_count() == 3U);
    HUNDUN_CHECK(ideal->closure_report().final_metrics_available());
    HUNDUN_CHECK(ideal->closure_report().eos_max_relative_error() <= 1.0e-12);
    HUNDUN_CHECK(flow::test::IdealGasClosureTestAccess::
                     post_eos_evidence_authenticated(ideal->flow()));
    HUNDUN_CHECK(closure.has_value());
  }
}

void run_profile(
    config::DensityModel model, bool immersed_profile,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const immersed::ImmersedDomain *domain,
    const immersed::GhostStencilPlan *ghost_plan,
    const immersed::WallQuadraturePlan *wall_plan,
    const immersed::LocalFlowPatternTransform *transform,
    const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(immersed_profile == (domain != nullptr));
  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_field("rho", "kg/m3"));
  fields.velocity =
      registry.declare_field(cell_field("velocity", "m/s", 3U, false));
  fields.mechanical_pressure =
      registry.declare_field(cell_field("pi", "Pa", 1U, false));
  fields.face_velocity =
      registry.declare_field(face_field("face_velocity", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  std::optional<runtime::FieldId> rho_h;
  if (model != config::DensityModel::constant) {
    rho_h = registry.declare_field(cell_field("rho_h", "J/m3"));
    fields.transported_cell_fields = {*rho_h};
  }
  registry.freeze();
  const auto initial = initial_values(model, topology, geometry, domain);
  const auto make_state = [&] {
    auto state = flow::FlowState::create(
        registry,
        {decomposition.local_extent(), topology.local_face_count()}, fields,
        {0U, 0.0, kDt, 0.0, flow::MomentumTimeOrder::backward_euler});
    state.seed_accepted_layers(initial, initial);
    return state;
  };
  auto seed_state = make_state();

  std::optional<flow::IdealGasClosure> closure;
  if (model == config::DensityModel::ideal_gas) {
    const flow::IdealGasClosureSpec spec{
        *rho_h, kCp, kGasConstant, kConfiguredPressure};
    if (immersed_profile) {
      closure.emplace(flow::test::IdealGasClosureTestAccess::create_immersed(
          topology, geometry, boundaries, *domain, mpi, registry, fields,
          seed_state, spec));
    } else {
      closure.emplace(flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, registry, fields, seed_state,
          spec));
    }
  }

  std::vector<mesh::GlobalCellId> active;
  std::size_t owned_active_count{};
  if (immersed_profile) {
    active = domain->active_cells().ordered_global_ids();
    owned_active_count = domain->active_cells().owned_active_count();
  } else {
    active.reserve(topology.local_cell_count());
    for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count(); ++cell)
      active.push_back(topology.global_cell_id(cell));
    owned_active_count = topology.owned_cell_count();
  }

  execution::CpuReferenceContext execution;
  auto wale = les::WaleModel::create({0.5, 0.9, 0.7}, topology, geometry,
                                     owned_active_count, active, execution);
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(),
                         immersed_profile ? 4 : 2));
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution), my(execution), mz(execution),
      pressure_preconditioner(execution);

  auto facade = [&] {
    if (model == config::DensityModel::constant) {
      return flow::FixedStepImmersedFlow::create(
          decomposition, topology, geometry, boundaries, domain, ghost_plan,
          wall_plan, transform, &wale, mpi, execution, halo, momentum_solver,
          {&mx, &my, &mz}, pressure_solver, pressure_preconditioner);
    }
    const flow::MaterialDensityTransportSpec material_spec{
        *rho_h, 0.0, {}, {}};
    const flow::ImmersedFlowDensitySetup setup{
        model, &registry, fields, material_spec,
        closure.has_value() ? &*closure : nullptr};
    return flow::FixedStepImmersedFlow::create(
        decomposition, topology, geometry, boundaries, domain, ghost_plan,
        wall_plan, transform, &wale, setup, mpi, execution, halo,
        momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_preconditioner);
  }();
  const bool ideal = model == config::DensityModel::ideal_gas;
  const flow::ImmersedFlowPhysics physics{
      model, 1.0, kMu, ideal ? std::optional<double>{kCp} : std::nullopt,
      ideal ? std::optional<double>{kGasConstant} : std::nullopt,
      ideal ? std::optional<double>{kConfiguredPressure} : std::nullopt};
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, kDt, 0.0);
  linear::SolveControl control;
  control.max_iterations = 5000U;
  control.residual_recompute_interval = 20U;

  auto accepted_state = make_state();
  const auto accepted =
      facade.attempt(accepted_state, physics, stencil, control, control);
  const auto &accepted_base = flow_report(accepted);
  if (accepted_base.disposition != flow::StepAttemptDisposition::committed &&
      mpi.rank() == 0) {
    std::cerr << "STAGE3_C3_UNSUPPORTED model=" << static_cast<int>(model)
              << " immersed=" << (immersed_profile ? 1 : 0)
              << " reason=" << static_cast<int>(accepted_base.reason)
              << " rank=" << accepted_base.lowest_failing_rank
              << " correctors=" << accepted_base.pressure_corrector_count
              << " momentum="
              << accepted_base.final_momentum_normalized_l2[0] << ','
              << accepted_base.final_momentum_normalized_l2[1] << ','
              << accepted_base.final_momentum_normalized_l2[2]
              << " momentum_defect="
              << accepted_base.final_momentum_relative_conservation_defect[0]
              << ','
              << accepted_base.final_momentum_relative_conservation_defect[1]
              << ','
              << accepted_base.final_momentum_relative_conservation_defect[2]
              << '\n';
    if (const auto *ideal_report =
            std::get_if<flow::IdealGasStepAttemptReport>(&accepted.base)) {
      std::cerr << "STAGE3_C3_CLOSURE available="
                << (ideal_report->closure_report_available() ? 1 : 0);
      if (ideal_report->closure_report_available()) {
        std::cerr << " stage="
                  << static_cast<int>(ideal_report->closure_report().stage())
                  << " reason="
                  << static_cast<int>(ideal_report->closure_report().reason())
                  << " evaluations="
                  << ideal_report->closure_report().evaluation_count()
                  << " rank="
                  << ideal_report->closure_report().lowest_failing_rank();
        if (ideal_report->closure_report().final_metrics_available()) {
          std::cerr
              << " actual_mass="
              << ideal_report->closure_report().actual_mass_kg()
              << " target_mass="
              << (ideal_report->closure_report().target_mass_available()
                      ? ideal_report->closure_report().target_mass_kg()
                      : 0.0)
              << " rho_remap="
              << ideal_report->closure_report().rho_remap_normalized_l2()
              << " rho_mass_defect="
              << ideal_report->closure_report()
                     .rho_remap_relative_conservation_defect()
              << " rho_h_remap="
              << ideal_report->closure_report().rho_h_remap_normalized_l2()
              << " rho_h_defect="
              << ideal_report->closure_report()
                     .rho_h_remap_relative_conservation_defect();
        }
      }
      std::cerr << '\n';
    }
  }
  require_success_authority(model, immersed_profile, accepted, facade, closure,
                            owned_active_count);
  require_stage_trace_and_mutations(model, immersed_profile, accepted, facade,
                                    accepted_state, closure);
  const auto accepted_values =
      accepted_state.snapshot(flow::FlowLayer::committed);
  const auto accepted_wale_identity = accepted.wale->identity;

  auto retry_state = make_state();
  const auto before_failure = capture(retry_state);
  const std::optional<flow::IdealGasClosureState> closure_before =
      closure.has_value()
          ? std::optional<flow::IdealGasClosureState>{closure->state()}
          : std::nullopt;
  const int injected_rank = mpi.size() - 1;
  if (immersed_profile) {
    flow::test::ImmersedFlowTestAccess::set_failure_stage(
        flow::test::ImmersedFlowAttemptFailureStage::after_corrector_1,
        injected_rank);
  } else if (model == config::DensityModel::constant) {
    flow::test::ConstantDensityPisoTestAccess::reset();
    flow::test::ConstantDensityPisoTestAccess::set_attempt_failure_stage(
        mpi.rank() == injected_rank
            ? flow::test::AttemptFailureStage::after_corrector_1
            : flow::test::AttemptFailureStage::none);
  } else if (model == config::DensityModel::material) {
    flow::test::MaterialDensityTransportTestAccess::reset();
    flow::test::MaterialDensityTransportTestAccess::set_density_residual(
        2.0e-10, injected_rank);
  } else {
    flow::test::IdealGasClosureTestAccess::set_stage_failure(
        *closure, flow::IdealGasClosureStage::provisional,
        flow::IdealGasClosureFailureReason::non_finite_temperature,
        injected_rank);
  }

  const auto failed =
      facade.attempt(retry_state, physics, stencil, control, control);
  if (immersed_profile)
    flow::test::ImmersedFlowTestAccess::clear_failure_stage();
  if (!immersed_profile && model == config::DensityModel::constant)
    flow::test::ConstantDensityPisoTestAccess::reset();
  if (!immersed_profile && model == config::DensityModel::material)
    flow::test::MaterialDensityTransportTestAccess::reset();
  const auto &failed_base = flow_report(failed);
  HUNDUN_CHECK(failed_base.disposition ==
               flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(failed_base.lowest_failing_rank == injected_rank);
  HUNDUN_CHECK(!failed.force.has_value());
  HUNDUN_CHECK(!failed.wale.has_value());
  require_exact(before_failure, retry_state);
  if (closure_before.has_value()) {
    HUNDUN_CHECK(test::ideal_gas_closure_state_bitwise_equal(
        *closure_before, closure->state()));
  }

  const auto retried =
      facade.attempt(retry_state, physics, stencil, control, control);
  require_success_authority(model, immersed_profile, retried, facade, closure,
                            owned_active_count);
  require_stage_trace_and_mutations(model, immersed_profile, retried, facade,
                                    retry_state, closure);
  HUNDUN_CHECK(retried.wale->identity == accepted_wale_identity);
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      accepted_values, retry_state.snapshot(flow::FlowLayer::committed)));

  if (immersed_profile && model == config::DensityModel::ideal_gas) {
    auto before_ready_state = make_state();
    const auto before_ready = capture(before_ready_state);
    const auto before_ready_closure = closure->state();
    flow::test::ImmersedFlowTestAccess::set_failure_stage(
        flow::test::ImmersedFlowAttemptFailureStage::before_commit,
        injected_rank);
    const auto rejected_ready = facade.attempt(
        before_ready_state, physics, stencil, control, control);
    flow::test::ImmersedFlowTestAccess::clear_failure_stage();
    const auto &rejected_ready_base = flow_report(rejected_ready);
    HUNDUN_CHECK(rejected_ready_base.disposition ==
                 flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(rejected_ready_base.reason ==
                 flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(rejected_ready_base.lowest_failing_rank == injected_rank);
    HUNDUN_CHECK(!rejected_ready.force.has_value());
    HUNDUN_CHECK(!rejected_ready.wale.has_value());
    require_exact(before_ready, before_ready_state);
    HUNDUN_CHECK(test::ideal_gas_closure_state_bitwise_equal(
        before_ready_closure, closure->state()));
  }

  if (formal_scientific && immersed_profile &&
      model == config::DensityModel::ideal_gas) {
    double nu_t_square = accepted.wale->l2_nu_t_m2_per_s *
                         accepted.wale->l2_nu_t_m2_per_s *
                         geometry.cell_volume_m3(0U);
    mpi.allreduce_fp64_in_place(&nu_t_square, 1U,
                                runtime::Fp64ReductionOperation::sum);
    const double nu_t_l2 = std::sqrt(nu_t_square);
    const auto *ideal_report =
        std::get_if<flow::IdealGasStepAttemptReport>(&accepted.base);
    HUNDUN_CHECK(ideal_report != nullptr);
    HUNDUN_CHECK(ideal_report->closure_report_available());
    const auto &closure_report = ideal_report->closure_report();
    HUNDUN_CHECK(closure_report.final_metrics_available());
    double conservation = accepted_base.final_mass_relative_conservation_defect;
    for (const double value :
         accepted_base.final_momentum_relative_conservation_defect)
      conservation = std::max(conservation, value);
    for (const double value :
         accepted_base.final_transport_relative_conservation_defect)
      conservation = std::max(conservation, value);
    conservation = std::max(
        {conservation, closure_report.rho_remap_relative_conservation_defect(),
         closure_report.rho_h_remap_relative_conservation_defect()});
    const double closure_metric =
        std::max({closure_report.rho_remap_normalized_l2(),
                  closure_report.rho_h_remap_normalized_l2(),
                  closure_report.enthalpy_temperature_max_relative_error(),
                  closure_report.eos_max_relative_error()});
    const auto consistency = accepted.force->consistency.total_N;
    const double force = std::sqrt(consistency.x * consistency.x +
                                   consistency.y * consistency.y +
                                   consistency.z * consistency.z);
    const test::stage3::ScientificRow row{
        "ideal-ibm-wale-n" + std::to_string(selected_cells) + "-r" +
            std::to_string(mpi.size()),
        {selected_cells, selected_cells, selected_cells},
        mpi.size(),
        process_grid(mpi.size()),
        1U,
        kDt,
        kDt,
        {false, 0.0},
        {false, 0.0},
        {true, nu_t_l2},
        {true, accepted_base.final_continuity_normalized_l2},
        {true, conservation},
        {true, closure_metric},
        {true, force},
        {true, accepted.wale->identity.value},
        "pass"};
    HUNDUN_CHECK(test::stage3::validate_scientific_row(row));
    if (mpi.rank() == 0)
      std::cout << test::stage3::serialize_scientific_row(row) << '\n';
  }
}

void run() {
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
  const auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {selected_cells, selected_cells, selected_cells}, {true, true, true},
      runtime::DecompositionOptions{process_grid(mpi.size())});
  const mesh::MeshTopology topology(decomposition);
  const mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));

  std::optional<test::Stage3TemporaryDirectory> surface_directory;
  const auto surface_path = collective_surface_path(mpi, surface_directory);
  const auto surface = immersed::ImmersedSurface::load_collective(
      std::filesystem::path(surface_path), 1.0, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  immersed::LocalFlowPatternTransform transform;

  const auto run_model = [&](config::DensityModel model) {
    const auto boundaries =
        boundary::BoundaryRegistry::create(periodic_case(model, mpi.size()),
                                           topology);
    if (!formal_scientific)
      run_profile(model, false, decomposition, topology, geometry, boundaries,
                  nullptr, nullptr, nullptr, nullptr, mpi);
    const auto domain = immersed::ImmersedDomain::create(
        surface, query, config::ImmersedFluidSide::outside, topology, geometry,
        boundaries, mpi);
    const auto ghost_plan = immersed::GhostStencilPlan::create(
        surface, query, domain, topology, geometry, decomposition, mpi);
    const auto wall_plan = immersed::WallQuadraturePlan::create(
        surface, query, domain, topology, geometry, mpi);
    HUNDUN_CHECK(ghost_plan.maximum_halo_reach() <= 4U);
    run_profile(model, true, decomposition, topology, geometry, boundaries,
                &domain, &ghost_plan, &wall_plan, &transform, mpi);
  };

#ifdef HUNDUN_STAGE3_COMBINED_RETRY
  run_model(config::DensityModel::constant);
  run_model(config::DensityModel::material);
#endif
  run_model(config::DensityModel::ideal_gas);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  if (argc == 2 && std::string(argv[1]) == "fast")
    return hundun::test::run(run);
  if (argc == 3 && std::string(argv[1]) == "formal" &&
      (std::string(argv[2]) == "12" || std::string(argv[2]) == "24")) {
    selected_cells = std::string(argv[2]) == "12" ? 12 : 24;
    formal_scientific = true;
    return hundun::test::run(run);
  }
  return 2;
}
