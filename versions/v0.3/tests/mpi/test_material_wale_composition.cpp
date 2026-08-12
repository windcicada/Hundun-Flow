// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/flow_material_density_piso_test_access.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/stage3_mms.hpp"
#include "tests/support/stage3_scientific_row.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/exec_execution.hpp"
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
#include <vector>

namespace {

using namespace hundun;

int selected_cells = 8;
bool formal_scientific{};
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kDt = 1.0e-4;
constexpr double kMu = 0.01;

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("unsupported material-WALE rank count");
}

config::FlowCaseConfig periodic_case(int ranks) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = "stage3-c2-material-wale";
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::material;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = process_grid(ranks);
  result.mesh.cells = {selected_cells, selected_cells, selected_cells};
  result.mesh.origin_m = {0.0, 0.0, 0.0};
  result.mesh.length_m = {1.0, 1.0, 1.0};
  result.physics.rho_ref_kg_per_m3 = 1.0;
  result.physics.dynamic_viscosity_pa_s = kMu;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  result.scalars.push_back({"alpha", 0.0});
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
          "stage3_c2_material_wale",
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
          "stage3_c2_material_wale",
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
    directory.emplace("stage3-c2-material-wale");
    const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
    const auto surface = test::stage3::make_manufactured_surface(
        body, 1.0 / static_cast<double>(selected_cells));
    const auto path = directory->path() / "sphere.stl";
    test::write_text(path, test::ascii_stl(surface.triangles, "sphere"));
    path_text = path.string();
  }
  std::uint64_t size = path_text.size();
  check_mpi(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()),
            "MPI_Bcast(material-WALE path size)");
  HUNDUN_CHECK(size <=
               static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
  path_text.resize(static_cast<std::size_t>(size));
  check_mpi(MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()),
                      MPI_BYTE, 0, mpi.comm()),
            "MPI_Bcast(material-WALE path)");
  return path_text;
}

double density(runtime::Real3 point) noexcept {
  return 1.0 + 0.10 * point.x + 0.05 * point.y - 0.03 * point.z;
}

runtime::Real3 velocity(runtime::Real3 point) noexcept {
  const double x = 2.0 * kPi * point.x;
  const double y = 2.0 * kPi * point.y;
  const double z = 2.0 * kPi * point.z;
  return {0.08 * std::sin(x) * std::cos(y) * std::cos(z),
          -0.08 * std::cos(x) * std::sin(y) * std::cos(z),
          0.0};
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

bool same(runtime::Int3 left, runtime::Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

flow::FlowLayerValues initial_values(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::ImmersedDomain *domain, double seeded_mass_flux,
    const runtime::MpiContext &mpi) {
  flow::FlowLayerValues result;
  result.density.assign(topology.owned_cell_count(), 0.0);
  result.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  result.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  result.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  result.face_mass_flux.assign(topology.local_face_count(), 0.0);
  result.transported_cell_fields.resize(2U);
  for (auto &field : result.transported_cell_fields)
    field.assign(topology.owned_cell_count(), 0.0);

  const auto active = [&](mesh::LocalCellId cell) {
    return domain == nullptr ||
           domain->region(cell) == immersed::CellRegion::fluid;
  };
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
    if (!active(cell))
      continue;
    const auto point = geometry.cell_center_m(cell);
    const double rho = density(point);
    const auto u = velocity(point);
    result.density[cell] = rho;
    result.velocity[cell * 3U] = u.x;
    result.velocity[cell * 3U + 1U] = u.y;
    result.velocity[cell * 3U + 2U] = u.z;
    result.transported_cell_fields[0][cell] = rho;
    result.transported_cell_fields[1][cell] = 0.2 * rho;
  }

  constexpr runtime::Int3 first{0, 0, 0};
  constexpr runtime::Int3 second{1, 0, 0};
  std::uint64_t selected_wall_face =
      std::numeric_limits<std::uint64_t>::max();
  if (domain != nullptr) {
    std::uint64_t local_wall_face =
        std::numeric_limits<std::uint64_t>::max();
    const auto linked = [&](mesh::LocalCellId cell) {
      const auto global = topology.global_cell_id(cell);
      return std::any_of(domain->links().begin(), domain->links().end(),
                         [global](const auto &link) {
                           return link.fluid_cell == global;
                         });
    };
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      const auto neighbour = topology.neighbour(face);
      if (!neighbour.has_value() ||
          topology.face_ownership(face) != mesh::EntityOwnership::owned ||
          !active(topology.owner(face)) || !active(*neighbour) ||
          (!linked(topology.owner(face)) && !linked(*neighbour)))
        continue;
      local_wall_face =
          std::min(local_wall_face,
                   static_cast<std::uint64_t>(topology.global_face_id(face)));
    }
    check_mpi(MPI_Allreduce(&local_wall_face, &selected_wall_face, 1U,
                            MPI_UINT64_T, MPI_MIN, mpi.comm()),
              "MPI_Allreduce(material-WALE wall face)");
    HUNDUN_CHECK(selected_wall_face !=
                 std::numeric_limits<std::uint64_t>::max());
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count(); ++face) {
    const auto owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    const bool owner_active = active(owner);
    const bool neighbour_active = neighbour.has_value() && active(*neighbour);
    if (!owner_active || (neighbour.has_value() && !neighbour_active))
      continue;
    const auto u = velocity(geometry.face_center_m(face));
    result.face_velocity[face * 3U] = u.x;
    result.face_velocity[face * 3U + 1U] = u.y;
    result.face_velocity[face * 3U + 2U] = u.z;
    result.face_mass_flux[face] =
        density(geometry.face_center_m(face)) *
        dot(u, geometry.face_area_vector_m2(face, mesh::FaceSide::owner));
    const bool body_fitted_seed =
        domain == nullptr && neighbour.has_value() &&
        !topology.periodic_pair(face).has_value() &&
        topology.cell_ownership(owner) == mesh::EntityOwnership::owned &&
        same(topology.global_cell(owner), first) &&
        same(topology.global_cell(*neighbour), second);
    const bool immersed_seed =
        domain != nullptr &&
        static_cast<std::uint64_t>(topology.global_face_id(face)) ==
            selected_wall_face;
    if (body_fitted_seed || immersed_seed)
      result.face_mass_flux[face] += seeded_mass_flux;
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

double maximum_difference(const std::vector<double> &left,
                          const std::vector<double> &right,
                          const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(left.size() == right.size());
  double result = 0.0;
  for (std::size_t index = 0U; index < left.size(); ++index)
    result = std::max(result, std::abs(left[index] - right[index]));
  mpi.allreduce_fp64_in_place(&result, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  return result;
}

void require_collective_any(bool local, const runtime::MpiContext &mpi) {
  double count = local ? 1.0 : 0.0;
  mpi.allreduce_fp64_in_place(&count, 1U,
                              runtime::Fp64ReductionOperation::sum);
  HUNDUN_CHECK(count > 0.0);
}

bool run_scenario(
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const immersed::ImmersedDomain *domain,
    const immersed::GhostStencilPlan *ghost_plan,
    const immersed::WallQuadraturePlan *wall_plan,
    const immersed::LocalFlowPatternTransform *transform,
    const runtime::MpiContext &mpi) {
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
  const auto rho_h = registry.declare_field(cell_field("rho_h", "J/m3"));
  const auto rho_alpha =
      registry.declare_field(cell_field("rho_alpha", "kg/m3"));
  fields.transported_cell_fields = {rho_h, rho_alpha};
  registry.freeze();
  const flow::MaterialDensityTransportSpec material_spec{
      rho_h, 0.0, {rho_alpha}, {0.0}};
  const auto make_state = [&](double seeded_mass_flux = 1.0e-4) {
    auto state = flow::FlowState::create(
        registry,
        {decomposition.local_extent(), topology.local_face_count()}, fields,
        {0U, 0.0, kDt, 0.0, flow::MomentumTimeOrder::backward_euler});
    const auto initial = initial_values(topology, geometry, domain,
                                        seeded_mass_flux, mpi);
    state.seed_accepted_layers(initial, initial);
    return state;
  };

  std::vector<mesh::GlobalCellId> active;
  std::size_t owned_active_count{};
  if (domain != nullptr) {
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
                         domain == nullptr ? 2 : 4));
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution), my(execution), mz(execution),
      pressure_preconditioner(execution);
  const flow::ImmersedFlowDensitySetup density_setup{
      config::DensityModel::material, &registry, fields, material_spec,
      nullptr};
  auto facade = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, domain, ghost_plan,
      wall_plan, transform, &wale, density_setup, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner);
  const flow::ImmersedFlowPhysics physics{
      config::DensityModel::material, 1.0, kMu, std::nullopt, std::nullopt,
      std::nullopt};
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, kDt, 0.0);
  linear::SolveControl control;
  control.max_iterations = 5000U;
  control.residual_recompute_interval = 20U;

  auto state = make_state();
  const auto accepted = facade.attempt(state, physics, stencil, control, control);
  const auto *material =
      std::get_if<flow::MaterialDensityStepAttemptReport>(&accepted.base);
  const bool supported =
      material != nullptr &&
      material->flow().disposition == flow::StepAttemptDisposition::committed &&
      accepted.wale.has_value() &&
      (domain == nullptr || accepted.force.has_value());
  if (!supported) {
    if (mpi.rank() == 0) {
      const auto reason = material == nullptr
                              ? static_cast<int>(
                                    std::get<flow::StepAttemptReport>(accepted.base)
                                        .reason)
                              : static_cast<int>(material->flow().reason);
      std::cerr << (domain == nullptr ? "BODY_MATERIAL_WALE_UNSUPPORTED"
                                      : "IBM_MATERIAL_WALE_UNSUPPORTED")
                << " reason=" << reason;
      if (material != nullptr) {
        std::cerr << " correctors=" << material->flow().pressure_corrector_count
                  << " momentum_residual="
                  << material->flow().final_momentum_normalized_l2[0] << ','
                  << material->flow().final_momentum_normalized_l2[1] << ','
                  << material->flow().final_momentum_normalized_l2[2]
                  << " momentum_defect="
                  << material->flow()
                         .final_momentum_relative_conservation_defect[0]
                  << ','
                  << material->flow()
                         .final_momentum_relative_conservation_defect[1]
                  << ','
                  << material->flow()
                         .final_momentum_relative_conservation_defect[2];
      }
      std::cerr << '\n';
    }
    return false;
  }

  HUNDUN_CHECK(material->flow().pressure_corrector_count == 2U);
  HUNDUN_CHECK(material->material_report_available());
  HUNDUN_CHECK(material->material_report().flux_provenance() ==
               flow::MaterialFluxProvenance::final_corrected);
  HUNDUN_CHECK(accepted.wale->identity.value != 0U);
  HUNDUN_CHECK(accepted.wale->owned_active_count == owned_active_count);
  double maximum_nu_t = accepted.wale->maximum_nu_t_m2_per_s;
  mpi.allreduce_fp64_in_place(&maximum_nu_t, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  HUNDUN_CHECK(maximum_nu_t > 0.0);

  HUNDUN_CHECK(
      flow::test::ImmersedFlowTestAccess::wale_evaluation_count(facade) == 1U);
  HUNDUN_CHECK(
      flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(facade) ==
      accepted.wale->identity);
  const auto accepted_values =
      state.snapshot(flow::FlowLayer::committed);
  const std::uint64_t accepted_wall_viscosity =
      flow::test::ImmersedFlowTestAccess::
          wall_effective_viscosity_fingerprint(facade);

  auto flux_variant_state = make_state(2.0e-4);
  const auto flux_variant_input =
      flux_variant_state.snapshot(flow::FlowLayer::committed);
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      state.snapshot(flow::FlowLayer::history),
      [&] {
        auto values = flux_variant_input;
        values.face_mass_flux = state.snapshot(flow::FlowLayer::history)
                                    .face_mass_flux;
        return values;
      }()));
  HUNDUN_CHECK(maximum_difference(
                   state.snapshot(flow::FlowLayer::history).face_mass_flux,
                   flux_variant_input.face_mass_flux, mpi) > 0.0);
  const auto flux_variant = facade.attempt(
      flux_variant_state, physics, stencil, control, control);
  const auto *flux_variant_material =
      std::get_if<flow::MaterialDensityStepAttemptReport>(&flux_variant.base);
  HUNDUN_CHECK(flux_variant_material != nullptr);
  HUNDUN_CHECK(flux_variant_material->flow().disposition ==
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(flux_variant.wale.has_value());
  require_collective_any(flux_variant.wale->identity != accepted.wale->identity,
                         mpi);
  if (domain != nullptr) {
    require_collective_any(
        flow::test::ImmersedFlowTestAccess::
                wall_effective_viscosity_fingerprint(facade) !=
            accepted_wall_viscosity,
        mpi);
  }

  auto viscosity_variant_state = make_state();
  auto viscosity_variant_physics = physics;
  viscosity_variant_physics.dynamic_viscosity_pa_s = 2.0 * kMu;
  const auto viscosity_variant = facade.attempt(
      viscosity_variant_state, viscosity_variant_physics, stencil, control,
      control);
  const auto *viscosity_variant_material =
      std::get_if<flow::MaterialDensityStepAttemptReport>(
          &viscosity_variant.base);
  HUNDUN_CHECK(viscosity_variant_material != nullptr);
  HUNDUN_CHECK(viscosity_variant_material->flow().disposition ==
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(viscosity_variant.wale.has_value());
  HUNDUN_CHECK(viscosity_variant.wale->identity == accepted.wale->identity);
  HUNDUN_CHECK(maximum_difference(
                   accepted_values.velocity,
                   viscosity_variant_state
                       .snapshot(flow::FlowLayer::committed)
                       .velocity,
                   mpi) > 1.0e-14);

  if (domain != nullptr) {
    HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::
                     wall_effective_viscosity_fingerprint(facade) != 0U);
  }

  auto retry_state = make_state();
  const auto before = capture(retry_state);
  flow::ImmersedFlowStepAttemptReport failed;
  if (domain != nullptr) {
    flow::test::ImmersedFlowTestAccess::set_failure_stage(
        flow::test::ImmersedFlowAttemptFailureStage::after_corrector_1);
    failed = facade.attempt(retry_state, physics, stencil, control, control);
    flow::test::ImmersedFlowTestAccess::clear_failure_stage();
  } else {
    auto failing_control = control;
    failing_control.max_iterations = 0U;
    failed = facade.attempt(retry_state, physics, stencil, failing_control,
                            control);
  }
  const auto *failed_material =
      std::get_if<flow::MaterialDensityStepAttemptReport>(&failed.base);
  HUNDUN_CHECK(failed_material != nullptr);
  HUNDUN_CHECK(failed_material->flow().disposition ==
               flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(!failed.force.has_value());
  HUNDUN_CHECK(!failed.wale.has_value());
  require_exact(before, retry_state);
  HUNDUN_CHECK(
      flow::test::ImmersedFlowTestAccess::wale_evaluation_count(facade) == 0U);

  const auto half_stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, 0.5 * kDt, 0.0);
  const auto retried = facade.attempt(retry_state, physics, half_stencil,
                                      control, control);
  const auto *retried_material =
      std::get_if<flow::MaterialDensityStepAttemptReport>(&retried.base);
  HUNDUN_CHECK(retried_material != nullptr);
  HUNDUN_CHECK(retried_material->flow().disposition ==
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(domain == nullptr ? !retried.force.has_value()
                                 : retried.force.has_value());
  HUNDUN_CHECK(retried.wale.has_value());
  HUNDUN_CHECK(retried.wale->identity != accepted.wale->identity);
  HUNDUN_CHECK(
      flow::test::ImmersedFlowTestAccess::wale_evaluation_count(facade) == 1U);
  HUNDUN_CHECK(
      flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(facade) ==
      retried.wale->identity);
  if (formal_scientific && domain != nullptr) {
    double nu_t_square = accepted.wale->l2_nu_t_m2_per_s *
                         accepted.wale->l2_nu_t_m2_per_s *
                         geometry.cell_volume_m3(0U);
    mpi.allreduce_fp64_in_place(&nu_t_square, 1U,
                                runtime::Fp64ReductionOperation::sum);
    const double nu_t_l2 = std::sqrt(nu_t_square);
    const auto &flow_report = material->flow();
    double conservation = flow_report.final_mass_relative_conservation_defect;
    for (const double value :
         flow_report.final_momentum_relative_conservation_defect)
      conservation = std::max(conservation, value);
    for (const double value :
         flow_report.final_transport_relative_conservation_defect)
      conservation = std::max(conservation, value);
    const auto &material_report = material->material_report();
    conservation = std::max(
        conservation, material_report.mass_relative_conservation_defect());
    for (const double value :
         material_report.transport_relative_conservation_defect())
      conservation = std::max(conservation, value);
    double closure = 0.0;
    for (const double value : material_report.transport_normalized_l2())
      closure = std::max(closure, value);
    const auto consistency = accepted.force->consistency.total_N;
    const double force = std::sqrt(consistency.x * consistency.x +
                                   consistency.y * consistency.y +
                                   consistency.z * consistency.z);
    const test::stage3::ScientificRow row{
        "material-ibm-wale-n" + std::to_string(selected_cells) + "-r" +
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
        {true, flow_report.final_continuity_normalized_l2},
        {true, conservation},
        {true, closure},
        {true, force},
        {true, accepted.wale->identity.value},
        "pass"};
    HUNDUN_CHECK(test::stage3::validate_scientific_row(row));
    if (mpi.rank() == 0)
      std::cout << test::stage3::serialize_scientific_row(row) << '\n';
  }
  return true;
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
  const auto boundaries =
      boundary::BoundaryRegistry::create(periodic_case(mpi.size()), topology);

  const bool body_fitted =
      formal_scientific ||
      run_scenario(decomposition, topology, geometry, boundaries, nullptr,
                   nullptr, nullptr, nullptr, mpi);

  std::optional<test::Stage3TemporaryDirectory> surface_directory;
  const auto surface_path = collective_surface_path(mpi, surface_directory);
  const auto surface = immersed::ImmersedSurface::load_collective(
      std::filesystem::path(surface_path), 1.0, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  const auto ghost_plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  const auto wall_plan = immersed::WallQuadraturePlan::create(
      surface, query, domain, topology, geometry, mpi);
  immersed::LocalFlowPatternTransform transform;
  HUNDUN_CHECK(ghost_plan.maximum_halo_reach() <= 4U);
  const bool immersed_flow =
      run_scenario(decomposition, topology, geometry, boundaries, &domain,
                   &ghost_plan, &wall_plan, &transform, mpi);

  HUNDUN_CHECK(body_fitted);
  HUNDUN_CHECK(immersed_flow);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  if (argc == 2 && std::string(argv[1]) == "fast") {
    return hundun::test::run(run);
  }
  if (argc == 3 && std::string(argv[1]) == "formal" &&
      (std::string(argv[2]) == "12" || std::string(argv[2]) == "24")) {
    selected_cells = std::string(argv[2]) == "12" ? 12 : 24;
    formal_scientific = true;
    return hundun::test::run(run);
  }
  return 2;
}
