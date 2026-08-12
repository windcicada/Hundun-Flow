// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/stage3_mms.hpp"
#include "tests/support/stage3_scientific_row.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/diag_immersed_module.hpp"
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

namespace {

using namespace hundun;

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kRho = 1.0;
constexpr double kMolecularMu = 0.01;
constexpr double kDt = 1.0e-4;
bool formal_scientific{};

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("unsupported immersed-WALE rank count");
}

config::FlowCaseConfig periodic_case(int cells, int ranks) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = "stage3-c1-immersed-wale-constant";
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::constant;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = process_grid(ranks);
  result.mesh.cells = {cells, cells, cells};
  result.mesh.origin_m = {0.0, 0.0, 0.0};
  result.mesh.length_m = {1.0, 1.0, 1.0};
  result.physics.rho_ref_kg_per_m3 = kRho;
  result.physics.dynamic_viscosity_pa_s = kMolecularMu;
  result.physics.inlet_consistency_rtol = 1.0e-12;
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

runtime::FieldDescriptor cell_field(const char *name,
                                    std::uint32_t components) {
  return {name,
          "1",
          "stage3_c1_immersed_wale_constant",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          4,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_field(const char *name,
                                    std::uint32_t components) {
  return {name,
          "1",
          "stage3_c1_immersed_wale_constant",
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
    std::optional<test::Stage3TemporaryDirectory> &directory, int cells) {
  std::string path_text;
  if (mpi.rank() == 0) {
    directory.emplace("stage3-c1-immersed-wale");
    const auto body = test::stage3::approved_body(
        test::stage3::BodyKind::sphere);
    const auto surface = test::stage3::make_manufactured_surface(
        body, 1.0 / static_cast<double>(cells));
    const auto path = directory->path() / "sphere.stl";
    test::write_text(path, test::ascii_stl(surface.triangles, "sphere"));
    path_text = path.string();
  }
  std::uint64_t size = path_text.size();
  check_mpi(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()),
            "MPI_Bcast(immersed-WALE path size)");
  HUNDUN_CHECK(size <=
               static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
  path_text.resize(static_cast<std::size_t>(size));
  check_mpi(MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()),
                      MPI_BYTE, 0, mpi.comm()),
            "MPI_Bcast(immersed-WALE path)");
  return path_text;
}

runtime::Real3 velocity(runtime::Real3 point) {
  constexpr double amplitude = 0.05;
  const double x = 2.0 * kPi * point.x;
  const double y = 2.0 * kPi * point.y;
  const double z = 2.0 * kPi * point.z;
  return {amplitude * std::sin(x) * std::cos(y) * std::cos(z),
          -amplitude * std::cos(x) * std::sin(y) * std::cos(z), 0.0};
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

flow::FlowLayerValues initial_values(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::ImmersedDomain &domain) {
  flow::FlowLayerValues result;
  result.density.assign(topology.owned_cell_count(), 0.0);
  result.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  result.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  result.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  result.face_mass_flux.assign(topology.local_face_count(), 0.0);
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const auto value = velocity(geometry.cell_center_m(cell));
    result.density[cell] = kRho;
    result.velocity[cell * 3U] = value.x;
    result.velocity[cell * 3U + 1U] = value.y;
    result.velocity[cell * 3U + 2U] = value.z;
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const bool owner_active =
        domain.region(topology.owner(face)) == immersed::CellRegion::fluid;
    const auto neighbour = topology.neighbour(face);
    const bool neighbour_active =
        neighbour.has_value() &&
        domain.region(*neighbour) == immersed::CellRegion::fluid;
    if (!owner_active || (neighbour.has_value() && !neighbour_active))
      continue;
    const auto value = velocity(geometry.face_center_m(face));
    result.face_velocity[face * 3U] = value.x;
    result.face_velocity[face * 3U + 1U] = value.y;
    result.face_velocity[face * 3U + 2U] = value.z;
    result.face_mass_flux[face] = dot(
        value, geometry.face_area_vector_m2(face, mesh::FaceSide::owner));
  }
  return result;
}

template <class Mutation>
void require_report_mutation_rejected(
    const flow::FixedStepImmersedFlow &facade, const flow::FlowState &state,
    const flow::ImmersedFlowStepAttemptReport &accepted,
    Mutation mutate) {
  auto changed = accepted;
  mutate(changed);
  bool rejected = false;
  try {
    static_cast<void>(facade.diagnostic_source(state, changed));
  } catch (const runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

template <class Mutation>
void require_diagnostic_mutation_rejected(
    const flow::FixedStepImmersedFlow &facade, const flow::FlowState &state,
    const flow::ImmersedFlowStepAttemptReport &accepted,
    Mutation mutate) {
  require_report_mutation_rejected(
      facade, state, accepted, [&](auto &changed) { mutate(*changed.wale); });
}

void run_case(int cells) {
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
  if (formal_scientific)
    HUNDUN_CHECK((cells == 24 &&
                  (mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4)) ||
                 (cells == 48 && mpi.size() == 1));
  const auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {cells, cells, cells}, {true, true, true},
      runtime::DecompositionOptions{process_grid(mpi.size())});
  const mesh::MeshTopology topology(decomposition);
  const mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  const auto boundaries = boundary::BoundaryRegistry::create(
      periodic_case(cells, mpi.size()), topology);

  std::optional<test::Stage3TemporaryDirectory> surface_directory;
  const auto surface_path =
      collective_surface_path(mpi, surface_directory, cells);
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
  HUNDUN_CHECK(ghost_plan.maximum_halo_reach() <= 4U);

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_field("rho", 1U));
  fields.velocity = registry.declare_field(cell_field("velocity", 3U));
  fields.mechanical_pressure =
      registry.declare_field(cell_field("pi", 1U));
  fields.face_velocity =
      registry.declare_field(face_field("face_velocity", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  registry.freeze();
  const auto make_state = [&] {
    auto state = flow::FlowState::create(
        registry,
        {decomposition.local_extent(), topology.local_face_count()}, fields,
        {0U, 0.0, kDt, 0.0, flow::MomentumTimeOrder::backward_euler});
    const auto initial = initial_values(topology, geometry, domain);
    state.seed_accepted_layers(initial, initial);
    return state;
  };

  execution::CpuReferenceContext execution;
  auto wale = les::WaleModel::create(
      {0.5, 0.9, 0.7}, topology, geometry,
      domain.active_cells().owned_active_count(),
      domain.active_cells().ordered_global_ids(), execution);
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 4));
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution), my(execution), mz(execution),
      pressure_preconditioner(execution);
  immersed::LocalFlowPatternTransform transform;
  auto facade = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, &domain, &ghost_plan,
      &wall_plan, &transform, &wale, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_preconditioner);
  const flow::ImmersedFlowPhysics physics{
      config::DensityModel::constant, kRho, kMolecularMu, std::nullopt,
      std::nullopt, std::nullopt};
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, kDt, 0.0);
  linear::SolveControl control;
  control.max_iterations = 5000U;
  control.residual_recompute_interval = 20U;

  auto state = make_state();
  const auto before_history = state.snapshot(flow::FlowLayer::history);
  const auto before_committed = state.snapshot(flow::FlowLayer::committed);
  const auto before_metadata = state.metadata();
  flow::test::ImmersedFlowTestAccess::set_failure_stage(
      flow::test::ImmersedFlowAttemptFailureStage::after_corrector_1);
  const auto failed = facade.attempt(state, physics, stencil, control, control);
  flow::test::ImmersedFlowTestAccess::clear_failure_stage();
  const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
  if (failed_base.disposition !=
          flow::StepAttemptDisposition::recoverable_failure &&
      mpi.rank() == 0) {
    std::cerr << "COMBINED_IBM_WALE_UNSUPPORTED reason="
              << static_cast<int>(failed_base.reason) << '\n';
  }
  HUNDUN_CHECK(failed_base.disposition ==
               flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(!failed.wale.has_value());
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      before_history, state.snapshot(flow::FlowLayer::history)));
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      before_committed, state.snapshot(flow::FlowLayer::committed)));
  HUNDUN_CHECK(test::accepted_step_metadata_bitwise_equal(before_metadata,
                                                          state.metadata()));
  HUNDUN_CHECK(
      flow::test::ImmersedFlowTestAccess::wale_evaluation_count(facade) == 0U);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::
                   wale_coefficient_identity(facade)
                   .value == 0U);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::
                   wall_effective_viscosity_fingerprint(facade) == 0U);

  const auto accepted = facade.attempt(state, physics, stencil, control, control);
  const auto &base = std::get<flow::StepAttemptReport>(accepted.base);
  if (base.disposition != flow::StepAttemptDisposition::committed &&
      mpi.rank() == 0) {
    std::cerr << "IMMERSED_WALE_CONSTANT_FAILURE reason="
              << static_cast<int>(base.reason)
              << " rank=" << base.lowest_failing_rank
              << " correctors=" << base.pressure_corrector_count << '\n';
  }
  HUNDUN_CHECK(base.disposition == flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(accepted.wale.has_value());
  HUNDUN_CHECK(base.pressure_corrector_count == 2U);
  HUNDUN_CHECK(
      flow::test::ImmersedFlowTestAccess::wale_evaluation_count(facade) == 1U);
  HUNDUN_CHECK(
      flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(facade) ==
      accepted.wale->identity);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::
                   wall_effective_viscosity_fingerprint(facade) != 0U);
  HUNDUN_CHECK(accepted.force.has_value());
  HUNDUN_CHECK(accepted.wale->identity.value != 0U);
  HUNDUN_CHECK(accepted.wale->owned_active_count ==
               domain.active_cells().owned_active_count());
  double maximum_nu_t = accepted.wale->maximum_nu_t_m2_per_s;
  mpi.allreduce_fp64_in_place(&maximum_nu_t, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  HUNDUN_CHECK(maximum_nu_t > 0.0);
  static_cast<void>(facade.diagnostic_source(state, accepted));
  require_diagnostic_mutation_rejected(
      facade, state, accepted,
      [](les::WaleSummary &summary) { ++summary.identity.value; });
  require_diagnostic_mutation_rejected(
      facade, state, accepted,
      [](les::WaleSummary &summary) { summary.minimum_nu_t_m2_per_s += 1.0; });
  require_diagnostic_mutation_rejected(
      facade, state, accepted,
      [](les::WaleSummary &summary) { summary.maximum_nu_t_m2_per_s += 1.0; });
  require_diagnostic_mutation_rejected(
      facade, state, accepted,
      [](les::WaleSummary &summary) { summary.l2_nu_t_m2_per_s += 1.0; });
  require_diagnostic_mutation_rejected(
      facade, state, accepted,
      [](les::WaleSummary &summary) { ++summary.exact_zero_count; });
  require_diagnostic_mutation_rejected(
      facade, state, accepted,
      [](les::WaleSummary &summary) { ++summary.owned_active_count; });
  require_report_mutation_rejected(
      facade, state, accepted, [](auto &changed) {
        std::get<flow::StepAttemptReport>(changed.base).attempted_dt_s += 1.0;
      });
  require_report_mutation_rejected(
      facade, state, accepted, [](auto &changed) {
        changed.force->operator_force.total_N.x += 1.0;
      });
  require_report_mutation_rejected(
      facade, state, accepted, [](auto &changed) {
        changed.force->budget_reaction.total_N.x += 1.0;
      });
  require_report_mutation_rejected(
      facade, state, accepted, [](auto &changed) {
        changed.force->surface_traction.total_N.x += 1.0;
      });
  require_report_mutation_rejected(
      facade, state, accepted, [](auto &changed) {
        changed.force->consistency.total_N.x += 1.0;
      });

  if (cells == 8) {
    const auto bdf2_stencil = flow::make_momentum_time_stencil(
        flow::MomentumTimeOrder::bdf2, kDt, kDt);
    const auto bdf2 =
        facade.attempt(state, physics, bdf2_stencil, control, control);
    const auto &bdf2_base = std::get<flow::StepAttemptReport>(bdf2.base);
    HUNDUN_CHECK(bdf2_base.disposition ==
                 flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(bdf2_base.pressure_corrector_count == 2U);
    HUNDUN_CHECK(bdf2.wale.has_value());
    HUNDUN_CHECK(bdf2.force.has_value());
    HUNDUN_CHECK(
        flow::test::ImmersedFlowTestAccess::wale_evaluation_count(facade) ==
        1U);
    HUNDUN_CHECK(
        flow::test::ImmersedFlowTestAccess::wale_coefficient_identity(facade) ==
        bdf2.wale->identity);
    HUNDUN_CHECK(bdf2.wale->identity != accepted.wale->identity);
    HUNDUN_CHECK(state.metadata().step == 2U);
    HUNDUN_CHECK(state.metadata().order == flow::MomentumTimeOrder::bdf2);
  }
  if (formal_scientific) {
    double nu_t_square = accepted.wale->l2_nu_t_m2_per_s *
                         accepted.wale->l2_nu_t_m2_per_s *
                         geometry.cell_volume_m3(0U);
    mpi.allreduce_fp64_in_place(&nu_t_square, 1U,
                                runtime::Fp64ReductionOperation::sum);
    const double nu_t_l2 = std::sqrt(nu_t_square);
    double conservation = base.final_mass_relative_conservation_defect;
    for (const double value : base.final_momentum_relative_conservation_defect)
      conservation = std::max(conservation, value);
    const auto consistency = accepted.force->consistency.total_N;
    const double force = std::sqrt(consistency.x * consistency.x +
                                   consistency.y * consistency.y +
                                   consistency.z * consistency.z);
    const test::stage3::ScientificRow row{
        "constant-ibm-wale-n" + std::to_string(cells) + "-r" +
            std::to_string(mpi.size()),
        {cells, cells, cells},
        mpi.size(),
        process_grid(mpi.size()),
        1U,
        kDt,
        kDt,
        {false, 0.0},
        {false, 0.0},
        {true, nu_t_l2},
        {true, base.final_continuity_normalized_l2},
        {true, conservation},
        {false, 0.0},
        {true, force},
        {true, accepted.wale->identity.value},
        "pass"};
    HUNDUN_CHECK(test::stage3::validate_scientific_row(row));
    if (mpi.rank() == 0)
      std::cout << test::stage3::serialize_scientific_row(row) << '\n';
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  int cells = 8;
  if (argc == 2 && std::string(argv[1]) == "fast") {
    cells = 12;
  } else if (argc == 3 && std::string(argv[1]) == "formal" &&
             (std::string(argv[2]) == "24" || std::string(argv[2]) == "48")) {
    cells = std::string(argv[2]) == "24" ? 24 : 48;
    formal_scientific = true;
  } else if (argc != 1) {
    return 2;
  }
  return hundun::test::run([&] { run_case(cells); });
}
