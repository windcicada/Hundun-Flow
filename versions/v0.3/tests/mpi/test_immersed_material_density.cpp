// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/flow_material_density_piso_test_access.hpp"
#include "tests/support/flow_material_density_transport_test_access.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/stage3_mms.hpp"
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
#include "hundun/lin_bicgstab.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_access_plan.hpp"
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

constexpr int kCells = 8;
constexpr double kDt = 1.0e-4;
constexpr double kMu = 0.01;

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  throw runtime::Error("unsupported immersed-material rank count");
}

config::FlowCaseConfig periodic_case(int ranks) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = "stage3-d1-immersed-material";
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::material;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = process_grid(ranks);
  result.mesh.cells = {kCells, kCells, kCells};
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
          "stage3_d1_immersed_material",
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
          "stage3_d1_immersed_material",
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
    directory.emplace("stage3-d1-immersed-material");
    const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
    const auto surface = test::stage3::make_manufactured_surface(
        body, 1.0 / static_cast<double>(kCells));
    const auto path = directory->path() / "sphere.stl";
    test::write_text(path, test::ascii_stl(surface.triangles, "sphere"));
    path_text = path.string();
  }
  std::uint64_t size = path_text.size();
  check_mpi(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()),
            "MPI_Bcast(immersed-material path size)");
  HUNDUN_CHECK(size <=
               static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
  path_text.resize(static_cast<std::size_t>(size));
  check_mpi(MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()),
                      MPI_BYTE, 0, mpi.comm()),
            "MPI_Bcast(immersed-material path)");
  return path_text;
}

double density(runtime::Real3 point) noexcept {
  return 1.0 + 0.10 * point.x + 0.05 * point.y - 0.03 * point.z;
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
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

flow::FlowLayerValues initial_values(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::ImmersedDomain &domain) {
  flow::FlowLayerValues result;
  result.density.assign(topology.owned_cell_count(), 0.0);
  result.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  result.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  result.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  result.face_mass_flux.assign(topology.local_face_count(), 0.0);
  result.transported_cell_fields.resize(2U);
  for (auto &field : result.transported_cell_fields)
    field.assign(topology.owned_cell_count(), 0.0);
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const double rho = density(geometry.cell_center_m(cell));
    result.density[cell] = rho;
    result.transported_cell_fields[0][cell] = rho;
    result.transported_cell_fields[1][cell] = 0.2 * rho;
  }
  constexpr runtime::Int3 first{0, 0, 0};
  constexpr runtime::Int3 second{1, 0, 0};
  constexpr double seeded_mass_flux = 1.0e-4;
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto neighbour = topology.neighbour(face);
    if (!neighbour.has_value() || topology.periodic_pair(face).has_value() ||
        topology.cell_ownership(topology.owner(face)) !=
            mesh::EntityOwnership::owned ||
        domain.region(topology.owner(face)) != immersed::CellRegion::fluid ||
        domain.region(*neighbour) != immersed::CellRegion::fluid)
      continue;
    const auto owner = topology.global_cell(topology.owner(face));
    const auto other = topology.global_cell(*neighbour);
    const auto same = [](runtime::Int3 left, runtime::Int3 right) {
      return left.x == right.x && left.y == right.y && left.z == right.z;
    };
    if (same(owner, first) && same(other, second))
      result.face_mass_flux[face] = seeded_mass_flux;
    else if (same(owner, second) && same(other, first))
      result.face_mass_flux[face] = -seeded_mass_flux;
  }
  return result;
}

template <class Function>
flow::MaterialDensityStepAttemptReport require_failed_attempt_neutral(
    flow::FixedStepImmersedFlow &facade, flow::FlowState &state,
    const flow::ImmersedFlowPhysics &physics,
    const flow::MomentumTimeStencil &stencil,
    const linear::SolveControl &control, Function inject) {
  const auto before = capture(state);
  inject(true);
  const auto failed = facade.attempt(state, physics, stencil, control, control);
  inject(false);
  HUNDUN_CHECK(std::holds_alternative<flow::MaterialDensityStepAttemptReport>(
      failed.base));
  const auto &material =
      std::get<flow::MaterialDensityStepAttemptReport>(failed.base);
  HUNDUN_CHECK(material.flow().disposition !=
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(flow::test::MaterialDensityPisoTestAccess::report_authenticated(
      material));
  HUNDUN_CHECK(!failed.force.has_value());
  HUNDUN_CHECK(!failed.wale.has_value());
  require_exact(before, state);
  return material;
}

void run_case(bool transaction_mode) {
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2);
  const auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {kCells, kCells, kCells}, {true, true, true},
      runtime::DecompositionOptions{process_grid(mpi.size())});
  const mesh::MeshTopology topology(decomposition);
  const mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  const auto boundaries =
      boundary::BoundaryRegistry::create(periodic_case(mpi.size()), topology);

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
  HUNDUN_CHECK(ghost_plan.maximum_halo_reach() <= 4U);

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_field("rho", "kg/m3"));
  fields.velocity = registry.declare_field(
      cell_field("velocity", "m/s", 3U, false));
  fields.mechanical_pressure =
      registry.declare_field(cell_field("pi", "Pa", 1U, false));
  fields.face_velocity =
      registry.declare_field(face_field("face_velocity", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  const auto rho_h =
      registry.declare_field(cell_field("rho_h", "J/m3"));
  const auto rho_alpha =
      registry.declare_field(cell_field("rho_alpha", "kg/m3"));
  fields.transported_cell_fields = {rho_h, rho_alpha};
  registry.freeze();
  const flow::MaterialDensityTransportSpec material_spec{
      rho_h, 0.0, {rho_alpha}, {0.0}};
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
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 4));
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution), my(execution), mz(execution),
      pressure_preconditioner(execution);
  immersed::LocalFlowPatternTransform transform;
  const flow::ImmersedFlowDensitySetup density_setup{
      config::DensityModel::material, &registry, fields, material_spec,
      nullptr};
  auto facade = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, &domain, &ghost_plan,
      &wall_plan, &transform, nullptr, density_setup, mpi, execution, halo,
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
  const auto initial = state.snapshot(flow::FlowLayer::committed);
  std::uint64_t seeded_faces = 0U;
  for (const double value : initial.face_mass_flux)
    seeded_faces += value != 0.0 ? 1U : 0U;
  std::uint64_t global_seeded_faces{};
  check_mpi(MPI_Allreduce(&seeded_faces, &global_seeded_faces, 1U,
                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(immersed-material seeded faces)");
  HUNDUN_CHECK(global_seeded_faces == 1U);
  const auto accepted = facade.attempt(state, physics, stencil, control, control);
  if (!std::holds_alternative<flow::MaterialDensityStepAttemptReport>(
          accepted.base)) {
    if (mpi.rank() == 0)
      std::cerr << "IMMERSED_MATERIAL_UNSUPPORTED\n";
    HUNDUN_CHECK(false);
  }
  const auto &material =
      std::get<flow::MaterialDensityStepAttemptReport>(accepted.base);
  if (material.flow().disposition != flow::StepAttemptDisposition::committed &&
      mpi.rank() == 0) {
    std::cerr << "IMMERSED_MATERIAL_FAILURE reason="
              << static_cast<int>(material.flow().reason)
              << " material_reason="
              << static_cast<int>(material.material_failure_reason())
              << " rank=" << material.flow().lowest_failing_rank
              << " correctors=" << material.flow().pressure_corrector_count
              << '\n';
  }
  HUNDUN_CHECK(material.flow().disposition ==
               flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(material.flow().pressure_corrector_count == 2U);
  HUNDUN_CHECK(material.material_report_available());
  HUNDUN_CHECK(material.material_report().flux_provenance() ==
               flow::MaterialFluxProvenance::final_corrected);
  HUNDUN_CHECK(material.material_report().mass_conservation_available());
  HUNDUN_CHECK(
      material.material_report().mass_relative_conservation_defect() <=
      5.0e-12);
  for (const double defect :
       material.flow().final_momentum_relative_conservation_defect) {
    HUNDUN_CHECK(std::isfinite(defect));
    HUNDUN_CHECK(defect <= 5.0e-11);
  }
  HUNDUN_CHECK(flow::test::MaterialDensityPisoTestAccess::report_authenticated(
      material));
  HUNDUN_CHECK(accepted.force.has_value());
  HUNDUN_CHECK(!accepted.wale.has_value());

  auto wall_density_probe_state = make_state();
  const auto wall_density_probe =
      flow::test::ImmersedFlowTestAccess::wall_predictor_values(
          facade, wall_density_probe_state, 1.0, stencil, 1.0);
  double maximum_product_wall_density_error = 0.0;
  double maximum_product_wall_gradient_error = 0.0;
  std::uint64_t product_wall_density_rows = 0U;
  constexpr runtime::Real3 product_analytic_gradient{0.10, 0.05, -0.03};
  for (const auto &value : wall_density_probe) {
    const auto link = std::find_if(
        domain.links().begin(), domain.links().end(),
        [&](const auto &candidate) { return candidate.id == value.link; });
    HUNDUN_CHECK(link != domain.links().end());
    maximum_product_wall_density_error =
        std::max(maximum_product_wall_density_error,
                 std::abs(value.wall_density_kg_per_m3 -
                          density(link->wall_intercept_m)));
    maximum_product_wall_gradient_error =
        std::max(maximum_product_wall_gradient_error,
                 std::abs(value.wall_density_normal_derivative_kg_per_m4 -
                          dot(product_analytic_gradient,
                              link->solid_to_fluid_normal)));
    ++product_wall_density_rows;
  }
  double maximum_product_wall_errors[2]{
      maximum_product_wall_density_error,
      maximum_product_wall_gradient_error};
  mpi.allreduce_fp64_in_place(maximum_product_wall_errors, 2U,
                              runtime::Fp64ReductionOperation::maximum);
  double product_wall_density_row_count =
      static_cast<double>(product_wall_density_rows);
  mpi.allreduce_fp64_in_place(&product_wall_density_row_count, 1U,
                              runtime::Fp64ReductionOperation::sum);
  HUNDUN_CHECK(product_wall_density_row_count > 0.0);
  HUNDUN_CHECK(maximum_product_wall_errors[0] <= 1.0e-12);
  HUNDUN_CHECK(maximum_product_wall_errors[1] <= 1.0e-10);

  if (transaction_mode) {
    static_cast<void>(facade.diagnostic_source(state, accepted));
    constexpr auto last_corruption =
        flow::test::MaterialReportCorruptionForTest::
            nested_conservation_outer_size;
    for (std::uint8_t value = 0U;
         value <= static_cast<std::uint8_t>(last_corruption); ++value) {
      auto changed = accepted;
      auto &changed_material =
          std::get<flow::MaterialDensityStepAttemptReport>(changed.base);
      flow::test::MaterialDensityPisoTestAccess::corrupt_report(
          changed_material,
          static_cast<flow::test::MaterialReportCorruptionForTest>(value));
      HUNDUN_CHECK(
          !flow::test::MaterialDensityPisoTestAccess::report_authenticated(
              changed_material));
      bool rejected = false;
      try {
        static_cast<void>(facade.diagnostic_source(state, changed));
      } catch (const runtime::Error &) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
    }

    const int injected_rank = mpi.size() > 1 ? 1 : 0;
    auto negative_wall_state = make_state();
    const auto negative_wall_report = require_failed_attempt_neutral(
        facade, negative_wall_state, physics, stencil, control,
        [&](bool enabled) {
          if (enabled)
            flow::test::ImmersedFlowTestAccess::set_wall_input_failure(
                flow::test::ImmersedFlowWallInputFailure::non_positive_density,
                injected_rank);
          else
            flow::test::ImmersedFlowTestAccess::clear_wall_input_failure();
        });
    HUNDUN_CHECK(negative_wall_report.material_failure_reason() ==
                 flow::MaterialTransportFailureReason::non_positive_density);
    HUNDUN_CHECK(negative_wall_report.flow().lowest_failing_rank ==
                 injected_rank);

    auto non_finite_wall_state = make_state();
    const auto non_finite_wall_report = require_failed_attempt_neutral(
        facade, non_finite_wall_state, physics, stencil, control,
        [&](bool enabled) {
          if (enabled)
            flow::test::ImmersedFlowTestAccess::set_wall_input_failure(
                flow::test::ImmersedFlowWallInputFailure::non_finite_density,
                injected_rank);
          else
            flow::test::ImmersedFlowTestAccess::clear_wall_input_failure();
        });
    HUNDUN_CHECK(non_finite_wall_report.material_failure_reason() ==
                 flow::MaterialTransportFailureReason::non_finite_state);
    HUNDUN_CHECK(non_finite_wall_report.flow().lowest_failing_rank ==
                 injected_rank);

    auto nan_state = make_state();
    const int nan_rank = injected_rank;
    if (mpi.rank() == nan_rank) {
      std::optional<std::size_t> donor;
      for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
           ++cell) {
        if (domain.region(cell) == immersed::CellRegion::fluid) {
          donor = static_cast<std::size_t>(cell);
          break;
        }
      }
      HUNDUN_CHECK(donor.has_value());
      flow::test::MaterialDensityTransportTestAccess::
          set_accepted_transport_value(
              nan_state, false, 0U, *donor,
              std::numeric_limits<double>::quiet_NaN());
    }
    const auto nan_before = capture(nan_state);
    const auto nan_report =
        facade.attempt(nan_state, physics, stencil, control, control);
    HUNDUN_CHECK(std::holds_alternative<
                 flow::MaterialDensityStepAttemptReport>(nan_report.base));
    HUNDUN_CHECK(std::get<flow::MaterialDensityStepAttemptReport>(nan_report.base)
                     .flow()
                     .disposition != flow::StepAttemptDisposition::committed);
    const auto &nan_material =
        std::get<flow::MaterialDensityStepAttemptReport>(nan_report.base);
    HUNDUN_CHECK(nan_material.material_failure_reason() ==
                 flow::MaterialTransportFailureReason::non_finite_state);
    HUNDUN_CHECK(nan_material.flow().lowest_failing_rank == nan_rank);
    require_exact(nan_before, nan_state);

    auto rank_failure_state = make_state();
    const auto rank_failure_report = require_failed_attempt_neutral(
        facade, rank_failure_state, physics, stencil, control,
        [&](bool enabled) {
          if (enabled)
            flow::test::ImmersedFlowTestAccess::set_failure_stage(
                flow::test::ImmersedFlowAttemptFailureStage::after_corrector_1,
                injected_rank);
          else
            flow::test::ImmersedFlowTestAccess::clear_failure_stage();
        });
    HUNDUN_CHECK(rank_failure_report.material_failure_reason() ==
                 flow::MaterialTransportFailureReason::none);
    HUNDUN_CHECK(rank_failure_report.flow().lowest_failing_rank ==
                 injected_rank);

    for (const auto stage :
         {flow::test::ImmersedFlowAttemptFailureStage::after_final_transport,
          flow::test::ImmersedFlowAttemptFailureStage::before_commit}) {
      auto late_failure_state = make_state();
      const auto late_failure_report = require_failed_attempt_neutral(
          facade, late_failure_state, physics, stencil, control,
          [&](bool enabled) {
            if (enabled)
              flow::test::ImmersedFlowTestAccess::set_failure_stage(
                  stage, injected_rank);
            else
              flow::test::ImmersedFlowTestAccess::clear_failure_stage();
          });
      HUNDUN_CHECK(late_failure_report.flow().lowest_failing_rank ==
                   injected_rank);
      HUNDUN_CHECK(
          late_failure_report.material_failure_reason() ==
          (stage == flow::test::ImmersedFlowAttemptFailureStage::
                        after_final_transport
               ? flow::MaterialTransportFailureReason::non_finite_state
               : flow::MaterialTransportFailureReason::none));
    }

    return;
  }

  const auto committed = state.snapshot(flow::FlowLayer::committed);
  double initial_active_mass = 0.0;
  double final_active_mass = 0.0;
  double maximum_density_change = 0.0;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
    if (domain.region(cell) == immersed::CellRegion::fluid) {
      HUNDUN_CHECK(committed.density[cell] > 0.0);
      HUNDUN_CHECK(std::isfinite(committed.density[cell]));
      const double volume = geometry.cell_volume_m3(cell);
      initial_active_mass += volume * initial.density[cell];
      final_active_mass += volume * committed.density[cell];
      maximum_density_change =
          std::max(maximum_density_change,
                   std::abs(committed.density[cell] - initial.density[cell]));
    } else {
      HUNDUN_CHECK(test::fp64_bits(committed.density[cell]) ==
                   test::fp64_bits(0.0));
      for (const auto &field : committed.transported_cell_fields)
        HUNDUN_CHECK(test::fp64_bits(field[cell]) == test::fp64_bits(0.0));
    }
  }
  double masses[2]{initial_active_mass, final_active_mass};
  mpi.allreduce_fp64_in_place(masses, 2U,
                              runtime::Fp64ReductionOperation::sum);
  const double relative_mass_error =
      std::abs(masses[1] - masses[0]) / std::max(std::abs(masses[0]), 1.0);
  HUNDUN_CHECK(relative_mass_error <= 5.0e-12);
  mpi.allreduce_fp64_in_place(&maximum_density_change, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  HUNDUN_CHECK(maximum_density_change > 1.0e-8);

  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count(); ++face) {
    const auto neighbour = topology.neighbour(face);
    if (!neighbour.has_value())
      continue;
    const bool owner_active =
        domain.region(topology.owner(face)) == immersed::CellRegion::fluid;
    const bool neighbour_active =
        domain.region(*neighbour) == immersed::CellRegion::fluid;
    if (owner_active != neighbour_active)
      HUNDUN_CHECK(test::fp64_bits(committed.face_mass_flux[face]) ==
                   test::fp64_bits(0.0));
  }

  runtime::FieldAccessPlan density_access(registry);
  constexpr runtime::PhaseId phase = 301U;
  constexpr runtime::ActorId actor = 301U;
  density_access.declare_access(phase, actor, fields.density,
                                runtime::AccessMode::read);
  density_access.freeze();
  const auto density_view = state.layer(flow::FlowLayer::committed)
                                .acquire_read<double>(density_access, phase,
                                                      actor, fields.density);
  double maximum_wall_value_error = 0.0;
  double maximum_gradient_error = 0.0;
  double maximum_zero_mutation_error = 0.0;
  std::uint64_t reconstructed_rows = 0U;
  constexpr runtime::Real3 analytic_gradient{0.10, 0.05, -0.03};
  for (const auto &link : domain.links()) {
    const auto &row = ghost_plan.reconstruction(link.id);
    const double rho_wall = row.value(link.wall_intercept_m, density_view, 0U);
    const auto gradient =
        row.gradient(link.wall_intercept_m, density_view, 0U);
    const double reconstructed_normal =
        dot(gradient, link.solid_to_fluid_normal);
    const double expected_normal =
        dot(analytic_gradient, link.solid_to_fluid_normal);
    const double expected_wall = density(link.wall_intercept_m);
    HUNDUN_CHECK(rho_wall > 0.0 && std::isfinite(rho_wall));
    HUNDUN_CHECK(std::isfinite(reconstructed_normal));
    maximum_wall_value_error =
        std::max(maximum_wall_value_error,
                 std::abs(rho_wall - expected_wall));
    maximum_gradient_error =
        std::max(maximum_gradient_error,
                 std::abs(reconstructed_normal - expected_normal));
    maximum_zero_mutation_error =
        std::max(maximum_zero_mutation_error, std::abs(expected_normal));
    ++reconstructed_rows;
  }
  double reconstruction_values[3]{maximum_wall_value_error,
                                  maximum_gradient_error,
                                  maximum_zero_mutation_error};
  mpi.allreduce_fp64_in_place(reconstruction_values, 3U,
                              runtime::Fp64ReductionOperation::maximum);
  double row_count = static_cast<double>(reconstructed_rows);
  mpi.allreduce_fp64_in_place(&row_count, 1U,
                              runtime::Fp64ReductionOperation::sum);
  HUNDUN_CHECK(row_count > 0.0);
  HUNDUN_CHECK(reconstruction_values[0] <= 1.0e-12);
  HUNDUN_CHECK(reconstruction_values[1] <= 1.0e-10);
  HUNDUN_CHECK(reconstruction_values[2] > 1.0e-3);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  try {
#ifdef HUNDUN_STAGE3_IMMERSED_MATERIAL_TRANSACTION
    run_case(true);
#else
    run_case(false);
#endif
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 8;
  }
}
