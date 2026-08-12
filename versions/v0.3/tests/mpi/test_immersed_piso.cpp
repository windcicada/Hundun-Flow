// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/flow_material_density_piso_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/ib_surface_query.hpp"
#include "hundun/lin_bicgstab.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/flow_state_equality.hpp"

#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"
#include "src/fvm_immersed_boundary_authority_detail.hpp"
#include "tests/support/fvm_immersed_operator_test_access.hpp"
#include "src/ib_deterministic_qr_detail.hpp"
#include "src/ib_periodic_surface_window_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <variant>
#include <vector>

namespace {

using namespace hundun;

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("unsupported Task 9 rank count");
}

int shared_row_oracle_field_coordinate(int coordinate, int begin,
                                       int local_size, int global_size,
                                       int available_ghost_width) {
  int result = coordinate - begin;
  const int lower_image = result - global_size;
  const int upper_image = result + global_size;
  const auto inside = [&](int value) {
    return value >= -available_ghost_width &&
           value < local_size + available_ghost_width;
  };
  if (!inside(result) && inside(lower_image))
    result = lower_image;
  else if (!inside(result) && inside(upper_image))
    result = upper_image;
  if (!inside(result))
    throw runtime::Error("shared-row oracle donor exceeds field halo");
  return result;
}

int legacy_shared_row_oracle_field_coordinate(int coordinate, int begin,
                                              int end, int global_size) {
  if (coordinate < begin && begin == 0)
    coordinate += global_size;
  if (coordinate >= end && end == global_size)
    coordinate -= global_size;
  return coordinate - begin;
}

config::FlowCaseConfig periodic_case(int ranks) {
  config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type = config::SimulationType::variable_density_flow;
  config.density_model = config::DensityModel::constant;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = process_grid(ranks);
  config.mesh.cells = {12, 12, 12};
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.dynamic_viscosity_pa_s = 0.01;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t index = 0U; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type = config::BoundaryType::periodic;
  }
  return config;
}

config::FlowCaseConfig open_case(int ranks, double inlet_velocity_x = 0.0) {
  auto config = periodic_case(ranks);
  config.boundaries[0].type = config::BoundaryType::velocity_inlet;
  config.boundaries[0].velocity_m_per_s =
      runtime::Real3{inlet_velocity_x, 0.0, 0.0};
  config.boundaries[0].thermal_authority =
      config::InletThermalAuthority::enthalpy;
  config.boundaries[0].enthalpy_J_per_kg = 1.0;
  config.boundaries[0].scalar_values = std::vector<config::InletScalarValue>{};
  config.boundaries[1].type = config::BoundaryType::pressure_outlet;
  config.boundaries[1].pressure_perturbation_pa = 0.0;
  for (std::size_t patch = 2U; patch < config.boundaries.size(); ++patch)
    config.boundaries[patch].type = config::BoundaryType::symmetry;
  return config;
}

runtime::FieldDescriptor cell_field(const char *name,
                                    std::uint32_t components,
                                    int ghost_width) {
  return {name,
          "1",
          "stage3_task9",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost_width,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_field(const char *name,
                                    std::uint32_t components) {
  return {name,
          "1",
          "stage3_task9",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

std::vector<test::StlFixtureTriangle> internal_sphere() {
  return test::projected_octahedral_sphere({0.5, 0.5, 0.5}, 0.2, 3U);
}

std::vector<test::StlFixtureTriangle> internal_tetrahedron() {
  return test::projected_octahedral_sphere({0.48, 0.51, 0.49}, 0.22, 3U);
}

enum class FixtureKind : std::uint8_t {
  outside_periodic,
  outside_open,
  outside_tetrahedron
};

struct OutletFaceHistory final {
  double committed_velocity_x_m_per_s{};
  double history_velocity_x_m_per_s{};
};

struct WallCellHistory final {
  double committed_scale{};
  double history_scale{};
};

struct Fixture final {
  explicit Fixture(const runtime::MpiContext &mpi,
                   FixtureKind kind = FixtureKind::outside_periodic,
                   int ghost_width_reduction = 0,
                   double inlet_velocity_x = 0.0)
      : directory(kind == FixtureKind::outside_open ? "task9-piso-open"
                  : kind == FixtureKind::outside_tetrahedron
                      ? "task11-piso-shared-row"
                      : "task9-piso-periodic"),
        decomposition(runtime::StructuredDecomposition::create(
            mpi, {12, 12, 12},
            kind == FixtureKind::outside_open
                ? std::array<bool, 3>{false, false, false}
                : std::array<bool, 3>{true, true, true},
            runtime::DecompositionOptions{process_grid(mpi.size())})),
        topology(decomposition),
        geometry(topology,
                 mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0})),
        boundaries(boundary::BoundaryRegistry::create(
            kind == FixtureKind::outside_open
                ? open_case(mpi.size(), inlet_velocity_x)
                : periodic_case(mpi.size()),
            topology)) {
    std::string path_text;
    if (mpi.rank() == 0) {
      const auto root_path = directory.path() / "body.stl";
      auto triangles = kind == FixtureKind::outside_tetrahedron
                           ? internal_tetrahedron()
                           : internal_sphere();
      test::write_text(root_path, test::ascii_stl(triangles, "body"));
      path_text = root_path.string();
    }
    std::uint64_t path_size = path_text.size();
    HUNDUN_CHECK(MPI_Bcast(&path_size, 1, MPI_UINT64_T, 0, mpi.comm()) ==
                 MPI_SUCCESS);
    path_text.resize(static_cast<std::size_t>(path_size));
    HUNDUN_CHECK(MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()),
                           MPI_BYTE, 0, mpi.comm()) == MPI_SUCCESS);
    const std::filesystem::path path(path_text);
    surface.emplace(
        immersed::ImmersedSurface::load_collective(path, 1.0, mpi, 0));
    query.emplace(immersed::SurfaceQuery::create(*surface));
    domain.emplace(immersed::ImmersedDomain::create(
        *surface, *query, config::ImmersedFluidSide::outside, topology,
        geometry, boundaries, mpi));
    ghost_plan.emplace(immersed::GhostStencilPlan::create(
        *surface, *query, *domain, topology, geometry, decomposition, mpi));
    const auto required_ghost_width =
        static_cast<int>(ghost_plan->maximum_halo_reach());
    HUNDUN_CHECK(ghost_width_reduction >= 0);
    HUNDUN_CHECK(ghost_width_reduction <= required_ghost_width);
    const int ghost_width = required_ghost_width - ghost_width_reduction;

    fields.density = registry.declare_field(cell_field("rho", 1U, ghost_width));
    fields.velocity =
        registry.declare_field(cell_field("velocity", 3U, ghost_width));
    fields.mechanical_pressure =
        registry.declare_field(cell_field("pi", 1U, ghost_width));
    fields.face_velocity =
        registry.declare_field(face_field("face_velocity", 3U));
    fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
    registry.freeze();
  }

  flow::FlowState make_state(bool pressure_checkerboard = false,
                             bool corrupt_inactive_cell = false,
                             bool corrupt_inactive_face = false,
                             std::optional<OutletFaceHistory> outlet_history =
                                 std::nullopt,
                             std::optional<WallCellHistory> wall_history =
                                 std::nullopt) const {
    const auto local = decomposition.local_extent();
    auto state = flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0, flow::MomentumTimeOrder::backward_euler});
    flow::FlowLayerValues values;
    values.density.resize(topology.owned_cell_count(), 0.0);
    values.velocity.resize(topology.owned_cell_count() * 3U, 0.0);
    values.mechanical_pressure.resize(topology.owned_cell_count(), 0.0);
    values.face_velocity.resize(topology.local_face_count() * 3U, 0.0);
    values.face_mass_flux.resize(topology.local_face_count(), 0.0);
    for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
         ++cell)
      if (domain->region(cell) == immersed::CellRegion::fluid) {
        values.density[cell] = 1.0;
        if (pressure_checkerboard) {
          const auto index = topology.global_cell(cell);
          values.mechanical_pressure[cell] =
              ((index.x + index.y + index.z) % 2 == 0) ? 1.0 : -1.0;
        }
      }
    if (corrupt_inactive_cell) {
      bool changed = false;
      for (mesh::LocalCellId cell = 0U;
           cell < topology.owned_cell_count() && !changed; ++cell) {
        if (domain->region(cell) != immersed::CellRegion::fluid) {
          values.density[cell] = 2.0;
          changed = true;
        }
      }
      HUNDUN_CHECK(changed);
    }
    if (corrupt_inactive_face) {
      bool changed = false;
      for (mesh::LocalFaceId face = 0U;
           face < topology.local_face_count() && !changed; ++face) {
        const auto neighbour = topology.neighbour(face);
        const bool owner_active =
            domain->region(topology.owner(face)) == immersed::CellRegion::fluid;
        const bool neighbour_active =
            neighbour.has_value() &&
            domain->region(*neighbour) == immersed::CellRegion::fluid;
        if (!(owner_active && (!neighbour.has_value() || neighbour_active))) {
          values.face_mass_flux[face] = -0.0;
          changed = true;
        }
      }
      HUNDUN_CHECK(changed);
    }
    HUNDUN_CHECK(!outlet_history.has_value() || !wall_history.has_value());
    if (wall_history.has_value()) {
      auto history = values;
      for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
           ++cell) {
        if (domain->region(cell) != immersed::CellRegion::fluid)
          continue;
        const auto center = geometry.cell_center_m(cell);
        constexpr double pi = 3.141592653589793238462643383279502884;
        const std::array<double, 3> pattern{
            1.0 + 0.25 * center.x + 0.10 * center.y +
                0.4 * std::sin(2.0 * pi * center.x) *
                    std::sin(2.0 * pi * center.y) *
                    std::sin(2.0 * pi * center.z),
            -0.4 + 0.15 * center.y + 0.07 * center.z -
                0.3 * std::cos(2.0 * pi * center.x) *
                    std::sin(2.0 * pi * center.z),
            0.3 - 0.11 * center.x + 0.09 * center.z +
                0.25 * std::sin(2.0 * pi * center.x) *
                    std::cos(2.0 * pi * center.y)};
        for (std::size_t component = 0U; component < 3U; ++component) {
          values.velocity[static_cast<std::size_t>(cell) * 3U + component] =
              wall_history->committed_scale * pattern[component];
          history.velocity[static_cast<std::size_t>(cell) * 3U + component] =
              wall_history->history_scale * pattern[component];
        }
      }
      state.seed_accepted_layers(history, values);
      return state;
    }
    if (!outlet_history.has_value()) {
      state.seed_accepted_layers(values, values);
      return state;
    }
    auto history = values;
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      if (topology.neighbour(face).has_value())
        continue;
      const auto patch = topology.patch_id(face);
      if (!patch.has_value() ||
          boundaries.patch(*patch).pressure_rule() !=
              boundary::PressureRule::prescribed_value ||
          domain->region(topology.owner(face)) != immersed::CellRegion::fluid)
        continue;
      const auto area =
          geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
      values.face_velocity[static_cast<std::size_t>(face) * 3U] =
          outlet_history->committed_velocity_x_m_per_s;
      values.face_mass_flux[face] =
          outlet_history->committed_velocity_x_m_per_s * area.x;
      history.face_velocity[static_cast<std::size_t>(face) * 3U] =
          outlet_history->history_velocity_x_m_per_s;
      history.face_mass_flux[face] =
          outlet_history->history_velocity_x_m_per_s * area.x;
    }
    state.seed_accepted_layers(history, values);
    return state;
  }

  flow::FlowState make_body_fitted_state() const {
    const auto local = decomposition.local_extent();
    auto state = flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0, flow::MomentumTimeOrder::backward_euler});
    flow::FlowLayerValues values;
    values.density.resize(topology.owned_cell_count(), 1.0);
    values.velocity.resize(topology.owned_cell_count() * 3U, 0.0);
    values.mechanical_pressure.resize(topology.owned_cell_count(), 0.0);
    values.face_velocity.resize(topology.local_face_count() * 3U, 0.0);
    values.face_mass_flux.resize(topology.local_face_count(), 0.0);
    state.seed_accepted_layers(values, values);
    return state;
  }

  flow::FlowState make_energy_probe_state(double phase) const {
    const auto local = decomposition.local_extent();
    auto state = flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0, flow::MomentumTimeOrder::backward_euler});
    flow::FlowLayerValues values;
    values.density.resize(topology.owned_cell_count(), 0.0);
    values.velocity.resize(topology.owned_cell_count() * 3U, 0.0);
    values.mechanical_pressure.resize(topology.owned_cell_count(), 0.0);
    values.face_velocity.resize(topology.local_face_count() * 3U, 0.0);
    values.face_mass_flux.resize(topology.local_face_count(), 0.0);
    constexpr double two_pi =
        6.2831853071795864769252867665590057683943387987502;
    const auto probe_velocity = [&](mesh::LocalCellId cell) {
      std::array<double, 3> result{};
      if (domain->region(cell) != immersed::CellRegion::fluid)
        return result;
      const auto center = geometry.cell_center_m(cell);
      const bool in_wall_support_shell =
          center.x >= 0.20 && center.x <= 0.80 && center.y >= 0.20 &&
          center.y <= 0.80 && center.z >= 0.20 && center.z <= 0.80;
      if (in_wall_support_shell)
        return result;
      const auto gid = topology.global_cell_id(cell);
      const double id_phase =
          0.00037 * static_cast<double>((gid % 104729U) + 1U);
      result[0] = 0.31 * std::sin(two_pi * (center.y + phase)) +
                  0.07 * std::cos(two_pi * center.z + id_phase);
      result[1] = -0.23 * std::cos(two_pi * (center.x - 0.5 * phase)) +
                  0.05 * std::sin(two_pi * center.z - id_phase);
      result[2] =
          0.19 * std::sin(two_pi * (center.x + center.y + phase)) -
          0.03 * std::cos(two_pi * center.z + 0.5 * id_phase);
      return result;
    };
    for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
         ++cell) {
      if (domain->region(cell) != immersed::CellRegion::fluid)
        continue;
      values.density[cell] = 1.0;
      const auto center = geometry.cell_center_m(cell);
      const auto gid = topology.global_cell_id(cell);
      const double id_phase =
          0.00037 * static_cast<double>((gid % 104729U) + 1U);
      const auto velocity = probe_velocity(cell);
      for (std::size_t component = 0U; component < velocity.size(); ++component)
        values.velocity[cell * 3U + component] = velocity[component];
      values.mechanical_pressure[cell] =
          0.83 * std::sin(two_pi * (center.x + 0.25 * phase)) +
          0.29 * std::cos(two_pi * (center.y - center.z)) + id_phase;
    }
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      const auto neighbour = topology.neighbour(face);
      const bool owner_active =
          domain->region(topology.owner(face)) == immersed::CellRegion::fluid;
      const bool neighbour_active =
          neighbour.has_value() &&
          domain->region(*neighbour) == immersed::CellRegion::fluid;
      if (!owner_active || !neighbour_active)
        continue;
      const auto owner_velocity = probe_velocity(topology.owner(face));
      const auto neighbour_velocity = probe_velocity(*neighbour);
      std::array<double, 3> centered_velocity{};
      for (std::size_t component = 0U; component < 3U; ++component) {
        centered_velocity[component] =
            0.5 * (owner_velocity[component] + neighbour_velocity[component]);
        values.face_velocity[static_cast<std::size_t>(face) * 3U + component] =
            centered_velocity[component];
      }
      const auto area =
          geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
      values.face_mass_flux[face] = centered_velocity[0] * area.x +
                                    centered_velocity[1] * area.y +
                                    centered_velocity[2] * area.z;
    }
    state.seed_accepted_layers(values, values);
    return state;
  }

  test::Stage3TemporaryDirectory directory;
  runtime::StructuredDecomposition decomposition;
  mesh::MeshTopology topology;
  mesh::MeshGeometry geometry;
  boundary::BoundaryRegistry boundaries;
  std::optional<immersed::ImmersedSurface> surface;
  std::optional<immersed::SurfaceQuery> query;
  std::optional<immersed::ImmersedDomain> domain;
  std::optional<immersed::GhostStencilPlan> ghost_plan;
  immersed::LocalFlowPatternTransform transform;
  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
};

std::size_t neighbour_occurrence(runtime::Int3 fluid, runtime::Int3 solid) {
  const runtime::Int3 offset{solid.x - fluid.x, solid.y - fluid.y,
                             solid.z - fluid.z};
  if (offset.x == -1 && offset.y == 0 && offset.z == 0)
    return 0U;
  if (offset.x == 1 && offset.y == 0 && offset.z == 0)
    return 1U;
  if (offset.x == 0 && offset.y == -1 && offset.z == 0)
    return 2U;
  if (offset.x == 0 && offset.y == 1 && offset.z == 0)
    return 3U;
  if (offset.x == 0 && offset.y == 0 && offset.z == -1)
    return 4U;
  if (offset.x == 0 && offset.y == 0 && offset.z == 1)
    return 5U;
  throw runtime::Error("Task 9 wall link is not a logical neighbour");
}

double distance(runtime::Real3 lhs, runtime::Real3 rhs) {
  const double dx = lhs.x - rhs.x;
  const double dy = lhs.y - rhs.y;
  const double dz = lhs.z - rhs.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void check_wall_effective_measures(const flow::FixedStepImmersedFlow &immersed_flow,
                                   const Fixture &fixture) {
  const auto snapshots =
      flow::test::ImmersedFlowTestAccess::wall_effective_measures(immersed_flow);
  std::vector<const immersed::ImmersedLink *> links;
  links.reserve(fixture.domain->links().size());
  for (const auto &link : fixture.domain->links())
    links.push_back(&link);
  std::sort(links.begin(), links.end(),
            [](const auto *lhs, const auto *rhs) { return lhs->id < rhs->id; });
  std::size_t expected_index = 0U;
  double legacy_l1_measure_difference = 0.0;
  for (const auto *link_pointer : links) {
    const auto &link = *link_pointer;
    const auto fluid = fixture.topology.find_local_cell(link.fluid_cell);
    if (!fluid.has_value() ||
        fixture.topology.cell_ownership(*fluid) != mesh::EntityOwnership::owned)
      continue;
    const auto solid = fixture.topology.find_local_cell(link.solid_cell);
    HUNDUN_CHECK(solid.has_value());
    mesh::LocalFaceId link_face = std::numeric_limits<mesh::LocalFaceId>::max();
    for (mesh::LocalFaceId face = 0U;
         face < fixture.topology.local_face_count(); ++face) {
      const auto neighbour = fixture.topology.neighbour(face);
      if (!neighbour.has_value())
        continue;
      const auto owner_id =
          fixture.topology.global_cell_id(fixture.topology.owner(face));
      const auto neighbour_id = fixture.topology.global_cell_id(*neighbour);
      if ((owner_id == link.fluid_cell && neighbour_id == link.solid_cell) ||
          (owner_id == link.solid_cell && neighbour_id == link.fluid_cell)) {
        link_face = face;
        break;
      }
    }
    HUNDUN_CHECK(link_face != std::numeric_limits<mesh::LocalFaceId>::max());
    const auto occurrence =
        neighbour_occurrence(fixture.topology.global_cell(*fluid),
                             fixture.topology.global_cell(*solid));
    const auto fluid_center = fixture.geometry.cell_center_m(*fluid);
    const auto solid_center = fixture.geometry.cell_center_m(*solid);
    const auto face_center = fixture.geometry.face_center_m(link_face);
    const double center_distance = distance(fluid_center, solid_center);
    const double owner_to_face = distance(fluid_center, face_center);
    HUNDUN_CHECK(center_distance > 0.0);
    HUNDUN_CHECK(owner_to_face > 0.0 && owner_to_face < center_distance);
    auto area = fixture.geometry.face_area_vector_m2(
        link_face, mesh::FaceSide::owner);
    if (fixture.topology.owner(link_face) != *fluid)
      area = {-area.x, -area.y, -area.z};
    const double expected =
        -(area.x * link.solid_to_fluid_normal.x +
          area.y * link.solid_to_fluid_normal.y +
          area.z * link.solid_to_fluid_normal.z);
    HUNDUN_CHECK(expected > 0.0 && std::isfinite(expected));
    immersed::LocalCoefficientRow background{};
    background.neighbour[occurrence] = owner_to_face / center_distance;
    background.diagonal = 1.0 - background.neighbour[occurrence];
    const auto transformed = fixture.transform.transform_full(
        background, link.fluid_to_wall_fraction, link.solid_to_fluid_normal);
    double legacy_l1_factor =
        std::abs(transformed.diagonal) + std::abs(transformed.source);
    for (const double coefficient : transformed.neighbour)
      legacy_l1_factor += std::abs(coefficient);
    const double legacy_l1_measure =
        legacy_l1_factor * fixture.geometry.face_area_m2(link_face);
    legacy_l1_measure_difference = std::max(
        legacy_l1_measure_difference,
        std::abs(expected - legacy_l1_measure));

    HUNDUN_CHECK(expected_index < snapshots.size());
    HUNDUN_CHECK(snapshots[expected_index].link == link.id);
    HUNDUN_CHECK_NEAR(snapshots[expected_index].effective_measure_m2, expected,
                      64.0 * std::numeric_limits<double>::epsilon() *
                          std::max(1.0, std::abs(expected)));
    ++expected_index;
  }
  HUNDUN_CHECK(expected_index == snapshots.size());
  HUNDUN_CHECK(legacy_l1_measure_difference > 1.0e-10);
}

void check_state_equal(const flow::FlowState &state,
                       const flow::FlowLayerValues &history,
                       const flow::FlowLayerValues &committed,
                       const flow::FlowLayerValues &trial,
                       const flow::AcceptedStepMetadata &metadata) {
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      state.snapshot(flow::FlowLayer::history), history));
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      state.snapshot(flow::FlowLayer::committed), committed));
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      state.snapshot(flow::FlowLayer::trial), trial));
  HUNDUN_CHECK(
      test::accepted_step_metadata_bitwise_equal(state.metadata(), metadata));
}

bool final_residual_reports_bitwise_equal(
    const flow::test::ImmersedFlowFinalMomentumPressureResidualReport &left,
    const flow::test::ImmersedFlowFinalMomentumPressureResidualReport &right) {
  return left.active_global_cell_ids == right.active_global_cell_ids &&
         test::fp64_vector_bitwise_equal(left.used_residual,
                                         right.used_residual) &&
         test::fp64_vector_bitwise_equal(left.direct_residual,
                                         right.direct_residual) &&
         test::fp64_bits(left.used_l2) == test::fp64_bits(right.used_l2) &&
         test::fp64_bits(left.direct_l2) == test::fp64_bits(right.direct_l2) &&
         test::fp64_bits(left.difference_l2) ==
             test::fp64_bits(right.difference_l2) &&
         test::fp64_bits(left.used_linf) == test::fp64_bits(right.used_linf) &&
         test::fp64_bits(left.direct_linf) ==
             test::fp64_bits(right.direct_linf) &&
         test::fp64_bits(left.difference_linf) ==
             test::fp64_bits(right.difference_linf) &&
         left.maximum_difference_global_cell_id ==
             right.maximum_difference_global_cell_id &&
         left.maximum_difference_component ==
             right.maximum_difference_component &&
         test::fp64_bits(left.maximum_difference_used_value) ==
             test::fp64_bits(right.maximum_difference_used_value) &&
         test::fp64_bits(left.maximum_difference_direct_value) ==
             test::fp64_bits(right.maximum_difference_direct_value);
}

bool accepted_inner_work_equal(
    const flow::test::ImmersedFlowExactPredictorWorkspaceSnapshot &left,
    const flow::test::ImmersedFlowExactPredictorWorkspaceSnapshot &right) noexcept {
  return left.accepted_inner_response_count ==
             right.accepted_inner_response_count &&
         left.accepted_inner_solve_count == right.accepted_inner_solve_count &&
         left.accepted_inner_iteration_count ==
             right.accepted_inner_iteration_count &&
         left.accepted_inner_matvec_count ==
             right.accepted_inner_matvec_count &&
         left.accepted_inner_preconditioner_apply_count ==
             right.accepted_inner_preconditioner_apply_count &&
         left.accepted_inner_global_reduction_count ==
             right.accepted_inner_global_reduction_count;
}

bool final_residual_vectors_agree(const std::vector<double> &used,
                                  const std::vector<double> &direct) noexcept {
  if (used.empty() || used.size() != direct.size() || used.size() % 3U != 0U)
    return false;
  double used_square = 0.0;
  double direct_square = 0.0;
  double difference_square = 0.0;
  double used_linf = 0.0;
  double direct_linf = 0.0;
  double difference_linf = 0.0;
  for (std::size_t index = 0U; index < used.size(); ++index) {
    if (!std::isfinite(used[index]) || !std::isfinite(direct[index]))
      return false;
    const double difference = used[index] - direct[index];
    used_square += used[index] * used[index];
    direct_square += direct[index] * direct[index];
    difference_square += difference * difference;
    used_linf = std::max(used_linf, std::abs(used[index]));
    direct_linf = std::max(direct_linf, std::abs(direct[index]));
    difference_linf = std::max(difference_linf, std::abs(difference));
  }
  const double used_l2 = std::sqrt(used_square);
  const double direct_l2 = std::sqrt(direct_square);
  const double difference_l2 = std::sqrt(difference_square);
  if (!std::isfinite(used_l2) || !std::isfinite(direct_l2) ||
      !std::isfinite(difference_l2))
    return false;
  const double factor =
      512.0 * std::numeric_limits<double>::epsilon();
  return difference_l2 <= factor * std::max({1.0, used_l2, direct_l2}) &&
         difference_linf <=
             factor * std::max({1.0, used_linf, direct_linf});
}

bool routes_agree(
    const flow::test::ImmersedFlowFinalMomentumPressureResidualReport &report) {
  if (report.active_global_cell_ids.empty() ||
      report.used_residual.size() != report.active_global_cell_ids.size() * 3U ||
      report.direct_residual.size() != report.used_residual.size() ||
      !std::isfinite(report.used_l2) || !std::isfinite(report.direct_l2) ||
      !std::isfinite(report.difference_l2) ||
      !std::isfinite(report.used_linf) || !std::isfinite(report.direct_linf) ||
      !std::isfinite(report.difference_linf))
    return false;
  const double factor =
      512.0 * std::numeric_limits<double>::epsilon();
  return report.difference_l2 <=
             factor * std::max({1.0, report.used_l2, report.direct_l2}) &&
         report.difference_linf <=
             factor * std::max({1.0, report.used_linf, report.direct_linf});
}

bool pressure_flux_reports_bitwise_equal(
    const flow::test::ImmersedFlowPressureFluxIdentityReport &left,
    const flow::test::ImmersedFlowPressureFluxIdentityReport &right) {
  return left.active_global_cell_ids == right.active_global_cell_ids &&
         test::fp64_vector_bitwise_equal(
             left.operator_residual_per_volume,
             right.operator_residual_per_volume) &&
         test::fp64_vector_bitwise_equal(
             left.routed_flux_divergence_per_volume,
             right.routed_flux_divergence_per_volume) &&
         test::fp64_bits(left.operator_l2) ==
             test::fp64_bits(right.operator_l2) &&
         test::fp64_bits(left.routed_l2) ==
             test::fp64_bits(right.routed_l2) &&
         test::fp64_bits(left.difference_l2) ==
             test::fp64_bits(right.difference_l2) &&
         test::fp64_bits(left.operator_linf) ==
             test::fp64_bits(right.operator_linf) &&
         test::fp64_bits(left.routed_linf) ==
             test::fp64_bits(right.routed_linf) &&
         test::fp64_bits(left.difference_linf) ==
             test::fp64_bits(right.difference_linf) &&
         left.maximum_difference_global_cell_id ==
             right.maximum_difference_global_cell_id &&
         test::fp64_bits(left.maximum_difference_operator_value) ==
             test::fp64_bits(right.maximum_difference_operator_value) &&
         test::fp64_bits(left.maximum_difference_routed_value) ==
             test::fp64_bits(right.maximum_difference_routed_value);
}

bool pressure_flux_routes_agree(
    const flow::test::ImmersedFlowPressureFluxIdentityReport &report) noexcept {
  if (report.active_global_cell_ids.empty() ||
      report.operator_residual_per_volume.size() !=
          report.active_global_cell_ids.size() ||
      report.routed_flux_divergence_per_volume.size() !=
          report.operator_residual_per_volume.size())
    return false;
  double operator_square = 0.0;
  double routed_square = 0.0;
  double difference_square = 0.0;
  double operator_linf = 0.0;
  double routed_linf = 0.0;
  double difference_linf = 0.0;
  double maximum_difference = -1.0;
  std::uint64_t maximum_gid = 0U;
  double maximum_operator_value = 0.0;
  double maximum_routed_value = 0.0;
  for (std::size_t index = 0U;
       index < report.operator_residual_per_volume.size(); ++index) {
    const double operator_value = report.operator_residual_per_volume[index];
    const double routed_value =
        report.routed_flux_divergence_per_volume[index];
    if (!std::isfinite(operator_value) || !std::isfinite(routed_value))
      return false;
    const double difference = operator_value - routed_value;
    operator_square += operator_value * operator_value;
    routed_square += routed_value * routed_value;
    difference_square += difference * difference;
    operator_linf = std::max(operator_linf, std::abs(operator_value));
    routed_linf = std::max(routed_linf, std::abs(routed_value));
    difference_linf = std::max(difference_linf, std::abs(difference));
    if (std::abs(difference) > maximum_difference) {
      maximum_difference = std::abs(difference);
      maximum_gid = report.active_global_cell_ids[index];
      maximum_operator_value = operator_value;
      maximum_routed_value = routed_value;
    }
  }
  const double operator_l2 = std::sqrt(operator_square);
  const double routed_l2 = std::sqrt(routed_square);
  const double difference_l2 = std::sqrt(difference_square);
  const double factor = 64.0 * std::numeric_limits<double>::epsilon();
  return std::isfinite(operator_l2) && std::isfinite(routed_l2) &&
         std::isfinite(difference_l2) &&
         test::fp64_bits(operator_l2) ==
             test::fp64_bits(report.operator_l2) &&
         test::fp64_bits(routed_l2) == test::fp64_bits(report.routed_l2) &&
         test::fp64_bits(difference_l2) ==
             test::fp64_bits(report.difference_l2) &&
         test::fp64_bits(operator_linf) ==
             test::fp64_bits(report.operator_linf) &&
         test::fp64_bits(routed_linf) ==
             test::fp64_bits(report.routed_linf) &&
         test::fp64_bits(difference_linf) ==
             test::fp64_bits(report.difference_linf) &&
         maximum_gid == report.maximum_difference_global_cell_id &&
         test::fp64_bits(maximum_operator_value) ==
             test::fp64_bits(report.maximum_difference_operator_value) &&
         test::fp64_bits(maximum_routed_value) ==
             test::fp64_bits(report.maximum_difference_routed_value) &&
         difference_l2 <=
             factor * std::max({1.0, operator_l2, routed_l2}) &&
         difference_linf <=
             factor * std::max({1.0, operator_linf, routed_linf});
}

bool cell_pressure_correction_records_bitwise_equal(
    const flow::test::ImmersedFlowCellPressureCorrectionRecord &left,
    const flow::test::ImmersedFlowCellPressureCorrectionRecord &right) {
  return left.active_global_cell_ids == right.active_global_cell_ids &&
         test::fp64_vector_bitwise_equal(left.exact_velocity_change,
                                         right.exact_velocity_change) &&
         test::fp64_vector_bitwise_equal(
             left.momentum_operator_velocity_change,
             right.momentum_operator_velocity_change) &&
         test::fp64_vector_bitwise_equal(left.lfp_pressure_residual_change,
                                         right.lfp_pressure_residual_change) &&
         test::fp64_vector_bitwise_equal(left.pressure_before_pa,
                                         right.pressure_before_pa) &&
         test::fp64_vector_bitwise_equal(left.pressure_correction_pa,
                                         right.pressure_correction_pa) &&
         test::fp64_vector_bitwise_equal(left.pressure_after_pa,
                                         right.pressure_after_pa) &&
         test::fp64_vector_bitwise_equal(left.pressure_rhs_per_volume,
                                         right.pressure_rhs_per_volume) &&
         test::fp64_vector_bitwise_equal(
             left.compact_pressure_action_per_volume,
             right.compact_pressure_action_per_volume) &&
         test::fp64_vector_bitwise_equal(
             left.interface_pressure_action_per_volume,
             right.interface_pressure_action_per_volume) &&
         test::fp64_vector_bitwise_equal(
             left.hybrid_pressure_action_per_volume,
             right.hybrid_pressure_action_per_volume) &&
         test::fp64_vector_bitwise_equal(left.affine_wall_source_per_volume,
                                         right.affine_wall_source_per_volume) &&
         test::fp64_vector_bitwise_equal(
             left.compact_predictor_defect_per_volume,
             right.compact_predictor_defect_per_volume) &&
         test::fp64_vector_bitwise_equal(
             left.hybrid_predictor_defect_per_volume,
             right.hybrid_predictor_defect_per_volume) &&
         test::fp64_bits(left.momentum_l2) ==
             test::fp64_bits(right.momentum_l2) &&
         test::fp64_bits(left.pressure_l2) ==
             test::fp64_bits(right.pressure_l2) &&
         test::fp64_bits(left.closure_l2) ==
             test::fp64_bits(right.closure_l2) &&
         test::fp64_bits(left.momentum_linf) ==
             test::fp64_bits(right.momentum_linf) &&
         test::fp64_bits(left.pressure_linf) ==
             test::fp64_bits(right.pressure_linf) &&
         test::fp64_bits(left.closure_linf) ==
             test::fp64_bits(right.closure_linf) &&
         left.maximum_closure_global_cell_id ==
             right.maximum_closure_global_cell_id &&
         left.maximum_closure_component == right.maximum_closure_component &&
         test::fp64_bits(left.maximum_closure_momentum_value) ==
             test::fp64_bits(right.maximum_closure_momentum_value) &&
         test::fp64_bits(left.maximum_closure_pressure_value) ==
             test::fp64_bits(right.maximum_closure_pressure_value) &&
         test::fp64_bits(left.wall_predictor_mass_flux_l2_kg_per_s) ==
             test::fp64_bits(right.wall_predictor_mass_flux_l2_kg_per_s) &&
         test::fp64_bits(left.wall_predictor_mass_flux_linf_kg_per_s) ==
             test::fp64_bits(right.wall_predictor_mass_flux_linf_kg_per_s) &&
         test::fp64_bits(left.wall_correction_gradient_l2_pa_per_m) ==
             test::fp64_bits(right.wall_correction_gradient_l2_pa_per_m) &&
         test::fp64_bits(left.wall_correction_gradient_linf_pa_per_m) ==
             test::fp64_bits(right.wall_correction_gradient_linf_pa_per_m);
}

bool cell_pressure_correction_reports_bitwise_equal(
    const flow::test::ImmersedFlowCellPressureCorrectionReport &left,
    const flow::test::ImmersedFlowCellPressureCorrectionReport &right) {
  if (left.correctors.size() != right.correctors.size())
    return false;
  for (std::size_t index = 0U; index < left.correctors.size(); ++index)
    if (!cell_pressure_correction_records_bitwise_equal(
            left.correctors[index], right.correctors[index]))
      return false;
  return test::fp64_bits(left.inter_corrector_authority_difference_l2) ==
         test::fp64_bits(right.inter_corrector_authority_difference_l2);
}

bool cell_pressure_correction_routes_agree(
    const flow::test::ImmersedFlowCellPressureCorrectionRecord &record) noexcept {
  if (record.active_global_cell_ids.empty() ||
      record.active_global_cell_ids.size() >
          std::numeric_limits<std::size_t>::max() / 3U ||
      record.exact_velocity_change.size() !=
          record.active_global_cell_ids.size() * 3U ||
      record.momentum_operator_velocity_change.size() !=
          record.active_global_cell_ids.size() * 3U ||
      record.lfp_pressure_residual_change.size() !=
          record.momentum_operator_velocity_change.size() ||
      record.pressure_before_pa.size() !=
          record.active_global_cell_ids.size() ||
      record.pressure_correction_pa.size() !=
          record.active_global_cell_ids.size() ||
      record.pressure_after_pa.size() != record.active_global_cell_ids.size() ||
      record.pressure_rhs_per_volume.size() !=
          record.active_global_cell_ids.size() ||
      record.compact_pressure_action_per_volume.size() !=
          record.active_global_cell_ids.size() ||
      record.interface_pressure_action_per_volume.size() !=
          record.active_global_cell_ids.size() ||
      record.hybrid_pressure_action_per_volume.size() !=
          record.active_global_cell_ids.size() ||
      record.affine_wall_source_per_volume.size() !=
          record.active_global_cell_ids.size() ||
      record.compact_predictor_defect_per_volume.size() !=
          record.active_global_cell_ids.size() ||
      record.hybrid_predictor_defect_per_volume.size() !=
          record.active_global_cell_ids.size())
    return false;
  double momentum_square = 0.0;
  double pressure_square = 0.0;
  double closure_square = 0.0;
  double momentum_linf = 0.0;
  double pressure_linf = 0.0;
  double closure_linf = 0.0;
  double maximum_closure = -1.0;
  std::uint64_t maximum_gid = 0U;
  std::uint32_t maximum_component = 0U;
  double maximum_momentum_value = 0.0;
  double maximum_pressure_value = 0.0;
  double predictor_identity_scale = 1.0;
  double predictor_identity_linf = 0.0;
  for (std::size_t row = 0U; row < record.active_global_cell_ids.size();
       ++row) {
    const double rhs = record.pressure_rhs_per_volume[row];
    const double compact = record.compact_pressure_action_per_volume[row];
    const double interface = record.interface_pressure_action_per_volume[row];
    const double hybrid = record.hybrid_pressure_action_per_volume[row];
    const double compact_defect =
        record.compact_predictor_defect_per_volume[row];
    const double hybrid_defect =
        record.hybrid_predictor_defect_per_volume[row];
    const double hybrid_identity = hybrid - compact - interface;
    const double compact_defect_identity = compact_defect - rhs - compact;
    const double hybrid_defect_identity =
        hybrid_defect - compact_defect - interface;
    if (!std::isfinite(rhs) || !std::isfinite(compact) ||
        !std::isfinite(interface) || !std::isfinite(hybrid) ||
        !std::isfinite(record.affine_wall_source_per_volume[row]) ||
        !std::isfinite(compact_defect) || !std::isfinite(hybrid_defect) ||
        !std::isfinite(hybrid_identity) ||
        !std::isfinite(compact_defect_identity) ||
        !std::isfinite(hybrid_defect_identity))
      return false;
    predictor_identity_scale =
        std::max({predictor_identity_scale, std::abs(rhs), std::abs(compact),
                  std::abs(interface), std::abs(hybrid),
                  std::abs(compact_defect), std::abs(hybrid_defect)});
    predictor_identity_linf =
        std::max({predictor_identity_linf, std::abs(hybrid_identity),
                  std::abs(compact_defect_identity),
                  std::abs(hybrid_defect_identity)});
  }
  for (std::size_t index = 0U;
       index < record.momentum_operator_velocity_change.size(); ++index) {
    if (!std::isfinite(record.exact_velocity_change[index]))
      return false;
    const double momentum = record.momentum_operator_velocity_change[index];
    const double pressure = record.lfp_pressure_residual_change[index];
    const double closure = momentum + pressure;
    if (!std::isfinite(momentum) || !std::isfinite(pressure) ||
        !std::isfinite(closure))
      return false;
    momentum_square += momentum * momentum;
    pressure_square += pressure * pressure;
    closure_square += closure * closure;
    momentum_linf = std::max(momentum_linf, std::abs(momentum));
    pressure_linf = std::max(pressure_linf, std::abs(pressure));
    closure_linf = std::max(closure_linf, std::abs(closure));
    if (std::abs(closure) > maximum_closure) {
      maximum_closure = std::abs(closure);
      maximum_gid = record.active_global_cell_ids[index / 3U];
      maximum_component = static_cast<std::uint32_t>(index % 3U);
      maximum_momentum_value = momentum;
      maximum_pressure_value = pressure;
    }
  }
  const double momentum_l2 = std::sqrt(momentum_square);
  const double pressure_l2 = std::sqrt(pressure_square);
  const double closure_l2 = std::sqrt(closure_square);
  const double factor = 512.0 * std::numeric_limits<double>::epsilon();
  return test::fp64_bits(momentum_l2) == test::fp64_bits(record.momentum_l2) &&
         test::fp64_bits(pressure_l2) == test::fp64_bits(record.pressure_l2) &&
         test::fp64_bits(closure_l2) == test::fp64_bits(record.closure_l2) &&
         test::fp64_bits(momentum_linf) ==
             test::fp64_bits(record.momentum_linf) &&
         test::fp64_bits(pressure_linf) ==
             test::fp64_bits(record.pressure_linf) &&
         test::fp64_bits(closure_linf) ==
             test::fp64_bits(record.closure_linf) &&
         maximum_gid == record.maximum_closure_global_cell_id &&
         maximum_component == record.maximum_closure_component &&
         test::fp64_bits(maximum_momentum_value) ==
             test::fp64_bits(record.maximum_closure_momentum_value) &&
         test::fp64_bits(maximum_pressure_value) ==
             test::fp64_bits(record.maximum_closure_pressure_value) &&
         predictor_identity_linf <= factor * predictor_identity_scale &&
         closure_l2 <=
             factor * std::max({1.0, momentum_l2, pressure_l2}) &&
         closure_linf <=
             factor * std::max({1.0, momentum_linf, pressure_linf});
}

std::array<double, 18> spatial_energy_scalars(
    const flow::test::ImmersedFlowSpatialEnergyReport &report) {
  return {report.kinetic_energy_J,
          report.total_power_W,
          report.convective_power_W,
          report.centered_convective_power_W,
          report.reconstruction_power_W,
          report.pressure_power_W,
          report.pressure_continuity_power_W,
          report.pressure_adjoint_defect_W,
          report.viscous_power_W,
          report.implicit_viscous_reference_power_W,
          report.residual_closure_l2_N,
          report.residual_closure_linf_N,
          report.mass_divergence_l2_kg_per_s,
          report.mass_divergence_linf_kg_per_s,
          report.stationary_wall_flux_linf_kg_per_s,
          report.physical_boundary_flux_linf_kg_per_s,
          report.maximum_closure_total_value_N,
          report.maximum_closure_parts_value_N};
}

bool fp64_arrays_bitwise_equal(const std::array<double, 18> &left,
                               const std::array<double, 18> &right) noexcept {
  for (std::size_t index = 0U; index < left.size(); ++index)
    if (test::fp64_bits(left[index]) != test::fp64_bits(right[index]))
      return false;
  return true;
}

bool spatial_energy_reports_bitwise_equal(
    const flow::test::ImmersedFlowSpatialEnergyReport &left,
    const flow::test::ImmersedFlowSpatialEnergyReport &right) noexcept {
  return left.active_global_cell_ids == right.active_global_cell_ids &&
         test::fp64_vector_bitwise_equal(left.cell_volume_m3,
                                         right.cell_volume_m3) &&
         test::fp64_vector_bitwise_equal(left.velocity_m_per_s,
                                         right.velocity_m_per_s) &&
         test::fp64_vector_bitwise_equal(left.pressure_pa,
                                         right.pressure_pa) &&
         test::fp64_vector_bitwise_equal(left.mass_divergence_kg_per_s,
                                         right.mass_divergence_kg_per_s) &&
         test::fp64_vector_bitwise_equal(left.total_residual_N,
                                         right.total_residual_N) &&
         test::fp64_vector_bitwise_equal(left.convective_residual_N,
                                         right.convective_residual_N) &&
         test::fp64_vector_bitwise_equal(left.pressure_residual_N,
                                         right.pressure_residual_N) &&
         test::fp64_vector_bitwise_equal(left.viscous_residual_N,
                                         right.viscous_residual_N) &&
         test::fp64_vector_bitwise_equal(
             left.implicit_viscous_reference_residual_N,
             right.implicit_viscous_reference_residual_N) &&
         fp64_arrays_bitwise_equal(spatial_energy_scalars(left),
                                   spatial_energy_scalars(right)) &&
         left.maximum_closure_global_cell_id ==
             right.maximum_closure_global_cell_id &&
         left.maximum_closure_component == right.maximum_closure_component;
}

bool spatial_energy_report_is_self_consistent(
    const flow::test::ImmersedFlowSpatialEnergyReport &report,
    double rho_ref_kg_per_m3) noexcept {
  const std::size_t count = report.active_global_cell_ids.size();
  if (count == 0U || !(rho_ref_kg_per_m3 > 0.0) ||
      report.cell_volume_m3.size() != count ||
      report.velocity_m_per_s.size() != count * 3U ||
      report.pressure_pa.size() != count ||
      report.mass_divergence_kg_per_s.size() != count ||
      report.total_residual_N.size() != count * 3U ||
      report.convective_residual_N.size() != count * 3U ||
      report.pressure_residual_N.size() != count * 3U ||
      report.viscous_residual_N.size() != count * 3U ||
      report.implicit_viscous_reference_residual_N.size() != count * 3U ||
      !std::is_sorted(report.active_global_cell_ids.begin(),
                      report.active_global_cell_ids.end()) ||
      std::adjacent_find(report.active_global_cell_ids.begin(),
                         report.active_global_cell_ids.end()) !=
          report.active_global_cell_ids.end())
    return false;

  flow::test::ImmersedFlowSpatialEnergyReport recomputed;
  double closure_square = 0.0;
  double divergence_square = 0.0;
  double closure_scale = 0.0;
  double maximum_closure = -1.0;
  for (std::size_t row = 0U; row < count; ++row) {
    const double volume = report.cell_volume_m3[row];
    const double pressure = report.pressure_pa[row];
    const double divergence = report.mass_divergence_kg_per_s[row];
    if (!(volume > 0.0) || !std::isfinite(volume) ||
        !std::isfinite(pressure) || !std::isfinite(divergence))
      return false;
    divergence_square += divergence * divergence;
    recomputed.mass_divergence_linf_kg_per_s =
        std::max(recomputed.mass_divergence_linf_kg_per_s,
                 std::abs(divergence));
    double speed_squared = 0.0;
    for (std::size_t component = 0U; component < 3U; ++component) {
      const auto offset = row * 3U + component;
      const double velocity = report.velocity_m_per_s[offset];
      const double total = report.total_residual_N[offset];
      const double convective = report.convective_residual_N[offset];
      const double pressure_residual = report.pressure_residual_N[offset];
      const double viscous = report.viscous_residual_N[offset];
      const double implicit_viscous_reference =
          report.implicit_viscous_reference_residual_N[offset];
      const double parts = convective + pressure_residual + viscous;
      const double closure = total - parts;
      if (!std::isfinite(velocity) || !std::isfinite(total) ||
          !std::isfinite(convective) || !std::isfinite(pressure_residual) ||
          !std::isfinite(viscous) ||
          !std::isfinite(implicit_viscous_reference) ||
          !std::isfinite(closure))
        return false;
      speed_squared += velocity * velocity;
      recomputed.total_power_W += velocity * total;
      recomputed.convective_power_W += velocity * convective;
      recomputed.pressure_power_W += velocity * pressure_residual;
      recomputed.viscous_power_W += velocity * viscous;
      recomputed.implicit_viscous_reference_power_W +=
          velocity * implicit_viscous_reference;
      closure_square += closure * closure;
      closure_scale += std::abs(total) + std::abs(parts);
      recomputed.residual_closure_linf_N =
          std::max(recomputed.residual_closure_linf_N, std::abs(closure));
      if (std::abs(closure) > maximum_closure) {
        maximum_closure = std::abs(closure);
        recomputed.maximum_closure_global_cell_id =
            report.active_global_cell_ids[row];
        recomputed.maximum_closure_component =
            static_cast<std::uint32_t>(component);
        recomputed.maximum_closure_total_value_N = total;
        recomputed.maximum_closure_parts_value_N = parts;
      }
    }
    recomputed.kinetic_energy_J +=
        0.5 * rho_ref_kg_per_m3 * volume * speed_squared;
    recomputed.centered_convective_power_W +=
        0.5 * speed_squared * divergence;
    recomputed.pressure_continuity_power_W +=
        pressure * divergence / rho_ref_kg_per_m3;
  }
  recomputed.reconstruction_power_W =
      recomputed.convective_power_W - recomputed.centered_convective_power_W;
  recomputed.pressure_adjoint_defect_W =
      recomputed.pressure_power_W + recomputed.pressure_continuity_power_W;
  recomputed.residual_closure_l2_N = std::sqrt(closure_square);
  recomputed.mass_divergence_l2_kg_per_s = std::sqrt(divergence_square);
  recomputed.stationary_wall_flux_linf_kg_per_s =
      report.stationary_wall_flux_linf_kg_per_s;
  recomputed.physical_boundary_flux_linf_kg_per_s =
      report.physical_boundary_flux_linf_kg_per_s;
  const double closure_bound =
      2048.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, closure_scale);
  return fp64_arrays_bitwise_equal(spatial_energy_scalars(recomputed),
                                   spatial_energy_scalars(report)) &&
         recomputed.maximum_closure_global_cell_id ==
             report.maximum_closure_global_cell_id &&
         recomputed.maximum_closure_component ==
             report.maximum_closure_component &&
         report.residual_closure_l2_N <= closure_bound &&
         report.residual_closure_linf_N <= closure_bound;
}

struct SpatialEnergyClassification final {
  bool pressure_adjoint{};
  bool reconstruction_dissipative{};
  bool viscous_dissipative{};
  bool implicit_viscous_reference_dissipative{};
  bool stationary_boundary_flux{};
};

SpatialEnergyClassification classify_spatial_energy(
    const flow::test::ImmersedFlowSpatialEnergyReport &report,
    double rho_ref_kg_per_m3) noexcept {
  double pressure_scale = 0.0;
  double reconstruction_scale = 0.0;
  double viscous_scale = 0.0;
  double implicit_viscous_reference_scale = 0.0;
  for (std::size_t row = 0U; row < report.active_global_cell_ids.size();
       ++row) {
    double speed_squared = 0.0;
    for (std::size_t component = 0U; component < 3U; ++component) {
      const auto offset = row * 3U + component;
      const double velocity = report.velocity_m_per_s[offset];
      speed_squared += velocity * velocity;
      pressure_scale +=
          std::abs(velocity * report.pressure_residual_N[offset]);
      reconstruction_scale +=
          std::abs(velocity * report.convective_residual_N[offset]);
      viscous_scale +=
          std::abs(velocity * report.viscous_residual_N[offset]);
      implicit_viscous_reference_scale +=
          std::abs(velocity *
                   report.implicit_viscous_reference_residual_N[offset]);
    }
    pressure_scale +=
        std::abs(report.pressure_pa[row] *
                 report.mass_divergence_kg_per_s[row] /
                 rho_ref_kg_per_m3);
    reconstruction_scale +=
        std::abs(0.5 * speed_squared *
                 report.mass_divergence_kg_per_s[row]);
  }
  const auto bound = [](double scale) {
    return 4096.0 * std::numeric_limits<double>::epsilon() *
           std::max(1.0, scale);
  };
  return {std::abs(report.pressure_adjoint_defect_W) <=
              bound(pressure_scale),
          report.reconstruction_power_W >= -bound(reconstruction_scale),
          report.viscous_power_W >= -bound(viscous_scale),
          report.implicit_viscous_reference_power_W >=
              -bound(implicit_viscous_reference_scale),
          test::fp64_bits(report.stationary_wall_flux_linf_kg_per_s) ==
                  test::fp64_bits(0.0) &&
              test::fp64_bits(report.physical_boundary_flux_linf_kg_per_s) ==
                  test::fp64_bits(0.0)};
}

void print_cell_pressure_correction_record(
    const runtime::MpiContext &mpi, runtime::Int3 global_extent,
    std::size_t corrector,
    const flow::test::ImmersedFlowCellPressureCorrectionRecord &record) {
  if (mpi.rank() != 0)
    return;
  double parity[10]{};
  double squares[6]{};
  double maxima[6]{};
  const auto nx = static_cast<std::uint64_t>(global_extent.x);
  const auto ny = static_cast<std::uint64_t>(global_extent.y);
  const auto plane = nx * ny;
  for (std::size_t row = 0U; row < record.active_global_cell_ids.size();
       ++row) {
    const auto gid = record.active_global_cell_ids[row];
    const auto i = gid % nx;
    const auto j = (gid / nx) % ny;
    const auto k = gid / plane;
    const double sign = ((i + j + k) & 1U) == 0U ? 1.0 : -1.0;
    parity[0] += sign * record.pressure_before_pa[row];
    parity[1] += sign * record.pressure_correction_pa[row];
    parity[2] += sign * record.pressure_after_pa[row];
    parity[3] += sign * record.pressure_rhs_per_volume[row];
    const double diagnostics[6]{
        record.compact_pressure_action_per_volume[row],
        record.interface_pressure_action_per_volume[row],
        record.hybrid_pressure_action_per_volume[row],
        record.affine_wall_source_per_volume[row],
        record.compact_predictor_defect_per_volume[row],
        record.hybrid_predictor_defect_per_volume[row]};
    for (std::size_t index = 0U; index < 6U; ++index) {
      parity[index + 4U] += sign * diagnostics[index];
      squares[index] += diagnostics[index] * diagnostics[index];
      maxima[index] = std::max(maxima[index], std::abs(diagnostics[index]));
    }
  }
  const double count =
      static_cast<double>(record.active_global_cell_ids.size());
  std::cerr << std::setprecision(17)
            << "cell_pressure_correction_authority corrector=" << corrector
            << " momentum_l2=" << record.momentum_l2
            << " pressure_l2=" << record.pressure_l2
            << " closure_l2=" << record.closure_l2
            << " momentum_linf=" << record.momentum_linf
            << " pressure_linf=" << record.pressure_linf
            << " closure_linf=" << record.closure_linf
            << " maximum_gid=" << record.maximum_closure_global_cell_id
            << " component=" << record.maximum_closure_component
            << " momentum=" << record.maximum_closure_momentum_value
            << " pressure=" << record.maximum_closure_pressure_value
            << " parity_before=" << parity[0] / count
            << " parity_correction=" << parity[1] / count
            << " parity_after=" << parity[2] / count
            << " rhs_parity=" << parity[3] / count
            << " compact_parity=" << parity[4] / count
            << " interface_parity=" << parity[5] / count
            << " hybrid_parity=" << parity[6] / count
            << " affine_parity=" << parity[7] / count
            << " compact_defect_parity=" << parity[8] / count
            << " hybrid_defect_parity=" << parity[9] / count
            << " compact_defect_l2=" << std::sqrt(squares[4])
            << " compact_defect_linf=" << maxima[4]
            << " hybrid_defect_l2=" << std::sqrt(squares[5])
            << " hybrid_defect_linf=" << maxima[5]
            << " wall_flux_l2=" << record.wall_predictor_mass_flux_l2_kg_per_s
            << " wall_flux_linf="
            << record.wall_predictor_mass_flux_linf_kg_per_s
            << " wall_gradient_l2="
            << record.wall_correction_gradient_l2_pa_per_m
            << " wall_gradient_linf="
            << record.wall_correction_gradient_linf_pa_per_m << '\n';
}

void print_pressure_flux_report(
    const runtime::MpiContext &mpi, const char *kind,
    const flow::test::ImmersedFlowPressureFluxIdentityReport &report) {
  if (mpi.rank() != 0)
    return;
  std::cerr << std::setprecision(17)
            << "pressure_flux_identity kind=" << kind
            << " operator_l2=" << report.operator_l2
            << " routed_l2=" << report.routed_l2
            << " difference_l2=" << report.difference_l2
            << " operator_linf=" << report.operator_linf
            << " routed_linf=" << report.routed_linf
            << " difference_linf=" << report.difference_linf
            << " maximum_gid="
            << report.maximum_difference_global_cell_id
            << " operator=" << report.maximum_difference_operator_value
            << " routed=" << report.maximum_difference_routed_value << '\n';
}

void print_final_residual_report(
    const runtime::MpiContext &mpi, const char *order,
    const flow::test::ImmersedFlowFinalMomentumPressureResidualReport &report) {
  if (mpi.rank() != 0)
    return;
  std::cerr << std::setprecision(17)
            << "final_momentum_pressure_residual order=" << order
            << " used_l2=" << report.used_l2
            << " direct_l2=" << report.direct_l2
            << " difference_l2=" << report.difference_l2
            << " used_linf=" << report.used_linf
            << " direct_linf=" << report.direct_linf
            << " difference_linf=" << report.difference_linf
            << " maximum_gid="
            << report.maximum_difference_global_cell_id
            << " component=" << report.maximum_difference_component
            << " used=" << report.maximum_difference_used_value
            << " direct=" << report.maximum_difference_direct_value << '\n';
}

void exercise_spatial_energy_diagnostics(
    const runtime::MpiContext &mpi, const Fixture &fixture,
    flow::FixedStepImmersedFlow &immersed_flow) {
  auto check_probe = [&](double phase) {
    auto probe = fixture.make_energy_probe_state(phase);
    const auto history_before = probe.snapshot(flow::FlowLayer::history);
    const auto committed_before = probe.snapshot(flow::FlowLayer::committed);
    const auto trial_before = probe.snapshot(flow::FlowLayer::trial);
    const auto metadata_before = probe.metadata();
    const std::array<std::uint64_t, 3> allocations_before{
        flow::test::MaterialDensityPisoTestAccess::state_allocation_identity(
            probe, flow::FlowLayer::history),
        flow::test::MaterialDensityPisoTestAccess::state_allocation_identity(
            probe, flow::FlowLayer::committed),
        flow::test::MaterialDensityPisoTestAccess::state_allocation_identity(
            probe, flow::FlowLayer::trial)};
    const auto revision_before =
        flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow);
    const auto report =
        flow::test::ImmersedFlowTestAccess::spatial_energy_terms(
            immersed_flow, probe, flow::FlowLayer::committed, 1.0, 0.01);
    const auto repeated =
        flow::test::ImmersedFlowTestAccess::spatial_energy_terms(
            immersed_flow, probe, flow::FlowLayer::committed, 1.0, 0.01);
    check_state_equal(probe, history_before, committed_before, trial_before,
                      metadata_before);
    const std::array<std::uint64_t, 3> allocations_after{
        flow::test::MaterialDensityPisoTestAccess::state_allocation_identity(
            probe, flow::FlowLayer::history),
        flow::test::MaterialDensityPisoTestAccess::state_allocation_identity(
            probe, flow::FlowLayer::committed),
        flow::test::MaterialDensityPisoTestAccess::state_allocation_identity(
            probe, flow::FlowLayer::trial)};
    HUNDUN_CHECK(allocations_after == allocations_before);
    HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow) ==
                 revision_before);
    HUNDUN_CHECK(spatial_energy_reports_bitwise_equal(report, repeated));
    HUNDUN_CHECK(spatial_energy_report_is_self_consistent(report, 1.0));
    auto ordinary_mutation = report;
    ordinary_mutation.velocity_m_per_s.front() += 1.0;
    HUNDUN_CHECK(
        !spatial_energy_report_is_self_consistent(ordinary_mutation, 1.0));
    auto nested_size_mutation = report;
    nested_size_mutation.cell_volume_m3.pop_back();
    HUNDUN_CHECK(
        !spatial_energy_report_is_self_consistent(nested_size_mutation, 1.0));
    auto aggregate_mutation = report;
    aggregate_mutation.total_power_W += 1.0;
    HUNDUN_CHECK(
        !spatial_energy_report_is_self_consistent(aggregate_mutation, 1.0));
    const auto exact_copy = report;
    HUNDUN_CHECK(spatial_energy_reports_bitwise_equal(report, exact_copy));
    if (mpi.rank() == 0)
      std::cerr << std::setprecision(17)
                << "spatial_energy phase=" << phase
                << " total=" << report.total_power_W
                << " convective=" << report.convective_power_W
                << " centered=" << report.centered_convective_power_W
                << " reconstruction=" << report.reconstruction_power_W
                << " pressure=" << report.pressure_power_W
                << " pressure_continuity="
                << report.pressure_continuity_power_W
                << " pressure_adjoint=" << report.pressure_adjoint_defect_W
                << " viscous=" << report.viscous_power_W
                << " implicit_viscous_reference="
                << report.implicit_viscous_reference_power_W
                << " closure_l2=" << report.residual_closure_l2_N << '\n';
    return report;
  };
  const auto first = check_probe(0.17);
  const auto second = check_probe(0.43);
  HUNDUN_CHECK(first.active_global_cell_ids ==
               second.active_global_cell_ids);
  const auto first_classification = classify_spatial_energy(first, 1.0);
  const auto second_classification = classify_spatial_energy(second, 1.0);
  const bool pressure_consistent =
      first_classification.pressure_adjoint &&
      second_classification.pressure_adjoint;
  const bool reconstruction_dissipative =
      first_classification.reconstruction_dissipative &&
      second_classification.reconstruction_dissipative;
  const bool viscous_dissipative = first_classification.viscous_dissipative &&
                                   second_classification.viscous_dissipative;
  const bool implicit_viscous_reference_dissipative =
      first_classification.implicit_viscous_reference_dissipative &&
      second_classification.implicit_viscous_reference_dissipative;
  const bool stationary_boundary_flux =
      first_classification.stationary_boundary_flux &&
      second_classification.stationary_boundary_flux;

  if (mpi.rank() == 0)
    std::cerr << std::setprecision(17)
              << "spatial_energy_classification pressure_boundary_projection="
              << (pressure_consistent ? "not_observed" : "confirmed")
              << " reconstruction_energy_production="
              << (reconstruction_dissipative ? "not_observed" : "confirmed")
              << " viscous_energy_sign="
              << (viscous_dissipative ? "not_observed" : "confirmed")
              << " implicit_viscous_reference_energy_sign="
              << (implicit_viscous_reference_dissipative ? "not_observed"
                                                         : "confirmed")
              << " stationary_boundary_flux="
              << (stationary_boundary_flux ? "not_observed" : "confirmed")
              << '\n';
  HUNDUN_CHECK(pressure_consistent);
  HUNDUN_CHECK(reconstruction_dissipative);
  HUNDUN_CHECK(viscous_dissipative);
  HUNDUN_CHECK(implicit_viscous_reference_dissipative);
  HUNDUN_CHECK(stationary_boundary_flux);
}

flow::test::ImmersedFlowFinalMomentumPressureResidualReport exercise_constraint_mode(
    const runtime::MpiContext &mpi, FixtureKind kind, bool expected_reference,
    bool checkerboard = false,
    flow::MomentumTimeOrder order = flow::MomentumTimeOrder::backward_euler,
    double dt_s = 0.01, double previous_dt_s = 0.0) {
  Fixture fixture(mpi, kind);
  auto halo = runtime::HaloExchange::create(
      fixture.decomposition,
      runtime::ExchangePlan::create(fixture.decomposition,
                                    fixture.decomposition.local_extent(),
                                    static_cast<int>(
                                        fixture.ghost_plan
                                            ->maximum_halo_reach())));
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan, nullptr,
      &fixture.transform, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);
  check_wall_effective_measures(immersed_flow, fixture);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::has_active_pressure_reference(
                   immersed_flow) == expected_reference);
  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                    1.0,
                                    0.01,
                                    std::nullopt,
                                    std::nullopt,
                                    std::nullopt};
  const auto stencil =
      flow::make_momentum_time_stencil(order, dt_s, previous_dt_s);
  auto state = fixture.make_state(checkerboard);
  linear::SolveControl pressure_control;
  if (checkerboard) {
    pressure_control.atol = 1.0e-14;
    pressure_control.rtol = 1.0e-13;
  }
  const auto history_before_attempt = state.snapshot(flow::FlowLayer::history);
  const auto committed_before_attempt =
      state.snapshot(flow::FlowLayer::committed);
  const auto trial_before_attempt = state.snapshot(flow::FlowLayer::trial);
  const auto metadata_before_attempt = state.metadata();
  const auto result =
      immersed_flow.attempt(state, physics, stencil, {}, pressure_control);
  const auto &base = std::get<flow::StepAttemptReport>(result.base);
  if (base.disposition != flow::StepAttemptDisposition::committed) {
    if (mpi.rank() == 0)
      std::cerr << std::setprecision(17)
                << "projection_attempt_rejected disposition="
                << static_cast<int>(base.disposition)
                << " reason=" << static_cast<int>(base.reason)
                << " correctors=" << base.pressure_corrector_count
                << " residuals=" << base.final_momentum_normalized_l2[0]
                << ',' << base.final_momentum_normalized_l2[1] << ','
                << base.final_momentum_normalized_l2[2] << '\n';
    HUNDUN_CHECK(base.disposition ==
                 flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(base.reason ==
                 flow::StepFailureReason::final_momentum_residual);
    HUNDUN_CHECK(base.pressure_corrector_count == 2U);
    check_state_equal(state, history_before_attempt, committed_before_attempt,
                      trial_before_attempt, metadata_before_attempt);
  }
  if (checkerboard || expected_reference) {
    const std::size_t count =
        fixture.domain->active_cells().owned_active_count();
    std::vector<double> pressure_probe(count);
    for (std::size_t row = 0U; row < count; ++row) {
      const auto id = fixture.domain->active_cells().ordered_global_ids()[row];
      pressure_probe[row] =
          std::sin(0.017 * static_cast<double>((id % 104729U) + 1U)) +
          0.25 * std::cos(0.031 * static_cast<double>((id % 65537U) + 3U));
    }
    std::vector<flow::test::ImmersedFlowWallGradientSnapshot> zero_gradients;
    zero_gradients.reserve(fixture.domain->links().size());
    for (const auto &link : fixture.domain->links())
      zero_gradients.push_back({link.id, 0.0});
    std::sort(zero_gradients.begin(), zero_gradients.end(),
              [](const auto &left, const auto &right) {
                return left.link < right.link;
              });
    zero_gradients.erase(
        std::unique(zero_gradients.begin(), zero_gradients.end(),
                    [](const auto &left, const auto &right) {
                      return left.link == right.link;
                    }),
        zero_gradients.end());
    auto nonzero_gradients = zero_gradients;
    for (auto &gradient : nonzero_gradients)
      gradient.normal_gradient_pa_per_m =
          0.125 * std::sin(0.019 * static_cast<double>(gradient.link + 1U));
    const auto check_identity = [&](const char *kind,
                                    const auto &wall_gradients) {
      const auto history_before_query =
          state.snapshot(flow::FlowLayer::history);
      const auto committed_before_query =
          state.snapshot(flow::FlowLayer::committed);
      const auto trial_before_query = state.snapshot(flow::FlowLayer::trial);
      const auto metadata_before_query = state.metadata();
      const auto revision_before_query =
          flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow);
      const auto report =
          flow::test::ImmersedFlowTestAccess::exact_predictor_schur_identity(
              immersed_flow, state, physics.rho_ref_kg_per_m3, stencil,
              pressure_probe, wall_gradients);
      check_state_equal(state, history_before_query, committed_before_query,
                        trial_before_query, metadata_before_query);
      HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow) ==
                   revision_before_query);
      const auto repeated =
          flow::test::ImmersedFlowTestAccess::exact_predictor_schur_identity(
              immersed_flow, state, physics.rho_ref_kg_per_m3, stencil,
              pressure_probe, wall_gradients);
      check_state_equal(state, history_before_query, committed_before_query,
                        trial_before_query, metadata_before_query);
      HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow) ==
                   revision_before_query);
      HUNDUN_CHECK(pressure_flux_reports_bitwise_equal(report, repeated));
      print_pressure_flux_report(mpi, kind, report);
      HUNDUN_CHECK(pressure_flux_routes_agree(report));
      auto ordinary_mutation = report;
      ordinary_mutation.routed_flux_divergence_per_volume.front() += 1.0;
      HUNDUN_CHECK(!pressure_flux_routes_agree(ordinary_mutation));
      auto nested_size_mutation = report;
      nested_size_mutation.routed_flux_divergence_per_volume.pop_back();
      HUNDUN_CHECK(!pressure_flux_routes_agree(nested_size_mutation));
      auto aggregate_mutation = report;
      aggregate_mutation.difference_l2 += 1.0;
      HUNDUN_CHECK(!pressure_flux_routes_agree(aggregate_mutation));
      return report;
    };
    const auto homogeneous = check_identity(
        expected_reference ? "pressure_reference_homogeneous"
                           : "closed_homogeneous",
        zero_gradients);
    const auto affine = check_identity(
        expected_reference ? "pressure_reference_affine" : "closed_affine",
        nonzero_gradients);

    std::vector<double> second_pressure_probe(count);
    for (std::size_t row = 0U; row < count; ++row) {
      const auto id = fixture.domain->active_cells().ordered_global_ids()[row];
      second_pressure_probe[row] =
          std::cos(0.013 * static_cast<double>((id % 99991U) + 5U)) -
          0.375 * std::sin(0.029 * static_cast<double>((id % 49157U) + 7U));
    }
    const auto workspace_before =
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
    HUNDUN_CHECK(workspace_before.probe_state_creation_count == 1U);
    const auto first_action =
        flow::test::ImmersedFlowTestAccess::active_pressure_operator_values(
            immersed_flow, pressure_probe, zero_gradients);
    const auto second_action =
        flow::test::ImmersedFlowTestAccess::active_pressure_operator_values(
            immersed_flow, second_pressure_probe, zero_gradients);
    const auto workspace_after =
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
    HUNDUN_CHECK(workspace_after.probe_state_creation_count ==
                 workspace_before.probe_state_creation_count);
    HUNDUN_CHECK(workspace_after.response_capacity_growth_count ==
                 workspace_before.response_capacity_growth_count);
    HUNDUN_CHECK(first_action.size() == count);
    HUNDUN_CHECK(second_action.size() == count);
    double bilinear[3]{};
    for (std::size_t row = 0U; row < count; ++row) {
      const auto local = fixture.topology.find_local_cell(
          fixture.domain->active_cells().ordered_global_ids()[row]);
      HUNDUN_CHECK(local.has_value());
      const double volume = fixture.geometry.cell_volume_m3(*local);
      bilinear[0] += volume * pressure_probe[row] * second_action[row];
      bilinear[1] += volume * second_pressure_probe[row] * first_action[row];
      bilinear[2] +=
          volume * (std::abs(pressure_probe[row] * second_action[row]) +
                    std::abs(second_pressure_probe[row] * first_action[row]));
    }
    mpi.allreduce_fp64_in_place(bilinear, 3U,
                                runtime::Fp64ReductionOperation::sum);
    const double symmetry_tolerance =
        4096.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, bilinear[2]);
    if (mpi.rank() == 0)
      std::cerr << std::setprecision(17)
                << "exact_predictor_schur_symmetry xAy=" << bilinear[0]
                << " yAx=" << bilinear[1]
                << " defect=" << bilinear[0] - bilinear[1]
                << " tolerance=" << symmetry_tolerance << '\n';
    HUNDUN_CHECK(std::abs(bilinear[0] - bilinear[1]) > symmetry_tolerance);

    HUNDUN_CHECK(homogeneous.operator_residual_per_volume.size() ==
                 affine.operator_residual_per_volume.size());
    HUNDUN_CHECK(homogeneous.routed_flux_divergence_per_volume.size() ==
                 affine.routed_flux_divergence_per_volume.size());
    double affine_operator_response = 0.0;
    double affine_routed_response = 0.0;
    for (std::size_t row = 0U;
         row < homogeneous.operator_residual_per_volume.size(); ++row) {
      affine_operator_response = std::max(
          affine_operator_response,
          std::abs(affine.operator_residual_per_volume[row] -
                   homogeneous.operator_residual_per_volume[row]));
      affine_routed_response = std::max(
          affine_routed_response,
          std::abs(affine.routed_flux_divergence_per_volume[row] -
                   homogeneous.routed_flux_divergence_per_volume[row]));
    }
    mpi.allreduce_fp64_in_place(&affine_operator_response, 1U,
                                runtime::Fp64ReductionOperation::maximum);
    mpi.allreduce_fp64_in_place(&affine_routed_response, 1U,
                                runtime::Fp64ReductionOperation::maximum);
    HUNDUN_CHECK(
        affine_operator_response >
        256.0 * std::numeric_limits<double>::epsilon());
    HUNDUN_CHECK_NEAR(affine_routed_response, affine_operator_response,
                      64.0 * std::numeric_limits<double>::epsilon() *
                          std::max(1.0, affine_operator_response));
  }
  {
    const auto history_before_query = state.snapshot(flow::FlowLayer::history);
    const auto committed_before_query =
        state.snapshot(flow::FlowLayer::committed);
    const auto trial_before_query = state.snapshot(flow::FlowLayer::trial);
    const auto metadata_before_query = state.metadata();
    const auto pressure_revision_before_query =
        flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow);
    const auto correction_report =
        flow::test::ImmersedFlowTestAccess::cell_pressure_correction_authority(
            immersed_flow);
    check_state_equal(state, history_before_query, committed_before_query,
                      trial_before_query, metadata_before_query);
    HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow) ==
                 pressure_revision_before_query);
    const auto repeated =
        flow::test::ImmersedFlowTestAccess::cell_pressure_correction_authority(
            immersed_flow);
    HUNDUN_CHECK(cell_pressure_correction_reports_bitwise_equal(
        correction_report, repeated));
    HUNDUN_CHECK(correction_report.correctors.size() == 2U);
    HUNDUN_CHECK(
        test::fp64_bits(
            correction_report.inter_corrector_authority_difference_l2) ==
        test::fp64_bits(0.0));
    for (std::size_t corrector = 0U;
         corrector < correction_report.correctors.size(); ++corrector) {
      print_cell_pressure_correction_record(
          mpi, fixture.topology.global_extent(), corrector + 1U,
          correction_report.correctors[corrector]);
      HUNDUN_CHECK(cell_pressure_correction_routes_agree(
          correction_report.correctors[corrector]));
    }
    if (checkerboard) {
      auto ordinary_mutation = correction_report;
      ordinary_mutation.correctors.front()
          .momentum_operator_velocity_change.front() += 1.0;
      HUNDUN_CHECK(!cell_pressure_correction_routes_agree(
          ordinary_mutation.correctors.front()));
      auto nested_size_mutation = correction_report;
      nested_size_mutation.correctors.front()
          .lfp_pressure_residual_change.pop_back();
      HUNDUN_CHECK(!cell_pressure_correction_routes_agree(
          nested_size_mutation.correctors.front()));
      auto aggregate_mutation = correction_report;
      aggregate_mutation.correctors.front().closure_l2 += 1.0;
      HUNDUN_CHECK(!cell_pressure_correction_routes_agree(
          aggregate_mutation.correctors.front()));
      auto diagnostic_nested_mutation = correction_report;
      diagnostic_nested_mutation.correctors.front()
          .pressure_after_pa.pop_back();
      HUNDUN_CHECK(!cell_pressure_correction_routes_agree(
          diagnostic_nested_mutation.correctors.front()));
      auto diagnostic_scalar_mutation = correction_report;
      diagnostic_scalar_mutation.correctors.front()
          .wall_predictor_mass_flux_l2_kg_per_s += 1.0;
      HUNDUN_CHECK(!cell_pressure_correction_reports_bitwise_equal(
          correction_report, diagnostic_scalar_mutation));
      auto predictor_value_mutation = correction_report;
      predictor_value_mutation.correctors.front()
          .compact_pressure_action_per_volume.front() += 1.0;
      HUNDUN_CHECK(!cell_pressure_correction_routes_agree(
          predictor_value_mutation.correctors.front()));
      std::size_t interface_mutation_index = 0U;
      double interface_mutation_magnitude = 0.0;
      for (std::size_t row = 0U;
           row < correction_report.correctors.front()
                     .interface_pressure_action_per_volume.size();
           ++row) {
        const double magnitude = std::abs(
            correction_report.correctors.front()
                .interface_pressure_action_per_volume[row]);
        if (magnitude > interface_mutation_magnitude) {
          interface_mutation_magnitude = magnitude;
          interface_mutation_index = row;
        }
      }
      HUNDUN_CHECK(interface_mutation_magnitude >
                   256.0 * std::numeric_limits<double>::epsilon());
      auto predictor_sign_mutation = correction_report;
      predictor_sign_mutation.correctors.front()
          .interface_pressure_action_per_volume[interface_mutation_index] *=
          -1.0;
      HUNDUN_CHECK(!cell_pressure_correction_routes_agree(
          predictor_sign_mutation.correctors.front()));
      auto predictor_duplicate_mutation = correction_report;
      predictor_duplicate_mutation.correctors.front()
          .interface_pressure_action_per_volume[interface_mutation_index] *=
          2.0;
      HUNDUN_CHECK(!cell_pressure_correction_routes_agree(
          predictor_duplicate_mutation.correctors.front()));
      auto predictor_nested_size_mutation = correction_report;
      predictor_nested_size_mutation.correctors.front()
          .hybrid_predictor_defect_per_volume.pop_back();
      HUNDUN_CHECK(!cell_pressure_correction_routes_agree(
          predictor_nested_size_mutation.correctors.front()));
      auto authority_mutation = correction_report;
      authority_mutation.inter_corrector_authority_difference_l2 += 1.0;
      HUNDUN_CHECK(!cell_pressure_correction_reports_bitwise_equal(
          correction_report, authority_mutation));
    }
  }
  if (checkerboard)
    exercise_spatial_energy_diagnostics(mpi, fixture, immersed_flow);
  HUNDUN_CHECK(base.disposition == flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(base.pressure_corrector_count == 2U);
  HUNDUN_CHECK(base.final_continuity_normalized_l2 <= 1.0e-10);
  for (const double value : base.final_momentum_normalized_l2)
    HUNDUN_CHECK(value <= 1.0e-9);
  const auto history_before_query = state.snapshot(flow::FlowLayer::history);
  const auto committed_before_query =
      state.snapshot(flow::FlowLayer::committed);
  const auto trial_before_query = state.snapshot(flow::FlowLayer::trial);
  const auto metadata_before_query = state.metadata();
  const auto pressure_revision_before_query =
      flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow);
  const auto final_residual =
      flow::test::ImmersedFlowTestAccess::final_momentum_pressure_residual_routes(
          immersed_flow, state);
  check_state_equal(state, history_before_query, committed_before_query,
                    trial_before_query, metadata_before_query);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow) ==
               pressure_revision_before_query);
  const auto repeated =
      flow::test::ImmersedFlowTestAccess::final_momentum_pressure_residual_routes(
          immersed_flow, state);
  check_state_equal(state, history_before_query, committed_before_query,
                    trial_before_query, metadata_before_query);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow) ==
               pressure_revision_before_query);
  HUNDUN_CHECK(final_residual_reports_bitwise_equal(final_residual, repeated));
  print_final_residual_report(
      mpi, order == flow::MomentumTimeOrder::backward_euler ? "BE" : "BDF2",
      final_residual);
  if (checkerboard) {
    const auto values = state.snapshot(flow::FlowLayer::committed);
    double parity[2]{};
    for (mesh::LocalCellId cell = 0U;
         cell < fixture.topology.owned_cell_count(); ++cell) {
      if (fixture.domain->region(cell) != immersed::CellRegion::fluid)
        continue;
      const auto index = fixture.topology.global_cell(cell);
      const double sign = ((index.x + index.y + index.z) % 2 == 0) ? 1.0 : -1.0;
      parity[0] += sign * values.mechanical_pressure[cell];
      parity[1] += 1.0;
    }
    mpi.allreduce_fp64_in_place(parity, 2U,
                                runtime::Fp64ReductionOperation::sum);
    const double parity_amplitude = std::abs(parity[0] / parity[1]);
    std::vector<mesh::GlobalCellId> wall_band_ids;
    wall_band_ids.reserve(fixture.domain->links().size());
    for (const auto &link : fixture.domain->links())
      wall_band_ids.push_back(link.fluid_cell);
    std::sort(wall_band_ids.begin(), wall_band_ids.end());
    wall_band_ids.erase(std::unique(wall_band_ids.begin(), wall_band_ids.end()),
                        wall_band_ids.end());
    double band_sums[6]{};
    double band_maxima[2]{};
    const auto owned_active_count =
        fixture.domain->active_cells().owned_active_count();
    std::vector<double> final_active_pressure(owned_active_count);
    for (std::size_t row = 0U; row < owned_active_count; ++row) {
      const auto gid = fixture.domain->active_cells().ordered_global_ids()[row];
      const auto local = fixture.topology.find_local_cell(gid);
      HUNDUN_CHECK(local.has_value());
      const auto logical = fixture.topology.global_cell(*local);
      const double sign =
          ((logical.x + logical.y + logical.z) % 2 == 0) ? 1.0 : -1.0;
      const double pressure = values.mechanical_pressure[*local];
      final_active_pressure[row] = pressure;
      const std::size_t band =
          std::binary_search(wall_band_ids.begin(), wall_band_ids.end(), gid)
              ? 0U
              : 1U;
      band_sums[band * 3U] += sign * pressure;
      band_sums[band * 3U + 1U] += 1.0;
      band_sums[band * 3U + 2U] += pressure * pressure;
      band_maxima[band] = std::max(band_maxima[band], std::abs(pressure));
    }
    mpi.allreduce_fp64_in_place(band_sums, 6U,
                                runtime::Fp64ReductionOperation::sum);
    mpi.allreduce_fp64_in_place(band_maxima, 2U,
                                runtime::Fp64ReductionOperation::maximum);
    std::vector<flow::test::ImmersedFlowWallGradientSnapshot> zero_gradients;
    zero_gradients.reserve(fixture.domain->links().size());
    for (const auto &link : fixture.domain->links())
      zero_gradients.push_back({link.id, 0.0});
    std::sort(zero_gradients.begin(), zero_gradients.end(),
              [](const auto &left, const auto &right) {
                return left.link < right.link;
              });
    zero_gradients.erase(std::unique(zero_gradients.begin(),
                                     zero_gradients.end(),
                                     [](const auto &left, const auto &right) {
                                       return left.link == right.link;
                                     }),
                         zero_gradients.end());
    const auto operator_values =
        flow::test::ImmersedFlowTestAccess::active_pressure_operator_values(
            immersed_flow, final_active_pressure, zero_gradients);
    HUNDUN_CHECK(operator_values.size() == final_active_pressure.size());
    double operator_evidence[3]{};
    for (std::size_t row = 0U; row < operator_values.size(); ++row) {
      operator_evidence[0] += final_active_pressure[row] * operator_values[row];
      operator_evidence[1] += operator_values[row] * operator_values[row];
      operator_evidence[2] +=
          final_active_pressure[row] * final_active_pressure[row];
    }
    mpi.allreduce_fp64_in_place(operator_evidence, 3U,
                                runtime::Fp64ReductionOperation::sum);
    if (mpi.rank() == 0) {
      const auto amplitude = [&](std::size_t band) {
        return std::abs(band_sums[band * 3U] / band_sums[band * 3U + 1U]);
      };
      std::cerr << std::setprecision(17)
                << "immersed_checkerboard parity=" << parity_amplitude
                << " wall_band_parity=" << amplitude(0U)
                << " bulk_parity=" << amplitude(1U)
                << " wall_band_l2=" << std::sqrt(band_sums[2U])
                << " bulk_l2=" << std::sqrt(band_sums[5U])
                << " wall_band_linf=" << band_maxima[0]
                << " bulk_linf=" << band_maxima[1]
                << " pressure_operator_energy=" << operator_evidence[0]
                << " pressure_operator_l2=" << std::sqrt(operator_evidence[1])
                << " pressure_l2=" << std::sqrt(operator_evidence[2]) << '\n';
    }
    HUNDUN_CHECK(parity_amplitude <= 1.0e-8);
  }
  return final_residual;
}

void exercise_pressure_outlet_face_history(const runtime::MpiContext &mpi) {
  Fixture fixture(mpi, FixtureKind::outside_open);
  auto halo = runtime::HaloExchange::create(
      fixture.decomposition,
      runtime::ExchangePlan::create(
          fixture.decomposition, fixture.decomposition.local_extent(),
          static_cast<int>(fixture.ghost_plan->maximum_halo_reach())));
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan, nullptr,
      &fixture.transform, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::bdf2, 0.01, 0.02);
  constexpr double diagonal = 1.0;

  const auto evaluate = [&](OutletFaceHistory face_history) {
    auto state = fixture.make_state(false, false, false, face_history);
    const auto values =
        flow::test::ImmersedFlowTestAccess::physical_outlet_predictor_values(
            immersed_flow, state, 1.0, stencil, diagonal);
    for (const auto &value : values) {
      const auto face = fixture.topology.find_local_face(value.global_face_id);
      HUNDUN_CHECK(face.has_value());
      const auto owner = fixture.topology.owner(*face);
      const double mobility =
          fixture.geometry.cell_volume_m3(owner) / diagonal;
      const double expected_x =
          mobility / stencil.dt_s *
          ((-stencil.alpha1) * face_history.committed_velocity_x_m_per_s +
           (-stencil.alpha2) * face_history.history_velocity_x_m_per_s);
      const double tolerance =
          64.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, std::abs(expected_x));
      HUNDUN_CHECK_NEAR(value.velocity_m_per_s.x, expected_x, tolerance);
      HUNDUN_CHECK(value.velocity_m_per_s.y == 0.0);
      HUNDUN_CHECK(value.velocity_m_per_s.z == 0.0);
      const auto area =
          fixture.geometry.face_area_vector_m2(*face, mesh::FaceSide::owner);
      HUNDUN_CHECK_NEAR(value.mass_flux_kg_per_s, expected_x * area.x,
                        tolerance);
    }
    return values;
  };

  const auto baseline = evaluate({0.375, -0.125});
  const auto current_mutation = evaluate({0.5, -0.125});
  const auto nested_history_mutation = evaluate({0.375, 0.25});
  std::uint64_t global_outlet_faces = baseline.size();
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &global_outlet_faces, 1,
                             MPI_UINT64_T, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_outlet_faces > 0U);
  HUNDUN_CHECK(baseline.size() == current_mutation.size());
  HUNDUN_CHECK(baseline.size() == nested_history_mutation.size());
  for (std::size_t index = 0U; index < baseline.size(); ++index) {
    HUNDUN_CHECK(baseline[index].global_face_id ==
                 current_mutation[index].global_face_id);
    HUNDUN_CHECK(baseline[index].global_face_id ==
                 nested_history_mutation[index].global_face_id);
    HUNDUN_CHECK(baseline[index].velocity_m_per_s.x !=
                 current_mutation[index].velocity_m_per_s.x);
    HUNDUN_CHECK(baseline[index].velocity_m_per_s.x !=
                 nested_history_mutation[index].velocity_m_per_s.x);
  }
}

void exercise_immersed_wall_face_history(const runtime::MpiContext &mpi) {
  Fixture fixture(mpi, FixtureKind::outside_tetrahedron);
  auto halo = runtime::HaloExchange::create(
      fixture.decomposition,
      runtime::ExchangePlan::create(
          fixture.decomposition, fixture.decomposition.local_extent(),
          static_cast<int>(fixture.ghost_plan->maximum_halo_reach())));
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan, nullptr,
      &fixture.transform, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::bdf2, 0.01, 0.02);
  const double volume = fixture.geometry.cell_volume_m3(0U);
  const double diagonal = stencil.alpha0 * volume / stencil.dt_s;
  const auto measures =
      flow::test::ImmersedFlowTestAccess::wall_effective_measures(immersed_flow);
  const auto operator_rows = finite_volume::test::ImmersedOperatorTestAccess::
      rows(flow::test::ImmersedFlowTestAccess::immersed_operator(immersed_flow));

  struct SharedRowReconstruction final {
    mesh::GlobalCellId fluid_cell{};
    immersed::detail::QuadraticFrame frame;
    std::vector<mesh::GlobalCellId> donor_ids;
    immersed::detail::DeterministicQr factorization;
  };
  std::vector<SharedRowReconstruction> shared_rows;
  std::vector<mesh::GlobalCellId> row_ids;
  for (const auto &row : operator_rows)
    if (!row.links.empty())
      row_ids.push_back(row.active_cell);
  std::sort(row_ids.begin(), row_ids.end());
  row_ids.erase(std::unique(row_ids.begin(), row_ids.end()), row_ids.end());
  const auto global_extent = fixture.topology.global_extent();
  const auto logical_cell = [&](mesh::GlobalCellId id) {
    const auto nx = static_cast<std::uint64_t>(global_extent.x);
    const auto ny = static_cast<std::uint64_t>(global_extent.y);
    const auto plane = nx * ny;
    return runtime::Int3{static_cast<int>(id % nx),
                         static_cast<int>((id / nx) % ny),
                         static_cast<int>(id / plane)};
  };
  const immersed::detail::PeriodicCellMapper periodic_cells(
      fixture.topology.global_extent(), fixture.topology.periodicity(),
      fixture.geometry.length_m());
  std::uint64_t periodic_oracle_image_count = 0U;
  for (const auto row_id : row_ids) {
    const auto local = fixture.topology.find_local_cell(row_id);
    HUNDUN_CHECK(local.has_value());
    std::vector<mesh::GlobalCellId> donor_ids;
    const auto operator_row =
        std::find_if(operator_rows.begin(), operator_rows.end(),
                     [&](const auto &row) { return row.active_cell == row_id; });
    HUNDUN_CHECK(operator_row != operator_rows.end());
    for (const auto &link : operator_row->links) {
      const auto &link_donors =
          immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
              fixture.ghost_plan->reconstruction(link.id));
      donor_ids.insert(donor_ids.end(), link_donors.begin(),
                       link_donors.end());
    }
    for (mesh::LocalFaceId face = 0U;
         face < fixture.topology.local_face_count(); ++face) {
      const auto neighbour = fixture.topology.neighbour(face);
      if (!neighbour.has_value() ||
          fixture.domain->region(fixture.topology.owner(face)) !=
              immersed::CellRegion::fluid ||
          fixture.domain->region(*neighbour) != immersed::CellRegion::fluid ||
          (fixture.topology.owner(face) != *local && *neighbour != *local))
        continue;
      donor_ids.push_back(
          fixture.topology.global_cell_id(fixture.topology.owner(face)));
      donor_ids.push_back(fixture.topology.global_cell_id(*neighbour));
    }
    std::sort(donor_ids.begin(), donor_ids.end());
    donor_ids.erase(std::unique(donor_ids.begin(), donor_ids.end()),
                    donor_ids.end());
    const double scale =
        std::cbrt(fixture.geometry.cell_volume_m3(*local));
    const immersed::detail::QuadraticFrame frame{
        fixture.geometry.cell_center_m(*local), {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, scale};
    std::vector<double> design(
        donor_ids.size() * immersed::detail::kQuadraticBasisSize, 0.0);
    const auto row_logical = logical_cell(row_id);
    for (std::size_t donor = 0U; donor < donor_ids.size(); ++donor) {
      const auto image =
          periodic_cells.nearest_image(logical_cell(donor_ids[donor]),
                                       row_logical);
      if (image.image.x != image.canonical.x ||
          image.image.y != image.canonical.y ||
          image.image.z != image.canonical.z)
        ++periodic_oracle_image_count;
      auto image_frame = frame;
      image_frame.origin_m = {
          frame.origin_m.x - image.shift_m.x,
          frame.origin_m.y - image.shift_m.y,
          frame.origin_m.z - image.shift_m.z};
      const auto moments = immersed::detail::quadratic_cell_average_basis(
          image.canonical, image_frame, fixture.topology,
          fixture.geometry);
      std::copy(moments.begin(), moments.end(),
                design.begin() + static_cast<std::ptrdiff_t>(
                                     donor * moments.size()));
    }
    shared_rows.push_back(
        {row_id, frame, std::move(donor_ids),
         immersed::detail::factorize_design_matrix(
             design, design.size() / immersed::detail::kQuadraticBasisSize,
             immersed::detail::kQuadraticBasisSize)});
  }
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &periodic_oracle_image_count, 1,
                             MPI_UINT64_T, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(periodic_oracle_image_count > 0U);

  std::uint64_t maximum_links_per_row = 0U;
  std::uint64_t multi_link_rows = 0U;
  std::uint64_t maximum_donor_symmetric_difference = 0U;
  for (const auto row_id : row_ids) {
    const auto operator_row =
        std::find_if(operator_rows.begin(), operator_rows.end(),
                     [&](const auto &row) { return row.active_cell == row_id; });
    HUNDUN_CHECK(operator_row != operator_rows.end());
    std::vector<const finite_volume::test::ImmersedWallLinkSnapshot *>
        row_links;
    for (const auto &link : operator_row->links)
      row_links.push_back(&link);
    maximum_links_per_row =
        std::max(maximum_links_per_row,
                 static_cast<std::uint64_t>(row_links.size()));
    if (row_links.size() <= 1U)
      continue;
    ++multi_link_rows;
    for (std::size_t left = 0U; left < row_links.size(); ++left) {
      const auto &left_donors =
          immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
              fixture.ghost_plan->reconstruction(row_links[left]->id));
      for (std::size_t right = left + 1U; right < row_links.size(); ++right) {
        const auto &right_donors =
            immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
                fixture.ghost_plan->reconstruction(row_links[right]->id));
        std::vector<mesh::GlobalCellId> difference;
        std::set_symmetric_difference(
            left_donors.begin(), left_donors.end(), right_donors.begin(),
            right_donors.end(), std::back_inserter(difference));
        maximum_donor_symmetric_difference =
            std::max(maximum_donor_symmetric_difference,
                     static_cast<std::uint64_t>(difference.size()));
      }
    }
  }
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &maximum_links_per_row, 1,
                             MPI_UINT64_T, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &multi_link_rows, 1, MPI_UINT64_T,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE,
                             &maximum_donor_symmetric_difference, 1,
                             MPI_UINT64_T, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  if (mpi.rank() == 0)
    std::cerr << "shared-row diagnostic: rows=" << row_ids.size()
              << " multi_link_rows=" << multi_link_rows
              << " max_links_per_row=" << maximum_links_per_row
              << " max_donor_symmetric_difference="
              << maximum_donor_symmetric_difference << '\n';

  const auto evaluate = [&](WallCellHistory history_scales,
                            bool require_independent_link_mutation) {
    auto state = fixture.make_state(false, false, false, std::nullopt,
                                    history_scales);
    const auto values =
        flow::test::ImmersedFlowTestAccess::wall_predictor_values(
            immersed_flow, state, 1.0, stencil, diagonal);
    const auto committed =
        state.layer(flow::FlowLayer::committed).view<double>(
            fixture.fields.velocity);
    const auto history = state.layer(flow::FlowLayer::history).view<double>(
        fixture.fields.velocity);
    double maximum_owner_mismatch = 0.0;
    double maximum_independent_link_mismatch = 0.0;
    for (const auto &value : values) {
      const finite_volume::test::ImmersedWallLinkSnapshot *link = nullptr;
      mesh::GlobalCellId fluid_global{};
      for (const auto &row : operator_rows) {
        const auto found = std::find_if(
            row.links.begin(), row.links.end(), [&](const auto &candidate) {
              return candidate.id == value.link;
            });
        if (found == row.links.end())
          continue;
        link = &*found;
        fluid_global = row.active_cell;
        break;
      }
      const auto measure =
          std::find_if(measures.begin(), measures.end(), [&](const auto &item) {
            return item.link == value.link;
          });
      HUNDUN_CHECK(link != nullptr);
      HUNDUN_CHECK(measure != measures.end());
      const auto fluid = fixture.topology.find_local_cell(fluid_global);
      HUNDUN_CHECK(fluid.has_value());
      const auto shared =
          std::find_if(shared_rows.begin(), shared_rows.end(),
                       [&](const auto &row) {
                         return row.fluid_cell == fluid_global;
                       });
      HUNDUN_CHECK(shared != shared_rows.end());
      const auto shared_value = [&](runtime::Real3 point,
                                    const runtime::FieldView<const double> &field,
                                    std::size_t component) {
        const auto functional =
            immersed::detail::quadratic_basis_at(point, shared->frame);
        const auto weights = shared->factorization.functional_weights(
            std::vector<double>(functional.begin(), functional.end()));
        HUNDUN_CHECK(weights.size() == shared->donor_ids.size());
        const auto box = fixture.topology.owned_global_box();
        const auto local_extent = field.interior_extent();
        const auto global_extent = fixture.topology.global_extent();
        const int ghost_width = field.ghost_width();
        double value = 0.0;
        for (std::size_t donor = 0U; donor < weights.size(); ++donor) {
          const auto logical = logical_cell(shared->donor_ids[donor]);
          const int i = shared_row_oracle_field_coordinate(
              logical.x, box.begin.x, local_extent.x, global_extent.x,
              ghost_width);
          const int j = shared_row_oracle_field_coordinate(
              logical.y, box.begin.y, local_extent.y, global_extent.y,
              ghost_width);
          const int k = shared_row_oracle_field_coordinate(
              logical.z, box.begin.z, local_extent.z, global_extent.z,
              ghost_width);
          value += weights[donor] *
                   field(i, j, k, static_cast<int>(component));
        }
        return value;
      };
      runtime::Real3 expected_velocity{};
      runtime::Real3 wrong_owner_velocity{};
      runtime::Real3 wrong_independent_link_velocity{};
      for (std::size_t component = 0U; component < 3U; ++component) {
        const double shared_n =
            shared_value(link->wall_intercept_m, committed, component);
        const double shared_nm1 =
            shared_value(link->wall_intercept_m, history, component);
        const double independent_n =
            fixture.ghost_plan->reconstruction(link->id)
                .value(link->wall_intercept_m, committed, component);
        const double independent_nm1 =
            fixture.ghost_plan->reconstruction(link->id)
                .value(link->wall_intercept_m, history, component);
        const double factor = volume / diagonal / stencil.dt_s;
        const double expected =
            shared_n + factor * (stencil.alpha1 * shared_n +
                                 stencil.alpha2 * shared_nm1);
        const double wrong_independent =
            independent_n + factor * (stencil.alpha1 * independent_n +
                                      stencil.alpha2 * independent_nm1);
        const auto global = fixture.topology.global_cell(*fluid);
        const auto box = fixture.topology.owned_global_box();
        const runtime::Int3 index{global.x - box.begin.x,
                                 global.y - box.begin.y,
                                 global.z - box.begin.z};
        const double owner_n =
            committed(index.x, index.y, index.z, static_cast<int>(component));
        const double owner_nm1 =
            history(index.x, index.y, index.z, static_cast<int>(component));
        const double wrong_owner =
            shared_n +
            factor * (stencil.alpha1 * owner_n +
                      stencil.alpha2 * owner_nm1);
        if (component == 0U) {
          expected_velocity.x = expected;
          wrong_owner_velocity.x = wrong_owner;
          wrong_independent_link_velocity.x = wrong_independent;
        } else if (component == 1U) {
          expected_velocity.y = expected;
          wrong_owner_velocity.y = wrong_owner;
          wrong_independent_link_velocity.y = wrong_independent;
        } else {
          expected_velocity.z = expected;
          wrong_owner_velocity.z = wrong_owner;
          wrong_independent_link_velocity.z = wrong_independent;
        }
      }
      const double expected_flux =
          measure->effective_measure_m2 *
          (expected_velocity.x * link->solid_to_fluid_normal.x +
           expected_velocity.y * link->solid_to_fluid_normal.y +
           expected_velocity.z * link->solid_to_fluid_normal.z);
      const double wrong_owner_flux =
          measure->effective_measure_m2 *
          (wrong_owner_velocity.x * link->solid_to_fluid_normal.x +
           wrong_owner_velocity.y * link->solid_to_fluid_normal.y +
           wrong_owner_velocity.z * link->solid_to_fluid_normal.z);
      const double wrong_independent_link_flux =
          measure->effective_measure_m2 *
          (wrong_independent_link_velocity.x *
               link->solid_to_fluid_normal.x +
           wrong_independent_link_velocity.y *
               link->solid_to_fluid_normal.y +
           wrong_independent_link_velocity.z *
               link->solid_to_fluid_normal.z);
      const double tolerance =
          512.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, std::abs(expected_flux));
      if (std::abs(value.predictor_mass_flux_kg_per_s - expected_flux) >
          tolerance) {
        const auto &product_authority =
            finite_volume::detail::ImmersedBoundaryAuthorityAccess::
                row_reconstruction(
                    flow::test::ImmersedFlowTestAccess::immersed_operator(immersed_flow),
                    value.link);
        const auto &product_donors =
            immersed::detail::QuadraticReconstructionWeights::
                donor_global_ids(product_authority);
        std::vector<mesh::GlobalCellId> donor_difference;
        std::set_symmetric_difference(
            shared->donor_ids.begin(), shared->donor_ids.end(),
            product_donors.begin(), product_donors.end(),
            std::back_inserter(donor_difference));
        std::cerr << "shared-row mismatch: rank=" << mpi.rank()
                  << " row=" << fluid_global << " link=" << value.link
                  << " product=" << std::setprecision(17)
                  << value.predictor_mass_flux_kg_per_s
                  << " expected=" << expected_flux
                  << " difference="
                  << value.predictor_mass_flux_kg_per_s - expected_flux
                  << " tolerance=" << tolerance
                  << " oracle_donors=" << shared->donor_ids.size()
                  << " product_donors=" << product_donors.size()
                  << " donor_symmetric_difference=" << donor_difference.size()
                  << '\n';
        HUNDUN_CHECK(false);
      }
      maximum_owner_mismatch =
          std::max(maximum_owner_mismatch,
                   std::abs(wrong_owner_flux - expected_flux));
      maximum_independent_link_mismatch =
          std::max(maximum_independent_link_mismatch,
                   std::abs(wrong_independent_link_flux - expected_flux));
    }
    mpi.allreduce_fp64_in_place(&maximum_owner_mismatch, 1U,
                                runtime::Fp64ReductionOperation::maximum);
    mpi.allreduce_fp64_in_place(
        &maximum_independent_link_mismatch, 1U,
        runtime::Fp64ReductionOperation::maximum);
    if (mpi.rank() == 0)
      std::cerr << "shared-row mutation diagnostic: owner="
                << maximum_owner_mismatch << " independent_link="
                << maximum_independent_link_mismatch << '\n';
    HUNDUN_CHECK(maximum_owner_mismatch > 1.0e-8);
    if (require_independent_link_mutation)
      HUNDUN_CHECK(maximum_independent_link_mismatch > 1.0e-10);
    std::uint64_t lookup_probe_count = 0U;
    if (!values.empty()) {
      static_cast<void>(
          finite_volume::detail::ImmersedBoundaryAuthorityAccess::
              row_reconstruction(
                  flow::test::ImmersedFlowTestAccess::immersed_operator(immersed_flow),
                  values.front().link));
      lookup_probe_count =
          finite_volume::test::ImmersedOperatorTestAccess::
              last_boundary_authority_lookup_probe_count(
                  flow::test::ImmersedFlowTestAccess::immersed_operator(immersed_flow));
    }
    HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &lookup_probe_count, 1,
                               MPI_UINT64_T, MPI_MAX, mpi.comm()) ==
                 MPI_SUCCESS);
    HUNDUN_CHECK(lookup_probe_count == 1U);
    bool missing_link_rejected = false;
    try {
      static_cast<void>(
          finite_volume::detail::ImmersedBoundaryAuthorityAccess::
              row_reconstruction(
                  flow::test::ImmersedFlowTestAccess::immersed_operator(immersed_flow),
                  std::numeric_limits<immersed::ImmersedLinkId>::max()));
    } catch (const runtime::Error &error) {
      missing_link_rejected =
          std::string(error.what()).find("outside the operator plan") !=
          std::string::npos;
    }
    HUNDUN_CHECK(missing_link_rejected);
    HUNDUN_CHECK(
        finite_volume::test::ImmersedOperatorTestAccess::
            last_boundary_authority_lookup_probe_count(
                flow::test::ImmersedFlowTestAccess::immersed_operator(immersed_flow)) ==
        0U);
    return values;
  };

  const auto current_mutation = evaluate({1.25, 1.0}, true);
  const auto steady = evaluate({1.0, 1.0}, false);
  const auto nested_history_mutation = evaluate({1.0, 0.5}, true);
  HUNDUN_CHECK(steady.size() == current_mutation.size());
  HUNDUN_CHECK(steady.size() == nested_history_mutation.size());
  double current_change = 0.0;
  double nested_change = 0.0;
  for (std::size_t index = 0U; index < steady.size(); ++index) {
    HUNDUN_CHECK(steady[index].link == current_mutation[index].link);
    HUNDUN_CHECK(steady[index].link == nested_history_mutation[index].link);
    current_change =
        std::max(current_change,
                 std::abs(steady[index].predictor_mass_flux_kg_per_s -
                          current_mutation[index].predictor_mass_flux_kg_per_s));
    nested_change = std::max(
        nested_change,
        std::abs(steady[index].predictor_mass_flux_kg_per_s -
                 nested_history_mutation[index].predictor_mass_flux_kg_per_s));
  }
  mpi.allreduce_fp64_in_place(&current_change, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  mpi.allreduce_fp64_in_place(&nested_change, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  HUNDUN_CHECK(current_change > 1.0e-8);
  HUNDUN_CHECK(nested_change > 1.0e-8);
}

void exercise_underprovisioned_wall_halo_rejection(
    const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(shared_row_oracle_field_coordinate(3, 4, 4, 12, 1) == -1);
  HUNDUN_CHECK(shared_row_oracle_field_coordinate(11, 0, 4, 12, 1) == -1);
  HUNDUN_CHECK(shared_row_oracle_field_coordinate(0, 6, 6, 12, 1) == 6);
  HUNDUN_CHECK(
      legacy_shared_row_oracle_field_coordinate(0, 6, 12, 12) == -6);
  bool oracle_rejected = false;
  try {
    static_cast<void>(
        shared_row_oracle_field_coordinate(1, 4, 4, 12, 1));
  } catch (const runtime::Error &) {
    oracle_rejected = true;
  }
  HUNDUN_CHECK(oracle_rejected);
  HUNDUN_CHECK(
      legacy_shared_row_oracle_field_coordinate(1, 4, 8, 12) == -3);

  Fixture fixture(mpi, FixtureKind::outside_tetrahedron, 1);
  HUNDUN_CHECK(fixture.ghost_plan->maximum_halo_reach() > 0U);
  auto halo = runtime::HaloExchange::create(
      fixture.decomposition,
      runtime::ExchangePlan::create(
          fixture.decomposition, fixture.decomposition.local_extent(),
          static_cast<int>(fixture.ghost_plan->maximum_halo_reach())));
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan, nullptr,
      &fixture.transform, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::bdf2, 0.01, 0.02);
  const double diagonal = stencil.alpha0 *
                          fixture.geometry.cell_volume_m3(0U) / stencil.dt_s;
  auto state = fixture.make_state(false, false, false, std::nullopt,
                                  WallCellHistory{1.0, 1.0});
  bool rejected = false;
  std::string message;
  try {
    static_cast<void>(
        flow::test::ImmersedFlowTestAccess::wall_predictor_values(
            immersed_flow, state, 1.0, stencil, diagonal));
  } catch (const runtime::Error &error) {
    rejected = true;
    message = error.what();
  }
  double all_rejected = rejected ? 1.0 : 0.0;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &all_rejected, 1, MPI_DOUBLE,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(all_rejected == 1.0);
  HUNDUN_CHECK(message ==
               "halo field layout is incompatible with the exchange plan");
}

void exercise_nonzero_inlet_exact_schur_linearity(
    const runtime::MpiContext &mpi) {
  Fixture fixture(mpi, FixtureKind::outside_open, 0, 1.0);
  auto halo = runtime::HaloExchange::create(
      fixture.decomposition,
      runtime::ExchangePlan::create(
          fixture.decomposition, fixture.decomposition.local_extent(),
          static_cast<int>(fixture.ghost_plan->maximum_halo_reach())));
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan, nullptr,
      &fixture.transform, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);
  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                          1.0,
                                          0.01,
                                          std::nullopt,
                                          std::nullopt,
                                          std::nullopt};
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);
  auto state = fixture.make_state();
  linear::SolveControl pressure_control;
  pressure_control.max_iterations = 1U;
  static_cast<void>(
      immersed_flow.attempt(state, physics, stencil, {}, pressure_control));

  const std::size_t count = fixture.domain->active_cells().owned_active_count();
  std::vector<flow::test::ImmersedFlowWallGradientSnapshot> zero_gradients;
  zero_gradients.reserve(fixture.domain->links().size());
  for (const auto &link : fixture.domain->links())
    zero_gradients.push_back({link.id, 0.0});
  std::sort(zero_gradients.begin(), zero_gradients.end(),
            [](const auto &left, const auto &right) {
              return left.link < right.link;
            });
  zero_gradients.erase(
      std::unique(zero_gradients.begin(), zero_gradients.end(),
                  [](const auto &left, const auto &right) {
                    return left.link == right.link;
                  }),
      zero_gradients.end());

  std::vector<double> zero(count, 0.0);
  std::vector<double> first(count);
  std::vector<double> second(count);
  std::vector<double> combined(count);
  constexpr double alpha = 0.375;
  constexpr double beta = -0.625;
  for (std::size_t row = 0U; row < count; ++row) {
    const auto id = fixture.domain->active_cells().ordered_global_ids()[row];
    first[row] =
        std::sin(0.017 * static_cast<double>((id % 104729U) + 1U));
    second[row] =
        std::cos(0.029 * static_cast<double>((id % 65537U) + 3U));
    combined[row] = alpha * first[row] + beta * second[row];
  }
  const auto zero_action =
      flow::test::ImmersedFlowTestAccess::active_pressure_operator_values(
          immersed_flow, zero, zero_gradients);
  const auto first_action =
      flow::test::ImmersedFlowTestAccess::active_pressure_operator_values(
          immersed_flow, first, zero_gradients);
  const auto second_action =
      flow::test::ImmersedFlowTestAccess::active_pressure_operator_values(
          immersed_flow, second, zero_gradients);
  const auto combined_action =
      flow::test::ImmersedFlowTestAccess::active_pressure_operator_values(
          immersed_flow, combined, zero_gradients);
  HUNDUN_CHECK(zero_action.size() == count);
  HUNDUN_CHECK(first_action.size() == count);
  HUNDUN_CHECK(second_action.size() == count);
  HUNDUN_CHECK(combined_action.size() == count);

  double norms[3]{};
  for (std::size_t row = 0U; row < count; ++row) {
    const double expected =
        alpha * first_action[row] + beta * second_action[row];
    const double defect = combined_action[row] - expected;
    norms[0] += zero_action[row] * zero_action[row];
    norms[1] += defect * defect;
    norms[2] += expected * expected;
  }
  mpi.allreduce_fp64_in_place(norms, 3U,
                              runtime::Fp64ReductionOperation::sum);
  const double zero_l2 = std::sqrt(norms[0]);
  const double defect_l2 = std::sqrt(norms[1]);
  const double scale_l2 = std::sqrt(norms[2]);
  const double tolerance = 1.0e-10 * std::max(1.0, scale_l2);
  if (mpi.rank() == 0)
    std::cerr << std::setprecision(17)
              << "nonzero_inlet_exact_schur_linearity zero_l2=" << zero_l2
              << " defect_l2=" << defect_l2
              << " scale_l2=" << scale_l2
              << " tolerance=" << tolerance << '\n';
  HUNDUN_CHECK(zero_l2 <= tolerance);
  HUNDUN_CHECK(defect_l2 <= tolerance);

  auto mutation = zero_action;
  HUNDUN_CHECK(!mutation.empty());
  mutation.front() += 1.0;
  double mutation_l2_squared = 0.0;
  for (const double value : mutation)
    mutation_l2_squared += value * value;
  mpi.allreduce_fp64_in_place(&mutation_l2_squared, 1U,
                              runtime::Fp64ReductionOperation::sum);
  HUNDUN_CHECK(std::sqrt(mutation_l2_squared) > tolerance);
}

void run() {
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  if (std::getenv("HUNDUN_EXACT_SCHUR_LINEARITY_ONLY") != nullptr) {
    exercise_nonzero_inlet_exact_schur_linearity(mpi);
    return;
  }
  if (std::getenv("HUNDUN_TASK11_SHARED_ROW_ONLY") != nullptr) {
    exercise_immersed_wall_face_history(mpi);
    exercise_underprovisioned_wall_halo_rejection(mpi);
    return;
  }
  if (std::getenv("HUNDUN_TASK11_CHECKERBOARD_ONLY") != nullptr) {
    static_cast<void>(exercise_constraint_mode(
        mpi, FixtureKind::outside_periodic, false, true,
        flow::MomentumTimeOrder::backward_euler, 0.001, 0.0));
    return;
  }
  exercise_pressure_outlet_face_history(mpi);
  exercise_immersed_wall_face_history(mpi);
  Fixture fixture(mpi);
  auto body_halo = runtime::HaloExchange::create(
      fixture.decomposition,
      runtime::ExchangePlan::create(fixture.decomposition,
                                    fixture.decomposition.local_extent(), 2));
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::ConjugateGradientSolver body_fitted_pressure_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  {
    linear::JacobiPreconditioner direct_mx(execution);
    linear::JacobiPreconditioner direct_my(execution);
    linear::JacobiPreconditioner direct_mz(execution);
    linear::JacobiPreconditioner direct_pressure_pc(execution);
    linear::JacobiPreconditioner wrapped_mx(execution);
    linear::JacobiPreconditioner wrapped_my(execution);
    linear::JacobiPreconditioner wrapped_mz(execution);
    linear::JacobiPreconditioner wrapped_pressure_pc(execution);
    auto direct = flow::FixedStepConstantDensityFlow::create(
        fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, mpi, execution, body_halo, momentum_solver,
        {&direct_mx, &direct_my, &direct_mz}, body_fitted_pressure_solver,
        direct_pressure_pc);
    auto wrapped = flow::FixedStepImmersedFlow::create(
        fixture.decomposition, fixture.topology, fixture.geometry,
        fixture.boundaries, nullptr, nullptr, nullptr, nullptr, nullptr, mpi,
        execution, body_halo, momentum_solver,
        {&wrapped_mx, &wrapped_my, &wrapped_mz},
        body_fitted_pressure_solver,
        wrapped_pressure_pc);
    auto direct_state = fixture.make_body_fitted_state();
    auto wrapped_state = fixture.make_body_fitted_state();
    const auto dispatch_stencil = flow::make_momentum_time_stencil(
        flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);
    const auto direct_report =
        direct.attempt(direct_state, 1.0, 0.01, dispatch_stencil, {}, {});
    const flow::ImmersedFlowPhysics dispatch_physics{config::DensityModel::constant,
                                               1.0,
                                               0.01,
                                               std::nullopt,
                                               std::nullopt,
                                               std::nullopt};
    const auto wrapped_result = wrapped.attempt(wrapped_state, dispatch_physics,
                                                dispatch_stencil, {}, {});
    const auto &wrapped_report =
        std::get<flow::StepAttemptReport>(wrapped_result.base);
    HUNDUN_CHECK(direct_report.disposition == wrapped_report.disposition);
    HUNDUN_CHECK(direct_report.reason == wrapped_report.reason);
    HUNDUN_CHECK(direct_report.pressure_corrector_count ==
                 wrapped_report.pressure_corrector_count);
    HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
        direct_state.snapshot(flow::FlowLayer::history),
        wrapped_state.snapshot(flow::FlowLayer::history)));
    HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
        direct_state.snapshot(flow::FlowLayer::committed),
        wrapped_state.snapshot(flow::FlowLayer::committed)));
    HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
        direct_state.snapshot(flow::FlowLayer::trial),
        wrapped_state.snapshot(flow::FlowLayer::trial)));
    HUNDUN_CHECK(test::accepted_step_metadata_bitwise_equal(
        direct_state.metadata(), wrapped_state.metadata()));
  }
  auto halo = runtime::HaloExchange::create(
      fixture.decomposition,
      runtime::ExchangePlan::create(
          fixture.decomposition, fixture.decomposition.local_extent(),
          static_cast<int>(fixture.ghost_plan->maximum_halo_reach())));
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan, nullptr,
      &fixture.transform, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);
  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                    1.0,
                                    0.01,
                                    std::nullopt,
                                    std::nullopt,
                                    std::nullopt};
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);
  const auto reference_stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, 0.1, 0.0);
  const auto momentum_algebra =
      flow::test::ImmersedFlowTestAccess::active_momentum_reference_probe(
          immersed_flow, physics.rho_ref_kg_per_m3, physics.dynamic_viscosity_pa_s,
          reference_stencil);
  HUNDUN_CHECK(momentum_algebra[0] <=
               512.0 * std::numeric_limits<double>::epsilon());
  HUNDUN_CHECK(momentum_algebra[1] > 0.0);
  HUNDUN_CHECK(momentum_algebra[2] > reference_stencil.alpha0 *
                                         physics.rho_ref_kg_per_m3 /
                                         reference_stencil.dt_s);
  HUNDUN_CHECK(momentum_algebra[3] > 0.0);
  HUNDUN_CHECK(momentum_algebra[6] > 1.0e-6);
  HUNDUN_CHECK(momentum_algebra[7] <=
               512.0 * std::numeric_limits<double>::epsilon());
  HUNDUN_CHECK(momentum_algebra[8] > 0.0);

  const auto zero_viscosity_algebra =
      flow::test::ImmersedFlowTestAccess::active_momentum_reference_probe(
          immersed_flow, physics.rho_ref_kg_per_m3, 0.0, reference_stencil);
  const double reference_time_rate = reference_stencil.alpha0 *
                                     physics.rho_ref_kg_per_m3 /
                                     reference_stencil.dt_s;
  HUNDUN_CHECK_NEAR(
      zero_viscosity_algebra[2], reference_time_rate,
      64.0 * std::numeric_limits<double>::epsilon() * reference_time_rate);
  HUNDUN_CHECK(zero_viscosity_algebra[3] == 0.0);
  HUNDUN_CHECK(zero_viscosity_algebra[6] == 0.0);
  HUNDUN_CHECK(zero_viscosity_algebra[7] == 0.0);
  HUNDUN_CHECK(zero_viscosity_algebra[8] == 0.0);

  const auto explicitly_stable_stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, 0.001, 0.0);
  const auto explicitly_stable_algebra =
      flow::test::ImmersedFlowTestAccess::active_momentum_reference_probe(
          immersed_flow, physics.rho_ref_kg_per_m3, physics.dynamic_viscosity_pa_s,
          explicitly_stable_stencil);
  HUNDUN_CHECK(explicitly_stable_algebra[0] <=
               512.0 * std::numeric_limits<double>::epsilon());
  HUNDUN_CHECK(explicitly_stable_algebra[1] > 0.0);
  const double explicitly_stable_time_rate = explicitly_stable_stencil.alpha0 *
                                             physics.rho_ref_kg_per_m3 /
                                             explicitly_stable_stencil.dt_s;
  HUNDUN_CHECK(explicitly_stable_algebra[2] > explicitly_stable_time_rate);
  HUNDUN_CHECK(explicitly_stable_algebra[3] > 0.0);
  HUNDUN_CHECK(explicitly_stable_algebra[6] > 1.0e-6);
  HUNDUN_CHECK(explicitly_stable_algebra[7] <=
               512.0 * std::numeric_limits<double>::epsilon());
  HUNDUN_CHECK(explicitly_stable_algebra[8] > 0.0);

  const auto explicitly_stable_bdf2 = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::bdf2, 0.001, 0.001);
  const auto bdf2_algebra =
      flow::test::ImmersedFlowTestAccess::active_momentum_reference_probe(
          immersed_flow, physics.rho_ref_kg_per_m3,
          physics.dynamic_viscosity_pa_s, explicitly_stable_bdf2);
  const double bdf2_time_rate = explicitly_stable_bdf2.alpha0 *
                                physics.rho_ref_kg_per_m3 /
                                explicitly_stable_bdf2.dt_s;
  HUNDUN_CHECK(bdf2_algebra[2] > bdf2_time_rate);
  HUNDUN_CHECK(bdf2_algebra[3] > 0.0);
  HUNDUN_CHECK(bdf2_algebra[4] <=
               512.0 * std::numeric_limits<double>::epsilon());
  HUNDUN_CHECK(bdf2_algebra[5] > 1.0e-6);
  HUNDUN_CHECK(bdf2_algebra[6] > 1.0e-6);
  HUNDUN_CHECK(bdf2_algebra[7] <=
               512.0 * std::numeric_limits<double>::epsilon());
  HUNDUN_CHECK(bdf2_algebra[8] > 0.0);

  auto state = fixture.make_state();
  auto accepted = immersed_flow.attempt(state, physics, stencil, {}, {});
  const auto &base = std::get<flow::StepAttemptReport>(accepted.base);
  HUNDUN_CHECK(base.disposition == flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(base.reason == flow::StepFailureReason::none);
  HUNDUN_CHECK(base.pressure_corrector_count == 2U);
  HUNDUN_CHECK(base.final_continuity_normalized_l2 <= 1.0e-10);
  for (const double value : base.final_momentum_normalized_l2)
    HUNDUN_CHECK(value <= 1.0e-9);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::last_corrector_count(immersed_flow) ==
               2U);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow) !=
               0U);
  HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::pressure_apply_schedule(
                   immersed_flow) == 1234U);
  const auto current_pressure_algebra =
      flow::test::ImmersedFlowTestAccess::active_pressure_current_algebra_probe(
          immersed_flow);
  HUNDUN_CHECK(std::isfinite(current_pressure_algebra[0]));
  HUNDUN_CHECK(std::isfinite(current_pressure_algebra[1]));
  HUNDUN_CHECK(std::isfinite(current_pressure_algebra[2]));
  HUNDUN_CHECK(current_pressure_algebra[0] <=
               4096.0 * std::numeric_limits<double>::epsilon());
  HUNDUN_CHECK(current_pressure_algebra[1] > 0.0);
  HUNDUN_CHECK(current_pressure_algebra[2] > 0.0);
  if (mpi.rank() == 0)
    std::cerr << std::setprecision(17)
              << "current_compact_pressure_algebra symmetry="
              << current_pressure_algebra[0]
              << " rayleigh=" << current_pressure_algebra[1] << ','
              << current_pressure_algebra[2] << '\n';
  const auto active_exchange =
      flow::test::ImmersedFlowTestAccess::active_exchange_performance(
          immersed_flow);
  const auto check_sparse_exchange = [&](const auto &counters) {
    HUNDUN_CHECK(counters.completed_exchanges > 0U);
    HUNDUN_CHECK(counters.receive_payload_bytes ==
                 counters.completed_exchanges *
                     active_exchange.ghost_active_count * sizeof(double));
  };
  check_sparse_exchange(active_exchange.pressure);
  for (const auto &counters : active_exchange.momentum)
    check_sparse_exchange(counters);
  HUNDUN_CHECK(!accepted.linear_solve_failure.has_value());
  auto nested_failure_flow = flow::FixedStepImmersedFlow::create(
      fixture.decomposition, fixture.topology, fixture.geometry,
      fixture.boundaries, &*fixture.domain, &*fixture.ghost_plan, nullptr,
      &fixture.transform, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);
  const auto exercise_nested_failure = [&](auto phase,
                                           std::uint32_t component) {
    auto failed_state = fixture.make_state();
    const auto history_before =
        failed_state.snapshot(flow::FlowLayer::history);
    const auto committed_before =
        failed_state.snapshot(flow::FlowLayer::committed);
    const auto metadata_before = failed_state.metadata();
    flow::test::ImmersedFlowTestAccess::set_linear_solve_failure(
        phase, component);
    const auto failed =
        nested_failure_flow.attempt(failed_state, physics, stencil, {}, {});
    flow::test::ImmersedFlowTestAccess::clear_linear_solve_failure();
    const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
    HUNDUN_CHECK(failed_base.disposition ==
                 flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(failed_base.reason ==
                 flow::StepFailureReason::pressure_linear_solve);
    HUNDUN_CHECK(failed.linear_solve_failure.has_value());
    HUNDUN_CHECK(failed.linear_solve_failure->phase == phase);
    HUNDUN_CHECK(failed.linear_solve_failure->pressure_corrector_index == 1U);
    HUNDUN_CHECK(failed.linear_solve_failure->component_index == component);
    if (phase ==
        flow::ImmersedLinearSolvePhase::pressure_independent_residual) {
      HUNDUN_CHECK(
          failed.linear_solve_failure->solve.reason ==
              linear::SolveTerminationReason::converged ||
          failed.linear_solve_failure->solve.reason ==
              linear::SolveTerminationReason::zero_right_hand_side);
    } else {
      HUNDUN_CHECK(failed.linear_solve_failure->solve.reason ==
                   linear::SolveTerminationReason::maximum_iterations);
    }
    if (phase !=
        flow::ImmersedLinearSolvePhase::pressure_independent_residual)
      HUNDUN_CHECK(failed.linear_solve_failure->solve.iterations > 0U);
    HUNDUN_CHECK(std::isfinite(
        failed.linear_solve_failure->solve.final_residual));
    HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
        failed_state.snapshot(flow::FlowLayer::history), history_before));
    HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
        failed_state.snapshot(flow::FlowLayer::committed), committed_before));
    HUNDUN_CHECK(test::accepted_step_metadata_bitwise_equal(
        failed_state.metadata(), metadata_before));
    return *failed.linear_solve_failure;
  };
  static_cast<void>(exercise_nested_failure(
      flow::ImmersedLinearSolvePhase::pressure_affine_momentum, 1U));
  const auto outer_failure = exercise_nested_failure(
      flow::ImmersedLinearSolvePhase::pressure_outer,
      std::numeric_limits<std::uint32_t>::max());
  const auto residual_failure = exercise_nested_failure(
      flow::ImmersedLinearSolvePhase::pressure_independent_residual,
      std::numeric_limits<std::uint32_t>::max());
  HUNDUN_CHECK(outer_failure.phase != residual_failure.phase);
  HUNDUN_CHECK(std::isfinite(residual_failure.independent_residual_l2));
  HUNDUN_CHECK(std::isfinite(residual_failure.acceptance_threshold));
  HUNDUN_CHECK(residual_failure.independent_residual_l2 >
               residual_failure.acceptance_threshold);
  HUNDUN_CHECK(
      !flow::test::ImmersedFlowTestAccess::has_active_pressure_reference(immersed_flow));
  const auto exact_diagonal_contract =
      flow::test::ImmersedFlowTestAccess::exact_predictor_diagonal_contract(
          immersed_flow);
  if (exact_diagonal_contract.advertised &&
      exact_diagonal_contract.relative_difference >
          512.0 * std::numeric_limits<double>::epsilon() &&
      mpi.rank() == 0)
    std::cerr << std::setprecision(17)
              << "exact_schur_diagonal_contract declared="
              << exact_diagonal_contract.declared_value
              << " applied_basis="
              << exact_diagonal_contract.applied_basis_value
              << " relative_difference="
              << exact_diagonal_contract.relative_difference << '\n';
  HUNDUN_CHECK(!exact_diagonal_contract.advertised ||
               exact_diagonal_contract.relative_difference <=
                   512.0 * std::numeric_limits<double>::epsilon());
  const auto exact_workspace_after_attempt =
      flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
  HUNDUN_CHECK(
      exact_workspace_after_attempt.cached_response_consumption_count == 2U);
  HUNDUN_CHECK(
      exact_workspace_after_attempt.cached_input_mutation_rejected);
  HUNDUN_CHECK(exact_workspace_after_attempt.accepted_inner_response_count >
               0U);
  HUNDUN_CHECK(exact_workspace_after_attempt.accepted_inner_solve_count ==
               3U *
                   exact_workspace_after_attempt.accepted_inner_response_count);
  HUNDUN_CHECK(exact_workspace_after_attempt.accepted_inner_matvec_count >=
               exact_workspace_after_attempt.accepted_inner_solve_count);
  HUNDUN_CHECK(exact_workspace_after_attempt.accepted_inner_iteration_count <=
               exact_workspace_after_attempt.accepted_inner_matvec_count);
  HUNDUN_CHECK(
      exact_workspace_after_attempt
          .accepted_inner_preconditioner_apply_count <=
      exact_workspace_after_attempt.accepted_inner_matvec_count);
  HUNDUN_CHECK(
      exact_workspace_after_attempt.accepted_inner_global_reduction_count >
      exact_workspace_after_attempt.accepted_inner_solve_count);
  HUNDUN_CHECK(
      exact_workspace_after_attempt.last_response_allocation_event_count ==
      0U);
  {
    const auto exact_copy = exact_workspace_after_attempt;
    HUNDUN_CHECK(
        accepted_inner_work_equal(exact_workspace_after_attempt, exact_copy));
    auto ordinary_mutation = exact_copy;
    ++ordinary_mutation.accepted_inner_matvec_count;
    HUNDUN_CHECK(!accepted_inner_work_equal(exact_workspace_after_attempt,
                                            ordinary_mutation));
  }
  if (std::getenv("HUNDUN_TASK11_EXACT_CACHE_ONLY") != nullptr)
    return;
  const auto local_revision =
      flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow);
  std::uint64_t revision_min{};
  std::uint64_t revision_max{};
  HUNDUN_CHECK(MPI_Allreduce(&local_revision, &revision_min, 1, MPI_UINT64_T,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&local_revision, &revision_max, 1, MPI_UINT64_T,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(revision_min == revision_max);
  {
    const std::size_t count =
        fixture.domain->active_cells().owned_active_count();
    std::vector<double> pressure_probe(count);
    for (std::size_t row = 0U; row < count; ++row) {
      const auto id = fixture.domain->active_cells().ordered_global_ids()[row];
      pressure_probe[row] =
          std::sin(0.031 * static_cast<double>((id % 104729U) + 1U));
    }
    std::vector<flow::test::ImmersedFlowWallGradientSnapshot> zero_wall_gradients;
    zero_wall_gradients.reserve(fixture.domain->links().size());
    for (const auto &link : fixture.domain->links())
      zero_wall_gradients.push_back({link.id, 0.0});
    std::sort(zero_wall_gradients.begin(), zero_wall_gradients.end(),
              [](const auto &left, const auto &right) {
                return left.link < right.link;
              });
    zero_wall_gradients.erase(
        std::unique(zero_wall_gradients.begin(), zero_wall_gradients.end(),
                    [](const auto &left, const auto &right) {
                      return left.link == right.link;
                    }),
        zero_wall_gradients.end());
    const auto baseline_operator =
        flow::test::ImmersedFlowTestAccess::active_pressure_operator_values(
            immersed_flow, pressure_probe, zero_wall_gradients);
    auto changed_state = fixture.make_state();
    flow::test::ImmersedFlowTestAccess::set_momentum_time_diagonal_scale(1.5);
    const auto changed =
        immersed_flow.attempt(changed_state, physics, stencil, {}, {});
    flow::test::ImmersedFlowTestAccess::clear_momentum_time_diagonal_scale();
    const auto &changed_base = std::get<flow::StepAttemptReport>(changed.base);
    HUNDUN_CHECK(changed_base.disposition ==
                 flow::StepAttemptDisposition::committed);
    const auto changed_operator =
        flow::test::ImmersedFlowTestAccess::active_pressure_operator_values(
            immersed_flow, pressure_probe, zero_wall_gradients);
    HUNDUN_CHECK(changed_operator.size() == baseline_operator.size());
    double maximum_change = 0.0;
    for (std::size_t row = 0U; row < baseline_operator.size(); ++row)
      maximum_change =
          std::max(maximum_change,
                   std::abs(changed_operator[row] - baseline_operator[row]));
    mpi.allreduce_fp64_in_place(&maximum_change, 1U,
                                runtime::Fp64ReductionOperation::maximum);
    HUNDUN_CHECK(maximum_change >
                 256.0 * std::numeric_limits<double>::epsilon());
  }
  {
    auto changed_state = fixture.make_state();
    const auto changed_stencil = flow::make_momentum_time_stencil(
        flow::MomentumTimeOrder::backward_euler, 0.005, 0.0);
    const auto changed =
        immersed_flow.attempt(changed_state, physics, changed_stencil, {}, {});
    const auto &changed_base = std::get<flow::StepAttemptReport>(changed.base);
    HUNDUN_CHECK(changed_base.disposition ==
                 flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::pressure_revision(immersed_flow) >
                 local_revision);
  }
  const auto accepted_values = state.snapshot(flow::FlowLayer::committed);
  for (mesh::LocalFaceId face = 0U; face < fixture.topology.local_face_count();
       ++face) {
    const auto neighbour = fixture.topology.neighbour(face);
    if (!neighbour.has_value())
      continue;
    const bool owner_active =
        fixture.domain->region(fixture.topology.owner(face)) ==
        immersed::CellRegion::fluid;
    const bool neighbour_active =
        fixture.domain->region(*neighbour) == immersed::CellRegion::fluid;
    if (owner_active != neighbour_active) {
      HUNDUN_CHECK(accepted_values.face_mass_flux[face] == 0.0);
      HUNDUN_CHECK(!std::signbit(accepted_values.face_mass_flux[face]));
    }
  }

  {
    auto inter_step_state = fixture.make_state(true);
    linear::SolveControl inter_step_pressure_control;
    inter_step_pressure_control.atol = 1.0e-14;
    inter_step_pressure_control.rtol = 1.0e-13;
    linear::SolveControl inter_step_momentum_control;
    inter_step_momentum_control.atol = 1.0e-18;
    inter_step_momentum_control.rtol = 1.0e-13;
    const auto first = immersed_flow.attempt(inter_step_state, physics,
                                     explicitly_stable_stencil,
                                     inter_step_momentum_control,
                                     inter_step_pressure_control);
    const auto &first_base = std::get<flow::StepAttemptReport>(first.base);
    if (first_base.disposition != flow::StepAttemptDisposition::committed &&
        mpi.rank() == 0) {
      std::cerr << std::setprecision(17)
                << "a22_inter_step_first_rejected disposition="
                << static_cast<int>(first_base.disposition)
                << " reason=" << static_cast<int>(first_base.reason)
                << " correctors=" << first_base.pressure_corrector_count
                << " continuity="
                << first_base.final_continuity_normalized_l2
                << " pressure=" << first_base.final_pressure_residual_l2
                << " momentum="
                << first_base.final_momentum_normalized_l2[0] << ','
                << first_base.final_momentum_normalized_l2[1] << ','
                << first_base.final_momentum_normalized_l2[2]
                << " pressure0_reason="
                << static_cast<int>(first_base.pressure[0].reason)
                << " pressure0_iterations="
                << first_base.pressure[0].iterations
                << " pressure0_initial="
                << first_base.pressure[0].initial_residual
                << " pressure0_recursive="
                << first_base.pressure[0].recursive_residual
                << " pressure0_final="
                << first_base.pressure[0].final_residual << '\n';
    }
    if (first_base.disposition != flow::StepAttemptDisposition::committed) {
      const auto correction_report =
          flow::test::ImmersedFlowTestAccess::cell_pressure_correction_authority(
              immersed_flow);
      for (std::size_t corrector = 0U;
           corrector < correction_report.correctors.size(); ++corrector)
        print_cell_pressure_correction_record(
            mpi, fixture.topology.global_extent(), corrector + 1U,
            correction_report.correctors[corrector]);
    }
    HUNDUN_CHECK(first_base.disposition ==
                 flow::StepAttemptDisposition::committed);
    const auto second =
        immersed_flow.attempt(inter_step_state, physics, explicitly_stable_stencil,
                       inter_step_momentum_control,
                       inter_step_pressure_control);
    const auto &second_base = std::get<flow::StepAttemptReport>(second.base);
    const double inter_step_difference =
        flow::test::ImmersedFlowTestAccess::
            inter_step_pressure_authority_difference_l2(immersed_flow);
    HUNDUN_CHECK(second_base.disposition ==
                 flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(test::fp64_bits(inter_step_difference) ==
                 test::fp64_bits(0.0));

    const auto bdf2_stencil = flow::make_momentum_time_stencil(
        flow::MomentumTimeOrder::bdf2, 0.001, 0.001);
    const auto third = immersed_flow.attempt(
        inter_step_state, physics, bdf2_stencil, inter_step_momentum_control,
        inter_step_pressure_control);
    const auto &third_base = std::get<flow::StepAttemptReport>(third.base);
    HUNDUN_CHECK(third_base.disposition ==
                 flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(
        test::fp64_bits(flow::test::ImmersedFlowTestAccess::
                            bdf2_history_pressure_authority_difference_l2(
                                immersed_flow)) == test::fp64_bits(0.0));

    const auto authority_before_failure =
        flow::test::ImmersedFlowTestAccess::pressure_authority_fingerprints(
            immersed_flow);
    const auto history_before_failure =
        inter_step_state.snapshot(flow::FlowLayer::history);
    const auto committed_before_failure =
        inter_step_state.snapshot(flow::FlowLayer::committed);
    const auto trial_before_failure =
        inter_step_state.snapshot(flow::FlowLayer::trial);
    const auto metadata_before_failure = inter_step_state.metadata();
    const auto inner_work_before_failure =
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
    flow::test::ImmersedFlowTestAccess::set_failure_stage(
        flow::test::ImmersedFlowAttemptFailureStage::after_corrector_1,
        mpi.size() > 1 ? 1 : 0);
    const auto failed = immersed_flow.attempt(
        inter_step_state, physics, bdf2_stencil, inter_step_momentum_control,
        inter_step_pressure_control);
    flow::test::ImmersedFlowTestAccess::clear_failure_stage();
    const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
    HUNDUN_CHECK(failed_base.disposition ==
                 flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(failed_base.reason ==
                 flow::StepFailureReason::non_finite_trial);
    check_state_equal(inter_step_state, history_before_failure,
                      committed_before_failure, trial_before_failure,
                      metadata_before_failure);
    HUNDUN_CHECK(
        flow::test::ImmersedFlowTestAccess::pressure_authority_fingerprints(
            immersed_flow) == authority_before_failure);
    HUNDUN_CHECK(accepted_inner_work_equal(
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow),
        inner_work_before_failure));

    const auto retry = immersed_flow.attempt(
        inter_step_state, physics, bdf2_stencil, inter_step_momentum_control,
        inter_step_pressure_control);
    const auto &retry_base = std::get<flow::StepAttemptReport>(retry.base);
    HUNDUN_CHECK(retry_base.disposition ==
                 flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(
        test::fp64_bits(flow::test::ImmersedFlowTestAccess::
                            inter_step_pressure_authority_difference_l2(
                                immersed_flow)) == test::fp64_bits(0.0));
    HUNDUN_CHECK(
        test::fp64_bits(flow::test::ImmersedFlowTestAccess::
                            bdf2_history_pressure_authority_difference_l2(
                                immersed_flow)) == test::fp64_bits(0.0));
  }

  for (const auto &[failure_stage, expected_reason] :
       {std::pair{flow::test::ImmersedFlowAttemptFailureStage::after_corrector_1,
                  flow::StepFailureReason::non_finite_trial},
        std::pair{flow::test::ImmersedFlowAttemptFailureStage::after_corrector_2,
                  flow::StepFailureReason::non_finite_trial},
        std::pair{
            flow::test::ImmersedFlowAttemptFailureStage::final_momentum_residual,
            flow::StepFailureReason::final_momentum_residual},
        std::pair{flow::test::ImmersedFlowAttemptFailureStage::after_final_transport,
                  flow::StepFailureReason::transport_failure},
        std::pair{flow::test::ImmersedFlowAttemptFailureStage::final_wall_penetration,
                  flow::StepFailureReason::final_continuity_residual},
        std::pair{
            flow::test::ImmersedFlowAttemptFailureStage::final_continuity_residual,
            flow::StepFailureReason::final_continuity_residual},
        std::pair{
            flow::test::ImmersedFlowAttemptFailureStage::final_pressure_residual,
            flow::StepFailureReason::final_pressure_residual}}) {
    auto rollback = fixture.make_state();
    const auto history = rollback.snapshot(flow::FlowLayer::history);
    const auto committed = rollback.snapshot(flow::FlowLayer::committed);
    const auto trial = rollback.snapshot(flow::FlowLayer::trial);
    const auto metadata = rollback.metadata();
    const auto inner_work_before_failure =
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
    flow::test::ImmersedFlowTestAccess::set_failure_stage(failure_stage);
    const auto failed = immersed_flow.attempt(rollback, physics, stencil, {}, {});
    flow::test::ImmersedFlowTestAccess::clear_failure_stage();
    const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
    HUNDUN_CHECK(failed_base.disposition ==
                 flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(failed_base.reason == expected_reason);
    HUNDUN_CHECK(failed_base.lowest_failing_rank == (mpi.size() > 1 ? 1 : 0));
    HUNDUN_CHECK(
        failed_base.pressure_corrector_count ==
        (failure_stage ==
                 flow::test::ImmersedFlowAttemptFailureStage::after_corrector_1
             ? 1U
             : 2U));
    check_state_equal(rollback, history, committed, trial, metadata);
    HUNDUN_CHECK(accepted_inner_work_equal(
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow),
        inner_work_before_failure));
  }

  const int selected_rank = mpi.size() > 1 ? 1 : 0;
  {
    auto rollback = fixture.make_state();
    const auto history = rollback.snapshot(flow::FlowLayer::history);
    const auto committed = rollback.snapshot(flow::FlowLayer::committed);
    const auto trial = rollback.snapshot(flow::FlowLayer::trial);
    const auto metadata = rollback.metadata();
    const auto inner_work_before_failure =
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
    auto rank_local_physics = physics;
    if (mpi.rank() == selected_rank)
      rank_local_physics.cp_J_per_kg_K = 1.0;
    const auto failed =
        immersed_flow.attempt(rollback, rank_local_physics, stencil, {}, {});
    const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
    HUNDUN_CHECK(failed_base.disposition ==
                 flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(failed_base.reason == flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(failed_base.lowest_failing_rank == selected_rank);
    check_state_equal(rollback, history, committed, trial, metadata);
    HUNDUN_CHECK(accepted_inner_work_equal(
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow),
        inner_work_before_failure));
  }
  {
    auto rollback = fixture.make_state();
    const auto history = rollback.snapshot(flow::FlowLayer::history);
    const auto committed = rollback.snapshot(flow::FlowLayer::committed);
    const auto trial = rollback.snapshot(flow::FlowLayer::trial);
    const auto metadata = rollback.metadata();
    const auto inner_work_before_failure =
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
    flow::test::ImmersedFlowTestAccess::set_momentum_time_diagonal_scale(
        mpi.rank() == selected_rank ? std::numeric_limits<double>::quiet_NaN()
                                    : 1.0);
    const auto failed = immersed_flow.attempt(rollback, physics, stencil, {}, {});
    flow::test::ImmersedFlowTestAccess::clear_momentum_time_diagonal_scale();
    const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
    HUNDUN_CHECK(failed_base.disposition ==
                 flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(failed_base.reason == flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(failed_base.lowest_failing_rank == selected_rank);
    check_state_equal(rollback, history, committed, trial, metadata);
    HUNDUN_CHECK(accepted_inner_work_equal(
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow),
        inner_work_before_failure));
  }
  {
    auto rollback = fixture.make_state();
    const auto history = rollback.snapshot(flow::FlowLayer::history);
    const auto committed = rollback.snapshot(flow::FlowLayer::committed);
    const auto trial = rollback.snapshot(flow::FlowLayer::trial);
    const auto metadata = rollback.metadata();
    const auto inner_work_before_failure =
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
    auto rank_local_stencil = stencil;
    if (mpi.rank() == selected_rank)
      rank_local_stencil.alpha0 += 1.0;
    const auto failed =
        immersed_flow.attempt(rollback, physics, rank_local_stencil, {}, {});
    const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
    HUNDUN_CHECK(failed_base.disposition ==
                 flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(failed_base.reason == flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(failed_base.lowest_failing_rank == selected_rank);
    check_state_equal(rollback, history, committed, trial, metadata);
    HUNDUN_CHECK(accepted_inner_work_equal(
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow),
        inner_work_before_failure));
  }
  for (const auto &[failure, expected_reason] :
       {std::pair{flow::test::ImmersedFlowWallInputFailure::non_positive_density,
                  flow::StepFailureReason::non_finite_trial},
        std::pair{flow::test::ImmersedFlowWallInputFailure::non_finite_density,
                  flow::StepFailureReason::non_finite_trial},
        std::pair{flow::test::ImmersedFlowWallInputFailure::non_positive_coefficient,
                  flow::StepFailureReason::non_finite_trial},
        std::pair{flow::test::ImmersedFlowWallInputFailure::non_finite_coefficient,
                  flow::StepFailureReason::non_finite_trial}}) {
    auto rollback = fixture.make_state();
    const auto history = rollback.snapshot(flow::FlowLayer::history);
    const auto committed = rollback.snapshot(flow::FlowLayer::committed);
    const auto trial = rollback.snapshot(flow::FlowLayer::trial);
    const auto metadata = rollback.metadata();
    const auto inner_work_before_failure =
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
    flow::test::ImmersedFlowTestAccess::set_wall_input_failure(failure,
                                                             selected_rank);
    const auto failed = immersed_flow.attempt(rollback, physics, stencil, {}, {});
    flow::test::ImmersedFlowTestAccess::clear_wall_input_failure();
    const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
    HUNDUN_CHECK(failed_base.disposition ==
                 flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(failed_base.reason == expected_reason);
    HUNDUN_CHECK(failed_base.lowest_failing_rank == selected_rank);
    check_state_equal(rollback, history, committed, trial, metadata);
    HUNDUN_CHECK(accepted_inner_work_equal(
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow),
        inner_work_before_failure));
  }

  {
    auto rollback = fixture.make_state();
    const auto history = rollback.snapshot(flow::FlowLayer::history);
    const auto committed = rollback.snapshot(flow::FlowLayer::committed);
    const auto trial = rollback.snapshot(flow::FlowLayer::trial);
    const auto metadata = rollback.metadata();
    const auto inner_work_before_failure =
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
    flow::test::ImmersedFlowTestAccess::set_wall_input_failure(
        flow::test::ImmersedFlowWallInputFailure::stale_preconditioner_revision,
        selected_rank);
    const auto failed = immersed_flow.attempt(rollback, physics, stencil, {}, {});
    flow::test::ImmersedFlowTestAccess::clear_wall_input_failure();
    const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
    HUNDUN_CHECK(failed_base.disposition ==
                 flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(failed_base.reason == flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(failed_base.lowest_failing_rank == selected_rank);
    check_state_equal(rollback, history, committed, trial, metadata);
    HUNDUN_CHECK(accepted_inner_work_equal(
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow),
        inner_work_before_failure));
  }

  for (const auto &corruption :
       {std::pair{true, false}, std::pair{false, true}}) {
    const bool selected = mpi.rank() == 0;
    auto rollback = fixture.make_state(false, selected && corruption.first,
                                       selected && corruption.second);
    const auto history = rollback.snapshot(flow::FlowLayer::history);
    const auto committed = rollback.snapshot(flow::FlowLayer::committed);
    const auto trial = rollback.snapshot(flow::FlowLayer::trial);
    const auto metadata = rollback.metadata();
    const auto inner_work_before_failure =
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow);
    const auto failed = immersed_flow.attempt(rollback, physics, stencil, {}, {});
    const auto &failed_base = std::get<flow::StepAttemptReport>(failed.base);
    HUNDUN_CHECK(failed_base.disposition ==
                 flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(failed_base.reason == flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(failed_base.lowest_failing_rank == 0);
    check_state_equal(rollback, history, committed, trial, metadata);
    HUNDUN_CHECK(accepted_inner_work_equal(
        flow::test::ImmersedFlowTestAccess::exact_predictor_workspace(immersed_flow),
        inner_work_before_failure));
  }

  const auto bdf2_residual = exercise_constraint_mode(
      mpi, FixtureKind::outside_periodic, false, false,
      flow::MomentumTimeOrder::bdf2, 0.01, 0.02);
  static_cast<void>(
      exercise_constraint_mode(mpi, FixtureKind::outside_open, true));
  const auto backward_euler_residual = exercise_constraint_mode(
      mpi, FixtureKind::outside_periodic, false, true,
      flow::MomentumTimeOrder::backward_euler, 0.001, 0.0);

  const auto exact = backward_euler_residual.used_residual;
  HUNDUN_CHECK(final_residual_vectors_agree(exact, exact));
  auto ordinary_mutation = exact;
  HUNDUN_CHECK(!ordinary_mutation.empty());
  ordinary_mutation.front() += 1.0;
  HUNDUN_CHECK(!final_residual_vectors_agree(exact, ordinary_mutation));
  auto size_mismatch = exact;
  size_mismatch.pop_back();
  HUNDUN_CHECK(!final_residual_vectors_agree(exact, size_mismatch));
  auto non_finite = exact;
  non_finite.front() = std::numeric_limits<double>::quiet_NaN();
  HUNDUN_CHECK(!final_residual_vectors_agree(exact, non_finite));

  HUNDUN_CHECK(routes_agree(backward_euler_residual));
  HUNDUN_CHECK(routes_agree(bdf2_residual));
}

} // namespace

int main(int argc, char **argv) {
  runtime::MpiEnvironment environment(argc, argv);
  return test::run(run);
}
