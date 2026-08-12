// SPDX-License-Identifier: Apache-2.0

#include "hundun/ib_ghost_stencil_plan.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "src/ib_deterministic_qr_detail.hpp"
#include "src/ib_ghost_stencil_plan_detail.hpp"
#include "src/ib_quadratic_reconstruction_detail.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <map>
#include <string>
#include <system_error>
#include <vector>

namespace {

using namespace hundun;
constexpr runtime::Int3 kExtent{12, 12, 12};

std::uint64_t bits(double value) {
  std::uint64_t result = 0U;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

double polynomial(runtime::Real3 point) {
  return 1.0 + 2.0 * point.x - 3.0 * point.y + 0.5 * point.z +
         0.7 * point.x * point.x - 0.2 * point.x * point.y +
         0.4 * point.x * point.z + 0.3 * point.y * point.y -
         0.6 * point.y * point.z + 0.9 * point.z * point.z;
}

runtime::Real3 polynomial_gradient(runtime::Real3 point) {
  return {2.0 + 1.4 * point.x - 0.2 * point.y + 0.4 * point.z,
          -3.0 - 0.2 * point.x + 0.6 * point.y - 0.6 * point.z,
          0.5 + 0.4 * point.x - 0.6 * point.y + 1.8 * point.z};
}

double cell_average(runtime::Int3 cell) {
  constexpr double h = 1.0 / static_cast<double>(kExtent.x);
  const runtime::Real3 center{(static_cast<double>(cell.x) + 0.5) * h,
                              (static_cast<double>(cell.y) + 0.5) * h,
                              (static_cast<double>(cell.z) + 0.5) * h};
  return polynomial(center) + (0.7 + 0.3 + 0.9) * h * h / 12.0;
}

double wall_normal_gradient_amplification(
    const immersed::QuadraticReconstruction &reconstruction,
    runtime::Real3 point, runtime::Real3 normal) {
  double result = 0.0;
  for (const auto &donor : immersed::detail::QuadraticReconstructionWeights::
           origin_constrained_directional_gradient_weights(reconstruction,
                                                           point, normal)) {
    result += std::abs(donor.weight);
  }
  return result;
}

bool donor_authority_equal(const std::vector<mesh::GlobalCellId> &lhs,
                           const std::vector<mesh::GlobalCellId> &rhs) {
  if (lhs.size() != rhs.size())
    return false;
  for (std::size_t index = 0U; index < lhs.size(); ++index)
    if (lhs[index] != rhs[index])
      return false;
  return true;
}

void prove_donor_authority_oracle_is_mutation_sensitive(
    const std::vector<mesh::GlobalCellId> &authority) {
  HUNDUN_CHECK(!authority.empty());
  HUNDUN_CHECK(donor_authority_equal(authority, authority));
  auto ordinary_mutation = authority;
  ordinary_mutation.front() ^= UINT64_C(1);
  HUNDUN_CHECK(!donor_authority_equal(authority, ordinary_mutation));
  auto nested_size_mutation = authority;
  nested_size_mutation.pop_back();
  HUNDUN_CHECK(!donor_authority_equal(authority, nested_size_mutation));
}

void check_feasible_owner_selector(const runtime::MpiContext &mpi) {
  constexpr int reach = 4;
  const std::vector<runtime::Box3> divergent_boxes{
      {{0, 0, 0}, {2, 2, 2}},
      {{8, 0, 0}, {10, 2, 2}},
      {{0, 0, 0}, {10, 2, 2}},
  };
  const std::vector<runtime::Int3> point_donors{{1, 0, 0}};
  const std::vector<runtime::Int3> pressure_donors{{9, 0, 0}};

  // Ranks 0 and 1 are the test-annotated association owners.  The helper has
  // no association-owner input, and only rank 2 covers the complete union.
  HUNDUN_CHECK(immersed::detail::select_wall_quadrature_execution_owner(
                   divergent_boxes, point_donors, pressure_donors, reach) ==
               2);
  HUNDUN_CHECK(immersed::detail::select_wall_quadrature_execution_owner(
                   divergent_boxes, point_donors, {}, reach) == 0);
  HUNDUN_CHECK(immersed::detail::select_wall_quadrature_execution_owner(
                   divergent_boxes, {}, pressure_donors, reach) == 1);

  const std::vector<runtime::Box3> tied_boxes{
      {{0, 0, 0}, {4, 4, 4}},
      {{0, 0, 0}, {4, 4, 4}},
  };
  HUNDUN_CHECK(immersed::detail::select_wall_quadrature_execution_owner(
                   tied_boxes, {{1, 1, 1}}, {{2, 2, 2}}, reach) == 0);

  bool local_ok = true;
  std::string local_message;
  try {
    static_cast<void>(
        immersed::detail::select_wall_quadrature_execution_owner(
            tied_boxes, {{100, 1, 1}}, {}, reach));
  } catch (const runtime::Error &error) {
    local_ok = false;
    local_message = error.what();
  }
  const auto status =
      runtime::collective_status(mpi, local_ok, local_message);
  HUNDUN_CHECK(!status.ok);
  HUNDUN_CHECK(status.message ==
               "wall quadrature donor exceeds owner Halo reach");
  HUNDUN_CHECK(status.failing_rank == 0);
}

std::vector<int>
triangle_owner_map(const immersed::WallQuadraturePlan &plan,
                   std::size_t triangle_count,
                   const runtime::MpiContext &mpi) {
  std::vector<int> local(triangle_count, mpi.size());
  for (const auto &point : plan.local_points()) {
    HUNDUN_CHECK(point.triangle < triangle_count);
    HUNDUN_CHECK(point.owner_rank == mpi.rank());
    auto &owner = local[static_cast<std::size_t>(point.triangle)];
    HUNDUN_CHECK(owner == mpi.size() || owner == point.owner_rank);
    owner = point.owner_rank;
  }
  std::vector<int> global(local.size(), mpi.size());
  HUNDUN_CHECK(MPI_Allreduce(local.data(), global.data(),
                             static_cast<int>(global.size()), MPI_INT, MPI_MIN,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(std::all_of(global.begin(), global.end(), [&](int owner) {
    return owner >= 0 && owner < mpi.size();
  }));
  return global;
}

runtime::Real3 midpoint(runtime::Real3 first, runtime::Real3 second) {
  return {(first.x + second.x) * 0.5, (first.y + second.y) * 0.5,
          (first.z + second.z) * 0.5};
}

std::vector<test::StlFixtureTriangle> refined_cube() {
  const auto coarse =
      test::translated(test::outward_cube(), {0.75, 0.75, 0.75});
  std::vector<test::StlFixtureTriangle> refined;
  refined.reserve(4U * coarse.size());
  for (const auto &triangle : coarse) {
    const auto ab = midpoint(triangle.vertices[0], triangle.vertices[1]);
    const auto bc = midpoint(triangle.vertices[1], triangle.vertices[2]);
    const auto ca = midpoint(triangle.vertices[2], triangle.vertices[0]);
    refined.push_back({triangle.file_normal, {triangle.vertices[0], ab, ca}});
    refined.push_back({triangle.file_normal, {ab, triangle.vertices[1], bc}});
    refined.push_back({triangle.file_normal, {ca, bc, triangle.vertices[2]}});
    refined.push_back({triangle.file_normal, {ab, bc, ca}});
  }
  return refined;
}

runtime::Int3 process_grid(int ranks, const std::string &name) {
  if (ranks == 1) {
    return {1, 1, 1};
  }
  if (ranks == 2) {
    return {2, 1, 1};
  }
  return name == "4x1x1" ? runtime::Int3{4, 1, 1} : runtime::Int3{2, 2, 1};
}

config::FlowCaseConfig config_for(int ranks, runtime::Int3 grid,
                                  bool warped = false) {
  config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "task5-wall";
  config.simulation_type = config::SimulationType::variable_density_flow;
  config.density_model = config::DensityModel::constant;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = grid;
  config.mesh.cells = kExtent;
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.mesh.mapping = warped ? config::MeshMapping::analytic_warped_box
                               : config::MeshMapping::uniform_box;
  if (warped) {
    config.mesh.warp_amplitude = runtime::Real3{0.02, -0.015, 0.01};
  }
  config.time.mode = config::TimeMode::fixed;
  config.time.steps = 1;
  config.time.initial_dt_s = 0.01;
  config.time.min_dt_s = 0.01;
  config.time.max_dt_s = 0.01;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.dynamic_viscosity_pa_s = 1.0e-3;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t patch = 0U; patch < names.size(); ++patch) {
    config::FlowBoundaryConfig boundary{};
    boundary.patch = names[patch];
    boundary.type = config::BoundaryType::no_slip_wall;
    config.boundaries[patch] = boundary;
  }
  return config;
}

class FixtureFile final {
public:
  explicit FixtureFile(const runtime::MpiContext &mpi) : mpi_(&mpi) {
    std::string text;
    if (mpi.rank() == 0) {
      path_ =
          std::filesystem::temp_directory_path() /
          ("hundun-task5-wall-" +
           std::to_string(
               std::chrono::steady_clock::now().time_since_epoch().count()) +
           ".stl");
      test::write_text(path_, test::ascii_stl(refined_cube(), "cube"));
      text = path_.string();
    }
    std::uint64_t size = text.size();
    HUNDUN_CHECK(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()) ==
                 MPI_SUCCESS);
    HUNDUN_CHECK(size <=
                 static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
    text.resize(static_cast<std::size_t>(size));
    HUNDUN_CHECK(MPI_Bcast(text.data(), static_cast<int>(text.size()), MPI_BYTE,
                           0, mpi.comm()) == MPI_SUCCESS);
    path_ = text;
    HUNDUN_CHECK(MPI_Barrier(mpi.comm()) == MPI_SUCCESS);
  }
  ~FixtureFile() {
    MPI_Barrier(mpi_->comm());
    if (mpi_->rank() == 0) {
      std::error_code error;
      std::filesystem::remove(path_, error);
    }
  }
  const std::filesystem::path &path() const { return path_; }

private:
  const runtime::MpiContext *mpi_;
  std::filesystem::path path_;
};

void run(const runtime::MpiContext &mpi, const std::string &mode,
         const std::string &grid_name) {
  HUNDUN_CHECK(mode == "success" || mode == "warped");
  const bool warped = mode == "warped";
  check_feasible_owner_selector(mpi);
  const auto grid = process_grid(mpi.size(), grid_name);
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, kExtent, {false, false, false}, runtime::DecompositionOptions{grid});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry =
      warped
          ? mesh::MeshGeometry(
                topology, mesh::AnalyticWarpedBoxMapping{{0.0, 0.0, 0.0},
                                                         {1.0, 1.0, 1.0},
                                                         {0.02, -0.015, 0.01}})
          : mesh::MeshGeometry(topology, mesh::UniformBoxMapping{
                                             {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}});
  auto boundaries = boundary::BoundaryRegistry::create(
      config_for(mpi.size(), grid, warped), topology);
  FixtureFile file(mpi);
  const auto surface =
      immersed::ImmersedSurface::load_collective(file.path(), 0.4, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  const auto ghost_plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  const auto plan = immersed::WallQuadraturePlan::create(
      surface, query, domain, topology, geometry, mpi);
  HUNDUN_CHECK(plan.maximum_halo_reach() == 4U);
  const auto owners =
      triangle_owner_map(plan, surface.triangle_count(), mpi);
  runtime::FieldRegistry registry;
  const auto field = registry.declare_field(runtime::FieldDescriptor{
      "wall_q", "1", "task5", runtime::FunctionSpace::cell_average,
      runtime::ScalarType::float64, 1U, 4, false,
      runtime::RestartPolicy::transient, runtime::OutputPolicy::never});
  registry.freeze();
  runtime::FieldStorage storage(registry, decomposition.local_extent());
  auto write = storage.view<double>(field);
  const auto box = decomposition.owned_box();
  const auto local_extent = decomposition.local_extent();
  for (int k = -4; k < local_extent.z + 4; ++k) {
    for (int j = -4; j < local_extent.y + 4; ++j) {
      for (int i = -4; i < local_extent.x + 4; ++i) {
        const runtime::Int3 global{box.begin.x + i, box.begin.y + j,
                                   box.begin.z + k};
        if (global.x < 0 || global.y < 0 || global.z < 0 ||
            global.x >= kExtent.x || global.y >= kExtent.y ||
            global.z >= kExtent.z) {
          write(i, j, k, 0U) = 0.0;
        } else {
          write(i, j, k, 0U) = cell_average(global);
        }
      }
    }
  }
  const auto read =
      static_cast<const runtime::FieldStorage &>(storage).view<double>(field);
  std::uint64_t local_count = plan.local_points().size();
  std::uint64_t global_count = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&local_count, &global_count, 1, MPI_UINT64_T,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_count == 3U * surface.triangle_count());
  double local_weight = 0.0;
  const auto &points = plan.local_points();
  std::map<immersed::ImmersedLinkId, std::vector<mesh::GlobalCellId>>
      donor_authority_by_link;
  for (std::size_t point_index = 0U; point_index < points.size();
       ++point_index) {
    const auto &point = points[point_index];
    const auto &velocity_reconstruction = point.reconstruction;
    const auto pressure_authority_reconstruction =
        immersed::detail::boundary_authority_reconstruction(
            point.reconstruction);
    const auto pressure_link =
        immersed::detail::boundary_authority_link(point.reconstruction);
    const auto &velocity_donors =
        immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
            velocity_reconstruction);
    const auto &pressure_donors =
        immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
            pressure_authority_reconstruction);
    const auto &ghost_donors =
        immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
            ghost_plan.reconstruction(pressure_link));
    HUNDUN_CHECK(donor_authority_equal(velocity_donors, pressure_donors));
    HUNDUN_CHECK(donor_authority_equal(pressure_donors, ghost_donors));
    const auto [authority, inserted] =
        donor_authority_by_link.emplace(pressure_link, velocity_donors);
    if (!inserted)
      HUNDUN_CHECK(donor_authority_equal(authority->second, velocity_donors));
    HUNDUN_CHECK(point.point_index < 3U);
    HUNDUN_CHECK(point.owner_rank == mpi.rank());
    HUNDUN_CHECK(point.reconstruction.quality().rank == 10U);
    HUNDUN_CHECK(point.reconstruction.quality().condition_estimate <= 1.0e8);
    HUNDUN_CHECK(point.reconstruction.quality().halo_reach <= 4U);
    const auto &triangle = surface.triangle(point.triangle);
    constexpr std::array<std::array<double, 3>, 3> barycentric{{
        {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
        {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
        {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
    }};
    runtime::Real3 expected_position{};
    for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
      expected_position.x += barycentric[point.point_index][vertex] *
                             triangle.vertices_m[vertex].x;
      expected_position.y += barycentric[point.point_index][vertex] *
                             triangle.vertices_m[vertex].y;
      expected_position.z += barycentric[point.point_index][vertex] *
                             triangle.vertices_m[vertex].z;
    }
    HUNDUN_CHECK(bits(point.position_m.x) == bits(expected_position.x));
    HUNDUN_CHECK(bits(point.position_m.y) == bits(expected_position.y));
    HUNDUN_CHECK(bits(point.position_m.z) == bits(expected_position.z));
    HUNDUN_CHECK(bits(point.solid_to_fluid_normal.x) ==
                 bits(triangle.geometric_outward_normal.x));
    HUNDUN_CHECK(bits(point.solid_to_fluid_normal.y) ==
                 bits(triangle.geometric_outward_normal.y));
    HUNDUN_CHECK(bits(point.solid_to_fluid_normal.z) ==
                 bits(triangle.geometric_outward_normal.z));
    HUNDUN_CHECK(bits(point.weight_m2) == bits(triangle.area_m2 / 3.0));
    const double selected_amplification = wall_normal_gradient_amplification(
        velocity_reconstruction, point.position_m, point.solid_to_fluid_normal);
    const double alternate_amplification = wall_normal_gradient_amplification(
        point.reconstruction, point.position_m, point.solid_to_fluid_normal);
    const double selection_tolerance = 512.0 *
                                       std::numeric_limits<double>::epsilon() *
                                       std::max(1.0, alternate_amplification);
    HUNDUN_CHECK(selected_amplification <=
                 alternate_amplification + selection_tolerance);
    if (!warped) {
      const double expected_value = polynomial(point.position_m);
      const auto expected_gradient = polynomial_gradient(point.position_m);
      const double tolerance = 4096.0 * std::numeric_limits<double>::epsilon() *
                               std::max(1.0, std::abs(expected_value));
      HUNDUN_CHECK_NEAR(point.reconstruction.value(point.position_m, read, 0U),
                        expected_value, tolerance);
      const auto gradient =
          point.reconstruction.gradient(point.position_m, read, 0U);
      HUNDUN_CHECK_NEAR(gradient.x, expected_gradient.x, 8.0 * tolerance);
      HUNDUN_CHECK_NEAR(gradient.y, expected_gradient.y, 8.0 * tolerance);
      HUNDUN_CHECK_NEAR(gradient.z, expected_gradient.z, 8.0 * tolerance);
    }
    local_weight += point.weight_m2;
  }
  int ranks_with_local_authority = donor_authority_by_link.empty() ? 0 : 1;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &ranks_with_local_authority, 1,
                             MPI_INT, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(ranks_with_local_authority > 0);
  if (!donor_authority_by_link.empty())
    prove_donor_authority_oracle_is_mutation_sensitive(
        donor_authority_by_link.begin()->second);
  double global_weight = local_weight;
  mpi.allreduce_fp64_in_place(&global_weight, 1U,
                              runtime::Fp64ReductionOperation::sum);
  double area = 0.0;
  for (immersed::TriangleId id = 0U; id < surface.triangle_count(); ++id) {
    area += surface.triangle(id).area_m2;
  }
  HUNDUN_CHECK_NEAR(global_weight, area,
                    512.0 * std::numeric_limits<double>::epsilon() *
                        std::max(1.0, std::abs(area)));
  unsigned long long value = plan.fingerprint();
  unsigned long long minimum = 0U;
  unsigned long long maximum = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&value, &minimum, 1, MPI_UNSIGNED_LONG_LONG,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&value, &maximum, 1, MPI_UNSIGNED_LONG_LONG,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(minimum == maximum);
  if (!warped && mpi.size() == 4 && grid_name == "default") {
    const runtime::Int3 alternate_grid{4, 1, 1};
    auto alternate_decomposition = runtime::StructuredDecomposition::create(
        mpi, kExtent, {false, false, false},
        runtime::DecompositionOptions{alternate_grid});
    mesh::MeshTopology alternate_topology(alternate_decomposition);
    mesh::MeshGeometry alternate_geometry(
        alternate_topology,
        mesh::UniformBoxMapping{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}});
    auto alternate_boundaries = boundary::BoundaryRegistry::create(
        config_for(mpi.size(), alternate_grid), alternate_topology);
    const auto alternate_domain = immersed::ImmersedDomain::create(
        surface, query, config::ImmersedFluidSide::outside, alternate_topology,
        alternate_geometry, alternate_boundaries, mpi);
    const auto alternate = immersed::WallQuadraturePlan::create(
        surface, query, alternate_domain, alternate_topology,
        alternate_geometry, mpi);
    const auto alternate_owners =
        triangle_owner_map(alternate, surface.triangle_count(), mpi);
    const auto first_owner_difference = std::mismatch(
        owners.begin(), owners.end(), alternate_owners.begin(),
        alternate_owners.end());
    HUNDUN_CHECK(first_owner_difference.first != owners.end());
    if (mpi.rank() == 0) {
      const auto triangle = static_cast<std::size_t>(
          std::distance(owners.begin(), first_owner_difference.first));
      std::fprintf(stdout,
                   "TASK11_A4_OWNER_MAP first_different_triangle=%zu "
                   "default_owner=%d alternate_owner=%d\n",
                   triangle, owners[triangle], alternate_owners[triangle]);
      std::fflush(stdout);
    }
    HUNDUN_CHECK(alternate.fingerprint() == plan.fingerprint());
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    run(mpi, argc > 1 ? argv[1] : "success", argc > 2 ? argv[2] : "default");
  });
}
