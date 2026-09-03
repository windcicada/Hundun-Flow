// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_ghost_stencil_plan.hpp"

#include "ib_deterministic_qr_detail.hpp"
#include "ib_ghost_stencil_plan_detail.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <system_error>

namespace {

using namespace hundun;
constexpr runtime::Int3 kExtent{12, 12, 12};

std::uint64_t bits(double value) {
  std::uint64_t result = 0U;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

void check_same_constraint(const immersed::AffineGhostConstraint &first,
                           const immersed::AffineGhostConstraint &second) {
  HUNDUN_CHECK(first.link == second.link);
  HUNDUN_CHECK(first.donors.size() == second.donors.size());
  for (std::size_t index = 0U; index < first.donors.size(); ++index) {
    HUNDUN_CHECK(first.donors[index].global_cell ==
                 second.donors[index].global_cell);
    HUNDUN_CHECK(bits(first.donors[index].weight) ==
                 bits(second.donors[index].weight));
  }
  HUNDUN_CHECK(bits(first.wall_value_weight) == bits(second.wall_value_weight));
  HUNDUN_CHECK(bits(first.wall_normal_gradient_weight_m) ==
               bits(second.wall_normal_gradient_weight_m));
}

runtime::Int3 logical_cell(mesh::GlobalCellId id) {
  const auto plane = static_cast<std::uint64_t>(kExtent.x) *
                     static_cast<std::uint64_t>(kExtent.y);
  return {static_cast<int>(id % static_cast<std::uint64_t>(kExtent.x)),
          static_cast<int>((id / static_cast<std::uint64_t>(kExtent.x)) %
                           static_cast<std::uint64_t>(kExtent.y)),
          static_cast<int>(id / plane)};
}

immersed::AffineGhostConstraint cell_average_zero_normal_constraint(
    const immersed::GhostStencilPlan &plan, const immersed::ImmersedLink &link,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry) {
  const auto &reconstruction = plan.reconstruction(link.id);
  const auto average =
      immersed::detail::QuadraticReconstructionWeights::cell_average_weights(
          reconstruction, logical_cell(link.solid_cell), topology, geometry);
  const auto normal = immersed::detail::QuadraticReconstructionWeights::
      directional_gradient_weights(reconstruction, link.wall_intercept_m,
                                   link.solid_to_fluid_normal);
  HUNDUN_CHECK(average.size() == normal.size());
  std::vector<immersed::WeightedDonor> donors(average.size());
  const auto solid = topology.find_local_cell(link.solid_cell);
  HUNDUN_CHECK(solid.has_value());
  const auto center = geometry.cell_center_m(*solid);
  const double distance =
      (center.x - link.wall_intercept_m.x) * link.solid_to_fluid_normal.x +
      (center.y - link.wall_intercept_m.y) * link.solid_to_fluid_normal.y +
      (center.z - link.wall_intercept_m.z) * link.solid_to_fluid_normal.z;
  for (std::size_t index = 0U; index < donors.size(); ++index) {
    HUNDUN_CHECK(average[index].global_cell == normal[index].global_cell);
    donors[index] = {average[index].global_cell,
                     average[index].weight - distance * normal[index].weight};
  }
  return {link.id, std::move(donors), 0.0, distance};
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

double cell_average(mesh::GlobalCellId id) {
  const auto cell = logical_cell(id);
  constexpr double h = 1.0 / static_cast<double>(kExtent.x);
  const runtime::Real3 center{(static_cast<double>(cell.x) + 0.5) * h,
                              (static_cast<double>(cell.y) + 0.5) * h,
                              (static_cast<double>(cell.z) + 0.5) * h};
  return polynomial(center) + (0.7 + 0.3 + 0.9) * h * h / 12.0;
}

double apply_constraint(const immersed::AffineGhostConstraint &constraint,
                        double wall_value, double wall_gradient) {
  double result = constraint.wall_value_weight * wall_value +
                  constraint.wall_normal_gradient_weight_m * wall_gradient;
  for (const auto donor : constraint.donors) {
    result += donor.weight * cell_average(donor.global_cell);
  }
  return result;
}

double apply_extrapolation(const immersed::FluidExtrapolation &extrapolation) {
  double result = 0.0;
  for (const auto donor : extrapolation.donors) {
    result += donor.weight * cell_average(donor.global_cell);
  }
  return result;
}

void check_functional_selection_oracle() {
  using immersed::detail::ReconstructionFunctionalScore;
  const std::vector<ReconstructionFunctionalScore> candidates{
      {3.0, 4.0, 20U, 40U},
      {2.0, 9.0, 30U, 30U},
      {2.0, 3.0, 31U, 20U},
      {2.0, 3.0, 29U, 50U},
      {2.0, 3.0, 29U, 10U}};
  HUNDUN_CHECK(immersed::detail::select_minimum_functional_score(candidates) ==
               4U);

  auto amplification_mutation = candidates;
  amplification_mutation[0].amplification = 1.0;
  HUNDUN_CHECK(immersed::detail::select_minimum_functional_score(
                   amplification_mutation) == 0U);

  auto fingerprint_mutation = candidates;
  fingerprint_mutation[3].pivot_fingerprint = 5U;
  HUNDUN_CHECK(immersed::detail::select_minimum_functional_score(
                   fingerprint_mutation) == 3U);

  for (const auto invalid : {std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN(), -1.0}) {
    auto invalid_candidates = candidates;
    invalid_candidates[2].amplification = invalid;
    bool rejected = false;
    try {
      static_cast<void>(immersed::detail::select_minimum_functional_score(
          invalid_candidates));
    } catch (const runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }

  for (const auto invalid : {std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN(), -1.0}) {
    auto invalid_candidates = candidates;
    invalid_candidates[2].condition_estimate = invalid;
    bool rejected = false;
    try {
      static_cast<void>(immersed::detail::select_minimum_functional_score(
          invalid_candidates));
    } catch (const runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }

  auto zero_donor_count = candidates;
  zero_donor_count[2].donor_count = 0U;
  bool rejected = false;
  try {
    static_cast<void>(
        immersed::detail::select_minimum_functional_score(zero_donor_count));
  } catch (const runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
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

runtime::Real3 midpoint(runtime::Real3 first, runtime::Real3 second) {
  return {(first.x + second.x) * 0.5, (first.y + second.y) * 0.5,
          (first.z + second.z) * 0.5};
}

std::vector<test::StlFixtureTriangle>
partition_conforming_cube(runtime::Real3 translation) {
  auto refined = test::translated(test::outward_cube(), translation);
  for (unsigned level = 0U; level < 2U; ++level) {
    std::vector<test::StlFixtureTriangle> next;
    next.reserve(4U * refined.size());
    for (const auto &triangle : refined) {
      const auto ab = midpoint(triangle.vertices[0], triangle.vertices[1]);
      const auto bc = midpoint(triangle.vertices[1], triangle.vertices[2]);
      const auto ca = midpoint(triangle.vertices[2], triangle.vertices[0]);
      next.push_back({triangle.file_normal, {triangle.vertices[0], ab, ca}});
      next.push_back({triangle.file_normal, {ab, triangle.vertices[1], bc}});
      next.push_back({triangle.file_normal, {ca, bc, triangle.vertices[2]}});
      next.push_back({triangle.file_normal, {ab, bc, ca}});
    }
    refined = std::move(next);
  }
  return refined;
}

config::FlowCaseConfig config_for(int ranks, runtime::Int3 grid,
                                  bool warped = false,
                                  bool periodic_x = false) {
  config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "task5-ghost";
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
    boundary.type = periodic_x && patch < 2U
                        ? config::BoundaryType::periodic
                        : config::BoundaryType::no_slip_wall;
    config.boundaries[patch] = boundary;
  }
  return config;
}

class FixtureFile final {
public:
  FixtureFile(const runtime::MpiContext &mpi, bool aligned,
              bool crosses_periodic_x = false)
      : mpi_(&mpi), scale_(aligned ? 2.0 / 3.0 : 0.4) {
    std::string text;
    if (mpi.rank() == 0) {
      path_ =
          std::filesystem::temp_directory_path() /
          ("hundun-task5-ghost-" +
           std::to_string(
               std::chrono::steady_clock::now().time_since_epoch().count()) +
           ".stl");
      test::write_text(
          path_,
          test::ascii_stl(
              partition_conforming_cube(
                  crosses_periodic_x
                      ? runtime::Real3{-0.25, 0.75, 0.75}
                      : aligned ? runtime::Real3{0.25, 0.25, 0.25}
                                : runtime::Real3{0.75, 0.75, 0.75}),
              "cube"));
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
  double scale() const noexcept { return scale_; }

private:
  const runtime::MpiContext *mpi_;
  std::filesystem::path path_;
  double scale_{};
};

std::uint64_t hash_text(const std::string &text) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const char raw_byte : text) {
    const auto byte = static_cast<unsigned char>(raw_byte);
    hash ^= byte;
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

void check_collective_failure(const runtime::MpiContext &mpi,
                              const std::string &message,
                              const std::string &marker,
                              int expected_failing_rank) {
  HUNDUN_CHECK(message.find(marker) != std::string::npos);
  HUNDUN_CHECK(message.find("lowest failing rank " +
                            std::to_string(expected_failing_rank)) !=
               std::string::npos);
  unsigned long long local = hash_text(message);
  unsigned long long minimum = 0U;
  unsigned long long maximum = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&local, &minimum, 1, MPI_UNSIGNED_LONG_LONG,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&local, &maximum, 1, MPI_UNSIGNED_LONG_LONG,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(minimum == maximum);
}

void prove_periodic_image_moment_oracle_is_mutation_sensitive() {
  constexpr double h = 1.0 / static_cast<double>(kExtent.x);
  constexpr int wrapped_global_cell = kExtent.x - 1;
  const double base_center =
      (static_cast<double>(wrapped_global_cell) + 0.5) * h;
  const double correct_image_center = base_center - 1.0;
  const double unwrapped_center = base_center;
  HUNDUN_CHECK(correct_image_center < 0.0);
  HUNDUN_CHECK(unwrapped_center > 0.0);

  const auto quadratic_cell_moment = [](double center) {
    return center * center + h * h / 12.0;
  };
  const double wrapped_moment = quadratic_cell_moment(correct_image_center);
  const double unwrapped_mutation = quadratic_cell_moment(unwrapped_center);
  const double expected_wrapped_moment = h * h / 3.0;
  HUNDUN_CHECK_NEAR(
      wrapped_moment, expected_wrapped_moment,
      64.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, std::abs(expected_wrapped_moment)));
  HUNDUN_CHECK(std::abs(unwrapped_mutation - wrapped_moment) > 0.5);
}

void run(const runtime::MpiContext &mpi, const std::string &mode,
         const std::string &grid_name) {
  check_functional_selection_oracle();
  const auto grid = process_grid(mpi.size(), grid_name);
  const bool warped = mode == "warped";
  const bool periodic_boundary = mode == "periodic_boundary";
  HUNDUN_CHECK(mode == "success" || mode == "warped" || mode == "failures" ||
               mode == "rank_disagreement" || periodic_boundary);
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, kExtent, {periodic_boundary, false, false},
      runtime::DecompositionOptions{grid});
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
      config_for(mpi.size(), grid, warped, periodic_boundary), topology);
  const bool failures = mode == "failures";
  FixtureFile file(mpi, failures, periodic_boundary);
  const auto surface = immersed::ImmersedSurface::load_collective(
      file.path(), file.scale(), mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  if (periodic_boundary) {
    prove_periodic_image_moment_oracle_is_mutation_sensitive();
    HUNDUN_CHECK(surface.bounding_box_min_m().x < 0.0);
    HUNDUN_CHECK(surface.bounding_box_max_m().x > 0.0);
    std::string message;
    try {
      static_cast<void>(immersed::ImmersedDomain::create(
          surface, query, config::ImmersedFluidSide::outside, topology,
          geometry, boundaries, mpi));
    } catch (const runtime::Error &error) {
      message = error.what();
    }
    check_collective_failure(mpi, message, "active_periodic_pairing failed", 0);

    FixtureFile centred_file(mpi, false);
    const auto centred_surface = immersed::ImmersedSurface::load_collective(
        centred_file.path(), centred_file.scale(), mpi, 0);
    HUNDUN_CHECK(centred_surface.bounding_box_min_m().x > 0.0);
    HUNDUN_CHECK(centred_surface.bounding_box_max_m().x < 1.0);
    const auto centred_query =
        immersed::SurfaceQuery::create(centred_surface);
    const auto centred_domain = immersed::ImmersedDomain::create(
        centred_surface, centred_query, config::ImmersedFluidSide::outside,
        topology, geometry, boundaries, mpi);
    std::uint64_t local_link_count = centred_domain.links().size();
    HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &local_link_count, 1,
                               MPI_UINT64_T, MPI_SUM, mpi.comm()) ==
                 MPI_SUCCESS);
    HUNDUN_CHECK(local_link_count > 0U);
    return;
  }
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  if (mode == "rank_disagreement") {
    HUNDUN_CHECK(mpi.size() >= 2);
    const auto alternate_surface = immersed::ImmersedSurface::load_collective(
        file.path(), file.scale() * 0.99, mpi, 0);
    const auto alternate_query =
        immersed::SurfaceQuery::create(alternate_surface);
    const auto &selected_query = mpi.rank() == 1 ? alternate_query : query;
    std::string message;
    try {
      static_cast<void>(immersed::GhostStencilPlan::create(
          surface, selected_query, domain, topology, geometry, decomposition,
          mpi));
    } catch (const runtime::Error &error) {
      message = error.what();
    }
    check_collective_failure(mpi, message, "rank input signature disagrees", 1);
    return;
  }
  if (failures) {
    std::string message;
    try {
      static_cast<void>(immersed::GhostStencilPlan::create(
          surface, query, domain, topology, geometry, decomposition, mpi));
    } catch (const runtime::Error &error) {
      message = error.what();
    }
    check_collective_failure(mpi, message, "donor selection failed", 0);
    return;
  }
  const auto plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  immersed::detail::ValidatedGeometryScope validated_geometry(topology,
                                                              geometry);
  std::uint64_t local_link_count = domain.links().size();
  std::uint64_t global_link_count = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&local_link_count, &global_link_count, 1,
                             MPI_UINT64_T, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_link_count > 0U);
  for (immersed::ImmersedLinkId id = 0U; id < global_link_count; ++id) {
    HUNDUN_CHECK(plan.velocity_constraint(id, 0U).link == id);
    HUNDUN_CHECK(plan.zero_normal_constraint(id).link == id);
    HUNDUN_CHECK(plan.density_extrapolation(id).link == id);
    const auto quality = plan.reconstruction(id).quality();
    HUNDUN_CHECK(quality.rank == 10U);
    HUNDUN_CHECK(quality.condition_estimate <= 1.0e8);
    HUNDUN_CHECK(quality.halo_reach <= 4U);
  }
  int local_nonpositive_density_sample = 0;
  std::uint32_t local_maximum_reach = 0U;
  for (const auto &link : domain.links()) {
    const auto &velocity = plan.velocity_constraint(link.id, 0U);
    HUNDUN_CHECK(velocity.link == link.id);
    HUNDUN_CHECK(velocity.donors.size() >= 14U);
    HUNDUN_CHECK(velocity.donors.size() <= 32U);
    HUNDUN_CHECK(velocity.wall_value_weight == 1.0);
    HUNDUN_CHECK(velocity.wall_normal_gradient_weight_m == 0.0);
    const auto anchor = logical_cell(link.fluid_cell);
    for (const auto donor : velocity.donors) {
      const auto logical = logical_cell(donor.global_cell);
      const auto distance = [](int first, int second) {
        const auto difference = static_cast<std::int64_t>(first) -
                                static_cast<std::int64_t>(second);
        return static_cast<std::uint32_t>(difference < 0 ? -difference
                                                         : difference);
      };
      local_maximum_reach = std::max(local_maximum_reach,
                                     std::max({distance(anchor.x, logical.x),
                                               distance(anchor.y, logical.y),
                                               distance(anchor.z, logical.z)}));
      for (const auto other : velocity.donors) {
        const auto other_logical = logical_cell(other.global_cell);
        local_maximum_reach =
            std::max(local_maximum_reach,
                     std::max({distance(other_logical.x, logical.x),
                               distance(other_logical.y, logical.y),
                               distance(other_logical.z, logical.z)}));
      }
    }
    check_same_constraint(velocity, plan.velocity_constraint(link.id, 1U));
    check_same_constraint(velocity, plan.velocity_constraint(link.id, 2U));
    const auto &normal = plan.zero_normal_constraint(link.id);
    const auto average_normal =
        cell_average_zero_normal_constraint(plan, link, topology, geometry);
    HUNDUN_CHECK(normal.wall_value_weight == 0.0);
    HUNDUN_CHECK(average_normal.wall_value_weight == 0.0);
    double pressure_row_constraint_coefficient = 0.0;
    for (const auto &donor : average_normal.donors)
      if (donor.global_cell != link.fluid_cell)
        pressure_row_constraint_coefficient += donor.weight;
    HUNDUN_CHECK(std::isfinite(pressure_row_constraint_coefficient));
    HUNDUN_CHECK(pressure_row_constraint_coefficient > 0.0);
    const auto solid = topology.find_local_cell(link.solid_cell);
    HUNDUN_CHECK(solid.has_value());
    const auto solid_center = geometry.cell_center_m(*solid);
    const double expected_d_g = (solid_center.x - link.wall_intercept_m.x) *
                                    link.solid_to_fluid_normal.x +
                                (solid_center.y - link.wall_intercept_m.y) *
                                    link.solid_to_fluid_normal.y +
                                (solid_center.z - link.wall_intercept_m.z) *
                                    link.solid_to_fluid_normal.z;
    HUNDUN_CHECK(expected_d_g < 0.0);
    HUNDUN_CHECK_NEAR(normal.wall_normal_gradient_weight_m, expected_d_g,
                      64.0 * std::numeric_limits<double>::epsilon() *
                          std::max(1.0, std::abs(expected_d_g)));
    HUNDUN_CHECK_NEAR(average_normal.wall_normal_gradient_weight_m,
                      expected_d_g,
                      64.0 * std::numeric_limits<double>::epsilon() *
                          std::max(1.0, std::abs(expected_d_g)));
    if (!warped) {
      const auto wall_gradient = polynomial_gradient(link.wall_intercept_m);
      const double wall_normal_gradient =
          wall_gradient.x * link.solid_to_fluid_normal.x +
          wall_gradient.y * link.solid_to_fluid_normal.y +
          wall_gradient.z * link.solid_to_fluid_normal.z;
      const double expected_ghost = polynomial(solid_center);
      double weight_l1 = 1.0;
      for (const auto donor : velocity.donors) {
        weight_l1 += std::abs(donor.weight);
      }
      const double tolerance =
          4096.0 * std::numeric_limits<double>::epsilon() *
          std::max({1.0, std::abs(expected_ghost), weight_l1});
      HUNDUN_CHECK_NEAR(
          apply_constraint(velocity, polynomial(link.wall_intercept_m), 0.0),
          expected_ghost, tolerance);
      HUNDUN_CHECK_NEAR(apply_constraint(normal, 0.0, wall_normal_gradient),
                        expected_ghost, tolerance);
      const double expected_ghost_average = cell_average(link.solid_cell);
      HUNDUN_CHECK_NEAR(
          apply_constraint(average_normal, 0.0, wall_normal_gradient),
          expected_ghost_average,
          4096.0 * std::numeric_limits<double>::epsilon() *
              std::max({1.0, std::abs(expected_ghost_average), weight_l1}));
      HUNDUN_CHECK(expected_ghost_average != expected_ghost);
      HUNDUN_CHECK_NEAR(
          apply_extrapolation(plan.density_extrapolation(link.id)),
          expected_ghost, tolerance);
    }
    HUNDUN_CHECK(plan.density_extrapolation(link.id).donors.size() ==
                 velocity.donors.size());
    HUNDUN_CHECK(plan.density_extrapolation(link.id).donors.size() <= 32U);
    const auto &density = plan.density_extrapolation(link.id);
    const auto minimum_weight =
        std::min_element(density.donors.begin(), density.donors.end(),
                         [](const auto &first, const auto &second) {
                           return first.weight < second.weight;
                         });
    if (minimum_weight != density.donors.end() &&
        minimum_weight->weight < 0.0) {
      double unit_result = 0.0;
      for (const auto donor : density.donors) {
        unit_result += donor.weight;
      }
      const double selected_value =
          std::max(1.0, (unit_result + 1.0) / -minimum_weight->weight) + 1.0;
      double result = 0.0;
      for (const auto &donor : density.donors) {
        const double positive_value =
            &donor == &*minimum_weight ? selected_value : 1.0;
        HUNDUN_CHECK(positive_value > 0.0);
        result += donor.weight * positive_value;
      }
      HUNDUN_CHECK(result < 0.0);
      local_nonpositive_density_sample = 1;
    }
    HUNDUN_CHECK(plan.reconstruction(link.id).quality().rank == 10U);
  }
  int global_nonpositive_density_sample = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_nonpositive_density_sample,
                             &global_nonpositive_density_sample, 1, MPI_INT,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_nonpositive_density_sample == 1);
  std::uint32_t global_maximum_reach = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&local_maximum_reach, &global_maximum_reach, 1,
                             MPI_UINT32_T, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(plan.maximum_halo_reach() == global_maximum_reach);
  HUNDUN_CHECK(plan.maximum_halo_reach() <= 4U);
  unsigned long long value = plan.fingerprint();
  unsigned long long minimum = 0U;
  unsigned long long maximum = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&value, &minimum, 1, MPI_UNSIGNED_LONG_LONG,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&value, &maximum, 1, MPI_UNSIGNED_LONG_LONG,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(minimum == maximum);
  bool rejected_component = false;
  try {
    static_cast<void>(plan.velocity_constraint(0U, 3U));
  } catch (const runtime::Error &) {
    rejected_component = true;
  }
  HUNDUN_CHECK(rejected_component);
  bool rejected_link = false;
  try {
    static_cast<void>(plan.reconstruction(global_link_count));
  } catch (const runtime::Error &) {
    rejected_link = true;
  }
  HUNDUN_CHECK(rejected_link);
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
    const auto alternate = immersed::GhostStencilPlan::create(
        surface, query, alternate_domain, alternate_topology,
        alternate_geometry, alternate_decomposition, mpi);
    HUNDUN_CHECK(alternate.fingerprint() == plan.fingerprint());
    HUNDUN_CHECK(alternate.maximum_halo_reach() == plan.maximum_halo_reach());
    for (immersed::ImmersedLinkId id = 0U; id < global_link_count; ++id) {
      check_same_constraint(plan.velocity_constraint(id, 0U),
                            alternate.velocity_constraint(id, 0U));
      check_same_constraint(plan.zero_normal_constraint(id),
                            alternate.zero_normal_constraint(id));
    }
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
