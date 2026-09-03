// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_wall_force.hpp"

#include "src/ib_quadratic_reconstruction_detail.hpp"
#include "src/ib_wall_force_detail.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"
#include "src/ib_deterministic_qr_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using namespace hundun;

constexpr runtime::Int3 kExtent{12, 12, 12};
constexpr int kGhost = 4;
constexpr std::array<double, 9> kGradient{1.0,  2.0,  -1.0, 0.5, -0.25,
                                          0.75, -1.5, 0.2,  0.1};
constexpr runtime::Real3 kPressureGradient{0.7, -0.4, 0.3};
constexpr runtime::Real3 kViscosityGradient{0.4, -0.2, 0.1};

std::uint64_t bits(double value) {
  std::uint64_t encoded{};
  std::memcpy(&encoded, &value, sizeof(encoded));
  return encoded;
}

runtime::Real3 add(runtime::Real3 lhs, runtime::Real3 rhs) {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

runtime::Real3 scale(double factor, runtime::Real3 value) {
  return {factor * value.x, factor * value.y, factor * value.z};
}

runtime::Real3 subtract(runtime::Real3 lhs, runtime::Real3 rhs) {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

runtime::Real3 cross(runtime::Real3 lhs, runtime::Real3 rhs) {
  return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

double dot(runtime::Real3 lhs, runtime::Real3 rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

runtime::Real3 midpoint(runtime::Real3 lhs, runtime::Real3 rhs) {
  return scale(0.5, add(lhs, rhs));
}

runtime::Real3 project_to_sphere(runtime::Real3 point, runtime::Real3 center,
                                 double radius) {
  const auto radial = subtract(point, center);
  const double length = std::sqrt(dot(radial, radial));
  HUNDUN_CHECK(length > 0.0);
  return add(center, scale(radius / length, radial));
}

test::StlFixtureTriangle sphere_triangle(runtime::Real3 first,
                                         runtime::Real3 second,
                                         runtime::Real3 third,
                                         runtime::Real3 center) {
  auto normal = cross(subtract(second, first), subtract(third, first));
  const auto centroid = scale(1.0 / 3.0, add(add(first, second), third));
  if (dot(normal, subtract(centroid, center)) < 0.0) {
    std::swap(second, third);
    normal = cross(subtract(second, first), subtract(third, first));
  }
  const double length = std::sqrt(dot(normal, normal));
  HUNDUN_CHECK(length > 0.0);
  return {scale(1.0 / length, normal), {first, second, third}};
}

std::vector<test::StlFixtureTriangle> faceted_sphere() {
  constexpr runtime::Real3 center{1.25, 1.25, 1.25};
  constexpr double radius = 0.45;
  const runtime::Real3 px{center.x + radius, center.y, center.z};
  const runtime::Real3 nx{center.x - radius, center.y, center.z};
  const runtime::Real3 py{center.x, center.y + radius, center.z};
  const runtime::Real3 ny{center.x, center.y - radius, center.z};
  const runtime::Real3 pz{center.x, center.y, center.z + radius};
  const runtime::Real3 nz{center.x, center.y, center.z - radius};
  std::vector<test::StlFixtureTriangle> triangles{
      sphere_triangle(pz, px, py, center), sphere_triangle(pz, py, nx, center),
      sphere_triangle(pz, nx, ny, center), sphere_triangle(pz, ny, px, center),
      sphere_triangle(nz, py, px, center), sphere_triangle(nz, nx, py, center),
      sphere_triangle(nz, ny, nx, center), sphere_triangle(nz, px, ny, center)};
  for (int level = 0; level < 2; ++level) {
    std::vector<test::StlFixtureTriangle> refined;
    refined.reserve(4U * triangles.size());
    for (const auto &triangle : triangles) {
      const auto ab = project_to_sphere(
          midpoint(triangle.vertices[0], triangle.vertices[1]), center, radius);
      const auto bc = project_to_sphere(
          midpoint(triangle.vertices[1], triangle.vertices[2]), center, radius);
      const auto ca = project_to_sphere(
          midpoint(triangle.vertices[2], triangle.vertices[0]), center, radius);
      refined.push_back(sphere_triangle(triangle.vertices[0], ab, ca, center));
      refined.push_back(sphere_triangle(ab, triangle.vertices[1], bc, center));
      refined.push_back(sphere_triangle(ca, bc, triangle.vertices[2], center));
      refined.push_back(sphere_triangle(ab, bc, ca, center));
    }
    triangles = std::move(refined);
  }
  return triangles;
}

std::vector<test::StlFixtureTriangle>
refined_cube(runtime::Real3 translation = {0.75, 0.75, 0.75}) {
  const auto coarse = test::translated(test::outward_cube(), translation);
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

runtime::Int3 process_grid(int ranks, bool alternate) {
  if (ranks == 1) {
    return {1, 1, 1};
  }
  if (ranks == 2) {
    return alternate ? runtime::Int3{1, 2, 1} : runtime::Int3{2, 1, 1};
  }
  return alternate ? runtime::Int3{4, 1, 1} : runtime::Int3{2, 2, 1};
}

config::FlowCaseConfig config_for(int ranks, runtime::Int3 grid) {
  config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "task10-wall-force";
  config.simulation_type = config::SimulationType::variable_density_flow;
  config.density_model = config::DensityModel::constant;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = grid;
  config.mesh.cells = kExtent;
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.mesh.mapping = config::MeshMapping::uniform_box;
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
  explicit FixtureFile(const runtime::MpiContext &mpi, bool permuted = false,
                       runtime::Real3 translation = {0.75, 0.75, 0.75},
                       std::vector<test::StlFixtureTriangle> triangles = {})
      : mpi_(&mpi) {
    std::string text;
    if (mpi.rank() == 0) {
      path_ =
          std::filesystem::temp_directory_path() /
          ("hundun-task10-wall-force-" +
           std::to_string(
               std::chrono::steady_clock::now().time_since_epoch().count()) +
           ".stl");
      if (triangles.empty())
        triangles = refined_cube(translation);
      if (permuted) {
        std::reverse(triangles.begin(), triangles.end());
      }
      test::write_text(path_, test::ascii_stl(triangles, "cube"));
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

runtime::FieldDescriptor descriptor(std::string name,
                                    std::uint32_t components,
                                    int ghost_width = kGhost) {
  return {std::move(name),
          "1",
          "task10",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost_width,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

struct FieldIds final {
  runtime::FieldId pressure{};
  runtime::FieldId velocity{};
  runtime::FieldId gradient{};
  runtime::FieldId viscosity{};
  runtime::FieldId bad_pressure{};
};

FieldIds declare_fields(runtime::FieldRegistry &registry,
                        int ghost_width = kGhost) {
  return {registry.declare_field(descriptor("pi", 1U, ghost_width)),
          registry.declare_field(descriptor("u", 3U, ghost_width)),
          registry.declare_field(descriptor("grad_u", 9U, ghost_width)),
          registry.declare_field(descriptor("mu_eff", 1U, ghost_width)),
          registry.declare_field(descriptor("bad_pi", 2U, ghost_width))};
}

runtime::Real3 cell_center(runtime::Int3 global) {
  constexpr double h = 1.0 / static_cast<double>(kExtent.x);
  return {(static_cast<double>(global.x) + 0.5) * h,
          (static_cast<double>(global.y) + 0.5) * h,
          (static_cast<double>(global.z) + 0.5) * h};
}

double pressure(runtime::Real3 point, double offset = 5.0) {
  return offset + dot(kPressureGradient, point);
}

double quadratic_pressure(runtime::Real3 point) {
  return pressure(point) + 0.2 * point.x * point.x - 0.15 * point.y * point.z +
         0.1 * point.z * point.z;
}

runtime::Real3 quadratic_pressure_gradient(runtime::Real3 point) {
  return {kPressureGradient.x + 0.4 * point.x,
          kPressureGradient.y - 0.15 * point.z,
          kPressureGradient.z - 0.15 * point.y + 0.2 * point.z};
}

double pressure_cell_average(runtime::Int3 global, bool cubic) {
  constexpr double h = 1.0 / static_cast<double>(kExtent.x);
  const auto point = cell_center(global);
  const double square_offset = h * h / 12.0;
  if (!cubic)
    return pressure(point) + 0.2 * (point.x * point.x + square_offset) -
           0.15 * point.y * point.z + 0.1 * (point.z * point.z + square_offset);
  const double cube_offset = h * h / 4.0;
  return pressure(point) +
         0.8 * (point.x * point.x * point.x + point.x * cube_offset) -
         0.5 * (point.y * point.y * point.y + point.y * cube_offset) +
         0.3 * (point.z * point.z * point.z + point.z * cube_offset) +
         0.2 * point.x * point.y * point.z;
}

double viscosity(runtime::Real3 point) {
  return 2.0 + dot(kViscosityGradient, point);
}

double velocity(runtime::Real3 point, std::size_t component) {
  const std::array<double, 3> coordinate{point.x, point.y, point.z};
  double result = 0.0;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    result += kGradient[component * 3U + axis] * coordinate[axis];
  }
  return result;
}

bool inside(runtime::Int3 global) {
  return global.x >= 0 && global.y >= 0 && global.z >= 0 &&
         global.x < kExtent.x && global.y < kExtent.y && global.z < kExtent.z;
}

runtime::Int3 logical_cell(mesh::GlobalCellId id) {
  const auto nx = static_cast<std::uint64_t>(kExtent.x);
  const auto ny = static_cast<std::uint64_t>(kExtent.y);
  const auto plane = nx * ny;
  HUNDUN_CHECK(id < plane * static_cast<std::uint64_t>(kExtent.z));
  return {static_cast<int>(id % nx), static_cast<int>((id / nx) % ny),
          static_cast<int>(id / plane)};
}

bool inside(runtime::Box3 box, runtime::Int3 cell) {
  return cell.x >= box.begin.x && cell.x < box.end.x &&
         cell.y >= box.begin.y && cell.y < box.end.y &&
         cell.z >= box.begin.z && cell.z < box.end.z;
}

void fill_fields(runtime::FieldStorage &storage, const FieldIds &ids,
                 const runtime::StructuredDecomposition &decomposition,
                 double pressure_offset = 5.0,
                 double pressure_gradient_scale = 1.0,
                 int ghost_width = kGhost) {
  auto pi = storage.view<double>(ids.pressure);
  auto u = storage.view<double>(ids.velocity);
  auto grad = storage.view<double>(ids.gradient);
  auto mu = storage.view<double>(ids.viscosity);
  auto bad = storage.view<double>(ids.bad_pressure);
  const auto local = decomposition.local_extent();
  const auto box = decomposition.owned_box();
  for (int k = -ghost_width; k < local.z + ghost_width; ++k) {
    for (int j = -ghost_width; j < local.y + ghost_width; ++j) {
      for (int i = -ghost_width; i < local.x + ghost_width; ++i) {
        const runtime::Int3 global{box.begin.x + i, box.begin.y + j,
                                   box.begin.z + k};
        const auto point = cell_center(global);
        const bool valid = inside(global);
        pi(i, j, k, 0) =
            valid ? pressure_offset +
                        pressure_gradient_scale * dot(kPressureGradient, point)
                  : 0.0;
        mu(i, j, k, 0) = valid ? viscosity(point) : 0.0;
        for (int component = 0; component < 3; ++component) {
          u(i, j, k, component) =
              valid ? velocity(point, static_cast<std::size_t>(component))
                    : 0.0;
        }
        for (int component = 0; component < 9; ++component) {
          grad(i, j, k, component) =
              valid ? kGradient[static_cast<std::size_t>(component)] : 0.0;
        }
        bad(i, j, k, 0) = pi(i, j, k, 0);
        bad(i, j, k, 1) = pi(i, j, k, 0);
      }
    }
  }
}

void fill_polynomial_pressure(
    runtime::FieldStorage &storage, runtime::FieldId pressure_id,
    const runtime::StructuredDecomposition &decomposition, bool cubic) {
  auto field = storage.view<double>(pressure_id);
  const auto local = decomposition.local_extent();
  const auto box = decomposition.owned_box();
  for (int k = -kGhost; k < local.z + kGhost; ++k)
    for (int j = -kGhost; j < local.y + kGhost; ++j)
      for (int i = -kGhost; i < local.x + kGhost; ++i) {
        const runtime::Int3 global{box.begin.x + i, box.begin.y + j,
                                   box.begin.z + k};
        field(i, j, k, 0) =
            inside(global) ? pressure_cell_average(global, cubic) : 0.0;
      }
}

template <class Mutation>
void mutate_component(runtime::FieldStorage &storage, runtime::FieldId id,
                      const runtime::StructuredDecomposition &decomposition,
                      int component, Mutation mutation) {
  auto field = storage.view<double>(id);
  const auto local = decomposition.local_extent();
  for (int k = -kGhost; k < local.z + kGhost; ++k) {
    for (int j = -kGhost; j < local.y + kGhost; ++j) {
      for (int i = -kGhost; i < local.x + kGhost; ++i) {
        field(i, j, k, component) = mutation(field(i, j, k, component));
      }
    }
  }
}

struct Views final {
  runtime::FieldView<const double> pressure;
  runtime::FieldView<const double> velocity;
  runtime::FieldView<const double> gradient;
  runtime::FieldView<const double> viscosity;
  runtime::FieldView<const double> bad_pressure;
};

Views views(const runtime::FieldStorage &storage, const FieldIds &ids) {
  return {
      storage.view<double>(ids.pressure), storage.view<double>(ids.velocity),
      storage.view<double>(ids.gradient), storage.view<double>(ids.viscosity),
      storage.view<double>(ids.bad_pressure)};
}

std::vector<std::uint64_t> snapshot(const runtime::FieldStorage &storage,
                                    const FieldIds &ids, runtime::Int3 local) {
  const auto read = views(storage, ids);
  std::vector<std::uint64_t> result;
  const std::array<std::pair<runtime::FieldView<const double>, int>, 5> fields{
      std::pair{read.pressure, 1}, std::pair{read.velocity, 3},
      std::pair{read.gradient, 9}, std::pair{read.viscosity, 1},
      std::pair{read.bad_pressure, 2}};
  for (const auto &field : fields) {
    for (int k = -kGhost; k < local.z + kGhost; ++k) {
      for (int j = -kGhost; j < local.y + kGhost; ++j) {
        for (int i = -kGhost; i < local.x + kGhost; ++i) {
          for (int component = 0; component < field.second; ++component) {
            result.push_back(bits(field.first(i, j, k, component)));
          }
        }
      }
    }
  }
  return result;
}

runtime::Real3 viscous_traction(runtime::Real3 normal, double mu,
                                const std::array<double, 9> &gradient) {
  const double divergence = gradient[0] + gradient[4] + gradient[8];
  const std::array<double, 3> n{normal.x, normal.y, normal.z};
  std::array<double, 3> result{};
  for (std::size_t row = 0U; row < 3U; ++row) {
    for (std::size_t column = 0U; column < 3U; ++column) {
      const double symmetric = gradient[row * 3U + column] +
                               gradient[column * 3U + row] -
                               (row == column ? (2.0 / 3.0) * divergence : 0.0);
      result[row] += mu * symmetric * n[column];
    }
  }
  return {result[0], result[1], result[2]};
}

struct Oracle final {
  immersed::ForceComponents force;
  immersed::MomentComponents moment;
  runtime::Real3 unprojected_viscous_N{};
  runtime::Real3 closure{};
  std::uint64_t count{};
};

Oracle oracle(const immersed::WallQuadraturePlan &plan,
              const runtime::MpiContext &mpi,
              const runtime::FieldView<const double> &velocity_field,
              double pressure_offset = 5.0) {
  Oracle result;
  const auto &points = plan.local_points();
  for (std::size_t index = 0U; index < points.size(); ++index) {
    const auto &point = points[index];
    const auto &velocity_reconstruction = point.reconstruction;
    const auto pressure_force =
        scale(-pressure(point.position_m, pressure_offset) * point.weight_m2,
              point.solid_to_fluid_normal);
    std::array<double, 9> gradient{};
    std::array<double, 9> unprojected_gradient{};
    for (std::size_t component = 0U; component < 3U; ++component) {
      const auto value = immersed::detail::gradient_with_origin_constraint(
          velocity_reconstruction, point.position_m, velocity_field, component,
          0.0);
      unprojected_gradient[component * 3U] = value.x;
      unprojected_gradient[component * 3U + 1U] = value.y;
      unprojected_gradient[component * 3U + 2U] = value.z;
      const double normal_derivative = value.x * point.solid_to_fluid_normal.x +
                                       value.y * point.solid_to_fluid_normal.y +
                                       value.z * point.solid_to_fluid_normal.z;
      gradient[component * 3U] =
          normal_derivative * point.solid_to_fluid_normal.x;
      gradient[component * 3U + 1U] =
          normal_derivative * point.solid_to_fluid_normal.y;
      gradient[component * 3U + 2U] =
          normal_derivative * point.solid_to_fluid_normal.z;
    }
    const auto viscous_force =
        scale(point.weight_m2,
              viscous_traction(point.solid_to_fluid_normal,
                               viscosity(point.position_m), gradient));
    const auto unprojected_viscous_force =
        scale(point.weight_m2, viscous_traction(point.solid_to_fluid_normal,
                                                viscosity(point.position_m),
                                                unprojected_gradient));
    const auto total_force = add(pressure_force, viscous_force);
    result.force.pressure_N = add(result.force.pressure_N, pressure_force);
    result.force.viscous_N = add(result.force.viscous_N, viscous_force);
    result.force.total_N = add(result.force.total_N, total_force);
    result.unprojected_viscous_N =
        add(result.unprojected_viscous_N, unprojected_viscous_force);
    result.moment.pressure_N_m = add(result.moment.pressure_N_m,
                                     cross(point.position_m, pressure_force));
    result.moment.viscous_N_m =
        add(result.moment.viscous_N_m, cross(point.position_m, viscous_force));
    result.moment.total_N_m =
        add(result.moment.total_N_m, cross(point.position_m, total_force));
    result.closure = add(result.closure,
                         scale(point.weight_m2, point.solid_to_fluid_normal));
    ++result.count;
  }
  std::array<double, 24> values{result.force.pressure_N.x,
                                result.force.pressure_N.y,
                                result.force.pressure_N.z,
                                result.force.viscous_N.x,
                                result.force.viscous_N.y,
                                result.force.viscous_N.z,
                                result.force.total_N.x,
                                result.force.total_N.y,
                                result.force.total_N.z,
                                result.moment.pressure_N_m.x,
                                result.moment.pressure_N_m.y,
                                result.moment.pressure_N_m.z,
                                result.moment.viscous_N_m.x,
                                result.moment.viscous_N_m.y,
                                result.moment.viscous_N_m.z,
                                result.moment.total_N_m.x,
                                result.moment.total_N_m.y,
                                result.moment.total_N_m.z,
                                result.unprojected_viscous_N.x,
                                result.unprojected_viscous_N.y,
                                result.unprojected_viscous_N.z,
                                result.closure.x,
                                result.closure.y,
                                result.closure.z};
  mpi.allreduce_fp64_in_place(values.data(), values.size(),
                              runtime::Fp64ReductionOperation::sum);
  std::uint64_t count = result.count;
  HUNDUN_CHECK(MPI_Allreduce(&count, &result.count, 1, MPI_UINT64_T, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  const auto assign = [&](runtime::Real3 &value, std::size_t index) {
    value = {values[index], values[index + 1U], values[index + 2U]};
  };
  assign(result.force.pressure_N, 0U);
  assign(result.force.viscous_N, 3U);
  assign(result.force.total_N, 6U);
  assign(result.moment.pressure_N_m, 9U);
  assign(result.moment.viscous_N_m, 12U);
  assign(result.moment.total_N_m, 15U);
  assign(result.unprojected_viscous_N, 18U);
  assign(result.closure, 21U);
  return result;
}

runtime::Real3 reconstructed_pressure_force(
    const immersed::WallQuadraturePlan &plan, const runtime::MpiContext &mpi,
    const runtime::FieldView<const double> &pressure_field,
    bool coupled_normal_constraint) {
  runtime::Real3 result{};
  const auto &points = plan.local_points();
  for (std::size_t index = 0U; index < points.size(); ++index) {
    const auto &point = points[index];
    double value =
        point.reconstruction.value(point.position_m, pressure_field, 0U);
    if (coupled_normal_constraint) {
      const auto gradient =
          point.reconstruction.gradient(point.position_m, pressure_field, 0U);
      const double normal_gradient = dot(gradient, point.solid_to_fluid_normal);
      value = immersed::detail::value_with_origin_normal_gradient(
          point.reconstruction, point.position_m, pressure_field, 0U,
          normal_gradient);
    }
    result = add(result,
                 scale(-value * point.weight_m2, point.solid_to_fluid_normal));
  }
  std::array<double, 3> reduced{result.x, result.y, result.z};
  mpi.allreduce_fp64_in_place(reduced.data(), reduced.size(),
                              runtime::Fp64ReductionOperation::sum);
  return {reduced[0], reduced[1], reduced[2]};
}

template <class Pressure>
runtime::Real3 sampled_pressure_force(const immersed::WallQuadraturePlan &plan,
                                      const runtime::MpiContext &mpi,
                                      Pressure evaluate) {
  runtime::Real3 result{};
  for (const auto &point : plan.local_points())
    result = add(result, scale(-evaluate(point.position_m) * point.weight_m2,
                               point.solid_to_fluid_normal));
  std::array<double, 3> reduced{result.x, result.y, result.z};
  mpi.allreduce_fp64_in_place(reduced.data(), reduced.size(),
                              runtime::Fp64ReductionOperation::sum);
  return {reduced[0], reduced[1], reduced[2]};
}

void check_near(runtime::Real3 actual, runtime::Real3 expected,
                double tolerance) {
  HUNDUN_CHECK_NEAR(actual.x, expected.x, tolerance);
  HUNDUN_CHECK_NEAR(actual.y, expected.y, tolerance);
  HUNDUN_CHECK_NEAR(actual.z, expected.z, tolerance);
}

void check_bits(runtime::Real3 lhs, runtime::Real3 rhs) {
  HUNDUN_CHECK(bits(lhs.x) == bits(rhs.x));
  HUNDUN_CHECK(bits(lhs.y) == bits(rhs.y));
  HUNDUN_CHECK(bits(lhs.z) == bits(rhs.z));
}

void check_same(const immersed::WallForceSample &lhs,
                const immersed::WallForceSample &rhs) {
  check_bits(lhs.surface_traction.pressure_N, rhs.surface_traction.pressure_N);
  check_bits(lhs.surface_traction.viscous_N, rhs.surface_traction.viscous_N);
  check_bits(lhs.surface_traction.total_N, rhs.surface_traction.total_N);
  check_bits(lhs.moment_about_global_origin.pressure_N_m,
             rhs.moment_about_global_origin.pressure_N_m);
  check_bits(lhs.moment_about_global_origin.viscous_N_m,
             rhs.moment_about_global_origin.viscous_N_m);
  check_bits(lhs.moment_about_global_origin.total_N_m,
             rhs.moment_about_global_origin.total_N_m);
  check_bits(lhs.area_vector_closure_m2, rhs.area_vector_closure_m2);
  HUNDUN_CHECK(lhs.quadrature_point_count == rhs.quadrature_point_count);
  HUNDUN_CHECK(lhs.lowest_failing_rank == rhs.lowest_failing_rank);
}

void check_failed(const immersed::WallForceSample &sample, int rank) {
  HUNDUN_CHECK(sample.lowest_failing_rank == rank);
  HUNDUN_CHECK(sample.quadrature_point_count == 0U);
  const std::array<runtime::Real3, 7> values{
      sample.surface_traction.pressure_N,
      sample.surface_traction.viscous_N,
      sample.surface_traction.total_N,
      sample.moment_about_global_origin.pressure_N_m,
      sample.moment_about_global_origin.viscous_N_m,
      sample.moment_about_global_origin.total_N_m,
      sample.area_vector_closure_m2};
  for (const auto value : values) {
    HUNDUN_CHECK(bits(value.x) == 0U);
    HUNDUN_CHECK(bits(value.y) == 0U);
    HUNDUN_CHECK(bits(value.z) == 0U);
  }
}

double max_norm(runtime::Real3 value) {
  return std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
}

bool finite(runtime::Real3 value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

bool finite(const immersed::WallForceSample &sample) {
  return finite(sample.surface_traction.pressure_N) &&
         finite(sample.surface_traction.viscous_N) &&
         finite(sample.surface_traction.total_N) &&
         finite(sample.moment_about_global_origin.pressure_N_m) &&
         finite(sample.moment_about_global_origin.viscous_N_m) &&
         finite(sample.moment_about_global_origin.total_N_m) &&
         finite(sample.area_vector_closure_m2);
}

void run(const runtime::MpiContext &mpi, bool alternate) {
  const auto grid = process_grid(mpi.size(), alternate);
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, kExtent, {false, false, false}, runtime::DecompositionOptions{grid});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}});
  auto boundaries = boundary::BoundaryRegistry::create(
      config_for(mpi.size(), grid), topology);
  FixtureFile file(mpi);
  const auto surface =
      immersed::ImmersedSurface::load_collective(file.path(), 0.4, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  const auto plan = immersed::WallQuadraturePlan::create(
      surface, query, domain, topology, geometry, mpi);
  const auto integrator = immersed::WallForceIntegrator::create(plan, mpi);
  const auto plan_fingerprint = plan.fingerprint();

  runtime::FieldRegistry registry;
  const auto ids = declare_fields(registry);
  registry.freeze();
  runtime::FieldStorage storage(registry, decomposition.local_extent());
  fill_fields(storage, ids, decomposition);
  auto read = views(storage, ids);
  const auto before = snapshot(storage, ids, decomposition.local_extent());
  const auto sample = integrator.integrate(read.pressure, read.velocity,
                                           read.gradient, read.viscosity);
  HUNDUN_CHECK(sample.lowest_failing_rank == -1);
  HUNDUN_CHECK(snapshot(storage, ids, decomposition.local_extent()) == before);

  if (mpi.size() == 1 && !alternate) {
    const immersed::WallQuadraturePoint *selected_point = nullptr;
    immersed::WeightedDonor selected_donor{};
    for (const auto &point : plan.local_points()) {
      const auto weights =
          immersed::detail::QuadraticReconstructionWeights::value_weights(
              point.reconstruction, point.position_m);
      const auto found = std::min_element(
          weights.begin(), weights.end(),
          [](const auto &left, const auto &right) {
            return left.weight < right.weight;
          });
      if (found != weights.end() && found->weight < -1.0e-6) {
        selected_point = &point;
        selected_donor = *found;
        break;
      }
    }
    HUNDUN_CHECK(selected_point != nullptr);
    auto mu = storage.view<double>(ids.viscosity);
    const auto local = decomposition.local_extent();
    const auto box = decomposition.owned_box();
    for (int k = -kGhost; k < local.z + kGhost; ++k)
      for (int j = -kGhost; j < local.y + kGhost; ++j)
        for (int i = -kGhost; i < local.x + kGhost; ++i)
          mu(i, j, k, 0) = inside({box.begin.x + i, box.begin.y + j,
                                   box.begin.z + k})
                               ? 1.0
                               : 0.0;
    const auto donor = logical_cell(selected_donor.global_cell);
    mu(donor.x - box.begin.x, donor.y - box.begin.y,
       donor.z - box.begin.z, 0) =
        1.0 + 2.0 / std::abs(selected_donor.weight);
    read = views(storage, ids);
    const double unbounded = selected_point->reconstruction.value(
        selected_point->position_m, read.viscosity, 0U);
    HUNDUN_CHECK(unbounded < 0.0);
    const auto bounded = integrator.integrate(
        read.pressure, read.velocity, read.gradient, read.viscosity);
    HUNDUN_CHECK(bounded.lowest_failing_rank == -1);
    HUNDUN_CHECK(finite(bounded));
    fill_fields(storage, ids, decomposition);
    read = views(storage, ids);
  }

  bool three_layer_contract_checked = false;
  bool three_layer_collectively_rejected = false;
  if (mpi.size() == 1 && !alternate) {
    std::uint64_t point_donor_count = 0U;
    std::uint64_t pressure_authority_donor_count = 0U;
    const auto owned = decomposition.owned_box();
    for (const auto &point : plan.local_points()) {
      const auto &point_donors =
          immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
              point.reconstruction);
      const auto pressure_authority =
          immersed::detail::boundary_authority_reconstruction(
              point.reconstruction);
      const auto &pressure_authority_donors =
          immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
              pressure_authority);
      HUNDUN_CHECK(!point_donors.empty());
      HUNDUN_CHECK(!pressure_authority_donors.empty());
      for (const auto donor : point_donors)
        HUNDUN_CHECK(inside(owned, logical_cell(donor)));
      for (const auto donor : pressure_authority_donors)
        HUNDUN_CHECK(inside(owned, logical_cell(donor)));
      point_donor_count += point_donors.size();
      pressure_authority_donor_count += pressure_authority_donors.size();
    }

    constexpr int three_layer_width = 3;
    runtime::FieldRegistry three_layer_registry;
    const auto three_layer_ids =
        declare_fields(three_layer_registry, three_layer_width);
    three_layer_registry.freeze();
    runtime::FieldStorage three_layer_storage(
        three_layer_registry, decomposition.local_extent());
    fill_fields(three_layer_storage, three_layer_ids, decomposition, 5.0, 1.0,
                three_layer_width);
    const auto three_layer_read =
        views(three_layer_storage, three_layer_ids);
    HUNDUN_CHECK(three_layer_read.pressure.ghost_width() == three_layer_width);
    HUNDUN_CHECK(three_layer_read.velocity.ghost_width() == three_layer_width);
    HUNDUN_CHECK(three_layer_read.gradient.ghost_width() == three_layer_width);
    HUNDUN_CHECK(three_layer_read.viscosity.ghost_width() == three_layer_width);
    const auto three_layer_sample = integrator.integrate(
        three_layer_read.pressure, three_layer_read.velocity,
        three_layer_read.gradient, three_layer_read.viscosity);
    const bool legacy_success = three_layer_sample.lowest_failing_rank == -1;
    three_layer_collectively_rejected =
        three_layer_sample.lowest_failing_rank == 0;
    HUNDUN_CHECK(legacy_success || three_layer_collectively_rejected);
    if (legacy_success) {
      HUNDUN_CHECK(finite(three_layer_sample));
      check_same(three_layer_sample, sample);
    } else {
      check_failed(three_layer_sample, 0);
    }
    std::fprintf(
        stdout,
        "TASK11_A5_WALL_FORCE point_donors=%llu "
        "pressure_authority_donors=%llu all_donors_owned=1 "
        "field_ghost_width=3 legacy_success=%d "
        "future_collective_reject_rank0=%d\n",
        static_cast<unsigned long long>(point_donor_count),
        static_cast<unsigned long long>(pressure_authority_donor_count),
        legacy_success ? 1 : 0,
        three_layer_collectively_rejected ? 1 : 0);
    std::fflush(stdout);
    three_layer_contract_checked = true;
  }

  std::vector<immersed::detail::WallPressureNormalGradient> wall_gradients;
  wall_gradients.reserve(domain.links().size());
  for (const auto &link : domain.links())
    wall_gradients.push_back(
        {link.id, dot(kPressureGradient, link.solid_to_fluid_normal)});
  const auto constrained_sample =
      immersed::detail::integrate_with_wall_pressure_authority(
          integrator, read.pressure, read.velocity, read.gradient,
          read.viscosity, wall_gradients);
  HUNDUN_CHECK(constrained_sample.lowest_failing_rank == -1);
  HUNDUN_CHECK(snapshot(storage, ids, decomposition.local_extent()) == before);
  check_near(constrained_sample.surface_traction.pressure_N,
             sample.surface_traction.pressure_N, 5.0e-12);

  {
    FixtureFile curved_file(mpi, false, {}, faceted_sphere());
    const auto curved_surface = immersed::ImmersedSurface::load_collective(
        curved_file.path(), 0.4, mpi, 0);
    const auto curved_query = immersed::SurfaceQuery::create(curved_surface);
    const auto curved_domain = immersed::ImmersedDomain::create(
        curved_surface, curved_query, config::ImmersedFluidSide::outside,
        topology, geometry, boundaries, mpi);
    const auto curved_plan = immersed::WallQuadraturePlan::create(
        curved_surface, curved_query, curved_domain, topology, geometry, mpi);
    const auto curved_integrator =
        immersed::WallForceIntegrator::create(curved_plan, mpi);
    if (mpi.size() == 1 && !alternate) {
      const auto curved_ghost_plan = immersed::GhostStencilPlan::create(
          curved_surface, curved_query, curved_domain, topology, geometry,
          decomposition, mpi);
      fill_polynomial_pressure(storage, ids.pressure, decomposition, true);
      const auto cubic_read = views(storage, ids);
      bool found_discriminating_row = false;
      bool independent_fit_mutation_detected = false;
      for (std::size_t left = 0U; left < curved_domain.links().size(); ++left) {
        const auto &left_link = curved_domain.links()[left];
        for (std::size_t right = left + 1U;
             right < curved_domain.links().size(); ++right) {
          const auto &right_link = curved_domain.links()[right];
          if (left_link.fluid_cell != right_link.fluid_cell)
            continue;
          const auto &left_independent =
              immersed::detail::QuadraticReconstructionWeights::
                  donor_global_ids(
                      curved_ghost_plan.reconstruction(left_link.id));
          const auto &right_independent =
              immersed::detail::QuadraticReconstructionWeights::
                  donor_global_ids(
                      curved_ghost_plan.reconstruction(right_link.id));
          if (left_independent == right_independent)
            continue;
          const immersed::QuadraticReconstruction *left_force = nullptr;
          const immersed::QuadraticReconstruction *right_force = nullptr;
          for (const auto &point : curved_plan.local_points()) {
            const auto link = immersed::detail::boundary_authority_link(
                point.reconstruction);
            if (link == left_link.id)
              left_force = &point.reconstruction;
            if (link == right_link.id)
              right_force = &point.reconstruction;
          }
          if (left_force == nullptr || right_force == nullptr)
            continue;
          const auto left_authority =
              immersed::detail::boundary_authority_reconstruction(*left_force);
          const auto right_authority =
              immersed::detail::boundary_authority_reconstruction(
                  *right_force);
          const auto &left_force_donors =
              immersed::detail::QuadraticReconstructionWeights::
                  donor_global_ids(left_authority);
          const auto &right_force_donors =
              immersed::detail::QuadraticReconstructionWeights::
                  donor_global_ids(right_authority);
          std::vector<mesh::GlobalCellId> expected_donors;
          for (const auto &row_link : curved_domain.links()) {
            if (row_link.fluid_cell != left_link.fluid_cell)
              continue;
            const auto &donors =
                immersed::detail::QuadraticReconstructionWeights::
                    donor_global_ids(
                        curved_ghost_plan.reconstruction(row_link.id));
            expected_donors.insert(expected_donors.end(), donors.begin(),
                                   donors.end());
          }
          const auto fluid =
              topology.find_local_cell(left_link.fluid_cell);
          HUNDUN_CHECK(fluid.has_value());
          const auto logical = topology.global_cell(*fluid);
          constexpr std::array<runtime::Int3, 6> offsets{
              runtime::Int3{-1, 0, 0}, runtime::Int3{1, 0, 0},
              runtime::Int3{0, -1, 0}, runtime::Int3{0, 1, 0},
              runtime::Int3{0, 0, -1}, runtime::Int3{0, 0, 1}};
          for (const auto offset : offsets) {
            const runtime::Int3 neighbour{logical.x + offset.x,
                                           logical.y + offset.y,
                                           logical.z + offset.z};
            if (!inside(neighbour))
              continue;
            const auto neighbour_id = topology.global_cell_id(neighbour);
            const auto local = topology.find_local_cell(neighbour_id);
            if (local.has_value() &&
                curved_domain.region(*local) == immersed::CellRegion::fluid)
              expected_donors.push_back(neighbour_id);
          }
          std::sort(expected_donors.begin(), expected_donors.end());
          expected_donors.erase(
              std::unique(expected_donors.begin(), expected_donors.end()),
              expected_donors.end());
          HUNDUN_CHECK(left_force_donors == expected_donors);
          HUNDUN_CHECK(left_force_donors == right_force_donors);
          const auto common_point =
              midpoint(left_link.wall_intercept_m, right_link.wall_intercept_m);
          const double left_shared_value =
              left_authority.value(common_point, cubic_read.pressure, 0U);
          const double right_shared_value =
              right_authority.value(common_point, cubic_read.pressure, 0U);
          const auto left_shared_gradient =
              left_authority.gradient(common_point, cubic_read.pressure, 0U);
          const auto right_shared_gradient =
              right_authority.gradient(common_point, cubic_read.pressure, 0U);
          const double basis_scale = std::max(
              {1.0, std::abs(left_shared_value), std::abs(right_shared_value),
               std::abs(left_shared_gradient.x),
               std::abs(left_shared_gradient.y),
               std::abs(left_shared_gradient.z),
               std::abs(right_shared_gradient.x),
               std::abs(right_shared_gradient.y),
               std::abs(right_shared_gradient.z)});
          const double basis_tolerance =
              8192.0 * std::numeric_limits<double>::epsilon() * basis_scale *
              (1.0 + std::max(left_authority.quality().condition_estimate,
                              right_authority.quality().condition_estimate));
          HUNDUN_CHECK(std::abs(left_shared_value - right_shared_value) <=
                       basis_tolerance);
          HUNDUN_CHECK(max_norm(subtract(left_shared_gradient,
                                         right_shared_gradient)) <=
                       basis_tolerance);
          const double left_independent_value =
              curved_ghost_plan.reconstruction(left_link.id)
                  .value(common_point, cubic_read.pressure, 0U);
          const double right_independent_value =
              curved_ghost_plan.reconstruction(right_link.id)
                  .value(common_point, cubic_read.pressure, 0U);
          independent_fit_mutation_detected =
              independent_fit_mutation_detected ||
              std::abs(left_independent_value - right_independent_value) >
                  basis_tolerance;
          found_discriminating_row = true;
        }
      }
      HUNDUN_CHECK(found_discriminating_row);
      HUNDUN_CHECK(independent_fit_mutation_detected);
      fill_fields(storage, ids, decomposition);
      read = views(storage, ids);
    }
    std::vector<immersed::detail::WallPressureNormalGradient> curved_gradients;
    curved_gradients.reserve(curved_domain.links().size());
    for (const auto &link : curved_domain.links())
      curved_gradients.push_back(
          {link.id, dot(kPressureGradient, link.solid_to_fluid_normal)});
    const auto curved_unconstrained = curved_integrator.integrate(
        read.pressure, read.velocity, read.gradient, read.viscosity);
    const auto curved_constrained =
        immersed::detail::integrate_with_wall_pressure_authority(
            curved_integrator, read.pressure, read.velocity, read.gradient,
            read.viscosity, curved_gradients);
    HUNDUN_CHECK(curved_unconstrained.lowest_failing_rank == -1);
    HUNDUN_CHECK(curved_constrained.lowest_failing_rank == -1);
    check_near(curved_constrained.surface_traction.pressure_N,
               curved_unconstrained.surface_traction.pressure_N, 5.0e-12);
    check_near(curved_constrained.moment_about_global_origin.pressure_N_m,
               curved_unconstrained.moment_about_global_origin.pressure_N_m,
               5.0e-12);

    fill_polynomial_pressure(storage, ids.pressure, decomposition, false);
    read = views(storage, ids);
    std::vector<immersed::detail::WallPressureNormalGradient>
        curved_quadratic_gradients;
    curved_quadratic_gradients.reserve(curved_domain.links().size());
    for (const auto &link : curved_domain.links())
      curved_quadratic_gradients.push_back(
          {link.id, dot(quadratic_pressure_gradient(link.wall_intercept_m),
                        link.solid_to_fluid_normal)});
    const auto curved_quadratic =
        immersed::detail::integrate_with_wall_pressure_authority(
            curved_integrator, read.pressure, read.velocity, read.gradient,
            read.viscosity, curved_quadratic_gradients);
    HUNDUN_CHECK(curved_quadratic.lowest_failing_rank == -1);
    check_near(curved_quadratic.surface_traction.pressure_N,
               sampled_pressure_force(curved_plan, mpi, quadratic_pressure),
               5.0e-12);
    fill_fields(storage, ids, decomposition);
    read = views(storage, ids);
  }

  auto mutated_wall_gradients = wall_gradients;
  for (auto &gradient : mutated_wall_gradients)
    gradient.value += 0.25 * (1.0 + static_cast<double>(gradient.link % 7U));
  const auto mutated_constrained_sample =
      immersed::detail::integrate_with_wall_pressure_authority(
          integrator, read.pressure, read.velocity, read.gradient,
          read.viscosity, mutated_wall_gradients);
  HUNDUN_CHECK(mutated_constrained_sample.lowest_failing_rank == -1);
  HUNDUN_CHECK(
      max_norm(
          add(mutated_constrained_sample.surface_traction.pressure_N,
              scale(-1.0, constrained_sample.surface_traction.pressure_N))) >
          1.0e-8 ||
      max_norm(add(
          mutated_constrained_sample.moment_about_global_origin.pressure_N_m,
          scale(-1.0,
                constrained_sample.moment_about_global_origin.pressure_N_m))) >
          1.0e-8);

  int local_gradient_rank = wall_gradients.empty() ? mpi.size() : mpi.rank();
  int lowest_gradient_rank = mpi.size();
  HUNDUN_CHECK(MPI_Allreduce(&local_gradient_rank, &lowest_gradient_rank, 1,
                             MPI_INT, MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(lowest_gradient_rank < mpi.size());
  check_failed(immersed::detail::integrate_with_wall_pressure_authority(
                   integrator, read.pressure, read.velocity, read.gradient,
                   read.viscosity, {}),
               lowest_gradient_rank);
  auto duplicate_wall_gradients = wall_gradients;
  if (mpi.rank() == lowest_gradient_rank)
    duplicate_wall_gradients.push_back(duplicate_wall_gradients.front());
  check_failed(immersed::detail::integrate_with_wall_pressure_authority(
                   integrator, read.pressure, read.velocity, read.gradient,
                   read.viscosity, duplicate_wall_gradients),
               lowest_gradient_rank);
  if (mpi.size() > 1) {
    constexpr immersed::ImmersedLinkId moved_link = 0U;
    const auto owns_moved_link = std::find_if(
        wall_gradients.begin(), wall_gradients.end(),
        [](const auto &gradient) { return gradient.link == moved_link; });
    int local_owner =
        owns_moved_link == wall_gradients.end() ? mpi.size() : mpi.rank();
    int owner = mpi.size();
    HUNDUN_CHECK(MPI_Allreduce(&local_owner, &owner, 1, MPI_INT, MPI_MIN,
                               mpi.comm()) == MPI_SUCCESS);
    HUNDUN_CHECK(owner < mpi.size());
    const int wrong_provider = (owner + 1) % mpi.size();
    auto wrong_provider_gradients = wall_gradients;
    if (mpi.rank() == owner) {
      wrong_provider_gradients.erase(std::find_if(
          wrong_provider_gradients.begin(), wrong_provider_gradients.end(),
          [](const auto &gradient) { return gradient.link == moved_link; }));
    }
    if (mpi.rank() == wrong_provider)
      wrong_provider_gradients.push_back({moved_link, 0.0});
    check_failed(immersed::detail::integrate_with_wall_pressure_authority(
                     integrator, read.pressure, read.velocity, read.gradient,
                     read.viscosity, wrong_provider_gradients),
                 std::min(owner, wrong_provider));
  }
  auto extra_wall_gradients = wall_gradients;
  if (mpi.rank() == lowest_gradient_rank)
    extra_wall_gradients.push_back(
        {std::numeric_limits<immersed::ImmersedLinkId>::max(), 0.0});
  check_failed(immersed::detail::integrate_with_wall_pressure_authority(
                   integrator, read.pressure, read.velocity, read.gradient,
                   read.viscosity, extra_wall_gradients),
               lowest_gradient_rank);
  auto nonfinite_wall_gradients = wall_gradients;
  if (mpi.rank() == lowest_gradient_rank)
    nonfinite_wall_gradients.front().value =
        std::numeric_limits<double>::quiet_NaN();
  check_failed(immersed::detail::integrate_with_wall_pressure_authority(
                   integrator, read.pressure, read.velocity, read.gradient,
                   read.viscosity, nonfinite_wall_gradients),
               lowest_gradient_rank);

  const auto expected = oracle(plan, mpi, read.velocity);
  HUNDUN_CHECK(sample.quadrature_point_count == 3U * surface.triangle_count());
  HUNDUN_CHECK(sample.quadrature_point_count == expected.count);
  const double tolerance =
      5.0e-12 *
      std::max(
          {1.0,
           std::hypot(expected.force.pressure_N.x, expected.force.pressure_N.y,
                      expected.force.pressure_N.z),
           std::hypot(expected.force.viscous_N.x, expected.force.viscous_N.y,
                      expected.force.viscous_N.z)});
  check_near(sample.surface_traction.pressure_N, expected.force.pressure_N,
             tolerance);
  check_near(sample.surface_traction.viscous_N, expected.force.viscous_N,
             tolerance);
  HUNDUN_CHECK(max_norm(add(expected.force.viscous_N,
                            scale(-1.0, expected.unprojected_viscous_N))) >
               100.0 * tolerance);
  check_near(sample.surface_traction.total_N, expected.force.total_N,
             tolerance);
  check_near(sample.moment_about_global_origin.pressure_N_m,
             expected.moment.pressure_N_m, tolerance);
  check_near(sample.moment_about_global_origin.viscous_N_m,
             expected.moment.viscous_N_m, tolerance);
  check_near(sample.moment_about_global_origin.total_N_m,
             expected.moment.total_N_m, tolerance);
  check_near(sample.area_vector_closure_m2, expected.closure, tolerance);
  check_near(sample.surface_traction.total_N,
             add(sample.surface_traction.pressure_N,
                 sample.surface_traction.viscous_N),
             tolerance);
  check_near(sample.moment_about_global_origin.total_N_m,
             add(sample.moment_about_global_origin.pressure_N_m,
                 sample.moment_about_global_origin.viscous_N_m),
             tolerance);

  constexpr double volume = 0.4 * 0.4 * 0.4;
  check_near(sample.surface_traction.pressure_N,
             scale(-volume, kPressureGradient), tolerance);
  HUNDUN_CHECK(max_norm(sample.surface_traction.pressure_N) > 1.0e-6);
  HUNDUN_CHECK(max_norm(sample.surface_traction.viscous_N) > 1.0e-6);
  HUNDUN_CHECK(max_norm(sample.moment_about_global_origin.pressure_N_m) >
               1.0e-6);

  fill_polynomial_pressure(storage, ids.pressure, decomposition, false);
  read = views(storage, ids);
  const auto quadratic_sample = integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  HUNDUN_CHECK(quadratic_sample.lowest_failing_rank == -1);
  check_near(quadratic_sample.surface_traction.pressure_N,
             sampled_pressure_force(plan, mpi, quadratic_pressure), tolerance);

  fill_polynomial_pressure(storage, ids.pressure, decomposition, true);
  read = views(storage, ids);
  const auto nonlinear_before =
      snapshot(storage, ids, decomposition.local_extent());
  const auto constrained_pressure_force =
      reconstructed_pressure_force(plan, mpi, read.pressure, true);
  const auto unconstrained_pressure_force =
      reconstructed_pressure_force(plan, mpi, read.pressure, false);
  check_near(constrained_pressure_force, unconstrained_pressure_force,
             tolerance);
  const auto nonlinear_sample = integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  HUNDUN_CHECK(nonlinear_sample.lowest_failing_rank == -1);
  HUNDUN_CHECK(snapshot(storage, ids, decomposition.local_extent()) ==
               nonlinear_before);
  check_near(nonlinear_sample.surface_traction.pressure_N,
             constrained_pressure_force, tolerance);
  fill_fields(storage, ids, decomposition);
  read = views(storage, ids);

  for (int velocity_component = 0; velocity_component < 3; ++velocity_component)
    mutate_component(storage, ids.velocity, decomposition, velocity_component,
                     [](double) { return 0.0; });
  read = views(storage, ids);
  const auto zero_velocity = integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  HUNDUN_CHECK(zero_velocity.lowest_failing_rank == -1);
  check_near(zero_velocity.surface_traction.viscous_N, {}, tolerance);
  fill_fields(storage, ids, decomposition);
  read = views(storage, ids);

  const auto counters_before = mpi.fp64_reduction_counters();
  const auto repeated = integrator.integrate(read.pressure, read.velocity,
                                             read.gradient, read.viscosity);
  const auto counters_after = mpi.fp64_reduction_counters();
  check_same(sample, repeated);
  HUNDUN_CHECK(counters_after.collective_calls ==
               counters_before.collective_calls);
  HUNDUN_CHECK(counters_after.reduced_scalars ==
               counters_before.reduced_scalars);
  HUNDUN_CHECK(counters_after.logical_payload_bytes ==
               counters_before.logical_payload_bytes);
  HUNDUN_CHECK(plan.fingerprint() == plan_fingerprint);

  FixtureFile permuted_file(mpi, true);
  const auto permuted_surface = immersed::ImmersedSurface::load_collective(
      permuted_file.path(), 0.4, mpi, 0);
  const auto permuted_query = immersed::SurfaceQuery::create(permuted_surface);
  const auto permuted_domain = immersed::ImmersedDomain::create(
      permuted_surface, permuted_query, config::ImmersedFluidSide::outside,
      topology, geometry, boundaries, mpi);
  const auto permuted_plan = immersed::WallQuadraturePlan::create(
      permuted_surface, permuted_query, permuted_domain, topology, geometry,
      mpi);
  const auto permuted_integrator =
      immersed::WallForceIntegrator::create(permuted_plan, mpi);
  const auto permuted_sample = permuted_integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  check_same(sample, permuted_sample);

  FixtureFile translated_file(mpi, false, {0.76, 0.75, 0.75});
  const auto translated_surface = immersed::ImmersedSurface::load_collective(
      translated_file.path(), 0.4, mpi, 0);
  const auto translated_query =
      immersed::SurfaceQuery::create(translated_surface);
  const auto translated_domain = immersed::ImmersedDomain::create(
      translated_surface, translated_query, config::ImmersedFluidSide::outside,
      topology, geometry, boundaries, mpi);
  const auto translated_plan = immersed::WallQuadraturePlan::create(
      translated_surface, translated_query, translated_domain, topology,
      geometry, mpi);
  const auto translated_integrator =
      immersed::WallForceIntegrator::create(translated_plan, mpi);
  const auto translated_sample = translated_integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  HUNDUN_CHECK(
      bits(translated_sample.moment_about_global_origin.pressure_N_m.y) !=
          bits(sample.moment_about_global_origin.pressure_N_m.y) ||
      bits(translated_sample.moment_about_global_origin.pressure_N_m.z) !=
          bits(sample.moment_about_global_origin.pressure_N_m.z));

  fill_fields(storage, ids, decomposition, 42.0, 0.0);
  read = views(storage, ids);
  const auto constant_pressure = integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  check_near(constant_pressure.surface_traction.pressure_N,
             scale(-42.0, constant_pressure.area_vector_closure_m2), tolerance);
  check_near(constant_pressure.surface_traction.pressure_N, {}, tolerance);

  fill_fields(storage, ids, decomposition, 42.0);
  read = views(storage, ids);
  const auto shifted = integrator.integrate(read.pressure, read.velocity,
                                            read.gradient, read.viscosity);
  check_near(shifted.surface_traction.pressure_N,
             sample.surface_traction.pressure_N, tolerance);
  check_near(shifted.moment_about_global_origin.pressure_N_m,
             sample.moment_about_global_origin.pressure_N_m, tolerance);
  fill_fields(storage, ids, decomposition);

  mutate_component(storage, ids.viscosity, decomposition, 0,
                   [](double value) { return 1.25 * value; });
  read = views(storage, ids);
  const auto viscosity_mutated = integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  HUNDUN_CHECK(bits(viscosity_mutated.surface_traction.viscous_N.x) !=
                   bits(sample.surface_traction.viscous_N.x) ||
               bits(viscosity_mutated.surface_traction.viscous_N.y) !=
                   bits(sample.surface_traction.viscous_N.y) ||
               bits(viscosity_mutated.surface_traction.viscous_N.z) !=
                   bits(sample.surface_traction.viscous_N.z));
  fill_fields(storage, ids, decomposition);

  mutate_component(storage, ids.gradient, decomposition, 1,
                   [](double value) { return 1.25 * value; });
  read = views(storage, ids);
  const auto gradient_mutated = integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  check_same(sample, gradient_mutated);
  fill_fields(storage, ids, decomposition);

  fill_fields(storage, ids, decomposition, 5.0, 1.25);
  read = views(storage, ids);
  const auto mutated = integrator.integrate(read.pressure, read.velocity,
                                            read.gradient, read.viscosity);
  HUNDUN_CHECK(bits(mutated.surface_traction.pressure_N.x) !=
                   bits(sample.surface_traction.pressure_N.x) ||
               bits(mutated.surface_traction.pressure_N.y) !=
                   bits(sample.surface_traction.pressure_N.y) ||
               bits(mutated.surface_traction.pressure_N.z) !=
                   bits(sample.surface_traction.pressure_N.z));
  fill_fields(storage, ids, decomposition);

  read = views(storage, ids);
  const auto bad_components = integrator.integrate(
      read.bad_pressure, read.velocity, read.gradient, read.viscosity);
  check_failed(bad_components, 0);

  runtime::FieldStorage mismatched(
      registry, runtime::Int3{decomposition.local_extent().x + 1,
                              decomposition.local_extent().y,
                              decomposition.local_extent().z});
  const auto mismatched_read = views(mismatched, ids);
  const auto mismatched_layout = integrator.integrate(
      mismatched_read.pressure, read.velocity, read.gradient, read.viscosity);
  check_failed(mismatched_layout, 0);

  int local_owner = plan.local_points().empty() ? -1 : mpi.rank();
  int failing_rank = -1;
  HUNDUN_CHECK(MPI_Allreduce(&local_owner, &failing_rank, 1, MPI_INT, MPI_MAX,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(failing_rank >= 0);

  if (mpi.rank() == failing_rank) {
    auto write = storage.view<double>(ids.pressure);
    const auto local = decomposition.local_extent();
    for (int k = -kGhost; k < local.z + kGhost; ++k) {
      for (int j = -kGhost; j < local.y + kGhost; ++j) {
        for (int i = -kGhost; i < local.x + kGhost; ++i) {
          write(i, j, k, 0) = std::numeric_limits<double>::quiet_NaN();
        }
      }
    }
  }
  read = views(storage, ids);
  const auto nonfinite = integrator.integrate(read.pressure, read.velocity,
                                              read.gradient, read.viscosity);
  check_failed(nonfinite, failing_rank);
  fill_fields(storage, ids, decomposition);

  if (mpi.rank() == failing_rank) {
    mutate_component(storage, ids.velocity, decomposition, 2, [](double) {
      return std::numeric_limits<double>::quiet_NaN();
    });
  }
  read = views(storage, ids);
  const auto before_nonfinite_velocity =
      snapshot(storage, ids, decomposition.local_extent());
  const auto nonfinite_velocity = integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  check_failed(nonfinite_velocity, failing_rank);
  HUNDUN_CHECK(snapshot(storage, ids, decomposition.local_extent()) ==
               before_nonfinite_velocity);
  fill_fields(storage, ids, decomposition);

  if (mpi.rank() == failing_rank) {
    auto write = storage.view<double>(ids.gradient);
    const auto local = decomposition.local_extent();
    for (int k = -kGhost; k < local.z + kGhost; ++k) {
      for (int j = -kGhost; j < local.y + kGhost; ++j) {
        for (int i = -kGhost; i < local.x + kGhost; ++i) {
          write(i, j, k, 5) = std::numeric_limits<double>::infinity();
        }
      }
    }
  }
  read = views(storage, ids);
  const auto nonfinite_gradient = integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  check_failed(nonfinite_gradient, failing_rank);
  fill_fields(storage, ids, decomposition);

  if (mpi.rank() == failing_rank) {
    mutate_component(storage, ids.viscosity, decomposition, 0, [](double) {
      return std::numeric_limits<double>::infinity();
    });
  }
  read = views(storage, ids);
  const auto before_nonfinite_viscosity =
      snapshot(storage, ids, decomposition.local_extent());
  const auto nonfinite_viscosity = integrator.integrate(
      read.pressure, read.velocity, read.gradient, read.viscosity);
  check_failed(nonfinite_viscosity, failing_rank);
  HUNDUN_CHECK(snapshot(storage, ids, decomposition.local_extent()) ==
               before_nonfinite_viscosity);
  fill_fields(storage, ids, decomposition);

  read = views(storage, ids);
  if (mpi.rank() == failing_rank) {
    storage.begin_rebuild();
  }
  const auto stale = integrator.integrate(read.pressure, read.velocity,
                                          read.gradient, read.viscosity);
  check_failed(stale, failing_rank);
  fill_fields(storage, ids, decomposition);
  read = views(storage, ids);

  if (mpi.rank() == failing_rank) {
    auto write = storage.view<double>(ids.viscosity);
    const auto local = decomposition.local_extent();
    for (int k = -kGhost; k < local.z + kGhost; ++k) {
      for (int j = -kGhost; j < local.y + kGhost; ++j) {
        for (int i = -kGhost; i < local.x + kGhost; ++i) {
          write(i, j, k, 0) = -1.0;
        }
      }
    }
  }
  read = views(storage, ids);
  const auto negative_mu = integrator.integrate(read.pressure, read.velocity,
                                                read.gradient, read.viscosity);
  check_failed(negative_mu, failing_rank);
  HUNDUN_CHECK(plan.fingerprint() == plan_fingerprint);

  if (mpi.size() == 1 && !alternate) {
    std::fprintf(stdout,
                 "TASK11_A5_WALL_FORCE expected="
                 "collective_width3_reject_before_donor_read_rank0\n");
    std::fflush(stdout);
    HUNDUN_CHECK(three_layer_contract_checked);
    HUNDUN_CHECK(three_layer_collectively_rejected);
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    const bool alternate = argc == 2 && std::string(argv[1]) == "alternate";
    HUNDUN_CHECK(argc == 1 || alternate);
    HUNDUN_CHECK(!alternate || mpi.size() == 2 || mpi.size() == 4);
    run(mpi, alternate);
  });
}
