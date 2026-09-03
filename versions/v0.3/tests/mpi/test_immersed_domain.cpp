// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_domain.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using hundun::boundary::BoundaryRegistry;
using hundun::config::BoundaryType;
using hundun::config::DensityModel;
using hundun::config::FlowBoundaryConfig;
using hundun::config::FlowCaseConfig;
using hundun::config::ImmersedFluidSide;
using hundun::config::InletScalarValue;
using hundun::config::InletThermalAuthority;
using hundun::config::MeshMapping;
using hundun::config::PatchName;
using hundun::config::SimulationType;
using hundun::config::TimeMode;
using hundun::immersed::CellRegion;
using hundun::immersed::ImmersedDomain;
using hundun::immersed::ImmersedLink;
using hundun::immersed::ImmersedSurface;
using hundun::immersed::SurfaceQuery;
using hundun::mesh::AnalyticWarpedBoxMapping;
using hundun::mesh::GlobalCellId;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::mesh::UniformBoxMapping;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::Real3;
using hundun::runtime::StructuredDecomposition;

constexpr Int3 kExtent{12, 12, 12};
constexpr Real3 kOrigin{0.0, 0.0, 0.0};
constexpr Real3 kLength{1.0, 1.0, 1.0};
constexpr Real3 kWarp{0.02, -0.015, 0.01};
constexpr std::array<PatchName, 6> kPatchNames{
    PatchName::x_min, PatchName::x_max, PatchName::y_min,
    PatchName::y_max, PatchName::z_min, PatchName::z_max};

class SharedDirectory final {
public:
  explicit SharedDirectory(const MpiContext &mpi)
      : mpi_(&mpi), rank_(mpi.rank()) {
    std::string text;
    if (rank_ == 0) {
      const auto stamp =
          std::chrono::steady_clock::now().time_since_epoch().count();
      path_ = std::filesystem::temp_directory_path() /
              ("hundun-stage3-domain-" + std::to_string(stamp));
      std::filesystem::create_directories(path_);
      text = path_.generic_string();
    }
    std::uint64_t length = text.size();
    HUNDUN_CHECK(MPI_Bcast(&length, 1, MPI_UINT64_T, 0, mpi_->comm()) ==
                 MPI_SUCCESS);
    HUNDUN_CHECK(length <=
                 static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
    text.resize(static_cast<std::size_t>(length));
    HUNDUN_CHECK(MPI_Bcast(text.data(), static_cast<int>(text.size()), MPI_BYTE,
                           0, mpi_->comm()) == MPI_SUCCESS);
    path_ = text;
    HUNDUN_CHECK(MPI_Barrier(mpi_->comm()) == MPI_SUCCESS);
  }

  ~SharedDirectory() {
    MPI_Barrier(mpi_->comm());
    if (rank_ == 0) {
      std::error_code error;
      std::filesystem::remove_all(path_, error);
    }
  }

  const std::filesystem::path &path() const noexcept { return path_; }

private:
  const MpiContext *mpi_{};
  int rank_{};
  std::filesystem::path path_;
};

std::uint64_t bits(double value) {
  std::uint64_t result = 0U;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

Int3 process_grid(int ranks, bool alternate) {
  switch (ranks) {
  case 1:
    return {1, 1, 1};
  case 2:
    return alternate ? Int3{1, 2, 1} : Int3{2, 1, 1};
  case 4:
    return alternate ? Int3{4, 1, 1} : Int3{2, 2, 1};
  default:
    throw hundun::runtime::Error("unsupported Task 3 test rank count");
  }
}

FlowBoundaryConfig wall_boundary(std::size_t stable_id) {
  FlowBoundaryConfig boundary{};
  boundary.patch = kPatchNames.at(stable_id);
  boundary.type = BoundaryType::no_slip_wall;
  return boundary;
}

FlowCaseConfig closed_config(int ranks, Int3 grid, Int3 extent = kExtent) {
  FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "stage3_task3_domain";
  config.simulation_type = SimulationType::variable_density_flow;
  config.density_model = DensityModel::constant;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = grid;
  config.mesh.cells = extent;
  config.mesh.origin_m = kOrigin;
  config.mesh.length_m = kLength;
  config.mesh.mapping = MeshMapping::uniform_box;
  config.time.mode = TimeMode::fixed;
  config.time.steps = 1;
  config.time.initial_dt_s = 0.01;
  config.time.min_dt_s = 0.01;
  config.time.max_dt_s = 0.01;
  config.time.cfl_target = 0.5;
  config.time.diffusion_number_target = 0.25;
  config.time.growth_factor = 1.25;
  config.time.retry_factor = 0.5;
  config.time.max_retries = 8;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.dynamic_viscosity_pa_s = 1.0e-3;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.boundaries = {wall_boundary(0U), wall_boundary(1U), wall_boundary(2U),
                       wall_boundary(3U), wall_boundary(4U), wall_boundary(5U)};
  return config;
}

FlowCaseConfig open_config(int ranks, Int3 grid, Int3 extent) {
  FlowCaseConfig config = closed_config(ranks, grid, extent);
  config.boundaries[0].type = BoundaryType::velocity_inlet;
  config.boundaries[0].velocity_m_per_s = Real3{1.0, 0.0, 0.0};
  config.boundaries[0].thermal_authority = InletThermalAuthority::enthalpy;
  config.boundaries[0].enthalpy_J_per_kg = 1.0;
  config.boundaries[0].scalar_values = std::vector<InletScalarValue>{};
  config.boundaries[1].type = BoundaryType::pressure_outlet;
  config.boundaries[1].pressure_perturbation_pa = 0.0;
  return config;
}

std::vector<hundun::test::StlFixtureTriangle> cube_between(double minimum,
                                                           double maximum) {
  auto triangles = hundun::test::outward_cube();
  const double scale = maximum - minimum;
  for (auto &triangle : triangles) {
    for (Real3 &vertex : triangle.vertices) {
      vertex.x = minimum + scale * vertex.x;
      vertex.y = minimum + scale * vertex.y;
      vertex.z = minimum + scale * vertex.z;
    }
  }
  return triangles;
}

std::vector<hundun::test::StlFixtureTriangle> box_between(Real3 minimum,
                                                          Real3 maximum) {
  auto triangles = hundun::test::outward_cube();
  for (auto &triangle_value : triangles) {
    for (Real3 &vertex : triangle_value.vertices) {
      vertex.x = minimum.x + (maximum.x - minimum.x) * vertex.x;
      vertex.y = minimum.y + (maximum.y - minimum.y) * vertex.y;
      vertex.z = minimum.z + (maximum.z - minimum.z) * vertex.z;
    }
  }
  return triangles;
}

hundun::test::StlFixtureTriangle triangle(Real3 a, Real3 b, Real3 c) {
  const Real3 first{b.x - a.x, b.y - a.y, b.z - a.z};
  const Real3 second{c.x - a.x, c.y - a.y, c.z - a.z};
  Real3 normal{first.y * second.z - first.z * second.y,
               first.z * second.x - first.x * second.z,
               first.x * second.y - first.y * second.x};
  const double magnitude = std::hypot(normal.x, normal.y, normal.z);
  HUNDUN_CHECK(std::isfinite(magnitude) && magnitude > 0.0);
  normal = {normal.x / magnitude, normal.y / magnitude, normal.z / magnitude};
  return {normal, {a, b, c}};
}

std::vector<hundun::test::StlFixtureTriangle> connected_concave_extrusion() {
  const Real3 v4{0.25, 0.25, 0.75};
  const Real3 v5{0.75, 0.25, 0.75};
  const Real3 v6{0.75, 0.75, 0.75};
  const Real3 v7{0.25, 0.75, 0.75};
  const Real3 dent{0.35, 0.50, 0.45};
  auto cube = cube_between(0.25, 0.75);
  std::vector<hundun::test::StlFixtureTriangle> triangles;
  triangles.reserve(cube.size() + 2U);
  for (std::size_t index = 0U; index < cube.size(); ++index) {
    if (index != 2U && index != 3U) {
      triangles.push_back(cube[index]);
    }
  }
  triangles.push_back(triangle(v4, v5, dent));
  triangles.push_back(triangle(v5, v6, dent));
  triangles.push_back(triangle(v6, v7, dent));
  triangles.push_back(triangle(v7, v4, dent));
  return triangles;
}

std::filesystem::path
write_surface(const SharedDirectory &directory, const MpiContext &mpi,
              std::string name,
              const std::vector<hundun::test::StlFixtureTriangle> &triangles) {
  const auto path = directory.path() / (name + ".stl");
  if (mpi.rank() == 0) {
    hundun::test::write_text(
        path, hundun::test::ascii_stl(triangles, std::move(name)));
  }
  HUNDUN_CHECK(MPI_Barrier(mpi.comm()) == MPI_SUCCESS);
  return path;
}

std::filesystem::path write_cube(const SharedDirectory &directory,
                                 const MpiContext &mpi, double minimum = 0.25,
                                 double maximum = 0.75) {
  return write_surface(directory, mpi, "task3-domain-cube",
                       cube_between(minimum, maximum));
}

template <class Operation>
void expect_error(Operation operation, const std::string &marker) {
  bool rejected = false;
  try {
    operation();
  } catch (const hundun::runtime::Error &error) {
    rejected = std::string(error.what()).find(marker) != std::string::npos;
  }
  HUNDUN_CHECK(rejected);
}

template <class Operation>
void expect_collective_error(Operation operation, const MpiContext &mpi,
                             const std::string &marker,
                             int expected_failing_rank) {
  bool threw = false;
  std::string message;
  try {
    operation();
  } catch (const hundun::runtime::Error &error) {
    threw = true;
    message = error.what();
  }
  int every_threw = threw ? 1 : 0;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &every_threw, 1, MPI_INT, MPI_MIN,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(every_threw == 1);
  std::string reference = mpi.rank() == 0 ? message : std::string{};
  std::uint64_t length = reference.size();
  HUNDUN_CHECK(MPI_Bcast(&length, 1, MPI_UINT64_T, 0, mpi.comm()) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(length <=
               static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
  reference.resize(static_cast<std::size_t>(length));
  HUNDUN_CHECK(MPI_Bcast(reference.data(), static_cast<int>(reference.size()),
                         MPI_BYTE, 0, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(message == reference);
  if (message.find(marker) == std::string::npos) {
    throw std::runtime_error("expected collective marker '" + marker +
                             "', received '" + message + "'");
  }
  HUNDUN_CHECK(message.find("lowest failing rank " +
                            std::to_string(expected_failing_rank)) !=
               std::string::npos);
}

std::vector<std::uint64_t> gather_u64(const std::vector<std::uint64_t> &local,
                                      const MpiContext &mpi) {
  HUNDUN_CHECK(local.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()));
  const int local_count = static_cast<int>(local.size());
  std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
  HUNDUN_CHECK(MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1,
                             MPI_INT, mpi.comm()) == MPI_SUCCESS);
  std::vector<int> displacements(counts.size());
  int total = 0;
  for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
    HUNDUN_CHECK(counts[rank] >= 0 &&
                 total <= std::numeric_limits<int>::max() - counts[rank]);
    displacements[rank] = total;
    total += counts[rank];
  }
  std::vector<std::uint64_t> gathered(static_cast<std::size_t>(total));
  HUNDUN_CHECK(MPI_Allgatherv(local.data(), local_count, MPI_UINT64_T,
                              gathered.data(), counts.data(),
                              displacements.data(), MPI_UINT64_T,
                              mpi.comm()) == MPI_SUCCESS);
  return gathered;
}

double value_from_bits(std::uint64_t value) {
  double result = 0.0;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::vector<Real3> gather_owned_cell_centers(const MeshTopology &topology,
                                             const MeshGeometry &geometry,
                                             const MpiContext &mpi) {
  constexpr std::size_t record_width = 4U;
  std::vector<std::uint64_t> local;
  local.reserve(topology.owned_cell_count() * record_width);
  for (std::size_t cell = 0U; cell < topology.owned_cell_count(); ++cell) {
    const Real3 center = geometry.cell_center_m(cell);
    local.insert(local.end(), {topology.global_cell_id(cell), bits(center.x),
                               bits(center.y), bits(center.z)});
  }
  const auto gathered = gather_u64(local, mpi);
  HUNDUN_CHECK(gathered.size() % record_width == 0U);
  std::vector<std::array<std::uint64_t, record_width>> records(gathered.size() /
                                                               record_width);
  for (std::size_t record = 0U; record < records.size(); ++record) {
    std::copy_n(gathered.begin() +
                    static_cast<std::ptrdiff_t>(record * record_width),
                record_width, records[record].begin());
  }
  std::sort(records.begin(), records.end(),
            [](const auto &first, const auto &second) {
              return first[0] < second[0];
            });
  HUNDUN_CHECK(records.size() == topology.global_cell_count());
  std::vector<Real3> result;
  result.reserve(records.size());
  for (std::size_t record = 0U; record < records.size(); ++record) {
    HUNDUN_CHECK(records[record][0] == record);
    result.push_back({value_from_bits(records[record][1]),
                      value_from_bits(records[record][2]),
                      value_from_bits(records[record][3])});
  }
  return result;
}

std::vector<double> gather_owned_cell_sizes(const MeshTopology &topology,
                                            const MeshGeometry &geometry,
                                            const MpiContext &mpi) {
  constexpr std::size_t record_width = 2U;
  std::vector<std::uint64_t> local;
  local.reserve(topology.owned_cell_count() * record_width);
  for (std::size_t cell = 0U; cell < topology.owned_cell_count(); ++cell) {
    local.insert(local.end(), {topology.global_cell_id(cell),
                               bits(std::cbrt(geometry.cell_volume_m3(cell)))});
  }
  const auto gathered = gather_u64(local, mpi);
  HUNDUN_CHECK(gathered.size() % record_width == 0U);
  std::vector<std::array<std::uint64_t, record_width>> records(gathered.size() /
                                                               record_width);
  for (std::size_t record = 0U; record < records.size(); ++record) {
    std::copy_n(gathered.begin() +
                    static_cast<std::ptrdiff_t>(record * record_width),
                record_width, records[record].begin());
  }
  std::sort(records.begin(), records.end(),
            [](const auto &first, const auto &second) {
              return first[0] < second[0];
            });
  HUNDUN_CHECK(records.size() == topology.global_cell_count());
  std::vector<double> result;
  result.reserve(records.size());
  for (std::size_t record = 0U; record < records.size(); ++record) {
    HUNDUN_CHECK(records[record][0] == record);
    const double h = value_from_bits(records[record][1]);
    HUNDUN_CHECK(std::isfinite(h) && h > 0.0);
    result.push_back(h);
  }
  return result;
}

using PackedLink = std::array<std::uint64_t, 11>;

PackedLink pack(const ImmersedLink &link) {
  return {link.id,
          link.fluid_cell,
          link.solid_cell,
          link.triangle,
          bits(link.wall_intercept_m.x),
          bits(link.wall_intercept_m.y),
          bits(link.wall_intercept_m.z),
          bits(link.solid_to_fluid_normal.x),
          bits(link.solid_to_fluid_normal.y),
          bits(link.solid_to_fluid_normal.z),
          bits(link.fluid_to_wall_fraction)};
}

struct LinkKey final {
  GlobalCellId first{};
  GlobalCellId second{};
  std::uint64_t triangle{};
  std::array<std::uint64_t, 3> intercept{};
};

bool less_key(const LinkKey &first, const LinkKey &second) {
  return std::tie(first.first, first.second, first.triangle, first.intercept) <
         std::tie(second.first, second.second, second.triangle,
                  second.intercept);
}

struct OracleLink final {
  LinkKey key{};
  PackedLink packed{};
};

std::vector<OracleLink> oracle_links(const MeshTopology &topology,
                                     const ImmersedSurface &surface,
                                     const SurfaceQuery &query,
                                     ImmersedFluidSide fluid_side,
                                     const std::vector<Real3> &global_centers) {
  std::vector<OracleLink> result;
  const auto append = [&](Int3 lower, Int3 upper) {
    const GlobalCellId lower_id = topology.global_cell_id(lower);
    const GlobalCellId upper_id = topology.global_cell_id(upper);
    const Real3 lower_center = global_centers.at(lower_id);
    const Real3 upper_center = global_centers.at(upper_id);
    const CellRegion lower_region = query.classify(lower_center, fluid_side);
    const CellRegion upper_region = query.classify(upper_center, fluid_side);
    if (lower_region == upper_region) {
      return;
    }
    const auto hits = query.segment_intersections(lower_center, upper_center);
    HUNDUN_CHECK(hits.size() == 1U);
    const auto &hit = hits.front();
    const bool lower_fluid = lower_region == CellRegion::fluid;
    const auto &triangle = surface.triangle(hit.triangle);
    const double sign = fluid_side == ImmersedFluidSide::outside ? 1.0 : -1.0;
    PackedLink packed{
        0U,
        lower_fluid ? lower_id : upper_id,
        lower_fluid ? upper_id : lower_id,
        hit.triangle,
        bits(hit.point_m.x),
        bits(hit.point_m.y),
        bits(hit.point_m.z),
        bits(sign * triangle.geometric_outward_normal.x),
        bits(sign * triangle.geometric_outward_normal.y),
        bits(sign * triangle.geometric_outward_normal.z),
        bits(lower_fluid ? hit.segment_fraction : 1.0 - hit.segment_fraction)};
    result.push_back(
        {{std::min(lower_id, upper_id),
          std::max(lower_id, upper_id),
          hit.triangle,
          {bits(hit.point_m.x), bits(hit.point_m.y), bits(hit.point_m.z)}},
         packed});
  };

  for (int k = 0; k < kExtent.z; ++k) {
    for (int j = 0; j < kExtent.y; ++j) {
      for (int i = 1; i < kExtent.x; ++i) {
        append({i - 1, j, k}, {i, j, k});
      }
    }
  }
  for (int k = 0; k < kExtent.z; ++k) {
    for (int j = 1; j < kExtent.y; ++j) {
      for (int i = 0; i < kExtent.x; ++i) {
        append({i, j - 1, k}, {i, j, k});
      }
    }
  }
  for (int k = 1; k < kExtent.z; ++k) {
    for (int j = 0; j < kExtent.y; ++j) {
      for (int i = 0; i < kExtent.x; ++i) {
        append({i, j, k - 1}, {i, j, k});
      }
    }
  }
  std::sort(result.begin(), result.end(),
            [](const OracleLink &first, const OracleLink &second) {
              return less_key(first.key, second.key);
            });
  for (std::size_t index = 0U; index < result.size(); ++index) {
    result[index].packed[0] = index;
  }
  return result;
}

void hash_u64(std::uint64_t &hash, std::uint64_t value) {
  constexpr std::uint64_t prime = UINT64_C(1099511628211);
  for (unsigned byte = 0U; byte < 8U; ++byte) {
    hash ^= (value >> (8U * byte)) & UINT64_C(0xff);
    hash *= prime;
  }
}

std::uint64_t oracle_classification_fingerprint(
    const ImmersedSurface &surface, const SurfaceQuery &query,
    ImmersedFluidSide fluid_side, const std::vector<Real3> &global_centers) {
  constexpr std::uint64_t offset = UINT64_C(14695981039346656037);
  std::uint64_t hash = offset;
  hash_u64(hash, UINT64_C(0x48554e44434c5333));
  hash_u64(hash, surface.fingerprint());
  hash_u64(hash, query.fingerprint());
  hash_u64(hash, static_cast<std::uint64_t>(fluid_side));
  hash_u64(hash, global_centers.size());
  for (std::size_t id = 0U; id < global_centers.size(); ++id) {
    hash_u64(hash, id);
    hash_u64(hash, static_cast<std::uint64_t>(
                       query.classify(global_centers[id], fluid_side)));
  }
  return hash;
}

std::uint64_t
oracle_coverage_fingerprint(const ImmersedSurface &surface,
                            const std::vector<PackedLink> &links,
                            const std::vector<Real3> &global_centers,
                            const std::vector<double> &global_h,
                            const std::vector<std::uint64_t> &global_active) {
  constexpr std::array<std::array<double, 3>, 3> barycentric{{
      {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
  }};
  constexpr double witness_factor = 2.0 * 1.7320508075688772935;
  constexpr std::uint64_t offset = UINT64_C(14695981039346656037);
  const auto distance = [](Real3 first, Real3 second) {
    return std::hypot(first.x - second.x, first.y - second.y,
                      first.z - second.z);
  };

  HUNDUN_CHECK(global_centers.size() == global_h.size());
  std::uint64_t hash = offset;
  hash_u64(hash, UINT64_C(0x48554e4443565233));
  hash_u64(hash, surface.fingerprint());
  hash_u64(hash, surface.triangle_count());
  double maximum_distance = 0.0;
  for (std::uint64_t triangle_id = 0U; triangle_id < surface.triangle_count();
       ++triangle_id) {
    const auto &triangle_value = surface.triangle(triangle_id);
    for (std::size_t point_index = 0U; point_index < barycentric.size();
         ++point_index) {
      Real3 point{};
      for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
        point.x += barycentric[point_index][vertex] *
                   triangle_value.vertices_m[vertex].x;
        point.y += barycentric[point_index][vertex] *
                   triangle_value.vertices_m[vertex].y;
        point.z += barycentric[point_index][vertex] *
                   triangle_value.vertices_m[vertex].z;
      }

      const PackedLink *selected = nullptr;
      double selected_distance = std::numeric_limits<double>::infinity();
      for (const PackedLink &link : links) {
        HUNDUN_CHECK(link[1] < global_h.size());
        const Real3 intercept{value_from_bits(link[4]),
                              value_from_bits(link[5]),
                              value_from_bits(link[6])};
        const double candidate_distance = distance(point, intercept);
        if (candidate_distance <= 2.0 * global_h[link[1]] &&
            (candidate_distance < selected_distance ||
             (candidate_distance == selected_distance &&
              (selected == nullptr || link[0] < (*selected)[0])))) {
          selected = &link;
          selected_distance = candidate_distance;
        }
      }
      HUNDUN_CHECK(selected != nullptr);
      const std::uint64_t fluid = (*selected)[1];
      const std::uint64_t solid = (*selected)[2];
      HUNDUN_CHECK(fluid < global_centers.size());
      HUNDUN_CHECK(solid < global_centers.size());
      HUNDUN_CHECK(std::binary_search(global_active.begin(),
                                      global_active.end(), fluid));
      const double fluid_distance = distance(point, global_centers[fluid]);
      const double solid_distance = distance(point, global_centers[solid]);
      HUNDUN_CHECK(fluid_distance <= witness_factor * global_h[fluid]);
      HUNDUN_CHECK(solid_distance <= witness_factor * global_h[solid]);
      maximum_distance = std::max(maximum_distance, selected_distance);
      hash_u64(hash, triangle_id);
      hash_u64(hash, point_index);
      hash_u64(hash, (*selected)[0]);
      hash_u64(hash, fluid);
      hash_u64(hash, solid);
      hash_u64(hash, bits(selected_distance));
      hash_u64(hash, bits(fluid_distance));
      hash_u64(hash, bits(solid_distance));
    }
  }
  hash_u64(hash, bits(maximum_distance));
  return hash;
}

struct DomainEvidence final {
  std::uint64_t classification{};
  std::uint64_t coverage{};
  std::uint64_t active_fingerprint{};
  std::uint64_t owned_active{};
  std::vector<std::uint64_t> ordered_active;
  std::vector<PackedLink> links;
  std::vector<std::uint64_t> global_active;
  std::vector<PackedLink> global_links;
  std::vector<std::vector<std::uint64_t>> nested;
};

bool same_vector(const std::vector<std::uint64_t> &first,
                 const std::vector<std::uint64_t> &second) {
  if (first.size() != second.size()) {
    return false;
  }
  for (std::size_t index = 0U; index < first.size(); ++index) {
    if (first[index] != second[index]) {
      return false;
    }
  }
  return true;
}

bool same_links(const std::vector<PackedLink> &first,
                const std::vector<PackedLink> &second) {
  if (first.size() != second.size()) {
    return false;
  }
  for (std::size_t link = 0U; link < first.size(); ++link) {
    for (std::size_t field = 0U; field < first[link].size(); ++field) {
      if (first[link][field] != second[link][field]) {
        return false;
      }
    }
  }
  return true;
}

bool same_evidence(const DomainEvidence &first, const DomainEvidence &second) {
  if (first.classification != second.classification ||
      first.coverage != second.coverage ||
      first.active_fingerprint != second.active_fingerprint ||
      first.owned_active != second.owned_active ||
      !same_vector(first.ordered_active, second.ordered_active) ||
      !same_links(first.links, second.links) ||
      !same_vector(first.global_active, second.global_active) ||
      !same_links(first.global_links, second.global_links) ||
      first.nested.size() != second.nested.size()) {
    return false;
  }
  for (std::size_t outer = 0U; outer < first.nested.size(); ++outer) {
    if (!same_vector(first.nested[outer], second.nested[outer])) {
      return false;
    }
  }
  return true;
}

DomainEvidence observe(const ImmersedDomain &domain) {
  DomainEvidence result;
  result.classification = domain.classification_fingerprint();
  result.coverage = domain.surface_coverage_fingerprint();
  result.active_fingerprint = domain.active_cells().fingerprint();
  result.owned_active = domain.active_cells().owned_active_count();
  result.ordered_active = domain.active_cells().ordered_global_ids();
  for (const ImmersedLink &link : domain.links()) {
    result.links.push_back(pack(link));
  }
  result.nested.push_back(result.ordered_active);
  result.nested.emplace_back();
  for (const PackedLink &link : result.links) {
    result.nested.back().push_back(link[0]);
  }
  return result;
}

void check_equality_oracle(const DomainEvidence &source) {
  DomainEvidence exact = source;
  HUNDUN_CHECK(same_evidence(source, exact));
  DomainEvidence ordinary = source;
  ordinary.classification ^= UINT64_C(1);
  HUNDUN_CHECK(!same_evidence(source, ordinary));
  DomainEvidence nested = source;
  HUNDUN_CHECK(!nested.nested.empty() && !nested.nested.front().empty());
  nested.nested.front().front() ^= UINT64_C(1);
  HUNDUN_CHECK(!same_evidence(source, nested));
  DomainEvidence inner_size = source;
  inner_size.nested.emplace_back();
  HUNDUN_CHECK(!same_evidence(source, inner_size));
}

void require_same_on_ranks(std::uint64_t value, const MpiContext &mpi) {
  std::uint64_t minimum = value;
  std::uint64_t maximum = value;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(minimum == maximum);
}

DomainEvidence run_success_case(const MpiContext &mpi,
                                const std::filesystem::path &surface_path,
                                bool warped, ImmersedFluidSide fluid_side,
                                Int3 grid) {
  const auto config = closed_config(mpi.size(), grid);
  auto decomposition = StructuredDecomposition::create(
      mpi, kExtent, {false, false, false}, DecompositionOptions{grid});
  MeshTopology topology(decomposition);
  std::unique_ptr<MeshGeometry> geometry;
  if (warped) {
    geometry = std::make_unique<MeshGeometry>(
        topology, AnalyticWarpedBoxMapping{kOrigin, kLength, kWarp});
  } else {
    geometry = std::make_unique<MeshGeometry>(
        topology, UniformBoxMapping{kOrigin, kLength});
  }
  BoundaryRegistry boundaries = BoundaryRegistry::create(config, topology);
  const ImmersedSurface surface =
      ImmersedSurface::load_collective(surface_path, 1.0, mpi, 0);
  const SurfaceQuery query = SurfaceQuery::create(surface);
  const ImmersedDomain domain = ImmersedDomain::create(
      surface, query, fluid_side, topology, *geometry, boundaries, mpi);
  const auto global_centers =
      gather_owned_cell_centers(topology, *geometry, mpi);
  const auto global_h = gather_owned_cell_sizes(topology, *geometry, mpi);

  require_same_on_ranks(domain.classification_fingerprint(), mpi);
  require_same_on_ranks(domain.surface_coverage_fingerprint(), mpi);
  HUNDUN_CHECK(domain.classification_fingerprint() != 0U);
  HUNDUN_CHECK(domain.surface_coverage_fingerprint() != 0U);
  HUNDUN_CHECK(domain.active_cells().fingerprint() != 0U);
  HUNDUN_CHECK(domain.classification_fingerprint() ==
               oracle_classification_fingerprint(surface, query, fluid_side,
                                                 global_centers));

  std::vector<std::uint64_t> expected_local_ids;
  std::size_t expected_owned_active_count = 0U;
  for (std::size_t local = 0U; local < topology.local_cell_count(); ++local) {
    const CellRegion expected =
        query.classify(geometry->cell_center_m(local), fluid_side);
    HUNDUN_CHECK(domain.region(local) == expected);
    const bool active = expected == CellRegion::fluid;
    HUNDUN_CHECK(domain.active_cells().active(local) == active);
    const auto active_index = domain.active_cells().active_index(local);
    HUNDUN_CHECK(active_index.has_value() == active);
    if (active) {
      expected_owned_active_count +=
          local < topology.owned_cell_count() ? 1U : 0U;
      HUNDUN_CHECK(*active_index <
                   domain.active_cells().ordered_global_ids().size());
      HUNDUN_CHECK(domain.active_cells().ordered_global_ids()[*active_index] ==
                   topology.global_cell_id(local));
      expected_local_ids.push_back(topology.global_cell_id(local));
    }
  }
  HUNDUN_CHECK(domain.active_cells().owned_active_count() ==
               expected_owned_active_count);
  const auto owned_end =
      expected_local_ids.begin() +
      static_cast<std::ptrdiff_t>(domain.active_cells().owned_active_count());
  HUNDUN_CHECK(std::is_sorted(expected_local_ids.begin(), owned_end));
  HUNDUN_CHECK(std::is_sorted(owned_end, expected_local_ids.end()));
  HUNDUN_CHECK(same_vector(expected_local_ids,
                           domain.active_cells().ordered_global_ids()));
  HUNDUN_CHECK(domain.active_cells().local_active_count() ==
               expected_local_ids.size());
  expect_error(
      [&] { static_cast<void>(domain.region(topology.local_cell_count())); },
      "local cell");
  expect_error(
      [&] {
        static_cast<void>(
            domain.active_cells().active(topology.local_cell_count()));
      },
      "local cell");
  std::uint64_t expected_active_fingerprint = UINT64_C(14695981039346656037);
  hash_u64(expected_active_fingerprint, UINT64_C(0x48554e4441435433));
  hash_u64(expected_active_fingerprint, topology.local_cell_count());
  hash_u64(expected_active_fingerprint,
           domain.active_cells().owned_active_count());
  for (const std::uint64_t id : expected_local_ids) {
    hash_u64(expected_active_fingerprint, id);
  }
  HUNDUN_CHECK(domain.active_cells().fingerprint() ==
               expected_active_fingerprint);

  std::vector<std::uint64_t> local_owned_active;
  for (std::size_t row = 0U; row < domain.active_cells().owned_active_count();
       ++row) {
    local_owned_active.push_back(
        domain.active_cells().ordered_global_ids()[row]);
  }
  auto global_active = gather_u64(local_owned_active, mpi);
  std::sort(global_active.begin(), global_active.end());
  std::vector<std::uint64_t> expected_global_active;
  for (int k = 0; k < kExtent.z; ++k) {
    for (int j = 0; j < kExtent.y; ++j) {
      for (int i = 0; i < kExtent.x; ++i) {
        const Int3 cell{i, j, k};
        if (query.classify(global_centers.at(topology.global_cell_id(cell)),
                           fluid_side) == CellRegion::fluid) {
          expected_global_active.push_back(topology.global_cell_id(cell));
        }
      }
    }
  }
  HUNDUN_CHECK(same_vector(global_active, expected_global_active));

  std::vector<std::uint64_t> local_packed;
  HUNDUN_CHECK(
      std::is_sorted(domain.links().begin(), domain.links().end(),
                     [](const ImmersedLink &first, const ImmersedLink &second) {
                       return first.id < second.id;
                     }));
  for (const ImmersedLink &link : domain.links()) {
    const PackedLink fields = pack(link);
    local_packed.insert(local_packed.end(), fields.begin(), fields.end());
    const auto fluid = topology.find_local_cell(link.fluid_cell);
    const auto solid = topology.find_local_cell(link.solid_cell);
    HUNDUN_CHECK(fluid.has_value() && solid.has_value());
    bool published_by_face_owner = false;
    for (std::size_t face = 0U; face < topology.local_face_count(); ++face) {
      const auto neighbour = topology.neighbour(face);
      if (!neighbour.has_value()) {
        continue;
      }
      const std::uint64_t owner_global =
          topology.global_cell_id(topology.owner(face));
      const std::uint64_t neighbour_global =
          topology.global_cell_id(*neighbour);
      if ((owner_global == link.fluid_cell &&
           neighbour_global == link.solid_cell) ||
          (owner_global == link.solid_cell &&
           neighbour_global == link.fluid_cell)) {
        HUNDUN_CHECK(!topology.patch_id(face).has_value());
        HUNDUN_CHECK(topology.face_ownership(face) ==
                     hundun::mesh::EntityOwnership::owned);
        published_by_face_owner = true;
        break;
      }
    }
    HUNDUN_CHECK(published_by_face_owner);
    const auto hits = query.segment_intersections(
        geometry->cell_center_m(*fluid), geometry->cell_center_m(*solid));
    HUNDUN_CHECK(hits.size() == 1U);
    HUNDUN_CHECK(hits.front().triangle == link.triangle);
    HUNDUN_CHECK(bits(hits.front().point_m.x) == bits(link.wall_intercept_m.x));
    HUNDUN_CHECK(bits(hits.front().point_m.y) == bits(link.wall_intercept_m.y));
    HUNDUN_CHECK(bits(hits.front().point_m.z) == bits(link.wall_intercept_m.z));
    HUNDUN_CHECK(link.fluid_to_wall_fraction > 0.0 &&
                 link.fluid_to_wall_fraction < 1.0);
  }
  auto gathered_packed = gather_u64(local_packed, mpi);
  HUNDUN_CHECK(gathered_packed.size() % PackedLink{}.size() == 0U);
  std::vector<PackedLink> actual_links(gathered_packed.size() /
                                       PackedLink{}.size());
  for (std::size_t link = 0U; link < actual_links.size(); ++link) {
    std::copy_n(gathered_packed.begin() +
                    static_cast<std::ptrdiff_t>(link * PackedLink{}.size()),
                PackedLink{}.size(), actual_links[link].begin());
  }
  std::sort(actual_links.begin(), actual_links.end(),
            [](const PackedLink &first, const PackedLink &second) {
              return first[0] < second[0];
            });
  const auto expected_links =
      oracle_links(topology, surface, query, fluid_side, global_centers);
  HUNDUN_CHECK(actual_links.size() == expected_links.size());
  for (std::size_t link = 0U; link < actual_links.size(); ++link) {
    HUNDUN_CHECK(actual_links[link][0] == link);
    for (std::size_t field = 0U; field < actual_links[link].size(); ++field) {
      if (actual_links[link][field] != expected_links[link].packed[field]) {
        throw std::runtime_error(
            "immersed link oracle mismatch at link " + std::to_string(link) +
            " field " + std::to_string(field) +
            ": actual=" + std::to_string(actual_links[link][field]) +
            " expected=" + std::to_string(expected_links[link].packed[field]));
      }
    }
  }
  HUNDUN_CHECK(domain.surface_coverage_fingerprint() ==
               oracle_coverage_fingerprint(surface, actual_links,
                                           global_centers, global_h,
                                           global_active));

  const ImmersedDomain repeated = ImmersedDomain::create(
      surface, query, fluid_side, topology, *geometry, boundaries, mpi);
  const DomainEvidence first = observe(domain);
  const DomainEvidence second = observe(repeated);
  DomainEvidence first_with_global = first;
  DomainEvidence second_with_global = second;
  first_with_global.global_active = global_active;
  second_with_global.global_active = global_active;
  first_with_global.global_links = actual_links;
  second_with_global.global_links = actual_links;
  HUNDUN_CHECK(same_evidence(first_with_global, second_with_global));
  check_equality_oracle(first_with_global);
  return first_with_global;
}

void success_cases(const MpiContext &mpi, const SharedDirectory &directory,
                   bool alternate_grid) {
  const auto path = write_cube(directory, mpi);
  const Int3 first_grid = process_grid(mpi.size(), alternate_grid);
  for (const bool warped : {false, true}) {
    for (const ImmersedFluidSide side :
         {ImmersedFluidSide::outside, ImmersedFluidSide::inside}) {
      static_cast<void>(run_success_case(mpi, path, warped, side, first_grid));
    }
  }
  if (mpi.size() == 4 && !alternate_grid) {
    const DomainEvidence first = run_success_case(
        mpi, path, false, ImmersedFluidSide::outside, process_grid(4, false));
    const DomainEvidence second = run_success_case(
        mpi, path, false, ImmersedFluidSide::outside, process_grid(4, true));
    HUNDUN_CHECK(first.classification == second.classification);
    HUNDUN_CHECK(first.coverage == second.coverage);
    HUNDUN_CHECK(same_vector(first.global_active, second.global_active));
    HUNDUN_CHECK(same_links(first.global_links, second.global_links));
  }
}

void failure_cases(const MpiContext &mpi, const SharedDirectory &directory) {
  HUNDUN_CHECK(mpi.size() >= 2);
  const auto path = write_cube(directory, mpi);
  const Int3 grid = process_grid(mpi.size(), false);
  const auto config = closed_config(mpi.size(), grid);
  auto decomposition = StructuredDecomposition::create(
      mpi, kExtent, {false, false, false}, DecompositionOptions{grid});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(topology, UniformBoxMapping{kOrigin, kLength});
  BoundaryRegistry boundaries = BoundaryRegistry::create(config, topology);
  const ImmersedSurface surface =
      ImmersedSurface::load_collective(path, 1.0, mpi, 0);
  const SurfaceQuery query = SurfaceQuery::create(surface);
  const auto other_query_path = write_surface(
      directory, mpi, "task3-other-query-surface", cube_between(0.20, 0.80));
  const ImmersedSurface other_query_surface =
      ImmersedSurface::load_collective(other_query_path, 1.0, mpi, 0);
  const SurfaceQuery other_query = SurfaceQuery::create(other_query_surface);
  const std::uint64_t surface_before = surface.fingerprint();
  const std::uint64_t query_before = query.fingerprint();
  const std::uint64_t other_surface_before = other_query_surface.fingerprint();
  const std::uint64_t other_query_before = other_query.fingerprint();
  const std::size_t cells_before = topology.local_cell_count();
  const std::size_t faces_before = topology.local_face_count();
  const Real3 origin_before = geometry.origin_m();
  const bool open_before = boundaries.open_domain();
  const auto require_unchanged = [&] {
    HUNDUN_CHECK(surface.fingerprint() == surface_before);
    HUNDUN_CHECK(query.fingerprint() == query_before);
    HUNDUN_CHECK(other_query_surface.fingerprint() == other_surface_before);
    HUNDUN_CHECK(other_query.fingerprint() == other_query_before);
    HUNDUN_CHECK(topology.local_cell_count() == cells_before);
    HUNDUN_CHECK(topology.local_face_count() == faces_before);
    HUNDUN_CHECK(bits(geometry.origin_m().x) == bits(origin_before.x));
    HUNDUN_CHECK(bits(geometry.origin_m().y) == bits(origin_before.y));
    HUNDUN_CHECK(bits(geometry.origin_m().z) == bits(origin_before.z));
    HUNDUN_CHECK(boundaries.open_domain() == open_before);
  };

  expect_collective_error(
      [&] {
        const auto side = mpi.rank() == 1 ? static_cast<ImmersedFluidSide>(255)
                                          : ImmersedFluidSide::outside;
        static_cast<void>(ImmersedDomain::create(surface, query, side, topology,
                                                 geometry, boundaries, mpi));
      },
      mpi, "fluid_side", 1);
  require_unchanged();

  expect_collective_error(
      [&] {
        const SurfaceQuery &selected = mpi.rank() == 1 ? other_query : query;
        static_cast<void>(ImmersedDomain::create(
            surface, selected, ImmersedFluidSide::outside, topology, geometry,
            boundaries, mpi));
      },
      mpi, "query_surface_compatible", 1);
  require_unchanged();

  constexpr Int3 other_extent{13, 12, 12};
  auto other_decomposition = StructuredDecomposition::create(
      mpi, other_extent, {false, false, false}, DecompositionOptions{grid});
  MeshTopology other_topology(other_decomposition);
  MeshGeometry other_geometry(other_topology,
                              UniformBoxMapping{kOrigin, kLength});
  expect_collective_error(
      [&] {
        const MeshGeometry &selected =
            mpi.rank() == 1 ? other_geometry : geometry;
        static_cast<void>(
            ImmersedDomain::create(surface, query, ImmersedFluidSide::outside,
                                   topology, selected, boundaries, mpi));
      },
      mpi, "geometry", 1);
  require_unchanged();

  const auto expect_fixture_failure = [&](const std::filesystem::path &fixture,
                                          ImmersedFluidSide side, Int3 extent,
                                          const std::string &marker) {
    auto fixture_decomposition = StructuredDecomposition::create(
        mpi, extent, {false, false, false}, DecompositionOptions{grid});
    MeshTopology fixture_topology(fixture_decomposition);
    MeshGeometry fixture_geometry(fixture_topology,
                                  UniformBoxMapping{kOrigin, kLength});
    BoundaryRegistry fixture_boundaries = BoundaryRegistry::create(
        closed_config(mpi.size(), grid, extent), fixture_topology);
    ImmersedSurface fixture_surface = [&] {
      try {
        return ImmersedSurface::load_collective(fixture, 1.0, mpi, 0);
      } catch (const hundun::runtime::Error &error) {
        throw std::runtime_error("fixture '" + fixture.filename().string() +
                                 "' failed to load: " + error.what());
      }
    }();
    const SurfaceQuery fixture_query = SurfaceQuery::create(fixture_surface);
    const std::uint64_t surface_before = fixture_surface.fingerprint();
    const std::uint64_t query_before = fixture_query.fingerprint();
    const std::size_t cells_before = fixture_topology.local_cell_count();
    const std::size_t faces_before = fixture_topology.local_face_count();
    const Real3 origin_before = fixture_geometry.origin_m();
    const bool open_before = fixture_boundaries.open_domain();
    expect_collective_error(
        [&] {
          static_cast<void>(ImmersedDomain::create(
              fixture_surface, fixture_query, side, fixture_topology,
              fixture_geometry, fixture_boundaries, mpi));
        },
        mpi, marker, 0);
    HUNDUN_CHECK(fixture_surface.fingerprint() == surface_before);
    HUNDUN_CHECK(fixture_query.fingerprint() == query_before);
    HUNDUN_CHECK(fixture_topology.local_cell_count() == cells_before);
    HUNDUN_CHECK(fixture_topology.local_face_count() == faces_before);
    HUNDUN_CHECK(bits(fixture_geometry.origin_m().x) == bits(origin_before.x));
    HUNDUN_CHECK(bits(fixture_geometry.origin_m().y) == bits(origin_before.y));
    HUNDUN_CHECK(bits(fixture_geometry.origin_m().z) == bits(origin_before.z));
    HUNDUN_CHECK(fixture_boundaries.open_domain() == open_before);
  };

  expect_fixture_failure(write_surface(directory, mpi, "task3-zero-fluid",
                                       cube_between(-0.2, 1.2)),
                         ImmersedFluidSide::outside, kExtent,
                         "fluid_cell_nonzero");
  expect_fixture_failure(write_surface(directory, mpi, "task3-zero-solid",
                                       cube_between(0.07, 0.095)),
                         ImmersedFluidSide::outside, kExtent,
                         "solid_cell_nonzero");
  expect_fixture_failure(write_surface(directory, mpi, "task3-single-cell",
                                       cube_between(0.52, 0.56)),
                         ImmersedFluidSide::outside, kExtent,
                         "surface_resolved_by_multiple_centres");
  expect_fixture_failure(
      write_surface(directory, mpi, "task3-domain-separation",
                    cube_between(0.10, 0.75)),
      ImmersedFluidSide::outside, kExtent, "surface_domain_separation");
  expect_fixture_failure(write_surface(directory, mpi,
                                       "task3-coverage-distance",
                                       cube_between(0.20, 0.80)),
                         ImmersedFluidSide::outside, Int3{64, 3, 64},
                         "coverage_point_link_distance");
  const auto witness_fixture = write_surface(
      directory, mpi, "task3-coverage-witness",
      box_between(Real3{0.45, 0.25, 0.25}, Real3{0.75, 0.75, 0.75}));
  expect_fixture_failure(witness_fixture, ImmersedFluidSide::outside,
                         Int3{3, 64, 64}, "coverage_solid_witness");
  expect_fixture_failure(witness_fixture, ImmersedFluidSide::inside,
                         Int3{3, 64, 64}, "coverage_fluid_witness");
  const auto concave_fixture =
      write_surface(directory, mpi, "task3-multiple-centre-intersections",
                    connected_concave_extrusion());
  expect_fixture_failure(concave_fixture, ImmersedFluidSide::outside,
                         Int3{3, 11, 11}, "link_single_intersection");

  constexpr Int3 concave_extent{3, 11, 11};
  auto concave_decomposition = StructuredDecomposition::create(
      mpi, concave_extent, {false, false, false}, DecompositionOptions{grid});
  MeshTopology concave_topology(concave_decomposition);
  MeshGeometry concave_geometry(concave_topology,
                                UniformBoxMapping{kOrigin, kLength});
  BoundaryRegistry concave_boundaries = BoundaryRegistry::create(
      open_config(mpi.size(), grid, concave_extent), concave_topology);
  const ImmersedSurface concave_surface =
      ImmersedSurface::load_collective(concave_fixture, 1.0, mpi, 0);
  const SurfaceQuery concave_query = SurfaceQuery::create(concave_surface);
  expect_collective_error(
      [&] {
        static_cast<void>(ImmersedDomain::create(
            concave_surface, concave_query, ImmersedFluidSide::inside,
            concave_topology, concave_geometry, concave_boundaries, mpi));
      },
      mpi, "link_single_intersection", 0);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = MpiContext::duplicate(MPI_COMM_WORLD);
  SharedDirectory directory(mpi);
  const std::string mode = argc > 1 ? argv[1] : "success";
  if (mode == "success") {
    const bool alternate = argc > 2 && std::string(argv[2]) == "alternate";
    success_cases(mpi, directory, alternate);
  } else if (mode == "failures") {
    failure_cases(mpi, directory);
  } else {
    HUNDUN_CHECK(false);
  }
  return 0;
}
