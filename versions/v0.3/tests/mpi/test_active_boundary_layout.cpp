// SPDX-License-Identifier: Apache-2.0

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_surface_query.hpp"
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
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using hundun::boundary::BoundaryKind;
using hundun::boundary::BoundaryRegistry;
using hundun::config::BoundaryType;
using hundun::config::DensityModel;
using hundun::config::FlowBoundaryConfig;
using hundun::config::FlowCaseConfig;
using hundun::config::ImmersedFluidSide;
using hundun::config::InletThermalAuthority;
using hundun::config::MeshMapping;
using hundun::config::PatchName;
using hundun::config::SimulationType;
using hundun::config::TimeMode;
using hundun::immersed::ActiveBoundaryLayout;
using hundun::immersed::ImmersedDomain;
using hundun::immersed::ImmersedSurface;
using hundun::immersed::SurfaceQuery;
using hundun::mesh::GlobalFaceId;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::Real3;
using hundun::runtime::StructuredDecomposition;

constexpr Int3 kExtent{12, 12, 12};
constexpr Real3 kOrigin{0.0, 0.0, 0.0};
constexpr Real3 kLength{1.0, 1.0, 1.0};
constexpr std::array<PatchName, 6> kPatchNames{
    PatchName::x_min, PatchName::x_max, PatchName::y_min,
    PatchName::y_max, PatchName::z_min, PatchName::z_max};

Int3 process_grid(int ranks, bool alternate) {
  switch (ranks) {
  case 1:
    return {1, 1, 1};
  case 2:
    return alternate ? Int3{1, 2, 1} : Int3{2, 1, 1};
  case 4:
    return alternate ? Int3{4, 1, 1} : Int3{2, 2, 1};
  default:
    throw hundun::runtime::Error("unsupported active-boundary rank count");
  }
}

FlowBoundaryConfig boundary(std::size_t stable_id, BoundaryType type) {
  FlowBoundaryConfig result{};
  result.patch = kPatchNames.at(stable_id);
  result.type = type;
  if (type == BoundaryType::velocity_inlet) {
    result.velocity_m_per_s = Real3{1.0, 0.0, 0.0};
    result.thermal_authority = InletThermalAuthority::enthalpy;
    result.enthalpy_J_per_kg = 1.0;
    result.scalar_values = std::vector<hundun::config::InletScalarValue>{};
  } else if (type == BoundaryType::pressure_outlet) {
    result.pressure_perturbation_pa = 0.0;
  }
  return result;
}

FlowCaseConfig make_config(int ranks, Int3 grid,
                           const std::array<BoundaryType, 6> &kinds) {
  FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "stage3_task3_active_boundary";
  config.simulation_type = SimulationType::variable_density_flow;
  config.density_model = DensityModel::constant;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = grid;
  config.mesh.cells = kExtent;
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
  for (std::size_t patch = 0U; patch < kinds.size(); ++patch) {
    config.boundaries[patch] = boundary(patch, kinds[patch]);
  }
  return config;
}

std::array<bool, 3> periodic_axes(const std::array<BoundaryType, 6> &kinds) {
  return {
      kinds[0] == BoundaryType::periodic && kinds[1] == BoundaryType::periodic,
      kinds[2] == BoundaryType::periodic && kinds[3] == BoundaryType::periodic,
      kinds[4] == BoundaryType::periodic && kinds[5] == BoundaryType::periodic};
}

class SharedDirectory final {
public:
  explicit SharedDirectory(const MpiContext &mpi) : mpi_(&mpi) {
    std::array<char, 256> storage{};
    if (mpi.rank() == 0) {
      const auto stamp =
          std::chrono::steady_clock::now().time_since_epoch().count();
      const std::string value =
          (std::filesystem::temp_directory_path() /
           ("hundun-stage3-active-boundary-" + std::to_string(stamp)))
              .string();
      HUNDUN_CHECK(value.size() + 1U <= storage.size());
      std::copy(value.begin(), value.end(), storage.begin());
      std::filesystem::create_directories(value);
    }
    HUNDUN_CHECK(MPI_Bcast(storage.data(), static_cast<int>(storage.size()),
                           MPI_CHAR, 0, mpi.comm()) == MPI_SUCCESS);
    path_ = storage.data();
    HUNDUN_CHECK(MPI_Barrier(mpi.comm()) == MPI_SUCCESS);
  }

  ~SharedDirectory() {
    MPI_Barrier(mpi_->comm());
    if (mpi_->rank() == 0) {
      std::error_code error;
      std::filesystem::remove_all(path_, error);
    }
  }

  const std::filesystem::path &path() const noexcept { return path_; }

private:
  const MpiContext *mpi_{};
  std::filesystem::path path_;
};

std::vector<hundun::test::StlFixtureTriangle> box_triangles(Real3 minimum,
                                                            Real3 maximum) {
  auto triangles = hundun::test::outward_cube();
  for (auto &triangle : triangles) {
    for (Real3 &vertex : triangle.vertices) {
      vertex.x = minimum.x + (maximum.x - minimum.x) * vertex.x;
      vertex.y = minimum.y + (maximum.y - minimum.y) * vertex.y;
      vertex.z = minimum.z + (maximum.z - minimum.z) * vertex.z;
    }
  }
  return triangles;
}

std::filesystem::path write_box(const SharedDirectory &directory,
                                const MpiContext &mpi, std::string name,
                                Real3 minimum, Real3 maximum) {
  const auto path = directory.path() / (name + ".stl");
  if (mpi.rank() == 0) {
    hundun::test::write_text(
        path, hundun::test::ascii_stl(box_triangles(minimum, maximum), name));
  }
  HUNDUN_CHECK(MPI_Barrier(mpi.comm()) == MPI_SUCCESS);
  return path;
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
  std::vector<std::uint64_t> result(static_cast<std::size_t>(total));
  HUNDUN_CHECK(MPI_Allgatherv(local.data(), local_count, MPI_UINT64_T,
                              result.data(), counts.data(),
                              displacements.data(), MPI_UINT64_T,
                              mpi.comm()) == MPI_SUCCESS);
  return result;
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
  std::string message;
  try {
    operation();
  } catch (const hundun::runtime::Error &error) {
    message = error.what();
  }
  HUNDUN_CHECK(!message.empty());
  HUNDUN_CHECK(message.find(marker) != std::string::npos);
  HUNDUN_CHECK(message.find("lowest failing rank " +
                            std::to_string(expected_failing_rank)) !=
               std::string::npos);
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
}

struct BoundaryEvidence final {
  std::array<std::vector<GlobalFaceId>, 6> patch_faces;
  bool open{};
  bool pressure_reference{};
  std::uint64_t fingerprint{};
};

bool same_faces(const std::vector<GlobalFaceId> &first,
                const std::vector<GlobalFaceId> &second) {
  if (first.size() != second.size()) {
    return false;
  }
  for (std::size_t face = 0U; face < first.size(); ++face) {
    if (first[face] != second[face]) {
      return false;
    }
  }
  return true;
}

bool same_evidence(const BoundaryEvidence &first,
                   const BoundaryEvidence &second) {
  if (first.open != second.open ||
      first.pressure_reference != second.pressure_reference ||
      first.fingerprint != second.fingerprint) {
    return false;
  }
  for (std::size_t patch = 0U; patch < first.patch_faces.size(); ++patch) {
    if (!same_faces(first.patch_faces[patch], second.patch_faces[patch])) {
      return false;
    }
  }
  return true;
}

BoundaryEvidence observe(const ActiveBoundaryLayout &layout) {
  BoundaryEvidence result{};
  for (std::uint32_t patch = 0U; patch < 6U; ++patch) {
    result.patch_faces[patch] = layout.patch_faces(patch);
  }
  result.open = layout.open_domain();
  result.pressure_reference = layout.has_pressure_reference();
  result.fingerprint = layout.fingerprint();
  return result;
}

BoundaryEvidence run_case(const MpiContext &mpi,
                          const std::filesystem::path &surface_path,
                          const std::array<BoundaryType, 6> &kinds,
                          ImmersedFluidSide side, Int3 grid) {
  const auto config = make_config(mpi.size(), grid, kinds);
  auto decomposition = StructuredDecomposition::create(
      mpi, kExtent, periodic_axes(kinds), DecompositionOptions{grid});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(topology,
                        hundun::mesh::UniformBoxMapping{kOrigin, kLength});
  BoundaryRegistry boundaries = BoundaryRegistry::create(config, topology);
  const ImmersedSurface surface =
      ImmersedSurface::load_collective(surface_path, 1.0, mpi, 0);
  const SurfaceQuery query = SurfaceQuery::create(surface);
  const ImmersedDomain domain = ImmersedDomain::create(
      surface, query, side, topology, geometry, boundaries, mpi);
  const ActiveBoundaryLayout &layout = domain.active_boundaries();

  const bool expected_open =
      std::find(kinds.begin(), kinds.end(), BoundaryType::velocity_inlet) !=
      kinds.end();
  const bool expected_reference =
      std::find(kinds.begin(), kinds.end(), BoundaryType::pressure_outlet) !=
      kinds.end();
  HUNDUN_CHECK(layout.open_domain() == expected_open);
  HUNDUN_CHECK(layout.has_pressure_reference() == expected_reference);
  require_same_on_ranks(layout.fingerprint(), mpi);
  HUNDUN_CHECK(layout.fingerprint() != 0U);

  for (std::uint32_t patch = 0U; patch < 6U; ++patch) {
    std::vector<std::uint64_t> local_expected;
    for (const auto face : topology.patch(patch).local_faces()) {
      if (domain.region(topology.owner(face)) ==
          hundun::immersed::CellRegion::fluid) {
        local_expected.push_back(topology.global_face_id(face));
      }
    }
    auto expected = gather_u64(local_expected, mpi);
    std::sort(expected.begin(), expected.end());
    HUNDUN_CHECK(same_faces(expected, layout.patch_faces(patch)));
    HUNDUN_CHECK(std::is_sorted(layout.patch_faces(patch).begin(),
                                layout.patch_faces(patch).end()));
    if (kinds[patch] == BoundaryType::periodic) {
      const std::uint32_t paired_patch =
          patch % 2U == 0U ? patch + 1U : patch - 1U;
      for (const auto face : topology.patch(patch).local_faces()) {
        if (domain.region(topology.owner(face)) !=
            hundun::immersed::CellRegion::fluid) {
          continue;
        }
        const auto pair = topology.periodic_pair(face);
        HUNDUN_CHECK(pair.has_value());
        HUNDUN_CHECK(
            std::binary_search(layout.patch_faces(paired_patch).begin(),
                               layout.patch_faces(paired_patch).end(), *pair));
      }
    }
  }
  expect_error([&] { static_cast<void>(layout.patch_faces(6U)); },
               "stable patch ID");

  const ImmersedDomain repeated = ImmersedDomain::create(
      surface, query, side, topology, geometry, boundaries, mpi);
  const BoundaryEvidence first = observe(layout);
  const BoundaryEvidence second = observe(repeated.active_boundaries());
  HUNDUN_CHECK(same_evidence(first, second));
  BoundaryEvidence mutated = first;
  mutated.fingerprint ^= UINT64_C(1);
  HUNDUN_CHECK(!same_evidence(first, mutated));
  BoundaryEvidence nested_mutation = first;
  nested_mutation.patch_faces[0].push_back(
      std::numeric_limits<GlobalFaceId>::max());
  HUNDUN_CHECK(!same_evidence(first, nested_mutation));
  return first;
}

void success_cases(const MpiContext &mpi, const SharedDirectory &directory,
                   bool alternate) {
  const auto cube = write_box(directory, mpi, "active-boundary-cube",
                              Real3{0.25, 0.25, 0.25}, Real3{0.75, 0.75, 0.75});
  const Int3 grid = process_grid(mpi.size(), alternate);
  const std::array<BoundaryType, 6> walls{
      BoundaryType::no_slip_wall, BoundaryType::no_slip_wall,
      BoundaryType::no_slip_wall, BoundaryType::no_slip_wall,
      BoundaryType::no_slip_wall, BoundaryType::no_slip_wall};
  const std::array<BoundaryType, 6> open{
      BoundaryType::velocity_inlet, BoundaryType::pressure_outlet,
      BoundaryType::no_slip_wall,   BoundaryType::no_slip_wall,
      BoundaryType::no_slip_wall,   BoundaryType::no_slip_wall};
  const std::array<BoundaryType, 6> x_periodic{
      BoundaryType::periodic,     BoundaryType::periodic,
      BoundaryType::no_slip_wall, BoundaryType::no_slip_wall,
      BoundaryType::no_slip_wall, BoundaryType::no_slip_wall};

  const BoundaryEvidence closed_outside =
      run_case(mpi, cube, walls, ImmersedFluidSide::outside, grid);
  for (const auto &faces : closed_outside.patch_faces) {
    HUNDUN_CHECK(!faces.empty());
  }
  HUNDUN_CHECK(!closed_outside.open && !closed_outside.pressure_reference);

  const BoundaryEvidence open_outside =
      run_case(mpi, cube, open, ImmersedFluidSide::outside, grid);
  HUNDUN_CHECK(open_outside.open && open_outside.pressure_reference);
  HUNDUN_CHECK(!open_outside.patch_faces[0].empty());
  HUNDUN_CHECK(!open_outside.patch_faces[1].empty());

  const BoundaryEvidence inside =
      run_case(mpi, cube, walls, ImmersedFluidSide::inside, grid);
  for (const auto &faces : inside.patch_faces) {
    HUNDUN_CHECK(faces.empty());
  }
  HUNDUN_CHECK(!inside.open && !inside.pressure_reference);

  const BoundaryEvidence inactive_periodic =
      run_case(mpi, cube, x_periodic, ImmersedFluidSide::inside, grid);
  HUNDUN_CHECK(inactive_periodic.patch_faces[0].empty());
  HUNDUN_CHECK(inactive_periodic.patch_faces[1].empty());
  HUNDUN_CHECK(inside.fingerprint != inactive_periodic.fingerprint);

  const BoundaryEvidence active_periodic =
      run_case(mpi, cube, x_periodic, ImmersedFluidSide::outside, grid);
  HUNDUN_CHECK(!active_periodic.patch_faces[0].empty());
  HUNDUN_CHECK(!active_periodic.patch_faces[1].empty());
  HUNDUN_CHECK(closed_outside.fingerprint != active_periodic.fingerprint);

  if (mpi.size() == 4 && !alternate) {
    const BoundaryEvidence other = run_case(
        mpi, cube, walls, ImmersedFluidSide::outside, process_grid(4, true));
    HUNDUN_CHECK(closed_outside.fingerprint == other.fingerprint);
    for (std::size_t patch = 0U; patch < 6U; ++patch) {
      HUNDUN_CHECK(same_faces(closed_outside.patch_faces[patch],
                              other.patch_faces[patch]));
    }
  }
}

void failure_cases(const MpiContext &mpi, const SharedDirectory &directory) {
  HUNDUN_CHECK(mpi.size() >= 2);
  const Int3 grid = process_grid(mpi.size(), false);
  const auto cube = write_box(directory, mpi, "active-boundary-failure-cube",
                              Real3{0.25, 0.25, 0.25}, Real3{0.75, 0.75, 0.75});
  const auto near_x_min =
      write_box(directory, mpi, "active-boundary-near-x-min",
                Real3{-0.1, 0.25, 0.25}, Real3{0.4, 0.75, 0.75});
  const auto x_min_slab =
      write_box(directory, mpi, "active-boundary-x-min-slab",
                Real3{-0.1, -0.1, -0.1}, Real3{0.4, 1.1, 1.1});
  const std::array<BoundaryType, 6> open{
      BoundaryType::velocity_inlet, BoundaryType::pressure_outlet,
      BoundaryType::no_slip_wall,   BoundaryType::no_slip_wall,
      BoundaryType::no_slip_wall,   BoundaryType::no_slip_wall};
  const std::array<BoundaryType, 6> walls{
      BoundaryType::no_slip_wall, BoundaryType::no_slip_wall,
      BoundaryType::no_slip_wall, BoundaryType::no_slip_wall,
      BoundaryType::no_slip_wall, BoundaryType::no_slip_wall};
  const std::array<BoundaryType, 6> x_periodic{
      BoundaryType::periodic,     BoundaryType::periodic,
      BoundaryType::no_slip_wall, BoundaryType::no_slip_wall,
      BoundaryType::no_slip_wall, BoundaryType::no_slip_wall};

  const auto run_failure = [&](const std::filesystem::path &path,
                               const std::array<BoundaryType, 6> &kinds,
                               ImmersedFluidSide side,
                               const std::string &marker) {
    auto decomposition = StructuredDecomposition::create(
        mpi, kExtent, periodic_axes(kinds), DecompositionOptions{grid});
    MeshTopology topology(decomposition);
    MeshGeometry geometry(topology,
                          hundun::mesh::UniformBoxMapping{kOrigin, kLength});
    BoundaryRegistry boundaries = BoundaryRegistry::create(
        make_config(mpi.size(), grid, kinds), topology);
    const ImmersedSurface surface =
        ImmersedSurface::load_collective(path, 1.0, mpi, 0);
    const SurfaceQuery query = SurfaceQuery::create(surface);
    const std::uint64_t surface_before = surface.fingerprint();
    const std::uint64_t query_before = query.fingerprint();
    const std::size_t cells_before = topology.local_cell_count();
    const std::size_t faces_before = topology.local_face_count();
    const bool open_before = boundaries.open_domain();
    std::array<BoundaryKind, 6> kinds_before{};
    for (std::uint32_t patch = 0U; patch < 6U; ++patch) {
      kinds_before[patch] = boundaries.patch(patch).kind();
    }
    expect_collective_error(
        [&] {
          static_cast<void>(ImmersedDomain::create(
              surface, query, side, topology, geometry, boundaries, mpi));
        },
        mpi, marker, 0);
    HUNDUN_CHECK(surface.fingerprint() == surface_before);
    HUNDUN_CHECK(query.fingerprint() == query_before);
    HUNDUN_CHECK(topology.local_cell_count() == cells_before);
    HUNDUN_CHECK(topology.local_face_count() == faces_before);
    HUNDUN_CHECK(boundaries.open_domain() == open_before);
    for (std::uint32_t patch = 0U; patch < 6U; ++patch) {
      HUNDUN_CHECK(boundaries.patch(patch).kind() == kinds_before[patch]);
    }
  };

  {
    auto decomposition = StructuredDecomposition::create(
        mpi, kExtent, periodic_axes(walls), DecompositionOptions{grid});
    MeshTopology topology(decomposition);
    MeshGeometry geometry(topology,
                          hundun::mesh::UniformBoxMapping{kOrigin, kLength});
    const auto &local_kinds = mpi.rank() == 1 ? open : walls;
    BoundaryRegistry boundaries = BoundaryRegistry::create(
        make_config(mpi.size(), grid, local_kinds), topology);
    const ImmersedSurface surface =
        ImmersedSurface::load_collective(cube, 1.0, mpi, 0);
    const SurfaceQuery query = SurfaceQuery::create(surface);
    const std::uint64_t surface_before = surface.fingerprint();
    const std::uint64_t query_before = query.fingerprint();
    const bool open_before = boundaries.open_domain();
    expect_collective_error(
        [&] {
          static_cast<void>(
              ImmersedDomain::create(surface, query, ImmersedFluidSide::outside,
                                     topology, geometry, boundaries, mpi));
        },
        mpi, "boundary_registry_consistent", 1);
    HUNDUN_CHECK(surface.fingerprint() == surface_before);
    HUNDUN_CHECK(query.fingerprint() == query_before);
    HUNDUN_CHECK(boundaries.open_domain() == open_before);
  }

  run_failure(cube, open, ImmersedFluidSide::inside,
              "active_velocity_inlet_nonzero");
  run_failure(near_x_min, open, ImmersedFluidSide::inside,
              "active_pressure_outlet_nonzero");
  run_failure(near_x_min, x_periodic, ImmersedFluidSide::outside,
              "active_periodic_pairing");
  run_failure(x_min_slab, x_periodic, ImmersedFluidSide::outside,
              "active_periodic_pairing");
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = MpiContext::duplicate(MPI_COMM_WORLD);
  SharedDirectory directory(mpi);
  const std::string mode = argc > 1 ? argv[1] : "success";
  if (mode == "success") {
    success_cases(mpi, directory,
                  argc > 2 && std::string(argv[2]) == "alternate");
  } else if (mode == "failures") {
    failure_cases(mpi, directory);
  } else {
    HUNDUN_CHECK(false);
  }
  return 0;
}
