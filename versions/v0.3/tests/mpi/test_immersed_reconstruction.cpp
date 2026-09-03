// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/fvm_immersed_reconstruction.hpp"

#include "src/fvm_immersed_boundary_authority_detail.hpp"
#include "tests/support/fvm_immersed_reconstruction_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
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
#include <vector>

namespace {

using namespace hundun;
using finite_volume::test::ImmersedReconstructionTestAccess;
using finite_volume::test::ReconstructionExecutionStage;

constexpr runtime::Int3 kExtent{12, 12, 12};
constexpr runtime::PhaseId kPhase = 71U;
constexpr runtime::ActorId kActor = 19U;

enum class BoundaryMode : std::uint8_t {
  no_slip_wall,
  periodic,
  symmetry,
  open
};

std::uint64_t bits(double value) {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

runtime::Real3 add(runtime::Real3 a, runtime::Real3 b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

runtime::Real3 subtract(runtime::Real3 a, runtime::Real3 b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

runtime::Real3 multiply(runtime::Real3 a, double s) {
  return {a.x * s, a.y * s, a.z * s};
}

double dot(runtime::Real3 a, runtime::Real3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

runtime::Real3 cross(runtime::Real3 a, runtime::Real3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

runtime::Real3 polynomial_gradient(runtime::Real3 p) {
  return {2.0 + 1.4 * p.x - 0.2 * p.y + 0.4 * p.z,
          -3.0 - 0.2 * p.x + 0.6 * p.y - 0.6 * p.z,
          0.5 + 0.4 * p.x - 0.6 * p.y + 1.8 * p.z};
}

double independent_cell_average(runtime::Int3 cell,
                                const mesh::MeshGeometry &geometry) {
  constexpr std::array<runtime::Int3, 8> offsets{
      runtime::Int3{0, 0, 0}, runtime::Int3{1, 0, 0}, runtime::Int3{1, 1, 0},
      runtime::Int3{0, 1, 0}, runtime::Int3{0, 0, 1}, runtime::Int3{1, 0, 1},
      runtime::Int3{1, 1, 1}, runtime::Int3{0, 1, 1}};
  std::array<runtime::Real3, 8> vertices{};
  runtime::Real3 reference{};
  for (std::size_t i = 0U; i < vertices.size(); ++i) {
    vertices[i] = geometry.vertex_position_m(
        {cell.x + offsets[i].x, cell.y + offsets[i].y, cell.z + offsets[i].z});
    reference = add(reference, multiply(vertices[i], 0.125));
  }
  constexpr std::array<std::array<std::size_t, 3>, 12> triangles{{
      {{0, 4, 7}},
      {{0, 7, 3}},
      {{1, 2, 6}},
      {{1, 6, 5}},
      {{0, 1, 5}},
      {{0, 5, 4}},
      {{3, 7, 6}},
      {{3, 6, 2}},
      {{0, 3, 2}},
      {{0, 2, 1}},
      {{4, 5, 6}},
      {{4, 6, 7}},
  }};
  double integral = 0.0;
  double volume = 0.0;
  for (const auto triangle : triangles) {
    const std::array<runtime::Real3, 4> tetra{reference, vertices[triangle[0]],
                                              vertices[triangle[1]],
                                              vertices[triangle[2]]};
    const double tetra_volume =
        dot(subtract(tetra[1], tetra[0]),
            cross(subtract(tetra[2], tetra[0]), subtract(tetra[3], tetra[0]))) /
        6.0;
    HUNDUN_CHECK(tetra_volume > 0.0);
    // Exact analytic quadratic tetra average, independent of product moments.
    std::array<runtime::Real3, 4> p = tetra;
    runtime::Real3 mean{};
    for (const auto v : p)
      mean = add(mean, multiply(v, 0.25));
    double exact = 1.0 + 2.0 * mean.x - 3.0 * mean.y + 0.5 * mean.z;
    const auto product_mean = [&](int lhs, int rhs) {
      double sum_lhs = 0.0;
      double sum_rhs = 0.0;
      double diagonal = 0.0;
      for (const auto v : p) {
        const std::array<double, 3> q{v.x, v.y, v.z};
        sum_lhs += q[static_cast<std::size_t>(lhs)];
        sum_rhs += q[static_cast<std::size_t>(rhs)];
        diagonal +=
            q[static_cast<std::size_t>(lhs)] * q[static_cast<std::size_t>(rhs)];
      }
      return (sum_lhs * sum_rhs + diagonal) / 20.0;
    };
    exact += 0.7 * product_mean(0, 0) - 0.2 * product_mean(0, 1) +
             0.4 * product_mean(0, 2) + 0.3 * product_mean(1, 1) -
             0.6 * product_mean(1, 2) + 0.9 * product_mean(2, 2);
    integral += tetra_volume * exact;
    volume += tetra_volume;
  }
  return integral / volume;
}

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  return {2, 2, 1};
}

config::FlowCaseConfig case_config(int ranks, runtime::Int3 grid, bool warped,
                                   BoundaryMode mode) {
  config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "task7-immersed-reconstruction";
  config.simulation_type = config::SimulationType::variable_density_flow;
  config.density_model = config::DensityModel::constant;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = grid;
  config.mesh.cells = kExtent;
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.mesh.mapping = warped ? config::MeshMapping::analytic_warped_box
                               : config::MeshMapping::uniform_box;
  if (warped)
    config.mesh.warp_amplitude = runtime::Real3{0.02, -0.015, 0.01};
  config.time.mode = config::TimeMode::fixed;
  config.time.steps = 1;
  config.time.initial_dt_s = 0.01;
  config.time.min_dt_s = 0.01;
  config.time.max_dt_s = 0.01;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.dynamic_viscosity_pa_s = 1.0e-3;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back(config::FlowScalarConfig{"alpha", 0.01});
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t p = 0U; p < names.size(); ++p) {
    config::FlowBoundaryConfig boundary{};
    boundary.patch = names[p];
    boundary.type =
        mode == BoundaryMode::periodic   ? config::BoundaryType::periodic
        : mode == BoundaryMode::symmetry ? config::BoundaryType::symmetry
                                         : config::BoundaryType::no_slip_wall;
    config.boundaries[p] = boundary;
  }
  if (mode == BoundaryMode::open) {
    auto &inlet = config.boundaries[0];
    inlet.type = config::BoundaryType::velocity_inlet;
    inlet.velocity_m_per_s = runtime::Real3{1.0, 0.2, -0.1};
    inlet.thermal_authority = config::InletThermalAuthority::enthalpy;
    inlet.enthalpy_J_per_kg = 12.5;
    inlet.scalar_values =
        std::vector<config::InletScalarValue>{{"alpha", 0.375}};
    auto &outlet = config.boundaries[1];
    outlet.type = config::BoundaryType::pressure_outlet;
    outlet.pressure_perturbation_pa = 17.0;
    for (std::size_t p = 2U; p < config.boundaries.size(); ++p)
      config.boundaries[p].type = config::BoundaryType::symmetry;
  }
  return config;
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

class FixtureFile final {
public:
  explicit FixtureFile(const runtime::MpiContext &mpi) : mpi_(&mpi) {
    std::string text;
    if (mpi.rank() == 0) {
      path_ =
          std::filesystem::temp_directory_path() /
          ("hundun-task7-" +
           std::to_string(
               std::chrono::steady_clock::now().time_since_epoch().count()) +
           ".stl");
      test::write_text(
          path_, test::ascii_stl(
                     partition_conforming_cube({0.75, 0.75, 0.75}),
                     "task7-cube"));
      text = path_.string();
    }
    std::uint64_t size = text.size();
    HUNDUN_CHECK(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()) ==
                 MPI_SUCCESS);
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
  const std::filesystem::path &path() const noexcept { return path_; }

private:
  const runtime::MpiContext *mpi_;
  std::filesystem::path path_;
};

runtime::FieldDescriptor cell_descriptor(const char *name,
                                         std::uint32_t components, int ghost) {
  return {name,
          "1",
          "task7",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_descriptor(const char *name,
                                         std::uint32_t components) {
  return {name,
          "1",
          "task7",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

struct Fields final {
  runtime::FieldId scalar{};
  runtime::FieldId velocity{};
  runtime::FieldId scalar_gradient{};
  runtime::FieldId velocity_gradient{};
  runtime::FieldId scalar_faces{};
  runtime::FieldId velocity_faces{};
  runtime::FieldId mass_flux{};
};

Fields declare_fields(runtime::FieldRegistry &registry, int ghost) {
  Fields result{};
  result.scalar = registry.declare_field(cell_descriptor("q", 1U, ghost));
  result.velocity = registry.declare_field(cell_descriptor("u", 3U, ghost));
  result.scalar_gradient =
      registry.declare_field(cell_descriptor("grad_q", 3U, 0));
  result.velocity_gradient =
      registry.declare_field(cell_descriptor("grad_u", 9U, 0));
  result.scalar_faces = registry.declare_field(face_descriptor("q_face", 1U));
  result.velocity_faces = registry.declare_field(face_descriptor("u_face", 3U));
  result.mass_flux = finite_volume::declare_face_mass_flux(registry);
  return result;
}

finite_volume::ReconstructionFieldBinding
binding(runtime::FieldStorage &storage, const runtime::FieldAccessPlan &access,
        runtime::FieldId field) {
  return {storage, access, kPhase, kActor, field};
}

void check_trace() {
  const std::vector<ReconstructionExecutionStage> expected{
      ReconstructionExecutionStage::halo_begin,
      ReconstructionExecutionStage::active_interior,
      ReconstructionExecutionStage::halo_wait,
      ReconstructionExecutionStage::remote_ghost_symbols,
      ReconstructionExecutionStage::partition_boundary,
      ReconstructionExecutionStage::immersed_interface};
  HUNDUN_CHECK(ImmersedReconstructionTestAccess::trace() == expected);
  HUNDUN_CHECK(ImmersedReconstructionTestAccess::inactive_read_attempts() ==
               0U);
}

std::size_t local_owned_index(const mesh::MeshTopology &topology,
                              mesh::GlobalCellId global) {
  const auto local = topology.find_local_cell(global);
  HUNDUN_CHECK(local.has_value());
  HUNDUN_CHECK(*local < topology.owned_cell_count());
  return *local;
}

runtime::Int3 wrap_cell(runtime::Int3 cell, bool periodic) {
  if (!periodic)
    return cell;
  int *coordinates[3]{&cell.x, &cell.y, &cell.z};
  const int extents[3]{kExtent.x, kExtent.y, kExtent.z};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    *coordinates[axis] %= extents[axis];
    if (*coordinates[axis] < 0)
      *coordinates[axis] += extents[axis];
  }
  return cell;
}

struct CanonicalFaceLine final {
  mesh::LocalCellId p{};
  mesh::LocalCellId n{};
  runtime::Int3 pm1{};
  runtime::Int3 np1{};
  double mass_flux{};
};

struct LinkPair final {
  std::uint64_t id{};
  mesh::GlobalCellId fluid{};
  mesh::GlobalCellId solid{};
};

std::vector<LinkPair> gather_link_pairs(const immersed::ImmersedDomain &domain,
                                        const runtime::MpiContext &mpi) {
  std::vector<std::uint64_t> local;
  local.reserve(domain.links().size() * 3U);
  for (const auto &link : domain.links()) {
    local.push_back(link.id);
    local.push_back(link.fluid_cell);
    local.push_back(link.solid_cell);
  }
  const int local_count = static_cast<int>(local.size());
  std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
  HUNDUN_CHECK(MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1,
                             MPI_INT, mpi.comm()) == MPI_SUCCESS);
  std::vector<int> displacements(counts.size(), 0);
  int total = 0;
  for (std::size_t rank = 0U; rank < counts.size(); ++rank) {
    displacements[rank] = total;
    total += counts[rank];
  }
  std::vector<std::uint64_t> gathered(static_cast<std::size_t>(total));
  HUNDUN_CHECK(MPI_Allgatherv(local.data(), local_count, MPI_UINT64_T,
                              gathered.data(), counts.data(),
                              displacements.data(), MPI_UINT64_T,
                              mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(gathered.size() % 3U == 0U);
  std::vector<LinkPair> result;
  result.reserve(gathered.size() / 3U);
  for (std::size_t index = 0U; index < gathered.size(); index += 3U)
    result.push_back(
        {gathered[index], gathered[index + 1U], gathered[index + 2U]});
  std::sort(result.begin(), result.end(),
            [](const auto &lhs, const auto &rhs) { return lhs.id < rhs.id; });
  return result;
}

CanonicalFaceLine canonical_face_line(const mesh::MeshTopology &topology,
                                      mesh::LocalFaceId face,
                                      double local_mass_flux) {
  const auto neighbour = topology.neighbour(face);
  HUNDUN_CHECK(neighbour.has_value());
  const auto periodic_pair = topology.periodic_pair(face);
  const bool reversed = periodic_pair.has_value() &&
                        topology.global_face_id(face) > *periodic_pair;
  const auto p = reversed ? *neighbour : topology.owner(face);
  const auto n = reversed ? topology.owner(face) : *neighbour;
  const auto p_global = topology.global_cell(p);
  const auto n_global = topology.global_cell(n);
  runtime::Int3 delta{n_global.x - p_global.x, n_global.y - p_global.y,
                      n_global.z - p_global.z};
  if (periodic_pair.has_value()) {
    if (std::abs(delta.x) > 1)
      delta.x = delta.x > 0 ? -1 : 1;
    if (std::abs(delta.y) > 1)
      delta.y = delta.y > 0 ? -1 : 1;
    if (std::abs(delta.z) > 1)
      delta.z = delta.z > 0 ? -1 : 1;
  }
  return {p,
          n,
          {p_global.x - delta.x, p_global.y - delta.y, p_global.z - delta.z},
          {n_global.x + delta.x, n_global.y + delta.y, n_global.z + delta.z},
          reversed ? -local_mass_flux : local_mass_flux};
}

double stencil_oracle(const mesh::MeshTopology &topology,
                      const mesh::MeshGeometry &geometry,
                      const immersed::GhostStencilPlan &ghost_plan,
                      const std::vector<LinkPair> &links,
                      mesh::GlobalCellId fluid, runtime::Int3 cell,
                      bool periodic, double fallback, double scale,
                      bool velocity, std::size_t component) {
  const bool inside = cell.x >= 0 && cell.y >= 0 && cell.z >= 0 &&
                      cell.x < kExtent.x && cell.y < kExtent.y &&
                      cell.z < kExtent.z;
  if (!inside && !periodic)
    return fallback;
  cell = wrap_cell(cell, periodic);
  const auto target = topology.global_cell_id(cell);
  const auto link =
      std::find_if(links.begin(), links.end(), [&](const auto &candidate) {
        return candidate.fluid == fluid && candidate.solid == target;
      });
  if (link == links.end())
    return scale * independent_cell_average(cell, geometry);
  const auto &donors =
      velocity ? ghost_plan.velocity_constraint(link->id, component).donors
               : ghost_plan.zero_normal_constraint(link->id).donors;
  double result = 0.0;
  for (const auto &donor : donors) {
    const auto nx = static_cast<std::uint64_t>(kExtent.x);
    const auto ny = static_cast<std::uint64_t>(kExtent.y);
    const auto plane = nx * ny;
    const runtime::Int3 logical{static_cast<int>(donor.global_cell % nx),
                                static_cast<int>((donor.global_cell / nx) % ny),
                                static_cast<int>(donor.global_cell / plane)};
    result +=
        donor.weight * scale * independent_cell_average(logical, geometry);
  }
  return result;
}

double mc_oracle(const mesh::MeshTopology &topology,
                 const mesh::MeshGeometry &geometry,
                 const immersed::GhostStencilPlan &ghost_plan,
                 const std::vector<LinkPair> &links, mesh::LocalFaceId face,
                 double local_mass_flux, bool periodic, double scale,
                 bool velocity = false, std::size_t component = 0U) {
  const auto line = canonical_face_line(topology, face, local_mass_flux);
  const double qp =
      scale * independent_cell_average(topology.global_cell(line.p), geometry);
  const double qn =
      scale * independent_cell_average(topology.global_cell(line.n), geometry);
  const double qm = stencil_oracle(topology, geometry, ghost_plan, links,
                                   topology.global_cell_id(line.p), line.pm1,
                                   periodic, qp, scale, velocity, component);
  const double qplus = stencil_oracle(topology, geometry, ghost_plan, links,
                                      topology.global_cell_id(line.n), line.np1,
                                      periodic, qn, scale, velocity, component);
  const double slope_p = finite_volume::monotonized_central(qp - qm, qn - qp);
  const double slope_n =
      finite_volume::monotonized_central(qn - qp, qplus - qn);
  return line.mass_flux >= 0.0 ? qp + 0.5 * slope_p : qn - 0.5 * slope_n;
}

void run_success(const runtime::MpiContext &mpi, bool warped,
                 BoundaryMode boundary_mode) {
  const auto grid = process_grid(mpi.size());
  const bool periodic = boundary_mode == BoundaryMode::periodic;
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, kExtent, {periodic, periodic, periodic},
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
      case_config(mpi.size(), grid, warped, boundary_mode), topology);
  FixtureFile fixture(mpi);
  const auto surface =
      immersed::ImmersedSurface::load_collective(fixture.path(), 0.4, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  const auto ghost_plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  const auto global_links = gather_link_pairs(domain, mpi);
  const int ghost = static_cast<int>(ghost_plan.maximum_halo_reach());
  HUNDUN_CHECK(ghost >= 1);
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), ghost));
  auto provider = finite_volume::ImmersedReconstruction::create(
      topology, geometry, boundaries, domain, ghost_plan, decomposition, mpi,
      halo);

  runtime::FieldRegistry registry;
  const Fields ids = declare_fields(registry, ghost);
  registry.freeze();
  runtime::FieldStorage storage(
      registry, runtime::FieldLayoutSet{decomposition.local_extent(),
                                        topology.local_face_count()});
  runtime::FieldAccessPlan access(registry);
  access.declare_access(kPhase, kActor, ids.scalar, runtime::AccessMode::read);
  access.declare_access(kPhase, kActor, ids.velocity,
                        runtime::AccessMode::read);
  access.declare_access(kPhase, kActor, ids.scalar_gradient,
                        runtime::AccessMode::read_write);
  access.declare_access(kPhase, kActor, ids.velocity_gradient,
                        runtime::AccessMode::read_write);
  access.declare_access(kPhase, kActor, ids.scalar_faces,
                        runtime::AccessMode::read_write);
  access.declare_access(kPhase, kActor, ids.velocity_faces,
                        runtime::AccessMode::read_write);
  access.declare_access(kPhase, kActor, ids.mass_flux,
                        runtime::AccessMode::read);
  access.freeze();

  auto scalar = storage.view<double>(ids.scalar);
  auto velocity = storage.view<double>(ids.velocity);
  const auto box = decomposition.owned_box();
  const auto local = decomposition.local_extent();
  for (int k = -ghost; k < local.z + ghost; ++k) {
    for (int j = -ghost; j < local.y + ghost; ++j) {
      for (int i = -ghost; i < local.x + ghost; ++i) {
        const runtime::Int3 global{box.begin.x + i, box.begin.y + j,
                                   box.begin.z + k};
        double value = 0.0;
        if (global.x >= 0 && global.y >= 0 && global.z >= 0 &&
            global.x < kExtent.x && global.y < kExtent.y &&
            global.z < kExtent.z) {
          value = independent_cell_average(global, geometry);
        }
        scalar(i, j, k, 0) = value;
        velocity(i, j, k, 0) = value;
        velocity(i, j, k, 1) = 2.0 * value;
        velocity(i, j, k, 2) = -0.5 * value;
      }
    }
  }
  for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::solid)
      continue;
    const auto global = topology.global_cell(cell);
    const int i = global.x - box.begin.x;
    const int j = global.y - box.begin.y;
    const int k = global.z - box.begin.z;
    scalar(i, j, k, 0) = 12345.0;
    velocity(i, j, k, 0) = -23456.0;
    velocity(i, j, k, 1) = 34567.0;
    velocity(i, j, k, 2) = -45678.0;
  }
  auto face_flux = storage.acquire_face_write<double>(
      [&] {
        runtime::FieldAccessPlan plan(registry);
        plan.declare_access(kPhase + 1U, kActor, ids.mass_flux,
                            runtime::AccessMode::write);
        plan.freeze();
        return plan;
      }(),
      kPhase + 1U, kActor, ids.mass_flux);
  std::array<int, 3> local_momentum_policy_counts{};
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    switch (topology.global_face_id(face) % 3U) {
    case 0U:
      face_flux(face, 0) = 0.75;
      break;
    case 1U:
      face_flux(face, 0) = -0.5;
      break;
    default:
      face_flux(face, 0) = 0.0;
      break;
    }
  }

  auto mass_flux = finite_volume::FaceMassFlux::acquire(
      registry, storage, access, kPhase, kActor, ids.mass_flux, topology);
  const auto snapshot = ImmersedReconstructionTestAccess::snapshot(provider);
  HUNDUN_CHECK(snapshot.owned_active_count ==
               domain.active_cells().owned_active_count());
  HUNDUN_CHECK(snapshot.interior_rows.size() + snapshot.partition_rows.size() +
                   snapshot.interface_rows.size() ==
               snapshot.owned_active_count);
  HUNDUN_CHECK(!snapshot.interface_rows.empty());
  HUNDUN_CHECK(snapshot.fingerprint == provider.dependency_fingerprint());
  HUNDUN_CHECK(std::is_sorted(snapshot.interior_rows.begin(),
                              snapshot.interior_rows.end()));
  HUNDUN_CHECK(std::is_sorted(snapshot.partition_rows.begin(),
                              snapshot.partition_rows.end()));
  HUNDUN_CHECK(std::is_sorted(snapshot.interface_rows.begin(),
                              snapshot.interface_rows.end()));
  std::vector<mesh::GlobalCellId> planned_rows = snapshot.interior_rows;
  planned_rows.insert(planned_rows.end(), snapshot.partition_rows.begin(),
                      snapshot.partition_rows.end());
  planned_rows.insert(planned_rows.end(), snapshot.interface_rows.begin(),
                      snapshot.interface_rows.end());
  std::sort(planned_rows.begin(), planned_rows.end());
  HUNDUN_CHECK(std::adjacent_find(planned_rows.begin(), planned_rows.end()) ==
               planned_rows.end());
  for (const auto global_id : planned_rows) {
    const auto cell = local_owned_index(topology, global_id);
    HUNDUN_CHECK(domain.region(cell) == immersed::CellRegion::fluid);
  }
  std::uint64_t fingerprint_min = 0U;
  std::uint64_t fingerprint_max = 0U;
  const std::uint64_t fingerprint = provider.dependency_fingerprint();
  HUNDUN_CHECK(MPI_Allreduce(&fingerprint, &fingerprint_min, 1, MPI_UINT64_T,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&fingerprint, &fingerprint_max, 1, MPI_UINT64_T,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(fingerprint_min == fingerprint_max);

  ImmersedReconstructionTestAccess::reset_trace();
  ImmersedReconstructionTestAccess::fail_on_inactive_read(true);
  provider.compute_gradient(
      finite_volume::GradientScheme::weighted_least_squares,
      finite_volume::FiniteVolumeQuantity::pressure(),
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.scalar_gradient));
  check_trace();
  const auto gradient =
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.scalar_gradient);
  const double tolerance = 32768.0 * std::numeric_limits<double>::epsilon();
  for (const auto global_id : snapshot.interface_rows) {
    const auto cell = local_owned_index(topology, global_id);
    const auto logical = topology.global_cell(cell);
    const auto expected = polynomial_gradient(geometry.cell_center_m(cell));
    HUNDUN_CHECK_NEAR(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                               logical.z - box.begin.z, 0),
                      expected.x,
                      tolerance * std::max(1.0, std::abs(expected.x)));
    HUNDUN_CHECK_NEAR(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                               logical.z - box.begin.z, 1),
                      expected.y,
                      tolerance * std::max(1.0, std::abs(expected.y)));
    HUNDUN_CHECK_NEAR(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                               logical.z - box.begin.z, 2),
                      expected.z,
                      tolerance * std::max(1.0, std::abs(expected.z)));
  }
  std::vector<finite_volume::detail::ImmersedWallNormalGradient>
      wall_normal_conditions;
  wall_normal_conditions.reserve(global_links.size());
  for (const auto &link : global_links)
    wall_normal_conditions.push_back({link.id, 0.75});
  finite_volume::detail::compute_pressure_gradient_with_wall_normal_constraints(
      provider, finite_volume::GradientScheme::weighted_least_squares,
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.scalar_gradient), wall_normal_conditions);
  std::vector<std::uint64_t> constrained_gradient_bits;
  for (const auto global_id : snapshot.interface_rows) {
    const auto cell = local_owned_index(topology, global_id);
    const auto logical = topology.global_cell(cell);
    for (int component = 0; component < 3; ++component)
      constrained_gradient_bits.push_back(
          bits(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                        logical.z - box.begin.z, component)));
  }
  for (auto &condition : wall_normal_conditions)
    condition.value = 0.0;
  finite_volume::detail::compute_pressure_gradient_with_wall_normal_constraints(
      provider, finite_volume::GradientScheme::weighted_least_squares,
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.scalar_gradient), wall_normal_conditions);
  std::size_t constrained_apply_allocations = 0U;
  {
    test::allocation_probe::AllocationAttemptGuard allocation_guard;
    finite_volume::detail::
        compute_pressure_gradient_with_wall_normal_constraints(
            provider, finite_volume::GradientScheme::weighted_least_squares,
            binding(storage, access, ids.scalar),
            binding(storage, access, ids.scalar_gradient),
            wall_normal_conditions);
    constrained_apply_allocations = allocation_guard.attempts();
  }
  HUNDUN_CHECK(constrained_apply_allocations == 0U);
  bool every_interface_gradient_changed = true;
  std::size_t constrained_index = 0U;
  for (const auto global_id : snapshot.interface_rows) {
    const auto cell = local_owned_index(topology, global_id);
    const auto logical = topology.global_cell(cell);
    bool row_changed = false;
    for (int component = 0; component < 3; ++component)
      row_changed =
          row_changed ||
          constrained_gradient_bits[constrained_index++] !=
              bits(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                            logical.z - box.begin.z, component));
    every_interface_gradient_changed =
        every_interface_gradient_changed && row_changed;
  }
  HUNDUN_CHECK(every_interface_gradient_changed);
  {
    auto missing = wall_normal_conditions;
    missing.clear();
    std::string message;
    try {
      finite_volume::detail::
          compute_pressure_gradient_with_wall_normal_constraints(
              provider, finite_volume::GradientScheme::weighted_least_squares,
              binding(storage, access, ids.scalar),
              binding(storage, access, ids.scalar_gradient), missing);
    } catch (const runtime::Error &error) {
      message = error.what();
    }
    HUNDUN_CHECK(message.find("condition is missing") != std::string::npos);
  }
  {
    auto duplicated = wall_normal_conditions;
    duplicated.push_back(duplicated.front());
    std::string message;
    try {
      finite_volume::detail::
          compute_pressure_gradient_with_wall_normal_constraints(
              provider, finite_volume::GradientScheme::weighted_least_squares,
              binding(storage, access, ids.scalar),
              binding(storage, access, ids.scalar_gradient), duplicated);
    } catch (const runtime::Error &error) {
      message = error.what();
    }
    HUNDUN_CHECK(message.find("condition is duplicated") != std::string::npos);
  }
  provider.compute_gradient(
      finite_volume::GradientScheme::weighted_least_squares,
      finite_volume::FiniteVolumeQuantity::density(),
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.scalar_gradient));

  ImmersedReconstructionTestAccess::reset_trace();
  provider.compute_gradient(finite_volume::GradientScheme::green_gauss,
                            finite_volume::FiniteVolumeQuantity::velocity(),
                            binding(storage, access, ids.velocity),
                            binding(storage, access, ids.velocity_gradient));
  check_trace();
  const auto velocity_gradient =
      static_cast<const runtime::FieldStorage &>(storage).view<double>(
          ids.velocity_gradient);
  int local_constrained_difference = 0;
  for (const auto global_id : snapshot.interface_rows) {
    const auto cell = local_owned_index(topology, global_id);
    const auto logical = topology.global_cell(cell);
    const auto expected = polynomial_gradient(geometry.cell_center_m(cell));
    const int i = logical.x - box.begin.x;
    const int j = logical.y - box.begin.y;
    const int k = logical.z - box.begin.z;
    const double expected_values[9]{
        expected.x,        expected.y,        expected.z,
        2.0 * expected.x,  2.0 * expected.y,  2.0 * expected.z,
        -0.5 * expected.x, -0.5 * expected.y, -0.5 * expected.z};
    for (int component = 0; component < 9; ++component) {
      const double observed = velocity_gradient(i, j, k, component);
      HUNDUN_CHECK(std::isfinite(observed));
      if (std::abs(observed - expected_values[component]) >
          2.0 * tolerance * std::max(1.0, std::abs(expected_values[component])))
        local_constrained_difference = 1;
    }
  }
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &local_constrained_difference, 1,
                             MPI_INT, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(local_constrained_difference == 1);

  for (int k = -ghost; k < local.z + ghost; ++k)
    for (int j = -ghost; j < local.y + ghost; ++j)
      for (int i = -ghost; i < local.x + ghost; ++i) {
        velocity(i, j, k, 0) = 1.0;
        velocity(i, j, k, 1) = -2.0;
        velocity(i, j, k, 2) = 0.5;
      }
  provider.compute_gradient(
      finite_volume::GradientScheme::weighted_least_squares,
      finite_volume::FiniteVolumeQuantity::velocity(),
      binding(storage, access, ids.velocity),
      binding(storage, access, ids.velocity_gradient));
  double local_restoring_gradient = 0.0;
  for (const auto global_id : snapshot.interface_rows) {
    const auto cell = local_owned_index(topology, global_id);
    const auto logical = topology.global_cell(cell);
    const int i = logical.x - box.begin.x;
    const int j = logical.y - box.begin.y;
    const int k = logical.z - box.begin.z;
    for (int component = 0; component < 9; ++component)
      local_restoring_gradient =
          std::max(local_restoring_gradient,
                   std::abs(velocity_gradient(i, j, k, component)));
  }
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &local_restoring_gradient, 1,
                             MPI_DOUBLE, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(local_restoring_gradient > 1.0e-6);

  for (int k = -ghost; k < local.z + ghost; ++k)
    for (int j = -ghost; j < local.y + ghost; ++j)
      for (int i = -ghost; i < local.x + ghost; ++i) {
        const runtime::Int3 global{box.begin.x + i, box.begin.y + j,
                                   box.begin.z + k};
        double value = 0.0;
        if (global.x >= 0 && global.y >= 0 && global.z >= 0 &&
            global.x < kExtent.x && global.y < kExtent.y &&
            global.z < kExtent.z)
          value = independent_cell_average(global, geometry);
        velocity(i, j, k, 0) = value;
        velocity(i, j, k, 1) = 2.0 * value;
        velocity(i, j, k, 2) = -0.5 * value;
      }
  for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::solid)
      continue;
    const auto global = topology.global_cell(cell);
    velocity(global.x - box.begin.x, global.y - box.begin.y,
             global.z - box.begin.z, 0) = -23456.0;
    velocity(global.x - box.begin.x, global.y - box.begin.y,
             global.z - box.begin.z, 1) = 34567.0;
    velocity(global.x - box.begin.x, global.y - box.begin.y,
             global.z - box.begin.z, 2) = -45678.0;
  }

  provider.compute_gradient(
      finite_volume::GradientScheme::weighted_least_squares,
      finite_volume::FiniteVolumeQuantity::scalar(0U),
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.scalar_gradient));

  std::vector<std::uint64_t> first_active_bits;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const auto logical = topology.global_cell(cell);
    for (int component = 0; component < 3; ++component)
      first_active_bits.push_back(
          bits(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                        logical.z - box.begin.z, component)));
  }
  for (mesh::LocalCellId cell = 0U; cell < topology.local_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::solid)
      continue;
    const auto global = topology.global_cell(cell);
    scalar(global.x - box.begin.x, global.y - box.begin.y,
           global.z - box.begin.z, 0) = -98765.0;
  }
  provider.compute_gradient(
      finite_volume::GradientScheme::weighted_least_squares,
      finite_volume::FiniteVolumeQuantity::scalar(0U),
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.scalar_gradient));
  std::size_t bit_index = 0U;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const auto logical = topology.global_cell(cell);
    for (int component = 0; component < 3; ++component)
      HUNDUN_CHECK(
          first_active_bits[bit_index++] ==
          bits(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                        logical.z - box.begin.z, component)));
  }

  const auto poison_ghost_cells = [&](auto &field, int components,
                                      double value) {
    if (mpi.size() <= 1)
      return;
    for (int k = -ghost; k < local.z + ghost; ++k)
      for (int j = -ghost; j < local.y + ghost; ++j)
        for (int i = -ghost; i < local.x + ghost; ++i)
          if (i < 0 || i >= local.x || j < 0 || j >= local.y || k < 0 ||
              k >= local.z)
            for (int component = 0; component < components; ++component)
              field(i, j, k, component) = value;
  };
  poison_ghost_cells(scalar, 1, -77777.0);
  ImmersedReconstructionTestAccess::reset_trace();
  provider.reconstruct_transport_faces(
      finite_volume::FiniteVolumeQuantity::enthalpy(), mass_flux,
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.scalar_faces));
  check_trace();
  const auto scalar_faces =
      static_cast<const runtime::FieldStorage &>(storage)
          .acquire_face_read<double>(access, kPhase, kActor, ids.scalar_faces);
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto neighbour = topology.neighbour(face);
    if (neighbour.has_value() &&
        domain.region(topology.owner(face)) != domain.region(*neighbour))
      HUNDUN_CHECK(bits(scalar_faces(face, 0)) == bits(0.0));
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto patch = topology.patch_id(face);
    if (!patch.has_value() || topology.periodic_pair(face).has_value() ||
        domain.region(topology.owner(face)) != immersed::CellRegion::fluid)
      continue;
    const auto logical = topology.global_cell(topology.owner(face));
    const double interior =
        scalar(logical.x - box.begin.x, logical.y - box.begin.y,
               logical.z - box.begin.z, 0);
    const double expected = boundaries.evaluate_enthalpy(*patch, interior).face;
    HUNDUN_CHECK(bits(scalar_faces(face, 0)) == bits(expected));
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    if (!neighbour.has_value() ||
        domain.region(owner) != immersed::CellRegion::fluid ||
        domain.region(*neighbour) != immersed::CellRegion::fluid)
      continue;
    const double expected =
        mc_oracle(topology, geometry, ghost_plan, global_links, face,
                  face_flux(face, 0), periodic, 1.0);
    HUNDUN_CHECK_NEAR(scalar_faces(face, 0), expected,
                      tolerance * std::max(1.0, std::abs(expected)));
  }
  std::array<int, 3> local_flux_sign_counts{};
  for (const auto global_face : snapshot.interface_faces) {
    const auto face = topology.find_local_face(global_face);
    if (!face.has_value())
      continue;
    const auto neighbour = topology.neighbour(*face);
    if (!neighbour.has_value() ||
        domain.region(topology.owner(*face)) != immersed::CellRegion::fluid ||
        domain.region(*neighbour) != immersed::CellRegion::fluid)
      continue;
    const double canonical_flux =
        canonical_face_line(topology, *face, face_flux(*face, 0)).mass_flux;
    ++local_flux_sign_counts[canonical_flux > 0.0   ? 0U
                             : canonical_flux < 0.0 ? 1U
                                                    : 2U];
  }
  std::array<int, 3> global_flux_sign_counts{};
  HUNDUN_CHECK(MPI_Allreduce(local_flux_sign_counts.data(),
                             global_flux_sign_counts.data(),
                             static_cast<int>(global_flux_sign_counts.size()),
                             MPI_INT, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  for (const int count : global_flux_sign_counts)
    HUNDUN_CHECK(count > 0);
  provider.reconstruct_transport_faces(
      finite_volume::FiniteVolumeQuantity::density(), mass_flux,
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.scalar_faces));
  provider.reconstruct_transport_faces(
      finite_volume::FiniteVolumeQuantity::scalar(0U), mass_flux,
      binding(storage, access, ids.scalar),
      binding(storage, access, ids.scalar_faces));
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto patch = topology.patch_id(face);
    if (!patch.has_value() || topology.periodic_pair(face).has_value() ||
        domain.region(topology.owner(face)) != immersed::CellRegion::fluid)
      continue;
    const auto logical = topology.global_cell(topology.owner(face));
    const double interior =
        scalar(logical.x - box.begin.x, logical.y - box.begin.y,
               logical.z - box.begin.z, 0);
    const double expected =
        boundaries.evaluate_scalar(*patch, 0U, interior).face;
    HUNDUN_CHECK(bits(scalar_faces(face, 0)) == bits(expected));
  }

  poison_ghost_cells(velocity, 3, -88888.0);
  ImmersedReconstructionTestAccess::reset_trace();
  provider.reconstruct_momentum_faces(
      mass_flux, binding(storage, access, ids.velocity),
      binding(storage, access, ids.velocity_faces));
  check_trace();
  const auto velocity_faces =
      static_cast<const runtime::FieldStorage &>(storage)
          .acquire_face_read<double>(access, kPhase, kActor,
                                     ids.velocity_faces);
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto neighbour = topology.neighbour(face);
    if (!neighbour.has_value() ||
        domain.region(topology.owner(face)) == domain.region(*neighbour))
      continue;
    HUNDUN_CHECK(bits(velocity_faces(face, 0)) == bits(0.0));
    HUNDUN_CHECK(bits(velocity_faces(face, 1)) == bits(0.0));
    HUNDUN_CHECK(bits(velocity_faces(face, 2)) == bits(0.0));
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto patch = topology.patch_id(face);
    if (!patch.has_value() || topology.periodic_pair(face).has_value() ||
        domain.region(topology.owner(face)) != immersed::CellRegion::fluid)
      continue;
    const auto logical = topology.global_cell(topology.owner(face));
    const runtime::Real3 interior{
        velocity(logical.x - box.begin.x, logical.y - box.begin.y,
                 logical.z - box.begin.z, 0),
        velocity(logical.x - box.begin.x, logical.y - box.begin.y,
                 logical.z - box.begin.z, 1),
        velocity(logical.x - box.begin.x, logical.y - box.begin.y,
                 logical.z - box.begin.z, 2)};
    const auto expected = boundaries.evaluate_velocity(
        *patch, interior,
        geometry.face_area_vector_m2(face, mesh::FaceSide::owner));
    HUNDUN_CHECK(bits(velocity_faces(face, 0)) == bits(expected.face.x));
    HUNDUN_CHECK(bits(velocity_faces(face, 1)) == bits(expected.face.y));
    HUNDUN_CHECK(bits(velocity_faces(face, 2)) == bits(expected.face.z));
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto neighbour = topology.neighbour(face);
    if (!neighbour.has_value() ||
        domain.region(topology.owner(face)) != immersed::CellRegion::fluid ||
        domain.region(*neighbour) != immersed::CellRegion::fluid)
      continue;
    const auto line = canonical_face_line(topology, face, face_flux(face, 0));
    const bool eligible_for_centered =
        geometry.face_skewness(face) <= 0.25 &&
        geometry.face_non_orthogonality_degrees(face) <= 70.0;
    const bool in_interface_band = std::find(snapshot.interface_faces.begin(),
                                             snapshot.interface_faces.end(),
                                             topology.global_face_id(face)) !=
                                   snapshot.interface_faces.end();
    const auto centered_component = [&](double scale) {
      return 0.5 * scale *
             (independent_cell_average(topology.global_cell(line.p), geometry) +
              independent_cell_average(topology.global_cell(line.n), geometry));
    };
    const auto expected_component = [&](double scale) {
      if (eligible_for_centered && !in_interface_band)
        return centered_component(scale);
      const std::size_t component = scale == 1.0 ? 0U : scale == 2.0 ? 1U : 2U;
      return mc_oracle(topology, geometry, ghost_plan, global_links, face,
                       face_flux(face, 0), periodic, scale, true, component);
    };
    const double expected_x = expected_component(1.0);
    const double expected_y = expected_component(2.0);
    const double expected_z = expected_component(-0.5);
    HUNDUN_CHECK_NEAR(velocity_faces(face, 0), expected_x,
                      tolerance * std::max(1.0, std::abs(expected_x)));
    HUNDUN_CHECK_NEAR(velocity_faces(face, 1), expected_y,
                      tolerance * std::max(1.0, std::abs(expected_y)));
    HUNDUN_CHECK_NEAR(velocity_faces(face, 2), expected_z,
                      tolerance * std::max(1.0, std::abs(expected_z)));
    if (eligible_for_centered && in_interface_band) {
      ++local_momentum_policy_counts[0];
      if (std::abs(expected_x - centered_component(1.0)) > 32.0 * tolerance)
        ++local_momentum_policy_counts[1];
    } else if (eligible_for_centered) {
      ++local_momentum_policy_counts[2];
    }
  }
  std::array<int, 3> global_momentum_policy_counts{};
  HUNDUN_CHECK(
      MPI_Allreduce(local_momentum_policy_counts.data(),
                    global_momentum_policy_counts.data(),
                    static_cast<int>(global_momentum_policy_counts.size()),
                    MPI_INT, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_momentum_policy_counts[0] > 0);
  HUNDUN_CHECK(global_momentum_policy_counts[1] > 0);
  HUNDUN_CHECK(global_momentum_policy_counts[2] > 0);

  poison_ghost_cells(velocity, 3, -99999.0);
  ImmersedReconstructionTestAccess::force_momentum_fallback(true);
  ImmersedReconstructionTestAccess::reset_trace();
  provider.reconstruct_momentum_faces(
      mass_flux, binding(storage, access, ids.velocity),
      binding(storage, access, ids.velocity_faces));
  check_trace();
  ImmersedReconstructionTestAccess::force_momentum_fallback(false);
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const auto neighbour = topology.neighbour(face);
    if (!neighbour.has_value() ||
        domain.region(topology.owner(face)) != immersed::CellRegion::fluid ||
        domain.region(*neighbour) != immersed::CellRegion::fluid)
      continue;
    constexpr std::array<double, 3> scales{1.0, 2.0, -0.5};
    for (std::size_t component = 0U; component < scales.size(); ++component) {
      const double expected = mc_oracle(
          topology, geometry, ghost_plan, global_links, face,
          face_flux(face, 0), periodic, scales[component], true, component);
      HUNDUN_CHECK_NEAR(velocity_faces(face, static_cast<int>(component)),
                        expected,
                        2.0 * tolerance * std::max(1.0, std::abs(expected)));
    }
  }

  std::vector<std::uint64_t> before_failure;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    const auto logical = topology.global_cell(cell);
    for (int component = 0; component < 3; ++component)
      before_failure.push_back(
          bits(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                        logical.z - box.begin.z, component)));
  }
  std::string alias_message;
  try {
    provider.compute_gradient(finite_volume::GradientScheme::green_gauss,
                              finite_volume::FiniteVolumeQuantity::density(),
                              binding(storage, access, ids.scalar),
                              binding(storage, access, ids.scalar));
  } catch (const runtime::Error &error) {
    alias_message = error.what();
  }
  HUNDUN_CHECK(alias_message.find("alias") != std::string::npos);
  std::size_t before_index = 0U;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    const auto logical = topology.global_cell(cell);
    for (int component = 0; component < 3; ++component)
      HUNDUN_CHECK(
          before_failure[before_index++] ==
          bits(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                        logical.z - box.begin.z, component)));
  }

  mesh::MeshTopology equivalent_but_distinct_topology(decomposition);
  auto distinct_mass_flux = finite_volume::FaceMassFlux::acquire(
      registry, storage, access, kPhase, kActor, ids.mass_flux,
      equivalent_but_distinct_topology);
  std::string topology_identity_message;
  try {
    provider.reconstruct_momentum_faces(
        distinct_mass_flux, binding(storage, access, ids.velocity),
        binding(storage, access, ids.velocity_faces));
  } catch (const runtime::Error &error) {
    topology_identity_message = error.what();
  }
  HUNDUN_CHECK(topology_identity_message.find("topology identity") !=
               std::string::npos);

  std::optional<std::pair<runtime::Int3, double>> poisoned_value;
  if (mpi.rank() == 0) {
    for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
         ++cell) {
      if (domain.region(cell) != immersed::CellRegion::fluid)
        continue;
      const auto logical = topology.global_cell(cell);
      const runtime::Int3 index{logical.x - box.begin.x,
                                logical.y - box.begin.y,
                                logical.z - box.begin.z};
      poisoned_value = std::pair{index, scalar(index.x, index.y, index.z, 0)};
      scalar(index.x, index.y, index.z, 0) =
          std::numeric_limits<double>::quiet_NaN();
      break;
    }
    HUNDUN_CHECK(poisoned_value.has_value());
  }
  std::string post_halo_message;
  try {
    provider.compute_gradient(
        finite_volume::GradientScheme::weighted_least_squares,
        finite_volume::FiniteVolumeQuantity::density(),
        binding(storage, access, ids.scalar),
        binding(storage, access, ids.scalar_gradient));
  } catch (const runtime::Error &error) {
    post_halo_message = error.what();
  }
  if (poisoned_value.has_value()) {
    const auto [index, value] = *poisoned_value;
    scalar(index.x, index.y, index.z, 0) = value;
  }
  HUNDUN_CHECK(post_halo_message.find("lowest failing rank 0") !=
               std::string::npos);
  std::uint64_t failure_hash = UINT64_C(14695981039346656037);
  for (const char character : post_halo_message) {
    failure_hash ^= static_cast<unsigned char>(character);
    failure_hash *= UINT64_C(1099511628211);
  }
  std::uint64_t failure_hash_min = 0U;
  std::uint64_t failure_hash_max = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&failure_hash, &failure_hash_min, 1, MPI_UINT64_T,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&failure_hash, &failure_hash_max, 1, MPI_UINT64_T,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(failure_hash_min == failure_hash_max);
  before_index = 0U;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    const auto logical = topology.global_cell(cell);
    for (int component = 0; component < 3; ++component)
      HUNDUN_CHECK(
          before_failure[before_index++] ==
          bits(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                        logical.z - box.begin.z, component)));
  }

  if (mpi.size() > 1) {
    std::string collective_message;
    try {
      const auto selected_output =
          mpi.rank() == 1 ? binding(storage, access, ids.scalar)
                          : binding(storage, access, ids.scalar_gradient);
      provider.compute_gradient(finite_volume::GradientScheme::green_gauss,
                                finite_volume::FiniteVolumeQuantity::density(),
                                binding(storage, access, ids.scalar),
                                selected_output);
    } catch (const runtime::Error &error) {
      collective_message = error.what();
    }
    HUNDUN_CHECK(collective_message.find("lowest failing rank 1") !=
                 std::string::npos);
    std::uint64_t message_hash = UINT64_C(14695981039346656037);
    for (const char character : collective_message) {
      message_hash ^= static_cast<unsigned char>(character);
      message_hash *= UINT64_C(1099511628211);
    }
    std::uint64_t minimum = 0U;
    std::uint64_t maximum = 0U;
    HUNDUN_CHECK(MPI_Allreduce(&message_hash, &minimum, 1, MPI_UINT64_T,
                               MPI_MIN, mpi.comm()) == MPI_SUCCESS);
    HUNDUN_CHECK(MPI_Allreduce(&message_hash, &maximum, 1, MPI_UINT64_T,
                               MPI_MAX, mpi.comm()) == MPI_SUCCESS);
    HUNDUN_CHECK(minimum == maximum);

    std::string scheme_message;
    try {
      const auto selected_scheme =
          mpi.rank() == 1
              ? static_cast<finite_volume::GradientScheme>(UINT8_C(255))
              : finite_volume::GradientScheme::green_gauss;
      provider.compute_gradient(selected_scheme,
                                finite_volume::FiniteVolumeQuantity::density(),
                                binding(storage, access, ids.scalar),
                                binding(storage, access, ids.scalar_gradient));
    } catch (const runtime::Error &error) {
      scheme_message = error.what();
    }
    HUNDUN_CHECK(scheme_message.find("lowest failing rank 1") !=
                 std::string::npos);

    std::optional<finite_volume::FaceMassFlux> moved_flux;
    if (mpi.rank() == 1)
      moved_flux.emplace(std::move(mass_flux));
    std::string moved_message;
    try {
      provider.reconstruct_momentum_faces(
          mass_flux, binding(storage, access, ids.velocity),
          binding(storage, access, ids.velocity_faces));
    } catch (const runtime::Error &error) {
      moved_message = error.what();
    }
    HUNDUN_CHECK(moved_message.find("lowest failing rank 1") !=
                 std::string::npos);
  }

  if (ghost > 1) {
    auto narrow_halo = runtime::HaloExchange::create(
        decomposition,
        runtime::ExchangePlan::create(decomposition,
                                      decomposition.local_extent(), ghost - 1));
    std::string narrow_message;
    try {
      static_cast<void>(finite_volume::ImmersedReconstruction::create(
          topology, geometry, boundaries, domain, ghost_plan, decomposition,
          mpi, narrow_halo));
    } catch (const runtime::Error &error) {
      narrow_message = error.what();
    }
    HUNDUN_CHECK(narrow_message.find("insufficient") != std::string::npos);
  }
  if (mpi.size() > 1) {
    MPI_Comm reordered = MPI_COMM_NULL;
    HUNDUN_CHECK(MPI_Comm_split(mpi.comm(), 0, mpi.size() - mpi.rank(),
                                &reordered) == MPI_SUCCESS);
    {
      auto reordered_context = runtime::MpiContext::duplicate(reordered);
      std::string communicator_message;
      try {
        static_cast<void>(finite_volume::ImmersedReconstruction::create(
            topology, geometry, boundaries, domain, ghost_plan, decomposition,
            reordered_context, halo));
      } catch (const runtime::Error &error) {
        communicator_message = error.what();
      }
      HUNDUN_CHECK(communicator_message.find("communicator") !=
                   std::string::npos);
    }
    HUNDUN_CHECK(MPI_Comm_free(&reordered) == MPI_SUCCESS);
  }
  if (boundary_mode == BoundaryMode::open) {
    for (int k = -ghost; k < local.z + ghost; ++k)
      for (int j = -ghost; j < local.y + ghost; ++j)
        for (int i = -ghost; i < local.x + ghost; ++i)
          scalar(i, j, k, 0) = 0.0;
    std::vector<finite_volume::detail::ImmersedWallNormalGradient>
        homogeneous_wall_conditions;
    homogeneous_wall_conditions.reserve(global_links.size());
    for (const auto &link : global_links)
      homogeneous_wall_conditions.push_back({link.id, 0.0});
    finite_volume::detail::
        compute_pressure_gradient_with_wall_normal_constraints(
            provider, finite_volume::GradientScheme::weighted_least_squares,
            binding(storage, access, ids.scalar),
            binding(storage, access, ids.scalar_gradient),
            homogeneous_wall_conditions);
    for (const auto global_id : planned_rows) {
      const auto cell = local_owned_index(topology, global_id);
      const auto logical = topology.global_cell(cell);
      for (int component = 0; component < 3; ++component)
        HUNDUN_CHECK(
            bits(gradient(logical.x - box.begin.x, logical.y - box.begin.y,
                          logical.z - box.begin.z, component)) == bits(0.0));
    }
  }
  HUNDUN_CHECK(provider.dependency_fingerprint() == fingerprint);
}

} // namespace

int main(int argc, char **argv) {
  runtime::MpiEnvironment environment(argc, argv);
  auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return test::run([&] {
    const std::string mode = argc > 1 ? argv[1] : "uniform";
    if (mode == "uniform")
      run_success(mpi, false, BoundaryMode::no_slip_wall);
    else if (mode == "warped")
      run_success(mpi, true, BoundaryMode::no_slip_wall);
    else if (mode == "periodic")
      run_success(mpi, false, BoundaryMode::periodic);
    else if (mode == "symmetry")
      run_success(mpi, false, BoundaryMode::symmetry);
    else if (mode == "open")
      run_success(mpi, false, BoundaryMode::open);
    else
      throw runtime::Error("unknown immersed reconstruction test mode");
  });
}
