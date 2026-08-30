// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_initialization.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string_view>
#include <vector>

namespace hundun::v04::detail {

enum class FreshProjectionFailurePoint : std::uint8_t {
  none,
  exact_operator_boundary,
  prepare_mg_boundary,
  audit_boundary
};

void set_fresh_projection_failure_for_test(FreshProjectionFailurePoint point,
                                           int failing_rank) noexcept;
void clear_fresh_projection_failure_for_test() noexcept;

} // namespace hundun::v04::detail

namespace {

using namespace hundun::v04;

bool collective(bool local, MPI_Comm communicator = MPI_COMM_WORLD) {
  const int value = local ? 1 : 0;
  int result = 0;
  return MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MIN, communicator) ==
             MPI_SUCCESS &&
         result != 0;
}

bool expect(bool local, int rank, std::string_view message) {
  const bool passed = collective(local);
  if (!passed && !local)
    std::cerr << "rank " << rank << " FAIL: " << message << '\n';
  return passed;
}

bool same_status(Status status) {
  const std::uint64_t packed =
      (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&packed, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(&packed, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

bool same_u64(std::uint64_t value) {
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField make_field(FieldId id, Int3 cells, std::uint8_t components,
                      std::uint8_t ghosts, StorageIdentity identity) {
  OwnedField result;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.storage.assign(nx * ny * nz * components, 0.0);
  result.view.base =
      result.storage.data() + ghosts + ghosts * nx + ghosts * nx * ny;
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = components;
  result.view.stride_y = nx;
  result.view.stride_z = nx * ny;
  result.view.component_stride = nx * ny * nz;
  result.view.field = id;
  result.view.revision = 1U;
  result.view.storage_identity = identity;
  result.view.revision_domain = identity + 100000U;
  return result;
}

struct OwnedFace {
  std::vector<double> storage;
  FaceFieldView view{};
};

OwnedFace make_face(CartesianAxis axis, Int3 cells, StorageIdentity identity) {
  OwnedFace result;
  Int3 extents = cells;
  if (axis == CartesianAxis::x)
    ++extents.x;
  if (axis == CartesianAxis::y)
    ++extents.y;
  if (axis == CartesianAxis::z)
    ++extents.z;
  result.storage.assign(
      static_cast<std::size_t>(extents.x) * extents.y * extents.z, 0.0);
  result.view = {result.storage.data(),
                 extents,
                 static_cast<std::size_t>(extents.x),
                 static_cast<std::size_t>(extents.x) * extents.y,
                 axis,
                 identity,
                 identity + 200000U};
  return result;
}

struct OwnedFlux {
  OwnedFace x;
  OwnedFace y;
  OwnedFace z;
  FaceFluxView view{};
};

OwnedFlux make_flux(Int3 cells, StorageIdentity identity,
                    RevisionToken revision) {
  OwnedFlux result;
  result.x = make_face(CartesianAxis::x, cells, identity);
  result.y = make_face(CartesianAxis::y, cells, identity);
  result.z = make_face(CartesianAxis::z, cells, identity);
  result.view = {result.x.view, result.y.view, result.z.view, revision, {}};
  return result;
}

CartesianMeshSpec mesh_spec(Int3 cells) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = cells;
  mesh.minimum_spacing = {1.0e-12, 1.0e-12, 1.0e-12};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {1000000U, 1U << 29U};
  return mesh;
}

std::size_t flat(Int3 extents, Int3 index) {
  return static_cast<std::size_t>(index.x) +
         static_cast<std::size_t>(extents.x) *
             (static_cast<std::size_t>(index.y) +
              static_cast<std::size_t>(extents.y) * index.z);
}

bool owns(const MeshPatch &patch, Int3 global) {
  return global.x >= patch.begin.x &&
         global.x < patch.begin.x + patch.cells.x &&
         global.y >= patch.begin.y &&
         global.y < patch.begin.y + patch.cells.y &&
         global.z >= patch.begin.z && global.z < patch.begin.z + patch.cells.z;
}

struct Fixture {
  int rank{};
  int size{};
  bool fragmented{};
  bool x_only_periodic{};
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  CartesianKernelPlan kernels;
  OwnedField density;
  OwnedField velocity;
  OwnedField velocity_previous;
  OwnedField velocity_trial;
  OwnedField chi;
  OwnedField rhs;
  OwnedField diagonal;
  OwnedField candidate_velocity;
  OwnedFace physical_x;
  OwnedFace physical_y;
  OwnedFace physical_z;
  OwnedFace solver_x;
  OwnedFace solver_y;
  OwnedFace solver_z;
  OwnedFlux flux;
  OwnedFlux candidate_flux;
  std::vector<std::uint8_t> active_cells;
  std::vector<std::uint8_t> active_x;
  std::vector<std::uint8_t> active_y;
  std::vector<std::uint8_t> active_z;
  LinearWorkspaceRequirements linear_requirements{};
  OwnedField linear_vectors;
  OwnedField linear_scalars;
  SolverWorkspace linear_workspace;
  ReductionEngine reductions;
  HaloEngine operator_halo;
  HaloEngine correction_halo;
  MgWorkspaceRequirements mg_requirements{};
  OwnedField mg_vectors;
  MgWorkspace mg_workspace;
  HaloEngine mg_halo;
  std::vector<HaloEngine> coarse_halos;
  std::vector<HaloEngine *> coarse_halo_pointers;
  FreshStartKinematicProjectionPlan plan;
  LinearSolveResult last_solve{};
  Status compile_status{StatusCode::invalid_plan, 0U};
  const char *failure_stage{"none"};
  double fixed_flux_sentinel{-3.75};

  bool global_active(Int3 global) const noexcept {
    if (fragmented)
      return (global.x & 1) == 0;
    const int nx = geometry.global_cells().x;
    return global.x != nx / 3 && global.x != 2 * nx / 3;
  }

  bool initialize(FreshStartProjectionLinearRoute route, Int3 global_cells,
                  bool use_fragmented = false, bool use_x_only_periodic = false,
                  bool mismatch_shared_mask = false) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    fragmented = use_fragmented;
    x_only_periodic = use_x_only_periodic;
    ValidatedModel model;
    model.mesh = mesh_spec(global_cells);
    model.fingerprint = 0x97260001U;
    model.pressure_reference = PressureReferenceKind::closed_mass;
    for (std::size_t index = 0U; index < model.boundaries.size(); ++index) {
      model.boundaries[index].flow_kind = !x_only_periodic || index < 2U
                                              ? BoundaryKind::periodic
                                              : BoundaryKind::no_slip_wall;
      model.boundaries[index].thermal_kind = BoundaryKind::none;
    }
    model.schemes.momentum = ConvectionScheme::central2;
    model.schemes.enthalpy = ConvectionScheme::central2;
    model.schemes.species = ConvectionScheme::central2;
    model.schemes.passive_scalar = ConvectionScheme::central2;
    FieldRegistry registry;
    FieldId id = 0U;
    const bool dependencies =
        registry.require_field("rho", 1U, 2U, id) &&
        registry.require_field("U", 3U, 2U, id) &&
        registry.require_field("pi", 1U, 2U, id) &&
        registry.require_field("h", 1U, 2U, id) &&
        registry.require_field("T", 1U, 2U, id) &&
        CartesianGeometryCompiler::compile(MPI_COMM_WORLD, model.mesh, {},
                                           geometry, patch) &&
        BoundaryCompiler::compile(MPI_COMM_WORLD, model, geometry, patch,
                                  registry, boundary, schemes, time) &&
        CartesianKernelPlan::compile(schemes, geometry, patch, boundary,
                                     kernels);
    if (!collective(dependencies)) {
      failure_stage = "dependencies";
      if (rank == 0)
        std::cerr << "fixture setup failed: " << failure_stage << '\n';
      return false;
    }

    const Int3 cells = patch.cells;
    const StorageIdentity salt =
        9726000U + static_cast<StorageIdentity>(route) * 1000U;
    density = make_field(10U, cells, 1U, 1U, salt + 1U);
    velocity = make_field(11U, cells, 3U, 0U, salt + 2U);
    velocity_previous = make_field(16U, cells, 3U, 0U, salt + 3U);
    velocity_trial = make_field(17U, cells, 3U, 0U, salt + 4U);
    chi = make_field(12U, cells, 1U, 1U, salt + 5U);
    rhs = make_field(13U, cells, 1U, 0U, salt + 6U);
    diagonal = make_field(14U, cells, 1U, 0U, salt + 7U);
    candidate_velocity = make_field(15U, cells, 3U, 0U, salt + 8U);
    physical_x = make_face(CartesianAxis::x, cells, salt + 9U);
    physical_y = make_face(CartesianAxis::y, cells, salt + 10U);
    physical_z = make_face(CartesianAxis::z, cells, salt + 11U);
    solver_x = make_face(CartesianAxis::x, cells, salt + 12U);
    solver_y = make_face(CartesianAxis::y, cells, salt + 13U);
    solver_z = make_face(CartesianAxis::z, cells, salt + 14U);
    flux = make_flux(cells, salt + 20U, 101U);
    candidate_flux = make_flux(cells, salt + 30U, 102U);
    std::fill(density.storage.begin(), density.storage.end(), 1.0);
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          velocity.view.unchecked({x, y, z}, 0U) = 1.0;
          velocity_previous.view.unchecked({x, y, z}, 0U) = 1.0;
          velocity_trial.view.unchecked({x, y, z}, 0U) = 1.0;
        }

    const Int3 xe{cells.x + 1, cells.y, cells.z};
    const Int3 ye{cells.x, cells.y + 1, cells.z};
    const Int3 ze{cells.x, cells.y, cells.z + 1};
    active_cells.assign(static_cast<std::size_t>(cells.x) * cells.y * cells.z,
                        0U);
    active_x.assign(static_cast<std::size_t>(xe.x) * xe.y * xe.z, 0U);
    active_y.assign(static_cast<std::size_t>(ye.x) * ye.y * ye.z, 0U);
    active_z.assign(static_cast<std::size_t>(ze.x) * ze.y * ze.z, 0U);
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 global{patch.begin.x + x, patch.begin.y + y,
                            patch.begin.z + z};
          active_cells[flat(cells, {x, y, z})] =
              global_active(global) ? 1U : 0U;
        }
    for (std::int32_t z = 0; z < xe.z; ++z)
      for (std::int32_t y = 0; y < xe.y; ++y)
        for (std::int32_t x = 0; x < xe.x; ++x) {
          const int face = patch.begin.x + x;
          const Int3 left{(face - 1 + global_cells.x) % global_cells.x,
                          patch.begin.y + y, patch.begin.z + z};
          const Int3 right{face % global_cells.x, patch.begin.y + y,
                           patch.begin.z + z};
          const bool active = global_active(left) && global_active(right);
          active_x[flat(xe, {x, y, z})] = active ? 1U : 0U;
          flux.x.view.unchecked({x, y, z}) = active && !fragmented ? 1.0 : 0.0;
        }
    for (std::int32_t z = 0; z < ye.z; ++z)
      for (std::int32_t y = 0; y < ye.y; ++y)
        for (std::int32_t x = 0; x < ye.x; ++x) {
          const int face = patch.begin.y + y;
          const Int3 low{patch.begin.x + x,
                         (face - 1 + global_cells.y) % global_cells.y,
                         patch.begin.z + z};
          const Int3 high{patch.begin.x + x, face % global_cells.y,
                          patch.begin.z + z};
          active_y[flat(ye, {x, y, z})] =
              global_active(low) && global_active(high) ? 1U : 0U;
        }
    for (std::int32_t z = 0; z < ze.z; ++z)
      for (std::int32_t y = 0; y < ze.y; ++y)
        for (std::int32_t x = 0; x < ze.x; ++x) {
          const int face = patch.begin.z + z;
          const Int3 low{patch.begin.x + x, patch.begin.y + y,
                         (face - 1 + global_cells.z) % global_cells.z};
          const Int3 high{patch.begin.x + x, patch.begin.y + y,
                          face % global_cells.z};
          active_z[flat(ze, {x, y, z})] =
              global_active(low) && global_active(high) ? 1U : 0U;
        }

    if (!fragmented) {
      const Int3 solid_global{global_cells.x / 3, 0, 0};
      if (owns(patch, solid_global)) {
        const Int3 local{solid_global.x - patch.begin.x,
                         solid_global.y - patch.begin.y,
                         solid_global.z - patch.begin.z};
        flux.y.view.unchecked(local) = fixed_flux_sentinel;
      }
    }
    if (mismatch_shared_mask && size > 1 && rank == 0) {
      if (patch.process_grid.x > 1)
        active_x[flat(xe, {cells.x, 0, 0})] ^= 1U;
      else if (patch.process_grid.y > 1)
        active_y[flat(ye, {0, cells.y, 0})] ^= 1U;
      else
        active_z[flat(ze, {0, 0, cells.z})] ^= 1U;
    }

    const LinearAlgorithm algorithm =
        route == FreshStartProjectionLinearRoute::native_mg_fgmres
            ? LinearAlgorithm::fgmres
            : LinearAlgorithm::pcg;
    const std::uint32_t restart =
        algorithm == LinearAlgorithm::fgmres ? 20U : 0U;
    bool services_ready = static_cast<bool>(make_linear_workspace_requirements(
        algorithm, cells, 1U, restart, ReductionMode::mpi_allreduce, 201U,
        linear_requirements));
    if (!collective(services_ready)) {
      failure_stage = "linear requirements";
      if (rank == 0)
        std::cerr << "fixture setup failed: " << failure_stage << '\n';
      return false;
    }
    linear_vectors = make_field(30U, cells, linear_requirements.vector_slots,
                                1U, salt + 40U);
    linear_scalars = make_field(
        31U,
        {static_cast<std::int32_t>(linear_requirements.scalar_doubles), 1, 1},
        1U, 0U, salt + 40U);
    services_ready =
        services_ready &&
        SolverWorkspace::bind(linear_requirements, linear_vectors.view,
                              linear_scalars.view, linear_workspace) &&
        ReductionEngine::compile(MPI_COMM_WORLD, ReductionMode::mpi_allreduce,
                                 linear_requirements.reduction_capacity,
                                 reductions);
    if (!collective(services_ready)) {
      failure_stage = "linear workspace/reductions";
      if (rank == 0)
        std::cerr << "fixture setup failed: " << failure_stage << '\n';
      return false;
    }
    const std::array<HaloFieldSpec, 1U> operator_fields{
        {{linear_vectors.view.field, 1U, 1U}}};
    const std::array<HaloFieldSpec, 1U> correction_fields{
        {{chi.view.field, 1U, 1U}}};
    services_ready =
        services_ready &&
        operator_halo.reserve(MPI_COMM_WORLD, patch,
                              {operator_fields.data(), operator_fields.size()},
                              boundary.halo_topology()) &&
        correction_halo.reserve(
            MPI_COMM_WORLD, patch,
            {correction_fields.data(), correction_fields.size()},
            boundary.halo_topology());
    if (!collective(services_ready)) {
      failure_stage = "fine halos";
      if (rank == 0)
        std::cerr << "fixture setup failed: " << failure_stage << '\n';
      return false;
    }
    MgRuntimeServices mg_services{};
    if (route == FreshStartProjectionLinearRoute::native_mg_fgmres) {
      MgHierarchyPolicy policy;
      policy.pre_sweeps = 1U;
      policy.post_sweeps = 2U;
      policy.point_smoother = MgPointSmootherKind::chebyshev_jacobi;
      policy.cycle = MgCycleKind::f_cycle;
      policy.maximum_levels = 2U;
      services_ready = static_cast<bool>(make_mg_workspace_requirements(
          MPI_COMM_WORLD, geometry, patch, policy, 202U, mg_requirements));
      if (!collective(services_ready) || mg_requirements.level_count == 0U) {
        failure_stage = "MG requirements";
        if (rank == 0)
          std::cerr << "fixture setup failed: " << failure_stage << '\n';
        return false;
      }
      mg_vectors =
          make_field(32U, mg_requirements.arena_shape, 1U, 1U, salt + 42U);
      services_ready =
          services_ready &&
          MgWorkspace::bind(mg_requirements, mg_vectors.view, mg_workspace);
      if (!collective(services_ready)) {
        failure_stage = "MG workspace";
        if (rank == 0)
          std::cerr << "fixture setup failed: " << failure_stage << '\n';
        return false;
      }
      const std::array<HaloFieldSpec, 1U> mg_fields{
          {{mg_vectors.view.field, 1U, 1U}}};
      services_ready =
          services_ready &&
          mg_halo.reserve(MPI_COMM_WORLD, mg_requirements.levels[0U].patch,
                          {mg_fields.data(), mg_fields.size()},
                          boundary.halo_topology());
      if (!collective(services_ready)) {
        failure_stage = "MG finest halo";
        if (rank == 0)
          std::cerr << "fixture setup failed: " << failure_stage << '\n';
        return false;
      }
      coarse_halos.resize(mg_requirements.level_count - 1U);
      coarse_halo_pointers.resize(coarse_halos.size());
      for (std::size_t level = 1U; level < mg_requirements.level_count;
           ++level) {
        services_ready =
            services_ready &&
            coarse_halos[level - 1U].reserve(
                MPI_COMM_WORLD, mg_requirements.levels[level].patch,
                {mg_fields.data(), mg_fields.size()}, boundary.halo_topology());
        coarse_halo_pointers[level - 1U] = &coarse_halos[level - 1U];
      }
      mg_services = {
          &mg_halo,
          &reductions,
          &mg_workspace,
          {coarse_halo_pointers.data(), coarse_halo_pointers.size()}};
    }
    if (!collective(services_ready))
      return false;

    FreshStartKinematicProjectionSpec spec;
    spec.communicator = MPI_COMM_WORLD;
    spec.geometry = &geometry;
    spec.kernels = &kernels;
    spec.boundary = &boundary;
    spec.patch = patch;
    spec.activity = {{active_cells.data(), active_cells.size()},
                     {active_x.data(), active_x.size()},
                     {active_y.data(), active_y.size()},
                     {active_z.data(), active_z.size()},
                     301U,
                     302U};
    spec.route = route;
    spec.solve = route == FreshStartProjectionLinearRoute::native_mg_fgmres
                     ? LinearSolveControl{1.0e-13, 1.0e-12, 500U, 4U, restart}
                     : LinearSolveControl{1.0e-11, 5.0e-10, 500U, 4U, restart};
    if (route == FreshStartProjectionLinearRoute::jacobi_pcg_oracle)
      spec.continuity_absolute_tolerance = 5.0e-8;
    spec.mg_policy.pre_sweeps = 1U;
    spec.mg_policy.post_sweeps = 2U;
    spec.mg_policy.point_smoother = MgPointSmootherKind::chebyshev_jacobi;
    spec.mg_policy.cycle = MgCycleKind::f_cycle;
    spec.mg_policy.maximum_levels = 2U;
    const FreshStartKinematicProjectionServices services{
        &operator_halo,   401U,       linear_vectors.view.field,
        &correction_halo, 402U,       &linear_workspace,
        &reductions,      mg_services};
    const FreshStartKinematicProjectionWorkspace workspace{
        chi.view,           rhs.view,
        diagonal.view,      physical_x.view,
        physical_y.view,    physical_z.view,
        solver_x.view,      solver_y.view,
        solver_z.view,      candidate_velocity.view,
        candidate_flux.view};
    compile_status = FreshStartKinematicProjectionPlan::compile(
        spec, services, workspace, plan);
    if (!compile_status)
      failure_stage = "projection compile";
    return static_cast<bool>(compile_status);
  }

  Status project(FreshStartKinematicProjectionCandidateCertificate &out) {
    FreshStartKinematicProjectionPreparedCertificate prepared;
    Status status =
        plan.prepare({as_const(density.view), as_const(velocity.view),
                      as_const(velocity_previous.view),
                      as_const(velocity_trial.view), as_const(flux.view), 501U},
                     prepared);
    if (!status)
      return status;
    FreshStartKinematicProjectionSolvedCertificate solved;
    status = plan.solve(prepared, solved);
    if (!status)
      return status;
    last_solve = solved.result();
    return plan.audit(solved, out);
  }
};

bool test_distributed_projection(int rank, int size) {
  Fixture fixture;
  const Int3 global{8 * size, 6, 4};
  bool passed =
      expect(fixture.initialize(
                 FreshStartProjectionLinearRoute::native_mg_fgmres, global),
             rank, "distributed production projection compiles");
  if (!passed)
    return false;
  const FreshStartKinematicProjectionRedCertificate &red = fixture.plan.red();
  passed &= expect(red.valid() && red.component_count == 2U &&
                       red.anchored_component_count == 2U &&
                       red.component_collective_payload_u64 > 0U &&
                       red.component_collective_global_payload_u64 >=
                           red.component_collective_payload_u64 &&
                       red.component_collective_global_payload_u64 <=
                           red.component_collective_payload_u64 *
                               static_cast<std::uint64_t>(size) &&
                       same_u64(red.component_graph) &&
                       same_u64(red.component_collective_payload_u64) &&
                       same_u64(red.component_collective_global_payload_u64),
                   rank,
                   "distributed graph has two compatible components and "
                   "collective payload bounds");
  FreshStartKinematicProjectionCandidateCertificate candidate;
  const Status projected = fixture.project(candidate);
  if (!projected && rank == 0)
    std::cerr << "distributed projection status "
              << static_cast<unsigned>(projected.code) << ':'
              << projected.detail
              << " solve residual=" << fixture.last_solve.final_true_residual
              << " iterations=" << fixture.last_solve.iterations << '\n';
  if (!projected && projected.detail == 9725U) {
    double local_maximum = 0.0;
    for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y)
        for (std::int32_t x = 0; x < fixture.patch.cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (fixture.active_cells[flat(fixture.patch.cells, cell)] == 0U)
            continue;
          const double residual =
              fixture.candidate_flux.x.view.unchecked({x + 1, y, z}) -
              fixture.candidate_flux.x.view.unchecked(cell) +
              fixture.candidate_flux.y.view.unchecked({x, y + 1, z}) -
              fixture.candidate_flux.y.view.unchecked(cell) +
              fixture.candidate_flux.z.view.unchecked({x, y, z + 1}) -
              fixture.candidate_flux.z.view.unchecked(cell);
          local_maximum = std::max(local_maximum, std::abs(residual));
        }
    double global_maximum = 0.0;
    MPI_Allreduce(&local_maximum, &global_maximum, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    if (rank == 0)
      std::cerr << "distributed rejected continuity=" << global_maximum << '\n';
  }
  if (projected && !(candidate.final_continuity_maximum() < 2.0e-10) &&
      rank == 0)
    std::cerr << "distributed final continuity "
              << candidate.final_continuity_maximum() << '\n';
  passed &= expect(
      same_status(projected) && projected && candidate.valid() &&
          candidate.final_continuity_maximum() < 2.0e-10,
      rank, "FGMRES/Native-MG closes distributed continuity collectively");
  if (!passed)
    return false;

  int local_anchor_count = 0;
  bool local_anchor_zero = true;
  for (const Int3 global_anchor :
       {Int3{0, 0, 0}, Int3{global.x / 3 + 1, 0, 0}}) {
    if (!owns(fixture.patch, global_anchor))
      continue;
    const Int3 local{global_anchor.x - fixture.patch.begin.x,
                     global_anchor.y - fixture.patch.begin.y,
                     global_anchor.z - fixture.patch.begin.z};
    ++local_anchor_count;
    const double value = fixture.chi.view.unchecked(local, 0U);
    const bool incident_solver_faces_zero =
        fixture.solver_x.view.unchecked(local) == 0.0 &&
        fixture.solver_x.view.unchecked({local.x + 1, local.y, local.z}) ==
            0.0 &&
        fixture.solver_y.view.unchecked(local) == 0.0 &&
        fixture.solver_y.view.unchecked({local.x, local.y + 1, local.z}) ==
            0.0 &&
        fixture.solver_z.view.unchecked(local) == 0.0 &&
        fixture.solver_z.view.unchecked({local.x, local.y, local.z + 1}) == 0.0;
    const bool physical_graph_survives =
        fixture.physical_x.view.unchecked(local) > 0.0 ||
        fixture.physical_x.view.unchecked({local.x + 1, local.y, local.z}) >
            0.0 ||
        fixture.physical_y.view.unchecked(local) > 0.0 ||
        fixture.physical_y.view.unchecked({local.x, local.y + 1, local.z}) >
            0.0 ||
        fixture.physical_z.view.unchecked(local) > 0.0 ||
        fixture.physical_z.view.unchecked({local.x, local.y, local.z + 1}) >
            0.0;
    local_anchor_zero = local_anchor_zero && std::abs(value) < 1.0e-10 &&
                        fixture.diagonal.view.unchecked(local, 0U) == 1.0 &&
                        fixture.rhs.view.unchecked(local, 0U) == 0.0 &&
                        incident_solver_faces_zero && physical_graph_survives;
  }
  int global_anchor_count = 0;
  MPI_Allreduce(&local_anchor_count, &global_anchor_count, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  passed &=
      expect(local_anchor_zero && global_anchor_count == 2, rank,
             "each disconnected component uses its global-minimum GID anchor");

  bool solid_candidate_zero = true;
  const Int3 cells = fixture.patch.cells;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (fixture.active_cells[flat(cells, cell)] != 0U)
          continue;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double value =
              fixture.candidate_velocity.view.unchecked(cell, component);
          const double positive_zero = 0.0;
          solid_candidate_zero =
              solid_candidate_zero &&
              std::memcmp(&value, &positive_zero, sizeof(value)) == 0;
        }
      }
  bool fixed_flux_unchanged = true;
  const Int3 sentinel_global{global.x / 3, 0, 0};
  if (owns(fixture.patch, sentinel_global)) {
    const Int3 local{sentinel_global.x - fixture.patch.begin.x,
                     sentinel_global.y - fixture.patch.begin.y,
                     sentinel_global.z - fixture.patch.begin.z};
    fixed_flux_unchanged =
        std::memcmp(&fixture.candidate_flux.y.view.unchecked(local),
                    &fixture.fixed_flux_sentinel,
                    sizeof(fixture.fixed_flux_sentinel)) == 0;
  }
  passed &= expect(solid_candidate_zero && fixed_flux_unchanged, rank,
                   "IBM candidate writes stationary-wall +0 and preserves "
                   "fixed cut flux bitwise");

  std::array<FieldView, 3U> layers{fixture.velocity.view,
                                   fixture.velocity_previous.view,
                                   fixture.velocity_trial.view};
  const Status committed = fixture.plan.commit(
      candidate, {layers.data(), layers.size()}, fixture.flux.view);
  bool solid_history_zero = true;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (fixture.active_cells[flat(cells, cell)] != 0U)
          continue;
        for (FieldView layer : layers)
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            const double value = layer.unchecked(cell, component);
            const double positive_zero = 0.0;
            solid_history_zero =
                solid_history_zero &&
                std::memcmp(&value, &positive_zero, sizeof(value)) == 0;
          }
      }
  passed &= expect(
      same_status(committed) && committed && solid_history_zero, rank,
      "joint commit publishes stationary-wall velocity to all three histories");
  return passed;
}

bool test_rank_local_incompatibility(int rank, int size) {
  Fixture fixture;
  bool passed = expect(
      fixture.initialize(FreshStartProjectionLinearRoute::jacobi_pcg_oracle,
                         {8 * size, 6, 4}),
      rank, "incompatible distributed fixture compiles");
  if (!passed)
    return false;
  if (rank == 0)
    fixture.flux.x.view.unchecked({0, 0, 0}) += 0.25;
  const auto candidate_velocity_before = fixture.candidate_velocity.storage;
  const auto candidate_x_before = fixture.candidate_flux.x.storage;
  const auto candidate_y_before = fixture.candidate_flux.y.storage;
  const auto candidate_z_before = fixture.candidate_flux.z.storage;
  FreshStartKinematicProjectionPreparedCertificate prepared;
  const Status status = fixture.plan.prepare(
      {as_const(fixture.density.view), as_const(fixture.velocity.view),
       as_const(fixture.velocity_previous.view),
       as_const(fixture.velocity_trial.view), as_const(fixture.flux.view),
       601U},
      prepared);
  passed &= expect(
      same_status(status) && status.code == StatusCode::rejected_step &&
          status.detail == 9723U && !prepared.valid() &&
          fixture.candidate_velocity.storage == candidate_velocity_before &&
          fixture.candidate_flux.x.storage == candidate_x_before &&
          fixture.candidate_flux.y.storage == candidate_y_before &&
          fixture.candidate_flux.z.storage == candidate_z_before,
      rank,
      "one-rank compatibility poison is collectively rejected with zero "
      "candidate writes");
  return passed;
}

bool test_collective_commit_transaction(int rank, int size) {
  Fixture fixture;
  bool passed = expect(
      fixture.initialize(FreshStartProjectionLinearRoute::native_mg_fgmres,
                         {8 * size, 6, 4}),
      rank, "collective commit transaction fixture compiles");
  if (!passed)
    return false;
  FreshStartKinematicProjectionCandidateCertificate candidate;
  Status status = fixture.project(candidate);
  passed &= expect(same_status(status) && status && candidate.valid(), rank,
                   "collective commit transaction produces a candidate");
  if (!passed)
    return false;

  const auto velocity_before = fixture.velocity.storage;
  const auto previous_before = fixture.velocity_previous.storage;
  const auto trial_before = fixture.velocity_trial.storage;
  const auto x_before = fixture.flux.x.storage;
  const auto y_before = fixture.flux.y.storage;
  const auto z_before = fixture.flux.z.storage;
  if (rank == 0)
    fixture.candidate_velocity.view.unchecked({0, 0, 0}, 0U) += 0.125;
  std::array<FieldView, 3U> layers{fixture.velocity.view,
                                   fixture.velocity_previous.view,
                                   fixture.velocity_trial.view};
  status = fixture.plan.commit(candidate, {layers.data(), layers.size()},
                               fixture.flux.view);
  passed &= expect(
      same_status(status) && status.code == StatusCode::invalid_plan &&
          status.detail == 9726U &&
          fixture.velocity.storage == velocity_before &&
          fixture.velocity_previous.storage == previous_before &&
          fixture.velocity_trial.storage == trial_before &&
          fixture.flux.x.storage == x_before &&
          fixture.flux.y.storage == y_before &&
          fixture.flux.z.storage == z_before,
      rank,
      "one-rank candidate tamper collectively rejects before every live store");
  if (!passed)
    return false;

  FreshStartKinematicProjectionCandidateCertificate recovered;
  status = fixture.project(recovered);
  passed &= expect(same_status(status) && status && recovered.valid(), rank,
                   "a fresh prepare/solve/audit recovers commit authority");
  if (!passed)
    return false;
  status = fixture.plan.commit(recovered, {layers.data(), layers.size()},
                               fixture.flux.view);
  passed &= expect(same_status(status) && status, rank,
                   "recovered collective authority commits jointly");
  return passed;
}

bool test_neighbor_mask_mismatch(int rank, int size) {
  if (size == 1)
    return true;
  Fixture fixture;
  const bool compiled =
      fixture.initialize(FreshStartProjectionLinearRoute::jacobi_pcg_oracle,
                         {8 * size, 6, 4}, false, false, true);
  return expect(!compiled && same_status(fixture.compile_status) &&
                    fixture.compile_status.code == StatusCode::invalid_plan,
                rank,
                "one-rank shared-face mask mismatch is collectively rejected");
}

bool test_payload_is_surface_bounded(int rank, int size) {
  Fixture small;
  Fixture large;
  bool passed = expect(
      small.initialize(FreshStartProjectionLinearRoute::jacobi_pcg_oracle,
                       {8 * size, 6, 4}, true, true),
      rank, "small fragmented MPI fixture compiles");
  passed &= expect(
      large.initialize(FreshStartProjectionLinearRoute::jacobi_pcg_oracle,
                       {16 * size, 6, 4}, true, true),
      rank, "large fragmented MPI fixture compiles");
  if (!passed)
    return false;
  const auto &small_red = small.plan.red();
  const auto &large_red = large.plan.red();
  passed &= expect(
      small.patch.process_grid.y == 1 && small.patch.process_grid.z == 1 &&
          large.patch.process_grid.y == 1 && large.patch.process_grid.z == 1 &&
          large_red.active_cells == 2U * small_red.active_cells &&
          large_red.component_count == 2U * small_red.component_count &&
          large_red.component_collective_payload_u64 ==
              small_red.component_collective_payload_u64 &&
          large_red.component_collective_global_payload_u64 ==
              small_red.component_collective_global_payload_u64,
      rank,
      "fixed partition surface has constant communication peak when interior "
      "fragments double");
  return passed;
}

bool test_rank_local_stage_failure_rendezvous(int rank, int size) {
  const int failing_rank = size - 1;
  const Int3 global{8 * size, 6, 4};
  bool passed = true;

  Fixture prepare_fixture;
  passed &=
      expect(prepare_fixture.initialize(
                 FreshStartProjectionLinearRoute::native_mg_fgmres, global),
             rank, "MG-boundary failure fixture compiles");
  if (!passed)
    return false;
  detail::set_fresh_projection_failure_for_test(
      detail::FreshProjectionFailurePoint::prepare_mg_boundary, failing_rank);
  FreshStartKinematicProjectionPreparedCertificate rejected_prepare;
  Status status = prepare_fixture.plan.prepare(
      {as_const(prepare_fixture.density.view),
       as_const(prepare_fixture.velocity.view),
       as_const(prepare_fixture.velocity_previous.view),
       as_const(prepare_fixture.velocity_trial.view),
       as_const(prepare_fixture.flux.view), 701U},
      rejected_prepare);
  detail::clear_fresh_projection_failure_for_test();
  passed &= expect(
      same_status(status) && status.code == StatusCode::invalid_plan &&
          status.detail == 9722U && !rejected_prepare.valid(),
      rank,
      "one-rank MG-boundary preparation failure rendezvous rejects every rank");
  if (!passed)
    return false;

  FreshStartKinematicProjectionPreparedCertificate prepared;
  status = prepare_fixture.plan.prepare(
      {as_const(prepare_fixture.density.view),
       as_const(prepare_fixture.velocity.view),
       as_const(prepare_fixture.velocity_previous.view),
       as_const(prepare_fixture.velocity_trial.view),
       as_const(prepare_fixture.flux.view), 702U},
      prepared);
  passed &= expect(same_status(status) && status && prepared.valid(), rank,
                   "MG-boundary failure leaves prepare retryable");
  if (!passed)
    return false;

  const LinearReductionCounters failure_reductions_before =
      prepare_fixture.reductions.counters();
  detail::set_fresh_projection_failure_for_test(
      detail::FreshProjectionFailurePoint::exact_operator_boundary,
      failing_rank);
  FreshStartKinematicProjectionSolvedCertificate rejected_solve;
  status = prepare_fixture.plan.solve(prepared, rejected_solve);
  detail::clear_fresh_projection_failure_for_test();
  const LinearReductionCounters failure_reductions_after =
      prepare_fixture.reductions.counters();
  const std::uint64_t failure_calls =
      failure_reductions_after.calls - failure_reductions_before.calls;
  const std::uint64_t failure_blocking_operations =
      failure_reductions_after.blocking_operations -
      failure_reductions_before.blocking_operations;
  const int observed_failing_rank =
      prepare_fixture.reductions.lowest_failing_rank();
  passed &=
      expect(same_status(status) && status.code == StatusCode::invalid_plan &&
                 status.detail == 9724U && !rejected_solve.valid(),
             rank,
             "one-rank exact-operator boundary failure rendezvous aborts the "
             "shared Krylov path");
  const bool collective_provenance = observed_failing_rank == failing_rank &&
                                     failure_calls == 2U &&
                                     failure_blocking_operations == 10U;
  if (!collective_provenance)
    std::cerr << "rank " << rank
              << " exact-operator collective provenance lowest="
              << observed_failing_rank << " expected=" << failing_rank
              << " calls=" << failure_calls
              << " blocking=" << failure_blocking_operations
              << " expected-blocking=10\n";
  passed &=
      expect(collective_provenance, rank,
             "exact-operator consensus preserves the lowest failing rank and "
             "Krylov consumes collective scope without a redundant reduction");
  if (!passed)
    return false;

  FreshStartKinematicProjectionSolvedCertificate recovered_solve;
  status = prepare_fixture.plan.solve(prepared, recovered_solve);
  const bool clean_solve_replay =
      same_status(status) && status && recovered_solve.valid() &&
      (recovered_solve.result().termination == LinearTermination::converged ||
       recovered_solve.result().termination == LinearTermination::zero_rhs) &&
      recovered_solve.result().lowest_failing_rank == -1 &&
      prepare_fixture.reductions.lowest_failing_rank() == -1;
  passed &= expect(
      clean_solve_replay, rank,
      "exact-operator failure leaves solve replayable with cleared provenance");
  if (!passed)
    return false;

  Fixture audit_fixture;
  passed &=
      expect(audit_fixture.initialize(
                 FreshStartProjectionLinearRoute::native_mg_fgmres, global),
             rank, "audit-boundary failure fixture compiles");
  if (!passed)
    return false;
  FreshStartKinematicProjectionPreparedCertificate audit_prepared;
  status = audit_fixture.plan.prepare(
      {as_const(audit_fixture.density.view),
       as_const(audit_fixture.velocity.view),
       as_const(audit_fixture.velocity_previous.view),
       as_const(audit_fixture.velocity_trial.view),
       as_const(audit_fixture.flux.view), 703U},
      audit_prepared);
  FreshStartKinematicProjectionSolvedCertificate audit_solved;
  if (status)
    status = audit_fixture.plan.solve(audit_prepared, audit_solved);
  passed &= expect(same_status(status) && status && audit_solved.valid(), rank,
                   "audit-boundary fixture reaches a solved projection");
  if (!passed)
    return false;
  detail::set_fresh_projection_failure_for_test(
      detail::FreshProjectionFailurePoint::audit_boundary, failing_rank);
  FreshStartKinematicProjectionCandidateCertificate rejected_audit;
  status = audit_fixture.plan.audit(audit_solved, rejected_audit);
  detail::clear_fresh_projection_failure_for_test();
  passed &=
      expect(same_status(status) && status.code == StatusCode::invalid_plan &&
                 status.detail == 9725U && !rejected_audit.valid(),
             rank,
             "one-rank post-correction-halo boundary failure rendezvous clears "
             "candidate authority");
  if (!passed)
    return false;
  FreshStartKinematicProjectionCandidateCertificate recovered;
  status = audit_fixture.plan.audit(audit_solved, recovered);
  passed &= expect(same_status(status) && status && recovered.valid(), rank,
                   "post-halo boundary failure leaves audit replayable");
  return passed;
}

} // namespace

int main(int argc, char **argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    return 2;
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  bool passed = test_distributed_projection(rank, size);
  passed &= test_collective_commit_transaction(rank, size);
  passed &= test_rank_local_incompatibility(rank, size);
  passed &= test_neighbor_mask_mismatch(rank, size);
  passed &= test_payload_is_surface_bounded(rank, size);
  passed &= test_rank_local_stage_failure_rendezvous(rank, size);
  const bool global_passed = collective(passed);
  if (rank == 0 && global_passed)
    std::cout << "fresh-start kinematic projection MPI tests passed on " << size
              << " rank(s)\n";
  MPI_Finalize();
  return global_passed ? 0 : 1;
}
