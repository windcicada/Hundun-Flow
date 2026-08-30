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

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view message) {
  if (!condition)
    std::cerr << "FAIL: " << message << '\n';
  return condition;
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

struct Fixture {
  FreshStartProjectionLinearRoute route{};
  bool fragmented_activity{};
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

  bool global_active(Int3 global) const noexcept {
    if (fragmented_activity)
      return (global.x & 1) == 0;
    const int nx = geometry.global_cells().x;
    const int first_solid = nx / 3;
    const int second_solid = 2 * nx / 3;
    return global.x != first_solid && global.x != second_solid;
  }

  bool initialize(FreshStartProjectionLinearRoute selected,
                  bool incompatible = false, bool no_ibm = false,
                  bool asymmetric_activity = false,
                  bool alias_workspace = false, Int3 global_cells = {8, 6, 4},
                  bool x_only_periodic = false, bool fragmented = false,
                  bool alias_face_with_linear_workspace = false) {
    route = selected;
    fragmented_activity = fragmented;
    ValidatedModel model;
    model.mesh = mesh_spec(global_cells);
    model.fingerprint = 0x97210001U;
    model.pressure_reference = PressureReferenceKind::closed_mass;
    for (std::size_t face_index = 0U; face_index < model.boundaries.size();
         ++face_index) {
      BoundaryFaceSpec &face = model.boundaries[face_index];
      face.flow_kind = !x_only_periodic || face_index < 2U
                           ? BoundaryKind::periodic
                           : BoundaryKind::no_slip_wall;
      face.thermal_kind = BoundaryKind::none;
    }
    model.schemes.momentum = ConvectionScheme::central2;
    model.schemes.enthalpy = ConvectionScheme::central2;
    model.schemes.species = ConvectionScheme::central2;
    model.schemes.passive_scalar = ConvectionScheme::central2;
    FieldRegistry registry;
    FieldId id = 0U;
    if (!registry.require_field("rho", 1U, 2U, id) ||
        !registry.require_field("U", 3U, 2U, id) ||
        !registry.require_field("pi", 1U, 2U, id) ||
        !registry.require_field("h", 1U, 2U, id) ||
        !registry.require_field("T", 1U, 2U, id) ||
        !CartesianGeometryCompiler::compile(MPI_COMM_SELF, model.mesh, {},
                                            geometry, patch) ||
        !BoundaryCompiler::compile(MPI_COMM_SELF, model, geometry, patch,
                                   registry, boundary, schemes, time) ||
        !CartesianKernelPlan::compile(schemes, geometry, patch, boundary,
                                      kernels)) {
      std::cerr << "fixture dependency compile failed\n";
      return false;
    }

    const Int3 cells = patch.cells;
    const StorageIdentity salt =
        selected == FreshStartProjectionLinearRoute::native_mg_fgmres
            ? 9721000U
            : 9722000U;
    density = make_field(10U, cells, 1U, 1U, salt + 1U);
    velocity = make_field(11U, cells, 3U, 0U, salt + 2U);
    velocity_previous = make_field(16U, cells, 3U, 0U, salt + 13U);
    velocity_trial = make_field(17U, cells, 3U, 0U, salt + 14U);
    chi = make_field(12U, cells, 1U, 1U, salt + 3U);
    rhs = make_field(13U, cells, 1U, 0U, salt + 4U);
    diagonal = make_field(14U, cells, 1U, 0U, salt + 5U);
    candidate_velocity = make_field(15U, cells, 3U, 0U, salt + 6U);
    physical_x = make_face(CartesianAxis::x, cells, salt + 7U);
    physical_y = make_face(CartesianAxis::y, cells, salt + 8U);
    physical_z = make_face(CartesianAxis::z, cells, salt + 9U);
    solver_x = make_face(CartesianAxis::x, cells, salt + 10U);
    solver_y = make_face(CartesianAxis::y, cells, salt + 11U);
    solver_z = make_face(CartesianAxis::z, cells, salt + 12U);
    flux = make_flux(cells, salt + 20U, 101U);
    candidate_flux = make_flux(cells, salt + 30U, 102U);
    std::fill(density.storage.begin(), density.storage.end(), 1.0);
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          velocity.view.unchecked({x, y, z}, 0U) = 1.0;
          velocity.view.unchecked({x, y, z}, 1U) = 0.0;
          velocity.view.unchecked({x, y, z}, 2U) = 0.0;
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            velocity_previous.view.unchecked({x, y, z}, component) =
                velocity.view.unchecked({x, y, z}, component);
            velocity_trial.view.unchecked({x, y, z}, component) =
                velocity.view.unchecked({x, y, z}, component);
          }
        }

    active_cells.assign(static_cast<std::size_t>(cells.x) * cells.y * cells.z,
                        0U);
    const Int3 xe{cells.x + 1, cells.y, cells.z};
    const Int3 ye{cells.x, cells.y + 1, cells.z};
    const Int3 ze{cells.x, cells.y, cells.z + 1};
    active_x.assign(static_cast<std::size_t>(xe.x) * xe.y * xe.z, 0U);
    active_y.assign(static_cast<std::size_t>(ye.x) * ye.y * ye.z, 0U);
    active_z.assign(static_cast<std::size_t>(ze.x) * ze.y * ze.z, 0U);
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 local{x, y, z};
          const Int3 global{patch.begin.x + x, patch.begin.y + y,
                            patch.begin.z + z};
          active_cells[flat(cells, local)] = global_active(global) ? 1U : 0U;
        }
    for (std::int32_t z = 0; z < xe.z; ++z)
      for (std::int32_t y = 0; y < xe.y; ++y)
        for (std::int32_t x = 0; x < xe.x; ++x) {
          const int global_face = patch.begin.x + x;
          const Int3 left{(global_face - 1 + global_cells.x) % global_cells.x,
                          patch.begin.y + y, patch.begin.z + z};
          const Int3 right{global_face % global_cells.x, patch.begin.y + y,
                           patch.begin.z + z};
          active_x[flat(xe, {x, y, z})] =
              global_active(left) && global_active(right) ? 1U : 0U;
          flux.x.view.unchecked({x, y, z}) =
              active_x[flat(xe, {x, y, z})] != 0U ? 1.0 : 0.0;
        }
    for (std::int32_t z = 0; z < ye.z; ++z)
      for (std::int32_t y = 0; y < ye.y; ++y)
        for (std::int32_t x = 0; x < ye.x; ++x) {
          const int gy = patch.begin.y + y;
          const Int3 low{patch.begin.x + x,
                         (gy - 1 + global_cells.y) % global_cells.y,
                         patch.begin.z + z};
          const Int3 high{patch.begin.x + x, gy % global_cells.y,
                          patch.begin.z + z};
          active_y[flat(ye, {x, y, z})] =
              global_active(low) && global_active(high) ? 1U : 0U;
        }
    for (std::int32_t z = 0; z < ze.z; ++z)
      for (std::int32_t y = 0; y < ze.y; ++y)
        for (std::int32_t x = 0; x < ze.x; ++x) {
          const int gz = patch.begin.z + z;
          const Int3 low{patch.begin.x + x, patch.begin.y + y,
                         (gz - 1 + global_cells.z) % global_cells.z};
          const Int3 high{patch.begin.x + x, patch.begin.y + y,
                          gz % global_cells.z};
          active_z[flat(ze, {x, y, z})] =
              global_active(low) && global_active(high) ? 1U : 0U;
        }
    if (incompatible)
      flux.x.view.unchecked({0, 0, 0}) += 0.25;
    if (asymmetric_activity) {
      const std::int32_t face_x = global_cells.x / 3 + 1;
      active_x[flat(xe, {face_x, 0, 0})] = 1U;
    }

    const LinearAlgorithm algorithm =
        selected == FreshStartProjectionLinearRoute::native_mg_fgmres
            ? LinearAlgorithm::fgmres
            : LinearAlgorithm::pcg;
    const std::uint32_t restart =
        algorithm == LinearAlgorithm::fgmres ? 20U : 0U;
    if (!make_linear_workspace_requirements(algorithm, cells, 1U, restart,
                                            ReductionMode::mpi_allreduce, 201U,
                                            linear_requirements)) {
      std::cerr << "linear requirements failed\n";
      return false;
    }
    linear_vectors = make_field(30U, cells, linear_requirements.vector_slots,
                                1U, salt + 40U);
    linear_scalars = make_field(
        31U,
        {static_cast<std::int32_t>(linear_requirements.scalar_doubles), 1, 1},
        1U, 0U, salt + 40U);
    if (!SolverWorkspace::bind(linear_requirements, linear_vectors.view,
                               linear_scalars.view, linear_workspace) ||
        !ReductionEngine::compile(MPI_COMM_SELF, ReductionMode::mpi_allreduce,
                                  linear_requirements.reduction_capacity,
                                  reductions)) {
      std::cerr << "linear workspace/reduction failed\n";
      return false;
    }
    const std::array<HaloFieldSpec, 1U> operator_fields{
        {{linear_vectors.view.field, 1U, 1U}}};
    const std::array<HaloFieldSpec, 1U> correction_fields{
        {{chi.view.field, 1U, 1U}}};
    if (!operator_halo.reserve(MPI_COMM_SELF, patch,
                               {operator_fields.data(), operator_fields.size()},
                               boundary.halo_topology()) ||
        !correction_halo.reserve(
            MPI_COMM_SELF, patch,
            {correction_fields.data(), correction_fields.size()},
            boundary.halo_topology())) {
      std::cerr << "operator/correction halo failed\n";
      return false;
    }

    MgRuntimeServices mg_services{};
    if (selected == FreshStartProjectionLinearRoute::native_mg_fgmres) {
      MgHierarchyPolicy policy;
      policy.pre_sweeps = 1U;
      policy.post_sweeps = 2U;
      policy.point_smoother = MgPointSmootherKind::chebyshev_jacobi;
      policy.cycle = MgCycleKind::f_cycle;
      policy.maximum_levels = 2U;
      const Status mg_required = make_mg_workspace_requirements(
          MPI_COMM_SELF, geometry, patch, policy, 202U, mg_requirements);
      if (!mg_required) {
        std::cerr << "mg requirements failed "
                  << static_cast<unsigned>(mg_required.code) << ':'
                  << mg_required.detail
                  << " levels=" << static_cast<unsigned>(policy.maximum_levels)
                  << " min="
                  << static_cast<unsigned>(policy.minimum_coarse_extent)
                  << " global=" << geometry.global_cells().x << ','
                  << geometry.global_cells().y << ','
                  << geometry.global_cells().z << " patch=" << patch.cells.x
                  << ',' << patch.cells.y << ',' << patch.cells.z << '\n';
        return false;
      }
      mg_vectors =
          make_field(32U, mg_requirements.arena_shape, 1U, 1U, salt + 42U);
      if (!MgWorkspace::bind(mg_requirements, mg_vectors.view, mg_workspace)) {
        std::cerr << "mg workspace failed\n";
        return false;
      }
      const std::array<HaloFieldSpec, 1U> mg_fields{
          {{mg_vectors.view.field, 1U, 1U}}};
      if (!mg_halo.reserve(MPI_COMM_SELF, mg_requirements.levels[0U].patch,
                           {mg_fields.data(), mg_fields.size()},
                           boundary.halo_topology())) {
        std::cerr << "mg finest halo failed\n";
        return false;
      }
      coarse_halos.resize(mg_requirements.level_count - 1U);
      coarse_halo_pointers.resize(coarse_halos.size());
      for (std::size_t level = 1U; level < mg_requirements.level_count;
           ++level) {
        if (!coarse_halos[level - 1U].reserve(
                MPI_COMM_SELF, mg_requirements.levels[level].patch,
                {mg_fields.data(), mg_fields.size()}, boundary.halo_topology()))
          return false;
        coarse_halo_pointers[level - 1U] = &coarse_halos[level - 1U];
      }
      mg_services = {
          &mg_halo,
          &reductions,
          &mg_workspace,
          {coarse_halo_pointers.data(), coarse_halo_pointers.size()}};
    }

    FreshStartKinematicProjectionSpec spec;
    spec.communicator = MPI_COMM_SELF;
    spec.geometry = &geometry;
    spec.kernels = &kernels;
    spec.boundary = &boundary;
    spec.patch = patch;
    if (!no_ibm)
      spec.activity = {{active_cells.data(), active_cells.size()},
                       {active_x.data(), active_x.size()},
                       {active_y.data(), active_y.size()},
                       {active_z.data(), active_z.size()},
                       301U,
                       302U};
    spec.route = selected;
    spec.solve = selected == FreshStartProjectionLinearRoute::jacobi_pcg_oracle
                     ? LinearSolveControl{1.0e-11, 5.0e-10, 500U, 4U, restart}
                     : LinearSolveControl{1.0e-11, 1.0e-10, 500U, 4U, restart};
    if (selected == FreshStartProjectionLinearRoute::jacobi_pcg_oracle)
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
    FreshStartKinematicProjectionWorkspace workspace{
        chi.view,           rhs.view,
        diagonal.view,      physical_x.view,
        physical_y.view,    physical_z.view,
        solver_x.view,      solver_y.view,
        solver_z.view,      candidate_velocity.view,
        candidate_flux.view};
    if (alias_workspace)
      workspace.diagonal = rhs.view;
    if (alias_face_with_linear_workspace)
      workspace.x_physical_mobility.base = linear_vectors.view.base;
    const Status compiled = FreshStartKinematicProjectionPlan::compile(
        spec, services, workspace, plan);
    if (!compiled && !alias_workspace && !asymmetric_activity &&
        !alias_face_with_linear_workspace)
      std::cerr << "projection compile status "
                << static_cast<unsigned>(compiled.code) << ':'
                << compiled.detail
                << " route=" << static_cast<unsigned>(selected) << '\n';
    return static_cast<bool>(compiled);
  }

  Status project(FreshStartKinematicProjectionCandidateCertificate &candidate) {
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
    return plan.audit(solved, candidate);
  }
};

std::vector<double> owned_values(ConstFieldView field) {
  std::vector<double> values;
  for (std::int32_t z = 0; z < field.interior.z; ++z)
    for (std::int32_t y = 0; y < field.interior.y; ++y)
      for (std::int32_t x = 0; x < field.interior.x; ++x)
        for (std::uint8_t component = 0U; component < field.components;
             ++component)
          values.push_back(field.unchecked({x, y, z}, component));
  return values;
}

bool test_red_projection_and_oracle() {
  Fixture production;
  Fixture oracle;
  bool passed = expect(
      production.initialize(FreshStartProjectionLinearRoute::native_mg_fgmres),
      "production FGMRES/Native-MG projection compiles");
  passed &= expect(
      oracle.initialize(FreshStartProjectionLinearRoute::jacobi_pcg_oracle),
      "fixed-SPD Jacobi/PCG oracle compiles");
  if (!passed)
    return false;
  const auto &red = production.plan.red();
  passed &= expect(
      red.valid() && red.component_count == 2U &&
          red.anchored_component_count == 2U &&
          red.physical_mobility_separate && red.per_component_compatibility &&
          red.global_minimum_gid_anchors && red.mg_null_space_none &&
          red.component_collective_payload_u64 > 0U &&
          red.component_collective_global_payload_u64 >=
              red.component_collective_payload_u64 &&
          red.runtime_workspace_preallocated &&
          red.no_immersed_bitwise_bypass && red.three_layer_joint_commit &&
          !red.chi_writes_pressure,
      "RED binds two disconnected compatible components, their anchors, and no "
      "pressure write");
  FreshStartKinematicProjectionCandidateCertificate production_candidate;
  FreshStartKinematicProjectionCandidateCertificate oracle_candidate;
  const Status production_status = production.project(production_candidate);
  const Status oracle_status = oracle.project(oracle_candidate);
  if (!production_status || !oracle_status)
    std::cerr << "projection statuses production="
              << static_cast<unsigned>(production_status.code) << ':'
              << production_status.detail
              << " oracle=" << static_cast<unsigned>(oracle_status.code) << ':'
              << oracle_status.detail << '\n';
  passed &=
      expect(production_status && oracle_status &&
                 production_candidate.valid() && oracle_candidate.valid() &&
                 production_candidate.final_continuity_maximum() < 2.0e-10 &&
                 oracle_candidate.final_continuity_maximum() < 5.0e-8,
             "production and independent PCG oracle both close continuity");
  if (!production_status || !oracle_status)
    return false;

  bool cut_faces_zero = true;
  for (std::size_t index = 0U; index < production.active_x.size(); ++index)
    if (production.active_x[index] == 0U)
      cut_faces_zero =
          cut_faces_zero && production.physical_x.storage[index] == 0.0;
  bool solid_velocity_is_wall_zero = true;
  const Int3 cells = production.patch.cells;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (production.active_cells[flat(cells, cell)] != 0U)
          continue;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double after =
              production.candidate_velocity.view.unchecked(cell, component);
          const double positive_zero = 0.0;
          solid_velocity_is_wall_zero =
              solid_velocity_is_wall_zero &&
              std::memcmp(&positive_zero, &after, sizeof(after)) == 0;
        }
      }
  passed &= expect(cut_faces_zero && solid_velocity_is_wall_zero &&
                       red.ibm_gradient_from_face_activity &&
                       red.cut_face_zero_mobility &&
                       red.inactive_velocity_stationary_wall_zero,
                   "G_IBM consumes the face activity: cut mobility is zero and "
                   "solid velocity is exact stationary-wall +0");

  const std::vector<double> production_chi =
      owned_values(as_const(production.chi.view));
  const std::vector<double> oracle_chi =
      owned_values(as_const(oracle.chi.view));
  double maximum_difference = 0.0;
  for (std::size_t index = 0U; index < production_chi.size(); ++index)
    maximum_difference =
        std::max(maximum_difference,
                 std::abs(production_chi[index] - oracle_chi[index]));
  passed &= expect(
      maximum_difference < 5.0e-8,
      "FGMRES/Native-MG and Jacobi/PCG solve the same anchored Schur system");

  bool separated_anchor_face = false;
  for (std::size_t index = 0U; index < production.physical_x.storage.size();
       ++index)
    separated_anchor_face =
        separated_anchor_face || (production.physical_x.storage[index] > 0.0 &&
                                  production.solver_x.storage[index] == 0.0);
  passed &=
      expect(separated_anchor_face, "physical mobility survives where the "
                                    "anchored solver off-diagonal is removed");

  const std::vector<double> velocity_before = production.velocity.storage;
  const std::vector<double> velocity_previous_before =
      production.velocity_previous.storage;
  const std::vector<double> velocity_trial_before =
      production.velocity_trial.storage;
  const std::vector<double> flux_before = production.flux.x.storage;
  std::array<FieldView, 3U> velocity_layers{production.velocity.view,
                                            production.velocity_previous.view,
                                            production.velocity_trial.view};
  passed &= expect(static_cast<bool>(production.plan.commit(
                       production_candidate,
                       {velocity_layers.data(), velocity_layers.size()},
                       production.flux.view)),
                   "joint no-fail U/phi candidate commits");
  passed &= expect(
      production.velocity.storage != velocity_before &&
          production.velocity_previous.storage != velocity_previous_before &&
          production.velocity_trial.storage != velocity_trial_before &&
          production.velocity.storage == production.velocity_previous.storage &&
          production.velocity.storage == production.velocity_trial.storage &&
          production.flux.x.storage != flux_before,
      "compatible cold start atomically publishes U to all three layers and "
      "mass flux");
  bool committed_solid_history_zero = true;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (production.active_cells[flat(cells, cell)] != 0U)
          continue;
        for (FieldView layer : velocity_layers)
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            const double value = layer.unchecked(cell, component);
            const double positive_zero = 0.0;
            committed_solid_history_zero =
                committed_solid_history_zero &&
                std::memcmp(&value, &positive_zero, sizeof(value)) == 0;
          }
      }
  passed &= expect(committed_solid_history_zero,
                   "all three committed solid velocity histories are exact "
                   "stationary-wall +0");
  return passed;
}

bool test_component_compatibility_rejects_without_candidate_write() {
  Fixture fixture;
  bool passed =
      expect(fixture.initialize(
                 FreshStartProjectionLinearRoute::native_mg_fgmres, true),
             "incompatible component fixture compiles");
  if (!passed)
    return false;
  const std::vector<double> candidate_velocity_before =
      fixture.candidate_velocity.storage;
  const std::vector<double> candidate_flux_before =
      fixture.candidate_flux.x.storage;
  FreshStartKinematicProjectionPreparedCertificate prepared;
  const Status rejected = fixture.plan.prepare(
      {as_const(fixture.density.view), as_const(fixture.velocity.view),
       as_const(fixture.velocity_previous.view),
       as_const(fixture.velocity_trial.view), as_const(fixture.flux.view),
       601U},
      prepared);
  if (rejected.code != StatusCode::rejected_step || rejected.detail != 9723U)
    std::cerr << "compatibility status " << static_cast<unsigned>(rejected.code)
              << ':' << rejected.detail << '\n';
  passed &= expect(
      rejected.code == StatusCode::rejected_step && rejected.detail == 9723U &&
          !prepared.valid() &&
          fixture.candidate_velocity.storage == candidate_velocity_before &&
          fixture.candidate_flux.x.storage == candidate_flux_before,
      "per-component incompatibility rejects before candidate U/phi writes");
  return passed;
}

bool test_tampered_candidate_cannot_partially_commit() {
  Fixture fixture;
  bool passed = expect(
      fixture.initialize(FreshStartProjectionLinearRoute::native_mg_fgmres),
      "tamper fixture compiles");
  FreshStartKinematicProjectionCandidateCertificate candidate;
  const Status projected = fixture.project(candidate);
  if (!projected)
    std::cerr << "tamper projection status "
              << static_cast<unsigned>(projected.code) << ':'
              << projected.detail << '\n';
  passed &= expect(static_cast<bool>(projected), "tamper fixture projects");
  if (!passed)
    return false;
  fixture.candidate_velocity.view.unchecked({0, 0, 0}, 0U) += 0.5;
  const std::vector<double> live_velocity = fixture.velocity.storage;
  const std::vector<double> live_previous = fixture.velocity_previous.storage;
  const std::vector<double> live_trial = fixture.velocity_trial.storage;
  const std::vector<double> live_x = fixture.flux.x.storage;
  const std::vector<double> live_y = fixture.flux.y.storage;
  const std::vector<double> live_z = fixture.flux.z.storage;
  std::array<FieldView, 3U> velocity_layers{fixture.velocity.view,
                                            fixture.velocity_previous.view,
                                            fixture.velocity_trial.view};
  const Status rejected = fixture.plan.commit(
      candidate, {velocity_layers.data(), velocity_layers.size()},
      fixture.flux.view);
  passed &= expect(
      rejected.code == StatusCode::invalid_plan &&
          fixture.velocity.storage == live_velocity &&
          fixture.velocity_previous.storage == live_previous &&
          fixture.velocity_trial.storage == live_trial &&
          fixture.flux.x.storage == live_x &&
          fixture.flux.y.storage == live_y && fixture.flux.z.storage == live_z,
      "tampered joint authority rejects before the first live U/phi store");
  return passed;
}

bool test_stale_base_and_single_use_authority() {
  Fixture audit_fixture;
  bool passed = expect(audit_fixture.initialize(
                           FreshStartProjectionLinearRoute::native_mg_fgmres),
                       "stale-audit fixture compiles");
  if (!passed)
    return false;
  FreshStartKinematicProjectionPreparedCertificate prepared;
  Status status = audit_fixture.plan.prepare(
      {as_const(audit_fixture.density.view),
       as_const(audit_fixture.velocity.view),
       as_const(audit_fixture.velocity_previous.view),
       as_const(audit_fixture.velocity_trial.view),
       as_const(audit_fixture.flux.view), 701U},
      prepared);
  FreshStartKinematicProjectionSolvedCertificate solved;
  if (status)
    status = audit_fixture.plan.solve(prepared, solved);
  passed &= expect(static_cast<bool>(status),
                   "stale-audit fixture prepares and solves");
  if (!passed)
    return false;
  audit_fixture.velocity.view.unchecked({0, 0, 0}, 0U) += 0.125;
  const auto candidate_velocity_before =
      audit_fixture.candidate_velocity.storage;
  const auto candidate_x_before = audit_fixture.candidate_flux.x.storage;
  FreshStartKinematicProjectionCandidateCertificate stale_candidate;
  const Status stale_audit = audit_fixture.plan.audit(solved, stale_candidate);
  passed &= expect(
      stale_audit.code == StatusCode::invalid_plan &&
          !stale_candidate.valid() &&
          audit_fixture.candidate_velocity.storage ==
              candidate_velocity_before &&
          audit_fixture.candidate_flux.x.storage == candidate_x_before,
      "base U mutation after solve rejects before audit candidate writes");

  Fixture post_audit_fixture;
  passed &= expect(post_audit_fixture.initialize(
                       FreshStartProjectionLinearRoute::native_mg_fgmres),
                   "post-audit stale-base fixture compiles");
  FreshStartKinematicProjectionCandidateCertificate post_audit_candidate;
  if (passed)
    status = post_audit_fixture.project(post_audit_candidate);
  passed &= expect(static_cast<bool>(status),
                   "post-audit stale-base fixture projects");
  if (!passed)
    return false;
  post_audit_fixture.flux.z.view.unchecked({0, 0, 0}) += 0.25;
  const auto stale_velocity = post_audit_fixture.velocity.storage;
  const auto stale_previous = post_audit_fixture.velocity_previous.storage;
  const auto stale_trial = post_audit_fixture.velocity_trial.storage;
  const auto stale_x = post_audit_fixture.flux.x.storage;
  const auto stale_y = post_audit_fixture.flux.y.storage;
  const auto stale_z = post_audit_fixture.flux.z.storage;
  std::array<FieldView, 3U> stale_layers{
      post_audit_fixture.velocity.view,
      post_audit_fixture.velocity_previous.view,
      post_audit_fixture.velocity_trial.view};
  const Status stale_commit = post_audit_fixture.plan.commit(
      post_audit_candidate, {stale_layers.data(), stale_layers.size()},
      post_audit_fixture.flux.view);
  passed &= expect(
      stale_commit.code == StatusCode::invalid_plan &&
          post_audit_fixture.velocity.storage == stale_velocity &&
          post_audit_fixture.velocity_previous.storage == stale_previous &&
          post_audit_fixture.velocity_trial.storage == stale_trial &&
          post_audit_fixture.flux.x.storage == stale_x &&
          post_audit_fixture.flux.y.storage == stale_y &&
          post_audit_fixture.flux.z.storage == stale_z,
      "base phi mutation after audit rejects before every joint live store");

  Fixture commit_fixture;
  passed &= expect(commit_fixture.initialize(
                       FreshStartProjectionLinearRoute::native_mg_fgmres),
                   "single-use commit fixture compiles");
  FreshStartKinematicProjectionCandidateCertificate candidate;
  if (passed)
    status = commit_fixture.project(candidate);
  passed &=
      expect(static_cast<bool>(status), "single-use commit fixture projects");
  if (!passed)
    return false;
  std::array<FieldView, 3U> layers{commit_fixture.velocity.view,
                                   commit_fixture.velocity_previous.view,
                                   commit_fixture.velocity_trial.view};
  status = commit_fixture.plan.commit(candidate, {layers.data(), layers.size()},
                                      commit_fixture.flux.view);
  passed &=
      expect(static_cast<bool>(status), "first use of joint authority commits");
  const auto velocity_after = commit_fixture.velocity.storage;
  const auto previous_after = commit_fixture.velocity_previous.storage;
  const auto trial_after = commit_fixture.velocity_trial.storage;
  const auto flux_x_after = commit_fixture.flux.x.storage;
  const auto flux_y_after = commit_fixture.flux.y.storage;
  const auto flux_z_after = commit_fixture.flux.z.storage;
  const Status double_commit = commit_fixture.plan.commit(
      candidate, {layers.data(), layers.size()}, commit_fixture.flux.view);
  passed &=
      expect(double_commit.code == StatusCode::invalid_plan &&
                 commit_fixture.velocity.storage == velocity_after &&
                 commit_fixture.velocity_previous.storage == previous_after &&
                 commit_fixture.velocity_trial.storage == trial_after &&
                 commit_fixture.flux.x.storage == flux_x_after &&
                 commit_fixture.flux.y.storage == flux_y_after &&
                 commit_fixture.flux.z.storage == flux_z_after,
             "consumed authority rejects a double commit with zero writes");
  return passed;
}

bool test_alias_and_activity_contracts() {
  Fixture workspace_alias;
  bool passed =
      expect(!workspace_alias.initialize(
                 FreshStartProjectionLinearRoute::jacobi_pcg_oracle, false,
                 false, false, true),
             "compile rejects overlapping rhs/diagonal scratch bindings");
  Fixture face_linear_alias;
  passed &=
      expect(!face_linear_alias.initialize(
                 FreshStartProjectionLinearRoute::jacobi_pcg_oracle, false,
                 false, false, false, {8, 6, 4}, false, false, true),
             "compile rejects face scratch overlapping Krylov storage");
  Fixture asymmetric;
  passed &= expect(
      !asymmetric.initialize(FreshStartProjectionLinearRoute::jacobi_pcg_oracle,
                             false, false, true),
      "compile collectively rejects an active face with an inactive endpoint");
  Fixture input_alias;
  passed &= expect(input_alias.initialize(
                       FreshStartProjectionLinearRoute::jacobi_pcg_oracle),
                   "input-alias fixture compiles");
  if (!passed)
    return false;
  const auto rhs_before = input_alias.rhs.storage;
  FreshStartKinematicProjectionPreparedCertificate prepared;
  const Status rejected = input_alias.plan.prepare(
      {as_const(input_alias.density.view), as_const(input_alias.velocity.view),
       as_const(input_alias.velocity.view),
       as_const(input_alias.velocity_trial.view),
       as_const(input_alias.flux.view), 702U},
      prepared);
  passed &=
      expect(rejected.code == StatusCode::invalid_plan && !prepared.valid() &&
                 input_alias.rhs.storage == rhs_before,
             "prepare rejects overlapping Fresh velocity layers before scratch "
             "writes");
  FaceFluxView aliased_flux = input_alias.flux.view;
  aliased_flux.x.base = input_alias.linear_vectors.view.base;
  FreshStartKinematicProjectionPreparedCertificate face_prepared;
  const Status face_rejected = input_alias.plan.prepare(
      {as_const(input_alias.density.view), as_const(input_alias.velocity.view),
       as_const(input_alias.velocity_previous.view),
       as_const(input_alias.velocity_trial.view), as_const(aliased_flux), 703U},
      face_prepared);
  passed &= expect(face_rejected.code == StatusCode::invalid_plan &&
                       !face_prepared.valid() &&
                       input_alias.rhs.storage == rhs_before,
                   "prepare rejects live face flux overlapping Krylov storage "
                   "before scratch writes");
  return passed;
}

bool test_no_ibm_bitwise_bypass() {
  Fixture fixture;
  bool passed = expect(
      fixture.initialize(FreshStartProjectionLinearRoute::native_mg_fgmres,
                         false, true),
      "no-IBM Fresh bypass fixture compiles");
  if (!passed)
    return false;
  fixture.velocity_previous.view.unchecked({0, 0, 0}, 0U) = -7.25;
  fixture.velocity_trial.view.unchecked({1, 0, 0}, 1U) = 3.5;
  fixture.flux.y.view.unchecked({0, 1, 0}) = -0.75;
  const auto velocity_before = fixture.velocity.storage;
  const auto previous_before = fixture.velocity_previous.storage;
  const auto trial_before = fixture.velocity_trial.storage;
  const auto x_before = fixture.flux.x.storage;
  const auto y_before = fixture.flux.y.storage;
  const auto z_before = fixture.flux.z.storage;
  FreshStartKinematicProjectionCandidateCertificate candidate;
  const Status projected = fixture.project(candidate);
  passed &= expect(static_cast<bool>(projected) && candidate.valid(),
                   "no-IBM bypass still produces auditable joint authority");
  if (!passed)
    return false;
  std::array<FieldView, 3U> layers{fixture.velocity.view,
                                   fixture.velocity_previous.view,
                                   fixture.velocity_trial.view};
  const Status committed = fixture.plan.commit(
      candidate, {layers.data(), layers.size()}, fixture.flux.view);
  passed &= expect(committed && fixture.velocity.storage == velocity_before &&
                       fixture.velocity_previous.storage == previous_before &&
                       fixture.velocity_trial.storage == trial_before &&
                       fixture.flux.x.storage == x_before &&
                       fixture.flux.y.storage == y_before &&
                       fixture.flux.z.storage == z_before,
                   "empty activity takes a literal zero-store bitwise bypass");
  return passed;
}

bool test_component_payload_is_surface_bounded() {
  Fixture small;
  Fixture large;
  bool passed = expect(
      small.initialize(FreshStartProjectionLinearRoute::jacobi_pcg_oracle,
                       false, false, false, false, {8, 6, 4}, true, true),
      "small fragmented local-only component fixture compiles");
  passed &= expect(
      large.initialize(FreshStartProjectionLinearRoute::jacobi_pcg_oracle,
                       false, false, false, false, {16, 6, 4}, true, true),
      "large fragmented local-only component fixture compiles");
  if (!passed)
    return false;
  const auto &small_red = small.plan.red();
  const auto &large_red = large.plan.red();
  passed &= expect(large_red.active_cells > small_red.active_cells &&
                       large_red.component_count > small_red.component_count &&
                       large_red.component_collective_payload_u64 ==
                           small_red.component_collective_payload_u64 &&
                       large_red.component_collective_global_payload_u64 ==
                           small_red.component_collective_global_payload_u64,
                   "communication peak follows the fixed x-boundary area, not "
                   "added interior volume/local fragments");
  return passed;
}

} // namespace

int main(int argc, char **argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    return 2;
  bool passed = test_red_projection_and_oracle();
  passed &= test_component_compatibility_rejects_without_candidate_write();
  passed &= test_tampered_candidate_cannot_partially_commit();
  passed &= test_stale_base_and_single_use_authority();
  passed &= test_alias_and_activity_contracts();
  passed &= test_no_ibm_bitwise_bypass();
  passed &= test_component_payload_is_surface_bounded();
  MPI_Finalize();
  if (!passed)
    return 1;
  std::cout << "fresh-start kinematic projection unit tests passed\n";
  return 0;
}
