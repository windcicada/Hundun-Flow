// SPDX-License-Identifier: Apache-2.0

#include "../support/candidate_boundary_fixture.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition)
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  return condition;
}

bool all_true(bool local, MPI_Comm communicator) {
  const int input = local ? 1 : 0;
  int result = 0;
  return MPI_Allreduce(&input, &result, 1, MPI_INT, MPI_MIN, communicator) ==
             MPI_SUCCESS &&
         result != 0;
}

std::size_t cell_offset(Int3 cells, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(cells.y) *
                  static_cast<std::size_t>(cell.z));
}

std::array<TriangleInput, 12U> domain_clipped_slab() noexcept {
  const Real3 a{0.4, -0.25, -0.25};
  const Real3 b{0.4, -0.25, 1.25};
  const Real3 c{0.4, 1.25, -0.25};
  const Real3 d{0.4, 1.25, 1.25};
  const Real3 e{1.25, -0.25, -0.25};
  const Real3 f{1.25, -0.25, 1.25};
  const Real3 g{1.25, 1.25, -0.25};
  const Real3 h{1.25, 1.25, 1.25};
  return {TriangleInput{a, d, b}, TriangleInput{a, c, d},
          TriangleInput{e, f, h}, TriangleInput{e, h, g},
          TriangleInput{a, b, f}, TriangleInput{a, f, e},
          TriangleInput{c, g, h}, TriangleInput{c, h, d},
          TriangleInput{a, e, g}, TriangleInput{a, g, c},
          TriangleInput{b, d, h}, TriangleInput{b, h, f}};
}

bool test_domain_clipped_solid_owner_flux(MPI_Comm world, int rank) {
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.immersed = true;
  spec.immersed_touches_x_min = true;
  spec.cells_per_axis = 24;
  spec.inlet_velocity = 1.25;
  CandidateBoundaryScratch candidate;
  const bool initialized = fixture.initialize(world, spec);
  const bool staged =
      initialized && fixture.stage(0.0, 16.0, 0.0, 41000U, candidate);
  if (!initialized || !staged)
    std::cerr << "rank " << rank
              << " clipped-inlet step=" << fixture.diagnostic_step << " status="
              << static_cast<unsigned>(fixture.diagnostic_status.code) << '/'
              << fixture.diagnostic_status.detail << '\n';
  bool passed = expect(
      initialized && staged && candidate.final_boundary.valid(), rank,
      "domain-clipped IBM candidate boundary chain compiles and finalizes");
  if (!all_true(passed, world))
    return false;

  const Span<const std::uint8_t> region = fixture.immersed_topology.region();
  std::uint64_t local_solid_owner_faces = 0U;
  std::uint64_t local_nonzero_solid_owner_faces = 0U;
  if (fixture.local_face_owner(CartesianFace::x_min)) {
    for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z) {
      for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
        const Int3 owner{0, y, z};
        if (region.data[cell_offset(fixture.patch.cells, owner)] !=
            static_cast<std::uint8_t>(RegionFlag::solid))
          continue;
        ++local_solid_owner_faces;
        const double flux = candidate.final_flux.x.unchecked({0, y, z});
        if (!std::isfinite(flux) || flux != 0.0)
          ++local_nonzero_solid_owner_faces;
      }
    }
  }
  std::uint64_t global_solid_owner_faces = local_solid_owner_faces;
  std::uint64_t global_nonzero_solid_owner_faces =
      local_nonzero_solid_owner_faces;
  MPI_Allreduce(MPI_IN_PLACE, &global_solid_owner_faces, 1, MPI_UINT64_T,
                MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &global_nonzero_solid_owner_faces, 1,
                MPI_UINT64_T, MPI_SUM, world);
  passed &= expect(global_solid_owner_faces > 0U &&
                       global_nonzero_solid_owner_faces == 0U,
                   rank,
                   "solid-owner physical inlet faces are exact zero after "
                   "final EOS flux closure");
  return all_true(passed, world);
}

bool test_mass_flow_counts_only_active_faces(MPI_Comm world, int rank) {
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.immersed = true;
  spec.immersed_touches_x_min = true;
  spec.cells_per_axis = 24;
  spec.inlet = CandidateBoundaryInlet::mass_flow;
  spec.mass_flow_rate = 0.375;
  CandidateBoundaryScratch candidate;
  bool passed =
      expect(fixture.initialize(world, spec) &&
                 fixture.stage(0.0, 8.0, 0.0, 42000U, candidate) &&
                 candidate.final_boundary.valid(),
             rank, "domain-clipped IBM mass-flow candidate finalizes");
  if (!all_true(passed, world))
    return false;

  const Span<const std::uint8_t> region = fixture.immersed_topology.region();
  std::uint64_t local_active = 0U;
  std::uint64_t local_inactive = 0U;
  double local_achieved = 0.0;
  double local_max_active_error = 0.0;
  if (fixture.local_face_owner(CartesianFace::x_min)) {
    for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
        const bool active =
            region.data[cell_offset(fixture.patch.cells, {0, y, z})] ==
            static_cast<std::uint8_t>(RegionFlag::fluid);
        const double flux = candidate.final_flux.x.unchecked({0, y, z});
        local_achieved += flux;
        local_active += active ? 1U : 0U;
        local_inactive += active ? 0U : 1U;
        if (!active && flux != 0.0)
          local_max_active_error = std::numeric_limits<double>::infinity();
      }
  }
  std::uint64_t global_active = local_active;
  std::uint64_t global_inactive = local_inactive;
  double global_achieved = local_achieved;
  MPI_Allreduce(MPI_IN_PLACE, &global_active, 1, MPI_UINT64_T, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &global_inactive, 1, MPI_UINT64_T, MPI_SUM,
                world);
  MPI_Allreduce(MPI_IN_PLACE, &global_achieved, 1, MPI_DOUBLE, MPI_SUM, world);
  const double expected_active_flux =
      global_active == 0U
          ? 0.0
          : spec.mass_flow_rate / static_cast<double>(global_active);
  if (fixture.local_face_owner(CartesianFace::x_min)) {
    for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
        const bool active =
            region.data[cell_offset(fixture.patch.cells, {0, y, z})] ==
            static_cast<std::uint8_t>(RegionFlag::fluid);
        if (!active)
          continue;
        local_max_active_error =
            std::max(local_max_active_error,
                     std::abs(candidate.final_flux.x.unchecked({0, y, z}) -
                              expected_active_flux));
      }
  }
  double global_max_active_error = local_max_active_error;
  MPI_Allreduce(MPI_IN_PLACE, &global_max_active_error, 1, MPI_DOUBLE, MPI_MAX,
                world);
  passed &= expect(
      global_active > 0U && global_inactive > 0U &&
          std::abs(global_achieved - spec.mass_flow_rate) <=
              128.0 * std::numeric_limits<double>::epsilon() &&
          global_max_active_error <=
              128.0 * std::numeric_limits<double>::epsilon(),
      rank,
      "mass-flow capacity and achieved target contain only active inlet faces");
  return all_true(passed, world);
}

bool test_backflow_skips_inactive_outlet_faces(MPI_Comm world, int rank) {
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.immersed = true;
  spec.immersed_touches_x_min = true;
  spec.reverse_open_boundaries = true;
  spec.cells_per_axis = 24;
  spec.allow_backflow = true;
  spec.backflow_velocity = -0.4;
  CandidateBoundaryScratch candidate;
  const bool initialized = fixture.initialize(world, spec);
  const bool staged =
      initialized && fixture.stage(1.0, -1000.0, 0.0, 43000U, candidate);
  if (!initialized || !staged)
    std::cerr << "rank " << rank
              << " clipped-outlet step=" << fixture.diagnostic_step
              << " status="
              << static_cast<unsigned>(fixture.diagnostic_status.code) << '/'
              << fixture.diagnostic_status.detail << '\n';
  bool passed =
      expect(initialized && staged && candidate.final_boundary.valid(), rank,
             "domain-clipped IBM outlet backflow candidate finalizes");
  if (!all_true(passed, world))
    return false;

  const Span<const std::uint8_t> region = fixture.immersed_topology.region();
  std::uint64_t local_inactive = 0U;
  std::uint64_t local_active_backflow = 0U;
  std::uint64_t local_inactive_nonzero = 0U;
  if (fixture.local_face_owner(CartesianFace::x_min)) {
    const std::int32_t face_x = 0;
    const std::int32_t owner_x = 0;
    for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
        const bool active =
            region.data[cell_offset(fixture.patch.cells, {owner_x, y, z})] ==
            static_cast<std::uint8_t>(RegionFlag::fluid);
        const double flux = candidate.final_flux.x.unchecked({face_x, y, z});
        if (active)
          local_active_backflow += flux > 0.0 ? 1U : 0U;
        else {
          ++local_inactive;
          local_inactive_nonzero += flux == 0.0 ? 0U : 1U;
        }
      }
  }
  std::uint64_t global_inactive = local_inactive;
  std::uint64_t global_active_backflow = local_active_backflow;
  std::uint64_t global_inactive_nonzero = local_inactive_nonzero;
  MPI_Allreduce(MPI_IN_PLACE, &global_inactive, 1, MPI_UINT64_T, MPI_SUM,
                world);
  MPI_Allreduce(MPI_IN_PLACE, &global_active_backflow, 1, MPI_UINT64_T, MPI_SUM,
                world);
  MPI_Allreduce(MPI_IN_PLACE, &global_inactive_nonzero, 1, MPI_UINT64_T,
                MPI_SUM, world);
  if (!(global_inactive > 0U && global_active_backflow > 0U &&
        global_inactive_nonzero == 0U &&
        candidate.final_boundary.outlet_fixed_point_iterations() == 1U))
    std::cerr << "rank " << rank << " outlet inactive=" << global_inactive
              << " active-backflow=" << global_active_backflow
              << " inactive-nonzero=" << global_inactive_nonzero << " fixed="
              << candidate.final_boundary.outlet_fixed_point_iterations()
              << '\n';
  passed &=
      expect(global_inactive > 0U && global_active_backflow > 0U &&
                 global_inactive_nonzero == 0U &&
                 candidate.final_boundary.outlet_fixed_point_iterations() == 1U,
             rank,
             "outlet backflow EOS closure skips inactive faces and preserves "
             "active fixed point");
  return all_true(passed, world);
}

std::vector<double> face_values(ConstFaceFluxView flux) {
  std::vector<double> values;
  for (ConstFaceFieldView face : {flux.x, flux.y, flux.z})
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          values.push_back(face.unchecked({x, y, z}));
  return values;
}

bool test_activity_authority_collective_fail_closed(
    MPI_Comm world, int rank,
    IbmPhysicalBoundaryFluxAuthority& finalize_lifetime) {
  int size = 0;
  MPI_Comm_size(world, &size);
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.immersed = true;
  spec.immersed_touches_x_min = true;
  spec.cells_per_axis = 24;
  CandidateBoundaryScratch candidate;
  bool passed = expect(fixture.initialize(world, spec) &&
                           fixture.stage(0.0, 16.0, 0.0, 44000U, candidate),
                       rank, "collective activity-authority fixture finalizes");
  if (!all_true(passed, world))
    return false;
  const std::vector<double> before =
      face_values(as_const(candidate.final_flux));

  IbmPhysicalBoundaryFluxAuthority preserved_authority;
  passed &=
      expect(static_cast<bool>(IbmPhysicalBoundaryFluxAuthority::compile(
                 world, fixture.geometry, fixture.patch,
                 fixture.immersed_topology, fixture.immersed_interface,
                 preserved_authority)),
             rank, "preserved physical-boundary activity authority compiles");
  if (!all_true(passed, world))
    return false;

  CandidateBoundaryFixture foreign;
  CandidateBoundaryFixtureSpec foreign_spec = spec;
  foreign_spec.immersed_touches_x_min = false;
  passed &= expect(foreign.initialize(world, foreign_spec), rank,
                   "foreign IBM activity fixture compiles");
  if (!all_true(passed, world))
    return false;

  const PlanFingerprint preserved_local =
      preserved_authority.local_fingerprint();
  const PlanFingerprint preserved_collective =
      preserved_authority.collective_fingerprint();
  const Status crossed_plan_compile = IbmPhysicalBoundaryFluxAuthority::compile(
      world, fixture.geometry, fixture.patch, fixture.immersed_topology,
      foreign.immersed_interface, preserved_authority);
  passed &= expect(
      crossed_plan_compile.code == StatusCode::invalid_plan &&
          preserved_authority.local_fingerprint() == preserved_local &&
          preserved_authority.collective_fingerprint() ==
              preserved_collective,
      rank,
      "topology/interface cross-pair rejects and preserves sealed authority");
  if (!all_true(passed, world))
    return false;

  MeshPatch forged_patch = fixture.patch;
  forged_patch.process_coord.x = forged_patch.process_grid.x;
  const Status forged_patch_compile =
      IbmPhysicalBoundaryFluxAuthority::compile(
          world, fixture.geometry, forged_patch, fixture.immersed_topology,
          fixture.immersed_interface, preserved_authority);
  passed &= expect(
      forged_patch_compile.code == StatusCode::invalid_plan &&
          preserved_authority.local_fingerprint() == preserved_local &&
          preserved_authority.collective_fingerprint() ==
              preserved_collective,
      rank,
      "noncanonical patch provenance rejects and preserves sealed authority");
  if (!all_true(passed, world))
    return false;

  if (size > 1) {
    const std::array<TriangleInput, 12U> slab = domain_clipped_slab();
    const StlScanBudget scan_budget{UINT64_C(268435456),
                                    UINT64_C(536870912), UINT64_C(4000000),
                                    UINT64_C(10000), 1U};
    StlScanPlan slab_scan;
    ImmersedSurfacePlan slab_surface;
    passed &= expect(
        static_cast<bool>(StlScanCompiler::compile_triangles(
            fixture.geometry, fixture.patch, {slab.data(), slab.size()},
            CartesianAxis::y, scan_budget, slab_scan)) &&
            static_cast<bool>(
                ImmersedSurfaceCompiler::compile(slab_scan, slab_surface)),
        rank, "domain-clipped slab surface compiles");
    if (!all_true(passed, world))
      return false;

    ImmersedPlanLimits slab_limits;
    ImmersedDomainBoundaryPolicy slab_policy;
    slab_policy.allow_one_sided_quadratic.fill(true);
    EBTopology outside_topology;
    BoundaryStencilPlan outside_boundary;
    IbmEquationInterfacePlan outside_interface;
    EBTopology inside_topology;
    BoundaryStencilPlan inside_boundary;
    IbmEquationInterfacePlan inside_interface;
    const Status outside_topology_status = EBTopologyCompiler::compile(
        world, fixture.geometry, fixture.patch, slab_scan, slab_surface,
        ImmersedFluidSide::outside, slab_limits, outside_topology);
    const Status outside_boundary_status =
        outside_topology_status
            ? BoundaryStencilCompiler::compile(
                  world, fixture.geometry, fixture.patch, slab_surface,
                  outside_topology, slab_policy, slab_limits,
                  outside_boundary)
            : outside_topology_status;
    const Status outside_interface_status =
        outside_boundary_status
            ? IbmEquationInterfacePlan::compile(
                  fixture.equations.kernels(), outside_topology,
                  outside_boundary, outside_topology.interface_metric(),
                  outside_interface)
            : outside_boundary_status;
    const Status inside_topology_status = EBTopologyCompiler::compile(
        world, fixture.geometry, fixture.patch, slab_scan, slab_surface,
        ImmersedFluidSide::inside, slab_limits, inside_topology);
    const Status inside_boundary_status =
        inside_topology_status
            ? BoundaryStencilCompiler::compile(
                  world, fixture.geometry, fixture.patch, slab_surface,
                  inside_topology, slab_policy, slab_limits, inside_boundary)
            : inside_topology_status;
    const Status inside_interface_status =
        inside_boundary_status
            ? IbmEquationInterfacePlan::compile(
                  fixture.equations.kernels(), inside_topology,
                  inside_boundary, inside_topology.interface_metric(),
                  inside_interface)
            : inside_boundary_status;
    passed &= expect(
        outside_interface_status && inside_interface_status &&
            outside_topology.geometry_fingerprint() ==
                inside_topology.geometry_fingerprint() &&
            outside_topology.surface_fingerprint() ==
                inside_topology.surface_fingerprint() &&
            outside_topology.fingerprint() != inside_topology.fingerprint(),
        rank, "same geometry/surface inside/outside slab interfaces compile");
    if (!passed)
      std::cerr << "rank " << rank << " slab status outside="
                << outside_topology_status.detail << '/'
                << outside_boundary_status.detail << '/'
                << outside_interface_status.detail << " inside="
                << inside_topology_status.detail << '/'
                << inside_boundary_status.detail << '/'
                << inside_interface_status.detail << '\n';
    if (!all_true(passed, world))
      return false;

    const bool use_inside = rank == size - 1;
    const Status mixed_side_compile =
        IbmPhysicalBoundaryFluxAuthority::compile(
            world, fixture.geometry, fixture.patch,
            use_inside ? inside_topology : outside_topology,
            use_inside ? inside_interface : outside_interface,
            preserved_authority);
    passed &= expect(
        mixed_side_compile.code == StatusCode::invalid_plan &&
            preserved_authority.local_fingerprint() == preserved_local &&
            preserved_authority.collective_fingerprint() ==
                preserved_collective,
        rank,
        "mixed inside/outside global topology rejects collectively without "
        "overwriting sealed authority");
    if (!all_true(passed, world))
      return false;
  }

  PressureEnergyCandidateBoundaryFinalizerBinding binding =
      fixture.finalizer_binding();
  PressureEnergyCandidateBoundaryFinalizer rejected;
  MPI_Comm congruent = MPI_COMM_NULL;
  MPI_Comm reversed = MPI_COMM_NULL;
  int reversed_relation = MPI_UNEQUAL;
  const int congruent_status = MPI_Comm_dup(world, &congruent);
  const int reversed_status =
      MPI_Comm_split(world, 0, size - 1 - rank, &reversed);
  const int relation_status =
      reversed_status == MPI_SUCCESS
          ? MPI_Comm_compare(world, reversed, &reversed_relation)
          : MPI_ERR_COMM;
  const bool communicator_setup =
      congruent_status == MPI_SUCCESS && reversed_status == MPI_SUCCESS &&
      relation_status == MPI_SUCCESS;
  bool communicator_passed =
      expect(communicator_setup, rank,
             "congruent and reversed communicators compile");
  IbmPhysicalBoundaryFluxAuthority duplicated_authority;
  if (communicator_setup) {
    const Status duplicated_authority_status =
        IbmPhysicalBoundaryFluxAuthority::compile(
            congruent, fixture.geometry, fixture.patch,
            fixture.immersed_topology, fixture.immersed_interface,
            duplicated_authority);
    communicator_passed &= expect(
        duplicated_authority_status &&
            duplicated_authority.matches(
                congruent, &fixture.geometry, fixture.patch,
                &fixture.immersed_interface) &&
            duplicated_authority.matches(
                world, &fixture.geometry, fixture.patch,
                &fixture.immersed_interface) &&
            (size == 1 ||
             !duplicated_authority.matches(
                 reversed, &fixture.geometry, fixture.patch,
                 &fixture.immersed_interface)),
        rank,
        "authority accepts IDENT/CONGRUENT order and rejects SIMILAR order");

    binding.communicator = congruent;
    binding.immersed_physical_boundary_flux = &duplicated_authority;
    PressureEnergyCandidateBoundaryFinalizer congruent_finalizer;
    const Status congruent_bind_status =
        PressureEnergyCandidateBoundaryFinalizer::bind(
            binding, congruent_finalizer);
    communicator_passed &=
        expect(congruent_bind_status && congruent_finalizer.ready(), rank,
               "CONGRUENT communicator preserves authority binding");

    if (size > 1) {
      binding = fixture.finalizer_binding();
      binding.communicator = reversed;
      const Status reversed_bind_status =
          PressureEnergyCandidateBoundaryFinalizer::bind(binding, rejected);
      communicator_passed &= expect(
          reversed_relation == MPI_SIMILAR &&
              reversed_bind_status.code == StatusCode::invalid_plan &&
              !rejected.ready() &&
              face_values(as_const(candidate.final_flux)) == before,
          rank,
          "SIMILAR reversed communicator rejects binding without flux write");

      const Status reversed_compile =
          IbmPhysicalBoundaryFluxAuthority::compile(
              reversed, fixture.geometry, fixture.patch,
              fixture.immersed_topology, fixture.immersed_interface,
              preserved_authority);
      communicator_passed &= expect(
          reversed_compile.code == StatusCode::invalid_plan &&
              preserved_authority.local_fingerprint() == preserved_local &&
              preserved_authority.collective_fingerprint() ==
                  preserved_collective,
          rank,
          "reversed communicator patch provenance preserves old authority");
    }
  }
  if (reversed != MPI_COMM_NULL)
    MPI_Comm_free(&reversed);
  if (congruent != MPI_COMM_NULL)
    MPI_Comm_free(&congruent);
  communicator_passed &= expect(
      duplicated_authority.matches(
          world, &fixture.geometry, fixture.patch,
          &fixture.immersed_interface),
      rank, "authority-owned duplicate survives caller communicator free");
  passed &= communicator_passed;
  if (!all_true(passed, world))
    return false;

  binding = fixture.finalizer_binding();
  binding.immersed_physical_boundary_flux = nullptr;
  Status rejected_status =
      PressureEnergyCandidateBoundaryFinalizer::bind(binding, rejected);
  passed &= expect(rejected_status.code == StatusCode::invalid_plan &&
                       !rejected.ready() &&
                       face_values(as_const(candidate.final_flux)) == before,
                   rank,
                   "IBM interface without activity authority rejects "
                   "collectively without final-flux write");
  if (!all_true(passed, world))
    return false;

  binding = fixture.finalizer_binding();
  binding.immersed_interface = nullptr;
  rejected_status =
      PressureEnergyCandidateBoundaryFinalizer::bind(binding, rejected);
  passed &= expect(rejected_status.code == StatusCode::invalid_plan &&
                       !rejected.ready() &&
                       face_values(as_const(candidate.final_flux)) == before,
                   rank,
                   "IBM activity authority without interface rejects "
                   "collectively without final-flux write");
  if (!all_true(passed, world))
    return false;

  binding = fixture.finalizer_binding();
  binding.immersed_physical_boundary_flux =
      &foreign.immersed_physical_boundary_flux;
  rejected_status =
      PressureEnergyCandidateBoundaryFinalizer::bind(binding, rejected);
  passed &= expect(rejected_status.code == StatusCode::invalid_plan &&
                       !rejected.ready() &&
                       face_values(as_const(candidate.final_flux)) == before,
                   rank,
                   "foreign activity authority rejects collectively without "
                   "final-flux write");
  if (!all_true(passed, world))
    return false;

  if (size > 1) {
    binding = fixture.finalizer_binding();
    if (rank == size - 1)
      binding.immersed_physical_boundary_flux =
          &foreign.immersed_physical_boundary_flux;
    rejected_status =
        PressureEnergyCandidateBoundaryFinalizer::bind(binding, rejected);
    passed &= expect(
        rejected_status.code == StatusCode::invalid_plan && !rejected.ready() &&
            face_values(as_const(candidate.final_flux)) == before,
        rank,
        "cross-rank activity fingerprint mismatch rejects collectively "
        "without final-flux write");
    if (!all_true(passed, world))
      return false;

    const bool use_foreign = rank == size - 1;
    const Status mixed_compile = IbmPhysicalBoundaryFluxAuthority::compile(
        world, use_foreign ? foreign.geometry : fixture.geometry,
        use_foreign ? foreign.patch : fixture.patch,
        use_foreign ? foreign.immersed_topology : fixture.immersed_topology,
        use_foreign ? foreign.immersed_interface : fixture.immersed_interface,
        preserved_authority);
    passed &= expect(
        mixed_compile.code == StatusCode::invalid_plan &&
            preserved_authority.local_fingerprint() == preserved_local &&
            preserved_authority.collective_fingerprint() ==
                preserved_collective,
        rank,
        "cross-rank topology mismatch leaves the prior sealed authority "
        "unchanged");
  }
  const bool globally_passed = all_true(passed, world);
  if (globally_passed)
    finalize_lifetime = std::move(duplicated_authority);
  return globally_passed;
}

} // namespace

int main(int argc, char **argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    return 2;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  IbmPhysicalBoundaryFluxAuthority finalize_lifetime;
  const bool passed =
      test_domain_clipped_solid_owner_flux(MPI_COMM_WORLD, rank) &&
      test_mass_flow_counts_only_active_faces(MPI_COMM_WORLD, rank) &&
      test_backflow_skips_inactive_outlet_faces(MPI_COMM_WORLD, rank) &&
      test_activity_authority_collective_fail_closed(
          MPI_COMM_WORLD, rank, finalize_lifetime);
  MPI_Finalize();
  return passed ? 0 : 1;
}
