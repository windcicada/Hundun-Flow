// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::v04::test {

inline constexpr StlScanBudget kForceScanBudget{
    UINT64_C(268435456), UINT64_C(536870912), UINT64_C(4000000),
    UINT64_C(10000), 1U};

inline CartesianMeshSpec force_mesh(std::int32_t cells = 24) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {-1.5, -1.5, -1.5};
  mesh.upper = {1.5, 1.5, 1.5};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {cells, cells, cells};
  mesh.minimum_spacing = {1.0e-9, 1.0e-9, 1.0e-9};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {UINT64_C(1000000), UINT64_C(1073741824)};
  return mesh;
}

inline std::array<TriangleInput, 12U> force_cube() {
  const Real3 a{-0.5, -0.5, -0.5};
  const Real3 b{-0.5, -0.5, 0.5};
  const Real3 c{-0.5, 0.5, -0.5};
  const Real3 d{-0.5, 0.5, 0.5};
  const Real3 e{0.5, -0.5, -0.5};
  const Real3 f{0.5, -0.5, 0.5};
  const Real3 g{0.5, 0.5, -0.5};
  const Real3 h{0.5, 0.5, 0.5};
  return {TriangleInput{a, d, b}, TriangleInput{a, c, d},
          TriangleInput{e, f, h}, TriangleInput{e, h, g},
          TriangleInput{a, b, f}, TriangleInput{a, f, e},
          TriangleInput{c, g, h}, TriangleInput{c, h, d},
          TriangleInput{a, e, g}, TriangleInput{a, g, c},
          TriangleInput{b, d, h}, TriangleInput{b, h, f}};
}

struct ForceOwnedField {
  std::vector<double> storage;
  FieldView view{};
};

inline ForceOwnedField make_force_field(FieldId field, Int3 cells,
                                        std::uint8_t components,
                                        std::uint8_t ghosts,
                                        RevisionToken revision,
                                        StorageIdentity identity) {
  ForceOwnedField result;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.storage.assign(nx * ny * nz * components, 0.0);
  result.view.base = result.storage.data() + ghosts + ghosts * nx +
                     ghosts * nx * ny;
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = components;
  result.view.stride_y = nx;
  result.view.stride_z = nx * ny;
  result.view.component_stride = nx * ny * nz;
  result.view.field = field;
  result.view.revision = revision;
  result.view.storage_identity = identity;
  result.view.revision_domain = 93001U;
  return result;
}

class CommittedFluxFixture {
 public:
  bool initialize(Int3 cells) {
    FieldId dependency = 0U;
    if (!registry_.declare_field("force.final", 1U, 0U, dependency) ||
        !registry_.freeze(schema_)) {
      return false;
    }
    dependency_ = dependency;
    const std::array requests{ArenaFieldRequest{
        dependency_, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
    if (!ArenaLayout::compile(schema_, {requests.data(), requests.size()},
                              layout_) ||
        !StateLayers::allocate(layout_, layers_) ||
        !AttemptTransaction::create(layers_.field_count(), 1U,
                                    layers_.field_count(), transaction_) ||
        !FaceFluxStorage::allocate_final(cells, storage_) ||
        !authority_.claim(77U, 0U, transaction_, writer_) ||
        !transaction_.begin(layers_) ||
        !transaction_.revise_trial(dependency_)) {
      return false;
    }
    const RevisionDependency dependency_stamp{
        AttemptTransaction::field_revision_source(dependency_),
        transaction_.trial_revision(dependency_)};
    const std::array dependencies{dependency_stamp};
    return writer_.begin_pending(transaction_, storage_, pending_) &&
           writer_.publish_pending(
               {dependencies.data(), dependencies.size()}, pending_) &&
           transaction_.collective_finish(MPI_COMM_SELF, Status{}) &&
           writer_.committed(storage_, flux_);
  }

  ConstFaceFluxView flux() const noexcept { return flux_; }

 private:
  FieldRegistry registry_;
  FieldSchema schema_;
  ArenaLayout layout_;
  StateLayers layers_;
  AttemptTransaction transaction_;
  FaceFluxStorage storage_;
  FinalFaceFluxAuthority authority_;
  FinalFaceFluxWriter writer_;
  PendingFaceFluxView pending_;
  ConstFaceFluxView flux_{};
  FieldId dependency_{};
};

struct IbmForceFixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  StlScanPlan scan;
  ImmersedSurfacePlan surface;
  EBTopology topology;
  BoundaryStencilPlan boundary;
  SurfaceQuadraturePlan quadrature;
  ForceOwnedField velocity;
  ForceOwnedField pressure;
  ForceOwnedField gradient;
  ForceOwnedField viscosity;
  CommittedFluxFixture committed_flux;

  bool initialize(MPI_Comm communicator = MPI_COMM_SELF,
                  std::int32_t cells = 24) {
    const auto triangles = force_cube();
    ImmersedPlanLimits limits;
    if (!CartesianGeometryCompiler::compile(communicator, force_mesh(cells), {},
                                            geometry, patch) ||
        !StlScanCompiler::compile_triangles(
            geometry, patch, {triangles.data(), triangles.size()},
            CartesianAxis::y, kForceScanBudget, scan) ||
        !ImmersedSurfaceCompiler::compile(scan, surface) ||
        !EBTopologyCompiler::compile(
            communicator, geometry, patch, scan, surface,
            ImmersedFluidSide::outside, limits, topology) ||
        !BoundaryStencilCompiler::compile(communicator, geometry, patch,
                                          surface, topology, limits,
                                          boundary) ||
        !SurfaceQuadratureCompiler::compile(communicator, geometry, patch,
                                            surface, topology, limits,
                                            quadrature)) {
      return false;
    }
    const std::uint8_t ghosts =
        quadrature.reconstruction().maximum_halo_reach();
    velocity = make_force_field(0U, patch.cells, 3U, ghosts, 301U, 401U);
    pressure = make_force_field(1U, patch.cells, 1U, ghosts, 302U, 402U);
    gradient = make_force_field(2U, patch.cells, 9U, ghosts, 303U, 403U);
    viscosity = make_force_field(3U, patch.cells, 1U, ghosts, 304U, 404U);
    return committed_flux.initialize(patch.cells);
  }

  double extrapolated_centre(const AxisMetrics& axis,
                             std::int32_t global) const {
    if (global >= 0 &&
        global < static_cast<std::int32_t>(axis.centres().size)) {
      return axis.centres().data[static_cast<std::size_t>(global)];
    }
    const double width = axis.uniform_width();
    return global < 0
               ? axis.centres().data[0U] + width * global
               : axis.centres().data[axis.centres().size - 1U] +
                     width * (global -
                              static_cast<std::int32_t>(axis.centres().size) +
                              1);
  }

  void fill_analytic(double pressure_reference = 101000.0) {
    constexpr VelocityGradient exact_gradient{{
        0.2, 0.3, -0.1, -0.4, -0.05, 0.2, 0.1, -0.15, -0.15}};
    const std::uint8_t ghosts = pressure.view.ghosts.x;
    for (std::int32_t z = -ghosts; z < patch.cells.z + ghosts; ++z) {
      for (std::int32_t y = -ghosts; y < patch.cells.y + ghosts; ++y) {
        for (std::int32_t x = -ghosts; x < patch.cells.x + ghosts; ++x) {
          const Int3 local{x, y, z};
          const Int3 global{patch.begin.x + x, patch.begin.y + y,
                            patch.begin.z + z};
          const double px = extrapolated_centre(geometry.x(), global.x);
          const double py = extrapolated_centre(geometry.y(), global.y);
          const double pz = extrapolated_centre(geometry.z(), global.z);
          const double width = geometry.x().uniform_width();
          const double perturbation =
              2.0 + 0.7 * px - 0.4 * py + 0.2 * pz +
              0.1 * (px * px + py * py + pz * pz +
                     3.0 * width * width / 12.0);
          pressure.view.unchecked(local, 0U) = perturbation;
          viscosity.view.unchecked(local, 0U) = 1.9e-5;
          for (std::uint8_t component = 0U; component < 9U; ++component) {
            gradient.view.unchecked(local, component) =
                exact_gradient.value[component];
          }
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            velocity.view.unchecked(local, component) =
                exact_gradient.value[component * 3U] * px +
                exact_gradient.value[component * 3U + 1U] * py +
                exact_gradient.value[component * 3U + 2U] * pz;
          }
        }
      }
    }
    (void)pressure_reference;
  }

  FinalSurfaceState state(ConstFaceFluxView flux,
                          double pressure_reference = 101000.0) const {
    constexpr PlanFingerprint turbulence_plan = 901U;
    FinalSurfaceState result;
    result.terminal_plan = 902U;
    result.terminal_state = 903U;
    result.final_flux = flux.revision;
    result.face_flux = flux;
    result.final_velocity = as_const(velocity.view);
    result.pressure_perturbation = as_const(pressure.view);
    result.velocity_gradient = as_const(gradient.view);
    result.effective_viscosity = as_const(viscosity.view);
    result.gradient_authority = make_derived_revision_tuple(
        result.final_velocity, geometry, patch, 904U, 905U, 906U,
        turbulence_plan);
    result.turbulence = {turbulence_plan, 907U, 908U,
                         gradient.view.revision, viscosity.view.revision,
                         909U};
    result.geometry = geometry.topology_revision();
    result.pressure_reference = pressure_reference;
    return result;
  }
};

}  // namespace hundun::v04::test
