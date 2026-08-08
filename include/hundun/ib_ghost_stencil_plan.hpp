// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/ib_domain.hpp"
#include "hundun/ib_quadratic_reconstruction.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace hundun::mesh {
class MeshGeometry;
}

namespace hundun::runtime {
class MpiContext;
class StructuredDecomposition;
} // namespace hundun::runtime

namespace hundun::finite_volume {
class ImmersedOperatorAdapter;
}

namespace hundun::immersed {

struct AffineGhostConstraint final {
  ImmersedLinkId link{};
  std::vector<WeightedDonor> donors;
  double wall_value_weight{};
  double wall_normal_gradient_weight_m{};
};

struct FluidExtrapolation final {
  ImmersedLinkId link{};
  std::vector<WeightedDonor> donors;
};

namespace detail {
struct GhostStencilPlanStorage;
struct WallQuadraturePlanStorage;
} // namespace detail

class GhostStencilPlan final {
public:
  static GhostStencilPlan
  create(const ImmersedSurface &, const SurfaceQuery &, const ImmersedDomain &,
         const mesh::MeshTopology &, const mesh::MeshGeometry &,
         const runtime::StructuredDecomposition &, const runtime::MpiContext &);

  const AffineGhostConstraint &velocity_constraint(ImmersedLinkId,
                                                   std::size_t component) const;
  const AffineGhostConstraint &zero_normal_constraint(ImmersedLinkId) const;
  const FluidExtrapolation &density_extrapolation(ImmersedLinkId) const;
  const QuadraticReconstruction &reconstruction(ImmersedLinkId) const;
  std::uint32_t maximum_halo_reach() const noexcept;
  std::uint64_t fingerprint() const noexcept;

private:
  std::size_t immersed_operator_link_count() const noexcept;
  const ImmersedLink &link_for_immersed_operator(ImmersedLinkId) const;
  runtime::Real3
  surface_measure_vector_m2_for_immersed_operator(ImmersedLinkId) const;
  runtime::Real3
  surface_patch_centroid_m_for_immersed_operator(ImmersedLinkId) const;

  explicit GhostStencilPlan(
      std::shared_ptr<const detail::GhostStencilPlanStorage> storage)
      : storage_(std::move(storage)) {}

  std::shared_ptr<const detail::GhostStencilPlanStorage> storage_;
  friend class finite_volume::ImmersedOperatorAdapter;
};

struct WallQuadraturePoint final {
  TriangleId triangle{};
  std::uint32_t point_index{};
  runtime::Real3 position_m{};
  runtime::Real3 solid_to_fluid_normal{};
  double weight_m2{};
  int owner_rank{};
  QuadraticReconstruction reconstruction;
};

class WallQuadraturePlan final {
public:
  static WallQuadraturePlan create(const ImmersedSurface &,
                                   const SurfaceQuery &, const ImmersedDomain &,
                                   const mesh::MeshTopology &,
                                   const mesh::MeshGeometry &,
                                   const runtime::MpiContext &);
  const std::vector<WallQuadraturePoint> &local_points() const noexcept;
  std::uint32_t maximum_halo_reach() const noexcept;
  std::uint64_t fingerprint() const noexcept;

private:
  explicit WallQuadraturePlan(
      std::shared_ptr<const detail::WallQuadraturePlanStorage> storage)
      : storage_(std::move(storage)) {}

  std::shared_ptr<const detail::WallQuadraturePlanStorage> storage_;
};

} // namespace hundun::immersed
