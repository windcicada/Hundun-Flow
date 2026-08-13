// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_execution.hpp"
#include "hundun/v04_mesh.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::v04 {

using GlobalCellId = std::uint64_t;
using SurfaceTriangleId = std::uint64_t;

inline constexpr std::uint32_t kInvalidIbmIndex = UINT32_MAX;
inline constexpr SurfaceTriangleId kInvalidSurfaceTriangle = UINT64_MAX;

enum class QuadraticConstraint : std::uint8_t {
  none,
  origin_value,
  origin_normal_gradient
};
enum class QuadraticFunctionalKind : std::uint8_t {
  value,
  directional_derivative
};

struct QuadraticFrame {
  Real3 origin{};
  Real3 normal{};
  Real3 tangent1{};
  Real3 tangent2{};
  double scale{};
  Int3 anchor_global_cell{};
};

struct QuadraticDonorCell {
  GlobalCellId global_cell{};
  Int3 global_index{};
  // Index relative to the owning rank's first owned cell.  Negative and
  // beyond-interior coordinates address the already registered halo.
  Int3 local_index{};
  Real3 centre{};
  Real3 widths{};
};

struct QuadraticFunctionalRequest {
  QuadraticFunctionalKind kind{QuadraticFunctionalKind::value};
  QuadraticConstraint constraint{QuadraticConstraint::none};
  Real3 point{};
  Real3 direction{};
};

struct QuadraticStencilRequest {
  QuadraticFrame frame{};
  Span<const QuadraticDonorCell> donors{};
  Span<const QuadraticFunctionalRequest> functionals{};
};

struct QuadraticStencilLimits {
  std::uint8_t minimum_donors{14U};
  std::uint8_t maximum_donors{32U};
  std::uint8_t maximum_reach{4U};
  std::uint8_t minimum_normal_bands{3U};
  double condition_limit{1.0e8};
};

struct QuadraticStencilQuality {
  std::uint8_t donor_count{};
  std::uint8_t normal_band_count{};
  std::uint8_t quadrant_mask{};
  std::uint8_t reach{};
  std::uint8_t rank{};
  double condition_estimate{};
  double functional_l1{};
  PlanFingerprint pivot_fingerprint{};
};

struct QuadraticStencilGroup {
  std::uint32_t donor_begin{};
  std::uint32_t row_begin{};
  std::uint8_t row_count{};
  QuadraticStencilQuality quality{};
  PlanFingerprint fingerprint{};
};

struct QuadraticAffineRow {
  std::uint32_t group{};
  std::uint32_t weight_begin{};
  double wall_value_weight{};
  double wall_normal_gradient_weight_m{};
};

class QuadraticStencilPlan {
 public:
  QuadraticStencilPlan() noexcept = default;
  QuadraticStencilPlan(const QuadraticStencilPlan&) = delete;
  QuadraticStencilPlan& operator=(const QuadraticStencilPlan&) = delete;
  QuadraticStencilPlan(QuadraticStencilPlan&&) noexcept = default;
  QuadraticStencilPlan& operator=(QuadraticStencilPlan&&) noexcept = default;

  Span<const QuadraticStencilGroup> groups() const noexcept {
    return {groups_.data(), groups_.size()};
  }
  Span<const QuadraticAffineRow> rows() const noexcept {
    return {rows_.data(), rows_.size()};
  }
  Span<const GlobalCellId> donor_global_cells() const noexcept {
    return {donor_global_cells_.data(), donor_global_cells_.size()};
  }
  Span<const Int3> donor_local_indices() const noexcept {
    return {donor_local_indices_.data(), donor_local_indices_.size()};
  }
  Span<const double> weights() const noexcept {
    return {weights_.data(), weights_.size()};
  }
  std::uint8_t maximum_halo_reach() const noexcept {
    return maximum_halo_reach_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  friend class QuadraticStencilCompiler;
  friend class BoundaryStencilCompiler;
  friend class SurfaceQuadratureCompiler;
  void refresh_fingerprint() noexcept;
  std::vector<QuadraticStencilGroup> groups_;
  std::vector<QuadraticAffineRow> rows_;
  std::vector<GlobalCellId> donor_global_cells_;
  std::vector<Int3> donor_local_indices_;
  std::vector<double> weights_;
  std::uint8_t maximum_halo_reach_{};
  PlanFingerprint fingerprint_{};
};

class QuadraticStencilCompiler {
 public:
  static Status compile(Span<const QuadraticStencilRequest> requests,
                        QuadraticStencilLimits limits,
                        QuadraticStencilPlan& out) noexcept;
};

Status evaluate_quadratic_row(const QuadraticStencilPlan& plan,
                              std::uint32_t row, ConstFieldView field,
                              std::uint8_t component, double wall_value,
                              double wall_normal_gradient,
                              double& out) noexcept;

struct SurfaceTriangle {
  SurfaceTriangleId id{};
  std::array<Real3, 3U> vertices{};
  Real3 geometric_outward_normal{};
  Real3 centroid{};
  double area{};
};

struct ClosestSurfacePoint {
  SurfaceTriangleId triangle{kInvalidSurfaceTriangle};
  Real3 point{};
  Real3 geometric_outward_normal{};
  double squared_distance{};
};

struct SurfaceSegmentIntersection {
  SurfaceTriangleId triangle{kInvalidSurfaceTriangle};
  Real3 point{};
  double segment_fraction{};
};

class ImmersedSurfacePlan {
 public:
  ImmersedSurfacePlan() noexcept = default;
  ImmersedSurfacePlan(const ImmersedSurfacePlan&) = delete;
  ImmersedSurfacePlan& operator=(const ImmersedSurfacePlan&) = delete;
  ImmersedSurfacePlan(ImmersedSurfacePlan&&) noexcept = default;
  ImmersedSurfacePlan& operator=(ImmersedSurfacePlan&&) noexcept = default;

  Span<const SurfaceTriangle> triangles() const noexcept {
    return {triangles_.data(), triangles_.size()};
  }
  Real3 bounding_box_min() const noexcept { return bounding_box_min_; }
  Real3 bounding_box_max() const noexcept { return bounding_box_max_; }
  double closed_volume() const noexcept { return closed_volume_; }
  PlanFingerprint source_triangle_fingerprint() const noexcept {
    return source_triangle_fingerprint_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

  Status closest_point(Real3 point, ClosestSurfacePoint& out) const noexcept;
  Status first_segment_intersection(
      Real3 begin, Real3 end,
      SurfaceSegmentIntersection& out) const noexcept;

 private:
  struct BvhNode {
    Real3 lower{};
    Real3 upper{};
    std::uint32_t begin{};
    std::uint32_t count{};
    std::uint32_t left{kInvalidIbmIndex};
    std::uint32_t right{kInvalidIbmIndex};
  };
  friend class ImmersedSurfaceCompiler;
  std::vector<SurfaceTriangle> triangles_;
  std::vector<std::uint32_t> bvh_triangle_ids_;
  std::vector<BvhNode> bvh_nodes_;
  Real3 bounding_box_min_{};
  Real3 bounding_box_max_{};
  double closed_volume_{};
  PlanFingerprint source_triangle_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

class ImmersedSurfaceCompiler {
 public:
  static Status compile(const StlScanPlan& scan,
                        ImmersedSurfacePlan& out) noexcept;
};

enum class ImmersedFaceDirection : std::uint8_t {
  x_negative,
  x_positive,
  y_negative,
  y_positive,
  z_negative,
  z_positive
};

struct ImmersedLink {
  std::uint64_t global_link{};
  GlobalCellId fluid_cell{};
  GlobalCellId solid_cell{};
  Int3 fluid_global_index{};
  Int3 solid_global_index{};
  Int3 fluid_local_index{};
  Int3 solid_local_index{};
  ImmersedFaceDirection direction{ImmersedFaceDirection::x_negative};
  SurfaceTriangleId triangle{kInvalidSurfaceTriangle};
  Real3 wall_point{};
  Real3 solid_to_fluid_normal{};
  Real3 surface_measure_vector{};
  Real3 surface_patch_centroid{};
};

struct ImmersedPlanLimits {
  QuadraticStencilLimits stencil{};
  std::uint64_t maximum_persistent_bytes_per_rank{UINT64_C(1073741824)};
  std::uint64_t maximum_peak_bytes_per_rank{UINT64_C(2147483648)};
  std::uint64_t maximum_local_links{UINT64_C(16777216)};
  std::uint64_t maximum_local_quadrature_points{UINT64_C(16777216)};
};

class EBTopology {
 public:
  EBTopology() noexcept = default;
  EBTopology(const EBTopology&) = delete;
  EBTopology& operator=(const EBTopology&) = delete;
  EBTopology(EBTopology&&) noexcept = default;
  EBTopology& operator=(EBTopology&&) noexcept = default;

  Span<const std::uint8_t> region() const noexcept {
    return {region_.data(), region_.size()};
  }
  Span<const std::uint32_t> interface_cells() const noexcept {
    return {interface_cells_.data(), interface_cells_.size()};
  }
  Span<const ImmersedLink> links() const noexcept {
    return {links_.data(), links_.size()};
  }
  ImmersedFluidSide fluid_side() const noexcept { return fluid_side_; }
  std::uint8_t region_halo_width() const noexcept { return region_halo_width_; }
  bool is_fluid_global(Int3 global_index) const noexcept;
  RevisionToken geometry_revision() const noexcept {
    return geometry_revision_;
  }
  PlanFingerprint geometry_fingerprint() const noexcept {
    return geometry_fingerprint_;
  }
  PlanFingerprint surface_fingerprint() const noexcept {
    return surface_fingerprint_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  int lowest_failing_rank() const noexcept { return lowest_failing_rank_; }

 private:
  friend class EBTopologyCompiler;
  friend class BoundaryStencilCompiler;
  friend class SurfaceQuadratureCompiler;
  std::vector<std::uint8_t> region_;
  std::vector<std::uint8_t> halo_region_;
  std::vector<std::uint32_t> interface_cells_;
  std::vector<ImmersedLink> links_;
  Int3 global_cells_{};
  MeshPatch patch_{};
  std::size_t halo_stride_y_{};
  std::size_t halo_stride_z_{};
  ImmersedFluidSide fluid_side_{ImmersedFluidSide::outside};
  std::uint8_t region_halo_width_{};
  RevisionToken geometry_revision_{};
  PlanFingerprint geometry_fingerprint_{};
  PlanFingerprint surface_fingerprint_{};
  PlanFingerprint fingerprint_{};
  int lowest_failing_rank_{-1};
};

class EBTopologyCompiler {
 public:
  static Status compile(MPI_Comm communicator,
                        const CartesianGeometryPlan& geometry,
                        const MeshPatch& patch, const StlScanPlan& scan,
                        const ImmersedSurfacePlan& surface,
                        ImmersedFluidSide fluid_side,
                        ImmersedPlanLimits limits,
                        EBTopology& out) noexcept;
};

struct BoundaryStencilLink {
  std::uint32_t topology_link{};
  std::uint32_t reconstruction_group{};
  std::uint32_t dirichlet_value_row{};
  std::uint32_t zero_normal_value_row{};
  std::uint32_t wall_value_row{};
  std::uint32_t wall_normal_gradient_row{};
};

class BoundaryStencilPlan {
 public:
  BoundaryStencilPlan() noexcept = default;
  BoundaryStencilPlan(const BoundaryStencilPlan&) = delete;
  BoundaryStencilPlan& operator=(const BoundaryStencilPlan&) = delete;
  BoundaryStencilPlan(BoundaryStencilPlan&&) noexcept = default;
  BoundaryStencilPlan& operator=(BoundaryStencilPlan&&) noexcept = default;

  Span<const BoundaryStencilLink> links() const noexcept {
    return {links_.data(), links_.size()};
  }
  const QuadraticStencilPlan& reconstruction() const noexcept {
    return reconstruction_;
  }
  std::uint8_t maximum_halo_reach() const noexcept {
    return reconstruction_.maximum_halo_reach();
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  int lowest_failing_rank() const noexcept { return lowest_failing_rank_; }

 private:
  friend class BoundaryStencilCompiler;
  std::vector<BoundaryStencilLink> links_;
  QuadraticStencilPlan reconstruction_;
  PlanFingerprint fingerprint_{};
  int lowest_failing_rank_{-1};
};

class BoundaryStencilCompiler {
 public:
  static Status compile(MPI_Comm communicator,
                        const CartesianGeometryPlan& geometry,
                        const MeshPatch& patch,
                        const ImmersedSurfacePlan& surface,
                        const EBTopology& topology,
                        ImmersedPlanLimits limits,
                        BoundaryStencilPlan& out) noexcept;
};

struct SurfaceQuadraturePoint {
  SurfaceTriangleId triangle{kInvalidSurfaceTriangle};
  std::uint8_t point_index{};
  Real3 position{};
  Real3 solid_to_fluid_normal{};
  double weight{};
  int owner_rank{-1};
  GlobalCellId owner_cell{};
  std::uint32_t reconstruction_group{kInvalidIbmIndex};
  std::uint32_t wall_value_row{kInvalidIbmIndex};
  std::uint32_t wall_normal_gradient_row{kInvalidIbmIndex};
};

class SurfaceQuadraturePlan {
 public:
  SurfaceQuadraturePlan() noexcept = default;
  SurfaceQuadraturePlan(const SurfaceQuadraturePlan&) = delete;
  SurfaceQuadraturePlan& operator=(const SurfaceQuadraturePlan&) = delete;
  SurfaceQuadraturePlan(SurfaceQuadraturePlan&&) noexcept = default;
  SurfaceQuadraturePlan& operator=(SurfaceQuadraturePlan&&) noexcept = default;

  Span<const SurfaceQuadraturePoint> local_points() const noexcept {
    return {local_points_.data(), local_points_.size()};
  }
  const QuadraticStencilPlan& reconstruction() const noexcept {
    return reconstruction_;
  }
  std::uint64_t global_point_count() const noexcept {
    return global_point_count_;
  }
  PlanFingerprint physical_fingerprint() const noexcept {
    return physical_fingerprint_;
  }
  PlanFingerprint local_layout_fingerprint() const noexcept {
    return local_layout_fingerprint_;
  }
  int lowest_failing_rank() const noexcept { return lowest_failing_rank_; }

 private:
  friend class SurfaceQuadratureCompiler;
  std::vector<SurfaceQuadraturePoint> local_points_;
  QuadraticStencilPlan reconstruction_;
  std::uint64_t global_point_count_{};
  PlanFingerprint physical_fingerprint_{};
  PlanFingerprint local_layout_fingerprint_{};
  int lowest_failing_rank_{-1};
};

class SurfaceQuadratureCompiler {
 public:
  static Status compile(MPI_Comm communicator,
                        const CartesianGeometryPlan& geometry,
                        const MeshPatch& patch,
                        const ImmersedSurfacePlan& surface,
                        const EBTopology& topology,
                        ImmersedPlanLimits limits,
                        SurfaceQuadraturePlan& out) noexcept;
};

}  // namespace hundun::v04
