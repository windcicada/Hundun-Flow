// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_linear.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_physics.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::v04 {

namespace detail { class MixtureEnthalpyDiffusion; class MixtureEnthalpyConvection; class IbmScalarTransport; }

class IbmPhysicalBoundaryFluxAuthority;
class ImportedIbmCompiler;
class BoundaryStencilPlan;
class PressureEnergyPressureFluxOperator;
class PressureEnergySchurOperator;
class PressureEnergySharedPressureCertificate;
class PressureEnergySharedPressureInputCertificate;
class SurfaceQuadraturePlan;

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
enum class IbmReconstructionOrder : std::uint8_t {
  linear = 1U,
  quadratic = 2U
};
enum class IbmReconstructionFallbackReason : std::uint8_t {
  none,
  quadratic_donors,
  quadratic_coverage,
  quadratic_rank,
  quadratic_condition
};

struct QuadraticFrame {
  Real3 origin{};
  Real3 normal{};
  Real3 tangent1{};
  Real3 tangent2{};
  double scale{};
  Int3 anchor_global_cell{};
  // Standard interior stencils require all four bits.  A separately
  // certified physical-domain/EB intersection may require the exact subset
  // geometrically reachable inside the domain; it never changes the
  // positive-normal, normal-band, rank, or condition contracts.
  std::uint8_t required_quadrant_mask{0x0fU};
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
  IbmReconstructionPolicy policy{
      IbmReconstructionPolicy::strict_quadratic};
  std::uint8_t minimum_linear_donors{6U};
  std::uint8_t minimum_linear_normal_bands{2U};
  std::uint8_t standard_reach{4U};
};

struct QuadraticStencilQuality {
  std::uint8_t donor_count{};
  std::uint8_t normal_band_count{};
  std::uint8_t quadrant_mask{};
  std::uint8_t required_quadrant_mask{0x0fU};
  std::uint8_t reach{};
  std::uint8_t rank{};
  double condition_estimate{};
  double functional_l1{};
  PlanFingerprint pivot_fingerprint{};
  IbmReconstructionOrder order{IbmReconstructionOrder::quadratic};
  IbmReconstructionFallbackReason fallback_reason{
      IbmReconstructionFallbackReason::none};
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

struct IbmReconstructionAudit {
  bool valid{};
  IbmReconstructionPolicy policy{
      IbmReconstructionPolicy::strict_quadratic};
  std::uint8_t standard_reach{};
  std::uint64_t group_count{};
  std::uint64_t quadratic_groups{};
  std::uint64_t linear_groups{};
  std::uint64_t expanded_search_groups{};
  std::uint64_t rank_fallback_groups{};
  std::uint64_t condition_fallback_groups{};
  std::uint64_t coverage_fallback_groups{};
  std::uint64_t donor_fallback_groups{};
  double maximum_condition_estimate{};
  double maximum_functional_l1{};
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
  // Cold-plan authority for periodic image donors.  A true axis entry means
  // both corresponding physical faces were explicitly validated as a
  // periodic pair by BoundaryStencilCompiler.  Generic standalone quadratic
  // plans never acquire this authority implicitly.
  bool periodic_axis(CartesianAxis axis) const noexcept {
    const std::size_t index = static_cast<std::size_t>(axis);
    return index < periodic_axes_.size() && periodic_axes_[index];
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  const IbmReconstructionAudit& audit() const noexcept { return audit_; }

 private:
  friend class ImportedIbmCompiler;
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
  std::array<bool, 3U> periodic_axes_{};
  IbmReconstructionAudit audit_{};
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

// High-order property reconstruction with a hard positive maximum principle.
// The quadratic value is retained inside the positive donor envelope and is
// projected to that envelope only when the polynomial overshoots.
Status evaluate_positive_bounded_quadratic_row(
    const QuadraticStencilPlan& plan, std::uint32_t row,
    ConstFieldView field, std::uint8_t component, double& out) noexcept;

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
  // Canonical closed-manifold adjacency, sorted by neighbouring triangle id.
  // Sorting makes the graph independent of accepted input winding.
  Span<const std::array<std::uint32_t, 3U>> triangle_neighbours() const
      noexcept {
    return {triangle_neighbours_.data(), triangle_neighbours_.size()};
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
  std::vector<std::array<std::uint32_t, 3U>> triangle_neighbours_;
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
  // Cartesian finite-volume control-face measure.  This is deliberately
  // distinct from the physical STL quadrature assigned by
  // IbmInterfaceMetricPlan.
  double cartesian_control_face_area{};
  Real3 surface_patch_centroid{};
};

// Rank-local binding for one prescribed Cartesian mass flux across an
// immersed fluid--solid link.  face_mass_flux uses FaceFluxView orientation:
// positive values point from the negative-index cell toward the positive one.
// A cold adapter may resolve a global source description into the subset of
// global_link values owned by this rank without exposing case/schema details
// to the equation layer.
struct IbmInterfaceMassFluxSource {
  std::uint64_t global_link{};
  double face_mass_flux{};
};

struct IbmInterfaceInletState {
  std::uint64_t global_link{};
  double face_mass_flux{};
  Real3 velocity{};
  double enthalpy{};
  Span<const double> independent_species{};
};

enum class IbmInterfaceInletFieldKind : std::uint8_t {
  velocity,
  enthalpy,
  independent_species,
  kinetic_energy
};

struct IbmInterfaceInletField {
  IbmInterfaceInletFieldKind kind{IbmInterfaceInletFieldKind::enthalpy};
  std::size_t component{};
};

enum class IbmInterfaceInletEvaluation : std::uint8_t {
  value,
  fixed_state_variation
};

// Physical STL quadrature assigned to one Cartesian fluid--solid link.
// Matrix entries use row-major [spatial_coordinate][normal_component].
// normal_first_moment therefore stores integral(x_i n_j dA), while
// normal_second_moment stores integral(n_i n_j dA).
struct IbmInterfaceLinkMetric {
  std::uint64_t global_link{};
  // Active control-link triangle that receives the nearest-triangle BFS
  // aggregate; this is the partition seed, not a claim that all quadrature
  // originated on one STL triangle.
  SurfaceTriangleId source_triangle{kInvalidSurfaceTriangle};
  double physical_quadrature_area{};
  Real3 physical_area_vector{};   // integral(n dA)
  Real3 physical_first_moment{};  // integral(x dA)
  std::array<double, 9U> normal_first_moment{};
  std::array<double, 9U> normal_second_moment{};
};

struct IbmInterfaceMetricConservation {
  double cartesian_control_area{};
  double physical_quadrature_area{};
  Real3 physical_area_vector{};
  Real3 physical_first_moment{};
  std::array<double, 9U> normal_first_moment{};
  std::array<double, 9U> normal_second_moment{};
};

struct IbmInterfaceMetricResources {
  // Metric-plan storage only.  EBTopologyCompiler separately accounts for
  // the topology storage that remains live while this cold plan is built.
  std::uint64_t persistent_bytes_per_rank{};
  std::uint64_t peak_bytes_per_rank{};
  // Sum of the MPI_DOUBLE element counts passed to metric Allreduces on one
  // rank; status/contract collectives are intentionally excluded.
  std::uint64_t collective_doubles_per_rank{};
};

class IbmInterfaceMetricPlan {
 public:
  IbmInterfaceMetricPlan() noexcept = default;
  IbmInterfaceMetricPlan(const IbmInterfaceMetricPlan&) = delete;
  IbmInterfaceMetricPlan& operator=(const IbmInterfaceMetricPlan&) = delete;
  IbmInterfaceMetricPlan(IbmInterfaceMetricPlan&&) noexcept = default;
  IbmInterfaceMetricPlan& operator=(IbmInterfaceMetricPlan&&) noexcept =
      default;

  Span<const IbmInterfaceLinkMetric> links() const noexcept {
    return {links_.data(), links_.size()};
  }
  const IbmInterfaceMetricConservation& conservation() const noexcept {
    return conservation_;
  }
  IbmInterfaceMetricResources resources() const noexcept { return resources_; }
  RevisionToken geometry_revision() const noexcept {
    return geometry_revision_;
  }
  PlanFingerprint geometry_fingerprint() const noexcept {
    return geometry_fingerprint_;
  }
  PlanFingerprint surface_fingerprint() const noexcept {
    return surface_fingerprint_;
  }
  PlanFingerprint physical_fingerprint() const noexcept {
    return physical_fingerprint_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  friend class ImportedIbmCompiler;
  friend class IbmInterfaceMetricCompiler;
  std::vector<IbmInterfaceLinkMetric> links_;
  IbmInterfaceMetricConservation conservation_{};
  IbmInterfaceMetricResources resources_{};
  RevisionToken geometry_revision_{};
  PlanFingerprint geometry_fingerprint_{};
  PlanFingerprint surface_fingerprint_{};
  PlanFingerprint physical_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

struct ImmersedPlanLimits {
  QuadraticStencilLimits stencil{};
  std::uint64_t maximum_persistent_bytes_per_rank{UINT64_C(1073741824)};
  std::uint64_t maximum_peak_bytes_per_rank{UINT64_C(2147483648)};
  std::uint64_t maximum_local_links{UINT64_C(16777216)};
  std::uint64_t maximum_local_quadrature_points{UINT64_C(16777216)};
};

struct ImmersedDomainBoundaryPolicy {
  // Face order is x_min, x_max, y_min, y_max, z_min, z_max.  A true entry
  // permits a distinct one-sided, full-quadratic intersection stencil at
  // that physical boundary.  Standard interior links remain four-quadrant.
  std::array<bool, 6U> allow_one_sided_quadratic{};
  // Face order is x_min, x_max, y_min, y_max, z_min, z_max.  A true pair
  // authorizes cold periodic image donor wrapping for that axis.  The two
  // faces of an axis must agree; Product derives this array only from a
  // validated paired BoundaryKind::periodic input.
  std::array<bool, 6U> allow_periodic_images{};
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
  const IbmInterfaceMetricPlan& interface_metric() const noexcept {
    return interface_metric_;
  }
  ImmersedFluidSide fluid_side() const noexcept { return fluid_side_; }
  std::uint8_t region_halo_width() const noexcept { return region_halo_width_; }
  bool is_fluid_global(Int3 global_index) const noexcept;
  Int3 global_cells() const noexcept { return global_cells_; }
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
  friend class IbmPhysicalBoundaryFluxAuthority;
  friend class ImportedIbmCompiler;
  friend class EBTopologyCompiler;
  friend class BoundaryStencilCompiler;
  friend class SurfaceQuadratureCompiler;
  friend class IbmInterfaceMetricCompiler;
  std::vector<std::uint8_t> region_;
  std::vector<std::uint8_t> halo_region_;
  std::vector<std::uint32_t> interface_cells_;
  std::vector<ImmersedLink> links_;
  IbmInterfaceMetricPlan interface_metric_;
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

// Compile a Cartesian immersed-boundary authority directly from an immutable,
// globally replicated cell marker (RegionFlag values in x-fastest order).
// This path intentionally does not reinterpret the marker through an STL.
class ImportedIbmCompiler {
 public:
  static Status compile(MPI_Comm communicator,
                        const CartesianGeometryPlan& geometry,
                        const MeshPatch& patch,
                        Span<const std::uint8_t> global_marker,
                        PlanFingerprint marker_source,
                        ImmersedPlanLimits limits, EBTopology& topology,
                        BoundaryStencilPlan& boundary,
                        SurfaceQuadraturePlan& quadrature) noexcept;
};

class IbmInterfaceMetricCompiler {
 public:
  static Status compile(MPI_Comm communicator,
                        const CartesianGeometryPlan& geometry,
                        const MeshPatch& patch,
                        const ImmersedSurfacePlan& surface,
                        const EBTopology& topology,
                        ImmersedPlanLimits limits,
                        IbmInterfaceMetricPlan& out) noexcept;

 private:
  friend class EBTopologyCompiler;
  // Embedded topology compilation keeps its Cartesian arrays resident while
  // building the metric.  Charging that storage here prevents the cold
  // allocation itself from crossing the caller's per-rank peak budget.
  static Status compile_with_resident_storage(
      MPI_Comm communicator, const CartesianGeometryPlan& geometry,
      const MeshPatch& patch, const ImmersedSurfacePlan& surface,
      const EBTopology& topology, ImmersedPlanLimits limits,
      std::uint64_t resident_persistent_bytes,
      IbmInterfaceMetricPlan& out) noexcept;
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
  friend class ImportedIbmCompiler;
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
  static Status compile(MPI_Comm communicator,
                        const CartesianGeometryPlan& geometry,
                        const MeshPatch& patch,
                        const ImmersedSurfacePlan& surface,
                        const EBTopology& topology,
                        ImmersedDomainBoundaryPolicy boundary_policy,
                        ImmersedPlanLimits limits,
                        BoundaryStencilPlan& out) noexcept;
};

struct IbmCellEquationView {
  FieldView diagonal{};
  FieldView rhs{};
  FieldView residual{};
};

// Immutable per-link replacement schedule for Cartesian equations cut by a
// stationary immersed wall.  All hot methods are allocation-free and operate
// on the owning fluid row; shared solid ghost values are never used as
// boundary authority.
// Cold COAST coupling uses fluid-side pressure differences at cut faces;
// the quadratic Neumann reconstruction remains the default for other methods.
enum class IbmPressureGradientKind : std::uint8_t {
  quadratic_neumann, coast_fluid_delta
};

class IbmEquationInterfacePlan {
 public:
  // An unbound/default object rejects Status-returning hot calls with
  // invalid_plan before touching outputs. These methods do not do MPI consensus.
  static Status compile(const CartesianKernelPlan& kernels,
                        const EBTopology& topology,
                        const BoundaryStencilPlan& boundary,
                        const IbmInterfaceMetricPlan& metric,
                        IbmEquationInterfacePlan& out,
                        IbmPressureGradientKind pressure_gradient =
                            IbmPressureGradientKind::quadratic_neumann) noexcept;
  static Status compile(
      const CartesianKernelPlan& kernels, const EBTopology& topology,
      const BoundaryStencilPlan& boundary,
      const IbmInterfaceMetricPlan& metric,
      Span<const IbmInterfaceMassFluxSource> mass_flux_sources,
      IbmEquationInterfacePlan& out,
                        IbmPressureGradientKind pressure_gradient =
                            IbmPressureGradientKind::quadratic_neumann) noexcept;
  static Status compile(
      const CartesianKernelPlan& kernels, const EBTopology& topology,
      const BoundaryStencilPlan& boundary,
      const IbmInterfaceMetricPlan& metric,
      Span<const IbmInterfaceInletState> inlet_states,
      std::size_t independent_species_count,
      IbmEquationInterfacePlan& out,
                        IbmPressureGradientKind pressure_gradient =
                            IbmPressureGradientKind::quadratic_neumann) noexcept;
  static Status compile(const CartesianKernelPlan& kernels,
                        const EBTopology& topology,
                        const BoundaryStencilPlan& boundary,
                        IbmEquationInterfacePlan& out,
                        IbmPressureGradientKind pressure_gradient =
                            IbmPressureGradientKind::quadratic_neumann) noexcept;
  static Status compile(
      const CartesianKernelPlan& kernels, const EBTopology& topology,
      const BoundaryStencilPlan& boundary,
      Span<const IbmInterfaceMassFluxSource> mass_flux_sources,
      IbmEquationInterfacePlan& out,
                        IbmPressureGradientKind pressure_gradient =
                            IbmPressureGradientKind::quadratic_neumann) noexcept;
  static Status compile(
      const CartesianKernelPlan& kernels, const EBTopology& topology,
      const BoundaryStencilPlan& boundary,
      Span<const IbmInterfaceInletState> inlet_states,
      std::size_t independent_species_count,
      IbmEquationInterfacePlan& out,
                        IbmPressureGradientKind pressure_gradient =
                            IbmPressureGradientKind::quadratic_neumann) noexcept;

  Status constrain_interface_flux(FaceFluxView flux) const noexcept;
  // Compatibility spelling for the sealed-wall path.  Source-aware plans
  // enforce their prescribed values after clearing all inactive faces.
  Status zero_interface_flux(FaceFluxView flux) const noexcept;
  Span<const std::uint8_t> cell_activity() const noexcept {
    return topology_ == nullptr ? Span<const std::uint8_t>{} : topology_->region();
  }
  Status validate_interface_flux(ConstFaceFluxView flux,
                                 double absolute_tolerance = 0.0) const
      noexcept;
  // Overwrite only prescribed source links in an unsealed face-value field.
  // Fixed inlet data has exactly zero directional variation.  Callers remain
  // responsible for sealing the resulting reconstruction into their own
  // numeric certificate after this operation succeeds.
  Status override_source_face_values(
      IbmInterfaceInletField field,
      IbmInterfaceInletEvaluation evaluation,
      FrozenConvectionFaceOutput values) const noexcept;
  // Freeze the Cartesian target reconstruction, replace prescribed source
  // links by the fixed inlet state, and reseal that sparse ownership into the
  // numeric certificate.  The resulting E_h face derivative is exactly +0 on
  // every listed source in both generic and compiled paths.
  Status freeze_source_convection_faces(
      IbmInterfaceInletField field, ConvectionScheme scheme,
      ConstFaceFluxView target_flux, ConstFieldView transported,
      std::uint8_t component, FrozenConvectionContext context,
      FrozenConvectionFaceOutput output,
      FrozenConvectionFaceField& frozen) const noexcept;
  Status validate_frozen_source_face_values(
      IbmInterfaceInletField field,
      const FrozenConvectionFaceField& frozen) const noexcept;
  // Returns true only for an exact prescribed source face.  phi is assigned
  // only on success; sealed immersed faces remain absent from this authority.
  bool prescribed_face_flux(CartesianAxis axis, Int3 face,
                            double& phi) const noexcept;
  // Adds only the sparse difference between ordinary Cartesian face
  // reconstruction and the prescribed inlet state.  The Cartesian kernel
  // remains the sole bulk divergence implementation.
  Status add_source_convection_correction(
      IbmInterfaceInletField field, ConvectionScheme scheme,
      ConstFieldView transported, double scale, FieldView output,
      KernelBox box = {}) const noexcept;
  Status add_source_first_order_upwind_correction(
      IbmInterfaceInletField field, ConstFieldView transported, double scale,
      FieldView output, KernelBox box = {}) const noexcept;
  // K is sampled directly from the three-component U stencil; no persisted
  // kinetic-energy field or second reconstruction implementation is needed.
  Status add_source_kinetic_convection_correction(
      ConvectionScheme scheme, ConstFieldView velocity, double scale,
      FieldView output, KernelBox box = {}) const noexcept;
  bool has_inlet_sources() const noexcept {
    return inlet_state_bound_ && !prescribed_interface_fluxes_.empty();
  }
  // Physical work into the fluid at fixed-velocity inlet links, using the
  // same resolved traction as constrain_momentum (stationary walls do no work).
  Status inlet_viscous_work_input(ConstFieldView velocity,
                                 ConstFieldView effective_viscosity,
                                 double& input) const noexcept;
  Status constrain_pressure_predictor(FieldView h_by_a,
                                      FaceFluxView phi_h_by_a) const noexcept;
  Status constrain_corrected_state(FieldView velocity,
                                   FaceFluxView flux) const noexcept;
  Status constrain_momentum(ConstFieldView velocity,
                            ConstFieldView velocity_gradient,
                            ConstFieldView pressure_perturbation,
                            ConstFieldView density,
                            ConstFieldView molecular_viscosity,
                            ConstFieldView effective_viscosity,
                            const TurbulencePlan* wall_treatment,
                            IbmCellEquationView system) const noexcept;
  Status correct_pressure_gradient(ConstFieldView pressure,
                                   FieldView gradient) const noexcept;
  Status correct_pressure_work(ConstFieldView pressure,
                               ConstFieldView velocity,
                               FieldView rate) const noexcept;
  Status correct_velocity_gradient(ConstFieldView velocity,
                                   FieldView velocity_gradient,
                                   Real3 wall_velocity = {}) const noexcept;
  // Add the sparse cut-face difference to Cartesian compatible viscous
  // heating [W/m^3], using exactly the traction in constrain_momentum.
  // total_energy_work selects div(U_face dot tau) instead of internal
  // heating: a stationary wall contributes zero physical boundary work.
  Status correct_viscous_heating(ConstFieldView velocity,
                                 ConstFieldView velocity_gradient,
                                 ConstFieldView density,
                                 ConstFieldView molecular_viscosity,
                                 ConstFieldView effective_viscosity,
                                 const TurbulencePlan* wall_treatment,
                                 FieldView rate, KernelBox box = {},
                                 bool total_energy_work = false) const noexcept;
  // All three diffusion corrections require rate storage to be disjoint
  // from transported and diffusivity, including partial/ghost overlap.
  // The two read-only inputs may share storage.
  Status correct_zero_normal_diffusion(ConstFieldView transported,
                                      ConstFieldView diffusivity,
                                      FieldView rate) const noexcept;
  // Correct an already assembled Cartesian diffusion rate on the binary
  // fluid control volumes: impermeable scalar cut faces carry exactly zero
  // flux. This is distinct from the reconstructed thermal Neumann operator.
  Status correct_impermeable_scalar_diffusion(ConstFieldView transported,
                                             ConstFieldView diffusivity,
                                             FieldView rate) const noexcept;
  Status correct_positive_bounded_zero_normal_diffusion(
      ConstFieldView transported, ConstFieldView diffusivity,
      FieldView rate) const noexcept;
  RevisionToken constrain_certificate(RevisionToken prior,
                                      RevisionToken field,
                                      RevisionToken coefficient) const
      noexcept;
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  friend class detail::MixtureEnthalpyDiffusion;
  friend class detail::MixtureEnthalpyConvection;
  friend class detail::IbmScalarTransport;
  friend class IbmPhysicalBoundaryFluxAuthority;
  struct WallLinearization {
    double distance{};
    double gradient_majorant{};
    double solid_pressure_derivative_weight{};
  };
  struct ViscousBoundaryTraction {
    std::array<double, 3U> residual{};
    Real3 boundary_velocity{};
    double viscosity{};
    double wall_drag_coefficient{};
    bool equilibrium_wall{};
  };
  Status viscous_boundary_traction(
      std::size_t row_index, ConstFieldView velocity, ConstFieldView density,
      ConstFieldView molecular_viscosity, ConstFieldView effective_viscosity,
      const TurbulencePlan* wall_treatment,
      ViscousBoundaryTraction& out) const noexcept;
  struct PrescribedInterfaceFlux {
    std::uint32_t topology_link{};
    double face_mass_flux{};
    Real3 velocity{};
    double enthalpy{};
    std::size_t independent_species_begin{};
    bool has_inlet_state{};
  };
  const PrescribedInterfaceFlux* inlet_for_link(
      std::uint32_t topology_link) const noexcept;
  static Status compile_sources(
      const CartesianKernelPlan& kernels, const EBTopology& topology,
      const BoundaryStencilPlan& boundary,
      const IbmInterfaceMetricPlan& metric,
      Span<const IbmInterfaceMassFluxSource> mass_flux_sources,
      Span<const IbmInterfaceInletState> inlet_states,
      std::size_t independent_species_count, bool inlet_state_bound,
      IbmEquationInterfacePlan& out,
      IbmPressureGradientKind pressure_gradient) noexcept;
  Status add_source_convection_correction_impl(
      IbmInterfaceInletField field, const ConvectionScheme* scheme,
      ConstFieldView transported, double scale, FieldView output,
      KernelBox box, bool kinetic_from_velocity = false) const noexcept;
  const CartesianKernelPlan* kernels_{};
  const EBTopology* topology_{};
  const BoundaryStencilPlan* boundary_{};
  const IbmInterfaceMetricPlan* metric_{};
  std::vector<WallLinearization> wall_linearization_;
  std::vector<std::uint8_t> scalar_blocked_faces_;
  std::vector<PrescribedInterfaceFlux> prescribed_interface_fluxes_;
  std::vector<FrozenConvectionFixedFace> prescribed_source_faces_;
  std::vector<double> prescribed_independent_species_;
  std::size_t independent_species_count_{};
  bool inlet_state_bound_{};
  IbmPressureGradientKind pressure_gradient_{IbmPressureGradientKind::quadratic_neumann};
  PlanFingerprint fingerprint_{};
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
  friend class ImportedIbmCompiler;
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
  static Status compile(MPI_Comm communicator,
                        const CartesianGeometryPlan& geometry,
                        const MeshPatch& patch,
                        const ImmersedSurfacePlan& surface,
                        const EBTopology& topology,
                        ImmersedDomainBoundaryPolicy boundary_policy,
                        ImmersedPlanLimits limits,
                        SurfaceQuadraturePlan& out) noexcept;
};

struct RemoteDonorFieldSpec {
  FieldId field{};
  std::uint8_t components{};
};

// Explicit geometry-certified ghost targets for non-quadratic consumers,
// including tensor parcel stencils. Values are borrowed only during analyze.
struct RemoteDonorTargets {
  PlanFingerprint geometry{};
  std::uint8_t reach{};
  std::array<bool, 3> periodic{};
  Span<const GlobalCellId> global_cells{};
  Span<const Int3> local_indices{};
};

struct RemoteDonorExchangeStats {
  std::uint64_t received_cells{};
  std::uint64_t supplied_cells{};
  std::uint64_t bytes_per_exchange{};
  std::uint32_t peer_messages{};
};

struct RemoteDonorExchangeCounters {
  std::uint64_t exchange_calls{};
  std::uint64_t peer_messages{};
  std::uint64_t bytes{};
  std::uint64_t control_consensus_calls{};
};

// Immutable compact gather for off-rank quadratic donors.  It exchanges only
// the deduplicated global cells referenced by the sealed reconstruction and
// writes them into their certified local ghost coordinates.
class RemoteDonorExchangePlan {
 public:
  RemoteDonorExchangePlan() noexcept = default;
  ~RemoteDonorExchangePlan() noexcept;
  RemoteDonorExchangePlan(const RemoteDonorExchangePlan&) = delete;
  RemoteDonorExchangePlan& operator=(const RemoteDonorExchangePlan&) = delete;
  RemoteDonorExchangePlan(RemoteDonorExchangePlan&&) noexcept;
  RemoteDonorExchangePlan& operator=(RemoteDonorExchangePlan&&) noexcept;

  static Status compile(MPI_Comm communicator,
                        Int3 global_cells,
                        MeshPatch patch,
                        const QuadraticStencilPlan& reconstruction,
                        Span<const RemoteDonorFieldSpec> fields,
                        StageId stage,
                        RemoteDonorExchangePlan& out) noexcept;
  static Status analyze(MPI_Comm communicator,
                        Int3 global_cells,
                        MeshPatch patch,
                        const QuadraticStencilPlan& reconstruction,
                        Span<const RemoteDonorFieldSpec> fields,
                        StageId stage,
                        RemoteDonorExchangePlan& out) noexcept;
  static Status analyze_cells(MPI_Comm communicator, Int3 global_cells,
                              MeshPatch patch, RemoteDonorTargets targets,
                              Span<const RemoteDonorFieldSpec> fields,
                              StageId stage,
                              RemoteDonorExchangePlan &out) noexcept;
  Status bind(MPI_Comm communicator) noexcept;
  // Pure rank-local validation.  Callers that may present rank-dependent
  // views must globalize this status before any rank enters exchange().
  Status preflight_exchange(StageId stage,
                            Span<const FieldView> fields) const noexcept;
  Status exchange(StageId stage, Span<FieldView> fields) noexcept;
  // Collective maximum reach sealed by analyze().  This is deliberately not
  // the rank-local reconstruction reach: consumers that bind a collective
  // candidate path must present one identical capability on every rank.
  std::uint8_t reach() const noexcept;
  RemoteDonorExchangeStats stats() const noexcept;
  RemoteDonorExchangeCounters runtime_counters() const noexcept;
  PlanFingerprint fingerprint() const noexcept;
  bool ready() const noexcept;

 private:
   static Status analyze_values(MPI_Comm communicator, Int3 global_cells,
                                MeshPatch patch, RemoteDonorTargets targets,
                                Span<const RemoteDonorFieldSpec> fields,
                                StageId stage,
                                RemoteDonorExchangePlan &out) noexcept;
   struct Impl;
   void release() noexcept;
   Impl *implementation_{};
};

struct SurfaceTractionPoint {
  Real3 position{};
  Real3 solid_to_fluid_normal{};
  double weight{};
  double absolute_pressure{};
  VelocityGradient velocity_gradient{};
  double effective_viscosity{};
};

struct SurfaceForce {
  Real3 pressure{};
  Real3 viscous{};
  Real3 total{};
  Real3 moment{};
  RevisionToken revision{};
};

struct FinalSurfaceState {
  PlanFingerprint terminal_plan{};
  RevisionToken terminal_state{};
  RevisionToken final_flux{};
  ConstFaceFluxView face_flux{};
  ConstFieldView final_velocity{};
  ConstFieldView pressure_perturbation{};
  ConstFieldView velocity_gradient{};
  ConstFieldView effective_viscosity{};
  DerivedRevisionTuple gradient_authority{};
  TurbulenceCertificate turbulence{};
  RevisionToken geometry{};
  double pressure_reference{};
  Real3 moment_origin{};
};

struct FinalForceCertificate {
  PlanFingerprint plan{};
  RevisionToken terminal_state{};
  RevisionToken final_flux{};
  RevisionToken force{};
  RevisionToken state{};

  bool valid() const noexcept {
    return plan != 0U && terminal_state != 0U && final_flux != 0U &&
           force != 0U && state != 0U;
  }
};

Status integrate_surface_traction(Span<const SurfaceTractionPoint> points,
                                  Real3 moment_origin,
                                  SurfaceForce& out) noexcept;
Status evaluate_surface_force(MPI_Comm communicator,
                              const SurfaceQuadraturePlan& quadrature,
                              const FinalSurfaceState& state,
                              SurfaceForce& out) noexcept;

class FinalForceCache {
 public:
  FinalForceCache() noexcept = default;
  ~FinalForceCache() noexcept;
  FinalForceCache(const FinalForceCache&) = delete;
  FinalForceCache& operator=(const FinalForceCache&) = delete;
  FinalForceCache(FinalForceCache&&) noexcept;
  FinalForceCache& operator=(FinalForceCache&&) noexcept;

  static Status bind(MPI_Comm communicator,
                     const SurfaceQuadraturePlan& quadrature,
                     RevisionToken geometry_revision, StageId stage,
                     RevisionSlotId cache_slot,
                     FinalForceCache& out) noexcept;
  Status prepare(const FinalSurfaceState& state,
                 Span<const RevisionDependency> dependencies,
                 AttemptTransaction& transaction,
                 FinalForceCertificate& certificate) noexcept;
  Status finalize(const AttemptTransaction& transaction) noexcept;
  Status committed(SurfaceForce& force,
                   FinalForceCertificate& certificate) const noexcept;
  PlanFingerprint fingerprint() const noexcept;

 private:
  struct Impl;
  void release() noexcept;
  Impl* implementation_{};
};

class IbmPressureOperator final : public LinearOperator {
 public:
  IbmPressureOperator() noexcept = default;
  ~IbmPressureOperator() noexcept override;
  IbmPressureOperator(const IbmPressureOperator&) = delete;
  IbmPressureOperator& operator=(const IbmPressureOperator&) = delete;
  IbmPressureOperator(IbmPressureOperator&&) noexcept;
  IbmPressureOperator& operator=(IbmPressureOperator&&) noexcept;

  static Status bind(const LinearOperator& regular,
                     const EBTopology& topology,
                     const BoundaryStencilPlan& boundary,
                     ConstFaceFieldView x_coefficient,
                     ConstFaceFieldView y_coefficient,
                     ConstFaceFieldView z_coefficient,
                     RevisionToken geometry_revision,
                     IbmPressureOperator& out) noexcept;
  // Cold-bind against a persistent regular operator before its first numeric
  // refresh.  The exact certificate remains delegated to that operator and
  // therefore becomes current with each pressure-correction lifecycle.
  static Status bind_lifecycle(const LinearOperator& regular,
                               Int3 local_shape,
                               const EBTopology& topology,
                               const BoundaryStencilPlan& boundary,
                               ConstFaceFieldView x_coefficient,
                               ConstFaceFieldView y_coefficient,
                               ConstFaceFieldView z_coefficient,
                               RevisionToken geometry_revision,
                               RemoteDonorExchangePlan* donor_exchange,
                               StageId donor_stage,
                               IbmPressureOperator& out) noexcept;
  LinearOperatorCertificate certificate() const noexcept override;
  Status apply(FieldView x, FieldView y) const noexcept override;
  LinearOperatorFailureProvenance failure_provenance() const
      noexcept override;
  Status mask_solid_rhs(FieldView rhs) const noexcept;
  PlanFingerprint fingerprint() const noexcept;

 private:
  friend class PressureEnergySchurOperator;
  Status certify_pressure_energy_shared_halo(
      const PressureEnergyPressureFluxOperator& energy_pressure,
      PressureEnergySharedPressureCertificate& certificate) const noexcept;
  Status exchange_pressure_energy_shared_input(
      FieldView input,
      const PressureEnergySharedPressureCertificate& certificate,
      PressureEnergySharedPressureInputCertificate& input_certificate) const
      noexcept;
  Status apply_pressure_energy_shared_input(
      FieldView input, FieldView output,
      const PressureEnergySharedPressureCertificate& certificate,
      const PressureEnergySharedPressureInputCertificate& input_certificate)
      const noexcept;
  Status decorate_pressure_action(FieldView input, FieldView output) const
      noexcept;
  static Status bind_internal(const LinearOperator& regular,
                              Int3 local_shape,
                              const EBTopology& topology,
                              const BoundaryStencilPlan& boundary,
                              ConstFaceFieldView x_coefficient,
                              ConstFaceFieldView y_coefficient,
                              ConstFaceFieldView z_coefficient,
                              RevisionToken geometry_revision,
                              RemoteDonorExchangePlan* donor_exchange,
                              StageId donor_stage,
                              IbmPressureOperator& out) noexcept;
  struct Impl;
  void release() noexcept;
  Impl* implementation_{};
};

}  // namespace hundun::v04
