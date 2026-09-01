// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_field.hpp"

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace hundun::v04 {

enum class CartesianAxis : std::uint8_t;
enum class ConvectionScheme : std::uint8_t;
enum class GeometryKind : std::uint8_t;
class BoundaryPlan;
class CartesianGeometryPlan;
class FaceFluxStorage;
class FinalFaceFluxWriter;
class PressureVelocityCoupler;
class ProductDriver;
class SchemePlan;
struct MeshPatch;

using ArenaOwnerId = std::uint16_t;
using RevisionSlotId = std::uint16_t;
using StorageIdentity = std::uint64_t;
using RevisionDomainIdentity = std::uint64_t;
using RevisionSourceId = std::uint32_t;

enum class FieldLifetime : std::uint8_t {
  state_layer,
  pending_cache,
  persistent_workspace,
  step_scratch
};

struct FieldPlacement {
  ArenaOwnerId owner{};
};

struct ArenaFieldRequest {
  FieldId id{};
  Int3 interior{};
  FieldPlacement placement{};
  FieldLifetime lifetime{FieldLifetime::state_layer};
};

struct ArenaFieldLayout {
  FieldId id{};
  Int3 interior{};
  Int3 ghosts{};
  std::uint8_t components{};
  std::size_t stride_y{};
  std::size_t stride_z{};
  std::size_t component_stride{};
  std::size_t replicas{};
  std::size_t replica_stride_doubles{};
  ArenaOwnerId owner{};
  std::size_t owner_index{};
  std::size_t raw_offset_doubles{};
  std::size_t offset_doubles{};
  std::size_t span_doubles{};
  FieldLifetime lifetime{FieldLifetime::state_layer};
};

struct ArenaOwnerLayout {
  ArenaOwnerId owner{};
  std::size_t offset_doubles{};
  std::size_t span_doubles{};
};

class ArenaLayout {
 public:
  static Status compile(const FieldSchema& schema,
                        Span<const ArenaFieldRequest> requests,
                        ArenaLayout& out);

  std::size_t field_count() const noexcept { return fields_.size(); }
  std::size_t owner_count() const noexcept { return owners_.size(); }
  std::size_t total_doubles() const noexcept { return total_doubles_; }
  const ArenaFieldLayout* field(FieldId id) const noexcept;
  const ArenaOwnerLayout* owner(std::size_t index) const noexcept;

 private:
  friend class FieldStorage;
  friend class StateLayers;
  std::vector<ArenaFieldLayout> fields_;
  std::vector<ArenaOwnerLayout> owners_;
  std::size_t total_doubles_{};
};

class RevisionSet {
 public:
  RevisionSet() noexcept = default;
  RevisionSet(const RevisionSet&) = delete;
  RevisionSet& operator=(const RevisionSet&) = delete;
  RevisionSet(RevisionSet&&) noexcept = default;
  RevisionSet& operator=(RevisionSet&&) noexcept = default;

  static Status create(const FieldSchema& schema, RevisionSet& out);
  RevisionToken token(FieldId field) const noexcept;
  Status revise(FieldId field) noexcept;

 private:
  friend class FieldStorage;
  friend class StateLayers;
  friend class AttemptTransaction;
  std::vector<RevisionToken> tokens_;
  RevisionDomainIdentity identity_{};
};

template <class T>
struct BasicFieldView {
  static_assert(std::is_same_v<T, double> ||
                    std::is_same_v<T, const double>,
                "BasicFieldView supports FP64 fields only");

  T* base{};
  Int3 interior{};
  Int3 ghosts{};
  std::uint8_t components{};
  std::size_t stride_y{};
  std::size_t stride_z{};
  std::size_t component_stride{};
  std::size_t replica{};
  FieldId field{};
  RevisionToken revision{};
  StorageIdentity storage_identity{};
  RevisionDomainIdentity revision_domain{};

  T& unchecked(Int3 index, std::uint8_t component) const noexcept {
    const auto x = static_cast<std::ptrdiff_t>(index.x);
    const auto y = static_cast<std::ptrdiff_t>(index.y);
    const auto z = static_cast<std::ptrdiff_t>(index.z);
    const auto sy = static_cast<std::ptrdiff_t>(stride_y);
    const auto sz = static_cast<std::ptrdiff_t>(stride_z);
    const auto sc = static_cast<std::ptrdiff_t>(component_stride);
    const auto c = static_cast<std::ptrdiff_t>(component);
    return base[x + y * sy + z * sz + c * sc];
  }
};

using FieldView = BasicFieldView<double>;
using ConstFieldView = BasicFieldView<const double>;

ConstFieldView as_const(FieldView view) noexcept;

template <class T>
struct BasicFaceFieldView {
  static_assert(std::is_same_v<T, double> ||
                    std::is_same_v<T, const double>,
                "BasicFaceFieldView supports FP64 fields only");

  T* base{};
  Int3 extents{};
  std::size_t stride_y{};
  std::size_t stride_z{};
  CartesianAxis axis{};
  StorageIdentity storage_identity{};
  RevisionDomainIdentity revision_domain{};

  T& unchecked(Int3 index) const noexcept {
    const auto x = static_cast<std::ptrdiff_t>(index.x);
    const auto y = static_cast<std::ptrdiff_t>(index.y);
    const auto z = static_cast<std::ptrdiff_t>(index.z);
    const auto sy = static_cast<std::ptrdiff_t>(stride_y);
    const auto sz = static_cast<std::ptrdiff_t>(stride_z);
    return base[x + y * sy + z * sz];
  }
};

using FaceFieldView = BasicFaceFieldView<double>;
using ConstFaceFieldView = BasicFaceFieldView<const double>;

namespace detail {
struct PendingFaceFluxAccess;
}

struct ConstFaceFluxView;

class FaceFluxCertificate {
 public:
  bool valid() const noexcept {
    return revision_ != 0U && authority_ != 0U && storage_ != 0U &&
           revision_domain_ != 0U && x_base_ != nullptr &&
           y_base_ != nullptr && z_base_ != nullptr && cells_.x > 0 &&
           cells_.y > 0 && cells_.z > 0;
  }
  RevisionToken revision() const noexcept { return revision_; }
  PlanFingerprint authority() const noexcept { return authority_; }
  StorageIdentity storage() const noexcept { return storage_; }
  RevisionDomainIdentity revision_domain() const noexcept {
    return revision_domain_;
  }
  bool matches(ConstFaceFluxView flux) const noexcept;
  friend bool operator==(FaceFluxCertificate left,
                         FaceFluxCertificate right) noexcept {
    return left.revision_ == right.revision_ &&
           left.authority_ == right.authority_ &&
           left.storage_ == right.storage_ &&
           left.revision_domain_ == right.revision_domain_ &&
           left.x_base_ == right.x_base_ && left.y_base_ == right.y_base_ &&
           left.z_base_ == right.z_base_ &&
           left.x_stride_y_ == right.x_stride_y_ &&
           left.x_stride_z_ == right.x_stride_z_ &&
           left.y_stride_y_ == right.y_stride_y_ &&
           left.y_stride_z_ == right.y_stride_z_ &&
           left.z_stride_y_ == right.z_stride_y_ &&
           left.z_stride_z_ == right.z_stride_z_ &&
           left.cells_.x == right.cells_.x &&
           left.cells_.y == right.cells_.y &&
           left.cells_.z == right.cells_.z;
  }

 private:
  friend class FaceFluxConsumer;
  friend class FinalFaceFluxWriter;
  friend class PressureVelocityCoupler;
  RevisionToken revision_{};
  PlanFingerprint authority_{};
  StorageIdentity storage_{};
  RevisionDomainIdentity revision_domain_{};
  const double* x_base_{};
  const double* y_base_{};
  const double* z_base_{};
  std::size_t x_stride_y_{};
  std::size_t x_stride_z_{};
  std::size_t y_stride_y_{};
  std::size_t y_stride_z_{};
  std::size_t z_stride_y_{};
  std::size_t z_stride_z_{};
  Int3 cells_{};
};

struct FaceFluxView {
  FaceFieldView x{};
  FaceFieldView y{};
  FaceFieldView z{};
  RevisionToken revision{};
  FaceFluxCertificate certificate{};
};

struct ConstFaceFluxView {
  ConstFaceFieldView x{};
  ConstFaceFieldView y{};
  ConstFaceFieldView z{};
  RevisionToken revision{};
  FaceFluxCertificate certificate{};
};

class PendingFaceFluxView {
 public:
  PendingFaceFluxView() noexcept = default;
  ~PendingFaceFluxView() noexcept;
  PendingFaceFluxView(const PendingFaceFluxView&) = delete;
  PendingFaceFluxView& operator=(const PendingFaceFluxView&) = delete;
  PendingFaceFluxView(PendingFaceFluxView&&) = delete;
  PendingFaceFluxView& operator=(PendingFaceFluxView&&) = delete;

  RevisionToken revision() const noexcept { return revision_; }
  bool valid() const noexcept { return revision_ != 0U; }

 private:
  friend struct detail::PendingFaceFluxAccess;
  friend class FaceFluxStorage;
  friend class FinalFaceFluxWriter;
  friend class PressureVelocityCoupler;
  FaceFieldView x_{};
  FaceFieldView y_{};
  FaceFieldView z_{};
  RevisionToken revision_{};
  std::uint64_t writer_identity_{};
  std::uint64_t attempt_identity_{};
  FaceFluxStorage* storage_{};
  FinalFaceFluxWriter* writer_{};
};

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace detail {
Status overwrite_pending_face_flux_for_test(PendingFaceFluxView& pending,
                                            double value) noexcept;
}  // namespace detail
#endif

ConstFaceFieldView as_const(FaceFieldView view) noexcept;
ConstFaceFluxView as_const(FaceFluxView view) noexcept;

struct FaceFluxStorageCounters {
  std::uint64_t aligned_payload_allocations{};
  std::uint64_t aligned_payload_bytes{};
  std::uint64_t replicas{};
  std::uint64_t directional_blocks{};
};

class CartesianKernelPlan;

class FaceFluxStorage {
 public:
  FaceFluxStorage() noexcept = default;
  ~FaceFluxStorage() noexcept;
  FaceFluxStorage(const FaceFluxStorage&) = delete;
  FaceFluxStorage& operator=(const FaceFluxStorage&) = delete;
  FaceFluxStorage(FaceFluxStorage&& other) = delete;
  FaceFluxStorage& operator=(FaceFluxStorage&& other) = delete;

  static Status allocate_workspace(Int3 cells, std::size_t replicas,
                                   FaceFluxStorage& out);
  static Status allocate_final(Int3 cells, FaceFluxStorage& out);
  Status workspace_view(std::size_t replica, RevisionToken revision,
                        FaceFluxView& out) noexcept;
  Status view(std::size_t replica, RevisionToken revision,
              ConstFaceFluxView& out) const noexcept;
  FaceFluxStorageCounters counters() const noexcept { return counters_; }

 private:
  friend struct detail::PendingFaceFluxAccess;
  friend class FinalFaceFluxWriter;
  static Status allocate_impl(Int3 cells, std::size_t replicas,
                              bool final_storage, FaceFluxStorage& out);
  void release() noexcept;
  void swap(FaceFluxStorage& other) noexcept;
  Status view_impl(std::size_t replica, RevisionToken revision,
                   FaceFluxView& out) noexcept;
  Status pending_view_impl(std::size_t replica, RevisionToken revision,
                           std::uint64_t writer_identity,
                           std::uint64_t attempt_identity,
                           PendingFaceFluxView& out) noexcept;
  Status view_impl(std::size_t replica, RevisionToken revision,
                   ConstFaceFluxView& out) const noexcept;
  double* data_{};
  Int3 cells_{};
  Int3 extents_[3]{};
  std::size_t stride_y_[3]{};
  std::size_t stride_z_[3]{};
  std::size_t offsets_[3]{};
  std::size_t replica_stride_{};
  std::size_t replicas_{};
  StorageIdentity identity_{};
  RevisionDomainIdentity revision_domain_{};
  FaceFluxStorageCounters counters_{};
  PlanFingerprint authority_identity_{};
  std::uint64_t pending_writer_identity_{};
  std::uint64_t pending_attempt_identity_{};
  FinalFaceFluxWriter* pending_writer_{};
  bool final_storage_{};
};

struct KernelBox {
  Int3 begin{};
  Int3 cells{};
};

struct KernelCounters {
  std::uint64_t invocations{};
  std::uint64_t cells{};
  std::uint64_t faces{};
  std::uint64_t logical_bytes_read{};
  std::uint64_t logical_bytes_written{};
};

struct KernelInvocation {
  Span<const ConstFieldView> reads{};
  Span<const FieldView> writes{};
  KernelBox box{};
  std::uint8_t read_component_begin{};
  std::uint8_t write_component_begin{};
  std::uint8_t component_count{1U};
  RevisionToken required_face_flux_revision{};
  KernelCounters* counters{};
};

namespace detail {

struct CartesianMetricPacket {
  const double* faces{};
  const double* centres{};
  const double* widths{};
  const double* inverse_widths{};
  std::size_t cells{};
  std::int32_t global_begin{};
  double local_face_origin{};
  double local_centre_origin{};
  double uniform_width{};
  double uniform_inverse_width{};
};

}  // namespace detail

class CartesianKernelPlan {
 public:
  CartesianKernelPlan() noexcept = default;
  CartesianKernelPlan(const CartesianKernelPlan&) = delete;
  CartesianKernelPlan& operator=(const CartesianKernelPlan&) = delete;
  CartesianKernelPlan(CartesianKernelPlan&& other) noexcept;
  CartesianKernelPlan& operator=(CartesianKernelPlan&& other) noexcept;

  static Status compile(const SchemePlan& schemes,
                        const CartesianGeometryPlan& geometry,
                        const MeshPatch& patch,
                        const BoundaryPlan& boundary,
                        CartesianKernelPlan& out) noexcept;

  GeometryKind geometry_kind() const noexcept { return geometry_kind_; }
  Int3 cells() const noexcept { return cells_; }
  std::uint8_t reach() const noexcept { return reach_; }
  double limiter() const noexcept { return limiter_; }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  const detail::CartesianMetricPacket& metric(
      std::size_t axis) const noexcept {
    return metrics_[axis];
  }

 private:
  friend Status cartesian_gradient(const CartesianKernelPlan&,
                                   const KernelInvocation&) noexcept;
  friend Status reconstruct_mass_flux(const CartesianKernelPlan&,
                                      const KernelInvocation&,
                                      FaceFluxView&) noexcept;
  friend Status cartesian_face_divergence(const CartesianKernelPlan&,
                                          ConstFaceFluxView,
                                          const KernelInvocation&) noexcept;
  friend Status cartesian_convection(const CartesianKernelPlan&,
                                     ConvectionScheme, ConstFaceFluxView,
                                     const KernelInvocation&) noexcept;
  friend Status cartesian_diffusion(const CartesianKernelPlan&,
                                    ConstFieldView,
                                    const KernelInvocation&) noexcept;
  void reset() noexcept;
  void move_from(CartesianKernelPlan&& other) noexcept;
  void rebind_metrics() noexcept;
  Int3 patch_begin_{};
  Int3 cells_{};
  Int3 process_grid_{};
  Int3 process_coord_{};
  PlanFingerprint boundary_identity_{};
  PlanFingerprint scheme_identity_{};
  PlanFingerprint fingerprint_{};
  GeometryKind geometry_kind_{};
  double limiter_{1.0};
  std::uint8_t reach_{};
  std::vector<double> metric_faces_[3];
  std::vector<double> metric_centres_[3];
  std::vector<double> metric_widths_[3];
  std::vector<double> metric_inverse_widths_[3];
  detail::CartesianMetricPacket metrics_[3]{};
};

Status cartesian_gradient(const CartesianKernelPlan& plan,
                          const KernelInvocation& invocation) noexcept;
Status reconstruct_mass_flux(const CartesianKernelPlan& plan,
                             const KernelInvocation& invocation,
                             FaceFluxView& flux) noexcept;
Status reconstruct_mass_flux(const CartesianKernelPlan& plan,
                             const KernelInvocation& invocation,
                             PendingFaceFluxView& flux) noexcept;
Status cartesian_face_divergence(const CartesianKernelPlan& plan,
                                 ConstFaceFluxView flux,
                                 const KernelInvocation& invocation) noexcept;
Status cartesian_provisional_face_divergence(
    const CartesianKernelPlan& plan, ConstFaceFluxView flux,
    const KernelInvocation& invocation) noexcept;
Status cartesian_convection(const CartesianKernelPlan& plan,
                            ConvectionScheme scheme, ConstFaceFluxView flux,
                            const KernelInvocation& invocation) noexcept;
// Conservative target-layer divergence using an attempt-local, deliberately
// uncommitted flux.  This is distinct from both predictor semantics and the
// final published-flux authority boundary.
Status cartesian_target_convection(
    const CartesianKernelPlan& plan, ConvectionScheme scheme,
    ConstFaceFluxView flux, const KernelInvocation& invocation) noexcept;
Status cartesian_provisional_convection(
    const CartesianKernelPlan& plan, ConvectionScheme scheme,
    ConstFaceFluxView flux, const KernelInvocation& invocation) noexcept;

// Exact face values frozen from the same reconstruction used by
// cartesian_target_convection.  reconstruction is the rank-invariant semantic
// token; revision seals the exact numeric target and may be rank-local.
// local_binding additionally seals the rank-local input and output storage
// identities used to produce these bytes.
struct FrozenConvectionFaceField {
  ConstFaceFieldView x{};
  ConstFaceFieldView y{};
  ConstFaceFieldView z{};
  RevisionToken revision{};
  PlanFingerprint reconstruction{};
  PlanFingerprint local_binding{};
  bool exact_target_reconstruction{};

  bool valid() const noexcept {
    return x.base != nullptr && y.base != nullptr && z.base != nullptr &&
           revision != 0U && reconstruction != 0U && local_binding != 0U &&
           exact_target_reconstruction;
  }
};

struct FrozenConvectionContext {
  // Rank-invariant authority of the equation/scheme/boundary family.
  PlanFingerprint collective_semantics{};
  // Exact numeric closure revision for the current target state; this may be
  // rank-local and therefore never enters a collective contract verbatim.
  RevisionToken closure{};
};

struct FrozenConvectionFaceOutput {
  FaceFieldView x{};
  FaceFieldView y{};
  FaceFieldView z{};
};

enum class FrozenConvectionLinearizationPolicy : std::uint8_t {
  classical_active_branch,
  semismooth_generalized_zero_slope
};

// Opaque, cold-compiled limiter authority for repeated directional actions.
// Each uint16 entry encodes the exact left/right limited-central branch and a
// semismooth-generalized bit.  It is metadata, not a transported FP32 field;
// all numerical values and arithmetic remain FP64.
struct FrozenConvectionBranchOutput {
  Span<std::uint16_t> values{};
};

struct FrozenConvectionBranchPlan {
  Span<const std::uint16_t> values{};
  Int3 cells{};
  PlanFingerprint kernels{};
  RevisionToken revision{};
  PlanFingerprint reconstruction{};
  PlanFingerprint branch_authority{};
  PlanFingerprint local_binding{};
  FrozenConvectionLinearizationPolicy policy{
      FrozenConvectionLinearizationPolicy::classical_active_branch};
  std::uint64_t generalized_face_count{};
  bool classical_everywhere{};

  bool valid() const noexcept {
    return values.data != nullptr && values.size != 0U && cells.x > 0 &&
           cells.y > 0 && cells.z > 0 && kernels != 0U && revision != 0U &&
           reconstruction != 0U &&
           branch_authority != 0U && local_binding != 0U &&
           classical_everywhere == (generalized_face_count == 0U) &&
           (policy == FrozenConvectionLinearizationPolicy::classical_active_branch
                ? classical_everywhere
                : policy == FrozenConvectionLinearizationPolicy::
                                semismooth_generalized_zero_slope);
  }
};

// limited_central2-only optimized route.  Other schemes retain the generic
// derivative path.  Compilation bitwise-checks the frozen nonlinear face field
// before publishing any branch byte.
Status compile_frozen_limited_central2_branches(
    const CartesianKernelPlan& plan, ConstFaceFluxView target_flux,
    ConstFieldView target, std::uint8_t target_component,
    FrozenConvectionContext context,
    FrozenConvectionLinearizationPolicy policy,
    const FrozenConvectionFaceField& frozen,
    FrozenConvectionBranchOutput output,
    FrozenConvectionBranchPlan& branches) noexcept;

Status validate_frozen_limited_central2_branches(
    const CartesianKernelPlan& plan, ConstFaceFluxView target_flux,
    ConstFieldView target, std::uint8_t target_component,
    FrozenConvectionContext context,
    const FrozenConvectionFaceField& frozen,
    const FrozenConvectionBranchPlan& branches) noexcept;

// Applies the original FP64 directional_limited_slope_values arithmetic in
// its original order while consuming the already-certified branch selection.
Status apply_frozen_limited_central2_branches(
    const CartesianKernelPlan& plan,
    const FrozenConvectionBranchPlan& branches, ConstFieldView variation,
    std::uint8_t variation_component,
    FrozenConvectionFaceOutput output) noexcept;

// Two-pass and allocation-free: all face values are preflighted before any
// output byte is changed.  The input flux must be an uncommitted target-layer
// view, matching cartesian_target_convection semantics.
Status freeze_cartesian_target_convection_faces(
    const CartesianKernelPlan& plan, ConvectionScheme scheme,
    ConstFaceFluxView target_flux, ConstFieldView transported,
    std::uint8_t component, FrozenConvectionContext context,
    FrozenConvectionFaceOutput output,
    FrozenConvectionFaceField& frozen) noexcept;

// Directional face values obtained by applying one explicit linearization
// policy to a FrozenConvectionFaceField and a cell-centred variation.
// reconstruction seals only rank-invariant scheme/context/policy semantics.
// revision, branch_authority, local_binding, counts, and all views are
// rank-local.  This certificate intentionally makes no blanket "exact"
// claim: classical_everywhere identifies the ordinary derivative case, while
// generalized_face_count records semismooth generalized faces.
struct FrozenConvectionFaceDirectionalDerivative {
  ConstFaceFieldView x{};
  ConstFaceFieldView y{};
  ConstFaceFieldView z{};
  RevisionToken revision{};
  PlanFingerprint reconstruction{};
  PlanFingerprint branch_authority{};
  PlanFingerprint local_binding{};
  FrozenConvectionLinearizationPolicy policy{
      FrozenConvectionLinearizationPolicy::classical_active_branch};
  std::uint64_t generalized_face_count{};
  bool classical_everywhere{};

  bool valid() const noexcept {
    return x.base != nullptr && y.base != nullptr && z.base != nullptr &&
           revision != 0U && reconstruction != 0U &&
           branch_authority != 0U && local_binding != 0U &&
           classical_everywhere == (generalized_face_count == 0U) &&
           (policy ==
                FrozenConvectionLinearizationPolicy::classical_active_branch
                ? classical_everywhere
                : policy == FrozenConvectionLinearizationPolicy::
                                semismooth_generalized_zero_slope);
  }
};

// Exactly two face traversals and allocation-free.  The first traversal
// preflights values, bitwise-compares the nonlinear target reconstruction with
// frozen, and seals branch authority/counts; the second commits output.  The
// caller's revision discipline must prevent concurrent raw mutation between
// traversals because no face-field lease exists at this seam.  A repeated E_h
// Krylov hot path that cannot afford branch discovery on every application
// must introduce separately owned, pre-registered branch-mask storage rather
// than weakening this numeric comparison.
Status differentiate_frozen_cartesian_target_convection_faces(
    const CartesianKernelPlan& plan, ConvectionScheme scheme,
    ConstFaceFluxView target_flux, ConstFieldView target,
    std::uint8_t target_component, FrozenConvectionContext context,
    FrozenConvectionLinearizationPolicy policy,
    const FrozenConvectionFaceField& frozen, ConstFieldView variation,
    std::uint8_t variation_component, FrozenConvectionFaceOutput output,
    FrozenConvectionFaceDirectionalDerivative& derivative) noexcept;

struct ConvectionPointDiagnostic {
  double divergence{};
  double mass_divergence{};
  double maximum_face_envelope_violation{};
  double selected_face_value{};
  double selected_donor_minimum{};
  double selected_donor_maximum{};
  bool face_envelope_checked{};
  bool face_envelope_valid{};
};

// Failure-path observer for one already-computed Cartesian convection cell.
// It executes no collective and allocates no storage. The returned divergence
// uses the production reconstruction; its adjacent-cell donor envelope is an
// independent check on the TVD face value.
Status diagnose_cartesian_convection_point(
    const CartesianKernelPlan& plan, ConvectionScheme scheme,
    ConstFaceFluxView flux, ConstFieldView transported,
    std::uint8_t component, Int3 cell,
    ConvectionPointDiagnostic& diagnostic) noexcept;

Status cartesian_diffusion(const CartesianKernelPlan& plan,
                           ConstFieldView diffusivity,
                           const KernelInvocation& invocation) noexcept;

enum class StateVisibility : std::uint8_t {
  accepted,
  trial,
  pending,
  committed_snapshot,
  workspace
};

enum class StageKind : std::uint8_t { compute, service, commit };

enum class GraphNodeKind : std::uint8_t {
  compute,
  halo_begin,
  compute_interior,
  halo_finish,
  compute_boundary,
  service,
  collective_consensus,
  commit
};

struct FieldAccessSpec {
  FieldId field{};
  StateVisibility visibility{StateVisibility::accepted};
  friend bool operator==(FieldAccessSpec left,
                         FieldAccessSpec right) noexcept {
    return left.field == right.field &&
           left.visibility == right.visibility;
  }
};

struct GraphFieldSpec {
  FieldAccessSpec access{};
  std::uint8_t ghost_capacity{};
  bool initially_available{};
};

struct StageResourceSpec {
  std::uint64_t merged_halo_messages{};
  std::uint64_t merged_halo_bytes{};
  std::uint64_t numeric_refills{};
  std::uint64_t hierarchy_rebuilds{};
  std::uint64_t cache_publishes{};
  std::uint64_t linear_iterations{};
  std::uint64_t stage_wall_nanoseconds{};
};

struct StageSpec {
  StageId id{};
  Span<const FieldAccessSpec> reads{};
  Span<const FieldAccessSpec> writes{};
  Span<const FieldAccessSpec> ghosts{};
  Span<const std::uint8_t> ghost_widths{};
  Span<const FieldAccessSpec> invalidates{};
  std::size_t workspace_bytes{};
  std::size_t workspace_alignment{64U};
  StageId workspace_live_through{};
  std::size_t fixed_workspace_offset{};
  StageResourceSpec resources{};
  StageKind kind{StageKind::compute};
  bool has_fixed_workspace_offset{};
  bool collective_consensus{};
};

struct ResourceContract {
  std::uint64_t max_live_workspace_bytes{};
  std::uint64_t allocation_allowance{};
  std::uint64_t merged_halo_messages{};
  std::uint64_t merged_halo_bytes{};
  std::uint64_t numeric_refills{};
  std::uint64_t hierarchy_rebuilds{};
  std::uint64_t cache_publishes{};
  std::uint64_t linear_iterations{};
  std::uint64_t stage_wall_nanoseconds{};
};

struct ResourceCounters {
  std::uint64_t peak_workspace_bytes{};
  std::uint64_t allocations{};
  std::uint64_t merged_halo_messages{};
  std::uint64_t merged_halo_bytes{};
  std::uint64_t numeric_refills{};
  std::uint64_t hierarchy_rebuilds{};
  std::uint64_t cache_publishes{};
  std::uint64_t linear_iterations{};
  std::uint64_t stage_wall_nanoseconds{};
};

Status add_resource_counters(ResourceCounters& counters,
                             ResourceCounters increment) noexcept;
Status validate_resource_counters(const ResourceContract& contract,
                                  const ResourceCounters& counters) noexcept;

struct FrozenStage {
  StageId id{};
  StageKind kind{StageKind::compute};
  std::uint32_t registration_ordinal{};
  std::uint32_t node_begin{};
  std::uint16_t node_count{};
  std::uint32_t read_begin{};
  std::uint16_t read_count{};
  std::uint32_t write_begin{};
  std::uint16_t write_count{};
  std::uint32_t ghost_begin{};
  std::uint16_t ghost_count{};
  std::uint32_t invalidation_begin{};
  std::uint16_t invalidation_count{};
  std::size_t workspace_offset{};
  std::size_t workspace_bytes{};
  std::size_t workspace_alignment{64U};
  StageId workspace_live_through{};
  StageResourceSpec resources{};
  bool collective_consensus{};
};

struct GraphNode {
  StageId stage{};
  GraphNodeKind kind{GraphNodeKind::compute};
  std::uint32_t ordinal{};
};

struct GraphEdge {
  std::uint32_t from{};
  std::uint32_t to{};
  friend bool operator==(GraphEdge left, GraphEdge right) noexcept {
    return left.from == right.from && left.to == right.to;
  }
};

class FrozenExecutionGraph {
 public:
  FrozenExecutionGraph() noexcept = default;
  FrozenExecutionGraph(const FrozenExecutionGraph&) = delete;
  FrozenExecutionGraph& operator=(const FrozenExecutionGraph&) = delete;
  FrozenExecutionGraph(FrozenExecutionGraph&&) noexcept = default;
  FrozenExecutionGraph& operator=(FrozenExecutionGraph&&) noexcept = default;

  Span<const FrozenStage> stages() const noexcept {
    return {stages_.data(), stages_.size()};
  }
  Span<const GraphNode> nodes() const noexcept {
    return {nodes_.data(), nodes_.size()};
  }
  Span<const GraphEdge> edges() const noexcept {
    return {edges_.data(), edges_.size()};
  }
  Span<const FieldAccessSpec> reads(StageId stage) const noexcept;
  Span<const FieldAccessSpec> writes(StageId stage) const noexcept;
  Span<const FieldAccessSpec> ghosts(StageId stage) const noexcept;
  Span<const std::uint8_t> ghost_widths(StageId stage) const noexcept;
  Span<const FieldAccessSpec> invalidations(StageId stage) const noexcept;
  const FrozenStage* stage(StageId id) const noexcept;
  const ResourceContract& resources() const noexcept { return resources_; }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  friend class ExecutionGraphCompiler;
  std::vector<FrozenStage> stages_;
  std::vector<GraphNode> nodes_;
  std::vector<GraphEdge> edges_;
  std::vector<FieldAccessSpec> reads_;
  std::vector<FieldAccessSpec> writes_;
  std::vector<FieldAccessSpec> ghosts_;
  std::vector<std::uint8_t> ghost_widths_;
  std::vector<FieldAccessSpec> invalidations_;
  ResourceContract resources_{};
  PlanFingerprint fingerprint_{};
};

class ExecutionGraphCompiler {
 public:
  ExecutionGraphCompiler() noexcept = default;
  ExecutionGraphCompiler(const ExecutionGraphCompiler&) = delete;
  ExecutionGraphCompiler& operator=(const ExecutionGraphCompiler&) = delete;
  ExecutionGraphCompiler(ExecutionGraphCompiler&&) noexcept = default;
  ExecutionGraphCompiler& operator=(ExecutionGraphCompiler&&) noexcept =
      default;

  Status configure(Span<const GraphFieldSpec> fields) noexcept;
  Status register_stage(const StageSpec& stage) noexcept;
  Status freeze(FrozenExecutionGraph& out) noexcept;
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
  Status freeze_for_test(FrozenExecutionGraph& out) noexcept;
#endif
  bool frozen() const noexcept { return frozen_; }

 private:
  struct OwnedStageSpec {
    StageId id{};
    std::vector<FieldAccessSpec> reads;
    std::vector<FieldAccessSpec> writes;
    std::vector<FieldAccessSpec> ghosts;
    std::vector<std::uint8_t> ghost_widths;
    std::vector<FieldAccessSpec> invalidates;
    std::size_t workspace_bytes{};
    std::size_t workspace_alignment{64U};
    StageId workspace_live_through{};
    std::size_t fixed_workspace_offset{};
    StageResourceSpec resources{};
    StageKind kind{StageKind::compute};
    bool has_fixed_workspace_offset{};
    bool collective_consensus{};
  };

  std::vector<GraphFieldSpec> fields_;
  std::vector<OwnedStageSpec> stages_;
  bool configured_{};
  bool frozen_{};
};

struct ExecutionCounters {
  std::uint64_t aligned_payload_allocations{};
  std::uint64_t aligned_payload_bytes{};
  std::uint64_t whole_field_copies{};
};

class FieldStorage {
 public:
  FieldStorage() noexcept = default;
  ~FieldStorage() noexcept;
  FieldStorage(const FieldStorage&) = delete;
  FieldStorage& operator=(const FieldStorage&) = delete;
  FieldStorage(FieldStorage&& other) noexcept;
  FieldStorage& operator=(FieldStorage&& other) noexcept;

  static Status allocate(const ArenaLayout& layout, FieldStorage& out);
  Status view(FieldId field, std::size_t replica,
              const RevisionSet& revisions, FieldView& out) noexcept;
  Status view(FieldId field, std::size_t replica,
              const RevisionSet& revisions, ConstFieldView& out) const noexcept;
  double* checked_ptr(const FieldView& view, Int3 index,
                      std::uint8_t component,
                      const RevisionSet& revisions) noexcept;
  const double* checked_ptr(const ConstFieldView& view, Int3 index,
                            std::uint8_t component,
                            const RevisionSet& revisions) const noexcept;
  ExecutionCounters counters() const noexcept { return counters_; }

 private:
  friend class StateLayers;
  struct OwnerAllocation {
    ArenaOwnerId owner{};
    double* data{};
    std::size_t doubles{};
  };

  void release() noexcept;
  const ArenaFieldLayout* field_layout(FieldId field) const noexcept;
  Status view(FieldId field, std::size_t replica, RevisionToken revision,
              RevisionDomainIdentity revision_domain,
              FieldView& out) noexcept;
  Status view(FieldId field, std::size_t replica, RevisionToken revision,
              RevisionDomainIdentity revision_domain,
              ConstFieldView& out) const noexcept;
  double* checked_ptr(const FieldView& view, Int3 index,
                      std::uint8_t component,
                      RevisionToken expected_revision,
                      RevisionDomainIdentity expected_domain) noexcept;
  const double* checked_ptr(const ConstFieldView& view, Int3 index,
                            std::uint8_t component,
                            RevisionToken expected_revision,
                            RevisionDomainIdentity expected_domain) const noexcept;
  template <class T>
  T* checked_ptr_impl(const BasicFieldView<T>& view, Int3 index,
                      std::uint8_t component,
                      RevisionToken expected_revision,
                      RevisionDomainIdentity expected_domain) const noexcept;

  std::vector<ArenaFieldLayout> fields_;
  std::vector<OwnerAllocation> owners_;
  ExecutionCounters counters_{};
  StorageIdentity identity_{};
};

enum class StateRole : std::uint8_t {
  accepted_n,
  accepted_n_minus_one,
  trial
};

struct PendingCacheStamp {
  RevisionToken cache_revision{};
};

struct RevisionDependency {
  RevisionSourceId source{};
  RevisionToken revision{};
};

class StateLayers {
 public:
  StateLayers() noexcept = default;
  StateLayers(const StateLayers&) = delete;
  StateLayers& operator=(const StateLayers&) = delete;
  StateLayers(StateLayers&& other) noexcept;
  StateLayers& operator=(StateLayers&& other) noexcept;

  static Status allocate(const ArenaLayout& layout, StateLayers& out);
  Status view(StateRole role, FieldId field, FieldView& out) noexcept;
  Status view(StateRole role, FieldId field, ConstFieldView& out) const noexcept;
  Status runtime_view(FieldLifetime lifetime, FieldId field,
                      FieldView& out) noexcept;
  Status runtime_view(FieldLifetime lifetime, FieldId field,
                      ConstFieldView& out) const noexcept;
  RevisionToken runtime_revision(FieldLifetime lifetime,
                                 FieldId field) const noexcept;
  Status revise_runtime(FieldLifetime lifetime, FieldId field) noexcept;
  std::size_t handle(StateRole role) const noexcept;
  RevisionToken revision(StateRole role, FieldId field) const noexcept;
  RevisionToken state_revision(StateRole role) const noexcept;
  std::size_t field_count() const noexcept { return field_count_; }
  double* checked_ptr(StateRole role, const FieldView& view, Int3 index,
                      std::uint8_t component) noexcept;
  const double* checked_ptr(StateRole role, const ConstFieldView& view,
                            Int3 index,
                            std::uint8_t component) const noexcept;
  ExecutionCounters counters() const noexcept { return storage_.counters(); }

 private:
  friend class AttemptTransaction;
  static constexpr std::size_t kRoleCount = 3U;
  std::size_t role_index(StateRole role) const noexcept;
  StorageIdentity storage_identity() const noexcept {
    return storage_.identity_;
  }
  void rotate_commit() noexcept;
  Status issue_revision(std::size_t replica, FieldId field) noexcept;

  FieldStorage storage_;
  std::vector<RevisionToken> revisions_;
  RevisionToken state_revisions_[kRoleCount]{};
  std::vector<std::uint8_t> state_fields_;
  std::size_t field_count_{};
  std::size_t role_handles_[kRoleCount]{0U, 1U, 2U};
  RevisionToken next_revision_{1U};
};

class AttemptTransaction;

enum class AttemptFinishDecision : std::uint8_t { accept, reject };

// The public default constructor only creates an invalid placeholder.  A
// consumable certificate can be produced solely by collective_prepare and is
// bound to one exact transaction attempt and mutation revision.
class PreparedAttemptFinish {
 public:
  PreparedAttemptFinish() noexcept = default;
  PreparedAttemptFinish(const PreparedAttemptFinish&) = delete;
  PreparedAttemptFinish& operator=(const PreparedAttemptFinish&) = delete;
  PreparedAttemptFinish(PreparedAttemptFinish&&) = delete;
  PreparedAttemptFinish& operator=(PreparedAttemptFinish&&) = delete;

  bool valid() const noexcept { return valid_ && !consumed_; }
  AttemptFinishDecision decision() const noexcept { return decision_; }
  Status outcome() const noexcept { return outcome_; }
  int lowest_failing_rank() const noexcept { return lowest_failing_rank_; }

 private:
  friend class AttemptTransaction;
  AttemptTransaction* owner_{};
  StateLayers* layers_{};
  FinalFaceFluxWriter* final_face_flux_writer_{};
  StorageIdentity bound_layers_identity_{};
  std::uint64_t attempt_identity_{};
  std::uint64_t mutation_revision_{};
  std::uint64_t final_face_flux_writer_identity_{};
  Status outcome_{StatusCode::invalid_plan, 0U};
  int lowest_failing_rank_{-1};
  AttemptFinishDecision decision_{AttemptFinishDecision::reject};
  bool locally_active_{};
  bool valid_{};
  bool consumed_{};
};

class AttemptTransaction {
 public:
  AttemptTransaction() noexcept = default;
  ~AttemptTransaction() noexcept;
  AttemptTransaction(const AttemptTransaction&) = delete;
  AttemptTransaction& operator=(const AttemptTransaction&) = delete;
  AttemptTransaction(AttemptTransaction&&) = delete;
  AttemptTransaction& operator=(AttemptTransaction&&) = delete;

  static Status create(std::size_t field_capacity,
                       std::size_t revision_slot_capacity,
                       std::size_t revision_source_capacity,
                       AttemptTransaction& out);
  Status begin(StateLayers& layers) noexcept;
  Status revise_trial(FieldId field) noexcept;
  RevisionToken trial_revision(FieldId field) const noexcept;
  static RevisionSourceId field_revision_source(FieldId field) noexcept {
    return static_cast<RevisionSourceId>(field) + 1U;
  }
  Status bind_dependency(RevisionDependency dependency) noexcept;
  Status publish_pending_cache(RevisionSlotId slot,
                               Span<const RevisionDependency> dependencies,
                               PendingCacheStamp stamp) noexcept;
  Status collective_prepare(MPI_Comm communicator, Status local_status,
                            PreparedAttemptFinish& out) noexcept;
  void commit_accept(PreparedAttemptFinish& prepared) noexcept;
  // commit_reject also accepts a prepared accept certificate.  This is the
  // no-fail abort path when another module's prepare phase cannot commit.
  void commit_reject(PreparedAttemptFinish& prepared) noexcept;
  Status collective_finish(MPI_Comm communicator,
                           Status local_status) noexcept;

  bool active() const noexcept { return active_; }
  bool finished() const noexcept { return finished_; }
  bool committed() const noexcept { return committed_; }
  std::uint64_t attempt_identity() const noexcept {
    return attempt_identity_;
  }
  int lowest_failing_rank() const noexcept { return lowest_failing_rank_; }
  bool pending_cache_valid(RevisionSlotId slot) const noexcept;
  RevisionToken pending_cache(RevisionSlotId slot) const noexcept;

 private:
  friend class FinalFaceFluxAuthority;
  friend class FinalFaceFluxWriter;
  Status claim_final_face_flux_writer(RevisionSlotId slot,
                                      std::uint64_t writer_identity,
                                      FinalFaceFluxWriter& writer) noexcept;
  Status publish_final_face_flux_cache(
      RevisionSlotId slot, Span<const RevisionDependency> dependencies,
      PendingCacheStamp stamp, std::uint64_t writer_identity) noexcept;
  Status preflight_final_face_flux_cache(
      RevisionSlotId slot, Span<const RevisionDependency> dependencies,
      PendingCacheStamp stamp,
      std::uint64_t writer_identity) const noexcept;
  Status publish_pending_cache_impl(
      RevisionSlotId slot, Span<const RevisionDependency> dependencies,
      PendingCacheStamp stamp, std::uint64_t writer_identity) noexcept;
  Status preflight_pending_cache_impl(
      RevisionSlotId slot, Span<const RevisionDependency> dependencies,
      PendingCacheStamp stamp,
      std::uint64_t writer_identity) const noexcept;
  void detach_final_face_flux_writer(FinalFaceFluxWriter& writer,
                                     bool fail_active) noexcept;
  void invalidate_final_face_flux_writer(
      const FinalFaceFluxWriter& writer) noexcept;
  void advance_mutation_revision() noexcept;
  bool matches_prepared(const PreparedAttemptFinish& prepared) const noexcept;
  void discard_attempt() noexcept;
  StateLayers* layers_{};
  std::vector<RevisionToken> active_caches_;
  std::vector<RevisionToken> pending_caches_;
  std::vector<RevisionToken> pending_cache_dependency_revisions_;
  std::vector<RevisionSourceId> pending_cache_dependency_sources_;
  std::vector<std::size_t> pending_cache_dependency_counts_;
  std::vector<RevisionDependency> dependency_catalog_;
  std::size_t dependency_count_{};
  std::size_t dependency_capacity_per_slot_{};
  std::vector<RevisionToken> cache_commit_buffer_;
  std::vector<std::uint8_t> required_pending_caches_;
  std::vector<std::uint64_t> pending_cache_writer_identities_;
  std::vector<std::uint8_t> revised_fields_;
  StorageIdentity bound_layers_identity_{};
  Status attempt_status_{};
  std::uint64_t attempt_identity_{};
  std::uint64_t mutation_revision_{};
  std::uint64_t final_face_flux_writer_identity_{};
  FinalFaceFluxWriter* final_face_flux_writer_{};
  bool active_{};
  bool finished_{};
  bool committed_{};
  int lowest_failing_rank_{-1};
};

class FinalFaceFluxWriter {
 public:
  FinalFaceFluxWriter() noexcept = default;
  ~FinalFaceFluxWriter() noexcept;
  FinalFaceFluxWriter(const FinalFaceFluxWriter&) = delete;
  FinalFaceFluxWriter& operator=(const FinalFaceFluxWriter&) = delete;

  Status begin_pending(AttemptTransaction& transaction,
                       FaceFluxStorage& storage,
                       PendingFaceFluxView& out) noexcept;
  Status initialize_committed(FaceFluxStorage& storage,
                              ConstFaceFluxView source) noexcept;
  // Restart images already represent both accepted time levels while the
  // controller is forced through backward-Euler recovery.  Initialize the
  // active and previous handles directly from those bytes; no synthetic
  // fresh-start flux is inserted into the restored lineage.
  Status initialize_restored(FaceFluxStorage& storage,
                             ConstFaceFluxView source) noexcept;
  // Exact Restart installs distinct t_n/t_{n-1} payloads and preserves their
  // committed revision lineage.  The one-source overload is the legacy BE
  // recovery path and deliberately aliases the two history handles.
  Status initialize_restored_history(
      FaceFluxStorage& storage, ConstFaceFluxView accepted,
      ConstFaceFluxView previous) noexcept;
  // Pure, allocation-free validation for a later publish_pending call.  A
  // failed preflight neither poisons the attempt nor consumes the pending
  // lease, so callers may first reach WORLD consensus and then publish.
  Status preflight_publish_pending(
      Span<const RevisionDependency> dependencies,
      const PendingFaceFluxView& pending) const noexcept;
  Status publish_pending(Span<const RevisionDependency> dependencies,
                         PendingFaceFluxView& pending) noexcept;
  Status committed(const FaceFluxStorage& storage,
                   ConstFaceFluxView& out) const noexcept;
  Status committed_previous(const FaceFluxStorage& storage,
                            ConstFaceFluxView& out) const noexcept;

 private:
  friend class AttemptTransaction;
  friend class FinalFaceFluxAuthority;
  friend class FaceFluxStorage;
  friend class PendingFaceFluxView;
  friend class ProductDriver;
  Status restore_committed(FaceFluxStorage& storage,
                           ConstFaceFluxView source) noexcept;
  bool ready_for_collective(const AttemptTransaction& transaction) const
      noexcept;
  void complete_from_transaction(const AttemptTransaction& transaction,
                                 bool committed) noexcept;
  void detach_from_transaction(AttemptTransaction& transaction) noexcept;
  void detach_from_storage(FaceFluxStorage& storage) noexcept;
  void abandon_pending_view(PendingFaceFluxView& pending) noexcept;
  void invalidate_pending_view() noexcept;
  void clear_storage_lease() noexcept;
  StageId stage_{};
  RevisionSlotId cache_slot_{};
  PlanFingerprint authority_fingerprint_{};
  RevisionToken active_revision_{};
  RevisionToken previous_revision_{};
  RevisionToken pending_revision_{};
  std::uint64_t pending_attempt_identity_{};
  std::size_t active_replica_{};
  std::size_t previous_replica_{};
  std::size_t pending_replica_{1U};
  AttemptTransaction* authority_transaction_{};
  AttemptTransaction* transaction_{};
  FaceFluxStorage* storage_{};
  PendingFaceFluxView* pending_view_{};
  StorageIdentity bound_storage_identity_{};
  RevisionDomainIdentity bound_revision_domain_{};
  Int3 bound_cells_{};
  bool issued_{};
  bool pending_{};
  bool pending_published_{};
};

class FinalFaceFluxAuthority {
 public:
  Status claim(StageId stage, RevisionSlotId cache_slot,
               AttemptTransaction& transaction,
               FinalFaceFluxWriter& out) noexcept;
  bool claimed() const noexcept { return claimed_; }

 private:
  PlanFingerprint fingerprint_{};
  bool claimed_{};
};

class FaceFluxConsumer {
 public:
  Status bind(FaceFluxCertificate certificate) noexcept;
  Status consume(ConstFaceFluxView flux) noexcept;
  RevisionToken consumed_revision() const noexcept {
    return consumed_revision_;
  }

 private:
  FaceFluxCertificate required_certificate_{};
  RevisionToken consumed_revision_{};
};

}  // namespace hundun::v04
