// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_execution.hpp"

#include "hundun/v04_boundary.hpp"
#include "hundun/v04_mesh.hpp"

#include "core_arena_detail.hpp"
#include "solver_cartesian_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kFaceStorage = 911U;
constexpr std::uint32_t kFaceExtent = 912U;
constexpr std::uint32_t kFaceAuthority = 913U;
constexpr std::uint32_t kFaceTransaction = 914U;
constexpr std::uint32_t kFaceKernel = 915U;
constexpr std::uint32_t kFaceNumerical = 916U;
constexpr std::size_t kFinalFaceFluxReplicas = 3U;

bool checked_add(std::size_t left, std::size_t right,
                 std::size_t& out) noexcept {
  return detail::checked_add(left, right, out);
}

bool checked_multiply(std::size_t left, std::size_t right,
                      std::size_t& out) noexcept {
  return detail::checked_multiply(left, right, out);
}

bool checked_align(std::size_t value, std::size_t& out) noexcept {
  return detail::checked_align(value, detail::kDoublesPerCacheLine, out);
}

template <class T>
BasicFaceFieldView<T> make_axis_view(
    T* data, Int3 extents, std::size_t stride_y, std::size_t stride_z,
    CartesianAxis axis, StorageIdentity storage_identity,
    RevisionDomainIdentity revision_domain) noexcept {
  return {data, extents, stride_y, stride_z, axis, storage_identity,
          revision_domain};
}

inline double cell_product(ConstFieldView density, ConstFieldView velocity,
                           Int3 index,
                           std::uint8_t velocity_component) noexcept {
  return density.unchecked(index, 0U) *
         velocity.unchecked(index, velocity_component);
}

template <std::size_t Axis>
inline std::int32_t axis_index(Int3 value) noexcept {
  if constexpr (Axis == 0U) {
    return value.x;
  } else if constexpr (Axis == 1U) {
    return value.y;
  } else {
    return value.z;
  }
}

template <std::size_t Axis>
inline Int3 set_axis(Int3 value, std::int32_t selected) noexcept {
  if constexpr (Axis == 0U) {
    value.x = selected;
  } else if constexpr (Axis == 1U) {
    value.y = selected;
  } else {
    value.z = selected;
  }
  return value;
}

template <bool Uniform, std::size_t Axis>
Status reconstruct_axis(const CartesianKernelPlan& plan,
                        ConstFieldView density, ConstFieldView velocity,
                        KernelBox box, FaceFieldView output) noexcept {
  constexpr auto velocity_component = static_cast<std::uint8_t>(Axis);
  const Int3 cell_end{box.begin.x + box.cells.x,
                      box.begin.y + box.cells.y,
                      box.begin.z + box.cells.z};
  Int3 face_begin = box.begin;
  Int3 face_end = cell_end;
  if constexpr (Axis == 0U) {
    face_begin.x = box.begin.x == 0 ? 0 : box.begin.x + 1;
    face_end.x = cell_end.x + 1;
  } else if constexpr (Axis == 1U) {
    face_begin.y = box.begin.y == 0 ? 0 : box.begin.y + 1;
    face_end.y = cell_end.y + 1;
  } else {
    face_begin.z = box.begin.z == 0 ? 0 : box.begin.z + 1;
    face_end.z = cell_end.z + 1;
  }
  for (std::int32_t z = face_begin.z; z < face_end.z; ++z) {
    for (std::int32_t y = face_begin.y; y < face_end.y; ++y) {
      for (std::int32_t x = face_begin.x; x < face_end.x; ++x) {
        const Int3 face{x, y, z};
        const double rate =
            detail::metric_interpolate_face<Uniform>(
                plan, Axis, axis_index<Axis>(face),
                cell_product(density, velocity,
                             set_axis<Axis>(face,
                                            axis_index<Axis>(face) - 1),
                             velocity_component),
                cell_product(density, velocity, face, velocity_component)) *
            detail::metric_face_area<Uniform>(plan, Axis, face);
        if (!std::isfinite(rate)) {
          return {StatusCode::numerical_failure, kFaceNumerical};
        }
        output.unchecked(face) = rate;
      }
    }
  }
  return {};
}

bool checked_face_count(KernelBox box, std::uint64_t& out) noexcept {
  const auto cx = static_cast<std::uint64_t>(box.cells.x);
  const auto cy = static_cast<std::uint64_t>(box.cells.y);
  const auto cz = static_cast<std::uint64_t>(box.cells.z);
  const std::uint64_t bx = box.begin.x == 0 ? 1U : 0U;
  const std::uint64_t by = box.begin.y == 0 ? 1U : 0U;
  const std::uint64_t bz = box.begin.z == 0 ? 1U : 0U;
  std::size_t x_faces = 0U;
  std::size_t y_faces = 0U;
  std::size_t z_faces = 0U;
  std::size_t sum = 0U;
  return detail::checked_multiply(static_cast<std::size_t>(cx + bx),
                                  static_cast<std::size_t>(cy), x_faces) &&
         detail::checked_multiply(x_faces, static_cast<std::size_t>(cz),
                                  x_faces) &&
         detail::checked_multiply(static_cast<std::size_t>(cx),
                                  static_cast<std::size_t>(cy + by), y_faces) &&
         detail::checked_multiply(y_faces, static_cast<std::size_t>(cz),
                                  y_faces) &&
         detail::checked_multiply(static_cast<std::size_t>(cx),
                                  static_cast<std::size_t>(cy), z_faces) &&
         detail::checked_multiply(z_faces,
                                  static_cast<std::size_t>(cz + bz), z_faces) &&
         detail::checked_add(x_faces, y_faces, sum) &&
         detail::checked_add(sum, z_faces, sum) &&
         (out = static_cast<std::uint64_t>(sum), true);
}

template <bool Uniform>
Status reconstruct_all_axes(const CartesianKernelPlan& plan,
                            ConstFieldView density, ConstFieldView velocity,
                            KernelBox box, FaceFluxView flux) noexcept {
  Status status = reconstruct_axis<Uniform, 0U>(plan, density, velocity, box,
                                                flux.x);
  if (status) {
    status = reconstruct_axis<Uniform, 1U>(plan, density, velocity, box,
                                           flux.y);
  }
  if (status) {
    status = reconstruct_axis<Uniform, 2U>(plan, density, velocity, box,
                                           flux.z);
  }
  return status;
}

template <bool Uniform>
Status divergence_kernel(const CartesianKernelPlan& plan,
                         ConstFaceFluxView flux, KernelBox box,
                         FieldView output,
                         std::uint8_t output_component) noexcept {
  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double raw = flux.x.unchecked({x + 1, y, z}) -
                               flux.x.unchecked({x, y, z}) +
                           flux.y.unchecked({x, y + 1, z}) -
                               flux.y.unchecked({x, y, z}) +
                           flux.z.unchecked({x, y, z + 1}) -
                               flux.z.unchecked({x, y, z});
        const double value =
            raw * detail::metric_inverse_volume<Uniform>(plan, cell);
        if (!std::isfinite(value)) {
          return {StatusCode::numerical_failure, kFaceNumerical};
        }
        output.unchecked(cell, output_component) = value;
      }
    }
  }
  return {};
}

}  // namespace

ConstFaceFieldView as_const(FaceFieldView view) noexcept {
  return {view.base,
          view.extents,
          view.stride_y,
          view.stride_z,
          view.axis,
          view.storage_identity,
          view.revision_domain};
}

ConstFaceFluxView as_const(FaceFluxView view) noexcept {
  return {as_const(view.x), as_const(view.y), as_const(view.z), view.revision,
          view.certificate};
}

bool FaceFluxCertificate::matches(ConstFaceFluxView flux) const noexcept {
  return valid() && flux.certificate == *this && flux.revision == revision_ &&
         flux.x.storage_identity == storage_ &&
         flux.y.storage_identity == storage_ &&
         flux.z.storage_identity == storage_ &&
         flux.x.revision_domain == revision_domain_ &&
         flux.y.revision_domain == revision_domain_ &&
         flux.z.revision_domain == revision_domain_ &&
         flux.x.base == x_base_ && flux.y.base == y_base_ &&
         flux.z.base == z_base_ && flux.x.axis == CartesianAxis::x &&
         flux.y.axis == CartesianAxis::y && flux.z.axis == CartesianAxis::z &&
         flux.x.stride_y == x_stride_y_ && flux.x.stride_z == x_stride_z_ &&
         flux.y.stride_y == y_stride_y_ && flux.y.stride_z == y_stride_z_ &&
         flux.z.stride_y == z_stride_y_ && flux.z.stride_z == z_stride_z_ &&
         flux.x.extents.x == cells_.x + 1 &&
         flux.x.extents.y == cells_.y && flux.x.extents.z == cells_.z &&
         flux.y.extents.x == cells_.x &&
         flux.y.extents.y == cells_.y + 1 &&
         flux.y.extents.z == cells_.z && flux.z.extents.x == cells_.x &&
         flux.z.extents.y == cells_.y &&
         flux.z.extents.z == cells_.z + 1;
}

namespace detail {

struct PendingFaceFluxAccess {
  static FaceFluxView raw(PendingFaceFluxView& pending) noexcept {
    return {pending.x_, pending.y_, pending.z_, pending.revision_, {}};
  }
  static bool lease_valid(const PendingFaceFluxView& pending) noexcept {
    return pending.storage_ != nullptr &&
           pending.storage_->pending_writer_identity_ ==
               pending.writer_identity_ &&
           pending.storage_->pending_attempt_identity_ ==
               pending.attempt_identity_;
  }
};

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
Status overwrite_pending_face_flux_for_test(PendingFaceFluxView& pending,
                                            double value) noexcept {
  if (!pending.valid() || !std::isfinite(value) ||
      !PendingFaceFluxAccess::lease_valid(pending)) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  const FaceFluxView raw = PendingFaceFluxAccess::raw(pending);
  const FaceFieldView faces[3]{raw.x, raw.y, raw.z};
  for (const FaceFieldView face : faces)
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          face.unchecked({x, y, z}) = value;
  return {};
}
#endif

}  // namespace detail

PendingFaceFluxView::~PendingFaceFluxView() noexcept {
  if (writer_ != nullptr) {
    writer_->abandon_pending_view(*this);
  }
}

FaceFluxStorage::~FaceFluxStorage() noexcept { release(); }

void FaceFluxStorage::swap(FaceFluxStorage& other) noexcept {
  using std::swap;
  swap(data_, other.data_);
  swap(cells_, other.cells_);
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    swap(extents_[axis], other.extents_[axis]);
    swap(stride_y_[axis], other.stride_y_[axis]);
    swap(stride_z_[axis], other.stride_z_[axis]);
    swap(offsets_[axis], other.offsets_[axis]);
  }
  swap(replica_stride_, other.replica_stride_);
  swap(replicas_, other.replicas_);
  swap(identity_, other.identity_);
  swap(revision_domain_, other.revision_domain_);
  swap(counters_, other.counters_);
  swap(authority_identity_, other.authority_identity_);
  swap(pending_writer_identity_, other.pending_writer_identity_);
  swap(pending_attempt_identity_, other.pending_attempt_identity_);
  swap(pending_writer_, other.pending_writer_);
  swap(final_storage_, other.final_storage_);
}

void FaceFluxStorage::release() noexcept {
  if (pending_writer_ != nullptr) {
    pending_writer_->detach_from_storage(*this);
  }
  ::operator delete(data_, std::align_val_t{detail::kCacheLineBytes});
  data_ = nullptr;
  cells_ = {};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    extents_[axis] = {};
    stride_y_[axis] = 0U;
    stride_z_[axis] = 0U;
    offsets_[axis] = 0U;
  }
  replica_stride_ = 0U;
  replicas_ = 0U;
  identity_ = 0U;
  revision_domain_ = 0U;
  counters_ = {};
  authority_identity_ = 0U;
  pending_writer_identity_ = 0U;
  pending_attempt_identity_ = 0U;
  pending_writer_ = nullptr;
  final_storage_ = false;
}

Status FaceFluxStorage::allocate_impl(Int3 cells, std::size_t replicas,
                                      bool final_storage,
                                      FaceFluxStorage& out) {
  if (out.data_ != nullptr || out.authority_identity_ != 0U ||
      !detail::valid_cells(cells) || replicas == 0U ||
      cells.x == std::numeric_limits<std::int32_t>::max() ||
      cells.y == std::numeric_limits<std::int32_t>::max() ||
      cells.z == std::numeric_limits<std::int32_t>::max()) {
    return {StatusCode::invalid_plan, kFaceExtent};
  }
  const Int3 extents[3]{{cells.x + 1, cells.y, cells.z},
                        {cells.x, cells.y + 1, cells.z},
                        {cells.x, cells.y, cells.z + 1}};
  if (!detail::valid_cells(extents[0]) || !detail::valid_cells(extents[1]) ||
      !detail::valid_cells(extents[2])) {
    return {StatusCode::invalid_plan, kFaceExtent};
  }
  try {
    FaceFluxStorage candidate;
    candidate.cells_ = cells;
    candidate.final_storage_ = final_storage;
    std::size_t cursor = 0U;
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      candidate.extents_[axis] = extents[axis];
      if (!checked_align(static_cast<std::size_t>(extents[axis].x),
                         candidate.stride_y_[axis]) ||
          !checked_multiply(candidate.stride_y_[axis],
                            static_cast<std::size_t>(extents[axis].y),
                            candidate.stride_z_[axis]) ||
          !checked_align(cursor, candidate.offsets_[axis])) {
        return {StatusCode::invalid_plan, kFaceExtent};
      }
      std::size_t block = 0U;
      if (!checked_multiply(candidate.stride_z_[axis],
                            static_cast<std::size_t>(extents[axis].z), block) ||
          !checked_add(candidate.offsets_[axis], block, cursor)) {
        return {StatusCode::invalid_plan, kFaceExtent};
      }
    }
    if (!checked_align(cursor, candidate.replica_stride_)) {
      return {StatusCode::invalid_plan, kFaceExtent};
    }
    std::size_t doubles = 0U;
    std::size_t bytes = 0U;
    if (!checked_multiply(candidate.replica_stride_, replicas, doubles) ||
        !checked_multiply(doubles, sizeof(double), bytes) || bytes == 0U) {
      return {StatusCode::invalid_plan, kFaceExtent};
    }
    candidate.data_ = static_cast<double*>(::operator new(
        bytes, std::align_val_t{detail::kCacheLineBytes}));
    std::uninitialized_fill_n(candidate.data_, doubles, 0.0);
    candidate.replicas_ = replicas;
    candidate.identity_ = detail::issue_identity();
    candidate.revision_domain_ = detail::issue_identity();
    if (candidate.identity_ == 0U || candidate.revision_domain_ == 0U) {
      return {StatusCode::invalid_plan, kFaceStorage};
    }
    if (replicas > std::numeric_limits<std::uint64_t>::max() / 3U) {
      return {StatusCode::invalid_plan, kFaceExtent};
    }
    candidate.counters_ = {
        1U, static_cast<std::uint64_t>(bytes),
        static_cast<std::uint64_t>(replicas),
        static_cast<std::uint64_t>(replicas) * 3U};
    out.swap(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kFaceStorage};
  }
}

Status FaceFluxStorage::allocate_workspace(Int3 cells, std::size_t replicas,
                                           FaceFluxStorage& out) {
  return allocate_impl(cells, replicas, false, out);
}

Status FaceFluxStorage::allocate_final(Int3 cells, FaceFluxStorage& out) {
  return allocate_impl(cells, kFinalFaceFluxReplicas, true, out);
}

Status FaceFluxStorage::view_impl(std::size_t replica,
                                  RevisionToken revision,
                                  FaceFluxView& out) noexcept {
  if (data_ == nullptr || replica >= replicas_ || revision == 0U ||
      identity_ == 0U || revision_domain_ == 0U) {
    return {StatusCode::invalid_plan, kFaceStorage};
  }
  double* const replica_base = data_ + replica * replica_stride_;
  FaceFluxView candidate;
  candidate.x = make_axis_view(replica_base + offsets_[0U], extents_[0U],
                               stride_y_[0U], stride_z_[0U], CartesianAxis::x,
                               identity_, revision_domain_);
  candidate.y = make_axis_view(replica_base + offsets_[1U], extents_[1U],
                               stride_y_[1U], stride_z_[1U], CartesianAxis::y,
                               identity_, revision_domain_);
  candidate.z = make_axis_view(replica_base + offsets_[2U], extents_[2U],
                               stride_y_[2U], stride_z_[2U], CartesianAxis::z,
                               identity_, revision_domain_);
  candidate.revision = revision;
  candidate.certificate = {};
  out = candidate;
  return {};
}

Status FaceFluxStorage::pending_view_impl(
    std::size_t replica, RevisionToken revision, std::uint64_t writer_identity,
    std::uint64_t attempt_identity, PendingFaceFluxView& out) noexcept {
  FaceFluxView raw;
  const Status viewed = view_impl(replica, revision, raw);
  if (!viewed) {
    return viewed;
  }
  if (writer_identity == 0U || attempt_identity == 0U || !final_storage_) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  if (out.valid()) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  out.x_ = raw.x;
  out.y_ = raw.y;
  out.z_ = raw.z;
  out.revision_ = revision;
  out.writer_identity_ = writer_identity;
  out.attempt_identity_ = attempt_identity;
  out.storage_ = this;
  return {};
}

Status FaceFluxStorage::view_impl(std::size_t replica,
                                  RevisionToken revision,
                                  ConstFaceFluxView& out) const noexcept {
  if (data_ == nullptr || replica >= replicas_ || revision == 0U ||
      identity_ == 0U || revision_domain_ == 0U) {
    return {StatusCode::invalid_plan, kFaceStorage};
  }
  const double* const replica_base = data_ + replica * replica_stride_;
  ConstFaceFluxView candidate;
  candidate.x = make_axis_view(replica_base + offsets_[0U], extents_[0U],
                               stride_y_[0U], stride_z_[0U], CartesianAxis::x,
                               identity_, revision_domain_);
  candidate.y = make_axis_view(replica_base + offsets_[1U], extents_[1U],
                               stride_y_[1U], stride_z_[1U], CartesianAxis::y,
                               identity_, revision_domain_);
  candidate.z = make_axis_view(replica_base + offsets_[2U], extents_[2U],
                               stride_y_[2U], stride_z_[2U], CartesianAxis::z,
                               identity_, revision_domain_);
  candidate.revision = revision;
  candidate.certificate = {};
  out = candidate;
  return {};
}

Status FaceFluxStorage::workspace_view(std::size_t replica,
                                       RevisionToken revision,
                                       FaceFluxView& out) noexcept {
  if (final_storage_) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  return view_impl(replica, revision, out);
}

Status FaceFluxStorage::view(std::size_t replica, RevisionToken revision,
                             ConstFaceFluxView& out) const noexcept {
  if (final_storage_) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  return view_impl(replica, revision, out);
}

Status FinalFaceFluxAuthority::claim(StageId stage, RevisionSlotId cache_slot,
                                     AttemptTransaction& transaction,
                                     FinalFaceFluxWriter& out) noexcept {
  if (claimed_ || out.issued_ || stage == 0U) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  const PlanFingerprint identity = detail::issue_identity();
  if (identity == 0U) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  const Status mandatory =
      transaction.claim_final_face_flux_writer(cache_slot, identity, out);
  if (!mandatory) {
    return mandatory;
  }
  out.stage_ = stage;
  out.cache_slot_ = cache_slot;
  out.authority_fingerprint_ = identity;
  out.authority_transaction_ = &transaction;
  out.issued_ = true;
  fingerprint_ = identity;
  claimed_ = true;
  return {};
}

Status FinalFaceFluxWriter::begin_pending(AttemptTransaction& transaction,
                                          FaceFluxStorage& storage,
                                          PendingFaceFluxView& out) noexcept {
  if (!issued_ || pending_ || authority_transaction_ != &transaction ||
      !transaction.active() ||
      transaction.attempt_identity() == 0U ||
      storage.replicas_ != kFinalFaceFluxReplicas ||
      !storage.final_storage_ ||
      storage_ != nullptr || transaction_ != nullptr ||
      (storage.authority_identity_ != 0U &&
       storage.authority_identity_ != authority_fingerprint_) ||
      (bound_storage_identity_ != 0U &&
       (bound_storage_identity_ != storage.identity_ ||
        bound_revision_domain_ != storage.revision_domain_ ||
       !detail::same_cells(bound_cells_, storage.cells_))) ||
      pending_replica_ >= kFinalFaceFluxReplicas ||
      pending_replica_ == active_replica_ ||
      (previous_revision_ != 0U &&
       pending_replica_ == previous_replica_)) {
    return {StatusCode::invalid_plan, kFaceTransaction};
  }
  const RevisionToken candidate =
      active_revision_ == 0U ? RevisionToken{1U} : active_revision_ + 1U;
  if (candidate == 0U || candidate <= active_revision_) {
    return {StatusCode::invalid_plan, kFaceTransaction};
  }
  const Status viewed = storage.pending_view_impl(
      pending_replica_, candidate, authority_fingerprint_,
      transaction.attempt_identity(), out);
  if (!viewed) {
    return viewed;
  }
  storage.pending_writer_identity_ = authority_fingerprint_;
  storage.pending_attempt_identity_ = transaction.attempt_identity();
  storage.pending_writer_ = this;
  out.writer_ = this;
  pending_revision_ = candidate;
  pending_attempt_identity_ = transaction.attempt_identity();
  transaction_ = &transaction;
  storage_ = &storage;
  pending_view_ = &out;
  if (bound_storage_identity_ == 0U) {
    bound_storage_identity_ = storage.identity_;
    bound_revision_domain_ = storage.revision_domain_;
    bound_cells_ = storage.cells_;
  }
  storage.authority_identity_ = authority_fingerprint_;
  pending_ = true;
  pending_published_ = false;
  return {};
}

Status FinalFaceFluxWriter::initialize_committed(
    FaceFluxStorage& storage, ConstFaceFluxView source) noexcept {
  if (!issued_ || active_revision_ != 0U || previous_revision_ != 0U ||
      pending_ ||
      authority_transaction_ == nullptr || authority_transaction_->active() ||
      storage.replicas_ != kFinalFaceFluxReplicas ||
      !storage.final_storage_ ||
      storage.authority_identity_ != 0U || storage.pending_writer_ != nullptr ||
      storage_ != nullptr || transaction_ != nullptr ||
      !detail::valid_flux_view(source, storage.cells_, source.revision)) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  const ConstFaceFieldView inputs[3]{source.x, source.y, source.z};
  for (const ConstFaceFieldView& input : inputs) {
    for (std::int32_t z = 0; z < input.extents.z; ++z)
      for (std::int32_t y = 0; y < input.extents.y; ++y)
        for (std::int32_t x = 0; x < input.extents.x; ++x)
          if (!std::isfinite(input.unchecked({x, y, z})))
            return {StatusCode::numerical_failure, kFaceAuthority};
  }
  const RevisionToken initial_revision = 1U;
  constexpr std::size_t initial_replica = 0U;
  FaceFluxView destination;
  Status status =
      storage.view_impl(initial_replica, initial_revision, destination);
  if (!status) return status;
  const FaceFieldView outputs[3]{destination.x, destination.y, destination.z};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    const ConstFaceFieldView input = inputs[axis];
    const FaceFieldView output = outputs[axis];
    for (std::int32_t z = 0; z < input.extents.z; ++z)
      for (std::int32_t y = 0; y < input.extents.y; ++y)
        for (std::int32_t x = 0; x < input.extents.x; ++x)
          output.unchecked({x, y, z}) = input.unchecked({x, y, z});
  }
  active_replica_ = initial_replica;
  previous_replica_ = initial_replica;
  pending_replica_ = 1U;
  active_revision_ = initial_revision;
  // A fresh state has two mathematically identical history levels.  Publish
  // both authorities at once so BDF history consumers never observe an
  // address without a corresponding revision/certificate lineage.
  previous_revision_ = initial_revision;
  bound_storage_identity_ = storage.identity_;
  bound_revision_domain_ = storage.revision_domain_;
  bound_cells_ = storage.cells_;
  storage.authority_identity_ = authority_fingerprint_;
  return {};
}

Status FinalFaceFluxWriter::initialize_restored(
    FaceFluxStorage& storage, ConstFaceFluxView source) noexcept {
  if (!issued_ || active_revision_ != 0U || previous_revision_ != 0U ||
      pending_ || authority_transaction_ == nullptr ||
      authority_transaction_->active() ||
      storage.replicas_ != kFinalFaceFluxReplicas ||
      !storage.final_storage_ || storage.authority_identity_ != 0U ||
      storage.pending_writer_ != nullptr || storage_ != nullptr ||
      transaction_ != nullptr || bound_storage_identity_ != 0U ||
      bound_revision_domain_ != 0U ||
      !detail::valid_flux_view(source, storage.cells_, source.revision)) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  const ConstFaceFieldView inputs[3]{source.x, source.y, source.z};
  for (const ConstFaceFieldView input : inputs) {
    for (std::int32_t z = 0; z < input.extents.z; ++z)
      for (std::int32_t y = 0; y < input.extents.y; ++y)
        for (std::int32_t x = 0; x < input.extents.x; ++x)
          if (!std::isfinite(input.unchecked({x, y, z})))
            return {StatusCode::numerical_failure, kFaceAuthority};
  }

  // Revision one is reserved for a genuine fresh-start baseline.  A restart
  // enters at revision two, matching the lineage previously produced by
  // initialize_committed()+restore_committed() without materializing that
  // incompatible intermediate state.
  constexpr RevisionToken restored_revision = 2U;
  constexpr std::size_t restored_replica = 0U;
  FaceFluxView destination;
  Status status =
      storage.view_impl(restored_replica, restored_revision, destination);
  if (!status) return status;
  const FaceFieldView outputs[3]{destination.x, destination.y, destination.z};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    const ConstFaceFieldView input = inputs[axis];
    const FaceFieldView output = outputs[axis];
    for (std::int32_t z = 0; z < input.extents.z; ++z)
      for (std::int32_t y = 0; y < input.extents.y; ++y)
        for (std::int32_t x = 0; x < input.extents.x; ++x)
          output.unchecked({x, y, z}) = input.unchecked({x, y, z});
  }
  active_replica_ = restored_replica;
  previous_replica_ = restored_replica;
  pending_replica_ = 1U;
  active_revision_ = restored_revision;
  previous_revision_ = restored_revision;
  bound_storage_identity_ = storage.identity_;
  bound_revision_domain_ = storage.revision_domain_;
  bound_cells_ = storage.cells_;
  storage.authority_identity_ = authority_fingerprint_;
  return {};
}

Status FinalFaceFluxWriter::initialize_restored_history(
    FaceFluxStorage& storage, ConstFaceFluxView accepted,
    ConstFaceFluxView previous) noexcept {
  if (!issued_ || active_revision_ != 0U || previous_revision_ != 0U ||
      pending_ || authority_transaction_ == nullptr ||
      authority_transaction_->active() ||
      storage.replicas_ != kFinalFaceFluxReplicas ||
      !storage.final_storage_ || storage.authority_identity_ != 0U ||
      storage.pending_writer_ != nullptr || storage_ != nullptr ||
      transaction_ != nullptr || bound_storage_identity_ != 0U ||
      bound_revision_domain_ != 0U ||
      !detail::valid_flux_view(accepted, storage.cells_, accepted.revision) ||
      !detail::valid_flux_view(previous, storage.cells_, previous.revision) ||
      previous.revision >= accepted.revision ||
      accepted.revision == std::numeric_limits<RevisionToken>::max()) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  const std::array<ConstFaceFluxView, 2U> sources{{previous, accepted}};
  for (const ConstFaceFluxView source : sources) {
    const std::array<ConstFaceFieldView, 3U> inputs{
        source.x, source.y, source.z};
    for (const ConstFaceFieldView input : inputs)
      for (std::int32_t z = 0; z < input.extents.z; ++z)
        for (std::int32_t y = 0; y < input.extents.y; ++y)
          for (std::int32_t x = 0; x < input.extents.x; ++x)
            if (!std::isfinite(input.unchecked({x, y, z})))
              return {StatusCode::numerical_failure, kFaceAuthority};
  }

  constexpr std::size_t previous_replica = 0U;
  constexpr std::size_t accepted_replica = 1U;
  const auto copy_to_replica = [&](ConstFaceFluxView source,
                                   std::size_t replica) noexcept {
    FaceFluxView destination;
    Status status = storage.view_impl(replica, source.revision, destination);
    if (!status) return status;
    const std::array<ConstFaceFieldView, 3U> inputs{
        source.x, source.y, source.z};
    const std::array<FaceFieldView, 3U> outputs{
        destination.x, destination.y, destination.z};
    for (std::size_t axis = 0U; axis < inputs.size(); ++axis)
      for (std::int32_t z = 0; z < inputs[axis].extents.z; ++z)
        for (std::int32_t y = 0; y < inputs[axis].extents.y; ++y)
          for (std::int32_t x = 0; x < inputs[axis].extents.x; ++x)
            outputs[axis].unchecked({x, y, z}) =
                inputs[axis].unchecked({x, y, z});
    return Status{};
  };
  Status status = copy_to_replica(previous, previous_replica);
  if (status) status = copy_to_replica(accepted, accepted_replica);
  if (!status) return status;

  active_replica_ = accepted_replica;
  previous_replica_ = previous_replica;
  pending_replica_ = 2U;
  active_revision_ = accepted.revision;
  previous_revision_ = previous.revision;
  bound_storage_identity_ = storage.identity_;
  bound_revision_domain_ = storage.revision_domain_;
  bound_cells_ = storage.cells_;
  storage.authority_identity_ = authority_fingerprint_;
  return {};
}

Status FinalFaceFluxWriter::restore_committed(
    FaceFluxStorage& storage, ConstFaceFluxView source) noexcept {
  if (!issued_ || active_revision_ == 0U || pending_ ||
      authority_transaction_ == nullptr || authority_transaction_->active() ||
      storage.replicas_ != kFinalFaceFluxReplicas ||
      !storage.final_storage_ ||
      storage.authority_identity_ != authority_fingerprint_ ||
      storage.pending_writer_ != nullptr || storage_ != nullptr ||
      transaction_ != nullptr || bound_storage_identity_ != storage.identity_ ||
      bound_revision_domain_ != storage.revision_domain_ ||
      !detail::same_cells(bound_cells_, storage.cells_) ||
      !detail::valid_flux_view(source, storage.cells_, source.revision))
    return {StatusCode::invalid_plan, kFaceAuthority};
  const ConstFaceFieldView inputs[3]{source.x, source.y, source.z};
  for (const ConstFaceFieldView input : inputs)
    for (std::int32_t z = 0; z < input.extents.z; ++z)
      for (std::int32_t y = 0; y < input.extents.y; ++y)
        for (std::int32_t x = 0; x < input.extents.x; ++x)
          if (!std::isfinite(input.unchecked({x, y, z})))
            return {StatusCode::numerical_failure, kFaceAuthority};
  const RevisionToken restored_revision = active_revision_ + 1U;
  if (restored_revision == 0U)
    return {StatusCode::invalid_plan, kFaceAuthority};
  FaceFluxView destination;
  Status status =
      storage.view_impl(active_replica_, restored_revision, destination);
  if (!status) return status;
  const FaceFieldView outputs[3]{destination.x, destination.y, destination.z};
  for (std::size_t axis = 0U; axis < 3U; ++axis)
    for (std::int32_t z = 0; z < inputs[axis].extents.z; ++z)
      for (std::int32_t y = 0; y < inputs[axis].extents.y; ++y)
        for (std::int32_t x = 0; x < inputs[axis].extents.x; ++x)
          outputs[axis].unchecked({x, y, z}) =
              inputs[axis].unchecked({x, y, z});
  previous_replica_ = active_replica_;
  previous_revision_ = restored_revision;
  active_revision_ = restored_revision;
  return {};
}

Status FinalFaceFluxWriter::publish_pending(
    Span<const RevisionDependency> dependencies,
    PendingFaceFluxView& pending) noexcept {
  const Status preflight = preflight_publish_pending(dependencies, pending);
  if (!preflight) {
    // Dependency/cache failures historically poison a direct publication.
    // Re-enter the transaction publisher only for those failures so it can
    // latch the exact status; writer/view association failures retain their
    // existing local rejection semantics.
    if (preflight.detail != kFaceTransaction && transaction_ != nullptr) {
      return transaction_->publish_final_face_flux_cache(
          cache_slot_, dependencies, PendingCacheStamp{pending_revision_},
          authority_fingerprint_);
    }
    return preflight;
  }
  const Status published = transaction_->publish_final_face_flux_cache(
      cache_slot_, dependencies, PendingCacheStamp{pending_revision_},
      authority_fingerprint_);
  if (published) {
    pending_published_ = true;
    pending.x_ = {};
    pending.y_ = {};
    pending.z_ = {};
    pending.revision_ = 0U;
    pending.writer_identity_ = 0U;
    pending.attempt_identity_ = 0U;
    pending.storage_ = nullptr;
    pending.writer_ = nullptr;
    pending_view_ = nullptr;
  }
  return published;
}

Status FinalFaceFluxWriter::preflight_publish_pending(
    Span<const RevisionDependency> dependencies,
    const PendingFaceFluxView& pending) const noexcept {
  const FaceFluxView pending_flux{pending.x_, pending.y_, pending.z_,
                                  pending.revision_, {}};
  if (!issued_ || !pending_ || pending_published_ || transaction_ == nullptr ||
      storage_ == nullptr || authority_transaction_ != transaction_ ||
      pending_view_ != &pending || pending.writer_ != this ||
      transaction_->final_face_flux_writer_ != this ||
      transaction_->final_face_flux_writer_identity_ !=
          authority_fingerprint_ ||
      !transaction_->active() ||
      transaction_->attempt_identity() != pending_attempt_identity_ ||
      pending.writer_identity_ != authority_fingerprint_ ||
      pending.attempt_identity_ != pending_attempt_identity_ ||
      pending.storage_ != storage_ ||
      storage_->pending_writer_ != this ||
      storage_->pending_writer_identity_ != authority_fingerprint_ ||
      storage_->pending_attempt_identity_ != pending_attempt_identity_ ||
      storage_->authority_identity_ != authority_fingerprint_ ||
      bound_storage_identity_ != storage_->identity_ ||
      bound_revision_domain_ != storage_->revision_domain_ ||
      !detail::same_cells(bound_cells_, storage_->cells_) ||
      pending.revision_ != pending_revision_ || pending_revision_ == 0U ||
      !detail::valid_flux_view(pending_flux, storage_->cells_) ||
      pending.x_.storage_identity != storage_->identity_ ||
      pending.x_.revision_domain != storage_->revision_domain_) {
    return {StatusCode::invalid_plan, kFaceTransaction};
  }
  return transaction_->preflight_final_face_flux_cache(
      cache_slot_, dependencies, PendingCacheStamp{pending_revision_},
      authority_fingerprint_);
}

bool FinalFaceFluxWriter::ready_for_collective(
    const AttemptTransaction& transaction) const noexcept {
  return pending_ && pending_published_ && transaction_ == &transaction &&
         storage_ != nullptr && transaction.active() &&
         transaction.attempt_identity() == pending_attempt_identity_ &&
         bound_storage_identity_ == storage_->identity_ &&
         bound_revision_domain_ == storage_->revision_domain_ &&
         detail::same_cells(bound_cells_, storage_->cells_) &&
         transaction.pending_caches_[cache_slot_] == pending_revision_;
}

void FinalFaceFluxWriter::complete_from_transaction(
    const AttemptTransaction& transaction, bool committed) noexcept {
  if (transaction_ != &transaction) {
    return;
  }
  if (pending_ && committed) {
    previous_replica_ = active_replica_;
    previous_revision_ = active_revision_;
    active_replica_ = pending_replica_;
    active_revision_ = pending_revision_;
    for (std::size_t replica = 0U; replica < kFinalFaceFluxReplicas;
         ++replica) {
      if (replica != active_replica_ && replica != previous_replica_) {
        pending_replica_ = replica;
        break;
      }
    }
  } else if (pending_) {
    // The rejected pending replica is already disjoint from both committed
    // time levels. Reuse it without touching any committed byte or handle.
  }
  clear_storage_lease();
  invalidate_pending_view();
  pending_revision_ = 0U;
  pending_attempt_identity_ = 0U;
  transaction_ = nullptr;
  storage_ = nullptr;
  pending_ = false;
  pending_published_ = false;
}

FinalFaceFluxWriter::~FinalFaceFluxWriter() noexcept {
  AttemptTransaction* const transaction = authority_transaction_;
  if (transaction != nullptr) {
    transaction->detach_final_face_flux_writer(*this, true);
  } else {
    clear_storage_lease();
    invalidate_pending_view();
    transaction_ = nullptr;
    storage_ = nullptr;
    pending_ = false;
    pending_published_ = false;
  }
}

void FinalFaceFluxWriter::clear_storage_lease() noexcept {
  if (storage_ != nullptr && storage_->pending_writer_ == this) {
    storage_->pending_writer_identity_ = 0U;
    storage_->pending_attempt_identity_ = 0U;
    storage_->pending_writer_ = nullptr;
  }
}

void FinalFaceFluxWriter::invalidate_pending_view() noexcept {
  if (pending_view_ != nullptr) {
    pending_view_->x_ = {};
    pending_view_->y_ = {};
    pending_view_->z_ = {};
    pending_view_->revision_ = 0U;
    pending_view_->writer_identity_ = 0U;
    pending_view_->attempt_identity_ = 0U;
    pending_view_->storage_ = nullptr;
    pending_view_->writer_ = nullptr;
    pending_view_ = nullptr;
  }
}

void FinalFaceFluxWriter::detach_from_transaction(
    AttemptTransaction& transaction) noexcept {
  if (transaction_ == &transaction) {
    clear_storage_lease();
    invalidate_pending_view();
    pending_revision_ = 0U;
    pending_attempt_identity_ = 0U;
    transaction_ = nullptr;
    storage_ = nullptr;
    pending_ = false;
    pending_published_ = false;
  }
  if (authority_transaction_ == &transaction) {
    authority_transaction_ = nullptr;
  }
}

void FinalFaceFluxWriter::detach_from_storage(
    FaceFluxStorage& storage) noexcept {
  if (storage_ != &storage) {
    return;
  }
  if (transaction_ != nullptr) {
    transaction_->invalidate_final_face_flux_writer(*this);
  }
  clear_storage_lease();
  invalidate_pending_view();
  pending_revision_ = 0U;
  pending_attempt_identity_ = 0U;
  storage_ = nullptr;
  pending_ = false;
  pending_published_ = false;
}

void FinalFaceFluxWriter::abandon_pending_view(
    PendingFaceFluxView& pending) noexcept {
  if (pending_view_ != &pending) {
    return;
  }
  if (transaction_ != nullptr) {
    transaction_->invalidate_final_face_flux_writer(*this);
  }
  clear_storage_lease();
  pending_view_ = nullptr;
  pending_revision_ = 0U;
  pending_attempt_identity_ = 0U;
  transaction_ = nullptr;
  storage_ = nullptr;
  pending_ = false;
  pending_published_ = false;
  pending.x_ = {};
  pending.y_ = {};
  pending.z_ = {};
  pending.revision_ = 0U;
  pending.writer_identity_ = 0U;
  pending.attempt_identity_ = 0U;
  pending.storage_ = nullptr;
  pending.writer_ = nullptr;
}

Status FinalFaceFluxWriter::committed(const FaceFluxStorage& storage,
                                      ConstFaceFluxView& out) const noexcept {
  if (!issued_ || active_revision_ == 0U ||
      storage.replicas_ != kFinalFaceFluxReplicas ||
      !storage.final_storage_ ||
      bound_storage_identity_ != storage.identity_ ||
      bound_revision_domain_ != storage.revision_domain_ ||
      !detail::same_cells(bound_cells_, storage.cells_)) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  const Status viewed = storage.view_impl(active_replica_, active_revision_,
                                          out);
  if (viewed) {
    out.certificate.revision_ = active_revision_;
    out.certificate.authority_ = authority_fingerprint_;
    out.certificate.storage_ = storage.identity_;
    out.certificate.revision_domain_ = storage.revision_domain_;
    out.certificate.x_base_ = out.x.base;
    out.certificate.y_base_ = out.y.base;
    out.certificate.z_base_ = out.z.base;
    out.certificate.x_stride_y_ = out.x.stride_y;
    out.certificate.x_stride_z_ = out.x.stride_z;
    out.certificate.y_stride_y_ = out.y.stride_y;
    out.certificate.y_stride_z_ = out.y.stride_z;
    out.certificate.z_stride_y_ = out.z.stride_y;
    out.certificate.z_stride_z_ = out.z.stride_z;
    out.certificate.cells_ = storage.cells_;
  }
  return viewed;
}

Status FinalFaceFluxWriter::committed_previous(
    const FaceFluxStorage& storage, ConstFaceFluxView& out) const noexcept {
  if (!issued_ || previous_revision_ == 0U ||
      storage.replicas_ != kFinalFaceFluxReplicas ||
      !storage.final_storage_ ||
      bound_storage_identity_ != storage.identity_ ||
      bound_revision_domain_ != storage.revision_domain_ ||
      !detail::same_cells(bound_cells_, storage.cells_)) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  const Status viewed =
      storage.view_impl(previous_replica_, previous_revision_, out);
  if (viewed) {
    out.certificate.revision_ = previous_revision_;
    out.certificate.authority_ = authority_fingerprint_;
    out.certificate.storage_ = storage.identity_;
    out.certificate.revision_domain_ = storage.revision_domain_;
    out.certificate.x_base_ = out.x.base;
    out.certificate.y_base_ = out.y.base;
    out.certificate.z_base_ = out.z.base;
    out.certificate.x_stride_y_ = out.x.stride_y;
    out.certificate.x_stride_z_ = out.x.stride_z;
    out.certificate.y_stride_y_ = out.y.stride_y;
    out.certificate.y_stride_z_ = out.y.stride_z;
    out.certificate.z_stride_y_ = out.z.stride_y;
    out.certificate.z_stride_z_ = out.z.stride_z;
    out.certificate.cells_ = storage.cells_;
  }
  return viewed;
}

Status FaceFluxConsumer::bind(FaceFluxCertificate certificate) noexcept {
  if (!certificate.valid() ||
      (required_certificate_.valid() && consumed_revision_ == 0U)) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  required_certificate_ = certificate;
  consumed_revision_ = 0U;
  return {};
}

Status FaceFluxConsumer::consume(ConstFaceFluxView flux) noexcept {
  if (!required_certificate_.matches(flux) || consumed_revision_ != 0U) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  consumed_revision_ = flux.revision;
  return {};
}

Status reconstruct_mass_flux(const CartesianKernelPlan& plan,
                             const KernelInvocation& invocation,
                             FaceFluxView& flux) noexcept {
  if (plan.fingerprint() == 0U || invocation.reads.data == nullptr ||
      invocation.reads.size != 2U || invocation.writes.size != 0U ||
      !detail::valid_kernel_box(invocation.box, plan.cells()) ||
      !detail::valid_cell_view(invocation.reads.data[0U], plan.cells(), 0U,
                               1U, 1U) ||
      !detail::valid_cell_view(invocation.reads.data[1U], plan.cells(), 0U,
                               3U, 1U) ||
      !detail::valid_flux_view(flux, plan.cells()) ||
      flux.certificate.valid() ||
      invocation.component_count != 1U ||
      invocation.required_face_flux_revision != 0U) {
    return {StatusCode::invalid_plan, kFaceKernel};
  }
  const ConstFieldView density = invocation.reads.data[0U];
  const ConstFieldView velocity = invocation.reads.data[1U];
  std::uint64_t faces = 0U;
  if (!checked_face_count(invocation.box, faces) ||
      faces > std::numeric_limits<std::uint64_t>::max() / 4U) {
    return {StatusCode::invalid_plan, kFaceKernel};
  }
  KernelCounters prepared{};
  if (invocation.counters != nullptr) {
    prepared = *invocation.counters;
    const Status counted = detail::add_kernel_counter_totals(
        &prepared, detail::box_cell_count(invocation.box), faces, 4U * faces,
        faces);
    if (!counted) {
      return counted;
    }
  }
  const Status status =
      plan.geometry_kind() == GeometryKind::uniform
          ? reconstruct_all_axes<true>(plan, density, velocity, invocation.box,
                                       flux)
          : reconstruct_all_axes<false>(plan, density, velocity,
                                        invocation.box, flux);
  if (!status) {
    return status;
  }
  if (invocation.counters != nullptr) {
    *invocation.counters = prepared;
  }
  return {};
}

Status reconstruct_mass_flux(const CartesianKernelPlan& plan,
                             const KernelInvocation& invocation,
                             PendingFaceFluxView& pending) noexcept {
  if (!pending.valid()) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  if (!detail::PendingFaceFluxAccess::lease_valid(pending)) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  FaceFluxView raw = detail::PendingFaceFluxAccess::raw(pending);
  return reconstruct_mass_flux(plan, invocation, raw);
}

Status cartesian_provisional_face_divergence(
    const CartesianKernelPlan& plan, ConstFaceFluxView flux,
    const KernelInvocation& invocation) noexcept {
  if (plan.fingerprint() == 0U || invocation.reads.size != 0U ||
      invocation.writes.data == nullptr || invocation.writes.size != 1U ||
      invocation.component_count != 1U ||
      invocation.required_face_flux_revision == 0U ||
      !detail::valid_kernel_box(invocation.box, plan.cells()) ||
      !detail::valid_flux_view(flux, plan.cells(),
                               invocation.required_face_flux_revision) ||
      !detail::valid_cell_view(invocation.writes.data[0U], plan.cells(),
                               invocation.write_component_begin, 1U)) {
    return {StatusCode::invalid_plan, kFaceKernel};
  }
  KernelCounters prepared{};
  if (invocation.counters != nullptr) {
    prepared = *invocation.counters;
    const Status counted = detail::add_kernel_counters(
        &prepared, invocation.box, 6U, 6U, 1U);
    if (!counted) {
      return counted;
    }
  }
  const FieldView output = invocation.writes.data[0U];
  const Status evaluated =
      plan.geometry_kind() == GeometryKind::uniform
          ? divergence_kernel<true>(plan, flux, invocation.box, output,
                                    invocation.write_component_begin)
          : divergence_kernel<false>(plan, flux, invocation.box, output,
                                     invocation.write_component_begin);
  if (!evaluated) {
    return evaluated;
  }
  if (invocation.counters != nullptr) {
    *invocation.counters = prepared;
  }
  return {};
}

Status cartesian_face_divergence(const CartesianKernelPlan& plan,
                                 ConstFaceFluxView flux,
                                 const KernelInvocation& invocation) noexcept {
  if (!flux.certificate.matches(flux)) {
    return {StatusCode::invalid_plan, kFaceAuthority};
  }
  return cartesian_provisional_face_divergence(plan, flux, invocation);
}

}  // namespace hundun::v04
