// SPDX-License-Identifier: Apache-2.0

#include "hundun/linear/ghosted_vector.hpp"

#include "ghosted_vector_halo_detail.hpp"
#include "hundun/runtime/error.hpp"

#include <algorithm>
#include <memory>
#include <new>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace hundun::linear {

struct VectorLayout::Data final {
  Data(std::size_t owned,
       std::vector<mesh::GlobalCellId> ordered_ids) noexcept
      : owned_count(owned), global_ids(std::move(ordered_ids)) {}

  std::size_t owned_count{};
  std::vector<mesh::GlobalCellId> global_ids;
};

namespace {

const std::vector<mesh::GlobalCellId>& empty_ids() noexcept {
  static const std::vector<mesh::GlobalCellId> value;
  return value;
}

}  // namespace

VectorLayout::VectorLayout(std::shared_ptr<const Data> data) noexcept
    : data_(std::move(data)) {}

VectorLayout VectorLayout::from_topology(
    const mesh::MeshTopology& topology) {
  std::vector<mesh::GlobalCellId> ids;
  try {
    ids.reserve(topology.local_cell_count());
    for (std::size_t local = 0; local < topology.local_cell_count(); ++local) {
      ids.push_back(topology.global_cell_id(local));
    }
  } catch (const runtime::Error&) {
    throw;
  } catch (const std::bad_alloc&) {
    throw runtime::Error("vector layout topology capture allocation failed");
  } catch (const std::length_error&) {
    throw runtime::Error("vector layout topology capture is unsupported");
  }
  return VectorLayout(topology.owned_cell_count(), std::move(ids));
}

VectorLayout::VectorLayout(
    std::size_t owned_count,
    std::vector<mesh::GlobalCellId> ordered_global_ids) {
  if (owned_count > ordered_global_ids.size()) {
    throw runtime::Error(
        "vector layout owned count exceeds its ordered global IDs");
  }
  try {
    std::unordered_set<mesh::GlobalCellId> unique;
    unique.reserve(ordered_global_ids.size());
    for (const mesh::GlobalCellId id : ordered_global_ids) {
      if (!unique.insert(id).second) {
        throw runtime::Error("vector layout contains a duplicate global ID");
      }
    }
    data_ = std::make_shared<Data>(owned_count,
                                   std::move(ordered_global_ids));
  } catch (const runtime::Error&) {
    throw;
  } catch (const std::bad_alloc&) {
    throw runtime::Error("vector layout allocation failed");
  } catch (const std::length_error&) {
    throw runtime::Error("vector layout size is unsupported");
  }
}

std::size_t VectorLayout::owned_count() const noexcept {
  return data_ ? data_->owned_count : 0U;
}

std::size_t VectorLayout::ghost_count() const noexcept {
  return local_count() - owned_count();
}

std::size_t VectorLayout::local_count() const noexcept {
  return data_ ? data_->global_ids.size() : 0U;
}

const std::vector<mesh::GlobalCellId>& VectorLayout::global_ids() const
    noexcept {
  return data_ ? data_->global_ids : empty_ids();
}

bool operator==(const VectorLayout& left, const VectorLayout& right) noexcept {
  if (left.data_ == right.data_) {
    return true;
  }
  return left.owned_count() == right.owned_count() &&
         left.global_ids() == right.global_ids();
}

bool operator!=(const VectorLayout& left, const VectorLayout& right) noexcept {
  return !(left == right);
}

GhostedVector::GhostedVector(execution::ExecutionContext& context,
                             VectorLayout layout)
    : layout_(std::move(layout)),
      buffer_(context, detail::checked_vector_bytes(layout_.local_count())) {}

GhostedVector::GhostedVector(GhostedVector&&) noexcept = default;
GhostedVector& GhostedVector::operator=(GhostedVector&&) noexcept = default;

std::size_t GhostedVector::owned_count() const noexcept {
  return layout_.owned_count();
}

std::size_t GhostedVector::ghost_count() const noexcept {
  return layout_.ghost_count();
}

std::size_t GhostedVector::local_count() const noexcept {
  return layout_.local_count();
}

const VectorLayout& GhostedVector::layout() const noexcept { return layout_; }

execution::AllocationIdentity GhostedVector::allocation_identity() const
    noexcept {
  return buffer_.allocation_identity();
}

std::uint64_t GhostedVector::epoch() const noexcept { return buffer_.epoch(); }

execution::BackendIdentity GhostedVector::backend_identity() const noexcept {
  return buffer_.backend_identity();
}

execution::ExecutionSpace GhostedVector::space() const noexcept {
  return buffer_.space();
}

execution::VectorView<double> GhostedVector::local_view() {
  return buffer_.view(0U, local_count());
}

execution::VectorView<const double> GhostedVector::local_view() const {
  return buffer_.view(0U, local_count());
}

execution::VectorView<double> GhostedVector::owned_view() {
  return buffer_.view(0U, owned_count());
}

execution::VectorView<const double> GhostedVector::owned_view() const {
  return buffer_.view(0U, owned_count());
}

execution::VectorView<double> GhostedVector::ghost_view() {
  return buffer_.view(owned_count(), ghost_count());
}

execution::VectorView<const double> GhostedVector::ghost_view() const {
  return buffer_.view(owned_count(), ghost_count());
}

}  // namespace hundun::linear
