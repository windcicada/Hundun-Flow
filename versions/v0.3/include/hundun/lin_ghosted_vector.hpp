// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/exec_execution.hpp"
#include "hundun/mesh_topology.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace hundun::linear {

class VectorLayout final {
 public:
  VectorLayout() noexcept = default;
  static VectorLayout from_topology(const mesh::MeshTopology& topology);
  VectorLayout(std::size_t owned_count,
               std::vector<mesh::GlobalCellId> ordered_global_ids);

  std::size_t owned_count() const noexcept;
  std::size_t ghost_count() const noexcept;
  std::size_t local_count() const noexcept;
  const std::vector<mesh::GlobalCellId>& global_ids() const noexcept;

  friend bool operator==(const VectorLayout& left,
                         const VectorLayout& right) noexcept;
  friend bool operator!=(const VectorLayout& left,
                         const VectorLayout& right) noexcept;

 private:
  struct Data;
  explicit VectorLayout(std::shared_ptr<const Data> data) noexcept;
  std::shared_ptr<const Data> data_;
};

class GhostedVector final {
 public:
  GhostedVector(execution::ExecutionContext& context, VectorLayout layout);
  ~GhostedVector() = default;

  GhostedVector(GhostedVector&&) noexcept;
  GhostedVector& operator=(GhostedVector&&) noexcept;
  GhostedVector(const GhostedVector&) = delete;
  GhostedVector& operator=(const GhostedVector&) = delete;

  std::size_t owned_count() const noexcept;
  std::size_t ghost_count() const noexcept;
  std::size_t local_count() const noexcept;
  const VectorLayout& layout() const noexcept;
  execution::AllocationIdentity allocation_identity() const noexcept;
  std::uint64_t epoch() const noexcept;
  execution::BackendIdentity backend_identity() const noexcept;
  execution::ExecutionSpace space() const noexcept;

  execution::VectorView<double> local_view();
  execution::VectorView<const double> local_view() const;
  execution::VectorView<double> owned_view();
  execution::VectorView<const double> owned_view() const;
  execution::VectorView<double> ghost_view();
  execution::VectorView<const double> ghost_view() const;

 private:
  VectorLayout layout_;
  execution::Buffer buffer_;
};

}  // namespace hundun::linear
