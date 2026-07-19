// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/linear/ghosted_vector.hpp"

#include <cstddef>
#include <memory>

namespace hundun::runtime {
class StructuredDecomposition;
}

namespace hundun::linear {

enum class BufferHaloPath {
  host_direct,
  device_direct,
  device_host_staged
};

class GhostedVectorHalo final {
 public:
  static GhostedVectorHalo create(
      const runtime::StructuredDecomposition& decomposition,
      const mesh::MeshTopology& topology,
      execution::ExecutionContext& context);

  ~GhostedVectorHalo() noexcept;
  GhostedVectorHalo(GhostedVectorHalo&&) noexcept;
  GhostedVectorHalo& operator=(GhostedVectorHalo&&) = delete;
  GhostedVectorHalo(const GhostedVectorHalo&) = delete;
  GhostedVectorHalo& operator=(const GhostedVectorHalo&) = delete;

  BufferHaloPath path() const;
  std::size_t owned_count() const;
  std::size_t ghost_count() const;
  std::size_t send_value_count() const;
  std::size_t receive_value_count() const;

  void exchange(GhostedVector& vector);
  void begin(const GhostedVector& vector);
  void wait(GhostedVector& vector);

 private:
  class Impl;
  explicit GhostedVectorHalo(std::unique_ptr<Impl> implementation) noexcept;
  std::unique_ptr<Impl> implementation_;
};

}  // namespace hundun::linear
