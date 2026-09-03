// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_types.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string_view>
#include <vector>

namespace hundun::runtime {
class StructuredDecomposition;
}

namespace hundun::mesh {

using GlobalCellId = std::uint64_t;
using GlobalFaceId = std::uint64_t;
using LocalCellId = std::size_t;
using LocalFaceId = std::size_t;

enum class FaceAxis : std::uint8_t { x = 0, y = 1, z = 2 };
enum class EntityOwnership : std::uint8_t { owned, ghost };
enum class PatchPairingKind : std::uint8_t { none, periodic };

struct LogicalFace {
  FaceAxis axis{};
  runtime::Int3 coordinate{};
};

class MeshTopology;

class BoundaryPatch final {
 public:
  BoundaryPatch(const BoundaryPatch&) = delete;
  BoundaryPatch& operator=(const BoundaryPatch&) = delete;
  BoundaryPatch(BoundaryPatch&&) noexcept = default;
  BoundaryPatch& operator=(BoundaryPatch&&) noexcept = default;

  std::uint32_t stable_id() const noexcept;
  std::string_view name() const noexcept;
  PatchPairingKind pairing_kind() const noexcept;
  std::optional<std::uint32_t> paired_patch_id() const noexcept;
  const std::vector<LocalFaceId>& local_faces() const noexcept;
  bool contains(LocalFaceId local_face) const noexcept;

 private:
  friend class MeshTopology;

  BoundaryPatch(std::uint32_t stable_id, std::string_view name,
                PatchPairingKind pairing_kind,
                std::optional<std::uint32_t> paired_patch_id,
                std::vector<LocalFaceId> local_faces);

  std::uint32_t stable_id_{};
  std::string_view name_{};
  PatchPairingKind pairing_kind_{PatchPairingKind::none};
  std::optional<std::uint32_t> paired_patch_id_{};
  std::vector<LocalFaceId> local_faces_{};
};

class MeshTopology final {
 public:
  explicit MeshTopology(const runtime::StructuredDecomposition& decomposition);
  ~MeshTopology();
  MeshTopology(const MeshTopology&) = delete;
  MeshTopology& operator=(const MeshTopology&) = delete;
  MeshTopology(MeshTopology&&) noexcept;
  MeshTopology& operator=(MeshTopology&&) noexcept;

  runtime::Int3 global_extent() const noexcept;
  runtime::Box3 owned_global_box() const noexcept;
  std::uint64_t global_cell_count() const noexcept;
  std::uint64_t global_face_count() const noexcept;
  std::uint64_t global_face_count(FaceAxis axis) const;

  std::size_t owned_cell_count() const noexcept;
  std::size_t ghost_cell_count() const noexcept;
  std::size_t local_cell_count() const noexcept;
  std::size_t owned_face_count() const noexcept;
  std::size_t ghost_face_count() const noexcept;
  std::size_t local_face_count() const noexcept;

  GlobalCellId global_cell_id(runtime::Int3 global_cell) const;
  GlobalCellId global_cell_id(LocalCellId local_cell) const;
  runtime::Int3 global_cell(LocalCellId local_cell) const;
  EntityOwnership cell_ownership(LocalCellId local_cell) const;
  std::optional<LocalCellId> find_local_cell(GlobalCellId global_id) const;

  GlobalFaceId global_face_id(LogicalFace face) const;
  GlobalFaceId global_face_id(LocalFaceId local_face) const;
  LogicalFace logical_face(LocalFaceId local_face) const;
  EntityOwnership face_ownership(LocalFaceId local_face) const;
  LocalCellId owner(LocalFaceId local_face) const;
  std::optional<LocalCellId> neighbour(LocalFaceId local_face) const;
  std::optional<LocalFaceId> find_local_face(GlobalFaceId global_id) const;
  std::optional<std::uint32_t> patch_id(LocalFaceId local_face) const;
  std::optional<GlobalFaceId> periodic_pair(LocalFaceId local_face) const;

  const std::array<BoundaryPatch, 6>& patches() const noexcept;
  const BoundaryPatch& patch(std::uint32_t stable_id) const;

 private:
  struct Impl;
  static std::unique_ptr<Impl> create_impl(
      const runtime::StructuredDecomposition& decomposition);

  std::unique_ptr<Impl> impl_;
};

}  // namespace hundun::mesh
