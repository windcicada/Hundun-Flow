// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/mesh_topology.hpp"
#include "hundun/rt_types.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>

namespace hundun::mesh {

enum class MappingKind : std::uint8_t {
  uniform_box,
  analytic_warped_box,
};

enum class FaceSide : std::uint8_t { owner, neighbour };

struct MappingJacobian {
  runtime::Real3 d_xi_m{};
  runtime::Real3 d_eta_m{};
  runtime::Real3 d_zeta_m{};

  double determinant_m3() const noexcept;
};

class UniformBoxMapping final {
 public:
  UniformBoxMapping(runtime::Real3 origin_m, runtime::Real3 length_m);

  runtime::Real3 origin_m() const noexcept;
  runtime::Real3 length_m() const noexcept;
  runtime::Real3 map(runtime::Real3 logical) const;
  MappingJacobian jacobian(runtime::Real3 logical) const;

 private:
  runtime::Real3 origin_m_{};
  runtime::Real3 length_m_{};
};

class AnalyticWarpedBoxMapping final {
 public:
  AnalyticWarpedBoxMapping(runtime::Real3 origin_m, runtime::Real3 length_m,
                           runtime::Real3 amplitude);

  runtime::Real3 origin_m() const noexcept;
  runtime::Real3 length_m() const noexcept;
  runtime::Real3 amplitude() const noexcept;
  runtime::Real3 map(runtime::Real3 logical) const;
  MappingJacobian jacobian(runtime::Real3 logical) const;

 private:
  runtime::Real3 origin_m_{};
  runtime::Real3 length_m_{};
  runtime::Real3 amplitude_{};
};

class MeshGeometry final {
 public:
  MeshGeometry(const MeshTopology&, UniformBoxMapping);
  MeshGeometry(const MeshTopology&, AnalyticWarpedBoxMapping);
  ~MeshGeometry();
  MeshGeometry(const MeshGeometry&) = delete;
  MeshGeometry& operator=(const MeshGeometry&) = delete;
  MeshGeometry(MeshGeometry&&) noexcept;
  MeshGeometry& operator=(MeshGeometry&&) noexcept;

  MappingKind mapping_kind() const noexcept;
  runtime::Int3 global_extent() const noexcept;
  runtime::Box3 owned_global_box() const noexcept;
  runtime::Real3 origin_m() const noexcept;
  runtime::Real3 length_m() const noexcept;
  std::optional<runtime::Real3> uniform_spacing_m() const noexcept;
  std::size_t local_face_count() const noexcept;

  bool compatible(const MeshTopology&) const;
  void require_compatible(const MeshTopology&) const;

  runtime::Real3 vertex_position_m(runtime::Int3 global_vertex) const;
  runtime::Real3 cell_center_m(LocalCellId) const;
  double cell_volume_m3(LocalCellId) const;
  double minimum_jacobian_determinant_m3(LocalCellId) const;
  runtime::Real3 face_center_m(LocalFaceId) const;
  runtime::Real3 face_displacement_m(LocalFaceId) const;
  runtime::Real3 face_area_vector_m2(LocalFaceId, FaceSide) const;
  double face_area_m2(LocalFaceId) const;
  double face_skewness(LocalFaceId) const;
  double face_non_orthogonality_degrees(LocalFaceId) const;
  runtime::Real3 cell_closure_m2(LocalCellId owned_cell) const;

 private:
  struct Impl;
  static std::unique_ptr<Impl> create_impl(const MeshTopology&,
                                           UniformBoxMapping);
  static std::unique_ptr<Impl> create_impl(const MeshTopology&,
                                           AnalyticWarpedBoxMapping);

  std::unique_ptr<Impl> impl_;
};

}  // namespace hundun::mesh
