// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/fvm_cell_centered.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_storage.hpp"

#include <cstdint>
#include <memory>

namespace hundun::boundary {
class BoundaryRegistry;
}

namespace hundun::immersed {
class GhostStencilPlan;
class ImmersedDomain;
} // namespace hundun::immersed

namespace hundun::mesh {
class MeshGeometry;
class MeshTopology;
} // namespace hundun::mesh

namespace hundun::runtime {
class HaloExchange;
class MpiContext;
class StructuredDecomposition;
} // namespace hundun::runtime

#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
namespace hundun::finite_volume::test {
class ImmersedReconstructionTestAccess;
}
#endif

namespace hundun::finite_volume {

class ImmersedOperatorAdapter;

struct ReconstructionFieldBinding final {
  runtime::FieldStorage &storage;
  const runtime::FieldAccessPlan &access_plan;
  runtime::PhaseId phase;
  runtime::ActorId actor;
  runtime::FieldId field;
};

class ImmersedReconstruction final {
public:
  static ImmersedReconstruction
  create(const mesh::MeshTopology &, const mesh::MeshGeometry &,
         const boundary::BoundaryRegistry &, const immersed::ImmersedDomain &,
         const immersed::GhostStencilPlan &,
         const runtime::StructuredDecomposition &, const runtime::MpiContext &,
         runtime::HaloExchange &);

  ~ImmersedReconstruction() noexcept;
  ImmersedReconstruction(ImmersedReconstruction &&) noexcept;
  ImmersedReconstruction &operator=(ImmersedReconstruction &&) = delete;
  ImmersedReconstruction(const ImmersedReconstruction &) = delete;
  ImmersedReconstruction &operator=(const ImmersedReconstruction &) = delete;

  void compute_gradient(GradientScheme, FiniteVolumeQuantity,
                        const ReconstructionFieldBinding &cell_values,
                        const ReconstructionFieldBinding &cell_gradients) const;
  void reconstruct_transport_faces(
      FiniteVolumeQuantity, const FaceMassFlux &,
      const ReconstructionFieldBinding &cell_values,
      const ReconstructionFieldBinding &face_values) const;
  void reconstruct_momentum_faces(
      const FaceMassFlux &, const ReconstructionFieldBinding &velocity,
      const ReconstructionFieldBinding &face_velocity) const;
  std::uint64_t dependency_fingerprint() const noexcept;

private:
  struct Impl;
  explicit ImmersedReconstruction(std::unique_ptr<Impl>) noexcept;
  Impl &require_impl() const;
  void require_immersed_operator_compatible(
      const mesh::MeshTopology &, const mesh::MeshGeometry &,
      const immersed::ImmersedDomain &,
      const immersed::GhostStencilPlan &) const;
  const runtime::MpiContext &mpi_for_immersed_operator() const;

  std::unique_ptr<Impl> impl_;
  friend class ImmersedOperatorAdapter;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  friend class test::ImmersedReconstructionTestAccess;
#endif
};

} // namespace hundun::finite_volume
