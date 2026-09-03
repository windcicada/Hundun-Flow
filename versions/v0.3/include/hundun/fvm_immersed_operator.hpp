// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/fvm_cell_centered.hpp"

#include <array>
#include <cstdint>
#include <memory>

namespace hundun::immersed {
class GhostStencilPlan;
class ImmersedDomain;
class LocalFlowPatternTransform;
} // namespace hundun::immersed

namespace hundun::mesh {
class MeshGeometry;
class MeshTopology;
} // namespace hundun::mesh

#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
namespace hundun::finite_volume::test {
class ImmersedOperatorTestAccess;
}
#endif

namespace hundun::finite_volume {

namespace detail {
class ImmersedBoundaryAuthorityAccess;
}

class ImmersedReconstruction;

struct ImmersedResidualParts final {
  std::array<double, 3> pressure{};
  std::array<double, 3> viscous{};
  std::array<double, 3> convective{};
};

struct ImmersedOperatorReport final {
  std::uint64_t active_row_count{};
  std::uint64_t replacement_group_count{};
  std::uint64_t simultaneous_substitution_count{};
  std::uint64_t row_fingerprint{};
  std::uint64_t algebraic_occurrence_count{};
  double replacement_coefficient_l2{};
  std::uint64_t limiting_case_status{};
  ImmersedResidualParts budget_reaction_N{};
};

class ImmersedOperatorAdapter final {
public:
  static ImmersedOperatorAdapter
  create(const mesh::MeshTopology &, const mesh::MeshGeometry &,
         const immersed::ImmersedDomain &, const immersed::GhostStencilPlan &,
         const immersed::LocalFlowPatternTransform &,
         const ImmersedReconstruction &);

  ~ImmersedOperatorAdapter() noexcept;
  ImmersedOperatorAdapter(ImmersedOperatorAdapter &&) noexcept;
  ImmersedOperatorAdapter &operator=(ImmersedOperatorAdapter &&) = delete;
  ImmersedOperatorAdapter(const ImmersedOperatorAdapter &) = delete;
  ImmersedOperatorAdapter &operator=(const ImmersedOperatorAdapter &) = delete;

  void accumulate_momentum(
      const FaceMassFlux &,
      const runtime::FaceFieldView<const double> &face_velocity,
      const runtime::FieldView<const double> &velocity,
      const runtime::FieldView<const double> &pressure,
      const runtime::FieldView<const double> &velocity_gradient,
      const runtime::FaceFieldView<const double> &dynamic_viscosity_by_face,
      const runtime::FieldView<double> &residual) const;

  void accumulate_transport(
      FiniteVolumeQuantity, const FaceMassFlux &,
      const runtime::FaceFieldView<const double> &face_values,
      const runtime::FieldView<const double> &values,
      const runtime::FieldView<const double> &gradients,
      const runtime::FaceFieldView<const double> &gamma_by_face,
      const runtime::FieldView<double> &residual) const;

  ImmersedOperatorReport report() const;

private:
  struct Impl;
  explicit ImmersedOperatorAdapter(std::unique_ptr<Impl>) noexcept;
  Impl &require_impl() const;
  std::unique_ptr<Impl> impl_;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  friend class test::ImmersedOperatorTestAccess;
#endif
  friend class detail::ImmersedBoundaryAuthorityAccess;
};

} // namespace hundun::finite_volume
