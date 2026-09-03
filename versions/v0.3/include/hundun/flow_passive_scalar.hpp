// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_types.hpp"

namespace hundun::mesh {
class UniformStructuredMesh;
}

namespace hundun::runtime {
class FieldStorage;
class HaloExchange;
class MpiContext;
class StructuredDecomposition;
} // namespace hundun::runtime

namespace hundun::solver {

[[nodiscard]] double mc_limiter(double left, double right) noexcept;

class PassiveScalarSolver final {
public:
  // Borrows context and halo for the solver lifetime; neither object may be
  // destroyed or moved first. Construction validates halo against the
  // supplied decomposition and retains no field view.
  PassiveScalarSolver(const runtime::MpiContext &context,
                      const runtime::StructuredDecomposition &decomposition,
                      const mesh::UniformStructuredMesh &mesh,
                      runtime::HaloExchange &halo,
                      runtime::Real3 velocity_m_per_s,
                      double diffusivity_m2_per_s);
  void advance_ssprk2(runtime::FieldStorage &storage, runtime::FieldId scalar,
                      runtime::FieldId stage, double dt_s);

private:
  const runtime::MpiContext *context_{};
  runtime::HaloExchange *halo_{};
  runtime::Int3 local_extent_{};
  runtime::Real3 spacing_m_{};
  runtime::Real3 velocity_m_per_s_{};
  int halo_ghost_width_{};
};

[[nodiscard]] double global_mass(const runtime::MpiContext &context,
                                 const mesh::UniformStructuredMesh &mesh,
                                 const runtime::FieldStorage &storage,
                                 runtime::FieldId scalar);
[[nodiscard]] double global_l1_error(const runtime::MpiContext &context,
                                     const mesh::UniformStructuredMesh &mesh,
                                     const runtime::FieldStorage &storage,
                                     runtime::FieldId actual,
                                     runtime::FieldId reference);

} // namespace hundun::solver
