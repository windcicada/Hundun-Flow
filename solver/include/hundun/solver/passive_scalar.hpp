// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/field_descriptor.hpp"
#include "hundun/runtime/types.hpp"

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
