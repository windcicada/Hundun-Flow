// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/diagnostics/structured_diagnostics.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/finite_volume/matrix_free_poisson.hpp"
#include "hundun/flow/adaptive_time_control.hpp"
#include "hundun/linear/ghosted_vector.hpp"
#include "hundun/linear/ghosted_vector_halo.hpp"
#include "hundun/linear/linear_system.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/structured_decomposition.hpp"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::diagnostics {

struct FieldRoleDiagnostic final {
  std::string role;
  runtime::FieldId field{};
};

struct FieldLayoutDiagnosticSource final {
  const runtime::FieldRegistry* registry{};
  runtime::FieldLayoutSet layout{};
  std::vector<FieldRoleDiagnostic> roles;
};

struct LinearSolveDiagnosticSource final {
  std::string_view instance_id;
  const linear::SolveReport* report{};
};

struct SharedFluxDiagnosticSource final {
  runtime::FieldId field{};
  std::size_t face_count{};
  bool final_flux{};
};

struct ConstantDensityPisoDiagnosticSource final {
  const flow::StepAttemptReport* report{};
};

struct FlowDriverDiagnosticSource final {
  config::DensityModel density_model{config::DensityModel::constant};
  std::uint64_t step{};
  double time_s{};
  std::size_t attempt_count{};
  flow::TimeAdvanceDisposition disposition{
      flow::TimeAdvanceDisposition::non_retryable_failure};
  flow::StepFailureReason reason{flow::StepFailureReason::invalid_input};
  int lowest_failing_rank{-1};
};

#define HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(Type)                       \
  DiagnosticDescriptor describe_diagnostics(const Type&) noexcept;         \
  std::vector<std::string_view> diagnostic_fingerprint_field_ids(           \
      const Type&);                                                         \
  void collect_diagnostics(const Type&, const DiagnosticRequest&,           \
                           DiagnosticSink&)

HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(runtime::MpiContext);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(runtime::StructuredDecomposition);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(FieldLayoutDiagnosticSource);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(runtime::ExchangePlan);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(mesh::MeshTopology);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(mesh::MeshGeometry);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(execution::ExecutionContext);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(execution::Buffer);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(execution::VectorView<const double>);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(linear::GhostedVector);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(linear::GhostedVectorHalo);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(linear::LinearOperator);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(
    finite_volume::MatrixFreePoissonOperator);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(LinearSolveDiagnosticSource);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(SharedFluxDiagnosticSource);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(boundary::BoundaryRegistry);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(ConstantDensityPisoDiagnosticSource);
HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER(FlowDriverDiagnosticSource);

#undef HUNDUN_DECLARE_LOCAL_DIAGNOSTIC_ADAPTER

void collect_diagnostics(const runtime::MpiContext&,
                         const runtime::MpiContext& collective_context,
                         const DiagnosticRequest&, DiagnosticSink&);

}  // namespace hundun::diagnostics
