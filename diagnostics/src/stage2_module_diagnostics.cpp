// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/stage2_module_diagnostics.hpp"

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace hundun::diagnostics {
namespace {

constexpr DiagnosticCapabilityFlags kSummary =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary);
constexpr DiagnosticCapabilityFlags kSummaryCounters =
    kSummary |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters);
constexpr DiagnosticCapabilityFlags kSummaryInvariants =
    kSummary |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants);

DiagnosticDescriptor descriptor(DiagnosticModuleKind kind,
                                std::string_view module_id,
                                DiagnosticCapabilityFlags capabilities) noexcept {
  return {kDiagnosticRecordSchemaV1, kind, module_id, module_id, capabilities};
}

DiagnosticRecord base_record(const DiagnosticDescriptor& descriptor,
                             const DiagnosticRequest& request) {
  DiagnosticRecord record;
  record.module_kind = descriptor.module_kind;
  record.module_id = std::string(descriptor.module_id);
  record.instance_id = std::string(descriptor.instance_id);
  record.level = request.level;
  record.scope = request.scope;
  record.rank = request.frame.rank;
  record.step = request.frame.step;
  record.time_s = describe_fp64(request.frame.time_s);
  record.phase = std::string(request.frame.phase);
  return record;
}

DiagnosticMetric metric(std::string id, std::string unit, double value) {
  return {std::move(id), DiagnosticMetricKind::state_summary, std::move(unit),
          describe_fp64(value)};
}

DiagnosticCounter counter(std::string id, std::uint64_t value,
                          std::string unit = "count") {
  return {std::move(id), std::move(unit), value};
}

DiagnosticInvariant finite_invariant(std::string id, std::string unit,
                                     double value) {
  DiagnosticInvariant result{std::move(id),
                             std::move(unit),
                             describe_fp64(value),
                             {},
                             InvariantRelation::finite,
                             false};
  result.passed = evaluate_invariant(result);
  return result;
}

DiagnosticInvariant positive_invariant(std::string id, std::string unit,
                                       double value) {
  DiagnosticInvariant result{std::move(id),
                             std::move(unit),
                             describe_fp64(value),
                             {},
                             InvariantRelation::positive,
                             false};
  result.passed = evaluate_invariant(result);
  return result;
}

void finish_record(DiagnosticRecord& record) {
  DiagnosticFingerprintAccumulator fingerprint;
  std::uint64_t ordinal = 0;
  for (const auto& invariant : record.invariants)
    fingerprint.add(invariant.id, ordinal++, 0U, invariant.observed);
  for (const auto& value : record.metrics)
    fingerprint.add(value.id, ordinal++, 0U, value.value);
  for (const auto& value : record.counters)
    fingerprint.add(value.id, ordinal++, 0U,
                    describe_fp64(static_cast<double>(value.value)));
  for (const auto& value : record.identities) {
    fingerprint.add(value.subject_id, ordinal++, 0U,
                    describe_fp64(static_cast<double>(
                        value.revision.value_or(value.generation.value_or(
                            value.allocation_identity.value_or(0U))))));
  }
  record.state_fingerprint = fingerprint.finish();
}

template <class Fill>
void collect_common(const DiagnosticDescriptor& descriptor,
                    const DiagnosticRequest& request, DiagnosticSink& sink,
                    Fill&& fill) {
  validate(request, descriptor);
  DiagnosticRecord record = base_record(descriptor, request);
  std::forward<Fill>(fill)(record);
  finish_record(record);
  validate(record, descriptor, request);
  try {
    sink.submit(record);
  } catch (const DiagnosticCollectionError&) {
    throw;
  } catch (const std::exception& error) {
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", -1,
                                    error.what());
  }
}

void require_summary(const DiagnosticRequest& request) {
  if (request.level != DiagnosticLevel::summary)
    throw DiagnosticCollectionError(DiagnosticFailureClass::capability,
                                    "diagnostics.request.capability", -1,
                                    "unsupported diagnostic request");
}

template <class Integer>
double as_double(Integer value) {
  return static_cast<double>(value);
}

std::uint64_t checked_u64(std::size_t value) {
  return static_cast<std::uint64_t>(value);
}

std::string failure_code(flow::StepFailureReason reason) {
  switch (reason) {
  case flow::StepFailureReason::none:
    return "none";
  case flow::StepFailureReason::invalid_input:
    return "flow.invalid-input";
  case flow::StepFailureReason::momentum_linear_solve:
    return "flow.momentum-linear-solve";
  case flow::StepFailureReason::pressure_linear_solve:
    return "flow.pressure-linear-solve";
  case flow::StepFailureReason::non_finite_trial:
    return "flow.non-finite-trial";
  case flow::StepFailureReason::boundary_backflow:
    return "flow.boundary-backflow";
  case flow::StepFailureReason::transport_failure:
    return "flow.transport-failure";
  case flow::StepFailureReason::final_momentum_residual:
    return "flow.final-momentum-residual";
  case flow::StepFailureReason::final_transport_residual:
    return "flow.final-transport-residual";
  case flow::StepFailureReason::final_conservation_defect:
    return "flow.final-conservation-defect";
  case flow::StepFailureReason::final_continuity_residual:
    return "flow.final-continuity-residual";
  case flow::StepFailureReason::final_pressure_residual:
    return "flow.final-pressure-residual";
  case flow::StepFailureReason::collective_operation:
    return "flow.collective-operation";
  case flow::StepFailureReason::density_closure_failure:
    return "flow.density-closure-failure";
  }
  return "flow.invalid-input";
}

DiagnosticFailureClass failure_class(flow::StepFailureReason reason) {
  switch (reason) {
  case flow::StepFailureReason::none:
    return DiagnosticFailureClass::none;
  case flow::StepFailureReason::invalid_input:
    return DiagnosticFailureClass::invalid_input;
  case flow::StepFailureReason::non_finite_trial:
    return DiagnosticFailureClass::non_finite_state;
  case flow::StepFailureReason::boundary_backflow:
    return DiagnosticFailureClass::boundary;
  case flow::StepFailureReason::collective_operation:
    return DiagnosticFailureClass::collective_operation;
  case flow::StepFailureReason::momentum_linear_solve:
  case flow::StepFailureReason::pressure_linear_solve:
  case flow::StepFailureReason::transport_failure:
  case flow::StepFailureReason::density_closure_failure:
    return DiagnosticFailureClass::numerical_breakdown;
  case flow::StepFailureReason::final_momentum_residual:
  case flow::StepFailureReason::final_transport_residual:
  case flow::StepFailureReason::final_continuity_residual:
  case flow::StepFailureReason::final_pressure_residual:
    return DiagnosticFailureClass::non_convergence;
  case flow::StepFailureReason::final_conservation_defect:
    return DiagnosticFailureClass::conservation;
  }
  return DiagnosticFailureClass::invalid_input;
}

void set_failure(DiagnosticRecord& record, flow::StepFailureReason reason,
                 int lowest_rank) {
  if (reason == flow::StepFailureReason::none)
    return;
  record.status = DiagnosticStatus::failed;
  record.failure = {failure_class(reason), failure_code(reason),
                    record.scope == DiagnosticScope::collective
                        ? std::max(0, lowest_rank)
                        : -1};
}

}  // namespace

DiagnosticDescriptor describe_diagnostics(const runtime::MpiContext&) noexcept {
  return descriptor(DiagnosticModuleKind::mpi, "hundun.runtime.mpi_context",
                    kSummaryCounters |
                        static_cast<DiagnosticCapabilityFlags>(
                            DiagnosticCapability::collective));
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const runtime::MpiContext&) {
  return {"rank", "size", "thread_level"};
}
void collect_diagnostics(const runtime::MpiContext& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(describe_diagnostics(source), request, sink,
                 [&](DiagnosticRecord& record) {
                   if (request.scope != DiagnosticScope::local)
                     throw DiagnosticCollectionError(
                         DiagnosticFailureClass::capability,
                         "diagnostics.request.collective-context", -1,
                         "collective context is required");
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("rank", "count", as_double(source.rank())),
                         metric("size", "count", as_double(source.size())),
                         metric("thread_level", "count",
                                as_double(source.thread_level()))};
                   } else {
                     const auto values = source.fp64_reduction_counters();
                     record.counters = {
                         counter("collective_calls", values.collective_calls),
                         counter("logical_payload_bytes",
                                 values.logical_payload_bytes, "byte"),
                         counter("reduced_scalars", values.reduced_scalars)};
                   }
                 });
}
void collect_diagnostics(const runtime::MpiContext& source,
                         const runtime::MpiContext& collective_context,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  validate(request, describe_diagnostics(source));
  if (request.scope != DiagnosticScope::collective)
    throw DiagnosticCollectionError(DiagnosticFailureClass::invalid_request,
                                    "diagnostics.request.scope", -1,
                                    "collective overload requires collective scope");
  const std::array<int, 3> local{source.rank(), source.size(),
                                 source.thread_level()};
  std::array<int, 3> minimum = local;
  std::array<int, 3> maximum = local;
  runtime::check_mpi_result(
      MPI_Allreduce(local.data(), minimum.data(), 3, MPI_INT, MPI_MIN,
                    collective_context.comm()),
      "MPI_Allreduce MPI diagnostic minimum");
  runtime::check_mpi_result(
      MPI_Allreduce(local.data(), maximum.data(), 3, MPI_INT, MPI_MAX,
                    collective_context.comm()),
      "MPI_Allreduce MPI diagnostic maximum");
  collect_common(describe_diagnostics(source), request, sink,
                 [&](DiagnosticRecord& record) {
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("rank_max", "count", as_double(maximum[0])),
                         metric("rank_min", "count", as_double(minimum[0])),
                         metric("size", "count", as_double(source.size()))};
                   } else {
                     const auto values = source.fp64_reduction_counters();
                     record.counters = {
                         counter("collective_calls", values.collective_calls),
                         counter("logical_payload_bytes",
                                 values.logical_payload_bytes, "byte"),
                         counter("reduced_scalars", values.reduced_scalars)};
                   }
                 });
}

#define HUNDUN_SIMPLE_SUMMARY_ADAPTER(Type, Kind, Id, Body)                 \
  DiagnosticDescriptor describe_diagnostics(const Type&) noexcept {         \
    return descriptor(DiagnosticModuleKind::Kind, Id, kSummary);            \
  }                                                                         \
  std::vector<std::string_view> diagnostic_fingerprint_field_ids(           \
      const Type&) {                                                        \
    return {};                                                              \
  }                                                                         \
  void collect_diagnostics(const Type& source,                              \
                           const DiagnosticRequest& request,                \
                           DiagnosticSink& sink) {                           \
    collect_common(describe_diagnostics(source), request, sink,             \
                   [&](DiagnosticRecord& record) {                           \
                     require_summary(request);                              \
                     record.metrics = Body;                                 \
                   });                                                      \
  }

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    runtime::StructuredDecomposition, runtime,
    "hundun.runtime.structured_decomposition",
    (std::vector<DiagnosticMetric>{
        metric("global_cell_count", "count",
               as_double(static_cast<std::uint64_t>(source.global_extent().x) *
                         static_cast<std::uint64_t>(source.global_extent().y) *
                         static_cast<std::uint64_t>(source.global_extent().z))),
        metric("local_cell_count", "count",
               as_double(static_cast<std::uint64_t>(source.local_extent().x) *
                         static_cast<std::uint64_t>(source.local_extent().y) *
                         static_cast<std::uint64_t>(source.local_extent().z))),
        metric("process_grid_x", "count",
               as_double(source.process_grid().x)),
        metric("process_grid_y", "count",
               as_double(source.process_grid().y)),
        metric("process_grid_z", "count",
               as_double(source.process_grid().z))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    runtime::ExchangePlan, halo, "hundun.runtime.halo",
    (std::vector<DiagnosticMetric>{
        metric("ghost_width", "count", as_double(source.ghost_width())),
        metric("region_count", "count", as_double(source.regions().size()))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    mesh::MeshTopology, mesh_topology, "hundun.mesh.topology",
    (std::vector<DiagnosticMetric>{
        metric("global_cell_count", "count",
               as_double(source.global_cell_count())),
        metric("global_face_count", "count",
               as_double(source.global_face_count())),
        metric("local_cell_count", "count",
               as_double(source.local_cell_count())),
        metric("local_face_count", "count",
               as_double(source.local_face_count())),
        metric("owned_cell_count", "count",
               as_double(source.owned_cell_count())),
        metric("owned_face_count", "count",
               as_double(source.owned_face_count()))}))

DiagnosticDescriptor
describe_diagnostics(const FieldLayoutDiagnosticSource&) noexcept {
  return descriptor(DiagnosticModuleKind::field, "hundun.runtime.field_layout",
                    kSummary);
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const FieldLayoutDiagnosticSource&) {
  return {"cell_count", "face_count", "field_count"};
}
void collect_diagnostics(const FieldLayoutDiagnosticSource& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(describe_diagnostics(source), request, sink,
                 [&](DiagnosticRecord& record) {
                   require_summary(request);
                   if (source.registry == nullptr ||
                       !source.registry->frozen())
                     throw DiagnosticCollectionError(
                         DiagnosticFailureClass::invalid_input,
                         "field-layout.invalid", -1,
                         "field layout source is not frozen");
                   record.metrics = {
                       metric("cell_count", "count",
                              as_double(static_cast<std::uint64_t>(
                                  source.layout.cell_interior_extent.x) *
                                        static_cast<std::uint64_t>(
                                            source.layout.cell_interior_extent.y) *
                                        static_cast<std::uint64_t>(
                                            source.layout.cell_interior_extent.z))),
                       metric("face_count", "count",
                              as_double(source.layout.face_count)),
                       metric("field_count", "count",
                              as_double(source.registry->size()))};
                 });
}

DiagnosticDescriptor
describe_diagnostics(const mesh::MeshGeometry&) noexcept {
  return descriptor(DiagnosticModuleKind::mesh_geometry,
                    "hundun.mesh.geometry", kSummaryInvariants);
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const mesh::MeshGeometry&) {
  return {"cell_volume_sum", "maximum_closure_norm",
          "minimum_jacobian"};
}
void collect_diagnostics(const mesh::MeshGeometry& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(describe_diagnostics(source), request, sink,
                 [&](DiagnosticRecord& record) {
                   const auto box = source.owned_global_box();
                   const std::size_t count =
                       static_cast<std::size_t>(box.end.x - box.begin.x) *
                       static_cast<std::size_t>(box.end.y - box.begin.y) *
                       static_cast<std::size_t>(box.end.z - box.begin.z);
                   double volume_sum = 0.0;
                   double minimum_jacobian =
                       std::numeric_limits<double>::infinity();
                   double maximum_closure = 0.0;
                   for (std::size_t cell = 0; cell < count; ++cell) {
                     volume_sum += source.cell_volume_m3(cell);
                     minimum_jacobian =
                         std::min(minimum_jacobian,
                                  source.minimum_jacobian_determinant_m3(cell));
                     const auto closure = source.cell_closure_m2(cell);
                     maximum_closure =
                         std::max(maximum_closure,
                                  std::hypot(closure.x,
                                             std::hypot(closure.y, closure.z)));
                   }
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("cell_volume_sum", "m3", volume_sum),
                         metric("maximum_closure_norm", "m2",
                                maximum_closure),
                         metric("minimum_jacobian", "m3",
                                minimum_jacobian)};
                   } else {
                     record.invariants = {
                         finite_invariant("cell_closure_finite", "m2",
                                          maximum_closure),
                         positive_invariant("cell_volume_sum_positive", "m3",
                                            volume_sum),
                         positive_invariant("jacobian_positive", "m3",
                                            minimum_jacobian)};
                   }
                 });
}

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    execution::ExecutionContext, execution, "hundun.execution.context",
    (std::vector<DiagnosticMetric>{
        metric("backend_identity", "count",
               as_double(source.backend_identity())),
        metric("ordered", "count", source.ordered() ? 1.0 : 0.0),
        metric("space", "count",
               source.space() == execution::ExecutionSpace::host ? 0.0
                                                                 : 1.0)}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    execution::Buffer, execution, "hundun.execution.buffer",
    (std::vector<DiagnosticMetric>{
        metric("allocation_identity", "count",
               as_double(source.allocation_identity())),
        metric("byte_size", "byte", as_double(source.byte_size())),
        metric("epoch", "count", as_double(source.epoch()))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    execution::VectorView<const double>, execution,
    "hundun.execution.vector_view",
    (std::vector<DiagnosticMetric>{
        metric("allocation_identity", "count",
               as_double(source.allocation_identity())),
        metric("element_count", "count", as_double(source.size())),
        metric("epoch", "count", as_double(source.epoch())),
        metric("offset_bytes", "byte", as_double(source.offset_bytes())),
        metric("stride", "count", as_double(source.stride()))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    linear::GhostedVector, execution, "hundun.linear.ghosted_vector",
    (std::vector<DiagnosticMetric>{
        metric("allocation_identity", "count",
               as_double(source.allocation_identity())),
        metric("epoch", "count", as_double(source.epoch())),
        metric("ghost_count", "count", as_double(source.ghost_count())),
        metric("local_count", "count", as_double(source.local_count())),
        metric("owned_count", "count", as_double(source.owned_count()))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    linear::GhostedVectorHalo, halo, "hundun.linear.ghosted_vector_halo",
    (std::vector<DiagnosticMetric>{
        metric("ghost_count", "count", as_double(source.ghost_count())),
        metric("owned_count", "count", as_double(source.owned_count())),
        metric("path", "count", as_double(static_cast<int>(source.path()))),
        metric("receive_value_count", "count",
               as_double(source.receive_value_count())),
        metric("send_value_count", "count",
               as_double(source.send_value_count()))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    linear::LinearOperator, linear_operator, "hundun.linear.operator",
    (std::vector<DiagnosticMetric>{
        metric("diagonal_available", "count",
               source.has_diagonal() ? 1.0 : 0.0),
        metric("domain_local_count", "count",
               as_double(source.domain_layout().local_count())),
        metric("range_local_count", "count",
               as_double(source.range_layout().local_count())),
        metric("revision", "count", as_double(source.revision()))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    finite_volume::MatrixFreePoissonOperator, linear_operator,
    "hundun.finite_volume.poisson",
    (std::vector<DiagnosticMetric>{
        metric("constraint_mode", "count",
               as_double(static_cast<int>(source.constraint_mode()))),
        metric("diagonal_available", "count",
               source.has_diagonal() ? 1.0 : 0.0),
        metric("revision", "count", as_double(source.revision())),
        metric("solver_family", "count",
               as_double(static_cast<int>(source.solver_family())))}))

DiagnosticDescriptor
describe_diagnostics(const LinearSolveDiagnosticSource& source) noexcept {
  return {kDiagnosticRecordSchemaV1,
          DiagnosticModuleKind::linear_solver,
          "hundun.linear.solve",
          source.instance_id,
          kSummaryCounters};
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const LinearSolveDiagnosticSource&) {
  return {"iterations", "termination"};
}
void collect_diagnostics(const LinearSolveDiagnosticSource& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(describe_diagnostics(source), request, sink,
                 [&](DiagnosticRecord& record) {
                   if (source.report == nullptr || source.instance_id.empty())
                     throw DiagnosticCollectionError(
                         DiagnosticFailureClass::invalid_input,
                         "linear-solve.invalid", -1,
                         "linear solve source is invalid");
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("iterations", "count",
                                as_double(source.report->iterations)),
                         metric("termination", "count",
                                as_double(static_cast<int>(
                                    source.report->reason)))};
                   } else {
                     record.counters = {
                         counter("global_reduction_count",
                                 source.report->global_reduction_count),
                         counter("matvec_count",
                                 source.report->matvec_count),
                         counter("preconditioner_apply_count",
                                 source.report->preconditioner_apply_count)};
                   }
                 });
}

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    SharedFluxDiagnosticSource, finite_volume,
    "hundun.finite_volume.shared_flux",
    (std::vector<DiagnosticMetric>{
        metric("face_count", "count", as_double(source.face_count)),
        metric("field_id", "count",
               as_double(static_cast<std::uint64_t>(source.field))),
        metric("final_flux", "count", source.final_flux ? 1.0 : 0.0)}))

DiagnosticDescriptor
describe_diagnostics(const boundary::BoundaryRegistry&) noexcept {
  return descriptor(DiagnosticModuleKind::boundary,
                    "hundun.boundary.registry", kSummaryCounters);
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const boundary::BoundaryRegistry&) {
  return {"open_domain", "scalar_count"};
}
void collect_diagnostics(const boundary::BoundaryRegistry& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(describe_diagnostics(source), request, sink,
                 [&](DiagnosticRecord& record) {
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("open_domain", "count",
                                source.open_domain() ? 1.0 : 0.0),
                         metric("scalar_count", "count",
                                as_double(source.scalar_count()))};
                   } else {
                     std::uint64_t periodic = 0;
                     std::uint64_t walls = 0;
                     for (std::uint32_t patch = 0; patch < 6U; ++patch) {
                       const auto kind = source.patch(patch).kind();
                       periodic += kind == boundary::BoundaryKind::periodic;
                       walls += kind == boundary::BoundaryKind::no_slip_wall;
                     }
                     record.counters = {
                         counter("no_slip_wall_count", walls),
                         counter("patch_count", 6U),
                         counter("periodic_patch_count", periodic)};
                   }
                 });
}

DiagnosticDescriptor describe_diagnostics(
    const ConstantDensityPisoDiagnosticSource&) noexcept {
  return descriptor(DiagnosticModuleKind::piso,
                    "hundun.flow.constant_density_piso",
                    kSummaryCounters);
}
std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const ConstantDensityPisoDiagnosticSource&) {
  return {"corrector_count", "disposition", "reason"};
}
void collect_diagnostics(const ConstantDensityPisoDiagnosticSource& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(describe_diagnostics(source), request, sink,
                 [&](DiagnosticRecord& record) {
                   if (source.report == nullptr)
                     throw DiagnosticCollectionError(
                         DiagnosticFailureClass::invalid_input,
                         "constant-piso.invalid", -1,
                         "constant density PISO source is invalid");
                   const auto& report = *source.report;
                   set_failure(record, report.reason,
                               report.lowest_failing_rank);
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("continuity_normalized_l2", "1",
                                report.final_continuity_normalized_l2),
                         metric("pressure_residual_l2", "1",
                                report.final_pressure_residual_l2)};
                   } else {
                     record.counters = {
                         counter("pressure_corrector_count",
                                 report.pressure_corrector_count),
                         counter("transported_field_count",
                                 checked_u64(
                                     report.final_transport_normalized_l2
                                         .size()))};
                   }
                 });
}

DiagnosticDescriptor
describe_diagnostics(const FlowDriverDiagnosticSource&) noexcept {
  return descriptor(DiagnosticModuleKind::flow_driver,
                    "hundun.application.flow_driver", kSummaryCounters);
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const FlowDriverDiagnosticSource&) {
  return {"attempt_count", "density_model", "step", "time"};
}
void collect_diagnostics(const FlowDriverDiagnosticSource& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(describe_diagnostics(source), request, sink,
                 [&](DiagnosticRecord& record) {
                   set_failure(record, source.reason,
                               source.lowest_failing_rank);
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("density_model", "count",
                                as_double(static_cast<int>(
                                    source.density_model))),
                         metric("time", "s", source.time_s)};
                   } else {
                     record.counters = {
                         counter("attempt_count",
                                 checked_u64(source.attempt_count)),
                         counter("step", source.step)};
                   }
                 });
}

#undef HUNDUN_SIMPLE_SUMMARY_ADAPTER

}  // namespace hundun::diagnostics
