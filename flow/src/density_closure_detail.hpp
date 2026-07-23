// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/material_density_piso.hpp"
#include "hundun/runtime/field_view.hpp"

#include <cmath>
#include <cstdint>
#include <limits>

namespace hundun::flow {

class IdealGasClosure;
class IdealGasClosureDiagnosticSource;

namespace detail {

enum class DensityClosureStage : std::uint8_t { predictor, provisional, final };
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
enum class DensityClosureOuterPoint : std::uint8_t {
  momentum_after_predictor,
  pressure_after_first_corrector
};
#endif

struct DensityClosureEvaluation final {
  bool accepted{};
  bool recoverable{};
  int lowest_failing_rank{-1};
};

struct DensityClosureHooks final {
  void *object{};
  runtime::FieldId enthalpy_density{};
  double enthalpy_rate_J_per_kg_s{};
  void (*begin)(void *, const FlowState &, std::uint64_t){};
  DensityClosureEvaluation (*evaluate)(void *, FlowState &,
                                       DensityClosureStage){};
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  int (*outer_failure)(void *, DensityClosureOuterPoint){};
#endif
  void (*before_outlet)(void *, FlowState &){};
  int (*before_prepare)(void *, FlowState &, AcceptedStepMetadata){};
  int (*prepare)(void *){};
  void (*publish)(void *) noexcept {};
  void (*rollback)(void *) noexcept {};

  explicit operator bool() const noexcept { return object != nullptr; }
};

struct DensityClosureBridge final {
  static FixedStepMaterialDensityFlow create_open_capable(
      const runtime::StructuredDecomposition &, const mesh::MeshTopology &,
      const mesh::MeshGeometry &, const boundary::BoundaryRegistry &,
      const runtime::MpiContext &, execution::ExecutionContext &,
      runtime::HaloExchange &, const linear::LinearSolver &,
      std::array<linear::Preconditioner *, 3>, const linear::LinearSolver &,
      linear::Preconditioner &, const runtime::FieldRegistry &, FlowFieldIds,
      MaterialDensityTransportSpec);

  static MaterialDensityStepAttemptReport
  attempt(const FixedStepMaterialDensityFlow &, FlowState &, double,
          const MomentumTimeStencil &, const linear::SolveControl &,
          const linear::SolveControl &, const DensityClosureHooks &);

  static bool
  report_authenticated(const MaterialDensityStepAttemptReport &) noexcept;
  static std::uint64_t
  report_seal(const MaterialDensityStepAttemptReport &) noexcept;
  static MaterialDensityStepAttemptReport make_report();
  static bool post_eos_evidence_authenticated(
      const MaterialDensityStepAttemptReport &) noexcept;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  static bool post_evidence_mutation_rejected(
      const MaterialDensityStepAttemptReport &report, std::uint8_t mutation) {
    auto copy = report;
    auto &post = *copy.post_closure_report_;
    switch (mutation) {
    case 0U:
      ++post.attempt_identity_;
      break;
    case 1U:
      ++post.shared_face_mass_flux_field_;
      break;
    case 2U:
      post.flux_provenance_ = MaterialFluxProvenance::predictor;
      break;
    case 3U:
      post.stencil_.dt_s = std::nextafter(
          post.stencil_.dt_s, std::numeric_limits<double>::infinity());
      break;
    case 4U:
      post.finalization_identity_ = 0U;
      break;
    case 5U:
      post.finalization_identity_ = copy.material_finalization_identity_;
      break;
    case 6U:
      post.transport_residual_available_.pop_back();
      break;
    case 7U:
      post.transport_normalized_l2_.front() = std::nextafter(
          post.transport_normalized_l2_.front(),
          std::numeric_limits<double>::infinity());
      break;
    default:
      return false;
    }
    post.seal();
    copy.seal_ = copy.compute_seal();
    return !report_authenticated(copy);
  }
#endif
};

struct DensityClosureAdapter final {
  static DensityClosureHooks bind(IdealGasClosure &, runtime::FieldId,
                                  double enthalpy_rate_J_per_kg_s);
  static void begin(void *, const FlowState &, std::uint64_t);
  static DensityClosureEvaluation evaluate(void *, FlowState &,
                                           DensityClosureStage);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  static int outer_failure(void *, DensityClosureOuterPoint);
#endif
  static void before_outlet(void *, FlowState &);
  static int before_prepare(void *, FlowState &, AcceptedStepMetadata);
  static int prepare(void *);
  static void publish(void *) noexcept;
  static void rollback(void *) noexcept;
  static double gas_constant_J_per_kg_K(const IdealGasClosure &) noexcept;
  static double cp_J_per_kg_K(const IdealGasClosure &) noexcept;
};

struct DensityClosureReadSession final {
  runtime::FieldView<const double> density;
  runtime::FieldView<const double> enthalpy_density;
  runtime::Box3 owned_box{};
  runtime::Int3 global_extent{};
};

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
struct DensityClosureDiagnosticTestState final {
  std::uint8_t fault{};
  int fault_rank{-1};
  bool consumed{};
  std::uint64_t observations{};
  std::uint64_t fingerprint_items{};
  std::uint64_t summary_items{};
  std::uint64_t invariant_items{};
  std::uint64_t sample_items{};
  std::uint64_t retained_sample_items{};
  std::uint64_t allocation_events{};
  std::uint64_t full_field_copy_attempts{};
  std::uint64_t collective_calls{};
};
#endif

struct DensityClosureDiagnosticAccess final {
  static DensityClosureReadSession
  acquire_committed(const IdealGasClosureDiagnosticSource &);
  static double
  gas_constant_J_per_kg_K(const IdealGasClosureDiagnosticSource &);
  static double cp_J_per_kg_K(const IdealGasClosureDiagnosticSource &);
  static const runtime::MpiContext &
  mpi(const IdealGasClosureDiagnosticSource &);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  static DensityClosureDiagnosticTestState &
  test_state(const IdealGasClosureDiagnosticSource &);
#endif
};

} // namespace detail
} // namespace hundun::flow
