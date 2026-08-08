// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow_ideal_gas_closure.hpp"
#include "hundun/flow_material_density_piso.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_field_view.hpp"

#include <cmath>
#include <cstdint>
#include <limits>

namespace hundun::flow {

class IdealGasClosure;
class IdealGasClosureDiagnosticSource;

namespace detail {

class DensityClosurePreflightFailure final : public runtime::Error {
public:
  explicit DensityClosurePreflightFailure(int failing_rank)
      : runtime::Error("density closure local preflight failed"),
        failing_rank_(failing_rank) {}

  int failing_rank() const noexcept { return failing_rank_; }

private:
  int failing_rank_;
};

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
  void (*prepare_attempt)(void *){};
  void (*begin_allocation_observation)(void *) noexcept {};
  void (*end_allocation_observation)(void *) noexcept {};
  int (*outer_failure)(void *, DensityClosureOuterPoint){};
  void (*after_halo)(void *, DensityClosureStage, runtime::FieldId,
                     runtime::FieldId){};
  void (*before_post_assessment)(void *, FlowState &){};
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
      MaterialDensityTransportSpec
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      ,
      int material_factory_construction_fault_rank
#endif
  );

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
  static bool material_report_authenticated_raw(
      const MaterialDensityStepAttemptReport &) noexcept;
  static void corrupt_material_report_raw(MaterialDensityStepAttemptReport &,
                                          std::uint8_t);
  static void force_finalization_identity_wrap(FixedStepMaterialDensityFlow &);
  static void force_transport_finalization_identity_wrap_raw(
      MaterialDensityTransport &) noexcept;
  static void set_post_assessment_fault(FixedStepMaterialDensityFlow &,
                                        std::uint8_t kind, int rank);
#endif
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  static bool post_evidence_mutation_rejected(
      const MaterialDensityStepAttemptReport &report, std::uint8_t mutation) {
    auto copy = report;
    constexpr std::uint8_t report_field_count = 29U;
    constexpr std::uint8_t post_authority_begin = report_field_count;
    constexpr std::uint8_t post_coherent_begin = 2U * report_field_count;
    constexpr std::uint8_t pre_report_begin = 3U * report_field_count;
    constexpr std::uint8_t pre_authority_begin = 4U * report_field_count;
    constexpr std::uint8_t pre_coherent_begin = 5U * report_field_count;
    constexpr std::uint8_t outer_begin = 6U * report_field_count;
    constexpr std::uint8_t outer_count = 20U;
    if (!copy.material_report_ || !copy.pre_closure_authority_ ||
        !copy.post_closure_report_ || !copy.post_closure_authority_ ||
        mutation >= outer_begin + outer_count)
      return false;

    const auto next = [](double value) noexcept {
      return std::nextafter(value, std::numeric_limits<double>::infinity());
    };
    const auto change_size = [](auto &values) {
      if (values.empty())
        values.push_back({});
      else
        values.pop_back();
    };
    const auto mutate_report = [&](MaterialDensityTransportReport &target,
                                   std::uint8_t field) {
      switch (field) {
      case 0U:
        target.disposition_ =
            target.disposition_ == MaterialTransportDisposition::finalized
                ? MaterialTransportDisposition::recoverable_failure
                : MaterialTransportDisposition::finalized;
        break;
      case 1U:
        target.reason_ = target.reason_ == MaterialTransportFailureReason::none
                             ? MaterialTransportFailureReason::invalid_input
                             : MaterialTransportFailureReason::none;
        break;
      case 2U:
        ++target.lowest_failing_rank_;
        break;
      case 3U:
        target.stencil_.order =
            target.stencil_.order == MomentumTimeOrder::backward_euler
                ? MomentumTimeOrder::bdf2
                : MomentumTimeOrder::backward_euler;
        break;
      case 4U:
        target.stencil_.dt_s = next(target.stencil_.dt_s);
        break;
      case 5U:
        target.stencil_.previous_dt_s = next(target.stencil_.previous_dt_s);
        break;
      case 6U:
        target.stencil_.alpha0 = next(target.stencil_.alpha0);
        break;
      case 7U:
        target.stencil_.alpha1 = next(target.stencil_.alpha1);
        break;
      case 8U:
        target.stencil_.alpha2 = next(target.stencil_.alpha2);
        break;
      case 9U:
        target.flux_provenance_ =
            target.flux_provenance_ == MaterialFluxProvenance::final_corrected
                ? MaterialFluxProvenance::predictor
                : MaterialFluxProvenance::final_corrected;
        break;
      case 10U:
        ++target.attempt_identity_;
        break;
      case 11U:
        ++target.finalization_identity_;
        break;
      case 12U:
        ++target.shared_face_mass_flux_field_;
        break;
      case 13U:
        target.density_residual_available_ =
            !target.density_residual_available_;
        break;
      case 14U:
        target.density_normalized_l2_ = next(target.density_normalized_l2_);
        break;
      case 15U:
        change_size(target.transport_residual_available_);
        break;
      case 16U:
        target.transport_residual_available_.front() ^= 1U;
        break;
      case 17U:
        change_size(target.transport_normalized_l2_);
        break;
      case 18U:
        target.transport_normalized_l2_.front() =
            next(target.transport_normalized_l2_.front());
        break;
      case 19U:
        target.mass_conservation_available_ =
            !target.mass_conservation_available_;
        break;
      case 20U:
        target.mass_relative_conservation_defect_ =
            next(target.mass_relative_conservation_defect_);
        break;
      case 21U:
        change_size(target.transport_conservation_available_);
        break;
      case 22U:
        target.transport_conservation_available_.front() ^= 1U;
        break;
      case 23U:
        change_size(target.transport_relative_conservation_defect_);
        break;
      case 24U:
        target.transport_relative_conservation_defect_.front() =
            next(target.transport_relative_conservation_defect_.front());
        break;
      case 25U:
        target.minimum_density_available_ = !target.minimum_density_available_;
        break;
      case 26U:
        target.minimum_density_kg_per_m3_ =
            next(target.minimum_density_kg_per_m3_);
        break;
      case 27U:
        ++target.minimum_density_global_cell_;
        break;
      case 28U:
        ++target.minimum_density_rank_;
        break;
      }
      target.seal();
    };

    if (mutation < report_field_count) {
      mutate_report(*copy.post_closure_report_, mutation);
    } else if (mutation < post_coherent_begin) {
      mutate_report(*copy.post_closure_authority_,
                    static_cast<std::uint8_t>(mutation - post_authority_begin));
    } else if (mutation < pre_report_begin) {
      const auto field =
          static_cast<std::uint8_t>(mutation - post_coherent_begin);
      mutate_report(*copy.post_closure_report_, field);
      mutate_report(*copy.post_closure_authority_, field);
    } else if (mutation < pre_authority_begin) {
      mutate_report(*copy.material_report_,
                    static_cast<std::uint8_t>(mutation - pre_report_begin));
    } else if (mutation < pre_coherent_begin) {
      mutate_report(*copy.pre_closure_authority_,
                    static_cast<std::uint8_t>(mutation - pre_authority_begin));
    } else if (mutation < outer_begin) {
      const auto field =
          static_cast<std::uint8_t>(mutation - pre_coherent_begin);
      mutate_report(*copy.material_report_, field);
      mutate_report(*copy.pre_closure_authority_, field);
    } else {
      const auto field = static_cast<std::uint8_t>(mutation - outer_begin);
      switch (field) {
      case 0U:
        copy.material_failure_reason_ =
            copy.material_failure_reason_ ==
                    MaterialTransportFailureReason::none
                ? MaterialTransportFailureReason::invalid_input
                : MaterialTransportFailureReason::none;
        break;
      case 1U:
        ++copy.material_field_count_;
        break;
      case 2U:
        ++copy.shared_face_mass_flux_field_;
        break;
      case 3U:
        copy.flux_provenance_ =
            copy.flux_provenance_ == MaterialFluxProvenance::final_corrected
                ? MaterialFluxProvenance::predictor
                : MaterialFluxProvenance::final_corrected;
        break;
      case 4U:
        ++copy.attempt_identity_;
        break;
      case 5U:
        ++copy.material_attempt_identity_;
        break;
      case 6U:
        ++copy.material_finalization_identity_;
        break;
      case 7U:
        change_size(copy.flow_.final_transport_normalized_l2);
        break;
      case 8U:
        copy.flow_.final_transport_normalized_l2.front() =
            next(copy.flow_.final_transport_normalized_l2.front());
        break;
      case 9U:
        change_size(copy.flow_.final_transport_relative_conservation_defect);
        break;
      case 10U:
        copy.flow_.final_transport_relative_conservation_defect.front() = next(
            copy.flow_.final_transport_relative_conservation_defect.front());
        break;
      case 11U:
        copy.mass_conservation_available_ = !copy.mass_conservation_available_;
        break;
      case 12U:
        copy.flow_.final_mass_relative_conservation_defect =
            next(copy.flow_.final_mass_relative_conservation_defect);
        break;
      case 13U:
        copy.closure_origin_ = !copy.closure_origin_;
        break;
      case 14U:
        copy.pre_closure_authority_.reset();
        break;
      case 15U:
        ++copy.pre_closure_report_seal_authority_;
        break;
      case 16U:
        copy.post_closure_evidence_available_ =
            !copy.post_closure_evidence_available_;
        break;
      case 17U:
        copy.post_closure_report_.reset();
        break;
      case 18U:
        copy.post_closure_authority_.reset();
        break;
      case 19U:
        ++copy.post_closure_report_seal_authority_;
        break;
      }
    }
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
  static void prepare_attempt(void *);
  static void begin_allocation_observation(void *) noexcept;
  static void end_allocation_observation(void *) noexcept;
  static int outer_failure(void *, DensityClosureOuterPoint);
  static void after_halo(void *, DensityClosureStage, runtime::FieldId,
                         runtime::FieldId);
  static void before_post_assessment(void *, FlowState &);
#endif
  static void before_outlet(void *, FlowState &);
  static int before_prepare(void *, FlowState &, AcceptedStepMetadata);
  static int prepare(void *);
  static void publish(void *) noexcept;
  static void rollback(void *) noexcept;
  static double gas_constant_J_per_kg_K(const IdealGasClosure &) noexcept;
  static double cp_J_per_kg_K(const IdealGasClosure &) noexcept;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  static bool same_rank_reason_precedence_raw() noexcept;
  static bool post_store_rank_marker_collision_free_raw(int) noexcept;
  static bool final_gate_rank_marker_collision_free_raw(int) noexcept;
  static bool candidate_pressure_change_rejected_raw(
      const IdealGasClosureReport &, bool) noexcept;
  static bool closure_report_authenticated_raw(
      const IdealGasClosureReport &) noexcept;
  static IdealGasClosure create_raw(
      const mesh::MeshTopology &, const mesh::MeshGeometry &,
      const boundary::BoundaryRegistry &, const runtime::MpiContext &,
      const runtime::FieldRegistry &, const FlowFieldIds &, const FlowState &,
      IdealGasClosureSpec, std::uint8_t, int);
  static int preflight_failure_rank_raw(const runtime::Error &) noexcept;
  static std::uint64_t preflight_wire_exchange_count_raw(
      const IdealGasClosure &) noexcept;
  static IdealGasClosureFailureReason create_validation_failure_reason_raw(
      const runtime::Error &) noexcept;
  static void set_restore_preparation_fault_raw(int) noexcept;
  static void set_restore_shape_fault_raw(int) noexcept;
  static IdealGasClosureReport evaluate_raw(IdealGasClosure &, FlowState &,
                                            IdealGasClosureStage);
  static void set_stage_failure_raw(IdealGasClosure &, IdealGasClosureStage,
                                    IdealGasClosureFailureReason, int);
  static void set_metric_gate_failure_raw(IdealGasClosure &, std::uint8_t,
                                          int);
  static void set_post_store_corruption_raw(IdealGasClosure &, int, bool);
  static std::vector<std::array<std::uint64_t, 3>> halo_trace_raw(
      const IdealGasClosure &);
  static void set_candidate_precedence_fault_raw(IdealGasClosure &, int);
  static void set_outer_failure_raw(IdealGasClosure &, std::uint8_t, int);
  static void set_prepare_fault_raw(IdealGasClosure &, bool, int);
  static void set_post_store_mpi_fault_raw(IdealGasClosure &, int);
  static void set_post_assessment_fault_raw(IdealGasClosure &, std::uint8_t,
                                            int);
  static void set_attempt_layout_fault_raw(IdealGasClosure &, int);
  static void set_outlet_backflow_fault_raw(IdealGasClosure &);
  static void set_facade_create_fault_raw(IdealGasClosure &, int);
  static void set_material_factory_create_fault_raw(IdealGasClosure &, int);
  static bool consume_facade_create_fault_raw(IdealGasClosure &,
                                              int) noexcept;
  static int consume_material_factory_create_fault_raw(
      IdealGasClosure &) noexcept;
  static void set_attempt_preparation_fault_raw(IdealGasClosure &,
                                                std::uint8_t, int);
  static void set_controlled_allocation_raw(IdealGasClosure &, int);
  static bool allocation_observation_active_raw(
      const IdealGasClosure &) noexcept;
  static void consume_attempt_preparation_fault_raw(IdealGasClosure &);
#endif
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
  static std::uint64_t source_generation_raw(
      const IdealGasClosureDiagnosticSource &);
  static DensityClosureDiagnosticTestState &
  test_state(const IdealGasClosureDiagnosticSource &);
#endif
};

} // namespace detail
} // namespace hundun::flow
