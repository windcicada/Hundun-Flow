// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/flow_state.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/mpi_context.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include <array>
#include <vector>
#endif

namespace hundun::boundary {
class BoundaryRegistry;
}

namespace hundun::flow {

namespace detail {
struct DensityClosureAdapter;
}

class FixedStepMaterialDensityFlow;
class FixedStepIdealGasFlow;
class IdealGasClosureDiagnosticSource;

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace test {
class IdealGasClosureTestAccess;
}
#endif

enum class IdealGasPressureMode : std::uint8_t { closed_dynamic, open_fixed };

struct IdealGasClosureSpec final {
  runtime::FieldId enthalpy_density{};
  double cp_J_per_kg_K{};
  double gas_constant_J_per_kg_K{};
  double configured_thermodynamic_pressure_pa{};
};

struct IdealGasClosureState final {
  IdealGasPressureMode mode{IdealGasPressureMode::closed_dynamic};
  double thermodynamic_pressure_pa{};
  std::optional<double> target_mass_kg;
  std::uint64_t revision{};
};

enum class IdealGasClosureDisposition : std::uint8_t {
  closed,
  recoverable_failure,
  non_retryable_failure
};

enum class IdealGasClosureStage : std::uint8_t {
  none,
  predictor,
  provisional,
  final
};

enum class IdealGasClosureFailureReason : std::uint8_t {
  none,
  invalid_input,
  non_finite_enthalpy,
  non_positive_enthalpy,
  non_finite_temperature,
  non_positive_temperature,
  denominator_breakdown,
  non_finite_pressure,
  non_positive_pressure,
  non_finite_density,
  non_positive_density,
  eos_residual,
  remap_residual,
  mass_conservation,
  enthalpy_conservation,
  collective_operation
};

class IdealGasClosureReport final {
public:
  IdealGasClosureReport(const IdealGasClosureReport &) = default;
  IdealGasClosureReport(IdealGasClosureReport &&) noexcept = default;
  IdealGasClosureReport &operator=(const IdealGasClosureReport &) = default;
  IdealGasClosureReport &operator=(IdealGasClosureReport &&) noexcept = default;

  IdealGasClosureDisposition disposition() const noexcept;
  IdealGasClosureFailureReason reason() const noexcept;
  IdealGasClosureStage stage() const noexcept;
  int lowest_failing_rank() const noexcept;
  std::uint64_t attempt_identity() const noexcept;
  std::uint32_t evaluation_count() const noexcept;
  std::uint64_t collective_count() const noexcept;
  IdealGasPressureMode pressure_mode() const noexcept;
  double configured_pressure_pa() const noexcept;
  bool candidate_pressure_available() const noexcept;
  double candidate_pressure_pa() const;
  bool target_mass_available() const noexcept;
  double target_mass_kg() const;
  bool final_metrics_available() const noexcept;
  double actual_mass_kg() const;
  double temperature_min_K() const;
  double temperature_max_K() const;
  double enthalpy_min_J_per_kg() const;
  double enthalpy_max_J_per_kg() const;
  double density_min_kg_per_m3() const;
  double density_max_kg_per_m3() const;
  double rho_remap_normalized_l2() const;
  double rho_h_remap_normalized_l2() const;
  double rho_remap_relative_conservation_defect() const;
  double rho_h_remap_relative_conservation_defect() const;
  double enthalpy_temperature_max_relative_error() const;
  double eos_max_relative_error() const;

private:
  IdealGasClosureReport() = default;
  IdealGasClosureDisposition disposition_{
      IdealGasClosureDisposition::non_retryable_failure};
  IdealGasClosureFailureReason reason_{
      IdealGasClosureFailureReason::invalid_input};
  IdealGasClosureStage stage_{IdealGasClosureStage::none};
  int lowest_failing_rank_{-1};
  std::uint64_t attempt_identity_{};
  std::uint32_t evaluation_count_{};
  std::uint64_t collective_count_{};
  IdealGasPressureMode pressure_mode_{IdealGasPressureMode::closed_dynamic};
  double configured_pressure_pa_{};
  bool candidate_pressure_available_{};
  double candidate_pressure_pa_{};
  std::optional<double> target_mass_kg_;
  bool final_metrics_available_{};
  double actual_mass_kg_{};
  double temperature_min_K_{};
  double temperature_max_K_{};
  double enthalpy_min_J_per_kg_{};
  double enthalpy_max_J_per_kg_{};
  double density_min_kg_per_m3_{};
  double density_max_kg_per_m3_{};
  double rho_remap_normalized_l2_{};
  double rho_h_remap_normalized_l2_{};
  double rho_remap_relative_conservation_defect_{};
  double rho_h_remap_relative_conservation_defect_{};
  double enthalpy_temperature_max_relative_error_{};
  double eos_max_relative_error_{};
  std::uint64_t seal_{};

  std::uint64_t compute_seal() const noexcept;
  bool semantic_valid() const noexcept;
  void seal() noexcept;
  bool authenticated() const noexcept;
  friend class IdealGasClosure;
  friend class FixedStepIdealGasFlow;
  friend class IdealGasStepAttemptReport;
  friend class IdealGasClosureDiagnosticSource;
  friend struct detail::DensityClosureAdapter;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::IdealGasClosureTestAccess;
#endif
};

class IdealGasClosure final {
public:
  static IdealGasClosure
  create(const mesh::MeshTopology &, const mesh::MeshGeometry &,
         const boundary::BoundaryRegistry &, const runtime::MpiContext &,
         const runtime::FieldRegistry &, const FlowFieldIds &,
         const FlowState &initialized_state, IdealGasClosureSpec);

  ~IdealGasClosure() noexcept;
  IdealGasClosure(IdealGasClosure &&) noexcept;
  IdealGasClosure &operator=(IdealGasClosure &&) = delete;
  IdealGasClosure(const IdealGasClosure &) = delete;
  IdealGasClosure &operator=(const IdealGasClosure &) = delete;

  IdealGasClosureState state() const;

private:
  struct Impl;
  explicit IdealGasClosure(std::unique_ptr<Impl>) noexcept;
  static IdealGasClosure
  create_internal(const mesh::MeshTopology &, const mesh::MeshGeometry &,
                  const boundary::BoundaryRegistry &,
                  const runtime::MpiContext &, const runtime::FieldRegistry &,
                  const FlowFieldIds &, const FlowState &, IdealGasClosureSpec
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
                  ,
                  int, int
#endif
  );
  void begin_attempt(const FlowState &, std::uint64_t);
  const IdealGasClosureReport &evaluate(FlowState &, IdealGasClosureStage);
  int prepare_commit();
  void publish_commit() noexcept;
  void rollback() noexcept;
  void before_outlet(FlowState &);
  int before_prepare(FlowState &, AcceptedStepMetadata);
  const IdealGasClosureReport &latest_report() const;
  bool matches(const mesh::MeshTopology &, const mesh::MeshGeometry &,
               const boundary::BoundaryRegistry &, const runtime::MpiContext &,
               const runtime::FieldRegistry &,
               const FlowFieldIds &) const noexcept;
  double cp_J_per_kg_K() const noexcept;
  double gas_constant_J_per_kg_K() const noexcept;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  void set_post_store_corruption_for_test(int rank, bool enthalpy_density);
  void set_candidate_precedence_fault_for_test(int rank);
  void set_stage_failure_for_test(IdealGasClosureStage,
                                  IdealGasClosureFailureReason, int rank);
  void set_metric_gate_failure_for_test(std::uint8_t kind, int rank);
  void set_post_assessment_fault_for_test(std::uint8_t kind, int rank);
  void before_post_assessment_for_test(FlowState &);
  void record_halo_for_test(std::uint8_t stage, runtime::FieldId density,
                            runtime::FieldId enthalpy_density);
  std::vector<std::array<std::uint64_t, 3>> halo_trace_for_test() const;
  void set_outer_failure_for_test(std::uint8_t point, int rank);
  int outer_failure_for_test(std::uint8_t point);
  void set_prepare_fault_for_test(bool state_prepare, int rank);
  void set_post_store_mpi_fault_for_test(int rank);
  void set_attempt_layout_fault_for_test(int rank);
  void set_outlet_backflow_fault_for_test();
#endif
  std::unique_ptr<Impl> impl_;

  friend class FixedStepMaterialDensityFlow;
  friend class FixedStepIdealGasFlow;
  friend class IdealGasStepAttemptReport;
  friend class IdealGasClosureDiagnosticSource;
  friend struct detail::DensityClosureAdapter;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::IdealGasClosureTestAccess;
#endif
};

} // namespace hundun::flow
