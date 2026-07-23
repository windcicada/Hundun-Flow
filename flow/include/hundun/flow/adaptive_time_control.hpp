// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/config/resolved_case.hpp"
#include "hundun/flow/ideal_gas_piso.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <variant>

namespace hundun::diagnostics::detail {
struct TimeControlAdapter;
}

namespace hundun::flow {

namespace detail {
struct AdaptiveTimeControlAccess;
struct AdaptiveTimeControlEngine;
}
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace test {
class AdaptiveTimeControlTestAccess;
}
#endif

enum class TimeAdvanceDisposition : std::uint8_t {
  committed,
  retry_limit_reached,
  minimum_dt_failure,
  non_retryable_failure
};

struct TimeControlState final {
  std::uint32_t schema_version{1};
  std::uint64_t accepted_step{};
  double proposed_next_dt_s{};
  double last_accepted_dt_s{};
  MomentumTimeOrder last_accepted_order{MomentumTimeOrder::backward_euler};
  bool history_ready{};
  bool last_all_linear_solves_within_half_limit{};
  double last_convective_rate_per_s{};
  double last_diffusive_rate_per_s{};
  bool last_stability_metrics_available{};
  std::uint32_t last_retry_count{};
  std::uint64_t revision{};
  std::uint64_t state_seal{};
};

struct TimeAttemptSummary final {
  double attempted_dt_s{};
  MomentumTimeOrder order{MomentumTimeOrder::backward_euler};
  StepAttemptDisposition disposition{
      StepAttemptDisposition::non_retryable_failure};
  StepFailureReason reason{StepFailureReason::invalid_input};
  int lowest_failing_rank{-1};
  bool all_linear_solves_within_half_limit{};
};

using DensityStepAttemptReport =
    std::variant<StepAttemptReport, MaterialDensityStepAttemptReport,
                 IdealGasStepAttemptReport>;

class TimeAdvanceReport final {
public:
  ~TimeAdvanceReport() noexcept;
  TimeAdvanceReport(TimeAdvanceReport &&) noexcept;
  TimeAdvanceReport &operator=(TimeAdvanceReport &&) = delete;
  TimeAdvanceReport(const TimeAdvanceReport &) = delete;
  TimeAdvanceReport &operator=(const TimeAdvanceReport &) = delete;

  TimeAdvanceDisposition disposition() const noexcept;
  StepFailureReason reason() const noexcept;
  int lowest_failing_rank() const noexcept;
  std::size_t attempt_count() const noexcept;
  const TimeAttemptSummary &attempt(std::size_t) const;
  bool final_attempt_available() const noexcept;
  const DensityStepAttemptReport &final_attempt() const;
  double accepted_dt_s() const noexcept;
  double proposed_next_dt_s() const noexcept;
  double convective_rate_per_s() const noexcept;
  double diffusive_rate_per_s() const noexcept;
  bool stability_metrics_available() const noexcept;
  bool limited_by_min_dt() const noexcept;

private:
  enum class PreflightCategory : std::uint8_t {
    none,
    config,
    identity,
    layout,
    capability,
    state,
    transport_authority,
    preparation,
    report
  };
  struct ConstantReportTag final {};
  struct MaterialReportTag final {};
  struct IdealGasReportTag final {};
  explicit TimeAdvanceReport(ConstantReportTag) noexcept;
  explicit TimeAdvanceReport(MaterialReportTag) noexcept;
  explicit TimeAdvanceReport(IdealGasReportTag) noexcept;
  static bool report_authenticated(
      const StepAttemptReport &) noexcept;
  static bool report_authenticated(
      const MaterialDensityStepAttemptReport &) noexcept;
  static bool report_authenticated(
      const IdealGasStepAttemptReport &) noexcept;
  void reset_moved_from() noexcept;

  std::array<TimeAttemptSummary, 9> attempts_{};
  std::size_t attempt_count_{};
  DensityStepAttemptReport final_attempt_;
  bool final_attempt_available_{};
  TimeAdvanceDisposition disposition_{
      TimeAdvanceDisposition::non_retryable_failure};
  StepFailureReason reason_{StepFailureReason::invalid_input};
  int lowest_failing_rank_{-1};
  double accepted_dt_s_{};
  double proposed_next_dt_s_{};
  double convective_rate_per_s_{};
  double diffusive_rate_per_s_{};
  bool stability_metrics_available_{};
  bool limited_by_min_dt_{};
  std::uint64_t controller_identity_{};
  std::uint64_t report_identity_{};
  std::uint64_t observed_flow_state_identity_{};
  std::uint64_t observed_step_{};
  double observed_time_s_{};
  AcceptedStepMetadata observed_metadata_{};
  PreflightCategory preflight_category_{PreflightCategory::none};
  friend class Bdf2RetryController;
  friend struct detail::AdaptiveTimeControlEngine;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::AdaptiveTimeControlTestAccess;
#endif
  friend struct ::hundun::diagnostics::detail::TimeControlAdapter;
};

class TimeControlDiagnosticSource final {
public:
  ~TimeControlDiagnosticSource() noexcept;
  TimeControlDiagnosticSource(TimeControlDiagnosticSource &&) noexcept;
  TimeControlDiagnosticSource &operator=(TimeControlDiagnosticSource &&) =
      delete;
  TimeControlDiagnosticSource(const TimeControlDiagnosticSource &) = delete;
  TimeControlDiagnosticSource &
  operator=(const TimeControlDiagnosticSource &) = delete;

private:
  struct Impl;
  explicit TimeControlDiagnosticSource(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend class Bdf2RetryController;
  friend struct ::hundun::diagnostics::detail::TimeControlAdapter;
};

class Bdf2RetryController final {
public:
  static Bdf2RetryController
  create(const config::FlowTimeConfig &, config::DensityModel,
         const mesh::MeshTopology &, const mesh::MeshGeometry &,
         const runtime::MpiContext &, const FlowState &);
  static Bdf2RetryController
  restore(const config::FlowTimeConfig &, config::DensityModel,
          const mesh::MeshTopology &, const mesh::MeshGeometry &,
          const runtime::MpiContext &, const FlowState &,
          const TimeControlState &);

  ~Bdf2RetryController() noexcept;
  Bdf2RetryController(Bdf2RetryController &&) noexcept;
  Bdf2RetryController &operator=(Bdf2RetryController &&) = delete;
  Bdf2RetryController(const Bdf2RetryController &) = delete;
  Bdf2RetryController &operator=(const Bdf2RetryController &) = delete;

  TimeControlState state() const noexcept;
  TimeAdvanceReport
  advance(FlowState &, FixedStepConstantDensityFlow &, double rho_ref,
          double mu, const linear::SolveControl &momentum,
          const linear::SolveControl &pressure);
  TimeAdvanceReport
  advance(FlowState &, FixedStepMaterialDensityFlow &, double mu,
          const linear::SolveControl &momentum,
          const linear::SolveControl &pressure);
  TimeAdvanceReport
  advance(FlowState &, FixedStepIdealGasFlow &, double mu,
          const linear::SolveControl &momentum,
          const linear::SolveControl &pressure);
  TimeControlDiagnosticSource
  diagnostic_source(const FlowState &, const TimeAdvanceReport &) const;

private:
  struct Impl;
  explicit Bdf2RetryController(std::unique_ptr<Impl>,
                               TimeControlState observer_state) noexcept;
  std::unique_ptr<Impl> impl_;
  TimeControlState observer_state_{};
  friend struct detail::AdaptiveTimeControlEngine;
};

} // namespace hundun::flow
