// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/adaptive_time_control.hpp"

#include "adaptive_time_control_detail.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <type_traits>
#include <utility>

namespace hundun::flow {
namespace {

constexpr std::uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
constexpr char kSealDomain[] = "hundun-time-control-state-seal-v1";

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool positive_zero(double value) noexcept { return bits(value) == 0U; }

void byte(std::uint64_t &hash, std::uint8_t value) noexcept {
  hash = (hash ^ value) * kFnvPrime;
}
template <class UInt>
void little_endian(std::uint64_t &hash, UInt value) noexcept {
  for (std::size_t i = 0; i < sizeof(UInt); ++i)
    byte(hash, static_cast<std::uint8_t>(value >> (8U * i)));
}
void fp64(std::uint64_t &hash, double value) noexcept {
  little_endian(hash, bits(value));
}

bool valid_config(const config::FlowTimeConfig &config) noexcept {
  return config.steps >= 0 && std::isfinite(config.initial_dt_s) &&
         std::isfinite(config.min_dt_s) && std::isfinite(config.max_dt_s) &&
         config.initial_dt_s > 0.0 && config.min_dt_s > 0.0 &&
         config.max_dt_s >= config.min_dt_s &&
         config.initial_dt_s >= config.min_dt_s &&
         config.initial_dt_s <= config.max_dt_s &&
         config.cfl_target == 0.5 &&
         config.diffusion_number_target == 0.25 &&
         config.growth_factor == 1.25 && config.retry_factor == 0.5 &&
         config.max_retries == 8;
}

double expected_proposal(const config::FlowTimeConfig &config,
                         const TimeControlState &state) noexcept {
  if (config.mode == config::TimeMode::fixed)
    return config.initial_dt_s;
  double result =
      std::min(config.max_dt_s,
               state.last_accepted_dt_s *
                   (state.last_all_linear_solves_within_half_limit
                        ? config.growth_factor
                        : 1.0));
  if (state.last_convective_rate_per_s > 0.0)
    result =
        std::min(result, config.cfl_target / state.last_convective_rate_per_s);
  if (state.last_diffusive_rate_per_s > 0.0)
    result = std::min(
        result,
        config.diffusion_number_target / state.last_diffusive_rate_per_s);
  return std::max(result, config.min_dt_s);
}

StepAttemptReport const &base_report(const StepAttemptReport &report) {
  return report;
}
StepAttemptReport const &
base_report(const MaterialDensityStepAttemptReport &report) {
  return report.flow();
}
StepAttemptReport const &base_report(const IdealGasStepAttemptReport &report) {
  return report.flow().flow();
}

bool half_work_gate(const StepAttemptReport &report,
                    const linear::SolveControl &momentum,
                    const linear::SolveControl &pressure) noexcept {
  if (report.disposition != StepAttemptDisposition::committed)
    return false;
  const auto momentum_limit = momentum.max_iterations / 2U;
  const auto pressure_limit = pressure.max_iterations / 2U;
  return std::all_of(report.momentum.components.begin(),
                     report.momentum.components.end(),
                     [=](const auto &solve) {
                       return solve.iterations <= momentum_limit;
                     }) &&
         std::all_of(report.pressure.begin(), report.pressure.end(),
                     [=](const auto &solve) {
                       return solve.iterations <= pressure_limit;
                     });
}

std::string ranked_message(std::string_view prefix, int rank) {
  return std::string(prefix) + ": lowest failing rank " +
         std::to_string(rank);
}

} // namespace

detail::TransportDiffusivityAuthority
detail::make_transport_diffusivity_authority(
    const std::vector<double> &coefficients) {
  TransportDiffusivityAuthority result;
  result.count = coefficients.size();
  for (double value : coefficients) {
    std::uint64_t value_bits = bits(value);
    for (std::size_t i = 0; i < sizeof(value_bits); ++i) {
      result.ordered_fingerprint ^=
          static_cast<std::uint8_t>(value_bits >> (8U * i));
      result.ordered_fingerprint *= kFnvPrime;
    }
    if (!std::isfinite(value) || value < 0.0)
      result.maximum = std::numeric_limits<double>::quiet_NaN();
    else if (std::isfinite(result.maximum))
      result.maximum = std::max(result.maximum, value);
  }
  return result;
}

std::uint64_t detail::TimeControlStateCodec::seal(
    const config::FlowTimeConfig &config, config::DensityModel model,
    const TimeControlState &state) noexcept {
  std::uint64_t hash = kFnvOffset;
  for (char ch : std::string_view(kSealDomain))
    byte(hash, static_cast<std::uint8_t>(ch));
  little_endian(hash, static_cast<std::uint32_t>(config.mode));
  fp64(hash, config.initial_dt_s);
  fp64(hash, config.min_dt_s);
  fp64(hash, config.max_dt_s);
  fp64(hash, config.cfl_target);
  fp64(hash, config.diffusion_number_target);
  fp64(hash, config.growth_factor);
  fp64(hash, config.retry_factor);
  little_endian(hash, static_cast<std::uint32_t>(config.max_retries));
  little_endian(hash, static_cast<std::uint32_t>(model));
  little_endian(hash, state.schema_version);
  little_endian(hash, state.accepted_step);
  fp64(hash, state.proposed_next_dt_s);
  fp64(hash, state.last_accepted_dt_s);
  byte(hash, static_cast<std::uint8_t>(state.last_accepted_order));
  byte(hash, state.history_ready ? 1U : 0U);
  byte(hash, state.last_all_linear_solves_within_half_limit ? 1U : 0U);
  fp64(hash, state.last_convective_rate_per_s);
  fp64(hash, state.last_diffusive_rate_per_s);
  byte(hash, state.last_stability_metrics_available ? 1U : 0U);
  little_endian(hash, state.last_retry_count);
  little_endian(hash, state.revision);
  return hash;
}

bool detail::TimeControlStateCodec::semantically_valid(
    const config::FlowTimeConfig &config, config::DensityModel model,
    const AcceptedStepMetadata &metadata,
    const TimeControlState &state) noexcept {
  if (!valid_config(config) || state.schema_version != 1U ||
      state.accepted_step != metadata.step ||
      state.revision != state.accepted_step ||
      state.history_ready != (state.accepted_step != 0U) ||
      state.last_retry_count > 8U ||
      !std::isfinite(state.proposed_next_dt_s) ||
      state.proposed_next_dt_s < config.min_dt_s ||
      state.proposed_next_dt_s > config.max_dt_s ||
      !std::isfinite(metadata.time_s) || metadata.time_s < 0.0 ||
      state.state_seal != seal(config, model, state))
    return false;
  if (state.accepted_step == 0U)
    return positive_zero(metadata.time_s) && positive_zero(metadata.previous_dt_s) &&
           metadata.order == MomentumTimeOrder::backward_euler &&
           bits(metadata.dt_s) == bits(config.initial_dt_s) &&
           positive_zero(state.last_accepted_dt_s) &&
           state.last_accepted_order == MomentumTimeOrder::backward_euler &&
           !state.last_all_linear_solves_within_half_limit &&
           !state.last_stability_metrics_available &&
           positive_zero(state.last_convective_rate_per_s) &&
           positive_zero(state.last_diffusive_rate_per_s) &&
           state.last_retry_count == 0U &&
           bits(state.proposed_next_dt_s) == bits(config.initial_dt_s);
  if (!(metadata.time_s > 0.0) || !(metadata.dt_s > 0.0) ||
      !(metadata.previous_dt_s > 0.0) ||
      bits(state.last_accepted_dt_s) != bits(metadata.dt_s) ||
      state.last_accepted_order != metadata.order)
    return false;
  const double ratio = metadata.dt_s / metadata.previous_dt_s;
  const auto expected_order =
      state.accepted_step == 1U || ratio < 0.5 || ratio > 2.0
          ? MomentumTimeOrder::backward_euler
          : MomentumTimeOrder::bdf2;
  if (metadata.order != expected_order)
    return false;
  if (config.mode == config::TimeMode::fixed)
    return !state.last_stability_metrics_available &&
           positive_zero(state.last_convective_rate_per_s) &&
           positive_zero(state.last_diffusive_rate_per_s) &&
           bits(state.proposed_next_dt_s) == bits(config.initial_dt_s);
  return state.last_stability_metrics_available &&
         std::isfinite(state.last_convective_rate_per_s) &&
         state.last_convective_rate_per_s >= 0.0 &&
         std::isfinite(state.last_diffusive_rate_per_s) &&
         state.last_diffusive_rate_per_s >= 0.0 &&
         bits(state.proposed_next_dt_s) == bits(expected_proposal(config, state));
}

std::string detail::render_owned_cells(runtime::Box3 box) {
  return "cell.f64.owned." + std::to_string(box.begin.x) + "." +
         std::to_string(box.begin.y) + "." + std::to_string(box.begin.z) +
         "." + std::to_string(box.end.x) + "." + std::to_string(box.end.y) +
         "." + std::to_string(box.end.z);
}
std::string detail::render_global_cells(runtime::Int3 extent) {
  return "cell.f64.global." + std::to_string(extent.x) + "." +
         std::to_string(extent.y) + "." + std::to_string(extent.z);
}
std::string detail::render_owned_faces(std::size_t count) {
  return "face.f64.owned." + std::to_string(count);
}
std::string detail::render_global_faces(std::uint64_t count) {
  return "face.f64.global." + std::to_string(count);
}

TimeAdvanceReport::TimeAdvanceReport(ConstantReportTag) noexcept
    : final_attempt_(StepAttemptReport{}) {}
TimeAdvanceReport::TimeAdvanceReport(MaterialReportTag) noexcept
    : final_attempt_(MaterialDensityStepAttemptReport{}) {}
TimeAdvanceReport::TimeAdvanceReport(IdealGasReportTag) noexcept
    : final_attempt_(IdealGasStepAttemptReport{}) {}
TimeAdvanceReport::~TimeAdvanceReport() noexcept = default;
TimeAdvanceReport::TimeAdvanceReport(TimeAdvanceReport &&other) noexcept
    : attempts_(other.attempts_), attempt_count_(other.attempt_count_),
      final_attempt_(std::move(other.final_attempt_)),
      final_attempt_available_(other.final_attempt_available_),
      disposition_(other.disposition_), reason_(other.reason_),
      lowest_failing_rank_(other.lowest_failing_rank_),
      accepted_dt_s_(other.accepted_dt_s_),
      proposed_next_dt_s_(other.proposed_next_dt_s_),
      convective_rate_per_s_(other.convective_rate_per_s_),
      diffusive_rate_per_s_(other.diffusive_rate_per_s_),
      stability_metrics_available_(other.stability_metrics_available_),
      limited_by_min_dt_(other.limited_by_min_dt_),
      controller_identity_(other.controller_identity_),
      report_identity_(other.report_identity_),
      preflight_category_(other.preflight_category_) {
  other.reset_moved_from();
}
void TimeAdvanceReport::reset_moved_from() noexcept {
  attempts_ = {};
  attempt_count_ = 0U;
  final_attempt_available_ = false;
  disposition_ = TimeAdvanceDisposition::non_retryable_failure;
  reason_ = StepFailureReason::invalid_input;
  lowest_failing_rank_ = -1;
  accepted_dt_s_ = proposed_next_dt_s_ = convective_rate_per_s_ =
      diffusive_rate_per_s_ = 0.0;
  stability_metrics_available_ = limited_by_min_dt_ = false;
  controller_identity_ = report_identity_ = 0U;
  preflight_category_ = PreflightCategory::none;
}
TimeAdvanceDisposition TimeAdvanceReport::disposition() const noexcept {
  return disposition_;
}
StepFailureReason TimeAdvanceReport::reason() const noexcept { return reason_; }
int TimeAdvanceReport::lowest_failing_rank() const noexcept {
  return lowest_failing_rank_;
}
std::size_t TimeAdvanceReport::attempt_count() const noexcept {
  return attempt_count_;
}
const TimeAttemptSummary &TimeAdvanceReport::attempt(std::size_t index) const {
  if (index >= attempt_count_)
    throw runtime::Error("time advance attempt index is out of range");
  return attempts_[index];
}
bool TimeAdvanceReport::final_attempt_available() const noexcept {
  return final_attempt_available_;
}
const DensityStepAttemptReport &TimeAdvanceReport::final_attempt() const {
  if (!final_attempt_available_)
    throw runtime::Error("time advance final attempt is unavailable");
  return final_attempt_;
}
double TimeAdvanceReport::accepted_dt_s() const noexcept {
  return accepted_dt_s_;
}
double TimeAdvanceReport::proposed_next_dt_s() const noexcept {
  return proposed_next_dt_s_;
}
double TimeAdvanceReport::convective_rate_per_s() const noexcept {
  return convective_rate_per_s_;
}
double TimeAdvanceReport::diffusive_rate_per_s() const noexcept {
  return diffusive_rate_per_s_;
}
bool TimeAdvanceReport::stability_metrics_available() const noexcept {
  return stability_metrics_available_;
}
bool TimeAdvanceReport::limited_by_min_dt() const noexcept {
  return limited_by_min_dt_;
}

TimeControlDiagnosticSource::TimeControlDiagnosticSource(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
TimeControlDiagnosticSource::~TimeControlDiagnosticSource() noexcept = default;
TimeControlDiagnosticSource::TimeControlDiagnosticSource(
    TimeControlDiagnosticSource &&) noexcept = default;

struct Bdf2RetryController::Impl final {
  config::FlowTimeConfig config;
  config::DensityModel model;
  const mesh::MeshTopology *topology{};
  const mesh::MeshGeometry *geometry{};
  const runtime::MpiContext *mpi{};
  const FlowState *state{};
  std::uint64_t identity{1U};
  std::uint64_t report_identity{};
  bool active{};
};

namespace {

TimeControlState fresh_state(const config::FlowTimeConfig &config,
                             config::DensityModel model) noexcept {
  TimeControlState state;
  state.proposed_next_dt_s = config.initial_dt_s;
  state.state_seal = detail::TimeControlStateCodec::seal(config, model, state);
  return state;
}

void collective_factory_check(const runtime::MpiContext &mpi, bool local_ok,
                              std::string_view prefix) {
  const auto status = runtime::collective_status(mpi, local_ok, prefix);
  if (!status.ok)
    throw runtime::Error(ranked_message(prefix, status.failing_rank));
}

bool controls_valid(const linear::SolveControl &control) noexcept {
  return std::isfinite(control.atol) && control.atol >= 0.0 &&
         std::isfinite(control.rtol) && control.rtol >= 0.0 &&
         control.max_iterations > 0U;
}

bool agrees_u64(const runtime::MpiContext &mpi, std::uint64_t local,
                std::string_view operation) {
  std::uint64_t root = local;
  runtime::check_mpi_result(
      MPI_Bcast(&root, 1, MPI_UINT64_T, 0, mpi.comm()), operation);
  return root == local;
}

struct Stability final {
  double convection{};
  double diffusion{};
};

template <class ControllerImpl>
Stability stability_rates(const ControllerImpl &impl,
                          const FlowState &state, double density_constant,
                          bool variable_density, double gamma) {
  const auto density = detail::AdaptiveTimeControlAccess::committed_density(state);
  const auto flux =
      detail::AdaptiveTimeControlAccess::committed_face_mass_flux(state);
  const auto cells = impl.topology->owned_cell_count();
  if (density.size() != cells ||
      flux.size() != impl.topology->local_face_count())
    throw runtime::Error("time-control stability layout is invalid");
  std::vector<double> flux_sum(cells);
  std::vector<double> geometry_sum(cells);
  for (mesh::LocalFaceId face = 0; face < impl.topology->local_face_count();
       ++face) {
    const auto owner = impl.topology->owner(face);
    const double mass_flux = flux[face];
    const auto area = impl.geometry->face_area_m2(face);
    const auto displacement = impl.geometry->face_displacement_m(face);
    const auto area_vector =
        impl.geometry->face_area_vector_m2(face, mesh::FaceSide::owner);
    const double denominator =
        std::abs(area_vector.x * displacement.x +
                 area_vector.y * displacement.y +
                 area_vector.z * displacement.z);
    const double factor = area * area / denominator;
    if (!std::isfinite(mass_flux) || !std::isfinite(factor) || factor <= 0.0)
      throw runtime::Error("time-control committed face state is invalid");
    if (owner < cells) {
      flux_sum[owner] += std::abs(mass_flux);
      geometry_sum[owner] += factor;
    }
    const auto neighbour = impl.topology->neighbour(face);
    if (neighbour && *neighbour < cells) {
      flux_sum[*neighbour] += std::abs(mass_flux);
      geometry_sum[*neighbour] += factor;
    }
  }
  Stability result;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const double rho = variable_density ? density[cell] : density_constant;
    const double volume = impl.geometry->cell_volume_m3(cell);
    if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(volume) ||
        volume <= 0.0)
      throw runtime::Error("time-control committed cell state is invalid");
    result.convection =
        std::max(result.convection, flux_sum[cell] / (2.0 * rho * volume));
    result.diffusion =
        std::max(result.diffusion,
                 gamma * geometry_sum[cell] / (2.0 * rho * volume));
  }
  return result;
}

} // namespace

struct detail::AdaptiveTimeControlEngine final {
template <class Report>
static void store_final(TimeAdvanceReport &destination,
                        Report &&source) noexcept {
  using Concrete = std::decay_t<Report>;
  std::get<Concrete>(destination.final_attempt_) =
      std::forward<Report>(source);
  destination.final_attempt_available_ = true;
}

template <class Facade, class Report, class Invoke>
static TimeAdvanceReport run(Bdf2RetryController &controller, FlowState &state,
                              Facade &facade, double density_constant,
                              bool variable_density, double mu,
                              const linear::SolveControl &momentum,
                              const linear::SolveControl &pressure,
                              TimeAdvanceReport report, Invoke &&invoke) {
  auto &impl = *controller.impl_;
  struct ActiveGuard final {
    bool &active;
    ~ActiveGuard() { active = false; }
  } guard{impl.active};
  report.controller_identity_ = impl.identity;
  report.report_identity_ = impl.report_identity + 1U;
  report.proposed_next_dt_s_ = controller.observer_state_.proposed_next_dt_s;

  bool local_ok = detail::AdaptiveTimeControlAccess::state_identity(state) ==
                      impl.state &&
                  controls_valid(momentum) && controls_valid(pressure) &&
                  std::isfinite(mu) && mu >= 0.0;
  if (!variable_density)
    local_ok = local_ok && std::isfinite(density_constant) &&
               density_constant > 0.0;
  auto status =
      runtime::collective_status(*impl.mpi, local_ok,
                                 "time-control.preflight.config");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ = TimeAdvanceReport::PreflightCategory::config;
    impl.report_identity = report.report_identity_;
    return report;
  }
  const bool state_agrees =
      agrees_u64(*impl.mpi, controller.observer_state_.state_seal,
                 "MPI_Bcast(time-control state agreement)");
  status = runtime::collective_status(
      *impl.mpi, state_agrees, "time-control.preflight.state");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ = TimeAdvanceReport::PreflightCategory::state;
    impl.report_identity = report.report_identity_;
    return report;
  }
  if (controller.observer_state_.revision ==
          std::numeric_limits<std::uint64_t>::max() ||
      state.metadata().step == std::numeric_limits<std::uint64_t>::max()) {
    report.preflight_category_ = TimeAdvanceReport::PreflightCategory::state;
    impl.report_identity = report.report_identity_;
    return report;
  }
  const auto &authority = detail::AdaptiveTimeControlAccess::authority(facade);
  const bool authority_ok =
      std::isfinite(authority.maximum) && authority.maximum >= 0.0;
  status = runtime::collective_status(
      *impl.mpi, authority_ok,
      "time-control.preflight.transport-authority");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ =
        TimeAdvanceReport::PreflightCategory::transport_authority;
    impl.report_identity = report.report_identity_;
    return report;
  }
  const bool authority_agrees =
      agrees_u64(*impl.mpi, authority.count,
                 "MPI_Bcast(time-control transport count)") &&
      agrees_u64(*impl.mpi, authority.ordered_fingerprint,
                 "MPI_Bcast(time-control transport fingerprint)") &&
      agrees_u64(*impl.mpi, bits(authority.maximum),
                 "MPI_Bcast(time-control transport maximum)") &&
      agrees_u64(*impl.mpi, bits(mu), "MPI_Bcast(time-control viscosity)") &&
      (variable_density ||
       agrees_u64(*impl.mpi, bits(density_constant),
                  "MPI_Bcast(time-control density)"));
  status = runtime::collective_status(
      *impl.mpi, authority_agrees,
      "time-control.preflight.transport-authority");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ =
        TimeAdvanceReport::PreflightCategory::transport_authority;
    impl.report_identity = report.report_identity_;
    return report;
  }

  Stability rates{};
  if (impl.config.mode == config::TimeMode::adaptive) {
    bool rates_ok = true;
    try {
      rates = stability_rates(impl, state, density_constant, variable_density,
                              std::max(mu, authority.maximum));
    } catch (...) {
      rates_ok = false;
    }
    status = runtime::collective_status(
        *impl.mpi, rates_ok, "time-control.preflight.state");
    if (!status.ok) {
      report.lowest_failing_rank_ = status.failing_rank;
      report.preflight_category_ = TimeAdvanceReport::PreflightCategory::state;
      impl.report_identity = report.report_identity_;
      return report;
    }
    double pair[]{rates.convection, rates.diffusion};
    impl.mpi->allreduce_fp64_in_place(
        pair, 2U, runtime::Fp64ReductionOperation::maximum);
    rates = {pair[0], pair[1]};
    report.stability_metrics_available_ = true;
    report.convective_rate_per_s_ = rates.convection;
    report.diffusive_rate_per_s_ = rates.diffusion;
  }

  double dt = impl.config.mode == config::TimeMode::fixed
                  ? impl.config.initial_dt_s
                  : std::clamp(controller.observer_state_.proposed_next_dt_s,
                               impl.config.min_dt_s, impl.config.max_dt_s);
  if (impl.config.mode == config::TimeMode::adaptive) {
    double candidate = dt;
    if (rates.convection > 0.0)
      candidate =
          std::min(candidate, impl.config.cfl_target / rates.convection);
    if (rates.diffusion > 0.0)
      candidate = std::min(candidate,
                           impl.config.diffusion_number_target / rates.diffusion);
    report.limited_by_min_dt_ = candidate <= impl.config.min_dt_s;
    dt = std::max(candidate, impl.config.min_dt_s);
  }

  std::uint32_t retry_count{};
  while (true) {
    const auto metadata = state.metadata();
    const double ratio = metadata.step == 0U ? 0.0 : dt / metadata.dt_s;
    const auto order =
        metadata.step != 0U && ratio >= 0.5 && ratio <= 2.0
            ? MomentumTimeOrder::bdf2
            : MomentumTimeOrder::backward_euler;
    const auto stencil =
        make_momentum_time_stencil(order, dt, metadata.dt_s);
    Report attempt = invoke(stencil);
    const auto &base = base_report(attempt);
    const bool work_gate = half_work_gate(base, momentum, pressure);
    report.attempts_[report.attempt_count_++] =
        {dt, order, base.disposition, base.reason, base.lowest_failing_rank,
         work_gate};
    report.reason_ = base.reason;
    report.lowest_failing_rank_ = base.lowest_failing_rank;
    store_final(report, std::move(attempt));

    if (base.disposition == StepAttemptDisposition::committed) {
      TimeControlState next = controller.observer_state_;
      next.accepted_step = state.metadata().step;
      next.last_accepted_dt_s = dt;
      next.last_accepted_order = order;
      next.history_ready = true;
      next.last_all_linear_solves_within_half_limit = work_gate;
      next.last_convective_rate_per_s =
          impl.config.mode == config::TimeMode::adaptive ? rates.convection
                                                        : 0.0;
      next.last_diffusive_rate_per_s =
          impl.config.mode == config::TimeMode::adaptive ? rates.diffusion
                                                        : 0.0;
      next.last_stability_metrics_available =
          impl.config.mode == config::TimeMode::adaptive;
      next.last_retry_count = retry_count;
      ++next.revision;
      next.proposed_next_dt_s =
          impl.config.mode == config::TimeMode::fixed
              ? impl.config.initial_dt_s
              : expected_proposal(impl.config, next);
      next.state_seal =
          detail::TimeControlStateCodec::seal(impl.config, impl.model, next);
      controller.observer_state_ = next;
      report.disposition_ = TimeAdvanceDisposition::committed;
      report.reason_ = StepFailureReason::none;
      report.lowest_failing_rank_ = -1;
      report.accepted_dt_s_ = dt;
      report.proposed_next_dt_s_ = next.proposed_next_dt_s;
      impl.report_identity = report.report_identity_;
      return report;
    }
    if (base.disposition != StepAttemptDisposition::recoverable_failure) {
      report.disposition_ = TimeAdvanceDisposition::non_retryable_failure;
      impl.report_identity = report.report_identity_;
      return report;
    }
    if (bits(dt) == bits(impl.config.min_dt_s)) {
      report.disposition_ = TimeAdvanceDisposition::minimum_dt_failure;
      report.limited_by_min_dt_ = true;
      impl.report_identity = report.report_identity_;
      return report;
    }
    if (retry_count == 8U) {
      report.disposition_ = TimeAdvanceDisposition::retry_limit_reached;
      impl.report_identity = report.report_identity_;
      return report;
    }
    ++retry_count;
    const double candidate = impl.config.retry_factor * dt;
    report.limited_by_min_dt_ = candidate <= impl.config.min_dt_s;
    dt = std::max(impl.config.min_dt_s, candidate);
  }
}
};

Bdf2RetryController Bdf2RetryController::create(
    const config::FlowTimeConfig &config, config::DensityModel model,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const runtime::MpiContext &mpi, const FlowState &state) {
  bool local_ok = valid_config(config);
  collective_factory_check(mpi, local_ok,
                           "time-control.create.invalid-config");
  local_ok = geometry.compatible(topology);
  collective_factory_check(mpi, local_ok, "time-control.create.layout");
  auto observer = fresh_state(config, model);
  local_ok = detail::TimeControlStateCodec::semantically_valid(
      config, model, state.metadata(), observer);
  collective_factory_check(mpi, local_ok, "time-control.create.state");
  const auto metadata = state.metadata();
  const bool agreement =
      agrees_u64(mpi, observer.state_seal,
                 "MPI_Bcast(time-control create state)") &&
      agrees_u64(mpi, metadata.step,
                 "MPI_Bcast(time-control create step)") &&
      agrees_u64(mpi, bits(metadata.time_s),
                 "MPI_Bcast(time-control create time)");
  collective_factory_check(mpi, agreement, "time-control.create.agreement");
  auto impl = std::make_unique<Impl>();
  impl->config = config;
  impl->model = model;
  impl->topology = &topology;
  impl->geometry = &geometry;
  impl->mpi = &mpi;
  impl->state = &state;
  return Bdf2RetryController(std::move(impl), observer);
}

Bdf2RetryController Bdf2RetryController::restore(
    const config::FlowTimeConfig &config, config::DensityModel model,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const runtime::MpiContext &mpi, const FlowState &state,
    const TimeControlState &observer) {
  collective_factory_check(mpi, valid_config(config),
                           "time-control.restore.invalid-config");
  collective_factory_check(mpi, geometry.compatible(topology),
                           "time-control.restore.layout");
  collective_factory_check(
      mpi, detail::TimeControlStateCodec::semantically_valid(
               config, model, state.metadata(), observer),
      "time-control.restore.state");
  const auto metadata = state.metadata();
  const bool agreement =
      agrees_u64(mpi, observer.state_seal,
                 "MPI_Bcast(time-control restore state)") &&
      agrees_u64(mpi, metadata.step,
                 "MPI_Bcast(time-control restore step)") &&
      agrees_u64(mpi, bits(metadata.time_s),
                 "MPI_Bcast(time-control restore time)");
  collective_factory_check(mpi, agreement, "time-control.restore.agreement");
  auto impl = std::make_unique<Impl>();
  impl->config = config;
  impl->model = model;
  impl->topology = &topology;
  impl->geometry = &geometry;
  impl->mpi = &mpi;
  impl->state = &state;
  return Bdf2RetryController(std::move(impl), observer);
}

Bdf2RetryController::Bdf2RetryController(
    std::unique_ptr<Impl> impl, TimeControlState observer_state) noexcept
    : impl_(std::move(impl)), observer_state_(observer_state) {}
Bdf2RetryController::~Bdf2RetryController() noexcept = default;
Bdf2RetryController::Bdf2RetryController(Bdf2RetryController &&other) noexcept
    : impl_(std::move(other.impl_)), observer_state_(other.observer_state_) {
  other.observer_state_ = {};
}
TimeControlState Bdf2RetryController::state() const noexcept {
  return observer_state_;
}

TimeAdvanceReport Bdf2RetryController::advance(
    FlowState &state, FixedStepConstantDensityFlow &facade, double rho_ref,
    double mu, const linear::SolveControl &momentum,
    const linear::SolveControl &pressure) {
  if (!impl_ || impl_->active)
    throw runtime::Error("time-control controller is moved-from or active");
  impl_->active = true;
  TimeAdvanceReport result(TimeAdvanceReport::ConstantReportTag{});
  return detail::AdaptiveTimeControlEngine::run<
      FixedStepConstantDensityFlow, StepAttemptReport>(
      *this, state, facade, rho_ref, false, mu, momentum, pressure,
      std::move(result), [&](const auto &stencil) {
        return facade.attempt(state, rho_ref, mu, stencil, momentum, pressure);
      });
}
TimeAdvanceReport Bdf2RetryController::advance(
    FlowState &state, FixedStepMaterialDensityFlow &facade, double mu,
    const linear::SolveControl &momentum,
    const linear::SolveControl &pressure) {
  if (!impl_ || impl_->active)
    throw runtime::Error("time-control controller is moved-from or active");
  impl_->active = true;
  TimeAdvanceReport result(TimeAdvanceReport::MaterialReportTag{});
  return detail::AdaptiveTimeControlEngine::run<
      FixedStepMaterialDensityFlow, MaterialDensityStepAttemptReport>(
      *this, state, facade, 0.0, true, mu, momentum, pressure,
      std::move(result), [&](const auto &stencil) {
        return facade.attempt(state, mu, stencil, momentum, pressure);
      });
}
TimeAdvanceReport Bdf2RetryController::advance(
    FlowState &state, FixedStepIdealGasFlow &facade, double mu,
    const linear::SolveControl &momentum,
    const linear::SolveControl &pressure) {
  if (!impl_ || impl_->active)
    throw runtime::Error("time-control controller is moved-from or active");
  impl_->active = true;
  TimeAdvanceReport result(TimeAdvanceReport::IdealGasReportTag{});
  return detail::AdaptiveTimeControlEngine::run<
      FixedStepIdealGasFlow, IdealGasStepAttemptReport>(
      *this, state, facade, 0.0, true, mu, momentum, pressure,
      std::move(result), [&](const auto &stencil) {
        return facade.attempt(state, mu, stencil, momentum, pressure);
      });
}

TimeControlDiagnosticSource Bdf2RetryController::diagnostic_source(
    const FlowState &state, const TimeAdvanceReport &report) const {
  if (!impl_)
    throw runtime::Error("time-control diagnostic source is stale");
  if (&state != impl_->state || report.controller_identity_ != impl_->identity ||
      report.report_identity_ == 0U ||
      report.report_identity_ != impl_->report_identity)
    throw runtime::Error("time-control diagnostic source identity mismatch");
  auto source = std::make_unique<TimeControlDiagnosticSource::Impl>();
  auto &snapshot = source->snapshot;
  snapshot.state = observer_state_;
  snapshot.attempts = report.attempts_;
  snapshot.attempt_count = report.attempt_count_;
  snapshot.disposition = report.disposition_;
  snapshot.reason = report.reason_;
  snapshot.lowest_failing_rank = report.lowest_failing_rank_;
  snapshot.accepted_dt_s = report.accepted_dt_s_;
  snapshot.proposed_next_dt_s = report.proposed_next_dt_s_;
  snapshot.convective_rate_per_s = report.convective_rate_per_s_;
  snapshot.diffusive_rate_per_s = report.diffusive_rate_per_s_;
  snapshot.stability_metrics_available = report.stability_metrics_available_;
  snapshot.limited_by_min_dt = report.limited_by_min_dt_;
  snapshot.config = impl_->config;
  snapshot.model = impl_->model;
  snapshot.controller_identity = impl_->identity;
  snapshot.report_identity = report.report_identity_;
  snapshot.flow_state_identity =
      detail::AdaptiveTimeControlAccess::diagnostic_identity(state);
  snapshot.relative_rank = impl_->mpi->rank();
  snapshot.observed_step = state.metadata().step;
  snapshot.observed_time_s = state.metadata().time_s;
  snapshot.observed_metadata = state.metadata();
  snapshot.local_box = impl_->topology->owned_global_box();
  snapshot.global_extent = impl_->topology->global_extent();
  snapshot.local_faces = impl_->topology->local_face_count();
  snapshot.global_faces = impl_->topology->global_face_count();
  for (mesh::LocalFaceId face = 0; face < snapshot.local_faces; ++face)
    if (impl_->topology->cell_ownership(impl_->topology->owner(face)) ==
        mesh::EntityOwnership::owned)
      ++snapshot.canonical_owned_faces;
  snapshot.local_cell_layout = detail::render_owned_cells(snapshot.local_box);
  snapshot.global_cell_layout =
      detail::render_global_cells(snapshot.global_extent);
  snapshot.local_face_layout =
      detail::render_owned_faces(snapshot.canonical_owned_faces);
  snapshot.global_face_layout =
      detail::render_global_faces(snapshot.global_faces);
  return TimeControlDiagnosticSource(std::move(source));
}

static_assert(
    std::is_nothrow_move_constructible_v<TimeAdvanceReport>);
static_assert(
    std::is_nothrow_move_constructible_v<TimeControlDiagnosticSource>);
static_assert(
    std::is_nothrow_move_constructible_v<Bdf2RetryController>);
static_assert(
    std::is_nothrow_move_assignable_v<MaterialDensityStepAttemptReport>);
static_assert(std::is_nothrow_move_assignable_v<IdealGasStepAttemptReport>);

} // namespace hundun::flow
