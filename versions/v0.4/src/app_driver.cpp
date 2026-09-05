// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include <fcntl.h>
#include <sys/resource.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <sstream>
#include <string_view>
#include <vector>

#include "app_driver_detail.hpp"
#include "app_evidence_detail.hpp"
#include "app_identity_detail.hpp"
#include "hundun/v04_app.hpp"
#include "hundun/v04_case.hpp"
#include "io_output_detail.hpp"

namespace hundun::v04 {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace detail {
namespace {
thread_local ApplicationFailurePhase failure_phase =
    ApplicationFailurePhase::none;
thread_local int failure_rank = -1;
}  // namespace
void arm_application_local_allocation_failure_for_test(
    ApplicationFailurePhase phase, int rank) noexcept {
  failure_phase = phase;
  failure_rank = rank;
}
}  // namespace detail
#endif
namespace {

namespace fs = std::filesystem;

constexpr std::uint32_t kApplicationInput = 10501U;
constexpr std::uint32_t kApplicationPath = 10502U;
constexpr std::uint32_t kApplicationTemplate = 10503U;

void local_allocation_checkpoint(ApplicationFailurePhase phase, int rank) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (detail::failure_phase == phase && detail::failure_rank == rank) {
    detail::failure_phase = ApplicationFailurePhase::none;
    throw std::bad_alloc{};
  }
#else
  (void)phase;
  (void)rank;
#endif
}

static_assert(kRuntimePressureEnergyRefinementCapacity ==
              kPressureEnergyRefinementCapacity);

RuntimePressureEnergyRefinementTermination
runtime_pressure_energy_refinement_termination(
    PressureEnergyRefinementTermination termination) noexcept {
  switch (termination) {
    case PressureEnergyRefinementTermination::component_residuals_converged:
      return RuntimePressureEnergyRefinementTermination::
          component_residuals_converged;
    case PressureEnergyRefinementTermination::iteration_capacity_exhausted:
      return RuntimePressureEnergyRefinementTermination::
          iteration_capacity_exhausted;
    case PressureEnergyRefinementTermination::rejected_candidate:
      return RuntimePressureEnergyRefinementTermination::rejected_candidate;
    case PressureEnergyRefinementTermination::none:
      return RuntimePressureEnergyRefinementTermination::none;
  }
  return RuntimePressureEnergyRefinementTermination::none;
}

Status collect_evidence_resources(
    MPI_Comm communicator, const DriverResourceReport& local,
    DriverResourceReport& global) noexcept {
  // Job-wide traffic is a sum.  Logical solver work and elapsed reduction
  // time are reported at the slowest rank so every rank emits one identical
  // committed-step record without hiding load imbalance.
  std::array<std::uint64_t, 16U> maxima{{
      local.structured_messages,
      local.structured_bytes,
      local.ibm_messages,
      local.ibm_bytes,
      local.reduction_collectives,
      local.reduction_nanoseconds,
      local.linear_iterations,
      local.exact_numeric_refills,
      local.hierarchy_rebuilds,
      local.preconditioner_applications,
      local.structured_exchanges,
      local.ibm_exchanges,
      local.reduction_logical_bytes,
      local.mg_blocking_collectives,
      local.mg_collective_logical_bytes,
      local.predictor_blocking_collectives,
  }};
  if (MPI_Allreduce(MPI_IN_PLACE, maxima.data(),
                    static_cast<int>(maxima.size()), MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kApplicationInput};
  int ranks = 0;
  if (MPI_Comm_size(communicator, &ranks) != MPI_SUCCESS || ranks <= 0)
    return {StatusCode::mpi_failure, kApplicationInput};
  for (std::size_t index = 0U; index < 4U; ++index)
    if (maxima[index] >
        UINT64_MAX / static_cast<std::uint64_t>(ranks))
      return {StatusCode::invalid_plan, kApplicationInput};
  std::array<std::uint64_t, 4U> sums{{
      local.structured_messages,
      local.structured_bytes,
      local.ibm_messages,
      local.ibm_bytes,
  }};
  if (MPI_Allreduce(MPI_IN_PLACE, sums.data(), static_cast<int>(sums.size()),
                    MPI_UINT64_T, MPI_SUM, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kApplicationInput};
  global = local;
  global.structured_messages = sums[0U];
  global.structured_bytes = sums[1U];
  global.ibm_messages = sums[2U];
  global.ibm_bytes = sums[3U];
  global.reduction_collectives = maxima[4U];
  global.reduction_nanoseconds = maxima[5U];
  global.linear_iterations = maxima[6U];
  global.exact_numeric_refills = maxima[7U];
  global.hierarchy_rebuilds = maxima[8U];
  global.preconditioner_applications = maxima[9U];
  global.structured_exchanges = maxima[10U];
  global.ibm_exchanges = maxima[11U];
  global.reduction_logical_bytes = maxima[12U];
  global.mg_blocking_collectives = maxima[13U];
  global.mg_collective_logical_bytes = maxima[14U];
  global.predictor_blocking_collectives = maxima[15U];
  return {};
}

Status collect_peak_rss(MPI_Comm communicator, std::uint64_t& maximum_rank,
                        std::uint64_t& maximum_node) noexcept {
  rusage usage{};
  Status status;
  if (::getrusage(RUSAGE_SELF, &usage) != 0 || usage.ru_maxrss < 0)
    status = {StatusCode::io_failure, kApplicationInput};
  const auto kib = static_cast<std::uint64_t>(usage.ru_maxrss);
  if (status && kib > UINT64_MAX / 1024U)
    status = {StatusCode::invalid_plan, kApplicationInput};
  status = detail::output_collective_status(communicator, status);
  if (!status) return status;
  const std::uint64_t local = kib * 1024U;
  if (MPI_Allreduce(&local, &maximum_rank, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kApplicationInput};

  MPI_Comm shared = MPI_COMM_NULL;
  if (MPI_Comm_split_type(communicator, MPI_COMM_TYPE_SHARED, 0,
                          MPI_INFO_NULL, &shared) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kApplicationInput};
  std::uint64_t node = 0U;
  const int node_sum = MPI_Allreduce(&local, &node, 1, MPI_UINT64_T, MPI_SUM,
                                     shared);
  const int free_status = MPI_Comm_free(&shared);
  if (node_sum != MPI_SUCCESS || free_status != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kApplicationInput};
  if (MPI_Allreduce(&node, &maximum_node, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kApplicationInput};
  return {};
}

Status collect_stage_timings(
    MPI_Comm communicator, const DriverStepReport& step,
    std::array<StageTimingRecord, kDriverTimedStageCapacity>& stages) noexcept {
  constexpr std::array<StageId, kDriverTimedStageCapacity> kStages{
      10U, 12U, 15U, 20U, 30U, 40U, 45U, 50U, 60U, 70U};
  std::array<std::uint64_t, kDriverTimedStageCapacity> local{};
  for (std::size_t index = 0U; index < step.stage_timing_count; ++index) {
    const DriverStageTiming timing = step.stage_timings[index];
    const auto found = std::find(kStages.begin(), kStages.end(), timing.stage);
    if (timing.stage == 0U || found == kStages.end())
      return {StatusCode::invalid_plan, kApplicationInput};
    const std::size_t target =
        static_cast<std::size_t>(found - kStages.begin());
    if (local[target] != 0U)
      return {StatusCode::invalid_plan, kApplicationInput};
    local[target] = timing.nanoseconds;
  }
  auto minimum = local;
  auto maximum = local;
  auto sum = local;
  if (MPI_Allreduce(MPI_IN_PLACE, minimum.data(),
                    static_cast<int>(minimum.size()), MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, maximum.data(),
                    static_cast<int>(maximum.size()), MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kApplicationInput};
  int ranks = 0;
  if (MPI_Comm_size(communicator, &ranks) != MPI_SUCCESS || ranks <= 0)
    return {StatusCode::mpi_failure, kApplicationInput};
  for (std::uint64_t value : maximum)
    if (value > UINT64_MAX / static_cast<std::uint64_t>(ranks))
      return {StatusCode::invalid_plan, kApplicationInput};
  if (MPI_Allreduce(MPI_IN_PLACE, sum.data(), static_cast<int>(sum.size()),
                    MPI_UINT64_T, MPI_SUM, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kApplicationInput};
  for (std::size_t index = 0U; index < stages.size(); ++index) {
    stages[index] = {kStages[index], minimum[index],
                     sum[index] / static_cast<std::uint64_t>(ranks),
                     maximum[index]};
  }
  return {};
}

constexpr std::string_view kCaseJson = R"json({
  "schema_version": 1,
  "units": "SI",
  "mesh": {
    "kind": "uniform",
    "domain": {"lower": [0, 0, 0], "upper": [1, 1, 1]},
    "exact_cells": [8, 8, 8],
    "base_spacing": null,
    "minimum_spacing": [0.125, 0.125, 0.125],
    "max_growth_ratio": 1.0,
    "focus_regions": [],
    "limits": {"max_global_cells": 512,
               "max_memory_bytes_per_rank": 268435456},
    "data_files": [],
    "immersed_boundary": null
  },
  "flow": {"model": "single_phase_low_mach_compressible",
           "pressure_reference": "closed_mass", "reacting": false},
  "solver": {
    "coupling": "PISO",
    "pressure_correctors": 2,
    "pressure_linear": {
      "absolute_tolerance": 1e-13,
      "relative_tolerance": 1e-13,
      "maximum_iterations": 400,
      "true_residual_interval": 4,
      "krylov_restart": 12
    },
    "terminal_tolerances": {
      "eos": 1e-10,
      "continuity": 1e-10,
      "closed_mass": 1e-10,
      "gauge": 1e-10
    }
  },
  "turbulence": {"model": "vreman_wall_function"},
  "thermophysics": {"data_file": "thermophysics.d"},
  "transported_scalars": [],
  "boundaries": {
    "x_min": {"flow_kind":"periodic","thermal_kind":"none","velocity":[0,0,0],"direction":[1,0,0],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]},
    "x_max": {"flow_kind":"periodic","thermal_kind":"none","velocity":[0,0,0],"direction":[1,0,0],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]},
    "y_min": {"flow_kind":"periodic","thermal_kind":"none","velocity":[0,0,0],"direction":[0,1,0],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]},
    "y_max": {"flow_kind":"periodic","thermal_kind":"none","velocity":[0,0,0],"direction":[0,1,0],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]},
    "z_min": {"flow_kind":"periodic","thermal_kind":"none","velocity":[0,0,0],"direction":[0,0,1],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]},
    "z_max": {"flow_kind":"periodic","thermal_kind":"none","velocity":[0,0,0],"direction":[0,0,1],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]}
  },
  "schemes": {"momentum":"limited_central2","enthalpy":"limited_central2","species":"tvd2","passive_scalar":"tvd2","diffusion":"central2","limiter":1.0},
  "time": {"control":"adaptive_flow","scheme":"variable_bdf2","initial_dt":0.001,"minimum_dt":1e-8,"maximum_dt":0.1,"convective_cfl":0.8,"viscous_cfl":0.5,"thermal_cfl":0.5,"species_cfl":0.5,"acoustic_cfl":0.8,"maximum_growth":1.2,"retry_factor":0.5,"maximum_retries":6,"minimum_bdf_ratio":0.25,"maximum_bdf_ratio":4.0}
})json";

constexpr std::string_view kThermophysics = R"data(HUNDUN_THERMOPHYSICS_V1
temperature_bounds 200 2000
temperature_inversion 1e-12 64
closed_mass_newton 1e-12 32 0.2
species_count 1
species air
molecular_weight 28.96546
temperature_switch 1000
nasa7_low 3.5 0 0 0 0 0 0
nasa7_high 3.5 0 0 0 0 0 0
transport_constant 1.8e-5 0.026
end_species
end
)data";

bool write_exclusive(const fs::path& path, std::string_view bytes) noexcept {
  const int descriptor =
      ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0644);
  if (descriptor < 0) return false;
  std::size_t cursor = 0U;
  bool okay = true;
  while (cursor < bytes.size()) {
    const ssize_t count =
        ::write(descriptor, bytes.data() + cursor, bytes.size() - cursor);
    if (count <= 0) {
      okay = false;
      break;
    }
    cursor += static_cast<std::size_t>(count);
  }
  if (okay && ::fsync(descriptor) != 0) okay = false;
  if (::close(descriptor) != 0) okay = false;
  return okay;
}

bool sync_directory(const fs::path& path) noexcept {
  const int descriptor = ::open(path.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
  if (descriptor < 0) return false;
  const bool okay = ::fsync(descriptor) == 0;
  return ::close(descriptor) == 0 && okay;
}

bool within(const fs::path& candidate, const fs::path& parent) noexcept {
  // Both operands have already been made absolute and canonical.  Comparing
  // native strings at a separator boundary avoids allocating path components
  // (some libc++ path iterators terminate on component-allocation failure).
  const auto& value = candidate.native();
  const auto& base = parent.native();
  std::size_t length = base.size();
  while (length > 1U && base[length - 1U] == fs::path::preferred_separator)
    --length;
  return length != 0U && value.size() >= length &&
         value.compare(0U, length, base, 0U, length) == 0 &&
         (value.size() == length ||
          base[length - 1U] == fs::path::preferred_separator ||
          value[length] == fs::path::preferred_separator);
}

bool absolute_normalized(const fs::path& path, fs::path& out) {
  std::error_code error;
  out = fs::weakly_canonical(path, error);
  if (!error) return true;
  error.clear();
  out = fs::absolute(path, error).lexically_normal();
  return !error;
}

}  // namespace

Status ApplicationService::validate(MPI_Comm communicator,
                                    const std::filesystem::path& case_root,
                                    CaseValidationReport& out) noexcept try {
  if (communicator == MPI_COMM_NULL || case_root.empty())
    return {StatusCode::invalid_case, kApplicationInput};
  ValidatedModel model;
  Status status = CaseCompiler::load_and_compile(communicator, case_root, model);
  if (!status) return status;
  CompiledCasePlan product;
  status = ProductCompiler::compile(communicator, model, case_root, product);
  if (!status) return status;
  const CaseValidationReport candidate{model.fingerprint, product.fingerprint(),
                                       product.summary()};
  out = candidate;
  return {};
} catch (...) {
  return {StatusCode::allocation_failure, kApplicationInput};
}

Status ApplicationService::initialize_case_directory(
    const std::filesystem::path& output_directory) noexcept try {
  if (output_directory.empty())
    return {StatusCode::invalid_case, kApplicationTemplate};
  std::error_code error;
  if (fs::exists(output_directory, error) || error)
    return {StatusCode::invalid_case, kApplicationTemplate};
  if (!fs::create_directories(output_directory, error) || error)
    return {StatusCode::io_failure, kApplicationTemplate};
  if (!write_exclusive(output_directory / "case.json", kCaseJson) ||
      !write_exclusive(output_directory / "thermophysics.d",
                       kThermophysics) ||
      !sync_directory(output_directory)) {
    return {StatusCode::io_failure, kApplicationTemplate};
  }
  return {};
} catch (...) {
  return {StatusCode::io_failure, kApplicationTemplate};
}

Status ApplicationService::validate_run_directories(
    const std::filesystem::path& case_root,
    const std::filesystem::path& run_directory,
    const std::filesystem::path& source_root) noexcept try {
  fs::path canonical_case;
  fs::path canonical_run;
  fs::path canonical_source;
  if (case_root.empty() || run_directory.empty() || source_root.empty() ||
      !absolute_normalized(case_root, canonical_case) ||
      !absolute_normalized(run_directory, canonical_run) ||
      !absolute_normalized(source_root, canonical_source) ||
      canonical_case == canonical_run || within(canonical_run, canonical_case) ||
      within(canonical_run, canonical_source)) {
    return {StatusCode::invalid_case, kApplicationPath};
  }
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kApplicationPath};
} catch (...) {
  return {StatusCode::invalid_case, kApplicationPath};
}

Status ApplicationService::run(MPI_Comm communicator,
                               const ApplicationRunOptions& options,
                               ApplicationRunReport& report) noexcept try {
  if (communicator == MPI_COMM_NULL || options.case_root.empty() ||
      options.run_directory.empty() || options.source_root.empty() ||
      options.steps == 0U)
    return {StatusCode::invalid_case, kApplicationInput};
  Status status = detail::output_collective_status(
      communicator,
      validate_run_directories(options.case_root, options.run_directory,
                               options.source_root));
  ValidatedModel model;
  if (status)
    status = CaseCompiler::load_and_compile(communicator, options.case_root,
                                            model);
  CompiledCasePlan plan;
  if (status)
    status = ProductCompiler::compile(communicator, model, options.case_root,
                                      plan);
  if (!status) return status;
  const PlanFingerprint product_fingerprint = plan.fingerprint();
  const PlanFingerprint cpu_plan_fingerprint = plan.cpu_plan_fingerprint();
  const PlanFingerprint stl_fingerprint = plan.stl_fingerprint();
  const IoServicePlan* sealed_services = plan.io_services();
  if (sealed_services == nullptr || sealed_services->fingerprint() == 0U)
    return {StatusCode::invalid_plan, kApplicationInput};
  const IoServicePlan& services = *sealed_services;
  ProductDriver driver;
  status = ProductDriver::create(communicator, std::move(plan), driver);
  std::uint64_t starting_step = 0U;
  RuntimeRunStartAnchor run_start;
  bool restart_backward_euler_recovery = false;
  if (status && !options.restart_directory.empty()) {
    RestartExpected expected;
    status = driver.restart_expected(expected);
    RestartImage image;
    if (status)
      status = RestartReader::load(communicator, options.restart_directory,
                                   expected, image);
    if (status) {
      starting_step = image.step;
      run_start.kind = RuntimeRunStartKind::restart;
      run_start.previous_step = image.step;
      run_start.previous_time = image.time;
      run_start.restart_manifest_sha256 = image.source_manifest_sha256;
      restart_backward_euler_recovery = image.backward_euler_recovery;
      status = driver.initialize_restart(image);
    }
  } else if (status) {
    std::vector<double> initial_scalars(model.transported_scalars.size(), 0.0);
    DriverInitialState initial;
    initial.transported_scalars =
        {initial_scalars.data(), initial_scalars.size()};
    for (const BoundaryFaceSpec& boundary : model.boundaries) {
      if (std::isfinite(boundary.temperature) && boundary.temperature > 0.0)
        initial.temperature = boundary.temperature;
      if (model.pressure_reference ==
              PressureReferenceKind::boundary_absolute &&
          boundary.flow_kind == BoundaryKind::pressure_outlet &&
          std::isfinite(boundary.pressure) && boundary.pressure > 0.0)
        initial.pressure_reference = boundary.pressure;
      if (boundary.flow_kind == BoundaryKind::velocity_inlet ||
          boundary.flow_kind == BoundaryKind::static_state_inlet ||
          boundary.flow_kind == BoundaryKind::total_state_inlet)
        initial.velocity = boundary.velocity;
    }
    status = driver.initialize(initial);
  }
  if (!status) return status;
  RuntimeCandidateIdentity candidate_identity;
  status = detail::runtime_candidate_identity(communicator,
                                              candidate_identity);
  if (!status) return status;
  const PlanFingerprint build_identity =
      detail::runtime_sha256_fingerprint(
          candidate_identity.build_manifest);
  const PlanFingerprint binary_identity =
      detail::runtime_sha256_fingerprint(candidate_identity.executable);
  if (build_identity == 0U || binary_identity == 0U)
    return {StatusCode::invalid_plan, kApplicationInput};
  if (options.steps > UINT64_MAX - starting_step)
    return {StatusCode::invalid_case, kApplicationInput};
  const std::uint64_t target_step = starting_step + options.steps;

  int rank = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kApplicationInput};
  report = {};
  report.case_model = model.fingerprint;
  report.product = product_fingerprint;
  report.accepted_steps = starting_step;
  report.final_time = run_start.previous_time;
  for (std::uint64_t step_index = 0U;
       step_index < options.steps && status; ++step_index) {
    report.failure_phase = ApplicationFailurePhase::advance;
    const auto begin = std::chrono::steady_clock::now();
    DriverStepReport step;
    LocalTimeLimits time_limits = options.time_limits;
    status = driver.constrain_convective_time_limit(time_limits);
    if (status) status = driver.advance(time_limits, step);
    report.attempts = step.attempts;
    // Commit is authoritative even if a later resource collective or writer
    // fails.  Never report an accepted step as a failed numerical attempt.
    if (step.accepted) {
      report.accepted_steps = step.accepted_step;
      report.final_time = step.accepted_time;
    }
    const auto end = std::chrono::steady_clock::now();
    const auto local_nanoseconds = static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin)
            .count());
    std::uint64_t maximum_nanoseconds = 0U;
    if (MPI_Allreduce(&local_nanoseconds, &maximum_nanoseconds, 1,
                      MPI_UINT64_T, MPI_MAX, communicator) != MPI_SUCCESS)
      status = {StatusCode::mpi_failure, kApplicationInput};
    if (!status) {
      report.failed_stage = step.failed_stage;
      report.attempts = step.attempts;
      report.failure = status;
      report.requested_bdf = step.proposal.bdf;
      report.effective_bdf = step.effective_bdf;
      report.thermophysical_predictor_calls =
          step.thermophysical_predictor_calls;
      report.temporal_method_fallback = step.temporal_method_fallback;
      report.piso = step.piso;
      report.momentum_predictor_limiter =
          step.momentum_predictor_limiter;
      report.momentum_predictor_solve = step.momentum_predictor_solve;
      report.numerical_failure = step.numerical_failure;
      report.thermophysical_predictor = step.thermophysical_predictor;
      report.pressure_energy_globalization =
          step.pressure_energy_globalization;
      break;
    }
    report.requested_bdf = step.proposal.bdf;
    report.effective_bdf = step.effective_bdf;
    report.thermophysical_predictor_calls =
        step.thermophysical_predictor_calls;
    report.temporal_method_fallback = step.temporal_method_fallback;
    report.piso = step.piso;
    report.momentum_predictor_limiter = step.momentum_predictor_limiter;
    report.momentum_predictor_solve = step.momentum_predictor_solve;
    report.thermophysical_predictor = step.thermophysical_predictor;
    report.pressure_energy_globalization =
        step.pressure_energy_globalization;
    report.maximum_advective_convective_cfl_out = std::max(
        report.maximum_advective_convective_cfl_out,
        step.momentum_predictor_limiter.advective_cfl.out_max);
    report.maximum_advective_convective_cfl_abs = std::max(
        report.maximum_advective_convective_cfl_abs,
        step.momentum_predictor_limiter.advective_cfl.absolute_max);
    report.advective_convective_cfl_limit =
        step.momentum_predictor_limiter.advective_cfl.limit;
    report.maximum_committed_convective_cfl_out = std::max(
        report.maximum_committed_convective_cfl_out,
        step.piso.committed_convective_cfl_out_max);
    report.maximum_committed_convective_cfl_abs = std::max(
        report.maximum_committed_convective_cfl_abs,
        step.piso.committed_convective_cfl_abs_max);
    report.committed_convective_cfl_limit =
        step.piso.committed_convective_cfl_limit;
    report.minimum_predictor_theta = std::min(
        report.minimum_predictor_theta,
        step.thermophysical_predictor.theta);
    if (step.thermophysical_predictor.limited)
      ++report.predictor_limiter_activations;
    report.predictor_low_order_transport_passes +=
        step.thermophysical_predictor.low_order_transport_passes;
    report.predictor_low_order_halo_exchanges +=
        step.thermophysical_predictor.low_order_halo_exchanges;

    const bool output =
        options.output_interval != 0U &&
        (step.accepted_step % options.output_interval == 0U ||
         step.accepted_step == target_step);
    const bool restart =
        options.restart_interval != 0U &&
        (step.accepted_step % options.restart_interval == 0U ||
         step.accepted_step == target_step);
    CommittedOutputSnapshot snapshot;
    fs::path output_path;
    if (output) {
      report.failure_phase = ApplicationFailurePhase::visit;
      status = detail::output_collective_stage(communicator, [&] {
        local_allocation_checkpoint(report.failure_phase, rank);
        output_path = options.run_directory / "Visit";
        return driver.committed_output_snapshot(snapshot);
      });
    }
    if (status && output)
      status =
          VisitWriter::write(communicator, output_path, services, snapshot);
    if (status && output) {
      report.failure_phase = ApplicationFailurePhase::screen;
      std::string screen_text;
      status = detail::output_collective_stage(communicator, [&] {
        local_allocation_checkpoint(report.failure_phase, rank);
        output_path = options.run_directory / "screen.log";
        std::ostringstream summary;
        summary << "attempts=" << step.attempts
                << " continuity=" << step.piso.continuity_residual
                << " energy=" << step.piso.energy_residual
                << " eos=" << step.piso.eos_residual << " advective_cfl_out="
                << step.momentum_predictor_limiter.advective_cfl.out_max
                << " advective_cfl_abs="
                << step.momentum_predictor_limiter.advective_cfl.absolute_max
                << " advective_cfl_limit="
                << step.momentum_predictor_limiter.advective_cfl.limit
                << " committed_cfl_out="
                << step.piso.committed_convective_cfl_out_max
                << " committed_cfl_abs="
                << step.piso.committed_convective_cfl_abs_max
                << " committed_cfl_limit="
                << step.piso.committed_convective_cfl_limit
                << " predictor_limited="
                << (step.thermophysical_predictor.limited ? 1 : 0)
                << " predictor_theta=" << step.thermophysical_predictor.theta;
        screen_text = summary.str();
        return Status{};
      });
      if (status)
        status = ScreenWriter::append(communicator, output_path, services,
                                      snapshot, screen_text);
    }
    if (status && output) {
      report.failure_phase = ApplicationFailurePhase::monitor;
      std::string monitor_text;
      status = detail::output_collective_stage(communicator, [&] {
        local_allocation_checkpoint(report.failure_phase, rank);
        output_path = options.run_directory / "monitor.jsonl";
        std::ostringstream payload;
        payload << "{\"attempts\":" << step.attempts
                << ",\"continuity\":" << step.piso.continuity_residual
                << ",\"energy\":" << step.piso.energy_residual
                << ",\"eos\":" << step.piso.eos_residual
                << ",\"advective_convective_cfl_out\":"
                << step.momentum_predictor_limiter.advective_cfl.out_max
                << ",\"advective_convective_cfl_abs\":"
                << step.momentum_predictor_limiter.advective_cfl.absolute_max
                << ",\"advective_convective_cfl_limit\":"
                << step.momentum_predictor_limiter.advective_cfl.limit
                << ",\"committed_convective_cfl_out\":"
                << step.piso.committed_convective_cfl_out_max
                << ",\"committed_convective_cfl_abs\":"
                << step.piso.committed_convective_cfl_abs_max
                << ",\"committed_convective_cfl_limit\":"
                << step.piso.committed_convective_cfl_limit
                << ",\"predictor_limited\":"
                << (step.thermophysical_predictor.limited ? "true" : "false")
                << ",\"predictor_theta\":"
                << step.thermophysical_predictor.theta << '}';
        monitor_text = payload.str();
        return Status{};
      });
      if (status)
        status = MonitorWriter::append(communicator, output_path, services,
                                       snapshot, monitor_text);
    }
    RestartSnapshot restart_snapshot;
    if (status && restart) {
      report.failure_phase = ApplicationFailurePhase::restart;
      status = detail::output_collective_stage(communicator, [&] {
        local_allocation_checkpoint(report.failure_phase, rank);
        output_path = options.run_directory / "Restart";
        return driver.committed_restart_snapshot(restart_snapshot);
      });
    }
    if (status && restart)
      status = RestartWriter::write(communicator, output_path, restart_snapshot,
                                    {1U});
    if (status) {
      report.failure_phase = ApplicationFailurePhase::resources;
      status = detail::output_collective_stage(communicator, [&] {
        local_allocation_checkpoint(report.failure_phase, rank);
        return Status{};
      });
      if (!status) break;
      DriverResourceReport evidence_resources;
      status = collect_evidence_resources(communicator, step.resources,
                                          evidence_resources);
      if (!status) break;
      std::uint64_t maximum_rank_rss = 0U;
      std::uint64_t maximum_node_rss = 0U;
      status = collect_peak_rss(communicator, maximum_rank_rss,
                                maximum_node_rss);
      if (!status) break;
      std::array<StageTimingRecord, kDriverTimedStageCapacity>
          evidence_stages{};
      status = collect_stage_timings(communicator, step, evidence_stages);
      if (!status) break;
      RuntimeEvidenceRecord evidence;
      evidence.build = build_identity;
      evidence.binary = binary_identity;
      evidence.candidate_identity = candidate_identity;
      evidence.case_model = model.fingerprint;
      evidence.stl = stl_fingerprint;
      evidence.product = product_fingerprint;
      evidence.cpu_plan = cpu_plan_fingerprint;
      evidence.run_start = run_start;
      evidence.step = step.accepted_step;
      evidence.previous_committed_time = step.proposal.time;
      evidence.time = step.accepted_time;
      evidence.requested_bdf_order = step.proposal.bdf.order;
      evidence.bdf_order = step.effective_bdf.order;
      evidence.coupling =
          model.solver.coupling == CouplingKind::simple
              ? RuntimeCouplingKind::simple
              : RuntimeCouplingKind::piso;
      evidence.thermophysical_predictor_calls =
          step.thermophysical_predictor_calls;
      evidence.temporal_method_fallback = step.temporal_method_fallback;
      evidence.maximum_rank_step_nanoseconds = maximum_nanoseconds;
      evidence.maximum_rank_rss_bytes = maximum_rank_rss;
      evidence.maximum_node_rss_bytes = maximum_node_rss;
      evidence.structured_messages = evidence_resources.structured_messages;
      evidence.structured_bytes = evidence_resources.structured_bytes;
      evidence.ibm_messages = evidence_resources.ibm_messages;
      evidence.ibm_bytes = evidence_resources.ibm_bytes;
      // ReductionEngine and the native MG direct collectives are blocking in
      // v0.4.  Keep their time/byte classes separate, but include both call
      // classes in the committed blocking-collective authority.
      if (evidence_resources.mg_blocking_collectives >
          UINT64_MAX - evidence_resources.reduction_collectives) {
        status = {StatusCode::invalid_plan, kApplicationInput};
        break;
      }
      evidence.blocking_collectives =
          evidence_resources.reduction_collectives +
          evidence_resources.mg_blocking_collectives;
      if (evidence_resources.predictor_blocking_collectives >
          UINT64_MAX - evidence.blocking_collectives) {
        status = {StatusCode::invalid_plan, kApplicationInput};
        break;
      }
      evidence.blocking_collectives +=
          evidence_resources.predictor_blocking_collectives;
      evidence.reduction_nanoseconds =
          evidence_resources.reduction_nanoseconds;
      evidence.linear_iterations = evidence_resources.linear_iterations;
      evidence.exact_numeric_refills =
          evidence_resources.exact_numeric_refills;
      evidence.coarse_numeric_refills =
          evidence_resources.exact_numeric_refills;
      evidence.preconditioner_setups =
          evidence_resources.hierarchy_rebuilds;
      evidence.preconditioner_reuses =
          evidence_resources.exact_numeric_refills >=
                  evidence_resources.hierarchy_rebuilds
              ? evidence_resources.exact_numeric_refills -
                    evidence_resources.hierarchy_rebuilds
              : 0U;
      evidence.pressure = step.piso.pressure;
      evidence.pressure_solve_calls = step.piso.pressure_solve_calls;
      // ProductDriver accepts the complete attempt against the coupled
      // continuity--energy terminal contract.  SIMPLE may use a pressure-only
      // C1 direction internally; C2 and the exact candidate audit still close
      // the same coupled attempt.  Termination alone cannot identify this
      // contract because zero-RHS exits before any optional linear audit.
      evidence.pressure_solve_contract =
          RuntimePressureSolveContract::continuity_energy_coupled;
      evidence.pressure_energy_refinement_solve_calls =
          step.piso.pressure_energy_refinement_solve_calls;
      evidence.pressure_energy_refinement_termination =
          runtime_pressure_energy_refinement_termination(
              step.piso.pressure_energy_refinement_termination);
      for (std::uint8_t index = 0U;
           index < step.piso.pressure_energy_refinement_solve_calls &&
           index < evidence.pressure_energy_refinement.size();
           ++index) {
        const PisoPressureEnergyRefinementSolveReport& source =
            step.piso.pressure_energy_refinement[index];
        RuntimePressureEnergyRefinementSolve& destination =
            evidence.pressure_energy_refinement[index];
        destination.solve = source.solve;
        destination.target_generation = source.target_generation;
        destination.collective_lineage = source.collective_lineage;
        destination.ordinal = source.ordinal;
      }
      evidence.terminal_physical_audit.present =
          step.piso.final_flux_revision != 0U;
      evidence.terminal_physical_audit.final_flux_revision =
          step.piso.final_flux_revision;
      evidence.terminal_physical_audit.eos_residual =
          step.piso.eos_residual;
      evidence.terminal_physical_audit.eos_tolerance =
          model.solver.terminal.eos;
      evidence.terminal_physical_audit.continuity_residual =
          step.piso.continuity_residual;
      evidence.terminal_physical_audit.continuity_tolerance =
          model.solver.terminal.continuity;
      evidence.terminal_physical_audit.energy_residual =
          step.piso.energy_residual;
      // The v0.4 coupled energy terminal gate deliberately shares the
      // configured normalized continuity tolerance.
      evidence.terminal_physical_audit.energy_tolerance =
          model.solver.terminal.continuity;
      evidence.momentum_predictor_passes =
          step.momentum_predictor_solve.predictor_passes;
      evidence.terminal_physical_audit.closed_mass_residual =
          step.piso.closed_mass_residual;
      evidence.terminal_physical_audit.closed_mass_tolerance =
          model.solver.terminal.closed_mass;
      evidence.terminal_physical_audit.gauge_residual =
          step.piso.gauge_residual;
      evidence.terminal_physical_audit.gauge_tolerance =
          model.solver.terminal.gauge;
      status = detail::runtime_committed_cfl(
          communicator, step.piso.committed_convective_cfl,
          evidence.committed_convective_cfl);
      if (!status) break;
      evidence.momentum_predictor = step.momentum_predictor_solve.components;
      evidence.momentum_predictor_solve_calls =
          step.momentum_predictor_solve.solve_calls;
      evidence.predictor_enthalpy_endpoint =
          step.thermophysical_predictor.enthalpy_endpoint;
      evidence.predictor_enthalpy_endpoint_alpha =
          step.thermophysical_predictor.enthalpy_endpoint_alpha;
      evidence.predictor_bdf_endpoint_alpha =
          step.thermophysical_predictor.bdf_endpoint_alpha;
      evidence.predictor_source_endpoint_alpha =
          step.thermophysical_predictor.source_endpoint_alpha;
      evidence.predictor_enthalpy_solve_calls =
          step.thermophysical_predictor.enthalpy_solve_calls;
      evidence.momentum_predictor_theta =
          step.momentum_predictor_limiter.theta;
      evidence.momentum_predictor_activations =
          step.momentum_predictor_limiter.activations;
      evidence.momentum_correction_metrics_applicable =
          step.momentum_predictor_limiter.correction_metrics_applicable;
      evidence.momentum_minimum_face_alpha =
          step.momentum_predictor_limiter.minimum_face_alpha;
      evidence.momentum_active_correction_faces =
          step.momentum_predictor_limiter.active_correction_faces;
      evidence.momentum_limited_face_fraction =
          step.momentum_predictor_limiter.limited_face_fraction;
      evidence.momentum_predictor_limited =
          step.momentum_predictor_limiter.limited;
      status = detail::runtime_advective_cfl(
          communicator, step.momentum_predictor_limiter.advective_cfl,
          evidence.momentum_advective_cfl);
      if (!status) break;
      evidence.predictor_theta = step.thermophysical_predictor.theta;
      evidence.predictor_mass_flux_scale =
          step.thermophysical_predictor.mass_flux_scale;
      evidence.predictor_low_margin =
          step.thermophysical_predictor.low_margin;
      evidence.predictor_high_margin =
          step.thermophysical_predictor.high_margin;
      evidence.predictor_low_order_transport_passes =
          step.thermophysical_predictor.low_order_transport_passes;
      evidence.predictor_low_order_substeps =
          step.thermophysical_predictor.low_order_substeps;
      evidence.predictor_low_order_halo_exchanges =
          step.thermophysical_predictor.low_order_halo_exchanges;
      evidence.predictor_blocking_collectives =
          step.thermophysical_predictor.blocking_collectives;
      evidence.predictor_limiting_cell_x =
          step.thermophysical_predictor.limiting_cell.x;
      evidence.predictor_limiting_cell_y =
          step.thermophysical_predictor.limiting_cell.y;
      evidence.predictor_limiting_cell_z =
          step.thermophysical_predictor.limiting_cell.z;
      evidence.predictor_limiting_rank =
          step.thermophysical_predictor.limiting_rank;
      evidence.predictor_constraint = static_cast<std::uint8_t>(
          step.thermophysical_predictor.constraint);
      evidence.predictor_low_state = static_cast<std::uint8_t>(
          step.thermophysical_predictor.low_state);
      evidence.predictor_limited = step.thermophysical_predictor.limited;
      evidence.stages = {evidence_stages.data(), evidence_stages.size()};
      evidence.startup = step.accepted_step == 1U;
      evidence.retry = step.attempts != 1U;
      evidence.restart_recovery = restart_backward_euler_recovery &&
                                  step_index == 0U;
      // Task 20 owns candidate identity and eligibility.  Ordinary product
      // runs are deliberately incapable of masquerading as benchmark data.
      evidence.statistics_eligible = false;
      report.failure_phase = ApplicationFailurePhase::evidence;
      status = detail::output_collective_stage(communicator, [&] {
        local_allocation_checkpoint(report.failure_phase, rank);
        output_path = options.run_directory / "evidence.jsonl";
        return Status{};
      });
      if (status)
        status = EvidenceWriter::append(communicator, output_path, services,
                                        evidence);
    }
  }
  if (!status) {
    report.failure = status;
    return status;
  }
  CommittedOutputSnapshot final_snapshot;
  status = driver.committed_output_snapshot(final_snapshot);
  if (!status) return status;
  report.accepted_steps = final_snapshot.step;
  report.final_time = final_snapshot.time;
  report.failure_phase = ApplicationFailurePhase::none;
  (void)rank;
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kApplicationInput};
} catch (...) {
  return {StatusCode::invalid_case, kApplicationInput};
}

}  // namespace hundun::v04
