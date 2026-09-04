// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_product.hpp"
#include "hundun/v04_app.hpp"
#include "hundun/v04_initialization.hpp"

#include "common_terminal_audit.h"
#include "core_product_freeze_detail.hpp"
#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_piso_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
#include <atomic>
#endif
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kProductInput = 10201U;
constexpr std::uint32_t kProductRegistration = 10202U;
constexpr std::uint32_t kProductAnalysis = 10203U;
constexpr std::uint32_t kProductAllocation = 10204U;
constexpr std::uint32_t kProductInstantiation = 10205U;
constexpr std::uint32_t kProductCapacity = 10206U;
constexpr std::uint32_t kProductCommunication = 10207U;
constexpr std::uint32_t kProductBinding = 10208U;
constexpr std::uint32_t kProductCollective = 10209U;
constexpr std::uint32_t kProductPressureEnergy = 10210U;
constexpr std::uint32_t kProductConvectiveCfl = 10211U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr StageId kIbmGradientDonorStage = 131U;
constexpr StageId kIbmPressureCorrectionDonorStage = 150U;
constexpr StageId kIbmCandidatePressureCorrectionDonorStage = 151U;
constexpr StageId kIbmCandidateVelocityDonorStage = 153U;
constexpr StageId kIbmCandidateEnergyRateDonorStage = 154U;
constexpr StageId kCoupledStateC1Stage = 147U;
constexpr StageId kCoupledStateC2Stage = 157U;
constexpr StageId kCoupledStateTerminalStage = 167U;
constexpr StageId kCandidateCorrectionC1Stage = 148U;
constexpr StageId kCandidateStateC1Stage = 149U;
constexpr StageId kCandidateCorrectionC2Stage = 158U;
constexpr StageId kCandidateStateC2Stage = 159U;
constexpr StageId kCandidateFinalizerStateC1Stage = 152U;
constexpr StageId kCandidateFinalizerStateC2Stage = 162U;
constexpr StageId kFreshProjectionOperatorStage = 171U;
constexpr StageId kFreshProjectionCorrectionStage = 172U;
constexpr StageId kFreshVelocityAcceptedHaloStage = 173U;
constexpr StageId kFreshVelocityPreviousHaloStage = 174U;
constexpr StageId kFreshVelocityTrialHaloStage = 175U;
constexpr std::uint32_t kAuxiliaryFgmresMaximumRestart = 12U;

std::uint64_t product_double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

std::uint64_t product_ulp_distance(double left, double right) noexcept {
  if (!std::isfinite(left) || !std::isfinite(right))
    return std::numeric_limits<std::uint64_t>::max();
  const auto ordered = [](double value) noexcept {
    const std::uint64_t bits = product_double_bits(value);
    constexpr std::uint64_t sign = UINT64_C(1) << 63U;
    return (bits & sign) != 0U ? ~bits : bits | sign;
  };
  const std::uint64_t lhs = ordered(left);
  const std::uint64_t rhs = ordered(right);
  return lhs >= rhs ? lhs - rhs : rhs - lhs;
}

Status convective_cfl_acceptance_status(TimeControlKind control,
                                        double outward_max,
                                        double limit) noexcept {
  if (!std::isfinite(outward_max) || outward_max < 0.0 ||
      !std::isfinite(limit) || !(limit > 0.0) ||
      (control != TimeControlKind::fixed &&
       control != TimeControlKind::adaptive_flow &&
       control != TimeControlKind::adaptive_acoustic))
    return {StatusCode::invalid_plan, kProductConvectiveCfl};
  constexpr double kCflRoundoffSlack =
      64.0 * std::numeric_limits<double>::epsilon();
  if (outward_max <= limit * (1.0 + kCflRoundoffSlack)) return {};
  return control == TimeControlKind::fixed
             ? Status{StatusCode::invalid_case, kProductConvectiveCfl}
             : Status{StatusCode::rejected_step, kProductConvectiveCfl};
}

Status refresh_coast_native_air_effective_thermal_transport(
    const TransportPlan& transport, ConstFieldView molecular_viscosity,
    ConstFieldView effective_viscosity, ConstFieldView heat_capacity,
    FieldView conductivity, FieldView enthalpy_diffusivity) noexcept {
  if (transport.kernel() != TransportKernel::coast_native_air) return {};
  const Int3 cells = molecular_viscosity.interior;
  const auto valid_scalar = [&](ConstFieldView view) noexcept {
    return view.base != nullptr && view.components == 1U &&
           view.interior.x == cells.x && view.interior.y == cells.y &&
           view.interior.z == cells.z;
  };
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0 ||
      !valid_scalar(effective_viscosity) || !valid_scalar(heat_capacity) ||
      !valid_scalar(as_const(conductivity)) ||
      !valid_scalar(as_const(enthalpy_diffusivity))) {
    return {StatusCode::invalid_plan, kProductBinding};
  }
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        double lambda = 0.0;
        double gamma_h = 0.0;
        const Status status = transport.effective_enthalpy_transport(
            molecular_viscosity.unchecked(cell, 0U),
            effective_viscosity.unchecked(cell, 0U),
            heat_capacity.unchecked(cell, 0U), lambda, gamma_h);
        if (!status) return status;
        conductivity.unchecked(cell, 0U) = lambda;
        enthalpy_diffusivity.unchecked(cell, 0U) = gamma_h;
      }
    }
  }
  return {};
}

Status exchange_effective_thermal_ghosts(
    HaloEngine& halo, StageId stage, const BoundaryPlan& boundary,
    FieldView& conductivity, FieldView& enthalpy_diffusivity,
    Status prerequisite) noexcept {
  std::array<FieldView, 2U> fields{conductivity, enthalpy_diffusivity};
  HaloTicket ticket;
  Status status =
      halo.begin(stage, {fields.data(), fields.size()}, prerequisite, ticket);
  if (status)
    status = halo.finish(ticket, {fields.data(), fields.size()});
  if (status) {
    conductivity = fields[0U];
    enthalpy_diffusivity = fields[1U];
    status = apply_physical_zero_gradient(
        boundary, {fields.data(), fields.size()});
    conductivity = fields[0U];
    enthalpy_diffusivity = fields[1U];
  }
  return status;
}

bool product_entirely_periodic(const BoundaryPlan& boundary) noexcept {
  for (std::size_t index = 0U; index < 6U; ++index) {
    const BoundaryFacePlan* face = nullptr;
    if (!boundary.face(static_cast<CartesianFace>(index), face) ||
        face == nullptr || !face->periodic)
      return false;
  }
  return true;
}

bool product_candidate_boundary_supported(
    const BoundaryPlan& boundary) noexcept {
  if (boundary.pressure_reference() !=
      PressureReferenceKind::boundary_absolute)
    return false;
  for (std::size_t index = 0U; index < 6U; ++index) {
    const BoundaryFacePlan* face = nullptr;
    if (!boundary.face(static_cast<CartesianFace>(index), face) ||
        face == nullptr)
      return false;
    switch (face->flow_kind) {
      case BoundaryKind::velocity_inlet:
      case BoundaryKind::mass_flow_inlet:
      case BoundaryKind::pressure_outlet:
      case BoundaryKind::no_slip_wall:
      case BoundaryKind::moving_wall:
      case BoundaryKind::slip:
      case BoundaryKind::symmetry:
      case BoundaryKind::periodic:
        break;
      case BoundaryKind::none:
      case BoundaryKind::static_state_inlet:
      case BoundaryKind::total_state_inlet:
      case BoundaryKind::nscbc_inlet:
      case BoundaryKind::nscbc_outlet:
      case BoundaryKind::adiabatic_wall:
      case BoundaryKind::isothermal_wall:
      case BoundaryKind::heat_flux_wall:
        return false;
    }
  }
  return true;
}

Status reduce_ibm_reconstruction_audit(
    MPI_Comm communicator, const IbmReconstructionAudit& local,
    IbmReconstructionAudit& global) noexcept {
  if (!local.valid ||
      local.policy > IbmReconstructionPolicy::adaptive_order ||
      local.standard_reach == 0U ||
      local.group_count != local.quadratic_groups + local.linear_groups) {
    return {StatusCode::invalid_plan, kProductCollective};
  }
  std::array<std::uint64_t, 9U> counts{{
      local.group_count,
      local.quadratic_groups,
      local.linear_groups,
      local.expanded_search_groups,
      local.rank_fallback_groups,
      local.condition_fallback_groups,
      local.coverage_fallback_groups,
      local.donor_fallback_groups,
      local.valid ? 1U : 0U,
  }};
  if (MPI_Allreduce(MPI_IN_PLACE, counts.data(),
                    static_cast<int>(counts.size()), MPI_UINT64_T, MPI_SUM,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kProductCollective};
  }
  std::array<unsigned int, 2U> identity{{
      static_cast<unsigned int>(local.policy),
      static_cast<unsigned int>(local.standard_reach),
  }};
  std::array<unsigned int, 2U> minimum = identity;
  std::array<unsigned int, 2U> maximum = identity;
  if (MPI_Allreduce(MPI_IN_PLACE, minimum.data(),
                    static_cast<int>(minimum.size()), MPI_UNSIGNED, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, maximum.data(),
                    static_cast<int>(maximum.size()), MPI_UNSIGNED, MPI_MAX,
                    communicator) != MPI_SUCCESS ||
      minimum != maximum) {
    return {StatusCode::mpi_failure, kProductCollective};
  }
  std::array<double, 2U> maxima{{local.maximum_condition_estimate,
                                 local.maximum_functional_l1}};
  if (MPI_Allreduce(MPI_IN_PLACE, maxima.data(),
                    static_cast<int>(maxima.size()), MPI_DOUBLE, MPI_MAX,
                    communicator) != MPI_SUCCESS ||
      !std::isfinite(maxima[0U]) || !std::isfinite(maxima[1U])) {
    return {StatusCode::mpi_failure, kProductCollective};
  }
  global = {};
  global.valid = true;
  global.policy = local.policy;
  global.standard_reach = local.standard_reach;
  global.group_count = counts[0U];
  global.quadratic_groups = counts[1U];
  global.linear_groups = counts[2U];
  global.expanded_search_groups = counts[3U];
  global.rank_fallback_groups = counts[4U];
  global.condition_fallback_groups = counts[5U];
  global.coverage_fallback_groups = counts[6U];
  global.donor_fallback_groups = counts[7U];
  global.maximum_condition_estimate = maxima[0U];
  global.maximum_functional_l1 = maxima[1U];
  return {};
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
std::atomic<bool> g_candidate_globalization_armed{false};
std::atomic<bool> g_candidate_globalization_published{false};
detail::PressureEnergyCandidateGlobalizationDiagnostic
    g_candidate_globalization_diagnostic{};
std::atomic<bool> g_candidate_storage_published{false};
detail::PressureEnergyCandidateStorageDiagnostic g_candidate_storage_diagnostic{};
std::atomic<bool> g_product_piso_bind_failure_armed{false};
std::atomic<int> g_product_piso_bind_failure_rank{-1};
std::atomic<bool> g_candidate_poison_armed{false};
std::atomic<int> g_candidate_poison_rank{-1};
std::atomic<std::uint8_t> g_candidate_poison_kind{0U};
std::atomic<std::uint8_t> g_candidate_poison_corrector{0U};
std::atomic<bool> g_final_flux_history_published{false};
detail::ProductFinalFluxHistoryDiagnostic g_final_flux_history_diagnostic{};
std::atomic<bool> g_pressure_correction_warm_start_published{false};
std::atomic<bool> g_pressure_correction_warm_start_suppressed{false};
detail::PressureCorrectionWarmStartDiagnostic
    g_pressure_correction_warm_start_diagnostic{};
std::atomic<bool> g_fresh_initialization_published{false};
std::atomic<bool> g_fresh_initialization_poison_armed{false};
std::atomic<int> g_fresh_initialization_poison_rank{-1};
std::atomic<std::uint8_t> g_fresh_initialization_poison_kind{0U};
std::atomic<RevisionToken> g_fresh_initialization_generation{0U};
detail::FreshInitializationDiagnostic g_fresh_initialization_diagnostic{};
std::atomic<bool> g_cold_velocity_dependents_published{false};
detail::ColdVelocityDependentsDiagnostic
    g_cold_velocity_dependents_diagnostic{};
std::atomic<bool> g_restart_restore_allocation_failure_armed{false};
std::atomic<int> g_restart_restore_allocation_failure_rank{-1};
std::atomic<std::uint8_t> g_restart_restore_allocation_failure_point{0U};
std::atomic<int> g_restart_restore_allocation_lowest_failing_rank{-1};

bool consume_candidate_globalization_arm() noexcept {
  return g_candidate_globalization_armed.exchange(false,
                                                   std::memory_order_acq_rel);
}

void publish_candidate_globalization_diagnostic(
    const detail::PressureEnergyCandidateGlobalizationDiagnostic&
        diagnostic) noexcept {
  g_candidate_globalization_diagnostic = diagnostic;
  g_candidate_globalization_published.store(diagnostic.valid,
                                             std::memory_order_release);
}

void publish_candidate_storage_diagnostic(
    const detail::PressureEnergyCandidateStorageDiagnostic& diagnostic) noexcept {
  g_candidate_storage_diagnostic = diagnostic;
  g_candidate_storage_published.store(diagnostic.valid,
                                      std::memory_order_release);
}

void publish_fresh_initialization_diagnostic(
    const detail::FreshInitializationDiagnostic& diagnostic) noexcept {
  g_fresh_initialization_diagnostic = diagnostic;
  g_fresh_initialization_published.store(diagnostic.valid,
                                         std::memory_order_release);
}

bool consume_restart_restore_allocation_failure(
    detail::RestartRestoreAllocationPoint point, int rank) noexcept {
  if (!g_restart_restore_allocation_failure_armed.load(
          std::memory_order_acquire) ||
      g_restart_restore_allocation_failure_point.load(
          std::memory_order_acquire) != static_cast<std::uint8_t>(point))
    return false;
  const bool armed = g_restart_restore_allocation_failure_armed.exchange(
      false, std::memory_order_acq_rel);
  return armed &&
         rank == g_restart_restore_allocation_failure_rank.load(
                     std::memory_order_acquire);
}

bool consume_product_piso_bind_failure(int rank) noexcept {
  if (!g_product_piso_bind_failure_armed.load(std::memory_order_acquire) ||
      rank !=
          g_product_piso_bind_failure_rank.load(std::memory_order_acquire))
    return false;
  return g_product_piso_bind_failure_armed.exchange(
      false, std::memory_order_acq_rel);
}
#endif

std::uint64_t diagnostic_global_cell(Int3 cell, Int3 global) noexcept {
  return static_cast<std::uint64_t>(cell.x) +
         static_cast<std::uint64_t>(global.x) *
             (static_cast<std::uint64_t>(cell.y) +
              static_cast<std::uint64_t>(global.y) *
                  static_cast<std::uint64_t>(cell.z));
}

bool diagnostic_inside(double value, double minimum, double maximum) noexcept {
  return std::isfinite(value) && std::isfinite(minimum) &&
         std::isfinite(maximum) && value >= minimum && value <= maximum;
}

Status publish_numerical_failure_context(
    MPI_Comm communicator, const NumericalFailureContext& local,
    NumericalFailureContext& published) noexcept {
  static_assert(std::is_trivially_copyable<NumericalFailureContext>::value,
                "failure context must remain a fixed-size MPI payload");
  const std::uint64_t absent = std::numeric_limits<std::uint64_t>::max();
  const std::uint64_t local_cell = local.valid ? local.global_cell : absent;
  std::uint64_t selected_cell = absent;
  if (MPI_Allreduce(&local_cell, &selected_cell, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kProductCollective};
  }
  if (selected_cell == absent) {
    published = {};
    return {};
  }
  int rank = -1;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kProductCollective};
  }
  const int local_owner =
      local.valid && local.global_cell == selected_cell
          ? rank
          : std::numeric_limits<int>::max();
  int owner = std::numeric_limits<int>::max();
  if (MPI_Allreduce(&local_owner, &owner, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      owner == std::numeric_limits<int>::max() ||
      sizeof(NumericalFailureContext) >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    return {StatusCode::mpi_failure, kProductCollective};
  }
  NumericalFailureContext candidate = rank == owner ? local
                                                     : NumericalFailureContext{};
  if (MPI_Bcast(&candidate, static_cast<int>(sizeof(candidate)), MPI_BYTE,
                owner, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kProductCollective};
  }
  if (!candidate.valid || candidate.global_cell != selected_cell ||
      candidate.rank != owner) {
    return {StatusCode::invalid_plan, kProductCollective};
  }
  published = candidate;
  return {};
}

struct ProductFields {
  FieldId rho{}, velocity{}, pressure{}, enthalpy{}, temperature{};
  FieldId molecular_viscosity{}, effective_viscosity{}, compressibility{};
  FieldId enthalpy_compressibility{};
  FieldId velocity_gradient{}, thermal_conductivity{}, heat_capacity{};
  FieldId enthalpy_diffusivity{};
  FieldId momentum_diagonal{}, momentum_rhs{}, momentum_residual{};
  FieldId pressure_diagonal{}, pressure_rhs{}, pressure_correction{};
  FieldId r_au{}, h_by_a{}, pressure_gradient{};
  FieldId predictor_density{}, predictor_accepted_advection{};
  FieldId predictor_previous_advection{}, predictor_low_bundle{}, eos_density{};
  FieldId enthalpy_nonadvective_rate{}, scalar_diffusivity{};
  FieldId pressure_energy_c_h{}, pressure_energy_c_h_row_scale{};
  FieldId pressure_energy_e_p{}, pressure_energy_e_h{};
  FieldId pressure_energy_continuity_residual{};
  FieldId pressure_energy_energy_residual{};
  FieldId enthalpy_correction{};
  FieldId pressure_energy_delta_temperature{};
  FieldId pressure_energy_candidate_pressure{};
  FieldId pressure_energy_candidate_pressure_correction{};
  FieldId pressure_energy_candidate_enthalpy{};
  FieldId pressure_energy_candidate_density{};
  FieldId pressure_energy_candidate_temperature{};
  FieldId pressure_energy_candidate_velocity{};
  FieldId pressure_energy_candidate_molecular_viscosity{};
  FieldId pressure_energy_candidate_effective_viscosity{};
  FieldId pressure_energy_candidate_velocity_gradient{};
  FieldId pressure_energy_candidate_compressibility{};
  FieldId pressure_energy_candidate_enthalpy_compressibility{};
  FieldId pressure_energy_candidate_thermal_conductivity{};
  FieldId pressure_energy_candidate_heat_capacity{};
  FieldId pressure_energy_candidate_enthalpy_diffusivity{};
  FieldId schur_continuity_response{}, schur_eliminated_enthalpy{};
  FieldId schur_energy_response{};
  FieldId krylov_vectors{}, krylov_scalars{}, mg_arena{};
  std::vector<FieldId> scalars;
  std::vector<FieldId> pressure_energy_candidate_species;
  std::vector<FieldId> scalar_nonadvective_rates;
  std::vector<TransportedScalarRole> scalar_roles;
};

Status require(FieldRegistry& registry, const char* name,
               std::uint8_t components, std::uint8_t ghosts,
               FieldId& field) noexcept {
  return registry.require_field(name, components, ghosts, field);
}

Status register_fields(FieldRegistry& registry, ProductFields& fields,
                       const ValidatedModel& model, std::uint8_t ghosts,
                       const LinearWorkspaceRequirements& krylov,
                       const LinearWorkspaceRequirements& auxiliary_krylov,
                       const MgWorkspaceRequirements& mg) {
  Status status = require(registry, "rho", 1U, ghosts, fields.rho);
  if (status) status = require(registry, "U", 3U, ghosts, fields.velocity);
  if (status) status = require(registry, "pi", 1U, ghosts, fields.pressure);
  if (status) status = require(registry, "h", 1U, ghosts, fields.enthalpy);
  if (status)
    status = require(registry, "T", 1U, ghosts, fields.temperature);
  if (status)
    status = require(registry, "mu", 1U, ghosts,
                     fields.molecular_viscosity);
  if (status)
    status = require(registry, "mu_eff", 1U, ghosts,
                     fields.effective_viscosity);
  if (status)
    status = require(registry, "drho_dp_hY", 1U, ghosts,
                     fields.compressibility);
  if (status)
    status = require(registry, "drho_dh_pY", 1U, ghosts,
                     fields.enthalpy_compressibility);
  if (status)
    status = require(registry, "grad_U", 9U, ghosts,
                     fields.velocity_gradient);
  if (status)
    status = require(registry, "lambda", 1U, ghosts,
                     fields.thermal_conductivity);
  if (status)
    status = require(registry, "cp", 1U, ghosts, fields.heat_capacity);
  if (status)
    status = require(registry, "lambda_over_cp", 1U, ghosts,
                     fields.enthalpy_diffusivity);
  if (status)
    status = require(registry, "momentum_diagonal", 3U, 0U,
                     fields.momentum_diagonal);
  if (status)
    status = require(registry, "momentum_rhs", 3U, 0U,
                     fields.momentum_rhs);
  if (status)
    status = require(registry, "momentum_residual", 3U, 0U,
                     fields.momentum_residual);
  if (status)
    status = require(registry, "pressure_diagonal", 1U, 0U,
                     fields.pressure_diagonal);
  if (status)
    status = require(registry, "pressure_rhs", 1U, 0U,
                     fields.pressure_rhs);
  if (status)
    status = require(registry, "delta_pi", 1U, ghosts,
                     fields.pressure_correction);
  if (status) status = require(registry, "rAU", 3U, 1U, fields.r_au);
  if (status) status = require(registry, "HbyA", 3U, ghosts, fields.h_by_a);
  if (status)
    status = require(registry, "grad_delta_pi", 3U, 1U,
                     fields.pressure_gradient);
  if (status)
    status = require(registry, "predictor_rho_star", 1U, 0U,
                     fields.predictor_density);
  if (status)
    status = require(registry, "predictor_advection_n", 1U, 0U,
                     fields.predictor_accepted_advection);
  if (status)
    status = require(registry, "predictor_advection_nm1", 1U, 0U,
                     fields.predictor_previous_advection);
  if (status && model.transported_scalars.size() > 253U)
    status = {StatusCode::invalid_case, kProductRegistration};
  if (status)
    status = require(
        registry, "predictor_low_bundle",
        static_cast<std::uint8_t>(2U + model.transported_scalars.size()),
        ghosts, fields.predictor_low_bundle);
  if (status)
    status = require(registry, "terminal_eos_rho", 1U, 0U,
                     fields.eos_density);
  if (status)
    status = require(registry, "rhs_nonadv_h", 1U, 0U,
                     fields.enthalpy_nonadvective_rate);
  if (status)
    status = require(registry, "scalar_mass_diffusivity", 1U, 1U,
                     fields.scalar_diffusivity);
  if (status)
    status = require(registry, "pressure_energy_C_h", 1U, 0U,
                     fields.pressure_energy_c_h);
  if (status)
    status = require(registry, "pressure_energy_C_h_row_scale", 1U, 0U,
                     fields.pressure_energy_c_h_row_scale);
  if (status)
    status = require(registry, "pressure_energy_E_p", 1U, 0U,
                     fields.pressure_energy_e_p);
  if (status)
    status = require(registry, "pressure_energy_E_h", 1U, 0U,
                     fields.pressure_energy_e_h);
  if (status)
    status = require(registry, "pressure_energy_R_C", 1U, 0U,
                     fields.pressure_energy_continuity_residual);
  if (status)
    status = require(registry, "pressure_energy_R_E", 1U, 0U,
                     fields.pressure_energy_energy_residual);
  if (status)
    status = require(registry, "delta_h", 1U, 0U,
                     fields.enthalpy_correction);
  if (status)
    status = require(registry, "pressure_energy_delta_T", 1U, 1U,
                     fields.pressure_energy_delta_temperature);
  if (status)
    status = require(registry, "pressure_energy_candidate_pi", 1U, ghosts,
                     fields.pressure_energy_candidate_pressure);
  if (status)
    status = require(registry, "pressure_energy_candidate_delta_pi", 1U,
                     ghosts,
                     fields.pressure_energy_candidate_pressure_correction);
  if (status)
    status = require(registry, "pressure_energy_candidate_h", 1U, ghosts,
                     fields.pressure_energy_candidate_enthalpy);
  if (status)
    status = require(registry, "pressure_energy_candidate_rho", 1U, ghosts,
                     fields.pressure_energy_candidate_density);
  if (status)
    status = require(registry, "pressure_energy_candidate_T", 1U, ghosts,
                     fields.pressure_energy_candidate_temperature);
  if (status)
    status = require(registry, "pressure_energy_candidate_U", 3U, ghosts,
                     fields.pressure_energy_candidate_velocity);
  if (status)
    status = require(registry, "pressure_energy_candidate_mu", 1U, ghosts,
                     fields.pressure_energy_candidate_molecular_viscosity);
  if (status)
    status = require(registry, "pressure_energy_candidate_mu_eff", 1U, ghosts,
                     fields.pressure_energy_candidate_effective_viscosity);
  if (status)
    status = require(registry, "pressure_energy_candidate_grad_U", 9U, ghosts,
                     fields.pressure_energy_candidate_velocity_gradient);
  if (status)
    status = require(registry, "pressure_energy_candidate_drho_dp_hY", 1U,
                     ghosts,
                     fields.pressure_energy_candidate_compressibility);
  if (status)
    status = require(
        registry, "pressure_energy_candidate_drho_dh_pY", 1U, ghosts,
        fields.pressure_energy_candidate_enthalpy_compressibility);
  if (status)
    status = require(registry, "pressure_energy_candidate_lambda", 1U, ghosts,
                     fields.pressure_energy_candidate_thermal_conductivity);
  if (status)
    status = require(registry, "pressure_energy_candidate_cp", 1U, ghosts,
                     fields.pressure_energy_candidate_heat_capacity);
  if (status)
    status = require(registry, "pressure_energy_candidate_lambda_over_cp", 1U,
                     ghosts,
                     fields.pressure_energy_candidate_enthalpy_diffusivity);
  if (status)
    status = require(registry,
                     "pressure_energy_schur_continuity_response", 1U, 0U,
                     fields.schur_continuity_response);
  if (status)
    status = require(registry, "pressure_energy_schur_eliminated_h", 1U, 2U,
                     fields.schur_eliminated_enthalpy);
  if (status)
    status = require(registry, "pressure_energy_schur_energy_response", 1U,
                     0U, fields.schur_energy_response);
  if (status)
    status = require(
        registry, "krylov_vectors",
        std::max(krylov.vector_slots, auxiliary_krylov.vector_slots),
        std::max(krylov.ghost_width, auxiliary_krylov.ghost_width),
        fields.krylov_vectors);
  if (status)
    status = require(registry, "krylov_scalars", 1U, 0U,
                     fields.krylov_scalars);
  if (status)
    status = require(registry, "mg_arena", 1U, 1U, fields.mg_arena);
  if (!status) return status;
  try {
    fields.scalars.reserve(model.transported_scalars.size());
    fields.pressure_energy_candidate_species.reserve(
        model.transported_scalars.size());
    fields.scalar_nonadvective_rates.reserve(
        model.transported_scalars.size());
    fields.scalar_roles.reserve(model.transported_scalars.size());
    for (const TransportedScalarSpec& scalar : model.transported_scalars) {
      FieldId id = 0U;
      status = registry.require_field(scalar.stable_name, 1U, ghosts, id);
      if (!status) return status;
      fields.scalars.push_back(id);
      fields.scalar_roles.push_back(scalar.role);
      if (scalar.role == TransportedScalarRole::species) {
        FieldId candidate_species = 0U;
        const std::string candidate_name =
            "pressure_energy_candidate_" + scalar.stable_name;
        status = registry.require_field(candidate_name, 1U, ghosts,
                                        candidate_species);
        if (!status) return status;
        fields.pressure_energy_candidate_species.push_back(
            candidate_species);
      }
      FieldId rate = 0U;
      const std::string rate_name = "rhs_nonadv_" + scalar.stable_name;
      status = registry.require_field(rate_name, 1U, 0U, rate);
      if (!status) return status;
      fields.scalar_nonadvective_rates.push_back(rate);
    }
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kProductRegistration};
  }
  (void)mg;
  return {};
}

bool state_field(FieldId field, const ProductFields& fields) noexcept {
  if (field == fields.rho || field == fields.velocity ||
      field == fields.pressure || field == fields.enthalpy ||
      field == fields.temperature ||
      field == fields.enthalpy_nonadvective_rate) {
    return true;
  }
  return std::find(fields.scalars.begin(), fields.scalars.end(), field) !=
             fields.scalars.end() ||
         std::find(fields.scalar_nonadvective_rates.begin(),
                   fields.scalar_nonadvective_rates.end(), field) !=
             fields.scalar_nonadvective_rates.end();
}

Status compile_graph(const ProductFields& fields, std::uint8_t ghosts,
                     Int3 local_shape, std::size_t local_cells,
                     RemoteDonorExchangeStats pressure_donors,
                     RemoteDonorExchangeStats candidate_pressure_donors,
                     RemoteDonorExchangeStats candidate_velocity_donors,
                     RemoteDonorExchangeStats candidate_rate_donors,
                     RemoteDonorExchangeStats gradient_donors,
                     RemoteDonorExchangeStats momentum_donors,
                     RemoteDonorExchangeStats rate_donors,
                     RemoteDonorExchangeStats force_donors,
                     FrozenExecutionGraph& graph) noexcept {
  std::size_t default_workspace = 0U;
  std::size_t momentum_workspace = 0U;
  std::size_t predictor_halo_bytes = 0U;
  std::size_t predictor_donor_halo_bytes = 0U;
  std::size_t thermo_halo_bytes = 0U;
  std::size_t turbulence_halo_bytes = 0U;
  std::size_t momentum_halo_bytes = 0U;
  std::size_t momentum_limiter_halo_bytes = 0U;
  std::size_t pressure_halo_bytes = 0U;
  std::size_t candidate_correction_halo_bytes = 0U;
  std::size_t candidate_state_full_halo_bytes = 0U;
  std::size_t candidate_state_face_halo_bytes = 0U;
  std::size_t candidate_state_halo_bytes = 0U;
  std::size_t thermal_halo_bytes = 0U;
  std::size_t candidate_finalizer_halo_bytes = 0U;
  std::size_t force_halo_bytes = 0U;
  if (!detail::product_field_bytes(local_cells, 8U, default_workspace) ||
      !detail::product_field_bytes(local_cells, 12U, momentum_workspace) ||
      !detail::product_halo_bytes(local_shape,
                                  1U + fields.scalars.size(), ghosts,
                                  predictor_halo_bytes) ||
      !detail::product_halo_bytes(local_shape, 1U, 1U,
                                  predictor_donor_halo_bytes) ||
      !detail::product_halo_bytes(local_shape,
                                  1U + fields.scalars.size(), ghosts,
                                  thermo_halo_bytes) ||
      !detail::product_halo_bytes(local_shape, 6U, ghosts,
                                  turbulence_halo_bytes) ||
      !detail::product_halo_bytes(local_shape, 15U, ghosts,
                                  momentum_halo_bytes) ||
      !detail::product_halo_bytes(local_shape, 3U, 1U,
                                  momentum_limiter_halo_bytes) ||
      !detail::product_halo_bytes(local_shape, 8U, ghosts,
                                  pressure_halo_bytes) ||
      !detail::product_halo_bytes(local_shape, 1U, 1U,
                                  candidate_correction_halo_bytes) ||
      !detail::product_halo_bytes(
          local_shape, 6U + fields.pressure_energy_candidate_species.size(),
          ghosts, candidate_state_full_halo_bytes) ||
      !detail::product_halo_bytes(local_shape, 1U, 1U,
                                  candidate_state_face_halo_bytes) ||
      !detail::product_halo_bytes(local_shape, 2U, 1U,
                                  thermal_halo_bytes) ||
      !detail::product_halo_bytes(
          local_shape, 7U + fields.pressure_energy_candidate_species.size(),
          1U, candidate_finalizer_halo_bytes) ||
      !detail::product_halo_bytes(local_shape,
                                  15U + fields.scalars.size(), ghosts,
                                  force_halo_bytes)) {
    return {StatusCode::invalid_plan, kProductAnalysis};
  }
  if (!detail::product_checked_add(candidate_state_full_halo_bytes,
                                   candidate_state_face_halo_bytes,
                                   candidate_state_halo_bytes))
    return {StatusCode::invalid_plan, kProductAnalysis};
  if (predictor_donor_halo_bytes >
      std::numeric_limits<std::size_t>::max() - predictor_halo_bytes) {
    return {StatusCode::invalid_plan, kProductAnalysis};
  }
  const std::size_t predictor_route_halo_bytes =
      predictor_halo_bytes + predictor_donor_halo_bytes;
  if (pressure_donors.bytes_per_exchange >
          UINT64_MAX - static_cast<std::uint64_t>(pressure_halo_bytes) ||
      gradient_donors.bytes_per_exchange >
          UINT64_MAX - static_cast<std::uint64_t>(turbulence_halo_bytes) ||
      momentum_donors.bytes_per_exchange >
          UINT64_MAX - static_cast<std::uint64_t>(momentum_halo_bytes) ||
      force_donors.bytes_per_exchange >
          UINT64_MAX - static_cast<std::uint64_t>(force_halo_bytes) ||
      rate_donors.bytes_per_exchange >
          UINT64_MAX - static_cast<std::uint64_t>(force_halo_bytes) -
                           force_donors.bytes_per_exchange ||
      gradient_donors.bytes_per_exchange >
          UINT64_MAX - static_cast<std::uint64_t>(force_halo_bytes) -
                           force_donors.bytes_per_exchange -
                           rate_donors.bytes_per_exchange)
    return {StatusCode::invalid_plan, kProductAnalysis};
  constexpr std::uint64_t kCandidateEvaluationsPerCorrector =
      static_cast<std::uint64_t>(
          kPressureEnergyGlobalizationCandidateCount + 2U);
  std::uint64_t candidate_bytes_per_evaluation =
      static_cast<std::uint64_t>(candidate_correction_halo_bytes);
  const auto checked_add_u64 = [](std::uint64_t left, std::uint64_t right,
                                  std::uint64_t& out) noexcept {
    if (right > UINT64_MAX - left) return false;
    out = left + right;
    return true;
  };
  const auto checked_multiply_u64 = [](std::uint64_t left,
                                       std::uint64_t right,
                                       std::uint64_t& out) noexcept {
    if (left != 0U && right > UINT64_MAX / left) return false;
    out = left * right;
    return true;
  };
  const std::uint64_t thermal_halo_messages =
      static_cast<std::uint64_t>(thermal_halo_bytes / sizeof(double));
  const auto add_live_thermal_halo_budget =
      [&](StageResourceSpec& resources,
          std::uint64_t exchange_count) noexcept {
        std::uint64_t added_messages = 0U;
        std::uint64_t added_bytes = 0U;
        std::uint64_t total_messages = 0U;
        std::uint64_t total_bytes = 0U;
        if (!checked_multiply_u64(thermal_halo_messages, exchange_count,
                                  added_messages) ||
            !checked_multiply_u64(
                static_cast<std::uint64_t>(thermal_halo_bytes),
                exchange_count, added_bytes) ||
            !checked_add_u64(resources.merged_halo_messages, added_messages,
                             total_messages) ||
            !checked_add_u64(resources.merged_halo_bytes, added_bytes,
                             total_bytes))
          return false;
        resources.merged_halo_messages = total_messages;
        resources.merged_halo_bytes = total_bytes;
        return true;
      };
  // ResourceContract is deliberately an upper bound.  Reserve all four
  // independent candidate halo routes for every evaluation, even though the
  // finalizer route is inactive for a closed-periodic case.  A transport
  // chunk contains at least one double, so the total payload-double count is
  // also a conservative message cap independent of MPI chunking policy.
  std::uint64_t candidate_messages_per_evaluation = 0U;
  std::uint64_t candidate_halo_bytes = 0U;
  std::uint64_t candidate_halo_messages = 0U;
  if (!checked_add_u64(
          candidate_messages_per_evaluation,
          static_cast<std::uint64_t>(candidate_correction_halo_bytes /
                                     sizeof(double)),
          candidate_messages_per_evaluation) ||
      !checked_add_u64(
          candidate_messages_per_evaluation,
          static_cast<std::uint64_t>(candidate_state_halo_bytes /
                                     sizeof(double)),
          candidate_messages_per_evaluation) ||
      !checked_add_u64(
          candidate_messages_per_evaluation,
          thermal_halo_messages,
          candidate_messages_per_evaluation) ||
      !checked_add_u64(
          candidate_messages_per_evaluation,
          static_cast<std::uint64_t>(candidate_finalizer_halo_bytes /
                                     sizeof(double)),
          candidate_messages_per_evaluation) ||
      !checked_add_u64(candidate_bytes_per_evaluation,
                       static_cast<std::uint64_t>(candidate_state_halo_bytes),
                       candidate_bytes_per_evaluation) ||
      !checked_add_u64(candidate_bytes_per_evaluation,
                       static_cast<std::uint64_t>(thermal_halo_bytes),
                       candidate_bytes_per_evaluation) ||
      !checked_add_u64(
          candidate_bytes_per_evaluation,
          static_cast<std::uint64_t>(candidate_finalizer_halo_bytes),
          candidate_bytes_per_evaluation) ||
      !checked_add_u64(candidate_bytes_per_evaluation,
                       candidate_pressure_donors.bytes_per_exchange,
                       candidate_bytes_per_evaluation) ||
      !checked_add_u64(candidate_bytes_per_evaluation,
                       candidate_velocity_donors.bytes_per_exchange,
                       candidate_bytes_per_evaluation) ||
      !checked_add_u64(candidate_bytes_per_evaluation,
                       candidate_rate_donors.bytes_per_exchange,
                       candidate_bytes_per_evaluation) ||
      !checked_add_u64(candidate_messages_per_evaluation,
                       candidate_pressure_donors.peer_messages,
                       candidate_messages_per_evaluation) ||
      !checked_add_u64(candidate_messages_per_evaluation,
                       candidate_velocity_donors.peer_messages,
                       candidate_messages_per_evaluation) ||
      !checked_add_u64(candidate_messages_per_evaluation,
                       candidate_rate_donors.peer_messages,
                       candidate_messages_per_evaluation) ||
      !checked_multiply_u64(candidate_bytes_per_evaluation,
                            kCandidateEvaluationsPerCorrector,
                            candidate_halo_bytes) ||
      !checked_multiply_u64(candidate_messages_per_evaluation,
                            kCandidateEvaluationsPerCorrector,
                            candidate_halo_messages))
    return {StatusCode::invalid_plan, kProductAnalysis};
  const FieldAccessSpec trial_rho{fields.rho, StateVisibility::trial};
  const FieldAccessSpec trial_u{fields.velocity, StateVisibility::trial};
  const FieldAccessSpec trial_pi{fields.pressure, StateVisibility::trial};
  const FieldAccessSpec trial_h{fields.enthalpy, StateVisibility::trial};
  const FieldAccessSpec trial_t{fields.temperature, StateVisibility::trial};
  const FieldAccessSpec accepted_rho{fields.rho, StateVisibility::accepted};
  const FieldAccessSpec accepted_u{fields.velocity, StateVisibility::accepted};
  const FieldAccessSpec accepted_h{fields.enthalpy,
                                   StateVisibility::accepted};
  const FieldAccessSpec accepted_h_rate{
      fields.enthalpy_nonadvective_rate, StateVisibility::accepted};
  const FieldAccessSpec trial_h_rate{fields.enthalpy_nonadvective_rate,
                                     StateVisibility::trial};
  const FieldAccessSpec trial_mu{fields.molecular_viscosity,
                                 StateVisibility::workspace};
  const FieldAccessSpec trial_mu_eff{fields.effective_viscosity,
                                     StateVisibility::workspace};
  const FieldAccessSpec trial_grad{fields.velocity_gradient,
                                   StateVisibility::workspace};
  const FieldAccessSpec work_lambda{fields.thermal_conductivity,
                                    StateVisibility::workspace};
  const FieldAccessSpec work_cp{fields.heat_capacity,
                                StateVisibility::workspace};
  const FieldAccessSpec work_drho_dp{fields.compressibility,
                                     StateVisibility::workspace};
  const FieldAccessSpec work_drho_dh{fields.enthalpy_compressibility,
                                     StateVisibility::workspace};
  const FieldAccessSpec work_lambda_cp{fields.enthalpy_diffusivity,
                                       StateVisibility::workspace};
  const FieldAccessSpec work_predictor_rho{fields.predictor_density,
                                           StateVisibility::workspace};
  const FieldAccessSpec work_predictor_n{
      fields.predictor_accepted_advection, StateVisibility::workspace};
  const FieldAccessSpec work_predictor_nm1{
      fields.predictor_previous_advection, StateVisibility::workspace};
  const FieldAccessSpec work_predictor_low{
      fields.predictor_low_bundle, StateVisibility::workspace};
  const FieldAccessSpec work_scalar_diffusivity{
      fields.scalar_diffusivity, StateVisibility::workspace};
  const FieldAccessSpec work_r_au{fields.r_au, StateVisibility::workspace};
  const FieldAccessSpec work_h_by_a{fields.h_by_a,
                                    StateVisibility::workspace};
  const FieldAccessSpec work_pressure_gradient{
      fields.pressure_gradient, StateVisibility::workspace};
  const FieldAccessSpec work_pressure_correction{
      fields.pressure_correction, StateVisibility::workspace};
  const FieldAccessSpec work_momentum_diagonal{
      fields.momentum_diagonal, StateVisibility::workspace};
  const FieldAccessSpec work_momentum_rhs{fields.momentum_rhs,
                                          StateVisibility::workspace};
  const FieldAccessSpec work_momentum_residual{
      fields.momentum_residual, StateVisibility::workspace};
  const FieldAccessSpec work_pressure_diagonal{
      fields.pressure_diagonal, StateVisibility::workspace};
  const FieldAccessSpec work_pressure_rhs{fields.pressure_rhs,
                                          StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_c_h{
      fields.pressure_energy_c_h, StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_c_h_row_scale{
      fields.pressure_energy_c_h_row_scale, StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_e_p{
      fields.pressure_energy_e_p, StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_e_h{
      fields.pressure_energy_e_h, StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_r_c{
      fields.pressure_energy_continuity_residual,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_r_e{
      fields.pressure_energy_energy_residual, StateVisibility::workspace};
  const FieldAccessSpec work_enthalpy_correction{
      fields.enthalpy_correction, StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_delta_temperature{
      fields.pressure_energy_delta_temperature,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_pi{
      fields.pressure_energy_candidate_pressure,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_delta_pi{
      fields.pressure_energy_candidate_pressure_correction,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_h{
      fields.pressure_energy_candidate_enthalpy,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_rho{
      fields.pressure_energy_candidate_density,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_t{
      fields.pressure_energy_candidate_temperature,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_u{
      fields.pressure_energy_candidate_velocity,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_mu{
      fields.pressure_energy_candidate_molecular_viscosity,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_mu_eff{
      fields.pressure_energy_candidate_effective_viscosity,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_grad_u{
      fields.pressure_energy_candidate_velocity_gradient,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_drho_dp{
      fields.pressure_energy_candidate_compressibility,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_drho_dh{
      fields.pressure_energy_candidate_enthalpy_compressibility,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_lambda{
      fields.pressure_energy_candidate_thermal_conductivity,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_cp{
      fields.pressure_energy_candidate_heat_capacity,
      StateVisibility::workspace};
  const FieldAccessSpec work_pressure_energy_candidate_lambda_cp{
      fields.pressure_energy_candidate_enthalpy_diffusivity,
      StateVisibility::workspace};
  const FieldAccessSpec work_schur_continuity_response{
      fields.schur_continuity_response, StateVisibility::workspace};
  const FieldAccessSpec work_schur_eliminated_enthalpy{
      fields.schur_eliminated_enthalpy, StateVisibility::workspace};
  const FieldAccessSpec work_schur_energy_response{
      fields.schur_energy_response, StateVisibility::workspace};
  std::vector<GraphFieldSpec> declarations{
      {trial_rho, ghosts, false}, {trial_u, ghosts, true},
      {trial_pi, ghosts, true},  {trial_h, ghosts, false},
      {trial_t, ghosts, false},   {accepted_rho, 0U, true},
      {accepted_u, 0U, true},     {accepted_h, ghosts, true},
      {accepted_h_rate, 0U, true},
      {trial_h_rate, 0U, false},
      {trial_mu, 1U, false},
      {trial_mu_eff, ghosts, false}, {trial_grad, ghosts, false},
      {work_lambda, 1U, false}, {work_cp, 0U, false},
      {work_drho_dp, 0U, false}, {work_drho_dh, 0U, false},
      {work_lambda_cp, 1U, false},
      {work_predictor_rho, 0U, false}, {work_predictor_n, 0U, false},
      {work_predictor_nm1, 1U, false},
      {work_predictor_low, ghosts, false},
      {work_scalar_diffusivity, 1U, false},
      {work_r_au, 1U, true}, {work_h_by_a, ghosts, true},
      {work_pressure_gradient, 1U, true},
      {work_pressure_correction, ghosts, true},
      {work_momentum_diagonal, 0U, false},
      {work_momentum_rhs, 0U, false},
      {work_momentum_residual, 0U, false},
      {work_pressure_diagonal, 0U, false},
      {work_pressure_rhs, 0U, false},
      {work_pressure_energy_c_h, 0U, false},
      {work_pressure_energy_c_h_row_scale, 0U, false},
      {work_pressure_energy_e_p, 0U, true},
      {work_pressure_energy_e_h, 0U, true},
      {work_pressure_energy_r_c, 0U, false},
      {work_pressure_energy_r_e, 0U, false},
      {work_enthalpy_correction, 0U, false},
      {work_pressure_energy_delta_temperature, 1U, false},
      {work_pressure_energy_candidate_pi, ghosts, false},
      {work_pressure_energy_candidate_delta_pi, ghosts, false},
      {work_pressure_energy_candidate_h, ghosts, false},
      {work_pressure_energy_candidate_rho, ghosts, false},
      {work_pressure_energy_candidate_t, ghosts, false},
      {work_pressure_energy_candidate_u, ghosts, false},
      {work_pressure_energy_candidate_mu, ghosts, false},
      {work_pressure_energy_candidate_mu_eff, ghosts, false},
      {work_pressure_energy_candidate_grad_u, ghosts, false},
      {work_pressure_energy_candidate_drho_dp, ghosts, false},
      {work_pressure_energy_candidate_drho_dh, ghosts, false},
      {work_pressure_energy_candidate_lambda, ghosts, false},
      {work_pressure_energy_candidate_cp, ghosts, false},
      {work_pressure_energy_candidate_lambda_cp, ghosts, false},
      {work_schur_continuity_response, 0U, false},
      {work_schur_eliminated_enthalpy, 2U, false},
      {work_schur_energy_response, 0U, false},
      {{fields.velocity, StateVisibility::committed_snapshot}, 0U, true},
      {{fields.pressure, StateVisibility::committed_snapshot}, 0U, true},
      {{fields.enthalpy, StateVisibility::committed_snapshot}, 0U, true}};
  for (FieldId scalar : fields.scalars) {
    declarations.push_back(
        {{scalar, StateVisibility::trial}, ghosts, false});
    declarations.push_back(
        {{scalar, StateVisibility::accepted}, ghosts, true});
    declarations.push_back(
        {{scalar, StateVisibility::committed_snapshot}, 0U, true});
  }
  for (FieldId scalar : fields.pressure_energy_candidate_species)
    declarations.push_back(
        {{scalar, StateVisibility::workspace}, ghosts, false});
  for (FieldId rate : fields.scalar_nonadvective_rates) {
    declarations.push_back({{rate, StateVisibility::trial}, 0U, false});
    declarations.push_back({{rate, StateVisibility::accepted}, 0U, true});
  }
  ExecutionGraphCompiler compiler;
  Status status = compiler.configure({declarations.data(), declarations.size()});
  if (!status) return status;
  const auto register_mutating = [&](StageId id,
                                     std::vector<FieldAccessSpec> reads,
                                     std::vector<FieldAccessSpec> writes,
                                     std::vector<FieldAccessSpec> ghost_fields,
                                     StageResourceSpec resources,
                                     bool consensus) {
    std::vector<std::uint8_t> widths(ghost_fields.size(), ghosts);
    StageSpec stage;
    stage.id = id;
    stage.reads = {reads.data(), reads.size()};
    stage.writes = {writes.data(), writes.size()};
    stage.ghosts = {ghost_fields.data(), ghost_fields.size()};
    stage.ghost_widths = {widths.data(), widths.size()};
    std::vector<FieldAccessSpec> invalidates;
    for (FieldAccessSpec write : writes) {
      if (std::find(reads.begin(), reads.end(), write) != reads.end())
        invalidates.push_back(write);
    }
    stage.invalidates = {invalidates.data(), invalidates.size()};
    stage.workspace_bytes = default_workspace;
    stage.workspace_alignment = 64U;
    stage.resources = resources;
    stage.collective_consensus = consensus;
    return compiler.register_stage(stage);
  };
  std::vector<FieldAccessSpec> predictor_reads{
      accepted_rho, accepted_h, accepted_h_rate};
  std::vector<FieldAccessSpec> predictor_writes{
      trial_rho, trial_h, work_predictor_rho, work_predictor_n,
      work_predictor_nm1, work_predictor_low};
  std::vector<FieldAccessSpec> predictor_ghosts{accepted_h};
  for (std::size_t index = 0U; index < fields.scalars.size(); ++index) {
    predictor_reads.push_back(
        {fields.scalars[index], StateVisibility::accepted});
    predictor_reads.push_back(
        {fields.scalar_nonadvective_rates[index], StateVisibility::accepted});
    predictor_writes.push_back(
        {fields.scalars[index], StateVisibility::trial});
    predictor_ghosts.push_back(
        {fields.scalars[index], StateVisibility::accepted});
  }
  status = register_mutating(
      10U, std::move(predictor_reads), std::move(predictor_writes),
      std::move(predictor_ghosts),
      {12U, predictor_route_halo_bytes, 0U, 0U, 0U, 128U, 0U}, false);
  if (status) {
    StageSpec contribution_one;
    contribution_one.id = 12U;
    status = compiler.register_stage(contribution_one);
  }
  std::vector<FieldAccessSpec> thermo_reads{trial_rho, trial_h, trial_pi,
                                            trial_u};
  std::vector<FieldAccessSpec> thermo_ghosts{trial_h};
  for (FieldId scalar : fields.scalars) {
    thermo_reads.push_back({scalar, StateVisibility::trial});
    thermo_ghosts.push_back({scalar, StateVisibility::trial});
  }
  if (status)
    status = register_mutating(
        15U, std::move(thermo_reads),
        {trial_rho, trial_pi, trial_t, trial_mu, trial_mu_eff, work_lambda,
         work_cp, work_drho_dp, work_drho_dh, work_lambda_cp},
        std::move(thermo_ghosts),
        {6U, thermo_halo_bytes, 0U, 0U, 0U, 0U, 0U}, true);
  if (status) {
    StageResourceSpec turbulence_resources{
        6U + gradient_donors.peer_messages,
        turbulence_halo_bytes + gradient_donors.bytes_per_exchange,
        0U, 0U, 1U, 0U, 0U};
    if (!add_live_thermal_halo_budget(turbulence_resources, 1U))
      return {StatusCode::invalid_plan, kProductAnalysis};
    status = register_mutating(
        20U,
        {trial_rho, trial_u, trial_mu, trial_mu_eff, accepted_rho,
         accepted_u, work_h_by_a},
        {trial_grad, trial_mu_eff, work_h_by_a}, {trial_u, work_h_by_a},
        turbulence_resources, false);
  }
  if (status) {
    StageSpec momentum;
    const std::array reads{trial_rho, trial_u, trial_pi, trial_mu_eff,
                           trial_grad, work_h_by_a};
    const std::array writes{trial_u, work_momentum_diagonal,
                            work_momentum_rhs, work_momentum_residual,
                            work_h_by_a};
    const std::array invalidates{trial_u, work_h_by_a};
    const std::array ghost_fields{trial_rho, trial_u, trial_pi,
                                  trial_mu_eff, trial_grad};
    const std::array<std::uint8_t, 5U> widths{ghosts, ghosts, ghosts, ghosts,
                                              ghosts};
    momentum.id = 30U;
    momentum.reads = {reads.data(), reads.size()};
    momentum.writes = {writes.data(), writes.size()};
    momentum.invalidates = {invalidates.data(), invalidates.size()};
    momentum.ghosts = {ghost_fields.data(), ghost_fields.size()};
    momentum.ghost_widths = {widths.data(), widths.size()};
    momentum.workspace_bytes = momentum_workspace;
    std::uint64_t limiter_messages = 0U;
    std::uint64_t limiter_bytes = 0U;
    std::uint64_t momentum_messages = 0U;
    std::uint64_t momentum_bytes = 0U;
    if (!checked_multiply_u64(6U, 4U, limiter_messages) ||
        !checked_multiply_u64(
            static_cast<std::uint64_t>(momentum_limiter_halo_bytes), 4U,
            limiter_bytes) ||
        !checked_add_u64(6U, momentum_donors.peer_messages,
                         momentum_messages) ||
        !checked_add_u64(momentum_messages, limiter_messages,
                         momentum_messages) ||
        !checked_add_u64(static_cast<std::uint64_t>(momentum_halo_bytes),
                         momentum_donors.bytes_per_exchange,
                         momentum_bytes) ||
        !checked_add_u64(momentum_bytes, limiter_bytes, momentum_bytes))
      return {StatusCode::invalid_plan, kProductAnalysis};
    momentum.resources.merged_halo_messages = momentum_messages;
    momentum.resources.merged_halo_bytes = momentum_bytes;
    momentum.resources.numeric_refills = 1U;
    momentum.resources.linear_iterations = 192U;
    status = compiler.register_stage(momentum);
  }
  StageResourceSpec pressure_resources;
  if (!checked_add_u64(6U + pressure_donors.peer_messages,
                       candidate_halo_messages,
                       pressure_resources.merged_halo_messages) ||
      !checked_add_u64(
          static_cast<std::uint64_t>(pressure_halo_bytes) +
              pressure_donors.bytes_per_exchange,
          candidate_halo_bytes,
          pressure_resources.merged_halo_bytes) ||
      !add_live_thermal_halo_budget(pressure_resources, 1U))
    return {StatusCode::invalid_plan, kProductAnalysis};
  pressure_resources.numeric_refills = 3U;
  pressure_resources.linear_iterations = 400U;
  const auto register_pressure = [&](StageId id,
                                     StageResourceSpec resources) {
    std::vector<FieldAccessSpec> reads{
        trial_rho, trial_u, trial_pi, trial_h, trial_t,
        work_momentum_diagonal, work_momentum_rhs, work_r_au, work_h_by_a,
        work_pressure_gradient, work_pressure_correction, work_drho_dp,
        work_drho_dh,
        work_pressure_energy_e_p, work_pressure_energy_e_h};
    std::vector<FieldAccessSpec> writes{
        trial_rho, trial_u, trial_pi, trial_h, trial_t,
        work_r_au, work_h_by_a, work_drho_dp, work_drho_dh,
        work_pressure_gradient, work_pressure_diagonal, work_pressure_rhs,
        work_pressure_correction, work_pressure_energy_c_h,
        work_pressure_energy_c_h_row_scale, work_pressure_energy_e_p,
        work_pressure_energy_e_h, work_pressure_energy_r_c,
        work_pressure_energy_r_e, work_enthalpy_correction,
        work_pressure_energy_delta_temperature,
        work_pressure_energy_candidate_pi,
        work_pressure_energy_candidate_delta_pi,
        work_pressure_energy_candidate_h,
        work_pressure_energy_candidate_rho,
        work_pressure_energy_candidate_t,
        work_pressure_energy_candidate_u,
        work_pressure_energy_candidate_mu,
        work_pressure_energy_candidate_mu_eff,
        work_pressure_energy_candidate_grad_u,
        work_pressure_energy_candidate_drho_dp,
        work_pressure_energy_candidate_drho_dh,
        work_pressure_energy_candidate_lambda,
        work_pressure_energy_candidate_cp,
        work_pressure_energy_candidate_lambda_cp,
        work_schur_continuity_response, work_schur_eliminated_enthalpy,
        work_schur_energy_response};
    for (FieldId species : fields.pressure_energy_candidate_species)
      writes.push_back({species, StateVisibility::workspace});
    if (id == 50U) {
      reads.push_back(work_pressure_diagonal);
      reads.push_back(work_pressure_rhs);
      reads.push_back(work_pressure_energy_c_h);
      reads.push_back(work_pressure_energy_c_h_row_scale);
      reads.push_back(work_pressure_energy_r_c);
      reads.push_back(work_pressure_energy_r_e);
      reads.push_back(work_enthalpy_correction);
      reads.push_back(work_pressure_energy_delta_temperature);
      reads.push_back(work_pressure_energy_candidate_pi);
      reads.push_back(work_pressure_energy_candidate_delta_pi);
      reads.push_back(work_pressure_energy_candidate_h);
      reads.push_back(work_pressure_energy_candidate_rho);
      reads.push_back(work_pressure_energy_candidate_t);
      reads.push_back(work_pressure_energy_candidate_u);
      reads.push_back(work_pressure_energy_candidate_mu);
      reads.push_back(work_pressure_energy_candidate_mu_eff);
      reads.push_back(work_pressure_energy_candidate_grad_u);
      reads.push_back(work_pressure_energy_candidate_drho_dp);
      reads.push_back(work_pressure_energy_candidate_drho_dh);
      reads.push_back(work_pressure_energy_candidate_lambda);
      reads.push_back(work_pressure_energy_candidate_cp);
      reads.push_back(work_pressure_energy_candidate_lambda_cp);
      reads.push_back(work_schur_continuity_response);
      reads.push_back(work_schur_eliminated_enthalpy);
      reads.push_back(work_schur_energy_response);
      for (FieldId species : fields.pressure_energy_candidate_species)
        reads.push_back({species, StateVisibility::workspace});
    }
    const std::array ghost_fields{trial_rho, work_r_au, work_h_by_a,
                                  work_pressure_gradient};
    // The pressure reconstruction consumes the Cartesian convection reach,
    // not the larger IBM donor-search capacity carried by state fields.
    // Exchanging all four IBM layers here both violated the coupler's sealed
    // halo contract and paid for unused communication on every corrector.
    constexpr std::uint8_t kCartesianEquationReach = 2U;
    const std::array<std::uint8_t, 4U> widths{
        1U, 1U, std::min(ghosts, kCartesianEquationReach), 1U};
    StageSpec stage;
    stage.id = id;
    stage.reads = {reads.data(), reads.size()};
    stage.writes = {writes.data(), writes.size()};
    stage.ghosts = {ghost_fields.data(), ghost_fields.size()};
    stage.ghost_widths = {widths.data(), widths.size()};
    std::vector<FieldAccessSpec> invalidates{
        trial_rho, trial_u, trial_pi, trial_h, trial_t, work_r_au,
        work_h_by_a, work_pressure_gradient, work_pressure_correction,
        work_drho_dp, work_drho_dh,
        work_pressure_energy_e_p, work_pressure_energy_e_h};
    if (id == 50U) {
      invalidates.push_back(work_pressure_diagonal);
      invalidates.push_back(work_pressure_rhs);
      invalidates.push_back(work_pressure_energy_c_h);
      invalidates.push_back(work_pressure_energy_c_h_row_scale);
      invalidates.push_back(work_pressure_energy_r_c);
      invalidates.push_back(work_pressure_energy_r_e);
      invalidates.push_back(work_enthalpy_correction);
      invalidates.push_back(work_pressure_energy_delta_temperature);
      invalidates.push_back(work_pressure_energy_candidate_pi);
      invalidates.push_back(work_pressure_energy_candidate_delta_pi);
      invalidates.push_back(work_pressure_energy_candidate_h);
      invalidates.push_back(work_pressure_energy_candidate_rho);
      invalidates.push_back(work_pressure_energy_candidate_t);
      invalidates.push_back(work_pressure_energy_candidate_u);
      invalidates.push_back(work_pressure_energy_candidate_mu);
      invalidates.push_back(work_pressure_energy_candidate_mu_eff);
      invalidates.push_back(work_pressure_energy_candidate_grad_u);
      invalidates.push_back(work_pressure_energy_candidate_drho_dp);
      invalidates.push_back(work_pressure_energy_candidate_drho_dh);
      invalidates.push_back(work_pressure_energy_candidate_lambda);
      invalidates.push_back(work_pressure_energy_candidate_cp);
      invalidates.push_back(work_pressure_energy_candidate_lambda_cp);
      invalidates.push_back(work_schur_continuity_response);
      invalidates.push_back(work_schur_eliminated_enthalpy);
      invalidates.push_back(work_schur_energy_response);
      for (FieldId species : fields.pressure_energy_candidate_species)
        invalidates.push_back({species, StateVisibility::workspace});
    }
    stage.invalidates = {invalidates.data(), invalidates.size()};
    stage.workspace_bytes = default_workspace;
    stage.workspace_alignment = 64U;
    stage.resources = resources;
    stage.collective_consensus = true;
    return compiler.register_stage(stage);
  };
  if (status) status = register_pressure(40U, pressure_resources);
  if (status) {
    StageSpec contribution_two;
    contribution_two.id = 45U;
    if (!add_live_thermal_halo_budget(contribution_two.resources, 2U))
      return {StatusCode::invalid_plan, kProductAnalysis};
    status = compiler.register_stage(contribution_two);
  }
  if (status) {
    StageResourceSpec second = pressure_resources;
    // The inherited budget covers the initial C2 solve, its complete bounded
    // candidate-globalization loop, and the terminal coupled-state refresh.
    // Every bounded refinement may repeat both the complete candidate loop
    // and one live thermal refresh.
    constexpr std::uint64_t refinement_count =
        static_cast<std::uint64_t>(kPressureEnergyRefinementCapacity);
    std::uint64_t refinement_candidate_bytes = 0U;
    std::uint64_t refinement_candidate_messages = 0U;
    std::uint64_t second_bytes = 0U;
    std::uint64_t second_messages = 0U;
    if (!checked_multiply_u64(candidate_halo_bytes, refinement_count,
                              refinement_candidate_bytes) ||
        !checked_multiply_u64(candidate_halo_messages, refinement_count,
                              refinement_candidate_messages) ||
        !checked_add_u64(second.merged_halo_bytes,
                         refinement_candidate_bytes, second_bytes) ||
        !checked_add_u64(second.merged_halo_messages,
                         refinement_candidate_messages, second_messages))
      return {StatusCode::invalid_plan, kProductAnalysis};
    second.merged_halo_bytes = second_bytes;
    second.merged_halo_messages = second_messages;
    if (!add_live_thermal_halo_budget(second, refinement_count))
      return {StatusCode::invalid_plan, kProductAnalysis};
    second.cache_publishes = 1U;
    status = register_pressure(50U, second);
  }
  if (status) {
    StageSpec force;
    std::vector<FieldAccessSpec> reads{
        trial_u, trial_pi, trial_h, trial_t, trial_grad, trial_mu,
        trial_mu_eff, work_lambda, work_predictor_n};
    std::vector<FieldAccessSpec> writes{
        trial_h_rate, work_predictor_n, work_scalar_diffusivity};
    std::vector<FieldAccessSpec> ghosts_force{
        trial_pi, trial_h, trial_t, trial_grad, trial_mu, trial_mu_eff,
        work_lambda};
    for (std::size_t index = 0U; index < fields.scalars.size(); ++index) {
      reads.push_back({fields.scalars[index], StateVisibility::trial});
      writes.push_back({fields.scalar_nonadvective_rates[index],
                        StateVisibility::trial});
      ghosts_force.push_back(
          {fields.scalars[index], StateVisibility::trial});
    }
    std::vector<std::uint8_t> widths(ghosts_force.size(), ghosts);
    widths[4U] = 1U;
    widths[6U] = 1U;
    force.id = 60U;
    force.reads = {reads.data(), reads.size()};
    force.writes = {writes.data(), writes.size()};
    force.ghosts = {ghosts_force.data(), ghosts_force.size()};
    force.ghost_widths = {widths.data(), widths.size()};
    const std::array<FieldAccessSpec, 1U> force_invalidates{
        work_predictor_n};
    force.invalidates = {force_invalidates.data(), force_invalidates.size()};
    force.resources.merged_halo_messages =
        6U + gradient_donors.peer_messages + force_donors.peer_messages +
        rate_donors.peer_messages;
    force.resources.merged_halo_bytes =
        force_halo_bytes + gradient_donors.bytes_per_exchange +
        force_donors.bytes_per_exchange + rate_donors.bytes_per_exchange;
    if (!add_live_thermal_halo_budget(force.resources, 1U))
      return {StatusCode::invalid_plan, kProductAnalysis};
    force.resources.cache_publishes = 1U;
    force.collective_consensus = true;
    status = compiler.register_stage(force);
  }
  if (status) {
    StageSpec commit;
    commit.id = 70U;
    commit.kind = StageKind::commit;
    status = compiler.register_stage(commit);
  }
  const std::array<FieldAccessSpec, 3U> base_snapshots{{
      {fields.velocity, StateVisibility::committed_snapshot},
      {fields.pressure, StateVisibility::committed_snapshot},
      {fields.enthalpy, StateVisibility::committed_snapshot}}};
  for (std::size_t service = 0U; service < 5U && status; ++service) {
    StageSpec stage;
    stage.id = static_cast<StageId>(200U + service);
    stage.kind = StageKind::service;
    stage.reads = {base_snapshots.data(), base_snapshots.size()};
    stage.collective_consensus = true;
    status = compiler.register_stage(stage);
  }
  return status ? compiler.freeze(graph) : status;
}

Status collective_semantic(MPI_Comm communicator,
  PlanFingerprint fingerprint) noexcept {
  PlanFingerprint minimum = 0U;
  PlanFingerprint maximum = 0U;
  const int minimum_rc = MPI_Allreduce(&fingerprint, &minimum, 1,
                                       MPI_UINT64_T, MPI_MIN, communicator);
  const int maximum_rc = MPI_Allreduce(&fingerprint, &maximum, 1,
                                       MPI_UINT64_T, MPI_MAX, communicator);
  if (minimum_rc != MPI_SUCCESS || maximum_rc != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kProductCollective};
  }
  return minimum == maximum && minimum != 0U
             ? Status{}
             : Status{StatusCode::invalid_plan, kProductCollective};
}

const FieldDescriptor* find_descriptor(const FieldSchema& schema,
                                       FieldId field) noexcept {
  for (const FieldDescriptor& descriptor : schema) {
    if (descriptor.id == field) return &descriptor;
  }
  return nullptr;
}

Status reserve_stage_halo(MPI_Comm communicator, const MeshPatch& patch,
                          HaloTopology topology,
                          const FieldSchema& schema,
                          const FrozenExecutionGraph& graph, StageId stage_id,
                          HaloEngine& out) noexcept {
  const FrozenStage* stage = graph.stage(stage_id);
  if (stage == nullptr) return {StatusCode::invalid_plan, kProductCommunication};
  const Span<const FieldAccessSpec> ghosts = graph.ghosts(stage_id);
  const Span<const std::uint8_t> widths = graph.ghost_widths(stage_id);
  if (ghosts.size == 0U || ghosts.size != widths.size) {
    return {StatusCode::invalid_plan, kProductCommunication};
  }
  try {
    std::vector<HaloFieldSpec> fields;
    fields.reserve(ghosts.size);
    for (std::size_t index = 0U; index < ghosts.size; ++index) {
      const FieldDescriptor* descriptor =
          find_descriptor(schema, ghosts.data[index].field);
      if (descriptor == nullptr || widths.data[index] == 0U ||
          widths.data[index] > descriptor->ghost_width) {
        return {StatusCode::invalid_plan, kProductCommunication};
      }
      fields.push_back({descriptor->id, widths.data[index],
                        descriptor->components});
    }
    Status status = out.reserve(communicator, patch,
                                {fields.data(), fields.size()}, topology);
    if (!status) return status;
    const HaloPlanStats stats = out.plan_stats();
    if (stats.maximum_messages_per_exchange >
            stage->resources.merged_halo_messages ||
        stats.maximum_bytes_per_exchange >
            stage->resources.merged_halo_bytes) {
      return {StatusCode::invalid_plan,
              static_cast<std::uint32_t>(kProductCommunication + stage_id)};
    }
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kProductCommunication};
  } catch (...) {
    return {StatusCode::invalid_plan, kProductCommunication};
  }
}

Status make_pressure_face_views(std::vector<double>& storage, Int3 cells,
                                FaceFieldView& x, FaceFieldView& y,
                                FaceFieldView& z) noexcept {
  std::size_t expected = 0U;
  if (!detail::product_face_doubles(cells, expected) ||
      storage.size() != expected || storage.data() == nullptr)
    return {StatusCode::invalid_plan, kProductBinding};
  const Int3 x_extent{cells.x + 1, cells.y, cells.z};
  const Int3 y_extent{cells.x, cells.y + 1, cells.z};
  const Int3 z_extent{cells.x, cells.y, cells.z + 1};
  const std::size_t x_values = static_cast<std::size_t>(x_extent.x) *
                               x_extent.y * x_extent.z;
  const std::size_t y_values = static_cast<std::size_t>(y_extent.x) *
                               y_extent.y * y_extent.z;
  const StorageIdentity identity =
      static_cast<StorageIdentity>(reinterpret_cast<std::uintptr_t>(
          storage.data()));
  const RevisionDomainIdentity domain = identity ^ UINT64_C(0xa04face0);
  if (identity == 0U || domain == 0U) return {StatusCode::invalid_plan,
                                             kProductBinding};
  x = {storage.data(), x_extent, static_cast<std::size_t>(x_extent.x),
       static_cast<std::size_t>(x_extent.x) * x_extent.y, CartesianAxis::x,
       identity, domain};
  y = {storage.data() + x_values, y_extent,
       static_cast<std::size_t>(y_extent.x),
       static_cast<std::size_t>(y_extent.x) * y_extent.y, CartesianAxis::y,
       identity, domain};
  z = {storage.data() + x_values + y_values, z_extent,
       static_cast<std::size_t>(z_extent.x),
       static_cast<std::size_t>(z_extent.x) * z_extent.y, CartesianAxis::z,
       identity, domain};
  return {};
}

Status make_pressure_energy_compiled_cell_views(
    std::vector<double>& storage, Int3 cells, FieldId semantic_field,
    RevisionToken revision, FieldView& local_diagonal,
    FieldView& response_stage) noexcept {
  std::size_t cell_values = 0U;
  std::size_t expected = 0U;
  if (!detail::product_cell_count(cells, cell_values) ||
      !detail::product_checked_multiply(cell_values, 2U, expected) ||
      storage.size() != expected || storage.data() == nullptr ||
      semantic_field == 0U || revision == 0U) {
    return {StatusCode::invalid_plan, kProductBinding};
  }
  const StorageIdentity identity =
      static_cast<StorageIdentity>(reinterpret_cast<std::uintptr_t>(
          storage.data()));
  const RevisionDomainIdentity domain = identity ^ UINT64_C(0xa04e4c0d);
  if (identity == 0U || domain == 0U) {
    return {StatusCode::invalid_plan, kProductBinding};
  }
  const auto make = [&](double* base) noexcept {
    FieldView view;
    view.base = base;
    view.interior = cells;
    view.ghosts = {};
    view.components = 1U;
    view.stride_y = static_cast<std::size_t>(cells.x);
    view.stride_z = view.stride_y * static_cast<std::size_t>(cells.y);
    view.component_stride = cell_values;
    view.replica = 0U;
    view.field = semantic_field;
    view.revision = revision;
    view.storage_identity = identity;
    view.revision_domain = domain;
    return view;
  };
  local_diagonal = make(storage.data());
  response_stage = make(storage.data() + cell_values);
  return {};
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
bool candidate_field_range(const FieldView& view,
                           const ArenaFieldLayout& layout,
                           std::uintptr_t& begin,
                           std::uintptr_t& end) noexcept {
  if (view.base == nullptr || layout.offset_doubles < layout.raw_offset_doubles ||
      layout.replica_stride_doubles == 0U ||
      view.replica >= layout.replicas)
    return false;
  const std::size_t interior_offset =
      layout.offset_doubles - layout.raw_offset_doubles;
  const std::uintptr_t base =
      reinterpret_cast<std::uintptr_t>(view.base);
  const std::size_t prefix_bytes = interior_offset * sizeof(double);
  const std::size_t span_bytes =
      layout.replica_stride_doubles * sizeof(double);
  if (base < prefix_bytes ||
      base - prefix_bytes >
          std::numeric_limits<std::uintptr_t>::max() - span_bytes)
    return false;
  begin = base - prefix_bytes;
  end = begin + span_bytes;
  return begin < end;
}

Status capture_candidate_field_storage(
    StateLayers& layers, const ArenaLayout& layout, FieldId candidate_field,
    FieldId live_field, bool live_is_trial,
    detail::PressureEnergyCandidateFieldStorageDiagnostic& out) noexcept {
  FieldView candidate;
  FieldView live;
  Status status = layers.runtime_view(FieldLifetime::persistent_workspace,
                                      candidate_field, candidate);
  if (status) {
    status = live_is_trial
                 ? layers.view(StateRole::trial, live_field, live)
                 : layers.runtime_view(FieldLifetime::persistent_workspace,
                                       live_field, live);
  }
  const ArenaFieldLayout* candidate_layout = layout.field(candidate_field);
  const ArenaFieldLayout* live_layout = layout.field(live_field);
  if (!status || candidate_layout == nullptr || live_layout == nullptr ||
      candidate.ghosts.x != candidate.ghosts.y ||
      candidate.ghosts.x != candidate.ghosts.z ||
      !candidate_field_range(candidate, *candidate_layout,
                             out.candidate_begin, out.candidate_end) ||
      !candidate_field_range(live, *live_layout, out.trial_begin,
                             out.trial_end)) {
    return status ? Status{StatusCode::invalid_plan, kProductBinding} : status;
  }
  out.candidate_field = candidate.field;
  out.trial_field = live.field;
  out.candidate_base = reinterpret_cast<std::uintptr_t>(candidate.base);
  out.trial_base = reinterpret_cast<std::uintptr_t>(live.base);
  out.candidate_revision = candidate.revision;
  out.trial_revision = live.revision;
  out.candidate_storage = candidate.storage_identity;
  out.trial_storage = live.storage_identity;
  out.candidate_revision_domain = candidate.revision_domain;
  out.trial_revision_domain = live.revision_domain;
  out.components = candidate.components;
  out.ghost_width = static_cast<std::uint8_t>(candidate.ghosts.x);
  const bool disjoint = out.candidate_end <= out.trial_begin ||
                        out.trial_end <= out.candidate_begin;
  if (out.candidate_field == out.trial_field ||
      out.candidate_revision == 0U || out.trial_revision == 0U ||
      out.candidate_revision == out.trial_revision ||
      out.candidate_storage == 0U ||
      out.candidate_storage != out.trial_storage ||
      out.candidate_revision_domain == 0U ||
      out.candidate_revision_domain != out.trial_revision_domain ||
      !disjoint)
    return {StatusCode::invalid_plan, kProductBinding};
  return {};
}

bool candidate_ranges_disjoint(
    const detail::PressureEnergyCandidateFieldStorageDiagnostic& left,
    const detail::PressureEnergyCandidateFieldStorageDiagnostic& right) noexcept {
  return left.candidate_begin != 0U && left.candidate_begin < left.candidate_end &&
         right.candidate_begin != 0U &&
         right.candidate_begin < right.candidate_end &&
         (left.candidate_end <= right.candidate_begin ||
          right.candidate_end <= left.candidate_begin);
}
#endif

Status make_pressure_mg_activity(
    const EBTopology& topology, Int3 cells,
    std::vector<std::uint8_t>& cell_activity,
    std::vector<std::uint8_t>& x_activity,
    std::vector<std::uint8_t>& y_activity,
    std::vector<std::uint8_t>& z_activity,
    PlanFingerprint& local_fingerprint,
    PlanFingerprint& collective_fingerprint) {
  std::size_t cell_values = 0U;
  if (!detail::product_cell_count(cells, cell_values))
    return {StatusCode::invalid_plan, kProductBinding};
  const Int3 xe{cells.x + 1, cells.y, cells.z};
  const Int3 ye{cells.x, cells.y + 1, cells.z};
  const Int3 ze{cells.x, cells.y, cells.z + 1};
  std::size_t x_values = 0U;
  std::size_t y_values = 0U;
  std::size_t z_values = 0U;
  if (!detail::product_cell_count(xe, x_values) ||
      !detail::product_cell_count(ye, y_values) ||
      !detail::product_cell_count(ze, z_values) ||
      topology.region().data == nullptr ||
      topology.region().size != cell_values) {
    return {StatusCode::invalid_plan, kProductBinding};
  }
  cell_activity.assign(topology.region().data,
                       topology.region().data + topology.region().size);
  x_activity.assign(x_values, 1U);
  y_activity.assign(y_values, 1U);
  z_activity.assign(z_values, 1U);
  const auto offset = [](Int3 shape, Int3 value) noexcept {
    return static_cast<std::size_t>(value.x) +
           static_cast<std::size_t>(shape.x) *
               (static_cast<std::size_t>(value.y) +
                static_cast<std::size_t>(shape.y) *
                    static_cast<std::size_t>(value.z));
  };
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        if (cell_activity[offset(cells, {i, j, k})] != 0U) continue;
        x_activity[offset(xe, {i, j, k})] = 0U;
        x_activity[offset(xe, {i + 1, j, k})] = 0U;
        y_activity[offset(ye, {i, j, k})] = 0U;
        y_activity[offset(ye, {i, j + 1, k})] = 0U;
        z_activity[offset(ze, {i, j, k})] = 0U;
        z_activity[offset(ze, {i, j, k + 1})] = 0U;
      }
    }
  }
  const Span<const ImmersedLink> links = topology.links();
  for (std::size_t index = 0U; index < links.size; ++index) {
    const ImmersedLink& link = links.data[index];
    const Int3 cell = link.fluid_local_index;
    if (cell.x < 0 || cell.y < 0 || cell.z < 0 || cell.x >= cells.x ||
        cell.y >= cells.y || cell.z >= cells.z) {
      return {StatusCode::invalid_plan, kProductBinding};
    }
    switch (link.direction) {
      case ImmersedFaceDirection::x_negative:
        x_activity[offset(xe, cell)] = 0U;
        break;
      case ImmersedFaceDirection::x_positive:
        x_activity[offset(xe, {cell.x + 1, cell.y, cell.z})] = 0U;
        break;
      case ImmersedFaceDirection::y_negative:
        y_activity[offset(ye, cell)] = 0U;
        break;
      case ImmersedFaceDirection::y_positive:
        y_activity[offset(ye, {cell.x, cell.y + 1, cell.z})] = 0U;
        break;
      case ImmersedFaceDirection::z_negative:
        z_activity[offset(ze, cell)] = 0U;
        break;
      case ImmersedFaceDirection::z_positive:
        z_activity[offset(ze, {cell.x, cell.y, cell.z + 1})] = 0U;
        break;
    }
  }
  std::uint64_t local = kFnvOffset;
  const auto hash_values = [&](const std::vector<std::uint8_t>& values) {
    local = detail::product_mix(local, values.size());
    for (const std::uint8_t value : values)
      local = detail::product_mix(local, value);
  };
  hash_values(cell_activity);
  hash_values(x_activity);
  hash_values(y_activity);
  hash_values(z_activity);
  local = detail::product_mix(local, topology.fingerprint());
  local_fingerprint = local == 0U ? 1U : local;
  std::uint64_t collective = kFnvOffset;
  collective = detail::product_mix(collective,
                                   topology.geometry_fingerprint());
  collective = detail::product_mix(collective,
                                   topology.surface_fingerprint());
  collective = detail::product_mix(collective,
                                   topology.geometry_revision());
  collective_fingerprint = collective == 0U ? 1U : collective;
  return {};
}

void fill_field(FieldView view, Span<const double> components) noexcept {
  for (std::uint8_t component = 0U; component < view.components; ++component)
    for (std::int32_t z = -view.ghosts.z;
         z < view.interior.z + view.ghosts.z; ++z)
      for (std::int32_t y = -view.ghosts.y;
           y < view.interior.y + view.ghosts.y; ++y)
        for (std::int32_t x = -view.ghosts.x;
             x < view.interior.x + view.ghosts.x; ++x)
          view.unchecked({x, y, z}, component) = components.data[component];
}

void fill_field(FieldView view, double value) noexcept {
  for (std::uint8_t component = 0U; component < view.components; ++component)
    for (std::int32_t z = -view.ghosts.z;
         z < view.interior.z + view.ghosts.z; ++z)
      for (std::int32_t y = -view.ghosts.y;
           y < view.interior.y + view.ghosts.y; ++y)
        for (std::int32_t x = -view.ghosts.x;
             x < view.interior.x + view.ghosts.x; ++x)
          view.unchecked({x, y, z}, component) = value;
}

Int3 boundary_owner_cell(CartesianFace face, Int3 cells,
                         std::int32_t inner,
                         std::int32_t outer) noexcept {
  switch (face) {
    case CartesianFace::x_min:
      return {0, inner, outer};
    case CartesianFace::x_max:
      return {cells.x - 1, inner, outer};
    case CartesianFace::y_min:
      return {inner, 0, outer};
    case CartesianFace::y_max:
      return {inner, cells.y - 1, outer};
    case CartesianFace::z_min:
      return {inner, outer, 0};
    case CartesianFace::z_max:
      return {inner, outer, cells.z - 1};
  }
  return {};
}

Int3 boundary_first_ghost_cell(CartesianFace face, Int3 cells,
                               std::int32_t inner,
                               std::int32_t outer) noexcept {
  switch (face) {
    case CartesianFace::x_min:
      return {-1, inner, outer};
    case CartesianFace::x_max:
      return {cells.x, inner, outer};
    case CartesianFace::y_min:
      return {inner, -1, outer};
    case CartesianFace::y_max:
      return {inner, cells.y, outer};
    case CartesianFace::z_min:
      return {inner, outer, -1};
    case CartesianFace::z_max:
      return {inner, outer, cells.z};
  }
  return {};
}

Status local_pressure_outlet_closure(
    const BoundaryPlan& boundary,
    const std::array<BoundaryFaceSpec, 6U>& boundary_specs,
    ConstFieldView pressure_perturbation, double pressure_reference,
    double& maximum, std::uint64_t& samples) noexcept {
  maximum = 0.0;
  samples = 0U;
  const Int3 cells = boundary.local_cells();
  if (pressure_perturbation.base == nullptr ||
      pressure_perturbation.field != boundary.pressure_field() ||
      pressure_perturbation.components != 1U ||
      pressure_perturbation.interior.x != cells.x ||
      pressure_perturbation.interior.y != cells.y ||
      pressure_perturbation.interior.z != cells.z ||
      pressure_perturbation.ghosts.x < 1 ||
      pressure_perturbation.ghosts.y < 1 ||
      pressure_perturbation.ghosts.z < 1 ||
      !std::isfinite(pressure_reference) || pressure_reference <= 0.0) {
    return {StatusCode::invalid_plan, kProductBinding};
  }
  const Span<const BoundaryIndexSpan> spans = boundary.spans();
  for (std::size_t span_index = 0U; span_index < spans.size; ++span_index) {
    const BoundaryIndexSpan& span = spans.data[span_index];
    if (span.stage != BoundaryStage::pressure) continue;
    const std::size_t face_index = static_cast<std::size_t>(span.face);
    if (face_index >= boundary_specs.size())
      return {StatusCode::invalid_plan, kProductBinding};
    const BoundaryFaceSpec& spec = boundary_specs[face_index];
    if (spec.flow_kind != BoundaryKind::pressure_outlet) continue;
    if (span.relation != BoundaryRelation::dirichlet ||
        span.value_source != BoundaryValueSource::resolved_scalar ||
        span.field != pressure_perturbation.field ||
        span.component_begin != 0U || span.component_count != 1U ||
        !std::isfinite(spec.pressure) || spec.pressure <= 0.0) {
      return {StatusCode::invalid_plan, kProductBinding};
    }
    for (std::uint32_t outer = 0U; outer < span.tangent_outer_count;
         ++outer) {
      for (std::uint32_t inner = 0U; inner < span.tangent_inner_count;
           ++inner) {
        const Int3 owner = boundary_owner_cell(
            span.face, cells, static_cast<std::int32_t>(inner),
            static_cast<std::int32_t>(outer));
        const Int3 ghost = boundary_first_ghost_cell(
            span.face, cells, static_cast<std::int32_t>(inner),
            static_cast<std::int32_t>(outer));
        const double face_absolute_pressure =
            pressure_reference +
            0.5 * (pressure_perturbation.unchecked(owner, 0U) +
                   pressure_perturbation.unchecked(ghost, 0U));
        double residual = 0.0;
        if (hf_coast_common_terminal_outlet_v1(
                face_absolute_pressure, spec.pressure, &residual) != 0)
          return {StatusCode::numerical_failure, kProductBinding};
        maximum = std::max(maximum, residual);
        if (samples == std::numeric_limits<std::uint64_t>::max())
          return {StatusCode::invalid_plan, kProductBinding};
        ++samples;
      }
    }
  }
  return {};
}

const BoundaryIndexSpan* scalar_boundary_span(
    const BoundaryPlan& boundary, CartesianFace face,
    FieldId field) noexcept {
  const Span<const BoundaryIndexSpan> spans = boundary.spans();
  for (std::size_t index = 0U; index < spans.size; ++index) {
    const BoundaryIndexSpan& span = spans.data[index];
    if (span.stage == BoundaryStage::scalar && span.face == face &&
        span.field == field)
      return &span;
  }
  return nullptr;
}

Real3 boundary_outward_normal(CartesianFace face) noexcept {
  switch (face) {
    case CartesianFace::x_min:
      return {-1.0, 0.0, 0.0};
    case CartesianFace::x_max:
      return {1.0, 0.0, 0.0};
    case CartesianFace::y_min:
      return {0.0, -1.0, 0.0};
    case CartesianFace::y_max:
      return {0.0, 1.0, 0.0};
    case CartesianFace::z_min:
      return {0.0, 0.0, -1.0};
    case CartesianFace::z_max:
      return {0.0, 0.0, 1.0};
  }
  return {};
}

double boundary_face_area(const CartesianGeometryPlan& geometry,
                          const MeshPatch& patch, CartesianFace face,
                          Int3 owner) noexcept {
  const Int3 global{patch.begin.x + owner.x, patch.begin.y + owner.y,
                    patch.begin.z + owner.z};
  const Span<const double> dx = geometry.x().widths();
  const Span<const double> dy = geometry.y().widths();
  const Span<const double> dz = geometry.z().widths();
  if (global.x < 0 || global.y < 0 || global.z < 0 ||
      static_cast<std::size_t>(global.x) >= dx.size ||
      static_cast<std::size_t>(global.y) >= dy.size ||
      static_cast<std::size_t>(global.z) >= dz.size)
    return std::numeric_limits<double>::quiet_NaN();
  if (face == CartesianFace::x_min || face == CartesianFace::x_max)
    return dy.data[global.y] * dz.data[global.z];
  if (face == CartesianFace::y_min || face == CartesianFace::y_max)
    return dx.data[global.x] * dz.data[global.z];
  return dx.data[global.x] * dy.data[global.y];
}

double vector_dot(Real3 left, Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

Status resolve_static_boundary_values(
    MPI_Comm communicator, const BoundaryPlan& boundary,
    const std::array<BoundaryFaceSpec, 6U>& boundary_specs,
    const CartesianGeometryPlan& geometry, MeshPatch patch,
    const ThermodynamicsPlan& thermodynamics, double pressure_reference,
    ConstFieldView density, ConstFieldView velocity, ConstFieldView enthalpy,
    Span<const ConstFieldView> species,
    Span<const ConstFieldView> passive,
    Span<double> composition, std::vector<double>& scalars,
    std::vector<Real3>& vectors,
    const std::vector<double>& normal_gradients,
    bool resolve_vectors, Status prerequisite = {}) noexcept {
  if (communicator == MPI_COMM_NULL)
    return {StatusCode::invalid_plan, kProductBinding};
  Status local = prerequisite;
  if (local &&
      (!std::isfinite(pressure_reference) || pressure_reference <= 0.0 ||
       species.size != thermodynamics.independent_species_count() ||
       composition.size != species.size ||
       scalars.size() != boundary.resolved_scalar_count() ||
       vectors.size() != boundary.resolved_vector_count() ||
       normal_gradients.size() !=
           boundary.resolved_normal_gradient_count()))
    local = {StatusCode::invalid_plan, kProductBinding};

  // Flux-derived thermal/species gradients still require material-field
  // resolution. Fail closed instead of publishing zero gradients.
  if (local && !normal_gradients.empty())
    local = {StatusCode::invalid_plan, kProductBinding};

  const Span<const BoundaryIndexSpan> spans = boundary.spans();
  const Span<const double> temperatures = boundary.temperature_targets();
  const Span<const double> backflow_temperatures =
      boundary.backflow_temperature_targets();
  const Span<const double> scalar_targets = boundary.scalar_targets();
  const Span<const double> scalar_backflow_targets =
      boundary.scalar_backflow_targets();
  const auto resolve_scalars = [&]() noexcept {
    if (!local) return;
    for (std::size_t span_index = 0U;
         span_index < spans.size && local; ++span_index) {
      const BoundaryIndexSpan& span = spans.data[span_index];
      if (span.value_source != BoundaryValueSource::resolved_scalar)
        continue;
      const std::size_t begin = span.resolved_begin;
      const std::size_t count = span.resolved_stride;
      if (begin > scalars.size() || count > scalars.size() - begin ||
          span.parameter >= boundary.parameter_count()) {
        local = {StatusCode::invalid_plan, kProductBinding};
        break;
      }
      if (span.stage == BoundaryStage::pressure) {
        double perturbation = 0.0;
        const Status resolved = boundary.pressure_perturbation_target(
            span.face, pressure_reference, perturbation);
        if (!resolved) {
          local = resolved;
          break;
        }
        std::fill_n(scalars.data() + begin, count, perturbation);
        continue;
      }
      if (span.stage == BoundaryStage::scalar) {
        const BoundaryFacePlan* face_plan = nullptr;
        const std::size_t face_index =
            static_cast<std::size_t>(span.face);
        const bool conditional_backflow =
            face_index < boundary_specs.size() &&
            boundary.face(span.face, face_plan) && face_plan != nullptr &&
            face_plan->flow_kind == BoundaryKind::pressure_outlet &&
            boundary_specs[face_index].allow_backflow;
        ConstFieldView transported_scalar;
        for (std::size_t candidate = 0U; candidate < species.size;
             ++candidate)
          if (species.data[candidate].field == span.field) {
            transported_scalar = species.data[candidate];
            break;
          }
        if (transported_scalar.base == nullptr)
          for (std::size_t candidate = 0U; candidate < passive.size;
               ++candidate)
            if (passive.data[candidate].field == span.field) {
              transported_scalar = passive.data[candidate];
              break;
            }
        if (!conditional_backflow || transported_scalar.base == nullptr ||
            span.parameter >= scalar_backflow_targets.size) {
          local = {StatusCode::invalid_plan, kProductBinding};
          break;
        }
        const Real3 outward = boundary_outward_normal(span.face);
        for (std::uint32_t outer = 0U;
             outer < span.tangent_outer_count && local; ++outer) {
          for (std::uint32_t inner = 0U;
               inner < span.tangent_inner_count; ++inner) {
            const std::size_t face_cell =
                static_cast<std::size_t>(outer) *
                    span.tangent_inner_count +
                inner;
            const Int3 owner = boundary_owner_cell(
                span.face, boundary.local_cells(),
                static_cast<std::int32_t>(inner),
                static_cast<std::int32_t>(outer));
            const Real3 owner_velocity{
                velocity.unchecked(owner, 0U),
                velocity.unchecked(owner, 1U),
                velocity.unchecked(owner, 2U)};
            const bool backflow =
                vector_dot(owner_velocity, outward) < 0.0;
            scalars[begin + face_cell] =
                backflow
                    ? scalar_backflow_targets.data[span.parameter]
                    : transported_scalar.unchecked(owner, 0U);
          }
        }
        continue;
      }
      if (span.stage != BoundaryStage::enthalpy ||
          span.parameter >= temperatures.size) {
        local = {StatusCode::invalid_plan, kProductBinding};
        break;
      }

      const BoundaryFacePlan* face_plan = nullptr;
      if (!boundary.face(span.face, face_plan) || face_plan == nullptr) {
        local = {StatusCode::invalid_plan, kProductBinding};
        break;
      }
      const bool conditional_backflow =
          face_plan->flow_kind == BoundaryKind::pressure_outlet &&
          boundary_specs[static_cast<std::size_t>(span.face)].allow_backflow;
      const Real3 outward = boundary_outward_normal(span.face);

      for (std::uint32_t outer = 0U;
           outer < span.tangent_outer_count && local; ++outer) {
        for (std::uint32_t inner = 0U;
             inner < span.tangent_inner_count; ++inner) {
          const std::size_t face_cell =
              static_cast<std::size_t>(outer) * span.tangent_inner_count +
              inner;
          const Int3 owner = boundary_owner_cell(
              span.face, boundary.local_cells(),
              static_cast<std::int32_t>(inner),
              static_cast<std::int32_t>(outer));
          bool backflow = false;
          if (conditional_backflow) {
            const Real3 owner_velocity{
                velocity.unchecked(owner, 0U), velocity.unchecked(owner, 1U),
                velocity.unchecked(owner, 2U)};
            backflow = vector_dot(owner_velocity, outward) < 0.0;
            if (!backflow) {
              scalars[begin + face_cell] = enthalpy.unchecked(owner, 0U);
              continue;
            }
          }
          for (std::size_t species_index = 0U;
               species_index < species.size; ++species_index) {
            double value =
                species.data[species_index].unchecked(owner, 0U);
            const BoundaryIndexSpan* scalar_span = scalar_boundary_span(
                boundary, span.face, species.data[species_index].field);
            if (scalar_span != nullptr) {
              if (backflow &&
                  scalar_span->parameter < scalar_backflow_targets.size) {
                value = scalar_backflow_targets.data[scalar_span->parameter];
              } else if (scalar_span->relation ==
                             BoundaryRelation::dirichlet &&
                         scalar_span->parameter < scalar_targets.size) {
                value = scalar_targets.data[scalar_span->parameter];
              }
            }
            composition.data[species_index] = value;
          }
          double resolved_enthalpy = 0.0;
          double cp = 0.0;
          double gas = 0.0;
          const double temperature =
              backflow
                  ? backflow_temperatures.data[span.parameter]
                  : temperatures.data[span.parameter];
          const Status resolved = thermodynamics.mixture_enthalpy(
              temperature,
              {composition.data, composition.size}, resolved_enthalpy, cp,
              gas);
          if (!resolved) {
            local = resolved;
            break;
          }
          scalars[begin + face_cell] = resolved_enthalpy;
        }
      }
    }
  };

  if (!resolve_vectors) {
    resolve_scalars();
    return local;
  }
  std::array<double, 6U> mass_flow_scales{};
  mass_flow_scales.fill(1.0);
  bool mass_flow_collective_failed = false;
  for (std::size_t face_index = 0U; face_index < 6U; ++face_index) {
    const auto face = static_cast<CartesianFace>(face_index);
    const BoundaryFaceSpec& runtime_spec = boundary_specs[face_index];
    const BoundaryFacePlan* face_plan = nullptr;
    const bool valid_face = boundary.face(face, face_plan) &&
                            face_plan != nullptr;
    if (!valid_face && local)
      local = {StatusCode::invalid_plan, kProductBinding};
    if (runtime_spec.flow_kind != BoundaryKind::mass_flow_inlet &&
        !(runtime_spec.flow_kind == BoundaryKind::pressure_outlet &&
          runtime_spec.allow_backflow))
      continue;
    if (runtime_spec.flow_kind == BoundaryKind::mass_flow_inlet) {
      const Real3 resolved_direction = runtime_spec.direction;
      const Real3 inward = [&]() noexcept {
        const Real3 outward = boundary_outward_normal(face);
        return Real3{-outward.x, -outward.y, -outward.z};
      }();
      const double inward_component = vector_dot(resolved_direction, inward);
      double local_capacity = 0.0;
      if (local && valid_face && face_plan->local_owner) {
        const std::int32_t inner_count =
            face_index < 2U ? patch.cells.y : patch.cells.x;
        const std::int32_t outer_count =
            face_index >= 4U ? patch.cells.y : patch.cells.z;
        for (std::int32_t outer = 0; outer < outer_count; ++outer)
          for (std::int32_t inner = 0; inner < inner_count; ++inner) {
            const Int3 owner =
                boundary_owner_cell(face, patch.cells, inner, outer);
            local_capacity += density.unchecked(owner, 0U) *
                              boundary_face_area(geometry, patch, face, owner) *
                              inward_component;
          }
      }
      double global_capacity = 0.0;
      if (MPI_Allreduce(&local_capacity, &global_capacity, 1, MPI_DOUBLE,
                        MPI_SUM, communicator) != MPI_SUCCESS) {
        // Keep the per-face collective sequence fixed.  A rank-local MPI
        // error must not let that rank skip a later configured inlet while
        // peers continue into its Allreduce.
        mass_flow_collective_failed = true;
      }
      if (local &&
          (!(global_capacity > 0.0) || !std::isfinite(global_capacity) ||
           !std::isfinite(runtime_spec.mass_flow_rate)))
        local = {StatusCode::numerical_failure, kProductBinding};
      if (local)
        mass_flow_scales[face_index] =
            runtime_spec.mass_flow_rate / global_capacity;
    }
  }
  const int local_collective_failure =
      mass_flow_collective_failed ? 1 : 0;
  int collective_failure = 0;
  const int collective_failure_result = MPI_Allreduce(
      &local_collective_failure, &collective_failure, 1, MPI_INT, MPI_MAX,
      communicator);
  if (collective_failure_result != MPI_SUCCESS || collective_failure != 0)
    return {StatusCode::mpi_failure, kProductCollective};

  resolve_scalars();
  if (local) {
    for (std::size_t face_index = 0U; face_index < 6U; ++face_index) {
      const auto face = static_cast<CartesianFace>(face_index);
      const BoundaryFaceSpec& runtime_spec = boundary_specs[face_index];
      const BoundaryFacePlan* face_plan = nullptr;
      if (!boundary.face(face, face_plan) || face_plan == nullptr) {
        local = {StatusCode::invalid_plan, kProductBinding};
        break;
      }
      if (runtime_spec.flow_kind != BoundaryKind::mass_flow_inlet &&
          !(runtime_spec.flow_kind == BoundaryKind::pressure_outlet &&
            runtime_spec.allow_backflow))
        continue;
      const Real3 resolved_direction = runtime_spec.direction;
      const double scale =
          runtime_spec.flow_kind == BoundaryKind::mass_flow_inlet
              ? mass_flow_scales[face_index]
              : 1.0;
      for (std::size_t span_index = 0U; span_index < spans.size;
           ++span_index) {
        const BoundaryIndexSpan& span = spans.data[span_index];
        if (span.stage != BoundaryStage::momentum || span.face != face ||
            span.value_source != BoundaryValueSource::resolved_vector)
          continue;
        const std::size_t begin = span.resolved_begin;
        if (begin > vectors.size() ||
            span.resolved_stride > vectors.size() - begin) {
          local = {StatusCode::invalid_plan, kProductBinding};
          break;
        }
        for (std::uint32_t outer = 0U; outer < span.tangent_outer_count;
             ++outer)
          for (std::uint32_t inner = 0U; inner < span.tangent_inner_count;
               ++inner) {
            const std::size_t face_cell =
                static_cast<std::size_t>(outer) * span.tangent_inner_count +
                inner;
            if (runtime_spec.flow_kind == BoundaryKind::mass_flow_inlet) {
              vectors[begin + face_cell] = {scale * resolved_direction.x,
                                            scale * resolved_direction.y,
                                            scale * resolved_direction.z};
            } else {
              const Int3 owner = boundary_owner_cell(
                  face, patch.cells, static_cast<std::int32_t>(inner),
                  static_cast<std::int32_t>(outer));
              const Real3 owner_velocity{
                  velocity.unchecked(owner, 0U),
                  velocity.unchecked(owner, 1U),
                  velocity.unchecked(owner, 2U)};
              if (vector_dot(owner_velocity, boundary_outward_normal(face)) >=
                  0.0)
                vectors[begin + face_cell] = owner_velocity;
              else
                vectors[begin + face_cell] = runtime_spec.backflow_velocity;
            }
          }
      }
    }
  }

  // Static/total state targets need the characteristic state resolver; any
  // remaining vector slice is therefore unsupported at this node.
  if (local) {
    for (std::size_t span_index = 0U;
         span_index < spans.size && local; ++span_index) {
      const BoundaryIndexSpan& span = spans.data[span_index];
      if (span.value_source == BoundaryValueSource::resolved_vector) {
        const BoundaryFacePlan* face_plan = nullptr;
        if (!boundary.face(span.face, face_plan) || face_plan == nullptr ||
            (boundary_specs[static_cast<std::size_t>(span.face)].flow_kind !=
                 BoundaryKind::mass_flow_inlet &&
             boundary_specs[static_cast<std::size_t>(span.face)].flow_kind !=
                 BoundaryKind::pressure_outlet))
          local = {StatusCode::invalid_plan, kProductBinding};
      }
    }
  }
  return local;
}

}  // namespace

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace detail {

void arm_pressure_energy_candidate_globalization_once_for_test() noexcept {
  g_candidate_globalization_published.store(false, std::memory_order_release);
  g_candidate_globalization_diagnostic = {};
  g_candidate_globalization_armed.store(true, std::memory_order_release);
}

void clear_pressure_energy_candidate_globalization_for_test() noexcept {
  g_candidate_globalization_armed.store(false, std::memory_order_release);
  g_candidate_globalization_published.store(false, std::memory_order_release);
  g_candidate_globalization_diagnostic = {};
}

bool pressure_energy_candidate_globalization_diagnostic_for_test(
    PressureEnergyCandidateGlobalizationDiagnostic& out) noexcept {
  if (!g_candidate_globalization_published.load(std::memory_order_acquire)) {
    out = {};
    return false;
  }
  out = g_candidate_globalization_diagnostic;
  return out.valid;
}

bool pressure_energy_candidate_storage_diagnostic_for_test(
    PressureEnergyCandidateStorageDiagnostic& out) noexcept {
  if (!g_candidate_storage_published.load(std::memory_order_acquire)) {
    out = {};
    return false;
  }
  out = g_candidate_storage_diagnostic;
  return out.valid;
}

void arm_product_piso_bind_failure_once_for_test(int failing_rank) noexcept {
  g_product_piso_bind_failure_rank.store(failing_rank,
                                         std::memory_order_release);
  g_product_piso_bind_failure_armed.store(true, std::memory_order_release);
}

void clear_product_piso_bind_failure_for_test() noexcept {
  g_product_piso_bind_failure_armed.store(false, std::memory_order_release);
  g_product_piso_bind_failure_rank.store(-1, std::memory_order_release);
}

void arm_pressure_energy_candidate_poison_once_for_test(
    int poison_rank, PressureEnergyCandidatePoisonKind kind,
    std::uint8_t corrector) noexcept {
  g_candidate_poison_rank.store(poison_rank, std::memory_order_release);
  g_candidate_poison_kind.store(static_cast<std::uint8_t>(kind),
                                std::memory_order_release);
  g_candidate_poison_corrector.store(corrector, std::memory_order_release);
  g_candidate_poison_armed.store(true, std::memory_order_release);
}

void clear_pressure_energy_candidate_poison_for_test() noexcept {
  g_candidate_poison_armed.store(false, std::memory_order_release);
  g_candidate_poison_rank.store(-1, std::memory_order_release);
  g_candidate_poison_corrector.store(0U, std::memory_order_release);
}

bool product_final_flux_history_diagnostic_for_test(
    ProductFinalFluxHistoryDiagnostic& out) noexcept {
  if (!g_final_flux_history_published.load(std::memory_order_acquire)) {
    out = {};
    return false;
  }
  out = g_final_flux_history_diagnostic;
  return out.valid;
}

bool pressure_correction_warm_start_diagnostic_for_test(
    PressureCorrectionWarmStartDiagnostic& out) noexcept {
  if (!g_pressure_correction_warm_start_published.load(
          std::memory_order_acquire)) {
    out = {};
    return false;
  }
  out = g_pressure_correction_warm_start_diagnostic;
  return out.valid;
}

void suppress_pressure_correction_warm_start_once_for_test() noexcept {
  g_pressure_correction_warm_start_suppressed.store(
      true, std::memory_order_release);
}

void arm_fresh_initialization_candidate_poison_once_for_test(
    int poison_rank, FreshInitializationPoisonKind kind) noexcept {
  g_fresh_initialization_poison_rank.store(poison_rank,
                                           std::memory_order_release);
  g_fresh_initialization_poison_kind.store(
      static_cast<std::uint8_t>(kind), std::memory_order_release);
  g_fresh_initialization_poison_armed.store(true,
                                            std::memory_order_release);
  g_fresh_initialization_published.store(false, std::memory_order_release);
  g_fresh_initialization_diagnostic = {};
}

void clear_fresh_initialization_diagnostic_for_test() noexcept {
  g_fresh_initialization_poison_armed.store(false,
                                            std::memory_order_release);
  g_fresh_initialization_poison_rank.store(-1, std::memory_order_release);
  g_fresh_initialization_published.store(false, std::memory_order_release);
  g_fresh_initialization_diagnostic = {};
  g_cold_velocity_dependents_published.store(false,
                                             std::memory_order_release);
  g_cold_velocity_dependents_diagnostic = {};
}

bool fresh_initialization_diagnostic_for_test(
    FreshInitializationDiagnostic& out) noexcept {
  if (!g_fresh_initialization_published.load(std::memory_order_acquire)) {
    out = {};
    return false;
  }
  out = g_fresh_initialization_diagnostic;
  return out.valid;
}

bool cold_velocity_dependents_diagnostic_for_test(
    ColdVelocityDependentsDiagnostic& out) noexcept {
  if (!g_cold_velocity_dependents_published.load(std::memory_order_acquire)) {
    out = {};
    return false;
  }
  out = g_cold_velocity_dependents_diagnostic;
  return out.valid;
}

void arm_restart_restore_allocation_failure_once_for_test(
    int failing_rank, RestartRestoreAllocationPoint point) noexcept {
  g_restart_restore_allocation_failure_rank.store(failing_rank,
                                                   std::memory_order_release);
  g_restart_restore_allocation_failure_point.store(
      static_cast<std::uint8_t>(point), std::memory_order_release);
  g_restart_restore_allocation_lowest_failing_rank.store(
      -1, std::memory_order_release);
  g_restart_restore_allocation_failure_armed.store(
      true, std::memory_order_release);
}

void clear_restart_restore_allocation_failure_for_test() noexcept {
  g_restart_restore_allocation_failure_armed.store(
      false, std::memory_order_release);
  g_restart_restore_allocation_failure_rank.store(-1,
                                                   std::memory_order_release);
  g_restart_restore_allocation_lowest_failing_rank.store(
      -1, std::memory_order_release);
}

int restart_restore_allocation_lowest_failing_rank_for_test() noexcept {
  return g_restart_restore_allocation_lowest_failing_rank.load(
      std::memory_order_acquire);
}

Status convective_cfl_acceptance_status_for_test(
    TimeControlKind control, double outward_max, double limit) noexcept {
  return convective_cfl_acceptance_status(control, outward_max, limit);
}

}  // namespace detail
#endif

struct CompiledCasePlan::Impl {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  CpuExecutionPlan cpu;
  std::optional<StlScanPlan> scan;
  std::optional<ImmersedSurfacePlan> surface;
  std::optional<EBTopology> topology;
  std::optional<BoundaryStencilPlan> ibm_boundary;
  std::optional<IbmEquationInterfacePlan> ibm_equations;
  std::optional<IbmPhysicalBoundaryFluxAuthority>
      ibm_physical_boundary_flux;
  std::optional<SurfaceQuadraturePlan> quadrature;
  IbmReconstructionAudit ibm_boundary_reconstruction{};
  IbmReconstructionAudit ibm_surface_reconstruction{};
  BoundaryPlan boundary;
  std::array<BoundaryFaceSpec, 6U> boundary_specs{};
  SchemePlan schemes;
  TimeSchemePlan time;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  ClosedMassPlan closed_mass;
  ContributionRegistry contributions;
  TurbulencePlan turbulence;
  EquationPlanSet equations;
  PisoPlan piso;
  PressureVelocityCoupler coupler;
  PressureEnergyCandidateBoundaryFinalizer candidate_boundary_finalizer;
  std::optional<FreshStartKinematicProjectionPlan> fresh_projection;
  ReductionEngine reductions;
  FieldSchema schema;
  ArenaLayout layout;
  StateLayers layers;
  FrozenExecutionGraph graph;
  IoServicePlan io;
  LinearWorkspaceRequirements krylov_requirements{};
  LinearWorkspaceRequirements auxiliary_krylov_requirements{};
  MgWorkspaceRequirements mg_requirements{};
  SolverWorkspace krylov_workspace;
  SolverWorkspace auxiliary_krylov_workspace;
  MgWorkspace mg_workspace;
  std::array<HaloEngine, 6U> stage_halos;
  HaloEngine predictor_donor_halo;
  HaloEngine final_velocity_halo;
  HaloEngine coupled_state_halo;
  HaloEngine candidate_state_halo;
  HaloEngine coupled_thermal_halo;
  HaloEngine candidate_thermal_halo;
  HaloEngine candidate_finalizer_state_halo;
  HaloEngine correction_halo;
  HaloEngine candidate_correction_halo;
  HaloEngine pressure_energy_enthalpy_halo;
  HaloEngine momentum_limiter_halo;
  HaloEngine krylov_halo;
  HaloEngine mg_halo;
  std::vector<HaloEngine> coarse_halos;
  std::vector<HaloEngine*> coarse_halo_pointers;
  std::optional<RemoteDonorExchangePlan> ibm_pressure_donors;
  std::optional<RemoteDonorExchangePlan> ibm_candidate_pressure_donors;
  std::optional<RemoteDonorExchangePlan> ibm_candidate_velocity_donors;
  std::optional<RemoteDonorExchangePlan> ibm_candidate_rate_donors;
  std::optional<RemoteDonorExchangePlan> ibm_gradient_donors;
  std::optional<RemoteDonorExchangePlan> ibm_momentum_donors;
  std::optional<RemoteDonorExchangePlan> ibm_rate_donors;
  std::optional<RemoteDonorExchangePlan> ibm_force_donors;
  FaceFluxStorage final_flux;
  FaceFluxStorage phi_workspace;
  std::vector<double> pressure_face_storage;
  // The nonlinear energy assembly coefficients, frozen target h_f values,
  // and directional delta-h_f values are simultaneously live inside a
  // Schur action and therefore require disjoint cold storage.
  std::vector<double> energy_assembly_face_storage;
  std::vector<double> energy_frozen_enthalpy_face_storage;
  std::vector<double> energy_directional_enthalpy_face_storage;
  // Cold-compiled E_h factors.  The two cell spans hold the immutable local
  // diagonal and the transactional response staging field.  Limiter branch
  // codes are compact metadata; neither allocation participates in a halo.
  std::vector<double> energy_compiled_enthalpy_cell_storage;
  std::vector<std::uint16_t> energy_compiled_enthalpy_branch_storage;
  std::vector<std::uint8_t> pressure_mg_cell_activity;
  std::vector<std::uint8_t> pressure_mg_x_activity;
  std::vector<std::uint8_t> pressure_mg_y_activity;
  std::vector<std::uint8_t> pressure_mg_z_activity;
  PlanFingerprint pressure_mg_activity_fingerprint{};
  PlanFingerprint pressure_mg_activity_collective{};
  std::array<ProductFreezePhase, 10U> phases{};
  PlanSummary summary{};
  ProductFields fields;
  PlanFingerprint fingerprint{};
  PlanFingerprint schema_fingerprint{};
  PlanFingerprint cpu_fingerprint{};
  PlanFingerprint stl_fingerprint{};
  std::uintptr_t state_address{};
};

struct ProductDriver::Impl {
  CompiledCasePlan plan;
  AttemptTransaction transaction;
  FinalFaceFluxAuthority final_flux_authority;
  FinalFaceFluxWriter final_flux_writer;
  PendingFaceFluxView pending_flux;
  std::optional<FinalForceCache> final_force_cache;
  TimeControllerState time;
  PressureLinearOperator pressure_operator;
  std::optional<IbmPressureOperator> ibm_pressure_operator;
  NativeCartesianMgPlan pressure_mg;
  ConservativeEnthalpyEndpoint enthalpy_endpoint;
  MgPlanCounters pressure_mg_counters{};
  ResourceCounters resources{};
  MPI_Comm communicator{MPI_COMM_NULL};
  std::vector<SnapshotFieldView> output_fields;
  std::vector<RestartFieldView> restart_fields;
  std::vector<RestartFieldView> restart_previous_fields;
  std::vector<RestartFieldView> restart_rate_fields;
  std::vector<RestartFieldView> restart_previous_rate_fields;
  std::vector<RestartExpectedField> restart_expected_fields;
  std::vector<RestartExpectedField> restart_expected_rate_fields;
  std::vector<double> species_values;
  std::vector<ConstFieldView> species_accepted;
  std::vector<ConstFieldView> species_previous;
  std::vector<ConstFieldView> passive_accepted;
  std::vector<ConstFieldView> passive_previous;
  std::vector<FieldView> species_trial;
  std::vector<FieldView> pressure_energy_candidate_species;
  std::vector<ConstFieldView> pressure_energy_candidate_species_const;
  std::vector<ConstFieldView>
      pressure_energy_candidate_species_boundary_aliases;
  std::vector<PrimitiveHistory> pressure_energy_candidate_species_history;
  // Owned-cell snapshot of the live same-target energy residual.  Candidate
  // evaluation overwrites the shared residual workspace, so this cold buffer
  // is the allocation-free alpha-zero IBM equivalence oracle.
  std::vector<double> pressure_energy_baseline_energy_residual;
  std::vector<FieldView> passive_trial;
  std::vector<FieldView> species_low;
  std::vector<FieldView> passive_low;
  ThermophysicalGhostHistory enthalpy_ghosts{};
  std::vector<ThermophysicalGhostHistory> species_ghosts;
  std::vector<ThermophysicalGhostHistory> passive_ghosts;
  ThermophysicalGhostAuthority trial_enthalpy_ghost{};
  std::vector<ThermophysicalGhostAuthority> trial_species_ghosts;
  std::vector<ThermophysicalGhostAuthority> trial_passive_ghosts;
  std::vector<PrimitiveHistory> species_history;
  std::vector<PrimitiveHistory> passive_history;
  std::vector<PredictorRateHistory> species_rate_history;
  std::vector<PredictorRateHistory> passive_rate_history;
  std::vector<FieldView> species_rate_output;
  std::vector<FieldView> passive_rate_output;
  std::vector<double> boundary_scalars;
  std::vector<Real3> boundary_vectors;
  std::vector<double> boundary_normal_gradients;
  std::vector<FieldView> halo_views;
  // Cold-sized publication dependencies: five primitive fields, every
  // independent species that determines EOS boundary density/flux, and the
  // pressure-reference authority.
  std::vector<RevisionDependency> final_dependencies;
  // Preallocated [pressure, enthalpy, independent species...] authority used
  // by the allocation-free physical thermophysical ghost closure.
  std::vector<BoundaryGhostFieldAuthority> boundary_thermo_authority_fields;
  BoundaryThermophysicalGhostCertificate boundary_thermo_certificate{};
  double pressure_reference{};
  double previous_pressure_reference{};
  struct PendingPressureReference {
    double value{};
    PressureReferenceCertificate certificate{};
    RevisionToken terminal_audit{};

    bool valid() const noexcept {
      return std::isfinite(value) && value > 0.0 && certificate.valid() &&
             terminal_audit != 0U;
    }
  } pending_pressure_reference{};
  double closed_mass_target{};
  bool initialized{};
  bool pressure_mg_initialized{};
  bool pressure_correction_warm_start_valid{};
  bool pending_force_cache{};
  StageId attempt_stage{};
  BdfCoefficients effective_bdf{};
  std::uint8_t thermophysical_predictor_calls{};
  bool temporal_method_fallback{};
  NumericalFailureContext numerical_failure{};
  ThermophysicalPredictorDiagnostics predictor_diagnostics{};
  MomentumPredictorLimiterReport momentum_predictor_limiter{};
  MomentumPredictorSolveReport momentum_predictor_solve{};
  PressureEnergyGlobalizationAttemptReport pressure_energy_globalization{};
  std::array<DriverStageTiming, kDriverTimedStageCapacity> attempt_timings{};
  std::size_t attempt_timing_count{};
  StageId timed_stage{};
  std::chrono::steady_clock::time_point timed_stage_begin{};

  BoundaryResolvedValues resolved_boundary_values() const noexcept {
    return {{boundary_scalars.data(), boundary_scalars.size()},
            {boundary_vectors.data(), boundary_vectors.size()},
            {boundary_normal_gradients.data(),
             boundary_normal_gradients.size()}};
  }
  Status execute_attempt(const StepTime& step,
                         PisoAttemptReport& report,
                         PreparedAttemptFinish& prepared) noexcept;
  void clear_pending_attempt_side_state() noexcept;
  void discard_pending_attempt_side_state() noexcept;
  void finalize_pending_force_cache() noexcept;
  void commit_pending_attempt_side_state() noexcept;
  DriverResourceReport resource_snapshot() const noexcept;
  void reset_stage_timings(StageId stage) noexcept;
  void begin_timed_stage(StageId stage) noexcept;
  void finish_stage_timings() noexcept;
  void accumulate_stage_timings(DriverStepReport& report) const noexcept;
  Status initialize_common_fields(
      const DriverInitialState& initial, double enthalpy,
      const ThermoState& thermo, const MolecularTransportState& transport,
      double heat_capacity) noexcept;
  Status rebuild_cold_velocity_dependents(double pressure_reference,
                                          bool restart) noexcept;

  ~Impl() noexcept {
    if (communicator != MPI_COMM_NULL) {
      int finalized = 0;
      if (MPI_Finalized(&finalized) == MPI_SUCCESS && finalized == 0)
        MPI_Comm_free(&communicator);
      else
        communicator = MPI_COMM_NULL;
    }
  }
};

void ProductDriver::Impl::clear_pending_attempt_side_state() noexcept {
  pending_pressure_reference = {};
  pending_force_cache = false;
  trial_enthalpy_ghost = {};
  boundary_thermo_certificate = {};
  std::fill(trial_species_ghosts.begin(), trial_species_ghosts.end(),
            ThermophysicalGhostAuthority{});
  std::fill(trial_passive_ghosts.begin(), trial_passive_ghosts.end(),
            ThermophysicalGhostAuthority{});
}

void ProductDriver::Impl::discard_pending_attempt_side_state() noexcept {
  // A failed solve may have overwritten the persistent pressure-correction
  // workspace.  Its warm-start bit is therefore an accepted-state authority,
  // not an attempt-local convergence hint, and must be retired on every
  // rejected transaction.
  pressure_correction_warm_start_valid = false;
  clear_pending_attempt_side_state();
}

void ProductDriver::Impl::finalize_pending_force_cache() noexcept {
  if (!pending_force_cache) return;
  const Status finalized =
      final_force_cache.has_value()
          ? final_force_cache->finalize(transaction)
          : Status{StatusCode::invalid_plan, kProductBinding};
  // A prepared transaction freezes the pending cache identity.  Once its
  // no-fail commit/reject has run, FinalForceCache::finalize is therefore a
  // mechanical publication/discard and cannot legitimately fail.
  assert(static_cast<bool>(finalized));
  static_cast<void>(finalized);
  pending_force_cache = false;
}

void ProductDriver::Impl::commit_pending_attempt_side_state() noexcept {
  // The pressure-reference scalar is part of the same accepted physical
  // state as p, h, rho, T, U, and the final mass flux.  A successful
  // transaction must therefore have carried the terminal audit authority,
  // rather than rotating an unbound side-channel double.
  assert(pending_pressure_reference.valid());
  enthalpy_ghosts.previous = enthalpy_ghosts.accepted;
  enthalpy_ghosts.accepted = trial_enthalpy_ghost;
  for (std::size_t index = 0U; index < species_ghosts.size(); ++index) {
    species_ghosts[index].previous = species_ghosts[index].accepted;
    species_ghosts[index].accepted = trial_species_ghosts[index];
  }
  for (std::size_t index = 0U; index < passive_ghosts.size(); ++index) {
    passive_ghosts[index].previous = passive_ghosts[index].accepted;
    passive_ghosts[index].accepted = trial_passive_ghosts[index];
  }
  previous_pressure_reference = pressure_reference;
  pressure_reference = pending_pressure_reference.value;
  pressure_correction_warm_start_valid = true;
  clear_pending_attempt_side_state();
}

void ProductDriver::Impl::reset_stage_timings(StageId stage) noexcept {
  attempt_timings = {};
  attempt_timing_count = 0U;
  timed_stage = stage;
  timed_stage_begin = std::chrono::steady_clock::now();
}

void ProductDriver::Impl::begin_timed_stage(StageId stage) noexcept {
  const auto now = std::chrono::steady_clock::now();
  if (timed_stage != 0U &&
      attempt_timing_count < attempt_timings.size()) {
    attempt_timings[attempt_timing_count++] = {
        timed_stage,
        static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                now - timed_stage_begin)
                .count())};
  }
  timed_stage = stage;
  timed_stage_begin = now;
}

void ProductDriver::Impl::finish_stage_timings() noexcept {
  if (timed_stage == 0U) return;
  begin_timed_stage(0U);
}

void ProductDriver::Impl::accumulate_stage_timings(
    DriverStepReport& report) const noexcept {
  for (std::size_t source = 0U; source < attempt_timing_count; ++source) {
    const DriverStageTiming timing = attempt_timings[source];
    std::size_t target = 0U;
    while (target < report.stage_timing_count &&
           report.stage_timings[target].stage != timing.stage)
      ++target;
    if (target == report.stage_timing_count) {
      if (target == report.stage_timings.size()) continue;
      report.stage_timings[target].stage = timing.stage;
      ++report.stage_timing_count;
    }
    std::uint64_t& accumulated = report.stage_timings[target].nanoseconds;
    accumulated = timing.nanoseconds > UINT64_MAX - accumulated
                      ? UINT64_MAX
                      : accumulated + timing.nanoseconds;
  }
}

DriverResourceReport ProductDriver::Impl::resource_snapshot() const noexcept {
  DriverResourceReport result;
  const auto add = [](std::uint64_t& target, std::uint64_t value) {
    target = value > UINT64_MAX - target ? UINT64_MAX : target + value;
  };
  const auto halo = [&](const HaloEngine& engine) {
    const HaloRuntimeCounters counters = engine.runtime_counters();
    add(result.structured_exchanges, counters.begin_calls);
    add(result.structured_messages, counters.messages_started);
    add(result.structured_bytes, counters.bytes_packed);
    add(result.structured_bytes, counters.bytes_unpacked);
  };
  const CompiledCasePlan::Impl& product = *plan.implementation_;
  for (const HaloEngine& engine : product.stage_halos) halo(engine);
  halo(product.predictor_donor_halo);
  halo(product.final_velocity_halo);
  halo(product.coupled_state_halo);
  halo(product.candidate_state_halo);
  halo(product.coupled_thermal_halo);
  halo(product.candidate_thermal_halo);
  halo(product.candidate_finalizer_state_halo);
  halo(product.correction_halo);
  halo(product.candidate_correction_halo);
  halo(product.pressure_energy_enthalpy_halo);
  halo(product.momentum_limiter_halo);
  halo(product.krylov_halo);
  halo(product.mg_halo);
  for (const HaloEngine& engine : product.coarse_halos) halo(engine);
  const auto donor = [&](const std::optional<RemoteDonorExchangePlan>& plan) {
    if (!plan.has_value()) return;
    const RemoteDonorExchangeCounters counters = plan->runtime_counters();
    add(result.ibm_exchanges, counters.exchange_calls);
    add(result.ibm_messages, counters.peer_messages);
    add(result.ibm_bytes, counters.bytes);
  };
  donor(product.ibm_pressure_donors);
  donor(product.ibm_candidate_pressure_donors);
  donor(product.ibm_candidate_velocity_donors);
  donor(product.ibm_candidate_rate_donors);
  donor(product.ibm_gradient_donors);
  donor(product.ibm_momentum_donors);
  donor(product.ibm_rate_donors);
  donor(product.ibm_force_donors);
  const LinearReductionCounters reduction = product.reductions.counters();
  result.reduction_collectives = reduction.blocking_operations;
  result.reduction_nanoseconds = reduction.wall_nanoseconds;
  result.reduction_logical_bytes = reduction.logical_bytes;
  result.reduction_tree_messages = reduction.tree_messages;
  result.linear_iterations = resources.linear_iterations;
  const MgPlanCounters mg = pressure_mg.counters();
  add(result.structured_messages, mg.point_to_point_messages);
  add(result.structured_bytes, mg.point_to_point_bytes);
  result.mg_blocking_collectives = mg.blocking_collectives;
  result.mg_collective_logical_bytes = mg.collective_logical_bytes;
  result.exact_numeric_refills = mg.numeric_refreshes;
  result.hierarchy_rebuilds = mg.hierarchy_rebuilds;
  result.preconditioner_applications = mg.applications;
  return result;
}

DriverResourceReport resource_difference(DriverResourceReport after,
                                         DriverResourceReport before) noexcept {
  DriverResourceReport result;
  const auto difference = [](std::uint64_t newer,
                             std::uint64_t older) noexcept {
    return newer >= older ? newer - older : UINT64_MAX;
  };
  result.structured_exchanges = difference(
      after.structured_exchanges, before.structured_exchanges);
  result.structured_messages = difference(
      after.structured_messages, before.structured_messages);
  result.structured_bytes =
      difference(after.structured_bytes, before.structured_bytes);
  result.ibm_exchanges = difference(after.ibm_exchanges, before.ibm_exchanges);
  result.ibm_messages = difference(after.ibm_messages, before.ibm_messages);
  result.ibm_bytes = difference(after.ibm_bytes, before.ibm_bytes);
  result.reduction_collectives = difference(
      after.reduction_collectives, before.reduction_collectives);
  result.reduction_nanoseconds = difference(
      after.reduction_nanoseconds, before.reduction_nanoseconds);
  result.reduction_logical_bytes = difference(
      after.reduction_logical_bytes, before.reduction_logical_bytes);
  result.reduction_tree_messages = difference(
      after.reduction_tree_messages, before.reduction_tree_messages);
  result.mg_blocking_collectives = difference(
      after.mg_blocking_collectives, before.mg_blocking_collectives);
  result.mg_collective_logical_bytes = difference(
      after.mg_collective_logical_bytes,
      before.mg_collective_logical_bytes);
  result.linear_iterations =
      difference(after.linear_iterations, before.linear_iterations);
  result.exact_numeric_refills = difference(
      after.exact_numeric_refills, before.exact_numeric_refills);
  result.hierarchy_rebuilds = difference(
      after.hierarchy_rebuilds, before.hierarchy_rebuilds);
  result.preconditioner_applications = difference(
      after.preconditioner_applications,
      before.preconditioner_applications);
  return result;
}

CompiledCasePlan::~CompiledCasePlan() noexcept { release(); }
CompiledCasePlan::CompiledCasePlan(CompiledCasePlan&& other) noexcept
    : implementation_(std::exchange(other.implementation_, nullptr)) {}
CompiledCasePlan& CompiledCasePlan::operator=(CompiledCasePlan&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = std::exchange(other.implementation_, nullptr);
  }
  return *this;
}
void CompiledCasePlan::release() noexcept {
  delete std::exchange(implementation_, nullptr);
}
PlanFingerprint CompiledCasePlan::fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->fingerprint;
}
PlanSummary CompiledCasePlan::summary() const noexcept {
  return implementation_ == nullptr ? PlanSummary{} : implementation_->summary;
}
Span<const ProductFreezePhase> CompiledCasePlan::freeze_order() const noexcept {
  return implementation_ == nullptr
             ? Span<const ProductFreezePhase>{}
             : Span<const ProductFreezePhase>{implementation_->phases.data(),
                                              implementation_->phases.size()};
}
const FieldSchema* CompiledCasePlan::field_schema() const noexcept {
  return implementation_ == nullptr ? nullptr : &implementation_->schema;
}
const ArenaLayout* CompiledCasePlan::arena_layout() const noexcept {
  return implementation_ == nullptr ? nullptr : &implementation_->layout;
}
const FrozenExecutionGraph* CompiledCasePlan::execution_graph() const noexcept {
  return implementation_ == nullptr ? nullptr : &implementation_->graph;
}
const IoServicePlan* CompiledCasePlan::io_services() const noexcept {
  return implementation_ == nullptr ? nullptr : &implementation_->io;
}
const ContributionPlan* CompiledCasePlan::contribution_plan() const noexcept {
  return implementation_ == nullptr ? nullptr
                                    : implementation_->contributions.plan();
}
const CpuExecutionPlan* CompiledCasePlan::cpu_execution_plan() const noexcept {
  return implementation_ == nullptr ? nullptr : &implementation_->cpu;
}
PlanFingerprint CompiledCasePlan::cpu_plan_fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->cpu_fingerprint;
}
PlanFingerprint CompiledCasePlan::stl_fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->stl_fingerprint;
}
std::uintptr_t CompiledCasePlan::state_storage_address() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->state_address;
}
std::uintptr_t CompiledCasePlan::krylov_storage_address() const noexcept {
  return implementation_ == nullptr
             ? 0U
             : implementation_->krylov_workspace.vector_storage_address();
}
std::uintptr_t CompiledCasePlan::mg_storage_address() const noexcept {
  return implementation_ == nullptr
             ? 0U
             : implementation_->mg_workspace.storage_address();
}

Status ProductCompiler::compile(MPI_Comm communicator,
                                const ValidatedModel& model,
                                const std::filesystem::path& case_root,
                                CompiledCasePlan& out) noexcept try {
  if (communicator == MPI_COMM_NULL || model.fingerprint == 0U ||
      out.implementation_ != nullptr) {
    return {StatusCode::invalid_plan, kProductInput};
  }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  g_candidate_storage_published.store(false, std::memory_order_release);
  g_candidate_storage_diagnostic = {};
#endif
  std::unique_ptr<CompiledCasePlan::Impl> candidate{
      new (std::nothrow) CompiledCasePlan::Impl};
  if (!candidate) {
    return {StatusCode::allocation_failure, kProductAllocation};
  }
  candidate->boundary_specs = model.boundaries;
  Status status = CartesianGeometryCompiler::compile(
      communicator, model.mesh, {}, candidate->geometry, candidate->patch);
  CpuExecutionRequest cpu_request;
  cpu_request.threads_per_rank = 1U;
  cpu_request.pure_mpi = true;
  cpu_request.bind_threads = false;
  if (status)
    status = CpuExecutionPlan::compile(communicator, cpu_request,
                                       candidate->cpu);
  if (status)
    candidate->cpu_fingerprint = candidate->cpu.semantic_fingerprint();
  const bool immersed = model.immersed_boundary.has_value();
  if (status && immersed) {
    candidate->scan.emplace();
    StlScanBudget scan_budget{model.mesh.limits.max_memory_bytes_per_rank / 2U,
                              model.mesh.limits.max_memory_bytes_per_rank,
                              UINT64_C(16000000), UINT64_C(100000), 1U};
    status = StlScanCompiler::compile(
        communicator, case_root, model.immersed_boundary->stl_file,
        candidate->geometry, candidate->patch, CartesianAxis::y, scan_budget,
        *candidate->scan);
    if (status) {
      candidate->surface.emplace();
      status = ImmersedSurfaceCompiler::compile(*candidate->scan,
                                                *candidate->surface);
      if (status)
        candidate->stl_fingerprint =
            candidate->surface->source_triangle_fingerprint();
    }
    if (status) {
      candidate->topology.emplace();
      ImmersedPlanLimits limits;
      status = EBTopologyCompiler::compile(
          communicator, candidate->geometry, candidate->patch,
          *candidate->scan, *candidate->surface,
          model.immersed_boundary->fluid_side, limits, *candidate->topology);
    }
    if (status) {
      status = make_pressure_mg_activity(
          *candidate->topology, candidate->patch.cells,
          candidate->pressure_mg_cell_activity,
          candidate->pressure_mg_x_activity,
          candidate->pressure_mg_y_activity,
          candidate->pressure_mg_z_activity,
          candidate->pressure_mg_activity_fingerprint,
          candidate->pressure_mg_activity_collective);
    }
    if (status) {
      ImmersedPlanLimits limits;
      limits.stencil.policy =
          model.immersed_boundary->reconstruction_policy;
      if (limits.stencil.policy ==
          IbmReconstructionPolicy::adaptive_order) {
        limits.stencil.standard_reach = 2U;
      }
      ImmersedDomainBoundaryPolicy boundary_policy;
      for (std::size_t face = 0U; face < model.boundaries.size(); ++face) {
        const BoundaryKind kind = model.boundaries[face].flow_kind;
        boundary_policy.allow_one_sided_quadratic[face] =
            kind == BoundaryKind::symmetry || kind == BoundaryKind::slip;
      }
      // Case validation already requires paired periodic faces.  Re-derive
      // the authority here from both faces so an IBM policy can never acquire
      // a one-sided or caller-invented periodic image permission.
      for (int axis = 0; axis < 3 && status; ++axis) {
        const std::size_t lower = static_cast<std::size_t>(2 * axis);
        const std::size_t upper = lower + 1U;
        const bool lower_periodic =
            model.boundaries[lower].flow_kind == BoundaryKind::periodic;
        const bool upper_periodic =
            model.boundaries[upper].flow_kind == BoundaryKind::periodic;
        if (lower_periodic != upper_periodic) {
          status = {StatusCode::invalid_plan, kProductInput};
          break;
        }
        boundary_policy.allow_periodic_images[lower] =
            lower_periodic && upper_periodic;
        boundary_policy.allow_periodic_images[upper] =
            lower_periodic && upper_periodic;
      }
      candidate->ibm_boundary.emplace();
      candidate->quadrature.emplace();
      if (status)
        status = BoundaryStencilCompiler::compile(
            communicator, candidate->geometry, candidate->patch,
            *candidate->surface, *candidate->topology, boundary_policy, limits,
            *candidate->ibm_boundary);
      if (status)
        status = SurfaceQuadratureCompiler::compile(
            communicator, candidate->geometry, candidate->patch,
            *candidate->surface, *candidate->topology, boundary_policy, limits,
            *candidate->quadrature);
      if (status)
        status = reduce_ibm_reconstruction_audit(
            communicator, candidate->ibm_boundary->reconstruction().audit(),
            candidate->ibm_boundary_reconstruction);
      if (status)
        status = reduce_ibm_reconstruction_audit(
            communicator, candidate->quadrature->reconstruction().audit(),
            candidate->ibm_surface_reconstruction);
    }
    if (status) {
      // Scan triangles and line intersections are cold construction
      // workspace.  Surface/topology plans are the sealed runtime authority.
      candidate->scan.reset();
    }
  }
  if (!status) return status;
  candidate->phases[0U] = ProductFreezePhase::geometry_and_decomposition;
  std::size_t local_cells = 0U;
  if (!detail::product_cell_count(candidate->patch.cells, local_cells)) {
    return {StatusCode::invalid_plan, kProductInput};
  }
  std::size_t maximum_cells_per_rank = 0U;
  if (!detail::product_cell_count(candidate->geometry.global_cells(),
                                  maximum_cells_per_rank)) {
    return {StatusCode::invalid_plan, kProductInput};
  }
  PisoPlanSpec piso_spec;
  piso_spec.coupling = model.solver.coupling;
  piso_spec.pressure_correctors = 2U;
  piso_spec.pressure_stage = 50U;
  piso_spec.final_flux_slot = 0U;
  piso_spec.pressure_algorithm = model.solver.pressure.algorithm;
  piso_spec.mg_correction_scaling = model.solver.pressure.mg_correction_scaling;
  piso_spec.pressure_solve = {
      model.solver.pressure.absolute_tolerance,
      model.solver.pressure.relative_tolerance,
      model.solver.pressure.maximum_iterations,
      model.solver.pressure.true_residual_interval,
      model.solver.pressure.krylov_restart};
  piso_spec.eos_tolerance = model.solver.terminal.eos;
  piso_spec.continuity_tolerance = model.solver.terminal.continuity;
  // Energy is normalized by the target-time temporal energy scale in the
  // product driver, so the existing dimensionless terminal conservation
  // tolerance is the compatible authority until a separate case key is
  // introduced.
  piso_spec.energy_tolerance = model.solver.terminal.continuity;
  piso_spec.closed_mass_tolerance = model.solver.terminal.closed_mass;
  piso_spec.gauge_tolerance = model.solver.terminal.gauge;
  constexpr std::uint8_t krylov_ghosts = 1U;
  status = make_linear_workspace_requirements(
      piso_spec.pressure_algorithm, candidate->patch.cells, krylov_ghosts,
      piso_spec.pressure_solve.restart, ReductionMode::mpi_allreduce, 1U,
      candidate->krylov_requirements);
  if (status) {
    if (piso_spec.pressure_algorithm == LinearAlgorithm::fgmres) {
      // Preserve the accepted default identity and layout.  Momentum and the
      // enthalpy endpoint have always shared this storage sequentially with
      // pressure and use at most twelve configured FGMRES directions.
      candidate->auxiliary_krylov_requirements =
          candidate->krylov_requirements;
    } else {
      // Pressure may select a non-FGMRES algorithm, while momentum and the
      // conservative enthalpy endpoint remain FGMRES/Jacobi solves.  Bind an
      // auxiliary logical workspace to the same preallocated storage; these
      // consumers execute in disjoint stages and add no hot allocation.
      status = make_linear_workspace_requirements(
          LinearAlgorithm::fgmres, candidate->patch.cells, krylov_ghosts,
          kAuxiliaryFgmresMaximumRestart, ReductionMode::mpi_allreduce, 1U,
          candidate->auxiliary_krylov_requirements);
    }
  }
  NativeCartesianMgSpec mg_spec;
  mg_spec.policy = detail::production_pressure_mg_policy();
  if (status)
    status = make_mg_workspace_requirements(
        communicator, candidate->geometry, candidate->patch, mg_spec.policy,
        1U, candidate->mg_requirements);
  FieldRegistry registry;
  const std::uint8_t ghosts = immersed ? 4U : 2U;
  if (status)
    status = register_fields(registry, candidate->fields, model, ghosts,
                             candidate->krylov_requirements,
                             candidate->auxiliary_krylov_requirements,
                             candidate->mg_requirements);
  if (status && immersed) {
    const std::array<RemoteDonorFieldSpec, 1U> pressure_fields{{
        {candidate->fields.pressure_correction, 1U}}};
    candidate->ibm_pressure_donors.emplace();
    status = RemoteDonorExchangePlan::analyze(
        communicator, candidate->geometry.global_cells(), candidate->patch,
        candidate->ibm_boundary->reconstruction(),
        {pressure_fields.data(), pressure_fields.size()},
        kIbmPressureCorrectionDonorStage,
        *candidate->ibm_pressure_donors);
    const std::array<RemoteDonorFieldSpec, 1U> candidate_pressure_fields{{
        {candidate->fields.pressure_energy_candidate_pressure_correction,
         1U}}};
    if (status) {
      candidate->ibm_candidate_pressure_donors.emplace();
      status = RemoteDonorExchangePlan::analyze(
          communicator, candidate->geometry.global_cells(), candidate->patch,
          candidate->ibm_boundary->reconstruction(),
          {candidate_pressure_fields.data(), candidate_pressure_fields.size()},
          kIbmCandidatePressureCorrectionDonorStage,
          *candidate->ibm_candidate_pressure_donors);
    }
    const std::array<RemoteDonorFieldSpec, 1U> candidate_velocity_fields{{
        {candidate->fields.pressure_energy_candidate_velocity, 3U}}};
    if (status) {
      candidate->ibm_candidate_velocity_donors.emplace();
      status = RemoteDonorExchangePlan::analyze(
          communicator, candidate->geometry.global_cells(), candidate->patch,
          candidate->ibm_boundary->reconstruction(),
          {candidate_velocity_fields.data(), candidate_velocity_fields.size()},
          kIbmCandidateVelocityDonorStage,
          *candidate->ibm_candidate_velocity_donors);
    }
    const std::array<RemoteDonorFieldSpec, 2U> candidate_rate_fields{{
        {candidate->fields.pressure_energy_candidate_pressure, 1U},
        {candidate->fields.pressure_energy_candidate_temperature, 1U}}};
    if (status) {
      candidate->ibm_candidate_rate_donors.emplace();
      status = RemoteDonorExchangePlan::analyze(
          communicator, candidate->geometry.global_cells(), candidate->patch,
          candidate->ibm_boundary->reconstruction(),
          {candidate_rate_fields.data(), candidate_rate_fields.size()},
          kIbmCandidateEnergyRateDonorStage,
          *candidate->ibm_candidate_rate_donors);
    }
    const std::array<RemoteDonorFieldSpec, 1U> gradient_fields{{
        {candidate->fields.velocity, 3U}}};
    if (status) {
      candidate->ibm_gradient_donors.emplace();
      status = RemoteDonorExchangePlan::analyze(
          communicator, candidate->geometry.global_cells(), candidate->patch,
          candidate->ibm_boundary->reconstruction(),
          {gradient_fields.data(), gradient_fields.size()},
          kIbmGradientDonorStage, *candidate->ibm_gradient_donors);
    }
    const std::array<RemoteDonorFieldSpec, 2U> momentum_fields{{
        {candidate->fields.effective_viscosity, 1U},
        {candidate->fields.pressure, 1U}}};
    if (status) {
      candidate->ibm_momentum_donors.emplace();
      status = RemoteDonorExchangePlan::analyze(
          communicator, candidate->geometry.global_cells(), candidate->patch,
          candidate->ibm_boundary->reconstruction(),
          {momentum_fields.data(), momentum_fields.size()}, 130U,
          *candidate->ibm_momentum_donors);
    }
    std::vector<RemoteDonorFieldSpec> rate_fields;
    if (status) {
      rate_fields.reserve(2U + candidate->fields.scalars.size());
      rate_fields.push_back({candidate->fields.pressure, 1U});
      rate_fields.push_back({candidate->fields.temperature, 1U});
      for (FieldId scalar : candidate->fields.scalars)
        rate_fields.push_back({scalar, 1U});
      candidate->ibm_rate_donors.emplace();
      status = RemoteDonorExchangePlan::analyze(
          communicator, candidate->geometry.global_cells(), candidate->patch,
          candidate->ibm_boundary->reconstruction(),
          {rate_fields.data(), rate_fields.size()}, 161U,
          *candidate->ibm_rate_donors);
    }
    const std::array<RemoteDonorFieldSpec, 4U> force_fields{{
        {candidate->fields.velocity, 3U},
        {candidate->fields.pressure, 1U},
        {candidate->fields.velocity_gradient, 9U},
        {candidate->fields.effective_viscosity, 1U}}};
    if (status) {
      candidate->ibm_force_donors.emplace();
      status = RemoteDonorExchangePlan::analyze(
          communicator, candidate->geometry.global_cells(), candidate->patch,
          candidate->quadrature->reconstruction(),
          {force_fields.data(), force_fields.size()}, 160U,
          *candidate->ibm_force_donors);
    }
  }
  if (!status) return status;
  candidate->phases[1U] = ProductFreezePhase::capability_registration;
  candidate->schema_fingerprint = registry.fingerprint();
  status = compile_graph(
      candidate->fields, ghosts, candidate->patch.cells, local_cells,
      candidate->ibm_pressure_donors.has_value()
          ? candidate->ibm_pressure_donors->stats()
          : RemoteDonorExchangeStats{},
      candidate->ibm_candidate_pressure_donors.has_value()
          ? candidate->ibm_candidate_pressure_donors->stats()
          : RemoteDonorExchangeStats{},
      candidate->ibm_candidate_velocity_donors.has_value()
          ? candidate->ibm_candidate_velocity_donors->stats()
          : RemoteDonorExchangeStats{},
      candidate->ibm_candidate_rate_donors.has_value()
          ? candidate->ibm_candidate_rate_donors->stats()
          : RemoteDonorExchangeStats{},
      candidate->ibm_gradient_donors.has_value()
          ? candidate->ibm_gradient_donors->stats()
          : RemoteDonorExchangeStats{},
      candidate->ibm_momentum_donors.has_value()
          ? candidate->ibm_momentum_donors->stats()
          : RemoteDonorExchangeStats{},
      candidate->ibm_rate_donors.has_value()
          ? candidate->ibm_rate_donors->stats()
          : RemoteDonorExchangeStats{},
      candidate->ibm_force_donors.has_value()
          ? candidate->ibm_force_donors->stats()
          : RemoteDonorExchangeStats{},
      candidate->graph);
  std::vector<SnapshotFieldSpec> snapshots{
      {candidate->fields.velocity, 3U}, {candidate->fields.pressure, 1U},
      {candidate->fields.enthalpy, 1U}};
  for (FieldId scalar : candidate->fields.scalars) snapshots.push_back({scalar, 1U});
  std::size_t snapshot_components = 0U;
  std::size_t snapshot_bytes = 0U;
  std::size_t coordinate_values = 0U;
  std::size_t coordinate_bytes = 0U;
  std::size_t staging_bytes = 0U;
  if (!detail::product_checked_add(5U, candidate->fields.scalars.size(),
                                   snapshot_components) ||
      !detail::product_field_bytes(local_cells, snapshot_components,
                                   snapshot_bytes) ||
      !detail::product_checked_add(
          static_cast<std::size_t>(candidate->patch.cells.x),
          static_cast<std::size_t>(candidate->patch.cells.y),
          coordinate_values) ||
      !detail::product_checked_add(
          coordinate_values,
          static_cast<std::size_t>(candidate->patch.cells.z),
          coordinate_values) ||
      !detail::product_checked_add(coordinate_values, 3U,
                                   coordinate_values) ||
      !detail::product_checked_multiply(coordinate_values, sizeof(double),
                                        coordinate_bytes) ||
      !detail::product_checked_add(snapshot_bytes, coordinate_bytes,
                                   staging_bytes) ||
      !detail::product_checked_add(staging_bytes, 65536U, staging_bytes)) {
    return {StatusCode::invalid_plan, kProductAnalysis};
  }
  std::array<RuntimeServiceCapacity, 5U> services{};
  for (std::size_t index = 0U; index < services.size(); ++index) {
    services[index] = {static_cast<RuntimeServiceKind>(index),
                       static_cast<StageId>(200U + index), snapshot_bytes,
                       staging_bytes, 8U};
  }
  if (status)
    status = IoServicePlan::compile({snapshots.data(), snapshots.size()},
                                    {services.data(), services.size()},
                                    local_cells, candidate->io);
  std::uint32_t registered_capabilities =
      detail::product_fields | detail::product_structured_ghosts |
      detail::product_cache_dependencies | detail::product_exact_numeric |
      detail::product_coarse_numeric | detail::product_preconditioner |
      detail::product_workspace | detail::product_service_snapshot |
      detail::product_collective_epochs;
  if (immersed) registered_capabilities |= detail::product_ibm_donors;
  if (status)
    status = detail::validate_product_capabilities(registered_capabilities,
                                                   immersed);
  if (!status) return status;
  candidate->phases[2U] = ProductFreezePhase::logical_analysis;
  FieldRegistry schema_registry = registry;
  status = schema_registry.freeze(candidate->schema);
  std::vector<ArenaFieldRequest> requests;
  try {
    requests.reserve(candidate->schema.size());
    for (const FieldDescriptor& descriptor : candidate->schema) {
      Int3 interior = candidate->patch.cells;
      if (descriptor.id == candidate->fields.krylov_scalars) {
        const std::size_t scalar_doubles = std::max(
            candidate->krylov_requirements.scalar_doubles,
            candidate->auxiliary_krylov_requirements.scalar_doubles);
        if (scalar_doubles >
            static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
          status = {StatusCode::invalid_plan, kProductCapacity};
          break;
        }
        interior = {static_cast<std::int32_t>(scalar_doubles), 1, 1};
      } else if (descriptor.id == candidate->fields.mg_arena) {
        interior = candidate->mg_requirements.arena_shape;
      }
      requests.push_back({descriptor.id, interior, {0U},
                          state_field(descriptor.id, candidate->fields)
                              ? FieldLifetime::state_layer
                              : FieldLifetime::persistent_workspace});
    }
  } catch (const std::bad_alloc&) {
    status = {StatusCode::allocation_failure, kProductAllocation};
  }
  if (status)
    status = ArenaLayout::compile(candidate->schema,
                                  {requests.data(), requests.size()},
                                  candidate->layout);
  if (status) status = StateLayers::allocate(candidate->layout, candidate->layers);
  if (!status) {
    return status;
  }
  candidate->phases[3U] = ProductFreezePhase::schema_and_allocation;
  const PlanFingerprint registry_before = registry.fingerprint();
  status = BoundaryCompiler::compile(
      communicator, model, candidate->geometry, candidate->patch, registry,
      candidate->boundary, candidate->schemes, candidate->time);
  if (status && registry.fingerprint() != registry_before)
    status = {StatusCode::invalid_plan, kProductInstantiation};
  if (status)
    status = ThermodynamicsPlan::compile(model.thermophysics,
                                         {model.transported_scalars.data(),
                                          model.transported_scalars.size()},
                                         candidate->thermodynamics);
  if (status)
    status = TransportPlan::compile(model.thermophysics,
                                    candidate->thermodynamics,
                                    candidate->transport);
  if (status)
    status = ClosedMassPlan::compile(model.pressure_reference,
                                     model.thermophysics,
                                     candidate->closed_mass);
  std::vector<FieldId> declared;
  if (status) {
    declared.reserve(candidate->schema.size());
    for (const FieldDescriptor& descriptor : candidate->schema)
      declared.push_back(descriptor.id);
    status = candidate->contributions.configure(
        {declared.data(), declared.size()});
  }
  TurbulencePlanSpec turbulence_spec;
  turbulence_spec.kind = model.turbulence;
  if (status)
    status = TurbulencePlan::compile(
        communicator, turbulence_spec, candidate->geometry, candidate->patch,
        candidate->fields.effective_viscosity, 20U,
        candidate->contributions, candidate->turbulence);
  if (status) status = candidate->contributions.freeze();
  std::vector<ScalarEquationSpec> scalar_specs;
  if (status) {
    for (std::size_t index = 0U; index < candidate->fields.scalars.size();
         ++index) {
      scalar_specs.push_back({candidate->fields.scalars[index],
                              model.transported_scalars[index].role,
                              model.transported_scalars[index].molecular_schmidt,
                              model.transported_scalars[index].turbulent_schmidt});
    }
    EquationPlanSpec equation_spec;
    equation_spec.density = candidate->fields.rho;
    equation_spec.velocity = candidate->fields.velocity;
    equation_spec.pressure_perturbation = candidate->fields.pressure;
    equation_spec.enthalpy = candidate->fields.enthalpy;
    equation_spec.temperature = candidate->fields.temperature;
    equation_spec.effective_viscosity = candidate->fields.effective_viscosity;
    equation_spec.velocity_gradient = candidate->fields.velocity_gradient;
    equation_spec.pressure_compressibility = candidate->fields.compressibility;
    equation_spec.pressure_reference = model.pressure_reference;
    equation_spec.closed_mass_service_stage =
        model.pressure_reference == PressureReferenceKind::closed_mass ? 11U : 0U;
    equation_spec.scalars = {scalar_specs.data(), scalar_specs.size()};
    // This value participates in the semantic equation fingerprint, so use a
    // rank-invariant conservative ceiling. It is a validation bound only; all
    // numeric storage below remains sized to the actual local patch.
    equation_spec.maximum_cells_per_rank = maximum_cells_per_rank;
    status = EquationPlanSet::compile(
        communicator, candidate->schemes, candidate->geometry,
        candidate->patch, candidate->boundary, candidate->contributions,
        candidate->thermodynamics, candidate->transport, equation_spec,
        candidate->equations);
  }
  if (status && immersed) {
    candidate->ibm_equations.emplace();
    status = IbmEquationInterfacePlan::compile(
        candidate->equations.kernels(), *candidate->topology,
        *candidate->ibm_boundary, *candidate->ibm_equations);
    if (status) {
      candidate->ibm_physical_boundary_flux.emplace();
      status = IbmPhysicalBoundaryFluxAuthority::compile(
          communicator, candidate->geometry, candidate->patch,
          *candidate->topology, *candidate->ibm_equations,
          *candidate->ibm_physical_boundary_flux);
    }
  }
  if (status)
    status = PisoPlan::compile(communicator, candidate->equations, piso_spec,
                               candidate->piso);
  if (!status) {
    return status;
  }
  candidate->phases[4U] = ProductFreezePhase::plan_instantiation;
  std::size_t face_doubles = 0U;
  std::size_t compiled_enthalpy_cell_doubles = 0U;
  if (!detail::product_face_doubles(candidate->patch.cells, face_doubles)) {
    return {StatusCode::invalid_plan, kProductCapacity};
  }
  if (!detail::product_checked_multiply(
          local_cells, 2U, compiled_enthalpy_cell_doubles)) {
    return {StatusCode::invalid_plan, kProductCapacity};
  }
  try {
    candidate->pressure_face_storage.assign(face_doubles, 0.0);
    candidate->energy_assembly_face_storage.assign(face_doubles, 0.0);
    candidate->energy_frozen_enthalpy_face_storage.assign(face_doubles, 0.0);
    candidate->energy_directional_enthalpy_face_storage.assign(face_doubles,
                                                               0.0);
    candidate->energy_compiled_enthalpy_cell_storage.assign(
        compiled_enthalpy_cell_doubles, 0.0);
    candidate->energy_compiled_enthalpy_branch_storage.assign(face_doubles,
                                                              0U);
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kProductCapacity};
  }
  status = FaceFluxStorage::allocate_final(candidate->patch.cells,
                                           candidate->final_flux);
  if (status)
    // Four pressure/candidate fluxes can be live together.  Two additional
    // replicas retain all three momentum AFC reconstructions while their
    // common face alpha is formed, avoiding a second reconstruction pass.
    status = FaceFluxStorage::allocate_workspace(candidate->patch.cells, 6U,
                                                 candidate->phi_workspace);
  if (!status) {
    return status;
  }
  candidate->phases[5U] = ProductFreezePhase::numeric_capacity;
  constexpr std::array<StageId, 6U> kHaloStages{
      {10U, 15U, 20U, 30U, 40U, 60U}};
  for (std::size_t index = 0U; index < kHaloStages.size() && status; ++index) {
    status = reserve_stage_halo(
        communicator, candidate->patch, candidate->boundary.halo_topology(),
        candidate->schema, candidate->graph, kHaloStages[index],
        candidate->stage_halos[index]);
  }
  const std::array<HaloFieldSpec, 1U> predictor_donor_halo{{
      {candidate->fields.rho, 1U, 1U}}};
  if (status)
    status = candidate->predictor_donor_halo.reserve(
        communicator, candidate->patch,
        {predictor_donor_halo.data(), predictor_donor_halo.size()},
        candidate->boundary.halo_topology());
  // Stage 20 now carries U and the effective-method momentum history.  The
  // terminal rate refresh still needs only U, so retain its original narrow
  // communication payload in a separately frozen plan.
  const std::array<HaloFieldSpec, 1U> final_velocity_halo{{
      {candidate->fields.velocity, ghosts, 3U}}};
  if (status)
    status = candidate->final_velocity_halo.reserve(
        communicator, candidate->patch,
        {final_velocity_halo.data(), final_velocity_halo.size()},
        candidate->boundary.halo_topology());
  // The same-target pressure-energy loop first publishes primitive authority.
  // Effective thermal material depends on the resulting velocity ghosts and
  // therefore owns a later, narrow exchange below.
  const std::array<HaloFieldSpec, 5U> coupled_state_halo{{
      {candidate->fields.velocity, ghosts, 3U},
      {candidate->fields.rho, 1U, 1U},
      {candidate->fields.pressure, ghosts, 1U},
      {candidate->fields.enthalpy, ghosts, 1U},
      {candidate->fields.temperature, ghosts, 1U},
  }};
  if (status)
    status = candidate->coupled_state_halo.reserve(
        communicator, candidate->patch,
        {coupled_state_halo.data(), coupled_state_halo.size()},
        candidate->boundary.halo_topology());
  // Candidate residual evaluation owns an independent communication lineage.
  // It must never issue candidate field revisions through the live trial
  // coupled-state halo, even though both routes have the same geometric reach.
  std::vector<HaloFieldSpec> candidate_state_halo;
  if (status) {
    candidate_state_halo.reserve(
        5U + candidate->fields.pressure_energy_candidate_species.size());
    candidate_state_halo.push_back(
        {candidate->fields.pressure_energy_candidate_velocity, ghosts, 3U});
    candidate_state_halo.push_back(
        {candidate->fields.pressure_energy_candidate_density, 1U, 1U});
    candidate_state_halo.push_back(
        {candidate->fields.pressure_energy_candidate_pressure, ghosts, 1U});
    candidate_state_halo.push_back(
        {candidate->fields.pressure_energy_candidate_enthalpy, ghosts, 1U});
    candidate_state_halo.push_back(
        {candidate->fields.pressure_energy_candidate_temperature, ghosts, 1U});
    for (FieldId species :
         candidate->fields.pressure_energy_candidate_species)
      candidate_state_halo.push_back({species, ghosts, 1U});
    status = candidate->candidate_state_halo.reserve(
        communicator, candidate->patch,
        {candidate_state_halo.data(), candidate_state_halo.size()},
        candidate->boundary.halo_topology());
  }
  const std::array<HaloFieldSpec, 2U> coupled_thermal_halo{{
      {candidate->fields.thermal_conductivity, 1U, 1U},
      {candidate->fields.enthalpy_diffusivity, 1U, 1U},
  }};
  if (status)
    status = candidate->coupled_thermal_halo.reserve(
        communicator, candidate->patch,
        {coupled_thermal_halo.data(), coupled_thermal_halo.size()},
        candidate->boundary.halo_topology());
  const std::array<HaloFieldSpec, 2U> candidate_thermal_halo{{
      {candidate->fields.pressure_energy_candidate_thermal_conductivity, 1U,
       1U},
      {candidate->fields.pressure_energy_candidate_enthalpy_diffusivity, 1U,
       1U},
  }};
  if (status)
    status = candidate->candidate_thermal_halo.reserve(
        communicator, candidate->patch,
        {candidate_thermal_halo.data(), candidate_thermal_halo.size()},
        candidate->boundary.halo_topology());
  // Physical boundary finalization consumes an exact, narrow candidate
  // primitive contract.  It must not accept the residual-state halo as
  // a superset because that would erase U/Y lineage and permit stale ghosts.
  std::vector<HaloFieldSpec> candidate_finalizer_state_halo;
  if (status) {
    candidate_finalizer_state_halo.reserve(
        5U + candidate->fields.pressure_energy_candidate_species.size());
    candidate_finalizer_state_halo.push_back(
        {candidate->fields.pressure_energy_candidate_pressure, 1U, 1U});
    candidate_finalizer_state_halo.push_back(
        {candidate->fields.pressure_energy_candidate_enthalpy, 1U, 1U});
    candidate_finalizer_state_halo.push_back(
        {candidate->fields.pressure_energy_candidate_density, 1U, 1U});
    candidate_finalizer_state_halo.push_back(
        {candidate->fields.pressure_energy_candidate_temperature, 1U, 1U});
    candidate_finalizer_state_halo.push_back(
        {candidate->fields.pressure_energy_candidate_velocity, 1U, 3U});
    for (FieldId species :
         candidate->fields.pressure_energy_candidate_species)
      candidate_finalizer_state_halo.push_back({species, 1U, 1U});
    status = candidate->candidate_finalizer_state_halo.reserve(
        communicator, candidate->patch,
        {candidate_finalizer_state_halo.data(),
         candidate_finalizer_state_halo.size()},
        candidate->boundary.halo_topology());
  }
  const std::array<HaloFieldSpec, 1U> correction_halo{{
      {candidate->fields.pressure_correction, 1U, 1U}}};
  if (status)
    status = candidate->correction_halo.reserve(
        communicator, candidate->patch,
        {correction_halo.data(), correction_halo.size()},
        candidate->boundary.halo_topology());
  const std::array<HaloFieldSpec, 1U> candidate_correction_halo{{
      {candidate->fields.pressure_energy_candidate_pressure_correction, 1U,
       1U}}};
  if (status)
    status = candidate->candidate_correction_halo.reserve(
        communicator, candidate->patch,
        {candidate_correction_halo.data(), candidate_correction_halo.size()},
        candidate->boundary.halo_topology());
  const std::array<HaloFieldSpec, 2U> pressure_energy_enthalpy_halo{{
      {candidate->fields.schur_eliminated_enthalpy, 2U, 1U},
      {candidate->fields.pressure_energy_delta_temperature, 1U, 1U},
  }};
  if (status)
    status = candidate->pressure_energy_enthalpy_halo.reserve(
        communicator, candidate->patch,
        {pressure_energy_enthalpy_halo.data(),
         pressure_energy_enthalpy_halo.size()},
        candidate->boundary.halo_topology());
  const std::array<HaloFieldSpec, 1U> momentum_limiter_halo{{
      {candidate->fields.h_by_a, 1U, 3U}}};
  if (status)
    status = candidate->momentum_limiter_halo.reserve(
        communicator, candidate->patch,
        {momentum_limiter_halo.data(), momentum_limiter_halo.size()},
        candidate->boundary.halo_topology());
  const std::array<HaloFieldSpec, 1U> krylov_halo{{
      {candidate->fields.krylov_vectors, 1U, 1U}}};
  if (status)
    status = candidate->krylov_halo.reserve(
        communicator, candidate->patch,
        {krylov_halo.data(), krylov_halo.size()},
        candidate->boundary.halo_topology());
  const std::array<HaloFieldSpec, 1U> mg_halo{{
      {candidate->fields.mg_arena, 1U, 1U}}};
  if (status)
    status = candidate->mg_halo.reserve(
        communicator, candidate->mg_requirements.levels[0U].patch,
        {mg_halo.data(), mg_halo.size()},
        candidate->boundary.halo_topology());
  if (status) {
    candidate->coarse_halos.resize(candidate->mg_requirements.level_count - 1U);
    candidate->coarse_halo_pointers.resize(candidate->coarse_halos.size());
    for (std::size_t level = 1U;
         level < candidate->mg_requirements.level_count; ++level) {
      status = candidate->coarse_halos[level - 1U].reserve(
          communicator, candidate->mg_requirements.levels[level].patch,
          {mg_halo.data(), mg_halo.size()},
          candidate->boundary.halo_topology());
      if (!status) break;
      candidate->coarse_halo_pointers[level - 1U] =
          &candidate->coarse_halos[level - 1U];
    }
  }
  if (status && candidate->ibm_pressure_donors.has_value())
    status = candidate->ibm_pressure_donors->bind(communicator);
  if (status && candidate->ibm_candidate_pressure_donors.has_value())
    status = candidate->ibm_candidate_pressure_donors->bind(communicator);
  if (status && candidate->ibm_candidate_velocity_donors.has_value())
    status = candidate->ibm_candidate_velocity_donors->bind(communicator);
  if (status && candidate->ibm_candidate_rate_donors.has_value())
    status = candidate->ibm_candidate_rate_donors->bind(communicator);
  if (status && candidate->ibm_gradient_donors.has_value())
    status = candidate->ibm_gradient_donors->bind(communicator);
  if (status && candidate->ibm_momentum_donors.has_value())
    status = candidate->ibm_momentum_donors->bind(communicator);
  if (status && candidate->ibm_rate_donors.has_value())
    status = candidate->ibm_rate_donors->bind(communicator);
  if (status && candidate->ibm_force_donors.has_value())
    status = candidate->ibm_force_donors->bind(communicator);
  if (!status) {
    return status;
  }
  candidate->phases[6U] = ProductFreezePhase::communication_binding;
  FieldView krylov_vectors;
  FieldView krylov_scalars;
  FieldView mg_arena;
  status = candidate->layers.runtime_view(
      FieldLifetime::persistent_workspace, candidate->fields.krylov_vectors,
      krylov_vectors);
  if (status)
    status = candidate->layers.runtime_view(
        FieldLifetime::persistent_workspace, candidate->fields.krylov_scalars,
        krylov_scalars);
  if (status)
    status = candidate->layers.runtime_view(
        FieldLifetime::persistent_workspace, candidate->fields.mg_arena,
        mg_arena);
  if (status)
    status = SolverWorkspace::bind(candidate->krylov_requirements,
                                   krylov_vectors, krylov_scalars,
                                   candidate->krylov_workspace);
  if (status && piso_spec.pressure_algorithm != LinearAlgorithm::fgmres)
    status = SolverWorkspace::bind(candidate->auxiliary_krylov_requirements,
                                   krylov_vectors, krylov_scalars,
                                   candidate->auxiliary_krylov_workspace);
  if (status)
    status = MgWorkspace::bind(candidate->mg_requirements, mg_arena,
                               candidate->mg_workspace);
  FieldView r_au;
  FieldView h_by_a;
  FieldView pressure_gradient;
  if (status)
    status = candidate->layers.runtime_view(
        FieldLifetime::persistent_workspace, candidate->fields.r_au, r_au);
  if (status)
    status = candidate->layers.runtime_view(
        FieldLifetime::persistent_workspace, candidate->fields.h_by_a,
        h_by_a);
  if (status)
    status = candidate->layers.runtime_view(
        FieldLifetime::persistent_workspace,
        candidate->fields.pressure_gradient, pressure_gradient);
  FaceFieldView x_pressure;
  FaceFieldView y_pressure;
  FaceFieldView z_pressure;
  if (status)
    status = make_pressure_face_views(candidate->pressure_face_storage,
                                      candidate->patch.cells, x_pressure,
                                      y_pressure, z_pressure);
  FaceFieldView x_energy_assembly;
  FaceFieldView y_energy_assembly;
  FaceFieldView z_energy_assembly;
  FaceFieldView x_energy_frozen;
  FaceFieldView y_energy_frozen;
  FaceFieldView z_energy_frozen;
  FaceFieldView x_energy_directional;
  FaceFieldView y_energy_directional;
  FaceFieldView z_energy_directional;
  if (status)
    status = make_pressure_face_views(
        candidate->energy_assembly_face_storage, candidate->patch.cells,
        x_energy_assembly, y_energy_assembly, z_energy_assembly);
  if (status)
    status = make_pressure_face_views(
        candidate->energy_frozen_enthalpy_face_storage,
        candidate->patch.cells, x_energy_frozen, y_energy_frozen,
        z_energy_frozen);
  if (status)
    status = make_pressure_face_views(
        candidate->energy_directional_enthalpy_face_storage,
        candidate->patch.cells, x_energy_directional, y_energy_directional,
        z_energy_directional);
  if (status) {
    const std::array<StorageIdentity, 4U> identities{{
        x_pressure.storage_identity, x_energy_assembly.storage_identity,
        x_energy_frozen.storage_identity,
        x_energy_directional.storage_identity}};
    const bool coherent_axes =
        x_pressure.storage_identity == y_pressure.storage_identity &&
        x_pressure.storage_identity == z_pressure.storage_identity &&
        x_energy_assembly.storage_identity ==
            y_energy_assembly.storage_identity &&
        x_energy_assembly.storage_identity ==
            z_energy_assembly.storage_identity &&
        x_energy_frozen.storage_identity == y_energy_frozen.storage_identity &&
        x_energy_frozen.storage_identity == z_energy_frozen.storage_identity &&
        x_energy_directional.storage_identity ==
            y_energy_directional.storage_identity &&
        x_energy_directional.storage_identity ==
            z_energy_directional.storage_identity;
    bool independent = coherent_axes;
    for (std::size_t left = 0U; left < identities.size(); ++left)
      for (std::size_t right = left + 1U; right < identities.size(); ++right)
        independent = independent && identities[left] != identities[right];
    if (!independent)
      status = {StatusCode::invalid_plan, kProductBinding};
  }
  FaceFluxView phi_h_by_a;
  if (status)
    status = candidate->phi_workspace.workspace_view(0U, 1U, phi_h_by_a);
  FaceFluxView provisional_flux_capacity;
  FaceFluxView candidate_flux_capacity;
  FaceFluxView candidate_final_flux_capacity;
  if (status)
    status = candidate->phi_workspace.workspace_view(
        1U, 1U, provisional_flux_capacity);
  if (status)
    status = candidate->phi_workspace.workspace_view(
        2U, 1U, candidate_flux_capacity);
  if (status)
    status = candidate->phi_workspace.workspace_view(
        3U, 1U, candidate_final_flux_capacity);
  if (status) {
    const bool same_storage =
        phi_h_by_a.x.storage_identity ==
            provisional_flux_capacity.x.storage_identity &&
        phi_h_by_a.x.storage_identity ==
            candidate_flux_capacity.x.storage_identity &&
        phi_h_by_a.x.storage_identity ==
            candidate_final_flux_capacity.x.storage_identity;
    const bool disjoint_replicas =
        phi_h_by_a.x.base != provisional_flux_capacity.x.base &&
        phi_h_by_a.x.base != candidate_flux_capacity.x.base &&
        provisional_flux_capacity.x.base != candidate_flux_capacity.x.base &&
        phi_h_by_a.x.base != candidate_final_flux_capacity.x.base &&
        provisional_flux_capacity.x.base !=
            candidate_final_flux_capacity.x.base &&
        candidate_flux_capacity.x.base != candidate_final_flux_capacity.x.base;
    if (!same_storage || !disjoint_replicas)
      status = {StatusCode::invalid_plan, kProductBinding};
  }
  if (status)
    status = ReductionEngine::compile(
        communicator, ReductionMode::mpi_allreduce,
        std::max<std::size_t>(
            std::max<std::size_t>(
                candidate->krylov_requirements.reduction_capacity,
                candidate->auxiliary_krylov_requirements.reduction_capacity),
            8U),
        candidate->reductions);
  // Freeze the incompatible-cold-start projection only after every borrowed
  // address (halos, reductions, linear/MG workspaces, cell fields, and face
  // replicas) has reached its final location.  CompiledCasePlan::Impl itself
  // is heap-owned, so these addresses remain stable when the public plan is
  // moved into ProductDriver.
  if (status) {
    FieldView fresh_chi;
    FieldView fresh_rhs;
    FieldView fresh_diagonal;
    FieldView fresh_candidate_velocity;
    status = candidate->layers.runtime_view(
        FieldLifetime::persistent_workspace,
        candidate->fields.pressure_correction, fresh_chi);
    if (status)
      status = candidate->layers.runtime_view(
          FieldLifetime::persistent_workspace, candidate->fields.pressure_rhs,
          fresh_rhs);
    if (status)
      status = candidate->layers.runtime_view(
          FieldLifetime::persistent_workspace,
          candidate->fields.pressure_diagonal, fresh_diagonal);
    if (status)
      status = candidate->layers.runtime_view(
          FieldLifetime::persistent_workspace,
          candidate->fields.pressure_energy_candidate_velocity,
          fresh_candidate_velocity);
    if (status) {
      const MgDomainActivityView activity =
          candidate->topology.has_value()
              ? MgDomainActivityView{
                    {candidate->pressure_mg_cell_activity.data(),
                     candidate->pressure_mg_cell_activity.size()},
                    {candidate->pressure_mg_x_activity.data(),
                     candidate->pressure_mg_x_activity.size()},
                    {candidate->pressure_mg_y_activity.data(),
                     candidate->pressure_mg_y_activity.size()},
                    {candidate->pressure_mg_z_activity.data(),
                     candidate->pressure_mg_z_activity.size()},
                    candidate->pressure_mg_activity_fingerprint,
                    candidate->pressure_mg_activity_collective}
              : MgDomainActivityView{};
      SolverWorkspace* fresh_solver_workspace =
          piso_spec.pressure_algorithm == LinearAlgorithm::fgmres
              ? &candidate->krylov_workspace
              : &candidate->auxiliary_krylov_workspace;
      const std::uint32_t fresh_restart =
          piso_spec.pressure_algorithm == LinearAlgorithm::fgmres
              ? piso_spec.pressure_solve.restart
              : kAuxiliaryFgmresMaximumRestart;
      const bool fresh_capacity_valid =
          piso_spec.pressure_solve.maximum_iterations <=
          std::numeric_limits<std::uint32_t>::max() - fresh_restart;
      if (!fresh_capacity_valid)
        status = {StatusCode::invalid_plan, kProductCapacity};
      // The cold-start projection solve gets exactly one Krylov restart of
      // headroom; the regular pressure-solve budget remains the case authority.
      const std::uint32_t fresh_maximum_iterations =
          fresh_capacity_valid
              ? piso_spec.pressure_solve.maximum_iterations + fresh_restart
              : 0U;
      FreshStartKinematicProjectionSpec fresh_spec;
      fresh_spec.communicator = communicator;
      fresh_spec.geometry = &candidate->geometry;
      fresh_spec.kernels = &candidate->equations.kernels();
      fresh_spec.boundary = &candidate->boundary;
      fresh_spec.patch = candidate->patch;
      fresh_spec.activity = activity;
      fresh_spec.route =
          FreshStartProjectionLinearRoute::native_mg_fgmres;
      fresh_spec.solve = {1.0e-13,
                          1.0e-13,
                          fresh_maximum_iterations,
                          piso_spec.pressure_solve.true_residual_interval,
                          fresh_restart};
      fresh_spec.mg_policy = mg_spec.policy;
      fresh_spec.compatibility_absolute_tolerance = 1.0e-13;
      fresh_spec.compatibility_relative_tolerance = 1.0e-12;
      // Krylov convergence uses the stricter 1e-13/1e-13 gate above, while
      // compatibility retains its separately certified 1e-13/1e-12 gate.
      // Preserve the Fresh module's independently certified physical
      // continuity gate; making it tighter here rejects the exact anchored
      // residual at the one-cell gauge row without improving the solved
      // operator.
      fresh_spec.continuity_absolute_tolerance = 1.0e-12;
      fresh_spec.continuity_relative_tolerance = 1.0e-10;
      const FreshStartKinematicProjectionServices fresh_services{
          &candidate->krylov_halo,
          kFreshProjectionOperatorStage,
          candidate->fields.krylov_vectors,
          &candidate->correction_halo,
          kFreshProjectionCorrectionStage,
          fresh_solver_workspace,
          &candidate->reductions,
          {&candidate->mg_halo,
           &candidate->reductions,
           &candidate->mg_workspace,
           {candidate->coarse_halo_pointers.data(),
            candidate->coarse_halo_pointers.size()}}};
      const FreshStartKinematicProjectionWorkspace fresh_workspace{
          fresh_chi,
          fresh_rhs,
          fresh_diagonal,
          x_energy_assembly,
          y_energy_assembly,
          z_energy_assembly,
          x_pressure,
          y_pressure,
          z_pressure,
          fresh_candidate_velocity,
          candidate_flux_capacity};
      if (status) {
        candidate->fresh_projection.emplace();
        status = FreshStartKinematicProjectionPlan::compile(
            fresh_spec, fresh_services, fresh_workspace,
            *candidate->fresh_projection);
        if (!status) candidate->fresh_projection.reset();
      }
    }
  }
  if (status) {
    const PisoCouplerWorkspace coupler_workspace{
        r_au, h_by_a, pressure_gradient, x_pressure, y_pressure, z_pressure,
        phi_h_by_a};
    PisoCouplerServices coupler_services{
        communicator,
        &candidate->geometry,
        candidate->patch,
        &candidate->boundary,
        &candidate->thermodynamics,
        &candidate->stage_halos[4U],
        40U,
        candidate->fields.rho,
        &candidate->correction_halo,
        41U,
        candidate->fields.pressure_correction,
        candidate->topology.has_value()
            ? PressureContinuityActivityView{
                  {candidate->pressure_mg_cell_activity.data(),
                   candidate->pressure_mg_cell_activity.size()},
                  {candidate->pressure_mg_x_activity.data(),
                   candidate->pressure_mg_x_activity.size()},
                  {candidate->pressure_mg_y_activity.data(),
                   candidate->pressure_mg_y_activity.size()},
                  {candidate->pressure_mg_z_activity.data(),
                   candidate->pressure_mg_z_activity.size()},
                  candidate->pressure_mg_activity_fingerprint,
                  candidate->pressure_mg_activity_collective}
            : PressureContinuityActivityView{}};
    coupler_services.pressure_correction_donors =
        candidate->ibm_pressure_donors.has_value()
            ? &*candidate->ibm_pressure_donors
            : nullptr;
    coupler_services.pressure_correction_donor_stage =
        candidate->ibm_pressure_donors.has_value()
            ? kIbmPressureCorrectionDonorStage
            : 0U;
    coupler_services.candidate_pressure_correction_donors =
        candidate->ibm_candidate_pressure_donors.has_value()
            ? &*candidate->ibm_candidate_pressure_donors
            : nullptr;
    coupler_services.candidate_pressure_correction_donor_stage =
        candidate->ibm_candidate_pressure_donors.has_value()
            ? kIbmCandidatePressureCorrectionDonorStage
            : 0U;
    coupler_services.candidate_pressure_correction_field =
        candidate->ibm_candidate_pressure_donors.has_value()
            ? candidate->fields.pressure_energy_candidate_pressure_correction
            : 0U;
    coupler_services.candidate_pressure_correction_donor_reach =
        candidate->ibm_candidate_pressure_donors.has_value()
            ? candidate->ibm_candidate_pressure_donors->reach()
            : 0U;
    coupler_services.candidate_pressure_correction_donor_fingerprint =
        candidate->ibm_candidate_pressure_donors.has_value()
            ? candidate->ibm_candidate_pressure_donors->fingerprint()
            : 0U;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    int compiler_rank = -1;
    if (MPI_Comm_rank(communicator, &compiler_rank) == MPI_SUCCESS &&
        consume_product_piso_bind_failure(compiler_rank))
      ++coupler_services.candidate_pressure_correction_donor_reach;
#endif
    coupler_services.immersed_interface =
        candidate->ibm_equations.has_value() ? &*candidate->ibm_equations
                                             : nullptr;
    status = PressureVelocityCoupler::bind(
        candidate->piso, candidate->equations, coupler_services,
        coupler_workspace, candidate->coupler);
    // PressureVelocityCoupler::bind is rank-local after its services have
    // been sealed.  A localized IBM can therefore expose a bad local
    // capability on only one rank.  Every rank that entered this block must
    // globalize the result before any successful rank enters the collective
    // candidate-finalizer bind below.
    status = candidate->reductions.consensus(status);
  }
  if (status && product_candidate_boundary_supported(candidate->boundary)) {
    const PressureEnergyCandidateBoundaryFinalizerBinding finalizer_binding{
        communicator,
        &candidate->geometry,
        candidate->patch,
        &candidate->boundary,
        &candidate->equations.kernels(),
        &candidate->thermodynamics,
        &candidate->transport,
        &candidate->coupler,
        candidate->ibm_equations.has_value() ? &*candidate->ibm_equations
                                             : nullptr,
        candidate->ibm_physical_boundary_flux.has_value()
            ? &*candidate->ibm_physical_boundary_flux
            : nullptr,
        candidate->ibm_candidate_pressure_donors.has_value()
            ? &*candidate->ibm_candidate_pressure_donors
            : nullptr,
        candidate->ibm_candidate_pressure_donors.has_value()
            ? kIbmCandidatePressureCorrectionDonorStage
            : StageId{0U},
        candidate->ibm_candidate_pressure_donors.has_value()
            ? candidate->fields.pressure_energy_candidate_pressure_correction
            : FieldId{0U},
        candidate->ibm_candidate_pressure_donors.has_value()
            ? candidate->ibm_candidate_pressure_donors->reach()
            : std::uint8_t{0U}};
    status = PressureEnergyCandidateBoundaryFinalizer::bind(
        finalizer_binding, candidate->candidate_boundary_finalizer);
  }
  FieldView state_probe;
  if (status)
    status = candidate->layers.view(StateRole::trial, candidate->fields.rho,
                                    state_probe);
  if (!status) {
    // Preserve the failing subsystem detail.  Replacing it with the generic
    // product-binding tag made sealed-product diagnostics unable to identify
    // which view/coupler contract was violated.
    return status;
  }
  candidate->state_address =
      reinterpret_cast<std::uintptr_t>(state_probe.base);
  candidate->phases[7U] = ProductFreezePhase::view_and_graph_binding;
  const FaceFluxStorageCounters workspace_flux_capacity =
      candidate->phi_workspace.counters();
  const FaceFluxStorageCounters final_flux_capacity =
      candidate->final_flux.counters();
  const std::array<const HaloEngine*, 7U> independent_halos{{
      &candidate->coupled_state_halo,
      &candidate->candidate_state_halo,
      &candidate->coupled_thermal_halo,
      &candidate->candidate_thermal_halo,
      &candidate->candidate_finalizer_state_halo,
      &candidate->correction_halo,
      &candidate->candidate_correction_halo,
  }};
  bool independent_halo_lineage = true;
  for (std::size_t index = 0U;
       index < independent_halos.size() && independent_halo_lineage; ++index) {
    const HaloEngine& halo = *independent_halos[index];
    independent_halo_lineage = halo.ready() && halo.instance_identity() != 0U;
    for (std::size_t other = index + 1U;
         other < independent_halos.size() && independent_halo_lineage; ++other)
      independent_halo_lineage =
          halo.instance_identity() !=
          independent_halos[other]->instance_identity();
  }
  const bool independent_candidate_donor_lineage =
      !immersed ||
      (candidate->ibm_candidate_pressure_donors.has_value() &&
       candidate->ibm_candidate_pressure_donors->ready() &&
       candidate->ibm_candidate_pressure_donors->fingerprint() != 0U &&
       candidate->ibm_candidate_velocity_donors.has_value() &&
       candidate->ibm_candidate_velocity_donors->ready() &&
       candidate->ibm_candidate_velocity_donors->fingerprint() != 0U &&
       candidate->ibm_candidate_rate_donors.has_value() &&
       candidate->ibm_candidate_rate_donors->ready() &&
       candidate->ibm_candidate_rate_donors->fingerprint() != 0U &&
       candidate->ibm_gradient_donors.has_value() &&
       candidate->ibm_rate_donors.has_value() &&
       candidate->ibm_pressure_donors.has_value() &&
       candidate->ibm_candidate_pressure_donors->fingerprint() !=
           candidate->ibm_pressure_donors->fingerprint() &&
       candidate->ibm_candidate_velocity_donors->fingerprint() !=
           candidate->ibm_gradient_donors->fingerprint() &&
       candidate->ibm_candidate_rate_donors->fingerprint() !=
           candidate->ibm_rate_donors->fingerprint());
  if (!independent_halo_lineage || !independent_candidate_donor_lineage ||
      workspace_flux_capacity.aligned_payload_allocations != 1U ||
      workspace_flux_capacity.replicas != 6U ||
      workspace_flux_capacity.directional_blocks != 18U ||
      final_flux_capacity.aligned_payload_allocations != 1U ||
      final_flux_capacity.replicas != 3U ||
      final_flux_capacity.directional_blocks != 9U)
    return {StatusCode::invalid_plan, kProductBinding};
  const std::array<FieldId, 14U> candidate_lineage_fields{{
      candidate->fields.pressure_energy_candidate_pressure,
      candidate->fields.pressure_energy_candidate_enthalpy,
      candidate->fields.pressure_energy_candidate_density,
      candidate->fields.pressure_energy_candidate_temperature,
      candidate->fields.pressure_energy_candidate_velocity,
      candidate->fields.pressure_energy_candidate_pressure_correction,
      candidate->fields.pressure_energy_candidate_molecular_viscosity,
      candidate->fields.pressure_energy_candidate_effective_viscosity,
      candidate->fields.pressure_energy_candidate_velocity_gradient,
      candidate->fields.pressure_energy_candidate_compressibility,
      candidate->fields.pressure_energy_candidate_enthalpy_compressibility,
      candidate->fields.pressure_energy_candidate_thermal_conductivity,
      candidate->fields.pressure_energy_candidate_heat_capacity,
      candidate->fields.pressure_energy_candidate_enthalpy_diffusivity,
  }};
  std::uint64_t candidate_lineage = kFnvOffset;
  for (FieldId field : candidate_lineage_fields) {
    const FieldDescriptor* descriptor = find_descriptor(candidate->schema, field);
    if (descriptor == nullptr)
      return {StatusCode::invalid_plan, kProductBinding};
    candidate_lineage = detail::product_mix(candidate_lineage, field);
    candidate_lineage =
        detail::product_mix(candidate_lineage, descriptor->components);
    candidate_lineage =
        detail::product_mix(candidate_lineage, descriptor->ghost_width);
  }
  // Only the four pressure/candidate flux roles belong to checkpoint plan
  // semantics.  The two extra replicas are transient AFC cache capacity and
  // must not invalidate an otherwise compatible restart.
  candidate_lineage = detail::product_mix(candidate_lineage, 4U);
  candidate_lineage = detail::product_mix(
      candidate_lineage,
      7U + candidate->fields.pressure_energy_candidate_species.size());
  candidate_lineage = detail::product_mix(candidate_lineage, 1U);
  for (FieldId species :
       candidate->fields.pressure_energy_candidate_species)
    candidate_lineage = detail::product_mix(candidate_lineage, species);
  if (immersed) {
    // Remote-donor fingerprints include rank-local receive/supply layouts.
    // Candidate publication certificates bind those exact local plans, while
    // the product semantic lineage must remain rank invariant.
    candidate_lineage = detail::product_mix(
        candidate_lineage, kIbmCandidatePressureCorrectionDonorStage);
    candidate_lineage = detail::product_mix(
        candidate_lineage, kIbmCandidateVelocityDonorStage);
    candidate_lineage = detail::product_mix(
        candidate_lineage, kIbmCandidateEnergyRateDonorStage);
  }
  if (candidate_lineage == 0U) candidate_lineage = 1U;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  detail::PressureEnergyCandidateStorageDiagnostic storage_diagnostic;
  const std::array<std::pair<FieldId, FieldId>,
                   detail::kPressureEnergyCandidateFieldCount>
      primitive_pairs{{
          {candidate->fields.pressure_energy_candidate_pressure,
           candidate->fields.pressure},
          {candidate->fields.pressure_energy_candidate_enthalpy,
           candidate->fields.enthalpy},
          {candidate->fields.pressure_energy_candidate_density,
           candidate->fields.rho},
          {candidate->fields.pressure_energy_candidate_temperature,
           candidate->fields.temperature},
          {candidate->fields.pressure_energy_candidate_velocity,
           candidate->fields.velocity},
      }};
  for (std::size_t index = 0U;
       index < primitive_pairs.size() && status; ++index)
    status = capture_candidate_field_storage(
        candidate->layers, candidate->layout, primitive_pairs[index].first,
        primitive_pairs[index].second, true,
        storage_diagnostic.fields[index]);
  const std::array<std::pair<FieldId, FieldId>,
                   detail::kPressureEnergyCandidateScratchFieldCount>
      scratch_pairs{{
          {candidate->fields.pressure_energy_candidate_pressure_correction,
           candidate->fields.pressure_correction},
          {candidate->fields.pressure_energy_candidate_molecular_viscosity,
           candidate->fields.molecular_viscosity},
          {candidate->fields.pressure_energy_candidate_effective_viscosity,
           candidate->fields.effective_viscosity},
          {candidate->fields.pressure_energy_candidate_velocity_gradient,
           candidate->fields.velocity_gradient},
          {candidate->fields.pressure_energy_candidate_compressibility,
           candidate->fields.compressibility},
          {candidate->fields.pressure_energy_candidate_enthalpy_compressibility,
           candidate->fields.enthalpy_compressibility},
          {candidate->fields.pressure_energy_candidate_thermal_conductivity,
           candidate->fields.thermal_conductivity},
          {candidate->fields.pressure_energy_candidate_heat_capacity,
           candidate->fields.heat_capacity},
          {candidate->fields.pressure_energy_candidate_enthalpy_diffusivity,
           candidate->fields.enthalpy_diffusivity},
      }};
  for (std::size_t index = 0U; index < scratch_pairs.size() && status;
       ++index)
    status = capture_candidate_field_storage(
        candidate->layers, candidate->layout, scratch_pairs[index].first,
        scratch_pairs[index].second, false,
        storage_diagnostic.scratch_fields[index]);
  for (std::size_t left = 0U;
       left < storage_diagnostic.fields.size() && status; ++left) {
    for (std::size_t right = left + 1U;
         right < storage_diagnostic.fields.size(); ++right)
      if (!candidate_ranges_disjoint(storage_diagnostic.fields[left],
                                     storage_diagnostic.fields[right])) {
        status = {StatusCode::invalid_plan, kProductBinding};
        break;
      }
    for (const auto& scratch : storage_diagnostic.scratch_fields)
      if (!candidate_ranges_disjoint(storage_diagnostic.fields[left],
                                     scratch)) {
        status = {StatusCode::invalid_plan, kProductBinding};
        break;
      }
  }
  for (std::size_t left = 0U;
       left < storage_diagnostic.scratch_fields.size() && status; ++left)
    for (std::size_t right = left + 1U;
         right < storage_diagnostic.scratch_fields.size(); ++right)
      if (!candidate_ranges_disjoint(
              storage_diagnostic.scratch_fields[left],
              storage_diagnostic.scratch_fields[right])) {
        status = {StatusCode::invalid_plan, kProductBinding};
        break;
      }
  std::array<FaceFluxView,
             detail::kPressureEnergyCandidateFluxReplicaCount>
      flux_replicas{};
  for (std::size_t replica = 0U;
       replica < flux_replicas.size() && status; ++replica)
    status = candidate->phi_workspace.workspace_view(
        replica, static_cast<RevisionToken>(replica + 1U),
        flux_replicas[replica]);
  std::uintptr_t flux_stride_bytes = 0U;
  if (status) {
    const std::uintptr_t first =
        reinterpret_cast<std::uintptr_t>(flux_replicas[0U].x.base);
    const std::uintptr_t second =
        reinterpret_cast<std::uintptr_t>(flux_replicas[1U].x.base);
    const std::uintptr_t third =
        reinterpret_cast<std::uintptr_t>(flux_replicas[2U].x.base);
    const std::uintptr_t fourth =
        reinterpret_cast<std::uintptr_t>(flux_replicas[3U].x.base);
    if (first == 0U || second <= first || third <= second ||
        fourth <= third || second - first != third - second ||
        second - first != fourth - third)
      status = {StatusCode::invalid_plan, kProductBinding};
    else
      flux_stride_bytes = second - first;
  }
  for (std::size_t replica = 0U;
       replica < flux_replicas.size() && status; ++replica) {
    const FaceFluxView& flux = flux_replicas[replica];
    auto& captured = storage_diagnostic.flux_replicas[replica];
    captured.bases = {{
        reinterpret_cast<std::uintptr_t>(flux.x.base),
        reinterpret_cast<std::uintptr_t>(flux.y.base),
        reinterpret_cast<std::uintptr_t>(flux.z.base),
    }};
    captured.replica_begin = captured.bases[0U];
    if (captured.replica_begin >
        std::numeric_limits<std::uintptr_t>::max() - flux_stride_bytes) {
      status = {StatusCode::invalid_plan, kProductBinding};
      break;
    }
    captured.replica_end = captured.replica_begin + flux_stride_bytes;
    captured.storage = flux.x.storage_identity;
    captured.revision_domain = flux.x.revision_domain;
    captured.revision = flux.revision;
    if (captured.storage == 0U || captured.revision_domain == 0U ||
        captured.revision == 0U ||
        flux.y.storage_identity != captured.storage ||
        flux.z.storage_identity != captured.storage ||
        flux.y.revision_domain != captured.revision_domain ||
        flux.z.revision_domain != captured.revision_domain ||
        captured.bases[1U] < captured.replica_begin ||
        captured.bases[1U] >= captured.replica_end ||
        captured.bases[2U] < captured.replica_begin ||
        captured.bases[2U] >= captured.replica_end)
      status = {StatusCode::invalid_plan, kProductBinding};
  }
  if (!status) return status;
  storage_diagnostic.lineage_fingerprint = candidate_lineage;
  storage_diagnostic.coupled_state_halo =
      candidate->coupled_state_halo.instance_identity();
  storage_diagnostic.candidate_state_halo =
      candidate->candidate_state_halo.instance_identity();
  storage_diagnostic.coupled_thermal_halo =
      candidate->coupled_thermal_halo.instance_identity();
  storage_diagnostic.candidate_thermal_halo =
      candidate->candidate_thermal_halo.instance_identity();
  storage_diagnostic.candidate_finalizer_state_halo =
      candidate->candidate_finalizer_state_halo.instance_identity();
  storage_diagnostic.correction_halo =
      candidate->correction_halo.instance_identity();
  storage_diagnostic.candidate_correction_halo =
      candidate->candidate_correction_halo.instance_identity();
  storage_diagnostic.coupled_state_halo_plan =
      candidate->coupled_state_halo.plan_stats();
  storage_diagnostic.candidate_state_halo_plan =
      candidate->candidate_state_halo.plan_stats();
  storage_diagnostic.coupled_thermal_halo_plan =
      candidate->coupled_thermal_halo.plan_stats();
  storage_diagnostic.candidate_thermal_halo_plan =
      candidate->candidate_thermal_halo.plan_stats();
  storage_diagnostic.candidate_finalizer_state_halo_plan =
      candidate->candidate_finalizer_state_halo.plan_stats();
  storage_diagnostic.correction_halo_plan =
      candidate->correction_halo.plan_stats();
  storage_diagnostic.candidate_correction_halo_plan =
      candidate->candidate_correction_halo.plan_stats();
  storage_diagnostic.local_ibm_reconstruction_reach =
      candidate->ibm_boundary.has_value()
          ? candidate->ibm_boundary->maximum_halo_reach()
          : std::uint8_t{0U};
  if (candidate->ibm_candidate_pressure_donors.has_value()) {
    storage_diagnostic.candidate_pressure_donor_plan =
        candidate->ibm_candidate_pressure_donors->fingerprint();
    storage_diagnostic.candidate_pressure_donor_stats =
        candidate->ibm_candidate_pressure_donors->stats();
    storage_diagnostic.candidate_pressure_donor_reach =
        candidate->ibm_candidate_pressure_donors->reach();
  }
  if (candidate->ibm_candidate_velocity_donors.has_value()) {
    storage_diagnostic.candidate_velocity_donor_plan =
        candidate->ibm_candidate_velocity_donors->fingerprint();
    storage_diagnostic.candidate_velocity_donor_stats =
        candidate->ibm_candidate_velocity_donors->stats();
  }
  if (candidate->ibm_candidate_rate_donors.has_value()) {
    storage_diagnostic.candidate_rate_donor_plan =
        candidate->ibm_candidate_rate_donors->fingerprint();
    storage_diagnostic.candidate_rate_donor_stats =
        candidate->ibm_candidate_rate_donors->stats();
  }
  storage_diagnostic.workspace_flux_capacity = workspace_flux_capacity;
  storage_diagnostic.final_flux_capacity = final_flux_capacity;
  storage_diagnostic.field_storage_capacity = candidate->layers.counters();
  storage_diagnostic.arena_doubles = candidate->layout.total_doubles();
  storage_diagnostic.execution_graph = candidate->graph.fingerprint();
  const FrozenStage* const corrector_one_stage =
      candidate->graph.stage(40U);
  const FrozenStage* const corrector_two_stage =
      candidate->graph.stage(50U);
  if (corrector_one_stage == nullptr || corrector_two_stage == nullptr)
    return {StatusCode::invalid_plan, kProductBinding};
  storage_diagnostic.corrector_one_resource_contract =
      corrector_one_stage->resources;
  storage_diagnostic.corrector_two_resource_contract =
      corrector_two_stage->resources;
#endif
  std::uint64_t semantic = kFnvOffset;
  semantic = detail::product_mix(semantic, model.fingerprint);
  semantic = detail::product_mix(semantic, candidate->geometry.fingerprint());
  semantic = detail::product_mix(semantic, candidate->cpu_fingerprint);
  semantic = detail::product_mix(semantic,
                                 candidate->boundary.semantic_fingerprint());
  semantic = detail::product_mix(semantic,
                                 candidate->equations.semantic_fingerprint());
  semantic = detail::product_mix(semantic, candidate->piso.fingerprint());
  semantic = detail::product_mix(semantic, candidate->turbulence.fingerprint());
  semantic = detail::product_mix(semantic, candidate->schema_fingerprint);
  semantic = detail::product_mix(semantic, candidate_lineage);
  semantic = detail::product_mix(semantic, immersed ? 1U : 0U);
  if (immersed) {
    semantic = detail::product_mix(
        semantic, candidate->quadrature->physical_fingerprint());
  }
  candidate->fingerprint = semantic == 0U ? 1U : semantic;
  status = collective_semantic(communicator, candidate->fingerprint);
  const bool auxiliary_workspace_valid =
      piso_spec.pressure_algorithm == LinearAlgorithm::fgmres ||
      candidate->auxiliary_krylov_workspace.vector_storage_address() != 0U;
  if (!status || candidate->graph.fingerprint() == 0U ||
      candidate->io.fingerprint() == 0U ||
      candidate->cpu_fingerprint == 0U ||
      candidate->state_address == 0U ||
      candidate->krylov_workspace.vector_storage_address() == 0U ||
      !auxiliary_workspace_valid ||
      candidate->mg_workspace.storage_address() == 0U) {
    return status ? Status{StatusCode::invalid_plan, kProductBinding} : status;
  }
  candidate->phases[8U] = ProductFreezePhase::validation;
  candidate->phases[9U] = ProductFreezePhase::sealed;
  candidate->summary.coupling = candidate->piso.coupling();
  candidate->summary.global_cells = candidate->geometry.global_cells();
  candidate->summary.local_cells = candidate->patch.cells;
  candidate->summary.field_count = candidate->schema.size();
  candidate->summary.arena_doubles = candidate->layout.total_doubles();
  candidate->summary.graph_stage_count = candidate->graph.stages().size;
  candidate->summary.graph_node_count = candidate->graph.nodes().size;
  candidate->summary.maximum_workspace_bytes =
      candidate->graph.resources().max_live_workspace_bytes;
  candidate->summary.service_staging_bytes =
      candidate->io.maximum_staging_bytes();
  candidate->summary.pressure_correctors =
      candidate->piso.pressure_correctors();
  candidate->summary.pressure_absolute_tolerance =
      model.solver.pressure.absolute_tolerance;
  candidate->summary.pressure_relative_tolerance =
      model.solver.pressure.relative_tolerance;
  candidate->summary.pressure_maximum_iterations =
      model.solver.pressure.maximum_iterations;
  candidate->summary.pressure_true_residual_interval =
      model.solver.pressure.true_residual_interval;
  candidate->summary.pressure_krylov_restart =
      model.solver.pressure.krylov_restart;
  candidate->summary.terminal_eos_tolerance = model.solver.terminal.eos;
  candidate->summary.terminal_continuity_tolerance =
      model.solver.terminal.continuity;
  candidate->summary.terminal_closed_mass_tolerance =
      model.solver.terminal.closed_mass;
  candidate->summary.terminal_gauge_tolerance = model.solver.terminal.gauge;
  candidate->summary.immersed = immersed;
  candidate->summary.ibm_boundary_reconstruction =
      candidate->ibm_boundary_reconstruction;
  candidate->summary.ibm_surface_reconstruction =
      candidate->ibm_surface_reconstruction;
  candidate->summary.exact_numeric_certified = false;
  candidate->summary.preconditioner_setup_certified = false;
  candidate->summary.sealed = true;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  storage_diagnostic.plan = candidate->fingerprint;
  storage_diagnostic.valid = true;
  publish_candidate_storage_diagnostic(storage_diagnostic);
#endif
  out.release();
  out.implementation_ = candidate.release();
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kProductAllocation};
} catch (...) {
  return {StatusCode::invalid_plan, kProductInput};
}

ProductDriver::~ProductDriver() noexcept { release(); }
ProductDriver::ProductDriver(ProductDriver&& other) noexcept
    : implementation_(std::exchange(other.implementation_, nullptr)) {}
ProductDriver& ProductDriver::operator=(ProductDriver&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = std::exchange(other.implementation_, nullptr);
  }
  return *this;
}
void ProductDriver::release() noexcept {
  delete std::exchange(implementation_, nullptr);
}

Status ProductDriver::create(MPI_Comm communicator, CompiledCasePlan&& plan,
                             ProductDriver& out) noexcept try {
  if (communicator == MPI_COMM_NULL || out.implementation_ != nullptr ||
      plan.implementation_ == nullptr || !plan.summary().sealed) {
    return {StatusCode::invalid_plan, kProductInput};
  }
  std::unique_ptr<ProductDriver::Impl> candidate{
      new (std::nothrow) ProductDriver::Impl};
  if (!candidate)
    return {StatusCode::allocation_failure, kProductAllocation};
  if (MPI_Comm_dup(communicator, &candidate->communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kProductCollective};
  CompiledCasePlan::Impl& product = *plan.implementation_;
  Status status = AttemptTransaction::create(
      product.layers.field_count(), 8U, product.layers.field_count() + 8U,
      candidate->transaction);
  if (status)
    status = candidate->final_flux_authority.claim(
        product.piso.pressure_stage(), product.piso.final_flux_slot(),
        candidate->transaction, candidate->final_flux_writer);
  if (status && product.quadrature.has_value()) {
    candidate->final_force_cache.emplace();
    status = FinalForceCache::bind(
        communicator, *product.quadrature,
        product.geometry.topology_revision(), 60U, 1U,
        *candidate->final_force_cache);
  }
  if (status)
    status = TimeControllerState::start(product.time, 0.0, candidate->time);
  if (status) {
    candidate->output_fields.resize(product.io.snapshot_fields().size);
    candidate->restart_fields.resize(3U + product.fields.scalars.size());
    candidate->restart_previous_fields.resize(
        3U + product.fields.scalars.size());
    candidate->restart_expected_fields.resize(
        3U + product.fields.scalars.size());
    candidate->restart_rate_fields.resize(
        1U + product.fields.scalar_nonadvective_rates.size());
    candidate->restart_previous_rate_fields.resize(
        1U + product.fields.scalar_nonadvective_rates.size());
    candidate->restart_expected_rate_fields.resize(
        1U + product.fields.scalar_nonadvective_rates.size());
    candidate->restart_expected_fields[0U] = {
        RestartFieldRole::velocity, product.fields.velocity, 3U};
    candidate->restart_expected_fields[1U] = {
        RestartFieldRole::pressure_perturbation, product.fields.pressure, 1U};
    candidate->restart_expected_fields[2U] = {
        RestartFieldRole::enthalpy, product.fields.enthalpy, 1U};
    for (std::size_t index = 0U; index < product.fields.scalars.size(); ++index)
      candidate->restart_expected_fields[index + 3U] = {
          product.fields.scalar_roles[index] == TransportedScalarRole::species
              ? RestartFieldRole::independent_species
              : RestartFieldRole::transported_scalar,
          product.fields.scalars[index], 1U};
    candidate->restart_expected_rate_fields[0U] = {
        RestartFieldRole::enthalpy_nonadvective_rate,
        product.fields.enthalpy_nonadvective_rate, 1U};
    for (std::size_t index = 0U;
         index < product.fields.scalar_nonadvective_rates.size(); ++index)
      candidate->restart_expected_rate_fields[index + 1U] = {
          RestartFieldRole::scalar_nonadvective_rate,
          product.fields.scalar_nonadvective_rates[index], 1U};
    candidate->species_values.resize(
        product.thermodynamics.independent_species_count());
    const std::size_t species =
        product.thermodynamics.independent_species_count();
    const std::size_t passive =
        product.fields.scalars.size() - species;
    candidate->species_accepted.resize(species);
    candidate->species_previous.resize(species);
    candidate->passive_accepted.resize(passive);
    candidate->passive_previous.resize(passive);
    candidate->species_trial.resize(species);
    candidate->pressure_energy_candidate_species.resize(species);
    candidate->pressure_energy_candidate_species_const.resize(species);
    candidate->pressure_energy_candidate_species_boundary_aliases.resize(
        species);
    candidate->pressure_energy_candidate_species_history.resize(species);
    std::size_t candidate_cell_count = 0U;
    if (!detail::product_cell_count(product.patch.cells,
                                    candidate_cell_count))
      return {StatusCode::invalid_plan, kProductBinding};
    candidate->pressure_energy_baseline_energy_residual.resize(
        candidate_cell_count);
    candidate->passive_trial.resize(passive);
    candidate->species_low.resize(species);
    candidate->passive_low.resize(passive);
    candidate->species_ghosts.resize(species);
    candidate->passive_ghosts.resize(passive);
    candidate->trial_species_ghosts.resize(species);
    candidate->trial_passive_ghosts.resize(passive);
    candidate->species_history.resize(species);
    candidate->passive_history.resize(passive);
    candidate->species_rate_history.resize(species);
    candidate->passive_rate_history.resize(passive);
    candidate->species_rate_output.resize(species);
    candidate->passive_rate_output.resize(passive);
    candidate->boundary_scalars.resize(
        product.boundary.resolved_scalar_count());
    candidate->boundary_vectors.resize(
        product.boundary.resolved_vector_count());
    candidate->boundary_normal_gradients.resize(
        product.boundary.resolved_normal_gradient_count());
    candidate->boundary_thermo_authority_fields.resize(species + 2U);
    // Reserve the current terminal-rate payload plus the coupled scalar views
    // required by the same-target pressure-energy correction. The coupled
    // halo is bound separately; this vector is the allocation-free hot-path
    // staging capacity shared by both routes.
    candidate->halo_views.resize(16U + product.fields.scalars.size());
    candidate->final_dependencies.resize(6U + species);
  }
  FieldView pressure_diagonal;
  FieldView pressure_rhs;
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace,
        product.fields.pressure_diagonal, pressure_diagonal);
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace, product.fields.pressure_rhs,
        pressure_rhs);
  const std::array<FieldId, 26U> pressure_energy_workspace_fields{{
      product.fields.enthalpy_compressibility,
      product.fields.pressure_energy_c_h,
      product.fields.pressure_energy_c_h_row_scale,
      product.fields.pressure_energy_e_p,
      product.fields.pressure_energy_e_h,
      product.fields.pressure_energy_continuity_residual,
      product.fields.pressure_energy_energy_residual,
      product.fields.enthalpy_correction,
      product.fields.pressure_energy_delta_temperature,
      product.fields.pressure_energy_candidate_pressure,
      product.fields.pressure_energy_candidate_pressure_correction,
      product.fields.pressure_energy_candidate_enthalpy,
      product.fields.pressure_energy_candidate_density,
      product.fields.pressure_energy_candidate_temperature,
      product.fields.pressure_energy_candidate_velocity,
      product.fields.pressure_energy_candidate_molecular_viscosity,
      product.fields.pressure_energy_candidate_effective_viscosity,
      product.fields.pressure_energy_candidate_velocity_gradient,
      product.fields.pressure_energy_candidate_compressibility,
      product.fields.pressure_energy_candidate_enthalpy_compressibility,
      product.fields.pressure_energy_candidate_thermal_conductivity,
      product.fields.pressure_energy_candidate_heat_capacity,
      product.fields.pressure_energy_candidate_enthalpy_diffusivity,
      product.fields.schur_continuity_response,
      product.fields.schur_eliminated_enthalpy,
      product.fields.schur_energy_response,
  }};
  std::array<FieldView, pressure_energy_workspace_fields.size()>
      pressure_energy_workspace_views{};
  for (std::size_t index = 0U;
       index < pressure_energy_workspace_fields.size() && status; ++index) {
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace,
        pressure_energy_workspace_fields[index],
        pressure_energy_workspace_views[index]);
  }
  if (status)
    status = product.coupler.bind_pressure_operator(
        {candidate->communicator, &product.krylov_halo, 140U,
         product.fields.krylov_vectors, 1U},
        {pressure_diagonal, pressure_rhs}, candidate->pressure_operator);
  if (status && product.ibm_boundary.has_value() &&
      product.topology.has_value()) {
    FaceFieldView x_pressure;
    FaceFieldView y_pressure;
    FaceFieldView z_pressure;
    status = make_pressure_face_views(product.pressure_face_storage,
                                      product.patch.cells, x_pressure,
                                      y_pressure, z_pressure);
    if (status) {
      candidate->ibm_pressure_operator.emplace();
      status = IbmPressureOperator::bind_lifecycle(
          candidate->pressure_operator, product.patch.cells,
          *product.topology, *product.ibm_boundary, as_const(x_pressure),
          as_const(y_pressure), as_const(z_pressure),
          product.geometry.topology_revision(), nullptr, 0U,
          *candidate->ibm_pressure_operator);
    }
  }
  FieldView enthalpy_krylov_vectors;
  FieldView enthalpy_krylov_scalars;
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace, product.fields.krylov_vectors,
        enthalpy_krylov_vectors);
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace, product.fields.krylov_scalars,
        enthalpy_krylov_scalars);
  if (status) {
    const MgDomainActivityView activity =
        product.topology.has_value()
            ? MgDomainActivityView{
                  {product.pressure_mg_cell_activity.data(),
                   product.pressure_mg_cell_activity.size()},
                  {product.pressure_mg_x_activity.data(),
                   product.pressure_mg_x_activity.size()},
                  {product.pressure_mg_y_activity.data(),
                   product.pressure_mg_y_activity.size()},
                  {product.pressure_mg_z_activity.data(),
                   product.pressure_mg_z_activity.size()},
                  product.pressure_mg_activity_fingerprint,
                  product.pressure_mg_activity_collective}
            : MgDomainActivityView{};
    const ConservativeEnthalpyEndpointServices endpoint_services{
        candidate->communicator,
        &product.equations.kernels(),
        &product.boundary,
        product.patch,
        activity,
        &product.krylov_halo,
        &product.reductions,
        15U,
        product.auxiliary_krylov_requirements,
        enthalpy_krylov_vectors,
        enthalpy_krylov_scalars,
        product.fields.enthalpy};
    status = ConservativeEnthalpyEndpoint::bind(
        endpoint_services, candidate->enthalpy_endpoint);
  }
  if (!status) return status;
  candidate->plan = std::move(plan);
  out.implementation_ = candidate.release();
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kProductAllocation};
} catch (...) {
  return {StatusCode::invalid_plan, kProductInput};
}

Status ProductDriver::Impl::initialize_common_fields(
    const DriverInitialState& initial, double enthalpy,
    const ThermoState& thermo, const MolecularTransportState& transport,
    double heat_capacity) noexcept {
  CompiledCasePlan::Impl& product = *plan.implementation_;
  Status status;
  const std::array<StateRole, 3U> roles{{
      StateRole::accepted_n, StateRole::accepted_n_minus_one,
      StateRole::trial}};
  const std::array<double, 3U> velocity{{
      initial.velocity.x, initial.velocity.y, initial.velocity.z}};
  for (StateRole role : roles) {
    FieldView field;
    if (!(status = product.layers.view(role, product.fields.rho, field))) break;
    fill_field(field, thermo.rho);
    if (!(status = product.layers.view(role, product.fields.velocity, field)))
      break;
    fill_field(field, {velocity.data(), velocity.size()});
    if (!(status = product.layers.view(role, product.fields.pressure, field)))
      break;
    fill_field(field, 0.0);
    if (!(status = product.layers.view(role, product.fields.enthalpy, field)))
      break;
    fill_field(field, enthalpy);
    if (!(status = product.layers.view(role, product.fields.temperature,
                                       field)))
      break;
    fill_field(field, initial.temperature);
    for (std::size_t index = 0U; index < product.fields.scalars.size(); ++index) {
      if (!(status = product.layers.view(role, product.fields.scalars[index],
                                         field)))
        break;
      fill_field(field, initial.transported_scalars.data[index]);
    }
    if (status) {
      status = product.layers.view(
          role, product.fields.enthalpy_nonadvective_rate, field);
      if (status) fill_field(field, 0.0);
    }
    for (std::size_t index = 0U;
         index < product.fields.scalar_nonadvective_rates.size() && status;
         ++index) {
      status = product.layers.view(
          role, product.fields.scalar_nonadvective_rates[index], field);
      if (status) fill_field(field, 0.0);
    }
    if (!status) break;
  }
  const auto fill_runtime = [&](FieldId id, double value) {
    FieldView field;
    Status viewed = product.layers.revise_runtime(
        FieldLifetime::persistent_workspace, id);
    if (viewed)
      viewed = product.layers.runtime_view(
          FieldLifetime::persistent_workspace, id, field);
    if (viewed) fill_field(field, value);
    return viewed;
  };
  if (status)
    status = fill_runtime(product.fields.molecular_viscosity,
                          transport.viscosity);
  if (status)
    status = fill_runtime(product.fields.effective_viscosity,
                          transport.viscosity);
  if (status)
    status = fill_runtime(product.fields.compressibility,
                          thermo.drho_dp_hY);
  if (status)
    status = fill_runtime(product.fields.enthalpy_compressibility,
                          thermo.drho_dh_pY);
  if (status)
    status = fill_runtime(product.fields.thermal_conductivity,
                          transport.conductivity);
  if (status)
    status = fill_runtime(product.fields.heat_capacity, heat_capacity);
  if (status)
    status = fill_runtime(product.fields.enthalpy_diffusivity,
                          transport.conductivity / heat_capacity);
  if (status) status = fill_runtime(product.fields.predictor_density, 0.0);
  if (status)
    status = fill_runtime(product.fields.predictor_accepted_advection, 0.0);
  if (status)
    status = fill_runtime(product.fields.predictor_previous_advection, 0.0);
  if (status)
    status = fill_runtime(product.fields.predictor_low_bundle, 0.0);
  if (status) status = fill_runtime(product.fields.eos_density, thermo.rho);
  if (status) status = fill_runtime(product.fields.scalar_diffusivity, 0.0);
  if (status) status = fill_runtime(product.fields.pressure_diagonal, 0.0);
  if (status) status = fill_runtime(product.fields.pressure_rhs, 0.0);
  if (status) status = fill_runtime(product.fields.pressure_correction, 0.0);
  if (status)
    status = fill_runtime(product.fields.pressure_energy_c_h, 0.0);
  if (status)
    status = fill_runtime(product.fields.pressure_energy_c_h_row_scale, 0.0);
  if (status)
    status = fill_runtime(product.fields.pressure_energy_e_p, 0.0);
  if (status)
    status = fill_runtime(product.fields.pressure_energy_e_h, 0.0);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_continuity_residual, 0.0);
  if (status)
    status = fill_runtime(product.fields.pressure_energy_energy_residual, 0.0);
  if (status)
    status = fill_runtime(product.fields.enthalpy_correction, 0.0);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_delta_temperature, 0.0);
  if (status)
    status = fill_runtime(product.fields.pressure_energy_candidate_pressure,
                          0.0);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_candidate_pressure_correction, 0.0);
  if (status)
    status = fill_runtime(product.fields.pressure_energy_candidate_enthalpy,
                          enthalpy);
  if (status)
    status = fill_runtime(product.fields.pressure_energy_candidate_density,
                          thermo.rho);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_candidate_temperature,
        thermo.temperature);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_candidate_molecular_viscosity,
        transport.viscosity);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_candidate_effective_viscosity,
        transport.viscosity);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_candidate_compressibility,
        thermo.drho_dp_hY);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_candidate_enthalpy_compressibility,
        thermo.drho_dh_pY);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_candidate_thermal_conductivity,
        transport.conductivity);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_candidate_heat_capacity,
        heat_capacity);
  if (status)
    status = fill_runtime(
        product.fields.pressure_energy_candidate_enthalpy_diffusivity,
        transport.conductivity / heat_capacity);
  std::size_t candidate_species_index = 0U;
  for (std::size_t scalar = 0U;
       scalar < product.fields.scalars.size() && status; ++scalar) {
    if (product.fields.scalar_roles[scalar] !=
        TransportedScalarRole::species)
      continue;
    if (candidate_species_index >=
        product.fields.pressure_energy_candidate_species.size()) {
      status = {StatusCode::invalid_plan, kProductBinding};
      break;
    }
    status = fill_runtime(
        product.fields.pressure_energy_candidate_species
            [candidate_species_index++],
        initial.transported_scalars.data[scalar]);
  }
  if (status && candidate_species_index !=
                    product.fields.pressure_energy_candidate_species.size())
    status = {StatusCode::invalid_plan, kProductBinding};
  if (status)
    status = fill_runtime(product.fields.schur_continuity_response, 0.0);
  if (status)
    status = fill_runtime(product.fields.schur_eliminated_enthalpy, 0.0);
  if (status)
    status = fill_runtime(product.fields.schur_energy_response, 0.0);
  FieldView vector_field;
  const std::array<double, 3U> zero_vector{{0.0, 0.0, 0.0}};
  const std::array<double, 3U> unit_vector{{1.0, 1.0, 1.0}};
  const auto vector_runtime_view_for_write = [&](FieldId id, FieldView& out) {
    Status revised = product.layers.revise_runtime(
        FieldLifetime::persistent_workspace, id);
    return revised ? product.layers.runtime_view(
                         FieldLifetime::persistent_workspace, id, out)
                   : revised;
  };
  if (status) {
    status = vector_runtime_view_for_write(product.fields.momentum_diagonal,
                                           vector_field);
    if (status)
      fill_field(vector_field, {unit_vector.data(), unit_vector.size()});
  }
  if (status) {
    status = vector_runtime_view_for_write(product.fields.momentum_rhs,
                                           vector_field);
    if (status)
      fill_field(vector_field, {zero_vector.data(), zero_vector.size()});
  }
  if (status) {
    status = vector_runtime_view_for_write(product.fields.momentum_residual,
                                           vector_field);
    if (status)
      fill_field(vector_field, {zero_vector.data(), zero_vector.size()});
  }
  if (status) {
    status = vector_runtime_view_for_write(product.fields.r_au, vector_field);
    if (status)
      fill_field(vector_field, {unit_vector.data(), unit_vector.size()});
  }
  if (status) {
    status = vector_runtime_view_for_write(product.fields.h_by_a,
                                           vector_field);
    if (status) fill_field(vector_field, {velocity.data(), velocity.size()});
  }
  if (status) {
    status = vector_runtime_view_for_write(product.fields.pressure_gradient,
                                           vector_field);
    if (status)
      fill_field(vector_field, {zero_vector.data(), zero_vector.size()});
  }
  if (status) {
    status = vector_runtime_view_for_write(
        product.fields.pressure_energy_candidate_velocity, vector_field);
    if (status) fill_field(vector_field, {velocity.data(), velocity.size()});
  }
  if (status) {
    status = vector_runtime_view_for_write(
        product.fields.pressure_energy_candidate_velocity_gradient,
        vector_field);
    if (status) {
      const std::array<double, 9U> zero_gradient{};
      fill_field(vector_field,
                 {zero_gradient.data(), zero_gradient.size()});
    }
  }
  if (status) {
    status = vector_runtime_view_for_write(product.fields.velocity_gradient,
                                           vector_field);
    if (status) {
      const std::array<double, 9U> zero_gradient{};
      fill_field(vector_field, {zero_gradient.data(), zero_gradient.size()});
    }
  }
  return status;
}

Status ProductDriver::Impl::rebuild_cold_velocity_dependents(
    double cold_pressure_reference, bool restart) noexcept {
  CompiledCasePlan::Impl& product = *plan.implementation_;
  const Int3 cells = product.patch.cells;
  const KernelBox full_box{{0, 0, 0}, cells};
  const std::array<StateRole, 3U> roles{{
      StateRole::accepted_n, StateRole::accepted_n_minus_one,
      StateRole::trial}};
  const std::array<StageId, 3U> velocity_halo_stages{{
      kFreshVelocityAcceptedHaloStage, kFreshVelocityPreviousHaloStage,
      kFreshVelocityTrialHaloStage}};
  Status status;
  std::array<FieldView, 3U> velocity_layers{};
  for (std::size_t layer = 0U; layer < roles.size() && status; ++layer)
    status = product.layers.view(roles[layer], product.fields.velocity,
                                 velocity_layers[layer]);

  const StateLayers& const_layers = product.layers;
  ConstFieldView accepted_density;
  ConstFieldView accepted_enthalpy;
  if (status)
    status = const_layers.view(StateRole::accepted_n, product.fields.rho,
                               accepted_density);
  if (status)
    status = const_layers.view(StateRole::accepted_n,
                               product.fields.enthalpy, accepted_enthalpy);
  std::size_t species_index = 0U;
  for (std::size_t index = 0U;
       index < product.fields.scalars.size() && status; ++index) {
    if (product.fields.scalar_roles[index] !=
        TransportedScalarRole::species)
      continue;
    status = const_layers.view(StateRole::accepted_n,
                               product.fields.scalars[index],
                               species_accepted[species_index]);
    if (status && species_accepted[species_index].base != nullptr)
      species_values[species_index] =
          species_accepted[species_index].unchecked({0, 0, 0}, 0U);
    ++species_index;
  }
  if (status && species_index != species_accepted.size())
    status = {StatusCode::invalid_plan, kProductBinding};
  // The boundary resolver can enter mass-flow collectives.  A rank-local
  // state-layer/view failure must be made uniform before any rank calls it.
  status = product.reductions.consensus(status);
  if (!status) return status;
  if (status)
    status = resolve_static_boundary_values(
        communicator, product.boundary, product.boundary_specs,
        product.geometry, product.patch, product.thermodynamics,
        cold_pressure_reference, accepted_density,
        as_const(velocity_layers[0U]), accepted_enthalpy,
        {species_accepted.data(), species_accepted.size()},
        {passive_accepted.data(), passive_accepted.size()},
        {species_values.data(), species_values.size()}, boundary_scalars,
        boundary_vectors, boundary_normal_gradients, true);
  status = product.reductions.consensus(status);
  if (!status) return status;

  // Synchronize all three identical cold velocity histories.  Remote IBM
  // donor data is refreshed independently for each layer; only ghosts are
  // touched, so Restart checkpoint interiors remain byte-for-byte intact.
  for (std::size_t layer = 0U; layer < velocity_layers.size(); ++layer) {
    std::array<FieldView, 1U> fields{velocity_layers[layer]};
    HaloTicket ticket;
    status = product.final_velocity_halo.begin(
        velocity_halo_stages[layer], {fields.data(), fields.size()}, ticket);
    if (status)
      status = product.final_velocity_halo.finish(
          ticket, {fields.data(), fields.size()});
    if (status) {
      velocity_layers[layer] = fields[0U];
      status = apply_boundary_ghosts(
          BoundaryStage::momentum, product.boundary,
          {&velocity_layers[layer], 1U}, resolved_boundary_values());
    }
    status = product.reductions.consensus(status);
    if (!status) return status;
    if (product.ibm_gradient_donors.has_value()) {
      fields[0U] = velocity_layers[layer];
      status = product.ibm_gradient_donors->exchange(
          kIbmGradientDonorStage, {fields.data(), fields.size()});
      if (status) velocity_layers[layer] = fields[0U];
      status = product.reductions.consensus(status);
      if (!status) return status;
    }
  }

  const auto runtime_view_for_write = [&](FieldId field, FieldView& view) {
    Status revised = product.layers.revise_runtime(
        FieldLifetime::persistent_workspace, field);
    return revised ? product.layers.runtime_view(
                         FieldLifetime::persistent_workspace, field, view)
                   : revised;
  };
  FieldView velocity_gradient;
  FieldView molecular_viscosity;
  FieldView effective_viscosity;
  FieldView conductivity;
  FieldView heat_capacity;
  FieldView enthalpy_diffusivity;
  status = runtime_view_for_write(product.fields.velocity_gradient,
                                  velocity_gradient);
  if (status) {
    const std::array<ConstFieldView, 1U> reads{
        as_const(velocity_layers[0U])};
    const std::array<FieldView, 1U> writes{velocity_gradient};
    status = cartesian_gradient(
        product.equations.kernels(),
        {{reads.data(), reads.size()}, {writes.data(), writes.size()},
         full_box, 0U, 0U, 3U, 0U, nullptr});
  }
  if (status && product.ibm_equations.has_value())
    status = product.ibm_equations->correct_velocity_gradient(
        as_const(velocity_layers[0U]), velocity_gradient);
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace,
        product.fields.molecular_viscosity, molecular_viscosity);
  if (status)
    status = runtime_view_for_write(product.fields.effective_viscosity,
                                    effective_viscosity);
  TurbulenceCertificate turbulence_certificate;
  if (status) {
    const TurbulenceUpdateInput turbulence_input{
        accepted_density, as_const(molecular_viscosity), {},
        velocity_gradient.revision, as_const(velocity_gradient)};
    status = product.turbulence.update(turbulence_input, effective_viscosity,
                                       turbulence_certificate);
  }
  if (status && product.transport.kernel() ==
                    TransportKernel::coast_native_air)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace, product.fields.heat_capacity,
        heat_capacity);
  if (status && product.transport.kernel() ==
                    TransportKernel::coast_native_air)
    status = runtime_view_for_write(product.fields.thermal_conductivity,
                                    conductivity);
  if (status && product.transport.kernel() ==
                    TransportKernel::coast_native_air)
    status = runtime_view_for_write(product.fields.enthalpy_diffusivity,
                                    enthalpy_diffusivity);
  if (status && product.transport.kernel() ==
                    TransportKernel::coast_native_air)
    status = refresh_coast_native_air_effective_thermal_transport(
        product.transport, as_const(molecular_viscosity),
        as_const(effective_viscosity), as_const(heat_capacity), conductivity,
        enthalpy_diffusivity);
  if (status && product.transport.kernel() !=
                    TransportKernel::coast_native_air)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace,
        product.fields.thermal_conductivity, conductivity);
  if (product.transport.kernel() == TransportKernel::coast_native_air)
    status = exchange_effective_thermal_ghosts(
        product.coupled_thermal_halo, 60U, product.boundary, conductivity,
        enthalpy_diffusivity, status);
  status = product.reductions.consensus(status);
  if (!status) return status;

  FieldView trial_density;
  FieldView trial_pressure;
  if (status)
    status = product.layers.view(StateRole::trial, product.fields.rho,
                                 trial_density);
  if (status)
    status = product.layers.view(StateRole::trial, product.fields.pressure,
                                 trial_pressure);
  // stage_halos[3] is a communication stage; do not let a local view failure
  // decide whether a rank enters it.
  status = product.reductions.consensus(status);
  if (!status) return status;
  if (status) {
    std::array<FieldView, 5U> fields{
        trial_density, velocity_layers[2U], trial_pressure,
        effective_viscosity, velocity_gradient};
    HaloTicket ticket;
    status = product.stage_halos[3U].begin(
        30U, {fields.data(), fields.size()}, ticket);
    if (status)
      status = product.stage_halos[3U].finish(
          ticket, {fields.data(), fields.size()});
    if (status) {
      trial_density = fields[0U];
      velocity_layers[2U] = fields[1U];
      trial_pressure = fields[2U];
      effective_viscosity = fields[3U];
      velocity_gradient = fields[4U];
      std::array<FieldView, 3U> zero_gradient{
          trial_density, effective_viscosity, velocity_gradient};
      status = apply_physical_zero_gradient(
          product.boundary,
          {zero_gradient.data(), zero_gradient.size()});
      if (status)
        status = apply_boundary_ghosts(
            BoundaryStage::pressure, product.boundary, {&trial_pressure, 1U},
            resolved_boundary_values());
    }
  }
  status = product.reductions.consensus(status);
  if (!status) return status;

  FieldView trial_enthalpy;
  FieldView trial_temperature;
  if (status)
    status = product.layers.view(StateRole::trial, product.fields.enthalpy,
                                 trial_enthalpy);
  if (status)
    status = product.layers.view(StateRole::trial, product.fields.temperature,
                                 trial_temperature);
  species_index = 0U;
  std::size_t passive_index = 0U;
  for (std::size_t index = 0U;
       index < product.fields.scalars.size() && status; ++index) {
    FieldView scalar;
    status = product.layers.view(StateRole::trial,
                                 product.fields.scalars[index], scalar);
    if (!status) break;
    if (product.fields.scalar_roles[index] ==
        TransportedScalarRole::species)
      species_trial[species_index++] = scalar;
    else
      passive_trial[passive_index++] = scalar;
  }
  const std::size_t rate_halo_count = 7U + product.fields.scalars.size();
  if (status && rate_halo_count > halo_views.size())
    status = {StatusCode::invalid_plan, kProductCommunication};
  // Likewise, all local scalar/view checks precede one common stage-60 halo
  // decision.
  status = product.reductions.consensus(status);
  if (!status) return status;
  if (status) {
    std::size_t halo_count = 0U;
    halo_views[halo_count++] = trial_pressure;
    halo_views[halo_count++] = trial_enthalpy;
    halo_views[halo_count++] = trial_temperature;
    halo_views[halo_count++] = velocity_gradient;
    halo_views[halo_count++] = molecular_viscosity;
    halo_views[halo_count++] = effective_viscosity;
    halo_views[halo_count++] = conductivity;
    for (FieldView value : species_trial) halo_views[halo_count++] = value;
    for (FieldView value : passive_trial) halo_views[halo_count++] = value;
    HaloTicket ticket;
    status = product.stage_halos[5U].begin(
        60U, {halo_views.data(), halo_count}, ticket);
    if (status)
      status = product.stage_halos[5U].finish(
          ticket, {halo_views.data(), halo_count});
    if (status) {
      trial_pressure = halo_views[0U];
      trial_enthalpy = halo_views[1U];
      trial_temperature = halo_views[2U];
      velocity_gradient = halo_views[3U];
      molecular_viscosity = halo_views[4U];
      effective_viscosity = halo_views[5U];
      conductivity = halo_views[6U];
      species_index = 0U;
      passive_index = 0U;
      for (std::size_t index = 0U; index < product.fields.scalars.size();
           ++index) {
        if (product.fields.scalar_roles[index] ==
            TransportedScalarRole::species)
          species_trial[species_index++] = halo_views[index + 7U];
        else
          passive_trial[passive_index++] = halo_views[index + 7U];
      }
      std::array<FieldView, 5U> derived{
          trial_temperature, velocity_gradient, molecular_viscosity,
          effective_viscosity, conductivity};
      status = apply_physical_zero_gradient(
          product.boundary, {derived.data(), derived.size()});
      if (status)
        status = apply_boundary_ghosts(
            BoundaryStage::pressure, product.boundary, {&trial_pressure, 1U},
            resolved_boundary_values());
      if (status)
        status = apply_boundary_ghosts(
            BoundaryStage::enthalpy, product.boundary,
            {&trial_enthalpy, 1U}, resolved_boundary_values());
      if (status && !product.fields.scalars.empty()) {
        // Rebuild the scalar list in schema order; the compact species/passive
        // arrays above are retained for the equation-state binding below.
        species_index = 0U;
        passive_index = 0U;
        for (std::size_t index = 0U; index < product.fields.scalars.size();
             ++index)
          halo_views[index] =
              product.fields.scalar_roles[index] ==
                      TransportedScalarRole::species
                  ? species_trial[species_index++]
                  : passive_trial[passive_index++];
        if (status)
          status = apply_boundary_ghosts(
              BoundaryStage::scalar, product.boundary,
              {halo_views.data(), product.fields.scalars.size()},
              resolved_boundary_values());
      }
    }
  }
  status = product.reductions.consensus(status);
  if (!status) return status;

  if (product.ibm_rate_donors.has_value()) {
    std::size_t halo_count = 0U;
    halo_views[halo_count++] = trial_pressure;
    halo_views[halo_count++] = trial_temperature;
    for (FieldView value : species_trial) halo_views[halo_count++] = value;
    for (FieldView value : passive_trial) halo_views[halo_count++] = value;
    status = product.ibm_rate_donors->exchange(
        161U, {halo_views.data(), halo_count});
    if (status) {
      trial_pressure = halo_views[0U];
      trial_temperature = halo_views[1U];
      species_index = 0U;
      passive_index = 0U;
      for (std::size_t index = 0U; index < product.fields.scalars.size();
           ++index) {
        if (product.fields.scalar_roles[index] ==
            TransportedScalarRole::species)
          species_trial[species_index++] = halo_views[index + 2U];
        else
          passive_trial[passive_index++] = halo_views[index + 2U];
      }
    }
    status = product.reductions.consensus(status);
    if (!status) return status;
  }

  const auto history = [&](FieldId field, PrimitiveHistory& out) {
    return const_layers.view(StateRole::trial, field, out.trial) &&
                   const_layers.view(StateRole::accepted_n, field,
                                     out.accepted) &&
                   const_layers.view(StateRole::accepted_n_minus_one, field,
                                     out.previous)
               ? Status{}
               : Status{StatusCode::invalid_plan, kProductBinding};
  };
  EquationStateView equation_state;
  if (status) status = history(product.fields.rho, equation_state.density);
  if (status)
    status = history(product.fields.velocity, equation_state.velocity);
  if (status)
    status = history(product.fields.pressure,
                     equation_state.pressure_perturbation);
  if (status)
    status = history(product.fields.enthalpy, equation_state.enthalpy);
  if (status)
    status = history(product.fields.temperature, equation_state.temperature);
  species_index = 0U;
  passive_index = 0U;
  for (std::size_t index = 0U;
       index < product.fields.scalars.size() && status; ++index) {
    PrimitiveHistory scalar_history;
    status = history(product.fields.scalars[index], scalar_history);
    if (!status) break;
    if (product.fields.scalar_roles[index] ==
        TransportedScalarRole::species)
      species_history[species_index++] = scalar_history;
    else
      passive_history[passive_index++] = scalar_history;
  }
  equation_state.independent_species = {species_history.data(),
                                        species_history.size()};
  equation_state.passive_scalars = {passive_history.data(),
                                    passive_history.size()};
  equation_state.pressure_reference = cold_pressure_reference;
  equation_state.accepted_pressure_reference = cold_pressure_reference;
  equation_state.previous_pressure_reference = cold_pressure_reference;

  FieldView enthalpy_rate_trial;
  if (status)
    status = product.layers.view(StateRole::trial,
                                 product.fields.enthalpy_nonadvective_rate,
                                 enthalpy_rate_trial);
  species_index = 0U;
  passive_index = 0U;
  for (std::size_t index = 0U;
       index < product.fields.scalar_nonadvective_rates.size() && status;
       ++index) {
    FieldView rate;
    status = product.layers.view(
        StateRole::trial, product.fields.scalar_nonadvective_rates[index],
        rate);
    if (!status) break;
    if (product.fields.scalar_roles[index] ==
        TransportedScalarRole::species)
      species_rate_output[species_index++] = rate;
    else
      passive_rate_output[passive_index++] = rate;
  }
  FieldView diffusion_scratch;
  FieldView scalar_diffusivity;
  if (status)
    status = runtime_view_for_write(
        product.fields.predictor_accepted_advection, diffusion_scratch);
  if (status)
    status = runtime_view_for_write(product.fields.scalar_diffusivity,
                                    scalar_diffusivity);
  EquationMaterialView material;
  material.molecular_viscosity = as_const(molecular_viscosity);
  material.effective_viscosity = as_const(effective_viscosity);
  material.thermal_conductivity = as_const(conductivity);
  ThermophysicalRateCertificate rate_certificate;
  if (status)
    status = evaluate_thermophysical_rates(
        product.equations,
        {equation_state, material, as_const(velocity_gradient),
         {1.0, -1.0, 0.0, 1U}, 1U, {},
         product.ibm_equations.has_value() ? &*product.ibm_equations
                                           : nullptr},
        {enthalpy_rate_trial,
         {species_rate_output.data(), species_rate_output.size()},
         {passive_rate_output.data(), passive_rate_output.size()},
         diffusion_scratch, scalar_diffusivity},
        rate_certificate);
  status = product.reductions.consensus(status);
  if (!status) return status;

  const auto copy_rate_to_history = [&](FieldId field,
                                        ConstFieldView source) noexcept {
    for (StateRole role : {StateRole::accepted_n,
                           StateRole::accepted_n_minus_one}) {
      FieldView destination;
      Status copied = product.layers.view(role, field, destination);
      if (!copied) return copied;
      for (std::int32_t z = 0; z < cells.z; ++z)
        for (std::int32_t y = 0; y < cells.y; ++y)
          for (std::int32_t x = 0; x < cells.x; ++x)
            destination.unchecked({x, y, z}, 0U) =
                source.unchecked({x, y, z}, 0U);
    }
    return Status{};
  };
  status = copy_rate_to_history(product.fields.enthalpy_nonadvective_rate,
                                as_const(enthalpy_rate_trial));
  for (std::size_t index = 0U;
       index < product.fields.scalar_nonadvective_rates.size() && status;
       ++index) {
    ConstFieldView rate;
    status = const_layers.view(
        StateRole::trial, product.fields.scalar_nonadvective_rates[index],
        rate);
    if (status)
      status = copy_rate_to_history(
          product.fields.scalar_nonadvective_rates[index], rate);
  }
  FieldView h_by_a;
  if (status)
    status = runtime_view_for_write(product.fields.h_by_a, h_by_a);
  if (status) {
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          for (std::uint8_t component = 0U; component < 3U; ++component)
            h_by_a.unchecked({x, y, z}, component) =
                velocity_layers[0U].unchecked({x, y, z}, component);
  }
  status = product.reductions.consensus(status);
  if (!status) return status;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  detail::ColdVelocityDependentsDiagnostic cold_diagnostic;
  cold_diagnostic.valid = true;
  cold_diagnostic.restart = restart;
  cold_diagnostic.rebuilt = true;
  cold_diagnostic.driver_plan = product.fingerprint;
  double local_h_by_a_difference = 0.0;
  std::uint64_t local_gradient_nonfinite = 0U;
  std::uint64_t local_effective_nonpositive = 0U;
  std::uint64_t local_solid_gradient_nonfinite = 0U;
  bool local_rate_layers_equal = true;
  std::uint64_t local_lineage = kFnvOffset;
  const Span<const std::uint8_t> cold_activity =
      product.topology.has_value() ? product.topology->region()
                                   : Span<const std::uint8_t>{};
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const std::size_t flat =
            static_cast<std::size_t>(x) +
            static_cast<std::size_t>(cells.x) *
                (static_cast<std::size_t>(y) +
                 static_cast<std::size_t>(cells.y) *
                     static_cast<std::size_t>(z));
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          local_h_by_a_difference = std::max(
              local_h_by_a_difference,
              std::abs(h_by_a.unchecked(cell, component) -
                       velocity_layers[0U].unchecked(cell, component)));
          local_lineage = detail::product_mix(
              local_lineage,
              product_double_bits(h_by_a.unchecked(cell, component)));
        }
        for (std::uint8_t component = 0U; component < 9U; ++component) {
          const double value = velocity_gradient.unchecked(cell, component);
          if (!std::isfinite(value)) {
            ++local_gradient_nonfinite;
            if (cold_activity.size != 0U &&
                cold_activity.data[flat] == 0U)
              ++local_solid_gradient_nonfinite;
          }
          local_lineage =
              detail::product_mix(local_lineage, product_double_bits(value));
        }
        const double effective =
            effective_viscosity.unchecked(cell, 0U);
        if (!std::isfinite(effective) || !(effective > 0.0))
          ++local_effective_nonpositive;
        local_lineage = detail::product_mix(
            local_lineage, product_double_bits(effective));
      }
  const auto inspect_rate_layers = [&](FieldId field) noexcept {
    std::array<ConstFieldView, 3U> layers{};
    for (std::size_t layer = 0U; layer < roles.size(); ++layer)
      if (!const_layers.view(roles[layer], field, layers[layer]))
        return false;
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const std::uint64_t bits =
              product_double_bits(layers[0U].unchecked(cell, 0U));
          local_lineage = detail::product_mix(local_lineage, bits);
          if (bits != product_double_bits(layers[1U].unchecked(cell, 0U)) ||
              bits != product_double_bits(layers[2U].unchecked(cell, 0U)))
            return false;
        }
    return true;
  };
  local_rate_layers_equal =
      inspect_rate_layers(product.fields.enthalpy_nonadvective_rate);
  for (FieldId field : product.fields.scalar_nonadvective_rates)
    local_rate_layers_equal = inspect_rate_layers(field) &&
                              local_rate_layers_equal;
  std::array<std::uint64_t, 3U> local_counts{
      local_gradient_nonfinite, local_effective_nonpositive,
      local_solid_gradient_nonfinite};
  std::array<std::uint64_t, 3U> global_counts{};
  double global_h_by_a_difference = 0.0;
  int local_rates_equal = local_rate_layers_equal ? 1 : 0;
  int global_rates_equal = 0;
  std::uint64_t lineage_xor = 0U;
  std::uint64_t lineage_sum = 0U;
  const int envelope_rc = MPI_Allreduce(
      &local_h_by_a_difference, &global_h_by_a_difference, 1, MPI_DOUBLE,
      MPI_MAX, communicator);
  const int counts_rc = MPI_Allreduce(
      local_counts.data(), global_counts.data(),
      static_cast<int>(local_counts.size()), MPI_UINT64_T, MPI_SUM,
      communicator);
  const int rates_rc = MPI_Allreduce(&local_rates_equal, &global_rates_equal,
                                     1, MPI_INT, MPI_MIN, communicator);
  const int xor_rc = MPI_Allreduce(&local_lineage, &lineage_xor, 1,
                                   MPI_UINT64_T, MPI_BXOR, communicator);
  const int sum_rc = MPI_Allreduce(&local_lineage, &lineage_sum, 1,
                                   MPI_UINT64_T, MPI_SUM, communicator);
  const bool diagnostic_collective =
      envelope_rc == MPI_SUCCESS && counts_rc == MPI_SUCCESS &&
      rates_rc == MPI_SUCCESS && xor_rc == MPI_SUCCESS &&
      sum_rc == MPI_SUCCESS;
  cold_diagnostic.valid = diagnostic_collective;
  cold_diagnostic.maximum_h_by_a_velocity_difference =
      global_h_by_a_difference;
  cold_diagnostic.velocity_gradient_nonfinite_count = global_counts[0U];
  cold_diagnostic.effective_viscosity_nonpositive_count = global_counts[1U];
  cold_diagnostic.solid_velocity_gradient_nonfinite_count = global_counts[2U];
  cold_diagnostic.rate_layers_bitwise_equal = global_rates_equal != 0;
  cold_diagnostic.lineage = detail::product_mix(
      detail::product_mix(product.fingerprint, lineage_xor), lineage_sum);
  if (cold_diagnostic.lineage == 0U) cold_diagnostic.lineage = 1U;
  g_cold_velocity_dependents_diagnostic = cold_diagnostic;
  g_cold_velocity_dependents_published.store(
      cold_diagnostic.valid, std::memory_order_release);
#else
  static_cast<void>(restart);
#endif
  return {};
}

Status ProductDriver::restart_expected(RestartExpected& out) noexcept {
  if (implementation_ == nullptr || implementation_->initialized)
    return {StatusCode::invalid_plan, kProductInput};
  const ProductDriver::Impl& runtime = *implementation_;
  const CompiledCasePlan::Impl& product = *runtime.plan.implementation_;
  out = {product.geometry.global_cells(),
         product.patch,
         runtime.plan.fingerprint(),
         product.schema_fingerprint,
         product.geometry.fingerprint(),
         {runtime.restart_expected_fields.data(),
          runtime.restart_expected_fields.size()},
         {runtime.restart_expected_rate_fields.data(),
          runtime.restart_expected_rate_fields.size()}};
  return {};
}

Status ProductDriver::initialize(const DriverInitialState& initial) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  // Each initialization attempt owns a fresh diagnostic epoch.  In
  // particular, an empty-activity Fresh bypass must never inherit an IBM or
  // Restart derived-state certificate from an earlier driver with the same
  // compiled-plan fingerprint.
  g_cold_velocity_dependents_published.store(false,
                                             std::memory_order_release);
  g_cold_velocity_dependents_diagnostic = {};
  g_fresh_initialization_published.store(false, std::memory_order_release);
  g_fresh_initialization_diagnostic = {};
#endif
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kProductInput};
  }
  ProductDriver::Impl& runtime = *implementation_;
  CompiledCasePlan::Impl& product = *runtime.plan.implementation_;
  Status status =
      !runtime.initialized && std::isfinite(initial.pressure_reference) &&
              initial.pressure_reference > 0.0 &&
              std::isfinite(initial.temperature) &&
              std::isfinite(initial.velocity.x) &&
              std::isfinite(initial.velocity.y) &&
              std::isfinite(initial.velocity.z) &&
              std::isfinite(initial.start_time) &&
              initial.transported_scalars.size ==
                  product.fields.scalars.size() &&
              (initial.transported_scalars.size == 0U ||
               initial.transported_scalars.data != nullptr)
          ? Status{}
          : Status{StatusCode::invalid_plan, kProductInput};
  std::size_t species_index = 0U;
  for (std::size_t index = 0U;
       index < product.fields.scalars.size() && status; ++index) {
    const double value = initial.transported_scalars.data[index];
    if (!std::isfinite(value) || value < 0.0 || value > 1.0)
      status = {StatusCode::numerical_failure, kProductInput};
    if (status && product.fields.scalar_roles[index] ==
                      TransportedScalarRole::species)
      ++species_index;
  }
  if (status && species_index != runtime.species_values.size())
    status = {StatusCode::invalid_plan, kProductInput};
  status = product.reductions.consensus(status);
  if (!status) return status;

  std::uint64_t initial_contract =
      detail::product_mix(kFnvOffset, UINT64_C(0x6672657368696e69));
  initial_contract = detail::product_mix(
      initial_contract, product_double_bits(initial.pressure_reference));
  initial_contract = detail::product_mix(
      initial_contract, product_double_bits(initial.temperature));
  initial_contract = detail::product_mix(
      initial_contract, product_double_bits(initial.velocity.x));
  initial_contract = detail::product_mix(
      initial_contract, product_double_bits(initial.velocity.y));
  initial_contract = detail::product_mix(
      initial_contract, product_double_bits(initial.velocity.z));
  initial_contract = detail::product_mix(
      initial_contract, product_double_bits(initial.start_time));
  initial_contract = detail::product_mix(
      initial_contract, initial.transported_scalars.size);
  for (std::size_t index = 0U; index < initial.transported_scalars.size;
       ++index)
    initial_contract = detail::product_mix(
        initial_contract,
        product_double_bits(initial.transported_scalars.data[index]));
  if (initial_contract == 0U) initial_contract = 1U;
  status = product.reductions.consensus_contract(initial_contract);
  if (!status) return status;

  species_index = 0U;
  for (std::size_t index = 0U; index < product.fields.scalars.size(); ++index)
    if (product.fields.scalar_roles[index] ==
        TransportedScalarRole::species)
      runtime.species_values[species_index++] =
          initial.transported_scalars.data[index];
  double enthalpy = 0.0;
  double heat_capacity = 0.0;
  double gas_constant = 0.0;
  status = product.thermodynamics.mixture_enthalpy(
      initial.temperature,
      {runtime.species_values.data(), runtime.species_values.size()}, enthalpy,
      heat_capacity, gas_constant);
  ThermoState thermo;
  if (status)
    status = product.thermodynamics.evaluate(
        initial.pressure_reference, enthalpy,
        {runtime.species_values.data(), runtime.species_values.size()},
        initial.velocity, thermo, initial.temperature);
  MolecularTransportState transport;
  if (status)
    status = product.transport.evaluate(
        initial.temperature,
        {runtime.species_values.data(), runtime.species_values.size()},
        transport);
  if (!status || !std::isfinite(thermo.rho) || thermo.rho <= 0.0 ||
      !std::isfinite(thermo.drho_dp_hY) || thermo.drho_dp_hY <= 0.0 ||
      !std::isfinite(thermo.drho_dh_pY) || thermo.drho_dh_pY >= 0.0 ||
      !std::isfinite(transport.viscosity) || transport.viscosity <= 0.0 ||
      !std::isfinite(transport.conductivity) ||
      transport.conductivity <= 0.0 || heat_capacity <= 0.0) {
    if (status)
      status = {StatusCode::numerical_failure, kProductInput};
  }
  status = product.reductions.consensus(status);
  if (!status) return status;
  double local_volume = 0.0;
  const Span<const double> dx = product.geometry.x().widths();
  const Span<const double> dy = product.geometry.y().widths();
  const Span<const double> dz = product.geometry.z().widths();
  const Span<const std::uint8_t> active =
      product.topology.has_value() ? product.topology->region()
                                   : Span<const std::uint8_t>{};
  for (std::int32_t z = 0; z < product.patch.cells.z; ++z)
    for (std::int32_t y = 0; y < product.patch.cells.y; ++y)
      for (std::int32_t x = 0; x < product.patch.cells.x; ++x)
        if (active.size == 0U ||
            active.data[static_cast<std::size_t>(x) +
                        static_cast<std::size_t>(product.patch.cells.x) *
                            (static_cast<std::size_t>(y) +
                             static_cast<std::size_t>(product.patch.cells.y) *
                                 z)] != 0U)
          local_volume +=
              dx.data[static_cast<std::size_t>(product.patch.begin.x + x)] *
              dy.data[static_cast<std::size_t>(product.patch.begin.y + y)] *
              dz.data[static_cast<std::size_t>(product.patch.begin.z + z)];
  const double local_mass = thermo.rho * local_volume;
  double global_mass = 0.0;
  status = product.reductions.checked_sum({&local_mass, 1U},
                                          {&global_mass, 1U});
  if (status && (!std::isfinite(global_mass) || global_mass <= 0.0))
    status = {StatusCode::numerical_failure, kProductInput};
  if (!status) return status;
  TimeControllerState controller;
  status = TimeControllerState::start(product.time, initial.start_time,
                                      controller);
  status = product.reductions.consensus(status);
  if (!status) return status;

  status = runtime.initialize_common_fields(initial, enthalpy, thermo,
                                            transport, heat_capacity);
  // Every rank must reach the same boundary-resolution sequence.  Field
  // binding is rank-local, while mass-flow boundary resolution may execute
  // collectives, so a local initialization failure cannot be allowed to
  // decide whether the next collective is entered.
  status = product.reductions.consensus(status);
  if (!status) return status;

  // The predictor treats the committed face flux as the fixed-flux boundary
  // authority.  Seed that history from the resolved physical velocity state,
  // rather than from the unconstrained interior initial value.  This is a
  // cold-start service; mass-flow normalization adds no hot-step collective.
  const StateLayers& const_layers = product.layers;
  ConstFieldView initial_density;
  ConstFieldView initial_velocity;
  ConstFieldView initial_enthalpy;
  if (status)
    status = const_layers.view(StateRole::accepted_n, product.fields.rho,
                               initial_density);
  if (status)
    status = const_layers.view(StateRole::accepted_n, product.fields.velocity,
                               initial_velocity);
  if (status)
    status = const_layers.view(StateRole::accepted_n, product.fields.enthalpy,
                               initial_enthalpy);
  species_index = 0U;
  std::size_t passive_index = 0U;
  for (std::size_t index = 0U;
       index < product.fields.scalars.size() && status; ++index) {
    if (product.fields.scalar_roles[index] ==
        TransportedScalarRole::species) {
      status = const_layers.view(StateRole::accepted_n,
                                 product.fields.scalars[index],
                                 runtime.species_accepted[species_index++]);
    } else {
      status = const_layers.view(StateRole::accepted_n,
                                 product.fields.scalars[index],
                                 runtime.passive_accepted[passive_index++]);
    }
  }
  if (status &&
      (species_index != runtime.species_accepted.size() ||
       passive_index != runtime.passive_accepted.size()))
    status = {StatusCode::invalid_plan, kProductBinding};
  status = product.reductions.consensus(status);
  if (status)
    status = resolve_static_boundary_values(
        runtime.communicator, product.boundary, product.boundary_specs,
        product.geometry, product.patch, product.thermodynamics,
        initial.pressure_reference, initial_density, initial_velocity,
        initial_enthalpy,
        {runtime.species_accepted.data(), runtime.species_accepted.size()},
        {runtime.passive_accepted.data(), runtime.passive_accepted.size()},
        {runtime.species_values.data(), runtime.species_values.size()},
        runtime.boundary_scalars, runtime.boundary_vectors,
        runtime.boundary_normal_gradients, true);
  FieldView initial_density_with_ghosts;
  FieldView initial_pressure_with_ghosts;
  FieldView initial_enthalpy_with_ghosts;
  FieldView initial_temperature_with_ghosts;
  FieldView initial_velocity_with_ghosts;
  if (status)
    status = product.layers.view(StateRole::accepted_n, product.fields.rho,
                                 initial_density_with_ghosts);
  if (status)
    status = product.layers.view(StateRole::accepted_n,
                                 product.fields.pressure,
                                 initial_pressure_with_ghosts);
  if (status)
    status = product.layers.view(StateRole::accepted_n,
                                 product.fields.enthalpy,
                                 initial_enthalpy_with_ghosts);
  if (status)
    status = product.layers.view(StateRole::accepted_n,
                                 product.fields.temperature,
                                 initial_temperature_with_ghosts);
  if (status)
    status = product.layers.view(StateRole::accepted_n,
                                 product.fields.velocity,
                                 initial_velocity_with_ghosts);
  if (status)
    status = apply_boundary_ghosts(
        BoundaryStage::momentum, product.boundary,
        {&initial_velocity_with_ghosts, 1U},
        runtime.resolved_boundary_values());
  if (status)
    status = apply_boundary_ghosts(
        BoundaryStage::pressure, product.boundary,
        {&initial_pressure_with_ghosts, 1U},
        runtime.resolved_boundary_values());
  if (status)
    status = apply_boundary_ghosts(
        BoundaryStage::enthalpy, product.boundary,
        {&initial_enthalpy_with_ghosts, 1U},
        runtime.resolved_boundary_values());
  if (status && !product.fields.scalars.empty()) {
    species_index = 0U;
    passive_index = 0U;
    for (std::size_t scalar = 0U;
         scalar < product.fields.scalars.size() && status; ++scalar) {
      status = product.layers.view(StateRole::accepted_n,
                                   product.fields.scalars[scalar],
                                   runtime.halo_views[scalar]);
      if (!status) break;
      if (product.fields.scalar_roles[scalar] ==
          TransportedScalarRole::species) {
        runtime.species_accepted[species_index++] =
            as_const(runtime.halo_views[scalar]);
      } else {
        runtime.passive_accepted[passive_index++] =
            as_const(runtime.halo_views[scalar]);
      }
    }
    if (status)
      status = apply_boundary_ghosts(
          BoundaryStage::scalar, product.boundary,
          {runtime.halo_views.data(), product.fields.scalars.size()},
          runtime.resolved_boundary_values());
    if (status) {
      species_index = 0U;
      passive_index = 0U;
      for (std::size_t scalar = 0U; scalar < product.fields.scalars.size();
           ++scalar) {
        if (product.fields.scalar_roles[scalar] ==
            TransportedScalarRole::species)
          runtime.species_accepted[species_index++] =
              as_const(runtime.halo_views[scalar]);
        else
          runtime.passive_accepted[passive_index++] =
              as_const(runtime.halo_views[scalar]);
      }
    }
  }
  FieldView initial_heat_capacity;
  FieldView initial_compressibility;
  FieldView initial_enthalpy_compressibility;
  FieldView initial_molecular_viscosity;
  FieldView initial_conductivity;
  FieldView initial_enthalpy_diffusivity;
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace, product.fields.heat_capacity,
        initial_heat_capacity);
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace, product.fields.compressibility,
        initial_compressibility);
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace,
        product.fields.enthalpy_compressibility,
        initial_enthalpy_compressibility);
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace,
        product.fields.molecular_viscosity, initial_molecular_viscosity);
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace,
        product.fields.thermal_conductivity, initial_conductivity);
  if (status)
    status = product.layers.runtime_view(
        FieldLifetime::persistent_workspace,
        product.fields.enthalpy_diffusivity, initial_enthalpy_diffusivity);
  if (status && runtime.boundary_thermo_authority_fields.size() !=
                    runtime.species_accepted.size() + 2U)
    status = {StatusCode::invalid_plan, kProductBinding};
  if (status) {
    runtime.boundary_thermo_authority_fields[0U] =
        make_boundary_ghost_field_authority(
            as_const(initial_pressure_with_ghosts));
    runtime.boundary_thermo_authority_fields[1U] =
        make_boundary_ghost_field_authority(
            as_const(initial_enthalpy_with_ghosts));
    for (std::size_t species = 0U;
         species < runtime.species_accepted.size(); ++species)
      runtime.boundary_thermo_authority_fields[species + 2U] =
          make_boundary_ghost_field_authority(
              runtime.species_accepted[species]);
    const std::int32_t reach_value = std::max(
        {initial_pressure_with_ghosts.ghosts.x,
         initial_pressure_with_ghosts.ghosts.y,
         initial_pressure_with_ghosts.ghosts.z});
    if (reach_value <= 0 || reach_value > UINT8_MAX) {
      status = {StatusCode::invalid_plan, kProductBinding};
    } else {
      const BoundaryThermophysicalGhostAuthority authority{
          static_cast<std::uintptr_t>(product.fingerprint),
          product.boundary.revision(),
          product.boundary.local_layout_fingerprint(), product.patch.cells,
          initial_pressure_with_ghosts.ghosts,
          static_cast<std::uint8_t>(reach_value),
          {runtime.boundary_thermo_authority_fields.data(),
           runtime.boundary_thermo_authority_fields.size()}};
      status = BoundaryThermophysicalFaceClosure::close(
          product.boundary, product.thermodynamics, product.transport,
          {initial.pressure_reference,
           as_const(initial_pressure_with_ghosts),
           as_const(initial_enthalpy_with_ghosts),
           {runtime.species_accepted.data(),
            runtime.species_accepted.size()},
           authority},
          {initial_density_with_ghosts, initial_temperature_with_ghosts,
           initial_heat_capacity, initial_compressibility,
           initial_enthalpy_compressibility, initial_molecular_viscosity,
           initial_conductivity, initial_enthalpy_diffusivity});
    }
  }
  initial_density = as_const(initial_density_with_ghosts);
  initial_velocity = as_const(initial_velocity_with_ghosts);
  initial_enthalpy = as_const(initial_enthalpy_with_ghosts);
  status = product.reductions.consensus(status);
  if (!status) return status;

  FaceFluxView mechanical_flux;
  status = product.phi_workspace.workspace_view(0U, 1U, mechanical_flux);
  status = product.reductions.consensus(status);
  if (!status) return status;
  for (std::int32_t z = 0; z < mechanical_flux.x.extents.z; ++z)
    for (std::int32_t y = 0; y < mechanical_flux.x.extents.y; ++y)
      for (std::int32_t x = 0; x < mechanical_flux.x.extents.x; ++x)
        mechanical_flux.x.unchecked({x, y, z}) =
            thermo.rho * initial.velocity.x *
            dy.data[static_cast<std::size_t>(product.patch.begin.y + y)] *
            dz.data[static_cast<std::size_t>(product.patch.begin.z + z)];
  for (std::int32_t z = 0; z < mechanical_flux.y.extents.z; ++z)
    for (std::int32_t y = 0; y < mechanical_flux.y.extents.y; ++y)
      for (std::int32_t x = 0; x < mechanical_flux.y.extents.x; ++x)
        mechanical_flux.y.unchecked({x, y, z}) =
            thermo.rho * initial.velocity.y *
            dx.data[static_cast<std::size_t>(product.patch.begin.x + x)] *
            dz.data[static_cast<std::size_t>(product.patch.begin.z + z)];
  for (std::int32_t z = 0; z < mechanical_flux.z.extents.z; ++z)
    for (std::int32_t y = 0; y < mechanical_flux.z.extents.y; ++y)
      for (std::int32_t x = 0; x < mechanical_flux.z.extents.x; ++x)
        mechanical_flux.z.unchecked({x, y, z}) =
            thermo.rho * initial.velocity.z *
            dx.data[static_cast<std::size_t>(product.patch.begin.x + x)] *
            dy.data[static_cast<std::size_t>(product.patch.begin.y + y)];
  const std::array<FaceFieldView, 3U> mechanical_flux_fields{
      mechanical_flux.x, mechanical_flux.y, mechanical_flux.z};
  for (std::size_t face_index = 0U; face_index < 6U && status;
       ++face_index) {
    const auto face = static_cast<CartesianFace>(face_index);
    const BoundaryFacePlan* face_plan = nullptr;
    if (!product.boundary.face(face, face_plan) || face_plan == nullptr) {
      status = {StatusCode::invalid_plan, kProductBinding};
      break;
    }
    if (!face_plan->local_owner || face_plan->periodic) continue;
    const auto axis = face_index < 2U
                          ? CartesianAxis::x
                          : (face_index < 4U ? CartesianAxis::y
                                             : CartesianAxis::z);
    const bool high = (face_index & 1U) != 0U;
    const std::int32_t normal_extent =
        axis == CartesianAxis::x
            ? product.patch.cells.x
            : (axis == CartesianAxis::y ? product.patch.cells.y
                                         : product.patch.cells.z);
    const std::int32_t normal = high ? normal_extent : 0;
    const std::int32_t inner_count =
        axis == CartesianAxis::x ? product.patch.cells.y
                                 : product.patch.cells.x;
    const std::int32_t outer_count =
        axis == CartesianAxis::z ? product.patch.cells.y
                                 : product.patch.cells.z;
    const std::uint8_t component =
        static_cast<std::uint8_t>(face_index / 2U);
    FaceFieldView flux = mechanical_flux_fields[component];
    for (std::int32_t outer = 0; outer < outer_count; ++outer) {
      for (std::int32_t inner = 0; inner < inner_count; ++inner) {
        const Int3 owner =
            boundary_owner_cell(face, product.patch.cells, inner, outer);
        Int3 ghost = owner;
        Int3 face_cell{};
        if (axis == CartesianAxis::x) {
          ghost.x = high ? product.patch.cells.x : -1;
          face_cell = {normal, inner, outer};
        } else if (axis == CartesianAxis::y) {
          ghost.y = high ? product.patch.cells.y : -1;
          face_cell = {inner, normal, outer};
        } else {
          ghost.z = high ? product.patch.cells.z : -1;
          face_cell = {inner, outer, normal};
        }
        const Int3 left = high ? owner : ghost;
        const Int3 right = high ? ghost : owner;
        const double velocity_face = detail::interpolate_face(
            product.equations.kernels(), axis, normal,
            initial_velocity_with_ghosts.unchecked(left, component),
            initial_velocity_with_ghosts.unchecked(right, component));
        const double value =
            thermo.rho * velocity_face *
            detail::face_area(product.equations.kernels(), axis, face_cell);
        if (!std::isfinite(velocity_face) || !std::isfinite(value)) {
          status = {StatusCode::numerical_failure, kProductBinding};
          break;
        }
        flux.unchecked(face_cell) = value;
      }
      if (!status) break;
    }
  }
  FaceFluxView initial_flux = mechanical_flux;
  const bool fresh_open_physical_flux =
      product.boundary.pressure_reference() ==
      PressureReferenceKind::boundary_absolute;
  if (status && fresh_open_physical_flux)
    status = product.phi_workspace.workspace_view(1U, 1U, initial_flux);
  if (status && fresh_open_physical_flux) {
    status = product.candidate_boundary_finalizer.ready()
                 ? product.candidate_boundary_finalizer
                       .close_fresh_physical_flux(
                           {initial.pressure_reference,
                            as_const(initial_pressure_with_ghosts),
                            as_const(initial_velocity_with_ghosts),
                            {runtime.species_accepted.data(),
                             runtime.species_accepted.size()},
                            as_const(mechanical_flux), initial_flux},
                           product.reductions)
                 : Status{StatusCode::invalid_plan, kProductBinding};
  }
  if (status && !fresh_open_physical_flux &&
      product.ibm_equations.has_value())
    status = product.ibm_equations->zero_interface_flux(initial_flux);
  // A fresh IBM state is not generally compatible with the conservative face
  // flux assembled above.  Project U and phi as one certified transaction
  // before the final-flux writer can publish either value as accepted history.
  // The projection API deliberately has no pressure destination.
  std::array<FieldView, 3U> initial_velocity_layers{};
  if (status) {
    const std::array<StateRole, 3U> roles{{
        StateRole::accepted_n, StateRole::accepted_n_minus_one,
        StateRole::trial}};
    for (std::size_t role = 0U; role < roles.size() && status; ++role)
      status = product.layers.view(roles[role], product.fields.velocity,
                                   initial_velocity_layers[role]);
  }
  if (status && !product.fresh_projection.has_value())
    status = {StatusCode::invalid_plan, kProductBinding};
  // Flux construction, IBM interface closure, and layer lookup are all
  // local.  Rendezvous once more before prepare(), whose first operation is
  // collective, so every rank either enters it or returns together.
  status = product.reductions.consensus(status);
  if (!status) return status;
  FreshStartKinematicProjectionPreparedCertificate fresh_prepared;
  FreshStartKinematicProjectionSolvedCertificate fresh_solved;
  FreshStartKinematicProjectionCandidateCertificate fresh_candidate;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  detail::FreshInitializationDiagnostic fresh_diagnostic;
  fresh_diagnostic.valid = true;
  fresh_diagnostic.driver_plan = product.fingerprint;
  fresh_diagnostic.generation =
      g_fresh_initialization_generation.fetch_add(
          RevisionToken{1U}, std::memory_order_acq_rel) +
      RevisionToken{1U};
  fresh_diagnostic.immersed = product.topology.has_value();
  fresh_diagnostic.projection_attempted = product.topology.has_value();
  fresh_diagnostic.no_ibm_bypassed = !product.topology.has_value();
  if (product.fresh_projection.has_value())
    fresh_diagnostic.red_plan = product.fresh_projection->red().plan;
#endif
  if (status) {
    RevisionToken fresh_state = detail::product_mix(
        product.fingerprint, initial_density.revision);
    for (const FieldView layer : initial_velocity_layers)
      fresh_state = detail::product_mix(fresh_state, layer.revision);
    fresh_state = detail::product_mix(fresh_state, initial_flux.revision);
    if (fresh_state == 0U) fresh_state = 1U;
    const FreshStartKinematicProjectionInput fresh_input{
        initial_density,
        as_const(initial_velocity_layers[0U]),
        as_const(initial_velocity_layers[1U]),
        as_const(initial_velocity_layers[2U]),
        as_const(initial_flux),
        fresh_state};
    status = product.fresh_projection->prepare(fresh_input, fresh_prepared);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    if (status)
      fresh_diagnostic.prepared_lineage = fresh_prepared.lineage();
#endif
  }
  if (status)
    status = product.fresh_projection->solve(fresh_prepared, fresh_solved,
                                             &runtime.resources);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (status) {
    fresh_diagnostic.solve = fresh_solved.result();
    fresh_diagnostic.solved_lineage = fresh_solved.lineage();
  }
#endif
  if (status)
    status = product.fresh_projection->audit(fresh_solved, fresh_candidate);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (status) {
    fresh_diagnostic.audited = true;
    fresh_diagnostic.candidate_lineage = fresh_candidate.lineage();
    fresh_diagnostic.initial_continuity_maximum =
        fresh_candidate.initial_continuity_maximum();
    fresh_diagnostic.final_continuity_maximum =
        fresh_candidate.final_continuity_maximum();
    fresh_diagnostic.continuity_limit =
        1.0e-12 + 1.0e-10 *
                       fresh_candidate.initial_continuity_maximum();
  }
  if (status && g_fresh_initialization_poison_armed.exchange(
                    false, std::memory_order_acq_rel)) {
    int rank = -1;
    if (MPI_Comm_rank(runtime.communicator, &rank) != MPI_SUCCESS) {
      status = {StatusCode::mpi_failure, kProductCollective};
    } else if (rank == g_fresh_initialization_poison_rank.load(
                           std::memory_order_acquire)) {
      const auto poison_kind = static_cast<detail::FreshInitializationPoisonKind>(
          g_fresh_initialization_poison_kind.load(std::memory_order_acquire));
      if (poison_kind == detail::FreshInitializationPoisonKind::velocity) {
        FieldView poisoned;
        status = product.layers.runtime_view(
            FieldLifetime::persistent_workspace,
            product.fields.pressure_energy_candidate_velocity, poisoned);
        if (status) {
          double& value = poisoned.unchecked({0, 0, 0}, 0U);
          value = std::nextafter(value,
                                 std::numeric_limits<double>::infinity());
        }
      } else {
        FaceFluxView poisoned;
        status = product.phi_workspace.workspace_view(
            2U, initial_flux.revision, poisoned);
        if (status) {
          double& value = poisoned.x.unchecked({0, 0, 0});
          value = std::nextafter(value,
                                 std::numeric_limits<double>::infinity());
        }
      }
    }
  }
#endif
  status = product.reductions.consensus(status);
  if (status)
    status = product.fresh_projection->commit(
        fresh_candidate,
        {initial_velocity_layers.data(), initial_velocity_layers.size()},
        initial_flux);
  // Projection may update pressure-outlet faces and their owner velocity.
  // Replay the exact same physical closure against that committed U/phi
  // candidate and require a bit-stable fixed point before either becomes
  // accepted history.  The replay uses replica 3; projection owns replica 2.
  if (status && fresh_open_physical_flux) {
    FaceFluxView replayed_flux;
    status = product.phi_workspace.workspace_view(
        3U, initial_flux.revision, replayed_flux);
    if (status)
      status = product.candidate_boundary_finalizer
                   .close_fresh_physical_flux(
                       {initial.pressure_reference,
                        as_const(initial_pressure_with_ghosts),
                        as_const(initial_velocity_layers[0U]),
                        {runtime.species_accepted.data(),
                         runtime.species_accepted.size()},
                        as_const(initial_flux), replayed_flux},
                       product.reductions);
    bool bitwise_sealed = static_cast<bool>(status);
    if (status) {
      const std::array<ConstFaceFieldView, 3U> accepted_faces{
          as_const(initial_flux.x), as_const(initial_flux.y),
          as_const(initial_flux.z)};
      const std::array<ConstFaceFieldView, 3U> replayed_faces{
          as_const(replayed_flux.x), as_const(replayed_flux.y),
          as_const(replayed_flux.z)};
      for (std::size_t axis = 0U;
           axis < accepted_faces.size() && bitwise_sealed; ++axis)
        for (std::int32_t z = 0;
             z < accepted_faces[axis].extents.z && bitwise_sealed; ++z)
          for (std::int32_t y = 0;
               y < accepted_faces[axis].extents.y && bitwise_sealed; ++y)
            for (std::int32_t x = 0;
                 x < accepted_faces[axis].extents.x; ++x)
              if (std::memcmp(
                      &accepted_faces[axis].unchecked({x, y, z}),
                      &replayed_faces[axis].unchecked({x, y, z}),
                      sizeof(double)) != 0) {
                bitwise_sealed = false;
                break;
              }
    }
    status = product.reductions.consensus(
        status && bitwise_sealed
            ? Status{}
            : (status ? Status{StatusCode::rejected_step, kProductBinding}
                      : status));
  }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (status) fresh_diagnostic.committed = true;
#endif
  // The empty-activity route is a literal no-store bypass.  The original
  // uniform cold histories are already kinematically compatible, so even
  // derived-workspace revisions must remain untouched in this branch.
  if (status && product.topology.has_value())
    status = runtime.rebuild_cold_velocity_dependents(
        initial.pressure_reference, false);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (status && g_cold_velocity_dependents_published.load(
                    std::memory_order_acquire)) {
    const detail::ColdVelocityDependentsDiagnostic& derived =
        g_cold_velocity_dependents_diagnostic;
    if (derived.valid && !derived.restart &&
        derived.driver_plan == product.fingerprint) {
      fresh_diagnostic.derived_velocity_dependents_rebuilt = derived.rebuilt;
      fresh_diagnostic.derived_velocity_lineage = derived.lineage;
      fresh_diagnostic.maximum_h_by_a_velocity_difference =
          derived.maximum_h_by_a_velocity_difference;
      fresh_diagnostic.velocity_gradient_nonfinite_count =
          derived.velocity_gradient_nonfinite_count;
      fresh_diagnostic.effective_viscosity_nonpositive_count =
          derived.effective_viscosity_nonpositive_count;
      fresh_diagnostic.solid_velocity_gradient_nonfinite_count =
          derived.solid_velocity_gradient_nonfinite_count;
      fresh_diagnostic.velocity_dependent_rate_layers_bitwise_equal =
          derived.rate_layers_bitwise_equal;
    }
  }
#endif
  if (status)
    status = runtime.final_flux_writer.initialize_committed(
        product.final_flux, as_const(initial_flux));
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  {
    double local_face_envelope = 0.0;
    const Span<const std::uint8_t> diagnostic_activity =
        product.topology.has_value() ? product.topology->region()
                                     : Span<const std::uint8_t>{};
    for (std::int32_t z = 0; z < product.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < product.patch.cells.y; ++y)
        for (std::int32_t x = 0; x < product.patch.cells.x; ++x) {
          const std::size_t flat =
              static_cast<std::size_t>(x) +
              static_cast<std::size_t>(product.patch.cells.x) *
                  (static_cast<std::size_t>(y) +
                   static_cast<std::size_t>(product.patch.cells.y) *
                       static_cast<std::size_t>(z));
          if (diagnostic_activity.size != 0U &&
              diagnostic_activity.data[flat] == 0U)
            continue;
          const double envelope =
              std::abs(initial_flux.x.unchecked({x, y, z})) +
              std::abs(initial_flux.x.unchecked({x + 1, y, z})) +
              std::abs(initial_flux.y.unchecked({x, y, z})) +
              std::abs(initial_flux.y.unchecked({x, y + 1, z})) +
              std::abs(initial_flux.z.unchecked({x, y, z})) +
              std::abs(initial_flux.z.unchecked({x, y, z + 1}));
          local_face_envelope = std::max(local_face_envelope, envelope);
        }
    double global_face_envelope = 0.0;
    if (MPI_Allreduce(&local_face_envelope, &global_face_envelope, 1,
                      MPI_DOUBLE, MPI_MAX, runtime.communicator) == MPI_SUCCESS)
      fresh_diagnostic.maximum_face_envelope = global_face_envelope;

    std::uint64_t local_cut_nonzero = 0U;
    std::uint64_t local_cut_negative_zero = 0U;
    if (product.topology.has_value()) {
      const Span<const ImmersedLink> links = product.topology->links();
      for (std::size_t link_index = 0U; link_index < links.size; ++link_index) {
        const ImmersedLink& link = links.data[link_index];
        Int3 face = link.fluid_local_index;
        const FaceFieldView* selected = &initial_flux.x;
        switch (link.direction) {
          case ImmersedFaceDirection::x_negative:
            selected = &initial_flux.x;
            break;
          case ImmersedFaceDirection::x_positive:
            selected = &initial_flux.x;
            ++face.x;
            break;
          case ImmersedFaceDirection::y_negative:
            selected = &initial_flux.y;
            break;
          case ImmersedFaceDirection::y_positive:
            selected = &initial_flux.y;
            ++face.y;
            break;
          case ImmersedFaceDirection::z_negative:
            selected = &initial_flux.z;
            break;
          case ImmersedFaceDirection::z_positive:
            selected = &initial_flux.z;
            ++face.z;
            break;
        }
        const double value = selected->unchecked(face);
        // Both IEEE zero signs are physically zero flux.  Keep the sign count
        // separately so the diagnostic can still certify the wire image
        // without misclassifying a restored -0.0 face as leakage.
        if (value != 0.0) ++local_cut_nonzero;
        if (value == 0.0 && std::signbit(value))
          ++local_cut_negative_zero;
      }
    }
    std::array<std::uint64_t, 2U> local_cut{
        local_cut_nonzero, local_cut_negative_zero};
    std::array<std::uint64_t, 2U> global_cut{};
    if (MPI_Allreduce(local_cut.data(), global_cut.data(),
                      static_cast<int>(local_cut.size()), MPI_UINT64_T,
                      MPI_SUM, runtime.communicator) == MPI_SUCCESS) {
      fresh_diagnostic.cut_face_nonzero_count = global_cut[0U];
      fresh_diagnostic.cut_face_negative_zero_count = global_cut[1U];
    }

    std::uint64_t local_solid_nonzero = 0U;
    bool local_layers_equal = true;
    for (std::int32_t z = 0; z < product.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < product.patch.cells.y; ++y)
        for (std::int32_t x = 0; x < product.patch.cells.x; ++x) {
          const std::size_t flat =
              static_cast<std::size_t>(x) +
              static_cast<std::size_t>(product.patch.cells.x) *
                  (static_cast<std::size_t>(y) +
                   static_cast<std::size_t>(product.patch.cells.y) *
                       static_cast<std::size_t>(z));
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            const std::uint64_t accepted_bits = product_double_bits(
                initial_velocity_layers[0U].unchecked({x, y, z}, component));
            for (std::size_t layer = 1U;
                 layer < initial_velocity_layers.size(); ++layer)
              local_layers_equal =
                  local_layers_equal &&
                  accepted_bits == product_double_bits(
                                       initial_velocity_layers[layer].unchecked(
                                           {x, y, z}, component));
            if (diagnostic_activity.size != 0U &&
                diagnostic_activity.data[flat] == 0U &&
                accepted_bits != product_double_bits(0.0))
              ++local_solid_nonzero;
          }
        }
    std::uint64_t global_solid_nonzero = 0U;
    int local_equal = local_layers_equal ? 1 : 0;
    int global_equal = 0;
    if (MPI_Allreduce(&local_solid_nonzero, &global_solid_nonzero, 1,
                      MPI_UINT64_T, MPI_SUM,
                      runtime.communicator) == MPI_SUCCESS)
      fresh_diagnostic.solid_nonpositive_zero_component_count =
          global_solid_nonzero;
    if (MPI_Allreduce(&local_equal, &global_equal, 1, MPI_INT, MPI_MIN,
                      runtime.communicator) == MPI_SUCCESS)
      fresh_diagnostic.velocity_layers_bitwise_equal = global_equal != 0;
  }
  fresh_diagnostic.terminal_status = status;
  publish_fresh_initialization_diagnostic(fresh_diagnostic);
#endif
  if (!status) return status;
  runtime.time = std::move(controller);
  runtime.pressure_reference = initial.pressure_reference;
  runtime.previous_pressure_reference = initial.pressure_reference;
  runtime.closed_mass_target = global_mass;
  runtime.pressure_correction_warm_start_valid = false;
  runtime.enthalpy_ghosts = {};
  std::fill(runtime.species_ghosts.begin(), runtime.species_ghosts.end(),
            ThermophysicalGhostHistory{});
  std::fill(runtime.passive_ghosts.begin(), runtime.passive_ghosts.end(),
            ThermophysicalGhostHistory{});
  runtime.predictor_diagnostics = {};
  runtime.initialized = true;
  return {};
}

Status ProductDriver::initialize_restart(const RestartImage& image) noexcept
    try {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  g_cold_velocity_dependents_published.store(false,
                                             std::memory_order_release);
  g_cold_velocity_dependents_diagnostic = {};
  g_fresh_initialization_published.store(false, std::memory_order_release);
  g_fresh_initialization_diagnostic = {};
#endif
  if (implementation_ == nullptr)
    return {StatusCode::invalid_plan, kProductInput};
  ProductDriver::Impl& runtime = *implementation_;
  CompiledCasePlan::Impl& product = *runtime.plan.implementation_;
  const bool exact_history = !image.backward_euler_recovery;
  Status status =
      !runtime.initialized &&
              std::isfinite(image.time) && std::isfinite(image.dt) &&
              image.dt > 0.0 && std::isfinite(image.pressure_reference) &&
              image.pressure_reference > 0.0 && image.step != 0U &&
              (!exact_history ||
               (image.controller_state != 0U &&
                std::isfinite(image.previous_pressure_reference) &&
                image.previous_pressure_reference > 0.0 &&
                std::isfinite(image.closed_mass_target) &&
                image.closed_mass_target > 0.0 &&
                image.previous_mass_flux_revision != 0U &&
                image.previous_mass_flux_revision <
                    image.final_mass_flux_revision &&
                image.final_mass_flux_revision !=
                    std::numeric_limits<RevisionToken>::max()))
          ? Status{}
          : Status{StatusCode::invalid_plan, kProductInput};
  const auto same_int3 = [](Int3 left, Int3 right) noexcept {
    return left.x == right.x && left.y == right.y && left.z == right.z;
  };
  const auto same_patch = [&](MeshPatch left, MeshPatch right) noexcept {
    return same_int3(left.begin, right.begin) &&
           same_int3(left.cells, right.cells) &&
           same_int3(left.process_grid, right.process_grid) &&
           same_int3(left.process_coord, right.process_coord);
  };
  if (status &&
      (!same_int3(image.global_cells, product.geometry.global_cells()) ||
       !same_patch(image.patch, product.patch) ||
       image.plan != runtime.plan.fingerprint() ||
       image.schema != product.schema_fingerprint ||
       image.geometry != product.geometry.fingerprint() ||
       image.fields.size() != runtime.restart_expected_fields.size()))
    status = {StatusCode::invalid_plan, kProductInput};
  const std::size_t local_cells =
      static_cast<std::size_t>(product.patch.cells.x) *
      static_cast<std::size_t>(product.patch.cells.y) *
      static_cast<std::size_t>(product.patch.cells.z);
  const auto validate_fields = [&](const std::vector<RestartImageField>& fields,
                                   const std::vector<RestartExpectedField>&
                                       expected_fields) {
    if (fields.size() != expected_fields.size())
      return Status{StatusCode::invalid_plan, kProductInput};
    for (std::size_t index = 0U; index < fields.size(); ++index) {
      const RestartImageField& field = fields[index];
      const RestartExpectedField expected = expected_fields[index];
      if (field.role != expected.role || field.field != expected.field ||
          field.components != expected.components ||
          field.values.size() != local_cells * field.components)
        return Status{StatusCode::invalid_plan, kProductInput};
      for (double value : field.values)
        if (!std::isfinite(value))
          return Status{StatusCode::numerical_failure, kProductInput};
    }
    return Status{};
  };
  if (status)
    status = validate_fields(image.fields, runtime.restart_expected_fields);
  if (status && exact_history)
    status = validate_fields(image.previous_fields,
                             runtime.restart_expected_fields);
  if (status && exact_history)
    status = validate_fields(image.accepted_rate_fields,
                             runtime.restart_expected_rate_fields);
  if (status && exact_history)
    status = validate_fields(image.previous_rate_fields,
                             runtime.restart_expected_rate_fields);
  const Int3 cells = product.patch.cells;
  const std::size_t x_faces =
      static_cast<std::size_t>(cells.x + 1) * cells.y * cells.z;
  const std::size_t y_faces =
      static_cast<std::size_t>(cells.x) * (cells.y + 1) * cells.z;
  const std::size_t z_faces =
      static_cast<std::size_t>(cells.x) * cells.y * (cells.z + 1);
  if (status &&
      (image.final_mass_flux[0U].size() != x_faces ||
       image.final_mass_flux[1U].size() != y_faces ||
       image.final_mass_flux[2U].size() != z_faces ||
       (exact_history &&
        (image.previous_mass_flux[0U].size() != x_faces ||
         image.previous_mass_flux[1U].size() != y_faces ||
         image.previous_mass_flux[2U].size() != z_faces))))
    status = {StatusCode::invalid_plan, kProductInput};
  for (const std::vector<double>& face : image.final_mass_flux)
    for (double value : face)
      if (status && !std::isfinite(value))
        status = {StatusCode::numerical_failure, kProductInput};
  if (exact_history)
    for (const std::vector<double>& face : image.previous_mass_flux)
      for (double value : face)
        if (status && !std::isfinite(value))
          status = {StatusCode::numerical_failure, kProductInput};
  // All malformed-checkpoint decisions become collective before any rank
  // constructs a view into the checkpoint or enters a later mass reduction.
  status = product.reductions.consensus(status);
  if (!status) return status;
  TimeControllerState controller;
  status = exact_history
               ? TimeControllerState::restart_exact(
                     product.time, image.time, image.dt, image.step,
                     image.controller_state, controller)
               : TimeControllerState::restart(
                     product.time, image.time, image.dt, image.step,
                     controller);
  status = product.reductions.consensus(status);
  if (!status) return status;
  const auto face = [](const std::vector<double>& values, Int3 extents,
                       CartesianAxis axis, StorageIdentity identity) {
    ConstFaceFieldView view;
    view.base = values.data();
    view.extents = extents;
    view.stride_y = static_cast<std::size_t>(extents.x);
    view.stride_z = view.stride_y * extents.y;
    view.axis = axis;
    view.storage_identity = identity;
    view.revision_domain = identity ^ UINT64_C(0x4f04);
    return view;
  };
  const RevisionToken flux_revision =
      exact_history ? image.final_mass_flux_revision : RevisionToken{17U};
  ConstFaceFluxView restored_flux{
      face(image.final_mass_flux[0U], {cells.x + 1, cells.y, cells.z},
           CartesianAxis::x, 401U),
      face(image.final_mass_flux[1U], {cells.x, cells.y + 1, cells.z},
           CartesianAxis::y, 401U),
      face(image.final_mass_flux[2U], {cells.x, cells.y, cells.z + 1},
           CartesianAxis::z, 401U),
      flux_revision,
      {}};
  ConstFaceFluxView restored_previous_flux{
      face(image.previous_mass_flux[0U], {cells.x + 1, cells.y, cells.z},
           CartesianAxis::x, 402U),
      face(image.previous_mass_flux[1U], {cells.x, cells.y + 1, cells.z},
           CartesianAxis::y, 402U),
      face(image.previous_mass_flux[2U], {cells.x, cells.y, cells.z + 1},
           CartesianAxis::z, 402U),
      image.previous_mass_flux_revision,
      {}};
  // Restart owns a shared face from one patch, while an immersed link can
  // belong to its neighbor.  Restore the payload first, then reapply the
  // mandatory zero-flux constraint to both accepted history levels.
  const auto canonicalize_ibm_flux =
      [&](ConstFaceFluxView source, std::size_t replica,
          ConstFaceFluxView& out) noexcept {
        FaceFluxView destination;
        Status local = product.phi_workspace.workspace_view(
            replica, source.revision, destination);
        if (!local) return local;
        const std::array<ConstFaceFieldView, 3U> inputs{
            source.x, source.y, source.z};
        const std::array<FaceFieldView, 3U> outputs{
            destination.x, destination.y, destination.z};
        for (std::size_t axis = 0U; axis < inputs.size(); ++axis)
          for (std::int32_t z = 0; z < inputs[axis].extents.z; ++z)
            for (std::int32_t y = 0; y < inputs[axis].extents.y; ++y)
              for (std::int32_t x = 0; x < inputs[axis].extents.x; ++x)
                outputs[axis].unchecked({x, y, z}) =
                    inputs[axis].unchecked({x, y, z});
        local = product.ibm_equations->zero_interface_flux(destination);
        if (local) out = as_const(destination);
        return local;
      };
  if (product.ibm_equations.has_value())
    status = canonicalize_ibm_flux(restored_flux, 0U, restored_flux);
  if (status && exact_history && product.ibm_equations.has_value())
    status = canonicalize_ibm_flux(restored_previous_flux, 1U,
                                   restored_previous_flux);
  if (status && product.ibm_equations.has_value())
    status = product.ibm_equations->validate_interface_flux(restored_flux);
  if (status && exact_history && product.ibm_equations.has_value())
    status =
        product.ibm_equations->validate_interface_flux(restored_previous_flux);
  status = product.reductions.consensus(status);
  if (!status) return status;

  const RestartImageField& velocity = image.fields[0U];
  const RestartImageField& pressure = image.fields[1U];
  const RestartImageField& enthalpy = image.fields[2U];
  const Span<const double> dx = product.geometry.x().widths();
  const Span<const double> dy = product.geometry.y().widths();
  const Span<const double> dz = product.geometry.z().widths();
  // rho, T, mu, lambda, cp, drho/dp, drho/dh, lambda/cp. This cold staging makes
  // thermodynamic validation atomic with respect to the driver state.
  std::vector<std::array<double, 8U>> restored_thermo;
  std::vector<std::array<double, 8U>> restored_previous_thermo;
  Status allocation_status;
  try {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    int restart_rank = -1;
    (void)MPI_Comm_rank(runtime.communicator, &restart_rank);
    if (consume_restart_restore_allocation_failure(
            detail::RestartRestoreAllocationPoint::thermophysical_staging,
            restart_rank))
      throw std::bad_alloc{};
#endif
    restored_thermo.resize(local_cells);
    if (exact_history) restored_previous_thermo.resize(local_cells);
  } catch (const std::bad_alloc&) {
    allocation_status = {StatusCode::allocation_failure, kProductAllocation};
  } catch (...) {
    allocation_status = {StatusCode::invalid_plan, kProductInput};
  }
  status = product.reductions.consensus(allocation_status);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  g_restart_restore_allocation_lowest_failing_rank.store(
      product.reductions.lowest_failing_rank(), std::memory_order_release);
#endif
  if (!status) return status;
  double local_mass = 0.0;
  const auto evaluate_level =
      [&](const std::vector<RestartImageField>& fields,
          double reference_pressure,
          std::vector<std::array<double, 8U>>& states,
          bool accumulate_mass) {
        const RestartImageField& level_velocity = fields[0U];
        const RestartImageField& level_pressure = fields[1U];
        const RestartImageField& level_enthalpy = fields[2U];
        Status level_status;
        for (std::int32_t z = 0; z < cells.z && level_status; ++z)
          for (std::int32_t y = 0; y < cells.y && level_status; ++y)
            for (std::int32_t x = 0; x < cells.x; ++x) {
              const std::size_t cell =
                  static_cast<std::size_t>(x) +
                  static_cast<std::size_t>(cells.x) *
                      (static_cast<std::size_t>(y) +
                       static_cast<std::size_t>(cells.y) * z);
              std::size_t species = 0U;
              for (std::size_t scalar = 0U;
                   scalar < product.fields.scalars.size(); ++scalar)
                if (product.fields.scalar_roles[scalar] ==
                    TransportedScalarRole::species)
                  runtime.species_values[species++] =
                      fields[scalar + 3U].values[cell];
              const Real3 value{
                  level_velocity.values[cell * 3U],
                  level_velocity.values[cell * 3U + 1U],
                  level_velocity.values[cell * 3U + 2U]};
              ThermoState thermo;
              level_status =
                  product.thermodynamics.evaluate_from_reference_pressure(
                      reference_pressure, level_pressure.values[cell],
                      level_enthalpy.values[cell],
                      {runtime.species_values.data(),
                       runtime.species_values.size()},
                      value, thermo, 300.0);
              MolecularTransportState transport;
              if (level_status)
                level_status = product.transport.evaluate(
                    thermo.temperature,
                    {runtime.species_values.data(),
                     runtime.species_values.size()},
                    transport);
              if (!level_status || !(thermo.rho > 0.0) ||
                  !(thermo.cp > 0.0) || !(thermo.drho_dp_hY > 0.0) ||
                  !(transport.viscosity > 0.0) ||
                  !(thermo.drho_dh_pY < 0.0) ||
                  !(transport.conductivity > 0.0)) {
                if (level_status)
                  level_status = {StatusCode::numerical_failure,
                                  kProductInput};
                break;
              }
              states[cell] = {
                  thermo.rho, thermo.temperature, transport.viscosity,
                  transport.conductivity, thermo.cp, thermo.drho_dp_hY,
                  thermo.drho_dh_pY, transport.conductivity / thermo.cp};
              if (accumulate_mass &&
                  (!product.topology.has_value() ||
                   product.topology->region().data[cell] != 0U)) {
                const std::size_t gx = static_cast<std::size_t>(
                    product.patch.begin.x + x);
                const std::size_t gy = static_cast<std::size_t>(
                    product.patch.begin.y + y);
                const std::size_t gz = static_cast<std::size_t>(
                    product.patch.begin.z + z);
                local_mass +=
                    thermo.rho * dx.data[gx] * dy.data[gy] * dz.data[gz];
              }
            }
        return level_status;
      };
  status = evaluate_level(image.fields, image.pressure_reference,
                          restored_thermo, true);
  if (status && exact_history)
    status = evaluate_level(image.previous_fields,
                            image.previous_pressure_reference,
                            restored_previous_thermo, false);
  status = product.reductions.consensus(status);
  if (!status) return status;
  double global_mass = 0.0;
  status = product.reductions.checked_sum({&local_mass, 1U},
                                          {&global_mass, 1U});
  if (status && (!(global_mass > 0.0) || !std::isfinite(global_mass)))
    status = {StatusCode::numerical_failure, kProductInput};
  status = product.reductions.consensus(status);
  if (!status) return status;

  // Reset state-layer histories and persistent workspaces without entering
  // the fresh-start path.  In particular, a restart must not evaluate the
  // default 300 K seed, resolve fresh boundary fluxes, or run the future IBM
  // compatible-start projection.  The first restored cell is only a finite
  // ghost/workspace fill value; every restored interior authority is copied
  // from the image below before the flux lineage is published.
  std::vector<double> restored_scalars;
  allocation_status = {};
  try {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    int restart_rank = -1;
    (void)MPI_Comm_rank(runtime.communicator, &restart_rank);
    if (consume_restart_restore_allocation_failure(
            detail::RestartRestoreAllocationPoint::scalar_seed,
            restart_rank))
      throw std::bad_alloc{};
#endif
    restored_scalars.resize(product.fields.scalars.size(), 0.0);
  } catch (const std::bad_alloc&) {
    allocation_status = {StatusCode::allocation_failure, kProductAllocation};
  } catch (...) {
    allocation_status = {StatusCode::invalid_plan, kProductInput};
  }
  status = product.reductions.consensus(allocation_status);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  g_restart_restore_allocation_lowest_failing_rank.store(
      product.reductions.lowest_failing_rank(), std::memory_order_release);
#endif
  if (!status) return status;
  for (std::size_t scalar = 0U; scalar < restored_scalars.size(); ++scalar)
    restored_scalars[scalar] = image.fields[scalar + 3U].values[0U];
  DriverInitialState restored_seed;
  restored_seed.pressure_reference = image.pressure_reference;
  restored_seed.temperature = restored_thermo[0U][1U];
  restored_seed.velocity = {velocity.values[0U], velocity.values[1U],
                            velocity.values[2U]};
  restored_seed.transported_scalars = {restored_scalars.data(),
                                       restored_scalars.size()};
  restored_seed.start_time = image.time;
  ThermoState common_thermo;
  common_thermo.rho = restored_thermo[0U][0U];
  common_thermo.temperature = restored_thermo[0U][1U];
  common_thermo.cp = restored_thermo[0U][4U];
  common_thermo.drho_dp_hY = restored_thermo[0U][5U];
  common_thermo.drho_dh_pY = restored_thermo[0U][6U];
  MolecularTransportState common_transport;
  common_transport.viscosity = restored_thermo[0U][2U];
  common_transport.conductivity = restored_thermo[0U][3U];
  status = runtime.initialize_common_fields(
      restored_seed, enthalpy.values[0U], common_thermo, common_transport,
      restored_thermo[0U][4U]);
  status = product.reductions.consensus(status);
  if (!status) return status;

  const std::array<StateRole, 3U> roles{{
      StateRole::accepted_n, StateRole::accepted_n_minus_one,
      StateRole::trial}};
  const auto copy_image = [&](const RestartImageField& source,
                              FieldView destination) noexcept {
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const std::size_t cell = static_cast<std::size_t>(x) +
                                   static_cast<std::size_t>(cells.x) *
                                       (static_cast<std::size_t>(y) +
                                        static_cast<std::size_t>(cells.y) * z);
          for (std::uint8_t component = 0U;
               component < source.components; ++component)
            destination.unchecked({x, y, z}, component) =
                source.values[cell * source.components + component];
        }
  };
  for (std::size_t field_index = 0U; field_index < image.fields.size();
       ++field_index) {
    for (StateRole role : roles) {
      FieldView destination;
      status = product.layers.view(
          role, runtime.restart_expected_fields[field_index].field,
          destination);
      const RestartImageField& source =
          exact_history && role == StateRole::accepted_n_minus_one
              ? image.previous_fields[field_index]
              : image.fields[field_index];
      if (status) copy_image(source, destination);
    }
  }

  std::array<FieldView, 3U> rho_views;
  std::array<FieldView, 3U> temperature_views;
  for (std::size_t role = 0U; role < roles.size(); ++role) {
    status = product.layers.view(roles[role], product.fields.rho,
                                 rho_views[role]);
    if (status)
      status = product.layers.view(roles[role], product.fields.temperature,
                                   temperature_views[role]);
  }
  const auto runtime_view_for_write = [&](FieldId field, FieldView& view) {
    Status revised = product.layers.revise_runtime(
        FieldLifetime::persistent_workspace, field);
    return revised ? product.layers.runtime_view(
                         FieldLifetime::persistent_workspace, field, view)
                   : revised;
  };
  FieldView viscosity;
  FieldView effective_viscosity;
  FieldView conductivity;
  FieldView heat_capacity;
  FieldView compressibility;
  FieldView enthalpy_compressibility;
  FieldView enthalpy_diffusivity;
  FieldView eos_density;
  status = runtime_view_for_write(product.fields.molecular_viscosity,
                                  viscosity);
  if (status)
    status = runtime_view_for_write(product.fields.effective_viscosity,
                                    effective_viscosity);
  if (status)
    status = runtime_view_for_write(product.fields.thermal_conductivity,
                                    conductivity);
  if (status)
    status = runtime_view_for_write(product.fields.heat_capacity,
                                    heat_capacity);
  if (status)
    status = runtime_view_for_write(product.fields.compressibility,
                                    compressibility);
  if (status)
    status = runtime_view_for_write(product.fields.enthalpy_compressibility,
                                    enthalpy_compressibility);
  if (status)
    status = runtime_view_for_write(product.fields.enthalpy_diffusivity,
                                    enthalpy_diffusivity);
  if (status)
    status = runtime_view_for_write(product.fields.eos_density, eos_density);
  for (std::int32_t z = 0; z < cells.z && status; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const std::size_t cell = static_cast<std::size_t>(x) +
                                 static_cast<std::size_t>(cells.x) *
                                     (static_cast<std::size_t>(y) +
                                      static_cast<std::size_t>(cells.y) * z);
        const std::array<double, 8U>& restored = restored_thermo[cell];
        const std::array<double, 8U>& restored_previous =
            exact_history ? restored_previous_thermo[cell] : restored;
        const Int3 index{x, y, z};
        for (std::size_t role = 0U; role < roles.size(); ++role) {
          const std::array<double, 8U>& selected =
              exact_history &&
                      roles[role] == StateRole::accepted_n_minus_one
                  ? restored_previous
                  : restored;
          rho_views[role].unchecked(index, 0U) = selected[0U];
          temperature_views[role].unchecked(index, 0U) = selected[1U];
        }
        viscosity.unchecked(index, 0U) = restored[2U];
        effective_viscosity.unchecked(index, 0U) = restored[2U];
        conductivity.unchecked(index, 0U) = restored[3U];
        heat_capacity.unchecked(index, 0U) = restored[4U];
        compressibility.unchecked(index, 0U) = restored[5U];
        enthalpy_compressibility.unchecked(index, 0U) = restored[6U];
        enthalpy_diffusivity.unchecked(index, 0U) = restored[7U];
        eos_density.unchecked(index, 0U) = restored[0U];
      }

  // The cold-dependent rebuild immediately enters boundary/halo collectives.
  // No rank may skip it due to a local state-layer or runtime-view failure.
  status = product.reductions.consensus(status);
  if (!status) return status;

  // Restart bypasses the Fresh projection, but it must not retain the uniform
  // grad(U), turbulence, or rate placeholders installed by
  // initialize_common_fields.  Rebuild those derived histories from the
  // restored U/state while leaving checkpoint U and phi interiors untouched.
  status = runtime.rebuild_cold_velocity_dependents(
      image.pressure_reference, true);
  if (!status) return status;

  if (exact_history) {
    const auto restore_rates = [&](StateRole role,
                                   const std::vector<RestartImageField>&
                                       sources) {
      for (std::size_t rate = 0U; rate < sources.size(); ++rate) {
        FieldView destination;
        Status copied = product.layers.view(
            role, runtime.restart_expected_rate_fields[rate].field,
            destination);
        if (!copied) return copied;
        copy_image(sources[rate], destination);
      }
      return Status{};
    };
    status = restore_rates(StateRole::accepted_n,
                           image.accepted_rate_fields);
    if (status)
      status = restore_rates(StateRole::accepted_n_minus_one,
                             image.previous_rate_fields);
    if (status)
      status = restore_rates(StateRole::trial, image.accepted_rate_fields);
    status = product.reductions.consensus(status);
    if (!status) return status;
  }

  runtime.enthalpy_ghosts = {};
  std::fill(runtime.species_ghosts.begin(), runtime.species_ghosts.end(),
            ThermophysicalGhostHistory{});
  std::fill(runtime.passive_ghosts.begin(), runtime.passive_ghosts.end(),
            ThermophysicalGhostHistory{});
  if (exact_history) {
    std::size_t halo_count = 0U;
    const StateLayers& restored_layers = product.layers;
    ConstFieldView previous_density;
    ConstFieldView previous_velocity;
    status = restored_layers.view(StateRole::accepted_n_minus_one,
                                  product.fields.rho, previous_density);
    if (status)
      status = restored_layers.view(StateRole::accepted_n_minus_one,
                                    product.fields.velocity,
                                    previous_velocity);
    FieldView previous_enthalpy;
    if (status)
      status = product.layers.view(StateRole::accepted_n_minus_one,
                                   product.fields.enthalpy,
                                   previous_enthalpy);
    if (status) runtime.halo_views[halo_count++] = previous_enthalpy;
    for (std::size_t scalar = 0U;
         scalar < product.fields.scalars.size() && status; ++scalar) {
      FieldView previous_scalar;
      status = product.layers.view(StateRole::accepted_n_minus_one,
                                   product.fields.scalars[scalar],
                                   previous_scalar);
      if (status) runtime.halo_views[halo_count++] = previous_scalar;
    }
    status = product.reductions.consensus(status);
    if (!status) return status;
    HaloTicket ticket;
    status = product.stage_halos[0U].begin(
        10U, {runtime.halo_views.data(), halo_count}, ticket);
    if (status)
      status = product.stage_halos[0U].finish(
          ticket, {runtime.halo_views.data(), halo_count});
    if (status) {
      previous_enthalpy = runtime.halo_views[0U];
      std::size_t species = 0U;
      std::size_t passive = 0U;
      for (std::size_t scalar = 0U;
           scalar < product.fields.scalars.size(); ++scalar) {
        const ConstFieldView value =
            as_const(runtime.halo_views[scalar + 1U]);
        if (product.fields.scalar_roles[scalar] ==
            TransportedScalarRole::species)
          runtime.species_previous[species++] = value;
        else
          runtime.passive_previous[passive++] = value;
      }
      if (species != runtime.species_previous.size() ||
          passive != runtime.passive_previous.size())
        status = {StatusCode::invalid_plan, kProductBinding};
    }
    if (status)
      status = resolve_static_boundary_values(
          runtime.communicator, product.boundary, product.boundary_specs,
          product.geometry, product.patch, product.thermodynamics,
          image.previous_pressure_reference, previous_density,
          previous_velocity, as_const(previous_enthalpy),
          {runtime.species_previous.data(),
           runtime.species_previous.size()},
          {runtime.passive_previous.data(),
           runtime.passive_previous.size()},
          {runtime.species_values.data(), runtime.species_values.size()},
          runtime.boundary_scalars, runtime.boundary_vectors,
          runtime.boundary_normal_gradients, false);
    if (status) {
      status = apply_boundary_ghosts(
          BoundaryStage::enthalpy, product.boundary,
          {&previous_enthalpy, 1U}, runtime.resolved_boundary_values());
    }
    if (status && product.fields.scalars.size() != 0U)
      status = apply_boundary_ghosts(
          BoundaryStage::scalar, product.boundary,
          {runtime.halo_views.data() + 1U,
           product.fields.scalars.size()},
          runtime.resolved_boundary_values());
    status = product.reductions.consensus(status);
    if (!status) return status;

    const auto ghost_authority = [&](ConstFieldView view,
                                     std::uint8_t reach,
                                     ThermophysicalGhostAuthority& out) {
      if (product.stage_halos[0U].ghost_revision(view.field) !=
              view.revision ||
          product.stage_halos[0U].instance_identity() == 0U || reach == 0U ||
          view.ghosts.x < reach || view.ghosts.y < reach ||
          view.ghosts.z < reach)
        return Status{StatusCode::invalid_plan, kProductCommunication};
      out = {product.stage_halos[0U].instance_identity(),
             view.field,
             view.revision,
             view.storage_identity,
             view.revision_domain,
             product.geometry.topology_revision(),
             product.boundary.revision(),
             reach};
      return Status{};
    };
    status = ghost_authority(
        as_const(previous_enthalpy),
        product.equations.thermophysical_predictor().enthalpy_reach(),
        runtime.enthalpy_ghosts.previous);
    std::size_t species = 0U;
    std::size_t passive = 0U;
    for (std::size_t scalar = 0U;
         scalar < product.fields.scalars.size() && status; ++scalar) {
      const ConstFieldView value =
          as_const(runtime.halo_views[scalar + 1U]);
      if (product.fields.scalar_roles[scalar] ==
          TransportedScalarRole::species) {
        status = ghost_authority(
            value,
            product.equations.thermophysical_predictor().species_reach(),
            runtime.species_ghosts[species++].previous);
      } else {
        status = ghost_authority(
            value,
            product.equations.thermophysical_predictor()
                .passive_scalar_reach(),
            runtime.passive_ghosts[passive++].previous);
      }
    }
    status = product.reductions.consensus(status);
    if (!status) return status;
  }

  // This is the final fallible publication on the restart path.  All
  // thermophysical and ghost-history validation has completed, so a failure
  // cannot leave a visible time/controller authority paired with a partial
  // face-flux lineage.  Legacy images synthesize revision two; exact images
  // retain their distinct accepted/previous revision tokens.
  status = exact_history
               ? runtime.final_flux_writer.initialize_restored_history(
                     product.final_flux, restored_flux,
                     restored_previous_flux)
               : runtime.final_flux_writer.initialize_restored(
                     product.final_flux, restored_flux);
  if (!status) return status;
  runtime.time = std::move(controller);
  runtime.pressure_reference = image.pressure_reference;
  runtime.previous_pressure_reference =
      exact_history ? image.previous_pressure_reference
                    : image.pressure_reference;
  runtime.closed_mass_target =
      exact_history ? image.closed_mass_target : global_mass;
  runtime.pressure_correction_warm_start_valid = false;
  runtime.predictor_diagnostics = {};
  runtime.initialized = true;
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kProductAllocation};
} catch (...) {
  return {StatusCode::invalid_plan, kProductInput};
}

Status ProductDriver::Impl::execute_attempt(
    const StepTime& step, PisoAttemptReport& report,
    PreparedAttemptFinish& prepared) noexcept {
  CompiledCasePlan::Impl& product = *plan.implementation_;
  effective_bdf = step.bdf;
  thermophysical_predictor_calls = 0U;
  temporal_method_fallback = false;
  numerical_failure = {};
  predictor_diagnostics = {};
  momentum_predictor_limiter = {};
  momentum_predictor_solve = {};
  pressure_energy_globalization = {};
  assert(!pending_force_cache);
  assert(!pending_flux.valid());
  clear_pending_attempt_side_state();
  const bool suppress_pressure_correction_warm_start =
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
      g_pressure_correction_warm_start_suppressed.exchange(
          false, std::memory_order_acq_rel);
#else
      false;
#endif
  const bool use_pressure_correction_warm_start =
      pressure_correction_warm_start_valid && step.attempt == 0U &&
      step.origin == StepOrigin::accepted &&
      !suppress_pressure_correction_warm_start;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  g_pressure_correction_warm_start_diagnostic = {
      true, step.origin, step.attempt,
      pressure_correction_warm_start_valid,
      use_pressure_correction_warm_start};
  g_pressure_correction_warm_start_published.store(
      true, std::memory_order_release);
#endif
  const Int3 cells = product.patch.cells;
  const KernelBox full_box{{0, 0, 0}, cells};
  const StateLayers& const_layers = product.layers;
  const BoundaryResolvedValues boundary_values = resolved_boundary_values();
  const auto history = [&](FieldId field, PrimitiveHistory& out) {
    return const_layers.view(StateRole::trial, field, out.trial) &&
                   const_layers.view(StateRole::accepted_n, field,
                                     out.accepted) &&
                   const_layers.view(StateRole::accepted_n_minus_one, field,
                                     out.previous)
               ? Status{}
               : Status{StatusCode::invalid_plan, kProductBinding};
  };
  const auto runtime_write_view = [&](FieldId field, FieldView& out) {
    Status status = product.layers.revise_runtime(
        FieldLifetime::persistent_workspace, field);
    return status ? product.layers.runtime_view(
                        FieldLifetime::persistent_workspace, field, out)
                  : status;
  };
  const auto copy_interior = [&](ConstFieldView source, FieldView target) {
    if (source.interior.x != target.interior.x ||
        source.interior.y != target.interior.y ||
        source.interior.z != target.interior.z ||
        source.components != target.components)
      return Status{StatusCode::invalid_plan, kProductBinding};
    for (std::uint8_t component = 0U; component < source.components;
         ++component)
      for (std::int32_t z = 0; z < cells.z; ++z)
        for (std::int32_t y = 0; y < cells.y; ++y)
          for (std::int32_t x = 0; x < cells.x; ++x)
            target.unchecked({x, y, z}, component) =
                source.unchecked({x, y, z}, component);
    return Status{};
  };
  const auto same_interior_bits = [&](ConstFieldView left,
                                      ConstFieldView right) noexcept {
    if (left.interior.x != right.interior.x ||
        left.interior.y != right.interior.y ||
        left.interior.z != right.interior.z ||
        left.components != right.components)
      return false;
    for (std::uint8_t component = 0U; component < left.components;
         ++component)
      for (std::int32_t z = 0; z < cells.z; ++z)
        for (std::int32_t y = 0; y < cells.y; ++y)
          for (std::int32_t x = 0; x < cells.x; ++x)
            if (product_double_bits(
                    left.unchecked({x, y, z}, component)) !=
                product_double_bits(
                    right.unchecked({x, y, z}, component)))
              return false;
    return true;
  };
  const auto mix_field_numeric = [&](std::uint64_t hash,
                                     ConstFieldView field) noexcept {
    hash = detail::product_mix(hash, field.components);
    for (std::uint8_t component = 0U; component < field.components;
         ++component)
      for (std::int32_t z = 0; z < cells.z; ++z)
        for (std::int32_t y = 0; y < cells.y; ++y)
          for (std::int32_t x = 0; x < cells.x; ++x)
            hash = detail::product_mix(
                hash, product_double_bits(
                          field.unchecked({x, y, z}, component)));
    return hash;
  };
  const auto mix_flux_numeric = [](std::uint64_t hash,
                                   ConstFaceFluxView flux) noexcept {
    const std::array<ConstFaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
    for (ConstFaceFieldView face : faces)
      for (std::int32_t z = 0; z < face.extents.z; ++z)
        for (std::int32_t y = 0; y < face.extents.y; ++y)
          for (std::int32_t x = 0; x < face.extents.x; ++x)
            hash = detail::product_mix(
                hash, product_double_bits(face.unchecked({x, y, z})));
    return hash;
  };
  const auto exchange = [&](HaloEngine& halo, StageId stage,
                            std::size_t count, Status prerequisite = {}) {
    if (prerequisite && (count == 0U || count > halo_views.size()))
      prerequisite = {StatusCode::invalid_plan, kProductCommunication};
    HaloTicket ticket;
    Status status =
        halo.begin(stage, {halo_views.data(), count}, prerequisite, ticket);
    return status ? halo.finish(ticket, {halo_views.data(), count}) : status;
  };
  const auto zero_field = [&](FieldView field) {
    for (std::uint8_t component = 0U; component < field.components;
         ++component)
      for (std::int32_t z = -field.ghosts.z;
           z < field.interior.z + field.ghosts.z; ++z)
        for (std::int32_t y = -field.ghosts.y;
             y < field.interior.y + field.ghosts.y; ++y)
          for (std::int32_t x = -field.ghosts.x;
               x < field.interior.x + field.ghosts.x; ++x)
            field.unchecked({x, y, z}, component) = 0.0;
  };

  reset_stage_timings(10U);
  Status status = transaction.begin(product.layers);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  // Diagnostic atomics are process-local, but every diagnostic branch below
  // may enter operators or reductions.  Freeze an all-rank consistency
  // contract before any such branch can run; tests must arm/clear on every
  // rank.  A partial-rank arm fails closed instead of changing collective
  // call counts.
  const int local_diagnostic_flags[2U]{
      g_candidate_globalization_armed.load(std::memory_order_acquire) ? 1 : 0,
      g_candidate_globalization_published.load(std::memory_order_acquire)
          ? 1
          : 0};
  int minimum_diagnostic_flags[2U]{};
  int maximum_diagnostic_flags[2U]{};
  const int minimum_flags_rc =
      MPI_Allreduce(local_diagnostic_flags, minimum_diagnostic_flags, 2,
                    MPI_INT, MPI_MIN, communicator);
  const int maximum_flags_rc =
      MPI_Allreduce(local_diagnostic_flags, maximum_diagnostic_flags, 2,
                    MPI_INT, MPI_MAX, communicator);
  if (status &&
      (minimum_flags_rc != MPI_SUCCESS || maximum_flags_rc != MPI_SUCCESS))
    status = {StatusCode::mpi_failure, kProductCollective};
  if (status &&
      (minimum_diagnostic_flags[0U] != maximum_diagnostic_flags[0U] ||
       minimum_diagnostic_flags[1U] != maximum_diagnostic_flags[1U]))
    status = {StatusCode::invalid_plan, kProductCollective};
#endif
  status = product.reductions.consensus(status);
  attempt_stage = 10U;
  // Every state-layer field is attempt-owned, even when its final value is
  // produced late in the graph.  Revise all trial authorities before the
  // first write so commit completeness is mechanical and rollback restores
  // the whole coupled state.
  const std::array<FieldId, 6U> primitive_state_fields{
      product.fields.rho,
      product.fields.velocity,
      product.fields.pressure,
      product.fields.enthalpy,
      product.fields.temperature,
      product.fields.enthalpy_nonadvective_rate};
  for (FieldId field : primitive_state_fields) {
    if (status) status = transaction.revise_trial(field);
  }
  for (FieldId field : product.fields.scalars) {
    if (status) status = transaction.revise_trial(field);
  }
  for (FieldId field : product.fields.scalar_nonadvective_rates) {
    if (status) status = transaction.revise_trial(field);
  }
  PrimitiveHistory rho_history;
  PrimitiveHistory velocity_history;
  PrimitiveHistory pressure_history;
  PrimitiveHistory enthalpy_history;
  PrimitiveHistory temperature_history;
  if (status) status = history(product.fields.rho, rho_history);
  if (status) status = history(product.fields.velocity, velocity_history);
  if (status) status = history(product.fields.pressure, pressure_history);
  if (status) status = history(product.fields.enthalpy, enthalpy_history);
  if (status) status = history(product.fields.temperature, temperature_history);
  FieldView trial_velocity;
  FieldView trial_pressure;
  FieldView trial_density;
  if (status)
    status = product.layers.view(StateRole::trial, product.fields.rho,
                                 trial_density);
  if (status)
    status = product.layers.view(StateRole::trial, product.fields.velocity,
                                 trial_velocity);
  if (status)
    status = product.layers.view(StateRole::trial, product.fields.pressure,
                                 trial_pressure);
  if (status) status = copy_interior(velocity_history.accepted, trial_velocity);
  if (status) status = copy_interior(pressure_history.accepted, trial_pressure);

  std::size_t species_index = 0U;
  std::size_t passive_index = 0U;
  for (std::size_t index = 0U;
       index < product.fields.scalars.size() && status; ++index) {
    PrimitiveHistory scalar_history;
    status = history(product.fields.scalars[index], scalar_history);
    PrimitiveHistory rate_history;
    if (status)
      status = history(product.fields.scalar_nonadvective_rates[index],
                       rate_history);
    FieldView scalar_trial;
    if (status)
      status = product.layers.view(StateRole::trial,
                                   product.fields.scalars[index],
                                   scalar_trial);
    if (!status) break;
    if (product.fields.scalar_roles[index] ==
        TransportedScalarRole::species) {
      species_history[species_index] = scalar_history;
      species_accepted[species_index] = scalar_history.accepted;
      species_previous[species_index] =
          step.bdf.order == 2U ? scalar_history.previous : ConstFieldView{};
      species_trial[species_index] = scalar_trial;
      species_rate_history[species_index] = {
          rate_history.accepted,
          step.bdf.order == 2U ? rate_history.previous : ConstFieldView{}};
      ++species_index;
    } else {
      passive_history[passive_index] = scalar_history;
      passive_accepted[passive_index] = scalar_history.accepted;
      passive_previous[passive_index] =
          step.bdf.order == 2U ? scalar_history.previous : ConstFieldView{};
      passive_trial[passive_index] = scalar_trial;
      passive_rate_history[passive_index] = {
          rate_history.accepted,
          step.bdf.order == 2U ? rate_history.previous : ConstFieldView{}};
      ++passive_index;
    }
  }
  PrimitiveHistory enthalpy_rate_history;
  if (status)
    status = history(product.fields.enthalpy_nonadvective_rate,
                     enthalpy_rate_history);

  // Predictor halo: accepted h and all accepted transported scalars.
  std::size_t halo_count = 0U;
  FieldView accepted_enthalpy_mutable;
  if (status)
    status = product.layers.view(StateRole::accepted_n,
                                 product.fields.enthalpy,
                                 accepted_enthalpy_mutable);
  if (status) halo_views[halo_count++] = accepted_enthalpy_mutable;
  for (std::size_t index = 0U;
       index < product.fields.scalars.size() && status; ++index) {
    FieldView value;
    status = product.layers.view(StateRole::accepted_n,
                                 product.fields.scalars[index], value);
    if (status) halo_views[halo_count++] = value;
  }
  if (status) status = exchange(product.stage_halos[0U], 10U, halo_count);
  if (status)
    status = resolve_static_boundary_values(
        communicator, product.boundary, product.boundary_specs,
        product.geometry, product.patch,
        product.thermodynamics, pressure_reference, rho_history.accepted,
        velocity_history.accepted, enthalpy_history.accepted,
        {species_accepted.data(), species_accepted.size()},
        {passive_accepted.data(), passive_accepted.size()},
        {species_values.data(), species_values.size()}, boundary_scalars,
        boundary_vectors, boundary_normal_gradients, false);
  if (status) {
    accepted_enthalpy_mutable = halo_views[0U];
    status = apply_boundary_ghosts(BoundaryStage::enthalpy, product.boundary,
                                   {&accepted_enthalpy_mutable, 1U},
                                   boundary_values);
  }
  if (status && product.fields.scalars.size() != 0U) {
    status = apply_boundary_ghosts(
        BoundaryStage::scalar, product.boundary,
        {halo_views.data() + 1U, product.fields.scalars.size()},
        boundary_values);
  }
  if (status) enthalpy_history.accepted = as_const(accepted_enthalpy_mutable);
  species_index = 0U;
  passive_index = 0U;
  for (std::size_t index = 0U;
       index < product.fields.scalars.size() && status; ++index) {
    const ConstFieldView accepted = as_const(halo_views[index + 1U]);
    if (product.fields.scalar_roles[index] ==
        TransportedScalarRole::species)
      species_accepted[species_index++] = accepted;
    else
      passive_accepted[passive_index++] = accepted;
  }
  const auto publish_predictor_ghost_authority =
      [&](HaloEngine& halo, ConstFieldView view, std::uint8_t reach,
          ThermophysicalGhostAuthority& authority) noexcept {
        if (halo.ghost_revision(view.field) != view.revision ||
            halo.instance_identity() == 0U || reach == 0U ||
            view.ghosts.x < reach || view.ghosts.y < reach ||
            view.ghosts.z < reach) {
          return Status{StatusCode::invalid_plan, kProductCommunication};
        }
        authority = {halo.instance_identity(),
                     view.field,
                     view.revision,
                     view.storage_identity,
                     view.revision_domain,
                     product.geometry.topology_revision(),
                     product.boundary.revision(),
                     reach};
        return Status{};
      };
  if (status)
    status = publish_predictor_ghost_authority(
        product.stage_halos[0U], enthalpy_history.accepted,
        product.equations.thermophysical_predictor().enthalpy_reach(),
        enthalpy_ghosts.accepted);
  for (std::size_t index = 0U; index < species_accepted.size() && status;
       ++index) {
    status = publish_predictor_ghost_authority(
        product.stage_halos[0U], species_accepted[index],
        product.equations.thermophysical_predictor().species_reach(),
        species_ghosts[index].accepted);
  }
  for (std::size_t index = 0U; index < passive_accepted.size() && status;
       ++index) {
    status = publish_predictor_ghost_authority(
        product.stage_halos[0U], passive_accepted[index],
        product.equations.thermophysical_predictor().passive_scalar_reach(),
        passive_ghosts[index].accepted);
  }

  ConstFaceFluxView accepted_flux;
  ConstFaceFluxView previous_flux;
  if (status)
    status = final_flux_writer.committed(product.final_flux, accepted_flux);
  if (status && step.bdf.order == 2U)
    status = final_flux_writer.committed_previous(product.final_flux,
                                                  previous_flux);
  FaceFluxView provisional_flux;
  if (status)
    status = product.phi_workspace.workspace_view(
        1U, transaction.attempt_identity(), provisional_flux);
  FieldView predictor_n;
  FieldView predictor_nm1;
  FieldView predictor_low_bundle;
  if (status)
    status = runtime_write_view(product.fields.predictor_accepted_advection,
                                predictor_n);
  if (status)
    status = runtime_write_view(product.fields.predictor_previous_advection,
                                predictor_nm1);
  if (status)
    status = runtime_write_view(product.fields.predictor_low_bundle,
                                predictor_low_bundle);
  FieldView trial_enthalpy;
  if (status)
    status = product.layers.view(StateRole::trial, product.fields.enthalpy,
                                 trial_enthalpy);
  const auto bundle_component = [](FieldView bundle, std::size_t component,
                                   FieldId semantic_field) noexcept {
    FieldView view = bundle;
    view.base += component * bundle.component_stride;
    view.components = 1U;
    view.field = semantic_field;
    return view;
  };
  FieldView predictor_low_density;
  FieldView predictor_low_enthalpy;
  if (status && predictor_low_bundle.components !=
                    2U + product.fields.scalars.size()) {
    status = {StatusCode::invalid_plan, kProductBinding};
  }
  if (status) {
    predictor_low_density =
        bundle_component(predictor_low_bundle, 0U, product.fields.rho);
    predictor_low_enthalpy =
        bundle_component(predictor_low_bundle, 1U, product.fields.enthalpy);
    species_index = 0U;
    passive_index = 0U;
    for (std::size_t index = 0U; index < product.fields.scalars.size();
         ++index) {
      const FieldView component = bundle_component(
          predictor_low_bundle, index + 2U, product.fields.scalars[index]);
      if (product.fields.scalar_roles[index] ==
          TransportedScalarRole::species) {
        species_low[species_index++] = component;
      } else {
        passive_low[passive_index++] = component;
      }
    }
  }
  ThermophysicalPredictorInput predictor_input;
  predictor_input.dt = step.dt;
  predictor_input.bdf = effective_bdf;
  predictor_input.time = step.generation;
  predictor_input.geometry = product.geometry.topology_revision();
  predictor_input.boundary = product.boundary.revision();
  predictor_input.transport = product.transport.fingerprint();
  predictor_input.density_accepted = rho_history.accepted;
  predictor_input.density_previous =
      effective_bdf.order == 2U ? rho_history.previous : ConstFieldView{};
  predictor_input.enthalpy_accepted = enthalpy_history.accepted;
  predictor_input.enthalpy_previous =
      effective_bdf.order == 2U ? enthalpy_history.previous
                                : ConstFieldView{};
  predictor_input.species_accepted = {species_accepted.data(),
                                      species_accepted.size()};
  predictor_input.species_previous = {species_previous.data(),
                                      species_previous.size()};
  predictor_input.passive_scalars_accepted = {passive_accepted.data(),
                                              passive_accepted.size()};
  predictor_input.passive_scalars_previous = {passive_previous.data(),
                                              passive_previous.size()};
  predictor_input.mass_flux_accepted = accepted_flux;
  predictor_input.mass_flux_previous = previous_flux;
  predictor_input.enthalpy_nonadvective_rhs = {
      enthalpy_rate_history.accepted,
      effective_bdf.order == 2U ? enthalpy_rate_history.previous
                                : ConstFieldView{}};
  predictor_input.species_nonadvective_rhs = {
      species_rate_history.data(), species_rate_history.size()};
  predictor_input.passive_scalar_nonadvective_rhs = {
      passive_rate_history.data(), passive_rate_history.size()};
  predictor_input.enthalpy_ghosts = enthalpy_ghosts;
  predictor_input.species_ghosts = {species_ghosts.data(),
                                    species_ghosts.size()};
  predictor_input.passive_scalar_ghosts = {passive_ghosts.data(),
                                           passive_ghosts.size()};
  ThermophysicalPredictorOutput predictor_output{
      trial_enthalpy,
      {species_trial.data(), species_trial.size()},
      trial_density,
      predictor_n,
      predictor_nm1,
      {passive_trial.data(), passive_trial.size()},
      predictor_low_density,
      predictor_low_enthalpy,
      {species_low.data(), species_low.size()},
      {passive_low.data(), passive_low.size()},
      provisional_flux};
  ThermophysicalPredictorCertificate predictor_certificate;
  halo_count = 0U;
  if (status) halo_views[halo_count++] = predictor_low_enthalpy;
  for (FieldView value : species_low) halo_views[halo_count++] = value;
  for (FieldView value : passive_low) halo_views[halo_count++] = value;
  const ThermophysicalPredictorSlowPath predictor_slow_path{
      &product.predictor_donor_halo, 16U,
      {halo_views.data(), halo_count}, boundary_values,
      &enthalpy_endpoint, &resources,
      product.ibm_equations.has_value() ? &*product.ibm_equations : nullptr};
  thermophysical_predictor_calls = 1U;
  status = product.equations.thermophysical_predictor().predict(
      communicator, status, predictor_input, predictor_output,
      predictor_slow_path, predictor_diagnostics, predictor_certificate);
  const bool fallback_eligible =
      !status && status.code == StatusCode::numerical_failure &&
      step.bdf.order == 2U && predictor_diagnostics.failure.valid &&
      predictor_diagnostics.failure.reason ==
          ThermophysicalPredictorFailureReason::
              low_bdf_source_base_admissibility;
  if (fallback_eligible) {
    const ThermophysicalPredictorDiagnostics requested_diagnostics =
        predictor_diagnostics;
    effective_bdf = {1.0 / step.dt, -1.0 / step.dt, 0.0, 1U};
    temporal_method_fallback = true;
    thermophysical_predictor_calls = 2U;
    predictor_input.bdf = effective_bdf;
    predictor_input.density_previous = {};
    predictor_input.enthalpy_previous = {};
    predictor_input.mass_flux_previous = {};
    predictor_input.enthalpy_nonadvective_rhs.previous = {};
    predictor_input.enthalpy_ghosts.previous = {};
    for (std::size_t index = 0U; index < species_previous.size(); ++index) {
      species_previous[index] = {};
      species_rate_history[index].previous = {};
    }
    for (std::size_t index = 0U; index < passive_previous.size(); ++index) {
      passive_previous[index] = {};
      passive_rate_history[index].previous = {};
    }
    predictor_diagnostics = {};
    predictor_certificate = {};
    status = product.equations.thermophysical_predictor().predict(
        communicator, Status{}, predictor_input, predictor_output,
        predictor_slow_path, predictor_diagnostics, predictor_certificate);
    const std::uint32_t requested_collectives =
        requested_diagnostics.blocking_collectives;
    const std::uint32_t effective_collectives =
        predictor_diagnostics.blocking_collectives;
    predictor_diagnostics.blocking_collectives =
        effective_collectives >
                std::numeric_limits<std::uint32_t>::max() -
                    requested_collectives
            ? std::numeric_limits<std::uint32_t>::max()
            : requested_collectives + effective_collectives;
  }

  // Publish predicted h/scalar ghosts before EOS/transport evaluation.
  halo_count = 0U;
  if (status) halo_views[halo_count++] = trial_enthalpy;
  for (FieldView value : species_trial) halo_views[halo_count++] = value;
  for (FieldView value : passive_trial) halo_views[halo_count++] = value;
  if (status) status = exchange(product.stage_halos[1U], 15U, halo_count);
  if (status) {
    trial_enthalpy = halo_views[0U];
    status = apply_boundary_ghosts(BoundaryStage::enthalpy, product.boundary,
                                   {&trial_enthalpy, 1U}, boundary_values);
  }
  if (status && product.fields.scalars.size() != 0U)
    status = apply_boundary_ghosts(
        BoundaryStage::scalar, product.boundary,
        {halo_views.data() + 1U, product.fields.scalars.size()},
        boundary_values);
  species_index = 0U;
  passive_index = 0U;
  for (std::size_t index = 0U;
       index < product.fields.scalars.size() && status; ++index) {
    if (product.fields.scalar_roles[index] ==
        TransportedScalarRole::species)
      species_trial[species_index++] = halo_views[index + 1U];
    else
      passive_trial[passive_index++] = halo_views[index + 1U];
  }
  if (status)
    status = publish_predictor_ghost_authority(
        product.stage_halos[1U], as_const(trial_enthalpy),
        product.equations.thermophysical_predictor().enthalpy_reach(),
        trial_enthalpy_ghost);
  for (std::size_t index = 0U; index < species_trial.size() && status;
       ++index) {
    status = publish_predictor_ghost_authority(
        product.stage_halos[1U], as_const(species_trial[index]),
        product.equations.thermophysical_predictor().species_reach(),
        trial_species_ghosts[index]);
  }
  for (std::size_t index = 0U; index < passive_trial.size() && status;
       ++index) {
    status = publish_predictor_ghost_authority(
        product.stage_halos[1U], as_const(passive_trial[index]),
        product.equations.thermophysical_predictor().passive_scalar_reach(),
        trial_passive_ghosts[index]);
  }

  FieldView trial_temperature;
  FieldView molecular_viscosity;
  FieldView effective_viscosity;
  FieldView compressibility;
  FieldView enthalpy_compressibility;
  FieldView conductivity;
  FieldView heat_capacity;
  FieldView enthalpy_diffusivity;
  FieldView eos_density;
  if (status)
    status = product.layers.view(StateRole::trial, product.fields.temperature,
                                 trial_temperature);
  if (status)
    status = runtime_write_view(product.fields.molecular_viscosity,
                                molecular_viscosity);
  if (status)
    status = runtime_write_view(product.fields.effective_viscosity,
                                effective_viscosity);
  if (status)
    status = runtime_write_view(product.fields.compressibility,
                                compressibility);
  if (status)
    status = runtime_write_view(product.fields.enthalpy_compressibility,
                                enthalpy_compressibility);
  if (status)
    status = runtime_write_view(product.fields.thermal_conductivity,
                                conductivity);
  if (status)
    status = runtime_write_view(product.fields.heat_capacity,
                                heat_capacity);
  if (status)
    status = runtime_write_view(product.fields.enthalpy_diffusivity,
                                enthalpy_diffusivity);
  if (status)
    status = runtime_write_view(product.fields.eos_density, eos_density);

  const auto refresh_live_effective_thermal_ghosts =
      [&](StageId stage, Status prerequisite) {
        if (product.transport.kernel() != TransportKernel::coast_native_air)
          return prerequisite;
        if (prerequisite)
          prerequisite = runtime_write_view(
              product.fields.thermal_conductivity, conductivity);
        if (prerequisite)
          prerequisite = runtime_write_view(
              product.fields.enthalpy_diffusivity, enthalpy_diffusivity);
        if (prerequisite)
          prerequisite = refresh_coast_native_air_effective_thermal_transport(
              product.transport, as_const(molecular_viscosity),
              as_const(effective_viscosity), as_const(heat_capacity),
              conductivity, enthalpy_diffusivity);
        prerequisite = exchange_effective_thermal_ghosts(
            product.coupled_thermal_halo, stage, product.boundary, conductivity,
            enthalpy_diffusivity, prerequisite);
        return product.reductions.consensus(prerequisite);
      };

  const auto capture_thermo_failure = [&](Int3 cell, double p_ref,
                                          Status failure) noexcept {
    NumericalFailureContext context;
    context.valid = true;
    context.field = NumericalFailureField::enthalpy;
    context.first_bad_contributor =
        NumericalFailureContributor::thermodynamic_inversion;
    context.field_id = product.fields.enthalpy;
    context.failure = failure;
    context.stage = 15U;
    context.attempted_step =
        step.accepted_step == std::numeric_limits<std::uint64_t>::max()
            ? step.accepted_step
            : step.accepted_step + 1U;
    context.generation = step.generation;
    context.time = step.time;
    context.target_time = step.time + step.dt;
    context.dt_before = step.dt;
    context.method_before = effective_bdf;
    context.origin_before = step.origin;
    context.dt_after = std::numeric_limits<double>::quiet_NaN();
    context.failed_value = trial_enthalpy.unchecked(cell, 0U);
    context.temperature_before =
        temperature_history.accepted.unchecked(cell, 0U);
    context.density_before = rho_history.accepted.unchecked(cell, 0U);
    context.density_previous =
        effective_bdf.order == 2U
            ? rho_history.previous.unchecked(cell, 0U)
            : std::numeric_limits<double>::quiet_NaN();
    context.density_predicted = trial_density.unchecked(cell, 0U);
    context.enthalpy_before =
        enthalpy_history.accepted.unchecked(cell, 0U);
    context.enthalpy_previous =
        effective_bdf.order == 2U
            ? enthalpy_history.previous.unchecked(cell, 0U)
            : std::numeric_limits<double>::quiet_NaN();
    context.cp_before = heat_capacity.unchecked(cell, 0U);
    context.pressure_absolute =
        p_ref + trial_pressure.unchecked(cell, 0U);
    context.enthalpy_accepted_revision = enthalpy_history.accepted.revision;
    context.enthalpy_previous_revision =
        effective_bdf.order == 2U ? enthalpy_history.previous.revision : 0U;
    context.nonadvective_accepted_revision =
        enthalpy_rate_history.accepted.revision;
    context.nonadvective_previous_revision =
        effective_bdf.order == 2U
            ? enthalpy_rate_history.previous.revision
            : 0U;
    context.mass_flux_accepted_revision = accepted_flux.revision;
    context.mass_flux_previous_revision =
        effective_bdf.order == 2U ? previous_flux.revision : 0U;
    context.temperature_accepted_revision =
        temperature_history.accepted.revision;
    context.conductivity_revision = conductivity.revision;
    context.effective_viscosity_revision = effective_viscosity.revision;

    int rank = -1;
    if (MPI_Comm_rank(communicator, &rank) == MPI_SUCCESS) {
      context.rank = rank;
    }
    context.global_index = {product.patch.begin.x + cell.x,
                            product.patch.begin.y + cell.y,
                            product.patch.begin.z + cell.z};
    const Int3 global_cells = product.geometry.global_cells();
    context.global_cell =
        diagnostic_global_cell(context.global_index, global_cells);
    const auto near_lower = [](std::int32_t index) { return index < 2; };
    const auto near_upper = [](std::int32_t index, std::int32_t extent) {
      return index >= extent - 2;
    };
    context.physical_boundary_stencil =
        near_lower(context.global_index.x) ||
        near_upper(context.global_index.x, global_cells.x) ||
        near_lower(context.global_index.y) ||
        near_upper(context.global_index.y, global_cells.y) ||
        near_lower(context.global_index.z) ||
        near_upper(context.global_index.z, global_cells.z);
    context.mpi_boundary_stencil =
        ((near_lower(cell.x) && product.patch.begin.x != 0) ||
         (near_upper(cell.x, cells.x) &&
          product.patch.begin.x + cells.x != global_cells.x) ||
         (near_lower(cell.y) && product.patch.begin.y != 0) ||
         (near_upper(cell.y, cells.y) &&
          product.patch.begin.y + cells.y != global_cells.y) ||
         (near_lower(cell.z) && product.patch.begin.z != 0) ||
         (near_upper(cell.z, cells.z) &&
          product.patch.begin.z + cells.z != global_cells.z));

    double independent_sum = 0.0;
    const std::size_t full_species_count = species_trial.size() + 1U;
    context.mass_fraction_count = std::min(
        full_species_count, kNumericalFailureMassFractionCapacity);
    context.mass_fractions_truncated =
        full_species_count > kNumericalFailureMassFractionCapacity;
    for (std::size_t species = 0U; species < species_trial.size(); ++species) {
      const double value = species_trial[species].unchecked(cell, 0U);
      species_values[species] = value;
      independent_sum += value;
      if (species < context.mass_fraction_count)
        context.mass_fractions[species] = value;
    }
    if (species_trial.size() < context.mass_fraction_count)
      context.mass_fractions[species_trial.size()] = 1.0 - independent_sum;

    double lower_cp = 0.0;
    double lower_gas = 0.0;
    double upper_cp = 0.0;
    double upper_gas = 0.0;
    const Status lower = product.thermodynamics.mixture_enthalpy(
        product.thermodynamics.minimum_temperature(),
        {species_values.data(), species_values.size()},
        context.allowed_minimum, lower_cp, lower_gas);
    const Status upper = product.thermodynamics.mixture_enthalpy(
        product.thermodynamics.maximum_temperature(),
        {species_values.data(), species_values.size()},
        context.allowed_maximum, upper_cp, upper_gas);
    if (lower && upper) {
      if (context.failed_value < context.allowed_minimum && lower_cp > 0.0) {
        context.temperature_estimate =
            product.thermodynamics.minimum_temperature() +
            (context.failed_value - context.allowed_minimum) / lower_cp;
      } else if (context.failed_value > context.allowed_maximum &&
                 upper_cp > 0.0) {
        context.temperature_estimate =
            product.thermodynamics.maximum_temperature() +
            (context.failed_value - context.allowed_maximum) / upper_cp;
      } else {
        context.temperature_estimate = context.temperature_before;
      }
    } else {
      context.allowed_minimum = std::numeric_limits<double>::quiet_NaN();
      context.allowed_maximum = std::numeric_limits<double>::quiet_NaN();
      context.temperature_estimate =
          std::numeric_limits<double>::quiet_NaN();
    }

    const ConvectionScheme convection =
        product.equations.thermophysical_predictor().enthalpy_convection();
    ConvectionPointDiagnostic accepted_diagnostic;
    Status diagnostic = diagnose_cartesian_convection_point(
        product.equations.kernels(), convection, accepted_flux,
        enthalpy_history.accepted, 0U, cell, accepted_diagnostic);
    if (diagnostic) {
      context.mass_divergence_accepted =
          accepted_diagnostic.mass_divergence;
      context.advection_accepted = accepted_diagnostic.divergence;
      context.face_envelope_checked =
          accepted_diagnostic.face_envelope_checked;
      context.face_envelope_valid = accepted_diagnostic.face_envelope_valid;
      context.maximum_face_envelope_violation =
          accepted_diagnostic.maximum_face_envelope_violation;
      context.selected_face_value = accepted_diagnostic.selected_face_value;
      context.selected_donor_minimum =
          accepted_diagnostic.selected_donor_minimum;
      context.selected_donor_maximum =
          accepted_diagnostic.selected_donor_maximum;
    } else {
      context.mass_divergence_accepted =
          std::numeric_limits<double>::quiet_NaN();
      context.advection_accepted =
          std::numeric_limits<double>::quiet_NaN();
    }
    if (effective_bdf.order == 2U) {
      ConvectionPointDiagnostic previous_diagnostic;
      diagnostic = diagnose_cartesian_convection_point(
          product.equations.kernels(), convection, previous_flux,
          enthalpy_history.previous, 0U, cell, previous_diagnostic);
      if (diagnostic) {
        context.mass_divergence_previous =
            previous_diagnostic.mass_divergence;
        context.advection_previous = previous_diagnostic.divergence;
        context.face_envelope_checked =
            context.face_envelope_checked &&
            previous_diagnostic.face_envelope_checked;
        context.face_envelope_valid =
            context.face_envelope_valid &&
            previous_diagnostic.face_envelope_valid;
        if (previous_diagnostic.maximum_face_envelope_violation >
            context.maximum_face_envelope_violation) {
          context.maximum_face_envelope_violation =
              previous_diagnostic.maximum_face_envelope_violation;
          context.selected_face_value =
              previous_diagnostic.selected_face_value;
          context.selected_donor_minimum =
              previous_diagnostic.selected_donor_minimum;
          context.selected_donor_maximum =
              previous_diagnostic.selected_donor_maximum;
        }
      } else {
        context.mass_divergence_previous =
            std::numeric_limits<double>::quiet_NaN();
        context.advection_previous =
            std::numeric_limits<double>::quiet_NaN();
      }
    }
    context.nonadvective_accepted =
        enthalpy_rate_history.accepted.base == nullptr
            ? 0.0
            : enthalpy_rate_history.accepted.unchecked(cell, 0U);
    context.nonadvective_previous =
        effective_bdf.order != 2U ||
                enthalpy_rate_history.previous.base == nullptr
            ? 0.0
            : enthalpy_rate_history.previous.unchecked(cell, 0U);

    const double ratio = effective_bdf.order == 2U
                             ? -effective_bdf.a1 * step.dt - 1.0
                             : 0.0;
    const double extrapolate_accepted =
        effective_bdf.order == 2U ? 1.0 + ratio : 1.0;
    const double extrapolate_previous =
        effective_bdf.order == 2U ? -ratio : 0.0;
    const double denominator = effective_bdf.a0 * context.density_predicted;
    if (std::isfinite(denominator) && denominator > 0.0) {
      const double conservative_history =
          -effective_bdf.a1 * context.density_before *
          context.enthalpy_before;
      const double conservative_previous =
          effective_bdf.order == 2U
              ? -effective_bdf.a2 * context.density_previous *
                    context.enthalpy_previous
              : 0.0;
      context.conservative_history_value =
          (conservative_history + conservative_previous) / denominator;
      context.accepted_advection_delta =
          -extrapolate_accepted * context.advection_accepted / denominator;
      context.accepted_nonadvective_delta =
          extrapolate_accepted * context.nonadvective_accepted / denominator;
      context.previous_advection_delta =
          -extrapolate_previous * context.advection_previous / denominator;
      context.previous_nonadvective_delta =
          extrapolate_previous * context.nonadvective_previous / denominator;
      context.reconstructed_value =
          context.conservative_history_value +
          context.accepted_advection_delta +
          context.accepted_nonadvective_delta +
          context.previous_advection_delta +
          context.previous_nonadvective_delta;
      double cumulative = context.conservative_history_value;
      if (!diagnostic_inside(cumulative, context.allowed_minimum,
                             context.allowed_maximum)) {
        context.first_bad_contributor =
            NumericalFailureContributor::conservative_history;
      }
      const auto add_contributor = [&](double value,
                                       NumericalFailureContributor which) {
        if (context.first_bad_contributor ==
            NumericalFailureContributor::thermodynamic_inversion) {
          const double next = cumulative + value;
          if (diagnostic_inside(cumulative, context.allowed_minimum,
                                context.allowed_maximum) &&
              !diagnostic_inside(next, context.allowed_minimum,
                                 context.allowed_maximum)) {
            context.first_bad_contributor = which;
          }
          cumulative = next;
        } else {
          cumulative += value;
        }
      };
      add_contributor(context.accepted_advection_delta,
                      NumericalFailureContributor::accepted_advection);
      add_contributor(context.accepted_nonadvective_delta,
                      NumericalFailureContributor::accepted_nonadvective);
      add_contributor(context.previous_advection_delta,
                      NumericalFailureContributor::previous_advection);
      add_contributor(context.previous_nonadvective_delta,
                      NumericalFailureContributor::previous_nonadvective);
    }

    context.counterfactual_be_density =
        context.density_before -
        step.dt * context.mass_divergence_accepted;
    const double be_conservative =
        context.density_before * context.enthalpy_before -
        step.dt *
            (context.advection_accepted - context.nonadvective_accepted);
    context.counterfactual_be_enthalpy =
        std::isfinite(context.counterfactual_be_density) &&
                context.counterfactual_be_density > 0.0
            ? be_conservative / context.counterfactual_be_density
            : std::numeric_limits<double>::quiet_NaN();
    context.counterfactual_be_admissible =
        context.counterfactual_be_density > 0.0 &&
        diagnostic_inside(context.counterfactual_be_enthalpy,
                          context.allowed_minimum,
                          context.allowed_maximum);

    ConstFieldView accepted_gradient;
    if (const_layers.runtime_view(FieldLifetime::persistent_workspace,
                                  product.fields.velocity_gradient,
                                  accepted_gradient)) {
      context.velocity_gradient_revision = accepted_gradient.revision;
      VelocityGradient gradient;
      for (std::uint8_t component = 0U; component < 9U; ++component)
        gradient.value[component] =
            accepted_gradient.unchecked(cell, component);
      const Status dissipated = newtonian_viscous_dissipation(
          gradient, effective_viscosity.unchecked(cell, 0U),
          context.viscous_dissipation_accepted);
      if (!dissipated)
        context.viscous_dissipation_accepted =
            std::numeric_limits<double>::quiet_NaN();
    } else {
      context.viscous_dissipation_accepted =
          std::numeric_limits<double>::quiet_NaN();
    }

    const auto accepted_thermal_point = [&](Int3 point, double& temperature,
                                            double& lambda) noexcept {
      for (std::size_t species = 0U; species < species_accepted.size();
           ++species)
        species_values[species] =
            species_accepted[species].unchecked(point, 0U);
      ThermoState thermo;
      Status evaluated = product.thermodynamics.evaluate(
          1.0, enthalpy_history.accepted.unchecked(point, 0U),
          {species_values.data(), species_values.size()}, {}, thermo);
      MolecularTransportState transport;
      if (evaluated)
        evaluated = product.transport.evaluate(
            thermo.temperature,
            {species_values.data(), species_values.size()}, transport);
      if (!evaluated) return evaluated;
      temperature = thermo.temperature;
      lambda = transport.conductivity;
      return Status{};
    };
    const std::array<Int3, 7U> points{{
        cell,
        {cell.x - 1, cell.y, cell.z},
        {cell.x + 1, cell.y, cell.z},
        {cell.x, cell.y - 1, cell.z},
        {cell.x, cell.y + 1, cell.z},
        {cell.x, cell.y, cell.z - 1},
        {cell.x, cell.y, cell.z + 1},
    }};
    std::array<double, 7U> point_temperature{};
    std::array<double, 7U> point_conductivity{};
    Status conduction_status;
    for (std::size_t point = 0U; point < points.size() && conduction_status;
         ++point) {
      conduction_status = accepted_thermal_point(
          points[point], point_temperature[point], point_conductivity[point]);
    }
    const auto transmissibility = [&](CartesianAxis axis, Int3 face,
                                      double left, double right) {
      const std::int32_t normal =
          axis == CartesianAxis::x
              ? face.x
              : (axis == CartesianAxis::y ? face.y : face.z);
      const double face_coordinate = detail::face_coordinate(
          product.equations.kernels(), axis, normal);
      const double left_coordinate = detail::centre_coordinate(
          product.equations.kernels(), axis, normal - 1);
      const double right_coordinate = detail::centre_coordinate(
          product.equations.kernels(), axis, normal);
      return detail::face_area(product.equations.kernels(), axis, face) /
             ((face_coordinate - left_coordinate) / left +
              (right_coordinate - face_coordinate) / right);
    };
    if (conduction_status) {
      const double txm = transmissibility(
          CartesianAxis::x, cell, point_conductivity[1U],
          point_conductivity[0U]);
      const double txp = transmissibility(
          CartesianAxis::x, {cell.x + 1, cell.y, cell.z},
          point_conductivity[0U], point_conductivity[2U]);
      const double tym = transmissibility(
          CartesianAxis::y, cell, point_conductivity[3U],
          point_conductivity[0U]);
      const double typ = transmissibility(
          CartesianAxis::y, {cell.x, cell.y + 1, cell.z},
          point_conductivity[0U], point_conductivity[4U]);
      const double tzm = transmissibility(
          CartesianAxis::z, cell, point_conductivity[5U],
          point_conductivity[0U]);
      const double tzp = transmissibility(
          CartesianAxis::z, {cell.x, cell.y, cell.z + 1},
          point_conductivity[0U], point_conductivity[6U]);
      context.diffusion_accepted =
          (txp * (point_temperature[2U] - point_temperature[0U]) -
           txm * (point_temperature[0U] - point_temperature[1U]) +
           typ * (point_temperature[4U] - point_temperature[0U]) -
           tym * (point_temperature[0U] - point_temperature[3U]) +
           tzp * (point_temperature[6U] - point_temperature[0U]) -
           tzm * (point_temperature[0U] - point_temperature[5U])) /
          detail::cell_volume(product.equations.kernels(), cell);
    } else {
      context.diffusion_accepted =
          std::numeric_limits<double>::quiet_NaN();
    }

    context.immersed_interface_cell = false;
    if (product.topology.has_value()) {
      const Span<const ImmersedLink> links = product.topology->links();
      for (std::size_t link = 0U; link < links.size; ++link) {
        const Int3 fluid = links.data[link].fluid_local_index;
        if (fluid.x == cell.x && fluid.y == cell.y && fluid.z == cell.z) {
          context.immersed_interface_cell = true;
          break;
        }
      }
    }
    context.explicit_source_accepted = 0.0;
    context.implicit_sink_accepted = 0.0;
    const ContributionPlan* contribution_plan = product.contributions.plan();
    const bool no_contributions = contribution_plan != nullptr &&
                                  contribution_plan->contributions().size == 0U;
    context.rate_breakdown_complete =
        conduction_status && no_contributions &&
        !context.immersed_interface_cell &&
        std::isfinite(context.diffusion_accepted) &&
        std::isfinite(context.viscous_dissipation_accepted);
    if (context.rate_breakdown_complete) {
      context.pressure_work_accepted =
          context.nonadvective_accepted - context.diffusion_accepted -
          context.viscous_dissipation_accepted;
    } else {
      context.pressure_work_accepted =
          std::numeric_limits<double>::quiet_NaN();
    }
    numerical_failure = context;
  };

  const PressureReferenceKind pressure_reference_kind =
      product.equations.pressure_reference().kind();
  const auto update_thermo_from_authority =
      [&](double p_ref, bool coupled_pressure_authority) {
    for (std::int32_t z = 0; z < cells.z; ++z) {
      for (std::int32_t y = 0; y < cells.y; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          for (std::size_t species = 0U; species < species_trial.size();
               ++species)
            species_values[species] =
                species_trial[species].unchecked(cell, 0U);
          const Real3 velocity{
              trial_velocity.unchecked(cell, 0U),
              trial_velocity.unchecked(cell, 1U),
              trial_velocity.unchecked(cell, 2U)};
          double density_value = trial_density.unchecked(cell, 0U);
          double absolute_pressure = 0.0;
          ThermoState thermo;
          Status evaluated;
          const bool use_density_authority =
              !coupled_pressure_authority &&
              pressure_reference_kind == PressureReferenceKind::closed_mass;
          if (use_density_authority) {
            evaluated = product.thermodynamics.evaluate_from_density(
                density_value, trial_enthalpy.unchecked(cell, 0U),
                {species_values.data(), species_values.size()}, velocity,
                absolute_pressure, thermo,
                temperature_history.accepted.unchecked(cell, 0U));
          } else {
            evaluated =
                product.thermodynamics.evaluate_from_reference_pressure(
                    p_ref, trial_pressure.unchecked(cell, 0U),
                    trial_enthalpy.unchecked(cell, 0U),
                    {species_values.data(), species_values.size()}, velocity,
                    thermo,
                    temperature_history.accepted.unchecked(cell, 0U));
            if (evaluated) {
              absolute_pressure =
                  p_ref + trial_pressure.unchecked(cell, 0U);
              density_value = thermo.rho;
            }
          }
          MolecularTransportState transport;
          if (evaluated)
            evaluated = product.transport.evaluate(
                thermo.temperature,
                {species_values.data(), species_values.size()}, transport);
          if (!evaluated) {
            capture_thermo_failure(cell, p_ref, evaluated);
            return evaluated;
          }
          const double pi = use_density_authority
                                ? absolute_pressure - p_ref
                                : trial_pressure.unchecked(cell, 0U);
          const double oracle_density = thermo.rho;
          if (!std::isfinite(pi) || !std::isfinite(oracle_density) ||
              oracle_density <= 0.0) {
            const Status failure{StatusCode::numerical_failure,
                                 kProductBinding};
            capture_thermo_failure(cell, p_ref, failure);
            return failure;
          }
          if (use_density_authority)
            trial_pressure.unchecked(cell, 0U) = pi;
          if (use_density_authority || coupled_pressure_authority)
            trial_density.unchecked(cell, 0U) = density_value;
          eos_density.unchecked(cell, 0U) = oracle_density;
          trial_temperature.unchecked(cell, 0U) = thermo.temperature;
          compressibility.unchecked(cell, 0U) = thermo.drho_dp_hY;
          enthalpy_compressibility.unchecked(cell, 0U) = thermo.drho_dh_pY;
          heat_capacity.unchecked(cell, 0U) = thermo.cp;
          molecular_viscosity.unchecked(cell, 0U) = transport.viscosity;
          effective_viscosity.unchecked(cell, 0U) = transport.viscosity;
          conductivity.unchecked(cell, 0U) = transport.conductivity;
          enthalpy_diffusivity.unchecked(cell, 0U) =
              transport.conductivity / thermo.cp;
        }
      }
    }
    return Status{};
  };
  double attempt_pressure_reference = pressure_reference;
  const auto pressure_energy_temporal_scale =
      [&](double target_pressure_reference, ConstFieldView target_density,
          ConstFieldView target_enthalpy, ConstFieldView target_pressure,
          Int3 cell, double volume, double& scale) noexcept {
        detail::ProductPressureEnergyTemporalState previous{};
        if (effective_bdf.order == 2U) {
          previous = {
              rho_history.previous.unchecked(cell, 0U),
              enthalpy_history.previous.unchecked(cell, 0U),
              previous_pressure_reference +
                  pressure_history.previous.unchecked(cell, 0U)};
        }
        return detail::product_pressure_energy_temporal_operand_scale(
            effective_bdf, volume,
            {target_density.unchecked(cell, 0U),
             target_enthalpy.unchecked(cell, 0U),
             target_pressure_reference +
                 target_pressure.unchecked(cell, 0U)},
            {rho_history.accepted.unchecked(cell, 0U),
             enthalpy_history.accepted.unchecked(cell, 0U),
             pressure_reference +
                 pressure_history.accepted.unchecked(cell, 0U)},
            previous, scale);
      };
  if (status) {
    begin_timed_stage(12U);
    attempt_stage = 12U;
  }
  if (status && product.contributions.plan() == nullptr)
    status = {StatusCode::invalid_plan, kProductBinding};
  if (status) {
    begin_timed_stage(15U);
    attempt_stage = 15U;
  }
  // A fixed absolute-pressure boundary evaluates its EOS density now but
  // keeps the conservative predictor density paired with the provisional
  // face flux through turbulence and the momentum predictor.  The explicit
  // handoff to the EOS density occurs only after momentum has converged.
  if (status)
    status = update_thermo_from_authority(attempt_pressure_reference, false);
  // normalize_closed_gauge performs a checked_sum after local structural and
  // positivity checks.  Freeze those checks into one rank-consistent
  // prerequisite so no peer can return while another enters the reduction.
  status = product.reductions.consensus(status);
  PressureReferenceCertificate pressure_reference_certificate;
  if (status && product.equations.pressure_reference().kind() ==
                    PressureReferenceKind::closed_mass) {
    status = product.equations.pressure_reference().normalize_closed_gauge(
        trial_pressure, as_const(compressibility), product.reductions,
        attempt_pressure_reference);
    ClosedMassResult closed_result;
    if (status) {
      const ClosedMassDensityFieldView closed_cells{
          as_const(trial_pressure),
          as_const(trial_density),
          as_const(compressibility),
          &product.geometry,
          product.patch,
          product.topology.has_value()
              ? product.topology->region()
              : Span<const std::uint8_t>{},
          predictor_certificate.state};
      status = product.equations.pressure_reference()
                   .certify_closed_mass_density_fields(
                       communicator, product.closed_mass,
                       product.thermodynamics, predictor_certificate,
                       product.equations.pressure_reference().service_stage(),
                       closed_cells, closed_mass_target,
                       attempt_pressure_reference, closed_result,
                       pressure_reference_certificate);
    }
    if (status) {
      attempt_pressure_reference = closed_result.pressure_reference;
      if (!std::isfinite(attempt_pressure_reference) ||
          attempt_pressure_reference <= 0.0) {
        status = {StatusCode::numerical_failure, kProductBinding};
      }
    }
  } else if (status) {
    status = product.equations.pressure_reference().certify_open_boundary(
        attempt_pressure_reference, product.boundary.revision(),
        predictor_certificate, pressure_reference_certificate);
  }
  if (status) {
    for (std::size_t species = 0U; species < species_trial.size(); ++species)
      species_accepted[species] = as_const(species_trial[species]);
    for (std::size_t passive = 0U; passive < passive_trial.size(); ++passive)
      passive_accepted[passive] = as_const(passive_trial[passive]);
  }
  status = resolve_static_boundary_values(
      communicator, product.boundary, product.boundary_specs,
      product.geometry, product.patch,
      product.thermodynamics, attempt_pressure_reference,
      pressure_reference_kind == PressureReferenceKind::boundary_absolute
          ? as_const(eos_density)
          : as_const(trial_density),
      as_const(trial_velocity),
      as_const(trial_enthalpy),
      {species_accepted.data(), species_accepted.size()},
      {passive_accepted.data(), passive_accepted.size()},
      {species_values.data(), species_values.size()}, boundary_scalars,
      boundary_vectors, boundary_normal_gradients, true, status);

  // The effective temporal method is now fixed.  Form its momentum history
  // before the existing stage-20 exchange so predictor, momentum and PISO
  // cannot consume different BDF authorities.  Moving this three-component
  // payload from stage 10 keeps the ordinary-path exchange count and total
  // structured bytes unchanged.
  FieldView temporal_momentum;
  if (status) {
    begin_timed_stage(20U);
    attempt_stage = 20U;
  }
  if (status)
    status = runtime_write_view(product.fields.h_by_a, temporal_momentum);
  for (std::int32_t z = 0; z < cells.z && status; ++z) {
    for (std::int32_t y = 0; y < cells.y && status; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho_n = rho_history.accepted.unchecked(cell, 0U);
        const double rho_nm1 = effective_bdf.order == 2U
                                   ? rho_history.previous.unchecked(cell, 0U)
                                   : 0.0;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double momentum =
              (-effective_bdf.a1 * rho_n *
                   velocity_history.accepted.unchecked(cell, component) -
               effective_bdf.a2 * rho_nm1 *
                   (effective_bdf.order == 2U
                        ? velocity_history.previous.unchecked(cell, component)
                        : 0.0)) /
              effective_bdf.a0;
          if (!std::isfinite(momentum)) {
            status = {StatusCode::numerical_failure, kProductBinding};
            break;
          }
          temporal_momentum.unchecked(cell, component) = momentum;
        }
        if (!status) break;
      }
    }
  }
  if (status) {
    halo_views[0U] = trial_velocity;
    halo_views[1U] = temporal_momentum;
  }
  // The predictor/thermophysical path immediately before stage 20 can report
  // a rank-local numerical status.  Fuse it into the halo's existing
  // preflight so no rank can skip this collective boundary.
  status = exchange(product.stage_halos[2U], 20U, 2U, status);
  if (!status && status.code == StatusCode::numerical_failure) {
    const NumericalFailureContext local_failure = numerical_failure;
    NumericalFailureContext selected_failure;
    const Status published = publish_numerical_failure_context(
        communicator, local_failure, selected_failure);
    if (!published) {
      status = published;
    } else if (selected_failure.valid) {
      numerical_failure = selected_failure;
    }
  }
  const bool momentum_path_active = static_cast<bool>(status);
  if (status) {
    trial_velocity = halo_views[0U];
    temporal_momentum = halo_views[1U];
    status = apply_boundary_ghosts(BoundaryStage::momentum, product.boundary,
                                   {&trial_velocity, 1U}, boundary_values);
  }
  status = product.reductions.consensus(status);
  if (status && product.ibm_gradient_donors.has_value()) {
    std::array<FieldView, 1U> donor_fields{trial_velocity};
    status = product.ibm_gradient_donors->exchange(
        kIbmGradientDonorStage, {donor_fields.data(), donor_fields.size()});
    trial_velocity = donor_fields[0U];
  }
  status = product.reductions.consensus(status);
  FieldView velocity_gradient;
  if (status)
    status = runtime_write_view(product.fields.velocity_gradient,
                                velocity_gradient);
  if (status) {
    const std::array<ConstFieldView, 1U> reads{as_const(trial_velocity)};
    const std::array<FieldView, 1U> writes{velocity_gradient};
    const KernelInvocation gradient_call{
        {reads.data(), reads.size()}, {writes.data(), writes.size()},
        full_box, 0U, 0U, 3U, 0U, nullptr};
    status = cartesian_gradient(product.equations.kernels(), gradient_call);
  }
  if (status && product.ibm_equations.has_value())
    status = product.ibm_equations->correct_velocity_gradient(
        as_const(trial_velocity), velocity_gradient);
  TurbulenceCertificate turbulence_certificate;
  if (status) {
    const TurbulenceUpdateInput turbulence_input{
        as_const(trial_density), as_const(molecular_viscosity), {},
        velocity_gradient.revision, as_const(velocity_gradient)};
    status = product.turbulence.update(turbulence_input, effective_viscosity,
                                       turbulence_certificate);
  }
  if (momentum_path_active)
    status = refresh_live_effective_thermal_ghosts(20U, status);
  const auto publish_momentum_state = [&](StageId halo_stage,
                                           Status prerequisite,
                                           bool enter_collective) {
    if (prerequisite) {
      halo_views[0U] = trial_density;
      halo_views[1U] = trial_velocity;
      halo_views[2U] = trial_pressure;
      halo_views[3U] = effective_viscosity;
      halo_views[4U] = velocity_gradient;
    }
    // grad(U) and turbulence are rank-local, while the momentum halo is an
    // existing collective preflight. Feed a local failure into that preflight
    // whenever peers may already be entering it.
    if (enter_collective)
      prerequisite =
          exchange(product.stage_halos[3U], halo_stage, 5U, prerequisite);
    if (prerequisite) {
      trial_density = halo_views[0U];
      trial_velocity = halo_views[1U];
      trial_pressure = halo_views[2U];
      effective_viscosity = halo_views[3U];
      velocity_gradient = halo_views[4U];
      std::array<FieldView, 4U> zero_gradient_fields{
          trial_density, molecular_viscosity, effective_viscosity,
          velocity_gradient};
      prerequisite = apply_physical_zero_gradient(
          product.boundary,
          {zero_gradient_fields.data(), zero_gradient_fields.size()});
      if (prerequisite)
        prerequisite = apply_boundary_ghosts(
            BoundaryStage::pressure, product.boundary,
            {&trial_pressure, 1U}, boundary_values);
    }
    prerequisite = product.reductions.consensus(prerequisite);
    if (prerequisite && product.ibm_momentum_donors.has_value()) {
      std::array<FieldView, 2U> donor_fields{effective_viscosity,
                                            trial_pressure};
      prerequisite = product.ibm_momentum_donors->exchange(
          130U, {donor_fields.data(), donor_fields.size()});
      effective_viscosity = donor_fields[0U];
      trial_pressure = donor_fields[1U];
    }
    return product.reductions.consensus(prerequisite);
  };
  status = publish_momentum_state(30U, status, momentum_path_active);

  FaceFluxView temporal_flux;
  if (status) {
    const auto divide_temporal = [&](Int3 cell) noexcept {
      const double rho = trial_density.unchecked(cell, 0U);
      if (!std::isfinite(rho) || rho <= 0.0) return false;
      for (std::uint8_t component = 0U; component < 3U; ++component) {
        const double value =
            temporal_momentum.unchecked(cell, component) / rho;
        if (!std::isfinite(value)) return false;
        temporal_momentum.unchecked(cell, component) = value;
      }
      return true;
    };
    for (std::int32_t z = 0; z < cells.z && status; ++z)
      for (std::int32_t y = 0; y < cells.y && status; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          if (!divide_temporal({x, y, z})) {
            status = {StatusCode::numerical_failure, kProductBinding};
            break;
          }
    for (std::size_t face_index = 0U; face_index < 6U && status;
         ++face_index) {
      const auto face = static_cast<CartesianFace>(face_index);
      const BoundaryFacePlan* plan = nullptr;
      if (!product.boundary.face(face, plan) || plan == nullptr) {
        status = {StatusCode::invalid_plan, kProductBinding};
        break;
      }
      if (plan->local_owner && !plan->periodic) continue;
      const bool high = (face_index & 1U) != 0U;
      const CartesianAxis axis =
          face_index < 2U
              ? CartesianAxis::x
              : (face_index < 4U ? CartesianAxis::y : CartesianAxis::z);
      const std::int32_t inner_count =
          axis == CartesianAxis::x ? cells.y : cells.x;
      const std::int32_t outer_count =
          axis == CartesianAxis::z ? cells.y : cells.z;
      for (std::int32_t outer = 0; outer < outer_count && status; ++outer) {
        for (std::int32_t inner = 0; inner < inner_count; ++inner) {
          Int3 ghost{};
          if (axis == CartesianAxis::x)
            ghost = {high ? cells.x : -1, inner, outer};
          else if (axis == CartesianAxis::y)
            ghost = {inner, high ? cells.y : -1, outer};
          else
            ghost = {inner, outer, high ? cells.z : -1};
          if (!divide_temporal(ghost)) {
            status = {StatusCode::numerical_failure, kProductBinding};
            break;
          }
        }
      }
    }
  }
  if (status) {
    FieldView temporal_velocity = temporal_momentum;
    temporal_velocity.field = product.fields.velocity;
    status = apply_boundary_ghosts(BoundaryStage::momentum,
                                   product.boundary,
                                   {&temporal_velocity, 1U}, boundary_values);
    temporal_momentum = temporal_velocity;
    temporal_momentum.field = product.fields.h_by_a;
  }
  if (status) {
    const RevisionToken paired_revision = provisional_flux.revision;
    const RevisionToken temporal_revision =
        paired_revision == std::numeric_limits<RevisionToken>::max()
            ? paired_revision - 1U
            : paired_revision + 1U;
    status = product.phi_workspace.workspace_view(0U, temporal_revision,
                                                  temporal_flux);
  }
  if (status) {
    const std::array<ConstFieldView, 2U> reads{
        as_const(trial_density), as_const(temporal_momentum)};
    const KernelInvocation reconstruct{
        {reads.data(), reads.size()}, {}, full_box, 0U, 0U, 1U, 0U,
        nullptr};
    status = reconstruct_mass_flux(product.equations.kernels(), reconstruct,
                                   temporal_flux);
  }
  if (status && product.ibm_equations.has_value())
    status = product.ibm_equations->zero_interface_flux(temporal_flux);

  FieldView momentum_diagonal;
  FieldView momentum_rhs;
  FieldView momentum_residual;
  FieldView momentum_low_order_rhs_delta;
  FaceFieldView x_coefficient;
  FaceFieldView y_coefficient;
  FaceFieldView z_coefficient;
  EquationStateView equation_state;
  EquationMaterialView material;
  EquationAssemblyContext assembly;
  EquationSystemView momentum_system;
  EquationAssemblyCertificate momentum_certificate;
  const MgDomainActivityView momentum_activity =
      product.topology.has_value()
          ? MgDomainActivityView{
                {product.pressure_mg_cell_activity.data(),
                 product.pressure_mg_cell_activity.size()},
                {product.pressure_mg_x_activity.data(),
                 product.pressure_mg_x_activity.size()},
                {product.pressure_mg_y_activity.data(),
                 product.pressure_mg_y_activity.size()},
                {product.pressure_mg_z_activity.data(),
                 product.pressure_mg_z_activity.size()},
                product.pressure_mg_activity_fingerprint,
                product.pressure_mg_activity_collective}
          : MgDomainActivityView{};
  const auto run_momentum_predictor =
      [&](Status prerequisite, bool time_stage, StageId solve_stage,
          std::uint8_t predictor_pass) {
        if (prerequisite && time_stage) {
          begin_timed_stage(30U);
          attempt_stage = 30U;
        }
        if (prerequisite)
          prerequisite = runtime_write_view(product.fields.momentum_diagonal,
                                            momentum_diagonal);
        if (prerequisite)
          prerequisite = runtime_write_view(product.fields.momentum_rhs,
                                            momentum_rhs);
        if (prerequisite)
          prerequisite = runtime_write_view(product.fields.momentum_residual,
                                            momentum_residual);
        if (prerequisite)
          prerequisite = runtime_write_view(product.fields.h_by_a,
                                            momentum_low_order_rhs_delta);
        if (prerequisite)
          prerequisite = make_pressure_face_views(
              product.pressure_face_storage, cells, x_coefficient,
              y_coefficient, z_coefficient);
        if (prerequisite) {
          prerequisite = history(product.fields.rho, equation_state.density);
          if (prerequisite)
            prerequisite =
                history(product.fields.velocity, equation_state.velocity);
          if (prerequisite)
            prerequisite = history(product.fields.pressure,
                                   equation_state.pressure_perturbation);
          if (prerequisite)
            prerequisite =
                history(product.fields.enthalpy, equation_state.enthalpy);
          if (prerequisite)
            prerequisite = history(product.fields.temperature,
                                   equation_state.temperature);
        }
        equation_state.pressure_reference = attempt_pressure_reference;
        equation_state.accepted_pressure_reference = pressure_reference;
        equation_state.previous_pressure_reference =
            previous_pressure_reference;
        equation_state.independent_species = {species_history.data(),
                                              species_history.size()};
        equation_state.passive_scalars = {passive_history.data(),
                                          passive_history.size()};
        material.molecular_viscosity = as_const(molecular_viscosity);
        material.effective_viscosity = as_const(effective_viscosity);
        material.thermal_conductivity = as_const(conductivity);
        material.heat_capacity = as_const(heat_capacity);
        material.enthalpy_diffusivity = as_const(enthalpy_diffusivity);
        material.pressure_compressibility = as_const(compressibility);
        assembly = {};
        assembly.dt = step.dt;
        assembly.bdf = effective_bdf;
        assembly.time = step.generation;
        assembly.geometry = product.geometry.topology_revision();
        assembly.boundary = product.boundary.revision();
        assembly.thermo = product.thermodynamics.fingerprint();
        assembly.transport = product.transport.fingerprint();
        assembly.face_flux = provisional_flux.revision;
        assembly.contribution_stage = 1U;
        assembly.scope = EquationAssemblyScope::momentum_predictor;
        assembly.mass_flux = as_const(provisional_flux);
        assembly.provisional_mass_flux = true;
        assembly.box = full_box;
        assembly.immersed_interface = product.ibm_equations.has_value()
                                          ? &*product.ibm_equations
                                          : nullptr;
        assembly.wall_treatment = product.ibm_equations.has_value()
                                      ? &product.turbulence
                                      : nullptr;
        momentum_system = {momentum_diagonal, momentum_rhs,
                           momentum_residual, x_coefficient, y_coefficient,
                           z_coefficient};
        momentum_certificate = {};
        if (prerequisite)
          prerequisite = assemble_momentum_predictor(
              product.equations.momentum(), equation_state, material,
              as_const(velocity_gradient), {}, assembly, momentum_system,
              momentum_low_order_rhs_delta, momentum_certificate);
        // Equation kernels are rank-local. Resolve a local failure before the
        // next halo collective so every rank follows the same rollback path.
        prerequisite = product.reductions.consensus(prerequisite);

        std::array<FaceFluxView, 3U> momentum_limiter_faces{};
        FaceFluxView momentum_limiter_alpha;
        const auto limiter_revision = [](RevisionToken base,
                                         RevisionToken offset) noexcept {
          return base <= std::numeric_limits<RevisionToken>::max() - offset
                     ? base + offset
                     : base - offset;
        };
        for (std::size_t component = 0U;
             component < momentum_limiter_faces.size() && prerequisite;
             ++component) {
          const RevisionToken revision = limiter_revision(
              momentum_certificate.state,
              static_cast<RevisionToken>(component + 1U));
          prerequisite = product.phi_workspace.workspace_view(
              2U + component, revision, momentum_limiter_faces[component]);
        }
        if (prerequisite) {
          const RevisionToken alpha_revision =
              limiter_revision(momentum_certificate.state, 4U);
          prerequisite = product.phi_workspace.workspace_view(
              5U, alpha_revision, momentum_limiter_alpha);
        }
        if (prerequisite)
          prerequisite = limit_momentum_predictor_correction(
              communicator, product.equations.momentum(), product.boundary,
              product.patch, momentum_certificate, as_const(trial_velocity),
              as_const(trial_density), step.dt,
              product.time.spec().convective_cfl, as_const(provisional_flux),
              momentum_activity, momentum_system,
              {momentum_low_order_rhs_delta, momentum_limiter_faces,
               momentum_limiter_alpha},
              product.momentum_limiter_halo, product.reductions,
              momentum_predictor_limiter);
        if (prerequisite) {
          const MomentumAdvectiveCflCertificate& cfl =
              momentum_predictor_limiter.advective_cfl;
          if (!cfl.valid() || cfl.plan != momentum_certificate.plan ||
              cfl.time != momentum_certificate.time ||
              cfl.density != trial_density.revision ||
              cfl.density_storage != trial_density.storage_identity ||
              cfl.density_revision_domain != trial_density.revision_domain ||
              !cfl.density_view_identity.matches(as_const(trial_density)) ||
              cfl.face_flux != provisional_flux.revision ||
              cfl.face_flux_storage != provisional_flux.x.storage_identity ||
              cfl.face_flux_revision_domain !=
                  provisional_flux.x.revision_domain ||
              !cfl.face_flux_view_identity.matches(as_const(provisional_flux)) ||
              cfl.activity_collective !=
                  momentum_activity.collective_fingerprint ||
              cfl.dt != momentum_certificate.dt ||
              momentum_certificate.dt != step.dt ||
              cfl.limit != product.time.spec().convective_cfl) {
            prerequisite = {StatusCode::invalid_plan,
                            kProductConvectiveCfl};
          } else {
            prerequisite = convective_cfl_acceptance_status(
                product.time.spec().control, cfl.out_max, cfl.limit);
          }
        }
        // The certificate carries rank-local storage and revision-domain
        // facts. Resolve any mismatch before another rank enters FGMRES.
        prerequisite = product.reductions.consensus(prerequisite);
        if (prerequisite) attempt_stage = solve_stage;
        if (prerequisite)
          prerequisite = solve_momentum_predictor(
              communicator, product.equations.momentum(), product.boundary,
              product.patch, momentum_certificate, as_const(provisional_flux),
              momentum_activity, momentum_system, trial_velocity,
              product.krylov_halo,
              product.piso.pressure_algorithm() == LinearAlgorithm::fgmres
                  ? product.krylov_workspace
                  : product.auxiliary_krylov_workspace,
              product.reductions, &resources, momentum_predictor_solve);
        if (prerequisite)
          prerequisite = transaction.revise_trial(product.fields.velocity);
        if (prerequisite)
          prerequisite = product.layers.view(
              StateRole::trial, product.fields.velocity, trial_velocity);
        if (prerequisite)
          momentum_predictor_solve.predictor_passes = predictor_pass;
        return prerequisite;
      };
  status = run_momentum_predictor(status, true, 31U, 1U);

  // The momentum predictor is complete.  Open domains now hand density
  // authority from the conservative predictor state to the pressure/EOS
  // state before C1.  C1's existing density halo publishes its ghosts, so
  // this handoff is one local scalar copy and no additional communication.
  if (status && pressure_reference_kind ==
                    PressureReferenceKind::boundary_absolute)
    status = transaction.revise_trial(product.fields.rho);
  if (status && pressure_reference_kind ==
                    PressureReferenceKind::boundary_absolute)
    status = product.layers.view(StateRole::trial, product.fields.rho,
                                 trial_density);
  if (status && pressure_reference_kind ==
                    PressureReferenceKind::boundary_absolute)
    status = copy_interior(as_const(eos_density), trial_density);

  // A pressure-energy corrector is a same-target nonlinear state, not a
  // pressure-only post-process.  Publish every primitive/material field that
  // the next continuity and energy residual consumes before asking the
  // coupler for HbyA/phiHbyA.  Effective thermal fields are published only
  // after the resulting velocity ghosts have refreshed turbulence.
  const auto refresh_coupled_state = [&](
      StageId halo_stage, BoundaryThermophysicalGhostPhase phase,
      Status prerequisite = {}) {
    boundary_thermo_certificate = {};
    Status refreshed = prerequisite;
    refreshed = resolve_static_boundary_values(
        communicator, product.boundary, product.boundary_specs,
        product.geometry, product.patch, product.thermodynamics,
        attempt_pressure_reference, as_const(trial_density),
        as_const(trial_velocity), as_const(trial_enthalpy),
        {species_accepted.data(), species_accepted.size()},
        {passive_accepted.data(), passive_accepted.size()},
        {species_values.data(), species_values.size()}, boundary_scalars,
        boundary_vectors, boundary_normal_gradients, true, refreshed);
    halo_views[0U] = trial_velocity;
    halo_views[1U] = trial_density;
    halo_views[2U] = trial_pressure;
    halo_views[3U] = trial_enthalpy;
    halo_views[4U] = trial_temperature;
    refreshed =
        exchange(product.coupled_state_halo, halo_stage, 5U, refreshed);
    if (refreshed) {
      trial_velocity = halo_views[0U];
      trial_density = halo_views[1U];
      trial_pressure = halo_views[2U];
      trial_enthalpy = halo_views[3U];
      trial_temperature = halo_views[4U];
      refreshed = apply_boundary_ghosts(
          BoundaryStage::momentum, product.boundary,
          {&trial_velocity, 1U}, boundary_values);
    }
    if (refreshed)
      refreshed = apply_boundary_ghosts(
          BoundaryStage::pressure, product.boundary,
          {&trial_pressure, 1U}, boundary_values);
    if (refreshed)
      refreshed = apply_boundary_ghosts(
          BoundaryStage::enthalpy, product.boundary,
          {&trial_enthalpy, 1U}, boundary_values);
    if (refreshed) {
      std::vector<BoundaryGhostFieldAuthority>& authorities =
          boundary_thermo_authority_fields;
      if (authorities.size() != species_accepted.size() + 2U) {
        refreshed = {StatusCode::invalid_plan, kProductBinding};
      } else {
        authorities[0U] =
            make_boundary_ghost_field_authority(as_const(trial_pressure));
        authorities[1U] =
            make_boundary_ghost_field_authority(as_const(trial_enthalpy));
        for (std::size_t species = 0U; species < species_accepted.size();
             ++species) {
          const ConstFieldView view = species_accepted[species];
          authorities[species + 2U] =
              make_boundary_ghost_field_authority(view);
        }
        std::uint64_t producer = detail::product_mix(
            product.coupled_state_halo.instance_identity(),
            product.stage_halos[1U].instance_identity());
        producer = detail::product_mix(producer, halo_stage);
        if (producer == 0U) producer = 1U;
        const std::int32_t reach_value = std::max(
            {trial_pressure.ghosts.x, trial_pressure.ghosts.y,
             trial_pressure.ghosts.z});
        if (reach_value <= 0 || reach_value > UINT8_MAX) {
          refreshed = {StatusCode::invalid_plan, kProductBinding};
        } else {
          const BoundaryThermophysicalGhostAuthority authority{
              static_cast<std::uintptr_t>(producer),
              product.boundary.revision(),
              product.boundary.local_layout_fingerprint(),
              cells,
              trial_pressure.ghosts,
              static_cast<std::uint8_t>(reach_value),
              {authorities.data(), authorities.size()}};
          BoundaryThermophysicalGhostCertificate candidate;
          refreshed = BoundaryThermophysicalFaceClosure::close(
              product.boundary, product.thermodynamics, product.transport,
              {attempt_pressure_reference, as_const(trial_pressure),
               as_const(trial_enthalpy),
               {species_accepted.data(), species_accepted.size()},
               authority},
              {trial_density, trial_temperature, heat_capacity,
               compressibility, enthalpy_compressibility,
               molecular_viscosity, conductivity,
               enthalpy_diffusivity},
              {step.generation, product.geometry.fingerprint(),
               pressure_reference_certificate.pressure_reference,
               product.boundary.revision(), phase},
              candidate);
          if (refreshed) boundary_thermo_certificate = candidate;
        }
      }
    }
    refreshed = product.reductions.consensus(refreshed);
    if (!refreshed) return refreshed;
    if (refreshed && product.ibm_gradient_donors.has_value()) {
      std::array<FieldView, 1U> donor_fields{trial_velocity};
      refreshed = product.ibm_gradient_donors->exchange(
          kIbmGradientDonorStage,
          {donor_fields.data(), donor_fields.size()});
      trial_velocity = donor_fields[0U];
    }
    refreshed = product.reductions.consensus(refreshed);
    if (!refreshed) return refreshed;
    if (refreshed && product.ibm_rate_donors.has_value()) {
      halo_count = 0U;
      halo_views[halo_count++] = trial_pressure;
      halo_views[halo_count++] = trial_temperature;
      for (FieldView value : species_trial) halo_views[halo_count++] = value;
      for (FieldView value : passive_trial) halo_views[halo_count++] = value;
      refreshed = product.ibm_rate_donors->exchange(
          161U, {halo_views.data(), halo_count});
      if (refreshed) {
        trial_pressure = halo_views[0U];
        trial_temperature = halo_views[1U];
        species_index = 0U;
        passive_index = 0U;
        for (std::size_t index = 0U;
             index < product.fields.scalars.size(); ++index) {
          if (product.fields.scalar_roles[index] ==
              TransportedScalarRole::species) {
            species_trial[species_index] = halo_views[index + 2U];
            species_accepted[species_index] =
                as_const(species_trial[species_index]);
            ++species_index;
          } else {
            passive_trial[passive_index++] = halo_views[index + 2U];
          }
        }
      }
    }
    refreshed = product.reductions.consensus(refreshed);
    if (!refreshed) return refreshed;
    if (refreshed)
      refreshed = runtime_write_view(product.fields.velocity_gradient,
                                     velocity_gradient);
    if (refreshed) {
      const std::array<ConstFieldView, 1U> reads{as_const(trial_velocity)};
      const std::array<FieldView, 1U> writes{velocity_gradient};
      refreshed = cartesian_gradient(
          product.equations.kernels(),
          {{reads.data(), reads.size()}, {writes.data(), writes.size()},
           full_box, 0U, 0U, 3U, 0U, nullptr});
    }
    if (refreshed && product.ibm_equations.has_value())
      refreshed = product.ibm_equations->correct_velocity_gradient(
          as_const(trial_velocity), velocity_gradient);
    if (refreshed)
      refreshed = runtime_write_view(product.fields.effective_viscosity,
                                     effective_viscosity);
    if (refreshed) {
      const TurbulenceUpdateInput turbulence_input{
          as_const(trial_density), as_const(molecular_viscosity), {},
          velocity_gradient.revision, as_const(velocity_gradient)};
      refreshed = product.turbulence.update(
          turbulence_input, effective_viscosity, turbulence_certificate);
    }
    if (product.transport.kernel() == TransportKernel::coast_native_air)
      return refresh_live_effective_thermal_ghosts(halo_stage, refreshed);
    return product.reductions.consensus(refreshed);
  };

  FieldView pressure_diagonal;
  FieldView pressure_rhs;
  FieldView pressure_correction;
  if (status)
    status = runtime_write_view(product.fields.pressure_diagonal,
                                pressure_diagonal);
  if (status)
    status = runtime_write_view(product.fields.pressure_rhs, pressure_rhs);
  if (status)
    status = runtime_write_view(product.fields.pressure_correction,
                                pressure_correction);
  if (status && !use_pressure_correction_warm_start)
    zero_field(pressure_correction);
  FieldView pressure_energy_c_h;
  FieldView pressure_energy_c_h_row_scale;
  FieldView pressure_energy_e_p;
  FieldView pressure_energy_e_h;
  FieldView pressure_energy_r_c;
  FieldView pressure_energy_r_e;
  FieldView enthalpy_correction;
  FieldView pressure_energy_candidate_pressure;
  FieldView pressure_energy_candidate_pressure_correction;
  FieldView pressure_energy_candidate_enthalpy;
  FieldView pressure_energy_candidate_density;
  FieldView pressure_energy_candidate_temperature;
  FieldView pressure_energy_candidate_velocity;
  FieldView pressure_energy_candidate_molecular_viscosity;
  FieldView pressure_energy_candidate_effective_viscosity;
  FieldView pressure_energy_candidate_velocity_gradient;
  FieldView pressure_energy_candidate_compressibility;
  FieldView pressure_energy_candidate_enthalpy_compressibility;
  FieldView pressure_energy_candidate_thermal_conductivity;
  FieldView pressure_energy_candidate_heat_capacity;
  FieldView pressure_energy_candidate_enthalpy_diffusivity;
  FieldView schur_continuity_response;
  FieldView schur_eliminated_enthalpy;
  FieldView schur_energy_response;
  FieldView pressure_energy_delta_temperature;
  FieldView pressure_energy_compiled_enthalpy_local;
  FieldView pressure_energy_compiled_enthalpy_response;
  if (status)
    status = runtime_write_view(product.fields.pressure_energy_c_h,
                                pressure_energy_c_h);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_c_h_row_scale,
        pressure_energy_c_h_row_scale);
  if (status)
    status = runtime_write_view(product.fields.pressure_energy_e_p,
                                pressure_energy_e_p);
  if (status)
    status = runtime_write_view(product.fields.pressure_energy_e_h,
                                pressure_energy_e_h);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_continuity_residual,
        pressure_energy_r_c);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_energy_residual,
        pressure_energy_r_e);
  if (status)
    status = runtime_write_view(product.fields.enthalpy_correction,
                                enthalpy_correction);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_pressure,
        pressure_energy_candidate_pressure);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_pressure_correction,
        pressure_energy_candidate_pressure_correction);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_enthalpy,
        pressure_energy_candidate_enthalpy);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_density,
        pressure_energy_candidate_density);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_temperature,
        pressure_energy_candidate_temperature);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_velocity,
        pressure_energy_candidate_velocity);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_molecular_viscosity,
        pressure_energy_candidate_molecular_viscosity);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_effective_viscosity,
        pressure_energy_candidate_effective_viscosity);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_velocity_gradient,
        pressure_energy_candidate_velocity_gradient);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_compressibility,
        pressure_energy_candidate_compressibility);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_enthalpy_compressibility,
        pressure_energy_candidate_enthalpy_compressibility);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_thermal_conductivity,
        pressure_energy_candidate_thermal_conductivity);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_heat_capacity,
        pressure_energy_candidate_heat_capacity);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_candidate_enthalpy_diffusivity,
        pressure_energy_candidate_enthalpy_diffusivity);
  if (status)
    status = runtime_write_view(product.fields.schur_continuity_response,
                                schur_continuity_response);
  if (status)
    status = runtime_write_view(product.fields.schur_eliminated_enthalpy,
                                schur_eliminated_enthalpy);
  if (status)
    status = runtime_write_view(product.fields.schur_energy_response,
                                schur_energy_response);
  if (status)
    status = runtime_write_view(
        product.fields.pressure_energy_delta_temperature,
        pressure_energy_delta_temperature);
  if (status)
    status = make_pressure_energy_compiled_cell_views(
        product.energy_compiled_enthalpy_cell_storage, cells,
        product.fields.pressure_energy_e_h, pressure_energy_e_h.revision,
        pressure_energy_compiled_enthalpy_local,
        pressure_energy_compiled_enthalpy_response);
  FaceFieldView energy_assembly_x_coefficient;
  FaceFieldView energy_assembly_y_coefficient;
  FaceFieldView energy_assembly_z_coefficient;
  FaceFieldView energy_frozen_x_enthalpy;
  FaceFieldView energy_frozen_y_enthalpy;
  FaceFieldView energy_frozen_z_enthalpy;
  FaceFieldView energy_directional_x_enthalpy;
  FaceFieldView energy_directional_y_enthalpy;
  FaceFieldView energy_directional_z_enthalpy;
  if (status)
    status = make_pressure_face_views(
        product.energy_assembly_face_storage, cells,
        energy_assembly_x_coefficient, energy_assembly_y_coefficient,
        energy_assembly_z_coefficient);
  if (status)
    status = make_pressure_face_views(
        product.energy_frozen_enthalpy_face_storage, cells,
        energy_frozen_x_enthalpy, energy_frozen_y_enthalpy,
        energy_frozen_z_enthalpy);
  if (status)
    status = make_pressure_face_views(
        product.energy_directional_enthalpy_face_storage, cells,
        energy_directional_x_enthalpy, energy_directional_y_enthalpy,
        energy_directional_z_enthalpy);
  PressureCorrectionSystemView pressure_system{pressure_diagonal,
                                                pressure_rhs};
  EquationSystemView pressure_energy_system{
      pressure_energy_e_h, pressure_energy_r_c, pressure_energy_r_e,
      energy_assembly_x_coefficient, energy_assembly_y_coefficient,
      energy_assembly_z_coefficient};
  const PressureEnergyCellActivity pressure_energy_activity =
      product.topology.has_value()
          ? PressureEnergyCellActivity{
                {product.pressure_mg_cell_activity.data(),
                 product.pressure_mg_cell_activity.size()},
                product.pressure_mg_activity_fingerprint,
                product.pressure_mg_activity_collective}
          : PressureEnergyCellActivity{};
  const PressureContinuityActivityView pressure_energy_continuity_activity =
      product.topology.has_value()
          ? PressureContinuityActivityView{
                {product.pressure_mg_cell_activity.data(),
                 product.pressure_mg_cell_activity.size()},
                {product.pressure_mg_x_activity.data(),
                 product.pressure_mg_x_activity.size()},
                {product.pressure_mg_y_activity.data(),
                 product.pressure_mg_y_activity.size()},
                {product.pressure_mg_z_activity.data(),
                 product.pressure_mg_z_activity.size()},
                product.pressure_mg_activity_fingerprint,
                product.pressure_mg_activity_collective}
          : PressureContinuityActivityView{};
  ConstFaceFluxView pressure_energy_target_flux;
  FrozenConvectionFaceField pressure_energy_frozen_enthalpy;
  const auto revise_pressure_energy_workspaces = [&]() {
    Status revised = runtime_write_view(product.fields.pressure_energy_c_h,
                                        pressure_energy_c_h);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_c_h_row_scale,
          pressure_energy_c_h_row_scale);
    if (revised)
      revised = runtime_write_view(product.fields.pressure_energy_e_p,
                                   pressure_energy_e_p);
    if (revised)
      revised = runtime_write_view(product.fields.pressure_energy_e_h,
                                   pressure_energy_e_h);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_continuity_residual,
          pressure_energy_r_c);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_energy_residual,
          pressure_energy_r_e);
    if (revised)
      revised = runtime_write_view(product.fields.enthalpy_correction,
                                   enthalpy_correction);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_pressure,
          pressure_energy_candidate_pressure);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_pressure_correction,
          pressure_energy_candidate_pressure_correction);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_enthalpy,
          pressure_energy_candidate_enthalpy);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_density,
          pressure_energy_candidate_density);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_temperature,
          pressure_energy_candidate_temperature);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_velocity,
          pressure_energy_candidate_velocity);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_molecular_viscosity,
          pressure_energy_candidate_molecular_viscosity);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_effective_viscosity,
          pressure_energy_candidate_effective_viscosity);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_velocity_gradient,
          pressure_energy_candidate_velocity_gradient);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_compressibility,
          pressure_energy_candidate_compressibility);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_enthalpy_compressibility,
          pressure_energy_candidate_enthalpy_compressibility);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_thermal_conductivity,
          pressure_energy_candidate_thermal_conductivity);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_heat_capacity,
          pressure_energy_candidate_heat_capacity);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_enthalpy_diffusivity,
          pressure_energy_candidate_enthalpy_diffusivity);
    for (std::size_t species = 0U;
         species < product.fields.pressure_energy_candidate_species.size() &&
         revised;
         ++species)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_species[species],
          pressure_energy_candidate_species[species]);
    if (revised)
      revised = runtime_write_view(
          product.fields.schur_continuity_response,
          schur_continuity_response);
    if (revised)
      revised = runtime_write_view(
          product.fields.schur_eliminated_enthalpy,
          schur_eliminated_enthalpy);
    if (revised)
      revised = runtime_write_view(product.fields.schur_energy_response,
                                   schur_energy_response);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_delta_temperature,
          pressure_energy_delta_temperature);
    if (revised)
      revised = make_pressure_energy_compiled_cell_views(
          product.energy_compiled_enthalpy_cell_storage, cells,
          product.fields.pressure_energy_e_h, pressure_energy_e_h.revision,
          pressure_energy_compiled_enthalpy_local,
          pressure_energy_compiled_enthalpy_response);
    if (revised)
      pressure_energy_system = {
          pressure_energy_e_h, pressure_energy_r_c, pressure_energy_r_e,
          energy_assembly_x_coefficient, energy_assembly_y_coefficient,
          energy_assembly_z_coefficient};
    return revised;
  };
  // A globalization sample is an independent numeric authority.  Revise
  // only its cold scratch and the reusable energy residual outputs; the
  // solved dp/dh direction and every live thermophysical workspace remain
  // untouched until a replay has passed the publication gate.
  const auto revise_pressure_energy_candidate_workspaces = [&]() {
    Status revised = runtime_write_view(
        product.fields.pressure_energy_candidate_pressure,
        pressure_energy_candidate_pressure);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_pressure_correction,
          pressure_energy_candidate_pressure_correction);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_enthalpy,
          pressure_energy_candidate_enthalpy);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_density,
          pressure_energy_candidate_density);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_temperature,
          pressure_energy_candidate_temperature);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_velocity,
          pressure_energy_candidate_velocity);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_molecular_viscosity,
          pressure_energy_candidate_molecular_viscosity);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_effective_viscosity,
          pressure_energy_candidate_effective_viscosity);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_velocity_gradient,
          pressure_energy_candidate_velocity_gradient);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_compressibility,
          pressure_energy_candidate_compressibility);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_enthalpy_compressibility,
          pressure_energy_candidate_enthalpy_compressibility);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_thermal_conductivity,
          pressure_energy_candidate_thermal_conductivity);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_heat_capacity,
          pressure_energy_candidate_heat_capacity);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_enthalpy_diffusivity,
          pressure_energy_candidate_enthalpy_diffusivity);
    for (std::size_t species = 0U;
         species < product.fields.pressure_energy_candidate_species.size() &&
         revised;
         ++species)
      revised = runtime_write_view(
          product.fields.pressure_energy_candidate_species[species],
          pressure_energy_candidate_species[species]);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_delta_temperature,
          pressure_energy_delta_temperature);
    if (revised)
      revised = runtime_write_view(product.fields.pressure_energy_e_h,
                                   pressure_energy_e_h);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_continuity_residual,
          pressure_energy_r_c);
    if (revised)
      revised = runtime_write_view(
          product.fields.pressure_energy_energy_residual,
          pressure_energy_r_e);
    if (revised)
      pressure_energy_system = {
          pressure_energy_e_h, pressure_energy_r_c, pressure_energy_r_e,
          energy_assembly_x_coefficient, energy_assembly_y_coefficient,
          energy_assembly_z_coefficient};
    return revised;
  };
  const auto revise_pressure_system_workspaces = [&]() {
    Status revised = runtime_write_view(product.fields.pressure_diagonal,
                                        pressure_diagonal);
    if (revised)
      revised = runtime_write_view(product.fields.pressure_rhs,
                                   pressure_rhs);
    if (revised)
      revised = runtime_write_view(product.fields.pressure_correction,
                                   pressure_correction);
    if (revised)
      pressure_system = {pressure_diagonal, pressure_rhs};
    return revised;
  };
  PisoPressureSolveEpoch solve_epoch;
  if (status) {
    begin_timed_stage(40U);
    attempt_stage = 40U;
  }
  if (status) status = solve_epoch.begin(product.piso);
  if (status)
    status = refresh_coupled_state(
        kCoupledStateC1Stage,
        BoundaryThermophysicalGhostPhase::corrector_one, status);
  const auto assemble_pressure_energy_residual_from_flux =
      [&](ConstFaceFluxView target_flux,
          EquationAssemblyScope scope) {
        Status assembled;
        if (assembled) {
          assembled = history(product.fields.rho, equation_state.density);
          if (assembled)
            assembled = history(product.fields.velocity,
                                equation_state.velocity);
          if (assembled)
            assembled = history(product.fields.pressure,
                                equation_state.pressure_perturbation);
          if (assembled)
            assembled = history(product.fields.enthalpy,
                                equation_state.enthalpy);
          if (assembled)
            assembled = history(product.fields.temperature,
                                equation_state.temperature);
          equation_state.pressure_reference = attempt_pressure_reference;
          equation_state.accepted_pressure_reference = pressure_reference;
          equation_state.previous_pressure_reference =
              previous_pressure_reference;
          equation_state.independent_species = {
              species_history.data(), species_history.size()};
          equation_state.passive_scalars = {
              passive_history.data(), passive_history.size()};
          material.molecular_viscosity = as_const(molecular_viscosity);
          material.effective_viscosity = as_const(effective_viscosity);
          material.thermal_conductivity = as_const(conductivity);
          material.heat_capacity = as_const(heat_capacity);
          material.enthalpy_diffusivity = as_const(enthalpy_diffusivity);
          material.pressure_compressibility = as_const(compressibility);
          assembly.scope = scope;
          assembly.face_flux = target_flux.revision;
          assembly.face_flux_authority =
              scope == EquationAssemblyScope::final_conservative
                  ? target_flux.certificate.authority()
                  : 0U;
          assembly.face_flux_storage =
              scope == EquationAssemblyScope::final_conservative
                  ? target_flux.certificate.storage()
                  : 0U;
          assembly.face_flux_revision_domain =
              scope == EquationAssemblyScope::final_conservative
                  ? target_flux.certificate.revision_domain()
                  : 0U;
          assembly.mass_flux = target_flux;
          assembly.provisional_mass_flux = false;
        }
        EquationAssemblyCertificate energy_certificate;
        if (assembled)
          assembled = assemble_enthalpy(
              product.equations.enthalpy(), equation_state, material,
              as_const(velocity_gradient), {}, assembly,
              pressure_energy_system, energy_certificate);
        if (assembled && product.ibm_equations.has_value()) {
          zero_field(pressure_energy_e_p);
          assembled = product.ibm_equations->correct_pressure_work(
              as_const(trial_pressure), as_const(trial_velocity),
              pressure_energy_e_p);
          for (std::int32_t z = 0; z < cells.z && assembled; ++z)
            for (std::int32_t y = 0; y < cells.y; ++y)
              for (std::int32_t x = 0; x < cells.x; ++x) {
                const Int3 cell{x, y, z};
                pressure_energy_r_e.unchecked(cell, 0U) -=
                    detail::cell_volume(product.equations.kernels(), cell) *
                    pressure_energy_e_p.unchecked(cell, 0U);
              }
          if (assembled) {
            zero_field(pressure_energy_e_p);
            assembled =
                product.ibm_equations->correct_zero_normal_diffusion(
                    as_const(trial_temperature), as_const(conductivity),
                    pressure_energy_e_p);
          }
          for (std::int32_t z = 0; z < cells.z && assembled; ++z)
            for (std::int32_t y = 0; y < cells.y; ++y)
              for (std::int32_t x = 0; x < cells.x; ++x) {
                const Int3 cell{x, y, z};
                pressure_energy_r_e.unchecked(cell, 0U) -=
                    detail::cell_volume(product.equations.kernels(), cell) *
                    pressure_energy_e_p.unchecked(cell, 0U);
              }
        }
        if (assembled && product.topology.has_value()) {
          const Span<const std::uint8_t> active =
              product.topology->region();
          std::size_t offset = 0U;
          for (std::int32_t z = 0; z < cells.z; ++z)
            for (std::int32_t y = 0; y < cells.y; ++y)
              for (std::int32_t x = 0; x < cells.x; ++x, ++offset)
                if (active.data[offset] == 0U)
                  pressure_energy_r_e.unchecked({x, y, z}, 0U) = 0.0;
        }
        return product.reductions.consensus(assembled);
      };
  const auto inspect_pressure_energy_target_flux =
      [&](const PisoIntermediateCertificate& intermediate) {
        ConstFaceFluxView target_flux;
        Status inspected = product.coupler.inspect_intermediate_flux(
            intermediate, target_flux);
        inspected = product.reductions.consensus(inspected);
        if (inspected) pressure_energy_target_flux = target_flux;
        return inspected;
      };
  const auto assemble_pressure_energy_residual =
      [&](const PisoIntermediateCertificate& intermediate) {
        Status inspected = inspect_pressure_energy_target_flux(intermediate);
        if (inspected)
          inspected = assemble_pressure_energy_residual_from_flux(
              pressure_energy_target_flux,
              EquationAssemblyScope::target_coupled);
        FrozenConvectionFaceField frozen;
        if (inspected)
          inspected = freeze_cartesian_target_convection_faces(
              product.equations.kernels(),
              product.equations.thermophysical_predictor()
                  .enthalpy_convection(),
              pressure_energy_target_flux, as_const(trial_enthalpy), 0U,
              {product.equations.semantic_fingerprint(),
               product.boundary.revision()},
              {energy_frozen_x_enthalpy, energy_frozen_y_enthalpy,
               energy_frozen_z_enthalpy},
              frozen);
        inspected = product.reductions.consensus(inspected);
        if (inspected) {
          pressure_energy_frozen_enthalpy = frozen;
        }
        return inspected;
      };
  PisoIntermediateInput intermediate_input;
  intermediate_input.momentum = momentum_certificate;
  intermediate_input.predictor = predictor_certificate;
  intermediate_input.pressure_reference = pressure_reference_certificate;
  intermediate_input.density = trial_density;
  intermediate_input.trial_velocity = as_const(trial_velocity);
  intermediate_input.trial_flux = as_const(provisional_flux);
  intermediate_input.momentum_system = momentum_system;
  intermediate_input.bdf = effective_bdf;
  intermediate_input.numeric_boundary = product.boundary.revision();
  intermediate_input.boundary_values = boundary_values;
  intermediate_input.corrector = 1U;
  intermediate_input.immersed_interface = product.ibm_equations.has_value()
                                              ? &*product.ibm_equations
                                              : nullptr;
  intermediate_input.temporal_reference = as_const(temporal_flux);
  intermediate_input.committed_face_history.accepted = accepted_flux;
  intermediate_input.committed_face_history.previous =
      effective_bdf.order == 2U ? previous_flux : ConstFaceFluxView{};
  intermediate_input.thermophysical_boundary = {
      boundary_thermo_certificate,
      {attempt_pressure_reference, as_const(trial_pressure),
       as_const(trial_enthalpy),
       {species_accepted.data(), species_accepted.size()},
       as_const(trial_density)}};
  PisoIntermediateCertificate intermediate_one;
  if (status) attempt_stage = 41U;
  if (status)
    status = product.coupler.refresh(intermediate_input, intermediate_one);
  if (status)
    status = product.piso.coupling() == CouplingKind::simple
                 ? inspect_pressure_energy_target_flux(intermediate_one)
                 : assemble_pressure_energy_residual(intermediate_one);
  PressureCorrectionInput pressure_input;
  pressure_input.intermediate = intermediate_one;
  pressure_input.pressure_reference = pressure_reference_certificate;
  pressure_input.density_trial = as_const(trial_density);
  pressure_input.density_accepted = rho_history.accepted;
  pressure_input.density_previous =
      effective_bdf.order == 2U ? rho_history.previous : ConstFieldView{};
  pressure_input.drho_dp_h_y = as_const(compressibility);
  pressure_input.bdf = effective_bdf;
  pressure_input.time = step.generation;
  pressure_input.geometry = product.geometry.topology_revision();
  pressure_input.numeric_boundary = product.boundary.revision();
  PressureCorrectionCertificate pressure_one;
  if (status) attempt_stage = 42U;
  if (status)
    status = product.coupler.assemble_pressure_system(
        pressure_input, pressure_system, pressure_one);
  const auto linear_identity = [&](const PressureCorrectionCertificate& p) {
    LinearIdentity identity;
    // Symbolic/hierarchy are rank-invariant semantic authorities. Numeric and
    // workspace identities remain rank-local and are validated locally by the
    // operator/preconditioner lifecycle.
    identity.symbolic = product.piso.fingerprint();
    identity.numeric = p.state;
    identity.hierarchy = product.geometry.fingerprint();
    identity.workspace = product.krylov_workspace.fingerprint();
    std::uint64_t fingerprint = kFnvOffset;
    fingerprint = detail::product_mix(fingerprint, identity.symbolic);
    fingerprint = detail::product_mix(fingerprint, identity.numeric);
    fingerprint = detail::product_mix(fingerprint, identity.hierarchy);
    fingerprint = detail::product_mix(fingerprint, identity.workspace);
    identity.fingerprint = fingerprint == 0U ? 1U : fingerprint;
    return identity;
  };
  const auto pressure_energy_active = [&](std::size_t offset) noexcept {
    return !product.topology.has_value() ||
           product.pressure_mg_cell_activity[offset] != 0U;
  };
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  struct PressureEnergyLinearPredictionDiagnostic {
    double normalized_continuity{std::numeric_limits<double>::infinity()};
    double normalized_energy{std::numeric_limits<double>::infinity()};
    bool valid{};
  };
  std::array<PressureEnergyLinearPredictionDiagnostic, 2U>
      pressure_energy_linear_predictions{};
#endif
  const auto form_pressure_energy_blocks = [&]() {
    Status formed;
    std::size_t offset = 0U;
    for (std::int32_t z = 0; z < cells.z && formed; ++z) {
      for (std::int32_t y = 0; y < cells.y && formed; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x, ++offset) {
          const Int3 cell{x, y, z};
          if (!pressure_energy_active(offset)) {
            pressure_energy_c_h.unchecked(cell, 0U) = 0.0;
            pressure_energy_c_h_row_scale.unchecked(cell, 0U) = 1.0;
            pressure_energy_e_p.unchecked(cell, 0U) = 0.0;
            pressure_energy_e_h.unchecked(cell, 0U) = 1.0;
            pressure_energy_r_c.unchecked(cell, 0U) = 0.0;
            pressure_energy_r_e.unchecked(cell, 0U) = 0.0;
            continue;
          }
          const double volume =
              detail::cell_volume(product.equations.kernels(), cell);
          const double rho = trial_density.unchecked(cell, 0U);
          const double h = trial_enthalpy.unchecked(cell, 0U);
          const double rho_p = compressibility.unchecked(cell, 0U);
          const double rho_h =
              enthalpy_compressibility.unchecked(cell, 0U);
          const double c_h = effective_bdf.a0 * volume * rho_h;
          const double temporal_reference =
              std::abs(effective_bdf.a0 * volume * rho /
                       std::max(std::abs(h), 1.0));
          const double row_scale =
              std::max(std::abs(c_h), temporal_reference);
          const double e_p =
              effective_bdf.a0 * volume * (h * rho_p - 1.0);
          const double e_h = pressure_energy_e_h.unchecked(cell, 0U);
          const double continuity_residual =
              -pressure_rhs.unchecked(cell, 0U);
          if (!std::isfinite(volume) || !(volume > 0.0) ||
              !std::isfinite(rho) || !(rho > 0.0) ||
              !std::isfinite(h) || !std::isfinite(rho_p) ||
              !std::isfinite(rho_h) || !std::isfinite(c_h) ||
              !std::isfinite(row_scale) || !(row_scale > 0.0) ||
              std::abs(c_h) / row_scale < 1.0e-12 ||
              !std::isfinite(e_p) || !std::isfinite(e_h) ||
              !std::isfinite(continuity_residual) ||
              !std::isfinite(pressure_energy_r_e.unchecked(cell, 0U))) {
            formed = {StatusCode::rejected_step,
                      kProductPressureEnergy};
            break;
          }
          pressure_energy_c_h.unchecked(cell, 0U) = c_h;
          pressure_energy_c_h_row_scale.unchecked(cell, 0U) = row_scale;
          pressure_energy_e_p.unchecked(cell, 0U) = e_p;
          pressure_energy_e_h.unchecked(cell, 0U) = e_h;
          pressure_energy_r_c.unchecked(cell, 0U) = continuity_residual;
        }
      }
    }
    return product.reductions.consensus(formed);
  };
  struct PressureEnergyJacobianObservation {
    bool valid{};
    PressureEnergyJacobianScope scope{
        PressureEnergyJacobianScope::generic_algebraic_quasi_newton};
  };
  const auto solve_pressure_energy =
      [&](std::uint8_t corrector,
          std::uint8_t refinement_iteration,
          const PisoIntermediateCertificate& intermediate,
          const PressureCorrectionCertificate& pressure,
          ConstFaceFluxView target_flux,
          const FrozenConvectionFaceField& frozen_enthalpy,
          LinearIdentity identity, MgCoefficientIdentity coefficients,
          const LinearSolveControl& solve_control,
          PressureEnergyJacobianObservation& jacobian_observation) {
        jacobian_observation = {};
        Status coupled =
            refinement_iteration == 0U
                ? solve_epoch.prepare_linear_lifecycle(
                      product.piso, corrector, pressure, identity,
                      coefficients, product.coupler, pressure_operator,
                      pressure_mg, pressure_system, pressure_correction,
                      product.krylov_workspace, &pressure_mg_counters)
                : solve_epoch.prepare_pressure_energy_refinement_lifecycle(
                      product.piso, refinement_iteration, pressure, identity,
                      coefficients, product.coupler, pressure_operator,
                      pressure_mg, pressure_system, pressure_correction,
                      product.krylov_workspace, &pressure_mg_counters);
        PressureEnergyPressureFluxOperator spatial_energy_pressure_operator;
        PressureEnergyEnthalpyOperator spatial_energy_enthalpy_operator;
        PressureEnergyDiagonalOperator diagonal_energy_pressure_operator;
        PressureEnergyDiagonalOperator diagonal_energy_enthalpy_operator;
        LinearOperator* energy_pressure_operator = nullptr;
        LinearOperator* energy_enthalpy_operator = nullptr;
        PressureEnergyPressureFluxCertificate energy_pressure_certificate;
        PressureEnergyEnthalpyCertificate energy_enthalpy_certificate;
        PressureEnergyDiagonalCertificate diagonal_energy_pressure_certificate;
        PressureEnergyDiagonalCertificate diagonal_energy_enthalpy_certificate;
        PressureEnergySchurBlockAuthority schur_block_authority;
        PressureEnergySchurOperator schur_operator;
        PressureEnergyJacobianCertificate jacobian_certificate;
        PressureEnergyFrozenFaceEnthalpy pressure_energy_enthalpy{
            frozen_enthalpy.x, frozen_enthalpy.y, frozen_enthalpy.z,
            frozen_enthalpy.revision, frozen_enthalpy.reconstruction, 0U};
        pressure_energy_enthalpy.local_binding =
            pressure_energy_frozen_face_enthalpy_local_binding(
                pressure_energy_enthalpy);
        PressureEnergyEnthalpyBinding energy_enthalpy_binding;
        energy_enthalpy_binding.geometry = &product.geometry;
        energy_enthalpy_binding.kernels = &product.equations.kernels();
        energy_enthalpy_binding.boundary = &product.boundary;
        energy_enthalpy_binding.patch = product.patch;
        energy_enthalpy_binding.convection =
            product.equations.thermophysical_predictor()
                .enthalpy_convection();
        energy_enthalpy_binding.services = {
            communicator, &product.pressure_energy_enthalpy_halo, 141U,
            product.fields.schur_eliminated_enthalpy,
            product.fields.pressure_energy_delta_temperature};
        energy_enthalpy_binding.authority = {
            effective_bdf,
            step.generation,
            product.geometry.topology_revision(),
            product.boundary.revision(),
            product.equations.thermodynamics_fingerprint(),
            product.equations.transport_fingerprint(),
            product.equations.semantic_fingerprint(),
            product.equations.thermodynamics_fingerprint(),
            product.equations.transport_fingerprint()};
        energy_enthalpy_binding.assembled_diagonal =
            as_const(pressure_energy_e_h);
        energy_enthalpy_binding.target_enthalpy =
            as_const(trial_enthalpy);
        energy_enthalpy_binding.density_enthalpy_derivative =
            as_const(enthalpy_compressibility);
        energy_enthalpy_binding.heat_capacity = as_const(heat_capacity);
        energy_enthalpy_binding.thermal_conductivity =
            as_const(conductivity);
        energy_enthalpy_binding.enthalpy_diffusivity =
            as_const(enthalpy_diffusivity);
        energy_enthalpy_binding.target_flux = target_flux;
        energy_enthalpy_binding.convection_context = {
            product.equations.semantic_fingerprint(),
            product.boundary.revision()};
        energy_enthalpy_binding.frozen_face_enthalpy = frozen_enthalpy;
        energy_enthalpy_binding.workspace = {
            pressure_energy_delta_temperature,
            {energy_directional_x_enthalpy,
             energy_directional_y_enthalpy,
             energy_directional_z_enthalpy},
            {// The enthalpy assembly face coefficients are dead throughout
             // this linear solve.  Reuse their frozen capacity for exact
             // lambda transmissibilities; every nonlinear candidate replay
             // reassembles its own coefficients before consuming them.
             pressure_energy_compiled_enthalpy_local,
             {energy_assembly_x_coefficient,
              energy_assembly_y_coefficient,
              energy_assembly_z_coefficient},
             {{product.energy_compiled_enthalpy_branch_storage.data(),
               product.energy_compiled_enthalpy_branch_storage.size()}},
             pressure_energy_compiled_enthalpy_response}};
        energy_enthalpy_binding.activity =
            pressure_energy_continuity_activity;
        energy_enthalpy_binding.identity = identity;
        energy_enthalpy_binding.linearization_policy =
            FrozenConvectionLinearizationPolicy::
                semismooth_generalized_zero_slope;
        if (coupled && product.piso.coupling() == CouplingKind::simple &&
            corrector == 1U && refinement_iteration == 0U) {
          LinearOperator& exact_pressure_operator =
              ibm_pressure_operator.has_value()
                  ? static_cast<LinearOperator&>(*ibm_pressure_operator)
                  : static_cast<LinearOperator&>(pressure_operator);
          zero_field(enthalpy_correction);
          LinearSolveControl simple_solve_control = solve_control;
          if (simple_solve_control.relative_tolerance >
              product.piso.pressure_solve().relative_tolerance)
            simple_solve_control.relative_tolerance =
                kSimplePressureInexactForcingRelativeToleranceCeiling;
          coupled = solve_epoch.solve_prepared(
              exact_pressure_operator,
              PisoPressureSolveContract::pressure_continuity,
              simple_solve_control, product.reductions, &resources);
          return product.reductions.consensus(coupled);
        }
        const bool simple_diagonal_schur =
            detail::product_use_simple_diagonal_schur(
                product.piso.coupling(), corrector, refinement_iteration,
                product.ibm_equations.has_value());
        PisoCartesianPressureWorkLinearization pressure_work;
        if (coupled && !simple_diagonal_schur &&
            !product.ibm_equations.has_value())
          coupled = product.coupler
                        .inspect_cartesian_pressure_work_linearization(
                            intermediate, pressure, as_const(trial_pressure),
                            as_const(trial_velocity), pressure_work);
        PressureEnergyPressureFluxBinding pressure_flux_binding;
        pressure_flux_binding.geometry = &product.geometry;
        pressure_flux_binding.boundary = &product.boundary;
        pressure_flux_binding.patch = product.patch;
        pressure_flux_binding.services = {
            communicator, &product.krylov_halo, 140U,
            product.fields.krylov_vectors, 1U};
        pressure_flux_binding.intermediate = intermediate;
        pressure_flux_binding.pressure = pressure;
        pressure_flux_binding.temporal_diagonal =
            as_const(pressure_energy_e_p);
        pressure_flux_binding.x_pressure_coefficient =
            as_const(x_coefficient);
        pressure_flux_binding.y_pressure_coefficient =
            as_const(y_coefficient);
        pressure_flux_binding.z_pressure_coefficient =
            as_const(z_coefficient);
        pressure_flux_binding.target_flux = target_flux;
        pressure_flux_binding.frozen_face_enthalpy =
            pressure_energy_enthalpy;
        pressure_flux_binding.activity =
            pressure_energy_continuity_activity;
        pressure_flux_binding.identity = identity;
        // IBM intentionally leaves pressure_work empty: its Cartesian
        // cell/face response is a flux-only spatial quasi-Newton block.  No
        // IBM donor-gradient derivative is claimed by either typed operator.
        pressure_flux_binding.pressure_work = pressure_work;
        if (simple_diagonal_schur) {
          if (coupled)
            coupled = PressureEnergyDiagonalOperator::bind(
                {as_const(pressure_energy_e_p), pressure_energy_activity,
                 identity, 0.0},
                diagonal_energy_pressure_operator,
                diagonal_energy_pressure_certificate);
          energy_pressure_operator = &diagonal_energy_pressure_operator;
          if (coupled)
            coupled = PressureEnergyDiagonalOperator::bind(
                {as_const(pressure_energy_e_h), pressure_energy_activity,
                 identity, 1.0},
                diagonal_energy_enthalpy_operator,
                diagonal_energy_enthalpy_certificate);
          energy_enthalpy_operator = &diagonal_energy_enthalpy_operator;
        } else {
          if (coupled)
            coupled = PressureEnergyPressureFluxOperator::bind(
                pressure_flux_binding, spatial_energy_pressure_operator,
                energy_pressure_certificate);
          energy_pressure_operator = &spatial_energy_pressure_operator;
          if (coupled)
            coupled = PressureEnergyEnthalpyOperator::bind(
                energy_enthalpy_binding, spatial_energy_enthalpy_operator,
                energy_enthalpy_certificate);
          energy_enthalpy_operator = &spatial_energy_enthalpy_operator;
        }
        if (coupled && simple_diagonal_schur)
          coupled = PressureEnergySchurBlockAuthority::ibm_double_diagonal(
              diagonal_energy_pressure_operator,
              diagonal_energy_enthalpy_operator, schur_block_authority);
        else if (coupled && product.ibm_equations.has_value())
          coupled = PressureEnergySchurBlockAuthority::
              ibm_cartesian_spatial_quasi_newton(
                  spatial_energy_pressure_operator,
                  spatial_energy_enthalpy_operator, schur_block_authority);
        else if (coupled)
          coupled = PressureEnergySchurBlockAuthority::exact_cartesian(
              spatial_energy_pressure_operator,
              spatial_energy_enthalpy_operator, schur_block_authority);
        LinearOperator* continuity_pressure_operator =
            ibm_pressure_operator.has_value()
                ? static_cast<LinearOperator*>(&*ibm_pressure_operator)
                : static_cast<LinearOperator*>(&pressure_operator);
        if (coupled)
          coupled = PressureEnergySchurOperator::bind(
              {continuity_pressure_operator,
               energy_pressure_operator,
               energy_enthalpy_operator,
               as_const(pressure_energy_c_h),
               as_const(pressure_energy_c_h_row_scale),
               {schur_continuity_response, schur_eliminated_enthalpy,
                schur_energy_response},
               pressure_energy_activity,
               1.0e-12,
               schur_block_authority},
              schur_operator, jacobian_certificate);
        if (coupled) {
          const PressureEnergyJacobianScope expected_scope =
              simple_diagonal_schur
                  ? PressureEnergyJacobianScope::
                        ibm_double_diagonal_quasi_newton
              : product.ibm_equations.has_value()
                  ? PressureEnergyJacobianScope::
                        ibm_cartesian_spatial_quasi_newton
                  : PressureEnergyJacobianScope::
                        exact_cartesian_frozen_spatial;
          if (!jacobian_certificate.valid() ||
              jacobian_certificate.jacobian_scope != expected_scope ||
              jacobian_certificate.full_nonlinear_jacobian)
            coupled = {StatusCode::invalid_plan, kProductPressureEnergy};
        }
        coupled = product.reductions.consensus(coupled);
        if (coupled)
          coupled = schur_operator.form_pressure_rhs(
              as_const(pressure_energy_r_c),
              as_const(pressure_energy_r_e), pressure_rhs);
        coupled = product.reductions.consensus(coupled);
        PressureEnergySchurPreparedApplyEpoch repeated_apply_epoch;
        if (coupled && !simple_diagonal_schur)
          coupled =
              schur_operator.prepare_repeated_apply(repeated_apply_epoch);
        coupled = product.reductions.consensus(coupled);
        const int prepare_lowest_failing_rank =
            product.reductions.lowest_failing_rank();
        const bool solve_invoked = static_cast<bool>(coupled);
        int close_lowest_failing_rank = prepare_lowest_failing_rank;
        if (solve_invoked) {
          jacobian_observation.scope = jacobian_certificate.jacobian_scope;
          jacobian_observation.valid = true;
          coupled = solve_epoch.solve_prepared(
              schur_operator,
              PisoPressureSolveContract::continuity_energy_coupled,
              solve_control, product.reductions, &resources);
          // This accessor is written directly from the LinearSolveResult and
          // avoids copying the full PisoAttemptReport on every hot solve.  An
          // invoked solve without a concrete rank deliberately closes with
          // -1 rather than manufacturing rank zero from a second consensus.
          close_lowest_failing_rank =
              solve_epoch.latest_solve_outcome_available()
                  ? solve_epoch.latest_solve_lowest_failing_rank()
                  : -1;
        }
        // Prepared Halo apply deliberately defers rank-local failures to the
        // Krylov checked reductions.  Publish one explicit all-rank solve
        // outcome before either Halo epoch can publish ghost revisions.
        if (solve_invoked)
          coupled = product.reductions.consensus(coupled);
        Status closed;
        if (repeated_apply_epoch.valid())
          closed = schur_operator.close_repeated_apply(
              repeated_apply_epoch, coupled, close_lowest_failing_rank);
        // Preserve a failed solve; otherwise publish even a rank-local close
        // anomaly before any rank can enter ordinary recovery.  This avoids a
        // conditional-collective split on the exceptional close path.
        coupled = product.reductions.consensus(coupled ? closed : coupled);
        if (coupled) {
          FieldView pressure_correction_for_operator = pressure_correction;
          pressure_correction_for_operator.field =
              product.fields.krylov_vectors;
          coupled = schur_operator.recover_enthalpy(
              as_const(pressure_energy_r_c),
              pressure_correction_for_operator, enthalpy_correction);
        }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        const bool diagnostic_active =
            g_candidate_globalization_armed.load(std::memory_order_acquire) ||
            g_candidate_globalization_published.load(
                std::memory_order_acquire);
        if (coupled && diagnostic_active &&
            !product.ibm_equations.has_value() && corrector >= 1U &&
            corrector <= pressure_energy_linear_predictions.size()) {
          // Evaluate the two frozen Jacobian rows on the recovered joint
          // direction.  This is diagnostic-only, but it deliberately uses the
          // bound production operators so a sign, field-identity, or target-
          // layer error cannot be hidden by a hand-written proxy.
          FieldView pressure_direction = pressure_correction;
          pressure_direction.field = product.fields.krylov_vectors;
          coupled = energy_pressure_operator->apply(
              pressure_direction,
              pressure_energy_candidate_pressure_correction);
          if (coupled) {
            for (std::int32_t z = 0; z < cells.z; ++z)
              for (std::int32_t y = 0; y < cells.y; ++y)
                for (std::int32_t x = 0; x < cells.x; ++x)
                  schur_eliminated_enthalpy.unchecked({x, y, z}, 0U) =
                      enthalpy_correction.unchecked({x, y, z}, 0U);
            coupled = energy_enthalpy_operator->apply(
                schur_eliminated_enthalpy, schur_energy_response);
          }
          double local_prediction[2U]{};
          std::size_t offset = 0U;
          for (std::int32_t z = 0; z < cells.z && coupled; ++z) {
            for (std::int32_t y = 0; y < cells.y && coupled; ++y) {
              for (std::int32_t x = 0; x < cells.x; ++x, ++offset) {
                if (!pressure_energy_active(offset)) continue;
                const Int3 cell{x, y, z};
                const double volume =
                    detail::cell_volume(product.equations.kernels(), cell);
                const double rho = trial_density.unchecked(cell, 0U);
                const double rho_n =
                    rho_history.accepted.unchecked(cell, 0U);
                const double rho_nm1 =
                    effective_bdf.order == 2U
                        ? rho_history.previous.unchecked(cell, 0U)
                        : 0.0;
                const double continuity_scale = std::max(
                    std::numeric_limits<double>::min(),
                    std::abs(volume * effective_bdf.a0 * rho) +
                        std::abs(volume * effective_bdf.a1 * rho_n) +
                        std::abs(volume * effective_bdf.a2 * rho_nm1) +
                        std::abs(target_flux.x.unchecked(cell)) +
                        std::abs(target_flux.x.unchecked({x + 1, y, z})) +
                        std::abs(target_flux.y.unchecked(cell)) +
                        std::abs(target_flux.y.unchecked({x, y + 1, z})) +
                        std::abs(target_flux.z.unchecked(cell)) +
                        std::abs(target_flux.z.unchecked({x, y, z + 1})));
                const double continuity_linear =
                    pressure_energy_r_c.unchecked(cell, 0U) +
                    schur_continuity_response.unchecked(cell, 0U) +
                    pressure_energy_c_h.unchecked(cell, 0U) *
                        enthalpy_correction.unchecked(cell, 0U);
                const double enthalpy =
                    trial_enthalpy.unchecked(cell, 0U);
                const double temperature =
                    trial_temperature.unchecked(cell, 0U);
                const double cp = heat_capacity.unchecked(cell, 0U);
                const double absolute_pressure = attempt_pressure_reference +
                                                 trial_pressure.unchecked(
                                                     cell, 0U);
                double energy_scale = 0.0;
                const bool energy_scale_valid =
                    pressure_energy_temporal_scale(
                        attempt_pressure_reference, as_const(trial_density),
                        as_const(trial_enthalpy), as_const(trial_pressure),
                        cell, volume, energy_scale);
                const double energy_linear =
                    pressure_energy_r_e.unchecked(cell, 0U) +
                    pressure_energy_candidate_pressure_correction.unchecked(
                        cell, 0U) +
                    schur_energy_response.unchecked(cell, 0U);
                if (!std::isfinite(continuity_scale) ||
                    !(continuity_scale > 0.0) ||
                    !std::isfinite(continuity_linear) ||
                    !std::isfinite(enthalpy) ||
                    !std::isfinite(temperature) || !(temperature > 0.0) ||
                    !std::isfinite(cp) || !(cp > 0.0) ||
                    !std::isfinite(absolute_pressure) ||
                    !(absolute_pressure > 0.0) || !energy_scale_valid ||
                    !std::isfinite(energy_linear)) {
                  coupled = {StatusCode::numerical_failure,
                             kProductPressureEnergy};
                  break;
                }
                local_prediction[0U] = std::max(
                    local_prediction[0U],
                    std::abs(continuity_linear) / continuity_scale);
                local_prediction[1U] =
                    std::max(local_prediction[1U],
                             std::abs(energy_linear) / energy_scale);
              }
            }
          }
          double global_prediction[2U]{};
          coupled = product.reductions.checked_max(
              {local_prediction, 2U}, {global_prediction, 2U}, coupled);
          if (coupled) {
            PressureEnergyLinearPredictionDiagnostic& prediction =
                pressure_energy_linear_predictions[corrector - 1U];
            prediction.normalized_continuity = global_prediction[0U];
            prediction.normalized_energy = global_prediction[1U];
            prediction.valid = true;
          }
        }
#endif
        return product.reductions.consensus(coupled);
      };
  struct PressureEnergyCandidateArtifacts {
    PressureEnergyGlobalizationSample sample{};
    PisoFrozenMomentumPressureStageCertificate pressure_stage{};
    PisoFrozenMomentumVelocityStageCertificate velocity_stage{};
    PisoFrozenMomentumFluxStageCertificate flux_stage{};
    PisoFrozenMomentumExactCandidateCertificate exact_certificate{};
    BoundaryThermophysicalGhostCertificate boundary_thermophysics{};
    FinalBoundaryFluxCertificate final_boundary{};
    EquationAssemblyCertificate energy{};
    PisoExactThermodynamicCandidateView exact{};
    ConstFaceFluxView flux{};
    double pressure_reference{};
    double alpha_zero_pressure_oracle_difference{};
    double alpha_zero_pressure_oracle_relative_error{};
    double alpha_zero_density_oracle_difference{};
    double alpha_zero_density_oracle_relative_error{};
    std::uint64_t alpha_zero_density_oracle_ulp{};
    double alpha_zero_temperature_oracle_difference{};
    double alpha_zero_temperature_oracle_relative_error{};
    std::uint64_t alpha_zero_temperature_oracle_ulp{};
    double alpha_zero_material_oracle_error{};
    double alpha_zero_gradient_oracle_error{};
    double alpha_zero_energy_residual_oracle_error{};
    double alpha_zero_effective_viscosity_oracle_error{};
    PlanFingerprint alpha_zero_oracle_numeric_lineage{};
    PlanFingerprint residual_replay_provenance{};
    bool alpha_zero_byte_equivalent{};
  };
  const auto semantic_field = [](ConstFieldView view, FieldId field) noexcept {
    view.field = field;
    return view;
  };
  const auto prepare_exact_thermodynamic_candidate =
      [&](const PressureCorrectionCertificate& pressure,
          ConstFieldView candidate_pressure_correction,
          ConstFieldView candidate_enthalpy_correction,
          ConstFieldView candidate_enthalpy,
          ConstFieldView candidate_density,
          ConstFieldView candidate_temperature,
          ConstFieldView candidate_pressure_compressibility,
          Span<const ConstFieldView> candidate_species,
          Span<const ConstFieldView> semantic_species,
          PisoExactThermodynamicCandidateView& candidate) {
        candidate = {};
        if (candidate_species.size != semantic_species.size ||
            (candidate_species.size != 0U &&
             (candidate_species.data == nullptr ||
              semantic_species.data == nullptr)))
          return Status{StatusCode::invalid_plan, kProductBinding};
        // Canonical composition is a physical EOS input.  Candidate scratch
        // FieldId/revision/storage differ on every replay and therefore belong
        // only to the exact certificate's rank-local scratch binding.
        std::uint64_t composition = detail::product_mix(
            kFnvOffset, UINT64_C(0x706563636f6d706f));
        composition = detail::product_mix(
            composition, product.thermodynamics.fingerprint());
        composition = detail::product_mix(
            composition,
            static_cast<std::uint64_t>(
                candidate_species.size));
        for (std::size_t species_index = 0U;
             species_index < candidate_species.size; ++species_index) {
          const ConstFieldView species = candidate_species.data[species_index];
          composition = mix_field_numeric(composition, species);
        }
        if (composition == 0U) composition = 1U;
        PisoExactEosClosureIdentity closure;
        closure.thermodynamics = pressure_reference_certificate.thermodynamics;
        closure.pressure_reference =
            pressure_reference_certificate.pressure_reference;
        closure.composition = composition;
        closure.pressure_state = make_piso_field_revision_identity(
            as_const(trial_pressure));
        closure.pressure_correction =
            make_piso_field_revision_identity(candidate_pressure_correction);
        closure.enthalpy_state = make_piso_field_revision_identity(
            as_const(trial_enthalpy));
        closure.enthalpy_correction =
            make_piso_field_revision_identity(candidate_enthalpy_correction);
        closure.candidate_enthalpy =
            make_piso_field_revision_identity(candidate_enthalpy);
        closure.candidate_density =
            make_piso_field_revision_identity(candidate_density);
        closure.candidate_temperature =
            make_piso_field_revision_identity(candidate_temperature);
        std::uint64_t closure_state = kFnvOffset;
        closure_state = detail::product_mix(
            closure_state, UINT64_C(0x7068656f73636c6f));
        closure_state = detail::product_mix(closure_state, pressure.state);
        closure_state = detail::product_mix(
            closure_state, pressure_reference_certificate.closure);
        closure_state = detail::product_mix(closure_state, composition);
        // This is the semantic EOS result, not its scratch allocation.  A
        // replay revises every candidate workspace, so revision tokens belong
        // in the field identities/scratch binding but must not perturb the
        // canonical closure carried by the closed-gauge result.
        closure_state = mix_field_numeric(
            closure_state, candidate_pressure_correction);
        closure_state = mix_field_numeric(
            closure_state, candidate_enthalpy_correction);
        closure_state = mix_field_numeric(closure_state,
                                          candidate_enthalpy);
        closure_state = mix_field_numeric(closure_state,
                                          candidate_density);
        closure_state = mix_field_numeric(closure_state,
                                          candidate_temperature);
        closure.closure = closure_state == 0U ? 1U : closure_state;
        candidate.enthalpy = candidate_enthalpy;
        candidate.density = candidate_density;
        candidate.temperature = candidate_temperature;
        candidate.closure = closure;
        candidate.independent_species = semantic_species;
        if (pressure_reference_kind == PressureReferenceKind::closed_mass) {
          candidate.pressure_compressibility =
              candidate_pressure_compressibility;
          ClosedGaugeCorrectionPrepareInput gauge;
          gauge.predecessor = pressure_reference_certificate;
          gauge.pressure_reference = attempt_pressure_reference;
          gauge.corrector = pressure.corrector;
          gauge.time = pressure.time;
          gauge.geometry = pressure.geometry;
          gauge.pressure_correction_authority = pressure.state;
          gauge.target_thermodynamic_closure = closure.closure;
          gauge.pressure_perturbation = as_const(trial_pressure);
          gauge.raw_pressure_correction = candidate_pressure_correction;
          gauge.candidate_pressure_compressibility =
              candidate.pressure_compressibility;
          gauge.activity = pressure_energy_activity;
          return product.equations.pressure_reference()
              .prepare_closed_gauge_correction(
                  gauge, product.reductions, candidate.closed_gauge);
        }
        return Status{};
      };
  const auto pressure_energy_direction_fingerprint =
      [&](std::uint8_t corrector,
          const PressureCorrectionCertificate& pressure,
          double (&local_maximum)[2U], Status& local) noexcept {
        local = {};
        local_maximum[0U] = 0.0;
        local_maximum[1U] = 0.0;
        std::uint64_t fingerprint = detail::product_mix(
            kFnvOffset, UINT64_C(0x7630347068646972));
        fingerprint = detail::product_mix(fingerprint, pressure.state);
        fingerprint = detail::product_mix(fingerprint, corrector);
        fingerprint = detail::product_mix(fingerprint, step.generation);
        const ConstFieldView pressure_direction =
            as_const(pressure_correction);
        const ConstFieldView enthalpy_direction =
            as_const(enthalpy_correction);
        const auto scalar_direction_shape = [&](ConstFieldView field) noexcept {
          return field.interior.x == cells.x &&
                 field.interior.y == cells.y &&
                 field.interior.z == cells.z && field.components == 1U;
        };
        if (!scalar_direction_shape(pressure_direction) ||
            !scalar_direction_shape(enthalpy_direction)) {
          local = {StatusCode::invalid_plan, kProductBinding};
          return fingerprint == 0U ? PlanFingerprint{1U} : fingerprint;
        }

        // Preserve the canonical fingerprint order: all dp bits precede all
        // dh bits.  The dh hash pass also consumes the corresponding dp value
        // for the existing paired finite/max certificate, removing the later
        // third traversal without changing its first-failure semantics.
        fingerprint = detail::product_mix(fingerprint,
                                          pressure_direction.components);
        for (std::int32_t z = 0; z < cells.z; ++z)
          for (std::int32_t y = 0; y < cells.y; ++y)
            for (std::int32_t x = 0; x < cells.x; ++x)
              fingerprint = detail::product_mix(
                  fingerprint,
                  product_double_bits(
                      pressure_direction.unchecked({x, y, z}, 0U)));
        fingerprint = detail::product_mix(fingerprint,
                                          enthalpy_direction.components);
        for (std::int32_t z = 0; z < cells.z; ++z)
          for (std::int32_t y = 0; y < cells.y; ++y)
            for (std::int32_t x = 0; x < cells.x; ++x) {
              const Int3 cell{x, y, z};
              const double dh = enthalpy_direction.unchecked(cell, 0U);
              fingerprint = detail::product_mix(
                  fingerprint, product_double_bits(dh));
              if (!local) continue;
              const double dp = pressure_direction.unchecked(cell, 0U);
              if (!std::isfinite(dp) || !std::isfinite(dh)) {
                local = {StatusCode::rejected_step,
                         kProductPressureEnergy};
                continue;
              }
              local_maximum[0U] =
                  std::max(local_maximum[0U], std::abs(dp));
              local_maximum[1U] =
                  std::max(local_maximum[1U], std::abs(dh));
            }
        return fingerprint == 0U ? PlanFingerprint{1U} : fingerprint;
      };
  const auto initialize_candidate_sample =
      [&](double alpha, std::size_t ordinal, std::uint8_t corrector,
          PlanFingerprint direction,
          PressureEnergyCandidateArtifacts& artifacts) noexcept {
        artifacts = {};
        PressureEnergyGlobalizationSample& sample = artifacts.sample;
        sample.alpha = alpha;
        sample.global_normalized_continuity =
            std::numeric_limits<double>::infinity();
        sample.global_normalized_energy =
            std::numeric_limits<double>::infinity();
        sample.corrector = corrector;
        sample.target_time = step.generation;
        sample.correction_direction = direction;
        std::uint64_t state = detail::product_mix(
            direction, UINT64_C(0x63616e6473746174));
        state = detail::product_mix(state, ordinal + 1U);
        state = detail::product_mix(state, product_double_bits(alpha));
        sample.state_provenance = state == 0U ? PlanFingerprint{1U} : state;
        std::uint64_t flux = detail::product_mix(
            direction, UINT64_C(0x63616e64666c7578));
        flux = detail::product_mix(flux, ordinal + 1U);
        flux = detail::product_mix(flux, product_double_bits(alpha));
        sample.mass_flux_provenance =
            flux == 0U ? PlanFingerprint{1U} : flux;
      };
  // Test diagnostics consume their one-shot arm before they replay the
  // candidate ladder.  Snapshot the request for the whole attempt so the
  // diagnostic cannot switch from the full oracle to the production fast
  // path halfway through that replay.
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const bool alpha_zero_diagnostic_oracle_requested =
      g_candidate_globalization_armed.load(std::memory_order_acquire) ||
      g_candidate_globalization_published.load(std::memory_order_acquire);
#endif
  const auto alpha_zero_diagnostic_oracle_armed = [&]() noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    return alpha_zero_diagnostic_oracle_requested;
#else
    return false;
#endif
  };
  const auto evaluate_pressure_energy_candidate =
      [&](std::uint8_t corrector,
          const PisoFrozenMomentumStageAuthority& frozen,
          const PressureCorrectionCertificate& pressure,
          const PisoFrozenMomentumExactCandidateCertificate* exact_baseline,
          PlanFingerprint direction, double alpha, std::size_t ordinal,
          ConstFaceFluxView baseline_flux,
          PressureEnergyCandidateArtifacts& artifacts) {
        initialize_candidate_sample(alpha, ordinal, corrector, direction,
                                    artifacts);
        const bool full_alpha_zero_oracle =
            alpha == 0.0 && alpha_zero_diagnostic_oracle_armed();
        Status evaluated = revise_pressure_energy_candidate_workspaces();
        if (evaluated &&
            pressure_energy_candidate_species.size() != species_trial.size())
          evaluated = {StatusCode::invalid_plan, kProductBinding};
        for (std::size_t species = 0U;
             species < species_trial.size() && evaluated; ++species)
          if ((evaluated = copy_interior(
                   as_const(species_trial[species]),
                   pressure_energy_candidate_species[species]))) {
            pressure_energy_candidate_species_const[species] =
                as_const(pressure_energy_candidate_species[species]);
            pressure_energy_candidate_species_history[species] =
                species_history[species];
            pressure_energy_candidate_species_history[species].trial =
                pressure_energy_candidate_species_const[species];
          }
        evaluated = product.reductions.consensus(evaluated);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 201U};

        const StageId correction_stage =
            corrector == 1U ? kCandidateCorrectionC1Stage
                            : kCandidateCorrectionC2Stage;
        const StageId state_stage =
            corrector == 1U ? kCandidateStateC1Stage
                            : kCandidateStateC2Stage;
        evaluated = product.coupler.form_frozen_momentum_scaled_pressure(
            frozen, as_const(pressure_correction),
            product.candidate_correction_halo, alpha,
            pressure_energy_candidate_pressure_correction,
            artifacts.pressure_stage);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 202U};
        halo_views[0U] = pressure_energy_candidate_pressure_correction;
        evaluated = exchange(product.candidate_correction_halo,
                             correction_stage, 1U, evaluated);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 203U};
        pressure_energy_candidate_pressure_correction = halo_views[0U];
        evaluated = product.coupler.stage_frozen_momentum_velocity(
            frozen, artifacts.pressure_stage,
            product.candidate_correction_halo,
            pressure_energy_candidate_pressure_correction,
            pressure_energy_candidate_velocity, artifacts.velocity_stage);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 204U};

        if (alpha == 0.0) {
          const bool velocity_matches = same_interior_bits(
              as_const(pressure_energy_candidate_velocity),
              as_const(trial_velocity));
          evaluated = product.reductions.consensus(
              velocity_matches
                  ? Status{}
                  : Status{StatusCode::invalid_plan,
                           kProductPressureEnergy});
          if (!evaluated) return evaluated;
        }

        const double alpha_zero_oracle_tolerance = std::max(
            product.summary.terminal_eos_tolerance,
            256.0 * std::numeric_limits<double>::epsilon());
        double local_alpha_zero_oracle[7U]{};
        std::uint64_t local_alpha_zero_oracle_lineage = detail::product_mix(
            kFnvOffset, UINT64_C(0x61306f7261636c65));
        Status local;
        for (std::int32_t z = 0; z < cells.z && local; ++z) {
          for (std::int32_t y = 0; y < cells.y && local; ++y) {
            for (std::int32_t x = 0; x < cells.x; ++x) {
              const Int3 cell{x, y, z};
              const double scaled_enthalpy =
                  alpha == 0.0
                      ? 0.0
                      : alpha * enthalpy_correction.unchecked(cell, 0U);
              const double candidate_pressure =
                  alpha == 0.0
                      ? trial_pressure.unchecked(cell, 0U)
                      : trial_pressure.unchecked(cell, 0U) +
                            pressure_energy_candidate_pressure_correction
                                .unchecked(cell, 0U);
              const double candidate_enthalpy =
                  alpha == 0.0
                      ? trial_enthalpy.unchecked(cell, 0U)
                      : trial_enthalpy.unchecked(cell, 0U) + scaled_enthalpy;
              pressure_energy_delta_temperature.unchecked(cell, 0U) =
                  scaled_enthalpy;
              if (!std::isfinite(candidate_pressure) ||
                  !std::isfinite(candidate_enthalpy) ||
                  !(attempt_pressure_reference + candidate_pressure > 0.0)) {
                local = {StatusCode::rejected_step,
                         kProductPressureEnergy};
                break;
              }
              pressure_energy_candidate_pressure.unchecked(cell, 0U) =
                  candidate_pressure;
              pressure_energy_candidate_enthalpy.unchecked(cell, 0U) =
                  candidate_enthalpy;
              if (alpha == 0.0) {
                // The live base is already an EOS-closed state.  In a closed
                // domain it was obtained through the density-authority inverse
                // map; applying the pressure-authority inverse map again can
                // move rho by one ulp even when p/h/T/U are bit identical.
                // Preserve that exact base closure so alpha zero and its
                // frozen total flux share the same density authority.
                const double base_density =
                    trial_density.unchecked(cell, 0U);
                const double base_temperature =
                    trial_temperature.unchecked(cell, 0U);
                if (full_alpha_zero_oracle) {
                  for (std::size_t species = 0U;
                       species < pressure_energy_candidate_species.size();
                       ++species)
                    species_values[species] =
                        pressure_energy_candidate_species[species].unchecked(
                            cell, 0U);
                  ThermoState oracle{};
                  Status oracle_status = product.thermodynamics
                                             .evaluate_from_reference_pressure(
                                                 attempt_pressure_reference,
                                                 candidate_pressure,
                                                 candidate_enthalpy,
                                                 {species_values.data(),
                                                  species_values.size()},
                                                 {pressure_energy_candidate_velocity
                                                      .unchecked(cell, 0U),
                                                  pressure_energy_candidate_velocity
                                                      .unchecked(cell, 1U),
                                                  pressure_energy_candidate_velocity
                                                      .unchecked(cell, 2U)},
                                                 oracle, base_temperature);
                  MolecularTransportState oracle_transport{};
                  if (oracle_status)
                    oracle_status = product.transport.evaluate(
                        oracle.temperature,
                        {species_values.data(), species_values.size()},
                        oracle_transport);
                  double oracle_conductivity = 0.0;
                  double oracle_enthalpy_diffusivity = 0.0;
                  if (oracle_status) {
                    oracle_conductivity = oracle_transport.conductivity;
                    oracle_enthalpy_diffusivity =
                        oracle_transport.conductivity / oracle.cp;
                    if (product.transport.kernel() ==
                        TransportKernel::coast_native_air)
                      oracle_status =
                          product.transport.effective_enthalpy_transport(
                              oracle_transport.viscosity,
                              effective_viscosity.unchecked(cell, 0U),
                              oracle.cp, oracle_conductivity,
                              oracle_enthalpy_diffusivity);
                  }
                  const auto normalized_difference =
                      [](double left, double right) noexcept {
                        return std::abs(left - right) /
                               std::max(
                                   {1.0, std::abs(left), std::abs(right)});
                      };
                  const double density_difference =
                      oracle_status
                          ? std::abs(oracle.rho - base_density)
                          : std::numeric_limits<double>::infinity();
                  const double temperature_difference =
                      oracle_status
                          ? std::abs(oracle.temperature - base_temperature)
                          : std::numeric_limits<double>::infinity();
                  double material_error =
                      std::numeric_limits<double>::infinity();
                  if (oracle_status) {
                    material_error = 0.0;
                    const std::array<std::pair<double, double>, 6U>
                        material_pairs{{
                            {oracle_transport.viscosity,
                             molecular_viscosity.unchecked(cell, 0U)},
                            {oracle.drho_dp_hY,
                             compressibility.unchecked(cell, 0U)},
                            {oracle.drho_dh_pY,
                             enthalpy_compressibility.unchecked(cell, 0U)},
                            {oracle_conductivity,
                             conductivity.unchecked(cell, 0U)},
                            {oracle.cp, heat_capacity.unchecked(cell, 0U)},
                            {oracle_enthalpy_diffusivity,
                             enthalpy_diffusivity.unchecked(cell, 0U)},
                        }};
                    for (const auto& pair : material_pairs)
                      material_error = std::max(
                          material_error,
                          normalized_difference(pair.first, pair.second));
                  }
                  const double density_error =
                      normalized_difference(oracle.rho, base_density);
                  const double temperature_error =
                      normalized_difference(oracle.temperature,
                                            base_temperature);
                  if (!oracle_status || !std::isfinite(density_difference) ||
                      !std::isfinite(temperature_difference) ||
                      !std::isfinite(material_error) ||
                      density_error > alpha_zero_oracle_tolerance ||
                      temperature_error > alpha_zero_oracle_tolerance ||
                      material_error > alpha_zero_oracle_tolerance) {
                    local = {StatusCode::rejected_step,
                             kProductPressureEnergy};
                    break;
                  }
                  local_alpha_zero_oracle[0U] = std::max(
                      local_alpha_zero_oracle[0U], density_difference);
                  local_alpha_zero_oracle[1U] = std::max(
                      local_alpha_zero_oracle[1U], temperature_difference);
                  local_alpha_zero_oracle[2U] = std::max(
                      local_alpha_zero_oracle[2U], material_error);
                  local_alpha_zero_oracle[3U] = std::max(
                      local_alpha_zero_oracle[3U], density_error);
                  local_alpha_zero_oracle[4U] = std::max(
                      local_alpha_zero_oracle[4U], temperature_error);
                  local_alpha_zero_oracle[5U] = std::max(
                      local_alpha_zero_oracle[5U],
                      static_cast<double>(
                          product_ulp_distance(oracle.rho, base_density)));
                  local_alpha_zero_oracle[6U] = std::max(
                      local_alpha_zero_oracle[6U],
                      static_cast<double>(product_ulp_distance(
                          oracle.temperature, base_temperature)));
                  for (double value :
                       {candidate_pressure, candidate_enthalpy, oracle.rho,
                        oracle.temperature, oracle.cp, oracle.drho_dp_hY,
                        oracle.drho_dh_pY, oracle_transport.viscosity,
                        oracle_conductivity, oracle_enthalpy_diffusivity})
                    local_alpha_zero_oracle_lineage = detail::product_mix(
                        local_alpha_zero_oracle_lineage,
                        product_double_bits(value));
                } else {
                  // The production baseline consumes the already-certified
                  // live closure.  Mix the actual bytes being reused so the
                  // fast path retains a deterministic, nonzero numeric oracle
                  // lineage without repeating EOS or transport evaluation.
                  for (double value :
                       {candidate_pressure, candidate_enthalpy, base_density,
                        base_temperature,
                        heat_capacity.unchecked(cell, 0U),
                        compressibility.unchecked(cell, 0U),
                        enthalpy_compressibility.unchecked(cell, 0U),
                        molecular_viscosity.unchecked(cell, 0U),
                        conductivity.unchecked(cell, 0U),
                        enthalpy_diffusivity.unchecked(cell, 0U)})
                    local_alpha_zero_oracle_lineage = detail::product_mix(
                        local_alpha_zero_oracle_lineage,
                        product_double_bits(value));
                }
                pressure_energy_candidate_density.unchecked(cell, 0U) =
                    base_density;
                pressure_energy_candidate_temperature.unchecked(cell, 0U) =
                    base_temperature;
                pressure_energy_candidate_molecular_viscosity.unchecked(
                    cell, 0U) = molecular_viscosity.unchecked(cell, 0U);
                pressure_energy_candidate_compressibility.unchecked(cell, 0U) =
                    compressibility.unchecked(cell, 0U);
                pressure_energy_candidate_enthalpy_compressibility.unchecked(
                    cell, 0U) =
                    enthalpy_compressibility.unchecked(cell, 0U);
                pressure_energy_candidate_thermal_conductivity.unchecked(
                    cell, 0U) = conductivity.unchecked(cell, 0U);
                pressure_energy_candidate_heat_capacity.unchecked(cell, 0U) =
                    heat_capacity.unchecked(cell, 0U);
                pressure_energy_candidate_enthalpy_diffusivity.unchecked(
                    cell, 0U) = enthalpy_diffusivity.unchecked(cell, 0U);
                continue;
              }
              for (std::size_t species = 0U;
                   species < pressure_energy_candidate_species.size();
                   ++species)
                species_values[species] =
                    pressure_energy_candidate_species[species].unchecked(
                        cell, 0U);
              ThermoState thermo;
              Status point =
                  product.thermodynamics.evaluate_from_reference_pressure(
                      attempt_pressure_reference, candidate_pressure,
                      candidate_enthalpy,
                      {species_values.data(), species_values.size()},
                      {pressure_energy_candidate_velocity.unchecked(cell,
                                                                      0U),
                       pressure_energy_candidate_velocity.unchecked(cell,
                                                                      1U),
                       pressure_energy_candidate_velocity.unchecked(cell,
                                                                      2U)},
                      thermo, trial_temperature.unchecked(cell, 0U));
              MolecularTransportState transport;
              if (point)
                point = product.transport.evaluate(
                    thermo.temperature,
                    {species_values.data(), species_values.size()},
                    transport);
              const bool admissible =
                  point && std::isfinite(thermo.rho) && thermo.rho > 0.0 &&
                  std::isfinite(thermo.temperature) &&
                  thermo.temperature > 0.0 && std::isfinite(thermo.cp) &&
                  thermo.cp > 0.0 &&
                  std::isfinite(thermo.drho_dp_hY) &&
                  thermo.drho_dp_hY > 0.0 &&
                  std::isfinite(thermo.drho_dh_pY) &&
                  thermo.drho_dh_pY < 0.0 &&
                  std::isfinite(transport.viscosity) &&
                  transport.viscosity > 0.0 &&
                  std::isfinite(transport.conductivity) &&
                  transport.conductivity > 0.0;
              if (!admissible) {
                local = {StatusCode::rejected_step,
                         kProductPressureEnergy};
                break;
              }
              pressure_energy_candidate_density.unchecked(cell, 0U) =
                  thermo.rho;
              pressure_energy_candidate_temperature.unchecked(cell, 0U) =
                  thermo.temperature;
              pressure_energy_candidate_molecular_viscosity.unchecked(
                  cell, 0U) = transport.viscosity;
              pressure_energy_candidate_compressibility.unchecked(cell,
                                                                    0U) =
                  thermo.drho_dp_hY;
              pressure_energy_candidate_enthalpy_compressibility.unchecked(
                  cell, 0U) = thermo.drho_dh_pY;
              pressure_energy_candidate_thermal_conductivity.unchecked(
                  cell, 0U) = transport.conductivity;
              pressure_energy_candidate_heat_capacity.unchecked(cell, 0U) =
                  thermo.cp;
              pressure_energy_candidate_enthalpy_diffusivity.unchecked(
                  cell, 0U) = transport.conductivity / thermo.cp;
            }
          }
        }
        evaluated = product.reductions.consensus(local);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 205U};
        if (full_alpha_zero_oracle) {
          double global_alpha_zero_oracle[7U]{};
          evaluated = product.reductions.checked_max(
              {local_alpha_zero_oracle, 7U},
              {global_alpha_zero_oracle, 7U});
          if (!evaluated) return evaluated;
          artifacts.alpha_zero_density_oracle_difference =
              global_alpha_zero_oracle[0U];
          artifacts.alpha_zero_temperature_oracle_difference =
              global_alpha_zero_oracle[1U];
          artifacts.alpha_zero_material_oracle_error =
              global_alpha_zero_oracle[2U];
          artifacts.alpha_zero_density_oracle_relative_error =
              global_alpha_zero_oracle[3U];
          artifacts.alpha_zero_temperature_oracle_relative_error =
              global_alpha_zero_oracle[4U];
          artifacts.alpha_zero_density_oracle_ulp =
              static_cast<std::uint64_t>(global_alpha_zero_oracle[5U]);
          artifacts.alpha_zero_temperature_oracle_ulp =
              static_cast<std::uint64_t>(global_alpha_zero_oracle[6U]);
        }

        ConstFieldView candidate_compressibility_semantic = semantic_field(
            as_const(pressure_energy_candidate_compressibility),
            product.fields.compressibility);
        evaluated = prepare_exact_thermodynamic_candidate(
            pressure,
            as_const(pressure_energy_candidate_pressure_correction),
            as_const(pressure_energy_delta_temperature),
            as_const(pressure_energy_candidate_enthalpy),
            as_const(pressure_energy_candidate_density),
            as_const(pressure_energy_candidate_temperature),
            candidate_compressibility_semantic,
            {pressure_energy_candidate_species_const.data(),
             pressure_energy_candidate_species_const.size()},
            {species_accepted.data(), species_accepted.size()},
            artifacts.exact);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 206U};
        artifacts.pressure_reference = attempt_pressure_reference;
        if (pressure_reference_kind == PressureReferenceKind::closed_mass) {
          if (!artifacts.exact.closed_gauge.valid())
            return Status{StatusCode::invalid_plan,
                          kProductPressureEnergy};
          artifacts.pressure_reference =
              artifacts.exact.closed_gauge.next_pressure_reference;
          const double shift = artifacts.exact.closed_gauge.shift;
          for (std::int32_t z = 0; z < cells.z; ++z)
            for (std::int32_t y = 0; y < cells.y; ++y)
              for (std::int32_t x = 0; x < cells.x; ++x)
                pressure_energy_candidate_pressure.unchecked({x, y, z},
                                                               0U) -= shift;
        }

        // Establish every candidate neighbour authority before any physical
        // boundary closure mutates face ghosts.  IBM reconstruction is a
        // separate remote-donor authority below; neither route is permitted
        // to borrow the live trial FieldIds.
        halo_views[0U] = pressure_energy_candidate_velocity;
        halo_views[1U] = pressure_energy_candidate_density;
        halo_views[2U] = pressure_energy_candidate_pressure;
        halo_views[3U] = pressure_energy_candidate_enthalpy;
        halo_views[4U] = pressure_energy_candidate_temperature;
        for (std::size_t species = 0U;
             species < pressure_energy_candidate_species.size(); ++species)
          halo_views[species + 5U] =
              pressure_energy_candidate_species[species];
        const std::size_t candidate_state_halo_count =
            5U + pressure_energy_candidate_species.size();
        evaluated =
            exchange(product.candidate_state_halo, state_stage,
                     candidate_state_halo_count, evaluated);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 207U};
        pressure_energy_candidate_velocity = halo_views[0U];
        pressure_energy_candidate_density = halo_views[1U];
        pressure_energy_candidate_pressure = halo_views[2U];
        pressure_energy_candidate_enthalpy = halo_views[3U];
        pressure_energy_candidate_temperature = halo_views[4U];
        for (std::size_t species = 0U;
             species < pressure_energy_candidate_species.size(); ++species) {
          pressure_energy_candidate_species[species] =
              halo_views[species + 5U];
          pressure_energy_candidate_species_const[species] =
              as_const(pressure_energy_candidate_species[species]);
          pressure_energy_candidate_species_history[species].trial =
              pressure_energy_candidate_species_const[species];
        }

        if (pressure_reference_kind ==
            PressureReferenceKind::boundary_absolute) {
          halo_count = 0U;
          halo_views[halo_count++] = pressure_energy_candidate_pressure;
          halo_views[halo_count++] = pressure_energy_candidate_enthalpy;
          halo_views[halo_count++] = pressure_energy_candidate_density;
          halo_views[halo_count++] = pressure_energy_candidate_temperature;
          halo_views[halo_count++] = pressure_energy_candidate_velocity;
          for (FieldView species : pressure_energy_candidate_species)
            halo_views[halo_count++] = species;
          evaluated = exchange(product.candidate_finalizer_state_halo,
                               corrector == 1U
                                   ? kCandidateFinalizerStateC1Stage
                                   : kCandidateFinalizerStateC2Stage,
                               halo_count, evaluated);
          if (!evaluated)
            return Status{evaluated.code,
                          kProductPressureEnergy + 252U};
          pressure_energy_candidate_pressure = halo_views[0U];
          pressure_energy_candidate_enthalpy = halo_views[1U];
          pressure_energy_candidate_density = halo_views[2U];
          pressure_energy_candidate_temperature = halo_views[3U];
          pressure_energy_candidate_velocity = halo_views[4U];
          for (std::size_t species = 0U;
               species < pressure_energy_candidate_species.size(); ++species) {
            pressure_energy_candidate_species[species] =
                halo_views[species + 5U];
            pressure_energy_candidate_species_const[species] =
                as_const(pressure_energy_candidate_species[species]);
            pressure_energy_candidate_species_history[species].trial =
                pressure_energy_candidate_species_const[species];
            FieldView semantic_species =
                pressure_energy_candidate_species[species];
            semantic_species.field = species_trial[species].field;
            pressure_energy_candidate_species_boundary_aliases[species] =
                as_const(semantic_species);
          }
          evaluated = resolve_static_boundary_values(
              communicator, product.boundary, product.boundary_specs,
              product.geometry, product.patch, product.thermodynamics,
              artifacts.pressure_reference,
              as_const(pressure_energy_candidate_density),
              as_const(pressure_energy_candidate_velocity),
              as_const(pressure_energy_candidate_enthalpy),
              {pressure_energy_candidate_species_boundary_aliases.data(),
               pressure_energy_candidate_species_boundary_aliases.size()},
              {passive_accepted.data(), passive_accepted.size()},
              {species_values.data(), species_values.size()},
              boundary_scalars, boundary_vectors,
              boundary_normal_gradients, true, evaluated);
          FieldView semantic_pressure = pressure_energy_candidate_pressure;
          semantic_pressure.field = product.fields.pressure;
          FieldView semantic_enthalpy = pressure_energy_candidate_enthalpy;
          semantic_enthalpy.field = product.fields.enthalpy;
          FieldView semantic_velocity = pressure_energy_candidate_velocity;
          semantic_velocity.field = product.fields.velocity;
          if (evaluated)
            evaluated = apply_boundary_ghosts(
                BoundaryStage::pressure, product.boundary,
                {&semantic_pressure, 1U}, boundary_values);
          if (evaluated)
            evaluated = apply_boundary_ghosts(
                BoundaryStage::enthalpy, product.boundary,
                {&semantic_enthalpy, 1U}, boundary_values);
          if (evaluated)
            evaluated = apply_boundary_ghosts(
                BoundaryStage::momentum, product.boundary,
                {&semantic_velocity, 1U}, boundary_values);
          if (evaluated && !product.fields.scalars.empty()) {
            std::size_t species_index = 0U;
            std::size_t passive_index = 0U;
            for (std::size_t scalar = 0U;
                 scalar < product.fields.scalars.size(); ++scalar) {
              if (product.fields.scalar_roles[scalar] ==
                  TransportedScalarRole::species) {
                if (species_index >=
                    pressure_energy_candidate_species.size()) {
                  evaluated = {StatusCode::invalid_plan, kProductBinding};
                  break;
                }
                halo_views[scalar] =
                    pressure_energy_candidate_species[species_index];
                halo_views[scalar].field =
                    product.fields.scalars[scalar];
                ++species_index;
              } else {
                if (passive_index >= passive_low.size()) {
                  evaluated = {StatusCode::invalid_plan, kProductBinding};
                  break;
                }
                // BoundaryStage::scalar validates the complete transported
                // field set.  Candidate energy/EOS does not consume passive
                // scalars, so use the already-preallocated low-order scratch
                // solely as a non-live carrier for those physical ghosts.
                halo_views[scalar] = passive_low[passive_index++];
              }
            }
            if (species_index != pressure_energy_candidate_species.size() ||
                passive_index != passive_trial.size())
              evaluated = {StatusCode::invalid_plan, kProductBinding};
            if (evaluated)
              evaluated = apply_boundary_ghosts(
                  BoundaryStage::scalar, product.boundary,
                  {halo_views.data(), product.fields.scalars.size()},
                  boundary_values);
            species_index = 0U;
            passive_index = 0U;
            for (std::size_t scalar = 0U;
                 scalar < product.fields.scalars.size() && evaluated;
                 ++scalar) {
              if (product.fields.scalar_roles[scalar] ==
                  TransportedScalarRole::species) {
                pressure_energy_candidate_species_boundary_aliases
                    [species_index++] = as_const(halo_views[scalar]);
              } else {
                ++passive_index;
              }
            }
          }
          if (evaluated) {
            std::vector<BoundaryGhostFieldAuthority>& authorities =
                boundary_thermo_authority_fields;
            if (authorities.size() !=
                pressure_energy_candidate_species.size() + 2U) {
              evaluated = {StatusCode::invalid_plan, kProductBinding};
            } else {
              authorities[0U] = make_boundary_ghost_field_authority(
                  as_const(semantic_pressure));
              authorities[1U] = make_boundary_ghost_field_authority(
                  as_const(semantic_enthalpy));
              for (std::size_t species = 0U;
                   species <
                   pressure_energy_candidate_species_boundary_aliases.size();
                   ++species)
                authorities[species + 2U] =
                    make_boundary_ghost_field_authority(
                        pressure_energy_candidate_species_boundary_aliases
                            [species]);
              std::uint64_t producer = detail::product_mix(
                  product.candidate_finalizer_state_halo.instance_identity(),
                  state_stage);
              producer = detail::product_mix(
                  producer, artifacts.pressure_stage.canonical_lineage());
              if (producer == 0U) producer = 1U;
              const std::int32_t reach_value = std::max(
                  {semantic_pressure.ghosts.x, semantic_pressure.ghosts.y,
                   semantic_pressure.ghosts.z});
              if (reach_value <= 0 || reach_value > UINT8_MAX) {
                evaluated = {StatusCode::invalid_plan, kProductBinding};
              } else {
                const BoundaryThermophysicalGhostAuthority authority{
                    static_cast<std::uintptr_t>(producer),
                    product.boundary.revision(),
                    product.boundary.local_layout_fingerprint(), cells,
                    semantic_pressure.ghosts,
                    static_cast<std::uint8_t>(reach_value),
                    {authorities.data(), authorities.size()}};
                evaluated = BoundaryThermophysicalFaceClosure::close(
                    product.boundary, product.thermodynamics,
                    product.transport,
                    {artifacts.pressure_reference,
                     as_const(semantic_pressure), as_const(semantic_enthalpy),
                     {pressure_energy_candidate_species_boundary_aliases
                          .data(),
                      pressure_energy_candidate_species_boundary_aliases
                          .size()},
                     authority},
                    {pressure_energy_candidate_density,
                     pressure_energy_candidate_temperature,
                     pressure_energy_candidate_heat_capacity,
                     pressure_energy_candidate_compressibility,
                     pressure_energy_candidate_enthalpy_compressibility,
                     pressure_energy_candidate_molecular_viscosity,
                     pressure_energy_candidate_thermal_conductivity,
                     pressure_energy_candidate_enthalpy_diffusivity},
                    {step.generation, product.geometry.fingerprint(),
                     pressure_reference_certificate.pressure_reference,
                     product.boundary.revision(),
                     corrector == 1U
                         ? BoundaryThermophysicalGhostPhase::corrector_one
                         : BoundaryThermophysicalGhostPhase::corrector_two},
                    artifacts.boundary_thermophysics);
              }
            }
          }
          evaluated = product.reductions.consensus(evaluated);
          if (!evaluated)
            return Status{evaluated.code,
                          kProductPressureEnergy + 251U};
        }

        if (evaluated && product.ibm_candidate_velocity_donors.has_value()) {
          std::array<FieldView, 1U> donor_fields{
              pressure_energy_candidate_velocity};
          evaluated = product.ibm_candidate_velocity_donors->exchange(
              kIbmCandidateVelocityDonorStage,
              {donor_fields.data(), donor_fields.size()});
          pressure_energy_candidate_velocity = donor_fields[0U];
        }
        evaluated = product.reductions.consensus(evaluated);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 255U};
        if (evaluated && product.ibm_candidate_rate_donors.has_value()) {
          std::array<FieldView, 2U> donor_fields{
              pressure_energy_candidate_pressure,
              pressure_energy_candidate_temperature};
          evaluated = product.ibm_candidate_rate_donors->exchange(
              kIbmCandidateEnergyRateDonorStage,
              {donor_fields.data(), donor_fields.size()});
          pressure_energy_candidate_pressure = donor_fields[0U];
          pressure_energy_candidate_temperature = donor_fields[1U];
        }
        evaluated = product.reductions.consensus(evaluated);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 255U};

        if (alpha == 0.0 && !full_alpha_zero_oracle) {
          // Alpha zero has exactly the live primitive state.  Reuse the
          // certified live gradient directly; boundary/IBM and final-flux
          // replay below remain candidate-local and are not bypassed.
          for (std::int32_t z = 0; z < cells.z; ++z)
            for (std::int32_t y = 0; y < cells.y; ++y)
              for (std::int32_t x = 0; x < cells.x; ++x)
                for (std::uint8_t component = 0U; component < 9U;
                     ++component) {
                  const double value =
                      velocity_gradient.unchecked({x, y, z}, component);
                  pressure_energy_candidate_velocity_gradient.unchecked(
                      {x, y, z}, component) = value;
                  local_alpha_zero_oracle_lineage = detail::product_mix(
                      local_alpha_zero_oracle_lineage,
                      product_double_bits(value));
                }
        } else {
          {
            const std::array<ConstFieldView, 1U> reads{
                as_const(pressure_energy_candidate_velocity)};
            const std::array<FieldView, 1U> writes{
                pressure_energy_candidate_velocity_gradient};
            local = cartesian_gradient(
                product.equations.kernels(),
                {{reads.data(), reads.size()}, {writes.data(), writes.size()},
                 full_box, 0U, 0U, 3U, 0U, nullptr});
          }
          if (local && product.ibm_equations.has_value())
            local = product.ibm_equations->correct_velocity_gradient(
                as_const(pressure_energy_candidate_velocity),
                pressure_energy_candidate_velocity_gradient);
          evaluated = product.reductions.consensus(local);
          if (!evaluated)
            return Status{evaluated.code,
                          kProductPressureEnergy + 208U};
        }
        if (full_alpha_zero_oracle) {
          double local_gradient_error = 0.0;
          for (std::int32_t z = 0; z < cells.z; ++z)
            for (std::int32_t y = 0; y < cells.y; ++y)
              for (std::int32_t x = 0; x < cells.x; ++x)
                for (std::uint8_t component = 0U; component < 9U;
                     ++component) {
                  const double oracle =
                      pressure_energy_candidate_velocity_gradient.unchecked(
                          {x, y, z}, component);
                  const double base = velocity_gradient.unchecked(
                      {x, y, z}, component);
                  const double scale =
                      std::max({1.0, std::abs(oracle), std::abs(base)});
                  local_gradient_error = std::max(
                      local_gradient_error, std::abs(oracle - base) / scale);
                  local_alpha_zero_oracle_lineage = detail::product_mix(
                      local_alpha_zero_oracle_lineage,
                      product_double_bits(oracle));
                }
          double global_gradient_error = 0.0;
          evaluated = product.reductions.checked_max(
              {&local_gradient_error, 1U}, {&global_gradient_error, 1U},
              std::isfinite(local_gradient_error) &&
                      local_gradient_error <= alpha_zero_oracle_tolerance
                  ? Status{}
                  : Status{StatusCode::rejected_step,
                           kProductPressureEnergy});
          if (!evaluated) return evaluated;
          artifacts.alpha_zero_gradient_oracle_error =
              global_gradient_error;
          evaluated = copy_interior(
              as_const(velocity_gradient),
              pressure_energy_candidate_velocity_gradient);
          if (!evaluated) return evaluated;
        }
        if (alpha == 0.0 && !full_alpha_zero_oracle) {
          for (std::int32_t z = 0; z < cells.z; ++z)
            for (std::int32_t y = 0; y < cells.y; ++y)
              for (std::int32_t x = 0; x < cells.x; ++x) {
                const double value =
                    effective_viscosity.unchecked({x, y, z}, 0U);
                pressure_energy_candidate_effective_viscosity.unchecked(
                    {x, y, z}, 0U) = value;
                local_alpha_zero_oracle_lineage = detail::product_mix(
                    local_alpha_zero_oracle_lineage,
                    product_double_bits(value));
              }
        } else {
          const TurbulenceCandidateInput turbulence_input{
              as_const(pressure_energy_candidate_density),
              as_const(pressure_energy_candidate_molecular_viscosity),
              as_const(pressure_energy_candidate_velocity_gradient),
              pressure_energy_candidate_velocity_gradient.revision};
          // This certificate only authenticates the independently recomputed
          // oracle output.  Alpha zero subsequently restores the conservative
          // live material bytes, so retaining the certificate in the published
          // artifacts would falsely authenticate an overwritten view.
          TurbulenceCandidateCertificate turbulence_oracle;
          evaluated =
              product.turbulence.evaluate_candidate_effective_viscosity(
                  turbulence_input,
                  pressure_energy_candidate_effective_viscosity,
                  turbulence_oracle);
          evaluated = product.reductions.consensus(evaluated);
          if (!evaluated)
            return Status{evaluated.code,
                          kProductPressureEnergy + 209U};
        }
        if (full_alpha_zero_oracle) {
          double local_effective_error = 0.0;
          for (std::int32_t z = 0; z < cells.z; ++z)
            for (std::int32_t y = 0; y < cells.y; ++y)
              for (std::int32_t x = 0; x < cells.x; ++x) {
                const double oracle =
                    pressure_energy_candidate_effective_viscosity.unchecked(
                        {x, y, z}, 0U);
                const double base =
                    effective_viscosity.unchecked({x, y, z}, 0U);
                const double scale =
                    std::max({1.0, std::abs(oracle), std::abs(base)});
                local_effective_error = std::max(
                    local_effective_error, std::abs(oracle - base) / scale);
                local_alpha_zero_oracle_lineage = detail::product_mix(
                    local_alpha_zero_oracle_lineage,
                    product_double_bits(oracle));
              }
          double global_effective_error = 0.0;
          evaluated = product.reductions.checked_max(
              {&local_effective_error, 1U}, {&global_effective_error, 1U},
              std::isfinite(local_effective_error) &&
                      local_effective_error <= alpha_zero_oracle_tolerance
                  ? Status{}
                  : Status{StatusCode::rejected_step,
                           kProductPressureEnergy});
          if (!evaluated) return evaluated;
          artifacts.alpha_zero_effective_viscosity_oracle_error =
              global_effective_error;
          evaluated = copy_interior(
              as_const(effective_viscosity),
              pressure_energy_candidate_effective_viscosity);
          if (!evaluated) return evaluated;
        }
        if (evaluated)
          evaluated = refresh_coast_native_air_effective_thermal_transport(
              product.transport,
              as_const(pressure_energy_candidate_molecular_viscosity),
              as_const(pressure_energy_candidate_effective_viscosity),
              as_const(pressure_energy_candidate_heat_capacity),
              pressure_energy_candidate_thermal_conductivity,
              pressure_energy_candidate_enthalpy_diffusivity);
        if (product.transport.kernel() == TransportKernel::coast_native_air) {
          evaluated = exchange_effective_thermal_ghosts(
              product.candidate_thermal_halo, state_stage, product.boundary,
              pressure_energy_candidate_thermal_conductivity,
              pressure_energy_candidate_enthalpy_diffusivity, evaluated);
          evaluated = product.reductions.consensus(evaluated);
        }
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 209U};

        RevisionToken flux_revision = detail::product_mix(
            direction, product_double_bits(alpha));
        flux_revision = detail::product_mix(flux_revision, ordinal + 1U);
        if (flux_revision == 0U) flux_revision = 1U;
        FaceFluxView candidate_flux;
        evaluated = product.phi_workspace.workspace_view(
            2U, flux_revision, candidate_flux);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 210U};
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        if (alpha == 0.0 && g_candidate_globalization_published.load(
                                std::memory_order_acquire)) {
          double local_difference[6U]{};
          for (std::int32_t z = 0; z < cells.z; ++z)
            for (std::int32_t y = 0; y < cells.y; ++y)
              for (std::int32_t x = 0; x < cells.x; ++x) {
                const Int3 cell{x, y, z};
                for (std::uint8_t component = 0U; component < 3U;
                     ++component) {
                  const double candidate_value =
                      pressure_energy_candidate_velocity.unchecked(
                          cell, component);
                  const double base_value =
                      trial_velocity.unchecked(cell, component);
                  local_difference[0U] = std::max(
                      local_difference[0U],
                      std::abs(candidate_value - base_value));
                  if (product_double_bits(candidate_value) !=
                      product_double_bits(base_value))
                    local_difference[3U] = 1.0;
                }
                const double candidate_density =
                    pressure_energy_candidate_density.unchecked(cell, 0U);
                const double base_density =
                    trial_density.unchecked(cell, 0U);
                local_difference[1U] = std::max(
                    local_difference[1U],
                    std::abs(candidate_density - base_density));
                if (product_double_bits(candidate_density) !=
                    product_double_bits(base_density))
                  local_difference[4U] = 1.0;
                const double candidate_temperature =
                    pressure_energy_candidate_temperature.unchecked(cell,
                                                                      0U);
                const double base_temperature =
                    trial_temperature.unchecked(cell, 0U);
                local_difference[2U] = std::max(
                    local_difference[2U],
                    std::abs(candidate_temperature - base_temperature));
                if (product_double_bits(candidate_temperature) !=
                    product_double_bits(base_temperature))
                  local_difference[5U] = 1.0;
              }
          double global_difference[6U]{};
          const Status diagnostic_status = product.reductions.checked_max(
              {local_difference, 6U}, {global_difference, 6U});
          if (diagnostic_status) {
            detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic =
                g_candidate_globalization_diagnostic;
            diagnostic.alpha_zero_maximum_velocity_difference =
                global_difference[0U];
            diagnostic.alpha_zero_maximum_density_difference =
                artifacts.alpha_zero_density_oracle_difference;
            diagnostic.alpha_zero_maximum_temperature_difference =
                artifacts.alpha_zero_temperature_oracle_difference;
            diagnostic.alpha_zero_maximum_density_relative_error =
                artifacts.alpha_zero_density_oracle_relative_error;
            diagnostic.alpha_zero_maximum_temperature_relative_error =
                artifacts.alpha_zero_temperature_oracle_relative_error;
            diagnostic.alpha_zero_maximum_density_ulp_difference =
                artifacts.alpha_zero_density_oracle_ulp;
            diagnostic.alpha_zero_maximum_temperature_ulp_difference =
                artifacts.alpha_zero_temperature_oracle_ulp;
            diagnostic.alpha_zero_material_oracle_error =
                artifacts.alpha_zero_material_oracle_error;
            diagnostic.alpha_zero_gradient_oracle_error =
                artifacts.alpha_zero_gradient_oracle_error;
            diagnostic.alpha_zero_energy_residual_oracle_error =
                artifacts.alpha_zero_energy_residual_oracle_error;
            diagnostic.alpha_zero_effective_viscosity_oracle_error =
                artifacts.alpha_zero_effective_viscosity_oracle_error;
            diagnostic.alpha_zero_oracle_numeric_lineage =
                artifacts.alpha_zero_oracle_numeric_lineage;
            diagnostic.alpha_zero_velocity_bitwise_equal =
                global_difference[3U] == 0.0;
            diagnostic.alpha_zero_density_bitwise_equal =
                global_difference[4U] == 0.0;
            diagnostic.alpha_zero_temperature_bitwise_equal =
                global_difference[5U] == 0.0;
            publish_candidate_globalization_diagnostic(diagnostic);
          }
        }
#endif
        // stage_frozen_momentum_flux is the single alpha-zero density
        // authority.  It checks the complete face-neighbour envelope against
        // the density snapshot captured by refresh() and publishes any
        // mismatch collectively before forming a face flux.  Do not repeat
        // the same full-field scan and consensus here.
        evaluated = product.coupler.stage_frozen_momentum_flux(
            frozen, artifacts.velocity_stage,
            as_const(pressure_energy_candidate_density),
            candidate_flux, artifacts.flux_stage);
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 211U};
        artifacts.flux = as_const(candidate_flux);
        if (pressure_reference_kind ==
            PressureReferenceKind::boundary_absolute) {
          RevisionToken final_flux_revision = detail::product_mix(
              flux_revision, UINT64_C(0x66696e616c706879));
          if (final_flux_revision == 0U) final_flux_revision = 1U;
          FaceFluxView candidate_final_flux;
          evaluated = product.phi_workspace.workspace_view(
              3U, final_flux_revision, candidate_final_flux);
          if (!evaluated)
            return Status{evaluated.code,
                          kProductPressureEnergy + 253U};
          ConstFieldView semantic_pressure =
              as_const(pressure_energy_candidate_pressure);
          semantic_pressure.field = product.fields.pressure;
          ConstFieldView semantic_enthalpy =
              as_const(pressure_energy_candidate_enthalpy);
          semantic_enthalpy.field = product.fields.enthalpy;
          const PressureEnergyCandidateBoundaryFinalizeInput finalizer_input{
              frozen,
              artifacts.pressure_stage,
              artifacts.velocity_stage,
              artifacts.flux_stage,
              pressure_reference_certificate,
              artifacts.pressure_reference,
              as_const(pressure_energy_candidate_pressure),
              as_const(pressure_energy_candidate_enthalpy),
              as_const(pressure_energy_candidate_density),
              as_const(pressure_energy_candidate_temperature),
              as_const(pressure_energy_candidate_velocity),
              {pressure_energy_candidate_species_const.data(),
               pressure_energy_candidate_species_const.size()},
              artifacts.exact.closure.composition,
              {artifacts.boundary_thermophysics,
               {artifacts.pressure_reference, semantic_pressure,
                semantic_enthalpy,
                {pressure_energy_candidate_species_boundary_aliases.data(),
                 pressure_energy_candidate_species_boundary_aliases.size()},
                as_const(pressure_energy_candidate_density)}},
              &product.candidate_finalizer_state_halo,
              as_const(candidate_flux),
              candidate_final_flux};
          evaluated = product.candidate_boundary_finalizer.finalize(
              finalizer_input, product.reductions,
              artifacts.final_boundary);
          if (!evaluated)
            return Status{evaluated.code,
                          kProductPressureEnergy + 254U};
          artifacts.flux = as_const(candidate_final_flux);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          // One-shot, rank-local mutation after the independent finalizer has
          // sealed the state/flux pair.  Exact certification must detect this
          // no-revision write collectively and leave the live transaction
          // replayable; production builds contain no mutation seam.
          const std::uint8_t poison_corrector =
              g_candidate_poison_corrector.load(std::memory_order_acquire);
          if ((poison_corrector == 0U || poison_corrector == corrector) &&
              g_candidate_poison_armed.exchange(
                  false, std::memory_order_acq_rel)) {
            int local_rank = -1;
            if (MPI_Comm_rank(communicator, &local_rank) == MPI_SUCCESS &&
                local_rank == g_candidate_poison_rank.load(
                                  std::memory_order_acquire) &&
                static_cast<detail::PressureEnergyCandidatePoisonKind>(
                    g_candidate_poison_kind.load(
                        std::memory_order_acquire)) ==
                    detail::PressureEnergyCandidatePoisonKind::
                        density_without_revision) {
              double& poisoned = pressure_energy_candidate_density.unchecked(
                  {0, 0, 0}, 0U);
              poisoned += 0.125;
            }
          }
#endif
        }
        const PisoFrozenMomentumExactCandidateInput exact_input{
            as_const(enthalpy_correction),
            as_const(pressure_energy_candidate_pressure_correction),
            as_const(pressure_energy_delta_temperature),
            {trial_velocity, trial_pressure, trial_enthalpy, trial_density,
             trial_temperature},
            as_const(pressure_energy_candidate_pressure),
            artifacts.exact,
            as_const(pressure_energy_candidate_velocity),
            artifacts.flux,
            artifacts.final_boundary,
            {pressure_energy_candidate_species_const.data(),
             pressure_energy_candidate_species_const.size()}};
        if (alpha == 0.0) {
          if (exact_baseline != nullptr)
            return Status{StatusCode::invalid_plan,
                          kProductPressureEnergy};
          evaluated = product.coupler.certify_frozen_momentum_exact_baseline(
              frozen, artifacts.pressure_stage, artifacts.velocity_stage,
              artifacts.flux_stage, exact_input, product.reductions,
              artifacts.exact_certificate);
        } else {
          if (exact_baseline == nullptr || !exact_baseline->valid())
            return Status{StatusCode::invalid_plan,
                          kProductPressureEnergy};
          evaluated = product.coupler.certify_frozen_momentum_exact_candidate(
              frozen, *exact_baseline, artifacts.pressure_stage,
              artifacts.velocity_stage, artifacts.flux_stage, exact_input,
              product.reductions, artifacts.exact_certificate);
        }
        if (!evaluated)
          return Status{evaluated.code,
                        kProductPressureEnergy + 212U};
        if (!artifacts.exact_certificate.valid() ||
            artifacts.exact_certificate.alpha() != alpha ||
            artifacts.exact_certificate.corrector() != corrector) {
          return Status{StatusCode::invalid_plan,
                        kProductPressureEnergy};
        }
        PressureEnergyGlobalizationSample& certified_sample =
            artifacts.sample;
        certified_sample.target_time =
            artifacts.exact_certificate.target_time();
        certified_sample.correction_direction =
            artifacts.exact_certificate.correction_direction();
        certified_sample.state_provenance =
            artifacts.exact_certificate.candidate_state_provenance();
        certified_sample.mass_flux_provenance =
            artifacts.exact_certificate.candidate_mass_flux_provenance();

        EquationStateView candidate_state = equation_state;
        ConstFieldView candidate_density_semantic =
            semantic_field(as_const(pressure_energy_candidate_density),
                           product.fields.rho);
        candidate_state.density.trial = candidate_density_semantic;
        candidate_state.velocity.trial = semantic_field(
            as_const(pressure_energy_candidate_velocity),
            product.fields.velocity);
        candidate_state.pressure_perturbation.trial = semantic_field(
            as_const(pressure_energy_candidate_pressure),
            product.fields.pressure);
        candidate_state.enthalpy.trial = semantic_field(
            as_const(pressure_energy_candidate_enthalpy),
            product.fields.enthalpy);
        candidate_state.temperature.trial = semantic_field(
            as_const(pressure_energy_candidate_temperature),
            product.fields.temperature);
        candidate_state.independent_species = {
            pressure_energy_candidate_species_history.data(),
            pressure_energy_candidate_species_history.size()};
        candidate_state.pressure_reference = artifacts.pressure_reference;
        EquationMaterialView candidate_material = material;
        candidate_material.molecular_viscosity = semantic_field(
            as_const(pressure_energy_candidate_molecular_viscosity),
            product.fields.molecular_viscosity);
        candidate_material.effective_viscosity = semantic_field(
            as_const(pressure_energy_candidate_effective_viscosity),
            product.fields.effective_viscosity);
        candidate_material.thermal_conductivity = semantic_field(
            as_const(pressure_energy_candidate_thermal_conductivity),
            product.fields.thermal_conductivity);
        candidate_material.heat_capacity = semantic_field(
            as_const(pressure_energy_candidate_heat_capacity),
            product.fields.heat_capacity);
        candidate_material.enthalpy_diffusivity = semantic_field(
            as_const(pressure_energy_candidate_enthalpy_diffusivity),
            product.fields.enthalpy_diffusivity);
        candidate_material.pressure_compressibility = semantic_field(
            as_const(pressure_energy_candidate_compressibility),
            product.fields.compressibility);
        ConstFieldView candidate_gradient = semantic_field(
            as_const(pressure_energy_candidate_velocity_gradient),
            product.fields.velocity_gradient);
        EquationAssemblyContext candidate_assembly = assembly;
        candidate_assembly.scope = EquationAssemblyScope::target_coupled;
        candidate_assembly.face_flux = artifacts.flux.revision;
        candidate_assembly.face_flux_authority = 0U;
        candidate_assembly.face_flux_storage = 0U;
        candidate_assembly.face_flux_revision_domain = 0U;
        candidate_assembly.mass_flux = artifacts.flux;
        candidate_assembly.provisional_mass_flux = false;
        candidate_assembly.immersed_interface =
            product.ibm_equations.has_value() ? &*product.ibm_equations
                                              : nullptr;
        candidate_assembly.wall_treatment = nullptr;
        const auto assemble_candidate_energy = [&]() noexcept {
          const TargetCoupledEnthalpyResidualWorkspace residual_workspace{
              pressure_energy_e_h, pressure_energy_r_c,
              pressure_energy_e_p};
          Status assembled = assemble_target_coupled_enthalpy_residual(
              product.equations.enthalpy(), candidate_state,
              candidate_material, candidate_gradient, candidate_assembly,
              pressure_energy_r_e, residual_workspace, artifacts.energy);
          if (assembled && product.ibm_equations.has_value()) {
            zero_field(pressure_energy_e_p);
            assembled = product.ibm_equations->correct_pressure_work(
                as_const(pressure_energy_candidate_pressure),
                as_const(pressure_energy_candidate_velocity),
                pressure_energy_e_p);
            for (std::int32_t z = 0; z < cells.z && assembled; ++z)
              for (std::int32_t y = 0; y < cells.y && assembled; ++y)
                for (std::int32_t x = 0; x < cells.x; ++x) {
                  const Int3 cell{x, y, z};
                  pressure_energy_r_e.unchecked(cell, 0U) -=
                      detail::cell_volume(product.equations.kernels(), cell) *
                      pressure_energy_e_p.unchecked(cell, 0U);
                }
            if (assembled) {
              zero_field(pressure_energy_e_p);
              assembled =
                  product.ibm_equations->correct_zero_normal_diffusion(
                      as_const(pressure_energy_candidate_temperature),
                      as_const(
                          pressure_energy_candidate_thermal_conductivity),
                      pressure_energy_e_p);
            }
            for (std::int32_t z = 0; z < cells.z && assembled; ++z)
              for (std::int32_t y = 0; y < cells.y && assembled; ++y)
                for (std::int32_t x = 0; x < cells.x; ++x) {
                  const Int3 cell{x, y, z};
                  pressure_energy_r_e.unchecked(cell, 0U) -=
                      detail::cell_volume(product.equations.kernels(), cell) *
                      pressure_energy_e_p.unchecked(cell, 0U);
                }
          }
          if (assembled && product.topology.has_value()) {
            const Span<const std::uint8_t> active =
                product.topology->region();
            std::size_t active_offset = 0U;
            for (std::int32_t z = 0; z < cells.z; ++z)
              for (std::int32_t y = 0; y < cells.y; ++y)
                for (std::int32_t x = 0; x < cells.x;
                     ++x, ++active_offset)
                  if (active.data[active_offset] == 0U)
                    pressure_energy_r_e.unchecked({x, y, z}, 0U) = 0.0;
          }
          return assembled;
        };
        local = assemble_candidate_energy();
        local = product.reductions.consensus(local);
        if (!local) return local;
        if (alpha == 0.0) {
          const std::size_t expected_energy_values =
              static_cast<std::size_t>(cells.x) *
              static_cast<std::size_t>(cells.y) *
              static_cast<std::size_t>(cells.z);
          double local_energy_error = 0.0;
          bool local_energy_bitwise_equal = true;
          std::size_t energy_offset = 0U;
          if (pressure_energy_baseline_energy_residual.size() !=
              expected_energy_values)
            local = {StatusCode::invalid_plan, kProductBinding};
          // Open-boundary alpha zero is defined by the independent finalizer:
          // the same-layer EOS-closed state plus its final physical boundary
          // flux.  The pre-finalizer residual was assembled with provisional
          // mechanical flux and is not that physical baseline.  Preserve the
          // first full assembly.  An armed diagnostic repeats the complete
          // energy path and demands deterministic equality; production uses
          // the first physical assembly as both baseline and oracle.
          if (local && pressure_reference_kind ==
                           PressureReferenceKind::boundary_absolute) {
            for (std::int32_t z = 0; z < cells.z; ++z)
              for (std::int32_t y = 0; y < cells.y; ++y)
                for (std::int32_t x = 0; x < cells.x;
                     ++x, ++energy_offset)
                  pressure_energy_baseline_energy_residual[energy_offset] =
                      pressure_energy_r_e.unchecked({x, y, z}, 0U);
            if (full_alpha_zero_oracle) {
              local = assemble_candidate_energy();
              local = product.reductions.consensus(local);
              if (!local) return local;
            }
            energy_offset = 0U;
          }
          for (std::int32_t z = 0; z < cells.z && local; ++z)
            for (std::int32_t y = 0; y < cells.y && local; ++y)
              for (std::int32_t x = 0; x < cells.x;
                   ++x, ++energy_offset) {
                const double assembled_residual =
                    pressure_energy_r_e.unchecked({x, y, z}, 0U);
                const double oracle =
                    !full_alpha_zero_oracle &&
                            pressure_reference_kind ==
                                PressureReferenceKind::boundary_absolute
                        ? pressure_energy_baseline_energy_residual[energy_offset]
                        : assembled_residual;
                if (full_alpha_zero_oracle) {
                  const double base =
                      pressure_energy_baseline_energy_residual[energy_offset];
                  local_energy_bitwise_equal &=
                      product_double_bits(oracle) ==
                      product_double_bits(base);
                  const Int3 cell{x, y, z};
                  const double volume =
                      detail::cell_volume(product.equations.kernels(), cell);
                  double temporal_operand_scale = 0.0;
                  if (!pressure_energy_temporal_scale(
                          artifacts.pressure_reference,
                          as_const(pressure_energy_candidate_density),
                          as_const(pressure_energy_candidate_enthalpy),
                          as_const(pressure_energy_candidate_pressure), cell,
                          volume, temporal_operand_scale)) {
                    local = {StatusCode::numerical_failure,
                             kProductPressureEnergy};
                    break;
                  }
                  // The residual may be a near-zero cancellation of large
                  // temporal and pressure operands.  Scale its replay error by
                  // those operands, not by the cancelled residual itself.
                  const double scale = std::max(
                      {1.0, std::abs(oracle), std::abs(base),
                       temporal_operand_scale});
                  local_energy_error = std::max(
                      local_energy_error, std::abs(oracle - base) / scale);
                }
                local_alpha_zero_oracle_lineage = detail::product_mix(
                    local_alpha_zero_oracle_lineage,
                    product_double_bits(oracle));
              }
          if (full_alpha_zero_oracle) {
            double global_energy_error = 0.0;
            const bool gauge_representation_shifted =
                pressure_reference_kind == PressureReferenceKind::closed_mass &&
                artifacts.exact.closed_gauge.shift != 0.0;
            // A closed-mass gauge shift preserves absolute pressure but can
            // change the rounded pressure-work cancellation when p_ref + pi
            // is reconstructed on a different rank decomposition.  The
            // pressure oracle below independently certifies that gauge
            // equivalence.  Keep the stronger byte oracle for an unshifted
            // representation; for a shifted gauge require the recorded
            // residual difference to remain inside an independent
            // machine-roundoff bound.  The physical EOS tolerance is much
            // too wide to certify replay equivalence.
            const bool local_energy_oracle_equivalent =
                detail::alpha_zero_energy_replay_equivalent(
                    local_energy_bitwise_equal,
                    gauge_representation_shifted, local_energy_error);
            evaluated = product.reductions.checked_max(
                {&local_energy_error, 1U}, {&global_energy_error, 1U},
                local && local_energy_oracle_equivalent &&
                        std::isfinite(local_energy_error) &&
                        local_energy_error <=
                            detail::kAlphaZeroEnergyReplayRoundoffTolerance
                    ? Status{}
                    : Status{StatusCode::rejected_step,
                             kProductPressureEnergy});
            if (!evaluated) return evaluated;
            artifacts.alpha_zero_energy_residual_oracle_error =
                global_energy_error;
          }
        }
        evaluated = product.reductions.consensus(local);
        if (!evaluated) return evaluated;

        double local_metrics[2U]{};
        std::size_t offset = 0U;
        for (std::int32_t z = 0; z < cells.z && local; ++z) {
          for (std::int32_t y = 0; y < cells.y && local; ++y) {
            for (std::int32_t x = 0; x < cells.x; ++x, ++offset) {
              if (!pressure_energy_active(offset)) continue;
              const Int3 cell{x, y, z};
              const double rho =
                  pressure_energy_candidate_density.unchecked(cell, 0U);
              const double rho_n =
                  rho_history.accepted.unchecked(cell, 0U);
              const double rho_nm1 =
                  effective_bdf.order == 2U
                      ? rho_history.previous.unchecked(cell, 0U)
                      : 0.0;
              const double volume =
                  detail::cell_volume(product.equations.kernels(), cell);
              double eos_residual = 0.0;
              double continuity_residual = 0.0;
              double mass = 0.0;
              double volume_sum = 0.0;
              double absolute_pi = 0.0;
              double compressibility_moment = 0.0;
              double compressibility_weight = 0.0;
              const int terminal = hf_coast_common_terminal_cell_v1(
                  rho, rho, rho_n, rho_nm1, volume, effective_bdf.a0,
                  effective_bdf.a1, effective_bdf.a2,
                  artifacts.flux.x.unchecked(cell),
                  artifacts.flux.x.unchecked({x + 1, y, z}),
                  artifacts.flux.y.unchecked(cell),
                  artifacts.flux.y.unchecked({x, y + 1, z}),
                  artifacts.flux.z.unchecked(cell),
                  artifacts.flux.z.unchecked({x, y, z + 1}),
                  pressure_energy_candidate_pressure.unchecked(cell, 0U),
                  pressure_reference_kind == PressureReferenceKind::closed_mass
                      ? 1
                      : 0,
                  pressure_energy_candidate_compressibility.unchecked(cell,
                                                                        0U),
                  &eos_residual, &continuity_residual, &mass, &volume_sum,
                  &absolute_pi, &compressibility_moment,
                  &compressibility_weight);
              const double h =
                  pressure_energy_candidate_enthalpy.unchecked(cell, 0U);
              const double temperature =
                  pressure_energy_candidate_temperature.unchecked(cell, 0U);
              const double cp =
                  pressure_energy_candidate_heat_capacity.unchecked(cell,
                                                                      0U);
              const double p_abs =
                  artifacts.pressure_reference +
                  pressure_energy_candidate_pressure.unchecked(cell, 0U);
              double energy_scale = 0.0;
              const bool energy_scale_valid =
                  pressure_energy_temporal_scale(
                      artifacts.pressure_reference,
                      as_const(pressure_energy_candidate_density),
                      as_const(pressure_energy_candidate_enthalpy),
                      as_const(pressure_energy_candidate_pressure), cell,
                      volume, energy_scale);
              const double energy_metric =
                  std::abs(pressure_energy_r_e.unchecked(cell, 0U)) /
                  energy_scale;
              if (terminal != 0 || !std::isfinite(continuity_residual) ||
                  !std::isfinite(h) || !std::isfinite(temperature) ||
                  !(temperature > 0.0) || !std::isfinite(cp) ||
                  !(cp > 0.0) || !std::isfinite(p_abs) ||
                  !(p_abs > 0.0) || !energy_scale_valid ||
                  !std::isfinite(energy_metric)) {
                local = {StatusCode::rejected_step,
                         kProductPressureEnergy};
                break;
              }
              local_metrics[0U] =
                  std::max(local_metrics[0U], continuity_residual);
              local_metrics[1U] =
                  std::max(local_metrics[1U], energy_metric);
            }
          }
        }
        double global_metrics[2U]{};
        evaluated = product.reductions.checked_max(
            {local_metrics, 2U}, {global_metrics, 2U}, local);
        if (!evaluated) return evaluated;

        PressureEnergyGlobalizationSample& sample = artifacts.sample;
        sample.global_normalized_continuity = global_metrics[0U];
        sample.global_normalized_energy = global_metrics[1U];
        sample.thermodynamically_admissible = true;
        sample.state_and_flux_finite = true;

        bool alpha_zero_pressure_equivalent = true;
        bool alpha_zero_gauge_certified = true;
        if (alpha == 0.0) {
          double local_pressure_error[2U]{};
          Status local_pressure_status;
          if (pressure_reference_kind == PressureReferenceKind::closed_mass) {
            alpha_zero_gauge_certified =
                artifacts.exact.closed_gauge.valid() &&
                artifacts.exact.closed_gauge.predecessor_pressure_reference ==
                    pressure_reference_certificate.pressure_reference &&
                product_double_bits(artifacts.pressure_reference) ==
                    product_double_bits(
                        artifacts.exact.closed_gauge.next_pressure_reference);
          }
          for (std::int32_t z = 0; z < cells.z && local_pressure_status; ++z) {
            for (std::int32_t y = 0; y < cells.y && local_pressure_status; ++y) {
              for (std::int32_t x = 0; x < cells.x; ++x) {
                const Int3 cell{x, y, z};
                const double candidate_pi =
                    pressure_energy_candidate_pressure.unchecked(cell, 0U);
                const double base_pi = trial_pressure.unchecked(cell, 0U);
                if (pressure_reference_kind !=
                    PressureReferenceKind::closed_mass) {
                  if (product_double_bits(candidate_pi) !=
                      product_double_bits(base_pi)) {
                    local_pressure_status = {
                        StatusCode::invalid_plan, kProductPressureEnergy};
                    break;
                  }
                  continue;
                }
                const double candidate_absolute =
                    artifacts.pressure_reference + candidate_pi;
                const double base_absolute =
                    attempt_pressure_reference + base_pi;
                const double difference =
                    std::abs(candidate_absolute - base_absolute);
                const double scale = std::max(
                    {1.0, std::abs(candidate_absolute),
                     std::abs(base_absolute)});
                const double relative = difference / scale;
                if (!std::isfinite(difference) || !std::isfinite(relative) ||
                    relative > alpha_zero_oracle_tolerance) {
                  local_pressure_status = {
                      StatusCode::invalid_plan, kProductPressureEnergy};
                  break;
                }
                if (full_alpha_zero_oracle) {
                  local_pressure_error[0U] =
                      std::max(local_pressure_error[0U], difference);
                  local_pressure_error[1U] =
                      std::max(local_pressure_error[1U], relative);
                }
                local_alpha_zero_oracle_lineage = detail::product_mix(
                    local_alpha_zero_oracle_lineage,
                    product_double_bits(candidate_absolute));
                local_alpha_zero_oracle_lineage = detail::product_mix(
                    local_alpha_zero_oracle_lineage,
                    product_double_bits(base_absolute));
              }
            }
          }
          if (full_alpha_zero_oracle) {
            double global_pressure_error[2U]{};
            evaluated = product.reductions.checked_max(
                {local_pressure_error, 2U}, {global_pressure_error, 2U},
                local_pressure_status && alpha_zero_gauge_certified
                    ? Status{}
                    : Status{StatusCode::invalid_plan,
                             kProductPressureEnergy});
            if (!evaluated) return evaluated;
            artifacts.alpha_zero_pressure_oracle_difference =
                global_pressure_error[0U];
            artifacts.alpha_zero_pressure_oracle_relative_error =
                global_pressure_error[1U];
          } else {
            evaluated = product.reductions.consensus(
                local_pressure_status && alpha_zero_gauge_certified
                    ? Status{}
                    : Status{StatusCode::invalid_plan,
                             kProductPressureEnergy});
            if (!evaluated) return evaluated;
          }
          std::uint64_t lineage_xor = 0U;
          std::uint64_t lineage_sum = 0U;
          int communicator_size = 0;
          const int xor_result = MPI_Allreduce(
              &local_alpha_zero_oracle_lineage, &lineage_xor, 1,
              MPI_UINT64_T, MPI_BXOR, communicator);
          const int sum_result = MPI_Allreduce(
              &local_alpha_zero_oracle_lineage, &lineage_sum, 1,
              MPI_UINT64_T, MPI_SUM, communicator);
          const int size_result =
              MPI_Comm_size(communicator, &communicator_size);
          Status oracle_collective_status =
              xor_result == MPI_SUCCESS && sum_result == MPI_SUCCESS &&
                      size_result == MPI_SUCCESS && communicator_size > 0
                  ? Status{}
                  : Status{StatusCode::mpi_failure, kProductCollective};
          oracle_collective_status =
              product.reductions.consensus(oracle_collective_status);
          if (!oracle_collective_status) return oracle_collective_status;
          std::uint64_t oracle_lineage = detail::product_mix(
              UINT64_C(0x61306f72636c676c), lineage_xor);
          oracle_lineage = detail::product_mix(oracle_lineage, lineage_sum);
          oracle_lineage = detail::product_mix(
              oracle_lineage, static_cast<std::uint64_t>(communicator_size));
          oracle_lineage = detail::product_mix(
              oracle_lineage,
              product_double_bits(
                  pressure_reference_kind == PressureReferenceKind::closed_mass
                      ? artifacts.exact.closed_gauge.post_shift_gauge_residual
                      : 0.0));
          artifacts.alpha_zero_oracle_numeric_lineage =
              oracle_lineage == 0U ? PlanFingerprint{1U} : oracle_lineage;
          alpha_zero_pressure_equivalent = true;
        }

        bool alpha_zero_state_bits_match = alpha_zero_pressure_equivalent;
        bool alpha_zero_flux_bits_match = true;
        const auto has_cell_shape = [&](ConstFieldView field) noexcept {
          return field.interior.x == cells.x &&
                 field.interior.y == cells.y &&
                 field.interior.z == cells.z && field.components != 0U;
        };
        const auto mix_candidate_field_numeric =
            [&](std::uint64_t hash, ConstFieldView candidate,
                ConstFieldView live,
                const ConstFieldView* second_live = nullptr) noexcept {
              if (alpha != 0.0) return mix_field_numeric(hash, candidate);
              hash = detail::product_mix(hash, candidate.components);
              const bool candidate_shape = has_cell_shape(candidate);
              const bool live_shape =
                  has_cell_shape(live) &&
                  candidate.components == live.components;
              const bool second_shape =
                  second_live == nullptr ||
                  (has_cell_shape(*second_live) &&
                   candidate.components == second_live->components);
              const bool comparable =
                  candidate_shape && live_shape && second_shape;
              if (!comparable) alpha_zero_state_bits_match = false;
              if (!candidate_shape) return hash;
              for (std::uint8_t component = 0U;
                   component < candidate.components; ++component)
                for (std::int32_t z = 0; z < cells.z; ++z)
                  for (std::int32_t y = 0; y < cells.y; ++y)
                    for (std::int32_t x = 0; x < cells.x; ++x) {
                      const Int3 cell{x, y, z};
                      const std::uint64_t candidate_bits =
                          product_double_bits(
                              candidate.unchecked(cell, component));
                      hash = detail::product_mix(hash, candidate_bits);
                      if (comparable &&
                          (candidate_bits !=
                               product_double_bits(
                                   live.unchecked(cell, component)) ||
                           (second_live != nullptr &&
                            candidate_bits !=
                                product_double_bits(second_live->unchecked(
                                    cell, component)))))
                        alpha_zero_state_bits_match = false;
                    }
              return hash;
            };
        const auto mix_candidate_flux_numeric =
            [&](std::uint64_t hash, ConstFaceFluxView candidate,
                ConstFaceFluxView live) noexcept {
              if (alpha != 0.0) return mix_flux_numeric(hash, candidate);
              const std::array<ConstFaceFieldView, 3U> candidate_faces{
                  candidate.x, candidate.y, candidate.z};
              const std::array<ConstFaceFieldView, 3U> live_faces{
                  live.x, live.y, live.z};
              const std::array<Int3, 3U> expected_extents{{
                  {cells.x + 1, cells.y, cells.z},
                  {cells.x, cells.y + 1, cells.z},
                  {cells.x, cells.y, cells.z + 1},
              }};
              for (std::size_t axis = 0U; axis < candidate_faces.size();
                   ++axis) {
                const ConstFaceFieldView candidate_face =
                    candidate_faces[axis];
                const ConstFaceFieldView live_face = live_faces[axis];
                const Int3 expected = expected_extents[axis];
                const bool candidate_shape =
                    candidate_face.extents.x == expected.x &&
                    candidate_face.extents.y == expected.y &&
                    candidate_face.extents.z == expected.z;
                const bool comparable =
                    candidate_shape && live_face.extents.x == expected.x &&
                    live_face.extents.y == expected.y &&
                    live_face.extents.z == expected.z;
                if (!comparable) alpha_zero_flux_bits_match = false;
                if (!candidate_shape) continue;
                for (std::int32_t z = 0; z < candidate_face.extents.z; ++z)
                  for (std::int32_t y = 0; y < candidate_face.extents.y; ++y)
                    for (std::int32_t x = 0; x < candidate_face.extents.x;
                         ++x) {
                      const Int3 face{x, y, z};
                      const std::uint64_t candidate_bits = product_double_bits(
                          candidate_face.unchecked(face));
                      hash = detail::product_mix(hash, candidate_bits);
                      if (comparable &&
                          candidate_bits != product_double_bits(
                                                live_face.unchecked(face)))
                        alpha_zero_flux_bits_match = false;
                    }
              }
              return hash;
            };

        std::uint64_t residual = detail::product_mix(
            artifacts.exact_certificate.canonical_lineage(),
            UINT64_C(0x726573696475616c));
        residual = detail::product_mix(
            residual, product.turbulence.fingerprint());
        residual = detail::product_mix(
            residual, product.equations.semantic_fingerprint());
        residual = detail::product_mix(
            residual, product_double_bits(artifacts.pressure_reference));
        if (pressure_reference_kind == PressureReferenceKind::closed_mass) {
          residual = detail::product_mix(
              residual,
              product_double_bits(artifacts.exact.closed_gauge.shift));
          residual = detail::product_mix(
              residual,
              product_double_bits(artifacts.exact.closed_gauge.global_moment));
          residual = detail::product_mix(
              residual,
              product_double_bits(artifacts.exact.closed_gauge.global_weight));
        }
        residual = mix_field_numeric(
            residual, as_const(pressure_energy_candidate_pressure));
        residual = mix_candidate_field_numeric(
            residual, as_const(pressure_energy_candidate_enthalpy),
            as_const(trial_enthalpy));
        const ConstFieldView eos_density_baseline = as_const(eos_density);
        residual = mix_candidate_field_numeric(
            residual, as_const(pressure_energy_candidate_density),
            as_const(trial_density), &eos_density_baseline);
        residual = mix_candidate_field_numeric(
            residual, as_const(pressure_energy_candidate_temperature),
            as_const(trial_temperature));
        residual = mix_candidate_field_numeric(
            residual, as_const(pressure_energy_candidate_velocity),
            as_const(trial_velocity));
        residual = mix_candidate_field_numeric(
            residual,
            as_const(pressure_energy_candidate_molecular_viscosity),
            as_const(molecular_viscosity));
        residual = mix_candidate_field_numeric(
            residual,
            as_const(pressure_energy_candidate_effective_viscosity),
            as_const(effective_viscosity));
        residual = mix_candidate_field_numeric(
            residual,
            as_const(pressure_energy_candidate_velocity_gradient),
            as_const(velocity_gradient));
        residual = mix_candidate_field_numeric(
            residual, as_const(pressure_energy_candidate_compressibility),
            as_const(compressibility));
        residual = mix_candidate_field_numeric(
            residual,
            as_const(pressure_energy_candidate_enthalpy_compressibility),
            as_const(enthalpy_compressibility));
        residual = mix_candidate_field_numeric(
            residual,
            as_const(pressure_energy_candidate_thermal_conductivity),
            as_const(conductivity));
        residual = mix_candidate_field_numeric(
            residual, as_const(pressure_energy_candidate_heat_capacity),
            as_const(heat_capacity));
        residual = mix_candidate_field_numeric(
            residual,
            as_const(pressure_energy_candidate_enthalpy_diffusivity),
            as_const(enthalpy_diffusivity));
        residual = mix_field_numeric(residual, as_const(pressure_energy_r_e));
        residual = mix_candidate_flux_numeric(residual, artifacts.flux,
                                              baseline_flux);
        residual = detail::product_mix(
            residual, product_double_bits(global_metrics[0U]));
        residual = detail::product_mix(
            residual, product_double_bits(global_metrics[1U]));
        if (alpha == 0.0) {
          residual = detail::product_mix(
              residual,
              product_double_bits(
                  artifacts.alpha_zero_pressure_oracle_difference));
          residual = detail::product_mix(
              residual,
              product_double_bits(
                  artifacts.alpha_zero_pressure_oracle_relative_error));
          residual = detail::product_mix(
              residual,
              product_double_bits(
                  artifacts.alpha_zero_density_oracle_difference));
          residual = detail::product_mix(
              residual,
              product_double_bits(
                  artifacts.alpha_zero_density_oracle_relative_error));
          residual = detail::product_mix(
              residual, artifacts.alpha_zero_density_oracle_ulp);
          residual = detail::product_mix(
              residual,
              product_double_bits(
                  artifacts.alpha_zero_temperature_oracle_difference));
          residual = detail::product_mix(
              residual,
              product_double_bits(
                  artifacts.alpha_zero_temperature_oracle_relative_error));
          residual = detail::product_mix(
              residual, artifacts.alpha_zero_temperature_oracle_ulp);
          residual = detail::product_mix(
              residual,
              product_double_bits(artifacts.alpha_zero_material_oracle_error));
          residual = detail::product_mix(
              residual,
              product_double_bits(artifacts.alpha_zero_gradient_oracle_error));
          residual = detail::product_mix(
              residual,
              product_double_bits(
                  artifacts.alpha_zero_energy_residual_oracle_error));
          residual = detail::product_mix(
              residual,
              product_double_bits(
                  artifacts.alpha_zero_effective_viscosity_oracle_error));
          residual = detail::product_mix(
              residual, artifacts.alpha_zero_oracle_numeric_lineage);
        }
        artifacts.residual_replay_provenance =
            residual == 0U ? PlanFingerprint{1U} : residual;

        if (alpha == 0.0) {
          const double local_alpha_zero_mismatch =
              alpha_zero_state_bits_match && alpha_zero_gauge_certified &&
                      alpha_zero_flux_bits_match
                  ? 0.0
                  : 1.0;
          double global_alpha_zero_mismatch = 1.0;
          evaluated = product.reductions.checked_max(
              {&local_alpha_zero_mismatch, 1U},
              {&global_alpha_zero_mismatch, 1U});
          if (!evaluated) return evaluated;
          artifacts.alpha_zero_byte_equivalent =
              global_alpha_zero_mismatch == 0.0;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          if (g_candidate_globalization_published.load(
                  std::memory_order_acquire)) {
            detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic =
                g_candidate_globalization_diagnostic;
            diagnostic.alpha_zero_maximum_pressure_difference =
                artifacts.alpha_zero_pressure_oracle_difference;
            diagnostic.alpha_zero_maximum_pressure_relative_error =
                artifacts.alpha_zero_pressure_oracle_relative_error;
            diagnostic.alpha_zero_maximum_density_difference =
                artifacts.alpha_zero_density_oracle_difference;
            diagnostic.alpha_zero_maximum_density_relative_error =
                artifacts.alpha_zero_density_oracle_relative_error;
            diagnostic.alpha_zero_maximum_density_ulp_difference =
                artifacts.alpha_zero_density_oracle_ulp;
            diagnostic.alpha_zero_maximum_temperature_difference =
                artifacts.alpha_zero_temperature_oracle_difference;
            diagnostic.alpha_zero_maximum_temperature_relative_error =
                artifacts.alpha_zero_temperature_oracle_relative_error;
            diagnostic.alpha_zero_maximum_temperature_ulp_difference =
                artifacts.alpha_zero_temperature_oracle_ulp;
            diagnostic.alpha_zero_material_oracle_error =
                artifacts.alpha_zero_material_oracle_error;
            diagnostic.alpha_zero_gradient_oracle_error =
                artifacts.alpha_zero_gradient_oracle_error;
            diagnostic.alpha_zero_energy_residual_oracle_error =
                artifacts.alpha_zero_energy_residual_oracle_error;
            diagnostic.alpha_zero_effective_viscosity_oracle_error =
                artifacts.alpha_zero_effective_viscosity_oracle_error;
            diagnostic.alpha_zero_oracle_numeric_lineage =
                artifacts.alpha_zero_oracle_numeric_lineage;
            publish_candidate_globalization_diagnostic(diagnostic);
          }
#endif
        }
        return evaluated;
      };
  const auto preflight_pressure_energy_candidate = [&]() {
    Status admissible;
    std::size_t offset = 0U;
    for (std::int32_t z = 0; z < cells.z && admissible; ++z) {
      for (std::int32_t y = 0; y < cells.y && admissible; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x, ++offset) {
          const Int3 cell{x, y, z};
          const double delta_pressure =
              pressure_correction.unchecked(cell, 0U);
          const double delta_enthalpy =
              enthalpy_correction.unchecked(cell, 0U);
          const double candidate_pressure =
              trial_pressure.unchecked(cell, 0U) +
              delta_pressure;
          const double candidate_enthalpy =
              trial_enthalpy.unchecked(cell, 0U) +
              delta_enthalpy;
          if (!std::isfinite(candidate_pressure) ||
              !std::isfinite(candidate_enthalpy) ||
              !(attempt_pressure_reference + candidate_pressure > 0.0)) {
            admissible = {StatusCode::rejected_step,
                          kProductPressureEnergy};
            break;
          }
          if (!pressure_energy_active(offset)) {
            if (delta_pressure != 0.0 || delta_enthalpy != 0.0 ||
                !std::isfinite(trial_density.unchecked(cell, 0U)) ||
                !(trial_density.unchecked(cell, 0U) > 0.0) ||
                !std::isfinite(trial_temperature.unchecked(cell, 0U)) ||
                !(trial_temperature.unchecked(cell, 0U) > 0.0)) {
              admissible = {StatusCode::rejected_step,
                            kProductPressureEnergy};
              break;
            }
            pressure_energy_candidate_enthalpy.unchecked(cell, 0U) =
                candidate_enthalpy;
            pressure_energy_candidate_density.unchecked(cell, 0U) =
                trial_density.unchecked(cell, 0U);
            pressure_energy_candidate_temperature.unchecked(cell, 0U) =
                trial_temperature.unchecked(cell, 0U);
            continue;
          }
          for (std::size_t species = 0U;
               species < species_trial.size(); ++species)
            species_values[species] =
                species_trial[species].unchecked(cell, 0U);
          ThermoState thermo;
          Status evaluated =
              product.thermodynamics.evaluate_from_reference_pressure(
                  attempt_pressure_reference, candidate_pressure,
                  candidate_enthalpy,
                  {species_values.data(), species_values.size()},
                  {trial_velocity.unchecked(cell, 0U),
                   trial_velocity.unchecked(cell, 1U),
                   trial_velocity.unchecked(cell, 2U)},
                  thermo, trial_temperature.unchecked(cell, 0U));
          MolecularTransportState transport;
          if (evaluated)
            evaluated = product.transport.evaluate(
                thermo.temperature,
                {species_values.data(), species_values.size()}, transport);
          if (!evaluated || !std::isfinite(thermo.rho) ||
              !(thermo.rho > 0.0) ||
              !std::isfinite(thermo.temperature) ||
              !(thermo.temperature > 0.0) || !std::isfinite(thermo.cp) ||
              !(thermo.cp > 0.0) ||
              !std::isfinite(thermo.drho_dp_hY) ||
              !(thermo.drho_dp_hY > 0.0) ||
              !std::isfinite(thermo.drho_dh_pY) ||
              !(thermo.drho_dh_pY < 0.0) ||
              !std::isfinite(transport.viscosity) ||
              !(transport.viscosity > 0.0) ||
              !std::isfinite(transport.conductivity) ||
              !(transport.conductivity > 0.0)) {
            admissible = {StatusCode::rejected_step,
                          kProductPressureEnergy};
            break;
          }
          pressure_energy_candidate_enthalpy.unchecked(cell, 0U) =
              candidate_enthalpy;
          pressure_energy_candidate_density.unchecked(cell, 0U) =
              thermo.rho;
          pressure_energy_candidate_temperature.unchecked(cell, 0U) =
              thermo.temperature;
          molecular_viscosity.unchecked(cell, 0U) = transport.viscosity;
          effective_viscosity.unchecked(cell, 0U) = transport.viscosity;
          compressibility.unchecked(cell, 0U) = thermo.drho_dp_hY;
          enthalpy_compressibility.unchecked(cell, 0U) =
              thermo.drho_dh_pY;
          conductivity.unchecked(cell, 0U) = transport.conductivity;
          heat_capacity.unchecked(cell, 0U) = thermo.cp;
          enthalpy_diffusivity.unchecked(cell, 0U) =
              transport.conductivity / thermo.cp;
          eos_density.unchecked(cell, 0U) = thermo.rho;
        }
      }
    }
    return product.reductions.consensus(admissible);
  };
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const auto observe_pressure_energy_candidate_globalization =
      [&](std::uint8_t corrector) noexcept {
        if (!consume_candidate_globalization_arm()) return;

        constexpr std::size_t samples =
            detail::kPressureEnergyCandidateGlobalizationSampleCapacity;
        constexpr std::uint64_t absent_cell =
            std::numeric_limits<std::uint64_t>::max();
        constexpr std::uint32_t absent_reason =
            std::numeric_limits<std::uint32_t>::max();
        const double infinity = std::numeric_limits<double>::infinity();
        detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic;
        diagnostic.attempted_step =
            step.accepted_step == std::numeric_limits<std::uint64_t>::max()
                ? step.accepted_step
                : step.accepted_step + 1U;
        diagnostic.generation = step.generation;
        diagnostic.attempt = step.attempt;
        diagnostic.corrector = corrector;
        diagnostic.sample_count = static_cast<std::uint8_t>(samples);

        std::array<std::uint64_t, samples> local_first_cell{};
        std::array<std::uint64_t, samples> global_first_cell{};
        std::array<std::uint32_t, samples> local_first_reason{};
        std::array<std::uint32_t, samples> global_first_reason{};
        std::array<double, samples> local_minimum_pressure{};
        std::array<double, samples> global_minimum_pressure{};
        std::array<double, samples> local_minimum_temperature{};
        std::array<double, samples> global_minimum_temperature{};
        std::array<double, samples> local_minimum_density{};
        std::array<double, samples> global_minimum_density{};
        local_first_cell.fill(absent_cell);
        local_first_reason.fill(absent_reason);
        local_minimum_pressure.fill(infinity);
        local_minimum_temperature.fill(infinity);
        local_minimum_density.fill(infinity);

        const auto record_failure =
            [&](std::size_t sample, Int3 cell,
                detail::PressureEnergyCandidateFailureReason reason) {
              const Int3 global{product.patch.begin.x + cell.x,
                                product.patch.begin.y + cell.y,
                                product.patch.begin.z + cell.z};
              const std::uint64_t gid = diagnostic_global_cell(
                  global, product.geometry.global_cells());
              if (gid < local_first_cell[sample]) {
                local_first_cell[sample] = gid;
                local_first_reason[sample] =
                    static_cast<std::uint32_t>(reason);
              }
            };

        double local_maximum_correction[2U]{0.0, 0.0};
        for (std::int32_t z = 0; z < cells.z; ++z) {
          for (std::int32_t y = 0; y < cells.y; ++y) {
            for (std::int32_t x = 0; x < cells.x; ++x) {
              const Int3 cell{x, y, z};
              const double delta_pressure =
                  pressure_correction.unchecked(cell, 0U);
              const double delta_enthalpy =
                  enthalpy_correction.unchecked(cell, 0U);
              if (std::isfinite(delta_pressure))
                local_maximum_correction[0U] =
                    std::max(local_maximum_correction[0U],
                             std::abs(delta_pressure));
              else
                local_maximum_correction[0U] = infinity;
              if (std::isfinite(delta_enthalpy))
                local_maximum_correction[1U] =
                    std::max(local_maximum_correction[1U],
                             std::abs(delta_enthalpy));
              else
                local_maximum_correction[1U] = infinity;
            }
          }
        }

        // The fixed certificate covers alpha = 1, 1/2, ..., 2^-24.  This is
        // small enough to distinguish a globalization opportunity from an
        // inadmissible base state without silently treating alpha=0 as a
        // candidate correction.
        for (std::size_t sample = 0U; sample < samples; ++sample) {
          const double alpha = std::ldexp(1.0, -static_cast<int>(sample));
          diagnostic.samples[sample].alpha = alpha;
          std::size_t offset = 0U;
          for (std::int32_t z = 0; z < cells.z; ++z) {
            for (std::int32_t y = 0; y < cells.y; ++y) {
              for (std::int32_t x = 0; x < cells.x; ++x, ++offset) {
                const Int3 cell{x, y, z};
                const double delta_pressure =
                    pressure_correction.unchecked(cell, 0U);
                const double delta_enthalpy =
                    enthalpy_correction.unchecked(cell, 0U);
                const double candidate_pressure =
                    trial_pressure.unchecked(cell, 0U) +
                    alpha * delta_pressure;
                const double candidate_enthalpy =
                    trial_enthalpy.unchecked(cell, 0U) +
                    alpha * delta_enthalpy;
                const double absolute_pressure =
                    attempt_pressure_reference + candidate_pressure;
                if (std::isfinite(absolute_pressure))
                  local_minimum_pressure[sample] =
                      std::min(local_minimum_pressure[sample],
                               absolute_pressure);

                using Failure =
                    detail::PressureEnergyCandidateFailureReason;
                if (!std::isfinite(candidate_pressure)) {
                  record_failure(sample, cell, Failure::nonfinite_pressure);
                  continue;
                }
                if (!std::isfinite(candidate_enthalpy)) {
                  record_failure(sample, cell, Failure::nonfinite_enthalpy);
                  continue;
                }
                if (!(absolute_pressure > 0.0)) {
                  record_failure(sample, cell,
                                 Failure::nonpositive_absolute_pressure);
                  continue;
                }
                if (!pressure_energy_active(offset)) {
                  if (delta_pressure != 0.0 || delta_enthalpy != 0.0) {
                    record_failure(sample, cell,
                                   Failure::inactive_nonzero_correction);
                  } else if (!std::isfinite(
                                 trial_density.unchecked(cell, 0U)) ||
                             !(trial_density.unchecked(cell, 0U) > 0.0)) {
                    record_failure(sample, cell, Failure::inactive_density);
                  } else if (!std::isfinite(
                                 trial_temperature.unchecked(cell, 0U)) ||
                             !(trial_temperature.unchecked(cell, 0U) > 0.0)) {
                    record_failure(sample, cell,
                                   Failure::inactive_temperature);
                  } else {
                    local_minimum_density[sample] =
                        std::min(local_minimum_density[sample],
                                 trial_density.unchecked(cell, 0U));
                    local_minimum_temperature[sample] =
                        std::min(local_minimum_temperature[sample],
                                 trial_temperature.unchecked(cell, 0U));
                  }
                  continue;
                }
                for (std::size_t species = 0U;
                     species < species_trial.size(); ++species)
                  species_values[species] =
                      species_trial[species].unchecked(cell, 0U);
                ThermoState thermo;
                const Status evaluated =
                    product.thermodynamics.evaluate_from_reference_pressure(
                        attempt_pressure_reference, candidate_pressure,
                        candidate_enthalpy,
                        {species_values.data(), species_values.size()},
                        {trial_velocity.unchecked(cell, 0U),
                         trial_velocity.unchecked(cell, 1U),
                         trial_velocity.unchecked(cell, 2U)},
                        thermo, trial_temperature.unchecked(cell, 0U));
                if (!evaluated) {
                  record_failure(sample, cell,
                                 Failure::thermodynamic_evaluation);
                  continue;
                }
                if (std::isfinite(thermo.rho))
                  local_minimum_density[sample] =
                      std::min(local_minimum_density[sample], thermo.rho);
                if (std::isfinite(thermo.temperature))
                  local_minimum_temperature[sample] =
                      std::min(local_minimum_temperature[sample],
                               thermo.temperature);
                if (!std::isfinite(thermo.rho) || !(thermo.rho > 0.0)) {
                  record_failure(sample, cell, Failure::density);
                  continue;
                }
                if (!std::isfinite(thermo.temperature) ||
                    !(thermo.temperature > 0.0)) {
                  record_failure(sample, cell, Failure::temperature);
                  continue;
                }
                if (!std::isfinite(thermo.cp) || !(thermo.cp > 0.0)) {
                  record_failure(sample, cell, Failure::heat_capacity);
                  continue;
                }
                if (!std::isfinite(thermo.drho_dp_hY) ||
                    !(thermo.drho_dp_hY > 0.0)) {
                  record_failure(sample, cell,
                                 Failure::pressure_compressibility);
                  continue;
                }
                if (!std::isfinite(thermo.drho_dh_pY) ||
                    !(thermo.drho_dh_pY < 0.0)) {
                  record_failure(sample, cell,
                                 Failure::enthalpy_compressibility);
                  continue;
                }
                MolecularTransportState transport;
                const Status transported = product.transport.evaluate(
                    thermo.temperature,
                    {species_values.data(), species_values.size()},
                    transport);
                if (!transported) {
                  record_failure(sample, cell,
                                 Failure::transport_evaluation);
                  continue;
                }
                if (!std::isfinite(transport.viscosity) ||
                    !(transport.viscosity > 0.0)) {
                  record_failure(sample, cell, Failure::viscosity);
                  continue;
                }
                if (!std::isfinite(transport.conductivity) ||
                    !(transport.conductivity > 0.0)) {
                  record_failure(sample, cell, Failure::conductivity);
                }
              }
            }
          }
        }

        double global_maximum_correction[2U]{};
        bool collected = true;
        collected =
            (MPI_Allreduce(local_maximum_correction,
                           global_maximum_correction, 2, MPI_DOUBLE, MPI_MAX,
                           communicator) == MPI_SUCCESS) &&
            collected;
        collected =
            (MPI_Allreduce(local_first_cell.data(), global_first_cell.data(),
                           static_cast<int>(samples), MPI_UINT64_T, MPI_MIN,
                           communicator) == MPI_SUCCESS) &&
            collected;
        collected =
            (MPI_Allreduce(local_minimum_pressure.data(),
                           global_minimum_pressure.data(),
                           static_cast<int>(samples), MPI_DOUBLE, MPI_MIN,
                           communicator) == MPI_SUCCESS) &&
            collected;
        collected =
            (MPI_Allreduce(local_minimum_temperature.data(),
                           global_minimum_temperature.data(),
                           static_cast<int>(samples), MPI_DOUBLE, MPI_MIN,
                           communicator) == MPI_SUCCESS) &&
            collected;
        collected =
            (MPI_Allreduce(local_minimum_density.data(),
                           global_minimum_density.data(),
                           static_cast<int>(samples), MPI_DOUBLE, MPI_MIN,
                           communicator) == MPI_SUCCESS) &&
            collected;
        for (std::size_t sample = 0U; sample < samples; ++sample)
          if (local_first_cell[sample] != global_first_cell[sample])
            local_first_reason[sample] = absent_reason;
        collected =
            (MPI_Allreduce(local_first_reason.data(),
                           global_first_reason.data(),
                           static_cast<int>(samples), MPI_UINT32_T, MPI_MIN,
                           communicator) == MPI_SUCCESS) &&
            collected;
        if (!collected) {
          publish_candidate_globalization_diagnostic(diagnostic);
          return;
        }

        diagnostic.maximum_absolute_pressure_correction =
            global_maximum_correction[0U];
        diagnostic.maximum_absolute_enthalpy_correction =
            global_maximum_correction[1U];
        for (std::size_t sample = 0U; sample < samples; ++sample) {
          auto& result = diagnostic.samples[sample];
          result.admissible = global_first_cell[sample] == absent_cell;
          result.first_failing_global_cell = global_first_cell[sample];
          result.first_failure_reason =
              result.admissible
                  ? detail::PressureEnergyCandidateFailureReason::none
                  : static_cast<
                        detail::PressureEnergyCandidateFailureReason>(
                        global_first_reason[sample]);
          result.minimum_absolute_pressure =
              global_minimum_pressure[sample];
          result.minimum_temperature = global_minimum_temperature[sample];
          result.minimum_density = global_minimum_density[sample];
          if (result.admissible &&
              diagnostic.first_admissible_sample ==
                  std::numeric_limits<std::uint8_t>::max())
            diagnostic.first_admissible_sample =
                static_cast<std::uint8_t>(sample);
        }
        diagnostic.valid = true;
        publish_candidate_globalization_diagnostic(diagnostic);
      };
#endif
  struct PressureEnergyCandidateLoopResult {
    PisoFrozenMomentumStageAuthority authority{};
    PisoFrozenMomentumExactCandidateCertificate exact_baseline{};
    PressureEnergyCandidateArtifacts replay{};
    PressureEnergyGlobalizationSelectionCertificate selection{};
    PressureEnergyStationaryCertificate stationary{};
    PlanFingerprint direction{};
    double maximum_pressure_direction{};
    double maximum_enthalpy_direction{};
    double selected_alpha{};
    std::uint8_t selected_ordinal{};
    std::uint8_t evaluated_candidates{};
    std::uint64_t runtime_halo_messages{};
    std::uint64_t runtime_halo_bytes{};
    std::uint64_t sealed_halo_messages{};
    std::uint64_t sealed_halo_bytes{};
    bool replay_valid{};
  };
  const auto pressure_energy_loop_merit =
      [](const PressureEnergyCandidateLoopResult& loop,
         double& merit) noexcept {
        merit = 0.0;
        return loop.replay_valid && detail::product_pressure_coupled_merit(
            loop.replay.sample.global_normalized_continuity,
            loop.replay.sample.global_normalized_energy, merit);
      };
  const auto pressure_energy_components_converged =
      [&](const PressureEnergyCandidateLoopResult& loop) noexcept {
        return loop.replay_valid &&
               loop.replay.sample.thermodynamically_admissible &&
               loop.replay.sample.state_and_flux_finite &&
               std::isfinite(
                   loop.replay.sample.global_normalized_continuity) &&
               std::isfinite(loop.replay.sample.global_normalized_energy) &&
               loop.replay.sample.global_normalized_continuity <=
                   product.summary.terminal_continuity_tolerance &&
               loop.replay.sample.global_normalized_energy <=
                   product.summary.terminal_continuity_tolerance;
      };
  const auto pressure_inexact_control =
      [&](const PressureEnergyCandidateLoopResult* loop,
          double previous_merit,
          bool previous_merit_available) noexcept {
        detail::ProductPressureInexactForcingState forcing;
        forcing.previous_merit = previous_merit;
        forcing.terminal_tolerance =
            product.summary.terminal_continuity_tolerance;
        forcing.previous_merit_available = previous_merit_available;
        double merit = 0.0;
        if (loop != nullptr && pressure_energy_loop_merit(*loop, merit)) {
          forcing.normalized_continuity =
              loop->replay.sample.global_normalized_continuity;
          forcing.normalized_energy =
              loop->replay.sample.global_normalized_energy;
          forcing.residual_available = true;
        }
        return detail::product_pressure_inexact_forcing_control(
            product.piso.pressure_solve(), forcing);
      };
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const auto publish_pressure_energy_candidate_loop_diagnostic =
      [&](const PressureEnergyGlobalizationSample& baseline,
          const std::array<PressureEnergyGlobalizationSample,
                           kPressureEnergyGlobalizationCandidateCount>&
              candidates,
          const PressureEnergyCandidateLoopResult& loop,
          bool selection_valid, bool committed) noexcept {
        if (!g_candidate_globalization_published.load(
                std::memory_order_acquire))
          return;
        detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic =
            g_candidate_globalization_diagnostic;
        diagnostic.production_candidate_loop = true;
        diagnostic.corrector = baseline.corrector;
        diagnostic.sample_count = loop.evaluated_candidates;
        diagnostic.baseline_commit = loop.stationary.valid();
        diagnostic.selection_valid = selection_valid;
        diagnostic.replay_valid = loop.replay_valid;
        diagnostic.committed = committed;
        diagnostic.baseline_normalized_continuity =
            baseline.global_normalized_continuity;
        diagnostic.baseline_normalized_energy =
            baseline.global_normalized_energy;
        diagnostic.selected_normalized_continuity =
            loop.replay.sample.global_normalized_continuity;
        diagnostic.selected_normalized_energy =
            loop.replay.sample.global_normalized_energy;
        diagnostic.selected_alpha = loop.selected_alpha;
        if (baseline.corrector == 1U) {
          diagnostic.corrector_one_baseline_normalized_continuity =
              baseline.global_normalized_continuity;
          diagnostic.corrector_one_baseline_normalized_energy =
              baseline.global_normalized_energy;
          diagnostic.corrector_one_selected_normalized_continuity =
              loop.replay.sample.global_normalized_continuity;
          diagnostic.corrector_one_selected_normalized_energy =
              loop.replay.sample.global_normalized_energy;
          diagnostic.corrector_one_selected_alpha = loop.selected_alpha;
        }
        if (pressure_energy_linear_predictions[0U].valid) {
          diagnostic.corrector_one_linear_predicted_normalized_continuity =
              pressure_energy_linear_predictions[0U].normalized_continuity;
          diagnostic.corrector_one_linear_predicted_normalized_energy =
              pressure_energy_linear_predictions[0U].normalized_energy;
        }
        if (pressure_energy_linear_predictions[1U].valid) {
          diagnostic.corrector_two_linear_predicted_normalized_continuity =
              pressure_energy_linear_predictions[1U].normalized_continuity;
          diagnostic.corrector_two_linear_predicted_normalized_energy =
              pressure_energy_linear_predictions[1U].normalized_energy;
        }
        diagnostic.correction_direction = loop.direction;
        diagnostic.selected_state_provenance =
            loop.replay.sample.state_provenance;
        diagnostic.selected_mass_flux_provenance =
            loop.replay.sample.mass_flux_provenance;
        diagnostic.final_boundary_flux_certified =
            loop.replay.final_boundary.valid() &&
            loop.replay.exact_certificate.final_boundary_flux_certified() &&
            loop.replay.exact_certificate.final_boundary_flux_provenance() ==
                loop.replay.final_boundary.final_flux_provenance();
        if (diagnostic.final_boundary_flux_certified) {
          diagnostic.final_boundary_canonical_lineage =
              loop.replay.final_boundary.canonical_lineage();
          diagnostic.final_physical_flux_provenance =
              loop.replay.final_boundary.final_flux_provenance();
          diagnostic.final_inlet_mass_target =
              loop.replay.final_boundary.inlet_mass_target();
        diagnostic.final_inlet_mass_achieved =
              loop.replay.final_boundary.inlet_mass_achieved();
        }
        diagnostic.candidate_runtime_halo_messages =
            loop.runtime_halo_messages;
        diagnostic.candidate_runtime_halo_bytes = loop.runtime_halo_bytes;
        diagnostic.candidate_sealed_halo_messages =
            loop.sealed_halo_messages;
        diagnostic.candidate_sealed_halo_bytes = loop.sealed_halo_bytes;
        if (product.topology.has_value()) {
          const Span<const ImmersedLink> links = product.topology->links();
          diagnostic.local_ibm_interface_links = links.size;
          for (std::size_t index = 0U; index < links.size; ++index) {
            const ImmersedLink& link = links.data[index];
            const Int3 cell = link.fluid_local_index;
            double value = 0.0;
            switch (link.direction) {
              case ImmersedFaceDirection::x_negative:
                value = loop.replay.flux.x.unchecked(cell);
                break;
              case ImmersedFaceDirection::x_positive:
                value = loop.replay.flux.x.unchecked(
                    {cell.x + 1, cell.y, cell.z});
                break;
              case ImmersedFaceDirection::y_negative:
                value = loop.replay.flux.y.unchecked(cell);
                break;
              case ImmersedFaceDirection::y_positive:
                value = loop.replay.flux.y.unchecked(
                    {cell.x, cell.y + 1, cell.z});
                break;
              case ImmersedFaceDirection::z_negative:
                value = loop.replay.flux.z.unchecked(cell);
                break;
              case ImmersedFaceDirection::z_positive:
                value = loop.replay.flux.z.unchecked(
                    {cell.x, cell.y, cell.z + 1});
                break;
            }
            if (product_double_bits(value) != product_double_bits(0.0)) {
              ++diagnostic.local_ibm_interface_nonzero_count;
              if (value == 0.0)
                ++diagnostic.local_ibm_interface_negative_zero_count;
            }
          }
        }
        diagnostic.first_admissible_sample =
            std::numeric_limits<std::uint8_t>::max();
        for (std::size_t index = 0U; index < candidates.size(); ++index) {
          auto& alpha_diagnostic = diagnostic.samples[index];
          if (index >= loop.evaluated_candidates) {
            alpha_diagnostic = {};
            continue;
          }
          alpha_diagnostic.alpha = candidates[index].alpha;
          alpha_diagnostic.admissible =
              candidates[index].thermodynamically_admissible &&
              candidates[index].state_and_flux_finite;
          if (baseline.corrector != 1U) {
            alpha_diagnostic.first_failing_global_cell =
                std::numeric_limits<std::uint64_t>::max();
            alpha_diagnostic.first_failure_reason =
                alpha_diagnostic.admissible
                    ? detail::PressureEnergyCandidateFailureReason::none
                    : detail::PressureEnergyCandidateFailureReason::
                          production_candidate_evaluation;
            alpha_diagnostic.minimum_absolute_pressure =
                std::numeric_limits<double>::infinity();
            alpha_diagnostic.minimum_temperature =
                std::numeric_limits<double>::infinity();
            alpha_diagnostic.minimum_density =
                std::numeric_limits<double>::infinity();
          }
          if (alpha_diagnostic.admissible &&
              diagnostic.first_admissible_sample ==
                  std::numeric_limits<std::uint8_t>::max())
            diagnostic.first_admissible_sample =
                static_cast<std::uint8_t>(index);
          alpha_diagnostic.normalized_continuity =
              candidates[index].global_normalized_continuity;
          alpha_diagnostic.normalized_energy =
              candidates[index].global_normalized_energy;
          alpha_diagnostic.state_and_flux_finite =
              candidates[index].state_and_flux_finite;
        }
        diagnostic.valid = true;
        publish_candidate_globalization_diagnostic(diagnostic);
      };
#endif
  const auto run_pressure_energy_candidate_loop =
      [&](std::uint8_t corrector,
          std::uint8_t refinement_iteration,
          const PressureEnergyJacobianObservation& jacobian_observation,
          const PisoIntermediateCertificate& intermediate,
          const PressureCorrectionCertificate& pressure,
          ConstFaceFluxView baseline_flux,
          double extrapolated_alpha,
          bool stationary_only,
          PressureEnergyCandidateLoopResult& loop) {
        loop = {};
        const DriverResourceReport candidate_resource_before =
            resource_snapshot();
        const std::size_t expected_energy_values =
            static_cast<std::size_t>(cells.x) *
            static_cast<std::size_t>(cells.y) *
            static_cast<std::size_t>(cells.z);
        if (pressure_energy_baseline_energy_residual.size() !=
            expected_energy_values)
          return Status{StatusCode::invalid_plan, kProductBinding};
        // The live residual is an independent alpha-zero oracle only for an
        // armed diagnostic in a periodic/closed domain.  Open boundaries must
        // first install the candidate's final physical flux, so their baseline
        // is captured after the first candidate energy assembly instead.
        if (alpha_zero_diagnostic_oracle_armed() &&
            pressure_reference_kind !=
                PressureReferenceKind::boundary_absolute) {
          std::size_t baseline_energy_offset = 0U;
          for (std::int32_t z = 0; z < cells.z; ++z)
            for (std::int32_t y = 0; y < cells.y; ++y)
              for (std::int32_t x = 0; x < cells.x;
                   ++x, ++baseline_energy_offset)
                pressure_energy_baseline_energy_residual
                    [baseline_energy_offset] =
                        pressure_energy_r_e.unchecked({x, y, z}, 0U);
        }
        double local_direction[2U]{};
        Status local;
        loop.direction = pressure_energy_direction_fingerprint(
            corrector, pressure, local_direction, local);
        if (local &&
            (!std::isfinite(extrapolated_alpha) ||
             extrapolated_alpha < 1.0 ||
             extrapolated_alpha > kPressureEnergyAitkenMaximumAlpha))
          local = {StatusCode::invalid_plan, kProductPressureEnergy};
        double global_direction[2U]{};
        Status candidate_status = product.reductions.checked_max(
            {local_direction, 2U}, {global_direction, 2U}, local);
        if (!candidate_status) return candidate_status;
        loop.maximum_pressure_direction = global_direction[0U];
        loop.maximum_enthalpy_direction = global_direction[1U];

        PisoFrozenMomentumStageAuthority frozen;
        candidate_status =
            product.coupler.make_frozen_momentum_stage_authority(
                intermediate, pressure, frozen);
        if (!candidate_status)
          return Status{candidate_status.code,
                        kProductPressureEnergy + 101U};
        loop.authority = frozen;
        constexpr std::size_t baseline_ordinal =
            kPressureEnergyGlobalizationCandidateCount;
        PressureEnergyCandidateArtifacts baseline_artifacts;
        candidate_status = evaluate_pressure_energy_candidate(
            corrector, frozen, pressure, nullptr, loop.direction, 0.0,
            baseline_ordinal,
            baseline_flux, baseline_artifacts);
        if (!candidate_status) return candidate_status;
        loop.exact_baseline = baseline_artifacts.exact_certificate;
        const bool exact_baseline_binding =
            loop.exact_baseline.valid() &&
            loop.exact_baseline.alpha() == 0.0 &&
            loop.exact_baseline.corrector() == corrector &&
            loop.exact_baseline.target_time() == step.generation &&
            loop.exact_baseline.baseline_state_provenance() ==
                loop.exact_baseline.candidate_state_provenance() &&
            loop.exact_baseline.baseline_mass_flux_provenance() ==
                loop.exact_baseline.candidate_mass_flux_provenance() &&
            baseline_artifacts.sample.state_provenance ==
                loop.exact_baseline.baseline_state_provenance() &&
            baseline_artifacts.sample.mass_flux_provenance ==
                loop.exact_baseline.baseline_mass_flux_provenance();
        candidate_status = product.reductions.consensus(
            exact_baseline_binding
                ? Status{}
                : Status{StatusCode::invalid_plan,
                         kProductPressureEnergy});
        if (!candidate_status) return candidate_status;
        loop.direction = loop.exact_baseline.correction_direction();
        const PressureEnergyGlobalizationSample baseline =
            baseline_artifacts.sample;

        std::array<PressureEnergyGlobalizationSample,
                   kPressureEnergyGlobalizationCandidateCount>
            candidates{};
        std::array<PlanFingerprint,
                   kPressureEnergyGlobalizationCandidateCount>
            residual_replay_provenance{};
        const bool baseline_terminal_compatible =
            baseline.thermodynamically_admissible &&
            baseline.state_and_flux_finite &&
            baseline_artifacts.alpha_zero_byte_equivalent &&
            baseline.global_normalized_continuity <=
                product.summary.terminal_continuity_tolerance &&
            baseline.global_normalized_energy <=
                product.summary.terminal_continuity_tolerance;
        if (stationary_only && !baseline_terminal_compatible) {
          loop.replay = baseline_artifacts;
          return Status{};
        }
        bool selection_valid = false;
        bool selected_artifact_valid = false;
        // A terminal exact baseline is a typed no-op authority regardless of
        // the (unconsumed) raw Newton direction.  Even a uniform converged
        // state can produce roundoff-sized nonzero dp/dh; sending that
        // direction through Armijo would require an impossible strict descent
        // from an already-terminal merit.  The stationary certificate itself
        // still proves that the alpha-zero scaled corrections are exactly
        // zero and binds the complete exact baseline residual/state/flux.
        if (baseline_terminal_compatible) {
          candidate_status = product.coupler.certify_frozen_momentum_stationary(
              frozen, loop.exact_baseline,
              baseline.global_normalized_continuity,
              baseline.global_normalized_energy, product.reductions,
              loop.stationary);
          if (!candidate_status || !loop.stationary.valid())
            return candidate_status
                       ? Status{StatusCode::invalid_plan,
                                kProductPressureEnergy}
                       : candidate_status;
          loop.selected_alpha = 0.0;
          loop.selected_ordinal =
              static_cast<std::uint8_t>(baseline_ordinal);
          // The exact alpha-zero evaluation already occupies the candidate
          // scratch and carries the complete state/flux certificate.  Keep
          // that immutable artifact for publication instead of replaying the
          // same no-op candidate through every full-field product kernel.
          loop.replay = baseline_artifacts;
          selected_artifact_valid = true;
        } else {
          Status selection_status{StatusCode::rejected_step,
                                  kProductPressureEnergy};
          const auto exact_candidate_matches =
              [&](const PressureEnergyCandidateArtifacts& artifacts,
                  const PressureEnergyGlobalizationSample& sample,
                  double alpha) noexcept {
                const PisoFrozenMomentumExactCandidateCertificate& exact =
                    artifacts.exact_certificate;
                return exact.valid() && exact.alpha() == alpha &&
                       exact.corrector() == corrector &&
                       exact.target_time() == baseline.target_time &&
                       exact.correction_direction() ==
                           baseline.correction_direction &&
                       exact.baseline_state_provenance() ==
                           baseline.state_provenance &&
                       exact.baseline_mass_flux_provenance() ==
                           baseline.mass_flux_provenance &&
                       exact.candidate_state_provenance() ==
                           sample.state_provenance &&
                       exact.candidate_mass_flux_provenance() ==
                           sample.mass_flux_provenance;
              };
          if (extrapolated_alpha > 1.0) {
            PressureEnergyCandidateArtifacts extrapolated;
            candidate_status = evaluate_pressure_energy_candidate(
                corrector, frozen, pressure, &loop.exact_baseline,
                loop.direction, extrapolated_alpha, 0U, baseline_flux,
                extrapolated);
            const bool complete_sample =
                static_cast<bool>(candidate_status);
            if (!candidate_status &&
                candidate_status.code != StatusCode::rejected_step &&
                candidate_status.code != StatusCode::numerical_failure)
              return candidate_status;
            if (complete_sample) {
              candidates[0U] = extrapolated.sample;
              candidate_status = product.reductions.consensus(
                  exact_candidate_matches(extrapolated, candidates[0U],
                                          extrapolated_alpha)
                      ? Status{}
                      : Status{StatusCode::invalid_plan,
                               kProductPressureEnergy});
              if (!candidate_status) return candidate_status;
              residual_replay_provenance[0U] =
                  extrapolated.residual_replay_provenance;
              loop.evaluated_candidates = 1U;
              selection_status = select_pressure_energy_extrapolation(
                  baseline, candidates[0U], loop.selection);
              if (selection_status) {
                selection_valid = true;
                loop.replay = extrapolated;
                selected_artifact_valid = true;
              } else if (selection_status.code !=
                         StatusCode::rejected_step) {
                return selection_status;
              }
            } else {
              candidate_status = {};
            }
          }
          if (!selection_valid) {
            for (std::size_t index = 0U; index < candidates.size(); ++index) {
              const double alpha =
                  std::ldexp(1.0, -static_cast<int>(index));
              PressureEnergyCandidateArtifacts artifacts;
              candidate_status = evaluate_pressure_energy_candidate(
                  corrector, frozen, pressure, &loop.exact_baseline,
                  loop.direction, alpha, index, baseline_flux, artifacts);
              candidates[index] = artifacts.sample;
              const bool complete_sample =
                  static_cast<bool>(candidate_status);
              if (!candidate_status) {
                if (candidate_status.code != StatusCode::rejected_step &&
                    candidate_status.code != StatusCode::numerical_failure)
                  return candidate_status;
                PressureEnergyGlobalizationSample& rejected =
                    candidates[index];
                rejected.alpha = alpha;
                rejected.corrector = baseline.corrector;
                rejected.target_time = baseline.target_time;
                rejected.correction_direction =
                    baseline.correction_direction;
                rejected.global_normalized_continuity =
                    std::numeric_limits<double>::infinity();
                rejected.global_normalized_energy =
                    std::numeric_limits<double>::infinity();
                rejected.thermodynamically_admissible = false;
                rejected.state_and_flux_finite = false;
                PlanFingerprint state = baseline.state_provenance;
                PlanFingerprint flux = baseline.mass_flux_provenance;
                do {
                  ++state;
                  if (state == 0U) ++state;
                } while (state == baseline.state_provenance ||
                         std::any_of(
                             candidates.begin(), candidates.begin() + index,
                             [state](const PressureEnergyGlobalizationSample&
                                         prior) noexcept {
                               return prior.state_provenance == state;
                             }));
                do {
                  ++flux;
                  if (flux == 0U) ++flux;
                } while (flux == baseline.mass_flux_provenance ||
                         std::any_of(
                             candidates.begin(), candidates.begin() + index,
                             [flux](const PressureEnergyGlobalizationSample&
                                        prior) noexcept {
                               return prior.mass_flux_provenance == flux;
                             }));
                rejected.state_provenance = state;
                rejected.mass_flux_provenance = flux;
                candidate_status = {};
              } else {
                candidate_status = product.reductions.consensus(
                    exact_candidate_matches(artifacts, candidates[index], alpha)
                        ? Status{}
                        : Status{StatusCode::invalid_plan,
                                 kProductPressureEnergy});
                if (!candidate_status) return candidate_status;
              }
              residual_replay_provenance[index] =
                  complete_sample
                      ? artifacts.residual_replay_provenance
                      : PlanFingerprint{};
              loop.evaluated_candidates =
                  static_cast<std::uint8_t>(index + 1U);
              selection_status = select_pressure_energy_globalization(
                  baseline,
                  {candidates.data(), loop.evaluated_candidates},
                  loop.selection);
              if (selection_status) {
                selection_valid = true;
                // Selection is attempted immediately after this candidate was
                // evaluated.  Earlier admissible candidates would already
                // have terminated the loop, so the selected artifact is the
                // one still resident in the preallocated candidate scratch.
                // Promote its typed certificate rather than evaluating it a
                // second time.
                loop.replay = artifacts;
                selected_artifact_valid = complete_sample;
                break;
              }
              if (selection_status.code != StatusCode::rejected_step)
                return selection_status;
            }
          }
          pressure_energy_globalization.valid = true;
          pressure_energy_globalization.corrector = corrector;
          pressure_energy_globalization.sample_count =
              loop.evaluated_candidates;
          pressure_energy_globalization.maximum_absolute_pressure_correction =
              loop.maximum_pressure_direction;
          pressure_energy_globalization.maximum_absolute_enthalpy_correction =
              loop.maximum_enthalpy_direction;
          pressure_energy_globalization.baseline = baseline;
          pressure_energy_globalization.candidates = candidates;
          if (!selection_valid) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
            publish_pressure_energy_candidate_loop_diagnostic(
                baseline, candidates, loop, false, false);
#endif
            return selection_status;
          }
          loop.selected_alpha = loop.selection.alpha;
          loop.selected_ordinal = loop.selection.selected_halvings;
        }

        candidate_status = product.reductions.consensus(
            selected_artifact_valid
                ? Status{}
                : Status{StatusCode::invalid_plan,
                         kProductPressureEnergy});
        if (!candidate_status) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          publish_pressure_energy_candidate_loop_diagnostic(
              baseline, candidates, loop, selection_valid, false);
#endif
          return candidate_status;
        }
        const PressureEnergyGlobalizationSample& expected =
            loop.stationary.valid()
                ? baseline
                : candidates[loop.selection.selected_halvings];
        const PlanFingerprint expected_residual_replay =
            loop.stationary.valid()
                ? baseline_artifacts.residual_replay_provenance
                : residual_replay_provenance[
                      loop.selection.selected_halvings];
        const PressureEnergyGlobalizationSample& replay = loop.replay.sample;
        const bool matching_replay =
            replay.alpha == expected.alpha &&
            replay.corrector == expected.corrector &&
            replay.target_time == expected.target_time &&
            replay.correction_direction == expected.correction_direction &&
            replay.state_provenance == expected.state_provenance &&
            replay.mass_flux_provenance == expected.mass_flux_provenance &&
            product_double_bits(replay.global_normalized_continuity) ==
                product_double_bits(
                    expected.global_normalized_continuity) &&
            product_double_bits(replay.global_normalized_energy) ==
                product_double_bits(expected.global_normalized_energy) &&
            replay.thermodynamically_admissible ==
                expected.thermodynamically_admissible &&
            replay.state_and_flux_finite == expected.state_and_flux_finite &&
            loop.replay.residual_replay_provenance ==
                expected_residual_replay &&
            (!loop.stationary.valid() ||
             loop.replay.alpha_zero_byte_equivalent) &&
            (loop.stationary.valid() ||
             (loop.selection.valid() &&
              replay.state_provenance ==
                  loop.selection.candidate_state_provenance &&
              replay.mass_flux_provenance ==
                  loop.selection.candidate_mass_flux_provenance));
        candidate_status = product.reductions.consensus(
            matching_replay
                ? Status{}
                : Status{StatusCode::invalid_plan,
                         kProductPressureEnergy});
        loop.replay_valid = static_cast<bool>(candidate_status);
        const DriverResourceReport candidate_resource_after =
            resource_snapshot();
        const FrozenStage* const candidate_stage =
            product.graph.stage(corrector == 1U ? 40U : 50U);
        const bool monotonic_candidate_resources =
            candidate_resource_after.structured_messages >=
                candidate_resource_before.structured_messages &&
            candidate_resource_after.structured_bytes >=
                candidate_resource_before.structured_bytes &&
            candidate_resource_after.ibm_messages >=
                candidate_resource_before.ibm_messages &&
            candidate_resource_after.ibm_bytes >=
                candidate_resource_before.ibm_bytes &&
            ((candidate_resource_after.structured_bytes -
              candidate_resource_before.structured_bytes) %
                 2U ==
             0U);
        if (monotonic_candidate_resources && candidate_stage != nullptr) {
          const std::uint64_t structured_messages =
              candidate_resource_after.structured_messages -
              candidate_resource_before.structured_messages;
          const std::uint64_t structured_bytes =
              candidate_resource_after.structured_bytes -
              candidate_resource_before.structured_bytes;
          const std::uint64_t ibm_messages =
              candidate_resource_after.ibm_messages -
              candidate_resource_before.ibm_messages;
          const std::uint64_t ibm_bytes =
              candidate_resource_after.ibm_bytes -
              candidate_resource_before.ibm_bytes;
          loop.runtime_halo_messages =
              ibm_messages > UINT64_MAX - structured_messages
                  ? UINT64_MAX
                  : structured_messages + ibm_messages;
          loop.runtime_halo_bytes =
              ibm_bytes > UINT64_MAX - structured_bytes / 2U
                  ? UINT64_MAX
                  : structured_bytes / 2U + ibm_bytes;
          loop.sealed_halo_messages =
              candidate_stage->resources.merged_halo_messages;
          loop.sealed_halo_bytes =
              candidate_stage->resources.merged_halo_bytes;
        }
        const bool resources_within_seal =
            monotonic_candidate_resources && candidate_stage != nullptr &&
            loop.runtime_halo_messages <= loop.sealed_halo_messages &&
            loop.runtime_halo_bytes <= loop.sealed_halo_bytes;
        candidate_status = product.reductions.consensus(
            candidate_status && resources_within_seal
                ? Status{}
                : Status{StatusCode::invalid_plan,
                         kProductPressureEnergy});
        loop.replay_valid = loop.replay_valid && resources_within_seal &&
                            static_cast<bool>(candidate_status);
        if (candidate_status) {
          const std::size_t trajectory_index =
              pressure_energy_globalization.trajectory_count;
          const bool trajectory_available =
              trajectory_index <
              pressure_energy_globalization.trajectory.size();
          // trajectory_count is advanced on every rank by the same selected
          // candidate loop.  Capacity is therefore a deterministic local
          // hot-resource contract and does not require another collective.
          if (!trajectory_available)
            candidate_status = {StatusCode::invalid_plan,
                                kProductPressureEnergy};
          if (candidate_status) {
            PressureEnergyGlobalizationIterationReport& iteration =
                pressure_energy_globalization
                    .trajectory[trajectory_index];
            iteration.valid = true;
            iteration.corrector = corrector;
            iteration.refinement_iteration = refinement_iteration;
            iteration.jacobian_scope_valid = jacobian_observation.valid;
            iteration.jacobian_scope = jacobian_observation.scope;
            iteration.maximum_absolute_pressure_correction =
                loop.maximum_pressure_direction;
            iteration.maximum_absolute_enthalpy_correction =
                loop.maximum_enthalpy_direction;
            iteration.baseline = baseline;
            iteration.selected = loop.replay.sample;
            ++pressure_energy_globalization.trajectory_count;
          }
        }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        publish_pressure_energy_candidate_loop_diagnostic(
            baseline, candidates, loop, selection_valid, false);
#endif
        return candidate_status;
      };
  const auto revise_coupled_trial_views = [&]() {
    Status revised;
    const std::array<FieldId, 5U> fields{
        product.fields.velocity, product.fields.pressure,
        product.fields.rho, product.fields.enthalpy,
        product.fields.temperature};
    for (FieldId field : fields)
      if (revised) revised = transaction.revise_trial(field);
    if (revised)
      revised = product.layers.view(StateRole::trial,
                                    product.fields.velocity,
                                    trial_velocity);
    if (revised)
      revised = product.layers.view(StateRole::trial,
                                    product.fields.pressure,
                                    trial_pressure);
    if (revised)
      revised = product.layers.view(StateRole::trial, product.fields.rho,
                                    trial_density);
    if (revised)
      revised = product.layers.view(StateRole::trial,
                                    product.fields.enthalpy,
                                    trial_enthalpy);
    if (revised)
      revised = product.layers.view(StateRole::trial,
                                    product.fields.temperature,
                                    trial_temperature);
    return revised;
  };
  const auto revise_exact_thermodynamic_workspaces = [&]() {
    Status revised = runtime_write_view(product.fields.molecular_viscosity,
                                        molecular_viscosity);
    if (revised)
      revised = runtime_write_view(product.fields.effective_viscosity,
                                   effective_viscosity);
    if (revised)
      revised = runtime_write_view(product.fields.compressibility,
                                   compressibility);
    if (revised)
      revised = runtime_write_view(
          product.fields.enthalpy_compressibility,
          enthalpy_compressibility);
    if (revised)
      revised = runtime_write_view(product.fields.thermal_conductivity,
                                   conductivity);
    if (revised)
      revised = runtime_write_view(product.fields.heat_capacity,
                                   heat_capacity);
    if (revised)
      revised = runtime_write_view(product.fields.enthalpy_diffusivity,
                                   enthalpy_diffusivity);
    if (revised)
      revised = runtime_write_view(product.fields.eos_density,
                                   eos_density);
    return revised;
  };
  const auto same_field_shape = [](ConstFieldView source,
                                   ConstFieldView target) noexcept {
    return source.base != nullptr && target.base != nullptr &&
           source.interior.x == target.interior.x &&
           source.interior.y == target.interior.y &&
           source.interior.z == target.interior.z &&
           source.components == target.components;
  };
  const auto candidate_material_sources = [&]() noexcept {
    return std::array<ConstFieldView, 9U>{
        as_const(pressure_energy_candidate_molecular_viscosity),
        as_const(pressure_energy_candidate_effective_viscosity),
        as_const(pressure_energy_candidate_velocity_gradient),
        as_const(pressure_energy_candidate_compressibility),
        as_const(pressure_energy_candidate_enthalpy_compressibility),
        as_const(pressure_energy_candidate_thermal_conductivity),
        as_const(pressure_energy_candidate_heat_capacity),
        as_const(pressure_energy_candidate_enthalpy_diffusivity),
        as_const(pressure_energy_candidate_density)};
  };
  const auto candidate_material_destinations = [&]() noexcept {
    return std::array<FieldView, 9U>{
        molecular_viscosity, effective_viscosity, velocity_gradient,
        compressibility, enthalpy_compressibility, conductivity,
        heat_capacity, enthalpy_diffusivity, eos_density};
  };
  const auto prepare_pressure_energy_candidate_material_publication = [&]() {
    // Reserve every destination revision and validate the complete copy graph
    // before the coupler reaches its collective no-fail publication boundary.
    // The matching publish step below performs only unchecked scalar stores.
    Status prepared = runtime_write_view(
        product.fields.molecular_viscosity, molecular_viscosity);
    if (prepared)
      prepared = runtime_write_view(product.fields.effective_viscosity,
                                    effective_viscosity);
    if (prepared)
      prepared = runtime_write_view(product.fields.velocity_gradient,
                                    velocity_gradient);
    if (prepared)
      prepared = runtime_write_view(product.fields.compressibility,
                                    compressibility);
    if (prepared)
      prepared = runtime_write_view(product.fields.enthalpy_compressibility,
                                    enthalpy_compressibility);
    if (prepared)
      prepared = runtime_write_view(product.fields.thermal_conductivity,
                                    conductivity);
    if (prepared)
      prepared = runtime_write_view(product.fields.heat_capacity,
                                    heat_capacity);
    if (prepared)
      prepared = runtime_write_view(product.fields.enthalpy_diffusivity,
                                    enthalpy_diffusivity);
    if (prepared)
      prepared = runtime_write_view(product.fields.eos_density, eos_density);
    if (!prepared) return product.reductions.consensus(prepared);

    const auto sources = candidate_material_sources();
    const auto destinations = candidate_material_destinations();
    bool valid = true;
    for (std::size_t index = 0U; index < sources.size(); ++index) {
      valid = valid && same_field_shape(
                           sources[index], as_const(destinations[index]));
      for (FieldView destination : destinations)
        valid = valid && !detail::field_views_overlap(
                             sources[index], as_const(destination));
      for (std::size_t other = index + 1U; other < sources.size(); ++other) {
        valid = valid && !detail::field_views_overlap(
                             sources[index], sources[other]);
        valid = valid && !detail::field_views_overlap(
                             as_const(destinations[index]),
                             as_const(destinations[other]));
      }
    }
    return product.reductions.consensus(
        valid ? Status{}
              : Status{StatusCode::invalid_plan, kProductBinding});
  };
  const auto publish_pressure_energy_candidate_materials = [&]() noexcept {
    const auto sources = candidate_material_sources();
    const auto destinations = candidate_material_destinations();
    for (std::size_t index = 0U; index < sources.size(); ++index) {
      const ConstFieldView source = sources[index];
      const FieldView destination = destinations[index];
      for (std::uint8_t component = 0U; component < source.components;
           ++component)
        for (std::int32_t z = 0; z < cells.z; ++z)
          for (std::int32_t y = 0; y < cells.y; ++y)
            for (std::int32_t x = 0; x < cells.x; ++x)
              destination.unchecked({x, y, z}, component) =
                  source.unchecked({x, y, z}, component);
    }
  };
  const auto same_pressure_reference =
      [](const PressureReferenceCertificate& left,
         const PressureReferenceCertificate& right) noexcept {
        return left.plan == right.plan &&
               left.predictor == right.predictor &&
               left.thermodynamics == right.thermodynamics &&
               left.closure == right.closure && left.time == right.time &&
               left.pressure_reference == right.pressure_reference &&
               left.kind == right.kind;
      };
  const auto same_pressure_reference_semantics =
      [](const PressureReferenceCertificate& left,
         const PressureReferenceCertificate& right) noexcept {
        return left.valid() && right.valid() && left.plan == right.plan &&
               left.predictor == right.predictor &&
               left.thermodynamics == right.thermodynamics &&
               left.closure == right.closure && left.time == right.time &&
               left.kind == right.kind;
      };
  const auto candidate_reference_transition_valid =
      [&](const PressureEnergyCandidateArtifacts& candidate) noexcept {
        if (pressure_reference_kind ==
            PressureReferenceKind::boundary_absolute)
          return !candidate.exact.closed_gauge.valid() &&
                 candidate.final_boundary.valid() &&
                 product_double_bits(candidate.pressure_reference) ==
                     product_double_bits(attempt_pressure_reference);
        if (pressure_reference_kind != PressureReferenceKind::closed_mass)
          return false;
        const ClosedGaugeCorrectionCertificate& gauge =
            candidate.exact.closed_gauge;
        return gauge.valid() &&
               gauge.predecessor_pressure_reference ==
                   pressure_reference_certificate.pressure_reference &&
               same_pressure_reference_semantics(
                   gauge.output_pressure_reference,
                   pressure_reference_certificate) &&
               product_double_bits(candidate.pressure_reference) ==
                   product_double_bits(gauge.next_pressure_reference);
      };
  const bool periodic_pressure_energy_candidate_scope =
      product_entirely_periodic(product.boundary);
  const bool open_pressure_energy_candidate_scope =
      product_candidate_boundary_supported(product.boundary) &&
      product.candidate_boundary_finalizer.ready();
  const bool pressure_energy_candidate_scope =
      periodic_pressure_energy_candidate_scope ||
      open_pressure_energy_candidate_scope;
  // Boundary-absolute pressure/enthalpy publication is legal only through
  // the typed final-boundary candidate transaction.  Static/total-state and
  // NSCBC faces do not yet have that issuer; allowing them to fall through
  // to correct_coupled_* would silently revive the pressure-only split that
  // this path was introduced to eliminate.
  if (status &&
      pressure_reference_kind == PressureReferenceKind::boundary_absolute &&
      !pressure_energy_candidate_scope) {
    attempt_stage = 44U;
    status = product.reductions.consensus(
        {StatusCode::invalid_plan, kProductPressureEnergy});
  }
  LinearIdentity identity_one = linear_identity(pressure_one);
  MgCoefficientIdentity coefficient_one{pressure_one.state,
                                        pressure_one.state, 0.0};
  if (status && !pressure_mg_initialized) {
    attempt_stage = 43U;
    NativeCartesianMgSpec spec;
    status = product.coupler.make_native_pressure_mg_spec(
        communicator, identity_one, coefficient_one, spec);
    if (status && product.topology.has_value()) {
      spec.activity = {
          {product.pressure_mg_cell_activity.data(),
           product.pressure_mg_cell_activity.size()},
          {product.pressure_mg_x_activity.data(),
           product.pressure_mg_x_activity.size()},
          {product.pressure_mg_y_activity.data(),
           product.pressure_mg_y_activity.size()},
          {product.pressure_mg_z_activity.data(),
           product.pressure_mg_z_activity.size()},
          product.pressure_mg_activity_fingerprint,
          product.pressure_mg_activity_collective};
    }
    if (status) {
      const MgRuntimeServices services{
          &product.mg_halo, &product.reductions, &product.mg_workspace,
          {product.coarse_halo_pointers.data(),
           product.coarse_halo_pointers.size()}};
      status = product.coupler.compile_native_pressure_mg(
          pressure_one, spec, services, pressure_system, pressure_mg,
          &pressure_mg_counters);
    }
    if (status) pressure_mg_initialized = true;
  }
  if (status) attempt_stage = 44U;
  if (status && ibm_pressure_operator.has_value())
    status = ibm_pressure_operator->mask_solid_rhs(pressure_system.rhs);
  if (status && product.piso.coupling() != CouplingKind::simple)
    status = form_pressure_energy_blocks();
  PressureEnergyJacobianObservation jacobian_one;
  if (status)
    status = solve_pressure_energy(
        1U, 0U, intermediate_one, pressure_one, pressure_energy_target_flux,
        pressure_energy_frozen_enthalpy, identity_one, coefficient_one,
        pressure_inexact_control(nullptr, 0.0, false), jacobian_one);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (status) observe_pressure_energy_candidate_globalization(1U);
#endif
  const bool terminal_path_active = static_cast<bool>(status);
  if (solve_epoch.solve_calls() != 0U) {
    const Status observed = solve_epoch.observe(report);
    if (status && !observed) status = observed;
  }
  PressureEnergyCandidateLoopResult candidate_loop_one;
  if (status && pressure_energy_candidate_scope) {
    // The exact certificate binds these final trial revisions as its base.
    // No state-layer revision is permitted between this point and commit.
    status = revise_coupled_trial_views();
    if (status)
      status = run_pressure_energy_candidate_loop(
          1U, 0U, jacobian_one, intermediate_one, pressure_one,
          pressure_energy_target_flux,
          1.0,
          false,
          candidate_loop_one);
  } else if (status) {
    status = revise_exact_thermodynamic_workspaces();
    if (status) status = preflight_pressure_energy_candidate();
    if (status) status = revise_coupled_trial_views();
  }
  if (status) attempt_stage = 45U;
  if (status) {
    RevisionToken corrected_flux_revision = detail::product_mix(
        pressure_one.state,
        transaction.trial_revision(product.fields.enthalpy));
    corrected_flux_revision = detail::product_mix(
        corrected_flux_revision,
        transaction.trial_revision(product.fields.temperature));
    corrected_flux_revision = detail::product_mix(
        corrected_flux_revision,
        transaction.trial_revision(product.fields.rho));
    corrected_flux_revision = detail::product_mix(
        corrected_flux_revision,
        pressure_energy_candidate_density.revision);
    if (corrected_flux_revision == 0U) corrected_flux_revision = 1U;
    if (corrected_flux_revision == provisional_flux.revision) {
      corrected_flux_revision =
          corrected_flux_revision ==
                  std::numeric_limits<RevisionToken>::max()
              ? RevisionToken{1U}
              : corrected_flux_revision + 1U;
    }
    status = product.phi_workspace.workspace_view(
        1U, corrected_flux_revision, provisional_flux);
  }
  PisoStateCorrectionCertificate corrected_one;
  PisoExactThermodynamicCandidateView candidate_one;
  if (status) attempt_stage = 46U;
  if (status && !pressure_energy_candidate_scope)
    status = prepare_exact_thermodynamic_candidate(
        pressure_one, as_const(pressure_correction),
        as_const(enthalpy_correction),
        as_const(pressure_energy_candidate_enthalpy),
        as_const(pressure_energy_candidate_density),
        as_const(pressure_energy_candidate_temperature),
        as_const(compressibility),
        {species_accepted.data(), species_accepted.size()},
        {species_accepted.data(), species_accepted.size()}, candidate_one);
  if (status && pressure_energy_candidate_scope) {
    const bool publishable =
        candidate_loop_one.replay_valid &&
        candidate_loop_one.replay.exact_certificate.valid() &&
        candidate_reference_transition_valid(candidate_loop_one.replay);
    status = product.reductions.consensus(
        publishable
            ? Status{}
            : Status{StatusCode::invalid_plan,
                     kProductPressureEnergy});
    if (status)
      status = prepare_pressure_energy_candidate_material_publication();
    // This commit owns the final collective validation and then crosses its
    // no-fail state/flux copy boundary.  Everything below a successful call
    // is therefore scalar assignment or a prevalidated direct material copy.
    if (status && candidate_loop_one.stationary.valid())
      status = product.coupler.commit_frozen_momentum_stationary_trial_state(
          candidate_loop_one.authority,
          candidate_loop_one.replay.exact_certificate,
          candidate_loop_one.stationary,
          {trial_velocity, trial_pressure, trial_enthalpy, trial_density,
           trial_temperature},
          provisional_flux, product.reductions, corrected_one);
    else if (status)
      status = product.coupler.commit_frozen_momentum_coupled_trial_state(
          candidate_loop_one.authority,
          candidate_loop_one.replay.exact_certificate,
          candidate_loop_one.selection,
          {trial_velocity, trial_pressure, trial_enthalpy, trial_density,
           trial_temperature},
          provisional_flux, product.reductions, corrected_one);
    if (status) {
      publish_pressure_energy_candidate_materials();
      if (pressure_reference_kind == PressureReferenceKind::closed_mass) {
        const ClosedGaugeCorrectionCertificate& gauge =
            candidate_loop_one.replay.exact.closed_gauge;
        attempt_pressure_reference = gauge.next_pressure_reference;
        pressure_reference_certificate = gauge.output_pressure_reference;
      }
      intermediate_input.pressure_reference = pressure_reference_certificate;
      pressure_input.pressure_reference = pressure_reference_certificate;
    }
  } else if (status) {
    status = product.coupler.correct_coupled_trial_state(
        pressure_one, pressure_correction, as_const(enthalpy_correction),
        {trial_velocity, trial_pressure, trial_enthalpy, trial_density,
         trial_temperature},
        candidate_one, provisional_flux, product.reductions, corrected_one);
  }
  if (status && !pressure_energy_candidate_scope) {
    const PisoExactThermodynamicCandidateView& committed_candidate =
        candidate_one;
    const bool reference_matches =
        corrected_one.valid() &&
        same_pressure_reference(corrected_one.input_pressure_reference,
                                pressure_reference_certificate) &&
        (pressure_reference_kind != PressureReferenceKind::closed_mass ||
         (committed_candidate.closed_gauge.valid() &&
          same_pressure_reference(
              corrected_one.output_pressure_reference,
              committed_candidate.closed_gauge.output_pressure_reference)));
    if (!reference_matches) {
      status = {StatusCode::invalid_plan, kProductBinding};
    } else {
      if (pressure_reference_kind == PressureReferenceKind::closed_mass)
        attempt_pressure_reference =
            committed_candidate.closed_gauge.next_pressure_reference;
      pressure_reference_certificate =
          corrected_one.output_pressure_reference;
      intermediate_input.pressure_reference = pressure_reference_certificate;
      pressure_input.pressure_reference = pressure_reference_certificate;
    }
  }
  if (!pressure_energy_candidate_scope)
    status = product.reductions.consensus(status);

  // Corrector one changes rho. Integral and conditional boundary targets
  // must therefore be refreshed before corrector two forms its HbyA/phiHbyA
  // state; the storage is preallocated and the common velocity/pressure
  // inlet path performs no extra collective.
  if (status) {
    begin_timed_stage(45U);
    attempt_stage = 45U;
  }
  if (status && product.contributions.plan() == nullptr)
    status = {StatusCode::invalid_plan, kProductBinding};
  if (terminal_path_active)
    status = refresh_coupled_state(
        kCoupledStateC2Stage,
        BoundaryThermophysicalGhostPhase::corrector_two, status);
  if (status && product.piso.coupling() == CouplingKind::simple)
    status = publish_momentum_state(45U, status, true);
  if (status && product.piso.coupling() == CouplingKind::simple)
    status = run_momentum_predictor(status, false, 46U, 2U);
  if (status && product.piso.coupling() == CouplingKind::simple)
    status = refresh_coupled_state(
        kCoupledStateC2Stage,
        BoundaryThermophysicalGhostPhase::corrector_two, status);
  if (status && product.piso.coupling() == CouplingKind::simple) {
    intermediate_input.momentum = momentum_certificate;
    intermediate_input.momentum_system = momentum_system;
  }

  intermediate_input.corrector = 2U;
  intermediate_input.temporal_reference = {};
  intermediate_input.committed_face_history = {};
  if (status) {
    begin_timed_stage(50U);
    attempt_stage = 50U;
  }
  intermediate_input.prior_corrector = corrected_one.state;
  intermediate_input.density = trial_density;
  intermediate_input.trial_velocity = as_const(trial_velocity);
  intermediate_input.trial_flux = as_const(provisional_flux);
  intermediate_input.thermophysical_boundary = {
      boundary_thermo_certificate,
      {attempt_pressure_reference, as_const(trial_pressure),
       as_const(trial_enthalpy),
       {species_accepted.data(), species_accepted.size()},
       as_const(trial_density)}};
  PisoIntermediateCertificate intermediate_two;
  if (status) attempt_stage = 51U;
  if (terminal_path_active)
    status = product.coupler.refresh(intermediate_input, intermediate_two,
                                     status);
  const bool probe_stationary_two =
      pressure_energy_candidate_scope &&
      product.piso.coupling() == CouplingKind::simple &&
      pressure_energy_components_converged(candidate_loop_one);
  if (status && !probe_stationary_two)
    status = revise_pressure_energy_workspaces();
  if (status)
    status = probe_stationary_two
                 ? inspect_pressure_energy_target_flux(intermediate_two)
                 : assemble_pressure_energy_residual(intermediate_two);
  // C2 is a new numeric authority over the same preallocated storage.  Its
  // system and solution views must therefore carry fresh revisions rather
  // than silently overwriting the C1 certificate's fields.
  if (status) status = revise_pressure_system_workspaces();
  pressure_input.intermediate = intermediate_two;
  pressure_input.density_trial = as_const(trial_density);
  PressureCorrectionCertificate pressure_two;
  if (status) attempt_stage = 52U;
  if (status)
    status = product.coupler.assemble_pressure_system(
        pressure_input, pressure_system, pressure_two);
  LinearIdentity identity_two = linear_identity(pressure_two);
  MgCoefficientIdentity coefficient_two{pressure_two.state,
                                        pressure_two.state, 0.0};
  PressureEnergyCandidateLoopResult candidate_loop_two;
  PressureEnergyJacobianObservation jacobian_two;
  bool stationary_two = false;
  if (status && probe_stationary_two) {
    zero_field(pressure_correction);
    zero_field(enthalpy_correction);
    status = revise_coupled_trial_views();
    if (status)
      status = run_pressure_energy_candidate_loop(
          2U, 0U, jacobian_two, intermediate_two, pressure_two,
          pressure_energy_target_flux, 1.0, true, candidate_loop_two);
    stationary_two =
        static_cast<bool>(status) && candidate_loop_two.stationary.valid();
    if (stationary_two)
      status = solve_epoch.record_stationary(
          product.piso, pressure_two, candidate_loop_two.stationary,
          product.coupler);
    else if (status) {
      // The exact physical-boundary baseline did not remain terminal after
      // the fresh SIMPLE momentum sweep. Candidate evaluation revised cold
      // scratch and live destination identities, so rebuild the original C2
      // system before taking its ordinary coupled direction.
      intermediate_input.density = trial_density;
      intermediate_input.trial_velocity = as_const(trial_velocity);
      intermediate_input.trial_flux = as_const(provisional_flux);
      intermediate_input.thermophysical_boundary = {
          boundary_thermo_certificate,
          {attempt_pressure_reference, as_const(trial_pressure),
           as_const(trial_enthalpy),
           {species_accepted.data(), species_accepted.size()},
           as_const(trial_density)}};
      status = product.coupler.refresh(intermediate_input, intermediate_two);
      if (status) status = revise_pressure_energy_workspaces();
      if (status)
        status = assemble_pressure_energy_residual(intermediate_two);
      if (status) status = revise_pressure_system_workspaces();
      pressure_input.intermediate = intermediate_two;
      pressure_input.density_trial = as_const(trial_density);
      if (status)
        status = product.coupler.assemble_pressure_system(
            pressure_input, pressure_system, pressure_two);
      identity_two = linear_identity(pressure_two);
      coefficient_two = {pressure_two.state, pressure_two.state, 0.0};
    }
  }
  if (status) attempt_stage = 53U;
  if (status && !stationary_two && ibm_pressure_operator.has_value())
    status = ibm_pressure_operator->mask_solid_rhs(pressure_system.rhs);
  if (status && !stationary_two)
    status = form_pressure_energy_blocks();
  if (status && !stationary_two)
    status = solve_pressure_energy(
        2U, 0U, intermediate_two, pressure_two, pressure_energy_target_flux,
        pressure_energy_frozen_enthalpy, identity_two, coefficient_two,
        // Both mandatory PISO correctors are inexact nonlinear directions.
        // Exact C2 replay below either accepts the five gates or enters a
        // residual-aware refinement; it never accepts this linear norm alone.
        pressure_inexact_control(nullptr, 0.0, false), jacobian_two);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (status) observe_pressure_energy_candidate_globalization(2U);
#endif
  if (solve_epoch.solve_calls() != 0U) {
    const Status observed = solve_epoch.observe(report);
    if (status && !observed) status = observed;
  }
  if (status && pressure_energy_candidate_scope && !stationary_two) {
    // Freeze the final C2 trial identities before the baseline is certified;
    // none of these five views may be revised before the exact pending commit.
    status = revise_coupled_trial_views();
    if (status)
      status = run_pressure_energy_candidate_loop(
          2U, 0U, jacobian_two, intermediate_two, pressure_two,
          pressure_energy_target_flux,
          1.0,
          false,
          candidate_loop_two);
  } else if (status) {
    status = revise_exact_thermodynamic_workspaces();
    if (status) status = preflight_pressure_energy_candidate();
  }

  // C2's Schur system is deliberately a frozen-spatial quasi-Newton
  // linearization.  The exact nonlinear replay above may therefore select a
  // joint descent state which is not yet inside the independent continuity
  // and energy terminal gates.  Re-enter C2 only through the typed
  // provisional-state authority: atomically publish the exact p/h/rho/T/U
  // state and final physical mass flux, refresh every state-dependent
  // closure, then form a new direction at the same target time.  These are
  // nonlinear refinements of corrector two, not extra PISO correctors.
  bool pressure_energy_refinement_converged =
      !pressure_energy_candidate_scope ||
      pressure_energy_components_converged(candidate_loop_two);
  double pressure_energy_previous_merit = 0.0;
  bool pressure_energy_previous_merit_available =
      pressure_energy_loop_merit(candidate_loop_one,
                                 pressure_energy_previous_merit);
  std::uint8_t refinement_iteration = 1U;
  while (status && pressure_energy_candidate_scope &&
         !pressure_energy_refinement_converged &&
         refinement_iteration <= kPressureEnergyRefinementCapacity) {
    const bool publishable =
        candidate_loop_two.replay_valid &&
        candidate_loop_two.replay.exact_certificate.valid() &&
        candidate_loop_two.selection.valid() &&
        candidate_reference_transition_valid(candidate_loop_two.replay);
    status = product.reductions.consensus(
        publishable
            ? Status{}
            : Status{StatusCode::invalid_plan, kProductPressureEnergy});
    if (status)
      status = prepare_pressure_energy_candidate_material_publication();

    // Replica one is the live provisional mass-flux carrier.  Give every
    // refinement publication a fresh revision before the no-fail atomic
    // copy; the selected candidate remains sealed in the candidate replica.
    RevisionToken refinement_flux_revision = detail::product_mix(
        candidate_loop_two.replay.sample.mass_flux_provenance,
        refinement_iteration);
    refinement_flux_revision = detail::product_mix(
        refinement_flux_revision,
        transaction.trial_revision(product.fields.rho));
    refinement_flux_revision = detail::product_mix(
        refinement_flux_revision,
        transaction.trial_revision(product.fields.enthalpy));
    if (refinement_flux_revision == 0U) refinement_flux_revision = 1U;
    if (refinement_flux_revision == provisional_flux.revision) {
      refinement_flux_revision =
          refinement_flux_revision ==
                  std::numeric_limits<RevisionToken>::max()
              ? RevisionToken{1U}
              : refinement_flux_revision + 1U;
    }
    if (status)
      status = product.phi_workspace.workspace_view(
          1U, refinement_flux_revision, provisional_flux);

    PisoStateCorrectionCertificate refinement_correction;
    PisoPressureEnergyRefinementStateCertificate refinement_state;
    if (status)
      status = product.coupler
                   .commit_frozen_momentum_coupled_refinement_trial_state(
                       candidate_loop_two.authority,
                       candidate_loop_two.replay.exact_certificate,
                       candidate_loop_two.selection,
                       {trial_velocity, trial_pressure, trial_enthalpy,
                        trial_density, trial_temperature},
                       provisional_flux, refinement_iteration,
                       product.reductions, refinement_correction,
                       refinement_state);
    if (status) {
      publish_pressure_energy_candidate_materials();
      if (pressure_reference_kind == PressureReferenceKind::closed_mass) {
        const ClosedGaugeCorrectionCertificate& gauge =
            candidate_loop_two.replay.exact.closed_gauge;
        attempt_pressure_reference = gauge.next_pressure_reference;
        pressure_reference_certificate = gauge.output_pressure_reference;
      }
      intermediate_input.pressure_reference =
          pressure_reference_certificate;
      pressure_input.pressure_reference = pressure_reference_certificate;
    }

    // Halo, physical-boundary, IBM donor, material and turbulence state must
    // all correspond to the just-published exact state before the typed
    // authority is consumed.  The authority seals interior primitive bytes;
    // this refresh changes only certified ghost/donor closures.
    if (status)
      status = refresh_coupled_state(
          kCoupledStateC2Stage,
          BoundaryThermophysicalGhostPhase::corrector_two, status);
    intermediate_input.prior_corrector = refinement_correction.state;
    intermediate_input.density = trial_density;
    intermediate_input.trial_velocity = as_const(trial_velocity);
    intermediate_input.trial_flux = as_const(provisional_flux);
    intermediate_input.thermophysical_boundary = {
        boundary_thermo_certificate,
        {attempt_pressure_reference, as_const(trial_pressure),
         as_const(trial_enthalpy),
         {species_accepted.data(), species_accepted.size()},
         as_const(trial_density)}};
    PisoIntermediateCertificate refined_intermediate;
    // Once issued, the authority must be consumed even if the intervening
    // ghost/material refresh fails.  Passing that failure as prerequisite
    // keeps every rank on the same collective fail-close path and prevents a
    // stale provisional-state capability from surviving the rejected round.
    if (refinement_state.valid())
      status = product.coupler.refresh_pressure_energy_refinement(
          refinement_state, intermediate_input,
          {trial_velocity, trial_pressure, trial_enthalpy, trial_density,
           trial_temperature},
          refined_intermediate, status);
    if (status) status = revise_pressure_energy_workspaces();
    if (status)
      status = assemble_pressure_energy_residual(refined_intermediate);
    if (status) status = revise_pressure_system_workspaces();
    pressure_input.intermediate = refined_intermediate;
    pressure_input.density_trial = as_const(trial_density);
    PressureCorrectionCertificate refined_pressure;
    if (status)
      status = product.coupler.assemble_pressure_system(
          pressure_input, pressure_system, refined_pressure);
    const LinearIdentity refinement_identity =
        linear_identity(refined_pressure);
    const MgCoefficientIdentity refinement_coefficient{
        refined_pressure.state, refined_pressure.state, 0.0};
    double pressure_energy_current_merit = 0.0;
    const bool pressure_energy_current_merit_available =
        pressure_energy_loop_merit(candidate_loop_two,
                                   pressure_energy_current_merit);
    const double refinement_extrapolated_alpha =
        pressure_energy_previous_merit_available &&
                pressure_energy_current_merit_available
            ? detail::product_pressure_aitken_initial_alpha(
                  pressure_energy_previous_merit,
                  pressure_energy_current_merit,
                  candidate_loop_two.selected_alpha)
            : 1.0;
    const LinearSolveControl refinement_solve_control =
        pressure_inexact_control(
            &candidate_loop_two, pressure_energy_previous_merit,
            pressure_energy_previous_merit_available);
    if (status && ibm_pressure_operator.has_value())
      status = ibm_pressure_operator->mask_solid_rhs(pressure_system.rhs);
    if (status) status = form_pressure_energy_blocks();
    PressureEnergyJacobianObservation refinement_jacobian;
    if (status)
      status = solve_pressure_energy(
          2U, refinement_iteration, refined_intermediate, refined_pressure,
          pressure_energy_target_flux, pressure_energy_frozen_enthalpy,
          refinement_identity, refinement_coefficient,
          refinement_solve_control, refinement_jacobian);
    if (solve_epoch.solve_calls() != 0U) {
      const Status observed = solve_epoch.observe(report);
      if (status && !observed) status = observed;
    }

    if (status) {
      intermediate_two = refined_intermediate;
      pressure_two = refined_pressure;
      // Freeze fresh destination revisions before certifying the next exact
      // baseline and its selected state/flux replay.
      status = revise_coupled_trial_views();
    }
    if (status)
      status = run_pressure_energy_candidate_loop(
          2U, refinement_iteration, refinement_jacobian, intermediate_two,
          pressure_two,
          pressure_energy_target_flux,
          refinement_extrapolated_alpha,
          false,
          candidate_loop_two);
    if (status)
      pressure_energy_refinement_converged =
          pressure_energy_components_converged(candidate_loop_two);
    pressure_energy_previous_merit = pressure_energy_current_merit;
    pressure_energy_previous_merit_available =
        pressure_energy_current_merit_available;
    ++refinement_iteration;
  }
  if (pressure_energy_candidate_scope) {
    if (!status) {
      report.pressure_energy_refinement_termination =
          PressureEnergyRefinementTermination::rejected_candidate;
    } else if (pressure_energy_refinement_converged) {
      report.pressure_energy_refinement_termination =
          PressureEnergyRefinementTermination::component_residuals_converged;
    } else {
      report.pressure_energy_refinement_termination =
          PressureEnergyRefinementTermination::iteration_capacity_exhausted;
      // The resource bound is itself an acceptance boundary.  Do not rely on
      // the later terminal audit to happen to reject the last replay: a
      // refreshed replay and the final assembly can straddle a component
      // tolerance by roundoff.  Publish the complete observed prefix, then
      // reject collectively before pending flux or terminal-audit authority
      // exists.
      attempt_stage = 54U;
      status = product.reductions.consensus(
          {StatusCode::rejected_step, kProductPressureEnergy});
    }
  }
  if (status) attempt_stage = 54U;
  if (status)
    status = final_flux_writer.begin_pending(transaction, product.final_flux,
                                             pending_flux);
  if (status && !pressure_energy_candidate_scope)
    status = revise_coupled_trial_views();
  PisoStateCorrectionCertificate corrected_two;
  PisoExactThermodynamicCandidateView candidate_two;
  if (status) attempt_stage = 55U;
  if (status && !pressure_energy_candidate_scope)
    status = prepare_exact_thermodynamic_candidate(
        pressure_two, as_const(pressure_correction),
        as_const(enthalpy_correction),
        as_const(pressure_energy_candidate_enthalpy),
        as_const(pressure_energy_candidate_density),
        as_const(pressure_energy_candidate_temperature),
        as_const(compressibility),
        {species_accepted.data(), species_accepted.size()},
        {species_accepted.data(), species_accepted.size()}, candidate_two);
  if (status && pressure_energy_candidate_scope) {
    const bool publishable =
        candidate_loop_two.replay_valid &&
        candidate_loop_two.replay.exact_certificate.valid() &&
        candidate_reference_transition_valid(candidate_loop_two.replay);
    status = product.reductions.consensus(
        publishable
            ? Status{}
            : Status{StatusCode::invalid_plan,
                     kProductPressureEnergy});
    if (status)
      status = prepare_pressure_energy_candidate_material_publication();
    // As in C1, this is the final fallible operation for the selected C2
    // state.  A successful return has already crossed the coupler's internal
    // collective gate and copied all five fields plus the pending flux.
    if (status && candidate_loop_two.stationary.valid())
      status =
          product.coupler.commit_frozen_momentum_stationary_pending_state(
              candidate_loop_two.authority,
              candidate_loop_two.replay.exact_certificate,
              candidate_loop_two.stationary,
              {trial_velocity, trial_pressure, trial_enthalpy, trial_density,
               trial_temperature},
              pending_flux, product.reductions, corrected_two);
    else if (status)
      status = product.coupler.commit_frozen_momentum_coupled_pending_state(
          candidate_loop_two.authority,
          candidate_loop_two.replay.exact_certificate,
          candidate_loop_two.selection,
          {trial_velocity, trial_pressure, trial_enthalpy, trial_density,
           trial_temperature},
          pending_flux, product.reductions, corrected_two);
    if (status) {
      publish_pressure_energy_candidate_materials();
      if (pressure_reference_kind == PressureReferenceKind::closed_mass) {
        const ClosedGaugeCorrectionCertificate& gauge =
            candidate_loop_two.replay.exact.closed_gauge;
        attempt_pressure_reference = gauge.next_pressure_reference;
        pressure_reference_certificate = gauge.output_pressure_reference;
      }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
      if (g_candidate_globalization_published.load(
              std::memory_order_acquire)) {
        detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic =
            g_candidate_globalization_diagnostic;
        diagnostic.committed = true;
        publish_candidate_globalization_diagnostic(diagnostic);
      }
#endif
    }
  } else if (status) {
    status = product.coupler.correct_coupled_pending_state(
        pressure_two, pressure_correction, as_const(enthalpy_correction),
        {trial_velocity, trial_pressure, trial_enthalpy, trial_density,
         trial_temperature},
        candidate_two, pending_flux, product.reductions, corrected_two);
  }
  if (status && !pressure_energy_candidate_scope) {
    const bool reference_matches =
        corrected_two.valid() &&
        same_pressure_reference(corrected_two.input_pressure_reference,
                                pressure_reference_certificate) &&
        (pressure_reference_kind != PressureReferenceKind::closed_mass ||
         (candidate_two.closed_gauge.valid() &&
          same_pressure_reference(
              corrected_two.output_pressure_reference,
              candidate_two.closed_gauge.output_pressure_reference)));
    if (!reference_matches) {
      status = {StatusCode::invalid_plan, kProductBinding};
    } else {
      if (pressure_reference_kind == PressureReferenceKind::closed_mass)
        attempt_pressure_reference =
            candidate_two.closed_gauge.next_pressure_reference;
      pressure_reference_certificate =
          corrected_two.output_pressure_reference;
    }
  }
  if (!pressure_energy_candidate_scope)
    status = product.reductions.consensus(status);

  // Re-evaluate the complete conservative energy equation against the exact
  // EOS state and the actual C2 pending flux.  This is the nonlinear closure
  // gate for the quasi-Newton Schur certificate; it is intentionally formed
  // before audit_pending_final so energy participates in the same terminal
  // accept/reject decision as continuity and EOS.
  double local_normalized_energy_residual = 0.0;
  ConstFaceFluxView terminal_energy_flux;
  if (status)
    status = refresh_coupled_state(
        kCoupledStateTerminalStage,
        BoundaryThermophysicalGhostPhase::terminal, status);
  if (status) status = revise_pressure_energy_workspaces();
  if (status)
    status = product.coupler.inspect_corrected_pending(
        corrected_two, pending_flux, terminal_energy_flux);
  status = product.reductions.consensus(status);
  if (status)
    status = assemble_pressure_energy_residual_from_flux(
        terminal_energy_flux,
        EquationAssemblyScope::final_conservative);
  if (status) {
    Status local_energy_status;
    std::size_t offset = 0U;
    for (std::int32_t z = 0; z < cells.z && local_energy_status; ++z) {
      for (std::int32_t y = 0; y < cells.y && local_energy_status; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x, ++offset) {
          if (!pressure_energy_active(offset)) continue;
          const Int3 cell{x, y, z};
          const double volume =
              detail::cell_volume(product.equations.kernels(), cell);
          const double rho = trial_density.unchecked(cell, 0U);
          const double h = trial_enthalpy.unchecked(cell, 0U);
          const double temperature =
              trial_temperature.unchecked(cell, 0U);
          const double cp = heat_capacity.unchecked(cell, 0U);
          const double p_abs = attempt_pressure_reference +
                               trial_pressure.unchecked(cell, 0U);
          double scale = 0.0;
          const bool scale_valid = pressure_energy_temporal_scale(
              attempt_pressure_reference, as_const(trial_density),
              as_const(trial_enthalpy), as_const(trial_pressure), cell,
              volume, scale);
          const double residual =
              pressure_energy_r_e.unchecked(cell, 0U);
          const double metric = std::abs(residual) / scale;
          if (!std::isfinite(volume) || !(volume > 0.0) ||
              !std::isfinite(rho) || !(rho > 0.0) ||
              !std::isfinite(h) || !std::isfinite(temperature) ||
              !(temperature > 0.0) || !std::isfinite(cp) ||
              !(cp > 0.0) || !std::isfinite(p_abs) || !(p_abs > 0.0) ||
              !scale_valid ||
              !std::isfinite(residual) || !std::isfinite(metric)) {
            local_energy_status = {StatusCode::rejected_step,
                                   kProductPressureEnergy};
            break;
          }
          local_normalized_energy_residual =
              std::max(local_normalized_energy_residual, metric);
        }
      }
    }
    status = product.reductions.consensus(local_energy_status);
  }

  double local_boundary_closure_residual = 0.0;
  std::uint64_t local_boundary_closure_samples = 0U;
  if (status && pressure_reference_certificate.kind !=
                    PressureReferenceKind::closed_mass) {
    status = apply_boundary_ghosts(BoundaryStage::pressure, product.boundary,
                                   {&trial_pressure, 1U}, boundary_values);
    if (status)
      status = local_pressure_outlet_closure(
          product.boundary, product.boundary_specs, as_const(trial_pressure),
          attempt_pressure_reference, local_boundary_closure_residual,
          local_boundary_closure_samples);
  }

  if (status) status = runtime_write_view(product.fields.eos_density, eos_density);
  if (status) {
    // At frozen h/Y, ideal-gas rho is affine in local absolute pressure;
    // re-evaluation remains the independent terminal EOS oracle.
    for (std::int32_t z = 0; z < cells.z && status; ++z) {
      for (std::int32_t y = 0; y < cells.y && status; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          for (std::size_t species = 0U; species < species_trial.size();
               ++species)
            species_values[species] =
                species_trial[species].unchecked(cell, 0U);
          ThermoState thermo;
          status = product.thermodynamics.evaluate_from_reference_pressure(
              attempt_pressure_reference,
              trial_pressure.unchecked(cell, 0U),
              trial_enthalpy.unchecked(cell, 0U),
              {species_values.data(), species_values.size()},
              {trial_velocity.unchecked(cell, 0U),
               trial_velocity.unchecked(cell, 1U),
               trial_velocity.unchecked(cell, 2U)},
              thermo, trial_temperature.unchecked(cell, 0U));
          if (!status) break;
          eos_density.unchecked(cell, 0U) = thermo.rho;
        }
      }
    }
  }
  PisoTerminalAuditInput audit;
  if (status) {
    begin_timed_stage(60U);
    attempt_stage = 60U;
  }
  audit.correction = corrected_two;
  audit.pressure_reference = pressure_reference_certificate;
  audit.density = as_const(trial_density);
  audit.eos_density = as_const(eos_density);
  audit.density_accepted = rho_history.accepted;
  audit.density_previous =
      effective_bdf.order == 2U ? rho_history.previous : ConstFieldView{};
  audit.pressure_perturbation = as_const(trial_pressure);
  audit.drho_dp_h_y =
      pressure_energy_candidate_scope
          ? candidate_loop_two.replay.exact.pressure_compressibility
          : as_const(compressibility);
  audit.bdf = effective_bdf;
  audit.step_dt = step.dt;
  audit.convective_cfl_limit = product.time.spec().convective_cfl;
  audit.closed_mass_target = closed_mass_target;
  audit.boundary_closure_residual = local_boundary_closure_residual;
  audit.boundary_closure_samples = local_boundary_closure_samples;
  audit.energy_residual = local_normalized_energy_residual;
  // Terminal audit must consume the exact activity storage certified by the
  // C2 correction, not an equal-valued copy owned by the topology plan.
  audit.active = pressure_energy_activity.cells;
  audit.thermophysical_boundary = {
      boundary_thermo_certificate,
      {attempt_pressure_reference, as_const(trial_pressure),
       as_const(trial_enthalpy),
       {species_accepted.data(), species_accepted.size()},
       as_const(trial_density)}};
  PisoTerminalCertificate terminal;
  FinalForceCertificate force_certificate;
  if (terminal_path_active)
    status = product.coupler.audit_pending_final(
        audit, pending_flux, product.reductions, report, terminal, status);
  // The terminal audit has reconstructed the committed value from the exact
  // C2 density and final mass flux.  Apply the ceiling before state/flux
  // publication: fixed scientific cases fail fatally instead of silently
  // changing dt=0.006, while adaptive cases reject this attempt so the time
  // controller can retry at its recorded smaller step.
  if (status) {
    status = convective_cfl_acceptance_status(
        product.time.spec().control,
        report.committed_convective_cfl_out_max,
        report.committed_convective_cfl_limit);
  }
  const RevisionSourceId pressure_reference_revision_source =
      static_cast<RevisionSourceId>(product.layers.field_count() + 1U);
  if (status &&
      (!terminal.valid() ||
       !same_pressure_reference(terminal.pressure_reference,
                                pressure_reference_certificate)))
    status = {StatusCode::invalid_plan, kProductBinding};
  if (status)
    status = transaction.bind_dependency(
        {pressure_reference_revision_source,
         terminal.pressure_reference.pressure_reference});
  status = product.reductions.consensus(status);

  // Final-state rate histories are committed with the state and feed the
  // next predictor. They do not consume or publish another face flux.
  if (status) attempt_stage = 61U;
  if (status) halo_views[0U] = trial_velocity;
  if (status) status = exchange(product.final_velocity_halo, 61U, 1U);
  if (status) trial_velocity = halo_views[0U];
  if (status)
    status = apply_boundary_ghosts(BoundaryStage::momentum, product.boundary,
                                   {&trial_velocity, 1U}, boundary_values);
  status = product.reductions.consensus(status);
  if (status && product.ibm_gradient_donors.has_value()) {
    std::array<FieldView, 1U> donor_fields{trial_velocity};
    status = product.ibm_gradient_donors->exchange(
        kIbmGradientDonorStage, {donor_fields.data(), donor_fields.size()});
    trial_velocity = donor_fields[0U];
  }
  status = product.reductions.consensus(status);
  const bool final_rate_path_active = static_cast<bool>(status);
  if (status) status = runtime_write_view(product.fields.velocity_gradient,
                                          velocity_gradient);
  if (status) {
    const std::array<ConstFieldView, 1U> reads{as_const(trial_velocity)};
    const std::array<FieldView, 1U> writes{velocity_gradient};
    status = cartesian_gradient(
        product.equations.kernels(),
        {{reads.data(), reads.size()}, {writes.data(), writes.size()},
         full_box, 0U, 0U, 3U, 0U, nullptr});
  }
  if (status && product.ibm_equations.has_value())
    status = product.ibm_equations->correct_velocity_gradient(
        as_const(trial_velocity), velocity_gradient);
  if (status)
    status = runtime_write_view(product.fields.effective_viscosity,
                                effective_viscosity);
  if (status) {
    const TurbulenceUpdateInput turbulence_input{
        as_const(trial_density), as_const(molecular_viscosity), {},
        velocity_gradient.revision, as_const(velocity_gradient)};
    status = product.turbulence.update(turbulence_input, effective_viscosity,
                                       turbulence_certificate);
  }
  if (final_rate_path_active && product.transport.kernel() ==
                                    TransportKernel::coast_native_air)
    status = refresh_live_effective_thermal_ghosts(60U, status);
  else
    status = product.reductions.consensus(status);
  const std::size_t final_rate_halo_count =
      7U + product.fields.scalars.size();
  if (status) {
    halo_count = 0U;
    halo_views[halo_count++] = trial_pressure;
    halo_views[halo_count++] = trial_enthalpy;
    halo_views[halo_count++] = trial_temperature;
    halo_views[halo_count++] = velocity_gradient;
    halo_views[halo_count++] = molecular_viscosity;
    halo_views[halo_count++] = effective_viscosity;
    halo_views[halo_count++] = conductivity;
    for (FieldView value : species_trial) halo_views[halo_count++] = value;
    for (FieldView value : passive_trial) halo_views[halo_count++] = value;
  }
  if (final_rate_path_active)
    status = exchange(product.stage_halos[5U], 60U,
                      final_rate_halo_count, status);
  if (status) {
    trial_pressure = halo_views[0U];
    trial_enthalpy = halo_views[1U];
    trial_temperature = halo_views[2U];
    velocity_gradient = halo_views[3U];
    molecular_viscosity = halo_views[4U];
    effective_viscosity = halo_views[5U];
    conductivity = halo_views[6U];
    species_index = 0U;
    passive_index = 0U;
    for (std::size_t index = 0U; index < product.fields.scalars.size();
         ++index) {
      if (product.fields.scalar_roles[index] ==
          TransportedScalarRole::species)
        species_trial[species_index++] = halo_views[index + 7U];
      else
        passive_trial[passive_index++] = halo_views[index + 7U];
    }
    std::array<FieldView, 4U> derived{velocity_gradient,
                                     molecular_viscosity,
                                     effective_viscosity, conductivity};
    status = apply_physical_zero_gradient(
        product.boundary, {derived.data(), derived.size()});
  }
  if (status) {
    for (std::size_t species = 0U; species < species_trial.size(); ++species)
      species_accepted[species] = as_const(species_trial[species]);
    for (std::size_t passive = 0U; passive < passive_trial.size(); ++passive)
      passive_accepted[passive] = as_const(passive_trial[passive]);
    status = resolve_static_boundary_values(
        communicator, product.boundary, product.boundary_specs,
        product.geometry, product.patch, product.thermodynamics,
        attempt_pressure_reference, as_const(trial_density),
        as_const(trial_velocity), as_const(trial_enthalpy),
        {species_accepted.data(), species_accepted.size()},
        {passive_accepted.data(), passive_accepted.size()},
        {species_values.data(), species_values.size()}, boundary_scalars,
        boundary_vectors, boundary_normal_gradients, false, status);
  }
  status = product.reductions.consensus(status);
  if (status && !product.fields.scalars.empty()) {
    species_index = 0U;
    passive_index = 0U;
    for (std::size_t scalar = 0U; scalar < product.fields.scalars.size();
         ++scalar) {
      if (product.fields.scalar_roles[scalar] ==
          TransportedScalarRole::species)
        halo_views[scalar] = species_trial[species_index++];
      else
        halo_views[scalar] = passive_trial[passive_index++];
    }
    status = apply_boundary_ghosts(
        BoundaryStage::scalar, product.boundary,
        {halo_views.data(), product.fields.scalars.size()}, boundary_values);
    if (status) {
      species_index = 0U;
      passive_index = 0U;
      for (std::size_t scalar = 0U; scalar < product.fields.scalars.size();
           ++scalar) {
        if (product.fields.scalar_roles[scalar] ==
            TransportedScalarRole::species) {
          species_trial[species_index] = halo_views[scalar];
          species_accepted[species_index] =
              as_const(species_trial[species_index]);
          ++species_index;
        } else {
          passive_trial[passive_index] = halo_views[scalar];
          passive_accepted[passive_index] =
              as_const(passive_trial[passive_index]);
          ++passive_index;
        }
      }
    }
  }
  status = product.reductions.consensus(status);
  if (status && product.ibm_rate_donors.has_value()) {
    halo_count = 0U;
    halo_views[halo_count++] = trial_pressure;
    halo_views[halo_count++] = trial_temperature;
    for (FieldView value : species_trial) halo_views[halo_count++] = value;
    for (FieldView value : passive_trial) halo_views[halo_count++] = value;
    status = product.ibm_rate_donors->exchange(
        161U, {halo_views.data(), halo_count});
    if (status) {
      trial_pressure = halo_views[0U];
      trial_temperature = halo_views[1U];
      species_index = 0U;
      passive_index = 0U;
      for (std::size_t index = 0U; index < product.fields.scalars.size();
           ++index) {
        if (product.fields.scalar_roles[index] ==
            TransportedScalarRole::species)
          species_trial[species_index++] = halo_views[index + 2U];
        else
          passive_trial[passive_index++] = halo_views[index + 2U];
      }
    }
  }
  status = product.reductions.consensus(status);
  if (status && final_force_cache.has_value()) {
    if (product.ibm_force_donors.has_value()) {
      std::array<FieldView, 4U> donor_fields{
          trial_velocity, trial_pressure, velocity_gradient,
          effective_viscosity};
      status = product.ibm_force_donors->exchange(
          160U, {donor_fields.data(), donor_fields.size()});
      trial_velocity = donor_fields[0U];
      trial_pressure = donor_fields[1U];
      velocity_gradient = donor_fields[2U];
      effective_viscosity = donor_fields[3U];
    }
    status = product.reductions.consensus(status);
    ConstFaceFluxView inspected_final_flux;
    if (status)
      status = product.coupler.inspect_pending_final(
          terminal, pending_flux, inspected_final_flux);
    if (status)
      status = transaction.bind_dependency(
          {AttemptTransaction::field_revision_source(
               product.fields.velocity_gradient),
           velocity_gradient.revision});
    if (status)
      status = transaction.bind_dependency(
          {AttemptTransaction::field_revision_source(
               product.fields.effective_viscosity),
           effective_viscosity.revision});
    const std::array<RevisionDependency, 5U> force_dependencies{{
        {AttemptTransaction::field_revision_source(product.fields.velocity),
         transaction.trial_revision(product.fields.velocity)},
        {AttemptTransaction::field_revision_source(product.fields.pressure),
         transaction.trial_revision(product.fields.pressure)},
        {AttemptTransaction::field_revision_source(
             product.fields.velocity_gradient),
         velocity_gradient.revision},
        {AttemptTransaction::field_revision_source(
             product.fields.effective_viscosity),
         effective_viscosity.revision},
        {pressure_reference_revision_source,
         terminal.pressure_reference.pressure_reference}}};
    if (status) {
      FinalSurfaceState surface_state;
      surface_state.terminal_plan = terminal.plan;
      surface_state.terminal_state = terminal.audit_state;
      surface_state.final_flux = terminal.final_flux;
      surface_state.face_flux = inspected_final_flux;
      surface_state.final_velocity = as_const(trial_velocity);
      surface_state.pressure_perturbation = as_const(trial_pressure);
      surface_state.velocity_gradient = as_const(velocity_gradient);
      surface_state.effective_viscosity = as_const(effective_viscosity);
      surface_state.gradient_authority = make_derived_revision_tuple(
          surface_state.final_velocity, product.geometry, product.patch,
          product.boundary.revision(),
          product.boundary.semantic_fingerprint(),
          turbulence_certificate.state, turbulence_certificate.plan);
      surface_state.turbulence = turbulence_certificate;
      surface_state.geometry = product.geometry.topology_revision();
      surface_state.pressure_reference = attempt_pressure_reference;
      status = final_force_cache->prepare(
          surface_state,
          {force_dependencies.data(), force_dependencies.size()}, transaction,
          force_certificate);
      if (status && force_certificate.valid()) pending_force_cache = true;
    }
  }
  FieldView enthalpy_rate_output;
  FieldView rate_scratch;
  FieldView scalar_diffusivity;
  if (status)
    status = product.layers.view(StateRole::trial,
                                 product.fields.enthalpy_nonadvective_rate,
                                 enthalpy_rate_output);
  species_index = 0U;
  passive_index = 0U;
  for (std::size_t index = 0U;
       index < product.fields.scalar_nonadvective_rates.size() && status;
       ++index) {
    FieldView rate;
    status = product.layers.view(
        StateRole::trial, product.fields.scalar_nonadvective_rates[index],
        rate);
    if (!status) break;
    if (product.fields.scalar_roles[index] ==
        TransportedScalarRole::species)
      species_rate_output[species_index++] = rate;
    else
      passive_rate_output[passive_index++] = rate;
  }
  if (status)
    status = runtime_write_view(product.fields.predictor_accepted_advection,
                                rate_scratch);
  if (status)
    status = runtime_write_view(product.fields.scalar_diffusivity,
                                scalar_diffusivity);
  if (status) {
    material.molecular_viscosity = as_const(molecular_viscosity);
    material.effective_viscosity = as_const(effective_viscosity);
    material.thermal_conductivity = as_const(conductivity);
    equation_state.pressure_reference = attempt_pressure_reference;
    status = history(product.fields.rho, equation_state.density);
    if (status) status = history(product.fields.velocity, equation_state.velocity);
    if (status)
      status = history(product.fields.pressure,
                       equation_state.pressure_perturbation);
    if (status) status = history(product.fields.enthalpy, equation_state.enthalpy);
    if (status)
      status = history(product.fields.temperature, equation_state.temperature);
    equation_state.independent_species = {species_history.data(),
                                          species_history.size()};
    equation_state.passive_scalars = {passive_history.data(),
                                      passive_history.size()};
  }
  ThermophysicalRateCertificate rate_certificate;
  if (status) attempt_stage = 62U;
  if (status)
    status = evaluate_thermophysical_rates(
        product.equations,
        {equation_state, material, as_const(velocity_gradient), effective_bdf,
         1U,
         {}, product.ibm_equations.has_value() ? &*product.ibm_equations
                                               : nullptr},
        {enthalpy_rate_output,
         {species_rate_output.data(), species_rate_output.size()},
         {passive_rate_output.data(), passive_rate_output.size()},
         rate_scratch, scalar_diffusivity},
        rate_certificate);
  // The final conservative mass flux is published only against the complete
  // same-target thermodynamic state.  Omitting h/T would allow a stale
  // pressure-only correction certificate to publish after EOS closure had
  // changed the face-density authority.
  const std::size_t expected_final_dependencies =
      6U + species_trial.size();
  if (status &&
      final_dependencies.size() != expected_final_dependencies)
    status = {StatusCode::invalid_plan, kProductBinding};
  if (status) {
    final_dependencies[0U] = {
        AttemptTransaction::field_revision_source(product.fields.rho),
        transaction.trial_revision(product.fields.rho)};
    final_dependencies[1U] = {
        AttemptTransaction::field_revision_source(product.fields.velocity),
        transaction.trial_revision(product.fields.velocity)};
    final_dependencies[2U] = {
        AttemptTransaction::field_revision_source(product.fields.pressure),
        transaction.trial_revision(product.fields.pressure)};
    final_dependencies[3U] = {
        AttemptTransaction::field_revision_source(product.fields.enthalpy),
        transaction.trial_revision(product.fields.enthalpy)};
    final_dependencies[4U] = {
        AttemptTransaction::field_revision_source(product.fields.temperature),
        transaction.trial_revision(product.fields.temperature)};
    std::size_t dependency = 5U;
    for (std::size_t scalar = 0U; scalar < product.fields.scalars.size();
         ++scalar) {
      if (product.fields.scalar_roles[scalar] !=
          TransportedScalarRole::species)
        continue;
      const FieldId field = product.fields.scalars[scalar];
      final_dependencies[dependency++] = {
          AttemptTransaction::field_revision_source(field),
          transaction.trial_revision(field)};
    }
    if (dependency != expected_final_dependencies - 1U) {
      status = {StatusCode::invalid_plan, kProductBinding};
    } else {
      final_dependencies[dependency] = {
          pressure_reference_revision_source,
          terminal.pressure_reference.pressure_reference};
    }
  }
  status = product.reductions.consensus(status);
  if (status)
    attempt_stage = 63U;
  if (status)
    status = product.coupler.publish_pending_final(
        terminal,
        {final_dependencies.data(), final_dependencies.size()},
        {species_accepted.data(), species_accepted.size()},
        product.reductions,
        final_flux_writer, pending_flux);
  if (status) {
    begin_timed_stage(70U);
    attempt_stage = 70U;
  }
  // Finalize the solve report while the epoch is still an attempt-local,
  // stateless validator.  A failure here must participate in transaction
  // consensus rather than being discovered after state rotation.
  if (status) status = solve_epoch.finalize(report);
  if (status)
    pending_pressure_reference = {
        attempt_pressure_reference, terminal.pressure_reference,
        terminal.audit_state};
  const Status prepare_status =
      transaction.collective_prepare(communicator, status, prepared);
  if (!prepare_status)
    status = prepare_status;
  else if (prepared.decision() == AttemptFinishDecision::reject)
    status = prepared.outcome();
  else
    status = {};
  finish_stage_timings();
  return status;
}

Status ProductDriver::advance(LocalTimeLimits limits,
                              DriverStepReport& report) noexcept {
  if (implementation_ == nullptr || !implementation_->initialized ||
      implementation_->time.has_active_proposal()) {
    return {StatusCode::invalid_plan, kProductInput};
  }
  DriverStepReport candidate;
  std::uint64_t predictor_blocking_collectives = 0U;
  const DriverResourceReport resources_before =
      implementation_->resource_snapshot();
  StepTime proposal;
  Status status = implementation_->time.propose(
      implementation_->communicator, limits, proposal);
  Status last_attempt_status;
  while (status) {
    ++candidate.attempts;
    candidate.proposal = proposal;
    PisoAttemptReport attempt_report;
    PreparedAttemptFinish prepared_attempt;
    const Status attempt_status =
        implementation_->execute_attempt(proposal, attempt_report,
                                         prepared_attempt);
    implementation_->accumulate_stage_timings(candidate);
    last_attempt_status = attempt_status;
    candidate.effective_bdf = implementation_->effective_bdf;
    candidate.thermophysical_predictor_calls =
        implementation_->thermophysical_predictor_calls;
    candidate.temporal_method_fallback =
        implementation_->temporal_method_fallback;
    candidate.piso = attempt_report;
    candidate.thermophysical_predictor =
        implementation_->predictor_diagnostics;
    candidate.pressure_energy_globalization =
        implementation_->pressure_energy_globalization;
    candidate.momentum_predictor_limiter =
        implementation_->momentum_predictor_limiter;
    candidate.momentum_predictor_solve =
        implementation_->momentum_predictor_solve;
    const std::uint64_t predictor_calls =
        implementation_->predictor_diagnostics.blocking_collectives;
    predictor_blocking_collectives =
        predictor_calls > UINT64_MAX - predictor_blocking_collectives
            ? UINT64_MAX
            : predictor_blocking_collectives + predictor_calls;
    if (!attempt_status) {
      candidate.failed_stage = implementation_->attempt_stage;
      candidate.failure = attempt_status;
      if (!candidate.numerical_failure.valid &&
          implementation_->numerical_failure.valid) {
        candidate.numerical_failure = implementation_->numerical_failure;
      }
    }
    if (!prepared_attempt.valid()) {
      candidate.resources = resource_difference(
          implementation_->resource_snapshot(), resources_before);
      candidate.resources.predictor_blocking_collectives =
          predictor_blocking_collectives;
      report = candidate;
      return attempt_status;
    }
    const AttemptFinishDecision attempt_decision =
        prepared_attempt.decision();

    StepTime next;
    PreparedTimeFinish prepared_time;
    const Status time_prepare = implementation_->time.prepare_finish(
        implementation_->communicator, proposal, attempt_status,
        prepared_time);
    if (!time_prepare) {
      implementation_->transaction.commit_reject(prepared_attempt);
      implementation_->finalize_pending_force_cache();
      implementation_->discard_pending_attempt_side_state();
      candidate.resources = resource_difference(
          implementation_->resource_snapshot(), resources_before);
      candidate.resources.predictor_blocking_collectives =
          predictor_blocking_collectives;
      report = candidate;
      return time_prepare.code == StatusCode::rejected_step &&
                     !last_attempt_status
                 ? last_attempt_status
                 : time_prepare;
    }

    if (prepared_time.decision() == TimeFinishDecision::accept) {
      const bool consistent_accept =
          attempt_status && attempt_decision == AttemptFinishDecision::accept;
      assert(consistent_accept);
      if (!consistent_accept) {
        implementation_->transaction.commit_reject(prepared_attempt);
        implementation_->finalize_pending_force_cache();
        implementation_->discard_pending_attempt_side_state();
        return {StatusCode::invalid_plan, kProductBinding};
      }
      implementation_->transaction.commit_accept(prepared_attempt);
      implementation_->time.commit_accept(prepared_time);
      implementation_->finalize_pending_force_cache();
      implementation_->commit_pending_attempt_side_state();
      candidate.piso = attempt_report;
      candidate.accepted_step = implementation_->time.accepted_step();
      candidate.accepted_time = implementation_->time.time();
      candidate.resources = resource_difference(
          implementation_->resource_snapshot(), resources_before);
      candidate.resources.predictor_blocking_collectives =
          predictor_blocking_collectives;
      candidate.accepted = true;
      report = candidate;
      return {};
    }

    implementation_->transaction.commit_reject(prepared_attempt);
    implementation_->finalize_pending_force_cache();
    implementation_->discard_pending_attempt_side_state();
    if (prepared_time.decision() == TimeFinishDecision::fatal) {
      const Status finish = prepared_time.outcome();
      // The attempt transaction and every attempt-local side cache have
      // already been rolled back.  Retire the fatal time proposal as the
      // final no-fail cleanup step so a later advance starts a fresh BE
      // recovery proposal instead of remaining permanently active.
      implementation_->time.commit_fatal(prepared_time);
      candidate.resources = resource_difference(
          implementation_->resource_snapshot(), resources_before);
      candidate.resources.predictor_blocking_collectives =
          predictor_blocking_collectives;
      report = candidate;
      return finish.code == StatusCode::rejected_step &&
                     !last_attempt_status
                 ? last_attempt_status
                 : finish;
    }

    const bool consistent_retry =
        !attempt_status && attempt_decision == AttemptFinishDecision::reject &&
        prepared_time.decision() == TimeFinishDecision::retry;
    assert(consistent_retry);
    if (!consistent_retry)
      return {StatusCode::invalid_plan, kProductBinding};
    implementation_->time.commit_retry(prepared_time, next);
    if (candidate.numerical_failure.valid) {
      candidate.numerical_failure.retry_proposed = true;
      candidate.numerical_failure.dt_after = next.dt;
      candidate.numerical_failure.method_after = next.bdf;
      candidate.numerical_failure.origin_after = next.origin;
    }
    proposal = next;
  }
  return status;
}

Status ProductDriver::committed_output_snapshot(
    CommittedOutputSnapshot& out) noexcept {
  if (implementation_ == nullptr || !implementation_->initialized)
    return {StatusCode::invalid_plan, kProductInput};
  ProductDriver::Impl& runtime = *implementation_;
  CompiledCasePlan::Impl& product = *runtime.plan.implementation_;
  const Span<const SnapshotFieldSpec> sealed = product.io.snapshot_fields();
  for (std::size_t index = 0U; index < sealed.size; ++index) {
    ConstFieldView view;
    Status status = product.layers.view(StateRole::accepted_n,
                                        sealed.data[index].field, view);
    if (!status) return status;
    const FieldDescriptor* descriptor =
        find_descriptor(product.schema, sealed.data[index].field);
    if (descriptor == nullptr)
      return {StatusCode::invalid_plan, kProductBinding};
    runtime.output_fields[index] = {descriptor->stable_name, view,
                                    view.revision};
  }
  out = {&product.geometry,
         product.patch,
         runtime.plan.fingerprint(),
         product.schema_fingerprint,
         runtime.time.time(),
         runtime.time.accepted_step(),
         {runtime.output_fields.data(), runtime.output_fields.size()},
         true};
  return {};
}

Status ProductDriver::committed_restart_snapshot(RestartSnapshot& out) noexcept {
  if (implementation_ == nullptr || !implementation_->initialized ||
      implementation_->time.accepted_step() == 0U)
    return {StatusCode::invalid_plan, kProductInput};
  ProductDriver::Impl& runtime = *implementation_;
  CompiledCasePlan::Impl& product = *runtime.plan.implementation_;
  std::size_t index = 0U;
  const auto append = [&](RestartFieldRole role, FieldId field,
                          StateRole state,
                          std::vector<RestartFieldView>& fields,
                          std::size_t& selected) {
    ConstFieldView view;
    const Status status = product.layers.view(state, field, view);
    if (status) fields[selected++] = {role, view};
    return status;
  };
  Status status = append(RestartFieldRole::velocity, product.fields.velocity,
                         StateRole::accepted_n, runtime.restart_fields,
                         index);
  if (status)
    status = append(RestartFieldRole::pressure_perturbation,
                    product.fields.pressure, StateRole::accepted_n,
                    runtime.restart_fields, index);
  if (status)
    status = append(RestartFieldRole::enthalpy, product.fields.enthalpy,
                    StateRole::accepted_n, runtime.restart_fields, index);
  for (std::size_t scalar = 0U;
       scalar < product.fields.scalars.size() && status; ++scalar) {
    status = append(product.fields.scalar_roles[scalar] ==
                            TransportedScalarRole::species
                        ? RestartFieldRole::independent_species
                        : RestartFieldRole::transported_scalar,
                    product.fields.scalars[scalar], StateRole::accepted_n,
                    runtime.restart_fields, index);
  }
  index = 0U;
  if (status)
    status = append(RestartFieldRole::velocity, product.fields.velocity,
                    StateRole::accepted_n_minus_one,
                    runtime.restart_previous_fields, index);
  if (status)
    status = append(RestartFieldRole::pressure_perturbation,
                    product.fields.pressure,
                    StateRole::accepted_n_minus_one,
                    runtime.restart_previous_fields, index);
  if (status)
    status = append(RestartFieldRole::enthalpy, product.fields.enthalpy,
                    StateRole::accepted_n_minus_one,
                    runtime.restart_previous_fields, index);
  for (std::size_t scalar = 0U;
       scalar < product.fields.scalars.size() && status; ++scalar) {
    status = append(product.fields.scalar_roles[scalar] ==
                            TransportedScalarRole::species
                        ? RestartFieldRole::independent_species
                        : RestartFieldRole::transported_scalar,
                    product.fields.scalars[scalar],
                    StateRole::accepted_n_minus_one,
                    runtime.restart_previous_fields, index);
  }
  index = 0U;
  if (status)
    status = append(RestartFieldRole::enthalpy_nonadvective_rate,
                    product.fields.enthalpy_nonadvective_rate,
                    StateRole::accepted_n, runtime.restart_rate_fields,
                    index);
  for (std::size_t scalar = 0U;
       scalar < product.fields.scalar_nonadvective_rates.size() && status;
       ++scalar)
    status = append(RestartFieldRole::scalar_nonadvective_rate,
                    product.fields.scalar_nonadvective_rates[scalar],
                    StateRole::accepted_n, runtime.restart_rate_fields,
                    index);
  index = 0U;
  if (status)
    status = append(RestartFieldRole::enthalpy_nonadvective_rate,
                    product.fields.enthalpy_nonadvective_rate,
                    StateRole::accepted_n_minus_one,
                    runtime.restart_previous_rate_fields, index);
  for (std::size_t scalar = 0U;
       scalar < product.fields.scalar_nonadvective_rates.size() && status;
       ++scalar)
    status = append(RestartFieldRole::scalar_nonadvective_rate,
                    product.fields.scalar_nonadvective_rates[scalar],
                    StateRole::accepted_n_minus_one,
                    runtime.restart_previous_rate_fields, index);
  ConstFaceFluxView flux;
  ConstFaceFluxView previous_flux;
  if (status)
    status = runtime.final_flux_writer.committed(product.final_flux, flux);
  if (status)
    status = runtime.final_flux_writer.committed_previous(product.final_flux,
                                                          previous_flux);
  if (!status) return status;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const auto payload = [](ConstFaceFluxView value) noexcept {
    std::uint64_t hash = detail::product_mix(
        UINT64_C(0xcbf29ce484222325), UINT64_C(0x666c757868697374));
    const std::array<ConstFaceFieldView, 3U> faces{
        value.x, value.y, value.z};
    for (ConstFaceFieldView face : faces) {
      hash = detail::product_mix(hash,
                                 static_cast<std::uint8_t>(face.axis));
      hash = detail::product_mix(hash, face.extents.x);
      hash = detail::product_mix(hash, face.extents.y);
      hash = detail::product_mix(hash, face.extents.z);
      for (std::int32_t z = 0; z < face.extents.z; ++z)
        for (std::int32_t y = 0; y < face.extents.y; ++y)
          for (std::int32_t x = 0; x < face.extents.x; ++x)
            hash = detail::product_mix(
                hash, product_double_bits(face.unchecked({x, y, z})));
    }
    return hash == 0U ? PlanFingerprint{1U} : hash;
  };
  detail::ProductFinalFluxHistoryDiagnostic flux_history;
  flux_history.valid = flux.certificate.valid() &&
                       previous_flux.certificate.valid() &&
                       flux.revision != 0U && previous_flux.revision != 0U;
  flux_history.accepted_revision = flux.revision;
  flux_history.previous_revision = previous_flux.revision;
  flux_history.accepted_certificate = flux.certificate;
  flux_history.previous_certificate = previous_flux.certificate;
  flux_history.accepted_payload = payload(flux);
  flux_history.previous_payload = payload(previous_flux);
  g_final_flux_history_diagnostic = flux_history;
  g_final_flux_history_published.store(flux_history.valid,
                                       std::memory_order_release);
#endif
  out = {product.geometry.global_cells(),
         product.patch,
         runtime.plan.fingerprint(),
         product.schema_fingerprint,
         product.geometry.fingerprint(),
         runtime.time.time(),
         runtime.time.last_accepted_dt(),
         runtime.pressure_reference,
         runtime.time.accepted_step(),
         runtime.time.next_generation(),
         {runtime.restart_fields.data(), runtime.restart_fields.size()},
         flux,
         {runtime.restart_previous_fields.data(),
          runtime.restart_previous_fields.size()},
         {runtime.restart_rate_fields.data(),
          runtime.restart_rate_fields.size()},
         {runtime.restart_previous_rate_fields.data(),
          runtime.restart_previous_rate_fields.size()},
         previous_flux,
         runtime.previous_pressure_reference,
         runtime.closed_mass_target};
  return {};
}

Status ProductDriver::committed_surface_force(
    SurfaceForce& force, FinalForceCertificate& certificate) const noexcept {
  if (implementation_ == nullptr || !implementation_->initialized ||
      !implementation_->final_force_cache.has_value())
    return {StatusCode::invalid_plan, kProductInput};
  return implementation_->final_force_cache->committed(force, certificate);
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
Status ProductDriver::committed_final_mass_flux_for_test(
    ConstFaceFluxView& out) const noexcept {
  out = {};
  if (implementation_ == nullptr || !implementation_->initialized)
    return {StatusCode::invalid_plan, kProductInput};
  const ProductDriver::Impl& runtime = *implementation_;
  const CompiledCasePlan::Impl& product = *runtime.plan.implementation_;
  return runtime.final_flux_writer.committed(product.final_flux, out);
}
#endif

bool ProductDriver::initialized() const noexcept {
  return implementation_ != nullptr && implementation_->initialized;
}
double ProductDriver::pressure_reference() const noexcept {
  return implementation_ == nullptr ? 0.0
                                    : implementation_->pressure_reference;
}
double ProductDriver::closed_mass_target() const noexcept {
  return implementation_ == nullptr ? 0.0
                                    : implementation_->closed_mass_target;
}

}  // namespace hundun::v04
