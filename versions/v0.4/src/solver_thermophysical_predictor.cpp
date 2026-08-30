// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"
#include "hundun/v04_parallel.hpp"

#include "solver_equation_detail.hpp"
#include "solver_thermophysical_predictor_detail.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kPredictorPlan = 984U;
constexpr std::uint32_t kPredictorNumerical = 985U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
std::atomic<bool> g_low_bdf_source_base_failure_once{false};
#endif

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value),
                "pressure-reference identity requires binary64");
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

double bits_double(std::uint64_t bits) noexcept {
  double value = 0.0;
  static_assert(sizeof(bits) == sizeof(value),
                "predictor failure provenance requires binary64");
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

void mark_failure(
    ThermophysicalPredictorFailure& failure,
    ThermophysicalPredictorFailureReason reason,
    ThermophysicalPredictorFailureField field, std::uint32_t field_index,
    ThermophysicalAdmissibilityConstraint constraint, std::int32_t rank,
    Int3 global_index, std::uint32_t substep, bool has_cell,
    bool has_substep, bool has_margins, double low_margin = 0.0,
    double high_margin = 0.0) noexcept {
  if (failure.valid) return;
  failure = ThermophysicalPredictorFailure{};
  failure.valid = true;
  failure.reason = reason;
  failure.field = field;
  failure.constraint = constraint;
  failure.field_index = field_index;
  failure.rank = rank;
  failure.global_index = global_index;
  failure.substep = substep;
  failure.has_cell = has_cell;
  failure.has_substep = has_substep;
  failure.has_margins = has_margins;
  failure.low_margin = low_margin;
  failure.high_margin = high_margin;
}

void record_scalar(ThermophysicalPredictorFailure& failure,
                   std::uint32_t bit, double& destination,
                   double value) noexcept {
  failure.scalar_mask |= bit;
  destination = value;
}

ThermophysicalPredictorFailureField failure_field_for_constraint(
    ThermophysicalAdmissibilityConstraint constraint) noexcept {
  switch (constraint) {
    case ThermophysicalAdmissibilityConstraint::density:
      return ThermophysicalPredictorFailureField::density;
    case ThermophysicalAdmissibilityConstraint::independent_species:
      return ThermophysicalPredictorFailureField::independent_species;
    case ThermophysicalAdmissibilityConstraint::enthalpy_lower:
    case ThermophysicalAdmissibilityConstraint::enthalpy_upper:
      return ThermophysicalPredictorFailureField::enthalpy;
    case ThermophysicalAdmissibilityConstraint::dependent_species:
    case ThermophysicalAdmissibilityConstraint::none:
      return ThermophysicalPredictorFailureField::none;
  }
  return ThermophysicalPredictorFailureField::none;
}

constexpr std::size_t kFailureWireWords = 23U;

std::uint64_t signed_bits(std::int32_t value) noexcept {
  return static_cast<std::uint64_t>(static_cast<std::int64_t>(value));
}

std::int32_t bits_signed(std::uint64_t value) noexcept {
  return static_cast<std::int32_t>(static_cast<std::int64_t>(value));
}

std::array<std::uint64_t, kFailureWireWords> pack_failure(
    const ThermophysicalPredictorFailure& failure) noexcept {
  std::array<std::uint64_t, kFailureWireWords> wire{};
  wire[0U] = (failure.valid ? 1U : 0U) |
             (static_cast<std::uint64_t>(failure.reason) << 8U) |
             (static_cast<std::uint64_t>(failure.field) << 16U) |
             (static_cast<std::uint64_t>(failure.constraint) << 24U) |
             (failure.has_cell ? UINT64_C(1) << 32U : 0U) |
             (failure.has_substep ? UINT64_C(1) << 33U : 0U) |
             (failure.has_margins ? UINT64_C(1) << 34U : 0U);
  wire[1U] = failure.field_index;
  wire[2U] = signed_bits(failure.rank);
  wire[3U] = signed_bits(failure.global_index.x);
  wire[4U] = signed_bits(failure.global_index.y);
  wire[5U] = signed_bits(failure.global_index.z);
  wire[6U] = failure.substep;
  wire[7U] = failure.scalar_mask;
  wire[8U] = double_bits(failure.low_margin);
  wire[9U] = double_bits(failure.high_margin);
  wire[10U] = double_bits(failure.dt_substep);
  wire[11U] = double_bits(failure.maximum_cfl);
  wire[12U] = double_bits(failure.local_cfl);
  wire[13U] = double_bits(failure.density_current);
  wire[14U] = double_bits(failure.density_next);
  wire[15U] = double_bits(failure.quantity_current);
  wire[16U] = double_bits(failure.quantity_previous);
  wire[17U] = double_bits(failure.divergence);
  wire[18U] = double_bits(failure.nonadvective_rhs);
  wire[19U] = double_bits(failure.conserved_next);
  wire[20U] = double_bits(failure.observed_value);
  wire[21U] = double_bits(failure.allowed_lower);
  wire[22U] = double_bits(failure.allowed_upper);
  return wire;
}

void unpack_failure(const std::array<std::uint64_t, kFailureWireWords>& wire,
                    ThermophysicalPredictorFailure& failure) noexcept {
  failure = ThermophysicalPredictorFailure{};
  failure.valid = (wire[0U] & 1U) != 0U;
  failure.reason = static_cast<ThermophysicalPredictorFailureReason>(
      (wire[0U] >> 8U) & 0xffU);
  failure.field = static_cast<ThermophysicalPredictorFailureField>(
      (wire[0U] >> 16U) & 0xffU);
  failure.constraint = static_cast<ThermophysicalAdmissibilityConstraint>(
      (wire[0U] >> 24U) & 0xffU);
  failure.has_cell = (wire[0U] & (UINT64_C(1) << 32U)) != 0U;
  failure.has_substep = (wire[0U] & (UINT64_C(1) << 33U)) != 0U;
  failure.has_margins = (wire[0U] & (UINT64_C(1) << 34U)) != 0U;
  failure.field_index = static_cast<std::uint32_t>(wire[1U]);
  failure.rank = bits_signed(wire[2U]);
  failure.global_index = {bits_signed(wire[3U]), bits_signed(wire[4U]),
                          bits_signed(wire[5U])};
  failure.substep = static_cast<std::uint32_t>(wire[6U]);
  failure.scalar_mask = static_cast<std::uint32_t>(wire[7U]);
  failure.low_margin = bits_double(wire[8U]);
  failure.high_margin = bits_double(wire[9U]);
  failure.dt_substep = bits_double(wire[10U]);
  failure.maximum_cfl = bits_double(wire[11U]);
  failure.local_cfl = bits_double(wire[12U]);
  failure.density_current = bits_double(wire[13U]);
  failure.density_next = bits_double(wire[14U]);
  failure.quantity_current = bits_double(wire[15U]);
  failure.quantity_previous = bits_double(wire[16U]);
  failure.divergence = bits_double(wire[17U]);
  failure.nonadvective_rhs = bits_double(wire[18U]);
  failure.conserved_next = bits_double(wire[19U]);
  failure.observed_value = bits_double(wire[20U]);
  failure.allowed_lower = bits_double(wire[21U]);
  failure.allowed_upper = bits_double(wire[22U]);
}

std::uint64_t predictor_state_hash(
    const ThermophysicalPredictorInput& input,
    ThermophysicalPredictorOutput output) noexcept {
  std::uint64_t hash = kFnvOffset;
  hash = hash_mix(hash, double_bits(input.bdf.a0));
  hash = hash_mix(hash, double_bits(input.bdf.a1));
  hash = hash_mix(hash, double_bits(input.bdf.a2));
  hash = hash_mix(hash, input.bdf.order);
  const auto mix_view = [&](ConstFieldView view) noexcept {
    hash = hash_mix(hash, view.field);
    hash = hash_mix(hash, view.revision);
    hash = hash_mix(hash, view.storage_identity);
    hash = hash_mix(hash, view.revision_domain);
  };
  mix_view(input.density_accepted);
  mix_view(input.enthalpy_accepted);
  if (input.bdf.order == 2U) {
    mix_view(input.density_previous);
    mix_view(input.enthalpy_previous);
  }
  for (std::size_t i = 0U; i < input.species_accepted.size; ++i) {
    mix_view(input.species_accepted.data[i]);
    if (input.bdf.order == 2U) {
      mix_view(input.species_previous.data[i]);
    }
  }
  for (std::size_t i = 0U; i < input.passive_scalars_accepted.size; ++i) {
    mix_view(input.passive_scalars_accepted.data[i]);
    if (input.bdf.order == 2U) {
      mix_view(input.passive_scalars_previous.data[i]);
    }
  }
  mix_view(input.enthalpy_nonadvective_rhs.accepted);
  if (input.bdf.order == 2U) {
    mix_view(input.enthalpy_nonadvective_rhs.previous);
  }
  for (std::size_t i = 0U; i < input.species_nonadvective_rhs.size; ++i) {
    mix_view(input.species_nonadvective_rhs.data[i].accepted);
    if (input.bdf.order == 2U) {
      mix_view(input.species_nonadvective_rhs.data[i].previous);
    }
  }
  for (std::size_t i = 0U;
       i < input.passive_scalar_nonadvective_rhs.size; ++i) {
    mix_view(input.passive_scalar_nonadvective_rhs.data[i].accepted);
    if (input.bdf.order == 2U) {
      mix_view(input.passive_scalar_nonadvective_rhs.data[i].previous);
    }
  }
  mix_view(as_const(output.enthalpy));
  mix_view(as_const(output.density_workspace));
  for (std::size_t i = 0U; i < output.independent_species.size; ++i) {
    mix_view(as_const(output.independent_species.data[i]));
  }
  for (std::size_t i = 0U; i < output.passive_scalars.size; ++i) {
    mix_view(as_const(output.passive_scalars.data[i]));
  }
  hash = hash_mix(hash, input.mass_flux_accepted.certificate.authority());
  hash = hash_mix(hash, input.mass_flux_accepted.certificate.storage());
  hash = hash_mix(
      hash, input.mass_flux_accepted.certificate.revision_domain());
  if (input.bdf.order == 2U) {
    hash = hash_mix(hash, input.mass_flux_previous.certificate.authority());
    hash = hash_mix(hash, input.mass_flux_previous.certificate.storage());
    hash = hash_mix(
        hash, input.mass_flux_previous.certificate.revision_domain());
  }
  const auto mix_face_view = [&](ConstFaceFieldView view) noexcept {
    hash = hash_mix(hash, view.storage_identity);
    hash = hash_mix(hash, view.revision_domain);
    hash = hash_mix(hash,
                    static_cast<std::uint64_t>(
                        reinterpret_cast<std::uintptr_t>(view.base)));
  };
  const ConstFaceFluxView paired = as_const(output.paired_mass_flux);
  hash = hash_mix(hash, paired.revision);
  mix_face_view(paired.x);
  mix_face_view(paired.y);
  mix_face_view(paired.z);
  return hash == 0U ? 1U : hash;
}

bool valid_rate(ConstFieldView view, Int3 cells, bool required) noexcept {
  if (view.base == nullptr) {
    return !required && view.field == 0U && view.revision == 0U &&
           view.storage_identity == 0U && view.revision_domain == 0U;
  }
  return detail::valid_cell_view(view, cells, 0U, 1U, 0U);
}

bool rate_is_zero(ConstFieldView view) noexcept { return view.base == nullptr; }

bool empty_field(ConstFieldView view) noexcept {
  return view.base == nullptr && view.field == 0U && view.revision == 0U &&
         view.storage_identity == 0U && view.revision_domain == 0U;
}

bool output_aliases_flux(FieldView output, ConstFaceFluxView flux) noexcept {
  return detail::cell_face_views_overlap(output, flux.x) ||
         detail::cell_face_views_overlap(output, flux.y) ||
         detail::cell_face_views_overlap(output, flux.z);
}

bool output_aliases_input(FieldView output, ConstFieldView input) noexcept {
  return input.base != nullptr &&
         detail::field_views_overlap(as_const(output), input);
}

bool bdf_matches_dt(double dt, BdfCoefficients bdf,
                    double& extrapolate_accepted,
                    double& extrapolate_previous) noexcept {
  if (!std::isfinite(dt) || dt <= 0.0 ||
      !detail::valid_bdf_coefficients(bdf)) {
    return false;
  }
  const double scale =
      std::max({1.0, std::abs(bdf.a0), std::abs(bdf.a1), std::abs(bdf.a2)});
  const double tolerance =
      256.0 * std::numeric_limits<double>::epsilon() * scale;
  if (bdf.order == 1U) {
    if (std::abs(bdf.a0 - 1.0 / dt) > tolerance ||
        std::abs(bdf.a1 + 1.0 / dt) > tolerance || bdf.a2 != 0.0) {
      return false;
    }
    extrapolate_accepted = 1.0;
    extrapolate_previous = 0.0;
    return true;
  }
  const double ratio = -bdf.a1 * dt - 1.0;
  if (!std::isfinite(ratio) || ratio <= 0.0) {
    return false;
  }
  const double denominator = dt * (1.0 + ratio);
  const double expected_a0 = (1.0 + 2.0 * ratio) / denominator;
  const double expected_a2 = ratio * ratio / denominator;
  if (std::abs(bdf.a0 - expected_a0) > tolerance ||
      std::abs(bdf.a2 - expected_a2) > tolerance) {
    return false;
  }
  extrapolate_accepted = 1.0 + ratio;
  extrapolate_previous = -ratio;
  return true;
}

double flux_divergence_density(const CartesianKernelPlan& kernels,
                               ConstFaceFluxView flux, Int3 cell) noexcept {
  const double integral =
      flux.x.unchecked({cell.x + 1, cell.y, cell.z}) -
      flux.x.unchecked(cell) +
      flux.y.unchecked({cell.x, cell.y + 1, cell.z}) -
      flux.y.unchecked(cell) +
      flux.z.unchecked({cell.x, cell.y, cell.z + 1}) -
      flux.z.unchecked(cell);
  return integral / detail::cell_volume(kernels, cell);
}

bool valid_flux(ConstFaceFluxView flux, Int3 cells) noexcept {
  return detail::valid_flux_view(flux, cells, flux.revision) &&
         flux.certificate.valid() && flux.certificate.matches(flux);
}

bool ghost_authority_matches(ThermophysicalGhostAuthority authority,
                             ConstFieldView view, RevisionToken geometry,
                             RevisionToken boundary,
                             std::uint8_t required_reach) noexcept {
  return authority.valid() && authority.field == view.field &&
         authority.state == view.revision &&
         authority.storage == view.storage_identity &&
         authority.revision_domain == view.revision_domain &&
         authority.geometry == geometry && authority.boundary == boundary &&
         authority.reach >= required_reach;
}

double first_order_upwind_divergence(const CartesianKernelPlan& kernels,
                                     ConstFaceFluxView flux,
                                     ConstFieldView quantity,
                                     Int3 cell) noexcept {
  const auto face_value = [&](double mass_flux, Int3 negative,
                              Int3 positive) noexcept {
    return mass_flux >= 0.0 ? quantity.unchecked(negative, 0U)
                            : quantity.unchecked(positive, 0U);
  };
  const Int3 xp{cell.x + 1, cell.y, cell.z};
  const Int3 xm{cell.x - 1, cell.y, cell.z};
  const Int3 yp{cell.x, cell.y + 1, cell.z};
  const Int3 ym{cell.x, cell.y - 1, cell.z};
  const Int3 zp{cell.x, cell.y, cell.z + 1};
  const Int3 zm{cell.x, cell.y, cell.z - 1};
  const double fxp = flux.x.unchecked(xp);
  const double fxm = flux.x.unchecked(cell);
  const double fyp = flux.y.unchecked(yp);
  const double fym = flux.y.unchecked(cell);
  const double fzp = flux.z.unchecked(zp);
  const double fzm = flux.z.unchecked(cell);
  const double integral =
      fxp * face_value(fxp, cell, xp) -
      fxm * face_value(fxm, xm, cell) +
      fyp * face_value(fyp, cell, yp) -
      fym * face_value(fym, ym, cell) +
      fzp * face_value(fzp, cell, zp) -
      fzm * face_value(fzm, zm, cell);
  return integral / detail::cell_volume(kernels, cell);
}

double outgoing_mass_flux(ConstFaceFluxView flux, Int3 cell) noexcept {
  const double right = flux.x.unchecked({cell.x + 1, cell.y, cell.z});
  const double left = flux.x.unchecked(cell);
  const double top = flux.y.unchecked({cell.x, cell.y + 1, cell.z});
  const double bottom = flux.y.unchecked(cell);
  const double front = flux.z.unchecked({cell.x, cell.y, cell.z + 1});
  const double back = flux.z.unchecked(cell);
  return std::max(right, 0.0) + std::max(-left, 0.0) +
         std::max(top, 0.0) + std::max(-bottom, 0.0) +
         std::max(front, 0.0) + std::max(-back, 0.0);
}

Status broadcast_failure(MPI_Comm communicator, int root, int rank,
                          ThermophysicalPredictorFailure& failure) noexcept {
  std::array<std::uint64_t, kFailureWireWords> wire{};
  if (rank == root) wire = pack_failure(failure);
  if (MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT64_T, root,
                communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPredictorPlan};
  }
  unpack_failure(wire, failure);
  return {};
}

Status collective_status(MPI_Comm communicator, Status local, int rank,
                         int size,
                         ThermophysicalPredictorFailure& local_failure,
                         ThermophysicalPredictorFailure& selected_failure) noexcept {
  const int candidate = local ? size : rank;
  int lowest = size;
  if (MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
      MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPredictorPlan};
  }
  if (lowest == size) return {};
  std::uint64_t encoded = 0U;
  if (rank == lowest) {
    encoded = (static_cast<std::uint64_t>(local.code) << 32U) | local.detail;
  }
  if (MPI_Bcast(&encoded, 1, MPI_UINT64_T, lowest, communicator) !=
      MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPredictorPlan};
  }
  const Status published{static_cast<StatusCode>(encoded >> 32U),
                         static_cast<std::uint32_t>(encoded)};
  if (!published && published.code == StatusCode::numerical_failure) {
    Status status = broadcast_failure(communicator, lowest, rank,
                                      local_failure);
    if (!status) return status;
    selected_failure = local_failure;
  }
  return published;
}

Status collective_high_state(MPI_Comm communicator, Status local,
                             bool locally_admissible, int rank, int size,
                             bool& globally_admissible,
                             ThermophysicalPredictorFailure& local_failure,
                             ThermophysicalPredictorFailure& selected_failure) noexcept {
  std::array<int, 2U> values{{local ? size : rank,
                             locally_admissible ? 1 : 0}};
  if (MPI_Allreduce(MPI_IN_PLACE, values.data(),
                    static_cast<int>(values.size()), MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPredictorPlan};
  }
  if (values[0U] != size) {
    std::uint64_t encoded = 0U;
    if (rank == values[0U]) {
      encoded = (static_cast<std::uint64_t>(local.code) << 32U) | local.detail;
    }
    if (MPI_Bcast(&encoded, 1, MPI_UINT64_T, values[0U], communicator) !=
        MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kPredictorPlan};
    }
    const Status published{static_cast<StatusCode>(encoded >> 32U),
                           static_cast<std::uint32_t>(encoded)};
    if (!published && published.code == StatusCode::numerical_failure) {
      Status status = broadcast_failure(communicator, values[0U], rank,
                                        local_failure);
      if (!status) return status;
      selected_failure = local_failure;
    }
    return published;
  }
  globally_admissible = values[1U] != 0;
  return {};
}

}  // namespace

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace detail {

void arm_low_bdf_source_base_failure_once_for_test() noexcept {
  g_low_bdf_source_base_failure_once.store(true, std::memory_order_relaxed);
}

void clear_low_bdf_source_base_failure_for_test() noexcept {
  g_low_bdf_source_base_failure_once.store(false, std::memory_order_relaxed);
}

}  // namespace detail
#endif

Status ThermophysicalPredictorPlan::predict(
    MPI_Comm communicator, Status prerequisite,
    const ThermophysicalPredictorInput& input,
    ThermophysicalPredictorOutput output,
    const ThermophysicalPredictorSlowPath& slow_path,
    ThermophysicalPredictorDiagnostics& diagnostics,
    ThermophysicalPredictorCertificate& certificate) const noexcept {
  int rank = -1;
  int size = 0;
  std::uint32_t blocking_collectives = 0U;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kPredictorPlan};
  }

  ThermophysicalPredictorFailure failure;
  LinearSolveResult enthalpy_endpoint_result;
  double enthalpy_endpoint_alpha = 1.0;
  double bdf_endpoint_alpha = 1.0;
  double source_endpoint_alpha = 1.0;
  std::uint8_t enthalpy_solve_calls = 0U;
  const auto publish_failure = [&](Status status) noexcept {
    diagnostics.enthalpy_endpoint = enthalpy_endpoint_result;
    diagnostics.enthalpy_endpoint_alpha = enthalpy_endpoint_alpha;
    diagnostics.bdf_endpoint_alpha = bdf_endpoint_alpha;
    diagnostics.source_endpoint_alpha = source_endpoint_alpha;
    diagnostics.enthalpy_solve_calls = enthalpy_solve_calls;
    diagnostics.blocking_collectives = blocking_collectives;
    if (!status && status.code == StatusCode::numerical_failure &&
        failure.valid) {
      diagnostics.failure = failure;
    }
    return status;
  };
  const bool second_order = input.bdf.order == 2U;
  Status local = prerequisite;
  const auto reject_plan = [&]() noexcept {
    if (local) local = {StatusCode::invalid_plan, kPredictorPlan};
  };
  const std::size_t species_count = species_.size();
  const std::size_t passive_count = passive_scalars_.size();
  if (boundary_ == nullptr || species_enthalpy_minimum_.size() != species_count ||
      species_enthalpy_maximum_.size() != species_count ||
      !std::isfinite(dependent_enthalpy_minimum_) ||
      !std::isfinite(dependent_enthalpy_maximum_) ||
      input.species_accepted.size != species_count ||
      input.species_previous.size != species_count ||
      input.passive_scalars_accepted.size != passive_count ||
      input.passive_scalars_previous.size != passive_count ||
      input.species_ghosts.size != species_count ||
      input.passive_scalar_ghosts.size != passive_count ||
      (species_count != 0U &&
       (input.species_accepted.data == nullptr ||
        input.species_previous.data == nullptr ||
        input.species_ghosts.data == nullptr)) ||
      (passive_count != 0U &&
       (input.passive_scalars_accepted.data == nullptr ||
        input.passive_scalars_previous.data == nullptr ||
        input.passive_scalar_ghosts.data == nullptr))) {
    reject_plan();
  }
  std::uintptr_t exchange_plan = 0U;
  if (local) {
    const ThermophysicalGhostHistory authority = input.enthalpy_ghosts;
    if (!ghost_authority_matches(authority.accepted,
                                 input.enthalpy_accepted, input.geometry,
                                 input.boundary, enthalpy_reach_) ||
        (second_order
             ? (!ghost_authority_matches(authority.previous,
                                         input.enthalpy_previous,
                                         input.geometry, input.boundary,
                                         enthalpy_reach_) ||
                authority.previous.exchange_plan !=
                    authority.accepted.exchange_plan)
             : false)) {
      reject_plan();
    } else {
      exchange_plan = authority.accepted.exchange_plan;
    }
  }
  for (std::size_t index = 0U; index < species_count && local; ++index) {
    const ThermophysicalGhostHistory authority =
        input.species_ghosts.data[index];
    if (!ghost_authority_matches(authority.accepted,
                                 input.species_accepted.data[index],
                                 input.geometry, input.boundary,
                                 species_reach_) ||
        authority.accepted.exchange_plan != exchange_plan ||
        (second_order
             ? (!ghost_authority_matches(
                    authority.previous, input.species_previous.data[index],
                    input.geometry, input.boundary, species_reach_) ||
                authority.previous.exchange_plan != exchange_plan)
             : false)) {
      reject_plan();
    }
  }
  for (std::size_t index = 0U; index < passive_count && local; ++index) {
    const ThermophysicalGhostHistory authority =
        input.passive_scalar_ghosts.data[index];
    if (!ghost_authority_matches(authority.accepted,
                                 input.passive_scalars_accepted.data[index],
                                 input.geometry, input.boundary,
                                 passive_scalar_reach_) ||
        authority.accepted.exchange_plan != exchange_plan ||
        (second_order
             ? (!ghost_authority_matches(
                    authority.previous,
                    input.passive_scalars_previous.data[index],
                    input.geometry, input.boundary,
                    passive_scalar_reach_) ||
                authority.previous.exchange_plan != exchange_plan)
             : false)) {
      reject_plan();
    }
  }

  ThermophysicalPredictorCertificate high_certificate;
  if (local) {
    local = predict_high_local(input, output, high_certificate, failure);
    if (failure.valid && failure.rank < 0) failure.rank = rank;
  }

  const auto audit_paired_mass = [&]() noexcept {
    for (std::int32_t z = 0; z < cells_.z; ++z) {
      for (std::int32_t y = 0; y < cells_.y; ++y) {
        for (std::int32_t x = 0; x < cells_.x; ++x) {
          const Int3 cell{x, y, z};
          const double volume = detail::cell_volume(*kernels_, cell);
          const double rho_star =
              output.density_workspace.unchecked(cell, 0U);
          const double rho_n = input.density_accepted.unchecked(cell, 0U);
          const double rho_nm1 =
              input.bdf.order == 2U
                  ? input.density_previous.unchecked(cell, 0U)
                  : 0.0;
          const ConstFaceFluxView paired = as_const(output.paired_mass_flux);
          const double flux_x_positive =
              paired.x.unchecked({cell.x + 1, cell.y, cell.z});
          const double flux_x_negative = paired.x.unchecked(cell);
          const double flux_y_positive =
              paired.y.unchecked({cell.x, cell.y + 1, cell.z});
          const double flux_y_negative = paired.y.unchecked(cell);
          const double flux_z_positive =
              paired.z.unchecked({cell.x, cell.y, cell.z + 1});
          const double flux_z_negative = paired.z.unchecked(cell);
          const double integral =
              flux_x_positive - flux_x_negative + flux_y_positive -
              flux_y_negative + flux_z_positive - flux_z_negative;
          const double storage =
              volume * (input.bdf.a0 * rho_star + input.bdf.a1 * rho_n +
                        input.bdf.a2 * rho_nm1);
          const double residual = storage + integral;
          // The paired identity is a cancellation of three BDF storage
          // history terms and six oriented face terms.  Its backward-error
          // scale must retain those operands: using the already-cancelled
          // storage/integral makes an exact small-dt identity look O(1) while
          // its individual BDF terms are O(1/dt).
          const double storage_term_magnitude =
              std::abs(volume * input.bdf.a0 * rho_star) +
              std::abs(volume * input.bdf.a1 * rho_n) +
              std::abs(volume * input.bdf.a2 * rho_nm1);
          const double flux_term_magnitude =
              std::abs(flux_x_positive) + std::abs(flux_x_negative) +
              std::abs(flux_y_positive) + std::abs(flux_y_negative) +
              std::abs(flux_z_positive) + std::abs(flux_z_negative);
          const double scale =
              std::max(1.0, storage_term_magnitude + flux_term_magnitude);
          if (!std::isfinite(volume) || !(volume > 0.0) ||
              !std::isfinite(residual) || !std::isfinite(scale) ||
              std::abs(residual) >
                  2048.0 * std::numeric_limits<double>::epsilon() * scale) {
            mark_failure(
                failure,
                ThermophysicalPredictorFailureReason::blended_density,
                ThermophysicalPredictorFailureField::mass_flux, UINT32_MAX,
                ThermophysicalAdmissibilityConstraint::density, rank,
                {cell.x + patch_begin_.x, cell.y + patch_begin_.y,
                 cell.z + patch_begin_.z},
                0U, true, false, false);
            if (failure.valid) {
              record_scalar(
                  failure, ThermophysicalPredictorFailure::scalar_divergence,
                  failure.divergence, integral / volume);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_observed_value,
                  failure.observed_value, residual);
            }
            return Status{StatusCode::numerical_failure,
                          kPredictorNumerical};
          }
        }
      }
    }
    return Status{};
  };
  if (local && slow_path.immersed_interface != nullptr)
    local = slow_path.immersed_interface->zero_interface_flux(
        output.paired_mass_flux);
  if (local) local = audit_paired_mass();

  const Int3 end{cells_.x, cells_.y, cells_.z};
  const auto tuple_admissible =
      [&](ConstFieldView density, ConstFieldView enthalpy,
          Span<FieldView> independent_species,
          Span<FieldView> passive_scalars,
          bool conserved_quantities = false) noexcept {
        if (independent_species.size != species_count ||
            passive_scalars.size != passive_count) {
          return false;
        }
        for (std::int32_t z = 0; z < end.z; ++z) {
          for (std::int32_t y = 0; y < end.y; ++y) {
            for (std::int32_t x = 0; x < end.x; ++x) {
              const Int3 cell{x, y, z};
              const double rho = density.unchecked(cell, 0U);
              const double stored_h = enthalpy.unchecked(cell, 0U);
              const double h = conserved_quantities ? stored_h / rho
                                                    : stored_h;
              if (!std::isfinite(rho) || !(rho > 0.0) ||
                  !std::isfinite(h)) {
                return false;
              }
              double dependent_density = rho;
              double lower = rho * dependent_enthalpy_minimum_;
              double upper = rho * dependent_enthalpy_maximum_;
              for (std::size_t index = 0U; index < species_count; ++index) {
                const double stored_value =
                    independent_species.data[index].unchecked(cell, 0U);
                const double value = conserved_quantities
                                         ? stored_value / rho
                                         : stored_value;
                const double species_density =
                    conserved_quantities ? stored_value : rho * value;
                if (!std::isfinite(value) ||
                    !std::isfinite(species_density) ||
                    species_density < 0.0) {
                  return false;
                }
                dependent_density -= species_density;
                lower += species_density *
                         (species_enthalpy_minimum_[index] -
                          dependent_enthalpy_minimum_);
                upper += species_density *
                         (species_enthalpy_maximum_[index] -
                          dependent_enthalpy_maximum_);
              }
              const double density_enthalpy =
                  conserved_quantities ? stored_h : rho * h;
              if (!std::isfinite(dependent_density) ||
                  dependent_density < 0.0 || !std::isfinite(lower) ||
                  !std::isfinite(upper) ||
                  !std::isfinite(density_enthalpy) ||
                  density_enthalpy < lower || density_enthalpy > upper) {
                return false;
              }
              for (std::size_t index = 0U; index < passive_count; ++index) {
                const double stored_value =
                    passive_scalars.data[index].unchecked(cell, 0U);
                const double value = conserved_quantities
                                         ? stored_value / rho
                                         : stored_value;
                if (!std::isfinite(value)) {
                  return false;
                }
              }
            }
          }
        }
        return true;
      };

  const auto tuple_nonenthalpy_admissible =
      [&](ConstFieldView density, Span<FieldView> independent_species,
          Span<FieldView> passive_scalars) noexcept {
        if (independent_species.size != species_count ||
            passive_scalars.size != passive_count) {
          return false;
        }
        for (std::int32_t z = 0; z < end.z; ++z) {
          for (std::int32_t y = 0; y < end.y; ++y) {
            for (std::int32_t x = 0; x < end.x; ++x) {
              const Int3 cell{x, y, z};
              const double rho = density.unchecked(cell, 0U);
              if (!std::isfinite(rho) || !(rho > 0.0)) return false;
              double dependent_density = rho;
              for (std::size_t index = 0U; index < species_count; ++index) {
                const double value =
                    independent_species.data[index].unchecked(cell, 0U);
                const double species_density = rho * value;
                if (!std::isfinite(value) ||
                    !std::isfinite(species_density) ||
                    species_density < 0.0) {
                  return false;
                }
                dependent_density -= species_density;
              }
              if (!std::isfinite(dependent_density) ||
                  dependent_density < 0.0) {
                return false;
              }
              for (std::size_t index = 0U; index < passive_count; ++index) {
                if (!std::isfinite(
                        passive_scalars.data[index].unchecked(cell, 0U))) {
                  return false;
                }
              }
            }
          }
        }
        return true;
      };

  // This is deliberately a failure-only second pass.  The admissibility
  // oracle above remains the sole normal-path check and keeps its original
  // traversal and arithmetic.  On a rejected endpoint this pass identifies
  // the first scalar/constraint in that same z/y/x, quantity order so that
  // the fixed provenance record is useful without paying for a scan on
  // successful calls.
  const auto capture_tuple_failure =
      [&](ConstFieldView density, ConstFieldView enthalpy,
          Span<FieldView> independent_species, Span<FieldView> passive_scalars,
          ThermophysicalPredictorFailureReason reason,
          std::uint32_t one_based_substep,
          bool conserved_quantities = false) noexcept {
        if (failure.valid) return;
        const bool has_substep = one_based_substep != 0U;
        if (independent_species.size != species_count ||
            passive_scalars.size != passive_count) {
          mark_failure(failure, reason,
                       ThermophysicalPredictorFailureField::none,
                       UINT32_MAX,
                       ThermophysicalAdmissibilityConstraint::none, rank,
                       {-1, -1, -1}, 0U, false, has_substep, false);
          return;
        }
        for (std::int32_t z = 0; z < end.z; ++z) {
          for (std::int32_t y = 0; y < end.y; ++y) {
            for (std::int32_t x = 0; x < end.x; ++x) {
              const Int3 cell{x, y, z};
              const Int3 global_cell{cell.x + patch_begin_.x,
                                     cell.y + patch_begin_.y,
                                     cell.z + patch_begin_.z};
              const double rho = density.unchecked(cell, 0U);
              const double stored_h = enthalpy.unchecked(cell, 0U);
              const double h = conserved_quantities ? stored_h / rho : stored_h;
              if (!std::isfinite(rho) || !(rho > 0.0)) {
                mark_failure(
                    failure, reason,
                    ThermophysicalPredictorFailureField::density, 0U,
                    ThermophysicalAdmissibilityConstraint::density, rank,
                    global_cell, one_based_substep, true, has_substep,
                    std::isfinite(rho), std::isfinite(rho) ? rho : 0.0,
                    0.0);
                if (failure.valid && std::isfinite(rho)) {
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_observed_value,
                                failure.observed_value, rho);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_allowed_lower,
                                failure.allowed_lower, 0.0);
                }
                return;
              }
              if (!std::isfinite(h)) {
                mark_failure(
                    failure, reason,
                    ThermophysicalPredictorFailureField::enthalpy, 0U,
                    ThermophysicalAdmissibilityConstraint::none, rank,
                    global_cell, one_based_substep, true, has_substep, false);
                if (failure.valid) {
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_observed_value,
                                failure.observed_value, h);
                }
                return;
              }
              double dependent_density = rho;
              double lower = rho * dependent_enthalpy_minimum_;
              double upper = rho * dependent_enthalpy_maximum_;
              for (std::size_t index = 0U; index < species_count; ++index) {
                const double stored_value =
                    independent_species.data[index].unchecked(cell, 0U);
                const double value = conserved_quantities
                                         ? stored_value / rho
                                         : stored_value;
                const double species_density =
                    conserved_quantities ? stored_value : rho * value;
                if (!std::isfinite(value) ||
                    !std::isfinite(species_density) ||
                    species_density < 0.0) {
                  mark_failure(
                      failure, reason,
                      ThermophysicalPredictorFailureField::independent_species,
                      static_cast<std::uint32_t>(index),
                      ThermophysicalAdmissibilityConstraint::
                          independent_species,
                      rank, global_cell, one_based_substep, true, has_substep,
                      false);
                  if (failure.valid) {
                    record_scalar(failure,
                                  ThermophysicalPredictorFailure::
                                      scalar_density_current,
                                  failure.density_current, rho);
                    record_scalar(failure,
                                  ThermophysicalPredictorFailure::
                                      scalar_observed_value,
                                  failure.observed_value, value);
                    record_scalar(failure,
                                  ThermophysicalPredictorFailure::
                                      scalar_allowed_lower,
                                  failure.allowed_lower, 0.0);
                  }
                  return;
                }
                dependent_density -= species_density;
                lower += species_density *
                         (species_enthalpy_minimum_[index] -
                          dependent_enthalpy_minimum_);
                upper += species_density *
                         (species_enthalpy_maximum_[index] -
                          dependent_enthalpy_maximum_);
              }
              if (!std::isfinite(dependent_density) ||
                  dependent_density < 0.0) {
                mark_failure(
                    failure, reason,
                    ThermophysicalPredictorFailureField::none, UINT32_MAX,
                    ThermophysicalAdmissibilityConstraint::dependent_species,
                    rank, global_cell, one_based_substep, true, has_substep,
                    std::isfinite(dependent_density),
                    std::isfinite(dependent_density) ? dependent_density
                                                      : 0.0,
                    0.0);
                if (failure.valid && std::isfinite(dependent_density)) {
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_observed_value,
                                failure.observed_value, dependent_density);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_allowed_lower,
                                failure.allowed_lower, 0.0);
                }
                return;
              }
              const double density_enthalpy =
                  conserved_quantities ? stored_h : rho * h;
              if (!std::isfinite(lower) || !std::isfinite(upper) ||
                  !std::isfinite(density_enthalpy)) {
                mark_failure(
                    failure, reason,
                    ThermophysicalPredictorFailureField::enthalpy, 0U,
                    ThermophysicalAdmissibilityConstraint::none, rank,
                    global_cell, one_based_substep, true, has_substep, false);
                if (failure.valid) {
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_observed_value,
                                failure.observed_value, density_enthalpy);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_conserved_next,
                                failure.conserved_next, density_enthalpy);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_quantity_current,
                                failure.quantity_current, h);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_allowed_lower,
                                failure.allowed_lower, lower);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_allowed_upper,
                                failure.allowed_upper, upper);
                }
                return;
              }
              if (density_enthalpy < lower) {
                mark_failure(
                    failure, reason,
                    ThermophysicalPredictorFailureField::enthalpy, 0U,
                    ThermophysicalAdmissibilityConstraint::enthalpy_lower,
                    rank, global_cell, one_based_substep, true, has_substep,
                    true, density_enthalpy - lower,
                    upper - density_enthalpy);
                if (failure.valid) {
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_density_current,
                                failure.density_current, rho);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_observed_value,
                                failure.observed_value, density_enthalpy);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_conserved_next,
                                failure.conserved_next, density_enthalpy);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_quantity_current,
                                failure.quantity_current, h);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_allowed_lower,
                                failure.allowed_lower, lower);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_allowed_upper,
                                failure.allowed_upper, upper);
                }
                return;
              }
              if (density_enthalpy > upper) {
                mark_failure(
                    failure, reason,
                    ThermophysicalPredictorFailureField::enthalpy, 0U,
                    ThermophysicalAdmissibilityConstraint::enthalpy_upper,
                    rank, global_cell, one_based_substep, true, has_substep,
                    true, density_enthalpy - lower,
                    upper - density_enthalpy);
                if (failure.valid) {
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_density_current,
                                failure.density_current, rho);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_observed_value,
                                failure.observed_value, density_enthalpy);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_allowed_lower,
                                failure.allowed_lower, lower);
                  record_scalar(failure,
                                ThermophysicalPredictorFailure::
                                    scalar_allowed_upper,
                                failure.allowed_upper, upper);
                }
                return;
              }
              for (std::size_t index = 0U; index < passive_count; ++index) {
                const double stored_value =
                    passive_scalars.data[index].unchecked(cell, 0U);
                const double value = conserved_quantities
                                         ? stored_value / rho
                                         : stored_value;
                if (!std::isfinite(value)) {
                  mark_failure(
                      failure, reason,
                      ThermophysicalPredictorFailureField::passive_scalar,
                      static_cast<std::uint32_t>(index),
                      ThermophysicalAdmissibilityConstraint::none, rank,
                      global_cell, one_based_substep, true, has_substep, false);
                  if (failure.valid) {
                    record_scalar(failure,
                                  ThermophysicalPredictorFailure::
                                      scalar_observed_value,
                                  failure.observed_value, value);
                  }
                  return;
                }
              }
            }
          }
        }
        mark_failure(failure, reason,
                     ThermophysicalPredictorFailureField::none, UINT32_MAX,
                     ThermophysicalAdmissibilityConstraint::none, rank,
                     {-1, -1, -1}, 0U, false, has_substep, false);
      };

  bool locally_high_admissible = false;
  if (local) {
    locally_high_admissible = tuple_admissible(
        as_const(output.density_workspace), as_const(output.enthalpy),
        output.independent_species, output.passive_scalars);
  }
  bool globally_high_admissible = false;
  Status consensus = collective_high_state(
      communicator, local, locally_high_admissible, rank, size,
      globally_high_admissible, failure, failure);
  if (!consensus) return publish_failure(consensus);
  ++blocking_collectives;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (second_order && g_low_bdf_source_base_failure_once.exchange(
                          false, std::memory_order_relaxed)) {
    mark_failure(
        failure,
        ThermophysicalPredictorFailureReason::
            low_bdf_source_base_admissibility,
        ThermophysicalPredictorFailureField::enthalpy, UINT32_MAX,
        ThermophysicalAdmissibilityConstraint::enthalpy_lower, rank,
        patch_begin_, 1U, true, true, false);
    local = {StatusCode::numerical_failure, kPredictorNumerical};
    consensus = collective_status(communicator, local, rank, size, failure,
                                  failure);
    ++blocking_collectives;
    if (!consensus) return publish_failure(consensus);
  }
#endif
  if (globally_high_admissible) {
    diagnostics = {};
    diagnostics.blocking_collectives = blocking_collectives;
    certificate = high_certificate;
    return {};
  }
  const int local_high_nonenthalpy =
      local && tuple_nonenthalpy_admissible(
                   as_const(output.density_workspace),
                   output.independent_species, output.passive_scalars)
          ? 1
          : 0;
  int global_high_nonenthalpy = 0;
  if (MPI_Allreduce(&local_high_nonenthalpy, &global_high_nonenthalpy, 1,
                    MPI_INT, MPI_MIN, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPredictorPlan};
  }
  ++blocking_collectives;
  double local_accepted_flux = 0.0;
  for (std::int32_t z = 0; z < end.z; ++z) {
    for (std::int32_t y = 0; y < end.y; ++y) {
      for (std::int32_t x = 0; x < end.x; ++x) {
        local_accepted_flux =
            std::max(local_accepted_flux,
                     outgoing_mass_flux(input.mass_flux_accepted,
                                        {x, y, z}));
      }
    }
  }
  double global_accepted_flux = 0.0;
  if (MPI_Allreduce(&local_accepted_flux, &global_accepted_flux, 1,
                    MPI_DOUBLE, MPI_MAX, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPredictorPlan};
  }
  ++blocking_collectives;

  const auto valid_low_field = [&](FieldView view, FieldId field,
                                   std::uint8_t reach) noexcept {
    return view.field == field &&
           detail::valid_cell_view(as_const(view), cells_, 0U, 1U, reach);
  };
  if (!valid_low_field(output.low_order_density_workspace, density_, 0U) ||
      !valid_low_field(output.low_order_enthalpy_workspace, enthalpy_, 1U) ||
      output.low_order_independent_species.size != species_count ||
      output.low_order_passive_scalars.size != passive_count ||
      (species_count != 0U &&
       output.low_order_independent_species.data == nullptr) ||
      (passive_count != 0U && output.low_order_passive_scalars.data == nullptr)) {
    local = {StatusCode::invalid_plan, kPredictorPlan};
  }
  for (std::size_t index = 0U; index < species_count && local; ++index) {
    if (!valid_low_field(output.low_order_independent_species.data[index],
                         species_[index], 1U)) {
      local = {StatusCode::invalid_plan, kPredictorPlan};
    }
  }
  for (std::size_t index = 0U; index < passive_count && local; ++index) {
    if (!valid_low_field(output.low_order_passive_scalars.data[index],
                         passive_scalars_[index], 1U)) {
      local = {StatusCode::invalid_plan, kPredictorPlan};
    }
  }
  const auto low_aliases = [&](ConstFieldView low, ConstFieldView high) {
    return detail::field_views_overlap(low, high);
  };
  if (local &&
      (low_aliases(as_const(output.low_order_density_workspace),
                   as_const(output.density_workspace)) ||
       low_aliases(as_const(output.low_order_enthalpy_workspace),
                   as_const(output.enthalpy)))) {
    local = {StatusCode::invalid_plan, kPredictorPlan};
  }
  for (std::size_t index = 0U; index < species_count && local; ++index) {
    if (low_aliases(
            as_const(output.low_order_independent_species.data[index]),
            as_const(output.independent_species.data[index]))) {
      local = {StatusCode::invalid_plan, kPredictorPlan};
    }
  }
  for (std::size_t index = 0U; index < passive_count && local; ++index) {
    if (low_aliases(as_const(output.low_order_passive_scalars.data[index]),
                    as_const(output.passive_scalars.data[index]))) {
      local = {StatusCode::invalid_plan, kPredictorPlan};
    }
  }
  const std::size_t low_field_count = 2U + species_count + passive_count;
  const auto low_field = [&](std::size_t index) noexcept -> FieldView {
    if (index == 0U) return output.low_order_density_workspace;
    if (index == 1U) return output.low_order_enthalpy_workspace;
    index -= 2U;
    if (index < species_count)
      return output.low_order_independent_species.data[index];
    return output.low_order_passive_scalars.data[index - species_count];
  };
  for (std::size_t left = 0U; left < low_field_count && local; ++left) {
    const FieldView candidate = low_field(left);
    for (std::size_t right = left + 1U; right < low_field_count; ++right) {
      if (detail::field_views_overlap(as_const(candidate),
                                      as_const(low_field(right)))) {
        local = {StatusCode::invalid_plan, kPredictorPlan};
        break;
      }
    }
    const std::array<FieldView, 4U> primary_outputs{
        output.enthalpy, output.density_workspace,
        output.accepted_advection_workspace,
        output.previous_advection_workspace};
    for (FieldView high : primary_outputs) {
      if (detail::field_views_overlap(as_const(candidate), as_const(high))) {
        local = {StatusCode::invalid_plan, kPredictorPlan};
        break;
      }
    }
    for (std::size_t index = 0U; index < species_count && local; ++index) {
      if (detail::field_views_overlap(
              as_const(candidate),
              as_const(output.independent_species.data[index]))) {
        local = {StatusCode::invalid_plan, kPredictorPlan};
      }
    }
    for (std::size_t index = 0U; index < passive_count && local; ++index) {
      if (detail::field_views_overlap(
              as_const(candidate), as_const(output.passive_scalars.data[index]))) {
        local = {StatusCode::invalid_plan, kPredictorPlan};
      }
    }
    const auto aliases = [&](ConstFieldView input_view) noexcept {
      return output_aliases_input(candidate, input_view);
    };
    if (local &&
        (aliases(input.density_accepted) ||
         aliases(input.enthalpy_accepted) ||
         aliases(input.enthalpy_nonadvective_rhs.accepted) ||
         (second_order &&
          (aliases(input.density_previous) ||
           aliases(input.enthalpy_previous) ||
           aliases(input.enthalpy_nonadvective_rhs.previous))) ||
         output_aliases_flux(candidate, input.mass_flux_accepted) ||
         (second_order &&
          output_aliases_flux(candidate, input.mass_flux_previous)))) {
      local = {StatusCode::invalid_plan, kPredictorPlan};
    }
    for (std::size_t index = 0U; index < species_count && local; ++index) {
      const PredictorRateHistory rate =
          input.species_nonadvective_rhs.data[index];
      if (aliases(input.species_accepted.data[index]) ||
          aliases(rate.accepted) ||
          (second_order &&
           (aliases(input.species_previous.data[index]) ||
            aliases(rate.previous)))) {
        local = {StatusCode::invalid_plan, kPredictorPlan};
      }
    }
    for (std::size_t index = 0U; index < passive_count && local; ++index) {
      const PredictorRateHistory rate =
          input.passive_scalar_nonadvective_rhs.data[index];
      if (aliases(input.passive_scalars_accepted.data[index]) ||
          aliases(rate.accepted) ||
          (second_order &&
           (aliases(input.passive_scalars_previous.data[index]) ||
            aliases(rate.previous)))) {
        local = {StatusCode::invalid_plan, kPredictorPlan};
      }
    }
  }
  consensus = collective_status(communicator, local, rank, size, failure,
                               failure);
  if (!consensus) return publish_failure(consensus);
  ++blocking_collectives;

  const std::size_t exchange_count = 1U + species_count + passive_count;
  if (exchange_count >
      std::numeric_limits<std::uint64_t>::max()) {
    mark_failure(
        failure, ThermophysicalPredictorFailureReason::transport_pass_overflow,
        ThermophysicalPredictorFailureField::none, UINT32_MAX,
        ThermophysicalAdmissibilityConstraint::none, rank, {-1, -1, -1}, 0U,
        false, false, false);
    return publish_failure(
        {StatusCode::numerical_failure, kPredictorNumerical});
  }
  const std::uint32_t substeps = 1U;
  std::uint32_t low_halo_exchanges = 0U;
  std::uint64_t low_transport_passes =
      static_cast<std::uint64_t>(exchange_count);
  ThermophysicalLowStateKind low_state =
      ThermophysicalLowStateKind::bdf_accepted_rate;
  double low_mass_flux_scale = 1.0;
  bool low_flux_is_high = false;
  bool low_flux_is_accepted = false;

  // The low endpoint is built in conserved form.  Keeping rho*q in the cold
  // workspace lets the scaled-Euler successor reuse the one accepted
  // transport divergence even when the BDF endpoint itself is inadmissible.
  bool bdf_candidate_finite = true;
  for (std::int32_t z = 0; z < end.z; ++z) {
    for (std::int32_t y = 0; y < end.y; ++y) {
      for (std::int32_t x = 0; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho_n =
            input.density_accepted.unchecked(cell, 0U);
        const double rho_previous =
            second_order ? input.density_previous.unchecked(cell, 0U) : 0.0;
        const double divergence =
            flux_divergence_density(*kernels_, input.mass_flux_accepted, cell);
        const double rho_bdf =
            (-input.bdf.a1 * rho_n -
             (second_order ? input.bdf.a2 * rho_previous : 0.0) -
             divergence) /
            input.bdf.a0;
        output.low_order_density_workspace.unchecked(cell, 0U) = rho_bdf;
        bdf_candidate_finite =
            bdf_candidate_finite && std::isfinite(rho_n) &&
            std::isfinite(rho_previous) && std::isfinite(divergence) &&
            std::isfinite(rho_bdf);
      }
    }
  }

  const auto build_bdf_quantity =
      [&](ConstFieldView accepted, ConstFieldView previous,
          PredictorRateHistory rate, FieldView candidate) noexcept {
        bool finite = true;
        for (std::int32_t z = 0; z < end.z; ++z) {
          for (std::int32_t y = 0; y < end.y; ++y) {
            for (std::int32_t x = 0; x < end.x; ++x) {
              const Int3 cell{x, y, z};
              const double q_n = accepted.unchecked(cell, 0U);
              const double q_previous =
                  second_order ? previous.unchecked(cell, 0U) : 0.0;
              const double divergence = first_order_upwind_divergence(
                  *kernels_, input.mass_flux_accepted, accepted, cell);
              const double rhs =
                  rate_is_zero(rate.accepted)
                      ? 0.0
                      : rate.accepted.unchecked(cell, 0U);
              const double rho_n =
                  input.density_accepted.unchecked(cell, 0U);
              const double rho_previous =
                  second_order
                      ? input.density_previous.unchecked(cell, 0U)
                      : 0.0;
              const double conserved_bdf =
                  (-input.bdf.a1 * rho_n * q_n -
                   (second_order ? input.bdf.a2 * rho_previous * q_previous
                                 : 0.0) -
                   (divergence - rhs)) /
                  input.bdf.a0;
              candidate.unchecked(cell, 0U) = conserved_bdf;
              finite =
                  finite && std::isfinite(q_n) &&
                  std::isfinite(q_previous) && std::isfinite(divergence) &&
                  std::isfinite(rhs) && std::isfinite(conserved_bdf);
            }
          }
        }
        return finite;
      };

  bdf_candidate_finite =
      build_bdf_quantity(input.enthalpy_accepted, input.enthalpy_previous,
                         input.enthalpy_nonadvective_rhs,
                         output.low_order_enthalpy_workspace) &&
      bdf_candidate_finite;
  for (std::size_t index = 0U; index < species_count; ++index) {
    bdf_candidate_finite =
        build_bdf_quantity(
            input.species_accepted.data[index],
            input.species_previous.data[index],
            input.species_nonadvective_rhs.data[index],
            output.low_order_independent_species.data[index]) &&
        bdf_candidate_finite;
  }
  for (std::size_t index = 0U; index < passive_count; ++index) {
    bdf_candidate_finite =
        build_bdf_quantity(
            input.passive_scalars_accepted.data[index],
            input.passive_scalars_previous.data[index],
            input.passive_scalar_nonadvective_rhs.data[index],
            output.low_order_passive_scalars.data[index]) &&
        bdf_candidate_finite;
  }

  const bool locally_bdf_admissible =
      bdf_candidate_finite &&
      tuple_admissible(as_const(output.low_order_density_workspace),
                       as_const(output.low_order_enthalpy_workspace),
                       output.low_order_independent_species,
                       output.low_order_passive_scalars, true);
  int local_bdf_route = locally_bdf_admissible ? 1 : 0;
  int global_bdf_route = 0;
  if (MPI_Allreduce(&local_bdf_route, &global_bdf_route, 1, MPI_INT,
                    MPI_MIN, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPredictorPlan};
  }
  ++blocking_collectives;

  const auto convert_bdf_to_intensive =
      [&]() noexcept {
        for (std::int32_t z = 0; z < end.z; ++z) {
          for (std::int32_t y = 0; y < end.y; ++y) {
            for (std::int32_t x = 0; x < end.x; ++x) {
              const Int3 cell{x, y, z};
              const double rho =
                  output.low_order_density_workspace.unchecked(cell, 0U);
              output.low_order_enthalpy_workspace.unchecked(cell, 0U) /=
                  rho;
              for (std::size_t index = 0U; index < species_count; ++index) {
                output.low_order_independent_species.data[index].unchecked(
                    cell, 0U) /= rho;
              }
              for (std::size_t index = 0U; index < passive_count; ++index) {
                output.low_order_passive_scalars.data[index].unchecked(
                    cell, 0U) /= rho;
              }
            }
          }
        }
      };

  constexpr double kAffineRoundoffSafety =
      1.0 - 64.0 * std::numeric_limits<double>::epsilon();
  constexpr double kRepresentationGuardFactor = 512.0;
  constexpr double kInvariantInteriorFraction = 0.1;
  const auto affine_admissibility_factor =
      [&](double anchor_margin, double trial_margin, double scale, bool strict,
          double& factor, double* selected_target) noexcept {
        if (!std::isfinite(anchor_margin) ||
            !std::isfinite(trial_margin) || !std::isfinite(scale) ||
            scale < 0.0 ||
            (strict ? !(anchor_margin > 0.0) : anchor_margin < 0.0)) {
          return false;
        }
        const double guard =
            kRepresentationGuardFactor *
            std::numeric_limits<double>::epsilon() * std::max(1.0, scale);
        if (!std::isfinite(guard) || !(guard > 0.0)) return false;
        // A merely representable positive endpoint leaves no room for the
        // downstream pressure correction.  Keep the low tuple inside the
        // same ten-percent admissible reserve used by the application-level
        // density globalization.  Zero species margins remain exactly zero.
        const double interior_margin =
            kInvariantInteriorFraction * anchor_margin;
        const double target_margin =
            std::min(anchor_margin, std::max(guard, interior_margin));
        if (selected_target != nullptr) *selected_target = target_margin;
        if (trial_margin >= target_margin) {
          factor = 1.0;
          return true;
        }
        const double denominator = anchor_margin - trial_margin;
        if (!std::isfinite(denominator) || !(denominator > 0.0)) {
          return false;
        }
        double candidate =
            (anchor_margin - target_margin) / denominator;
        if (!std::isfinite(candidate)) return false;
        candidate = std::clamp(candidate, 0.0, 1.0);
        if (candidate < 1.0 || strict) {
          candidate =
              std::nextafter(candidate * kAffineRoundoffSafety, 0.0);
        }
        if (!std::isfinite(candidate) || candidate < 0.0 ||
            candidate > 1.0) {
          return false;
        }
        factor = candidate;
        return true;
      };

  if (global_high_nonenthalpy != 0 &&
      std::isfinite(global_accepted_flux) && global_accepted_flux > 0.0) {
    low_state = ThermophysicalLowStateKind::implicit_upwind;
    low_transport_passes = 1U;
    low_mass_flux_scale = 1.0;
    low_flux_is_high = true;
    for (std::int32_t z = 0; z < end.z; ++z) {
      for (std::int32_t y = 0; y < end.y; ++y) {
        for (std::int32_t x = 0; x < end.x; ++x) {
          const Int3 cell{x, y, z};
          output.low_order_density_workspace.unchecked(cell, 0U) =
              output.density_workspace.unchecked(cell, 0U);
          for (std::size_t index = 0U; index < species_count; ++index) {
            output.low_order_independent_species.data[index].unchecked(
                cell, 0U) =
                output.independent_species.data[index].unchecked(cell, 0U);
          }
          for (std::size_t index = 0U; index < passive_count; ++index) {
            output.low_order_passive_scalars.data[index].unchecked(cell, 0U) =
                output.passive_scalars.data[index].unchecked(cell, 0U);
          }
        }
      }
    }
    if (slow_path.enthalpy_endpoint == nullptr) {
      mark_failure(
          failure,
          ThermophysicalPredictorFailureReason::implicit_enthalpy_endpoint,
          ThermophysicalPredictorFailureField::enthalpy, 0U,
          ThermophysicalAdmissibilityConstraint::none, rank,
          {-1, -1, -1}, 0U, false, false, false);
      local = {StatusCode::numerical_failure, kPredictorNumerical};
    }
    consensus = collective_status(communicator, local, rank, size, failure,
                                  failure);
    if (!consensus) return publish_failure(consensus);
    ++blocking_collectives;
    const ConservativeEnthalpyEndpointInput endpoint_input{
        input.bdf,
        input.time,
        as_const(output.density_workspace),
        input.enthalpy_accepted,
        as_const(output.enthalpy),
        input.mass_flux_accepted,
        output.accepted_advection_workspace,
        output.previous_advection_workspace,
        output.low_order_enthalpy_workspace,
        slow_path.resources};
    ++enthalpy_solve_calls;
    local = slow_path.enthalpy_endpoint->solve(endpoint_input,
                                               enthalpy_endpoint_result);
    if (!local) {
      mark_failure(
          failure,
          ThermophysicalPredictorFailureReason::implicit_enthalpy_endpoint,
          ThermophysicalPredictorFailureField::enthalpy, 0U,
          ThermophysicalAdmissibilityConstraint::none, rank,
          {-1, -1, -1}, 0U, false, false, false);
    }

    // Certify the converged implicit endpoint against the nonzero,
    // composition-dependent thermophysical interval.  The accepted
    // enthalpy is the common reference while density and composition remain
    // the already-admissible high values.  A communicator-wide alpha scales
    // the complete implicit defect; it does not clip a cell or damp mass
    // flux.
    double local_endpoint_alpha = local ? 1.0 : -1.0;
    for (std::int32_t z = 0; z < end.z && local; ++z) {
      for (std::int32_t y = 0; y < end.y && local; ++y) {
        for (std::int32_t x = 0; x < end.x && local; ++x) {
          const Int3 cell{x, y, z};
          const Int3 global_cell{cell.x + patch_begin_.x,
                                 cell.y + patch_begin_.y,
                                 cell.z + patch_begin_.z};
          const double rho =
              output.low_order_density_workspace.unchecked(cell, 0U);
          const double accepted = input.enthalpy_accepted.unchecked(cell, 0U);
          const double implicit =
              output.low_order_enthalpy_workspace.unchecked(cell, 0U);
          double dependent_density = rho;
          double lower = rho * dependent_enthalpy_minimum_;
          double upper = rho * dependent_enthalpy_maximum_;
          for (std::size_t index = 0U; index < species_count; ++index) {
            const double value =
                output.low_order_independent_species.data[index].unchecked(
                    cell, 0U);
            const double species_density = rho * value;
            dependent_density -= species_density;
            lower += species_density *
                     (species_enthalpy_minimum_[index] -
                      dependent_enthalpy_minimum_);
            upper += species_density *
                     (species_enthalpy_maximum_[index] -
                      dependent_enthalpy_maximum_);
          }
          const double accepted_conserved = rho * accepted;
          const double implicit_conserved = rho * implicit;
          ThermophysicalAdmissibilityConstraint constraint =
              ThermophysicalAdmissibilityConstraint::none;
          if (!std::isfinite(rho) || !(rho > 0.0) ||
              !std::isfinite(dependent_density) || dependent_density < 0.0 ||
              !std::isfinite(lower) || !std::isfinite(upper) ||
              !(lower <= upper) || !std::isfinite(accepted_conserved) ||
              !std::isfinite(implicit_conserved)) {
            mark_failure(
                failure,
                ThermophysicalPredictorFailureReason::
                    implicit_enthalpy_endpoint,
                ThermophysicalPredictorFailureField::enthalpy, 0U,
                constraint, rank, global_cell, 1U, true, true, false);
            local = {StatusCode::numerical_failure, kPredictorNumerical};
            local_endpoint_alpha = -1.0;
            break;
          }
          if (accepted_conserved < lower) {
            constraint =
                ThermophysicalAdmissibilityConstraint::enthalpy_lower;
          } else if (accepted_conserved > upper) {
            constraint =
                ThermophysicalAdmissibilityConstraint::enthalpy_upper;
          }
          if (constraint != ThermophysicalAdmissibilityConstraint::none) {
            mark_failure(
                failure,
                ThermophysicalPredictorFailureReason::
                    implicit_enthalpy_endpoint,
                ThermophysicalPredictorFailureField::enthalpy, 0U,
                constraint, rank, global_cell, 1U, true, true, true,
                accepted_conserved - lower, upper - accepted_conserved);
            if (failure.valid) {
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_density_current,
                  failure.density_current, rho);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_quantity_current,
                  failure.quantity_current, accepted);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_observed_value,
                  failure.observed_value, accepted_conserved);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_allowed_lower,
                  failure.allowed_lower, lower);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_allowed_upper,
                  failure.allowed_upper, upper);
            }
            local = {StatusCode::numerical_failure, kPredictorNumerical};
            local_endpoint_alpha = -1.0;
            break;
          }

          double candidate = 1.0;
          if (implicit_conserved < lower) {
            const double denominator =
                accepted_conserved - implicit_conserved;
            candidate = (accepted_conserved - lower) / denominator;
          } else if (implicit_conserved > upper) {
            const double denominator =
                implicit_conserved - accepted_conserved;
            candidate = (upper - accepted_conserved) / denominator;
          }
          if (!std::isfinite(candidate) || candidate < 0.0 ||
              candidate > 1.0) {
            mark_failure(
                failure,
                ThermophysicalPredictorFailureReason::
                    implicit_enthalpy_endpoint,
                ThermophysicalPredictorFailureField::enthalpy, 0U,
                implicit_conserved < lower
                    ? ThermophysicalAdmissibilityConstraint::enthalpy_lower
                    : ThermophysicalAdmissibilityConstraint::enthalpy_upper,
                rank, global_cell, 1U, true, true, false);
            local = {StatusCode::numerical_failure, kPredictorNumerical};
            local_endpoint_alpha = -1.0;
            break;
          }
          local_endpoint_alpha =
              std::min(local_endpoint_alpha, candidate);
        }
      }
    }
    double global_endpoint_alpha = -1.0;
    if (MPI_Allreduce(&local_endpoint_alpha, &global_endpoint_alpha, 1,
                      MPI_DOUBLE, MPI_MIN, communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kPredictorPlan};
    }
    ++blocking_collectives;
    if (local && std::isfinite(global_endpoint_alpha) &&
        global_endpoint_alpha >= 0.0 && global_endpoint_alpha <= 1.0) {
      if (global_endpoint_alpha < 1.0) {
        constexpr double kEndpointSafety =
            1.0 - 64.0 * std::numeric_limits<double>::epsilon();
        global_endpoint_alpha = std::nextafter(
            global_endpoint_alpha * kEndpointSafety, 0.0);
        low_state =
            ThermophysicalLowStateKind::implicit_upwind_source_limited;
        for (std::int32_t z = 0; z < end.z; ++z) {
          for (std::int32_t y = 0; y < end.y; ++y) {
            for (std::int32_t x = 0; x < end.x; ++x) {
              const Int3 cell{x, y, z};
              const double accepted =
                  input.enthalpy_accepted.unchecked(cell, 0U);
              const double implicit =
                  output.low_order_enthalpy_workspace.unchecked(cell, 0U);
              output.low_order_enthalpy_workspace.unchecked(cell, 0U) =
                  accepted + global_endpoint_alpha * (implicit - accepted);
            }
          }
        }
      }
      enthalpy_endpoint_alpha = global_endpoint_alpha;
      if (!tuple_admissible(
              as_const(output.low_order_density_workspace),
              as_const(output.low_order_enthalpy_workspace),
              output.low_order_independent_species,
              output.low_order_passive_scalars)) {
        capture_tuple_failure(
            as_const(output.low_order_density_workspace),
            as_const(output.low_order_enthalpy_workspace),
            output.low_order_independent_species,
            output.low_order_passive_scalars,
            ThermophysicalPredictorFailureReason::implicit_enthalpy_endpoint,
            1U);
        local = {StatusCode::numerical_failure, kPredictorNumerical};
      }
    }
    consensus = collective_status(communicator, local, rank, size, failure,
                                  failure);
    if (!consensus) return publish_failure(consensus);
    ++blocking_collectives;
  } else if (global_bdf_route != 0) {
    convert_bdf_to_intensive();
    low_flux_is_accepted = true;
  } else {
    // Build the zero-incoming BDF/source base in conserved form.  It must be
    // admissible before any outgoing transport is admitted; incoming upwind
    // donor tuples are then convex-cone additions.
    const auto donor_intensive_admissible = [&](Int3 donor) noexcept {
      const double h = input.enthalpy_accepted.unchecked(donor, 0U);
      double dependent = 1.0;
      double lower = dependent_enthalpy_minimum_;
      double upper = dependent_enthalpy_maximum_;
      for (std::size_t index = 0U; index < species_count; ++index) {
        const double value =
            input.species_accepted.data[index].unchecked(donor, 0U);
        if (!std::isfinite(value) || value < 0.0) return false;
        dependent -= value;
        lower += value * (species_enthalpy_minimum_[index] -
                          dependent_enthalpy_minimum_);
        upper += value * (species_enthalpy_maximum_[index] -
                          dependent_enthalpy_maximum_);
      }
      if (!std::isfinite(h) || !std::isfinite(dependent) ||
          dependent < 0.0 || !std::isfinite(lower) ||
          !std::isfinite(upper) || h < lower || h > upper) {
        return false;
      }
      for (std::size_t index = 0U; index < passive_count; ++index) {
        if (!std::isfinite(input.passive_scalars_accepted.data[index]
                               .unchecked(donor, 0U))) {
          return false;
        }
      }
      return true;
    };
    const auto source_free_base_quantity = [&](ConstFieldView accepted,
                                               ConstFieldView previous,
                                               Int3 cell) noexcept {
      const double rho_n = input.density_accepted.unchecked(cell, 0U);
      const double rho_nm1 =
          second_order ? input.density_previous.unchecked(cell, 0U) : 0.0;
      const double q_n = accepted.unchecked(cell, 0U);
      const double q_nm1 = second_order ? previous.unchecked(cell, 0U) : 0.0;
      return (-input.bdf.a1 * rho_n * q_n -
              input.bdf.a2 * rho_nm1 * q_nm1) /
             input.bdf.a0;
    };
    const auto full_source_base_quantity =
        [&](ConstFieldView accepted, ConstFieldView previous,
            PredictorRateHistory rate, Int3 cell) noexcept {
      const double rho_n = input.density_accepted.unchecked(cell, 0U);
      const double rho_nm1 =
          second_order ? input.density_previous.unchecked(cell, 0U) : 0.0;
      const double q_n = accepted.unchecked(cell, 0U);
      const double q_nm1 = second_order ? previous.unchecked(cell, 0U) : 0.0;
      const double source = rate_is_zero(rate.accepted)
                                ? 0.0
                                : rate.accepted.unchecked(cell, 0U);
      return (-input.bdf.a1 * rho_n * q_n -
              input.bdf.a2 * rho_nm1 * q_nm1 + source) /
             input.bdf.a0;
    };
    bool local_base_valid = true;
    for (std::int32_t z = 0; z < end.z; ++z) {
      for (std::int32_t y = 0; y < end.y; ++y) {
        for (std::int32_t x = 0; x < end.x; ++x) {
          const Int3 cell{x, y, z};
          const double rho_n = input.density_accepted.unchecked(cell, 0U);
          const double rho_nm1 =
              second_order ? input.density_previous.unchecked(cell, 0U) : 0.0;
          const double base_rho =
              (-input.bdf.a1 * rho_n - input.bdf.a2 * rho_nm1) /
              input.bdf.a0;
          output.low_order_density_workspace.unchecked(cell, 0U) = base_rho;
          local_base_valid =
              local_base_valid && std::isfinite(rho_n) && rho_n > 0.0 &&
              donor_intensive_admissible(cell);
          output.low_order_enthalpy_workspace.unchecked(cell, 0U) =
              source_free_base_quantity(input.enthalpy_accepted,
                                        input.enthalpy_previous, cell);
          for (std::size_t index = 0U; index < species_count; ++index) {
            output.low_order_independent_species.data[index].unchecked(
                cell, 0U) =
                source_free_base_quantity(input.species_accepted.data[index],
                                          input.species_previous.data[index],
                                          cell);
          }
          for (std::size_t index = 0U; index < passive_count; ++index) {
            output.low_order_passive_scalars.data[index].unchecked(cell, 0U) =
                source_free_base_quantity(
                    input.passive_scalars_accepted.data[index],
                    input.passive_scalars_previous.data[index], cell);
          }
        }
      }
    }
    local_base_valid = local_base_valid && tuple_admissible(
        as_const(output.low_order_density_workspace),
        as_const(output.low_order_enthalpy_workspace),
        output.low_order_independent_species,
        output.low_order_passive_scalars, true);
    if (!local_base_valid) {
      capture_tuple_failure(
          as_const(output.low_order_density_workspace),
          as_const(output.low_order_enthalpy_workspace),
          output.low_order_independent_species,
          output.low_order_passive_scalars,
          ThermophysicalPredictorFailureReason::
              low_bdf_source_base_admissibility,
          1U, true);
      local = {StatusCode::numerical_failure, kPredictorNumerical};
    }

    double local_source_endpoint_alpha = local_base_valid ? 1.0 : -1.0;
    for (std::int32_t z = 0; z < end.z && local; ++z) {
      for (std::int32_t y = 0; y < end.y && local; ++y) {
        for (std::int32_t x = 0; x < end.x; ++x) {
          const Int3 cell{x, y, z};
          const double base_rho =
              output.low_order_density_workspace.unchecked(cell, 0U);
          bool valid = std::isfinite(base_rho) && base_rho > 0.0;
          double cell_factor = 1.0;
          const auto constrain_source = [&](double anchor_margin,
                                            double trial_margin,
                                            bool strict) noexcept {
            double candidate = 1.0;
            valid = valid && affine_admissibility_factor(
                                 anchor_margin, trial_margin,
                                 std::max({1.0, std::abs(anchor_margin),
                                           std::abs(trial_margin)}),
                                 strict, candidate, nullptr);
            cell_factor = std::min(cell_factor, candidate);
          };

          constrain_source(base_rho, base_rho, true);
          double base_dependent = base_rho;
          double full_dependent = base_rho;
          double base_lower = base_rho * dependent_enthalpy_minimum_;
          double full_lower = base_lower;
          double base_upper = base_rho * dependent_enthalpy_maximum_;
          double full_upper = base_upper;
          for (std::size_t index = 0U; index < species_count; ++index) {
            const double base_species =
                output.low_order_independent_species.data[index].unchecked(
                    cell, 0U);
            const double full_species = full_source_base_quantity(
                input.species_accepted.data[index],
                input.species_previous.data[index],
                input.species_nonadvective_rhs.data[index], cell);
            constrain_source(base_species, full_species, false);
            base_dependent -= base_species;
            full_dependent -= full_species;
            const double lower_delta = species_enthalpy_minimum_[index] -
                                       dependent_enthalpy_minimum_;
            const double upper_delta = species_enthalpy_maximum_[index] -
                                       dependent_enthalpy_maximum_;
            base_lower += base_species * lower_delta;
            full_lower += full_species * lower_delta;
            base_upper += base_species * upper_delta;
            full_upper += full_species * upper_delta;
          }
          constrain_source(base_dependent, full_dependent, false);

          const double base_enthalpy =
              output.low_order_enthalpy_workspace.unchecked(cell, 0U);
          const double full_enthalpy = full_source_base_quantity(
              input.enthalpy_accepted, input.enthalpy_previous,
              input.enthalpy_nonadvective_rhs, cell);
          constrain_source(base_enthalpy - base_lower,
                           full_enthalpy - full_lower, false);
          constrain_source(base_upper - base_enthalpy,
                           full_upper - full_enthalpy, false);

          for (std::size_t index = 0U; index < passive_count; ++index) {
            const double base_passive =
                output.low_order_passive_scalars.data[index].unchecked(cell,
                                                                       0U);
            const double full_passive = full_source_base_quantity(
                input.passive_scalars_accepted.data[index],
                input.passive_scalars_previous.data[index],
                input.passive_scalar_nonadvective_rhs.data[index], cell);
            valid = valid && std::isfinite(base_passive) &&
                    std::isfinite(full_passive);
          }
          if (!valid || !std::isfinite(cell_factor) || cell_factor < 0.0 ||
              cell_factor > 1.0) {
            mark_failure(
                failure,
                ThermophysicalPredictorFailureReason::low_tuple_admissibility,
                ThermophysicalPredictorFailureField::nonadvective_rhs,
                UINT32_MAX, ThermophysicalAdmissibilityConstraint::none, rank,
                {cell.x + patch_begin_.x, cell.y + patch_begin_.y,
                 cell.z + patch_begin_.z},
                1U, true, true, false);
            local = {StatusCode::numerical_failure, kPredictorNumerical};
            local_source_endpoint_alpha = -1.0;
            break;
          }
          local_source_endpoint_alpha =
              std::min(local_source_endpoint_alpha, cell_factor);
        }
      }
    }

    double global_source_endpoint_alpha = -1.0;
    if (MPI_Allreduce(&local_source_endpoint_alpha,
                      &global_source_endpoint_alpha, 1, MPI_DOUBLE, MPI_MIN,
                      communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kPredictorPlan};
    }
    ++blocking_collectives;
    if (!std::isfinite(global_source_endpoint_alpha) ||
        global_source_endpoint_alpha < 0.0 ||
        global_source_endpoint_alpha > 1.0) {
      consensus = collective_status(communicator, local, rank, size, failure,
                                    failure);
      ++blocking_collectives;
      if (consensus) {
        mark_failure(
            failure,
            ThermophysicalPredictorFailureReason::low_tuple_admissibility,
            ThermophysicalPredictorFailureField::nonadvective_rhs,
            UINT32_MAX, ThermophysicalAdmissibilityConstraint::none, rank,
            {-1, -1, -1}, 1U, false, true, false);
        consensus = {StatusCode::numerical_failure, kPredictorNumerical};
      }
      return publish_failure(consensus);
    }
    source_endpoint_alpha = global_source_endpoint_alpha;

    const auto publish_source_quantity =
        [&](ConstFieldView accepted, ConstFieldView previous,
            PredictorRateHistory rate, FieldView destination) noexcept {
      for (std::int32_t z = 0; z < end.z; ++z) {
        for (std::int32_t y = 0; y < end.y; ++y) {
          for (std::int32_t x = 0; x < end.x; ++x) {
            const Int3 cell{x, y, z};
            const double full = full_source_base_quantity(
                accepted, previous, rate, cell);
            if (source_endpoint_alpha == 1.0) {
              destination.unchecked(cell, 0U) = full;
            } else if (source_endpoint_alpha != 0.0) {
              const double anchor = destination.unchecked(cell, 0U);
              destination.unchecked(cell, 0U) =
                  anchor + source_endpoint_alpha * (full - anchor);
            }
          }
        }
      }
    };
    publish_source_quantity(input.enthalpy_accepted,
                            input.enthalpy_previous,
                            input.enthalpy_nonadvective_rhs,
                            output.low_order_enthalpy_workspace);
    for (std::size_t index = 0U; index < species_count; ++index) {
      publish_source_quantity(
          input.species_accepted.data[index],
          input.species_previous.data[index],
          input.species_nonadvective_rhs.data[index],
          output.low_order_independent_species.data[index]);
    }
    for (std::size_t index = 0U; index < passive_count; ++index) {
      publish_source_quantity(
          input.passive_scalars_accepted.data[index],
          input.passive_scalars_previous.data[index],
          input.passive_scalar_nonadvective_rhs.data[index],
          output.low_order_passive_scalars.data[index]);
    }
    if (!tuple_admissible(
            as_const(output.low_order_density_workspace),
            as_const(output.low_order_enthalpy_workspace),
            output.low_order_independent_species,
            output.low_order_passive_scalars, true)) {
      capture_tuple_failure(
          as_const(output.low_order_density_workspace),
          as_const(output.low_order_enthalpy_workspace),
          output.low_order_independent_species,
          output.low_order_passive_scalars,
          ThermophysicalPredictorFailureReason::low_tuple_admissibility, 1U,
          true);
      local = {StatusCode::numerical_failure, kPredictorNumerical};
    }

    double local_minimum_factor = 1.0;
    for (std::int32_t z = 0; z < end.z && local; ++z) {
      for (std::int32_t y = 0; y < end.y && local; ++y) {
        for (std::int32_t x = 0; x < end.x; ++x) {
          const Int3 cell{x, y, z};
          const double volume = detail::cell_volume(*kernels_, cell);
          const double outgoing = outgoing_mass_flux(
              input.mass_flux_accepted, cell);
          const double amount = outgoing / (input.bdf.a0 * volume);
          const double rho_n = input.density_accepted.unchecked(cell, 0U);
          const double base_rho =
              output.low_order_density_workspace.unchecked(cell, 0U);
          bool valid = std::isfinite(volume) && volume > 0.0 &&
                       std::isfinite(outgoing) && outgoing >= 0.0 &&
                       std::isfinite(amount) && amount >= 0.0 &&
                       std::isfinite(rho_n) && rho_n > 0.0;
          double factor = 1.0;
          const auto constrain_factor = [&](double base_margin,
                                            double full_margin,
                                            bool strict) noexcept {
            double candidate = 1.0;
            valid = valid && affine_admissibility_factor(
                                 base_margin, full_margin,
                                 std::max({1.0, std::abs(base_margin),
                                           std::abs(full_margin)}),
                                 strict, candidate, nullptr);
            factor = std::min(factor, candidate);
          };
          constrain_factor(base_rho, base_rho - amount, true);
          double base_dependent = base_rho;
          double full_dependent = base_rho - amount;
          double base_lower = base_rho * dependent_enthalpy_minimum_;
          double full_lower =
              (base_rho - amount) * dependent_enthalpy_minimum_;
          double base_upper = base_rho * dependent_enthalpy_maximum_;
          double full_upper =
              (base_rho - amount) * dependent_enthalpy_maximum_;
          double accepted_sum = 0.0;
          for (std::size_t index = 0U; index < species_count; ++index) {
            const double y_n =
                input.species_accepted.data[index].unchecked(cell, 0U);
            const double base_species =
                output.low_order_independent_species.data[index].unchecked(
                    cell, 0U);
            const double full_species = base_species - amount * y_n;
            valid = valid && std::isfinite(y_n) && y_n >= 0.0;
            accepted_sum += y_n;
            constrain_factor(base_species, full_species, false);
            base_dependent -= base_species;
            full_dependent -= full_species;
            const double lower_delta = species_enthalpy_minimum_[index] -
                                       dependent_enthalpy_minimum_;
            const double upper_delta = species_enthalpy_maximum_[index] -
                                       dependent_enthalpy_maximum_;
            base_lower += base_species * lower_delta;
            full_lower += full_species * lower_delta;
            base_upper += base_species * upper_delta;
            full_upper += full_species * upper_delta;
          }
          valid = valid && std::isfinite(accepted_sum) &&
                  accepted_sum <= 1.0;
          constrain_factor(base_dependent, full_dependent, false);
          const double h_n = input.enthalpy_accepted.unchecked(cell, 0U);
          const double base_rho_h =
              output.low_order_enthalpy_workspace.unchecked(cell, 0U);
          const double full_rho_h = base_rho_h - amount * h_n;
          valid = valid && std::isfinite(h_n);
          constrain_factor(base_rho_h - base_lower,
                           full_rho_h - full_lower, false);
          constrain_factor(base_upper - base_rho_h,
                           full_upper - full_rho_h, false);
          for (std::size_t index = 0U; index < passive_count; ++index) {
            const double q_n = input.passive_scalars_accepted.data[index]
                                   .unchecked(cell, 0U);
            const double base = output.low_order_passive_scalars.data[index]
                                    .unchecked(cell, 0U);
            valid = valid && std::isfinite(q_n) && std::isfinite(base) &&
                    std::isfinite(base - amount * q_n);
          }
          if (!valid || !std::isfinite(factor) || factor < 0.0 ||
              factor > 1.0) {
            mark_failure(
                failure,
                ThermophysicalPredictorFailureReason::low_tuple_admissibility,
                ThermophysicalPredictorFailureField::mass_flux, UINT32_MAX,
                ThermophysicalAdmissibilityConstraint::none, rank,
                {cell.x + patch_begin_.x, cell.y + patch_begin_.y,
                 cell.z + patch_begin_.z},
                1U, true, true, false);
            local = {StatusCode::numerical_failure, kPredictorNumerical};
            break;
          }
          output.low_order_density_workspace.unchecked(cell, 0U) = factor;
          local_minimum_factor = std::min(local_minimum_factor, factor);
        }
      }
    }
    consensus = collective_status(communicator, local, rank, size, failure,
                                  failure);
    if (!consensus) return publish_failure(consensus);
    ++blocking_collectives;
    if (slow_path.halo == nullptr || !slow_path.halo->ready() ||
        slow_path.halo_stage == 0U ||
        output.low_order_density_workspace.ghosts.x < 1 ||
        output.low_order_density_workspace.ghosts.y < 1 ||
        output.low_order_density_workspace.ghosts.z < 1) {
      return publish_failure(
          {StatusCode::invalid_plan, kPredictorPlan});
    }
    FieldView donor_factor = output.low_order_density_workspace;
    HaloTicket donor_ticket;
    local = slow_path.halo->begin(slow_path.halo_stage,
                                  {&donor_factor, 1U}, donor_ticket);
    if (local)
      local = slow_path.halo->finish(donor_ticket, {&donor_factor, 1U});
    if (!local) return publish_failure(local);
    output.low_order_density_workspace = donor_factor;
    low_halo_exchanges = 1U;

    double global_minimum_factor = 1.0;
    if (MPI_Allreduce(&local_minimum_factor, &global_minimum_factor, 1,
                      MPI_DOUBLE, MPI_MIN, communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kPredictorPlan};
    }
    ++blocking_collectives;
    if (!std::isfinite(global_minimum_factor) ||
        global_minimum_factor < 0.0 || global_minimum_factor > 1.0) {
      return publish_failure(
          {StatusCode::numerical_failure, kPredictorNumerical});
    }
    low_mass_flux_scale = global_minimum_factor;
    low_state = source_endpoint_alpha < 1.0
                    ? ThermophysicalLowStateKind::
                          bdf_local_donor_flux_source_limited
                    : ThermophysicalLowStateKind::bdf_local_donor_flux;
    bdf_endpoint_alpha = 1.0;

    const auto physical_boundary = [&](CartesianAxis axis,
                                       bool high) noexcept {
      const CartesianFace face =
          axis == CartesianAxis::x
              ? (high ? CartesianFace::x_max : CartesianFace::x_min)
              : (axis == CartesianAxis::y
                     ? (high ? CartesianFace::y_max : CartesianFace::y_min)
                     : (high ? CartesianFace::z_max
                             : CartesianFace::z_min));
      const BoundaryFacePlan* plan = nullptr;
      return boundary_->face(face, plan) && plan != nullptr &&
             plan->local_owner && !plan->periodic;
    };
    const ConstFaceFieldView accepted_fields[]{input.mass_flux_accepted.x,
                                                input.mass_flux_accepted.y,
                                                input.mass_flux_accepted.z};
    const FaceFieldView limited_fields[]{output.paired_mass_flux.x,
                                         output.paired_mass_flux.y,
                                         output.paired_mass_flux.z};
    for (std::size_t axis_index = 0U; axis_index < 3U && local;
         ++axis_index) {
      const auto axis = static_cast<CartesianAxis>(axis_index);
      const FaceFieldView limited = limited_fields[axis_index];
      for (std::int32_t z = 0; z < limited.extents.z && local; ++z) {
        for (std::int32_t y = 0; y < limited.extents.y && local; ++y) {
          for (std::int32_t x = 0; x < limited.extents.x; ++x) {
            const Int3 face{x, y, z};
            const std::int32_t normal =
                axis == CartesianAxis::x
                    ? x
                    : (axis == CartesianAxis::y ? y : z);
            const std::int32_t extent =
                axis == CartesianAxis::x
                    ? end.x
                    : (axis == CartesianAxis::y ? end.y : end.z);
            Int3 negative = face;
            if (axis == CartesianAxis::x)
              --negative.x;
            else if (axis == CartesianAxis::y)
              --negative.y;
            else
              --negative.z;
            const Int3 positive = face;
            const double accepted = accepted_fields[axis_index].unchecked(face);
            const bool donor_is_negative = accepted >= 0.0;
            const bool donor_outside =
                (normal == 0 && donor_is_negative) ||
                (normal == extent && !donor_is_negative);
            const bool physical =
                donor_outside &&
                physical_boundary(axis, normal == extent);
            const Int3 donor = donor_is_negative ? negative : positive;
            if (physical && !donor_intensive_admissible(donor)) {
              mark_failure(
                  failure,
                  ThermophysicalPredictorFailureReason::low_tuple_admissibility,
                  ThermophysicalPredictorFailureField::mass_flux,
                  UINT32_MAX,
                  ThermophysicalAdmissibilityConstraint::none, rank,
                  {-1, -1, -1}, 1U, false, true, false);
              local = {StatusCode::numerical_failure, kPredictorNumerical};
              break;
            }
            const double factor =
                physical ? 1.0
                         : donor_factor.unchecked(donor, 0U);
            const double value = accepted * factor;
            if (!std::isfinite(accepted) || !std::isfinite(factor) ||
                factor < 0.0 || factor > 1.0 || !std::isfinite(value)) {
              mark_failure(
                  failure,
                  ThermophysicalPredictorFailureReason::low_divergence,
                  ThermophysicalPredictorFailureField::mass_flux,
                  UINT32_MAX, ThermophysicalAdmissibilityConstraint::none,
                  rank, {-1, -1, -1}, 1U, false, true, false);
              local = {StatusCode::numerical_failure, kPredictorNumerical};
              break;
            }
            limited.unchecked(face) = value;
          }
        }
      }
    }
    if (local && slow_path.immersed_interface != nullptr)
      local = slow_path.immersed_interface->zero_interface_flux(
          output.paired_mass_flux);

    const ConstFaceFluxView low_flux = as_const(output.paired_mass_flux);
    for (std::int32_t z = 0; z < end.z && local; ++z) {
      for (std::int32_t y = 0; y < end.y && local; ++y) {
        for (std::int32_t x = 0; x < end.x; ++x) {
          const Int3 cell{x, y, z};
          const double rho_n = input.density_accepted.unchecked(cell, 0U);
          const double rho_nm1 =
              second_order ? input.density_previous.unchecked(cell, 0U) : 0.0;
          output.low_order_density_workspace.unchecked(cell, 0U) =
              (-input.bdf.a1 * rho_n - input.bdf.a2 * rho_nm1 -
               flux_divergence_density(*kernels_, low_flux, cell)) /
              input.bdf.a0;
        }
      }
    }
    const auto build_limited_quantity =
        [&](ConstFieldView accepted, ConstFieldView previous,
            PredictorRateHistory rate, FieldView destination) noexcept {
          for (std::int32_t z = 0; z < end.z; ++z) {
            for (std::int32_t y = 0; y < end.y; ++y) {
              for (std::int32_t x = 0; x < end.x; ++x) {
                const Int3 cell{x, y, z};
                const double rho_n =
                    input.density_accepted.unchecked(cell, 0U);
                const double rho_nm1 =
                    second_order
                        ? input.density_previous.unchecked(cell, 0U)
                        : 0.0;
                const double q_nm1 =
                    second_order ? previous.unchecked(cell, 0U) : 0.0;
                const double source = rate_is_zero(rate.accepted)
                                          ? 0.0
                                          : rate.accepted.unchecked(cell, 0U);
                const double limited_source =
                    source_endpoint_alpha == 1.0
                        ? source
                        : source_endpoint_alpha * source;
                destination.unchecked(cell, 0U) =
                    (-input.bdf.a1 * rho_n * accepted.unchecked(cell, 0U) -
                     input.bdf.a2 * rho_nm1 * q_nm1 -
                     first_order_upwind_divergence(*kernels_, low_flux,
                                                   accepted, cell) +
                     limited_source) /
                    input.bdf.a0;
              }
            }
          }
        };
    if (local) {
      build_limited_quantity(input.enthalpy_accepted,
                             input.enthalpy_previous,
                             input.enthalpy_nonadvective_rhs,
                             output.low_order_enthalpy_workspace);
      for (std::size_t index = 0U; index < species_count; ++index)
        build_limited_quantity(
            input.species_accepted.data[index],
            input.species_previous.data[index],
            input.species_nonadvective_rhs.data[index],
            output.low_order_independent_species.data[index]);
      for (std::size_t index = 0U; index < passive_count; ++index)
        build_limited_quantity(
            input.passive_scalars_accepted.data[index],
            input.passive_scalars_previous.data[index],
            input.passive_scalar_nonadvective_rhs.data[index],
            output.low_order_passive_scalars.data[index]);
      if (!tuple_admissible(
              as_const(output.low_order_density_workspace),
              as_const(output.low_order_enthalpy_workspace),
              output.low_order_independent_species,
              output.low_order_passive_scalars, true)) {
        capture_tuple_failure(
            as_const(output.low_order_density_workspace),
            as_const(output.low_order_enthalpy_workspace),
            output.low_order_independent_species,
            output.low_order_passive_scalars,
            ThermophysicalPredictorFailureReason::low_tuple_admissibility,
            1U, true);
        local = {StatusCode::numerical_failure, kPredictorNumerical};
      } else {
        convert_bdf_to_intensive();
      }
    }
  }
  consensus = collective_status(communicator, local, rank, size, failure,
                               failure);
  if (!consensus) return publish_failure(consensus);
  ++blocking_collectives;

  double local_theta = 1.0;
  ThermophysicalPredictorDiagnostics local_diagnostics;
  local_diagnostics.theta = 1.0;
  local_diagnostics.mass_flux_scale = low_mass_flux_scale;
  local_diagnostics.low_order_substeps = substeps;
  local_diagnostics.low_order_halo_exchanges = low_halo_exchanges;
  local_diagnostics.low_state = low_state;
  local_diagnostics.enthalpy_endpoint = enthalpy_endpoint_result;
  local_diagnostics.enthalpy_endpoint_alpha = enthalpy_endpoint_alpha;
  local_diagnostics.bdf_endpoint_alpha = bdf_endpoint_alpha;
  local_diagnostics.source_endpoint_alpha = source_endpoint_alpha;
  local_diagnostics.enthalpy_solve_calls = enthalpy_solve_calls;
  const auto constrain = [&](double low_margin, double high_margin,
                             double scale, bool strict,
                             ThermophysicalAdmissibilityConstraint constraint,
                             Int3 cell) noexcept {
    double candidate = 1.0;
    double target_margin = 0.0;
    if (!affine_admissibility_factor(low_margin, high_margin, scale, strict,
                                     candidate, &target_margin)) {
      candidate = -1.0;
    }
    if (candidate < local_theta) {
      local_theta = candidate;
      local_diagnostics.theta = candidate;
      // Evidence margins are measured from the guarded interior target.  This
      // preserves the established nonnegative/negative limiter certificate
      // even when a physically positive high margin lies inside the reserved
      // binary64 representation band.
      local_diagnostics.low_margin = low_margin - target_margin;
      local_diagnostics.high_margin = high_margin - target_margin;
      local_diagnostics.constraint = constraint;
      local_diagnostics.limiting_cell = {
          cell.x + patch_begin_.x,
          cell.y + patch_begin_.y,
          cell.z + patch_begin_.z};
      local_diagnostics.limiting_rank = rank;
    }
  };

  for (std::int32_t z = 0; z < end.z && local; ++z) {
    for (std::int32_t y = 0; y < end.y && local; ++y) {
      for (std::int32_t x = 0; x < end.x && local; ++x) {
        const Int3 cell{x, y, z};
        const double rho_low =
            output.low_order_density_workspace.unchecked(cell, 0U);
        const double rho_high = output.density_workspace.unchecked(cell, 0U);
        // Density is the sole strict inequality.  A zero high endpoint must
        // therefore select an inward theta rather than masquerade as an
        // admissible zero-margin constraint.
        const double density_high_margin =
            rho_high > 0.0
                ? rho_high
                : (rho_high < 0.0
                       ? rho_high
                       : -std::numeric_limits<double>::epsilon() *
                             std::max(1.0, std::abs(rho_low)));
        constrain(rho_low, density_high_margin,
                  std::max({1.0, std::abs(rho_low), std::abs(rho_high)}),
                  true,
                  ThermophysicalAdmissibilityConstraint::density, cell);
        double dependent_low = rho_low;
        double dependent_high = rho_high;
        double lower_low = rho_low * dependent_enthalpy_minimum_;
        double lower_high = rho_high * dependent_enthalpy_minimum_;
        double upper_low = rho_low * dependent_enthalpy_maximum_;
        double upper_high = rho_high * dependent_enthalpy_maximum_;
        for (std::size_t index = 0U; index < species_count; ++index) {
          const double species_low =
              rho_low * output.low_order_independent_species.data[index]
                            .unchecked(cell, 0U);
          const double high_species_value =
              output.independent_species.data[index].unchecked(cell, 0U);
          // predict_high_local stores the conserved numerator directly when
          // rho_high is exactly zero; such a state can only reach this slow
          // path and is never published in intensive form.
          const double species_high =
              rho_high == 0.0 ? high_species_value
                              : rho_high * high_species_value;
          constrain(species_low, species_high,
                    std::max({1.0, std::abs(rho_low), std::abs(rho_high),
                              std::abs(species_low), std::abs(species_high)}),
                    false,
                    ThermophysicalAdmissibilityConstraint::independent_species,
                    cell);
          dependent_low -= species_low;
          dependent_high -= species_high;
          lower_low += species_low *
                       (species_enthalpy_minimum_[index] -
                        dependent_enthalpy_minimum_);
          lower_high += species_high *
                        (species_enthalpy_minimum_[index] -
                         dependent_enthalpy_minimum_);
          upper_low += species_low *
                       (species_enthalpy_maximum_[index] -
                        dependent_enthalpy_maximum_);
          upper_high += species_high *
                        (species_enthalpy_maximum_[index] -
                         dependent_enthalpy_maximum_);
        }
        constrain(dependent_low, dependent_high,
                  std::max({1.0, std::abs(rho_low), std::abs(rho_high),
                            std::abs(dependent_low),
                            std::abs(dependent_high)}),
                  false,
                  ThermophysicalAdmissibilityConstraint::dependent_species,
                  cell);
        const double density_enthalpy_low =
            rho_low *
            output.low_order_enthalpy_workspace.unchecked(cell, 0U);
        const double high_enthalpy_value =
            output.enthalpy.unchecked(cell, 0U);
        const double density_enthalpy_high =
            rho_high == 0.0 ? high_enthalpy_value
                            : rho_high * high_enthalpy_value;
        constrain(density_enthalpy_low - lower_low,
                  density_enthalpy_high - lower_high,
                  std::max({1.0, std::abs(density_enthalpy_low),
                            std::abs(density_enthalpy_high),
                            std::abs(lower_low), std::abs(lower_high)}),
                  false,
                  ThermophysicalAdmissibilityConstraint::enthalpy_lower,
                  cell);
        constrain(upper_low - density_enthalpy_low,
                  upper_high - density_enthalpy_high,
                  std::max({1.0, std::abs(density_enthalpy_low),
                            std::abs(density_enthalpy_high),
                            std::abs(upper_low), std::abs(upper_high)}),
                  false,
                  ThermophysicalAdmissibilityConstraint::enthalpy_upper,
                  cell);
      }
    }
  }

  struct ThetaRank {
    double theta;
    int rank;
  } local_theta_rank{local_theta, rank}, global_theta_rank{};
  if (MPI_Allreduce(&local_theta_rank, &global_theta_rank, 1, MPI_DOUBLE_INT,
                    MPI_MINLOC, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPredictorPlan};
  }
  ++blocking_collectives;
  if (!std::isfinite(global_theta_rank.theta) ||
      global_theta_rank.theta < 0.0 || global_theta_rank.theta > 1.0 ||
      global_theta_rank.rank < 0 || global_theta_rank.rank >= size) {
    const int failure_root = global_theta_rank.rank >= 0 &&
                                     global_theta_rank.rank < size
                                 ? global_theta_rank.rank
                                 : 0;
    if (rank == failure_root) {
      const auto constraint = local_diagnostics.constraint;
      const bool has_cell =
          local_diagnostics.limiting_rank >= 0 && constraint !=
                                                     ThermophysicalAdmissibilityConstraint::none;
      mark_failure(
          failure, ThermophysicalPredictorFailureReason::limiter_theta,
          failure_field_for_constraint(constraint), UINT32_MAX, constraint,
          failure_root,
          has_cell ? local_diagnostics.limiting_cell : Int3{-1, -1, -1}, 0U,
          has_cell, false,
          std::isfinite(local_diagnostics.low_margin) &&
              std::isfinite(local_diagnostics.high_margin),
          local_diagnostics.low_margin, local_diagnostics.high_margin);
      if (failure.valid) {
        record_scalar(failure,
                      ThermophysicalPredictorFailure::scalar_observed_value,
                      failure.observed_value, global_theta_rank.theta);
        record_scalar(failure,
                      ThermophysicalPredictorFailure::scalar_allowed_lower,
                      failure.allowed_lower, 0.0);
        record_scalar(failure,
                      ThermophysicalPredictorFailure::scalar_allowed_upper,
                      failure.allowed_upper, 1.0);
      }
    }
    const Status failure_status =
        broadcast_failure(communicator, failure_root, rank, failure);
    if (!failure_status) return failure_status;
    return publish_failure(
        {StatusCode::numerical_failure, kPredictorNumerical});
  }
  std::array<double, 3U> selected_reals{};
  std::array<std::uint64_t, 9U> selected_integers{};
  if (rank == global_theta_rank.rank) {
    selected_reals = {global_theta_rank.theta,
                      local_diagnostics.low_margin,
                      local_diagnostics.high_margin};
    selected_integers = {
        local_diagnostics.low_order_substeps,
        local_diagnostics.low_order_halo_exchanges,
        low_transport_passes,
        static_cast<std::uint64_t>(local_diagnostics.limiting_cell.x),
        static_cast<std::uint64_t>(local_diagnostics.limiting_cell.y),
        static_cast<std::uint64_t>(local_diagnostics.limiting_cell.z),
        static_cast<std::uint64_t>(local_diagnostics.limiting_rank),
        static_cast<std::uint64_t>(local_diagnostics.constraint),
        1U};
  }
  if (MPI_Bcast(selected_reals.data(),
                static_cast<int>(selected_reals.size()), MPI_DOUBLE,
                global_theta_rank.rank, communicator) != MPI_SUCCESS ||
      MPI_Bcast(selected_integers.data(),
                static_cast<int>(selected_integers.size()), MPI_UINT64_T,
                global_theta_rank.rank, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPredictorPlan};
  }
  blocking_collectives += 2U;
  if (!std::isfinite(selected_reals[0U]) ||
      !std::isfinite(selected_reals[1U]) ||
      !std::isfinite(selected_reals[2U]) ||
      selected_reals[0U] != global_theta_rank.theta ||
      !(selected_reals[1U] >= 0.0) || !(selected_reals[2U] < 0.0) ||
      selected_integers[0U] != substeps ||
      selected_integers[1U] != low_halo_exchanges ||
      selected_integers[2U] != low_transport_passes ||
      selected_integers[0U] >
          std::numeric_limits<std::uint32_t>::max() ||
      selected_integers[1U] >
          std::numeric_limits<std::uint32_t>::max() ||
      selected_integers[3U] >
          static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max()) ||
      selected_integers[4U] >
          static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max()) ||
      selected_integers[5U] >
          static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max()) ||
      selected_integers[6U] >= static_cast<std::uint64_t>(size) ||
      selected_integers[7U] < static_cast<std::uint64_t>(
                                  ThermophysicalAdmissibilityConstraint::
                                      density) ||
      selected_integers[7U] > static_cast<std::uint64_t>(
                                  ThermophysicalAdmissibilityConstraint::
                                      enthalpy_upper) ||
      selected_integers[8U] != 1U) {
    const bool metadata_constraint_valid =
        selected_integers[7U] >= static_cast<std::uint64_t>(
                                     ThermophysicalAdmissibilityConstraint::
                                         density) &&
        selected_integers[7U] <= static_cast<std::uint64_t>(
                                     ThermophysicalAdmissibilityConstraint::
                                         enthalpy_upper);
    const auto metadata_constraint = metadata_constraint_valid
                                         ? static_cast<
                                               ThermophysicalAdmissibilityConstraint>(
                                               selected_integers[7U])
                                         : ThermophysicalAdmissibilityConstraint::
                                               none;
    const bool metadata_cell_valid =
        selected_integers[3U] <= static_cast<std::uint64_t>(
                                     std::numeric_limits<std::int32_t>::max()) &&
        selected_integers[4U] <= static_cast<std::uint64_t>(
                                     std::numeric_limits<std::int32_t>::max()) &&
        selected_integers[5U] <= static_cast<std::uint64_t>(
                                     std::numeric_limits<std::int32_t>::max());
    mark_failure(
        failure, ThermophysicalPredictorFailureReason::limiter_metadata,
        failure_field_for_constraint(metadata_constraint), UINT32_MAX,
        metadata_constraint,
        selected_integers[6U] < static_cast<std::uint64_t>(size)
            ? static_cast<std::int32_t>(selected_integers[6U])
            : -1,
        metadata_cell_valid
            ? Int3{static_cast<std::int32_t>(selected_integers[3U]),
                   static_cast<std::int32_t>(selected_integers[4U]),
                   static_cast<std::int32_t>(selected_integers[5U])}
            : Int3{-1, -1, -1},
        0U, metadata_cell_valid, false,
        std::isfinite(selected_reals[1U]) &&
            std::isfinite(selected_reals[2U]),
        selected_reals[1U], selected_reals[2U]);
    if (failure.valid) {
      record_scalar(failure,
                    ThermophysicalPredictorFailure::scalar_observed_value,
                    failure.observed_value, selected_reals[0U]);
      record_scalar(failure,
                    ThermophysicalPredictorFailure::scalar_allowed_lower,
                    failure.allowed_lower, 0.0);
      record_scalar(failure,
                    ThermophysicalPredictorFailure::scalar_allowed_upper,
                    failure.allowed_upper, 1.0);
    }
    return publish_failure(
        {StatusCode::numerical_failure, kPredictorNumerical});
  }
  ThermophysicalPredictorDiagnostics selected;
  selected.theta = selected_reals[0U];
  selected.low_margin = selected_reals[1U];
  selected.high_margin = selected_reals[2U];
  selected.low_order_substeps =
      static_cast<std::uint32_t>(selected_integers[0U]);
  selected.low_order_halo_exchanges =
      static_cast<std::uint32_t>(selected_integers[1U]);
  selected.low_order_transport_passes = selected_integers[2U];
  selected.mass_flux_scale = low_mass_flux_scale;
  selected.limiting_cell = {
      static_cast<std::int32_t>(selected_integers[3U]),
      static_cast<std::int32_t>(selected_integers[4U]),
      static_cast<std::int32_t>(selected_integers[5U])};
  selected.limiting_rank =
      static_cast<std::int32_t>(selected_integers[6U]);
  selected.constraint = static_cast<ThermophysicalAdmissibilityConstraint>(
      selected_integers[7U]);
  selected.low_state = low_state;
  selected.enthalpy_endpoint = enthalpy_endpoint_result;
  selected.enthalpy_endpoint_alpha = enthalpy_endpoint_alpha;
  selected.bdf_endpoint_alpha = bdf_endpoint_alpha;
  selected.source_endpoint_alpha = source_endpoint_alpha;
  selected.enthalpy_solve_calls = enthalpy_solve_calls;
  selected.limited = true;

  const double theta = global_theta_rank.theta;
  for (std::int32_t z = 0; z < end.z && local; ++z) {
    for (std::int32_t y = 0; y < end.y && local; ++y) {
      for (std::int32_t x = 0; x < end.x && local; ++x) {
        const Int3 cell{x, y, z};
        const double rho_low =
            output.low_order_density_workspace.unchecked(cell, 0U);
        const double rho_high = output.density_workspace.unchecked(cell, 0U);
        if (theta == 0.0) {
          if (!std::isfinite(rho_low) || !(rho_low > 0.0)) {
            mark_failure(
                failure,
                ThermophysicalPredictorFailureReason::blended_density,
                ThermophysicalPredictorFailureField::density, 0U,
                ThermophysicalAdmissibilityConstraint::density, rank,
                {cell.x + patch_begin_.x, cell.y + patch_begin_.y,
                 cell.z + patch_begin_.z},
                0U, true, false, false);
            if (failure.valid) {
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_density_current,
                  failure.density_current, rho_low);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_density_next,
                  failure.density_next, rho_high);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_observed_value,
                  failure.observed_value, rho_low);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_allowed_lower,
                  failure.allowed_lower, 0.0);
            }
            local = {StatusCode::numerical_failure, kPredictorNumerical};
            break;
          }
          output.density_workspace.unchecked(cell, 0U) = rho_low;
          output.enthalpy.unchecked(cell, 0U) =
              output.low_order_enthalpy_workspace.unchecked(cell, 0U);
          for (std::size_t index = 0U; index < species_count; ++index) {
            output.independent_species.data[index].unchecked(cell, 0U) =
                output.low_order_independent_species.data[index].unchecked(
                    cell, 0U);
          }
          for (std::size_t index = 0U; index < passive_count; ++index) {
            output.passive_scalars.data[index].unchecked(cell, 0U) =
                output.low_order_passive_scalars.data[index].unchecked(cell,
                                                                       0U);
          }
          continue;
        }
        const double rho = rho_low + theta * (rho_high - rho_low);
        if (!std::isfinite(rho) || !(rho > 0.0)) {
          mark_failure(
              failure,
              ThermophysicalPredictorFailureReason::blended_density,
              ThermophysicalPredictorFailureField::density, 0U,
              ThermophysicalAdmissibilityConstraint::density, rank,
              {cell.x + patch_begin_.x, cell.y + patch_begin_.y,
               cell.z + patch_begin_.z},
              0U, true, false, false);
          if (failure.valid) {
            record_scalar(
                failure,
                ThermophysicalPredictorFailure::scalar_density_current,
                failure.density_current, rho_low);
            record_scalar(
                failure,
                ThermophysicalPredictorFailure::scalar_density_next,
                failure.density_next, rho_high);
            record_scalar(
                failure,
                ThermophysicalPredictorFailure::scalar_observed_value,
                failure.observed_value, rho);
            record_scalar(
                failure,
                ThermophysicalPredictorFailure::scalar_allowed_lower,
                failure.allowed_lower, 0.0);
          }
          local = {StatusCode::numerical_failure, kPredictorNumerical};
          break;
        }
        const auto blend = [&](double low, double high) noexcept {
          const double conserved_low = rho_low * low;
          const double conserved_high =
              rho_high == 0.0 ? high : rho_high * high;
          return (conserved_low + theta * (conserved_high - conserved_low)) /
                 rho;
        };
        output.enthalpy.unchecked(cell, 0U) = blend(
            output.low_order_enthalpy_workspace.unchecked(cell, 0U),
            output.enthalpy.unchecked(cell, 0U));
        for (std::size_t index = 0U; index < species_count; ++index) {
          output.independent_species.data[index].unchecked(cell, 0U) = blend(
              output.low_order_independent_species.data[index].unchecked(cell,
                                                                          0U),
              output.independent_species.data[index].unchecked(cell, 0U));
        }
        for (std::size_t index = 0U; index < passive_count; ++index) {
          output.passive_scalars.data[index].unchecked(cell, 0U) = blend(
              output.low_order_passive_scalars.data[index].unchecked(cell,
                                                                      0U),
              output.passive_scalars.data[index].unchecked(cell, 0U));
        }
        output.density_workspace.unchecked(cell, 0U) = rho;
      }
    }
  }
  if (local && !low_flux_is_high) {
    double extrapolate_accepted = 0.0;
    double extrapolate_previous = 0.0;
    if (!bdf_matches_dt(input.dt, input.bdf, extrapolate_accepted,
                        extrapolate_previous)) {
      local = {StatusCode::invalid_plan, kPredictorPlan};
    }
    const ConstFaceFieldView accepted_fields[]{
        input.mass_flux_accepted.x, input.mass_flux_accepted.y,
        input.mass_flux_accepted.z};
    const ConstFaceFieldView previous_fields[]{
        input.mass_flux_previous.x, input.mass_flux_previous.y,
        input.mass_flux_previous.z};
    const FaceFieldView paired_fields[]{output.paired_mass_flux.x,
                                        output.paired_mass_flux.y,
                                        output.paired_mass_flux.z};
    for (std::size_t axis = 0U; axis < 3U && local; ++axis) {
      const Int3 extents = paired_fields[axis].extents;
      for (std::int32_t z = 0; z < extents.z && local; ++z) {
        for (std::int32_t y = 0; y < extents.y && local; ++y) {
          for (std::int32_t x = 0; x < extents.x; ++x) {
            const Int3 face{x, y, z};
            const double accepted = accepted_fields[axis].unchecked(face);
            const double previous =
                second_order ? previous_fields[axis].unchecked(face) : 0.0;
            const double high = extrapolate_accepted * accepted +
                                extrapolate_previous * previous;
            const double low =
                low_flux_is_accepted
                    ? accepted
                    : paired_fields[axis].unchecked(face);
            const double blended = low + theta * (high - low);
            if (!std::isfinite(blended)) {
              local = {StatusCode::numerical_failure, kPredictorNumerical};
              break;
            }
            paired_fields[axis].unchecked(face) = blended;
          }
        }
      }
    }
  }
  if (local && slow_path.immersed_interface != nullptr)
    local = slow_path.immersed_interface->zero_interface_flux(
        output.paired_mass_flux);
  if (local &&
      !tuple_admissible(as_const(output.density_workspace),
                        as_const(output.enthalpy),
                        output.independent_species,
                        output.passive_scalars)) {
    capture_tuple_failure(
        as_const(output.density_workspace), as_const(output.enthalpy),
        output.independent_species, output.passive_scalars,
        ThermophysicalPredictorFailureReason::final_tuple_admissibility, 0U);
    local = {StatusCode::numerical_failure, kPredictorNumerical};
  }
  if (local) local = audit_paired_mass();
  consensus = collective_status(communicator, local, rank, size, failure,
                               failure);
  if (!consensus) return publish_failure(consensus);
  ++blocking_collectives;

  high_certificate.state =
      hash_mix(high_certificate.state, double_bits(theta));
  high_certificate.state = hash_mix(high_certificate.state, substeps);
  high_certificate.state =
      hash_mix(high_certificate.state, static_cast<std::uint64_t>(low_state));
  high_certificate.state =
      hash_mix(high_certificate.state, double_bits(low_mass_flux_scale));
  high_certificate.state =
      hash_mix(high_certificate.state, double_bits(enthalpy_endpoint_alpha));
  high_certificate.state =
      hash_mix(high_certificate.state, double_bits(bdf_endpoint_alpha));
  high_certificate.state =
      hash_mix(high_certificate.state, double_bits(source_endpoint_alpha));
  high_certificate.state =
      high_certificate.state == 0U ? 1U : high_certificate.state;
  diagnostics = selected;
  diagnostics.blocking_collectives = blocking_collectives;
  certificate = high_certificate;
  return {};
}

Status ThermophysicalPredictorPlan::predict_high_local(
    const ThermophysicalPredictorInput& input,
    ThermophysicalPredictorOutput output,
    ThermophysicalPredictorCertificate& certificate,
    ThermophysicalPredictorFailure& failure) const noexcept {
  double extrapolate_accepted = 0.0;
  double extrapolate_previous = 0.0;
  const bool second_order = input.bdf.order == 2U;
  if (kernels_ == nullptr || fingerprint_ == 0U || input.time == 0U ||
      input.geometry != geometry_revision_ ||
      input.boundary != boundary_revision_ ||
      input.transport != transport_fingerprint_ ||
      !bdf_matches_dt(input.dt, input.bdf, extrapolate_accepted,
                      extrapolate_previous) ||
      input.species_accepted.size != species_.size() ||
      input.species_previous.size != species_.size() ||
      output.independent_species.size != species_.size() ||
      input.species_nonadvective_rhs.size != species_.size() ||
      input.passive_scalars_accepted.size != passive_scalars_.size() ||
      input.passive_scalars_previous.size != passive_scalars_.size() ||
      output.passive_scalars.size != passive_scalars_.size() ||
      input.passive_scalar_nonadvective_rhs.size !=
          passive_scalars_.size() ||
      (species_.size() != 0U &&
       (input.species_accepted.data == nullptr ||
        input.species_previous.data == nullptr ||
        output.independent_species.data == nullptr ||
        input.species_nonadvective_rhs.data == nullptr)) ||
      (passive_scalars_.size() != 0U &&
       (input.passive_scalars_accepted.data == nullptr ||
        input.passive_scalars_previous.data == nullptr ||
        output.passive_scalars.data == nullptr ||
        input.passive_scalar_nonadvective_rhs.data == nullptr)) ||
      input.density_accepted.field != density_ ||
      (second_order && input.density_previous.field != density_) ||
      input.enthalpy_accepted.field != enthalpy_ ||
      (second_order && input.enthalpy_previous.field != enthalpy_) ||
      output.enthalpy.field != enthalpy_ ||
      !detail::valid_cell_view(input.density_accepted, cells_, 0U, 1U, 0U) ||
      (second_order
           ? !detail::valid_cell_view(input.density_previous, cells_, 0U, 1U,
                                      0U)
           : !empty_field(input.density_previous)) ||
      !detail::valid_cell_view(input.enthalpy_accepted, cells_, 0U, 1U,
                               enthalpy_reach_) ||
      (second_order
           ? !detail::valid_cell_view(input.enthalpy_previous, cells_, 0U, 1U,
                                      enthalpy_reach_)
           : !empty_field(input.enthalpy_previous)) ||
      !detail::valid_cell_view(output.enthalpy, cells_, 0U, 1U) ||
      !detail::valid_cell_view(output.density_workspace, cells_, 0U, 1U) ||
      !detail::valid_cell_view(output.accepted_advection_workspace, cells_,
                               0U, 1U) ||
      !detail::valid_cell_view(output.previous_advection_workspace, cells_,
                               0U, 1U) ||
      !detail::valid_flux_view(output.paired_mass_flux, cells_) ||
      output.paired_mass_flux.revision == 0U ||
      output.paired_mass_flux.certificate.valid() ||
      !valid_flux(input.mass_flux_accepted, cells_) ||
      (second_order &&
       (!valid_flux(input.mass_flux_previous, cells_) ||
        input.mass_flux_previous.certificate.authority() !=
            input.mass_flux_accepted.certificate.authority() ||
        input.mass_flux_previous.certificate.storage() !=
            input.mass_flux_accepted.certificate.storage() ||
        input.mass_flux_previous.certificate.revision_domain() !=
            input.mass_flux_accepted.certificate.revision_domain() ||
        input.mass_flux_previous.revision ==
            input.mass_flux_accepted.revision)) ||
      !valid_rate(input.enthalpy_nonadvective_rhs.accepted, cells_, false) ||
      (second_order
           ? !valid_rate(input.enthalpy_nonadvective_rhs.previous, cells_,
                         false)
           : !empty_field(input.enthalpy_nonadvective_rhs.previous))) {
    mark_failure(
        failure,
        ThermophysicalPredictorFailureReason::high_input_preflight,
        ThermophysicalPredictorFailureField::none, UINT32_MAX,
        ThermophysicalAdmissibilityConstraint::none, -1, {-1, -1, -1}, 0U,
        false, false, false);
    return {StatusCode::invalid_plan, kPredictorPlan};
  }

  for (std::size_t i = 0U; i < species_.size(); ++i) {
    const PredictorRateHistory rate =
        input.species_nonadvective_rhs.data[i];
    if (input.species_accepted.data[i].field != species_[i] ||
        (second_order &&
         input.species_previous.data[i].field != species_[i]) ||
        output.independent_species.data[i].field != species_[i] ||
        !detail::valid_cell_view(input.species_accepted.data[i], cells_, 0U,
                                 1U, species_reach_) ||
        (second_order
             ? !detail::valid_cell_view(input.species_previous.data[i],
                                        cells_, 0U, 1U, species_reach_)
             : !empty_field(input.species_previous.data[i])) ||
        !detail::valid_cell_view(output.independent_species.data[i], cells_,
                                 0U, 1U) ||
        !valid_rate(rate.accepted, cells_, false) ||
        (second_order ? !valid_rate(rate.previous, cells_, false)
                      : !empty_field(rate.previous))) {
      mark_failure(
          failure,
          ThermophysicalPredictorFailureReason::high_input_preflight,
          ThermophysicalPredictorFailureField::independent_species,
          static_cast<std::uint32_t>(i),
          ThermophysicalAdmissibilityConstraint::none, -1, {-1, -1, -1}, 0U,
          false, false, false);
      return {StatusCode::invalid_plan, kPredictorPlan};
    }
  }
  for (std::size_t i = 0U; i < passive_scalars_.size(); ++i) {
    const PredictorRateHistory rate =
        input.passive_scalar_nonadvective_rhs.data[i];
    if (input.passive_scalars_accepted.data[i].field !=
            passive_scalars_[i] ||
        (second_order &&
         input.passive_scalars_previous.data[i].field !=
             passive_scalars_[i]) ||
        output.passive_scalars.data[i].field != passive_scalars_[i] ||
        !detail::valid_cell_view(input.passive_scalars_accepted.data[i],
                                 cells_, 0U, 1U,
                                 passive_scalar_reach_) ||
        (second_order
             ? !detail::valid_cell_view(
                   input.passive_scalars_previous.data[i], cells_, 0U, 1U,
                   passive_scalar_reach_)
             : !empty_field(input.passive_scalars_previous.data[i])) ||
        !detail::valid_cell_view(output.passive_scalars.data[i], cells_, 0U,
                                 1U) ||
        !valid_rate(rate.accepted, cells_, false) ||
        (second_order ? !valid_rate(rate.previous, cells_, false)
                      : !empty_field(rate.previous))) {
      mark_failure(
          failure,
          ThermophysicalPredictorFailureReason::high_input_preflight,
          ThermophysicalPredictorFailureField::passive_scalar,
          static_cast<std::uint32_t>(i),
          ThermophysicalAdmissibilityConstraint::none, -1, {-1, -1, -1}, 0U,
          false, false, false);
      return {StatusCode::invalid_plan, kPredictorPlan};
    }
  }

  const std::array<FieldView, 4U> primary_outputs{
      output.enthalpy, output.density_workspace,
      output.accepted_advection_workspace,
      output.previous_advection_workspace};
  for (std::size_t left = 0U; left < primary_outputs.size(); ++left) {
    for (std::size_t right = left + 1U; right < primary_outputs.size();
         ++right) {
      if (detail::field_views_overlap(as_const(primary_outputs[left]),
                                      as_const(primary_outputs[right]))) {
        return {StatusCode::invalid_plan, kPredictorPlan};
      }
    }
  }
  for (std::size_t i = 0U; i < species_.size(); ++i) {
    const FieldView species_output = output.independent_species.data[i];
    for (FieldView primary : primary_outputs) {
      if (detail::field_views_overlap(as_const(species_output),
                                      as_const(primary))) {
        return {StatusCode::invalid_plan, kPredictorPlan};
      }
    }
    for (std::size_t j = 0U; j < i; ++j) {
      if (detail::field_views_overlap(
              as_const(species_output),
              as_const(output.independent_species.data[j]))) {
        return {StatusCode::invalid_plan, kPredictorPlan};
      }
    }
  }
  for (std::size_t i = 0U; i < passive_scalars_.size(); ++i) {
    const FieldView scalar_output = output.passive_scalars.data[i];
    for (FieldView primary : primary_outputs) {
      if (detail::field_views_overlap(as_const(scalar_output),
                                      as_const(primary))) {
        return {StatusCode::invalid_plan, kPredictorPlan};
      }
    }
    for (std::size_t j = 0U; j < species_.size(); ++j) {
      if (detail::field_views_overlap(
              as_const(scalar_output),
              as_const(output.independent_species.data[j]))) {
        return {StatusCode::invalid_plan, kPredictorPlan};
      }
    }
    for (std::size_t j = 0U; j < i; ++j) {
      if (detail::field_views_overlap(
              as_const(scalar_output),
              as_const(output.passive_scalars.data[j]))) {
        return {StatusCode::invalid_plan, kPredictorPlan};
      }
    }
  }

  const auto aliases_any_output = [&](ConstFieldView view) noexcept {
    if (view.base == nullptr) {
      return false;
    }
    for (FieldView primary : primary_outputs) {
      if (output_aliases_input(primary, view)) {
        return true;
      }
    }
    for (std::size_t i = 0U; i < species_.size(); ++i) {
      if (output_aliases_input(output.independent_species.data[i], view)) {
        return true;
      }
    }
    for (std::size_t i = 0U; i < passive_scalars_.size(); ++i) {
      if (output_aliases_input(output.passive_scalars.data[i], view)) {
        return true;
      }
    }
    return false;
  };
  if (aliases_any_output(input.density_accepted) ||
      aliases_any_output(input.enthalpy_accepted) ||
      aliases_any_output(input.enthalpy_nonadvective_rhs.accepted) ||
      (second_order &&
       (aliases_any_output(input.density_previous) ||
        aliases_any_output(input.enthalpy_previous) ||
        aliases_any_output(input.enthalpy_nonadvective_rhs.previous)))) {
    return {StatusCode::invalid_plan, kPredictorPlan};
  }
  for (std::size_t i = 0U; i < species_.size(); ++i) {
    if (aliases_any_output(input.species_accepted.data[i]) ||
        aliases_any_output(
            input.species_nonadvective_rhs.data[i].accepted) ||
        (second_order &&
         (aliases_any_output(input.species_previous.data[i]) ||
          aliases_any_output(
              input.species_nonadvective_rhs.data[i].previous)))) {
      return {StatusCode::invalid_plan, kPredictorPlan};
    }
  }
  for (std::size_t i = 0U; i < passive_scalars_.size(); ++i) {
    if (aliases_any_output(input.passive_scalars_accepted.data[i]) ||
        aliases_any_output(
            input.passive_scalar_nonadvective_rhs.data[i].accepted) ||
        (second_order &&
         (aliases_any_output(input.passive_scalars_previous.data[i]) ||
          aliases_any_output(
              input.passive_scalar_nonadvective_rhs.data[i].previous)))) {
      return {StatusCode::invalid_plan, kPredictorPlan};
    }
  }
  for (FieldView primary : primary_outputs) {
    if (output_aliases_flux(primary, input.mass_flux_accepted) ||
        (second_order &&
         output_aliases_flux(primary, input.mass_flux_previous))) {
      return {StatusCode::invalid_plan, kPredictorPlan};
    }
    const ConstFaceFluxView paired = as_const(output.paired_mass_flux);
    if (detail::cell_face_views_overlap(primary, paired.x) ||
        detail::cell_face_views_overlap(primary, paired.y) ||
        detail::cell_face_views_overlap(primary, paired.z)) {
      return {StatusCode::invalid_plan, kPredictorPlan};
    }
  }
  for (std::size_t i = 0U; i < species_.size(); ++i) {
    if (output_aliases_flux(output.independent_species.data[i],
                            input.mass_flux_accepted) ||
        (second_order &&
         output_aliases_flux(output.independent_species.data[i],
                             input.mass_flux_previous))) {
      return {StatusCode::invalid_plan, kPredictorPlan};
    }
  }
  for (std::size_t i = 0U; i < passive_scalars_.size(); ++i) {
    if (output_aliases_flux(output.passive_scalars.data[i],
                            input.mass_flux_accepted) ||
        (second_order &&
         output_aliases_flux(output.passive_scalars.data[i],
                             input.mass_flux_previous))) {
      return {StatusCode::invalid_plan, kPredictorPlan};
    }
  }
  const ConstFaceFluxView paired_read = as_const(output.paired_mass_flux);
  const ConstFaceFieldView paired_faces[]{paired_read.x, paired_read.y,
                                           paired_read.z};
  const ConstFaceFieldView accepted_faces[]{input.mass_flux_accepted.x,
                                             input.mass_flux_accepted.y,
                                             input.mass_flux_accepted.z};
  const ConstFaceFieldView previous_faces[]{input.mass_flux_previous.x,
                                             input.mass_flux_previous.y,
                                             input.mass_flux_previous.z};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    if (detail::face_views_overlap(paired_faces[axis],
                                   accepted_faces[axis]) ||
        (second_order &&
         detail::face_views_overlap(paired_faces[axis],
                                    previous_faces[axis]))) {
      return {StatusCode::invalid_plan, kPredictorPlan};
    }
  }

  const KernelBox box{{0, 0, 0}, cells_};
  if (!detail::finite_field_box(input.density_accepted, box, 0U, 1U) ||
      (second_order &&
       !detail::finite_field_box(input.density_previous, box, 0U, 1U)) ||
      !detail::finite_face_flux(input.mass_flux_accepted, box) ||
      (second_order &&
       !detail::finite_face_flux(input.mass_flux_previous, box)) ||
      !detail::finite_face_neighbour_slabs(
          input.enthalpy_accepted, box, 0U, 1U, enthalpy_reach_) ||
      (second_order &&
       !detail::finite_face_neighbour_slabs(
           input.enthalpy_previous, box, 0U, 1U, enthalpy_reach_)) ||
      (!rate_is_zero(input.enthalpy_nonadvective_rhs.accepted) &&
       !detail::finite_field_box(
           input.enthalpy_nonadvective_rhs.accepted, box, 0U, 1U)) ||
       (second_order &&
       !rate_is_zero(input.enthalpy_nonadvective_rhs.previous) &&
       !detail::finite_field_box(
           input.enthalpy_nonadvective_rhs.previous, box, 0U, 1U))) {
    ThermophysicalPredictorFailureField field_kind =
        ThermophysicalPredictorFailureField::none;
    std::uint32_t field_index = UINT32_MAX;
    if (!detail::finite_field_box(input.density_accepted, box, 0U, 1U) ||
        (second_order &&
         !detail::finite_field_box(input.density_previous, box, 0U, 1U))) {
      field_kind = ThermophysicalPredictorFailureField::density;
      field_index = 0U;
    } else if (!detail::finite_face_flux(input.mass_flux_accepted, box) ||
               (second_order &&
                !detail::finite_face_flux(input.mass_flux_previous, box))) {
      field_kind = ThermophysicalPredictorFailureField::mass_flux;
      field_index = UINT32_MAX;
    } else if (!detail::finite_face_neighbour_slabs(
                   input.enthalpy_accepted, box, 0U, 1U,
                   enthalpy_reach_) ||
               (second_order &&
                !detail::finite_face_neighbour_slabs(
                    input.enthalpy_previous, box, 0U, 1U,
                    enthalpy_reach_))) {
      field_kind = ThermophysicalPredictorFailureField::enthalpy;
      field_index = 0U;
    } else {
      field_kind = ThermophysicalPredictorFailureField::nonadvective_rhs;
      field_index = 0U;
    }
    mark_failure(
        failure,
        ThermophysicalPredictorFailureReason::high_input_preflight,
        field_kind, field_index,
        ThermophysicalAdmissibilityConstraint::none, -1, {-1, -1, -1}, 0U,
        false, false, false);
    return {StatusCode::numerical_failure, kPredictorNumerical};
  }
  for (std::size_t i = 0U; i < species_.size(); ++i) {
    const PredictorRateHistory rate =
        input.species_nonadvective_rhs.data[i];
    if (!detail::finite_face_neighbour_slabs(
            input.species_accepted.data[i], box, 0U, 1U, species_reach_) ||
        (second_order &&
         !detail::finite_face_neighbour_slabs(
             input.species_previous.data[i], box, 0U, 1U,
             species_reach_)) ||
        (!rate_is_zero(rate.accepted) &&
         !detail::finite_field_box(rate.accepted, box, 0U, 1U)) ||
        (second_order && !rate_is_zero(rate.previous) &&
         !detail::finite_field_box(rate.previous, box, 0U, 1U))) {
      const bool species_state_bad =
          !detail::finite_face_neighbour_slabs(
              input.species_accepted.data[i], box, 0U, 1U,
              species_reach_) ||
          (second_order &&
           !detail::finite_face_neighbour_slabs(
               input.species_previous.data[i], box, 0U, 1U,
               species_reach_));
      mark_failure(
          failure,
          ThermophysicalPredictorFailureReason::high_input_preflight,
          species_state_bad
              ? ThermophysicalPredictorFailureField::independent_species
              : ThermophysicalPredictorFailureField::nonadvective_rhs,
          static_cast<std::uint32_t>(i),
          ThermophysicalAdmissibilityConstraint::none, -1, {-1, -1, -1}, 0U,
          false, false, false);
      return {StatusCode::numerical_failure, kPredictorNumerical};
    }
  }
  for (std::size_t i = 0U; i < passive_scalars_.size(); ++i) {
    const PredictorRateHistory rate =
        input.passive_scalar_nonadvective_rhs.data[i];
    if (!detail::finite_face_neighbour_slabs(
            input.passive_scalars_accepted.data[i], box, 0U, 1U,
            passive_scalar_reach_) ||
        (second_order &&
         !detail::finite_face_neighbour_slabs(
             input.passive_scalars_previous.data[i], box, 0U, 1U,
             passive_scalar_reach_)) ||
        (!rate_is_zero(rate.accepted) &&
         !detail::finite_field_box(rate.accepted, box, 0U, 1U)) ||
        (second_order && !rate_is_zero(rate.previous) &&
         !detail::finite_field_box(rate.previous, box, 0U, 1U))) {
      const bool passive_state_bad =
          !detail::finite_face_neighbour_slabs(
              input.passive_scalars_accepted.data[i], box, 0U, 1U,
              passive_scalar_reach_) ||
          (second_order &&
           !detail::finite_face_neighbour_slabs(
               input.passive_scalars_previous.data[i], box, 0U, 1U,
               passive_scalar_reach_));
      mark_failure(
          failure,
          ThermophysicalPredictorFailureReason::high_input_preflight,
          passive_state_bad
              ? ThermophysicalPredictorFailureField::passive_scalar
              : ThermophysicalPredictorFailureField::nonadvective_rhs,
          static_cast<std::uint32_t>(i),
          ThermophysicalAdmissibilityConstraint::none, -1, {-1, -1, -1}, 0U,
          false, false, false);
      return {StatusCode::numerical_failure, kPredictorNumerical};
    }
  }

  const Int3 end{cells_.x, cells_.y, cells_.z};
  const ConstFaceFieldView accepted_face_fields[]{
      input.mass_flux_accepted.x, input.mass_flux_accepted.y,
      input.mass_flux_accepted.z};
  const ConstFaceFieldView previous_face_fields[]{
      input.mass_flux_previous.x, input.mass_flux_previous.y,
      input.mass_flux_previous.z};
  const FaceFieldView paired_face_fields[]{output.paired_mass_flux.x,
                                           output.paired_mass_flux.y,
                                           output.paired_mass_flux.z};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    const Int3 extents = paired_face_fields[axis].extents;
    for (std::int32_t z = 0; z < extents.z; ++z) {
      for (std::int32_t y = 0; y < extents.y; ++y) {
        for (std::int32_t x = 0; x < extents.x; ++x) {
          const Int3 face{x, y, z};
          const double value =
              extrapolate_accepted * accepted_face_fields[axis].unchecked(face) +
              extrapolate_previous *
                  (second_order ? previous_face_fields[axis].unchecked(face)
                                : 0.0);
          if (!std::isfinite(value)) {
            mark_failure(
                failure,
                ThermophysicalPredictorFailureReason::high_density_bdf,
                ThermophysicalPredictorFailureField::mass_flux, UINT32_MAX,
                ThermophysicalAdmissibilityConstraint::none, -1,
                {-1, -1, -1}, 0U, false, false, false);
            return {StatusCode::numerical_failure, kPredictorNumerical};
          }
          paired_face_fields[axis].unchecked(face) = value;
        }
      }
    }
  }
  for (std::int32_t z = 0; z < end.z; ++z) {
    for (std::int32_t y = 0; y < end.y; ++y) {
      for (std::int32_t x = 0; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho_n = input.density_accepted.unchecked(cell, 0U);
        const double rho_nm1 =
            second_order ? input.density_previous.unchecked(cell, 0U) : 0.0;
        const double mass_rate_n = flux_divergence_density(
            *kernels_, input.mass_flux_accepted, cell);
        const double mass_rate_nm1 =
            second_order
                ? flux_divergence_density(*kernels_,
                                          input.mass_flux_previous, cell)
                : 0.0;
        const double rho_star =
            (-input.bdf.a1 * rho_n - input.bdf.a2 * rho_nm1 -
             extrapolate_accepted * mass_rate_n -
             extrapolate_previous * mass_rate_nm1) /
            input.bdf.a0;
        if (!std::isfinite(rho_n) || rho_n <= 0.0 ||
            (second_order &&
             (!std::isfinite(rho_nm1) || rho_nm1 <= 0.0)) ||
            !std::isfinite(rho_star)) {
          mark_failure(
              failure,
              ThermophysicalPredictorFailureReason::high_density_bdf,
              ThermophysicalPredictorFailureField::density, 0U,
              ThermophysicalAdmissibilityConstraint::density, -1,
              {cell.x + patch_begin_.x, cell.y + patch_begin_.y,
               cell.z + patch_begin_.z},
              0U, true, false, false);
          if (failure.valid) {
            record_scalar(
                failure,
                ThermophysicalPredictorFailure::scalar_density_current,
                failure.density_current, rho_n);
            record_scalar(
                failure,
                ThermophysicalPredictorFailure::scalar_density_next,
                failure.density_next, rho_star);
            record_scalar(
                failure,
                ThermophysicalPredictorFailure::scalar_divergence,
                failure.divergence, mass_rate_n);
            record_scalar(
                failure,
                ThermophysicalPredictorFailure::scalar_observed_value,
                failure.observed_value, rho_star);
            record_scalar(
                failure,
                ThermophysicalPredictorFailure::scalar_allowed_lower,
                failure.allowed_lower, 0.0);
          }
          return {StatusCode::numerical_failure, kPredictorNumerical};
        }
        output.density_workspace.unchecked(cell, 0U) = rho_star;
      }
    }
  }

  const auto predict_quantity = [&](ConstFieldView accepted,
                                    ConstFieldView previous,
                                    PredictorRateHistory nonadvective_rhs,
                                    ConvectionScheme convection,
                                    FieldView predicted,
                                    ThermophysicalPredictorFailureField
                                        field_kind,
                                    std::uint32_t field_index) noexcept -> Status {
    const std::array<ConstFieldView, 1U> accepted_reads{accepted};
    const std::array<FieldView, 1U> accepted_writes{
        output.accepted_advection_workspace};
    const KernelInvocation accepted_call{
        {accepted_reads.data(), accepted_reads.size()},
        {accepted_writes.data(), accepted_writes.size()}, box, 0U, 0U, 1U,
        input.mass_flux_accepted.revision, input.counters};
    Status status = cartesian_convection(*kernels_, convection,
                                         input.mass_flux_accepted,
                                         accepted_call);
    if (!status) {
      if (status.code == StatusCode::numerical_failure) {
        mark_failure(
            failure,
            ThermophysicalPredictorFailureReason::high_quantity_transport,
            field_kind, field_index,
            ThermophysicalAdmissibilityConstraint::none, -1, {-1, -1, -1},
            0U, false, false, false);
      }
      return status;
    }
    if (second_order) {
      const std::array<ConstFieldView, 1U> previous_reads{previous};
      const std::array<FieldView, 1U> previous_writes{
          output.previous_advection_workspace};
      const KernelInvocation previous_call{
          {previous_reads.data(), previous_reads.size()},
          {previous_writes.data(), previous_writes.size()}, box, 0U, 0U, 1U,
          input.mass_flux_previous.revision, input.counters};
      status = cartesian_convection(*kernels_, convection,
                                    input.mass_flux_previous,
                                    previous_call);
      if (!status) {
        if (status.code == StatusCode::numerical_failure) {
          mark_failure(
              failure,
              ThermophysicalPredictorFailureReason::high_quantity_transport,
              field_kind, field_index,
              ThermophysicalAdmissibilityConstraint::none, -1,
              {-1, -1, -1}, 0U, false, false, false);
        }
        return status;
      }
    }

    for (std::int32_t z = 0; z < end.z; ++z) {
      for (std::int32_t y = 0; y < end.y; ++y) {
        for (std::int32_t x = 0; x < end.x; ++x) {
          const Int3 cell{x, y, z};
          const double rho_n = input.density_accepted.unchecked(cell, 0U);
          const double rho_nm1 =
              second_order ? input.density_previous.unchecked(cell, 0U) : 0.0;
          const double rho_star =
              output.density_workspace.unchecked(cell, 0U);
          const double rhs_n = rate_is_zero(nonadvective_rhs.accepted)
                                   ? 0.0
                                   : nonadvective_rhs.accepted.unchecked(cell,
                                                                          0U);
          const double rhs_nm1 =
              !second_order || rate_is_zero(nonadvective_rhs.previous)
                  ? 0.0
                  : nonadvective_rhs.previous.unchecked(cell, 0U);
          const double transport_n =
              output.accepted_advection_workspace.unchecked(cell, 0U) - rhs_n;
          const double transport_nm1 =
              second_order
                  ? output.previous_advection_workspace.unchecked(cell, 0U) -
                        rhs_nm1
                  : 0.0;
          const double previous_quantity =
              second_order ? previous.unchecked(cell, 0U) : 0.0;
          const double rho_q_star =
              (-input.bdf.a1 * rho_n * accepted.unchecked(cell, 0U) -
               input.bdf.a2 * rho_nm1 * previous_quantity -
               extrapolate_accepted * transport_n -
               extrapolate_previous * transport_nm1) /
              input.bdf.a0;
          // An exactly zero high density cannot represent an intensive
          // quantity.  Store its conserved numerator as a slow-path-only
          // marker; the wrapper reconstructs the common conservative blend.
          const double value =
              rho_star == 0.0 ? rho_q_star : rho_q_star / rho_star;
          if (!std::isfinite(rho_star) || !std::isfinite(rho_q_star) ||
              !std::isfinite(value)) {
            mark_failure(
                failure,
                ThermophysicalPredictorFailureReason::high_quantity_bdf,
                field_kind, field_index,
                ThermophysicalAdmissibilityConstraint::none, -1,
                {cell.x + patch_begin_.x, cell.y + patch_begin_.y,
                 cell.z + patch_begin_.z},
                0U, true, false, false);
            if (failure.valid) {
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_density_current,
                  failure.density_current, rho_n);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_density_next,
                  failure.density_next, rho_star);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_quantity_current,
                  failure.quantity_current,
                  accepted.unchecked(cell, 0U));
              if (second_order) {
                record_scalar(
                    failure,
                    ThermophysicalPredictorFailure::scalar_quantity_previous,
                    failure.quantity_previous, previous_quantity);
              }
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::
                      scalar_nonadvective_rhs,
                  failure.nonadvective_rhs, rhs_n);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_conserved_next,
                  failure.conserved_next, rho_q_star);
              record_scalar(
                  failure,
                  ThermophysicalPredictorFailure::scalar_observed_value,
                  failure.observed_value, value);
            }
            return {StatusCode::numerical_failure, kPredictorNumerical};
          }
          predicted.unchecked(cell, 0U) = value;
        }
      }
    }
    return {};
  };

  Status status = predict_quantity(
      input.enthalpy_accepted, input.enthalpy_previous,
      input.enthalpy_nonadvective_rhs, enthalpy_convection_, output.enthalpy,
      ThermophysicalPredictorFailureField::enthalpy, 0U);
  if (!status) {
    return status;
  }
  for (std::size_t i = 0U; i < species_.size(); ++i) {
    status = predict_quantity(input.species_accepted.data[i],
                              input.species_previous.data[i],
                              input.species_nonadvective_rhs.data[i],
                              species_convection_,
                              output.independent_species.data[i],
                              ThermophysicalPredictorFailureField::
                                  independent_species,
                              static_cast<std::uint32_t>(i));
    if (!status) {
      return status;
    }
  }
  for (std::size_t i = 0U; i < passive_scalars_.size(); ++i) {
    status = predict_quantity(
        input.passive_scalars_accepted.data[i],
        input.passive_scalars_previous.data[i],
        input.passive_scalar_nonadvective_rhs.data[i],
        passive_scalar_convection_, output.passive_scalars.data[i],
        ThermophysicalPredictorFailureField::passive_scalar,
        static_cast<std::uint32_t>(i));
    if (!status) {
      return status;
    }
  }

  ThermophysicalPredictorCertificate candidate;
  candidate.plan = fingerprint_;
  candidate.time = input.time;
  candidate.geometry = input.geometry;
  candidate.accepted_face_flux = input.mass_flux_accepted.revision;
  candidate.previous_face_flux =
      second_order ? input.mass_flux_previous.revision : 0U;
  candidate.committed_face_flux_authority =
      input.mass_flux_accepted.certificate.authority();
  candidate.committed_face_flux_storage =
      input.mass_flux_accepted.certificate.storage();
  candidate.committed_face_flux_revision_domain =
      input.mass_flux_accepted.certificate.revision_domain();
  candidate.predicted_density = output.density_workspace.revision;
  candidate.predicted_density_storage =
      output.density_workspace.storage_identity;
  candidate.predicted_density_revision_domain =
      output.density_workspace.revision_domain;
  candidate.paired_face_flux = output.paired_mass_flux.revision;
  candidate.paired_face_flux_storage = output.paired_mass_flux.x.storage_identity;
  candidate.paired_face_flux_revision_domain =
      output.paired_mass_flux.x.revision_domain;
  candidate.state = predictor_state_hash(input, output);
  candidate.order = input.bdf.order;
  certificate = candidate;
  return {};
}

namespace {

constexpr std::uint64_t kClosedGaugePrepareSchema =
    UINT64_C(0x7630346367707265);

std::uint64_t nonzero_hash(std::uint64_t hash) noexcept {
  return hash == 0U ? UINT64_C(0x9e3779b97f4a7c15) : hash;
}

bool same_pressure_reference_certificate(
    const PressureReferenceCertificate& left,
    const PressureReferenceCertificate& right) noexcept {
  return left.plan == right.plan && left.predictor == right.predictor &&
         left.thermodynamics == right.thermodynamics &&
         left.closure == right.closure && left.time == right.time &&
         left.pressure_reference == right.pressure_reference &&
         left.kind == right.kind;
}

std::size_t gauge_cell_count(Int3 cells) noexcept {
  return static_cast<std::size_t>(cells.x) *
         static_cast<std::size_t>(cells.y) *
         static_cast<std::size_t>(cells.z);
}

std::size_t gauge_cell_offset(Int3 cells, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(cells.y) *
                  static_cast<std::size_t>(cell.z));
}

bool valid_gauge_activity(PressureEnergyCellActivity activity,
                          Int3 cells) noexcept {
  if (activity.cells.size == 0U) {
    return activity.cells.data == nullptr &&
           activity.local_fingerprint == 0U &&
           activity.collective_fingerprint == 0U;
  }
  if (activity.cells.data == nullptr ||
      activity.cells.size != gauge_cell_count(cells) ||
      activity.local_fingerprint == 0U ||
      activity.collective_fingerprint == 0U) {
    return false;
  }
  for (std::size_t cell = 0U; cell < activity.cells.size; ++cell) {
    if (activity.cells.data[cell] > 1U) return false;
  }
  return true;
}

bool gauge_cell_is_active(PressureEnergyCellActivity activity, Int3 cells,
                          Int3 cell) noexcept {
  return activity.cells.size == 0U ||
         activity.cells.data[gauge_cell_offset(cells, cell)] == 1U;
}

std::uint64_t mix_gauge_collective_field(std::uint64_t hash,
                                         ConstFieldView field) noexcept {
  hash = hash_mix(hash, field.field);
  return hash_mix(hash, field.revision != 0U);
}

std::uint64_t mix_gauge_local_field(std::uint64_t hash,
                                    ConstFieldView field) noexcept {
  hash = hash_mix(
      hash, static_cast<std::uint64_t>(
                reinterpret_cast<std::uintptr_t>(field.base)));
  hash = hash_mix(hash, static_cast<std::uint32_t>(field.interior.x));
  hash = hash_mix(hash, static_cast<std::uint32_t>(field.interior.y));
  hash = hash_mix(hash, static_cast<std::uint32_t>(field.interior.z));
  hash = hash_mix(hash, static_cast<std::uint32_t>(field.ghosts.x));
  hash = hash_mix(hash, static_cast<std::uint32_t>(field.ghosts.y));
  hash = hash_mix(hash, static_cast<std::uint32_t>(field.ghosts.z));
  hash = hash_mix(hash, field.components);
  hash = hash_mix(hash, field.stride_y);
  hash = hash_mix(hash, field.stride_z);
  hash = hash_mix(hash, field.component_stride);
  hash = hash_mix(hash, field.replica);
  hash = hash_mix(hash, field.field);
  hash = hash_mix(hash, field.revision);
  hash = hash_mix(hash, field.storage_identity);
  return hash_mix(hash, field.revision_domain);
}

class GaugeCompensatedSum {
 public:
  void add(double value) noexcept {
    const double next = sum_ + value;
    correction_ += std::abs(sum_) >= std::abs(value)
                       ? (sum_ - next) + value
                       : (value - next) + sum_;
    sum_ = next;
  }

  double value() const noexcept { return sum_ + correction_; }

 private:
  double sum_{};
  double correction_{};
};

struct ClosedGaugeLocalEvaluation {
  double moment{};
  double weight{};
  std::size_t active_cells{};
  RevisionToken content{};
};

struct ClosedGaugePostShiftEvaluation {
  double moment{};
  double absolute_moment{};
};

Status evaluate_closed_gauge_local(
    const ClosedGaugeCorrectionPrepareInput& input,
    const CartesianKernelPlan& kernels, Int3 cells,
    ClosedGaugeLocalEvaluation& evaluation) noexcept {
  GaugeCompensatedSum moment;
  GaugeCompensatedSum weight;
  std::size_t active_cells = 0U;
  std::uint64_t content = hash_mix(kFnvOffset, kClosedGaugePrepareSchema);
  content = mix_gauge_local_field(content, input.pressure_perturbation);
  content = mix_gauge_local_field(content, input.raw_pressure_correction);
  content = mix_gauge_local_field(
      content, input.candidate_pressure_compressibility);
  content = hash_mix(content, input.activity.local_fingerprint);

  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double pi = input.pressure_perturbation.unchecked(cell, 0U);
        const double correction =
            input.raw_pressure_correction.unchecked(cell, 0U);
        const double compressibility =
            input.candidate_pressure_compressibility.unchecked(cell, 0U);
        const double corrected_pi = pi + correction;
        const double absolute_pressure =
            input.pressure_reference + corrected_pi;
        const bool active =
            gauge_cell_is_active(input.activity, cells, cell);
        content = hash_mix(content, double_bits(pi));
        content = hash_mix(content, double_bits(correction));
        content = hash_mix(content, double_bits(compressibility));
        content = hash_mix(content, active ? 1U : 0U);
        if (!std::isfinite(pi) || !std::isfinite(correction) ||
            !std::isfinite(corrected_pi) ||
            !std::isfinite(absolute_pressure) || absolute_pressure <= 0.0) {
          return {StatusCode::numerical_failure, kPredictorNumerical};
        }
        if (!active) continue;

        const double volume = detail::cell_volume(kernels, cell);
        const double cell_weight = volume * compressibility;
        const double cell_moment = cell_weight * corrected_pi;
        content = hash_mix(content, double_bits(volume));
        if (!std::isfinite(compressibility) || compressibility <= 0.0 ||
            !std::isfinite(volume) || volume <= 0.0 ||
            !std::isfinite(cell_weight) || cell_weight <= 0.0 ||
            !std::isfinite(cell_moment)) {
          return {StatusCode::numerical_failure, kPredictorNumerical};
        }
        weight.add(cell_weight);
        moment.add(cell_moment);
        ++active_cells;
      }
    }
  }
  const double local_moment = moment.value();
  const double local_weight = weight.value();
  if (!std::isfinite(local_moment) || !std::isfinite(local_weight) ||
      local_weight < 0.0 ||
      (active_cells != 0U && local_weight <= 0.0) ||
      (active_cells == 0U &&
       (local_moment != 0.0 || local_weight != 0.0))) {
    return {StatusCode::numerical_failure, kPredictorNumerical};
  }
  evaluation = {local_moment, local_weight, active_cells,
                nonzero_hash(content)};
  return {};
}

Status evaluate_closed_gauge_post_shift_local(
    const ClosedGaugeCorrectionPrepareInput& input,
    const CartesianKernelPlan& kernels, Int3 cells, double shift,
    double next_pressure_reference,
    ClosedGaugePostShiftEvaluation& evaluation) noexcept {
  GaugeCompensatedSum moment;
  GaugeCompensatedSum absolute_moment;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double pi = input.pressure_perturbation.unchecked(cell, 0U);
        const double correction =
            input.raw_pressure_correction.unchecked(cell, 0U);
        // Match the required no-fail write expression exactly.
        const double next_pi = (pi + correction) - shift;
        const double next_absolute_pressure =
            next_pressure_reference + next_pi;
        if (!std::isfinite(next_pi) ||
            !std::isfinite(next_absolute_pressure) ||
            next_absolute_pressure <= 0.0) {
          return {StatusCode::numerical_failure, kPredictorNumerical};
        }
        if (!gauge_cell_is_active(input.activity, cells, cell)) continue;
        const double compressibility =
            input.candidate_pressure_compressibility.unchecked(cell, 0U);
        const double volume = detail::cell_volume(kernels, cell);
        const double weighted_pi = volume * compressibility * next_pi;
        if (!std::isfinite(weighted_pi)) {
          return {StatusCode::numerical_failure, kPredictorNumerical};
        }
        moment.add(weighted_pi);
        absolute_moment.add(std::abs(weighted_pi));
      }
    }
  }
  const double local_moment = moment.value();
  const double local_absolute_moment = absolute_moment.value();
  if (!std::isfinite(local_moment) ||
      !std::isfinite(local_absolute_moment) ||
      local_absolute_moment < 0.0) {
    return {StatusCode::numerical_failure, kPredictorNumerical};
  }
  evaluation = {local_moment, local_absolute_moment};
  return {};
}

PlanFingerprint closed_gauge_collective_contract(
    PlanFingerprint plan_fingerprint,
    const ClosedGaugeCorrectionPrepareInput& input) noexcept {
  std::uint64_t hash = hash_mix(kFnvOffset, kClosedGaugePrepareSchema);
  hash = hash_mix(hash, plan_fingerprint);
  hash = hash_mix(hash, input.predecessor.plan);
  hash = hash_mix(hash, input.predecessor.predictor);
  hash = hash_mix(hash, input.predecessor.thermodynamics);
  hash = hash_mix(hash, input.predecessor.closure);
  // A predecessor pressure-reference revision is a rank-local state token.
  // Its exact value is mixed below into rank_local_transaction only.
  hash = hash_mix(hash, input.predecessor.pressure_reference != 0U);
  hash = hash_mix(hash, static_cast<std::uint8_t>(input.predecessor.kind));
  hash = hash_mix(hash, double_bits(input.pressure_reference));
  hash = hash_mix(hash, input.corrector);
  hash = hash_mix(hash, input.time);
  hash = hash_mix(hash, input.geometry);
  hash = hash_mix(hash, input.pressure_correction_authority != 0U);
  hash = hash_mix(hash, input.target_thermodynamic_closure != 0U);
  hash = mix_gauge_collective_field(hash, input.pressure_perturbation);
  hash = mix_gauge_collective_field(hash, input.raw_pressure_correction);
  hash = mix_gauge_collective_field(
      hash, input.candidate_pressure_compressibility);
  hash = hash_mix(hash, input.activity.collective_fingerprint);
  return nonzero_hash(hash);
}

PlanFingerprint closed_gauge_collective_transaction(
    PlanFingerprint contract, double shift, double next_pressure_reference,
    double global_moment, double global_weight,
    double global_post_shift_moment,
    double global_post_shift_absolute_moment,
    double post_shift_gauge_residual,
    double post_shift_gauge_tolerance) noexcept {
  std::uint64_t hash = hash_mix(contract, double_bits(global_moment));
  hash = hash_mix(hash, double_bits(global_weight));
  hash = hash_mix(hash, double_bits(shift));
  hash = hash_mix(hash, double_bits(next_pressure_reference));
  hash = hash_mix(hash, double_bits(global_post_shift_moment));
  hash = hash_mix(hash, double_bits(global_post_shift_absolute_moment));
  hash = hash_mix(hash, double_bits(post_shift_gauge_residual));
  hash = hash_mix(hash, double_bits(post_shift_gauge_tolerance));
  return nonzero_hash(hash);
}

double closed_gauge_post_shift_tolerance(
    double shift, double global_weight,
    double global_post_shift_absolute_moment) noexcept {
  const double absolute_mean =
      global_post_shift_absolute_moment / global_weight;
  const double scale = std::max({1.0, std::abs(shift), absolute_mean});
  return 512.0 * std::numeric_limits<double>::epsilon() * scale;
}

RevisionToken closed_gauge_rank_local_transaction(
    PlanFingerprint collective_transaction,
    const ClosedGaugeCorrectionPrepareInput& input,
    const ClosedGaugeLocalEvaluation& local,
    const ClosedGaugePostShiftEvaluation& post_shift) noexcept {
  std::uint64_t hash = hash_mix(collective_transaction,
                                input.predecessor.pressure_reference);
  hash = hash_mix(hash, input.predecessor.time);
  hash = hash_mix(hash, input.pressure_correction_authority);
  hash = hash_mix(hash, input.target_thermodynamic_closure);
  hash = hash_mix(hash, input.activity.local_fingerprint);
  hash = hash_mix(hash, local.content);
  hash = hash_mix(hash, double_bits(local.moment));
  hash = hash_mix(hash, double_bits(local.weight));
  hash = hash_mix(hash, local.active_cells);
  hash = hash_mix(hash, double_bits(post_shift.moment));
  hash = hash_mix(hash, double_bits(post_shift.absolute_moment));
  return nonzero_hash(hash);
}

bool closed_gauge_structural_input_valid(
    const ClosedGaugeCorrectionPrepareInput& input, Int3 cells,
    PlanFingerprint plan_fingerprint, PlanFingerprint predictor_fingerprint,
    PlanFingerprint thermodynamics_fingerprint,
    RevisionToken geometry_revision,
    FieldId pressure_perturbation_field,
    FieldId compressibility_field) noexcept {
  return input.predecessor.valid() &&
         input.predecessor.plan == plan_fingerprint &&
         input.predecessor.predictor == predictor_fingerprint &&
         input.predecessor.thermodynamics == thermodynamics_fingerprint &&
         input.predecessor.kind == PressureReferenceKind::closed_mass &&
         input.predecessor.time == input.time &&
         std::isfinite(input.pressure_reference) &&
         input.pressure_reference > 0.0 &&
         (input.corrector == 1U || input.corrector == 2U) &&
         input.time != 0U && input.geometry == geometry_revision &&
         input.pressure_correction_authority != 0U &&
         input.target_thermodynamic_closure != 0U &&
         detail::valid_cell_view(input.pressure_perturbation, cells, 0U, 1U,
                                 0U) &&
         input.pressure_perturbation.field == pressure_perturbation_field &&
         detail::valid_cell_view(input.raw_pressure_correction, cells, 0U, 1U,
                                 0U) &&
         detail::valid_cell_view(input.candidate_pressure_compressibility,
                                 cells, 0U, 1U, 0U) &&
         input.candidate_pressure_compressibility.field ==
             compressibility_field &&
         !detail::field_views_overlap(input.pressure_perturbation,
                                      input.raw_pressure_correction) &&
         !detail::field_views_overlap(
             input.pressure_perturbation,
             input.candidate_pressure_compressibility) &&
         !detail::field_views_overlap(
             input.raw_pressure_correction,
             input.candidate_pressure_compressibility) &&
         valid_gauge_activity(input.activity, cells);
}

}  // namespace

Status PressureReferencePlan::solve_closed_mass(
    MPI_Comm communicator, const ClosedMassPlan& closed_mass,
    const ThermodynamicsPlan& thermodynamics,
    const ThermophysicalPredictorCertificate& predictor, StageId stage,
    const ClosedMassCellView& cells, double target_mass,
    double current_pressure_reference, ClosedMassResult& result,
    PressureReferenceCertificate& certificate) const noexcept {
  if (kind_ != PressureReferenceKind::closed_mass ||
      gauge_ != PressureGauge::compressibility_weighted_zero_mean ||
      service_stage_ == 0U || stage != service_stage_ ||
      fingerprint_ == 0U || predictor_fingerprint_ == 0U ||
      !predictor.valid() || predictor.plan != predictor_fingerprint_ ||
      cells.predictor_state != predictor.state ||
      thermodynamics.fingerprint() != thermodynamics_fingerprint_ ||
      closed_mass.authority() != PressureReferenceKind::closed_mass ||
      closed_mass.fingerprint() == 0U) {
    return {StatusCode::invalid_plan, kPredictorPlan};
  }
  ClosedMassResult candidate_result;
  const Status status = closed_mass.solve(
      communicator, thermodynamics, cells, target_mass,
      current_pressure_reference, candidate_result);
  if (!status) {
    return status;
  }
  std::uint64_t revision = kFnvOffset;
  revision = hash_mix(revision, fingerprint_);
  revision = hash_mix(revision, predictor.state);
  revision = hash_mix(revision, predictor.time);
  revision = hash_mix(revision, closed_mass.fingerprint());
  revision = hash_mix(revision, double_bits(candidate_result.pressure_reference));
  revision = hash_mix(revision, double_bits(candidate_result.mass));
  revision = hash_mix(revision, double_bits(candidate_result.residual));
  revision = revision == 0U ? 1U : revision;
  PressureReferenceCertificate candidate_certificate;
  candidate_certificate.plan = fingerprint_;
  candidate_certificate.predictor = predictor.plan;
  candidate_certificate.thermodynamics = thermodynamics_fingerprint_;
  candidate_certificate.closure = closed_mass.fingerprint();
  candidate_certificate.time = predictor.time;
  candidate_certificate.pressure_reference = revision;
  candidate_certificate.kind = kind_;
  result = candidate_result;
  certificate = candidate_certificate;
  return {};
}

Status PressureReferencePlan::solve_closed_mass_fields(
    MPI_Comm communicator, const ClosedMassPlan& closed_mass,
    const ThermodynamicsPlan& thermodynamics,
    const ThermophysicalPredictorCertificate& predictor, StageId stage,
    const ClosedMassFieldView& cells, double target_mass,
    double current_pressure_reference, ClosedMassResult& result,
    PressureReferenceCertificate& certificate) const noexcept {
  if (kind_ != PressureReferenceKind::closed_mass ||
      gauge_ != PressureGauge::compressibility_weighted_zero_mean ||
      service_stage_ == 0U || stage != service_stage_ ||
      fingerprint_ == 0U || predictor_fingerprint_ == 0U ||
      !predictor.valid() || predictor.plan != predictor_fingerprint_ ||
      cells.predictor_state != predictor.state ||
      thermodynamics.fingerprint() != thermodynamics_fingerprint_ ||
      closed_mass.authority() != PressureReferenceKind::closed_mass ||
      closed_mass.fingerprint() == 0U) {
    return {StatusCode::invalid_plan, kPredictorPlan};
  }
  ClosedMassResult candidate_result;
  const Status status = closed_mass.solve_fields(
      communicator, thermodynamics, cells, target_mass,
      current_pressure_reference, candidate_result);
  if (!status) return status;
  std::uint64_t revision = kFnvOffset;
  revision = hash_mix(revision, fingerprint_);
  revision = hash_mix(revision, predictor.state);
  revision = hash_mix(revision, predictor.time);
  revision = hash_mix(revision, closed_mass.fingerprint());
  revision = hash_mix(revision,
                      double_bits(candidate_result.pressure_reference));
  revision = hash_mix(revision, double_bits(candidate_result.mass));
  revision = hash_mix(revision, double_bits(candidate_result.residual));
  revision = revision == 0U ? 1U : revision;
  const PressureReferenceCertificate candidate_certificate{
      fingerprint_, predictor.plan, thermodynamics_fingerprint_,
      closed_mass.fingerprint(), predictor.time, revision, kind_};
  result = candidate_result;
  certificate = candidate_certificate;
  return {};
}

Status PressureReferencePlan::certify_closed_mass_density_fields(
    MPI_Comm communicator, const ClosedMassPlan& closed_mass,
    const ThermodynamicsPlan& thermodynamics,
    const ThermophysicalPredictorCertificate& predictor, StageId stage,
    const ClosedMassDensityFieldView& cells, double target_mass,
    double current_pressure_reference, ClosedMassResult& result,
    PressureReferenceCertificate& certificate) const noexcept {
  if (kind_ != PressureReferenceKind::closed_mass ||
      gauge_ != PressureGauge::compressibility_weighted_zero_mean ||
      service_stage_ == 0U || stage != service_stage_ ||
      fingerprint_ == 0U || predictor_fingerprint_ == 0U ||
      !predictor.valid() || predictor.plan != predictor_fingerprint_ ||
      cells.predictor_state != predictor.state ||
      thermodynamics.fingerprint() != thermodynamics_fingerprint_ ||
      closed_mass.authority() != PressureReferenceKind::closed_mass ||
      closed_mass.fingerprint() == 0U) {
    return {StatusCode::invalid_plan, kPredictorPlan};
  }
  ClosedMassResult candidate_result;
  const Status status = closed_mass.certify_density_fields(
      communicator, cells, target_mass, current_pressure_reference,
      candidate_result);
  if (!status) return status;
  std::uint64_t revision = kFnvOffset;
  revision = hash_mix(revision, fingerprint_);
  revision = hash_mix(revision, predictor.state);
  revision = hash_mix(revision, predictor.time);
  revision = hash_mix(revision, closed_mass.fingerprint());
  revision = hash_mix(revision, cells.pressure_perturbation.revision);
  revision = hash_mix(revision, cells.density.revision);
  revision = hash_mix(revision, cells.pressure_compressibility.revision);
  revision = hash_mix(revision,
                      double_bits(candidate_result.pressure_reference));
  revision = hash_mix(revision, double_bits(candidate_result.mass));
  revision = hash_mix(revision, double_bits(candidate_result.residual));
  revision = revision == 0U ? 1U : revision;
  const PressureReferenceCertificate candidate_certificate{
      fingerprint_, predictor.plan, thermodynamics_fingerprint_,
      closed_mass.fingerprint(), predictor.time, revision, kind_};
  result = candidate_result;
  certificate = candidate_certificate;
  return {};
}

Status PressureReferencePlan::normalize_closed_gauge(
    FieldView pressure_perturbation, ConstFieldView drho_dp_h_y,
    ReductionEngine& reductions, double& pressure_reference) const noexcept {
  if (kind_ != PressureReferenceKind::closed_mass ||
      gauge_ != PressureGauge::compressibility_weighted_zero_mean ||
      kernels_ == nullptr || fingerprint_ == 0U ||
      !std::isfinite(pressure_reference) || pressure_reference <= 0.0 ||
      !detail::valid_cell_view(pressure_perturbation, cells_, 0U, 1U) ||
      !detail::valid_cell_view(drho_dp_h_y, cells_, 0U, 1U, 0U) ||
      detail::field_views_overlap(as_const(pressure_perturbation),
                                  drho_dp_h_y)) {
    return {StatusCode::invalid_plan, kPredictorPlan};
  }
  double local[2]{};
  for (std::int32_t z = 0; z < cells_.z; ++z) {
    for (std::int32_t y = 0; y < cells_.y; ++y) {
      for (std::int32_t x = 0; x < cells_.x; ++x) {
        const Int3 cell{x, y, z};
        const double pi = pressure_perturbation.unchecked(cell, 0U);
        const double compressibility = drho_dp_h_y.unchecked(cell, 0U);
        const double volume = detail::cell_volume(*kernels_, cell);
        if (!std::isfinite(pi) || !std::isfinite(compressibility) ||
            compressibility <= 0.0 || !std::isfinite(volume) ||
            volume <= 0.0 || !std::isfinite(pressure_reference + pi) ||
            pressure_reference + pi <= 0.0) {
          return {StatusCode::numerical_failure, kPredictorNumerical};
        }
        local[0U] += volume * compressibility * pi;
        local[1U] += volume * compressibility;
      }
    }
  }
  double global[2]{};
  Status status = reductions.checked_sum({local, 2U}, {global, 2U});
  if (!status) return status;
  if (!std::isfinite(global[0U]) || !std::isfinite(global[1U]) ||
      global[1U] <= 0.0) {
    return {StatusCode::numerical_failure, kPredictorNumerical};
  }
  const double shift = global[0U] / global[1U];
  const double next_reference = pressure_reference + shift;
  if (!std::isfinite(shift) || !std::isfinite(next_reference) ||
      next_reference <= 0.0) {
    return {StatusCode::numerical_failure, kPredictorNumerical};
  }
  for (std::int32_t z = 0; z < cells_.z; ++z) {
    for (std::int32_t y = 0; y < cells_.y; ++y) {
      for (std::int32_t x = 0; x < cells_.x; ++x) {
        const Int3 cell{x, y, z};
        pressure_perturbation.unchecked(cell, 0U) -= shift;
      }
    }
  }
  pressure_reference = next_reference;
  return {};
}

Status PressureReferencePlan::prepare_closed_gauge_correction(
    const ClosedGaugeCorrectionPrepareInput& input,
    ReductionEngine& reductions,
    ClosedGaugeCorrectionCertificate& certificate) const noexcept {
  const bool local_structure =
      kind_ == PressureReferenceKind::closed_mass &&
      gauge_ == PressureGauge::compressibility_weighted_zero_mean &&
      kernels_ != nullptr && fingerprint_ != 0U &&
      predictor_fingerprint_ != 0U && thermodynamics_fingerprint_ != 0U &&
      geometry_revision_ != 0U &&
      closed_gauge_structural_input_valid(
          input, cells_, fingerprint_, predictor_fingerprint_,
          thermodynamics_fingerprint_, geometry_revision_,
          pressure_perturbation_field_,
          compressibility_field_);
  const PlanFingerprint contract =
      local_structure
          ? closed_gauge_collective_contract(fingerprint_, input)
          : 0U;
  Status status = reductions.consensus_contract(contract);
  if (!status) return status;

  ClosedGaugeLocalEvaluation local_evaluation;
  Status local = evaluate_closed_gauge_local(input, *kernels_, cells_,
                                             local_evaluation);
  const double local_values[2U]{local_evaluation.moment,
                                local_evaluation.weight};
  double global_values[2U]{};
  status = reductions.checked_sum({local_values, 2U}, {global_values, 2U},
                                  local);
  if (!status) return status;

  const double global_moment = global_values[0U];
  const double global_weight = global_values[1U];
  if (!std::isfinite(global_moment) || !std::isfinite(global_weight) ||
      global_weight <= 0.0) {
    return {StatusCode::numerical_failure, kPredictorNumerical};
  }
  const double shift = global_moment / global_weight;
  const double next_pressure_reference = input.pressure_reference + shift;
  if (!std::isfinite(shift) || !std::isfinite(next_pressure_reference) ||
      next_pressure_reference <= 0.0) {
    return {StatusCode::numerical_failure, kPredictorNumerical};
  }

  ClosedGaugePostShiftEvaluation local_post_shift;
  local = evaluate_closed_gauge_post_shift_local(
      input, *kernels_, cells_, shift, next_pressure_reference,
      local_post_shift);
  const double local_post_shift_values[2U]{local_post_shift.moment,
                                           local_post_shift.absolute_moment};
  double global_post_shift_values[2U]{};
  status = reductions.checked_sum({local_post_shift_values, 2U},
                                  {global_post_shift_values, 2U}, local);
  if (!status) return status;
  const double global_post_shift_moment = global_post_shift_values[0U];
  const double global_post_shift_absolute_moment =
      global_post_shift_values[1U];
  const double post_shift_gauge_residual =
      std::abs(global_post_shift_moment) / global_weight;
  const double post_shift_gauge_tolerance =
      closed_gauge_post_shift_tolerance(
          shift, global_weight, global_post_shift_absolute_moment);
  if (!std::isfinite(global_post_shift_moment) ||
      !std::isfinite(global_post_shift_absolute_moment) ||
      global_post_shift_absolute_moment < 0.0 ||
      !std::isfinite(post_shift_gauge_residual) ||
      !std::isfinite(post_shift_gauge_tolerance) ||
      post_shift_gauge_tolerance <= 0.0 ||
      post_shift_gauge_residual > post_shift_gauge_tolerance) {
    return {StatusCode::numerical_failure, kPredictorNumerical};
  }

  const PlanFingerprint collective_transaction =
      closed_gauge_collective_transaction(
          contract, shift, next_pressure_reference, global_moment,
          global_weight, global_post_shift_moment,
          global_post_shift_absolute_moment, post_shift_gauge_residual,
          post_shift_gauge_tolerance);
  const RevisionToken rank_local_transaction =
      closed_gauge_rank_local_transaction(collective_transaction, input,
                                          local_evaluation,
                                          local_post_shift);
  PressureReferenceCertificate output = input.predecessor;
  output.time = input.time;
  output.pressure_reference = rank_local_transaction;

  ClosedGaugeCorrectionCertificate candidate;
  candidate.shift = shift;
  candidate.next_pressure_reference = next_pressure_reference;
  candidate.local_moment = local_evaluation.moment;
  candidate.local_weight = local_evaluation.weight;
  candidate.global_moment = global_moment;
  candidate.global_weight = global_weight;
  candidate.local_post_shift_moment = local_post_shift.moment;
  candidate.local_post_shift_absolute_moment =
      local_post_shift.absolute_moment;
  candidate.global_post_shift_moment = global_post_shift_moment;
  candidate.global_post_shift_absolute_moment =
      global_post_shift_absolute_moment;
  candidate.post_shift_gauge_residual = post_shift_gauge_residual;
  candidate.post_shift_gauge_tolerance = post_shift_gauge_tolerance;
  candidate.output_pressure_reference = output;
  candidate.predecessor_pressure_reference =
      input.predecessor.pressure_reference;
  candidate.time = input.time;
  candidate.geometry = input.geometry;
  candidate.pressure_correction_authority =
      input.pressure_correction_authority;
  candidate.target_thermodynamic_closure =
      input.target_thermodynamic_closure;
  candidate.activity_local_fingerprint = input.activity.local_fingerprint;
  candidate.activity_collective_fingerprint =
      input.activity.collective_fingerprint;
  candidate.collective_transaction = collective_transaction;
  candidate.rank_local_transaction = rank_local_transaction;
  candidate.local_active_cells = local_evaluation.active_cells;
  candidate.corrector = input.corrector;
  if (!candidate.valid()) {
    return {StatusCode::numerical_failure, kPredictorNumerical};
  }
  certificate = candidate;
  return {};
}

bool PressureReferencePlan::matches_closed_gauge_correction(
    const ClosedGaugeCorrectionPrepareInput& input,
    const ClosedGaugeCorrectionCertificate& certificate) const noexcept {
  if (kind_ != PressureReferenceKind::closed_mass ||
      gauge_ != PressureGauge::compressibility_weighted_zero_mean ||
      kernels_ == nullptr || fingerprint_ == 0U ||
      predictor_fingerprint_ == 0U || thermodynamics_fingerprint_ == 0U ||
      geometry_revision_ == 0U || !certificate.valid() ||
      !closed_gauge_structural_input_valid(
          input, cells_, fingerprint_, predictor_fingerprint_,
          thermodynamics_fingerprint_, geometry_revision_,
          pressure_perturbation_field_,
          compressibility_field_)) {
    return false;
  }
  ClosedGaugeLocalEvaluation local;
  if (!evaluate_closed_gauge_local(input, *kernels_, cells_, local)) {
    return false;
  }
  if (!std::isfinite(certificate.global_moment) ||
      !std::isfinite(certificate.global_weight) ||
      certificate.global_weight <= 0.0) {
    return false;
  }
  const double expected_shift =
      certificate.global_moment / certificate.global_weight;
  const double expected_reference =
      input.pressure_reference + expected_shift;
  if (!std::isfinite(expected_shift) || !std::isfinite(expected_reference) ||
      expected_reference <= 0.0 ||
      double_bits(certificate.shift) != double_bits(expected_shift) ||
      double_bits(certificate.next_pressure_reference) !=
          double_bits(expected_reference) ||
      double_bits(certificate.local_moment) != double_bits(local.moment) ||
      double_bits(certificate.local_weight) != double_bits(local.weight) ||
      certificate.local_active_cells != local.active_cells) {
    return false;
  }
  ClosedGaugePostShiftEvaluation local_post_shift;
  if (!evaluate_closed_gauge_post_shift_local(
          input, *kernels_, cells_, expected_shift, expected_reference,
          local_post_shift)) {
    return false;
  }
  const double expected_post_shift_gauge_residual =
      std::abs(certificate.global_post_shift_moment) /
      certificate.global_weight;
  const double expected_post_shift_gauge_tolerance =
      closed_gauge_post_shift_tolerance(
          expected_shift, certificate.global_weight,
          certificate.global_post_shift_absolute_moment);
  if (!std::isfinite(certificate.global_post_shift_moment) ||
      !std::isfinite(certificate.global_post_shift_absolute_moment) ||
      certificate.global_post_shift_absolute_moment < 0.0 ||
      !std::isfinite(expected_post_shift_gauge_residual) ||
      !std::isfinite(expected_post_shift_gauge_tolerance) ||
      expected_post_shift_gauge_tolerance <= 0.0 ||
      expected_post_shift_gauge_residual >
          expected_post_shift_gauge_tolerance ||
      double_bits(certificate.local_post_shift_moment) !=
          double_bits(local_post_shift.moment) ||
      double_bits(certificate.local_post_shift_absolute_moment) !=
          double_bits(local_post_shift.absolute_moment) ||
      double_bits(certificate.post_shift_gauge_residual) !=
          double_bits(expected_post_shift_gauge_residual) ||
      double_bits(certificate.post_shift_gauge_tolerance) !=
          double_bits(expected_post_shift_gauge_tolerance)) {
    return false;
  }
  const PlanFingerprint contract =
      closed_gauge_collective_contract(fingerprint_, input);
  const PlanFingerprint collective_transaction =
      closed_gauge_collective_transaction(
          contract, expected_shift, expected_reference,
          certificate.global_moment, certificate.global_weight,
          certificate.global_post_shift_moment,
          certificate.global_post_shift_absolute_moment,
          expected_post_shift_gauge_residual,
          expected_post_shift_gauge_tolerance);
  const RevisionToken rank_local_transaction =
      closed_gauge_rank_local_transaction(collective_transaction, input,
                                          local, local_post_shift);
  PressureReferenceCertificate expected_output = input.predecessor;
  expected_output.time = input.time;
  expected_output.pressure_reference = rank_local_transaction;
  return certificate.predecessor_pressure_reference ==
             input.predecessor.pressure_reference &&
         certificate.time == input.time &&
         certificate.geometry == input.geometry &&
         certificate.corrector == input.corrector &&
         certificate.pressure_correction_authority ==
             input.pressure_correction_authority &&
         certificate.target_thermodynamic_closure ==
             input.target_thermodynamic_closure &&
         certificate.activity_local_fingerprint ==
             input.activity.local_fingerprint &&
         certificate.activity_collective_fingerprint ==
             input.activity.collective_fingerprint &&
         certificate.collective_transaction == collective_transaction &&
         certificate.rank_local_transaction == rank_local_transaction &&
         same_pressure_reference_certificate(
             certificate.output_pressure_reference, expected_output);
}

Status PressureReferencePlan::certify_open_boundary(
    double pressure_reference, RevisionToken boundary_value_revision,
    const ThermophysicalPredictorCertificate& predictor,
    PressureReferenceCertificate& certificate) const noexcept {
  if (kind_ != PressureReferenceKind::boundary_absolute ||
      gauge_ != PressureGauge::absolute_boundary_dirichlet ||
      service_stage_ != 0U || fingerprint_ == 0U ||
      predictor_fingerprint_ == 0U || !predictor.valid() ||
      predictor.plan != predictor_fingerprint_ ||
      !std::isfinite(pressure_reference) || pressure_reference <= 0.0 ||
      boundary_value_revision == 0U) {
    return {StatusCode::invalid_plan, kPredictorPlan};
  }
  std::uint64_t revision = kFnvOffset;
  revision = hash_mix(revision, fingerprint_);
  revision = hash_mix(revision, predictor.state);
  revision = hash_mix(revision, predictor.time);
  revision = hash_mix(revision, boundary_value_revision);
  revision = hash_mix(revision, double_bits(pressure_reference));
  revision = revision == 0U ? 1U : revision;
  PressureReferenceCertificate candidate;
  candidate.plan = fingerprint_;
  candidate.predictor = predictor.plan;
  candidate.thermodynamics = thermodynamics_fingerprint_;
  candidate.closure = boundary_value_revision;
  candidate.time = predictor.time;
  candidate.pressure_reference = revision;
  candidate.kind = kind_;
  certificate = candidate;
  return {};
}

}  // namespace hundun::v04
