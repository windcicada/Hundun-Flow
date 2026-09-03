// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_flow.hpp"

#include <mpi.h>

#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <utility>

namespace allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void observe() noexcept {
  if (enabled.load(std::memory_order_acquire))
    count.fetch_add(1U, std::memory_order_relaxed);
}

void* allocate(std::size_t bytes) {
  observe();
  void* result = std::malloc(bytes == 0U ? 1U : bytes);
  if (result == nullptr) throw std::bad_alloc{};
  return result;
}

void* allocate_aligned(std::size_t bytes, std::size_t alignment) {
  observe();
  void* result = nullptr;
  const std::size_t requested = bytes == 0U ? alignment : bytes;
  if (posix_memalign(&result, alignment, requested) != 0 || result == nullptr)
    throw std::bad_alloc{};
  return result;
}

class Guard {
 public:
  Guard() noexcept {
    count.store(0U, std::memory_order_relaxed);
    enabled.store(true, std::memory_order_release);
  }
  ~Guard() { enabled.store(false, std::memory_order_release); }

  Guard(const Guard&) = delete;
  Guard& operator=(const Guard&) = delete;
};

}  // namespace allocation_observer

void* operator new(std::size_t bytes) {
  return allocation_observer::allocate(bytes);
}
void* operator new[](std::size_t bytes) {
  return allocation_observer::allocate(bytes);
}
void* operator new(std::size_t bytes, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      bytes, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t bytes, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      bytes, static_cast<std::size_t>(alignment));
}
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t,
                       std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;

static_assert(noexcept(select_pressure_energy_globalization(
    std::declval<const PressureEnergyGlobalizationSample&>(),
    std::declval<Span<const PressureEnergyGlobalizationSample>>(),
    std::declval<PressureEnergyGlobalizationSelectionCertificate&>())));
static_assert(noexcept(select_pressure_energy_extrapolation(
    std::declval<const PressureEnergyGlobalizationSample&>(),
    std::declval<const PressureEnergyGlobalizationSample&>(),
    std::declval<PressureEnergyGlobalizationSelectionCertificate&>())));

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

PressureEnergyGlobalizationSample make_baseline(std::uint8_t corrector) {
  PressureEnergyGlobalizationSample sample;
  sample.alpha = 0.0;
  sample.global_normalized_continuity = 1.0;
  sample.global_normalized_energy = 1.0;
  sample.thermodynamically_admissible = true;
  sample.state_and_flux_finite = true;
  sample.corrector = corrector;
  sample.target_time = 101U;
  sample.correction_direction = 201U;
  sample.state_provenance = 301U;
  sample.mass_flux_provenance = 401U;
  return sample;
}

std::array<PressureEnergyGlobalizationSample,
           kPressureEnergyGlobalizationCandidateCount>
make_candidates(const PressureEnergyGlobalizationSample& baseline) {
  std::array<PressureEnergyGlobalizationSample,
             kPressureEnergyGlobalizationCandidateCount>
      candidates{};
  for (std::size_t index = 0U; index < candidates.size(); ++index) {
    auto& sample = candidates[index];
    sample = baseline;
    sample.alpha = std::ldexp(1.0, -static_cast<int>(index));
    sample.global_normalized_continuity = 2.0;
    sample.global_normalized_energy = 2.0;
    sample.state_provenance = 1000U + index;
    sample.mass_flux_provenance = 2000U + index;
  }
  return candidates;
}

bool test_full_step_is_selected_first() {
  const auto baseline = make_baseline(1U);
  auto candidates = make_candidates(baseline);
  candidates[0U].global_normalized_continuity = 0.6;
  candidates[0U].global_normalized_energy = 0.7;
  candidates[1U].global_normalized_continuity = 0.2;
  candidates[1U].global_normalized_energy = 0.2;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  const Status status = select_pressure_energy_globalization(
      baseline, {candidates.data(), candidates.size()}, certificate);
  return expect(static_cast<bool>(status), "full-step selection succeeds") &&
         expect(certificate.valid(), "selection returns a typed certificate") &&
         expect(certificate.scope ==
                    PressureEnergyGlobalizationScope::
                        frozen_momentum_continuity_energy_globalization &&
                    certificate.corrector == 1U &&
                    certificate.target_time == baseline.target_time &&
                    certificate.correction_direction ==
                        baseline.correction_direction &&
                    certificate.alpha == 1.0 &&
                    certificate.selected_halvings == 0U &&
                    certificate.candidate_state_provenance ==
                        candidates[0U].state_provenance &&
                    certificate.candidate_mass_flux_provenance ==
                        candidates[0U].mass_flux_provenance,
                "certificate binds C1, target, direction, alpha, state and flux");
}

bool test_safeguarded_extrapolation_is_distinct_and_fail_closed() {
  const auto baseline = make_baseline(2U);
  auto extrapolated = baseline;
  extrapolated.alpha = 1.5;
  extrapolated.global_normalized_continuity = 0.4;
  extrapolated.global_normalized_energy = 0.5;
  extrapolated.state_provenance = 1001U;
  extrapolated.mass_flux_provenance = 2001U;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  Status status = select_pressure_energy_extrapolation(
      baseline, extrapolated, certificate);
  bool passed = expect(static_cast<bool>(status) && certificate.valid() &&
                           certificate.extrapolated &&
                           certificate.alpha == extrapolated.alpha &&
                           certificate.selected_halvings == 0U,
                       "an exact decreasing alpha>1 candidate is certified");

  auto forged = certificate;
  forged.extrapolated = false;
  passed &= expect(!forged.valid(),
                   "an extrapolated certificate cannot claim the legacy ladder");

  auto nondecreasing = extrapolated;
  nondecreasing.global_normalized_continuity = 1.1;
  nondecreasing.global_normalized_energy = 1.1;
  status = select_pressure_energy_extrapolation(
      baseline, nondecreasing, certificate);
  passed &= expect(status.code == StatusCode::rejected_step &&
                       !certificate.valid(),
                   "a non-decreasing extrapolation falls back without authority");

  auto out_of_range = extrapolated;
  out_of_range.alpha = 2.01;
  status = select_pressure_energy_extrapolation(
      baseline, out_of_range, certificate);
  passed &= expect(status.code == StatusCode::invalid_plan &&
                       !certificate.valid(),
                   "an extrapolation above the safeguarded cap is invalid");
  return passed;
}

bool test_evaluated_prefix_selects_without_requiring_the_frozen_tail() {
  const auto baseline = make_baseline(2U);
  auto candidates = make_candidates(baseline);
  candidates[0U].global_normalized_continuity = 1.1;
  candidates[0U].global_normalized_energy = 1.1;
  candidates[1U].global_normalized_continuity = 0.7;
  candidates[1U].global_normalized_energy = 0.6;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  const Status status = select_pressure_energy_globalization(
      baseline, {candidates.data(), 2U}, certificate);
  bool passed =
      expect(static_cast<bool>(status) && certificate.valid() &&
                 certificate.selected_halvings == 1U,
             "an evaluated prefix selects its first Armijo candidate");

  auto malformed_prefix = candidates;
  malformed_prefix[1U].alpha = malformed_prefix[0U].alpha;
  const Status malformed = select_pressure_energy_globalization(
      baseline, {malformed_prefix.data(), 2U}, certificate);
  passed &= expect(malformed.code == StatusCode::invalid_plan &&
                       !certificate.valid(),
                   "an evaluated prefix still validates every supplied rung");
  return passed;
}

bool test_non_decreasing_full_step_is_rejected_before_c2_half_step() {
  const auto baseline = make_baseline(2U);
  auto candidates = make_candidates(baseline);
  candidates[0U].global_normalized_continuity = 0.99;
  candidates[0U].global_normalized_energy = 1.01;
  candidates[1U].global_normalized_continuity = 0.8;
  candidates[1U].global_normalized_energy = 0.7;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  const Status status = select_pressure_energy_globalization(
      baseline, {candidates.data(), candidates.size()}, certificate);
  return expect(static_cast<bool>(status),
                "line search skips an admissible non-decreasing full step") &&
         expect(certificate.valid() && certificate.corrector == 2U &&
                    certificate.alpha == 0.5 &&
                    certificate.selected_halvings == 1U &&
                    certificate.candidate_state_provenance ==
                        candidates[1U].state_provenance &&
                    certificate.candidate_mass_flux_provenance ==
                        candidates[1U].mass_flux_provenance,
                "C2 half-step is the first signed decreasing candidate");
}

bool test_armijo_envelope_is_stricter_than_any_decrease() {
  const auto baseline = make_baseline(1U);
  auto candidates = make_candidates(baseline);
  candidates[0U].global_normalized_continuity = 0.99995;
  candidates[0U].global_normalized_energy = 0.99995;
  candidates[1U].global_normalized_continuity = 0.99994;
  candidates[1U].global_normalized_energy = 0.99994;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  const Status status = select_pressure_energy_globalization(
      baseline, {candidates.data(), candidates.size()}, certificate);
  return expect(static_cast<bool>(status) && certificate.valid() &&
                    certificate.selected_halvings == 1U &&
                    certificate.candidate_merit ==
                        std::hypot(0.99994, 0.99994) &&
                    certificate.candidate_merit <=
                        certificate.armijo_upper_bound,
                "strict decrease below Armijo is skipped before half step");
}

bool test_coupled_merit_accepts_physical_residual_exchange() {
  auto baseline = make_baseline(1U);
  baseline.global_normalized_continuity = 0.00417391;
  baseline.global_normalized_energy = 0.00369408;
  auto candidates = make_candidates(baseline);
  // Re3900 half-step witness: the same-target pressure--enthalpy direction
  // reduces the joint residual by about 24%. Continuity rises only 0.85%
  // while energy falls 86%; a componentwise max merit falsely rejects it.
  candidates[8U].global_normalized_continuity = 0.00420923;
  candidates[8U].global_normalized_energy = 0.000517943;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  const Status status = select_pressure_energy_globalization(
      baseline, {candidates.data(), candidates.size()}, certificate);
  return expect(static_cast<bool>(status) && certificate.valid(),
                "joint continuity-energy descent accepts residual exchange") &&
         expect(certificate.selected_halvings == 8U &&
                    certificate.alpha == std::ldexp(1.0, -8) &&
                    certificate.candidate_normalized_continuity >
                        certificate.baseline_normalized_continuity &&
                    certificate.candidate_normalized_energy <
                        certificate.baseline_normalized_energy &&
                    certificate.candidate_merit < certificate.baseline_merit,
                "certificate preserves the observed C-up/E-down coupled witness");
}

bool test_duplicate_candidate_provenance_is_not_a_ladder() {
  const auto baseline = make_baseline(1U);
  auto candidates = make_candidates(baseline);
  candidates[0U].global_normalized_continuity = 0.5;
  candidates[0U].global_normalized_energy = 0.5;
  candidates[7U].state_provenance = candidates[6U].state_provenance;
  candidates[7U].mass_flux_provenance =
      candidates[6U].mass_flux_provenance;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  const Status status = select_pressure_energy_globalization(
      baseline, {candidates.data(), candidates.size()}, certificate);
  return expect(status.code == StatusCode::invalid_plan,
                "duplicate candidate state/flux provenance is rejected") &&
         expect(!certificate.valid(),
                "duplicate sample cannot leave a selection certificate");
}

bool test_alpha_zero_baseline_and_sample_provenance_are_mandatory() {
  const auto baseline = make_baseline(1U);
  auto candidates = make_candidates(baseline);
  candidates[0U].global_normalized_continuity = 0.5;
  candidates[0U].global_normalized_energy = 0.5;
  PressureEnergyGlobalizationSelectionCertificate certificate;
  bool passed = expect(static_cast<bool>(select_pressure_energy_globalization(
                           baseline, {candidates.data(), candidates.size()},
                           certificate)) &&
                           certificate.valid(),
                       "alpha=0 baseline is accepted");

  auto fake_baseline = baseline;
  fake_baseline.alpha = std::ldexp(1.0, -24);
  Status status = select_pressure_energy_globalization(
      fake_baseline, {candidates.data(), candidates.size()}, certificate);
  passed &= expect(status.code == StatusCode::invalid_plan &&
                       !certificate.valid(),
                   "nonzero baseline alpha is rejected and clears output");

  auto malformed = candidates;
  malformed[9U].alpha = malformed[8U].alpha;
  status = select_pressure_energy_globalization(
      baseline, {malformed.data(), malformed.size()}, certificate);
  passed &= expect(status.code == StatusCode::invalid_plan &&
                       !certificate.valid(),
                   "duplicate or fake alpha is rejected");

  malformed = candidates;
  malformed[24U].target_time = baseline.target_time + 1U;
  status = select_pressure_energy_globalization(
      baseline, {malformed.data(), malformed.size()}, certificate);
  passed &= expect(status.code == StatusCode::invalid_plan &&
                       !certificate.valid(),
                   "stale target-time sample is rejected");

  malformed = candidates;
  malformed[24U].correction_direction = baseline.correction_direction + 1U;
  status = select_pressure_energy_globalization(
      baseline, {malformed.data(), malformed.size()}, certificate);
  passed &= expect(status.code == StatusCode::invalid_plan &&
                       !certificate.valid(),
                   "foreign correction direction is rejected");

  malformed = candidates;
  malformed[24U].corrector = 2U;
  status = select_pressure_energy_globalization(
      baseline, {malformed.data(), malformed.size()}, certificate);
  passed &= expect(status.code == StatusCode::invalid_plan &&
                       !certificate.valid(),
                   "foreign corrector sample is rejected");

  malformed = candidates;
  malformed[24U].state_provenance = baseline.state_provenance;
  status = select_pressure_energy_globalization(
      baseline, {malformed.data(), malformed.size()}, certificate);
  passed &= expect(status.code == StatusCode::invalid_plan &&
                       !certificate.valid(),
                   "candidate cannot reuse baseline state provenance");
  return passed;
}

bool test_positivity_alone_does_not_authorize_a_step() {
  const auto baseline = make_baseline(1U);
  auto candidates = make_candidates(baseline);
  candidates[0U].global_normalized_continuity = 0.99;
  candidates[0U].global_normalized_energy = 295.0;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  const Status status = select_pressure_energy_globalization(
      baseline, {candidates.data(), candidates.size()}, certificate);
  return expect(status.code == StatusCode::rejected_step,
                "positive finite candidate with energy merit 295 is rejected") &&
         expect(!certificate.valid(),
                "an exhausted ladder returns no clipping certificate");
}

bool test_nonfinite_and_inadmissible_candidates_are_skipped() {
  const auto baseline = make_baseline(2U);
  auto candidates = make_candidates(baseline);
  candidates[0U].global_normalized_energy =
      std::numeric_limits<double>::quiet_NaN();
  candidates[1U].global_normalized_continuity =
      std::numeric_limits<double>::infinity();
  candidates[2U].global_normalized_continuity = 0.2;
  candidates[2U].global_normalized_energy = 0.2;
  candidates[2U].state_and_flux_finite = false;
  candidates[3U].global_normalized_continuity = 0.1;
  candidates[3U].global_normalized_energy = 0.1;
  candidates[3U].thermodynamically_admissible = false;
  candidates[4U].global_normalized_continuity = 0.75;
  candidates[4U].global_normalized_energy = 0.8;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  Status status = select_pressure_energy_globalization(
      baseline, {candidates.data(), candidates.size()}, certificate);
  bool passed = expect(static_cast<bool>(status) && certificate.valid() &&
                           certificate.selected_halvings == 4U &&
                           certificate.alpha == std::ldexp(1.0, -4),
                       "NaN, Inf, nonfinite-state and nonpositive samples skip");

  auto invalid_baseline = baseline;
  invalid_baseline.global_normalized_energy =
      std::numeric_limits<double>::infinity();
  status = select_pressure_energy_globalization(
      invalid_baseline, {candidates.data(), candidates.size()}, certificate);
  passed &= expect(status.code == StatusCode::invalid_plan &&
                       !certificate.valid(),
                   "nonfinite baseline cannot define a merit policy");
  return passed;
}

bool test_last_ladder_point_and_no_full_newton_claim() {
  const auto baseline = make_baseline(1U);
  auto candidates = make_candidates(baseline);
  candidates.back().global_normalized_continuity = 0.8;
  candidates.back().global_normalized_energy = 0.9;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  const Status status = select_pressure_energy_globalization(
      baseline, {candidates.data(), candidates.size()}, certificate);
  bool passed = expect(static_cast<bool>(status) && certificate.valid() &&
                           certificate.selected_halvings == 24U &&
                           certificate.alpha == std::ldexp(1.0, -24),
                       "the fixed ladder includes 2^-24 exactly");
  auto forged = certificate;
  forged.selection_provenance ^= UINT64_C(0x100);
  passed &= expect(!forged.valid(),
                   "selection provenance is recomputed, not just nonzero");
  certificate.full_nonlinear_newton = true;
  passed &= expect(!certificate.valid(),
                   "frozen-momentum certificate cannot claim full Newton");
  return passed;
}

bool test_joint_l2_policy_has_distinct_frozen_provenance() {
  // Frozen from the original f5e0fcd L-infinity policy.  Keep this literal:
  // it is the mutation oracle proving that a merit-policy change cannot reuse
  // the prior schema identity.
  constexpr PlanFingerprint kLegacyLInfinityProvenance =
      UINT64_C(0x9ef260037a8ea7d7);
  // Frozen for the versioned "v04pegl2" policy schema and the same public
  // fixture.  A future policy change must deliberately version this oracle.
  constexpr PlanFingerprint kJointL2Provenance = UINT64_C(0xe6d46d99f1321347);

  const auto baseline = make_baseline(1U);
  auto candidates = make_candidates(baseline);
  candidates[0U].global_normalized_continuity = 0.6;
  candidates[0U].global_normalized_energy = 0.7;

  PressureEnergyGlobalizationSelectionCertificate certificate;
  const Status status = select_pressure_energy_globalization(
      baseline, {candidates.data(), candidates.size()}, certificate);
  bool passed =
      expect(static_cast<bool>(status) && certificate.valid(),
             "joint-L2 provenance fixture selects a valid full step") &&
      expect(certificate.selection_provenance == kJointL2Provenance,
             "joint-L2 policy is pinned to its versioned provenance") &&
      expect(certificate.selection_provenance != kLegacyLInfinityProvenance,
             "joint-L2 policy cannot reuse the frozen L-infinity identity");

  auto legacy_mutation = certificate;
  legacy_mutation.selection_provenance = kLegacyLInfinityProvenance;
  passed &= expect(!legacy_mutation.valid(),
                   "a frozen L-infinity provenance mutation is rejected");
  return passed;
}

bool test_selection_is_allocation_free_and_rank_independent() {
  const auto baseline = make_baseline(1U);
  auto candidates = make_candidates(baseline);
  candidates[3U].global_normalized_continuity = 0.7;
  candidates[3U].global_normalized_energy = 0.8;
  PressureEnergyGlobalizationSelectionCertificate certificate;
  Status status;
  std::size_t allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    status = select_pressure_energy_globalization(
        baseline, {candidates.data(), candidates.size()}, certificate);
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  bool passed = expect(static_cast<bool>(status) && certificate.valid() &&
                           allocations == 0U,
                       "pure selection performs no dynamic allocation");

  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  const std::uint64_t local = certificate.selection_provenance;
  const int min_status = MPI_Allreduce(&local, &minimum, 1, MPI_UINT64_T,
                                       MPI_MIN, MPI_COMM_WORLD);
  const int max_status = MPI_Allreduce(&local, &maximum, 1, MPI_UINT64_T,
                                       MPI_MAX, MPI_COMM_WORLD);
  passed &= expect(min_status == MPI_SUCCESS && max_status == MPI_SUCCESS &&
                       minimum != 0U && minimum == maximum,
                   "identical global inputs sign the same certificate on "
                   "every rank");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  bool passed = test_full_step_is_selected_first();
  passed &= test_safeguarded_extrapolation_is_distinct_and_fail_closed();
  passed &= test_evaluated_prefix_selects_without_requiring_the_frozen_tail();
  passed &= test_non_decreasing_full_step_is_rejected_before_c2_half_step();
  passed &= test_armijo_envelope_is_stricter_than_any_decrease();
  passed &= test_coupled_merit_accepts_physical_residual_exchange();
  passed &= test_duplicate_candidate_provenance_is_not_a_ladder();
  passed &= test_alpha_zero_baseline_and_sample_provenance_are_mandatory();
  passed &= test_positivity_alone_does_not_authorize_a_step();
  passed &= test_nonfinite_and_inadmissible_candidates_are_skipped();
  passed &= test_last_ladder_point_and_no_full_newton_claim();
  passed &= test_joint_l2_policy_has_distinct_frozen_provenance();
  passed &= test_selection_is_allocation_free_and_rank_independent();
  const int local = passed ? 1 : 0;
  int global = 0;
  const int reduce_status =
      MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  const int finalize_status = MPI_Finalize();
  return reduce_status == MPI_SUCCESS && finalize_status == MPI_SUCCESS &&
                 global == 1
             ? 0
             : 1;
}
