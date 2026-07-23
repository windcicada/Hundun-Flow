// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "material-density PISO test access requires tests-on build"
#endif

#include <cmath>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace hundun::flow {
class FlowState;
class FixedStepMaterialDensityFlow;
class MaterialDensityStepAttemptReport;
class PisoCoupler;
enum class FlowLayer : std::uint8_t;
struct AcceptedStepMetadata;
}
namespace hundun::runtime {
class MpiContext;
}

namespace hundun::flow::test {

struct MaterialDensityVortexSource final {
  double x{};
  double y{};
  double z{};
};

struct MaterialMomentumConservationInput final {
  double momentum_n_minus_1{};
  double momentum_n{};
  double momentum_n_plus_1{};
  double boundary_n_minus_1{};
  double boundary_n{};
  double source_n_plus_1{};
  double momentum_abs_n_minus_1{};
  double momentum_abs_n{};
  double momentum_abs_n_plus_1{};
  double boundary_abs_n_minus_1{};
  double boundary_abs_n{};
  double source_abs_n_plus_1{};
  double dt_s{};
  double alpha0{};
  double alpha2{};
  bool bdf2{};
};

struct MaterialPhaseFailureForTest final {
  bool failed{};
  std::uint8_t reason{};
  bool recoverable{};
};

struct MaterialPhaseSelectionForTest final {
  std::uint8_t reason{};
  int lowest_failing_rank{-1};
  bool recoverable{};
};

struct MaterialPressureEvidenceForTest final {
  std::vector<double> rhs_raw;
  std::vector<double> rhs_solve;
  std::vector<double> correction;
  std::vector<double> final_face_density;
  double rhs_l2{};
  std::uint32_t corrector_ordinal{};
  bool token_available{};
  bool final_operator_available{};
  std::uint64_t final_operator_revision{};
};

struct FacadeCacheSnapshot final {
  struct Workspace final {
    std::uintptr_t identity{};
    std::size_t capacity{};
  };
  struct MomentumOperator final {
    std::uintptr_t identity{};
    std::uint64_t revision{};
    std::vector<double> diagonal;
  };

  std::vector<Workspace> workspaces;
  std::array<MomentumOperator, 3> operators{};
  std::size_t operator_count{};
  bool delegated{};
};

FacadeCacheSnapshot material_facade_cache_values_for_ideal(
    const FixedStepMaterialDensityFlow &);

enum class MaterialReportCorruptionForTest : std::uint8_t {
  success_corrector_count,
  success_provenance,
  success_shared_field,
  reliable_collective_rank,
  material_reason_mapping,
  parent_transport_size,
  material_transport_size,
  unavailable_numeric_value,
  material_count_zero,
  material_count_plus_two,
  material_count_five,
  parent_transport_residual_value,
  nested_transport_residual_value,
  parent_transport_conservation_value,
  nested_transport_conservation_value,
  nested_density_residual_availability,
  nested_transport_residual_availability,
  nested_mass_conservation_availability,
  nested_transport_conservation_availability,
  nested_minimum_density_availability,
  nested_attempt_identity,
  nested_finalization_identity,
  nested_shared_field,
  nested_provenance,
  nested_residual_outer_size,
  nested_conservation_outer_size
};

enum class MaterialTerminalModeForTest : std::uint8_t {
  none,
  returned_rankless,
  thrown_operation,
  returned_reliable
};

enum class MaterialTransportAuthorityMutation : std::uint8_t {
  omitted,
  reordered,
  same_maximum_different_sequence
};

enum class MaterialTerminalPointForTest : std::uint8_t {
  predictor_stage,
  momentum_x,
  momentum_y,
  momentum_z,
  pressure_corrector_one,
  provisional_stage,
  pressure_corrector_two,
  public_finalizer,
  final_continuity_reduction,
  final_continuity_status,
  final_pressure_entry,
  final_pressure_gamma_sum,
  final_pressure_gamma_count,
  final_pressure_residual_reduction,
  final_momentum_residual_reduction,
  final_momentum_conservation_reduction,
  final_momentum_status,
  final_conservation_status,
  count
};

namespace detail {

inline std::atomic<int> material_terminal_point{-1};
inline std::atomic<int> material_terminal_mode{
    static_cast<int>(MaterialTerminalModeForTest::none)};
inline std::array<std::atomic<std::uint64_t>,
                  static_cast<std::size_t>(
                      MaterialTerminalPointForTest::count)>
    material_terminal_calls{};
inline std::atomic<std::uint32_t> face_flux_path_observation_depth{};

class FaceFluxPathObservation final {
public:
  explicit FaceFluxPathObservation(bool enabled = true) noexcept
      : enabled_(enabled) {
    if (enabled_)
      face_flux_path_observation_depth.fetch_add(1U,
                                                  std::memory_order_relaxed);
  }
  ~FaceFluxPathObservation() noexcept {
    if (enabled_)
      face_flux_path_observation_depth.fetch_sub(1U,
                                                  std::memory_order_relaxed);
  }
  FaceFluxPathObservation(const FaceFluxPathObservation &) = delete;
  FaceFluxPathObservation &operator=(const FaceFluxPathObservation &) = delete;

private:
  bool enabled_{};
};

inline MaterialTerminalModeForTest
reach_material_terminal_point(MaterialTerminalPointForTest point) noexcept {
  material_terminal_calls[static_cast<std::size_t>(point)].fetch_add(
      1U, std::memory_order_relaxed);
  if (material_terminal_point.load(std::memory_order_relaxed) !=
      static_cast<int>(point))
    return MaterialTerminalModeForTest::none;
  return static_cast<MaterialTerminalModeForTest>(
      material_terminal_mode.load(std::memory_order_relaxed));
}

} // namespace detail

class MaterialDensityPisoTestAccess final {
public:
  static constexpr std::uint32_t contract_version() noexcept { return 1U; }

  static double momentum_conservation_defect(
      const MaterialMomentumConservationInput &) noexcept;
  static bool report_authenticated(
      const MaterialDensityStepAttemptReport &) noexcept;
  static void corrupt_report(MaterialDensityStepAttemptReport &,
                             MaterialReportCorruptionForTest);
  static void set_preflight_allocation_failure_rank(int) noexcept;
  static void reset_preflight_allocation_failure() noexcept;
  static void set_terminal_fault(MaterialTerminalPointForTest point,
                                 MaterialTerminalModeForTest mode) noexcept {
    detail::material_terminal_point.store(static_cast<int>(point),
                                          std::memory_order_relaxed);
    detail::material_terminal_mode.store(static_cast<int>(mode),
                                         std::memory_order_relaxed);
  }
  static void reset_terminal_fault() noexcept {
    detail::material_terminal_point.store(-1, std::memory_order_relaxed);
    detail::material_terminal_mode.store(
        static_cast<int>(MaterialTerminalModeForTest::none),
        std::memory_order_relaxed);
    for (auto &calls : detail::material_terminal_calls)
      calls.store(0U, std::memory_order_relaxed);
  }
  static std::uint64_t
  terminal_point_calls(MaterialTerminalPointForTest point) noexcept {
    return detail::material_terminal_calls[static_cast<std::size_t>(point)]
        .load(std::memory_order_relaxed);
  }
  static MaterialPhaseSelectionForTest select_phase_failure(
      const runtime::MpiContext &, MaterialPhaseFailureForTest);
  static void require_reliable_collective_result(std::uint8_t reason,
                                                 int lowest_failing_rank);
  static MaterialPressureEvidenceForTest
  material_pressure_evidence(const FixedStepMaterialDensityFlow &);
  static FacadeCacheSnapshot
  facade_cache_snapshot(const FixedStepMaterialDensityFlow &);
  static MaterialPressureEvidenceForTest
  material_pressure_evidence(const PisoCoupler &);
  static const std::vector<double> &
  finalizer_flux_evidence(const FixedStepMaterialDensityFlow &);

  static MaterialDensityVortexSource vortex_source(double x, double y,
                                                    double mu) noexcept {
    const double density = 1.0 + 0.1 * std::sin(x) * std::sin(y);
    return {density * std::sin(x) * std::cos(x) +
                2.0 * mu * std::sin(x) * std::cos(y),
            density * std::sin(y) * std::cos(y) -
                2.0 * mu * std::cos(x) * std::sin(y),
            0.0};
  }

  static void force_state_diagnostic_identity(FlowState &,
                                               std::uint64_t) noexcept;
  static std::uint64_t
  state_diagnostic_identity(const FlowState &) noexcept;
  static bool state_attempt_active(const FlowState &) noexcept;
  static std::uint64_t state_attempt_identity(const FlowState &) noexcept;
  static std::uint64_t state_allocation_identity(const FlowState &,
                                                 FlowLayer);
  static void set_accepted_face_mass_flux(FlowState &, std::size_t face,
                                          double value);
  static void force_state_metadata(FlowState &,
                                   AcceptedStepMetadata) noexcept;
  static bool face_flux_path_observation_active() noexcept;
  static void force_flow_attempt_identity(FixedStepMaterialDensityFlow &,
                                          std::uint64_t) noexcept;
  static bool has_diagnostic_report(
      const FixedStepMaterialDensityFlow &) noexcept;
  static void enable_vortex_source(FixedStepMaterialDensityFlow &,
                                   bool) noexcept;
  static void mutate_transport_authority(
      FixedStepMaterialDensityFlow &,
      MaterialTransportAuthorityMutation) noexcept;
  static double material_face_value(
      double predictor_velocity_interpolation,
      double predictor_momentum_interpolation, double face_density,
      double cell_momentum_n, double face_momentum_n,
      double cell_momentum_n_minus_1, double face_momentum_n_minus_1,
      double mobility, double dt_s, double alpha1, double alpha2,
      double pressure_correction) noexcept;
};

} // namespace hundun::flow::test
