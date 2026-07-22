// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/ideal_gas_closure.hpp"

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include "ideal_gas_closure_test_access.hpp"
#endif

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace hundun::flow {
namespace {

constexpr runtime::PhaseId kStatePhase = 1800U;
constexpr runtime::ActorId kStateActor = 1800U;
constexpr std::uint64_t kClosureReportSeed = 0x696465616c676173ULL;

std::uint64_t fp_bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

void mix(std::uint64_t &hash, std::uint64_t value) noexcept {
  hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6U) + (hash >> 2U);
}

bool same(runtime::Int3 left, runtime::Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

double relative_error(double actual, double expected) noexcept {
  return std::abs(actual - expected) /
         std::max({std::abs(actual), std::abs(expected), DBL_MIN});
}

double relative_product_error(double a, double b, double c,
                              double expected) noexcept {
  const long double actual_value = static_cast<long double>(a) * b * c;
  const long double expected_value = expected;
  const long double denominator =
      std::max({std::abs(actual_value), std::abs(expected_value),
                std::numeric_limits<long double>::min()});
  const long double value =
      std::abs(actual_value - expected_value) / denominator;
  return std::isfinite(value) &&
                 value <= static_cast<long double>(
                              std::numeric_limits<double>::max())
             ? static_cast<double>(value)
             : std::numeric_limits<double>::infinity();
}

double checked_ratio(long double numerator, long double denominator,
                     const char *message) {
  if (!(denominator > 0.0L) || !std::isfinite(denominator))
    throw runtime::Error(message);
  const long double value = numerator / denominator;
  if (!(value > 0.0L) || !std::isfinite(value) ||
      value > static_cast<long double>(std::numeric_limits<double>::max()) ||
      value <
          static_cast<long double>(std::numeric_limits<double>::denorm_min()))
    throw runtime::Error(message);
  return static_cast<double>(value);
}

class CompensatedSum final {
public:
  void add(double value) noexcept {
    const double adjusted = value - correction_;
    const double next = sum_ + adjusted;
    correction_ = (next - sum_) - adjusted;
    sum_ = next;
  }
  double value() const noexcept { return sum_; }

private:
  double sum_{};
  double correction_{};
};

struct FailureSelection final {
  IdealGasClosureFailureReason reason{IdealGasClosureFailureReason::none};
  int rank{-1};
};

FailureSelection agree_failure(const runtime::MpiContext &mpi,
                               IdealGasClosureFailureReason local) {
  const int candidate =
      local == IdealGasClosureFailureReason::none ? mpi.size() : mpi.rank();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(ideal-gas closure failure rank)");
  if (lowest == mpi.size())
    return {};
  int encoded = mpi.rank() == lowest ? static_cast<int>(local) : 0;
  runtime::check_mpi_result(MPI_Bcast(&encoded, 1, MPI_INT, lowest, mpi.comm()),
                            "MPI_Bcast(ideal-gas closure failure reason)");
  if (encoded <= static_cast<int>(IdealGasClosureFailureReason::none) ||
      encoded >
          static_cast<int>(IdealGasClosureFailureReason::collective_operation))
    throw runtime::Error("ideal-gas closure failure payload is invalid");
  return {static_cast<IdealGasClosureFailureReason>(encoded), lowest};
}

struct PreflightWire final {
  std::uint64_t mode{};
  std::uint64_t enthalpy{};
  std::uint64_t density{};
  std::uint64_t face_flux{};
  std::uint64_t cp{};
  std::uint64_t gas_constant{};
  std::uint64_t configured_pressure{};
  std::uint64_t global_cells[3]{};
  std::uint64_t attempt_identity{};
  std::uint64_t local_cells[3]{};
  std::int64_t owned_begin[3]{};
  std::int64_t owned_end[3]{};
  std::uint64_t local_faces{};
};

void agree_preflight(const runtime::MpiContext &mpi, const PreflightWire &local,
                     const char *operation) {
  std::vector<PreflightWire> all(static_cast<std::size_t>(mpi.size()));
  runtime::check_mpi_result(
      MPI_Allgather(&local, static_cast<int>(sizeof(local)), MPI_BYTE,
                    all.data(), static_cast<int>(sizeof(local)), MPI_BYTE,
                    mpi.comm()),
      operation);
  constexpr std::size_t common_bytes = offsetof(PreflightWire, local_cells);
  if (std::any_of(all.begin() + 1, all.end(), [&](const auto &candidate) {
        return std::memcmp(&candidate, &all.front(), common_bytes) != 0;
      }))
    throw runtime::Error("ideal-gas closure rank preflight disagrees");
  const auto &common = all.front();
  if (common.global_cells[0] == 0U || common.global_cells[1] == 0U ||
      common.global_cells[2] == 0U)
    throw runtime::Error("ideal-gas closure global extent is invalid");
  std::uint64_t covered = 0U;
  for (std::size_t left = 0; left < all.size(); ++left) {
    const auto &candidate = all[left];
    std::uint64_t volume = 1U;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if (candidate.owned_begin[axis] < 0 ||
          candidate.owned_end[axis] <= candidate.owned_begin[axis] ||
          static_cast<std::uint64_t>(candidate.owned_end[axis]) >
              common.global_cells[axis] ||
          static_cast<std::uint64_t>(candidate.owned_end[axis] -
                                     candidate.owned_begin[axis]) !=
              candidate.local_cells[axis] ||
          candidate.local_cells[axis] == 0U ||
          volume > std::numeric_limits<std::uint64_t>::max() /
                       candidate.local_cells[axis])
        throw runtime::Error("ideal-gas closure ownership is invalid");
      volume *= candidate.local_cells[axis];
    }
    if (candidate.local_faces == 0U ||
        covered > std::numeric_limits<std::uint64_t>::max() - volume)
      throw runtime::Error("ideal-gas closure local topology is invalid");
    covered += volume;
    for (std::size_t right = 0; right < left; ++right) {
      bool overlaps = true;
      for (std::size_t axis = 0; axis < 3U; ++axis)
        overlaps = overlaps &&
                   candidate.owned_begin[axis] < all[right].owned_end[axis] &&
                   all[right].owned_begin[axis] < candidate.owned_end[axis];
      if (overlaps)
        throw runtime::Error("ideal-gas closure ownership overlaps");
    }
  }
  std::uint64_t global_volume = 1U;
  for (const auto extent : common.global_cells) {
    if (global_volume > std::numeric_limits<std::uint64_t>::max() / extent)
      throw runtime::Error("ideal-gas closure global extent overflows");
    global_volume *= extent;
  }
  if (covered != global_volume)
    throw runtime::Error("ideal-gas closure ownership does not exactly cover");
}

std::size_t cell_count(runtime::Int3 extent) {
  return static_cast<std::size_t>(extent.x) *
         static_cast<std::size_t>(extent.y) *
         static_cast<std::size_t>(extent.z);
}

template <class Function>
void for_each_cell(runtime::Int3 extent, Function &&function) {
  std::size_t offset = 0U;
  for (int k = 0; k < extent.z; ++k)
    for (int j = 0; j < extent.y; ++j)
      for (int i = 0; i < extent.x; ++i)
        function(i, j, k, offset++);
}

IdealGasClosureFailureReason first_state_failure(double rho, double q,
                                                 double cp) noexcept {
  if (!std::isfinite(q))
    return IdealGasClosureFailureReason::non_finite_enthalpy;
  const double h = q / rho;
  if (!std::isfinite(h))
    return IdealGasClosureFailureReason::non_finite_enthalpy;
  if (!(h > 0.0))
    return IdealGasClosureFailureReason::non_positive_enthalpy;
  const double temperature = h / cp;
  if (!std::isfinite(temperature))
    return IdealGasClosureFailureReason::non_finite_temperature;
  if (!(temperature > 0.0))
    return IdealGasClosureFailureReason::non_positive_temperature;
  if (!std::isfinite(rho))
    return IdealGasClosureFailureReason::non_finite_density;
  if (!(rho > 0.0))
    return IdealGasClosureFailureReason::non_positive_density;
  return IdealGasClosureFailureReason::none;
}

IdealGasClosureFailureReason
lower_failure(IdealGasClosureFailureReason left,
              IdealGasClosureFailureReason right) noexcept {
  if (left == IdealGasClosureFailureReason::none)
    return right;
  if (right == IdealGasClosureFailureReason::none)
    return left;
  return static_cast<std::uint8_t>(left) < static_cast<std::uint8_t>(right)
             ? left
             : right;
}

} // namespace

struct IdealGasClosure::Impl final {
  const mesh::MeshTopology *topology{};
  const mesh::MeshGeometry *geometry{};
  const boundary::BoundaryRegistry *boundaries{};
  const runtime::MpiContext *mpi{};
  const runtime::FieldRegistry *registry{};
  FlowFieldIds fields;
  IdealGasClosureSpec spec;
  IdealGasClosureState committed;
  IdealGasClosureState trial;
  IdealGasClosureState prepared;
  bool active{};
  bool prepared_valid{};
  IdealGasClosureReport latest;
  bool latest_available{};
  std::uint64_t attempt_identity{};
  std::uint64_t source_generation{1U};
};

IdealGasClosureDisposition IdealGasClosureReport::disposition() const noexcept {
  return disposition_;
}
IdealGasClosureFailureReason IdealGasClosureReport::reason() const noexcept {
  return reason_;
}
IdealGasClosureStage IdealGasClosureReport::stage() const noexcept {
  return stage_;
}
int IdealGasClosureReport::lowest_failing_rank() const noexcept {
  return lowest_failing_rank_;
}
std::uint64_t IdealGasClosureReport::attempt_identity() const noexcept {
  return attempt_identity_;
}
std::uint32_t IdealGasClosureReport::evaluation_count() const noexcept {
  return evaluation_count_;
}
std::uint64_t IdealGasClosureReport::collective_count() const noexcept {
  return collective_count_;
}
IdealGasPressureMode IdealGasClosureReport::pressure_mode() const noexcept {
  return pressure_mode_;
}
double IdealGasClosureReport::configured_pressure_pa() const noexcept {
  return configured_pressure_pa_;
}
bool IdealGasClosureReport::candidate_pressure_available() const noexcept {
  return candidate_pressure_available_;
}
double IdealGasClosureReport::candidate_pressure_pa() const {
  if (!candidate_pressure_available_)
    throw runtime::Error("ideal-gas candidate pressure is unavailable");
  return candidate_pressure_pa_;
}
bool IdealGasClosureReport::target_mass_available() const noexcept {
  return target_mass_kg_.has_value();
}
double IdealGasClosureReport::target_mass_kg() const {
  if (!target_mass_kg_)
    throw runtime::Error("ideal-gas target mass is unavailable");
  return *target_mass_kg_;
}
bool IdealGasClosureReport::final_metrics_available() const noexcept {
  return final_metrics_available_;
}

#define HUNDUN_CLOSURE_METRIC_ACCESSOR(method, member)                         \
  double IdealGasClosureReport::method() const {                               \
    if (!final_metrics_available_)                                             \
      throw runtime::Error("ideal-gas final closure metrics are unavailable"); \
    return member;                                                             \
  }
HUNDUN_CLOSURE_METRIC_ACCESSOR(actual_mass_kg, actual_mass_kg_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(temperature_min_K, temperature_min_K_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(temperature_max_K, temperature_max_K_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(enthalpy_min_J_per_kg, enthalpy_min_J_per_kg_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(enthalpy_max_J_per_kg, enthalpy_max_J_per_kg_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(density_min_kg_per_m3, density_min_kg_per_m3_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(density_max_kg_per_m3, density_max_kg_per_m3_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(rho_remap_normalized_l2,
                               rho_remap_normalized_l2_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(rho_h_remap_normalized_l2,
                               rho_h_remap_normalized_l2_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(rho_remap_relative_conservation_defect,
                               rho_remap_relative_conservation_defect_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(rho_h_remap_relative_conservation_defect,
                               rho_h_remap_relative_conservation_defect_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(enthalpy_temperature_max_relative_error,
                               enthalpy_temperature_max_relative_error_)
HUNDUN_CLOSURE_METRIC_ACCESSOR(eos_max_relative_error, eos_max_relative_error_)
#undef HUNDUN_CLOSURE_METRIC_ACCESSOR

std::uint64_t IdealGasClosureReport::compute_seal() const noexcept {
  std::uint64_t hash = kClosureReportSeed;
  mix(hash, static_cast<std::uint64_t>(disposition_));
  mix(hash, static_cast<std::uint64_t>(reason_));
  mix(hash, static_cast<std::uint64_t>(stage_));
  mix(hash, static_cast<std::uint64_t>(lowest_failing_rank_ + 1));
  mix(hash, attempt_identity_);
  mix(hash, evaluation_count_);
  mix(hash, collective_count_);
  mix(hash, static_cast<std::uint64_t>(pressure_mode_));
  mix(hash, fp_bits(configured_pressure_pa_));
  mix(hash, candidate_pressure_available_ ? 1U : 0U);
  mix(hash, fp_bits(candidate_pressure_pa_));
  mix(hash, target_mass_kg_.has_value() ? 1U : 0U);
  if (target_mass_kg_)
    mix(hash, fp_bits(*target_mass_kg_));
  mix(hash, final_metrics_available_ ? 1U : 0U);
  for (const double value :
       {actual_mass_kg_, temperature_min_K_, temperature_max_K_,
        enthalpy_min_J_per_kg_, enthalpy_max_J_per_kg_, density_min_kg_per_m3_,
        density_max_kg_per_m3_, rho_remap_normalized_l2_,
        rho_h_remap_normalized_l2_, rho_remap_relative_conservation_defect_,
        rho_h_remap_relative_conservation_defect_,
        enthalpy_temperature_max_relative_error_, eos_max_relative_error_})
    mix(hash, fp_bits(value));
  return hash;
}
bool IdealGasClosureReport::semantic_valid() const noexcept {
  const auto finite_nonnegative = [](double value) noexcept {
    return std::isfinite(value) && value >= 0.0;
  };
  if (static_cast<std::uint8_t>(disposition_) >
          static_cast<std::uint8_t>(
              IdealGasClosureDisposition::non_retryable_failure) ||
      static_cast<std::uint8_t>(reason_) >
          static_cast<std::uint8_t>(
              IdealGasClosureFailureReason::collective_operation) ||
      static_cast<std::uint8_t>(stage_) >
          static_cast<std::uint8_t>(IdealGasClosureStage::final) ||
      static_cast<std::uint8_t>(pressure_mode_) >
          static_cast<std::uint8_t>(IdealGasPressureMode::open_fixed) ||
      attempt_identity_ == 0U || evaluation_count_ == 0U ||
      evaluation_count_ > 3U ||
      static_cast<std::uint8_t>(stage_) != evaluation_count_ ||
      !(configured_pressure_pa_ > 0.0) ||
      !std::isfinite(configured_pressure_pa_) ||
      target_mass_kg_.has_value() !=
          (pressure_mode_ == IdealGasPressureMode::closed_dynamic) ||
      (target_mass_kg_ &&
       (!(*target_mass_kg_ > 0.0) || !std::isfinite(*target_mass_kg_))) ||
      (candidate_pressure_available_ &&
       (!(candidate_pressure_pa_ > 0.0) ||
        !std::isfinite(candidate_pressure_pa_))) ||
      (!candidate_pressure_available_ && fp_bits(candidate_pressure_pa_) != 0U))
    return false;
  const bool prefix_valid =
      (evaluation_count_ == 1U &&
       (collective_count_ == 3U || collective_count_ == 5U)) ||
      (evaluation_count_ == 2U &&
       (collective_count_ == 7U || collective_count_ == 9U)) ||
      (evaluation_count_ == 3U &&
       (collective_count_ == 11U || collective_count_ == 13U ||
        collective_count_ == 14U));
  if (!prefix_valid ||
      (candidate_pressure_available_ != (collective_count_ >= 5U)) ||
      (final_metrics_available_ && collective_count_ != 14U))
    return false;
  if (final_metrics_available_) {
    if (!(actual_mass_kg_ > 0.0) || !std::isfinite(actual_mass_kg_) ||
        !(temperature_min_K_ > 0.0) || !std::isfinite(temperature_min_K_) ||
        temperature_max_K_ < temperature_min_K_ ||
        !std::isfinite(temperature_max_K_) || !(enthalpy_min_J_per_kg_ > 0.0) ||
        !std::isfinite(enthalpy_min_J_per_kg_) ||
        enthalpy_max_J_per_kg_ < enthalpy_min_J_per_kg_ ||
        !std::isfinite(enthalpy_max_J_per_kg_) ||
        !(density_min_kg_per_m3_ > 0.0) ||
        !std::isfinite(density_min_kg_per_m3_) ||
        density_max_kg_per_m3_ < density_min_kg_per_m3_ ||
        !std::isfinite(density_max_kg_per_m3_) ||
        !finite_nonnegative(rho_remap_normalized_l2_) ||
        !finite_nonnegative(rho_h_remap_normalized_l2_) ||
        !finite_nonnegative(rho_remap_relative_conservation_defect_) ||
        !finite_nonnegative(rho_h_remap_relative_conservation_defect_) ||
        !finite_nonnegative(enthalpy_temperature_max_relative_error_) ||
        !finite_nonnegative(eos_max_relative_error_))
      return false;
  }
  if (disposition_ == IdealGasClosureDisposition::closed)
    return reason_ == IdealGasClosureFailureReason::none &&
           lowest_failing_rank_ == -1 &&
           (stage_ != IdealGasClosureStage::final || final_metrics_available_);
  return reason_ != IdealGasClosureFailureReason::none &&
         lowest_failing_rank_ >= 0;
}
void IdealGasClosureReport::seal() noexcept { seal_ = compute_seal(); }
bool IdealGasClosureReport::authenticated() const noexcept {
  return seal_ != 0U && seal_ == compute_seal() && semantic_valid();
}

IdealGasClosure::IdealGasClosure(std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
IdealGasClosure::~IdealGasClosure() noexcept = default;
IdealGasClosure::IdealGasClosure(IdealGasClosure &&other) noexcept
    : impl_(std::move(other.impl_)) {
  if (impl_) {
    if (impl_->source_generation == std::numeric_limits<std::uint64_t>::max())
      impl_->source_generation = 0U;
    else
      ++impl_->source_generation;
  }
}

IdealGasClosure IdealGasClosure::create(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::MpiContext &mpi, const runtime::FieldRegistry &registry,
    const FlowFieldIds &fields, const FlowState &state,
    IdealGasClosureSpec spec) {
  bool geometry_valid = true;
  try {
    geometry.require_compatible(topology);
  } catch (...) {
    geometry_valid = false;
  }
  const auto mode = boundaries.open_domain()
                        ? IdealGasPressureMode::open_fixed
                        : IdealGasPressureMode::closed_dynamic;
  const auto local_extent = state.layer(FlowLayer::committed).interior_extent();
  const auto global = topology.global_extent();
  const auto owned = topology.owned_global_box();
  const PreflightWire wire{static_cast<std::uint64_t>(mode),
                           spec.enthalpy_density,
                           fields.density,
                           fields.face_mass_flux,
                           fp_bits(spec.cp_J_per_kg_K),
                           fp_bits(spec.gas_constant_J_per_kg_K),
                           fp_bits(spec.configured_thermodynamic_pressure_pa),
                           {static_cast<std::uint64_t>(global.x),
                            static_cast<std::uint64_t>(global.y),
                            static_cast<std::uint64_t>(global.z)},
                           0U,
                           {static_cast<std::uint64_t>(local_extent.x),
                            static_cast<std::uint64_t>(local_extent.y),
                            static_cast<std::uint64_t>(local_extent.z)},
                           {owned.begin.x, owned.begin.y, owned.begin.z},
                           {owned.end.x, owned.end.y, owned.end.z},
                           topology.local_face_count()};
  agree_preflight(mpi, wire,
                  "MPI_Allgather(ideal-gas closure create preflight)");

  bool valid =
      geometry_valid && !state.attempt_active() &&
      &state.solver_registry() == &registry &&
      state.fields().density == fields.density &&
      state.fields().face_mass_flux == fields.face_mass_flux &&
      state.fields().transported_cell_fields ==
          fields.transported_cell_fields &&
      !fields.transported_cell_fields.empty() &&
      fields.transported_cell_fields.front() == spec.enthalpy_density &&
      spec.cp_J_per_kg_K > 0.0 && std::isfinite(spec.cp_J_per_kg_K) &&
      spec.gas_constant_J_per_kg_K > 0.0 &&
      std::isfinite(spec.gas_constant_J_per_kg_K) &&
      spec.configured_thermodynamic_pressure_pa > 0.0 &&
      std::isfinite(spec.configured_thermodynamic_pressure_pa) &&
      same(local_extent, {topology.owned_global_box().end.x -
                              topology.owned_global_box().begin.x,
                          topology.owned_global_box().end.y -
                              topology.owned_global_box().begin.y,
                          topology.owned_global_box().end.z -
                              topology.owned_global_box().begin.z});
  try {
    const auto &rho_descriptor = registry.descriptor(fields.density);
    const auto &h_descriptor = registry.descriptor(spec.enthalpy_density);
    valid = valid && registry.frozen() &&
            rho_descriptor.space == runtime::FunctionSpace::cell_average &&
            rho_descriptor.scalar_type == runtime::ScalarType::float64 &&
            rho_descriptor.components == 1U &&
            rho_descriptor.ghost_width >= 2 && rho_descriptor.conservative &&
            rho_descriptor.unit == "kg/m3" &&
            h_descriptor.space == runtime::FunctionSpace::cell_average &&
            h_descriptor.scalar_type == runtime::ScalarType::float64 &&
            h_descriptor.components == 1U && h_descriptor.ghost_width >= 2 &&
            h_descriptor.conservative && h_descriptor.unit == "J/m3";
  } catch (...) {
    valid = false;
  }
  const auto committed = state.snapshot(FlowLayer::committed);
  const auto history = state.snapshot(FlowLayer::history);
  const std::size_t enthalpy_index = 0U;
  valid = valid && committed.density.size() == topology.owned_cell_count() &&
          history.density.size() == topology.owned_cell_count() &&
          committed.transported_cell_fields.size() > enthalpy_index &&
          history.transported_cell_fields.size() > enthalpy_index;
  std::array<CompensatedSum, 4> local_sums{};
  if (valid) {
    for (std::size_t cell = 0; cell < topology.owned_cell_count(); ++cell) {
      const double volume = geometry.cell_volume_m3(cell);
      const double rho = committed.density[cell];
      const double rho_history = history.density[cell];
      const double q = committed.transported_cell_fields[enthalpy_index][cell];
      const double q_history =
          history.transported_cell_fields[enthalpy_index][cell];
      const auto current_failure =
          first_state_failure(rho, q, spec.cp_J_per_kg_K);
      const auto history_failure =
          first_state_failure(rho_history, q_history, spec.cp_J_per_kg_K);
      if (current_failure != IdealGasClosureFailureReason::none ||
          history_failure != IdealGasClosureFailureReason::none ||
          !(volume > 0.0) || !std::isfinite(volume)) {
        valid = false;
        break;
      }
      const double temperature = (q / rho) / spec.cp_J_per_kg_K;
      const double temperature_history =
          (q_history / rho_history) / spec.cp_J_per_kg_K;
      if (mode == IdealGasPressureMode::closed_dynamic)
        local_sums[0].add(volume * rho);
      local_sums[1].add(volume / temperature);
      if (mode == IdealGasPressureMode::closed_dynamic)
        local_sums[2].add(volume * rho_history);
      local_sums[3].add(volume / temperature_history);
    }
  }
  if (valid && mode == IdealGasPressureMode::open_fixed) {
    const auto inlet = boundaries.velocity_inlet_patch_id();
    valid = inlet.has_value();
    if (valid) {
      const auto &inlet_state = boundaries.patch(*inlet).inlet_state();
      valid = inlet_state.has_value() && inlet_state->temperature_K.has_value();
      const double h = valid ? inlet_state->enthalpy_J_per_kg : 0.0;
      const double rho = valid ? inlet_state->density_kg_per_m3 : 0.0;
      const double temperature = valid ? *inlet_state->temperature_K : 0.0;
      valid = valid && h > 0.0 && rho > 0.0 && temperature > 0.0 &&
              std::isfinite(h) && std::isfinite(rho) &&
              std::isfinite(temperature) &&
              relative_error(h, spec.cp_J_per_kg_K * temperature) <= 1.0e-12 &&
              relative_product_error(
                  rho, spec.gas_constant_J_per_kg_K, temperature,
                  spec.configured_thermodynamic_pressure_pa) <= 1.0e-12;
    }
  }
  const auto validation =
      agree_failure(mpi, valid ? IdealGasClosureFailureReason::none
                               : IdealGasClosureFailureReason::invalid_input);
  if (validation.reason != IdealGasClosureFailureReason::none)
    throw runtime::Error("ideal-gas closure initial state is invalid");

  std::array<double, 4> sums{};
  std::transform(local_sums.begin(), local_sums.end(), sums.begin(),
                 [](const auto &sum) { return sum.value(); });
  mpi.allreduce_fp64_in_place(sums.data(), sums.size(),
                              runtime::Fp64ReductionOperation::sum);
  const double target_mass = sums[0];
  double current_pressure = spec.configured_thermodynamic_pressure_pa;
  double history_pressure = spec.configured_thermodynamic_pressure_pa;
  if (mode == IdealGasPressureMode::closed_dynamic) {
    current_pressure = checked_ratio(
        static_cast<long double>(target_mass) * spec.gas_constant_J_per_kg_K,
        sums[1], "ideal-gas initial pressure is invalid");
    history_pressure = checked_ratio(
        static_cast<long double>(target_mass) * spec.gas_constant_J_per_kg_K,
        sums[3], "ideal-gas historical pressure is invalid");
  }
  std::array<double, 10> maxima{};
  std::fill_n(maxima.begin(), 6U, -std::numeric_limits<double>::infinity());
  for (std::size_t cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const double rho = committed.density[cell];
    const double rho_history = history.density[cell];
    const double h = committed.transported_cell_fields[0][cell] / rho;
    const double hh = history.transported_cell_fields[0][cell] / rho_history;
    const double temperature = h / spec.cp_J_per_kg_K;
    const double th = hh / spec.cp_J_per_kg_K;
    maxima[0] = std::max(maxima[0], -rho);
    maxima[1] = std::max(maxima[1], -rho_history);
    maxima[2] = std::max(maxima[2], -h);
    maxima[3] = std::max(maxima[3], -hh);
    maxima[4] = std::max(maxima[4], -temperature);
    maxima[5] = std::max(maxima[5], -th);
    maxima[6] = std::max(maxima[6],
                         relative_error(h, spec.cp_J_per_kg_K * temperature));
    maxima[7] =
        std::max(maxima[7], relative_error(hh, spec.cp_J_per_kg_K * th));
    maxima[8] = std::max(
        maxima[8], relative_product_error(rho, spec.gas_constant_J_per_kg_K,
                                          temperature, current_pressure));
    maxima[9] =
        std::max(maxima[9], relative_product_error(rho_history,
                                                   spec.gas_constant_J_per_kg_K,
                                                   th, history_pressure));
  }
  mpi.allreduce_fp64_in_place(maxima.data(), maxima.size(),
                              runtime::Fp64ReductionOperation::maximum);
  const bool consistent =
      sums[1] > 0.0 && sums[3] > 0.0 &&
      (mode == IdealGasPressureMode::open_fixed ||
       (target_mass > 0.0 && std::isfinite(target_mass) &&
        relative_error(sums[2], target_mass) <= 5.0e-12)) &&
      relative_error(current_pressure,
                     spec.configured_thermodynamic_pressure_pa) <= 1.0e-12 &&
      maxima[6] <= 1.0e-12 && maxima[7] <= 1.0e-12 && maxima[8] <= 1.0e-12 &&
      maxima[9] <= 1.0e-12;
  if (!consistent)
    throw runtime::Error("ideal-gas closure initial EOS is inconsistent");

  auto impl = std::make_unique<Impl>();
  impl->topology = &topology;
  impl->geometry = &geometry;
  impl->boundaries = &boundaries;
  impl->mpi = &mpi;
  impl->registry = &registry;
  impl->fields = fields;
  impl->spec = spec;
  impl->committed = {mode, current_pressure,
                     mode == IdealGasPressureMode::closed_dynamic
                         ? std::optional<double>(target_mass)
                         : std::nullopt,
                     0U};
  impl->trial = impl->committed;
  impl->prepared = impl->committed;
  return IdealGasClosure(std::move(impl));
}

IdealGasClosureState IdealGasClosure::state() const {
  if (!impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  return impl_->committed;
}

void IdealGasClosure::begin_attempt(const FlowState &state,
                                    std::uint64_t identity) {
  if (!impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  if (impl_->active)
    throw runtime::Error("ideal-gas closure call overlaps another call");
  if (!state.attempt_active() || identity == 0U)
    throw runtime::Error("ideal-gas closure attempt is invalid");
  const auto local_extent = state.layer(FlowLayer::committed).interior_extent();
  const auto global = impl_->topology->global_extent();
  const auto owned = impl_->topology->owned_global_box();
  const PreflightWire wire{
      static_cast<std::uint64_t>(impl_->committed.mode),
      impl_->spec.enthalpy_density,
      impl_->fields.density,
      impl_->fields.face_mass_flux,
      fp_bits(impl_->spec.cp_J_per_kg_K),
      fp_bits(impl_->spec.gas_constant_J_per_kg_K),
      fp_bits(impl_->spec.configured_thermodynamic_pressure_pa),
      {static_cast<std::uint64_t>(global.x),
       static_cast<std::uint64_t>(global.y),
       static_cast<std::uint64_t>(global.z)},
      identity,
      {static_cast<std::uint64_t>(local_extent.x),
       static_cast<std::uint64_t>(local_extent.y),
       static_cast<std::uint64_t>(local_extent.z)},
      {owned.begin.x, owned.begin.y, owned.begin.z},
      {owned.end.x, owned.end.y, owned.end.z},
      impl_->topology->local_face_count()};
  agree_preflight(*impl_->mpi, wire,
                  "MPI_Allgather(ideal-gas closure attempt preflight)");
  impl_->active = true;
  impl_->prepared_valid = false;
  impl_->attempt_identity = identity;
  impl_->trial = impl_->committed;
  impl_->latest = IdealGasClosureReport();
  impl_->latest.attempt_identity_ = identity;
  impl_->latest.pressure_mode_ = impl_->committed.mode;
  impl_->latest.configured_pressure_pa_ =
      impl_->spec.configured_thermodynamic_pressure_pa;
  impl_->latest.target_mass_kg_ = impl_->committed.target_mass_kg;
  impl_->latest.collective_count_ = 1U;
  impl_->latest_available = false;
}

const IdealGasClosureReport &
IdealGasClosure::evaluate(FlowState &state, IdealGasClosureStage stage) {
  if (!impl_ || !impl_->active || !state.attempt_active() ||
      stage == IdealGasClosureStage::none)
    throw runtime::Error("ideal-gas closure evaluation is invalid");
  auto &report = impl_->latest;
  report.stage_ = stage;
  ++report.evaluation_count_;
  report.disposition_ = IdealGasClosureDisposition::recoverable_failure;
  report.reason_ = IdealGasClosureFailureReason::none;
  report.lowest_failing_rank_ = -1;
  report.final_metrics_available_ = false;
  impl_->latest_available = true;

  auto &trial = state.solver_layer(FlowLayer::trial);
  const auto &access = state.solver_access_plan();
  const auto extent = trial.interior_extent();
  const std::size_t count = cell_count(extent);
  std::vector<double> rho_tilde(count), q_tilde(count), h(count),
      temperature(count), rho_eos(count), q_eos(count);
  CompensatedSum denominator;
  IdealGasClosureFailureReason local = IdealGasClosureFailureReason::none;
  {
    const auto rho_read = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, impl_->fields.density);
    const auto q_read = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, impl_->spec.enthalpy_density);
    for_each_cell(extent, [&](int i, int j, int k, std::size_t offset) {
      rho_tilde[offset] = rho_read(i, j, k, 0);
      q_tilde[offset] = q_read(i, j, k, 0);
      const auto cell_failure = first_state_failure(
          rho_tilde[offset], q_tilde[offset], impl_->spec.cp_J_per_kg_K);
      local = lower_failure(local, cell_failure);
      if (cell_failure != IdealGasClosureFailureReason::none)
        return;
      h[offset] = q_tilde[offset] / rho_tilde[offset];
      temperature[offset] = h[offset] / impl_->spec.cp_J_per_kg_K;
      denominator.add(impl_->geometry->cell_volume_m3(offset) /
                      temperature[offset]);
    });
  }
  double global_denominator =
      local == IdealGasClosureFailureReason::none ? denominator.value() : 0.0;
  impl_->mpi->allreduce_fp64_in_place(&global_denominator, 1U,
                                      runtime::Fp64ReductionOperation::sum);
  ++report.collective_count_;
  double pressure = impl_->spec.configured_thermodynamic_pressure_pa;
  if (local == IdealGasClosureFailureReason::none) {
    if (!(global_denominator > 0.0) || !std::isfinite(global_denominator)) {
      local = IdealGasClosureFailureReason::denominator_breakdown;
    } else if (impl_->committed.mode == IdealGasPressureMode::closed_dynamic) {
      try {
        pressure = checked_ratio(
            static_cast<long double>(*impl_->committed.target_mass_kg) *
                impl_->spec.gas_constant_J_per_kg_K,
            global_denominator, "ideal-gas candidate pressure is invalid");
      } catch (...) {
        local = IdealGasClosureFailureReason::non_finite_pressure;
      }
    }
  }
  if (local == IdealGasClosureFailureReason::none) {
    if (!std::isfinite(pressure))
      local = IdealGasClosureFailureReason::non_finite_pressure;
    else if (!(pressure > 0.0))
      local = IdealGasClosureFailureReason::non_positive_pressure;
  }
  if (local == IdealGasClosureFailureReason::none) {
    for (std::size_t cell = 0; cell < count; ++cell) {
      IdealGasClosureFailureReason cell_failure =
          IdealGasClosureFailureReason::none;
      try {
        rho_eos[cell] = checked_ratio(
            pressure,
            static_cast<long double>(impl_->spec.gas_constant_J_per_kg_K) *
                temperature[cell],
            "ideal-gas candidate density is invalid");
        const long double product =
            static_cast<long double>(rho_eos[cell]) * h[cell];
        if (!std::isfinite(product) ||
            std::abs(product) >
                static_cast<long double>(std::numeric_limits<double>::max()))
          throw runtime::Error("ideal-gas candidate rho-h is invalid");
        q_eos[cell] = static_cast<double>(product);
      } catch (...) {
        cell_failure = IdealGasClosureFailureReason::non_finite_density;
      }
      if (cell_failure == IdealGasClosureFailureReason::none &&
          !(rho_eos[cell] > 0.0))
        cell_failure = IdealGasClosureFailureReason::non_positive_density;
      local = lower_failure(local, cell_failure);
    }
  }
  std::array<CompensatedSum, 10> local_sums{};
  std::array<double, 8> maxima{};
  std::fill_n(maxima.begin(), 3U, -std::numeric_limits<double>::infinity());
  if (local == IdealGasClosureFailureReason::none) {
    for (std::size_t cell = 0; cell < count; ++cell) {
      const double volume = impl_->geometry->cell_volume_m3(cell);
      const double drho = rho_eos[cell] - rho_tilde[cell];
      const double dq = q_eos[cell] - q_tilde[cell];
      local_sums[0].add(volume * drho * drho);
      local_sums[1].add(volume * rho_eos[cell] * rho_eos[cell]);
      local_sums[2].add(volume * rho_tilde[cell] * rho_tilde[cell]);
      local_sums[3].add(volume * dq * dq);
      local_sums[4].add(volume * q_eos[cell] * q_eos[cell]);
      local_sums[5].add(volume * q_tilde[cell] * q_tilde[cell]);
      local_sums[6].add(volume * drho);
      local_sums[7].add(volume * dq);
      local_sums[8].add(volume * rho_eos[cell]);
      local_sums[9].add(volume * q_eos[cell]);
      maxima[0] = std::max(maxima[0], -h[cell]);
      maxima[1] = std::max(maxima[1], -temperature[cell]);
      maxima[2] = std::max(maxima[2], -rho_eos[cell]);
      maxima[3] = std::max(maxima[3], h[cell]);
      maxima[4] = std::max(maxima[4], temperature[cell]);
      maxima[5] = std::max(maxima[5], rho_eos[cell]);
      maxima[6] = std::max(maxima[6],
                           relative_error(h[cell], impl_->spec.cp_J_per_kg_K *
                                                       temperature[cell]));
      maxima[7] = std::max(
          maxima[7], relative_product_error(rho_eos[cell],
                                            impl_->spec.gas_constant_J_per_kg_K,
                                            temperature[cell], pressure));
    }
    std::array<double, 10> local_values{};
    std::transform(local_sums.begin(), local_sums.end(), local_values.begin(),
                   [](const auto &sum) { return sum.value(); });
    if (!std::all_of(maxima.begin(), maxima.end(),
                     [](double value) { return std::isfinite(value); }))
      local = IdealGasClosureFailureReason::eos_residual;
    else if (!std::all_of(local_values.begin(), local_values.end(),
                          [](double value) { return std::isfinite(value); }))
      local = IdealGasClosureFailureReason::remap_residual;
  }
  const auto selected = agree_failure(*impl_->mpi, local);
  ++report.collective_count_;
  if (selected.reason != IdealGasClosureFailureReason::none) {
    report.reason_ = selected.reason;
    report.lowest_failing_rank_ = selected.rank;
    report.seal();
    return report;
  }
  report.candidate_pressure_available_ = true;
  report.candidate_pressure_pa_ = pressure;
  impl_->trial.thermodynamic_pressure_pa = pressure;

  std::array<double, 10> sums{};
  std::transform(local_sums.begin(), local_sums.end(), sums.begin(),
                 [](const auto &sum) { return sum.value(); });
  impl_->mpi->allreduce_fp64_in_place(sums.data(), sums.size(),
                                      runtime::Fp64ReductionOperation::sum);
  ++report.collective_count_;
  impl_->mpi->allreduce_fp64_in_place(maxima.data(), maxima.size(),
                                      runtime::Fp64ReductionOperation::maximum);
  ++report.collective_count_;
  report.actual_mass_kg_ = sums[8];
  report.enthalpy_min_J_per_kg_ = -maxima[0];
  report.temperature_min_K_ = -maxima[1];
  report.density_min_kg_per_m3_ = -maxima[2];
  report.enthalpy_max_J_per_kg_ = maxima[3];
  report.temperature_max_K_ = maxima[4];
  report.density_max_kg_per_m3_ = maxima[5];
  report.enthalpy_temperature_max_relative_error_ = maxima[6];
  report.eos_max_relative_error_ = maxima[7];
  report.rho_remap_normalized_l2_ =
      std::sqrt(sums[0]) /
      std::max({std::sqrt(sums[1]), std::sqrt(sums[2]), DBL_MIN});
  report.rho_h_remap_normalized_l2_ =
      std::sqrt(sums[3]) /
      std::max({std::sqrt(sums[4]), std::sqrt(sums[5]), DBL_MIN});
  report.rho_remap_relative_conservation_defect_ =
      std::abs(sums[6]) / std::max(std::abs(sums[8] - sums[6]), DBL_MIN);
  report.rho_h_remap_relative_conservation_defect_ =
      std::abs(sums[7]) / std::max(std::abs(sums[9] - sums[7]), DBL_MIN);

  const std::array reported_metrics{
      report.actual_mass_kg_,
      report.enthalpy_min_J_per_kg_,
      report.temperature_min_K_,
      report.density_min_kg_per_m3_,
      report.enthalpy_max_J_per_kg_,
      report.temperature_max_K_,
      report.density_max_kg_per_m3_,
      report.enthalpy_temperature_max_relative_error_,
      report.eos_max_relative_error_,
      report.rho_remap_normalized_l2_,
      report.rho_h_remap_normalized_l2_,
      report.rho_remap_relative_conservation_defect_,
      report.rho_h_remap_relative_conservation_defect_};
  if (!(report.actual_mass_kg_ > 0.0) ||
      !std::all_of(reported_metrics.begin(), reported_metrics.end(),
                   [](double value) { return std::isfinite(value); })) {
    report.reason_ = IdealGasClosureFailureReason::remap_residual;
    report.lowest_failing_rank_ = 0;
    report.seal();
    return report;
  }

  if (stage == IdealGasClosureStage::final) {
    IdealGasClosureFailureReason gate = IdealGasClosureFailureReason::none;
    if (report.enthalpy_temperature_max_relative_error_ > 1.0e-12 ||
        report.eos_max_relative_error_ > 1.0e-12)
      gate = IdealGasClosureFailureReason::eos_residual;
    else if (report.rho_remap_normalized_l2_ > 1.0e-10 ||
             report.rho_h_remap_normalized_l2_ > 1.0e-9)
      gate = IdealGasClosureFailureReason::remap_residual;
    else if (report.rho_remap_relative_conservation_defect_ > 5.0e-11 ||
             (impl_->committed.target_mass_kg &&
              relative_error(report.actual_mass_kg_,
                             *impl_->committed.target_mass_kg) > 5.0e-12))
      gate = IdealGasClosureFailureReason::mass_conservation;
    else if (report.rho_h_remap_relative_conservation_defect_ > 5.0e-11)
      gate = IdealGasClosureFailureReason::enthalpy_conservation;
    if (gate != IdealGasClosureFailureReason::none) {
      report.reason_ = gate;
      report.lowest_failing_rank_ = 0;
      report.seal();
      return report;
    }
  }

  auto rho_write = trial.acquire_write<double>(access, kStatePhase, kStateActor,
                                               impl_->fields.density);
  auto q_write = trial.acquire_write<double>(access, kStatePhase, kStateActor,
                                             impl_->spec.enthalpy_density);
  for_each_cell(extent, [&](int i, int j, int k, std::size_t offset) {
    rho_write(i, j, k, 0) = rho_eos[offset];
    q_write(i, j, k, 0) = q_eos[offset];
  });
  if (stage == IdealGasClosureStage::final) {
    report.final_metrics_available_ = true;
    std::array<double, 2> independent{};
    for (std::size_t cell = 0; cell < count; ++cell) {
      const double stored_h = q_eos[cell] / rho_eos[cell];
      const double stored_t = stored_h / impl_->spec.cp_J_per_kg_K;
      independent[0] = std::max(
          independent[0],
          relative_error(stored_h, impl_->spec.cp_J_per_kg_K * stored_t));
      independent[1] =
          std::max(independent[1],
                   relative_product_error(rho_eos[cell],
                                          impl_->spec.gas_constant_J_per_kg_K,
                                          stored_t, pressure));
    }
    impl_->mpi->allreduce_fp64_in_place(
        independent.data(), independent.size(),
        runtime::Fp64ReductionOperation::maximum);
    ++report.collective_count_;
    report.enthalpy_temperature_max_relative_error_ = independent[0];
    report.eos_max_relative_error_ = independent[1];
    if (independent[0] > 1.0e-12 || independent[1] > 1.0e-12) {
      report.reason_ = IdealGasClosureFailureReason::eos_residual;
      report.lowest_failing_rank_ = 0;
      report.seal();
      return report;
    }
  }
  report.disposition_ = IdealGasClosureDisposition::closed;
  report.reason_ = IdealGasClosureFailureReason::none;
  report.lowest_failing_rank_ = -1;
  report.seal();
  return report;
}

void IdealGasClosure::prepare_commit() {
  if (!impl_ || !impl_->active || !impl_->latest.authenticated() ||
      impl_->latest.disposition() != IdealGasClosureDisposition::closed ||
      impl_->latest.stage() != IdealGasClosureStage::final ||
      !impl_->latest.final_metrics_available())
    throw runtime::Error("ideal-gas closure commit preparation is invalid");
  impl_->prepared = impl_->trial;
  if (impl_->prepared.revision == std::numeric_limits<std::uint64_t>::max())
    throw runtime::Error("ideal-gas closure revision would wrap");
  ++impl_->prepared.revision;
  impl_->prepared_valid = true;
}

void IdealGasClosure::publish_commit() noexcept {
  if (!impl_ || !impl_->active || !impl_->prepared_valid)
    return;
  impl_->committed = impl_->prepared;
  impl_->active = false;
  impl_->prepared_valid = false;
}

void IdealGasClosure::rollback() noexcept {
  if (!impl_)
    return;
  impl_->trial = impl_->committed;
  impl_->prepared = impl_->committed;
  impl_->active = false;
  impl_->prepared_valid = false;
}

const IdealGasClosureReport &IdealGasClosure::latest_report() const {
  if (!impl_ || !impl_->latest_available || !impl_->latest.authenticated())
    throw runtime::Error("ideal-gas closure report is unavailable");
  return impl_->latest;
}

bool IdealGasClosure::matches(const mesh::MeshTopology &topology,
                              const mesh::MeshGeometry &geometry,
                              const boundary::BoundaryRegistry &boundaries,
                              const runtime::MpiContext &mpi,
                              const runtime::FieldRegistry &registry,
                              const FlowFieldIds &fields) const noexcept {
  return impl_ && impl_->topology == &topology &&
         impl_->geometry == &geometry && impl_->boundaries == &boundaries &&
         impl_->mpi == &mpi && impl_->registry == &registry &&
         impl_->fields.density == fields.density &&
         impl_->fields.face_mass_flux == fields.face_mass_flux &&
         impl_->fields.transported_cell_fields ==
             fields.transported_cell_fields;
}

double IdealGasClosure::cp_J_per_kg_K() const noexcept {
  return impl_ ? impl_->spec.cp_J_per_kg_K : 0.0;
}

double IdealGasClosure::gas_constant_J_per_kg_K() const noexcept {
  return impl_ ? impl_->spec.gas_constant_J_per_kg_K : 0.0;
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
bool test::IdealGasClosureTestAccess::
    same_rank_reason_precedence_is_enum_order() noexcept {
  return lower_failure(IdealGasClosureFailureReason::non_positive_density,
                       IdealGasClosureFailureReason::non_finite_enthalpy) ==
             IdealGasClosureFailureReason::non_finite_enthalpy &&
         lower_failure(IdealGasClosureFailureReason::non_positive_enthalpy,
                       IdealGasClosureFailureReason::non_finite_density) ==
             IdealGasClosureFailureReason::non_positive_enthalpy;
}
#endif

} // namespace hundun::flow
