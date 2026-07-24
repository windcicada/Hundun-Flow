// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/ideal_gas_closure.hpp"

#include "checkpoint_v2_detail.hpp"
#include "density_closure_detail.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include "checkpoint_v2_test_access.hpp"
#include "ideal_gas_closure_test_access.hpp"
#endif

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <limits>
#include <new>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::flow {
namespace {

constexpr runtime::PhaseId kStatePhase = 1800U;
constexpr runtime::ActorId kStateActor = 1800U;
constexpr std::uint64_t kClosureReportSeed = 0x696465616c676173ULL;

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
int restore_snapshot_preparation_fault_rank{-1};

bool create_fault(int selected_kind, int selected_rank,
                  test::IdealGasCreateFault fault, int rank) noexcept {
  return selected_kind == static_cast<int>(fault) && selected_rank == rank;
}

void require_create_reduction(test::IdealGasCreateFault fault,
                              const runtime::MpiContext &mpi,
                              const char *operation, int selected_kind,
                              int selected_rank) {
  const auto prepared = runtime::collective_status(
      mpi, !create_fault(selected_kind, selected_rank, fault, mpi.rank()),
      operation);
  if (!prepared.ok)
    runtime::check_mpi_result(MPI_ERR_OTHER, operation);
}
#endif

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

class DensityClosureCreateValidationFailure final : public runtime::Error {
public:
  DensityClosureCreateValidationFailure(IdealGasClosureFailureReason reason,
                                        int rank)
      : runtime::Error("ideal-gas closure initial state is invalid"),
        reason_(reason), rank_(rank) {}

  IdealGasClosureFailureReason reason() const noexcept { return reason_; }
  int failing_rank() const noexcept { return rank_; }

private:
  IdealGasClosureFailureReason reason_;
  int rank_;
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
  std::uint64_t first_owned_global_id{};
  std::uint64_t last_owned_global_id{};
  std::uint64_t local_layout_valid{};
};

static_assert(std::is_trivially_copyable_v<PreflightWire>);

enum class PreflightResultCode : int {
  ok = 0,
  invalid_layout = 1,
  invalid_contract = 2
};

struct PreflightResult final {
  PreflightResultCode code{PreflightResultCode::ok};
  int rank{-1};
};

bool checked_wire_volume(const PreflightWire &candidate,
                         const PreflightWire &common,
                         std::uint64_t &volume) noexcept {
  volume = 1U;
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
      return false;
    volume *= candidate.local_cells[axis];
  }
  return true;
}

bool wires_overlap(const PreflightWire &left,
                   const PreflightWire &right) noexcept {
  bool overlaps = true;
  for (std::size_t axis = 0; axis < 3U; ++axis)
    overlaps = overlaps && left.owned_begin[axis] < right.owned_end[axis] &&
               right.owned_begin[axis] < left.owned_end[axis];
  return overlaps;
}

int agree_preflight(const runtime::MpiContext &mpi, const PreflightWire &local,
                    std::vector<PreflightWire> &preflight_workspace,
                    std::uint64_t &preflight_wire_exchange_count,
                    const char *operation) {
  if (preflight_workspace.size() != static_cast<std::size_t>(mpi.size()))
    throw runtime::Error("ideal-gas closure preflight workspace is invalid");
  static_assert(sizeof(PreflightWire) <=
                static_cast<std::size_t>(std::numeric_limits<int>::max()));
  runtime::check_mpi_result(
      MPI_Allgather(&local, static_cast<int>(sizeof(local)), MPI_BYTE,
                    preflight_workspace.data(),
                    static_cast<int>(sizeof(PreflightWire)), MPI_BYTE,
                    mpi.comm()),
      operation);
  ++preflight_wire_exchange_count;

  PreflightResult result;
  constexpr std::size_t common_bytes = offsetof(PreflightWire, local_cells);
  const auto &common = preflight_workspace.front();
  std::uint64_t global_volume = 1U;
  bool contract_valid =
      common.mode <=
          static_cast<std::uint64_t>(IdealGasPressureMode::open_fixed) &&
      common.global_cells[0] != 0U && common.global_cells[1] != 0U &&
      common.global_cells[2] != 0U;
  for (const auto extent : common.global_cells) {
    if (extent == 0U ||
        global_volume > std::numeric_limits<std::uint64_t>::max() / extent) {
      contract_valid = false;
      break;
    }
    global_volume *= extent;
  }
  std::uint64_t covered = 0U;
  for (int rank = 0; rank < mpi.size(); ++rank) {
    const auto &candidate =
        preflight_workspace[static_cast<std::size_t>(rank)];
    if (candidate.local_layout_valid == 0U && result.rank < 0) {
      result.code = PreflightResultCode::invalid_layout;
      result.rank = rank;
    }
    if (std::memcmp(&candidate, &common, common_bytes) != 0)
      contract_valid = false;
    std::uint64_t volume{};
    if (!checked_wire_volume(candidate, common, volume) ||
        candidate.local_faces == 0U ||
        covered > std::numeric_limits<std::uint64_t>::max() - volume) {
      contract_valid = false;
    } else {
      const auto global_cell_id = [&](std::int64_t i, std::int64_t j,
                                      std::int64_t k) noexcept {
        return ((static_cast<std::uint64_t>(k) * common.global_cells[1] +
                 static_cast<std::uint64_t>(j)) *
                    common.global_cells[0] +
                static_cast<std::uint64_t>(i));
      };
      if (candidate.first_owned_global_id !=
              global_cell_id(candidate.owned_begin[0],
                             candidate.owned_begin[1],
                             candidate.owned_begin[2]) ||
          candidate.last_owned_global_id !=
              global_cell_id(candidate.owned_end[0] - 1,
                             candidate.owned_end[1] - 1,
                             candidate.owned_end[2] - 1))
        contract_valid = false;
      covered += volume;
    }
  }
  if (covered != global_volume)
    contract_valid = false;
  for (int left_rank = 0; left_rank < mpi.size(); ++left_rank)
    for (int right_rank = 0; right_rank < left_rank; ++right_rank)
      if (wires_overlap(
              preflight_workspace[static_cast<std::size_t>(left_rank)],
              preflight_workspace[static_cast<std::size_t>(right_rank)]))
        contract_valid = false;
  if (result.code == PreflightResultCode::ok && !contract_valid)
    result.code = PreflightResultCode::invalid_contract;
  if (result.code == PreflightResultCode::invalid_layout)
    return result.rank;
  if (result.code != PreflightResultCode::ok)
    throw runtime::Error("ideal-gas closure rank preflight disagrees");
  return -1;
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
  auto earlier = [](IdealGasClosureFailureReason left,
                    IdealGasClosureFailureReason right) noexcept {
    if (left == IdealGasClosureFailureReason::none)
      return right;
    if (right == IdealGasClosureFailureReason::none)
      return left;
    return static_cast<std::uint8_t>(left) < static_cast<std::uint8_t>(right)
               ? left
               : right;
  };
  IdealGasClosureFailureReason result = IdealGasClosureFailureReason::none;
  if (!std::isfinite(q))
    result = earlier(result,
                     IdealGasClosureFailureReason::non_finite_enthalpy);
  else if (!(q > 0.0))
    result = earlier(result,
                     IdealGasClosureFailureReason::non_positive_enthalpy);
  if (!std::isfinite(rho))
    result = earlier(result, IdealGasClosureFailureReason::non_finite_density);
  else if (!(rho > 0.0))
    result = earlier(result,
                     IdealGasClosureFailureReason::non_positive_density);
  // Every safely observable candidate enters enum-order selection before any
  // division, so invalid values cannot trigger floating-point exceptions.
  if (result != IdealGasClosureFailureReason::none)
    return result;
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
  return IdealGasClosureFailureReason::none;
}

bool same_layout(runtime::FieldLayoutSet left,
                 runtime::FieldLayoutSet right) noexcept {
  return same(left.cell_interior_extent, right.cell_interior_extent) &&
         left.face_count == right.face_count;
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

double post_store_rank_marker(int ranks, int rank) noexcept {
  return 2.0 * static_cast<double>(ranks - rank);
}

int decode_post_store_rank_marker(int ranks, double marker) noexcept {
  return ranks - static_cast<int>(marker / 2.0);
}

enum class FinalGateOrigin : std::uint8_t {
  none,
  eos,
  rho_remap,
  rho_h_remap,
  mass,
  enthalpy
};

double final_gate_rank_marker(int ranks, FinalGateOrigin origin,
                              int rank) noexcept {
  if (ranks <= 0 || rank < 0 || rank >= ranks ||
      origin == FinalGateOrigin::none)
    return 0.0;
  constexpr int origin_count = 5;
  const auto priority = origin_count + 1 - static_cast<int>(origin);
  return static_cast<double>(static_cast<std::uint64_t>(priority) *
                                 (static_cast<std::uint64_t>(ranks) + 1U) +
                             static_cast<std::uint64_t>(ranks - rank));
}

struct FinalGateSelection final {
  FinalGateOrigin origin{FinalGateOrigin::none};
  int rank{-1};
};

FinalGateSelection decode_final_gate_rank_marker(int ranks,
                                                 double marker) noexcept {
  if (ranks <= 0 || !(marker > 0.0) || !std::isfinite(marker) ||
      marker != std::floor(marker))
    return {};
  const auto encoded = static_cast<std::uint64_t>(marker);
  const auto stride = static_cast<std::uint64_t>(ranks) + 1U;
  const auto priority = encoded / stride;
  const auto remainder = encoded % stride;
  if (priority < 1U || priority > 5U || remainder < 1U ||
      remainder > static_cast<std::uint64_t>(ranks))
    return {};
  const int rank = ranks - static_cast<int>(remainder);
  const auto origin =
      static_cast<FinalGateOrigin>(6U - static_cast<std::uint8_t>(priority));
  if (rank < 0 || rank >= ranks || origin == FinalGateOrigin::none ||
      final_gate_rank_marker(ranks, origin, rank) != marker)
    return {};
  return {origin, rank};
}

IdealGasClosureFailureReason
final_gate_reason(FinalGateOrigin origin) noexcept {
  switch (origin) {
  case FinalGateOrigin::eos:
    return IdealGasClosureFailureReason::eos_residual;
  case FinalGateOrigin::rho_remap:
  case FinalGateOrigin::rho_h_remap:
    return IdealGasClosureFailureReason::remap_residual;
  case FinalGateOrigin::mass:
    return IdealGasClosureFailureReason::mass_conservation;
  case FinalGateOrigin::enthalpy:
    return IdealGasClosureFailureReason::enthalpy_conservation;
  case FinalGateOrigin::none:
    return IdealGasClosureFailureReason::none;
  }
  return IdealGasClosureFailureReason::none;
}

} // namespace

bool detail::validate_ideal_gas_restore_state(
    const runtime::MpiContext &mpi, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    double cp, double gas_constant, double configured_pressure,
    const FlowLayerValues &history,
    const FlowLayerValues &committed, const IdealGasClosureState &restored,
    std::uint64_t &collective_count) {
  if (!(cp > 0.0) || !(gas_constant > 0.0) ||
      !(configured_pressure > 0.0) || !std::isfinite(cp) ||
      !std::isfinite(gas_constant) || !std::isfinite(configured_pressure))
    return false;
  const bool open = boundaries.open_domain();
  bool valid =
      restored.mode == (open ? IdealGasPressureMode::open_fixed
                             : IdealGasPressureMode::closed_dynamic) &&
      restored.thermodynamic_pressure_pa > 0.0 &&
      std::isfinite(restored.thermodynamic_pressure_pa) &&
      restored.revision != std::numeric_limits<std::uint64_t>::max() &&
      restored.target_mass_kg.has_value() == !open;
  if (open)
    valid =
        valid &&
        fp_bits(restored.thermodynamic_pressure_pa) ==
            fp_bits(configured_pressure);
  const double target_mass = restored.target_mass_kg.value_or(
      std::numeric_limits<double>::quiet_NaN());
  if (!open)
    valid = valid && target_mass > 0.0 && std::isfinite(target_mass);

  std::array<double, 3> sums{};
  bool shapes_valid = true;
  const auto accumulate = [&](const FlowLayerValues &layer,
                              std::size_t mass_index,
                              bool inverse_temperature) {
    if (layer.transported_cell_fields.empty() ||
        layer.density.size() != topology.owned_cell_count() ||
        layer.transported_cell_fields.front().size() !=
            topology.owned_cell_count()) {
      valid = false;
      shapes_valid = false;
      return;
    }
    for (std::size_t cell = 0; cell < topology.owned_cell_count(); ++cell) {
      const double rho = layer.density[cell];
      const double q = layer.transported_cell_fields.front()[cell];
      const double enthalpy = q / rho;
      const double temperature = enthalpy / cp;
      const double volume = geometry.cell_volume_m3(cell);
      if (!(rho > 0.0) || !(q > 0.0) || !(enthalpy > 0.0) ||
          !(temperature > 0.0) || !(volume > 0.0) ||
          !std::isfinite(rho) || !std::isfinite(q) ||
          !std::isfinite(enthalpy) || !std::isfinite(temperature) ||
          !std::isfinite(volume)) {
        valid = false;
        continue;
      }
      sums[mass_index] += volume * rho;
      if (inverse_temperature)
        sums[2] += volume / temperature;
    }
  };
  accumulate(history, 0U, !open);
  accumulate(committed, 1U, false);
  if (!open) {
    mpi.allreduce_fp64_in_place(sums.data(), sums.size(),
                                runtime::Fp64ReductionOperation::sum);
    ++collective_count;
    valid = valid && sums[2] > 0.0 && std::isfinite(sums[2]) &&
            relative_error(sums[0], target_mass) <= 5.0e-12 &&
            relative_error(sums[1], target_mass) <= 5.0e-12;
  }
  if (!shapes_valid)
    return false;
  const double history_pressure =
      open ? restored.thermodynamic_pressure_pa
           : target_mass * gas_constant / sums[2];
  const auto eos_valid = [&](const FlowLayerValues &layer, double pressure) {
    if (!(pressure > 0.0) || !std::isfinite(pressure))
      return false;
    for (std::size_t cell = 0; cell < topology.owned_cell_count(); ++cell) {
      const double rho = layer.density[cell];
      const double q = layer.transported_cell_fields.front()[cell];
      const double temperature = (q / rho) / cp;
      const double ratio = rho * gas_constant * temperature / pressure;
      if (!std::isfinite(ratio) || std::abs(ratio - 1.0) > 1.0e-12)
        return false;
    }
    return true;
  };
  return valid && eos_valid(history, history_pressure) &&
         eos_valid(committed, restored.thermodynamic_pressure_pa);
}

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
  std::array<std::vector<double>, 6> cell_workspace;
  std::vector<PreflightWire> preflight_workspace;
  std::uint64_t preflight_wire_exchange_count{};
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  int facade_create_fault_rank{-1};
  int material_factory_create_fault_rank{-1};
  int controlled_allocation_rank{-1};
  bool allocation_observation_active{};
  int attempt_preparation_fault_kind{-1};
  int attempt_preparation_fault_rank{-1};
  int post_store_fault_rank{-1};
  bool post_store_fault_enthalpy{};
  int candidate_precedence_fault_rank{-1};
  IdealGasClosureStage stage_failure_stage{IdealGasClosureStage::none};
  IdealGasClosureFailureReason stage_failure_reason{
      IdealGasClosureFailureReason::none};
  int stage_failure_rank{-1};
  int metric_gate_failure_kind{-1};
  int metric_gate_failure_rank{-1};
  int post_assessment_fault_kind{-1};
  int post_assessment_fault_rank{-1};
  std::uint8_t outer_failure_point{std::numeric_limits<std::uint8_t>::max()};
  int outer_failure_rank{-1};
  int state_prepare_fault_rank{-1};
  int closure_prepare_fault_rank{-1};
  int post_store_mpi_fault_rank{-1};
  int attempt_layout_fault_rank{-1};
  bool outlet_backflow_fault{};
  std::vector<test::IdealGasHaloTraceEntry> halo_trace;
#endif
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
  return create_internal(topology, geometry, boundaries, mpi, registry, fields,
                         state, spec, nullptr
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
                         ,
                         -1, -1
#endif
  );
}

IdealGasClosure IdealGasClosure::restore(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::MpiContext &mpi, const runtime::FieldRegistry &registry,
    const FlowFieldIds &fields, const FlowState &state,
    IdealGasClosureSpec spec, IdealGasClosureState restored) {
  FlowLayerValues history;
  FlowLayerValues committed;
  bool snapshots_prepared = true;
  std::array<char, 96> snapshot_message{};
  std::string_view local_message;
  try {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (restore_snapshot_preparation_fault_rank == mpi.rank()) {
      restore_snapshot_preparation_fault_rank = -1;
      throw std::bad_alloc();
    }
#endif
    history = state.snapshot(FlowLayer::history);
    committed = state.snapshot(FlowLayer::committed);
  } catch (...) {
    snapshots_prepared = false;
    const int written = std::snprintf(
        snapshot_message.data(), snapshot_message.size(),
        "ideal-gas closure restore snapshot preparation failed on rank %d",
        mpi.rank());
    local_message =
        written > 0 && static_cast<std::size_t>(written) <
                           snapshot_message.size()
            ? std::string_view(snapshot_message.data(),
                               static_cast<std::size_t>(written))
            : std::string_view(
                  "ideal-gas closure restore snapshot preparation failed");
  }
  const auto snapshot_status =
      runtime::collective_status(mpi, snapshots_prepared, local_message);
  if (!snapshot_status.ok)
    throw runtime::Error(snapshot_status.message);
  std::uint64_t validation_collectives{};
  const bool valid = detail::validate_ideal_gas_restore_state(
      mpi, topology, geometry, boundaries, spec.cp_J_per_kg_K,
      spec.gas_constant_J_per_kg_K,
      spec.configured_thermodynamic_pressure_pa,
      history, committed, restored, validation_collectives);
  const auto status = runtime::collective_status(
      mpi, valid, "ideal-gas closure restore validation");
  if (!status.ok)
    throw runtime::Error(status.message);
  return create_internal(topology, geometry, boundaries, mpi, registry, fields,
                         state, spec, &restored
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
                         ,
                         -1, -1
#endif
  );
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void test::set_ideal_gas_restore_snapshot_preparation_fault(
    int rank) noexcept {
  restore_snapshot_preparation_fault_rank = rank;
}
#endif

IdealGasClosure IdealGasClosure::create_internal(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::MpiContext &mpi, const runtime::FieldRegistry &registry,
    const FlowFieldIds &fields, const FlowState &state, IdealGasClosureSpec spec,
    const IdealGasClosureState *restored_authority
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    ,
    int create_fault_kind, int create_fault_rank
#endif
) {
  IdealGasPressureMode mode{IdealGasPressureMode::closed_dynamic};
  runtime::FieldLayoutSet committed_layout{};
  runtime::FieldLayoutSet history_layout{};
  runtime::FieldLayoutSet trial_layout{};
  runtime::Int3 local_extent{};
  PreflightWire wire{};
  wire.enthalpy = spec.enthalpy_density;
  wire.density = fields.density;
  wire.face_flux = fields.face_mass_flux;
  wire.cp = fp_bits(spec.cp_J_per_kg_K);
  wire.gas_constant = fp_bits(spec.gas_constant_J_per_kg_K);
  wire.configured_pressure = fp_bits(spec.configured_thermodynamic_pressure_pa);
  bool preparation_valid = true;
  try {
    geometry.require_compatible(topology);
    if (!state.impl_)
      throw runtime::Error("ideal-gas closure state has been moved from");
    mode = boundaries.open_domain() ? IdealGasPressureMode::open_fixed
                                    : IdealGasPressureMode::closed_dynamic;
    committed_layout = state.layer(FlowLayer::committed).layout_set();
    history_layout = state.layer(FlowLayer::history).layout_set();
    trial_layout = state.layer(FlowLayer::trial).layout_set();
    local_extent = committed_layout.cell_interior_extent;
    const auto global = topology.global_extent();
    const auto owned = topology.owned_global_box();
    wire.mode = static_cast<std::uint64_t>(mode);
    wire.global_cells[0] = static_cast<std::uint64_t>(global.x);
    wire.global_cells[1] = static_cast<std::uint64_t>(global.y);
    wire.global_cells[2] = static_cast<std::uint64_t>(global.z);
    wire.local_cells[0] = static_cast<std::uint64_t>(local_extent.x);
    wire.local_cells[1] = static_cast<std::uint64_t>(local_extent.y);
    wire.local_cells[2] = static_cast<std::uint64_t>(local_extent.z);
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const auto begin =
          std::array{owned.begin.x, owned.begin.y, owned.begin.z};
      const auto end = std::array{owned.end.x, owned.end.y, owned.end.z};
      wire.owned_begin[axis] = begin[axis];
      wire.owned_end[axis] = end[axis];
    }
    wire.local_faces = committed_layout.face_count;
    wire.first_owned_global_id = topology.global_cell_id(owned.begin);
    wire.last_owned_global_id = topology.global_cell_id(
        {owned.end.x - 1, owned.end.y - 1, owned.end.z - 1});
    wire.local_layout_valid =
        same_layout(committed_layout, history_layout) &&
                same_layout(committed_layout, trial_layout) &&
                committed_layout.face_count == topology.local_face_count()
            ? 1U
            : 0U;
  } catch (...) {
    preparation_valid = false;
    wire.local_layout_valid = 0U;
  }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (create_fault(create_fault_kind, create_fault_rank,
                   test::IdealGasCreateFault::local_preparation, mpi.rank())) {
    preparation_valid = false;
    wire.local_layout_valid = 0U;
  } else if (create_fault(create_fault_kind, create_fault_rank,
                          test::IdealGasCreateFault::mode_disagreement,
                          mpi.rank())) {
    wire.mode =
        static_cast<std::uint64_t>(IdealGasPressureMode::open_fixed) + 1U;
  } else if (create_fault(create_fault_kind, create_fault_rank,
                          test::IdealGasCreateFault::ownership_gap,
                          mpi.rank())) {
    ++wire.owned_begin[0];
    --wire.local_cells[0];
  } else if (create_fault(create_fault_kind, create_fault_rank,
                          test::IdealGasCreateFault::ownership_overlap,
                          mpi.rank())) {
    if (wire.owned_begin[0] > 0) {
      --wire.owned_begin[0];
      ++wire.local_cells[0];
    } else {
      ++wire.owned_end[0];
      ++wire.local_cells[0];
    }
  } else if (create_fault(create_fault_kind, create_fault_rank,
                          test::IdealGasCreateFault::ownership_swap,
                          mpi.rank())) {
    std::swap(wire.first_owned_global_id, wire.last_owned_global_id);
  }
#endif
  std::vector<PreflightWire> preflight_workspace;
  bool preflight_workspace_ready = true;
  try {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (create_fault(
            create_fault_kind, create_fault_rank,
            test::IdealGasCreateFault::preflight_workspace_allocation,
            mpi.rank()))
      throw std::bad_alloc();
#endif
    preflight_workspace.resize(static_cast<std::size_t>(mpi.size()));
  } catch (...) {
    preflight_workspace_ready = false;
  }
  const int local_workspace_failure =
      preflight_workspace_ready ? mpi.size() : mpi.rank();
  int workspace_failure = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&local_workspace_failure, &workspace_failure, 1, MPI_INT,
                    MPI_MIN, mpi.comm()),
      "MPI_Allreduce(ideal-gas closure preflight workspace allocation)");
  if (workspace_failure != mpi.size())
    throw detail::DensityClosurePreflightFailure(workspace_failure);

  std::uint64_t preflight_wire_exchange_count = 0U;
  const int preparation_failure = agree_preflight(
      mpi, wire, preflight_workspace, preflight_wire_exchange_count,
      "MPI_Allgather(ideal-gas closure create preflight)");
  if (preparation_failure >= 0)
    throw detail::DensityClosurePreflightFailure(preparation_failure);

  bool valid = preparation_valid;
  IdealGasClosureFailureReason local_validation_failure =
      IdealGasClosureFailureReason::none;
  FlowLayerValues committed;
  FlowLayerValues history;
  std::array<CompensatedSum, 4> local_sums{};
  std::size_t owned_cell_count{};
  std::unique_ptr<Impl> impl;
  try {
    valid = valid && !state.attempt_active() &&
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
            (restored_authority == nullptr ||
             restored_authority->mode == mode);
    const auto owned = topology.owned_global_box();
    valid = valid && same(local_extent, {owned.end.x - owned.begin.x,
                                         owned.end.y - owned.begin.y,
                                         owned.end.z - owned.begin.z});
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
    committed = state.snapshot(FlowLayer::committed);
    history = state.snapshot(FlowLayer::history);
    constexpr std::size_t enthalpy_index = 0U;
    owned_cell_count = topology.owned_cell_count();
    valid = valid && committed.density.size() == owned_cell_count &&
            history.density.size() == owned_cell_count &&
            committed.transported_cell_fields.size() > enthalpy_index &&
            history.transported_cell_fields.size() > enthalpy_index &&
            committed.transported_cell_fields[enthalpy_index].size() ==
                owned_cell_count &&
            history.transported_cell_fields[enthalpy_index].size() ==
                owned_cell_count;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (create_fault(create_fault_kind, create_fault_rank,
                     test::IdealGasCreateFault::construction_allocation,
                     mpi.rank()))
      throw std::bad_alloc();
#endif
    impl = std::make_unique<Impl>();
    impl->topology = &topology;
    impl->geometry = &geometry;
    impl->boundaries = &boundaries;
    impl->mpi = &mpi;
    impl->registry = &registry;
    impl->fields = fields;
    impl->spec = spec;
    impl->preflight_workspace = std::move(preflight_workspace);
    impl->preflight_wire_exchange_count = preflight_wire_exchange_count;
    for (auto &workspace : impl->cell_workspace)
      workspace.assign(owned_cell_count, 0.0);
    if (valid) {
      for (std::size_t cell = 0; cell < owned_cell_count; ++cell) {
        const double volume = geometry.cell_volume_m3(cell);
        const double rho = committed.density[cell];
        const double rho_history = history.density[cell];
        const double q = committed.transported_cell_fields[0][cell];
        const double q_history = history.transported_cell_fields[0][cell];
        const auto current_failure =
            first_state_failure(rho, q, spec.cp_J_per_kg_K);
        const auto history_failure =
            first_state_failure(rho_history, q_history, spec.cp_J_per_kg_K);
        if (current_failure != IdealGasClosureFailureReason::none ||
            history_failure != IdealGasClosureFailureReason::none) {
          local_validation_failure =
              lower_failure(local_validation_failure,
                            lower_failure(current_failure, history_failure));
          valid = false;
          continue;
        }
        if (!(volume > 0.0) || !std::isfinite(volume)) {
          local_validation_failure = lower_failure(
              local_validation_failure,
              IdealGasClosureFailureReason::invalid_input);
          valid = false;
          continue;
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
        valid =
            inlet_state.has_value() && inlet_state->temperature_K.has_value();
        double h = valid ? inlet_state->enthalpy_J_per_kg : 0.0;
        double rho = valid ? inlet_state->density_kg_per_m3 : 0.0;
        double temperature = valid ? *inlet_state->temperature_K : 0.0;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
        if (create_fault(create_fault_kind, create_fault_rank,
                         test::IdealGasCreateFault::inlet_cp, mpi.rank()))
          h *= 1.01;
        else if (create_fault(create_fault_kind, create_fault_rank,
                              test::IdealGasCreateFault::inlet_gas_constant,
                              mpi.rank()))
          rho *= 1.01;
        else if (create_fault(create_fault_kind, create_fault_rank,
                              test::IdealGasCreateFault::inlet_pressure,
                              mpi.rank())) {
          h *= 1.01;
          temperature *= 1.01;
        }
#endif
        valid =
            valid && h > 0.0 && rho > 0.0 && temperature > 0.0 &&
            std::isfinite(h) && std::isfinite(rho) &&
            std::isfinite(temperature) &&
            relative_error(h, spec.cp_J_per_kg_K * temperature) <= 1.0e-12 &&
            relative_product_error(
                rho, spec.gas_constant_J_per_kg_K, temperature,
                spec.configured_thermodynamic_pressure_pa) <= 1.0e-12;
      }
    }
  } catch (...) {
    valid = false;
    local_validation_failure = IdealGasClosureFailureReason::invalid_input;
  }
  if (!valid &&
      local_validation_failure == IdealGasClosureFailureReason::none)
    local_validation_failure = IdealGasClosureFailureReason::invalid_input;
  const auto validation = agree_failure(mpi, local_validation_failure);
  if (validation.reason != IdealGasClosureFailureReason::none)
    throw DensityClosureCreateValidationFailure(validation.reason,
                                                validation.rank);

  std::array<double, 4> sums{};
  std::transform(local_sums.begin(), local_sums.end(), sums.begin(),
                 [](const auto &sum) { return sum.value(); });
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  require_create_reduction(test::IdealGasCreateFault::sum_reduction, mpi,
                           "ideal-gas create sum reduction", create_fault_kind,
                           create_fault_rank);
#endif
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
  for (std::size_t cell = 0; cell < owned_cell_count; ++cell) {
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
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  require_create_reduction(test::IdealGasCreateFault::maximum_reduction, mpi,
                           "ideal-gas create maximum reduction",
                           create_fault_kind, create_fault_rank);
#endif
  mpi.allreduce_fp64_in_place(maxima.data(), maxima.size(),
                              runtime::Fp64ReductionOperation::maximum);
  const bool consistent =
      sums[1] > 0.0 && sums[3] > 0.0 &&
      (mode == IdealGasPressureMode::open_fixed ||
       (target_mass > 0.0 && std::isfinite(target_mass) &&
        relative_error(sums[2], target_mass) <= 5.0e-12)) &&
      (restored_authority != nullptr ||
       relative_error(current_pressure,
                      spec.configured_thermodynamic_pressure_pa) <= 1.0e-12) &&
      maxima[6] <= 1.0e-12 && maxima[7] <= 1.0e-12 && maxima[8] <= 1.0e-12 &&
      maxima[9] <= 1.0e-12;
  if (!consistent)
    throw runtime::Error("ideal-gas closure initial EOS is inconsistent");

  impl->committed =
      restored_authority != nullptr
          ? *restored_authority
          : IdealGasClosureState{
                mode, current_pressure,
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
  runtime::FieldLayoutSet committed_layout{};
  runtime::FieldLayoutSet history_layout{};
  runtime::FieldLayoutSet trial_layout{};
  PreflightWire wire{};
  wire.mode = static_cast<std::uint64_t>(impl_->committed.mode);
  wire.enthalpy = impl_->spec.enthalpy_density;
  wire.density = impl_->fields.density;
  wire.face_flux = impl_->fields.face_mass_flux;
  wire.cp = fp_bits(impl_->spec.cp_J_per_kg_K);
  wire.gas_constant = fp_bits(impl_->spec.gas_constant_J_per_kg_K);
  wire.configured_pressure =
      fp_bits(impl_->spec.configured_thermodynamic_pressure_pa);
  wire.attempt_identity = identity;
  bool layout_valid = true;
  try {
    if (!state.impl_)
      throw runtime::Error("ideal-gas closure state has been moved from");
    committed_layout = state.layer(FlowLayer::committed).layout_set();
    history_layout = state.layer(FlowLayer::history).layout_set();
    trial_layout = state.layer(FlowLayer::trial).layout_set();
    const auto local_extent = committed_layout.cell_interior_extent;
    const auto global = impl_->topology->global_extent();
    const auto owned = impl_->topology->owned_global_box();
    wire.global_cells[0] = static_cast<std::uint64_t>(global.x);
    wire.global_cells[1] = static_cast<std::uint64_t>(global.y);
    wire.global_cells[2] = static_cast<std::uint64_t>(global.z);
    wire.local_cells[0] = static_cast<std::uint64_t>(local_extent.x);
    wire.local_cells[1] = static_cast<std::uint64_t>(local_extent.y);
    wire.local_cells[2] = static_cast<std::uint64_t>(local_extent.z);
    const auto begin = std::array{owned.begin.x, owned.begin.y, owned.begin.z};
    const auto end = std::array{owned.end.x, owned.end.y, owned.end.z};
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      wire.owned_begin[axis] = begin[axis];
      wire.owned_end[axis] = end[axis];
    }
    wire.local_faces = committed_layout.face_count;
    wire.first_owned_global_id = impl_->topology->global_cell_id(owned.begin);
    wire.last_owned_global_id = impl_->topology->global_cell_id(
        {owned.end.x - 1, owned.end.y - 1, owned.end.z - 1});
    layout_valid =
        !impl_->active && state.attempt_active() && identity != 0U &&
        same_layout(committed_layout, history_layout) &&
        same_layout(committed_layout, trial_layout) &&
        committed_layout.face_count == impl_->topology->local_face_count();
  } catch (...) {
    layout_valid = false;
  }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (impl_->attempt_layout_fault_rank >= 0) {
    const int target = impl_->attempt_layout_fault_rank;
    impl_->attempt_layout_fault_rank = -1;
    if (target == impl_->mpi->rank())
      layout_valid = false;
  }
#endif
  wire.local_layout_valid = layout_valid ? 1U : 0U;
  const int layout_failure = agree_preflight(
      *impl_->mpi, wire, impl_->preflight_workspace,
      impl_->preflight_wire_exchange_count,
      "MPI_Allgather(ideal-gas closure attempt preflight)");
  if (layout_failure >= 0)
    throw detail::DensityClosurePreflightFailure(layout_failure);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  impl_->halo_trace.clear();
#endif
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
  if (std::any_of(impl_->cell_workspace.begin(), impl_->cell_workspace.end(),
                  [count](const auto &workspace) {
                    return workspace.size() != count;
                  }))
    throw runtime::Error("ideal-gas closure workspace layout is invalid");
  auto &rho_tilde = impl_->cell_workspace[0];
  auto &q_tilde = impl_->cell_workspace[1];
  auto &h = impl_->cell_workspace[2];
  auto &temperature = impl_->cell_workspace[3];
  auto &rho_eos = impl_->cell_workspace[4];
  auto &q_eos = impl_->cell_workspace[5];
  std::optional<runtime::FieldView<double>> rho_write;
  std::optional<runtime::FieldView<double>> q_write;
  CompensatedSum denominator;
  IdealGasClosureFailureReason local = IdealGasClosureFailureReason::none;
  try {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    const auto local_fault =
        stage == IdealGasClosureStage::predictor
            ? test::IdealGasAttemptPreparationFault::predictor_local
        : stage == IdealGasClosureStage::provisional
            ? test::IdealGasAttemptPreparationFault::provisional_local
            : test::IdealGasAttemptPreparationFault::final_local;
    if (impl_->attempt_preparation_fault_kind ==
            static_cast<int>(local_fault) &&
        impl_->attempt_preparation_fault_rank == impl_->mpi->rank()) {
      impl_->attempt_preparation_fault_kind = -1;
      impl_->attempt_preparation_fault_rank = -1;
      throw runtime::Error("injected ideal-gas local preparation failure");
    }
    if (stage == IdealGasClosureStage::predictor &&
        impl_->controlled_allocation_rank == impl_->mpi->rank()) {
      impl_->controlled_allocation_rank = -1;
      auto controlled = std::make_unique<std::byte[]>(1U);
      static_cast<void>(controlled);
      throw runtime::Error("injected ideal-gas controlled allocation");
    }
    if (stage == IdealGasClosureStage::predictor &&
        impl_->candidate_precedence_fault_rank >= 0) {
      const int target = impl_->candidate_precedence_fault_rank;
      impl_->candidate_precedence_fault_rank = -1;
      if (target == impl_->mpi->rank()) {
        auto rho_fault = trial.acquire_write<double>(
            access, kStatePhase, kStateActor, impl_->fields.density);
        auto q_fault = trial.acquire_write<double>(
            access, kStatePhase, kStateActor, impl_->spec.enthalpy_density);
        rho_fault(0, 0, 0, 0) = 0.0;
        q_fault(0, 0, 0, 0) = std::numeric_limits<double>::quiet_NaN();
      }
    }
#endif
    const auto rho_read = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, impl_->fields.density);
    const auto q_read = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, impl_->spec.enthalpy_density);
    rho_write.emplace(trial.acquire_write<double>(
        access, kStatePhase, kStateActor, impl_->fields.density));
    q_write.emplace(trial.acquire_write<double>(
        access, kStatePhase, kStateActor, impl_->spec.enthalpy_density));
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    const auto write_fault =
        stage == IdealGasClosureStage::predictor
            ? test::IdealGasAttemptPreparationFault::
                  predictor_write_capability
        : stage == IdealGasClosureStage::provisional
            ? test::IdealGasAttemptPreparationFault::
                  provisional_write_capability
            : test::IdealGasAttemptPreparationFault::final_write_capability;
    if (impl_->attempt_preparation_fault_kind ==
            static_cast<int>(write_fault) &&
        impl_->attempt_preparation_fault_rank == impl_->mpi->rank()) {
      impl_->attempt_preparation_fault_kind = -1;
      impl_->attempt_preparation_fault_rank = -1;
      throw runtime::Error("injected ideal-gas write capability failure");
    }
#endif
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
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    local = IdealGasClosureFailureReason::invalid_input;
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
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (stage == impl_->stage_failure_stage) {
    if (impl_->stage_failure_rank == impl_->mpi->rank())
      local = lower_failure(local, impl_->stage_failure_reason);
    impl_->stage_failure_stage = IdealGasClosureStage::none;
    impl_->stage_failure_reason = IdealGasClosureFailureReason::none;
    impl_->stage_failure_rank = -1;
  }
#endif
  std::array<CompensatedSum, 10> local_sums{};
  std::array<double, 10> local_values{};
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
  }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (stage == IdealGasClosureStage::final &&
      impl_->metric_gate_failure_rank >= 0) {
    const int target = impl_->metric_gate_failure_rank;
    const int kind = impl_->metric_gate_failure_kind;
    impl_->metric_gate_failure_kind = -1;
    impl_->metric_gate_failure_rank = -1;
    if (target == impl_->mpi->rank()) {
      switch (static_cast<test::IdealGasMetricGateFault>(kind)) {
      case test::IdealGasMetricGateFault::eos:
        maxima[7] = std::max(maxima[7], 1.0e-6);
        break;
      case test::IdealGasMetricGateFault::rho_remap:
        local_sums[0].add(1.0e6);
        break;
      case test::IdealGasMetricGateFault::rho_h_remap:
        local_sums[3].add(1.0e12);
        break;
      case test::IdealGasMetricGateFault::mass:
        local_sums[6].add(1.0);
        break;
      case test::IdealGasMetricGateFault::enthalpy:
        local_sums[7].add(1.0e6);
        break;
      }
    }
  }
#endif
  std::transform(local_sums.begin(), local_sums.end(), local_values.begin(),
                 [](const auto &sum) { return sum.value(); });
  if (local == IdealGasClosureFailureReason::none) {
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

  auto sums = local_values;
  impl_->mpi->allreduce_fp64_in_place(sums.data(), sums.size(),
                                      runtime::Fp64ReductionOperation::sum);
  ++report.collective_count_;
  if (stage == IdealGasClosureStage::final) {
    const double rho_remap =
        std::sqrt(sums[0]) /
        std::max({std::sqrt(sums[1]), std::sqrt(sums[2]), DBL_MIN});
    const double rho_h_remap =
        std::sqrt(sums[3]) /
        std::max({std::sqrt(sums[4]), std::sqrt(sums[5]), DBL_MIN});
    const double mass_defect =
        std::abs(sums[6]) / std::max(std::abs(sums[8] - sums[6]), DBL_MIN);
    const double enthalpy_defect =
        std::abs(sums[7]) / std::max(std::abs(sums[9] - sums[7]), DBL_MIN);
    const bool target_mass_failure =
        impl_->committed.target_mass_kg &&
        relative_error(sums[8], *impl_->committed.target_mass_kg) > 5.0e-12;
    const double local_rho_remap =
        std::sqrt(local_values[0]) /
        std::max(
            {std::sqrt(local_values[1]), std::sqrt(local_values[2]), DBL_MIN});
    const double local_rho_h_remap =
        std::sqrt(local_values[3]) /
        std::max(
            {std::sqrt(local_values[4]), std::sqrt(local_values[5]), DBL_MIN});
    const double local_mass_defect =
        std::abs(local_values[6]) /
        std::max(std::abs(local_values[8] - local_values[6]), DBL_MIN);
    const double local_enthalpy_defect =
        std::abs(local_values[7]) /
        std::max(std::abs(local_values[9] - local_values[7]), DBL_MIN);
    FinalGateOrigin origin = FinalGateOrigin::none;
    if (maxima[6] > 1.0e-12 || maxima[7] > 1.0e-12)
      origin = FinalGateOrigin::eos;
    else if (rho_remap > 1.0e-10 && local_rho_remap > 1.0e-10)
      origin = FinalGateOrigin::rho_remap;
    else if (rho_h_remap > 1.0e-9 && local_rho_h_remap > 1.0e-9)
      origin = FinalGateOrigin::rho_h_remap;
    else if ((mass_defect > 5.0e-11 || target_mass_failure) &&
             (local_mass_defect > 5.0e-11 ||
              (target_mass_failure && impl_->mpi->rank() == 0)))
      origin = FinalGateOrigin::mass;
    else if (enthalpy_defect > 5.0e-11 && local_enthalpy_defect > 5.0e-11)
      origin = FinalGateOrigin::enthalpy;
    if (origin != FinalGateOrigin::none)
      maxima[0] = final_gate_rank_marker(impl_->mpi->size(), origin,
                                         impl_->mpi->rank());
  }
  impl_->mpi->allreduce_fp64_in_place(maxima.data(), maxima.size(),
                                      runtime::Fp64ReductionOperation::maximum);
  ++report.collective_count_;
  const auto final_gate_selection =
      stage == IdealGasClosureStage::final
          ? decode_final_gate_rank_marker(impl_->mpi->size(), maxima[0])
          : FinalGateSelection{};
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
      report.lowest_failing_rank_ =
          final_gate_reason(final_gate_selection.origin) == gate
              ? final_gate_selection.rank
              : 0;
      report.seal();
      return report;
    }
  }

  // Both checked write capabilities and the exact field/count/layout contract
  // were validated before the ordinary agreement. Publication after the
  // metric gate is deliberately a private no-allocation/no-fail operation;
  // no checked or kernel view is retained across a collective.
  static_assert(noexcept(trial.publish_validated_cell_interior_double(
      impl_->fields.density, rho_eos.data(), rho_eos.size())));
  trial.publish_validated_cell_interior_double(
      impl_->fields.density, rho_eos.data(), rho_eos.size());
  trial.publish_validated_cell_interior_double(
      impl_->spec.enthalpy_density, q_eos.data(), q_eos.size());

  if (stage == IdealGasClosureStage::final) {
    report.final_metrics_available_ = true;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (impl_->post_store_fault_rank >= 0) {
      const int target = impl_->post_store_fault_rank;
      impl_->post_store_fault_rank = -1;
      const auto field = impl_->post_store_fault_enthalpy
                             ? impl_->spec.enthalpy_density
                             : impl_->fields.density;
      if (target == impl_->mpi->rank()) {
        auto corrupt = trial.acquire_write<double>(access, kStatePhase,
                                                   kStateActor, field);
        if (impl_->post_store_fault_enthalpy)
          corrupt(0, 0, 0, 0) *= 1.01;
        else
          corrupt(0, 0, 0, 0) = 0.0;
      }
    }
#endif
    std::array<double, 2> independent{};
    try {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      if (impl_->attempt_preparation_fault_kind ==
              static_cast<int>(
                  test::IdealGasAttemptPreparationFault::final_readback) &&
          impl_->attempt_preparation_fault_rank == impl_->mpi->rank()) {
        impl_->attempt_preparation_fault_kind = -1;
        impl_->attempt_preparation_fault_rank = -1;
        throw runtime::Error("injected ideal-gas final readback failure");
      }
#endif
      const auto stored_rho = trial.acquire_read<double>(
          access, kStatePhase, kStateActor, impl_->fields.density);
      const auto stored_q = trial.acquire_read<double>(
          access, kStatePhase, kStateActor, impl_->spec.enthalpy_density);
      for_each_cell(extent, [&](int i, int j, int k, std::size_t) {
        const double rho = stored_rho(i, j, k, 0);
        const double q = stored_q(i, j, k, 0);
        const auto failure = first_state_failure(
            rho, q, impl_->spec.cp_J_per_kg_K);
        if (failure != IdealGasClosureFailureReason::none) {
          independent[0] = std::max(independent[0], 1.0);
          independent[1] = std::max(independent[1], 1.0);
          return;
        }
        const double stored_h = q / rho;
        const double stored_t = stored_h / impl_->spec.cp_J_per_kg_K;
        independent[0] = std::max(
            independent[0],
            relative_error(stored_h, impl_->spec.cp_J_per_kg_K * stored_t));
        independent[1] = std::max(
            independent[1], relative_product_error(
                                rho, impl_->spec.gas_constant_J_per_kg_K,
                                stored_t, pressure));
      });
    } catch (const runtime::MpiOperationError &) {
      throw;
    } catch (...) {
      independent[0] =
          post_store_rank_marker(impl_->mpi->size(), impl_->mpi->rank());
      independent[1] = 1.0;
    }
    const bool local_post_store_failure =
        !(independent[0] <= 1.0e-12 && independent[1] <= 1.0e-12);
    // On failure the h-T value is unavailable.  Its slot carries an exact
    // integer rank marker so MAX chooses the lowest failing rank without
    // changing the frozen two-scalar FP64 payload.
    if (local_post_store_failure)
      independent[0] =
          post_store_rank_marker(impl_->mpi->size(), impl_->mpi->rank());
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (impl_->post_store_mpi_fault_rank >= 0) {
      const int target = impl_->post_store_mpi_fault_rank;
      impl_->post_store_mpi_fault_rank = -1;
      const auto preparation = runtime::collective_status(
          *impl_->mpi, impl_->mpi->rank() != target,
          "ideal-gas post-store reduction test fault");
      if (preparation.ok)
        throw runtime::Error(
            "ideal-gas post-store MPI test target was not selected");
      runtime::check_mpi_result(
          MPI_ERR_OTHER,
          "MPI_Allreduce(ideal-gas post-store reduction test fault)");
    }
#endif
    impl_->mpi->allreduce_fp64_in_place(
        independent.data(), independent.size(),
        runtime::Fp64ReductionOperation::maximum);
    ++report.collective_count_;
    report.enthalpy_temperature_max_relative_error_ = independent[0];
    report.eos_max_relative_error_ = independent[1];
    if (independent[0] > 1.0e-12 || independent[1] > 1.0e-12) {
      report.reason_ = IdealGasClosureFailureReason::eos_residual;
      report.lowest_failing_rank_ = decode_post_store_rank_marker(
          impl_->mpi->size(), independent[0]);
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

int IdealGasClosure::prepare_commit() {
  if (!impl_ || !impl_->active || !impl_->latest.authenticated() ||
      impl_->latest.disposition() != IdealGasClosureDisposition::closed ||
      impl_->latest.stage() != IdealGasClosureStage::final ||
      !impl_->latest.final_metrics_available())
    throw runtime::Error("ideal-gas closure commit preparation is invalid");
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (impl_->closure_prepare_fault_rank >= 0) {
    const int target = impl_->closure_prepare_fault_rank;
    impl_->closure_prepare_fault_rank = -1;
    bool prepared = true;
    try {
      impl_->prepared = impl_->trial;
      if (target == impl_->mpi->rank())
        impl_->prepared.revision = std::numeric_limits<std::uint64_t>::max();
      if (impl_->prepared.revision ==
          std::numeric_limits<std::uint64_t>::max())
        throw runtime::Error("ideal-gas closure revision would wrap");
      ++impl_->prepared.revision;
      impl_->prepared_valid = true;
    } catch (...) {
      prepared = false;
    }
    const auto status = runtime::collective_status(
        *impl_->mpi, prepared, "ideal-gas closure commit preparation failed");
    if (!status.ok)
      return status.failing_rank;
    return -1;
  }
#endif
  impl_->prepared = impl_->trial;
  if (impl_->prepared.revision == std::numeric_limits<std::uint64_t>::max())
    throw runtime::Error("ideal-gas closure revision would wrap");
  ++impl_->prepared.revision;
  impl_->prepared_valid = true;
  return -1;
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

void IdealGasClosure::before_outlet(FlowState &state) {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (!impl_ || !impl_->outlet_backflow_fault)
    return;
  impl_->outlet_backflow_fault = false;
  const auto outlet = impl_->boundaries->pressure_outlet_patch_id();
  if (!outlet)
    return;
  const auto &faces = impl_->topology->patch(*outlet).local_faces();
  if (faces.empty())
    return;
  auto &trial = state.solver_layer(FlowLayer::trial);
  auto flux = trial.acquire_face_write<double>(
      state.solver_access_plan(), kStatePhase, kStateActor,
      impl_->fields.face_mass_flux);
  flux(faces.front(), 0) = -std::max(1.0, std::abs(flux(faces.front(), 0)));
#else
  static_cast<void>(state);
#endif
}

int IdealGasClosure::before_prepare(FlowState &state,
                                    AcceptedStepMetadata accepted) {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (impl_ && impl_->state_prepare_fault_rank >= 0) {
    const int target = impl_->state_prepare_fault_rank;
    impl_->state_prepare_fault_rank = -1;
    bool prepared = true;
    if (target == impl_->mpi->rank()) {
      try {
        state.prepare_commit_attempt(accepted);
      } catch (...) {
        prepared = false;
      }
      // A successful target preparation is deliberately rejected so the
      // collective rollback path must undo a genuinely prepared state.
      prepared = false;
    }
    const auto status = runtime::collective_status(
        *impl_->mpi, prepared, "FlowState commit preparation test fault");
    return status.ok ? -1 : status.failing_rank;
  }
#else
  static_cast<void>(state);
  static_cast<void>(accepted);
#endif
  return -1;
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
  return first_state_failure(
             0.0, std::numeric_limits<double>::quiet_NaN(), 1000.0) ==
             IdealGasClosureFailureReason::non_finite_enthalpy &&
         first_state_failure(std::numeric_limits<double>::quiet_NaN(), 0.0,
                             1000.0) ==
             IdealGasClosureFailureReason::non_positive_enthalpy;
}

bool test::IdealGasClosureTestAccess::post_store_rank_marker_is_collision_free(
    int ranks) noexcept {
  if (ranks <= 0)
    return false;
  double selected = 0.0;
  for (int rank = ranks - 1; rank >= 0; --rank) {
    const double marker = post_store_rank_marker(ranks, rank);
    if (!std::isfinite(marker) || marker <= 1.0e-12 ||
        decode_post_store_rank_marker(ranks, marker) != rank)
      return false;
    selected = std::max(selected, marker);
  }
  return decode_post_store_rank_marker(ranks, selected) == 0 &&
         std::max(1.0e-12, selected) == selected;
}

bool test::IdealGasClosureTestAccess::final_gate_rank_marker_is_collision_free(
    int ranks) noexcept {
  if (ranks <= 0)
    return false;
  double selected = 0.0;
  for (std::uint8_t encoded = 1U; encoded <= 5U; ++encoded) {
    const auto origin = static_cast<FinalGateOrigin>(encoded);
    for (int rank = ranks - 1; rank >= 0; --rank) {
      const double marker = final_gate_rank_marker(ranks, origin, rank);
      const auto decoded = decode_final_gate_rank_marker(ranks, marker);
      if (!std::isfinite(marker) || marker <= 0.0 || decoded.origin != origin ||
          decoded.rank != rank)
        return false;
      selected = std::max(selected, marker);
    }
  }
  const auto decoded = decode_final_gate_rank_marker(ranks, selected);
  return decoded.origin == FinalGateOrigin::eos && decoded.rank == 0 &&
         final_gate_rank_marker(ranks, decoded.origin, decoded.rank) ==
             selected;
}

bool test::IdealGasClosureTestAccess::candidate_pressure_mutation_rejected(
    const IdealGasClosureReport &report, bool availability) noexcept {
  auto copy = report;
  if (availability)
    copy.candidate_pressure_available_ =
        !copy.candidate_pressure_available_;
  else
    copy.candidate_pressure_pa_ =
        copy.candidate_pressure_pa_ == 0.0
            ? 1.0
            : std::nextafter(copy.candidate_pressure_pa_,
                             std::numeric_limits<double>::infinity());
  return !copy.authenticated();
}

IdealGasClosure test::IdealGasClosureTestAccess::create(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::MpiContext &mpi, const runtime::FieldRegistry &registry,
    const FlowFieldIds &fields, const FlowState &state,
    IdealGasClosureSpec spec, IdealGasCreateFault fault, int rank) {
  return IdealGasClosure::create_internal(topology, geometry, boundaries, mpi,
                                          registry, fields, state, spec, nullptr,
                                          static_cast<int>(fault), rank);
}

int test::IdealGasClosureTestAccess::preflight_failure_rank(
    const runtime::Error &error) noexcept {
  const auto *failure =
      dynamic_cast<const detail::DensityClosurePreflightFailure *>(&error);
  if (failure != nullptr)
    return failure->failing_rank();
  const auto *validation =
      dynamic_cast<const DensityClosureCreateValidationFailure *>(&error);
  return validation == nullptr ? -1 : validation->failing_rank();
}

std::uint64_t
test::IdealGasClosureTestAccess::preflight_wire_exchange_count(
    const IdealGasClosure &closure) noexcept {
  return closure.impl_ ? closure.impl_->preflight_wire_exchange_count : 0U;
}

IdealGasClosureFailureReason
test::IdealGasClosureTestAccess::create_validation_failure_reason(
    const runtime::Error &error) noexcept {
  const auto *validation =
      dynamic_cast<const DensityClosureCreateValidationFailure *>(&error);
  return validation == nullptr ? IdealGasClosureFailureReason::none
                               : validation->reason();
}

void test::IdealGasClosureTestAccess::begin_attempt(
    IdealGasClosure &closure, FlowState &state, std::uint64_t identity) {
  closure.begin_attempt(state, identity);
}

IdealGasClosureReport test::IdealGasClosureTestAccess::evaluate(
    IdealGasClosure &closure, FlowState &state, IdealGasClosureStage stage) {
  return closure.evaluate(state, stage);
}

void test::IdealGasClosureTestAccess::rollback(
    IdealGasClosure &closure) noexcept {
  closure.rollback();
}

void test::IdealGasClosureTestAccess::set_stage_failure(
    IdealGasClosure &closure, IdealGasClosureStage stage,
    IdealGasClosureFailureReason reason, int rank) {
  closure.set_stage_failure_for_test(stage, reason, rank);
}

void test::IdealGasClosureTestAccess::set_metric_gate_failure(
    IdealGasClosure &closure, IdealGasMetricGateFault kind, int rank) {
  closure.set_metric_gate_failure_for_test(static_cast<std::uint8_t>(kind),
                                           rank);
}

void test::IdealGasClosureTestAccess::set_post_store_corruption(
    IdealGasClosure &closure, int rank, bool enthalpy_density) {
  closure.set_post_store_corruption_for_test(rank, enthalpy_density);
}

void test::IdealGasClosureTestAccess::set_facade_create_fault(
    IdealGasClosure &closure, int rank) {
  if (!closure.impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  closure.impl_->facade_create_fault_rank = rank;
}

void test::IdealGasClosureTestAccess::set_material_factory_create_fault(
    IdealGasClosure &closure, int rank) {
  if (!closure.impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  closure.impl_->material_factory_create_fault_rank = rank;
}

bool test::IdealGasClosureTestAccess::consume_facade_create_fault(
    IdealGasClosure &closure, int rank) noexcept {
  if (!closure.impl_ || closure.impl_->facade_create_fault_rank != rank)
    return false;
  closure.impl_->facade_create_fault_rank = -1;
  return true;
}

int test::IdealGasClosureTestAccess::consume_material_factory_create_fault(
    IdealGasClosure &closure) noexcept {
  if (!closure.impl_)
    return -1;
  const int rank = closure.impl_->material_factory_create_fault_rank;
  closure.impl_->material_factory_create_fault_rank = -1;
  return rank;
}

void test::IdealGasClosureTestAccess::begin_allocation_observation(
    IdealGasClosure &closure) noexcept {
  if (closure.impl_)
    closure.impl_->allocation_observation_active = true;
}

void test::IdealGasClosureTestAccess::end_allocation_observation(
    IdealGasClosure &closure) noexcept {
  if (closure.impl_)
    closure.impl_->allocation_observation_active = false;
}

void test::IdealGasClosureTestAccess::set_attempt_preparation_fault(
    IdealGasClosure &closure, IdealGasAttemptPreparationFault fault,
    int rank) {
  if (!closure.impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  closure.impl_->attempt_preparation_fault_kind = static_cast<int>(fault);
  closure.impl_->attempt_preparation_fault_rank = rank;
}

void test::IdealGasClosureTestAccess::set_controlled_allocation(
    IdealGasClosure &closure, int rank) {
  if (!closure.impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  closure.impl_->controlled_allocation_rank = rank;
}

bool test::IdealGasClosureTestAccess::allocation_observation_active(
    const IdealGasClosure &closure) noexcept {
  return closure.impl_ && closure.impl_->allocation_observation_active;
}

void test::IdealGasClosureTestAccess::consume_attempt_preparation_fault(
    IdealGasClosure &closure) {
  if (!closure.impl_ ||
      closure.impl_->attempt_preparation_fault_rank !=
          closure.impl_->mpi->rank())
    return;
  const auto fault = static_cast<IdealGasAttemptPreparationFault>(
      closure.impl_->attempt_preparation_fault_kind);
  switch (fault) {
  case IdealGasAttemptPreparationFault::pre_authority_preparation:
  case IdealGasAttemptPreparationFault::post_authority_preparation:
    closure.impl_->attempt_preparation_fault_kind = -1;
    closure.impl_->attempt_preparation_fault_rank = -1;
    throw std::bad_alloc();
  default:
    return;
  }
}

void IdealGasClosure::set_post_store_corruption_for_test(
    int rank, bool enthalpy_density) {
  if (!impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  impl_->post_store_fault_rank = rank;
  impl_->post_store_fault_enthalpy = enthalpy_density;
}
void IdealGasClosure::record_halo_for_test(
    std::uint8_t stage, runtime::FieldId density,
    runtime::FieldId enthalpy_density) {
  if (!impl_ ||
      stage < static_cast<std::uint8_t>(IdealGasClosureStage::predictor) ||
      stage > static_cast<std::uint8_t>(IdealGasClosureStage::final))
    throw runtime::Error("ideal-gas Halo trace is invalid");
  impl_->halo_trace.push_back(
      {static_cast<IdealGasClosureStage>(stage), density,
       enthalpy_density});
}

std::vector<std::array<std::uint64_t, 3>>
IdealGasClosure::halo_trace_for_test() const {
  if (!impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  std::vector<std::array<std::uint64_t, 3>> result;
  result.reserve(impl_->halo_trace.size());
  for (const auto &entry : impl_->halo_trace)
    result.push_back(
        {static_cast<std::uint64_t>(entry.stage), entry.density,
         entry.enthalpy_density});
  return result;
}
void IdealGasClosure::set_candidate_precedence_fault_for_test(int rank) {
  if (!impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  impl_->candidate_precedence_fault_rank = rank;
}
void IdealGasClosure::set_stage_failure_for_test(
    IdealGasClosureStage stage, IdealGasClosureFailureReason reason, int rank) {
  if (!impl_ || stage == IdealGasClosureStage::none ||
      reason == IdealGasClosureFailureReason::none)
    throw runtime::Error("ideal-gas closure stage fault is invalid");
  impl_->stage_failure_stage = stage;
  impl_->stage_failure_reason = reason;
  impl_->stage_failure_rank = rank;
}
void IdealGasClosure::set_metric_gate_failure_for_test(std::uint8_t kind,
                                                       int rank) {
  if (!impl_ ||
      kind >
          static_cast<std::uint8_t>(test::IdealGasMetricGateFault::enthalpy) ||
      rank < 0 || rank >= impl_->mpi->size())
    throw runtime::Error("ideal-gas closure metric-gate fault is invalid");
  impl_->metric_gate_failure_kind = static_cast<int>(kind);
  impl_->metric_gate_failure_rank = rank;
}
void IdealGasClosure::set_post_assessment_fault_for_test(std::uint8_t kind,
                                                         int rank) {
  if (!impl_ ||
      kind > static_cast<std::uint8_t>(
                 test::IdealGasPostAssessmentFault::non_positive_density) ||
      rank < 0 || rank >= impl_->mpi->size())
    throw runtime::Error("ideal-gas post-assessment fault is invalid");
  impl_->post_assessment_fault_kind = static_cast<int>(kind);
  impl_->post_assessment_fault_rank = rank;
}
void IdealGasClosure::before_post_assessment_for_test(FlowState &state) {
  if (!impl_ || impl_->post_assessment_fault_rank < 0)
    return;
  const int kind = impl_->post_assessment_fault_kind;
  const int target = impl_->post_assessment_fault_rank;
  impl_->post_assessment_fault_kind = -1;
  impl_->post_assessment_fault_rank = -1;
  if (target != impl_->mpi->rank())
    return;
  auto &trial = state.solver_layer(FlowLayer::trial);
  if (kind == static_cast<int>(
                  test::IdealGasPostAssessmentFault::non_finite_state) ||
      kind == static_cast<int>(
                  test::IdealGasPostAssessmentFault::non_positive_density)) {
    auto density =
        trial.acquire_write<double>(state.solver_access_plan(), kStatePhase,
                                    kStateActor, impl_->fields.density);
    density(0, 0, 0, 0) =
        kind == static_cast<int>(
                    test::IdealGasPostAssessmentFault::non_finite_state)
            ? std::numeric_limits<double>::quiet_NaN()
            : 0.0;
    return;
  }
}
void IdealGasClosure::set_outer_failure_for_test(std::uint8_t point, int rank) {
  if (!impl_ || point > 1U || rank < 0 || rank >= impl_->mpi->size())
    throw runtime::Error("ideal-gas outer failure test control is invalid");
  impl_->outer_failure_point = point;
  impl_->outer_failure_rank = rank;
}
int IdealGasClosure::outer_failure_for_test(std::uint8_t point) {
  if (!impl_ || impl_->outer_failure_point != point)
    return -1;
  const int rank = impl_->outer_failure_rank;
  impl_->outer_failure_point = std::numeric_limits<std::uint8_t>::max();
  impl_->outer_failure_rank = -1;
  return rank;
}
void IdealGasClosure::set_prepare_fault_for_test(bool state_prepare, int rank) {
  if (!impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  (state_prepare ? impl_->state_prepare_fault_rank
                 : impl_->closure_prepare_fault_rank) = rank;
}
void IdealGasClosure::set_post_store_mpi_fault_for_test(int rank) {
  if (!impl_ || rank < 0 || rank >= impl_->mpi->size())
    throw runtime::Error("ideal-gas closure post-store MPI fault is invalid");
  impl_->post_store_mpi_fault_rank = rank;
}
void IdealGasClosure::set_attempt_layout_fault_for_test(int rank) {
  if (!impl_ || rank < 0 || rank >= impl_->mpi->size())
    throw runtime::Error("ideal-gas closure attempt layout fault is invalid");
  impl_->attempt_layout_fault_rank = rank;
}
void IdealGasClosure::set_outlet_backflow_fault_for_test() {
  if (!impl_)
    throw runtime::Error("ideal-gas closure has been moved from");
  impl_->outlet_backflow_fault = true;
}
#endif

} // namespace hundun::flow
