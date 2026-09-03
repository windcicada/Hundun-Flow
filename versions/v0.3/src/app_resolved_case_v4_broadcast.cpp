// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "app_resolved_case_v4_broadcast_detail.hpp"

#include "app_resolved_case_v3_broadcast_detail.hpp"
#include "cfg_resolved_case_v4_loader_detail.hpp"

#include "hundun/cfg_resolved_case_v4_loader.hpp"
#include "hundun/rt_error.hpp"
#include "rt_mpi_error_detail.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::config {
namespace {

using runtime::Error;

enum class BroadcastFault {
  none,
#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
  reverse_species,
  rank_local_error
#endif
};

struct Identity final {
  int rank{};
  int size{};
  int root{};
};

void broadcast_bytes(MPI_Comm comm, int root, char *data, std::uint64_t size,
                     std::string_view operation) {
  constexpr std::uint64_t chunk =
      static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  for (std::uint64_t offset = 0U; offset < size;) {
    const int count = static_cast<int>(std::min(size - offset, chunk));
    runtime::detail::check_mpi(
        MPI_Bcast(data + static_cast<std::size_t>(offset), count, MPI_BYTE,
                  root, comm),
        operation);
    offset += static_cast<std::uint64_t>(count);
  }
}

[[noreturn]] void uniform_error(MPI_Comm comm, const Identity &identity,
                                bool local_ok, std::string local_message) {
  const int candidate = local_ok ? identity.size : identity.rank;
  int failing_rank = identity.size;
  runtime::detail::check_mpi(
      MPI_Allreduce(&candidate, &failing_rank, 1, MPI_INT, MPI_MIN, comm),
      "MPI_Allreduce reacting resolved-case failing rank");
  std::uint64_t length = identity.rank == failing_rank
                             ? static_cast<std::uint64_t>(local_message.size())
                             : 0U;
  runtime::detail::check_mpi(
      MPI_Bcast(&length, 1, MPI_UINT64_T, failing_rank, comm),
      "MPI_Bcast reacting resolved-case error length");
  if (length == 0U || length > static_cast<std::uint64_t>(
                                   std::numeric_limits<std::size_t>::max())) {
    throw Error("reacting resolved-case collective error has invalid length");
  }
  local_message.resize(static_cast<std::size_t>(length));
  broadcast_bytes(comm, failing_rank, local_message.data(), length,
                  "MPI_Bcast reacting resolved-case error bytes");
  throw Error(std::move(local_message));
}

void converge(MPI_Comm comm, const Identity &identity, bool local_ok,
              std::string local_message) {
  const int candidate = local_ok ? identity.size : identity.rank;
  int failing_rank = identity.size;
  runtime::detail::check_mpi(
      MPI_Allreduce(&candidate, &failing_rank, 1, MPI_INT, MPI_MIN, comm),
      "MPI_Allreduce reacting resolved-case status");
  if (failing_rank != identity.size) {
    uniform_error(comm, identity, local_ok, std::move(local_message));
  }
}

Identity collective_identity(MPI_Comm comm, int root) {
  runtime::detail::require_mpi_active("broadcast reacting resolved case");
  if (comm == MPI_COMM_NULL) {
    throw Error(
        "reacting resolved-case broadcast requires an intracommunicator");
  }
  int intercommunicator = 0;
  runtime::detail::check_mpi(MPI_Comm_test_inter(comm, &intercommunicator),
                             "MPI_Comm_test_inter reacting resolved case");
  if (intercommunicator != 0) {
    throw Error(
        "reacting resolved-case broadcast requires an intracommunicator");
  }
  Identity identity{};
  runtime::detail::check_mpi(MPI_Comm_rank(comm, &identity.rank),
                             "MPI_Comm_rank reacting resolved case");
  runtime::detail::check_mpi(MPI_Comm_size(comm, &identity.size),
                             "MPI_Comm_size reacting resolved case");
  int minimum = root;
  int maximum = root;
  runtime::detail::check_mpi(
      MPI_Allreduce(&root, &minimum, 1, MPI_INT, MPI_MIN, comm),
      "MPI_Allreduce reacting resolved-case minimum root");
  runtime::detail::check_mpi(
      MPI_Allreduce(&root, &maximum, 1, MPI_INT, MPI_MAX, comm),
      "MPI_Allreduce reacting resolved-case maximum root");
  if (minimum != maximum) {
    throw Error("reacting resolved-case broadcast root differs across ranks");
  }
  if (minimum < 0 || minimum >= identity.size) {
    throw Error(
        "reacting resolved-case broadcast root is outside communicator");
  }
  identity.root = minimum;
  return identity;
}

void broadcast_string(MPI_Comm comm, const Identity &identity,
                      std::string &value, std::string_view label) {
  std::uint64_t length = identity.rank == identity.root
                             ? static_cast<std::uint64_t>(value.size())
                             : 0U;
  runtime::detail::check_mpi(
      MPI_Bcast(&length, 1, MPI_UINT64_T, identity.root, comm),
      std::string("MPI_Bcast reacting resolved-case ") + std::string(label) +
          " length");
  bool ok = length <=
            static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max());
  std::string message;
  if (ok) {
    try {
      value.resize(static_cast<std::size_t>(length));
    } catch (const std::exception &error) {
      ok = false;
      message = error.what();
    }
  } else {
    message = "reacting resolved-case string exceeds host size domain";
  }
  converge(comm, identity, ok, std::move(message));
  broadcast_bytes(comm, identity.root, value.data(), length,
                  std::string("MPI_Bcast reacting resolved-case ") +
                      std::string(label));
}

template <class T>
void broadcast_pod(MPI_Comm comm, const Identity &identity, T &value,
                   std::string_view label) {
  static_assert(std::is_trivially_copyable_v<T>);
  runtime::detail::check_mpi(MPI_Bcast(&value, static_cast<int>(sizeof(T)),
                                       MPI_BYTE, identity.root, comm),
                             std::string("MPI_Bcast reacting resolved-case ") +
                                 std::string(label));
}

std::uint64_t local_species_fingerprint(const std::vector<std::string> &names) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const auto &name : names) {
    for (const char character : name) {
      hash ^= static_cast<unsigned char>(character);
      hash *= UINT64_C(1099511628211);
    }
    hash ^= UINT64_C(255);
    hash *= UINT64_C(1099511628211);
  }
  return hash == 0U ? 1U : hash;
}

ResolvedReactingCaseV4 broadcast_impl(MPI_Comm comm, int root,
                                      const ResolvedReactingCaseV4 *root_case,
                                      BroadcastFault fault) {
  const Identity identity = collective_identity(comm, root);
  const bool pointer_ok = identity.rank == identity.root ? root_case != nullptr
                                                         : root_case == nullptr;
  converge(comm, identity, pointer_ok,
           identity.rank == identity.root
               ? "reacting resolved-case broadcast root requires case data"
               : "reacting resolved-case broadcast non-root data must be null");

  bool preparation_ok = true;
  std::string preparation_error;
  if (identity.rank == identity.root) {
    try {
      static_cast<void>(to_resolved_reacting_json_v4(*root_case));
    } catch (const std::exception &error) {
      preparation_ok = false;
      preparation_error = error.what();
    }
  }
  converge(comm, identity, preparation_ok, std::move(preparation_error));

  ResolvedReactingCaseV4 result{};
  ResolvedCaseV3 common_root;
  if (identity.rank == identity.root) {
    FlowCaseConfig common_projection = root_case->common_flow;
    common_projection.density_model = DensityModel::constant;
    if (root_case->immersed_boundary.model != ImmersedBoundaryModel::none ||
        root_case->les.model != LesModel::none) {
      ImmersedFlowCaseConfig stage3{};
      stage3.schema_version = 3;
      stage3.common_flow = std::move(common_projection);
      stage3.immersed_boundary = root_case->immersed_boundary;
      stage3.les = root_case->les;
      common_root = ResolvedCaseV3(std::move(stage3));
    } else {
      common_root = ResolvedCaseV3(std::move(common_projection));
    }
  }
  const ResolvedCaseV3 common = broadcast_resolved_case_v3(
      comm, identity.root,
      identity.rank == identity.root ? &common_root : nullptr);
  if (const auto *stage3 = std::get_if<ImmersedFlowCaseConfig>(&common)) {
    result.common_flow = stage3->common_flow;
    result.immersed_boundary = stage3->immersed_boundary;
    result.les = stage3->les;
  } else {
    result.common_flow = std::get<FlowCaseConfig>(common);
    result.immersed_boundary.model = ImmersedBoundaryModel::none;
    result.les.model = LesModel::none;
  }
  result.common_flow.density_model = DensityModel::ideal_gas;

  int schema_version =
      identity.rank == identity.root ? root_case->schema_version : 0;
  broadcast_pod(comm, identity, schema_version, "schema version");
  result.schema_version = schema_version;
  std::string mechanism_file = identity.rank == identity.root
                                   ? root_case->mechanism.file.generic_string()
                                   : std::string();
  std::string mechanism_sha = identity.rank == identity.root
                                  ? root_case->mechanism.sha256
                                  : std::string();
  std::string mechanism_phase = identity.rank == identity.root
                                    ? root_case->mechanism.phase
                                    : std::string();
  broadcast_string(comm, identity, mechanism_file, "mechanism file");
  broadcast_string(comm, identity, mechanism_sha, "mechanism hash");
  broadcast_string(comm, identity, mechanism_phase, "mechanism phase");
  result.mechanism = {mechanism_file, mechanism_sha, mechanism_phase};

  result.chemistry = identity.rank == identity.root ? root_case->chemistry
                                                    : ChemistrySolverConfig{};
  result.initial_p0_pa =
      identity.rank == identity.root ? root_case->initial_p0_pa : 0.0;
  result.initial_temperature_k =
      identity.rank == identity.root ? root_case->initial_temperature_k : 0.0;
  int pressure_mode = identity.rank == identity.root
                          ? static_cast<int>(root_case->pressure_mode)
                          : 0;
  broadcast_pod(comm, identity, result.chemistry.relative_tolerance,
                "chemistry relative tolerance");
  broadcast_pod(comm, identity, result.chemistry.absolute_tolerance,
                "chemistry absolute tolerance");
  broadcast_pod(comm, identity, result.chemistry.maximum_internal_steps,
                "chemistry maximum internal steps");
  broadcast_pod(comm, identity, result.initial_p0_pa, "initial p0");
  broadcast_pod(comm, identity, result.initial_temperature_k,
                "initial temperature");
  broadcast_pod(comm, identity, pressure_mode, "pressure mode");
  result.pressure_mode = static_cast<PressureConstraintMode>(pressure_mode);

  std::uint64_t species_count =
      identity.rank == identity.root
          ? static_cast<std::uint64_t>(root_case->species_names.size())
          : 0U;
  broadcast_pod(comm, identity, species_count, "species count");
  bool count_ok = species_count > 0U &&
                  species_count <= static_cast<std::uint64_t>(
                                       std::numeric_limits<std::size_t>::max());
  converge(comm, identity, count_ok,
           "reacting resolved-case invalid species count");
  result.species_names.resize(static_cast<std::size_t>(species_count));
  result.initial_mass_fractions.resize(static_cast<std::size_t>(species_count));
  for (std::size_t index = 0; index < result.species_names.size(); ++index) {
    if (identity.rank == identity.root) {
      result.species_names[index] = root_case->species_names[index];
      result.initial_mass_fractions[index] =
          root_case->initial_mass_fractions[index];
    }
    broadcast_string(comm, identity, result.species_names[index],
                     "species name");
    broadcast_pod(comm, identity, result.initial_mass_fractions[index],
                  "mass fraction");
  }
  std::uint64_t expected_fingerprint =
      identity.rank == identity.root ? root_case->composition_fingerprint : 0U;
  broadcast_pod(comm, identity, expected_fingerprint,
                "composition fingerprint");

  for (std::size_t index = 0; index < result.boundary_reacting.size();
       ++index) {
    int species_mode =
        identity.rank == identity.root &&
                root_case->boundary_reacting[index].non_catalytic_impermeable
            ? 1
            : 0;
    int thermal_present =
        identity.rank == identity.root &&
                root_case->boundary_reacting[index].thermal.has_value()
            ? 1
            : 0;
    int thermal_mode =
        identity.rank == identity.root && thermal_present != 0
            ? static_cast<int>(
                  root_case->boundary_reacting[index].thermal->mode)
            : 0;
    int temperature_present = identity.rank == identity.root &&
                                      thermal_present != 0 &&
                                      root_case->boundary_reacting[index]
                                          .thermal->temperature_k.has_value()
                                  ? 1
                                  : 0;
    double temperature =
        identity.rank == identity.root && temperature_present != 0
            ? *root_case->boundary_reacting[index].thermal->temperature_k
            : 0.0;
    broadcast_pod(comm, identity, species_mode, "boundary species mode");
    broadcast_pod(comm, identity, thermal_present, "boundary thermal presence");
    broadcast_pod(comm, identity, thermal_mode, "boundary thermal mode");
    broadcast_pod(comm, identity, temperature_present,
                  "boundary temperature presence");
    broadcast_pod(comm, identity, temperature, "boundary temperature");
    auto &boundary = result.boundary_reacting[index];
    boundary.non_catalytic_impermeable = species_mode != 0;
    if (thermal_present != 0) {
      boundary.thermal.emplace();
      boundary.thermal->mode = static_cast<ReactingThermalMode>(thermal_mode);
      if (temperature_present != 0) {
        boundary.thermal->temperature_k = temperature;
      }
    }
  }
#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
  const int mutation_rank = (identity.root + 1) % identity.size;
  if (fault == BroadcastFault::reverse_species &&
      identity.rank == mutation_rank && result.species_names.size() > 1U) {
    std::reverse(result.species_names.begin(), result.species_names.end());
  }
  const bool injected_ok = fault != BroadcastFault::rank_local_error ||
                           identity.rank != mutation_rank;
  converge(comm, identity, injected_ok,
           "injected rank-local reacting resolved-case error");
#else
  static_cast<void>(fault);
#endif

  result.composition_fingerprint =
      local_species_fingerprint(result.species_names);
  const bool species_order_ok =
      result.composition_fingerprint == expected_fingerprint;
  converge(
      comm, identity, species_order_ok,
      "reacting resolved-case species order does not match root fingerprint");

  bool validation_ok = true;
  std::string validation_error;
  try {
    const std::string canonical = to_resolved_reacting_json_v4(result);
    result = detail::parse_resolved_reacting_case_v4_json(
        canonical, "reacting-broadcast-case.json");
  } catch (const std::exception &error) {
    validation_ok = false;
    validation_error = error.what();
  }
  converge(comm, identity, validation_ok, std::move(validation_error));
  return result;
}

} // namespace

ResolvedReactingCaseV4
broadcast_resolved_reacting_case_v4(MPI_Comm comm, int root,
                                    const ResolvedReactingCaseV4 *root_case) {
  return broadcast_impl(comm, root, root_case, BroadcastFault::none);
}

#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
namespace detail {

ResolvedReactingCaseV4 broadcast_resolved_reacting_case_v4_with_fault(
    MPI_Comm comm, int root, const ResolvedReactingCaseV4 *root_case,
    ResolvedReactingCaseV4BroadcastFault fault) {
  BroadcastFault internal = BroadcastFault::none;
  switch (fault) {
  case ResolvedReactingCaseV4BroadcastFault::none:
    break;
  case ResolvedReactingCaseV4BroadcastFault::reverse_species:
    internal = BroadcastFault::reverse_species;
    break;
  case ResolvedReactingCaseV4BroadcastFault::rank_local_error:
    internal = BroadcastFault::rank_local_error;
    break;
  }
  return broadcast_impl(comm, root, root_case, internal);
}

} // namespace detail
#endif

} // namespace hundun::config
