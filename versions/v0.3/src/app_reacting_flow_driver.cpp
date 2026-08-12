// SPDX-License-Identifier: Apache-2.0

#include "app_reacting_flow_driver_detail.hpp"

#include "hundun/rt_error.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <set>
#include <stdexcept>

namespace hundun::application {
namespace {

bool lowercase_hex_sha256(const std::string &value) {
  return value.size() == 64U &&
         std::all_of(value.begin(), value.end(), [](char character) {
           return (character >= '0' && character <= '9') ||
                  (character >= 'a' && character <= 'f');
         });
}

void validate_reacting_case(const config::ResolvedReactingCaseV4 &resolved) {
  if (resolved.schema_version != 4 || resolved.common_flow.schema_version != 2 ||
      resolved.mechanism.file.empty() ||
      !lowercase_hex_sha256(resolved.mechanism.sha256) ||
      resolved.mechanism.phase.empty() || !std::isfinite(resolved.initial_p0_pa) ||
      resolved.initial_p0_pa <= 0.0 ||
      !std::isfinite(resolved.initial_temperature_k) ||
      resolved.initial_temperature_k <= 0.0 ||
      resolved.species_names.empty() ||
      resolved.species_names.size() != resolved.initial_mass_fractions.size() ||
      resolved.composition_fingerprint == 0U)
    throw std::invalid_argument("reacting driver case identity is invalid");
  std::set<std::string> species;
  double sum{};
  for (std::size_t k = 0; k < resolved.species_names.size(); ++k) {
    if (resolved.species_names[k].empty() ||
        !species.insert(resolved.species_names[k]).second ||
        !std::isfinite(resolved.initial_mass_fractions[k]) ||
        resolved.initial_mass_fractions[k] < 0.0)
      throw std::invalid_argument("reacting driver species identity is invalid");
    sum += resolved.initial_mass_fractions[k];
  }
  if (std::abs(sum - 1.0) > 1.0e-12)
    throw std::invalid_argument("reacting driver mass fractions are invalid");
  for (const auto &boundary : resolved.boundary_reacting)
    if (boundary.thermal && !boundary.non_catalytic_impermeable)
      throw std::invalid_argument("catalytic reacting boundary is unsupported");
  const bool immersed = resolved.immersed_boundary.model !=
                        config::ImmersedBoundaryModel::none;
  if (immersed != resolved.immersed_boundary.geometry.has_value())
    throw std::invalid_argument("reacting driver IBM configuration is incomplete");
  const bool wale = resolved.les.model == config::LesModel::wale;
  if (wale != resolved.les.wale.has_value())
    throw std::invalid_argument("reacting driver WALE configuration is incomplete");
}

} // namespace

ReactingDriverPlan
plan_reacting_flow_case(const config::ResolvedReactingCaseV4 &resolved) {
  validate_reacting_case(resolved);
  ReactingDriverPlan plan;
  plan.closed_pressure = resolved.pressure_mode !=
                         config::PressureConstraintMode::open_fixed_p0;
  plan.immersed_boundary = resolved.immersed_boundary.model !=
                           config::ImmersedBoundaryModel::none;
  plan.wale = resolved.les.model == config::LesModel::wale;
  plan.backend_runtime_count = 1U;
  plan.workspace_pool_count = 1U;
  plan.operator_order = {"chemistry-half-1", "scalar-transport", "piso-1",
                         "chemistry-half-2", "piso-2"};
  if (plan.closed_pressure)
    plan.operator_order.push_back("p0-constraint");
  return plan;
}

int run_reacting_flow_case(const config::ResolvedReactingCaseV4 &resolved,
                           MPI_Comm communicator) {
  if (communicator == MPI_COMM_NULL)
    throw runtime::Error("reacting driver communicator is null");
  const auto plan = plan_reacting_flow_case(resolved);
  int local = plan.backend_runtime_count == 1U &&
                      plan.workspace_pool_count == 1U
                  ? 1
                  : 0;
  int global{};
  if (MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, communicator) !=
      MPI_SUCCESS)
    throw runtime::Error("reacting driver collective validation failed");
  if (global != 1)
    throw runtime::Error("reacting driver construction is inconsistent");
  return EXIT_SUCCESS;
}

} // namespace hundun::application
