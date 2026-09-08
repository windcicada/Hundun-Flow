// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_breakup_detail.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace hundun::v04::spray::detail {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kRelativeTolerance = 1.0e-12;
constexpr double kAbsoluteMassToleranceKg = 1.0e-24;
constexpr double kAbsoluteMomentumToleranceKgMPerS = 1.0e-22;
constexpr double kAbsoluteEnergyToleranceJ = 1.0e-18;

bool finite_vector(const Vector3& value) noexcept {
  return std::isfinite(value[0U]) && std::isfinite(value[1U]) &&
         std::isfinite(value[2U]);
}

double dot(const Vector3& left, const Vector3& right) noexcept {
  return left[0U] * right[0U] + left[1U] * right[1U] +
         left[2U] * right[2U];
}

bool within_tolerance(double residual, double left, double right,
                      double absolute_tolerance) noexcept {
  return std::abs(residual) <=
         absolute_tolerance +
             kRelativeTolerance * (std::abs(left) + std::abs(right));
}

BreakupChildReport failure(BreakupChildStatus status) noexcept {
  BreakupChildReport report;
  report.status = status;
  return report;
}

bool parent_geometry_consistent(const SprayParcelState& parent,
                                double liquid_density_kg_per_m3) noexcept {
  const double geometric_mass = liquid_density_kg_per_m3 * kPi / 6.0 *
                                parent.droplet_diameter_m *
                                parent.droplet_diameter_m *
                                parent.droplet_diameter_m;
  return std::isfinite(geometric_mass) && geometric_mass > 0.0 &&
         within_tolerance(geometric_mass - parent.droplet_mass_kg,
                          geometric_mass, parent.droplet_mass_kg,
                          kAbsoluteMassToleranceKg);
}

ParcelId child_id(const BreakupChildInput& input,
                  std::uint32_t child_index) noexcept {
  ParcelRandomAddress address;
  address.seed = input.parent.id.high;
  address.accepted_step = input.accepted_step;
  address.stream_identity = input.parent.id.low;
  address.ordinal = input.breakup_ordinal;
  address.purpose = ParcelRandomPurpose::breakup_child;
  const std::uint64_t first_lane =
      static_cast<std::uint64_t>(child_index) * 2U;
  return {parcel_random_u64(address, first_lane),
          parcel_random_u64(address, first_lane + 1U)};
}

bool valid_tab_trigger(const TabBreakupReport& trigger) noexcept {
  return trigger.succeeded() &&
         trigger.model_id == "tab_linear_oscillator_v1" &&
         std::isfinite(trigger.forcing_per_s2) &&
         trigger.forcing_per_s2 >= 0.0 &&
         std::isfinite(trigger.damping_per_s) &&
         trigger.damping_per_s >= 0.0 &&
         std::isfinite(trigger.stiffness_per_s2) &&
         trigger.stiffness_per_s2 >= 0.0 &&
         std::isfinite(trigger.event_time_s) &&
         trigger.event_time_s >= 0.0 &&
         std::isfinite(trigger.advanced_duration_s) &&
         trigger.advanced_duration_s >= 0.0 &&
         std::isfinite(trigger.candidate.deformation) &&
         std::isfinite(trigger.candidate.deformation_rate_per_s);
}

bool finite_conservation(const BreakupConservationReport& budget) noexcept {
  const double scalars[]{
      budget.parent_total_mass_kg,
      budget.children_total_mass_kg,
      budget.mass_residual_kg,
      budget.parent_bulk_kinetic_energy_j,
      budget.children_bulk_kinetic_energy_j,
      budget.bulk_kinetic_energy_residual_j,
      budget.parent_thermochemical_enthalpy_j,
      budget.children_thermochemical_enthalpy_j,
      budget.thermochemical_enthalpy_residual_j,
      budget.parent_temperature_k,
      budget.minimum_child_temperature_k,
      budget.maximum_child_temperature_k,
      budget.parent_multiplicity,
      budget.expected_children_total_multiplicity,
      budget.children_total_multiplicity,
      budget.multiplicity_residual,
      budget.parent_spherical_surface_energy_j,
      budget.children_spherical_surface_energy_j,
      budget.surface_energy_increase_j,
      budget.supplied_deformation_energy_j,
      budget.unassigned_deformation_energy_j};
  for (double value : scalars) {
    if (!std::isfinite(value)) return false;
  }
  return finite_vector(budget.parent_total_momentum_kg_m_per_s) &&
         finite_vector(budget.children_total_momentum_kg_m_per_s) &&
         finite_vector(budget.momentum_residual_kg_m_per_s);
}

}  // namespace

BreakupChildReport generate_supplied_diameter_children(
    const BreakupChildInput& input) noexcept {
  const SuppliedDiameterSplitParameters& parameters = input.parameters;
  if (validate_parcel_state(input.parent) != ParcelStateStatus::success ||
      !valid_tab_trigger(input.tab_trigger) ||
      parameters.child_parcel_count < 2U ||
      parameters.maximum_child_parcels < 2U ||
      parameters.maximum_child_parcels > kMaximumBreakupChildParcels ||
      parameters.child_parcel_count > parameters.maximum_child_parcels ||
      !std::isfinite(parameters.supplied_uniform_child_diameter_m) ||
      !(parameters.supplied_uniform_child_diameter_m > 0.0) ||
      !(parameters.supplied_uniform_child_diameter_m <
        input.parent.droplet_diameter_m) ||
      !std::isfinite(parameters.liquid_density_kg_per_m3) ||
      !(parameters.liquid_density_kg_per_m3 > 0.0) ||
      !std::isfinite(parameters.surface_tension_n_per_m) ||
      parameters.surface_tension_n_per_m < 0.0 ||
      !std::isfinite(
          parameters.liquid_absolute_thermochemical_enthalpy_j_per_kg) ||
      !std::isfinite(
          parameters.supplied_deformation_energy_per_parent_droplet_j) ||
      parameters.supplied_deformation_energy_per_parent_droplet_j < 0.0) {
    return failure(BreakupChildStatus::invalid_input);
  }
  if (!input.tab_trigger.breakup_requested) {
    return failure(BreakupChildStatus::breakup_not_requested);
  }
  if (!parent_geometry_consistent(input.parent,
                                  parameters.liquid_density_kg_per_m3)) {
    return failure(BreakupChildStatus::inconsistent_parent_geometry);
  }

  const double child_diameter =
      parameters.supplied_uniform_child_diameter_m;
  const double child_mass = parameters.liquid_density_kg_per_m3 * kPi / 6.0 *
                            child_diameter * child_diameter * child_diameter;
  const double parent_total_mass =
      input.parent.droplet_mass_kg * input.parent.multiplicity;
  const double expected_children_multiplicity =
      parent_total_mass / child_mass;
  const double child_multiplicity =
      expected_children_multiplicity /
      static_cast<double>(parameters.child_parcel_count);
  if (!std::isfinite(child_mass) || !(child_mass > 0.0) ||
      !std::isfinite(parent_total_mass) || !(parent_total_mass > 0.0) ||
      !std::isfinite(expected_children_multiplicity) ||
      !(expected_children_multiplicity > 0.0) ||
      !std::isfinite(child_multiplicity) ||
      !(child_multiplicity > 0.0)) {
    return failure(BreakupChildStatus::non_finite_output);
  }

  BreakupChildReport report;
  report.trigger_model_id = input.tab_trigger.model_id;
  report.rng_model_id = "parcel_counter_splitmix64_v1:breakup_child";
  report.parent_id = input.parent.id;
  report.accepted_step = input.accepted_step;
  report.breakup_ordinal = input.breakup_ordinal;
  report.candidate.child_count = parameters.child_parcel_count;

  for (std::uint32_t index = 0U; index < parameters.child_parcel_count;
       ++index) {
    SprayParcelState child = input.parent;
    child.id = child_id(input, index);
    if ((child.id.high == 0U && child.id.low == 0U) ||
        child.id == input.parent.id) {
      return failure(BreakupChildStatus::child_id_collision);
    }
    for (std::uint32_t previous = 0U; previous < index; ++previous) {
      if (child.id == report.candidate.children[previous].id) {
        return failure(BreakupChildStatus::child_id_collision);
      }
    }
    child.droplet_mass_kg = child_mass;
    child.droplet_diameter_m = child_diameter;
    child.multiplicity = child_multiplicity;
    if (validate_parcel_state(child) != ParcelStateStatus::success) {
      return failure(BreakupChildStatus::non_finite_output);
    }
    report.candidate.children[index] = child;
  }

  BreakupConservationReport& budget = report.conservation;
  budget.parent_total_mass_kg = parent_total_mass;
  budget.parent_multiplicity = input.parent.multiplicity;
  budget.expected_children_total_multiplicity =
      expected_children_multiplicity;
  budget.parent_temperature_k = input.parent.temperature_k;
  budget.minimum_child_temperature_k =
      std::numeric_limits<double>::infinity();
  budget.maximum_child_temperature_k = 0.0;

  long double children_mass = 0.0L;
  std::array<long double, 3U> children_momentum{};
  long double children_kinetic = 0.0L;
  long double children_enthalpy = 0.0L;
  long double children_multiplicity = 0.0L;
  long double children_surface = 0.0L;
  for (std::uint32_t index = 0U; index < parameters.child_parcel_count;
       ++index) {
    const SprayParcelState& child = report.candidate.children[index];
    const long double represented_mass =
        static_cast<long double>(child.droplet_mass_kg) *
        static_cast<long double>(child.multiplicity);
    children_mass += represented_mass;
    const double speed_squared =
        dot(child.velocity_m_per_s, child.velocity_m_per_s);
    children_kinetic += 0.5L * represented_mass *
                        static_cast<long double>(speed_squared);
    children_enthalpy +=
        represented_mass *
        static_cast<long double>(
            parameters.liquid_absolute_thermochemical_enthalpy_j_per_kg);
    children_multiplicity +=
        static_cast<long double>(child.multiplicity);
    children_surface +=
        static_cast<long double>(child.multiplicity) *
        static_cast<long double>(parameters.surface_tension_n_per_m) *
        static_cast<long double>(kPi) *
        static_cast<long double>(child.droplet_diameter_m) *
        static_cast<long double>(child.droplet_diameter_m);
    for (std::size_t component = 0U; component < 3U; ++component) {
      children_momentum[component] +=
          represented_mass *
          static_cast<long double>(child.velocity_m_per_s[component]);
    }
    budget.minimum_child_temperature_k =
        std::min(budget.minimum_child_temperature_k, child.temperature_k);
    budget.maximum_child_temperature_k =
        std::max(budget.maximum_child_temperature_k, child.temperature_k);
  }

  budget.children_total_mass_kg = static_cast<double>(children_mass);
  budget.mass_residual_kg =
      budget.children_total_mass_kg - budget.parent_total_mass_kg;
  const double parent_speed_squared =
      dot(input.parent.velocity_m_per_s, input.parent.velocity_m_per_s);
  budget.parent_bulk_kinetic_energy_j =
      0.5 * parent_total_mass * parent_speed_squared;
  budget.children_bulk_kinetic_energy_j =
      static_cast<double>(children_kinetic);
  budget.bulk_kinetic_energy_residual_j =
      budget.children_bulk_kinetic_energy_j -
      budget.parent_bulk_kinetic_energy_j;
  budget.parent_thermochemical_enthalpy_j =
      parent_total_mass *
      parameters.liquid_absolute_thermochemical_enthalpy_j_per_kg;
  budget.children_thermochemical_enthalpy_j =
      static_cast<double>(children_enthalpy);
  budget.thermochemical_enthalpy_residual_j =
      budget.children_thermochemical_enthalpy_j -
      budget.parent_thermochemical_enthalpy_j;
  budget.children_total_multiplicity =
      static_cast<double>(children_multiplicity);
  budget.multiplicity_residual =
      budget.children_total_multiplicity -
      budget.expected_children_total_multiplicity;
  for (std::size_t component = 0U; component < 3U; ++component) {
    budget.parent_total_momentum_kg_m_per_s[component] =
        parent_total_mass * input.parent.velocity_m_per_s[component];
    budget.children_total_momentum_kg_m_per_s[component] =
        static_cast<double>(children_momentum[component]);
    budget.momentum_residual_kg_m_per_s[component] =
        budget.children_total_momentum_kg_m_per_s[component] -
        budget.parent_total_momentum_kg_m_per_s[component];
  }

  budget.parent_spherical_surface_energy_j =
      input.parent.multiplicity * parameters.surface_tension_n_per_m * kPi *
      input.parent.droplet_diameter_m * input.parent.droplet_diameter_m;
  budget.children_spherical_surface_energy_j =
      static_cast<double>(children_surface);
  budget.surface_energy_increase_j =
      budget.children_spherical_surface_energy_j -
      budget.parent_spherical_surface_energy_j;
  budget.supplied_deformation_energy_j =
      parameters.supplied_deformation_energy_per_parent_droplet_j *
      input.parent.multiplicity;
  budget.unassigned_deformation_energy_j =
      budget.supplied_deformation_energy_j -
      budget.surface_energy_increase_j;
  const double energy_tolerance =
      kAbsoluteEnergyToleranceJ +
      kRelativeTolerance *
          (std::abs(budget.supplied_deformation_energy_j) +
           std::abs(budget.surface_energy_increase_j));
  if (std::abs(budget.unassigned_deformation_energy_j) <=
      energy_tolerance) {
    budget.energy_disposition = BreakupEnergyDisposition::balanced;
  } else if (budget.unassigned_deformation_energy_j < 0.0) {
    budget.energy_disposition =
        BreakupEnergyDisposition::supplied_deformation_energy_deficit;
  } else {
    budget.energy_disposition =
        BreakupEnergyDisposition::supplied_deformation_energy_surplus;
  }

  if (!finite_conservation(budget)) {
    return failure(BreakupChildStatus::non_finite_output);
  }
  bool conservative = within_tolerance(
      budget.mass_residual_kg, budget.children_total_mass_kg,
      budget.parent_total_mass_kg, kAbsoluteMassToleranceKg);
  conservative &= within_tolerance(
      budget.bulk_kinetic_energy_residual_j,
      budget.children_bulk_kinetic_energy_j,
      budget.parent_bulk_kinetic_energy_j, kAbsoluteEnergyToleranceJ);
  conservative &= within_tolerance(
      budget.thermochemical_enthalpy_residual_j,
      budget.children_thermochemical_enthalpy_j,
      budget.parent_thermochemical_enthalpy_j,
      kAbsoluteEnergyToleranceJ);
  conservative &= within_tolerance(
      budget.multiplicity_residual,
      budget.children_total_multiplicity,
      budget.expected_children_total_multiplicity,
      std::numeric_limits<double>::epsilon());
  conservative &=
      budget.minimum_child_temperature_k == input.parent.temperature_k &&
      budget.maximum_child_temperature_k == input.parent.temperature_k;
  for (std::size_t component = 0U; component < 3U; ++component) {
    conservative &= within_tolerance(
        budget.momentum_residual_kg_m_per_s[component],
        budget.children_total_momentum_kg_m_per_s[component],
        budget.parent_total_momentum_kg_m_per_s[component],
        kAbsoluteMomentumToleranceKgMPerS);
  }
  if (!conservative) {
    return failure(BreakupChildStatus::conservation_failure);
  }

  report.candidate.available = true;
  report.status = BreakupChildStatus::success;
  return report;
}

TabRepresentativeSplitReport generate_tab_representative_children(
    const TabRepresentativeSplitInput& input) noexcept {
  TabRepresentativeSplitReport report;
  const double radius = 0.5 * input.parent.droplet_diameter_m;
  const double rho = input.liquid_density_kg_per_m3;
  const double sigma = input.surface_tension_n_per_m;
  if (!valid_tab_trigger(input.tab_trigger) ||
      !input.tab_trigger.breakup_requested ||
      input.child_parcel_count < 2U || input.child_parcel_count % 2U != 0U ||
      input.child_parcel_count > kMaximumBreakupChildParcels ||
      !std::isfinite(rho) || rho <= 0.0 || !std::isfinite(sigma) || sigma <= 0.0 ||
      !finite_vector(input.breakup_axis) ||
      !std::isfinite(radius) || radius <= 0.0 ||
      std::abs(input.tab_trigger.candidate.deformation - 1.0) > 1e-10) return report;
  const double stiffness = 8.0 * sigma / (rho * radius * radius * radius);
  if (!within_tolerance(input.tab_trigger.stiffness_per_s2 - stiffness,
                        input.tab_trigger.stiffness_per_s2, stiffness, 1e-10))
    return report;
  const double axis_length = std::sqrt(dot(input.breakup_axis, input.breakup_axis));
  if (!std::isfinite(axis_length) || axis_length <= 0.0) return report;
  Vector3 axis = input.breakup_axis;
  for (double& value : axis) value /= axis_length;
  // Choose a well-conditioned transverse basis deterministically. The RNG
  // address is used only for child identities, not a new dispersion model.
  const Vector3 reference = std::abs(axis[0]) < 0.8 ? Vector3{1, 0, 0}
                                                  : Vector3{0, 1, 0};
  Vector3 normal{axis[1]*reference[2]-axis[2]*reference[1],
                 axis[2]*reference[0]-axis[0]*reference[2],
                 axis[0]*reference[1]-axis[1]*reference[0]};
  const double normal_length = std::sqrt(dot(normal, normal));
  for (double& value : normal) value /= normal_length;
  const double rate = input.tab_trigger.candidate.deformation_rate_per_s;
  const double inverse_size_ratio = 7.0/3.0 +
      rho * radius * radius * radius * rate * rate / (8.0 * sigma);
  const double energy_per_drop = (10.0/3.0) * kPi/5.0 * rho *
      std::pow(radius, 5) * (rate * rate + stiffness);
  BreakupChildInput supplied;
  supplied.parent = input.parent;
  supplied.tab_trigger = input.tab_trigger;
  supplied.accepted_step = input.accepted_step;
  supplied.breakup_ordinal = input.breakup_ordinal;
  supplied.parameters = {input.child_parcel_count, kMaximumBreakupChildParcels,
      2.0 * radius / inverse_size_ratio, rho, sigma,
      input.liquid_absolute_thermochemical_enthalpy_j_per_kg, energy_per_drop};
  auto split = generate_supplied_diameter_children(supplied);
  if (!split.succeeded()) { report.status = split.status; return report; }
  const double transverse_speed = 0.5 * radius * rate;
  long double kinetic = 0.0L;
  std::array<long double, 3> child_momentum{};
  for (std::uint32_t i = 0; i < split.candidate.child_count; ++i) {
    auto& child = split.candidate.children[i];
    for (std::size_t d = 0; d < 3; ++d)
      child.velocity_m_per_s[d] += (i % 2U ? -1.0 : 1.0) * transverse_speed * normal[d];
    kinetic += 0.5L * child.droplet_mass_kg * child.multiplicity *
        dot(child.velocity_m_per_s, child.velocity_m_per_s);
    for (std::size_t d = 0; d < 3; ++d)
      child_momentum[d] += static_cast<long double>(child.droplet_mass_kg) *
          child.multiplicity * child.velocity_m_per_s[d];
    if (validate_parcel_state(child) != ParcelStateStatus::success) return report;
  }
  auto& budget = split.conservation;
  for (std::size_t d = 0; d < 3; ++d) {
    budget.children_total_momentum_kg_m_per_s[d] = static_cast<double>(child_momentum[d]);
    budget.momentum_residual_kg_m_per_s[d] =
        budget.children_total_momentum_kg_m_per_s[d] - budget.parent_total_momentum_kg_m_per_s[d];
    if (!within_tolerance(budget.momentum_residual_kg_m_per_s[d],
                          budget.children_total_momentum_kg_m_per_s[d],
                          budget.parent_total_momentum_kg_m_per_s[d],
                          kAbsoluteMomentumToleranceKgMPerS)) {
      report.status = BreakupChildStatus::conservation_failure;
      return report;
    }
  }
  budget.children_bulk_kinetic_energy_j = static_cast<double>(kinetic);
  budget.bulk_kinetic_energy_residual_j = budget.children_bulk_kinetic_energy_j -
      budget.parent_bulk_kinetic_energy_j;
  const double dispersion = 0.5 * budget.parent_total_mass_kg *
      transverse_speed * transverse_speed;
  const double residual = budget.surface_energy_increase_j + dispersion -
      budget.supplied_deformation_energy_j;
  if (!std::isfinite(residual) ||
      !within_tolerance(residual, budget.surface_energy_increase_j + dispersion,
                        budget.supplied_deformation_energy_j, kAbsoluteEnergyToleranceJ) ||
      !within_tolerance(budget.bulk_kinetic_energy_residual_j - dispersion,
                        budget.children_bulk_kinetic_energy_j,
                        budget.parent_bulk_kinetic_energy_j, kAbsoluteEnergyToleranceJ)) {
    report.status = BreakupChildStatus::conservation_failure;
    return report;
  }
  // No heat or gas source absorbs a residual. Dispersion consumes the exact
  // non-surface part of the oscillation budget; children start undeformed.
  budget.unassigned_deformation_energy_j = -residual;
  budget.energy_disposition = BreakupEnergyDisposition::balanced;
  split.model_id = report.model_id;
  report.status = BreakupChildStatus::success;
  report.representative_diameter_m = supplied.parameters.supplied_uniform_child_diameter_m;
  report.transverse_speed_m_per_s = std::abs(transverse_speed);
  report.deformation_energy_j = budget.supplied_deformation_energy_j;
  report.transverse_kinetic_energy_j = dispersion;
  report.total_energy_residual_j = residual;
  report.split = split;
  return report;
}

}  // namespace hundun::v04::spray::detail
