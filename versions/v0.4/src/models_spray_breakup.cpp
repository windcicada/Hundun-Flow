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

ParcelId child_id(const SprayParcelState& parent, std::uint64_t accepted_step,
                  std::uint64_t breakup_ordinal,
                  std::uint32_t child_index) noexcept {
  ParcelRandomAddress address;
  address.seed = parent.id.high;
  address.accepted_step = accepted_step;
  address.stream_identity = parent.id.low;
  address.ordinal = breakup_ordinal;
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

// Shared audit for equal and unequal children. Surface energy is an explicit
// ledger alongside the preserved bulk kinetic energy and absolute enthalpy.
BreakupChildReport audit_children(
    const SprayParcelState& parent, double surface_tension_n_per_m,
    double liquid_enthalpy_j_per_kg, double deformation_energy_per_droplet_j,
    double expected_children_multiplicity, BreakupChildReport report) noexcept {
  const double parent_total_mass = parent.droplet_mass_kg * parent.multiplicity;
  BreakupConservationReport& budget = report.conservation;
  budget.parent_total_mass_kg = parent_total_mass;
  budget.parent_multiplicity = parent.multiplicity;
  budget.expected_children_total_multiplicity =
      expected_children_multiplicity;
  budget.parent_temperature_k = parent.temperature_k;
  budget.minimum_child_temperature_k =
      std::numeric_limits<double>::infinity();
  budget.maximum_child_temperature_k = 0.0;

  long double children_mass = 0.0L;
  std::array<long double, 3U> children_momentum{};
  long double children_kinetic = 0.0L;
  long double children_enthalpy = 0.0L;
  long double children_multiplicity = 0.0L;
  long double children_surface = 0.0L;
  for (std::uint32_t index = 0U; index < report.candidate.child_count;
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
            liquid_enthalpy_j_per_kg);
    children_multiplicity +=
        static_cast<long double>(child.multiplicity);
    children_surface +=
        static_cast<long double>(child.multiplicity) *
        static_cast<long double>(surface_tension_n_per_m) *
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
      dot(parent.velocity_m_per_s, parent.velocity_m_per_s);
  budget.parent_bulk_kinetic_energy_j =
      0.5 * parent_total_mass * parent_speed_squared;
  budget.children_bulk_kinetic_energy_j =
      static_cast<double>(children_kinetic);
  budget.bulk_kinetic_energy_residual_j =
      budget.children_bulk_kinetic_energy_j -
      budget.parent_bulk_kinetic_energy_j;
  budget.parent_thermochemical_enthalpy_j =
      parent_total_mass *
      liquid_enthalpy_j_per_kg;
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
        parent_total_mass * parent.velocity_m_per_s[component];
    budget.children_total_momentum_kg_m_per_s[component] =
        static_cast<double>(children_momentum[component]);
    budget.momentum_residual_kg_m_per_s[component] =
        budget.children_total_momentum_kg_m_per_s[component] -
        budget.parent_total_momentum_kg_m_per_s[component];
  }

  budget.parent_spherical_surface_energy_j =
      parent.multiplicity * surface_tension_n_per_m * kPi *
      parent.droplet_diameter_m * parent.droplet_diameter_m;
  budget.children_spherical_surface_energy_j =
      static_cast<double>(children_surface);
  budget.surface_energy_increase_j =
      budget.children_spherical_surface_energy_j -
      budget.parent_spherical_surface_energy_j;
  budget.supplied_deformation_energy_j =
      deformation_energy_per_droplet_j *
      parent.multiplicity;
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
      budget.minimum_child_temperature_k == parent.temperature_k &&
      budget.maximum_child_temperature_k == parent.temperature_k;
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
    child.id = child_id(input.parent, input.accepted_step,
                        input.breakup_ordinal, index);
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

  return audit_children(input.parent, parameters.surface_tension_n_per_m,
      parameters.liquid_absolute_thermochemical_enthalpy_j_per_kg,
      parameters.supplied_deformation_energy_per_parent_droplet_j,
      expected_children_multiplicity, report);
}

BreakupChildReport generate_sgs_children(const SgsChildInput& input) noexcept {
  const auto fail = [](BreakupChildStatus status) {
    auto report = failure(status);
    report.model_id = "conservative_sgs_binary_split_v1";
    return report;
  };
  const auto& trigger = input.trigger;
  const auto& history = trigger.candidate;
  const double positive[]{input.liquid_density_kg_per_m3,
                          input.surface_tension_n_per_m};
  for (double value : positive) {
    if (!std::isfinite(value) || !(value > 0.0))
      return fail(BreakupChildStatus::invalid_input);
  }
  const double diagnostics[]{history.mean_dissipation_m2_per_s3,
      history.dissipation_age_s, history.mean_rate_per_s, history.rate_age_s,
      trigger.weber_number, trigger.critical_diameter_m,
      trigger.kolmogorov_length_m, trigger.deterministic_rate_per_s,
      trigger.stochastic_rate_per_s, trigger.minimum_diameter_ratio,
      trigger.maximum_diameter_ratio, trigger.distribution_lambda};
  for (double value : diagnostics) {
    if (!std::isfinite(value) || value < 0.0)
      return fail(BreakupChildStatus::invalid_input);
  }
  if (validate_parcel_state(input.parent) != ParcelStateStatus::success ||
      !std::isfinite(input.liquid_absolute_thermochemical_enthalpy_j_per_kg) ||
      !trigger.succeeded() ||
      trigger.model_id != "deterministic_sgs_martinez_bazan_20_v1" ||
      history.poisson_multiplier > 7U) {
    return fail(BreakupChildStatus::invalid_input);
  }
  if (!trigger.breakup_requested)
    return fail(BreakupChildStatus::breakup_not_requested);
  if (!trigger.rate_active || history.poisson_multiplier == 0U ||
      !(history.mean_rate_per_s > 0.0) ||
      !(history.rate_age_s >= 1.0 /
          (history.poisson_multiplier * history.mean_rate_per_s))) {
    return fail(BreakupChildStatus::invalid_input);
  }
  const auto ratios = trigger.daughter_diameter_ratios;
  for (double ratio : ratios) {
    if (!std::isfinite(ratio) || !(ratio > 0.0) || ratio > 1.0)
      return fail(BreakupChildStatus::invalid_input);
  }
  const double volume_sum = ratios[0] * ratios[0] * ratios[0] +
                            ratios[1] * ratios[1] * ratios[1];
  if (std::abs(volume_sum - 1.0) >
      64.0 * std::numeric_limits<double>::epsilon()) {
    return fail(BreakupChildStatus::invalid_input);
  }
  if (!parent_geometry_consistent(input.parent, input.liquid_density_kg_per_m3))
    return fail(BreakupChildStatus::inconsistent_parent_geometry);

  // Form the smaller mass first. At a CDF endpoint the larger diameter can
  // round to the parent diameter while the smaller daughter is representable.
  const std::size_t small = ratios[0] <= ratios[1] ? 0U : 1U;
  std::array<double, 2U> masses{};
  masses[small] = input.parent.droplet_mass_kg * ratios[small] *
                  ratios[small] * ratios[small];
  masses[1U - small] = input.parent.droplet_mass_kg - masses[small];
  BreakupChildReport report;
  report.model_id = "conservative_sgs_binary_split_v1";
  report.trigger_model_id = trigger.model_id;
  report.rng_model_id = "parcel_counter_splitmix64_v1:breakup_child";
  report.parent_id = input.parent.id;
  report.accepted_step = input.accepted_step;
  report.breakup_ordinal = input.breakup_ordinal;
  report.candidate.child_count = 2U;
  for (std::uint32_t index = 0U; index < 2U; ++index) {
    auto child = input.parent;
    child.id = child_id(input.parent, input.accepted_step, input.breakup_ordinal, index);
    if ((child.id.high == 0U && child.id.low == 0U) || child.id == input.parent.id ||
        (index == 1U && child.id == report.candidate.children[0].id)) {
      return fail(BreakupChildStatus::child_id_collision);
    }
    child.droplet_diameter_m = input.parent.droplet_diameter_m * ratios[index];
    child.droplet_mass_kg = masses[index];
    if (validate_parcel_state(child) != ParcelStateStatus::success)
      return fail(BreakupChildStatus::non_finite_output);
    report.candidate.children[index] = child;
  }
  report = audit_children(input.parent, input.surface_tension_n_per_m,
      input.liquid_absolute_thermochemical_enthalpy_j_per_kg, 0.0,
      2.0 * input.parent.multiplicity, report);
  if (!report.succeeded()) return fail(report.status);
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
