// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_physics.hpp"

#include "../../src/physics_input_detail.hpp"

#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <locale>
#include <new>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void observe() noexcept {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
}

void* allocate(std::size_t size) {
  observe();
  void* result = std::malloc(size == 0U ? 1U : size);
  if (result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}

void* allocate_aligned(std::size_t size, std::size_t alignment) {
  observe();
  void* result = nullptr;
  const std::size_t requested = size == 0U ? alignment : size;
  if (posix_memalign(&result, alignment, requested) != 0 || result == nullptr) {
    throw std::bad_alloc{};
  }
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

void* operator new(std::size_t size) {
  return allocation_observer::allocate(size);
}
void* operator new[](std::size_t size) {
  return allocation_observer::allocate(size);
}
void* operator new(std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
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

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool near(double actual, double expected, double tolerance) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <=
             tolerance * std::max(1.0, std::abs(expected));
}

class CommaDecimalPoint final : public std::numpunct<char> {
 protected:
  char do_decimal_point() const override { return ','; }
};

class ScopedGlobalLocale final {
 public:
  explicit ScopedGlobalLocale(const std::locale& locale)
      : previous_(std::locale::global(locale)) {}
  ~ScopedGlobalLocale() { std::locale::global(previous_); }

  ScopedGlobalLocale(const ScopedGlobalLocale&) = delete;
  ScopedGlobalLocale& operator=(const ScopedGlobalLocale&) = delete;

 private:
  std::locale previous_;
};

constexpr std::string_view kThermophysicalText = R"data(
# strict direct-root thermophysical table
HUNDUN_THERMOPHYSICS_V1
temperature_bounds 200 3000
temperature_inversion 1e-12 64
closed_mass_newton 1e-12 32 0.2
species_count 1
species air
molecular_weight 28.96546
temperature_switch 1000
nasa7_low 3.5 0 0 0 0 0 0
nasa7_high 3.5 0 0 0 0 0 0
transport_sutherland 1.716e-5 273.15 110.4 0.71
end_species
end
)data";

constexpr std::string_view kCoastNativeAirText = R"data(
HUNDUN_THERMOPHYSICS_V1
temperature_bounds 200 6000
temperature_inversion 1e-12 64
closed_mass_newton 1e-12 32 0.2
species_count 1
species air
molecular_weight 28.850334
temperature_switch 1000
nasa7_low 3.5838100068 -7.2700635412e-4 1.67056387003e-6 -1.091801341e-10 -4.317787988e-13 -1050.5394088 3.1124135035
nasa7_high 3.1013370688 1.24138813631e-3 -4.1882038804e-7 6.641656204e-11 -3.9127843272e-15 -985.27467132 5.3560174057
transport_coast_native_air
end_species
end
)data";

SpeciesThermophysicalSpec constant_species(std::string_view name, double mw,
                                           double cp) {
  SpeciesThermophysicalSpec species;
  species.stable_name.assign(name.data(), name.size());
  species.molecular_weight = mw;
  species.temperature_switch = 1000.0;
  const double r = kUniversalGasConstant / mw;
  species.nasa7_low[0U] = cp / r;
  species.nasa7_high[0U] = cp / r;
  species.viscosity_reference = 1.0e-5;
  species.conductivity = 0.02;
  return species;
}

SpeciesThermophysicalSpec varying_species(std::string_view name, double mw,
                                          double a1, double a2) {
  SpeciesThermophysicalSpec species = constant_species(name, mw, 1000.0);
  species.nasa7_low = {a1, a2, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double high_a2 = 0.5 * a2;
  const double switch_temperature = species.temperature_switch;
  const double high_a1 = a1 + (a2 - high_a2) * switch_temperature;
  const double low_h_at_switch =
      a1 * switch_temperature + 0.5 * a2 * switch_temperature *
                                    switch_temperature;
  const double high_h_without_offset =
      high_a1 * switch_temperature +
      0.5 * high_a2 * switch_temperature * switch_temperature;
  species.nasa7_high = {high_a1, high_a2, 0.0, 0.0, 0.0,
                        low_h_at_switch - high_h_without_offset, 0.0};
  return species;
}

ThermophysicalSpec base_spec() {
  ThermophysicalSpec spec;
  spec.data_file = "thermo.d";
  spec.minimum_temperature = 200.0;
  spec.maximum_temperature = 3000.0;
  spec.temperature_relative_tolerance = 1.0e-12;
  spec.maximum_temperature_iterations = 80U;
  spec.closed_mass_relative_tolerance = 1.0e-12;
  spec.maximum_closed_mass_iterations = 20U;
  spec.maximum_closed_mass_relative_step = 0.5;
  return spec;
}

constexpr ThermalRevisionTuple thermal_revisions() noexcept {
  return {11U, 12U, 21U, 22U, 31U, 32U, 41U};
}

double species_cp_oracle(const SpeciesThermophysicalSpec& species,
                         const std::array<double, 7U>& coefficients,
                         double temperature) {
  const double cp_over_r =
      ((((coefficients[4U] * temperature + coefficients[3U]) * temperature +
         coefficients[2U]) *
            temperature +
        coefficients[1U]) *
           temperature +
       coefficients[0U]);
  return kUniversalGasConstant * cp_over_r / species.molecular_weight;
}

double species_h_oracle(const SpeciesThermophysicalSpec& species,
                        const std::array<double, 7U>& coefficients,
                        double temperature) {
  const double t2 = temperature * temperature;
  const double t3 = t2 * temperature;
  const double t4 = t3 * temperature;
  const double h_over_r =
      coefficients[0U] * temperature + coefficients[1U] * t2 * 0.5 +
      coefficients[2U] * t3 / 3.0 + coefficients[3U] * t4 * 0.25 +
      coefficients[4U] * t4 * temperature * 0.2 + coefficients[5U];
  return kUniversalGasConstant * h_over_r / species.molecular_weight;
}

bool test_constant_cp_path() {
  ThermophysicalSpec spec = base_spec();
  spec.species.push_back(constant_species("A", 20.0, 1000.0));
  spec.species.push_back(constant_species("B", 40.0, 1200.0));
  std::array<TransportedScalarSpec, 2U> catalog{
      TransportedScalarSpec{"A", TransportedScalarRole::species},
      TransportedScalarSpec{"tracer",
                            TransportedScalarRole::passive_scalar}};
  ThermodynamicsPlan plan;
  bool passed = expect(static_cast<bool>(ThermodynamicsPlan::compile(
                           spec, Span<const TransportedScalarSpec>{
                                     catalog.data(), catalog.size()},
                           plan)),
                       "constant-cp thermo plan compiles");
  passed &= expect(plan.kernel() == ThermodynamicsKernel::constant_cp &&
                       plan.species_count() == 2U &&
                       plan.independent_species_count() == 1U &&
                       plan.dependent_species_index() == 1U &&
                       plan.fingerprint() != 0U,
                   "constant-cp plan publishes deterministic N-1 metadata");

  const std::array<double, 1U> y{0.25};
  double h = 0.0;
  double cp = 0.0;
  double gas_constant = 0.0;
  passed &= expect(static_cast<bool>(plan.mixture_enthalpy(
                       600.0, Span<const double>{y.data(), y.size()}, h, cp,
                       gas_constant)),
                   "mixture enthalpy evaluates");
  constexpr double expected_cp = 1150.0;
  const double expected_r =
      0.25 * kUniversalGasConstant / 20.0 +
      0.75 * kUniversalGasConstant / 40.0;
  passed &= expect(near(cp, expected_cp, 1.0e-13) &&
                       near(h, expected_cp * 600.0, 1.0e-13) &&
                       near(gas_constant, expected_r, 1.0e-13),
                   "constant-cp mixture uses mass-fraction weighting");

  ThermoState state;
  const double speed = 12.5;
  const double p = 123456.0;
  passed &= expect(static_cast<bool>(plan.evaluate(
                       p, h, Span<const double>{y.data(), y.size()},
                       Real3{speed, 0.0, 0.0}, state)),
                   "constant-cp state evaluates from h");
  const double gamma = expected_cp / (expected_cp - expected_r);
  const double sound = std::sqrt(gamma * expected_r * 600.0);
  passed &= expect(near(state.temperature, 600.0, 1.0e-11) &&
                       near(state.rho, p / (expected_r * 600.0), 1.0e-12) &&
                       near(state.cp, expected_cp, 1.0e-12) &&
                       near(state.gas_constant, expected_r, 1.0e-12) &&
                       near(state.gamma, gamma, 1.0e-12) &&
                       near(state.drho_dp_hY,
                            1.0 / (expected_r * 600.0), 1.0e-12) &&
                       near(state.sound_speed, sound, 1.0e-12) &&
                       near(state.mach, speed / sound, 1.0e-12),
                   "EOS, analytic drho/dp, sound speed, and Mach are consistent");

  ThermoState referenced;
  passed &= expect(static_cast<bool>(plan.evaluate_from_reference_pressure(
                       120000.0, 3456.0, h,
                       Span<const double>{y.data(), y.size()},
                       Real3{speed, 0.0, 0.0}, referenced)) &&
                       near(referenced.rho, state.rho, 1.0e-14),
                   "p_abs is exactly p_ref plus pi");

  const double dp = 0.25;
  ThermoState plus;
  ThermoState minus;
  passed &= expect(static_cast<bool>(plan.evaluate(
                       p + dp, h, Span<const double>{y.data(), y.size()}, {},
                       plus)) &&
                       static_cast<bool>(plan.evaluate(
                           p - dp, h,
                           Span<const double>{y.data(), y.size()}, {}, minus)),
                   "finite-difference derivative states evaluate");
  passed &= expect(near((plus.rho - minus.rho) / (2.0 * dp),
                        state.drho_dp_hY, 1.0e-8),
                   "analytic drho/dp matches finite difference at fixed h,Y");

  const double dh = 0.05;
  passed &= expect(static_cast<bool>(plan.evaluate(
                       p, h + dh, Span<const double>{y.data(), y.size()}, {},
                       plus)) &&
                       static_cast<bool>(plan.evaluate(
                           p, h - dh,
                           Span<const double>{y.data(), y.size()}, {}, minus)),
                   "constant-cp enthalpy derivative states evaluate");
  passed &= expect(near((plus.rho - minus.rho) / (2.0 * dh),
                        state.drho_dh_pY, 1.0e-7) &&
                       state.drho_dh_pY < 0.0,
                   "analytic drho/dh matches finite difference at fixed p,Y");

  ThermophysicalSpec generic_spec = spec;
  constexpr double roundoff_slope = 1.0e-18;
  for (SpeciesThermophysicalSpec& species : generic_spec.species) {
    species.nasa7_low[1U] = roundoff_slope;
    species.nasa7_high[1U] = roundoff_slope;
  }
  ThermodynamicsPlan generic;
  passed &= expect(static_cast<bool>(ThermodynamicsPlan::compile(
                       generic_spec,
                       Span<const TransportedScalarSpec>{catalog.data(),
                                                         catalog.size()},
                       generic)) &&
                       generic.kernel() == ThermodynamicsKernel::nasa7,
                   "a roundoff-close polynomial selects the generic kernel");
  double generic_h = 0.0;
  double generic_cp = 0.0;
  double generic_r = 0.0;
  ThermoState generic_state;
  passed &= expect(
      static_cast<bool>(generic.mixture_enthalpy(
          600.0, Span<const double>{y.data(), y.size()}, generic_h,
          generic_cp, generic_r)) &&
          static_cast<bool>(generic.evaluate(
              p, generic_h, Span<const double>{y.data(), y.size()},
              Real3{speed, 0.0, 0.0}, generic_state)) &&
          near(generic_state.temperature, state.temperature, 1.0e-13) &&
          near(generic_state.rho, state.rho, 1.0e-13) &&
          near(generic_state.cp, state.cp, 1.0e-13) &&
          near(generic_state.gas_constant, state.gas_constant, 1.0e-13) &&
          near(generic_state.sound_speed, state.sound_speed, 1.0e-13),
      "constant-cp specialization agrees with a generic roundoff oracle");

  Status repeated_status;
  ThermoState repeated;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    for (std::size_t iteration = 0U; iteration < 1000U; ++iteration) {
      repeated_status = plan.evaluate(
          p, h, Span<const double>{y.data(), y.size()},
          Real3{speed, 0.0, 0.0}, repeated, 600.0);
    }
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(static_cast<bool>(repeated_status) &&
                       hot_allocations == 0U &&
                       repeated.temperature == state.temperature &&
                       repeated.rho == state.rho,
                   "constant-cp hot evaluation is allocation-free");
  return passed;
}

bool test_conserved_enthalpy_bounds() {
  bool passed = true;

  ThermophysicalSpec single_spec = base_spec();
  single_spec.maximum_temperature = 3200.0;
  single_spec.species.push_back(constant_species("air", 28.96546, 1005.0));
  ThermodynamicsPlan single;
  const Span<const TransportedScalarSpec> empty_catalog{nullptr, 0U};
  passed &= expect(static_cast<bool>(ThermodynamicsPlan::compile(
                           single_spec, empty_catalog, single)),
                   "single-species thermodynamics compiles for bounds");
  const double single_density = 1.7;
  double single_lower = 0.0;
  double single_upper = 0.0;
  passed &= expect(static_cast<bool>(single.species_enthalpy_bounds(
                       0U, single_lower, single_upper)),
                   "single-species endpoint cache is published");
  for (const double endpoint : {single_spec.minimum_temperature,
                                single_spec.maximum_temperature}) {
    double mixture_h = 0.0;
    double mixture_cp = 0.0;
    double mixture_r = 0.0;
    passed &= expect(static_cast<bool>(single.mixture_enthalpy(
                         endpoint, Span<const double>{nullptr, 0U}, mixture_h,
                         mixture_cp, mixture_r)),
                     "single-species endpoint mixture enthalpy evaluates");
    const double conserved = single_density * mixture_h;
    passed &= expect(
        near(single_density *
                 (endpoint == single_spec.minimum_temperature ? single_lower
                                                               : single_upper),
             conserved, 2.0e-14),
        "single-species bound equals rho times endpoint enthalpy");
  }

  const PlanFingerprint single_fingerprint = single.fingerprint();
  ThermodynamicsPlan single_moved = std::move(single);
  passed &= expect(single_moved.fingerprint() == single_fingerprint,
                   "moving a thermodynamics plan retains bound identity");
  ThermodynamicsPlan single_recompiled;
  passed &= expect(static_cast<bool>(ThermodynamicsPlan::compile(
                           single_spec, empty_catalog, single_recompiled)) &&
                       single_recompiled.fingerprint() == single_fingerprint,
                   "cached endpoint bounds preserve deterministic fingerprint");

  ThermophysicalSpec multi_spec = base_spec();
  multi_spec.maximum_temperature = 3200.0;
  multi_spec.species.push_back(varying_species("A", 24.0, 3.2, 4.0e-4));
  multi_spec.species.push_back(varying_species("B", 32.0, 3.8, 2.0e-4));
  multi_spec.species.push_back(varying_species("C", 40.0, 4.1, 1.0e-4));
  std::array<TransportedScalarSpec, 2U> multi_catalog{
      TransportedScalarSpec{"A", TransportedScalarRole::species},
      TransportedScalarSpec{"B", TransportedScalarRole::species}};
  ThermodynamicsPlan multi;
  passed &= expect(static_cast<bool>(ThermodynamicsPlan::compile(
                           multi_spec,
                           Span<const TransportedScalarSpec>{
                               multi_catalog.data(), multi_catalog.size()},
                           multi)),
                   "multi-species thermodynamics compiles for bounds");

  std::array<TransportedScalarSpec, 2U> reordered_catalog{
      TransportedScalarSpec{"B", TransportedScalarRole::species},
      TransportedScalarSpec{"A", TransportedScalarRole::species}};
  ThermodynamicsPlan reordered;
  passed &= expect(static_cast<bool>(ThermodynamicsPlan::compile(
                           multi_spec,
                           Span<const TransportedScalarSpec>{
                               reordered_catalog.data(),
                               reordered_catalog.size()},
                           reordered)),
                   "reordered species catalog compiles for bounds");
  double mapped_minimum = 0.0;
  double mapped_maximum = 0.0;
  double direct_minimum = 0.0;
  double direct_maximum = 0.0;
  passed &= expect(
      static_cast<bool>(reordered.independent_species_enthalpy_bounds(
          0U, mapped_minimum, mapped_maximum)) &&
          static_cast<bool>(reordered.species_enthalpy_bounds(
              1U, direct_minimum, direct_maximum)) &&
          mapped_minimum == direct_minimum && mapped_maximum == direct_maximum,
      "independent endpoint lookup follows the catalog-to-species mapping");

  // This vector is deliberately outside the composition simplex.  The API
  // is a linear conserved-state operation and must still produce both bounds.
  const double non_simplex_density = 1.0;
  const std::array<double, 2U> non_simplex_species{1.2, -0.4};
  const double non_simplex_dependent =
      non_simplex_density - non_simplex_species[0U] -
      non_simplex_species[1U];
  const std::array<double, 3U> non_simplex_full_species{
      non_simplex_species[0U], non_simplex_species[1U],
      non_simplex_dependent};
  for (const double endpoint : {multi_spec.minimum_temperature,
                                multi_spec.maximum_temperature}) {
    double expected = 0.0;
    double actual = 0.0;
    for (std::size_t species_index = 0U; species_index < 3U;
         ++species_index) {
      double at_minimum = 0.0;
      double at_maximum = 0.0;
      passed &= expect(static_cast<bool>(multi.species_enthalpy_bounds(
                           species_index, at_minimum, at_maximum)),
                       "multi-species endpoint cache covers every species");
      actual += non_simplex_full_species[species_index] *
                (endpoint == multi_spec.minimum_temperature ? at_minimum
                                                             : at_maximum);
      const SpeciesThermophysicalSpec& species =
          multi_spec.species[species_index];
      const auto& coefficients =
          endpoint <= species.temperature_switch ? species.nasa7_low
                                                 : species.nasa7_high;
      expected += non_simplex_full_species[species_index] *
                  species_h_oracle(species, coefficients, endpoint);
    }
    passed &= expect(near(actual, expected, 3.0e-14),
                     "multi-species bound is the complete linear species sum");
  }

  const double physical_density = 2.4;
  const std::array<double, 2U> physical_species{0.72, 0.48};
  const std::array<double, 2U> physical_fractions{
      physical_species[0U] / physical_density,
      physical_species[1U] / physical_density};
  const std::array<double, 3U> physical_full_species{
      physical_species[0U], physical_species[1U],
      physical_density - physical_species[0U] - physical_species[1U]};
  for (const double endpoint : {multi_spec.minimum_temperature,
                                multi_spec.maximum_temperature}) {
    double mixture_h = 0.0;
    double mixture_cp = 0.0;
    double mixture_r = 0.0;
    passed &= expect(static_cast<bool>(multi.mixture_enthalpy(
                         endpoint,
                         Span<const double>{physical_fractions.data(),
                                             physical_fractions.size()},
                         mixture_h, mixture_cp, mixture_r)),
                     "multi-species endpoint mixture enthalpy evaluates");
    double actual = 0.0;
    for (std::size_t species_index = 0U; species_index < 3U;
         ++species_index) {
      double at_minimum = 0.0;
      double at_maximum = 0.0;
      passed &= expect(static_cast<bool>(multi.species_enthalpy_bounds(
                           species_index, at_minimum, at_maximum)),
                       "physical endpoint cache lookup succeeds");
      actual += physical_full_species[species_index] *
                (endpoint == multi_spec.minimum_temperature ? at_minimum
                                                             : at_maximum);
    }
    passed &= expect(near(actual, physical_density * mixture_h, 3.0e-14),
                     "multi-species bound equals rho times mixture enthalpy");
  }

  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    for (std::size_t iteration = 0U; iteration < 1000U; ++iteration) {
      for (std::size_t species_index = 0U; species_index < 3U;
           ++species_index) {
        double lower = 0.0;
        double upper = 0.0;
        (void)multi.species_enthalpy_bounds(species_index, lower, upper);
      }
    }
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(hot_allocations == 0U,
                   "conserved endpoint bounds are allocation-free");

  double untouched_lower = 17.0;
  double untouched_upper = 19.0;
  passed &= expect(!multi.species_enthalpy_bounds(
                       multi.species_count(), untouched_lower, untouched_upper) &&
                       untouched_lower == 17.0 && untouched_upper == 19.0,
                   "out-of-range endpoint lookup is rejected without publication");
  return passed;
}

bool test_nasa7_inversion_and_validation() {
  ThermophysicalSpec spec = base_spec();
  spec.species.push_back(varying_species("A", 24.0, 3.2, 4.0e-4));
  spec.species.push_back(varying_species("B", 32.0, 3.8, 2.0e-4));
  std::array<TransportedScalarSpec, 1U> catalog{
      TransportedScalarSpec{"A", TransportedScalarRole::species}};
  ThermodynamicsPlan plan;
  bool passed = expect(static_cast<bool>(ThermodynamicsPlan::compile(
                           spec, Span<const TransportedScalarSpec>{
                                     catalog.data(), catalog.size()},
                           plan)) &&
                           plan.kernel() == ThermodynamicsKernel::nasa7,
                       "temperature-dependent NASA7 plan compiles");
  const std::array<double, 1U> y{0.35};
  double h = 0.0;
  double cp = 0.0;
  double r = 0.0;
  passed &= expect(static_cast<bool>(plan.mixture_enthalpy(
                       1375.0, Span<const double>{y.data(), y.size()}, h, cp,
                       r)),
                   "NASA7 forward enthalpy evaluates");
  const double expected_high_h =
      y[0U] * species_h_oracle(spec.species[0U],
                               spec.species[0U].nasa7_high, 1375.0) +
      (1.0 - y[0U]) *
          species_h_oracle(spec.species[1U], spec.species[1U].nasa7_high,
                           1375.0);
  const double expected_high_cp =
      y[0U] * species_cp_oracle(spec.species[0U],
                                spec.species[0U].nasa7_high, 1375.0) +
      (1.0 - y[0U]) *
          species_cp_oracle(spec.species[1U], spec.species[1U].nasa7_high,
                            1375.0);
  const double wrong_low_cp =
      y[0U] * species_cp_oracle(spec.species[0U],
                                spec.species[0U].nasa7_low, 1375.0) +
      (1.0 - y[0U]) *
          species_cp_oracle(spec.species[1U], spec.species[1U].nasa7_low,
                            1375.0);
  passed &= expect(near(h, expected_high_h, 2.0e-13) &&
                       near(cp, expected_high_cp, 2.0e-13) &&
                       !near(cp, wrong_low_cp, 1.0e-6),
                   "temperature above the switch uses the high NASA interval");
  ThermoState state;
  passed &= expect(static_cast<bool>(plan.evaluate(
                       101325.0, h,
                       Span<const double>{y.data(), y.size()},
                       Real3{3.0, 4.0, 12.0}, state, 210.0)),
                   "safeguarded inversion converges from a poor hint");
  passed &= expect(near(state.temperature, 1375.0, 2.0e-11) &&
                       near(state.cp, cp, 2.0e-11) &&
                       near(state.mach, 13.0 / state.sound_speed, 1.0e-13),
                   "NASA7 inversion returns the forward state and |U|/a Mach");

  const double dh = std::max(1.0e-3, std::abs(h) * 2.0e-7);
  ThermoState plus_h;
  ThermoState minus_h;
  passed &= expect(static_cast<bool>(plan.evaluate(
                       101325.0, h + dh,
                       Span<const double>{y.data(), y.size()}, {}, plus_h,
                       1375.0)) &&
                       static_cast<bool>(plan.evaluate(
                           101325.0, h - dh,
                           Span<const double>{y.data(), y.size()}, {}, minus_h,
                           1375.0)),
                   "NASA7 enthalpy derivative states evaluate");
  passed &= expect(near((plus_h.rho - minus_h.rho) / (2.0 * dh),
                        state.drho_dh_pY, 2.0e-7) &&
                       state.drho_dh_pY < 0.0,
                   "NASA7 analytic drho/dh matches finite difference");

  Status repeated_status;
  ThermoState repeated;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    for (std::size_t iteration = 0U; iteration < 1000U; ++iteration) {
      repeated_status = plan.evaluate(
          101325.0, h, Span<const double>{y.data(), y.size()},
          Real3{3.0, 4.0, 12.0}, repeated, 1375.0);
    }
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(static_cast<bool>(repeated_status) &&
                       hot_allocations == 0U &&
                       near(repeated.temperature, 1375.0, 2.0e-11),
                   "generic NASA7 hot evaluation is allocation-free");

  ThermalState thermal;
  constexpr ThermalRevisionTuple revisions = thermal_revisions();
  passed &= expect(static_cast<bool>(plan.evaluate_thermal(
                       h, Span<const double>{y.data(), y.size()}, revisions,
                       thermal, 210.0)),
                   "pressure-independent thermal state evaluates once");
  ThermoState low_pressure;
  ThermoState high_pressure;
  passed &= expect(static_cast<bool>(plan.complete_state(
                       90000.0, thermal, revisions,
                       Real3{0.0, 0.0, 0.0},
                       low_pressure)) &&
                       static_cast<bool>(plan.complete_state(
                           180000.0, thermal, revisions,
                           Real3{3.0, 4.0, 0.0},
                           high_pressure)),
                   "one thermal state completes multiple pressure/velocity states");
  passed &= expect(near(high_pressure.rho, 2.0 * low_pressure.rho,
                        1.0e-13) &&
                       low_pressure.temperature == thermal.temperature() &&
                       high_pressure.temperature == thermal.temperature() &&
                       low_pressure.mach == 0.0 &&
                       near(high_pressure.mach,
                            5.0 / high_pressure.sound_speed, 1.0e-13),
                   "pressure and velocity completion does not alter thermal state");

  const auto rejects_stale_thermal = [&](ThermalRevisionTuple stale,
                                         std::string_view description) {
    ThermoState untouched;
    untouched.rho = 77.0;
    return expect(!plan.complete_state(101325.0, thermal, stale, {}, untouched) &&
                      untouched.rho == 77.0,
                  description);
  };
  ThermalRevisionTuple stale = revisions;
  ++stale.enthalpy;
  passed &= rejects_stale_thermal(stale,
                                  "stale enthalpy revision rejects thermal reuse");
  stale = revisions;
  ++stale.composition;
  passed &= rejects_stale_thermal(
      stale, "stale composition revision rejects thermal reuse");
  stale = revisions;
  ++stale.enthalpy_storage;
  passed &= rejects_stale_thermal(
      stale, "stale enthalpy storage rejects thermal reuse");
  stale = revisions;
  ++stale.composition_storage;
  passed &= rejects_stale_thermal(
      stale, "stale composition storage rejects thermal reuse");
  stale = revisions;
  ++stale.enthalpy_revision_domain;
  passed &= rejects_stale_thermal(
      stale, "stale enthalpy domain rejects thermal reuse");
  stale = revisions;
  ++stale.composition_revision_domain;
  passed &= rejects_stale_thermal(
      stale, "stale composition domain rejects thermal reuse");
  stale = revisions;
  ++stale.patch_identity;
  passed &= rejects_stale_thermal(stale,
                                  "stale patch identity rejects thermal reuse");

  ThermophysicalSpec endpoint_spec = spec;
  endpoint_spec.maximum_temperature_iterations = 1U;
  ThermodynamicsPlan endpoint_plan;
  passed &= expect(static_cast<bool>(ThermodynamicsPlan::compile(
                       endpoint_spec, Span<const TransportedScalarSpec>{
                                          catalog.data(), catalog.size()},
                       endpoint_plan)),
                   "generic endpoint plan compiles with one allowed iteration");
  for (const double endpoint :
       {endpoint_spec.minimum_temperature, endpoint_spec.maximum_temperature}) {
    double endpoint_h = 0.0;
    double endpoint_cp = 0.0;
    double endpoint_r = 0.0;
    ThermalState endpoint_state;
    passed &= expect(static_cast<bool>(endpoint_plan.mixture_enthalpy(
                         endpoint, Span<const double>{y.data(), y.size()},
                         endpoint_h, endpoint_cp, endpoint_r)) &&
                         static_cast<bool>(endpoint_plan.evaluate_thermal(
                             endpoint_h,
                             Span<const double>{y.data(), y.size()},
                             revisions, endpoint_state)) &&
                         endpoint_state.temperature() == endpoint,
                     "exact Tmin/Tmax enthalpy bypasses iterative inversion");
  }

  ThermoState sentinel;
  sentinel.rho = 77.0;
  const std::array<double, 1U> negative{-0.01};
  const std::array<double, 1U> above_sum{1.01};
  passed &= expect(!plan.evaluate(101325.0, h,
                                  Span<const double>{negative.data(), 1U},
                                  {}, sentinel) &&
                       sentinel.rho == 77.0,
                   "negative independent species rejects without output publication");
  passed &= expect(!plan.evaluate(101325.0, h,
                                  Span<const double>{above_sum.data(), 1U},
                                  {}, sentinel),
                   "independent species sum above one is rejected");
  passed &= expect(!plan.evaluate(0.0, h,
                                  Span<const double>{y.data(), y.size()}, {},
                                  sentinel),
                   "nonpositive absolute pressure is rejected");
  passed &= expect(!plan.evaluate_from_reference_pressure(
                       std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::max(), h,
                       Span<const double>{y.data(), y.size()}, {}, sentinel),
                   "overflowing p_ref plus pi is rejected");
  passed &= expect(!plan.evaluate(
                       101325.0, std::numeric_limits<double>::infinity(),
                       Span<const double>{y.data(), y.size()}, {}, sentinel),
                   "nonfinite enthalpy is rejected");

  ThermophysicalSpec invalid = spec;
  invalid.species[1U].stable_name = "A";
  ThermodynamicsPlan rejected;
  passed &= expect(!ThermodynamicsPlan::compile(
                       invalid, Span<const TransportedScalarSpec>{
                                    catalog.data(), catalog.size()},
                       rejected),
                   "duplicate thermo species are rejected");
  invalid = spec;
  invalid.species[0U].stable_name = "wrong";
  passed &= expect(!ThermodynamicsPlan::compile(
                       invalid, Span<const TransportedScalarSpec>{
                                    catalog.data(), catalog.size()},
                       rejected),
                   "thermo species bind explicitly to catalog species names");
  ThermophysicalSpec reordered = spec;
  std::swap(reordered.species[0U], reordered.species[1U]);
  ThermodynamicsPlan reordered_plan;
  passed &= expect(static_cast<bool>(ThermodynamicsPlan::compile(
                       reordered, Span<const TransportedScalarSpec>{
                                      catalog.data(), catalog.size()},
                       reordered_plan)) &&
                       reordered_plan.dependent_species_index() == 0U,
                   "species bind by stable name independent of .d ordering");
  double reordered_h = 0.0;
  double reordered_cp = 0.0;
  double reordered_r = 0.0;
  passed &= expect(static_cast<bool>(reordered_plan.mixture_enthalpy(
                       1375.0, Span<const double>{y.data(), y.size()},
                       reordered_h, reordered_cp, reordered_r)) &&
                       near(reordered_h, h, 1.0e-13) &&
                       near(reordered_cp, cp, 1.0e-13) &&
                       near(reordered_r, r, 1.0e-13),
                   "reordered thermophysical species preserve mixture state");
  invalid = spec;
  invalid.species.pop_back();
  passed &= expect(!ThermodynamicsPlan::compile(
                       invalid, Span<const TransportedScalarSpec>{
                                    catalog.data(), catalog.size()},
                       rejected),
                   "strict N-1 requires exactly one dependent thermo species");
  return passed;
}

bool test_thermophysical_text_contract() {
  ThermophysicalSpec parsed;
  bool passed = expect(static_cast<bool>(
                           detail::parse_thermophysical_text(
                               kThermophysicalText, parsed)),
                       "versioned thermophysical .d parses strictly");
  passed &= expect(parsed.species.size() == 1U &&
                       parsed.species[0].stable_name == "air" &&
                       parsed.species[0].transport_law ==
                           TransportLaw::sutherland &&
                       parsed.maximum_temperature_iterations == 64U &&
                       parsed.maximum_closed_mass_iterations == 32U,
                   "typed .d values are published without filename guessing");
  passed &= expect(detail::valid_thermophysical_spec(parsed),
                   "parsed typed payload passes the shared validator");
  ThermophysicalSpec coast_native_air;
  passed &= expect(
      static_cast<bool>(detail::parse_thermophysical_text(
          kCoastNativeAirText, coast_native_air)) &&
          coast_native_air.species.size() == 1U &&
          coast_native_air.species.front().transport_law ==
              TransportLaw::coast_native_air &&
          coast_native_air.species.front().viscosity_reference == 0.0,
      "COAST-native air has a parameter-free strict transport record");
  passed &= expect(detail::thermophysical_spec_fingerprint(parsed) == 0U,
                   "fingerprint rejects a payload without source identity");
  parsed.data_file = "thermophysics.d";
  const PlanFingerprint fingerprint =
      detail::thermophysical_spec_fingerprint(parsed);
  passed &= expect(fingerprint != 0U &&
                       fingerprint ==
                           detail::thermophysical_spec_fingerprint(parsed),
                   "complete thermophysical fingerprint is deterministic");
  ThermophysicalSpec changed_fingerprint = parsed;
  changed_fingerprint.species[0U].prandtl = 0.72;
  passed &= expect(
      detail::thermophysical_spec_fingerprint(changed_fingerprint) !=
          fingerprint,
      "complete fingerprint covers species transport fields");
  ThermophysicalSpec sentinel;
  sentinel.minimum_temperature = 77.0;
  std::string trailing{kThermophysicalText};
  trailing += "unknown_record 1\n";
  passed &= expect(
      !detail::parse_thermophysical_text(trailing, sentinel) &&
          sentinel.minimum_temperature == 77.0,
      "unknown trailing records reject atomically");
  std::string discontinuous{kThermophysicalText};
  const std::string old_high = "nasa7_high 3.5 0 0 0 0 0 0";
  const std::size_t position = discontinuous.find(old_high);
  passed &= expect(position != std::string::npos,
                   "continuity mutation locates the high polynomial");
  if (position != std::string::npos) {
    discontinuous.replace(position, old_high.size(),
                          "nasa7_high 3.5 0 0 0 0 100 0");
    passed &= expect(!detail::parse_thermophysical_text(
                         discontinuous, sentinel),
                     "discontinuous NASA enthalpy intervals are rejected");
  }

  std::string hidden_cp_failure{kThermophysicalText};
  const std::string old_low = "nasa7_low 3.5 0 0 0 0 0 0";
  const std::string hidden_polynomial =
      "16.5 -0.08 0.0001 0 0 0 0";
  std::size_t low_position = hidden_cp_failure.find(old_low);
  std::size_t high_position = hidden_cp_failure.find(old_high);
  passed &= expect(low_position != std::string::npos &&
                       high_position != std::string::npos,
                   "interior cp mutation locates both NASA intervals");
  if (low_position != std::string::npos &&
      high_position != std::string::npos) {
    hidden_cp_failure.replace(low_position, old_low.size(),
                              "nasa7_low " + hidden_polynomial);
    high_position = hidden_cp_failure.find(old_high);
    hidden_cp_failure.replace(high_position, old_high.size(),
                              "nasa7_high " + hidden_polynomial);
    passed &= expect(!detail::parse_thermophysical_text(
                         hidden_cp_failure, sentinel),
                     "full-interval proof rejects cp<=R between old samples");
  }

  for (const std::string& endpoint_polynomial :
       {std::string{"0.21875 0.00390625 0 0 0 0 0"},
        std::string{"1.732421875 -0.000244140625 0 0 0 0 0"}}) {
    std::string endpoint_cp{kThermophysicalText};
    low_position = endpoint_cp.find(old_low);
    high_position = endpoint_cp.find(old_high);
    passed &= expect(low_position != std::string::npos &&
                         high_position != std::string::npos,
                     "cp endpoint mutation locates both NASA intervals");
    if (low_position != std::string::npos &&
        high_position != std::string::npos) {
      endpoint_cp.replace(low_position, old_low.size(),
                          "nasa7_low " + endpoint_polynomial);
      high_position = endpoint_cp.find(old_high);
      endpoint_cp.replace(high_position, old_high.size(),
                          "nasa7_high " + endpoint_polynomial);
      passed &= expect(!detail::parse_thermophysical_text(
                           endpoint_cp, sentinel),
                       "strict cp>R rejects equality at Tmin or Tmax");
    }
  }

  std::string cp_discontinuous{kThermophysicalText};
  high_position = cp_discontinuous.find(old_high);
  passed &= expect(high_position != std::string::npos,
                   "cp continuity mutation locates the high polynomial");
  if (high_position != std::string::npos) {
    cp_discontinuous.replace(high_position, old_high.size(),
                             "nasa7_high 3.6 0 0 0 0 -100 0");
    passed &= expect(!detail::parse_thermophysical_text(
                         cp_discontinuous, sentinel),
                     "cp jump rejects even when switch enthalpy is continuous");
  }

  {
    const std::locale comma_locale{
        std::locale::classic(), new CommaDecimalPoint};
    const ScopedGlobalLocale scoped_locale(comma_locale);
    ThermophysicalSpec locale_parsed;
    passed &= expect(static_cast<bool>(
                         detail::parse_thermophysical_text(
                             kThermophysicalText, locale_parsed)),
                     "dot-decimal input is independent of the global locale");
    std::string comma_decimal{kThermophysicalText};
    const std::string dot_bounds = "temperature_bounds 200 3000";
    const std::size_t bounds_position = comma_decimal.find(dot_bounds);
    passed &= expect(bounds_position != std::string::npos,
                     "locale mutation locates temperature bounds");
    if (bounds_position != std::string::npos) {
      comma_decimal.replace(bounds_position, dot_bounds.size(),
                            "temperature_bounds 200,0 3000");
      sentinel.minimum_temperature = 77.0;
      passed &= expect(!detail::parse_thermophysical_text(
                           comma_decimal, sentinel) &&
                           sentinel.minimum_temperature == 77.0,
                       "comma decimal is rejected atomically in every locale");
    }
  }

  for (const std::pair<std::string, std::string>& mutation :
       {std::pair<std::string, std::string>{
            "temperature_inversion 1e-12 64",
            "temperature_inversion 1e-12 257"},
        std::pair<std::string, std::string>{
            "closed_mass_newton 1e-12 32 0.2",
            "closed_mass_newton 1e-12 257 0.2"},
        std::pair<std::string, std::string>{"temperature_switch 1000",
                                            "temperature_switch 200"},
        std::pair<std::string, std::string>{"temperature_switch 1000",
                                            "temperature_switch 3000"}}) {
    std::string invalid_control{kThermophysicalText};
    const std::size_t control_position = invalid_control.find(mutation.first);
    passed &= expect(control_position != std::string::npos,
                     "bounded-control mutation locates its input record");
    if (control_position != std::string::npos) {
      invalid_control.replace(control_position, mutation.first.size(),
                              mutation.second);
      passed &= expect(!detail::parse_thermophysical_text(
                           invalid_control, sentinel),
                       "iteration overflow or degenerate NASA interval rejects");
    }
  }
  return passed;
}

}  // namespace

int main() {
  bool passed = test_constant_cp_path();
  passed &= test_conserved_enthalpy_bounds();
  passed &= test_nasa7_inversion_and_validation();
  passed &= test_thermophysical_text_contract();
  if (passed) {
    std::cout << "v0.4 thermodynamics tests passed\n";
  }
  return passed ? 0 : 1;
}
