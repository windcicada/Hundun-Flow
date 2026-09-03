// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/bc_basic_boundary.hpp"

#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_kernel_field_view.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::boundary {
namespace detail {
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
namespace {

thread_local int next_outlet_topology_observation = 0;

} // namespace
#endif

std::string_view fixed_preflight_message(const runtime::Error &) noexcept {
  return "final pressure-outlet preflight rejected local data";
}

std::string_view fixed_preflight_message(const std::exception &) noexcept {
  return "final pressure-outlet preflight failed locally";
}

std::string_view fixed_unknown_preflight_message() noexcept {
  return "final pressure-outlet preflight failed unexpectedly";
}

#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
void set_next_outlet_topology_observation_raw(int observation) {
  if (observation == 0) {
    throw runtime::Error("outlet topology test observation must change data");
  }
  if (next_outlet_topology_observation != 0) {
    throw runtime::Error("outlet topology test observation is already pending");
  }
  next_outlet_topology_observation = observation;
}

int consume_next_outlet_topology_observation() noexcept {
  return std::exchange(next_outlet_topology_observation, 0);
}
#endif

} // namespace detail

namespace {

constexpr std::array<std::string_view, 6> kPatchNames{
    "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};
constexpr double kLimbBase = 4294967296.0;

template <class Enum> bool enum_in_range(Enum value, Enum last) {
  using Underlying = std::underlying_type_t<Enum>;
  const auto converted = static_cast<Underlying>(value);
  return converted >= Underlying{} &&
         converted <= static_cast<Underlying>(last);
}

bool finite(double value) { return std::isfinite(value); }

void require_finite(double value, const char *message) {
  if (!finite(value)) {
    throw runtime::Error(message);
  }
}

void require_positive(double value, const char *message) {
  if (!finite(value) || value <= 0.0) {
    throw runtime::Error(message);
  }
}

void require_finite(runtime::Real3 value, const char *message) {
  if (!finite(value.x) || !finite(value.y) || !finite(value.z)) {
    throw runtime::Error(message);
  }
}

bool valid_scalar_name(std::string_view name) {
  if (name.empty()) {
    return false;
  }
  const auto alpha_or_underscore = [](char character) {
    return (character >= 'A' && character <= 'Z') ||
           (character >= 'a' && character <= 'z') || character == '_';
  };
  const auto alnum_or_underscore = [&](char character) {
    return alpha_or_underscore(character) ||
           (character >= '0' && character <= '9');
  };
  return alpha_or_underscore(name.front()) &&
         std::all_of(name.begin() + 1, name.end(), alnum_or_underscore);
}

bool relatively_equal(double expected, double actual, double tolerance) {
  const double denominator = std::max({std::abs(expected), std::abs(actual),
                                       std::numeric_limits<double>::min()});
  return std::abs(expected - actual) / denominator <= tolerance;
}

double scaled_positive_ratio(double numerator, double denominator_a,
                             double denominator_b) {
  int numerator_exponent = 0;
  int denominator_a_exponent = 0;
  int denominator_b_exponent = 0;
  const double numerator_fraction = std::frexp(numerator, &numerator_exponent);
  const double denominator_a_fraction =
      std::frexp(denominator_a, &denominator_a_exponent);
  const double denominator_b_fraction =
      std::frexp(denominator_b, &denominator_b_exponent);
  const double fraction =
      (numerator_fraction / denominator_a_fraction) / denominator_b_fraction;
  const std::int64_t exponent =
      static_cast<std::int64_t>(numerator_exponent) -
      static_cast<std::int64_t>(denominator_a_exponent) -
      static_cast<std::int64_t>(denominator_b_exponent);
  if (exponent > std::numeric_limits<int>::max()) {
    return std::numeric_limits<double>::infinity();
  }
  if (exponent < std::numeric_limits<int>::min()) {
    return 0.0;
  }
  return std::scalbn(fraction, static_cast<int>(exponent));
}

double checked_long_double(long double value, const char *message) {
  const double converted = static_cast<double>(value);
  if (!std::isfinite(converted)) {
    throw runtime::Error(message);
  }
  return converted;
}

struct Rules final {
  BoundaryKind kind;
  VelocityRule velocity;
  PressureRule pressure;
  TransportRule density;
  TransportRule transport;
  MassFluxRule mass_flux;
};

Rules rules_for(config::BoundaryType type) {
  switch (type) {
  case config::BoundaryType::periodic:
    return {BoundaryKind::periodic,       VelocityRule::periodic_pair,
            PressureRule::periodic_pair,  TransportRule::periodic_pair,
            TransportRule::periodic_pair, MassFluxRule::periodic_pair};
  case config::BoundaryType::no_slip_wall:
    return {BoundaryKind::no_slip_wall,
            VelocityRule::prescribed_zero,
            PressureRule::zero_normal_gradient,
            TransportRule::copy_interior,
            TransportRule::zero_normal_diffusive_flux,
            MassFluxRule::identically_zero};
  case config::BoundaryType::symmetry:
    return {BoundaryKind::symmetry,
            VelocityRule::reflect_normal_copy_tangential,
            PressureRule::zero_normal_gradient,
            TransportRule::copy_interior,
            TransportRule::zero_normal_diffusive_flux,
            MassFluxRule::identically_zero};
  case config::BoundaryType::velocity_inlet:
    return {BoundaryKind::velocity_inlet,
            VelocityRule::prescribed_inlet,
            PressureRule::zero_normal_gradient,
            TransportRule::prescribed_value,
            TransportRule::prescribed_value,
            MassFluxRule::prescribed_inlet_state};
  case config::BoundaryType::pressure_outlet:
    return {BoundaryKind::pressure_outlet,  VelocityRule::pure_outflow,
            PressureRule::prescribed_value, TransportRule::pure_outflow,
            TransportRule::pure_outflow,    MassFluxRule::outflow_only};
  default:
    throw runtime::Error("boundary kind is invalid");
  }
}

bool has_any_prescription(const config::FlowBoundaryConfig &boundary) {
  return boundary.velocity_m_per_s.has_value() ||
         boundary.thermal_authority.has_value() ||
         boundary.temperature_K.has_value() ||
         boundary.enthalpy_J_per_kg.has_value() ||
         boundary.density_kg_per_m3.has_value() ||
         boundary.scalar_values.has_value() ||
         boundary.pressure_perturbation_pa.has_value();
}

std::vector<double>
ordered_scalar_values(const config::FlowBoundaryConfig &boundary,
                      const std::vector<std::string> &scalar_names) {
  if (!boundary.scalar_values.has_value()) {
    throw runtime::Error("velocity inlet scalar prescriptions are missing");
  }
  if (boundary.scalar_values->size() != scalar_names.size()) {
    throw runtime::Error("velocity inlet must prescribe exactly every scalar");
  }
  std::vector<double> values(scalar_names.size());
  std::vector<bool> seen(scalar_names.size(), false);
  for (const auto &prescription : *boundary.scalar_values) {
    const auto found = std::lower_bound(scalar_names.begin(),
                                        scalar_names.end(), prescription.name);
    if (found == scalar_names.end() || *found != prescription.name) {
      throw runtime::Error("velocity inlet references an unknown scalar");
    }
    const auto index = static_cast<std::size_t>(found - scalar_names.begin());
    if (seen[index]) {
      throw runtime::Error("velocity inlet scalar prescription is duplicated");
    }
    require_finite(prescription.value,
                   "velocity inlet scalar prescription is not finite");
    seen[index] = true;
    values[index] = prescription.value;
  }
  if (!std::all_of(seen.begin(), seen.end(),
                   [](bool value) { return value; })) {
    throw runtime::Error("velocity inlet scalar prescription is incomplete");
  }
  return values;
}

ResolvedInletState
materialize_inlet(const config::FlowCaseConfig &resolved,
                  const config::FlowBoundaryConfig &boundary,
                  const std::vector<std::string> &scalar_names) {
  if (!boundary.velocity_m_per_s.has_value()) {
    throw runtime::Error("velocity inlet prescription is missing velocity");
  }
  require_finite(*boundary.velocity_m_per_s,
                 "velocity inlet prescription contains non-finite velocity");
  if (!boundary.thermal_authority.has_value() ||
      !enum_in_range(*boundary.thermal_authority,
                     config::InletThermalAuthority::enthalpy)) {
    throw runtime::Error("velocity inlet thermal authority is invalid");
  }
  if (boundary.pressure_perturbation_pa.has_value()) {
    throw runtime::Error(
        "velocity inlet must not prescribe mechanical pressure");
  }

  ResolvedInletState state{};
  state.velocity_m_per_s = *boundary.velocity_m_per_s;
  state.scalar_values = ordered_scalar_values(boundary, scalar_names);
  const double tolerance = resolved.physics.inlet_consistency_rtol;
  require_positive(tolerance,
                   "inlet consistency tolerance must be positive and finite");

  switch (resolved.density_model) {
  case config::DensityModel::constant: {
    if (*boundary.thermal_authority !=
            config::InletThermalAuthority::enthalpy ||
        !boundary.enthalpy_J_per_kg.has_value()) {
      throw runtime::Error(
          "constant-density inlet requires authoritative enthalpy");
    }
    require_finite(*boundary.enthalpy_J_per_kg,
                   "velocity inlet enthalpy is not finite");
    if (boundary.temperature_K.has_value()) {
      throw runtime::Error(
          "constant-density inlet must not prescribe temperature");
    }
    if (boundary.density_kg_per_m3.has_value()) {
      require_positive(*boundary.density_kg_per_m3,
                       "velocity inlet density must be positive and finite");
      if (!relatively_equal(resolved.physics.rho_ref_kg_per_m3,
                            *boundary.density_kg_per_m3, tolerance)) {
        throw runtime::Error(
            "velocity inlet density is inconsistent with rho_ref");
      }
    }
    state.density_kg_per_m3 = resolved.physics.rho_ref_kg_per_m3;
    state.enthalpy_J_per_kg = *boundary.enthalpy_J_per_kg;
    break;
  }
  case config::DensityModel::material: {
    if (*boundary.thermal_authority !=
            config::InletThermalAuthority::enthalpy ||
        !boundary.enthalpy_J_per_kg.has_value() ||
        !boundary.density_kg_per_m3.has_value()) {
      throw runtime::Error(
          "material-density inlet requires enthalpy and density");
    }
    require_finite(*boundary.enthalpy_J_per_kg,
                   "velocity inlet enthalpy is not finite");
    require_positive(*boundary.density_kg_per_m3,
                     "velocity inlet density must be positive and finite");
    if (boundary.temperature_K.has_value()) {
      throw runtime::Error(
          "material-density inlet must not prescribe temperature");
    }
    state.density_kg_per_m3 = *boundary.density_kg_per_m3;
    state.enthalpy_J_per_kg = *boundary.enthalpy_J_per_kg;
    break;
  }
  case config::DensityModel::ideal_gas: {
    if (!resolved.physics.cp_J_per_kg_K.has_value() ||
        !resolved.physics.gas_constant_J_per_kg_K.has_value() ||
        !resolved.physics.thermodynamic_pressure_pa.has_value()) {
      throw runtime::Error("ideal-gas inlet physics is incomplete");
    }
    const double cp = *resolved.physics.cp_J_per_kg_K;
    const double gas_constant = *resolved.physics.gas_constant_J_per_kg_K;
    const double thermodynamic_pressure =
        *resolved.physics.thermodynamic_pressure_pa;
    require_positive(cp, "ideal-gas cp must be positive and finite");
    require_positive(gas_constant,
                     "ideal-gas constant must be positive and finite");
    require_positive(thermodynamic_pressure,
                     "thermodynamic pressure must be positive and finite");

    double temperature = 0.0;
    double enthalpy = 0.0;
    if (*boundary.thermal_authority ==
        config::InletThermalAuthority::temperature) {
      if (!boundary.temperature_K.has_value()) {
        throw runtime::Error(
            "ideal-gas inlet authoritative temperature is missing");
      }
      require_positive(*boundary.temperature_K,
                       "ideal-gas inlet temperature must be positive");
      temperature = *boundary.temperature_K;
      enthalpy =
          checked_long_double(static_cast<long double>(cp) * temperature,
                              "ideal-gas inlet derived enthalpy is not finite");
      if (enthalpy <= 0.0) {
        throw runtime::Error(
            "ideal-gas inlet derived enthalpy is not positive");
      }
    } else {
      if (!boundary.enthalpy_J_per_kg.has_value()) {
        throw runtime::Error(
            "ideal-gas inlet authoritative enthalpy is missing");
      }
      require_positive(*boundary.enthalpy_J_per_kg,
                       "ideal-gas inlet enthalpy must be positive");
      enthalpy = *boundary.enthalpy_J_per_kg;
      temperature = enthalpy / cp;
      require_positive(
          temperature,
          "ideal-gas inlet derived temperature is not positive and finite");
    }
    const double density = scaled_positive_ratio(thermodynamic_pressure,
                                                 gas_constant, temperature);
    require_positive(
        density, "ideal-gas inlet derived density is not positive and finite");
    if (boundary.temperature_K.has_value() &&
        !relatively_equal(temperature, *boundary.temperature_K, tolerance)) {
      throw runtime::Error(
          "ideal-gas inlet temperature conflicts with thermal authority");
    }
    if (boundary.enthalpy_J_per_kg.has_value() &&
        !relatively_equal(enthalpy, *boundary.enthalpy_J_per_kg, tolerance)) {
      throw runtime::Error(
          "ideal-gas inlet enthalpy conflicts with thermal authority");
    }
    if (boundary.density_kg_per_m3.has_value()) {
      require_positive(*boundary.density_kg_per_m3,
                       "ideal-gas inlet density must be positive");
      if (!relatively_equal(density, *boundary.density_kg_per_m3, tolerance)) {
        throw runtime::Error("ideal-gas inlet density conflicts with closure");
      }
    }
    state.temperature_K = temperature;
    state.enthalpy_J_per_kg = enthalpy;
    state.density_kg_per_m3 = density;
    break;
  }
  default:
    throw runtime::Error("density model is invalid");
  }
  return state;
}

bool same(runtime::Int3 left, runtime::Int3 right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(runtime::Box3 left, runtime::Box3 right) {
  return same(left.begin, right.begin) && same(left.end, right.end);
}

struct OutletFace final {
  mesh::LocalFaceId local_id;
  mesh::GlobalFaceId global_id;
};

struct TopologySignature final {
  runtime::Int3 global_extent{};
  runtime::Box3 owned_box{};
  std::uint64_t global_face_count{};
  std::size_t local_face_count{};
  std::array<std::optional<std::uint32_t>, 6> paired_patch_ids{};
  std::array<std::size_t, 6> patch_face_counts{};
};

TopologySignature topology_signature(const mesh::MeshTopology &topology) {
  TopologySignature signature{};
  signature.global_extent = topology.global_extent();
  signature.owned_box = topology.owned_global_box();
  signature.global_face_count = topology.global_face_count();
  signature.local_face_count = topology.local_face_count();
  for (std::size_t id = 0; id < kPatchNames.size(); ++id) {
    const auto &patch = topology.patch(static_cast<std::uint32_t>(id));
    if (patch.stable_id() != id || patch.name() != kPatchNames[id]) {
      throw runtime::Error("mesh topology patch identity is incompatible");
    }
    signature.paired_patch_ids[id] = patch.paired_patch_id();
    signature.patch_face_counts[id] = patch.local_faces().size();
  }
  return signature;
}

bool compatible(const TopologySignature &expected,
                const mesh::MeshTopology &topology,
                const std::vector<OutletFace> &expected_owned_outlet_faces,
                std::uint32_t outlet_id
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
                ,
                int observation) {
#else
                ) {
#endif
  if (!same(expected.global_extent, topology.global_extent()) ||
      !same(expected.owned_box, topology.owned_global_box()) ||
      expected.global_face_count != topology.global_face_count() ||
      expected.local_face_count != topology.local_face_count()) {
    return false;
  }
  for (std::size_t id = 0; id < kPatchNames.size(); ++id) {
    const auto &patch = topology.patch(static_cast<std::uint32_t>(id));
    if (patch.stable_id() != id || patch.name() != kPatchNames[id] ||
        patch.paired_patch_id() != expected.paired_patch_ids[id] ||
        patch.local_faces().size() != expected.patch_face_counts[id]) {
      return false;
    }
  }

  const auto &outlet_patch = topology.patch(outlet_id);
  std::size_t current_owned_outlet_face_count = 0U;
  for (const auto local_face : outlet_patch.local_faces()) {
    if (local_face >= topology.local_face_count()) {
      return false;
    }
    if (topology.face_ownership(local_face) == mesh::EntityOwnership::owned) {
      ++current_owned_outlet_face_count;
    }
  }
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
  if (observation == 6) {
    if (current_owned_outlet_face_count ==
        std::numeric_limits<std::size_t>::max()) {
      return false;
    }
    ++current_owned_outlet_face_count;
  }
#endif
  if (current_owned_outlet_face_count != expected_owned_outlet_faces.size()) {
    return false;
  }

#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
  bool first = true;
#endif
  for (const auto &expected_face : expected_owned_outlet_faces) {
    const auto local_id = expected_face.local_id;
    bool local_id_in_range = local_id < topology.local_face_count();
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
    if (first &&
        observation == 1) {
      local_id_in_range = false;
    }
#endif
    if (!local_id_in_range) {
      return false;
    }

    auto ownership = topology.face_ownership(local_id);
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
    if (first &&
        observation == 2) {
      ownership = mesh::EntityOwnership::ghost;
    }
#endif
    if (ownership != mesh::EntityOwnership::owned) {
      return false;
    }

    auto patch_id = topology.patch_id(local_id);
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
    if (first &&
        observation == 3) {
      patch_id = outlet_id == 0U ? 1U : 0U;
    }
#endif
    if (patch_id != std::optional<std::uint32_t>{outlet_id}) {
      return false;
    }

    bool patch_contains = outlet_patch.contains(local_id);
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
    if (first && observation == 4) {
      patch_contains = false;
    }
#endif
    if (!patch_contains) {
      return false;
    }

    auto global_id = topology.global_face_id(local_id);
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
    if (first &&
        observation == 5) {
      global_id = global_id == std::numeric_limits<mesh::GlobalFaceId>::max()
                      ? global_id - 1U
                      : global_id + 1U;
    }
#endif
    if (global_id != expected_face.global_id) {
      return false;
    }
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
    first = false;
#endif
  }
  return true;
}

std::array<double, 3> normalized(runtime::Real3 area) {
  require_finite(area, "boundary face area vector is not finite");
  const double scale =
      std::max({std::abs(area.x), std::abs(area.y), std::abs(area.z)});
  if (!(scale > 0.0)) {
    throw runtime::Error("boundary face area vector has zero magnitude");
  }
  const long double x = static_cast<long double>(area.x) / scale;
  const long double y = static_cast<long double>(area.y) / scale;
  const long double z = static_cast<long double>(area.z) / scale;
  const long double length = std::sqrt(x * x + y * y + z * z);
  if (!(length > 0.0L) || !std::isfinite(length)) {
    throw runtime::Error("boundary face normal cannot be normalized");
  }
  return {static_cast<double>(x / length), static_cast<double>(y / length),
          static_cast<double>(z / length)};
}

ScalarBoundaryValues copy_interior(double value) { return {value, value}; }

ScalarBoundaryValues prescribed(double interior, double boundary) {
  const double exterior =
      checked_long_double(2.0L * static_cast<long double>(boundary) - interior,
                          "boundary mirror value is not finite");
  return {boundary, exterior};
}

} // namespace

BoundaryDescriptor::BoundaryDescriptor(
    std::uint32_t stable_id, std::string name, BoundaryKind kind,
    VelocityRule velocity_rule, PressureRule pressure_rule,
    TransportRule density_rule, TransportRule enthalpy_rule,
    TransportRule scalar_rule, MassFluxRule mass_flux_rule,
    std::optional<std::uint32_t> paired_patch_id,
    std::optional<ResolvedInletState> inlet_state,
    std::optional<double> pressure_value_pa)
    : stable_id_(stable_id), name_(std::move(name)), kind_(kind),
      velocity_rule_(velocity_rule), pressure_rule_(pressure_rule),
      density_rule_(density_rule), enthalpy_rule_(enthalpy_rule),
      scalar_rule_(scalar_rule), mass_flux_rule_(mass_flux_rule),
      paired_patch_id_(paired_patch_id), inlet_state_(std::move(inlet_state)),
      pressure_value_pa_(pressure_value_pa) {}

BoundaryDescriptor::BoundaryDescriptor(BoundaryDescriptor &&) noexcept =
    default;
BoundaryDescriptor &
BoundaryDescriptor::operator=(BoundaryDescriptor &&) noexcept = default;

std::uint32_t BoundaryDescriptor::stable_id() const noexcept {
  return stable_id_;
}
std::string_view BoundaryDescriptor::name() const noexcept { return name_; }
BoundaryKind BoundaryDescriptor::kind() const noexcept { return kind_; }
VelocityRule BoundaryDescriptor::velocity_rule() const noexcept {
  return velocity_rule_;
}
PressureRule BoundaryDescriptor::pressure_rule() const noexcept {
  return pressure_rule_;
}
TransportRule BoundaryDescriptor::density_rule() const noexcept {
  return density_rule_;
}
TransportRule BoundaryDescriptor::enthalpy_rule() const noexcept {
  return enthalpy_rule_;
}
TransportRule BoundaryDescriptor::scalar_rule() const noexcept {
  return scalar_rule_;
}
MassFluxRule BoundaryDescriptor::mass_flux_rule() const noexcept {
  return mass_flux_rule_;
}
std::optional<std::uint32_t>
BoundaryDescriptor::paired_patch_id() const noexcept {
  return paired_patch_id_;
}
const std::optional<ResolvedInletState> &
BoundaryDescriptor::inlet_state() const noexcept {
  return inlet_state_;
}
std::optional<double> BoundaryDescriptor::pressure_value_pa() const noexcept {
  return pressure_value_pa_;
}

struct BoundaryRegistry::Impl final {
  std::vector<BoundaryDescriptor> patches;
  std::vector<std::string> scalar_names;
  bool open_domain{};
  std::optional<std::uint32_t> inlet_id;
  std::optional<std::uint32_t> outlet_id;
  TopologySignature topology;
  std::vector<OutletFace> owned_outlet_faces;
};

BoundaryRegistry::BoundaryRegistry(std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
BoundaryRegistry::~BoundaryRegistry() noexcept = default;
BoundaryRegistry::BoundaryRegistry(BoundaryRegistry &&) noexcept = default;

BoundaryRegistry
BoundaryRegistry::create(const config::FlowCaseConfig &resolved,
                         const mesh::MeshTopology &topology) try {
  if (resolved.schema_version != 2 ||
      resolved.simulation_type !=
          config::SimulationType::variable_density_flow) {
    throw runtime::Error(
        "boundary registry requires a resolved schema-v2 flow case");
  }
  if (!enum_in_range(resolved.density_model, config::DensityModel::ideal_gas)) {
    throw runtime::Error("boundary registry density model is invalid");
  }
  require_positive(resolved.physics.rho_ref_kg_per_m3,
                   "reference density must be positive and finite");
  require_positive(resolved.physics.inlet_consistency_rtol,
                   "inlet consistency tolerance must be positive and finite");
  if (resolved.density_model == config::DensityModel::ideal_gas) {
    if (!resolved.physics.cp_J_per_kg_K.has_value() ||
        !resolved.physics.gas_constant_J_per_kg_K.has_value() ||
        !resolved.physics.thermodynamic_pressure_pa.has_value()) {
      throw runtime::Error("ideal-gas boundary physics is incomplete");
    }
    require_positive(*resolved.physics.cp_J_per_kg_K,
                     "ideal-gas cp must be positive and finite");
    require_positive(*resolved.physics.gas_constant_J_per_kg_K,
                     "ideal-gas constant must be positive and finite");
    require_positive(*resolved.physics.thermodynamic_pressure_pa,
                     "thermodynamic pressure must be positive and finite");
  } else if (resolved.physics.cp_J_per_kg_K.has_value() ||
             resolved.physics.gas_constant_J_per_kg_K.has_value() ||
             resolved.physics.thermodynamic_pressure_pa.has_value()) {
    throw runtime::Error("non-ideal boundary case contains ideal-gas physics");
  }

  auto impl = std::make_unique<Impl>();
  impl->topology = topology_signature(topology);
  impl->scalar_names.reserve(resolved.scalars.size());
  for (const auto &scalar : resolved.scalars) {
    if (!valid_scalar_name(scalar.name)) {
      throw runtime::Error("boundary registry scalar name is invalid");
    }
    if (!impl->scalar_names.empty() &&
        impl->scalar_names.back() >= scalar.name) {
      throw runtime::Error(
          "boundary registry scalar names are not canonical and unique");
    }
    impl->scalar_names.push_back(scalar.name);
  }

  std::array<const config::FlowBoundaryConfig *, 6> ordered{};
  for (const auto &boundary : resolved.boundaries) {
    if (!enum_in_range(boundary.patch, config::PatchName::z_max)) {
      throw runtime::Error("boundary patch stable ID is invalid");
    }
    if (!enum_in_range(boundary.type, config::BoundaryType::pressure_outlet)) {
      throw runtime::Error("boundary kind is invalid");
    }
    const auto id = static_cast<std::size_t>(boundary.patch);
    if (ordered[id] != nullptr) {
      throw runtime::Error("boundary patch stable ID is duplicated");
    }
    ordered[id] = &boundary;
  }
  if (std::any_of(ordered.begin(), ordered.end(),
                  [](const auto *value) { return value == nullptr; })) {
    throw runtime::Error("boundary registry requires all six stable patches");
  }

  int inlet_count = 0;
  int outlet_count = 0;
  impl->patches.reserve(6U);
  for (std::size_t id = 0; id < ordered.size(); ++id) {
    const auto &boundary = *ordered[id];
    const auto &topology_patch = topology.patch(static_cast<std::uint32_t>(id));
    const bool config_periodic =
        boundary.type == config::BoundaryType::periodic;
    const bool topology_periodic =
        topology_patch.pairing_kind() == mesh::PatchPairingKind::periodic;
    if (config_periodic != topology_periodic) {
      throw runtime::Error(
          "resolved boundary periodicity does not match mesh topology");
    }
    const std::optional<std::uint32_t> expected_pair =
        config_periodic
            ? std::optional<std::uint32_t>{static_cast<std::uint32_t>(
                  id % 2U == 0U ? id + 1U : id - 1U)}
            : std::nullopt;
    if (topology_patch.paired_patch_id() != expected_pair) {
      throw runtime::Error(
          "mesh topology periodic patch pair is not reciprocal");
    }

    std::optional<ResolvedInletState> inlet_state;
    std::optional<double> pressure_value;
    if (boundary.type == config::BoundaryType::velocity_inlet) {
      ++inlet_count;
      impl->inlet_id = static_cast<std::uint32_t>(id);
      inlet_state = materialize_inlet(resolved, boundary, impl->scalar_names);
    } else if (boundary.type == config::BoundaryType::pressure_outlet) {
      ++outlet_count;
      impl->outlet_id = static_cast<std::uint32_t>(id);
      if (!boundary.pressure_perturbation_pa.has_value()) {
        throw runtime::Error("pressure outlet prescription is missing");
      }
      require_finite(*boundary.pressure_perturbation_pa,
                     "pressure outlet prescription is not finite");
      if (boundary.velocity_m_per_s.has_value() ||
          boundary.thermal_authority.has_value() ||
          boundary.temperature_K.has_value() ||
          boundary.enthalpy_J_per_kg.has_value() ||
          boundary.density_kg_per_m3.has_value() ||
          boundary.scalar_values.has_value()) {
        throw runtime::Error(
            "pressure outlet contains an incompatible prescription");
      }
      pressure_value = *boundary.pressure_perturbation_pa;
    } else if (has_any_prescription(boundary)) {
      throw runtime::Error(
          "closed or periodic patch contains an incompatible prescription");
    }

    const Rules rules = rules_for(boundary.type);
    impl->patches.push_back(BoundaryDescriptor{
        static_cast<std::uint32_t>(id), std::string(kPatchNames[id]),
        rules.kind, rules.velocity, rules.pressure, rules.density,
        rules.transport, rules.transport, rules.mass_flux, expected_pair,
        std::move(inlet_state), pressure_value});
  }

  if (!((inlet_count == 0 && outlet_count == 0) ||
        (inlet_count == 1 && outlet_count == 1))) {
    throw runtime::Error("boundary registry requires either a closed case or "
                         "one inlet/outlet pair");
  }
  impl->open_domain = inlet_count == 1;
  if (impl->open_domain) {
    const auto outlet_id = *impl->outlet_id;
    const auto &outlet_patch = topology.patch(outlet_id);
    impl->owned_outlet_faces.reserve(outlet_patch.local_faces().size());
    for (const auto local_face : outlet_patch.local_faces()) {
      if (topology.face_ownership(local_face) == mesh::EntityOwnership::owned) {
        impl->owned_outlet_faces.push_back(
            OutletFace{local_face, topology.global_face_id(local_face)});
      }
    }
  }
  return BoundaryRegistry(std::move(impl));
} catch (const runtime::Error &) {
  throw;
} catch (const std::bad_alloc &) {
  throw runtime::Error("boundary registry allocation failed");
} catch (const std::length_error &) {
  throw runtime::Error("boundary registry allocation size is unsupported");
}

const BoundaryDescriptor &
BoundaryRegistry::patch(std::uint32_t stable_id) const {
  if (stable_id >= impl_->patches.size()) {
    throw runtime::Error("boundary patch stable ID is outside [0, 6)");
  }
  return impl_->patches[stable_id];
}

std::size_t BoundaryRegistry::scalar_count() const noexcept {
  return impl_->scalar_names.size();
}

std::string_view BoundaryRegistry::scalar_name(std::size_t scalar) const {
  if (scalar >= impl_->scalar_names.size()) {
    throw runtime::Error("boundary scalar index is out of bounds");
  }
  return impl_->scalar_names[scalar];
}

bool BoundaryRegistry::open_domain() const noexcept {
  return impl_->open_domain;
}

std::optional<std::uint32_t>
BoundaryRegistry::velocity_inlet_patch_id() const noexcept {
  return impl_->inlet_id;
}

std::optional<std::uint32_t>
BoundaryRegistry::pressure_outlet_patch_id() const noexcept {
  return impl_->outlet_id;
}

VelocityBoundaryValues BoundaryRegistry::evaluate_velocity(
    std::uint32_t patch_id, runtime::Real3 interior_velocity,
    runtime::Real3 owner_outward_area_vector) const {
  const auto &descriptor = patch(patch_id);
  if (descriptor.velocity_rule() == VelocityRule::periodic_pair) {
    throw runtime::Error(
        "periodic velocity requires the topology-paired neighbour");
  }
  require_finite(interior_velocity, "boundary interior velocity is not finite");
  const auto normal = normalized(owner_outward_area_vector);
  switch (descriptor.velocity_rule()) {
  case VelocityRule::prescribed_zero:
    return {runtime::Real3{0.0, 0.0, 0.0},
            runtime::Real3{-interior_velocity.x, -interior_velocity.y,
                           -interior_velocity.z}};
  case VelocityRule::reflect_normal_copy_tangential: {
    const long double projection =
        static_cast<long double>(interior_velocity.x) * normal[0] +
        static_cast<long double>(interior_velocity.y) * normal[1] +
        static_cast<long double>(interior_velocity.z) * normal[2];
    const runtime::Real3 exterior{
        checked_long_double(static_cast<long double>(interior_velocity.x) -
                                2.0L * projection * normal[0],
                            "symmetry exterior velocity is not finite"),
        checked_long_double(static_cast<long double>(interior_velocity.y) -
                                2.0L * projection * normal[1],
                            "symmetry exterior velocity is not finite"),
        checked_long_double(static_cast<long double>(interior_velocity.z) -
                                2.0L * projection * normal[2],
                            "symmetry exterior velocity is not finite")};
    const runtime::Real3 face{
        checked_long_double(
            0.5L * (static_cast<long double>(interior_velocity.x) + exterior.x),
            "symmetry face velocity is not finite"),
        checked_long_double(
            0.5L * (static_cast<long double>(interior_velocity.y) + exterior.y),
            "symmetry face velocity is not finite"),
        checked_long_double(
            0.5L * (static_cast<long double>(interior_velocity.z) + exterior.z),
            "symmetry face velocity is not finite")};
    return {face, exterior};
  }
  case VelocityRule::prescribed_inlet: {
    const auto &inlet = descriptor.inlet_state();
    if (!inlet.has_value()) {
      throw runtime::Error("velocity inlet state is unavailable");
    }
    const auto prescribed_velocity = inlet->velocity_m_per_s;
    const runtime::Real3 exterior{
        checked_long_double(2.0L * prescribed_velocity.x - interior_velocity.x,
                            "inlet exterior velocity is not finite"),
        checked_long_double(2.0L * prescribed_velocity.y - interior_velocity.y,
                            "inlet exterior velocity is not finite"),
        checked_long_double(2.0L * prescribed_velocity.z - interior_velocity.z,
                            "inlet exterior velocity is not finite")};
    return {prescribed_velocity, exterior};
  }
  case VelocityRule::pure_outflow:
    return {interior_velocity, interior_velocity};
  default:
    throw runtime::Error("boundary velocity rule is invalid");
  }
}

ScalarBoundaryValues
BoundaryRegistry::evaluate_pressure(std::uint32_t patch_id,
                                    double interior_pi) const {
  const auto &descriptor = patch(patch_id);
  require_finite(interior_pi, "boundary interior pressure is not finite");
  switch (descriptor.pressure_rule()) {
  case PressureRule::periodic_pair:
    throw runtime::Error(
        "periodic pressure requires the topology-paired neighbour");
  case PressureRule::zero_normal_gradient:
    return copy_interior(interior_pi);
  case PressureRule::prescribed_value:
    if (!descriptor.pressure_value_pa().has_value()) {
      throw runtime::Error("pressure outlet value is unavailable");
    }
    return prescribed(interior_pi, *descriptor.pressure_value_pa());
  default:
    throw runtime::Error("boundary pressure rule is invalid");
  }
}

ScalarBoundaryValues
BoundaryRegistry::evaluate_density(std::uint32_t patch_id,
                                   double interior_density) const {
  const auto &descriptor = patch(patch_id);
  require_finite(interior_density, "boundary interior density is not finite");
  switch (descriptor.density_rule()) {
  case TransportRule::periodic_pair:
    throw runtime::Error(
        "periodic density requires the topology-paired neighbour");
  case TransportRule::copy_interior:
  case TransportRule::zero_normal_diffusive_flux:
  case TransportRule::pure_outflow:
    return copy_interior(interior_density);
  case TransportRule::prescribed_value:
    if (!descriptor.inlet_state().has_value()) {
      throw runtime::Error("velocity inlet density is unavailable");
    }
    return prescribed(interior_density,
                      descriptor.inlet_state()->density_kg_per_m3);
  default:
    throw runtime::Error("boundary density rule is invalid");
  }
}

ScalarBoundaryValues
BoundaryRegistry::evaluate_enthalpy(std::uint32_t patch_id,
                                    double interior_enthalpy) const {
  const auto &descriptor = patch(patch_id);
  require_finite(interior_enthalpy, "boundary interior enthalpy is not finite");
  switch (descriptor.enthalpy_rule()) {
  case TransportRule::periodic_pair:
    throw runtime::Error(
        "periodic enthalpy requires the topology-paired neighbour");
  case TransportRule::copy_interior:
  case TransportRule::zero_normal_diffusive_flux:
  case TransportRule::pure_outflow:
    return copy_interior(interior_enthalpy);
  case TransportRule::prescribed_value:
    if (!descriptor.inlet_state().has_value()) {
      throw runtime::Error("velocity inlet enthalpy is unavailable");
    }
    return prescribed(interior_enthalpy,
                      descriptor.inlet_state()->enthalpy_J_per_kg);
  default:
    throw runtime::Error("boundary enthalpy rule is invalid");
  }
}

ScalarBoundaryValues
BoundaryRegistry::evaluate_scalar(std::uint32_t patch_id, std::size_t scalar,
                                  double interior_value) const {
  if (scalar >= impl_->scalar_names.size()) {
    throw runtime::Error("boundary scalar index is out of bounds");
  }
  const auto &descriptor = patch(patch_id);
  require_finite(interior_value, "boundary interior scalar is not finite");
  switch (descriptor.scalar_rule()) {
  case TransportRule::periodic_pair:
    throw runtime::Error(
        "periodic scalar requires the topology-paired neighbour");
  case TransportRule::copy_interior:
  case TransportRule::zero_normal_diffusive_flux:
  case TransportRule::pure_outflow:
    return copy_interior(interior_value);
  case TransportRule::prescribed_value:
    if (!descriptor.inlet_state().has_value() ||
        descriptor.inlet_state()->scalar_values.size() !=
            impl_->scalar_names.size()) {
      throw runtime::Error("velocity inlet scalar state is unavailable");
    }
    return prescribed(interior_value,
                      descriptor.inlet_state()->scalar_values[scalar]);
  default:
    throw runtime::Error("boundary scalar rule is invalid");
  }
}

FinalFluxAdmissibility BoundaryRegistry::assess_final_pressure_outlet_flux(
    const mesh::MeshTopology &topology, const runtime::MpiContext &mpi,
    const runtime::FaceFieldView<const double> &final_face_mass_flux,
    std::uint64_t step, double time_s) const {
  if (!impl_->open_domain) {
    return {FinalFluxDecision::admissible, std::nullopt};
  }

  bool local_ok = true;
  std::string_view local_message;
  double local_minimum = 0.0;
  mesh::GlobalFaceId local_minimum_face =
      std::numeric_limits<mesh::GlobalFaceId>::max();
  bool local_negative = false;
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
  const auto topology_observation =
      detail::consume_next_outlet_topology_observation();
#endif
  try {
    if (!finite(time_s) || time_s < 0.0) {
      throw runtime::Error("final pressure-outlet assessment time is invalid");
    }
    if (!compatible(impl_->topology, topology, impl_->owned_outlet_faces,
                    *impl_->outlet_id
#ifdef HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS
                    , topology_observation
#endif
                    )) {
      throw runtime::Error(
          "final pressure-outlet assessment topology is incompatible");
    }
    if (final_face_mass_flux.face_count() != topology.local_face_count() ||
        final_face_mass_flux.components() != 1U) {
      throw runtime::Error(
          "final pressure-outlet face field layout is incompatible");
    }
    runtime::with_kernel_face_view(final_face_mass_flux, [&](auto flux) {
      for (const auto &outlet_face : impl_->owned_outlet_faces) {
        const double value = flux(outlet_face.local_id, 0);
        if (!finite(value)) {
          throw runtime::Error("final pressure-outlet face flux is not finite");
        }
        if (value < local_minimum ||
            (value == local_minimum && value < 0.0 &&
             outlet_face.global_id < local_minimum_face)) {
          local_minimum = value;
          local_minimum_face = outlet_face.global_id;
        }
        if (value < 0.0) {
          local_negative = true;
        }
      }
    });
  } catch (const runtime::Error &error) {
    local_ok = false;
    local_message = detail::fixed_preflight_message(error);
  } catch (const std::exception &error) {
    local_ok = false;
    local_message = detail::fixed_preflight_message(error);
  } catch (...) {
    local_ok = false;
    local_message = detail::fixed_unknown_preflight_message();
  }
  const auto status = runtime::collective_status(mpi, local_ok, local_message);
  if (!status.ok) {
    throw runtime::Error(status.message);
  }

  double severity = local_negative ? -local_minimum : 0.0;
  mpi.allreduce_fp64_in_place(&severity, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  if (severity == 0.0) {
    return {FinalFluxDecision::admissible, std::nullopt};
  }

  const bool tied = local_negative && -local_minimum == severity;
  const std::uint32_t local_high =
      static_cast<std::uint32_t>(local_minimum_face >> 32U);
  double high_score = tied ? kLimbBase - local_high : 0.0;
  mpi.allreduce_fp64_in_place(&high_score, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  const auto selected_high = static_cast<std::uint32_t>(kLimbBase - high_score);

  const std::uint32_t local_low =
      static_cast<std::uint32_t>(local_minimum_face & 0xffffffffULL);
  double low_score =
      tied && local_high == selected_high ? kLimbBase - local_low : 0.0;
  mpi.allreduce_fp64_in_place(&low_score, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  const auto selected_low = static_cast<std::uint32_t>(kLimbBase - low_score);

  double rank_score =
      local_negative ? static_cast<double>(mpi.size() - mpi.rank()) : 0.0;
  mpi.allreduce_fp64_in_place(&rank_score, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  const int lowest_failing_rank = mpi.size() - static_cast<int>(rank_score);
  const mesh::GlobalFaceId selected_face =
      (static_cast<mesh::GlobalFaceId>(selected_high) << 32U) | selected_low;
  return {FinalFluxDecision::outlet_backflow,
          OutletBackflowEvidence{*impl_->outlet_id, step, time_s, -severity,
                                 selected_face, lowest_failing_rank}};
}

} // namespace hundun::boundary
