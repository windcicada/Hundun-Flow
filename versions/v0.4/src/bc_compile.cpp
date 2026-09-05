// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "bc_identity_detail.hpp"
#include "hundun/v04_boundary.hpp"
#include "mesh_focus_detail.hpp"

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kBoundaryInput = 401U;
constexpr std::uint32_t kBoundaryCollective = 402U;
constexpr std::uint32_t kBoundaryPeriodic = 403U;
constexpr std::uint32_t kBoundaryClosure = 404U;
constexpr std::uint32_t kBoundaryScalar = 405U;
constexpr std::uint32_t kBoundaryScheme = 406U;
constexpr std::uint32_t kBoundaryRegistry = 407U;
constexpr std::uint32_t kBoundaryOverflow = 408U;
constexpr std::size_t kMaxTransportedScalars = 64U;
constexpr std::size_t kMaxStableNameBytes = 255U;

constexpr std::array<CartesianFace, 6U> kFaces{
    CartesianFace::x_min, CartesianFace::x_max, CartesianFace::y_min,
    CartesianFace::y_max, CartesianFace::z_min, CartesianFace::z_max};

struct BoundaryBuildData {
  std::array<BoundaryFacePlan, 6U> faces{};
  std::vector<BoundaryIndexSpan> spans;
  std::vector<BoundaryKernelBatch> batches;
  std::vector<BoundaryTransportedField> transported_fields;
  std::vector<double> velocity_x, velocity_y, velocity_z;
  std::vector<double> direction_x, direction_y, direction_z;
  std::vector<double> backflow_velocity_x, backflow_velocity_y,
      backflow_velocity_z;
  std::vector<double> scalar_targets, scalar_backflow_targets,
      pressure_targets, temperature_targets;
  std::vector<double> total_pressure_targets, total_temperature_targets;
  std::vector<double> backflow_temperature_targets;
  std::vector<double> heat_flux_targets, mass_flow_targets;
  std::vector<double> relaxation_rates, mach_limits;
  std::vector<std::uint8_t> allow_backflow;
  std::vector<BoundaryParameterRole> parameter_roles;
  std::vector<ScalarBoundaryKind> scalar_kinds;
  std::vector<BoundaryScalarRole> scalar_roles;
  std::vector<double> normal_distance_1, normal_distance_2;
  HaloTopology halo_topology{};
  FieldId velocity_field{}, pressure_field{}, enthalpy_field{};
  PressureReferenceKind pressure_reference{
      PressureReferenceKind::boundary_absolute};
  Int3 local_cells{};
  std::size_t resolved_scalar_count{}, resolved_vector_count{};
  std::size_t resolved_normal_gradient_count{};
  std::uint8_t required_ghost_width{};
  RevisionToken revision{};
  PlanFingerprint semantic_fingerprint{}, local_layout_fingerprint{};
};

std::size_t face_index(CartesianFace face) noexcept {
  return static_cast<std::size_t>(face);
}

int face_axis(CartesianFace face) noexcept {
  const std::size_t index = face_index(face);
  return index < 2U ? 0 : (index < 4U ? 1 : 2);
}

CartesianFace opposite(CartesianFace face) noexcept {
  const std::size_t index = face_index(face);
  return static_cast<CartesianFace>(index ^ 1U);
}

std::int32_t component(Int3 value, int axis) noexcept {
  return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

bool finite(double value) noexcept { return std::isfinite(value); }

bool finite(Real3 value) noexcept {
  return finite(value.x) && finite(value.y) && finite(value.z);
}

bool valid_flow_kind(BoundaryKind kind) noexcept {
  return kind >= BoundaryKind::velocity_inlet &&
         kind <= BoundaryKind::periodic;
}

bool valid_thermal_kind(BoundaryKind kind) noexcept {
  return kind == BoundaryKind::none ||
         kind == BoundaryKind::adiabatic_wall ||
         kind == BoundaryKind::isothermal_wall ||
         kind == BoundaryKind::heat_flux_wall;
}

bool valid_scalar_kind(ScalarBoundaryKind kind) noexcept {
  return static_cast<std::uint8_t>(kind) <=
         static_cast<std::uint8_t>(ScalarBoundaryKind::convective);
}

bool valid_stable_name(std::string_view value) noexcept {
  if (value.empty() || value.size() > kMaxStableNameBytes) {
    return false;
  }
  for (const unsigned char character : value) {
    if (std::isalnum(character) == 0 && character != '_' &&
        character != '-') {
      return false;
    }
  }
  return true;
}

bool valid_scalar_catalog(const ValidatedModel& model) noexcept {
  if (model.transported_scalars.size() > kMaxTransportedScalars) {
    return false;
  }
  for (std::size_t index = 0U; index < model.transported_scalars.size();
       ++index) {
    const TransportedScalarSpec& scalar = model.transported_scalars[index];
    if (!valid_stable_name(scalar.stable_name) || scalar.stable_name == "U" ||
        scalar.stable_name == "pi" || scalar.stable_name == "h" ||
        static_cast<std::uint8_t>(scalar.role) >
            static_cast<std::uint8_t>(TransportedScalarRole::passive_scalar)) {
      return false;
    }
    for (std::size_t prior = 0U; prior < index; ++prior) {
      if (model.transported_scalars[prior].stable_name == scalar.stable_name) {
        return false;
      }
    }
  }
  for (const BoundaryFaceSpec& face : model.boundaries) {
    if (face.scalars.size() > kMaxTransportedScalars) {
      return false;
    }
    if (face.flow_kind == BoundaryKind::periodic) {
      if (!face.scalars.empty()) {
        return false;
      }
      continue;
    }
    if (face.scalars.size() != model.transported_scalars.size()) {
      return false;
    }
    for (const TransportedScalarSpec& declared : model.transported_scalars) {
      std::size_t matches = 0U;
      for (const ScalarBoundarySpec& closure : face.scalars) {
        matches += closure.stable_name == declared.stable_name ? 1U : 0U;
      }
      if (matches != 1U) {
        return false;
      }
    }
  }
  return true;
}

double outward_component(CartesianFace face, Real3 value) noexcept {
  switch (face) {
    case CartesianFace::x_min:
      return -value.x;
    case CartesianFace::x_max:
      return value.x;
    case CartesianFace::y_min:
      return -value.y;
    case CartesianFace::y_max:
      return value.y;
    case CartesianFace::z_min:
      return -value.z;
    case CartesianFace::z_max:
      return value.z;
  }
  return 0.0;
}

bool impermeable_boundary(BoundaryKind kind) noexcept {
  return kind == BoundaryKind::no_slip_wall ||
         kind == BoundaryKind::moving_wall || kind == BoundaryKind::slip ||
         kind == BoundaryKind::symmetry || kind == BoundaryKind::periodic;
}

bool valid_convection(ConvectionScheme scheme) noexcept {
  return static_cast<std::uint8_t>(scheme) <=
         static_cast<std::uint8_t>(ConvectionScheme::tvd2);
}

bool flow_wall(BoundaryKind kind) noexcept {
  return kind == BoundaryKind::no_slip_wall ||
         kind == BoundaryKind::moving_wall || kind == BoundaryKind::slip ||
         kind == BoundaryKind::symmetry;
}

bool pressure_authority(BoundaryKind kind) noexcept {
  return kind == BoundaryKind::pressure_outlet ||
         kind == BoundaryKind::nscbc_outlet ||
         kind == BoundaryKind::static_state_inlet ||
         kind == BoundaryKind::total_state_inlet;
}

bool inlet_composition_authority(BoundaryKind kind) noexcept {
  return kind == BoundaryKind::velocity_inlet ||
         kind == BoundaryKind::mass_flow_inlet ||
         kind == BoundaryKind::static_state_inlet ||
         kind == BoundaryKind::total_state_inlet ||
         kind == BoundaryKind::nscbc_inlet;
}

const TransportedScalarSpec* declared_scalar(
    const ValidatedModel& model, std::string_view stable_name) noexcept {
  const auto found = std::find_if(
      model.transported_scalars.begin(), model.transported_scalars.end(),
      [stable_name](const TransportedScalarSpec& candidate) {
        return candidate.stable_name == stable_name;
      });
  return found == model.transported_scalars.end() ? nullptr : &*found;
}

bool valid_species_composition(const ValidatedModel& model,
                               const BoundaryFaceSpec& face) noexcept {
  const bool complete_inlet = inlet_composition_authority(face.flow_kind);
  double dirichlet_sum = 0.0;
  for (const ScalarBoundarySpec& closure : face.scalars) {
    const TransportedScalarSpec* declared =
        declared_scalar(model, closure.stable_name);
    if (declared == nullptr ||
        declared->role != TransportedScalarRole::species) {
      continue;
    }
    if (complete_inlet && closure.kind != ScalarBoundaryKind::dirichlet) {
      return false;
    }
    if (closure.kind == ScalarBoundaryKind::dirichlet) {
      if (closure.value < 0.0 || closure.value > 1.0) {
        return false;
      }
      dirichlet_sum += closure.value;
    }
  }
  return dirichlet_sum <= 1.0;
}

bool resolved_pressure_ghost(BoundaryKind kind) noexcept {
  return kind == BoundaryKind::pressure_outlet ||
         kind == BoundaryKind::static_state_inlet ||
         kind == BoundaryKind::total_state_inlet;
}

BoundaryFlowKernel flow_kernel(BoundaryKind kind) noexcept {
  switch (kind) {
    case BoundaryKind::mass_flow_inlet:
      return BoundaryFlowKernel::mass_flow_constraint;
    case BoundaryKind::static_state_inlet:
      return BoundaryFlowKernel::static_state;
    case BoundaryKind::total_state_inlet:
      return BoundaryFlowKernel::total_state;
    case BoundaryKind::pressure_outlet:
      return BoundaryFlowKernel::pressure_backflow;
    case BoundaryKind::nscbc_inlet:
    case BoundaryKind::nscbc_outlet:
      return BoundaryFlowKernel::characteristic;
    default:
      return BoundaryFlowKernel::ordinary;
  }
}

bool local_owner(const MeshPatch& patch, CartesianFace face) noexcept {
  const int axis = face_axis(face);
  const bool minimum = (face_index(face) & 1U) == 0U;
  const std::int32_t coordinate = component(patch.process_coord, axis);
  const std::int32_t partitions = component(patch.process_grid, axis);
  return minimum ? coordinate == 0 : coordinate + 1 == partitions;
}

std::uint8_t scheme_reach(ConvectionScheme scheme) noexcept {
  return scheme == ConvectionScheme::central2 ? 1U : 2U;
}

std::uint8_t required_ghost_width(const SchemeSpec& schemes) noexcept {
  return std::max(
      {scheme_reach(schemes.momentum), scheme_reach(schemes.enthalpy),
       scheme_reach(schemes.species),
       scheme_reach(schemes.passive_scalar), std::uint8_t{1U}});
}

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= 1099511628211ULL;
  return hash;
}

std::uint64_t double_bits(double value) noexcept {
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

Status collective_status(MPI_Comm communicator, int rank, int size,
                         Status local, int& lowest) noexcept {
  const std::uint64_t candidate =
      local.code == StatusCode::ok
          ? std::numeric_limits<std::uint64_t>::max()
          : static_cast<std::uint64_t>(rank);
  std::uint64_t selected = std::numeric_limits<std::uint64_t>::max();
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kBoundaryCollective};
  }
  if (selected == std::numeric_limits<std::uint64_t>::max()) {
    lowest = -1;
    return {};
  }
  if (selected >= static_cast<std::uint64_t>(size)) {
    lowest = -1;
    return {StatusCode::mpi_failure, kBoundaryCollective};
  }
  const int failing_rank = static_cast<int>(selected);
  std::array<std::uint64_t, 2U> wire{};
  if (rank == failing_rank) {
    wire[0] = static_cast<std::uint64_t>(local.code);
    wire[1] = local.detail;
  }
  if (MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT64_T,
                failing_rank, communicator) != MPI_SUCCESS ||
      wire[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      wire[1] > std::numeric_limits<std::uint32_t>::max()) {
    lowest = -1;
    return {StatusCode::mpi_failure, kBoundaryCollective};
  }
  lowest = failing_rank;
  return {static_cast<StatusCode>(wire[0]),
          static_cast<std::uint32_t>(wire[1])};
}

Status validate_geometry_patch(const CartesianGeometryPlan& geometry,
                               const MeshPatch& patch, int rank, int size,
                               std::uint8_t ghost_width) noexcept {
  const Int3 global = geometry.global_cells();
  MeshPatch authoritative;
  const Status decomposition =
      detail::make_mesh_patch(rank, size, global, authoritative);
  if (geometry.fingerprint() == 0U || geometry.topology_revision() == 0U ||
      !decomposition || global.x <= 0 || global.y <= 0 || global.z <= 0 ||
      patch.begin.x != authoritative.begin.x ||
      patch.begin.y != authoritative.begin.y ||
      patch.begin.z != authoritative.begin.z ||
      patch.cells.x != authoritative.cells.x ||
      patch.cells.y != authoritative.cells.y ||
      patch.cells.z != authoritative.cells.z ||
      patch.process_grid.x != authoritative.process_grid.x ||
      patch.process_grid.y != authoritative.process_grid.y ||
      patch.process_grid.z != authoritative.process_grid.z ||
      patch.process_coord.x != authoritative.process_coord.x ||
      patch.process_coord.y != authoritative.process_coord.y ||
      patch.process_coord.z != authoritative.process_coord.z ||
      patch.cells.y <= 0 || patch.cells.z <= 0 || patch.process_grid.x <= 0 ||
      patch.process_grid.y <= 0 || patch.process_grid.z <= 0 ||
      patch.process_coord.x < 0 || patch.process_coord.y < 0 ||
      patch.process_coord.z < 0 ||
      patch.process_coord.x >= patch.process_grid.x ||
      patch.process_coord.y >= patch.process_grid.y ||
      patch.process_coord.z >= patch.process_grid.z || ghost_width == 0U ||
      patch.cells.x < ghost_width || patch.cells.y < ghost_width ||
      patch.cells.z < ghost_width) {
    return {StatusCode::invalid_plan, kBoundaryInput};
  }
  const std::int64_t patch_end_x =
      static_cast<std::int64_t>(patch.begin.x) + patch.cells.x;
  const std::int64_t patch_end_y =
      static_cast<std::int64_t>(patch.begin.y) + patch.cells.y;
  const std::int64_t patch_end_z =
      static_cast<std::int64_t>(patch.begin.z) + patch.cells.z;
  const std::int64_t process_count =
      static_cast<std::int64_t>(patch.process_grid.x) * patch.process_grid.y *
      patch.process_grid.z;
  const std::int64_t expected_rank =
      patch.process_coord.x +
      static_cast<std::int64_t>(patch.process_grid.x) *
          (patch.process_coord.y +
           static_cast<std::int64_t>(patch.process_grid.y) *
               patch.process_coord.z);
  if (patch_end_x > global.x || patch_end_y > global.y ||
      patch_end_z > global.z || process_count != size ||
      expected_rank != rank) {
    return {StatusCode::invalid_plan, kBoundaryInput};
  }
  return {};
}

Status validate_model(const ValidatedModel& model) noexcept {
  if (model.fingerprint == 0U ||
      !valid_convection(model.schemes.momentum) ||
      !valid_convection(model.schemes.enthalpy) ||
      !valid_convection(model.schemes.species) ||
      !valid_convection(model.schemes.passive_scalar) ||
      model.schemes.diffusion != DiffusionScheme::central2 ||
      !finite(model.schemes.limiter) || model.schemes.limiter < 0.0 ||
      model.schemes.limiter > 1.0 ||
      !valid_scalar_catalog(model) ||
      static_cast<std::uint8_t>(model.pressure_reference) >
          static_cast<std::uint8_t>(PressureReferenceKind::closed_mass)) {
    return {StatusCode::invalid_plan, kBoundaryScheme};
  }
  bool has_pressure_authority = false;
  for (std::size_t index = 0U; index < model.boundaries.size(); ++index) {
    const BoundaryFaceSpec& spec = model.boundaries[index];
    if (!valid_flow_kind(spec.flow_kind) ||
        !valid_thermal_kind(spec.thermal_kind) ||
        !finite(spec.velocity) || !finite(spec.direction) ||
        !finite(spec.backflow_velocity) || !finite(spec.mass_flow_rate) ||
        !finite(spec.pressure) || !finite(spec.temperature) ||
        !finite(spec.total_pressure) || !finite(spec.total_temperature) ||
        !finite(spec.backflow_temperature) || !finite(spec.heat_flux) ||
        !finite(spec.relaxation) || !finite(spec.mach_limit) ||
        spec.relaxation < 0.0 || spec.mach_limit <= 0.0 ||
        spec.mach_limit >= 1.0) {
      return {StatusCode::invalid_plan, kBoundaryInput};
    }
    const CartesianFace face = kFaces[index];
    if (spec.flow_kind == BoundaryKind::periodic) {
      const BoundaryFaceSpec& paired =
          model.boundaries[face_index(opposite(face))];
      if (paired.flow_kind != BoundaryKind::periodic ||
          spec.thermal_kind != BoundaryKind::none || !spec.scalars.empty()) {
        return {StatusCode::invalid_plan, kBoundaryPeriodic};
      }
    }
    if (spec.thermal_kind != BoundaryKind::none &&
        !flow_wall(spec.flow_kind)) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if (spec.thermal_kind == BoundaryKind::isothermal_wall &&
        spec.temperature <= 0.0) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if (spec.flow_kind == BoundaryKind::velocity_inlet &&
        (spec.temperature <= 0.0 ||
         !(outward_component(face, spec.velocity) < 0.0))) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if ((spec.flow_kind == BoundaryKind::pressure_outlet ||
         spec.flow_kind == BoundaryKind::nscbc_outlet) &&
        spec.pressure <= 0.0) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if (spec.flow_kind == BoundaryKind::mass_flow_inlet &&
        (spec.mass_flow_rate <= 0.0 || spec.temperature <= 0.0)) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if ((spec.flow_kind == BoundaryKind::static_state_inlet ||
         spec.flow_kind == BoundaryKind::nscbc_inlet) &&
        (spec.pressure <= 0.0 || spec.temperature <= 0.0)) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if (spec.flow_kind == BoundaryKind::total_state_inlet &&
        (spec.total_pressure <= 0.0 || spec.total_temperature <= 0.0)) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if (spec.flow_kind == BoundaryKind::no_slip_wall &&
        (spec.velocity.x != 0.0 || spec.velocity.y != 0.0 ||
         spec.velocity.z != 0.0)) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if ((spec.flow_kind == BoundaryKind::pressure_outlet ||
         spec.flow_kind == BoundaryKind::nscbc_outlet) &&
        spec.allow_backflow && spec.backflow_temperature <= 0.0) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if ((spec.flow_kind == BoundaryKind::mass_flow_inlet ||
         spec.flow_kind == BoundaryKind::static_state_inlet ||
         spec.flow_kind == BoundaryKind::total_state_inlet) &&
        !(outward_component(face, spec.direction) < 0.0)) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if (spec.flow_kind == BoundaryKind::nscbc_inlet &&
        !(outward_component(face, spec.velocity) < 0.0)) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    if ((spec.flow_kind == BoundaryKind::pressure_outlet ||
         spec.flow_kind == BoundaryKind::nscbc_outlet) &&
        spec.allow_backflow &&
        !(outward_component(face, spec.backflow_velocity) < 0.0)) {
      return {StatusCode::invalid_plan, kBoundaryClosure};
    }
    for (std::size_t scalar = 0U; scalar < spec.scalars.size(); ++scalar) {
      const ScalarBoundarySpec& value = spec.scalars[scalar];
      if (value.stable_name.empty() || !valid_scalar_kind(value.kind) ||
          !valid_scalar_kind(value.backflow_kind) || !finite(value.value) ||
          !finite(value.backflow_value)) {
        return {StatusCode::invalid_plan, kBoundaryScalar};
      }
      for (std::size_t prior = 0U; prior < scalar; ++prior) {
        if (spec.scalars[prior].stable_name == value.stable_name) {
          return {StatusCode::invalid_plan, kBoundaryScalar};
        }
      }
    }
    if (!valid_species_composition(model, spec)) {
      return {StatusCode::invalid_plan, kBoundaryScalar};
    }
    if ((spec.flow_kind == BoundaryKind::pressure_outlet ||
         spec.flow_kind == BoundaryKind::nscbc_outlet) &&
        spec.allow_backflow) {
      double independent_species_sum = 0.0;
      for (const ScalarBoundarySpec& value : spec.scalars) {
        const TransportedScalarSpec* declared =
            declared_scalar(model, value.stable_name);
        if (declared == nullptr ||
            value.backflow_kind != ScalarBoundaryKind::dirichlet) {
          return {StatusCode::invalid_plan, kBoundaryScalar};
        }
        if (declared->role == TransportedScalarRole::species) {
          if (value.backflow_value < 0.0 || value.backflow_value > 1.0) {
            return {StatusCode::invalid_plan, kBoundaryScalar};
          }
          independent_species_sum += value.backflow_value;
        }
      }
      if (independent_species_sum > 1.0) {
        return {StatusCode::invalid_plan, kBoundaryScalar};
      }
    }
    has_pressure_authority |= pressure_authority(spec.flow_kind);
  }
  if ((model.pressure_reference ==
           PressureReferenceKind::boundary_absolute &&
       !has_pressure_authority) ||
      (model.pressure_reference == PressureReferenceKind::closed_mass &&
       has_pressure_authority)) {
    return {StatusCode::invalid_plan, kBoundaryClosure};
  }
  if (model.pressure_reference == PressureReferenceKind::closed_mass) {
    for (std::size_t index = 0U; index < model.boundaries.size(); ++index) {
      const BoundaryFaceSpec& face = model.boundaries[index];
      if (!impermeable_boundary(face.flow_kind) ||
          (face.flow_kind == BoundaryKind::moving_wall &&
           outward_component(kFaces[index], face.velocity) != 0.0)) {
        return {StatusCode::invalid_plan, kBoundaryClosure};
      }
    }
  }
  return {};
}

BoundaryRelation flow_relation(BoundaryKind kind) noexcept {
  if (kind == BoundaryKind::velocity_inlet ||
      kind == BoundaryKind::moving_wall ||
      kind == BoundaryKind::no_slip_wall) {
    return BoundaryRelation::dirichlet;
  }
  if (kind == BoundaryKind::slip || kind == BoundaryKind::symmetry) {
    return BoundaryRelation::reflect_normal;
  }
  return BoundaryRelation::zero_gradient;
}

BoundaryRelation scalar_relation(ScalarBoundaryKind kind) noexcept {
  switch (kind) {
    case ScalarBoundaryKind::dirichlet:
      return BoundaryRelation::dirichlet;
    case ScalarBoundaryKind::normal_flux:
      return BoundaryRelation::normal_gradient;
    case ScalarBoundaryKind::zero_gradient:
      return BoundaryRelation::zero_gradient;
    case ScalarBoundaryKind::convective:
      return BoundaryRelation::convective;
  }
  return BoundaryRelation::zero_gradient;
}

Status tangent_counts(const MeshPatch& patch, CartesianFace face,
                      std::uint32_t& inner,
                      std::uint32_t& outer) noexcept {
  const auto checked = [](std::int32_t value,
                          std::uint32_t& out) noexcept -> bool {
    if (value <= 0) {
      return false;
    }
    out = static_cast<std::uint32_t>(value);
    return true;
  };
  if (face_axis(face) == 0) {
    return checked(patch.cells.y, inner) && checked(patch.cells.z, outer)
               ? Status{}
               : Status{StatusCode::invalid_plan, kBoundaryOverflow};
  }
  if (face_axis(face) == 1) {
    return checked(patch.cells.x, inner) && checked(patch.cells.z, outer)
               ? Status{}
               : Status{StatusCode::invalid_plan, kBoundaryOverflow};
  }
  return checked(patch.cells.x, inner) && checked(patch.cells.y, outer)
             ? Status{}
             : Status{StatusCode::invalid_plan, kBoundaryOverflow};
}

Status add_span(BoundaryBuildData& plan, BoundaryStage stage,
                BoundaryRelation relation, CartesianFace face, FieldId field,
                std::uint8_t component_begin, std::uint8_t component_count,
                std::uint8_t ghosts, std::uint32_t parameter,
                BoundaryValueSource value_source, const MeshPatch& patch) {
  std::uint32_t inner = 0U;
  std::uint32_t outer = 0U;
  const Status counts = tangent_counts(patch, face, inner, outer);
  if (!counts) {
    return counts;
  }
  const std::uint64_t face_elements =
      static_cast<std::uint64_t>(inner) * static_cast<std::uint64_t>(outer);
  if (face_elements == 0U ||
      face_elements > std::numeric_limits<std::uint32_t>::max()) {
    return {StatusCode::invalid_plan, kBoundaryOverflow};
  }
  std::uint32_t resolved_begin = kInvalidBoundaryParameter;
  std::uint32_t resolved_stride = 0U;
  std::size_t* resolved_count = nullptr;
  switch (value_source) {
    case BoundaryValueSource::resolved_scalar:
      resolved_count = &plan.resolved_scalar_count;
      break;
    case BoundaryValueSource::resolved_vector:
      resolved_count = &plan.resolved_vector_count;
      break;
    case BoundaryValueSource::resolved_normal_gradient:
      resolved_count = &plan.resolved_normal_gradient_count;
      break;
    case BoundaryValueSource::none:
    case BoundaryValueSource::compiled_scalar:
    case BoundaryValueSource::compiled_vector:
      break;
  }
  if (resolved_count != nullptr) {
    const std::size_t stride = static_cast<std::size_t>(face_elements);
    if (*resolved_count > std::numeric_limits<std::uint32_t>::max() ||
        stride > std::numeric_limits<std::uint32_t>::max() -
                     *resolved_count) {
      return {StatusCode::invalid_plan, kBoundaryOverflow};
    }
    resolved_begin = static_cast<std::uint32_t>(*resolved_count);
    resolved_stride = static_cast<std::uint32_t>(stride);
    *resolved_count += stride;
  }
  plan.spans.push_back(BoundaryIndexSpan{
      stage, relation, face, field, component_begin, component_count, ghosts,
      inner, outer, parameter, resolved_begin, value_source,
      resolved_stride});
  return {};
}

Status append_parameter(BoundaryBuildData& plan,
                        const BoundaryFaceSpec& spec,
                        BoundaryParameterRole role,
                        ScalarBoundaryKind scalar_kind,
                        BoundaryScalarRole scalar_role,
                        double scalar_target,
                        double distance_1, double distance_2,
                        std::uint32_t& out) {
  if (plan.scalar_targets.size() >=
      static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    return {StatusCode::invalid_plan, kBoundaryOverflow};
  }
  out = static_cast<std::uint32_t>(plan.scalar_targets.size());
  plan.velocity_x.push_back(spec.velocity.x);
  plan.velocity_y.push_back(spec.velocity.y);
  plan.velocity_z.push_back(spec.velocity.z);
  plan.direction_x.push_back(spec.direction.x);
  plan.direction_y.push_back(spec.direction.y);
  plan.direction_z.push_back(spec.direction.z);
  plan.backflow_velocity_x.push_back(spec.backflow_velocity.x);
  plan.backflow_velocity_y.push_back(spec.backflow_velocity.y);
  plan.backflow_velocity_z.push_back(spec.backflow_velocity.z);
  plan.scalar_targets.push_back(scalar_target);
  plan.scalar_backflow_targets.push_back(scalar_target);
  plan.pressure_targets.push_back(spec.pressure);
  plan.temperature_targets.push_back(spec.temperature);
  plan.total_pressure_targets.push_back(spec.total_pressure);
  plan.total_temperature_targets.push_back(spec.total_temperature);
  plan.backflow_temperature_targets.push_back(spec.backflow_temperature);
  plan.heat_flux_targets.push_back(spec.heat_flux);
  plan.mass_flow_targets.push_back(spec.mass_flow_rate);
  plan.relaxation_rates.push_back(spec.relaxation);
  plan.mach_limits.push_back(spec.mach_limit);
  plan.allow_backflow.push_back(spec.allow_backflow ? 1U : 0U);
  plan.parameter_roles.push_back(role);
  plan.scalar_kinds.push_back(scalar_kind);
  plan.scalar_roles.push_back(scalar_role);
  plan.normal_distance_1.push_back(distance_1);
  plan.normal_distance_2.push_back(distance_2);
  return {};
}

Status normal_distances(const CartesianGeometryPlan& geometry,
                        CartesianFace face, double& first,
                        double& second) noexcept {
  const CartesianAxis axis = face_axis(face) == 0
                                 ? CartesianAxis::x
                                 : (face_axis(face) == 1 ? CartesianAxis::y
                                                         : CartesianAxis::z);
  const AxisMetrics& metrics = geometry.axis(axis);
  const Span<const double> faces = metrics.faces();
  const Span<const double> centres = metrics.centres();
  if (faces.size < 2U || centres.size + 1U != faces.size) {
    return {StatusCode::invalid_plan, kBoundaryInput};
  }
  const bool minimum = (face_index(face) & 1U) == 0U;
  if (minimum) {
    first = 2.0 * (centres.data[0] - faces.data[0]);
    second = centres.size > 1U
                 ? 2.0 * (centres.data[1] - faces.data[0])
                 : 3.0 * first;
  } else {
    const std::size_t last = centres.size - 1U;
    first = 2.0 * (faces.data[faces.size - 1U] - centres.data[last]);
    second = centres.size > 1U
                 ? 2.0 * (faces.data[faces.size - 1U] -
                          centres.data[last - 1U])
                 : 3.0 * first;
  }
  return std::isfinite(first) && std::isfinite(second) && first > 0.0 &&
                 second > first
             ? Status{}
             : Status{StatusCode::invalid_plan, kBoundaryInput};
}

Status compile_local_plan(const ValidatedModel& model,
                          const CartesianGeometryPlan& geometry,
                          const MeshPatch& patch,
                          FieldId velocity, FieldId pressure, FieldId enthalpy,
                          const std::vector<std::pair<std::string, FieldId>>&
                              scalar_fields,
                          BoundaryBuildData& plan) {
  plan.velocity_field = velocity;
  plan.pressure_field = pressure;
  plan.enthalpy_field = enthalpy;
  plan.pressure_reference = model.pressure_reference;
  plan.local_cells = patch.cells;
  plan.required_ghost_width = required_ghost_width(model.schemes);
  plan.transported_fields.reserve(model.transported_scalars.size());
  for (const TransportedScalarSpec& declared : model.transported_scalars) {
    const auto found = std::lower_bound(
        scalar_fields.begin(), scalar_fields.end(), declared.stable_name,
        [](const auto& entry, const std::string& name) {
          return entry.first < name;
        });
    if (found == scalar_fields.end() ||
        found->first != declared.stable_name) {
      return {StatusCode::invalid_plan, kBoundaryScalar};
    }
    plan.transported_fields.push_back(BoundaryTransportedField{
        found->second,
        declared.role == TransportedScalarRole::species
            ? BoundaryScalarRole::species
            : BoundaryScalarRole::passive_scalar});
  }
  plan.halo_topology = {
      model.boundaries[face_index(CartesianFace::x_min)].flow_kind ==
          BoundaryKind::periodic,
      model.boundaries[face_index(CartesianFace::y_min)].flow_kind ==
          BoundaryKind::periodic,
      model.boundaries[face_index(CartesianFace::z_min)].flow_kind ==
          BoundaryKind::periodic};
  for (const CartesianFace face : kFaces) {
    const std::size_t index = face_index(face);
    const BoundaryFaceSpec& spec = model.boundaries[index];
    BoundaryFacePlan& face_plan = plan.faces[index];
    face_plan.flow_kind = spec.flow_kind;
    face_plan.thermal_kind = spec.thermal_kind;
    face_plan.flow_kernel = flow_kernel(spec.flow_kind);
    face_plan.periodic = spec.flow_kind == BoundaryKind::periodic;
    face_plan.local_owner = local_owner(patch, face);
    if (face_plan.periodic || !face_plan.local_owner) {
      continue;
    }
    double distance_1 = 0.0;
    double distance_2 = 0.0;
    Status status =
        normal_distances(geometry, face, distance_1, distance_2);
    if (!status) {
      return status;
    }
    std::uint32_t flow_parameter = 0U;
    status = append_parameter(plan, spec, BoundaryParameterRole::flow,
                              ScalarBoundaryKind::zero_gradient,
                              BoundaryScalarRole::none,
                              spec.pressure, distance_1, distance_2,
                              flow_parameter);
    if (!status) {
      return status;
    }
    face_plan.flow_parameter = flow_parameter;
    BoundaryRelation momentum_relation = flow_relation(spec.flow_kind);
    BoundaryValueSource momentum_source = BoundaryValueSource::none;
    if (spec.flow_kind == BoundaryKind::velocity_inlet ||
        spec.flow_kind == BoundaryKind::moving_wall ||
        spec.flow_kind == BoundaryKind::no_slip_wall) {
      momentum_source = BoundaryValueSource::compiled_vector;
    } else if (face_plan.flow_kernel ==
                   BoundaryFlowKernel::mass_flow_constraint ||
               face_plan.flow_kernel == BoundaryFlowKernel::static_state ||
               face_plan.flow_kernel == BoundaryFlowKernel::total_state) {
      // Thermo-, integral-, backflow-, and characteristic-dependent targets
      // are resolved into field space by their owning runtime service.
      momentum_relation = BoundaryRelation::dirichlet;
      momentum_source = BoundaryValueSource::resolved_vector;
    } else if (spec.flow_kind == BoundaryKind::pressure_outlet &&
               spec.allow_backflow) {
      momentum_relation = BoundaryRelation::dirichlet;
      momentum_source = BoundaryValueSource::resolved_vector;
    }
    // NSCBC alone owns characteristic residuals. Its stencil ghost fill is
    // extrapolation and may not impose an additional U Dirichlet condition.
    status = add_span(plan, BoundaryStage::momentum, momentum_relation, face,
                      velocity, 0U, 3U, plan.required_ghost_width,
                      flow_parameter, momentum_source, patch);
    if (!status) {
      return status;
    }
    const bool resolved_pressure = resolved_pressure_ghost(spec.flow_kind);
    status = add_span(
        plan, BoundaryStage::pressure,
        resolved_pressure ? BoundaryRelation::dirichlet
                          : BoundaryRelation::zero_gradient,
        face, pressure, 0U, 1U, plan.required_ghost_width, flow_parameter,
        resolved_pressure ? BoundaryValueSource::resolved_scalar
                          : BoundaryValueSource::none,
        patch);
    if (!status) {
      return status;
    }
    std::uint32_t thermal_parameter = flow_parameter;
    if (spec.thermal_kind != BoundaryKind::none) {
      status = append_parameter(plan, spec, BoundaryParameterRole::thermal,
                                ScalarBoundaryKind::zero_gradient,
                                BoundaryScalarRole::none,
                                spec.temperature, distance_1, distance_2,
                                thermal_parameter);
      if (!status) {
        return status;
      }
    }
    face_plan.thermal_parameter = thermal_parameter;
    const bool flow_resolved_enthalpy =
        spec.flow_kind == BoundaryKind::velocity_inlet ||
        spec.flow_kind == BoundaryKind::mass_flow_inlet ||
        spec.flow_kind == BoundaryKind::static_state_inlet ||
        spec.flow_kind == BoundaryKind::total_state_inlet ||
        (spec.flow_kind == BoundaryKind::pressure_outlet &&
         spec.allow_backflow);
    const BoundaryRelation thermal_relation =
        flow_resolved_enthalpy
            ? BoundaryRelation::dirichlet
            : (spec.thermal_kind == BoundaryKind::isothermal_wall
            ? BoundaryRelation::dirichlet
            : (spec.thermal_kind == BoundaryKind::heat_flux_wall
                   ? BoundaryRelation::normal_gradient
                   : BoundaryRelation::zero_gradient));
    BoundaryValueSource thermal_source = BoundaryValueSource::none;
    if (thermal_relation == BoundaryRelation::dirichlet) {
      thermal_source = BoundaryValueSource::resolved_scalar;
    } else if (thermal_relation == BoundaryRelation::normal_gradient) {
      thermal_source = BoundaryValueSource::resolved_normal_gradient;
    }
    status = add_span(plan, BoundaryStage::enthalpy, thermal_relation, face,
                      enthalpy, 0U, 1U, plan.required_ghost_width,
                      thermal_parameter, thermal_source, patch);
    if (!status) {
      return status;
    }
    if (plan.scalar_targets.size() >
        static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
      return {StatusCode::invalid_plan, kBoundaryOverflow};
    }
    face_plan.scalar_begin =
        static_cast<std::uint32_t>(plan.scalar_targets.size());
    for (const ScalarBoundarySpec& scalar : spec.scalars) {
      const auto found = std::lower_bound(
          scalar_fields.begin(), scalar_fields.end(), scalar.stable_name,
          [](const auto& entry, const std::string& name) {
            return entry.first < name;
          });
      if (found == scalar_fields.end() ||
          found->first != scalar.stable_name) {
        return {StatusCode::invalid_plan, kBoundaryScalar};
      }
      std::uint32_t parameter = 0U;
      const auto declared = std::find_if(
          model.transported_scalars.begin(),
          model.transported_scalars.end(),
          [&scalar](const TransportedScalarSpec& candidate) {
            return candidate.stable_name == scalar.stable_name;
          });
      if (declared == model.transported_scalars.end()) {
        return {StatusCode::invalid_plan, kBoundaryScalar};
      }
      const BoundaryScalarRole scalar_role =
          declared->role == TransportedScalarRole::species
              ? BoundaryScalarRole::species
              : BoundaryScalarRole::passive_scalar;
      status = append_parameter(
          plan, spec, BoundaryParameterRole::transported_scalar, scalar.kind,
          scalar_role, scalar.value, distance_1, distance_2, parameter);
      if (!status) {
        return status;
      }
      const bool conditional_backflow =
          (spec.flow_kind == BoundaryKind::pressure_outlet ||
           spec.flow_kind == BoundaryKind::nscbc_outlet) &&
          spec.allow_backflow;
      plan.scalar_backflow_targets[parameter] = scalar.backflow_value;
      const BoundaryRelation relation =
          conditional_backflow ? BoundaryRelation::convective
                               : scalar_relation(scalar.kind);
      BoundaryValueSource source = BoundaryValueSource::none;
      if (conditional_backflow) {
        source = BoundaryValueSource::resolved_scalar;
      } else if (relation == BoundaryRelation::dirichlet) {
        source = BoundaryValueSource::compiled_scalar;
      } else if (relation == BoundaryRelation::normal_gradient) {
        source = BoundaryValueSource::resolved_normal_gradient;
      } else if (relation == BoundaryRelation::convective) {
        source = BoundaryValueSource::resolved_scalar;
      }
      status = add_span(plan, BoundaryStage::scalar, relation, face,
                        found->second, 0U, 1U,
                        plan.required_ghost_width, parameter, source, patch);
      if (!status) {
        return status;
      }
    }
    face_plan.scalar_count = static_cast<std::uint32_t>(
        plan.scalar_targets.size() - face_plan.scalar_begin);
  }
  std::stable_sort(plan.spans.begin(), plan.spans.end(),
                   [](const BoundaryIndexSpan& left,
                      const BoundaryIndexSpan& right) {
                     return std::tie(left.stage, left.relation, left.face,
                                     left.field, left.component_begin) <
                            std::tie(right.stage, right.relation, right.face,
                                     right.field, right.component_begin);
                   });
  for (std::size_t begin = 0U; begin < plan.spans.size();) {
    std::size_t end = begin + 1U;
    while (end < plan.spans.size() &&
           plan.spans[end].stage == plan.spans[begin].stage &&
           plan.spans[end].relation == plan.spans[begin].relation) {
      ++end;
    }
    if (begin > std::numeric_limits<std::uint32_t>::max() ||
        end - begin > std::numeric_limits<std::uint32_t>::max()) {
      return {StatusCode::invalid_plan, kBoundaryOverflow};
    }
    plan.batches.push_back(BoundaryKernelBatch{
        plan.spans[begin].stage, plan.spans[begin].relation,
        static_cast<std::uint32_t>(begin),
        static_cast<std::uint32_t>(end - begin)});
    begin = end;
  }
  return {};
}

void hash_real3(std::uint64_t& hash, Real3 value) noexcept {
  hash = hash_mix(hash, double_bits(value.x));
  hash = hash_mix(hash, double_bits(value.y));
  hash = hash_mix(hash, double_bits(value.z));
}

void hash_text(std::uint64_t& hash, std::string_view value) noexcept {
  hash = hash_mix(hash, value.size());
  for (const unsigned char character : value) {
    hash = hash_mix(hash, character);
  }
}

std::uint64_t semantic_hash(const ValidatedModel& model,
                            PlanFingerprint geometry,
                            PlanFingerprint registry,
                            FieldId velocity, FieldId pressure,
                            FieldId enthalpy,
                            const std::vector<std::pair<std::string, FieldId>>&
                                scalars) noexcept {
  std::uint64_t hash = 1469598103934665603ULL;
  // Never trust the cached model fingerprint as the only contract: public
  // callers may construct a ValidatedModel directly. Hash every numerical
  // authority that can alter a compiled kernel or runtime target.
  hash = hash_mix(hash, model.fingerprint);
  hash = hash_mix(hash, geometry);
  hash = hash_mix(hash, registry);
  hash = hash_mix(hash, velocity);
  hash = hash_mix(hash, pressure);
  hash = hash_mix(hash, enthalpy);
  hash = hash_mix(hash, static_cast<std::uint8_t>(model.pressure_reference));
  hash = hash_mix(hash, model.transported_scalars.size());
  for (const TransportedScalarSpec& scalar : model.transported_scalars) {
    hash_text(hash, scalar.stable_name);
    hash = hash_mix(hash, static_cast<std::uint8_t>(scalar.role));
  }
  hash = hash_mix(hash, static_cast<std::uint8_t>(model.schemes.momentum));
  hash = hash_mix(hash, static_cast<std::uint8_t>(model.schemes.enthalpy));
  hash = hash_mix(hash, static_cast<std::uint8_t>(model.schemes.species));
  hash =
      hash_mix(hash, static_cast<std::uint8_t>(model.schemes.passive_scalar));
  hash = hash_mix(hash, static_cast<std::uint8_t>(model.schemes.diffusion));
  hash = hash_mix(hash, double_bits(model.schemes.limiter));
  hash = hash_mix(hash, static_cast<std::uint8_t>(model.time.control));
  hash = hash_mix(hash, static_cast<std::uint8_t>(model.time.scheme));
  hash = hash_mix(hash, double_bits(model.time.initial_dt));
  hash = hash_mix(hash, double_bits(model.time.minimum_dt));
  hash = hash_mix(hash, double_bits(model.time.maximum_dt));
  hash = hash_mix(hash, double_bits(model.time.convective_cfl));
  hash = hash_mix(hash, double_bits(model.time.viscous_cfl));
  hash = hash_mix(hash, double_bits(model.time.thermal_cfl));
  hash = hash_mix(hash, double_bits(model.time.species_cfl));
  hash = hash_mix(hash, double_bits(model.time.acoustic_cfl));
  hash = hash_mix(hash, double_bits(model.time.maximum_growth));
  hash = hash_mix(hash, double_bits(model.time.retry_factor));
  hash = hash_mix(hash, model.time.maximum_retries);
  hash = hash_mix(hash, double_bits(model.time.minimum_bdf_ratio));
  hash = hash_mix(hash, double_bits(model.time.maximum_bdf_ratio));
  for (const BoundaryFaceSpec& face : model.boundaries) {
    hash = hash_mix(hash, static_cast<std::uint8_t>(face.flow_kind));
    hash = hash_mix(hash, static_cast<std::uint8_t>(face.thermal_kind));
    hash_real3(hash, face.velocity);
    hash_real3(hash, face.direction);
    hash_real3(hash, face.backflow_velocity);
    hash = hash_mix(hash, double_bits(face.mass_flow_rate));
    hash = hash_mix(hash, double_bits(face.pressure));
    hash = hash_mix(hash, double_bits(face.temperature));
    hash = hash_mix(hash, double_bits(face.total_pressure));
    hash = hash_mix(hash, double_bits(face.total_temperature));
    hash = hash_mix(hash, double_bits(face.backflow_temperature));
    hash = hash_mix(hash, double_bits(face.heat_flux));
    hash = hash_mix(hash, double_bits(face.relaxation));
    hash = hash_mix(hash, double_bits(face.mach_limit));
    hash = hash_mix(hash, face.allow_backflow ? 1U : 0U);
    hash = hash_mix(hash, face.scalars.size());
    for (const ScalarBoundarySpec& scalar : face.scalars) {
      hash_text(hash, scalar.stable_name);
      hash = hash_mix(hash, static_cast<std::uint8_t>(scalar.kind));
      hash = hash_mix(hash, double_bits(scalar.value));
      hash = hash_mix(hash, static_cast<std::uint8_t>(scalar.backflow_kind));
      hash = hash_mix(hash, double_bits(scalar.backflow_value));
    }
  }
  for (const auto& scalar : scalars) {
    hash = hash_mix(hash, scalar.second);
    hash_text(hash, scalar.first);
  }
  return hash == 0U ? 1U : hash;
}

std::uint64_t local_hash(const BoundaryBuildData& plan,
                         const MeshPatch& patch) noexcept {
  std::uint64_t hash = plan.semantic_fingerprint;
  hash = hash_mix(hash, static_cast<std::uint32_t>(patch.begin.x));
  hash = hash_mix(hash, static_cast<std::uint32_t>(patch.begin.y));
  hash = hash_mix(hash, static_cast<std::uint32_t>(patch.begin.z));
  hash = hash_mix(hash, static_cast<std::uint32_t>(patch.cells.x));
  hash = hash_mix(hash, static_cast<std::uint32_t>(patch.cells.y));
  hash = hash_mix(hash, static_cast<std::uint32_t>(patch.cells.z));
  hash = hash_mix(hash, plan.spans.size());
  for (const BoundaryIndexSpan& span : plan.spans) {
    hash = hash_mix(hash, static_cast<std::uint8_t>(span.stage));
    hash = hash_mix(hash, static_cast<std::uint8_t>(span.relation));
    hash = hash_mix(hash, static_cast<std::uint8_t>(span.face));
    hash = hash_mix(hash, span.field);
    hash = hash_mix(hash, span.tangent_inner_count);
    hash = hash_mix(hash, span.tangent_outer_count);
    hash = hash_mix(hash, span.resolved_begin);
    hash = hash_mix(hash, static_cast<std::uint8_t>(span.value_source));
    hash = hash_mix(hash, span.resolved_stride);
  }
  return hash == 0U ? 1U : hash;
}

}  // namespace

Status detail::boundary_identity_for_registry(const ValidatedModel& model,
                                              const BoundaryPlan& boundary,
                                              PlanFingerprint registry,
                                              PlanFingerprint& out) noexcept
    try {
  const auto fields = boundary.transported_fields();
  if (registry == 0U || boundary.semantic_fingerprint() == 0U ||
      fields.size != model.transported_scalars.size())
    return {StatusCode::invalid_plan, kBoundaryRegistry};
  std::vector<std::pair<std::string, FieldId>> scalars;
  scalars.reserve(fields.size);
  for (std::size_t i = 0U; i < fields.size; ++i)
    scalars.emplace_back(model.transported_scalars[i].stable_name,
                         fields.data[i].field);
  std::sort(scalars.begin(), scalars.end());
  out = semantic_hash(model, boundary.geometry_fingerprint(), registry,
                      boundary.velocity_field(), boundary.pressure_field(),
                      boundary.enthalpy_field(), scalars);
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kBoundaryRegistry};
} catch (...) {
  return {StatusCode::invalid_plan, kBoundaryRegistry};
}

Status BoundaryCompiler::compile(MPI_Comm communicator,
                                 const ValidatedModel& model,
                                 const CartesianGeometryPlan& geometry,
                                 const MeshPatch& patch,
                                 FieldRegistry& registry,
                                 BoundaryPlan& boundary, SchemePlan& schemes,
                                 TimeSchemePlan& time,
                                 BoundaryCompileDiagnostics* diagnostics) noexcept {
  if (diagnostics != nullptr) {
    diagnostics->lowest_failing_rank = -1;
  }
  if (communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, kBoundaryInput};
  }
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kBoundaryCollective};
  }
  Status local = validate_geometry_patch(
      geometry, patch, rank, size, required_ghost_width(model.schemes));
  if (local) {
    local = validate_model(model);
  }
  int lowest = -1;
  Status consensus =
      collective_status(communicator, rank, size, local, lowest);
  if (!consensus) {
    if (diagnostics != nullptr) {
      diagnostics->lowest_failing_rank = lowest;
    }
    return consensus;
  }

  try {
    FieldRegistry registry_candidate = registry;
    BoundaryPlan boundary_candidate;
    BoundaryBuildData boundary_build;
    SchemePlan scheme_candidate;
    TimeSchemePlan time_candidate;
    scheme_candidate.momentum_ = model.schemes.momentum;
    scheme_candidate.enthalpy_ = model.schemes.enthalpy;
    scheme_candidate.species_ = model.schemes.species;
    scheme_candidate.passive_scalar_ = model.schemes.passive_scalar;
    scheme_candidate.diffusion_ = model.schemes.diffusion;
    scheme_candidate.limiter_ = model.schemes.limiter;
    scheme_candidate.required_ghost_width_ =
        required_ghost_width(model.schemes);

    FieldId velocity = 0U;
    FieldId pressure = 0U;
    FieldId enthalpy = 0U;
    local = registry_candidate.require_field(
        "U", 3U, scheme_candidate.required_ghost_width_, velocity);
    if (local) {
      local = registry_candidate.require_field(
          "pi", 1U, scheme_candidate.required_ghost_width_, pressure);
    }
    if (local) {
      local = registry_candidate.require_field(
          "h", 1U, scheme_candidate.required_ghost_width_, enthalpy);
    }
    std::vector<std::string> scalar_names;
    if (local) {
      scalar_names.reserve(model.transported_scalars.size());
      for (const TransportedScalarSpec& scalar : model.transported_scalars) {
        scalar_names.push_back(scalar.stable_name);
      }
      std::sort(scalar_names.begin(), scalar_names.end());
    }
    std::vector<std::pair<std::string, FieldId>> scalar_fields;
    if (local) {
      scalar_fields.reserve(scalar_names.size());
      for (const std::string& name : scalar_names) {
        FieldId id = 0U;
        local = registry_candidate.require_field(
            name, 1U, scheme_candidate.required_ghost_width_, id);
        if (!local) {
          break;
        }
        scalar_fields.emplace_back(name, id);
      }
    }
    if (local) {
      local = TimeSchemePlan::compile(model.time, time_candidate);
    }
    if (local) {
      local = compile_local_plan(model, geometry, patch, velocity, pressure,
                                 enthalpy, scalar_fields, boundary_build);
    }
    consensus = collective_status(communicator, rank, size, local, lowest);
    if (!consensus) {
      if (diagnostics != nullptr) {
        diagnostics->lowest_failing_rank = lowest;
      }
      return consensus;
    }
    const std::uint64_t semantic = semantic_hash(
        model, geometry.fingerprint(), registry_candidate.fingerprint(),
        velocity, pressure, enthalpy, scalar_fields);
    std::uint64_t minimum_semantic = semantic;
    std::uint64_t maximum_semantic = semantic;
    if (MPI_Allreduce(MPI_IN_PLACE, &minimum_semantic, 1, MPI_UINT64_T,
                      MPI_MIN, communicator) != MPI_SUCCESS ||
        MPI_Allreduce(MPI_IN_PLACE, &maximum_semantic, 1, MPI_UINT64_T,
                      MPI_MAX, communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kBoundaryCollective};
    }
    if (minimum_semantic != maximum_semantic) {
      // The semantic comparison proves divergence but intentionally does not
      // mutate a previously published plan. Attribute it to the first rank
      // whose semantic differs from rank zero's authority.
      std::uint64_t authority = semantic;
      if (MPI_Bcast(&authority, 1, MPI_UINT64_T, 0, communicator) !=
          MPI_SUCCESS) {
        return {StatusCode::mpi_failure, kBoundaryCollective};
      }
      const int candidate = semantic == authority ? size : rank;
      int first = size;
      if (MPI_Allreduce(&candidate, &first, 1, MPI_INT, MPI_MIN,
                        communicator) != MPI_SUCCESS) {
        return {StatusCode::mpi_failure, kBoundaryCollective};
      }
      if (diagnostics != nullptr) {
        diagnostics->lowest_failing_rank = first < size ? first : -1;
      }
      return {StatusCode::invalid_plan, kBoundaryRegistry};
    }
    boundary_build.semantic_fingerprint = semantic;
    boundary_build.local_layout_fingerprint =
        local_hash(boundary_build, patch);
    boundary_build.revision = semantic;
    boundary_candidate.faces_ = std::move(boundary_build.faces);
    boundary_candidate.spans_ = std::move(boundary_build.spans);
    boundary_candidate.batches_ = std::move(boundary_build.batches);
    boundary_candidate.transported_fields_ =
        std::move(boundary_build.transported_fields);
    boundary_candidate.velocity_x_ = std::move(boundary_build.velocity_x);
    boundary_candidate.velocity_y_ = std::move(boundary_build.velocity_y);
    boundary_candidate.velocity_z_ = std::move(boundary_build.velocity_z);
    boundary_candidate.direction_x_ = std::move(boundary_build.direction_x);
    boundary_candidate.direction_y_ = std::move(boundary_build.direction_y);
    boundary_candidate.direction_z_ = std::move(boundary_build.direction_z);
    boundary_candidate.backflow_velocity_x_ =
        std::move(boundary_build.backflow_velocity_x);
    boundary_candidate.backflow_velocity_y_ =
        std::move(boundary_build.backflow_velocity_y);
    boundary_candidate.backflow_velocity_z_ =
        std::move(boundary_build.backflow_velocity_z);
    boundary_candidate.scalar_targets_ =
        std::move(boundary_build.scalar_targets);
    boundary_candidate.scalar_backflow_targets_ =
        std::move(boundary_build.scalar_backflow_targets);
    boundary_candidate.pressure_targets_ =
        std::move(boundary_build.pressure_targets);
    boundary_candidate.temperature_targets_ =
        std::move(boundary_build.temperature_targets);
    boundary_candidate.total_pressure_targets_ =
        std::move(boundary_build.total_pressure_targets);
    boundary_candidate.total_temperature_targets_ =
        std::move(boundary_build.total_temperature_targets);
    boundary_candidate.backflow_temperature_targets_ =
        std::move(boundary_build.backflow_temperature_targets);
    boundary_candidate.heat_flux_targets_ =
        std::move(boundary_build.heat_flux_targets);
    boundary_candidate.mass_flow_targets_ =
        std::move(boundary_build.mass_flow_targets);
    boundary_candidate.relaxation_rates_ =
        std::move(boundary_build.relaxation_rates);
    boundary_candidate.mach_limits_ = std::move(boundary_build.mach_limits);
    boundary_candidate.allow_backflow_ =
        std::move(boundary_build.allow_backflow);
    boundary_candidate.parameter_roles_ =
        std::move(boundary_build.parameter_roles);
    boundary_candidate.scalar_kinds_ =
        std::move(boundary_build.scalar_kinds);
    boundary_candidate.scalar_roles_ =
        std::move(boundary_build.scalar_roles);
    boundary_candidate.normal_distance_1_ =
        std::move(boundary_build.normal_distance_1);
    boundary_candidate.normal_distance_2_ =
        std::move(boundary_build.normal_distance_2);
    boundary_candidate.halo_topology_ = boundary_build.halo_topology;
    boundary_candidate.velocity_field_ = boundary_build.velocity_field;
    boundary_candidate.pressure_field_ = boundary_build.pressure_field;
    boundary_candidate.enthalpy_field_ = boundary_build.enthalpy_field;
    boundary_candidate.pressure_reference_ =
        boundary_build.pressure_reference;
    boundary_candidate.local_cells_ = boundary_build.local_cells;
    boundary_candidate.resolved_scalar_count_ =
        boundary_build.resolved_scalar_count;
    boundary_candidate.resolved_vector_count_ =
        boundary_build.resolved_vector_count;
    boundary_candidate.resolved_normal_gradient_count_ =
        boundary_build.resolved_normal_gradient_count;
    boundary_candidate.required_ghost_width_ =
        boundary_build.required_ghost_width;
    boundary_candidate.revision_ = boundary_build.revision;
    boundary_candidate.geometry_fingerprint_ = geometry.fingerprint();
    boundary_candidate.semantic_fingerprint_ =
        boundary_build.semantic_fingerprint;
    boundary_candidate.local_layout_fingerprint_ =
        boundary_build.local_layout_fingerprint;
    scheme_candidate.fingerprint_ = hash_mix(
        semantic, static_cast<std::uint8_t>(model.schemes.momentum));
    if (scheme_candidate.fingerprint_ == 0U) {
      scheme_candidate.fingerprint_ = 1U;
    }
    // Registry and all plans publish together only after every rank has the
    // same semantic FieldId signature and complete local compact layout.
    registry = std::move(registry_candidate);
    boundary = std::move(boundary_candidate);
    schemes = std::move(scheme_candidate);
    time = std::move(time_candidate);
    return {};
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, 0U};
  } catch (...) {
    local = {StatusCode::invalid_plan, kBoundaryInput};
  }
  const Status result =
      collective_status(communicator, rank, size, local, lowest);
  if (!result && diagnostics != nullptr) {
    diagnostics->lowest_failing_rank = lowest;
  }
  return result;
}

}  // namespace hundun::v04
