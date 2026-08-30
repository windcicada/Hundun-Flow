// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_boundary.hpp"

#include "field_view_interval_detail.hpp"

#include "hundun/v04_physics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kClosureInput = 736U;
constexpr std::uint32_t kClosureState = 737U;
constexpr std::size_t kMaximumIndependentSpecies = 64U;
constexpr std::size_t kOutputCount = 8U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

std::uint64_t nonzero(std::uint64_t value) noexcept {
  return value == 0U ? 1U : value;
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

std::uint64_t mix_int3(std::uint64_t hash, Int3 value) noexcept {
  hash = hash_mix(hash,
                  static_cast<std::uint32_t>(value.x));
  hash = hash_mix(hash,
                  static_cast<std::uint32_t>(value.y));
  return hash_mix(hash, static_cast<std::uint32_t>(value.z));
}

std::uint64_t mix_view_identity(std::uint64_t hash,
                                ConstFieldView view) noexcept {
  hash = hash_mix(
      hash, static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(
                view.base)));
  hash = mix_int3(hash, view.interior);
  hash = mix_int3(hash, view.ghosts);
  hash = hash_mix(hash, view.components);
  hash = hash_mix(hash, view.stride_y);
  hash = hash_mix(hash, view.stride_z);
  hash = hash_mix(hash, view.component_stride);
  hash = hash_mix(hash, view.replica);
  hash = hash_mix(hash, view.field);
  hash = hash_mix(hash, view.revision);
  hash = hash_mix(hash, view.storage_identity);
  return hash_mix(hash, view.revision_domain);
}

std::uint64_t mix_view_storage_lineage(std::uint64_t hash,
                                       ConstFieldView view) noexcept {
  hash = hash_mix(
      hash, static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(
                view.base)));
  hash = mix_int3(hash, view.interior);
  hash = mix_int3(hash, view.ghosts);
  hash = hash_mix(hash, view.components);
  hash = hash_mix(hash, view.stride_y);
  hash = hash_mix(hash, view.stride_z);
  hash = hash_mix(hash, view.component_stride);
  hash = hash_mix(hash, view.replica);
  hash = hash_mix(hash, view.field);
  hash = hash_mix(hash, view.storage_identity);
  return hash_mix(hash, view.revision_domain);
}

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

template <class T>
bool valid_scalar_view(BasicFieldView<T> view, Int3 interior,
                       Int3 ghosts) noexcept {
  detail::FieldStorageInterval interval{};
  return view.base != nullptr && same(view.interior, interior) &&
         same(view.ghosts, ghosts) && view.components == 1U &&
         view.revision != 0U && view.storage_identity != 0U &&
         view.revision_domain != 0U &&
         detail::field_storage_interval(view, interval);
}

bool matches(BoundaryGhostFieldAuthority authority,
             ConstFieldView view) noexcept {
  return authority.base == view.base && authority.replica == view.replica &&
         authority.field == view.field && authority.revision == view.revision &&
         authority.storage == view.storage_identity &&
         authority.revision_domain == view.revision_domain &&
         authority.base != nullptr && authority.revision != 0U &&
         authority.storage != 0U &&
         authority.revision_domain != 0U;
}

bool matches_legacy(BoundaryGhostFieldAuthority authority,
                    ConstFieldView view) noexcept {
  return authority.field == view.field && authority.revision == view.revision &&
         authority.storage == view.storage_identity &&
         authority.revision_domain == view.revision_domain &&
         authority.revision != 0U && authority.storage != 0U &&
         authority.revision_domain != 0U;
}

bool distinct_field(FieldId field, ConstFieldView pressure,
                    ConstFieldView enthalpy,
                    Span<const ConstFieldView> species) noexcept {
  if (field == pressure.field || field == enthalpy.field) return false;
  for (std::size_t index = 0U; index < species.size; ++index) {
    if (field == species.data[index].field) return false;
  }
  return true;
}

bool finite_positive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

bool physical_faces(const BoundaryPlan &boundary,
                    std::array<bool, 6U> &physical) noexcept {
  for (std::size_t face_index = 0U; face_index < physical.size();
       ++face_index) {
    const BoundaryFacePlan *face = nullptr;
    if (!boundary.face(static_cast<CartesianFace>(face_index), face) ||
        face == nullptr) {
      return false;
    }
    physical[face_index] = face->local_owner && !face->periodic;
  }
  return true;
}

template <class Visitor>
Status visit_box(Int3 begin, Int3 end, Visitor &visitor) noexcept {
  for (std::int32_t z = begin.z; z < end.z; ++z) {
    for (std::int32_t y = begin.y; y < end.y; ++y) {
      for (std::int32_t x = begin.x; x < end.x; ++x) {
        const Status status = visitor(Int3{x, y, z});
        if (!status) return status;
      }
    }
  }
  return {};
}

// Visit each locally owned physical ghost exactly once.  X shells own every
// admissible x corner, Y shells own only x-interior cells, and Z shells own
// only x/y-interior cells.  Thus MPI/periodic corners are excluded and the
// work is O(surface * reach), not O(padded volume).
template <class Visitor>
Status for_each_physical_ghost(
    Int3 interior, Int3 ghosts, const std::array<bool, 6U> &physical,
    Visitor &&visitor) noexcept {
  const std::int32_t y_begin = physical[2U] ? -ghosts.y : 0;
  const std::int32_t y_end =
      physical[3U] ? interior.y + ghosts.y : interior.y;
  const std::int32_t z_begin = physical[4U] ? -ghosts.z : 0;
  const std::int32_t z_end =
      physical[5U] ? interior.z + ghosts.z : interior.z;
  if (physical[0U]) {
    const Status status = visit_box(
        {-ghosts.x, y_begin, z_begin}, {0, y_end, z_end}, visitor);
    if (!status) return status;
  }
  if (physical[1U]) {
    const Status status = visit_box(
        {interior.x, y_begin, z_begin},
        {interior.x + ghosts.x, y_end, z_end}, visitor);
    if (!status) return status;
  }
  if (physical[2U]) {
    const Status status = visit_box(
        {0, -ghosts.y, z_begin}, {interior.x, 0, z_end}, visitor);
    if (!status) return status;
  }
  if (physical[3U]) {
    const Status status = visit_box(
        {0, interior.y, z_begin},
        {interior.x, interior.y + ghosts.y, z_end}, visitor);
    if (!status) return status;
  }
  if (physical[4U]) {
    const Status status = visit_box(
        {0, 0, -ghosts.z}, {interior.x, interior.y, 0}, visitor);
    if (!status) return status;
  }
  if (physical[5U]) {
    const Status status = visit_box(
        {0, 0, interior.z},
        {interior.x, interior.y, interior.z + ghosts.z}, visitor);
    if (!status) return status;
  }
  return {};
}

struct ThermophysicalGhostDigests {
  std::uint64_t primitive{kFnvOffset};
  std::uint64_t density{kFnvOffset};
};

void mix_ghost_digests(ThermophysicalGhostDigests &digests,
                       const BoundaryThermophysicalGhostBinding &binding,
                       Int3 cell) noexcept {
  digests.primitive = mix_int3(digests.primitive, cell);
  digests.primitive = hash_mix(
      digests.primitive,
      double_bits(binding.pressure_perturbation.unchecked(cell, 0U)));
  digests.primitive = hash_mix(
      digests.primitive, double_bits(binding.enthalpy.unchecked(cell, 0U)));
  for (std::size_t species = 0U;
       species < binding.independent_species.size; ++species) {
    digests.primitive = hash_mix(
        digests.primitive,
        double_bits(
            binding.independent_species.data[species].unchecked(cell, 0U)));
  }
  digests.density = mix_int3(digests.density, cell);
  digests.density = hash_mix(
      digests.density, double_bits(binding.density.unchecked(cell, 0U)));
}

ThermophysicalGhostDigests thermophysical_ghost_digests(
    const BoundaryThermophysicalGhostBinding &binding, Int3 interior,
    Int3 ghosts, const std::array<bool, 6U> &physical) noexcept {
  ThermophysicalGhostDigests digests;
  const Status visited = for_each_physical_ghost(
      interior, ghosts, physical, [&](Int3 cell) noexcept {
        mix_ghost_digests(digests, binding, cell);
        return Status{};
      });
  (void)visited;
  digests.primitive = nonzero(digests.primitive);
  digests.density = nonzero(digests.density);
  return digests;
}

std::uint64_t consumer_binding(
    const BoundaryPlan &boundary,
    const BoundaryThermophysicalGhostBinding &binding) noexcept {
  std::uint64_t hash = kFnvOffset;
  hash = hash_mix(hash, boundary.local_layout_fingerprint());
  hash = hash_mix(hash, double_bits(binding.pressure_reference));
  hash = mix_view_identity(hash, binding.pressure_perturbation);
  hash = mix_view_identity(hash, binding.enthalpy);
  hash = hash_mix(hash, binding.independent_species.size);
  for (std::size_t species = 0U;
       species < binding.independent_species.size; ++species) {
    hash = mix_view_identity(hash,
                             binding.independent_species.data[species]);
  }
  hash = mix_view_identity(hash, binding.density);
  return nonzero(hash);
}

Status evaluate_cell(
    const ThermodynamicsPlan &thermodynamics, const TransportPlan &transport,
    const BoundaryThermophysicalGhostInput &input, Int3 cell,
    std::array<double, kMaximumIndependentSpecies> &composition,
    std::array<double, kOutputCount> &values) noexcept {
  for (std::size_t species = 0U; species < input.independent_species.size;
       ++species) {
    composition[species] =
        input.independent_species.data[species].unchecked(cell, 0U);
  }
  const double perturbation = input.pressure_perturbation.unchecked(cell, 0U);
  const double enthalpy = input.enthalpy.unchecked(cell, 0U);
  const double pressure_absolute = input.pressure_reference + perturbation;
  if (!std::isfinite(input.pressure_reference) ||
      !std::isfinite(perturbation) || !finite_positive(pressure_absolute) ||
      !std::isfinite(enthalpy)) {
    return {StatusCode::numerical_failure, kClosureState};
  }

  const Span<const double> fractions{composition.data(),
                                     input.independent_species.size};
  ThermoState thermo;
  Status status = thermodynamics.evaluate(pressure_absolute, enthalpy,
                                          fractions, Real3{}, thermo);
  if (!status) return status;
  MolecularTransportState molecular;
  status = transport.evaluate(thermo.temperature, fractions, molecular);
  if (!status) return status;
  const double conductivity_over_cp = molecular.conductivity / thermo.cp;
  if (!finite_positive(thermo.rho) || !finite_positive(thermo.temperature) ||
      !finite_positive(thermo.cp) || !finite_positive(thermo.drho_dp_hY) ||
      !std::isfinite(thermo.drho_dh_pY) || thermo.drho_dh_pY >= 0.0 ||
      !finite_positive(molecular.viscosity) ||
      !finite_positive(molecular.conductivity) ||
      !finite_positive(conductivity_over_cp)) {
    return {StatusCode::numerical_failure, kClosureState};
  }
  values = {thermo.rho,
            thermo.temperature,
            thermo.cp,
            thermo.drho_dp_hY,
            thermo.drho_dh_pY,
            molecular.viscosity,
            molecular.conductivity,
            conductivity_over_cp};
  return {};
}

}  // namespace

bool BoundaryThermophysicalGhostCertificate::valid() const noexcept {
  return collective_semantics_ != 0U && collective_lineage_ != 0U &&
         collective_target_ != 0U && rank_local_lineage_ != 0U &&
         rank_local_binding_ != 0U && consumer_binding_ != 0U &&
         boundary_semantic_ != 0U && boundary_layout_ != 0U &&
         thermodynamics_ != 0U && transport_ != 0U &&
         boundary_revision_ != 0U && target_time_ != 0U && geometry_ != 0U &&
         pressure_reference_ != 0U && numeric_boundary_ != 0U &&
         primitive_digest_ != 0U && density_digest_ != 0U &&
         phase_ != BoundaryThermophysicalGhostPhase::invalid;
}

bool BoundaryThermophysicalGhostCertificate::matches(
    const BoundaryPlan &boundary, BoundaryThermophysicalGhostContext context,
    const BoundaryThermophysicalGhostBinding &binding) const noexcept {
  if (!valid() || !context.valid() ||
      context.numeric_boundary != boundary.revision() ||
      context.geometry != boundary.geometry_fingerprint() ||
      boundary.semantic_fingerprint() != boundary_semantic_ ||
      boundary.local_layout_fingerprint() != boundary_layout_ ||
      boundary.revision() != boundary_revision_ ||
      context.target_time != target_time_ || context.geometry != geometry_ ||
      context.pressure_reference != pressure_reference_ ||
      context.numeric_boundary != numeric_boundary_ || context.phase != phase_ ||
      double_bits(binding.pressure_reference) != pressure_reference_bits_ ||
      binding.pressure_perturbation.field != boundary.pressure_field() ||
      binding.enthalpy.field != boundary.enthalpy_field() ||
      binding.independent_species.size > kMaximumIndependentSpecies ||
      (binding.independent_species.size != 0U &&
       binding.independent_species.data == nullptr)) {
    return false;
  }
  const Int3 interior = boundary.local_cells();
  const Int3 ghosts = binding.pressure_perturbation.ghosts;
  if (!valid_scalar_view(binding.pressure_perturbation, interior, ghosts) ||
      !valid_scalar_view(binding.enthalpy, interior, ghosts) ||
      !valid_scalar_view(binding.density, interior, ghosts)) {
    return false;
  }
  for (std::size_t species = 0U;
       species < binding.independent_species.size; ++species) {
    if (!valid_scalar_view(binding.independent_species.data[species], interior,
                           ghosts)) {
      return false;
    }
  }
  std::array<bool, 6U> physical{};
  if (!physical_faces(boundary, physical)) return false;
  const ThermophysicalGhostDigests digests = thermophysical_ghost_digests(
      binding, interior, ghosts, physical);
  return consumer_binding(boundary, binding) == consumer_binding_ &&
         digests.primitive == primitive_digest_ &&
         digests.density == density_digest_;
}

Status BoundaryThermophysicalFaceClosure::close(
    const BoundaryPlan &boundary, const ThermodynamicsPlan &thermodynamics,
    const TransportPlan &transport,
    const BoundaryThermophysicalGhostInput &input,
    const BoundaryThermophysicalGhostOutput &output) noexcept {
  constexpr std::size_t kMaximumAuthorityFields =
      kMaximumIndependentSpecies + 2U;
  const BoundaryThermophysicalGhostAuthority authority = input.authority;
  if (!authority.valid() ||
      input.independent_species.size > kMaximumIndependentSpecies ||
      (input.independent_species.size != 0U &&
       input.independent_species.data == nullptr) ||
      authority.fields.size != input.independent_species.size + 2U ||
      !matches_legacy(authority.fields.data[0U],
                      input.pressure_perturbation) ||
      !matches_legacy(authority.fields.data[1U], input.enthalpy)) {
    return {StatusCode::invalid_plan, kClosureInput};
  }
  std::array<BoundaryGhostFieldAuthority, kMaximumAuthorityFields>
      exact_authorities{};
  exact_authorities[0U] =
      make_boundary_ghost_field_authority(input.pressure_perturbation);
  exact_authorities[1U] =
      make_boundary_ghost_field_authority(input.enthalpy);
  for (std::size_t species = 0U; species < input.independent_species.size;
       ++species) {
    const ConstFieldView view = input.independent_species.data[species];
    if (!matches_legacy(authority.fields.data[species + 2U], view)) {
      return {StatusCode::invalid_plan, kClosureInput};
    }
    exact_authorities[species + 2U] =
        make_boundary_ghost_field_authority(view);
  }
  BoundaryThermophysicalGhostInput certified_input = input;
  certified_input.authority.fields = {
      exact_authorities.data(), input.independent_species.size + 2U};
  BoundaryThermophysicalGhostCertificate ignored;
  return close(boundary, thermodynamics, transport, certified_input, output,
               {input.pressure_perturbation.revision,
                boundary.geometry_fingerprint(), authority.producer,
                boundary.revision(),
                BoundaryThermophysicalGhostPhase::terminal},
               ignored);
}

Status BoundaryThermophysicalFaceClosure::close(
    const BoundaryPlan &boundary, const ThermodynamicsPlan &thermodynamics,
    const TransportPlan &transport,
    const BoundaryThermophysicalGhostInput &input,
    const BoundaryThermophysicalGhostOutput &output,
    BoundaryThermophysicalGhostContext context,
    BoundaryThermophysicalGhostCertificate &certificate) noexcept {
  certificate = {};
  const BoundaryThermophysicalGhostAuthority authority = input.authority;
  const Int3 interior = boundary.local_cells();
  const Int3 ghosts = authority.ghosts;
  if (boundary.semantic_fingerprint() == 0U ||
      boundary.local_layout_fingerprint() == 0U ||
      thermodynamics.fingerprint() == 0U || transport.fingerprint() == 0U ||
      !context.valid() || context.numeric_boundary != boundary.revision() ||
      context.geometry != boundary.geometry_fingerprint() ||
      !authority.valid() || authority.boundary != boundary.revision() ||
      authority.boundary_layout != boundary.local_layout_fingerprint() ||
      !same(authority.interior, interior) ||
      authority.fields.size != input.independent_species.size + 2U ||
      input.independent_species.size > kMaximumIndependentSpecies ||
      (input.independent_species.size != 0U &&
       input.independent_species.data == nullptr) ||
      input.independent_species.size !=
          thermodynamics.independent_species_count() ||
      input.pressure_perturbation.field != boundary.pressure_field() ||
      input.enthalpy.field != boundary.enthalpy_field() ||
      ghosts.x < static_cast<std::int32_t>(boundary.required_ghost_width()) ||
      ghosts.y < static_cast<std::int32_t>(boundary.required_ghost_width()) ||
      ghosts.z < static_cast<std::int32_t>(boundary.required_ghost_width()) ||
      ghosts.x > UINT8_MAX || ghosts.y > UINT8_MAX || ghosts.z > UINT8_MAX ||
      authority.reach <
          static_cast<std::uint8_t>(std::max({ghosts.x, ghosts.y, ghosts.z})) ||
      !valid_scalar_view(input.pressure_perturbation, interior, ghosts) ||
      !valid_scalar_view(input.enthalpy, interior, ghosts) ||
      !matches(authority.fields.data[0U], input.pressure_perturbation) ||
      !matches(authority.fields.data[1U], input.enthalpy)) {
    return {StatusCode::invalid_plan, kClosureInput};
  }

  std::size_t expected_species = 0U;
  const Span<const BoundaryTransportedField> transported =
      boundary.transported_fields();
  for (std::size_t field = 0U; field < transported.size; ++field) {
    if (transported.data[field].role != BoundaryScalarRole::species) continue;
    if (expected_species >= input.independent_species.size ||
        transported.data[field].field !=
            input.independent_species.data[expected_species].field) {
      return {StatusCode::invalid_plan, kClosureInput};
    }
    ++expected_species;
  }
  if (expected_species != input.independent_species.size) {
    return {StatusCode::invalid_plan, kClosureInput};
  }

  for (std::size_t species = 0U; species < input.independent_species.size;
       ++species) {
    const ConstFieldView view = input.independent_species.data[species];
    if (!valid_scalar_view(view, interior, ghosts) ||
        !matches(authority.fields.data[species + 2U], view)) {
      return {StatusCode::invalid_plan, kClosureInput};
    }
    if (view.field == input.pressure_perturbation.field ||
        view.field == input.enthalpy.field ||
        detail::field_views_overlap(view, input.pressure_perturbation) ||
        detail::field_views_overlap(view, input.enthalpy)) {
      return {StatusCode::invalid_plan, kClosureInput};
    }
    for (std::size_t prior = 0U; prior < species; ++prior) {
      if (view.field == input.independent_species.data[prior].field ||
          detail::field_views_overlap(view,
                                      input.independent_species.data[prior])) {
        return {StatusCode::invalid_plan, kClosureInput};
      }
    }
  }
  if (input.pressure_perturbation.field == input.enthalpy.field ||
      detail::field_views_overlap(input.pressure_perturbation,
                                  input.enthalpy)) {
    return {StatusCode::invalid_plan, kClosureInput};
  }

  const std::array<FieldView, kOutputCount> outputs{
      output.density,
      output.temperature,
      output.heat_capacity_cp,
      output.density_pressure_derivative,
      output.density_enthalpy_derivative,
      output.molecular_viscosity,
      output.thermal_conductivity,
      output.conductivity_over_cp};
  for (std::size_t index = 0U; index < outputs.size(); ++index) {
    const FieldView view = outputs[index];
    if (!valid_scalar_view(view, interior, ghosts) ||
        !distinct_field(view.field, input.pressure_perturbation, input.enthalpy,
                        input.independent_species) ||
        detail::field_views_overlap(as_const(view),
                                    input.pressure_perturbation) ||
        detail::field_views_overlap(as_const(view), input.enthalpy)) {
      return {StatusCode::invalid_plan, kClosureInput};
    }
    for (std::size_t species = 0U; species < input.independent_species.size;
         ++species) {
      if (detail::field_views_overlap(
              as_const(view), input.independent_species.data[species])) {
        return {StatusCode::invalid_plan, kClosureInput};
      }
    }
    for (std::size_t prior = 0U; prior < index; ++prior) {
      if (view.field == outputs[prior].field ||
          detail::field_views_overlap(as_const(view),
                                      as_const(outputs[prior]))) {
        return {StatusCode::invalid_plan, kClosureInput};
      }
    }
  }

  std::array<bool, 6U> physical{};
  if (!physical_faces(boundary, physical))
    return {StatusCode::invalid_plan, kClosureInput};

  std::array<double, kMaximumIndependentSpecies> composition{};
  std::array<double, kOutputCount> values{};
  Status pass = for_each_physical_ghost(
      interior, ghosts, physical, [&](Int3 cell) noexcept {
        return evaluate_cell(thermodynamics, transport, input, cell,
                             composition, values);
      });
  if (!pass) return pass;
  const BoundaryThermophysicalGhostBinding binding{
      input.pressure_reference, input.pressure_perturbation, input.enthalpy,
      input.independent_species, as_const(output.density)};
  ThermophysicalGhostDigests digests;
  pass = for_each_physical_ghost(
      interior, ghosts, physical, [&](Int3 cell) noexcept {
        const Status status = evaluate_cell(thermodynamics, transport, input,
                                            cell, composition, values);
        if (!status) return status;
        for (std::size_t field = 0U; field < outputs.size(); ++field) {
          outputs[field].unchecked(cell, 0U) = values[field];
        }
        mix_ghost_digests(digests, binding, cell);
        return Status{};
      });
  if (!pass) return pass;
  const std::uint64_t primitive = nonzero(digests.primitive);
  const std::uint64_t density = nonzero(digests.density);
  std::uint64_t semantics = kFnvOffset;
  semantics = hash_mix(semantics, UINT64_C(0x425447484f535431));
  semantics = hash_mix(semantics, boundary.semantic_fingerprint());
  semantics = hash_mix(semantics, thermodynamics.fingerprint());
  semantics = hash_mix(semantics, transport.fingerprint());
  semantics = hash_mix(semantics, boundary.pressure_field());
  semantics = hash_mix(semantics, boundary.enthalpy_field());
  semantics = hash_mix(semantics, input.independent_species.size);
  semantics = hash_mix(semantics, boundary.required_ghost_width());
  for (std::size_t species = 0U;
       species < input.independent_species.size; ++species) {
    semantics =
        hash_mix(semantics, input.independent_species.data[species].field);
  }
  semantics = nonzero(semantics);
  std::uint64_t lineage = kFnvOffset;
  lineage = hash_mix(lineage, semantics);
  lineage = hash_mix(lineage, boundary.revision());
  lineage = hash_mix(lineage, context.target_time);
  lineage = hash_mix(lineage, context.geometry);
  // PressureReferenceCertificate::pressure_reference includes the
  // rank-local predictor state.  Its exact value belongs to the local
  // binding below; the collective target certifies only its presence and the
  // globally shared binary64 pressure_reference value.
  lineage = hash_mix(lineage, context.pressure_reference != 0U);
  lineage = hash_mix(lineage, context.numeric_boundary);
  lineage = hash_mix(lineage, double_bits(input.pressure_reference));
  lineage = nonzero(lineage);
  std::uint64_t target = kFnvOffset;
  target = hash_mix(target, lineage);
  target = hash_mix(target, static_cast<std::uint8_t>(context.phase));
  target = nonzero(target);
  const std::uint64_t consumer = consumer_binding(boundary, binding);
  std::uint64_t local_lineage = kFnvOffset;
  local_lineage = hash_mix(local_lineage, lineage);
  local_lineage =
      hash_mix(local_lineage, boundary.local_layout_fingerprint());
  local_lineage = mix_int3(local_lineage, authority.interior);
  local_lineage = mix_int3(local_lineage, authority.ghosts);
  local_lineage = hash_mix(local_lineage, authority.reach);
  local_lineage = hash_mix(local_lineage, context.pressure_reference);
  local_lineage =
      mix_view_storage_lineage(local_lineage, binding.pressure_perturbation);
  local_lineage =
      mix_view_storage_lineage(local_lineage, binding.enthalpy);
  for (std::size_t species = 0U;
       species < binding.independent_species.size; ++species) {
    local_lineage = mix_view_storage_lineage(
        local_lineage, binding.independent_species.data[species]);
  }
  local_lineage = mix_view_storage_lineage(local_lineage, binding.density);
  local_lineage = nonzero(local_lineage);
  std::uint64_t local = kFnvOffset;
  local = hash_mix(local, local_lineage);
  local = hash_mix(local, target);
  local = hash_mix(local, authority.producer);
  local = hash_mix(local, consumer);
  local = hash_mix(local, primitive);
  local = hash_mix(local, density);
  local = nonzero(local);
  BoundaryThermophysicalGhostCertificate candidate;
  candidate.collective_semantics_ = semantics;
  candidate.collective_lineage_ = lineage;
  candidate.collective_target_ = target;
  candidate.rank_local_lineage_ = local_lineage;
  candidate.rank_local_binding_ = local;
  candidate.consumer_binding_ = consumer;
  candidate.boundary_semantic_ = boundary.semantic_fingerprint();
  candidate.boundary_layout_ = boundary.local_layout_fingerprint();
  candidate.thermodynamics_ = thermodynamics.fingerprint();
  candidate.transport_ = transport.fingerprint();
  candidate.boundary_revision_ = boundary.revision();
  candidate.target_time_ = context.target_time;
  candidate.geometry_ = context.geometry;
  candidate.pressure_reference_ = context.pressure_reference;
  candidate.numeric_boundary_ = context.numeric_boundary;
  candidate.pressure_reference_bits_ = double_bits(input.pressure_reference);
  candidate.primitive_digest_ = primitive;
  candidate.density_digest_ = density;
  candidate.phase_ = context.phase;
  certificate = candidate;
  return {};
}

}  // namespace hundun::v04
