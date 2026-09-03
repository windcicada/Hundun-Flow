// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_equation_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kPressureEnergyThermo = 5781U;
constexpr std::uint32_t kPressureEnergyTemporal = 5785U;
constexpr std::uint32_t kPressureEnergySchurBinding = 5782U;
constexpr std::uint32_t kPressureEnergySchurPivot = 5783U;
constexpr std::uint32_t kPressureEnergySchurApply = 5784U;
constexpr std::uint32_t kPressureEnergyDiagonalBinding = 5786U;
constexpr std::uint32_t kPressureEnergyDiagonalApply = 5787U;
constexpr std::uint32_t kPressureEnergyPressureFluxBinding = 5788U;
constexpr std::uint32_t kPressureEnergyPressureFluxApply = 5789U;
constexpr std::uint32_t kPressureEnergyEnthalpyBinding = 5790U;
constexpr std::uint32_t kPressureEnergyEnthalpyApply = 5791U;
constexpr std::uint32_t kPressureEnergyGlobalization = 5792U;
constexpr std::uint64_t kPressureEnergySchurSchema =
    UINT64_C(0x7630347068657363);
constexpr std::uint64_t kPressureEnergyDiagonalSchema =
    UINT64_C(0x7630347068656467);
constexpr std::uint64_t kPressureEnergyPressureFluxSchema =
    UINT64_C(0x7630347065666c78);
constexpr std::uint64_t kPressureEnergyEnthalpySchema =
    UINT64_C(0x7630347065656e68);
// "v04pegl2": the joint Euclidean merit is a different policy from the
// original componentwise L-infinity globalization and must sign a distinct
// provenance lineage.
constexpr std::uint64_t kPressureEnergyGlobalizationSchema =
    UINT64_C(0x7630347065676c32);
// "v04peax1": a bounded alpha>1 exact candidate is deliberately distinct
// from the legacy alpha=1,1/2,... globalization ladder.
constexpr std::uint64_t kPressureEnergyExtrapolationSchema =
    UINT64_C(0x7630347065617831);
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value),
                "pressure-energy identity requires binary64");
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool same_bits(double left, double right) noexcept {
  return double_bits(left) == double_bits(right);
}

std::uint64_t nonzero_hash(std::uint64_t value) noexcept {
  return value == 0U ? UINT64_C(0x9e3779b97f4a7c15) : value;
}

void retain_first_failure(Status candidate, Status& deferred) noexcept {
  if (deferred && !candidate) deferred = candidate;
}

double pressure_energy_globalization_merit(
    const PressureEnergyGlobalizationSample& sample) noexcept {
  // A coupled Newton direction may exchange residual between the continuity
  // and energy blocks while reducing their joint norm.  A componentwise max
  // is nonsmooth at the active-block switch and falsely rejects that valid
  // descent.  hypot is the overflow-safe Euclidean merit for the two already
  // normalized residuals.
  return std::hypot(sample.global_normalized_continuity,
                    sample.global_normalized_energy);
}

PlanFingerprint pressure_energy_globalization_selection_provenance(
    const PressureEnergyGlobalizationSample& baseline,
    const PressureEnergyGlobalizationSample& selected,
    std::uint8_t selected_halvings) noexcept {
  std::uint64_t hash = hash_mix(kFnvOffset, kPressureEnergyGlobalizationSchema);
  hash = hash_mix(hash, baseline.corrector);
  hash = hash_mix(hash, baseline.target_time);
  hash = hash_mix(hash, baseline.correction_direction);
  hash = hash_mix(hash, baseline.state_provenance);
  hash = hash_mix(hash, baseline.mass_flux_provenance);
  hash = hash_mix(hash, selected_halvings);
  hash = hash_mix(hash, double_bits(selected.alpha));
  hash = hash_mix(hash, selected.state_provenance);
  hash = hash_mix(hash, selected.mass_flux_provenance);
  hash = hash_mix(hash,
                  double_bits(baseline.global_normalized_continuity));
  hash = hash_mix(hash, double_bits(baseline.global_normalized_energy));
  hash = hash_mix(hash,
                  double_bits(selected.global_normalized_continuity));
  hash = hash_mix(hash, double_bits(selected.global_normalized_energy));
  return nonzero_hash(hash);
}

PlanFingerprint pressure_energy_extrapolation_selection_provenance(
    const PressureEnergyGlobalizationSample& baseline,
    const PressureEnergyGlobalizationSample& selected) noexcept {
  std::uint64_t hash =
      hash_mix(kFnvOffset, kPressureEnergyExtrapolationSchema);
  hash = hash_mix(hash, baseline.corrector);
  hash = hash_mix(hash, baseline.target_time);
  hash = hash_mix(hash, baseline.correction_direction);
  hash = hash_mix(hash, baseline.state_provenance);
  hash = hash_mix(hash, baseline.mass_flux_provenance);
  hash = hash_mix(hash, double_bits(selected.alpha));
  hash = hash_mix(hash, selected.state_provenance);
  hash = hash_mix(hash, selected.mass_flux_provenance);
  hash = hash_mix(hash,
                  double_bits(baseline.global_normalized_continuity));
  hash = hash_mix(hash, double_bits(baseline.global_normalized_energy));
  hash = hash_mix(hash,
                  double_bits(selected.global_normalized_continuity));
  hash = hash_mix(hash, double_bits(selected.global_normalized_energy));
  return nonzero_hash(hash);
}

Status certify_pressure_energy_globalization_selection(
    const PressureEnergyGlobalizationSample& baseline,
    const PressureEnergyGlobalizationSample& selected,
    std::uint8_t selected_halvings, bool extrapolated,
    PressureEnergyGlobalizationSelectionCertificate& certificate) noexcept {
  const double baseline_merit = pressure_energy_globalization_merit(baseline);
  const double candidate_merit = pressure_energy_globalization_merit(selected);
  const double armijo_upper_bound =
      (1.0 - kPressureEnergyGlobalizationArmijoCoefficient * selected.alpha) *
      baseline_merit;
  const bool acceptable =
      selected.thermodynamically_admissible &&
      selected.state_and_flux_finite && std::isfinite(candidate_merit) &&
      selected.global_normalized_continuity >= 0.0 &&
      selected.global_normalized_energy >= 0.0 &&
      candidate_merit < baseline_merit &&
      candidate_merit <= armijo_upper_bound;
  if (!acceptable)
    return {StatusCode::rejected_step, kPressureEnergyGlobalization};

  PressureEnergyGlobalizationSelectionCertificate selected_certificate;
  selected_certificate.scope = PressureEnergyGlobalizationScope::
      frozen_momentum_continuity_energy_globalization;
  selected_certificate.alpha = selected.alpha;
  selected_certificate.baseline_normalized_continuity =
      baseline.global_normalized_continuity;
  selected_certificate.baseline_normalized_energy =
      baseline.global_normalized_energy;
  selected_certificate.candidate_normalized_continuity =
      selected.global_normalized_continuity;
  selected_certificate.candidate_normalized_energy =
      selected.global_normalized_energy;
  selected_certificate.baseline_merit = baseline_merit;
  selected_certificate.candidate_merit = candidate_merit;
  selected_certificate.armijo_upper_bound = armijo_upper_bound;
  selected_certificate.target_time = baseline.target_time;
  selected_certificate.correction_direction = baseline.correction_direction;
  selected_certificate.baseline_state_provenance =
      baseline.state_provenance;
  selected_certificate.baseline_mass_flux_provenance =
      baseline.mass_flux_provenance;
  selected_certificate.candidate_state_provenance =
      selected.state_provenance;
  selected_certificate.candidate_mass_flux_provenance =
      selected.mass_flux_provenance;
  selected_certificate.corrector = baseline.corrector;
  selected_certificate.selected_halvings = selected_halvings;
  selected_certificate.thermodynamically_admissible = true;
  selected_certificate.state_and_flux_finite = true;
  selected_certificate.strict_merit_decrease = true;
  selected_certificate.armijo_sufficient_decrease = true;
  selected_certificate.full_nonlinear_newton = false;
  selected_certificate.extrapolated = extrapolated;
  selected_certificate.selection_provenance =
      extrapolated
          ? pressure_energy_extrapolation_selection_provenance(baseline,
                                                                selected)
          : pressure_energy_globalization_selection_provenance(
                baseline, selected, selected_halvings);
  if (!selected_certificate.valid())
    return {StatusCode::invalid_plan, kPressureEnergyGlobalization};
  certificate = selected_certificate;
  return {};
}

bool same_shape(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool valid_identity(const LinearIdentity& identity) noexcept {
  return identity.symbolic != 0U && identity.numeric != 0U &&
         identity.hierarchy != 0U && identity.workspace != 0U &&
         identity.fingerprint != 0U;
}

bool same_identity(LinearIdentity left, LinearIdentity right) noexcept {
  return left.symbolic == right.symbolic && left.numeric == right.numeric &&
         left.hierarchy == right.hierarchy &&
         left.workspace == right.workspace &&
         left.fingerprint == right.fingerprint;
}

bool same_certificate(LinearOperatorCertificate left,
                      LinearOperatorCertificate right) noexcept {
  return same_identity(left.identity, right.identity) &&
         left.collective_fingerprint == right.collective_fingerprint &&
         same_shape(left.local_shape, right.local_shape) &&
         left.operator_class == right.operator_class;
}

std::size_t cell_count(Int3 shape) noexcept {
  return static_cast<std::size_t>(shape.x) *
         static_cast<std::size_t>(shape.y) *
         static_cast<std::size_t>(shape.z);
}

std::size_t cell_offset(Int3 shape, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(shape.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(shape.y) *
                  static_cast<std::size_t>(cell.z));
}

std::size_t face_offset(Int3 shape, Int3 face) noexcept {
  return static_cast<std::size_t>(face.x) +
         static_cast<std::size_t>(shape.x) *
             (static_cast<std::size_t>(face.y) +
              static_cast<std::size_t>(shape.y) *
                  static_cast<std::size_t>(face.z));
}

Int3 face_shape(Int3 cells, CartesianAxis axis) noexcept {
  if (axis == CartesianAxis::x) {
    ++cells.x;
  } else if (axis == CartesianAxis::y) {
    ++cells.y;
  } else {
    ++cells.z;
  }
  return cells;
}

ConstFaceFieldView select_face(ConstFaceFieldView x, ConstFaceFieldView y,
                               ConstFaceFieldView z,
                               CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x ? x : (axis == CartesianAxis::y ? y : z);
}

Span<const std::uint8_t> select_face_activity(
    PressureContinuityActivityView activity, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? activity.x_faces
             : (axis == CartesianAxis::y ? activity.y_faces : activity.z_faces);
}

bool active_face(PressureContinuityActivityView activity, CartesianAxis axis,
                 Int3 cells, Int3 face) noexcept {
  const Span<const std::uint8_t> selected =
      select_face_activity(activity, axis);
  return selected.size == 0U ||
         selected.data[face_offset(face_shape(cells, axis), face)] == 1U;
}

bool valid_face_view(ConstFaceFieldView view, CartesianAxis axis,
                     Int3 cells) noexcept {
  detail::FieldStorageInterval interval;
  return detail::face_storage_interval(view, interval) && view.axis == axis &&
         same_shape(view.extents, face_shape(cells, axis)) &&
         view.storage_identity != 0U && view.revision_domain != 0U;
}

bool valid_pressure_flux_activity(PressureContinuityActivityView activity,
                                  Int3 cells) noexcept {
  const bool empty = activity.cells.size == 0U && activity.x_faces.size == 0U &&
                     activity.y_faces.size == 0U && activity.z_faces.size == 0U;
  if (empty) {
    return activity.cells.data == nullptr && activity.x_faces.data == nullptr &&
           activity.y_faces.data == nullptr &&
           activity.z_faces.data == nullptr &&
           activity.local_fingerprint == 0U &&
           activity.collective_fingerprint == 0U;
  }
  const Int3 x_shape = face_shape(cells, CartesianAxis::x);
  const Int3 y_shape = face_shape(cells, CartesianAxis::y);
  const Int3 z_shape = face_shape(cells, CartesianAxis::z);
  const auto count = [](Int3 shape) noexcept {
    return static_cast<std::size_t>(shape.x) *
           static_cast<std::size_t>(shape.y) *
           static_cast<std::size_t>(shape.z);
  };
  if (activity.cells.data == nullptr ||
      activity.cells.size != cell_count(cells) ||
      activity.x_faces.data == nullptr ||
      activity.x_faces.size != count(x_shape) ||
      activity.y_faces.data == nullptr ||
      activity.y_faces.size != count(y_shape) ||
      activity.z_faces.data == nullptr ||
      activity.z_faces.size != count(z_shape) ||
      activity.local_fingerprint == 0U ||
      activity.collective_fingerprint == 0U) {
    return false;
  }
  const Span<const std::uint8_t> spans[]{activity.cells, activity.x_faces,
                                         activity.y_faces, activity.z_faces};
  for (Span<const std::uint8_t> span : spans) {
    for (std::size_t index = 0U; index < span.size; ++index) {
      if (span.data[index] > 1U) return false;
    }
  }
  // A solid control volume cannot retain a live interface in the conservative
  // pressure/energy block.  Fluid-fluid faces may still be disabled by a
  // stricter topology authority.
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (activity.cells.data[cell_offset(cells, cell)] != 0U) continue;
        const Int3 x_plus{x + 1, y, z};
        const Int3 y_plus{x, y + 1, z};
        const Int3 z_plus{x, y, z + 1};
        if (active_face(activity, CartesianAxis::x, cells, cell) ||
            active_face(activity, CartesianAxis::x, cells, x_plus) ||
            active_face(activity, CartesianAxis::y, cells, cell) ||
            active_face(activity, CartesianAxis::y, cells, y_plus) ||
            active_face(activity, CartesianAxis::z, cells, cell) ||
            active_face(activity, CartesianAxis::z, cells, z_plus)) {
          return false;
        }
      }
    }
  }
  return true;
}

bool valid_activity(PressureEnergyCellActivity activity,
                    Int3 shape) noexcept {
  if (activity.cells.size == 0U) {
    return activity.cells.data == nullptr &&
           activity.local_fingerprint == 0U &&
           activity.collective_fingerprint == 0U;
  }
  if (activity.cells.data == nullptr ||
      activity.cells.size != cell_count(shape) ||
      activity.local_fingerprint == 0U ||
      activity.collective_fingerprint == 0U) {
    return false;
  }
  for (std::size_t index = 0U; index < activity.cells.size; ++index) {
    if (activity.cells.data[index] > 1U) return false;
  }
  return true;
}

bool active_cell(PressureEnergyCellActivity activity, Int3 shape,
                 Int3 cell) noexcept {
  return activity.cells.size == 0U ||
         activity.cells.data[cell_offset(shape, cell)] == 1U;
}

std::uint64_t mix_identity(std::uint64_t hash,
                           LinearIdentity identity) noexcept {
  hash = hash_mix(hash, identity.symbolic);
  hash = hash_mix(hash, identity.hierarchy);
  // Numeric coefficients and workspace replicas are deliberately rank-local
  // identities.  Match the Krylov solve-contract convention: their presence
  // is collective, their values are not.  Symbolic and hierarchy authorities
  // remain exact collective values.
  hash = hash_mix(hash, identity.numeric != 0U);
  hash = hash_mix(hash, identity.workspace != 0U);
  return hash_mix(hash, identity.fingerprint != 0U);
}

std::uint64_t mix_collective_view(std::uint64_t hash,
                                  ConstFieldView view) noexcept {
  hash = hash_mix(hash, view.field);
  // Revision and replica are producer-local authorities.  The exact values
  // belong to the rank-local binding certificate; only a valid revision's
  // presence may enter the collective solve contract.
  return hash_mix(hash, view.revision != 0U);
}

std::uint64_t mix_local_face(std::uint64_t hash,
                             ConstFaceFieldView view) noexcept {
  hash = hash_mix(hash, view.storage_identity);
  hash = hash_mix(hash, view.revision_domain);
  hash = hash_mix(hash, reinterpret_cast<std::uintptr_t>(view.base));
  hash = hash_mix(hash, static_cast<std::uint32_t>(view.extents.x));
  hash = hash_mix(hash, static_cast<std::uint32_t>(view.extents.y));
  hash = hash_mix(hash, static_cast<std::uint32_t>(view.extents.z));
  hash = hash_mix(hash, view.stride_y);
  hash = hash_mix(hash, view.stride_z);
  return hash_mix(hash, static_cast<std::uint8_t>(view.axis));
}

std::uint64_t mix_local_field(std::uint64_t hash,
                              ConstFieldView view) noexcept {
  hash = hash_mix(hash, reinterpret_cast<std::uintptr_t>(view.base));
  hash = hash_mix(hash, static_cast<std::uint32_t>(view.interior.x));
  hash = hash_mix(hash, static_cast<std::uint32_t>(view.interior.y));
  hash = hash_mix(hash, static_cast<std::uint32_t>(view.interior.z));
  hash = hash_mix(hash, static_cast<std::uint32_t>(view.ghosts.x));
  hash = hash_mix(hash, static_cast<std::uint32_t>(view.ghosts.y));
  hash = hash_mix(hash, static_cast<std::uint32_t>(view.ghosts.z));
  hash = hash_mix(hash, view.components);
  hash = hash_mix(hash, view.stride_y);
  hash = hash_mix(hash, view.stride_z);
  hash = hash_mix(hash, view.component_stride);
  hash = hash_mix(hash, view.field);
  hash = hash_mix(hash, view.revision);
  hash = hash_mix(hash, view.storage_identity);
  hash = hash_mix(hash, view.revision_domain);
  return hash_mix(hash, view.replica);
}

bool valid_enthalpy_authority(
    const PressureEnergyEnthalpyAuthority& authority) noexcept {
  return detail::valid_bdf_coefficients(authority.bdf) &&
         authority.target_time != 0U && authority.geometry != 0U &&
         authority.numeric_boundary != 0U && authority.thermodynamics != 0U &&
         authority.transport != 0U && authority.equation_semantics != 0U &&
         authority.thermodynamics_semantics != 0U &&
         authority.transport_semantics != 0U;
}

PlanFingerprint enthalpy_collective_fingerprint(
    const PressureEnergyEnthalpyBinding& binding) noexcept {
  std::uint64_t hash = hash_mix(kFnvOffset, kPressureEnergyEnthalpySchema);
  hash = mix_identity(hash, binding.identity);
  hash = hash_mix(hash, static_cast<std::uint8_t>(binding.convection));
  hash =
      hash_mix(hash, static_cast<std::uint8_t>(binding.linearization_policy));
  hash = hash_mix(hash, double_bits(binding.authority.bdf.a0));
  hash = hash_mix(hash, double_bits(binding.authority.bdf.a1));
  hash = hash_mix(hash, double_bits(binding.authority.bdf.a2));
  hash = hash_mix(hash, binding.authority.bdf.order);
  hash = hash_mix(hash, binding.authority.target_time);
  hash = hash_mix(hash, binding.authority.geometry);
  hash = hash_mix(hash, binding.authority.numeric_boundary != 0U);
  hash = hash_mix(hash, binding.authority.thermodynamics != 0U);
  hash = hash_mix(hash, binding.authority.transport != 0U);
  hash = hash_mix(hash, binding.authority.equation_semantics);
  hash = hash_mix(hash, binding.authority.thermodynamics_semantics);
  hash = hash_mix(hash, binding.authority.transport_semantics);
  hash = hash_mix(hash, binding.geometry->fingerprint());
  hash = hash_mix(hash, binding.boundary->semantic_fingerprint());
  hash = hash_mix(hash, binding.convection_context.collective_semantics);
  hash = hash_mix(hash, binding.frozen_face_enthalpy.reconstruction);
  hash = hash_mix(hash, binding.target_flux.revision != 0U);
  const ConstFieldView fields[]{binding.assembled_diagonal,
                                binding.target_enthalpy,
                                binding.density_enthalpy_derivative,
                                binding.heat_capacity,
                                binding.thermal_conductivity,
                                binding.enthalpy_diffusivity,
                                as_const(binding.workspace.delta_temperature)};
  for (ConstFieldView field : fields)
    hash = mix_collective_view(hash, field);
  const bool compiled =
      binding.workspace.compiled.local_diagonal.base != nullptr ||
      binding.workspace.compiled.response_stage.base != nullptr ||
      binding.workspace.compiled.directional_branches.values.data != nullptr ||
      binding.workspace.compiled.thermal_conductance.x.base != nullptr ||
      binding.workspace.compiled.thermal_conductance.y.base != nullptr ||
      binding.workspace.compiled.thermal_conductance.z.base != nullptr;
  hash = hash_mix(hash, compiled ? 1U : 0U);
  if (compiled) {
    hash = mix_collective_view(
        hash, as_const(binding.workspace.compiled.local_diagonal));
    hash = mix_collective_view(
        hash, as_const(binding.workspace.compiled.response_stage));
  }
  hash = hash_mix(hash, binding.services.enthalpy_variation_field);
  hash = hash_mix(hash, binding.services.temperature_variation_field);
  hash = hash_mix(hash, binding.activity.collective_fingerprint);
  return nonzero_hash(hash);
}

RevisionToken enthalpy_binding_revision(
    const PressureEnergyEnthalpyBinding& binding,
    PlanFingerprint directional_branch_authority) noexcept {
  std::uint64_t hash = enthalpy_collective_fingerprint(binding);
  hash = hash_mix(hash, binding.geometry->topology_revision());
  hash = hash_mix(hash, binding.kernels->fingerprint());
  hash = hash_mix(hash, binding.boundary->revision());
  hash = hash_mix(hash, binding.boundary->local_layout_fingerprint());
  hash = hash_mix(hash, binding.convection_context.closure);
  hash = hash_mix(hash, binding.frozen_face_enthalpy.revision);
  hash = hash_mix(hash, binding.frozen_face_enthalpy.local_binding);
  hash = hash_mix(hash, binding.target_flux.revision);
  hash = hash_mix(hash, binding.services.halo->instance_identity());
  hash = hash_mix(hash, binding.services.halo_stage);
  hash = hash_mix(hash, binding.activity.local_fingerprint);
  hash = hash_mix(hash, directional_branch_authority);
  const ConstFieldView fields[]{binding.assembled_diagonal,
                                binding.target_enthalpy,
                                binding.density_enthalpy_derivative,
                                binding.heat_capacity,
                                binding.thermal_conductivity,
                                binding.enthalpy_diffusivity,
                                as_const(binding.workspace.delta_temperature)};
  for (ConstFieldView field : fields)
    hash = mix_local_field(hash, field);
  const bool compiled =
      binding.workspace.compiled.local_diagonal.base != nullptr ||
      binding.workspace.compiled.response_stage.base != nullptr ||
      binding.workspace.compiled.directional_branches.values.data != nullptr ||
      binding.workspace.compiled.thermal_conductance.x.base != nullptr ||
      binding.workspace.compiled.thermal_conductance.y.base != nullptr ||
      binding.workspace.compiled.thermal_conductance.z.base != nullptr;
  hash = hash_mix(hash, compiled ? 1U : 0U);
  if (compiled) {
    hash = mix_local_field(
        hash, as_const(binding.workspace.compiled.local_diagonal));
    hash = mix_local_field(
        hash, as_const(binding.workspace.compiled.response_stage));
    hash = hash_mix(
        hash, reinterpret_cast<std::uintptr_t>(
                  binding.workspace.compiled.directional_branches.values.data));
    hash = hash_mix(
        hash, binding.workspace.compiled.directional_branches.values.size);
  }
  const ConstFaceFieldView faces[]{
      binding.target_flux.x,
      binding.target_flux.y,
      binding.target_flux.z,
      binding.frozen_face_enthalpy.x,
      binding.frozen_face_enthalpy.y,
      binding.frozen_face_enthalpy.z,
      as_const(binding.workspace.directional_enthalpy.x),
      as_const(binding.workspace.directional_enthalpy.y),
      as_const(binding.workspace.directional_enthalpy.z)};
  for (ConstFaceFieldView face : faces)
    hash = mix_local_face(hash, face);
  if (compiled) {
    const ConstFaceFieldView compiled_faces[]{
        as_const(binding.workspace.compiled.thermal_conductance.x),
        as_const(binding.workspace.compiled.thermal_conductance.y),
        as_const(binding.workspace.compiled.thermal_conductance.z)};
    for (ConstFaceFieldView face : compiled_faces)
      hash = mix_local_face(hash, face);
  }
  return nonzero_hash(hash);
}

PlanFingerprint pressure_flux_collective_fingerprint(
    const PressureEnergyPressureFluxBinding& binding,
    const PressureCorrectionBoundaryCertificate& boundary) noexcept {
  std::uint64_t hash = hash_mix(kFnvOffset, kPressureEnergyPressureFluxSchema);
  hash = mix_identity(hash, binding.identity);
  hash = hash_mix(hash, boundary.semantic);
  hash = hash_mix(hash, binding.pressure.plan != 0U);
  hash = hash_mix(hash, binding.intermediate.corrector);
  // These authorities may contain rank-local storage in their value.  Their
  // exact values enter binding_revision below; only their collective presence
  // enters the Krylov certificate.
  hash = hash_mix(hash, binding.intermediate.pressure_face_coefficient != 0U);
  hash = hash_mix(hash, binding.intermediate.phi_h_by_a != 0U);
  hash = hash_mix(hash, binding.intermediate.dependency != 0U);
  hash = mix_collective_view(hash, binding.temporal_diagonal);
  hash = hash_mix(hash, binding.target_flux.revision != 0U);
  hash = hash_mix(hash, binding.frozen_face_enthalpy.revision != 0U);
  hash = hash_mix(hash, binding.frozen_face_enthalpy.reconstruction);
  hash = hash_mix(hash, binding.pressure.geometry);
  hash = hash_mix(hash, binding.pressure.numeric_boundary);
  hash = hash_mix(
      hash, binding.pressure.thermophysical_boundary_semantics);
  hash = hash_mix(
      hash, binding.pressure.thermophysical_boundary_collective_lineage);
  hash =
      hash_mix(hash, binding.pressure.thermophysical_boundary_target);
  hash = hash_mix(hash, binding.activity.collective_fingerprint);
  hash = hash_mix(hash, binding.pressure_work.valid());
  return nonzero_hash(hash);
}

RevisionToken pressure_flux_binding_revision(
    const PressureEnergyPressureFluxBinding& binding,
    const PressureCorrectionBoundaryCertificate& boundary) noexcept {
  std::uint64_t hash =
      pressure_flux_collective_fingerprint(binding, boundary);
  hash = hash_mix(hash, boundary.rank_local_layout);
  hash = hash_mix(hash, binding.pressure.plan);
  hash = hash_mix(hash, binding.intermediate.plan);
  hash = hash_mix(hash, binding.pressure.state);
  hash = hash_mix(hash, binding.intermediate.dependency);
  hash = hash_mix(hash, binding.intermediate.pressure_face_coefficient);
  hash = hash_mix(hash, binding.intermediate.phi_h_by_a);
  hash = hash_mix(
      hash,
      binding.pressure.thermophysical_boundary_rank_local_binding);
  hash = hash_mix(
      hash, binding.pressure.thermophysical_boundary_rank_local_lineage);
  hash = hash_mix(hash, binding.activity.local_fingerprint);
  hash = hash_mix(hash, binding.temporal_diagonal.storage_identity);
  hash = hash_mix(hash, binding.temporal_diagonal.revision_domain);
  hash = hash_mix(hash, binding.temporal_diagonal.revision);
  hash = hash_mix(hash, binding.temporal_diagonal.replica);
  hash = hash_mix(
      hash, reinterpret_cast<std::uintptr_t>(binding.temporal_diagonal.base));
  hash = hash_mix(hash, binding.target_flux.revision);
  hash = hash_mix(hash, binding.frozen_face_enthalpy.revision);
  hash = hash_mix(hash, binding.frozen_face_enthalpy.local_binding);
  const ConstFaceFieldView faces[]{
      binding.x_pressure_coefficient, binding.y_pressure_coefficient,
      binding.z_pressure_coefficient, binding.target_flux.x,
      binding.target_flux.y,          binding.target_flux.z,
      binding.frozen_face_enthalpy.x, binding.frozen_face_enthalpy.y,
      binding.frozen_face_enthalpy.z};
  for (ConstFaceFieldView face : faces) hash = mix_local_face(hash, face);
  return nonzero_hash(hash);
}

PlanFingerprint diagonal_collective_fingerprint(
    const PressureEnergyDiagonalBinding& binding) noexcept {
  std::uint64_t hash = hash_mix(kFnvOffset, kPressureEnergyDiagonalSchema);
  hash = mix_identity(hash, binding.identity);
  hash = mix_collective_view(hash, binding.diagonal);
  hash = hash_mix(hash, double_bits(binding.inactive_diagonal));
  hash = hash_mix(hash, binding.activity.collective_fingerprint);
  return nonzero_hash(hash);
}

PlanFingerprint schur_collective_fingerprint(
    const LinearOperatorCertificate& cp,
    const LinearOperatorCertificate& ep,
    const LinearOperatorCertificate& eh,
    const PressureEnergySchurBinding& binding) noexcept {
  std::uint64_t hash = hash_mix(kFnvOffset, kPressureEnergySchurSchema);
  hash = mix_identity(hash, cp.identity);
  hash = hash_mix(hash, cp.collective_fingerprint);
  hash = hash_mix(hash, ep.collective_fingerprint);
  hash = hash_mix(hash, eh.collective_fingerprint);
  hash = hash_mix(hash, static_cast<std::uint8_t>(cp.operator_class));
  hash = hash_mix(hash, static_cast<std::uint8_t>(ep.operator_class));
  hash = hash_mix(hash, static_cast<std::uint8_t>(eh.operator_class));
  hash = mix_collective_view(hash,
                             binding.continuity_enthalpy_diagonal);
  hash = mix_collective_view(hash,
                             binding.continuity_enthalpy_row_scale);
  hash = hash_mix(hash, binding.activity.collective_fingerprint);
  hash = hash_mix(hash, double_bits(binding.scaled_pivot_floor));
  return nonzero_hash(hash);
}

RevisionToken schur_block_revision(
    const LinearOperatorCertificate& cp,
    const LinearOperatorCertificate& ep,
    const LinearOperatorCertificate& eh,
    const PressureEnergySchurBinding& binding) noexcept {
  std::uint64_t hash = schur_collective_fingerprint(cp, ep, eh, binding);
  hash = hash_mix(hash, binding.activity.local_fingerprint);
  hash = hash_mix(hash,
                  binding.continuity_enthalpy_diagonal.storage_identity);
  hash = hash_mix(hash,
                  binding.continuity_enthalpy_diagonal.revision_domain);
  hash = hash_mix(hash, binding.continuity_enthalpy_diagonal.revision);
  hash = hash_mix(hash, binding.continuity_enthalpy_diagonal.replica);
  hash = hash_mix(hash,
                  binding.continuity_enthalpy_row_scale.storage_identity);
  hash = hash_mix(hash,
                  binding.continuity_enthalpy_row_scale.revision_domain);
  hash = hash_mix(hash, binding.continuity_enthalpy_row_scale.revision);
  hash = hash_mix(hash, binding.continuity_enthalpy_row_scale.replica);
  hash = hash_mix(
      hash, reinterpret_cast<std::uintptr_t>(
                binding.continuity_enthalpy_diagonal.base));
  hash = hash_mix(
      hash, reinterpret_cast<std::uintptr_t>(
                binding.continuity_enthalpy_row_scale.base));
  return nonzero_hash(hash);
}

template <class T>
bool valid_scalar_view(BasicFieldView<T> view, Int3 shape) noexcept {
  detail::FieldStorageInterval interval;
  return same_shape(view.interior, shape) && view.components == 1U &&
         view.revision != 0U && view.storage_identity != 0U &&
         view.revision_domain != 0U &&
         detail::field_storage_interval(view, interval);
}

std::size_t pressure_energy_face_count(Int3 cells,
                                       CartesianAxis axis) noexcept {
  const Int3 shape = face_shape(cells, axis);
  return static_cast<std::size_t>(shape.x) *
         static_cast<std::size_t>(shape.y) *
         static_cast<std::size_t>(shape.z);
}

std::size_t pressure_energy_branch_count(Int3 cells) noexcept {
  return pressure_energy_face_count(cells, CartesianAxis::x) +
         pressure_energy_face_count(cells, CartesianAxis::y) +
         pressure_energy_face_count(cells, CartesianAxis::z);
}

bool compiled_enthalpy_workspace_present(
    const PressureEnergyEnthalpyCompiledWorkspace& workspace) noexcept {
  return workspace.local_diagonal.base != nullptr ||
         workspace.response_stage.base != nullptr ||
         workspace.directional_branches.values.data != nullptr ||
         workspace.directional_branches.values.size != 0U ||
         workspace.thermal_conductance.x.base != nullptr ||
         workspace.thermal_conductance.y.base != nullptr ||
         workspace.thermal_conductance.z.base != nullptr;
}

bool valid_compiled_enthalpy_workspace(
    const PressureEnergyEnthalpyCompiledWorkspace& workspace,
    Int3 cells) noexcept {
  return valid_scalar_view(workspace.local_diagonal, cells) &&
         valid_scalar_view(workspace.response_stage, cells) &&
         workspace.directional_branches.values.data != nullptr &&
         workspace.directional_branches.values.size ==
             pressure_energy_branch_count(cells) &&
         valid_face_view(as_const(workspace.thermal_conductance.x),
                         CartesianAxis::x, cells) &&
         valid_face_view(as_const(workspace.thermal_conductance.y),
                         CartesianAxis::y, cells) &&
         valid_face_view(as_const(workspace.thermal_conductance.z),
                         CartesianAxis::z, cells);
}

PlanFingerprint compiled_enthalpy_local_binding(
    const PressureEnergyEnthalpyBinding& binding,
    const FrozenConvectionBranchPlan& branches) noexcept {
  std::uint64_t hash = hash_mix(kFnvOffset, kPressureEnergyEnthalpySchema);
  hash = hash_mix(hash, branches.revision);
  hash = hash_mix(hash, branches.branch_authority);
  hash = hash_mix(hash, branches.local_binding);
  hash = mix_local_field(
      hash, as_const(binding.workspace.compiled.local_diagonal));
  hash = mix_local_field(hash,
                         as_const(binding.workspace.compiled.response_stage));
  hash = hash_mix(hash,
                  reinterpret_cast<std::uintptr_t>(branches.values.data));
  hash = hash_mix(hash, branches.values.size);
  const ConstFaceFieldView faces[]{
      as_const(binding.workspace.compiled.thermal_conductance.x),
      as_const(binding.workspace.compiled.thermal_conductance.y),
      as_const(binding.workspace.compiled.thermal_conductance.z)};
  for (ConstFaceFieldView face : faces) hash = mix_local_face(hash, face);
  return nonzero_hash(hash);
}

RevisionToken compiled_enthalpy_numeric_revision(
    const PressureEnergyEnthalpyBinding& binding,
    const FrozenConvectionBranchPlan& branches) noexcept {
  std::uint64_t hash = compiled_enthalpy_local_binding(binding, branches);
  hash = hash_mix(hash, binding.assembled_diagonal.revision);
  hash = hash_mix(hash, binding.target_enthalpy.revision);
  hash = hash_mix(hash, binding.density_enthalpy_derivative.revision);
  hash = hash_mix(hash, binding.heat_capacity.revision);
  hash = hash_mix(hash, binding.thermal_conductivity.revision);
  hash = hash_mix(hash, binding.enthalpy_diffusivity.revision);
  hash = hash_mix(hash, binding.authority.target_time);
  return nonzero_hash(hash);
}

template <class Left, class Right>
bool overlaps(BasicFieldView<Left> left,
              BasicFieldView<Right> right) noexcept {
  return detail::field_views_overlap(left, right);
}

bool valid_component_certificate(LinearOperatorCertificate certificate,
                                 Int3 shape) noexcept {
  return valid_identity(certificate.identity) &&
         certificate.collective_fingerprint != 0U &&
         same_shape(certificate.local_shape, shape) &&
         (certificate.operator_class == LinearOperatorClass::spd ||
          certificate.operator_class == LinearOperatorClass::nonsymmetric);
}

bool components_current(
    const LinearOperator* continuity_pressure,
    LinearOperatorCertificate continuity_pressure_certificate,
    const LinearOperator* energy_pressure,
    LinearOperatorCertificate energy_pressure_certificate,
    const LinearOperator* energy_enthalpy,
    LinearOperatorCertificate energy_enthalpy_certificate) noexcept {
  return continuity_pressure != nullptr && energy_pressure != nullptr &&
         energy_enthalpy != nullptr &&
         same_certificate(continuity_pressure->certificate(),
                          continuity_pressure_certificate) &&
         same_certificate(energy_pressure->certificate(),
                          energy_pressure_certificate) &&
         same_certificate(energy_enthalpy->certificate(),
                          energy_enthalpy_certificate);
}

Status component_apply(const LinearOperator& component, FieldView input,
                       FieldView output,
                       LinearOperatorFailureProvenance& failure) noexcept {
  const Status status = component.apply(input, output);
  if (status) return status;
  const LinearOperatorFailureProvenance component_failure =
      component.failure_provenance();
  if (component_failure.status.code == status.code &&
      component_failure.status.detail == status.detail) {
    failure = component_failure;
  } else {
    failure = {status, LinearOperatorStatusScope::rank_local, -1};
  }
  return status;
}

Status capture_component_failure(
    const LinearOperator& component, Status status,
    LinearOperatorFailureProvenance& failure) noexcept {
  if (status) return status;
  const LinearOperatorFailureProvenance component_failure =
      component.failure_provenance();
  if (component_failure.status.code == status.code &&
      component_failure.status.detail == status.detail) {
    failure = component_failure;
  } else {
    failure = {status, LinearOperatorStatusScope::rank_local, -1};
  }
  return status;
}

template <class Function>
void for_each_cell(Int3 cells, Function&& function) noexcept {
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        function(Int3{x, y, z});
}

double pressure_gradient_component(const CartesianKernelPlan& kernels,
                                   ConstFieldView pressure, Int3 cell,
                                   std::size_t axis) noexcept {
  const std::int32_t normal =
      axis == 0U ? cell.x : (axis == 1U ? cell.y : cell.z);
  const detail::DerivativeWeights weights =
      kernels.geometry_kind() == GeometryKind::uniform
          ? detail::metric_derivative_weights<true>(kernels, axis, normal)
          : detail::metric_derivative_weights<false>(kernels, axis, normal);
  Int3 minus = cell;
  Int3 plus = cell;
  if (axis == 0U) {
    --minus.x;
    ++plus.x;
  } else if (axis == 1U) {
    --minus.y;
    ++plus.y;
  } else {
    --minus.z;
    ++plus.z;
  }
  return weights.minus * pressure.unchecked(minus, 0U) +
         weights.centre * pressure.unchecked(cell, 0U) +
         weights.plus * pressure.unchecked(plus, 0U);
}

bool runtime_views_valid(ConstFieldView first, ConstFieldView second,
                         FieldView output, Int3 shape,
                         const PressureEnergySchurWorkspace& workspace,
                         ConstFieldView diagonal,
                         ConstFieldView row_scale) noexcept {
  if (!valid_scalar_view(first, shape) ||
      !valid_scalar_view(second, shape) ||
      !valid_scalar_view(output, shape) || overlaps(first, second) ||
      overlaps(first, output) || overlaps(second, output)) {
    return false;
  }
  const ConstFieldView workspaces[]{as_const(workspace.continuity_response),
                                    as_const(workspace.eliminated_enthalpy),
                                    as_const(workspace.energy_response)};
  for (ConstFieldView value : workspaces) {
    if (overlaps(first, value) || overlaps(second, value) ||
        overlaps(as_const(output), value) || overlaps(diagonal, value) ||
        overlaps(row_scale, value)) {
      return false;
    }
  }
  return !overlaps(first, diagonal) && !overlaps(second, diagonal) &&
         !overlaps(as_const(output), diagonal) &&
         !overlaps(first, row_scale) && !overlaps(second, row_scale) &&
         !overlaps(as_const(output), row_scale);
}

}  // namespace

bool PressureEnergyGlobalizationSelectionCertificate::valid() const noexcept {
  const double expected_alpha =
      std::ldexp(1.0, -static_cast<int>(selected_halvings));
  const bool valid_alpha =
      extrapolated
          ? selected_halvings == 0U && std::isfinite(alpha) && alpha > 1.0 &&
                alpha <= kPressureEnergyAitkenMaximumAlpha
          : std::isfinite(alpha) && alpha == expected_alpha;
  const double expected_baseline_merit = std::hypot(
      baseline_normalized_continuity, baseline_normalized_energy);
  const double expected_candidate_merit = std::hypot(
      candidate_normalized_continuity, candidate_normalized_energy);
  const double expected_armijo =
      (1.0 - kPressureEnergyGlobalizationArmijoCoefficient * alpha) *
      baseline_merit;
  if (scope != PressureEnergyGlobalizationScope::
                   frozen_momentum_continuity_energy_globalization ||
      (corrector != 1U && corrector != 2U) || target_time == 0U ||
      correction_direction == 0U || baseline_state_provenance == 0U ||
      baseline_mass_flux_provenance == 0U ||
      candidate_state_provenance == 0U ||
      candidate_mass_flux_provenance == 0U ||
      candidate_state_provenance == baseline_state_provenance ||
      candidate_mass_flux_provenance == baseline_mass_flux_provenance ||
      selected_halvings >= kPressureEnergyGlobalizationCandidateCount ||
      !valid_alpha ||
      !std::isfinite(baseline_normalized_continuity) ||
      baseline_normalized_continuity < 0.0 ||
      !std::isfinite(baseline_normalized_energy) ||
      baseline_normalized_energy < 0.0 ||
      !std::isfinite(candidate_normalized_continuity) ||
      candidate_normalized_continuity < 0.0 ||
      !std::isfinite(candidate_normalized_energy) ||
      candidate_normalized_energy < 0.0 ||
      !std::isfinite(baseline_merit) ||
      baseline_merit != expected_baseline_merit ||
      !std::isfinite(candidate_merit) ||
      candidate_merit != expected_candidate_merit ||
      !std::isfinite(armijo_upper_bound) ||
      armijo_upper_bound != expected_armijo ||
      !thermodynamically_admissible || !state_and_flux_finite ||
      !strict_merit_decrease || !(candidate_merit < baseline_merit) ||
      !armijo_sufficient_decrease ||
      !(candidate_merit <= armijo_upper_bound) || full_nonlinear_newton) {
    return false;
  }

  PressureEnergyGlobalizationSample baseline;
  baseline.alpha = 0.0;
  baseline.global_normalized_continuity = baseline_normalized_continuity;
  baseline.global_normalized_energy = baseline_normalized_energy;
  baseline.thermodynamically_admissible = true;
  baseline.state_and_flux_finite = true;
  baseline.corrector = corrector;
  baseline.target_time = target_time;
  baseline.correction_direction = correction_direction;
  baseline.state_provenance = baseline_state_provenance;
  baseline.mass_flux_provenance = baseline_mass_flux_provenance;
  PressureEnergyGlobalizationSample selected = baseline;
  selected.alpha = alpha;
  selected.global_normalized_continuity = candidate_normalized_continuity;
  selected.global_normalized_energy = candidate_normalized_energy;
  selected.state_provenance = candidate_state_provenance;
  selected.mass_flux_provenance = candidate_mass_flux_provenance;
  const PlanFingerprint expected_provenance =
      extrapolated
          ? pressure_energy_extrapolation_selection_provenance(baseline,
                                                                selected)
          : pressure_energy_globalization_selection_provenance(
                baseline, selected, selected_halvings);
  return selection_provenance == expected_provenance;
}

PlanFingerprint pressure_energy_frozen_face_enthalpy_local_binding(
    const PressureEnergyFrozenFaceEnthalpy& frozen) noexcept {
  std::uint64_t hash = hash_mix(kFnvOffset, kPressureEnergyPressureFluxSchema);
  hash = hash_mix(hash, frozen.revision);
  hash = hash_mix(hash, frozen.reconstruction);
  hash = mix_local_face(hash, frozen.x);
  hash = mix_local_face(hash, frozen.y);
  hash = mix_local_face(hash, frozen.z);
  return nonzero_hash(hash);
}

bool PressureEnergyPressureFluxOperator::exact_pressure_work_current()
    const noexcept {
  return certificate_.valid() &&
         certificate_.pressure_work_scope ==
             PressureEnergyPressureWorkScope::exact_cartesian &&
         pressure_work_.valid() && pressure_work_.issuer_ != nullptr &&
         pressure_work_.issuer_->matches_cartesian_pressure_work_linearization(
             pressure_work_, false);
}

bool PressureEnergySchurBlockAuthority::valid() const noexcept {
  const bool components =
      valid_component_certificate(energy_pressure_,
                                  energy_pressure_.local_shape) &&
      valid_component_certificate(energy_enthalpy_,
                                  energy_pressure_.local_shape) &&
      same_identity(energy_pressure_.identity,
                    energy_enthalpy_.identity) &&
      same_shape(energy_pressure_.local_shape,
                 energy_enthalpy_.local_shape) &&
      collective_fingerprint_ != 0U && rank_local_revision_ != 0U &&
      ((activity_local_fingerprint_ == 0U) ==
       (activity_collective_fingerprint_ == 0U));
  if (!components) return false;
  return scope_ ==
             PressureEnergySchurBlockScope::exact_cartesian_frozen_spatial
             ? activity_local_fingerprint_ == 0U &&
                   activity_collective_fingerprint_ == 0U
         : scope_ == PressureEnergySchurBlockScope::
                         ibm_cartesian_spatial_quasi_newton
             ? activity_local_fingerprint_ != 0U &&
                   activity_collective_fingerprint_ != 0U
         : scope_ == PressureEnergySchurBlockScope::
                         ibm_double_diagonal_quasi_newton
             ? activity_local_fingerprint_ != 0U &&
                   activity_collective_fingerprint_ != 0U
         : scope_ == PressureEnergySchurBlockScope::
                         generic_algebraic_quasi_newton;
}

bool PressureEnergySchurBlockAuthority::matches_operators(
    const LinearOperator* energy_pressure,
    const LinearOperator* energy_enthalpy) const noexcept {
  if (!valid() || energy_pressure == nullptr || energy_enthalpy == nullptr ||
      !same_certificate(energy_pressure_, energy_pressure->certificate()) ||
      !same_certificate(energy_enthalpy_, energy_enthalpy->certificate())) {
    return false;
  }
  if (scope_ ==
      PressureEnergySchurBlockScope::exact_cartesian_frozen_spatial) {
    return energy_pressure == energy_pressure_operator_ &&
           energy_enthalpy == energy_enthalpy_operator_ &&
           static_cast<const PressureEnergyPressureFluxOperator*>(
               energy_pressure_operator_)
               ->exact_pressure_work_current();
  }
  if (scope_ == PressureEnergySchurBlockScope::
                    ibm_cartesian_spatial_quasi_newton ||
      scope_ == PressureEnergySchurBlockScope::
                    ibm_double_diagonal_quasi_newton) {
    return energy_pressure == energy_pressure_operator_ &&
           energy_enthalpy == energy_enthalpy_operator_;
  }
  return energy_pressure_operator_ == nullptr &&
         energy_enthalpy_operator_ == nullptr;
}

Status PressureEnergySchurBlockAuthority::exact_cartesian(
    const PressureEnergyPressureFluxOperator& energy_pressure_operator,
    const PressureEnergyEnthalpyOperator& energy_enthalpy_operator,
    PressureEnergySchurBlockAuthority& out) noexcept {
  const PressureEnergyPressureFluxCertificate& energy_pressure =
      energy_pressure_operator.pressure_flux_certificate();
  const PressureEnergyEnthalpyCertificate& energy_enthalpy =
      energy_enthalpy_operator.enthalpy_certificate();
  const bool valid =
      energy_pressure_operator.exact_pressure_work_current() &&
      energy_pressure.valid() && energy_enthalpy.valid() &&
      energy_pressure.pressure_work_scope ==
          PressureEnergyPressureWorkScope::exact_cartesian &&
      energy_pressure.full_cartesian_pressure_work &&
      !energy_pressure.flux_only_quasi_newton &&
      energy_pressure.activity_local_fingerprint == 0U &&
      energy_pressure.activity_collective_fingerprint == 0U &&
      energy_pressure.inactive_cells == 0U &&
      energy_enthalpy.activity_local_fingerprint == 0U &&
      energy_enthalpy.activity_collective_fingerprint == 0U &&
      energy_enthalpy.inactive_cells == 0U &&
      energy_enthalpy.exact_cartesian_spatial_response &&
      energy_enthalpy.exact_temperature_space_conduction &&
      !energy_enthalpy.ibm_spatial_derivative &&
      same_identity(energy_pressure.linear.identity,
                    energy_enthalpy.linear.identity) &&
      same_shape(energy_pressure.linear.local_shape,
                 energy_enthalpy.linear.local_shape);
  if (!valid) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }
  PressureEnergySchurBlockAuthority candidate;
  candidate.scope_ =
      PressureEnergySchurBlockScope::exact_cartesian_frozen_spatial;
  candidate.energy_pressure_operator_ = &energy_pressure_operator;
  candidate.energy_enthalpy_operator_ = &energy_enthalpy_operator;
  candidate.energy_pressure_ = energy_pressure.linear;
  candidate.energy_enthalpy_ = energy_enthalpy.linear;
  std::uint64_t collective =
      hash_mix(kFnvOffset, UINT64_C(0x657863617274626c));
  collective = hash_mix(collective,
                        energy_pressure.linear.collective_fingerprint);
  collective = hash_mix(collective,
                        energy_enthalpy.linear.collective_fingerprint);
  collective = hash_mix(collective,
                        energy_pressure.pressure_work_linearization != 0U);
  candidate.collective_fingerprint_ = nonzero_hash(collective);
  std::uint64_t local = hash_mix(candidate.collective_fingerprint_,
                                 energy_pressure.binding_revision);
  local = hash_mix(local, energy_enthalpy.binding_revision);
  local = hash_mix(local, energy_pressure.pressure_work_linearization);
  candidate.rank_local_revision_ = nonzero_hash(local);
  if (!candidate.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }
  out = candidate;
  return {};
}

Status PressureEnergySchurBlockAuthority::
    ibm_cartesian_spatial_quasi_newton(
        const PressureEnergyPressureFluxOperator& energy_pressure_operator,
        const PressureEnergyEnthalpyOperator& energy_enthalpy_operator,
        PressureEnergySchurBlockAuthority& out) noexcept {
  const PressureEnergyPressureFluxCertificate& energy_pressure =
      energy_pressure_operator.pressure_flux_certificate();
  const PressureEnergyEnthalpyCertificate& energy_enthalpy =
      energy_enthalpy_operator.enthalpy_certificate();
  const bool valid =
      energy_pressure.valid() && energy_enthalpy.valid() &&
      energy_pressure.pressure_work_scope ==
          PressureEnergyPressureWorkScope::flux_only_quasi_newton &&
      energy_pressure.pressure_work_linearization == 0U &&
      !energy_pressure.full_cartesian_pressure_work &&
      energy_pressure.flux_only_quasi_newton &&
      energy_pressure.activity_local_fingerprint != 0U &&
      energy_pressure.activity_collective_fingerprint != 0U &&
      energy_pressure.activity_local_fingerprint ==
          energy_enthalpy.activity_local_fingerprint &&
      energy_pressure.activity_collective_fingerprint ==
          energy_enthalpy.activity_collective_fingerprint &&
      energy_pressure.active_cells == energy_enthalpy.active_cells &&
      energy_pressure.inactive_cells == energy_enthalpy.inactive_cells &&
      energy_enthalpy.exact_cartesian_spatial_response &&
      energy_enthalpy.exact_temperature_space_conduction &&
      !energy_enthalpy.ibm_spatial_derivative &&
      energy_enthalpy.inactive_rows_identity &&
      energy_enthalpy.inactive_interfaces_zero &&
      same_identity(energy_pressure.linear.identity,
                    energy_enthalpy.linear.identity) &&
      same_shape(energy_pressure.linear.local_shape,
                 energy_enthalpy.linear.local_shape);
  if (!valid) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }
  PressureEnergySchurBlockAuthority candidate;
  candidate.scope_ =
      PressureEnergySchurBlockScope::ibm_cartesian_spatial_quasi_newton;
  candidate.energy_pressure_operator_ = &energy_pressure_operator;
  candidate.energy_enthalpy_operator_ = &energy_enthalpy_operator;
  candidate.energy_pressure_ = energy_pressure.linear;
  candidate.energy_enthalpy_ = energy_enthalpy.linear;
  candidate.activity_local_fingerprint_ =
      energy_pressure.activity_local_fingerprint;
  candidate.activity_collective_fingerprint_ =
      energy_pressure.activity_collective_fingerprint;
  std::uint64_t collective =
      hash_mix(kFnvOffset, UINT64_C(0x69626d737061746c));
  collective = hash_mix(collective,
                        energy_pressure.linear.collective_fingerprint);
  collective = hash_mix(collective,
                        energy_enthalpy.linear.collective_fingerprint);
  collective = hash_mix(collective,
                        candidate.activity_collective_fingerprint_);
  collective = hash_mix(collective,
                        energy_pressure.flux_only_quasi_newton);
  collective = hash_mix(collective,
                        energy_enthalpy.exact_cartesian_spatial_response);
  collective = hash_mix(collective,
                        energy_enthalpy.ibm_spatial_derivative);
  candidate.collective_fingerprint_ = nonzero_hash(collective);
  std::uint64_t local = hash_mix(candidate.collective_fingerprint_,
                                 energy_pressure.binding_revision);
  local = hash_mix(local, energy_enthalpy.binding_revision);
  local = hash_mix(local, candidate.activity_local_fingerprint_);
  candidate.rank_local_revision_ = nonzero_hash(local);
  if (!candidate.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }
  out = candidate;
  return {};
}

Status PressureEnergySchurBlockAuthority::ibm_double_diagonal(
    const PressureEnergyDiagonalOperator& energy_pressure_operator,
    const PressureEnergyDiagonalOperator& energy_enthalpy_operator,
    PressureEnergySchurBlockAuthority& out) noexcept {
  const PressureEnergyDiagonalCertificate& energy_pressure =
      energy_pressure_operator.diagonal_certificate();
  const PressureEnergyDiagonalCertificate& energy_enthalpy =
      energy_enthalpy_operator.diagonal_certificate();
  const bool valid =
      energy_pressure.valid() && energy_enthalpy.valid() &&
      energy_pressure.inactive_diagonal == 0.0 &&
      energy_enthalpy.inactive_diagonal == 1.0 &&
      energy_pressure.activity_local_fingerprint != 0U &&
      energy_pressure.activity_collective_fingerprint != 0U &&
      energy_pressure.activity_local_fingerprint ==
          energy_enthalpy.activity_local_fingerprint &&
      energy_pressure.activity_collective_fingerprint ==
          energy_enthalpy.activity_collective_fingerprint &&
      energy_pressure.active_cells == energy_enthalpy.active_cells &&
      energy_pressure.inactive_cells == energy_enthalpy.inactive_cells &&
      same_identity(energy_pressure.linear.identity,
                    energy_enthalpy.linear.identity) &&
      same_shape(energy_pressure.linear.local_shape,
                 energy_enthalpy.linear.local_shape);
  if (!valid) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }
  PressureEnergySchurBlockAuthority candidate;
  candidate.scope_ =
      PressureEnergySchurBlockScope::ibm_double_diagonal_quasi_newton;
  candidate.energy_pressure_operator_ = &energy_pressure_operator;
  candidate.energy_enthalpy_operator_ = &energy_enthalpy_operator;
  candidate.energy_pressure_ = energy_pressure.linear;
  candidate.energy_enthalpy_ = energy_enthalpy.linear;
  candidate.activity_local_fingerprint_ =
      energy_pressure.activity_local_fingerprint;
  candidate.activity_collective_fingerprint_ =
      energy_pressure.activity_collective_fingerprint;
  std::uint64_t collective =
      hash_mix(kFnvOffset, UINT64_C(0x69626d6464696167));
  collective = hash_mix(collective,
                        energy_pressure.linear.collective_fingerprint);
  collective = hash_mix(collective,
                        energy_enthalpy.linear.collective_fingerprint);
  collective = hash_mix(collective,
                        candidate.activity_collective_fingerprint_);
  collective = hash_mix(collective,
                        double_bits(energy_pressure.inactive_diagonal));
  collective = hash_mix(collective,
                        double_bits(energy_enthalpy.inactive_diagonal));
  candidate.collective_fingerprint_ = nonzero_hash(collective);
  std::uint64_t local = hash_mix(candidate.collective_fingerprint_,
                                 energy_pressure.diagonal_revision);
  local = hash_mix(local, energy_enthalpy.diagonal_revision);
  local = hash_mix(local, candidate.activity_local_fingerprint_);
  candidate.rank_local_revision_ = nonzero_hash(local);
  if (!candidate.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }
  out = candidate;
  return {};
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
Status PressureEnergySchurBlockAuthority::
    generic_algebraic_quasi_newton_for_test(
        LinearOperatorCertificate energy_pressure,
        LinearOperatorCertificate energy_enthalpy,
        PressureEnergyCellActivity activity,
        PressureEnergySchurBlockAuthority& out) noexcept {
  if (!valid_component_certificate(energy_pressure,
                                   energy_pressure.local_shape) ||
      !valid_component_certificate(energy_enthalpy,
                                   energy_pressure.local_shape) ||
      !same_identity(energy_pressure.identity, energy_enthalpy.identity) ||
      !same_shape(energy_pressure.local_shape,
                  energy_enthalpy.local_shape) ||
      !valid_activity(activity, energy_pressure.local_shape)) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }
  PressureEnergySchurBlockAuthority candidate;
  candidate.scope_ =
      PressureEnergySchurBlockScope::generic_algebraic_quasi_newton;
  candidate.energy_pressure_ = energy_pressure;
  candidate.energy_enthalpy_ = energy_enthalpy;
  candidate.activity_local_fingerprint_ = activity.local_fingerprint;
  candidate.activity_collective_fingerprint_ =
      activity.collective_fingerprint;
  std::uint64_t collective =
      hash_mix(kFnvOffset, UINT64_C(0x7465737467656e71));
  collective = hash_mix(collective,
                        energy_pressure.collective_fingerprint);
  collective = hash_mix(collective,
                        energy_enthalpy.collective_fingerprint);
  collective = hash_mix(collective,
                        activity.collective_fingerprint);
  candidate.collective_fingerprint_ = nonzero_hash(collective);
  std::uint64_t local = hash_mix(candidate.collective_fingerprint_,
                                 energy_pressure.identity.numeric);
  local = hash_mix(local, energy_enthalpy.identity.numeric);
  local = hash_mix(local, activity.local_fingerprint);
  candidate.rank_local_revision_ = nonzero_hash(local);
  if (!candidate.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }
  out = candidate;
  return {};
}
#endif

double pressure_correction_jump(ConstFieldView correction,
                                const CartesianGeometryPlan& geometry,
                                MeshPatch patch, const BoundaryPlan& boundary,
                                CartesianAxis axis, Int3 face) noexcept {
  PressureCorrectionBoundaryPlan plan;
  if (!PressureCorrectionBoundaryPlan::compile(geometry, patch, boundary,
                                               plan)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return plan.jump(correction, axis, face);
}

double pressure_correction_mass_flux_response(
    ConstFieldView correction, const CartesianGeometryPlan& geometry,
    MeshPatch patch, const BoundaryPlan& boundary, CartesianAxis axis,
    Int3 face, double pressure_face_coefficient,
    double correction_scale) noexcept {
  PressureCorrectionBoundaryPlan plan;
  if (!PressureCorrectionBoundaryPlan::compile(geometry, patch, boundary,
                                               plan)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return plan.mass_flux_response(correction, axis, face,
                                 pressure_face_coefficient,
                                 correction_scale);
}

Status PressureEnergyPressureFluxOperator::bind(
    const PressureEnergyPressureFluxBinding& binding,
    PressureEnergyPressureFluxOperator& out,
    PressureEnergyPressureFluxCertificate& certificate) noexcept {
  certificate = {};
  const Int3 cells = binding.patch.cells;
  const PisoCartesianPressureWorkLinearization& pressure_work =
      binding.pressure_work;
  const bool empty_pressure_work =
      pressure_work.issuer_ == nullptr && pressure_work.kernels_ == nullptr &&
      pressure_work.geometry_ == nullptr &&
      pressure_work.patch_.cells.x == 0 &&
      pressure_work.patch_.cells.y == 0 &&
      pressure_work.patch_.cells.z == 0 &&
      !pressure_work.intermediate_.valid() &&
      !pressure_work.pressure_.valid() &&
      pressure_work.target_pressure_perturbation_.base == nullptr &&
      pressure_work.target_velocity_.base == nullptr &&
      pressure_work.h_by_a_.base == nullptr &&
      pressure_work.r_au_.base == nullptr &&
      pressure_work.authority_ == 0U &&
      pressure_work.numeric_fingerprint_ == 0U;
  const bool exact_pressure_work = pressure_work.valid();
  PressureCorrectionBoundaryPlan pressure_boundary;
  const Status pressure_boundary_status =
      binding.geometry != nullptr && binding.boundary != nullptr
          ? PressureCorrectionBoundaryPlan::compile(
                *binding.geometry, binding.patch, *binding.boundary,
                pressure_boundary)
          : Status{StatusCode::invalid_plan,
                   kPressureEnergyPressureFluxBinding};
  const bool plans_match =
      static_cast<bool>(pressure_boundary_status) &&
      binding.geometry != nullptr && binding.boundary != nullptr &&
      binding.geometry->fingerprint() != 0U &&
      binding.boundary->semantic_fingerprint() != 0U &&
      same_shape(binding.boundary->local_cells(), cells) &&
      binding.pressure.valid() && binding.intermediate.valid() &&
      binding.pressure.plan == binding.intermediate.plan &&
      binding.pressure.intermediate == binding.intermediate.dependency &&
      binding.pressure.corrector == binding.intermediate.corrector &&
      binding.pressure.geometry == binding.geometry->topology_revision() &&
      binding.pressure.numeric_boundary ==
          pressure_boundary.certificate().source_revision &&
      binding.pressure.thermophysical_boundary_semantics ==
          binding.intermediate.thermophysical_boundary_semantics &&
      binding.pressure.thermophysical_boundary_target ==
          binding.intermediate.thermophysical_boundary_target &&
      binding.pressure.thermophysical_boundary_rank_local_binding ==
          binding.intermediate
              .thermophysical_boundary_rank_local_binding &&
      binding.pressure.thermophysical_boundary_collective_lineage ==
          binding.intermediate
              .thermophysical_boundary_collective_lineage &&
      binding.pressure.thermophysical_boundary_rank_local_lineage ==
          binding.intermediate
              .thermophysical_boundary_rank_local_lineage &&
      (empty_pressure_work ||
       (exact_pressure_work && pressure_work.geometry_ == binding.geometry &&
        same_shape(pressure_work.patch_.cells, binding.patch.cells) &&
        pressure_work.patch_.begin.x == binding.patch.begin.x &&
        pressure_work.patch_.begin.y == binding.patch.begin.y &&
        pressure_work.patch_.begin.z == binding.patch.begin.z &&
        pressure_work.intermediate_.dependency ==
            binding.intermediate.dependency &&
        pressure_work.pressure_.state == binding.pressure.state &&
        pressure_work.issuer_ != nullptr &&
        pressure_work.issuer_->matches_cartesian_pressure_work_linearization(
            pressure_work, true)));
  const std::array<HaloFieldSpec, 1U> halo_fields{
      {{binding.services.solution_field, binding.services.halo_width, 1U}}};
  const bool services_valid =
      binding.services.communicator != MPI_COMM_NULL &&
      binding.services.halo != nullptr && binding.services.halo->ready() &&
      binding.services.halo_stage != 0U && binding.services.halo_width >= 1U &&
      static_cast<bool>(binding.services.halo->validate_contract(
          binding.services.communicator, binding.patch,
          {halo_fields.data(), halo_fields.size()},
          binding.boundary == nullptr ? HaloTopology{}
                                      : binding.boundary->halo_topology()));
  const PressureEnergyFrozenFaceEnthalpy frozen = binding.frozen_face_enthalpy;
  const bool views_valid =
      valid_identity(binding.identity) &&
      valid_scalar_view(binding.temporal_diagonal, cells) &&
      valid_face_view(binding.x_pressure_coefficient, CartesianAxis::x,
                      cells) &&
      valid_face_view(binding.y_pressure_coefficient, CartesianAxis::y,
                      cells) &&
      valid_face_view(binding.z_pressure_coefficient, CartesianAxis::z,
                      cells) &&
      detail::valid_flux_view(binding.target_flux, cells,
                              binding.target_flux.revision) &&
      !binding.target_flux.certificate.valid() &&
      valid_face_view(frozen.x, CartesianAxis::x, cells) &&
      valid_face_view(frozen.y, CartesianAxis::y, cells) &&
      valid_face_view(frozen.z, CartesianAxis::z, cells) &&
      frozen.revision != 0U && frozen.reconstruction != 0U &&
      frozen.local_binding != 0U &&
      frozen.local_binding ==
          pressure_energy_frozen_face_enthalpy_local_binding(frozen) &&
      valid_pressure_flux_activity(binding.activity, cells) &&
      (!exact_pressure_work ||
       (binding.activity.cells.size == 0U &&
        binding.activity.x_faces.size == 0U &&
        binding.activity.y_faces.size == 0U &&
        binding.activity.z_faces.size == 0U &&
        binding.activity.local_fingerprint == 0U &&
        binding.activity.collective_fingerprint == 0U));
  if (!plans_match || !services_valid || !views_valid) {
    return {StatusCode::invalid_plan, kPressureEnergyPressureFluxBinding};
  }

  const ConstFaceFieldView faces[]{binding.x_pressure_coefficient,
                                   binding.y_pressure_coefficient,
                                   binding.z_pressure_coefficient,
                                   binding.target_flux.x,
                                   binding.target_flux.y,
                                   binding.target_flux.z,
                                   frozen.x,
                                   frozen.y,
                                   frozen.z};
  bool aliases = false;
  for (ConstFaceFieldView face : faces) {
    aliases = aliases ||
              detail::cell_face_views_overlap(binding.temporal_diagonal, face);
  }
  for (std::size_t left = 0U; left < std::size(faces); ++left) {
    for (std::size_t right = left + 1U; right < std::size(faces); ++right) {
      aliases =
          aliases || detail::face_views_overlap(faces[left], faces[right]);
    }
  }
  if (exact_pressure_work) {
    const ConstFieldView pressure_work_fields[]{
        pressure_work.target_pressure_perturbation_,
        pressure_work.target_velocity_, pressure_work.h_by_a_,
        pressure_work.r_au_};
    for (ConstFieldView field : pressure_work_fields) {
      aliases = aliases || overlaps(binding.temporal_diagonal, field);
      for (ConstFaceFieldView face : faces) {
        aliases = aliases || detail::cell_face_views_overlap(field, face);
      }
    }
    for (std::size_t left = 0U; left < std::size(pressure_work_fields);
         ++left) {
      for (std::size_t right = left + 1U;
           right < std::size(pressure_work_fields); ++right) {
        aliases = aliases ||
                  overlaps(pressure_work_fields[left],
                           pressure_work_fields[right]);
      }
    }
  }
  if (aliases) {
    return {StatusCode::invalid_plan, kPressureEnergyPressureFluxBinding};
  }

  std::size_t active_cells = 0U;
  bool admissible = true;
  for_each_cell(cells, [&](Int3 cell) {
    const bool active =
        binding.activity.cells.size == 0U ||
        binding.activity.cells.data[cell_offset(cells, cell)] == 1U;
    if (!active) return;
    ++active_cells;
    if (!std::isfinite(binding.temporal_diagonal.unchecked(cell, 0U))) {
      admissible = false;
    }
  });
  const CartesianAxis axes[]{CartesianAxis::x, CartesianAxis::y,
                             CartesianAxis::z};
  for (CartesianAxis axis : axes) {
    const ConstFaceFieldView coefficient = select_face(
        binding.x_pressure_coefficient, binding.y_pressure_coefficient,
        binding.z_pressure_coefficient, axis);
    const ConstFaceFieldView target =
        select_face(binding.target_flux.x, binding.target_flux.y,
                    binding.target_flux.z, axis);
    const ConstFaceFieldView enthalpy =
        select_face(frozen.x, frozen.y, frozen.z, axis);
    for (std::int32_t z = 0; z < coefficient.extents.z; ++z) {
      for (std::int32_t y = 0; y < coefficient.extents.y; ++y) {
        for (std::int32_t x = 0; x < coefficient.extents.x; ++x) {
          const Int3 face{x, y, z};
          if (!active_face(binding.activity, axis, cells, face)) continue;
          const double a = coefficient.unchecked(face);
          const double phi = target.unchecked(face);
          const double h = enthalpy.unchecked(face);
          if (!std::isfinite(a) || a < 0.0 || !std::isfinite(phi) ||
              !std::isfinite(h)) {
            admissible = false;
          }
        }
      }
    }
  }
  if (!admissible) {
    return {StatusCode::rejected_step, kPressureEnergyPressureFluxBinding};
  }

  PressureEnergyPressureFluxOperator candidate;
  candidate.pressure_boundary_ = std::move(pressure_boundary);
  candidate.patch_ = binding.patch;
  candidate.services_ = binding.services;
  candidate.temporal_diagonal_ = binding.temporal_diagonal;
  candidate.x_pressure_coefficient_ = binding.x_pressure_coefficient;
  candidate.y_pressure_coefficient_ = binding.y_pressure_coefficient;
  candidate.z_pressure_coefficient_ = binding.z_pressure_coefficient;
  candidate.target_flux_ = binding.target_flux;
  candidate.frozen_face_enthalpy_ = frozen;
  candidate.activity_ = binding.activity;
  candidate.pressure_work_ = pressure_work;
  candidate.certificate_.linear = {
      binding.identity,
      pressure_flux_collective_fingerprint(
          binding, candidate.pressure_boundary_.certificate()),
      cells,
      LinearOperatorClass::nonsymmetric};
  candidate.certificate_.binding_revision = pressure_flux_binding_revision(
      binding, candidate.pressure_boundary_.certificate());
  if (exact_pressure_work) {
    candidate.certificate_.binding_revision = nonzero_hash(hash_mix(
        candidate.certificate_.binding_revision, pressure_work.authority_));
  }
  candidate.certificate_.pressure_system = binding.pressure.state;
  candidate.certificate_.intermediate_dependency =
      binding.intermediate.dependency;
  candidate.certificate_.pressure_face_coefficients =
      binding.intermediate.pressure_face_coefficient;
  candidate.certificate_.target_flux_identity = binding.intermediate.phi_h_by_a;
  candidate.certificate_.target_flux_revision = binding.target_flux.revision;
  candidate.certificate_.temporal_diagonal = binding.temporal_diagonal.revision;
  candidate.certificate_.frozen_face_enthalpy = frozen.revision;
  candidate.certificate_.frozen_face_enthalpy_local_binding =
      frozen.local_binding;
  candidate.certificate_.geometry = binding.pressure.geometry;
  candidate.certificate_.numeric_boundary = binding.pressure.numeric_boundary;
  candidate.certificate_.pressure_boundary_semantic =
      candidate.pressure_boundary_.certificate().semantic;
  candidate.certificate_.pressure_boundary_rank_local_layout =
      candidate.pressure_boundary_.certificate().rank_local_layout;
  candidate.certificate_.frozen_reconstruction = frozen.reconstruction;
  candidate.certificate_.activity_local_fingerprint =
      binding.activity.local_fingerprint;
  candidate.certificate_.activity_collective_fingerprint =
      binding.activity.collective_fingerprint;
  candidate.certificate_.active_cells = active_cells;
  candidate.certificate_.inactive_cells = cell_count(cells) - active_cells;
  candidate.certificate_.pressure_work_scope =
      exact_pressure_work
          ? PressureEnergyPressureWorkScope::exact_cartesian
          : PressureEnergyPressureWorkScope::flux_only_quasi_newton;
  candidate.certificate_.pressure_work_linearization =
      exact_pressure_work ? pressure_work.authority_ : 0U;
  candidate.certificate_.conservative_face_response = true;
  candidate.certificate_.exact_piso_flux_jump = true;
  candidate.certificate_.frozen_target_enthalpy = true;
  candidate.certificate_.full_cartesian_pressure_work = exact_pressure_work;
  candidate.certificate_.flux_only_quasi_newton = !exact_pressure_work;
  candidate.certificate_.allocation_free_apply = true;
  candidate.certificate_.thermophysical_boundary_semantics =
      binding.pressure.thermophysical_boundary_semantics;
  candidate.certificate_.thermophysical_boundary_target =
      binding.pressure.thermophysical_boundary_target;
  candidate.certificate_.thermophysical_boundary_rank_local_binding =
      binding.pressure.thermophysical_boundary_rank_local_binding;
  candidate.certificate_.thermophysical_boundary_collective_lineage =
      binding.pressure.thermophysical_boundary_collective_lineage;
  candidate.certificate_.thermophysical_boundary_rank_local_lineage =
      binding.pressure.thermophysical_boundary_rank_local_lineage;
  if (!candidate.certificate_.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergyPressureFluxBinding};
  }
  out = std::move(candidate);
  certificate = out.certificate_;
  return {};
}

Status PressureEnergyPressureFluxOperator::validate_apply_views(
    FieldView input, FieldView output) const noexcept {
  const Int3 cells = certificate_.linear.local_shape;
  bool aliases =
      detail::field_views_overlap(as_const(input), as_const(output)) ||
      detail::field_views_overlap(as_const(input), temporal_diagonal_) ||
      detail::field_views_overlap(as_const(output), temporal_diagonal_);
  const ConstFaceFieldView faces[]{
      x_pressure_coefficient_, y_pressure_coefficient_,
      z_pressure_coefficient_, target_flux_.x,
      target_flux_.y,          target_flux_.z,
      frozen_face_enthalpy_.x, frozen_face_enthalpy_.y,
      frozen_face_enthalpy_.z};
  for (ConstFaceFieldView face : faces) {
    aliases = aliases ||
              detail::cell_face_views_overlap(as_const(input), face) ||
              detail::cell_face_views_overlap(as_const(output), face);
  }
  const bool exact_pressure_work =
      certificate_.pressure_work_scope ==
      PressureEnergyPressureWorkScope::exact_cartesian;
  if (exact_pressure_work) {
    const ConstFieldView pressure_work_fields[]{
        pressure_work_.target_pressure_perturbation_,
        pressure_work_.target_velocity_, pressure_work_.h_by_a_,
        pressure_work_.r_au_};
    for (ConstFieldView field : pressure_work_fields) {
      aliases = aliases ||
                detail::field_views_overlap(as_const(input), field) ||
                detail::field_views_overlap(as_const(output), field);
    }
  }
  const bool current =
      certificate_.valid() && pressure_boundary_.current() &&
      services_.halo != nullptr && services_.halo->ready() &&
      pressure_boundary_.certificate().geometry == certificate_.geometry &&
      (!exact_pressure_work ||
       (pressure_work_.issuer_ != nullptr &&
        pressure_work_.issuer_
            ->matches_cartesian_pressure_work_linearization(pressure_work_,
                                                             false)));
  const bool valid = current && input.field == services_.solution_field &&
                     detail::valid_cell_view(as_const(input), cells, 0U, 1U,
                                             services_.halo_width) &&
                     detail::valid_cell_view(output, cells, 0U, 1U) && !aliases;
  if (!valid) {
    const Status status{StatusCode::invalid_plan,
                        kPressureEnergyPressureFluxApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  return {};
}

Status PressureEnergyPressureFluxOperator::apply(
    FieldView input, FieldView output) const noexcept {
  failure_ = {};
  Status status = validate_apply_views(input, output);
  if (!status) return status;

  std::array<FieldView, 1U> fields{input};
  HaloTicket ticket;
  status = services_.halo->begin(services_.halo_stage,
                                 {fields.data(), fields.size()}, ticket);
  if (status) {
    status = services_.halo->finish(ticket, {fields.data(), fields.size()});
  }
  if (!status) {
    failure_ = {status, LinearOperatorStatusScope::collective,
                services_.halo->lowest_failing_rank()};
    return status;
  }
  input = fields[0U];
  return apply_after_exchange(input, output);
}

Status PressureEnergyPressureFluxOperator::apply_pressure_energy_shared_input(
    FieldView input, FieldView output,
    const PressureEnergySharedPressureCertificate& certificate,
    const PressureEnergySharedPressureInputCertificate& input_certificate)
    const noexcept {
  failure_ = {};
  Status status = validate_apply_views(input, output);
  if (!status) return status;
  const bool current =
      certificate.valid() && input_certificate.valid() &&
      certificate.energy_consumer_ == this &&
      same_certificate(certificate.energy_pressure_, this->certificate()) &&
      certificate.halo_instance_ == services_.halo->instance_identity() &&
      certificate.halo_stage_ == services_.halo_stage &&
      certificate.pressure_field_ == services_.solution_field &&
      certificate.halo_width_ == services_.halo_width &&
      input_certificate.issuer_ == certificate.regular_issuer_ &&
      input_certificate.shared_contract_ == certificate.collective_contract_ &&
      input_certificate.shared_rank_local_contract_ ==
          certificate.rank_local_contract_ &&
      input_certificate.halo_ghost_revision_ ==
          services_.halo->ghost_revision(input.field) &&
      input_certificate.input_rank_local_binding_ ==
          PressureEnergySharedPressureInputCertificate::local_binding(input) &&
      input_certificate.input_storage_ == input.storage_identity &&
      input_certificate.input_revision_domain_ == input.revision_domain &&
      input_certificate.input_revision_ == input.revision &&
      input_certificate.input_field_ == input.field;
  if (!current) {
    status = {StatusCode::invalid_plan,
              kPressureEnergyPressureFluxApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  return apply_after_exchange(input, output);
}

Status PressureEnergyPressureFluxOperator::apply_after_exchange(
    FieldView input, FieldView output) const noexcept {
  Status status = pressure_boundary_.fill_ghosts(input);
  if (!status) {
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  const Int3 cells = certificate_.linear.local_shape;
  const ConstFieldView pressure = as_const(input);
  const bool exact_pressure_work =
      certificate_.pressure_work_scope ==
      PressureEnergyPressureWorkScope::exact_cartesian;

  Status arithmetic;
  for_each_cell(cells, [&](Int3 cell) {
    const bool active = activity_.cells.size == 0U ||
                        activity_.cells.data[cell_offset(cells, cell)] == 1U;
    if (!active) {
      output.unchecked(cell, 0U) = 0.0;
      return;
    }
    double value =
        temporal_diagonal_.unchecked(cell, 0U) * pressure.unchecked(cell, 0U);
    for (CartesianAxis axis :
         {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
      const ConstFaceFieldView coefficient =
          select_face(x_pressure_coefficient_, y_pressure_coefficient_,
                      z_pressure_coefficient_, axis);
      const ConstFaceFieldView enthalpy =
          select_face(frozen_face_enthalpy_.x, frozen_face_enthalpy_.y,
                      frozen_face_enthalpy_.z, axis);
      Int3 plus = cell;
      Int3 minus_left = cell;
      bool minus_is_local_interior = false;
      bool plus_is_local_interior = false;
      if (axis == CartesianAxis::x) {
        ++plus.x;
        --minus_left.x;
        minus_is_local_interior = cell.x > 0;
        plus_is_local_interior = plus.x < cells.x;
      } else if (axis == CartesianAxis::y) {
        ++plus.y;
        --minus_left.y;
        minus_is_local_interior = cell.y > 0;
        plus_is_local_interior = plus.y < cells.y;
      } else {
        ++plus.z;
        --minus_left.z;
        minus_is_local_interior = cell.z > 0;
        plus_is_local_interior = plus.z < cells.z;
      }
      const bool minus_active = active_face(activity_, axis, cells, cell);
      const bool plus_active = active_face(activity_, axis, cells, plus);
      const double minus_response =
          minus_active
              ? (minus_is_local_interior
                     ? -coefficient.unchecked(cell) *
                           (pressure.unchecked(cell, 0U) -
                            pressure.unchecked(minus_left, 0U))
                     : pressure_boundary_.mass_flux_response_unchecked(
                           pressure, axis, cell,
                           coefficient.unchecked(cell)))
              : 0.0;
      const double plus_response =
          plus_active
              ? (plus_is_local_interior
                     ? -coefficient.unchecked(plus) *
                           (pressure.unchecked(plus, 0U) -
                            pressure.unchecked(cell, 0U))
                     : pressure_boundary_.mass_flux_response_unchecked(
                           pressure, axis, plus,
                           coefficient.unchecked(plus)))
              : 0.0;
      if (plus_active) {
        value += enthalpy.unchecked(plus) * plus_response;
      }
      if (minus_active) {
        value -= enthalpy.unchecked(cell) * minus_response;
      }
    }
    if (exact_pressure_work) {
      const double volume =
          detail::cell_volume(*pressure_work_.kernels_, cell);
      for (std::size_t component = 0U; component < 3U; ++component) {
        const double correction_gradient = pressure_gradient_component(
            *pressure_work_.kernels_, pressure, cell, component);
        const double target_gradient = pressure_gradient_component(
            *pressure_work_.kernels_,
            pressure_work_.target_pressure_perturbation_, cell, component);
        const double velocity =
            pressure_work_.h_by_a_.unchecked(
                cell, static_cast<std::uint8_t>(component));
        const double reciprocal =
            pressure_work_.r_au_.unchecked(
                cell, static_cast<std::uint8_t>(component));
        value += volume *
                 (-velocity * correction_gradient +
                  reciprocal * correction_gradient * target_gradient);
      }
    }
    output.unchecked(cell, 0U) = value;
    if (!std::isfinite(value)) {
      arithmetic = {StatusCode::numerical_failure,
                    kPressureEnergyPressureFluxApply};
    }
  });
  if (!arithmetic) {
    failure_ = {arithmetic, LinearOperatorStatusScope::rank_local, -1};
  }
  return arithmetic;
}

Status form_pressure_energy_thermo_jacobian(
    double pressure_absolute, double enthalpy, const ThermoState& state,
    PressureEnergyThermoJacobian& jacobian) noexcept {
  if (!std::isfinite(pressure_absolute) || pressure_absolute <= 0.0 ||
      !std::isfinite(enthalpy) || !std::isfinite(state.rho) ||
      state.rho <= 0.0 || !std::isfinite(state.temperature) ||
      state.temperature <= 0.0 || !std::isfinite(state.cp) ||
      state.cp <= 0.0 || !std::isfinite(state.drho_dp_hY) ||
      state.drho_dp_hY <= 0.0 || !std::isfinite(state.drho_dh_pY) ||
      state.drho_dh_pY >= 0.0) {
    return {StatusCode::numerical_failure, kPressureEnergyThermo};
  }
  const double eos_density = pressure_absolute * state.drho_dp_hY;
  const double eos_scale =
      std::max({1.0, std::abs(state.rho), std::abs(eos_density)});
  const double eos_tolerance =
      128.0 * std::numeric_limits<double>::epsilon() * eos_scale;
  if (!std::isfinite(eos_density) || eos_density <= 0.0 ||
      std::abs(eos_density - state.rho) > eos_tolerance) {
    return {StatusCode::numerical_failure, kPressureEnergyThermo};
  }
  PressureEnergyThermoJacobian candidate;
  candidate.density = state.rho;
  candidate.drho_dp_hY = state.drho_dp_hY;
  candidate.drho_dh_pY = state.drho_dh_pY;
  candidate.dq_dp_hY = enthalpy * state.drho_dp_hY - 1.0;
  candidate.dq_dh_pY = state.rho + enthalpy * state.drho_dh_pY;
  if (!std::isfinite(candidate.dq_dp_hY) ||
      !std::isfinite(candidate.dq_dh_pY)) {
    return {StatusCode::numerical_failure, kPressureEnergyThermo};
  }
  jacobian = candidate;
  return {};
}

Status linearize_pressure_energy_temporal(
    const PressureEnergyTemporalPoint& point,
    PressureEnergyTemporalLinearization& linearization) noexcept {
  const bool bdf2 = point.bdf.order == 2U;
  const bool valid_current =
      detail::valid_bdf_coefficients(point.bdf) &&
      std::isfinite(point.cell_volume) && point.cell_volume > 0.0 &&
      std::isfinite(point.pressure_absolute) &&
      point.pressure_absolute > 0.0 && std::isfinite(point.density) &&
      point.density > 0.0 && std::isfinite(point.enthalpy) &&
      std::isfinite(point.accepted_pressure_absolute) &&
      point.accepted_pressure_absolute > 0.0 &&
      std::isfinite(point.accepted_density) &&
      point.accepted_density > 0.0 &&
      std::isfinite(point.accepted_enthalpy) &&
      std::isfinite(point.target_thermo.density) &&
      point.target_thermo.density > 0.0 &&
      std::isfinite(point.target_thermo.drho_dp_hY) &&
      point.target_thermo.drho_dp_hY > 0.0 &&
      std::isfinite(point.target_thermo.drho_dh_pY) &&
      point.target_thermo.drho_dh_pY < 0.0 &&
      std::isfinite(point.target_thermo.dq_dp_hY) &&
      std::isfinite(point.target_thermo.dq_dh_pY);
  const bool valid_previous =
      !bdf2 ||
      (std::isfinite(point.previous_pressure_absolute) &&
       point.previous_pressure_absolute > 0.0 &&
       std::isfinite(point.previous_density) &&
       point.previous_density > 0.0 &&
       std::isfinite(point.previous_enthalpy));
  const double density_scale =
      std::max({1.0, std::abs(point.density),
                std::abs(point.target_thermo.density)});
  const bool same_target_density =
      valid_current &&
      std::abs(point.density - point.target_thermo.density) <=
          128.0 * std::numeric_limits<double>::epsilon() * density_scale;
  if (!valid_current || !valid_previous || !same_target_density) {
    return {StatusCode::numerical_failure, kPressureEnergyTemporal};
  }

  PressureEnergyTemporalLinearization candidate;
  candidate.target_q =
      point.density * point.enthalpy - point.pressure_absolute;
  candidate.accepted_q = point.accepted_density * point.accepted_enthalpy -
                         point.accepted_pressure_absolute;
  candidate.previous_q =
      bdf2 ? point.previous_density * point.previous_enthalpy -
                 point.previous_pressure_absolute
           : 0.0;
  candidate.continuity_residual =
      point.cell_volume *
      (point.bdf.a0 * point.density +
       point.bdf.a1 * point.accepted_density +
       point.bdf.a2 * (bdf2 ? point.previous_density : 0.0));
  candidate.energy_residual =
      point.cell_volume *
      (point.bdf.a0 * candidate.target_q +
       point.bdf.a1 * candidate.accepted_q +
       point.bdf.a2 * candidate.previous_q);
  const double temporal_scale = point.cell_volume * point.bdf.a0;
  candidate.continuity_pressure =
      temporal_scale * point.target_thermo.drho_dp_hY;
  candidate.continuity_enthalpy =
      temporal_scale * point.target_thermo.drho_dh_pY;
  candidate.energy_pressure =
      temporal_scale * point.target_thermo.dq_dp_hY;
  candidate.energy_enthalpy =
      temporal_scale * point.target_thermo.dq_dh_pY;
  const double outputs[]{
      candidate.target_q,
      candidate.accepted_q,
      candidate.previous_q,
      candidate.continuity_residual,
      candidate.energy_residual,
      candidate.continuity_pressure,
      candidate.continuity_enthalpy,
      candidate.energy_pressure,
      candidate.energy_enthalpy,
  };
  for (double value : outputs) {
    if (!std::isfinite(value)) {
      return {StatusCode::numerical_failure, kPressureEnergyTemporal};
    }
  }
  linearization = candidate;
  return {};
}

Status select_pressure_energy_globalization(
    const PressureEnergyGlobalizationSample& baseline,
    Span<const PressureEnergyGlobalizationSample> candidates,
    PressureEnergyGlobalizationSelectionCertificate& certificate) noexcept {
  certificate = {};
  const bool valid_baseline =
      baseline.alpha == 0.0 &&
      std::isfinite(baseline.global_normalized_continuity) &&
      baseline.global_normalized_continuity >= 0.0 &&
      std::isfinite(baseline.global_normalized_energy) &&
      baseline.global_normalized_energy >= 0.0 &&
      baseline.thermodynamically_admissible &&
      baseline.state_and_flux_finite &&
      (baseline.corrector == 1U || baseline.corrector == 2U) &&
      baseline.target_time != 0U && baseline.correction_direction != 0U &&
      baseline.state_provenance != 0U &&
      baseline.mass_flux_provenance != 0U;
  if (!valid_baseline || candidates.data == nullptr ||
      candidates.size == 0U ||
      candidates.size > kPressureEnergyGlobalizationCandidateCount) {
    return {StatusCode::invalid_plan, kPressureEnergyGlobalization};
  }

  for (std::size_t index = 0U; index < candidates.size; ++index) {
    const PressureEnergyGlobalizationSample& sample = candidates.data[index];
    const double expected_alpha =
        std::ldexp(1.0, -static_cast<int>(index));
    const bool valid_continuity_norm =
        !std::isfinite(sample.global_normalized_continuity) ||
        sample.global_normalized_continuity >= 0.0;
    const bool valid_energy_norm =
        !std::isfinite(sample.global_normalized_energy) ||
        sample.global_normalized_energy >= 0.0;
    const bool matching_context =
        sample.alpha == expected_alpha &&
        sample.corrector == baseline.corrector &&
        sample.target_time == baseline.target_time &&
        sample.correction_direction == baseline.correction_direction &&
        sample.state_provenance != 0U &&
        sample.mass_flux_provenance != 0U &&
        sample.state_provenance != baseline.state_provenance &&
        sample.mass_flux_provenance != baseline.mass_flux_provenance &&
        valid_continuity_norm && valid_energy_norm;
    if (!matching_context) {
      return {StatusCode::invalid_plan, kPressureEnergyGlobalization};
    }
    for (std::size_t prior = 0U; prior < index; ++prior) {
      const PressureEnergyGlobalizationSample& previous =
          candidates.data[prior];
      if (sample.state_provenance == previous.state_provenance ||
          sample.mass_flux_provenance == previous.mass_flux_provenance) {
        return {StatusCode::invalid_plan, kPressureEnergyGlobalization};
      }
    }
  }

  const double baseline_merit =
      pressure_energy_globalization_merit(baseline);
  const PressureEnergyGlobalizationSample* selected = nullptr;
  double candidate_merit = 0.0;
  double armijo_upper_bound = 0.0;
  std::uint8_t selected_halvings = 0U;
  for (std::size_t index = 0U; index < candidates.size; ++index) {
    const PressureEnergyGlobalizationSample& sample = candidates.data[index];
    const double sample_merit = pressure_energy_globalization_merit(sample);
    const double sample_armijo =
        (1.0 - kPressureEnergyGlobalizationArmijoCoefficient * sample.alpha) *
        baseline_merit;
    const bool acceptable =
        sample.thermodynamically_admissible &&
        sample.state_and_flux_finite && std::isfinite(sample_merit) &&
        sample.global_normalized_continuity >= 0.0 &&
        sample.global_normalized_energy >= 0.0 &&
        sample_merit < baseline_merit && sample_merit <= sample_armijo;
    if (acceptable) {
      selected = &sample;
      candidate_merit = sample_merit;
      armijo_upper_bound = sample_armijo;
      selected_halvings = static_cast<std::uint8_t>(index);
      break;
    }
  }
  if (selected == nullptr) {
    return {StatusCode::rejected_step, kPressureEnergyGlobalization};
  }

  PressureEnergyGlobalizationSelectionCertificate selected_certificate;
  selected_certificate.scope = PressureEnergyGlobalizationScope::
      frozen_momentum_continuity_energy_globalization;
  selected_certificate.alpha = selected->alpha;
  selected_certificate.baseline_normalized_continuity =
      baseline.global_normalized_continuity;
  selected_certificate.baseline_normalized_energy =
      baseline.global_normalized_energy;
  selected_certificate.candidate_normalized_continuity =
      selected->global_normalized_continuity;
  selected_certificate.candidate_normalized_energy =
      selected->global_normalized_energy;
  selected_certificate.baseline_merit = baseline_merit;
  selected_certificate.candidate_merit = candidate_merit;
  selected_certificate.armijo_upper_bound = armijo_upper_bound;
  selected_certificate.target_time = baseline.target_time;
  selected_certificate.correction_direction =
      baseline.correction_direction;
  selected_certificate.baseline_state_provenance =
      baseline.state_provenance;
  selected_certificate.baseline_mass_flux_provenance =
      baseline.mass_flux_provenance;
  selected_certificate.candidate_state_provenance =
      selected->state_provenance;
  selected_certificate.candidate_mass_flux_provenance =
      selected->mass_flux_provenance;
  selected_certificate.corrector = baseline.corrector;
  selected_certificate.selected_halvings = selected_halvings;
  selected_certificate.thermodynamically_admissible = true;
  selected_certificate.state_and_flux_finite = true;
  selected_certificate.strict_merit_decrease = true;
  selected_certificate.armijo_sufficient_decrease = true;
  selected_certificate.full_nonlinear_newton = false;
  selected_certificate.selection_provenance =
      pressure_energy_globalization_selection_provenance(
          baseline, *selected, selected_certificate.selected_halvings);
  if (!selected_certificate.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergyGlobalization};
  }
  certificate = selected_certificate;
  return {};
}

Status select_pressure_energy_extrapolation(
    const PressureEnergyGlobalizationSample& baseline,
    const PressureEnergyGlobalizationSample& candidate,
    PressureEnergyGlobalizationSelectionCertificate& certificate) noexcept {
  certificate = {};
  const bool valid_baseline =
      baseline.alpha == 0.0 &&
      std::isfinite(baseline.global_normalized_continuity) &&
      baseline.global_normalized_continuity >= 0.0 &&
      std::isfinite(baseline.global_normalized_energy) &&
      baseline.global_normalized_energy >= 0.0 &&
      baseline.thermodynamically_admissible &&
      baseline.state_and_flux_finite &&
      (baseline.corrector == 1U || baseline.corrector == 2U) &&
      baseline.target_time != 0U && baseline.correction_direction != 0U &&
      baseline.state_provenance != 0U &&
      baseline.mass_flux_provenance != 0U;
  const bool valid_continuity_norm =
      !std::isfinite(candidate.global_normalized_continuity) ||
      candidate.global_normalized_continuity >= 0.0;
  const bool valid_energy_norm =
      !std::isfinite(candidate.global_normalized_energy) ||
      candidate.global_normalized_energy >= 0.0;
  const bool matching_candidate =
      std::isfinite(candidate.alpha) && candidate.alpha > 1.0 &&
      candidate.alpha <= kPressureEnergyAitkenMaximumAlpha &&
      candidate.corrector == baseline.corrector &&
      candidate.target_time == baseline.target_time &&
      candidate.correction_direction == baseline.correction_direction &&
      candidate.state_provenance != 0U &&
      candidate.mass_flux_provenance != 0U &&
      candidate.state_provenance != baseline.state_provenance &&
      candidate.mass_flux_provenance != baseline.mass_flux_provenance &&
      valid_continuity_norm && valid_energy_norm;
  if (!valid_baseline || !matching_candidate)
    return {StatusCode::invalid_plan, kPressureEnergyGlobalization};
  return certify_pressure_energy_globalization_selection(
      baseline, candidate, 0U, true, certificate);
}

Status PressureEnergyDiagonalOperator::bind(
    const PressureEnergyDiagonalBinding& binding,
    PressureEnergyDiagonalOperator& out,
    PressureEnergyDiagonalCertificate& certificate) noexcept {
  certificate = {};
  const Int3 shape = binding.diagonal.interior;
  if (shape.x <= 0 || shape.y <= 0 || shape.z <= 0 ||
      !valid_identity(binding.identity) ||
      !valid_scalar_view(binding.diagonal, shape) ||
      !valid_activity(binding.activity, shape) ||
      !std::isfinite(binding.inactive_diagonal)) {
    return {StatusCode::invalid_plan, kPressureEnergyDiagonalBinding};
  }

  std::size_t active_cells = 0U;
  bool admissible = true;
  for_each_cell(shape, [&](Int3 cell) {
    if (!active_cell(binding.activity, shape, cell)) return;
    ++active_cells;
    if (!std::isfinite(binding.diagonal.unchecked(cell, 0U))) {
      admissible = false;
    }
  });
  if (!admissible) {
    return {StatusCode::rejected_step, kPressureEnergyDiagonalBinding};
  }

  PressureEnergyDiagonalOperator candidate;
  candidate.diagonal_ = binding.diagonal;
  candidate.activity_ = binding.activity;
  candidate.inactive_diagonal_ = binding.inactive_diagonal;
  candidate.certificate_.linear = {
      binding.identity, diagonal_collective_fingerprint(binding), shape,
      LinearOperatorClass::nonsymmetric};
  candidate.certificate_.diagonal_revision = binding.diagonal.revision;
  candidate.certificate_.inactive_diagonal = binding.inactive_diagonal;
  candidate.certificate_.activity_local_fingerprint =
      binding.activity.local_fingerprint;
  candidate.certificate_.activity_collective_fingerprint =
      binding.activity.collective_fingerprint;
  candidate.certificate_.active_cells = active_cells;
  candidate.certificate_.inactive_cells = cell_count(shape) - active_cells;
  candidate.certificate_.cell_local = true;
  if (!candidate.certificate_.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergyDiagonalBinding};
  }
  out = std::move(candidate);
  certificate = out.certificate_;
  return {};
}

Status PressureEnergyDiagonalOperator::apply(FieldView input,
                                             FieldView output) const noexcept {
  failure_ = {};
  const Int3 shape = certificate_.linear.local_shape;
  if (!certificate_.valid() || !valid_scalar_view(input, shape) ||
      !valid_scalar_view(output, shape) ||
      overlaps(as_const(input), as_const(output)) ||
      overlaps(as_const(input), diagonal_) ||
      overlaps(as_const(output), diagonal_)) {
    const Status status{StatusCode::invalid_plan,
                        kPressureEnergyDiagonalApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  for_each_cell(shape, [&](Int3 cell) {
    const double coefficient = active_cell(activity_, shape, cell)
                                   ? diagonal_.unchecked(cell, 0U)
                                   : inactive_diagonal_;
    output.unchecked(cell, 0U) =
        coefficient * input.unchecked(cell, 0U);
  });
  return {};
}

bool PressureEnergyEnthalpyCertificate::valid() const noexcept {
  return linear.identity.symbolic != 0U && linear.identity.numeric != 0U &&
         linear.identity.hierarchy != 0U && linear.identity.workspace != 0U &&
         linear.identity.fingerprint != 0U &&
         linear.collective_fingerprint != 0U && linear.local_shape.x > 0 &&
         linear.local_shape.y > 0 && linear.local_shape.z > 0 &&
         linear.operator_class == LinearOperatorClass::nonsymmetric &&
         binding_revision != 0U && authority.target_time != 0U &&
         authority.geometry != 0U && authority.numeric_boundary != 0U &&
         authority.thermodynamics != 0U && authority.transport != 0U &&
         authority.equation_semantics != 0U &&
         authority.thermodynamics_semantics != 0U &&
         authority.transport_semantics != 0U && assembled_diagonal != 0U &&
         target_enthalpy != 0U && density_enthalpy_derivative != 0U &&
         heat_capacity != 0U && thermal_conductivity != 0U &&
         enthalpy_diffusivity != 0U && target_flux != 0U &&
         frozen_face_enthalpy != 0U && frozen_reconstruction != 0U &&
         frozen_local_binding != 0U && directional_reconstruction != 0U &&
         directional_branch_authority != 0U && kernels != 0U &&
         geometry_fingerprint != 0U && boundary_semantics != 0U &&
         boundary_rank_local_layout != 0U && halo_stage != 0U &&
         enthalpy_variation_field != 0U && temperature_variation_field != 0U &&
         enthalpy_variation_field != temperature_variation_field &&
         delta_temperature_revision != 0U && delta_temperature_storage != 0U &&
         delta_temperature_revision_domain != 0U &&
         directional_enthalpy_storage != 0U &&
         directional_enthalpy_revision_domain != 0U && halo_instance != 0U &&
         active_cells + inactive_cells ==
             static_cast<std::size_t>(linear.local_shape.x) *
                 static_cast<std::size_t>(linear.local_shape.y) *
                 static_cast<std::size_t>(linear.local_shape.z) &&
         ((activity_local_fingerprint == 0U) ==
          (activity_collective_fingerprint == 0U)) &&
         ((linearization_policy ==
               FrozenConvectionLinearizationPolicy::classical_active_branch &&
           generalized_face_count == 0U) ||
          linearization_policy == FrozenConvectionLinearizationPolicy::
                                      semismooth_generalized_zero_slope) &&
         detail::valid_bdf_coefficients(authority.bdf) &&
         exact_cartesian_spatial_response &&
         exact_temperature_space_conduction && !ibm_spatial_derivative &&
         inactive_rows_identity && inactive_interfaces_zero &&
         allocation_free_apply &&
         (compiled_factored_apply
              ? compiled_numeric_revision != 0U &&
                    compiled_local_binding != 0U
              : compiled_numeric_revision == 0U &&
                    compiled_local_binding == 0U);
}

Status PressureEnergyEnthalpyOperator::bind(
    const PressureEnergyEnthalpyBinding& binding,
    PressureEnergyEnthalpyOperator& out,
    PressureEnergyEnthalpyCertificate& certificate) noexcept {
  certificate = {};
  const Int3 cells = binding.patch.cells;
  std::array<HaloFieldSpec, 2U> halo_fields{{
      {binding.services.enthalpy_variation_field, 2U, 1U},
      {binding.services.temperature_variation_field, 1U, 1U},
  }};
  if (halo_fields[1U].field < halo_fields[0U].field) {
    std::swap(halo_fields[0U], halo_fields[1U]);
  }
  const bool valid_scheme = static_cast<std::uint8_t>(binding.convection) <=
                            static_cast<std::uint8_t>(ConvectionScheme::tvd2);
  const bool valid_policy =
      binding.linearization_policy ==
          FrozenConvectionLinearizationPolicy::classical_active_branch ||
      binding.linearization_policy == FrozenConvectionLinearizationPolicy::
                                          semismooth_generalized_zero_slope;
  const bool pointers_valid =
      binding.geometry != nullptr && binding.kernels != nullptr &&
      binding.boundary != nullptr && binding.services.halo != nullptr;
  const bool geometry_valid =
      pointers_valid && cells.x >= 2 && cells.y >= 2 && cells.z >= 2 &&
      binding.geometry->fingerprint() != 0U &&
      binding.geometry->topology_revision() == binding.authority.geometry &&
      same_shape(binding.kernels->cells(), cells) &&
      same_shape(binding.boundary->local_cells(), cells) &&
      binding.kernels->fingerprint() != 0U;
  const bool boundary_valid =
      pointers_valid && binding.boundary->revision() != 0U &&
      binding.boundary->semantic_fingerprint() != 0U &&
      binding.boundary->local_layout_fingerprint() != 0U &&
      binding.authority.numeric_boundary == binding.boundary->revision() &&
      binding.convection_context.closure == binding.authority.numeric_boundary;
  const bool semantics_valid =
      valid_scheme && valid_policy && valid_identity(binding.identity) &&
      valid_enthalpy_authority(binding.authority) &&
      binding.convection_context.collective_semantics ==
          binding.authority.equation_semantics;
  const bool services_valid =
      pointers_valid && binding.services.communicator != MPI_COMM_NULL &&
      binding.services.halo->ready() && binding.services.halo_stage != 0U &&
      binding.services.enthalpy_variation_field != 0U &&
      binding.services.temperature_variation_field != 0U &&
      binding.services.enthalpy_variation_field !=
          binding.services.temperature_variation_field;
  const bool halo_valid =
      pointers_valid && services_valid && boundary_valid &&
      static_cast<bool>(binding.services.halo->validate_contract(
          binding.services.communicator, binding.patch,
          {halo_fields.data(), halo_fields.size()},
          binding.boundary->halo_topology()));
  const bool cell_views_valid =
      valid_scalar_view(binding.assembled_diagonal, cells) &&
      detail::valid_cell_view(binding.target_enthalpy, cells, 0U, 1U, 2U) &&
      binding.target_enthalpy.field == binding.boundary->enthalpy_field() &&
      valid_scalar_view(binding.density_enthalpy_derivative, cells) &&
      detail::valid_cell_view(binding.heat_capacity, cells, 0U, 1U, 1U) &&
      detail::valid_cell_view(binding.thermal_conductivity, cells, 0U, 1U,
                              1U) &&
      detail::valid_cell_view(binding.enthalpy_diffusivity, cells, 0U, 1U,
                              1U) &&
      detail::valid_cell_view(as_const(binding.workspace.delta_temperature),
                              cells, 0U, 1U, 1U) &&
      binding.workspace.delta_temperature.field ==
          binding.services.temperature_variation_field;
  const bool face_views_valid =
      detail::valid_flux_view(binding.target_flux, cells,
                              binding.target_flux.revision) &&
      !binding.target_flux.certificate.valid() &&
      binding.frozen_face_enthalpy.valid() &&
      valid_face_view(binding.frozen_face_enthalpy.x, CartesianAxis::x,
                      cells) &&
      valid_face_view(binding.frozen_face_enthalpy.y, CartesianAxis::y,
                      cells) &&
      valid_face_view(binding.frozen_face_enthalpy.z, CartesianAxis::z,
                      cells) &&
      valid_face_view(as_const(binding.workspace.directional_enthalpy.x),
                      CartesianAxis::x, cells) &&
      valid_face_view(as_const(binding.workspace.directional_enthalpy.y),
                      CartesianAxis::y, cells) &&
      valid_face_view(as_const(binding.workspace.directional_enthalpy.z),
                      CartesianAxis::z, cells) &&
      valid_pressure_flux_activity(binding.activity, cells);
  const bool compiled_present =
      compiled_enthalpy_workspace_present(binding.workspace.compiled);
  const bool compiled_views_valid =
      !compiled_present ||
      (binding.convection == ConvectionScheme::limited_central2 &&
       valid_compiled_enthalpy_workspace(binding.workspace.compiled, cells));
  if (!pointers_valid || !geometry_valid || !boundary_valid ||
      !semantics_valid || !services_valid || !halo_valid || !cell_views_valid ||
      !face_views_valid || !compiled_views_valid) {
    return {StatusCode::invalid_plan, kPressureEnergyEnthalpyBinding};
  }

  const ConstFieldView cell_views[]{
      binding.assembled_diagonal,
      binding.target_enthalpy,
      binding.density_enthalpy_derivative,
      binding.heat_capacity,
      binding.thermal_conductivity,
      binding.enthalpy_diffusivity,
      as_const(binding.workspace.delta_temperature)};
  const ConstFaceFieldView face_views[]{
      binding.target_flux.x,
      binding.target_flux.y,
      binding.target_flux.z,
      binding.frozen_face_enthalpy.x,
      binding.frozen_face_enthalpy.y,
      binding.frozen_face_enthalpy.z,
      as_const(binding.workspace.directional_enthalpy.x),
      as_const(binding.workspace.directional_enthalpy.y),
      as_const(binding.workspace.directional_enthalpy.z)};
  bool aliases = false;
  for (std::size_t left = 0U; left < std::size(cell_views); ++left) {
    for (std::size_t right = left + 1U; right < std::size(cell_views);
         ++right) {
      aliases = aliases || overlaps(cell_views[left], cell_views[right]);
    }
    for (ConstFaceFieldView face : face_views) {
      aliases =
          aliases || detail::cell_face_views_overlap(cell_views[left], face);
    }
  }
  for (std::size_t left = 0U; left < std::size(face_views); ++left) {
    for (std::size_t right = left + 1U; right < std::size(face_views);
         ++right) {
      aliases = aliases ||
                detail::face_views_overlap(face_views[left], face_views[right]);
    }
  }
  if (compiled_present) {
    const ConstFieldView compiled_cells[]{
        as_const(binding.workspace.compiled.local_diagonal),
        as_const(binding.workspace.compiled.response_stage)};
    const ConstFaceFieldView compiled_faces[]{
        as_const(binding.workspace.compiled.thermal_conductance.x),
        as_const(binding.workspace.compiled.thermal_conductance.y),
        as_const(binding.workspace.compiled.thermal_conductance.z)};
    for (ConstFieldView compiled_cell : compiled_cells) {
      for (ConstFieldView cell : cell_views)
        aliases = aliases || overlaps(compiled_cell, cell);
      for (ConstFaceFieldView face : face_views)
        aliases = aliases ||
                  detail::cell_face_views_overlap(compiled_cell, face);
      for (ConstFaceFieldView face : compiled_faces)
        aliases = aliases ||
                  detail::cell_face_views_overlap(compiled_cell, face);
    }
    aliases = aliases || overlaps(compiled_cells[0U], compiled_cells[1U]);
    for (std::size_t left = 0U; left < std::size(compiled_faces); ++left) {
      for (ConstFieldView cell : cell_views)
        aliases = aliases ||
                  detail::cell_face_views_overlap(cell, compiled_faces[left]);
      for (std::size_t right = left + 1U;
           right < std::size(compiled_faces); ++right) {
        aliases = aliases || detail::face_views_overlap(compiled_faces[left],
                                                        compiled_faces[right]);
      }
      for (ConstFaceFieldView face : face_views)
        aliases = aliases ||
                  detail::face_views_overlap(compiled_faces[left], face);
    }
  }
  if (aliases) {
    return {StatusCode::invalid_plan, kPressureEnergyEnthalpyBinding};
  }

  std::size_t active_cells = 0U;
  bool admissible = true;
  for_each_cell(cells, [&](Int3 cell) {
    if (binding.activity.cells.size != 0U &&
        binding.activity.cells.data[cell_offset(cells, cell)] == 0U) {
      return;
    }
    ++active_cells;
    const double diagonal = binding.assembled_diagonal.unchecked(cell, 0U);
    const double enthalpy = binding.target_enthalpy.unchecked(cell, 0U);
    const double rho_h =
        binding.density_enthalpy_derivative.unchecked(cell, 0U);
    const double cp = binding.heat_capacity.unchecked(cell, 0U);
    const double volume = detail::cell_volume(*binding.kernels, cell);
    const double proxy_diagonal = detail::diffusion_diagonal(
        *binding.kernels, binding.enthalpy_diffusivity, cell);
    const double local = diagonal - proxy_diagonal +
                         binding.authority.bdf.a0 * volume * enthalpy * rho_h;
    if (!std::isfinite(diagonal) || !(diagonal > 0.0) ||
        !std::isfinite(enthalpy) || !std::isfinite(rho_h) || !(rho_h < 0.0) ||
        !std::isfinite(cp) || !(cp > 0.0) || !std::isfinite(volume) ||
        !(volume > 0.0) || !std::isfinite(proxy_diagonal) ||
        !(proxy_diagonal > 0.0) || !std::isfinite(local)) {
      admissible = false;
    }
  });
  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    ConstFaceFieldView selected_flux =
        select_face(binding.target_flux.x, binding.target_flux.y,
                    binding.target_flux.z, axis);
    for (std::int32_t z = 0; z < selected_flux.extents.z; ++z) {
      for (std::int32_t y = 0; y < selected_flux.extents.y; ++y) {
        for (std::int32_t x = 0; x < selected_flux.extents.x; ++x) {
          const Int3 face{x, y, z};
          if (!std::isfinite(selected_flux.unchecked(face))) {
            admissible = false;
          }
          const double conductivity = detail::positive_transmissibility(
              *binding.kernels, binding.thermal_conductivity, axis, face);
          const double proxy = detail::positive_transmissibility(
              *binding.kernels, binding.enthalpy_diffusivity, axis, face);
          if (!std::isfinite(conductivity) || !(conductivity > 0.0) ||
              !std::isfinite(proxy) || !(proxy > 0.0)) {
            admissible = false;
          }
        }
      }
    }
  }
  if (!admissible) {
    return {StatusCode::rejected_step, kPressureEnergyEnthalpyBinding};
  }

  FrozenConvectionFaceDirectionalDerivative directional;
  const Status directional_status =
      differentiate_frozen_cartesian_target_convection_faces(
          *binding.kernels, binding.convection, binding.target_flux,
          binding.target_enthalpy, 0U, binding.convection_context,
          binding.linearization_policy, binding.frozen_face_enthalpy,
          binding.target_enthalpy, 0U, binding.workspace.directional_enthalpy,
          directional);
  if (!directional_status)
    return directional_status;
  if (!directional.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergyEnthalpyBinding};
  }

  FrozenConvectionBranchPlan compiled_branches;
  if (compiled_present) {
    const Status compiled_status = compile_frozen_limited_central2_branches(
        *binding.kernels, binding.target_flux, binding.target_enthalpy, 0U,
        binding.convection_context, binding.linearization_policy,
        binding.frozen_face_enthalpy,
        binding.workspace.compiled.directional_branches, compiled_branches);
    if (!compiled_status) return compiled_status;
    if (!compiled_branches.valid() ||
        compiled_branches.reconstruction != directional.reconstruction ||
        compiled_branches.branch_authority != directional.branch_authority ||
        compiled_branches.generalized_face_count !=
            directional.generalized_face_count ||
        compiled_branches.policy != directional.policy) {
      return {StatusCode::invalid_plan, kPressureEnergyEnthalpyBinding};
    }
    for_each_cell(cells, [&](Int3 cell) {
      double local = 1.0;
      if (binding.activity.cells.size == 0U ||
          binding.activity.cells.data[cell_offset(cells, cell)] != 0U) {
        const double diagonal =
            binding.assembled_diagonal.unchecked(cell, 0U);
        const double enthalpy = binding.target_enthalpy.unchecked(cell, 0U);
        const double rho_h =
            binding.density_enthalpy_derivative.unchecked(cell, 0U);
        const double volume = detail::cell_volume(*binding.kernels, cell);
        const double proxy_diagonal = detail::diffusion_diagonal(
            *binding.kernels, binding.enthalpy_diffusivity, cell);
        local = diagonal - proxy_diagonal +
                binding.authority.bdf.a0 * volume * enthalpy * rho_h;
      }
      binding.workspace.compiled.local_diagonal.unchecked(cell, 0U) = local;
      binding.workspace.compiled.response_stage.unchecked(cell, 0U) = 0.0;
    });
    for (CartesianAxis axis :
         {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
      FaceFieldView output =
          axis == CartesianAxis::x
              ? binding.workspace.compiled.thermal_conductance.x
              : (axis == CartesianAxis::y
                     ? binding.workspace.compiled.thermal_conductance.y
                     : binding.workspace.compiled.thermal_conductance.z);
      for (std::int32_t z = 0; z < output.extents.z; ++z) {
        for (std::int32_t y = 0; y < output.extents.y; ++y) {
          for (std::int32_t x = 0; x < output.extents.x; ++x) {
            const Int3 face{x, y, z};
            output.unchecked(face) = detail::positive_transmissibility(
                *binding.kernels, binding.thermal_conductivity, axis, face);
          }
        }
      }
    }
  }

  PressureEnergyEnthalpyOperator candidate;
  candidate.geometry_ = binding.geometry;
  candidate.kernels_ = binding.kernels;
  candidate.boundary_ = binding.boundary;
  candidate.patch_ = binding.patch;
  candidate.convection_ = binding.convection;
  candidate.services_ = binding.services;
  candidate.assembled_diagonal_ = binding.assembled_diagonal;
  candidate.target_enthalpy_ = binding.target_enthalpy;
  candidate.density_enthalpy_derivative_ = binding.density_enthalpy_derivative;
  candidate.heat_capacity_ = binding.heat_capacity;
  candidate.thermal_conductivity_ = binding.thermal_conductivity;
  candidate.enthalpy_diffusivity_ = binding.enthalpy_diffusivity;
  candidate.target_flux_ = binding.target_flux;
  candidate.convection_context_ = binding.convection_context;
  candidate.frozen_face_enthalpy_ = binding.frozen_face_enthalpy;
  candidate.workspace_ = binding.workspace;
  if (compiled_present) {
    candidate.compiled_directional_branches_ = compiled_branches;
    candidate.compiled_local_diagonal_ =
        as_const(binding.workspace.compiled.local_diagonal);
    candidate.thermal_conductance_x_ =
        as_const(binding.workspace.compiled.thermal_conductance.x);
    candidate.thermal_conductance_y_ =
        as_const(binding.workspace.compiled.thermal_conductance.y);
    candidate.thermal_conductance_z_ =
        as_const(binding.workspace.compiled.thermal_conductance.z);
    candidate.response_stage_ = binding.workspace.compiled.response_stage;
  }
  candidate.activity_ = binding.activity;
  candidate.certificate_.linear = {binding.identity,
                                   enthalpy_collective_fingerprint(binding),
                                   cells, LinearOperatorClass::nonsymmetric};
  candidate.certificate_.binding_revision =
      enthalpy_binding_revision(binding, directional.branch_authority);
  candidate.certificate_.authority = binding.authority;
  candidate.certificate_.assembled_diagonal =
      binding.assembled_diagonal.revision;
  candidate.certificate_.target_enthalpy = binding.target_enthalpy.revision;
  candidate.certificate_.density_enthalpy_derivative =
      binding.density_enthalpy_derivative.revision;
  candidate.certificate_.heat_capacity = binding.heat_capacity.revision;
  candidate.certificate_.thermal_conductivity =
      binding.thermal_conductivity.revision;
  candidate.certificate_.enthalpy_diffusivity =
      binding.enthalpy_diffusivity.revision;
  candidate.certificate_.target_flux = binding.target_flux.revision;
  candidate.certificate_.frozen_face_enthalpy =
      binding.frozen_face_enthalpy.revision;
  candidate.certificate_.frozen_reconstruction =
      binding.frozen_face_enthalpy.reconstruction;
  candidate.certificate_.frozen_local_binding =
      binding.frozen_face_enthalpy.local_binding;
  candidate.certificate_.directional_reconstruction =
      directional.reconstruction;
  candidate.certificate_.directional_branch_authority =
      directional.branch_authority;
  candidate.certificate_.geometry_fingerprint = binding.geometry->fingerprint();
  candidate.certificate_.kernels = binding.kernels->fingerprint();
  candidate.certificate_.boundary_semantics =
      binding.boundary->semantic_fingerprint();
  candidate.certificate_.boundary_rank_local_layout =
      binding.boundary->local_layout_fingerprint();
  candidate.certificate_.halo_stage = binding.services.halo_stage;
  candidate.certificate_.enthalpy_variation_field =
      binding.services.enthalpy_variation_field;
  candidate.certificate_.temperature_variation_field =
      binding.services.temperature_variation_field;
  candidate.certificate_.delta_temperature_revision =
      binding.workspace.delta_temperature.revision;
  candidate.certificate_.delta_temperature_storage =
      binding.workspace.delta_temperature.storage_identity;
  candidate.certificate_.delta_temperature_revision_domain =
      binding.workspace.delta_temperature.revision_domain;
  candidate.certificate_.directional_enthalpy_storage =
      binding.workspace.directional_enthalpy.x.storage_identity;
  candidate.certificate_.directional_enthalpy_revision_domain =
      binding.workspace.directional_enthalpy.x.revision_domain;
  candidate.certificate_.halo_instance =
      binding.services.halo->instance_identity();
  candidate.certificate_.activity_local_fingerprint =
      binding.activity.local_fingerprint;
  candidate.certificate_.activity_collective_fingerprint =
      binding.activity.collective_fingerprint;
  candidate.certificate_.active_cells = active_cells;
  candidate.certificate_.inactive_cells = cell_count(cells) - active_cells;
  candidate.certificate_.generalized_face_count =
      directional.generalized_face_count;
  candidate.certificate_.linearization_policy = binding.linearization_policy;
  candidate.certificate_.exact_cartesian_spatial_response = true;
  candidate.certificate_.exact_temperature_space_conduction = true;
  candidate.certificate_.ibm_spatial_derivative = false;
  candidate.certificate_.inactive_rows_identity = true;
  candidate.certificate_.inactive_interfaces_zero = true;
  candidate.certificate_.allocation_free_apply = true;
  candidate.certificate_.compiled_factored_apply = compiled_present;
  if (compiled_present) {
    candidate.certificate_.compiled_local_binding =
        compiled_enthalpy_local_binding(binding, compiled_branches);
    candidate.certificate_.compiled_numeric_revision =
        compiled_enthalpy_numeric_revision(binding, compiled_branches);
  }
  if (!candidate.certificate_.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergyEnthalpyBinding};
  }
  out = std::move(candidate);
  certificate = out.certificate_;
  return {};
}

Status PressureEnergyEnthalpyOperator::validate_compiled_snapshot()
    const noexcept {
  if (!certificate_.valid() || !certificate_.compiled_factored_apply ||
      geometry_ == nullptr || kernels_ == nullptr || boundary_ == nullptr ||
      !compiled_directional_branches_.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergyEnthalpyApply};
  }
  PressureEnergyEnthalpyBinding current_binding;
  current_binding.geometry = geometry_;
  current_binding.kernels = kernels_;
  current_binding.boundary = boundary_;
  current_binding.patch = patch_;
  current_binding.convection = convection_;
  current_binding.services = services_;
  current_binding.authority = certificate_.authority;
  current_binding.assembled_diagonal = assembled_diagonal_;
  current_binding.target_enthalpy = target_enthalpy_;
  current_binding.density_enthalpy_derivative = density_enthalpy_derivative_;
  current_binding.heat_capacity = heat_capacity_;
  current_binding.thermal_conductivity = thermal_conductivity_;
  current_binding.enthalpy_diffusivity = enthalpy_diffusivity_;
  current_binding.target_flux = target_flux_;
  current_binding.convection_context = convection_context_;
  current_binding.frozen_face_enthalpy = frozen_face_enthalpy_;
  current_binding.workspace = workspace_;
  current_binding.activity = activity_;
  current_binding.identity = certificate_.linear.identity;
  current_binding.linearization_policy = certificate_.linearization_policy;
  if (compiled_enthalpy_local_binding(current_binding,
                                       compiled_directional_branches_) !=
          certificate_.compiled_local_binding ||
      compiled_enthalpy_numeric_revision(current_binding,
                                          compiled_directional_branches_) !=
          certificate_.compiled_numeric_revision) {
    return {StatusCode::invalid_plan, kPressureEnergyEnthalpyApply};
  }

  Status status = validate_frozen_limited_central2_branches(
      *kernels_, target_flux_, target_enthalpy_, 0U, convection_context_,
      frozen_face_enthalpy_, compiled_directional_branches_);
  if (!status) return status;

  const Int3 cells = certificate_.linear.local_shape;
  bool finite = true;
  bool exact = true;
  for_each_cell(cells, [&](Int3 cell) {
    double expected = 1.0;
    if (activity_.cells.size == 0U ||
        activity_.cells.data[cell_offset(cells, cell)] != 0U) {
      const double diagonal = assembled_diagonal_.unchecked(cell, 0U);
      const double enthalpy = target_enthalpy_.unchecked(cell, 0U);
      const double rho_h =
          density_enthalpy_derivative_.unchecked(cell, 0U);
      const double volume = detail::cell_volume(*kernels_, cell);
      const double proxy_diagonal =
          detail::diffusion_diagonal(*kernels_, enthalpy_diffusivity_, cell);
      expected = diagonal - proxy_diagonal +
                 certificate_.authority.bdf.a0 * volume * enthalpy * rho_h;
    }
    const double compiled = compiled_local_diagonal_.unchecked(cell, 0U);
    finite = finite && std::isfinite(expected) && std::isfinite(compiled);
    exact = exact && same_bits(expected, compiled);
  });
  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    const ConstFaceFieldView compiled =
        select_face(thermal_conductance_x_, thermal_conductance_y_,
                    thermal_conductance_z_, axis);
    for (std::int32_t z = 0; z < compiled.extents.z; ++z) {
      for (std::int32_t y = 0; y < compiled.extents.y; ++y) {
        for (std::int32_t x = 0; x < compiled.extents.x; ++x) {
          const Int3 face{x, y, z};
          const double expected = detail::positive_transmissibility(
              *kernels_, thermal_conductivity_, axis, face);
          const double observed = compiled.unchecked(face);
          finite =
              finite && std::isfinite(expected) && std::isfinite(observed);
          exact = exact && same_bits(expected, observed);
        }
      }
    }
  }
  if (!finite)
    return {StatusCode::numerical_failure, kPressureEnergyEnthalpyApply};
  return exact ? Status{}
               : Status{StatusCode::invalid_plan,
                        kPressureEnergyEnthalpyApply};
}

Status PressureEnergyEnthalpyOperator::prepare_repeated_apply(
    PressureEnergyEnthalpyPreparedEpoch& epoch) const noexcept {
  failure_ = {};
  epoch = {};
  const Status status = validate_compiled_snapshot();
  if (!status) {
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  epoch.issuer_ = this;
  epoch.binding_revision_ = certificate_.binding_revision;
  epoch.compiled_numeric_revision_ = certificate_.compiled_numeric_revision;
  return {};
}

Status PressureEnergyEnthalpyOperator::close_repeated_apply(
    PressureEnergyEnthalpyPreparedEpoch& epoch) const noexcept {
  if (!epoch.valid() || epoch.issuer_ != this ||
      epoch.binding_revision_ != certificate_.binding_revision ||
      epoch.compiled_numeric_revision_ !=
          certificate_.compiled_numeric_revision) {
    return {StatusCode::invalid_plan, kPressureEnergyEnthalpyApply};
  }
  epoch = {};
  return {};
}

Status PressureEnergyEnthalpyOperator::apply_prepared(
    FieldView input, FieldView output,
    const PressureEnergyEnthalpyPreparedEpoch& epoch) const noexcept {
  if (!epoch.valid() || epoch.issuer_ != this ||
      epoch.binding_revision_ != certificate_.binding_revision ||
      epoch.compiled_numeric_revision_ !=
          certificate_.compiled_numeric_revision) {
    const Status status{StatusCode::invalid_plan,
                        kPressureEnergyEnthalpyApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  return apply_impl(input, output, true, nullptr);
}

Status PressureEnergyEnthalpyOperator::apply(FieldView input,
                                             FieldView output) const noexcept {
  return apply_impl(input, output, false, nullptr);
}

Status PressureEnergyEnthalpyOperator::enter_schur_prepared_halo()
    const noexcept {
  const bool current =
      certificate_.valid() && certificate_.compiled_factored_apply &&
      services_.halo != nullptr && services_.halo->ready() &&
      services_.halo->instance_identity() == certificate_.halo_instance;
  if (!current) {
    const Status status{StatusCode::invalid_plan,
                        kPressureEnergyEnthalpyApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  const Status status = services_.halo->enter_prepared_epoch();
  if (!status)
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
  return status;
}

Status PressureEnergyEnthalpyOperator::apply_schur_prepared(
    FieldView input, FieldView output, Status& deferred) const noexcept {
  return apply_impl(input, output, true, &deferred);
}

void PressureEnergyEnthalpyOperator::close_schur_prepared_halo(
    bool publish, int lowest_failing_rank) const noexcept {
  if (services_.halo == nullptr) return;
  services_.halo->close_prepared_epoch(publish, lowest_failing_rank);
}

Status PressureEnergyEnthalpyOperator::apply_impl(
    FieldView input, FieldView output, bool snapshot_validated,
    Status* prepared_deferred) const noexcept {
  failure_ = {};
  if (services_.halo == nullptr) {
    const Status status{StatusCode::invalid_plan, kPressureEnergyEnthalpyApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  const Int3 cells = certificate_.linear.local_shape;
  const ConstFieldView bound_fields[]{assembled_diagonal_,
                                      target_enthalpy_,
                                      density_enthalpy_derivative_,
                                      heat_capacity_,
                                      thermal_conductivity_,
                                      enthalpy_diffusivity_,
                                      as_const(workspace_.delta_temperature)};
  const ConstFaceFieldView bound_faces[]{
      target_flux_.x,
      target_flux_.y,
      target_flux_.z,
      frozen_face_enthalpy_.x,
      frozen_face_enthalpy_.y,
      frozen_face_enthalpy_.z,
      as_const(workspace_.directional_enthalpy.x),
      as_const(workspace_.directional_enthalpy.y),
      as_const(workspace_.directional_enthalpy.z)};
  bool aliases = detail::field_views_overlap(as_const(input), as_const(output));
  for (ConstFieldView field : bound_fields) {
    aliases = aliases || detail::field_views_overlap(as_const(input), field) ||
              detail::field_views_overlap(as_const(output), field);
  }
  for (ConstFaceFieldView face : bound_faces) {
    aliases = aliases ||
              detail::cell_face_views_overlap(as_const(input), face) ||
              detail::cell_face_views_overlap(as_const(output), face);
  }
  if (certificate_.compiled_factored_apply) {
    const ConstFieldView compiled_fields[]{compiled_local_diagonal_,
                                           as_const(response_stage_)};
    const ConstFaceFieldView compiled_faces[]{thermal_conductance_x_,
                                              thermal_conductance_y_,
                                              thermal_conductance_z_};
    for (ConstFieldView field : compiled_fields) {
      aliases =
          aliases || detail::field_views_overlap(as_const(input), field) ||
          detail::field_views_overlap(as_const(output), field);
    }
    for (ConstFaceFieldView face : compiled_faces) {
      aliases = aliases ||
                detail::cell_face_views_overlap(as_const(input), face) ||
                detail::cell_face_views_overlap(as_const(output), face);
    }
  }

  PressureEnergyEnthalpyBinding current_binding;
  current_binding.geometry = geometry_;
  current_binding.kernels = kernels_;
  current_binding.boundary = boundary_;
  current_binding.patch = patch_;
  current_binding.convection = convection_;
  current_binding.services = services_;
  current_binding.authority = certificate_.authority;
  current_binding.assembled_diagonal = assembled_diagonal_;
  current_binding.target_enthalpy = target_enthalpy_;
  current_binding.density_enthalpy_derivative = density_enthalpy_derivative_;
  current_binding.heat_capacity = heat_capacity_;
  current_binding.thermal_conductivity = thermal_conductivity_;
  current_binding.enthalpy_diffusivity = enthalpy_diffusivity_;
  current_binding.target_flux = target_flux_;
  current_binding.convection_context = convection_context_;
  current_binding.frozen_face_enthalpy = frozen_face_enthalpy_;
  current_binding.workspace = workspace_;
  current_binding.activity = activity_;
  current_binding.identity = certificate_.linear.identity;
  current_binding.linearization_policy = certificate_.linearization_policy;
  const bool compiled_current =
      !certificate_.compiled_factored_apply ||
      (compiled_directional_branches_.valid() &&
       compiled_directional_branches_.reconstruction ==
           certificate_.directional_reconstruction &&
       compiled_directional_branches_.branch_authority ==
           certificate_.directional_branch_authority &&
       compiled_directional_branches_.generalized_face_count ==
           certificate_.generalized_face_count &&
       compiled_directional_branches_.policy ==
           certificate_.linearization_policy &&
       compiled_enthalpy_local_binding(current_binding,
                                        compiled_directional_branches_) ==
           certificate_.compiled_local_binding &&
       compiled_enthalpy_numeric_revision(current_binding,
                                           compiled_directional_branches_) ==
           certificate_.compiled_numeric_revision);
  const bool current =
      certificate_.valid() && geometry_ != nullptr && kernels_ != nullptr &&
      boundary_ != nullptr && services_.halo->ready() &&
      geometry_->fingerprint() == certificate_.geometry_fingerprint &&
      geometry_->topology_revision() == certificate_.authority.geometry &&
      kernels_->fingerprint() == certificate_.kernels &&
      boundary_->revision() == certificate_.authority.numeric_boundary &&
      boundary_->semantic_fingerprint() == certificate_.boundary_semantics &&
      boundary_->local_layout_fingerprint() ==
          certificate_.boundary_rank_local_layout &&
      services_.halo->instance_identity() == certificate_.halo_instance &&
      services_.halo_stage == certificate_.halo_stage &&
      services_.enthalpy_variation_field ==
          certificate_.enthalpy_variation_field &&
      services_.temperature_variation_field ==
          certificate_.temperature_variation_field &&
      assembled_diagonal_.revision == certificate_.assembled_diagonal &&
      target_enthalpy_.revision == certificate_.target_enthalpy &&
      density_enthalpy_derivative_.revision ==
          certificate_.density_enthalpy_derivative &&
      heat_capacity_.revision == certificate_.heat_capacity &&
      thermal_conductivity_.revision == certificate_.thermal_conductivity &&
      enthalpy_diffusivity_.revision == certificate_.enthalpy_diffusivity &&
      target_flux_.revision == certificate_.target_flux &&
      frozen_face_enthalpy_.revision == certificate_.frozen_face_enthalpy &&
      frozen_face_enthalpy_.reconstruction ==
          certificate_.frozen_reconstruction &&
      frozen_face_enthalpy_.local_binding ==
          certificate_.frozen_local_binding &&
      workspace_.delta_temperature.field ==
          certificate_.temperature_variation_field &&
      workspace_.delta_temperature.revision ==
          certificate_.delta_temperature_revision &&
      workspace_.delta_temperature.storage_identity ==
          certificate_.delta_temperature_storage &&
      workspace_.delta_temperature.revision_domain ==
          certificate_.delta_temperature_revision_domain &&
      workspace_.directional_enthalpy.x.storage_identity ==
          certificate_.directional_enthalpy_storage &&
      workspace_.directional_enthalpy.y.storage_identity ==
          certificate_.directional_enthalpy_storage &&
      workspace_.directional_enthalpy.z.storage_identity ==
          certificate_.directional_enthalpy_storage &&
      workspace_.directional_enthalpy.x.revision_domain ==
          certificate_.directional_enthalpy_revision_domain &&
      workspace_.directional_enthalpy.y.revision_domain ==
          certificate_.directional_enthalpy_revision_domain &&
      workspace_.directional_enthalpy.z.revision_domain ==
          certificate_.directional_enthalpy_revision_domain &&
      activity_.local_fingerprint == certificate_.activity_local_fingerprint &&
      activity_.collective_fingerprint ==
          certificate_.activity_collective_fingerprint &&
      enthalpy_collective_fingerprint(current_binding) ==
          certificate_.linear.collective_fingerprint &&
      enthalpy_binding_revision(current_binding,
                                certificate_.directional_branch_authority) ==
          certificate_.binding_revision &&
      compiled_current;
  const bool views_valid =
      current && input.field == services_.enthalpy_variation_field &&
      detail::valid_cell_view(as_const(input), cells, 0U, 1U, 2U) &&
      valid_scalar_view(output, cells) && !aliases;
  Status prerequisite = views_valid ? Status{}
                                    : Status{StatusCode::invalid_plan,
                                             kPressureEnergyEnthalpyApply};
  if (prerequisite && certificate_.compiled_factored_apply &&
      !snapshot_validated) {
    prerequisite = validate_compiled_snapshot();
  }

  FieldView delta_temperature = workspace_.delta_temperature;
  if (prerequisite) {
    for_each_cell(cells, [&](Int3 cell) {
      const double direction = input.unchecked(cell, 0U);
      const double cp = heat_capacity_.unchecked(cell, 0U);
      const double delta_t = direction / cp;
      if (!std::isfinite(direction) || !std::isfinite(cp) || !(cp > 0.0) ||
          !std::isfinite(delta_t)) {
        prerequisite = {StatusCode::numerical_failure,
                        kPressureEnergyEnthalpyApply};
        return;
      }
      delta_temperature.unchecked(cell, 0U) = delta_t;
    });
  } else if (prepared_deferred != nullptr) {
    // The prepared epoch must still execute the certified payload schedule on
    // every rank.  Publish a deterministic benign payload when the local
    // arithmetic preflight has already failed; the deferred failure prevents
    // any response from being consumed.
    for_each_cell(cells, [&](Int3 cell) {
      delta_temperature.unchecked(cell, 0U) = 0.0;
    });
  }

  std::array<FieldView, 2U> fields{input, delta_temperature};
  HaloTicket ticket;
  Status status;
  if (prepared_deferred != nullptr) {
    retain_first_failure(prerequisite, *prepared_deferred);
    status = services_.halo->begin_prepared(
        services_.halo_stage, {fields.data(), fields.size()},
        *prepared_deferred, ticket);
    if (status) {
      status = services_.halo->finish_prepared(
          ticket, {fields.data(), fields.size()}, *prepared_deferred);
    }
    retain_first_failure(status, *prepared_deferred);
    if (!status || !*prepared_deferred) {
      const Status failure = !status ? status : *prepared_deferred;
      failure_ = {failure, LinearOperatorStatusScope::rank_local, -1};
      return failure;
    }
  } else {
    status = services_.halo->begin(services_.halo_stage,
                                   {fields.data(), fields.size()},
                                   prerequisite, ticket);
    if (status) {
      status = services_.halo->finish(ticket, {fields.data(), fields.size()});
    }
    if (!status) {
      failure_ = {status, LinearOperatorStatusScope::collective,
                  services_.halo->lowest_failing_rank()};
      return status;
    }
  }
  input = fields[0U];
  delta_temperature = fields[1U];

  status = apply_homogeneous_scalar_boundary_ghosts(
      BoundaryStage::enthalpy, *boundary_, target_enthalpy_.field, input, 2U);
  if (status) {
    status = apply_homogeneous_scalar_boundary_ghosts(
        BoundaryStage::enthalpy, *boundary_, target_enthalpy_.field,
        delta_temperature, 1U);
  }
  if (!status) {
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }

  FrozenConvectionFaceDirectionalDerivative directional;
  if (certificate_.compiled_factored_apply) {
    status = apply_frozen_limited_central2_branches(
        *kernels_, compiled_directional_branches_, as_const(input), 0U,
        workspace_.directional_enthalpy);
  } else {
    status = differentiate_frozen_cartesian_target_convection_faces(
        *kernels_, convection_, target_flux_, target_enthalpy_, 0U,
        convection_context_, certificate_.linearization_policy,
        frozen_face_enthalpy_, as_const(input), 0U,
        workspace_.directional_enthalpy, directional);
  }
  if (!status) {
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  if (!certificate_.compiled_factored_apply &&
      (!directional.valid() ||
       directional.reconstruction != certificate_.directional_reconstruction ||
       directional.branch_authority !=
           certificate_.directional_branch_authority ||
       directional.generalized_face_count !=
           certificate_.generalized_face_count ||
       directional.policy != certificate_.linearization_policy)) {
    status = {StatusCode::invalid_plan, kPressureEnergyEnthalpyApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }

  const ConstFieldView direction = as_const(input);
  const ConstFieldView delta_t = as_const(delta_temperature);
  const auto evaluate = [&](Int3 cell) noexcept {
    if (activity_.cells.size != 0U &&
        activity_.cells.data[cell_offset(cells, cell)] == 0U) {
      return direction.unchecked(cell, 0U);
    }
    const double local =
        certificate_.compiled_factored_apply
            ? compiled_local_diagonal_.unchecked(cell, 0U)
            : assembled_diagonal_.unchecked(cell, 0U) -
                  detail::diffusion_diagonal(*kernels_,
                                             enthalpy_diffusivity_, cell) +
                  certificate_.authority.bdf.a0 *
                      detail::cell_volume(*kernels_, cell) *
                      target_enthalpy_.unchecked(cell, 0U) *
                      density_enthalpy_derivative_.unchecked(cell, 0U);
    double value = local * direction.unchecked(cell, 0U);
    const double centre_temperature = delta_t.unchecked(cell, 0U);
    for (CartesianAxis axis :
         {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
      Int3 plus = cell;
      if (axis == CartesianAxis::x)
        ++plus.x;
      else if (axis == CartesianAxis::y)
        ++plus.y;
      else
        ++plus.z;
      const ConstFaceFieldView flux =
          select_face(target_flux_.x, target_flux_.y, target_flux_.z, axis);
      const ConstFaceFieldView face_direction =
          select_face(as_const(workspace_.directional_enthalpy.x),
                      as_const(workspace_.directional_enthalpy.y),
                      as_const(workspace_.directional_enthalpy.z), axis);
      if (active_face(activity_, axis, cells, plus))
        value += flux.unchecked(plus) * face_direction.unchecked(plus);
      if (active_face(activity_, axis, cells, cell))
        value -= flux.unchecked(cell) * face_direction.unchecked(cell);
      for (int side : {-1, 1}) {
        const Int3 face = side < 0 ? cell : plus;
        if (!active_face(activity_, axis, cells, face)) continue;
        Int3 neighbour = cell;
        if (axis == CartesianAxis::x)
          neighbour.x += side;
        else if (axis == CartesianAxis::y)
          neighbour.y += side;
        else
          neighbour.z += side;
        const double conductance =
            certificate_.compiled_factored_apply
                ? select_face(thermal_conductance_x_, thermal_conductance_y_,
                              thermal_conductance_z_, axis)
                      .unchecked(face)
                : detail::positive_transmissibility(
                      *kernels_, thermal_conductivity_, axis, face);
        value += conductance *
                 (centre_temperature - delta_t.unchecked(neighbour, 0U));
      }
    }
    return value;
  };

  bool finite = true;
  if (certificate_.compiled_factored_apply) {
    for_each_cell(cells, [&](Int3 cell) {
      const double value = evaluate(cell);
      if (!std::isfinite(value)) finite = false;
      response_stage_.unchecked(cell, 0U) = value;
    });
  } else {
    for_each_cell(cells, [&](Int3 cell) {
      if (!std::isfinite(evaluate(cell))) finite = false;
    });
  }
  if (!finite) {
    status = {StatusCode::numerical_failure, kPressureEnergyEnthalpyApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  if (certificate_.compiled_factored_apply) {
    for_each_cell(cells, [&](Int3 cell) {
      output.unchecked(cell, 0U) = response_stage_.unchecked(cell, 0U);
    });
  } else {
    for_each_cell(cells, [&](Int3 cell) {
      output.unchecked(cell, 0U) = evaluate(cell);
    });
  }
  return {};
}

Status PressureEnergySchurOperator::bind(
    const PressureEnergySchurBinding& binding,
    PressureEnergySchurOperator& out,
    PressureEnergyJacobianCertificate& certificate) noexcept {
  certificate = {};
  const Int3 shape =
      binding.continuity_enthalpy_diagonal.interior;
  const PressureEnergySchurWorkspace& workspace = binding.workspace;
  const LinearOperatorCertificate cp_certificate =
      binding.continuity_pressure == nullptr
          ? LinearOperatorCertificate{}
          : binding.continuity_pressure->certificate();
  const LinearOperatorCertificate ep_certificate =
      binding.energy_pressure == nullptr
          ? LinearOperatorCertificate{}
          : binding.energy_pressure->certificate();
  const LinearOperatorCertificate eh_certificate =
      binding.energy_enthalpy == nullptr
          ? LinearOperatorCertificate{}
          : binding.energy_enthalpy->certificate();
  const PressureEnergySchurBlockAuthority& block_authority =
      binding.block_authority;
  const bool authority_matches_components =
      block_authority.matches_operators(binding.energy_pressure,
                                        binding.energy_enthalpy) &&
      block_authority.activity_local_fingerprint_ ==
          binding.activity.local_fingerprint &&
      block_authority.activity_collective_fingerprint_ ==
          binding.activity.collective_fingerprint;
  const bool valid_header =
      shape.x > 0 && shape.y > 0 && shape.z > 0 &&
      std::isfinite(binding.scaled_pivot_floor) &&
      binding.scaled_pivot_floor > 0.0 &&
      binding.scaled_pivot_floor < 1.0 &&
      binding.continuity_pressure != nullptr &&
      binding.energy_pressure != nullptr &&
      binding.energy_enthalpy != nullptr &&
      valid_component_certificate(cp_certificate, shape) &&
      valid_component_certificate(ep_certificate, shape) &&
      valid_component_certificate(eh_certificate, shape) &&
      same_identity(cp_certificate.identity, ep_certificate.identity) &&
      same_identity(cp_certificate.identity, eh_certificate.identity) &&
      valid_scalar_view(binding.continuity_enthalpy_diagonal, shape) &&
      valid_scalar_view(binding.continuity_enthalpy_row_scale, shape) &&
      valid_scalar_view(workspace.continuity_response, shape) &&
      valid_scalar_view(workspace.eliminated_enthalpy, shape) &&
      valid_scalar_view(workspace.energy_response, shape) &&
      valid_activity(binding.activity, shape) &&
      authority_matches_components;
  if (!valid_header) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }

  const ConstFieldView diagonal = binding.continuity_enthalpy_diagonal;
  const ConstFieldView scale = binding.continuity_enthalpy_row_scale;
  const ConstFieldView continuity_workspace =
      as_const(workspace.continuity_response);
  const ConstFieldView eliminated_workspace =
      as_const(workspace.eliminated_enthalpy);
  const ConstFieldView energy_workspace = as_const(workspace.energy_response);
  const bool aliases =
      overlaps(diagonal, scale) || overlaps(diagonal, continuity_workspace) ||
      overlaps(diagonal, eliminated_workspace) ||
      overlaps(diagonal, energy_workspace) ||
      overlaps(scale, continuity_workspace) ||
      overlaps(scale, eliminated_workspace) ||
      overlaps(scale, energy_workspace) ||
      overlaps(continuity_workspace, eliminated_workspace) ||
      overlaps(continuity_workspace, energy_workspace) ||
      overlaps(eliminated_workspace, energy_workspace);
  if (aliases) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }

  double minimum_scaled = std::numeric_limits<double>::infinity();
  double maximum_scaled = 0.0;
  std::size_t active_cells = 0U;
  bool admissible = true;
  for_each_cell(shape, [&](Int3 cell) {
    if (!active_cell(binding.activity, shape, cell)) return;
    ++active_cells;
    const double pivot = diagonal.unchecked(cell, 0U);
    const double row_scale = scale.unchecked(cell, 0U);
    const double scaled = std::abs(pivot) / row_scale;
    if (!std::isfinite(pivot) || !std::isfinite(row_scale) ||
        !(row_scale > 0.0) || !std::isfinite(scaled) ||
        scaled < binding.scaled_pivot_floor) {
      admissible = false;
      return;
    }
    minimum_scaled = std::min(minimum_scaled, scaled);
    maximum_scaled = std::max(maximum_scaled, scaled);
  });
  if (!admissible ||
      (active_cells > 0U &&
       (!std::isfinite(minimum_scaled) || !(minimum_scaled > 0.0)))) {
    return {StatusCode::rejected_step, kPressureEnergySchurPivot};
  }
  if (active_cells == 0U) {
    minimum_scaled = 0.0;
    maximum_scaled = 0.0;
  }

  PressureEnergySchurOperator candidate;
  candidate.continuity_pressure_ = binding.continuity_pressure;
  candidate.energy_pressure_ = binding.energy_pressure;
  candidate.energy_enthalpy_ = binding.energy_enthalpy;
  const auto* typed_enthalpy =
      dynamic_cast<const PressureEnergyEnthalpyOperator*>(
          binding.energy_enthalpy);
  if (typed_enthalpy != nullptr &&
      typed_enthalpy->enthalpy_certificate().compiled_factored_apply) {
    candidate.compiled_energy_enthalpy_ = typed_enthalpy;
  }
  candidate.continuity_pressure_certificate_ = cp_certificate;
  candidate.energy_pressure_certificate_ = ep_certificate;
  candidate.energy_enthalpy_certificate_ = eh_certificate;
  candidate.continuity_enthalpy_diagonal_ = diagonal;
  candidate.continuity_enthalpy_row_scale_ = scale;
  candidate.workspace_ = workspace;
  candidate.activity_ = binding.activity;
  candidate.certificate_.schur = {
      cp_certificate.identity,
      nonzero_hash(hash_mix(
          schur_collective_fingerprint(cp_certificate, ep_certificate,
                                       eh_certificate, binding),
          block_authority.collective_fingerprint_)),
      shape,
      LinearOperatorClass::nonsymmetric};
  candidate.certificate_.block_jacobian = nonzero_hash(hash_mix(
      schur_block_revision(cp_certificate, ep_certificate, eh_certificate,
                           binding),
      block_authority.rank_local_revision_));
  candidate.certificate_.block_scope_authority =
      block_authority.rank_local_revision_;
  candidate.certificate_.continuity_enthalpy = diagonal.revision;
  candidate.certificate_.minimum_scaled_abs_c_h = minimum_scaled;
  candidate.certificate_.maximum_scaled_abs_c_h = maximum_scaled;
  candidate.certificate_.active_cells = active_cells;
  candidate.certificate_.inactive_cells = cell_count(shape) - active_cells;
  candidate.certificate_.activity_local_fingerprint =
      binding.activity.local_fingerprint;
  candidate.certificate_.activity_collective_fingerprint =
      binding.activity.collective_fingerprint;
  const auto* typed_energy_pressure =
      dynamic_cast<const PressureEnergyPressureFluxOperator*>(
          binding.energy_pressure);
  if (typed_energy_pressure != nullptr &&
      block_authority.scope_ ==
          PressureEnergySchurBlockScope::exact_cartesian_frozen_spatial) {
    const auto* typed_continuity =
        dynamic_cast<const PressureLinearOperator*>(
            binding.continuity_pressure);
    PressureEnergySharedPressureCertificate shared;
    if (typed_continuity != nullptr &&
        typed_continuity->certify_pressure_energy_shared_halo(
            *typed_continuity,
            PressureEnergySharedPressureScope::cartesian,
            *typed_energy_pressure, shared)) {
      candidate.shared_cartesian_pressure_ = typed_continuity;
      candidate.shared_energy_pressure_ = typed_energy_pressure;
      candidate.certificate_.shared_pressure = shared;
    }
  } else if (
      typed_energy_pressure != nullptr &&
      block_authority.scope_ == PressureEnergySchurBlockScope::
                                      ibm_cartesian_spatial_quasi_newton) {
    const auto* typed_continuity =
        dynamic_cast<const IbmPressureOperator*>(binding.continuity_pressure);
    PressureEnergySharedPressureCertificate shared;
    if (typed_continuity != nullptr &&
        typed_continuity->certify_pressure_energy_shared_halo(
            *typed_energy_pressure, shared)) {
      candidate.shared_cartesian_pressure_ = shared.regular_issuer_;
      candidate.shared_ibm_pressure_ = typed_continuity;
      candidate.shared_energy_pressure_ = typed_energy_pressure;
      candidate.certificate_.shared_pressure = shared;
    }
  }
  // The shared-halo policy changes the executable Schur route even though it
  // preserves the algebraic block identity.  Bind that policy into both the
  // collective and rank-local certificates so a generic three-halo operator
  // cannot be mistaken for the typed two-halo implementation.
  if (candidate.certificate_.shared_pressure.valid()) {
    candidate.certificate_.schur.collective_fingerprint = nonzero_hash(
        hash_mix(candidate.certificate_.schur.collective_fingerprint,
                 candidate.certificate_.shared_pressure
                     .collective_contract_));
    candidate.certificate_.block_jacobian = nonzero_hash(
        hash_mix(candidate.certificate_.block_jacobian,
                 candidate.certificate_.shared_pressure
                     .rank_local_contract_));
  }
  candidate.certificate_.sign_class = PressureEnergySchurSignClass::general;
  if (block_authority.scope_ ==
      PressureEnergySchurBlockScope::exact_cartesian_frozen_spatial) {
    candidate.certificate_.jacobian_scope =
        PressureEnergyJacobianScope::exact_cartesian_frozen_spatial;
  } else if (block_authority.scope_ ==
             PressureEnergySchurBlockScope::
                 ibm_cartesian_spatial_quasi_newton) {
    candidate.certificate_.jacobian_scope =
        PressureEnergyJacobianScope::ibm_cartesian_spatial_quasi_newton;
  } else if (block_authority.scope_ ==
             PressureEnergySchurBlockScope::
                 ibm_double_diagonal_quasi_newton) {
    candidate.certificate_.jacobian_scope =
        PressureEnergyJacobianScope::ibm_double_diagonal_quasi_newton;
  } else {
    candidate.certificate_.jacobian_scope =
        PressureEnergyJacobianScope::generic_algebraic_quasi_newton;
  }
  candidate.certificate_.exact_algebraic_schur = true;
  candidate.certificate_.full_nonlinear_jacobian = false;
  candidate.certificate_.exact_block_equivalent = true;
  candidate.certificate_.cell_local_continuity_enthalpy = true;
  candidate.certificate_.native_mg_preconditioner_only = true;
  if (!candidate.certificate_.valid()) {
    return {StatusCode::invalid_plan, kPressureEnergySchurBinding};
  }
  out = std::move(candidate);
  certificate = out.certificate_;
  return {};
}

Status PressureEnergySchurOperator::prepare_repeated_apply(
    PressureEnergySchurPreparedApplyEpoch& epoch) const noexcept {
  failure_ = {};
  epoch = {};
  if (repeated_apply_active_ || compiled_energy_enthalpy_ == nullptr ||
      shared_cartesian_pressure_ == nullptr ||
      shared_energy_pressure_ == nullptr ||
      !certificate_.shared_pressure.valid() ||
      !certificate_.valid() ||
      !components_current(continuity_pressure_,
                          continuity_pressure_certificate_, energy_pressure_,
                          energy_pressure_certificate_, energy_enthalpy_,
                          energy_enthalpy_certificate_)) {
    const Status status{StatusCode::invalid_plan,
                        kPressureEnergySchurApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  Status status =
      compiled_energy_enthalpy_->prepare_repeated_apply(enthalpy_epoch_);
  if (status) {
    status = shared_cartesian_pressure_
                 ->enter_pressure_energy_shared_prepared_halo(
                     certificate_.shared_pressure);
    pressure_halo_epoch_active_ = static_cast<bool>(status);
  }
  if (status) {
    status = compiled_energy_enthalpy_->enter_schur_prepared_halo();
    enthalpy_halo_epoch_active_ = static_cast<bool>(status);
  }
  if (!status) {
    if (enthalpy_halo_epoch_active_)
      compiled_energy_enthalpy_->close_schur_prepared_halo(false, -1);
    if (pressure_halo_epoch_active_)
      shared_cartesian_pressure_
          ->close_pressure_energy_shared_prepared_halo(false, -1);
    if (enthalpy_epoch_.valid()) {
      (void)compiled_energy_enthalpy_->close_repeated_apply(enthalpy_epoch_);
    }
    pressure_halo_epoch_active_ = false;
    enthalpy_halo_epoch_active_ = false;
    prepared_halo_epoch_ = 0U;
    prepared_pressure_exchange_ = 0U;
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  ++prepared_pressure_exchange_ordinal_;
  if (prepared_pressure_exchange_ordinal_ == 0U)
    ++prepared_pressure_exchange_ordinal_;
  prepared_halo_epoch_ = nonzero_hash(hash_mix(
      certificate_.block_jacobian, prepared_pressure_exchange_ordinal_));
  repeated_apply_active_ = true;
  epoch.issuer_ = this;
  epoch.block_jacobian_ = certificate_.block_jacobian;
  epoch.halo_epoch_ = prepared_halo_epoch_;
  return {};
}

Status PressureEnergySchurOperator::close_repeated_apply(
    PressureEnergySchurPreparedApplyEpoch& epoch) const noexcept {
  return close_repeated_apply(epoch, {}, -1);
}

Status PressureEnergySchurOperator::close_repeated_apply(
    PressureEnergySchurPreparedApplyEpoch& epoch,
    Status globally_consistent_solve_status,
    int lowest_failing_rank) const noexcept {
  const bool valid =
      repeated_apply_active_ && epoch.valid() && epoch.issuer_ == this &&
      epoch.block_jacobian_ == certificate_.block_jacobian &&
      epoch.halo_epoch_ == prepared_halo_epoch_ &&
      compiled_energy_enthalpy_ != nullptr &&
      shared_cartesian_pressure_ != nullptr &&
      pressure_halo_epoch_active_ && enthalpy_halo_epoch_active_ &&
      enthalpy_epoch_.valid();
  if (!valid) {
    if (enthalpy_halo_epoch_active_ && compiled_energy_enthalpy_ != nullptr)
      compiled_energy_enthalpy_->close_schur_prepared_halo(false, -1);
    if (pressure_halo_epoch_active_ && shared_cartesian_pressure_ != nullptr)
      shared_cartesian_pressure_
          ->close_pressure_energy_shared_prepared_halo(false, -1);
    if (enthalpy_epoch_.valid() && compiled_energy_enthalpy_ != nullptr)
      (void)compiled_energy_enthalpy_->close_repeated_apply(enthalpy_epoch_);
    repeated_apply_active_ = false;
    pressure_halo_epoch_active_ = false;
    enthalpy_halo_epoch_active_ = false;
    prepared_halo_epoch_ = 0U;
    prepared_pressure_exchange_ = 0U;
    epoch = {};
    return {StatusCode::invalid_plan, kPressureEnergySchurApply};
  }

  const Status enthalpy_close =
      compiled_energy_enthalpy_->close_repeated_apply(enthalpy_epoch_);
  const bool publish = static_cast<bool>(globally_consistent_solve_status) &&
                       static_cast<bool>(enthalpy_close);
  compiled_energy_enthalpy_->close_schur_prepared_halo(
      publish, publish ? -1 : lowest_failing_rank);
  shared_cartesian_pressure_->close_pressure_energy_shared_prepared_halo(
      publish, publish ? -1 : lowest_failing_rank);
  repeated_apply_active_ = false;
  pressure_halo_epoch_active_ = false;
  enthalpy_halo_epoch_active_ = false;
  prepared_halo_epoch_ = 0U;
  prepared_pressure_exchange_ = 0U;
  epoch = {};
  if (!enthalpy_close)
    return capture_component_failure(*energy_enthalpy_, enthalpy_close,
                                     failure_);
  return {};
}

Status PressureEnergySchurOperator::apply_energy_enthalpy(
    FieldView input, FieldView output) const noexcept {
  if (energy_enthalpy_ == nullptr) {
    const Status status{StatusCode::invalid_plan,
                        kPressureEnergySchurApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  if (!repeated_apply_active_)
    return component_apply(*energy_enthalpy_, input, output, failure_);
  if (compiled_energy_enthalpy_ == nullptr || !enthalpy_epoch_.valid()) {
    const Status status{StatusCode::invalid_plan,
                        kPressureEnergySchurApply};
    failure_ = {status, LinearOperatorStatusScope::rank_local, -1};
    return status;
  }
  const Status status = compiled_energy_enthalpy_->apply_prepared(
      input, output, enthalpy_epoch_);
  if (!status)
    return capture_component_failure(*energy_enthalpy_, status, failure_);
  return {};
}

Status PressureEnergySchurOperator::apply(FieldView pressure,
                                          FieldView output) const noexcept {
  failure_ = {};
  const Int3 shape = certificate_.schur.local_shape;
  const bool components_valid =
      certificate_.valid() && continuity_pressure_ != nullptr &&
      energy_pressure_ != nullptr && energy_enthalpy_ != nullptr &&
      components_current(continuity_pressure_,
                         continuity_pressure_certificate_, energy_pressure_,
                         energy_pressure_certificate_, energy_enthalpy_,
                         energy_enthalpy_certificate_);
  const bool views_valid =
      valid_scalar_view(pressure, shape) &&
      valid_scalar_view(output, shape) &&
      !overlaps(as_const(pressure), as_const(output)) &&
      !overlaps(as_const(pressure), continuity_enthalpy_diagonal_) &&
      !overlaps(as_const(output), continuity_enthalpy_diagonal_) &&
      !overlaps(as_const(pressure), continuity_enthalpy_row_scale_) &&
      !overlaps(as_const(output), continuity_enthalpy_row_scale_) &&
      !overlaps(as_const(pressure),
                as_const(workspace_.continuity_response)) &&
      !overlaps(as_const(pressure),
                as_const(workspace_.eliminated_enthalpy)) &&
      !overlaps(as_const(pressure), as_const(workspace_.energy_response)) &&
      !overlaps(as_const(output),
                as_const(workspace_.continuity_response)) &&
      !overlaps(as_const(output),
                as_const(workspace_.eliminated_enthalpy)) &&
      !overlaps(as_const(output), as_const(workspace_.energy_response));
  Status validation =
      components_valid && views_valid
          ? Status{}
          : Status{StatusCode::invalid_plan, kPressureEnergySchurApply};

  if (repeated_apply_active_) {
    const bool prepared_current =
        prepared_halo_epoch_ != 0U && pressure_halo_epoch_active_ &&
        enthalpy_halo_epoch_active_ && enthalpy_epoch_.valid() &&
        shared_cartesian_pressure_ != nullptr &&
        shared_energy_pressure_ != nullptr &&
        compiled_energy_enthalpy_ != nullptr &&
        certificate_.shared_pressure.valid();
    if (validation && !prepared_current) {
      validation = {StatusCode::invalid_plan, kPressureEnergySchurApply};
    }
    Status deferred = validation;

    ++prepared_pressure_exchange_ordinal_;
    if (prepared_pressure_exchange_ordinal_ == 0U)
      ++prepared_pressure_exchange_ordinal_;
    const RevisionToken exchange_token = nonzero_hash(hash_mix(
        prepared_halo_epoch_, prepared_pressure_exchange_ordinal_));
    FieldView pressure_payload = pressure;
    const bool pressure_exchange_view_valid =
        certificate_.shared_pressure.valid() &&
        pressure.field == certificate_.shared_pressure.pressure_field_ &&
        detail::valid_cell_view(
            as_const(pressure), shape, 0U, 1U,
            certificate_.shared_pressure.halo_width_);
    if (deferred && !pressure_exchange_view_valid) {
      deferred = {StatusCode::invalid_plan, kPressureEnergySchurApply};
    }
    if (!pressure_exchange_view_valid) {
      // Keep the fixed pressure-Halo schedule even for a rank-local invalid
      // Krylov view.  The Schur-owned reach-two eliminated field is a valid,
      // disjoint benign payload after relabelling; deferred prevents either
      // pressure block from consuming it.
      pressure_payload = workspace_.eliminated_enthalpy;
      pressure_payload.field = certificate_.shared_pressure.pressure_field_;
      if (valid_scalar_view(workspace_.eliminated_enthalpy, shape)) {
        for_each_cell(shape, [&](Int3 cell) {
          pressure_payload.unchecked(cell, 0U) = 0.0;
        });
      }
    }
    Status immediate{StatusCode::invalid_plan, kPressureEnergySchurApply};
    if (shared_cartesian_pressure_ != nullptr) {
      immediate = shared_cartesian_pressure_
                      ->exchange_pressure_energy_shared_input_prepared(
                          pressure_payload, certificate_.shared_pressure,
                          deferred);
    }
    retain_first_failure(immediate, deferred);
    if (immediate && deferred)
      prepared_pressure_exchange_ = exchange_token;

    if (deferred && prepared_pressure_exchange_ == exchange_token) {
      Status local = shared_cartesian_pressure_->apply_after_exchange(
          pressure_payload, workspace_.continuity_response);
      if (local && shared_ibm_pressure_ != nullptr) {
        local = shared_ibm_pressure_->decorate_pressure_action(
            pressure_payload, workspace_.continuity_response);
      }
      retain_first_failure(local, deferred);
    }

    if (deferred) {
      for_each_cell(shape, [&](Int3 cell) {
        workspace_.eliminated_enthalpy.unchecked(cell, 0U) =
            active_cell(activity_, shape, cell)
                ? workspace_.continuity_response.unchecked(cell, 0U) /
                      continuity_enthalpy_diagonal_.unchecked(cell, 0U)
                : 0.0;
      });
    } else if (valid_scalar_view(workspace_.eliminated_enthalpy, shape)) {
      for_each_cell(shape, [&](Int3 cell) {
        workspace_.eliminated_enthalpy.unchecked(cell, 0U) = 0.0;
      });
    }

    if (compiled_energy_enthalpy_ != nullptr) {
      immediate = compiled_energy_enthalpy_->apply_schur_prepared(
          workspace_.eliminated_enthalpy, workspace_.energy_response,
          deferred);
    } else {
      immediate = {StatusCode::invalid_plan, kPressureEnergySchurApply};
    }
    retain_first_failure(immediate, deferred);
    if (deferred && prepared_pressure_exchange_ == exchange_token) {
      // E_h uses a disjoint HaloEngine, so the completed pressure exchange is
      // still the exact shared C_p/E_p input.  Stage E_p over the now-dead C_p
      // response and commit the caller output only after both blocks succeed.
      const Status local = shared_energy_pressure_->apply_after_exchange(
          pressure_payload, workspace_.continuity_response);
      retain_first_failure(local, deferred);
    }
    prepared_pressure_exchange_ = 0U;
    if (!deferred) {
      failure_ = {deferred, LinearOperatorStatusScope::rank_local, -1};
      return deferred;
    }

    for_each_cell(shape, [&](Int3 cell) {
      if (active_cell(activity_, shape, cell)) {
        output.unchecked(cell, 0U) =
            workspace_.continuity_response.unchecked(cell, 0U) -
            workspace_.energy_response.unchecked(cell, 0U);
      } else {
        output.unchecked(cell, 0U) = pressure.unchecked(cell, 0U);
      }
    });
    return {};
  }

  if (!validation) {
    failure_ = {validation, LinearOperatorStatusScope::rank_local, -1};
    return validation;
  }

  Status status;
  const bool shared_pressure = certificate_.shared_pressure.valid();
  if (shared_pressure) {
    // Consume E_p immediately after the shared C_p exchange.  Besides saving
    // the second payload halo, this keeps the per-input certificate fail-closed:
    // no intervening exchange can advance the certified pressure ghost view.
    PressureEnergySharedPressureInputCertificate shared_input;
    if (shared_cartesian_pressure_ != nullptr) {
      status = shared_cartesian_pressure_
                   ->exchange_pressure_energy_shared_input(
                       pressure, certificate_.shared_pressure, shared_input);
      if (!status)
        return capture_component_failure(*continuity_pressure_, status,
                                         failure_);
      status = shared_cartesian_pressure_
                   ->apply_pressure_energy_shared_input(
                       pressure, workspace_.continuity_response,
                       certificate_.shared_pressure, shared_input);
    } else if (shared_ibm_pressure_ != nullptr) {
      status = shared_ibm_pressure_->exchange_pressure_energy_shared_input(
          pressure, certificate_.shared_pressure, shared_input);
      if (!status)
        return capture_component_failure(*continuity_pressure_, status,
                                         failure_);
      status = shared_ibm_pressure_->apply_pressure_energy_shared_input(
          pressure, workspace_.continuity_response,
          certificate_.shared_pressure, shared_input);
    } else {
      status = {StatusCode::invalid_plan, kPressureEnergySchurApply};
    }
    if (!status)
      return capture_component_failure(*continuity_pressure_, status,
                                       failure_);
    status = shared_energy_pressure_->apply_pressure_energy_shared_input(
        pressure, output, certificate_.shared_pressure, shared_input);
    if (!status)
      return capture_component_failure(*energy_pressure_, status, failure_);
  } else {
    status = component_apply(*continuity_pressure_, pressure,
                             workspace_.continuity_response, failure_);
    if (!status) return status;
  }
  for_each_cell(shape, [&](Int3 cell) {
    workspace_.eliminated_enthalpy.unchecked(cell, 0U) =
        active_cell(activity_, shape, cell)
            ? workspace_.continuity_response.unchecked(cell, 0U) /
                  continuity_enthalpy_diagonal_.unchecked(cell, 0U)
            : 0.0;
  });
  status = apply_energy_enthalpy(workspace_.eliminated_enthalpy,
                                 workspace_.energy_response);
  if (!status) return status;
  if (!shared_pressure) {
    status = component_apply(*energy_pressure_, pressure, output, failure_);
    if (!status) return status;
  }
  for_each_cell(shape, [&](Int3 cell) {
    if (active_cell(activity_, shape, cell)) {
      output.unchecked(cell, 0U) -=
          workspace_.energy_response.unchecked(cell, 0U);
    } else {
      output.unchecked(cell, 0U) = pressure.unchecked(cell, 0U);
    }
  });
  return {};
}

Status PressureEnergySchurOperator::form_pressure_rhs(
    ConstFieldView continuity_residual, ConstFieldView energy_residual,
    FieldView output) const noexcept {
  failure_ = {};
  const Int3 shape = certificate_.schur.local_shape;
  if (!certificate_.valid() || energy_enthalpy_ == nullptr ||
      !components_current(continuity_pressure_,
                          continuity_pressure_certificate_,
                          energy_pressure_, energy_pressure_certificate_,
                          energy_enthalpy_,
                          energy_enthalpy_certificate_) ||
      !runtime_views_valid(continuity_residual, energy_residual, output,
                           shape, workspace_,
                           continuity_enthalpy_diagonal_,
                           continuity_enthalpy_row_scale_)) {
    return {StatusCode::invalid_plan, kPressureEnergySchurApply};
  }
  for_each_cell(shape, [&](Int3 cell) {
    workspace_.eliminated_enthalpy.unchecked(cell, 0U) =
        active_cell(activity_, shape, cell)
            ? continuity_residual.unchecked(cell, 0U) /
                  continuity_enthalpy_diagonal_.unchecked(cell, 0U)
            : 0.0;
  });
  Status status = apply_energy_enthalpy(workspace_.eliminated_enthalpy,
                                        workspace_.energy_response);
  if (!status) return status;
  for_each_cell(shape, [&](Int3 cell) {
    output.unchecked(cell, 0U) = active_cell(activity_, shape, cell)
                                    ? -energy_residual.unchecked(cell, 0U) +
                                          workspace_.energy_response.unchecked(
                                              cell, 0U)
                                    : 0.0;
  });
  return {};
}

Status PressureEnergySchurOperator::recover_enthalpy(
    ConstFieldView continuity_residual, FieldView pressure_correction,
    FieldView enthalpy_correction) const noexcept {
  failure_ = {};
  const Int3 shape = certificate_.schur.local_shape;
  if (!certificate_.valid() || continuity_pressure_ == nullptr ||
      !components_current(continuity_pressure_,
                          continuity_pressure_certificate_,
                          energy_pressure_, energy_pressure_certificate_,
                          energy_enthalpy_,
                          energy_enthalpy_certificate_) ||
      !runtime_views_valid(continuity_residual,
                           as_const(pressure_correction),
                           enthalpy_correction, shape, workspace_,
                           continuity_enthalpy_diagonal_,
                           continuity_enthalpy_row_scale_)) {
    return {StatusCode::invalid_plan, kPressureEnergySchurApply};
  }
  Status status = component_apply(*continuity_pressure_, pressure_correction,
                                  workspace_.continuity_response, failure_);
  if (!status) return status;
  for_each_cell(shape, [&](Int3 cell) {
    enthalpy_correction.unchecked(cell, 0U) =
        active_cell(activity_, shape, cell)
            ? -(continuity_residual.unchecked(cell, 0U) +
                workspace_.continuity_response.unchecked(cell, 0U)) /
                  continuity_enthalpy_diagonal_.unchecked(cell, 0U)
            : 0.0;
  });
  return {};
}

Status PressureEnergySchurOperator::
    form_pressure_rhs_from_continuity_system_rhs(
        ConstFieldView continuity_system_rhs,
        ConstFieldView energy_residual, FieldView output) const noexcept {
  failure_ = {};
  const Int3 shape = certificate_.schur.local_shape;
  if (!certificate_.valid() || energy_enthalpy_ == nullptr ||
      !components_current(continuity_pressure_,
                          continuity_pressure_certificate_,
                          energy_pressure_, energy_pressure_certificate_,
                          energy_enthalpy_,
                          energy_enthalpy_certificate_) ||
      !runtime_views_valid(continuity_system_rhs, energy_residual, output,
                           shape, workspace_,
                           continuity_enthalpy_diagonal_,
                           continuity_enthalpy_row_scale_)) {
    return {StatusCode::invalid_plan, kPressureEnergySchurApply};
  }
  for_each_cell(shape, [&](Int3 cell) {
    workspace_.eliminated_enthalpy.unchecked(cell, 0U) =
        active_cell(activity_, shape, cell)
            ? continuity_system_rhs.unchecked(cell, 0U) /
                  continuity_enthalpy_diagonal_.unchecked(cell, 0U)
            : 0.0;
  });
  Status status = apply_energy_enthalpy(workspace_.eliminated_enthalpy,
                                        workspace_.energy_response);
  if (!status) return status;
  for_each_cell(shape, [&](Int3 cell) {
    output.unchecked(cell, 0U) =
        active_cell(activity_, shape, cell)
            ? -energy_residual.unchecked(cell, 0U) -
                  workspace_.energy_response.unchecked(cell, 0U)
            : 0.0;
  });
  return {};
}

Status PressureEnergySchurOperator::
    recover_enthalpy_from_continuity_system_rhs(
        ConstFieldView continuity_system_rhs,
        FieldView pressure_correction,
        FieldView enthalpy_correction) const noexcept {
  failure_ = {};
  const Int3 shape = certificate_.schur.local_shape;
  if (!certificate_.valid() || continuity_pressure_ == nullptr ||
      !components_current(continuity_pressure_,
                          continuity_pressure_certificate_,
                          energy_pressure_, energy_pressure_certificate_,
                          energy_enthalpy_,
                          energy_enthalpy_certificate_) ||
      !runtime_views_valid(continuity_system_rhs,
                           as_const(pressure_correction),
                           enthalpy_correction, shape, workspace_,
                           continuity_enthalpy_diagonal_,
                           continuity_enthalpy_row_scale_)) {
    return {StatusCode::invalid_plan, kPressureEnergySchurApply};
  }
  Status status = component_apply(*continuity_pressure_, pressure_correction,
                                  workspace_.continuity_response, failure_);
  if (!status) return status;
  for_each_cell(shape, [&](Int3 cell) {
    enthalpy_correction.unchecked(cell, 0U) =
        active_cell(activity_, shape, cell)
            ? (continuity_system_rhs.unchecked(cell, 0U) -
               workspace_.continuity_response.unchecked(cell, 0U)) /
                  continuity_enthalpy_diagonal_.unchecked(cell, 0U)
            : 0.0;
  });
  return {};
}

}  // namespace hundun::v04
