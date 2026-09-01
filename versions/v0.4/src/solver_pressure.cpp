// SPDX-License-Identifier: Apache-2.0

#include "solver_piso_detail.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_equation_detail.hpp"

#include <cmath>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04::detail {
namespace {

constexpr std::uint32_t kPressureAssembly = 1511U;
constexpr std::uint32_t kPressureNumerical = 1512U;
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
                "pressure identity requires binary64");
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

std::uint64_t mix_view(std::uint64_t hash, ConstFieldView view) noexcept {
  hash = hash_mix(hash, view.field);
  hash = hash_mix(hash, view.revision);
  hash = hash_mix(hash, view.storage_identity);
  return hash_mix(hash, view.revision_domain);
}

bool empty_view(ConstFieldView view) noexcept {
  return view.base == nullptr && view.revision == 0U &&
         view.storage_identity == 0U && view.revision_domain == 0U;
}

bool finite_nonnegative_faces(ConstFaceFieldView view) noexcept {
  for (std::int32_t z = 0; z < view.extents.z; ++z) {
    for (std::int32_t y = 0; y < view.extents.y; ++y) {
      for (std::int32_t x = 0; x < view.extents.x; ++x) {
        const double value = view.unchecked({x, y, z});
        if (!std::isfinite(value) || value < 0.0) {
          return false;
        }
      }
    }
  }
  return true;
}

bool finite_faces(ConstFaceFieldView view) noexcept {
  for (std::int32_t z = 0; z < view.extents.z; ++z) {
    for (std::int32_t y = 0; y < view.extents.y; ++y) {
      for (std::int32_t x = 0; x < view.extents.x; ++x) {
        if (!std::isfinite(view.unchecked({x, y, z}))) {
          return false;
        }
      }
    }
  }
  return true;
}

double face_sum(ConstFaceFieldView x, ConstFaceFieldView y,
                ConstFaceFieldView z, Int3 cell) noexcept {
  return x.unchecked(cell) + x.unchecked({cell.x + 1, cell.y, cell.z}) +
         y.unchecked(cell) + y.unchecked({cell.x, cell.y + 1, cell.z}) +
         z.unchecked(cell) + z.unchecked({cell.x, cell.y, cell.z + 1});
}

double outward_flux(ConstFaceFluxView flux, Int3 cell) noexcept {
  return flux.x.unchecked({cell.x + 1, cell.y, cell.z}) -
             flux.x.unchecked(cell) +
         flux.y.unchecked({cell.x, cell.y + 1, cell.z}) -
             flux.y.unchecked(cell) +
         flux.z.unchecked({cell.x, cell.y, cell.z + 1}) -
             flux.z.unchecked(cell);
}

bool aliases_output(ConstFieldView input,
                    PressureCorrectionSystemView system) noexcept {
  return field_views_overlap(input, as_const(system.diagonal)) ||
         field_views_overlap(input, as_const(system.rhs));
}

bool face_aliases_output(ConstFaceFieldView face,
                         PressureCorrectionSystemView system) noexcept {
  return cell_face_views_overlap(system.diagonal, face) ||
         cell_face_views_overlap(system.rhs, face);
}

}  // namespace

Status assemble_pressure_system_impl(
    const PressureAssemblyBinding& binding,
    const PressureCorrectionInput& input,
    PressureCorrectionSystemView system,
    PressureCorrectionCertificate& certificate) noexcept {
  const Int3 cells = binding.cells;
  const bool bdf2 = input.bdf.order == 2U;
  const ConstFaceFluxView phi = as_const(binding.workspace.phi_h_by_a);
  const ConstFaceFieldView x =
      as_const(binding.workspace.x_pressure_coefficient);
  const ConstFaceFieldView y =
      as_const(binding.workspace.y_pressure_coefficient);
  const ConstFaceFieldView z =
      as_const(binding.workspace.z_pressure_coefficient);
  const bool valid_history =
      valid_cell_view(input.density_accepted, cells, 0U, 1U, 0U) &&
      (bdf2 ? valid_cell_view(input.density_previous, cells, 0U, 1U, 0U)
            : empty_view(input.density_previous));
  const bool valid =
      binding.kernels != nullptr && binding.coupler != 0U &&
      binding.current.valid() && input.intermediate.valid() &&
      input.intermediate.plan == binding.coupler &&
      input.intermediate.plan == binding.current.plan &&
      input.intermediate.dependency == binding.current.dependency &&
      input.intermediate.corrector == binding.current.corrector &&
      input.intermediate.thermophysical_boundary_semantics ==
          binding.current.thermophysical_boundary_semantics &&
      input.intermediate.thermophysical_boundary_target ==
          binding.current.thermophysical_boundary_target &&
      input.intermediate.thermophysical_boundary_rank_local_binding ==
          binding.current.thermophysical_boundary_rank_local_binding &&
      input.intermediate.thermophysical_boundary_collective_lineage ==
          binding.current.thermophysical_boundary_collective_lineage &&
      input.intermediate.thermophysical_boundary_rank_local_lineage ==
          binding.current.thermophysical_boundary_rank_local_lineage &&
      input.intermediate.pressure_energy_refinement ==
          binding.current.pressure_energy_refinement &&
      input.intermediate.pressure_energy_refinement_collective_lineage ==
          binding.current.pressure_energy_refinement_collective_lineage &&
      input.intermediate.pressure_energy_refinement_lineage ==
          binding.current.pressure_energy_refinement_lineage &&
      input.pressure_reference.valid() &&
      input.pressure_reference.plan == binding.pressure_reference_plan &&
      input.pressure_reference.time == input.time && input.time != 0U &&
      binding.thermophysical_context.valid() &&
      binding.thermophysical_context.target_time == input.time &&
      binding.thermophysical_context.pressure_reference ==
          input.pressure_reference.pressure_reference &&
      binding.thermophysical_context.geometry ==
          binding.geometry_fingerprint &&
      binding.thermophysical_context.numeric_boundary ==
          input.numeric_boundary &&
      input.geometry != 0U && input.numeric_boundary != 0U &&
      valid_bdf_coefficients(input.bdf) && input.bdf.a0 > 0.0 &&
      valid_cell_view(input.density_trial, cells, 0U, 1U, 0U) &&
      valid_history &&
      valid_cell_view(input.drho_dp_h_y, cells, 0U, 1U, 0U) &&
      valid_cell_view(system.diagonal, cells, 0U, 1U) &&
      valid_cell_view(system.rhs, cells, 0U, 1U) &&
      !field_views_overlap(as_const(system.diagonal), as_const(system.rhs)) &&
      valid_flux_view(phi, cells, binding.workspace.phi_h_by_a.revision) &&
      valid_equation_face_view(binding.workspace.x_pressure_coefficient,
                               CartesianAxis::x, cells) &&
      valid_equation_face_view(binding.workspace.y_pressure_coefficient,
                               CartesianAxis::y, cells) &&
      valid_equation_face_view(binding.workspace.z_pressure_coefficient,
                               CartesianAxis::z, cells);
  if (!valid) {
    return {StatusCode::invalid_plan, kPressureAssembly};
  }

  const ConstFieldView cell_inputs[]{input.density_trial,
                                     input.density_accepted,
                                     input.drho_dp_h_y};
  for (ConstFieldView view : cell_inputs) {
    if (aliases_output(view, system)) {
      return {StatusCode::invalid_plan, kPressureAssembly};
    }
  }
  if (bdf2 && aliases_output(input.density_previous, system)) {
    return {StatusCode::invalid_plan, kPressureAssembly};
  }
  const ConstFaceFieldView faces[]{x, y, z, phi.x, phi.y, phi.z};
  for (ConstFaceFieldView face : faces) {
    if (face_aliases_output(face, system)) {
      return {StatusCode::invalid_plan, kPressureAssembly};
    }
  }

  const KernelBox box{{0, 0, 0}, cells};
  if (!finite_field_box(input.density_trial, box, 0U, 1U) ||
      !finite_field_box(input.density_accepted, box, 0U, 1U) ||
      (bdf2 &&
       !finite_field_box(input.density_previous, box, 0U, 1U)) ||
      !finite_field_box(input.drho_dp_h_y, box, 0U, 1U) ||
      !finite_nonnegative_faces(x) || !finite_nonnegative_faces(y) ||
      !finite_nonnegative_faces(z) || !finite_faces(phi.x) ||
      !finite_faces(phi.y) || !finite_faces(phi.z)) {
    return {StatusCode::numerical_failure, kPressureNumerical};
  }

  // Preflight the complete algebra so arithmetic overflow cannot partially
  // publish a pressure system.
  for (std::int32_t iz = 0; iz < cells.z; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        const double rho = input.density_trial.unchecked(cell, 0U);
        const double rho_n = input.density_accepted.unchecked(cell, 0U);
        const double rho_nm1 =
            bdf2 ? input.density_previous.unchecked(cell, 0U) : 0.0;
        const double derivative = input.drho_dp_h_y.unchecked(cell, 0U);
        const double volume = cell_volume(*binding.kernels, cell);
        const double density_defect =
            volume * (input.bdf.a0 * rho + input.bdf.a1 * rho_n +
                      input.bdf.a2 * rho_nm1);
        const double diagonal = input.bdf.a0 * volume * derivative +
                                face_sum(x, y, z, cell);
        const double rhs = -(density_defect + outward_flux(phi, cell));
        if (!(rho > 0.0) || !(rho_n > 0.0) ||
            (bdf2 && !(rho_nm1 > 0.0)) || !(derivative > 0.0) ||
            !(volume > 0.0) || !(diagonal > 0.0) ||
            !std::isfinite(density_defect) || !std::isfinite(diagonal) ||
            !std::isfinite(rhs)) {
          return {StatusCode::numerical_failure, kPressureNumerical};
        }
      }
    }
  }

  for (std::int32_t iz = 0; iz < cells.z; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        const double rho_nm1 =
            bdf2 ? input.density_previous.unchecked(cell, 0U) : 0.0;
        const double volume = cell_volume(*binding.kernels, cell);
        const double density_defect =
            volume *
            (input.bdf.a0 * input.density_trial.unchecked(cell, 0U) +
             input.bdf.a1 * input.density_accepted.unchecked(cell, 0U) +
             input.bdf.a2 * rho_nm1);
        system.diagonal.unchecked(cell, 0U) =
            input.bdf.a0 * volume *
                input.drho_dp_h_y.unchecked(cell, 0U) +
            face_sum(x, y, z, cell);
        system.rhs.unchecked(cell, 0U) =
            -(density_defect + outward_flux(phi, cell));
      }
    }
  }

  std::uint64_t state = kFnvOffset;
  state = hash_mix(state, input.intermediate.dependency);
  state = hash_mix(state, input.pressure_reference.pressure_reference);
  state = mix_view(state, input.density_trial);
  state = mix_view(state, input.density_accepted);
  if (bdf2) state = mix_view(state, input.density_previous);
  state = mix_view(state, input.drho_dp_h_y);
  state = hash_mix(state, double_bits(input.bdf.a0));
  state = hash_mix(state, double_bits(input.bdf.a1));
  state = hash_mix(state, double_bits(input.bdf.a2));
  state = hash_mix(state, input.bdf.order);
  state = hash_mix(state, input.time);
  state = hash_mix(state, input.geometry);
  state = hash_mix(state, input.numeric_boundary);
  state = hash_mix(
      state, input.intermediate.thermophysical_boundary_semantics);
  state = hash_mix(state,
                   input.intermediate.thermophysical_boundary_target);
  state = hash_mix(
      state,
      input.intermediate.thermophysical_boundary_rank_local_binding);
  state = hash_mix(
      state,
      input.intermediate.thermophysical_boundary_collective_lineage);
  state = hash_mix(
      state,
      input.intermediate.thermophysical_boundary_rank_local_lineage);
  state = hash_mix(state,
                   input.intermediate.pressure_energy_refinement);
  state = hash_mix(
      state,
      input.intermediate.pressure_energy_refinement_collective_lineage);
  state = hash_mix(
      state, input.intermediate.pressure_energy_refinement_lineage);
  state = mix_view(state, as_const(system.diagonal));
  state = mix_view(state, as_const(system.rhs));
  PressureCorrectionCertificate candidate;
  candidate.plan = binding.coupler;
  candidate.intermediate = input.intermediate.dependency;
  candidate.time = input.time;
  candidate.geometry = input.geometry;
  candidate.numeric_boundary = input.numeric_boundary;
  candidate.state = state == 0U ? 1U : state;
  candidate.corrector = input.intermediate.corrector;
  candidate.thermophysical_boundary_semantics =
      input.intermediate.thermophysical_boundary_semantics;
  candidate.thermophysical_boundary_target =
      input.intermediate.thermophysical_boundary_target;
  candidate.thermophysical_boundary_rank_local_binding =
      input.intermediate.thermophysical_boundary_rank_local_binding;
  candidate.thermophysical_boundary_collective_lineage =
      input.intermediate.thermophysical_boundary_collective_lineage;
  candidate.thermophysical_boundary_rank_local_lineage =
      input.intermediate.thermophysical_boundary_rank_local_lineage;
  candidate.pressure_energy_refinement =
      input.intermediate.pressure_energy_refinement;
  candidate.pressure_energy_refinement_collective_lineage =
      input.intermediate.pressure_energy_refinement_collective_lineage;
  candidate.pressure_energy_refinement_lineage =
      input.intermediate.pressure_energy_refinement_lineage;
  certificate = candidate;
  return {};
}

}  // namespace hundun::v04::detail

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kPressureOperator = 1513U;
constexpr std::uint32_t kPressureOperatorNumerical = 1514U;
constexpr std::uint64_t kOperatorFnvOffset =
    UINT64_C(1469598103934665603);
constexpr std::uint64_t kOperatorFnvPrime = UINT64_C(1099511628211);

std::uint64_t operator_hash(std::uint64_t hash,
                            std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kOperatorFnvPrime;
  return hash;
}

bool same_cells(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
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
         same_cells(left.local_shape, right.local_shape) &&
         left.operator_class == right.operator_class;
}

bool valid_linear_certificate(LinearOperatorCertificate certificate) noexcept {
  return certificate.identity.symbolic != 0U &&
         certificate.identity.numeric != 0U &&
         certificate.identity.hierarchy != 0U &&
         certificate.identity.workspace != 0U &&
         certificate.identity.fingerprint != 0U &&
         certificate.collective_fingerprint != 0U &&
         certificate.local_shape.x > 0 && certificate.local_shape.y > 0 &&
         certificate.local_shape.z > 0;
}

bool same_boundary(PressureCorrectionBoundaryCertificate left,
                   PressureCorrectionBoundaryCertificate right) noexcept {
  return left.valid() && right.valid() && left.semantic == right.semantic &&
         left.rank_local_layout == right.rank_local_layout &&
         left.geometry == right.geometry &&
         left.source_revision == right.source_revision &&
         same_cells(left.global_cells, right.global_cells) &&
         same_cells(left.patch_begin, right.patch_begin) &&
         same_cells(left.local_cells, right.local_cells) &&
         left.geometry_fingerprint == right.geometry_fingerprint;
}

bool valid_coefficient(ConstFaceFieldView view, CartesianAxis axis,
                       Int3 cells) noexcept {
  Int3 expected = cells;
  if (axis == CartesianAxis::x) {
    ++expected.x;
  } else if (axis == CartesianAxis::y) {
    ++expected.y;
  } else {
    ++expected.z;
  }
  detail::FieldStorageInterval interval{};
  return detail::face_storage_interval(view, interval) &&
         view.axis == axis && same_cells(view.extents, expected) &&
         view.storage_identity != 0U && view.revision_domain != 0U;
}

ConstFaceFieldView coefficient_for(ConstFaceFieldView x,
                                   ConstFaceFieldView y,
                                   ConstFaceFieldView z,
                                   CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? x
             : (axis == CartesianAxis::y ? y : z);
}

}  // namespace

PlanFingerprint PressureEnergySharedPressureInputCertificate::local_binding(
    FieldView input) noexcept {
  std::uint64_t hash = kOperatorFnvOffset;
  hash = operator_hash(
      hash, static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(
                input.base)));
  hash = operator_hash(hash, static_cast<std::uint64_t>(input.interior.x));
  hash = operator_hash(hash, static_cast<std::uint64_t>(input.interior.y));
  hash = operator_hash(hash, static_cast<std::uint64_t>(input.interior.z));
  hash = operator_hash(hash, static_cast<std::uint64_t>(input.ghosts.x));
  hash = operator_hash(hash, static_cast<std::uint64_t>(input.ghosts.y));
  hash = operator_hash(hash, static_cast<std::uint64_t>(input.ghosts.z));
  hash = operator_hash(hash, input.components);
  hash = operator_hash(hash, input.stride_y);
  hash = operator_hash(hash, input.stride_z);
  hash = operator_hash(hash, input.component_stride);
  hash = operator_hash(hash, input.replica);
  hash = operator_hash(hash, input.field);
  hash = operator_hash(hash, input.revision);
  hash = operator_hash(hash, input.storage_identity);
  hash = operator_hash(hash, input.revision_domain);
  return hash == 0U ? 1U : hash;
}

struct PressureLinearOperator::Impl {
  PressureCorrectionBoundaryPlan pressure_boundary{};
  HaloEngine* halo{};
  MeshPatch patch{};
  StageId halo_stage{};
  FieldId solution_field{};
  std::uint8_t halo_width{};
  ConstFieldView diagonal{};
  ConstFaceFieldView x{};
  ConstFaceFieldView y{};
  ConstFaceFieldView z{};
  PlanFingerprint coupler{};
  PlanFingerprint fingerprint{};
  LinearOperatorCertificate certificate{};
  LinearOperatorFailureProvenance failure_provenance{};
};

bool PressureEnergySharedPressureCertificate::valid() const noexcept {
  return (scope_ == PressureEnergySharedPressureScope::cartesian ||
          scope_ == PressureEnergySharedPressureScope::ibm_decorated) &&
         valid_linear_certificate(regular_pressure_) &&
         valid_linear_certificate(continuity_pressure_) &&
         valid_linear_certificate(energy_pressure_) &&
         same_identity(regular_pressure_.identity,
                       continuity_pressure_.identity) &&
         same_identity(regular_pressure_.identity, energy_pressure_.identity) &&
         same_cells(regular_pressure_.local_shape,
                    continuity_pressure_.local_shape) &&
         same_cells(regular_pressure_.local_shape,
                    energy_pressure_.local_shape) &&
         collective_contract_ != 0U && rank_local_contract_ != 0U &&
         halo_instance_ != 0U && halo_stage_ != 0U && halo_width_ != 0U &&
         geometry_ != 0U && boundary_semantics_ != 0U &&
         boundary_rank_local_layout_ != 0U && regular_issuer_ != nullptr &&
         continuity_owner_ != nullptr && energy_consumer_ != nullptr &&
         ((scope_ == PressureEnergySharedPressureScope::cartesian &&
           continuity_owner_ == regular_issuer_) ||
          (scope_ == PressureEnergySharedPressureScope::ibm_decorated &&
           continuity_owner_ != regular_issuer_));
}

bool PressureEnergySharedPressureInputCertificate::valid() const noexcept {
  return shared_contract_ != 0U && shared_rank_local_contract_ != 0U &&
         halo_ghost_revision_ != 0U && input_rank_local_binding_ != 0U &&
         input_storage_ != 0U && input_revision_domain_ != 0U &&
         input_revision_ != 0U && issuer_ != nullptr;
}

Status PressureLinearOperator::certify_pressure_energy_shared_halo(
    const LinearOperator& continuity_owner,
    PressureEnergySharedPressureScope scope,
    const PressureEnergyPressureFluxOperator& energy_pressure,
    PressureEnergySharedPressureCertificate& certificate) const noexcept {
  certificate = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPressureOperator};
  }
  const LinearOperatorCertificate regular = this->certificate();
  const LinearOperatorCertificate continuity = continuity_owner.certificate();
  const LinearOperatorCertificate energy = energy_pressure.certificate();
  const PressureCorrectionBoundaryCertificate regular_boundary =
      implementation_->pressure_boundary.certificate();
  const PressureCorrectionBoundaryCertificate energy_boundary =
      energy_pressure.pressure_boundary_.certificate();
  const bool owner_matches_scope =
      (scope == PressureEnergySharedPressureScope::cartesian &&
       &continuity_owner == this) ||
      (scope == PressureEnergySharedPressureScope::ibm_decorated &&
       &continuity_owner != this);
  const bool compatible =
      owner_matches_scope && valid_linear_certificate(regular) &&
      valid_linear_certificate(continuity) &&
      valid_linear_certificate(energy) &&
      regular.operator_class == LinearOperatorClass::spd &&
      continuity.operator_class == LinearOperatorClass::spd &&
      energy.operator_class == LinearOperatorClass::nonsymmetric &&
      same_identity(regular.identity, continuity.identity) &&
      same_identity(regular.identity, energy.identity) &&
      same_cells(regular.local_shape, continuity.local_shape) &&
      same_cells(regular.local_shape, energy.local_shape) &&
      implementation_->pressure_boundary.current() &&
      energy_pressure.certificate_.valid() &&
      energy_pressure.pressure_boundary_.current() &&
      same_boundary(regular_boundary, energy_boundary) &&
      implementation_->halo != nullptr && implementation_->halo->ready() &&
      energy_pressure.services_.halo == implementation_->halo &&
      energy_pressure.services_.halo != nullptr &&
      energy_pressure.services_.halo->ready() &&
      energy_pressure.services_.halo_stage == implementation_->halo_stage &&
      energy_pressure.services_.solution_field ==
          implementation_->solution_field &&
      energy_pressure.services_.halo_width == implementation_->halo_width &&
      same_cells(energy_pressure.patch_.cells, implementation_->patch.cells) &&
      same_cells(energy_pressure.patch_.begin, implementation_->patch.begin);
  if (!compatible) {
    return {StatusCode::invalid_plan, kPressureOperator};
  }

  PressureEnergySharedPressureCertificate candidate;
  candidate.scope_ = scope;
  candidate.regular_pressure_ = regular;
  candidate.continuity_pressure_ = continuity;
  candidate.energy_pressure_ = energy;
  candidate.halo_instance_ = implementation_->halo->instance_identity();
  candidate.halo_stage_ = implementation_->halo_stage;
  candidate.pressure_field_ = implementation_->solution_field;
  candidate.halo_width_ = implementation_->halo_width;
  candidate.geometry_ = regular_boundary.geometry;
  candidate.boundary_semantics_ = regular_boundary.semantic;
  candidate.boundary_rank_local_layout_ = regular_boundary.rank_local_layout;
  candidate.regular_issuer_ = this;
  candidate.continuity_owner_ = &continuity_owner;
  candidate.energy_consumer_ = &energy_pressure;

  std::uint64_t collective = kOperatorFnvOffset;
  collective = operator_hash(collective,
                             static_cast<std::uint8_t>(candidate.scope_));
  collective = operator_hash(collective, regular.collective_fingerprint);
  collective = operator_hash(collective, continuity.collective_fingerprint);
  collective = operator_hash(collective, energy.collective_fingerprint);
  collective = operator_hash(collective, candidate.halo_stage_);
  collective = operator_hash(collective, candidate.pressure_field_);
  collective = operator_hash(collective, candidate.halo_width_);
  collective = operator_hash(collective, candidate.geometry_);
  collective = operator_hash(collective, candidate.boundary_semantics_);
  candidate.collective_contract_ = collective == 0U ? 1U : collective;

  std::uint64_t local = kOperatorFnvOffset;
  local = operator_hash(local, candidate.collective_contract_);
  local = operator_hash(local, candidate.halo_instance_);
  local = operator_hash(local, candidate.boundary_rank_local_layout_);
  local = operator_hash(local, regular.identity.numeric);
  local = operator_hash(local, regular.identity.workspace);
  local = operator_hash(local, energy_pressure.certificate_.binding_revision);
  local = operator_hash(
      local, reinterpret_cast<std::uintptr_t>(&continuity_owner));
  local = operator_hash(
      local, reinterpret_cast<std::uintptr_t>(&energy_pressure));
  candidate.rank_local_contract_ = local == 0U ? 1U : local;
  if (!candidate.valid()) {
    return {StatusCode::invalid_plan, kPressureOperator};
  }
  certificate = candidate;
  return {};
}

Status PressureLinearOperator::exchange_pressure_energy_shared_input(
    FieldView input,
    const PressureEnergySharedPressureCertificate& certificate,
    PressureEnergySharedPressureInputCertificate& input_certificate) const
    noexcept {
  input_certificate = {};
  if (implementation_ != nullptr) implementation_->failure_provenance = {};
  const bool current =
      implementation_ != nullptr && certificate.valid() &&
      certificate.regular_issuer_ == this &&
      same_certificate(certificate.regular_pressure_, this->certificate()) &&
      certificate.halo_instance_ ==
          implementation_->halo->instance_identity() &&
      certificate.halo_stage_ == implementation_->halo_stage &&
      certificate.pressure_field_ == implementation_->solution_field &&
      certificate.halo_width_ == implementation_->halo_width &&
      implementation_->pressure_boundary.current() &&
      input.field == implementation_->solution_field &&
      detail::valid_cell_view(as_const(input), implementation_->patch.cells,
                              0U, 1U, implementation_->halo_width);
  if (!current) {
    return {StatusCode::invalid_plan, kPressureOperator};
  }

  std::array<FieldView, 1U> fields{input};
  HaloTicket ticket;
  Status status = implementation_->halo->begin(
      implementation_->halo_stage, {fields.data(), fields.size()}, ticket);
  if (status) {
    status = implementation_->halo->finish(ticket,
                                            {fields.data(), fields.size()});
  }
  if (!status) {
    implementation_->failure_provenance = {
        status, LinearOperatorStatusScope::collective,
        implementation_->halo->lowest_failing_rank()};
    return status;
  }
  input = fields[0U];
  PressureEnergySharedPressureInputCertificate issued;
  issued.shared_contract_ = certificate.collective_contract_;
  issued.shared_rank_local_contract_ = certificate.rank_local_contract_;
  issued.halo_ghost_revision_ =
      implementation_->halo->ghost_revision(input.field);
  issued.input_rank_local_binding_ =
      PressureEnergySharedPressureInputCertificate::local_binding(input);
  issued.input_storage_ = input.storage_identity;
  issued.input_revision_domain_ = input.revision_domain;
  issued.input_revision_ = input.revision;
  issued.input_field_ = input.field;
  issued.issuer_ = this;
  if (!issued.valid()) {
    return {StatusCode::invalid_plan, kPressureOperator};
  }
  input_certificate = issued;
  return {};
}

Status PressureLinearOperator::enter_pressure_energy_shared_prepared_halo(
    const PressureEnergySharedPressureCertificate& certificate) const
    noexcept {
  if (implementation_ != nullptr) implementation_->failure_provenance = {};
  const bool current =
      implementation_ != nullptr && certificate.valid() &&
      certificate.regular_issuer_ == this &&
      same_certificate(certificate.regular_pressure_, this->certificate()) &&
      certificate.halo_instance_ ==
          implementation_->halo->instance_identity() &&
      certificate.halo_stage_ == implementation_->halo_stage &&
      certificate.pressure_field_ == implementation_->solution_field &&
      certificate.halo_width_ == implementation_->halo_width &&
      implementation_->pressure_boundary.current();
  if (!current) {
    const Status status{StatusCode::invalid_plan, kPressureOperator};
    if (implementation_ != nullptr) {
      implementation_->failure_provenance = {
          status, LinearOperatorStatusScope::rank_local, -1};
    }
    return status;
  }
  const Status status = implementation_->halo->enter_prepared_epoch();
  if (!status) {
    implementation_->failure_provenance = {
        status, LinearOperatorStatusScope::rank_local, -1};
  }
  return status;
}

Status PressureLinearOperator::exchange_pressure_energy_shared_input_prepared(
    FieldView input,
    const PressureEnergySharedPressureCertificate& certificate,
    Status& deferred) const noexcept {
  if (implementation_ != nullptr) implementation_->failure_provenance = {};
  const bool current =
      implementation_ != nullptr && certificate.valid() &&
      certificate.regular_issuer_ == this &&
      same_certificate(certificate.regular_pressure_, this->certificate()) &&
      certificate.halo_instance_ ==
          implementation_->halo->instance_identity() &&
      certificate.halo_stage_ == implementation_->halo_stage &&
      certificate.pressure_field_ == implementation_->solution_field &&
      certificate.halo_width_ == implementation_->halo_width &&
      implementation_->pressure_boundary.current() &&
      input.field == implementation_->solution_field &&
      detail::valid_cell_view(as_const(input), implementation_->patch.cells,
                              0U, 1U, implementation_->halo_width);
  if (deferred && !current) {
    deferred = {StatusCode::invalid_plan, kPressureOperator};
  }

  std::array<FieldView, 1U> fields{input};
  HaloTicket ticket;
  Status immediate = implementation_ == nullptr
                         ? Status{StatusCode::invalid_plan, kPressureOperator}
                         : implementation_->halo->begin_prepared(
                               implementation_->halo_stage,
                               {fields.data(), fields.size()}, deferred,
                               ticket);
  if (immediate) {
    immediate = implementation_->halo->finish_prepared(
        ticket, {fields.data(), fields.size()}, deferred);
  }
  if (!immediate) {
    if (deferred) deferred = immediate;
    if (implementation_ != nullptr) {
      implementation_->failure_provenance = {
          immediate, LinearOperatorStatusScope::rank_local, -1};
    }
  }
  return immediate;
}

void PressureLinearOperator::close_pressure_energy_shared_prepared_halo(
    bool publish, int lowest_failing_rank) const noexcept {
  if (implementation_ == nullptr || implementation_->halo == nullptr) return;
  implementation_->halo->close_prepared_epoch(publish,
                                               lowest_failing_rank);
}

Status PressureLinearOperator::validate_apply_views(FieldView input,
                                                    FieldView output) const
    noexcept {
  if (implementation_ == nullptr ||
      implementation_->certificate.identity.fingerprint == 0U ||
      !implementation_->pressure_boundary.current() ||
      input.field != implementation_->solution_field ||
      !detail::valid_cell_view(as_const(input), implementation_->patch.cells,
                               0U, 1U, implementation_->halo_width) ||
      !detail::valid_cell_view(output, implementation_->patch.cells, 0U, 1U) ||
      detail::field_views_overlap(as_const(input), as_const(output))) {
    return {StatusCode::invalid_plan, kPressureOperator};
  }
  return {};
}

Status PressureLinearOperator::apply_after_exchange(FieldView input,
                                                    FieldView output) const
    noexcept {
  const Status valid = validate_apply_views(input, output);
  if (!valid) return valid;
  const ConstFieldView x = as_const(input);
  const Int3 cells = implementation_->patch.cells;
  Status arithmetic{};
  for (std::int32_t iz = 0; iz < cells.z; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        const double centre = x.unchecked(cell, 0U);
        double value = implementation_->diagonal.unchecked(cell, 0U) * centre;
        for (CartesianAxis axis : {CartesianAxis::x, CartesianAxis::y,
                                   CartesianAxis::z}) {
          const ConstFaceFieldView face = coefficient_for(
              implementation_->x, implementation_->y, implementation_->z,
              axis);
          Int3 plus_face = cell;
          Int3 minus_cell = cell;
          bool minus_is_local_interior = false;
          bool plus_is_local_interior = false;
          if (axis == CartesianAxis::x) {
            ++plus_face.x;
            --minus_cell.x;
            minus_is_local_interior = cell.x > 0;
            plus_is_local_interior = plus_face.x < cells.x;
          } else if (axis == CartesianAxis::y) {
            ++plus_face.y;
            --minus_cell.y;
            minus_is_local_interior = cell.y > 0;
            plus_is_local_interior = plus_face.y < cells.y;
          } else {
            ++plus_face.z;
            --minus_cell.z;
            minus_is_local_interior = cell.z > 0;
            plus_is_local_interior = plus_face.z < cells.z;
          }
          const double minus_neighbor =
              minus_is_local_interior
                  ? x.unchecked(minus_cell, 0U)
                  : implementation_->pressure_boundary
                        .neighbor_value_unchecked(x, cell, axis, -1);
          const double plus_neighbor =
              plus_is_local_interior
                  ? x.unchecked(plus_face, 0U)
                  : implementation_->pressure_boundary
                        .neighbor_value_unchecked(x, cell, axis, 1);
          value -= face.unchecked(cell) * minus_neighbor;
          value -= face.unchecked(plus_face) * plus_neighbor;
        }
        output.unchecked(cell, 0U) = value;
        if (!std::isfinite(value)) {
          arithmetic = {StatusCode::numerical_failure,
                        kPressureOperatorNumerical};
        }
      }
    }
  }
  if (!arithmetic) {
    implementation_->failure_provenance = {
        arithmetic, LinearOperatorStatusScope::rank_local, -1};
  }
  return arithmetic;
}

Status PressureLinearOperator::apply_pressure_energy_shared_input(
    FieldView input, FieldView output,
    const PressureEnergySharedPressureCertificate& certificate,
    const PressureEnergySharedPressureInputCertificate& input_certificate)
    const noexcept {
  if (implementation_ != nullptr) implementation_->failure_provenance = {};
  const bool current =
      implementation_ != nullptr && certificate.valid() &&
      input_certificate.valid() && certificate.regular_issuer_ == this &&
      input_certificate.issuer_ == this &&
      input_certificate.shared_contract_ == certificate.collective_contract_ &&
      input_certificate.shared_rank_local_contract_ ==
          certificate.rank_local_contract_ &&
      input_certificate.halo_ghost_revision_ ==
          implementation_->halo->ghost_revision(input.field) &&
      input_certificate.input_rank_local_binding_ ==
          PressureEnergySharedPressureInputCertificate::local_binding(input) &&
      input_certificate.input_storage_ == input.storage_identity &&
      input_certificate.input_revision_domain_ == input.revision_domain &&
      input_certificate.input_revision_ == input.revision &&
      input_certificate.input_field_ == input.field;
  if (!current) {
    return {StatusCode::invalid_plan, kPressureOperator};
  }
  return apply_after_exchange(input, output);
}

PressureLinearOperator::~PressureLinearOperator() noexcept { release(); }

PressureLinearOperator::PressureLinearOperator(
    PressureLinearOperator&& other) noexcept
    : implementation_(other.implementation_) {
  other.implementation_ = nullptr;
}

PressureLinearOperator& PressureLinearOperator::operator=(
    PressureLinearOperator&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = other.implementation_;
    other.implementation_ = nullptr;
  }
  return *this;
}

void PressureLinearOperator::release() noexcept {
  delete implementation_;
  implementation_ = nullptr;
}

Status PressureLinearOperator::bind_internal(
    const CartesianGeometryPlan& geometry, MeshPatch patch,
    const BoundaryPlan& boundary, PlanFingerprint coupler,
    ConstFaceFieldView x_coefficient,
    ConstFaceFieldView y_coefficient,
    ConstFaceFieldView z_coefficient,
    PressureOperatorServices services,
    PressureCorrectionSystemView system,
    PressureLinearOperator& out) noexcept {
  const Int3 cells = patch.cells;
  const std::array<HaloFieldSpec, 1U> halo_fields{{
      {services.solution_field, services.halo_width, 1U}}};
  const bool valid =
      geometry.fingerprint() != 0U && boundary.semantic_fingerprint() != 0U &&
      coupler != 0U && services.communicator != MPI_COMM_NULL &&
      services.halo != nullptr && services.halo->ready() &&
      services.halo_stage != 0U && services.halo_width != 0U &&
      same_cells(boundary.local_cells(), cells) &&
      detail::valid_cell_view(system.diagonal, cells, 0U, 1U) &&
      detail::valid_cell_view(system.rhs, cells, 0U, 1U) &&
      !detail::field_views_overlap(as_const(system.diagonal),
                                   as_const(system.rhs)) &&
      valid_coefficient(x_coefficient, CartesianAxis::x, cells) &&
      valid_coefficient(y_coefficient, CartesianAxis::y, cells) &&
      valid_coefficient(z_coefficient, CartesianAxis::z, cells) &&
      !detail::cell_face_views_overlap(system.diagonal, x_coefficient) &&
      !detail::cell_face_views_overlap(system.diagonal, y_coefficient) &&
      !detail::cell_face_views_overlap(system.diagonal, z_coefficient) &&
      static_cast<bool>(services.halo->validate_contract(
          services.communicator, patch,
          {halo_fields.data(), halo_fields.size()},
          boundary.halo_topology()));
  if (!valid) {
    return {StatusCode::invalid_plan, kPressureOperator};
  }
  PressureCorrectionBoundaryPlan pressure_boundary;
  const Status pressure_boundary_status =
      PressureCorrectionBoundaryPlan::compile(geometry, patch, boundary,
                                              pressure_boundary);
  if (!pressure_boundary_status) return pressure_boundary_status;
  auto* candidate = new (std::nothrow) Impl;
  if (candidate == nullptr) {
    return {StatusCode::allocation_failure, kPressureOperator};
  }
  candidate->pressure_boundary = std::move(pressure_boundary);
  candidate->halo = services.halo;
  candidate->patch = patch;
  candidate->halo_stage = services.halo_stage;
  candidate->solution_field = services.solution_field;
  candidate->halo_width = services.halo_width;
  candidate->diagonal = as_const(system.diagonal);
  candidate->x = x_coefficient;
  candidate->y = y_coefficient;
  candidate->z = z_coefficient;
  candidate->coupler = coupler;
  std::uint64_t fingerprint = kOperatorFnvOffset;
  fingerprint = operator_hash(
      fingerprint, candidate->pressure_boundary.certificate().semantic);
  fingerprint = operator_hash(
      fingerprint,
      candidate->pressure_boundary.certificate().rank_local_layout);
  fingerprint = operator_hash(fingerprint, coupler);
  fingerprint = operator_hash(fingerprint, system.diagonal.storage_identity);
  fingerprint = operator_hash(fingerprint, x_coefficient.storage_identity);
  fingerprint = operator_hash(fingerprint, y_coefficient.storage_identity);
  fingerprint = operator_hash(fingerprint, z_coefficient.storage_identity);
  fingerprint = operator_hash(fingerprint, services.solution_field);
  candidate->fingerprint = fingerprint == 0U ? 1U : fingerprint;
  out.release();
  out.implementation_ = candidate;
  return {};
}

Status PressureLinearOperator::refresh(
    PressureOperatorRevision revision) noexcept {
  if (implementation_ == nullptr || !revision.pressure.valid() ||
      !implementation_->pressure_boundary.current() ||
      revision.pressure.plan != implementation_->coupler ||
      revision.pressure.geometry !=
          implementation_->pressure_boundary.certificate().geometry ||
      revision.identity.fingerprint == 0U ||
      revision.collective_fingerprint == 0U) {
    return {StatusCode::invalid_plan, kPressureOperator};
  }
  LinearOperatorCertificate candidate;
  candidate.identity = revision.identity;
  candidate.collective_fingerprint = revision.collective_fingerprint;
  candidate.local_shape = implementation_->patch.cells;
  candidate.operator_class = LinearOperatorClass::spd;
  implementation_->certificate = candidate;
  return {};
}

LinearOperatorCertificate PressureLinearOperator::certificate() const
    noexcept {
  return implementation_ == nullptr ? LinearOperatorCertificate{}
                                    : implementation_->certificate;
}

LinearOperatorFailureProvenance
PressureLinearOperator::failure_provenance() const noexcept {
  return implementation_ == nullptr
             ? LinearOperatorFailureProvenance{}
             : implementation_->failure_provenance;
}

Status PressureLinearOperator::apply(FieldView x, FieldView y) const noexcept {
  if (implementation_ != nullptr) {
    implementation_->failure_provenance = {};
  }
  const Status valid = validate_apply_views(x, y);
  if (!valid) return valid;
  std::array<FieldView, 1U> fields{x};
  HaloTicket ticket;
  Status status = implementation_->halo->begin(
      implementation_->halo_stage, {fields.data(), fields.size()}, ticket);
  if (status) {
    status = implementation_->halo->finish(
        ticket, {fields.data(), fields.size()});
  }
  if (!status) {
    implementation_->failure_provenance = {
        status, LinearOperatorStatusScope::collective,
        implementation_->halo->lowest_failing_rank()};
    return status;
  }
  x = fields[0U];
  return apply_after_exchange(x, y);
}

PlanFingerprint PressureLinearOperator::fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->fingerprint;
}

std::uintptr_t PressureLinearOperator::coefficient_storage_address() const
    noexcept {
  return implementation_ == nullptr
             ? 0U
             : reinterpret_cast<std::uintptr_t>(implementation_->diagonal.base);
}

}  // namespace hundun::v04
