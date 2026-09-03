// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"

#include "common_terminal_audit.h"
#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_equation_detail.hpp"
#include "solver_piso_detail.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kPisoPlan = 1501U;
constexpr std::uint32_t kPisoCollective = 1502U;
constexpr std::uint32_t kPisoCoupler = 1503U;
constexpr std::uint32_t kPisoNumerical = 1504U;
constexpr std::uint32_t kPisoSolve = 1505U;
constexpr std::uint32_t kPisoPressureBoundary = 1506U;
constexpr double kPressureCorrectionSafetyFraction = 0.9;
constexpr std::uint64_t kPressureCorrectionAuditSemanticVersion = 4U;
constexpr std::uint64_t kPressureCorrectionBoundarySemanticSchema =
    UINT64_C(0x763034706362736d);
constexpr std::uint64_t kPressureCorrectionBoundaryLayoutSchema =
    UINT64_C(0x7630347063626c79);
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
std::atomic<int> g_post_halo_revalidation_failure_rank{-1};
#endif

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value),
                "PISO identity requires binary64");
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool finite_positive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

bool same_bdf_coefficients(BdfCoefficients left,
                           BdfCoefficients right) noexcept {
  return left.order == right.order &&
         double_bits(left.a0) == double_bits(right.a0) &&
         double_bits(left.a1) == double_bits(right.a1) &&
         double_bits(left.a2) == double_bits(right.a2);
}

PressureCorrectionFaceKind pressure_correction_kind(
    BoundaryKind kind) noexcept {
  if (kind == BoundaryKind::periodic) {
    return PressureCorrectionFaceKind::periodic;
  }
  if (kind == BoundaryKind::pressure_outlet ||
      kind == BoundaryKind::nscbc_outlet ||
      kind == BoundaryKind::static_state_inlet ||
      kind == BoundaryKind::total_state_inlet) {
    return PressureCorrectionFaceKind::homogeneous_dirichlet;
  }
  return PressureCorrectionFaceKind::homogeneous_neumann;
}

MgBoundaryKind mg_boundary_kind(PressureCorrectionFaceKind kind) noexcept {
  if (kind == PressureCorrectionFaceKind::periodic) {
    return MgBoundaryKind::periodic;
  }
  return kind == PressureCorrectionFaceKind::homogeneous_dirichlet
             ? MgBoundaryKind::dirichlet
             : MgBoundaryKind::neumann;
}

bool same_boundaries(MgBoundarySet left, MgBoundarySet right) noexcept {
  return left.x_min == right.x_min && left.x_max == right.x_max &&
         left.y_min == right.y_min && left.y_max == right.y_max &&
         left.z_min == right.z_min && left.z_max == right.z_max;
}

bool same_patch(MeshPatch left, MeshPatch right) noexcept {
  return left.begin.x == right.begin.x && left.begin.y == right.begin.y &&
         left.begin.z == right.begin.z && left.cells.x == right.cells.x &&
         left.cells.y == right.cells.y && left.cells.z == right.cells.z &&
         left.process_grid.x == right.process_grid.x &&
         left.process_grid.y == right.process_grid.y &&
         left.process_grid.z == right.process_grid.z &&
         left.process_coord.x == right.process_coord.x &&
         left.process_coord.y == right.process_coord.y &&
         left.process_coord.z == right.process_coord.z;
}

bool patch_rank(MeshPatch patch, int& rank, int& count) noexcept {
  const std::int64_t gx = patch.process_grid.x;
  const std::int64_t gy = patch.process_grid.y;
  const std::int64_t gz = patch.process_grid.z;
  if (gx <= 0 || gy <= 0 || gz <= 0 || patch.process_coord.x < 0 ||
      patch.process_coord.y < 0 || patch.process_coord.z < 0 ||
      patch.process_coord.x >= gx || patch.process_coord.y >= gy ||
      patch.process_coord.z >= gz ||
      gx > std::numeric_limits<int>::max() / gy ||
      gx * gy > std::numeric_limits<int>::max() / gz) {
    return false;
  }
  const std::int64_t process_count = gx * gy * gz;
  // CartesianGeometryCompiler assigns ranks in x-fastest order.  The halo
  // contract validated by bind() carries this exact patch/rank association.
  const std::int64_t process_rank =
      patch.process_coord.x +
      gx * (patch.process_coord.y + gy * patch.process_coord.z);
  if (process_rank < 0 || process_rank >= process_count) return false;
  rank = static_cast<int>(process_rank);
  count = static_cast<int>(process_count);
  return true;
}

bool encode_global_cell(Int3 index, Int3 cells,
                        std::uint64_t& encoded) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0 || index.x < 0 ||
      index.y < 0 || index.z < 0 || index.x >= cells.x ||
      index.y >= cells.y || index.z >= cells.z) {
    return false;
  }
  const std::uint64_t x = static_cast<std::uint64_t>(index.x);
  const std::uint64_t y = static_cast<std::uint64_t>(index.y);
  const std::uint64_t z = static_cast<std::uint64_t>(index.z);
  const std::uint64_t nx = static_cast<std::uint64_t>(cells.x);
  const std::uint64_t ny = static_cast<std::uint64_t>(cells.y);
  if (z > (std::numeric_limits<std::uint64_t>::max() - y) / ny) return false;
  const std::uint64_t yz = y + ny * z;
  if (yz > (std::numeric_limits<std::uint64_t>::max() - x) / nx)
    return false;
  encoded = x + nx * yz;
  return true;
}

std::int32_t axis_value(Int3 value, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? value.x
             : (axis == CartesianAxis::y ? value.y : value.z);
}

bool valid_pressure_algorithm(LinearAlgorithm algorithm) noexcept {
  return algorithm == LinearAlgorithm::fgmres ||
         algorithm == LinearAlgorithm::bicgstab;
}

bool valid_control(LinearAlgorithm algorithm, MgCorrectionScaling scaling,
                   LinearSolveControl control) noexcept {
  const bool valid_restart =
      algorithm == LinearAlgorithm::fgmres
          ? control.restart >= 2U && control.restart <= 64U
          : control.restart == 0U;
  const bool valid_pair =
      (algorithm == LinearAlgorithm::fgmres &&
       scaling == MgCorrectionScaling::residual_minimizing) ||
      (algorithm == LinearAlgorithm::bicgstab &&
       scaling == MgCorrectionScaling::unit_linear);
  return finite_positive(control.absolute_tolerance) &&
         finite_positive(control.relative_tolerance) &&
         control.relative_tolerance < 1.0 && control.maximum_iterations > 0U &&
         control.true_residual_interval > 0U && valid_restart && valid_pair;
}

PlanFingerprint spec_fingerprint(const EquationPlanSet& equations,
                                 const PisoPlanSpec& spec) noexcept {
  std::uint64_t hash = kFnvOffset;
  hash = hash_mix(hash, UINT64_C(0x7630347069736f32));
  // The cold PISO identity is collective and therefore uses the
  // decomposition-independent equation semantics.  The rank-local equation
  // fingerprint is retained separately for bind-time layout validation.
  hash = hash_mix(hash, equations.semantic_fingerprint());
  hash = hash_mix(hash, equations.thermophysical_predictor().fingerprint());
  hash = hash_mix(hash, equations.pressure_reference().fingerprint());
  if (spec.coupling != CouplingKind::piso)
    hash = hash_mix(hash, static_cast<std::uint8_t>(spec.coupling));
  hash = hash_mix(hash, spec.pressure_correctors);
  hash = hash_mix(hash, spec.pressure_stage);
  hash = hash_mix(hash, spec.final_flux_slot);
  // Keep the legacy default FGMRES tag in the plan identity.  The old
  // implementation hard-coded this value, so omitted-default plans must keep
  // their fingerprint byte-for-byte stable while explicit BiCGStab remains
  // distinguishable below.
  hash = hash_mix(hash, static_cast<std::uint8_t>(spec.pressure_algorithm));
  if (spec.mg_correction_scaling != MgCorrectionScaling::residual_minimizing) {
    hash =
        hash_mix(hash, static_cast<std::uint8_t>(spec.mg_correction_scaling));
  }
  hash = hash_mix(hash, double_bits(spec.pressure_solve.absolute_tolerance));
  hash = hash_mix(hash, double_bits(spec.pressure_solve.relative_tolerance));
  hash = hash_mix(hash, spec.pressure_solve.maximum_iterations);
  hash = hash_mix(hash, spec.pressure_solve.true_residual_interval);
  hash = hash_mix(hash, spec.pressure_solve.restart);
  hash = hash_mix(hash, double_bits(spec.eos_tolerance));
  hash = hash_mix(hash, double_bits(spec.continuity_tolerance));
  if (spec.energy_tolerance != 0.0)
    hash = hash_mix(hash, double_bits(spec.energy_tolerance));
  hash = hash_mix(hash, double_bits(spec.closed_mass_tolerance));
  hash = hash_mix(hash, double_bits(spec.gauge_tolerance));
  return hash == 0U ? 1U : hash;
}

Status collective_status(MPI_Comm communicator, Status local, int rank,
                         int size, int& lowest) noexcept {
  const int candidate = local ? size : rank;
  if (MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
      MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kPisoCollective};
  }
  if (lowest == size) {
    lowest = -1;
    return {};
  }
  std::uint64_t encoded = 0U;
  if (rank == lowest) {
    encoded = (static_cast<std::uint64_t>(local.code) << 32U) | local.detail;
  }
  if (MPI_Bcast(&encoded, 1, MPI_UINT64_T, lowest, communicator) !=
      MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kPisoCollective};
  }
  return {static_cast<StatusCode>(encoded >> 32U),
          static_cast<std::uint32_t>(encoded)};
}

bool same_cells(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool valid_face_workspace(FaceFieldView view, CartesianAxis axis,
                          Int3 cells) noexcept {
  return detail::valid_equation_face_view(view, axis, cells);
}

bool cell_aliases_face(ConstFieldView cell, FaceFieldView face) noexcept {
  return detail::cell_face_views_overlap(cell, face);
}

std::uint64_t mix_view(std::uint64_t hash, ConstFieldView view) noexcept {
  hash = hash_mix(hash, view.field);
  hash = hash_mix(hash, view.revision);
  hash = hash_mix(hash, view.storage_identity);
  return hash_mix(hash, view.revision_domain);
}

std::uint64_t mix_face(std::uint64_t hash,
                       ConstFaceFieldView view) noexcept {
  hash = hash_mix(hash, view.storage_identity);
  hash = hash_mix(hash, view.revision_domain);
  hash = hash_mix(hash, reinterpret_cast<std::uintptr_t>(view.base));
  hash = hash_mix(hash, view.stride_y);
  return hash_mix(hash, view.stride_z);
}

bool same_field_identity(ConstFieldView left, ConstFieldView right) noexcept {
  return left.base == right.base && same_cells(left.interior, right.interior) &&
         same_cells(left.ghosts, right.ghosts) &&
         left.components == right.components && left.stride_y == right.stride_y &&
         left.stride_z == right.stride_z &&
         left.component_stride == right.component_stride &&
         left.replica == right.replica && left.field == right.field &&
         left.revision == right.revision &&
         left.storage_identity == right.storage_identity &&
         left.revision_domain == right.revision_domain;
}

bool same_intermediate_certificate(
    const PisoIntermediateCertificate& left,
    const PisoIntermediateCertificate& right) noexcept {
  return left.plan == right.plan && left.r_au == right.r_au &&
         left.h_by_a == right.h_by_a &&
         left.pressure_face_coefficient == right.pressure_face_coefficient &&
         left.phi_h_by_a == right.phi_h_by_a &&
         left.trial_face_flux == right.trial_face_flux &&
         left.temporal_face_flux == right.temporal_face_flux &&
         left.committed_face_history == right.committed_face_history &&
         left.dependency == right.dependency &&
         left.corrector == right.corrector &&
         left.thermophysical_boundary_semantics ==
             right.thermophysical_boundary_semantics &&
         left.thermophysical_boundary_target ==
             right.thermophysical_boundary_target &&
         left.thermophysical_boundary_rank_local_binding ==
             right.thermophysical_boundary_rank_local_binding &&
         left.thermophysical_boundary_collective_lineage ==
             right.thermophysical_boundary_collective_lineage &&
         left.thermophysical_boundary_rank_local_lineage ==
             right.thermophysical_boundary_rank_local_lineage &&
         left.pressure_energy_refinement ==
             right.pressure_energy_refinement &&
         left.pressure_energy_refinement_collective_lineage ==
             right.pressure_energy_refinement_collective_lineage &&
         left.pressure_energy_refinement_lineage ==
             right.pressure_energy_refinement_lineage;
}

std::uint64_t mix_complete_view_identity(std::uint64_t hash,
                                         ConstFieldView view) noexcept;

PlanFingerprint scalar_reach_one_numeric_fingerprint(
    ConstFieldView field, Int3 cells) noexcept {
  std::uint64_t hash = hash_mix(kFnvOffset, UINT64_C(0x7069726561636831));
  hash = mix_complete_view_identity(hash, field);
  for (std::int32_t z = -1; z <= cells.z; ++z) {
    for (std::int32_t y = -1; y <= cells.y; ++y) {
      for (std::int32_t x = -1; x <= cells.x; ++x) {
        hash = hash_mix(
            hash, double_bits(field.unchecked({x, y, z}, 0U)));
      }
    }
  }
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

PlanFingerprint cartesian_pressure_work_numeric_fingerprint(
    ConstFieldView pressure, ConstFieldView velocity,
    ConstFieldView h_by_a, ConstFieldView r_au, Int3 cells,
    const PisoIntermediateCertificate& intermediate,
    const PressureCorrectionCertificate& pressure_system,
    PlanFingerprint kernels) noexcept {
  std::uint64_t hash = hash_mix(kFnvOffset, UINT64_C(0x65707870776f726b));
  hash = hash_mix(hash,
                  scalar_reach_one_numeric_fingerprint(pressure, cells));
  hash = mix_complete_view_identity(hash, velocity);
  hash = mix_complete_view_identity(hash, h_by_a);
  hash = mix_complete_view_identity(hash, r_au);
  hash = hash_mix(hash, intermediate.dependency);
  hash = hash_mix(hash, intermediate.h_by_a);
  hash = hash_mix(hash, intermediate.r_au);
  hash = hash_mix(hash, pressure_system.state);
  hash = hash_mix(hash, kernels);
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          hash = hash_mix(
              hash, double_bits(velocity.unchecked(cell, component)));
          hash = hash_mix(
              hash, double_bits(h_by_a.unchecked(cell, component)));
          hash = hash_mix(hash,
                          double_bits(r_au.unchecked(cell, component)));
        }
      }
    }
  }
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

struct ThermophysicalBoundaryTokens {
  PlanFingerprint semantics{};
  PlanFingerprint collective_lineage{};
  RevisionToken target{};
  PlanFingerprint rank_local_lineage{};
  PlanFingerprint rank_local_binding{};

  bool valid() const noexcept {
    return semantics != 0U && collective_lineage != 0U && target != 0U &&
           rank_local_lineage != 0U && rank_local_binding != 0U;
  }
};

bool entirely_periodic(const BoundaryPlan& boundary) noexcept {
  for (std::size_t index = 0U; index < 6U; ++index) {
    const BoundaryFacePlan* face = nullptr;
    if (!boundary.face(static_cast<CartesianFace>(index), face) ||
        face == nullptr || !face->periodic) {
      return false;
    }
  }
  return true;
}

bool supported_frozen_open_boundary(const BoundaryPlan& boundary) noexcept {
  for (std::size_t index = 0U; index < 6U; ++index) {
    const BoundaryFacePlan* face = nullptr;
    if (!boundary.face(static_cast<CartesianFace>(index), face) ||
        face == nullptr) {
      return false;
    }
    switch (face->flow_kind) {
      case BoundaryKind::velocity_inlet:
      case BoundaryKind::mass_flow_inlet:
      case BoundaryKind::pressure_outlet:
      case BoundaryKind::no_slip_wall:
      case BoundaryKind::moving_wall:
      case BoundaryKind::slip:
      case BoundaryKind::symmetry:
      case BoundaryKind::periodic:
        break;
      default:
        return false;
    }
  }
  return true;
}

CartesianFace cartesian_face(CartesianAxis axis, bool high) noexcept {
  const std::size_t axis_index = static_cast<std::size_t>(axis);
  return static_cast<CartesianFace>(2U * axis_index + (high ? 1U : 0U));
}

std::uint64_t mix_complete_view_identity(std::uint64_t hash,
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
  hash = hash_mix(hash, view.replica);
  hash = hash_mix(hash, view.field);
  hash = hash_mix(hash, view.revision);
  hash = hash_mix(hash, view.storage_identity);
  return hash_mix(hash, view.revision_domain);
}

std::uint64_t mix_stable_view_lineage(std::uint64_t hash,
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
  hash = hash_mix(hash, view.replica);
  hash = hash_mix(hash, view.field);
  hash = hash_mix(hash, view.storage_identity);
  return hash_mix(hash, view.revision_domain);
}

bool thermophysical_boundary_tokens(
    const BoundaryPlan& boundary, BoundaryThermophysicalGhostContext context,
    const BoundaryThermophysicalGhostUse& use, ConstFieldView density,
    PlanFingerprint expected_thermodynamics,
    PlanFingerprint expected_transport,
    ThermophysicalBoundaryTokens& tokens) noexcept {
  tokens = {};
  if (!context.valid() || context.numeric_boundary != boundary.revision() ||
      expected_thermodynamics == 0U || expected_transport == 0U)
    return false;
  if (entirely_periodic(boundary)) {
    std::uint64_t semantics = kFnvOffset;
    semantics = hash_mix(semantics, UINT64_C(0x706572696f646963));
    semantics = hash_mix(semantics, boundary.semantic_fingerprint());
    semantics = hash_mix(semantics, boundary.pressure_field());
    semantics = hash_mix(semantics, boundary.enthalpy_field());
    semantics = hash_mix(semantics, expected_thermodynamics);
    semantics = hash_mix(semantics, expected_transport);
    semantics = semantics == 0U ? 1U : semantics;
    std::uint64_t lineage = kFnvOffset;
    lineage = hash_mix(lineage, semantics);
    lineage = hash_mix(lineage, boundary.revision());
    lineage = hash_mix(lineage, context.target_time);
    lineage = hash_mix(lineage, context.geometry);
    lineage = hash_mix(lineage, context.pressure_reference != 0U);
    lineage = hash_mix(lineage, context.numeric_boundary);
    lineage = lineage == 0U ? 1U : lineage;
    std::uint64_t target = kFnvOffset;
    target = hash_mix(target, lineage);
    target = hash_mix(target, static_cast<std::uint8_t>(context.phase));
    target = target == 0U ? 1U : target;
    std::uint64_t local_lineage = kFnvOffset;
    local_lineage = hash_mix(local_lineage, lineage);
    local_lineage =
        hash_mix(local_lineage, boundary.local_layout_fingerprint());
    local_lineage = hash_mix(local_lineage, context.pressure_reference);
    local_lineage = mix_stable_view_lineage(local_lineage, density);
    local_lineage = local_lineage == 0U ? 1U : local_lineage;
    std::uint64_t local = kFnvOffset;
    local = hash_mix(local, local_lineage);
    local = hash_mix(local, target);
    tokens = {semantics, lineage, target, local_lineage,
              local == 0U ? 1U : local};
    return true;
  }
  const BoundaryThermophysicalGhostBinding& binding = use.binding;
  if (use.certificate.thermodynamics() != expected_thermodynamics ||
      use.certificate.transport() != expected_transport ||
      !same_field_identity(binding.density, density) ||
      !use.certificate.matches(boundary, context, binding)) {
    return false;
  }
  tokens = {use.certificate.collective_semantics(),
            use.certificate.collective_lineage(),
            use.certificate.collective_target(),
            use.certificate.rank_local_lineage(),
            use.certificate.rank_local_binding()};
  return tokens.valid();
}

std::uint64_t mix_candidate_field_values(std::uint64_t hash,
                                         ConstFieldView field, Int3 cells,
                                         std::uint8_t components,
                                         std::int32_t ghost_reach) noexcept {
  for (std::uint8_t component = 0U; component < components; ++component) {
    for (std::int32_t z = -ghost_reach; z < cells.z + ghost_reach; ++z) {
      for (std::int32_t y = -ghost_reach; y < cells.y + ghost_reach; ++y) {
        for (std::int32_t x = -ghost_reach; x < cells.x + ghost_reach; ++x) {
          const unsigned outside =
              static_cast<unsigned>(x < 0 || x >= cells.x) +
              static_cast<unsigned>(y < 0 || y >= cells.y) +
              static_cast<unsigned>(z < 0 || z >= cells.z);
          if (outside > 1U) continue;
          hash = hash_mix(
              hash, double_bits(field.unchecked({x, y, z}, component)));
        }
      }
    }
  }
  return hash;
}

void mix_candidate_field_values_pair(std::uint64_t& first,
                                     std::uint64_t& second,
                                     ConstFieldView field, Int3 cells,
                                     std::uint8_t components,
                                     std::int32_t ghost_reach) noexcept {
  for (std::uint8_t component = 0U; component < components; ++component) {
    for (std::int32_t z = -ghost_reach; z < cells.z + ghost_reach; ++z) {
      for (std::int32_t y = -ghost_reach; y < cells.y + ghost_reach; ++y) {
        for (std::int32_t x = -ghost_reach; x < cells.x + ghost_reach; ++x) {
          const unsigned outside =
              static_cast<unsigned>(x < 0 || x >= cells.x) +
              static_cast<unsigned>(y < 0 || y >= cells.y) +
              static_cast<unsigned>(z < 0 || z >= cells.z);
          if (outside > 1U) continue;
          const std::uint64_t bits =
              double_bits(field.unchecked({x, y, z}, component));
          first = hash_mix(first, bits);
          second = hash_mix(second, bits);
        }
      }
    }
  }
}

std::uint64_t mix_candidate_flux_values(std::uint64_t hash,
                                        ConstFaceFluxView flux) noexcept {
  const std::array<ConstFaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (ConstFaceFieldView face : faces) {
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          hash = hash_mix(hash,
                          double_bits(face.unchecked({x, y, z})));
  }
  return hash;
}

PlanFingerprint pressure_energy_refinement_state_numeric(
    PisoCoupledStateView state, Int3 cells) noexcept {
  std::uint64_t hash =
      hash_mix(kFnvOffset, UINT64_C(0x7265667374617465));
  hash = mix_candidate_field_values(
      hash, as_const(state.pressure_perturbation), cells, 1U, 0);
  hash = mix_candidate_field_values(hash, as_const(state.enthalpy), cells,
                                    1U, 0);
  hash = mix_candidate_field_values(hash, as_const(state.density), cells, 1U,
                                    0);
  hash = mix_candidate_field_values(hash, as_const(state.temperature), cells,
                                    1U, 0);
  hash = mix_candidate_field_values(hash, as_const(state.velocity), cells, 3U,
                                    0);
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

PlanFingerprint pressure_energy_refinement_flux_numeric(
    ConstFaceFluxView flux) noexcept {
  std::uint64_t hash =
      hash_mix(kFnvOffset, UINT64_C(0x726566666c7578));
  hash = mix_candidate_flux_values(hash, flux);
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

std::uint64_t mix_candidate_flux_binding(std::uint64_t hash,
                                         ConstFaceFluxView flux) noexcept {
  hash = mix_face(hash, flux.x);
  hash = mix_face(hash, flux.y);
  hash = mix_face(hash, flux.z);
  return hash_mix(hash, flux.revision);
}

std::uint64_t candidate_composition_numeric_local_provenance(
    Span<const ConstFieldView> independent_species, Int3 cells) noexcept {
  std::uint64_t hash =
      hash_mix(kFnvOffset, UINT64_C(0x636f6d706e756d76));
  hash = hash_mix(hash,
                  static_cast<std::uint64_t>(independent_species.size));
  for (std::size_t species = 0U; species < independent_species.size;
       ++species) {
    const ConstFieldView field = independent_species.data[species];
    hash = hash_mix(hash, field.components);
    for (std::uint8_t component = 0U; component < field.components;
         ++component)
      for (std::int32_t z = 0; z < cells.z; ++z)
        for (std::int32_t y = 0; y < cells.y; ++y)
          for (std::int32_t x = 0; x < cells.x; ++x)
            hash = hash_mix(
                hash,
                double_bits(field.unchecked({x, y, z}, component)));
  }
  return hash == 0U ? 1U : hash;
}

PlanFingerprint semantic_composition_rank_local_binding(
    Span<const ConstFieldView> independent_species) noexcept {
  std::uint64_t hash =
      hash_mix(kFnvOffset, UINT64_C(0x7465726d636f6d70));
  for (std::size_t species = 0U; species < independent_species.size;
       ++species)
    hash = mix_complete_view_identity(hash,
                                      independent_species.data[species]);
  return hash == 0U ? 1U : hash;
}

std::uint64_t exact_composition_identity_local(
    PlanFingerprint thermodynamics,
    Span<const ConstFieldView> independent_species, Int3 cells) noexcept {
  std::uint64_t hash =
      hash_mix(kFnvOffset, UINT64_C(0x706563636f6d706f));
  hash = hash_mix(hash, thermodynamics);
  hash = hash_mix(hash,
                  static_cast<std::uint64_t>(independent_species.size));
  for (std::size_t species = 0U; species < independent_species.size;
       ++species) {
    const ConstFieldView field = independent_species.data[species];
    hash = hash_mix(hash, field.components);
    for (std::uint8_t component = 0U; component < field.components;
         ++component)
      for (std::int32_t z = 0; z < cells.z; ++z)
        for (std::int32_t y = 0; y < cells.y; ++y)
          for (std::int32_t x = 0; x < cells.x; ++x)
            hash = hash_mix(
                hash,
                double_bits(field.unchecked({x, y, z}, component)));
  }
  return hash == 0U ? 1U : hash;
}

bool same_candidate_face_identity(ConstFaceFieldView left,
                                  ConstFaceFieldView right) noexcept {
  return left.base == right.base && same_cells(left.extents, right.extents) &&
         left.stride_y == right.stride_y &&
         left.stride_z == right.stride_z && left.axis == right.axis &&
         left.storage_identity == right.storage_identity &&
         left.revision_domain == right.revision_domain;
}

bool same_candidate_flux_identity(ConstFaceFluxView left,
                                  ConstFaceFluxView right) noexcept {
  return left.revision == right.revision &&
         same_candidate_face_identity(left.x, right.x) &&
         same_candidate_face_identity(left.y, right.y) &&
         same_candidate_face_identity(left.z, right.z);
}

bool valid_pressure_energy_refinement_state_views(
    PisoCoupledStateView state, Int3 cells) noexcept {
  const std::array<ConstFieldView, 5U> fields{
      as_const(state.velocity), as_const(state.pressure_perturbation),
      as_const(state.enthalpy), as_const(state.density),
      as_const(state.temperature)};
  if (!detail::valid_cell_view(fields[0U], cells, 0U, 3U, 0U)) return false;
  for (std::size_t index = 1U; index < fields.size(); ++index)
    if (!detail::valid_cell_view(fields[index], cells, 0U, 1U, 0U))
      return false;
  for (std::size_t index = 0U; index < fields.size(); ++index)
    for (std::size_t other = index + 1U; other < fields.size(); ++other)
      if (detail::field_views_overlap(fields[index], fields[other]))
        return false;
  return true;
}

bool refinement_state_flux_alias_free(PisoCoupledStateView state,
                                      ConstFaceFluxView flux) noexcept {
  const std::array<ConstFieldView, 5U> fields{
      as_const(state.velocity), as_const(state.pressure_perturbation),
      as_const(state.enthalpy), as_const(state.density),
      as_const(state.temperature)};
  const std::array<ConstFaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (ConstFieldView field : fields)
    for (ConstFaceFieldView face : faces)
      if (detail::cell_face_views_overlap(field, face)) return false;
  return true;
}

std::uint64_t frozen_exact_candidate_state_local_provenance(
    ConstFieldView pressure, ConstFieldView enthalpy, ConstFieldView density,
    ConstFieldView temperature, ConstFieldView velocity, Int3 cells,
    const PisoExactThermodynamicCandidateView& thermodynamic) noexcept {
  std::uint64_t hash =
      hash_mix(kFnvOffset, UINT64_C(0x66726f7a65786163));
  hash = hash_mix(hash, thermodynamic.closure.thermodynamics);
  hash = hash_mix(hash, thermodynamic.closure.pressure_reference);
  hash = hash_mix(hash, thermodynamic.closure.composition);
  hash = hash_mix(hash,
                  double_bits(thermodynamic.closed_gauge
                                  .next_pressure_reference));
  hash = mix_candidate_field_values(hash, pressure, cells, 1U, 0);
  hash = mix_candidate_field_values(hash, enthalpy, cells, 1U, 0);
  hash = mix_candidate_field_values(hash, density, cells, 1U, 0);
  hash = mix_candidate_field_values(hash, temperature, cells, 1U, 0);
  return mix_candidate_field_values(hash, velocity, cells, 3U, 0);
}

std::uint64_t frozen_exact_base_state_local_provenance(
    ConstFieldView pressure, ConstFieldView enthalpy, ConstFieldView density,
    ConstFieldView temperature, ConstFieldView velocity, Int3 cells,
    double absolute_pressure_reference) noexcept {
  std::uint64_t hash =
      hash_mix(kFnvOffset, UINT64_C(0x66726f7a65626173));
  hash = hash_mix(hash, double_bits(absolute_pressure_reference));
  hash = mix_candidate_field_values(hash, pressure, cells, 1U, 0);
  hash = mix_candidate_field_values(hash, enthalpy, cells, 1U, 0);
  hash = mix_candidate_field_values(hash, density, cells, 1U, 0);
  hash = mix_candidate_field_values(hash, temperature, cells, 1U, 0);
  return mix_candidate_field_values(hash, velocity, cells, 3U, 0);
}

std::uint64_t final_boundary_state_local_provenance(
    ConstFieldView pressure, ConstFieldView enthalpy, ConstFieldView density,
    ConstFieldView temperature, ConstFieldView velocity, Int3 cells,
    PlanFingerprint thermodynamics, RevisionToken pressure_reference,
    PlanFingerprint composition, double absolute_pressure_reference,
    PlanFingerprint composition_numeric_provenance) noexcept {
  std::uint64_t hash =
      hash_mix(kFnvOffset, UINT64_C(0x7065636266737461));
  hash = hash_mix(hash, thermodynamics);
  hash = hash_mix(hash, pressure_reference);
  hash = hash_mix(hash, composition);
  hash = hash_mix(hash, composition_numeric_provenance);
  hash = hash_mix(hash, double_bits(absolute_pressure_reference));
  const auto mix_with_coordinates =
      [&](std::uint64_t value, ConstFieldView field,
          std::uint8_t components) noexcept {
        for (std::int32_t z = 0; z < cells.z; ++z)
          for (std::int32_t y = 0; y < cells.y; ++y)
            for (std::int32_t x = 0; x < cells.x; ++x) {
              value = hash_mix(value, static_cast<std::uint32_t>(x));
              value = hash_mix(value, static_cast<std::uint32_t>(y));
              value = hash_mix(value, static_cast<std::uint32_t>(z));
              for (std::uint8_t component = 0U; component < components;
                   ++component)
                value = hash_mix(
                    value,
                    double_bits(field.unchecked({x, y, z}, component)));
            }
        return value;
      };
  hash = mix_with_coordinates(hash, pressure, 1U);
  hash = mix_with_coordinates(hash, enthalpy, 1U);
  hash = mix_with_coordinates(hash, density, 1U);
  hash = mix_with_coordinates(hash, temperature, 1U);
  return mix_with_coordinates(hash, velocity, 3U);
}

struct ExactCandidateFluxLocalAudit {
  std::uint64_t final_boundary_provenance{};
  std::uint64_t exact_candidate_provenance{};
  bool finite{};
};

ExactCandidateFluxLocalAudit exact_candidate_flux_local_audit(
    ConstFaceFluxView flux) noexcept {
  std::uint64_t final_boundary =
      hash_mix(kFnvOffset, UINT64_C(0x7065636266666c78));
  std::uint64_t exact_candidate =
      hash_mix(kFnvOffset, UINT64_C(0x6578616374666c78));
  bool finite = true;
  for (ConstFaceFieldView face :
       {flux.x, flux.y, flux.z}) {
    final_boundary =
        hash_mix(final_boundary, static_cast<std::uint8_t>(face.axis));
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x) {
          const double value = face.unchecked({x, y, z});
          const std::uint64_t bits = double_bits(value);
          final_boundary = hash_mix(final_boundary, bits);
          exact_candidate = hash_mix(exact_candidate, bits);
          finite = finite && std::isfinite(value);
        }
  }
  return {final_boundary, exact_candidate, finite};
}

bool same_piso_field_revision_identity(
    const PisoFieldRevisionIdentity& identity,
    ConstFieldView field) noexcept {
  return identity.valid() && identity.base == field.base &&
         identity.replica == field.replica && identity.field == field.field &&
         identity.revision == field.revision &&
         identity.storage == field.storage_identity &&
         identity.revision_domain == field.revision_domain;
}

bool empty_field_view(ConstFieldView field) noexcept {
  return field.base == nullptr && field.interior.x == 0 &&
         field.interior.y == 0 && field.interior.z == 0 &&
         field.ghosts.x == 0 && field.ghosts.y == 0 && field.ghosts.z == 0 &&
         field.components == 0U && field.stride_y == 0U &&
         field.stride_z == 0U && field.component_stride == 0U &&
         field.replica == 0U && field.field == 0U && field.revision == 0U &&
         field.storage_identity == 0U && field.revision_domain == 0U;
}

bool same_pressure_reference_certificate(
    const PressureReferenceCertificate& left,
    const PressureReferenceCertificate& right) noexcept {
  return left.plan == right.plan && left.predictor == right.predictor &&
         left.thermodynamics == right.thermodynamics &&
         left.closure == right.closure && left.time == right.time &&
         left.pressure_reference == right.pressure_reference &&
         left.kind == right.kind;
}

bool empty_closed_gauge_certificate(
    const ClosedGaugeCorrectionCertificate& certificate) noexcept {
  const PressureReferenceCertificate& output =
      certificate.output_pressure_reference;
  return certificate.shift == 0.0 &&
         certificate.next_pressure_reference == 0.0 &&
         certificate.local_moment == 0.0 &&
         certificate.local_weight == 0.0 &&
         certificate.global_moment == 0.0 &&
         certificate.global_weight == 0.0 &&
         certificate.local_post_shift_moment == 0.0 &&
         certificate.local_post_shift_absolute_moment == 0.0 &&
         certificate.global_post_shift_moment == 0.0 &&
         certificate.global_post_shift_absolute_moment == 0.0 &&
         certificate.post_shift_gauge_residual == 0.0 &&
         certificate.post_shift_gauge_tolerance == 0.0 &&
         output.plan == 0U && output.predictor == 0U &&
         output.thermodynamics == 0U && output.closure == 0U &&
         output.time == 0U && output.pressure_reference == 0U &&
         output.kind == PressureReferenceKind::boundary_absolute &&
         certificate.predecessor_pressure_reference == 0U &&
         certificate.time == 0U && certificate.geometry == 0U &&
         certificate.pressure_correction_authority == 0U &&
         certificate.target_thermodynamic_closure == 0U &&
         certificate.activity_local_fingerprint == 0U &&
         certificate.activity_collective_fingerprint == 0U &&
         certificate.collective_transaction == 0U &&
         certificate.rank_local_transaction == 0U &&
         certificate.local_active_cells == 0U &&
         certificate.corrector == 0U;
}

std::uint64_t mix_piso_field_revision_identity(
    std::uint64_t hash,
    const PisoFieldRevisionIdentity& identity) noexcept {
  hash = hash_mix(hash, reinterpret_cast<std::uintptr_t>(identity.base));
  hash = hash_mix(hash, identity.replica);
  hash = hash_mix(hash, identity.field);
  hash = hash_mix(hash, identity.revision);
  hash = hash_mix(hash, identity.storage);
  return hash_mix(hash, identity.revision_domain);
}

PlanFingerprint exact_eos_closure_fingerprint(
    const PisoExactEosClosureIdentity& identity) noexcept {
  std::uint64_t hash = kFnvOffset;
  hash = hash_mix(hash, UINT64_C(0x6578616374656f73));
  hash = hash_mix(hash, identity.thermodynamics);
  hash = hash_mix(hash, identity.pressure_reference);
  hash = hash_mix(hash, identity.composition);
  hash = mix_piso_field_revision_identity(hash, identity.pressure_state);
  hash = mix_piso_field_revision_identity(hash, identity.pressure_correction);
  hash = mix_piso_field_revision_identity(hash, identity.enthalpy_state);
  hash =
      mix_piso_field_revision_identity(hash, identity.enthalpy_correction);
  hash = mix_piso_field_revision_identity(hash, identity.candidate_enthalpy);
  hash = mix_piso_field_revision_identity(hash, identity.candidate_density);
  hash =
      mix_piso_field_revision_identity(hash, identity.candidate_temperature);
  hash = hash_mix(hash, identity.closure);
  return hash == 0U ? 1U : hash;
}

bool same_pressure_certificate(const PressureCorrectionCertificate& left,
                               const PressureCorrectionCertificate& right) noexcept {
  return left.plan == right.plan && left.intermediate == right.intermediate &&
         left.time == right.time && left.geometry == right.geometry &&
         left.numeric_boundary == right.numeric_boundary &&
         left.state == right.state && left.corrector == right.corrector &&
         left.thermophysical_boundary_semantics ==
             right.thermophysical_boundary_semantics &&
         left.thermophysical_boundary_target ==
             right.thermophysical_boundary_target &&
         left.thermophysical_boundary_rank_local_binding ==
             right.thermophysical_boundary_rank_local_binding &&
         left.thermophysical_boundary_collective_lineage ==
             right.thermophysical_boundary_collective_lineage &&
         left.thermophysical_boundary_rank_local_lineage ==
             right.thermophysical_boundary_rank_local_lineage &&
         left.pressure_energy_refinement ==
             right.pressure_energy_refinement &&
         left.pressure_energy_refinement_collective_lineage ==
             right.pressure_energy_refinement_collective_lineage &&
         left.pressure_energy_refinement_lineage ==
             right.pressure_energy_refinement_lineage;
}

bool same_state_correction_certificate(
    const PisoStateCorrectionCertificate& left,
    const PisoStateCorrectionCertificate& right) noexcept {
  return left.plan == right.plan &&
         left.pressure_system == right.pressure_system &&
         left.correction == right.correction &&
         left.enthalpy_correction == right.enthalpy_correction &&
         left.velocity == right.velocity && left.pressure == right.pressure &&
         left.enthalpy == right.enthalpy && left.density == right.density &&
         left.temperature == right.temperature &&
         left.face_flux == right.face_flux &&
         left.exact_eos_closure == right.exact_eos_closure &&
         left.state == right.state && left.corrector == right.corrector &&
         left.closure == right.closure &&
         left.thermophysical_boundary_semantics ==
             right.thermophysical_boundary_semantics &&
         left.thermophysical_boundary_target ==
             right.thermophysical_boundary_target &&
         left.thermophysical_boundary_rank_local_binding ==
             right.thermophysical_boundary_rank_local_binding &&
         left.thermophysical_boundary_collective_lineage ==
             right.thermophysical_boundary_collective_lineage &&
         left.thermophysical_boundary_rank_local_lineage ==
             right.thermophysical_boundary_rank_local_lineage &&
         same_pressure_reference_certificate(
             left.input_pressure_reference,
             right.input_pressure_reference) &&
         same_pressure_reference_certificate(
             left.output_pressure_reference,
             right.output_pressure_reference) &&
         left.closed_gauge_collective_transaction ==
             right.closed_gauge_collective_transaction &&
         left.closed_gauge_rank_local_transaction ==
             right.closed_gauge_rank_local_transaction &&
         left.composition_numeric_provenance ==
             right.composition_numeric_provenance &&
         left.composition_rank_local_binding ==
             right.composition_rank_local_binding &&
         left.independent_species_count ==
             right.independent_species_count;
}

bool same_linear_identity(const LinearIdentity& left,
                          const LinearIdentity& right) noexcept {
  return left.symbolic == right.symbolic && left.numeric == right.numeric &&
         left.hierarchy == right.hierarchy && left.workspace == right.workspace &&
         left.fingerprint == right.fingerprint;
}

bool empty_pressure_energy_refinement_report(
    const PisoPressureEnergyRefinementSolveReport& report) noexcept {
  return report.target_generation == 0U && report.collective_lineage == 0U &&
         report.pressure_state == 0U && report.linear_identity.symbolic == 0U &&
         report.linear_identity.numeric == 0U &&
         report.linear_identity.hierarchy == 0U &&
         report.linear_identity.workspace == 0U &&
         report.linear_identity.fingerprint == 0U && report.ordinal == 0U &&
         report.solve.status.code == StatusCode::ok &&
         report.solve.status.detail == 0U &&
         report.solve.termination == LinearTermination::invalid_plan;
}

bool valid_pressure_energy_refinement_prefix(
    const std::array<PisoPressureEnergyRefinementSolveReport,
                     kPressureEnergyRefinementCapacity>& reports,
    std::uint8_t count) noexcept {
  if (count > kPressureEnergyRefinementCapacity) return false;
  RevisionToken target_generation = 0U;
  for (std::size_t index = 0U; index < reports.size(); ++index) {
    if (index >= count) {
      if (!empty_pressure_energy_refinement_report(reports[index]))
        return false;
      continue;
    }
    const PisoPressureEnergyRefinementSolveReport& report = reports[index];
    if (!report.valid() || report.ordinal != index + 1U) {
      return false;
    }
    if (index == 0U)
      target_generation = report.target_generation;
    else if (report.target_generation != target_generation)
      return false;
    for (std::size_t prior = 0U; prior < index; ++prior)
      if (reports[prior].collective_lineage == report.collective_lineage)
        return false;
  }
  return true;
}

bool valid_pressure_energy_refinement_termination(
    PressureEnergyRefinementTermination termination,
    std::uint8_t count) noexcept {
  switch (termination) {
    case PressureEnergyRefinementTermination::none:
      return count == 0U;
    case PressureEnergyRefinementTermination::component_residuals_converged:
      return count <= kPressureEnergyRefinementCapacity;
    case PressureEnergyRefinementTermination::iteration_capacity_exhausted:
      return count == kPressureEnergyRefinementCapacity;
    case PressureEnergyRefinementTermination::rejected_candidate:
      return false;
  }
  return false;
}

bool same_linear_operator_certificate(
    const LinearOperatorCertificate& left,
    const LinearOperatorCertificate& right) noexcept {
  return same_linear_identity(left.identity, right.identity) &&
         left.collective_fingerprint == right.collective_fingerprint &&
         same_cells(left.local_shape, right.local_shape) &&
         left.operator_class == right.operator_class;
}

double face_pressure_coefficient(
    const CartesianKernelPlan& kernels, ConstFieldView density,
    ConstFieldView r_au, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch,
    const PressureCorrectionBoundaryPlan& pressure_boundary,
    CartesianAxis axis, Int3 face) noexcept {
  const Int3 cells = kernels.cells();
  const std::int32_t normal = axis == CartesianAxis::x
                                  ? face.x
                                  : (axis == CartesianAxis::y ? face.y
                                                              : face.z);
  const std::int32_t extent = axis == CartesianAxis::x
                                  ? cells.x
                                  : (axis == CartesianAxis::y ? cells.y
                                                              : cells.z);
  const std::uint8_t component =
      axis == CartesianAxis::x ? 0U : (axis == CartesianAxis::y ? 1U : 2U);
  (void)patch;
  PressureCorrectionFaceRule rule;
  if (!pressure_boundary.face_rule(axis, face, rule)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const bool physical = rule.physical;
  const bool periodic =
      rule.kind == PressureCorrectionFaceKind::periodic;
  if (rule.kind == PressureCorrectionFaceKind::homogeneous_neumann) {
    return 0.0;
  }
  Int3 left = face;
  if (axis == CartesianAxis::x) {
    --left.x;
  } else if (axis == CartesianAxis::y) {
    --left.y;
  } else {
    --left.z;
  }
  if (physical && !periodic) {
    const Int3 owner = normal == 0 ? face : left;
    const Int3 ghost = normal == 0 ? left : face;
    const double owner_coefficient =
        density.unchecked(owner, 0U) * r_au.unchecked(owner, component);
    const double ghost_coefficient =
        density.unchecked(ghost, 0U) * r_au.unchecked(ghost, component);
    const double distance =
        normal == 0
            ? detail::centre_coordinate(kernels, axis, 0) -
                  detail::face_coordinate(kernels, axis, 0)
            : detail::face_coordinate(kernels, axis, extent) -
                  detail::centre_coordinate(kernels, axis, extent - 1);
    // p'=0 is imposed at the boundary face, so the correction gradient spans
    // one owner-to-face distance.  The thermodynamic mobility at that face is
    // the harmonic interpolation of the owner and EOS-closed ghost states.
    return 2.0 * detail::face_area(kernels, axis, face) /
           (distance / owner_coefficient + distance / ghost_coefficient);
  }
  const double left_coefficient =
      density.unchecked(left, 0U) * r_au.unchecked(left, component);
  const double right_coefficient =
      density.unchecked(face, 0U) * r_au.unchecked(face, component);
  double left_distance = 0.0;
  double right_distance = 0.0;
  if (physical && periodic) {
    const Span<const double> widths = geometry.axis(axis).widths();
    left_distance = 0.5 * widths.data[widths.size - 1U];
    right_distance = 0.5 * widths.data[0U];
  } else {
    const double location = detail::face_coordinate(kernels, axis, normal);
    left_distance =
        location - detail::centre_coordinate(kernels, axis, normal - 1);
    right_distance =
        detail::centre_coordinate(kernels, axis, normal) - location;
  }
  return detail::face_area(kernels, axis, face) /
         (left_distance / left_coefficient +
          right_distance / right_coefficient);
}

template <CartesianAxis Axis>
Status fill_pressure_coefficients(const CartesianKernelPlan& kernels,
                                  ConstFieldView density,
                                  ConstFieldView r_au,
                                  const CartesianGeometryPlan& geometry,
                                  const MeshPatch& patch,
                                  const PressureCorrectionBoundaryPlan&
                                      pressure_boundary,
                                  FaceFieldView output) noexcept {
  const Int3 extents = output.extents;
  for (std::int32_t z = 0; z < extents.z; ++z) {
    for (std::int32_t y = 0; y < extents.y; ++y) {
      for (std::int32_t x = 0; x < extents.x; ++x) {
        const Int3 face{x, y, z};
        const double coefficient = face_pressure_coefficient(
            kernels, density, r_au, geometry, patch, pressure_boundary, Axis,
            face);
        if (!std::isfinite(coefficient) || coefficient < 0.0) {
          return {StatusCode::numerical_failure, kPisoNumerical};
        }
        output.unchecked(face) = coefficient;
      }
    }
  }
  return {};
}

void fill_physical_zero_gradient(FieldView field, Int3 cells,
                                 std::uint8_t components,
                                 std::uint8_t reach,
                                 const BoundaryPlan& boundary) noexcept {
  for (std::uint8_t layer = 1U; layer <= reach; ++layer) {
    const std::int32_t offset = static_cast<std::int32_t>(layer);
    for (std::size_t face_index = 0U; face_index < 6U; ++face_index) {
      const auto selected = static_cast<CartesianFace>(face_index);
      const BoundaryFacePlan* plan = nullptr;
      if (!boundary.face(selected, plan) || plan == nullptr ||
          !plan->local_owner || plan->periodic) {
        continue;
      }
      const bool high = (face_index & 1U) != 0U;
      const CartesianAxis axis =
          face_index < 2U ? CartesianAxis::x
                          : (face_index < 4U ? CartesianAxis::y
                                             : CartesianAxis::z);
      const std::int32_t inner_count =
          axis == CartesianAxis::x ? cells.y : cells.x;
      const std::int32_t outer_count =
          axis == CartesianAxis::z ? cells.y : cells.z;
      for (std::int32_t outer = 0; outer < outer_count; ++outer) {
        for (std::int32_t inner = 0; inner < inner_count; ++inner) {
          Int3 source{};
          Int3 destination{};
          if (axis == CartesianAxis::x) {
            source = {high ? cells.x - 1 : 0, inner, outer};
            destination = {high ? cells.x - 1 + offset : -offset, inner,
                           outer};
          } else if (axis == CartesianAxis::y) {
            source = {inner, high ? cells.y - 1 : 0, outer};
            destination = {inner, high ? cells.y - 1 + offset : -offset,
                           outer};
          } else {
            source = {inner, outer, high ? cells.z - 1 : 0};
            destination = {inner, outer,
                           high ? cells.z - 1 + offset : -offset};
          }
          for (std::uint8_t component = 0U; component < components;
               ++component) {
            field.unchecked(destination, component) =
                field.unchecked(source, component);
          }
        }
      }
    }
  }
}

}  // namespace

bool PisoCartesianPressureWorkLinearization::valid() const noexcept {
  return issuer_ != nullptr && kernels_ != nullptr && geometry_ != nullptr &&
         patch_.cells.x > 0 && patch_.cells.y > 0 && patch_.cells.z > 0 &&
         intermediate_.valid() && pressure_.valid() &&
         target_pressure_perturbation_.base != nullptr &&
         target_velocity_.base != nullptr && h_by_a_.base != nullptr &&
         r_au_.base != nullptr && authority_ != 0U &&
         numeric_fingerprint_ != 0U;
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace detail {

void arm_piso_post_halo_revalidation_failure_once_for_test(
    int failing_rank) noexcept {
  g_post_halo_revalidation_failure_rank.store(failing_rank,
                                               std::memory_order_relaxed);
}

void clear_piso_post_halo_revalidation_failure_for_test() noexcept {
  g_post_halo_revalidation_failure_rank.store(-1,
                                               std::memory_order_relaxed);
}

}  // namespace detail

FinalBoundaryFluxCertificate
FinalBoundaryFluxCertificateTestAccess::make_foreign(
    const PressureVelocityCoupler* claimed_coupler,
    ConstFaceFluxView final_flux) noexcept {
  FinalBoundaryFluxCertificate candidate;
  candidate.issuer_ =
      reinterpret_cast<const PressureEnergyCandidateBoundaryFinalizer*>(
          std::uintptr_t{1U});
  candidate.coupler_ = claimed_coupler;
  candidate.scope_ =
      PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm;
  candidate.stage_lineage_ = 1U;
  candidate.velocity_lineage_ = 1U;
  candidate.mechanical_flux_lineage_ = 1U;
  candidate.nonphysical_flux_provenance_ = 1U;
  candidate.pressure_outlet_provisional_provenance_ = 1U;
  candidate.boundary_semantic_ = 1U;
  candidate.boundary_layout_ = 1U;
  candidate.pressure_reference_ =
      {1U, 1U, 1U, 1U, 1U, 1U,
       PressureReferenceKind::boundary_absolute};
  candidate.absolute_pressure_reference_ = 101325.0;
  candidate.target_time_ = 1U;
  candidate.alpha_ = 0.5;
  const PisoFieldRevisionIdentity field_identity{
      final_flux.x.base, 0U, 1U, 1U, 1U, 1U};
  candidate.candidate_pressure_ = field_identity;
  candidate.candidate_enthalpy_ = field_identity;
  candidate.candidate_density_ = field_identity;
  candidate.candidate_temperature_ = field_identity;
  candidate.candidate_velocity_ = field_identity;
  candidate.composition_identity_ = 1U;
  candidate.composition_numeric_provenance_ = 1U;
  candidate.independent_species_views_ = {};
  candidate.independent_species_count_ = 0U;
  candidate.candidate_state_provenance_ = 1U;
  candidate.candidate_state_binding_ = 1U;
  candidate.state_halo_instance_ = 1U;
  candidate.state_halo_lineage_ = 1U;
  candidate.thermophysical_boundary_semantics_ = 1U;
  candidate.thermophysical_boundary_target_ = 1U;
  candidate.thermophysical_boundary_collective_lineage_ = 1U;
  candidate.thermophysical_boundary_rank_local_lineage_ = 1U;
  candidate.thermophysical_boundary_rank_local_binding_ = 1U;
  candidate.physical_state_halo_lineage_ = 1U;
  candidate.face_closure_lineage_ = 1U;
  candidate.inlet_flux_provenance_ = 1U;
  candidate.outlet_backflow_provenance_ = 1U;
  candidate.ibm_donor_lineage_ = 1U;
  candidate.ibm_geometry_lineage_ = 1U;
  candidate.ibm_zero_interface_lineage_ = 1U;
  candidate.final_flux_view_ = final_flux;
  candidate.final_flux_revision_ = final_flux.revision;
  candidate.final_flux_storage_ = final_flux.x.storage_identity;
  candidate.final_flux_revision_domain_ = final_flux.x.revision_domain;
  candidate.final_flux_provenance_ = 1U;
  candidate.canonical_lineage_ = 1U;
  candidate.scratch_binding_ = 1U;
  candidate.corrector_ = 1U;
  return candidate;
}
#endif

Status PressureCorrectionBoundaryPlan::compile(
    const CartesianGeometryPlan& geometry, MeshPatch patch,
    const BoundaryPlan& boundary,
    PressureCorrectionBoundaryPlan& out) noexcept {
  const Int3 global_cells = geometry.global_cells();
  int patch_rank_value = 0;
  int patch_rank_count = 0;
  const bool valid_patch =
      patch_rank(patch, patch_rank_value, patch_rank_count) &&
      global_cells.x > 0 && global_cells.y > 0 && global_cells.z > 0 &&
      patch.cells.x > 0 && patch.cells.y > 0 && patch.cells.z > 0 &&
      patch.begin.x >= 0 && patch.begin.y >= 0 && patch.begin.z >= 0 &&
      patch.begin.x <= global_cells.x - patch.cells.x &&
      patch.begin.y <= global_cells.y - patch.cells.y &&
      patch.begin.z <= global_cells.z - patch.cells.z &&
      same_cells(boundary.local_cells(), patch.cells);
  if (geometry.fingerprint() == 0U || geometry.topology_revision() == 0U ||
      boundary.semantic_fingerprint() == 0U ||
      boundary.local_layout_fingerprint() == 0U ||
      boundary.revision() == 0U ||
      boundary.geometry_fingerprint() != geometry.fingerprint() ||
      !valid_patch) {
    return {StatusCode::invalid_plan, kPisoPressureBoundary};
  }

  PressureCorrectionBoundaryPlan candidate;
  candidate.geometry_ = &geometry;
  candidate.boundary_ = &boundary;
  candidate.patch_ = patch;
  candidate.geometry_fingerprint_ = geometry.fingerprint();
  candidate.boundary_semantic_ = boundary.semantic_fingerprint();
  candidate.boundary_local_layout_ = boundary.local_layout_fingerprint();

  std::uint64_t semantic = hash_mix(
      kFnvOffset, kPressureCorrectionBoundarySemanticSchema);
  semantic = hash_mix(semantic, geometry.fingerprint());
  semantic = hash_mix(semantic, boundary.semantic_fingerprint());
  std::uint64_t local =
      hash_mix(kFnvOffset, kPressureCorrectionBoundaryLayoutSchema);
  bool valid_faces = true;
  for (std::size_t index = 0U; index < candidate.kinds_.size(); ++index) {
    const auto face = static_cast<CartesianFace>(index);
    const BoundaryFacePlan* source = nullptr;
    if (!boundary.face(face, source) || source == nullptr) {
      valid_faces = false;
      break;
    }
    const PressureCorrectionFaceKind kind =
        pressure_correction_kind(source->flow_kind);
    const bool periodic = kind == PressureCorrectionFaceKind::periodic;
    const CartesianAxis axis =
        index < 2U ? CartesianAxis::x
                   : (index < 4U ? CartesianAxis::y : CartesianAxis::z);
    const bool high = (index & 1U) != 0U;
    const std::int32_t begin = axis_value(patch.begin, axis);
    const std::int32_t cells = axis_value(patch.cells, axis);
    const std::int32_t global = axis_value(global_cells, axis);
    const bool expected_owner = high ? begin + cells == global : begin == 0;
    if (source->periodic != periodic ||
        source->local_owner != expected_owner) {
      valid_faces = false;
      break;
    }
    candidate.kinds_[index] = kind;
    candidate.local_physical_[index] = expected_owner;
    semantic = hash_mix(semantic, static_cast<std::uint8_t>(kind));
  }
  if (!valid_faces) {
    return {StatusCode::invalid_plan, kPisoPressureBoundary};
  }
  if (semantic == 0U) semantic = 1U;
  local = hash_mix(local, semantic);
  local = hash_mix(local, boundary.local_layout_fingerprint());
  const Int3 local_values[]{patch.begin, patch.cells, patch.process_grid,
                            patch.process_coord};
  for (Int3 value : local_values) {
    local = hash_mix(local, static_cast<std::uint32_t>(value.x));
    local = hash_mix(local, static_cast<std::uint32_t>(value.y));
    local = hash_mix(local, static_cast<std::uint32_t>(value.z));
  }
  for (bool physical : candidate.local_physical_) {
    local = hash_mix(local, physical ? 1U : 0U);
  }
  if (local == 0U) local = 1U;
  candidate.certificate_ = {semantic,
                            local,
                            geometry.topology_revision(),
                            boundary.revision(),
                            global_cells,
                            patch.begin,
                            patch.cells,
                            geometry.fingerprint()};
  if (!candidate.certificate_.valid()) {
    return {StatusCode::invalid_plan, kPisoPressureBoundary};
  }
  out = std::move(candidate);
  return {};
}

bool PressureCorrectionBoundaryPlan::current() const noexcept {
  return certificate_.valid() && geometry_ != nullptr && boundary_ != nullptr &&
         geometry_->fingerprint() == geometry_fingerprint_ &&
         geometry_->fingerprint() == certificate_.geometry_fingerprint &&
         geometry_->topology_revision() == certificate_.geometry &&
         same_cells(geometry_->global_cells(), certificate_.global_cells) &&
         boundary_->semantic_fingerprint() == boundary_semantic_ &&
         boundary_->geometry_fingerprint() ==
             certificate_.geometry_fingerprint &&
         boundary_->local_layout_fingerprint() == boundary_local_layout_ &&
         boundary_->revision() == certificate_.source_revision &&
         same_cells(boundary_->local_cells(), certificate_.local_cells);
}

PressureCorrectionFaceKind PressureCorrectionBoundaryPlan::kind(
    CartesianFace face) const noexcept {
  const std::size_t index = static_cast<std::size_t>(face);
  return current() && index < kinds_.size()
             ? kinds_[index]
             : PressureCorrectionFaceKind::invalid;
}

Status PressureCorrectionBoundaryPlan::face_rule(
    CartesianAxis axis, Int3 face, PressureCorrectionFaceRule& out) const
    noexcept {
  if (!current()) {
    return {StatusCode::invalid_plan, kPisoPressureBoundary};
  }
  return face_rule_unchecked(axis, face, out);
}

Status PressureCorrectionBoundaryPlan::face_rule_unchecked(
    CartesianAxis axis, Int3 face, PressureCorrectionFaceRule& out) const
    noexcept {
  const std::size_t axis_index = static_cast<std::size_t>(axis);
  if (axis_index >= 3U) {
    return {StatusCode::invalid_plan, kPisoPressureBoundary};
  }
  const Int3 cells = certificate_.local_cells;
  const std::int32_t normal = axis_value(face, axis);
  const std::int32_t extent = axis_value(cells, axis);
  const bool tangential_valid =
      axis == CartesianAxis::x
          ? face.y >= 0 && face.y < cells.y && face.z >= 0 &&
                face.z < cells.z
          : (axis == CartesianAxis::y
                 ? face.x >= 0 && face.x < cells.x && face.z >= 0 &&
                       face.z < cells.z
                 : face.x >= 0 && face.x < cells.x && face.y >= 0 &&
                       face.y < cells.y);
  if (!tangential_valid || normal < 0 || normal > extent) {
    return {StatusCode::invalid_plan, kPisoPressureBoundary};
  }
  const bool local_edge = normal == 0 || normal == extent;
  if (!local_edge) {
    out = {PressureCorrectionFaceKind::exchanged, false, false};
    return {};
  }
  const bool high = normal == extent;
  const std::int32_t global_face =
      axis_value(certificate_.patch_begin, axis) + normal;
  const bool physical =
      global_face == 0 ||
      global_face == axis_value(certificate_.global_cells, axis);
  if (!physical) {
    out = {PressureCorrectionFaceKind::exchanged, false, high};
    return {};
  }
  const std::size_t index = 2U * axis_index + (high ? 1U : 0U);
  if (!local_physical_[index]) {
    return {StatusCode::invalid_plan, kPisoPressureBoundary};
  }
  out = {kinds_[index], true, high};
  return {};
}

Status PressureCorrectionBoundaryPlan::mg_boundaries(
    MgBoundarySet& out) const noexcept {
  if (!current()) {
    return {StatusCode::invalid_plan, kPisoPressureBoundary};
  }
  const MgBoundarySet candidate{
      mg_boundary_kind(kinds_[0U]), mg_boundary_kind(kinds_[1U]),
      mg_boundary_kind(kinds_[2U]), mg_boundary_kind(kinds_[3U]),
      mg_boundary_kind(kinds_[4U]), mg_boundary_kind(kinds_[5U])};
  out = candidate;
  return {};
}

Status PressureCorrectionBoundaryPlan::fill_ghosts(
    FieldView correction) const noexcept {
  const Int3 cells = certificate_.local_cells;
  if (!current() ||
      !detail::valid_cell_view(as_const(correction), cells, 0U, 1U, 1U) ||
      correction.components != 1U) {
    return {StatusCode::invalid_plan, kPisoPressureBoundary};
  }
  for (std::size_t index = 0U; index < kinds_.size(); ++index) {
    const PressureCorrectionFaceKind face_kind = kinds_[index];
    if (!local_physical_[index] ||
        face_kind == PressureCorrectionFaceKind::periodic) {
      continue;
    }
    const CartesianAxis axis =
        index < 2U ? CartesianAxis::x
                   : (index < 4U ? CartesianAxis::y : CartesianAxis::z);
    const bool high = (index & 1U) != 0U;
    const std::int32_t inner_count =
        axis == CartesianAxis::x ? cells.y : cells.x;
    const std::int32_t outer_count =
        axis == CartesianAxis::z ? cells.y : cells.z;
    for (std::int32_t outer = 0; outer < outer_count; ++outer) {
      for (std::int32_t inner = 0; inner < inner_count; ++inner) {
        Int3 owner{};
        Int3 ghost{};
        if (axis == CartesianAxis::x) {
          owner = {high ? cells.x - 1 : 0, inner, outer};
          ghost = {high ? cells.x : -1, inner, outer};
        } else if (axis == CartesianAxis::y) {
          owner = {inner, high ? cells.y - 1 : 0, outer};
          ghost = {inner, high ? cells.y : -1, outer};
        } else {
          owner = {inner, outer, high ? cells.z - 1 : 0};
          ghost = {inner, outer, high ? cells.z : -1};
        }
        const double value = correction.unchecked(owner, 0U);
        correction.unchecked(ghost, 0U) =
            face_kind == PressureCorrectionFaceKind::homogeneous_dirichlet
                ? -value
                : value;
      }
    }
  }
  return {};
}

double PressureCorrectionBoundaryPlan::neighbor_value(
    ConstFieldView correction, Int3 cell, CartesianAxis axis,
    int direction) const noexcept {
  const Int3 cells = certificate_.local_cells;
  if (!current() || (direction != -1 && direction != 1) ||
      !detail::valid_cell_view(correction, cells, 0U, 1U, 1U) || cell.x < 0 ||
      cell.y < 0 || cell.z < 0 || cell.x >= cells.x || cell.y >= cells.y ||
      cell.z >= cells.z) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return neighbor_value_unchecked(correction, cell, axis, direction);
}

double PressureCorrectionBoundaryPlan::neighbor_value_unchecked(
    ConstFieldView correction, Int3 cell, CartesianAxis axis,
    int direction) const noexcept {
  const Int3 cells = certificate_.local_cells;
  Int3 selected = cell;
  if (axis == CartesianAxis::x) {
    selected.x += direction;
  } else if (axis == CartesianAxis::y) {
    selected.y += direction;
  } else if (axis == CartesianAxis::z) {
    selected.z += direction;
  } else {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const std::int32_t coordinate = axis_value(selected, axis);
  const std::int32_t extent = axis_value(cells, axis);
  if (coordinate >= 0 && coordinate < extent) {
    return correction.unchecked(selected, 0U);
  }
  Int3 face = cell;
  if (direction > 0) {
    if (axis == CartesianAxis::x) {
      ++face.x;
    } else if (axis == CartesianAxis::y) {
      ++face.y;
    } else {
      ++face.z;
    }
  }
  PressureCorrectionFaceRule rule;
  if (!face_rule_unchecked(axis, face, rule)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (rule.kind == PressureCorrectionFaceKind::homogeneous_dirichlet) {
    return 0.0;
  }
  if (rule.kind == PressureCorrectionFaceKind::homogeneous_neumann) {
    return correction.unchecked(cell, 0U);
  }
  return correction.unchecked(selected, 0U);
}

double PressureCorrectionBoundaryPlan::jump(ConstFieldView correction,
                                             CartesianAxis axis,
                                             Int3 face) const noexcept {
  if (!current() || !detail::valid_cell_view(
                        correction, certificate_.local_cells, 0U, 1U, 1U)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return jump_unchecked(correction, axis, face);
}

double PressureCorrectionBoundaryPlan::jump_unchecked(
    ConstFieldView correction, CartesianAxis axis, Int3 face) const noexcept {
  PressureCorrectionFaceRule rule;
  if (!face_rule_unchecked(axis, face, rule)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  Int3 left = face;
  if (axis == CartesianAxis::x) {
    --left.x;
  } else if (axis == CartesianAxis::y) {
    --left.y;
  } else {
    --left.z;
  }
  if (rule.kind == PressureCorrectionFaceKind::homogeneous_neumann) {
    return 0.0;
  }
  if (rule.kind == PressureCorrectionFaceKind::homogeneous_dirichlet) {
    const Int3 owner = rule.high ? left : face;
    const double value = correction.unchecked(owner, 0U);
    return rule.high ? -value : value;
  }
  return correction.unchecked(face, 0U) - correction.unchecked(left, 0U);
}

double PressureCorrectionBoundaryPlan::mass_flux_response(
    ConstFieldView correction, CartesianAxis axis, Int3 face,
    double pressure_face_coefficient, double correction_scale) const noexcept {
  if (!current() || !detail::valid_cell_view(
                        correction, certificate_.local_cells, 0U, 1U, 1U)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return mass_flux_response_unchecked(correction, axis, face,
                                      pressure_face_coefficient,
                                      correction_scale);
}

double PressureCorrectionBoundaryPlan::mass_flux_response_unchecked(
    ConstFieldView correction, CartesianAxis axis, Int3 face,
    double pressure_face_coefficient, double correction_scale) const noexcept {
  if (!std::isfinite(pressure_face_coefficient) ||
      !std::isfinite(correction_scale)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double correction_jump = jump_unchecked(correction, axis, face);
  const double applied_jump = correction_scale == 1.0
                                  ? correction_jump
                                  : correction_scale * correction_jump;
  return -pressure_face_coefficient * applied_jump;
}

Status PisoPlan::compile(MPI_Comm communicator,
                         const EquationPlanSet& equations,
                         const PisoPlanSpec& spec, PisoPlan& out,
                         PisoCompileDiagnostics* diagnostics) noexcept {
  if (diagnostics != nullptr) {
    diagnostics->lowest_failing_rank = -1;
  }
  int rank = 0;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kPisoCollective};
  }
  const Int3 cells = equations.pressure_reference().cells();
  Status local;
  if (equations.fingerprint() == 0U ||
      equations.thermophysical_predictor().fingerprint() == 0U ||
      equations.pressure_reference().fingerprint() == 0U || cells.x <= 0 ||
      cells.y <= 0 || cells.z <= 0 ||
      (spec.coupling != CouplingKind::piso &&
       spec.coupling != CouplingKind::simple) ||
      spec.pressure_correctors != 2U ||
      spec.pressure_stage == 0U ||
      !valid_pressure_algorithm(spec.pressure_algorithm) ||
      !valid_control(spec.pressure_algorithm, spec.mg_correction_scaling,
                     spec.pressure_solve) ||
      !finite_positive(spec.eos_tolerance) ||
      !finite_positive(spec.continuity_tolerance) ||
      (!finite_positive(spec.energy_tolerance) &&
       spec.energy_tolerance != 0.0) ||
      !finite_positive(spec.closed_mass_tolerance) ||
      !finite_positive(spec.gauge_tolerance) ||
      spec.pressure_solve.relative_tolerance >
          std::min(spec.continuity_tolerance, spec.gauge_tolerance)) {
    local = {StatusCode::invalid_plan, kPisoPlan};
  }
  int lowest = -1;
  Status consensus = collective_status(communicator, local, rank, size, lowest);
  if (!consensus) {
    if (diagnostics != nullptr) {
      diagnostics->lowest_failing_rank = lowest;
    }
    return consensus;
  }
  const PlanFingerprint fingerprint = spec_fingerprint(equations, spec);
  PlanFingerprint minimum = fingerprint;
  PlanFingerprint maximum = fingerprint;
  const int minimum_rc = MPI_Allreduce(MPI_IN_PLACE, &minimum, 1,
                                       MPI_UINT64_T, MPI_MIN, communicator);
  const int maximum_rc = MPI_Allreduce(MPI_IN_PLACE, &maximum, 1,
                                       MPI_UINT64_T, MPI_MAX, communicator);
  if (minimum_rc != MPI_SUCCESS || maximum_rc != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPisoCollective};
  }
  if (minimum != maximum) {
    if (diagnostics != nullptr) {
      diagnostics->lowest_failing_rank = 0;
    }
    return {StatusCode::invalid_plan, kPisoCollective};
  }
  PisoPlan candidate;
  candidate.cells_ = cells;
  candidate.equations_fingerprint_ = equations.fingerprint();
  candidate.predictor_fingerprint_ =
      equations.thermophysical_predictor().fingerprint();
  candidate.pressure_reference_fingerprint_ =
      equations.pressure_reference().fingerprint();
  candidate.coupling_ = spec.coupling;
  candidate.pressure_stage_ = spec.pressure_stage;
  candidate.final_flux_slot_ = spec.final_flux_slot;
  candidate.pressure_algorithm_ = spec.pressure_algorithm;
  candidate.mg_correction_scaling_ = spec.mg_correction_scaling;
  candidate.pressure_solve_ = spec.pressure_solve;
  candidate.eos_tolerance_ = spec.eos_tolerance;
  candidate.continuity_tolerance_ = spec.continuity_tolerance;
  candidate.energy_tolerance_ = spec.energy_tolerance;
  candidate.closed_mass_tolerance_ = spec.closed_mass_tolerance;
  candidate.gauge_tolerance_ = spec.gauge_tolerance;
  candidate.fingerprint_ = fingerprint;
  out = std::move(candidate);
  return {};
}

struct PressureVelocityCoupler::Impl {
  struct SealedPressureCorrectionAuthority {
    bool valid{};
    PressureCorrectionCertificate pressure{};
    LinearOperatorCertificate exact_operator{};
    ConstFieldView correction{};
    double alpha{};
  };

  struct PredecessorPressureApplication {
    bool valid{};
    double alpha{1.0};
    double maximum_depletion{};
  };

  struct CorrectedPendingAuthority {
    PisoStateCorrectionCertificate correction{};
    const PendingFaceFluxView* view{};
    FinalFaceFluxWriter* writer{};
    FaceFluxStorage* storage{};
    RevisionToken revision{};
    std::uint64_t writer_identity{};
    std::uint64_t attempt_identity{};
    PisoFieldRevisionIdentity pressure_compressibility{};
    Span<const std::uint8_t> active_cells{};
    PlanFingerprint activity_local_fingerprint{};
    PlanFingerprint activity_collective_fingerprint{};
  };

  const CartesianKernelPlan* kernels{};
  const CartesianGeometryPlan* geometry{};
  const BoundaryPlan* boundary{};
  const ThermodynamicsPlan* thermodynamics_plan{};
  const PressureReferencePlan* pressure_reference{};
  MPI_Comm communicator{MPI_COMM_NULL};
  int rank{-1};
  int size{};
  PressureCorrectionBoundaryPlan pressure_boundary{};
  HaloEngine* halo{};
  HaloEngine* correction_halo{};
  MeshPatch patch{};
  StageId halo_stage{};
  StageId correction_halo_stage{};
  FieldId density_field{};
  FieldId correction_field{};
  PisoCouplerWorkspace workspace{};
  Int3 cells{};
  PlanFingerprint plan{};
  PlanFingerprint momentum_plan{};
  PlanFingerprint predictor_plan{};
  PlanFingerprint pressure_reference_plan{};
  PlanFingerprint thermodynamics{};
  PlanFingerprint transport{};
  PressureReferenceKind pressure_reference_kind{
      PressureReferenceKind::boundary_absolute};
  PressureGauge pressure_gauge{PressureGauge::absolute_boundary_dirichlet};
  PlanFingerprint equations{};
  PlanFingerprint fingerprint{};
  LinearAlgorithm pressure_algorithm{LinearAlgorithm::fgmres};
  MgCorrectionScaling mg_correction_scaling{
      MgCorrectionScaling::residual_minimizing};
  double eos_tolerance{};
  double continuity_tolerance{};
  double energy_tolerance{};
  double closed_mass_tolerance{};
  double gauge_tolerance{};
  PisoIntermediateCertificate corrector_one{};
  PisoStateCorrectionCertificate corrected_c1{};
  PisoStateCorrectionCertificate current_corrected_c1{};
  PisoIntermediateCertificate refinement_root_c1{};
  PisoStateCorrectionCertificate refinement_predecessor{};
  PisoPressureEnergyRefinementStateCertificate refinement_authority{};
  PisoIntermediateCertificate current{};
  PisoIntermediateCertificate terminal_lineage_source{};
  ConstFieldView current_trial_velocity{};
  ConstFieldView current_pressure_perturbation{};
  PlanFingerprint current_pressure_perturbation_numeric_fingerprint{};
  PressureCorrectionCertificate current_pressure_work{};
  BoundaryThermophysicalGhostContext current_thermophysical_context{};
  PressureReferenceCertificate current_pressure_reference{};
  double current_absolute_pressure_reference{};
  PisoTerminalCertificate terminal{};
  const IbmEquationInterfacePlan* immersed_interface{};
  PressureContinuityActivityView continuity_activity{};
  RemoteDonorExchangePlan* pressure_correction_donors{};
  StageId pressure_correction_donor_stage{};
  RemoteDonorExchangePlan* candidate_pressure_correction_donors{};
  StageId candidate_pressure_correction_donor_stage{};
  FieldId candidate_pressure_correction_field{};
  std::uint8_t candidate_pressure_correction_donor_reach{};
  PlanFingerprint candidate_pressure_correction_donor_fingerprint{};
  PressureCorrectionInput pressure_input{};
  PressureCorrectionCertificate pressure_correction{};
  SealedPressureCorrectionAuthority sealed{};
  PredecessorPressureApplication predecessor_c1{};
  CorrectedPendingAuthority corrected_pending{};
  std::unique_ptr<double[]> frozen_candidate_density_storage{};
  std::unique_ptr<double[]> frozen_candidate_face_aux_storage{};
  mutable std::unique_ptr<std::uint64_t[]>
      frozen_candidate_collective_hashes{};
  mutable std::unique_ptr<double[]> frozen_candidate_composition_values{};
  std::unique_ptr<FieldId[]> independent_species_semantic_fields{};
  std::size_t independent_species_count{};
  FieldView frozen_candidate_density{};
  FaceFluxView frozen_candidate_face_aux{};
  BdfCoefficients frozen_candidate_bdf{};
  RevisionToken frozen_candidate_baseline{};
  PlanFingerprint frozen_candidate_numeric{};
  bool frozen_candidate_cartesian{};
  mutable PlanFingerprint current_frozen_exact_lineage{};
  mutable PlanFingerprint current_frozen_exact_scratch{};

  static constexpr std::size_t kCandidateProvenanceBatchCapacity = 4U;
  static constexpr std::uint8_t kCandidateProvenanceNoDependency =
      std::numeric_limits<std::uint8_t>::max();

  struct CandidateProvenanceBatchEntry {
    std::uint64_t local{};
    std::uint64_t domain{};
    std::uint8_t dependency{kCandidateProvenanceNoDependency};
    std::uint64_t suffix{};
  };

  Status collective_candidate_provenance_batch(
      const CandidateProvenanceBatchEntry* entries, std::size_t entry_count,
      PlanFingerprint* provenances) const noexcept {
    if (communicator == MPI_COMM_NULL || size <= 0 || entries == nullptr ||
        provenances == nullptr || entry_count == 0U ||
        entry_count > kCandidateProvenanceBatchCapacity ||
        !frozen_candidate_collective_hashes) {
      return {StatusCode::invalid_plan, kPisoCoupler};
    }
    for (std::size_t entry = 0U; entry < entry_count; ++entry)
      provenances[entry] = 0U;
    bool has_dependencies = false;
    std::array<std::uint64_t, 2U * kCandidateProvenanceBatchCapacity>
        local_payload{};
    for (std::size_t entry = 0U; entry < entry_count; ++entry) {
      if (entries[entry].dependency != kCandidateProvenanceNoDependency &&
          entries[entry].dependency >= entry) {
        return {StatusCode::invalid_plan, kPisoCoupler};
      }
      has_dependencies =
          has_dependencies ||
          entries[entry].dependency != kCandidateProvenanceNoDependency;
    }
    // A dependent rank-local fingerprint is mixed in the original order:
    // raw local bytes, the already-folded collective dependency, then that
    // rank's suffix.  Gather the suffix beside the raw hash so even a rejected
    // rank-divergent alpha produces exactly the same fold as the scalar route.
    const std::size_t payload_stride = has_dependencies ? 2U : 1U;
    for (std::size_t entry = 0U; entry < entry_count; ++entry) {
      local_payload[payload_stride * entry] = entries[entry].local;
      if (has_dependencies)
        local_payload[payload_stride * entry + 1U] = entries[entry].suffix;
    }
    const std::size_t payload_count = payload_stride * entry_count;
    const int batch_count = static_cast<int>(payload_count);
    const int gather_result = MPI_Allgather(
        local_payload.data(), batch_count, MPI_UINT64_T,
        frozen_candidate_collective_hashes.get(), batch_count, MPI_UINT64_T,
        communicator);
    const int local_failure = gather_result == MPI_SUCCESS ? 0 : 1;
    int global_failure = 0;
    const int consensus_result = MPI_Allreduce(
        &local_failure, &global_failure, 1, MPI_INT, MPI_MAX, communicator);
    if (consensus_result != MPI_SUCCESS || global_failure != 0) {
      return {StatusCode::mpi_failure, kPisoCollective};
    }
    for (std::size_t entry = 0U; entry < entry_count; ++entry) {
      std::uint64_t collective = hash_mix(kFnvOffset, entries[entry].domain);
      collective = hash_mix(collective, static_cast<std::uint64_t>(size));
      for (int rank_index = 0; rank_index < size; ++rank_index) {
        const std::size_t rank_offset =
            static_cast<std::size_t>(rank_index) * payload_count;
        const std::size_t entry_offset = payload_stride * entry;
        std::uint64_t rank_local = frozen_candidate_collective_hashes[
            rank_offset + entry_offset];
        if (entries[entry].dependency !=
            kCandidateProvenanceNoDependency) {
          rank_local = hash_mix(
              rank_local, provenances[entries[entry].dependency]);
          rank_local = hash_mix(
              rank_local, frozen_candidate_collective_hashes[
                              rank_offset + entry_offset + 1U]);
        }
        collective = hash_mix(collective, rank_local);
      }
      provenances[entry] = collective == 0U ? 1U : collective;
    }
    return {};
  }

  Status collective_candidate_provenance(
      std::uint64_t local, std::uint64_t domain,
      PlanFingerprint& provenance) const noexcept {
    const CandidateProvenanceBatchEntry entry{local, domain};
    return collective_candidate_provenance_batch(&entry, 1U, &provenance);
  }

  PlanFingerprint frozen_candidate_numeric_fingerprint() const noexcept {
    if (!current.valid() || frozen_candidate_baseline == 0U ||
        !detail::valid_bdf_coefficients(frozen_candidate_bdf))
      return 0U;
    std::uint64_t hash =
        hash_mix(kFnvOffset, UINT64_C(0x76303466726f7a63));
    hash = hash_mix(hash, fingerprint);
    hash = hash_mix(hash, current.dependency);
    hash = hash_mix(hash, current.corrector);
    hash = hash_mix(hash, double_bits(frozen_candidate_bdf.a0));
    hash = hash_mix(hash, double_bits(frozen_candidate_bdf.a1));
    hash = hash_mix(hash, double_bits(frozen_candidate_bdf.a2));
    hash = hash_mix(hash, frozen_candidate_bdf.order);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    // Keep the diagnostic core's byte-exact mutation detector.  Production
    // fields are mutated only through revision-issuing authorities, so the
    // complete binding below provides the same frozen-candidate guard without
    // rescanning these large arrays at every consumer boundary.
    hash = mix_candidate_field_values(
        hash, as_const(frozen_candidate_density), cells, 1U, 1);
    hash = mix_candidate_field_values(hash, as_const(workspace.h_by_a), cells,
                                      3U, 1);
    hash = mix_candidate_field_values(hash, as_const(workspace.r_au), cells,
                                      3U, 1);
    if (current.corrector == 1U) {
      if (!detail::valid_cell_view(current_pressure_perturbation, cells, 0U,
                                   1U, 1U))
        return 0U;
      hash = hash_mix(
          hash, scalar_reach_one_numeric_fingerprint(
                    current_pressure_perturbation, cells));
      hash = mix_candidate_field_values(
          hash, as_const(workspace.pressure_gradient), cells, 3U, 0);
    }
    hash = mix_candidate_flux_values(hash, as_const(workspace.phi_h_by_a));
    hash = mix_candidate_flux_values(
        hash, as_const(frozen_candidate_face_aux));
#else
    hash = mix_complete_view_identity(
        hash, as_const(frozen_candidate_density));
    hash = mix_complete_view_identity(hash, as_const(workspace.h_by_a));
    hash = mix_complete_view_identity(hash, as_const(workspace.r_au));
    if (current.corrector == 1U) {
      if (!detail::valid_cell_view(current_pressure_perturbation, cells, 0U,
                                   1U, 1U))
        return 0U;
      hash = mix_complete_view_identity(hash, current_pressure_perturbation);
      hash = mix_complete_view_identity(
          hash, as_const(workspace.pressure_gradient));
    }
    hash = mix_candidate_flux_binding(
        hash, as_const(workspace.phi_h_by_a));
    hash = mix_candidate_flux_binding(
        hash, as_const(frozen_candidate_face_aux));
#endif
    if (immersed_interface != nullptr ||
        candidate_pressure_correction_donors != nullptr) {
      hash = hash_mix(
          hash, immersed_interface == nullptr
                    ? PlanFingerprint{0U}
                    : immersed_interface->fingerprint());
      hash = hash_mix(hash, continuity_activity.local_fingerprint);
      hash = hash_mix(hash, continuity_activity.collective_fingerprint);
      hash = hash_mix(hash,
                      candidate_pressure_correction_donor_fingerprint);
      hash = hash_mix(hash, candidate_pressure_correction_donor_stage);
      hash = hash_mix(hash, candidate_pressure_correction_field);
      hash = hash_mix(hash, candidate_pressure_correction_donor_reach);
    }
    return hash == 0U ? PlanFingerprint{1U} : hash;
  }
};

PressureVelocityCoupler::~PressureVelocityCoupler() noexcept { release(); }

PressureVelocityCoupler::PressureVelocityCoupler(
    PressureVelocityCoupler&& other) noexcept
    : implementation_(other.implementation_) {
  other.implementation_ = nullptr;
}

PressureVelocityCoupler& PressureVelocityCoupler::operator=(
    PressureVelocityCoupler&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = other.implementation_;
    other.implementation_ = nullptr;
  }
  return *this;
}

void PressureVelocityCoupler::release() noexcept {
  delete implementation_;
  implementation_ = nullptr;
}

Status PressureVelocityCoupler::bind(
    const PisoPlan& plan, const EquationPlanSet& equations,
    PisoCouplerServices services,
    PisoCouplerWorkspace workspace, PressureVelocityCoupler& out) noexcept {
  int communicator_rank = -1;
  int communicator_size = 0;
  const bool valid_communicator =
      services.communicator != MPI_COMM_NULL &&
      MPI_Comm_rank(services.communicator, &communicator_rank) == MPI_SUCCESS &&
      MPI_Comm_size(services.communicator, &communicator_size) == MPI_SUCCESS &&
      communicator_rank >= 0 && communicator_rank < communicator_size &&
      communicator_size > 0;
  const Int3 cells = equations.pressure_reference().cells();
  const bool valid_cells = same_cells(cells, plan.cells_);
  const bool valid_cell_workspace =
      detail::valid_cell_view(workspace.r_au, cells, 0U, 3U) &&
      workspace.r_au.ghosts.x >= 1 && workspace.r_au.ghosts.y >= 1 &&
      workspace.r_au.ghosts.z >= 1 &&
      detail::valid_cell_view(workspace.h_by_a, cells, 0U, 3U) &&
      workspace.h_by_a.ghosts.x >=
          static_cast<std::int32_t>(equations.kernels().reach()) &&
      workspace.h_by_a.ghosts.y >=
          static_cast<std::int32_t>(equations.kernels().reach()) &&
      workspace.h_by_a.ghosts.z >=
          static_cast<std::int32_t>(equations.kernels().reach()) &&
      detail::valid_cell_view(workspace.pressure_gradient, cells, 0U, 3U) &&
      workspace.pressure_gradient.ghosts.x >= 1 &&
      workspace.pressure_gradient.ghosts.y >= 1 &&
      workspace.pressure_gradient.ghosts.z >= 1 &&
      !detail::field_views_overlap(as_const(workspace.r_au),
                                   as_const(workspace.h_by_a)) &&
      !detail::field_views_overlap(as_const(workspace.r_au),
                                   as_const(workspace.pressure_gradient)) &&
      !detail::field_views_overlap(as_const(workspace.h_by_a),
                                   as_const(workspace.pressure_gradient));
  const bool valid_faces =
      valid_face_workspace(workspace.x_pressure_coefficient,
                           CartesianAxis::x, cells) &&
      valid_face_workspace(workspace.y_pressure_coefficient,
                           CartesianAxis::y, cells) &&
      valid_face_workspace(workspace.z_pressure_coefficient,
                           CartesianAxis::z, cells) &&
      !detail::face_views_overlap(workspace.x_pressure_coefficient,
                                  workspace.y_pressure_coefficient) &&
      !detail::face_views_overlap(workspace.x_pressure_coefficient,
                                  workspace.z_pressure_coefficient) &&
      !detail::face_views_overlap(workspace.y_pressure_coefficient,
                                  workspace.z_pressure_coefficient);
  const ConstFaceFluxView phi = as_const(workspace.phi_h_by_a);
  const bool valid_phi =
      detail::valid_flux_view(phi, cells, phi.revision) &&
      !phi.certificate.valid();
  const FaceFieldView pressure_faces[]{workspace.x_pressure_coefficient,
                                       workspace.y_pressure_coefficient,
                                       workspace.z_pressure_coefficient};
  const ConstFaceFieldView phi_faces[]{phi.x, phi.y, phi.z};
  bool disjoint = valid_cell_workspace && valid_faces && valid_phi;
  for (FaceFieldView pressure_face : pressure_faces) {
    disjoint = disjoint &&
               !cell_aliases_face(as_const(workspace.r_au), pressure_face) &&
               !cell_aliases_face(as_const(workspace.h_by_a), pressure_face) &&
               !cell_aliases_face(as_const(workspace.pressure_gradient),
                                  pressure_face);
    for (ConstFaceFieldView phi_face : phi_faces) {
      disjoint = disjoint &&
                 !detail::face_views_overlap(pressure_face, phi_face);
    }
  }
  for (ConstFaceFieldView phi_face : phi_faces) {
    disjoint = disjoint &&
               !detail::cell_face_views_overlap(workspace.r_au, phi_face) &&
               !detail::cell_face_views_overlap(workspace.h_by_a, phi_face) &&
               !detail::cell_face_views_overlap(workspace.pressure_gradient,
                                                phi_face);
  }
  const std::array<HaloFieldSpec, 4U> halo_fields{{
      {services.density_field, 1U, 1U},
      {workspace.r_au.field, 1U, 3U},
      {workspace.h_by_a.field, equations.kernels().reach(), 3U},
      {workspace.pressure_gradient.field, 1U, 3U}}};
  const std::array<HaloFieldSpec, 1U> correction_halo_fields{{
      {services.correction_field, 1U, 1U}}};
  const bool valid_services =
      valid_communicator && services.geometry != nullptr &&
      services.geometry->fingerprint() != 0U && services.boundary != nullptr &&
      services.boundary->semantic_fingerprint() != 0U &&
      services.thermodynamics != nullptr &&
      services.thermodynamics->fingerprint() ==
          equations.thermodynamics_fingerprint() &&
      services.halo != nullptr && services.halo->ready() &&
      services.halo_stage != 0U &&
      services.density_field != workspace.r_au.field &&
      services.density_field != workspace.h_by_a.field &&
      services.density_field != workspace.pressure_gradient.field &&
      workspace.r_au.field != workspace.pressure_gradient.field &&
      workspace.h_by_a.field != workspace.pressure_gradient.field &&
      services.correction_halo != nullptr &&
      services.correction_halo->ready() &&
      services.correction_halo_stage != 0U &&
      same_cells(services.patch.cells, cells) &&
      same_cells(services.geometry->global_cells(), equations.global_cells()) &&
      same_cells(services.boundary->local_cells(), cells) &&
      static_cast<bool>(services.halo->validate_contract(
          services.communicator, services.patch,
          {halo_fields.data(), halo_fields.size()},
          services.boundary->halo_topology())) &&
      static_cast<bool>(services.correction_halo->validate_contract(
          services.communicator, services.patch,
          {correction_halo_fields.data(), correction_halo_fields.size()},
          services.boundary->halo_topology()));
  const std::size_t cell_count = static_cast<std::size_t>(cells.x) *
                                 static_cast<std::size_t>(cells.y) *
                                 static_cast<std::size_t>(cells.z);
  const std::size_t x_face_count =
      static_cast<std::size_t>(cells.x + 1) *
      static_cast<std::size_t>(cells.y) * static_cast<std::size_t>(cells.z);
  const std::size_t y_face_count = static_cast<std::size_t>(cells.x) *
                                   static_cast<std::size_t>(cells.y + 1) *
                                   static_cast<std::size_t>(cells.z);
  const std::size_t z_face_count = static_cast<std::size_t>(cells.x) *
                                   static_cast<std::size_t>(cells.y) *
                                   static_cast<std::size_t>(cells.z + 1);
  const auto binary_span = [](Span<const std::uint8_t> values,
                              std::size_t expected) noexcept {
    if (values.data == nullptr || values.size != expected) return false;
    for (std::size_t index = 0U; index < values.size; ++index)
      if (values.data[index] > 1U) return false;
    return true;
  };
  const PressureContinuityActivityView activity =
      services.continuity_activity;
  const bool empty_activity =
      activity.cells.data == nullptr && activity.cells.size == 0U &&
      activity.x_faces.data == nullptr && activity.x_faces.size == 0U &&
      activity.y_faces.data == nullptr && activity.y_faces.size == 0U &&
      activity.z_faces.data == nullptr && activity.z_faces.size == 0U &&
      activity.local_fingerprint == 0U &&
      activity.collective_fingerprint == 0U;
  const bool complete_activity =
      binary_span(activity.cells, cell_count) &&
      binary_span(activity.x_faces, x_face_count) &&
      binary_span(activity.y_faces, y_face_count) &&
      binary_span(activity.z_faces, z_face_count) &&
      activity.local_fingerprint != 0U &&
      activity.collective_fingerprint != 0U;
  const bool valid_pressure_donors =
      (services.pressure_correction_donors == nullptr &&
       services.pressure_correction_donor_stage == 0U) ||
      (services.pressure_correction_donors != nullptr &&
       services.pressure_correction_donors->ready() &&
       services.pressure_correction_donor_stage != 0U);
  const bool empty_candidate_pressure_donors =
      services.candidate_pressure_correction_donors == nullptr &&
      services.candidate_pressure_correction_donor_stage == 0U &&
      services.candidate_pressure_correction_field == 0U &&
      services.candidate_pressure_correction_donor_reach == 0U &&
      services.candidate_pressure_correction_donor_fingerprint == 0U;
  const bool complete_candidate_pressure_donors =
      services.candidate_pressure_correction_donors != nullptr &&
      services.candidate_pressure_correction_donors->ready() &&
      services.candidate_pressure_correction_donors !=
          services.pressure_correction_donors &&
      services.candidate_pressure_correction_donor_stage != 0U &&
      services.candidate_pressure_correction_donor_stage !=
          services.pressure_correction_donor_stage &&
      services.candidate_pressure_correction_field != 0U &&
      services.candidate_pressure_correction_field !=
          services.correction_field &&
      services.candidate_pressure_correction_donor_reach >= 1U &&
      services.candidate_pressure_correction_donor_reach ==
          services.candidate_pressure_correction_donors->reach() &&
      services.candidate_pressure_correction_donor_fingerprint != 0U &&
      services.candidate_pressure_correction_donor_fingerprint ==
          services.candidate_pressure_correction_donors->fingerprint();
  const bool valid_candidate_pressure_donors =
      empty_candidate_pressure_donors ||
      (complete_candidate_pressure_donors &&
       services.pressure_correction_donors != nullptr &&
       services.immersed_interface != nullptr &&
       services.immersed_interface->fingerprint() != 0U);
  if (plan.fingerprint_ == 0U || equations.fingerprint() == 0U ||
      equations.fingerprint() != plan.equations_fingerprint_ || !valid_cells ||
      !disjoint || !valid_services ||
      (!empty_activity && !complete_activity) || !valid_pressure_donors ||
      !valid_candidate_pressure_donors) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  PressureCorrectionBoundaryPlan pressure_boundary;
  const Status pressure_boundary_status =
      PressureCorrectionBoundaryPlan::compile(
          *services.geometry, services.patch, *services.boundary,
          pressure_boundary);
  if (!pressure_boundary_status) return pressure_boundary_status;
  auto* candidate = new (std::nothrow) Impl;
  if (candidate == nullptr) {
    return {StatusCode::allocation_failure, kPisoCoupler};
  }
  const std::size_t frozen_nx = static_cast<std::size_t>(cells.x + 2);
  const std::size_t frozen_ny = static_cast<std::size_t>(cells.y + 2);
  const std::size_t frozen_nz = static_cast<std::size_t>(cells.z + 2);
  const std::size_t frozen_density_count = frozen_nx * frozen_ny * frozen_nz;
  const std::size_t frozen_face_count = x_face_count + y_face_count +
                                        z_face_count;
  candidate->frozen_candidate_density_storage.reset(
      new (std::nothrow) double[frozen_density_count]);
  candidate->frozen_candidate_face_aux_storage.reset(
      new (std::nothrow) double[frozen_face_count]);
  candidate->frozen_candidate_collective_hashes.reset(
      new (std::nothrow) std::uint64_t[
          static_cast<std::size_t>(communicator_size) *
          2U * Impl::kCandidateProvenanceBatchCapacity]);
  candidate->independent_species_count =
      services.thermodynamics->independent_species_count();
  if (candidate->independent_species_count != 0U) {
    candidate->frozen_candidate_composition_values.reset(
        new (std::nothrow) double[candidate->independent_species_count]);
    candidate->independent_species_semantic_fields.reset(
        new (std::nothrow) FieldId[candidate->independent_species_count]);
  }
  if (!candidate->frozen_candidate_density_storage ||
      !candidate->frozen_candidate_face_aux_storage ||
      !candidate->frozen_candidate_collective_hashes ||
      (candidate->independent_species_count != 0U &&
       (!candidate->frozen_candidate_composition_values ||
        !candidate->independent_species_semantic_fields))) {
    delete candidate;
    return {StatusCode::allocation_failure, kPisoCoupler};
  }
  std::size_t semantic_species = 0U;
  const Span<const BoundaryTransportedField> transported =
      services.boundary->transported_fields();
  if (equations.species().size() != candidate->independent_species_count) {
    delete candidate;
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  for (std::size_t index = 0U; index < transported.size; ++index) {
    if (transported.data[index].role != BoundaryScalarRole::species) continue;
    const ScalarEquationSpec* expected =
        equations.species().spec(semantic_species);
    if (semantic_species >= candidate->independent_species_count ||
        expected == nullptr ||
        expected->role != TransportedScalarRole::species ||
        expected->field == 0U ||
        transported.data[index].field != expected->field) {
      delete candidate;
      return {StatusCode::invalid_plan, kPisoCoupler};
    }
    candidate->independent_species_semantic_fields[semantic_species++] =
        expected->field;
  }
  if (semantic_species != candidate->independent_species_count) {
    delete candidate;
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  double* const frozen_density_base =
      candidate->frozen_candidate_density_storage.get();
  const StorageIdentity frozen_density_storage =
      reinterpret_cast<std::uintptr_t>(frozen_density_base) == 0U
          ? StorageIdentity{1U}
          : static_cast<StorageIdentity>(
                reinterpret_cast<std::uintptr_t>(frozen_density_base));
  candidate->frozen_candidate_density = {
      frozen_density_base + 1U + frozen_nx + frozen_nx * frozen_ny,
      cells,
      {1, 1, 1},
      1U,
      frozen_nx,
      frozen_nx * frozen_ny,
      frozen_density_count,
      0U,
      services.density_field,
      1U,
      frozen_density_storage,
      static_cast<RevisionDomainIdentity>(frozen_density_storage ^
                                          UINT64_C(0x66726f7a656e7268))};
  if (candidate->frozen_candidate_density.revision_domain == 0U)
    candidate->frozen_candidate_density.revision_domain = 1U;
  double* const frozen_face_base =
      candidate->frozen_candidate_face_aux_storage.get();
  const StorageIdentity frozen_face_storage =
      reinterpret_cast<std::uintptr_t>(frozen_face_base) == 0U
          ? StorageIdentity{1U}
          : static_cast<StorageIdentity>(
                reinterpret_cast<std::uintptr_t>(frozen_face_base));
  RevisionDomainIdentity frozen_face_domain =
      static_cast<RevisionDomainIdentity>(frozen_face_storage ^
                                          UINT64_C(0x66726f7a656e6661));
  if (frozen_face_domain == 0U) frozen_face_domain = 1U;
  candidate->frozen_candidate_face_aux = {
      {frozen_face_base,
       {cells.x + 1, cells.y, cells.z},
       static_cast<std::size_t>(cells.x + 1),
       static_cast<std::size_t>(cells.x + 1) * cells.y,
       CartesianAxis::x,
       frozen_face_storage,
       frozen_face_domain},
      {frozen_face_base + x_face_count,
       {cells.x, cells.y + 1, cells.z},
       static_cast<std::size_t>(cells.x),
       static_cast<std::size_t>(cells.x) * (cells.y + 1),
       CartesianAxis::y,
       frozen_face_storage,
       frozen_face_domain},
      {frozen_face_base + x_face_count + y_face_count,
       {cells.x, cells.y, cells.z + 1},
       static_cast<std::size_t>(cells.x),
       static_cast<std::size_t>(cells.x) * cells.y,
       CartesianAxis::z,
       frozen_face_storage,
       frozen_face_domain},
      1U,
      {}};
  candidate->kernels = &equations.kernels();
  candidate->geometry = services.geometry;
  candidate->boundary = services.boundary;
  candidate->thermodynamics_plan = services.thermodynamics;
  candidate->pressure_reference = &equations.pressure_reference();
  candidate->communicator = services.communicator;
  candidate->rank = communicator_rank;
  candidate->size = communicator_size;
  candidate->pressure_boundary = std::move(pressure_boundary);
  candidate->halo = services.halo;
  candidate->correction_halo = services.correction_halo;
  candidate->patch = services.patch;
  candidate->halo_stage = services.halo_stage;
  candidate->correction_halo_stage = services.correction_halo_stage;
  candidate->density_field = services.density_field;
  candidate->correction_field = services.correction_field;
  candidate->continuity_activity = activity;
  candidate->pressure_correction_donors =
      services.pressure_correction_donors;
  candidate->pressure_correction_donor_stage =
      services.pressure_correction_donor_stage;
  candidate->candidate_pressure_correction_donors =
      services.candidate_pressure_correction_donors;
  candidate->candidate_pressure_correction_donor_stage =
      services.candidate_pressure_correction_donor_stage;
  candidate->candidate_pressure_correction_field =
      services.candidate_pressure_correction_field;
  candidate->candidate_pressure_correction_donor_reach =
      services.candidate_pressure_correction_donor_reach;
  candidate->candidate_pressure_correction_donor_fingerprint =
      services.candidate_pressure_correction_donor_fingerprint;
  candidate->immersed_interface = services.immersed_interface;
  candidate->workspace = workspace;
  candidate->cells = cells;
  candidate->plan = plan.fingerprint_;
  candidate->momentum_plan = equations.momentum().fingerprint();
  candidate->predictor_plan = plan.predictor_fingerprint_;
  candidate->pressure_reference_plan =
      plan.pressure_reference_fingerprint_;
  candidate->thermodynamics = equations.thermodynamics_fingerprint();
  candidate->transport = equations.transport_fingerprint();
  candidate->pressure_reference_kind =
      equations.pressure_reference().kind();
  candidate->pressure_gauge = equations.pressure_reference().gauge();
  candidate->equations = equations.fingerprint();
  candidate->pressure_algorithm = plan.pressure_algorithm_;
  candidate->mg_correction_scaling = plan.mg_correction_scaling_;
  candidate->eos_tolerance = plan.eos_tolerance_;
  candidate->continuity_tolerance = plan.continuity_tolerance_;
  candidate->energy_tolerance = plan.energy_tolerance_;
  candidate->closed_mass_tolerance = plan.closed_mass_tolerance_;
  candidate->gauge_tolerance = plan.gauge_tolerance_;
  std::uint64_t hash = kFnvOffset;
  hash = hash_mix(hash, candidate->plan);
  hash = mix_view(hash, as_const(workspace.r_au));
  hash = mix_view(hash, as_const(workspace.h_by_a));
  hash = mix_view(hash, as_const(workspace.pressure_gradient));
  hash = mix_face(hash, as_const(workspace.x_pressure_coefficient));
  hash = mix_face(hash, as_const(workspace.y_pressure_coefficient));
  hash = mix_face(hash, as_const(workspace.z_pressure_coefficient));
  hash = mix_face(hash, phi.x);
  hash = mix_face(hash, phi.y);
  hash = mix_face(hash, phi.z);
  hash = hash_mix(
      hash, candidate->pressure_boundary.certificate().semantic);
  hash = hash_mix(
      hash, candidate->pressure_boundary.certificate().rank_local_layout);
  if (candidate->pressure_correction_donors != nullptr) {
    hash = hash_mix(hash,
                    candidate->pressure_correction_donors->fingerprint());
    hash = hash_mix(hash, candidate->pressure_correction_donor_stage);
  }
  if (candidate->candidate_pressure_correction_donors != nullptr) {
    hash = hash_mix(
        hash, candidate->candidate_pressure_correction_donor_fingerprint);
    hash = hash_mix(
        hash, candidate->candidate_pressure_correction_donor_stage);
    hash = hash_mix(hash, candidate->candidate_pressure_correction_field);
    hash = hash_mix(
        hash, candidate->candidate_pressure_correction_donor_reach);
  }
  candidate->fingerprint = hash == 0U ? 1U : hash;
  out.release();
  out.implementation_ = candidate;
  return {};
}

Status PressureVelocityCoupler::refresh(
    const PisoIntermediateInput& input,
    PisoIntermediateCertificate& certificate,
    Status prerequisite) noexcept {
  return refresh_impl(nullptr, input, nullptr, certificate, prerequisite);
}

Status PressureVelocityCoupler::refresh_pressure_energy_refinement(
    const PisoPressureEnergyRefinementStateCertificate& refinement,
    const PisoIntermediateInput& input, PisoCoupledStateView state,
    PisoIntermediateCertificate& certificate,
    Status prerequisite) noexcept {
  return refresh_impl(&refinement, input, &state, certificate, prerequisite);
}

Status PressureVelocityCoupler::refresh_impl(
    const PisoPressureEnergyRefinementStateCertificate* refinement,
    const PisoIntermediateInput& input, const PisoCoupledStateView* state,
    PisoIntermediateCertificate& certificate,
    Status prerequisite) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  Impl& impl = *implementation_;
  impl.current = {};
  impl.current_trial_velocity = {};
  impl.current_pressure_perturbation = {};
  impl.current_pressure_perturbation_numeric_fingerprint = 0U;
  impl.current_pressure_work = {};
  impl.current_thermophysical_context = {};
  impl.current_pressure_reference = {};
  impl.current_absolute_pressure_reference = 0.0;
  impl.current_corrected_c1 = {};
  impl.sealed = {};
  impl.corrected_pending = {};
  impl.frozen_candidate_bdf = {};
  impl.frozen_candidate_baseline = 0U;
  impl.frozen_candidate_numeric = 0U;
  impl.current_frozen_exact_lineage = 0U;
  impl.current_frozen_exact_scratch = 0U;
  impl.frozen_candidate_cartesian = false;
  Status local = prerequisite;
  PisoIntermediateCertificate corrector_one_authorization{};
  PisoStateCorrectionCertificate corrected_c1_authorization{};
  const bool refinement_entry = refinement != nullptr;
  std::uint8_t refinement_iteration = 0U;
  PlanFingerprint refinement_collective_lineage = 0U;
  PlanFingerprint refinement_lineage = 0U;
  bool paired_flux_consumed = false;
  if (input.corrector == 1U) {
    impl.corrector_one = {};
    impl.corrected_c1 = {};
    impl.refinement_root_c1 = {};
    impl.refinement_predecessor = {};
    impl.refinement_authority = {};
    impl.terminal_lineage_source = {};
    impl.predecessor_c1 = {};
    if (local &&
        (refinement_entry ||
         input.immersed_interface != impl.immersed_interface))
      local = {StatusCode::invalid_plan, kPisoCoupler};
  } else if (input.corrector == 2U) {
    if (refinement_entry) {
      const bool state_views_valid =
          state != nullptr &&
          valid_pressure_energy_refinement_state_views(*state, impl.cells) &&
          same_piso_field_revision_identity(refinement->velocity_,
                                            as_const(state->velocity)) &&
          same_piso_field_revision_identity(
              refinement->pressure_, as_const(state->pressure_perturbation)) &&
          same_piso_field_revision_identity(refinement->enthalpy_,
                                            as_const(state->enthalpy)) &&
          same_piso_field_revision_identity(refinement->density_,
                                            as_const(state->density)) &&
          same_piso_field_revision_identity(refinement->temperature_,
                                            as_const(state->temperature)) &&
          same_field_identity(as_const(state->velocity),
                              input.trial_velocity) &&
          same_field_identity(as_const(state->density),
                              as_const(input.density)) &&
          same_field_identity(
              as_const(state->pressure_perturbation),
              input.thermophysical_boundary.binding.pressure_perturbation) &&
          same_field_identity(
              as_const(state->enthalpy),
              input.thermophysical_boundary.binding.enthalpy) &&
          same_field_identity(
              as_const(state->density),
              input.thermophysical_boundary.binding.density);
      const bool flux_views_valid =
          detail::valid_flux_view(refinement->flux_, impl.cells,
                                  refinement->flux_.revision) &&
          detail::valid_flux_view(input.trial_flux, impl.cells,
                                  input.trial_flux.revision) &&
          same_candidate_flux_identity(refinement->flux_, input.trial_flux) &&
          state_views_valid &&
          refinement_state_flux_alias_free(*state, input.trial_flux);
      const PlanFingerprint local_state_numeric =
          local && state_views_valid
              ? pressure_energy_refinement_state_numeric(*state, impl.cells)
              : PlanFingerprint{};
      const PlanFingerprint local_flux_numeric =
          local && flux_views_valid
              ? pressure_energy_refinement_flux_numeric(input.trial_flux)
              : PlanFingerprint{};
      const bool authority_matches =
          refinement->valid() && refinement->issuer_ == this &&
          impl.refinement_authority.valid() &&
          impl.refinement_authority.collective_lineage_ ==
              refinement->collective_lineage_ &&
          impl.refinement_authority.lineage_ == refinement->lineage_ &&
          same_state_correction_certificate(
              impl.refinement_authority.predecessor_,
              refinement->predecessor_) &&
          same_state_correction_certificate(impl.refinement_predecessor,
                                            refinement->predecessor_) &&
          state_views_valid &&
          flux_views_valid &&
          input.prior_corrector == refinement->predecessor_.state &&
          input.pressure_reference.pressure_reference ==
              refinement->predecessor_.output_pressure_reference
                  .pressure_reference &&
          local_state_numeric == refinement->rank_local_state_numeric_ &&
          local_flux_numeric == refinement->rank_local_flux_numeric_;
      corrector_one_authorization = impl.refinement_root_c1;
      corrected_c1_authorization = refinement->predecessor_;
      refinement_iteration = refinement->iteration_;
      refinement_collective_lineage = refinement->collective_lineage_;
      refinement_lineage = refinement->lineage_;
      impl.refinement_authority = {};
      impl.refinement_predecessor = {};
      if (local && !authority_matches)
        local = {StatusCode::invalid_plan, kPisoCoupler};
    } else {
      // The ordinary corrector-one authorization remains single-use.  Only a
      // typed provisional C2 commit may authorize a later refinement entry.
      corrector_one_authorization = impl.corrector_one;
      corrected_c1_authorization = impl.corrected_c1;
      impl.refinement_root_c1 = corrector_one_authorization;
    }
    impl.corrector_one = {};
    impl.corrected_c1 = {};
    impl.terminal_lineage_source = {};
    if (local &&
        (!corrector_one_authorization.valid() ||
         input.prior_corrector !=
             (corrected_c1_authorization.valid()
                  ? corrected_c1_authorization.state
                  : corrector_one_authorization.dependency) ||
         input.immersed_interface != impl.immersed_interface))
      local = {StatusCode::invalid_plan, kPisoCoupler};
  } else if (local || refinement_entry) {
    local = {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Int3 cells = impl.cells;
  const bool valid_certificates =
      impl.pressure_boundary.current() &&
      input.momentum.valid() &&
      input.momentum.plan == impl.momentum_plan &&
      input.momentum.scope == EquationAssemblyScope::momentum_predictor &&
      input.predictor.valid() && input.predictor.plan == impl.predictor_plan &&
      input.pressure_reference.valid() &&
      input.pressure_reference.plan == impl.pressure_reference_plan &&
      input.pressure_reference.kind == impl.pressure_reference_kind &&
      input.momentum.time == input.predictor.time &&
      input.momentum.time == input.pressure_reference.time &&
      input.momentum.geometry == input.predictor.geometry &&
      input.momentum.geometry ==
          impl.pressure_boundary.certificate().geometry;
  const auto empty_flux = [](ConstFaceFluxView flux) noexcept {
    return flux.x.base == nullptr && flux.y.base == nullptr &&
           flux.z.base == nullptr && flux.revision == 0U &&
           !flux.certificate.valid();
  };
  const ConstFaceFluxView phi_output = as_const(impl.workspace.phi_h_by_a);
  const bool valid_temporal_reference =
      input.corrector == 1U
          ? detail::valid_flux_view(input.temporal_reference, cells,
                                    input.temporal_reference.revision) &&
                !input.temporal_reference.certificate.valid() &&
                input.temporal_reference.revision !=
                    input.trial_flux.revision &&
                input.temporal_reference.x.base == phi_output.x.base &&
                input.temporal_reference.y.base == phi_output.y.base &&
                input.temporal_reference.z.base == phi_output.z.base &&
                input.temporal_reference.x.storage_identity ==
                    phi_output.x.storage_identity &&
                input.temporal_reference.x.revision_domain ==
                    phi_output.x.revision_domain
          : empty_flux(input.temporal_reference);
  const auto valid_committed_flux = [&](ConstFaceFluxView flux,
                                        RevisionToken revision) noexcept {
    return revision != 0U && flux.revision == revision &&
           detail::valid_flux_view(flux, cells, revision) &&
           flux.certificate.valid() && flux.certificate.matches(flux);
  };
  const ConstFaceFluxView accepted_history =
      input.committed_face_history.accepted;
  const ConstFaceFluxView previous_history =
      input.committed_face_history.previous;
  const bool valid_committed_face_history =
      input.corrector == 1U
          ? valid_committed_flux(accepted_history,
                                 input.predictor.accepted_face_flux) &&
                accepted_history.certificate.authority() ==
                    input.predictor.committed_face_flux_authority &&
                accepted_history.certificate.storage() ==
                    input.predictor.committed_face_flux_storage &&
                accepted_history.certificate.revision_domain() ==
                    input.predictor.committed_face_flux_revision_domain &&
                (input.bdf.order == 2U
                     ? valid_committed_flux(
                           previous_history,
                           input.predictor.previous_face_flux) &&
                           accepted_history.revision !=
                               previous_history.revision &&
                           accepted_history.certificate.authority() ==
                               previous_history.certificate.authority() &&
                           accepted_history.certificate.storage() ==
                               previous_history.certificate.storage() &&
                           accepted_history.certificate.revision_domain() ==
                               previous_history.certificate.revision_domain()
                     : empty_flux(previous_history))
          : empty_flux(accepted_history) && empty_flux(previous_history);
  const bool same_density_lineage =
      input.density.storage_identity ==
          input.predictor.predicted_density_storage &&
      input.density.revision_domain ==
          input.predictor.predicted_density_revision_domain;
  const bool valid_density_authority =
      input.corrector != 1U ||
      (same_density_lineage &&
       ((impl.pressure_reference_kind == PressureReferenceKind::closed_mass &&
         input.density.revision == input.predictor.predicted_density) ||
        (impl.pressure_reference_kind ==
             PressureReferenceKind::boundary_absolute &&
         input.density.revision != input.predictor.predicted_density)));
  const BoundaryThermophysicalGhostPhase thermophysical_phase =
      input.corrector == 1U
          ? BoundaryThermophysicalGhostPhase::corrector_one
          : (input.corrector == 2U
                 ? BoundaryThermophysicalGhostPhase::corrector_two
                 : BoundaryThermophysicalGhostPhase::invalid);
  const BoundaryThermophysicalGhostContext thermophysical_context{
      input.momentum.time, impl.geometry->fingerprint(),
      input.pressure_reference.pressure_reference, input.numeric_boundary,
      thermophysical_phase};
  ThermophysicalBoundaryTokens thermophysical_boundary;
  const bool valid_thermophysical_boundary =
      thermophysical_boundary_tokens(
          *impl.boundary, thermophysical_context,
          input.thermophysical_boundary, as_const(input.density),
          impl.thermodynamics, impl.transport,
          thermophysical_boundary);
  const bool closed_reference =
      impl.pressure_reference_kind == PressureReferenceKind::closed_mass;
  const bool valid_c2_pressure_reference =
      input.corrector == 1U || !corrected_c1_authorization.valid() ||
      (closed_reference
           ? same_pressure_reference_certificate(
                 input.pressure_reference,
                 corrected_c1_authorization.output_pressure_reference)
           : same_pressure_reference_certificate(
                 input.pressure_reference,
                 corrected_c1_authorization.input_pressure_reference));
  const bool valid_thermophysical_lineage =
      input.corrector == 1U ||
      (closed_reference && corrected_c1_authorization.valid()
           ? thermophysical_boundary.semantics ==
                 corrected_c1_authorization
                     .thermophysical_boundary_semantics
           : thermophysical_boundary.collective_lineage ==
                     corrector_one_authorization
                         .thermophysical_boundary_collective_lineage &&
                 thermophysical_boundary.rank_local_lineage ==
                     corrector_one_authorization
                         .thermophysical_boundary_rank_local_lineage);
  const ConstFieldView thermophysical_pressure =
      input.thermophysical_boundary.binding.pressure_perturbation;
  const bool has_thermophysical_pressure =
      !empty_field_view(thermophysical_pressure);
  const bool valid_thermophysical_pressure =
      has_thermophysical_pressure &&
      thermophysical_pressure.field == impl.boundary->pressure_field() &&
      detail::valid_cell_view(thermophysical_pressure, cells, 0U, 1U, 1U);
  const bool valid_inputs =
      detail::valid_bdf_coefficients(input.bdf) &&
      input.bdf.order == input.predictor.order &&
      input.numeric_boundary != 0U &&
      (input.corrector == 1U || input.corrector == 2U) &&
      (input.corrector == 1U ? input.prior_corrector == 0U : true) &&
      input.density.field == impl.density_field &&
      detail::valid_cell_view(as_const(input.density), cells, 0U, 1U, 1U) &&
      detail::valid_cell_view(input.trial_velocity, cells, 0U, 3U, 0U) &&
      detail::valid_flux_view(input.trial_flux, cells,
                              input.trial_flux.revision) &&
      !input.trial_flux.certificate.valid() &&
      valid_temporal_reference &&
      valid_committed_face_history &&
      valid_density_authority &&
      valid_thermophysical_boundary &&
      valid_thermophysical_pressure &&
      valid_c2_pressure_reference &&
      valid_thermophysical_lineage &&
      (input.corrector == 1U
           ? input.trial_flux.revision == input.momentum.face_flux &&
                 input.trial_flux.revision ==
                     input.predictor.paired_face_flux &&
                 input.trial_flux.x.storage_identity ==
                     input.predictor.paired_face_flux_storage &&
                 input.trial_flux.x.revision_domain ==
                     input.predictor.paired_face_flux_revision_domain
           : input.trial_flux.revision !=
                 corrector_one_authorization.trial_face_flux) &&
      detail::valid_cell_view(input.momentum_system.diagonal, cells, 0U, 3U) &&
      detail::valid_cell_view(input.momentum_system.rhs, cells, 0U, 3U);
  const ConstFieldView cell_inputs[]{
      as_const(input.density), input.trial_velocity, thermophysical_pressure,
      as_const(input.momentum_system.diagonal),
      as_const(input.momentum_system.rhs)};
  bool aliases = false;
  const ConstFaceFieldView trial_faces[]{input.trial_flux.x,
                                         input.trial_flux.y,
                                         input.trial_flux.z};
  const ConstFaceFieldView output_faces[]{phi_output.x, phi_output.y,
                                          phi_output.z};
  const FaceFieldView pressure_faces[]{
      impl.workspace.x_pressure_coefficient,
      impl.workspace.y_pressure_coefficient,
      impl.workspace.z_pressure_coefficient};
  if (local && (!valid_certificates || !valid_inputs))
    local = {StatusCode::invalid_plan, kPisoCoupler};
  if (local) {
    for (ConstFieldView view : cell_inputs) {
      aliases = aliases ||
                detail::field_views_overlap(view,
                                            as_const(impl.workspace.r_au)) ||
                detail::field_views_overlap(view,
                                            as_const(impl.workspace.h_by_a));
      const FaceFieldView pressure_faces[]{
          impl.workspace.x_pressure_coefficient,
          impl.workspace.y_pressure_coefficient,
          impl.workspace.z_pressure_coefficient};
      for (FaceFieldView face : pressure_faces) {
        aliases = aliases || detail::cell_face_views_overlap(view, face);
      }
      const ConstFaceFluxView phi = as_const(impl.workspace.phi_h_by_a);
      aliases = aliases || detail::cell_face_views_overlap(view, phi.x) ||
                detail::cell_face_views_overlap(view, phi.y) ||
                detail::cell_face_views_overlap(view, phi.z);
    }
    for (ConstFaceFieldView trial_face : trial_faces) {
      for (ConstFaceFieldView output_face : output_faces) {
        aliases = aliases ||
                  detail::face_views_overlap(trial_face, output_face);
      }
      for (FaceFieldView pressure_face : pressure_faces) {
        aliases = aliases || detail::face_views_overlap(trial_face,
                                                        pressure_face);
      }
    }
    if (input.corrector == 1U) {
      const ConstFaceFieldView accepted_faces[]{accepted_history.x,
                                                accepted_history.y,
                                                accepted_history.z};
      const ConstFaceFieldView previous_faces[]{previous_history.x,
                                                previous_history.y,
                                                previous_history.z};
      const ConstFaceFieldView history_faces[]{
          accepted_faces[0U], accepted_faces[1U], accepted_faces[2U],
          previous_faces[0U], previous_faces[1U], previous_faces[2U]};
      const std::size_t history_count = input.bdf.order == 2U ? 6U : 3U;
      for (std::size_t history_index = 0U;
           history_index < history_count; ++history_index) {
        const ConstFaceFieldView history = history_faces[history_index];
        for (ConstFieldView cell : cell_inputs) {
          aliases = aliases ||
                    detail::cell_face_views_overlap(cell, history);
        }
        aliases = aliases || detail::cell_face_views_overlap(
                                 as_const(impl.workspace.r_au), history) ||
                  detail::cell_face_views_overlap(
                      as_const(impl.workspace.h_by_a), history) ||
                  detail::cell_face_views_overlap(
                      as_const(impl.workspace.pressure_gradient), history);
        for (ConstFaceFieldView output_face : output_faces) {
          aliases = aliases ||
                    detail::face_views_overlap(history, output_face);
        }
        for (ConstFaceFieldView trial_face : trial_faces) {
          aliases = aliases ||
                    detail::face_views_overlap(history, trial_face);
        }
        for (FaceFieldView pressure_face : pressure_faces) {
          aliases = aliases ||
                    detail::face_views_overlap(history, pressure_face);
        }
      }
      if (input.bdf.order == 2U) {
        for (ConstFaceFieldView accepted_face : accepted_faces) {
          for (ConstFaceFieldView previous_face : previous_faces) {
            aliases = aliases || detail::face_views_overlap(
                                     accepted_face, previous_face);
          }
        }
      }
    }
  }
  if (local && aliases)
    local = {StatusCode::invalid_plan, kPisoCoupler};

  const KernelBox box{{0, 0, 0}, cells};
  if (local &&
      (!detail::finite_field_box(as_const(input.density), box, 0U, 1U) ||
       !detail::finite_field_box(input.trial_velocity, box, 0U, 3U) ||
       !detail::finite_field_box(as_const(input.momentum_system.diagonal), box,
                                 0U, 3U) ||
       !detail::finite_field_box(as_const(input.momentum_system.rhs), box, 0U,
                                 3U))) {
    local = {StatusCode::numerical_failure, kPisoNumerical};
  }
  if (local) {
    for (std::int32_t z = -1; z <= cells.z && local; ++z)
      for (std::int32_t y = -1; y <= cells.y && local; ++y)
        for (std::int32_t x = -1; x <= cells.x; ++x) {
          const unsigned outside =
              static_cast<unsigned>(x < 0 || x >= cells.x) +
              static_cast<unsigned>(y < 0 || y >= cells.y) +
              static_cast<unsigned>(z < 0 || z >= cells.z);
          if (outside <= 1U &&
              !std::isfinite(
                  thermophysical_pressure.unchecked({x, y, z}, 0U))) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
        }
  }
  if (local) {
    for (std::int32_t z = 0; z < cells.z && local; ++z) {
      for (std::int32_t y = 0; y < cells.y && local; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (!finite_positive(input.density.unchecked(cell, 0U))) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            const double diagonal =
                input.momentum_system.diagonal.unchecked(cell, component);
            const double rhs =
                input.momentum_system.rhs.unchecked(cell, component);
            const double volume = detail::cell_volume(*impl.kernels, cell);
            const double r_au = volume / diagonal;
            const double h_by_a = rhs / diagonal;
            if (!finite_positive(diagonal) || !finite_positive(volume) ||
                !finite_positive(r_au) || !std::isfinite(h_by_a)) {
              local = {StatusCode::numerical_failure, kPisoNumerical};
              break;
            }
          }
        }
      }
    }
  }
  if (local) {
    for (std::int32_t z = 0; z < cells.z; ++z) {
      for (std::int32_t y = 0; y < cells.y; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const double volume = detail::cell_volume(*impl.kernels, cell);
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            const double diagonal =
                input.momentum_system.diagonal.unchecked(cell, component);
            impl.workspace.r_au.unchecked(cell, component) = volume / diagonal;
            // The momentum stage now publishes a true converged predictor.
            // C1 consumes that velocity directly; C2 remains incremental and
            // consumes the already pressure-corrected C1 velocity.  rAU alone
            // is reconstructed from the unchanged assembled diagonal.
            impl.workspace.h_by_a.unchecked(cell, component) =
                input.trial_velocity.unchecked(cell, component);
          }
        }
      }
    }
  }
  if (local) {
    const std::array<ConstFieldView, 1U> pressure_reads{
        thermophysical_pressure};
    const std::array<FieldView, 1U> gradient_writes{
        impl.workspace.pressure_gradient};
    const KernelInvocation gradient_invocation{
        {pressure_reads.data(), pressure_reads.size()},
        {gradient_writes.data(), gradient_writes.size()}, box, 0U, 0U, 1U,
        0U, nullptr};
    local = cartesian_gradient(*impl.kernels, gradient_invocation);
    if (local && input.immersed_interface != nullptr) {
      local = input.immersed_interface->correct_pressure_gradient(
          thermophysical_pressure, impl.workspace.pressure_gradient);
    }
  }
  std::array<FieldView, 4U> halo_views{
      input.density, impl.workspace.r_au, impl.workspace.h_by_a,
      impl.workspace.pressure_gradient};
  HaloTicket ticket;
  Status status = impl.halo->begin(
      impl.halo_stage, {halo_views.data(), halo_views.size()}, local, ticket);
  if (status) {
    status = impl.halo->finish(ticket,
                               {halo_views.data(), halo_views.size()});
  }
  if (!status) {
    return status;
  }
  const FieldView density_with_ghosts = halo_views[0U];
  impl.workspace.r_au = halo_views[1U];
  impl.workspace.h_by_a = halo_views[2U];
  impl.workspace.pressure_gradient = halo_views[3U];
  const RevisionToken density_ghost_revision =
      impl.halo->ghost_revision(density_with_ghosts.field);
  ThermophysicalBoundaryTokens revalidated_thermophysical_boundary;
  Status post_halo_local;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const bool injected_post_halo_failure =
      g_post_halo_revalidation_failure_rank.exchange(
          -1, std::memory_order_relaxed) == impl.rank;
#else
  constexpr bool injected_post_halo_failure = false;
#endif
  if (density_ghost_revision != density_with_ghosts.revision ||
      injected_post_halo_failure ||
      !thermophysical_boundary_tokens(
          *impl.boundary, thermophysical_context,
          input.thermophysical_boundary, as_const(density_with_ghosts),
          impl.thermodynamics, impl.transport,
          revalidated_thermophysical_boundary) ||
      revalidated_thermophysical_boundary.semantics !=
          thermophysical_boundary.semantics ||
      revalidated_thermophysical_boundary.collective_lineage !=
          thermophysical_boundary.collective_lineage ||
      revalidated_thermophysical_boundary.target !=
          thermophysical_boundary.target ||
      revalidated_thermophysical_boundary.rank_local_lineage !=
          thermophysical_boundary.rank_local_lineage ||
      revalidated_thermophysical_boundary.rank_local_binding !=
          thermophysical_boundary.rank_local_binding) {
    post_halo_local = {StatusCode::invalid_plan, kPisoCoupler};
  }
  int post_halo_lowest = -1;
  const Status post_halo_consensus =
      collective_status(impl.communicator, post_halo_local, impl.rank,
                        impl.size, post_halo_lowest);
  if (!post_halo_consensus) return post_halo_consensus;
  // MPI exchange owns inter-rank ghosts only.  Physical rho ghosts arrive
  // already closed by the thermophysical boundary authority from the same
  // p/h/Y revision.  Replacing them with owner rho here breaks EOS closure
  // and changes the pressure-face coefficient at Dirichlet thermodynamic
  // boundaries.  rAU is algebraic and retains its zero-normal-gradient
  // closure below.
  fill_physical_zero_gradient(impl.workspace.r_au, cells, 3U, 1U,
                              *impl.boundary);
  fill_physical_zero_gradient(impl.workspace.pressure_gradient, cells, 3U,
                              1U, *impl.boundary);
  FieldView constrained_h_by_a = impl.workspace.h_by_a;
  constrained_h_by_a.field = impl.boundary->velocity_field();
  status = apply_boundary_ghosts(
      BoundaryStage::momentum, *impl.boundary,
      {&constrained_h_by_a, 1U}, input.boundary_values);
  if (!status) {
    return status;
  }
  status = fill_pressure_coefficients<CartesianAxis::x>(
      *impl.kernels, as_const(density_with_ghosts),
      as_const(impl.workspace.r_au),
      *impl.geometry, impl.patch, impl.pressure_boundary,
      impl.workspace.x_pressure_coefficient);
  if (status) {
    status = fill_pressure_coefficients<CartesianAxis::y>(
        *impl.kernels, as_const(density_with_ghosts),
        as_const(impl.workspace.r_au),
        *impl.geometry, impl.patch, impl.pressure_boundary,
        impl.workspace.y_pressure_coefficient);
  }
  if (status) {
    status = fill_pressure_coefficients<CartesianAxis::z>(
        *impl.kernels, as_const(density_with_ghosts),
        as_const(impl.workspace.r_au),
        *impl.geometry, impl.patch, impl.pressure_boundary,
        impl.workspace.z_pressure_coefficient);
  }
  if (!status) {
    return status;
  }
  if (input.corrector == 1U) {
    const ConstFaceFieldView temporal_faces[]{
        input.temporal_reference.x, input.temporal_reference.y,
        input.temporal_reference.z};
    const ConstFaceFieldView paired_faces[]{input.trial_flux.x,
                                            input.trial_flux.y,
                                            input.trial_flux.z};
    const ConstFaceFieldView accepted_faces[]{accepted_history.x,
                                              accepted_history.y,
                                              accepted_history.z};
    const ConstFaceFieldView previous_faces[]{previous_history.x,
                                              previous_history.y,
                                              previous_history.z};
    const FaceFieldView output_faces[]{impl.workspace.phi_h_by_a.x,
                                       impl.workspace.phi_h_by_a.y,
                                       impl.workspace.phi_h_by_a.z};
    const FaceFieldView auxiliary_faces[]{impl.frozen_candidate_face_aux.x,
                                          impl.frozen_candidate_face_aux.y,
                                          impl.frozen_candidate_face_aux.z};
    const ConstFaceFieldView coefficient_faces[]{
        as_const(impl.workspace.x_pressure_coefficient),
        as_const(impl.workspace.y_pressure_coefficient),
        as_const(impl.workspace.z_pressure_coefficient)};
    for (std::size_t axis_index = 0U; axis_index < 3U && status;
         ++axis_index) {
      const auto axis = static_cast<CartesianAxis>(axis_index);
      const FaceFieldView output = output_faces[axis_index];
      for (std::int32_t z = 0; z < output.extents.z && status; ++z) {
        for (std::int32_t y = 0; y < output.extents.y && status; ++y) {
          for (std::int32_t x = 0; x < output.extents.x; ++x) {
            const Int3 face{x, y, z};
            const std::int32_t normal = axis_value(face, axis);
            const std::int32_t extent = axis_value(cells, axis);
            const bool boundary_face = normal == 0 || normal == extent;
            bool fixed_flux = false;
            if (boundary_face) {
              PressureCorrectionFaceRule rule;
              status =
                  impl.pressure_boundary.face_rule_unchecked(axis, face, rule);
              if (!status) {
                break;
              }
              fixed_flux =
                  rule.kind ==
                  PressureCorrectionFaceKind::homogeneous_neumann;
            }
            if (fixed_flux) {
              const double paired =
                  paired_faces[axis_index].unchecked(face);
              if (!std::isfinite(paired)) {
                status = {StatusCode::numerical_failure, kPisoNumerical};
                break;
              }
              paired_flux_consumed = true;
              output.unchecked(face) = paired;
              auxiliary_faces[axis_index].unchecked(face) = 0.0;
              continue;
            }
            Int3 left = face;
            if (axis == CartesianAxis::x)
              --left.x;
            else if (axis == CartesianAxis::y)
              --left.y;
            else
              --left.z;
            const std::uint8_t component =
                static_cast<std::uint8_t>(axis_index);
            const double left_rho = density_with_ghosts.unchecked(left, 0U);
            const double right_rho = density_with_ghosts.unchecked(face, 0U);
            const double current =
                detail::interpolate_face(
                    *impl.kernels, axis, normal,
                    left_rho * impl.workspace.h_by_a.unchecked(left,
                                                               component),
                    right_rho * impl.workspace.h_by_a.unchecked(face,
                                                                component)) *
                detail::face_area(*impl.kernels, axis, face);
            const double rho_r_au = detail::interpolate_face(
                *impl.kernels, axis, normal,
                left_rho * impl.workspace.r_au.unchecked(left, component),
                right_rho * impl.workspace.r_au.unchecked(face, component));
            const double weight = input.bdf.a0 * rho_r_au;
            const double temporal =
                temporal_faces[axis_index].unchecked(face);
            const double accepted =
                accepted_faces[axis_index].unchecked(face);
            const double previous =
                input.bdf.order == 2U
                    ? previous_faces[axis_index].unchecked(face)
                    : 0.0;
            const double normalized_history =
                input.bdf.order == 2U
                    ? (-input.bdf.a1 * accepted -
                       input.bdf.a2 * previous) /
                          input.bdf.a0
                    : accepted;
            const double pressure_gradient_flux =
                detail::interpolate_face(
                    *impl.kernels, axis, normal,
                    left_rho *
                        impl.workspace.r_au.unchecked(left, component) *
                        impl.workspace.pressure_gradient.unchecked(
                            left, component),
                    right_rho *
                        impl.workspace.r_au.unchecked(face, component) *
                        impl.workspace.pressure_gradient.unchecked(
                            face, component)) *
                detail::face_area(*impl.kernels, axis, face);
            const double pressure_jump_flux =
                impl.pressure_boundary.mass_flux_response_unchecked(
                    thermophysical_pressure, axis, face,
                    coefficient_faces[axis_index].unchecked(face));
            const double value =
                current + weight * (normalized_history - temporal) +
                pressure_gradient_flux + pressure_jump_flux;
            if (!std::isfinite(current) || !std::isfinite(rho_r_au) ||
                !std::isfinite(weight) || weight < 0.0 ||
                !std::isfinite(temporal) ||
                !std::isfinite(accepted) || !std::isfinite(previous) ||
                !std::isfinite(normalized_history) ||
                !std::isfinite(pressure_gradient_flux) ||
                !std::isfinite(pressure_jump_flux) ||
                !std::isfinite(value)) {
              status = {StatusCode::numerical_failure, kPisoNumerical};
              break;
            }
            output.unchecked(face) = value;
            // C1 candidate replay keeps only the BDF/history velocity offset
            // theta frozen.  rho and rho*rAU are refreshed for every alpha.
            auxiliary_faces[axis_index].unchecked(face) =
                normalized_history - temporal;
          }
        }
      }
    }
    if (status) {
      RevisionToken corrected_revision = kFnvOffset;
      corrected_revision =
          hash_mix(corrected_revision, input.temporal_reference.revision);
      corrected_revision =
          hash_mix(corrected_revision, accepted_history.revision);
      corrected_revision =
          hash_mix(corrected_revision, previous_history.revision);
      if (paired_flux_consumed) {
        corrected_revision =
            hash_mix(corrected_revision, input.trial_flux.revision);
      }
      corrected_revision = hash_mix(corrected_revision, input.momentum.state);
      corrected_revision =
          hash_mix(corrected_revision, thermophysical_pressure.revision);
      impl.workspace.phi_h_by_a.revision =
          corrected_revision == 0U ? RevisionToken{1U}
                                   : corrected_revision;
    }
  } else {
    // Corrector two must use the already corrected C1 face flux as its
    // incremental base.  Keep the existing workspace/revision authority and
    // only copy/check its preallocated directional payload.
    const ConstFaceFieldView input_faces[]{
        input.trial_flux.x, input.trial_flux.y, input.trial_flux.z};
    const FaceFieldView output_faces[]{
        impl.workspace.phi_h_by_a.x, impl.workspace.phi_h_by_a.y,
        impl.workspace.phi_h_by_a.z};
    impl.workspace.phi_h_by_a.revision = input.trial_flux.revision;
    for (std::uint8_t axis = 0U; axis < 3U && status; ++axis) {
      const Int3 extents = output_faces[axis].extents;
      for (std::int32_t z = 0; z < extents.z && status; ++z) {
        for (std::int32_t y = 0; y < extents.y && status; ++y) {
          for (std::int32_t x = 0; x < extents.x; ++x) {
            const double value = input_faces[axis].unchecked({x, y, z});
            if (!std::isfinite(value)) {
              status = {StatusCode::numerical_failure, kPisoNumerical};
              break;
            }
            output_faces[axis].unchecked({x, y, z}) = value;
          }
        }
      }
    }
    // C2 starts from the already corrected C1 total flux, but its candidate
    // density response must not be lost by a byte-copy base.  Freeze the
    // explicit offset Phi0-Q(rho0), with Q=A I(rho HbyA), so replay can
    // refresh Q(rho_alpha) while preserving the exact C1 base at alpha zero.
    const FaceFieldView auxiliary_faces[]{impl.frozen_candidate_face_aux.x,
                                          impl.frozen_candidate_face_aux.y,
                                          impl.frozen_candidate_face_aux.z};
    for (std::uint8_t axis_index = 0U; axis_index < 3U && status;
         ++axis_index) {
      const auto axis = static_cast<CartesianAxis>(axis_index);
      const ConstFaceFieldView base = input_faces[axis_index];
      const FaceFieldView auxiliary = auxiliary_faces[axis_index];
      for (std::int32_t z = 0; z < auxiliary.extents.z && status; ++z) {
        for (std::int32_t y = 0; y < auxiliary.extents.y && status; ++y) {
          for (std::int32_t x = 0; x < auxiliary.extents.x; ++x) {
            const Int3 face{x, y, z};
            Int3 left = face;
            if (axis == CartesianAxis::x)
              --left.x;
            else if (axis == CartesianAxis::y)
              --left.y;
            else
              --left.z;
            const std::uint8_t component = axis_index;
            const double q =
                detail::interpolate_face(
                    *impl.kernels, axis, axis_value(face, axis),
                    density_with_ghosts.unchecked(left, 0U) *
                        impl.workspace.h_by_a.unchecked(left, component),
                    density_with_ghosts.unchecked(face, 0U) *
                        impl.workspace.h_by_a.unchecked(face, component)) *
                detail::face_area(*impl.kernels, axis, face);
            const double offset = base.unchecked(face) - q;
            if (!std::isfinite(q) || !std::isfinite(offset)) {
              status = {StatusCode::numerical_failure, kPisoNumerical};
              break;
            }
            auxiliary.unchecked(face) = offset;
          }
        }
      }
    }
  }
  if (!status) {
    return status;
  }
  if (input.immersed_interface != nullptr) {
    status = input.immersed_interface->constrain_pressure_predictor(
        impl.workspace.h_by_a, impl.workspace.phi_h_by_a);
    if (!status) return status;
  }

  std::uint64_t r_au_dependency = kFnvOffset;
  r_au_dependency = hash_mix(r_au_dependency, input.momentum.state);
  r_au_dependency = mix_view(
      r_au_dependency, as_const(input.momentum_system.diagonal));
  // Density is intentionally explicit rather than assumed to be represented
  // by the diagonal certificate: source/constraint implementations may keep
  // the same matrix storage while changing its thermophysical authority.
  r_au_dependency = mix_view(r_au_dependency, as_const(input.density));
  r_au_dependency = hash_mix(r_au_dependency, double_bits(input.bdf.a0));
  r_au_dependency = hash_mix(r_au_dependency, input.bdf.order);

  std::uint64_t h_by_a_dependency = kFnvOffset;
  h_by_a_dependency = hash_mix(h_by_a_dependency, r_au_dependency);
  h_by_a_dependency = mix_view(
      h_by_a_dependency, as_const(input.momentum_system.rhs));
  h_by_a_dependency = mix_view(h_by_a_dependency, input.trial_velocity);

  std::uint64_t pressure_face_dependency = kFnvOffset;
  pressure_face_dependency =
      hash_mix(pressure_face_dependency, r_au_dependency);
  pressure_face_dependency =
      mix_view(pressure_face_dependency, as_const(input.density));
  pressure_face_dependency = hash_mix(
      pressure_face_dependency,
      input.pressure_reference.pressure_reference);
  pressure_face_dependency =
      hash_mix(pressure_face_dependency, input.numeric_boundary);
  pressure_face_dependency =
      hash_mix(pressure_face_dependency, density_ghost_revision);
  pressure_face_dependency = hash_mix(
      pressure_face_dependency, thermophysical_boundary.semantics);
  pressure_face_dependency = hash_mix(
      pressure_face_dependency, thermophysical_boundary.collective_lineage);
  pressure_face_dependency =
      hash_mix(pressure_face_dependency, thermophysical_boundary.target);
  pressure_face_dependency = hash_mix(
      pressure_face_dependency, thermophysical_boundary.rank_local_binding);
  pressure_face_dependency = hash_mix(
      pressure_face_dependency, thermophysical_boundary.rank_local_lineage);

  RevisionToken committed_face_history_dependency = 0U;
  if (input.corrector == 1U) {
    const auto mix_committed_flux = [&](std::uint64_t hash,
                                        ConstFaceFluxView flux) noexcept {
      hash = mix_face(hash, flux.x);
      hash = mix_face(hash, flux.y);
      hash = mix_face(hash, flux.z);
      hash = hash_mix(hash, flux.revision);
      hash = hash_mix(hash, flux.certificate.revision());
      hash = hash_mix(hash, flux.certificate.authority());
      hash = hash_mix(hash, flux.certificate.storage());
      return hash_mix(hash, flux.certificate.revision_domain());
    };
    std::uint64_t history_dependency = kFnvOffset;
    history_dependency =
        mix_committed_flux(history_dependency, accepted_history);
    if (input.bdf.order == 2U) {
      history_dependency =
          mix_committed_flux(history_dependency, previous_history);
    }
    committed_face_history_dependency =
        history_dependency == 0U ? RevisionToken{1U} : history_dependency;
  }

  std::uint64_t phi_h_by_a_dependency = kFnvOffset;
  phi_h_by_a_dependency =
      hash_mix(phi_h_by_a_dependency, h_by_a_dependency);
  phi_h_by_a_dependency =
      hash_mix(phi_h_by_a_dependency, pressure_face_dependency);
  phi_h_by_a_dependency =
      hash_mix(phi_h_by_a_dependency, input.predictor.state);
  phi_h_by_a_dependency =
      hash_mix(phi_h_by_a_dependency, input.momentum.face_flux);
  if (input.corrector == 2U || paired_flux_consumed) {
    phi_h_by_a_dependency =
        mix_face(phi_h_by_a_dependency, input.trial_flux.x);
    phi_h_by_a_dependency =
        mix_face(phi_h_by_a_dependency, input.trial_flux.y);
    phi_h_by_a_dependency =
        mix_face(phi_h_by_a_dependency, input.trial_flux.z);
    phi_h_by_a_dependency =
        hash_mix(phi_h_by_a_dependency, input.trial_flux.revision);
  }
  if (input.corrector == 1U) {
    phi_h_by_a_dependency =
        mix_view(phi_h_by_a_dependency, thermophysical_pressure);
    phi_h_by_a_dependency = hash_mix(
        phi_h_by_a_dependency, committed_face_history_dependency);
    phi_h_by_a_dependency = mix_face(phi_h_by_a_dependency,
                                     input.temporal_reference.x);
    phi_h_by_a_dependency = mix_face(phi_h_by_a_dependency,
                                     input.temporal_reference.y);
    phi_h_by_a_dependency = mix_face(phi_h_by_a_dependency,
                                     input.temporal_reference.z);
    phi_h_by_a_dependency = hash_mix(
        phi_h_by_a_dependency, input.temporal_reference.revision);
  }
  phi_h_by_a_dependency =
      hash_mix(phi_h_by_a_dependency, input.numeric_boundary);
  phi_h_by_a_dependency =
      hash_mix(phi_h_by_a_dependency, double_bits(input.bdf.a0));
  phi_h_by_a_dependency =
      hash_mix(phi_h_by_a_dependency, double_bits(input.bdf.a1));
  phi_h_by_a_dependency =
      hash_mix(phi_h_by_a_dependency, double_bits(input.bdf.a2));
  phi_h_by_a_dependency =
      hash_mix(phi_h_by_a_dependency, input.bdf.order);
  if (input.immersed_interface != nullptr)
    phi_h_by_a_dependency = hash_mix(
        phi_h_by_a_dependency, input.immersed_interface->fingerprint());

  std::uint64_t dependency = kFnvOffset;
  dependency = hash_mix(dependency, phi_h_by_a_dependency);
  dependency = hash_mix(dependency, input.prior_corrector);
  dependency = hash_mix(dependency, input.corrector);
  dependency = hash_mix(dependency, thermophysical_boundary.semantics);
  dependency =
      hash_mix(dependency, thermophysical_boundary.collective_lineage);
  dependency = hash_mix(dependency, thermophysical_boundary.target);
  dependency =
      hash_mix(dependency, thermophysical_boundary.rank_local_binding);
  dependency =
      hash_mix(dependency, thermophysical_boundary.rank_local_lineage);
  dependency = hash_mix(dependency, refinement_iteration);
  dependency = hash_mix(dependency, refinement_collective_lineage);
  dependency = hash_mix(dependency, refinement_lineage);
  dependency = dependency == 0U ? 1U : dependency;
  PisoIntermediateCertificate candidate;
  candidate.plan = impl.fingerprint;
  candidate.r_au =
      hash_mix(r_au_dependency, impl.workspace.r_au.revision);
  candidate.h_by_a =
      hash_mix(h_by_a_dependency, impl.workspace.h_by_a.revision);
  candidate.pressure_face_coefficient = hash_mix(
      pressure_face_dependency,
      impl.workspace.x_pressure_coefficient.storage_identity ^
          impl.workspace.y_pressure_coefficient.storage_identity ^
          impl.workspace.z_pressure_coefficient.storage_identity);
  candidate.phi_h_by_a = hash_mix(
      phi_h_by_a_dependency, impl.workspace.phi_h_by_a.revision);
  candidate.trial_face_flux = input.trial_flux.revision;
  candidate.temporal_face_flux =
      input.corrector == 1U ? input.temporal_reference.revision : 0U;
  candidate.committed_face_history =
      input.corrector == 1U ? committed_face_history_dependency : 0U;
  candidate.dependency = dependency;
  candidate.corrector = input.corrector;
  candidate.thermophysical_boundary_semantics =
      thermophysical_boundary.semantics;
  candidate.thermophysical_boundary_target = thermophysical_boundary.target;
  candidate.thermophysical_boundary_rank_local_binding =
      thermophysical_boundary.rank_local_binding;
  candidate.thermophysical_boundary_collective_lineage =
      thermophysical_boundary.collective_lineage;
  candidate.thermophysical_boundary_rank_local_lineage =
      thermophysical_boundary.rank_local_lineage;
  candidate.pressure_energy_refinement = refinement_iteration;
  candidate.pressure_energy_refinement_collective_lineage =
      refinement_collective_lineage;
  candidate.pressure_energy_refinement_lineage = refinement_lineage;
  if (input.corrector == 1U) {
    impl.corrector_one = candidate;
  } else {
    impl.terminal_lineage_source = candidate;
    impl.current_corrected_c1 = corrected_c1_authorization;
  }
  impl.current = candidate;
  impl.current_trial_velocity = input.trial_velocity;
  if (has_thermophysical_pressure) {
    impl.current_pressure_perturbation = thermophysical_pressure;
    impl.current_pressure_perturbation_numeric_fingerprint =
        scalar_reach_one_numeric_fingerprint(thermophysical_pressure, cells);
  }
  impl.current_thermophysical_context = thermophysical_context;
  impl.current_absolute_pressure_reference =
      input.thermophysical_boundary.binding.pressure_reference;
  for (std::int32_t z = -1; z <= cells.z; ++z)
    for (std::int32_t y = -1; y <= cells.y; ++y)
      for (std::int32_t x = -1; x <= cells.x; ++x) {
        const unsigned outside =
            static_cast<unsigned>(x < 0 || x >= cells.x) +
            static_cast<unsigned>(y < 0 || y >= cells.y) +
            static_cast<unsigned>(z < 0 || z >= cells.z);
        if (outside > 1U) continue;
        impl.frozen_candidate_density.unchecked({x, y, z}, 0U) =
            density_with_ghosts.unchecked({x, y, z}, 0U);
      }
  RevisionToken frozen_revision = hash_mix(
      hash_mix(candidate.dependency, input.density.revision),
      density_ghost_revision);
  if (frozen_revision == 0U) frozen_revision = 1U;
  impl.frozen_candidate_density.revision = frozen_revision;
  impl.frozen_candidate_face_aux.revision =
      hash_mix(frozen_revision, candidate.phi_h_by_a);
  if (impl.frozen_candidate_face_aux.revision == 0U)
    impl.frozen_candidate_face_aux.revision = 1U;
  impl.frozen_candidate_bdf = input.bdf;
  impl.frozen_candidate_baseline =
      hash_mix(candidate.dependency, UINT64_C(0x66726f7a656e6261));
  if (impl.frozen_candidate_baseline == 0U)
    impl.frozen_candidate_baseline = 1U;
  const bool empty_activity =
      impl.continuity_activity.cells.data == nullptr &&
      impl.continuity_activity.cells.size == 0U &&
      impl.continuity_activity.x_faces.data == nullptr &&
      impl.continuity_activity.x_faces.size == 0U &&
      impl.continuity_activity.y_faces.data == nullptr &&
      impl.continuity_activity.y_faces.size == 0U &&
      impl.continuity_activity.z_faces.data == nullptr &&
      impl.continuity_activity.z_faces.size == 0U &&
      impl.continuity_activity.local_fingerprint == 0U &&
      impl.continuity_activity.collective_fingerprint == 0U;
  const bool complete_activity =
      impl.continuity_activity.cells.data != nullptr &&
      impl.continuity_activity.cells.size != 0U &&
      impl.continuity_activity.x_faces.data != nullptr &&
      impl.continuity_activity.x_faces.size != 0U &&
      impl.continuity_activity.y_faces.data != nullptr &&
      impl.continuity_activity.y_faces.size != 0U &&
      impl.continuity_activity.z_faces.data != nullptr &&
      impl.continuity_activity.z_faces.size != 0U &&
      impl.continuity_activity.local_fingerprint != 0U &&
      impl.continuity_activity.collective_fingerprint != 0U;
  const bool periodic_candidate_scope =
      entirely_periodic(*impl.boundary) &&
      impl.pressure_reference_kind == PressureReferenceKind::closed_mass &&
      ((input.immersed_interface == nullptr &&
        impl.pressure_correction_donors == nullptr &&
        impl.candidate_pressure_correction_donors == nullptr &&
        empty_activity) ||
       (input.immersed_interface != nullptr &&
        impl.pressure_correction_donors != nullptr &&
        impl.candidate_pressure_correction_donors != nullptr &&
        impl.candidate_pressure_correction_donor_fingerprint ==
            impl.candidate_pressure_correction_donors->fingerprint() &&
        complete_activity));
  const bool open_candidate_scope =
      !entirely_periodic(*impl.boundary) &&
      impl.pressure_reference_kind ==
          PressureReferenceKind::boundary_absolute &&
      input.immersed_interface == nullptr &&
      impl.pressure_correction_donors == nullptr &&
      impl.candidate_pressure_correction_donors == nullptr && empty_activity;
  const bool open_ibm_candidate_scope =
      !entirely_periodic(*impl.boundary) &&
      impl.pressure_reference_kind ==
          PressureReferenceKind::boundary_absolute &&
      input.immersed_interface != nullptr &&
      impl.pressure_correction_donors != nullptr &&
      impl.candidate_pressure_correction_donors != nullptr &&
      impl.candidate_pressure_correction_donor_fingerprint ==
          impl.candidate_pressure_correction_donors->fingerprint() &&
      complete_activity;
  impl.frozen_candidate_cartesian =
      periodic_candidate_scope || open_candidate_scope ||
      open_ibm_candidate_scope;
  impl.frozen_candidate_numeric =
      impl.frozen_candidate_numeric_fingerprint();
  if (impl.frozen_candidate_numeric == 0U) {
    impl.frozen_candidate_baseline = 0U;
    impl.frozen_candidate_cartesian = false;
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  impl.frozen_candidate_baseline = hash_mix(
      impl.frozen_candidate_baseline, impl.frozen_candidate_numeric);
  if (impl.frozen_candidate_baseline == 0U)
    impl.frozen_candidate_baseline = 1U;
  certificate = candidate;
  return {};
}

Status PressureVelocityCoupler::inspect_intermediate_flux(
    const PisoIntermediateCertificate& intermediate,
    ConstFaceFluxView& flux) const noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  const PisoIntermediateCertificate& current = impl.current;
  const bool same_current =
      intermediate.valid() && current.valid() &&
      intermediate.plan == current.plan &&
      intermediate.r_au == current.r_au &&
      intermediate.h_by_a == current.h_by_a &&
      intermediate.pressure_face_coefficient ==
          current.pressure_face_coefficient &&
      intermediate.phi_h_by_a == current.phi_h_by_a &&
      intermediate.trial_face_flux == current.trial_face_flux &&
      intermediate.temporal_face_flux == current.temporal_face_flux &&
      intermediate.committed_face_history ==
          current.committed_face_history &&
      intermediate.dependency == current.dependency &&
      intermediate.corrector == current.corrector &&
      intermediate.thermophysical_boundary_semantics ==
          current.thermophysical_boundary_semantics &&
      intermediate.thermophysical_boundary_target ==
          current.thermophysical_boundary_target &&
      intermediate.thermophysical_boundary_rank_local_binding ==
          current.thermophysical_boundary_rank_local_binding &&
      intermediate.thermophysical_boundary_collective_lineage ==
          current.thermophysical_boundary_collective_lineage &&
      intermediate.thermophysical_boundary_rank_local_lineage ==
          current.thermophysical_boundary_rank_local_lineage;
  ConstFaceFluxView candidate = as_const(impl.workspace.phi_h_by_a);
  if (!same_current ||
      !detail::valid_flux_view(candidate, impl.cells, candidate.revision) ||
      candidate.certificate.valid()) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  candidate.certificate = {};
  flux = candidate;
  return {};
}

Status PressureVelocityCoupler::inspect_cartesian_pressure_work_linearization(
    const PisoIntermediateCertificate& intermediate,
    const PressureCorrectionCertificate& pressure,
    ConstFieldView target_pressure_perturbation,
    ConstFieldView target_velocity,
    PisoCartesianPressureWorkLinearization& linearization) const noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  const bool empty_activity =
      impl.continuity_activity.cells.data == nullptr &&
      impl.continuity_activity.cells.size == 0U &&
      impl.continuity_activity.x_faces.data == nullptr &&
      impl.continuity_activity.x_faces.size == 0U &&
      impl.continuity_activity.y_faces.data == nullptr &&
      impl.continuity_activity.y_faces.size == 0U &&
      impl.continuity_activity.z_faces.data == nullptr &&
      impl.continuity_activity.z_faces.size == 0U &&
      impl.continuity_activity.local_fingerprint == 0U &&
      impl.continuity_activity.collective_fingerprint == 0U;
  const bool identities =
      impl.immersed_interface == nullptr &&
      impl.pressure_correction_donors == nullptr && empty_activity &&
      impl.current.valid() && impl.current_pressure_work.valid() &&
      same_intermediate_certificate(intermediate, impl.current) &&
      same_pressure_certificate(pressure, impl.current_pressure_work) &&
      pressure.intermediate == intermediate.dependency &&
      same_field_identity(target_pressure_perturbation,
                          impl.current_pressure_perturbation) &&
      same_field_identity(target_velocity, impl.current_trial_velocity) &&
      same_cells(target_pressure_perturbation.interior, impl.cells) &&
      target_pressure_perturbation.field == impl.boundary->pressure_field() &&
      target_velocity.field == impl.boundary->velocity_field() &&
      detail::valid_cell_view(target_pressure_perturbation, impl.cells, 0U,
                              1U, 1U) &&
      detail::valid_cell_view(target_velocity, impl.cells, 0U, 3U, 0U) &&
      detail::valid_cell_view(as_const(impl.workspace.h_by_a), impl.cells, 0U,
                              3U, 1U) &&
      detail::valid_cell_view(as_const(impl.workspace.r_au), impl.cells, 0U,
                              3U, 1U) &&
      impl.current_pressure_perturbation_numeric_fingerprint != 0U &&
      scalar_reach_one_numeric_fingerprint(target_pressure_perturbation,
                                           impl.cells) ==
          impl.current_pressure_perturbation_numeric_fingerprint;
  bool admissible = identities;
  if (admissible) {
    const ConstFieldView h_by_a = as_const(impl.workspace.h_by_a);
    const ConstFieldView r_au = as_const(impl.workspace.r_au);
    admissible =
        !detail::field_views_overlap(target_pressure_perturbation,
                                     target_velocity) &&
        !detail::field_views_overlap(target_pressure_perturbation, h_by_a) &&
        !detail::field_views_overlap(target_pressure_perturbation, r_au) &&
        !detail::field_views_overlap(target_velocity, h_by_a) &&
        !detail::field_views_overlap(target_velocity, r_au) &&
        !detail::field_views_overlap(h_by_a, r_au);
    for (std::int32_t z = -1; z <= impl.cells.z && admissible; ++z) {
      for (std::int32_t y = -1; y <= impl.cells.y && admissible; ++y) {
        for (std::int32_t x = -1; x <= impl.cells.x; ++x) {
          if (!std::isfinite(
                  target_pressure_perturbation.unchecked({x, y, z}, 0U))) {
            admissible = false;
            break;
          }
        }
      }
    }
    for (std::int32_t z = 0; z < impl.cells.z && admissible; ++z) {
      for (std::int32_t y = 0; y < impl.cells.y && admissible; ++y) {
        for (std::int32_t x = 0; x < impl.cells.x; ++x) {
          const Int3 cell{x, y, z};
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            const double velocity =
                target_velocity.unchecked(cell, component);
            const double predictor = h_by_a.unchecked(cell, component);
            const double reciprocal = r_au.unchecked(cell, component);
            if (!std::isfinite(velocity) || !std::isfinite(predictor) ||
                !finite_positive(reciprocal) ||
                double_bits(velocity) != double_bits(predictor)) {
              admissible = false;
              break;
            }
          }
          if (!admissible) break;
        }
      }
    }
  }
  if (!admissible) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }

  PisoCartesianPressureWorkLinearization candidate;
  candidate.issuer_ = this;
  candidate.kernels_ = impl.kernels;
  candidate.geometry_ = impl.geometry;
  candidate.patch_ = impl.patch;
  candidate.intermediate_ = intermediate;
  candidate.pressure_ = pressure;
  candidate.target_pressure_perturbation_ = target_pressure_perturbation;
  candidate.target_velocity_ = target_velocity;
  candidate.h_by_a_ = as_const(impl.workspace.h_by_a);
  candidate.r_au_ = as_const(impl.workspace.r_au);
  candidate.numeric_fingerprint_ =
      cartesian_pressure_work_numeric_fingerprint(
          candidate.target_pressure_perturbation_, candidate.target_velocity_,
          candidate.h_by_a_, candidate.r_au_, impl.cells, intermediate,
          pressure, impl.kernels->fingerprint());
  std::uint64_t authority = hash_mix(
      kFnvOffset, UINT64_C(0x7069736f6570776c));
  authority = hash_mix(authority, impl.fingerprint);
  authority = hash_mix(authority, intermediate.dependency);
  authority = hash_mix(authority, pressure.state);
  authority = hash_mix(authority, candidate.numeric_fingerprint_);
  authority = hash_mix(authority,
                       reinterpret_cast<std::uintptr_t>(this));
  candidate.authority_ = authority == 0U ? RevisionToken{1U} : authority;
  if (!candidate.valid() ||
      !matches_cartesian_pressure_work_linearization(candidate, true)) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  linearization = candidate;
  return {};
}

bool PressureVelocityCoupler::matches_cartesian_pressure_work_linearization(
    const PisoCartesianPressureWorkLinearization& linearization,
    bool verify_numeric_fingerprint) const noexcept {
  if (implementation_ == nullptr || !linearization.valid() ||
      linearization.issuer_ != this) {
    return false;
  }
  const Impl& impl = *implementation_;
  const bool current =
      linearization.kernels_ == impl.kernels &&
      linearization.geometry_ == impl.geometry &&
      same_patch(linearization.patch_, impl.patch) &&
      impl.current.valid() && impl.current_pressure_work.valid() &&
      same_intermediate_certificate(linearization.intermediate_,
                                    impl.current) &&
      same_pressure_certificate(linearization.pressure_,
                                impl.current_pressure_work) &&
      same_field_identity(linearization.target_pressure_perturbation_,
                          impl.current_pressure_perturbation) &&
      same_field_identity(linearization.target_velocity_,
                          impl.current_trial_velocity) &&
      same_field_identity(linearization.h_by_a_,
                          as_const(impl.workspace.h_by_a)) &&
      same_field_identity(linearization.r_au_,
                          as_const(impl.workspace.r_au)) &&
      linearization.numeric_fingerprint_ != 0U &&
      linearization.authority_ != 0U;
  return current &&
         (!verify_numeric_fingerprint ||
          cartesian_pressure_work_numeric_fingerprint(
              linearization.target_pressure_perturbation_,
              linearization.target_velocity_, linearization.h_by_a_,
              linearization.r_au_, impl.cells, linearization.intermediate_,
              linearization.pressure_, impl.kernels->fingerprint()) ==
              linearization.numeric_fingerprint_);
}

bool PressureVelocityCoupler::matches_candidate_boundary_finalizer_binding(
    const CartesianGeometryPlan* geometry, MeshPatch patch,
    const BoundaryPlan* boundary, const CartesianKernelPlan* kernels,
    const ThermodynamicsPlan* thermodynamics,
    const TransportPlan* transport,
    const IbmEquationInterfacePlan* immersed_interface,
    const RemoteDonorExchangePlan* candidate_donors,
    StageId candidate_donor_stage, FieldId candidate_field,
    std::uint8_t candidate_reach) const noexcept {
  if (implementation_ == nullptr || geometry == nullptr ||
      boundary == nullptr || kernels == nullptr || thermodynamics == nullptr ||
      transport == nullptr) {
    return false;
  }
  const Impl& impl = *implementation_;
  return geometry == impl.geometry && boundary == impl.boundary &&
         kernels == impl.kernels && same_patch(patch, impl.patch) &&
         thermodynamics->fingerprint() == impl.thermodynamics &&
         transport->fingerprint() == impl.transport &&
         immersed_interface == impl.immersed_interface &&
         candidate_donors ==
             impl.candidate_pressure_correction_donors &&
         candidate_donor_stage ==
             impl.candidate_pressure_correction_donor_stage &&
         candidate_field == impl.candidate_pressure_correction_field &&
         candidate_reach ==
             impl.candidate_pressure_correction_donor_reach &&
         impl.pressure_reference_kind ==
             PressureReferenceKind::boundary_absolute;
}

Status PressureVelocityCoupler::assemble_pressure_system(
    const PressureCorrectionInput& input,
    PressureCorrectionSystemView system,
    PressureCorrectionCertificate& certificate) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  Impl& impl = *implementation_;
  impl.sealed = {};
  impl.current_pressure_work = {};
  const PressureCorrectionBoundaryCertificate& boundary =
      impl.pressure_boundary.certificate();
  if (!impl.pressure_boundary.current() ||
      input.geometry != boundary.geometry ||
      input.numeric_boundary != boundary.source_revision) {
    impl.current_pressure_reference = {};
    impl.current_absolute_pressure_reference = 0.0;
    impl.pressure_input = {};
    impl.pressure_correction = {};
    return {StatusCode::invalid_plan, kPisoPressureBoundary};
  }
  const detail::PressureAssemblyBinding binding{
      impl.kernels, impl.workspace, impl.cells, impl.fingerprint,
      impl.pressure_reference_plan, impl.current,
      impl.current_thermophysical_context,
      impl.pressure_boundary.certificate().geometry_fingerprint};
  PressureCorrectionCertificate candidate;
  const Status status = detail::assemble_pressure_system_impl(
      binding, input, system, candidate);
  if (!status) {
    impl.current_pressure_reference = {};
    impl.current_absolute_pressure_reference = 0.0;
    impl.pressure_input = {};
    impl.pressure_correction = {};
    return status;
  }
  impl.current_pressure_reference = input.pressure_reference;
  impl.current_pressure_work = candidate;
  if (candidate.corrector == 2U) {
    impl.pressure_input = input;
    impl.pressure_correction = candidate;
  } else {
    impl.pressure_input = {};
    impl.pressure_correction = {};
  }
  certificate = candidate;
  return {};
}

Status PressureVelocityCoupler::make_frozen_momentum_stage_authority(
    const PisoIntermediateCertificate& intermediate,
    const PressureCorrectionCertificate& pressure,
    PisoFrozenMomentumStageAuthority& authority) const noexcept {
  authority = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  const bool periodic_scope = entirely_periodic(*impl.boundary);
  const bool supported_scope =
      periodic_scope || supported_frozen_open_boundary(*impl.boundary);
  Status local;
  if (!supported_scope || !impl.frozen_candidate_cartesian ||
      impl.frozen_candidate_baseline == 0U ||
      impl.frozen_candidate_numeric == 0U || !impl.current.valid() ||
      !impl.current_pressure_work.valid() || !intermediate.valid() ||
      !pressure.valid() ||
      !same_intermediate_certificate(intermediate, impl.current) ||
      !same_pressure_certificate(pressure, impl.current_pressure_work) ||
      pressure.intermediate != intermediate.dependency ||
      impl.frozen_candidate_numeric_fingerprint() !=
          impl.frozen_candidate_numeric) {
    local = {StatusCode::invalid_plan, kPisoCoupler};
  }
  int lowest = -1;
  const Status consensus = collective_status(
      impl.communicator, local, impl.rank, impl.size, lowest);
  if (!consensus) return consensus;

  PisoFrozenMomentumStageAuthority candidate;
  candidate.issuer_ = this;
  candidate.intermediate_ = intermediate;
  candidate.pressure_ = pressure;
  candidate.baseline_ = impl.frozen_candidate_baseline;
  candidate.scope_ = periodic_scope
                         ? PisoFrozenMomentumStageScope::cartesian_periodic
                         : PisoFrozenMomentumStageScope::
                               cartesian_open_boundary_ibm;
  candidate.corrector_ = intermediate.corrector;
  std::uint64_t lineage =
      hash_mix(kFnvOffset, UINT64_C(0x76303466726f7a61));
  lineage = hash_mix(lineage, intermediate.dependency);
  lineage = hash_mix(lineage, pressure.state);
  lineage = hash_mix(lineage, impl.frozen_candidate_baseline);
  lineage = hash_mix(lineage, impl.frozen_candidate_numeric);
  if (candidate.scope_ ==
          PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm ||
      impl.immersed_interface != nullptr) {
    lineage = hash_mix(lineage, static_cast<std::uint8_t>(candidate.scope_));
    lineage = hash_mix(
        lineage, impl.pressure_boundary.certificate().semantic);
    lineage = hash_mix(
        lineage, impl.pressure_boundary.certificate().rank_local_layout);
    lineage = hash_mix(
        lineage, static_cast<std::uint8_t>(impl.pressure_reference_kind));
    lineage = hash_mix(
        lineage, double_bits(impl.current_absolute_pressure_reference));
    lineage = hash_mix(
        lineage, impl.candidate_pressure_correction_donor_fingerprint);
    lineage = hash_mix(
        lineage, impl.continuity_activity.collective_fingerprint);
    lineage = hash_mix(
        lineage, impl.immersed_interface == nullptr
                     ? PlanFingerprint{0U}
                     : impl.immersed_interface->fingerprint());
  }
  lineage = hash_mix(lineage, intermediate.corrector);
  candidate.canonical_lineage_ = lineage == 0U ? 1U : lineage;
  Status issuance = candidate.valid()
                        ? Status{}
                        : Status{StatusCode::invalid_plan, kPisoCoupler};
  int issuance_lowest = -1;
  issuance = collective_status(impl.communicator, issuance, impl.rank,
                               impl.size, issuance_lowest);
  if (!issuance) return issuance;
  authority = candidate;
  return {};
}

Status PressureVelocityCoupler::form_frozen_momentum_scaled_pressure(
    const PisoFrozenMomentumStageAuthority& authority,
    ConstFieldView pressure_direction,
    const HaloEngine& candidate_correction_halo,
    double alpha, FieldView& scaled_pressure_correction,
    PisoFrozenMomentumPressureStageCertificate& certificate) const noexcept {
  certificate = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  const Int3 cells = impl.cells;
  const bool current_authority =
      authority.valid() && authority.issuer_ == this && impl.current.valid() &&
      impl.current_pressure_work.valid() &&
      same_intermediate_certificate(authority.intermediate_, impl.current) &&
      same_pressure_certificate(authority.pressure_,
                                impl.current_pressure_work) &&
      authority.baseline_ == impl.frozen_candidate_baseline &&
      authority.canonical_lineage_ != 0U && impl.frozen_candidate_cartesian &&
      impl.frozen_candidate_numeric_fingerprint() ==
          impl.frozen_candidate_numeric;
  const std::array<HaloFieldSpec, 1U> candidate_halo_contract{{
      {scaled_pressure_correction.field, 1U, 1U}}};
  const bool valid_views =
      detail::valid_cell_view(pressure_direction, cells, 0U, 1U, 0U) &&
      pressure_direction.field == impl.correction_field &&
      detail::valid_cell_view(as_const(scaled_pressure_correction), cells, 0U,
                              1U, 1U) &&
      scaled_pressure_correction.field != impl.correction_field &&
      pressure_direction.revision != 0U &&
      !detail::field_views_overlap(pressure_direction,
                                   as_const(scaled_pressure_correction)) &&
      !detail::field_views_overlap(as_const(scaled_pressure_correction),
                                   as_const(impl.workspace.h_by_a)) &&
      !detail::field_views_overlap(as_const(scaled_pressure_correction),
                                   as_const(impl.workspace.r_au)) &&
      candidate_correction_halo.ready() &&
      candidate_correction_halo.instance_identity() != 0U &&
      candidate_correction_halo.instance_identity() !=
          impl.correction_halo->instance_identity() &&
      candidate_correction_halo.ghost_revision(
          scaled_pressure_correction.field) !=
          scaled_pressure_correction.revision &&
      static_cast<bool>(candidate_correction_halo.validate_contract(
          impl.communicator, impl.patch,
          {candidate_halo_contract.data(), candidate_halo_contract.size()},
          impl.boundary->halo_topology()));
  Status local;
  if (!current_authority || !valid_views || !std::isfinite(alpha) ||
      alpha < 0.0 || alpha > kPressureEnergyAitkenMaximumAlpha) {
    local = {StatusCode::invalid_plan, kPisoCoupler};
  }
  if (local) {
    for (std::int32_t z = 0; z < cells.z && local; ++z)
      for (std::int32_t y = 0; y < cells.y && local; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          if (!std::isfinite(
                  pressure_direction.unchecked({x, y, z}, 0U))) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
  }
  int lowest = -1;
  Status consensus = collective_status(
      impl.communicator, local, impl.rank, impl.size, lowest);
  if (!consensus) return consensus;

  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        scaled_pressure_correction.unchecked({x, y, z}, 0U) =
            alpha == 0.0
                ? 0.0
                : alpha * pressure_direction.unchecked({x, y, z}, 0U);

  std::uint64_t direction_numeric =
      hash_mix(kFnvOffset, UINT64_C(0x76303466726f7a64));
  direction_numeric = mix_candidate_field_values(
      direction_numeric, pressure_direction, cells, 1U, 0);
  std::uint64_t correction_direction = hash_mix(
      authority.canonical_lineage_, UINT64_C(0x636f727264697265));
  correction_direction = hash_mix(correction_direction,
                                  pressure_direction.field);
  correction_direction = hash_mix(correction_direction, direction_numeric);
  PlanFingerprint collective_correction_direction = 0U;
  const Status direction_status = impl.collective_candidate_provenance(
      correction_direction, UINT64_C(0x636f727264697267),
      collective_correction_direction);
  if (!direction_status) return direction_status;
  std::uint64_t canonical = hash_mix(
      authority.canonical_lineage_, UINT64_C(0x7363616c65646470));
  canonical = hash_mix(canonical, double_bits(alpha));
  canonical = hash_mix(canonical, direction_numeric);
  canonical = mix_candidate_field_values(
      canonical, as_const(scaled_pressure_correction), cells, 1U, 0);
  std::uint64_t scratch =
      hash_mix(kFnvOffset, UINT64_C(0x7630347363616c65));
  scratch = mix_complete_view_identity(scratch, pressure_direction);
  scratch = mix_complete_view_identity(
      scratch, as_const(scaled_pressure_correction));
  scratch = hash_mix(scratch,
                     candidate_correction_halo.instance_identity());

  PisoFrozenMomentumPressureStageCertificate candidate;
  candidate.issuer_ = this;
  candidate.stage_lineage_ = authority.canonical_lineage_;
  candidate.pressure_direction_ =
      make_piso_field_revision_identity(pressure_direction);
  candidate.scaled_pressure_correction_ = make_piso_field_revision_identity(
      as_const(scaled_pressure_correction));
  candidate.pressure_direction_view_ = pressure_direction;
  candidate.scaled_pressure_correction_view_ =
      as_const(scaled_pressure_correction);
  candidate.halo_instance_ = candidate_correction_halo.instance_identity();
  candidate.alpha_ = alpha;
  candidate.correction_direction_ = collective_correction_direction;
  candidate.canonical_lineage_ = canonical == 0U ? 1U : canonical;
  candidate.scratch_binding_ = scratch == 0U ? 1U : scratch;
  candidate.corrector_ = authority.corrector_;
  Status issuance = candidate.valid()
                        ? Status{}
                        : Status{StatusCode::invalid_plan, kPisoCoupler};
  int issuance_lowest = -1;
  issuance = collective_status(impl.communicator, issuance, impl.rank,
                               impl.size, issuance_lowest);
  if (!issuance) return issuance;
  certificate = candidate;
  return {};
}

Status PressureVelocityCoupler::stage_frozen_momentum_velocity(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
    const HaloEngine& candidate_correction_halo,
    ConstFieldView scaled_pressure_correction, FieldView candidate_velocity,
    PisoFrozenMomentumVelocityStageCertificate& certificate) const noexcept {
  return stage_frozen_momentum_velocity_impl(
      authority, pressure_stage, candidate_correction_halo,
      scaled_pressure_correction, nullptr, candidate_velocity, certificate);
}

Status PressureVelocityCoupler::stage_frozen_momentum_velocity(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
    const HaloEngine& candidate_correction_halo,
    FieldView& scaled_pressure_correction, FieldView candidate_velocity,
    PisoFrozenMomentumVelocityStageCertificate& certificate) const noexcept {
  return stage_frozen_momentum_velocity_impl(
      authority, pressure_stage, candidate_correction_halo,
      as_const(scaled_pressure_correction), &scaled_pressure_correction,
      candidate_velocity, certificate);
}

Status PressureVelocityCoupler::stage_frozen_momentum_velocity_impl(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
    const HaloEngine& candidate_correction_halo,
    ConstFieldView scaled_pressure_correction,
    FieldView* mutable_scaled_pressure_correction,
    FieldView candidate_velocity,
    PisoFrozenMomentumVelocityStageCertificate& certificate) const noexcept {
  certificate = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  const Int3 cells = impl.cells;
  std::uint64_t direction_numeric =
      hash_mix(kFnvOffset, UINT64_C(0x76303466726f7a64));
  direction_numeric = mix_candidate_field_values(
      direction_numeric, pressure_stage.pressure_direction_view_, cells, 1U,
      0);
  std::uint64_t pressure_replay = hash_mix(
      authority.canonical_lineage_, UINT64_C(0x7363616c65646470));
  pressure_replay = hash_mix(pressure_replay,
                             double_bits(pressure_stage.alpha_));
  pressure_replay = hash_mix(pressure_replay, direction_numeric);
  pressure_replay = mix_candidate_field_values(
      pressure_replay, scaled_pressure_correction, cells, 1U, 0);
  const bool current_authority =
      authority.valid() && authority.issuer_ == this &&
      pressure_stage.valid() && pressure_stage.issuer_ == this &&
      pressure_stage.stage_lineage_ == authority.canonical_lineage_ &&
      pressure_stage.corrector_ == authority.corrector_ &&
      pressure_stage.halo_instance_ ==
          candidate_correction_halo.instance_identity() &&
      same_field_identity(scaled_pressure_correction,
                          pressure_stage.scaled_pressure_correction_view_) &&
      (pressure_replay == 0U ? PlanFingerprint{1U} : pressure_replay) ==
          pressure_stage.canonical_lineage_ &&
      impl.current.valid() && impl.current_pressure_work.valid() &&
      same_intermediate_certificate(authority.intermediate_, impl.current) &&
      same_pressure_certificate(authority.pressure_,
                                impl.current_pressure_work) &&
      authority.baseline_ == impl.frozen_candidate_baseline &&
      impl.frozen_candidate_cartesian &&
      impl.frozen_candidate_numeric_fingerprint() ==
          impl.frozen_candidate_numeric;
  const bool open_scope =
      authority.scope_ ==
      PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm;
  const bool ibm_scope = impl.immersed_interface != nullptr;
  const std::array<HaloFieldSpec, 1U> candidate_halo_contract{{
      {scaled_pressure_correction.field, 1U, 1U}}};
  bool valid_views =
      detail::valid_cell_view(scaled_pressure_correction, cells, 0U, 1U, 1U) &&
      scaled_pressure_correction.field != impl.correction_field &&
      detail::valid_cell_view(as_const(candidate_velocity), cells, 0U, 3U,
                              0U) &&
      candidate_velocity.field != impl.boundary->velocity_field() &&
      candidate_velocity.revision != 0U &&
      candidate_correction_halo.ready() &&
      candidate_correction_halo.instance_identity() != 0U &&
      candidate_correction_halo.instance_identity() !=
          impl.correction_halo->instance_identity() &&
      static_cast<bool>(candidate_correction_halo.validate_contract(
          impl.communicator, impl.patch,
          {candidate_halo_contract.data(), candidate_halo_contract.size()},
          impl.boundary->halo_topology())) &&
      candidate_correction_halo.ghost_revision(
          scaled_pressure_correction.field) ==
          scaled_pressure_correction.revision &&
      !detail::field_views_overlap(scaled_pressure_correction,
                                   as_const(candidate_velocity)) &&
      !detail::field_views_overlap(as_const(candidate_velocity),
                                   as_const(impl.workspace.h_by_a)) &&
      !detail::field_views_overlap(as_const(candidate_velocity),
                                   as_const(impl.workspace.r_au));
  if (open_scope || ibm_scope) {
    valid_views =
        valid_views && mutable_scaled_pressure_correction != nullptr &&
        same_field_identity(
            scaled_pressure_correction,
            as_const(*mutable_scaled_pressure_correction));
  }
  if (ibm_scope) {
    valid_views =
        valid_views &&
        impl.candidate_pressure_correction_donors != nullptr &&
        impl.candidate_pressure_correction_donors->ready() &&
        impl.candidate_pressure_correction_donors !=
            impl.pressure_correction_donors &&
        impl.candidate_pressure_correction_donor_stage != 0U &&
        impl.candidate_pressure_correction_field != 0U &&
        scaled_pressure_correction.field ==
            impl.candidate_pressure_correction_field &&
        scaled_pressure_correction.ghosts.x >=
            impl.candidate_pressure_correction_donor_reach &&
        scaled_pressure_correction.ghosts.y >=
            impl.candidate_pressure_correction_donor_reach &&
        scaled_pressure_correction.ghosts.z >=
            impl.candidate_pressure_correction_donor_reach &&
        impl.candidate_pressure_correction_donor_fingerprint != 0U &&
        impl.candidate_pressure_correction_donor_fingerprint ==
            impl.candidate_pressure_correction_donors->fingerprint();
  }
  Status local;
  if (!current_authority || !valid_views)
    local = {StatusCode::invalid_plan, kPisoCoupler};
  int lowest = -1;
  Status consensus = collective_status(
      impl.communicator, local, impl.rank, impl.size, lowest);
  if (!consensus) return consensus;

  if (open_scope) {
    local = impl.pressure_boundary.fill_ghosts(
        *mutable_scaled_pressure_correction);
    consensus = collective_status(
        impl.communicator, local, impl.rank, impl.size, lowest);
    if (!consensus) return consensus;
    scaled_pressure_correction =
        as_const(*mutable_scaled_pressure_correction);
  }

  if (ibm_scope) {
    std::array<FieldView, 1U> donor_views{
        *mutable_scaled_pressure_correction};
    local = impl.candidate_pressure_correction_donors->preflight_exchange(
        impl.candidate_pressure_correction_donor_stage,
        {donor_views.data(), donor_views.size()});
    consensus = collective_status(
        impl.communicator, local, impl.rank, impl.size, lowest);
    if (!consensus) return consensus;
    local = impl.candidate_pressure_correction_donors->exchange(
        impl.candidate_pressure_correction_donor_stage,
        {donor_views.data(), donor_views.size()});
    consensus = collective_status(
        impl.communicator, local, impl.rank, impl.size, lowest);
    if (!consensus) return consensus;
    *mutable_scaled_pressure_correction = donor_views[0U];
    scaled_pressure_correction = as_const(donor_views[0U]);
  }

  if (local) {
    for (std::int32_t z = -1; z <= cells.z && local; ++z)
      for (std::int32_t y = -1; y <= cells.y && local; ++y)
        for (std::int32_t x = -1; x <= cells.x; ++x) {
          const unsigned outside =
              static_cast<unsigned>(x < 0 || x >= cells.x) +
              static_cast<unsigned>(y < 0 || y >= cells.y) +
              static_cast<unsigned>(z < 0 || z >= cells.z);
          if (outside > 1U) continue;
          if (!std::isfinite(
                  scaled_pressure_correction.unchecked({x, y, z}, 0U))) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
        }
  }
  consensus = collective_status(
      impl.communicator, local, impl.rank, impl.size, lowest);
  if (!consensus) return consensus;

  if (!open_scope && pressure_stage.alpha_ == 0.0) {
    for (std::int32_t z = 0; z < cells.z && local; ++z)
      for (std::int32_t y = 0; y < cells.y && local; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          for (std::uint8_t component = 0U; component < 3U; ++component)
            candidate_velocity.unchecked(cell, component) =
                impl.workspace.h_by_a.unchecked(cell, component);
        }
  } else {
    const std::array<ConstFieldView, 1U> reads{scaled_pressure_correction};
    const std::array<FieldView, 1U> writes{candidate_velocity};
    const KernelInvocation invocation{{reads.data(), reads.size()},
                                      {writes.data(), writes.size()},
                                      {{0, 0, 0}, cells},
                                      0U,
                                      0U,
                                      1U,
                                      0U,
                                      nullptr};
    local = cartesian_gradient(*impl.kernels, invocation);
    if (local && ibm_scope)
      local = impl.immersed_interface->correct_pressure_gradient(
          scaled_pressure_correction, candidate_velocity);
    if (local) {
      for (std::int32_t z = 0; z < cells.z && local; ++z)
        for (std::int32_t y = 0; y < cells.y && local; ++y)
          for (std::int32_t x = 0; x < cells.x; ++x) {
            const Int3 cell{x, y, z};
            const std::size_t activity_offset =
                static_cast<std::size_t>(x) +
                static_cast<std::size_t>(cells.x) *
                    (static_cast<std::size_t>(y) +
                     static_cast<std::size_t>(cells.y) *
                         static_cast<std::size_t>(z));
            const bool active = !ibm_scope ||
                                impl.continuity_activity.cells
                                        .data[activity_offset] != 0U;
            for (std::uint8_t component = 0U; component < 3U; ++component) {
              const double value =
                  !active
                      ? 0.0
                      : (pressure_stage.alpha_ == 0.0
                             ? impl.workspace.h_by_a.unchecked(cell,
                                                               component)
                             : impl.workspace.h_by_a.unchecked(
                                   cell, component) -
                                   impl.workspace.r_au.unchecked(
                                       cell, component) *
                                       candidate_velocity.unchecked(
                                           cell, component));
              if (!std::isfinite(value)) {
                local = {StatusCode::numerical_failure, kPisoNumerical};
                break;
              }
              candidate_velocity.unchecked(cell, component) = value;
            }
            if (!local) break;
          }
    }
  }
  consensus = collective_status(
      impl.communicator, local, impl.rank, impl.size, lowest);
  if (!consensus) return consensus;

  std::uint64_t canonical = hash_mix(
      pressure_stage.canonical_lineage_, UINT64_C(0x76656c6f63697479));
  if (open_scope || ibm_scope) {
    canonical = hash_mix(
        canonical, impl.candidate_pressure_correction_donor_fingerprint);
    canonical = hash_mix(
        canonical, impl.candidate_pressure_correction_donor_stage);
    canonical = hash_mix(
        canonical, impl.candidate_pressure_correction_donor_reach);
    canonical = hash_mix(canonical,
                         ibm_scope ? impl.immersed_interface->fingerprint()
                                   : PlanFingerprint{0U});
    canonical = hash_mix(
        canonical, impl.continuity_activity.collective_fingerprint);
  }
  canonical = mix_candidate_field_values(
      canonical, scaled_pressure_correction, cells, 1U, 1);
  canonical = mix_candidate_field_values(
      canonical, as_const(candidate_velocity), cells, 3U, 0);
  std::uint64_t candidate_velocity_numeric =
      hash_mix(kFnvOffset, UINT64_C(0x63616e6476656c6e));
  candidate_velocity_numeric = mix_candidate_field_values(
      candidate_velocity_numeric, as_const(candidate_velocity), cells, 3U,
      0);
  std::uint64_t scratch =
      hash_mix(kFnvOffset, UINT64_C(0x76303466726f7a73));
  scratch = mix_complete_view_identity(scratch, scaled_pressure_correction);
  scratch = hash_mix(scratch,
                     candidate_correction_halo.instance_identity());
  if (open_scope || ibm_scope)
    scratch = hash_mix(
        scratch, impl.candidate_pressure_correction_donor_fingerprint);
  scratch =
      mix_complete_view_identity(scratch, as_const(candidate_velocity));

  PisoFrozenMomentumVelocityStageCertificate candidate;
  candidate.issuer_ = this;
  candidate.stage_lineage_ = authority.canonical_lineage_;
  candidate.pressure_stage_lineage_ = pressure_stage.canonical_lineage_;
  candidate.scaled_pressure_correction_ =
      make_piso_field_revision_identity(scaled_pressure_correction);
  candidate.candidate_velocity_ =
      make_piso_field_revision_identity(as_const(candidate_velocity));
  candidate.scaled_pressure_correction_view_ = scaled_pressure_correction;
  candidate.candidate_velocity_view_ = as_const(candidate_velocity);
  candidate.correction_halo_ = &candidate_correction_halo;
  candidate.halo_instance_ = candidate_correction_halo.instance_identity();
  candidate.halo_ghost_revision_ = candidate_correction_halo.ghost_revision(
      scaled_pressure_correction.field);
  candidate.candidate_donor_fingerprint_ =
      open_scope ? impl.candidate_pressure_correction_donor_fingerprint : 0U;
  candidate.candidate_donor_stage_ =
      open_scope ? impl.candidate_pressure_correction_donor_stage : 0U;
  candidate.candidate_donor_field_ =
      open_scope ? impl.candidate_pressure_correction_field : 0U;
  candidate.candidate_donor_reach_ =
      open_scope ? impl.candidate_pressure_correction_donor_reach : 0U;
  candidate.ibm_geometry_fingerprint_ =
      open_scope && impl.immersed_interface != nullptr
          ? impl.immersed_interface->fingerprint()
          : 0U;
  candidate.candidate_velocity_numeric_ =
      candidate_velocity_numeric == 0U ? 1U : candidate_velocity_numeric;
  candidate.alpha_ = pressure_stage.alpha_;
  candidate.canonical_lineage_ = canonical == 0U ? 1U : canonical;
  candidate.scratch_binding_ = scratch == 0U ? 1U : scratch;
  candidate.corrector_ = authority.corrector_;
  Status issuance = candidate.valid()
                        ? Status{}
                        : Status{StatusCode::invalid_plan, kPisoCoupler};
  int issuance_lowest = -1;
  issuance = collective_status(impl.communicator, issuance, impl.rank,
                               impl.size, issuance_lowest);
  if (!issuance) return issuance;
  certificate = candidate;
  return {};
}

Status PressureVelocityCoupler::stage_frozen_momentum_flux(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumVelocityStageCertificate& velocity,
    ConstFieldView candidate_density, FaceFluxView candidate_flux,
    PisoFrozenMomentumFluxStageCertificate& certificate) const noexcept {
  certificate = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  const Int3 cells = impl.cells;
  const bool current_authority =
      authority.valid() && authority.issuer_ == this && impl.current.valid() &&
      impl.current_pressure_work.valid() &&
      same_intermediate_certificate(authority.intermediate_, impl.current) &&
      same_pressure_certificate(authority.pressure_,
                                impl.current_pressure_work) &&
      authority.baseline_ == impl.frozen_candidate_baseline &&
      authority.canonical_lineage_ == velocity.stage_lineage_ &&
      velocity.valid() && velocity.issuer_ == this &&
      velocity.corrector_ == authority.corrector_ &&
      impl.frozen_candidate_cartesian &&
      impl.frozen_candidate_numeric_fingerprint() ==
          impl.frozen_candidate_numeric;
  const ConstFaceFluxView flux_read = as_const(candidate_flux);
  const std::array<HaloFieldSpec, 1U> candidate_halo_contract{{
      {velocity.scaled_pressure_correction_view_.field, 1U, 1U}}};
  const HaloEngine* const candidate_halo = velocity.correction_halo_;
  const bool matching_halo_authority =
      candidate_halo != nullptr && candidate_halo->ready() &&
      candidate_halo->instance_identity() == velocity.halo_instance_ &&
      candidate_halo->ghost_revision(
          velocity.scaled_pressure_correction_view_.field) ==
          velocity.halo_ghost_revision_ &&
      velocity.halo_ghost_revision_ ==
          velocity.scaled_pressure_correction_view_.revision &&
      static_cast<bool>(candidate_halo->validate_contract(
          impl.communicator, impl.patch,
          {candidate_halo_contract.data(), candidate_halo_contract.size()},
          impl.boundary->halo_topology()));
  bool valid_views =
      detail::valid_cell_view(candidate_density, cells, 0U, 1U, 1U) &&
      candidate_density.field != impl.density_field &&
      detail::valid_flux_view(flux_read, cells, candidate_flux.revision) &&
      !candidate_flux.certificate.valid() && candidate_density.revision != 0U &&
      matching_halo_authority;
  if (valid_views) {
    const std::array<ConstFaceFieldView, 3U> output_faces{
        flux_read.x, flux_read.y, flux_read.z};
    for (ConstFaceFieldView face : output_faces) {
      valid_views =
          valid_views &&
          !detail::cell_face_views_overlap(candidate_density, face) &&
          !detail::cell_face_views_overlap(
              velocity.scaled_pressure_correction_view_, face) &&
          !detail::cell_face_views_overlap(velocity.candidate_velocity_view_,
                                           face) &&
          !detail::face_views_overlap(face,
                                      as_const(impl.workspace.phi_h_by_a.x)) &&
          !detail::face_views_overlap(face,
                                      as_const(impl.workspace.phi_h_by_a.y)) &&
          !detail::face_views_overlap(face,
                                      as_const(impl.workspace.phi_h_by_a.z));
    }
  }

  std::uint64_t replay = hash_mix(
      velocity.pressure_stage_lineage_, UINT64_C(0x76656c6f63697479));
  if (authority.scope_ ==
          PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm ||
      impl.immersed_interface != nullptr) {
    replay = hash_mix(
        replay, impl.candidate_pressure_correction_donor_fingerprint);
    replay = hash_mix(
        replay, impl.candidate_pressure_correction_donor_stage);
    replay = hash_mix(
        replay, impl.candidate_pressure_correction_donor_reach);
    replay = hash_mix(
        replay, impl.immersed_interface == nullptr
                    ? PlanFingerprint{0U}
                    : impl.immersed_interface->fingerprint());
    replay = hash_mix(
        replay, impl.continuity_activity.collective_fingerprint);
  }
  replay = mix_candidate_field_values(
      replay, velocity.scaled_pressure_correction_view_, cells, 1U, 1);
  replay = mix_candidate_field_values(
      replay, velocity.candidate_velocity_view_, cells, 3U, 0);
  Status local;
  if (!current_authority || !valid_views ||
      (replay == 0U ? PlanFingerprint{1U} : replay) !=
          velocity.canonical_lineage_) {
    local = {StatusCode::invalid_plan, kPisoCoupler};
  }
  if (local) {
    for (std::int32_t z = -1; z <= cells.z && local; ++z)
      for (std::int32_t y = -1; y <= cells.y && local; ++y)
        for (std::int32_t x = -1; x <= cells.x; ++x) {
          const unsigned outside =
              static_cast<unsigned>(x < 0 || x >= cells.x) +
              static_cast<unsigned>(y < 0 || y >= cells.y) +
              static_cast<unsigned>(z < 0 || z >= cells.z);
          if (outside > 1U) continue;
          const double value = candidate_density.unchecked({x, y, z}, 0U);
          if (!finite_positive(value)) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
          if (velocity.alpha_ == 0.0 &&
              double_bits(value) !=
                  double_bits(impl.frozen_candidate_density.unchecked(
                      {x, y, z}, 0U))) {
            local = {StatusCode::invalid_plan, kPisoCoupler};
            break;
          }
        }
  }
  int lowest = -1;
  Status consensus = collective_status(
      impl.communicator, local, impl.rank, impl.size, lowest);
  if (!consensus) return consensus;

  const std::array<CartesianAxis, 3U> axes{
      CartesianAxis::x, CartesianAxis::y, CartesianAxis::z};
  const std::array<FaceFieldView, 3U> outputs{
      candidate_flux.x, candidate_flux.y, candidate_flux.z};
  const std::array<ConstFaceFieldView, 3U> bases{
      as_const(impl.workspace.phi_h_by_a.x),
      as_const(impl.workspace.phi_h_by_a.y),
      as_const(impl.workspace.phi_h_by_a.z)};
  const std::array<ConstFaceFieldView, 3U> auxiliaries{
      as_const(impl.frozen_candidate_face_aux.x),
      as_const(impl.frozen_candidate_face_aux.y),
      as_const(impl.frozen_candidate_face_aux.z)};
  const bool open_scope =
      authority.scope_ ==
      PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm;
  for (std::uint8_t axis_index = 0U; axis_index < 3U && local;
       ++axis_index) {
    const CartesianAxis axis = axes[axis_index];
    const FaceFieldView output = outputs[axis_index];
    for (std::int32_t z = 0; z < output.extents.z && local; ++z) {
      for (std::int32_t y = 0; y < output.extents.y && local; ++y) {
        for (std::int32_t x = 0; x < output.extents.x; ++x) {
          const Int3 face{x, y, z};
          PressureCorrectionFaceRule face_rule;
          local = impl.pressure_boundary.face_rule_unchecked(
              axis, face, face_rule);
          bool pressure_outlet = false;
          if (local && face_rule.physical) {
            const BoundaryFacePlan* boundary_face = nullptr;
            local = impl.boundary->face(
                cartesian_face(axis, face_rule.high), boundary_face);
            if (local && boundary_face == nullptr)
              local = {StatusCode::invalid_plan, kPisoCoupler};
            pressure_outlet =
                local && boundary_face->flow_kind ==
                             BoundaryKind::pressure_outlet;
          }
          if (!local) break;
          if (open_scope && face_rule.physical && !pressure_outlet) {
            // Velocity/mass-flow inlets, walls and symmetry faces are owned
            // exclusively by the boundary finalizer.  Zero is a defined
            // scratch placeholder, never publication authority.
            output.unchecked(face) = 0.0;
            continue;
          }
          Int3 left = face;
          if (axis == CartesianAxis::x)
            --left.x;
          else if (axis == CartesianAxis::y)
            --left.y;
          else
            --left.z;
          const std::uint8_t component = axis_index;
          const double coefficient = face_pressure_coefficient(
              *impl.kernels, candidate_density,
              as_const(impl.workspace.r_au), *impl.geometry, impl.patch,
              impl.pressure_boundary, axis, face);
          const double current =
              detail::interpolate_face(
                  *impl.kernels, axis, axis_value(face, axis),
                  candidate_density.unchecked(left, 0U) *
                      impl.workspace.h_by_a.unchecked(left, component),
                  candidate_density.unchecked(face, 0U) *
                      impl.workspace.h_by_a.unchecked(face, component)) *
              detail::face_area(*impl.kernels, axis, face);
          double q = current;
          if (authority.corrector_ == 1U) {
            const double rho_r_au = detail::interpolate_face(
                *impl.kernels, axis, axis_value(face, axis),
                candidate_density.unchecked(left, 0U) *
                    impl.workspace.r_au.unchecked(left, component),
                candidate_density.unchecked(face, 0U) *
                    impl.workspace.r_au.unchecked(face, component));
            q += impl.frozen_candidate_bdf.a0 * rho_r_au *
                 auxiliaries[axis_index].unchecked(face);
            const ConstFieldView predictor_density =
                as_const(impl.frozen_candidate_density);
            const double predictor_coefficient = face_pressure_coefficient(
                *impl.kernels, predictor_density,
                as_const(impl.workspace.r_au), *impl.geometry, impl.patch,
                impl.pressure_boundary, axis, face);
            const double pressure_gradient_flux =
                detail::interpolate_face(
                    *impl.kernels, axis, axis_value(face, axis),
                    predictor_density.unchecked(left, 0U) *
                        impl.workspace.r_au.unchecked(left, component) *
                        impl.workspace.pressure_gradient.unchecked(
                            left, component),
                    predictor_density.unchecked(face, 0U) *
                        impl.workspace.r_au.unchecked(face, component) *
                        impl.workspace.pressure_gradient.unchecked(
                            face, component)) *
                detail::face_area(*impl.kernels, axis, face);
            const double pressure_jump_flux =
                impl.pressure_boundary.mass_flux_response_unchecked(
                    impl.current_pressure_perturbation, axis, face,
                    predictor_coefficient);
            q += pressure_gradient_flux + pressure_jump_flux;
            if (!std::isfinite(predictor_coefficient) ||
                predictor_coefficient < 0.0 ||
                !std::isfinite(pressure_gradient_flux) ||
                !std::isfinite(pressure_jump_flux)) {
              local = {StatusCode::numerical_failure, kPisoNumerical};
              break;
            }
          } else {
            q += auxiliaries[axis_index].unchecked(face);
          }
          const double response =
              impl.pressure_boundary.mass_flux_response_unchecked(
                  velocity.scaled_pressure_correction_view_, axis, face,
                  coefficient);
          const double value = velocity.alpha_ == 0.0
                                   ? bases[axis_index].unchecked(face)
                                   : q + response;
          if (!std::isfinite(coefficient) || coefficient < 0.0 ||
              !std::isfinite(current) || !std::isfinite(q) ||
              !std::isfinite(response) || !std::isfinite(value)) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
          output.unchecked(face) = value;
        }
      }
    }
  }
  if (local && impl.immersed_interface != nullptr) {
    local = impl.immersed_interface->zero_interface_flux(candidate_flux);
    if (local)
      local = impl.immersed_interface->validate_interface_flux(
          as_const(candidate_flux), 0.0);
  }
  consensus = collective_status(
      impl.communicator, local, impl.rank, impl.size, lowest);
  if (!consensus) return consensus;

  std::uint64_t canonical = hash_mix(
      velocity.canonical_lineage_, UINT64_C(0x6d617373666c7578));
  PlanFingerprint nonphysical_flux_provenance = 0U;
  PlanFingerprint pressure_outlet_provisional_provenance = 0U;
  PlanFingerprint ibm_interface_provenance = 0U;
  if (!open_scope) {
    canonical = mix_candidate_field_values(canonical, candidate_density,
                                           cells, 1U, 1);
    canonical = mix_candidate_flux_values(canonical,
                                          as_const(candidate_flux));
  } else {
    std::uint64_t local_nonphysical =
        hash_mix(kFnvOffset, UINT64_C(0x6e6f6e7068797366));
    std::uint64_t local_pressure_outlet =
        hash_mix(kFnvOffset, UINT64_C(0x70726f766f757466));
    for (std::uint8_t axis_index = 0U; axis_index < 3U && local;
         ++axis_index) {
      const CartesianAxis axis = axes[axis_index];
      const ConstFaceFieldView output = as_const(outputs[axis_index]);
      local_nonphysical = hash_mix(local_nonphysical, axis_index);
      local_pressure_outlet = hash_mix(local_pressure_outlet, axis_index);
      for (std::int32_t z = 0; z < output.extents.z && local; ++z)
        for (std::int32_t y = 0; y < output.extents.y && local; ++y)
          for (std::int32_t x = 0; x < output.extents.x; ++x) {
            const Int3 face{x, y, z};
            PressureCorrectionFaceRule face_rule;
            local = impl.pressure_boundary.face_rule_unchecked(
                axis, face, face_rule);
            bool pressure_outlet = false;
            if (local && face_rule.physical) {
              const BoundaryFacePlan* boundary_face = nullptr;
              local = impl.boundary->face(
                  cartesian_face(axis, face_rule.high), boundary_face);
              if (local && boundary_face == nullptr)
                local = {StatusCode::invalid_plan, kPisoCoupler};
              pressure_outlet =
                  local && boundary_face->flow_kind ==
                               BoundaryKind::pressure_outlet;
            }
            if (!local) break;
            const std::uint64_t offset =
                static_cast<std::uint64_t>(x) +
                static_cast<std::uint64_t>(output.extents.x) *
                    (static_cast<std::uint64_t>(y) +
                     static_cast<std::uint64_t>(output.extents.y) *
                         static_cast<std::uint64_t>(z));
            if (!face_rule.physical) {
              local_nonphysical = hash_mix(local_nonphysical, offset);
              local_nonphysical = hash_mix(
                  local_nonphysical, double_bits(output.unchecked(face)));
            } else if (pressure_outlet) {
              local_pressure_outlet =
                  hash_mix(local_pressure_outlet, offset);
              local_pressure_outlet = hash_mix(
                  local_pressure_outlet, double_bits(output.unchecked(face)));
            }
          }
    }
    consensus = collective_status(
        impl.communicator, local, impl.rank, impl.size, lowest);
    if (!consensus) return consensus;
    Status provenance_status = impl.collective_candidate_provenance(
        local_nonphysical, UINT64_C(0x6e6f6e7068797367),
        nonphysical_flux_provenance);
    if (!provenance_status) return provenance_status;
    provenance_status = impl.collective_candidate_provenance(
        local_pressure_outlet, UINT64_C(0x70726f766f757467),
        pressure_outlet_provisional_provenance);
    if (!provenance_status) return provenance_status;
    std::uint64_t local_ibm =
        hash_mix(kFnvOffset, UINT64_C(0x69626d696e746572));
    local_ibm = hash_mix(
        local_ibm, impl.immersed_interface == nullptr
                       ? PlanFingerprint{0U}
                       : impl.immersed_interface->fingerprint());
    local_ibm = hash_mix(
        local_ibm, impl.continuity_activity.local_fingerprint);
    provenance_status = impl.collective_candidate_provenance(
        local_ibm, UINT64_C(0x69626d696e746567),
        ibm_interface_provenance);
    if (!provenance_status) return provenance_status;
    canonical = mix_candidate_field_values(canonical, candidate_density,
                                           cells, 1U, 1);
    canonical = hash_mix(canonical, nonphysical_flux_provenance);
    canonical = hash_mix(canonical,
                         pressure_outlet_provisional_provenance);
    canonical = hash_mix(canonical, ibm_interface_provenance);
  }
  std::uint64_t scratch =
      hash_mix(kFnvOffset, UINT64_C(0x763034666c757873));
  scratch = mix_complete_view_identity(scratch, candidate_density);
  scratch = mix_candidate_flux_binding(scratch, as_const(candidate_flux));
  std::uint64_t candidate_density_numeric =
      hash_mix(kFnvOffset, UINT64_C(0x63616e6464656e73));
  candidate_density_numeric = mix_candidate_field_values(
      candidate_density_numeric, candidate_density, cells, 1U, 1);
  std::uint64_t mechanical_flux_numeric =
      hash_mix(kFnvOffset, UINT64_C(0x6d656368666c786e));
  mechanical_flux_numeric = mix_candidate_flux_values(
      mechanical_flux_numeric, as_const(candidate_flux));

  PisoFrozenMomentumFluxStageCertificate candidate;
  candidate.issuer_ = this;
  candidate.stage_lineage_ = authority.canonical_lineage_;
  candidate.velocity_lineage_ = velocity.canonical_lineage_;
  candidate.candidate_density_ =
      make_piso_field_revision_identity(candidate_density);
  candidate.candidate_density_view_ = candidate_density;
  candidate.candidate_flux_view_ = as_const(candidate_flux);
  candidate.face_flux_revision_ = candidate_flux.revision;
  candidate.face_flux_storage_ = candidate_flux.x.storage_identity;
  candidate.face_flux_revision_domain_ =
      candidate_flux.x.revision_domain;
  candidate.alpha_ = velocity.alpha_;
  candidate.nonphysical_flux_provenance_ =
      nonphysical_flux_provenance;
  candidate.pressure_outlet_provisional_provenance_ =
      pressure_outlet_provisional_provenance;
  candidate.ibm_interface_provenance_ = ibm_interface_provenance;
  candidate.ibm_geometry_fingerprint_ =
      impl.immersed_interface != nullptr
          ? impl.immersed_interface->fingerprint()
          : 0U;
  candidate.candidate_density_numeric_ =
      candidate_density_numeric == 0U ? 1U : candidate_density_numeric;
  candidate.mechanical_flux_numeric_ =
      mechanical_flux_numeric == 0U ? 1U : mechanical_flux_numeric;
  candidate.canonical_lineage_ = canonical == 0U ? 1U : canonical;
  candidate.scratch_binding_ = scratch == 0U ? 1U : scratch;
  candidate.scope_ = authority.scope_;
  candidate.corrector_ = authority.corrector_;
  Status issuance = candidate.valid()
                        ? Status{}
                        : Status{StatusCode::invalid_plan, kPisoCoupler};
  int issuance_lowest = -1;
  issuance = collective_status(impl.communicator, issuance, impl.rank,
                               impl.size, issuance_lowest);
  if (!issuance) return issuance;
  certificate = candidate;
  return {};
}

Status PressureVelocityCoupler::certify_frozen_momentum_exact_baseline(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
    const PisoFrozenMomentumVelocityStageCertificate& velocity_stage,
    const PisoFrozenMomentumFluxStageCertificate& flux_stage,
    PisoFrozenMomentumExactCandidateInput input, ReductionEngine& reductions,
    PisoFrozenMomentumExactCandidateCertificate& certificate) const noexcept {
  return certify_frozen_momentum_exact_candidate_impl(
      authority, nullptr, pressure_stage, velocity_stage, flux_stage, input,
      reductions, certificate);
}

Status PressureVelocityCoupler::certify_frozen_momentum_exact_candidate(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumExactCandidateCertificate& baseline,
    const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
    const PisoFrozenMomentumVelocityStageCertificate& velocity_stage,
    const PisoFrozenMomentumFluxStageCertificate& flux_stage,
    PisoFrozenMomentumExactCandidateInput input, ReductionEngine& reductions,
    PisoFrozenMomentumExactCandidateCertificate& certificate) const noexcept {
  return certify_frozen_momentum_exact_candidate_impl(
      authority, &baseline, pressure_stage, velocity_stage, flux_stage, input,
      reductions, certificate);
}

Status PressureVelocityCoupler::certify_frozen_momentum_exact_candidate_impl(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumExactCandidateCertificate* baseline,
    const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
    const PisoFrozenMomentumVelocityStageCertificate& velocity_stage,
    const PisoFrozenMomentumFluxStageCertificate& flux_stage,
    PisoFrozenMomentumExactCandidateInput input, ReductionEngine& reductions,
    PisoFrozenMomentumExactCandidateCertificate& certificate) const noexcept {
  certificate = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  const Int3 cells = impl.cells;
  const bool open_scope =
      authority.scope_ ==
      PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm;
  bool independent_species_valid =
      impl.thermodynamics_plan != nullptr &&
      input.independent_species.size == impl.independent_species_count &&
      input.thermodynamic.independent_species.size ==
          impl.independent_species_count &&
      (input.independent_species.size == 0U ||
       (input.independent_species.data != nullptr &&
        input.thermodynamic.independent_species.data != nullptr));
  for (std::size_t species = 0U;
       species < input.independent_species.size && independent_species_valid;
       ++species) {
    const ConstFieldView scratch = input.independent_species.data[species];
    const ConstFieldView semantic =
        input.thermodynamic.independent_species.data[species];
    independent_species_valid =
        detail::valid_cell_view(scratch, cells, 0U, 1U, 0U) &&
        detail::valid_cell_view(semantic, cells, 0U, 1U, 0U) &&
        impl.independent_species_semantic_fields != nullptr &&
        semantic.field == impl.independent_species_semantic_fields[species];
  }
  const PlanFingerprint composition_numeric_provenance =
      independent_species_valid
          ? candidate_composition_numeric_local_provenance(
                input.independent_species, cells)
          : PlanFingerprint{};
  const PlanFingerprint semantic_composition_numeric_provenance =
      independent_species_valid
          ? candidate_composition_numeric_local_provenance(
                input.thermodynamic.independent_species, cells)
          : PlanFingerprint{};
  const PlanFingerprint expected_composition_identity =
      independent_species_valid
          ? exact_composition_identity_local(
                impl.thermodynamics, input.independent_species, cells)
          : PlanFingerprint{};
  independent_species_valid =
      independent_species_valid && composition_numeric_provenance != 0U &&
      composition_numeric_provenance ==
          semantic_composition_numeric_provenance &&
      expected_composition_identity != 0U &&
      input.thermodynamic.closure.composition ==
          expected_composition_identity;
  std::uint64_t composition_binding =
      hash_mix(kFnvOffset, UINT64_C(0x636f6d7062696e64));
  if (independent_species_valid)
    for (std::size_t species = 0U; species < input.independent_species.size;
         ++species) {
      composition_binding = mix_complete_view_identity(
          composition_binding, input.independent_species.data[species]);
      composition_binding = mix_complete_view_identity(
          composition_binding,
          input.thermodynamic.independent_species.data[species]);
    }
  if (composition_binding == 0U) composition_binding = 1U;
  int species_lowest = -1;
  Status species_status = collective_status(
      impl.communicator,
      independent_species_valid
          ? Status{}
          : Status{StatusCode::invalid_plan, kPisoCoupler},
      impl.rank, impl.size, species_lowest);
  if (!species_status) return species_status;
  ExactCandidateFluxLocalAudit candidate_flux_audit{};
  bool candidate_flux_audit_valid = false;
  if (open_scope) {
    const FinalBoundaryFluxCertificate& final_boundary =
        input.final_boundary_flux;
    bool matching_species = independent_species_valid &&
                            final_boundary.independent_species_count_ ==
                                input.independent_species.size &&
                            final_boundary.independent_species_views_.size ==
                                input.independent_species.size &&
                            (input.independent_species.size == 0U ||
                             final_boundary.independent_species_views_.data !=
                                 nullptr);
    for (std::size_t species = 0U;
         species < input.independent_species.size && matching_species;
         ++species)
      matching_species = same_field_identity(
          input.independent_species.data[species],
          final_boundary.independent_species_views_.data[species]);
    std::uint64_t final_state_binding =
        hash_mix(kFnvOffset, UINT64_C(0x7065636266736372));
    for (ConstFieldView field :
         {input.candidate_pressure, input.thermodynamic.enthalpy,
          input.thermodynamic.density, input.thermodynamic.temperature,
          input.candidate_velocity})
      final_state_binding = mix_piso_field_revision_identity(
          final_state_binding, make_piso_field_revision_identity(field));
    for (std::size_t species = 0U;
         species < input.independent_species.size; ++species)
      final_state_binding = mix_piso_field_revision_identity(
          final_state_binding,
          make_piso_field_revision_identity(
              input.independent_species.data[species]));
    if (final_state_binding == 0U) final_state_binding = 1U;
    const bool matching_final_boundary =
        final_boundary.valid() && final_boundary.coupler_ == this &&
        final_boundary.scope_ == authority.scope_ &&
        final_boundary.stage_lineage_ == authority.canonical_lineage_ &&
        final_boundary.velocity_lineage_ ==
            velocity_stage.canonical_lineage_ &&
        final_boundary.mechanical_flux_lineage_ ==
            flux_stage.canonical_lineage_ &&
        final_boundary.nonphysical_flux_provenance_ ==
            flux_stage.nonphysical_flux_provenance_ &&
        final_boundary.pressure_outlet_provisional_provenance_ ==
            flux_stage.pressure_outlet_provisional_provenance_ &&
        final_boundary.boundary_semantic_ ==
            impl.boundary->semantic_fingerprint() &&
        final_boundary.boundary_layout_ ==
            impl.boundary->local_layout_fingerprint() &&
        same_pressure_reference_certificate(
            final_boundary.pressure_reference_,
            impl.current_pressure_reference) &&
        double_bits(final_boundary.absolute_pressure_reference_) ==
            double_bits(impl.current_absolute_pressure_reference) &&
        final_boundary.target_time_ == authority.pressure_.time &&
        final_boundary.corrector_ == authority.corrector_ &&
        final_boundary.alpha_ == pressure_stage.alpha_ &&
        same_piso_field_revision_identity(
            final_boundary.candidate_pressure_, input.candidate_pressure) &&
        same_piso_field_revision_identity(
            final_boundary.candidate_enthalpy_,
            input.thermodynamic.enthalpy) &&
        same_piso_field_revision_identity(
            final_boundary.candidate_density_,
            input.thermodynamic.density) &&
        same_piso_field_revision_identity(
            final_boundary.candidate_temperature_,
            input.thermodynamic.temperature) &&
        same_piso_field_revision_identity(
            final_boundary.candidate_velocity_, input.candidate_velocity) &&
        input.thermodynamic.closure.valid() &&
        final_boundary.composition_identity_ ==
            input.thermodynamic.closure.composition &&
        final_boundary.composition_numeric_provenance_ ==
            composition_numeric_provenance &&
        final_boundary.candidate_state_binding_ == final_state_binding &&
        matching_species &&
        same_candidate_flux_identity(
            final_boundary.final_flux_view_, input.candidate_flux);
    Status final_boundary_local =
        matching_final_boundary
            ? Status{}
            : Status{StatusCode::invalid_plan, kPisoCoupler};
    int lowest = -1;
    const Status final_boundary_consensus = collective_status(
        impl.communicator, final_boundary_local, impl.rank, impl.size,
        lowest);
    if (!final_boundary_consensus) return final_boundary_consensus;
    const std::uint64_t local_final_state =
        final_boundary_state_local_provenance(
            input.candidate_pressure, input.thermodynamic.enthalpy,
            input.thermodynamic.density,
            input.thermodynamic.temperature, input.candidate_velocity,
            cells, input.thermodynamic.closure.thermodynamics,
            final_boundary.pressure_reference_.pressure_reference,
            input.thermodynamic.closure.composition,
            final_boundary.absolute_pressure_reference_,
            composition_numeric_provenance);
    candidate_flux_audit =
        exact_candidate_flux_local_audit(input.candidate_flux);
    candidate_flux_audit_valid = true;
    const std::array<Impl::CandidateProvenanceBatchEntry, 2U> replay_entries{{
        {local_final_state, UINT64_C(0x7065636266636f6c)},
        {candidate_flux_audit.final_boundary_provenance,
         UINT64_C(0x7065636266636f6d)},
    }};
    std::array<PlanFingerprint, 2U> replayed{};
    const Status replay_status = impl.collective_candidate_provenance_batch(
        replay_entries.data(), replay_entries.size(), replayed.data());
    if (!replay_status) return replay_status;
    final_boundary_local =
        replayed[0U] ==
                final_boundary.candidate_state_provenance_ &&
            replayed[1U] == final_boundary.final_flux_provenance_
            ? Status{}
            : Status{StatusCode::invalid_plan, kPisoCoupler};
    const Status replay_consensus = collective_status(
        impl.communicator, final_boundary_local, impl.rank, impl.size,
        lowest);
    if (!replay_consensus) return replay_consensus;
  } else {
    const Status scope_status = collective_status(
        impl.communicator,
        input.final_boundary_flux.valid()
            ? Status{StatusCode::invalid_plan, kPisoCoupler}
            : Status{},
        impl.rank, impl.size, species_lowest);
    if (!scope_status) return scope_status;
  }
  const ConstFieldView scaled_pressure_correction =
      input.scaled_pressure_correction;
  const ConstFieldView raw_enthalpy_direction = input.raw_enthalpy_direction;
  const ConstFieldView enthalpy_correction =
      input.scaled_enthalpy_correction;
  const PisoCoupledStateView state = input.base_state;
  const ConstFieldView candidate_pressure = input.candidate_pressure;
  const PisoExactThermodynamicCandidateView thermodynamic_candidate =
      input.thermodynamic;
  const ConstFieldView candidate_velocity = input.candidate_velocity;
  const ConstFaceFluxView candidate_flux = input.candidate_flux;
  const bool baseline_request = baseline == nullptr;
  const ConstFieldView base_velocity = as_const(state.velocity);
  const ConstFieldView base_pressure = as_const(state.pressure_perturbation);
  const ConstFieldView base_enthalpy = as_const(state.enthalpy);
  const ConstFieldView base_density = as_const(state.density);
  const ConstFieldView base_temperature = as_const(state.temperature);
  const PisoExactEosClosureIdentity& closure =
      thermodynamic_candidate.closure;
  const PressureReferenceCertificate& reference =
      impl.current_pressure_reference;
  const bool closed =
      impl.pressure_reference_kind == PressureReferenceKind::closed_mass;
  const PressureEnergyCellActivity gauge_activity{
      impl.continuity_activity.cells,
      impl.continuity_activity.local_fingerprint,
      impl.continuity_activity.collective_fingerprint};
  const ClosedGaugeCorrectionPrepareInput gauge_input{
      reference,
      impl.current_absolute_pressure_reference,
      authority.corrector_,
      authority.pressure_.time,
      authority.pressure_.geometry,
      authority.pressure_.state,
      closure.closure,
      base_pressure,
      scaled_pressure_correction,
      thermodynamic_candidate.pressure_compressibility,
      gauge_activity};
  const bool matching_gauge =
      closed
          ? impl.pressure_reference != nullptr &&
                impl.pressure_reference->matches_closed_gauge_correction(
                    gauge_input, thermodynamic_candidate.closed_gauge)
          : empty_closed_gauge_certificate(
                thermodynamic_candidate.closed_gauge);
  const double gauge_shift =
      closed ? thermodynamic_candidate.closed_gauge.shift : 0.0;
  const double next_absolute_pressure_reference =
      closed ? thermodynamic_candidate.closed_gauge.next_pressure_reference
             : impl.current_absolute_pressure_reference;

  std::uint64_t direction_numeric =
      hash_mix(kFnvOffset, UINT64_C(0x76303466726f7a64));
  std::uint64_t local_direction =
      hash_mix(kFnvOffset, UINT64_C(0x7068646972656374));
  local_direction = hash_mix(local_direction, authority.pressure_.state);
  local_direction = hash_mix(local_direction, authority.pressure_.time);
  mix_candidate_field_values_pair(
      direction_numeric, local_direction,
      pressure_stage.pressure_direction_view_, cells, 1U, 0);
  std::uint64_t pressure_replay = hash_mix(
      authority.canonical_lineage_, UINT64_C(0x7363616c65646470));
  pressure_replay =
      hash_mix(pressure_replay, double_bits(pressure_stage.alpha_));
  pressure_replay = hash_mix(pressure_replay, direction_numeric);
  pressure_replay = mix_candidate_field_values(
      pressure_replay, scaled_pressure_correction, cells, 1U, 0);

  std::uint64_t velocity_replay = hash_mix(
      pressure_stage.canonical_lineage_, UINT64_C(0x76656c6f63697479));
  if (authority.scope_ ==
          PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm ||
      impl.immersed_interface != nullptr) {
    velocity_replay = hash_mix(
        velocity_replay,
        impl.candidate_pressure_correction_donor_fingerprint);
    velocity_replay = hash_mix(
        velocity_replay, impl.candidate_pressure_correction_donor_stage);
    velocity_replay = hash_mix(
        velocity_replay, impl.candidate_pressure_correction_donor_reach);
    velocity_replay = hash_mix(
        velocity_replay, impl.immersed_interface == nullptr
                             ? PlanFingerprint{0U}
                             : impl.immersed_interface->fingerprint());
    velocity_replay = hash_mix(
        velocity_replay, impl.continuity_activity.collective_fingerprint);
  }
  velocity_replay = mix_candidate_field_values(
      velocity_replay, scaled_pressure_correction, cells, 1U, 1);
  velocity_replay = mix_candidate_field_values(
      velocity_replay, candidate_velocity, cells, 3U, 0);

  std::uint64_t flux_replay = hash_mix(
      velocity_stage.canonical_lineage_, UINT64_C(0x6d617373666c7578));
  flux_replay = mix_candidate_field_values(
      flux_replay, thermodynamic_candidate.density, cells, 1U, 1);
  if (open_scope) {
    flux_replay = hash_mix(
        flux_replay, flux_stage.nonphysical_flux_provenance_);
    flux_replay = hash_mix(
        flux_replay,
        flux_stage.pressure_outlet_provisional_provenance_);
    flux_replay = hash_mix(
        flux_replay, flux_stage.ibm_interface_provenance_);
  } else {
    flux_replay = mix_candidate_flux_values(flux_replay, candidate_flux);
  }

  const std::array<HaloFieldSpec, 1U> candidate_halo_contract{{
      {velocity_stage.scaled_pressure_correction_view_.field, 1U, 1U}}};
  const HaloEngine* const candidate_halo = velocity_stage.correction_halo_;
  const bool matching_halo_authority =
      candidate_halo != nullptr && candidate_halo->ready() &&
      candidate_halo->instance_identity() == velocity_stage.halo_instance_ &&
      candidate_halo->ghost_revision(
          velocity_stage.scaled_pressure_correction_view_.field) ==
          velocity_stage.halo_ghost_revision_ &&
      velocity_stage.halo_ghost_revision_ ==
          velocity_stage.scaled_pressure_correction_view_.revision &&
      static_cast<bool>(candidate_halo->validate_contract(
          impl.communicator, impl.patch,
          {candidate_halo_contract.data(), candidate_halo_contract.size()},
          impl.boundary->halo_topology()));

  const bool current_authority =
      authority.valid() && authority.issuer_ == this &&
      pressure_stage.valid() && pressure_stage.issuer_ == this &&
      velocity_stage.valid() && velocity_stage.issuer_ == this &&
      flux_stage.valid() && flux_stage.issuer_ == this &&
      pressure_stage.stage_lineage_ == authority.canonical_lineage_ &&
      velocity_stage.stage_lineage_ == authority.canonical_lineage_ &&
      flux_stage.stage_lineage_ == authority.canonical_lineage_ &&
      velocity_stage.pressure_stage_lineage_ ==
          pressure_stage.canonical_lineage_ &&
      flux_stage.velocity_lineage_ == velocity_stage.canonical_lineage_ &&
      pressure_stage.corrector_ == authority.corrector_ &&
      velocity_stage.corrector_ == authority.corrector_ &&
      flux_stage.corrector_ == authority.corrector_ &&
      pressure_stage.alpha_ == velocity_stage.alpha_ &&
      pressure_stage.alpha_ == flux_stage.alpha_ && impl.current.valid() &&
      impl.current_pressure_work.valid() &&
      same_intermediate_certificate(authority.intermediate_, impl.current) &&
      same_pressure_certificate(authority.pressure_,
                                impl.current_pressure_work) &&
      authority.baseline_ == impl.frozen_candidate_baseline &&
      impl.frozen_candidate_cartesian &&
      impl.frozen_candidate_numeric_fingerprint() ==
          impl.frozen_candidate_numeric &&
      matching_halo_authority &&
      (pressure_replay == 0U ? PlanFingerprint{1U} : pressure_replay) ==
          pressure_stage.canonical_lineage_ &&
      (velocity_replay == 0U ? PlanFingerprint{1U} : velocity_replay) ==
          velocity_stage.canonical_lineage_ &&
      (flux_replay == 0U ? PlanFingerprint{1U} : flux_replay) ==
          flux_stage.canonical_lineage_;

  const bool matching_views =
      same_field_identity(scaled_pressure_correction,
                          pressure_stage.scaled_pressure_correction_view_) &&
      same_field_identity(scaled_pressure_correction,
                          velocity_stage.scaled_pressure_correction_view_) &&
      same_field_identity(candidate_velocity,
                          velocity_stage.candidate_velocity_view_) &&
      same_field_identity(thermodynamic_candidate.density,
                          flux_stage.candidate_density_view_) &&
      (open_scope ||
       same_candidate_flux_identity(candidate_flux,
                                    flux_stage.candidate_flux_view_));
  const bool matching_closure =
      closure.valid() && reference.valid() &&
      closure.thermodynamics == reference.thermodynamics &&
      closure.pressure_reference == reference.pressure_reference &&
      same_piso_field_revision_identity(closure.pressure_state,
                                        base_pressure) &&
      same_piso_field_revision_identity(closure.pressure_correction,
                                        scaled_pressure_correction) &&
      same_piso_field_revision_identity(closure.enthalpy_state,
                                        base_enthalpy) &&
      same_piso_field_revision_identity(closure.enthalpy_correction,
                                        enthalpy_correction) &&
      same_piso_field_revision_identity(closure.candidate_enthalpy,
                                        thermodynamic_candidate.enthalpy) &&
      same_piso_field_revision_identity(closure.candidate_density,
                                        thermodynamic_candidate.density) &&
      same_piso_field_revision_identity(closure.candidate_temperature,
                                        thermodynamic_candidate.temperature);
  const bool matching_c1_output =
      authority.corrector_ == 1U ||
      (impl.current_corrected_c1.valid() &&
       (impl.current_corrected_c1.corrector == 1U ||
        impl.current_corrected_c1.corrector == 2U) &&
       same_pressure_reference_certificate(
           reference,
           impl.current_corrected_c1.output_pressure_reference));

  Status local = reductions.validate_communicator(impl.communicator);
  if (local &&
      (!current_authority || !matching_views || !matching_closure ||
       !matching_gauge || !matching_c1_output ||
       !detail::valid_cell_view(base_velocity, cells, 0U, 3U, 0U) ||
       !detail::valid_cell_view(base_pressure, cells, 0U, 1U, 0U) ||
       !detail::valid_cell_view(base_enthalpy, cells, 0U, 1U, 0U) ||
       !detail::valid_cell_view(base_density, cells, 0U, 1U, 0U) ||
       !detail::valid_cell_view(base_temperature, cells, 0U, 1U, 0U) ||
       !detail::valid_cell_view(scaled_pressure_correction, cells, 0U, 1U,
                               1U) ||
       !detail::valid_cell_view(raw_enthalpy_direction, cells, 0U, 1U, 0U) ||
       !detail::valid_cell_view(enthalpy_correction, cells, 0U, 1U, 0U) ||
       !detail::valid_cell_view(candidate_pressure, cells, 0U, 1U, 0U) ||
       !detail::valid_cell_view(thermodynamic_candidate.enthalpy, cells, 0U,
                               1U, 0U) ||
       !detail::valid_cell_view(thermodynamic_candidate.density, cells, 0U,
                               1U, 1U) ||
       !detail::valid_cell_view(thermodynamic_candidate.temperature, cells,
                               0U, 1U, 0U) ||
       !detail::valid_cell_view(candidate_velocity, cells, 0U, 3U, 0U) ||
       !(closed
             ? detail::valid_cell_view(
                   thermodynamic_candidate.pressure_compressibility, cells,
                   0U, 1U, 0U)
             : empty_field_view(
                   thermodynamic_candidate.pressure_compressibility)) ||
       !detail::valid_flux_view(candidate_flux, cells,
                                candidate_flux.revision))) {
    local = {StatusCode::invalid_plan, kPisoCoupler};
  }

  const std::array<ConstFieldView, 5U> base_fields{
      base_velocity, base_pressure, base_enthalpy, base_density,
      base_temperature};
  std::array<ConstFieldView, 9U> candidate_fields{
      scaled_pressure_correction,
      raw_enthalpy_direction,
      enthalpy_correction,
      candidate_pressure,
      thermodynamic_candidate.enthalpy,
      thermodynamic_candidate.density,
      thermodynamic_candidate.temperature,
      candidate_velocity,
      {}};
  std::size_t candidate_field_count = 8U;
  if (closed) {
    candidate_fields[candidate_field_count++] =
        thermodynamic_candidate.pressure_compressibility;
  }
  bool aliases = false;
  for (std::size_t left = 0U; local && left < base_fields.size(); ++left) {
    for (std::size_t right = left + 1U; right < base_fields.size(); ++right)
      aliases = aliases || detail::field_views_overlap(base_fields[left],
                                                       base_fields[right]);
    for (std::size_t candidate = 0U; candidate < candidate_field_count;
         ++candidate) {
      aliases = aliases || detail::field_views_overlap(base_fields[left],
                                                       candidate_fields[candidate]);
    }
  }
  for (std::size_t left = 0U; local && left < candidate_field_count; ++left)
    for (std::size_t right = left + 1U; right < candidate_field_count;
         ++right)
      aliases = aliases || detail::field_views_overlap(candidate_fields[left],
                                                       candidate_fields[right]);
  const std::array<ConstFaceFieldView, 3U> candidate_faces{
      candidate_flux.x, candidate_flux.y, candidate_flux.z};
  for (ConstFieldView field : base_fields)
    for (ConstFaceFieldView face : candidate_faces)
      aliases = aliases || detail::cell_face_views_overlap(field, face);
  for (std::size_t candidate = 0U; candidate < candidate_field_count;
       ++candidate) {
    const ConstFieldView field = candidate_fields[candidate];
    for (ConstFaceFieldView face : candidate_faces)
      aliases = aliases || detail::cell_face_views_overlap(field, face);
  }
  for (std::size_t left = 0U; local && left < candidate_faces.size(); ++left)
    for (std::size_t right = left + 1U; right < candidate_faces.size();
         ++right)
      aliases = aliases || detail::face_views_overlap(candidate_faces[left],
                                                      candidate_faces[right]);
  for (std::size_t species = 0U;
       local && species < input.independent_species.size; ++species) {
    const ConstFieldView raw = input.independent_species.data[species];
    const ConstFieldView semantic =
        thermodynamic_candidate.independent_species.data[species];
    aliases = aliases || detail::field_views_overlap(raw, semantic);
    for (ConstFieldView base : base_fields) {
      aliases = aliases || detail::field_views_overlap(raw, base) ||
                detail::field_views_overlap(semantic, base);
    }
    for (std::size_t candidate_index = 0U;
         candidate_index < candidate_field_count; ++candidate_index) {
      const ConstFieldView candidate = candidate_fields[candidate_index];
      aliases = aliases || detail::field_views_overlap(raw, candidate) ||
                detail::field_views_overlap(semantic, candidate);
    }
    for (ConstFaceFieldView face : candidate_faces) {
      aliases = aliases || detail::cell_face_views_overlap(raw, face) ||
                detail::cell_face_views_overlap(semantic, face);
    }
    for (std::size_t prior = 0U; prior < species; ++prior) {
      const ConstFieldView prior_raw =
          input.independent_species.data[prior];
      const ConstFieldView prior_semantic =
          thermodynamic_candidate.independent_species.data[prior];
      aliases = aliases || detail::field_views_overlap(raw, prior_raw) ||
                detail::field_views_overlap(raw, prior_semantic) ||
                detail::field_views_overlap(semantic, prior_raw) ||
                detail::field_views_overlap(semantic, prior_semantic);
    }
  }
  if (local && aliases)
    local = {StatusCode::invalid_plan, kPisoCoupler};

  for (std::int32_t z = 0; z < cells.z && local; ++z) {
    for (std::int32_t y = 0; y < cells.y && local; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double old_pressure = base_pressure.unchecked(cell, 0U);
        const double delta_pressure =
            scaled_pressure_correction.unchecked(cell, 0U);
        const double expected_pressure =
            (old_pressure + delta_pressure) - gauge_shift;
        const double selected_pressure =
            candidate_pressure.unchecked(cell, 0U);
        const double pressure_scale =
            std::max({1.0, std::abs(old_pressure),
                      std::abs(delta_pressure), std::abs(expected_pressure),
                      std::abs(selected_pressure)});
        const double pressure_error =
            std::abs(selected_pressure - expected_pressure);
        const double absolute_pressure =
            next_absolute_pressure_reference + selected_pressure;
        const double old_enthalpy = base_enthalpy.unchecked(cell, 0U);
        const double delta_enthalpy =
            enthalpy_correction.unchecked(cell, 0U);
        const double raw_delta_enthalpy =
            raw_enthalpy_direction.unchecked(cell, 0U);
        const double expected_delta_enthalpy =
            pressure_stage.alpha_ == 0.0
                ? 0.0
                : pressure_stage.alpha_ * raw_delta_enthalpy;
        const double delta_enthalpy_scale =
            std::max({1.0, std::abs(delta_enthalpy),
                      std::abs(raw_delta_enthalpy),
                      std::abs(expected_delta_enthalpy)});
        const double delta_enthalpy_error =
            std::abs(delta_enthalpy - expected_delta_enthalpy);
        const double selected_enthalpy =
            thermodynamic_candidate.enthalpy.unchecked(cell, 0U);
        const double enthalpy_scale =
            std::max({1.0, std::abs(old_enthalpy),
                      std::abs(delta_enthalpy),
                      std::abs(selected_enthalpy)});
        const double enthalpy_error =
            std::abs(selected_enthalpy -
                     (old_enthalpy + delta_enthalpy));
        const double selected_density =
            thermodynamic_candidate.density.unchecked(cell, 0U);
        const double selected_temperature =
            thermodynamic_candidate.temperature.unchecked(cell, 0U);
        bool finite_velocity = true;
        for (std::uint8_t component = 0U; component < 3U; ++component)
          finite_velocity =
              finite_velocity &&
              std::isfinite(candidate_velocity.unchecked(cell, component));
        if (!std::isfinite(old_pressure) ||
            !std::isfinite(delta_pressure) ||
            !std::isfinite(expected_pressure) ||
            !std::isfinite(selected_pressure) ||
            !std::isfinite(pressure_error) ||
            pressure_error > 16.0 * std::numeric_limits<double>::epsilon() *
                                 pressure_scale ||
            !std::isfinite(absolute_pressure) ||
            !(absolute_pressure > 0.0) || !std::isfinite(old_enthalpy) ||
            !std::isfinite(raw_delta_enthalpy) ||
            !std::isfinite(delta_enthalpy) ||
            !std::isfinite(expected_delta_enthalpy) ||
            !std::isfinite(delta_enthalpy_error) ||
            delta_enthalpy_error >
                16.0 * std::numeric_limits<double>::epsilon() *
                    delta_enthalpy_scale ||
            !std::isfinite(selected_enthalpy) ||
            !std::isfinite(enthalpy_error) ||
            enthalpy_error > 16.0 * std::numeric_limits<double>::epsilon() *
                                 enthalpy_scale ||
            !finite_positive(selected_density) ||
            !finite_positive(selected_temperature) ||
            !finite_positive(base_density.unchecked(cell, 0U)) ||
            !finite_positive(base_temperature.unchecked(cell, 0U)) ||
            !finite_velocity) {
          local = {StatusCode::numerical_failure, kPisoNumerical};
          break;
        }
        if (!open_scope) {
          for (std::size_t species = 0U;
               species < input.independent_species.size; ++species)
            impl.frozen_candidate_composition_values[species] =
                input.independent_species.data[species].unchecked(cell, 0U);
          const Real3 selected_velocity{
              candidate_velocity.unchecked(cell, 0U),
              candidate_velocity.unchecked(cell, 1U),
              candidate_velocity.unchecked(cell, 2U)};
          ThermoState replayed_thermo;
          const Status replayed =
              impl.thermodynamics_plan->evaluate_from_reference_pressure(
                  next_absolute_pressure_reference, selected_pressure,
                  selected_enthalpy,
                  {impl.frozen_candidate_composition_values.get(),
                   impl.independent_species_count},
                  selected_velocity, replayed_thermo,
                  selected_temperature);
          const double density_scale =
              std::max({1.0, std::abs(selected_density),
                        std::abs(replayed_thermo.rho)});
          const double temperature_scale =
              std::max({1.0, std::abs(selected_temperature),
                        std::abs(replayed_thermo.temperature)});
          const double selected_compressibility =
              thermodynamic_candidate.pressure_compressibility.unchecked(
                  cell, 0U);
          const double compressibility_scale =
              std::max({1.0, std::abs(selected_compressibility),
                        std::abs(replayed_thermo.drho_dp_hY)});
          constexpr double exact_eos_factor =
              64.0 * std::numeric_limits<double>::epsilon();
          if (!replayed || !finite_positive(replayed_thermo.rho) ||
              !finite_positive(replayed_thermo.temperature) ||
              !finite_positive(replayed_thermo.drho_dp_hY) ||
              std::abs(selected_density - replayed_thermo.rho) >
                  exact_eos_factor * density_scale ||
              std::abs(selected_temperature - replayed_thermo.temperature) >
                  exact_eos_factor * temperature_scale ||
              !finite_positive(selected_compressibility) ||
              std::abs(selected_compressibility -
                       replayed_thermo.drho_dp_hY) >
                  exact_eos_factor * compressibility_scale) {
            local = {StatusCode::invalid_plan, kPisoCoupler};
            break;
          }
        }
        if (baseline_request) {
          const double base_absolute_pressure =
              impl.current_absolute_pressure_reference + old_pressure;
          const double candidate_absolute_pressure =
              next_absolute_pressure_reference + selected_pressure;
          const double absolute_pressure_scale =
              std::max({1.0, std::abs(base_absolute_pressure),
                        std::abs(candidate_absolute_pressure)});
          bool exact_base =
              pressure_stage.alpha_ == 0.0 &&
              double_bits(delta_pressure) == 0U &&
              double_bits(delta_enthalpy) == 0U &&
              std::abs(candidate_absolute_pressure - base_absolute_pressure) <=
                  16.0 * std::numeric_limits<double>::epsilon() *
                      absolute_pressure_scale &&
              double_bits(selected_enthalpy) == double_bits(old_enthalpy) &&
              double_bits(selected_density) ==
                  double_bits(base_density.unchecked(cell, 0U)) &&
              double_bits(selected_temperature) ==
                  double_bits(base_temperature.unchecked(cell, 0U));
          for (std::uint8_t component = 0U; component < 3U; ++component)
            exact_base =
                exact_base &&
                double_bits(candidate_velocity.unchecked(cell, component)) ==
                    double_bits(base_velocity.unchecked(cell, component));
          if (!exact_base) {
            local = {StatusCode::invalid_plan, kPisoCoupler};
            break;
          }
        }
      }
    }
  }
  if (local && !candidate_flux_audit_valid) {
    candidate_flux_audit = exact_candidate_flux_local_audit(candidate_flux);
    candidate_flux_audit_valid = true;
  }
  if (local && !candidate_flux_audit.finite)
    local = {StatusCode::numerical_failure, kPisoNumerical};
  if (baseline_request && !open_scope) {
    const std::array<ConstFaceFieldView, 3U> base_faces{
        as_const(impl.workspace.phi_h_by_a.x),
        as_const(impl.workspace.phi_h_by_a.y),
        as_const(impl.workspace.phi_h_by_a.z)};
    for (std::size_t axis = 0U; axis < candidate_faces.size() && local;
         ++axis) {
      const ConstFaceFieldView candidate_face = candidate_faces[axis];
      const ConstFaceFieldView base_face = base_faces[axis];
      for (std::int32_t z = 0; z < candidate_face.extents.z && local; ++z)
        for (std::int32_t y = 0; y < candidate_face.extents.y && local; ++y)
          for (std::int32_t x = 0; x < candidate_face.extents.x; ++x)
            if (double_bits(candidate_face.unchecked({x, y, z})) !=
                double_bits(base_face.unchecked({x, y, z}))) {
              local = {StatusCode::invalid_plan, kPisoCoupler};
              break;
            }
    }
  }

  double local_guard = 0.0;
  double global_guard = 0.0;
  Status status = reductions.checked_max(
      {&local_guard, 1U}, {&global_guard, 1U}, local);
  if (!status) return status;

  const std::uint64_t local_base_state =
      frozen_exact_base_state_local_provenance(
          base_pressure, base_enthalpy, base_density, base_temperature,
          base_velocity, cells, impl.current_absolute_pressure_reference);

  local_direction = mix_candidate_field_values(
      local_direction, raw_enthalpy_direction, cells, 1U, 0);

  std::uint64_t local_state = frozen_exact_candidate_state_local_provenance(
      candidate_pressure, thermodynamic_candidate.enthalpy,
      thermodynamic_candidate.density,
      thermodynamic_candidate.temperature, candidate_velocity, cells,
      thermodynamic_candidate);
  local_state = hash_mix(local_state, composition_numeric_provenance);
  const std::uint64_t local_flux =
      candidate_flux_audit.exact_candidate_provenance;
  const std::uint64_t alpha_bits = double_bits(pressure_stage.alpha_);
  const std::array<Impl::CandidateProvenanceBatchEntry, 4U>
      provenance_entries{{
          {local_base_state, UINT64_C(0x6578616362617367)},
          {local_direction, UINT64_C(0x7068646972676c62)},
          {local_state, UINT64_C(0x6578616374737467), 1U, alpha_bits},
          {local_flux, UINT64_C(0x6578616374666c67), 1U, alpha_bits},
      }};
  std::array<PlanFingerprint, 4U> provenances{};
  status = impl.collective_candidate_provenance_batch(
      provenance_entries.data(), provenance_entries.size(),
      provenances.data());
  if (!status) return status;
  const PlanFingerprint base_state_provenance = provenances[0U];
  const PlanFingerprint correction_direction = provenances[1U];
  const PlanFingerprint candidate_state_provenance = provenances[2U];
  const PlanFingerprint candidate_flux_provenance = provenances[3U];

  const bool baseline_certificate = baseline_request;
  const bool matching_baseline =
      baseline_certificate
          ? pressure_stage.alpha_ == 0.0
          : baseline->valid() && baseline->issuer_ == this &&
                baseline->stage_lineage_ == authority.canonical_lineage_ &&
                baseline->corrector_ == authority.corrector_ &&
                baseline->target_time_ == authority.pressure_.time &&
                baseline->alpha_ == 0.0 && pressure_stage.alpha_ > 0.0 &&
                same_piso_field_revision_identity(baseline->base_velocity_,
                                                  base_velocity) &&
                same_piso_field_revision_identity(baseline->base_pressure_,
                                                  base_pressure) &&
                same_piso_field_revision_identity(baseline->base_enthalpy_,
                                                  base_enthalpy) &&
                same_piso_field_revision_identity(baseline->base_density_,
                                                  base_density) &&
                same_piso_field_revision_identity(baseline->base_temperature_,
                                                  base_temperature) &&
                baseline->base_state_provenance_ == base_state_provenance &&
                baseline->composition_numeric_provenance_ ==
                    composition_numeric_provenance &&
                baseline->independent_species_count_ ==
                    input.independent_species.size &&
                baseline->correction_direction_ == correction_direction &&
                baseline->baseline_state_provenance_ ==
                    baseline->candidate_state_provenance_ &&
                baseline->baseline_mass_flux_provenance_ ==
                    baseline->candidate_mass_flux_provenance_;
  local = matching_baseline
              ? Status{}
              : Status{StatusCode::invalid_plan, kPisoCoupler};
  status = reductions.checked_max(
      {&local_guard, 1U}, {&global_guard, 1U}, local);
  if (!status) return status;
  const PlanFingerprint baseline_state_provenance =
      baseline_certificate ? candidate_state_provenance
                           : baseline->candidate_state_provenance_;
  const PlanFingerprint baseline_flux_provenance =
      baseline_certificate ? candidate_flux_provenance
                           : baseline->candidate_mass_flux_provenance_;

  std::uint64_t exact_lineage = hash_mix(
      authority.canonical_lineage_, UINT64_C(0x657861637463616e));
  exact_lineage = hash_mix(exact_lineage, authority.pressure_.time);
  exact_lineage = hash_mix(exact_lineage, authority.corrector_);
  exact_lineage = hash_mix(exact_lineage, correction_direction);
  exact_lineage = hash_mix(exact_lineage,
                           double_bits(pressure_stage.alpha_));
  exact_lineage = hash_mix(exact_lineage, baseline_state_provenance);
  exact_lineage = hash_mix(exact_lineage, baseline_flux_provenance);
  exact_lineage = hash_mix(exact_lineage, candidate_state_provenance);
  exact_lineage = hash_mix(exact_lineage, candidate_flux_provenance);
  exact_lineage = hash_mix(exact_lineage, base_state_provenance);
  exact_lineage = hash_mix(exact_lineage,
                           composition_numeric_provenance);
  if (open_scope)
    exact_lineage = hash_mix(
        exact_lineage, input.final_boundary_flux.canonical_lineage_);
  std::uint64_t scratch =
      hash_mix(kFnvOffset, UINT64_C(0x6578616374736372));
  for (ConstFieldView field : base_fields)
    scratch = mix_complete_view_identity(scratch, field);
  for (ConstFieldView field : candidate_fields)
    scratch = mix_complete_view_identity(scratch, field);
  scratch = mix_complete_view_identity(
      scratch, pressure_stage.pressure_direction_view_);
  scratch = mix_complete_view_identity(
      scratch, thermodynamic_candidate.pressure_compressibility);
  scratch = mix_candidate_flux_binding(scratch, candidate_flux);
  scratch = hash_mix(scratch, exact_eos_closure_fingerprint(closure));
  for (std::size_t species = 0U; species < input.independent_species.size;
       ++species) {
    scratch = mix_complete_view_identity(
        scratch, input.independent_species.data[species]);
    scratch = mix_complete_view_identity(
        scratch, thermodynamic_candidate.independent_species.data[species]);
  }
  scratch = hash_mix(
      scratch, thermodynamic_candidate.closed_gauge.rank_local_transaction);
  if (open_scope)
    scratch = hash_mix(scratch,
                       input.final_boundary_flux.scratch_binding_);

  PisoFrozenMomentumExactCandidateCertificate issued;
  issued.issuer_ = this;
  issued.stage_lineage_ = authority.canonical_lineage_;
  issued.pressure_lineage_ = pressure_stage.canonical_lineage_;
  issued.velocity_lineage_ = velocity_stage.canonical_lineage_;
  issued.flux_lineage_ = flux_stage.canonical_lineage_;
  issued.pressure_stage_ = pressure_stage;
  issued.velocity_stage_ = velocity_stage;
  issued.flux_stage_ = flux_stage;
  issued.final_boundary_flux_ = input.final_boundary_flux;
  issued.base_velocity_ = make_piso_field_revision_identity(base_velocity);
  issued.base_pressure_ = make_piso_field_revision_identity(base_pressure);
  issued.base_enthalpy_ = make_piso_field_revision_identity(base_enthalpy);
  issued.base_density_ = make_piso_field_revision_identity(base_density);
  issued.base_temperature_ =
      make_piso_field_revision_identity(base_temperature);
  issued.raw_pressure_direction_ = make_piso_field_revision_identity(
      pressure_stage.pressure_direction_view_);
  issued.scaled_pressure_correction_ =
      make_piso_field_revision_identity(scaled_pressure_correction);
  issued.raw_enthalpy_direction_ =
      make_piso_field_revision_identity(raw_enthalpy_direction);
  issued.enthalpy_correction_ =
      make_piso_field_revision_identity(enthalpy_correction);
  issued.candidate_pressure_ =
      make_piso_field_revision_identity(candidate_pressure);
  issued.candidate_velocity_ =
      make_piso_field_revision_identity(candidate_velocity);
  issued.base_velocity_view_ = base_velocity;
  issued.base_pressure_view_ = base_pressure;
  issued.base_enthalpy_view_ = base_enthalpy;
  issued.base_density_view_ = base_density;
  issued.base_temperature_view_ = base_temperature;
  issued.raw_pressure_direction_view_ =
      pressure_stage.pressure_direction_view_;
  issued.scaled_pressure_correction_view_ = scaled_pressure_correction;
  issued.raw_enthalpy_direction_view_ = raw_enthalpy_direction;
  issued.enthalpy_correction_view_ = enthalpy_correction;
  issued.candidate_pressure_view_ = candidate_pressure;
  issued.candidate_velocity_view_ = candidate_velocity;
  issued.thermodynamic_candidate_ = thermodynamic_candidate;
  issued.candidate_flux_view_ = candidate_flux;
  issued.independent_species_views_ = input.independent_species;
  issued.alpha_ = pressure_stage.alpha_;
  issued.target_time_ = authority.pressure_.time;
  issued.correction_direction_ = correction_direction;
  issued.baseline_state_provenance_ = baseline_state_provenance;
  issued.baseline_mass_flux_provenance_ = baseline_flux_provenance;
  issued.candidate_state_provenance_ = candidate_state_provenance;
  issued.candidate_mass_flux_provenance_ = candidate_flux_provenance;
  issued.base_state_provenance_ = base_state_provenance;
  issued.composition_numeric_provenance_ =
      composition_numeric_provenance;
  issued.composition_binding_ = composition_binding;
  issued.independent_species_count_ = input.independent_species.size;
  issued.exact_lineage_ = exact_lineage == 0U ? 1U : exact_lineage;
  issued.scratch_binding_ = scratch == 0U ? 1U : scratch;
  issued.corrector_ = authority.corrector_;
  Status issuance = reductions.consensus(
      issued.valid() ? Status{}
                     : Status{StatusCode::invalid_plan, kPisoCoupler});
  if (!issuance) return issuance;
  impl.current_frozen_exact_lineage = issued.exact_lineage_;
  impl.current_frozen_exact_scratch = issued.scratch_binding_;
  certificate = issued;
  return {};
}

Status PressureVelocityCoupler::certify_frozen_momentum_stationary(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumExactCandidateCertificate& baseline,
    double global_normalized_continuity, double global_normalized_energy,
    ReductionEngine& reductions,
    PressureEnergyStationaryCertificate& certificate) const noexcept {
  certificate = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  const Int3 cells = impl.cells;
  Status local = reductions.validate_communicator(impl.communicator);
  const bool current =
      authority.valid() && authority.issuer_ == this && baseline.valid() &&
      baseline.issuer_ == this && baseline.alpha_ == 0.0 &&
      baseline.stage_lineage_ == authority.canonical_lineage_ &&
      baseline.corrector_ == authority.corrector_ &&
      baseline.target_time_ == authority.pressure_.time &&
      baseline.baseline_state_provenance_ ==
          baseline.candidate_state_provenance_ &&
      baseline.baseline_mass_flux_provenance_ ==
          baseline.candidate_mass_flux_provenance_ && impl.current.valid() &&
      baseline.exact_lineage_ == impl.current_frozen_exact_lineage &&
      baseline.scratch_binding_ == impl.current_frozen_exact_scratch &&
      impl.current_pressure_work.valid() &&
      same_intermediate_certificate(authority.intermediate_, impl.current) &&
      same_pressure_certificate(authority.pressure_,
                                impl.current_pressure_work) &&
      authority.baseline_ == impl.frozen_candidate_baseline &&
      impl.frozen_candidate_numeric_fingerprint() ==
          impl.frozen_candidate_numeric;
  if (local &&
      (!current || !std::isfinite(global_normalized_continuity) ||
       global_normalized_continuity < 0.0 ||
       !std::isfinite(global_normalized_energy) ||
       global_normalized_energy < 0.0 ||
       !std::isfinite(impl.continuity_tolerance) ||
       !(impl.continuity_tolerance > 0.0) ||
       !std::isfinite(impl.energy_tolerance) ||
       !(impl.energy_tolerance > 0.0))) {
    local = {StatusCode::invalid_plan, kPisoCoupler};
  }

  std::uint64_t contract =
      hash_mix(kFnvOffset, UINT64_C(0x73746174696f6e63));
  contract = hash_mix(contract, baseline.target_time_);
  contract = hash_mix(contract, baseline.corrector_);
  contract = hash_mix(contract, baseline.correction_direction_);
  contract = hash_mix(contract, baseline.base_state_provenance_);
  contract = hash_mix(contract, baseline.baseline_state_provenance_);
  contract = hash_mix(contract, baseline.baseline_mass_flux_provenance_);
  contract = hash_mix(contract, baseline.candidate_state_provenance_);
  contract = hash_mix(contract, baseline.candidate_mass_flux_provenance_);
  contract = hash_mix(contract, double_bits(global_normalized_continuity));
  contract = hash_mix(contract, double_bits(global_normalized_energy));
  contract = hash_mix(contract, double_bits(impl.continuity_tolerance));
  contract = hash_mix(contract, double_bits(impl.energy_tolerance));
  double local_maximum[2]{};
  if (local) {
    for (std::int32_t z = 0; z < cells.z && local; ++z) {
      for (std::int32_t y = 0; y < cells.y && local; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const double dp = baseline.scaled_pressure_correction_view_.unchecked(
              cell, 0U);
          const double dh = baseline.enthalpy_correction_view_.unchecked(
              cell, 0U);
          local_maximum[0U] = std::max(local_maximum[0U], std::abs(dp));
          local_maximum[1U] = std::max(local_maximum[1U], std::abs(dh));
          const double base_absolute_pressure =
              impl.current_absolute_pressure_reference +
              baseline.base_pressure_view_.unchecked(cell, 0U);
          const double candidate_absolute_pressure =
              (baseline.flux_stage_.scope_ ==
                       PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm
                   ? impl.current_absolute_pressure_reference
                   : baseline.thermodynamic_candidate_.closed_gauge
                         .next_pressure_reference) +
              baseline.candidate_pressure_view_.unchecked(cell, 0U);
          const double absolute_pressure_scale =
              std::max({1.0, std::abs(base_absolute_pressure),
                        std::abs(candidate_absolute_pressure)});
          bool same_baseline_state =
              std::isfinite(base_absolute_pressure) &&
              std::isfinite(candidate_absolute_pressure) &&
              std::abs(candidate_absolute_pressure - base_absolute_pressure) <=
                  16.0 * std::numeric_limits<double>::epsilon() *
                      absolute_pressure_scale &&
              double_bits(baseline.thermodynamic_candidate_
                              .enthalpy.unchecked(cell, 0U)) ==
                  double_bits(baseline.base_enthalpy_view_.unchecked(cell,
                                                                      0U)) &&
              double_bits(baseline.thermodynamic_candidate_
                              .density.unchecked(cell, 0U)) ==
                  double_bits(baseline.base_density_view_.unchecked(cell,
                                                                     0U)) &&
              double_bits(baseline.thermodynamic_candidate_
                              .temperature.unchecked(cell, 0U)) ==
                  double_bits(baseline.base_temperature_view_.unchecked(
                      cell, 0U));
          for (std::uint8_t component = 0U; component < 3U; ++component)
            same_baseline_state =
                same_baseline_state &&
                double_bits(baseline.candidate_velocity_view_.unchecked(
                    cell, component)) ==
                    double_bits(baseline.base_velocity_view_.unchecked(
                        cell, component));
          if (double_bits(dp) != 0U || double_bits(dh) != 0U ||
              !same_baseline_state) {
            local = {StatusCode::invalid_plan, kPisoCoupler};
            break;
          }
        }
      }
    }
  }
  double global_maximum[2]{};
  Status status = reductions.checked_max(
      {local_maximum, 2U}, {global_maximum, 2U}, local);
  if (!status) return status;
  status = reductions.consensus_contract(contract);
  if (!status) return status;
  if (double_bits(global_maximum[0U]) != 0U ||
      double_bits(global_maximum[1U]) != 0U ||
      global_normalized_continuity > impl.continuity_tolerance ||
      global_normalized_energy > impl.energy_tolerance) {
    return {StatusCode::rejected_step, kPisoCoupler};
  }

  if (baseline.flux_stage_.scope_ ==
      PisoFrozenMomentumStageScope::cartesian_periodic) {
    const std::array<ConstFaceFieldView, 3U> stationary_candidate_faces{
        baseline.candidate_flux_view_.x, baseline.candidate_flux_view_.y,
        baseline.candidate_flux_view_.z};
    const std::array<ConstFaceFieldView, 3U> stationary_base_faces{
        as_const(impl.workspace.phi_h_by_a.x),
        as_const(impl.workspace.phi_h_by_a.y),
        as_const(impl.workspace.phi_h_by_a.z)};
    for (std::size_t axis = 0U;
         axis < stationary_candidate_faces.size() && local; ++axis) {
      const ConstFaceFieldView candidate_face =
          stationary_candidate_faces[axis];
      const ConstFaceFieldView base_face = stationary_base_faces[axis];
      for (std::int32_t z = 0; z < candidate_face.extents.z && local; ++z)
        for (std::int32_t y = 0; y < candidate_face.extents.y && local; ++y)
          for (std::int32_t x = 0; x < candidate_face.extents.x; ++x)
            if (double_bits(candidate_face.unchecked({x, y, z})) !=
                double_bits(base_face.unchecked({x, y, z}))) {
              local = {StatusCode::invalid_plan, kPisoCoupler};
              break;
            }
    }
  }
  status = reductions.checked_max(
      {local_maximum, 1U}, {global_maximum, 1U}, local);
  if (!status) return status;

  const std::uint64_t local_base_state =
      frozen_exact_base_state_local_provenance(
          baseline.base_pressure_view_, baseline.base_enthalpy_view_,
          baseline.base_density_view_, baseline.base_temperature_view_,
          baseline.base_velocity_view_, cells,
          impl.current_absolute_pressure_reference);
  PlanFingerprint replay_base_state = 0U;
  status = impl.collective_candidate_provenance(
      local_base_state, UINT64_C(0x6578616362617367), replay_base_state);
  if (!status) return status;

  std::uint64_t local_direction =
      hash_mix(kFnvOffset, UINT64_C(0x7068646972656374));
  local_direction = hash_mix(local_direction, authority.pressure_.state);
  local_direction = hash_mix(local_direction, authority.pressure_.time);
  local_direction = mix_candidate_field_values(
      local_direction, baseline.raw_pressure_direction_view_, cells, 1U, 0);
  local_direction = mix_candidate_field_values(
      local_direction, baseline.raw_enthalpy_direction_view_, cells, 1U, 0);
  PlanFingerprint replay_direction = 0U;
  status = impl.collective_candidate_provenance(
      local_direction, UINT64_C(0x7068646972676c62), replay_direction);
  if (!status) return status;
  std::uint64_t local_state = frozen_exact_candidate_state_local_provenance(
      baseline.candidate_pressure_view_,
      baseline.thermodynamic_candidate_.enthalpy,
      baseline.thermodynamic_candidate_.density,
      baseline.thermodynamic_candidate_.temperature,
      baseline.candidate_velocity_view_, cells,
      baseline.thermodynamic_candidate_);
  local_state = hash_mix(local_state,
                         baseline.composition_numeric_provenance_);
  local_state = hash_mix(local_state, replay_direction);
  local_state = hash_mix(local_state, double_bits(0.0));
  PlanFingerprint replay_state = 0U;
  status = impl.collective_candidate_provenance(
      local_state, UINT64_C(0x6578616374737467), replay_state);
  if (!status) return status;
  std::uint64_t local_flux =
      hash_mix(kFnvOffset, UINT64_C(0x6578616374666c78));
  local_flux = mix_candidate_flux_values(local_flux,
                                         baseline.candidate_flux_view_);
  local_flux = hash_mix(local_flux, replay_direction);
  local_flux = hash_mix(local_flux, double_bits(0.0));
  PlanFingerprint replay_flux = 0U;
  status = impl.collective_candidate_provenance(
      local_flux, UINT64_C(0x6578616374666c67), replay_flux);
  if (!status) return status;
  local = replay_direction == baseline.correction_direction_ &&
                  replay_base_state == baseline.base_state_provenance_ &&
                  replay_state == baseline.candidate_state_provenance_ &&
                  replay_flux == baseline.candidate_mass_flux_provenance_
              ? Status{}
              : Status{StatusCode::invalid_plan, kPisoCoupler};
  status = reductions.checked_max(
      {local_maximum, 1U}, {global_maximum, 1U}, local);
  if (!status) return status;

  std::uint64_t lineage =
      hash_mix(baseline.exact_lineage_, UINT64_C(0x73746174696f6e79));
  lineage = hash_mix(lineage, double_bits(global_normalized_continuity));
  lineage = hash_mix(lineage, double_bits(global_normalized_energy));
  lineage = hash_mix(lineage, double_bits(impl.continuity_tolerance));
  lineage = hash_mix(lineage, double_bits(impl.energy_tolerance));
  PressureEnergyStationaryCertificate issued;
  issued.issuer_ = this;
  issued.exact_lineage_ = baseline.exact_lineage_;
  issued.target_time_ = baseline.target_time_;
  issued.correction_direction_ = baseline.correction_direction_;
  issued.state_provenance_ = baseline.candidate_state_provenance_;
  issued.mass_flux_provenance_ =
      baseline.candidate_mass_flux_provenance_;
  issued.maximum_pressure_correction_ = global_maximum[0U];
  issued.maximum_enthalpy_correction_ = global_maximum[1U];
  issued.normalized_continuity_ = global_normalized_continuity;
  issued.normalized_energy_ = global_normalized_energy;
  issued.continuity_limit_ = impl.continuity_tolerance;
  issued.energy_limit_ = impl.energy_tolerance;
  issued.stationary_lineage_ = lineage == 0U ? 1U : lineage;
  issued.corrector_ = baseline.corrector_;
  Status issuance = reductions.consensus(
      issued.valid() ? Status{}
                     : Status{StatusCode::invalid_plan, kPisoCoupler});
  if (!issuance) return issuance;
  certificate = issued;
  return {};
}

Status PressureVelocityCoupler::commit_frozen_momentum_exact_state_impl(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumExactCandidateCertificate& exact,
    const PressureEnergyGlobalizationSelectionCertificate* selection,
    const PressureEnergyStationaryCertificate* stationary,
    PisoCoupledStateView state, FaceFluxView output_flux,
    std::uint8_t corrector, ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  certificate = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  Impl& impl = *implementation_;
  const Int3 cells = impl.cells;
  const ConstFieldView base_velocity = as_const(state.velocity);
  const ConstFieldView base_pressure = as_const(state.pressure_perturbation);
  const ConstFieldView base_enthalpy = as_const(state.enthalpy);
  const ConstFieldView base_density = as_const(state.density);
  const ConstFieldView base_temperature = as_const(state.temperature);
  const ConstFaceFluxView candidate_flux = exact.candidate_flux_view_;
  const ConstFaceFluxView output_read = as_const(output_flux);
  const PisoExactThermodynamicCandidateView& thermodynamic =
      exact.thermodynamic_candidate_;
  const PisoExactEosClosureIdentity& closure = thermodynamic.closure;
  const PressureReferenceCertificate& reference =
      impl.current_pressure_reference;
  const bool selected_publication = selection != nullptr;
  const bool stationary_publication = stationary != nullptr;
  bool composition_views_valid =
      exact.independent_species_count_ ==
          exact.independent_species_views_.size &&
      thermodynamic.independent_species.size ==
          exact.independent_species_count_ &&
      (exact.independent_species_count_ == 0U ||
       (exact.independent_species_views_.data != nullptr &&
        thermodynamic.independent_species.data != nullptr));
  for (std::size_t species = 0U;
       species < exact.independent_species_views_.size &&
       composition_views_valid;
       ++species) {
    composition_views_valid =
        detail::valid_cell_view(exact.independent_species_views_.data[species],
                                cells, 0U, 1U, 0U) &&
        detail::valid_cell_view(thermodynamic.independent_species.data[species],
                                cells, 0U, 1U, 0U) &&
        impl.independent_species_semantic_fields != nullptr &&
        thermodynamic.independent_species.data[species].field ==
            impl.independent_species_semantic_fields[species];
  }
  const PlanFingerprint composition_numeric_provenance =
      composition_views_valid
          ? candidate_composition_numeric_local_provenance(
                exact.independent_species_views_, cells)
          : PlanFingerprint{};
  const PlanFingerprint semantic_composition_numeric_provenance =
      composition_views_valid
          ? candidate_composition_numeric_local_provenance(
                thermodynamic.independent_species, cells)
          : PlanFingerprint{};
  const PlanFingerprint semantic_composition_binding =
      composition_views_valid
          ? semantic_composition_rank_local_binding(
                thermodynamic.independent_species)
          : PlanFingerprint{};
  composition_views_valid =
      composition_views_valid && composition_numeric_provenance != 0U &&
      composition_numeric_provenance ==
          semantic_composition_numeric_provenance &&
      semantic_composition_binding != 0U;
  std::uint64_t composition_binding =
      hash_mix(kFnvOffset, UINT64_C(0x636f6d7062696e64));
  if (composition_views_valid)
    for (std::size_t species = 0U;
         species < exact.independent_species_views_.size; ++species) {
      composition_binding = mix_complete_view_identity(
          composition_binding,
          exact.independent_species_views_.data[species]);
      composition_binding = mix_complete_view_identity(
          composition_binding,
          thermodynamic.independent_species.data[species]);
    }
  if (composition_binding == 0U) composition_binding = 1U;

  Status local = reductions.validate_communicator(impl.communicator);
  const bool current =
      authority.valid() && authority.issuer_ == this && exact.valid() &&
      exact.issuer_ == this && exact.stage_lineage_ == authority.canonical_lineage_ &&
      exact.corrector_ == corrector && authority.corrector_ == corrector &&
      exact.target_time_ == authority.pressure_.time && impl.current.valid() &&
      exact.exact_lineage_ == impl.current_frozen_exact_lineage &&
      exact.scratch_binding_ == impl.current_frozen_exact_scratch &&
      impl.current_pressure_work.valid() &&
      same_intermediate_certificate(authority.intermediate_, impl.current) &&
      same_pressure_certificate(authority.pressure_,
                                impl.current_pressure_work) &&
      authority.baseline_ == impl.frozen_candidate_baseline &&
      impl.frozen_candidate_cartesian &&
      impl.frozen_candidate_numeric_fingerprint() ==
          impl.frozen_candidate_numeric;
  const bool selected =
      selected_publication && !stationary_publication && selection->valid() &&
      exact.alpha_ > 0.0 && selection->alpha == exact.alpha_ &&
      selection->corrector == corrector &&
      selection->target_time == exact.target_time_ &&
      selection->correction_direction == exact.correction_direction_ &&
      selection->baseline_state_provenance ==
          exact.baseline_state_provenance_ &&
      selection->baseline_mass_flux_provenance ==
          exact.baseline_mass_flux_provenance_ &&
      selection->candidate_state_provenance ==
          exact.candidate_state_provenance_ &&
      selection->candidate_mass_flux_provenance ==
          exact.candidate_mass_flux_provenance_;
  const bool stationary_selected =
      stationary_publication && !selected_publication && stationary->valid() &&
      stationary->issuer_ == this && exact.alpha_ == 0.0 &&
      stationary->exact_lineage_ == exact.exact_lineage_ &&
      stationary->corrector_ == corrector &&
      stationary->target_time_ == exact.target_time_ &&
      stationary->correction_direction_ == exact.correction_direction_ &&
      stationary->state_provenance_ == exact.candidate_state_provenance_ &&
      stationary->mass_flux_provenance_ ==
          exact.candidate_mass_flux_provenance_ &&
      exact.baseline_state_provenance_ ==
          exact.candidate_state_provenance_ &&
      exact.baseline_mass_flux_provenance_ ==
          exact.candidate_mass_flux_provenance_;
  const bool matching_base =
      same_piso_field_revision_identity(exact.base_velocity_, base_velocity) &&
      same_piso_field_revision_identity(exact.base_pressure_, base_pressure) &&
      same_piso_field_revision_identity(exact.base_enthalpy_, base_enthalpy) &&
      same_piso_field_revision_identity(exact.base_density_, base_density) &&
      same_piso_field_revision_identity(exact.base_temperature_,
                                        base_temperature);
  const bool matching_composition =
      composition_views_valid && composition_numeric_provenance != 0U &&
      composition_numeric_provenance ==
          exact.composition_numeric_provenance_ &&
      composition_binding == exact.composition_binding_ &&
      (!exact.final_boundary_flux_.valid() ||
       (exact.final_boundary_flux_.composition_numeric_provenance_ ==
            composition_numeric_provenance &&
        exact.final_boundary_flux_.independent_species_count_ ==
            exact.independent_species_count_));
  const bool valid_output =
      detail::valid_flux_view(output_read, cells, output_flux.revision) &&
      !output_flux.certificate.valid();
  if (local &&
      (!current || (!selected && !stationary_selected) || !matching_base ||
       !matching_composition ||
       !valid_output || corrector < 1U || corrector > 2U)) {
    local = {StatusCode::invalid_plan, kPisoCoupler};
  }

  const PisoFrozenMomentumPressureStageCertificate& pressure_stage =
      exact.pressure_stage_;
  const PisoFrozenMomentumVelocityStageCertificate& velocity_stage =
      exact.velocity_stage_;
  const PisoFrozenMomentumFluxStageCertificate& flux_stage = exact.flux_stage_;
  std::uint64_t direction_numeric =
      hash_mix(kFnvOffset, UINT64_C(0x76303466726f7a64));
  if (local)
    direction_numeric = mix_candidate_field_values(
        direction_numeric, exact.raw_pressure_direction_view_, cells, 1U, 0);
  std::uint64_t pressure_replay = hash_mix(
      authority.canonical_lineage_, UINT64_C(0x7363616c65646470));
  pressure_replay = hash_mix(pressure_replay, double_bits(exact.alpha_));
  pressure_replay = hash_mix(pressure_replay, direction_numeric);
  if (local)
    pressure_replay = mix_candidate_field_values(
        pressure_replay, exact.scaled_pressure_correction_view_, cells, 1U,
        0);
  std::uint64_t velocity_replay = hash_mix(
      pressure_stage.canonical_lineage_, UINT64_C(0x76656c6f63697479));
  if (authority.scope_ ==
          PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm ||
      impl.immersed_interface != nullptr) {
    velocity_replay = hash_mix(
        velocity_replay,
        impl.candidate_pressure_correction_donor_fingerprint);
    velocity_replay = hash_mix(
        velocity_replay, impl.candidate_pressure_correction_donor_stage);
    velocity_replay = hash_mix(
        velocity_replay, impl.candidate_pressure_correction_donor_reach);
    velocity_replay = hash_mix(
        velocity_replay, impl.immersed_interface == nullptr
                             ? PlanFingerprint{0U}
                             : impl.immersed_interface->fingerprint());
    velocity_replay = hash_mix(
        velocity_replay, impl.continuity_activity.collective_fingerprint);
  }
  if (local) {
    velocity_replay = mix_candidate_field_values(
        velocity_replay, exact.scaled_pressure_correction_view_, cells, 1U,
        1);
    velocity_replay = mix_candidate_field_values(
        velocity_replay, exact.candidate_velocity_view_, cells, 3U, 0);
  }
  std::uint64_t flux_replay = hash_mix(
      velocity_stage.canonical_lineage_, UINT64_C(0x6d617373666c7578));
  if (local) {
    flux_replay = mix_candidate_field_values(
        flux_replay, thermodynamic.density, cells, 1U, 1);
    if (authority.scope_ ==
        PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm) {
      flux_replay = hash_mix(
          flux_replay, flux_stage.nonphysical_flux_provenance_);
      flux_replay = hash_mix(
          flux_replay,
          flux_stage.pressure_outlet_provisional_provenance_);
      flux_replay = hash_mix(
          flux_replay, flux_stage.ibm_interface_provenance_);
    } else {
      flux_replay = mix_candidate_flux_values(flux_replay, candidate_flux);
    }
  }
  const std::array<HaloFieldSpec, 1U> candidate_halo_contract{{
      {velocity_stage.scaled_pressure_correction_view_.field, 1U, 1U}}};
  const HaloEngine* const candidate_halo = velocity_stage.correction_halo_;
  const bool matching_halo_authority =
      candidate_halo != nullptr && candidate_halo->ready() &&
      candidate_halo->instance_identity() == velocity_stage.halo_instance_ &&
      candidate_halo->ghost_revision(
          velocity_stage.scaled_pressure_correction_view_.field) ==
          velocity_stage.halo_ghost_revision_ &&
      velocity_stage.halo_ghost_revision_ ==
          velocity_stage.scaled_pressure_correction_view_.revision &&
      static_cast<bool>(candidate_halo->validate_contract(
          impl.communicator, impl.patch,
          {candidate_halo_contract.data(), candidate_halo_contract.size()},
          impl.boundary->halo_topology()));
  if (local &&
      (!matching_halo_authority ||
       (pressure_replay == 0U ? PlanFingerprint{1U} : pressure_replay) !=
           pressure_stage.canonical_lineage_ ||
       (velocity_replay == 0U ? PlanFingerprint{1U} : velocity_replay) !=
           velocity_stage.canonical_lineage_ ||
       (flux_replay == 0U ? PlanFingerprint{1U} : flux_replay) !=
           flux_stage.canonical_lineage_)) {
    local = {StatusCode::invalid_plan, kPisoCoupler};
  }

  const PressureEnergyCellActivity gauge_activity{
      impl.continuity_activity.cells,
      impl.continuity_activity.local_fingerprint,
      impl.continuity_activity.collective_fingerprint};
  const ClosedGaugeCorrectionPrepareInput gauge_input{
      reference,
      impl.current_absolute_pressure_reference,
      corrector,
      authority.pressure_.time,
      authority.pressure_.geometry,
      authority.pressure_.state,
      closure.closure,
      base_pressure,
      exact.scaled_pressure_correction_view_,
      thermodynamic.pressure_compressibility,
      gauge_activity};
  const bool closed =
      impl.pressure_reference_kind == PressureReferenceKind::closed_mass;
  const bool matching_gauge =
      closed ? impl.pressure_reference != nullptr &&
                   impl.pressure_reference->matches_closed_gauge_correction(
                       gauge_input, thermodynamic.closed_gauge)
             : empty_closed_gauge_certificate(
                   thermodynamic.closed_gauge);
  if (local && !matching_gauge)
    local = {StatusCode::invalid_plan, kPisoCoupler};

  const double gauge_shift =
      closed ? thermodynamic.closed_gauge.shift : 0.0;
  const double next_pressure_reference =
      closed ? thermodynamic.closed_gauge.next_pressure_reference
             : impl.current_absolute_pressure_reference;
  for (std::int32_t z = 0; z < cells.z && local; ++z) {
    for (std::int32_t y = 0; y < cells.y && local; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double dp =
            exact.scaled_pressure_correction_view_.unchecked(cell, 0U);
        const double raw_dh =
            exact.raw_enthalpy_direction_view_.unchecked(cell, 0U);
        const double dh = exact.enthalpy_correction_view_.unchecked(cell, 0U);
        const double expected_dh =
            exact.alpha_ == 0.0 ? 0.0 : exact.alpha_ * raw_dh;
        const double expected_pressure =
            (base_pressure.unchecked(cell, 0U) + dp) - gauge_shift;
        const double candidate_pressure =
            exact.candidate_pressure_view_.unchecked(cell, 0U);
        const double expected_enthalpy =
            base_enthalpy.unchecked(cell, 0U) + dh;
        const double candidate_enthalpy =
            thermodynamic.enthalpy.unchecked(cell, 0U);
        const double pressure_scale =
            std::max({1.0, std::abs(expected_pressure),
                      std::abs(candidate_pressure)});
        const double enthalpy_scale =
            std::max({1.0, std::abs(expected_enthalpy),
                      std::abs(candidate_enthalpy)});
        const double dh_scale =
            std::max({1.0, std::abs(expected_dh), std::abs(dh)});
        bool finite_velocity = true;
        for (std::uint8_t component = 0U; component < 3U; ++component)
          finite_velocity = finite_velocity && std::isfinite(
              exact.candidate_velocity_view_.unchecked(cell, component));
        if (!std::isfinite(raw_dh) || !std::isfinite(dh) ||
            !std::isfinite(expected_dh) ||
            std::abs(dh - expected_dh) >
                16.0 * std::numeric_limits<double>::epsilon() * dh_scale ||
            !std::isfinite(candidate_pressure) ||
            std::abs(candidate_pressure - expected_pressure) >
                16.0 * std::numeric_limits<double>::epsilon() *
                    pressure_scale ||
            !finite_positive(next_pressure_reference + candidate_pressure) ||
            !std::isfinite(candidate_enthalpy) ||
            std::abs(candidate_enthalpy - expected_enthalpy) >
                16.0 * std::numeric_limits<double>::epsilon() *
                    enthalpy_scale ||
            !finite_positive(thermodynamic.density.unchecked(cell, 0U)) ||
            !finite_positive(thermodynamic.temperature.unchecked(cell, 0U)) ||
            !finite_velocity) {
          local = {StatusCode::numerical_failure, kPisoNumerical};
          break;
        }
      }
    }
  }

  const std::array<ConstFaceFieldView, 3U> candidate_faces{
      candidate_flux.x, candidate_flux.y, candidate_flux.z};
  const std::array<ConstFaceFieldView, 3U> output_faces{
      output_read.x, output_read.y, output_read.z};
  const bool same_flux =
      same_candidate_flux_identity(candidate_flux, output_read);
  bool aliases = false;
  for (std::size_t axis = 0U; local && axis < output_faces.size(); ++axis) {
    for (std::size_t other = axis + 1U; other < output_faces.size(); ++other)
      aliases = aliases || detail::face_views_overlap(output_faces[axis],
                                                      output_faces[other]);
    if (!same_flux)
      for (ConstFaceFieldView candidate_face : candidate_faces)
        aliases = aliases || detail::face_views_overlap(output_faces[axis],
                                                        candidate_face);
    for (ConstFieldView field :
         {base_velocity, base_pressure, base_enthalpy, base_density,
          base_temperature, exact.raw_pressure_direction_view_,
          exact.scaled_pressure_correction_view_,
          exact.raw_enthalpy_direction_view_,
          exact.enthalpy_correction_view_, exact.candidate_pressure_view_,
          thermodynamic.enthalpy, thermodynamic.density,
          thermodynamic.temperature, exact.candidate_velocity_view_})
      aliases = aliases || detail::cell_face_views_overlap(field,
                                                           output_faces[axis]);
    if (closed) {
      aliases = aliases || detail::cell_face_views_overlap(
                               thermodynamic.pressure_compressibility,
                               output_faces[axis]);
    }
  }
  const std::array<ConstFieldView, 14U> commit_cells{
      base_velocity,
      base_pressure,
      base_enthalpy,
      base_density,
      base_temperature,
      exact.raw_pressure_direction_view_,
      exact.scaled_pressure_correction_view_,
      exact.raw_enthalpy_direction_view_,
      exact.enthalpy_correction_view_,
      exact.candidate_pressure_view_,
      thermodynamic.enthalpy,
      thermodynamic.density,
      thermodynamic.temperature,
      exact.candidate_velocity_view_};
  if (closed) {
    for (ConstFieldView field : commit_cells) {
      aliases = aliases || detail::field_views_overlap(
                               field,
                               thermodynamic.pressure_compressibility);
    }
  }
  for (std::size_t species = 0U;
       local && species < exact.independent_species_views_.size; ++species) {
    const ConstFieldView raw = exact.independent_species_views_.data[species];
    const ConstFieldView semantic =
        thermodynamic.independent_species.data[species];
    aliases = aliases || detail::field_views_overlap(raw, semantic);
    for (ConstFieldView field : commit_cells) {
      aliases = aliases || detail::field_views_overlap(raw, field) ||
                detail::field_views_overlap(semantic, field);
    }
    if (closed) {
      aliases = aliases || detail::field_views_overlap(
                               raw,
                               thermodynamic.pressure_compressibility) ||
                detail::field_views_overlap(
                    semantic,
                    thermodynamic.pressure_compressibility);
    }
    if (!empty_field_view(thermodynamic.pressure_compressibility)) {
      aliases = aliases || detail::field_views_overlap(
                               raw,
                               thermodynamic.pressure_compressibility) ||
                detail::field_views_overlap(
                    semantic, thermodynamic.pressure_compressibility);
    }
    for (ConstFaceFieldView face : candidate_faces) {
      aliases = aliases || detail::cell_face_views_overlap(raw, face) ||
                detail::cell_face_views_overlap(semantic, face);
    }
    for (ConstFaceFieldView face : output_faces) {
      aliases = aliases || detail::cell_face_views_overlap(raw, face) ||
                detail::cell_face_views_overlap(semantic, face);
    }
    for (std::size_t prior = 0U; prior < species; ++prior) {
      const ConstFieldView prior_raw =
          exact.independent_species_views_.data[prior];
      const ConstFieldView prior_semantic =
          thermodynamic.independent_species.data[prior];
      aliases = aliases || detail::field_views_overlap(raw, prior_raw) ||
                detail::field_views_overlap(raw, prior_semantic) ||
                detail::field_views_overlap(semantic, prior_raw) ||
                detail::field_views_overlap(semantic, prior_semantic);
    }
  }
  if (!empty_field_view(thermodynamic.pressure_compressibility)) {
    for (ConstFieldView field : commit_cells)
      aliases = aliases || detail::field_views_overlap(
                               thermodynamic.pressure_compressibility, field);
    for (ConstFaceFieldView face : candidate_faces)
      aliases = aliases || detail::cell_face_views_overlap(
                               thermodynamic.pressure_compressibility, face);
    for (ConstFaceFieldView face : output_faces)
      aliases = aliases || detail::cell_face_views_overlap(
                               thermodynamic.pressure_compressibility, face);
  }
  if (local && aliases)
    local = {StatusCode::invalid_plan, kPisoCoupler};

  double local_guard = 0.0;
  double global_guard = 0.0;
  Status status = reductions.checked_max(
      {&local_guard, 1U}, {&global_guard, 1U}, local);
  if (!status) return status;

  const std::uint64_t local_base_state =
      frozen_exact_base_state_local_provenance(
          base_pressure, base_enthalpy, base_density, base_temperature,
          base_velocity, cells, impl.current_absolute_pressure_reference);
  PlanFingerprint replay_base_state = 0U;
  status = impl.collective_candidate_provenance(
      local_base_state, UINT64_C(0x6578616362617367), replay_base_state);
  if (!status) return status;

  std::uint64_t local_direction =
      hash_mix(kFnvOffset, UINT64_C(0x7068646972656374));
  local_direction = hash_mix(local_direction, authority.pressure_.state);
  local_direction = hash_mix(local_direction, authority.pressure_.time);
  local_direction = mix_candidate_field_values(
      local_direction, exact.raw_pressure_direction_view_, cells, 1U, 0);
  local_direction = mix_candidate_field_values(
      local_direction, exact.raw_enthalpy_direction_view_, cells, 1U, 0);
  PlanFingerprint replay_direction = 0U;
  status = impl.collective_candidate_provenance(
      local_direction, UINT64_C(0x7068646972676c62), replay_direction);
  if (!status) return status;
  std::uint64_t local_state = frozen_exact_candidate_state_local_provenance(
      exact.candidate_pressure_view_, thermodynamic.enthalpy,
      thermodynamic.density, thermodynamic.temperature,
      exact.candidate_velocity_view_, cells, thermodynamic);
  local_state = hash_mix(local_state, composition_numeric_provenance);
  local_state = hash_mix(local_state, replay_direction);
  local_state = hash_mix(local_state, double_bits(exact.alpha_));
  PlanFingerprint replay_state = 0U;
  status = impl.collective_candidate_provenance(
      local_state, UINT64_C(0x6578616374737467), replay_state);
  if (!status) return status;
  std::uint64_t local_flux =
      hash_mix(kFnvOffset, UINT64_C(0x6578616374666c78));
  local_flux = mix_candidate_flux_values(local_flux, candidate_flux);
  local_flux = hash_mix(local_flux, replay_direction);
  local_flux = hash_mix(local_flux, double_bits(exact.alpha_));
  PlanFingerprint replay_flux = 0U;
  status = impl.collective_candidate_provenance(
      local_flux, UINT64_C(0x6578616374666c67), replay_flux);
  if (!status) return status;
  local = replay_base_state == exact.base_state_provenance_ &&
                  replay_direction == exact.correction_direction_ &&
                  replay_state == exact.candidate_state_provenance_ &&
                  replay_flux == exact.candidate_mass_flux_provenance_
              ? Status{}
              : Status{StatusCode::invalid_plan, kPisoCoupler};
  status = reductions.checked_max(
      {&local_guard, 1U}, {&global_guard, 1U}, local);
  if (!status) return status;

  // All collective and numerical failure points precede this no-fail copy.
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        state.pressure_perturbation.unchecked(cell, 0U) =
            exact.candidate_pressure_view_.unchecked(cell, 0U);
        state.enthalpy.unchecked(cell, 0U) =
            thermodynamic.enthalpy.unchecked(cell, 0U);
        state.density.unchecked(cell, 0U) =
            thermodynamic.density.unchecked(cell, 0U);
        state.temperature.unchecked(cell, 0U) =
            thermodynamic.temperature.unchecked(cell, 0U);
        for (std::uint8_t component = 0U; component < 3U; ++component)
          state.velocity.unchecked(cell, component) =
              exact.candidate_velocity_view_.unchecked(cell, component);
      }
    }
  }
  const std::array<FaceFieldView, 3U> output_writes{
      output_flux.x, output_flux.y, output_flux.z};
  for (std::size_t axis = 0U; axis < output_writes.size(); ++axis) {
    const FaceFieldView output = output_writes[axis];
    const ConstFaceFieldView source = candidate_faces[axis];
    for (std::int32_t z = 0; z < output.extents.z; ++z)
      for (std::int32_t y = 0; y < output.extents.y; ++y)
        for (std::int32_t x = 0; x < output.extents.x; ++x)
          output.unchecked({x, y, z}) = source.unchecked({x, y, z});
  }

  const PressureReferenceCertificate& output_reference =
      closed ? thermodynamic.closed_gauge.output_pressure_reference
             : reference;
  const PlanFingerprint closure_fingerprint =
      exact_eos_closure_fingerprint(closure);
  std::uint64_t state_hash =
      hash_mix(exact.exact_lineage_, UINT64_C(0x7075626c69736865));
  state_hash = hash_mix(
      state_hash,
      selected_publication ? selection->selection_provenance
                           : stationary->stationary_lineage_);
  state_hash = hash_mix(state_hash, output_flux.revision);
  PisoStateCorrectionCertificate committed;
  committed.plan = impl.fingerprint;
  committed.pressure_system = authority.pressure_.state;
  committed.correction = exact.scaled_pressure_correction_view_.revision;
  committed.enthalpy_correction = exact.enthalpy_correction_view_.revision;
  committed.velocity = state.velocity.revision;
  committed.pressure = state.pressure_perturbation.revision;
  committed.enthalpy = state.enthalpy.revision;
  committed.density = state.density.revision;
  committed.temperature = state.temperature.revision;
  committed.face_flux = output_flux.revision;
  committed.exact_eos_closure = closure_fingerprint;
  committed.state = state_hash == 0U ? 1U : state_hash;
  committed.corrector = corrector;
  committed.closure = PisoStateClosure::exact_eos;
  // These tokens authorize the live terminal refresh and therefore retain
  // the C2 live-storage lineage.  The open-boundary candidate scratch and its
  // complete physical closure remain independently sealed by
  // exact.final_boundary_flux_; substituting that scratch-local lineage here
  // would make a valid no-fail copy impossible to audit in the live fields.
  committed.thermophysical_boundary_semantics =
      authority.pressure_.thermophysical_boundary_semantics;
  committed.thermophysical_boundary_target =
      authority.pressure_.thermophysical_boundary_target;
  committed.thermophysical_boundary_rank_local_binding =
      authority.pressure_.thermophysical_boundary_rank_local_binding;
  committed.thermophysical_boundary_collective_lineage =
      authority.pressure_.thermophysical_boundary_collective_lineage;
  committed.thermophysical_boundary_rank_local_lineage =
      authority.pressure_.thermophysical_boundary_rank_local_lineage;
  committed.input_pressure_reference = reference;
  committed.output_pressure_reference = output_reference;
  committed.closed_gauge_collective_transaction =
      closed ? thermodynamic.closed_gauge.collective_transaction : 0U;
  committed.closed_gauge_rank_local_transaction =
      closed ? thermodynamic.closed_gauge.rank_local_transaction : 0U;
  committed.composition_numeric_provenance =
      composition_numeric_provenance;
  committed.composition_rank_local_binding =
      semantic_composition_binding;
  committed.independent_species_count = exact.independent_species_count_;
  if (corrector == 1U) {
    impl.predecessor_c1.valid = true;
    impl.predecessor_c1.alpha = exact.alpha_;
    impl.predecessor_c1.maximum_depletion = 0.0;
    impl.corrected_c1 = committed;
  }
  impl.current_pressure_work = {};
  impl.current = {};
  impl.current_pressure_reference = {};
  impl.current_absolute_pressure_reference = 0.0;
  impl.sealed = {};
  certificate = committed;
  return {};
}

Status PressureVelocityCoupler::commit_frozen_momentum_coupled_trial_state(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumExactCandidateCertificate& candidate,
    const PressureEnergyGlobalizationSelectionCertificate& selection,
    PisoCoupledStateView state, FaceFluxView trial_flux,
    ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  certificate = {};
  if (authority.corrector() != 1U)
    return {StatusCode::invalid_plan, kPisoCoupler};
  return commit_frozen_momentum_exact_state_impl(
      authority, candidate, &selection, nullptr, state, trial_flux, 1U,
      reductions, certificate);
}

Status PressureVelocityCoupler::
    commit_frozen_momentum_coupled_refinement_trial_state(
        const PisoFrozenMomentumStageAuthority& authority,
        const PisoFrozenMomentumExactCandidateCertificate& candidate,
        const PressureEnergyGlobalizationSelectionCertificate& selection,
        PisoCoupledStateView state, FaceFluxView trial_flux,
        std::uint8_t refinement_iteration, ReductionEngine& reductions,
        PisoStateCorrectionCertificate& correction,
        PisoPressureEnergyRefinementStateCertificate& refinement) noexcept {
  correction = {};
  refinement = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  Impl& impl = *implementation_;
  const ConstFaceFluxView trial_flux_read = as_const(trial_flux);
  Status local = reductions.validate_communicator(impl.communicator);
  const bool predecessor_refinement_valid =
      (authority.pressure_.pressure_energy_refinement == 0U &&
       authority.pressure_
               .pressure_energy_refinement_collective_lineage == 0U &&
       authority.pressure_.pressure_energy_refinement_lineage == 0U) ||
      (authority.pressure_.pressure_energy_refinement != 0U &&
       authority.pressure_.pressure_energy_refinement <
           kPressureEnergyRefinementCapacity &&
       authority.pressure_
               .pressure_energy_refinement_collective_lineage != 0U &&
       authority.pressure_.pressure_energy_refinement_lineage != 0U);
  if (local &&
      (authority.corrector() != 2U || refinement_iteration == 0U ||
       refinement_iteration > kPressureEnergyRefinementCapacity ||
       authority.pressure_.pressure_energy_refinement + 1U !=
           refinement_iteration ||
       !predecessor_refinement_valid ||
       !valid_pressure_energy_refinement_state_views(state, impl.cells) ||
       !detail::valid_flux_view(trial_flux_read, impl.cells,
                                trial_flux.revision) ||
       trial_flux.certificate.valid() ||
       !refinement_state_flux_alias_free(state, trial_flux_read) ||
       candidate.target_time_ == 0U ||
       candidate.candidate_state_provenance_ == 0U ||
       candidate.candidate_mass_flux_provenance_ == 0U)) {
    local = {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Status preflight = reductions.consensus(local);
  if (!preflight) return preflight;

  impl.refinement_authority = {};
  impl.refinement_predecessor = {};
  PisoStateCorrectionCertificate committed;
  Status status = commit_frozen_momentum_exact_state_impl(
      authority, candidate, &selection, nullptr, state, trial_flux, 2U,
      reductions, committed);
  if (!status) return status;

  PisoPressureEnergyRefinementStateCertificate issued;
  issued.issuer_ = this;
  issued.predecessor_ = committed;
  issued.velocity_ =
      make_piso_field_revision_identity(as_const(state.velocity));
  issued.pressure_ = make_piso_field_revision_identity(
      as_const(state.pressure_perturbation));
  issued.enthalpy_ =
      make_piso_field_revision_identity(as_const(state.enthalpy));
  issued.density_ =
      make_piso_field_revision_identity(as_const(state.density));
  issued.temperature_ =
      make_piso_field_revision_identity(as_const(state.temperature));
  issued.flux_ = as_const(trial_flux);
  issued.target_time_ = candidate.target_time_;
  issued.state_provenance_ = candidate.candidate_state_provenance_;
  issued.mass_flux_provenance_ = candidate.candidate_mass_flux_provenance_;
  issued.rank_local_state_numeric_ =
      pressure_energy_refinement_state_numeric(state, impl.cells);
  issued.rank_local_flux_numeric_ =
      pressure_energy_refinement_flux_numeric(as_const(trial_flux));
  issued.iteration_ = refinement_iteration;
  // candidate_state_provenance_ and candidate_mass_flux_provenance_ were
  // already agreed collectively while certifying/replaying the exact
  // candidate.  Deriving this public lineage is therefore rank invariant and
  // requires no additional collective on the post-publication hot path.
  std::uint64_t collective_lineage =
      hash_mix(kFnvOffset, UINT64_C(0x707266636f6c6c31));
  collective_lineage = hash_mix(collective_lineage, 1U);
  collective_lineage = hash_mix(collective_lineage, issued.target_time_);
  collective_lineage = hash_mix(collective_lineage, issued.iteration_);
  collective_lineage =
      hash_mix(collective_lineage, issued.state_provenance_);
  collective_lineage =
      hash_mix(collective_lineage, issued.mass_flux_provenance_);
  issued.collective_lineage_ =
      collective_lineage == 0U ? PlanFingerprint{1U} : collective_lineage;
  std::uint64_t lineage =
      hash_mix(kFnvOffset, UINT64_C(0x7065726566696e65));
  lineage = hash_mix(lineage, impl.fingerprint);
  lineage = hash_mix(lineage, issued.target_time_);
  lineage = hash_mix(lineage, issued.iteration_);
  lineage = hash_mix(lineage, issued.predecessor_.state);
  lineage = hash_mix(lineage, issued.state_provenance_);
  lineage = hash_mix(lineage, issued.mass_flux_provenance_);
  issued.lineage_ = lineage == 0U ? PlanFingerprint{1U} : lineage;
  // The exact commit above crossed the no-fail live-state publication
  // boundary.  Every field/flux identity and issuance input was collectively
  // preflighted before that copy; certificate construction is now an
  // invariant-only operation and must not introduce a post-publication error
  // path.
  assert(issued.valid());
  impl.refinement_predecessor = committed;
  impl.refinement_authority = issued;
  correction = committed;
  refinement = issued;
  return {};
}

Status PressureVelocityCoupler::commit_frozen_momentum_stationary_trial_state(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumExactCandidateCertificate& baseline,
    const PressureEnergyStationaryCertificate& stationary,
    PisoCoupledStateView state, FaceFluxView trial_flux,
    ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  certificate = {};
  if (authority.corrector() != 1U)
    return {StatusCode::invalid_plan, kPisoCoupler};
  return commit_frozen_momentum_exact_state_impl(
      authority, baseline, nullptr, &stationary, state, trial_flux, 1U,
      reductions, certificate);
}

Status PressureVelocityCoupler::commit_frozen_momentum_coupled_pending_state(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumExactCandidateCertificate& candidate,
    const PressureEnergyGlobalizationSelectionCertificate& selection,
    PisoCoupledStateView state, PendingFaceFluxView& pending_flux,
    ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  certificate = {};
  if (implementation_ != nullptr) {
    Impl& impl = *implementation_;
    impl.corrected_pending = {};
    impl.terminal = {};
    impl.sealed = {};
    impl.refinement_authority = {};
    impl.refinement_predecessor = {};
    impl.refinement_root_c1 = {};
  }
  if (implementation_ == nullptr || authority.corrector() != 2U ||
      !pending_flux.valid() || pending_flux.writer_ == nullptr ||
      pending_flux.storage_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  FaceFluxView output{pending_flux.x_, pending_flux.y_, pending_flux.z_,
                      pending_flux.revision_, {}};
  const Status status = commit_frozen_momentum_exact_state_impl(
      authority, candidate, &selection, nullptr, state, output, 2U,
      reductions, certificate);
  if (status) {
    Impl& impl = *implementation_;
    impl.terminal = {};
    impl.corrected_pending = {certificate,
                              &pending_flux,
                              pending_flux.writer_,
                              pending_flux.storage_,
                              pending_flux.revision_,
                              pending_flux.writer_identity_,
                              pending_flux.attempt_identity_};
    impl.corrected_pending.pressure_compressibility =
        make_piso_field_revision_identity(
            candidate.thermodynamic_candidate_.pressure_compressibility);
    impl.corrected_pending.active_cells = impl.continuity_activity.cells;
    impl.corrected_pending.activity_local_fingerprint =
        impl.continuity_activity.local_fingerprint;
    impl.corrected_pending.activity_collective_fingerprint =
        impl.continuity_activity.collective_fingerprint;
  }
  return status;
}

Status PressureVelocityCoupler::commit_frozen_momentum_stationary_pending_state(
    const PisoFrozenMomentumStageAuthority& authority,
    const PisoFrozenMomentumExactCandidateCertificate& baseline,
    const PressureEnergyStationaryCertificate& stationary,
    PisoCoupledStateView state, PendingFaceFluxView& pending_flux,
    ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  certificate = {};
  if (implementation_ != nullptr) {
    Impl& impl = *implementation_;
    impl.corrected_pending = {};
    impl.terminal = {};
    impl.sealed = {};
    impl.refinement_authority = {};
    impl.refinement_predecessor = {};
    impl.refinement_root_c1 = {};
  }
  if (implementation_ == nullptr || authority.corrector() != 2U ||
      !pending_flux.valid() || pending_flux.writer_ == nullptr ||
      pending_flux.storage_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  FaceFluxView output{pending_flux.x_, pending_flux.y_, pending_flux.z_,
                      pending_flux.revision_, {}};
  const Status status = commit_frozen_momentum_exact_state_impl(
      authority, baseline, nullptr, &stationary, state, output, 2U,
      reductions, certificate);
  if (status) {
    Impl& impl = *implementation_;
    impl.terminal = {};
    impl.corrected_pending = {certificate,
                              &pending_flux,
                              pending_flux.writer_,
                              pending_flux.storage_,
                              pending_flux.revision_,
                              pending_flux.writer_identity_,
                              pending_flux.attempt_identity_};
    impl.corrected_pending.pressure_compressibility =
        make_piso_field_revision_identity(
            baseline.thermodynamic_candidate_.pressure_compressibility);
    impl.corrected_pending.active_cells = impl.continuity_activity.cells;
    impl.corrected_pending.activity_local_fingerprint =
        impl.continuity_activity.local_fingerprint;
    impl.corrected_pending.activity_collective_fingerprint =
        impl.continuity_activity.collective_fingerprint;
  }
  return status;
}

Status PressureVelocityCoupler::bind_pressure_operator(
    PressureOperatorServices services,
    PressureCorrectionSystemView system,
    PressureLinearOperator& out) const noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  return PressureLinearOperator::bind_internal(
      *impl.geometry, impl.patch, *impl.boundary, impl.fingerprint,
      as_const(impl.workspace.x_pressure_coefficient),
      as_const(impl.workspace.y_pressure_coefficient),
      as_const(impl.workspace.z_pressure_coefficient), services, system, out);
}

Status PressureVelocityCoupler::make_native_pressure_mg_spec(
    MPI_Comm communicator, LinearIdentity identity,
    MgCoefficientIdentity coefficients,
    NativeCartesianMgSpec& out) const noexcept {
  if (implementation_ == nullptr || communicator == MPI_COMM_NULL ||
      identity.fingerprint == 0U || coefficients.revision == 0U ||
      coefficients.fingerprint == 0U) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }
  const Impl& impl = *implementation_;
  MgBoundarySet boundaries;
  const Status boundary_status =
      impl.pressure_boundary.mg_boundaries(boundaries);
  if (!boundary_status) return boundary_status;
  NativeCartesianMgSpec candidate;
  candidate.communicator = communicator;
  candidate.geometry = impl.geometry;
  candidate.patch = impl.patch;
  candidate.boundaries = boundaries;
  // The low-Mach pressure increment has positive local EOS storage
  // a0*V*(drho/dp) even in a closed domain. The p_ref/pi representation still
  // needs a gauge normalization, but the increment operator itself has no
  // constant null space.
  candidate.null_space = MgNullSpace::none;
  candidate.operator_class =
      MgOperatorClass::symmetric_diagonally_dominant_m_matrix;
  candidate.policy.pre_sweeps = 1U;
  candidate.policy.post_sweeps = 2U;
  candidate.policy.point_smoother = MgPointSmootherKind::chebyshev_jacobi;
  candidate.policy.cycle = MgCycleKind::f_cycle;
  candidate.policy.chebyshev_lower_spectrum_fraction = 0.3;
  candidate.correction_scaling = impl.mg_correction_scaling;
  candidate.identity = identity;
  candidate.coefficients = coefficients;
  out = candidate;
  return {};
}

Status PressureVelocityCoupler::compile_native_pressure_mg(
    const PressureCorrectionCertificate& pressure,
    const NativeCartesianMgSpec& spec,
    MgRuntimeServices services,
    PressureCorrectionSystemView system,
    NativeCartesianMgPlan& out,
    MgPlanCounters* counters) const noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }
  const Impl& impl = *implementation_;
  MgBoundarySet expected_boundaries;
  const Status boundary_status =
      impl.pressure_boundary.mg_boundaries(expected_boundaries);
  if (!pressure.valid() || pressure.plan != impl.fingerprint ||
      pressure.intermediate != impl.current.dependency ||
      pressure.thermophysical_boundary_semantics !=
          impl.current.thermophysical_boundary_semantics ||
      pressure.thermophysical_boundary_target !=
          impl.current.thermophysical_boundary_target ||
      pressure.thermophysical_boundary_rank_local_binding !=
          impl.current.thermophysical_boundary_rank_local_binding ||
      pressure.thermophysical_boundary_collective_lineage !=
          impl.current.thermophysical_boundary_collective_lineage ||
      pressure.thermophysical_boundary_rank_local_lineage !=
          impl.current.thermophysical_boundary_rank_local_lineage ||
      spec.geometry != impl.geometry || !same_patch(spec.patch, impl.patch) ||
      !boundary_status ||
      !same_boundaries(spec.boundaries, expected_boundaries) ||
      spec.correction_scaling != impl.mg_correction_scaling ||
      spec.identity.fingerprint == 0U || spec.coefficients.revision == 0U ||
      spec.coefficients.fingerprint != pressure.state ||
      !detail::valid_cell_view(system.diagonal, impl.cells, 0U, 1U) ||
      !detail::valid_cell_view(system.rhs, impl.cells, 0U, 1U)) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }
  return NativeCartesianMgPlan::compile(
      spec, services,
      {as_const(system.diagonal),
       as_const(impl.workspace.x_pressure_coefficient),
       as_const(impl.workspace.y_pressure_coefficient),
       as_const(impl.workspace.z_pressure_coefficient)},
      out, counters);
}

Status PressureVelocityCoupler::refresh_pressure_linear_lifecycle(
    const PressureCorrectionCertificate& pressure,
    LinearIdentity identity,
    MgCoefficientIdentity coefficients,
    PressureCorrectionSystemView system,
    PressureLinearOperator& linear_operator,
    NativeCartesianMgPlan& preconditioner,
    MgPlanCounters* counters, bool collective_fail_close) const noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }
  Impl& impl = *implementation_;
  impl.sealed = {};
  Status local;
  if (!impl.pressure_boundary.current() || !pressure.valid() ||
      pressure.plan != impl.fingerprint ||
      pressure.intermediate != impl.current.dependency ||
      pressure.thermophysical_boundary_semantics !=
          impl.current.thermophysical_boundary_semantics ||
      pressure.thermophysical_boundary_target !=
          impl.current.thermophysical_boundary_target ||
      pressure.thermophysical_boundary_rank_local_binding !=
          impl.current.thermophysical_boundary_rank_local_binding ||
      pressure.thermophysical_boundary_collective_lineage !=
          impl.current.thermophysical_boundary_collective_lineage ||
      pressure.thermophysical_boundary_rank_local_lineage !=
          impl.current.thermophysical_boundary_rank_local_lineage ||
      identity.fingerprint == 0U || coefficients.revision == 0U ||
      coefficients.fingerprint != pressure.state ||
      !detail::valid_cell_view(system.diagonal, impl.cells, 0U, 1U) ||
      !detail::valid_cell_view(system.rhs, impl.cells, 0U, 1U) ||
      linear_operator.fingerprint() == 0U ||
      preconditioner.certificate().identity.fingerprint == 0U) {
    local = {StatusCode::invalid_plan, kPisoSolve};
  }
  Status status = local;
  if (collective_fail_close) {
    int validation_lowest = -1;
    status = collective_status(impl.communicator, local, impl.rank,
                               impl.size, validation_lowest);
  }
  if (!status) return status;
  const MgCoefficientViews views{
      as_const(system.diagonal),
      as_const(impl.workspace.x_pressure_coefficient),
      as_const(impl.workspace.y_pressure_coefficient),
      as_const(impl.workspace.z_pressure_coefficient)};
  status = preconditioner.update_coefficients(
      identity, coefficients, views, counters);
  if (!status) {
    return status;
  }
  const LinearPreconditionerCertificate preconditioner_certificate =
      preconditioner.certificate();
  const LinearPreconditionerClass expected_class =
      impl.pressure_algorithm == LinearAlgorithm::bicgstab
          ? LinearPreconditionerClass::fixed_general
          : LinearPreconditionerClass::flexible;
  local = {};
  if (preconditioner_certificate.identity.fingerprint != identity.fingerprint ||
      preconditioner_certificate.collective_fingerprint == 0U ||
      preconditioner_certificate.preconditioner_class != expected_class) {
    local = {StatusCode::invalid_plan, kPisoSolve};
  }
  status = local;
  if (collective_fail_close) {
    int certificate_lowest = -1;
    status = collective_status(impl.communicator, local, impl.rank, impl.size,
                               certificate_lowest);
  }
  if (!status) return status;
  local = linear_operator.refresh(
      {pressure, identity,
       preconditioner_certificate.collective_fingerprint});
  if (!collective_fail_close) return local;
  int refresh_lowest = -1;
  return collective_status(impl.communicator, local, impl.rank, impl.size,
                           refresh_lowest);
}

Status PressureVelocityCoupler::correct_state_impl(
    const PressureCorrectionCertificate& pressure,
    FieldView correction,
    PisoTrialStateView state,
    FaceFluxView trial_flux,
    std::uint8_t corrector,
    StateCorrectionContract contract,
    double sealed_alpha,
    ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  Impl& impl = *implementation_;
  impl.current_pressure_work = {};
  if (corrector == 2U) {
    impl.terminal = {};
    impl.corrected_pending = {};
  }
  const Int3 cells = impl.cells;
  const ConstFieldView correction_read = as_const(correction);
  const ConstFieldView velocity = as_const(state.velocity);
  const ConstFieldView perturbation = as_const(state.pressure_perturbation);
  const ConstFieldView density = as_const(state.density);
  const ConstFaceFluxView phi_h_by_a = as_const(impl.workspace.phi_h_by_a);
  const ConstFaceFluxView output_flux = as_const(trial_flux);
  const bool valid =
      impl.pressure_boundary.current() && pressure.valid() &&
      impl.pressure_reference_kind ==
          PressureReferenceKind::boundary_absolute &&
      pressure.plan == impl.fingerprint &&
      (corrector == 1U || corrector == 2U) &&
      pressure.corrector == corrector && impl.current.valid() &&
      impl.current.corrector == corrector &&
      pressure.intermediate == impl.current.dependency &&
      pressure.thermophysical_boundary_semantics ==
          impl.current.thermophysical_boundary_semantics &&
      pressure.thermophysical_boundary_target ==
          impl.current.thermophysical_boundary_target &&
      pressure.thermophysical_boundary_rank_local_binding ==
          impl.current.thermophysical_boundary_rank_local_binding &&
      pressure.thermophysical_boundary_collective_lineage ==
          impl.current.thermophysical_boundary_collective_lineage &&
      pressure.thermophysical_boundary_rank_local_lineage ==
          impl.current.thermophysical_boundary_rank_local_lineage &&
      correction.field == impl.correction_field &&
      detail::valid_cell_view(correction_read, cells, 0U, 1U, 1U) &&
      detail::valid_cell_view(state.velocity, cells, 0U, 3U) &&
      detail::valid_cell_view(state.pressure_perturbation, cells, 0U, 1U) &&
      detail::valid_cell_view(state.density, cells, 0U, 1U) &&
      detail::valid_cell_view(state.drho_dp_h_y, cells, 0U, 1U, 0U) &&
      detail::valid_flux_view(trial_flux, cells) &&
      !trial_flux.certificate.valid() &&
      !detail::field_views_overlap(velocity, perturbation) &&
      !detail::field_views_overlap(velocity, density) &&
      !detail::field_views_overlap(perturbation, density) &&
      !detail::field_views_overlap(correction_read, velocity) &&
      !detail::field_views_overlap(correction_read, perturbation) &&
      !detail::field_views_overlap(correction_read, density) &&
      !detail::field_views_overlap(state.drho_dp_h_y, velocity) &&
      !detail::field_views_overlap(state.drho_dp_h_y, perturbation) &&
      !detail::field_views_overlap(state.drho_dp_h_y, density) &&
      !detail::field_views_overlap(
          correction_read, as_const(impl.workspace.pressure_gradient));
  if (!valid) {
    impl.current = {};
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const ConstFieldView cell_views[]{correction_read, velocity, perturbation,
                                    density, state.drho_dp_h_y};
  const ConstFaceFieldView output_faces[]{output_flux.x, output_flux.y,
                                          output_flux.z};
  for (ConstFieldView cell : cell_views) {
    for (ConstFaceFieldView face : output_faces) {
      if (detail::cell_face_views_overlap(cell, face)) {
        impl.current = {};
        return {StatusCode::invalid_plan, kPisoCoupler};
      }
    }
  }
  const ConstFaceFieldView internal_faces[]{
      phi_h_by_a.x, phi_h_by_a.y, phi_h_by_a.z,
      as_const(impl.workspace.x_pressure_coefficient),
      as_const(impl.workspace.y_pressure_coefficient),
      as_const(impl.workspace.z_pressure_coefficient)};
  for (ConstFaceFieldView output : output_faces) {
    for (ConstFaceFieldView internal : internal_faces) {
      if (detail::face_views_overlap(output, internal)) {
        impl.current = {};
        return {StatusCode::invalid_plan, kPisoCoupler};
      }
    }
  }

  double alpha = 1.0;
  double maximum_depletion = 0.0;
  Status status{};
  if (contract == StateCorrectionContract::pressure_sealed) {
    alpha = sealed_alpha;
  } else {
    double local_max_depletion = 0.0;
    Status local_status{};
    for (std::int32_t iz = 0; iz < cells.z && local_status; ++iz) {
      for (std::int32_t iy = 0; iy < cells.y && local_status; ++iy) {
        for (std::int32_t ix = 0; ix < cells.x; ++ix) {
          const Int3 cell{ix, iy, iz};
          const double rho = density.unchecked(cell, 0U);
          const double psi = state.drho_dp_h_y.unchecked(cell, 0U);
          const double delta = correction_read.unchecked(cell, 0U);
          if (!finite_positive(rho) || !std::isfinite(psi) || psi < 0.0 ||
              !std::isfinite(delta)) {
            local_status = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
          const double psi_delta = psi * delta;
          if (!std::isfinite(psi_delta)) {
            local_status = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
          const double depletion = std::max(0.0, -psi_delta / rho);
          if (!std::isfinite(depletion)) {
            local_status = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
          local_max_depletion = std::max(local_max_depletion, depletion);
        }
      }
    }
    status = reductions.checked_max(
        {&local_max_depletion, 1U}, {&maximum_depletion, 1U},
        local_status);
    if (!status) {
      impl.current = {};
      return status;
    }
    if (contract == StateCorrectionContract::pressure_unsealed &&
        maximum_depletion > kPressureCorrectionSafetyFraction) {
      alpha = std::nextafter(
          kPressureCorrectionSafetyFraction / maximum_depletion, 0.0);
    }
  }
  if (!std::isfinite(alpha) || !(alpha > 0.0) || alpha > 1.0) {
    impl.current = {};
    return {StatusCode::numerical_failure, kPisoNumerical};
  }
  if (alpha < 1.0) {
    for (std::int32_t iz = 0; iz < cells.z; ++iz) {
      for (std::int32_t iy = 0; iy < cells.y; ++iy) {
        for (std::int32_t ix = 0; ix < cells.x; ++ix) {
          const Int3 cell{ix, iy, iz};
          const double scaled = alpha * correction.unchecked(cell, 0U);
          if (!std::isfinite(scaled)) {
            impl.current = {};
            return {StatusCode::numerical_failure, kPisoNumerical};
          }
          correction.unchecked(cell, 0U) = scaled;
        }
      }
    }
  }

  std::array<FieldView, 1U> correction_views{correction};
  HaloTicket ticket;
  status = impl.correction_halo->begin(
      impl.correction_halo_stage,
      {correction_views.data(), correction_views.size()}, ticket);
  if (status) {
    status = impl.correction_halo->finish(
        ticket, {correction_views.data(), correction_views.size()});
  }
  if (!status) {
    impl.current = {};
    return status;
  }
  correction = correction_views[0U];
  status = impl.pressure_boundary.fill_ghosts(correction);
  status = reductions.consensus(status);
  if (!status) {
    impl.current = {};
    return status;
  }
  // Quadratic IBM reconstruction is a post-solve correction only. Its
  // compact donor gather must never enter a Krylov/MG operator application.
  if (impl.immersed_interface != nullptr &&
      impl.pressure_correction_donors != nullptr) {
    std::array<FieldView, 1U> donor_views{correction};
    status = impl.pressure_correction_donors->preflight_exchange(
        impl.pressure_correction_donor_stage,
        {donor_views.data(), donor_views.size()});
    status = reductions.consensus(status);
    if (status) {
      status = impl.pressure_correction_donors->exchange(
        impl.pressure_correction_donor_stage,
        {donor_views.data(), donor_views.size()});
    }
    status = reductions.consensus(status);
    if (!status) {
      impl.current = {};
      return status;
    }
    correction = donor_views[0U];
  }
  const ConstFieldView corrected_read = as_const(correction);
  const std::array<ConstFieldView, 1U> gradient_reads{corrected_read};
  const std::array<FieldView, 1U> gradient_writes{
      impl.workspace.pressure_gradient};
  const KernelBox box{{0, 0, 0}, cells};
  const KernelInvocation gradient_invocation{
      {gradient_reads.data(), gradient_reads.size()},
      {gradient_writes.data(), gradient_writes.size()}, box, 0U, 0U, 1U,
      0U, nullptr};
  status = cartesian_gradient(*impl.kernels, gradient_invocation);
  if (status && impl.immersed_interface != nullptr) {
    status = impl.immersed_interface->correct_pressure_gradient(
        corrected_read, impl.workspace.pressure_gradient);
  }
  status = reductions.consensus(status);
  if (!status) {
    impl.current = {};
    return status;
  }

  for (std::int32_t iz = 0; iz < cells.z; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        const double delta = corrected_read.unchecked(cell, 0U);
        const double next_pressure =
            perturbation.unchecked(cell, 0U) + delta;
        const double next_density =
            density.unchecked(cell, 0U) +
            state.drho_dp_h_y.unchecked(cell, 0U) * delta;
        if (!std::isfinite(delta) || !std::isfinite(next_pressure) ||
            !(next_density > 0.0) || !std::isfinite(next_density)) {
          impl.current = {};
          return {StatusCode::numerical_failure, kPisoNumerical};
        }
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double next_velocity =
              impl.workspace.h_by_a.unchecked(cell, component) -
              impl.workspace.r_au.unchecked(cell, component) *
                  impl.workspace.pressure_gradient.unchecked(cell, component);
          if (!std::isfinite(next_velocity)) {
            impl.current = {};
            return {StatusCode::numerical_failure, kPisoNumerical};
          }
        }
      }
    }
  }
  const CartesianAxis axes[]{CartesianAxis::x, CartesianAxis::y,
                             CartesianAxis::z};
  for (CartesianAxis axis : axes) {
    const ConstFaceFieldView source = detail::select(phi_h_by_a, axis);
    const ConstFaceFieldView coefficient =
        axis == CartesianAxis::x
            ? as_const(impl.workspace.x_pressure_coefficient)
            : (axis == CartesianAxis::y
                   ? as_const(impl.workspace.y_pressure_coefficient)
                   : as_const(impl.workspace.z_pressure_coefficient));
    const FaceFieldView output = detail::select(trial_flux, axis);
    for (std::int32_t iz = 0; iz < output.extents.z; ++iz) {
      for (std::int32_t iy = 0; iy < output.extents.y; ++iy) {
        for (std::int32_t ix = 0; ix < output.extents.x; ++ix) {
          const Int3 face{ix, iy, iz};
          const double next =
              source.unchecked(face) +
              impl.pressure_boundary.mass_flux_response_unchecked(
                  corrected_read, axis, face, coefficient.unchecked(face));
          if (!std::isfinite(next)) {
            impl.current = {};
            return {StatusCode::numerical_failure, kPisoNumerical};
          }
        }
      }
    }
  }

  for (std::int32_t iz = 0; iz < cells.z; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        const double delta = corrected_read.unchecked(cell, 0U);
        state.pressure_perturbation.unchecked(cell, 0U) += delta;
        state.density.unchecked(cell, 0U) +=
            state.drho_dp_h_y.unchecked(cell, 0U) * delta;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          state.velocity.unchecked(cell, component) =
              impl.workspace.h_by_a.unchecked(cell, component) -
              impl.workspace.r_au.unchecked(cell, component) *
                  impl.workspace.pressure_gradient.unchecked(cell, component);
        }
      }
    }
  }
  for (CartesianAxis axis : axes) {
    const ConstFaceFieldView source = detail::select(phi_h_by_a, axis);
    const ConstFaceFieldView coefficient =
        axis == CartesianAxis::x
            ? as_const(impl.workspace.x_pressure_coefficient)
            : (axis == CartesianAxis::y
                   ? as_const(impl.workspace.y_pressure_coefficient)
                   : as_const(impl.workspace.z_pressure_coefficient));
    const FaceFieldView output = detail::select(trial_flux, axis);
    for (std::int32_t iz = 0; iz < output.extents.z; ++iz) {
      for (std::int32_t iy = 0; iy < output.extents.y; ++iy) {
        for (std::int32_t ix = 0; ix < output.extents.x; ++ix) {
          const Int3 face{ix, iy, iz};
          output.unchecked(face) =
              source.unchecked(face) +
              impl.pressure_boundary.mass_flux_response_unchecked(
                  corrected_read, axis, face, coefficient.unchecked(face));
        }
      }
    }
  }

  if (impl.immersed_interface != nullptr) {
    status = impl.immersed_interface->constrain_corrected_state(
        state.velocity, trial_flux);
    if (!status) {
      impl.current = {};
      return status;
    }
  }

  std::uint64_t state_hash = kFnvOffset;
  state_hash = hash_mix(state_hash, pressure.state);
  state_hash = mix_view(state_hash, corrected_read);
  state_hash = mix_view(state_hash, as_const(state.velocity));
  state_hash = mix_view(state_hash, as_const(state.pressure_perturbation));
  state_hash = mix_view(state_hash, as_const(state.density));
  state_hash = mix_face(state_hash, output_flux.x);
  state_hash = mix_face(state_hash, output_flux.y);
  state_hash = mix_face(state_hash, output_flux.z);
  state_hash = hash_mix(
      state_hash, pressure.thermophysical_boundary_semantics);
  state_hash =
      hash_mix(state_hash, pressure.thermophysical_boundary_target);
  state_hash = hash_mix(
      state_hash, pressure.thermophysical_boundary_rank_local_binding);
  state_hash = hash_mix(
      state_hash, pressure.thermophysical_boundary_collective_lineage);
  state_hash = hash_mix(
      state_hash, pressure.thermophysical_boundary_rank_local_lineage);
  state_hash = hash_mix(
      state_hash, impl.current_pressure_reference.pressure_reference);
  PisoStateCorrectionCertificate candidate;
  candidate.plan = impl.fingerprint;
  candidate.pressure_system = pressure.state;
  candidate.correction = corrected_read.revision;
  candidate.velocity = state.velocity.revision;
  candidate.pressure = state.pressure_perturbation.revision;
  candidate.density = state.density.revision;
  candidate.face_flux = trial_flux.revision;
  candidate.state = state_hash == 0U ? 1U : state_hash;
  candidate.corrector = corrector;
  candidate.closure = PisoStateClosure::pressure_affine;
  candidate.thermophysical_boundary_semantics =
      pressure.thermophysical_boundary_semantics;
  candidate.thermophysical_boundary_target =
      pressure.thermophysical_boundary_target;
  candidate.thermophysical_boundary_rank_local_binding =
      pressure.thermophysical_boundary_rank_local_binding;
  candidate.thermophysical_boundary_collective_lineage =
      pressure.thermophysical_boundary_collective_lineage;
  candidate.thermophysical_boundary_rank_local_lineage =
      pressure.thermophysical_boundary_rank_local_lineage;
  candidate.input_pressure_reference = impl.current_pressure_reference;
  candidate.output_pressure_reference = impl.current_pressure_reference;
  if (corrector == 1U) {
    impl.predecessor_c1.valid = true;
    impl.predecessor_c1.alpha = alpha;
    impl.predecessor_c1.maximum_depletion = maximum_depletion;
    impl.corrected_c1 = candidate;
  }
  impl.current = {};
  impl.current_pressure_reference = {};
  impl.current_absolute_pressure_reference = 0.0;
  certificate = candidate;
  return {};
}

Status PressureVelocityCoupler::correct_exact_coupled_state_impl(
    const PressureCorrectionCertificate& pressure,
    FieldView pressure_correction,
    ConstFieldView enthalpy_correction,
    PisoCoupledStateView state,
    PisoExactThermodynamicCandidateView exact_candidate,
    FaceFluxView output_flux,
    std::uint8_t corrector,
    ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  Impl& impl = *implementation_;
  impl.current_pressure_work = {};
  if (corrector == 2U) {
    impl.terminal = {};
    impl.corrected_pending = {};
  }
  const Int3 cells = impl.cells;
  const ConstFieldView dp = as_const(pressure_correction);
  const ConstFieldView velocity = as_const(state.velocity);
  const ConstFieldView perturbation = as_const(state.pressure_perturbation);
  const ConstFieldView enthalpy = as_const(state.enthalpy);
  const ConstFieldView density = as_const(state.density);
  const ConstFieldView temperature = as_const(state.temperature);
  const ConstFaceFluxView phi_h_by_a = as_const(impl.workspace.phi_h_by_a);
  const ConstFaceFluxView output_read = as_const(output_flux);
  const PisoExactEosClosureIdentity& closure = exact_candidate.closure;
  const PressureReferenceCertificate& reference =
      impl.current_pressure_reference;
  const bool closed = impl.pressure_reference_kind ==
                      PressureReferenceKind::closed_mass;
  const PressureEnergyCellActivity gauge_activity{
      impl.continuity_activity.cells,
      impl.continuity_activity.local_fingerprint,
      impl.continuity_activity.collective_fingerprint};
  const ClosedGaugeCorrectionPrepareInput gauge_input{
      reference,
      impl.current_absolute_pressure_reference,
      corrector,
      pressure.time,
      pressure.geometry,
      pressure.state,
      closure.closure,
      perturbation,
      dp,
      exact_candidate.pressure_compressibility,
      gauge_activity};
  const bool matching_gauge =
      closed
          ? impl.pressure_reference != nullptr &&
                impl.pressure_reference->matches_closed_gauge_correction(
                    gauge_input, exact_candidate.closed_gauge)
          : empty_field_view(exact_candidate.pressure_compressibility) &&
                empty_closed_gauge_certificate(exact_candidate.closed_gauge);
  const double gauge_shift =
      closed ? exact_candidate.closed_gauge.shift : 0.0;
  const double next_absolute_pressure_reference =
      closed ? exact_candidate.closed_gauge.next_pressure_reference
             : impl.current_absolute_pressure_reference;
  const PressureReferenceCertificate& output_reference =
      closed ? exact_candidate.closed_gauge.output_pressure_reference
             : reference;

  bool composition_valid =
      impl.thermodynamics_plan != nullptr &&
      exact_candidate.independent_species.size ==
          impl.independent_species_count &&
      (exact_candidate.independent_species.size == 0U ||
       exact_candidate.independent_species.data != nullptr);
  for (std::size_t species = 0U;
       species < exact_candidate.independent_species.size && composition_valid;
       ++species) {
    composition_valid =
        detail::valid_cell_view(
            exact_candidate.independent_species.data[species], cells, 0U, 1U,
            0U) &&
        impl.independent_species_semantic_fields != nullptr &&
        exact_candidate.independent_species.data[species].field ==
            impl.independent_species_semantic_fields[species];
  }
  const PlanFingerprint composition_numeric =
      composition_valid
          ? candidate_composition_numeric_local_provenance(
                exact_candidate.independent_species, cells)
          : PlanFingerprint{};
  const PlanFingerprint composition_binding =
      composition_valid
          ? semantic_composition_rank_local_binding(
                exact_candidate.independent_species)
          : PlanFingerprint{};
  const PlanFingerprint expected_composition_identity =
      composition_valid
          ? exact_composition_identity_local(
                impl.thermodynamics,
                exact_candidate.independent_species, cells)
          : PlanFingerprint{};

  const bool matching_closure =
      closure.valid() && reference.valid() &&
      closure.thermodynamics == reference.thermodynamics &&
      closure.pressure_reference == reference.pressure_reference &&
      closure.composition == expected_composition_identity &&
      same_piso_field_revision_identity(closure.pressure_state,
                                        perturbation) &&
      same_piso_field_revision_identity(closure.pressure_correction, dp) &&
      same_piso_field_revision_identity(closure.enthalpy_state, enthalpy) &&
      same_piso_field_revision_identity(closure.enthalpy_correction,
                                        enthalpy_correction) &&
      same_piso_field_revision_identity(closure.candidate_enthalpy,
                                        exact_candidate.enthalpy) &&
      same_piso_field_revision_identity(closure.candidate_density,
                                        exact_candidate.density) &&
      same_piso_field_revision_identity(closure.candidate_temperature,
                                        exact_candidate.temperature);
  const bool matching_c1_output =
      corrector == 1U ||
      (impl.current_corrected_c1.valid() &&
       (impl.current_corrected_c1.corrector == 1U ||
        impl.current_corrected_c1.corrector == 2U) &&
       same_pressure_reference_certificate(
           reference,
           impl.current_corrected_c1.output_pressure_reference));
  const bool valid =
      impl.pressure_boundary.current() && pressure.valid() &&
      pressure.plan == impl.fingerprint &&
      (corrector == 1U || corrector == 2U) &&
      pressure.corrector == corrector && impl.current.valid() &&
      impl.current.corrector == corrector &&
      pressure.intermediate == impl.current.dependency &&
      pressure.thermophysical_boundary_semantics ==
          impl.current.thermophysical_boundary_semantics &&
      pressure.thermophysical_boundary_target ==
          impl.current.thermophysical_boundary_target &&
      pressure.thermophysical_boundary_rank_local_binding ==
          impl.current.thermophysical_boundary_rank_local_binding &&
      pressure.thermophysical_boundary_collective_lineage ==
          impl.current.thermophysical_boundary_collective_lineage &&
      pressure.thermophysical_boundary_rank_local_lineage ==
          impl.current.thermophysical_boundary_rank_local_lineage &&
      pressure_correction.field == impl.correction_field &&
      detail::valid_cell_view(dp, cells, 0U, 1U, 1U) &&
      detail::valid_cell_view(enthalpy_correction, cells, 0U, 1U, 0U) &&
      detail::valid_cell_view(state.velocity, cells, 0U, 3U) &&
      detail::valid_cell_view(state.pressure_perturbation, cells, 0U, 1U) &&
      detail::valid_cell_view(state.enthalpy, cells, 0U, 1U) &&
      detail::valid_cell_view(state.density, cells, 0U, 1U) &&
      detail::valid_cell_view(state.temperature, cells, 0U, 1U) &&
      detail::valid_cell_view(exact_candidate.enthalpy, cells, 0U, 1U, 0U) &&
      detail::valid_cell_view(exact_candidate.density, cells, 0U, 1U, 0U) &&
      detail::valid_cell_view(exact_candidate.temperature, cells, 0U, 1U,
                              0U) &&
      (closed
           ? detail::valid_cell_view(
                 exact_candidate.pressure_compressibility, cells, 0U, 1U,
                 0U)
           : empty_field_view(exact_candidate.pressure_compressibility)) &&
      detail::valid_flux_view(output_flux, cells) &&
      !output_flux.certificate.valid() && matching_closure && matching_gauge &&
      matching_c1_output && composition_valid && composition_numeric != 0U &&
      composition_binding != 0U;
  Status local = valid ? Status{}
                       : Status{StatusCode::invalid_plan, kPisoCoupler};

  const std::array<ConstFieldView, 5U> writable_cells{
      velocity, perturbation, enthalpy, density, temperature};
  const std::array<ConstFieldView, 6U> immutable_cells{
      dp, enthalpy_correction, exact_candidate.enthalpy,
      exact_candidate.density, exact_candidate.temperature,
      exact_candidate.pressure_compressibility};
  const std::array<ConstFieldView, 3U> internal_cells{
      as_const(impl.workspace.r_au), as_const(impl.workspace.h_by_a),
      as_const(impl.workspace.pressure_gradient)};
  bool aliases = false;
  for (std::size_t left = 0U; local && left < writable_cells.size(); ++left) {
    for (std::size_t right = left + 1U; right < writable_cells.size();
         ++right) {
      aliases = aliases || detail::field_views_overlap(
                               writable_cells[left], writable_cells[right]);
    }
  }
  for (std::size_t left = 0U; local && left < immutable_cells.size(); ++left) {
    if (empty_field_view(immutable_cells[left])) continue;
    for (std::size_t right = left + 1U; right < immutable_cells.size();
         ++right) {
      if (empty_field_view(immutable_cells[right])) continue;
      aliases = aliases || detail::field_views_overlap(
                               immutable_cells[left], immutable_cells[right]);
    }
  }
  for (ConstFieldView writable : writable_cells) {
    if (!local) break;
    for (ConstFieldView immutable : immutable_cells) {
      if (empty_field_view(immutable)) continue;
      aliases = aliases ||
                detail::field_views_overlap(writable, immutable);
    }
    for (ConstFieldView internal : internal_cells) {
      aliases = aliases || detail::field_views_overlap(writable, internal);
    }
  }
  for (std::size_t species = 0U;
       local && species < exact_candidate.independent_species.size;
       ++species) {
    const ConstFieldView composition =
        exact_candidate.independent_species.data[species];
    for (ConstFieldView writable : writable_cells)
      aliases = aliases || detail::field_views_overlap(composition, writable);
    for (ConstFieldView immutable : immutable_cells) {
      if (!empty_field_view(immutable))
        aliases = aliases ||
                  detail::field_views_overlap(composition, immutable);
    }
    for (std::size_t prior = 0U; prior < species; ++prior)
      aliases = aliases || detail::field_views_overlap(
                               composition,
                               exact_candidate.independent_species
                                   .data[prior]);
  }
  for (ConstFieldView immutable : immutable_cells) {
    if (empty_field_view(immutable)) continue;
    if (!local) break;
    for (ConstFieldView internal : internal_cells) {
      aliases = aliases || detail::field_views_overlap(immutable, internal);
    }
  }
  const std::array<ConstFaceFieldView, 3U> output_faces{
      output_read.x, output_read.y, output_read.z};
  for (ConstFieldView cell : writable_cells) {
    if (!local) break;
    for (ConstFaceFieldView face : output_faces) {
      aliases = aliases || detail::cell_face_views_overlap(cell, face);
    }
  }
  for (ConstFieldView cell : immutable_cells) {
    if (empty_field_view(cell)) continue;
    if (!local) break;
    for (ConstFaceFieldView face : output_faces) {
      aliases = aliases || detail::cell_face_views_overlap(cell, face);
    }
  }
  for (std::size_t left = 0U; local && left < output_faces.size(); ++left) {
    for (std::size_t right = left + 1U; right < output_faces.size(); ++right) {
      aliases = aliases ||
                detail::face_views_overlap(output_faces[left],
                                           output_faces[right]);
    }
  }
  const std::array<ConstFaceFieldView, 6U> internal_faces{
      phi_h_by_a.x, phi_h_by_a.y, phi_h_by_a.z,
      as_const(impl.workspace.x_pressure_coefficient),
      as_const(impl.workspace.y_pressure_coefficient),
      as_const(impl.workspace.z_pressure_coefficient)};
  for (ConstFieldView cell : writable_cells) {
    if (!local) break;
    for (ConstFaceFieldView internal : internal_faces) {
      aliases = aliases ||
                detail::cell_face_views_overlap(cell, internal);
    }
  }
  for (ConstFieldView cell : immutable_cells) {
    if (empty_field_view(cell)) continue;
    if (!local) break;
    for (ConstFaceFieldView internal : internal_faces) {
      aliases = aliases ||
                detail::cell_face_views_overlap(cell, internal);
    }
  }
  for (ConstFaceFieldView output : output_faces) {
    if (!local) break;
    for (ConstFieldView internal : internal_cells) {
      aliases = aliases ||
                detail::cell_face_views_overlap(internal, output);
    }
    for (ConstFaceFieldView internal : internal_faces) {
      aliases = aliases || detail::face_views_overlap(output, internal);
    }
  }
  if (local && aliases)
    local = {StatusCode::invalid_plan, kPisoCoupler};

  // Validate the caller's complete thermodynamic candidate collectively
  // before pressure-correction halo/ghost storage or any state/flux output is
  // touched.  In particular, no affine rho+rho_p*dp positivity decision is
  // made on this exact-EOS path.
  for (std::int32_t iz = 0; iz < cells.z && local; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y && local; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        const double old_pressure = perturbation.unchecked(cell, 0U);
        const double delta_pressure = dp.unchecked(cell, 0U);
        const double old_enthalpy = enthalpy.unchecked(cell, 0U);
        const double delta_enthalpy =
            enthalpy_correction.unchecked(cell, 0U);
        const double candidate_enthalpy =
            exact_candidate.enthalpy.unchecked(cell, 0U);
        const double candidate_density =
            exact_candidate.density.unchecked(cell, 0U);
        const double candidate_temperature =
            exact_candidate.temperature.unchecked(cell, 0U);
        const double expected_enthalpy = old_enthalpy + delta_enthalpy;
        const double enthalpy_scale =
            std::max({1.0, std::abs(old_enthalpy),
                      std::abs(delta_enthalpy),
                      std::abs(candidate_enthalpy)});
        const double closure_error =
            std::abs(candidate_enthalpy - expected_enthalpy);
        const double next_pressure = old_pressure + delta_pressure;
        const double next_absolute_pressure =
            impl.current_absolute_pressure_reference + next_pressure;
        const double represented_pressure = next_pressure - gauge_shift;
        const double represented_absolute_pressure =
            next_absolute_pressure_reference + represented_pressure;
        const double pressure_scale =
            std::max({1.0, std::abs(next_absolute_pressure),
                      std::abs(represented_absolute_pressure)});
        const double representation_error =
            std::abs(represented_absolute_pressure - next_absolute_pressure);
        if (!std::isfinite(old_pressure) ||
            !std::isfinite(delta_pressure) ||
            !std::isfinite(next_pressure) ||
            !std::isfinite(next_absolute_pressure) ||
            !(next_absolute_pressure > 0.0) ||
            !std::isfinite(represented_pressure) ||
            !std::isfinite(represented_absolute_pressure) ||
            !(represented_absolute_pressure > 0.0) ||
            !std::isfinite(representation_error) ||
            representation_error >
                16.0 * std::numeric_limits<double>::epsilon() *
                    pressure_scale ||
            !std::isfinite(old_enthalpy) ||
            !std::isfinite(delta_enthalpy) ||
            !std::isfinite(candidate_enthalpy) ||
            !std::isfinite(candidate_density) ||
            !(candidate_density > 0.0) ||
            !std::isfinite(candidate_temperature) ||
            !(candidate_temperature > 0.0) ||
            !std::isfinite(closure_error) ||
            closure_error >
                16.0 * std::numeric_limits<double>::epsilon() *
                    enthalpy_scale ||
            !std::isfinite(density.unchecked(cell, 0U)) ||
            !(density.unchecked(cell, 0U) > 0.0) ||
            !std::isfinite(temperature.unchecked(cell, 0U)) ||
            !(temperature.unchecked(cell, 0U) > 0.0)) {
          local = {StatusCode::numerical_failure, kPisoNumerical};
          break;
        }
      }
    }
  }
  double local_guard = 0.0;
  double global_guard = 0.0;
  Status status = reductions.checked_max(
      {&local_guard, 1U}, {&global_guard, 1U}, local);
  if (!status) {
    impl.current = {};
    impl.current_pressure_reference = {};
    impl.current_absolute_pressure_reference = 0.0;
    return status;
  }

  std::array<FieldView, 1U> correction_views{pressure_correction};
  HaloTicket ticket;
  status = impl.correction_halo->begin(
      impl.correction_halo_stage,
      {correction_views.data(), correction_views.size()}, ticket);
  if (status) {
    status = impl.correction_halo->finish(
        ticket, {correction_views.data(), correction_views.size()});
  }
  if (!status) {
    impl.current = {};
    impl.current_pressure_reference = {};
    impl.current_absolute_pressure_reference = 0.0;
    return status;
  }
  pressure_correction = correction_views[0U];
  status = impl.pressure_boundary.fill_ghosts(pressure_correction);
  status = reductions.consensus(status);
  if (!status) {
    impl.current = {};
    impl.current_pressure_reference = {};
    impl.current_absolute_pressure_reference = 0.0;
    return status;
  }
  if (impl.immersed_interface != nullptr &&
      impl.pressure_correction_donors != nullptr) {
    std::array<FieldView, 1U> donor_views{pressure_correction};
    status = impl.pressure_correction_donors->preflight_exchange(
        impl.pressure_correction_donor_stage,
        {donor_views.data(), donor_views.size()});
    status = reductions.consensus(status);
    if (status) {
      status = impl.pressure_correction_donors->exchange(
          impl.pressure_correction_donor_stage,
          {donor_views.data(), donor_views.size()});
    }
    status = reductions.consensus(status);
    if (!status) {
      impl.current = {};
      impl.current_pressure_reference = {};
      impl.current_absolute_pressure_reference = 0.0;
      return status;
    }
    pressure_correction = donor_views[0U];
  }
  const ConstFieldView corrected_dp = as_const(pressure_correction);
  const std::array<ConstFieldView, 1U> gradient_reads{corrected_dp};
  const std::array<FieldView, 1U> gradient_writes{
      impl.workspace.pressure_gradient};
  const KernelInvocation gradient_invocation{
      {gradient_reads.data(), gradient_reads.size()},
      {gradient_writes.data(), gradient_writes.size()},
      {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U, nullptr};
  status = cartesian_gradient(*impl.kernels, gradient_invocation);
  if (status && impl.immersed_interface != nullptr) {
    status = impl.immersed_interface->correct_pressure_gradient(
        corrected_dp, impl.workspace.pressure_gradient);
  }
  status = reductions.consensus(status);
  if (!status) {
    impl.current = {};
    impl.current_pressure_reference = {};
    impl.current_absolute_pressure_reference = 0.0;
    return status;
  }

  local = {};
  for (std::int32_t iz = 0; iz < cells.z && local; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y && local; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double next_velocity =
              impl.workspace.h_by_a.unchecked(cell, component) -
              impl.workspace.r_au.unchecked(cell, component) *
                  impl.workspace.pressure_gradient.unchecked(cell,
                                                              component);
          if (!std::isfinite(next_velocity)) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
        }
        if (!local) break;
      }
    }
  }
  const CartesianAxis axes[]{CartesianAxis::x, CartesianAxis::y,
                             CartesianAxis::z};
  for (CartesianAxis axis : axes) {
    if (!local) break;
    const ConstFaceFieldView source = detail::select(phi_h_by_a, axis);
    const ConstFaceFieldView coefficient =
        axis == CartesianAxis::x
            ? as_const(impl.workspace.x_pressure_coefficient)
            : (axis == CartesianAxis::y
                   ? as_const(impl.workspace.y_pressure_coefficient)
                   : as_const(impl.workspace.z_pressure_coefficient));
    const FaceFieldView output = detail::select(output_flux, axis);
    for (std::int32_t iz = 0; iz < output.extents.z && local; ++iz) {
      for (std::int32_t iy = 0; iy < output.extents.y && local; ++iy) {
        for (std::int32_t ix = 0; ix < output.extents.x; ++ix) {
          const Int3 face{ix, iy, iz};
          const double next =
              source.unchecked(face) +
              impl.pressure_boundary.mass_flux_response_unchecked(
                  corrected_dp, axis, face, coefficient.unchecked(face));
          if (!std::isfinite(next)) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
        }
      }
    }
  }
  status = reductions.checked_max(
      {&local_guard, 1U}, {&global_guard, 1U}, local);
  if (!status) {
    impl.current = {};
    impl.current_pressure_reference = {};
    impl.current_absolute_pressure_reference = 0.0;
    return status;
  }

  const bool staged_ibm = impl.immersed_interface != nullptr;
  if (staged_ibm) {
    // The IBM candidate is staged in coupler-owned storage so a failed mask
    // cannot expose a partial caller state.  Cartesian keeps its predictor
    // workspace immutable and needs no additional collective.
    for (std::int32_t iz = 0; iz < cells.z; ++iz) {
      for (std::int32_t iy = 0; iy < cells.y; ++iy) {
        for (std::int32_t ix = 0; ix < cells.x; ++ix) {
          const Int3 cell{ix, iy, iz};
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            impl.workspace.h_by_a.unchecked(cell, component) =
                impl.workspace.h_by_a.unchecked(cell, component) -
                impl.workspace.r_au.unchecked(cell, component) *
                    impl.workspace.pressure_gradient.unchecked(cell,
                                                                component);
          }
        }
      }
    }
    for (CartesianAxis axis : axes) {
      const ConstFaceFieldView source = detail::select(phi_h_by_a, axis);
      const ConstFaceFieldView coefficient =
          axis == CartesianAxis::x
              ? as_const(impl.workspace.x_pressure_coefficient)
              : (axis == CartesianAxis::y
                     ? as_const(impl.workspace.y_pressure_coefficient)
                     : as_const(impl.workspace.z_pressure_coefficient));
      const FaceFieldView staged =
          detail::select(impl.workspace.phi_h_by_a, axis);
      for (std::int32_t iz = 0; iz < staged.extents.z; ++iz) {
        for (std::int32_t iy = 0; iy < staged.extents.y; ++iy) {
          for (std::int32_t ix = 0; ix < staged.extents.x; ++ix) {
            const Int3 face{ix, iy, iz};
            staged.unchecked(face) =
                source.unchecked(face) +
                impl.pressure_boundary.mass_flux_response_unchecked(
                    corrected_dp, axis, face, coefficient.unchecked(face));
          }
        }
      }
    }
    const Status staged_local =
        impl.immersed_interface->constrain_corrected_state(
            impl.workspace.h_by_a, impl.workspace.phi_h_by_a);
    status = reductions.checked_max(
        {&local_guard, 1U}, {&global_guard, 1U}, staged_local);
    if (!status) {
      impl.current = {};
      impl.current_corrected_c1 = {};
      impl.current_pressure_reference = {};
      impl.current_absolute_pressure_reference = 0.0;
      return status;
    }
  }

  // EOS closure belongs to the same target layer as the mechanical
  // correction.  Replay it only after correction ghosts, donor exchange,
  // Cartesian gradient, and IBM constraints have formed the candidate U.
  // The earlier validation intentionally checks only caller scratch and
  // positivity, so a base-layer velocity can never stand in for this U.
  local = {};
  for (std::int32_t iz = 0; iz < cells.z && local; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y && local; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        const double represented_pressure =
            perturbation.unchecked(cell, 0U) +
            corrected_dp.unchecked(cell, 0U) - gauge_shift;
        const double candidate_enthalpy =
            exact_candidate.enthalpy.unchecked(cell, 0U);
        const double candidate_density =
            exact_candidate.density.unchecked(cell, 0U);
        const double candidate_temperature =
            exact_candidate.temperature.unchecked(cell, 0U);
        for (std::size_t species = 0U;
             species < exact_candidate.independent_species.size; ++species) {
          impl.frozen_candidate_composition_values[species] =
              exact_candidate.independent_species.data[species].unchecked(
                  cell, 0U);
        }
        Real3 replay_velocity;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double corrected_velocity =
              staged_ibm
                  ? impl.workspace.h_by_a.unchecked(cell, component)
                  : impl.workspace.h_by_a.unchecked(cell, component) -
                        impl.workspace.r_au.unchecked(cell, component) *
                            impl.workspace.pressure_gradient.unchecked(
                                cell, component);
          if (component == 0U)
            replay_velocity.x = corrected_velocity;
          else if (component == 1U)
            replay_velocity.y = corrected_velocity;
          else
            replay_velocity.z = corrected_velocity;
        }
        ThermoState replayed_thermo;
        const Status replayed =
            impl.thermodynamics_plan->evaluate_from_reference_pressure(
                next_absolute_pressure_reference, represented_pressure,
                candidate_enthalpy,
                {impl.frozen_candidate_composition_values.get(),
                 impl.independent_species_count},
                replay_velocity, replayed_thermo, candidate_temperature);
        const double density_scale =
            std::max({1.0, std::abs(candidate_density),
                      std::abs(replayed_thermo.rho)});
        const double temperature_scale =
            std::max({1.0, std::abs(candidate_temperature),
                      std::abs(replayed_thermo.temperature)});
        constexpr double exact_eos_factor =
            64.0 * std::numeric_limits<double>::epsilon();
        bool exact_eos =
            replayed && finite_positive(replayed_thermo.rho) &&
            finite_positive(replayed_thermo.temperature) &&
            finite_positive(replayed_thermo.drho_dp_hY) &&
            std::abs(candidate_density - replayed_thermo.rho) <=
                exact_eos_factor * density_scale &&
            std::abs(candidate_temperature - replayed_thermo.temperature) <=
                exact_eos_factor * temperature_scale;
        if (closed) {
          const double selected_compressibility =
              exact_candidate.pressure_compressibility.unchecked(cell, 0U);
          const double compressibility_scale =
              std::max({1.0, std::abs(selected_compressibility),
                        std::abs(replayed_thermo.drho_dp_hY)});
          exact_eos =
              exact_eos && finite_positive(selected_compressibility) &&
              std::abs(selected_compressibility -
                       replayed_thermo.drho_dp_hY) <=
                  exact_eos_factor * compressibility_scale;
        }
        if (!exact_eos) {
          local = {StatusCode::invalid_plan, kPisoCoupler};
          break;
        }
      }
    }
  }
  status = reductions.checked_max(
      {&local_guard, 1U}, {&global_guard, 1U}, local);
  if (!status) {
    impl.current = {};
    impl.current_corrected_c1 = {};
    impl.current_pressure_reference = {};
    impl.current_absolute_pressure_reference = 0.0;
    return status;
  }

  // Every possibly failing operation has passed a collective gate.  The
  // following caller-visible commit only copies the staged U/flux candidate
  // alongside p/h/rho/T and cannot fail.
  for (std::int32_t iz = 0; iz < cells.z; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        state.pressure_perturbation.unchecked(cell, 0U) =
            (state.pressure_perturbation.unchecked(cell, 0U) +
             corrected_dp.unchecked(cell, 0U)) -
            gauge_shift;
        state.enthalpy.unchecked(cell, 0U) =
            exact_candidate.enthalpy.unchecked(cell, 0U);
        state.density.unchecked(cell, 0U) =
            exact_candidate.density.unchecked(cell, 0U);
        state.temperature.unchecked(cell, 0U) =
            exact_candidate.temperature.unchecked(cell, 0U);
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          state.velocity.unchecked(cell, component) =
              staged_ibm
                  ? impl.workspace.h_by_a.unchecked(cell, component)
                  : impl.workspace.h_by_a.unchecked(cell, component) -
                        impl.workspace.r_au.unchecked(cell, component) *
                            impl.workspace.pressure_gradient.unchecked(
                                cell, component);
        }
      }
    }
  }
  for (CartesianAxis axis : axes) {
    const ConstFaceFieldView source = detail::select(phi_h_by_a, axis);
    const ConstFaceFieldView coefficient =
        axis == CartesianAxis::x
            ? as_const(impl.workspace.x_pressure_coefficient)
            : (axis == CartesianAxis::y
                   ? as_const(impl.workspace.y_pressure_coefficient)
                   : as_const(impl.workspace.z_pressure_coefficient));
    const FaceFieldView output = detail::select(output_flux, axis);
    for (std::int32_t iz = 0; iz < output.extents.z; ++iz) {
      for (std::int32_t iy = 0; iy < output.extents.y; ++iy) {
        for (std::int32_t ix = 0; ix < output.extents.x; ++ix) {
          const Int3 face{ix, iy, iz};
          output.unchecked(face) =
              staged_ibm
                  ? source.unchecked(face)
                  : source.unchecked(face) +
                        impl.pressure_boundary.mass_flux_response_unchecked(
                            corrected_dp, axis, face,
                            coefficient.unchecked(face));
        }
      }
    }
  }

  const PlanFingerprint closure_fingerprint =
      exact_eos_closure_fingerprint(closure);
  std::uint64_t state_hash = kFnvOffset;
  state_hash = hash_mix(state_hash, pressure.state);
  state_hash = mix_view(state_hash, corrected_dp);
  state_hash = mix_view(state_hash, enthalpy_correction);
  state_hash = mix_view(state_hash, as_const(state.velocity));
  state_hash = mix_view(state_hash, as_const(state.pressure_perturbation));
  state_hash = mix_view(state_hash, as_const(state.enthalpy));
  state_hash = mix_view(state_hash, as_const(state.density));
  state_hash = mix_view(state_hash, as_const(state.temperature));
  state_hash = hash_mix(state_hash, closure_fingerprint);
  state_hash = mix_face(state_hash, output_read.x);
  state_hash = mix_face(state_hash, output_read.y);
  state_hash = mix_face(state_hash, output_read.z);
  state_hash = hash_mix(
      state_hash, pressure.thermophysical_boundary_semantics);
  state_hash =
      hash_mix(state_hash, pressure.thermophysical_boundary_target);
  state_hash = hash_mix(
      state_hash, pressure.thermophysical_boundary_rank_local_binding);
  state_hash = hash_mix(
      state_hash, pressure.thermophysical_boundary_collective_lineage);
  state_hash = hash_mix(
      state_hash, pressure.thermophysical_boundary_rank_local_lineage);
  state_hash = hash_mix(state_hash, reference.pressure_reference);
  state_hash = hash_mix(state_hash, output_reference.pressure_reference);
  state_hash = hash_mix(
      state_hash, exact_candidate.closed_gauge.collective_transaction);
  state_hash = hash_mix(
      state_hash, exact_candidate.closed_gauge.rank_local_transaction);
  state_hash = hash_mix(state_hash, double_bits(gauge_shift));
  state_hash = hash_mix(state_hash,
                        double_bits(next_absolute_pressure_reference));
  state_hash = hash_mix(state_hash, composition_numeric);
  state_hash = hash_mix(
      state_hash,
      static_cast<std::uint64_t>(exact_candidate.independent_species.size));
  PisoStateCorrectionCertificate committed;
  committed.plan = impl.fingerprint;
  committed.pressure_system = pressure.state;
  committed.correction = corrected_dp.revision;
  committed.enthalpy_correction = enthalpy_correction.revision;
  committed.velocity = state.velocity.revision;
  committed.pressure = state.pressure_perturbation.revision;
  committed.enthalpy = state.enthalpy.revision;
  committed.density = state.density.revision;
  committed.temperature = state.temperature.revision;
  committed.face_flux = output_flux.revision;
  committed.exact_eos_closure = closure_fingerprint;
  committed.state = state_hash == 0U ? 1U : state_hash;
  committed.corrector = corrector;
  committed.closure = PisoStateClosure::exact_eos;
  committed.thermophysical_boundary_semantics =
      pressure.thermophysical_boundary_semantics;
  committed.thermophysical_boundary_target =
      pressure.thermophysical_boundary_target;
  committed.thermophysical_boundary_rank_local_binding =
      pressure.thermophysical_boundary_rank_local_binding;
  committed.thermophysical_boundary_collective_lineage =
      pressure.thermophysical_boundary_collective_lineage;
  committed.thermophysical_boundary_rank_local_lineage =
      pressure.thermophysical_boundary_rank_local_lineage;
  committed.input_pressure_reference = reference;
  committed.output_pressure_reference = output_reference;
  committed.closed_gauge_collective_transaction =
      exact_candidate.closed_gauge.collective_transaction;
  committed.closed_gauge_rank_local_transaction =
      exact_candidate.closed_gauge.rank_local_transaction;
  committed.composition_numeric_provenance = composition_numeric;
  committed.composition_rank_local_binding = composition_binding;
  committed.independent_species_count =
      exact_candidate.independent_species.size;
  if (corrector == 1U) {
    impl.predecessor_c1.valid = true;
    impl.predecessor_c1.alpha = 1.0;
    impl.predecessor_c1.maximum_depletion = 0.0;
    impl.corrected_c1 = committed;
  }
  impl.current = {};
  impl.current_pressure_reference = {};
  impl.current_absolute_pressure_reference = 0.0;
  certificate = committed;
  return {};
}

Status PressureVelocityCoupler::correct_trial_state(
    const PressureCorrectionCertificate& pressure,
    FieldView correction,
    PisoTrialStateView state,
    FaceFluxView trial_flux,
    ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  if (pressure.corrector != 1U) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  return correct_state_impl(pressure, correction, state, trial_flux, 1U,
                            StateCorrectionContract::pressure_unsealed, 1.0,
                            reductions, certificate);
}

Status PressureVelocityCoupler::correct_coupled_trial_state(
    const PressureCorrectionCertificate& pressure,
    FieldView pressure_correction,
    ConstFieldView enthalpy_correction,
    PisoCoupledStateView state,
    PisoExactThermodynamicCandidateView candidate,
    FaceFluxView trial_flux,
    ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  if (implementation_ == nullptr || pressure.corrector != 1U) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  implementation_->sealed = {};
  return correct_exact_coupled_state_impl(
      pressure, pressure_correction, enthalpy_correction, state, candidate,
      trial_flux, 1U, reductions, certificate);
}

Status PressureVelocityCoupler::correct_pending_state(
    const PressureCorrectionCertificate& pressure,
    FieldView correction,
    PisoTrialStateView state,
    PendingFaceFluxView& pending_flux,
    ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  Impl& impl = *implementation_;
  impl.corrected_pending = {};
  if (pressure.corrector != 2U || !pending_flux.valid() ||
      pending_flux.writer_ == nullptr || pending_flux.storage_ == nullptr) {
    impl.sealed = {};
    impl.current = {};
    impl.current_corrected_c1 = {};
    impl.current_pressure_reference = {};
    impl.current_absolute_pressure_reference = 0.0;
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  bool has_sealed_alpha = false;
  double sealed_alpha = 1.0;
  if (impl.sealed.valid) {
    const LinearOperatorCertificate& exact = impl.sealed.exact_operator;
    const bool complete_exact =
        exact.identity.fingerprint != 0U &&
        exact.collective_fingerprint != 0U &&
        exact.local_shape.x > 0 && exact.local_shape.y > 0 &&
        exact.local_shape.z > 0;
    const bool matching =
        complete_exact &&
        same_pressure_certificate(impl.sealed.pressure, pressure) &&
        same_field_identity(impl.sealed.correction, as_const(correction)) &&
        std::isfinite(impl.sealed.alpha) && impl.sealed.alpha > 0.0 &&
        impl.sealed.alpha <= 1.0;
    if (!matching) {
      impl.sealed = {};
      impl.current = {};
      return {StatusCode::invalid_plan, kPisoCoupler};
    }
    has_sealed_alpha = true;
    sealed_alpha = impl.sealed.alpha;
    impl.sealed = {};
  }
  FaceFluxView output;
  output.x = pending_flux.x_;
  output.y = pending_flux.y_;
  output.z = pending_flux.z_;
  output.revision = pending_flux.revision_;
  const Status status = correct_state_impl(
      pressure, correction, state, output, 2U,
      has_sealed_alpha ? StateCorrectionContract::pressure_sealed
                       : StateCorrectionContract::pressure_unsealed,
      sealed_alpha, reductions, certificate);
  if (status) {
    impl.corrected_pending = {certificate,
                              &pending_flux,
                              pending_flux.writer_,
                              pending_flux.storage_,
                              pending_flux.revision_,
                              pending_flux.writer_identity_,
                              pending_flux.attempt_identity_};
    impl.corrected_pending.active_cells = impl.continuity_activity.cells;
    impl.corrected_pending.activity_local_fingerprint =
        impl.continuity_activity.local_fingerprint;
    impl.corrected_pending.activity_collective_fingerprint =
        impl.continuity_activity.collective_fingerprint;
  }
  return status;
}

Status PressureVelocityCoupler::correct_coupled_pending_state(
    const PressureCorrectionCertificate& pressure,
    FieldView pressure_correction,
    ConstFieldView enthalpy_correction,
    PisoCoupledStateView state,
    PisoExactThermodynamicCandidateView candidate,
    PendingFaceFluxView& pending_flux,
    ReductionEngine& reductions,
    PisoStateCorrectionCertificate& certificate) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  Impl& impl = *implementation_;
  impl.sealed = {};
  impl.corrected_pending = {};
  if (pressure.corrector != 2U || !pending_flux.valid() ||
      pending_flux.writer_ == nullptr || pending_flux.storage_ == nullptr) {
    impl.current = {};
    impl.current_corrected_c1 = {};
    impl.current_pressure_reference = {};
    impl.current_absolute_pressure_reference = 0.0;
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  FaceFluxView output;
  output.x = pending_flux.x_;
  output.y = pending_flux.y_;
  output.z = pending_flux.z_;
  output.revision = pending_flux.revision_;
  const Status status = correct_exact_coupled_state_impl(
      pressure, pressure_correction, enthalpy_correction, state, candidate,
      output, 2U, reductions, certificate);
  if (status) {
    impl.corrected_pending = {certificate,
                              &pending_flux,
                              pending_flux.writer_,
                              pending_flux.storage_,
                              pending_flux.revision_,
                              pending_flux.writer_identity_,
                              pending_flux.attempt_identity_};
    impl.corrected_pending.pressure_compressibility =
        make_piso_field_revision_identity(
            candidate.pressure_compressibility);
    impl.corrected_pending.active_cells = impl.continuity_activity.cells;
    impl.corrected_pending.activity_local_fingerprint =
        impl.continuity_activity.local_fingerprint;
    impl.corrected_pending.activity_collective_fingerprint =
        impl.continuity_activity.collective_fingerprint;
  }
  return status;
}

Status PressureVelocityCoupler::inspect_corrected_pending(
    const PisoStateCorrectionCertificate& correction,
    const PendingFaceFluxView& pending_flux,
    ConstFaceFluxView& view) const noexcept {
  view = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  const Impl::CorrectedPendingAuthority& authority = impl.corrected_pending;
  const PisoStateCorrectionCertificate& expected = authority.correction;
  const bool same_correction =
      correction.valid() && correction.corrector == 2U &&
      same_state_correction_certificate(correction, expected);
  const bool valid =
      same_correction && authority.view == &pending_flux &&
      pending_flux.valid() && pending_flux.writer_ != nullptr &&
      pending_flux.storage_ != nullptr &&
      pending_flux.writer_ == authority.writer &&
      pending_flux.storage_ == authority.storage &&
      pending_flux.revision_ == authority.revision &&
      pending_flux.revision_ == correction.face_flux &&
      pending_flux.writer_identity_ == authority.writer_identity &&
      pending_flux.attempt_identity_ == authority.attempt_identity;
  if (!valid) return {StatusCode::invalid_plan, kPisoCoupler};

  ConstFaceFluxView candidate{as_const(pending_flux.x_),
                              as_const(pending_flux.y_),
                              as_const(pending_flux.z_),
                              pending_flux.revision_, {}};
  FaceFluxCertificate& flux_certificate = candidate.certificate;
  flux_certificate.revision_ = pending_flux.revision_;
  flux_certificate.authority_ = pending_flux.writer_identity_;
  flux_certificate.storage_ = candidate.x.storage_identity;
  flux_certificate.revision_domain_ = candidate.x.revision_domain;
  flux_certificate.x_base_ = candidate.x.base;
  flux_certificate.y_base_ = candidate.y.base;
  flux_certificate.z_base_ = candidate.z.base;
  flux_certificate.x_stride_y_ = candidate.x.stride_y;
  flux_certificate.x_stride_z_ = candidate.x.stride_z;
  flux_certificate.y_stride_y_ = candidate.y.stride_y;
  flux_certificate.y_stride_z_ = candidate.y.stride_z;
  flux_certificate.z_stride_y_ = candidate.z.stride_y;
  flux_certificate.z_stride_z_ = candidate.z.stride_z;
  flux_certificate.cells_ = impl.cells;
  if (!flux_certificate.matches(candidate))
    return {StatusCode::invalid_plan, kPisoCoupler};
  view = candidate;
  return {};
}

Status PressureVelocityCoupler::audit_pending_final(
    const PisoTerminalAuditInput& input,
    const PendingFaceFluxView& pending_flux,
    ReductionEngine& reductions,
    PisoAttemptReport& report,
    PisoTerminalCertificate& certificate,
    Status prerequisite) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  Impl& impl = *implementation_;
  impl.terminal = {};
  Status reduction_binding =
      reductions.validate_communicator(impl.communicator);
  int reduction_binding_lowest = -1;
  reduction_binding = collective_status(
      impl.communicator, reduction_binding, impl.rank, impl.size,
      reduction_binding_lowest);
  if (!reduction_binding) return reduction_binding;
  const Int3 cells = impl.cells;
  ConstFaceFluxView flux;
  flux.x = as_const(pending_flux.x_);
  flux.y = as_const(pending_flux.y_);
  flux.z = as_const(pending_flux.z_);
  flux.revision = pending_flux.revision_;
  const bool bdf2 = input.bdf.order == 2U;
  const bool empty_previous =
      input.density_previous.base == nullptr &&
      input.density_previous.revision == 0U &&
      input.density_previous.storage_identity == 0U &&
      input.density_previous.revision_domain == 0U;
  const bool closed = input.pressure_reference.kind ==
                      PressureReferenceKind::closed_mass;
  const std::size_t local_cell_count =
      static_cast<std::size_t>(cells.x) * cells.y * cells.z;
  const bool valid_active =
      input.active.size == 0U ||
      (input.active.data != nullptr && input.active.size == local_cell_count);
  const BoundaryThermophysicalGhostContext thermophysical_context{
      input.pressure_reference.time,
      impl.pressure_boundary.certificate().geometry_fingerprint,
      input.pressure_reference.pressure_reference,
      impl.pressure_boundary.certificate().source_revision,
      BoundaryThermophysicalGhostPhase::terminal};
  ThermophysicalBoundaryTokens thermophysical_boundary;
  const bool valid_source_lineage =
      impl.terminal_lineage_source.valid() &&
      impl.terminal_lineage_source.corrector == 2U &&
      input.correction.thermophysical_boundary_semantics ==
          impl.terminal_lineage_source.thermophysical_boundary_semantics &&
      input.correction.thermophysical_boundary_target ==
          impl.terminal_lineage_source.thermophysical_boundary_target &&
      input.correction.thermophysical_boundary_rank_local_binding ==
          impl.terminal_lineage_source
              .thermophysical_boundary_rank_local_binding &&
      input.correction.thermophysical_boundary_collective_lineage ==
          impl.terminal_lineage_source
              .thermophysical_boundary_collective_lineage &&
      input.correction.thermophysical_boundary_rank_local_lineage ==
          impl.terminal_lineage_source
              .thermophysical_boundary_rank_local_lineage;
  const bool valid_thermophysical_boundary =
      thermophysical_boundary_tokens(
          *impl.boundary, thermophysical_context,
          input.thermophysical_boundary, input.density,
          impl.thermodynamics, impl.transport,
          thermophysical_boundary) &&
      (entirely_periodic(*impl.boundary) ||
       same_field_identity(
           input.thermophysical_boundary.binding.pressure_perturbation,
           input.pressure_perturbation)) &&
      thermophysical_boundary.semantics ==
          input.correction.thermophysical_boundary_semantics &&
      (closed ||
       (thermophysical_boundary.collective_lineage ==
            input.correction.thermophysical_boundary_collective_lineage &&
        thermophysical_boundary.rank_local_lineage ==
            input.correction.thermophysical_boundary_rank_local_lineage)) &&
      valid_source_lineage;
  const Span<const ConstFieldView> terminal_species =
      input.thermophysical_boundary.binding.independent_species;
  const bool correction_has_composition =
      input.correction.composition_numeric_provenance != 0U;
  bool valid_terminal_species =
      (terminal_species.size == 0U || terminal_species.data != nullptr);
  if (correction_has_composition)
    valid_terminal_species =
        valid_terminal_species &&
        input.correction.independent_species_count == terminal_species.size &&
        terminal_species.size == impl.independent_species_count &&
        (terminal_species.size == 0U ||
         impl.independent_species_semantic_fields != nullptr);
  for (std::size_t species = 0U;
       species < terminal_species.size && valid_terminal_species; ++species)
    valid_terminal_species =
        detail::valid_cell_view(terminal_species.data[species], cells, 0U,
                                1U, 0U) &&
        (!correction_has_composition ||
         terminal_species.data[species].field ==
             impl.independent_species_semantic_fields[species]);
  const PlanFingerprint terminal_composition_numeric =
      valid_terminal_species
          ? candidate_composition_numeric_local_provenance(terminal_species,
                                                           cells)
          : PlanFingerprint{};
  const PlanFingerprint terminal_composition_binding =
      valid_terminal_species
          ? semantic_composition_rank_local_binding(terminal_species)
          : PlanFingerprint{};
  const bool valid_terminal_composition =
      valid_terminal_species &&
      (!correction_has_composition ||
       terminal_composition_numeric ==
               input.correction.composition_numeric_provenance) &&
      (!correction_has_composition ||
       terminal_composition_binding ==
           input.correction.composition_rank_local_binding);
  const Impl::CorrectedPendingAuthority& corrected_pending =
      impl.corrected_pending;
  const bool valid_corrected_pending =
      same_state_correction_certificate(
          input.correction, corrected_pending.correction) &&
      corrected_pending.view == &pending_flux && pending_flux.valid() &&
      pending_flux.writer_ != nullptr && pending_flux.storage_ != nullptr &&
      pending_flux.writer_ == corrected_pending.writer &&
      pending_flux.storage_ == corrected_pending.storage &&
      pending_flux.revision_ == corrected_pending.revision &&
      pending_flux.writer_identity_ == corrected_pending.writer_identity &&
      pending_flux.attempt_identity_ == corrected_pending.attempt_identity &&
      input.active.data == corrected_pending.active_cells.data &&
      input.active.size == corrected_pending.active_cells.size &&
      corrected_pending.activity_local_fingerprint ==
          impl.continuity_activity.local_fingerprint &&
      corrected_pending.activity_collective_fingerprint ==
          impl.continuity_activity.collective_fingerprint &&
      (input.correction.closure != PisoStateClosure::exact_eos ||
       input.pressure_reference.kind != PressureReferenceKind::closed_mass ||
       same_piso_field_revision_identity(
           corrected_pending.pressure_compressibility,
           input.drho_dp_h_y));
  const bool valid =
      input.correction.valid() && input.correction.corrector == 2U &&
      input.correction.plan == impl.fingerprint &&
      input.correction.face_flux == pending_flux.revision_ &&
      input.pressure_reference.valid() &&
      input.pressure_reference.plan == impl.pressure_reference_plan &&
      input.pressure_reference.kind == impl.pressure_reference_kind &&
      same_pressure_reference_certificate(
          input.pressure_reference,
          input.correction.output_pressure_reference) &&
      ((closed && impl.pressure_gauge ==
                       PressureGauge::compressibility_weighted_zero_mean) ||
       (!closed && impl.pressure_gauge ==
                        PressureGauge::absolute_boundary_dirichlet)) &&
      pending_flux.valid() && pending_flux.writer_ != nullptr &&
      pending_flux.storage_ != nullptr &&
      detail::valid_flux_view(flux, cells, pending_flux.revision_) &&
      detail::valid_bdf_coefficients(input.bdf) && input.bdf.a0 > 0.0 &&
      valid_corrected_pending &&
      detail::valid_cell_view(input.density, cells, 0U, 1U, 0U) &&
      detail::valid_cell_view(input.eos_density, cells, 0U, 1U, 0U) &&
      detail::valid_cell_view(input.density_accepted, cells, 0U, 1U, 0U) &&
      (bdf2
           ? detail::valid_cell_view(input.density_previous, cells, 0U, 1U,
                                     0U)
           : empty_previous) &&
      detail::valid_cell_view(input.pressure_perturbation, cells, 0U, 1U,
                              0U) &&
      (!closed || detail::valid_cell_view(input.drho_dp_h_y, cells, 0U, 1U,
                                          0U)) &&
      input.density.revision == input.correction.density &&
      input.pressure_perturbation.revision == input.correction.pressure &&
      std::isfinite(input.boundary_closure_residual) &&
      input.boundary_closure_residual >= 0.0 &&
      std::isfinite(input.energy_residual) && input.energy_residual >= 0.0 &&
      std::isfinite(input.step_dt) && input.step_dt > 0.0 &&
      same_bdf_coefficients(input.bdf, impl.frozen_candidate_bdf) &&
      detail::bdf_matches_time_step(input.bdf, input.step_dt) &&
      std::isfinite(input.convective_cfl_limit) &&
      input.convective_cfl_limit > 0.0 &&
      valid_active &&
      valid_thermophysical_boundary &&
      valid_terminal_composition &&
      (!closed || (std::isfinite(input.closed_mass_target) &&
                   input.closed_mass_target > 0.0));
  const auto face_offset = [](Int3 shape, Int3 face) noexcept {
    return static_cast<std::size_t>(face.x) +
           static_cast<std::size_t>(shape.x) *
               (static_cast<std::size_t>(face.y) +
                static_cast<std::size_t>(shape.y) *
                    static_cast<std::size_t>(face.z));
  };
  double local_max[7]{};
  double local_sum[5]{};
  std::array<ReductionMaximumLocation, 2U> local_cfl_winners{};
  Int3 local_continuity_cell{};
  bool local_continuity_location_valid = false;
  local_max[3U] = input.boundary_closure_residual;
  local_max[4U] = input.energy_residual;
  local_sum[4U] = static_cast<double>(input.boundary_closure_samples);
  Status local = prerequisite;
  if (local && !valid)
    local = {StatusCode::invalid_plan, kPisoCoupler};
  if (local) {
    for (std::int32_t iz = 0; iz < cells.z; ++iz) {
      for (std::int32_t iy = 0; iy < cells.y; ++iy) {
        for (std::int32_t ix = 0; ix < cells.x; ++ix) {
          const Int3 cell{ix, iy, iz};
          const std::size_t flat =
              static_cast<std::size_t>(ix) +
              static_cast<std::size_t>(cells.x) *
                  (static_cast<std::size_t>(iy) +
                   static_cast<std::size_t>(cells.y) * iz);
          if (input.active.size != 0U) {
            if (input.active.data[flat] > 1U) {
              local = {StatusCode::invalid_plan, kPisoCoupler};
              continue;
            }
            if (input.active.data[flat] == 0U) continue;
          }
          const double rho = input.density.unchecked(cell, 0U);
          const double eos = input.eos_density.unchecked(cell, 0U);
          const double rho_n = input.density_accepted.unchecked(cell, 0U);
          const double rho_nm1 =
              bdf2 ? input.density_previous.unchecked(cell, 0U) : 0.0;
          const double pi = input.pressure_perturbation.unchecked(cell, 0U);
          const double volume = detail::cell_volume(*impl.kernels, cell);
          const double fxm = flux.x.unchecked(cell);
          const double fxp = flux.x.unchecked({ix + 1, iy, iz});
          const double fym = flux.y.unchecked(cell);
          const double fyp = flux.y.unchecked({ix, iy + 1, iz});
          const double fzm = flux.z.unchecked(cell);
          const double fzp = flux.z.unchecked({ix, iy, iz + 1});
          const double compressibility =
              closed ? input.drho_dp_h_y.unchecked(cell, 0U) : 0.0;
          double eos_residual = 0.0;
          double continuity_residual = 0.0;
          double mass_contribution = 0.0;
          double volume_contribution = 0.0;
          double absolute_pi = 0.0;
          double compressibility_moment = 0.0;
          double compressibility_weight = 0.0;
          const bool inactive_xm =
              impl.continuity_activity.x_faces.size != 0U &&
              impl.continuity_activity.x_faces
                      .data[face_offset(flux.x.extents, cell)] == 0U;
          const bool inactive_xp =
              impl.continuity_activity.x_faces.size != 0U &&
              impl.continuity_activity.x_faces
                      .data[face_offset(flux.x.extents,
                                        {ix + 1, iy, iz})] == 0U;
          const bool inactive_ym =
              impl.continuity_activity.y_faces.size != 0U &&
              impl.continuity_activity.y_faces
                      .data[face_offset(flux.y.extents, cell)] == 0U;
          const bool inactive_yp =
              impl.continuity_activity.y_faces.size != 0U &&
              impl.continuity_activity.y_faces
                      .data[face_offset(flux.y.extents,
                                        {ix, iy + 1, iz})] == 0U;
          const bool inactive_zm =
              impl.continuity_activity.z_faces.size != 0U &&
              impl.continuity_activity.z_faces
                      .data[face_offset(flux.z.extents, cell)] == 0U;
          const bool inactive_zp =
              impl.continuity_activity.z_faces.size != 0U &&
              impl.continuity_activity.z_faces
                      .data[face_offset(flux.z.extents,
                                        {ix, iy, iz + 1})] == 0U;
          if (!std::isfinite(rho) || !(rho > 0.0) ||
              !std::isfinite(volume) || !(volume > 0.0) ||
              !std::isfinite(fxm) || !std::isfinite(fxp) ||
              !std::isfinite(fym) || !std::isfinite(fyp) ||
              !std::isfinite(fzm) || !std::isfinite(fzp)) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            continue;
          }
          // The IBM activity authority excludes solid cells and seals every
          // inactive fluid/solid control face to zero mass flux.  A non-zero
          // value there is an authority violation, not a flux that may be
          // silently omitted from the CFL reconstruction.
          if ((inactive_xm && fxm != 0.0) ||
              (inactive_xp && fxp != 0.0) ||
              (inactive_ym && fym != 0.0) ||
              (inactive_yp && fyp != 0.0) ||
              (inactive_zm && fzm != 0.0) ||
              (inactive_zp && fzp != 0.0)) {
            local = {StatusCode::invalid_plan, kPisoCoupler};
            continue;
          }
          const std::array<double, 6U> cell_flux{{
              fxm, fxp, fym, fyp, fzm, fzp}};
          const std::array<std::uint8_t, 6U> cell_active{{
              static_cast<std::uint8_t>(!inactive_xm),
              static_cast<std::uint8_t>(!inactive_xp),
              static_cast<std::uint8_t>(!inactive_ym),
              static_cast<std::uint8_t>(!inactive_yp),
              static_cast<std::uint8_t>(!inactive_zm),
              static_cast<std::uint8_t>(!inactive_zp),
          }};
          detail::CellConvectiveCflResult cell_cfl;
          const detail::CellConvectiveCflStatus cell_cfl_status =
              detail::evaluate_cell_convective_cfl(
                  rho, volume, cell_flux, cell_active, input.step_dt,
                  cell_cfl);
          if (cell_cfl_status != detail::CellConvectiveCflStatus::success) {
            const bool authority_failure =
                cell_cfl_status == detail::CellConvectiveCflStatus::
                                       inactive_nonzero_flux ||
                cell_cfl_status ==
                    detail::CellConvectiveCflStatus::invalid_activity;
            local = {authority_failure ? StatusCode::invalid_plan
                                       : StatusCode::numerical_failure,
                     authority_failure ? kPisoCoupler : kPisoNumerical};
            continue;
          }
          if (hf_coast_common_terminal_cell_v1(
                  rho, eos, rho_n, rho_nm1, volume, input.bdf.a0,
                  input.bdf.a1, input.bdf.a2, fxm, fxp, fym, fyp, fzm, fzp,
                  pi, closed ? 1 : 0, compressibility, &eos_residual,
                  &continuity_residual, &mass_contribution,
                  &volume_contribution, &absolute_pi,
                  &compressibility_moment, &compressibility_weight) != 0) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            continue;
          }
          // Record only the local location during the normal audit pass.
          // The detailed continuity witness is reconstructed from this one
          // cell only after the global continuity gate has failed.
          if (!local_continuity_location_valid ||
              continuity_residual > local_max[1U]) {
            local_continuity_cell = cell;
            local_continuity_location_valid = true;
          }
          local_max[0U] = std::max(local_max[0U], eos_residual);
          local_max[1U] = std::max(local_max[1U], continuity_residual);
          local_max[2U] = std::max(local_max[2U], absolute_pi);
          local_max[5U] = std::max(local_max[5U], cell_cfl.out);
          local_max[6U] =
              std::max(local_max[6U], cell_cfl.absolute);
          const Int3 global_cell{impl.patch.begin.x + cell.x,
                                 impl.patch.begin.y + cell.y,
                                 impl.patch.begin.z + cell.z};
          const Int3 global_cells = impl.geometry->global_cells();
          const std::uint64_t global_linear =
              static_cast<std::uint64_t>(global_cell.x) +
              static_cast<std::uint64_t>(global_cells.x) *
                  (static_cast<std::uint64_t>(global_cell.y) +
                   static_cast<std::uint64_t>(global_cells.y) *
                       static_cast<std::uint64_t>(global_cell.z));
          const std::array<double, 5U> payload{{
              cell_cfl.out, cell_cfl.absolute, cell_cfl.density_volume,
              cell_cfl.outgoing_mass_flow, cell_cfl.absolute_mass_flow}};
          const auto select_winner = [&](std::size_t index,
                                         double value) noexcept {
            ReductionMaximumLocation& winner = local_cfl_winners[index];
            if (!winner.valid || value > winner.value ||
                (value == winner.value &&
                 global_linear < winner.global_location)) {
              winner.valid = true;
              winner.value = value;
              winner.global_location = global_linear;
              winner.rank = impl.rank;
              winner.payload = payload;
            }
          };
          select_winner(0U, cell_cfl.out);
          select_winner(1U, cell_cfl.absolute);
          local_sum[0U] += mass_contribution;
          local_sum[1U] += volume_contribution;
          local_sum[2U] += compressibility_moment;
          local_sum[3U] += compressibility_weight;
        }
      }
    }
  }
  double global_max[7]{};
  Status status = reductions.checked_max({local_max, 7U},
                                         {global_max, 7U}, local);
  if (!status) {
    return status;
  }
  std::array<ReductionMaximumLocation, 2U> global_cfl_winners{};
  status = reductions.checked_max_locations(
      {local_cfl_winners.data(), local_cfl_winners.size()},
      {global_cfl_winners.data(), global_cfl_winners.size()});
  if (!status) return status;
  if (global_cfl_winners[0U].value != global_max[5U] ||
      global_cfl_winners[1U].value != global_max[6U]) {
    return {StatusCode::numerical_failure, kPisoNumerical};
  }
  double global_sum[5]{};
  status = reductions.checked_sum({local_sum, 5U}, {global_sum, 5U});
  if (!status) {
    return status;
  }
  if (!closed && (!(global_sum[4U] >= 1.0) ||
                  !std::isfinite(global_sum[4U]))) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  double mass_residual = 0.0;
  double gauge_residual = 0.0;
  if (hf_coast_common_terminal_finalize_v1(
          closed ? 1 : 0, global_sum[0U], global_sum[1U],
          input.closed_mass_target, global_sum[2U], global_sum[3U],
          global_max[2U], global_max[3U], &mass_residual,
          &gauge_residual) != 0) {
    return {StatusCode::numerical_failure, kPisoNumerical};
  }
  PisoContinuityWitness global_continuity_witness;
  if (global_max[1U] > impl.continuity_tolerance) {
    PisoContinuityWitness local_continuity_witness;
    if (local_continuity_location_valid) {
      const Int3 cell = local_continuity_cell;
      const double rho = input.density.unchecked(cell, 0U);
      const double rho_n = input.density_accepted.unchecked(cell, 0U);
      const double rho_nm1 =
          bdf2 ? input.density_previous.unchecked(cell, 0U) : 0.0;
      const double volume = detail::cell_volume(*impl.kernels, cell);
      const double fxm = flux.x.unchecked(cell);
      const double fxp = flux.x.unchecked({cell.x + 1, cell.y, cell.z});
      const double fym = flux.y.unchecked(cell);
      const double fyp = flux.y.unchecked({cell.x, cell.y + 1, cell.z});
      const double fzm = flux.z.unchecked(cell);
      const double fzp = flux.z.unchecked({cell.x, cell.y, cell.z + 1});
      const double unsteady =
          volume * (input.bdf.a0 * rho + input.bdf.a1 * rho_n +
                    input.bdf.a2 * rho_nm1);
      const double divergence =
          (fxp - fxm) + (fyp - fym) + (fzp - fzm);
      const double scale =
          std::abs(volume * input.bdf.a0 * rho) +
          std::abs(volume * input.bdf.a1 * rho_n) +
          std::abs(volume * input.bdf.a2 * rho_nm1) + std::abs(fxm) +
          std::abs(fxp) + std::abs(fym) + std::abs(fyp) + std::abs(fzm) +
          std::abs(fzp);
      const Int3 global{impl.patch.begin.x + cell.x,
                        impl.patch.begin.y + cell.y,
                        impl.patch.begin.z + cell.z};
      const Int3 global_cells = impl.geometry->global_cells();
      const std::uint64_t global_cell =
          static_cast<std::uint64_t>(global.x) +
          static_cast<std::uint64_t>(global_cells.x) *
              (static_cast<std::uint64_t>(global.y) +
               static_cast<std::uint64_t>(global_cells.y) *
                   static_cast<std::uint64_t>(global.z));
      local_continuity_witness = {
          true,
          global_cell,
          global,
          impl.rank,
          local_max[1U],
          unsteady + divergence,
          unsteady,
          divergence,
          scale,
          rho,
          rho_n,
          rho_nm1,
          {fxm, fxp, fym, fyp, fzm, fzp},
          global.x == 0 || global.y == 0 || global.z == 0 ||
              global.x + 1 == global_cells.x ||
              global.y + 1 == global_cells.y ||
              global.z + 1 == global_cells.z};
    }
    struct MaxLocation {
      double value;
      int rank;
    } local_location{local_continuity_witness.valid
                         ? local_continuity_witness.normalized_residual
                         : -1.0,
                     impl.rank},
        global_location{};
    if (MPI_Allreduce(&local_location, &global_location, 1, MPI_DOUBLE_INT,
                      MPI_MAXLOC, impl.communicator) != MPI_SUCCESS ||
        global_location.rank < 0 || global_location.rank >= impl.size) {
      return {StatusCode::mpi_failure, kPisoCollective};
    }
    if (impl.rank == global_location.rank)
      global_continuity_witness = local_continuity_witness;
    if (MPI_Bcast(&global_continuity_witness,
                  static_cast<int>(sizeof(global_continuity_witness)),
                  MPI_BYTE, global_location.rank,
                  impl.communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kPisoCollective};
    }
  }
  const Int3 global_cells = impl.geometry->global_cells();
  const auto cfl_witness = [&](const ReductionMaximumLocation& winner) {
    ConvectiveCflFailureWitness witness;
    if (!winner.valid || global_cells.x <= 0 || global_cells.y <= 0 ||
        global_cells.z <= 0)
      return witness;
    std::uint64_t remainder = winner.global_location;
    const std::uint64_t plane =
        static_cast<std::uint64_t>(global_cells.x) *
        static_cast<std::uint64_t>(global_cells.y);
    const std::uint64_t z = remainder / plane;
    remainder -= z * plane;
    const std::uint64_t y =
        remainder / static_cast<std::uint64_t>(global_cells.x);
    const std::uint64_t x =
        remainder - y * static_cast<std::uint64_t>(global_cells.x);
    if (x >= static_cast<std::uint64_t>(global_cells.x) ||
        y >= static_cast<std::uint64_t>(global_cells.y) ||
        z >= static_cast<std::uint64_t>(global_cells.z))
      return witness;
    witness.valid = true;
    witness.global_cell = {static_cast<std::int32_t>(x),
                           static_cast<std::int32_t>(y),
                           static_cast<std::int32_t>(z)};
    witness.rank = winner.rank;
    witness.out = winner.payload[0U];
    witness.absolute = winner.payload[1U];
    witness.density_volume = winner.payload[2U];
    witness.outgoing_mass_flow = winner.payload[3U];
    witness.absolute_mass_flow = winner.payload[4U];
    return witness;
  };
  CommittedConvectiveCflCertificate committed_cfl;
  committed_cfl.plan = impl.fingerprint;
  committed_cfl.correction_state = input.correction.state;
  committed_cfl.density = input.density.revision;
  committed_cfl.final_flux = pending_flux.revision_;
  committed_cfl.density_storage = input.density.storage_identity;
  committed_cfl.density_revision_domain = input.density.revision_domain;
  committed_cfl.face_flux_storage = flux.x.storage_identity;
  committed_cfl.face_flux_revision_domain = flux.x.revision_domain;
  committed_cfl.density_view_identity = {input.density};
  committed_cfl.face_flux_view_identity = {flux};
  committed_cfl.activity_collective =
      impl.continuity_activity.collective_fingerprint;
  committed_cfl.dt = input.step_dt;
  committed_cfl.out_max = global_max[5U];
  committed_cfl.absolute_max = global_max[6U];
  committed_cfl.limit = input.convective_cfl_limit;
  committed_cfl.out_winner = cfl_witness(global_cfl_winners[0U]);
  committed_cfl.absolute_winner = cfl_witness(global_cfl_winners[1U]);
  constexpr double kCommittedCflRoundoffSlack =
      64.0 * std::numeric_limits<double>::epsilon();
  if (committed_cfl.out_max >
      committed_cfl.limit * (1.0 + kCommittedCflRoundoffSlack)) {
    committed_cfl.failure_witness.valid = true;
    committed_cfl.failure_witness.out_winner = committed_cfl.out_winner;
    committed_cfl.failure_witness.absolute_winner =
        committed_cfl.absolute_winner;
  }
  if (!committed_cfl.valid())
    return {StatusCode::numerical_failure, kPisoNumerical};
  PisoAttemptReport candidate_report = report;
  candidate_report.eos_residual = global_max[0U];
  candidate_report.continuity_residual = global_max[1U];
  candidate_report.energy_residual = global_max[4U];
  candidate_report.closed_mass_residual = mass_residual;
  candidate_report.gauge_residual = gauge_residual;
  candidate_report.committed_convective_cfl_out_max = global_max[5U];
  candidate_report.committed_convective_cfl_abs_max = global_max[6U];
  candidate_report.committed_convective_cfl_limit =
      input.convective_cfl_limit;
  candidate_report.committed_convective_cfl = committed_cfl;
  candidate_report.final_flux_revision = pending_flux.revision_;
  candidate_report.continuity_witness = global_continuity_witness;
  report = candidate_report;
  int accepted = hf_coast_common_terminal_accept_v1(
      global_max[0U], global_max[1U], mass_residual, gauge_residual,
      impl.eos_tolerance, impl.continuity_tolerance,
      impl.closed_mass_tolerance, impl.gauge_tolerance);
  if (accepted > 0 && impl.energy_tolerance > 0.0 &&
      global_max[4U] > impl.energy_tolerance)
    accepted = 0;
  if (accepted < 0) {
    return {StatusCode::numerical_failure, kPisoNumerical};
  }
  if (accepted == 0) {
    return {StatusCode::rejected_step, kPisoSolve};
  }
  std::uint64_t audit = kFnvOffset;
  audit = hash_mix(audit, input.correction.state);
  audit = hash_mix(audit, pending_flux.revision_);
  audit = hash_mix(audit, double_bits(global_max[0U]));
  audit = hash_mix(audit, double_bits(global_max[1U]));
  audit = hash_mix(audit, double_bits(global_max[4U]));
  audit = hash_mix(audit, double_bits(mass_residual));
  audit = hash_mix(audit, double_bits(gauge_residual));
  audit = hash_mix(audit, double_bits(input.step_dt));
  audit = hash_mix(audit, double_bits(global_max[5U]));
  audit = hash_mix(audit, double_bits(global_max[6U]));
  audit = hash_mix(audit, double_bits(input.convective_cfl_limit));
  audit = hash_mix(audit, committed_cfl.density);
  audit = hash_mix(audit, committed_cfl.final_flux);
  audit = hash_mix(audit, committed_cfl.density_storage);
  audit = hash_mix(audit, committed_cfl.density_revision_domain);
  audit = hash_mix(audit, committed_cfl.face_flux_storage);
  audit = hash_mix(audit, committed_cfl.face_flux_revision_domain);
  audit = hash_mix(audit, committed_cfl.activity_collective);
  audit = hash_mix(audit, global_cfl_winners[0U].global_location);
  audit = hash_mix(audit, global_cfl_winners[1U].global_location);
  audit = hash_mix(
      audit, static_cast<std::uint64_t>(committed_cfl.out_winner.rank));
  audit = hash_mix(
      audit, static_cast<std::uint64_t>(committed_cfl.absolute_winner.rank));
  audit = hash_mix(audit,
                   double_bits(committed_cfl.out_winner.density_volume));
  audit = hash_mix(
      audit, double_bits(committed_cfl.out_winner.outgoing_mass_flow));
  audit = hash_mix(
      audit, double_bits(committed_cfl.out_winner.absolute_mass_flow));
  audit = hash_mix(
      audit, double_bits(committed_cfl.absolute_winner.density_volume));
  audit = hash_mix(
      audit,
      double_bits(committed_cfl.absolute_winner.outgoing_mass_flow));
  audit = hash_mix(
      audit,
      double_bits(committed_cfl.absolute_winner.absolute_mass_flow));
  audit = hash_mix(audit, thermophysical_boundary.semantics);
  audit = hash_mix(audit, thermophysical_boundary.target);
  audit = hash_mix(audit, thermophysical_boundary.rank_local_binding);
  audit = hash_mix(audit, thermophysical_boundary.collective_lineage);
  audit = hash_mix(audit, thermophysical_boundary.rank_local_lineage);
  audit = hash_mix(audit, input.pressure_reference.pressure_reference);
  audit = hash_mix(
      audit, input.correction.closed_gauge_collective_transaction);
  audit = hash_mix(audit, terminal_composition_numeric);
  audit = hash_mix(audit, terminal_composition_binding);
  audit = hash_mix(
      audit,
      static_cast<std::uint64_t>(terminal_species.size));
  PisoTerminalCertificate candidate;
  candidate.plan = impl.fingerprint;
  candidate.correction_state = input.correction.state;
  candidate.final_flux = pending_flux.revision_;
  candidate.audit_state = audit == 0U ? 1U : audit;
  candidate.thermophysical_boundary_semantics =
      thermophysical_boundary.semantics;
  candidate.thermophysical_boundary_target = thermophysical_boundary.target;
  candidate.thermophysical_boundary_rank_local_binding =
      thermophysical_boundary.rank_local_binding;
  candidate.thermophysical_boundary_collective_lineage =
      thermophysical_boundary.collective_lineage;
  candidate.thermophysical_boundary_rank_local_lineage =
      thermophysical_boundary.rank_local_lineage;
  candidate.pressure_reference = input.pressure_reference;
  candidate.closed_gauge_collective_transaction =
      input.correction.closed_gauge_collective_transaction;
  candidate.composition_numeric_provenance =
      correction_has_composition
          ? input.correction.composition_numeric_provenance
          : terminal_composition_numeric;
  candidate.composition_rank_local_binding =
      correction_has_composition
          ? input.correction.composition_rank_local_binding
          : terminal_composition_binding;
  candidate.independent_species_count =
      correction_has_composition
          ? input.correction.independent_species_count
          : terminal_species.size;
  impl.terminal = candidate;
  certificate = candidate;
  return {};
}

Status PressureVelocityCoupler::publish_pending_final(
    const PisoTerminalCertificate& terminal,
    Span<const RevisionDependency> dependencies,
    Span<const ConstFieldView> independent_species,
    ReductionEngine& reductions,
    FinalFaceFluxWriter& writer,
    PendingFaceFluxView& pending_flux) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  Impl& impl = *implementation_;
  Status local = reductions.validate_communicator(impl.communicator);
  bool composition_views_valid =
      terminal.independent_species_count == independent_species.size &&
      independent_species.size == impl.independent_species_count &&
      (independent_species.size == 0U ||
       (independent_species.data != nullptr &&
        impl.independent_species_semantic_fields != nullptr));
  for (std::size_t species = 0U;
       species < independent_species.size && composition_views_valid;
       ++species)
    composition_views_valid =
        detail::valid_cell_view(independent_species.data[species], impl.cells,
                                0U, 1U, 0U) &&
        independent_species.data[species].field ==
            impl.independent_species_semantic_fields[species];
  const PlanFingerprint replayed_composition =
      composition_views_valid
          ? candidate_composition_numeric_local_provenance(
                independent_species, impl.cells)
          : PlanFingerprint{};
  const PlanFingerprint replayed_composition_binding =
      composition_views_valid
          ? semantic_composition_rank_local_binding(independent_species)
          : PlanFingerprint{};
  const bool composition_matches =
      terminal.composition_numeric_provenance == 0U
          ? independent_species.size == 0U &&
                terminal.independent_species_count == 0U
          : composition_views_valid &&
                replayed_composition ==
                    terminal.composition_numeric_provenance &&
                replayed_composition_binding ==
                    terminal.composition_rank_local_binding;
  const bool dependencies_valid =
      dependencies.size == 0U || dependencies.data != nullptr;
  bool composition_dependencies_complete = dependencies_valid;
  for (std::size_t species = 0U;
       species < independent_species.size &&
       composition_dependencies_complete;
       ++species) {
    const ConstFieldView field = independent_species.data[species];
    const RevisionSourceId expected_source =
        AttemptTransaction::field_revision_source(field.field);
    std::size_t matches = 0U;
    for (std::size_t dependency = 0U; dependency < dependencies.size;
         ++dependency)
      if (dependencies.data[dependency].source == expected_source &&
          dependencies.data[dependency].revision == field.revision)
        ++matches;
    composition_dependencies_complete = matches == 1U;
  }
  const bool valid = terminal.valid() && impl.terminal.valid() &&
                     terminal.plan == impl.terminal.plan &&
                     terminal.correction_state ==
                         impl.terminal.correction_state &&
                     terminal.final_flux == impl.terminal.final_flux &&
                     terminal.audit_state == impl.terminal.audit_state &&
                     terminal.thermophysical_boundary_semantics ==
                         impl.terminal.thermophysical_boundary_semantics &&
                     terminal.thermophysical_boundary_target ==
                         impl.terminal.thermophysical_boundary_target &&
                     terminal.thermophysical_boundary_rank_local_binding ==
                         impl.terminal
                             .thermophysical_boundary_rank_local_binding &&
                     terminal.thermophysical_boundary_collective_lineage ==
                         impl.terminal
                             .thermophysical_boundary_collective_lineage &&
                     terminal.thermophysical_boundary_rank_local_lineage ==
                         impl.terminal
                             .thermophysical_boundary_rank_local_lineage &&
                     same_pressure_reference_certificate(
                         terminal.pressure_reference,
                         impl.terminal.pressure_reference) &&
                     terminal.closed_gauge_collective_transaction ==
                         impl.terminal.closed_gauge_collective_transaction &&
                     terminal.composition_numeric_provenance ==
                         impl.terminal.composition_numeric_provenance &&
                     terminal.composition_rank_local_binding ==
                         impl.terminal.composition_rank_local_binding &&
                     terminal.independent_species_count ==
                         impl.terminal.independent_species_count &&
                     composition_matches &&
                     composition_dependencies_complete &&
                     pending_flux.valid() &&
                     pending_flux.revision_ == terminal.final_flux &&
                     pending_flux.writer_ == &writer;
  if (local && !valid)
    local = {StatusCode::invalid_plan, kPisoCoupler};
  local = reductions.consensus(local);
  if (!local) {
    impl.terminal = {};
    return local;
  }
  Status preflight =
      writer.preflight_publish_pending(dependencies, pending_flux);
  preflight = reductions.consensus(preflight);
  if (!preflight) {
    impl.terminal = {};
    return preflight;
  }
  const Status published = writer.publish_pending(dependencies, pending_flux);
  const Status status = reductions.consensus(published);
  impl.terminal = {};
  return status;
}

Status PressureVelocityCoupler::inspect_pending_final(
    const PisoTerminalCertificate& terminal,
    const PendingFaceFluxView& pending_flux,
    ConstFaceFluxView& view) const noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoCoupler};
  }
  const Impl& impl = *implementation_;
  const bool valid = terminal.valid() && impl.terminal.valid() &&
                     terminal.plan == impl.terminal.plan &&
                     terminal.correction_state ==
                         impl.terminal.correction_state &&
                     terminal.final_flux == impl.terminal.final_flux &&
                     terminal.audit_state == impl.terminal.audit_state &&
                     terminal.thermophysical_boundary_semantics ==
                         impl.terminal.thermophysical_boundary_semantics &&
                     terminal.thermophysical_boundary_target ==
                         impl.terminal.thermophysical_boundary_target &&
                     terminal.thermophysical_boundary_rank_local_binding ==
                         impl.terminal
                             .thermophysical_boundary_rank_local_binding &&
                     terminal.thermophysical_boundary_collective_lineage ==
                         impl.terminal
                             .thermophysical_boundary_collective_lineage &&
                     terminal.thermophysical_boundary_rank_local_lineage ==
                         impl.terminal
                             .thermophysical_boundary_rank_local_lineage &&
                     same_pressure_reference_certificate(
                         terminal.pressure_reference,
                         impl.terminal.pressure_reference) &&
                     terminal.closed_gauge_collective_transaction ==
                         impl.terminal.closed_gauge_collective_transaction &&
                     terminal.composition_numeric_provenance ==
                         impl.terminal.composition_numeric_provenance &&
                     terminal.composition_rank_local_binding ==
                         impl.terminal.composition_rank_local_binding &&
                     terminal.independent_species_count ==
                         impl.terminal.independent_species_count &&
                     pending_flux.valid() &&
                     pending_flux.revision_ == terminal.final_flux &&
                     pending_flux.writer_ != nullptr &&
                     pending_flux.storage_ != nullptr;
  if (!valid) return {StatusCode::invalid_plan, kPisoCoupler};
  ConstFaceFluxView candidate{as_const(pending_flux.x_),
                              as_const(pending_flux.y_),
                              as_const(pending_flux.z_),
                              pending_flux.revision_, {}};
  FaceFluxCertificate& certificate = candidate.certificate;
  certificate.revision_ = pending_flux.revision_;
  certificate.authority_ = pending_flux.writer_identity_;
  certificate.storage_ = candidate.x.storage_identity;
  certificate.revision_domain_ = candidate.x.revision_domain;
  certificate.x_base_ = candidate.x.base;
  certificate.y_base_ = candidate.y.base;
  certificate.z_base_ = candidate.z.base;
  certificate.x_stride_y_ = candidate.x.stride_y;
  certificate.x_stride_z_ = candidate.x.stride_z;
  certificate.y_stride_y_ = candidate.y.stride_y;
  certificate.y_stride_z_ = candidate.y.stride_z;
  certificate.z_stride_y_ = candidate.z.stride_y;
  certificate.z_stride_z_ = candidate.z.stride_z;
  certificate.cells_ = impl.cells;
  if (!certificate.matches(candidate))
    return {StatusCode::invalid_plan, kPisoCoupler};
  view = candidate;
  return {};
}

class PressureContinuityConvergenceAudit final
    : public LinearConvergenceAudit {
 public:
  struct Provisional {
    bool valid{};
    PressureCorrectionCertificate pressure{};
    double alpha{};
    ConstFieldView audited_workspace{};
  };

  PressureContinuityConvergenceAudit(
      PressureVelocityCoupler& coupler,
      PressureCorrectionCertificate pressure,
      ConstFieldView sealed_rhs) noexcept
      : coupler_(&coupler), pressure_(pressure), sealed_rhs_(sealed_rhs) {}

  LinearConvergenceAuditCertificate certificate() const noexcept override {
    return coupler_ == nullptr
               ? LinearConvergenceAuditCertificate{}
               : coupler_->pressure_convergence_certificate(pressure_);
  }

  Status evaluate(ConstFieldView solution, ConstFieldView true_residual,
                  ReductionEngine& reductions,
                  LinearConvergenceAuditResult& result) noexcept override {
    provisional_ = {};
    if (coupler_ == nullptr)
      return {StatusCode::invalid_plan, kPisoSolve};
    double selected_alpha = 1.0;
    Status status = coupler_->audit_pressure_convergence(
        pressure_, solution, true_residual, reductions, selected_alpha,
        result);
    if (status && result.terminal_rejection) {
      status = coupler_->capture_pressure_failure_provenance(
          pressure_, solution, sealed_rhs_, reductions, result);
    }
    if (status && result.accepted) {
      provisional_.valid = true;
      provisional_.pressure = pressure_;
      provisional_.alpha = selected_alpha;
      provisional_.audited_workspace = solution;
    }
    return status;
  }

  bool take_provisional(Provisional& out) noexcept {
    out = provisional_;
    provisional_ = {};
    return out.valid;
  }

 private:
  PressureVelocityCoupler* coupler_{};
  PressureCorrectionCertificate pressure_{};
  ConstFieldView sealed_rhs_{};
  Provisional provisional_{};
};

LinearConvergenceAuditCertificate
PressureVelocityCoupler::pressure_convergence_certificate(
    const PressureCorrectionCertificate& pressure) const noexcept {
  if (implementation_ == nullptr ||
      !implementation_->pressure_boundary.current() || !pressure.valid() ||
      pressure.corrector != 2U ||
      !same_pressure_certificate(
          pressure, implementation_->pressure_correction)) {
    return {};
  }
  const Impl& impl = *implementation_;
  std::uint64_t fingerprint = kFnvOffset;
  fingerprint = hash_mix(fingerprint, UINT64_C(0x7032636f6e74696e));
  fingerprint = hash_mix(fingerprint, impl.plan);
  fingerprint = hash_mix(fingerprint, pressure.time);
  fingerprint = hash_mix(fingerprint, pressure.geometry);
  fingerprint = hash_mix(fingerprint, pressure.numeric_boundary);
  fingerprint = hash_mix(
      fingerprint, pressure.thermophysical_boundary_semantics);
  fingerprint = hash_mix(
      fingerprint, pressure.thermophysical_boundary_collective_lineage);
  fingerprint =
      hash_mix(fingerprint, pressure.thermophysical_boundary_target);
  fingerprint = hash_mix(fingerprint,
                         kPressureCorrectionAuditSemanticVersion);
  fingerprint = hash_mix(fingerprint, double_bits(impl.continuity_tolerance));
  fingerprint = hash_mix(
      fingerprint, double_bits(kPressureCorrectionSafetyFraction));
  fingerprint = hash_mix(
      fingerprint,
      impl.continuity_activity.collective_fingerprint == 0U
          ? UINT64_C(0x66756c6c646f6d61)
          : impl.continuity_activity.collective_fingerprint);
  return {fingerprint == 0U ? 1U : fingerprint};
}

Status PressureVelocityCoupler::audit_pressure_convergence(
    const PressureCorrectionCertificate& pressure,
    ConstFieldView correction, ConstFieldView true_residual,
    ReductionEngine& reductions, double& selected_alpha,
    LinearConvergenceAuditResult& result) noexcept {
  if (implementation_ == nullptr ||
      pressure_convergence_certificate(pressure).collective_fingerprint ==
          0U) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }
  Impl& impl = *implementation_;
  selected_alpha = 1.0;
  const PressureCorrectionInput& input = impl.pressure_input;
  const Int3 cells = impl.cells;
  const bool bdf2 = input.bdf.order == 2U;
  if (!detail::valid_cell_view(correction, cells, 0U, 1U, 1U) ||
      !detail::valid_cell_view(true_residual, cells, 0U, 1U, 0U) ||
      input.intermediate.corrector != 2U ||
      input.intermediate.dependency != pressure.intermediate ||
      input.intermediate.thermophysical_boundary_semantics !=
          pressure.thermophysical_boundary_semantics ||
      input.intermediate.thermophysical_boundary_target !=
          pressure.thermophysical_boundary_target ||
      input.intermediate.thermophysical_boundary_rank_local_binding !=
          pressure.thermophysical_boundary_rank_local_binding ||
      input.intermediate.thermophysical_boundary_collective_lineage !=
          pressure.thermophysical_boundary_collective_lineage ||
      input.intermediate.thermophysical_boundary_rank_local_lineage !=
          pressure.thermophysical_boundary_rank_local_lineage ||
      !detail::valid_cell_view(input.density_trial, cells, 0U, 1U, 0U) ||
      !detail::valid_cell_view(input.density_accepted, cells, 0U, 1U, 0U) ||
      (bdf2 && !detail::valid_cell_view(input.density_previous, cells, 0U,
                                       1U, 0U)) ||
      !detail::valid_cell_view(input.drho_dp_h_y, cells, 0U, 1U, 0U)) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }

  const PressureContinuityActivityView activity = impl.continuity_activity;
  const auto offset = [](Int3 shape, Int3 value) noexcept {
    return static_cast<std::size_t>(value.x) +
           static_cast<std::size_t>(shape.x) *
               (static_cast<std::size_t>(value.y) +
                static_cast<std::size_t>(shape.y) *
                    static_cast<std::size_t>(value.z));
  };
  double local_max_depletion = 0.0;
  Status local;
  for (std::int32_t iz = 0; iz < cells.z && local; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y && local; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        const double rho = input.density_trial.unchecked(cell, 0U);
        const double psi = input.drho_dp_h_y.unchecked(cell, 0U);
        const double delta = correction.unchecked(cell, 0U);
        if (!finite_positive(rho) || !std::isfinite(psi) || psi < 0.0 ||
            !std::isfinite(delta)) {
          local = {StatusCode::numerical_failure, kPisoNumerical};
          break;
        }
        const double psi_delta = psi * delta;
        if (!std::isfinite(psi_delta)) {
          local = {StatusCode::numerical_failure, kPisoNumerical};
          break;
        }
        const double depletion = std::max(0.0, -psi_delta / rho);
        if (!std::isfinite(depletion)) {
          local = {StatusCode::numerical_failure, kPisoNumerical};
          break;
        }
        local_max_depletion = std::max(local_max_depletion, depletion);
      }
    }
  }
  double global_max_depletion = 0.0;
  Status status = reductions.checked_max(
      {&local_max_depletion, 1U}, {&global_max_depletion, 1U}, local);
  if (!status) return status;
  double alpha = 1.0;
  if (global_max_depletion > kPressureCorrectionSafetyFraction) {
    alpha = std::nextafter(
        kPressureCorrectionSafetyFraction / global_max_depletion, 0.0);
  }
  if (!std::isfinite(alpha) || !(alpha > 0.0) || alpha > 1.0)
    return {StatusCode::numerical_failure, kPisoNumerical};
  const auto corrected_face = [&](CartesianAxis axis, Int3 face,
                                  double application_scale) noexcept {
    const ConstFaceFieldView source =
        detail::select(as_const(impl.workspace.phi_h_by_a), axis);
    const ConstFaceFieldView coefficient =
        axis == CartesianAxis::x
            ? as_const(impl.workspace.x_pressure_coefficient)
            : (axis == CartesianAxis::y
                   ? as_const(impl.workspace.y_pressure_coefficient)
                   : as_const(impl.workspace.z_pressure_coefficient));
    Span<const std::uint8_t> active_faces;
    if (axis == CartesianAxis::x)
      active_faces = activity.x_faces;
    else if (axis == CartesianAxis::y)
      active_faces = activity.y_faces;
    else
      active_faces = activity.z_faces;
    if (active_faces.size != 0U &&
        active_faces.data[offset(source.extents, face)] == 0U)
      return 0.0;
    return source.unchecked(face) +
           impl.pressure_boundary.mass_flux_response_unchecked(
               correction, axis, face, coefficient.unchecked(face),
               application_scale);
  };

  double local_max = 0.0;
  for (std::int32_t iz = 0; iz < cells.z; ++iz) {
    for (std::int32_t iy = 0; iy < cells.y; ++iy) {
      for (std::int32_t ix = 0; ix < cells.x; ++ix) {
        const Int3 cell{ix, iy, iz};
        if (activity.cells.size != 0U &&
            activity.cells.data[offset(cells, cell)] == 0U)
          continue;
        const double raw_delta = correction.unchecked(cell, 0U);
        const double delta = alpha == 1.0 ? raw_delta : alpha * raw_delta;
        const double rho = input.density_trial.unchecked(cell, 0U) +
                           input.drho_dp_h_y.unchecked(cell, 0U) * delta;
        if (!std::isfinite(raw_delta) || !std::isfinite(delta) ||
            !std::isfinite(rho) || !(rho > 0.0)) {
          local = {StatusCode::numerical_failure, kPisoNumerical};
          continue;
        }
        const double rho_n = input.density_accepted.unchecked(cell, 0U);
        const double rho_nm1 =
            bdf2 ? input.density_previous.unchecked(cell, 0U) : 0.0;
        const double fxm = corrected_face(CartesianAxis::x, cell, alpha);
        const double fxp =
            corrected_face(CartesianAxis::x, {ix + 1, iy, iz}, alpha);
        const double fym = corrected_face(CartesianAxis::y, cell, alpha);
        const double fyp =
            corrected_face(CartesianAxis::y, {ix, iy + 1, iz}, alpha);
        const double fzm = corrected_face(CartesianAxis::z, cell, alpha);
        const double fzp =
            corrected_face(CartesianAxis::z, {ix, iy, iz + 1}, alpha);
        double eos_residual = 0.0;
        double continuity_residual = 0.0;
        double mass = 0.0;
        double volume_sum = 0.0;
        double absolute_pressure = 0.0;
        double pressure_moment = 0.0;
        double pressure_weight = 0.0;
        if (hf_coast_common_terminal_cell_v1(
                rho, rho, rho_n, rho_nm1,
                detail::cell_volume(*impl.kernels, cell), input.bdf.a0,
                input.bdf.a1, input.bdf.a2, fxm, fxp, fym, fyp, fzm, fzp,
                0.0, 0, 0.0, &eos_residual, &continuity_residual, &mass,
                &volume_sum, &absolute_pressure, &pressure_moment,
                &pressure_weight) != 0) {
          local = {StatusCode::numerical_failure, kPisoNumerical};
          continue;
        }
        local_max = std::max(local_max, continuity_residual);
      }
    }
  }
  double local_unscaled_max = local_max;
  double local_operator_parity_error = 0.0;
  if (alpha < 1.0 && local) {
    local_unscaled_max = 0.0;
    for (std::int32_t iz = 0; iz < cells.z; ++iz) {
      for (std::int32_t iy = 0; iy < cells.y; ++iy) {
        for (std::int32_t ix = 0; ix < cells.x; ++ix) {
          const Int3 cell{ix, iy, iz};
          if (activity.cells.size != 0U &&
              activity.cells.data[offset(cells, cell)] == 0U)
            continue;
          const double delta = correction.unchecked(cell, 0U);
          const double rho = input.density_trial.unchecked(cell, 0U) +
                             input.drho_dp_h_y.unchecked(cell, 0U) * delta;
          const double rho_n = input.density_accepted.unchecked(cell, 0U);
          const double rho_nm1 =
              bdf2 ? input.density_previous.unchecked(cell, 0U) : 0.0;
          const double fxm = corrected_face(CartesianAxis::x, cell, 1.0);
          const double fxp =
              corrected_face(CartesianAxis::x, {ix + 1, iy, iz}, 1.0);
          const double fym = corrected_face(CartesianAxis::y, cell, 1.0);
          const double fyp =
              corrected_face(CartesianAxis::y, {ix, iy + 1, iz}, 1.0);
          const double fzm = corrected_face(CartesianAxis::z, cell, 1.0);
          const double fzp =
              corrected_face(CartesianAxis::z, {ix, iy, iz + 1}, 1.0);
          const double volume = detail::cell_volume(*impl.kernels, cell);
          const double unsteady =
              volume * (input.bdf.a0 * rho + input.bdf.a1 * rho_n +
                        input.bdf.a2 * rho_nm1);
          const double flux_sum = (fxp - fxm) + (fyp - fym) + (fzp - fzm);
          const double raw_continuity = unsteady + flux_sum;
          const double scale =
              std::abs(volume * input.bdf.a0 * rho) +
              std::abs(volume * input.bdf.a1 * rho_n) +
              std::abs(volume * input.bdf.a2 * rho_nm1) + std::abs(fxm) +
              std::abs(fxp) + std::abs(fym) + std::abs(fyp) +
              std::abs(fzm) + std::abs(fzp);
          const double denominator =
              std::max(scale, std::numeric_limits<double>::min());
          const double unscaled_metric =
              std::abs(raw_continuity) / denominator;
          const double parity_error =
              std::abs(raw_continuity +
                       true_residual.unchecked(cell, 0U)) /
              denominator;
          if (!std::isfinite(rho) || !std::isfinite(unsteady) ||
              !std::isfinite(flux_sum) || !std::isfinite(scale) ||
              !std::isfinite(unscaled_metric) ||
              !std::isfinite(parity_error)) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            continue;
          }
          local_unscaled_max =
              std::max(local_unscaled_max, unscaled_metric);
          local_operator_parity_error =
              std::max(local_operator_parity_error, parity_error);
        }
      }
    }
  }
  const double local_metrics[3U]{local_max, local_unscaled_max,
                                 local_operator_parity_error};
  double global_metrics[3U]{};
  status = reductions.checked_max({local_metrics, 3U}, {global_metrics, 3U},
                                  local);
  if (!status) return status;
  LinearConvergenceAuditResult candidate;
  candidate.metric = global_metrics[0U];
  candidate.limit = impl.continuity_tolerance;
  candidate.application_scale = alpha;
  candidate.unscaled_metric = global_metrics[1U];
  candidate.maximum_depletion = global_max_depletion;
  candidate.operator_parity_error = global_metrics[2U];
  candidate.accepted = candidate.metric <= impl.continuity_tolerance;
  candidate.terminal_rejection =
      !candidate.accepted && alpha < 1.0 &&
      candidate.unscaled_metric <= impl.continuity_tolerance &&
      candidate.operator_parity_error <= impl.continuity_tolerance;
  result = candidate;
  if (candidate.accepted) selected_alpha = alpha;
  return {};
}

Status PressureVelocityCoupler::capture_pressure_failure_provenance(
    const PressureCorrectionCertificate& pressure,
    ConstFieldView correction, ConstFieldView sealed_rhs,
    ReductionEngine& reductions,
    LinearConvergenceAuditResult& result) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }
  Impl& impl = *implementation_;
  const PressureCorrectionInput& input = impl.pressure_input;
  const Int3 cells = impl.cells;
  const bool bdf2 = input.bdf.order == 2U;
  const bool valid =
      impl.kernels != nullptr && impl.geometry != nullptr &&
      result.terminal_rejection && !result.accepted &&
      std::isfinite(result.application_scale) &&
      result.application_scale > 0.0 && result.application_scale < 1.0 &&
      std::isfinite(result.maximum_depletion) &&
      result.maximum_depletion > kPressureCorrectionSafetyFraction &&
      std::isfinite(result.metric) && std::isfinite(result.limit) &&
      result.metric > result.limit &&
      std::isfinite(result.unscaled_metric) &&
      result.unscaled_metric <= result.limit &&
      std::isfinite(result.operator_parity_error) &&
      result.operator_parity_error <= result.limit &&
      !result.failure_provenance.valid &&
      same_pressure_certificate(pressure, impl.pressure_correction) &&
      input.intermediate.corrector == 2U &&
      input.intermediate.dependency == pressure.intermediate &&
      input.intermediate.thermophysical_boundary_semantics ==
          pressure.thermophysical_boundary_semantics &&
      input.intermediate.thermophysical_boundary_target ==
          pressure.thermophysical_boundary_target &&
      input.intermediate.thermophysical_boundary_rank_local_binding ==
          pressure.thermophysical_boundary_rank_local_binding &&
      input.intermediate.thermophysical_boundary_collective_lineage ==
          pressure.thermophysical_boundary_collective_lineage &&
      input.intermediate.thermophysical_boundary_rank_local_lineage ==
          pressure.thermophysical_boundary_rank_local_lineage &&
      detail::valid_cell_view(correction, cells, 0U, 1U, 1U) &&
      detail::valid_cell_view(sealed_rhs, cells, 0U, 1U, 0U) &&
      detail::valid_cell_view(input.density_trial, cells, 0U, 1U, 0U) &&
      detail::valid_cell_view(input.density_accepted, cells, 0U, 1U, 0U) &&
      (!bdf2 || detail::valid_cell_view(input.density_previous, cells, 0U,
                                       1U, 0U)) &&
      detail::valid_cell_view(input.drho_dp_h_y, cells, 0U, 1U, 0U) &&
      impl.predecessor_c1.valid &&
      std::isfinite(impl.predecessor_c1.alpha) &&
      impl.predecessor_c1.alpha > 0.0 &&
      impl.predecessor_c1.alpha <= 1.0 &&
      std::isfinite(impl.predecessor_c1.maximum_depletion) &&
      impl.predecessor_c1.maximum_depletion >= 0.0;
  Status local = valid ? Status{}
                       : Status{StatusCode::invalid_plan, kPisoSolve};

  int local_rank = -1;
  int rank_count = 0;
  if (local && !patch_rank(impl.patch, local_rank, rank_count))
    local = {StatusCode::invalid_plan, kPisoSolve};

  bool has_limiting_cell = false;
  Int3 limiting_cell{-1, -1, -1};
  if (local) {
    // z/y/x traversal is the local restriction of the global x-fastest
    // linear order, so the first equality is the deterministic local tie.
    for (std::int32_t iz = 0; iz < cells.z && !has_limiting_cell; ++iz) {
      for (std::int32_t iy = 0; iy < cells.y && !has_limiting_cell; ++iy) {
        for (std::int32_t ix = 0; ix < cells.x; ++ix) {
          const Int3 cell{ix, iy, iz};
          const double rho = input.density_trial.unchecked(cell, 0U);
          const double psi = input.drho_dp_h_y.unchecked(cell, 0U);
          const double delta = correction.unchecked(cell, 0U);
          if (!finite_positive(rho) || !std::isfinite(psi) || psi < 0.0 ||
              !std::isfinite(delta)) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
          const double psi_delta = psi * delta;
          const double depletion = std::max(0.0, -psi_delta / rho);
          if (!std::isfinite(psi_delta) || !std::isfinite(depletion)) {
            local = {StatusCode::numerical_failure, kPisoNumerical};
            break;
          }
          if (depletion == result.maximum_depletion) {
            limiting_cell = cell;
            has_limiting_cell = true;
            break;
          }
        }
      }
    }
  }

  const double local_owner_token =
      local && has_limiting_cell
          ? static_cast<double>(rank_count - local_rank)
          : 0.0;
  double global_owner_token = 0.0;
  Status status = reductions.checked_max(
      {&local_owner_token, 1U}, {&global_owner_token, 1U}, local);
  if (!status) return status;
  if (global_owner_token < 1.0 ||
      global_owner_token > static_cast<double>(rank_count) ||
      std::trunc(global_owner_token) != global_owner_token) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }
  const int owner_rank =
      rank_count - static_cast<int>(global_owner_token);

  double local_primary[8U]{};
  double local_secondary[4U]{};
  if (local_rank == owner_rank) {
    if (!has_limiting_cell) {
      local = {StatusCode::invalid_plan, kPisoSolve};
    } else {
      const Int3 global_index{impl.patch.begin.x + limiting_cell.x,
                              impl.patch.begin.y + limiting_cell.y,
                              impl.patch.begin.z + limiting_cell.z};
      const double rho =
          input.density_trial.unchecked(limiting_cell, 0U);
      const double psi =
          input.drho_dp_h_y.unchecked(limiting_cell, 0U);
      const double delta = correction.unchecked(limiting_cell, 0U);
      const double depletion = std::max(0.0, -(psi * delta) / rho);
      const double rho_n =
          input.density_accepted.unchecked(limiting_cell, 0U);
      const double rho_nm1 =
          bdf2 ? input.density_previous.unchecked(limiting_cell, 0U) : 0.0;
      const double volume =
          detail::cell_volume(*impl.kernels, limiting_cell);
      const double storage =
          volume * (input.bdf.a0 * rho + input.bdf.a1 * rho_n +
                    input.bdf.a2 * rho_nm1);
      const ConstFaceFluxView flux =
          as_const(impl.workspace.phi_h_by_a);
      const double fxm = flux.x.unchecked(limiting_cell);
      const double fxp = flux.x.unchecked(
          {limiting_cell.x + 1, limiting_cell.y, limiting_cell.z});
      const double fym = flux.y.unchecked(limiting_cell);
      const double fyp = flux.y.unchecked(
          {limiting_cell.x, limiting_cell.y + 1, limiting_cell.z});
      const double fzm = flux.z.unchecked(limiting_cell);
      const double fzp = flux.z.unchecked(
          {limiting_cell.x, limiting_cell.y, limiting_cell.z + 1});
      const double flux_divergence =
          (fxp - fxm) + (fyp - fym) + (fzp - fzm);
      const double reconstructed_rhs = -(storage + flux_divergence);
      const double rhs = sealed_rhs.unchecked(limiting_cell, 0U);
      const std::size_t activity_offset =
          static_cast<std::size_t>(limiting_cell.x) +
          static_cast<std::size_t>(cells.x) *
              (static_cast<std::size_t>(limiting_cell.y) +
               static_cast<std::size_t>(cells.y) *
                   static_cast<std::size_t>(limiting_cell.z));
      const double pressure_activity =
          impl.continuity_activity.cells.size == 0U
              ? 1.0
              : static_cast<double>(
                    impl.continuity_activity.cells.data[activity_offset]);
      const bool valid_pressure_activity =
          pressure_activity == 0.0 || pressure_activity == 1.0;
      const double expected_sealed_rhs =
          pressure_activity == 1.0 ? reconstructed_rhs : 0.0;
      const double sealed_contract_mismatch =
          std::abs(expected_sealed_rhs - rhs);
      const double roundoff_scale = std::max(
          {std::abs(storage) + std::abs(fxm) + std::abs(fxp) +
               std::abs(fym) + std::abs(fyp) + std::abs(fzm) +
               std::abs(fzp),
           std::abs(expected_sealed_rhs),
           std::abs(rhs),
           std::numeric_limits<double>::min()});
      const double roundoff_limit =
          64.0 * std::numeric_limits<double>::epsilon() * roundoff_scale;
      if (!std::isfinite(rho) || !std::isfinite(psi) ||
          !std::isfinite(delta) || !std::isfinite(depletion) ||
          !std::isfinite(storage) || !std::isfinite(flux_divergence) ||
          !std::isfinite(reconstructed_rhs) || !std::isfinite(rhs) ||
          !valid_pressure_activity ||
          !std::isfinite(sealed_contract_mismatch) ||
          sealed_contract_mismatch > roundoff_limit) {
        local = {StatusCode::numerical_failure, kPisoNumerical};
      } else {
        local_primary[0U] = static_cast<double>(global_index.x);
        local_primary[1U] = static_cast<double>(global_index.y);
        local_primary[2U] = static_cast<double>(global_index.z);
        local_primary[3U] = rho;
        local_primary[4U] = psi;
        local_primary[5U] = delta;
        local_primary[6U] = depletion;
        local_primary[7U] = storage;
        local_secondary[0U] = flux_divergence;
        local_secondary[1U] = reconstructed_rhs;
        local_secondary[2U] = rhs;
        local_secondary[3U] = pressure_activity;
      }
    }
  }

  double global_primary[8U]{};
  status = reductions.checked_sum({local_primary, 8U}, {global_primary, 8U},
                                  local);
  if (!status) return status;
  double global_secondary[4U]{};
  status = reductions.checked_sum({local_secondary, 4U},
                                  {global_secondary, 4U}, local);
  if (!status) return status;

  const auto coordinate = [](double value, std::int32_t& out) noexcept {
    if (!std::isfinite(value) || value < 0.0 ||
        value > static_cast<double>(std::numeric_limits<std::int32_t>::max()) ||
        std::trunc(value) != value) {
      return false;
    }
    out = static_cast<std::int32_t>(value);
    return true;
  };
  Int3 global_index;
  std::uint64_t global_cell = 0U;
  const Int3 global_cells = impl.geometry->global_cells();
  if (!coordinate(global_primary[0U], global_index.x) ||
      !coordinate(global_primary[1U], global_index.y) ||
      !coordinate(global_primary[2U], global_index.z) ||
      !encode_global_cell(global_index, global_cells, global_cell) ||
      global_primary[6U] != result.maximum_depletion ||
      (global_secondary[3U] != 0.0 && global_secondary[3U] != 1.0)) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }

  LinearConvergenceFailureProvenance provenance;
  provenance.valid = true;
  provenance.global_cell = global_cell;
  provenance.global_index = global_index;
  provenance.owner_rank = owner_rank;
  provenance.pressure_activity =
      static_cast<std::uint8_t>(global_secondary[3U]);
  provenance.predecessor_application_scale = impl.predecessor_c1.alpha;
  provenance.predecessor_maximum_depletion =
      impl.predecessor_c1.maximum_depletion;
  provenance.density = global_primary[3U];
  provenance.psi = global_primary[4U];
  provenance.raw_correction = global_primary[5U];
  provenance.depletion = global_primary[6U];
  provenance.bdf_storage = global_primary[7U];
  provenance.flux_divergence = global_secondary[0U];
  provenance.reconstructed_rhs = global_secondary[1U];
  provenance.sealed_rhs = global_secondary[2U];
  provenance.rhs_absolute_mismatch =
      std::abs(provenance.reconstructed_rhs - provenance.sealed_rhs);
  provenance.rhs_relative_mismatch =
      provenance.rhs_absolute_mismatch /
      std::max({std::abs(provenance.reconstructed_rhs),
                std::abs(provenance.sealed_rhs),
                std::numeric_limits<double>::min()});
  result.failure_provenance = provenance;
  return {};
}

PisoPressureSolveEpoch::~PisoPressureSolveEpoch() noexcept {
  discard_workspace();
}

void PisoPressureSolveEpoch::discard_workspace() noexcept {
  SolverWorkspace* const prepared_workspace = prepared_.workspace;
  SolverWorkspace* const epoch_workspace = epoch_workspace_;
  if (workspace_ != nullptr) workspace_->recycle_clear();
  if (prepared_workspace != nullptr && prepared_workspace != workspace_)
    prepared_workspace->recycle_clear();
  if (epoch_workspace != nullptr && epoch_workspace != workspace_ &&
      epoch_workspace != prepared_workspace)
    epoch_workspace->recycle_clear();
  epoch_coupler_ = nullptr;
  epoch_workspace_ = nullptr;
  workspace_ = nullptr;
  epoch_communicator_ = MPI_COMM_NULL;
  epoch_rank_ = -1;
  epoch_size_ = 0;
  refinement_target_generation_ = 0U;
  prepared_ = {};
}

Status PisoPressureSolveEpoch::begin(const PisoPlan& plan) noexcept {
  const bool valid_restart =
      plan.pressure_algorithm() == LinearAlgorithm::fgmres
          ? plan.pressure_solve().restart != 0U
          : plan.pressure_algorithm() == LinearAlgorithm::bicgstab &&
                plan.pressure_solve().restart == 0U;
  if (active_ || epoch_coupler_ != nullptr || epoch_workspace_ != nullptr ||
      workspace_ != nullptr || epoch_communicator_ != MPI_COMM_NULL ||
      epoch_rank_ != -1 || epoch_size_ != 0 || prepared_.valid ||
      refinement_target_generation_ != 0U || plan.fingerprint() == 0U ||
      !valid_restart) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }
  results_ = {};
  refinement_results_ = {};
  epoch_coupler_ = nullptr;
  epoch_workspace_ = nullptr;
  epoch_communicator_ = MPI_COMM_NULL;
  epoch_rank_ = -1;
  epoch_size_ = 0;
  plan_ = plan.fingerprint();
  refinement_target_generation_ = 0U;
  solve_calls_ = 0U;
  refinement_solve_calls_ = 0U;
  latest_solve_outcome_available_ = false;
  latest_solve_lowest_failing_rank_ = -1;
  failed_ = false;
  active_ = true;
  return {};
}

Status PisoPressureSolveEpoch::prepare_linear_lifecycle(
    const PisoPlan& plan, std::uint8_t corrector,
    const PressureCorrectionCertificate& pressure,
    LinearIdentity identity,
    MgCoefficientIdentity coefficients,
    PressureVelocityCoupler& coupler,
    PressureLinearOperator& lifecycle_operator,
    NativeCartesianMgPlan& preconditioner,
    PressureCorrectionSystemView system,
    FieldView correction,
    SolverWorkspace& workspace,
    MgPlanCounters* mg_counters) noexcept {
  const auto reject = [&](Status status) noexcept {
    if (coupler.implementation_ != nullptr)
      coupler.implementation_->sealed = {};
    discard_workspace();
    workspace.recycle_clear();
    active_ = false;
    failed_ = true;
    return status;
  };
  const std::uint8_t expected = static_cast<std::uint8_t>(solve_calls_ + 1U);
  if (!active_ || failed_ || prepared_.valid ||
      plan.fingerprint() != plan_ || corrector != expected ||
      corrector > 2U || !pressure.valid() ||
      pressure.corrector != corrector ||
      pressure.pressure_energy_refinement != 0U ||
      pressure.pressure_energy_refinement_collective_lineage != 0U ||
      pressure.pressure_energy_refinement_lineage != 0U ||
      (corrector == 1U
           ? epoch_coupler_ != nullptr || epoch_workspace_ != nullptr ||
                 epoch_communicator_ != MPI_COMM_NULL || epoch_rank_ != -1 ||
                 epoch_size_ != 0
           : epoch_coupler_ != &coupler || epoch_workspace_ != &workspace ||
                 workspace_ != &workspace) ||
      !detail::valid_cell_view(correction, system.rhs.interior, 0U, 1U)) {
    return reject({StatusCode::invalid_plan, kPisoSolve});
  }
  Status status = coupler.refresh_pressure_linear_lifecycle(
      pressure, identity, coefficients, system, lifecycle_operator,
      preconditioner, mg_counters);
  if (!status) return reject(status);

  const LinearOperatorCertificate lifecycle_certificate =
      lifecycle_operator.certificate();
  const LinearPreconditionerCertificate preconditioner_certificate =
      preconditioner.certificate();
  if (!same_linear_identity(lifecycle_certificate.identity, identity) ||
      lifecycle_certificate.collective_fingerprint == 0U ||
      !same_cells(lifecycle_certificate.local_shape, system.rhs.interior) ||
      !same_linear_identity(preconditioner_certificate.identity, identity) ||
      preconditioner_certificate.collective_fingerprint == 0U) {
    return reject({StatusCode::invalid_plan, kPisoSolve});
  }

  prepared_.plan = &plan;
  prepared_.pressure = pressure;
  prepared_.identity = identity;
  prepared_.coupler = &coupler;
  prepared_.lifecycle_operator = &lifecycle_operator;
  prepared_.preconditioner = &preconditioner;
  prepared_.system = system;
  prepared_.correction = correction;
  prepared_.workspace = &workspace;
  prepared_.lifecycle_certificate = lifecycle_certificate;
  prepared_.preconditioner_certificate = preconditioner_certificate;
  prepared_.corrector = corrector;
  if (corrector == 1U) {
    const PressureVelocityCoupler::Impl& impl = *coupler.implementation_;
    epoch_coupler_ = &coupler;
    epoch_workspace_ = &workspace;
    epoch_communicator_ = impl.communicator;
    epoch_rank_ = impl.rank;
    epoch_size_ = impl.size;
  } else {
    refinement_target_generation_ = pressure.time;
  }
  prepared_.valid = true;
  return {};
}

Status PisoPressureSolveEpoch::
    prepare_pressure_energy_refinement_lifecycle(
        const PisoPlan& plan, std::uint8_t refinement_iteration,
        const PressureCorrectionCertificate& pressure,
        LinearIdentity identity, MgCoefficientIdentity coefficients,
        PressureVelocityCoupler& coupler,
        PressureLinearOperator& lifecycle_operator,
        NativeCartesianMgPlan& preconditioner,
        PressureCorrectionSystemView system, FieldView correction,
        SolverWorkspace& workspace, MgPlanCounters* mg_counters) noexcept {
  const auto reject = [&](Status status) noexcept {
    if (coupler.implementation_ != nullptr)
      coupler.implementation_->sealed = {};
    discard_workspace();
    workspace.recycle_clear();
    active_ = false;
    failed_ = true;
    return status;
  };
  const std::uint8_t expected =
      static_cast<std::uint8_t>(refinement_solve_calls_ + 1U);
  const bool epoch_collective_available =
      epoch_communicator_ != MPI_COMM_NULL && epoch_rank_ >= 0 &&
      epoch_size_ > 0 && epoch_rank_ < epoch_size_;
  Status local;
  bool prefix_valid = valid_pressure_energy_refinement_prefix(
      refinement_results_, refinement_solve_calls_);
  bool lineage_unique = true;
  for (std::uint8_t index = 0U; index < refinement_solve_calls_; ++index)
    lineage_unique =
        lineage_unique &&
        refinement_results_[index].collective_lineage !=
            pressure.pressure_energy_refinement_collective_lineage;
  if (!active_ || failed_ || prepared_.valid || solve_calls_ != 2U ||
      plan.fingerprint() != plan_ || refinement_iteration != expected ||
      refinement_iteration > kPressureEnergyRefinementCapacity ||
      !pressure.valid() || pressure.corrector != 2U ||
      refinement_target_generation_ == 0U ||
      pressure.time != refinement_target_generation_ || !prefix_valid ||
      !lineage_unique ||
      pressure.pressure_energy_refinement != refinement_iteration ||
      pressure.pressure_energy_refinement_collective_lineage == 0U ||
      pressure.pressure_energy_refinement_lineage == 0U ||
      epoch_coupler_ != &coupler ||
      epoch_workspace_ != &workspace || workspace_ != nullptr ||
      coupler.implementation_ == nullptr ||
      !detail::valid_cell_view(correction, system.rhs.interior, 0U, 1U)) {
    local = {StatusCode::invalid_plan, kPisoSolve};
  }
  if (!epoch_collective_available)
    return reject({StatusCode::invalid_plan, kPisoSolve});
  if (coupler.implementation_ != nullptr &&
      (coupler.implementation_->communicator != epoch_communicator_ ||
       coupler.implementation_->rank != epoch_rank_ ||
       coupler.implementation_->size != epoch_size_)) {
    local = {StatusCode::invalid_plan, kPisoSolve};
  }
  int preflight_lowest = -1;
  Status status = collective_status(epoch_communicator_, local, epoch_rank_,
                                    epoch_size_, preflight_lowest);
  if (!status) return reject(status);
  PressureVelocityCoupler::Impl& impl = *coupler.implementation_;

  workspace.recycle_clear();
  status = coupler.refresh_pressure_linear_lifecycle(
      pressure, identity, coefficients, system, lifecycle_operator,
      preconditioner, mg_counters, true);
  if (!status) return reject(status);

  const LinearOperatorCertificate lifecycle_certificate =
      lifecycle_operator.certificate();
  const LinearPreconditionerCertificate preconditioner_certificate =
      preconditioner.certificate();
  local = {};
  if (!same_linear_identity(lifecycle_certificate.identity, identity) ||
      lifecycle_certificate.collective_fingerprint == 0U ||
      !same_cells(lifecycle_certificate.local_shape, system.rhs.interior) ||
      !same_linear_identity(preconditioner_certificate.identity, identity) ||
      preconditioner_certificate.collective_fingerprint == 0U) {
    local = {StatusCode::invalid_plan, kPisoSolve};
  }
  int certificate_lowest = -1;
  status = collective_status(impl.communicator, local, impl.rank, impl.size,
                             certificate_lowest);
  if (!status) return reject(status);

  prepared_.plan = &plan;
  prepared_.pressure = pressure;
  prepared_.identity = identity;
  prepared_.coupler = &coupler;
  prepared_.lifecycle_operator = &lifecycle_operator;
  prepared_.preconditioner = &preconditioner;
  prepared_.system = system;
  prepared_.correction = correction;
  prepared_.workspace = &workspace;
  prepared_.lifecycle_certificate = lifecycle_certificate;
  prepared_.preconditioner_certificate = preconditioner_certificate;
  prepared_.corrector = 2U;
  prepared_.refinement_iteration = refinement_iteration;
  prepared_.pressure_energy_refinement = true;
  prepared_.valid = true;
  return {};
}

Status PisoPressureSolveEpoch::solve_prepared(
    LinearOperator& exact_operator,
    PisoPressureSolveContract contract,
    ReductionEngine& reductions,
    ResourceCounters* resources) noexcept {
  const LinearSolveControl control =
      prepared_.plan == nullptr ? LinearSolveControl{}
                                : prepared_.plan->pressure_solve();
  return solve_prepared(exact_operator, contract, control, reductions,
                        resources);
}

Status PisoPressureSolveEpoch::solve_prepared(
    LinearOperator& exact_operator,
    PisoPressureSolveContract contract,
    const LinearSolveControl& solve_control,
    ReductionEngine& reductions,
    ResourceCounters* resources) noexcept {
  // This narrow authority describes only the LinearSolveResult produced by
  // the current invocation.  Early lifecycle rejection deliberately leaves
  // it unavailable; a recorded solve with no specific failing rank publishes
  // available=true, rank=-1.
  latest_solve_outcome_available_ = false;
  latest_solve_lowest_failing_rank_ = -1;
  const bool valid_contract =
      contract == PisoPressureSolveContract::pressure_continuity ||
      contract == PisoPressureSolveContract::continuity_energy_coupled;
  if (!active_ || failed_ || !prepared_.valid || !valid_contract) {
    discard_workspace();
    active_ = false;
    failed_ = true;
    return {StatusCode::invalid_plan, kPisoSolve};
  }

  const PreparedLinearLifecycle prepared = prepared_;
  const bool pressure_energy_refinement =
      prepared.pressure_energy_refinement;
  const auto reject = [&](Status status) noexcept {
    if (prepared.coupler != nullptr &&
        prepared.coupler->implementation_ != nullptr)
      prepared.coupler->implementation_->sealed = {};
    discard_workspace();
    active_ = false;
    failed_ = true;
    return status;
  };
  const bool sequence_valid =
      pressure_energy_refinement
          ? (solve_calls_ == 2U &&
             prepared.corrector == 2U &&
             prepared.refinement_iteration ==
                 static_cast<std::uint8_t>(refinement_solve_calls_ + 1U) &&
             prepared.refinement_iteration <=
                 kPressureEnergyRefinementCapacity &&
             prepared.pressure.pressure_energy_refinement ==
                 prepared.refinement_iteration &&
             prepared.pressure.time == refinement_target_generation_ &&
             prepared.pressure
                     .pressure_energy_refinement_collective_lineage != 0U &&
             contract ==
                 PisoPressureSolveContract::continuity_energy_coupled)
          : (prepared.corrector ==
                 static_cast<std::uint8_t>(solve_calls_ + 1U) &&
             prepared.corrector <= 2U &&
             prepared.refinement_iteration == 0U &&
             prepared.pressure.pressure_energy_refinement == 0U &&
             prepared.pressure
                     .pressure_energy_refinement_collective_lineage == 0U);
  bool solve_control_valid = false;
  if (prepared.plan != nullptr) {
    const LinearSolveControl& base = prepared.plan->pressure_solve();
    const bool fixed_control_matches =
        solve_control.absolute_tolerance == base.absolute_tolerance &&
        solve_control.maximum_iterations == base.maximum_iterations &&
        solve_control.true_residual_interval ==
            base.true_residual_interval &&
        solve_control.restart == base.restart;
    const bool exact_control =
        solve_control.relative_tolerance == base.relative_tolerance;
    const double guarded_ceiling =
        contract == PisoPressureSolveContract::pressure_continuity &&
                prepared.plan->coupling() == CouplingKind::simple &&
                prepared.corrector == 1U && !pressure_energy_refinement
            ? kSimplePressureInexactForcingRelativeToleranceCeiling
            : kPressureInexactForcingRelativeToleranceCeiling;
    const bool guarded_inexact_control =
        (contract == PisoPressureSolveContract::continuity_energy_coupled ||
         (contract == PisoPressureSolveContract::pressure_continuity &&
          prepared.plan->coupling() == CouplingKind::simple &&
          prepared.corrector == 1U && !pressure_energy_refinement)) &&
        base.relative_tolerance <
            guarded_ceiling &&
        solve_control.relative_tolerance >= base.relative_tolerance &&
        solve_control.relative_tolerance <= guarded_ceiling;
    solve_control_valid =
        fixed_control_matches && (exact_control || guarded_inexact_control) &&
        valid_control(prepared.plan->pressure_algorithm(),
                      prepared.plan->pressure_correction_scaling(),
                      solve_control);
  }
  if (prepared.plan == nullptr || prepared.coupler == nullptr ||
      prepared.coupler->implementation_ == nullptr ||
      prepared.lifecycle_operator == nullptr ||
      prepared.preconditioner == nullptr || prepared.workspace == nullptr ||
      epoch_workspace_ != prepared.workspace || !sequence_valid ||
      !solve_control_valid ||
      !same_pressure_certificate(prepared.pressure,
                                 prepared_.pressure) ||
      !same_linear_identity(prepared.identity, prepared_.identity) ||
      !same_field_identity(as_const(prepared.system.diagonal),
                           as_const(prepared_.system.diagonal)) ||
      !same_field_identity(as_const(prepared.system.rhs),
                           as_const(prepared_.system.rhs)) ||
      !same_field_identity(as_const(prepared.correction),
                           as_const(prepared_.correction))) {
    return reject({StatusCode::invalid_plan, kPisoSolve});
  }

  const LinearOperatorCertificate lifecycle_now =
      prepared.lifecycle_operator->certificate();
  const LinearPreconditionerCertificate preconditioner_now =
      prepared.preconditioner->certificate();
  const LinearOperatorCertificate exact = exact_operator.certificate();
  if (!same_linear_operator_certificate(lifecycle_now,
                                        prepared.lifecycle_certificate) ||
      !same_linear_identity(preconditioner_now.identity,
                            prepared.preconditioner_certificate.identity) ||
      preconditioner_now.collective_fingerprint !=
          prepared.preconditioner_certificate.collective_fingerprint ||
      !same_linear_identity(exact.identity, prepared.identity) ||
      exact.collective_fingerprint == 0U ||
      !same_cells(exact.local_shape, prepared.system.rhs.interior)) {
    return reject({StatusCode::invalid_plan, kPisoSolve});
  }

  PisoPlan const& plan = *prepared.plan;
  PressureVelocityCoupler& coupler = *prepared.coupler;
  NativeCartesianMgPlan& preconditioner = *prepared.preconditioner;
  PressureCorrectionSystemView system = prepared.system;
  FieldView correction = prepared.correction;
  SolverWorkspace& workspace = *prepared.workspace;
  const std::uint8_t corrector = prepared.corrector;
  const LinearIdentity identity = prepared.identity;
  const PressureCorrectionCertificate pressure = prepared.pressure;
  const bool bicgstab = plan.pressure_algorithm() == LinearAlgorithm::bicgstab;
  const bool independent_simple_sweep =
      plan.coupling() == CouplingKind::simple && !pressure_energy_refinement;
  const bool pressure_continuity_contract =
      contract == PisoPressureSolveContract::pressure_continuity;
  if (!pressure_continuity_contract && coupler.implementation_ != nullptr)
    coupler.implementation_->sealed = {};

  Status status{};
  if (bicgstab || independent_simple_sweep) {
    // BiCGStab has no Arnoldi recycle span. SIMPLE also deliberately breaks
    // the frozen-momentum C1->C2 projection lineage because a fresh momentum
    // system is assembled between its pressure solves. Clear any stale span
    // while retaining the prepared workspace identity across both solves.
    if (corrector == 1U)
      workspace_ = &workspace;
    workspace.recycle_clear();
  } else {
    status = pressure_energy_refinement
                 ? workspace.recycle_begin_capture(system.rhs.interior,
                                                   identity.fingerprint)
                 : (corrector == 1U
                 ? workspace.recycle_begin_capture(system.rhs.interior,
                                                   identity.fingerprint)
                 : workspace.recycle_begin_projection(system.rhs.interior,
                                                      identity.fingerprint));
    if (!status) return reject(status);
    if (corrector == 1U)
      workspace_ = &workspace;
  }
  PressureContinuityConvergenceAudit convergence_audit(
      coupler, pressure, as_const(system.rhs));
  LinearConvergenceAudit* const selected_audit =
      pressure_continuity_contract && corrector == 2U
          ? &convergence_audit
          : nullptr;
  const LinearSolveInvocation invocation{
      as_const(system.rhs), correction, identity, solve_control,
      selected_audit};
  LinearSolveResult result =
      bicgstab ? solve_bicgstab(exact_operator, preconditioner, invocation,
                                workspace, reductions, resources)
               : solve_fgmres(exact_operator, preconditioner, invocation,
                              workspace, reductions, resources);
  latest_solve_outcome_available_ = true;
  latest_solve_lowest_failing_rank_ = result.lowest_failing_rank;
  const auto record_refinement_result = [&](LinearSolveResult solve) noexcept {
    PisoPressureEnergyRefinementSolveReport recorded;
    recorded.solve = solve;
    recorded.target_generation = prepared.pressure.time;
    recorded.collective_lineage =
        prepared.pressure.pressure_energy_refinement_collective_lineage;
    recorded.pressure_state = prepared.pressure.state;
    recorded.linear_identity = prepared.identity;
    recorded.ordinal = prepared.refinement_iteration;
    assert(recorded.valid());
    refinement_results_[refinement_solve_calls_] = recorded;
    ++refinement_solve_calls_;
  };
  const bool accepted =
      static_cast<bool>(result.status) &&
      (result.termination == LinearTermination::converged ||
       result.termination == LinearTermination::zero_rhs);
  if (!accepted) {
    // Preserve the failed corrector's already-computed linear provenance.
    // The attempt remains rejected and no correction is applied, but observe()
    // can now distinguish iteration exhaustion, recurrence breakdown and a
    // convergence-audit rejection without another operator application.
    if (pressure_energy_refinement) {
      record_refinement_result(result);
    } else {
      results_[solve_calls_] = result;
      ++solve_calls_;
    }
    return reject(result.status.code == StatusCode::ok
                      ? Status{StatusCode::rejected_step, kPisoSolve}
                      : result.status);
  }
  if (corrector == 2U && pressure_continuity_contract) {
    PressureContinuityConvergenceAudit::Provisional provisional;
    const bool provisional_valid =
        convergence_audit.take_provisional(provisional);
    const LinearOperatorCertificate exact_after = exact_operator.certificate();
    const bool exact_complete =
        same_linear_operator_certificate(exact_after, exact) &&
        exact_after.identity.fingerprint != 0U &&
        exact_after.collective_fingerprint != 0U &&
        same_linear_identity(exact_after.identity, identity) &&
        same_cells(exact_after.local_shape, system.rhs.interior);
    const bool audited_workspace_valid =
        detail::valid_cell_view(provisional.audited_workspace,
                                system.rhs.interior, 0U, 1U, 0U) &&
        provisional.audited_workspace.base != nullptr;
    if (!provisional_valid ||
        !same_pressure_certificate(provisional.pressure, pressure) ||
        !std::isfinite(provisional.alpha) || provisional.alpha <= 0.0 ||
        provisional.alpha > 1.0 || !exact_complete ||
        !audited_workspace_valid ||
        !detail::valid_cell_view(as_const(correction), system.rhs.interior,
                                 0U, 1U, 0U)) {
      return reject({StatusCode::invalid_plan, kPisoSolve});
    }
    PressureVelocityCoupler::Impl& impl = *coupler.implementation_;
    impl.sealed.valid = true;
    impl.sealed.pressure = pressure;
    impl.sealed.exact_operator = exact_after;
    impl.sealed.correction = as_const(correction);
    impl.sealed.alpha = provisional.alpha;
  } else if (!pressure_continuity_contract) {
    const LinearOperatorCertificate exact_after = exact_operator.certificate();
    if (!same_linear_operator_certificate(exact_after, exact))
      return reject({StatusCode::invalid_plan, kPisoSolve});
    if (coupler.implementation_ != nullptr)
      coupler.implementation_->sealed = {};
  }
  if (pressure_energy_refinement) {
    result.recycle_cycle_corrections =
        workspace.recycle_capture_cycle_corrections();
    result.recycle_capture_vector_passes =
        workspace.recycle_capture_vector_passes();
    result.recycle_capture_cycle_attempts =
        workspace.recycle_capture_cycle_attempts();
    result.recycle_capture_reduction_calls =
        workspace.recycle_capture_reduction_calls();
    result.recycle_capture_blocking_operations =
        workspace.recycle_capture_blocking_operations();
    workspace.recycle_clear();
    workspace_ = nullptr;
  } else if (bicgstab || independent_simple_sweep) {
    result.recycle_cycle_corrections = 0U;
    result.recycle_capture_vector_passes = 0U;
    result.recycle_capture_cycle_attempts = 0U;
    result.recycle_capture_reduction_calls = 0U;
    result.recycle_capture_blocking_operations = 0U;
    workspace.recycle_clear();
    workspace_ = corrector == 2U ? nullptr : workspace_;
  } else if (corrector == 1U) {
    result.recycle_cycle_corrections =
        workspace.recycle_capture_cycle_corrections();
    result.recycle_capture_vector_passes =
        workspace.recycle_capture_vector_passes();
    result.recycle_capture_cycle_attempts =
        workspace.recycle_capture_cycle_attempts();
    result.recycle_capture_reduction_calls =
        workspace.recycle_capture_reduction_calls();
    result.recycle_capture_blocking_operations =
        workspace.recycle_capture_blocking_operations();
  } else {
    result.recycle_cycle_corrections = 0U;
    result.recycle_capture_vector_passes = 0U;
    result.recycle_capture_cycle_attempts = 0U;
    result.recycle_capture_reduction_calls = 0U;
    result.recycle_capture_blocking_operations = 0U;
    workspace.recycle_clear();
    workspace_ = nullptr;
  }
  if (pressure_energy_refinement) {
    record_refinement_result(result);
  } else {
    results_[solve_calls_] = result;
    ++solve_calls_;
  }
  prepared_ = {};
  return {};
}

Status PisoPressureSolveEpoch::record_stationary(
    const PisoPlan& plan,
    const PressureCorrectionCertificate& pressure,
    const PressureEnergyStationaryCertificate& stationary,
    PressureVelocityCoupler& coupler) noexcept {
  const auto reject = [&](Status status) noexcept {
    if (coupler.implementation_ != nullptr)
      coupler.implementation_->sealed = {};
    discard_workspace();
    active_ = false;
    failed_ = true;
    return status;
  };
  const bool current =
      active_ && !failed_ && !prepared_.valid &&
      plan.fingerprint() == plan_ &&
      plan.coupling() == CouplingKind::simple && solve_calls_ == 1U &&
      refinement_solve_calls_ == 0U && pressure.valid() &&
      pressure.corrector == 2U &&
      pressure.pressure_energy_refinement == 0U && stationary.valid() &&
      stationary.issuer_ == &coupler && stationary.corrector_ == 2U &&
      stationary.target_time_ == pressure.time &&
      coupler.implementation_ != nullptr && epoch_coupler_ == &coupler &&
      epoch_workspace_ != nullptr && workspace_ == epoch_workspace_ &&
      epoch_communicator_ == coupler.implementation_->communicator &&
      epoch_rank_ == coupler.implementation_->rank &&
      epoch_size_ == coupler.implementation_->size &&
      stationary.exact_lineage_ ==
          coupler.implementation_->current_frozen_exact_lineage &&
      same_pressure_certificate(
          pressure, coupler.implementation_->current_pressure_work);
  if (!current)
    return reject({StatusCode::invalid_plan, kPisoSolve});

  LinearSolveResult result;
  result.status = {};
  result.termination = LinearTermination::zero_rhs;
  result.lowest_failing_rank = -1;
  results_[solve_calls_] = result;
  ++solve_calls_;
  latest_solve_outcome_available_ = true;
  latest_solve_lowest_failing_rank_ = -1;
  workspace_->recycle_clear();
  workspace_ = nullptr;
  return {};
}

Status PisoPressureSolveEpoch::solve(
    const PisoPlan& plan, std::uint8_t corrector,
    const PressureCorrectionCertificate& pressure,
    LinearIdentity identity,
    MgCoefficientIdentity coefficients,
    PressureVelocityCoupler& coupler,
    PressureLinearOperator& linear_operator,
    NativeCartesianMgPlan& preconditioner,
    PressureCorrectionSystemView system,
    FieldView correction,
    SolverWorkspace& workspace,
    ReductionEngine& reductions,
    ResourceCounters* resources,
    MgPlanCounters* mg_counters) noexcept {
  return solve(plan, corrector, pressure, identity, coefficients, coupler,
               linear_operator, linear_operator, preconditioner, system,
               correction, workspace, reductions, resources, mg_counters);
}

Status PisoPressureSolveEpoch::solve(
    const PisoPlan& plan, std::uint8_t corrector,
    const PressureCorrectionCertificate& pressure,
    LinearIdentity identity,
    MgCoefficientIdentity coefficients,
    PressureVelocityCoupler& coupler,
    PressureLinearOperator& lifecycle_operator,
    LinearOperator& exact_operator,
    NativeCartesianMgPlan& preconditioner,
    PressureCorrectionSystemView system,
    FieldView correction,
    SolverWorkspace& workspace,
    ReductionEngine& reductions,
    ResourceCounters* resources,
    MgPlanCounters* mg_counters) noexcept {
  Status status = prepare_linear_lifecycle(
      plan, corrector, pressure, identity, coefficients, coupler,
      lifecycle_operator, preconditioner, system, correction, workspace,
      mg_counters);
  if (!status) return status;
  return solve_prepared(exact_operator,
                        PisoPressureSolveContract::pressure_continuity,
                        reductions, resources);
}

Status PisoPressureSolveEpoch::finalize(PisoAttemptReport& report) noexcept {
  const bool refinement_prefix_valid =
      valid_pressure_energy_refinement_prefix(refinement_results_,
                                              refinement_solve_calls_);
  const bool termination_valid = valid_pressure_energy_refinement_termination(
      report.pressure_energy_refinement_termination,
      refinement_solve_calls_);
  if (!active_ || failed_ || prepared_.valid || solve_calls_ != 2U ||
      !refinement_prefix_valid || !termination_valid ||
      report.pressure_energy_refinement_termination ==
          PressureEnergyRefinementTermination::iteration_capacity_exhausted ||
      report.pressure_energy_refinement_solve_calls !=
          refinement_solve_calls_) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }
  PisoAttemptReport candidate = report;
  candidate.pressure = results_;
  candidate.pressure_energy_refinement = refinement_results_;
  candidate.pressure_solve_calls = 2U;
  candidate.pressure_energy_refinement_solve_calls =
      refinement_solve_calls_;
  report = candidate;
  epoch_coupler_ = nullptr;
  epoch_workspace_ = nullptr;
  workspace_ = nullptr;
  epoch_communicator_ = MPI_COMM_NULL;
  epoch_rank_ = -1;
  epoch_size_ = 0;
  refinement_target_generation_ = 0U;
  latest_solve_outcome_available_ = false;
  latest_solve_lowest_failing_rank_ = -1;
  active_ = false;
  return {};
}

Status PisoPressureSolveEpoch::observe(PisoAttemptReport& report) const
    noexcept {
  if (plan_ == 0U || solve_calls_ == 0U || solve_calls_ > 2U ||
      refinement_solve_calls_ > kPressureEnergyRefinementCapacity ||
      !valid_pressure_energy_refinement_prefix(refinement_results_,
                                               refinement_solve_calls_)) {
    return {StatusCode::invalid_plan, kPisoSolve};
  }
  PisoAttemptReport candidate = report;
  for (std::uint8_t index = 0U; index < solve_calls_; ++index)
    candidate.pressure[index] = results_[index];
  // A caller may reuse a failure report from an earlier, longer refinement
  // attempt.  Clear the entire typed array before publishing this epoch's
  // active prefix so no stale suffix can survive observation.
  candidate.pressure_energy_refinement = {};
  for (std::uint8_t index = 0U; index < refinement_solve_calls_; ++index)
    candidate.pressure_energy_refinement[index] =
        refinement_results_[index];
  candidate.pressure_solve_calls = solve_calls_;
  candidate.pressure_energy_refinement_solve_calls =
      refinement_solve_calls_;
  report = candidate;
  return {};
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
bool PisoPressureSolveEpoch::
    validate_pressure_energy_refinement_report_for_test(
        const PisoAttemptReport& report) noexcept {
  return valid_pressure_energy_refinement_prefix(
             report.pressure_energy_refinement,
             report.pressure_energy_refinement_solve_calls) &&
         valid_pressure_energy_refinement_termination(
             report.pressure_energy_refinement_termination,
             report.pressure_energy_refinement_solve_calls);
}
#endif

PlanFingerprint PressureVelocityCoupler::fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->fingerprint;
}

std::uintptr_t PressureVelocityCoupler::workspace_storage_address() const
    noexcept {
  return implementation_ == nullptr
             ? 0U
             : reinterpret_cast<std::uintptr_t>(
                   implementation_->workspace.r_au.base);
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
void PressureVelocityCoupler::set_frozen_stationary_tolerances_for_test(
    double continuity_tolerance, double energy_tolerance) noexcept {
  if (implementation_ == nullptr) return;
  implementation_->continuity_tolerance = continuity_tolerance;
  implementation_->energy_tolerance = energy_tolerance;
}
#endif

}  // namespace hundun::v04
