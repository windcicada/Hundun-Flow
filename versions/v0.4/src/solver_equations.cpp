// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_flow.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_equation_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kEquationPlan = 1401U;
constexpr std::uint32_t kEquationCollective = 1402U;
constexpr std::uint32_t kEquationAssembly = 1403U;
constexpr std::uint32_t kEquationNumerical = 1404U;
constexpr std::uint32_t kAssemblyEpoch = 1405U;
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
                "equation fingerprints require binary64");
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

PlanFingerprint finish_hash(std::uint64_t hash) noexcept {
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

bool equal_int3(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool valid_field_ids(const EquationPlanSpec& spec) noexcept {
  const FieldId ids[]{spec.density, spec.velocity,
                      spec.pressure_perturbation, spec.enthalpy,
                      spec.temperature, spec.effective_viscosity,
                      spec.pressure_compressibility,
                      spec.velocity_gradient};
  for (std::size_t i = 0U; i < 8U; ++i) {
    for (std::size_t j = i + 1U; j < 8U; ++j) {
      if (ids[i] == ids[j]) {
        return false;
      }
    }
  }
  return true;
}

std::size_t checked_cell_count(Int3 cells) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0) {
    return 0U;
  }
  const std::size_t x = static_cast<std::size_t>(cells.x);
  const std::size_t y = static_cast<std::size_t>(cells.y);
  const std::size_t z = static_cast<std::size_t>(cells.z);
  if (x > std::numeric_limits<std::size_t>::max() / y ||
      x * y > std::numeric_limits<std::size_t>::max() / z) {
    return 0U;
  }
  return x * y * z;
}

bool valid_scalar_catalog(Span<const ScalarEquationSpec> scalars,
                          std::size_t species_expected) noexcept {
  if (scalars.size != 0U && scalars.data == nullptr) {
    return false;
  }
  std::size_t species = 0U;
  for (std::size_t i = 0U; i < scalars.size; ++i) {
    const ScalarEquationSpec& scalar = scalars.data[i];
    if (!std::isfinite(scalar.molecular_schmidt) ||
        !std::isfinite(scalar.turbulent_schmidt) ||
        scalar.molecular_schmidt <= 0.0 ||
        scalar.turbulent_schmidt <= 0.0) {
      return false;
    }
    for (std::size_t prior = 0U; prior < i; ++prior) {
      if (scalars.data[prior].field == scalar.field) {
        return false;
      }
    }
    if (scalar.role == TransportedScalarRole::species) {
      ++species;
    } else if (scalar.role == TransportedScalarRole::passive_scalar) {
    } else {
      return false;
    }
  }
  return species == species_expected;
}

bool valid_equation_contribution_units(
    const ContributionRegistry& contributions,
    const EquationPlanSpec& spec) noexcept;

PlanFingerprint compute_semantic_fingerprint(
    const SchemePlan& schemes, const CartesianGeometryPlan& geometry,
    const BoundaryPlan& boundary, const ContributionRegistry& contributions,
    const ThermodynamicsPlan& thermodynamics, const TransportPlan& transport,
    const EquationPlanSpec& spec) noexcept {
  std::uint64_t hash = kFnvOffset;
  hash = hash_mix(hash, UINT64_C(0x7630346571756174));
  hash = hash_mix(hash, kMomentumPredictorLimiterPolicySchema);
  hash = hash_mix(hash, static_cast<std::uint8_t>(geometry.kind()));
  hash = hash_mix(hash, static_cast<std::uint32_t>(geometry.global_cells().x));
  hash = hash_mix(hash, static_cast<std::uint32_t>(geometry.global_cells().y));
  hash = hash_mix(hash, static_cast<std::uint32_t>(geometry.global_cells().z));
  hash = hash_mix(hash, double_bits(geometry.lower().x));
  hash = hash_mix(hash, double_bits(geometry.lower().y));
  hash = hash_mix(hash, double_bits(geometry.lower().z));
  hash = hash_mix(hash, double_bits(geometry.upper().x));
  hash = hash_mix(hash, double_bits(geometry.upper().y));
  hash = hash_mix(hash, double_bits(geometry.upper().z));
  const AxisMetrics* metrics[]{&geometry.x(), &geometry.y(), &geometry.z()};
  for (const AxisMetrics* metric : metrics) {
    hash = hash_mix(hash, metric->faces().size);
    for (std::size_t i = 0U; i < metric->faces().size; ++i) {
      hash = hash_mix(hash, double_bits(metric->faces().data[i]));
    }
  }
  hash = hash_mix(hash, geometry.topology_revision());
  hash = hash_mix(hash, static_cast<std::uint8_t>(schemes.momentum()));
  hash = hash_mix(hash, static_cast<std::uint8_t>(schemes.enthalpy()));
  hash = hash_mix(hash, static_cast<std::uint8_t>(schemes.species()));
  hash = hash_mix(hash, static_cast<std::uint8_t>(schemes.passive_scalar()));
  hash = hash_mix(hash, static_cast<std::uint8_t>(schemes.diffusion()));
  hash = hash_mix(hash, double_bits(schemes.limiter()));
  hash = hash_mix(hash, schemes.required_ghost_width());
  hash = hash_mix(hash, static_cast<std::uint8_t>(boundary.pressure_reference()));
  hash = hash_mix(hash, contributions.fingerprint());
  hash = hash_mix(hash, thermodynamics.fingerprint());
  hash = hash_mix(hash, transport.fingerprint());
  hash = hash_mix(hash, spec.density);
  hash = hash_mix(hash, spec.velocity);
  hash = hash_mix(hash, spec.pressure_perturbation);
  hash = hash_mix(hash, spec.enthalpy);
  hash = hash_mix(hash, spec.temperature);
  hash = hash_mix(hash, spec.effective_viscosity);
  hash = hash_mix(hash, spec.velocity_gradient);
  hash = hash_mix(hash, spec.pressure_compressibility);
  hash = hash_mix(hash, static_cast<std::uint8_t>(spec.pressure_reference));
  hash = hash_mix(hash, spec.maximum_cells_per_rank);
  hash = hash_mix(hash, spec.closed_mass_service_stage);
  hash = hash_mix(hash, spec.scalars.size);
  for (std::size_t index = 0U; index < spec.scalars.size; ++index) {
    const ScalarEquationSpec& scalar = spec.scalars.data[index];
    hash = hash_mix(hash, scalar.field);
    hash = hash_mix(hash, static_cast<std::uint8_t>(scalar.role));
    hash = hash_mix(hash, double_bits(scalar.molecular_schmidt));
    hash = hash_mix(hash, double_bits(scalar.turbulent_schmidt));
  }
  return finish_hash(hash);
}

Status local_plan_status(
    const SchemePlan& schemes, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, const BoundaryPlan& boundary,
    const ContributionRegistry& contributions,
    const ThermodynamicsPlan& thermodynamics, const TransportPlan& transport,
    const EquationPlanSpec& spec, std::size_t& local_cells) noexcept {
  local_cells = checked_cell_count(patch.cells);
  const std::size_t expected_species =
      thermodynamics.independent_species_count();
  if (schemes.fingerprint() == 0U || geometry.fingerprint() == 0U ||
      geometry.topology_revision() == 0U ||
      boundary.semantic_fingerprint() == 0U ||
      boundary.local_layout_fingerprint() == 0U || !contributions.frozen() ||
      contributions.fingerprint() == 0U ||
      thermodynamics.fingerprint() == 0U || transport.fingerprint() == 0U ||
      local_cells == 0U ||
      !equal_int3(boundary.local_cells(), patch.cells) ||
      boundary.velocity_field() != spec.velocity ||
      boundary.pressure_field() != spec.pressure_perturbation ||
      boundary.enthalpy_field() != spec.enthalpy ||
      spec.maximum_cells_per_rank == 0U ||
      local_cells > spec.maximum_cells_per_rank || !valid_field_ids(spec) ||
      boundary.pressure_reference() != spec.pressure_reference ||
      (spec.pressure_reference == PressureReferenceKind::closed_mass
           ? spec.closed_mass_service_stage == 0U
           : spec.closed_mass_service_stage != 0U) ||
      !valid_scalar_catalog(spec.scalars, expected_species) ||
      !valid_equation_contribution_units(contributions, spec)) {
    return {StatusCode::invalid_plan, kEquationPlan};
  }
  const FieldId base[]{spec.density, spec.velocity,
                       spec.pressure_perturbation, spec.enthalpy,
                       spec.temperature, spec.effective_viscosity,
                       spec.pressure_compressibility,
                       spec.velocity_gradient};
  for (std::size_t scalar = 0U; scalar < spec.scalars.size; ++scalar) {
    for (FieldId field : base) {
      if (spec.scalars.data[scalar].field == field) {
        return {StatusCode::invalid_plan, kEquationPlan};
      }
    }
  }
  return {};
}

Status collective_status(MPI_Comm communicator, Status local, int rank,
                         int size, int& lowest) noexcept {
  const int candidate = local ? size : rank;
  if (MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
      MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kEquationCollective};
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
    return {StatusCode::mpi_failure, kEquationCollective};
  }
  return {static_cast<StatusCode>(encoded >> 32U),
          static_cast<std::uint32_t>(encoded)};
}

bool compatible_history(PrimitiveHistory history, Int3 cells,
                        FieldId field, std::uint8_t components,
                        std::uint8_t trial_ghosts) noexcept {
  return history.trial.field == field && history.accepted.field == field &&
         history.previous.field == field &&
         detail::valid_cell_view(history.trial, cells, 0U, components,
                                 trial_ghosts) &&
         detail::valid_cell_view(history.accepted, cells, 0U, components,
                                 0U) &&
         detail::valid_cell_view(history.previous, cells, 0U, components,
                                 0U);
}

KernelBox resolved_box(const EquationAssemblyContext& context,
                       Int3 cells) noexcept {
  if (context.box.cells.x == 0 && context.box.cells.y == 0 &&
      context.box.cells.z == 0 && context.box.begin.x == 0 &&
      context.box.begin.y == 0 && context.box.begin.z == 0) {
    return {{0, 0, 0}, cells};
  }
  return context.box;
}

bool valid_context(const CartesianKernelPlan& kernels,
                   const EquationAssemblyContext& context,
                   KernelBox& box) noexcept {
  box = resolved_box(context, kernels.cells());
  if (!detail::bdf_matches_time_step(context.bdf, context.dt) ||
      context.time == 0U || context.geometry == 0U ||
      context.face_flux == 0U ||
      context.face_flux != context.mass_flux.revision ||
      !detail::valid_kernel_box(box, kernels.cells()) ||
      !detail::valid_flux_view(context.mass_flux, kernels.cells(),
                               context.face_flux)) {
    return false;
  }
  if (context.scope == EquationAssemblyScope::final_conservative) {
    return !context.provisional_mass_flux &&
           context.face_flux_authority != 0U &&
           context.face_flux_storage != 0U &&
           context.face_flux_revision_domain != 0U &&
           context.face_flux_authority ==
               context.mass_flux.certificate.authority() &&
           context.face_flux_storage ==
               context.mass_flux.certificate.storage() &&
           context.face_flux_revision_domain ==
               context.mass_flux.certificate.revision_domain() &&
           context.mass_flux.certificate.matches(context.mass_flux);
  }
  return context.scope == EquationAssemblyScope::momentum_predictor &&
         context.provisional_mass_flux &&
         !context.mass_flux.certificate.valid();
}

RevisionToken state_revision(const EquationStateView& state) noexcept {
  std::uint64_t hash = kFnvOffset;
  const PrimitiveHistory histories[]{state.density, state.velocity,
                                     state.pressure_perturbation,
                                     state.enthalpy, state.temperature};
  for (const PrimitiveHistory& history : histories) {
    const ConstFieldView views[]{history.trial, history.accepted,
                                 history.previous};
    for (ConstFieldView view : views) {
      hash = hash_mix(hash, view.revision);
      hash = hash_mix(hash, view.storage_identity);
      hash = hash_mix(hash, view.revision_domain);
    }
  }
  hash = hash_mix(hash, double_bits(state.pressure_reference));
  hash = hash_mix(hash, double_bits(state.accepted_pressure_reference));
  hash = hash_mix(hash, double_bits(state.previous_pressure_reference));
  return finish_hash(hash);
}

EquationAssemblyCertificate make_certificate(
    PlanFingerprint plan, const EquationStateView& state,
    const EquationAssemblyContext& context) noexcept {
  return {plan, context.scope, context.time, context.geometry,
          context.face_flux, state_revision(state), context.dt};
}

PlanFingerprint child_fingerprint(PlanFingerprint semantic,
                                  std::uint64_t tag,
                                  FieldId field = 0U) noexcept {
  return finish_hash(hash_mix(hash_mix(kFnvOffset, semantic ^ tag), field));
}

std::uint8_t convection_reach(ConvectionScheme scheme) noexcept {
  return scheme == ConvectionScheme::central2 ? 1U : 2U;
}

UnitDimension expected_contribution_units(
    FieldId field, const EquationPlanSpec& spec) noexcept {
  UnitDimension units;
  if (field == spec.velocity) {
    units.si_exponents = {1, -2, -2, 0, 0, 0, 0};
  } else if (field == spec.enthalpy) {
    units.si_exponents = {1, -1, -3, 0, 0, 0, 0};
  } else {
    units.si_exponents = {1, -3, -1, 0, 0, 0, 0};
  }
  return units;
}

bool valid_equation_contribution_units(
    const ContributionRegistry& contributions,
    const EquationPlanSpec& spec) noexcept {
  for (std::size_t index = 0U;
       index < contributions.contributions().size; ++index) {
    const CompiledContribution& contribution =
        contributions.contributions().data[index];
    bool recognized = contribution.conserved_quantity == spec.velocity ||
                      contribution.conserved_quantity == spec.enthalpy;
    for (std::size_t scalar = 0U; scalar < spec.scalars.size; ++scalar) {
      recognized = recognized ||
                   contribution.conserved_quantity ==
                       spec.scalars.data[scalar].field;
    }
    if (recognized &&
        !(contribution.units == expected_contribution_units(
                                   contribution.conserved_quantity, spec))) {
      return false;
    }
  }
  return true;
}

}  // namespace

Status AssemblyEpoch::begin(const EquationAssemblyContext& context,
                            EquationSystemView system) noexcept {
  const bool empty_box = context.box.begin.x == 0 &&
                         context.box.begin.y == 0 &&
                         context.box.begin.z == 0 &&
                         context.box.cells.x == 0 &&
                         context.box.cells.y == 0 &&
                         context.box.cells.z == 0;
  if (active_ || !empty_box ||
      !detail::bdf_matches_time_step(context.bdf, context.dt) ||
      context.time == 0U ||
      context.geometry == 0U || context.face_flux == 0U ||
      context.contribution_stage == 0U ||
      (context.scope != EquationAssemblyScope::momentum_predictor &&
       context.scope != EquationAssemblyScope::target_coupled &&
       context.scope != EquationAssemblyScope::final_conservative)) {
    return {StatusCode::invalid_plan, kAssemblyEpoch};
  }
  context_ = context;
  system_ = system;
  candidate_ = {};
  cells_ = {};
  covered_cells_ = 0U;
  tile_count_ = 0U;
  failure_ = {};
  active_ = true;
  return {};
}

Status AssemblyEpoch::fail(Status status) noexcept {
  if (!active_) {
    return {StatusCode::invalid_plan, kAssemblyEpoch};
  }
  if (status) {
    status = {StatusCode::invalid_plan, kAssemblyEpoch};
  }
  if (failure_) {
    failure_ = status;
  }
  return status;
}

Status AssemblyEpoch::record(
    KernelBox box, Int3 cells,
    const EquationAssemblyCertificate& certificate) noexcept {
  if (!active_ || !failure_ || !certificate.valid() ||
      !detail::valid_kernel_box(box, cells)) {
    return fail({StatusCode::invalid_plan, kAssemblyEpoch});
  }
  if (tile_count_ == 0U) {
    cells_ = cells;
    candidate_ = certificate;
  } else {
    const bool same_cells = cells_.x == cells.x && cells_.y == cells.y &&
                            cells_.z == cells.z;
    const bool same_certificate =
        candidate_.plan == certificate.plan &&
        candidate_.scope == certificate.scope &&
        candidate_.time == certificate.time &&
        candidate_.geometry == certificate.geometry &&
        candidate_.face_flux == certificate.face_flux &&
        candidate_.state == certificate.state &&
        candidate_.dt == certificate.dt;
    if (!same_cells || !same_certificate) {
      return fail({StatusCode::invalid_plan, kAssemblyEpoch});
    }
  }
  if (tile_count_ == maximum_tiles) {
    return fail({StatusCode::allocation_failure, kAssemblyEpoch});
  }
  const auto overlaps = [](TileRecord left, KernelBox right) noexcept {
    const Int3 left_end{left.begin_x + left.cells_x,
                        left.begin_y + left.cells_y,
                        left.begin_z + left.cells_z};
    const Int3 right_end{right.begin.x + right.cells.x,
                         right.begin.y + right.cells.y,
                         right.begin.z + right.cells.z};
    return left.begin_x < right_end.x && right.begin.x < left_end.x &&
           left.begin_y < right_end.y && right.begin.y < left_end.y &&
           left.begin_z < right_end.z && right.begin.z < left_end.z;
  };
  for (std::size_t index = 0U; index < tile_count_; ++index) {
    if (overlaps(tiles_[index], box)) {
      return fail({StatusCode::invalid_plan, kAssemblyEpoch});
    }
  }
  const std::uint64_t count = detail::box_cell_count(box);
  if (count == 0U ||
      covered_cells_ > std::numeric_limits<std::uint64_t>::max() - count) {
    return fail({StatusCode::invalid_plan, kAssemblyEpoch});
  }
  tiles_[tile_count_++] = {box.begin.x, box.begin.y, box.begin.z,
                           box.cells.x, box.cells.y, box.cells.z};
  covered_cells_ += count;
  return {};
}

Status AssemblyEpoch::finalize(
    EquationAssemblyCertificate& certificate) noexcept {
  if (!active_) {
    return {StatusCode::invalid_plan, kAssemblyEpoch};
  }
  active_ = false;
  if (!failure_) {
    return failure_;
  }
  const KernelBox full{{0, 0, 0}, cells_};
  const std::uint64_t expected = detail::box_cell_count(full);
  if (tile_count_ == 0U || expected == 0U || covered_cells_ != expected ||
      !candidate_.valid()) {
    return {StatusCode::invalid_plan, kAssemblyEpoch};
  }
  certificate = candidate_;
  return {};
}

void EquationPlanSet::reset() noexcept {
  kernels_ = CartesianKernelPlan{};
  continuity_.kernels_ = nullptr;
  continuity_.fingerprint_ = 0U;
  momentum_.kernels_ = nullptr;
  momentum_.contributions_.clear();
  momentum_.fingerprint_ = 0U;
  enthalpy_.kernels_ = nullptr;
  enthalpy_.contributions_.clear();
  enthalpy_.fingerprint_ = 0U;
  species_.kernels_ = nullptr;
  species_.specs_.clear();
  species_.contribution_counts_.clear();
  species_.contribution_begins_.clear();
  species_.contributions_.clear();
  species_.fingerprint_ = 0U;
  scalars_.kernels_ = nullptr;
  scalars_.specs_.clear();
  scalars_.contribution_counts_.clear();
  scalars_.contribution_begins_.clear();
  scalars_.contributions_.clear();
  scalars_.fingerprint_ = 0U;
  pressure_reference_ = {};
  thermophysical_predictor_.kernels_ = nullptr;
  thermophysical_predictor_.cells_ = {};
  thermophysical_predictor_.patch_begin_ = {};
  thermophysical_predictor_.density_ = 0U;
  thermophysical_predictor_.enthalpy_ = 0U;
  thermophysical_predictor_.species_.clear();
  thermophysical_predictor_.passive_scalars_.clear();
  thermophysical_predictor_.species_enthalpy_minimum_.clear();
  thermophysical_predictor_.species_enthalpy_maximum_.clear();
  thermophysical_predictor_.dependent_enthalpy_minimum_ = 0.0;
  thermophysical_predictor_.dependent_enthalpy_maximum_ = 0.0;
  thermophysical_predictor_.enthalpy_convection_ =
      ConvectionScheme::limited_central2;
  thermophysical_predictor_.species_convection_ = ConvectionScheme::tvd2;
  thermophysical_predictor_.passive_scalar_convection_ =
      ConvectionScheme::tvd2;
  thermophysical_predictor_.enthalpy_reach_ = 2U;
  thermophysical_predictor_.species_reach_ = 2U;
  thermophysical_predictor_.passive_scalar_reach_ = 2U;
  thermophysical_predictor_.boundary_ = nullptr;
  thermophysical_predictor_.geometry_revision_ = 0U;
  thermophysical_predictor_.boundary_revision_ = 0U;
  thermophysical_predictor_.transport_fingerprint_ = 0U;
  thermophysical_predictor_.fingerprint_ = 0U;
  global_cells_ = {};
  patch_begin_ = {};
  local_cells_ = 0U;
  semantic_fingerprint_ = 0U;
  thermodynamics_fingerprint_ = 0U;
  transport_fingerprint_ = 0U;
  fingerprint_ = 0U;
}

void EquationPlanSet::rebind() noexcept {
  CartesianKernelPlan* const kernels =
      kernels_.fingerprint() == 0U ? nullptr : &kernels_;
  continuity_.kernels_ = kernels;
  momentum_.kernels_ = kernels;
  enthalpy_.kernels_ = kernels;
  species_.kernels_ = kernels;
  scalars_.kernels_ = kernels;
  pressure_reference_.kernels_ = kernels;
  thermophysical_predictor_.kernels_ = kernels;
}

void EquationPlanSet::move_from(EquationPlanSet&& other) noexcept {
  kernels_ = std::move(other.kernels_);
  continuity_.cells_ = other.continuity_.cells_;
  continuity_.density_ = other.continuity_.density_;
  continuity_.velocity_ = other.continuity_.velocity_;
  continuity_.pressure_reference_ = other.continuity_.pressure_reference_;
  continuity_.geometry_revision_ = other.continuity_.geometry_revision_;
  continuity_.fingerprint_ = other.continuity_.fingerprint_;
  momentum_.cells_ = other.momentum_.cells_;
  momentum_.density_ = other.momentum_.density_;
  momentum_.velocity_ = other.momentum_.velocity_;
  momentum_.pressure_ = other.momentum_.pressure_;
  momentum_.viscosity_ = other.momentum_.viscosity_;
  momentum_.velocity_gradient_ = other.momentum_.velocity_gradient_;
  momentum_.convection_ = other.momentum_.convection_;
  momentum_.convection_reach_ = other.momentum_.convection_reach_;
  momentum_.stress_form_ = other.momentum_.stress_form_;
  momentum_.contributions_ = std::move(other.momentum_.contributions_);
  momentum_.geometry_revision_ = other.momentum_.geometry_revision_;
  momentum_.boundary_revision_ = other.momentum_.boundary_revision_;
  momentum_.transport_fingerprint_ = other.momentum_.transport_fingerprint_;
  momentum_.fingerprint_ = other.momentum_.fingerprint_;
  enthalpy_.cells_ = other.enthalpy_.cells_;
  enthalpy_.density_ = other.enthalpy_.density_;
  enthalpy_.velocity_ = other.enthalpy_.velocity_;
  enthalpy_.pressure_ = other.enthalpy_.pressure_;
  enthalpy_.enthalpy_ = other.enthalpy_.enthalpy_;
  enthalpy_.temperature_ = other.enthalpy_.temperature_;
  enthalpy_.velocity_gradient_ = other.enthalpy_.velocity_gradient_;
  enthalpy_.convection_ = other.enthalpy_.convection_;
  enthalpy_.convection_reach_ = other.enthalpy_.convection_reach_;
  enthalpy_.contributions_ = std::move(other.enthalpy_.contributions_);
  enthalpy_.geometry_revision_ = other.enthalpy_.geometry_revision_;
  enthalpy_.boundary_revision_ = other.enthalpy_.boundary_revision_;
  enthalpy_.thermodynamics_fingerprint_ =
      other.enthalpy_.thermodynamics_fingerprint_;
  enthalpy_.transport_fingerprint_ = other.enthalpy_.transport_fingerprint_;
  enthalpy_.fingerprint_ = other.enthalpy_.fingerprint_;
  species_.cells_ = other.species_.cells_;
  species_.density_ = other.species_.density_;
  species_.convection_ = other.species_.convection_;
  species_.convection_reach_ = other.species_.convection_reach_;
  species_.specs_ = std::move(other.species_.specs_);
  species_.contribution_counts_ =
      std::move(other.species_.contribution_counts_);
  species_.contribution_begins_ =
      std::move(other.species_.contribution_begins_);
  species_.contributions_ = std::move(other.species_.contributions_);
  species_.geometry_revision_ = other.species_.geometry_revision_;
  species_.boundary_revision_ = other.species_.boundary_revision_;
  species_.thermodynamics_fingerprint_ =
      other.species_.thermodynamics_fingerprint_;
  species_.transport_fingerprint_ = other.species_.transport_fingerprint_;
  species_.fingerprint_ = other.species_.fingerprint_;
  scalars_.cells_ = other.scalars_.cells_;
  scalars_.density_ = other.scalars_.density_;
  scalars_.convection_ = other.scalars_.convection_;
  scalars_.convection_reach_ = other.scalars_.convection_reach_;
  scalars_.specs_ = std::move(other.scalars_.specs_);
  scalars_.contribution_counts_ =
      std::move(other.scalars_.contribution_counts_);
  scalars_.contribution_begins_ =
      std::move(other.scalars_.contribution_begins_);
  scalars_.contributions_ = std::move(other.scalars_.contributions_);
  scalars_.geometry_revision_ = other.scalars_.geometry_revision_;
  scalars_.boundary_revision_ = other.scalars_.boundary_revision_;
  scalars_.thermodynamics_fingerprint_ =
      other.scalars_.thermodynamics_fingerprint_;
  scalars_.transport_fingerprint_ = other.scalars_.transport_fingerprint_;
  scalars_.fingerprint_ = other.scalars_.fingerprint_;
  pressure_reference_ = other.pressure_reference_;
  thermophysical_predictor_.cells_ =
      other.thermophysical_predictor_.cells_;
  thermophysical_predictor_.patch_begin_ =
      other.thermophysical_predictor_.patch_begin_;
  thermophysical_predictor_.density_ =
      other.thermophysical_predictor_.density_;
  thermophysical_predictor_.enthalpy_ =
      other.thermophysical_predictor_.enthalpy_;
  thermophysical_predictor_.species_ =
      std::move(other.thermophysical_predictor_.species_);
  thermophysical_predictor_.passive_scalars_ =
      std::move(other.thermophysical_predictor_.passive_scalars_);
  thermophysical_predictor_.species_enthalpy_minimum_ =
      std::move(other.thermophysical_predictor_.species_enthalpy_minimum_);
  thermophysical_predictor_.species_enthalpy_maximum_ =
      std::move(other.thermophysical_predictor_.species_enthalpy_maximum_);
  thermophysical_predictor_.dependent_enthalpy_minimum_ =
      other.thermophysical_predictor_.dependent_enthalpy_minimum_;
  thermophysical_predictor_.dependent_enthalpy_maximum_ =
      other.thermophysical_predictor_.dependent_enthalpy_maximum_;
  thermophysical_predictor_.enthalpy_convection_ =
      other.thermophysical_predictor_.enthalpy_convection_;
  thermophysical_predictor_.species_convection_ =
      other.thermophysical_predictor_.species_convection_;
  thermophysical_predictor_.passive_scalar_convection_ =
      other.thermophysical_predictor_.passive_scalar_convection_;
  thermophysical_predictor_.enthalpy_reach_ =
      other.thermophysical_predictor_.enthalpy_reach_;
  thermophysical_predictor_.species_reach_ =
      other.thermophysical_predictor_.species_reach_;
  thermophysical_predictor_.passive_scalar_reach_ =
      other.thermophysical_predictor_.passive_scalar_reach_;
  thermophysical_predictor_.boundary_ =
      other.thermophysical_predictor_.boundary_;
  thermophysical_predictor_.geometry_revision_ =
      other.thermophysical_predictor_.geometry_revision_;
  thermophysical_predictor_.boundary_revision_ =
      other.thermophysical_predictor_.boundary_revision_;
  thermophysical_predictor_.transport_fingerprint_ =
      other.thermophysical_predictor_.transport_fingerprint_;
  thermophysical_predictor_.fingerprint_ =
      other.thermophysical_predictor_.fingerprint_;
  global_cells_ = other.global_cells_;
  patch_begin_ = other.patch_begin_;
  local_cells_ = other.local_cells_;
  semantic_fingerprint_ = other.semantic_fingerprint_;
  thermodynamics_fingerprint_ = other.thermodynamics_fingerprint_;
  transport_fingerprint_ = other.transport_fingerprint_;
  fingerprint_ = other.fingerprint_;
  rebind();
  other.reset();
}

EquationPlanSet::EquationPlanSet(EquationPlanSet&& other) noexcept {
  move_from(std::move(other));
}

EquationPlanSet& EquationPlanSet::operator=(EquationPlanSet&& other) noexcept {
  if (this != &other) {
    reset();
    move_from(std::move(other));
  }
  return *this;
}

Status EquationPlanSet::compile(
    MPI_Comm communicator, const SchemePlan& schemes,
    const CartesianGeometryPlan& geometry, const MeshPatch& patch,
    const BoundaryPlan& boundary, const ContributionRegistry& contributions,
    const ThermodynamicsPlan& thermodynamics, const TransportPlan& transport,
    const EquationPlanSpec& spec, EquationPlanSet& out,
    EquationCompileDiagnostics* diagnostics) noexcept {
  if (diagnostics != nullptr) {
    diagnostics->lowest_failing_rank = -1;
  }
  int rank = 0;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kEquationCollective};
  }
  std::size_t local_cells = 0U;
  Status local = local_plan_status(schemes, geometry, patch, boundary,
                                   contributions, thermodynamics, transport,
                                   spec, local_cells);
  int lowest = -1;
  Status consensus =
      collective_status(communicator, local, rank, size, lowest);
  if (!consensus) {
    if (diagnostics != nullptr) {
      diagnostics->lowest_failing_rank = lowest;
    }
    return consensus;
  }

  const PlanFingerprint semantic = compute_semantic_fingerprint(
      schemes, geometry, boundary, contributions, thermodynamics, transport,
      spec);
  PlanFingerprint semantic_min = semantic;
  PlanFingerprint semantic_max = semantic;
  if (MPI_Allreduce(MPI_IN_PLACE, &semantic_min, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kEquationCollective};
  }
  if (MPI_Allreduce(MPI_IN_PLACE, &semantic_max, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kEquationCollective};
  }
  if (semantic_min != semantic_max) {
    PlanFingerprint authority = semantic;
    if (MPI_Bcast(&authority, 1, MPI_UINT64_T, 0, communicator) !=
        MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kEquationCollective};
    }
    const int divergent = semantic == authority ? size : rank;
    if (MPI_Allreduce(&divergent, &lowest, 1, MPI_INT, MPI_MIN,
                      communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kEquationCollective};
    }
    if (diagnostics != nullptr) {
      diagnostics->lowest_failing_rank = lowest < size ? lowest : -1;
    }
    return {StatusCode::invalid_plan, kEquationCollective};
  }

  EquationPlanSet candidate;
  try {
    local = CartesianKernelPlan::compile(schemes, geometry, patch, boundary,
                                         candidate.kernels_);
    if (!local) {
      throw local;
    }
    candidate.global_cells_ = geometry.global_cells();
    candidate.patch_begin_ = patch.begin;
    candidate.local_cells_ = local_cells;
    candidate.semantic_fingerprint_ = semantic;
    candidate.thermodynamics_fingerprint_ = thermodynamics.fingerprint();
    candidate.transport_fingerprint_ = transport.fingerprint();
    std::uint64_t layout_hash = hash_mix(kFnvOffset, semantic);
    layout_hash = hash_mix(layout_hash, boundary.local_layout_fingerprint());
    layout_hash = hash_mix(layout_hash, candidate.kernels_.fingerprint());
    layout_hash = hash_mix(layout_hash, static_cast<std::uint32_t>(patch.begin.x));
    layout_hash = hash_mix(layout_hash, static_cast<std::uint32_t>(patch.begin.y));
    layout_hash = hash_mix(layout_hash, static_cast<std::uint32_t>(patch.begin.z));
    layout_hash = hash_mix(layout_hash, local_cells);
    candidate.fingerprint_ = finish_hash(layout_hash);

    candidate.continuity_.cells_ = patch.cells;
    candidate.continuity_.density_ = spec.density;
    candidate.continuity_.velocity_ = spec.velocity;
    candidate.continuity_.pressure_reference_ = spec.pressure_reference;
    candidate.continuity_.geometry_revision_ = geometry.topology_revision();
    candidate.continuity_.fingerprint_ =
        child_fingerprint(semantic, UINT64_C(0x636f6e74696e), spec.density);

    candidate.momentum_.cells_ = patch.cells;
    candidate.momentum_.density_ = spec.density;
    candidate.momentum_.velocity_ = spec.velocity;
    candidate.momentum_.pressure_ = spec.pressure_perturbation;
    candidate.momentum_.viscosity_ = spec.effective_viscosity;
    candidate.momentum_.velocity_gradient_ = spec.velocity_gradient;
    candidate.momentum_.convection_ = schemes.momentum();
    candidate.momentum_.convection_reach_ =
        convection_reach(schemes.momentum());
    candidate.momentum_.geometry_revision_ = geometry.topology_revision();
    candidate.momentum_.boundary_revision_ = boundary.revision();
    candidate.momentum_.transport_fingerprint_ = transport.fingerprint();
    for (std::size_t i = 0U; i < contributions.contributions().size; ++i) {
      if (contributions.contributions().data[i].conserved_quantity ==
          spec.velocity) {
        candidate.momentum_.contributions_.push_back(
            contributions.contributions().data[i]);
      }
    }
    candidate.momentum_.fingerprint_ =
        child_fingerprint(semantic, UINT64_C(0x6d6f6d656e74), spec.velocity);

    candidate.enthalpy_.cells_ = patch.cells;
    candidate.enthalpy_.density_ = spec.density;
    candidate.enthalpy_.velocity_ = spec.velocity;
    candidate.enthalpy_.pressure_ = spec.pressure_perturbation;
    candidate.enthalpy_.enthalpy_ = spec.enthalpy;
    candidate.enthalpy_.temperature_ = spec.temperature;
    candidate.enthalpy_.velocity_gradient_ = spec.velocity_gradient;
    candidate.enthalpy_.convection_ = schemes.enthalpy();
    candidate.enthalpy_.convection_reach_ =
        convection_reach(schemes.enthalpy());
    candidate.enthalpy_.geometry_revision_ = geometry.topology_revision();
    candidate.enthalpy_.boundary_revision_ = boundary.revision();
    candidate.enthalpy_.thermodynamics_fingerprint_ =
        thermodynamics.fingerprint();
    candidate.enthalpy_.transport_fingerprint_ = transport.fingerprint();
    for (std::size_t i = 0U; i < contributions.contributions().size; ++i) {
      if (contributions.contributions().data[i].conserved_quantity ==
          spec.enthalpy) {
        candidate.enthalpy_.contributions_.push_back(
            contributions.contributions().data[i]);
      }
    }
    candidate.enthalpy_.fingerprint_ =
        child_fingerprint(semantic, UINT64_C(0x656e7468616c), spec.enthalpy);

    for (std::size_t index = 0U; index < spec.scalars.size; ++index) {
      const ScalarEquationSpec scalar = spec.scalars.data[index];
      auto* plan = scalar.role == TransportedScalarRole::species
                       ? &candidate.species_.specs_
                       : &candidate.scalars_.specs_;
      auto* counts = scalar.role == TransportedScalarRole::species
                         ? &candidate.species_.contribution_counts_
                         : &candidate.scalars_.contribution_counts_;
      plan->push_back(scalar);
      auto* begins = scalar.role == TransportedScalarRole::species
                         ? &candidate.species_.contribution_begins_
                         : &candidate.scalars_.contribution_begins_;
      begins->push_back(0U);
      counts->push_back(0U);
    }
    candidate.species_.cells_ = patch.cells;
    candidate.species_.density_ = spec.density;
    candidate.species_.convection_ = schemes.species();
    candidate.species_.convection_reach_ =
        convection_reach(schemes.species());
    candidate.species_.geometry_revision_ = geometry.topology_revision();
    candidate.species_.boundary_revision_ = boundary.revision();
    candidate.species_.thermodynamics_fingerprint_ =
        thermodynamics.fingerprint();
    candidate.species_.transport_fingerprint_ = transport.fingerprint();
    candidate.species_.fingerprint_ =
        child_fingerprint(semantic, UINT64_C(0x737065636965));
    candidate.scalars_.cells_ = patch.cells;
    candidate.scalars_.density_ = spec.density;
    candidate.scalars_.convection_ = schemes.passive_scalar();
    candidate.scalars_.convection_reach_ =
        convection_reach(schemes.passive_scalar());
    candidate.scalars_.geometry_revision_ = geometry.topology_revision();
    candidate.scalars_.boundary_revision_ = boundary.revision();
    candidate.scalars_.thermodynamics_fingerprint_ =
        thermodynamics.fingerprint();
    candidate.scalars_.transport_fingerprint_ = transport.fingerprint();
    candidate.scalars_.fingerprint_ =
        child_fingerprint(semantic, UINT64_C(0x7363616c6172));
    candidate.pressure_reference_.cells_ = patch.cells;
    candidate.pressure_reference_.kind_ = spec.pressure_reference;
    candidate.pressure_reference_.gauge_ =
        spec.pressure_reference == PressureReferenceKind::closed_mass
            ? PressureGauge::compressibility_weighted_zero_mean
            : PressureGauge::absolute_boundary_dirichlet;
    candidate.pressure_reference_.pressure_perturbation_field_ =
        spec.pressure_perturbation;
    candidate.pressure_reference_.compressibility_field_ =
        spec.pressure_compressibility;
    candidate.pressure_reference_.service_stage_ =
        spec.closed_mass_service_stage;
    candidate.pressure_reference_.geometry_revision_ =
        geometry.topology_revision();
    candidate.pressure_reference_.thermodynamics_fingerprint_ =
        thermodynamics.fingerprint();
    candidate.pressure_reference_.fingerprint_ = child_fingerprint(
        semantic, UINT64_C(0x7072657373757265),
        spec.pressure_compressibility);
    candidate.thermophysical_predictor_.cells_ = patch.cells;
    candidate.thermophysical_predictor_.patch_begin_ = patch.begin;
    candidate.thermophysical_predictor_.density_ = spec.density;
    candidate.thermophysical_predictor_.enthalpy_ = spec.enthalpy;
    candidate.thermophysical_predictor_.enthalpy_convection_ =
        schemes.enthalpy();
    candidate.thermophysical_predictor_.species_convection_ =
        schemes.species();
    candidate.thermophysical_predictor_.passive_scalar_convection_ =
        schemes.passive_scalar();
    candidate.thermophysical_predictor_.enthalpy_reach_ =
        convection_reach(schemes.enthalpy());
    candidate.thermophysical_predictor_.species_reach_ =
        convection_reach(schemes.species());
    candidate.thermophysical_predictor_.passive_scalar_reach_ =
        convection_reach(schemes.passive_scalar());
    candidate.thermophysical_predictor_.boundary_ = &boundary;
    candidate.thermophysical_predictor_.geometry_revision_ =
        geometry.topology_revision();
    candidate.thermophysical_predictor_.boundary_revision_ =
        boundary.revision();
    candidate.thermophysical_predictor_.transport_fingerprint_ =
        transport.fingerprint();
    for (const ScalarEquationSpec& scalar :
         candidate.species_.specs_) {
      candidate.thermophysical_predictor_.species_.push_back(scalar.field);
    }
    candidate.thermophysical_predictor_.species_enthalpy_minimum_.resize(
        candidate.thermophysical_predictor_.species_.size());
    candidate.thermophysical_predictor_.species_enthalpy_maximum_.resize(
        candidate.thermophysical_predictor_.species_.size());
    for (std::size_t index = 0U;
         index < candidate.thermophysical_predictor_.species_.size();
         ++index) {
      local = thermodynamics.independent_species_enthalpy_bounds(
          index,
          candidate.thermophysical_predictor_
              .species_enthalpy_minimum_[index],
          candidate.thermophysical_predictor_
              .species_enthalpy_maximum_[index]);
      if (!local) throw local;
    }
    local = thermodynamics.species_enthalpy_bounds(
        thermodynamics.dependent_species_index(),
        candidate.thermophysical_predictor_.dependent_enthalpy_minimum_,
        candidate.thermophysical_predictor_.dependent_enthalpy_maximum_);
    if (!local) throw local;
    for (const ScalarEquationSpec& scalar : candidate.scalars_.specs_) {
      candidate.thermophysical_predictor_.passive_scalars_.push_back(
          scalar.field);
    }
    candidate.thermophysical_predictor_.fingerprint_ = child_fingerprint(
        semantic, UINT64_C(0x746865726d707265), spec.enthalpy);
    candidate.pressure_reference_.predictor_fingerprint_ =
        candidate.thermophysical_predictor_.fingerprint_;
    std::size_t begin = 0U;
    for (std::size_t i = 0U; i < candidate.species_.specs_.size(); ++i) {
      candidate.species_.contribution_begins_[i] = begin;
      for (std::size_t contribution = 0U;
           contribution < contributions.contributions().size;
           ++contribution) {
        const CompiledContribution descriptor =
            contributions.contributions().data[contribution];
        if (descriptor.conserved_quantity ==
            candidate.species_.specs_[i].field) {
          candidate.species_.contributions_.push_back(descriptor);
          ++candidate.species_.contribution_counts_[i];
          ++begin;
        }
      }
    }
    begin = 0U;
    for (std::size_t i = 0U; i < candidate.scalars_.specs_.size(); ++i) {
      candidate.scalars_.contribution_begins_[i] = begin;
      for (std::size_t contribution = 0U;
           contribution < contributions.contributions().size;
           ++contribution) {
        const CompiledContribution descriptor =
            contributions.contributions().data[contribution];
        if (descriptor.conserved_quantity ==
            candidate.scalars_.specs_[i].field) {
          candidate.scalars_.contributions_.push_back(descriptor);
          ++candidate.scalars_.contribution_counts_[i];
          ++begin;
        }
      }
    }
    candidate.rebind();
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kEquationPlan};
  } catch (Status status) {
    local = status;
  } catch (...) {
    local = {StatusCode::invalid_plan, kEquationPlan};
  }
  consensus = collective_status(communicator, local, rank, size, lowest);
  if (!consensus) {
    if (diagnostics != nullptr) {
      diagnostics->lowest_failing_rank = lowest;
    }
    return consensus;
  }
  out = std::move(candidate);
  if (diagnostics != nullptr) {
    diagnostics->lowest_failing_rank = -1;
  }
  return {};
}

Status assemble_continuity_impl(
    const ContinuityEquationPlan& plan, const EquationStateView& state,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate, bool allow_partial) noexcept {
  if (plan.kernels_ == nullptr || plan.fingerprint_ == 0U ||
      plan.kernels_->fingerprint() == 0U ||
      context.geometry != plan.geometry_revision_ ||
      !compatible_history(state.density, plan.cells_, plan.density_, 1U, 0U)) {
    return {StatusCode::invalid_plan, kEquationAssembly};
  }
  KernelBox box{};
  if (!valid_context(*plan.kernels_, context, box) ||
      (!allow_partial && !detail::full_equation_box(box, plan.cells_)) ||
      !detail::valid_cell_view(system.residual, plan.cells_, 0U, 1U) ||
      !detail::finite_face_flux(context.mass_flux, box) ||
      detail::output_aliases_flux(system, false, context.mass_flux) ||
      detail::field_views_overlap(state.density.trial,
                                  as_const(system.residual)) ||
      detail::field_views_overlap(state.density.accepted,
                                  as_const(system.residual)) ||
      detail::field_views_overlap(state.density.previous,
                                  as_const(system.residual))) {
    return {StatusCode::invalid_plan, kEquationAssembly};
  }

  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double trial = state.density.trial.unchecked(cell, 0U);
        const double accepted = state.density.accepted.unchecked(cell, 0U);
        const double previous = state.density.previous.unchecked(cell, 0U);
        const double volume = detail::cell_volume(*plan.kernels_, cell);
        const double value =
            (context.bdf.a0 * trial + context.bdf.a1 * accepted +
             context.bdf.a2 * previous) * volume;
        if (!std::isfinite(trial) || !std::isfinite(accepted) ||
            !std::isfinite(previous) || trial <= 0.0 || accepted <= 0.0 ||
            previous <= 0.0 || !std::isfinite(value)) {
          return {StatusCode::numerical_failure, kEquationNumerical};
        }
      }
    }
  }

  const std::array<FieldView, 1U> writes{system.residual};
  const KernelInvocation divergence{{}, {writes.data(), writes.size()}, box,
                                    0U, 0U, 1U, context.face_flux,
                                    context.counters};
  const Status flux_status =
      context.scope == EquationAssemblyScope::final_conservative
          ? cartesian_face_divergence(*plan.kernels_, context.mass_flux,
                                      divergence)
          : cartesian_provisional_face_divergence(
                *plan.kernels_, context.mass_flux, divergence);
  if (!flux_status) {
    return flux_status;
  }

  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double unsteady =
            context.bdf.a0 * state.density.trial.unchecked(cell, 0U) +
            context.bdf.a1 * state.density.accepted.unchecked(cell, 0U) +
            context.bdf.a2 * state.density.previous.unchecked(cell, 0U);
        const double volume = detail::cell_volume(*plan.kernels_, cell);
        const double value =
            (system.residual.unchecked(cell, 0U) + unsteady) * volume;
        if (!std::isfinite(value)) {
          return {StatusCode::numerical_failure, kEquationNumerical};
        }
        system.residual.unchecked(cell, 0U) = value;
      }
    }
  }
  certificate = make_certificate(plan.fingerprint_, state, context);
  return {};
}

Status assemble_tile(AssemblyEpoch& epoch,
                     const ContinuityEquationPlan& plan,
                     const EquationStateView& state, KernelBox box) noexcept {
  if (!epoch.active_) {
    return {StatusCode::invalid_plan, kAssemblyEpoch};
  }
  EquationAssemblyContext context = epoch.context_;
  context.box = box;
  EquationAssemblyCertificate candidate;
  const Status status = assemble_continuity_impl(
      plan, state, context, epoch.system_, candidate, true);
  if (!status) {
    return epoch.fail(status);
  }
  return epoch.record(box, plan.cells(), candidate);
}

Status assemble_continuity(
    const ContinuityEquationPlan& plan, const EquationStateView& state,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate) noexcept {
  const KernelBox box = resolved_box(context, plan.cells_);
  if (!detail::full_equation_box(box, plan.cells_)) {
    return {StatusCode::invalid_plan, kEquationAssembly};
  }
  EquationAssemblyContext epoch_context = context;
  epoch_context.box = {};
  AssemblyEpoch epoch;
  Status status = epoch.begin(epoch_context, system);
  if (status) {
    status = assemble_tile(epoch, plan, state, box);
  }
  if (!status) {
    return status;
  }
  return epoch.finalize(certificate);
}

Status assemble_pressure_storage(
    const PressureReferencePlan& plan, ConstFieldView drho_dp_hY,
    const PressureReferenceCertificate& pressure_reference,
    const EquationAssemblyContext& context, FieldView cell_integral_diagonal,
    EquationAssemblyCertificate& certificate) noexcept {
  if (plan.fingerprint_ == 0U || plan.kernels_ == nullptr ||
      !pressure_reference.valid() ||
      pressure_reference.plan != plan.fingerprint_ ||
      pressure_reference.predictor != plan.predictor_fingerprint_ ||
      pressure_reference.thermodynamics !=
          plan.thermodynamics_fingerprint_ ||
      pressure_reference.time != context.time ||
      pressure_reference.kind != plan.kind_ ||
      context.geometry != plan.geometry_revision_ ||
      context.thermo != plan.thermodynamics_fingerprint_ ||
      (plan.kind_ == PressureReferenceKind::closed_mass
           ? context.contribution_stage != plan.service_stage_
           : plan.service_stage_ != 0U) ||
      context.time == 0U ||
      !detail::bdf_matches_time_step(context.bdf, context.dt) ||
      drho_dp_hY.field != plan.compressibility_field_ ||
      drho_dp_hY.revision != context.thermo ||
      !detail::valid_cell_view(drho_dp_hY, plan.cells_, 0U, 1U, 0U) ||
      !detail::valid_cell_view(cell_integral_diagonal, plan.cells_, 0U, 1U) ||
      detail::component_ranges_overlap(drho_dp_hY, 0U, 1U,
                                       cell_integral_diagonal, 0U, 1U)) {
    return {StatusCode::invalid_plan, kEquationAssembly};
  }
  const KernelBox box = resolved_box(context, plan.cells_);
  if (!detail::valid_kernel_box(box, plan.cells_) ||
      !detail::full_equation_box(box, plan.cells_)) {
    return {StatusCode::invalid_plan, kEquationAssembly};
  }
  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  // Validate first so a bad EOS derivative leaves the caller's numeric state
  // unchanged.
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double psi = drho_dp_hY.unchecked(cell, 0U);
        const double volume = detail::cell_volume(*plan.kernels_, cell);
        const double diagonal = context.bdf.a0 * psi * volume;
        if (!std::isfinite(psi) || psi <= 0.0 ||
            !std::isfinite(volume) || volume <= 0.0 ||
            !std::isfinite(diagonal) || diagonal <= 0.0) {
          return {StatusCode::numerical_failure, kEquationNumerical};
        }
      }
    }
  }
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        cell_integral_diagonal.unchecked(cell, 0U) =
            context.bdf.a0 * drho_dp_hY.unchecked(cell, 0U) *
            detail::cell_volume(*plan.kernels_, cell);
      }
    }
  }
  std::uint64_t state = hash_mix(kFnvOffset, drho_dp_hY.revision);
  state = hash_mix(state, drho_dp_hY.storage_identity);
  state = hash_mix(state, drho_dp_hY.revision_domain);
  state = hash_mix(state, pressure_reference.pressure_reference);
  certificate = {plan.fingerprint_, context.scope, context.time,
                 context.geometry, context.face_flux, finish_hash(state),
                 context.dt};
  return {};
}

}  // namespace hundun::v04
