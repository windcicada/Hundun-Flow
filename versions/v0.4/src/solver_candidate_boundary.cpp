// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <limits>
#include <new>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kCandidateBoundaryFinalizer = 1711U;
constexpr std::uint32_t kCandidateBoundaryNumerical = 1712U;
constexpr std::uint32_t kCandidateBoundaryBackflow = 1713U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  return hash * kFnvPrime;
}

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool finite_positive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

bool same_view_identity(ConstFieldView left, ConstFieldView right) noexcept {
  return left.base == right.base && left.replica == right.replica &&
         same(left.interior, right.interior) && same(left.ghosts, right.ghosts) &&
         left.components == right.components &&
         left.stride_y == right.stride_y && left.stride_z == right.stride_z &&
         left.component_stride == right.component_stride &&
         left.field == right.field && left.revision == right.revision &&
         left.storage_identity == right.storage_identity &&
         left.revision_domain == right.revision_domain;
}

// Candidate p/h live in independent scratch FieldIds.  The boundary closure
// consumes semantic aliases carrying the BoundaryPlan p/h FieldIds, but the
// alias must still name the exact same immutable payload and runtime lineage.
bool same_view_payload_identity(ConstFieldView left,
                                ConstFieldView right) noexcept {
  return left.base == right.base && left.replica == right.replica &&
         same(left.interior, right.interior) && same(left.ghosts, right.ghosts) &&
         left.components == right.components &&
         left.stride_y == right.stride_y && left.stride_z == right.stride_z &&
         left.component_stride == right.component_stride &&
         left.revision == right.revision &&
         left.storage_identity == right.storage_identity &&
         left.revision_domain == right.revision_domain;
}

bool same_face_identity(ConstFaceFieldView left,
                        ConstFaceFieldView right) noexcept {
  return left.base == right.base && same(left.extents, right.extents) &&
         left.stride_y == right.stride_y &&
         left.stride_z == right.stride_z && left.axis == right.axis &&
         left.storage_identity == right.storage_identity &&
         left.revision_domain == right.revision_domain;
}

bool same_flux_identity(ConstFaceFluxView left,
                        ConstFaceFluxView right) noexcept {
  return left.revision == right.revision &&
         same_face_identity(left.x, right.x) &&
         same_face_identity(left.y, right.y) &&
         same_face_identity(left.z, right.z);
}

std::uint64_t mix_view_binding(std::uint64_t hash,
                               ConstFieldView view) noexcept {
  hash = mix(hash, reinterpret_cast<std::uintptr_t>(view.base));
  hash = mix(hash, view.replica);
  hash = mix(hash, view.field);
  hash = mix(hash, view.revision);
  hash = mix(hash, view.storage_identity);
  return mix(hash, view.revision_domain);
}

std::uint64_t mix_field_values(std::uint64_t hash, ConstFieldView field,
                               Int3 cells, std::uint8_t components) noexcept {
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        hash = mix(hash, static_cast<std::uint32_t>(x));
        hash = mix(hash, static_cast<std::uint32_t>(y));
        hash = mix(hash, static_cast<std::uint32_t>(z));
        for (std::uint8_t component = 0U; component < components;
             ++component) {
          hash = mix(hash,
                     double_bits(field.unchecked({x, y, z}, component)));
        }
      }
  return hash;
}

std::uint64_t mix_stage_field_values(std::uint64_t hash,
                                     ConstFieldView field, Int3 cells,
                                     std::uint8_t components,
                                     std::int32_t ghost_reach) noexcept {
  for (std::uint8_t component_index = 0U;
       component_index < components; ++component_index)
    for (std::int32_t z = -ghost_reach;
         z < cells.z + ghost_reach; ++z)
      for (std::int32_t y = -ghost_reach;
           y < cells.y + ghost_reach; ++y)
        for (std::int32_t x = -ghost_reach;
             x < cells.x + ghost_reach; ++x) {
          const unsigned outside =
              static_cast<unsigned>(x < 0 || x >= cells.x) +
              static_cast<unsigned>(y < 0 || y >= cells.y) +
              static_cast<unsigned>(z < 0 || z >= cells.z);
          if (outside > 1U) continue;
          hash = mix(hash,
                     double_bits(field.unchecked({x, y, z},
                                                  component_index)));
        }
  return hash;
}

std::uint64_t state_local_provenance(
    const PressureEnergyCandidateBoundaryFinalizeInput& input, Int3 cells,
    PlanFingerprint thermodynamics,
    PlanFingerprint composition_numeric_provenance) noexcept {
  std::uint64_t hash = mix(kFnvOffset, UINT64_C(0x7065636266737461));
  hash = mix(hash, thermodynamics);
  hash = mix(hash, input.pressure_reference.pressure_reference);
  hash = mix(hash, input.composition_identity);
  hash = mix(hash, composition_numeric_provenance);
  hash = mix(hash, double_bits(input.absolute_pressure_reference));
  hash = mix_field_values(hash, input.pressure_perturbation, cells, 1U);
  hash = mix_field_values(hash, input.enthalpy, cells, 1U);
  hash = mix_field_values(hash, input.density, cells, 1U);
  hash = mix_field_values(hash, input.temperature, cells, 1U);
  return mix_field_values(hash, input.velocity, cells, 3U);
}

std::uint64_t composition_numeric_local_provenance(
    Span<const ConstFieldView> independent_species, Int3 cells) noexcept {
  std::uint64_t hash = mix(kFnvOffset, UINT64_C(0x636f6d706e756d76));
  hash = mix(hash, static_cast<std::uint64_t>(independent_species.size));
  for (std::size_t species = 0U; species < independent_species.size;
       ++species) {
    const ConstFieldView field = independent_species.data[species];
    hash = mix(hash, field.components);
    for (std::uint8_t component = 0U; component < field.components;
         ++component)
      for (std::int32_t z = 0; z < cells.z; ++z)
        for (std::int32_t y = 0; y < cells.y; ++y)
          for (std::int32_t x = 0; x < cells.x; ++x)
            hash = mix(
                hash,
                double_bits(field.unchecked({x, y, z}, component)));
  }
  return hash == 0U ? 1U : hash;
}

std::uint64_t mix_face_values(std::uint64_t hash,
                              ConstFaceFieldView face) noexcept {
  hash = mix(hash, static_cast<std::uint8_t>(face.axis));
  for (std::int32_t z = 0; z < face.extents.z; ++z)
    for (std::int32_t y = 0; y < face.extents.y; ++y)
      for (std::int32_t x = 0; x < face.extents.x; ++x)
        hash = mix(hash, double_bits(face.unchecked({x, y, z})));
  return hash;
}

std::uint64_t flux_local_provenance(ConstFaceFluxView flux) noexcept {
  std::uint64_t hash = mix(kFnvOffset, UINT64_C(0x7065636266666c78));
  hash = mix_face_values(hash, flux.x);
  hash = mix_face_values(hash, flux.y);
  return mix_face_values(hash, flux.z);
}

std::uint64_t mix_stage_flux_values(std::uint64_t hash,
                                    ConstFaceFluxView flux) noexcept {
  const std::array<ConstFaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (ConstFaceFieldView face : faces)
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          hash = mix(hash, double_bits(face.unchecked({x, y, z})));
  return hash;
}

CartesianFace selected_face(CartesianAxis axis, bool high) noexcept {
  return static_cast<CartesianFace>(
      2U * static_cast<std::size_t>(axis) + (high ? 1U : 0U));
}

template <class Visit>
void for_each_boundary_face_storage_order(CartesianAxis axis, Int3 cells,
                                          Visit&& visit) noexcept {
  if (axis == CartesianAxis::x) {
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::uint8_t side = 0U; side < 2U; ++side) {
          const bool high = side != 0U;
          if (!visit({high ? cells.x : 0, y, z}, high)) return;
        }
    return;
  }
  if (axis == CartesianAxis::y) {
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::uint8_t side = 0U; side < 2U; ++side) {
        const bool high = side != 0U;
        for (std::int32_t x = 0; x < cells.x; ++x)
          if (!visit({x, high ? cells.y : 0, z}, high)) return;
      }
    return;
  }
  for (std::uint8_t side = 0U; side < 2U; ++side) {
    const bool high = side != 0U;
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        if (!visit({x, y, high ? cells.z : 0}, high)) return;
  }
}

double outward_sign(bool high) noexcept { return high ? 1.0 : -1.0; }

Int3 owner_cell(CartesianAxis axis, bool high, Int3 face,
                Int3 cells) noexcept {
  Int3 owner = face;
  if (axis == CartesianAxis::x)
    owner.x = high ? cells.x - 1 : 0;
  else if (axis == CartesianAxis::y)
    owner.y = high ? cells.y - 1 : 0;
  else
    owner.z = high ? cells.z - 1 : 0;
  return owner;
}

double area(const CartesianGeometryPlan& geometry, MeshPatch patch,
            CartesianAxis axis, Int3 face) noexcept {
  const std::int32_t gx = patch.begin.x + face.x;
  const std::int32_t gy = patch.begin.y + face.y;
  const std::int32_t gz = patch.begin.z + face.z;
  const Span<const double> dx = geometry.x().widths();
  const Span<const double> dy = geometry.y().widths();
  const Span<const double> dz = geometry.z().widths();
  if (axis == CartesianAxis::x) {
    if (gy < 0 || gz < 0 || static_cast<std::size_t>(gy) >= dy.size ||
        static_cast<std::size_t>(gz) >= dz.size)
      return 0.0;
    return dy.data[gy] * dz.data[gz];
  }
  if (axis == CartesianAxis::y) {
    if (gx < 0 || gz < 0 || static_cast<std::size_t>(gx) >= dx.size ||
        static_cast<std::size_t>(gz) >= dz.size)
      return 0.0;
    return dx.data[gx] * dz.data[gz];
  }
  if (gx < 0 || gy < 0 || static_cast<std::size_t>(gx) >= dx.size ||
      static_cast<std::size_t>(gy) >= dy.size)
    return 0.0;
  return dx.data[gx] * dy.data[gy];
}

Real3 vector_parameter(const BoundaryPlan& boundary, std::uint32_t parameter,
                       bool backflow, bool direction) noexcept {
  if (direction)
    return {boundary.direction_x().data[parameter],
            boundary.direction_y().data[parameter],
            boundary.direction_z().data[parameter]};
  if (backflow)
    return {boundary.backflow_velocity_x().data[parameter],
            boundary.backflow_velocity_y().data[parameter],
            boundary.backflow_velocity_z().data[parameter]};
  return {boundary.velocity_x().data[parameter],
          boundary.velocity_y().data[parameter],
          boundary.velocity_z().data[parameter]};
}

double component(Real3 value, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? value.x
             : (axis == CartesianAxis::y ? value.y : value.z);
}

bool parameter_storage_valid(const BoundaryPlan& boundary) noexcept {
  const std::size_t count = boundary.parameter_count();
  return boundary.velocity_x().size == count &&
         boundary.velocity_y().size == count &&
         boundary.velocity_z().size == count &&
         boundary.direction_x().size == count &&
         boundary.direction_y().size == count &&
         boundary.direction_z().size == count &&
         boundary.backflow_velocity_x().size == count &&
         boundary.backflow_velocity_y().size == count &&
         boundary.backflow_velocity_z().size == count &&
         boundary.scalar_targets().size == count &&
         boundary.scalar_backflow_targets().size == count &&
         boundary.pressure_targets().size == count &&
         boundary.temperature_targets().size == count &&
         boundary.backflow_temperature_targets().size == count &&
         boundary.mass_flow_targets().size == count &&
         boundary.allow_backflow().size == count;
}

bool scalar_parameter(const BoundaryPlan& boundary, CartesianFace face,
                      FieldId field, std::uint32_t& parameter) noexcept {
  const Span<const BoundaryIndexSpan> spans = boundary.spans();
  for (std::size_t index = 0U; index < spans.size; ++index) {
    const BoundaryIndexSpan& span = spans.data[index];
    if (span.stage == BoundaryStage::scalar && span.face == face &&
        span.field == field) {
      parameter = span.parameter;
      return parameter < boundary.parameter_count();
    }
  }
  return false;
}

bool supported(const BoundaryPlan& boundary) noexcept {
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

bool checked_face_count(Int3 cells, std::size_t& x, std::size_t& y,
                        std::size_t& z) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0) return false;
  const auto nx = static_cast<std::size_t>(cells.x);
  const auto ny = static_cast<std::size_t>(cells.y);
  const auto nz = static_cast<std::size_t>(cells.z);
  if (nx == std::numeric_limits<std::size_t>::max()) return false;
  const auto multiply = [](std::size_t left, std::size_t right,
                           std::size_t& value) noexcept {
    if (left != 0U &&
        right > std::numeric_limits<std::size_t>::max() / left) {
      return false;
    }
    value = left * right;
    return true;
  };
  std::size_t xy = 0U;
  return multiply(nx + 1U, ny, xy) && multiply(xy, nz, x) &&
         multiply(nx, ny + 1U, xy) && multiply(xy, nz, y) &&
         multiply(nx, ny, xy) && multiply(xy, nz + 1U, z);
}

}  // namespace

struct PressureEnergyCandidateBoundaryFinalizer::Impl {
  MPI_Comm communicator{MPI_COMM_NULL};
  int rank{};
  int size{};
  const CartesianGeometryPlan* geometry{};
  MeshPatch patch{};
  const BoundaryPlan* boundary{};
  const CartesianKernelPlan* kernels{};
  const ThermodynamicsPlan* thermodynamics{};
  const TransportPlan* transport{};
  const PressureVelocityCoupler* coupler{};
  const IbmEquationInterfacePlan* immersed_interface{};
  const IbmPhysicalBoundaryFluxAuthority* immersed_physical_boundary_flux{};
  RemoteDonorExchangePlan* candidate_pressure_correction_donors{};
  StageId candidate_pressure_correction_donor_stage{};
  FieldId candidate_pressure_correction_field{};
  std::uint8_t candidate_pressure_correction_donor_reach{};
  std::vector<double> x_scratch;
  std::vector<double> y_scratch;
  std::vector<double> z_scratch;
  std::vector<std::uint64_t> collective_hashes;
  std::vector<HaloFieldSpec> halo_contract;
  std::vector<double> composition;
  PlanFingerprint fingerprint{};
};

namespace {

struct PhysicalBoundaryFluxEngineInput {
  double absolute_pressure_reference{};
  ConstFieldView pressure_perturbation{};
  ConstFieldView velocity{};
  Span<const ConstFieldView> independent_species{};
  ConstFaceFluxView mechanical_flux{};
  RevisionToken final_flux_revision{};
};

struct PhysicalBoundaryFluxEngineResult {
  FaceFluxView prepared{};
  std::uint64_t local_inlet_hash{};
  std::uint64_t local_outlet_hash{};
  std::array<double, 8U> global_mass_target{};
  std::array<double, 6U> global_achieved{};
};

template <class Implementation>
Status prepare_physical_boundary_flux(
    Implementation& impl, const PhysicalBoundaryFluxEngineInput& input,
    ReductionEngine& reductions,
    PhysicalBoundaryFluxEngineResult& result) noexcept {
  result = {};
  const Int3 cells = impl.patch.cells;
  const auto make_scratch_face = [](std::vector<double>& storage,
                                    ConstFaceFieldView model) noexcept {
    const std::size_t stride_y =
        static_cast<std::size_t>(model.extents.x);
    const std::size_t stride_z =
        stride_y * static_cast<std::size_t>(model.extents.y);
    return FaceFieldView{storage.data(),
                         model.extents,
                         stride_y,
                         stride_z,
                         model.axis,
                         model.storage_identity,
                         model.revision_domain};
  };
  FaceFluxView prepared{
      make_scratch_face(impl.x_scratch, input.mechanical_flux.x),
      make_scratch_face(impl.y_scratch, input.mechanical_flux.y),
      make_scratch_face(impl.z_scratch, input.mechanical_flux.z),
      input.final_flux_revision,
      {}};
  const std::array<ConstFaceFieldView, 3U> mechanical_faces{
      input.mechanical_flux.x, input.mechanical_flux.y,
      input.mechanical_flux.z};
  const std::array<FaceFieldView, 3U> prepared_faces{
      prepared.x, prepared.y, prepared.z};
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const ConstFaceFieldView source = mechanical_faces[axis_index];
    const FaceFieldView destination = prepared_faces[axis_index];
    for (std::int32_t z = 0; z < source.extents.z; ++z)
      for (std::int32_t y = 0; y < source.extents.y; ++y)
        for (std::int32_t x = 0; x < source.extents.x; ++x)
          destination.unchecked({x, y, z}) = source.unchecked({x, y, z});
  }

  std::array<double, 6U> local_capacity{};
  std::array<double, 8U> local_mass_target{};
  std::array<double, 6U> global_capacity{};
  std::uint64_t local_inlet_hash =
      mix(kFnvOffset, UINT64_C(0x696e6c6574666c78));
  std::uint64_t local_outlet_hash =
      mix(kFnvOffset, UINT64_C(0x6f75746c65746678));
  std::uint32_t local_fixed_point_iterations = 0U;
  const Span<const double> mass_targets = impl.boundary->mass_flow_targets();
  const Span<const double> temperatures =
      impl.boundary->temperature_targets();
  const Span<const double> backflow_temperatures =
      impl.boundary->backflow_temperature_targets();
  const Span<const double> pressure_targets =
      impl.boundary->pressure_targets();
  const Span<const std::uint8_t> allow_backflow =
      impl.boundary->allow_backflow();

  const auto configured_thermo = [&](CartesianFace face,
                                     const BoundaryFacePlan& face_plan,
                                     Int3 owner, bool backflow,
                                     Real3 configured_velocity,
                                     ThermoState& thermo) noexcept {
    if (face_plan.flow_parameter >= impl.boundary->parameter_count())
      return Status{StatusCode::invalid_plan,
                    kCandidateBoundaryFinalizer};
    for (std::size_t species = 0U;
         species < input.independent_species.size; ++species) {
      std::uint32_t parameter = 0U;
      if (!scalar_parameter(*impl.boundary, face,
                            input.independent_species.data[species].field,
                            parameter))
        return Status{StatusCode::invalid_plan,
                      kCandidateBoundaryFinalizer};
      impl.composition[species] =
          backflow
              ? impl.boundary->scalar_backflow_targets().data[parameter]
              : impl.boundary->scalar_targets().data[parameter];
    }
    const double temperature =
        backflow ? backflow_temperatures.data[face_plan.flow_parameter]
                 : temperatures.data[face_plan.flow_parameter];
    double enthalpy = 0.0;
    double cp = 0.0;
    double gas = 0.0;
    Status status = impl.thermodynamics->mixture_enthalpy(
        temperature, {impl.composition.data(), impl.composition.size()},
        enthalpy, cp, gas);
    if (!status) return status;
    const double absolute_pressure =
        backflow
            ? pressure_targets.data[face_plan.flow_parameter]
            : input.absolute_pressure_reference +
                  input.pressure_perturbation.unchecked(owner, 0U);
    return impl.thermodynamics->evaluate(
        absolute_pressure, enthalpy,
        {impl.composition.data(), impl.composition.size()},
        configured_velocity, thermo, temperature);
  };

  Status local{};
  if (impl.immersed_physical_boundary_flux != nullptr) {
    local = impl.immersed_physical_boundary_flux
                ->zero_inactive_physical_boundary_flux(prepared);
    local = reductions.consensus(local);
    if (!local) return local;
  }
  for (std::size_t axis_index = 0U; axis_index < 3U && local; ++axis_index) {
    const auto axis = static_cast<CartesianAxis>(axis_index);
    const FaceFieldView destination = prepared_faces[axis_index];
    for_each_boundary_face_storage_order(
        axis, cells, [&](Int3 face_index, bool high) noexcept {
          const CartesianFace face = selected_face(axis, high);
          const BoundaryFacePlan* face_plan = nullptr;
          local = impl.boundary->face(face, face_plan);
          if (!local || face_plan == nullptr) {
            local = {StatusCode::invalid_plan,
                     kCandidateBoundaryFinalizer};
            return false;
          }
          if (!face_plan->local_owner || face_plan->periodic) return true;
          if (face_plan->flow_parameter >= impl.boundary->parameter_count()) {
            local = {StatusCode::invalid_plan,
                     kCandidateBoundaryFinalizer};
            return false;
          }
          bool control_face_active = true;
          if (impl.immersed_physical_boundary_flux != nullptr)
            local = impl.immersed_physical_boundary_flux->physical_face_active(
                axis, face_index, control_face_active);
          if (!local) return false;
          if (!control_face_active) {
            if (face_plan->flow_kind == BoundaryKind::mass_flow_inlet) {
              const double target =
                  mass_targets.data[face_plan->flow_parameter];
              if (!finite_positive(target)) {
                local = {StatusCode::numerical_failure,
                         kCandidateBoundaryNumerical};
                return false;
              }
              local_mass_target[static_cast<std::size_t>(face)] = target;
            }
            return true;
          }
          const Int3 owner = owner_cell(axis, high, face_index, cells);
          const double face_area =
              area(*impl.geometry, impl.patch, axis, face_index);
          if (!finite_positive(face_area)) {
            local = {StatusCode::numerical_failure,
                     kCandidateBoundaryNumerical};
            return false;
          }
          if (face_plan->flow_kind == BoundaryKind::velocity_inlet) {
            const Real3 prescribed = vector_parameter(
                *impl.boundary, face_plan->flow_parameter, false, false);
            ThermoState thermo;
            local = configured_thermo(face, *face_plan, owner, false,
                                      prescribed, thermo);
            const double value =
                local ? thermo.rho * component(prescribed, axis) * face_area
                      : 0.0;
            if (local && !std::isfinite(value))
              local = {StatusCode::numerical_failure,
                       kCandidateBoundaryNumerical};
            if (!local) return false;
            destination.unchecked(face_index) = value;
            local_inlet_hash = mix(local_inlet_hash, axis_index);
            local_inlet_hash = mix(local_inlet_hash, high ? 1U : 0U);
            local_inlet_hash = mix(local_inlet_hash, double_bits(value));
          } else if (face_plan->flow_kind ==
                     BoundaryKind::mass_flow_inlet) {
            const Real3 direction = vector_parameter(
                *impl.boundary, face_plan->flow_parameter, false, true);
            ThermoState thermo;
            local = configured_thermo(face, *face_plan, owner, false,
                                      direction, thermo);
            const double outward_component =
                outward_sign(high) * component(direction, axis);
            const double capacity =
                local ? -thermo.rho * outward_component * face_area : 0.0;
            if (local &&
                (!finite_positive(capacity) ||
                 !finite_positive(
                     mass_targets.data[face_plan->flow_parameter])))
              local = {StatusCode::numerical_failure,
                       kCandidateBoundaryNumerical};
            if (!local) return false;
            local_capacity[static_cast<std::size_t>(face)] += capacity;
            local_mass_target[static_cast<std::size_t>(face)] =
                mass_targets.data[face_plan->flow_parameter];
            destination.unchecked(face_index) =
                thermo.rho * component(direction, axis) * face_area;
          } else if (face_plan->flow_kind ==
                     BoundaryKind::pressure_outlet) {
            const double provisional =
                mechanical_faces[axis_index].unchecked(face_index);
            const double outward = outward_sign(high) * provisional;
            const double owner_outward_velocity =
                outward_sign(high) *
                input.velocity.unchecked(owner, axis_index);
            if (allow_backflow.data[face_plan->flow_parameter] != 0U &&
                ((owner_outward_velocity < 0.0) != (outward < 0.0))) {
              local = {StatusCode::rejected_step,
                       kCandidateBoundaryBackflow};
              return false;
            }
            local_outlet_hash = mix(local_outlet_hash, axis_index);
            local_outlet_hash = mix(local_outlet_hash, high ? 1U : 0U);
            local_outlet_hash = mix(local_outlet_hash,
                                    double_bits(provisional));
            if (outward < 0.0) {
              if (allow_backflow.data[face_plan->flow_parameter] == 0U) {
                local = {StatusCode::rejected_step,
                         kCandidateBoundaryBackflow};
                return false;
              }
              const Real3 prescribed = vector_parameter(
                  *impl.boundary, face_plan->flow_parameter, true, false);
              ThermoState thermo;
              local = configured_thermo(face, *face_plan, owner, true,
                                        prescribed, thermo);
              const double value =
                  local ? thermo.rho * component(prescribed, axis) * face_area
                        : 0.0;
              if (local &&
                  (!std::isfinite(value) ||
                   !(outward_sign(high) * value < 0.0)))
                local = {StatusCode::numerical_failure,
                         kCandidateBoundaryBackflow};
              ThermoState replayed_thermo;
              if (local)
                local = configured_thermo(face, *face_plan, owner, true,
                                          prescribed, replayed_thermo);
              const double replayed_value =
                  local ? replayed_thermo.rho * component(prescribed, axis) *
                              face_area
                        : 0.0;
              if (local &&
                  (double_bits(replayed_value) != double_bits(value) ||
                   !(outward_sign(high) * replayed_value < 0.0)))
                local = {StatusCode::numerical_failure,
                         kCandidateBoundaryBackflow};
              if (!local) return false;
              destination.unchecked(face_index) = value;
              local_fixed_point_iterations =
                  std::max(local_fixed_point_iterations, 1U);
              local_outlet_hash = mix(local_outlet_hash, UINT64_C(1));
              local_outlet_hash = mix(local_outlet_hash, double_bits(value));
            } else {
              local_outlet_hash = mix(local_outlet_hash, UINT64_C(0));
              local_outlet_hash = mix(local_outlet_hash,
                                      double_bits(provisional));
            }
          } else {
            destination.unchecked(face_index) = 0.0;
          }
          return true;
        });
  }
  local_mass_target[6U] =
      static_cast<double>(local_fixed_point_iterations);
  local = reductions.checked_sum(
      {local_capacity.data(), local_capacity.size()},
      {global_capacity.data(), global_capacity.size()}, local);
  if (!local) return local;
  local = reductions.checked_max(
      {local_mass_target.data(), local_mass_target.size()},
      {result.global_mass_target.data(), result.global_mass_target.size()});
  if (!local) return local;

  std::array<double, 6U> local_achieved{};
  for (std::size_t axis_index = 0U; axis_index < 3U && local; ++axis_index) {
    const auto axis = static_cast<CartesianAxis>(axis_index);
    const FaceFieldView destination = prepared_faces[axis_index];
    for_each_boundary_face_storage_order(
        axis, cells, [&](Int3 index, bool high) noexcept {
          const CartesianFace face = selected_face(axis, high);
          const BoundaryFacePlan* face_plan = nullptr;
          local = impl.boundary->face(face, face_plan);
          if (!local || face_plan == nullptr) {
            local = {StatusCode::invalid_plan,
                     kCandidateBoundaryFinalizer};
            return false;
          }
          if (!face_plan->local_owner || face_plan->periodic ||
              face_plan->flow_kind != BoundaryKind::mass_flow_inlet)
            return true;
          bool control_face_active = true;
          if (impl.immersed_physical_boundary_flux != nullptr)
            local = impl.immersed_physical_boundary_flux->physical_face_active(
                axis, index, control_face_active);
          if (!local) return false;
          if (!control_face_active) return true;
          const std::size_t face_slot = static_cast<std::size_t>(face);
          if (!finite_positive(global_capacity[face_slot]) ||
              !finite_positive(result.global_mass_target[face_slot])) {
            local = {StatusCode::numerical_failure,
                     kCandidateBoundaryNumerical};
            return false;
          }
          const double scale = result.global_mass_target[face_slot] /
                               global_capacity[face_slot];
          const double value = destination.unchecked(index) * scale;
          if (!std::isfinite(value)) {
            local = {StatusCode::numerical_failure,
                     kCandidateBoundaryNumerical};
            return false;
          }
          destination.unchecked(index) = value;
          local_inlet_hash = mix(local_inlet_hash, axis_index);
          local_inlet_hash = mix(local_inlet_hash, high ? 1U : 0U);
          local_inlet_hash = mix(local_inlet_hash, double_bits(value));
          local_achieved[face_slot] += -outward_sign(high) * value;
          return true;
        });
  }
  if (!local) local_achieved.fill(0.0);
  local = reductions.checked_sum(
      {local_achieved.data(), local_achieved.size()},
      {result.global_achieved.data(), result.global_achieved.size()}, local);
  if (!local) return local;
  for (std::size_t face = 0U; face < 6U && local; ++face) {
    if (result.global_mass_target[face] == 0.0) continue;
    const double scale = std::max(1.0, result.global_mass_target[face]);
    if (!std::isfinite(result.global_achieved[face]) ||
        std::abs(result.global_achieved[face] -
                 result.global_mass_target[face]) >
            64.0 * std::numeric_limits<double>::epsilon() * scale)
      local = {StatusCode::numerical_failure,
               kCandidateBoundaryNumerical};
  }
  if (local && impl.immersed_interface != nullptr) {
    local = impl.immersed_interface->zero_interface_flux(prepared);
    if (local)
      local = impl.immersed_interface->validate_interface_flux(
          as_const(prepared), 0.0);
    if (local)
      local = impl.immersed_physical_boundary_flux
                  ->validate_inactive_physical_boundary_flux(
                      as_const(prepared), 0.0);
  }
  local = reductions.consensus(local);
  if (!local) return local;

  result.prepared = prepared;
  result.local_inlet_hash = local_inlet_hash;
  result.local_outlet_hash = local_outlet_hash;
  return {};
}

void copy_flux_no_fail(ConstFaceFluxView source, FaceFluxView destination) noexcept {
  const std::array<ConstFaceFieldView, 3U> source_faces{
      source.x, source.y, source.z};
  const std::array<FaceFieldView, 3U> destination_faces{
      destination.x, destination.y, destination.z};
  for (std::size_t axis = 0U; axis < source_faces.size(); ++axis)
    for (std::int32_t z = 0; z < source_faces[axis].extents.z; ++z)
      for (std::int32_t y = 0; y < source_faces[axis].extents.y; ++y)
        for (std::int32_t x = 0; x < source_faces[axis].extents.x; ++x)
          destination_faces[axis].unchecked({x, y, z}) =
              source_faces[axis].unchecked({x, y, z});
}

}  // namespace

PressureEnergyCandidateBoundaryFinalizer::
    ~PressureEnergyCandidateBoundaryFinalizer() noexcept {
  release();
}

PressureEnergyCandidateBoundaryFinalizer::
    PressureEnergyCandidateBoundaryFinalizer(
        PressureEnergyCandidateBoundaryFinalizer&& other) noexcept
    : implementation_(other.implementation_) {
  other.implementation_ = nullptr;
}

PressureEnergyCandidateBoundaryFinalizer&
PressureEnergyCandidateBoundaryFinalizer::operator=(
    PressureEnergyCandidateBoundaryFinalizer&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = other.implementation_;
    other.implementation_ = nullptr;
  }
  return *this;
}

void PressureEnergyCandidateBoundaryFinalizer::release() noexcept {
  delete implementation_;
  implementation_ = nullptr;
}

Status PressureEnergyCandidateBoundaryFinalizer::bind(
    const PressureEnergyCandidateBoundaryFinalizerBinding& binding,
    PressureEnergyCandidateBoundaryFinalizer& out) noexcept {
  int rank = 0;
  int size = 0;
  if (binding.communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(binding.communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(binding.communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::invalid_plan, kCandidateBoundaryFinalizer};
  }
  const bool empty_candidate_donor =
      binding.candidate_pressure_correction_donors == nullptr &&
      binding.candidate_pressure_correction_donor_stage == 0U &&
      binding.candidate_pressure_correction_field == 0U &&
      binding.candidate_pressure_correction_donor_reach == 0U;
  const bool complete_candidate_donor =
      binding.candidate_pressure_correction_donors != nullptr &&
      binding.candidate_pressure_correction_donors->ready() &&
      binding.candidate_pressure_correction_donor_stage != 0U &&
      binding.candidate_pressure_correction_field != 0U &&
      binding.candidate_pressure_correction_donor_reach >= 1U &&
      binding.candidate_pressure_correction_donor_reach ==
          binding.candidate_pressure_correction_donors->reach();
  bool local_valid =
      binding.geometry != nullptr && binding.boundary != nullptr &&
      binding.kernels != nullptr && binding.thermodynamics != nullptr &&
      binding.transport != nullptr && binding.coupler != nullptr &&
      binding.geometry->fingerprint() != 0U &&
      binding.boundary->semantic_fingerprint() != 0U &&
      binding.boundary->local_layout_fingerprint() != 0U &&
      binding.kernels->fingerprint() != 0U &&
      binding.thermodynamics->fingerprint() != 0U &&
      binding.transport->fingerprint() != 0U &&
      binding.boundary->pressure_reference() ==
          PressureReferenceKind::boundary_absolute &&
      binding.boundary->geometry_fingerprint() ==
          binding.geometry->fingerprint() &&
      same(binding.patch.cells, binding.boundary->local_cells()) &&
      same(binding.patch.cells, binding.kernels->cells()) &&
      binding.coupler->matches_candidate_boundary_finalizer_binding(
          binding.geometry, binding.patch, binding.boundary,
          binding.kernels, binding.thermodynamics, binding.transport,
          binding.immersed_interface,
          binding.candidate_pressure_correction_donors,
          binding.candidate_pressure_correction_donor_stage,
          binding.candidate_pressure_correction_field,
          binding.candidate_pressure_correction_donor_reach) &&
      supported(*binding.boundary) &&
      ((binding.immersed_interface == nullptr &&
        binding.immersed_physical_boundary_flux == nullptr &&
        empty_candidate_donor) ||
       (binding.immersed_interface != nullptr &&
        binding.immersed_physical_boundary_flux != nullptr &&
        binding.immersed_interface->fingerprint() != 0U &&
        binding.immersed_physical_boundary_flux->matches(
            binding.communicator, binding.geometry, binding.patch,
            binding.immersed_interface) &&
        complete_candidate_donor));
  int local_flag = local_valid ? 1 : 0;
  int global_flag = 0;
  if (MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_MIN,
                    binding.communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCandidateBoundaryFinalizer};
  }
  if (global_flag == 0) {
    return {StatusCode::invalid_plan, kCandidateBoundaryFinalizer};
  }
  const std::uint64_t local_activity_collective =
      binding.immersed_physical_boundary_flux == nullptr
          ? PlanFingerprint{0U}
          : binding.immersed_physical_boundary_flux
                ->collective_fingerprint();
  std::uint64_t minimum_activity_collective = 0U;
  std::uint64_t maximum_activity_collective = 0U;
  const int minimum_activity_status = MPI_Allreduce(
      &local_activity_collective, &minimum_activity_collective, 1,
      MPI_UINT64_T, MPI_MIN, binding.communicator);
  const int maximum_activity_status = MPI_Allreduce(
      &local_activity_collective, &maximum_activity_collective, 1,
      MPI_UINT64_T, MPI_MAX, binding.communicator);
  if (minimum_activity_status != MPI_SUCCESS ||
      maximum_activity_status != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCandidateBoundaryFinalizer};
  if (minimum_activity_collective != maximum_activity_collective)
    return {StatusCode::invalid_plan, kCandidateBoundaryFinalizer};

  std::size_t x_count = 0U;
  std::size_t y_count = 0U;
  std::size_t z_count = 0U;
  int allocation_failure_kind = 0;
  if (!checked_face_count(binding.patch.cells, x_count, y_count, z_count))
    allocation_failure_kind = 1;
  Impl* candidate =
      allocation_failure_kind == 0 ? new (std::nothrow) Impl : nullptr;
  if (allocation_failure_kind == 0 && candidate == nullptr)
    allocation_failure_kind = 2;
  if (allocation_failure_kind == 0) {
    try {
      candidate->x_scratch.resize(x_count);
      candidate->y_scratch.resize(y_count);
      candidate->z_scratch.resize(z_count);
      candidate->collective_hashes.resize(
          static_cast<std::size_t>(size) * 8U);
      candidate->halo_contract.resize(
          5U + binding.thermodynamics->independent_species_count());
      candidate->composition.resize(
          binding.thermodynamics->independent_species_count());
    } catch (...) {
      allocation_failure_kind = 2;
    }
  }
  int collective_failure_kind = 0;
  const int allocation_consensus = MPI_Allreduce(
      &allocation_failure_kind, &collective_failure_kind, 1, MPI_INT, MPI_MAX,
      binding.communicator);
  if (allocation_consensus != MPI_SUCCESS || collective_failure_kind != 0) {
    delete candidate;
    if (allocation_consensus != MPI_SUCCESS)
      return {StatusCode::mpi_failure, kCandidateBoundaryFinalizer};
    return {collective_failure_kind == 1 ? StatusCode::invalid_plan
                                         : StatusCode::allocation_failure,
            kCandidateBoundaryFinalizer};
  }
  candidate->communicator = binding.communicator;
  candidate->rank = rank;
  candidate->size = size;
  candidate->geometry = binding.geometry;
  candidate->patch = binding.patch;
  candidate->boundary = binding.boundary;
  candidate->kernels = binding.kernels;
  candidate->thermodynamics = binding.thermodynamics;
  candidate->transport = binding.transport;
  candidate->coupler = binding.coupler;
  candidate->immersed_interface = binding.immersed_interface;
  candidate->immersed_physical_boundary_flux =
      binding.immersed_physical_boundary_flux;
  candidate->candidate_pressure_correction_donors =
      binding.candidate_pressure_correction_donors;
  candidate->candidate_pressure_correction_donor_stage =
      binding.candidate_pressure_correction_donor_stage;
  candidate->candidate_pressure_correction_field =
      binding.candidate_pressure_correction_field;
  candidate->candidate_pressure_correction_donor_reach =
      binding.candidate_pressure_correction_donor_reach;
  std::uint64_t fingerprint = mix(kFnvOffset, UINT64_C(0x70656362666e6472));
  fingerprint = mix(fingerprint, binding.geometry->fingerprint());
  fingerprint = mix(fingerprint, binding.boundary->semantic_fingerprint());
  fingerprint = mix(fingerprint, binding.kernels->fingerprint());
  fingerprint = mix(fingerprint, binding.thermodynamics->fingerprint());
  fingerprint = mix(fingerprint, binding.transport->fingerprint());
  fingerprint = mix(
      fingerprint, binding.immersed_interface == nullptr
                       ? PlanFingerprint{0U}
                       : binding.immersed_interface->fingerprint());
  fingerprint = mix(
      fingerprint,
      binding.immersed_physical_boundary_flux == nullptr
          ? PlanFingerprint{0U}
          : binding.immersed_physical_boundary_flux->local_fingerprint());
  fingerprint = mix(
      fingerprint,
      binding.immersed_physical_boundary_flux == nullptr
          ? PlanFingerprint{0U}
          : binding.immersed_physical_boundary_flux
                ->collective_fingerprint());
  fingerprint = mix(
      fingerprint,
      binding.candidate_pressure_correction_donors == nullptr
          ? PlanFingerprint{0U}
          : binding.candidate_pressure_correction_donors->fingerprint());
  fingerprint = mix(
      fingerprint, binding.candidate_pressure_correction_donor_stage);
  fingerprint = mix(fingerprint,
                    binding.candidate_pressure_correction_field);
  fingerprint = mix(
      fingerprint, binding.candidate_pressure_correction_donor_reach);
  candidate->fingerprint = fingerprint == 0U ? 1U : fingerprint;
  out.release();
  out.implementation_ = candidate;
  return {};
}

bool PressureEnergyCandidateBoundaryFinalizer::ready() const noexcept {
  return implementation_ != nullptr && implementation_->fingerprint != 0U;
}

PlanFingerprint PressureEnergyCandidateBoundaryFinalizer::fingerprint()
    const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->fingerprint;
}

Status PressureEnergyCandidateBoundaryFinalizer::close_fresh_physical_flux(
    const FreshPhysicalBoundaryFluxClosureInput& input,
    ReductionEngine& reductions) noexcept {
  if (implementation_ == nullptr)
    return {StatusCode::invalid_plan, kCandidateBoundaryFinalizer};
  Impl& impl = *implementation_;
  const Int3 cells = impl.patch.cells;
  Status local = reductions.validate_communicator(impl.communicator);
  bool valid =
      finite_positive(input.absolute_pressure_reference) &&
      detail::valid_cell_view(input.pressure_perturbation, cells, 0U, 1U,
                              0U) &&
      input.pressure_perturbation.field == impl.boundary->pressure_field() &&
      detail::valid_cell_view(input.velocity, cells, 0U, 3U, 0U) &&
      input.independent_species.size ==
          impl.thermodynamics->independent_species_count() &&
      (input.independent_species.size == 0U ||
       input.independent_species.data != nullptr) &&
      detail::valid_flux_view(input.mechanical_flux, cells,
                              input.mechanical_flux.revision) &&
      detail::valid_flux_view(as_const(input.final_flux), cells,
                              input.final_flux.revision) &&
      !input.final_flux.certificate.valid() &&
      parameter_storage_valid(*impl.boundary);

  std::size_t expected_species = 0U;
  const Span<const BoundaryTransportedField> transported =
      impl.boundary->transported_fields();
  for (std::size_t field = 0U; field < transported.size && valid; ++field) {
    if (transported.data[field].role != BoundaryScalarRole::species) continue;
    valid = expected_species < input.independent_species.size &&
            detail::valid_cell_view(
                input.independent_species.data[expected_species], cells, 0U,
                1U, 0U) &&
            transported.data[field].field ==
                input.independent_species.data[expected_species].field;
    ++expected_species;
  }
  valid = valid && expected_species == input.independent_species.size;

  const std::array<ConstFaceFieldView, 3U> mechanical_faces{
      input.mechanical_flux.x, input.mechanical_flux.y,
      input.mechanical_flux.z};
  const std::array<ConstFaceFieldView, 3U> final_faces{
      as_const(input.final_flux.x), as_const(input.final_flux.y),
      as_const(input.final_flux.z)};
  for (std::size_t axis = 0U; axis < final_faces.size() && valid; ++axis) {
    valid = valid &&
            !detail::cell_face_views_overlap(input.pressure_perturbation,
                                             final_faces[axis]) &&
            !detail::cell_face_views_overlap(input.velocity,
                                             final_faces[axis]);
    for (std::size_t species = 0U;
         species < input.independent_species.size && valid; ++species)
      valid = !detail::cell_face_views_overlap(
          input.independent_species.data[species], final_faces[axis]);
    for (ConstFaceFieldView mechanical : mechanical_faces)
      valid = valid &&
              !detail::face_views_overlap(final_faces[axis], mechanical);
    for (std::size_t other = axis + 1U;
         other < final_faces.size() && valid; ++other)
      valid = valid &&
              !detail::face_views_overlap(final_faces[axis],
                                          final_faces[other]);
  }
  for (std::int32_t z = 0; z < cells.z && valid; ++z)
    for (std::int32_t y = 0; y < cells.y && valid; ++y)
      for (std::int32_t x = 0; x < cells.x && valid; ++x) {
        const Int3 cell{x, y, z};
        valid = finite_positive(input.absolute_pressure_reference +
                                input.pressure_perturbation.unchecked(cell,
                                                                      0U));
        for (std::uint8_t component_index = 0U;
             component_index < 3U && valid; ++component_index)
          valid = std::isfinite(
              input.velocity.unchecked(cell, component_index));
      }
  for (ConstFaceFieldView face : mechanical_faces)
    for (std::int32_t z = 0; z < face.extents.z && valid; ++z)
      for (std::int32_t y = 0; y < face.extents.y && valid; ++y)
        for (std::int32_t x = 0; x < face.extents.x && valid; ++x)
          valid = std::isfinite(face.unchecked({x, y, z}));
  if (local && !valid)
    local = {StatusCode::invalid_plan, kCandidateBoundaryFinalizer};
  local = reductions.consensus(local);
  if (!local) return local;

  const PhysicalBoundaryFluxEngineInput physical_input{
      input.absolute_pressure_reference,
      input.pressure_perturbation,
      input.velocity,
      input.independent_species,
      input.mechanical_flux,
      input.final_flux.revision};
  PhysicalBoundaryFluxEngineResult physical;
  local = prepare_physical_boundary_flux(impl, physical_input, reductions,
                                         physical);
  if (!local) return local;
  copy_flux_no_fail(as_const(physical.prepared), input.final_flux);
  return {};
}

Status PressureEnergyCandidateBoundaryFinalizer::finalize(
    const PressureEnergyCandidateBoundaryFinalizeInput& input,
    ReductionEngine& reductions,
    FinalBoundaryFluxCertificate& certificate) noexcept {
  certificate = {};
  if (implementation_ == nullptr)
    return {StatusCode::invalid_plan, kCandidateBoundaryFinalizer};
  Impl& impl = *implementation_;
  const Int3 cells = impl.patch.cells;
  const auto& authority = input.authority;
  const auto& pressure = input.pressure_stage;
  const auto& velocity = input.velocity_stage;
  const auto& flux = input.flux_stage;

  Status local = reductions.validate_communicator(impl.communicator);
  const bool stage_chain =
      authority.valid() && authority.issuer_ == impl.coupler &&
      authority.scope_ ==
          PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm &&
      pressure.valid() && pressure.issuer_ == impl.coupler &&
      velocity.valid() && velocity.issuer_ == impl.coupler && flux.valid() &&
      flux.issuer_ == impl.coupler &&
      pressure.stage_lineage_ == authority.canonical_lineage_ &&
      velocity.stage_lineage_ == authority.canonical_lineage_ &&
      flux.stage_lineage_ == authority.canonical_lineage_ &&
      velocity.pressure_stage_lineage_ == pressure.canonical_lineage_ &&
      flux.velocity_lineage_ == velocity.canonical_lineage_ &&
      pressure.corrector_ == authority.corrector_ &&
      velocity.corrector_ == authority.corrector_ &&
      flux.corrector_ == authority.corrector_ &&
      pressure.alpha_ == velocity.alpha_ && pressure.alpha_ == flux.alpha_ &&
      flux.scope_ == authority.scope_;
  const bool reference_valid =
      input.pressure_reference.valid() &&
      input.pressure_reference.kind ==
          PressureReferenceKind::boundary_absolute &&
      input.pressure_reference.time == authority.pressure_.time &&
      input.pressure_reference.thermodynamics ==
          impl.thermodynamics->fingerprint() &&
      finite_positive(input.absolute_pressure_reference) &&
      double_bits(input.absolute_pressure_reference) ==
          double_bits(input.thermophysical_boundary.binding
                          .pressure_reference);
  const bool views_valid =
      detail::valid_cell_view(input.pressure_perturbation, cells, 0U, 1U,
                              1U) &&
      detail::valid_cell_view(input.enthalpy, cells, 0U, 1U, 1U) &&
      detail::valid_cell_view(input.density, cells, 0U, 1U, 1U) &&
      detail::valid_cell_view(input.temperature, cells, 0U, 1U, 1U) &&
      detail::valid_cell_view(input.velocity, cells, 0U, 3U, 1U) &&
      detail::valid_flux_view(input.mechanical_flux, cells,
                              input.mechanical_flux.revision) &&
      detail::valid_flux_view(as_const(input.final_flux), cells,
                              input.final_flux.revision) &&
      !input.final_flux.certificate.valid() &&
      input.composition_identity != 0U && input.state_halo != nullptr &&
      input.independent_species.size ==
          impl.thermodynamics->independent_species_count() &&
      (input.independent_species.size == 0U ||
       input.independent_species.data != nullptr);
  const bool stage_views_match =
      same_view_identity(input.velocity, velocity.candidate_velocity_view_) &&
      same_view_identity(input.density, flux.candidate_density_view_) &&
      same_flux_identity(input.mechanical_flux, flux.candidate_flux_view_);
  std::uint64_t candidate_velocity_numeric =
      mix(kFnvOffset, UINT64_C(0x63616e6476656c6e));
  candidate_velocity_numeric = mix_stage_field_values(
      candidate_velocity_numeric, input.velocity, cells, 3U, 0);
  std::uint64_t candidate_density_numeric =
      mix(kFnvOffset, UINT64_C(0x63616e6464656e73));
  candidate_density_numeric = mix_stage_field_values(
      candidate_density_numeric, input.density, cells, 1U, 1);
  std::uint64_t mechanical_flux_numeric =
      mix(kFnvOffset, UINT64_C(0x6d656368666c786e));
  mechanical_flux_numeric = mix_stage_flux_values(
      mechanical_flux_numeric, input.mechanical_flux);
  const bool stage_numeric_match =
      (candidate_velocity_numeric == 0U ? 1U
                                        : candidate_velocity_numeric) ==
          velocity.candidate_velocity_numeric_ &&
      (candidate_density_numeric == 0U ? 1U
                                       : candidate_density_numeric) ==
          flux.candidate_density_numeric_ &&
      (mechanical_flux_numeric == 0U ? 1U
                                     : mechanical_flux_numeric) ==
          flux.mechanical_flux_numeric_;
  const bool thermo_views_match =
      same_view_payload_identity(
          input.pressure_perturbation,
          input.thermophysical_boundary.binding.pressure_perturbation) &&
      input.thermophysical_boundary.binding.pressure_perturbation.field ==
          impl.boundary->pressure_field() &&
      same_view_payload_identity(
          input.enthalpy,
          input.thermophysical_boundary.binding.enthalpy) &&
      input.thermophysical_boundary.binding.enthalpy.field ==
          impl.boundary->enthalpy_field() &&
      same_view_identity(input.density,
                         input.thermophysical_boundary.binding.density) &&
      input.independent_species.size ==
          input.thermophysical_boundary.binding.independent_species.size;
  bool species_valid =
      input.independent_species.size ==
          input.thermophysical_boundary.binding.independent_species.size &&
      (input.independent_species.size == 0U ||
       (input.independent_species.data != nullptr &&
        input.thermophysical_boundary.binding.independent_species.data !=
            nullptr));
  for (std::size_t species = 0U;
       species < input.independent_species.size && species_valid; ++species) {
    const ConstFieldView field = input.independent_species.data[species];
    const ConstFieldView semantic =
        input.thermophysical_boundary.binding.independent_species
            .data[species];
    species_valid =
        detail::valid_cell_view(field, cells, 0U, 1U, 1U) &&
        semantic.field != 0U && same_view_payload_identity(field, semantic);
    for (std::size_t prior = 0U; prior < species && species_valid; ++prior) {
      species_valid =
          semantic.field !=
              input.thermophysical_boundary.binding.independent_species
                  .data[prior]
                  .field &&
          field.field != input.independent_species.data[prior].field;
    }
  }
  const std::array<ConstFaceFieldView, 3U> mechanical_identity_faces{
      input.mechanical_flux.x, input.mechanical_flux.y,
      input.mechanical_flux.z};
  const std::array<ConstFaceFieldView, 3U> final_identity_faces{
      as_const(input.final_flux.x), as_const(input.final_flux.y),
      as_const(input.final_flux.z)};
  bool distinct_flux = true;
  for (std::size_t axis = 0U; axis < final_identity_faces.size(); ++axis) {
    const std::array<ConstFieldView, 5U> state_inputs{{
        input.pressure_perturbation, input.enthalpy, input.density,
        input.temperature, input.velocity}};
    for (ConstFieldView state : state_inputs)
      distinct_flux =
          distinct_flux &&
          !detail::cell_face_views_overlap(state, final_identity_faces[axis]);
    for (std::size_t species = 0U;
         species < input.independent_species.size; ++species)
      distinct_flux =
          distinct_flux &&
          !detail::cell_face_views_overlap(
              input.independent_species.data[species],
              final_identity_faces[axis]);
    for (ConstFaceFieldView mechanical : mechanical_identity_faces)
      distinct_flux =
          distinct_flux &&
          !detail::face_views_overlap(final_identity_faces[axis],
                                      mechanical);
    for (std::size_t other = axis + 1U;
         other < final_identity_faces.size(); ++other)
      distinct_flux =
          distinct_flux &&
          !detail::face_views_overlap(final_identity_faces[axis],
                                      final_identity_faces[other]);
  }
  const bool donor_chain = [&]() noexcept {
    if (impl.immersed_interface == nullptr) {
      return impl.immersed_physical_boundary_flux == nullptr &&
             velocity.candidate_donor_fingerprint_ == 0U &&
             velocity.ibm_geometry_fingerprint_ == 0U &&
             flux.ibm_geometry_fingerprint_ == 0U;
    }
    return impl.candidate_pressure_correction_donors != nullptr &&
           impl.immersed_physical_boundary_flux != nullptr &&
           impl.immersed_physical_boundary_flux->matches(
               impl.communicator, impl.geometry, impl.patch,
               impl.immersed_interface) &&
           impl.candidate_pressure_correction_donors->ready() &&
           velocity.candidate_donor_fingerprint_ ==
               impl.candidate_pressure_correction_donors->fingerprint() &&
           velocity.candidate_donor_stage_ ==
               impl.candidate_pressure_correction_donor_stage &&
           velocity.candidate_donor_field_ ==
               impl.candidate_pressure_correction_field &&
           velocity.candidate_donor_reach_ ==
               impl.candidate_pressure_correction_donor_reach &&
           velocity.ibm_geometry_fingerprint_ ==
               impl.immersed_interface->fingerprint() &&
           flux.ibm_geometry_fingerprint_ ==
               impl.immersed_interface->fingerprint();
  }();
  if (local &&
      (!stage_chain || !reference_valid || !views_valid ||
       !stage_views_match || !stage_numeric_match || !thermo_views_match ||
       !species_valid ||
       !distinct_flux || !donor_chain || !parameter_storage_valid(*impl.boundary)))
    local = {StatusCode::invalid_plan, kCandidateBoundaryFinalizer};

  local = reductions.consensus(local);
  if (!local) return local;
  const PlanFingerprint composition_numeric_provenance =
      composition_numeric_local_provenance(input.independent_species, cells);

  const BoundaryThermophysicalGhostContext thermo_context{
      authority.pressure_.time,
      impl.geometry->fingerprint(),
      input.pressure_reference.pressure_reference,
      authority.pressure_.numeric_boundary,
      authority.corrector_ == 1U
          ? BoundaryThermophysicalGhostPhase::corrector_one
          : BoundaryThermophysicalGhostPhase::corrector_two};
  local = input.thermophysical_boundary.certificate.matches(
              *impl.boundary, thermo_context,
              input.thermophysical_boundary.binding)
              ? Status{}
              : Status{StatusCode::invalid_plan,
                       kCandidateBoundaryFinalizer};
  local = reductions.consensus(local);
  if (!local) return local;

  impl.halo_contract[0U] = {input.pressure_perturbation.field, 1U, 1U};
  impl.halo_contract[1U] = {input.enthalpy.field, 1U, 1U};
  impl.halo_contract[2U] = {input.density.field, 1U, 1U};
  impl.halo_contract[3U] = {input.temperature.field, 1U, 1U};
  impl.halo_contract[4U] = {input.velocity.field, 1U, 3U};
  for (std::size_t species = 0U; species < input.independent_species.size;
       ++species)
    impl.halo_contract[5U + species] = {
        input.independent_species.data[species].field, 1U, 1U};
  std::sort(impl.halo_contract.begin(), impl.halo_contract.end(),
            [](const HaloFieldSpec& left,
               const HaloFieldSpec& right) noexcept {
              return left.field < right.field;
            });
  local = input.state_halo->ready() &&
                  input.state_halo->instance_identity() != 0U &&
                  input.state_halo->instance_identity() !=
                      pressure.halo_instance_
              ? input.state_halo->validate_contract(
                    impl.communicator, impl.patch,
                    {impl.halo_contract.data(), impl.halo_contract.size()},
                    impl.boundary->halo_topology())
              : Status{StatusCode::invalid_plan,
                       kCandidateBoundaryFinalizer};
  local = reductions.consensus(local);
  if (!local) return local;
  const auto halo_current = [&](ConstFieldView field) noexcept {
    return input.state_halo->ghost_revision(field.field) == field.revision;
  };
  bool all_halos_current = halo_current(input.pressure_perturbation) &&
                           halo_current(input.enthalpy) &&
                           halo_current(input.density) &&
                           halo_current(input.temperature) &&
                           halo_current(input.velocity);
  for (std::size_t species = 0U;
       species < input.independent_species.size && all_halos_current;
       ++species)
    all_halos_current = halo_current(input.independent_species.data[species]);
  local = reductions.consensus(
      all_halos_current
          ? Status{}
          : Status{StatusCode::invalid_plan,
                   kCandidateBoundaryFinalizer});
  if (!local) return local;

  for (std::int32_t z = 0; z < cells.z && local; ++z)
    for (std::int32_t y = 0; y < cells.y && local; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double absolute_pressure =
            input.absolute_pressure_reference +
            input.pressure_perturbation.unchecked(cell, 0U);
        bool finite_velocity = true;
        for (std::uint8_t component_index = 0U; component_index < 3U;
             ++component_index)
          finite_velocity =
              finite_velocity &&
              std::isfinite(input.velocity.unchecked(cell, component_index));
        if (!finite_positive(absolute_pressure) ||
            !std::isfinite(input.enthalpy.unchecked(cell, 0U)) ||
            !finite_positive(input.density.unchecked(cell, 0U)) ||
            !finite_positive(input.temperature.unchecked(cell, 0U)) ||
            !finite_velocity) {
          local = {StatusCode::numerical_failure,
                   kCandidateBoundaryNumerical};
          break;
        }
        for (std::size_t species = 0U;
             species < input.independent_species.size; ++species)
          impl.composition[species] =
              input.independent_species.data[species].unchecked(cell, 0U);
        const Real3 cell_velocity{
            input.velocity.unchecked(cell, 0U),
            input.velocity.unchecked(cell, 1U),
            input.velocity.unchecked(cell, 2U)};
        ThermoState replayed;
        local = impl.thermodynamics->evaluate(
            absolute_pressure, input.enthalpy.unchecked(cell, 0U),
            {impl.composition.data(), impl.composition.size()},
            cell_velocity, replayed,
            input.temperature.unchecked(cell, 0U));
        if (!local) break;
        const double candidate_density =
            input.density.unchecked(cell, 0U);
        const double candidate_temperature =
            input.temperature.unchecked(cell, 0U);
        const double density_scale =
            std::max({1.0, std::abs(candidate_density),
                      std::abs(replayed.rho)});
        const double temperature_scale =
            std::max({1.0, std::abs(candidate_temperature),
                      std::abs(replayed.temperature)});
        if (std::abs(candidate_density - replayed.rho) >
                64.0 * std::numeric_limits<double>::epsilon() *
                    density_scale ||
            std::abs(candidate_temperature - replayed.temperature) >
                64.0 * std::numeric_limits<double>::epsilon() *
                    temperature_scale) {
          local = {StatusCode::invalid_plan,
                   kCandidateBoundaryFinalizer};
          break;
        }
      }
  local = reductions.consensus(local);
  if (!local) return local;

  const PhysicalBoundaryFluxEngineInput physical_input{
      input.absolute_pressure_reference,
      input.pressure_perturbation,
      input.velocity,
      input.thermophysical_boundary.binding.independent_species,
      input.mechanical_flux,
      input.final_flux.revision};
  PhysicalBoundaryFluxEngineResult physical;
  local = prepare_physical_boundary_flux(impl, physical_input, reductions,
                                         physical);
  if (!local) return local;
  const FaceFluxView prepared = physical.prepared;
  const std::array<FaceFieldView, 3U> prepared_faces{
      prepared.x, prepared.y, prepared.z};
  const std::uint64_t local_inlet_hash = physical.local_inlet_hash;
  const std::uint64_t local_outlet_hash = physical.local_outlet_hash;
  const std::array<double, 8U>& global_mass_target =
      physical.global_mass_target;
  const std::array<double, 6U>& global_achieved = physical.global_achieved;
  const std::uint64_t local_state = state_local_provenance(
      input, cells, impl.thermodynamics->fingerprint(),
      composition_numeric_provenance);
  const std::uint64_t local_final_flux =
      flux_local_provenance(as_const(prepared));
  const std::uint64_t local_ibm =
      mix(mix(kFnvOffset, UINT64_C(0x69626d66696e616c)),
          impl.immersed_interface == nullptr
              ? PlanFingerprint{0U}
              : impl.immersed_interface->fingerprint());
  const std::uint64_t local_ibm_activity =
      mix(mix(local_ibm, UINT64_C(0x7068797361637469)),
          impl.immersed_physical_boundary_flux == nullptr
              ? PlanFingerprint{0U}
              : impl.immersed_physical_boundary_flux->local_fingerprint());
  const std::uint64_t local_ibm_control =
      mix(local_ibm_activity,
          impl.immersed_physical_boundary_flux == nullptr
              ? PlanFingerprint{0U}
              : impl.immersed_physical_boundary_flux
                    ->collective_fingerprint());
  const std::array<std::uint64_t, 5U> local_hashes{
      local_state, local_final_flux, local_inlet_hash, local_outlet_hash,
      local_ibm_control};
  std::array<PlanFingerprint, 5U> collective_hashes{};
  Status provenance_collective_status;
  for (std::size_t item = 0U; item < local_hashes.size(); ++item) {
    std::uint64_t* gathered =
        impl.collective_hashes.data() + item * static_cast<std::size_t>(impl.size);
    const int gather_result = MPI_Allgather(
        local_hashes.data() + item, 1, MPI_UINT64_T, gathered, 1,
        MPI_UINT64_T, impl.communicator);
    if (gather_result != MPI_SUCCESS) {
      provenance_collective_status = {
          StatusCode::mpi_failure, kCandidateBoundaryFinalizer};
      continue;
    }
    std::uint64_t collective =
        mix(kFnvOffset, UINT64_C(0x7065636266636f6c) + item);
    collective = mix(collective, static_cast<std::uint64_t>(impl.size));
    for (int rank = 0; rank < impl.size; ++rank)
      collective = mix(collective, gathered[static_cast<std::size_t>(rank)]);
    collective_hashes[item] = collective == 0U ? 1U : collective;
  }
  provenance_collective_status =
      reductions.consensus(provenance_collective_status);
  if (!provenance_collective_status) return provenance_collective_status;

  std::uint64_t state_binding =
      mix(kFnvOffset, UINT64_C(0x7065636266736372));
  state_binding = mix_view_binding(state_binding,
                                   input.pressure_perturbation);
  state_binding = mix_view_binding(state_binding, input.enthalpy);
  state_binding = mix_view_binding(state_binding, input.density);
  state_binding = mix_view_binding(state_binding, input.temperature);
  state_binding = mix_view_binding(state_binding, input.velocity);
  for (std::size_t species = 0U; species < input.independent_species.size;
       ++species)
    state_binding = mix_view_binding(
        state_binding, input.independent_species.data[species]);
  std::uint64_t halo_lineage =
      mix(kFnvOffset, UINT64_C(0x706563626668616c));
  halo_lineage = mix(halo_lineage,
                     input.state_halo->instance_identity());
  for (const HaloFieldSpec& spec : impl.halo_contract) {
    halo_lineage = mix(halo_lineage, spec.field);
    halo_lineage = mix(halo_lineage,
                       input.state_halo->ghost_revision(spec.field));
  }
  std::uint64_t ibm_donor_lineage =
      mix(kFnvOffset, UINT64_C(0x7065636266646f6e));
  ibm_donor_lineage = mix(
      ibm_donor_lineage,
      impl.candidate_pressure_correction_donors == nullptr
          ? PlanFingerprint{0U}
          : impl.candidate_pressure_correction_donors->fingerprint());
  ibm_donor_lineage = mix(
      ibm_donor_lineage, impl.candidate_pressure_correction_donor_stage);
  ibm_donor_lineage = mix(
      ibm_donor_lineage, impl.candidate_pressure_correction_field);
  ibm_donor_lineage = mix(
      ibm_donor_lineage, impl.candidate_pressure_correction_donor_reach);
  std::uint64_t canonical =
      mix(authority.canonical_lineage_, UINT64_C(0x7065636266636572));
  canonical = mix(canonical, velocity.canonical_lineage_);
  canonical = mix(canonical, flux.canonical_lineage_);
  canonical = mix(canonical, collective_hashes[0U]);
  canonical = mix(canonical, collective_hashes[1U]);
  canonical = mix(canonical, collective_hashes[2U]);
  canonical = mix(canonical, collective_hashes[3U]);
  canonical = mix(canonical, collective_hashes[4U]);
  canonical = mix(canonical, input.pressure_reference.pressure_reference);
  canonical = mix(canonical, double_bits(pressure.alpha_));
  std::uint64_t scratch = mix(state_binding, input.final_flux.revision);
  scratch = mix(scratch,
                reinterpret_cast<std::uintptr_t>(input.final_flux.x.base));
  scratch = mix(scratch,
                reinterpret_cast<std::uintptr_t>(input.final_flux.y.base));
  scratch = mix(scratch,
                reinterpret_cast<std::uintptr_t>(input.final_flux.z.base));

  double total_target = 0.0;
  double total_achieved = 0.0;
  for (std::size_t face = 0U; face < 6U; ++face) {
    total_target += global_mass_target[face];
    total_achieved += global_achieved[face];
  }
  FinalBoundaryFluxCertificate issued;
  issued.issuer_ = this;
  issued.coupler_ = impl.coupler;
  issued.scope_ = authority.scope_;
  issued.stage_lineage_ = authority.canonical_lineage_;
  issued.velocity_lineage_ = velocity.canonical_lineage_;
  issued.mechanical_flux_lineage_ = flux.canonical_lineage_;
  issued.nonphysical_flux_provenance_ =
      flux.nonphysical_flux_provenance_;
  issued.pressure_outlet_provisional_provenance_ =
      flux.pressure_outlet_provisional_provenance_;
  issued.boundary_semantic_ = impl.boundary->semantic_fingerprint();
  issued.boundary_layout_ = impl.boundary->local_layout_fingerprint();
  issued.pressure_reference_ = input.pressure_reference;
  issued.absolute_pressure_reference_ = input.absolute_pressure_reference;
  issued.target_time_ = authority.pressure_.time;
  issued.alpha_ = pressure.alpha_;
  issued.candidate_pressure_ =
      make_piso_field_revision_identity(input.pressure_perturbation);
  issued.candidate_enthalpy_ =
      make_piso_field_revision_identity(input.enthalpy);
  issued.candidate_density_ =
      make_piso_field_revision_identity(input.density);
  issued.candidate_temperature_ =
      make_piso_field_revision_identity(input.temperature);
  issued.candidate_velocity_ =
      make_piso_field_revision_identity(input.velocity);
  issued.composition_identity_ = input.composition_identity;
  issued.composition_numeric_provenance_ =
      composition_numeric_provenance;
  issued.independent_species_views_ = input.independent_species;
  issued.independent_species_count_ = input.independent_species.size;
  issued.candidate_state_provenance_ = collective_hashes[0U];
  issued.candidate_state_binding_ = state_binding == 0U ? 1U : state_binding;
  issued.state_halo_instance_ = input.state_halo->instance_identity();
  issued.state_halo_lineage_ = halo_lineage == 0U ? 1U : halo_lineage;
  issued.thermophysical_boundary_semantics_ =
      input.thermophysical_boundary.certificate.collective_semantics();
  issued.thermophysical_boundary_target_ =
      input.thermophysical_boundary.certificate.collective_target();
  issued.thermophysical_boundary_collective_lineage_ =
      input.thermophysical_boundary.certificate.collective_lineage();
  issued.thermophysical_boundary_rank_local_lineage_ =
      input.thermophysical_boundary.certificate.rank_local_lineage();
  issued.thermophysical_boundary_rank_local_binding_ =
      input.thermophysical_boundary.certificate.rank_local_binding();
  issued.physical_state_halo_lineage_ =
      mix(issued.state_halo_lineage_,
          issued.thermophysical_boundary_rank_local_binding_);
  issued.face_closure_lineage_ =
      issued.thermophysical_boundary_collective_lineage_;
  issued.inlet_flux_provenance_ = collective_hashes[2U];
  issued.outlet_backflow_provenance_ = collective_hashes[3U];
  issued.ibm_donor_lineage_ =
      ibm_donor_lineage == 0U ? 1U : ibm_donor_lineage;
  issued.ibm_geometry_lineage_ = collective_hashes[4U];
  issued.ibm_zero_interface_lineage_ =
      mix(flux.ibm_interface_provenance_, collective_hashes[4U]);
  issued.inlet_mass_target_ = total_target;
  issued.inlet_mass_achieved_ = total_achieved;
  issued.outlet_fixed_point_iterations_ = static_cast<std::uint32_t>(
      global_mass_target[6U]);
  issued.final_flux_view_ = as_const(input.final_flux);
  issued.final_flux_revision_ = input.final_flux.revision;
  issued.final_flux_storage_ = input.final_flux.x.storage_identity;
  issued.final_flux_revision_domain_ = input.final_flux.x.revision_domain;
  issued.final_flux_provenance_ = collective_hashes[1U];
  issued.canonical_lineage_ = canonical == 0U ? 1U : canonical;
  issued.scratch_binding_ = scratch == 0U ? 1U : scratch;
  issued.corrector_ = authority.corrector_;
  local = reductions.consensus(
      issued.valid()
          ? Status{}
          : Status{StatusCode::invalid_plan,
                   kCandidateBoundaryFinalizer});
  if (!local) return local;

  // Every fallible local/collective check is complete.  The publication
  // below is a bounded no-fail copy into caller-owned candidate scratch.
  const std::array<FaceFieldView, 3U> output_faces{
      input.final_flux.x, input.final_flux.y, input.final_flux.z};
  for (std::size_t axis_index = 0U; axis_index < 3U; ++axis_index) {
    const FaceFieldView source = prepared_faces[axis_index];
    const FaceFieldView destination = output_faces[axis_index];
    for (std::int32_t z = 0; z < source.extents.z; ++z)
      for (std::int32_t y = 0; y < source.extents.y; ++y)
        for (std::int32_t x = 0; x < source.extents.x; ++x)
          destination.unchecked({x, y, z}) =
              source.unchecked({x, y, z});
  }
  certificate = issued;
  return {};
}

}  // namespace hundun::v04
