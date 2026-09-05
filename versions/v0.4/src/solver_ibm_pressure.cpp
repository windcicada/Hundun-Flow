// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <utility>
#include <vector>

#include "field_view_interval_detail.hpp"
#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kIbmPressurePlan = 9201U;
constexpr std::uint32_t kIbmPressureApply = 9202U;
constexpr std::uint32_t kIbmPressureNumerical = 9203U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

bool same_shape(Int3 left, Int3 right) noexcept {
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
         same_shape(left.local_shape, right.local_shape) &&
         left.operator_class == right.operator_class;
}

template <class T>
bool valid_scalar(BasicFieldView<T> field, Int3 cells,
                  std::uint8_t ghosts) noexcept {
  detail::FieldStorageInterval interval{};
  return field.base != nullptr && same_shape(field.interior, cells) &&
         field.components == 1U && field.ghosts.x >= ghosts &&
         field.ghosts.y >= ghosts && field.ghosts.z >= ghosts &&
         field.storage_identity != 0U && field.revision_domain != 0U &&
         detail::field_storage_interval(field, interval);
}

bool valid_face(ConstFaceFieldView field, CartesianAxis axis,
                Int3 cells) noexcept {
  detail::FieldStorageInterval interval{};
  Int3 extents = cells;
  if (axis == CartesianAxis::x) {
    ++extents.x;
  } else if (axis == CartesianAxis::y) {
    ++extents.y;
  } else {
    ++extents.z;
  }
  return field.axis == axis && same_shape(field.extents, extents) &&
         detail::face_storage_interval(field, interval);
}

double coefficient(ConstFaceFieldView x, ConstFaceFieldView y,
                   ConstFaceFieldView z, const ImmersedLink& link) noexcept {
  const Int3 cell = link.fluid_local_index;
  switch (link.direction) {
    case ImmersedFaceDirection::x_negative:
      return x.unchecked(cell);
    case ImmersedFaceDirection::x_positive:
      return x.unchecked({cell.x + 1, cell.y, cell.z});
    case ImmersedFaceDirection::y_negative:
      return y.unchecked(cell);
    case ImmersedFaceDirection::y_positive:
      return y.unchecked({cell.x, cell.y + 1, cell.z});
    case ImmersedFaceDirection::z_negative:
      return z.unchecked(cell);
    case ImmersedFaceDirection::z_positive:
      return z.unchecked({cell.x, cell.y, cell.z + 1});
  }
  return -1.0;
}

}  // namespace

struct IbmPressureOperator::Impl {
  const LinearOperator* regular{};
  const EBTopology* topology{};
  const BoundaryStencilPlan* boundary{};
  ConstFaceFieldView x{};
  ConstFaceFieldView y{};
  ConstFaceFieldView z{};
  RevisionToken geometry{};
  PlanFingerprint fingerprint{};
  PlanFingerprint collective_static{};
  Int3 cells{};
  std::vector<Int3> solid_cells;
  LinearOperatorFailureProvenance failure_provenance{};
};

Status IbmPressureOperator::bind_internal(
    const LinearOperator& regular, Int3 cells, const EBTopology& topology,
    const BoundaryStencilPlan& boundary, ConstFaceFieldView x_coefficient,
    ConstFaceFieldView y_coefficient, ConstFaceFieldView z_coefficient,
    RevisionToken geometry_revision, RemoteDonorExchangePlan* donor_exchange,
    StageId donor_stage, IbmPressureOperator& out) noexcept try {
  (void)donor_exchange;
  (void)donor_stage;
  const std::size_t expected_cells =
      cells.x > 0 && cells.y > 0 && cells.z > 0
          ? static_cast<std::size_t>(cells.x) *
                static_cast<std::size_t>(cells.y) *
                static_cast<std::size_t>(cells.z)
          : 0U;
  if (expected_cells == 0U || topology.region().size != expected_cells ||
      topology.fingerprint() == 0U || topology.geometry_revision() == 0U ||
      topology.geometry_revision() != geometry_revision ||
      boundary.fingerprint() == 0U ||
      boundary.links().size != topology.links().size ||
      boundary.reconstruction().fingerprint() == 0U ||
      !valid_face(x_coefficient, CartesianAxis::x, cells) ||
      !valid_face(y_coefficient, CartesianAxis::y, cells) ||
      !valid_face(z_coefficient, CartesianAxis::z, cells)) {
    return {StatusCode::invalid_plan, kIbmPressurePlan};
  }
  auto candidate = std::make_unique<Impl>();
  candidate->regular = &regular;
  candidate->topology = &topology;
  candidate->boundary = &boundary;
  candidate->x = x_coefficient;
  candidate->y = y_coefficient;
  candidate->z = z_coefficient;
  candidate->geometry = geometry_revision;
  candidate->cells = cells;
  const Span<const std::uint8_t> region = topology.region();
  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x, ++flat) {
        if (region.data[flat] ==
            static_cast<std::uint8_t>(RegionFlag::solid)) {
          candidate->solid_cells.push_back({x, y, z});
        }
      }
    }
  }
  std::uint64_t local = kFnvOffset;
  local = mix(local, topology.fingerprint());
  local = mix(local, boundary.fingerprint());
  local = mix(local, geometry_revision);
  local = mix(local, x_coefficient.storage_identity);
  local = mix(local, y_coefficient.storage_identity);
  local = mix(local, z_coefficient.storage_identity);
  candidate->fingerprint = local == 0U ? 1U : local;
  std::uint64_t collective = kFnvOffset;
  collective = mix(collective, topology.geometry_fingerprint());
  collective = mix(collective, topology.surface_fingerprint());
  collective = mix(collective, geometry_revision);
  candidate->collective_static = collective == 0U ? 1U : collective;
  out.release();
  out.implementation_ = candidate.release();
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kIbmPressurePlan};
} catch (...) {
  return {StatusCode::invalid_plan, kIbmPressurePlan};
}

IbmPressureOperator::~IbmPressureOperator() noexcept { release(); }

IbmPressureOperator::IbmPressureOperator(IbmPressureOperator&& other) noexcept
    : implementation_(std::exchange(other.implementation_, nullptr)) {}

IbmPressureOperator& IbmPressureOperator::operator=(
    IbmPressureOperator&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = std::exchange(other.implementation_, nullptr);
  }
  return *this;
}

void IbmPressureOperator::release() noexcept {
  delete std::exchange(implementation_, nullptr);
}

Status IbmPressureOperator::bind(
    const LinearOperator& regular, const EBTopology& topology,
    const BoundaryStencilPlan& boundary,
    ConstFaceFieldView x_coefficient, ConstFaceFieldView y_coefficient,
    ConstFaceFieldView z_coefficient, RevisionToken geometry_revision,
    IbmPressureOperator& out) noexcept {
  const LinearOperatorCertificate base = regular.certificate();
  const Int3 cells = base.local_shape;
  if (base.identity.fingerprint == 0U || base.collective_fingerprint == 0U ||
      cells.x <= 0 || cells.y <= 0 || cells.z <= 0) {
    return {StatusCode::invalid_plan, kIbmPressurePlan};
  }
  return bind_internal(regular, cells, topology, boundary, x_coefficient,
                       y_coefficient, z_coefficient, geometry_revision,
                       nullptr, 0U, out);
}

Status IbmPressureOperator::bind_lifecycle(
    const LinearOperator& regular, Int3 local_shape,
    const EBTopology& topology, const BoundaryStencilPlan& boundary,
    ConstFaceFieldView x_coefficient, ConstFaceFieldView y_coefficient,
    ConstFaceFieldView z_coefficient, RevisionToken geometry_revision,
    RemoteDonorExchangePlan* donor_exchange,
    StageId donor_stage,
    IbmPressureOperator& out) noexcept {
  return bind_internal(regular, local_shape, topology, boundary,
                       x_coefficient, y_coefficient, z_coefficient,
                       geometry_revision, donor_exchange, donor_stage, out);
}

LinearOperatorCertificate IbmPressureOperator::certificate() const noexcept {
  if (implementation_ == nullptr) {
    return {};
  }
  LinearOperatorCertificate result = implementation_->regular->certificate();
  if (result.identity.fingerprint == 0U) {
    return {};
  }
  std::uint64_t collective = kFnvOffset;
  collective = mix(collective, result.collective_fingerprint);
  collective = mix(collective, implementation_->collective_static);
  result.collective_fingerprint = collective == 0U ? 1U : collective;
  // An impermeable Cartesian IBM link is removed from the conservative
  // pressure-flux graph.  Solid rows are isolated identities, so the exact
  // active-domain operator retains the regular operator's SPD class.
  result.operator_class = LinearOperatorClass::spd;
  return result;
}

LinearOperatorFailureProvenance
IbmPressureOperator::failure_provenance() const noexcept {
  return implementation_ == nullptr
             ? LinearOperatorFailureProvenance{}
             : implementation_->failure_provenance;
}

Status IbmPressureOperator::apply(FieldView x, FieldView y) const noexcept {
  if (implementation_ != nullptr) {
    implementation_->failure_provenance = {};
  }
  if (implementation_ == nullptr || certificate().identity.fingerprint == 0U ||
      !valid_scalar(x, implementation_->cells, 1U) ||
      !valid_scalar(y, implementation_->cells, 0U) ||
      detail::field_views_overlap(as_const(x), as_const(y))) {
    return {StatusCode::invalid_plan, kIbmPressureApply};
  }
  Status status = implementation_->regular->apply(x, y);
  if (!status) {
    implementation_->failure_provenance =
        implementation_->regular->failure_provenance();
    return status;
  }
  return decorate_pressure_action(x, y);
}

Status IbmPressureOperator::certify_pressure_energy_shared_halo(
    const PressureEnergyPressureFluxOperator& energy_pressure,
    PressureEnergySharedPressureCertificate& certificate) const noexcept {
  certificate = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kIbmPressurePlan};
  }
  const auto* regular =
      dynamic_cast<const PressureLinearOperator*>(implementation_->regular);
  if (regular == nullptr) {
    return {StatusCode::invalid_plan, kIbmPressurePlan};
  }
  return regular->certify_pressure_energy_shared_halo(
      *this, PressureEnergySharedPressureScope::ibm_decorated,
      energy_pressure, certificate);
}

Status IbmPressureOperator::exchange_pressure_energy_shared_input(
    FieldView input,
    const PressureEnergySharedPressureCertificate& certificate,
    PressureEnergySharedPressureInputCertificate& input_certificate) const
    noexcept {
  input_certificate = {};
  if (implementation_ == nullptr || !certificate.valid() ||
      certificate.scope_ !=
          PressureEnergySharedPressureScope::ibm_decorated ||
      certificate.continuity_owner_ != this ||
      !same_certificate(certificate.continuity_pressure_,
                        this->certificate())) {
    return {StatusCode::invalid_plan, kIbmPressureApply};
  }
  const auto* regular =
      dynamic_cast<const PressureLinearOperator*>(implementation_->regular);
  if (regular == nullptr) {
    return {StatusCode::invalid_plan, kIbmPressureApply};
  }
  const Status status = regular->exchange_pressure_energy_shared_input(
      input, certificate, input_certificate);
  if (!status) {
    implementation_->failure_provenance = regular->failure_provenance();
  }
  return status;
}

Status IbmPressureOperator::apply_pressure_energy_shared_input(
    FieldView input, FieldView output,
    const PressureEnergySharedPressureCertificate& certificate,
    const PressureEnergySharedPressureInputCertificate& input_certificate)
    const noexcept {
  if (implementation_ != nullptr) implementation_->failure_provenance = {};
  if (implementation_ == nullptr || !certificate.valid() ||
      certificate.scope_ !=
          PressureEnergySharedPressureScope::ibm_decorated ||
      certificate.continuity_owner_ != this ||
      !same_certificate(certificate.continuity_pressure_,
                        this->certificate())) {
    return {StatusCode::invalid_plan, kIbmPressureApply};
  }
  const auto* regular =
      dynamic_cast<const PressureLinearOperator*>(implementation_->regular);
  if (regular == nullptr) {
    return {StatusCode::invalid_plan, kIbmPressureApply};
  }
  Status status = regular->apply_pressure_energy_shared_input(
      input, output, certificate, input_certificate);
  if (!status) {
    implementation_->failure_provenance = regular->failure_provenance();
    return status;
  }
  return decorate_pressure_action(input, output);
}

Status IbmPressureOperator::decorate_pressure_action(
    FieldView x, FieldView y) const noexcept {
  if (implementation_ == nullptr || certificate().identity.fingerprint == 0U ||
      !valid_scalar(x, implementation_->cells, 1U) ||
      !valid_scalar(y, implementation_->cells, 0U) ||
      detail::field_views_overlap(as_const(x), as_const(y))) {
    return {StatusCode::invalid_plan, kIbmPressureApply};
  }
  const Span<const std::uint8_t> region = implementation_->topology->region();
  const std::size_t expected =
      static_cast<std::size_t>(implementation_->cells.x) *
      implementation_->cells.y * implementation_->cells.z;
  if (region.size != expected) {
    return {StatusCode::invalid_plan, kIbmPressureApply};
  }
  for (const Int3 cell : implementation_->solid_cells) {
    y.unchecked(cell, 0U) = x.unchecked(cell, 0U);
  }
  const Span<const ImmersedLink> links = implementation_->topology->links();
  for (std::size_t index = 0U; index < links.size; ++index) {
    const ImmersedLink& link = links.data[index];
    const double face = coefficient(implementation_->x, implementation_->y,
                                    implementation_->z, link);
    const double solid = x.unchecked(link.solid_local_index, 0U);
    const double fluid = x.unchecked(link.fluid_local_index, 0U);
    const double correction = face * (solid - fluid);
    if (!std::isfinite(face) || face < 0.0 ||
        !std::isfinite(solid) || !std::isfinite(fluid) ||
        !std::isfinite(correction)) {
      const Status status{StatusCode::numerical_failure,
                          kIbmPressureNumerical};
      implementation_->failure_provenance = {
          status, LinearOperatorStatusScope::rank_local, -1};
      return status;
    }
    y.unchecked(link.fluid_local_index, 0U) += correction;
  }
  return {};
}

Status IbmPressureOperator::mask_solid_rhs(FieldView rhs) const noexcept {
  if (implementation_ == nullptr ||
      !valid_scalar(rhs, implementation_->cells, 0U)) {
    return {StatusCode::invalid_plan, kIbmPressureApply};
  }
  const Span<const std::uint8_t> region = implementation_->topology->region();
  const std::size_t expected =
      static_cast<std::size_t>(implementation_->cells.x) *
      implementation_->cells.y * implementation_->cells.z;
  if (region.size != expected) {
    return {StatusCode::invalid_plan, kIbmPressureApply};
  }
  for (const Int3 cell : implementation_->solid_cells) {
    rhs.unchecked(cell, 0U) = 0.0;
  }
  return {};
}

PlanFingerprint IbmPressureOperator::fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->fingerprint;
}

}  // namespace hundun::v04
