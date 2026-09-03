// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_ibm.hpp"

#include "solver_ibm_force_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kForceInput = 9101U;
constexpr std::uint32_t kForceReconstruction = 9102U;
constexpr std::uint32_t kForceCollective = 9103U;
constexpr std::uint32_t kForceCache = 9104U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

bool valid_gradient(const VelocityGradient& gradient) noexcept {
  return std::all_of(gradient.value.begin(), gradient.value.end(),
                     [](double value) { return std::isfinite(value); });
}

Status add_point(const SurfaceTractionPoint& point, Real3 origin,
                 SurfaceForce& force) noexcept {
  const double normal_squared =
      point.solid_to_fluid_normal.x * point.solid_to_fluid_normal.x +
      point.solid_to_fluid_normal.y * point.solid_to_fluid_normal.y +
      point.solid_to_fluid_normal.z * point.solid_to_fluid_normal.z;
  if (!detail::finite_force_vector(point.position) ||
      !detail::finite_force_vector(point.solid_to_fluid_normal) ||
      !std::isfinite(point.weight) || !(point.weight > 0.0) ||
      !std::isfinite(point.absolute_pressure) ||
      !std::isfinite(point.effective_viscosity) ||
      !(point.effective_viscosity > 0.0) ||
      !valid_gradient(point.velocity_gradient) ||
      std::abs(normal_squared - 1.0) > 1.0e-10) {
    return {StatusCode::numerical_failure, kForceInput};
  }
  const Real3 pressure = detail::force_scale(
      -point.absolute_pressure, point.solid_to_fluid_normal);
  const auto& g = point.velocity_gradient.value;
  const double divergence = g[0U] + g[4U] + g[8U];
  double tau[3][3]{};
  for (std::size_t i = 0U; i < 3U; ++i) {
    for (std::size_t j = 0U; j < 3U; ++j) {
      tau[i][j] = point.effective_viscosity *
                  (g[i * 3U + j] + g[j * 3U + i] -
                   (i == j ? (2.0 / 3.0) * divergence : 0.0));
    }
  }
  Real3 viscous{
      tau[0][0] * point.solid_to_fluid_normal.x +
          tau[0][1] * point.solid_to_fluid_normal.y +
          tau[0][2] * point.solid_to_fluid_normal.z,
      tau[1][0] * point.solid_to_fluid_normal.x +
          tau[1][1] * point.solid_to_fluid_normal.y +
          tau[1][2] * point.solid_to_fluid_normal.z,
      tau[2][0] * point.solid_to_fluid_normal.x +
          tau[2][1] * point.solid_to_fluid_normal.y +
          tau[2][2] * point.solid_to_fluid_normal.z};
  const Real3 total = detail::force_add(pressure, viscous);
  const Real3 arm{point.position.x - origin.x,
                  point.position.y - origin.y,
                  point.position.z - origin.z};
  force.pressure = detail::force_add(
      force.pressure, detail::force_scale(point.weight, pressure));
  force.viscous = detail::force_add(
      force.viscous, detail::force_scale(point.weight, viscous));
  force.moment = detail::force_add(
      force.moment,
      detail::force_scale(point.weight, detail::force_cross(arm, total)));
  return {};
}

bool valid_surface_state(const SurfaceQuadraturePlan& quadrature,
                         const FinalSurfaceState& state) noexcept {
  const Int3 cells = state.pressure_perturbation.interior;
  const std::uint8_t reach = quadrature.reconstruction().maximum_halo_reach();
  return quadrature.physical_fingerprint() != 0U &&
         quadrature.local_layout_fingerprint() != 0U &&
         state.terminal_plan != 0U && state.terminal_state != 0U &&
         state.final_flux != 0U && state.face_flux.revision == state.final_flux &&
         state.face_flux.certificate.valid() &&
         state.face_flux.certificate.matches(state.face_flux) &&
         detail::valid_force_field(state.final_velocity, cells, 3U, 0U) &&
         detail::valid_force_field(state.pressure_perturbation, cells, 1U,
                                   reach) &&
         detail::valid_force_field(state.velocity_gradient, cells, 9U,
                                   reach) &&
         detail::valid_force_field(state.effective_viscosity, cells, 1U,
                                   reach) &&
         state.gradient_authority.velocity == state.final_velocity.revision &&
         state.gradient_authority.velocity_storage ==
             state.final_velocity.storage_identity &&
         state.gradient_authority.velocity_revision_domain ==
             state.final_velocity.revision_domain &&
         state.gradient_authority.geometry == state.geometry &&
         state.gradient_authority.geometry_plan_identity != 0U &&
         state.gradient_authority.patch_identity != 0U &&
         state.turbulence.valid() &&
         state.turbulence.gradient == state.velocity_gradient.revision &&
         state.turbulence.effective_viscosity ==
             state.effective_viscosity.revision &&
         state.gradient_authority.turbulence_plan_identity ==
             state.turbulence.plan &&
         state.geometry != 0U && std::isfinite(state.pressure_reference) &&
         detail::finite_force_vector(state.moment_origin);
}

RevisionToken force_state_revision(const SurfaceQuadraturePlan& quadrature,
                                   const FinalSurfaceState& state) noexcept {
  std::uint64_t hash = kFnvOffset;
  hash = mix(hash, quadrature.physical_fingerprint());
  hash = mix(hash, state.terminal_plan);
  hash = mix(hash, state.terminal_state);
  hash = mix(hash, state.final_flux);
  hash = mix(hash, state.final_velocity.revision);
  hash = mix(hash, state.pressure_perturbation.revision);
  hash = mix(hash, state.velocity_gradient.revision);
  hash = mix(hash, state.effective_viscosity.revision);
  hash = mix(hash, state.geometry);
  return hash == 0U ? 1U : hash;
}

Status force_consensus(MPI_Comm communicator, Status local) noexcept {
  int rank = -1;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kForceCollective};
  }
  const int candidate = local ? size : rank;
  int selected = size;
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kForceCollective};
  }
  if (selected == size) return {};
  std::uint64_t wire = 0U;
  if (rank == selected) {
    wire = (static_cast<std::uint64_t>(local.detail) << 16U) |
           static_cast<std::uint16_t>(local.code);
  }
  if (MPI_Bcast(&wire, 1, MPI_UINT64_T, selected, communicator) !=
      MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kForceCollective};
  }
  return {static_cast<StatusCode>(wire & UINT64_C(0xffff)),
          static_cast<std::uint32_t>(wire >> 16U)};
}

}  // namespace

Status integrate_surface_traction(Span<const SurfaceTractionPoint> points,
                                  Real3 moment_origin,
                                  SurfaceForce& out) noexcept {
  if (points.data == nullptr || points.size == 0U ||
      !detail::finite_force_vector(moment_origin)) {
    return {StatusCode::invalid_plan, kForceInput};
  }
  SurfaceForce candidate;
  for (std::size_t index = 0U; index < points.size; ++index) {
    const Status status = add_point(points.data[index], moment_origin,
                                    candidate);
    if (!status) {
      return status;
    }
  }
  candidate.total = detail::force_add(candidate.pressure, candidate.viscous);
  if (!detail::finite_force_vector(candidate.pressure) ||
      !detail::finite_force_vector(candidate.viscous) ||
      !detail::finite_force_vector(candidate.total) ||
      !detail::finite_force_vector(candidate.moment)) {
    return {StatusCode::numerical_failure, kForceInput};
  }
  candidate.revision = 1U;
  out = candidate;
  return {};
}

Status evaluate_surface_force(MPI_Comm communicator,
                              const SurfaceQuadraturePlan& quadrature,
                              const FinalSurfaceState& state,
                              SurfaceForce& out) noexcept {
  if (communicator == MPI_COMM_NULL)
    return {StatusCode::invalid_plan, kForceInput};
  Status local_status = valid_surface_state(quadrature, state)
                            ? Status{}
                            : Status{StatusCode::invalid_plan, kForceInput};
  Status agreed = force_consensus(communicator, local_status);
  if (!agreed) return agreed;
  SurfaceForce local_force;
  const Span<const SurfaceQuadraturePoint> points = quadrature.local_points();
  const QuadraticStencilPlan& reconstruction = quadrature.reconstruction();
  for (std::size_t index = 0U; index < points.size && local_status; ++index) {
    const SurfaceQuadraturePoint& quadrature_point = points.data[index];
    SurfaceTractionPoint point;
    point.position = quadrature_point.position;
    point.solid_to_fluid_normal =
        quadrature_point.solid_to_fluid_normal;
    point.weight = quadrature_point.weight;
    double perturbation = 0.0;
    Status status = evaluate_quadratic_row(
        reconstruction, quadrature_point.wall_value_row,
        state.pressure_perturbation, 0U, 0.0, 0.0, perturbation);
    point.absolute_pressure = state.pressure_reference + perturbation;
    for (std::uint8_t component = 0U; component < 9U && status; ++component) {
      status = evaluate_quadratic_row(
          reconstruction, quadrature_point.wall_value_row,
          state.velocity_gradient, component, 0.0, 0.0,
          point.velocity_gradient.value[component]);
    }
    if (status) {
      status = evaluate_positive_bounded_quadratic_row(
          reconstruction, quadrature_point.wall_value_row,
          state.effective_viscosity, 0U, point.effective_viscosity);
    }
    if (status) {
      status = add_point(point, state.moment_origin, local_force);
    }
    if (!status) {
      local_status = {status.code, kForceReconstruction};
    }
  }
  agreed = force_consensus(communicator, local_status);
  if (!agreed) return agreed;
  std::array<double, 9U> local_values{
      local_force.pressure.x, local_force.pressure.y,
      local_force.pressure.z, local_force.viscous.x,
      local_force.viscous.y,  local_force.viscous.z,
      local_force.moment.x,   local_force.moment.y,
      local_force.moment.z};
  std::array<double, 9U> global{};
  if (MPI_Allreduce(local_values.data(), global.data(),
                    static_cast<int>(global.size()), MPI_DOUBLE, MPI_SUM,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kForceCollective};
  }
  SurfaceForce candidate;
  candidate.pressure = {global[0U], global[1U], global[2U]};
  candidate.viscous = {global[3U], global[4U], global[5U]};
  candidate.total = detail::force_add(candidate.pressure, candidate.viscous);
  candidate.moment = {global[6U], global[7U], global[8U]};
  candidate.revision = force_state_revision(quadrature, state);
  if (!detail::finite_force_vector(candidate.pressure) ||
      !detail::finite_force_vector(candidate.viscous) ||
      !detail::finite_force_vector(candidate.total) ||
      !detail::finite_force_vector(candidate.moment)) {
    return {StatusCode::numerical_failure, kForceCollective};
  }
  out = candidate;
  return {};
}

struct FinalForceCache::Impl {
  MPI_Comm communicator{MPI_COMM_NULL};
  const SurfaceQuadraturePlan* quadrature{};
  SurfaceForce active_force{};
  SurfaceForce pending_force{};
  FinalForceCertificate active{};
  FinalForceCertificate pending{};
  RevisionToken geometry{};
  RevisionToken next_revision{1U};
  StageId stage{};
  RevisionSlotId slot{};
  PlanFingerprint fingerprint{};
  std::uint64_t pending_attempt{};
};

FinalForceCache::~FinalForceCache() noexcept { release(); }

FinalForceCache::FinalForceCache(FinalForceCache&& other) noexcept
    : implementation_(std::exchange(other.implementation_, nullptr)) {}

FinalForceCache& FinalForceCache::operator=(FinalForceCache&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = std::exchange(other.implementation_, nullptr);
  }
  return *this;
}

void FinalForceCache::release() noexcept {
  Impl* implementation = std::exchange(implementation_, nullptr);
  if (implementation != nullptr) {
    if (implementation->communicator != MPI_COMM_NULL) {
      MPI_Comm_free(&implementation->communicator);
    }
    delete implementation;
  }
}

Status FinalForceCache::bind(MPI_Comm communicator,
                             const SurfaceQuadraturePlan& quadrature,
                             RevisionToken geometry_revision, StageId stage,
                             RevisionSlotId cache_slot,
                             FinalForceCache& out) noexcept {
  if (communicator == MPI_COMM_NULL ||
      quadrature.physical_fingerprint() == 0U ||
      quadrature.local_layout_fingerprint() == 0U ||
      geometry_revision == 0U || stage == 0U) {
    return {StatusCode::invalid_plan, kForceCache};
  }
  Impl* candidate = new (std::nothrow) Impl;
  if (candidate == nullptr) {
    return {StatusCode::allocation_failure, kForceCache};
  }
  if (MPI_Comm_dup(communicator, &candidate->communicator) != MPI_SUCCESS ||
      MPI_Comm_set_errhandler(candidate->communicator, MPI_ERRORS_RETURN) !=
          MPI_SUCCESS) {
    if (candidate->communicator != MPI_COMM_NULL) {
      MPI_Comm_free(&candidate->communicator);
    }
    delete candidate;
    return {StatusCode::mpi_failure, kForceCollective};
  }
  candidate->quadrature = &quadrature;
  candidate->geometry = geometry_revision;
  candidate->stage = stage;
  candidate->slot = cache_slot;
  std::uint64_t hash = kFnvOffset;
  hash = mix(hash, quadrature.physical_fingerprint());
  hash = mix(hash, geometry_revision);
  hash = mix(hash, stage);
  hash = mix(hash, cache_slot);
  candidate->fingerprint = hash == 0U ? 1U : hash;
  out.release();
  out.implementation_ = candidate;
  return {};
}

Status FinalForceCache::prepare(
    const FinalSurfaceState& state,
    Span<const RevisionDependency> dependencies,
    AttemptTransaction& transaction,
    FinalForceCertificate& certificate) noexcept {
  if (implementation_ == nullptr || !transaction.active() ||
      transaction.attempt_identity() == 0U ||
      state.geometry != implementation_->geometry ||
      implementation_->pending.valid() ||
      implementation_->next_revision == 0U) {
    return {StatusCode::invalid_plan, kForceCache};
  }
  SurfaceForce force;
  const Status evaluated = evaluate_surface_force(
      implementation_->communicator, *implementation_->quadrature, state,
      force);
  if (!evaluated) {
    return evaluated;
  }
  const RevisionToken cache_revision = implementation_->next_revision;
  const Status published = transaction.publish_pending_cache(
      implementation_->slot, dependencies, {cache_revision});
  if (!published) {
    return published;
  }
  force.revision = cache_revision;
  FinalForceCertificate next{implementation_->fingerprint,
                             state.terminal_state,
                             state.final_flux,
                             cache_revision,
                             force_state_revision(
                                 *implementation_->quadrature, state)};
  implementation_->pending_force = force;
  implementation_->pending = next;
  implementation_->pending_attempt = transaction.attempt_identity();
  if (cache_revision == std::numeric_limits<RevisionToken>::max()) {
    implementation_->next_revision = 0U;
  } else {
    ++implementation_->next_revision;
  }
  certificate = next;
  return {};
}

Status FinalForceCache::finalize(
    const AttemptTransaction& transaction) noexcept {
  if (implementation_ == nullptr || !transaction.finished() ||
      !implementation_->pending.valid() ||
      implementation_->pending_attempt != transaction.attempt_identity()) {
    return {StatusCode::invalid_plan, kForceCache};
  }
  if (transaction.committed() &&
      transaction.pending_cache(implementation_->slot) ==
          implementation_->pending.force) {
    implementation_->active_force = implementation_->pending_force;
    implementation_->active = implementation_->pending;
  }
  implementation_->pending_force = {};
  implementation_->pending = {};
  implementation_->pending_attempt = 0U;
  return {};
}

Status FinalForceCache::committed(
    SurfaceForce& force, FinalForceCertificate& certificate) const noexcept {
  if (implementation_ == nullptr || !implementation_->active.valid()) {
    return {StatusCode::invalid_plan, kForceCache};
  }
  force = implementation_->active_force;
  certificate = implementation_->active;
  return {};
}

PlanFingerprint FinalForceCache::fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->fingerprint;
}

}  // namespace hundun::v04
