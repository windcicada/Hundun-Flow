// SPDX-License-Identifier: Apache-2.0

#include "hundun/finite_volume/matrix_free_poisson.hpp"

#include "matrix_free_poisson_detail.hpp"
#include "hundun/linear/ghosted_vector.hpp"
#include "hundun/linear/ghosted_vector_halo.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace hundun::finite_volume {
namespace detail {
namespace {

thread_local PoissonTestOptions test_options;
thread_local PoissonApplyTrace apply_trace;

}  // namespace

void set_poisson_test_options(PoissonTestOptions options) noexcept {
  test_options = options;
}

void reset_poisson_apply_trace() noexcept { apply_trace = {}; }

PoissonApplyTrace poisson_apply_trace() noexcept { return apply_trace; }

void record_poisson_apply_phase(PoissonApplyPhase phase) {
  if (!test_options.observe_apply &&
      !test_options.injected_apply_failure.has_value()) {
    return;
  }
  if (apply_trace.count >= apply_trace.phases.size()) {
    throw runtime::Error("Poisson apply trace capacity exceeded");
  }
  apply_trace.phases[apply_trace.count++] = phase;
  if (test_options.injected_apply_failure == phase) {
    throw runtime::Error("injected Poisson apply phase failure");
  }
}

void record_poisson_cleanup_wait() noexcept {
  if (!test_options.observe_apply &&
      !test_options.injected_apply_failure.has_value()) {
    return;
  }
  if (apply_trace.count < apply_trace.phases.size()) {
    apply_trace.phases[apply_trace.count++] = PoissonApplyPhase::halo_wait;
  }
}

bool consume_force_revision_wrap() noexcept {
  const bool result = test_options.force_revision_wrap;
  test_options.force_revision_wrap = false;
  return result;
}

bool consume_injected_reentrant_apply() noexcept {
  const bool result = test_options.inject_reentrant_apply;
  test_options.inject_reentrant_apply = false;
  return result;
}

}  // namespace detail

namespace {

constexpr std::uint64_t kInitialRevision = 1U;

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

runtime::Real3 subtract(runtime::Real3 left,
                        runtime::Real3 right) noexcept {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

runtime::Real3 multiply(runtime::Real3 value, double factor) noexcept {
  return {value.x * factor, value.y * factor, value.z * factor};
}

void validate_reference_context(const execution::ExecutionContext& context,
                                const char* subject) {
  if (context.backend_identity() == 0U ||
      context.space() != execution::ExecutionSpace::host ||
      !context.supports(execution::ExecutionCapability::host_access) ||
      !context.supports(execution::ExecutionCapability::buffer_allocation) ||
      !context.supports(execution::ExecutionCapability::transfer)) {
    throw runtime::Error(std::string(subject) +
                         " requires a host-accessible reference context");
  }
}

template <class View>
const double* validate_const_view(View view,
                                  const execution::ExecutionContext& context,
                                  std::size_t expected_size,
                                  const char* subject) {
  if (view.scalar_format() != execution::ScalarFormat::float64 ||
      view.size() != expected_size ||
      view.backend_identity() != context.backend_identity() ||
      view.space() != execution::ExecutionSpace::host ||
      view.space() != context.space()) {
    throw runtime::Error(std::string(subject) +
                         " view metadata does not match");
  }
  return view.data();
}

double* validate_mutable_view(execution::VectorView<double> view,
                              const execution::ExecutionContext& context,
                              std::size_t expected_size,
                              const char* subject) {
  if (!view.writable() ||
      view.scalar_format() != execution::ScalarFormat::float64 ||
      view.size() != expected_size ||
      view.backend_identity() != context.backend_identity() ||
      view.space() != execution::ExecutionSpace::host ||
      view.space() != context.space()) {
    throw runtime::Error(std::string(subject) +
                         " view metadata does not match");
  }
  return view.data();
}

class OperationGuard final {
 public:
  explicit OperationGuard(std::atomic<bool>& active) : active_(&active) {
    bool expected = false;
    if (!active_->compare_exchange_strong(expected, true)) {
      throw runtime::Error("matrix-free Poisson object is already active");
    }
  }

  ~OperationGuard() noexcept { active_->store(false); }
  OperationGuard(const OperationGuard&) = delete;
  OperationGuard& operator=(const OperationGuard&) = delete;

 private:
  std::atomic<bool>* active_;
};

struct FaceStencil final {
  mesh::LocalFaceId face{};
  mesh::LocalCellId owner{};
  mesh::LocalCellId neighbour{};
  bool has_neighbour{};
  bool owner_owned{};
  bool neighbour_owned{};
  bool periodic{};
  bool reference{};
  double orthogonal_factor{};
  double owner_inverse_volume{};
  double neighbour_inverse_volume{};
  runtime::Real3 nonorthogonal_area{};
};

std::vector<double> copy_positive_coefficients(
    execution::VectorView<const double> coefficients,
    const execution::ExecutionContext& context, std::size_t expected_size) {
  const double* values = validate_const_view(
      coefficients, context, expected_size, "Poisson face coefficient");
  std::vector<double> result;
  try {
    result.resize(expected_size);
  } catch (const std::bad_alloc&) {
    throw runtime::Error("Poisson face coefficient allocation failed");
  } catch (const std::length_error&) {
    throw runtime::Error("Poisson face coefficient count is unsupported");
  }
  for (std::size_t face = 0; face < expected_size; ++face) {
    const double value = values[face * coefficients.stride()];
    if (!std::isfinite(value) || value <= 0.0) {
      throw runtime::Error(
          "Poisson face coefficients must be finite and positive");
    }
    result[face] = value;
  }
  return result;
}

void validate_boundary(const mesh::MeshTopology& topology,
                       PoissonBoundarySpec boundary) {
  switch (boundary.mode) {
    case PressureConstraintMode::constant_nullspace:
      if (boundary.pressure_reference_patch_id.has_value()) {
        throw runtime::Error(
            "constant-nullspace Poisson boundary cannot name a reference "
            "patch");
      }
      return;
    case PressureConstraintMode::pressure_reference_patch:
      if (!boundary.pressure_reference_patch_id.has_value() ||
          *boundary.pressure_reference_patch_id >= 6U) {
        throw runtime::Error(
            "pressure-reference Poisson boundary requires patch 0..5");
      }
      if (topology.patch(*boundary.pressure_reference_patch_id)
              .pairing_kind() == mesh::PatchPairingKind::periodic) {
        throw runtime::Error(
            "pressure-reference Poisson patch must be non-periodic");
      }
      return;
    default:
      throw runtime::Error("invalid Poisson pressure constraint mode");
  }
}

std::vector<FaceStencil> build_stencils(
    const mesh::MeshTopology& topology, const mesh::MeshGeometry& geometry,
    PoissonBoundarySpec boundary) {
  std::vector<FaceStencil> result;
  try {
    result.reserve(topology.local_face_count());
  } catch (const std::bad_alloc&) {
    throw runtime::Error("Poisson stencil allocation failed");
  } catch (const std::length_error&) {
    throw runtime::Error("Poisson stencil count is unsupported");
  }
  for (mesh::LocalFaceId face = 0; face < topology.local_face_count();
       ++face) {
    const mesh::LocalCellId owner = topology.owner(face);
    const std::optional<mesh::LocalCellId> neighbour =
        topology.neighbour(face);
    const bool owner_owned =
        topology.cell_ownership(owner) == mesh::EntityOwnership::owned;
    const bool neighbour_owned =
        neighbour.has_value() && topology.cell_ownership(*neighbour) ==
                                     mesh::EntityOwnership::owned;
    const bool periodic = topology.periodic_pair(face).has_value();
    const std::optional<std::uint32_t> patch = topology.patch_id(face);
    const bool reference =
        patch.has_value() &&
        boundary.mode == PressureConstraintMode::pressure_reference_patch &&
        boundary.pressure_reference_patch_id == patch;
    if ((!neighbour.has_value() && !reference) ||
        (periodic && !owner_owned) ||
        (!owner_owned && !neighbour_owned)) {
      continue;
    }
    const runtime::Real3 area =
        geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
    const runtime::Real3 displacement = geometry.face_displacement_m(face);
    const double area_squared = dot(area, area);
    const double projection = dot(area, displacement);
    const double factor = area_squared / projection;
    if (!finite(area) || !finite(displacement) ||
        !std::isfinite(area_squared) || area_squared <= 0.0 ||
        !std::isfinite(projection) || projection <= 0.0 ||
        !std::isfinite(factor) || factor <= 0.0) {
      throw runtime::Error("Poisson orthogonal face factor is invalid");
    }
    const runtime::Real3 nonorthogonal =
        subtract(area, multiply(displacement, factor));
    if (!finite(nonorthogonal)) {
      throw runtime::Error("Poisson non-orthogonal face area is invalid");
    }
    FaceStencil stencil;
    stencil.face = face;
    stencil.owner = owner;
    stencil.has_neighbour = neighbour.has_value();
    stencil.neighbour = neighbour.value_or(0U);
    stencil.owner_owned = owner_owned;
    stencil.neighbour_owned = neighbour_owned;
    stencil.periodic = periodic;
    stencil.reference = reference;
    stencil.orthogonal_factor = factor;
    stencil.nonorthogonal_area = nonorthogonal;
    if (owner_owned) {
      stencil.owner_inverse_volume = 1.0 / geometry.cell_volume_m3(owner);
    }
    if (neighbour_owned) {
      stencil.neighbour_inverse_volume =
          1.0 / geometry.cell_volume_m3(*neighbour);
    }
    if ((owner_owned &&
         (!std::isfinite(stencil.owner_inverse_volume) ||
          stencil.owner_inverse_volume <= 0.0)) ||
        (neighbour_owned &&
         (!std::isfinite(stencil.neighbour_inverse_volume) ||
          stencil.neighbour_inverse_volume <= 0.0))) {
      throw runtime::Error("Poisson cell inverse volume is invalid");
    }
    result.push_back(stencil);
  }
  return result;
}

bool stencil_needs_remote(const FaceStencil& stencil,
                          std::size_t owned_count) noexcept {
  if (!stencil.has_neighbour || stencil.reference) {
    return false;
  }
  if (stencil.periodic) {
    return stencil.neighbour >= owned_count;
  }
  return (stencil.owner_owned && stencil.neighbour >= owned_count) ||
         (stencil.neighbour_owned && stencil.owner >= owned_count);
}

void apply_stencil(const FaceStencil& stencil,
                   const std::vector<double>& gamma, const double* values,
                   double* output) {
  const double coefficient =
      gamma[stencil.face] * stencil.orthogonal_factor;
  if (stencil.reference) {
    output[stencil.owner] += coefficient * values[stencil.owner] *
                             stencil.owner_inverse_volume;
    return;
  }
  const double difference =
      values[stencil.owner] - values[stencil.neighbour];
  if (stencil.owner_owned) {
    output[stencil.owner] +=
        coefficient * difference * stencil.owner_inverse_volume;
  }
  if (!stencil.periodic && stencil.neighbour_owned) {
    output[stencil.neighbour] -=
        coefficient * difference * stencil.neighbour_inverse_volume;
  }
}

std::vector<double> build_diagonal(
    std::size_t owned_count, const std::vector<FaceStencil>& stencils,
    const std::vector<double>& gamma) {
  std::vector<double> diagonal;
  try {
    diagonal.assign(owned_count, 0.0);
  } catch (const std::bad_alloc&) {
    throw runtime::Error("Poisson diagonal allocation failed");
  } catch (const std::length_error&) {
    throw runtime::Error("Poisson diagonal count is unsupported");
  }
  for (const FaceStencil& stencil : stencils) {
    const double coefficient =
        gamma[stencil.face] * stencil.orthogonal_factor;
    if (!std::isfinite(coefficient) || coefficient <= 0.0) {
      throw runtime::Error("Poisson diagonal face coefficient is invalid");
    }
    const bool periodic_self =
        stencil.periodic && stencil.has_neighbour &&
        stencil.owner == stencil.neighbour;
    if (stencil.owner_owned && !periodic_self) {
      const double contribution =
          coefficient * stencil.owner_inverse_volume;
      if (!std::isfinite(contribution) || contribution <= 0.0 ||
          !std::isfinite(diagonal[stencil.owner] + contribution)) {
        throw runtime::Error("Poisson owner diagonal is invalid");
      }
      diagonal[stencil.owner] += contribution;
    }
    if (!stencil.periodic && !stencil.reference &&
        stencil.neighbour_owned) {
      const double contribution =
          coefficient * stencil.neighbour_inverse_volume;
      if (!std::isfinite(contribution) || contribution <= 0.0 ||
          !std::isfinite(diagonal[stencil.neighbour] + contribution)) {
        throw runtime::Error("Poisson neighbour diagonal is invalid");
      }
      diagonal[stencil.neighbour] += contribution;
    }
  }
  for (const double value : diagonal) {
    if (!std::isfinite(value) || value < 0.0) {
      throw runtime::Error("Poisson diagonal must be finite and nonnegative");
    }
  }
  return diagonal;
}

}  // namespace

struct MatrixFreePoissonOperator::Impl final {
  Impl(const mesh::MeshTopology& topology,
       mesh::MappingKind mapping_kind,
       execution::ExecutionContext& execution_context,
       PoissonBoundarySpec supplied_boundary,
       std::vector<FaceStencil> supplied_stencils,
       std::vector<double> supplied_gamma,
       std::vector<double> supplied_diagonal)
      : context(&execution_context),
        layout(linear::VectorLayout::from_topology(topology)),
        workspace(execution_context, layout),
        output_scratch(topology.owned_cell_count(), 0.0),
        boundary(supplied_boundary),
        solver(mapping_kind == mesh::MappingKind::uniform_box
                   ? PoissonSolverFamily::conjugate_gradient
                   : PoissonSolverFamily::bicgstab),
        stencils(std::move(supplied_stencils)),
        gamma(std::move(supplied_gamma)),
        diagonal_values(std::move(supplied_diagonal)),
        completed(execution::ExecutionEvent::completed()) {}

  execution::ExecutionContext* context;
  linear::VectorLayout layout;
  linear::GhostedVector workspace;
  std::optional<linear::GhostedVectorHalo> halo;
  mutable std::vector<double> output_scratch;
  PoissonBoundarySpec boundary;
  PoissonSolverFamily solver;
  std::vector<FaceStencil> stencils;
  std::vector<double> gamma;
  std::vector<double> diagonal_values;
  std::uint64_t revision{kInitialRevision};
  execution::ExecutionEvent completed;
  mutable std::atomic<bool> active{false};
};

MatrixFreePoissonOperator MatrixFreePoissonOperator::create(
    const runtime::StructuredDecomposition& decomposition,
    const mesh::MeshTopology& topology, const mesh::MeshGeometry& geometry,
    execution::ExecutionContext& execution_context,
    execution::VectorView<const double> gamma_by_local_face,
    PoissonBoundarySpec boundary) {
  runtime::MpiContext preflight =
      runtime::MpiContext::duplicate(decomposition.comm());
  std::unique_ptr<Impl> implementation;
  bool local_ok = true;
  std::string_view local_message;
  try {
    validate_reference_context(execution_context, "matrix-free Poisson");
    geometry.require_compatible(topology);
    if (decomposition.global_extent().x != topology.global_extent().x ||
        decomposition.global_extent().y != topology.global_extent().y ||
        decomposition.global_extent().z != topology.global_extent().z ||
        decomposition.owned_box().begin.x !=
            topology.owned_global_box().begin.x ||
        decomposition.owned_box().begin.y !=
            topology.owned_global_box().begin.y ||
        decomposition.owned_box().begin.z !=
            topology.owned_global_box().begin.z ||
        decomposition.owned_box().end.x != topology.owned_global_box().end.x ||
        decomposition.owned_box().end.y != topology.owned_global_box().end.y ||
        decomposition.owned_box().end.z != topology.owned_global_box().end.z) {
      throw runtime::Error(
          "matrix-free Poisson decomposition and topology do not match");
    }
    validate_boundary(topology, boundary);
    std::vector<double> gamma = copy_positive_coefficients(
        gamma_by_local_face, execution_context, topology.local_face_count());
    std::vector<FaceStencil> stencils =
        build_stencils(topology, geometry, boundary);
    std::vector<double> diagonal =
        build_diagonal(topology.owned_cell_count(), stencils, gamma);
    implementation = std::make_unique<Impl>(
        topology, geometry.mapping_kind(), execution_context, boundary,
        std::move(stencils), std::move(gamma), std::move(diagonal));
  } catch (const runtime::Error& error) {
    local_ok = false;
    static_cast<void>(error);
    local_message = "local Poisson validation failed";
  } catch (const std::bad_alloc&) {
    local_ok = false;
    local_message = "matrix-free Poisson allocation failed";
  } catch (const std::length_error&) {
    local_ok = false;
    local_message = "matrix-free Poisson allocation size is unsupported";
  }
  const runtime::CollectiveStatus status = runtime::collective_status(
      preflight, local_ok, local_ok ? std::string_view{} : local_message);
  if (!status.ok || !implementation) {
    throw runtime::Error("matrix-free Poisson collective construction "
                         "failure: rank=" +
                         std::to_string(status.failing_rank) + " " +
                         status.message);
  }
  implementation->halo.emplace(linear::GhostedVectorHalo::create(
      decomposition, topology, execution_context));
  return MatrixFreePoissonOperator(std::move(implementation));
}

MatrixFreePoissonOperator::MatrixFreePoissonOperator(
    std::unique_ptr<Impl> implementation) noexcept
    : impl_(std::move(implementation)) {}

MatrixFreePoissonOperator::~MatrixFreePoissonOperator() noexcept = default;
MatrixFreePoissonOperator::MatrixFreePoissonOperator(
    MatrixFreePoissonOperator&&) noexcept = default;

linear::VectorLayout MatrixFreePoissonOperator::domain_layout() const {
  return impl_->layout;
}

linear::VectorLayout MatrixFreePoissonOperator::range_layout() const {
  return impl_->layout;
}

const execution::ExecutionContext& MatrixFreePoissonOperator::context() const {
  return *impl_->context;
}

std::uint64_t MatrixFreePoissonOperator::revision() const {
  return impl_->revision;
}

execution::ExecutionEvent MatrixFreePoissonOperator::apply(
    execution::VectorView<const double> x,
    execution::VectorView<double> y) const {
  OperationGuard guard(impl_->active);
  if (detail::consume_injected_reentrant_apply()) {
    static_cast<void>(apply(x, y));
  }
  const std::size_t owned_count = impl_->layout.owned_count();
  const double* input =
      validate_const_view(x, *impl_->context, owned_count, "Poisson input");
  double* output =
      validate_mutable_view(y, *impl_->context, owned_count, "Poisson output");
  auto local_view = impl_->workspace.local_view();
  double* local_values = local_view.data();
  for (std::size_t cell = 0; cell < owned_count; ++cell) {
    const double value = input[cell * x.stride()];
    if (!std::isfinite(value)) {
      throw runtime::Error("Poisson input must be finite");
    }
    local_values[cell] = value;
  }
  bool halo_started = false;
  try {
    impl_->halo->begin(impl_->workspace);
    halo_started = true;
    detail::record_poisson_apply_phase(
        detail::PoissonApplyPhase::halo_begin);
    std::fill(impl_->output_scratch.begin(), impl_->output_scratch.end(), 0.0);
    for (const FaceStencil& stencil : impl_->stencils) {
      if (!stencil_needs_remote(stencil, owned_count)) {
        apply_stencil(stencil, impl_->gamma, local_values,
                      impl_->output_scratch.data());
      }
    }
    detail::record_poisson_apply_phase(
        detail::PoissonApplyPhase::local_rows);
    impl_->halo->wait(impl_->workspace);
    halo_started = false;
    detail::record_poisson_apply_phase(detail::PoissonApplyPhase::halo_wait);
    for (const FaceStencil& stencil : impl_->stencils) {
      if (stencil_needs_remote(stencil, owned_count)) {
        apply_stencil(stencil, impl_->gamma, local_values,
                      impl_->output_scratch.data());
      }
    }
    detail::record_poisson_apply_phase(
        detail::PoissonApplyPhase::partition_rows);
  } catch (...) {
    const std::exception_ptr failure = std::current_exception();
    if (halo_started) {
      try {
        impl_->halo->wait(impl_->workspace);
        detail::record_poisson_cleanup_wait();
      } catch (...) {
      }
    }
    std::rethrow_exception(failure);
  }

  for (const double value : impl_->output_scratch) {
    if (!std::isfinite(value)) {
      throw runtime::Error("Poisson operator produced a non-finite result");
    }
  }
  for (std::size_t cell = 0; cell < owned_count; ++cell) {
    output[cell * y.stride()] = impl_->output_scratch[cell];
  }
  return impl_->completed;
}

bool MatrixFreePoissonOperator::has_diagonal() const { return true; }

execution::ExecutionEvent MatrixFreePoissonOperator::diagonal(
    execution::VectorView<double> d) const {
  OperationGuard guard(impl_->active);
  double* output = validate_mutable_view(
      d, *impl_->context, impl_->layout.owned_count(), "Poisson diagonal");
  for (std::size_t cell = 0; cell < impl_->diagonal_values.size(); ++cell) {
    output[cell * d.stride()] = impl_->diagonal_values[cell];
  }
  return impl_->completed;
}

PressureConstraintMode MatrixFreePoissonOperator::constraint_mode() const
    noexcept {
  return impl_->boundary.mode;
}

std::optional<std::uint32_t>
MatrixFreePoissonOperator::pressure_reference_patch_id() const noexcept {
  return impl_->boundary.pressure_reference_patch_id;
}

PoissonSolverFamily MatrixFreePoissonOperator::solver_family() const noexcept {
  return impl_->solver;
}

std::uint64_t MatrixFreePoissonOperator::replace_face_coefficients(
    execution::VectorView<const double> gamma_by_local_face) {
  OperationGuard guard(impl_->active);
  std::vector<double> candidate = copy_positive_coefficients(
      gamma_by_local_face, *impl_->context, impl_->gamma.size());
  std::vector<double> candidate_diagonal = build_diagonal(
      impl_->layout.owned_count(), impl_->stencils, candidate);
  if (detail::consume_force_revision_wrap() ||
      impl_->revision == std::numeric_limits<std::uint64_t>::max()) {
    throw runtime::Error("matrix-free Poisson revision would wrap");
  }
  impl_->gamma = std::move(candidate);
  impl_->diagonal_values = std::move(candidate_diagonal);
  ++impl_->revision;
  return impl_->revision;
}

execution::ExecutionEvent
MatrixFreePoissonOperator::accumulate_explicit_nonorthogonal_rhs(
    execution::VectorView<const double> gradient_x_by_local_face,
    execution::VectorView<const double> gradient_y_by_local_face,
    execution::VectorView<const double> gradient_z_by_local_face,
    execution::VectorView<double> rhs) const {
  OperationGuard guard(impl_->active);
  const std::size_t face_count = impl_->gamma.size();
  const double* gradient_x = validate_const_view(
      gradient_x_by_local_face, *impl_->context, face_count,
      "Poisson x face gradient");
  const double* gradient_y = validate_const_view(
      gradient_y_by_local_face, *impl_->context, face_count,
      "Poisson y face gradient");
  const double* gradient_z = validate_const_view(
      gradient_z_by_local_face, *impl_->context, face_count,
      "Poisson z face gradient");
  double* output = validate_mutable_view(
      rhs, *impl_->context, impl_->layout.owned_count(), "Poisson RHS");
  for (std::size_t cell = 0; cell < impl_->layout.owned_count(); ++cell) {
    const double value = output[cell * rhs.stride()];
    if (!std::isfinite(value)) {
      throw runtime::Error("Poisson RHS must be finite");
    }
    impl_->output_scratch[cell] = value;
  }
  for (const FaceStencil& stencil : impl_->stencils) {
    const std::size_t face = stencil.face;
    const runtime::Real3 gradient{
        gradient_x[face * gradient_x_by_local_face.stride()],
        gradient_y[face * gradient_y_by_local_face.stride()],
        gradient_z[face * gradient_z_by_local_face.stride()]};
    if (!finite(gradient)) {
      throw runtime::Error("Poisson face gradient must be finite");
    }
    const double correction =
        impl_->gamma[face] * dot(gradient, stencil.nonorthogonal_area);
    if (!std::isfinite(correction)) {
      throw runtime::Error("Poisson explicit correction is non-finite");
    }
    if (stencil.owner_owned) {
      impl_->output_scratch[stencil.owner] +=
          correction * stencil.owner_inverse_volume;
    }
    if (!stencil.periodic && !stencil.reference &&
        stencil.neighbour_owned) {
      impl_->output_scratch[stencil.neighbour] -=
          correction * stencil.neighbour_inverse_volume;
    }
  }
  for (const double value : impl_->output_scratch) {
    if (!std::isfinite(value)) {
      throw runtime::Error("Poisson explicit RHS candidate is non-finite");
    }
  }
  for (std::size_t cell = 0; cell < impl_->layout.owned_count(); ++cell) {
    output[cell * rhs.stride()] = impl_->output_scratch[cell];
  }
  return impl_->completed;
}

struct PoissonConstraint::Impl final {
  Impl(execution::ExecutionContext& execution_context,
       const runtime::MpiContext& mpi_context,
       PressureConstraintMode supplied_mode,
       linear::VectorLayout supplied_layout,
       std::vector<double> supplied_volumes)
      : context(&execution_context), mpi(&mpi_context), mode(supplied_mode),
        layout(std::move(supplied_layout)),
        volumes(std::move(supplied_volumes)) {}

  execution::ExecutionContext* context;
  const runtime::MpiContext* mpi;
  PressureConstraintMode mode;
  linear::VectorLayout layout;
  std::vector<double> volumes;
};

PoissonConstraint PoissonConstraint::create(
    const mesh::MeshTopology& topology, const mesh::MeshGeometry& geometry,
    execution::ExecutionContext& execution_context,
    const runtime::MpiContext& mpi_context, PressureConstraintMode mode) {
  validate_reference_context(execution_context, "Poisson constraint");
  geometry.require_compatible(topology);
  if (mode != PressureConstraintMode::constant_nullspace &&
      mode != PressureConstraintMode::pressure_reference_patch) {
    throw runtime::Error("invalid Poisson constraint mode");
  }
  std::vector<double> volumes;
  try {
    volumes.resize(topology.owned_cell_count());
    for (mesh::LocalCellId cell = 0; cell < topology.owned_cell_count();
         ++cell) {
      const double value = geometry.cell_volume_m3(cell);
      if (!std::isfinite(value) || value <= 0.0) {
        throw runtime::Error("Poisson constraint volume is invalid");
      }
      volumes[cell] = value;
    }
    return PoissonConstraint(std::make_unique<Impl>(
        execution_context, mpi_context, mode,
        linear::VectorLayout::from_topology(topology), std::move(volumes)));
  } catch (const runtime::Error&) {
    throw;
  } catch (const std::bad_alloc&) {
    throw runtime::Error("Poisson constraint allocation failed");
  } catch (const std::length_error&) {
    throw runtime::Error("Poisson constraint size is unsupported");
  }
}

PoissonConstraint::PoissonConstraint(std::unique_ptr<Impl> implementation)
    noexcept
    : impl_(std::move(implementation)) {}

PoissonConstraint::~PoissonConstraint() noexcept = default;
PoissonConstraint::PoissonConstraint(PoissonConstraint&&) noexcept = default;

PressureConstraintMode PoissonConstraint::mode() const noexcept {
  return impl_->mode;
}

namespace {

template <class ConstraintImplementation>
void apply_constraint_operation(ConstraintImplementation& impl,
                                execution::VectorView<double> values,
                                const char* subject) {
  if (impl.mode == PressureConstraintMode::pressure_reference_patch) {
    return;
  }
  double* data = nullptr;
  double local_weighted_sum = 0.0;
  double local_volume_sum = 0.0;
  bool local_ok = true;
  try {
    data = validate_mutable_view(values, *impl.context,
                                 impl.layout.owned_count(), subject);
    for (std::size_t cell = 0; cell < impl.volumes.size(); ++cell) {
      const double value = data[cell * values.stride()];
      const double weighted = impl.volumes[cell] * value;
      if (!std::isfinite(value) || !std::isfinite(weighted) ||
          !std::isfinite(local_weighted_sum + weighted) ||
          !std::isfinite(local_volume_sum + impl.volumes[cell])) {
        local_ok = false;
        break;
      }
      local_weighted_sum += weighted;
      local_volume_sum += impl.volumes[cell];
    }
  } catch (const runtime::Error&) {
    local_ok = false;
  }
  double preflight = local_ok ? 0.0 : 1.0;
  impl.mpi->allreduce_fp64_in_place(
      &preflight, 1U, runtime::Fp64ReductionOperation::maximum);
  if (preflight != 0.0 || data == nullptr) {
    throw runtime::Error("Poisson constraint collective preflight failed");
  }
  double sums[2]{local_weighted_sum, local_volume_sum};
  impl.mpi->allreduce_fp64_in_place(sums, 2U,
                                    runtime::Fp64ReductionOperation::sum);
  const double mean = sums[0] / sums[1];
  if (!std::isfinite(sums[0]) || !std::isfinite(sums[1]) || sums[1] <= 0.0 ||
      !std::isfinite(mean)) {
    throw runtime::Error("Poisson constraint global mean is invalid");
  }
  for (std::size_t cell = 0; cell < impl.volumes.size(); ++cell) {
    const double candidate = data[cell * values.stride()] - mean;
    if (!std::isfinite(candidate)) {
      throw runtime::Error("Poisson constraint candidate is non-finite");
    }
  }
  for (std::size_t cell = 0; cell < impl.volumes.size(); ++cell) {
    data[cell * values.stride()] -= mean;
  }
}

}  // namespace

void PoissonConstraint::project_rhs(
    execution::VectorView<double> rhs) const {
  apply_constraint_operation(*impl_, rhs, "Poisson RHS projection");
}

void PoissonConstraint::normalize_solution(
    execution::VectorView<double> x) const {
  apply_constraint_operation(*impl_, x, "Poisson solution normalization");
}

}  // namespace hundun::finite_volume
