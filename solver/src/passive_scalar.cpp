// SPDX-License-Identifier: Apache-2.0

#include "hundun/solver/passive_scalar.hpp"

#include "hundun/mesh/uniform_structured_mesh.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/field_view.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "mpi_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::solver {
namespace {

constexpr std::string_view kMpiFormattingFallback =
    "passive scalar MPI result formatting failed";
constexpr std::string_view kHaloCompatibilityFallback =
    "passive scalar halo compatibility validation failed";
constexpr std::string_view kHaloExchangeFallback =
    "passive scalar halo exchange failed";
constexpr std::string_view kStageExchangeFallback =
    "passive scalar stage exchange failed";

bool same(runtime::Int3 left, runtime::Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(runtime::Box3 left, runtime::Box3 right) noexcept {
  return same(left.begin, right.begin) && same(left.end, right.end);
}

void converge(const runtime::MpiContext &context, bool local_ok,
              std::string_view local_message) {
  const runtime::CollectiveStatus status =
      runtime::collective_status(context, local_ok, local_message);
  if (!status.ok) {
    throw runtime::Error(status.message.empty()
                             ? std::string("passive scalar operation failed")
                             : status.message);
  }
}

void converge_mpi_result(const runtime::MpiContext &context, int result,
                         std::string_view operation) {
  bool local_ok = true;
  std::string retained_message;
  std::string_view local_message;
  try {
    runtime::detail::check_mpi(result, operation);
  } catch (const runtime::Error &error) {
    local_ok = false;
    try {
      retained_message = error.what();
      local_message = retained_message;
    } catch (...) {
      local_message = kMpiFormattingFallback;
    }
  } catch (...) {
    local_ok = false;
    local_message = kMpiFormattingFallback;
  }

  const runtime::CollectiveStatus status =
      runtime::collective_status(context, local_ok, local_message);
  if (!status.ok) {
    throw runtime::Error(status.message);
  }
}

void retain_project_error(const runtime::Error &error,
                          std::string &retained_message,
                          std::string_view &selected_message,
                          std::string_view fallback) noexcept {
  selected_message = fallback;
  try {
    const char *message = error.what();
    if (message != nullptr && message[0] != '\0') {
      retained_message = message;
      if (!retained_message.empty()) {
        selected_message = retained_message;
      }
    }
  } catch (...) {
  }
}

std::uint64_t bit_pattern(double value) noexcept {
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

void require_agreed_value(const runtime::MpiContext &context, double value,
                          bool locally_valid,
                          std::string_view failure_message) {
  const std::uint64_t local_bits = bit_pattern(value);
  std::uint64_t minimum_bits = 0;
  std::uint64_t maximum_bits = 0;
  converge_mpi_result(context,
                      MPI_Allreduce(&local_bits, &minimum_bits, 1, MPI_UINT64_T,
                                    MPI_MIN, context.comm()),
                      "MPI_Allreduce passive scalar bit-pattern minimum");
  converge_mpi_result(context,
                      MPI_Allreduce(&local_bits, &maximum_bits, 1, MPI_UINT64_T,
                                    MPI_MAX, context.comm()),
                      "MPI_Allreduce passive scalar bit-pattern maximum");
  converge(context, locally_valid && minimum_bits == maximum_bits,
           failure_message);
}

std::size_t checked_owned_count(runtime::Int3 extent) {
  constexpr std::size_t limit = std::numeric_limits<std::size_t>::max();
  const auto x = static_cast<std::size_t>(extent.x);
  const auto y = static_cast<std::size_t>(extent.y);
  const auto z = static_cast<std::size_t>(extent.z);
  if (y != 0U && x > limit / y) {
    throw runtime::Error("passive scalar owned-cell count overflow");
  }
  const std::size_t xy = x * y;
  if (z != 0U && xy > limit / z) {
    throw runtime::Error("passive scalar owned-cell count overflow");
  }
  return xy * z;
}

std::size_t owned_index(runtime::Int3 extent, int i, int j, int k) noexcept {
  return (static_cast<std::size_t>(k) * static_cast<std::size_t>(extent.y) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(extent.x) +
         static_cast<std::size_t>(i);
}

double value_at(const runtime::FieldView<const double> &field,
                runtime::Int3 cell) {
  return field(cell.x, cell.y, cell.z, 0);
}

double face_flux(const runtime::FieldView<const double> &field,
                 runtime::Int3 lower, runtime::Int3 direction,
                 double velocity) {
  if (velocity == 0.0) {
    return 0.0;
  }
  const auto shifted = [lower, direction](int amount) {
    return runtime::Int3{lower.x + amount * direction.x,
                         lower.y + amount * direction.y,
                         lower.z + amount * direction.z};
  };
  const double minus_one = value_at(field, shifted(-1));
  const double lower_value = value_at(field, shifted(0));
  const double upper_value = value_at(field, shifted(1));
  const double plus_two = value_at(field, shifted(2));
  const double lower_slope =
      mc_limiter(lower_value - minus_one, upper_value - lower_value);
  const double upper_slope =
      mc_limiter(upper_value - lower_value, plus_two - upper_value);
  const double left_state = lower_value + 0.5 * lower_slope;
  const double right_state = upper_value - 0.5 * upper_slope;
  return velocity * (velocity > 0.0 ? left_state : right_state);
}

double spatial_operator(const runtime::FieldView<const double> &field,
                        runtime::Int3 cell, runtime::Real3 spacing,
                        runtime::Real3 velocity) {
  const double x_plus =
      face_flux(field, cell, runtime::Int3{1, 0, 0}, velocity.x);
  const double x_minus =
      face_flux(field, runtime::Int3{cell.x - 1, cell.y, cell.z},
                runtime::Int3{1, 0, 0}, velocity.x);
  const double y_plus =
      face_flux(field, cell, runtime::Int3{0, 1, 0}, velocity.y);
  const double y_minus =
      face_flux(field, runtime::Int3{cell.x, cell.y - 1, cell.z},
                runtime::Int3{0, 1, 0}, velocity.y);
  const double z_plus =
      face_flux(field, cell, runtime::Int3{0, 0, 1}, velocity.z);
  const double z_minus =
      face_flux(field, runtime::Int3{cell.x, cell.y, cell.z - 1},
                runtime::Int3{0, 0, 1}, velocity.z);
  return -(x_plus - x_minus) / spacing.x - (y_plus - y_minus) / spacing.y -
         (z_plus - z_minus) / spacing.z;
}

void validate_diagnostic_view(const mesh::UniformStructuredMesh &mesh,
                              const runtime::FieldStorage &storage,
                              const runtime::FieldView<const double> &view) {
  if (!same(storage.interior_extent(), mesh.local_extent()) ||
      !same(view.interior_extent(), mesh.local_extent())) {
    throw runtime::Error(
        "passive scalar diagnostic extent does not match the mesh");
  }
  if (view.components() != 1U) {
    throw runtime::Error(
        "passive scalar diagnostic field must have one component");
  }
  if (!std::isfinite(mesh.cell_volume_m3()) || mesh.cell_volume_m3() <= 0.0) {
    throw runtime::Error(
        "passive scalar diagnostic cell volume must be finite and positive");
  }
}

} // namespace

double mc_limiter(double left, double right) noexcept {
  if (!std::isfinite(left) || !std::isfinite(right) || left == 0.0 ||
      right == 0.0 || std::signbit(left) != std::signbit(right)) {
    return 0.0;
  }
  const double sign = left > 0.0 ? 1.0 : -1.0;
  return sign * std::min({2.0 * std::abs(left), 0.5 * std::abs(left + right),
                          2.0 * std::abs(right)});
}

PassiveScalarSolver::PassiveScalarSolver(
    const runtime::MpiContext &context,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::UniformStructuredMesh &mesh, runtime::HaloExchange &halo,
    runtime::Real3 velocity_m_per_s, double diffusivity_m2_per_s) {
  runtime::detail::require_mpi_active("construct passive scalar solver");

  // The constructor's collective preflight order is part of the public
  // contract: communicator, velocity x/y/z, diffusion, periodicity, geometry.
  const bool local_communicators_present =
      context.comm() != MPI_COMM_NULL && decomposition.comm() != MPI_COMM_NULL;
  converge(context, local_communicators_present,
           "passive scalar solver requires valid communicators");
  int communicator_relation = MPI_UNEQUAL;
  converge_mpi_result(context,
                      MPI_Comm_compare(context.comm(), decomposition.comm(),
                                       &communicator_relation),
                      "MPI_Comm_compare passive scalar communicators");
  const bool communicators_compatible = communicator_relation == MPI_IDENT ||
                                        communicator_relation == MPI_CONGRUENT;
  converge(context, communicators_compatible,
           "passive scalar communicators must be identical or congruent");

  require_agreed_value(context, velocity_m_per_s.x,
                       std::isfinite(velocity_m_per_s.x),
                       "passive scalar velocity x is invalid or differs");
  require_agreed_value(context, velocity_m_per_s.y,
                       std::isfinite(velocity_m_per_s.y),
                       "passive scalar velocity y is invalid or differs");
  require_agreed_value(context, velocity_m_per_s.z,
                       std::isfinite(velocity_m_per_s.z),
                       "passive scalar velocity z is invalid or differs");
  require_agreed_value(
      context, diffusivity_m2_per_s,
      std::isfinite(diffusivity_m2_per_s) && diffusivity_m2_per_s == 0.0,
      "Stage 1 passive scalar diffusion must be finite and zero");

  const auto periodic = decomposition.periodic();
  converge(context, periodic[0] && periodic[1] && periodic[2],
           "Stage 1 passive scalar mesh must be periodic on every axis");

  const runtime::Real3 spacing = mesh.spacing_m();
  const bool geometry_ok =
      same(mesh.local_extent(), decomposition.local_extent()) &&
      same(mesh.owned_global_box(), decomposition.owned_box()) &&
      std::isfinite(spacing.x) && spacing.x > 0.0 && std::isfinite(spacing.y) &&
      spacing.y > 0.0 && std::isfinite(spacing.z) && spacing.z > 0.0 &&
      std::isfinite(mesh.cell_volume_m3()) && mesh.cell_volume_m3() > 0.0;
  converge(context, geometry_ok,
           "passive scalar mesh geometry is incompatible");

  int halo_ghost_width = 0;
  bool local_halo_ok = true;
  std::string retained_halo_message;
  std::string_view local_halo_message;
  try {
    halo_ghost_width = halo.ghost_width();
    if (!halo.is_compatible_with(decomposition)) {
      local_halo_ok = false;
      local_halo_message = kHaloCompatibilityFallback;
    }
  } catch (const runtime::Error &error) {
    local_halo_ok = false;
    retain_project_error(error, retained_halo_message, local_halo_message,
                         kHaloCompatibilityFallback);
  } catch (...) {
    local_halo_ok = false;
    local_halo_message = kHaloCompatibilityFallback;
  }
  converge(context, local_halo_ok, local_halo_message);
  require_agreed_value(context, static_cast<double>(halo_ghost_width),
                       halo_ghost_width >= 2,
                       "passive scalar halo width is invalid or differs");

  context_ = &context;
  halo_ = &halo;
  local_extent_ = decomposition.local_extent();
  spacing_m_ = spacing;
  velocity_m_per_s_ = velocity_m_per_s;
  halo_ghost_width_ = halo_ghost_width;
}

void PassiveScalarSolver::advance_ssprk2(runtime::FieldStorage &storage,
                                         runtime::FieldId scalar,
                                         runtime::FieldId stage, double dt_s) {
  runtime::detail::require_mpi_active("advance passive scalar SSPRK2 step");
  require_agreed_value(*context_, dt_s, std::isfinite(dt_s) && dt_s > 0.0,
                       "passive scalar time step is invalid or differs");

  std::optional<runtime::FieldView<double>> scalar_view;
  std::optional<runtime::FieldView<double>> stage_view;
  std::vector<double> old_values;
  std::vector<double> final_values;
  bool local_preflight_ok = true;
  std::string_view local_preflight_message;
  try {
    // Fixed local preflight order: IDs, storage extent, scalar, stage, buffers.
    if (scalar == stage) {
      throw runtime::Error("passive scalar and stage fields must differ");
    }
    if (!same(storage.interior_extent(), local_extent_)) {
      throw runtime::Error(
          "passive scalar storage extent does not match the decomposition");
    }
    scalar_view.emplace(storage.view<double>(scalar));
    if (scalar_view->components() != 1U) {
      throw runtime::Error("passive scalar field must have one component");
    }
    if (!same(scalar_view->interior_extent(), local_extent_)) {
      throw runtime::Error("passive scalar field extent is incompatible");
    }
    if (scalar_view->ghost_width() < halo_ghost_width_) {
      throw runtime::Error(
          "passive scalar field has fewer ghosts than the halo width");
    }
    stage_view.emplace(storage.view<double>(stage));
    if (stage_view->components() != 1U) {
      throw runtime::Error("passive scalar stage must have one component");
    }
    if (!same(stage_view->interior_extent(), local_extent_)) {
      throw runtime::Error("passive scalar stage extent is incompatible");
    }
    if (stage_view->ghost_width() < halo_ghost_width_) {
      throw runtime::Error(
          "passive scalar stage has fewer ghosts than the halo width");
    }

    const std::size_t count = checked_owned_count(local_extent_);
    old_values.resize(count);
    final_values.resize(count);
    for (int k = 0; k < local_extent_.z; ++k) {
      for (int j = 0; j < local_extent_.y; ++j) {
        for (int i = 0; i < local_extent_.x; ++i) {
          old_values[owned_index(local_extent_, i, j, k)] =
              (*scalar_view)(i, j, k, 0);
        }
      }
    }
  } catch (const std::exception &) {
    local_preflight_ok = false;
    local_preflight_message = "passive scalar field preflight failed";
  } catch (...) {
    local_preflight_ok = false;
    local_preflight_message = "passive scalar field preflight failed";
  }
  converge(*context_, local_preflight_ok, local_preflight_message);

  bool local_step_ok = true;
  std::string retained_step_message;
  std::string_view local_step_message;
  try {
    halo_->exchange(storage, scalar);
  } catch (const runtime::Error &error) {
    local_step_ok = false;
    retain_project_error(error, retained_step_message, local_step_message,
                         kHaloExchangeFallback);
  } catch (...) {
    local_step_ok = false;
    local_step_message = kHaloExchangeFallback;
  }
  converge(*context_, local_step_ok, local_step_message);

  local_step_ok = true;
  local_step_message = {};
  try {
    const auto source =
        static_cast<const runtime::FieldStorage &>(storage).view<double>(
            scalar);
    for (int k = 0; k < local_extent_.z; ++k) {
      for (int j = 0; j < local_extent_.y; ++j) {
        for (int i = 0; i < local_extent_.x; ++i) {
          const runtime::Int3 cell{i, j, k};
          const std::size_t index = owned_index(local_extent_, i, j, k);
          (*stage_view)(i, j, k, 0) =
              old_values[index] + dt_s * spatial_operator(source, cell,
                                                          spacing_m_,
                                                          velocity_m_per_s_);
        }
      }
    }
  } catch (const std::exception &) {
    local_step_ok = false;
    local_step_message = "passive scalar first stage failed";
  } catch (...) {
    local_step_ok = false;
    local_step_message = "passive scalar first stage failed";
  }
  converge(*context_, local_step_ok, local_step_message);

  local_step_ok = true;
  local_step_message = {};
  std::string retained_stage_message;
  try {
    halo_->exchange(storage, stage);
  } catch (const runtime::Error &error) {
    local_step_ok = false;
    retain_project_error(error, retained_stage_message, local_step_message,
                         kStageExchangeFallback);
  } catch (...) {
    local_step_ok = false;
    local_step_message = kStageExchangeFallback;
  }
  converge(*context_, local_step_ok, local_step_message);

  local_step_ok = true;
  local_step_message = {};
  try {
    const auto source =
        static_cast<const runtime::FieldStorage &>(storage).view<double>(stage);
    for (int k = 0; k < local_extent_.z; ++k) {
      for (int j = 0; j < local_extent_.y; ++j) {
        for (int i = 0; i < local_extent_.x; ++i) {
          const runtime::Int3 cell{i, j, k};
          const std::size_t index = owned_index(local_extent_, i, j, k);
          final_values[index] =
              0.5 * old_values[index] +
              0.5 * ((*stage_view)(i, j, k, 0) +
                     dt_s * spatial_operator(source, cell, spacing_m_,
                                             velocity_m_per_s_));
        }
      }
    }
  } catch (const std::exception &) {
    local_step_ok = false;
    local_step_message = "passive scalar second stage failed";
  } catch (...) {
    local_step_ok = false;
    local_step_message = "passive scalar second stage failed";
  }
  converge(*context_, local_step_ok, local_step_message);

  for (int k = 0; k < local_extent_.z; ++k) {
    for (int j = 0; j < local_extent_.y; ++j) {
      for (int i = 0; i < local_extent_.x; ++i) {
        (*scalar_view)(i, j, k, 0) =
            final_values[owned_index(local_extent_, i, j, k)];
      }
    }
  }
}

double global_mass(const runtime::MpiContext &context,
                   const mesh::UniformStructuredMesh &mesh,
                   const runtime::FieldStorage &storage,
                   runtime::FieldId scalar) {
  std::optional<runtime::FieldView<const double>> view;
  bool local_ok = true;
  std::string_view local_message;
  try {
    view.emplace(storage.view<double>(scalar));
    validate_diagnostic_view(mesh, storage, *view);
  } catch (const std::exception &) {
    local_ok = false;
    local_message = "passive scalar mass preflight failed";
  } catch (...) {
    local_ok = false;
    local_message = "passive scalar mass preflight failed";
  }
  converge(context, local_ok, local_message);

  double local_mass = 0.0;
  const runtime::Int3 extent = mesh.local_extent();
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        local_mass += (*view)(i, j, k, 0) * mesh.cell_volume_m3();
      }
    }
  }
  double result = 0.0;
  converge_mpi_result(context,
                      MPI_Allreduce(&local_mass, &result, 1, MPI_DOUBLE,
                                    MPI_SUM, context.comm()),
                      "MPI_Allreduce passive scalar global mass");
  return result;
}

double global_l1_error(const runtime::MpiContext &context,
                       const mesh::UniformStructuredMesh &mesh,
                       const runtime::FieldStorage &storage,
                       runtime::FieldId actual, runtime::FieldId reference) {
  std::optional<runtime::FieldView<const double>> actual_view;
  std::optional<runtime::FieldView<const double>> reference_view;
  bool local_ok = true;
  std::string_view local_message;
  try {
    actual_view.emplace(storage.view<double>(actual));
    validate_diagnostic_view(mesh, storage, *actual_view);
    reference_view.emplace(storage.view<double>(reference));
    validate_diagnostic_view(mesh, storage, *reference_view);
  } catch (const std::exception &) {
    local_ok = false;
    local_message = "passive scalar L1 preflight failed";
  } catch (...) {
    local_ok = false;
    local_message = "passive scalar L1 preflight failed";
  }
  converge(context, local_ok, local_message);

  const runtime::Int3 extent = mesh.local_extent();
  double local_error = 0.0;
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        local_error += std::abs((*actual_view)(i, j, k, 0) -
                                (*reference_view)(i, j, k, 0)) *
                       mesh.cell_volume_m3();
      }
    }
  }
  const double local_volume =
      static_cast<double>(runtime::volume(extent)) * mesh.cell_volume_m3();
  double global_error = 0.0;
  double global_volume = 0.0;
  converge_mpi_result(context,
                      MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE,
                                    MPI_SUM, context.comm()),
                      "MPI_Allreduce passive scalar global L1 numerator");
  converge_mpi_result(context,
                      MPI_Allreduce(&local_volume, &global_volume, 1,
                                    MPI_DOUBLE, MPI_SUM, context.comm()),
                      "MPI_Allreduce passive scalar global L1 volume");
  if (!std::isfinite(global_volume) || global_volume <= 0.0) {
    throw runtime::Error(
        "passive scalar global volume must be finite and positive");
  }
  return global_error / global_volume;
}

} // namespace hundun::solver
