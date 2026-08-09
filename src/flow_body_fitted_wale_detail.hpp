// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "flow_fixed_step_detail.hpp"

#include "hundun/fvm_cell_centered.hpp"
#include "hundun/les_wale.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_field_view.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>

namespace hundun::flow::detail {

struct BodyFittedWaleAttemptData final {
  std::optional<les::WaleAttemptCoefficients> coefficients;
  std::optional<runtime::FaceFieldView<const double>>
      sgs_dynamic_viscosity_by_face;
  std::optional<runtime::FaceFieldView<const double>>
      effective_dynamic_viscosity_by_face;
  les::WaleSummary summary;
};

namespace body_fitted_wale_detail {

struct StructuredIndex final {
  int i{};
  int j{};
  int k{};
};

inline StructuredIndex map_cell(runtime::Int3 global, runtime::Box3 owned,
                                runtime::Int3 global_extent) {
  const runtime::Int3 local{owned.end.x - owned.begin.x,
                            owned.end.y - owned.begin.y,
                            owned.end.z - owned.begin.z};
  const auto axis = [](int coordinate, int begin, int end, int global_n,
                       int local_n) {
    if (coordinate >= begin && coordinate < end)
      return coordinate - begin;
    if (coordinate == begin - 1 ||
        (begin == 0 && coordinate == global_n - 1))
      return -1;
    if (coordinate == end || (end == global_n && coordinate == 0))
      return local_n;
    throw runtime::Error("WALE cell has no structured field mapping");
  };
  return {axis(global.x, owned.begin.x, owned.end.x, global_extent.x,
               local.x),
          axis(global.y, owned.begin.y, owned.end.y, global_extent.y,
               local.y),
          axis(global.z, owned.begin.z, owned.end.z, global_extent.z,
               local.z)};
}

inline std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

inline std::uint64_t hash_value(std::uint64_t hash,
                                std::uint64_t value) noexcept {
  constexpr std::uint64_t prime = UINT64_C(1099511628211);
  for (std::size_t byte = 0U; byte < sizeof(value); ++byte) {
    hash ^= (value >> (8U * byte)) & UINT64_C(0xff);
    hash *= prime;
  }
  return hash;
}

} // namespace body_fitted_wale_detail

template <class FlowImplementation>
BodyFittedWaleAttemptData prepare_body_fitted_wale_attempt(
    FlowImplementation &impl, FlowState &state, const les::WaleModel &model,
    double molecular_mu, const MomentumTimeStencil &stencil,
    runtime::FieldView<const double> rho_attempt,
    std::uint64_t &evaluation_count,
    runtime::PhaseId state_phase, runtime::ActorId state_actor,
    runtime::PhaseId scratch_phase, runtime::ActorId scratch_actor) {
  if (evaluation_count != 0U)
    throw runtime::Error("body-fitted WALE was evaluated more than once");
  ++evaluation_count;
  auto &committed = FlowStateSolverAccess::layer(state, FlowLayer::committed);
  auto &history = FlowStateSolverAccess::layer(state, FlowLayer::history);
  const auto &access = FlowStateSolverAccess::access(state);
  const auto fields = state.fields();
  const runtime::Int3 extent = impl.decomposition->local_extent();
  auto lagged = impl.scratch.storage->template acquire_write<double>(
      *impl.scratch.access, scratch_phase, scratch_actor,
      impl.scratch.lagged_velocity);
  const auto velocity_n = committed.template acquire_read<double>(
      access, state_phase, state_actor, fields.velocity);
  const auto velocity_nm1 = history.template acquire_read<double>(
      access, state_phase, state_actor, fields.velocity);
  const double ratio = stencil.order == MomentumTimeOrder::bdf2
                           ? stencil.dt_s / stencil.previous_dt_s
                           : 0.0;
  for (int k = -2; k < extent.z + 2; ++k) {
    for (int j = -2; j < extent.y + 2; ++j) {
      for (int i = -2; i < extent.x + 2; ++i) {
        for (int component = 0; component < 3; ++component) {
          const double current = velocity_n(i, j, k, component);
          const double previous = velocity_nm1(i, j, k, component);
          const double value =
              stencil.order == MomentumTimeOrder::bdf2
                  ? current + ratio * (current - previous)
                  : current;
          if (!std::isfinite(value))
            throw runtime::Error("WALE lagged velocity is non-finite");
          lagged(i, j, k, component) = value;
        }
      }
    }
  }
  const auto lagged_read =
      impl.scratch.storage->template acquire_read<double>(
          *impl.scratch.access, scratch_phase, scratch_actor,
          impl.scratch.lagged_velocity);
  auto gradient = impl.scratch.storage->template acquire_write<double>(
      *impl.scratch.access, scratch_phase, scratch_actor,
      impl.scratch.velocity_gradient);
  impl.fvm.compute_gradient(finite_volume::GradientScheme::green_gauss,
                            finite_volume::FiniteVolumeQuantity::velocity(),
                            *impl.boundaries, lagged_read, gradient);
  impl.halo->exchange(*impl.scratch.storage,
                      impl.scratch.velocity_gradient);

  const auto gradient_read =
      impl.scratch.storage->template acquire_read<double>(
          *impl.scratch.access, scratch_phase, scratch_actor,
          impl.scratch.velocity_gradient);
  constexpr std::uint64_t offset = UINT64_C(1469598103934665603);
  std::uint64_t committed_hash = offset;
  std::uint64_t history_hash = offset;
  std::uint64_t gradient_hash = offset;
  std::uint64_t density_hash = offset;
  const auto owned = impl.topology->owned_global_box();
  for (mesh::LocalCellId cell = 0U; cell < impl.topology->local_cell_count();
       ++cell) {
    const auto index = body_fitted_wale_detail::map_cell(
        impl.topology->global_cell(cell), owned,
        impl.topology->global_extent());
    for (int component = 0; component < 3; ++component) {
      committed_hash = body_fitted_wale_detail::hash_value(
          committed_hash,
          body_fitted_wale_detail::bits(
              velocity_n(index.i, index.j, index.k, component)));
      history_hash = body_fitted_wale_detail::hash_value(
          history_hash,
          body_fitted_wale_detail::bits(
              velocity_nm1(index.i, index.j, index.k, component)));
    }
    for (int component = 0; component < 9; ++component)
      gradient_hash = body_fitted_wale_detail::hash_value(
          gradient_hash,
          body_fitted_wale_detail::bits(
              gradient_read(index.i, index.j, index.k, component)));
    density_hash = body_fitted_wale_detail::hash_value(
        density_hash,
        body_fitted_wale_detail::bits(
            rho_attempt(index.i, index.j, index.k, 0)));
  }

  BodyFittedWaleAttemptData result;
  result.coefficients.emplace(model.evaluate(
      {state.metadata().step + 1U,
       stencil.dt_s,
       stencil.order == MomentumTimeOrder::bdf2
           ? les::WaleTimeOrder::bdf2
           : les::WaleTimeOrder::backward_euler,
       committed_hash,
       history_hash,
       gradient_hash,
       density_hash,
       gradient_read,
       rho_attempt}));
  if (result.coefficients->owned_active_count() !=
          impl.topology->owned_cell_count() ||
      result.coefficients->local_active_count() !=
          impl.topology->local_cell_count())
    throw runtime::Error("WALE coefficient layout is not topology order");
  const auto mu_sgs = result.coefficients->mu_sgs_pa_s();
  {
    auto sgs_cell = impl.scratch.storage->template acquire_write<double>(
        *impl.scratch.access, scratch_phase, scratch_actor,
        impl.scratch.sgs_dynamic_viscosity_cell);
    for (mesh::LocalCellId cell = 0U;
         cell < impl.topology->local_cell_count(); ++cell) {
      const auto index = body_fitted_wale_detail::map_cell(
          impl.topology->global_cell(cell), owned,
          impl.topology->global_extent());
      const double value = mu_sgs[cell];
      if (!(value >= 0.0) || !std::isfinite(value))
        throw runtime::Error("WALE cell viscosity is invalid");
      sgs_cell(index.i, index.j, index.k, 0) = value;
    }
  }
  impl.halo->exchange(*impl.scratch.storage,
                      impl.scratch.sgs_dynamic_viscosity_cell);
  const auto sgs_cell_read =
      impl.scratch.storage->template acquire_read<double>(
          *impl.scratch.access, scratch_phase, scratch_actor,
          impl.scratch.sgs_dynamic_viscosity_cell);
  auto sgs_face = impl.scratch.storage->template acquire_face_write<double>(
      *impl.scratch.access, scratch_phase, scratch_actor,
      impl.scratch.sgs_dynamic_viscosity);
  impl.fvm.interpolate_cell_scalar_to_faces(sgs_cell_read, sgs_face);
  auto effective_face =
      impl.scratch.storage->template acquire_face_write<double>(
          *impl.scratch.access, scratch_phase, scratch_actor,
          impl.scratch.effective_dynamic_viscosity);
  for (mesh::LocalFaceId face = 0U;
       face < impl.topology->local_face_count(); ++face) {
    const double face_sgs = sgs_face(face, 0);
    const double effective = molecular_mu + face_sgs;
    if (!(face_sgs >= 0.0) || !std::isfinite(face_sgs) ||
        !(effective >= 0.0) || !std::isfinite(effective))
      throw runtime::Error("WALE face viscosity is invalid");
    sgs_face(face, 0) = face_sgs;
    effective_face(face, 0) = effective;
  }
  result.sgs_dynamic_viscosity_by_face.emplace(
      impl.scratch.storage->template acquire_face_read<double>(
          *impl.scratch.access, scratch_phase, scratch_actor,
          impl.scratch.sgs_dynamic_viscosity));
  result.effective_dynamic_viscosity_by_face.emplace(
      impl.scratch.storage->template acquire_face_read<double>(
          *impl.scratch.access, scratch_phase, scratch_actor,
          impl.scratch.effective_dynamic_viscosity));
  result.summary = result.coefficients->summary();
  return result;
}

} // namespace hundun::flow::detail
