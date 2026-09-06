// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_flow.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_scalar_boundary_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

namespace hundun::v04::detail {

// Conservative transport of the mass correction, not a cellwise density
// rescaling. The predictor's *actual paired* flux includes its BDF/EX weights
// and any common conservative limiter. For dm_f=(phi_final-phi_star)/a0:
//   M_final q + sum_f dm_f q_upwind = (Mq)_star.
// The same operator transports every species, including the implicit balance
// species. Its row sum is M_star to the continuity residual. Positive masses
// therefore give an M-matrix and preserve the predictor's bounds. Donor order
// applies only to the flux correction (O(dt^2) in a smooth BDF2 step), not to
// the configured predictor convection scheme.
class ScalarMassRemap {
 public:
  static constexpr std::uint32_t kInvalid = 10214U;
  static constexpr std::uint32_t kNonconverged = 10215U;
  static constexpr std::uint32_t kRecouple = 10216U;
  static constexpr unsigned maximum_coupling_sweeps = 16U;
  static constexpr double tolerance =
      32.0 * std::numeric_limits<double>::epsilon();
  // Absolute composition backward error, after a complete EOS/momentum/
  // limiter recoupling. Distinct from the inner matrix's relative residual:
  // requiring the latter of the entire floating-point nonlinear map causes
  // roundoff-level AFC branch chatter. Independent inventory tests retain
  // their finite-dt 1e-12 gate; no global inventory correction is applied.
  static constexpr double composition_tolerance =
      128.0 * std::numeric_limits<double>::epsilon();

  struct Report {
    unsigned iterations{};
    double initial_species_residual{};
    double residual{};
    double mass_pairing_residual{};
  };

  // Local allocation and collective communication binding are separate so
  // the caller can agree bad_alloc before any rank enters HaloEngine::reserve.
  Status allocate(const MeshPatch& patch, Span<const FieldId> fields,
                  Span<const TransportedScalarRole> roles,
                  std::uint8_t ghosts) {
    cells_ = patch.cells;
    patch_ = patch;
    ghosts_ = ghosts;
    if (!valid_cells(cells_) || fields.size != roles.size ||
        fields.size == 0U || ghosts == 0U)
      return {StatusCode::invalid_plan, kInvalid};
    const auto multiply = [](std::size_t a, std::size_t b,
                             std::size_t& out) {
      if (b && a > std::numeric_limits<std::size_t>::max() / b) return false;
      out = a * b;
      return true;
    };
    std::size_t sy = static_cast<std::size_t>(cells_.x) + 2U * ghosts;
    std::size_t sz{}, stride{}, count{};
    if (!multiply(sy, static_cast<std::size_t>(cells_.y) + 2U * ghosts, sz) ||
        !multiply(sz, static_cast<std::size_t>(cells_.z) + 2U * ghosts, stride) ||
        !multiply(static_cast<std::size_t>(cells_.x), cells_.y, count) ||
        !multiply(count, cells_.z, count) ||
        !multiply(count, fields.size, quantity_count_))
      return {StatusCode::invalid_plan, kInvalid};
    count_ = count;
    mass_.resize(count);
    quantity_.resize(quantity_count_);
    next_.resize(quantity_count_);
    storage_.resize(fields.size);
    views_.resize(fields.size);
    halo_specs_.resize(fields.size);
    roles_.assign(roles.data, roles.data + roles.size);
    const auto passive_count = static_cast<std::size_t>(std::count(
        roles_.begin(), roles_.end(), TransportedScalarRole::passive_scalar));
    intervals_.resize(passive_count);
    bounds_local_.resize(2U * passive_count);
    bounds_global_.resize(2U * passive_count);
    for (std::size_t s = 0U; s < fields.size; ++s) {
      storage_[s].assign(stride, 0.0);
      FieldView v;
      v.base = storage_[s].data() + ghosts + ghosts * sy + ghosts * sz;
      v.interior = cells_;
      v.ghosts = {ghosts, ghosts, ghosts};
      v.components = 1U;
      v.stride_y = sy;
      v.stride_z = sz;
      v.component_stride = stride;
      v.field = fields.data[s];
      v.revision = 1U;
      v.storage_identity = reinterpret_cast<StorageIdentity>(storage_[s].data());
      v.revision_domain = reinterpret_cast<RevisionDomainIdentity>(this);
      views_[s] = v;
      halo_specs_[s] = {v.field, ghosts, 1U};
    }
    Status status = FaceFluxStorage::allocate_workspace(cells_, 1U, flux_);
    if (status) status = flux_.workspace_view(0U, 1U, predictor_flux_);
    return status;
  }

  Status bind(MPI_Comm comm, const BoundaryPlan& boundary) noexcept {
    boundary_ = &boundary;
    return halo_.reserve(comm, patch_, {halo_specs_.data(), halo_specs_.size()},
                         boundary.halo_topology());
  }

  const HaloEngine& halo() const noexcept { return halo_; }

  // Owned capacities, not a process RSS/peak estimate. Halo/MPI bookkeeping
  // is reported by its own service and is deliberately not double-counted.
  std::uint64_t owned_payload_bytes() const noexcept {
    std::uint64_t bytes = flux_.counters().aligned_payload_bytes;
    const auto add = [&](const auto& values) {
      using Element = typename std::decay_t<decltype(values)>::value_type;
      const auto count = values.capacity();
      if (count > (UINT64_MAX - bytes) / sizeof(Element)) bytes = UINT64_MAX;
      else bytes += count * sizeof(Element);
    };
    add(mass_); add(quantity_); add(next_); add(storage_);
    for (const auto& values : storage_) add(values);
    add(views_); add(halo_specs_); add(roles_); add(intervals_);
    add(bounds_local_); add(bounds_global_);
    return bytes;
  }

  Status prepare_passive_intervals(ThermophysicalPredictorInput& input,
                                  ReductionEngine& reductions) noexcept {
    if (intervals_.empty()) return {};
    Status local;
    if (input.passive_scalars_accepted.size != intervals_.size() ||
        input.passive_scalar_nonadvective_rhs.size != intervals_.size())
      local = {StatusCode::invalid_plan, kInvalid};
    std::fill(bounds_local_.begin(), bounds_local_.end(),
              -std::numeric_limits<double>::max());
    for (std::size_t s = 0U; s < intervals_.size() && local; ++s) {
      const auto q = input.passive_scalars_accepted.data[s];
      const auto rate = input.passive_scalar_nonadvective_rhs.data[s].accepted;
      const auto include = [&](double value) {
        if (!std::isfinite(value)) { local = {StatusCode::rejected_step, kInvalid}; return; }
        bounds_local_[2U*s] = std::max(bounds_local_[2U*s], -value);
        bounds_local_[2U*s+1U] = std::max(bounds_local_[2U*s+1U], value);
      };
      for (int z = 0; z < cells_.z; ++z)
        for (int y = 0; y < cells_.y; ++y)
          for (int x = 0; x < cells_.x; ++x) {
            const Int3 c{x,y,z};
            const double value = q.unchecked(c,0U);
            include(value);
            // Include the explicit physical/source endpoint. This permits
            // signed/source-driven passives without imposing [0,1]. A stable
            // source-free diffusion step does not enlarge the global range.
            if (rate.base)
              include(value + input.dt * rate.unchecked(c,0U) /
                                input.density_accepted.unchecked(c,0U));
            for (int f = 0; f < 6; ++f) {
              const BoundaryFacePlan* face = nullptr;
              const int axis = f/2, n = axis==0 ? x : axis==1 ? y : z;
              const int extent = axis==0 ? cells_.x : axis==1 ? cells_.y : cells_.z;
              if (n != (f%2 ? extent-1 : 0) ||
                  !boundary_->face(static_cast<CartesianFace>(f),face) ||
                  !face->local_owner || face->periodic) continue;
              Int3 donor = c;
              (axis==0 ? donor.x : axis==1 ? donor.y : donor.z) += f%2 ? 1 : -1;
              include(0.5*(value + q.unchecked(donor,0U)));
            }
          }
    }
    Status status = reductions.checked_max(
        {bounds_local_.data(), bounds_local_.size()},
        {bounds_global_.data(), bounds_global_.size()}, local);
    if (!status) return status;
    for (std::size_t s=0U; s<intervals_.size(); ++s) {
      const double lower=-bounds_global_[2U*s], upper=bounds_global_[2U*s+1U];
      const double guard=128.0*std::numeric_limits<double>::epsilon()*
                         std::max(std::abs(lower),std::abs(upper));
      intervals_[s]={lower-guard,upper+guard};
    }
    input.passive_intervals={intervals_.data(),intervals_.size()};
    return {};
  }

  Status capture(const CartesianKernelPlan& kernels, ConstFieldView density,
                 Span<const FieldView> scalars, ConstFaceFluxView flux,
                 double a0) noexcept {
    if (!std::isfinite(a0) || a0 <= 0.0 ||
        scalars.size != views_.size() ||
        !valid_cell_view(density, cells_, 0U, 1U, 0U) ||
        !valid_flux_view(flux, cells_, flux.revision))
      return {StatusCode::invalid_plan, kInvalid};
    a0_ = a0;
    const std::array<ConstFaceFieldView, 3U> source{flux.x, flux.y, flux.z};
    const std::array<FaceFieldView, 3U> target{
        predictor_flux_.x, predictor_flux_.y, predictor_flux_.z};
    for (std::size_t a = 0U; a < 3U; ++a)
      for (int z = 0; z < source[a].extents.z; ++z)
        for (int y = 0; y < source[a].extents.y; ++y)
          for (int x = 0; x < source[a].extents.x; ++x)
            target[a].unchecked({x, y, z}) = source[a].unchecked({x, y, z});
    for (std::size_t s = 0U; s < scalars.size; ++s)
      if (scalars.data[s].field != views_[s].field ||
          !valid_cell_view(as_const(scalars.data[s]), cells_, 0U, 1U, 0U))
        return {StatusCode::invalid_plan, kInvalid};
    std::size_t i = 0U;
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x, ++i) {
          const Int3 c{x, y, z};
          mass_[i] = density.unchecked(c, 0U) * cell_volume(kernels, c);
          if (!std::isfinite(mass_[i]) || mass_[i] <= 0.0)
            return {StatusCode::rejected_step, kInvalid};
          for (std::size_t s = 0U; s < scalars.size; ++s) {
            const double q = scalars.data[s].unchecked(c, 0U);
            if (!std::isfinite(q)) return {StatusCode::rejected_step, kInvalid};
            quantity_[s * count_ + i] = mass_[i] * q;
          }
        }
    return {};
  }

  Status solve(const CartesianKernelPlan& kernels, ConstFieldView density,
               ConstFieldView velocity, ConstFaceFluxView flux, Span<const FieldView> scalars,
               Span<const std::uint8_t> activity,
               BoundaryResolvedValues boundary_values,
               ReductionEngine& reductions, Report& report) noexcept {
    report = {};
    Status local;
    if (boundary_ == nullptr || scalars.size != views_.size() ||
        !valid_cell_view(density, cells_, 0U, 1U, 0U) ||
        !valid_cell_view(velocity, cells_, 0U, 3U, 0U) ||
        !valid_flux_view(flux, cells_, flux.revision) ||
        (activity.size != 0U && (activity.size != count_ || !activity.data)))
      local = {StatusCode::invalid_plan, kInvalid};
    for (std::size_t s = 0U; s < scalars.size && local; ++s) {
      if (scalars.data[s].field != views_[s].field ||
          !valid_cell_view(as_const(scalars.data[s]), cells_, 0U, 1U, 0U)) {
        local = {StatusCode::invalid_plan, kInvalid};
        break;
      }
      std::size_t i = 0U;
      for (int z = 0; z < cells_.z; ++z)
        for (int y = 0; y < cells_.y; ++y)
          for (int x = 0; x < cells_.x; ++x, ++i)
          {
            views_[s].unchecked({x, y, z}, 0U) =
                scalars.data[s].unchecked({x, y, z}, 0U);
          }
    }
    Status status = reductions.consensus(local);
    if (!status) return status;
    const std::array<ConstFaceFieldView, 3U> final{flux.x, flux.y, flux.z};
    const ConstFaceFluxView prediction = as_const(predictor_flux_);
    const std::array<ConstFaceFieldView, 3U> base{
        prediction.x, prediction.y, prediction.z};
    for (unsigned iteration = 0U; iteration < 128U; ++iteration) {
      for (FieldView& v : views_) ++v.revision;
      HaloTicket ticket;
      status = halo_.begin(176U, {views_.data(), views_.size()}, {}, ticket);
      if (status) status = halo_.finish(ticket, {views_.data(), views_.size()});
      if (status)
        status = apply_boundary_ghosts(BoundaryStage::scalar, *boundary_,
            {views_.data(), views_.size()}, boundary_values);
      status = reductions.consensus(status);
      if (!status) return status;
      double maximum[3U]{};
      std::size_t i = 0U;
      for (int z = 0; z < cells_.z && local; ++z)
        for (int y = 0; y < cells_.y && local; ++y)
          for (int x = 0; x < cells_.x; ++x, ++i) {
            const Int3 c{x, y, z};
            const double mass = density.unchecked(c, 0U) * cell_volume(kernels, c);
            const bool active = activity.size == 0U || activity.data[i] != 0U;
            std::array<double, 6U> outward{};
            std::array<Int3, 6U> neighbours{};
            double diagonal = mass, signed_mass = 0.0;
            for (std::size_t f = 0U; f < 6U; ++f) {
              const std::size_t axis = f / 2U;
              Int3 face = c, neighbour = c;
              auto coordinate = [&](Int3& point) -> int& {
                return axis == 0U ? point.x : axis == 1U ? point.y : point.z;
              };
              const int sign = (f % 2U) == 0U ? -1 : 1;
              if (sign > 0) ++coordinate(face);
              coordinate(neighbour) += sign;
              neighbours[f] = neighbour;
              outward[f] = sign * (final[axis].unchecked(face) -
                                  base[axis].unchecked(face)) / a0_;
              diagonal += std::max(outward[f], 0.0);
              signed_mass += outward[f];
            }
            if (!std::isfinite(diagonal) || diagonal <= 0.0 ||
                !std::isfinite(mass) || mass <= 0.0) {
              local = {StatusCode::rejected_step, kInvalid};
              break;
            }
            if (active)
              maximum[2U] = std::max(maximum[2U],
                  std::abs((mass - mass_[i]) + signed_mass) /
                      (mass + mass_[i]));
            for (std::size_t s = 0U; s < views_.size(); ++s) {
              const double q = views_[s].unchecked(c, 0U);
              long double rhs = quantity_[s * count_ + i];
              for (std::size_t f = 0U; f < 6U; ++f)
                if (outward[f] < 0.0)
                  rhs -= static_cast<long double>(outward[f]) *
                         scalar_upwind_donor(*boundary_,as_const(views_[s]),neighbours[f],velocity);
              const long double residual = static_cast<long double>(diagonal) * q - rhs;
              const long double scale = std::abs(static_cast<long double>(diagonal) * q) +
                                        std::abs(rhs);
              const double norm = scale == 0.0L ? 0.0 :
                  static_cast<double>(std::abs(residual) / scale);
              const double value = active ? static_cast<double>(rhs / diagonal) : q;
              if (!std::isfinite(value) || !std::isfinite(norm)) {
                local = {StatusCode::rejected_step, kInvalid};
                break;
              }
              next_[s * count_ + i] = value;
              if (active) {
                maximum[0U] = std::max(maximum[0U], norm);
                if (roles_[s] == TransportedScalarRole::species)
                  // Composition coupling is normalized by mixture mass,
                  // i.e. absolute mass-fraction error, including trace species.
                  // The inner remap solve still uses each equation's relative
                  // residual; no inventory is clipped or renormalized.
                  maximum[1U] = std::max(maximum[1U],
                      static_cast<double>(std::abs(residual)/(mass+mass_[i])));
              }
            }
          }
      double global[3U]{};
      status = reductions.checked_max({maximum, 3U}, {global, 3U}, local);
      if (!status) return status;
      if (iteration == 0U) report.initial_species_residual = global[1U];
      report.iterations = iteration;
      report.residual = global[0U];
      report.mass_pairing_residual = global[2U];
      if (report.residual <= tolerance) return {};
      for (std::size_t s = 0U; s < views_.size(); ++s) {
        i = 0U;
        for (int z = 0; z < cells_.z; ++z)
          for (int y = 0; y < cells_.y; ++y)
            for (int x = 0; x < cells_.x; ++x, ++i)
              views_[s].unchecked({x, y, z}, 0U) = next_[s * count_ + i];
      }
    }
    return {StatusCode::rejected_step, kNonconverged};
  }

  void copy_solution(Span<const FieldView> target,
                     TransportedScalarRole role) const noexcept {
    for (std::size_t s = 0U; s < views_.size(); ++s) {
      if (roles_[s] != role) continue;
      for (int z = 0; z < cells_.z; ++z)
        for (int y = 0; y < cells_.y; ++y)
          for (int x = 0; x < cells_.x; ++x)
            target.data[s].unchecked({x, y, z}, 0U) =
                views_[s].unchecked({x, y, z}, 0U);
    }
  }

 private:
  Int3 cells_{};
  MeshPatch patch_{};
  std::uint8_t ghosts_{};
  std::size_t count_{}, quantity_count_{};
  double a0_{};
  const BoundaryPlan* boundary_{};
  std::vector<double> mass_, quantity_, next_;
  std::vector<double> bounds_local_, bounds_global_;
  std::vector<ScalarAdmissibleInterval> intervals_;
  std::vector<std::vector<double>> storage_;
  std::vector<FieldView> views_;
  std::vector<HaloFieldSpec> halo_specs_;
  std::vector<TransportedScalarRole> roles_;
  HaloEngine halo_;
  FaceFluxStorage flux_;
  FaceFluxView predictor_flux_{};
};

}  // namespace hundun::v04::detail
