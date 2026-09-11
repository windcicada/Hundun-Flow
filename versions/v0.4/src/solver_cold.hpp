// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"
#include "hundun/v04_linear.hpp"
#include "hundun/v04_mpi_runtime.hpp"
#include "solver_cartesian_detail.hpp"

namespace hundun::v04::detail {

// COAST cold-flow row operators, adapted to Hundun field views and linear
// interfaces. Momentum uses CN midpoint storage; mass/energy/species use BE.
// Spatial and temporal identities are verified separately against the COAST
// condif, gvctr, cmod and step routines. IBM solids receive isolated unit rows.

// COAST gvctr/gvctr3 interpolate the coefficient and midpoint velocity
// separately, then multiply. Interpolating their cellwise product differs
// across a density/composition interface.
inline double cold_coefficient_face_flux(
    double lower_coefficient, double upper_coefficient, double lower_weight,
    const std::array<double, 3>& lower_midpoint_velocity,
    const std::array<double, 3>& upper_midpoint_velocity,
    const std::array<double, 3>& oriented_area, bool blocked) noexcept {
  if (blocked) return 0;
  const double coefficient =
      lower_weight * lower_coefficient + (1 - lower_weight) * upper_coefficient;
  double flux = 0;
  for (unsigned i = 0; i < 3; ++i)
    flux += coefficient *
            (lower_weight * lower_midpoint_velocity[i] +
             (1 - lower_weight) * upper_midpoint_velocity[i]) *
            oriented_area[i];
  return flux;
}

// Return metric*gamma for compress. Gradients are (p_i-p_{i-1})/distance,
// (p_{i+1}-p_i)/distance and (p_{i+2}-p_{i+1})/distance. Preserve vls's
// positive-flow IBM fallback to the upper gradient. Thresholds are explicit:
// COAST uses small=1e-30 and epsilon of its REAL kind for uround.
inline double cold_pressure_limiter_diffusion(
    double lower_gradient, double face_gradient, double upper_gradient,
    double lower_weight, double density_flux, bool lower_stencil_allowed,
    double small, double roundoff) noexcept {
  const bool use_lower = density_flux > 0 && lower_stencil_allowed;
  const double upstream = use_lower ? lower_gradient : upper_gradient;
  const double epsilon = std::max(
      small, roundoff * std::max(std::abs(upstream), std::abs(face_gradient)));
  double ratio =
      upstream / (face_gradient + std::copysign(epsilon, face_gradient));
  if (!std::isfinite(ratio)) ratio = 0;
  const double limiter = 1 - std::max(std::min(2 * ratio, 1.0), 0.0);
  return (use_lower ? 1 - lower_weight : lower_weight) *
         std::abs(density_flux) * limiter;
}

// Prototype for the COAST cold pressure correction, in positive-neighbour
// storage: A*x = diagonal*x - sum(neighbour*x_neighbour).
// Faces use x-,x+,y-,y+,z-,z+ and positive Cartesian flux orientation.
struct ColdPressureFace {
  double projection{};          // metric * interpolated dt
  double owner_lower_weight{};  // COAST w at the lower-index cell
  double density_flux{};        // interpolated drho/dp * U_mid dot area
  double limiter_diffusion{};   // metric * pressure TVD gamma
  double mass_flux{};
  bool solid_neighbour{};
  bool prescribed_mass_flux{};  // fixed IBM inlet flux remains in the RHS
};

struct ColdPressureRow {
  std::array<double, 6> neighbour{};
  double diagonal{};
  double rhs{};
};

inline bool assemble_cold_pressure_row(
    double volume, double dt, double density_pressure_derivative,
    double density_rate, const std::array<ColdPressureFace, 6>& faces,
    bool fluid, ColdPressureRow& out) noexcept {
  if (!std::isfinite(volume) || volume <= 0 || !std::isfinite(dt) || dt <= 0)
    return false;
  ColdPressureRow row;
  if (!fluid) {
    row.diagonal = 1;
    out = row;
    return true;
  }
  if (!std::isfinite(density_pressure_derivative) ||
      density_pressure_derivative < 0 || !std::isfinite(density_rate))
    return false;
  row.diagonal = density_pressure_derivative / dt;
  row.rhs = -density_rate;
  for (unsigned i = 0; i < 6; ++i) {
    const auto& f = faces[i];
    const double w = f.owner_lower_weight;
    if (!std::isfinite(f.projection) || f.projection < 0 || !std::isfinite(w) ||
        w < 0 || w > 1 || !std::isfinite(f.density_flux) ||
        !std::isfinite(f.mass_flux) || !std::isfinite(f.limiter_diffusion) ||
        f.limiter_diffusion < 0)
      return false;
    const bool lower = i % 2 == 0;
    const double sign = lower ? -1 : 1;
    const double neighbour_advection = lower ? w : -(1 - w);
    const double diagonal_advection = lower ? -(1 - w) : w;
    row.neighbour[i] =
        (f.projection +
         0.5 * (f.limiter_diffusion + neighbour_advection * f.density_flux)) /
        volume;
    row.diagonal +=
        (f.projection +
         0.5 * (f.limiter_diffusion + diagonal_advection * f.density_flux)) /
        volume;
    row.rhs -= sign * f.mass_flux / volume;
  }
  // COAST applies the pressure IBM coefficient hook after press/compress/step.
  // Its scalar-flux hook must already have supplied zero blocked-face fluxes.
  for (unsigned i = 0; i < 6; ++i)
    if (faces[i].solid_neighbour) {
      if ((faces[i].mass_flux != 0 && !faces[i].prescribed_mass_flux) ||
          faces[i].density_flux != 0)
        return false;
      row.diagonal -= row.neighbour[i];
      row.neighbour[i] = 0;
    }
  if (!std::isfinite(row.diagonal) || row.diagonal <= 0 ||
      !std::isfinite(row.rhs))
    return false;
  for (double value : row.neighbour)
    if (!std::isfinite(value)) return false;
  out = row;
  return true;
}

inline double cold_pressure_flux_correction(double mass_flux, double projection,
                                            double lower_dp,
                                            double upper_dp) noexcept {
  return mass_flux - projection * (upper_dp - lower_dp);
}

// Keep density and drhodt relaxation identical. COAST update.F90 uses 0.5;
// a full fixed-temperature pressure correction uses 1 and preserves the
// ideal-gas EOS when derivative=rho/p. The caller selects the stage convention.
inline bool correct_cold_pressure_density(double density, double density_rate,
                                          double derivative, double dp,
                                          double dt, double relaxation,
                                          double& corrected_density,
                                          double& corrected_rate) noexcept {
  if (!std::isfinite(density) || density <= 0 || !std::isfinite(density_rate) ||
      !std::isfinite(derivative) || derivative < 0 || !std::isfinite(dp) ||
      !std::isfinite(dt) || dt <= 0 || !std::isfinite(relaxation) ||
      relaxation <= 0 || relaxation > 1)
    return false;
  const double delta = relaxation * derivative * dp;
  const double rho = density + delta;
  const double rate = density_rate + delta / dt;
  if (!std::isfinite(rho) || rho <= 0 || !std::isfinite(rate)) return false;
  corrected_density = rho;
  corrected_rate = rate;
  return true;
}

// COAST cmod followed by step: centre only the spatial matrix action.
// The source RHS is retained in full. storage is rho for per-volume rows,
// or rho*V for integrated rows; the caller supplies the intended time level.
inline bool time_centre_cold_row(const ColdPressureRow& spatial,
                                 double old_value,
                                 const std::array<double, 6>& old_neighbours,
                                 double storage, double dt,
                                 ColdPressureRow& output) noexcept {
  if (!std::isfinite(old_value) || !std::isfinite(storage) || storage <= 0 ||
      !std::isfinite(dt) || dt <= 0 || !std::isfinite(spatial.diagonal) ||
      !std::isfinite(spatial.rhs))
    return false;
  auto row = spatial;
  row.diagonal *= 0.5;
  for (unsigned f = 0; f < 6; ++f) {
    if (!std::isfinite(row.neighbour[f]) || !std::isfinite(old_neighbours[f]))
      return false;
    row.neighbour[f] *= 0.5;
    row.rhs += row.neighbour[f] * old_neighbours[f];
  }
  row.rhs -= row.diagonal * old_value;
  const double transient = storage / dt;
  row.diagonal += transient;
  row.rhs += transient * old_value;
  if (!std::isfinite(row.diagonal) || row.diagonal <= 0 ||
      !std::isfinite(row.rhs))
    return false;
  output = row;
  return true;
}

struct ColdTransportFace {
  double diffusion{}, lower_weight{}, mass_flux{};
};
// COAST condif on an orthogonal grid, TVD disabled (momentum reference).
// This is the advective row; source and boundary elimination follow separately.
inline bool assemble_cold_transport_row(
    double volume, const std::array<ColdTransportFace, 6>& faces, double source,
    ColdPressureRow& output) noexcept {
  if (!std::isfinite(volume) || volume <= 0 || !std::isfinite(source))
    return false;
  ColdPressureRow row;
  row.rhs = source;
  for (unsigned f = 0; f < 6; ++f) {
    const auto& face = faces[f];
    if (!std::isfinite(face.diffusion) || face.diffusion < 0 ||
        !std::isfinite(face.lower_weight) || face.lower_weight < 0 ||
        face.lower_weight > 1 || !std::isfinite(face.mass_flux))
      return false;
    const double weight = f % 2 ? face.lower_weight - 1 : face.lower_weight;
    row.neighbour[f] = (face.diffusion + weight * face.mass_flux) / volume;
    row.diagonal += row.neighbour[f];
  }
  if (!std::isfinite(row.diagonal)) return false;
  output = row;
  return true;
}

// COAST PDF path: step without cmod. An advective BE row uses old density
// for exact equivalence with conservative transport and BE continuity.
inline bool add_cold_backward_euler_storage(const ColdPressureRow& spatial,
                                            double old_value,
                                            double old_storage, double dt,
                                            ColdPressureRow& output) noexcept {
  if (!std::isfinite(old_value) || !std::isfinite(old_storage) ||
      old_storage <= 0 || !std::isfinite(dt) || dt <= 0)
    return false;
  auto row = spatial;
  const double storage = old_storage / dt;
  row.diagonal += storage;
  row.rhs += storage * old_value;
  for (double c : row.neighbour)
    if (!std::isfinite(c)) return false;
  if (!std::isfinite(row.diagonal) || row.diagonal <= 0 ||
      !std::isfinite(row.rhs))
    return false;
  output = row;
  return true;
}

struct ColdGridReport {
  std::uint64_t fluid{}, solid{}, cut_faces{}, inlet_faces{};
  double constant_error{}, signed_inlet_mass{};
};

inline Status eliminate_cold_boundaries(
    const PressureCorrectionBoundaryPlan& boundary, Int3 cells,
    std::vector<ColdPressureRow>& rows) {
  std::size_t i{};
  for (int z = 0; z < cells.z; ++z)
    for (int y = 0; y < cells.y; ++y)
      for (int x = 0; x < cells.x; ++x, ++i) {
        const Int3 c{x, y, z};
        const int coord[]{x, y, z}, extent[]{cells.x, cells.y, cells.z};
        for (unsigned a = 0; a < 3; ++a)
          for (unsigned side = 0; side < 2; ++side) {
            if (coord[a] != (side ? extent[a] - 1 : 0)) continue;
            Int3 face = c;
            if (side) {
              if (a == 0)
                ++face.x;
              else if (a == 1)
                ++face.y;
              else
                ++face.z;
            }
            PressureCorrectionFaceRule rule;
            auto status =
                boundary.face_rule(static_cast<CartesianAxis>(a), face, rule);
            if (!status) return status;
            if (!rule.is_nonperiodic_boundary()) continue;
            const double coefficient = rows[i].neighbour[2 * a + side];
            if (rule.kind == PressureCorrectionFaceKind::homogeneous_neumann)
              rows[i].diagonal -= coefficient;
            else if (rule.kind ==
                     PressureCorrectionFaceKind::homogeneous_dirichlet)
              rows[i].diagonal += coefficient;
            else
              return {StatusCode::invalid_plan, 17801};
            rows[i].neighbour[2 * a + side] = 0;
          }
        if (!std::isfinite(rows[i].diagonal) || rows[i].diagonal <= 0)
          return {StatusCode::numerical_failure, 17802};
      }
  return {};
}

// The exact pressure rows use residual per volume. Native MG stores an
// integrated finite-volume matrix, so apply its inverse to V*r.
class VolumeScaledPreconditioner final : public LinearPreconditioner {
 public:
  VolumeScaledPreconditioner(NativeCartesianMgPlan& mg,
                            const CartesianKernelPlan& kernels,
                            FieldView scratch)
      : mg_(mg), kernels_(kernels), scratch_(scratch) {}
  LinearPreconditionerCertificate certificate() const noexcept override {
    return mg_.certificate();
  }
  Status prepare_batch(const LinearPreconditionerBatchDescriptor& descriptor,
                       LinearPreconditionerBatchTicket& ticket) noexcept override {
    return mg_.prepare_batch(descriptor, ticket);
  }
  Status apply_prepared(ConstFieldView residual, FieldView correction,
                        std::uint32_t iteration,
                        const LinearPreconditionerBatchTicket& ticket) noexcept override {
    scale_residual(residual);
    return mg_.apply_prepared(as_const(scratch_), correction, iteration, ticket);
  }
  Status apply(ConstFieldView residual, FieldView correction,
               std::uint32_t iteration) noexcept override {
    scale_residual(residual);
    return mg_.apply(as_const(scratch_), correction, iteration);
  }
 private:
  void scale_residual(ConstFieldView residual) noexcept {
    const auto cells = residual.interior;
    for (int z=0; z<cells.z; ++z) for (int y=0; y<cells.y; ++y)
      for (int x=0; x<cells.x; ++x) {
        const Int3 c{x,y,z};
        scratch_.unchecked(c,0) = residual.unchecked(c,0)*cell_volume(kernels_,c);
      }
  }

  NativeCartesianMgPlan& mg_;
  const CartesianKernelPlan& kernels_;
  FieldView scratch_;
};

class ColdPressureOperator final : public LinearOperator {
 public:
  ColdPressureOperator(const std::vector<ColdPressureRow>& rows, Int3 cells,
                       HaloEngine& halo, LinearIdentity identity)
      : rows_(rows), cells_(cells), halo_(halo), identity_(identity) {}
  LinearOperatorCertificate certificate() const noexcept override {
    return {identity_, 0x434f4c44504f5031ULL, cells_,
            LinearOperatorClass::nonsymmetric};
  }
  LinearOperatorFailureProvenance failure_provenance() const noexcept override {
    return failure_;
  }
  Status apply(FieldView input, FieldView output) const noexcept override {
    failure_ = {};
    HaloTicket ticket;
    auto status = halo_.begin(140U, {&input, 1U}, ticket);
    if (status) status = halo_.finish(ticket, {&input, 1U});
    if (!status) {
      failure_ = {status, LinearOperatorStatusScope::collective,
                  halo_.lowest_failing_rank()};
      return status;
    }
    std::size_t i{};
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x, ++i) {
          const auto& row = rows_[i];
          const Int3 c{x, y, z};
          const Int3 neighbour[] = {{x - 1, y, z}, {x + 1, y, z},
                                    {x, y - 1, z}, {x, y + 1, z},
                                    {x, y, z - 1}, {x, y, z + 1}};
          double value = row.diagonal * input.unchecked(c, 0);
          for (unsigned f = 0; f < 6; ++f)
            if (row.neighbour[f] != 0)
              value -= row.neighbour[f] * input.unchecked(neighbour[f], 0);
          output.unchecked(c, 0) = value;
        }
    return {};
  }

 private:
  const std::vector<ColdPressureRow>& rows_;
  Int3 cells_;
  HaloEngine& halo_;
  LinearIdentity identity_;
  mutable LinearOperatorFailureProvenance failure_{};
};

// Rank-local DILU: COAST uses a local incomplete-factorization pressure
// preconditioner. Retain opposite-direction coefficient products here since
// compress and volume normalization make the C++ pressure rows nonsymmetric.
class ColdPressureDilu final : public LinearPreconditioner {
 public:
  std::uint64_t owned_payload_bytes() const noexcept {
    return inverse_.capacity() * sizeof(double);
  }
  ColdPressureDilu(const std::vector<ColdPressureRow>& rows, Int3 cells,
                   LinearIdentity identity)
      : rows_(rows),
        cells_(cells),
        identity_(identity),
        inverse_(rows.size()) {}
  Status prepare() {
    const std::size_t stride[]{1, std::size_t(cells_.x),
                               std::size_t(cells_.x) * cells_.y};
    std::size_t i{};
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x, ++i) {
          const int c[]{x, y, z};
          double diagonal = rows_[i].diagonal;
          for (unsigned a = 0; a < 3; ++a)
            if (c[a] > 0)
              diagonal -= rows_[i].neighbour[2 * a] *
                          rows_[i - stride[a]].neighbour[2 * a + 1] *
                          inverse_[i - stride[a]];
          if (!std::isfinite(diagonal) || diagonal <= 0)
            return {StatusCode::numerical_failure, 17803};
          inverse_[i] = 1 / diagonal;
        }
    return {};
  }
  LinearPreconditionerCertificate certificate() const noexcept override {
    return {identity_, 0x434f4c4444494c55ULL,
            LinearPreconditionerClass::fixed_general};
  }
  Status apply(ConstFieldView input, FieldView output,
               std::uint32_t) noexcept override {
    std::size_t i{};
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x, ++i) {
          const auto& r = rows_[i];
          double v = input.unchecked({x, y, z}, 0);
          if (x > 0) v += r.neighbour[0] * output.unchecked({x - 1, y, z}, 0);
          if (y > 0) v += r.neighbour[2] * output.unchecked({x, y - 1, z}, 0);
          if (z > 0) v += r.neighbour[4] * output.unchecked({x, y, z - 1}, 0);
          output.unchecked({x, y, z}, 0) = v * inverse_[i];
        }
    i = rows_.size();
    for (int z = cells_.z - 1; z >= 0; --z)
      for (int y = cells_.y - 1; y >= 0; --y)
        for (int x = cells_.x - 1; x >= 0; --x) {
          --i;
          const auto& r = rows_[i];
          double v{};
          if (x + 1 < cells_.x)
            v += r.neighbour[1] * output.unchecked({x + 1, y, z}, 0);
          if (y + 1 < cells_.y)
            v += r.neighbour[3] * output.unchecked({x, y + 1, z}, 0);
          if (z + 1 < cells_.z)
            v += r.neighbour[5] * output.unchecked({x, y, z + 1}, 0);
          output.unchecked({x, y, z}, 0) += v * inverse_[i];
        }
    return {};
  }

 private:
  const std::vector<ColdPressureRow>& rows_;
  Int3 cells_;
  LinearIdentity identity_;
  std::vector<double> inverse_;
};

// COAST enables momentum VLS above its isothermal Mach threshold.
// Freeze one globally agreed choice for the whole nonlinear time step.
inline bool cold_momentum_uses_tvd(double isothermal_mach) noexcept {
  return isothermal_mach > 0.60;
}

struct ColdMomentumClosureReport {
  std::uint64_t solid_rows{}, cut_links{}, physical_links{};
};
// Reuse the authoritative Hundun pressure, stress, IBM and fixed-face residual.
// Replace its BDF/convection splitting with a frozen advective CN matrix.
inline Status close_cold_momentum_rows(
    const CartesianKernelPlan& kernels, MeshPatch patch, Int3 global_cells,
    const EBTopology* topology, const BoundaryPlan& boundary,
    const EquationStateView& state, const EquationAssemblyContext& context,
    const EquationSystemView& reference,
    std::array<std::vector<ColdPressureRow>, 3>& rows,
    ColdMomentumClosureReport& report, ConvectionScheme reference_convection,
    bool coast_momentum_tvd) {
  const auto cells = patch.cells;
  const auto count = std::size_t(cells.x) * cells.y * cells.z;
  for (auto& r : rows) r.resize(count);
  const auto coord = [](Int3 c, unsigned a) {
    return a == 0 ? c.x : a == 1 ? c.y : c.z;
  };
  std::array<bool, 3> periodic{};
  for (unsigned axis = 0; axis < 3; ++axis) {
    const BoundaryFacePlan* rule{};
    const auto status = boundary.face(static_cast<CartesianFace>(2 * axis), rule);
    if (!status || !rule) return {StatusCode::invalid_plan, 17811};
    periodic[axis] = rule->periodic;
  }
  const auto fluid = [&](Int3 c) {
    Int3 g{c.x + patch.begin.x, c.y + patch.begin.y, c.z + patch.begin.z};
    for (unsigned axis = 0; axis < 3; ++axis) {
      auto& value = axis == 0 ? g.x : axis == 1 ? g.y : g.z;
      const auto extent = coord(global_cells, axis);
      if (value < 0 || value >= extent) {
        if (!periodic[axis]) return true;
        value = (value % extent + extent) % extent;
      }
    }
    return !topology || topology->is_fluid_global(g);
  };
  std::size_t index{};
  for (int z = 0; z < cells.z; ++z)
    for (int y = 0; y < cells.y; ++y)
      for (int x = 0; x < cells.x; ++x, ++index) {
        const Int3 c{x, y, z};
        if (!fluid(c)) {
          report.solid_rows += 3;
          for (auto& r : rows) {
            r[index] = {};
            r[index].diagonal = 1;
          }
          continue;
        }
        const double volume = cell_volume(kernels, c),
                     rho = state.density.trial.unchecked(c, 0),
                     old_rho = state.density.accepted.unchecked(c, 0);
        const Int3 nb[] = {{x - 1, y, z}, {x + 1, y, z}, {x, y - 1, z},
                           {x, y + 1, z}, {x, y, z - 1}, {x, y, z + 1}};
        std::array<ColdTransportFace, 6> faces;
        double div{}, outgoing{}, diffusion_sum{};
        for (unsigned f = 0; f < 6; ++f) {
          unsigned a = f / 2;
          auto axis = static_cast<CartesianAxis>(a);
          const Int3 lower = f % 2 ? c : nb[f], upper = f % 2 ? nb[f] : c;
          double ds = centre_coordinate(kernels, axis, coord(upper, a)) -
                      centre_coordinate(kernels, axis, coord(lower, a));
          double weight = (centre_coordinate(kernels, axis, coord(upper, a)) -
                           face_coordinate(kernels, axis, coord(upper, a))) /
                          ds;
          double diffusion = a == 0 ? reference.x_coefficient.unchecked(upper)
                             : a == 1
                                 ? reference.y_coefficient.unchecked(upper)
                                 : reference.z_coefficient.unchecked(upper);
          double flux = a == 0   ? context.mass_flux.x.unchecked(upper)
                        : a == 1 ? context.mass_flux.y.unchecked(upper)
                                 : context.mass_flux.z.unchecked(upper);
          faces[f] = {diffusion, weight, flux};
          diffusion_sum += diffusion;
          div += (f % 2 ? 1 : -1) * flux / volume;
          outgoing += std::max((f % 2 ? 1 : -1) * flux, 0.0);
          if (!fluid(nb[f])) ++report.cut_links;
        }
        ColdPressureRow base;
        if (!assemble_cold_transport_row(volume, faces, 0, base))
          return {StatusCode::invalid_plan, 17810};
        for (unsigned component = 0; component < 3; ++component) {
          ColdPressureRow spatial = base;
          double convection_correction{};
          if (reference_convection != ConvectionScheme::central2 || coast_momentum_tvd) {
            auto limited_faces = faces;
            for (unsigned f = 0; f < 6; ++f) {
              const unsigned a = f / 2;
              const auto axis = static_cast<CartesianAxis>(a);
              const Int3 lower = f % 2 ? c : nb[f];
              const Int3 upper = f % 2 ? nb[f] : c;
              const int global_face = coord(upper, a) + coord(patch.begin, a);
              // Keep prescribed physical/IBM source faces in the native
              // boundary operator. Periodic faces share the interior policy.
              if (!fluid(lower) || !fluid(upper) ||
                  (!periodic[a] && (global_face <= 0 ||
                   global_face >= coord(global_cells, a)))) continue;
              const auto offset = [a](Int3 p, int d) {
                (a == 0 ? p.x : a == 1 ? p.y : p.z) += d;
                return p;
              };
              const Int3 before = offset(lower, -1), after = offset(upper, 1);
              const auto value_at = [&](Int3 p) {
                return state.velocity.trial.unchecked(p, component);
              };
              const auto gradient = [&](Int3 lo, Int3 hi) {
                return (value_at(hi) - value_at(lo)) /
                    (centre_coordinate(kernels, axis, coord(hi, a)) -
                     centre_coordinate(kernels, axis, coord(lo, a)));
              };
              const double weight = faces[f].lower_weight;
              const double flux = faces[f].mass_flux;
              const double limiter = coast_momentum_tvd
                  ? cold_pressure_limiter_diffusion(
                        gradient(before, lower), gradient(lower, upper),
                        gradient(upper, after), weight, flux,
                        fluid(lower) && fluid(upper) && fluid(after),
                        1e-30, std::numeric_limits<double>::epsilon())
                  : 0.0;
              const double extra = std::max(0.0, limiter - faces[f].diffusion);
              limited_faces[f].diffusion += extra;
              double native_face{};
              auto status = reconstruct_cartesian_convection_face(
                  kernels, reference_convection, state.velocity.trial,
                  component, axis, upper, flux, native_face);
              if (!status) return status;
              const double central = weight * value_at(lower) +
                                     (1 - weight) * value_at(upper);
              convection_correction += (f % 2 ? 1 : -1) *
                  (flux * (central - native_face) -
                   extra * (value_at(upper) - value_at(lower))) / volume;
            }
            if (!assemble_cold_transport_row(volume, limited_faces, 0, spatial))
              return {StatusCode::invalid_plan, 17810};
          }
          const double value = state.velocity.trial.unchecked(c, component);
          const double unsteady =
              context.bdf.a0 * rho * value +
              context.bdf.a1 * old_rho *
                  state.velocity.accepted.unchecked(c, component) +
              context.bdf.a2 * state.density.previous.unchecked(c, 0) *
                  state.velocity.previous.unchecked(c, component);
          const double wall_source_diagonal =
              reference.diagonal.unchecked(c, component) -
              context.bdf.a0 * rho * volume - diffusion_sum - outgoing;
          spatial.diagonal += wall_source_diagonal / volume;
          for (unsigned f = 0; f < 6; ++f)
            if (!fluid(nb[f])) spatial.neighbour[f] = 0;
          double action = spatial.diagonal * value;
          for (unsigned f = 0; f < 6; ++f)
            if (spatial.neighbour[f])
              action -= spatial.neighbour[f] *
                        state.velocity.trial.unchecked(nb[f], component);
          const double advective_residual =
              reference.residual.unchecked(c, component) / volume - unsteady -
              value * div + convection_correction;
          spatial.rhs = action - advective_residual;
          for (unsigned f = 0; f < 6; ++f) {
            unsigned a = f / 2;
            Int3 face = f % 2 ? nb[f] : c;
            int global_face = coord(face, a) + coord(patch.begin, a);
            if (global_face != 0 && global_face != coord(global_cells, a))
              continue;
            const BoundaryFacePlan* rule{};
            auto status = boundary.face(static_cast<CartesianFace>(f), rule);
            if (!status || !rule) return {StatusCode::invalid_plan, 17811};
            if (rule->periodic) continue;
            bool dirichlet{};
            switch (rule->flow_kind) {
              case BoundaryKind::symmetry:
                dirichlet = component == a;
                break;
              case BoundaryKind::no_slip_wall:
              case BoundaryKind::velocity_inlet:
              case BoundaryKind::mass_flow_inlet:
                dirichlet = true;
                break;
              case BoundaryKind::pressure_outlet:
                dirichlet = component != a &&
                    boundary.allow_backflow().data[rule->flow_parameter] != 0 &&
                    (f % 2 ? 1.0 : -1.0) *
                        state.velocity.trial.unchecked(c, a) < 0.0;
                break;
              case BoundaryKind::zero_gradient_mass_outlet:
                dirichlet = false;
                break;
              default:
                return {StatusCode::invalid_plan, 17812};
            }
            const double coefficient = spatial.neighbour[f];
            if (dirichlet) {
              const double boundary_value =
                  0.5 *
                  (value + state.velocity.trial.unchecked(nb[f], component));
              spatial.diagonal += coefficient;
              spatial.rhs += 2 * coefficient * boundary_value;
            } else
              spatial.diagonal -= coefficient;
            spatial.neighbour[f] = 0;
            ++report.physical_links;
          }
          std::array<double, 6> old;
          for (unsigned f = 0; f < 6; ++f)
            old[f] = state.velocity.accepted.unchecked(nb[f], component);
          if (!time_centre_cold_row(
                  spatial, state.velocity.accepted.unchecked(c, component), old,
                  0.5 * (rho + old_rho), context.dt, rows[component][index]))
            return {StatusCode::numerical_failure, 17813};
        }
      }
  return {};
}

// Upwind implicit splitting; the retained full residual supplies the higher
// order convection, temperature conduction, pressure and viscous work.
inline Status close_cold_enthalpy_rows(
    const CartesianKernelPlan& kernels, MeshPatch patch, Int3 global_cells,
    const EBTopology* topology, const BoundaryPlan& boundary,
    const EquationStateView& state, const EquationAssemblyContext& context,
    const EquationSystemView& reference, std::vector<ColdPressureRow>& rows) {
  auto cells = patch.cells;
  rows.resize(std::size_t(cells.x) * cells.y * cells.z);
  const auto coord = [](Int3 c, unsigned a) {
    return a == 0 ? c.x : a == 1 ? c.y : c.z;
  };
  const auto fluid = [&](Int3 c) {
    Int3 g{c.x + patch.begin.x, c.y + patch.begin.y, c.z + patch.begin.z};
    if (g.x < 0 || g.y < 0 || g.z < 0 || g.x >= global_cells.x ||
        g.y >= global_cells.y || g.z >= global_cells.z)
      return true;
    return !topology || topology->is_fluid_global(g);
  };
  std::size_t index{};
  for (int z = 0; z < cells.z; ++z)
    for (int y = 0; y < cells.y; ++y)
      for (int x = 0; x < cells.x; ++x, ++index) {
        Int3 c{x, y, z};
        if (!fluid(c)) {
          rows[index] = {};
          rows[index].diagonal = 1;
          rows[index].rhs = state.enthalpy.trial.unchecked(c, 0);
          continue;
        }
        double volume = cell_volume(kernels, c),
               rho = state.density.trial.unchecked(c, 0),
               rho_old = state.density.accepted.unchecked(c, 0),
               h = state.enthalpy.trial.unchecked(c, 0),
               h_old = state.enthalpy.accepted.unchecked(c, 0);
        const Int3 nb[] = {{x - 1, y, z}, {x + 1, y, z}, {x, y - 1, z},
                           {x, y + 1, z}, {x, y, z - 1}, {x, y, z + 1}};
        std::array<ColdTransportFace, 6> faces;
        std::array<double, 6> physical_diffusion;
        double div{}, diffusion_sum{};
        for (unsigned f = 0; f < 6; ++f) {
          unsigned a = f / 2;
          auto axis = static_cast<CartesianAxis>(a);
          Int3 lower = f % 2 ? c : nb[f], upper = f % 2 ? nb[f] : c;
          double ds = centre_coordinate(kernels, axis, coord(upper, a)) -
                      centre_coordinate(kernels, axis, coord(lower, a));
          double weight = (centre_coordinate(kernels, axis, coord(upper, a)) -
                           face_coordinate(kernels, axis, coord(upper, a))) /
                          ds;
          double diffusion = a == 0 ? reference.x_coefficient.unchecked(upper)
                             : a == 1
                                 ? reference.y_coefficient.unchecked(upper)
                                 : reference.z_coefficient.unchecked(upper);
          double flux = a == 0   ? context.mass_flux.x.unchecked(upper)
                        : a == 1 ? context.mass_flux.y.unchecked(upper)
                                 : context.mass_flux.z.unchecked(upper);
          physical_diffusion[f] = diffusion;
          diffusion_sum += diffusion;
          div += (f % 2 ? 1 : -1) * flux / volume;
          double upwind_diffusion =
              (flux >= 0 ? 1 - weight : weight) * std::abs(flux);
          faces[f] = {diffusion + upwind_diffusion, weight, flux};
        }
        ColdPressureRow spatial;
        if (!assemble_cold_transport_row(volume, faces, 0, spatial))
          return {StatusCode::invalid_plan, 17814};
        spatial.diagonal += (reference.diagonal.unchecked(c, 0) -
                             context.bdf.a0 * rho * volume - diffusion_sum) /
                            volume;
        for (unsigned f = 0; f < 6; ++f) {
          if (!fluid(nb[f])) {
            spatial.diagonal -= physical_diffusion[f] / volume;
            spatial.neighbour[f] = 0;
            continue;
          }
          unsigned a = f / 2;
          Int3 face = f % 2 ? nb[f] : c;
          int global_face = coord(face, a) + coord(patch.begin, a);
          if (global_face != 0 && global_face != coord(global_cells, a))
            continue;
          const BoundaryFacePlan* rule{};
          auto status = boundary.face(static_cast<CartesianFace>(f), rule);
          if (!status || !rule) return {StatusCode::invalid_plan, 17815};
          if (rule->periodic) continue;
          bool dirichlet = rule->flow_kind == BoundaryKind::velocity_inlet ||
                           rule->flow_kind == BoundaryKind::mass_flow_inlet ||
                           rule->thermal_kind == BoundaryKind::isothermal_wall ||
                           (rule->flow_kind == BoundaryKind::pressure_outlet &&
                            boundary.allow_backflow().data[rule->flow_parameter] != 0 &&
                            (f % 2 ? 1.0 : -1.0) *
                                state.velocity.trial.unchecked(c, a) < 0.0);
          if (dirichlet)
            spatial.diagonal += physical_diffusion[f] / volume;
          else
            spatial.diagonal -= spatial.neighbour[f];
          spatial.neighbour[f] = 0;
        }
        double action = spatial.diagonal * h;
        for (unsigned f = 0; f < 6; ++f)
          if (spatial.neighbour[f])
            action -=
                spatial.neighbour[f] * state.enthalpy.trial.unchecked(nb[f], 0);
        const double unsteady = (rho * h - rho_old * h_old) / context.dt;
        spatial.rhs = action - (reference.residual.unchecked(c, 0) / volume -
                                unsteady - h * div);
        if (!add_cold_backward_euler_storage(spatial, h_old, rho_old,
                                             context.dt, rows[index]))
          return {StatusCode::numerical_failure, 17816};
      }
  return {};
}

inline bool assemble_midpoint_cold_grid(
    const CartesianKernelPlan& kernels, MeshPatch patch, Int3 global_cells,
    const EBTopology* topology, const BoundaryPlan& boundary,
    ConstFieldView rho, ConstFieldView old_rho, ConstFieldView velocity,
    ConstFieldView old_velocity, ConstFieldView pressure,
    double reference_pressure, double dt, ConstFaceFluxView prescribed_flux,
    std::vector<ColdPressureRow>& rows, ColdGridReport& report,
    std::vector<std::array<ColdPressureFace, 6>>* saved_faces = nullptr) {
  const auto cells = patch.cells;
  if (!valid_cell_view(rho, cells, 0, 1, 1) ||
      !valid_cell_view(old_rho, cells, 0, 1, 1) ||
      !valid_cell_view(velocity, cells, 0, 3, 1) ||
      !valid_cell_view(old_velocity, cells, 0, 3, 1) ||
      !valid_cell_view(pressure, cells, 0, 1, 2))
    return false;
  rows.resize(std::size_t(cells.x) * cells.y * cells.z);
  if (saved_faces) saved_faces->resize(rows.size());
  const auto index = [](Int3 c, unsigned axis) {
    return axis == 0 ? c.x : axis == 1 ? c.y : c.z;
  };
  const auto shift = [](Int3 c, unsigned axis, int amount) {
    if (axis == 0)
      c.x += amount;
    else if (axis == 1)
      c.y += amount;
    else
      c.z += amount;
    return c;
  };
  const auto fluid = [&](Int3 c) {
    const Int3 g{c.x + patch.begin.x, c.y + patch.begin.y, c.z + patch.begin.z};
    if (g.x < 0 || g.y < 0 || g.z < 0 || g.x >= global_cells.x ||
        g.y >= global_cells.y || g.z >= global_cells.z)
      return true;  // external boundary ghosts are handled separately from IBM
    return !topology || topology->is_fluid_global(g);
  };
  const auto derivative = [&](Int3 c) {
    return rho.unchecked(c, 0) /
           (reference_pressure + pressure.unchecked(c, 0));
  };
  const auto midpoint_velocity = [&](Int3 c) {
    std::array<double, 3> value{};
    for (unsigned a = 0; a < 3; ++a)
      value[a] =
          0.5 * (velocity.unchecked(c, a) + old_velocity.unchecked(c, a));
    return value;
  };
  const auto midpoint_rho = [&](Int3 c) {
    return 0.5 * (rho.unchecked(c, 0) + old_rho.unchecked(c, 0));
  };
  std::size_t ordinal{};
  for (int z = 0; z < cells.z; ++z)
    for (int y = 0; y < cells.y; ++y)
      for (int x = 0; x < cells.x; ++x, ++ordinal) {
        const Int3 cell{x, y, z};
        std::array<ColdPressureFace, 6> faces{};
        const bool active = fluid(cell);
        double density_divergence{};
        if (active) {
          ++report.fluid;
          for (unsigned a = 0; a < 3; ++a)
            for (unsigned side = 0; side < 2; ++side) {
              const auto axis = static_cast<CartesianAxis>(a);
              const Int3 lo = shift(cell, a, side == 0 ? -1 : 0),
                         hi = shift(lo, a, 1), face = hi;
              auto& f = faces[2 * a + side];
              const double distance =
                  centre_coordinate(kernels, axis, index(hi, a)) -
                  centre_coordinate(kernels, axis, index(lo, a));
              f.owner_lower_weight =
                  (centre_coordinate(kernels, axis, index(hi, a)) -
                   face_coordinate(kernels, axis, index(face, a))) /
                  distance;
              const double area = face_area(kernels, axis, face);
              f.projection = .5 * dt * area / distance;
              f.solid_neighbour = !fluid(side == 0 ? lo : hi);
              if (f.solid_neighbour) {
                ++report.cut_faces;
                f.mass_flux = a == 0   ? prescribed_flux.x.unchecked(face)
                              : a == 1 ? prescribed_flux.y.unchecked(face)
                                       : prescribed_flux.z.unchecked(face);
                if (f.mass_flux != 0) {
                  f.prescribed_mass_flux = true;
                  ++report.inlet_faces;
                  report.signed_inlet_mass += (side ? 1 : -1) * f.mass_flux;
                }
                continue;
              }
              const int global_face = index(face, a) + index(patch.begin, a);
              if (global_face == 0 || global_face == index(global_cells, a)) {
                const BoundaryFacePlan* rule{};
                if (!boundary.face(
                        static_cast<CartesianFace>(2 * a + (global_face != 0)),
                        rule) ||
                    !rule)
                  return false;
                if (rule->flow_kind == BoundaryKind::pressure_outlet) {
                  // The thermophysical closure stores the physical outlet
                  // material in the exterior slot. At fixed p and frozen h/Y
                  // its density has zero pressure-correction response.
                  const Int3 exterior = global_face == 0 ? lo : hi;
                  const auto ul = midpoint_velocity(lo), uh = midpoint_velocity(hi);
                  f.mass_flux = midpoint_rho(exterior) * area *
                      (f.owner_lower_weight * ul[a] +
                       (1 - f.owner_lower_weight) * uh[a]);
                  continue;
                }
                if (!rule->periodic &&
                    rule->flow_kind !=
                        BoundaryKind::zero_gradient_mass_outlet &&
                    rule->flow_kind != BoundaryKind::pressure_outlet) {
                  f.prescribed_mass_flux = true;
                  f.mass_flux = a == 0   ? prescribed_flux.x.unchecked(face)
                                : a == 1 ? prescribed_flux.y.unchecked(face)
                                         : prescribed_flux.z.unchecked(face);
                  continue;
                }
              }
              std::array<double, 3> oriented_area{};
              oriented_area[a] = area;
              const auto ul = midpoint_velocity(lo), uh = midpoint_velocity(hi);
              f.mass_flux = cold_coefficient_face_flux(
                  midpoint_rho(lo), midpoint_rho(hi), f.owner_lower_weight, ul,
                  uh, oriented_area, false);
              f.density_flux = cold_coefficient_face_flux(
                  derivative(lo), derivative(hi), f.owner_lower_weight, ul, uh,
                  oriented_area, false);
              const auto gradient = [&](Int3 p0, Int3 p1) {
                const double ds =
                    centre_coordinate(kernels, axis, index(p1, a)) -
                    centre_coordinate(kernels, axis, index(p0, a));
                return (pressure.unchecked(p1, 0) - pressure.unchecked(p0, 0)) /
                       ds;
              };
              const Int3 before = shift(lo, a, -1), after = shift(hi, a, 1);
              // Collocated pressure/velocity coupling (d2pds on an orthogonal
              // mesh). Centre gradients use fluid values only at IBM cuts.
              // The midpoint mass flux responds with dt/2 to an endpoint
              // pressure change, matching f.projection in the correction.
              const double face_gradient = gradient(lo, hi);
              const double lower_gradient = fluid(before)
                  ? gradient(before, hi) : face_gradient;
              const double upper_gradient = fluid(after)
                  ? gradient(lo, after) : face_gradient;
              f.mass_flux += 0.5 * dt * area *
                  (0.5 * (lower_gradient + upper_gradient) - face_gradient);
              // Match vls: the positive-flow branch tests lo,hi,after activity.
              f.limiter_diffusion = cold_pressure_limiter_diffusion(
                  gradient(before, lo), gradient(lo, hi), gradient(hi, after),
                  f.owner_lower_weight, f.density_flux,
                  fluid(lo) && fluid(hi) && fluid(after), 1e-30,
                  std::numeric_limits<double>::epsilon());
              density_divergence += (side ? 1 : -1) * f.density_flux;
            }
        } else
          ++report.solid;
        const double volume = cell_volume(kernels, cell);
        const double d = active ? derivative(cell) : 0;
        const double rate =
            active ? (rho.unchecked(cell, 0) - old_rho.unchecked(cell, 0)) / dt
                   : 0;
        if (!assemble_cold_pressure_row(volume, dt, d, rate, faces, active,
                                        rows[ordinal]))
          return false;
        if (saved_faces) (*saved_faces)[ordinal] = faces;
        if (active) {
          const auto& row = rows[ordinal];
          double response = row.diagonal, scale = std::abs(row.diagonal);
          for (double q : row.neighbour) {
            response -= q;
            scale += std::abs(q);
          }
          const double expected = d / dt + 0.5 * density_divergence / volume;
          report.constant_error = std::max(
              report.constant_error, std::abs(response - expected) /
                                         std::max(scale, std::abs(expected)));
        }
      }
  return true;
}

}  // namespace hundun::v04::detail
