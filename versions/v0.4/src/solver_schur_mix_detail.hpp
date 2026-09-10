// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#pragma once
#include "hundun/v04_flow.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_conservative_energy_detail.hpp"
#include <algorithm>
#include <array>
#include <vector>
namespace hundun::v04 {
namespace detail {
// Frozen internal-fluid-face mixture convection response. Storage is reserved
// before stepping; each binding replaces every face weight, including zeros.
class SchurMixture {
public:
  FieldView eliminated;
  LinearOperatorCertificate cert;
  PlanFingerprint target_fingerprint{};
  PressureContinuityActivityView activity;
  struct Face {
    std::array<double, 4> h{};
    double p{};
  };
  std::array<std::vector<Face>, 3> faces;
  static Int3 offset(Int3 c, int a, int d) {
    (a == 0 ? c.x : a == 1 ? c.y : c.z) += d;
    return c;
  }
  std::size_t index(int a, Int3 c) const {
    auto n = offset(cert.local_shape, a, 1);
    return c.x + std::size_t(n.x) * (c.y + std::size_t(n.y) * c.z);
  }
  template <class F> void cells(F f) const {
    auto n = cert.local_shape;
    std::size_t i = 0;
    for (int z = 0; z < n.z; ++z)
      for (int y = 0; y < n.y; ++y)
        for (int x = 0; x < n.x; ++x, ++i)
          if (!activity.cells.size || activity.cells.data[i])
            f(Int3{x, y, z});
  }
  std::vector<double> pure;
  void reserve(Int3 shape, std::size_t species) {
    cert.local_shape = shape;
    pure.resize(species);
    for (int a = 0; a < 3; ++a) {
      auto n = offset(shape, a, 1);
      faces[a].resize(std::size_t(n.x) * n.y * n.z);
    }
  }
  static bool valid_face(ConstFaceFieldView face, Int3 cells,
                         int axis) noexcept {
    auto e = offset(cells, axis, 1);
    return face.base && face.axis == static_cast<CartesianAxis>(axis) &&
           same_cells(face.extents, e) && face.stride_y >= std::size_t(e.x) &&
           face.stride_z >= face.stride_y * std::size_t(e.y);
  }
  Status bind(const PressureCorrectionBoundaryPlan &boundary,
              const CartesianKernelPlan &kernels,
              const ThermodynamicsPlan &thermo, ConvectionScheme hs,
              ConvectionScheme ys, ConstFieldView h, ConstFieldView t,
              ConstFieldView cp, Span<const PrimitiveHistory> species,
              ConstFaceFluxView flux, std::array<ConstFaceFieldView, 3> coeff) {
    Status status;
    if (pure.size() != species.size || (species.size && !species.data))
      return {StatusCode::invalid_plan, 1583U};
    const auto shape = cert.local_shape;
    const auto hreach = hs == ConvectionScheme::central2 ? 1U : 2U;
    const auto yreach =
        std::max(hreach, ys == ConvectionScheme::central2 ? 1U : 2U);
    const auto count = std::size_t(shape.x) * shape.y * shape.z;
    if (!same_cells(shape, kernels.cells()) ||
        !valid_cell_view(h, shape, 0U, 1U, hreach) ||
        !valid_cell_view(t, shape, 0U, 1U, 1U) ||
        !valid_cell_view(cp, shape, 0U, 1U, 1U) ||
        (activity.cells.size &&
         (!activity.cells.data || activity.cells.size != count)))
      return {StatusCode::invalid_plan, 1583U};
    for (std::size_t i = 0; i < species.size; ++i)
      if (!valid_cell_view(species.data[i].trial, shape, 0U, 1U, yreach))
        return {StatusCode::invalid_plan, 1583U};
    const int first = hs == ConvectionScheme::central2 ? 1 : 0;
    const int last = hs == ConvectionScheme::central2 ? 3 : 4;
    target_fingerprint = thermo.fingerprint();
    const auto mix = [&](std::uint64_t value) {
      target_fingerprint ^= value;
      target_fingerprint *= UINT64_C(1099511628211);
    };
    mix(static_cast<unsigned>(hs));
    mix(static_cast<unsigned>(ys));
    for (auto field : {h, t, cp}) {
      mix(field.field);
      mix(field.revision);
    }
    for (std::size_t s = 0; s < species.size; ++s) {
      mix(species.data[s].trial.field);
      mix(species.data[s].trial.revision);
    }
    std::fill(pure.begin(), pure.end(), 0.0);
    std::array<ConstFaceFieldView, 3> phi{flux.x, flux.y, flux.z};
    std::array<Span<const std::uint8_t>, 3> active{
        activity.x_faces, activity.y_faces, activity.z_faces};
    for (int a = 0; a < 3; ++a) {
      const auto axis = static_cast<CartesianAxis>(a);
      auto n = offset(cert.local_shape, a, 1);
      if (!valid_face(phi[a], shape, a) || !valid_face(coeff[a], shape, a) ||
          (active[a].size &&
           (!active[a].data || active[a].size != faces[a].size())) ||
          faces[a].size() != std::size_t(n.x) * n.y * n.z)
        return {StatusCode::invalid_plan, 1583U};
      std::fill(faces[a].begin(), faces[a].end(), Face{});
      for (int z = 0; z < n.z; ++z)
        for (int y = 0; y < n.y; ++y)
          for (int x = 0; x < n.x; ++x) {
            Int3 f{x, y, z};
            auto i = index(a, f);
            if (active[a].size && !active[a].data[i])
              continue;
            PressureCorrectionFaceRule rule;
            status = boundary.face_rule(axis, f, rule);
            if (!status)
              return status;
            if (rule.is_nonperiodic_boundary())
              continue;
            const double mass = phi[a].unchecked(f);
            if (mass == 0)
              continue;
            auto left = offset(f, a, -1);
            int normal = a == 0 ? x : a == 1 ? y : z;
            const double tf = interpolate_face(
                kernels, axis, normal, t.unchecked(left, 0), t.unchecked(f, 0));
            double dummy, cpdep, r;
            status = thermo.mixture_enthalpy(tf, {pure.data(), pure.size()},
                                             dummy, cpdep, r);
            if (!status)
              return status;
            std::array<long double, 4> thermal{};
            std::array<double, 4> original{}, anchor{};
            for (int j = first; j < last; ++j)
              thermal[j] = original[j] = h.unchecked(offset(f, a, j - 2), 0);
            long double yf_h = 0;
            double yf_cp = 0;
            for (std::size_t s = 0; s < species.size; ++s) {
              double dh;
              status =
                  thermo.independent_species_enthalpy_difference(s, tf, dh);
              if (!status)
                return status;
              pure[s] = 1;
              double cps;
              status = thermo.mixture_enthalpy(tf, {pure.data(), pure.size()},
                                               dummy, cps, r);
              pure[s] = 0;
              if (!status)
                return status;
              double yf;
              auto ysview = species.data[s].trial;
              status = reconstruct_cartesian_convection_face(
                  kernels, ys, ysview, 0, axis, f, mass, yf);
              if (!status)
                return status;
              yf_h += static_cast<long double>(dh) * yf;
              yf_cp += (cps - cpdep) * yf;
              for (int j = first; j < last; ++j) {
                double yj = ysview.unchecked(offset(f, a, j - 2), 0);
                thermal[j] -= static_cast<long double>(dh) * yj;
                anchor[j] += (cps - cpdep) * yj;
              }
            }
            std::array<double, 4> q{};
            for (int j = 0; j < 4; ++j)
              q[j] = static_cast<double>(thermal[j]);
            const double tf_sensitivity =
                yf_cp - sampled_convection_direction(kernels, hs, q, anchor,
                                                     axis, f, mass);
            auto &result = faces[a][i];
            for (int j = 0; j < 4; ++j) {
              std::array<double, 4> basis{};
              basis[j] = 1;
              result.h[j] =
                  mass * (sampled_convection_direction(kernels, hs, q, basis,
                                                       axis, f, mass) -
                          sampled_convection_direction(kernels, hs, original,
                                                       basis, axis, f, mass));
            }
            result.h[1] += mass * tf_sensitivity *
                           interpolate_face(kernels, axis, normal,
                                            1.0 / cp.unchecked(left, 0), 0);
            result.h[2] += mass * tf_sensitivity *
                           interpolate_face(kernels, axis, normal, 0,
                                            1.0 / cp.unchecked(f, 0));
            double ordinary =
                sampled_convection_face(kernels, hs, original, axis, f, mass);
            double correction = static_cast<double>(
                static_cast<long double>(
                    sampled_convection_face(kernels, hs, q, axis, f, mass)) +
                yf_h - ordinary);
            result.p = coeff[a].unchecked(f) * correction;
            if (!std::isfinite(result.p))
              return {StatusCode::numerical_failure, 1583U};
            for (double weight : result.h)
              if (!std::isfinite(weight))
                return {StatusCode::numerical_failure, 1583U};
          }
    }
    return {};
  }
  double face(int a, Int3 f, ConstFieldView p, bool pressure) const {
    const auto &w = faces[a][index(a, f)];
    double v = 0;
    for (int j = 0; j < 4; ++j)
      if (w.h[j] != 0)
        v += w.h[j] * eliminated.unchecked(offset(f, a, j - 2), 0);
    if (pressure) {
      v = -v;
      if (w.p != 0)
        v += w.p * (p.unchecked(offset(f, a, -1), 0) - p.unchecked(f, 0));
    }
    return v;
  }
  void add(FieldView out, ConstFieldView p, bool pressure) const {
    cells([&](Int3 c) {
      double v = 0;
      for (int a = 0; a < 3; ++a)
        v += face(a, offset(c, a, 1), p, pressure) - face(a, c, p, pressure);
      out.unchecked(c, 0) += v;
    });
  }
  void add_pressure(FieldView out, ConstFieldView p) const {
    cells([&](Int3 c) {
      double v = 0;
      for (int a = 0; a < 3; ++a) {
        auto value = [&](Int3 f) {
          auto w = faces[a][index(a, f)].p;
          return w == 0 ? 0
                        : w * (p.unchecked(offset(f, a, -1), 0) -
                               p.unchecked(f, 0));
        };
        v += value(offset(c, a, 1)) - value(c);
      }
      out.unchecked(c, 0) += v;
    });
  }
};
} // namespace detail
} // namespace hundun::v04
