// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#pragma once
#include "hundun/v04_flow.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_schur_kin_detail.hpp"
#include "solver_schur_mix_detail.hpp"
#include <array>
#include <cstring>
#include <type_traits>
#include <vector>
namespace hundun::v04::detail {
// Coupled quasi-Newton action with internal-fluid-face EOS advection,
// mixture convection and temporal kinetic response. Physical boundary and
// transport-material derivatives retain the frozen approximation.
class PressureEnergyCoupledSchur : public LinearOperator {
public:
  PressureEnergyCoupledSchur() = default;
  PressureEnergyCoupledSchur(const PressureEnergyCoupledSchur &) = delete;
  PressureEnergyCoupledSchur &
  operator=(const PressureEnergyCoupledSchur &) = delete;
  const LinearOperator *cp{}, *ep{}, *eh{};
  const PressureEnergySchurOperator *pressure_pair{};
  ReductionEngine *reductions{};
  FieldId pressure_input_field{};
  SchurMixture *mixture{};
  SchurKinetic *kinetic{};
  SchurMixture mixture_storage;
  SchurKinetic kinetic_storage;
  PressureContinuityActivityView continuity_activity;
  bool rho_epoch{};
  bool enthalpy_halo_epoch{};
  bool sealed{};
  mutable LinearOperatorFailureProvenance failure_{};
  const CartesianKernelPlan *kernels{};
  ConstFieldView rp, rh, velocity;
  PressureEnergyCellActivity active;
  std::array<ConstFaceFieldView, 3> coeff, hface;
  LinearOperatorCertificate cert;
  double a0{};
  const PressureEnergyEnthalpyOperator *spatial{};
  mutable PressureEnergyEnthalpyPreparedEpoch epoch;
  ~PressureEnergyCoupledSchur() override { close(); }
  Status close() noexcept {
    Status s;
    if (enthalpy_halo_epoch)
      spatial->close_schur_prepared_halo(false, -1);
    enthalpy_halo_epoch = false;
    if (epoch.valid())
      s = spatial->close_repeated_apply(epoch);
    if (rho_epoch)
      halo.close_prepared_epoch(false);
    rho_epoch = false;
    sealed = false;
    spatial = nullptr;
    return s;
  }
  LinearOperatorFailureProvenance failure_provenance() const noexcept override {
    return failure_;
  }
  Status action(const LinearOperator &op, FieldView in, FieldView out) const {
    auto s = op.apply(in, out);
    if (!s)
      failure_ = op.failure_provenance();
    return s;
  }

  Status apply_eh() const {
    Status deferred;
    auto s = enthalpy_halo_epoch
                 ? spatial->apply_schur_prepared(fields[5], fields[6], deferred)
                 : eh->apply(fields[5], fields[6]);
    if (s)
      s = deferred;
    if (!s)
      failure_ = eh->failure_provenance();
    if (s && mixture)
      mixture->add(fields[6], {}, false);
    return s;
  }
  std::array<std::vector<double>, 3> left_weight, right_weight;
  std::size_t face_index(int a, Int3 f) const {
    auto n = cert.local_shape;
    (a == 0 ? n.x : a == 1 ? n.y : n.z)++;
    return std::size_t(f.x) +
           std::size_t(n.x) * (std::size_t(f.y) + std::size_t(n.y) * f.z);
  }
  PressureCorrectionBoundaryPlan boundary;
  mutable HaloEngine halo;
  mutable std::array<std::vector<double>, 9> memory;
  mutable std::array<FieldView, 9> fields;
  template <class F> void cells(F fn) const {
    const auto n = cert.local_shape;
    std::size_t i = 0;
    for (int z = 0; z < n.z; ++z)
      for (int y = 0; y < n.y; ++y)
        for (int x = 0; x < n.x; ++x, ++i)
          fn(Int3{x, y, z}, i);
  }
  bool fluid(std::size_t i) const {
    return !active.cells.size || active.cells.data[i];
  }
  Status allocate(const CartesianGeometryPlan &geom, MeshPatch patch,
                  const BoundaryPlan &bc, FieldId rho_field, FieldId hfield,
                  FieldId pfield, FieldId gfield, std::size_t species,
                  int donor_reach) {
    if (!valid_cells(patch.cells) || halo.ready())
      return {StatusCode::invalid_plan, 1584U};
    cert.local_shape = patch.cells;
    const auto n = cert.local_shape;
    const int g = 2;
    const std::size_t sy = n.x + 2 * g, sz = sy * (n.y + 2 * g),
                      count = sz * (n.z + 2 * g);
    for (std::size_t i = 0; i < 9; ++i) {
      memory[i].resize(count);
      auto &v = fields[i];
      v.base = memory[i].data() + g + g * sy + g * sz;
      v.interior = n;
      v.ghosts = {g, g, g};
      v.components = 1;
      v.stride_y = sy;
      v.stride_z = sz;
      v.component_stride = count;
      v.field = i == 5 ? hfield : rho_field;
      v.revision = 1;
      v.storage_identity = reinterpret_cast<StorageIdentity>(memory[i].data());
      v.revision_domain = reinterpret_cast<RevisionDomainIdentity>(this);
    }
    auto status =
        PressureCorrectionBoundaryPlan::compile(geom, patch, bc, boundary);
    for (int a = 0; a < 3; ++a) {
      auto e = n;
      (a == 0 ? e.x : a == 1 ? e.y : e.z)++;
      const auto nf = std::size_t(e.x) * e.y * e.z;
      left_weight[a].resize(nf);
      right_weight[a].resize(nf);
    }
    mixture_storage.reserve(n, species);
    kinetic_storage.reserve(n, donor_reach, pfield, gfield);
    return status;
  }
  Status bind_halo(MPI_Comm comm, MeshPatch patch, const BoundaryPlan &bc) {
    const HaloFieldSpec spec{fields[0].field, 1, 1};
    return halo.reserve(comm, patch, {&spec, 1}, bc.halo_topology());
  }
  Status prepare() {
    if (rho_epoch || epoch.valid() || !cp || !ep || !eh || !kernels ||
        !halo.ready())
      return {StatusCode::invalid_plan, 1584U};
    Status status;
    const auto n = cert.local_shape;
    sealed = false;
    const std::size_t count = std::size_t(n.x) * n.y * n.z;
    if (!same_cells(n, kernels->cells()) || !std::isfinite(a0) || a0 <= 0 ||
        !valid_cell_view(rp, n, 0U, 1U, 0U) ||
        !valid_cell_view(rh, n, 0U, 1U, 0U) ||
        !valid_cell_view(velocity, n, 0U, 3U, 1U) ||
        (active.cells.size &&
         (active.cells.size != count || !active.cells.data)))
      return {StatusCode::invalid_plan, 1584U};
    for (int a = 0; a < 3; ++a) {
      const auto flags = a == 0   ? continuity_activity.x_faces
                         : a == 1 ? continuity_activity.y_faces
                                  : continuity_activity.z_faces;
      auto e = n;
      (a == 0 ? e.x : a == 1 ? e.y : e.z)++;
      const auto nf = std::size_t(e.x) * e.y * e.z;
      if (!SchurMixture::valid_face(coeff[a], n, a) ||
          !SchurMixture::valid_face(hface[a], n, a) ||
          (flags.size && (flags.size != nf || !flags.data)))
        return {StatusCode::invalid_plan, 1584U};
    }
    bool eos_valid = true;
    cells([&](Int3 c, std::size_t i) {
      if (fluid(i))
        eos_valid = eos_valid && std::isfinite(rp.unchecked(c, 0)) &&
                    rp.unchecked(c, 0) > 0 &&
                    std::isfinite(rh.unchecked(c, 0)) && rh.unchecked(c, 0) < 0;
    });
    if (!eos_valid)
      return {StatusCode::numerical_failure, 1584U};
    for (auto &f : fields)
      ++f.revision;
    std::fill(kinetic_storage.pressure_memory.begin(),
              kinetic_storage.pressure_memory.end(), 0.0);
    if (status)
      for (int a = 0; a < 3; ++a) {
        auto e = n;
        (a == 0 ? e.x : a == 1 ? e.y : e.z)++;
        const auto axis = static_cast<CartesianAxis>(a);
        const std::size_t nf = std::size_t(e.x) * e.y * e.z;
        if (left_weight[a].size() != nf || right_weight[a].size() != nf)
          return {StatusCode::invalid_plan, 1584U};
        std::fill(left_weight[a].begin(), left_weight[a].end(), 0.0);
        std::fill(right_weight[a].begin(), right_weight[a].end(), 0.0);
        for (int z = 0; z < e.z; ++z)
          for (int y = 0; y < e.y; ++y)
            for (int x = 0; x < e.x; ++x) {
              const Int3 f{x, y, z};
              PressureCorrectionFaceRule rule;
              status = boundary.face_rule(axis, f, rule);
              if (!status)
                return status;
              const auto flags = a == 0   ? continuity_activity.x_faces
                                 : a == 1 ? continuity_activity.y_faces
                                          : continuity_activity.z_faces;
              if ((flags.size && !flags.data[face_index(a, f)]) ||
                  rule.is_nonperiodic_boundary() ||
                  coeff[a].unchecked(f) == 0.0)
                continue;
              Int3 l = f;
              (a == 0 ? l.x : a == 1 ? l.y : l.z)--;
              const auto i = face_index(a, f);
              const int normal = a == 0 ? f.x : a == 1 ? f.y : f.z;
              const double area = detail::face_area(*kernels, axis, f);
              left_weight[a][i] = area * detail::interpolate_face(
                                             *kernels, axis, normal,
                                             velocity.unchecked(l, a), 0.0);
              right_weight[a][i] =
                  area * detail::interpolate_face(*kernels, axis, normal, 0.0,
                                                  velocity.unchecked(f, a));
              if (!std::isfinite(left_weight[a][i]) ||
                  !std::isfinite(right_weight[a][i]) ||
                  !std::isfinite(hface[a].unchecked(f)))
                return {StatusCode::numerical_failure, 1584U};
            }
      }
    if (status) {
      status = halo.enter_prepared_epoch();
      rho_epoch = bool(status);
    }
    if (status)
      status = prepare_enthalpy();
    return status;
  }
  Status prepare_enthalpy() {
    if (epoch.valid() || enthalpy_halo_epoch)
      return {StatusCode::invalid_plan, 1584U};
    spatial = dynamic_cast<const PressureEnergyEnthalpyOperator *>(eh);
    if (spatial && spatial->enthalpy_certificate().compiled_factored_apply) {
      auto status = spatial->prepare_repeated_apply(epoch);
      if (status) {
        status = spatial->enter_schur_prepared_halo();
        enthalpy_halo_epoch = bool(status);
      }
      return status;
    }
    return {};
  }
  std::uint64_t owned_payload_bytes() const noexcept {
    std::uint64_t total = 0;
    const auto add = [&](const auto &v) {
      using T = typename std::decay_t<decltype(v)>::value_type;
      total += v.capacity() * sizeof(T);
    };
    for (const auto &v : memory)
      add(v);
    for (const auto &v : left_weight)
      add(v);
    for (const auto &v : right_weight)
      add(v);
    for (const auto &v : mixture_storage.faces)
      add(v);
    add(mixture_storage.pure);
    add(kinetic_storage.memory);
    add(kinetic_storage.pressure_memory);
    return total;
  }
  Status seal() noexcept {
    if (!rho_epoch || !mixture || !kinetic)
      return {StatusCode::invalid_plan, 1584U};
    const auto kinetic_status = kinetic->validate();
    if (!kinetic_status)
      return kinetic_status;
    auto &hash = cert.collective_fingerprint;
    const auto mix = [&](std::uint64_t value) {
      hash ^= value;
      hash *= UINT64_C(1099511628211);
    };
    mix(UINT64_C(0x72686f7363687572));
    mix(mixture->target_fingerprint);
    mix(continuity_activity.collective_fingerprint);
    mix(pressure_pair &&
        pressure_pair->jacobian_certificate().shared_pressure.valid());
    std::uint64_t bits{};
    std::memcpy(&bits, &a0, sizeof(bits));
    mix(bits);
    for (auto field : {rp, rh, velocity, kinetic->density, kinetic->velocity,
                       kinetic->r_au}) {
      mix(field.field);
      mix(field.revision);
    }
    if (!hash)
      hash = 1;
    sealed = true;
    return {};
  }
  LinearOperatorCertificate certificate() const noexcept override {
    return sealed ? cert : LinearOperatorCertificate{};
  }
  // fields: rho, D(rho), DE(rho), Cp, Ep, h, Eh, rhs, reserved.
  Status advect(bool energy = true) const {
    HaloTicket ticket;
    Status deferred;
    auto status = halo.begin_prepared(181, {&fields[0], 1}, deferred, ticket);
    if (status)
      status = halo.finish_prepared(ticket, {&fields[0], 1}, deferred);
    if (!status)
      return status;
    if (!deferred)
      return deferred;
    const auto r = as_const(fields[0]);
    const auto face = [&](int axis, Int3 f) -> double {
      const auto i = face_index(axis, f);
      if (left_weight[axis][i] == 0 && right_weight[axis][i] == 0)
        return 0;
      Int3 l = f;
      (axis == 0 ? l.x : axis == 1 ? l.y : l.z)--;
      return left_weight[axis][i] * r.unchecked(l, 0) +
             right_weight[axis][i] * r.unchecked(f, 0);
    };
    cells([&](Int3 c, std::size_t i) {
      double dc = 0, de = 0;
      if (fluid(i))
        for (int a = 0; a < 3; ++a) {
          Int3 f = c;
          (a == 0 ? f.x : a == 1 ? f.y : f.z)++;
          const auto heat = [&](Int3 f) {
            double h = hface[a].unchecked(f);
            if (mixture && coeff[a].unchecked(f) != 0)
              h +=
                  mixture->faces[a][face_index(a, f)].p / coeff[a].unchecked(f);
            return h;
          };
          const double lo = face(a, c), hi = face(a, f);
          dc += hi - lo;
          if (energy)
            de += (hi == 0.0 ? 0.0 : hi * heat(f)) -
                  (lo == 0.0 ? 0.0 : lo * heat(c));
        }
      fields[1].unchecked(c, 0) = dc;
      fields[2].unchecked(c, 0) = de;
    });
    return {};
  }
  // Fixed polynomial inverse, hence a linear operator at every Krylov call.
  Status inverse_ch(ConstFieldView rhs) const {
    cells([&](Int3 c, std::size_t i) {
      fields[7].unchecked(c, 0) = fluid(i) ? rhs.unchecked(c, 0) : 0;
      fields[0].unchecked(c, 0) =
          fluid(i)
              ? rhs.unchecked(c, 0) / (a0 * detail::cell_volume(*kernels, c))
              : 0;
    });
    for (int it = 0; it < 2; ++it) {
      // Intermediate inverse iterations consume only D(rho). The energy
      // response is evaluated once, on the final density direction below.
      auto s = advect(false);
      if (!s)
        return s;
      cells([&](Int3 c, std::size_t i) {
        fields[0].unchecked(c, 0) =
            fluid(i) ? (fields[7].unchecked(c, 0) - fields[1].unchecked(c, 0)) /
                           (a0 * detail::cell_volume(*kernels, c))
                     : 0;
      });
    }
    cells([&](Int3 c, std::size_t i) {
      fields[5].unchecked(c, 0) =
          fluid(i) ? fields[0].unchecked(c, 0) / rh.unchecked(c, 0) : 0;
    });
    return advect();
  }
  Status pressure_blocks(FieldView p) const {
    failure_ = {};
    if (pressure_pair && reductions &&
        pressure_pair->jacobian_certificate().shared_pressure.valid()) {
      const auto status = pressure_pair->apply_pressure_pair(
          p, fields[3], fields[4], *reductions);
      if (!status)
        failure_ = pressure_pair->failure_provenance();
      return status;
    }
    // Generic components can reject metadata before entering their Halo.
    // Agree the common input/workspace contract before invoking either one.
    Status status =
        cp && ep && p.field == pressure_input_field &&
                valid_cell_view(as_const(p), cert.local_shape, 0U, 1U, 1U) &&
                valid_cell_view(fields[3], cert.local_shape, 0U, 1U) &&
                valid_cell_view(fields[4], cert.local_shape, 0U, 1U) &&
                !field_views_overlap(as_const(p), as_const(fields[3])) &&
                !field_views_overlap(as_const(p), as_const(fields[4])) &&
                !field_views_overlap(as_const(fields[3]), as_const(fields[4]))
            ? Status{}
            : Status{StatusCode::invalid_plan, 1584U};
    const auto agree = [&](Status local) {
      if (!reductions)
        return local;
      const auto prior = failure_;
      const auto global = reductions->consensus(local);
      if (!global) {
        const bool collective_origin =
            prior.status.code == global.code &&
            prior.status.detail == global.detail &&
            prior.status_scope == LinearOperatorStatusScope::collective;
        failure_ = {global, LinearOperatorStatusScope::collective,
                    collective_origin ? prior.lowest_failing_rank
                                      : reductions->lowest_failing_rank()};
      }
      return global;
    };
    status = agree(status);
    if (status)
      status = agree(action(*cp, p, fields[3]));
    if (status)
      status = agree(action(*ep, p, fields[4]));
    return status;
  }
  void remove_pressure_time(FieldView p) const {
    cells([&](Int3 c, std::size_t i) {
      if (fluid(i))
        fields[3].unchecked(c, 0) -= a0 * detail::cell_volume(*kernels, c) *
                                     rp.unchecked(c, 0) * p.unchecked(c, 0);
    });
  }
  Status full_cp(FieldView p) const {
    const auto status = action(*cp, p, fields[3]);
    if (status)
      remove_pressure_time(p);
    return status;
  }
  Status apply(FieldView p, FieldView out) const noexcept override {
    failure_ = {};
    if (!sealed)
      return {StatusCode::invalid_plan, 1584U};
    auto s = pressure_blocks(p);
    if (s)
      remove_pressure_time(p);
    if (s && kinetic) {
      s = kinetic->add(p, fields[4]);
      if (!s)
        failure_ = kinetic->failure;
    }
    if (!s)
      return s;
    s = inverse_ch(as_const(fields[3]));
    if (!s)
      return s;
    cells([&](Int3 c, std::size_t i) {
      if (fluid(i))
        fields[5].unchecked(c, 0) +=
            rp.unchecked(c, 0) * p.unchecked(c, 0) / rh.unchecked(c, 0);
    });
    s = apply_eh();
    if (!s)
      return s;
    cells([&](Int3 c, std::size_t i) {
      out.unchecked(c, 0) = fluid(i) ? fields[4].unchecked(c, 0) -
                                           fields[6].unchecked(c, 0) -
                                           fields[2].unchecked(c, 0)
                                     : p.unchecked(c, 0);
    });
    if (mixture)
      mixture->add_pressure(out, as_const(p));
    return {};
  }
  // Evaluate the unreduced rows on an arbitrary joint direction. This shares
  // the exact bound face actions used by elimination and recovery.
  Status joint_linear(FieldView p, ConstFieldView h, FieldView c_out,
                      FieldView e_out) const {
    failure_ = {};
    if (!sealed)
      return {StatusCode::invalid_plan, 1584U};
    auto s = pressure_blocks(p);
    if (s) {
      s = kinetic->add(p, fields[4]);
      if (!s)
        failure_ = kinetic->failure;
    }
    if (!s)
      return s;
    cells([&](Int3 c, std::size_t i) {
      fields[5].unchecked(c, 0) = fluid(i) ? h.unchecked(c, 0) : 0.0;
      fields[0].unchecked(c, 0) =
          fluid(i) ? rp.unchecked(c, 0) * p.unchecked(c, 0) +
                         rh.unchecked(c, 0) * h.unchecked(c, 0)
                   : 0.0;
    });
    s = advect();
    if (s)
      s = apply_eh();
    if (!s)
      return s;
    cells([&](Int3 c, std::size_t i) {
      c_out.unchecked(c, 0) = fluid(i) ? fields[3].unchecked(c, 0) +
                                             a0 * cell_volume(*kernels, c) *
                                                 rh.unchecked(c, 0) *
                                                 h.unchecked(c, 0) +
                                             fields[1].unchecked(c, 0)
                                       : p.unchecked(c, 0);
      e_out.unchecked(c, 0) = fluid(i) ? fields[4].unchecked(c, 0) +
                                             fields[6].unchecked(c, 0) +
                                             fields[2].unchecked(c, 0)
                                       : 0.0;
    });
    mixture->add_pressure(e_out, as_const(p));
    return {};
  }
  Status rhs(ConstFieldView rc, ConstFieldView re, FieldView out) const {
    failure_ = {};
    if (!sealed)
      return {StatusCode::invalid_plan, 1584U};
    auto s = inverse_ch(rc);
    if (s)
      s = apply_eh();
    if (!s)
      return s;
    cells([&](Int3 c, std::size_t i) {
      out.unchecked(c, 0) = fluid(i) ? -re.unchecked(c, 0) +
                                           fields[6].unchecked(c, 0) +
                                           fields[2].unchecked(c, 0)
                                     : 0;
    });
    return {};
  }
  Status recover(ConstFieldView rc, FieldView p, FieldView h) const {
    failure_ = {};
    if (!sealed)
      return {StatusCode::invalid_plan, 1584U};
    auto s = full_cp(p);
    if (!s)
      return s;
    cells([&](Int3 c, std::size_t i) {
      fields[3].unchecked(c, 0) += rc.unchecked(c, 0);
    });
    s = inverse_ch(as_const(fields[3]));
    if (!s)
      return s;
    cells([&](Int3 c, std::size_t i) {
      h.unchecked(c, 0) =
          fluid(i)
              ? -fields[5].unchecked(c, 0) -
                    rp.unchecked(c, 0) * p.unchecked(c, 0) / rh.unchecked(c, 0)
              : 0;
    });
    return {};
  }
};
} // namespace hundun::v04::detail
