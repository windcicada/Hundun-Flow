// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/ibm_force_fixture.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

class SevenPointOperator final : public LinearOperator {
 public:
  explicit SevenPointOperator(Int3 cells) : cells_(cells) {
    certificate_.identity = {1U, 2U, 3U, 4U, 5U};
    certificate_.collective_fingerprint = 6U;
    certificate_.local_shape = cells;
    certificate_.operator_class = LinearOperatorClass::spd;
  }
  LinearOperatorCertificate certificate() const noexcept override {
    return certificate_;
  }
  Status apply(FieldView x, FieldView y) const noexcept override {
    for (std::int32_t z = 0; z < cells_.z; ++z) {
      for (std::int32_t yy = 0; yy < cells_.y; ++yy) {
        for (std::int32_t xx = 0; xx < cells_.x; ++xx) {
          const Int3 cell{xx, yy, z};
          y.unchecked(cell, 0U) =
              7.0 * x.unchecked(cell, 0U) -
              x.unchecked({xx - 1, yy, z}, 0U) -
              x.unchecked({xx + 1, yy, z}, 0U) -
              x.unchecked({xx, yy - 1, z}, 0U) -
              x.unchecked({xx, yy + 1, z}, 0U) -
              x.unchecked({xx, yy, z - 1}, 0U) -
              x.unchecked({xx, yy, z + 1}, 0U);
        }
      }
    }
    return {};
  }

 private:
  Int3 cells_{};
  LinearOperatorCertificate certificate_{};
};

class ProvenanceOperator final : public LinearOperator {
 public:
  explicit ProvenanceOperator(Int3 cells) : cells_(cells) {
    certificate_.identity = {11U, 12U, 13U, 14U, 15U};
    certificate_.collective_fingerprint = 16U;
    certificate_.local_shape = cells;
    certificate_.operator_class = LinearOperatorClass::spd;
  }

  LinearOperatorCertificate certificate() const noexcept override {
    return certificate_;
  }

  Status apply(FieldView x, FieldView y) const noexcept override {
    failure_ = {};
    if (fail_) {
      const Status status{StatusCode::mpi_failure, 95001U};
      failure_ = {status, LinearOperatorStatusScope::collective, 0};
      return status;
    }
    for (std::int32_t z = 0; z < cells_.z; ++z) {
      for (std::int32_t yy = 0; yy < cells_.y; ++yy) {
        for (std::int32_t xx = 0; xx < cells_.x; ++xx) {
          y.unchecked({xx, yy, z}, 0U) =
              x.unchecked({xx, yy, z}, 0U);
        }
      }
    }
    return {};
  }

  LinearOperatorFailureProvenance failure_provenance() const
      noexcept override {
    return failure_;
  }

  void fail(bool value) noexcept { fail_ = value; }

 private:
  Int3 cells_{};
  LinearOperatorCertificate certificate_{};
  bool fail_{true};
  mutable LinearOperatorFailureProvenance failure_{};
};

struct OwnedFace {
  std::vector<double> storage;
  FaceFieldView view{};
};

OwnedFace face(CartesianAxis axis, Int3 cells, StorageIdentity identity) {
  OwnedFace result;
  Int3 extents = cells;
  if (axis == CartesianAxis::x) ++extents.x;
  if (axis == CartesianAxis::y) ++extents.y;
  if (axis == CartesianAxis::z) ++extents.z;
  const std::size_t stride_y = static_cast<std::size_t>(extents.x);
  const std::size_t stride_z = stride_y * extents.y;
  result.storage.assign(stride_z * extents.z, 1.0);
  result.view = {result.storage.data(), extents, stride_y, stride_z, axis,
                 identity, 94001U};
  return result;
}

std::size_t flat(Int3 cells, Int3 cell) {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(cells.y) * cell.z);
}

bool test_exact_compact_correction() {
  IbmForceFixture fixture;
  bool passed = expect(fixture.initialize(),
                       "IBM pressure operator fixture compiles");
  if (!passed) return false;
  const Int3 cells = fixture.patch.cells;
  const std::uint8_t ghosts = fixture.boundary.maximum_halo_reach();
  ForceOwnedField x = make_force_field(20U, cells, 1U, ghosts, 501U, 601U);
  ForceOwnedField product = make_force_field(21U, cells, 1U, 0U, 502U, 602U);
  ForceOwnedField regular_value =
      make_force_field(22U, cells, 1U, 0U, 503U, 603U);
  for (std::int32_t z = -ghosts; z < cells.z + ghosts; ++z) {
    for (std::int32_t y = -ghosts; y < cells.y + ghosts; ++y) {
      for (std::int32_t xx = -ghosts; xx < cells.x + ghosts; ++xx) {
        const Int3 global{fixture.patch.begin.x + xx,
                          fixture.patch.begin.y + y,
                          fixture.patch.begin.z + z};
        const double px = fixture.extrapolated_centre(fixture.geometry.x(),
                                                      global.x);
        const double py = fixture.extrapolated_centre(fixture.geometry.y(),
                                                      global.y);
        const double pz = fixture.extrapolated_centre(fixture.geometry.z(),
                                                      global.z);
        x.view.unchecked({xx, y, z}, 0U) =
            1.0 + 0.4 * px - 0.3 * py + 0.2 * pz +
            0.1 * (px * px + py * py + pz * pz);
      }
    }
  }
  OwnedFace fx = face(CartesianAxis::x, cells, 701U);
  OwnedFace fy = face(CartesianAxis::y, cells, 702U);
  OwnedFace fz = face(CartesianAxis::z, cells, 703U);
  SevenPointOperator regular(cells);
  passed &= expect(static_cast<bool>(regular.apply(x.view,
                                                   regular_value.view)),
                   "regular Cartesian oracle applies");
  IbmPressureOperator ibm;
  passed &= expect(IbmPressureOperator::bind(
                       regular, fixture.topology, fixture.boundary,
                       as_const(fx.view), as_const(fy.view), as_const(fz.view),
                       fixture.geometry.topology_revision(), ibm) &&
                       ibm.certificate().operator_class ==
                           LinearOperatorClass::spd &&
                       ibm.apply(x.view, product.view),
                   "exact conservative IBM pressure operator applies");
  if (!passed) return false;

  std::vector<double> oracle = regular_value.storage;
  const Span<const std::uint8_t> region = fixture.topology.region();
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t xx = 0; xx < cells.x; ++xx) {
        const Int3 cell{xx, y, z};
        if (region.data[flat(cells, cell)] ==
            static_cast<std::uint8_t>(RegionFlag::solid)) {
          oracle[flat(cells, cell)] = x.view.unchecked(cell, 0U);
        }
      }
    }
  }
  const auto bindings = fixture.boundary.links();
  const auto links = fixture.topology.links();
  for (std::size_t index = 0U; index < bindings.size; ++index) {
    const BoundaryStencilLink& binding = bindings.data[index];
    const ImmersedLink& link = links.data[binding.topology_link];
    oracle[flat(cells, link.fluid_local_index)] +=
        x.view.unchecked(link.solid_local_index, 0U) -
        x.view.unchecked(link.fluid_local_index, 0U);
  }
  double maximum_error = 0.0;
  double maximum_compact_change = 0.0;
  for (std::size_t index = 0U; index < oracle.size(); ++index) {
    maximum_error =
        std::max(maximum_error, std::abs(product.storage[index] - oracle[index]));
    maximum_compact_change = std::max(
        maximum_compact_change,
        std::abs(product.storage[index] - regular_value.storage[index]));
  }
  passed &= expect(maximum_error < 2.0e-12 && maximum_compact_change > 1.0e-4,
                   "product removes impermeable interface flux exactly");

  ForceOwnedField rhs = make_force_field(23U, cells, 1U, 0U, 504U, 604U);
  std::fill(rhs.storage.begin(), rhs.storage.end(), 3.0);
  passed &= expect(static_cast<bool>(ibm.mask_solid_rhs(rhs.view)),
                   "solid RHS rows are masked to the identity equation");
  for (std::size_t index = 0U; index < rhs.storage.size(); ++index) {
    const double expected =
        region.data[index] == static_cast<std::uint8_t>(RegionFlag::solid)
            ? 0.0
            : 3.0;
    passed &= expect(rhs.storage[index] == expected,
                     "fluid RHS is retained and solid RHS is zero");
  }

  IbmPressureOperator rejected;
  passed &= expect(IbmPressureOperator::bind(
                       regular, fixture.topology, fixture.boundary,
                       as_const(fx.view), as_const(fy.view), as_const(fz.view),
                       fixture.geometry.topology_revision() + 1U, rejected)
                           .code == StatusCode::invalid_plan &&
                       rejected.fingerprint() == 0U,
                   "stale geometry cannot bind compact pressure rows");
  return passed;
}

bool test_failure_provenance_delegation() {
  IbmForceFixture fixture;
  bool passed = expect(fixture.initialize(),
                       "IBM failure-provenance fixture compiles");
  if (!passed) return false;
  const Int3 cells = fixture.patch.cells;
  const std::uint8_t ghosts = fixture.boundary.maximum_halo_reach();
  ForceOwnedField input =
      make_force_field(30U, cells, 1U, ghosts, 95101U, 95201U);
  ForceOwnedField output =
      make_force_field(31U, cells, 1U, 0U, 95102U, 95202U);
  std::fill(input.storage.begin(), input.storage.end(), 1.0);
  OwnedFace fx = face(CartesianAxis::x, cells, 95301U);
  OwnedFace fy = face(CartesianAxis::y, cells, 95302U);
  OwnedFace fz = face(CartesianAxis::z, cells, 95303U);
  ProvenanceOperator regular(cells);
  IbmPressureOperator ibm;
  passed &= expect(
      static_cast<bool>(IbmPressureOperator::bind(
          regular, fixture.topology, fixture.boundary, as_const(fx.view),
          as_const(fy.view), as_const(fz.view),
          fixture.geometry.topology_revision(), ibm)),
      "IBM failure-provenance adapter binds");
  if (!passed) return false;

  const Status delegated = ibm.apply(input.view, output.view);
  const LinearOperatorFailureProvenance delegated_provenance =
      ibm.failure_provenance();
  passed &= expect(
      delegated.code == StatusCode::mpi_failure &&
          delegated.detail == 95001U &&
          delegated_provenance.status.code == delegated.code &&
          delegated_provenance.status.detail == delegated.detail &&
          delegated_provenance.status_scope ==
              LinearOperatorStatusScope::collective &&
          delegated_provenance.lowest_failing_rank == 0,
      "IBM snapshots exact delegated collective provenance");

  const Status local_invalid = ibm.apply(input.view, input.view);
  const LinearOperatorFailureProvenance after_local_invalid =
      ibm.failure_provenance();
  passed &= expect(
      local_invalid.code == StatusCode::invalid_plan &&
          after_local_invalid.status &&
          after_local_invalid.status_scope ==
              LinearOperatorStatusScope::rank_local &&
          after_local_invalid.lowest_failing_rank == -1,
      "IBM local validation clears stale delegated provenance");

  regular.fail(false);
  const Status retry = ibm.apply(input.view, output.view);
  const LinearOperatorFailureProvenance after_retry =
      ibm.failure_provenance();
  passed &= expect(
      static_cast<bool>(retry) && after_retry.status &&
          after_retry.status_scope ==
              LinearOperatorStatusScope::rank_local &&
          after_retry.lowest_failing_rank == -1,
      "IBM successful retry keeps provenance empty");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  const bool passed = test_exact_compact_correction() &&
                      test_failure_provenance_delegation();
  MPI_Finalize();
  return passed ? 0 : 1;
}
