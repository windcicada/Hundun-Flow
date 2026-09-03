// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/piso_fixture.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string_view>

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

double face_sum(const PeriodicPisoFixture& fixture, Int3 cell) {
  return fixture.x_pressure_coefficient.view.unchecked(cell) +
         fixture.x_pressure_coefficient.view.unchecked(
             {cell.x + 1, cell.y, cell.z}) +
         fixture.y_pressure_coefficient.view.unchecked(cell) +
         fixture.y_pressure_coefficient.view.unchecked(
             {cell.x, cell.y + 1, cell.z}) +
         fixture.z_pressure_coefficient.view.unchecked(cell) +
         fixture.z_pressure_coefficient.view.unchecked(
             {cell.x, cell.y, cell.z + 1});
}

bool test_production_checkerboard_mode() {
  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(8, MPI_COMM_WORLD),
                       "production periodic PISO fixture compiles");
  if (!passed) {
    return false;
  }
  const BdfCoefficients bdf{10.0, -10.0, 0.0, 1U};
  const Int3 cells = fixture.patch.cells;
  for (std::int32_t z = -1; z <= cells.z; ++z) {
    for (std::int32_t y = -1; y <= cells.y; ++y) {
      for (std::int32_t x = -1; x <= cells.x; ++x) {
        fixture.pressure.view.unchecked({x, y, z}, 0U) =
            ((fixture.patch.begin.x + x + fixture.patch.begin.y + y +
              fixture.patch.begin.z + z) &
             1) == 0
                ? 1.0
                : -1.0;
      }
    }
  }
  ++fixture.pressure.view.revision;
  PisoIntermediateInput intermediate_input =
      fixture.intermediate_input(bdf, 5001U);
  PisoIntermediateCertificate intermediate;
  passed &= expect(static_cast<bool>(fixture.coupler.refresh(
                       intermediate_input, intermediate)),
                   "production rAU/HbyA/pressure-face lifecycle refreshes");
  double predictor_response_error = 0.0;
  double predictor_response_norm = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x <= cells.x; ++x) {
        const Int3 face{x, y, z};
        const Int3 left{x - 1, y, z};
        const double jump =
            fixture.pressure.view.unchecked(face, 0U) -
            fixture.pressure.view.unchecked(left, 0U);
        const double expected =
            -fixture.x_pressure_coefficient.view.unchecked(face) * jump;
        const double actual = fixture.phi_h_by_a.x.unchecked(face);
        predictor_response_error =
            std::max(predictor_response_error, std::abs(actual - expected));
        predictor_response_norm =
            std::max(predictor_response_norm, std::abs(expected));
      }
    }
  }
  passed &= expect(predictor_response_norm > 0.0 &&
                       predictor_response_error <=
                           1.0e-12 * predictor_response_norm,
                   "C1 predictor face flux consumes the current pressure "
                   "checkerboard through the spatial Rhie-Chow response");

  OwnedField accepted = make_field(50U, fixture.patch.cells, 1U, 0U, 5002U,
                                   6002U);
  OwnedField drho_dp = make_field(52U, fixture.patch.cells, 1U, 0U, 5003U,
                                  6003U);
  OwnedField diagonal = make_field(53U, fixture.patch.cells, 1U, 0U, 5004U,
                                   6004U);
  OwnedField rhs = make_field(54U, fixture.patch.cells, 1U, 0U, 5005U,
                              6005U);
  fill(accepted, 1.0);
  fill(drho_dp, 0.01);
  PressureCorrectionInput pressure_input;
  pressure_input.intermediate = intermediate;
  pressure_input.pressure_reference = intermediate_input.pressure_reference;
  pressure_input.density_trial = as_const(fixture.density.view);
  pressure_input.density_accepted = as_const(accepted.view);
  pressure_input.drho_dp_h_y = as_const(drho_dp.view);
  pressure_input.bdf = bdf;
  pressure_input.time = intermediate_input.momentum.time;
  pressure_input.geometry = intermediate_input.momentum.geometry;
  pressure_input.numeric_boundary = intermediate_input.numeric_boundary;
  const PressureCorrectionSystemView system{diagonal.view, rhs.view};
  PressureCorrectionCertificate pressure;
  passed &= expect(static_cast<bool>(fixture.coupler.assemble_pressure_system(
                       pressure_input, system, pressure)),
                   "production pressure system assembles");

  constexpr FieldId correction_field = 90U;
  OwnedField mode = make_field(correction_field, fixture.patch.cells, 1U, 1U,
                               5006U, 6006U);
  OwnedField applied = make_field(91U, fixture.patch.cells, 1U, 0U, 5007U,
                                  6007U);
  const std::array<HaloFieldSpec, 1U> halo_fields{{
      {correction_field, 1U, 1U}}};
  HaloEngine operator_halo;
  passed &= expect(static_cast<bool>(operator_halo.reserve(
                       MPI_COMM_WORLD, fixture.patch,
                       {halo_fields.data(), halo_fields.size()},
                       fixture.boundary.halo_topology())),
                   "pressure-operator halo reserves once");
  PressureLinearOperator pressure_operator;
  passed &= expect(static_cast<bool>(fixture.coupler.bind_pressure_operator(
                       {MPI_COMM_WORLD, &operator_halo, 5008U,
                        correction_field},
                       system, pressure_operator)) &&
                       static_cast<bool>(pressure_operator.refresh(
                           {pressure,
                            {5009U, 5010U, 5011U, 5012U, 5013U},
                            5014U})),
                   "production pressure operator binds and refreshes");
  if (!passed) {
    return false;
  }

  fill(mode, 1.0);
  passed &= expect(static_cast<bool>(pressure_operator.apply(mode.view,
                                                             applied.view)),
                   "constant pressure mode applies");
  double constant_numerator = 0.0;
  double constant_denominator = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        constant_numerator += applied.view.unchecked({x, y, z}, 0U);
        constant_denominator += 1.0;
      }
    }
  }
  const double constant_rayleigh =
      constant_numerator / constant_denominator;

  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        mode.view.unchecked({x, y, z}, 0U) =
            ((fixture.patch.begin.x + x + fixture.patch.begin.y + y +
              fixture.patch.begin.z + z) &
             1) == 0
                ? 1.0
                : -1.0;
      }
    }
  }
  ++mode.view.revision;
  passed &= expect(static_cast<bool>(pressure_operator.apply(mode.view,
                                                             applied.view)),
                   "three-dimensional checkerboard pressure mode applies");
  double checkerboard_numerator = 0.0;
  double checkerboard_denominator = 0.0;
  double oracle_numerator = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double volume =
            fixture.geometry.x().widths().data[static_cast<std::size_t>(x)] *
            fixture.geometry.y().widths().data[static_cast<std::size_t>(y)] *
            fixture.geometry.z().widths().data[static_cast<std::size_t>(z)];
        const double storage = bdf.a0 * volume * 0.01;
        const double value = mode.view.unchecked(cell, 0U);
        checkerboard_numerator +=
            value * applied.view.unchecked(cell, 0U);
        checkerboard_denominator += value * value;
        oracle_numerator += (storage + 2.0 * face_sum(fixture, cell)) *
                            value * value;
      }
    }
  }
  const double checkerboard_rayleigh =
      checkerboard_numerator / checkerboard_denominator;
  const double oracle_rayleigh = oracle_numerator / checkerboard_denominator;
  const double tolerance = 1.0e-12 * std::max(1.0, oracle_rayleigh);
  const double reference_volume =
      fixture.geometry.x().widths().data[0U] *
      fixture.geometry.y().widths().data[0U] *
      fixture.geometry.z().widths().data[0U];
  const double reference_storage = bdf.a0 * reference_volume * 0.01;
  passed &= expect(std::abs(constant_rayleigh - reference_storage) <
                       tolerance,
                   "constant mode contains only compressibility storage");
  passed &= expect(std::abs(checkerboard_rayleigh - oracle_rayleigh) <
                       tolerance,
                   "checkerboard eigenvalue matches independent face stencil");
  const double checkerboard_face_contribution =
      oracle_rayleigh - constant_rayleigh;
  passed &= expect(checkerboard_face_contribution > 10.0 * constant_rayleigh,
                   "volume-weighted face coupling strongly controls checkerboard mode");
  if (!passed) {
    std::cerr << "constant=" << constant_rayleigh
              << " checkerboard=" << checkerboard_rayleigh
              << " oracle=" << oracle_rayleigh
              << " predictor-response-error=" << predictor_response_error
              << " predictor-response-norm=" << predictor_response_norm
              << '\n';
  }
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  const bool passed = test_production_checkerboard_mode();
  MPI_Finalize();
  return passed ? 0 : 1;
}
