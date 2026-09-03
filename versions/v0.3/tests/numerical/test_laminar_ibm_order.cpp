// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/flow_immersed_test_access.hpp"
#include "tests/support/fvm_immersed_operator_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/ib_surface_query.hpp"
#include "hundun/lin_bicgstab.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/stage3_dual3.hpp"
#include "tests/support/stage3_mms.hpp"
#include "tests/support/test_main.hpp"
#include "src/ib_deterministic_qr_detail.hpp"

#include <mpi.h>

#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <thread>
#include <variant>
#include <vector>

namespace {

using namespace hundun;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr std::size_t kParallelWorkerBudget = 96U;

test::stage3::MmsCellAverage
rescale_manufactured_average(const test::stage3::MmsCellAverage &reference,
                             double reference_time_s, double time_s) {
  const double reference_amplitude =
      1.0 + 0.1 * std::sin(2.0 * kPi * reference_time_s);
  const double amplitude = 1.0 + 0.1 * std::sin(2.0 * kPi * time_s);
  const double reference_pressure_scale = std::cos(kPi * reference_time_s);
  HUNDUN_CHECK(std::abs(reference_amplitude) > 0.0);
  HUNDUN_CHECK(std::abs(reference_pressure_scale) > 0.0);
  const double velocity_scale = amplitude / reference_amplitude;
  const double pressure_scale =
      std::cos(kPi * time_s) / reference_pressure_scale;
  return {{velocity_scale * reference.velocity_m_per_s.x,
           velocity_scale * reference.velocity_m_per_s.y,
           velocity_scale * reference.velocity_m_per_s.z},
          pressure_scale * reference.mechanical_pressure_pa,
          {},
          {pressure_scale *
               reference.mechanical_pressure_gradient_pa_per_m.x,
           pressure_scale *
               reference.mechanical_pressure_gradient_pa_per_m.y,
           pressure_scale *
               reference.mechanical_pressure_gradient_pa_per_m.z}};
}

void check_dual3() {
  using test::stage3::Dual3;
  const auto x = Dual3::variable(0.3, 0U);
  const auto y = Dual3::variable(-0.2, 1U);
  const auto value = test::stage3::sin(x * y) + x * x * x;
  const double xy = 0.3 * -0.2;
  HUNDUN_CHECK_NEAR(value.value, std::sin(xy) + 0.3 * 0.3 * 0.3, 1.0e-15);
  HUNDUN_CHECK_NEAR(value.gradient[0], std::cos(xy) * -0.2 + 3.0 * 0.3 * 0.3,
                    1.0e-15);
  HUNDUN_CHECK_NEAR(value.gradient[1], std::cos(xy) * 0.3, 1.0e-15);
  HUNDUN_CHECK_NEAR(value.hessian[0][0], -std::sin(xy) * 0.04 + 6.0 * 0.3,
                    2.0e-15);
  HUNDUN_CHECK_NEAR(value.hessian[0][1], std::cos(xy) - std::sin(xy) * xy,
                    2.0e-15);
  HUNDUN_CHECK_NEAR(value.hessian[1][0], value.hessian[0][1], 0.0);
  HUNDUN_CHECK_NEAR(value.hessian[1][1], -std::sin(xy) * 0.09, 2.0e-15);
}

void check_mms(test::stage3::BodyKind kind) {
  const auto body = test::stage3::approved_body(kind);
  const auto sample =
      test::stage3::evaluate_mms(body, {0.21, 0.28, 0.37}, 0.013);
  HUNDUN_CHECK(test::stage3::finite(sample));
  HUNDUN_CHECK(std::abs(sample.velocity_m_per_s.x) +
                   std::abs(sample.velocity_m_per_s.y) +
                   std::abs(sample.velocity_m_per_s.z) >
               1.0e-14);
  const double divergence = sample.velocity_gradient_per_s[0][0] +
                            sample.velocity_gradient_per_s[1][1] +
                            sample.velocity_gradient_per_s[2][2];
  HUNDUN_CHECK(std::abs(divergence) <= 2.0e-12);
}

void check_exact_wall_zeros() {
  using test::stage3::BodyKind;
  const auto sphere = test::stage3::approved_body(BodyKind::sphere);
  const auto sphere_wall =
      test::stage3::evaluate_mms(sphere, {0.68, 0.5, 0.5}, 0.017);
  HUNDUN_CHECK(test::stage3::max_abs(sphere_wall.velocity_m_per_s) <= 2.0e-13);

  for (const auto kind : {BodyKind::sphere, BodyKind::finite_cylinder,
                          BodyKind::oblique_rectangular_prism}) {
    const auto outer = test::stage3::evaluate_mms(
        test::stage3::approved_body(kind), {0.0, 0.37, 0.61}, 0.017);
    HUNDUN_CHECK(test::stage3::max_abs(outer.velocity_m_per_s) <= 2.0e-13);
  }
}

void check_outer_pressure_neumann_compatibility() {
  const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
  constexpr double tolerance = 2.0e-13;
  constexpr std::array<runtime::Real3, 6> points{
      runtime::Real3{0.0, 0.37, 0.61}, runtime::Real3{1.0, 0.37, 0.61},
      runtime::Real3{0.29, 0.0, 0.61}, runtime::Real3{0.29, 1.0, 0.61},
      runtime::Real3{0.29, 0.37, 0.0}, runtime::Real3{0.29, 0.37, 1.0}};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    for (std::size_t side = 0U; side < 2U; ++side) {
      const auto sample =
          test::stage3::evaluate_mms(body, points[2U * axis + side], 0.013);
      const std::array<double, 3> gradient{
          sample.mechanical_pressure_gradient_pa_per_m.x,
          sample.mechanical_pressure_gradient_pa_per_m.y,
          sample.mechanical_pressure_gradient_pa_per_m.z};
      HUNDUN_CHECK(std::abs(gradient[axis]) <= tolerance);
      HUNDUN_CHECK(
          std::abs(sample.mechanical_pressure_hessian_pa_per_m2[axis][axis]) <=
          tolerance);
    }
  }
}

void check_surface_and_force_oracle(test::stage3::BodyKind kind) {
  const auto body = test::stage3::approved_body(kind);
  constexpr double h = 1.0 / 24.0;
  const auto coarse = test::stage3::make_manufactured_surface(body, h);
  const auto fine = test::stage3::make_manufactured_surface(body, 0.5 * h);
  HUNDUN_CHECK(!coarse.triangles.empty());
  HUNDUN_CHECK(!fine.triangles.empty());
  HUNDUN_CHECK(coarse.maximum_edge_m <= 0.45 * h * (1.0 + 1.0e-12));
  HUNDUN_CHECK(fine.maximum_edge_m <= 0.45 * 0.5 * h * (1.0 + 1.0e-12));
  if (kind != test::stage3::BodyKind::oblique_rectangular_prism) {
    HUNDUN_CHECK(coarse.maximum_chord_error_m > 0.0);
    HUNDUN_CHECK(fine.maximum_chord_error_m > 0.0);
    HUNDUN_CHECK(
        std::log(coarse.maximum_chord_error_m / fine.maximum_chord_error_m) /
            std::log(2.0) >=
        1.8);
  } else {
    HUNDUN_CHECK(coarse.maximum_chord_error_m == 0.0);
    HUNDUN_CHECK(fine.maximum_chord_error_m == 0.0);
  }

  const auto reference =
      test::stage3::analytic_force_reference(body, 0.013, 48U);
  const auto refined = test::stage3::analytic_force_reference(body, 0.013, 96U);
  const double pressure_scale = reference.surface_area_m2;
  const double viscous_scale =
      test::stage3::kDynamicViscosityPaS * reference.surface_area_m2;
  const double total_scale = std::max(pressure_scale, viscous_scale);
  constexpr double discrimination =
      256.0 * std::numeric_limits<double>::epsilon();
  HUNDUN_CHECK(refined.pressure_traction_rms_force_N >
               discrimination * pressure_scale);
  HUNDUN_CHECK(refined.viscous_traction_rms_force_N >
               discrimination * viscous_scale);
  HUNDUN_CHECK(refined.total_traction_rms_force_N >
               discrimination * total_scale);
  HUNDUN_CHECK(test::stage3::max_abs_difference(reference.force.pressure_N,
                                                refined.force.pressure_N) <=
               1.0e-13 * std::max(1.0, pressure_scale));
  HUNDUN_CHECK(test::stage3::max_abs_difference(reference.force.viscous_N,
                                                refined.force.viscous_N) <=
               1.0e-13 * std::max(1.0, viscous_scale));
  HUNDUN_CHECK(test::stage3::max_abs_difference(reference.force.total_N,
                                                refined.force.total_N) <=
               1.0e-13 * std::max(1.0, total_scale));
  HUNDUN_CHECK(std::abs(reference.pressure_traction_rms_force_N -
                        refined.pressure_traction_rms_force_N) <=
               1.0e-13 * std::max(1.0, pressure_scale));
  HUNDUN_CHECK(std::abs(reference.viscous_traction_rms_force_N -
                        refined.viscous_traction_rms_force_N) <=
               1.0e-13 * std::max(1.0, viscous_scale));
  HUNDUN_CHECK(std::abs(reference.total_traction_rms_force_N -
                        refined.total_traction_rms_force_N) <=
               1.0e-13 * std::max(1.0, total_scale));
}

test::stage3::ManufacturedCase formal_case(const std::string &selector,
                                           int cells);

void check_signed_force_reference_preflight() {
  constexpr double time_s = 0.013;
  const double epsilon = std::numeric_limits<double>::epsilon();
  const double sqrt_epsilon = std::sqrt(epsilon);
  const double force_scale =
      test::stage3::kDynamicViscosityPaS *
      test::stage3::kReferenceVelocityMPerS /
      test::stage3::kReferenceLengthM;

  const auto increment =
      test::stage3::absolute_viscous_traction_force_increment(
          {-2.0, 3.0, -4.0}, 0.25);
  const bool pure_increment_exact =
      increment.x == 0.5 && increment.y == 0.75 && increment.z == 1.0;

  test::stage3::AnalyticForceReference valid;
  valid.surface_area_m2 = 1.0;
  valid.viscous_traction_rms_force_N = force_scale / std::sqrt(3.0);
  valid.viscous_absolute_component_traction_force_N =
      {0.125 * force_scale, 0.125 * force_scale, 0.25 * force_scale};
  valid.viscous_traction_l1_force_N = 0.5 * force_scale;
  const double valid_force =
      std::max(2.0 * 1024.0 * sqrt_epsilon * force_scale,
               8192.0 * epsilon * force_scale);
  valid.force.viscous_N = {0.0, valid_force, 0.0};
  const bool valid_control_accepted =
      test::stage3::signed_viscous_force_reference_preflight_accepts(valid,
                                                                     valid);

  auto nonfinite_coarse = valid;
  auto nonfinite_refined = valid;
  nonfinite_coarse.force.viscous_N =
      {0.0, std::numeric_limits<double>::infinity(), 0.0};
  nonfinite_refined.force.viscous_N = nonfinite_coarse.force.viscous_N;
  const bool s0_nf_rejected =
      !test::stage3::signed_viscous_force_reference_preflight_accepts(
          nonfinite_coarse, nonfinite_refined);

  auto missing_component_coarse = valid;
  auto missing_component_refined = valid;
  missing_component_coarse.viscous_absolute_component_traction_force_N =
      {0.0, 0.125 * force_scale, 0.25 * force_scale};
  missing_component_coarse.viscous_traction_l1_force_N =
      (0.0 + 0.125 * force_scale) + 0.25 * force_scale;
  missing_component_refined = missing_component_coarse;
  const bool s0_n1a_rejected =
      !test::stage3::signed_viscous_force_reference_preflight_accepts(
          missing_component_coarse, missing_component_refined);

  auto nonfinite_l1_coarse = valid;
  auto nonfinite_l1_refined = valid;
  const double large_finite =
      0.75 * std::numeric_limits<double>::max();
  nonfinite_l1_coarse.viscous_absolute_component_traction_force_N =
      {large_finite, large_finite, force_scale};
  nonfinite_l1_coarse.viscous_traction_l1_force_N =
      (large_finite + large_finite) + force_scale;
  nonfinite_l1_refined = nonfinite_l1_coarse;
  const bool s0_n1b_rejected =
      !test::stage3::signed_viscous_force_reference_preflight_accepts(
          nonfinite_l1_coarse, nonfinite_l1_refined);

  auto l2_composition_coarse = valid;
  auto l2_composition_refined = valid;
  const auto valid_absolute =
      valid.viscous_absolute_component_traction_force_N;
  l2_composition_coarse.viscous_traction_l1_force_N = std::sqrt(
      valid_absolute.x * valid_absolute.x +
      valid_absolute.y * valid_absolute.y +
      valid_absolute.z * valid_absolute.z);
  l2_composition_refined = l2_composition_coarse;
  const bool s0_n1c_l2_rejected =
      !test::stage3::signed_viscous_force_reference_preflight_accepts(
          l2_composition_coarse, l2_composition_refined);

  auto x_composition_coarse = valid;
  auto x_composition_refined = valid;
  x_composition_coarse.viscous_traction_l1_force_N = valid_absolute.x;
  x_composition_refined = x_composition_coarse;
  const bool s0_n1c_x_rejected =
      !test::stage3::signed_viscous_force_reference_preflight_accepts(
          x_composition_coarse, x_composition_refined);

  auto roundoff_floor_coarse = valid;
  auto roundoff_floor_refined = valid;
  const double roundoff_force = 2048.0 * epsilon * force_scale;
  const double roundoff_certificate =
      roundoff_force / (2.0 * 1024.0 * sqrt_epsilon);
  const double roundoff_l1 = 0.5 * roundoff_certificate;
  roundoff_floor_coarse.force.viscous_N = {0.0, roundoff_force, 0.0};
  roundoff_floor_coarse.viscous_traction_rms_force_N =
      roundoff_certificate / std::sqrt(3.0);
  roundoff_floor_coarse.viscous_absolute_component_traction_force_N =
      {0.25 * roundoff_l1, 0.25 * roundoff_l1, 0.5 * roundoff_l1};
  roundoff_floor_coarse.viscous_traction_l1_force_N = roundoff_l1;
  roundoff_floor_refined = roundoff_floor_coarse;
  const bool s0_n2_rejected =
      !test::stage3::signed_viscous_force_reference_preflight_accepts(
          roundoff_floor_coarse, roundoff_floor_refined);

  auto cancellation_coarse = valid;
  auto cancellation_refined = valid;
  const double cancellation_force =
      0.75 * 1024.0 * sqrt_epsilon * force_scale;
  cancellation_coarse.force.viscous_N =
      {0.0, cancellation_force, 0.0};
  cancellation_refined = cancellation_coarse;
  const bool s0_n3_rejected =
      !test::stage3::signed_viscous_force_reference_preflight_accepts(
          cancellation_coarse, cancellation_refined);

  auto vector_stability_coarse = valid;
  auto vector_stability_refined = valid;
  const double stability_bound = 1.0e-13 * std::max(1.0, force_scale);
  vector_stability_coarse.force.viscous_N.y += 8.0 * stability_bound;
  const bool s0_n4_rejected =
      !test::stage3::signed_viscous_force_reference_preflight_accepts(
          vector_stability_coarse, vector_stability_refined);

  auto certificate_stability_coarse = valid;
  auto certificate_stability_refined = valid;
  certificate_stability_coarse.viscous_traction_rms_force_N =
      (force_scale + 1.5 * stability_bound) / std::sqrt(3.0);
  const bool s0_n5_rejected =
      !test::stage3::signed_viscous_force_reference_preflight_accepts(
          certificate_stability_coarse, certificate_stability_refined);

  constexpr std::array<const char *, 5> body_names{
      "sphere", "finite_cylinder", "force_certified_oblique_prism",
      "translated_sphere", "inside_sphere_cavity"};
  const std::array bodies{
      test::stage3::approved_body(test::stage3::BodyKind::sphere),
      test::stage3::approved_body(test::stage3::BodyKind::finite_cylinder),
      test::stage3::force_certified_oblique_prism(),
      test::stage3::translated_sphere(),
      test::stage3::approved_body(
          test::stage3::BodyKind::inside_sphere_cavity)};
  struct RealBodyObservation final {
    test::stage3::AnalyticForceReference reference;
    bool wrapper_accepted{};
    bool real_l1_valid{};
    bool omitted_separated{};
    bool reversed_separated{};
    double force_scale{};
    double force_net{};
    double traction_certificate{};
    std::optional<double> cancellation_ratio;
    std::optional<double> certified_cancellation_ratio;
    std::string wrapper_error;
  };
  std::array<RealBodyObservation, 5> real_observations{};
  for (std::size_t body_index = 0U; body_index < bodies.size(); ++body_index) {
    auto &observation = real_observations[body_index];
    try {
      observation.reference =
          test::stage3::verified_force_reference_preflight(bodies[body_index],
                                                           time_s);
      observation.wrapper_accepted = true;
      const auto absolute =
          observation.reference.viscous_absolute_component_traction_force_N;
      const double traction_l1 =
          observation.reference.viscous_traction_l1_force_N;
      observation.real_l1_valid =
          std::isfinite(traction_l1) && traction_l1 > 0.0 &&
          traction_l1 == (absolute.x + absolute.y) + absolute.z;
      observation.force_scale =
          test::stage3::kDynamicViscosityPaS *
          test::stage3::kReferenceVelocityMPerS *
          observation.reference.surface_area_m2 /
          test::stage3::kReferenceLengthM;
      const auto force = observation.reference.force.viscous_N;
      const bool force_components_finite =
          std::isfinite(force.x) && std::isfinite(force.y) &&
          std::isfinite(force.z);
      const bool scale_valid = std::isfinite(observation.force_scale) &&
                               observation.force_scale > 0.0;
      if (force_components_finite && scale_valid) {
        const double omitted_error = test::stage3::max_abs_difference(
                                         {0.0, 0.0, 0.0}, force) /
                                     observation.force_scale;
        const runtime::Real3 reversed{-force.x, -force.y, -force.z};
        const double reversed_error =
            test::stage3::max_abs_difference(reversed, force) /
            observation.force_scale;
        observation.omitted_separated =
            omitted_error >= 4096.0 * epsilon;
        observation.reversed_separated =
            reversed_error >= 8192.0 * epsilon;
      }
      observation.force_net = test::stage3::max_abs(force);
      observation.traction_certificate =
          std::sqrt(3.0) *
          observation.reference.viscous_traction_rms_force_N;
      const bool ratio_primitives_valid =
          observation.real_l1_valid && force_components_finite && scale_valid &&
          std::isfinite(observation.force_net) &&
          std::isfinite(observation.traction_certificate) &&
          observation.traction_certificate > 0.0;
      if (ratio_primitives_valid) {
        observation.cancellation_ratio =
            observation.force_net /
            observation.reference.viscous_traction_l1_force_N;
        observation.certified_cancellation_ratio =
            observation.force_net / observation.traction_certificate;
      }
    } catch (const std::exception &error) {
      observation.wrapper_error = error.what();
    } catch (...) {
      observation.wrapper_error = "non-standard exception";
    }
  }

  const auto focused_prism = bodies[2];
  const auto formal_prism = formal_case("prism_uniform_acceptance", 24).body;
  const bool prism_fixture_route_exact =
      focused_prism.kind == formal_prism.kind &&
      focused_prism.centre_m.x == formal_prism.centre_m.x &&
      focused_prism.centre_m.y == formal_prism.centre_m.y &&
      focused_prism.centre_m.z == formal_prism.centre_m.z &&
      focused_prism.radius_m == formal_prism.radius_m &&
      focused_prism.length_m == formal_prism.length_m &&
      focused_prism.axis.x == formal_prism.axis.x &&
      focused_prism.axis.y == formal_prism.axis.y &&
      focused_prism.axis.z == formal_prism.axis.z &&
      focused_prism.half_lengths_m.x == formal_prism.half_lengths_m.x &&
      focused_prism.half_lengths_m.y == formal_prism.half_lengths_m.y &&
      focused_prism.half_lengths_m.z == formal_prism.half_lengths_m.z &&
      focused_prism.prism_mms_factor_multiplier ==
          formal_prism.prism_mms_factor_multiplier;

  for (std::size_t body_index = 0U; body_index < body_names.size();
       ++body_index) {
    const auto &observation = real_observations[body_index];
    const auto force = observation.reference.force.viscous_N;
    const auto absolute =
        observation.reference.viscous_absolute_component_traction_force_N;
    std::cerr << "signed_force_reference_preflight body="
              << body_names[body_index] << " force=" << force.x << ','
              << force.y << ',' << force.z << " A_abs=" << absolute.x << ','
              << absolute.y << ',' << absolute.z
              << " T_abs="
              << observation.reference.viscous_traction_l1_force_N
              << " RMS="
              << observation.reference.viscous_traction_rms_force_N
              << " T_cert=" << observation.traction_certificate
              << " scale=" << observation.force_scale << " D_net=";
    if (observation.cancellation_ratio.has_value())
      std::cerr << *observation.cancellation_ratio;
    else
      std::cerr << "not-evaluated";
    std::cerr << " D_cert=";
    if (observation.certified_cancellation_ratio.has_value())
      std::cerr << *observation.certified_cancellation_ratio;
    else
      std::cerr << "not-evaluated";
    if (!observation.wrapper_error.empty())
      std::cerr << " wrapper_error=" << observation.wrapper_error;
    std::cerr << '\n';
  }

  std::cerr << "signed_force_reference_preflight observations"
            << " pure_increment_exact=" << pure_increment_exact
            << " valid_control_accepted=" << valid_control_accepted
            << " S0-NF_rejected=" << s0_nf_rejected
            << " S0-N1a_rejected=" << s0_n1a_rejected
            << " S0-N1b_rejected=" << s0_n1b_rejected
            << " S0-N1c-L2_rejected=" << s0_n1c_l2_rejected
            << " S0-N1c-X_rejected=" << s0_n1c_x_rejected
            << " S0-N2_rejected=" << s0_n2_rejected
            << " S0-N3_rejected=" << s0_n3_rejected
            << " S0-N4_rejected=" << s0_n4_rejected
            << " S0-N5_rejected=" << s0_n5_rejected;
  bool accepted = pure_increment_exact && valid_control_accepted &&
                  s0_nf_rejected && s0_n1a_rejected && s0_n1b_rejected &&
                  s0_n1c_l2_rejected && s0_n1c_x_rejected &&
                  s0_n2_rejected && s0_n3_rejected && s0_n4_rejected &&
                  s0_n5_rejected;
  for (std::size_t body_index = 0U; body_index < body_names.size();
       ++body_index) {
    const auto &observation = real_observations[body_index];
    std::cerr << ' ' << body_names[body_index]
              << "_wrapper_accepted=" << observation.wrapper_accepted << ' '
              << body_names[body_index]
              << "_real_l1_valid=" << observation.real_l1_valid << ' '
              << body_names[body_index]
              << "_omitted_separated=" << observation.omitted_separated << ' '
              << body_names[body_index]
              << "_reversed_separated=" << observation.reversed_separated;
    accepted = accepted && observation.wrapper_accepted &&
               observation.real_l1_valid && observation.omitted_separated &&
               observation.reversed_separated;
  }
  std::cerr << " prism_fixture_route_exact=" << prism_fixture_route_exact;
  accepted = accepted && prism_fixture_route_exact;
  std::cerr << '\n';
  HUNDUN_CHECK(accepted);
}

void check_source_quadrature(const runtime::MpiContext &mpi) {
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {4, 4, 4}, {false, false, false},
      runtime::DecompositionOptions{runtime::Int3{1, 1, 1}});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}});
  const auto cell =
      topology.find_local_cell(topology.global_cell_id({1, 1, 1}));
  HUNDUN_CHECK(cell.has_value());
  const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
  const auto average =
      test::stage3::source_cell_average(topology, geometry, *cell, body, 0.013);
  HUNDUN_CHECK(std::isfinite(average.x));
  HUNDUN_CHECK(std::isfinite(average.y));
  HUNDUN_CHECK(std::isfinite(average.z));
  runtime::Real3 gauss_reference{};
  constexpr double root = 0.774596669241483377035853079956;
  constexpr std::array<double, 3> nodes{-root, 0.0, root};
  constexpr std::array<double, 3> weights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  for (std::size_t k = 0U; k < nodes.size(); ++k)
    for (std::size_t j = 0U; j < nodes.size(); ++j)
      for (std::size_t i = 0U; i < nodes.size(); ++i) {
        const runtime::Real3 point{0.375 + 0.125 * nodes[i],
                                   0.375 + 0.125 * nodes[j],
                                   0.375 + 0.125 * nodes[k]};
        const auto value =
            test::stage3::evaluate_mms(body, point, 0.013).body_source_N_per_m3;
        const double weight = 0.125 * weights[i] * weights[j] * weights[k];
        gauss_reference.x += weight * value.x;
        gauss_reference.y += weight * value.y;
        gauss_reference.z += weight * value.z;
      }
  HUNDUN_CHECK(test::stage3::max_abs_difference(average, gauss_reference) <=
               2.0e-14);
  const auto centre =
      test::stage3::evaluate_mms(body, {0.375, 0.375, 0.375}, 0.013)
          .body_source_N_per_m3;
  HUNDUN_CHECK(test::stage3::max_abs_difference(average, centre) > 1.0e-7);
  const auto reference = test::stage3::evaluate_cell_average(
      topology, geometry, *cell, body, 0.013);
  geometry.require_compatible(topology);
  const auto validated_reference =
      test::stage3::detail::evaluate_cell_average_validated(topology, geometry,
                                                            *cell, body, 0.013);
  HUNDUN_CHECK(reference.velocity_m_per_s.x ==
               validated_reference.velocity_m_per_s.x);
  HUNDUN_CHECK(reference.velocity_m_per_s.y ==
               validated_reference.velocity_m_per_s.y);
  HUNDUN_CHECK(reference.velocity_m_per_s.z ==
               validated_reference.velocity_m_per_s.z);
  HUNDUN_CHECK(reference.mechanical_pressure_pa ==
               validated_reference.mechanical_pressure_pa);
  HUNDUN_CHECK(reference.body_source_N_per_m3.x ==
               validated_reference.body_source_N_per_m3.x);
  HUNDUN_CHECK(reference.body_source_N_per_m3.y ==
               validated_reference.body_source_N_per_m3.y);
  HUNDUN_CHECK(reference.body_source_N_per_m3.z ==
               validated_reference.body_source_N_per_m3.z);
  const auto direct = test::stage3::evaluate_cell_average(topology, geometry,
                                                          *cell, body, -0.007);
  const auto rescaled = rescale_manufactured_average(reference, 0.013, -0.007);
  HUNDUN_CHECK(test::stage3::max_abs_difference(direct.velocity_m_per_s,
                                                rescaled.velocity_m_per_s) <=
               2.0e-15);
  HUNDUN_CHECK_NEAR(direct.mechanical_pressure_pa,
                    rescaled.mechanical_pressure_pa, 2.0e-15);

  mesh::MeshGeometry warped_geometry(
      topology, mesh::AnalyticWarpedBoxMapping{
                    {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {0.02, -0.015, 0.01}});
  const auto warped = test::stage3::evaluate_cell_average(
      topology, warped_geometry, *cell, body, 0.013);
  const auto warped_repeat = test::stage3::evaluate_cell_average(
      topology, warped_geometry, *cell, body, 0.013);
  HUNDUN_CHECK(test::stage3::finite(test::stage3::evaluate_mms(
      body, warped_geometry.cell_center_m(*cell), 0.013)));
  HUNDUN_CHECK(warped.velocity_m_per_s.x == warped_repeat.velocity_m_per_s.x);
  HUNDUN_CHECK(warped.velocity_m_per_s.y == warped_repeat.velocity_m_per_s.y);
  HUNDUN_CHECK(warped.velocity_m_per_s.z == warped_repeat.velocity_m_per_s.z);
  HUNDUN_CHECK(warped.mechanical_pressure_pa ==
               warped_repeat.mechanical_pressure_pa);
  HUNDUN_CHECK(warped.body_source_N_per_m3.x ==
               warped_repeat.body_source_N_per_m3.x);
  HUNDUN_CHECK(warped.body_source_N_per_m3.y ==
               warped_repeat.body_source_N_per_m3.y);
  HUNDUN_CHECK(warped.body_source_N_per_m3.z ==
               warped_repeat.body_source_N_per_m3.z);
  HUNDUN_CHECK(test::stage3::max_abs_difference(warped.body_source_N_per_m3,
                                                average) > 1.0e-7);
}

double independent_face_flux_reference(const mesh::MeshTopology &topology,
                                       const mesh::MeshGeometry &geometry,
                                       mesh::LocalFaceId face,
                                       const test::stage3::BodySpec &body,
                                       double time_s) {
  const auto logical = topology.logical_face(face);
  const auto c = logical.coordinate;
  std::array<runtime::Int3, 4> coordinates{};
  switch (logical.axis) {
  case mesh::FaceAxis::x:
    coordinates = {
        runtime::Int3{c.x, c.y, c.z}, runtime::Int3{c.x, c.y + 1, c.z},
        runtime::Int3{c.x, c.y + 1, c.z + 1}, runtime::Int3{c.x, c.y, c.z + 1}};
    break;
  case mesh::FaceAxis::y:
    coordinates = {
        runtime::Int3{c.x, c.y, c.z}, runtime::Int3{c.x, c.y, c.z + 1},
        runtime::Int3{c.x + 1, c.y, c.z + 1}, runtime::Int3{c.x + 1, c.y, c.z}};
    break;
  case mesh::FaceAxis::z:
    coordinates = {
        runtime::Int3{c.x, c.y, c.z}, runtime::Int3{c.x + 1, c.y, c.z},
        runtime::Int3{c.x + 1, c.y + 1, c.z}, runtime::Int3{c.x, c.y + 1, c.z}};
    break;
  }
  std::array<runtime::Real3, 4> vertices{};
  for (std::size_t vertex = 0U; vertex < vertices.size(); ++vertex)
    vertices[vertex] = geometry.vertex_position_m(coordinates[vertex]);
  const auto add = [](runtime::Real3 left, runtime::Real3 right) {
    return runtime::Real3{left.x + right.x, left.y + right.y, left.z + right.z};
  };
  const auto subtract = [](runtime::Real3 left, runtime::Real3 right) {
    return runtime::Real3{left.x - right.x, left.y - right.y, left.z - right.z};
  };
  const auto multiply = [](double scale, runtime::Real3 value) {
    return runtime::Real3{scale * value.x, scale * value.y, scale * value.z};
  };
  const auto dot = [](runtime::Real3 left, runtime::Real3 right) {
    return left.x * right.x + left.y * right.y + left.z * right.z;
  };
  const auto cross = [](runtime::Real3 left, runtime::Real3 right) {
    return runtime::Real3{left.y * right.z - left.z * right.y,
                          left.z * right.x - left.x * right.z,
                          left.x * right.y - left.y * right.x};
  };
  constexpr std::array<double, 8> nodes{
      0.019855071751231884, 0.10166676129318664, 0.2372337950418355,
      0.4082826787521751,   0.5917173212478249,  0.7627662049581645,
      0.8983332387068134,   0.9801449282487681};
  constexpr std::array<double, 8> weights{
      0.05061426814518815, 0.11119051722668725, 0.15685332293894365,
      0.181341891689181,   0.181341891689181,   0.15685332293894365,
      0.11119051722668725, 0.05061426814518815};
  double positive_flux = 0.0;
  runtime::Real3 integrated_area{};
  const auto accumulate = [&](runtime::Real3 point,
                              runtime::Real3 area_weight) {
    const auto sample = test::stage3::evaluate_mms(body, point, time_s);
    HUNDUN_CHECK(test::stage3::finite(sample));
    positive_flux += dot(sample.velocity_m_per_s, area_weight);
    integrated_area = add(integrated_area, area_weight);
  };
  if (geometry.mapping_kind() == mesh::MappingKind::uniform_box) {
    const auto edge_a = subtract(vertices[1], vertices[0]);
    const auto edge_b = subtract(vertices[3], vertices[0]);
    const auto area = cross(edge_a, edge_b);
    for (std::size_t a = 0U; a < nodes.size(); ++a)
      for (std::size_t b = 0U; b < nodes.size(); ++b) {
        const auto point = add(vertices[0], add(multiply(nodes[a], edge_a),
                                                multiply(nodes[b], edge_b)));
        accumulate(point, multiply(weights[a] * weights[b], area));
      }
  } else {
    constexpr std::array<std::array<std::size_t, 3>, 2> triangles{
        std::array<std::size_t, 3>{0U, 1U, 2U},
        std::array<std::size_t, 3>{0U, 2U, 3U}};
    for (const auto &triangle : triangles) {
      const auto origin = vertices[triangle[0]];
      const auto edge_a = subtract(vertices[triangle[1]], origin);
      const auto edge_b = subtract(vertices[triangle[2]], origin);
      const auto area = cross(edge_a, edge_b);
      for (std::size_t a = 0U; a < nodes.size(); ++a)
        for (std::size_t b = 0U; b < nodes.size(); ++b) {
          const double one_minus_a = 1.0 - nodes[a];
          const auto point =
              add(origin, add(multiply(nodes[a], edge_a),
                              multiply(one_minus_a * nodes[b], edge_b)));
          accumulate(point,
                     multiply(weights[a] * weights[b] * one_minus_a, area));
        }
    }
  }
  const auto owner_area =
      geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
  const double orientation = dot(integrated_area, owner_area);
  HUNDUN_CHECK(std::isfinite(positive_flux));
  HUNDUN_CHECK(std::isfinite(orientation));
  HUNDUN_CHECK(orientation != 0.0);
  return orientation > 0.0 ? positive_flux : -positive_flux;
}

void check_face_history_quadrature(const runtime::MpiContext &mpi) {
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {6, 6, 6}, {false, false, false},
      runtime::DecompositionOptions{runtime::Int3{1, 1, 1}});
  mesh::MeshTopology topology(decomposition);
  const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
  const runtime::Int3 global_cell{1, 2, 3};
  const auto local_cell =
      topology.find_local_cell(topology.global_cell_id(global_cell));
  HUNDUN_CHECK(local_cell.has_value());
  const std::array<mesh::LogicalFace, 6> logical_faces{
      mesh::LogicalFace{mesh::FaceAxis::x, global_cell},
      mesh::LogicalFace{mesh::FaceAxis::x,
                        {global_cell.x + 1, global_cell.y, global_cell.z}},
      mesh::LogicalFace{mesh::FaceAxis::y, global_cell},
      mesh::LogicalFace{mesh::FaceAxis::y,
                        {global_cell.x, global_cell.y + 1, global_cell.z}},
      mesh::LogicalFace{mesh::FaceAxis::z, global_cell},
      mesh::LogicalFace{mesh::FaceAxis::z,
                        {global_cell.x, global_cell.y, global_cell.z + 1}}};

  for (const auto mapping : {test::stage3::ManufacturedMapping::uniform,
                             test::stage3::ManufacturedMapping::warped}) {
    mesh::MeshGeometry geometry =
        mapping == test::stage3::ManufacturedMapping::uniform
            ? mesh::MeshGeometry(
                  topology,
                  mesh::UniformBoxMapping{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}})
            : mesh::MeshGeometry(topology,
                                 mesh::AnalyticWarpedBoxMapping{
                                     {0.0, 0.0, 0.0},
                                     {1.0, 1.0, 1.0},
                                     {test::stage3::kWarpAmplitude[0],
                                      test::stage3::kWarpAmplitude[1],
                                      test::stage3::kWarpAmplitude[2]}});
    double integrated_divergence = 0.0;
    double midpoint_divergence = 0.0;
    bool midpoint_mutation_observed = false;
    for (const auto logical_face : logical_faces) {
      const auto face =
          topology.find_local_face(topology.global_face_id(logical_face));
      HUNDUN_CHECK(face.has_value());
      const auto history = test::stage3::evaluate_face_history(
          topology, geometry, *face, body, 0.013);
      const auto repeat = test::stage3::evaluate_face_history(
          topology, geometry, *face, body, 0.013);
      HUNDUN_CHECK(test::fp64_bits(history.velocity_m_per_s.x) ==
                   test::fp64_bits(repeat.velocity_m_per_s.x));
      HUNDUN_CHECK(test::fp64_bits(history.velocity_m_per_s.y) ==
                   test::fp64_bits(repeat.velocity_m_per_s.y));
      HUNDUN_CHECK(test::fp64_bits(history.velocity_m_per_s.z) ==
                   test::fp64_bits(repeat.velocity_m_per_s.z));
      HUNDUN_CHECK(test::fp64_bits(history.owner_normal_volume_flux_m3_per_s) ==
                   test::fp64_bits(repeat.owner_normal_volume_flux_m3_per_s));

      const auto owner_area =
          geometry.face_area_vector_m2(*face, mesh::FaceSide::owner);
      const double represented_flux =
          history.velocity_m_per_s.x * owner_area.x +
          history.velocity_m_per_s.y * owner_area.y +
          history.velocity_m_per_s.z * owner_area.z;
      const double flux_tolerance =
          64.0 * std::numeric_limits<double>::epsilon() *
          std::max({1.0, std::abs(represented_flux),
                    std::abs(history.owner_normal_volume_flux_m3_per_s)});
      HUNDUN_CHECK(std::abs(represented_flux -
                            history.owner_normal_volume_flux_m3_per_s) <=
                   flux_tolerance);
      const double independent_flux = independent_face_flux_reference(
          topology, geometry, *face, body, 0.013);
      HUNDUN_CHECK(std::abs(independent_flux -
                            history.owner_normal_volume_flux_m3_per_s) <=
                   512.0 * std::numeric_limits<double>::epsilon() *
                       std::max(1.0, std::abs(independent_flux)));

      const auto neighbour = topology.neighbour(*face);
      HUNDUN_CHECK(neighbour.has_value());
      const auto neighbour_area =
          geometry.face_area_vector_m2(*face, mesh::FaceSide::neighbour);
      HUNDUN_CHECK(test::stage3::max_abs_difference(
                       neighbour_area,
                       {-owner_area.x, -owner_area.y, -owner_area.z}) == 0.0);
      const double neighbour_flux =
          history.velocity_m_per_s.x * neighbour_area.x +
          history.velocity_m_per_s.y * neighbour_area.y +
          history.velocity_m_per_s.z * neighbour_area.z;
      HUNDUN_CHECK(std::abs(neighbour_flux +
                            history.owner_normal_volume_flux_m3_per_s) <=
                   flux_tolerance);

      double outward_sign = 0.0;
      if (topology.owner(*face) == *local_cell)
        outward_sign = 1.0;
      else {
        HUNDUN_CHECK(*neighbour == *local_cell);
        outward_sign = -1.0;
      }
      integrated_divergence +=
          outward_sign * history.owner_normal_volume_flux_m3_per_s;
      const auto midpoint = test::stage3::evaluate_mms(
          body, geometry.face_center_m(*face), 0.013);
      const double midpoint_owner_flux =
          midpoint.velocity_m_per_s.x * owner_area.x +
          midpoint.velocity_m_per_s.y * owner_area.y +
          midpoint.velocity_m_per_s.z * owner_area.z;
      midpoint_mutation_observed =
          midpoint_mutation_observed ||
          std::abs(midpoint_owner_flux - independent_flux) >
              1024.0 * std::numeric_limits<double>::epsilon();
      const runtime::Real3 correction{
          history.velocity_m_per_s.x - midpoint.velocity_m_per_s.x,
          history.velocity_m_per_s.y - midpoint.velocity_m_per_s.y,
          history.velocity_m_per_s.z - midpoint.velocity_m_per_s.z};
      const runtime::Real3 tangential_change{
          correction.y * owner_area.z - correction.z * owner_area.y,
          correction.z * owner_area.x - correction.x * owner_area.z,
          correction.x * owner_area.y - correction.y * owner_area.x};
      HUNDUN_CHECK(test::stage3::max_abs(tangential_change) <=
                   128.0 * std::numeric_limits<double>::epsilon());
      midpoint_divergence += outward_sign * midpoint_owner_flux;
    }
    HUNDUN_CHECK(midpoint_mutation_observed);
    HUNDUN_CHECK(std::abs(integrated_divergence) <= 2.0e-14);
    HUNDUN_CHECK(std::abs(midpoint_divergence) > 1.0e-12);
    HUNDUN_CHECK(std::abs(integrated_divergence) <
                 1.0e-3 * std::abs(midpoint_divergence));
  }
}

config::FlowCaseConfig manufactured_config(int cells) {
  config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "stage3-task11-manufactured";
  config.simulation_type = config::SimulationType::variable_density_flow;
  config.density_model = config::DensityModel::constant;
  config.resources.expected_ranks = 1;
  config.resources.process_grid = runtime::Int3{1, 1, 1};
  config.mesh.cells = {cells, cells, cells};
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.physics.rho_ref_kg_per_m3 = test::stage3::kReferenceDensityKgPerM3;
  config.physics.dynamic_viscosity_pa_s = test::stage3::kDynamicViscosityPaS;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t patch = 0U; patch < names.size(); ++patch) {
    config.boundaries[patch].patch = names[patch];
    config.boundaries[patch].type = config::BoundaryType::no_slip_wall;
  }
  return config;
}

runtime::FieldDescriptor manufactured_cell_field(const char *name,
                                                 std::uint32_t components) {
  return {name,
          "1",
          "stage3_task11",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          4,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor manufactured_face_field(const char *name,
                                                 std::uint32_t components) {
  return {name,
          "1",
          "stage3_task11",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

double canonical_zero(double value) {
  return std::abs(value) <= 1.0e-30 ? 0.0 : value;
}

struct ManufacturedRunResult final {
  double velocity_l2{};
  double velocity_linf{};
  double near_wall_velocity_l2{};
  double bulk_velocity_l2{};
  double pressure_l2{};
  double pressure_linf{};
  double near_wall_pressure_l2{};
  double bulk_pressure_l2{};
  double exact_velocity_linf{};
  double numerical_velocity_linf{};
  double linf_wall_distance_m{};
  runtime::Int3 linf_global_cell{};
  std::uint32_t linf_component{};
  double linf_signed_error{};
  double linf_exact_value{};
  double linf_numerical_value{};
  double linf_pressure_error_pa{};
  std::size_t linf_incident_wall_links{};
  std::array<double, 3> component_velocity_linf{};
  std::array<double, 4> wall_band_velocity_linf{};
  std::optional<flow::ForceAttemptReport> force;
};

double observed_order(double coarse, double fine) {
  HUNDUN_CHECK(std::isfinite(coarse));
  HUNDUN_CHECK(std::isfinite(fine));
  HUNDUN_CHECK(coarse > fine);
  HUNDUN_CHECK(fine > 0.0);
  return std::log(coarse / fine) / std::log(2.0);
}

constexpr double kFormalErrorFloor =
    4096.0 * std::numeric_limits<double>::epsilon();

bool strict_two_segment_order_passes(
    const std::array<double, 3> &errors) noexcept {
  for (const double error : errors)
    if (!std::isfinite(error) || !(error > 0.0))
      return false;
  if (!(errors[0] > errors[1] && errors[1] > errors[2]))
    return false;
  return errors[2] > kFormalErrorFloor &&
         observed_order(errors[0], errors[1]) >= 1.8 &&
         observed_order(errors[1], errors[2]) >= 1.8;
}

bool positive_decreasing_fast_pair(
    const std::array<double, 2> &errors) noexcept {
  return std::isfinite(errors[0]) && std::isfinite(errors[1]) &&
         errors[0] > errors[1] && errors[1] > 0.0;
}

constexpr std::size_t kPressureForceRow = 8U;
constexpr std::size_t kViscousTractionL2Row = 9U;
constexpr std::size_t kViscousForceRow = 10U;
constexpr std::size_t kTotalForceRow = 11U;
constexpr std::size_t kPressureConsistencyRow = 12U;
constexpr std::size_t kViscousConsistencyRow = 13U;
constexpr std::size_t kTotalConsistencyRow = 14U;

constexpr bool formal_order_is_required(std::size_t row) noexcept {
  return row <= kTotalConsistencyRow;
}

bool order_oracle_is_mutation_sensitive() {
  const double h = 1.0 / 12.0;
  const double quadratic_coarse = h * h;
  const double quadratic_fine = 0.25 * quadratic_coarse;
  const double disabled_quadratic_coarse = h;
  const double disabled_quadratic_fine = 0.5 * disabled_quadratic_coarse;
  const std::array<double, 3> endpoint_only_pressure_linf{
      4.28627e-3, 9.02527e-4, 3.23720e-4};
  const std::array<double, 3> nondecreasing_pressure_linf{
      4.28627e-3, 9.02527e-4, 9.02527e-4};
  const std::array<double, 3> nonfinite_pressure_linf{
      4.28627e-3, std::numeric_limits<double>::infinity(), 3.23720e-4};
  const std::array<double, 3> nonpositive_pressure_linf{4.28627e-3, 9.02527e-4,
                                                        0.0};
  const std::array<double, 3> low_endpoint_order_pressure_linf{4.0e-3, 1.0e-3,
                                                               4.0e-4};
  const std::array<double, 3> rejected_r20_pressure_linf{2.02187e-3, 3.78413e-4,
                                                         2.20932e-4};
  const std::array<double, 3> rejected_fine_segment_pressure_linf{16.0, 2.0,
                                                                  1.0};
  const std::array<double, 3> accepted_two_segment{16.0, 4.0, 1.0};
  const std::array<double, 3> rejected_first_segment{8.0, 4.0, 1.0};
  const std::array<double, 3> rejected_second_segment{16.0, 4.0, 2.0};
  const std::array<double, 3> rejected_subfloor_two_segment{
      8.0 * kFormalErrorFloor, 2.0 * kFormalErrorFloor,
      0.5 * kFormalErrorFloor};
  const std::array<double, 2> accepted_fast_pair{4.0, 1.0};
  const std::array<double, 2> nondecreasing_fast_pair{1.0, 1.0};
  const std::array<double, 2> nonfinite_fast_pair{
      1.0, std::numeric_limits<double>::infinity()};
  const std::array<double, 2> nonpositive_fast_pair{1.0, 0.0};
  return observed_order(quadratic_coarse, quadratic_fine) >= 1.8 &&
         observed_order(disabled_quadratic_coarse, disabled_quadratic_fine) <
             1.8 &&
         !strict_two_segment_order_passes(endpoint_only_pressure_linf) &&
         !strict_two_segment_order_passes(nondecreasing_pressure_linf) &&
         !strict_two_segment_order_passes(nonfinite_pressure_linf) &&
         !strict_two_segment_order_passes(nonpositive_pressure_linf) &&
         !strict_two_segment_order_passes(low_endpoint_order_pressure_linf) &&
         !strict_two_segment_order_passes(rejected_r20_pressure_linf) &&
         !strict_two_segment_order_passes(
             rejected_fine_segment_pressure_linf) &&
         strict_two_segment_order_passes(accepted_two_segment) &&
         !strict_two_segment_order_passes(rejected_first_segment) &&
         !strict_two_segment_order_passes(rejected_second_segment) &&
         !strict_two_segment_order_passes(rejected_subfloor_two_segment) &&
         positive_decreasing_fast_pair(accepted_fast_pair) &&
         !positive_decreasing_fast_pair(nondecreasing_fast_pair) &&
         !positive_decreasing_fast_pair(nonfinite_fast_pair) &&
         !positive_decreasing_fast_pair(nonpositive_fast_pair) &&
         formal_order_is_required(kPressureForceRow) &&
         formal_order_is_required(kViscousTractionL2Row) &&
         formal_order_is_required(kViscousForceRow) &&
         formal_order_is_required(kTotalForceRow) &&
         formal_order_is_required(kPressureConsistencyRow) &&
         formal_order_is_required(kViscousConsistencyRow) &&
         formal_order_is_required(kTotalConsistencyRow) &&
         !formal_order_is_required(kTotalConsistencyRow + 1U);
}

bool pressure_error_extremum_oracle_is_mutation_sensitive() {
  using test::stage3::PressureErrorExtremum;
  if (test::stage3::ManufacturedCase{}.surface_policy !=
      test::stage3::ManufacturedSurfacePolicy::per_level)
    return false;
  const std::vector<PressureErrorExtremum> candidates{
      {-2.0, 2.0, 42U, {2, 1, 0}, 0.125}, {1.0, 1.0, 55U, {3, 1, 0}, 0.25}};
  const auto selected =
      test::stage3::select_pressure_error_extremum(candidates);
  const auto repeated =
      test::stage3::select_pressure_error_extremum(candidates);
  if (!selected.has_value() || !repeated.has_value())
    return false;
  const auto top_two = test::stage3::select_pressure_error_extrema(candidates,
                                                                    2U);
  if (top_two.size() != 2U || top_two[0].global_cell_id != 42U ||
      top_two[1].global_cell_id != 55U)
    return false;
  if (!test::stage3::select_pressure_error_extrema(candidates, 0U).empty())
    return false;
  if (!test::stage3::select_pressure_error_extrema({}, 2U).empty())
    return false;
  const auto same_bits = [](double left, double right) noexcept {
    return std::memcmp(&left, &right, sizeof(double)) == 0;
  };
  if (selected->global_cell_id != 42U || selected->logical_cell.x != 2 ||
      selected->logical_cell.y != 1 || selected->logical_cell.z != 0 ||
      !same_bits(selected->signed_error_pa, -2.0) ||
      !same_bits(selected->absolute_error_pa, 2.0) ||
      !same_bits(selected->wall_distance_m, 0.125) ||
      !same_bits(selected->signed_error_pa, repeated->signed_error_pa) ||
      !same_bits(selected->absolute_error_pa, repeated->absolute_error_pa) ||
      selected->global_cell_id != repeated->global_cell_id ||
      selected->logical_cell.x != repeated->logical_cell.x ||
      selected->logical_cell.y != repeated->logical_cell.y ||
      selected->logical_cell.z != repeated->logical_cell.z ||
      !same_bits(selected->wall_distance_m, repeated->wall_distance_m))
    return false;

  auto larger_error = candidates;
  larger_error[1] = {3.0, 3.0, 55U, {3, 1, 0}, 0.25};
  const auto larger_selected =
      test::stage3::select_pressure_error_extremum(larger_error);
  if (!larger_selected.has_value() || larger_selected->global_cell_id != 55U)
    return false;

  auto smaller_id_tie = candidates;
  smaller_id_tie[1] = {2.0, 2.0, 7U, {1, 0, 0}, 0.375};
  const auto tie_selected =
      test::stage3::select_pressure_error_extremum(smaller_id_tie);
  if (!tie_selected.has_value() || tie_selected->global_cell_id != 7U)
    return false;

  auto nonfinite_error = candidates;
  nonfinite_error[1].signed_error_pa = std::numeric_limits<double>::infinity();
  nonfinite_error[1].absolute_error_pa =
      std::numeric_limits<double>::infinity();
  auto negative_distance = candidates;
  negative_distance[1].wall_distance_m = -1.0;
  auto nonfinite_distance = candidates;
  nonfinite_distance[1].wall_distance_m =
      std::numeric_limits<double>::quiet_NaN();
  return !test::stage3::select_pressure_error_extremum({}).has_value() &&
         !test::stage3::select_pressure_error_extremum(nonfinite_error)
              .has_value() &&
         !test::stage3::select_pressure_error_extremum(negative_distance)
              .has_value() &&
         !test::stage3::select_pressure_error_extremum(nonfinite_distance)
              .has_value();
}

bool near_wall_pressure_band_oracle_is_mutation_sensitive() {
  using test::stage3::normalized_near_wall_pressure_band;
  const auto requires_band = [&](double distance_over_h,
                                 std::size_t expected) {
    const auto actual = normalized_near_wall_pressure_band(distance_over_h);
    return actual.has_value() && *actual == expected;
  };
  return !normalized_near_wall_pressure_band(-0.01).has_value() &&
         requires_band(0.0, 0U) && requires_band(0.5, 0U) &&
         requires_band(0.500001, 1U) && requires_band(1.0, 1U) &&
         requires_band(1.000001, 2U) && requires_band(1.5, 2U) &&
         requires_band(1.500001, 3U) && requires_band(2.0, 3U) &&
         !normalized_near_wall_pressure_band(2.000001).has_value() &&
         !normalized_near_wall_pressure_band(
              std::numeric_limits<double>::quiet_NaN())
              .has_value();
}

bool fixed_physical_near_wall_band_oracle_is_mutation_sensitive() {
  using test::stage3::manufactured_near_wall_band_contains;
  constexpr double thickness = test::stage3::kReferenceLengthM / 6.0;
  return !manufactured_near_wall_band_contains(-0.01) &&
         manufactured_near_wall_band_contains(0.0) &&
         manufactured_near_wall_band_contains(thickness) &&
         !manufactured_near_wall_band_contains(
             std::nextafter(thickness,
                            std::numeric_limits<double>::infinity())) &&
         !manufactured_near_wall_band_contains(
              std::numeric_limits<double>::quiet_NaN());
}

bool pressure_error_moments_oracle_is_mutation_sensitive() {
  test::stage3::PressureErrorMoments moments;
  if (!test::stage3::accumulate_pressure_error_moments(moments, 1.0, 1.0) ||
      !test::stage3::accumulate_pressure_error_moments(moments, 3.0, 3.0))
    return false;
  const auto statistics =
      test::stage3::summarize_pressure_error_moments(moments);
  if (!statistics.has_value())
    return false;
  const double tolerance = 32.0 * std::numeric_limits<double>::epsilon();
  if (std::abs(statistics->mean_error_pa - 2.5) > tolerance ||
      std::abs(statistics->rms_error_pa - std::sqrt(7.0)) > tolerance ||
      std::abs(statistics->centered_rms_error_pa - std::sqrt(0.75)) >
          tolerance ||
      statistics->volume_m3 != 4.0 || statistics->cell_count != 2U)
    return false;
  const auto before = moments;
  if (test::stage3::accumulate_pressure_error_moments(moments, 4.0, 0.0) ||
      moments.signed_error_volume_pa_m3 !=
          before.signed_error_volume_pa_m3 ||
      moments.squared_error_volume_pa2_m3 !=
          before.squared_error_volume_pa2_m3 ||
      moments.volume_m3 != before.volume_m3 ||
      moments.cell_count != before.cell_count)
    return false;
  return !test::stage3::summarize_pressure_error_moments({}).has_value();
}

std::array<double, 15>
formal_errors(const test::stage3::ManufacturedRunResult &result) {
  return {
      result.errors.velocity_l2,           result.errors.velocity_linf,
      result.errors.near_wall_velocity_l2, result.errors.pressure_l2,
      result.errors.pressure_linf,         result.errors.near_wall_pressure_l2,
      result.errors.penetration_l2,        result.errors.penetration_linf,
      result.errors.pressure_force,        result.errors.viscous_traction_l2,
      result.errors.viscous_force,         result.errors.total_force,
      result.errors.pressure_consistency,  result.errors.viscous_consistency,
      result.errors.total_consistency};
}

runtime::Real3 select_operator_viscous_force(
    const test::stage3::ManufacturedRunResult &result) noexcept {
  return result.operator_force.viscous_N;
}

constexpr std::array<const char *, 15> kFormalErrorNames{
    "velocity_l2",          "velocity_linf",       "near_wall_velocity_l2",
    "pressure_l2",          "pressure_linf",       "near_wall_pressure_l2",
    "penetration_l2",       "penetration_linf",    "pressure_force",
    "viscous_traction_l2",  "viscous_force",       "total_force",
    "pressure_consistency", "viscous_consistency", "total_consistency"};

test::stage3::ManufacturedCase formal_case(const std::string &selector,
                                           int cells) {
  using test::stage3::BodyKind;
  using test::stage3::ManufacturedMapping;
  test::stage3::ManufacturedCase result;
  result.cells = cells;
  result.process_grid = {1, 1, 1};
  result.mapping = selector.find("_warped_") != std::string::npos
                       ? ManufacturedMapping::warped
                       : ManufacturedMapping::uniform;
  result.fluid_side = config::ImmersedFluidSide::outside;
  if (selector.rfind("sphere_", 0U) == 0U)
    result.body = test::stage3::approved_body(BodyKind::sphere);
  else if (selector.rfind("cylinder_", 0U) == 0U)
    result.body = test::stage3::approved_body(BodyKind::finite_cylinder);
  else if (selector.rfind("prism_", 0U) == 0U)
    result.body = test::stage3::force_certified_oblique_prism();
  else if (selector.rfind("translated_sphere_", 0U) == 0U)
    result.body = test::stage3::translated_sphere();
  else if (selector.rfind("inside_cavity_", 0U) == 0U) {
    result.body = test::stage3::approved_body(BodyKind::inside_sphere_cavity);
    result.fluid_side = config::ImmersedFluidSide::inside;
  } else {
    throw std::runtime_error("unknown manufactured acceptance selector");
  }
  return result;
}

void print_formal_result(const std::string &selector, int cells,
                         const test::stage3::ManufacturedRunResult &result) {
  const auto errors = formal_errors(result);
  std::cerr << selector << " cells=" << cells;
  for (std::size_t row = 0U; row < errors.size(); ++row)
    std::cerr << ' ' << kFormalErrorNames[row] << '=' << errors[row];
  std::cerr << " gauge_mean=" << result.active_pressure_mean
            << " active_cells=" << result.active_cell_count
            << " links=" << result.immersed_link_count
            << " wall_points=" << result.wall_point_count
            << " pressure_linf_signed_error="
            << result.pressure_error_extremum.signed_error_pa
            << " pressure_linf_absolute_error="
            << result.pressure_error_extremum.absolute_error_pa
            << " pressure_linf_global_cell="
            << result.pressure_error_extremum.global_cell_id
            << " pressure_linf_logical_cell="
            << result.pressure_error_extremum.logical_cell.x << ','
            << result.pressure_error_extremum.logical_cell.y << ','
            << result.pressure_error_extremum.logical_cell.z
            << " pressure_linf_wall_distance="
            << result.pressure_error_extremum.wall_distance_m
            << " pressure_extremum_authority_link="
            << result.pressure_extremum_authority.nearest_link
            << " pressure_extremum_authority_owner="
            << result.pressure_extremum_authority.owner_rank
            << " pressure_extremum_authority_direct="
            << result.pressure_extremum_authority.direct_fluid_cell
            << " pressure_extremum_authority_distance="
            << result.pressure_extremum_authority.cell_to_wall_distance_m
            << " pressure_extremum_authority_condition="
            << result.pressure_extremum_authority.condition_estimate
            << " pressure_extremum_authority_donors="
            << result.pressure_extremum_authority.donor_count
            << " pressure_extremum_row_links="
            << result.pressure_extremum_authority.row_link_count
            << " pressure_extremum_authority_fingerprint="
            << result.pressure_extremum_authority.donor_fingerprint
            << " pressure_extremum_value_amplification="
            << result.pressure_extremum_authority.pressure_value_amplification
            << " pressure_extremum_gradient_amplification="
            << result.pressure_extremum_authority.normal_gradient_amplification
            << '\n';
}

void run_formal_sequence(const std::string &selector,
                         const runtime::MpiContext &mpi,
                         std::array<int, 3> grids) {
  HUNDUN_CHECK(mpi.size() == 1);
  std::array<test::stage3::ManufacturedRunResult, 3> results;
  for (std::size_t level = 0U; level < grids.size(); ++level) {
    results[level] = test::stage3::run_manufactured_case(
        mpi, formal_case(selector, grids[level]));
    print_formal_result(selector, grids[level], results[level]);
    const double h = 1.0 / static_cast<double>(grids[level]);
    HUNDUN_CHECK(results[level].surface_maximum_edge_m <=
                 0.45 * h * (1.0 + 1.0e-12));
    if (selector.rfind("inside_cavity_", 0U) == 0U)
      HUNDUN_CHECK(std::abs(results[level].active_pressure_mean) <= 1.0e-12);
  }

  const bool planar_surface = selector.rfind("prism_", 0U) == 0U;
  if (planar_surface) {
    for (const auto &result : results)
      HUNDUN_CHECK(result.surface_maximum_chord_error_m == 0.0);
  } else {
    for (const auto &result : results)
      HUNDUN_CHECK(result.surface_maximum_chord_error_m > 0.0);
    HUNDUN_CHECK(observed_order(results[0].surface_maximum_chord_error_m,
                                results[1].surface_maximum_chord_error_m) >=
                 1.8);
    HUNDUN_CHECK(observed_order(results[1].surface_maximum_chord_error_m,
                                results[2].surface_maximum_chord_error_m) >=
                 1.8);
  }

  const auto coarse = formal_errors(results[0]);
  const auto medium = formal_errors(results[1]);
  const auto fine = formal_errors(results[2]);
  constexpr double floor = kFormalErrorFloor;
  constexpr double penetration_exact =
      8192.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, test::stage3::kReferenceVelocityMPerS);
  for (std::size_t row = 0U; row < coarse.size(); ++row) {
    HUNDUN_CHECK(std::isfinite(coarse[row]));
    HUNDUN_CHECK(std::isfinite(medium[row]));
    HUNDUN_CHECK(std::isfinite(fine[row]));
    const bool penetration = row == 6U || row == 7U;
    if (penetration && coarse[row] <= penetration_exact &&
        medium[row] <= penetration_exact && fine[row] <= penetration_exact) {
      std::cerr << selector << ' ' << kFormalErrorNames[row]
                << " exact_enforcement\n";
      continue;
    }
    if (row == 4U) {
      HUNDUN_CHECK(fine[row] > floor);
      const double first = observed_order(coarse[row], medium[row]);
      const double second = observed_order(medium[row], fine[row]);
      const double endpoint = std::log(coarse[row] / fine[row]) / std::log(4.0);
      std::cerr << selector << ' ' << kFormalErrorNames[row]
                << " orders=" << first << ',' << second
                << " endpoint_order=" << endpoint << '\n';
      HUNDUN_CHECK(strict_two_segment_order_passes(
          {coarse[row], medium[row], fine[row]}));
      continue;
    }
    HUNDUN_CHECK(formal_order_is_required(row));
    if (row == kViscousTractionL2Row) {
      const double first = observed_order(coarse[row], medium[row]);
      const double second = observed_order(medium[row], fine[row]);
      std::cerr << selector << ' ' << kFormalErrorNames[row]
                << " orders=" << first << ',' << second << '\n';
      HUNDUN_CHECK(strict_two_segment_order_passes(
          {coarse[row], medium[row], fine[row]}));
      continue;
    }
    HUNDUN_CHECK(formal_order_is_required(row));
    HUNDUN_CHECK(coarse[row] > medium[row]);
    HUNDUN_CHECK(medium[row] > fine[row]);
    HUNDUN_CHECK(fine[row] > floor);
    const double first = observed_order(coarse[row], medium[row]);
    const double second = observed_order(medium[row], fine[row]);
    std::cerr << selector << ' ' << kFormalErrorNames[row]
              << " orders=" << first << ',' << second << '\n';
    HUNDUN_CHECK(first >= 1.8);
    HUNDUN_CHECK(second >= 1.8);
  }
}

void run_functional_selection_fast(const runtime::MpiContext &mpi) {
  std::array<test::stage3::ManufacturedRunResult, 2> results;
  constexpr std::array<int, 2> grids{12, 24};
  for (std::size_t level = 0U; level < grids.size(); ++level) {
    results[level] = test::stage3::run_manufactured_case(
        mpi, formal_case("sphere_warped_acceptance", grids[level]));
    print_formal_result("functional_selection_fast", grids[level],
                        results[level]);
  }
  const auto coarse = formal_errors(results[0]);
  const auto fine = formal_errors(results[1]);
  for (const std::size_t row :
       {kPressureForceRow, kViscousTractionL2Row, kViscousForceRow,
        kTotalForceRow, kPressureConsistencyRow, kViscousConsistencyRow,
        kTotalConsistencyRow}) {
    const double order = observed_order(coarse[row], fine[row]);
    std::cerr << "functional_selection_fast " << kFormalErrorNames[row]
              << " order=" << order << '\n';
    HUNDUN_CHECK(order >= 1.8);
  }
}

std::optional<double> diagnostic_order(double coarse, double fine,
                                       int coarse_cells, int fine_cells);

void run_exact_force_consistency_fast(const runtime::MpiContext &mpi) {
  std::array<test::stage3::ManufacturedRunResult, 2> results;
  constexpr std::array<int, 2> grids{12, 24};
  for (std::size_t level = 0U; level < grids.size(); ++level) {
    results[level] = test::stage3::run_manufactured_case(
        mpi, formal_case("sphere_uniform_acceptance", grids[level]));
    print_formal_result("exact_force_consistency_fast", grids[level],
                        results[level]);
  }
  const auto norm = [](runtime::Real3 value) {
    return test::stage3::max_abs(value);
  };
  for (std::size_t level = 0U; level < grids.size(); ++level) {
    const auto &result = results[level];
    const double scale = result.pressure_force_scale_N;
    const runtime::Real3 exact_consistency{
        -result.pressure_measure.exact_state_a22_consistency_N.x,
        -result.pressure_measure.exact_state_a22_consistency_N.y,
        -result.pressure_measure.exact_state_a22_consistency_N.z};
    const runtime::Real3 exact_donor_consistency{
        -result.pressure_measure.exact_state_a22_donor_reaction_N.x,
        -result.pressure_measure.exact_state_a22_donor_reaction_N.y,
        -result.pressure_measure.exact_state_a22_donor_reaction_N.z};
    const runtime::Real3 exact_wall_consistency{
        -result.pressure_measure.exact_state_a22_wall_reaction_N.x -
            result.analytic_force.pressure_N.x,
        -result.pressure_measure.exact_state_a22_wall_reaction_N.y -
            result.analytic_force.pressure_N.y,
        -result.pressure_measure.exact_state_a22_wall_reaction_N.z -
            result.analytic_force.pressure_N.z};
    const auto &solved_consistency = result.consistency.pressure_N;
    const runtime::Real3 state_consistency{
        solved_consistency.x - exact_consistency.x,
        solved_consistency.y - exact_consistency.y,
        solved_consistency.z - exact_consistency.z};
    std::cerr << "exact_force_consistency_fast cells=" << grids[level]
              << " pressure_C_exact_norm=" << norm(exact_consistency) / scale
              << " pressure_C_solved_norm=" << norm(solved_consistency) / scale
              << " pressure_C_state_norm=" << norm(state_consistency) / scale
              << " pressure_donor_term_norm="
              << norm(exact_donor_consistency) / scale
              << " pressure_wall_term_norm="
              << norm(exact_wall_consistency) / scale
              << " pressure_state_defect_norm="
              << norm(result.pressure_measure.numerical_state_a22_defect_N) /
                     scale
              << " surface_geometry_error_p="
              << norm(runtime::Real3{
                     result.wall_plan_analytic_force.pressure_N.x -
                         result.analytic_force.pressure_N.x,
                     result.wall_plan_analytic_force.pressure_N.y -
                         result.analytic_force.pressure_N.y,
                     result.wall_plan_analytic_force.pressure_N.z -
                         result.analytic_force.pressure_N.z}) /
                     scale
              << " surface_reconstruction_error_p="
              << norm(runtime::Real3{
                     result.surface_traction.pressure_N.x -
                         result.wall_plan_analytic_force.pressure_N.x,
                     result.surface_traction.pressure_N.y -
                         result.wall_plan_analytic_force.pressure_N.y,
                     result.surface_traction.pressure_N.z -
                         result.wall_plan_analytic_force.pressure_N.z}) /
                     scale
              << " viscous_C_solved_norm=" << norm(result.consistency.viscous_N)
              << " row_wall_value_linf="
              << result.pressure_measure.row_wall_value_linf
              << " row_wall_value_l2="
              << result.pressure_measure.row_wall_value_l2
              << " per_link_wall_value_linf="
              << result.pressure_measure.per_link_wall_value_linf
              << " per_link_wall_value_l2="
              << result.pressure_measure.per_link_wall_value_l2
              << " per_link_origin_g_linf="
              << result.pressure_measure.per_link_origin_g_wall_value_linf
              << " per_link_origin_g_l2="
              << result.pressure_measure.per_link_origin_g_wall_value_l2
              << " linear_wall_value_linf="
              << result.pressure_measure.linear_wall_value_linf
              << " linear_hg_linf="
              << result.pressure_measure.linear_wall_value_hg_linf
              << " plain_linear_linf="
              << result.pressure_measure.plain_linear_wall_value_linf
              << " linear_extrap_wall_linf="
              << result.pressure_measure.linear_extrap_wall_linf
              << " linear_extrap_wall_l2="
              << result.pressure_measure.linear_extrap_wall_l2
              << " constrained_origin_linear_linf="
              << result.pressure_measure.constrained_origin_linear_linf
              << " constrained_origin_hg_linf="
              << result.pressure_measure.constrained_origin_hg_linear_linf
              << " constant_field_linf="
              << result.pressure_measure.constant_field_linf
              << " abs_sum_w_t1="
              << result.pressure_measure.abs_sum_w_t1_linf
              << " abs_sum_w_t2="
              << result.pressure_measure.abs_sum_w_t2_linf
              << " abs_sum_w_n="
              << result.pressure_measure.abs_sum_w_n_linf
              << " authority_linear_linf="
              << result.pressure_measure.authority_linear_linf
              << " ghost_donor_wall_linf="
              << result.pressure_measure.ghost_donor_wall_linf
              << " ghost_donor_wall_l2="
              << result.pressure_measure.ghost_donor_wall_l2
              << " ghost_links="
              << result.pressure_measure.ghost_donor_link_count
              << " linear_row_residual_linf="
              << result.pressure_measure.linear_row_residual_linf
              << " linear_row_residual_l2="
              << result.pressure_measure.linear_row_residual_l2
              << " mms_row_residual_linf="
              << result.pressure_measure.mms_row_residual_linf
              << " mms_row_residual_l2="
              << result.pressure_measure.mms_row_residual_l2
              << " mms_rows="
              << result.pressure_measure.mms_row_count
              << " link_true_normal_mismatch_linf="
              << result.pressure_measure.link_true_normal_mismatch_linf
              << " mms_row_coherence="
              << result.pressure_measure.mms_row_coherence_ratio
              << " linear_row_sum="
              << result.pressure_measure.linear_row_sum_N.x << ','
              << result.pressure_measure.linear_row_sum_N.y << ','
              << result.pressure_measure.linear_row_sum_N.z
              << " single_link_sum="
              << result.pressure_measure.single_link_row_sum_N.x << ','
              << result.pressure_measure.single_link_row_sum_N.y << ','
              << result.pressure_measure.single_link_row_sum_N.z
              << " multi_link_sum="
              << result.pressure_measure.multi_link_row_sum_N.x << ','
              << result.pressure_measure.multi_link_row_sum_N.y << ','
              << result.pressure_measure.multi_link_row_sum_N.z
              << " single_link_donor_sum="
              << result.pressure_measure.single_link_donor_sum_N.x << ','
              << result.pressure_measure.single_link_donor_sum_N.y << ','
              << result.pressure_measure.single_link_donor_sum_N.z
              << " single_link_wall_sum="
              << result.pressure_measure.single_link_wall_sum_N.x << ','
              << result.pressure_measure.single_link_wall_sum_N.y << ','
              << result.pressure_measure.single_link_wall_sum_N.z
              << " single_link_linear_donor_sum="
              << result.pressure_measure.single_link_linear_donor_sum_N.x << ','
              << result.pressure_measure.single_link_linear_donor_sum_N.y << ','
              << result.pressure_measure.single_link_linear_donor_sum_N.z
              << " single_link_donor_constant_sum="
              << result.pressure_measure.single_link_donor_constant_sum_N.x
              << ',' << result.pressure_measure.single_link_donor_constant_sum_N.y
              << ',' << result.pressure_measure.single_link_donor_constant_sum_N.z
              << " single_link_donor_moment_x_sum="
              << result.pressure_measure.single_link_donor_moment_x_sum_N.x
              << ',' << result.pressure_measure.single_link_donor_moment_x_sum_N.y
              << ',' << result.pressure_measure.single_link_donor_moment_x_sum_N.z
              << " single_link_donor_moment_y_sum="
              << result.pressure_measure.single_link_donor_moment_y_sum_N.x
              << ',' << result.pressure_measure.single_link_donor_moment_y_sum_N.y
              << ',' << result.pressure_measure.single_link_donor_moment_y_sum_N.z
              << " single_link_donor_moment_z_sum="
              << result.pressure_measure.single_link_donor_moment_z_sum_N.x
              << ',' << result.pressure_measure.single_link_donor_moment_z_sum_N.y
              << ',' << result.pressure_measure.single_link_donor_moment_z_sum_N.z
              << " single_link_wall_coefficient_sum="
              << result.pressure_measure.single_link_wall_coefficient_sum_N.x
              << ',' << result.pressure_measure.single_link_wall_coefficient_sum_N.y
              << ',' << result.pressure_measure.single_link_wall_coefficient_sum_N.z
              << " single_link_donor_second_moment_sum="
              << result.pressure_measure.single_link_donor_second_moment_sum_N.x
              << ',' << result.pressure_measure.single_link_donor_second_moment_sum_N.y
              << ',' << result.pressure_measure.single_link_donor_second_moment_sum_N.z
              << " single_link_authority_wall_flux_sum="
              << result.pressure_measure.single_link_authority_wall_flux_sum_N.x
              << ',' << result.pressure_measure.single_link_authority_wall_flux_sum_N.y
              << ',' << result.pressure_measure.single_link_authority_wall_flux_sum_N.z
              << " single_link_exact_face_flux_sum="
              << result.pressure_measure.single_link_exact_face_flux_sum_N.x
              << ',' << result.pressure_measure.single_link_exact_face_flux_sum_N.y
              << ',' << result.pressure_measure.single_link_exact_face_flux_sum_N.z
              << " multi_link_exact_face_flux_sum="
              << result.pressure_measure.multi_link_exact_face_flux_sum_N.x
              << ',' << result.pressure_measure.multi_link_exact_face_flux_sum_N.y
              << ',' << result.pressure_measure.multi_link_exact_face_flux_sum_N.z
              << " single_link_authority_face_diff_l2="
              << result.pressure_measure.single_link_authority_value_difference_l2_N.x
              << ',' << result.pressure_measure.single_link_authority_value_difference_l2_N.y
              << ',' << result.pressure_measure.single_link_authority_value_difference_l2_N.z
              << " analytic_pressure_force="
              << result.analytic_force.pressure_N.x << ','
              << result.analytic_force.pressure_N.y << ','
              << result.analytic_force.pressure_N.z
              << " surface_pressure_traction="
              << result.surface_traction.pressure_N.x << ','
              << result.surface_traction.pressure_N.y << ','
              << result.surface_traction.pressure_N.z
              << " wall_plan_analytic_pressure="
              << result.wall_plan_analytic_force.pressure_N.x << ','
              << result.wall_plan_analytic_force.pressure_N.y << ','
              << result.wall_plan_analytic_force.pressure_N.z
              << " single_rows="
              << result.pressure_measure.single_link_row_count
              << " multi_rows="
              << result.pressure_measure.multi_link_row_count
              << " linear_wall_value_l2="
              << result.pressure_measure.linear_wall_value_l2
              << " linear_links="
              << result.pressure_measure.linear_wall_value_count
              << " authority_wall_value_linf="
              << result.pressure_measure.authority_wall_value_linf
              << " authority_wall_value_l2="
              << result.pressure_measure.authority_wall_value_l2
              << " authority_centroid_value_linf="
              << result.pressure_measure.authority_centroid_value_linf
              << " authority_centroid_value_l2="
              << result.pressure_measure.authority_centroid_value_l2
              << " exact_centroid_quadrature_sum="
              << result.pressure_measure.exact_centroid_quadrature_sum_N.x
              << ',' << result.pressure_measure.exact_centroid_quadrature_sum_N.y
              << ',' << result.pressure_measure.exact_centroid_quadrature_sum_N.z
              << " authority_centroid_quadrature_sum="
              << result.pressure_measure.authority_centroid_quadrature_sum_N.x
              << ',' << result.pressure_measure.authority_centroid_quadrature_sum_N.y
              << ',' << result.pressure_measure.authority_centroid_quadrature_sum_N.z
              << " surface_vector_sum="
              << result.pressure_measure.surface_vector_sum_N.x << ','
              << result.pressure_measure.surface_vector_sum_N.y << ','
              << result.pressure_measure.surface_vector_sum_N.z
              << " ghost_measure_sum="
              << result.pressure_measure.ghost_measure_sum_N.x << ','
              << result.pressure_measure.ghost_measure_sum_N.y << ','
              << result.pressure_measure.ghost_measure_sum_N.z
              << " ghost_authority_centroid_quadrature_sum="
              << result.pressure_measure.ghost_authority_centroid_quadrature_sum_N.x
              << ',' << result.pressure_measure.ghost_authority_centroid_quadrature_sum_N.y
              << ',' << result.pressure_measure.ghost_authority_centroid_quadrature_sum_N.z
              << " err_delta_pearson="
              << result.pressure_measure.error_delta_pearson
              << " err_oblique_pearson="
              << result.pressure_measure.error_obliqueness_pearson
              << " err_curvature_pearson="
              << result.pressure_measure.error_curvature_pearson
              << " err_delta_loglog_slope="
              << result.pressure_measure.error_delta_loglog_slope
              << " worst_link_obliqueness="
              << result.pressure_measure.worst_link_obliqueness
              << " wall_value_links="
              << result.pressure_measure.wall_value_link_count
              << '\n';
  }
  const auto exact_order = diagnostic_order(
      norm(results[0].pressure_measure.exact_state_a22_consistency_N) /
          results[0].pressure_force_scale_N,
      norm(results[1].pressure_measure.exact_state_a22_consistency_N) /
          results[1].pressure_force_scale_N,
      12, 24);
  const auto solved_order = diagnostic_order(
      norm(results[0].consistency.pressure_N) /
          results[0].pressure_force_scale_N,
      norm(results[1].consistency.pressure_N) /
          results[1].pressure_force_scale_N,
      12, 24);
  const auto defect_order = diagnostic_order(
      norm(results[0].pressure_measure.numerical_state_a22_defect_N) /
          results[0].pressure_force_scale_N,
      norm(results[1].pressure_measure.numerical_state_a22_defect_N) /
          results[1].pressure_force_scale_N,
      12, 24);
  const auto donor_term_order = diagnostic_order(
      norm(results[0].pressure_measure.exact_state_a22_donor_reaction_N) /
          results[0].pressure_force_scale_N,
      norm(results[1].pressure_measure.exact_state_a22_donor_reaction_N) /
          results[1].pressure_force_scale_N,
      12, 24);
  const auto wall_term_order = diagnostic_order(
      norm(results[0].pressure_measure.exact_state_a22_wall_reaction_N) /
          results[0].pressure_force_scale_N,
      norm(results[1].pressure_measure.exact_state_a22_wall_reaction_N) /
          results[1].pressure_force_scale_N,
      12, 24);
  const auto format = [](const std::optional<double> &order) {
    return order.has_value() ? std::to_string(*order)
                             : std::string("nonmonotone");
  };
  std::cerr << "exact_force_consistency_fast pressure C_exact_order="
            << format(exact_order) << " C_solved_order="
            << format(solved_order) << " state_defect_order="
            << format(defect_order) << " donor_term_order="
            << format(donor_term_order) << " wall_term_order="
            << format(wall_term_order) << '\n';
  const auto geometry_error_order = [&](std::size_t part) {
    const auto level_error = [&](std::size_t level) {
      const runtime::Real3 a = part == 0U
                                   ? results[level].wall_plan_analytic_force
                                         .pressure_N
                                   : results[level].wall_plan_analytic_force
                                         .viscous_N;
      const runtime::Real3 b =
          part == 0U ? results[level].analytic_force.pressure_N
                     : results[level].analytic_force.viscous_N;
      return norm(runtime::Real3{a.x - b.x, a.y - b.y, a.z - b.z});
    };
    return diagnostic_order(level_error(0U), level_error(1U), 12, 24);
  };
  const auto reconstruction_error_order = [&](std::size_t part) {
    const auto level_error = [&](std::size_t level) {
      const runtime::Real3 a = part == 0U
                                   ? results[level].surface_traction.pressure_N
                                   : results[level].surface_traction.viscous_N;
      const runtime::Real3 b =
          part == 0U
              ? results[level].wall_plan_analytic_force.pressure_N
              : results[level].wall_plan_analytic_force.viscous_N;
      return norm(runtime::Real3{a.x - b.x, a.y - b.y, a.z - b.z});
    };
    return diagnostic_order(level_error(0U), level_error(1U), 12, 24);
  };
  std::cerr << "exact_force_consistency_fast surface pressure_geometry_order="
            << format(geometry_error_order(0U))
            << " pressure_reconstruction_order="
            << format(reconstruction_error_order(0U))
            << " viscous_geometry_order=" << format(geometry_error_order(1U))
            << " viscous_reconstruction_order="
            << format(reconstruction_error_order(1U)) << '\n';
  const auto wall_value_order = [&](std::size_t which) {
    const auto level_error = [&](std::size_t level) {
      return which == 0U
                 ? results[level].pressure_measure.row_wall_value_l2
                 : which == 1U
                       ? results[level].pressure_measure.per_link_wall_value_l2
                       : results[level]
                             .pressure_measure.authority_wall_value_l2;
    };
    return diagnostic_order(level_error(0U), level_error(1U), 12, 24);
  };
  std::cerr << "exact_force_consistency_fast row_wall_value_l2_order="
            << format(wall_value_order(0U)) << " per_link_wall_value_l2_order="
            << format(wall_value_order(1U)) << " authority_wall_value_l2_order="
            << format(wall_value_order(2U)) << '\n';
}

std::optional<double> diagnostic_order(double coarse, double fine,
                                       int coarse_cells, int fine_cells) {
  if (!(coarse_cells > 0 && fine_cells > coarse_cells) ||
      !std::isfinite(coarse) || !std::isfinite(fine) || !(coarse > fine) ||
      !(fine > 0.0))
    return std::nullopt;
  const double value =
      std::log(coarse / fine) / std::log(static_cast<double>(fine_cells) /
                                         static_cast<double>(coarse_cells));
  return std::isfinite(value) ? std::optional<double>(value) : std::nullopt;
}

const char *surface_policy_name(
    test::stage3::ManufacturedSurfacePolicy policy) noexcept {
  return policy == test::stage3::ManufacturedSurfacePolicy::fixed_48
             ? "fixed_48_surface"
             : "per_level_surface";
}

void run_pressure_extrema_screen(const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(mpi.size() == 1);
  constexpr std::array<int, 3> grids{12, 24, 48};
  constexpr std::array<test::stage3::ManufacturedSurfacePolicy, 2> policies{
      test::stage3::ManufacturedSurfacePolicy::per_level,
      test::stage3::ManufacturedSurfacePolicy::fixed_48};
  for (const auto policy : policies) {
    std::array<test::stage3::ManufacturedRunResult, 3> results;
    for (std::size_t level = 0U; level < grids.size(); ++level) {
      auto definition = formal_case("sphere_uniform_acceptance", grids[level]);
      definition.collect_force = false;
      definition.surface_policy = policy;
      definition.collect_pressure_extrema_diagnostics = true;
      results[level] =
          test::stage3::run_manufactured_case(mpi, definition);
      const auto &result = results[level];
      std::cerr << "pressure_extrema_screen policy="
                << surface_policy_name(policy) << " cells=" << grids[level]
                << " pressure_l2=" << result.errors.pressure_l2
                << " near_wall_pressure_l2="
                << result.errors.near_wall_pressure_l2
                << " pressure_linf=" << result.errors.pressure_linf
                << " top_count="
                << result.pressure_error_extrema_top_k.size() << '\n';
      for (std::size_t rank = 0U;
           rank < result.pressure_error_extrema_top_k.size(); ++rank) {
        const auto &detail = result.pressure_error_extrema_top_k[rank];
        std::cerr << "pressure_extrema_screen policy="
                  << surface_policy_name(policy) << " cells=" << grids[level]
                  << " top=" << rank
                  << " cell=" << detail.error.logical_cell.x << ','
                  << detail.error.logical_cell.y << ','
                  << detail.error.logical_cell.z
                  << " global_id=" << detail.error.global_cell_id
                  << " signed_error=" << detail.error.signed_error_pa
                  << " absolute_error=" << detail.error.absolute_error_pa
                  << " wall_distance=" << detail.error.wall_distance_m
                  << " wall_distance_over_h="
                  << detail.wall_distance_over_h
                  << " incident_links=" << detail.incident_wall_link_count
                  << " authority_available=" << detail.authority_available;
        if (detail.authority_available)
          std::cerr << " authority_link=" << detail.authority.nearest_link
                    << " authority_direct="
                    << detail.authority.direct_fluid_cell
                    << " authority_distance="
                    << detail.authority.cell_to_wall_distance_m
                    << " authority_condition="
                    << detail.authority.condition_estimate
                    << " authority_donors=" << detail.authority.donor_count
                    << " authority_fingerprint="
                    << detail.authority.donor_fingerprint;
        std::cerr << '\n';
      }
    }
    const auto print_order = [&](const char *name, double coarse,
                                 double fine, int coarse_cells,
                                 int fine_cells) {
      const auto order = diagnostic_order(coarse, fine, coarse_cells,
                                           fine_cells);
      std::cerr << "pressure_extrema_screen policy="
                << surface_policy_name(policy) << " " << name
                << " order="
                << (order.has_value() ? std::to_string(*order)
                                      : std::string("nonmonotone"))
                << '\n';
    };
    print_order("pressure_l2_12_24", results[0].errors.pressure_l2,
                results[1].errors.pressure_l2, grids[0], grids[1]);
    print_order("pressure_l2_24_48", results[1].errors.pressure_l2,
                results[2].errors.pressure_l2, grids[1], grids[2]);
    print_order("near_wall_pressure_l2_12_24",
                results[0].errors.near_wall_pressure_l2,
                results[1].errors.near_wall_pressure_l2, grids[0], grids[1]);
    print_order("near_wall_pressure_l2_24_48",
                results[1].errors.near_wall_pressure_l2,
                results[2].errors.near_wall_pressure_l2, grids[1], grids[2]);
    print_order("pressure_linf_12_24", results[0].errors.pressure_linf,
                results[1].errors.pressure_linf, grids[0], grids[1]);
    print_order("pressure_linf_24_48", results[1].errors.pressure_linf,
                results[2].errors.pressure_linf, grids[1], grids[2]);
  }
}

void print_pressure_partition_statistics(
    int cells, const char *partition, std::size_t band,
    const test::stage3::PressureErrorMoments &moments) {
  const auto statistics =
      test::stage3::summarize_pressure_error_moments(moments);
  std::cerr << "near_wall_pressure_diagnostic cells=" << cells
            << " partition=" << partition << " band=" << band;
  if (!statistics.has_value()) {
    std::cerr << " empty\n";
    return;
  }
  std::cerr << " count=" << statistics->cell_count
            << " volume=" << statistics->volume_m3
            << " mean=" << statistics->mean_error_pa
            << " rms=" << statistics->rms_error_pa
            << " centered_rms=" << statistics->centered_rms_error_pa << '\n';
}

void run_near_wall_pressure_diagnostic(const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(mpi.size() == 1);
  constexpr std::array<int, 2> grids{12, 24};
  std::array<test::stage3::ManufacturedRunResult, 2> results;
  for (std::size_t level = 0U; level < grids.size(); ++level) {
    auto definition =
        formal_case("sphere_uniform_acceptance", grids[level]);
    definition.collect_force = false;
    results[level] =
        test::stage3::run_manufactured_case(mpi, definition);
    const auto &diagnostics = results[level].near_wall_pressure;
    HUNDUN_CHECK(diagnostics.available);
    std::cerr << "near_wall_pressure_diagnostic cells=" << grids[level]
              << " pressure_l2=" << results[level].errors.pressure_l2
              << " near_wall_pressure_l2="
              << results[level].errors.near_wall_pressure_l2 << '\n';
    print_pressure_partition_statistics(grids[level], "all", 4U,
                                        diagnostics.total);
    print_pressure_partition_statistics(grids[level], "incident", 4U,
                                        diagnostics.incident_total);
    print_pressure_partition_statistics(grids[level], "nonincident", 4U,
                                        diagnostics.nonincident_total);
    for (std::size_t band = 0U; band < 4U; ++band) {
      print_pressure_partition_statistics(grids[level], "all", band,
                                          diagnostics.distance_bands[band]);
      print_pressure_partition_statistics(
          grids[level], "incident", band,
          diagnostics.incident_distance_bands[band]);
      print_pressure_partition_statistics(
          grids[level], "nonincident", band,
          diagnostics.nonincident_distance_bands[band]);
    }
  }

  const auto print_partition_order = [&](const char *partition,
                                         std::size_t band,
                                         const auto &coarse,
                                         const auto &fine) {
    const auto coarse_statistics =
        test::stage3::summarize_pressure_error_moments(coarse);
    const auto fine_statistics =
        test::stage3::summarize_pressure_error_moments(fine);
    std::cerr << "near_wall_pressure_diagnostic partition=" << partition
              << " band=" << band;
    if (!coarse_statistics.has_value() || !fine_statistics.has_value()) {
      std::cerr << " order=unavailable\n";
      return;
    }
    const auto rms_order = diagnostic_order(
        coarse_statistics->rms_error_pa, fine_statistics->rms_error_pa,
        grids[0], grids[1]);
    const auto centered_order = diagnostic_order(
        coarse_statistics->centered_rms_error_pa,
        fine_statistics->centered_rms_error_pa, grids[0], grids[1]);
    const auto mean_order = diagnostic_order(
        std::abs(coarse_statistics->mean_error_pa),
        std::abs(fine_statistics->mean_error_pa), grids[0], grids[1]);
    const auto format = [](const std::optional<double> &value) {
      return value.has_value() ? std::to_string(*value)
                               : std::string("nonmonotone");
    };
    std::cerr << " rms_order=" << format(rms_order)
              << " centered_rms_order=" << format(centered_order)
              << " abs_mean_order=" << format(mean_order) << '\n';
  };

  const auto fixed_physical_order = diagnostic_order(
      results[0].errors.near_wall_pressure_l2,
      results[1].errors.near_wall_pressure_l2, grids[0], grids[1]);
  std::cerr << "near_wall_pressure_diagnostic fixed_physical_rms_order="
            << (fixed_physical_order.has_value()
                    ? std::to_string(*fixed_physical_order)
                    : std::string("nonmonotone"))
            << '\n';

  print_partition_order("all", 4U,
                        results[0].near_wall_pressure.total,
                        results[1].near_wall_pressure.total);
  print_partition_order("incident", 4U,
                        results[0].near_wall_pressure.incident_total,
                        results[1].near_wall_pressure.incident_total);
  print_partition_order("nonincident", 4U,
                        results[0].near_wall_pressure.nonincident_total,
                        results[1].near_wall_pressure.nonincident_total);
  for (std::size_t band = 0U; band < 4U; ++band) {
    print_partition_order("all", band,
                          results[0].near_wall_pressure.distance_bands[band],
                          results[1].near_wall_pressure.distance_bands[band]);
    print_partition_order(
        "incident", band,
        results[0].near_wall_pressure.incident_distance_bands[band],
        results[1].near_wall_pressure.incident_distance_bands[band]);
    print_partition_order(
        "nonincident", band,
        results[0].near_wall_pressure.nonincident_distance_bands[band],
        results[1].near_wall_pressure.nonincident_distance_bands[band]);
  }
}

void run_functional_closure_screen(const runtime::MpiContext &mpi,
                                   std::array<int, 2> grids,
                                   bool require_pressure_linf_order) {
  struct Probe final {
    const char *selector{};
    std::array<std::size_t, 8> rows{};
    std::size_t row_count{};
    bool pressure_linf{};
  };
  constexpr std::array<Probe, 2> probes{
      Probe{"prism_warped_acceptance",
            {kPressureForceRow, kViscousTractionL2Row, kViscousForceRow,
             kTotalForceRow, kPressureConsistencyRow, kViscousConsistencyRow,
             kTotalConsistencyRow, 0U},
            7U,
            false},
      Probe{"translated_sphere_warped_acceptance",
            {2U, kPressureForceRow, kViscousTractionL2Row, kViscousForceRow,
             kTotalForceRow, kPressureConsistencyRow, kViscousConsistencyRow,
             kTotalConsistencyRow},
            8U,
            true}};
  bool accepted = true;
  for (const auto &probe : probes) {
    std::array<test::stage3::ManufacturedRunResult, 2> results;
    for (std::size_t level = 0U; level < grids.size(); ++level) {
      results[level] = test::stage3::run_manufactured_case(
          mpi, formal_case(probe.selector, grids[level]));
      print_formal_result("functional_closure_fast", grids[level],
                          results[level]);
      const auto &surface = results[level].surface_traction.pressure_N;
      const auto &reference = results[level].analytic_force.pressure_N;
      std::cerr << probe.selector << " cells=" << grids[level]
                << " pressure_surface=" << surface.x << ',' << surface.y << ','
                << surface.z << " pressure_reference=" << reference.x << ','
                << reference.y << ',' << reference.z
                << " pressure_component_error=" << (surface.x - reference.x)
                << ',' << (surface.y - reference.y) << ','
                << (surface.z - reference.z)
                << " pressure_scale=" << results[level].pressure_force_scale_N
                << '\n';
      const auto &measure = results[level].pressure_measure;
      std::cerr << probe.selector << " cells=" << grids[level]
                << " pressure_measure_background_reaction="
                << measure.background_reaction_N.x << ','
                << measure.background_reaction_N.y << ','
                << measure.background_reaction_N.z
                << " pressure_measure_projected_reaction="
                << measure.projected_reaction_N.x << ','
                << measure.projected_reaction_N.y << ','
                << measure.projected_reaction_N.z
                << " pressure_measure_surface_partition_reaction="
                << measure.surface_partition_reaction_N.x << ','
                << measure.surface_partition_reaction_N.y << ','
                << measure.surface_partition_reaction_N.z
                << " pressure_measure_exact_state_a22_reaction="
                << measure.exact_state_a22_reaction_N.x << ','
                << measure.exact_state_a22_reaction_N.y << ','
                << measure.exact_state_a22_reaction_N.z
                << " pressure_measure_background_consistency="
                << measure.background_consistency_N.x << ','
                << measure.background_consistency_N.y << ','
                << measure.background_consistency_N.z
                << " pressure_measure_projected_consistency="
                << measure.projected_consistency_N.x << ','
                << measure.projected_consistency_N.y << ','
                << measure.projected_consistency_N.z
                << " pressure_measure_surface_partition_consistency="
                << measure.surface_partition_consistency_N.x << ','
                << measure.surface_partition_consistency_N.y << ','
                << measure.surface_partition_consistency_N.z
                << " pressure_measure_exact_state_a22_consistency="
                << measure.exact_state_a22_consistency_N.x << ','
                << measure.exact_state_a22_consistency_N.y << ','
                << measure.exact_state_a22_consistency_N.z
                << " pressure_measure_numerical_state_a22_defect="
                << measure.numerical_state_a22_defect_N.x << ','
                << measure.numerical_state_a22_defect_N.y << ','
                << measure.numerical_state_a22_defect_N.z
                << " pressure_measure_product_reproduction_linf="
                << measure.product_reproduction_linf_N
                << " pressure_measure_area_closure="
                << measure.background_area_closure << ','
                << measure.projected_area_closure << ','
                << measure.surface_partition_area_closure
                << " pressure_measure_links=" << measure.link_count << ','
                << measure.surface_supported_link_count << ','
                << measure.missing_surface_link_count << ','
                << measure.multiple_surface_point_link_count << '\n';
      const auto &scalar = results[level].pressure_scalar;
      std::cerr << probe.selector << " cells=" << grids[level]
                << " pressure_scalar_numerical_face_consistency="
                << scalar.numerical_face_consistency_N.x << ','
                << scalar.numerical_face_consistency_N.y << ','
                << scalar.numerical_face_consistency_N.z
                << " pressure_scalar_exact_face_consistency="
                << scalar.exact_face_consistency_N.x << ','
                << scalar.exact_face_consistency_N.y << ','
                << scalar.exact_face_consistency_N.z
                << " pressure_scalar_exact_wall_consistency="
                << scalar.exact_wall_consistency_N.x << ','
                << scalar.exact_wall_consistency_N.y << ','
                << scalar.exact_wall_consistency_N.z
                << " pressure_scalar_exact_hybrid_consistency="
                << scalar.exact_hybrid_consistency_N.x << ','
                << scalar.exact_hybrid_consistency_N.y << ','
                << scalar.exact_hybrid_consistency_N.z
                << " pressure_scalar_numerical_projected_defect_consistency="
                << scalar.numerical_projected_defect_consistency_N.x << ','
                << scalar.numerical_projected_defect_consistency_N.y << ','
                << scalar.numerical_projected_defect_consistency_N.z
                << " pressure_scalar_exact_projected_defect_consistency="
                << scalar.exact_projected_defect_consistency_N.x << ','
                << scalar.exact_projected_defect_consistency_N.y << ','
                << scalar.exact_projected_defect_consistency_N.z
                << " pressure_scalar_numerical_wall_projected_consistency="
                << scalar.numerical_wall_projected_consistency_N.x << ','
                << scalar.numerical_wall_projected_consistency_N.y << ','
                << scalar.numerical_wall_projected_consistency_N.z
                << " pressure_scalar_exact_wall_projected_consistency="
                << scalar.exact_wall_projected_consistency_N.x << ','
                << scalar.exact_wall_projected_consistency_N.y << ','
                << scalar.exact_wall_projected_consistency_N.z << '\n';
      const auto &viscous_surface = results[level].surface_traction.viscous_N;
      const auto &viscous_reference = results[level].analytic_force.viscous_N;
      const auto &viscous_plan_reference =
          results[level].wall_plan_analytic_force.viscous_N;
      const auto viscous_operator =
          select_operator_viscous_force(results[level]);
      std::cerr << probe.selector << " cells=" << grids[level]
                << " viscous_surface=" << viscous_surface.x << ','
                << viscous_surface.y << ',' << viscous_surface.z
                << " viscous_reference=" << viscous_reference.x << ','
                << viscous_reference.y << ',' << viscous_reference.z
                << " viscous_component_error="
                << (viscous_surface.x - viscous_reference.x) << ','
                << (viscous_surface.y - viscous_reference.y) << ','
                << (viscous_surface.z - viscous_reference.z)
                << " viscous_plan_reference=" << viscous_plan_reference.x << ','
                << viscous_plan_reference.y << ',' << viscous_plan_reference.z
                << " surface_plan_error="
                << (viscous_surface.x - viscous_plan_reference.x) << ','
                << (viscous_surface.y - viscous_plan_reference.y) << ','
                << (viscous_surface.z - viscous_plan_reference.z)
                << " plan_continuous_error="
                << (viscous_plan_reference.x - viscous_reference.x) << ','
                << (viscous_plan_reference.y - viscous_reference.y) << ','
                << (viscous_plan_reference.z - viscous_reference.z)
                << " plan_full_gradient="
                << results[level].wall_plan_full_gradient_viscous_force_N.x
                << ','
                << results[level].wall_plan_full_gradient_viscous_force_N.y
                << ','
                << results[level].wall_plan_full_gradient_viscous_force_N.z
                << " link_projected="
                << results[level].wall_link_projected_viscous_force_N.x << ','
                << results[level].wall_link_projected_viscous_force_N.y << ','
                << results[level].wall_link_projected_viscous_force_N.z
                << " link_full_gradient="
                << results[level].wall_link_full_gradient_viscous_force_N.x
                << ','
                << results[level].wall_link_full_gradient_viscous_force_N.y
                << ','
                << results[level].wall_link_full_gradient_viscous_force_N.z
                << " viscous_operator=" << viscous_operator.x << ','
                << viscous_operator.y << ',' << viscous_operator.z
                << " operator_surface_difference="
                << (viscous_operator.x - viscous_surface.x) << ','
                << (viscous_operator.y - viscous_surface.y) << ','
                << (viscous_operator.z - viscous_surface.z)
                << " traction_jump_same_link="
                << results[level].wall_traction_jump.same_link_rms << ','
                << results[level].wall_traction_jump.same_link_max << ','
                << results[level].wall_traction_jump.same_link_pairs
                << " traction_jump_cross_link_same_donor="
                << results[level]
                       .wall_traction_jump.cross_link_same_donor_rms
                << ','
                << results[level]
                       .wall_traction_jump.cross_link_same_donor_max
                << ','
                << results[level]
                       .wall_traction_jump.cross_link_same_donor_pairs
                << " traction_jump_cross_donor="
                << results[level].wall_traction_jump.cross_donor_rms << ','
                << results[level].wall_traction_jump.cross_donor_max << ','
                << results[level].wall_traction_jump.cross_donor_pairs << '\n';
      const auto &conditioning =
          results[level].wall_functional_conditioning;
      std::cerr << probe.selector << " cells=" << grids[level]
                << " wall_functional_amplification="
                << conditioning.normal_gradient_amplification_rms << ','
                << conditioning.normal_gradient_amplification_max
                << " wall_functional_condition="
                << conditioning.condition_estimate_max
                << " wall_signed_traction_error="
                << conditioning.signed_traction_error_N.x << ','
                << conditioning.signed_traction_error_N.y << ','
                << conditioning.signed_traction_error_N.z
                << " wall_absolute_component_error="
                << conditioning.absolute_component_traction_error_N.x << ','
                << conditioning.absolute_component_traction_error_N.y << ','
                << conditioning.absolute_component_traction_error_N.z
                << " wall_absolute_traction_error="
                << conditioning.absolute_traction_error_N
                << " wall_cancellation_ratio="
                << conditioning.cancellation_ratio
                << " wall_functional_points=" << conditioning.point_count
                << '\n';
      const auto &traction_decomposition =
          results[level].wall_traction_error_decomposition;
      std::cerr << probe.selector << " cells=" << grids[level]
                << " traction_exact_cell_reconstruction_l2="
                << traction_decomposition.exact_cell_reconstruction_l2
                << " traction_numerical_cell_contamination_l2="
                << traction_decomposition.numerical_cell_contamination_l2
                << " traction_total_error_l2="
                << traction_decomposition.total_error_l2
                << " traction_split_closure_linf="
                << traction_decomposition.pointwise_split_closure_linf
                << " traction_decomposition_points="
                << traction_decomposition.point_count << '\n';
      const auto &extremum = results[level].pressure_extremum_authority;
      std::cerr << probe.selector << " cells=" << grids[level]
                << " pressure_extremum_authority_link="
                << extremum.nearest_link
                << " pressure_extremum_owner=" << extremum.owner_rank
                << " pressure_extremum_direct=" << extremum.direct_fluid_cell
                << " pressure_extremum_link_distance="
                << extremum.cell_to_wall_distance_m
                << " pressure_extremum_condition="
                << extremum.condition_estimate
                << " pressure_extremum_donors=" << extremum.donor_count
                << " pressure_extremum_donor_fingerprint="
                << extremum.donor_fingerprint
                << " pressure_extremum_value_amplification="
                << extremum.pressure_value_amplification
                << " pressure_extremum_gradient_amplification="
                << extremum.normal_gradient_amplification << '\n';
    }
    const std::array errors{formal_errors(results[0]),
                            formal_errors(results[1])};
    for (std::size_t row_index = 0U; row_index < probe.row_count; ++row_index) {
      const auto row = probe.rows[row_index];
      const auto first =
          diagnostic_order(errors[0][row], errors[1][row], grids[0], grids[1]);
      std::cerr << probe.selector << ' ' << kFormalErrorNames[row]
                << " fast_order=";
      if (first.has_value())
        std::cerr << *first;
      else
        std::cerr << "nonmonotone";
      std::cerr << '\n';
      accepted = accepted && first.has_value() && *first >= 1.8;
    }
    if (probe.pressure_linf) {
      constexpr std::size_t pressure_linf_row = 4U;
      const auto endpoint =
          diagnostic_order(errors[0][pressure_linf_row],
                           errors[1][pressure_linf_row], grids[0], grids[1]);
      std::cerr << probe.selector << " pressure_linf fast_order="
                << (endpoint.has_value() ? std::to_string(*endpoint)
                                         : std::string("nonmonotone"))
                << '\n';
      accepted = accepted && positive_decreasing_fast_pair(
                                   {errors[0][pressure_linf_row],
                                    errors[1][pressure_linf_row]});
      if (require_pressure_linf_order)
        accepted = accepted && endpoint.has_value() && *endpoint >= 1.8;
    }
    const auto pressure_measure_error = [&](std::size_t level,
                                            std::size_t candidate) {
      const auto &measure = results[level].pressure_measure;
      const runtime::Real3 value =
          candidate == 0U   ? measure.background_consistency_N
          : candidate == 1U ? measure.projected_consistency_N
          : candidate == 2U ? measure.surface_partition_consistency_N
          : candidate == 3U ? measure.exact_state_a22_consistency_N
                            : measure.numerical_state_a22_defect_N;
      return test::stage3::max_abs(value) /
             results[level].pressure_force_scale_N;
    };
    constexpr std::array<const char *, 5> pressure_measure_names{
        "background", "projected", "surface_partition", "exact_state_a22",
        "numerical_state_a22_defect"};
    for (std::size_t candidate = 0U; candidate < pressure_measure_names.size();
         ++candidate) {
      const auto order = diagnostic_order(
          pressure_measure_error(0U, candidate),
          pressure_measure_error(1U, candidate), grids[0], grids[1]);
      std::cerr << probe.selector << " pressure_measure_"
                << pressure_measure_names[candidate] << " fast_order="
                << (order.has_value() ? std::to_string(*order)
                                      : std::string("nonmonotone"))
                << '\n';
    }
    const auto pressure_scalar_error = [&](std::size_t level,
                                           std::size_t candidate) {
      const auto &scalar = results[level].pressure_scalar;
      const runtime::Real3 value =
          candidate == 0U   ? scalar.numerical_face_consistency_N
          : candidate == 1U ? scalar.exact_face_consistency_N
          : candidate == 2U ? scalar.exact_wall_consistency_N
          : candidate == 3U ? scalar.exact_hybrid_consistency_N
          : candidate == 4U ? scalar.numerical_projected_defect_consistency_N
          : candidate == 5U ? scalar.exact_projected_defect_consistency_N
          : candidate == 6U ? scalar.numerical_wall_projected_consistency_N
                            : scalar.exact_wall_projected_consistency_N;
      return test::stage3::max_abs(value) /
             results[level].pressure_force_scale_N;
    };
    constexpr std::array<const char *, 8> pressure_scalar_names{
        "numerical_face", "exact_face", "exact_wall", "exact_hybrid",
        "numerical_projected_defect", "exact_projected_defect",
        "numerical_wall_projected", "exact_wall_projected"};
    for (std::size_t candidate = 0U; candidate < pressure_scalar_names.size();
         ++candidate) {
      const auto order = diagnostic_order(
          pressure_scalar_error(0U, candidate),
          pressure_scalar_error(1U, candidate), grids[0], grids[1]);
      std::cerr << probe.selector << " pressure_scalar_"
                << pressure_scalar_names[candidate] << " fast_order="
                << (order.has_value() ? std::to_string(*order)
                                      : std::string("nonmonotone"))
                << '\n';
    }
  }
  HUNDUN_CHECK(accepted);
}

void run_exact_momentum_residual_fast(const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(mpi.size() == 1);
  const auto diagnostic_order = [](double coarse, double fine) {
    HUNDUN_CHECK(coarse > 0.0);
    HUNDUN_CHECK(fine > 0.0);
    HUNDUN_CHECK(std::isfinite(coarse));
    HUNDUN_CHECK(std::isfinite(fine));
    return std::log(coarse / fine) / std::log(2.0);
  };
  constexpr std::array<int, 2> grids{12, 24};
  constexpr std::array<const char *, 7> term_names{
      "time", "convection", "viscous_remainder", "pressure",
      "implicit_reference", "source", "total"};
  for (const auto *case_selector : {"sphere_uniform_acceptance",
                                    "prism_warped_acceptance",
                                    "translated_sphere_warped_acceptance"}) {
    std::array<test::stage3::ExactMomentumResidualDiagnostics, 2> results;
    for (std::size_t level = 0U; level < grids.size(); ++level) {
      auto definition = formal_case(case_selector, grids[level]);
      definition.collect_force = false;
      definition.collect_exact_momentum_residual = true;
      const auto run = test::stage3::run_manufactured_case(mpi, definition);
      results[level] = run.exact_momentum_residual;
      HUNDUN_CHECK(results[level].available);
      HUNDUN_CHECK(results[level].interface_row_count > 0U);
      HUNDUN_CHECK(results[level].bulk_row_count > 0U);
      HUNDUN_CHECK(results[level].single_link_row_count > 0U);
      HUNDUN_CHECK(results[level].multi_link_row_count > 0U);
      HUNDUN_CHECK(std::isfinite(
          results[level].single_link_coherent_pressure_balance_rms_N_per_m3));
      HUNDUN_CHECK(
          results[level]
              .single_link_coherent_integral_pressure_balance_rms_N_per_m3 >
          0.0);
      HUNDUN_CHECK(
          results[level]
              .single_link_unconstrained_integral_pressure_balance_rms_N_per_m3 >
          0.0);
      HUNDUN_CHECK(
          std::isfinite(results[level].pointwise_split_closure_linf_N));
      std::cerr << "exact_momentum_residual_fast case=" << case_selector
                << " cells=" << grids[level]
                << " interface_rows="
                << results[level].interface_row_count
                << " bulk_rows=" << results[level].bulk_row_count
                << " split_closure_linf="
                << results[level].pointwise_split_closure_linf_N
                << " interface_pressure_balance="
                << results[level].interface_pressure_balance_rms_N_per_m3
                << " bulk_pressure_balance="
                << results[level].bulk_pressure_balance_rms_N_per_m3
                << " interface_background_pressure_balance="
                << results[level]
                       .interface_background_pressure_balance_rms_N_per_m3
                << " bulk_background_pressure_balance="
                << results[level]
                       .bulk_background_pressure_balance_rms_N_per_m3
                << " interface_reconstructed_face_pressure_balance="
                << results[level]
                       .interface_reconstructed_face_pressure_balance_rms_N_per_m3
                << " bulk_reconstructed_face_pressure_balance="
                << results[level]
                       .bulk_reconstructed_face_pressure_balance_rms_N_per_m3
                << " interface_analytic_face_pressure_balance="
                << results[level]
                       .interface_analytic_face_pressure_balance_rms_N_per_m3
                << " bulk_analytic_face_pressure_balance="
                << results[level]
                       .bulk_analytic_face_pressure_balance_rms_N_per_m3
                << " interface_reconstructed_to_analytic_face="
                << results[level]
                       .interface_reconstructed_to_analytic_face_rms_N_per_m3
                << " interface_unconstrained_center_pressure_balance="
                << results[level]
                       .interface_unconstrained_center_pressure_balance_rms_N_per_m3
                << " bulk_unconstrained_center_pressure_balance="
                << results[level]
                       .bulk_unconstrained_center_pressure_balance_rms_N_per_m3
                << " interface_analytic_center_pressure_balance="
                << results[level]
                       .interface_analytic_center_pressure_balance_rms_N_per_m3
                << " bulk_analytic_center_pressure_balance="
                << results[level]
                       .bulk_analytic_center_pressure_balance_rms_N_per_m3
                << " interface_constrained_to_unconstrained_center="
                << results[level]
                       .interface_constrained_to_unconstrained_center_rms_N_per_m3
                << " interface_exact_shared_center_pressure_balance="
                << results[level]
                       .interface_exact_shared_center_pressure_balance_rms_N_per_m3
                << " bulk_exact_shared_center_pressure_balance="
                << results[level]
                       .bulk_exact_shared_center_pressure_balance_rms_N_per_m3
                << " interface_shared_center_correction="
                << results[level]
                       .interface_shared_center_correction_rms_N_per_m3
                << " shared_center_correction_conservation_linf="
                << results[level]
                       .shared_center_correction_conservation_linf_N
                << " interface_full_exact_center_pressure_balance="
                << results[level]
                       .interface_full_exact_center_pressure_balance_rms_N_per_m3
                << " bulk_full_exact_center_pressure_balance="
                << results[level]
                       .bulk_full_exact_center_pressure_balance_rms_N_per_m3
                << " interface_full_exact_integral_pressure_balance="
                << results[level]
                       .interface_full_exact_integral_pressure_balance_rms_N_per_m3
                << " bulk_full_exact_integral_pressure_balance="
                << results[level]
                       .bulk_full_exact_integral_pressure_balance_rms_N_per_m3
                << " interface_background_to_full_exact_integral="
                << results[level]
                       .interface_background_to_full_exact_integral_rms_N_per_m3
                << " single_link_rows="
                << results[level].single_link_row_count
                << " multi_link_rows="
                << results[level].multi_link_row_count
                << " single_link_current_background_pressure_balance="
                << results[level]
                       .single_link_current_background_pressure_balance_rms_N_per_m3
                << " single_link_full_exact_center_pressure_balance="
                << results[level]
                       .single_link_full_exact_center_pressure_balance_rms_N_per_m3
                << " single_link_coherent_pressure_balance="
                << results[level]
                       .single_link_coherent_pressure_balance_rms_N_per_m3
                << " single_link_coherent_integral_pressure_balance="
                << results[level]
                       .single_link_coherent_integral_pressure_balance_rms_N_per_m3
                << " single_link_unconstrained_integral_pressure_balance="
                << results[level]
                       .single_link_unconstrained_integral_pressure_balance_rms_N_per_m3
                << " single_link_constraint_integral_difference="
                << results[level]
                       .single_link_constraint_integral_difference_rms_N_per_m3
                << " interface_pressure_wall_defect="
                << results[level]
                       .interface_pressure_wall_defect_rms_N_per_m3
                << " bulk_pressure_wall_defect="
                << results[level].bulk_pressure_wall_defect_rms_N_per_m3
                << " interface_pressure_linf="
                << results[level]
                       .interface_pressure_balance_linf_N_per_m3[0]
                << " interface_background_pressure_linf="
                << results[level]
                       .interface_pressure_balance_linf_N_per_m3[1]
                << " interface_full_exact_center_pressure_linf="
                << results[level]
                       .interface_pressure_balance_linf_N_per_m3[2]
                << " interface_wall_defect_linf="
                << results[level]
                       .interface_pressure_balance_linf_N_per_m3[3]
                << " face_area_closure_linf="
                << results[level].face_area_closure_linf_m2
                << " face_correction_closure_linf="
                << results[level].face_correction_closure_linf_N;
      for (std::size_t term = 0U; term < term_names.size(); ++term) {
        HUNDUN_CHECK(std::isfinite(
            results[level].interface_rms_N_per_m3[term]));
        HUNDUN_CHECK(
            std::isfinite(results[level].bulk_rms_N_per_m3[term]));
        std::cerr << " interface_" << term_names[term] << '='
                  << results[level].interface_rms_N_per_m3[term]
                  << " bulk_" << term_names[term] << '='
                  << results[level].bulk_rms_N_per_m3[term];
      }
      std::cerr << '\n';
    }
    for (std::size_t term = 0U; term < term_names.size(); ++term) {
      const double interface_order = diagnostic_order(
          results[0].interface_rms_N_per_m3[term],
          results[1].interface_rms_N_per_m3[term]);
      const double bulk_order = diagnostic_order(
          results[0].bulk_rms_N_per_m3[term],
          results[1].bulk_rms_N_per_m3[term]);
      std::cerr << "exact_momentum_residual_fast case=" << case_selector
                << " term=" << term_names[term]
                << " interface_order=" << interface_order
                << " bulk_order=" << bulk_order << '\n';
    }
    std::cerr << "exact_momentum_residual_fast case=" << case_selector
              << " pressure_balance_interface_order="
              << diagnostic_order(
                     results[0].interface_pressure_balance_rms_N_per_m3,
                     results[1].interface_pressure_balance_rms_N_per_m3)
              << " pressure_balance_bulk_order="
              << diagnostic_order(
                     results[0].bulk_pressure_balance_rms_N_per_m3,
                     results[1].bulk_pressure_balance_rms_N_per_m3)
              << " background_pressure_balance_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_background_pressure_balance_rms_N_per_m3,
                     results[1]
                         .interface_background_pressure_balance_rms_N_per_m3)
              << " background_pressure_balance_bulk_order="
              << diagnostic_order(
                     results[0]
                         .bulk_background_pressure_balance_rms_N_per_m3,
                     results[1]
                         .bulk_background_pressure_balance_rms_N_per_m3)
              << " reconstructed_face_pressure_balance_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_reconstructed_face_pressure_balance_rms_N_per_m3,
                     results[1]
                         .interface_reconstructed_face_pressure_balance_rms_N_per_m3)
              << " reconstructed_face_pressure_balance_bulk_order="
              << diagnostic_order(
                     results[0]
                         .bulk_reconstructed_face_pressure_balance_rms_N_per_m3,
                     results[1]
                         .bulk_reconstructed_face_pressure_balance_rms_N_per_m3)
              << " analytic_face_pressure_balance_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_analytic_face_pressure_balance_rms_N_per_m3,
                     results[1]
                         .interface_analytic_face_pressure_balance_rms_N_per_m3)
              << " analytic_face_pressure_balance_bulk_order="
              << diagnostic_order(
                     results[0]
                         .bulk_analytic_face_pressure_balance_rms_N_per_m3,
                     results[1]
                         .bulk_analytic_face_pressure_balance_rms_N_per_m3)
              << " reconstructed_to_analytic_face_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_reconstructed_to_analytic_face_rms_N_per_m3,
                     results[1]
                         .interface_reconstructed_to_analytic_face_rms_N_per_m3)
              << " unconstrained_center_pressure_balance_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_unconstrained_center_pressure_balance_rms_N_per_m3,
                     results[1]
                         .interface_unconstrained_center_pressure_balance_rms_N_per_m3)
              << " unconstrained_center_pressure_balance_bulk_order="
              << diagnostic_order(
                     results[0]
                         .bulk_unconstrained_center_pressure_balance_rms_N_per_m3,
                     results[1]
                         .bulk_unconstrained_center_pressure_balance_rms_N_per_m3)
              << " analytic_center_pressure_balance_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_analytic_center_pressure_balance_rms_N_per_m3,
                     results[1]
                         .interface_analytic_center_pressure_balance_rms_N_per_m3)
              << " analytic_center_pressure_balance_bulk_order="
              << diagnostic_order(
                     results[0]
                         .bulk_analytic_center_pressure_balance_rms_N_per_m3,
                     results[1]
                         .bulk_analytic_center_pressure_balance_rms_N_per_m3)
              << " constrained_to_unconstrained_center_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_constrained_to_unconstrained_center_rms_N_per_m3,
                     results[1]
                         .interface_constrained_to_unconstrained_center_rms_N_per_m3)
              << " exact_shared_center_pressure_balance_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_exact_shared_center_pressure_balance_rms_N_per_m3,
                     results[1]
                         .interface_exact_shared_center_pressure_balance_rms_N_per_m3)
              << " exact_shared_center_pressure_balance_bulk_order="
              << diagnostic_order(
                     results[0]
                         .bulk_exact_shared_center_pressure_balance_rms_N_per_m3,
                     results[1]
                         .bulk_exact_shared_center_pressure_balance_rms_N_per_m3)
              << " shared_center_correction_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_shared_center_correction_rms_N_per_m3,
                     results[1]
                         .interface_shared_center_correction_rms_N_per_m3)
              << " full_exact_center_pressure_balance_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_full_exact_center_pressure_balance_rms_N_per_m3,
                     results[1]
                         .interface_full_exact_center_pressure_balance_rms_N_per_m3)
              << " full_exact_center_pressure_balance_bulk_order="
              << diagnostic_order(
                     results[0]
                         .bulk_full_exact_center_pressure_balance_rms_N_per_m3,
                     results[1]
                         .bulk_full_exact_center_pressure_balance_rms_N_per_m3)
              << " full_exact_integral_pressure_balance_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_full_exact_integral_pressure_balance_rms_N_per_m3,
                     results[1]
                         .interface_full_exact_integral_pressure_balance_rms_N_per_m3)
              << " full_exact_integral_pressure_balance_bulk_order="
              << diagnostic_order(
                     results[0]
                         .bulk_full_exact_integral_pressure_balance_rms_N_per_m3,
                     results[1]
                         .bulk_full_exact_integral_pressure_balance_rms_N_per_m3)
              << " background_to_full_exact_integral_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_background_to_full_exact_integral_rms_N_per_m3,
                     results[1]
                         .interface_background_to_full_exact_integral_rms_N_per_m3)
              << " single_link_current_background_pressure_balance_order="
              << diagnostic_order(
                     results[0]
                         .single_link_current_background_pressure_balance_rms_N_per_m3,
                     results[1]
                         .single_link_current_background_pressure_balance_rms_N_per_m3)
              << " single_link_full_exact_center_pressure_balance_order="
              << diagnostic_order(
                     results[0]
                         .single_link_full_exact_center_pressure_balance_rms_N_per_m3,
                     results[1]
                         .single_link_full_exact_center_pressure_balance_rms_N_per_m3)
              << " single_link_coherent_pressure_balance_order="
              << diagnostic_order(
                     results[0]
                         .single_link_coherent_pressure_balance_rms_N_per_m3,
                     results[1]
                         .single_link_coherent_pressure_balance_rms_N_per_m3)
              << " single_link_coherent_integral_pressure_balance_order="
              << diagnostic_order(
                     results[0]
                         .single_link_coherent_integral_pressure_balance_rms_N_per_m3,
                     results[1]
                         .single_link_coherent_integral_pressure_balance_rms_N_per_m3)
              << " single_link_unconstrained_integral_pressure_balance_order="
              << diagnostic_order(
                     results[0]
                         .single_link_unconstrained_integral_pressure_balance_rms_N_per_m3,
                     results[1]
                         .single_link_unconstrained_integral_pressure_balance_rms_N_per_m3)
              << " single_link_constraint_integral_difference_order="
              << diagnostic_order(
                     results[0]
                         .single_link_constraint_integral_difference_rms_N_per_m3,
                     results[1]
                         .single_link_constraint_integral_difference_rms_N_per_m3)
              << " pressure_wall_defect_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_pressure_wall_defect_rms_N_per_m3,
                     results[1]
                         .interface_pressure_wall_defect_rms_N_per_m3)
              << " pressure_linf_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_pressure_balance_linf_N_per_m3[0],
                     results[1]
                         .interface_pressure_balance_linf_N_per_m3[0])
              << " background_pressure_linf_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_pressure_balance_linf_N_per_m3[1],
                     results[1]
                         .interface_pressure_balance_linf_N_per_m3[1])
              << " full_exact_center_pressure_linf_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_pressure_balance_linf_N_per_m3[2],
                     results[1]
                         .interface_pressure_balance_linf_N_per_m3[2])
              << " wall_defect_linf_interface_order="
              << diagnostic_order(
                     results[0]
                         .interface_pressure_balance_linf_N_per_m3[3],
                     results[1]
                         .interface_pressure_balance_linf_N_per_m3[3])
              << '\n';
  }
}

ManufacturedRunResult run_sphere_smoke(const runtime::MpiContext &mpi,
                                       int cells, bool collect_force) {
  HUNDUN_CHECK(mpi.size() == 1);
  const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
  const auto started = std::chrono::steady_clock::now();
  auto phase = started;
  const auto elapsed = [&] {
    const auto now = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(now - phase).count();
    phase = now;
    return seconds;
  };
  const auto report_phase = [&](const char *name, double seconds) {
    std::cerr << "sphere_phase cells=" << cells << " name=" << name
              << " seconds=" << seconds << " total="
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - started)
                     .count()
              << '\n';
  };
  const double h = 1.0 / static_cast<double>(cells);
  const double dt = 0.05 * h * h;
  const auto surface_mesh = test::stage3::make_manufactured_surface(body, h);
  test::Stage3TemporaryDirectory directory("task11-laminar-sphere");
  const auto path = directory.path() / "sphere.stl";
  test::write_text(path, test::ascii_stl(surface_mesh.triangles, "sphere"));

  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {cells, cells, cells}, {false, false, false},
      runtime::DecompositionOptions{runtime::Int3{1, 1, 1}});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}});
  auto boundaries =
      boundary::BoundaryRegistry::create(manufactured_config(cells), topology);
  const auto surface =
      immersed::ImmersedSurface::load_collective(path, 1.0, mpi, 0);
  const double surface_seconds = elapsed();
  report_phase("surface", surface_seconds);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  const double domain_seconds = elapsed();
  report_phase("domain", domain_seconds);
  const auto ghost_plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  const double ghost_seconds = elapsed();
  report_phase("ghost", ghost_seconds);
  std::optional<immersed::WallQuadraturePlan> wall_plan;
  if (collect_force)
    wall_plan.emplace(immersed::WallQuadraturePlan::create(
        surface, query, domain, topology, geometry, mpi));
  const double wall_seconds = elapsed();
  report_phase("wall", wall_seconds);
  HUNDUN_CHECK(ghost_plan.maximum_halo_reach() <= 4U);

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(manufactured_cell_field("rho", 1U));
  fields.velocity =
      registry.declare_field(manufactured_cell_field("velocity", 3U));
  fields.mechanical_pressure =
      registry.declare_field(manufactured_cell_field("pi", 1U));
  fields.face_velocity =
      registry.declare_field(manufactured_face_field("face_velocity", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  registry.freeze();

  geometry.require_compatible(topology);
  const auto make_cell_averages = [&](double time_s) {
    std::vector<test::stage3::MmsCellAverage> averages(
        topology.owned_cell_count());
    const auto hardware_threads =
        std::max(1U, std::thread::hardware_concurrency());
    const std::size_t worker_count = std::min<std::size_t>(
        std::min<std::size_t>(
            hardware_threads,
            std::max<std::size_t>(1U,
                                  kParallelWorkerBudget /
                                      static_cast<std::size_t>(mpi.size()))),
        topology.owned_cell_count());
    std::vector<std::exception_ptr> failures(worker_count);
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    std::atomic<bool> cancel{false};
    try {
      for (std::size_t worker = 0U; worker < worker_count; ++worker) {
        workers.emplace_back([&, worker] {
          try {
            const auto begin =
                topology.owned_cell_count() * worker / worker_count;
            const auto end =
                topology.owned_cell_count() * (worker + 1U) / worker_count;
            for (mesh::LocalCellId cell = begin;
                 cell < end && !cancel.load(std::memory_order_relaxed);
                 ++cell) {
              if (domain.region(cell) == immersed::CellRegion::fluid)
                averages[cell] =
                    test::stage3::detail::evaluate_cell_average_validated(
                        topology, geometry, cell, body, time_s);
            }
          } catch (...) {
            failures[worker] = std::current_exception();
            cancel.store(true, std::memory_order_relaxed);
          }
        });
      }
    } catch (...) {
      cancel.store(true, std::memory_order_relaxed);
      for (auto &worker : workers)
        worker.join();
      throw;
    }
    for (auto &worker : workers)
      worker.join();
    for (const auto &failure : failures)
      if (failure)
        std::rethrow_exception(failure);
    return averages;
  };
  const auto next_averages = make_cell_averages(dt);
  const double oracle_next_seconds = elapsed();
  report_phase("oracle_next", oracle_next_seconds);
  auto previous_averages = next_averages;
  auto current_averages = next_averages;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    previous_averages[cell] =
        rescale_manufactured_average(next_averages[cell], dt, -dt);
    current_averages[cell] =
        rescale_manufactured_average(next_averages[cell], dt, 0.0);
  }
  const double oracle_rescale_seconds = elapsed();
  const double oracle_seconds = oracle_next_seconds + oracle_rescale_seconds;
  report_phase("oracle_rescale", oracle_rescale_seconds);

  const auto make_face_histories = [&](double time_s) {
    std::vector<test::stage3::MmsFaceHistory> histories(
        topology.local_face_count());
    const auto hardware_threads =
        std::max(1U, std::thread::hardware_concurrency());
    const std::size_t worker_count = std::min<std::size_t>(
        std::min<std::size_t>(
            hardware_threads,
            std::max<std::size_t>(1U,
                                  kParallelWorkerBudget /
                                      static_cast<std::size_t>(mpi.size()))),
        topology.local_face_count());
    std::vector<std::exception_ptr> failures(worker_count);
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    std::atomic<bool> cancel{false};
    try {
      for (std::size_t worker = 0U; worker < worker_count; ++worker) {
        workers.emplace_back([&, worker] {
          try {
            const auto begin =
                topology.local_face_count() * worker / worker_count;
            const auto end =
                topology.local_face_count() * (worker + 1U) / worker_count;
            for (mesh::LocalFaceId face = begin;
                 face < end && !cancel.load(std::memory_order_relaxed);
                 ++face) {
              const bool owner_active = domain.region(topology.owner(face)) ==
                                        immersed::CellRegion::fluid;
              const auto neighbour = topology.neighbour(face);
              const bool neighbour_active =
                  neighbour.has_value() &&
                  domain.region(*neighbour) == immersed::CellRegion::fluid;
              if (owner_active && (!neighbour.has_value() || neighbour_active))
                histories[face] =
                    test::stage3::detail::evaluate_face_history_validated(
                        topology, geometry, face, body, time_s);
            }
          } catch (...) {
            failures[worker] = std::current_exception();
            cancel.store(true, std::memory_order_relaxed);
          }
        });
      }
    } catch (...) {
      cancel.store(true, std::memory_order_relaxed);
      for (auto &worker : workers)
        worker.join();
      throw;
    }
    for (auto &worker : workers)
      worker.join();
    for (const auto &failure : failures)
      if (failure)
        std::rethrow_exception(failure);
    return histories;
  };
  const auto previous_face_histories = make_face_histories(-dt);
  const auto current_face_histories = make_face_histories(0.0);
  const auto next_face_histories =
      collect_force
          ? std::optional<std::vector<test::stage3::MmsFaceHistory>>{}
          : std::optional<std::vector<test::stage3::MmsFaceHistory>>{
                make_face_histories(dt)};

  const auto make_layer =
      [&](const std::vector<test::stage3::MmsCellAverage> &averages,
          const std::vector<test::stage3::MmsFaceHistory> &face_histories) {
        flow::FlowLayerValues values;
        values.density.resize(topology.owned_cell_count(), 0.0);
        values.velocity.resize(topology.owned_cell_count() * 3U, 0.0);
        values.mechanical_pressure.resize(topology.owned_cell_count(), 0.0);
        values.face_velocity.resize(topology.local_face_count() * 3U, 0.0);
        values.face_mass_flux.resize(topology.local_face_count(), 0.0);
        for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
             ++cell) {
          if (domain.region(cell) != immersed::CellRegion::fluid)
            continue;
          const auto &exact = averages[cell];
          values.density[cell] = test::stage3::kReferenceDensityKgPerM3;
          values.velocity[cell * 3U] = canonical_zero(exact.velocity_m_per_s.x);
          values.velocity[cell * 3U + 1U] =
              canonical_zero(exact.velocity_m_per_s.y);
          values.velocity[cell * 3U + 2U] =
              canonical_zero(exact.velocity_m_per_s.z);
          values.mechanical_pressure[cell] =
              canonical_zero(exact.mechanical_pressure_pa);
        }
        for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
             ++face) {
          const bool owner_active = domain.region(topology.owner(face)) ==
                                    immersed::CellRegion::fluid;
          const auto neighbour = topology.neighbour(face);
          const bool neighbour_active =
              neighbour.has_value() &&
              domain.region(*neighbour) == immersed::CellRegion::fluid;
          if (!owner_active || (neighbour.has_value() && !neighbour_active))
            continue;
          const auto &history = face_histories[face];
          values.face_velocity[face * 3U] =
              canonical_zero(history.velocity_m_per_s.x);
          values.face_velocity[face * 3U + 1U] =
              canonical_zero(history.velocity_m_per_s.y);
          values.face_velocity[face * 3U + 2U] =
              canonical_zero(history.velocity_m_per_s.z);
          values.face_mass_flux[face] =
              canonical_zero(test::stage3::kReferenceDensityKgPerM3 *
                             history.owner_normal_volume_flux_m3_per_s);
        }
        return values;
      };

  auto state = flow::FlowState::create(
      registry, {decomposition.local_extent(), topology.local_face_count()},
      fields, {1U, 0.0, dt, dt, flow::MomentumTimeOrder::bdf2});
  state.seed_accepted_layers(
      make_layer(previous_averages, previous_face_histories),
      make_layer(current_averages, current_face_histories));
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 4));
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  const double state_seconds = elapsed();
  report_phase("state", state_seconds);
  immersed::LocalFlowPatternTransform transform;
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, &domain, &ghost_plan,
      wall_plan ? &*wall_plan : nullptr, &transform, nullptr, mpi, execution,
      halo, momentum_solver, {&mx, &my, &mz}, pressure_solver, pressure_pc);
  const double immersed_flow_seconds = elapsed();
  report_phase("immersed_flow", immersed_flow_seconds);
  std::vector<double> source(domain.active_cells().owned_active_count() * 3U,
                             0.0);
  const auto &active_ids = domain.active_cells().ordered_global_ids();
  for (std::size_t row = 0U; row < domain.active_cells().owned_active_count();
       ++row) {
    const auto local = topology.find_local_cell(active_ids[row]);
    HUNDUN_CHECK(local.has_value());
    const auto &average = next_averages[*local].body_source_N_per_m3;
    source[row * 3U] = average.x;
    source[row * 3U + 1U] = average.y;
    source[row * 3U + 2U] = average.z;
  }
  flow::test::ImmersedFlowTestAccess::set_manufactured_body_source(immersed_flow,
                                                                 source);
  std::optional<flow::test::ImmersedFlowExactMomentumResidualReport>
      exact_momentum_residual;
  if (!collect_force) {
    HUNDUN_CHECK(next_face_histories.has_value());
    std::array<std::vector<flow::test::ImmersedFlowWallGradientSnapshot>, 3>
        exact_wall_gradients;
    const std::array<double, 3> times{-dt, 0.0, dt};
    for (std::size_t time_index = 0U; time_index < times.size(); ++time_index) {
      auto &values = exact_wall_gradients[time_index];
      values.reserve(domain.links().size());
      for (const auto &link : domain.links()) {
        const auto fluid = topology.find_local_cell(link.fluid_cell);
        if (!fluid.has_value() ||
            topology.cell_ownership(*fluid) != mesh::EntityOwnership::owned)
          continue;
        const auto sample =
            test::stage3::evaluate_mms(body, link.wall_intercept_m,
                                       times[time_index]);
        values.push_back(
            {link.id,
             sample.mechanical_pressure_gradient_pa_per_m.x *
                     link.solid_to_fluid_normal.x +
                 sample.mechanical_pressure_gradient_pa_per_m.y *
                     link.solid_to_fluid_normal.y +
                 sample.mechanical_pressure_gradient_pa_per_m.z *
                     link.solid_to_fluid_normal.z});
      }
    }
    exact_momentum_residual =
        flow::test::ImmersedFlowTestAccess::exact_momentum_residual_terms(
            immersed_flow, state,
            make_layer(next_averages, *next_face_histories),
            flow::make_momentum_time_stencil(
                flow::MomentumTimeOrder::bdf2, dt, dt),
            test::stage3::kReferenceDensityKgPerM3,
            test::stage3::kDynamicViscosityPaS, source,
            std::move(exact_wall_gradients));
  }
  const double source_seconds = elapsed();
  report_phase("source", source_seconds);
  std::cerr << "sphere_setup cells=" << cells << " surface=" << surface_seconds
            << " domain=" << domain_seconds << " ghost=" << ghost_seconds
            << " wall=" << wall_seconds << " oracle=" << oracle_seconds
            << " state=" << state_seconds << " immersed_flow=" << immersed_flow_seconds
            << " source=" << source_seconds << '\n';
  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                    test::stage3::kReferenceDensityKgPerM3,
                                    test::stage3::kDynamicViscosityPaS,
                                    std::nullopt,
                                    std::nullopt,
                                    std::nullopt};
  const auto stencil =
      flow::make_momentum_time_stencil(flow::MomentumTimeOrder::bdf2, dt, dt);
  linear::SolveControl solve_control;
  solve_control.atol = test::stage3::kManufacturedSolveAtol;
  solve_control.rtol = test::stage3::kManufacturedSolveRtol;
  solve_control.max_iterations = test::stage3::kManufacturedSolveMaxIterations;
  solve_control.residual_recompute_interval =
      test::stage3::kManufacturedResidualRecomputeInterval;
  const auto attempt =
      immersed_flow.attempt(state, physics, stencil, solve_control, solve_control);
  const double attempt_seconds = elapsed();
  std::cerr << "sphere_phases cells=" << cells << " surface=" << surface_seconds
            << " domain=" << domain_seconds << " ghost=" << ghost_seconds
            << " wall=" << wall_seconds << " oracle=" << oracle_seconds
            << " state=" << state_seconds << " immersed_flow=" << immersed_flow_seconds
            << " source=" << source_seconds << " attempt=" << attempt_seconds
            << " total="
            << std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                             started)
                   .count()
            << '\n';
  flow::test::ImmersedFlowTestAccess::clear_manufactured_body_source(immersed_flow);
  const auto &base = std::get<flow::StepAttemptReport>(attempt.base);
  std::cerr << "sphere_solves cells=" << cells
            << " momentum_iterations=" << base.momentum.components[0].iterations
            << ',' << base.momentum.components[1].iterations << ','
            << base.momentum.components[2].iterations
            << " pressure_iterations=" << base.pressure[0].iterations << ','
            << base.pressure[1].iterations
            << " pressure_matvecs=" << base.pressure[0].matvec_count << ','
            << base.pressure[1].matvec_count
            << " pressure_initial=" << base.pressure[0].initial_residual
            << " pressure_final=" << base.pressure[0].final_residual
            << " pressure_reason=" << static_cast<int>(base.pressure[0].reason)
            << '\n';
  HUNDUN_CHECK(base.disposition == flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(base.reason == flow::StepFailureReason::none);
  HUNDUN_CHECK(base.pressure_corrector_count == 2U);
  HUNDUN_CHECK(attempt.force.has_value() == collect_force);

  const auto committed = state.snapshot(flow::FlowLayer::committed);
  double error_square = 0.0;
  double volume_sum = 0.0;
  double error_linf = 0.0;
  double near_error_square = 0.0;
  double near_volume_sum = 0.0;
  double bulk_error_square = 0.0;
  double bulk_volume_sum = 0.0;
  double exact_linf = 0.0;
  double numerical_linf = 0.0;
  double pressure_error_square = 0.0;
  double pressure_error_linf = 0.0;
  double near_pressure_error_square = 0.0;
  double bulk_pressure_error_square = 0.0;
  double linf_distance = 0.0;
  runtime::Int3 linf_cell{};
  std::uint32_t linf_component = 0U;
  double linf_signed_error = 0.0;
  double linf_exact_value = 0.0;
  double linf_numerical_value = 0.0;
  double linf_pressure_error = 0.0;
  std::size_t linf_incident_wall_links = 0U;
  std::array<double, 3> component_linf{};
  std::array<double, 4> wall_band_linf{};
  std::vector<std::size_t> incident_wall_link_counts(
      topology.owned_cell_count(), 0U);
  for (const auto &link : domain.links()) {
    const auto local = topology.find_local_cell(link.fluid_cell);
    if (local.has_value() && *local < topology.owned_cell_count())
      ++incident_wall_link_counts[*local];
  }
  double pressure_error_integral = 0.0;
  double active_volume = 0.0;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const auto &exact = next_averages[cell];
    const double volume = geometry.cell_volume_m3(cell);
    pressure_error_integral +=
        (committed.mechanical_pressure[cell] - exact.mechanical_pressure_pa) *
        volume;
    active_volume += volume;
  }
  HUNDUN_CHECK(active_volume > 0.0);
  const double pressure_error_mean = pressure_error_integral / active_volume;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const auto &exact = next_averages[cell];
    const std::array<double, 3> expected{exact.velocity_m_per_s.x,
                                         exact.velocity_m_per_s.y,
                                         exact.velocity_m_per_s.z};
    const double volume = geometry.cell_volume_m3(cell);
    const double wall_distance = std::sqrt(
        query.closest_point(geometry.cell_center_m(cell)).squared_distance_m2);
    const double pressure_error = committed.mechanical_pressure[cell] -
                                  exact.mechanical_pressure_pa -
                                  pressure_error_mean;
    const std::size_t wall_band = wall_distance <= h       ? 0U
                                  : wall_distance <= 2.0 * h ? 1U
                                  : wall_distance <= 4.0 * h ? 2U
                                                             : 3U;
    pressure_error_square += pressure_error * pressure_error * volume;
    pressure_error_linf =
        std::max(pressure_error_linf, std::abs(pressure_error));
    if (test::stage3::manufactured_near_wall_band_contains(wall_distance))
      near_pressure_error_square += pressure_error * pressure_error * volume;
    else
      bulk_pressure_error_square += pressure_error * pressure_error * volume;
    for (std::size_t component = 0U; component < 3U; ++component) {
      const double numerical = committed.velocity[cell * 3U + component];
      const double error = numerical - expected[component];
      error_square += error * error * volume;
      if (test::stage3::manufactured_near_wall_band_contains(wall_distance)) {
        near_error_square += error * error * volume;
        near_volume_sum += volume;
      } else {
        bulk_error_square += error * error * volume;
        bulk_volume_sum += volume;
      }
      exact_linf = std::max(exact_linf, std::abs(expected[component]));
      numerical_linf = std::max(numerical_linf, std::abs(numerical));
      component_linf[component] =
          std::max(component_linf[component], std::abs(error));
      wall_band_linf[wall_band] =
          std::max(wall_band_linf[wall_band], std::abs(error));
      if (std::abs(error) > error_linf) {
        error_linf = std::abs(error);
        linf_distance = wall_distance;
        linf_cell = topology.global_cell(cell);
        linf_component = static_cast<std::uint32_t>(component);
        linf_signed_error = error;
        linf_exact_value = expected[component];
        linf_numerical_value = numerical;
        linf_pressure_error = pressure_error;
        linf_incident_wall_links = incident_wall_link_counts[cell];
      }
    }
    volume_sum += 3.0 * volume;
  }
  HUNDUN_CHECK(error_square > 0.0);
  HUNDUN_CHECK(volume_sum > 0.0);
  if (exact_momentum_residual.has_value()) {
    const auto &report = *exact_momentum_residual;
    const auto peak_gid = topology.global_cell_id(linf_cell);
    const auto peak_position =
        std::find(report.active_global_cell_ids.begin(),
                  report.active_global_cell_ids.end(), peak_gid);
    HUNDUN_CHECK(peak_position != report.active_global_cell_ids.end());
    const std::size_t peak_row = static_cast<std::size_t>(
        std::distance(report.active_global_cell_ids.begin(), peak_position));
    const std::size_t peak_offset = peak_row * 3U + linf_component;
    HUNDUN_CHECK(peak_offset < report.total_residual_N.size());
    std::size_t residual_peak_offset = 0U;
    for (std::size_t offset = 1U; offset < report.total_residual_N.size();
         ++offset)
      if (std::abs(report.total_residual_N[offset]) >
          std::abs(report.total_residual_N[residual_peak_offset]))
        residual_peak_offset = offset;
    const std::size_t residual_peak_row = residual_peak_offset / 3U;
    const auto print_exact_residual = [&](const char *kind, std::size_t row,
                                          std::size_t offset) {
      const auto local =
          topology.find_local_cell(report.active_global_cell_ids[row]);
      HUNDUN_CHECK(local.has_value());
      const auto logical = topology.global_cell(*local);
      const double volume = report.cell_volume_m3[row];
      std::cerr << "sphere_exact_residual cells=" << cells << " kind=" << kind
                << " cell=" << logical.x << ',' << logical.y << ','
                << logical.z << " component=" << offset % 3U
                << " interface="
                << static_cast<unsigned>(report.immersed_interface_row[row])
                << " total_N=" << report.total_residual_N[offset]
                << " total_per_volume="
                << report.total_residual_N[offset] / volume
                << " time_N=" << report.time_residual_N[offset]
                << " convective_N=" << report.convective_residual_N[offset]
                << " viscous_remainder_N="
                << report.viscous_remainder_residual_N[offset]
                << " pressure_N=" << report.pressure_residual_N[offset]
                << " background_pressure_N="
                << report.background_pressure_residual_N[offset]
                << " pressure_wall_defect_N="
                << report.pressure_wall_defect_residual_N[offset]
                << " implicit_reference_N="
                << report.implicit_viscous_reference_residual_N[offset]
                << " source_N=" << report.source_residual_N[offset] << '\n';
    };
    print_exact_residual("velocity_peak", peak_row, peak_offset);
    print_exact_residual("residual_peak", residual_peak_row,
                         residual_peak_offset);
    std::vector<const immersed::ImmersedLink *> peak_links;
    for (const auto &link : domain.links())
      if (link.fluid_cell == peak_gid)
        peak_links.push_back(&link);
    HUNDUN_CHECK(peak_links.size() == linf_incident_wall_links);
    HUNDUN_CHECK(!peak_links.empty());
    const immersed::detail::QuadraticFrame constraint_frame{
        geometry.cell_center_m(*topology.find_local_cell(peak_gid)),
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        h};
    std::vector<double> constraint_transpose(
        immersed::detail::kQuadraticBasisSize * peak_links.size(), 0.0);
    for (std::size_t constraint = 0U; constraint < peak_links.size();
         ++constraint) {
      const auto basis =
          immersed::detail::quadratic_directional_derivative_basis(
              peak_links[constraint]->wall_intercept_m,
              peak_links[constraint]->solid_to_fluid_normal,
              constraint_frame);
      for (std::size_t coefficient = 0U; coefficient < basis.size();
           ++coefficient)
        constraint_transpose[coefficient * peak_links.size() + constraint] =
            basis[coefficient];
    }
    const auto constraint_qr = immersed::detail::factorize_design_matrix(
        constraint_transpose, immersed::detail::kQuadraticBasisSize,
        peak_links.size());
    double maximum_normal_alignment = 0.0;
    double minimum_intercept_separation_over_h =
        std::numeric_limits<double>::infinity();
    for (std::size_t first = 0U; first < peak_links.size(); ++first)
      for (std::size_t second = first + 1U; second < peak_links.size();
           ++second) {
        const auto &a = *peak_links[first];
        const auto &b = *peak_links[second];
        maximum_normal_alignment = std::max(
            maximum_normal_alignment,
            std::abs(a.solid_to_fluid_normal.x *
                         b.solid_to_fluid_normal.x +
                     a.solid_to_fluid_normal.y *
                         b.solid_to_fluid_normal.y +
                     a.solid_to_fluid_normal.z *
                         b.solid_to_fluid_normal.z));
        const double dx = a.wall_intercept_m.x - b.wall_intercept_m.x;
        const double dy = a.wall_intercept_m.y - b.wall_intercept_m.y;
        const double dz = a.wall_intercept_m.z - b.wall_intercept_m.z;
        minimum_intercept_separation_over_h =
            std::min(minimum_intercept_separation_over_h,
                     std::sqrt(dx * dx + dy * dy + dz * dz) / h);
      }
    const auto pressure_rows =
        finite_volume::test::ImmersedOperatorTestAccess::
            interface_pressure_rows(
                flow::test::ImmersedFlowTestAccess::immersed_operator(immersed_flow));
    const auto pressure_row =
        std::find_if(pressure_rows.begin(), pressure_rows.end(),
                     [&](const auto &candidate) {
                       return candidate.momentum_cell == peak_gid;
                     });
    HUNDUN_CHECK(pressure_row != pressure_rows.end());
    std::array<double, 3> donor_l1_over_h2{};
    std::array<double, 3> wall_l1_over_h3{};
    for (const auto &term : pressure_row->a22_donor_terms)
      donor_l1_over_h2[term.output_component] +=
          std::abs(term.coefficient) / (h * h);
    for (const auto &term : pressure_row->a22_wall_terms)
      wall_l1_over_h3[term.output_component] +=
          std::abs(term.coefficient) / (h * h * h);
    bool identical_link_donors = true;
    std::size_t donor_union_count = 0U;
    std::vector<mesh::GlobalCellId> donor_union;
    std::vector<std::uint64_t> link_pivots;
    std::vector<double> link_conditions;
    for (std::size_t index = 0U; index < peak_links.size(); ++index) {
      const auto &reconstruction =
          ghost_plan.reconstruction(peak_links[index]->id);
      const auto &donors =
          immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
              reconstruction);
      if (index > 0U) {
        const auto &first =
            immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
                ghost_plan.reconstruction(peak_links.front()->id));
        identical_link_donors = identical_link_donors && donors == first;
      }
      donor_union.insert(donor_union.end(), donors.begin(), donors.end());
      link_pivots.push_back(reconstruction.quality().pivot_fingerprint);
      link_conditions.push_back(reconstruction.quality().condition_estimate);
    }
    std::sort(donor_union.begin(), donor_union.end());
    donor_union.erase(std::unique(donor_union.begin(), donor_union.end()),
                      donor_union.end());
    donor_union_count = donor_union.size();
    std::cerr << "sphere_multilink_condition cells=" << cells << " cell="
              << linf_cell.x << ',' << linf_cell.y << ',' << linf_cell.z
              << " links=" << peak_links.size()
              << " constraint_rank=" << constraint_qr.rank
              << " constraint_condition=" << constraint_qr.condition_estimate
              << " max_normal_alignment=" << maximum_normal_alignment
              << " min_intercept_separation_over_h="
              << minimum_intercept_separation_over_h
              << " a22_donor_l1_over_h2=" << donor_l1_over_h2[0] << ','
              << donor_l1_over_h2[1] << ',' << donor_l1_over_h2[2]
              << " a22_wall_l1_over_h3=" << wall_l1_over_h3[0] << ','
              << wall_l1_over_h3[1] << ',' << wall_l1_over_h3[2]
              << " identical_link_donors=" << identical_link_donors
              << " donor_union_count=" << donor_union_count
              << " link_pivots=" << link_pivots[0] << ',' << link_pivots[1]
              << " link_conditions=" << link_conditions[0] << ','
              << link_conditions[1] << '\n';
    const auto corrections =
        flow::test::ImmersedFlowTestAccess::cell_pressure_correction_authority(
            immersed_flow);
    HUNDUN_CHECK(corrections.correctors.size() == 2U);
    double summed_velocity_change = 0.0;
    for (std::size_t corrector = 0U;
         corrector < corrections.correctors.size(); ++corrector) {
      const auto &record = corrections.correctors[corrector];
      const auto position =
          std::find(record.active_global_cell_ids.begin(),
                    record.active_global_cell_ids.end(), peak_gid);
      HUNDUN_CHECK(position != record.active_global_cell_ids.end());
      const std::size_t row = static_cast<std::size_t>(
          std::distance(record.active_global_cell_ids.begin(), position));
      const std::size_t offset = row * 3U + linf_component;
      HUNDUN_CHECK(offset < record.exact_velocity_change.size());
      summed_velocity_change += record.exact_velocity_change[offset];
      std::cerr << "sphere_pressure_correction cells=" << cells
                << " corrector=" << corrector + 1U << " cell="
                << linf_cell.x << ',' << linf_cell.y << ',' << linf_cell.z
                << " component=" << linf_component
                << " velocity_change="
                << record.exact_velocity_change[offset]
                << " momentum_action="
                << record.momentum_operator_velocity_change[offset]
                << " pressure_action="
                << record.lfp_pressure_residual_change[offset]
                << " pressure_before=" << record.pressure_before_pa[row]
                << " pressure_correction="
                << record.pressure_correction_pa[row]
                << " pressure_after=" << record.pressure_after_pa[row]
                << " pressure_rhs_per_volume="
                << record.pressure_rhs_per_volume[row]
                << " compact_action_per_volume="
                << record.compact_pressure_action_per_volume[row]
                << " interface_action_per_volume="
                << record.interface_pressure_action_per_volume[row]
                << " affine_wall_source_per_volume="
                << record.affine_wall_source_per_volume[row] << '\n';
    }
    std::cerr << "sphere_velocity_split cells=" << cells << " cell="
              << linf_cell.x << ',' << linf_cell.y << ',' << linf_cell.z
              << " component=" << linf_component
              << " predictor="
              << linf_numerical_value - summed_velocity_change
              << " summed_pressure_change=" << summed_velocity_change
              << " final=" << linf_numerical_value
              << " exact=" << linf_exact_value << '\n';
  }
  ManufacturedRunResult result;
  result.velocity_l2 = std::sqrt(error_square / volume_sum);
  result.velocity_linf = error_linf;
  result.near_wall_velocity_l2 =
      std::sqrt(near_error_square / near_volume_sum);
  result.bulk_velocity_l2 = std::sqrt(bulk_error_square / bulk_volume_sum);
  result.pressure_l2 = std::sqrt(pressure_error_square / active_volume);
  result.pressure_linf = pressure_error_linf;
  result.near_wall_pressure_l2 =
      std::sqrt(near_pressure_error_square / near_volume_sum);
  result.bulk_pressure_l2 =
      std::sqrt(bulk_pressure_error_square / bulk_volume_sum);
  result.exact_velocity_linf = exact_linf;
  result.numerical_velocity_linf = numerical_linf;
  result.linf_wall_distance_m = linf_distance;
  result.linf_global_cell = linf_cell;
  result.linf_component = linf_component;
  result.linf_signed_error = linf_signed_error;
  result.linf_exact_value = linf_exact_value;
  result.linf_numerical_value = linf_numerical_value;
  result.linf_pressure_error_pa = linf_pressure_error;
  result.linf_incident_wall_links = linf_incident_wall_links;
  result.component_velocity_linf = component_linf;
  result.wall_band_velocity_linf = wall_band_linf;
  result.force = attempt.force;
  return result;
}

void run(const std::string &selector, const runtime::MpiContext &mpi) {
  const bool formal_acceptance =
      selector == "sphere_uniform_acceptance" ||
      selector == "sphere_warped_acceptance" ||
      selector == "cylinder_uniform_acceptance" ||
      selector == "cylinder_warped_acceptance" ||
      selector == "prism_uniform_acceptance" ||
      selector == "prism_warped_acceptance" ||
      selector == "translated_sphere_uniform_acceptance" ||
      selector == "translated_sphere_warped_acceptance" ||
      selector == "inside_cavity_uniform_acceptance";
  HUNDUN_CHECK(selector == "oracle" || selector == "surface_acceptance" ||
               selector == "signed_force_reference_preflight" ||
               selector == "sphere_smoke" || selector == "sphere_smoke_24" ||
               selector == "sphere_smoke_48" ||
               selector == "inside_smoke_12" || selector == "inside_smoke_24" ||
               selector == "inside_smoke_48" || selector == "sphere_fast" ||
               selector == "sphere_probe_12" ||
               selector == "sphere_probe_24" ||
               selector == "pressure_force_repair_fast" ||
               selector == "functional_selection_fast" ||
               selector == "functional_closure_fast" ||
               selector == "functional_closure_medium" ||
               selector == "pressure_extrema_screen" ||
               selector == "near_wall_pressure_diagnostic" ||
               selector == "exact_momentum_residual_fast" ||
               selector == "exact_force_consistency_fast" ||
               selector == "manufactured_runner_smoke" || formal_acceptance);
  if (selector == "signed_force_reference_preflight") {
    check_signed_force_reference_preflight();
    return;
  }
  HUNDUN_CHECK(order_oracle_is_mutation_sensitive());
  HUNDUN_CHECK(pressure_error_extremum_oracle_is_mutation_sensitive());
  HUNDUN_CHECK(near_wall_pressure_band_oracle_is_mutation_sensitive());
  HUNDUN_CHECK(fixed_physical_near_wall_band_oracle_is_mutation_sensitive());
  HUNDUN_CHECK(pressure_error_moments_oracle_is_mutation_sensitive());
  HUNDUN_CHECK(
      test::stage3::pressure_face_quadrature_oracle_is_mutation_sensitive());
  if (selector == "surface_acceptance") {
    if (mpi.rank() == 0) {
      const auto self = runtime::MpiContext::duplicate(MPI_COMM_SELF);
      test::Stage3TemporaryDirectory temporary("task11-surface-acceptance");
      for (const auto kind :
           {test::stage3::BodyKind::finite_cylinder,
            test::stage3::BodyKind::oblique_rectangular_prism}) {
        for (const int cells : {12, 24, 48}) {
          const auto fixture = test::stage3::make_manufactured_surface(
              test::stage3::approved_body(kind),
              1.0 / static_cast<double>(cells));
          const auto path = temporary.path() /
                            ("body-" + std::to_string(static_cast<int>(kind)) +
                             "-" + std::to_string(cells) + ".stl");
          test::write_text(path, test::ascii_stl(fixture.triangles, "body"));
          const auto surface =
              immersed::ImmersedSurface::load_collective(path, 1.0, self, 0);
          HUNDUN_CHECK(surface.triangle_count() == fixture.triangles.size());
        }
      }
    }
    return;
  }
  if (selector == "manufactured_runner_smoke") {
    const auto result = test::stage3::run_manufactured_case(
        mpi, formal_case("sphere_uniform_acceptance", 12));
    print_formal_result(selector, 12, result);
    HUNDUN_CHECK(result.errors.velocity_l2 > 0.0);
    HUNDUN_CHECK(result.errors.pressure_l2 > 0.0);
    HUNDUN_CHECK(result.errors.pressure_force > 0.0);
    const auto viscous_operator = select_operator_viscous_force(result);
    HUNDUN_CHECK(test::stage3::max_abs(viscous_operator) > 0.0);
    HUNDUN_CHECK(viscous_operator.x == result.operator_force.viscous_N.x);
    HUNDUN_CHECK(viscous_operator.y == result.operator_force.viscous_N.y);
    HUNDUN_CHECK(viscous_operator.z == result.operator_force.viscous_N.z);
    HUNDUN_CHECK(viscous_operator.x == -result.budget_reaction.viscous_N.x);
    HUNDUN_CHECK(viscous_operator.y == -result.budget_reaction.viscous_N.y);
    HUNDUN_CHECK(viscous_operator.z == -result.budget_reaction.viscous_N.z);
    const std::array<runtime::Real3, 3> operator_components{
        result.operator_force.pressure_N, result.operator_force.viscous_N,
        result.operator_force.total_N};
    const std::array<runtime::Real3, 3> budget_components{
        result.budget_reaction.pressure_N, result.budget_reaction.viscous_N,
        result.budget_reaction.total_N};
    const std::array<runtime::Real3, 3> surface_components{
        result.surface_traction.pressure_N, result.surface_traction.viscous_N,
        result.surface_traction.total_N};
    const std::array<runtime::Real3, 3> consistency_components{
        result.consistency.pressure_N, result.consistency.viscous_N,
        result.consistency.total_N};
    for (std::size_t component = 0U; component < operator_components.size();
         ++component) {
      const auto operator_force = operator_components[component];
      const auto budget_reaction = budget_components[component];
      const auto surface_traction = surface_components[component];
      const auto consistency = consistency_components[component];
      HUNDUN_CHECK(operator_force.x == -budget_reaction.x);
      HUNDUN_CHECK(operator_force.y == -budget_reaction.y);
      HUNDUN_CHECK(operator_force.z == -budget_reaction.z);
      HUNDUN_CHECK(consistency.x == operator_force.x - surface_traction.x);
      HUNDUN_CHECK(consistency.y == operator_force.y - surface_traction.y);
      HUNDUN_CHECK(consistency.z == operator_force.z - surface_traction.z);
    }
    return;
  }
  if (selector == "functional_selection_fast") {
    run_functional_selection_fast(mpi);
    return;
  }
  if (selector == "exact_force_consistency_fast") {
    run_exact_force_consistency_fast(mpi);
    return;
  }
  if (selector == "pressure_force_repair_fast") {
    std::array<test::stage3::ManufacturedRunResult, 2> results;
    constexpr std::array<int, 2> grids{24, 48};
    for (std::size_t level = 0U; level < grids.size(); ++level) {
      results[level] = test::stage3::run_manufactured_case(
          mpi, formal_case("sphere_uniform_acceptance", grids[level]));
      print_formal_result("pressure_force_repair_fast", grids[level],
                          results[level]);
    }
    const auto coarse = formal_errors(results[0]);
    const auto fine = formal_errors(results[1]);
    for (const std::size_t row : {8U, 11U}) {
      const double order = observed_order(coarse[row], fine[row]);
      std::cerr << "pressure_force_repair_fast " << kFormalErrorNames[row]
                << " order=" << order << '\n';
      HUNDUN_CHECK(order >= 1.8);
    }
    return;
  }
  if (selector == "functional_closure_fast") {
    run_functional_closure_screen(mpi, {12, 24}, false);
    return;
  }
  if (selector == "functional_closure_medium") {
    run_functional_closure_screen(mpi, {24, 48}, true);
    return;
  }
  if (selector == "pressure_extrema_screen") {
    run_pressure_extrema_screen(mpi);
    return;
  }
  if (selector == "near_wall_pressure_diagnostic") {
    run_near_wall_pressure_diagnostic(mpi);
    return;
  }
  if (selector == "exact_momentum_residual_fast") {
    run_exact_momentum_residual_fast(mpi);
    return;
  }
  if (selector == "inside_smoke_12" || selector == "inside_smoke_24" ||
      selector == "inside_smoke_48") {
    const int cells = selector == "inside_smoke_12"   ? 12
                      : selector == "inside_smoke_24" ? 24
                                                      : 48;
    const auto result = test::stage3::run_manufactured_case(
        mpi, formal_case("inside_cavity_uniform_acceptance", cells));
    print_formal_result(selector, cells, result);
    HUNDUN_CHECK(result.errors.velocity_l2 > 0.0);
    HUNDUN_CHECK(result.errors.pressure_l2 > 0.0);
    HUNDUN_CHECK(result.errors.pressure_force > 0.0);
    return;
  }
  if (formal_acceptance) {
    run_formal_sequence(selector, mpi, {12, 24, 48});
    return;
  }
  if (selector == "sphere_smoke" || selector == "sphere_smoke_24" ||
      selector == "sphere_smoke_48") {
    const int cells = selector == "sphere_smoke"      ? 12
                      : selector == "sphere_smoke_24" ? 24
                                                      : 48;
    const auto result = run_sphere_smoke(mpi, cells, true);
    std::cerr << "sphere_smoke cells=" << cells << " l2=" << result.velocity_l2
              << " linf=" << result.velocity_linf
              << " near_l2=" << result.near_wall_velocity_l2
              << " bulk_l2=" << result.bulk_velocity_l2
              << " pressure_l2=" << result.pressure_l2
              << " pressure_linf=" << result.pressure_linf
              << " near_pressure_l2=" << result.near_wall_pressure_l2
              << " bulk_pressure_l2=" << result.bulk_pressure_l2
              << " exact_linf=" << result.exact_velocity_linf
              << " numerical_linf=" << result.numerical_velocity_linf
              << " linf_distance=" << result.linf_wall_distance_m
              << " linf_cell=" << result.linf_global_cell.x << ','
              << result.linf_global_cell.y << ',' << result.linf_global_cell.z
              << '\n';
    HUNDUN_CHECK(std::isfinite(result.velocity_l2));
    HUNDUN_CHECK(std::isfinite(result.velocity_linf));
    HUNDUN_CHECK(result.velocity_l2 > 0.0);
    HUNDUN_CHECK(result.velocity_linf > 0.0);
    return;
  }
  if (selector == "sphere_fast") {
    const auto coarse = run_sphere_smoke(mpi, 12, false);
    const auto fine = run_sphere_smoke(mpi, 24, false);
    const double l2_order =
        observed_order(coarse.velocity_l2, fine.velocity_l2);
    const double linf_order =
        observed_order(coarse.velocity_linf, fine.velocity_linf);
    std::cerr << "sphere_fast velocity_l2=" << coarse.velocity_l2 << ','
              << fine.velocity_l2 << " order=" << l2_order
              << " velocity_linf=" << coarse.velocity_linf << ','
              << fine.velocity_linf << " order=" << linf_order
              << " near_l2=" << coarse.near_wall_velocity_l2 << ','
              << fine.near_wall_velocity_l2
              << " bulk_l2=" << coarse.bulk_velocity_l2 << ','
              << fine.bulk_velocity_l2 << '\n';
    const auto print_peak = [](const char *resolution,
                               const ManufacturedRunResult &result) {
      std::cerr << "sphere_fast_peak resolution=" << resolution
                << " cell=" << result.linf_global_cell.x << ','
                << result.linf_global_cell.y << ','
                << result.linf_global_cell.z
                << " component=" << result.linf_component
                << " wall_distance=" << result.linf_wall_distance_m
                << " signed_error=" << result.linf_signed_error
                << " exact=" << result.linf_exact_value
                << " numerical=" << result.linf_numerical_value
                << " pressure_error=" << result.linf_pressure_error_pa
                << " incident_wall_links="
                << result.linf_incident_wall_links << " component_linf="
                << result.component_velocity_linf[0] << ','
                << result.component_velocity_linf[1] << ','
                << result.component_velocity_linf[2] << " wall_band_linf="
                << result.wall_band_velocity_linf[0] << ','
                << result.wall_band_velocity_linf[1] << ','
                << result.wall_band_velocity_linf[2] << ','
                << result.wall_band_velocity_linf[3] << '\n';
    };
    print_peak("coarse", coarse);
    print_peak("fine", fine);
    HUNDUN_CHECK(l2_order >= 1.8);
    HUNDUN_CHECK(linf_order >= 1.8);
    return;
  }
  if (selector == "sphere_probe_12" || selector == "sphere_probe_24") {
    const int cells = selector == "sphere_probe_12" ? 12 : 24;
    const auto result = run_sphere_smoke(mpi, cells, false);
    std::cerr << "sphere_probe cells=" << cells
              << " velocity_l2=" << result.velocity_l2
              << " velocity_linf=" << result.velocity_linf
              << " pressure_l2=" << result.pressure_l2
              << " pressure_linf=" << result.pressure_linf << '\n';
    HUNDUN_CHECK(std::isfinite(result.velocity_l2));
    HUNDUN_CHECK(std::isfinite(result.velocity_linf));
    return;
  }
  check_dual3();
  check_mms(test::stage3::BodyKind::sphere);
  check_mms(test::stage3::BodyKind::finite_cylinder);
  check_mms(test::stage3::BodyKind::oblique_rectangular_prism);
  check_exact_wall_zeros();
  check_outer_pressure_neumann_compatibility();
  check_surface_and_force_oracle(test::stage3::BodyKind::finite_cylinder);
  check_surface_and_force_oracle(
      test::stage3::BodyKind::oblique_rectangular_prism);
  check_surface_and_force_oracle(test::stage3::BodyKind::sphere);
  check_surface_and_force_oracle(test::stage3::BodyKind::inside_sphere_cavity);
  check_source_quadrature(mpi);
  check_face_history_quadrature(mpi);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  const auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  const std::string selector = argc == 2 ? argv[1] : "";
  return hundun::test::run([&] { run(selector, mpi); });
}
