// SPDX-License-Identifier: Apache-2.0

#include "../support/piso_fixture.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
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

bool accept(TimeControllerState& state, const StepTime& step) {
  StepTime ignored;
  return static_cast<bool>(state.finish(MPI_COMM_SELF, step, Status{}, ignored));
}

bool time_coefficients(double requested_dt, BdfCoefficients& bdf,
                       double& current_dt, double& previous_dt) {
  constexpr double ratio = 1.25;
  previous_dt = requested_dt / ratio;
  TimeControlSpec spec;
  spec.control = TimeControlKind::adaptive_flow;
  spec.scheme = TimeScheme::variable_bdf2;
  spec.initial_dt = previous_dt;
  spec.minimum_dt = requested_dt * 1.0e-4;
  spec.maximum_dt = requested_dt * 10.0;
  spec.convective_cfl = 1.0;
  spec.viscous_cfl = 1.0;
  spec.thermal_cfl = 1.0;
  spec.species_cfl = 1.0;
  spec.acoustic_cfl = 1.0;
  spec.maximum_growth = ratio;
  spec.retry_factor = 0.5;
  spec.maximum_retries = 2U;
  spec.minimum_bdf_ratio = 0.2;
  spec.maximum_bdf_ratio = 5.0;

  TimeSchemePlan plan;
  TimeControllerState state;
  const double start_time = -(requested_dt + previous_dt);
  if (!TimeSchemePlan::compile(spec, plan) ||
      !TimeControllerState::start(plan, start_time, state)) {
    return false;
  }
  const double infinity = std::numeric_limits<double>::infinity();
  const LocalTimeLimits loose{infinity, infinity, infinity, infinity,
                              infinity};
  StepTime startup;
  if (!state.propose(MPI_COMM_SELF, loose, startup) ||
      startup.bdf.order != 1U || !accept(state, startup)) {
    return false;
  }
  StepTime step;
  if (!state.propose(MPI_COMM_SELF, loose, step) || step.bdf.order != 2U) {
    return false;
  }
  current_dt = step.dt;
  previous_dt = startup.dt;
  bdf = step.bdf;
  return true;
}

struct ErrorPair {
  double production{};
  double missing_a2{};
};

struct CommittedFluxHistory {
  FieldId dependency{};
  StateLayers layers;
  AttemptTransaction transaction;
  FaceFluxStorage storage;
  FinalFaceFluxAuthority authority;
  FinalFaceFluxWriter writer;
};

bool initialize_flux_history(Int3 cells, CommittedFluxHistory& history) {
  FieldRegistry registry;
  FieldSchema schema;
  if (!registry.declare_field("piso.temporal_flux_history", 1U, 0U,
                              history.dependency) ||
      !registry.freeze(schema)) {
    return false;
  }
  const std::array requests{ArenaFieldRequest{
      history.dependency, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout layout;
  return static_cast<bool>(ArenaLayout::compile(
             schema, {requests.data(), requests.size()}, layout)) &&
         static_cast<bool>(StateLayers::allocate(layout, history.layers)) &&
         static_cast<bool>(AttemptTransaction::create(
             history.layers.field_count(), 1U, history.layers.field_count(),
             history.transaction)) &&
         static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                           history.storage)) &&
         static_cast<bool>(history.authority.claim(
             71U, 0U, history.transaction, history.writer));
}

bool commit_uniform_x_flux(const CartesianKernelPlan& kernels, Int3 cells,
                           double mass_flux, CommittedFluxHistory& history,
                           ConstFaceFluxView& committed) {
  if (!history.transaction.begin(history.layers) ||
      !history.transaction.revise_trial(history.dependency)) {
    return false;
  }
  const RevisionDependency dependency{
      AttemptTransaction::field_revision_source(history.dependency),
      history.transaction.trial_revision(history.dependency)};
  PendingFaceFluxView pending;
  if (!history.writer.begin_pending(history.transaction, history.storage,
                                    pending)) {
    return false;
  }
  OwnedField density = make_field(80U, cells, 1U, 2U, 7101U, 8101U);
  OwnedField velocity = make_field(81U, cells, 3U, 2U, 7102U, 8102U);
  fill(density, 1.0);
  fill(velocity, 0.0);
  const double face_velocity =
      mass_flux * static_cast<double>(cells.y * cells.z);
  for (std::int32_t z = -2; z < cells.z + 2; ++z) {
    for (std::int32_t y = -2; y < cells.y + 2; ++y) {
      for (std::int32_t x = -2; x < cells.x + 2; ++x) {
        velocity.view.unchecked({x, y, z}, 0U) = face_velocity;
      }
    }
  }
  const std::array<ConstFieldView, 2U> reads{
      as_const(density.view), as_const(velocity.view)};
  const KernelInvocation invocation{
      {reads.data(), reads.size()}, {}, {{0, 0, 0}, cells},
      0U, 0U, 1U, 0U, nullptr};
  const std::array dependencies{dependency};
  return static_cast<bool>(reconstruct_mass_flux(kernels, invocation,
                                                 pending)) &&
         static_cast<bool>(history.writer.publish_pending(
             {dependencies.data(), dependencies.size()}, pending)) &&
         static_cast<bool>(history.transaction.collective_finish(
             MPI_COMM_SELF, Status{})) &&
         static_cast<bool>(history.writer.committed(history.storage,
                                                     committed));
}

void fill_face_flux(FaceFluxView flux, double value) {
  const std::array<FaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (FaceFieldView face : faces) {
    for (std::int32_t z = 0; z < face.extents.z; ++z) {
      for (std::int32_t y = 0; y < face.extents.y; ++y) {
        for (std::int32_t x = 0; x < face.extents.x; ++x) {
          face.unchecked({x, y, z}) = value;
        }
      }
    }
  }
}

void bind_committed_flux_history(PisoIntermediateInput& input,
                                 ConstFaceFluxView accepted,
                                 ConstFaceFluxView previous = {}) {
  input.predictor.accepted_face_flux = accepted.revision;
  input.predictor.previous_face_flux = previous.revision;
  input.predictor.committed_face_flux_authority =
      accepted.certificate.authority();
  input.predictor.committed_face_flux_storage =
      accepted.certificate.storage();
  input.predictor.committed_face_flux_revision_domain =
      accepted.certificate.revision_domain();
  input.committed_face_history = {accepted, previous};
}

bool evaluate_error(PeriodicPisoFixture& fixture, double requested_dt,
                    std::uint64_t revision_seed,
                    ConstFaceFluxView accepted_flux,
                    ConstFaceFluxView previous_flux, ErrorPair& error) {
  BdfCoefficients bdf;
  double dt = 0.0;
  double previous_dt = 0.0;
  if (!time_coefficients(requested_dt, bdf, dt, previous_dt)) {
    return false;
  }

  const auto density = [](double time) { return 2.0 + std::exp(time); };
  fill(fixture.density, density(0.0));
  ++fixture.density.view.revision;
  PisoIntermediateInput intermediate_input =
      fixture.intermediate_input(bdf, 3000U + revision_seed);
  bind_committed_flux_history(intermediate_input, accepted_flux,
                              previous_flux);
  PisoIntermediateCertificate intermediate;
  if (!fixture.coupler.refresh(intermediate_input, intermediate)) {
    return false;
  }

  OwnedField accepted = make_field(50U, fixture.patch.cells, 1U, 0U,
                                   3100U + revision_seed,
                                   4100U + revision_seed);
  OwnedField previous = make_field(51U, fixture.patch.cells, 1U, 0U,
                                   3200U + revision_seed,
                                   4200U + revision_seed);
  OwnedField drho_dp = make_field(52U, fixture.patch.cells, 1U, 0U,
                                  3300U + revision_seed,
                                  4300U + revision_seed);
  OwnedField diagonal = make_field(53U, fixture.patch.cells, 1U, 0U,
                                   3400U + revision_seed,
                                   4400U + revision_seed);
  OwnedField rhs = make_field(54U, fixture.patch.cells, 1U, 0U,
                              3500U + revision_seed,
                              4500U + revision_seed);
  fill(accepted, density(-dt));
  fill(previous, density(-(dt + previous_dt)));
  fill(drho_dp, 0.01);

  PressureCorrectionInput pressure_input;
  pressure_input.intermediate = intermediate;
  pressure_input.pressure_reference = intermediate_input.pressure_reference;
  pressure_input.density_trial = as_const(fixture.density.view);
  pressure_input.density_accepted = as_const(accepted.view);
  pressure_input.density_previous = as_const(previous.view);
  pressure_input.drho_dp_h_y = as_const(drho_dp.view);
  pressure_input.bdf = bdf;
  pressure_input.time = intermediate_input.momentum.time;
  pressure_input.geometry = intermediate_input.momentum.geometry;
  pressure_input.numeric_boundary = intermediate_input.numeric_boundary;
  const PressureCorrectionSystemView system{diagonal.view, rhs.view};
  PressureCorrectionCertificate certificate;
  if (!fixture.coupler.assemble_pressure_system(pressure_input, system,
                                                certificate)) {
    return false;
  }
  const double volume = 1.0 /
      static_cast<double>(fixture.patch.cells.x * fixture.patch.cells.y *
                          fixture.patch.cells.z);
  const double approximation = -rhs.view.unchecked({0, 0, 0}, 0U) / volume;
  error.production = std::abs(approximation - 1.0);

  // This is a valid but deliberately a2-suppressed mutation: retain a tiny
  // positive a2 so it still enters the production BDF2 path, and fold the
  // removed history weight into a1 to preserve the constant-state identity.
  BdfCoefficients missing = bdf;
  missing.a2 = std::numeric_limits<double>::min();
  missing.a1 = -(missing.a0 + missing.a2);
  pressure_input.bdf = missing;
  ++pressure_input.time;
  pressure_input.pressure_reference.time = pressure_input.time;
  intermediate_input = fixture.intermediate_input(missing, pressure_input.time);
  bind_committed_flux_history(intermediate_input, accepted_flux,
                              previous_flux);
  if (!fixture.coupler.refresh(intermediate_input, intermediate)) {
    return false;
  }
  pressure_input.intermediate = intermediate;
  pressure_input.pressure_reference = intermediate_input.pressure_reference;
  certificate = {};
  if (!fixture.coupler.assemble_pressure_system(pressure_input, system,
                                                certificate)) {
    return false;
  }
  const double mutated = -rhs.view.unchecked({0, 0, 0}, 0U) / volume;
  error.missing_a2 = std::abs(mutated - 1.0);
  return std::isfinite(error.production) && std::isfinite(error.missing_a2);
}

bool test_variable_step_bdf2_order() {
  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(4),
                       "production periodic PISO fixture compiles");
  if (!passed) {
    return false;
  }
  CommittedFluxHistory history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(
      initialize_flux_history(fixture.patch.cells, history) &&
          commit_uniform_x_flux(fixture.equations.kernels(),
                                fixture.patch.cells, 0.0, history,
                                previous_flux) &&
          commit_uniform_x_flux(fixture.equations.kernels(),
                                fixture.patch.cells, 0.0, history,
                                accepted_flux),
      "temporal-order committed zero-flux history publishes");
  if (!passed) {
    return false;
  }
  const std::array<double, 3U> steps{0.08, 0.04, 0.02};
  std::array<ErrorPair, steps.size()> errors{};
  for (std::size_t level = 0U; level < steps.size(); ++level) {
    passed &= expect(evaluate_error(fixture, steps[level], level,
                                    accepted_flux, previous_flux,
                                    errors[level]),
                     "production pressure assembly evaluates temporal defect");
  }
  if (!passed) {
    return false;
  }
  const double order_0 = std::log(errors[0U].production /
                                  errors[1U].production) /
                         std::log(2.0);
  const double order_1 = std::log(errors[1U].production /
                                  errors[2U].production) /
                         std::log(2.0);
  passed &= expect(order_0 >= 1.8 && order_1 >= 1.8,
                   "variable-step BDF2 pressure defect converges at order >= 1.8");
  passed &= expect(errors[2U].missing_a2 > 100.0 * errors[2U].production &&
                       errors[1U].missing_a2 >=
                           0.5 * errors[0U].missing_a2,
                   "omitting the previous-history coefficient is mutation RED");
  if (!passed) {
    std::cerr << "errors=" << errors[0U].production << ','
              << errors[1U].production << ',' << errors[2U].production
              << " orders=" << order_0 << ',' << order_1
              << " missing_a2=" << errors[0U].missing_a2 << ','
              << errors[1U].missing_a2 << ',' << errors[2U].missing_a2
              << '\n';
  }
  return passed;
}

bool test_c1_uses_normalized_bdf_face_history_not_paired_ex2() {
  PeriodicPisoFixture fixture;
  bool passed = expect(fixture.initialize(2),
                       "normalized-history periodic PISO fixture compiles");
  if (!passed) {
    return false;
  }

  CommittedFluxHistory history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(
      initialize_flux_history(fixture.patch.cells, history) &&
          commit_uniform_x_flux(fixture.equations.kernels(),
                                fixture.patch.cells, -1.0, history,
                                previous_flux) &&
          commit_uniform_x_flux(fixture.equations.kernels(),
                                fixture.patch.cells, 2.0, history,
                                accepted_flux),
      "two certified committed face-flux history layers publish");
  if (!passed) {
    return false;
  }

  struct Case {
    BdfCoefficients bdf;
    double paired;
    double normalized;
  };
  const std::array<Case, 3U> cases{{
      {{4.0 / 3.0, -3.0 / 2.0, 1.0 / 6.0, 2U}, 7.0 / 2.0,
       19.0 / 8.0},
      {{3.0 / 2.0, -2.0, 1.0 / 2.0, 2U}, 5.0, 3.0},
      {{5.0 / 3.0, -3.0, 4.0 / 3.0, 2U}, 8.0, 22.0 / 5.0},
  }};
  const double volume = 1.0 / 8.0;
  for (std::size_t index = 0U; index < cases.size(); ++index) {
    const Case& selected = cases[index];
    fill(fixture.velocity, 0.0);
    fill(fixture.momentum_diagonal, selected.bdf.a0 * volume);
    fill(fixture.momentum_rhs, 0.0);
    fill_face_flux(fixture.phi_h_by_a, 0.0);
    fill_face_flux(fixture.trial_flux, selected.paired);

    PisoIntermediateInput input = fixture.intermediate_input(
        selected.bdf, 7201U + static_cast<RevisionToken>(index));
    bind_committed_flux_history(input, accepted_flux, previous_flux);
    PisoIntermediateCertificate certificate;
    passed &= expect(static_cast<bool>(fixture.coupler.refresh(
                         input, certificate)) &&
                         certificate.committed_face_history != 0U,
                     "C1 accepts exact committed BDF face-history authority");
    const double observed = fixture.phi_h_by_a.x.unchecked({1, 0, 0});
    passed &= expect(
        std::abs(observed - selected.normalized) < 1.0e-12 &&
            std::abs(selected.normalized - selected.paired) > 1.0e-3,
        "C1 internal face uses normalized BDF history, not paired EX2");
  }

  const Case& independent_authority = cases[1U];
  fill(fixture.velocity, 0.0);
  fill(fixture.momentum_diagonal,
       independent_authority.bdf.a0 * volume);
  fill(fixture.momentum_rhs, 0.0);
  fill_face_flux(fixture.phi_h_by_a, 0.0);
  fill_face_flux(fixture.trial_flux,
                 std::numeric_limits<double>::quiet_NaN());
  PisoIntermediateInput independent_input = fixture.intermediate_input(
      independent_authority.bdf, 7211U);
  bind_committed_flux_history(independent_input, accepted_flux,
                              previous_flux);
  PisoIntermediateCertificate independent_certificate;
  passed &= expect(
      static_cast<bool>(fixture.coupler.refresh(
          independent_input, independent_certificate)) &&
          independent_certificate.valid() &&
          std::abs(fixture.phi_h_by_a.x.unchecked({1, 0, 0}) -
                   independent_authority.normalized) < 1.0e-12,
      "periodic C1 does not consume attempt-local paired-flux payload");
  return passed;
}

bool test_c1_backward_euler_uses_only_accepted_face_history() {
  PeriodicPisoFixture fixture;
  bool passed = expect(
      fixture.initialize(2),
      "backward-Euler committed-history PISO fixture compiles");
  if (!passed) {
    return false;
  }
  CommittedFluxHistory history;
  ConstFaceFluxView accepted_flux;
  passed &= expect(
      initialize_flux_history(fixture.patch.cells, history) &&
          commit_uniform_x_flux(fixture.equations.kernels(),
                                fixture.patch.cells, 2.0, history,
                                accepted_flux),
      "backward-Euler certified accepted face flux publishes");
  if (!passed) {
    return false;
  }

  const BdfCoefficients bdf{1.0, -1.0, 0.0, 1U};
  fill(fixture.velocity, 0.0);
  fill(fixture.momentum_diagonal, 1.0 / 8.0);
  fill(fixture.momentum_rhs, 0.0);
  fill_face_flux(fixture.phi_h_by_a, 0.0);
  fill_face_flux(fixture.trial_flux, 9.0);
  PisoIntermediateInput input = fixture.intermediate_input(bdf, 7301U);
  bind_committed_flux_history(input, accepted_flux);
  PisoIntermediateCertificate first;
  passed &= expect(
      static_cast<bool>(fixture.coupler.refresh(input, first)) &&
          std::abs(fixture.phi_h_by_a.x.unchecked({1, 0, 0}) - 2.0) <
              1.0e-12,
      "backward Euler uses accepted history and ignores paired EX flux");

  input.corrector = 2U;
  input.prior_corrector = first.dependency;
  input.temporal_reference = {};
  ConstFaceFluxView corrected_trial = as_const(fixture.trial_flux);
  ++corrected_trial.revision;
  input.trial_flux = corrected_trial;
  PisoIntermediateCertificate rejected = first;
  passed &= expect(
      fixture.coupler.refresh(input, rejected).code ==
              StatusCode::invalid_plan &&
          rejected.committed_face_history == first.committed_face_history,
      "corrector two rejects retained committed face-history input");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  const bool passed = test_variable_step_bdf2_order() &&
                      test_c1_uses_normalized_bdf_face_history_not_paired_ex2() &&
                      test_c1_backward_euler_uses_only_accepted_face_history();
  MPI_Finalize();
  return passed ? 0 : 1;
}
