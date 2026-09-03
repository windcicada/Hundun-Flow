// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_flow.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr FieldId kPressure = 2U;
constexpr FieldId kCompressibility = 6U;
constexpr RevisionToken kTime = 501U;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <=
             1.0e-13 * std::max(1.0, std::abs(expected));
}

std::uint64_t double_bits(double value) {
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool same_collective_status(MPI_Comm world, Status status) {
  const std::array<std::uint64_t, 2U> local{
      static_cast<std::uint64_t>(status.code), status.detail};
  std::array<std::uint64_t, 2U> minimum{};
  std::array<std::uint64_t, 2U> maximum{};
  return MPI_Allreduce(local.data(), minimum.data(), 2, MPI_UINT64_T, MPI_MIN,
                       world) == MPI_SUCCESS &&
         MPI_Allreduce(local.data(), maximum.data(), 2, MPI_UINT64_T, MPI_MAX,
                       world) == MPI_SUCCESS &&
         minimum == maximum;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {8, 4, 4};
  mesh.minimum_spacing = {0.125, 0.25, 0.25};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = 8U * 4U * 4U;
  mesh.limits.max_memory_bytes_per_rank = 1U << 28U;
  return mesh;
}

ValidatedModel model_for(const CartesianMeshSpec& mesh) {
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x14c10501U;
  model.pressure_reference = PressureReferenceKind::closed_mass;
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::central2;
  model.schemes.species = ConvectionScheme::central2;
  model.schemes.passive_scalar = ConvectionScheme::central2;
  model.schemes.diffusion = DiffusionScheme::central2;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::no_slip_wall;
    face.thermal_kind = BoundaryKind::adiabatic_wall;
    face.mach_limit = 0.95;
  }
  return model;
}

ThermophysicalSpec thermo_spec() {
  ThermophysicalSpec spec;
  spec.data_file = "analytic.d";
  spec.minimum_temperature = 200.0;
  spec.maximum_temperature = 2000.0;
  spec.temperature_relative_tolerance = 1.0e-12;
  spec.maximum_temperature_iterations = 64U;
  spec.closed_mass_relative_tolerance = 1.0e-12;
  spec.maximum_closed_mass_iterations = 32U;
  spec.maximum_closed_mass_relative_step = 0.2;
  SpeciesThermophysicalSpec air;
  air.stable_name = "air";
  air.molecular_weight = 28.96546;
  air.temperature_switch = 1000.0;
  air.nasa7_low[0U] = 3.5;
  air.nasa7_high[0U] = 3.5;
  air.viscosity_reference = 1.8e-5;
  air.conductivity = 0.026;
  spec.species.push_back(air);
  return spec;
}

struct Fixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  ContributionRegistry contributions;
  EquationPlanSet equations;
  ReductionEngine reductions;
};

bool make_fixture(MPI_Comm world, Fixture& out) {
  const CartesianMeshSpec mesh = mesh_spec();
  const ValidatedModel model = model_for(mesh);
  FieldRegistry registry;
  FieldId density = 0U;
  FieldId velocity = 0U;
  FieldId pressure = 0U;
  FieldId enthalpy = 0U;
  FieldId temperature = 0U;
  if (!registry.require_field("rho", 1U, 2U, density) || density != 0U ||
      !registry.require_field("U", 3U, 2U, velocity) || velocity != 1U ||
      !registry.require_field("pi", 1U, 2U, pressure) ||
      pressure != kPressure ||
      !registry.require_field("h", 1U, 2U, enthalpy) || enthalpy != 3U ||
      !registry.require_field("T", 1U, 2U, temperature) ||
      temperature != 4U ||
      !CartesianGeometryCompiler::compile(world, mesh, {}, out.geometry,
                                          out.patch) ||
      !BoundaryCompiler::compile(world, model, out.geometry, out.patch,
                                 registry, out.boundary, out.schemes,
                                 out.time)) {
    return false;
  }
  const ThermophysicalSpec thermo = thermo_spec();
  if (!ThermodynamicsPlan::compile(thermo, {}, out.thermodynamics) ||
      !TransportPlan::compile(thermo, out.thermodynamics, out.transport)) {
    return false;
  }
  const std::array<FieldId, 8U> fields{0U, 1U, 2U, 3U,
                                      4U, 5U, 6U, 7U};
  if (!out.contributions.configure({fields.data(), fields.size()}) ||
      !out.contributions.freeze()) {
    return false;
  }
  EquationPlanSpec plan;
  plan.density = density;
  plan.velocity = velocity;
  plan.pressure_perturbation = pressure;
  plan.enthalpy = enthalpy;
  plan.temperature = temperature;
  plan.effective_viscosity = 5U;
  plan.velocity_gradient = 7U;
  plan.pressure_compressibility = kCompressibility;
  plan.pressure_reference = PressureReferenceKind::closed_mass;
  plan.closed_mass_service_stage = 37U;
  plan.maximum_cells_per_rank = 8U * 4U * 4U;
  return EquationPlanSet::compile(
             world, out.schemes, out.geometry, out.patch, out.boundary,
             out.contributions, out.thermodynamics, out.transport, plan,
             out.equations) &&
         ReductionEngine::compile(world, ReductionMode::mpi_allreduce, 2U,
                                  out.reductions);
}

struct OwnedField {
  std::vector<double> values;
  FieldView view{};
};

OwnedField make_field(FieldId field, Int3 cells, RevisionToken revision,
                      StorageIdentity storage,
                      RevisionDomainIdentity revision_domain) {
  OwnedField owned;
  const std::size_t count = static_cast<std::size_t>(cells.x) * cells.y *
                            static_cast<std::size_t>(cells.z);
  owned.values.assign(count, 0.0);
  owned.view.base = owned.values.data();
  owned.view.interior = cells;
  owned.view.components = 1U;
  owned.view.stride_y = static_cast<std::size_t>(cells.x);
  owned.view.stride_z = static_cast<std::size_t>(cells.x) *
                        static_cast<std::size_t>(cells.y);
  owned.view.component_stride = count;
  owned.view.field = field;
  owned.view.revision = revision;
  owned.view.storage_identity = storage;
  owned.view.revision_domain = revision_domain;
  return owned;
}

PressureReferenceCertificate predecessor(const Fixture& fixture) {
  const PressureReferencePlan& plan = fixture.equations.pressure_reference();
  return {plan.fingerprint(),
          fixture.equations.thermophysical_predictor().fingerprint(),
          fixture.thermodynamics.fingerprint(),
          601U,
          kTime,
          602U,
          PressureReferenceKind::closed_mass};
}

bool same_certificate(const PressureReferenceCertificate& left,
                      const PressureReferenceCertificate& right) {
  return left.plan == right.plan && left.predictor == right.predictor &&
         left.thermodynamics == right.thermodynamics &&
         left.closure == right.closure && left.time == right.time &&
         left.pressure_reference == right.pressure_reference &&
         left.kind == right.kind;
}

bool same_gauge_certificate(const ClosedGaugeCorrectionCertificate& left,
                            const ClosedGaugeCorrectionCertificate& right) {
  return double_bits(left.shift) == double_bits(right.shift) &&
         double_bits(left.next_pressure_reference) ==
             double_bits(right.next_pressure_reference) &&
         double_bits(left.local_moment) == double_bits(right.local_moment) &&
         double_bits(left.local_weight) == double_bits(right.local_weight) &&
         double_bits(left.global_moment) ==
             double_bits(right.global_moment) &&
         double_bits(left.global_weight) == double_bits(right.global_weight) &&
         double_bits(left.local_post_shift_moment) ==
             double_bits(right.local_post_shift_moment) &&
         double_bits(left.local_post_shift_absolute_moment) ==
             double_bits(right.local_post_shift_absolute_moment) &&
         double_bits(left.global_post_shift_moment) ==
             double_bits(right.global_post_shift_moment) &&
         double_bits(left.global_post_shift_absolute_moment) ==
             double_bits(right.global_post_shift_absolute_moment) &&
         double_bits(left.post_shift_gauge_residual) ==
             double_bits(right.post_shift_gauge_residual) &&
         double_bits(left.post_shift_gauge_tolerance) ==
             double_bits(right.post_shift_gauge_tolerance) &&
         same_certificate(left.output_pressure_reference,
                          right.output_pressure_reference) &&
         left.predecessor_pressure_reference ==
             right.predecessor_pressure_reference &&
         left.time == right.time && left.geometry == right.geometry &&
         left.pressure_correction_authority ==
             right.pressure_correction_authority &&
         left.target_thermodynamic_closure ==
             right.target_thermodynamic_closure &&
         left.activity_local_fingerprint ==
             right.activity_local_fingerprint &&
         left.activity_collective_fingerprint ==
             right.activity_collective_fingerprint &&
         left.collective_transaction == right.collective_transaction &&
         left.rank_local_transaction == right.rank_local_transaction &&
         left.local_active_cells == right.local_active_cells &&
         left.corrector == right.corrector;
}

ClosedGaugeCorrectionCertificate marker_certificate() {
  ClosedGaugeCorrectionCertificate marker;
  marker.shift = -11.0;
  marker.next_pressure_reference = 12.0;
  marker.local_moment = 13.0;
  marker.local_weight = 14.0;
  marker.global_moment = 15.0;
  marker.global_weight = 16.0;
  marker.local_post_shift_moment = 16.5;
  marker.local_post_shift_absolute_moment = 17.0;
  marker.global_post_shift_moment = 18.0;
  marker.global_post_shift_absolute_moment = 19.0;
  marker.post_shift_gauge_residual = 20.0;
  marker.post_shift_gauge_tolerance = 21.0;
  marker.output_pressure_reference = {22U, 23U, 24U, 25U, 26U, 27U,
                                      PressureReferenceKind::closed_mass};
  marker.predecessor_pressure_reference = 28U;
  marker.time = 29U;
  marker.geometry = 30U;
  marker.pressure_correction_authority = 31U;
  marker.target_thermodynamic_closure = 32U;
  marker.activity_local_fingerprint = 33U;
  marker.activity_collective_fingerprint = 34U;
  marker.collective_transaction = 35U;
  marker.rank_local_transaction = 36U;
  marker.local_active_cells = 37U;
  marker.corrector = 2U;
  return marker;
}

bool test_constant_correction_is_a_reference_shift(MPI_Comm world, int rank,
                                                   int size) {
  Fixture fixture;
  bool passed = expect(make_fixture(world, fixture), rank,
                       "closed-gauge prepare fixture compiles");
  if (!passed) return false;

  const Int3 cells = fixture.patch.cells;
  OwnedField pi = make_field(kPressure, cells, 701U + rank, 710U + rank,
                             720U + rank);
  OwnedField correction = make_field(20U, cells, 702U + rank, 730U + rank,
                                     740U + rank);
  OwnedField chi = make_field(kCompressibility, cells, 703U + rank,
                              750U + rank, 760U + rank);
  std::fill(pi.values.begin(), pi.values.end(), 0.0);
  std::fill(correction.values.begin(), correction.values.end(), 7.25);
  for (std::size_t cell = 0U; cell < chi.values.size(); ++cell) {
    chi.values[cell] = 1.0e-5 * (1.0 + static_cast<double>(cell % 7U));
  }
  std::vector<std::uint8_t> active(pi.values.size(), 1U);
  const PressureEnergyCellActivity activity{
      {active.data(), active.size()},
      static_cast<PlanFingerprint>(770U + rank), 780U};
  const ClosedGaugeCorrectionPrepareInput input{
      predecessor(fixture),
      100.0,
      1U,
      kTime,
      fixture.geometry.topology_revision(),
      static_cast<RevisionToken>(781U + rank),
      static_cast<PlanFingerprint>(782U + rank),
      as_const(pi.view),
      as_const(correction.view),
      as_const(chi.view),
      activity};
  const std::vector<double> pi_before = pi.values;
  const std::vector<double> correction_before = correction.values;
  const std::vector<double> chi_before = chi.values;
  ClosedGaugeCorrectionCertificate certificate;
  const PressureReferencePlan& pressure_reference =
      fixture.equations.pressure_reference();
  const Status status = pressure_reference.prepare_closed_gauge_correction(
      input, fixture.reductions, certificate);
  passed &= expect(
      pressure_reference.pressure_perturbation_field() == kPressure &&
          status && certificate.valid() && close(certificate.shift, 7.25) &&
          close(certificate.next_pressure_reference, 107.25) &&
          certificate.output_pressure_reference.valid() &&
          certificate.output_pressure_reference.pressure_reference ==
              certificate.rank_local_transaction &&
          fixture.equations.pressure_reference()
              .matches_closed_gauge_correction(input, certificate) &&
          pi.values == pi_before && correction.values == correction_before &&
          chi.values == chi_before,
      rank,
      "constant raw dp is prepared as a read-only absolute-pressure-preserving reference shift");

  unsigned long long local_collective =
      static_cast<unsigned long long>(certificate.collective_transaction);
  unsigned long long minimum_collective = 0U;
  unsigned long long maximum_collective = 0U;
  unsigned long long local_rank_token =
      static_cast<unsigned long long>(certificate.rank_local_transaction);
  unsigned long long minimum_rank_token = 0U;
  unsigned long long maximum_rank_token = 0U;
  MPI_Allreduce(&local_collective, &minimum_collective, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MIN, world);
  MPI_Allreduce(&local_collective, &maximum_collective, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MAX, world);
  MPI_Allreduce(&local_rank_token, &minimum_rank_token, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MIN, world);
  MPI_Allreduce(&local_rank_token, &maximum_rank_token, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MAX, world);
  passed &= expect(minimum_collective == maximum_collective &&
                       minimum_collective != 0U &&
                       (size == 1 || minimum_rank_token != maximum_rank_token),
                   rank,
                   "collective transaction excludes rank-local storage while local transaction binds it");
  return passed;
}

bool test_candidate_compressibility_weights_the_target_gauge(MPI_Comm world,
                                                              int rank) {
  Fixture fixture;
  bool passed = expect(make_fixture(world, fixture), rank,
                       "candidate-weight gauge fixture compiles");
  if (!passed) return false;

  const Int3 cells = fixture.patch.cells;
  OwnedField pi = make_field(kPressure, cells, 801U + rank, 810U + rank,
                             820U + rank);
  OwnedField correction = make_field(20U, cells, 802U + rank, 830U + rank,
                                     840U + rank);
  OwnedField candidate_chi = make_field(
      kCompressibility, cells, 803U + rank, 850U + rank, 860U + rank);
  OwnedField uniform_chi = make_field(kCompressibility, cells, 804U + rank,
                                      870U + rank, 880U + rank);
  bool decomposition_preserves_absolute_pressure = true;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const std::int32_t global_x = fixture.patch.begin.x + x;
        pi.view.unchecked(cell, 0U) = 0.0;
        correction.view.unchecked(cell, 0U) = global_x < 4 ? 2.0 : 8.0;
        candidate_chi.view.unchecked(cell, 0U) =
            global_x < 4 ? 1.0e-5 : 3.0e-5;
        uniform_chi.view.unchecked(cell, 0U) = 1.0e-5;
      }
    }
  }
  std::vector<std::uint8_t> active(pi.values.size(), 1U);
  const PressureEnergyCellActivity activity{
      {active.data(), active.size()},
      static_cast<PlanFingerprint>(890U + rank), 900U};
  ClosedGaugeCorrectionPrepareInput candidate_input{
      predecessor(fixture),
      100.0,
      2U,
      kTime,
      fixture.geometry.topology_revision(),
      static_cast<RevisionToken>(901U + rank),
      static_cast<PlanFingerprint>(902U + rank),
      as_const(pi.view),
      as_const(correction.view),
      as_const(candidate_chi.view),
      activity};
  ClosedGaugeCorrectionCertificate candidate_certificate;
  const Status candidate_status = fixture.equations.pressure_reference()
                                      .prepare_closed_gauge_correction(
                                          candidate_input,
                                          fixture.reductions,
                                          candidate_certificate);
  passed &= expect(
      candidate_status && candidate_certificate.valid() &&
          close(candidate_certificate.shift, 6.5) &&
          candidate_certificate.post_shift_gauge_residual <=
              candidate_certificate.post_shift_gauge_tolerance &&
          fixture.equations.pressure_reference()
              .matches_closed_gauge_correction(candidate_input,
                                                candidate_certificate),
      rank,
      "target candidate chi weights q={2,8} as (2+3*8)/(1+3)=6.5");
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const bool lower = fixture.patch.begin.x + x < 4;
        const double next_pi =
            (pi.view.unchecked(cell, 0U) +
             correction.view.unchecked(cell, 0U)) -
            candidate_certificate.shift;
        const double next_absolute =
            candidate_certificate.next_pressure_reference + next_pi;
        const double raw_absolute =
            candidate_input.pressure_reference +
            pi.view.unchecked(cell, 0U) +
            correction.view.unchecked(cell, 0U);
        decomposition_preserves_absolute_pressure &=
            close(next_pi, lower ? -4.5 : 1.5) &&
            close(next_absolute, raw_absolute);
      }
    }
  }
  passed &= expect(
      decomposition_preserves_absolute_pressure, rank,
      "prepared decomposition gives analytic zero-mean pi and preserves every absolute pressure");

  ClosedGaugeCorrectionPrepareInput uniform_input = candidate_input;
  uniform_input.candidate_pressure_compressibility =
      as_const(uniform_chi.view);
  ClosedGaugeCorrectionCertificate uniform_certificate;
  const Status uniform_status = fixture.equations.pressure_reference()
                                    .prepare_closed_gauge_correction(
                                        uniform_input, fixture.reductions,
                                        uniform_certificate);
  passed &= expect(
      uniform_status && uniform_certificate.valid() &&
          close(uniform_certificate.shift, 5.0) &&
          uniform_certificate.collective_transaction !=
              candidate_certificate.collective_transaction &&
          uniform_certificate.rank_local_transaction !=
              candidate_certificate.rank_local_transaction &&
          !fixture.equations.pressure_reference()
               .matches_closed_gauge_correction(uniform_input,
                                                 candidate_certificate),
      rank,
      "changing target chi changes the weighted shift and invalidates the prior transaction");
  return passed;
}

bool test_empty_local_activity_and_all_cell_positivity(MPI_Comm world,
                                                       int rank) {
  Fixture fixture;
  bool passed = expect(make_fixture(world, fixture), rank,
                       "empty-local-active gauge fixture compiles");
  if (!passed) return false;

  const Int3 cells = fixture.patch.cells;
  OwnedField pi = make_field(kPressure, cells, 901U + rank, 910U + rank,
                             920U + rank);
  OwnedField correction = make_field(20U, cells, 902U + rank, 930U + rank,
                                     940U + rank);
  OwnedField chi = make_field(kCompressibility, cells, 903U + rank,
                              950U + rank, 960U + rank);
  std::fill(pi.values.begin(), pi.values.end(), 0.0);
  std::fill(correction.values.begin(), correction.values.end(),
            rank == 0 ? 4.0 : 9.0);
  std::fill(chi.values.begin(), chi.values.end(), 2.0e-5);
  std::vector<std::uint8_t> active(pi.values.size(), rank == 0 ? 1U : 0U);
  const PressureEnergyCellActivity activity{
      {active.data(), active.size()},
      static_cast<PlanFingerprint>(970U + rank), 980U};
  ClosedGaugeCorrectionPrepareInput input{
      predecessor(fixture),
      100.0,
      1U,
      kTime,
      fixture.geometry.topology_revision(),
      static_cast<RevisionToken>(994U + rank),
      static_cast<PlanFingerprint>(995U + rank),
      as_const(pi.view),
      as_const(correction.view),
      as_const(chi.view),
      activity};
  ClosedGaugeCorrectionCertificate certificate;
  const Status status = fixture.equations.pressure_reference()
                            .prepare_closed_gauge_correction(
                                input, fixture.reductions, certificate);
  passed &= expect(
      status && certificate.valid() && close(certificate.shift, 4.0) &&
          ((rank == 0 && certificate.local_active_cells == pi.values.size() &&
            certificate.local_weight > 0.0) ||
           (rank != 0 && certificate.local_active_cells == 0U &&
            certificate.local_moment == 0.0 &&
            certificate.local_weight == 0.0)),
      rank,
      "an empty local fluid partition contributes exact zero while global weight remains positive");
  bool inactive_decomposition = true;
  if (rank != 0) {
    for (std::size_t cell = 0U; cell < pi.values.size(); ++cell) {
      const double next_pi =
          (pi.values[cell] + correction.values[cell]) - certificate.shift;
      inactive_decomposition &=
          close(next_pi, 5.0) &&
          close(certificate.next_pressure_reference + next_pi,
                input.pressure_reference + pi.values[cell] +
                    correction.values[cell]);
    }
  }
  passed &= expect(inactive_decomposition, rank,
                   "the same gauge shift is prepared for every inactive IBM cell");

  const ClosedGaugeCorrectionCertificate marker = marker_certificate();
  ClosedGaugeCorrectionCertificate rejected = marker;
  std::fill(active.begin(), active.end(), 0U);
  input.activity.local_fingerprint =
      static_cast<PlanFingerprint>(990U + rank);
  input.activity.collective_fingerprint = 991U;
  const Status zero_weight = fixture.equations.pressure_reference()
                                 .prepare_closed_gauge_correction(
                                     input, fixture.reductions, rejected);
  passed &= expect(!zero_weight &&
                       zero_weight.code == StatusCode::numerical_failure &&
                       same_collective_status(world, zero_weight) &&
                       same_gauge_certificate(rejected, marker),
                   rank,
                   "a globally empty activity map rejects atomically because gauge weight is zero");

  std::fill(active.begin(), active.end(), rank == 0 ? 1U : 0U);
  if (rank == 0 && !active.empty()) active[0U] = 0U;
  input.activity.local_fingerprint =
      static_cast<PlanFingerprint>(992U + rank);
  input.activity.collective_fingerprint = 993U;
  correction.values[0U] = -101.0;
  rejected = marker;
  const Status bad_absolute = fixture.equations.pressure_reference()
                                  .prepare_closed_gauge_correction(
                                      input, fixture.reductions, rejected);
  passed &= expect(!bad_absolute &&
                       bad_absolute.code == StatusCode::numerical_failure &&
                       same_collective_status(world, bad_absolute) &&
                       same_gauge_certificate(rejected, marker),
                   rank,
                   "non-positive candidate absolute pressure rejects even in an inactive cell");
  return passed;
}

bool test_mutation_and_collective_fault_atomicity(MPI_Comm world, int rank,
                                                  int size) {
  Fixture fixture;
  bool passed = expect(make_fixture(world, fixture), rank,
                       "gauge mutation/fault fixture compiles");
  if (!passed) return false;

  const Int3 cells = fixture.patch.cells;
  OwnedField pi = make_field(kPressure, cells, 1001U + rank, 1010U + rank,
                             1020U + rank);
  OwnedField correction = make_field(20U, cells, 1002U + rank,
                                     1030U + rank, 1040U + rank);
  OwnedField chi = make_field(kCompressibility, cells, 1003U + rank,
                              1050U + rank, 1060U + rank);
  std::fill(pi.values.begin(), pi.values.end(), 0.25);
  std::fill(correction.values.begin(), correction.values.end(), 1.75);
  std::fill(chi.values.begin(), chi.values.end(), 2.0e-5);
  std::vector<std::uint8_t> active(pi.values.size(), 1U);
  const PressureEnergyCellActivity activity{
      {active.data(), active.size()},
      static_cast<PlanFingerprint>(1070U + rank), 1080U};
  ClosedGaugeCorrectionPrepareInput input{
      predecessor(fixture),
      100.0,
      1U,
      kTime,
      fixture.geometry.topology_revision(),
      static_cast<RevisionToken>(1081U + rank),
      static_cast<PlanFingerprint>(1082U + rank),
      as_const(pi.view),
      as_const(correction.view),
      as_const(chi.view),
      activity};
  ClosedGaugeCorrectionCertificate accepted;
  const PressureReferencePlan& plan =
      fixture.equations.pressure_reference();
  const bool prepared = static_cast<bool>(plan.prepare_closed_gauge_correction(
      input, fixture.reductions, accepted));
  const LinearReductionCounters before_match = fixture.reductions.counters();
  const bool baseline_matches =
      plan.matches_closed_gauge_correction(input, accepted);
  const LinearReductionCounters after_match = fixture.reductions.counters();
  passed &= expect(prepared && baseline_matches &&
                       before_match.calls == after_match.calls &&
                       before_match.scalars == after_match.scalars &&
                       before_match.blocking_operations ==
                           after_match.blocking_operations,
                   rank, "baseline gauge transaction prepares and matches");
  const ClosedGaugeCorrectionCertificate marker = marker_certificate();

  const double correction_saved = correction.values[0U];
  correction.values[0U] += 0.125;
  passed &= expect(!plan.matches_closed_gauge_correction(input, accepted),
                   rank,
                   "raw-dp data mutation without a revision bump is detected");
  correction.values[0U] = correction_saved;

  const double chi_saved = chi.values[0U];
  chi.values[0U] *= 1.5;
  passed &= expect(!plan.matches_closed_gauge_correction(input, accepted),
                   rank,
                   "candidate-chi data mutation without a revision bump is detected");
  chi.values[0U] = chi_saved;

  active[0U] = 0U;
  passed &= expect(!plan.matches_closed_gauge_correction(input, accepted),
                   rank,
                   "IBM activity content mutation invalidates the transaction");
  active[0U] = 1U;

  ClosedGaugeCorrectionPrepareInput mutated_input = input;
  ++mutated_input.predecessor.pressure_reference;
  passed &= expect(
      !plan.matches_closed_gauge_correction(mutated_input, accepted), rank,
      "predecessor pressure-reference revision invalidates the transaction");
  mutated_input = input;
  mutated_input.corrector = 2U;
  passed &= expect(
      !plan.matches_closed_gauge_correction(mutated_input, accepted), rank,
      "corrector identity invalidates the transaction");
  mutated_input = input;
  ++mutated_input.raw_pressure_correction.storage_identity;
  passed &= expect(
      !plan.matches_closed_gauge_correction(mutated_input, accepted), rank,
      "rank-local raw-dp storage substitution invalidates the transaction");
  mutated_input = input;
  ++mutated_input.pressure_correction_authority;
  passed &= expect(
      !plan.matches_closed_gauge_correction(mutated_input, accepted), rank,
      "pressure-correction producer authority invalidates the transaction");
  mutated_input = input;
  ++mutated_input.target_thermodynamic_closure;
  passed &= expect(
      !plan.matches_closed_gauge_correction(mutated_input, accepted), rank,
      "same-target EOS/thermo closure invalidates the transaction");

  ClosedGaugeCorrectionCertificate structural_failure = marker;
  mutated_input = input;
  mutated_input.pressure_perturbation.field = kPressure + 1U;
  const Status wrong_pressure_field = plan.prepare_closed_gauge_correction(
      mutated_input, fixture.reductions, structural_failure);
  passed &= expect(
      !wrong_pressure_field &&
          wrong_pressure_field.code == StatusCode::invalid_plan &&
          same_collective_status(world, wrong_pressure_field) &&
          same_gauge_certificate(structural_failure, marker),
      rank,
      "a foreign scalar cannot substitute for the compiled pi field");
  structural_failure = marker;
  mutated_input = input;
  mutated_input.target_thermodynamic_closure = 0U;
  const Status missing_thermodynamic_lineage =
      plan.prepare_closed_gauge_correction(mutated_input, fixture.reductions,
                                           structural_failure);
  passed &= expect(
      !missing_thermodynamic_lineage &&
          missing_thermodynamic_lineage.code == StatusCode::invalid_plan &&
          same_collective_status(world, missing_thermodynamic_lineage) &&
          same_gauge_certificate(structural_failure, marker),
      rank, "a bare candidate chi without target EOS lineage is rejected");

  ClosedGaugeCorrectionCertificate mutated_certificate = accepted;
  mutated_certificate.shift += 0.25;
  passed &= expect(
      !plan.matches_closed_gauge_correction(input, mutated_certificate), rank,
      "shift mutation invalidates the transaction");
  mutated_certificate = accepted;
  ++mutated_certificate.output_pressure_reference.pressure_reference;
  passed &= expect(
      !plan.matches_closed_gauge_correction(input, mutated_certificate), rank,
      "output pressure-reference mutation invalidates the transaction");
  mutated_certificate = accepted;
  mutated_certificate.global_weight *= 2.0;
  passed &= expect(
      !plan.matches_closed_gauge_correction(input, mutated_certificate), rank,
      "global-weight mutation invalidates the transaction");

  ClosedGaugeCorrectionCertificate corrected_layer;
  mutated_input = input;
  mutated_input.corrector = 2U;
  passed &= expect(
      plan.prepare_closed_gauge_correction(mutated_input, fixture.reductions,
                                           corrected_layer) &&
          corrected_layer.collective_transaction !=
              accepted.collective_transaction &&
          corrected_layer.rank_local_transaction !=
              accepted.rank_local_transaction,
      rank, "corrector changes both collective and rank-local lineage");

  ClosedGaugeCorrectionCertificate revised_predecessor;
  mutated_input = input;
  ++mutated_input.predecessor.pressure_reference;
  passed &= expect(
      plan.prepare_closed_gauge_correction(mutated_input, fixture.reductions,
                                           revised_predecessor) &&
          revised_predecessor.collective_transaction ==
              accepted.collective_transaction &&
          revised_predecessor.rank_local_transaction !=
              accepted.rank_local_transaction &&
          revised_predecessor.output_pressure_reference.pressure_reference ==
              revised_predecessor.rank_local_transaction,
      rank,
      "exact predecessor revision is rank-local lineage and cannot alias the old output authority");

  ClosedGaugeCorrectionCertificate revised_pressure_authority;
  mutated_input = input;
  ++mutated_input.pressure_correction_authority;
  passed &= expect(
      plan.prepare_closed_gauge_correction(
          mutated_input, fixture.reductions, revised_pressure_authority) &&
          revised_pressure_authority.collective_transaction ==
              accepted.collective_transaction &&
          revised_pressure_authority.rank_local_transaction !=
              accepted.rank_local_transaction,
      rank,
      "exact pressure-correction authority is bound in rank-local lineage");

  ClosedGaugeCorrectionCertificate revised_thermodynamic_closure;
  mutated_input = input;
  ++mutated_input.target_thermodynamic_closure;
  passed &= expect(
      plan.prepare_closed_gauge_correction(
          mutated_input, fixture.reductions,
          revised_thermodynamic_closure) &&
          revised_thermodynamic_closure.collective_transaction ==
              accepted.collective_transaction &&
          revised_thermodynamic_closure.rank_local_transaction !=
              accepted.rank_local_transaction,
      rank, "same-target EOS/thermo closure is bound in rank-local lineage");

  ClosedGaugeCorrectionCertificate failed = marker;
  const std::vector<double> pi_before = pi.values;
  const std::vector<double> correction_before = correction.values;
  const std::vector<double> chi_before = chi.values;
  const int failing_rank = size - 1;
  const Status armed = fixture.reductions.arm_checked_sum_fault_for_test(
      1U, failing_rank);
  const Status fault = plan.prepare_closed_gauge_correction(
      input, fixture.reductions, failed);
  passed &= expect(
      armed && !fault && fault.code == StatusCode::numerical_failure &&
          same_collective_status(world, fault) &&
          fixture.reductions.lowest_failing_rank() == failing_rank &&
          same_gauge_certificate(failed, marker) && pi.values == pi_before &&
          correction.values == correction_before && chi.values == chi_before,
      rank,
      "rank-selective checked-sum fault rejects collectively without publishing certificate or mutating fields");

  ClosedGaugeCorrectionCertificate retry;
  passed &= expect(
      plan.prepare_closed_gauge_correction(input, fixture.reductions, retry) &&
          plan.matches_closed_gauge_correction(input, retry),
      rank, "one-shot collective failure can retry immediately");

  failed = marker;
  const Status post_shift_armed =
      fixture.reductions.arm_checked_sum_fault_for_test(2U, failing_rank);
  const Status post_shift_fault = plan.prepare_closed_gauge_correction(
      input, fixture.reductions, failed);
  passed &= expect(
      post_shift_armed && !post_shift_fault &&
          post_shift_fault.code == StatusCode::numerical_failure &&
          same_collective_status(world, post_shift_fault) &&
          fixture.reductions.lowest_failing_rank() == failing_rank &&
          same_gauge_certificate(failed, marker) && pi.values == pi_before &&
          correction.values == correction_before && chi.values == chi_before,
      rank,
      "rank-selective post-shift audit fault is collective and failure-atomic");
  ClosedGaugeCorrectionCertificate post_shift_retry;
  passed &= expect(
      plan.prepare_closed_gauge_correction(
          input, fixture.reductions, post_shift_retry) &&
          plan.matches_closed_gauge_correction(input, post_shift_retry),
      rank, "post-shift audit failure can retry immediately");

  if (size > 1) {
    ClosedGaugeCorrectionPrepareInput divergent = input;
    if (rank == size - 1) divergent.corrector = 2U;
    failed = marker;
    const Status mismatch = plan.prepare_closed_gauge_correction(
        divergent, fixture.reductions, failed);
    passed &= expect(
        !mismatch && mismatch.code == StatusCode::invalid_plan &&
            same_collective_status(world, mismatch) &&
            same_gauge_certificate(failed, marker),
        rank,
        "rank-divergent collective semantics reject before the gauge sum");
  }
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  bool passed = test_constant_correction_is_a_reference_shift(
      MPI_COMM_WORLD, rank, size);
  passed &= test_candidate_compressibility_weights_the_target_gauge(
      MPI_COMM_WORLD, rank);
  passed &= test_empty_local_activity_and_all_cell_positivity(
      MPI_COMM_WORLD, rank);
  passed &= test_mutation_and_collective_fault_atomicity(
      MPI_COMM_WORLD, rank, size);
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
