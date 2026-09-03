// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_flow.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr FieldId kDensity = 0U;
constexpr FieldId kVelocity = 1U;
constexpr FieldId kPressure = 2U;
constexpr FieldId kEnthalpy = 3U;
constexpr FieldId kTemperature = 4U;
constexpr FieldId kSpecies = 5U;
constexpr FieldId kEffectiveViscosity = 6U;
constexpr FieldId kCompressibility = 7U;
constexpr FieldId kVelocityGradient = 8U;
constexpr double kPi = 3.141592653589793238462643383279502884;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <=
             2.0e-12 * std::max(1.0, std::abs(expected));
}

struct OwnedField {
  std::vector<double> bytes;
  FieldView view{};
};

OwnedField make_field(FieldId id, Int3 cells, std::uint8_t ghosts,
                      RevisionToken revision) {
  OwnedField field;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  field.bytes.assign(nx * ny * nz, 0.0);
  field.view.base =
      field.bytes.data() + ghosts + ghosts * nx + ghosts * nx * ny;
  field.view.interior = cells;
  field.view.ghosts = {ghosts, ghosts, ghosts};
  field.view.components = 1U;
  field.view.stride_y = nx;
  field.view.stride_z = nx * ny;
  field.view.component_stride = nx * ny * nz;
  field.view.field = id;
  field.view.revision = revision;
  field.view.storage_identity = 1000U + revision;
  field.view.revision_domain = 7001U;
  return field;
}

void fill(OwnedField& field, double value) {
  std::fill(field.bytes.begin(), field.bytes.end(), value);
}

CartesianMeshSpec mesh_spec(std::int32_t n, bool stretched = false) {
  CartesianMeshSpec mesh;
  mesh.kind = stretched ? GeometryKind::tensor_stretched
                        : GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {n, n, n};
  const double inverse_n = 1.0 / static_cast<double>(n);
  if (stretched) {
    mesh.has_base_spacing = true;
    mesh.base_spacing = {1.08 * inverse_n, 1.08 * inverse_n,
                         1.08 * inverse_n};
    mesh.minimum_spacing = {0.92 * inverse_n, 0.92 * inverse_n,
                            0.92 * inverse_n};
    mesh.max_growth_ratio = 1.0 + 0.8 * inverse_n;
    mesh.focus_regions.push_back(
        {{0.35, 0.35, 0.35}, {0.65, 0.65, 0.65},
         {0.96 * inverse_n, 0.96 * inverse_n, 0.96 * inverse_n}});
  } else {
    mesh.minimum_spacing = {inverse_n, inverse_n, inverse_n};
    mesh.max_growth_ratio = 1.0;
  }
  mesh.limits.max_global_cells = static_cast<std::uint64_t>(n) * n * n;
  mesh.limits.max_memory_bytes_per_rank = 1U << 26U;
  return mesh;
}

SpeciesThermophysicalSpec species_spec(std::string_view name,
                                       double molecular_weight) {
  SpeciesThermophysicalSpec species;
  species.stable_name.assign(name.data(), name.size());
  species.molecular_weight = molecular_weight;
  species.temperature_switch = 1000.0;
  species.nasa7_low[0U] = 3.5;
  species.nasa7_high[0U] = 3.5;
  species.viscosity_reference = 1.8e-5;
  species.conductivity = 0.026;
  return species;
}

ThermophysicalSpec thermophysical_spec() {
  ThermophysicalSpec spec;
  spec.data_file = "analytic.d";
  spec.minimum_temperature = 200.0;
  spec.maximum_temperature = 2000.0;
  spec.temperature_relative_tolerance = 1.0e-12;
  spec.maximum_temperature_iterations = 64U;
  spec.closed_mass_relative_tolerance = 1.0e-12;
  spec.maximum_closed_mass_iterations = 32U;
  spec.maximum_closed_mass_relative_step = 0.2;
  spec.species.push_back(species_spec("species_a", 28.0));
  spec.species.push_back(species_spec("species_b", 32.0));
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
};

bool make_fixture(std::int32_t n, PressureReferenceKind pressure_reference,
                  Fixture& out, bool stretched = false) {
  const CartesianMeshSpec mesh = mesh_spec(n, stretched);
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x14a7U;
  model.pressure_reference = pressure_reference;
  model.transported_scalars = {
      {"species_a", TransportedScalarRole::species, 1.0, 1.0}};
  if (pressure_reference == PressureReferenceKind::closed_mass) {
    for (BoundaryFaceSpec& face : model.boundaries) {
      face.flow_kind = BoundaryKind::periodic;
      face.thermal_kind = BoundaryKind::none;
    }
  } else {
    for (BoundaryFaceSpec& face : model.boundaries) {
      face.flow_kind = BoundaryKind::no_slip_wall;
      face.thermal_kind = BoundaryKind::adiabatic_wall;
      face.scalars.push_back(ScalarBoundarySpec{
          "species_a", ScalarBoundaryKind::zero_gradient, 0.0,
          ScalarBoundaryKind::zero_gradient, 0.0});
    }
    model.boundaries[0U].flow_kind = BoundaryKind::velocity_inlet;
    model.boundaries[0U].thermal_kind = BoundaryKind::none;
    model.boundaries[0U].velocity = {0.3, 0.0, 0.0};
    model.boundaries[0U].temperature = 300.0;
    model.boundaries[0U].scalars[0U].kind =
        ScalarBoundaryKind::dirichlet;
    model.boundaries[0U].scalars[0U].value = 0.2;
    model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
    model.boundaries[1U].thermal_kind = BoundaryKind::none;
    model.boundaries[1U].pressure = 101325.0;
  }
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::central2;
  model.schemes.species = ConvectionScheme::central2;
  model.schemes.passive_scalar = ConvectionScheme::central2;

  FieldRegistry registry;
  FieldId id = 0U;
  if (!registry.require_field("rho", 1U, 2U, id) || id != kDensity ||
      !registry.require_field("U", 3U, 2U, id) || id != kVelocity ||
      !registry.require_field("pi", 1U, 2U, id) || id != kPressure ||
      !registry.require_field("h", 1U, 2U, id) || id != kEnthalpy ||
      !registry.require_field("T", 1U, 2U, id) || id != kTemperature ||
      !CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh, {},
                                          out.geometry, out.patch) ||
      !BoundaryCompiler::compile(MPI_COMM_SELF, model, out.geometry,
                                 out.patch, registry, out.boundary,
                                 out.schemes, out.time)) {
    return false;
  }

  const ThermophysicalSpec thermophysics = thermophysical_spec();
  if (!ThermodynamicsPlan::compile(
          thermophysics,
          {model.transported_scalars.data(), model.transported_scalars.size()},
          out.thermodynamics) ||
      !TransportPlan::compile(thermophysics, out.thermodynamics,
                              out.transport)) {
    return false;
  }
  const std::array<FieldId, 9U> declared{
      kDensity,          kVelocity,      kPressure,
      kEnthalpy,         kTemperature,   kSpecies,
      kEffectiveViscosity, kCompressibility, kVelocityGradient};
  if (!out.contributions.configure({declared.data(), declared.size()}) ||
      !out.contributions.freeze()) {
    return false;
  }
  const std::array<ScalarEquationSpec, 1U> scalars{{
      {kSpecies, TransportedScalarRole::species, 1.0, 1.0},
  }};
  EquationPlanSpec spec;
  spec.density = kDensity;
  spec.velocity = kVelocity;
  spec.pressure_perturbation = kPressure;
  spec.enthalpy = kEnthalpy;
  spec.temperature = kTemperature;
  spec.effective_viscosity = kEffectiveViscosity;
  spec.pressure_compressibility = kCompressibility;
  spec.velocity_gradient = kVelocityGradient;
  spec.pressure_reference = pressure_reference;
  spec.scalars = {scalars.data(), scalars.size()};
  spec.closed_mass_service_stage =
      pressure_reference == PressureReferenceKind::closed_mass ? 1U : 0U;
  spec.maximum_cells_per_rank = static_cast<std::size_t>(n) * n * n;
  return static_cast<bool>(EquationPlanSet::compile(
      MPI_COMM_SELF, out.schemes, out.geometry, out.patch, out.boundary,
      out.contributions, out.thermodynamics, out.transport, spec,
      out.equations));
}

struct FluxHistoryFixture {
  FieldId dependency{};
  StateLayers layers;
  AttemptTransaction transaction;
  FaceFluxStorage storage;
  FinalFaceFluxAuthority authority;
  FinalFaceFluxWriter writer;
};

bool make_flux_history(Int3 cells, FluxHistoryFixture& fixture) {
  FieldRegistry registry;
  FieldSchema schema;
  if (!registry.declare_field("flux_dependency", 1U, 0U,
                              fixture.dependency) ||
      !registry.freeze(schema)) {
    return false;
  }
  const std::array requests{ArenaFieldRequest{
      fixture.dependency, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout layout;
  return static_cast<bool>(ArenaLayout::compile(
             schema, {requests.data(), requests.size()}, layout)) &&
         static_cast<bool>(StateLayers::allocate(layout, fixture.layers)) &&
         static_cast<bool>(AttemptTransaction::create(
             fixture.layers.field_count(), 1U, fixture.layers.field_count(),
             fixture.transaction)) &&
         static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                           fixture.storage)) &&
         static_cast<bool>(fixture.authority.claim(
             41U, 0U, fixture.transaction, fixture.writer));
}

bool commit_uniform_flux(const CartesianKernelPlan& kernels, Int3 cells,
                         double velocity_x, FluxHistoryFixture& fixture,
                         ConstFaceFluxView& committed) {
  if (!fixture.transaction.begin(fixture.layers) ||
      !fixture.transaction.revise_trial(fixture.dependency)) {
    return false;
  }
  const RevisionDependency dependency{
      AttemptTransaction::field_revision_source(fixture.dependency),
      fixture.transaction.trial_revision(fixture.dependency)};
  PendingFaceFluxView pending;
  if (!fixture.writer.begin_pending(fixture.transaction, fixture.storage,
                                    pending)) {
    return false;
  }
  OwnedField rho = make_field(40U, cells, 2U, 501U);
  OwnedField velocity_x_field = make_field(41U, cells, 2U, 502U);
  // Reconstruct expects one three-component velocity field.  Reuse one
  // contiguous backing with explicit component strides.
  const std::size_t scalar_size = velocity_x_field.bytes.size();
  velocity_x_field.bytes.resize(3U * scalar_size, 0.0);
  velocity_x_field.view.base = velocity_x_field.bytes.data() + 2U +
                               2U * velocity_x_field.view.stride_y +
                               2U * velocity_x_field.view.stride_z;
  velocity_x_field.view.components = 3U;
  velocity_x_field.view.component_stride = scalar_size;
  fill(rho, 1.0);
  for (std::int32_t z = -2; z < cells.z + 2; ++z) {
    for (std::int32_t y = -2; y < cells.y + 2; ++y) {
      for (std::int32_t x = -2; x < cells.x + 2; ++x) {
        velocity_x_field.view.unchecked({x, y, z}, 0U) = velocity_x;
      }
    }
  }
  const std::array<ConstFieldView, 2U> reads{
      as_const(rho.view), as_const(velocity_x_field.view)};
  const KernelInvocation call{{reads.data(), reads.size()}, {},
                              {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U, nullptr};
  const std::array dependencies{dependency};
  return static_cast<bool>(reconstruct_mass_flux(kernels, call, pending)) &&
         static_cast<bool>(fixture.writer.publish_pending(
             {dependencies.data(), dependencies.size()}, pending)) &&
         static_cast<bool>(fixture.transaction.collective_finish(MPI_COMM_SELF,
                                                                  Status{})) &&
         static_cast<bool>(fixture.writer.committed(fixture.storage,
                                                     committed));
}

bool same_certificate(ThermophysicalPredictorCertificate left,
                      ThermophysicalPredictorCertificate right) {
  return left.plan == right.plan && left.time == right.time &&
         left.geometry == right.geometry &&
         left.accepted_face_flux == right.accepted_face_flux &&
         left.previous_face_flux == right.previous_face_flux &&
         left.committed_face_flux_authority ==
             right.committed_face_flux_authority &&
         left.committed_face_flux_storage ==
             right.committed_face_flux_storage &&
         left.committed_face_flux_revision_domain ==
             right.committed_face_flux_revision_domain &&
         left.state == right.state && left.order == right.order;
}

ThermophysicalGhostAuthority make_ghost_authority(
    ConstFieldView view, RevisionToken geometry, RevisionToken boundary,
    std::uintptr_t exchange_plan) {
  return {exchange_plan,
          view.field,
          view.revision,
          view.storage_identity,
          view.revision_domain,
          geometry,
          boundary,
          static_cast<std::uint8_t>(view.ghosts.x)};
}

void attach_ghost_authority(
    const Fixture& fixture, ThermophysicalPredictorInput& input,
    std::array<ThermophysicalGhostHistory, 1U>& species_ghosts,
    std::uintptr_t exchange_plan = static_cast<std::uintptr_t>(0xC0A57U)) {
  const RevisionToken geometry = fixture.geometry.topology_revision();
  const RevisionToken boundary = fixture.boundary.revision();
  input.enthalpy_ghosts.accepted = make_ghost_authority(
      input.enthalpy_accepted, geometry, boundary, exchange_plan);
  input.enthalpy_ghosts.previous =
      input.bdf.order == 2U
          ? make_ghost_authority(input.enthalpy_previous, geometry, boundary,
                                 exchange_plan)
          : ThermophysicalGhostAuthority{};
  if (input.species_accepted.size != 0U) {
    species_ghosts[0U].accepted = make_ghost_authority(
        input.species_accepted.data[0U], geometry, boundary, exchange_plan);
    species_ghosts[0U].previous =
        input.bdf.order == 2U
            ? make_ghost_authority(input.species_previous.data[0U], geometry,
                                   boundary, exchange_plan)
            : ThermophysicalGhostAuthority{};
  } else {
    species_ghosts[0U] = {};
  }
  input.species_ghosts = {species_ghosts.data(), input.species_accepted.size};
  input.passive_scalar_ghosts = {};
}

double periodic_center(const AxisMetrics& axis, std::int32_t index,
                       std::int32_t cells) {
  const std::int32_t wrapped = (index % cells + cells) % cells;
  const std::int32_t periods = (index - wrapped) / cells;
  return axis.centres().data[static_cast<std::size_t>(wrapped)] +
         static_cast<double>(periods) *
             (axis.faces().data[static_cast<std::size_t>(cells)] -
              axis.faces().data[0U]);
}

double run_predictor_mms(std::int32_t n, bool stretched) {
  Fixture fixture;
  if (!make_fixture(n, PressureReferenceKind::closed_mass, fixture,
                    stretched)) {
    return HUGE_VAL;
  }

  const Int3 cells = fixture.patch.cells;
  FluxHistoryFixture flux_history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  constexpr double velocity_x = 0.3;
  if (!make_flux_history(cells, flux_history) ||
      !commit_uniform_flux(fixture.equations.kernels(), cells, velocity_x,
                           flux_history, previous_flux) ||
      !commit_uniform_flux(fixture.equations.kernels(), cells, velocity_x,
                           flux_history, accepted_flux)) {
    return HUGE_VAL;
  }
  FaceFluxStorage paired_flux_storage;
  FaceFluxView paired_flux;
  if (!FaceFluxStorage::allocate_workspace(cells, 1U,
                                           paired_flux_storage) ||
      !paired_flux_storage.workspace_view(0U, 3016U, paired_flux)) {
    return HUGE_VAL;
  }
  OwnedField rho_n = make_field(kDensity, cells, 0U, 3001U);
  OwnedField rho_nm1 = make_field(kDensity, cells, 0U, 3002U);
  OwnedField h_n = make_field(kEnthalpy, cells, 2U, 3003U);
  OwnedField h_nm1 = make_field(kEnthalpy, cells, 2U, 3004U);
  OwnedField y_n = make_field(kSpecies, cells, 2U, 3005U);
  OwnedField y_nm1 = make_field(kSpecies, cells, 2U, 3006U);
  OwnedField h_star = make_field(kEnthalpy, cells, 0U, 3007U);
  OwnedField y_star = make_field(kSpecies, cells, 0U, 3008U);
  OwnedField rho_work = make_field(40U, cells, 0U, 3009U);
  OwnedField work_n = make_field(41U, cells, 0U, 3010U);
  OwnedField work_nm1 = make_field(42U, cells, 0U, 3011U);
  OwnedField low_rho = make_field(kDensity, cells, 0U, 3012U);
  OwnedField low_h = make_field(kEnthalpy, cells, 1U, 3013U);
  OwnedField low_y = make_field(kSpecies, cells, 1U, 3014U);
  fill(rho_n, 1.0);
  fill(rho_nm1, 1.0);
  fill(y_n, 0.2);
  fill(y_nm1, 0.2);
  constexpr double amplitude = 1000.0;
  constexpr double wave = 2.0 * kPi;
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const double x = periodic_center(fixture.geometry.x(), i, n);
        const double value = 300000.0 + amplitude * std::sin(wave * x);
        h_n.view.unchecked({i, j, k}, 0U) = value;
        h_nm1.view.unchecked({i, j, k}, 0U) = value;
      }
    }
  }
  const std::array<ConstFieldView, 1U> species_n{as_const(y_n.view)};
  const std::array<ConstFieldView, 1U> species_nm1{as_const(y_nm1.view)};
  const std::array<PredictorRateHistory, 1U> zero_species_rate{{{}}};
  std::array<FieldView, 1U> species_output{y_star.view};
  ThermophysicalPredictorInput input;
  input.dt = 0.1;
  input.bdf = {15.0, -20.0, 5.0, 2U};
  input.time = 3020U;
  input.geometry = fixture.geometry.topology_revision();
  input.boundary = fixture.boundary.revision();
  input.transport = fixture.transport.fingerprint();
  input.density_accepted = as_const(rho_n.view);
  input.density_previous = as_const(rho_nm1.view);
  input.enthalpy_accepted = as_const(h_n.view);
  input.enthalpy_previous = as_const(h_nm1.view);
  input.species_accepted = {species_n.data(), species_n.size()};
  input.species_previous = {species_nm1.data(), species_nm1.size()};
  input.mass_flux_accepted = accepted_flux;
  input.mass_flux_previous = previous_flux;
  input.species_nonadvective_rhs = {zero_species_rate.data(),
                                    zero_species_rate.size()};
  std::array<ThermophysicalGhostHistory, 1U> ghost_history{};
  attach_ghost_authority(fixture, input, ghost_history);
  ThermophysicalPredictorOutput output{
      h_star.view, {species_output.data(), species_output.size()},
      rho_work.view, work_n.view, work_nm1.view, {}, low_rho.view, low_h.view,
      {&low_y.view, 1U}, {}, paired_flux};
  ThermophysicalPredictorDiagnostics diagnostics;
  ThermophysicalPredictorCertificate certificate;
  if (!fixture.equations.thermophysical_predictor().predict(
          MPI_COMM_SELF, Status{}, input, output, {}, diagnostics,
          certificate)) {
    return HUGE_VAL;
  }
  long double error = 0.0L;
  long double measure = 0.0L;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const double x = fixture.geometry.x().centres().data[
            static_cast<std::size_t>(i)];
        const double correction =
            velocity_x * amplitude * wave * std::cos(wave * x) / 15.0;
        const double expected =
            300000.0 + amplitude * std::sin(wave * x) - correction;
        const double difference =
            h_star.view.unchecked({i, j, k}, 0U) - expected;
        error += difference * difference;
        measure += correction * correction;
      }
    }
  }
  return std::sqrt(static_cast<double>(error / measure));
}

double observed_order(double coarse, double fine) {
  return std::log(coarse / fine) / std::log(2.0);
}

bool test_predictor_mms_orders() {
  bool passed = true;
  for (const bool stretched : {false, true}) {
    const double e16 = run_predictor_mms(16, stretched);
    const double e32 = run_predictor_mms(32, stretched);
    const double e64 = run_predictor_mms(64, stretched);
    passed &= expect(std::isfinite(e16) && e16 > e32 && e32 > e64,
                     "committed-flux thermophysical predictor MMS converges");
    passed &= expect(observed_order(e16, e32) >= 1.8 &&
                         observed_order(e32, e64) >= 1.8,
                     "uniform/stretched predictor order >= 1.8");
  }
  return passed;
}

bool test_bdf2_predictor_and_mutations() {
  Fixture fixture;
  constexpr std::int32_t n = 4;
  bool passed = expect(make_fixture(n, PressureReferenceKind::closed_mass,
                                    fixture),
                       "predictor fixture compiles");
  if (!passed) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  FaceFluxStorage paired_flux_storage;
  FaceFluxView paired_flux;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, paired_flux_storage)) &&
          static_cast<bool>(paired_flux_storage.workspace_view(
              0U, 619U, paired_flux)),
      "attempt-local paired face-flux workspace allocates");
  if (!passed) {
    return false;
  }
  FluxHistoryFixture flux_history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(make_flux_history(cells, flux_history) &&
                       commit_uniform_flux(fixture.equations.kernels(), cells,
                                           0.2, flux_history,
                                           previous_flux) &&
                       commit_uniform_flux(fixture.equations.kernels(), cells,
                                           0.3, flux_history,
                                           accepted_flux),
                   "two committed flux revisions publish");
  if (!passed) {
    return false;
  }

  OwnedField rho_n = make_field(kDensity, cells, 0U, 601U);
  OwnedField rho_nm1 = make_field(kDensity, cells, 0U, 602U);
  OwnedField h_n = make_field(kEnthalpy, cells, 2U, 603U);
  OwnedField h_nm1 = make_field(kEnthalpy, cells, 2U, 604U);
  OwnedField y_n = make_field(kSpecies, cells, 2U, 605U);
  OwnedField y_nm1 = make_field(kSpecies, cells, 2U, 606U);
  OwnedField h_rhs_n = make_field(20U, cells, 0U, 607U);
  OwnedField h_rhs_nm1 = make_field(20U, cells, 0U, 608U);
  OwnedField y_rhs_n = make_field(21U, cells, 0U, 609U);
  OwnedField y_rhs_nm1 = make_field(21U, cells, 0U, 610U);
  OwnedField h_star = make_field(kEnthalpy, cells, 0U, 611U);
  OwnedField y_star = make_field(kSpecies, cells, 0U, 612U);
  OwnedField rho_work = make_field(30U, cells, 0U, 613U);
  OwnedField work_n = make_field(31U, cells, 0U, 614U);
  OwnedField work_nm1 = make_field(32U, cells, 0U, 615U);
  OwnedField low_rho = make_field(kDensity, cells, 0U, 616U);
  OwnedField low_h = make_field(kEnthalpy, cells, 1U, 617U);
  OwnedField low_y = make_field(kSpecies, cells, 1U, 618U);
  fill(rho_n, 1.0);
  fill(rho_nm1, 1.0);
  fill(h_n, 300000.0);
  fill(h_nm1, 299000.0);
  fill(y_n, 0.2);
  fill(y_nm1, 0.1);
  fill(h_rhs_n, 3000.0);
  fill(h_rhs_nm1, 1000.0);
  fill(y_rhs_n, 0.03);
  fill(y_rhs_nm1, 0.01);
  fill(h_star, -91.0);
  fill(y_star, -92.0);

  const std::array<ConstFieldView, 1U> species_n{as_const(y_n.view)};
  const std::array<ConstFieldView, 1U> species_nm1{as_const(y_nm1.view)};
  const std::array<PredictorRateHistory, 1U> species_rhs{{
      {as_const(y_rhs_n.view), as_const(y_rhs_nm1.view)},
  }};
  std::array<FieldView, 1U> species_output{y_star.view};
  ThermophysicalPredictorInput input;
  input.dt = 0.1;
  input.bdf = {15.0, -20.0, 5.0, 2U};
  input.time = 701U;
  input.geometry = fixture.geometry.topology_revision();
  input.boundary = fixture.boundary.revision();
  input.transport = fixture.transport.fingerprint();
  input.density_accepted = as_const(rho_n.view);
  input.density_previous = as_const(rho_nm1.view);
  input.enthalpy_accepted = as_const(h_n.view);
  input.enthalpy_previous = as_const(h_nm1.view);
  input.species_accepted = {species_n.data(), species_n.size()};
  input.species_previous = {species_nm1.data(), species_nm1.size()};
  input.mass_flux_accepted = accepted_flux;
  input.mass_flux_previous = previous_flux;
  input.enthalpy_nonadvective_rhs = {as_const(h_rhs_n.view),
                                     as_const(h_rhs_nm1.view)};
  input.species_nonadvective_rhs = {species_rhs.data(), species_rhs.size()};
  std::array<ThermophysicalGhostHistory, 1U> ghost_history{};
  attach_ghost_authority(fixture, input, ghost_history);
  ThermophysicalPredictorOutput output{
      h_star.view, {species_output.data(), species_output.size()},
      rho_work.view, work_n.view, work_nm1.view, {}, low_rho.view, low_h.view,
      {&low_y.view, 1U}, {}, paired_flux};
  ThermophysicalPredictorDiagnostics diagnostics;
  ThermophysicalPredictorCertificate certificate;
  passed &= expect(static_cast<bool>(
                       fixture.equations.thermophysical_predictor().predict(
                           MPI_COMM_SELF, Status{}, input, output, {},
                           diagnostics, certificate)) &&
                       certificate.valid(),
                   "BDF2 predictor consumes committed flux histories");
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        passed &= expect(close(h_star.view.unchecked({x, y, z}, 0U),
                               4510000.0 / 15.0) &&
                             close(y_star.view.unchecked({x, y, z}, 0U),
                                   3.55 / 15.0),
                         "BDF2 conservative predictor matches oracle");
      }
    }
  }

  const std::vector<double> h_snapshot = h_star.bytes;
  const std::vector<double> y_snapshot = y_star.bytes;
  const ThermophysicalPredictorCertificate marker = certificate;
  ConstFaceFluxView stale = accepted_flux;
  ++stale.revision;
  input.mass_flux_accepted = stale;
  passed &= expect(
      fixture.equations.thermophysical_predictor()
                  .predict(MPI_COMM_SELF, Status{}, input, output, {},
                           diagnostics, certificate)
                  .code == StatusCode::invalid_plan &&
          h_star.bytes == h_snapshot && y_star.bytes == y_snapshot &&
          same_certificate(certificate, marker),
      "stale committed-flux identity rejects before trial output writes");

  input.mass_flux_accepted = accepted_flux;
  ThermophysicalPredictorOutput aliased = output;
  aliased.enthalpy = h_n.view;
  passed &= expect(
      fixture.equations.thermophysical_predictor()
              .predict(MPI_COMM_SELF, Status{}, input, aliased, {}, diagnostics,
                       certificate)
              .code == StatusCode::invalid_plan &&
          same_certificate(certificate, marker),
      "predictor rejects accepted-history/output alias");

  input.dt = 0.2;
  passed &= expect(fixture.equations.thermophysical_predictor()
                           .predict(MPI_COMM_SELF, Status{}, input, output, {},
                                    diagnostics, certificate)
                           .code == StatusCode::invalid_plan &&
                       same_certificate(certificate, marker),
                   "BDF coefficients inconsistent with dt reject");
  input.dt = 0.1;

  fill(h_rhs_n, std::numeric_limits<double>::quiet_NaN());
  passed &= expect(fixture.equations.thermophysical_predictor()
                           .predict(MPI_COMM_SELF, Status{}, input, output, {},
                                    diagnostics, certificate)
                           .code == StatusCode::numerical_failure &&
                       h_star.bytes == h_snapshot &&
                       y_star.bytes == y_snapshot &&
                       same_certificate(certificate, marker),
                   "non-finite registered RHS rejects before outputs");
  fill(h_rhs_n, 3000.0);

  ConstFieldView revised_density = as_const(rho_n.view);
  ++revised_density.revision;
  input.density_accepted = revised_density;
  ThermophysicalPredictorCertificate revised_certificate;
  passed &= expect(static_cast<bool>(
                       fixture.equations.thermophysical_predictor().predict(
                           MPI_COMM_SELF, Status{}, input, output, {},
                           diagnostics, revised_certificate)) &&
                       revised_certificate.valid() &&
                       revised_certificate.state != marker.state,
                   "density-history revision changes predictor certificate");

  ClosedMassPlan closed_mass;
  const ThermophysicalSpec thermophysics = thermophysical_spec();
  passed &= expect(static_cast<bool>(ClosedMassPlan::compile(
                       PressureReferenceKind::closed_mass, thermophysics,
                       closed_mass)),
                   "closed-mass plan compiles");
  const std::size_t cell_count = static_cast<std::size_t>(cells.x) * cells.y *
                                 static_cast<std::size_t>(cells.z);
  std::vector<double> pi(cell_count, 0.0);
  std::vector<double> volume(cell_count, 1.0 / cell_count);
  std::vector<std::uint8_t> active(cell_count, 1U);
  const std::array<double, 1U> composition{
      y_star.view.unchecked({0, 0, 0}, 0U)};
  ThermoState thermo_state;
  constexpr double p_ref = 101325.0;
  passed &= expect(static_cast<bool>(fixture.thermodynamics.evaluate(
                       p_ref, h_star.view.unchecked({0, 0, 0}, 0U),
                       {composition.data(), composition.size()}, {},
                       thermo_state)),
                   "predicted h/Y evaluate through production EOS");
  const ClosedMassCellView closed_cells{
      {pi.data(), pi.size()},
      {h_star.bytes.data(), h_star.bytes.size()},
      {y_star.bytes.data(), y_star.bytes.size()},
      {volume.data(), volume.size()},
      {active.data(), active.size()},
      revised_certificate.state,
  };
  ClosedMassResult closed_result;
  PressureReferenceCertificate pressure_certificate;
  passed &= expect(static_cast<bool>(
                       fixture.equations.pressure_reference().solve_closed_mass(
                           MPI_COMM_SELF, closed_mass, fixture.thermodynamics,
                           revised_certificate, 1U, closed_cells,
                           thermo_state.rho, p_ref, closed_result,
                           pressure_certificate)) &&
                       pressure_certificate.valid() &&
                       close(closed_result.pressure_reference, p_ref) &&
                       close(closed_result.residual, 0.0),
                   "predictor certificate orders closed-mass p_ref solve");
  const ClosedMassResult accepted_result = closed_result;
  const PressureReferenceCertificate accepted_pressure_certificate =
      pressure_certificate;
  passed &= expect(
      fixture.equations.pressure_reference()
                  .solve_closed_mass(
                      MPI_COMM_SELF, closed_mass, fixture.thermodynamics,
                      revised_certificate, 2U, closed_cells, thermo_state.rho,
                      p_ref, closed_result, pressure_certificate)
                  .code == StatusCode::invalid_plan &&
          close(closed_result.pressure_reference,
                accepted_result.pressure_reference) &&
          pressure_certificate.pressure_reference ==
              accepted_pressure_certificate.pressure_reference,
      "wrong coupling stage cannot publish p_ref or its certificate");
  ClosedMassCellView stale_closed_cells = closed_cells;
  ++stale_closed_cells.predictor_state;
  passed &= expect(
      fixture.equations.pressure_reference()
                  .solve_closed_mass(
                      MPI_COMM_SELF, closed_mass, fixture.thermodynamics,
                      revised_certificate, 1U, stale_closed_cells,
                      thermo_state.rho, p_ref, closed_result,
                      pressure_certificate)
                  .code == StatusCode::invalid_plan &&
          pressure_certificate.pressure_reference ==
              accepted_pressure_certificate.pressure_reference,
      "closed-mass view/predictor identity mutation rejects atomically");

  Fixture open_fixture;
  passed &= expect(make_fixture(n, PressureReferenceKind::boundary_absolute,
                                open_fixture),
                   "open pressure-reference fixture compiles");
  input.geometry = open_fixture.geometry.topology_revision();
  input.boundary = open_fixture.boundary.revision();
  input.transport = open_fixture.transport.fingerprint();
  attach_ghost_authority(open_fixture, input, ghost_history);
  ThermophysicalPredictorCertificate open_predictor;
  passed &= expect(static_cast<bool>(
                       open_fixture.equations.thermophysical_predictor().predict(
                           MPI_COMM_SELF, Status{}, input, output, {},
                           diagnostics, open_predictor)) &&
                       open_predictor.valid(),
                   "open case preserves predictor-first ordering");
  PressureReferenceCertificate open_certificate;
  passed &= expect(static_cast<bool>(
                       open_fixture.equations.pressure_reference()
                           .certify_open_boundary(101325.0, 901U,
                                                  open_predictor,
                                                  open_certificate)) &&
                       open_certificate.valid() &&
                       open_certificate.kind ==
                           PressureReferenceKind::boundary_absolute,
                   "open absolute-pressure boundary publishes authority");
  const PressureReferenceCertificate open_marker = open_certificate;
  passed &= expect(
      open_fixture.equations.pressure_reference()
                  .certify_open_boundary(101325.0, 0U, open_predictor,
                                         open_certificate)
                  .code == StatusCode::invalid_plan &&
          open_certificate.pressure_reference ==
              open_marker.pressure_reference,
      "open boundary revision mutation rejects atomically");

  const std::array<ConstFieldView, 1U> empty_species_previous{{{}}};
  const std::array<PredictorRateHistory, 1U> be_species_rhs{{
      {as_const(y_rhs_n.view), {}},
  }};
  input.bdf = {10.0, -10.0, 0.0, 1U};
  input.density_previous = {};
  input.enthalpy_previous = {};
  input.species_previous = {empty_species_previous.data(),
                            empty_species_previous.size()};
  input.mass_flux_previous = {};
  input.enthalpy_nonadvective_rhs = {as_const(h_rhs_n.view), {}};
  input.species_nonadvective_rhs = {be_species_rhs.data(),
                                    be_species_rhs.size()};
  attach_ghost_authority(open_fixture, input, ghost_history);
  ThermophysicalPredictorCertificate be_certificate;
  passed &= expect(static_cast<bool>(
                       open_fixture.equations.thermophysical_predictor().predict(
                           MPI_COMM_SELF, Status{}, input, output, {},
                           diagnostics, be_certificate)) &&
                       be_certificate.valid() && be_certificate.order == 1U &&
                       be_certificate.previous_face_flux == 0U,
                   "BE startup consumes no previous state or flux authority");
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        passed &= expect(
            close(h_star.view.unchecked({x, y, z}, 0U), 300300.0) &&
                close(y_star.view.unchecked({x, y, z}, 0U), 0.203),
            "BE startup predictor matches accepted-only oracle");
      }
    }
  }
  return passed;
}

bool test_nonadvective_rate_path() {
  Fixture fixture;
  constexpr std::int32_t n = 4;
  bool passed = expect(make_fixture(n, PressureReferenceKind::closed_mass,
                                    fixture),
                       "nonadvective-rate fixture compiles");
  if (!passed) return false;
  const Int3 cells = fixture.patch.cells;
  const auto make_components = [&](FieldId field, std::uint8_t components,
                                   std::uint8_t ghosts,
                                   RevisionToken revision) {
    OwnedField value = make_field(field, cells, ghosts, revision);
    const std::size_t scalar_size = value.bytes.size();
    value.bytes.resize(scalar_size * components, 0.0);
    value.view.base = value.bytes.data() + ghosts +
                      ghosts * value.view.stride_y +
                      ghosts * value.view.stride_z;
    value.view.components = components;
    value.view.component_stride = scalar_size;
    return value;
  };

  OwnedField rho = make_field(kDensity, cells, 2U, 8001U);
  OwnedField rho_n = make_field(kDensity, cells, 2U, 8002U);
  OwnedField rho_nm1 = make_field(kDensity, cells, 2U, 8003U);
  OwnedField velocity = make_components(kVelocity, 3U, 2U, 8004U);
  OwnedField velocity_n = make_components(kVelocity, 3U, 2U, 8005U);
  OwnedField velocity_nm1 = make_components(kVelocity, 3U, 2U, 8006U);
  OwnedField pressure = make_field(kPressure, cells, 2U, 8007U);
  OwnedField pressure_n = make_field(kPressure, cells, 2U, 8008U);
  OwnedField pressure_nm1 = make_field(kPressure, cells, 2U, 8009U);
  OwnedField enthalpy = make_field(kEnthalpy, cells, 2U, 8010U);
  OwnedField enthalpy_n = make_field(kEnthalpy, cells, 2U, 8011U);
  OwnedField enthalpy_nm1 = make_field(kEnthalpy, cells, 2U, 8012U);
  OwnedField temperature = make_field(kTemperature, cells, 2U, 8013U);
  OwnedField temperature_n = make_field(kTemperature, cells, 2U, 8014U);
  OwnedField temperature_nm1 = make_field(kTemperature, cells, 2U, 8015U);
  OwnedField species = make_field(kSpecies, cells, 2U, 8016U);
  OwnedField species_n = make_field(kSpecies, cells, 2U, 8017U);
  OwnedField species_nm1 = make_field(kSpecies, cells, 2U, 8018U);
  OwnedField molecular = make_field(60U, cells, 1U, 8019U);
  OwnedField effective = make_field(kEffectiveViscosity, cells, 1U, 8020U);
  OwnedField conductivity = make_field(61U, cells, 1U, 8021U);
  OwnedField gradient = make_components(kVelocityGradient, 9U, 0U, 8022U);
  OwnedField enthalpy_rate = make_field(62U, cells, 0U, 8023U);
  OwnedField species_rate = make_field(63U, cells, 0U, 8024U);
  OwnedField diffusion_scratch = make_field(64U, cells, 0U, 8025U);
  OwnedField scalar_diffusivity = make_field(65U, cells, 1U, 8026U);
  for (OwnedField* field : {&rho, &rho_n, &rho_nm1}) fill(*field, 1.0);
  for (OwnedField* field : {&velocity, &velocity_n, &velocity_nm1,
                            &pressure, &pressure_n, &pressure_nm1,
                            &gradient})
    fill(*field, 0.0);
  for (OwnedField* field : {&enthalpy, &enthalpy_n, &enthalpy_nm1})
    fill(*field, 300000.0);
  for (OwnedField* field : {&temperature, &temperature_n, &temperature_nm1})
    fill(*field, 300.0);
  for (OwnedField* field : {&species, &species_n, &species_nm1})
    fill(*field, 0.2);
  fill(molecular, 1.8e-5);
  fill(effective, 2.4e-5);
  fill(conductivity, 0.026);
  fill(enthalpy_rate, -7.0);
  fill(species_rate, -8.0);

  const std::array<PrimitiveHistory, 1U> species_history{{
      {as_const(species.view), as_const(species_n.view),
       as_const(species_nm1.view)},
  }};
  EquationStateView state;
  state.density = {as_const(rho.view), as_const(rho_n.view),
                   as_const(rho_nm1.view)};
  state.velocity = {as_const(velocity.view), as_const(velocity_n.view),
                    as_const(velocity_nm1.view)};
  state.pressure_perturbation = {
      as_const(pressure.view), as_const(pressure_n.view),
      as_const(pressure_nm1.view)};
  state.enthalpy = {as_const(enthalpy.view), as_const(enthalpy_n.view),
                    as_const(enthalpy_nm1.view)};
  state.temperature = {
      as_const(temperature.view), as_const(temperature_n.view),
      as_const(temperature_nm1.view)};
  state.independent_species = {species_history.data(), species_history.size()};
  state.pressure_reference = 101325.0;
  state.accepted_pressure_reference = 101325.0;
  state.previous_pressure_reference = 101325.0;
  EquationMaterialView material;
  material.molecular_viscosity = as_const(molecular.view);
  material.effective_viscosity = as_const(effective.view);
  material.thermal_conductivity = as_const(conductivity.view);
  std::array<FieldView, 1U> species_rates{species_rate.view};
  ThermophysicalRateInput input{state, material, as_const(gradient.view),
                                {10.0, -10.0, 0.0, 1U}, 1U, {}};
  ThermophysicalRateOutput output{
      enthalpy_rate.view, {species_rates.data(), species_rates.size()}, {},
      diffusion_scratch.view, scalar_diffusivity.view};
  ThermophysicalRateCertificate certificate;
  Status status = evaluate_thermophysical_rates(fixture.equations, input,
                                                output, certificate);
  passed &= expect(status && certificate.valid(),
                   "production nonadvective-rate path evaluates");
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        passed &= expect(
            close(enthalpy_rate.view.unchecked({x, y, z}, 0U), 0.0) &&
                close(species_rate.view.unchecked({x, y, z}, 0U), 0.0),
            "constant state has zero nonadvective enthalpy/species rate");

  // The target-time pressure derivative belongs to the coupled
  // BDF(rho*h-p) residual.  Persisting it here would make the next predictor
  // EX2 a derivative history and recreate the step-count feedback loop.
  fill(pressure, 2.0);
  fill(pressure_n, 1.0);
  fill(pressure_nm1, 0.0);
  ++pressure.view.revision;
  ++pressure_n.view.revision;
  ++pressure_nm1.view.revision;
  input.state.pressure_perturbation = {
      as_const(pressure.view), as_const(pressure_n.view),
      as_const(pressure_nm1.view)};
  ThermophysicalRateCertificate temporal_pressure_certificate;
  status = evaluate_thermophysical_rates(
      fixture.equations, input, output, temporal_pressure_certificate);
  passed &= expect(status && temporal_pressure_certificate.valid() &&
                       close(enthalpy_rate.view.unchecked({1, 1, 1}, 0U),
                             0.0),
                   "stored spatial energy rate excludes BDF(p)");

  // Excluding BDF(p) must not remove the physical spatial pressure work.
  // On the uniform unit cube p=x/n has grad(p)=(1,0,0), so U=(1/4,0,0)
  // must leave U.grad(p)=1/4 in the persisted spatial rate.
  fill(velocity, 0.0);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        velocity.view.unchecked({x, y, z}, 0U) = 0.25;
  for (OwnedField* field : {&pressure, &pressure_n, &pressure_nm1}) {
    for (std::int32_t z = -2; z < cells.z + 2; ++z)
      for (std::int32_t y = -2; y < cells.y + 2; ++y)
        for (std::int32_t x = -2; x < cells.x + 2; ++x)
          field->view.unchecked({x, y, z}, 0U) =
              static_cast<double>(x) / static_cast<double>(n);
  }
  ++velocity.view.revision;
  ++pressure.view.revision;
  ++pressure_n.view.revision;
  ++pressure_nm1.view.revision;
  input.state.velocity.trial = as_const(velocity.view);
  input.state.pressure_perturbation = {
      as_const(pressure.view), as_const(pressure_n.view),
      as_const(pressure_nm1.view)};
  ThermophysicalRateCertificate spatial_pressure_certificate;
  status = evaluate_thermophysical_rates(
      fixture.equations, input, output, spatial_pressure_certificate);
  passed &= expect(status && spatial_pressure_certificate.valid() &&
                       close(enthalpy_rate.view.unchecked({1, 1, 1}, 0U),
                             0.25),
                   "stored spatial energy rate retains U dot grad(p)");

  fill(velocity, 0.0);
  fill(pressure, 0.0);
  fill(pressure_n, 0.0);
  fill(pressure_nm1, 0.0);
  ++velocity.view.revision;
  ++pressure.view.revision;
  ++pressure_n.view.revision;
  ++pressure_nm1.view.revision;
  input.state.velocity.trial = as_const(velocity.view);
  input.state.pressure_perturbation = {
      as_const(pressure.view), as_const(pressure_n.view),
      as_const(pressure_nm1.view)};

  temperature.view.unchecked({1, 1, 1}, 0U) = 310.0;
  ++temperature.view.revision;
  input.state.temperature.trial = as_const(temperature.view);
  ThermophysicalRateCertificate conduction_certificate;
  status = evaluate_thermophysical_rates(fixture.equations, input, output,
                                         conduction_certificate);
  passed &= expect(status && conduction_certificate.valid() &&
                       conduction_certificate.state != certificate.state &&
                       !close(enthalpy_rate.view.unchecked({1, 1, 1}, 0U),
                              0.0),
                   "temperature curvature enters the production rate path");

  ThermophysicalRateOutput aliased = output;
  aliased.diffusion_scratch = enthalpy_rate.view;
  const ThermophysicalRateCertificate marker = conduction_certificate;
  passed &= expect(
      evaluate_thermophysical_rates(fixture.equations, input, aliased,
                                    conduction_certificate)
              .code == StatusCode::invalid_plan &&
          conduction_certificate.equations == marker.equations &&
          conduction_certificate.state == marker.state &&
          conduction_certificate.rates == marker.rates,
      "rate workspace/output alias rejects without replacing certificate");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  const bool passed = test_bdf2_predictor_and_mutations();
  const bool mms_passed = test_predictor_mms_orders();
  const bool rates_passed = test_nonadvective_rate_path();
  MPI_Finalize();
  return passed && mms_passed && rates_passed ? 0 : 1;
}
