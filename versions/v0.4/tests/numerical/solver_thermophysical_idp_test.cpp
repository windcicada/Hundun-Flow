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
constexpr FieldId kPassive = 6U;
constexpr FieldId kEffectiveViscosity = 7U;
constexpr FieldId kCompressibility = 8U;
constexpr FieldId kVelocityGradient = 9U;
constexpr double kDt = 0.1;
constexpr BdfCoefficients kBdf{15.0, -20.0, 5.0, 2U};
constexpr double kSharedFaceDt = 0.006;
constexpr double kSharedFaceRho = 0.56445295871161083;
constexpr double kSharedFaceDivergence = 107.47787969771555;
constexpr BdfCoefficients kSharedFaceBdf{
    1.5 / kSharedFaceDt, -2.0 / kSharedFaceDt, 0.5 / kSharedFaceDt, 2U};

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected, double relative = 2.0e-12) {
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::abs(actual - expected) <=
             relative * std::max(1.0, std::abs(expected));
}

bool same_bits(double left, double right) noexcept {
  return std::memcmp(&left, &right, sizeof(left)) == 0;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField make_field(FieldId field, Int3 cells, std::uint8_t ghosts,
                      RevisionToken revision, StorageIdentity storage_id) {
  OwnedField result;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.storage.assign(nx * ny * nz, 0.0);
  result.view.base = result.storage.data() + ghosts + ghosts * nx +
                     ghosts * nx * ny;
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = 1U;
  result.view.stride_y = nx;
  result.view.stride_z = nx * ny;
  result.view.component_stride = nx * ny * nz;
  result.view.field = field;
  result.view.revision = revision;
  result.view.storage_identity = storage_id;
  result.view.revision_domain = UINT64_C(0x1d7001);
  return result;
}

OwnedField make_field_components(FieldId field, Int3 cells,
                                 std::uint8_t components,
                                 std::uint8_t ghosts,
                                 RevisionToken revision,
                                 StorageIdentity storage_id) {
  OwnedField result =
      make_field(field, cells, ghosts, revision, storage_id);
  const std::size_t scalar_size = result.storage.size();
  result.storage.resize(static_cast<std::size_t>(components) * scalar_size,
                        0.0);
  result.view.base =
      result.storage.data() + ghosts + ghosts * result.view.stride_y +
      ghosts * result.view.stride_z;
  result.view.components = components;
  result.view.component_stride = scalar_size;
  return result;
}

void fill(OwnedField& field, double value) {
  std::fill(field.storage.begin(), field.storage.end(), value);
}

CartesianMeshSpec mesh_spec(Int3 cells) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = cells;
  mesh.minimum_spacing = {
      1.0 / static_cast<double>(cells.x),
      1.0 / static_cast<double>(cells.y),
      1.0 / static_cast<double>(cells.z)};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells =
      static_cast<std::uint64_t>(cells.x) *
      static_cast<std::uint64_t>(cells.y) *
      static_cast<std::uint64_t>(cells.z);
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
  spec.data_file = "analytic-idp.d";
  spec.minimum_temperature = 200.0;
  spec.maximum_temperature = 3200.0;
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
  LinearWorkspaceRequirements endpoint_requirements;
  OwnedField endpoint_vectors;
  OwnedField endpoint_scalars;
  HaloEngine endpoint_halo;
  HaloEngine donor_halo;
  FaceFluxStorage paired_flux;
  ReductionEngine endpoint_reductions;
  std::vector<std::uint8_t> endpoint_cells;
  std::vector<std::uint8_t> endpoint_x_faces;
  std::vector<std::uint8_t> endpoint_y_faces;
  std::vector<std::uint8_t> endpoint_z_faces;
  ConservativeEnthalpyEndpoint endpoint;
};

bool make_fixture(Fixture& out, MPI_Comm communicator = MPI_COMM_SELF,
                  Int3 global_cells = {2, 2, 2},
                  bool isolate_internal_face = false,
                  bool physical_x = false) {
  const CartesianMeshSpec mesh = mesh_spec(global_cells);
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = UINT64_C(0x1d7001);
  model.pressure_reference = PressureReferenceKind::closed_mass;
  model.transported_scalars = {
      {"species_a", TransportedScalarRole::species, 1.0, 1.0},
      {"tracer", TransportedScalarRole::passive_scalar, 1.0, 1.0}};
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::periodic;
    face.thermal_kind = BoundaryKind::none;
  }
  if (physical_x) {
    for (std::size_t face = 0U; face < 2U; ++face) {
      model.boundaries[face].flow_kind = BoundaryKind::no_slip_wall;
      model.boundaries[face].thermal_kind = BoundaryKind::adiabatic_wall;
      model.boundaries[face].scalars.push_back(
          {"species_a", ScalarBoundaryKind::zero_gradient, 0.0,
           ScalarBoundaryKind::zero_gradient, 0.0});
      model.boundaries[face].scalars.push_back(
          {"tracer", ScalarBoundaryKind::zero_gradient, 0.0,
           ScalarBoundaryKind::zero_gradient, 0.0});
    }
  }
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::limited_central2;
  model.schemes.species = ConvectionScheme::tvd2;
  model.schemes.passive_scalar = ConvectionScheme::tvd2;

  FieldRegistry registry;
  FieldId id = 0U;
  if (!registry.require_field("rho", 1U, 2U, id) || id != kDensity ||
      !registry.require_field("U", 3U, 2U, id) || id != kVelocity ||
      !registry.require_field("pi", 1U, 2U, id) || id != kPressure ||
      !registry.require_field("h", 1U, 2U, id) || id != kEnthalpy ||
      !registry.require_field("T", 1U, 2U, id) || id != kTemperature ||
      !CartesianGeometryCompiler::compile(communicator, mesh, {},
                                          out.geometry, out.patch) ||
      !BoundaryCompiler::compile(communicator, model, out.geometry,
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
  const std::array<FieldId, 10U> declared{
      kDensity, kVelocity, kPressure, kEnthalpy, kTemperature, kSpecies,
      kPassive, kEffectiveViscosity, kCompressibility, kVelocityGradient};
  if (!out.contributions.configure({declared.data(), declared.size()}) ||
      !out.contributions.freeze()) {
    return false;
  }
  const std::array<ScalarEquationSpec, 2U> scalars{{
      {kSpecies, TransportedScalarRole::species, 1.0, 1.0},
      {kPassive, TransportedScalarRole::passive_scalar, 1.0, 1.0},
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
  spec.pressure_reference = PressureReferenceKind::closed_mass;
  spec.scalars = {scalars.data(), scalars.size()};
  spec.closed_mass_service_stage = 1U;
  spec.maximum_cells_per_rank =
      static_cast<std::size_t>(global_cells.x) *
      static_cast<std::size_t>(global_cells.y) *
      static_cast<std::size_t>(global_cells.z);
  if (!EquationPlanSet::compile(
          communicator, out.schemes, out.geometry, out.patch, out.boundary,
          out.contributions, out.thermodynamics, out.transport, spec,
          out.equations) ||
      !make_linear_workspace_requirements(
          LinearAlgorithm::fgmres, out.patch.cells, 1U, 12U,
          ReductionMode::mpi_allreduce, 920U,
          out.endpoint_requirements)) {
    return false;
  }
  out.endpoint_vectors = make_field_components(
      90U, out.patch.cells, out.endpoint_requirements.vector_slots, 1U, 921U,
      922U);
  out.endpoint_scalars = make_field(
      91U,
      {static_cast<std::int32_t>(out.endpoint_requirements.scalar_doubles),
       1, 1},
      0U, 923U, 922U);
  const std::array<HaloFieldSpec, 1U> halo_fields{{
      {out.endpoint_vectors.view.field, 1U, 1U},
  }};
  const std::array<HaloFieldSpec, 1U> donor_halo_fields{{
      {kDensity, 1U, 1U},
  }};
  if (!out.endpoint_halo.reserve(
          communicator, out.patch,
          {halo_fields.data(), halo_fields.size()},
          out.boundary.halo_topology()) ||
      !out.donor_halo.reserve(
          communicator, out.patch,
          {donor_halo_fields.data(), donor_halo_fields.size()},
          out.boundary.halo_topology()) ||
      !FaceFluxStorage::allocate_workspace(out.patch.cells, 1U,
                                           out.paired_flux) ||
      !ReductionEngine::compile(
          communicator, ReductionMode::mpi_allreduce,
          out.endpoint_requirements.reduction_capacity,
          out.endpoint_reductions)) {
    return false;
  }
  MgDomainActivityView activity;
  if (isolate_internal_face) {
    const auto count = [](Int3 shape) {
      return static_cast<std::size_t>(shape.x) *
             static_cast<std::size_t>(shape.y) *
             static_cast<std::size_t>(shape.z);
    };
    const Int3 x_faces{out.patch.cells.x + 1, out.patch.cells.y,
                       out.patch.cells.z};
    const Int3 y_faces{out.patch.cells.x, out.patch.cells.y + 1,
                       out.patch.cells.z};
    const Int3 z_faces{out.patch.cells.x, out.patch.cells.y,
                       out.patch.cells.z + 1};
    out.endpoint_cells.assign(count(out.patch.cells), 1U);
    out.endpoint_x_faces.assign(count(x_faces), 1U);
    out.endpoint_y_faces.assign(count(y_faces), 1U);
    out.endpoint_z_faces.assign(count(z_faces), 1U);
    out.endpoint_x_faces[1U] = 0U;
    activity = {{out.endpoint_cells.data(), out.endpoint_cells.size()},
                {out.endpoint_x_faces.data(), out.endpoint_x_faces.size()},
                {out.endpoint_y_faces.data(), out.endpoint_y_faces.size()},
                {out.endpoint_z_faces.data(), out.endpoint_z_faces.size()},
                UINT64_C(0x1d7021), UINT64_C(0x1d7022)};
  }
  const ConservativeEnthalpyEndpointServices endpoint_services{
      communicator,
      &out.equations.kernels(),
      &out.boundary,
      out.patch,
      activity,
      &out.endpoint_halo,
      &out.endpoint_reductions,
      15U,
      out.endpoint_requirements,
      out.endpoint_vectors.view,
      out.endpoint_scalars.view,
      kEnthalpy};
  return static_cast<bool>(
      ConservativeEnthalpyEndpoint::bind(endpoint_services, out.endpoint));
}

struct FluxHistory {
  FieldId dependency{};
  StateLayers layers;
  AttemptTransaction transaction;
  FaceFluxStorage storage;
  FinalFaceFluxAuthority authority;
  FinalFaceFluxWriter writer;
};

bool make_flux_history(Int3 cells, FluxHistory& result) {
  FieldRegistry registry;
  FieldSchema schema;
  if (!registry.declare_field("idp.flux_dependency", 1U, 0U,
                              result.dependency) ||
      !registry.freeze(schema)) {
    return false;
  }
  const std::array requests{ArenaFieldRequest{
      result.dependency, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout layout;
  return static_cast<bool>(ArenaLayout::compile(
             schema, {requests.data(), requests.size()}, layout)) &&
         static_cast<bool>(StateLayers::allocate(layout, result.layers)) &&
         static_cast<bool>(AttemptTransaction::create(
             result.layers.field_count(), 1U, result.layers.field_count(),
             result.transaction)) &&
         static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                            result.storage)) &&
         static_cast<bool>(result.authority.claim(
             41U, 0U, result.transaction, result.writer));
}

bool commit_zero_flux(const CartesianKernelPlan& kernels, Int3 cells,
                      FluxHistory& history, ConstFaceFluxView& committed) {
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

  OwnedField rho = make_field(kDensity, cells, 2U, 501U, 601U);
  OwnedField velocity = make_field(kVelocity, cells, 2U, 502U, 602U);
  const std::size_t scalar_size = velocity.storage.size();
  velocity.storage.resize(3U * scalar_size, 0.0);
  velocity.view.base = velocity.storage.data() + 2U +
                       2U * velocity.view.stride_y +
                       2U * velocity.view.stride_z;
  velocity.view.components = 3U;
  velocity.view.component_stride = scalar_size;
  fill(rho, 1.0);
  const std::array<ConstFieldView, 2U> reads{as_const(rho.view),
                                             as_const(velocity.view)};
  const KernelInvocation invocation{
      {reads.data(), reads.size()}, {}, {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U,
      nullptr};
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

ThermophysicalGhostAuthority ghost_authority(ConstFieldView view,
                                              RevisionToken geometry,
                                              RevisionToken boundary,
                                              std::uintptr_t exchange_plan) {
  return {exchange_plan, view.field, view.revision, view.storage_identity,
          view.revision_domain, geometry, boundary,
          static_cast<std::uint8_t>(view.ghosts.x)};
}

bool same_certificate(const ThermophysicalPredictorCertificate& left,
                      const ThermophysicalPredictorCertificate& right) {
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
         left.predicted_density == right.predicted_density &&
         left.predicted_density_storage ==
             right.predicted_density_storage &&
         left.predicted_density_revision_domain ==
             right.predicted_density_revision_domain &&
         left.paired_face_flux == right.paired_face_flux &&
         left.paired_face_flux_storage == right.paired_face_flux_storage &&
         left.paired_face_flux_revision_domain ==
             right.paired_face_flux_revision_domain &&
         left.state == right.state && left.order == right.order;
}

struct PredictorData {
  OwnedField rho;
  OwnedField rho_previous;
  OwnedField enthalpy;
  OwnedField enthalpy_previous;
  OwnedField species;
  OwnedField species_previous;
  OwnedField enthalpy_rhs;
  OwnedField enthalpy_rhs_previous;
  OwnedField species_rhs;
  OwnedField species_rhs_previous;
  OwnedField passive;
  OwnedField passive_previous;
  OwnedField passive_rhs;
  OwnedField passive_rhs_previous;
  OwnedField high_enthalpy;
  OwnedField high_species;
  OwnedField high_passive;
  OwnedField density_workspace;
  OwnedField accepted_advection_workspace;
  OwnedField previous_advection_workspace;
  OwnedField low_density_workspace;
  OwnedField low_enthalpy_workspace;
  OwnedField low_species_workspace;
  OwnedField low_passive_workspace;

  explicit PredictorData(Int3 cells, double accepted_h, double previous_h,
                         double accepted_rhs,
                         double previous_rhs = 0.0)
      : rho(make_field(kDensity, cells, 0U, 601U, 701U)),
        rho_previous(make_field(kDensity, cells, 0U, 602U, 702U)),
        enthalpy(make_field(kEnthalpy, cells, 2U, 603U, 703U)),
        enthalpy_previous(make_field(kEnthalpy, cells, 2U, 604U, 704U)),
        species(make_field(kSpecies, cells, 2U, 605U, 705U)),
        species_previous(make_field(kSpecies, cells, 2U, 606U, 706U)),
        enthalpy_rhs(make_field(20U, cells, 0U, 607U, 707U)),
        enthalpy_rhs_previous(make_field(20U, cells, 0U, 608U, 708U)),
        species_rhs(make_field(21U, cells, 0U, 609U, 709U)),
        species_rhs_previous(make_field(21U, cells, 0U, 610U, 710U)),
        passive(make_field(kPassive, cells, 2U, 619U, 719U)),
        passive_previous(make_field(kPassive, cells, 2U, 620U, 720U)),
        passive_rhs(make_field(22U, cells, 0U, 621U, 721U)),
        passive_rhs_previous(make_field(22U, cells, 0U, 622U, 722U)),
        high_enthalpy(make_field(kEnthalpy, cells, 0U, 611U, 711U)),
        high_species(make_field(kSpecies, cells, 0U, 612U, 712U)),
        high_passive(make_field(kPassive, cells, 0U, 623U, 723U)),
        density_workspace(make_field(30U, cells, 0U, 613U, 713U)),
        accepted_advection_workspace(
            make_field(31U, cells, 0U, 614U, 714U)),
        previous_advection_workspace(
            make_field(32U, cells, 0U, 615U, 715U)),
        low_density_workspace(make_field(kDensity, cells, 1U, 616U, 716U)),
        low_enthalpy_workspace(make_field(kEnthalpy, cells, 1U, 617U, 717U)),
        low_species_workspace(make_field(kSpecies, cells, 1U, 618U, 718U)),
        low_passive_workspace(make_field(kPassive, cells, 1U, 624U, 724U)) {
    fill(rho, 1.0);
    fill(rho_previous, 1.0);
    fill(enthalpy, accepted_h);
    fill(enthalpy_previous, previous_h);
    fill(species, 0.2);
    fill(species_previous, 0.2);
    fill(enthalpy_rhs, accepted_rhs);
    fill(enthalpy_rhs_previous, previous_rhs);
    fill(species_rhs, 0.0);
    fill(species_rhs_previous, 0.0);
    fill(passive, 0.4);
    fill(passive_previous, 0.4);
    fill(passive_rhs, 0.0);
    fill(passive_rhs_previous, 0.0);
    fill(high_enthalpy, -1.0);
    fill(high_species, -2.0);
    fill(high_passive, -2.5);
    fill(density_workspace, -3.0);
    fill(accepted_advection_workspace, -4.0);
    fill(previous_advection_workspace, -5.0);
    fill(low_density_workspace, -6.0);
    fill(low_enthalpy_workspace, -7.0);
    fill(low_species_workspace, -8.0);
    fill(low_passive_workspace, -9.0);
  }
};

struct PredictorCall {
  std::array<ConstFieldView, 1U> species_accepted;
  std::array<ConstFieldView, 1U> species_previous;
  std::array<PredictorRateHistory, 1U> species_rates;
  std::array<ThermophysicalGhostHistory, 1U> species_ghosts;
  std::array<FieldView, 1U> species_output;
  std::array<FieldView, 1U> low_species_output;
  std::array<ConstFieldView, 1U> passive_accepted;
  std::array<ConstFieldView, 1U> passive_previous;
  std::array<PredictorRateHistory, 1U> passive_rates;
  std::array<ThermophysicalGhostHistory, 1U> passive_ghosts;
  std::array<FieldView, 1U> passive_output;
  std::array<FieldView, 1U> low_passive_output;
  std::array<FieldView, 3U> slow_exchange_fields;
  ThermophysicalPredictorInput input{};
  ThermophysicalPredictorOutput output{};
  ThermophysicalPredictorSlowPath slow_path{};
  ThermophysicalPredictorDiagnostics diagnostics{};
  ThermophysicalPredictorCertificate certificate{};
};

PredictorCall make_call(Fixture& fixture, PredictorData& data,
                        ConstFaceFluxView accepted_flux,
                        ConstFaceFluxView previous_flux,
                        std::uintptr_t exchange_plan,
                        RevisionToken time) {
  PredictorCall call;
  call.species_accepted = {as_const(data.species.view)};
  call.species_previous = {as_const(data.species_previous.view)};
  call.species_rates[0U] = {as_const(data.species_rhs.view),
                            as_const(data.species_rhs_previous.view)};
  call.species_ghosts[0U].accepted =
      ghost_authority(as_const(data.species.view),
                      fixture.geometry.topology_revision(),
                      fixture.boundary.revision(), exchange_plan);
  call.species_ghosts[0U].previous =
      ghost_authority(as_const(data.species_previous.view),
                      fixture.geometry.topology_revision(),
                      fixture.boundary.revision(), exchange_plan);
  call.species_output = {data.high_species.view};
  call.low_species_output = {data.low_species_workspace.view};
  call.passive_accepted = {as_const(data.passive.view)};
  call.passive_previous = {as_const(data.passive_previous.view)};
  call.passive_rates[0U] = {as_const(data.passive_rhs.view),
                            as_const(data.passive_rhs_previous.view)};
  call.passive_ghosts[0U].accepted =
      ghost_authority(as_const(data.passive.view),
                      fixture.geometry.topology_revision(),
                      fixture.boundary.revision(), exchange_plan);
  call.passive_ghosts[0U].previous =
      ghost_authority(as_const(data.passive_previous.view),
                      fixture.geometry.topology_revision(),
                      fixture.boundary.revision(), exchange_plan);
  call.passive_output = {data.high_passive.view};
  call.low_passive_output = {data.low_passive_workspace.view};
  // The slow path exchanges the cold low-order bundle in place.  High-order
  // outputs are immutable candidate results and must never be used as halo
  // workspace when the IDP activates substeps.
  call.slow_exchange_fields = {data.low_enthalpy_workspace.view,
                               data.low_species_workspace.view,
                               data.low_passive_workspace.view};

  call.input.dt = kDt;
  call.input.bdf = kBdf;
  call.input.time = time;
  call.input.geometry = fixture.geometry.topology_revision();
  call.input.boundary = fixture.boundary.revision();
  call.input.transport = fixture.transport.fingerprint();
  call.input.density_accepted = as_const(data.rho.view);
  call.input.density_previous = as_const(data.rho_previous.view);
  call.input.enthalpy_accepted = as_const(data.enthalpy.view);
  call.input.enthalpy_previous = as_const(data.enthalpy_previous.view);
  call.input.species_accepted = {call.species_accepted.data(),
                                 call.species_accepted.size()};
  call.input.species_previous = {call.species_previous.data(),
                                 call.species_previous.size()};
  call.input.passive_scalars_accepted = {call.passive_accepted.data(),
                                         call.passive_accepted.size()};
  call.input.passive_scalars_previous = {call.passive_previous.data(),
                                         call.passive_previous.size()};
  call.input.mass_flux_accepted = accepted_flux;
  call.input.mass_flux_previous = previous_flux;
  call.input.enthalpy_nonadvective_rhs = {
      as_const(data.enthalpy_rhs.view),
      as_const(data.enthalpy_rhs_previous.view)};
  call.input.species_nonadvective_rhs = {call.species_rates.data(),
                                         call.species_rates.size()};
  call.input.passive_scalar_nonadvective_rhs = {call.passive_rates.data(),
                                                 call.passive_rates.size()};
  call.input.enthalpy_ghosts = {
      ghost_authority(as_const(data.enthalpy.view),
                      fixture.geometry.topology_revision(),
                      fixture.boundary.revision(), exchange_plan),
      ghost_authority(as_const(data.enthalpy_previous.view),
                      fixture.geometry.topology_revision(),
                      fixture.boundary.revision(), exchange_plan)};
  call.input.species_ghosts = {call.species_ghosts.data(),
                               call.species_ghosts.size()};
  call.input.passive_scalar_ghosts = {call.passive_ghosts.data(),
                                      call.passive_ghosts.size()};

  call.output = {
      data.high_enthalpy.view,
      {call.species_output.data(), call.species_output.size()},
      data.density_workspace.view,
      data.accepted_advection_workspace.view,
      data.previous_advection_workspace.view,
      {call.passive_output.data(), call.passive_output.size()},
      data.low_density_workspace.view,
      data.low_enthalpy_workspace.view,
      {call.low_species_output.data(), call.low_species_output.size()},
      {call.low_passive_output.data(), call.low_passive_output.size()}};
  FaceFluxView paired_flux;
  if (!fixture.paired_flux.workspace_view(
          0U, time == std::numeric_limits<RevisionToken>::max()
                  ? RevisionToken{1U}
                  : time + 1U,
          paired_flux)) {
    return call;
  }
  call.output.paired_mass_flux = paired_flux;
  call.slow_path.halo = &fixture.donor_halo;
  call.slow_path.halo_stage = 16U;
  call.slow_path.exchange_fields = {call.slow_exchange_fields.data(),
                                    call.slow_exchange_fields.size()};
  call.slow_path.enthalpy_endpoint = &fixture.endpoint;
  return call;
}

double high_rho(double rho_n, double rho_previous, double mass_n,
                double mass_previous) {
  const double extrapolate_accepted = 2.0;
  const double extrapolate_previous = -1.0;
  return (-kBdf.a1 * rho_n - kBdf.a2 * rho_previous -
          extrapolate_accepted * mass_n -
          extrapolate_previous * mass_previous) /
         kBdf.a0;
}

double high_rho_quantity(double rho_n, double q_n, double rho_previous,
                         double q_previous, double transport_n,
                         double transport_previous) {
  return (-kBdf.a1 * rho_n * q_n - kBdf.a2 * rho_previous * q_previous -
          2.0 * transport_n - (-1.0) * transport_previous) /
         kBdf.a0;
}

double bdf_accepted_rate_density(double rho_n, double rho_previous,
                                 double mass_divergence) {
  return (-kBdf.a1 * rho_n - kBdf.a2 * rho_previous - mass_divergence) /
         kBdf.a0;
}

double bdf_accepted_rate_quantity(double rho_n, double q_n,
                                  double rho_previous, double q_previous,
                                  double transport, double rhs) {
  return (-kBdf.a1 * rho_n * q_n - kBdf.a2 * rho_previous * q_previous -
          (transport - rhs)) /
         kBdf.a0;
}

bool bound_for_state(const ThermodynamicsPlan& thermodynamics, double rho,
                     double independent_species_density, double& lower,
                     double& upper) {
  double independent_min = 0.0;
  double independent_max = 0.0;
  double dependent_min = 0.0;
  double dependent_max = 0.0;
  if (!thermodynamics.independent_species_enthalpy_bounds(
          0U, independent_min, independent_max) ||
      !thermodynamics.species_enthalpy_bounds(
          thermodynamics.dependent_species_index(), dependent_min,
          dependent_max)) {
    return false;
  }
  const double dependent_density = rho - independent_species_density;
  lower = independent_species_density * independent_min +
          dependent_density * dependent_min;
  upper = independent_species_density * independent_max +
          dependent_density * dependent_max;
  return std::isfinite(lower) && std::isfinite(upper) && lower <= upper;
}

bool commit_shared_face_transfer(const CartesianKernelPlan& kernels,
                                 Int3 cells, FluxHistory& history,
                                 ConstFaceFluxView& committed,
                                 double target_divergence =
                                     kSharedFaceDivergence) {
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

  // The two active cells are (0,0,0) and (1,0,0).  Their only nonzero
  // transfer is the internal x face; the two external x-face values cancel
  // against ghost products, leaving one reciprocal shared-face contribution.
  OwnedField density = make_field(kDensity, cells, 2U, 801U, 901U);
  OwnedField velocity = make_field(kVelocity, cells, 2U, 802U, 902U);
  const std::size_t scalar_size = velocity.storage.size();
  velocity.storage.resize(3U * scalar_size, 0.0);
  velocity.view.base = velocity.storage.data() + 2U +
                       2U * velocity.view.stride_y +
                       2U * velocity.view.stride_z;
  velocity.view.components = 3U;
  velocity.view.component_stride = scalar_size;
  fill(density, kSharedFaceRho);

  constexpr double cell_volume = 1.0 / 8.0;
  constexpr double face_area = 1.0 / 4.0;
  const double face_flux = target_divergence * cell_volume;
  const double normal_velocity =
      face_flux / (kSharedFaceRho * face_area);
  for (std::int32_t x = 0; x < cells.x; ++x) {
    velocity.view.unchecked({x, 0, 0}, 0U) = normal_velocity;
  }
  velocity.view.unchecked({-1, 0, 0}, 0U) = -normal_velocity;
  velocity.view.unchecked({cells.x, 0, 0}, 0U) = -normal_velocity;

  const std::array<ConstFieldView, 2U> reads{as_const(density.view),
                                             as_const(velocity.view)};
  const KernelInvocation invocation{
      {reads.data(), reads.size()}, {}, {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U,
      nullptr};
  const std::array dependencies{dependency};
  return static_cast<bool>(reconstruct_mass_flux(kernels, invocation, pending)) &&
         static_cast<bool>(history.writer.publish_pending(
             {dependencies.data(), dependencies.size()}, pending)) &&
         static_cast<bool>(history.transaction.collective_finish(
             MPI_COMM_SELF, Status{})) &&
         static_cast<bool>(history.writer.committed(history.storage, committed));
}

bool commit_x_boundary_transfer(Fixture& fixture, FluxHistory& history,
                                double lower_face_flux,
                                double upper_face_flux,
                                ConstFaceFluxView& committed) {
  if (fixture.patch.cells.x != 2 ||
      !history.transaction.begin(history.layers) ||
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

  OwnedField density =
      make_field(kDensity, fixture.patch.cells, 2U, 811U, 911U);
  OwnedField velocity =
      make_field_components(kVelocity, fixture.patch.cells, 3U, 2U, 812U,
                            912U);
  fill(density, kSharedFaceRho);
  fill(velocity, 0.0);
  const Int3 global = fixture.geometry.global_cells();
  const double face_area =
      1.0 / (static_cast<double>(global.y) * global.z);
  const double lower_velocity =
      lower_face_flux / (kSharedFaceRho * face_area);
  const double upper_velocity =
      upper_face_flux / (kSharedFaceRho * face_area);
  // The two interior velocities are zero.  Mirrored ghost values therefore
  // prescribe exactly the requested boundary-face averages while the local
  // internal x face stays zero.
  velocity.view.unchecked({-1, 0, 0}, 0U) = 2.0 * lower_velocity;
  velocity.view.unchecked({fixture.patch.cells.x, 0, 0}, 0U) =
      2.0 * upper_velocity;

  const std::array<ConstFieldView, 2U> reads{as_const(density.view),
                                             as_const(velocity.view)};
  const KernelInvocation invocation{
      {reads.data(), reads.size()}, {}, {{0, 0, 0}, fixture.patch.cells},
      0U, 0U, 1U, 0U, nullptr};
  const std::array dependencies{dependency};
  return static_cast<bool>(reconstruct_mass_flux(
             fixture.equations.kernels(), invocation, pending)) &&
         static_cast<bool>(history.writer.publish_pending(
             {dependencies.data(), dependencies.size()}, pending)) &&
         static_cast<bool>(history.transaction.collective_finish(
             MPI_COMM_SELF, Status{})) &&
         static_cast<bool>(history.writer.committed(history.storage,
                                                    committed));
}

double mass_divergence(ConstFaceFluxView flux, Int3 cell) {
  constexpr double cell_volume = 1.0 / 8.0;
  const double integral =
      flux.x.unchecked({cell.x + 1, cell.y, cell.z}) -
      flux.x.unchecked(cell) + flux.y.unchecked({cell.x, cell.y + 1, cell.z}) -
      flux.y.unchecked(cell) + flux.z.unchecked({cell.x, cell.y, cell.z + 1}) -
      flux.z.unchecked(cell);
  return integral / cell_volume;
}

bool test_shared_face_density_red_regression(Fixture& fixture,
                                             std::uintptr_t exchange_plan) {
  FluxHistory history;
  if (!expect(make_flux_history(fixture.patch.cells, history),
              "shared-face red fixture allocates committed flux history")) {
    return false;
  }
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  if (!expect(commit_shared_face_transfer(fixture.equations.kernels(),
                                          fixture.patch.cells, history,
                                          previous_flux) &&
                  commit_shared_face_transfer(fixture.equations.kernels(),
                                              fixture.patch.cells, history,
                                              accepted_flux),
              "shared-face red fixture publishes two transfer histories")) {
    return false;
  }

  // A previous-history source makes only the EX2 high state inadmissible.
  // The accepted-level low source is exactly zero, so once density is
  // repaired the low enthalpy remains the constant admissible 300 kJ/kg.
  PredictorData data(fixture.patch.cells, 300000.0, 300000.0, 0.0,
                     4000000.0);
  fill(data.rho, kSharedFaceRho);
  fill(data.rho_previous, kSharedFaceRho);
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                 exchange_plan, 908U);
  call.input.dt = kSharedFaceDt;
  call.input.bdf = kSharedFaceBdf;
  const Status status = fixture.equations.thermophysical_predictor().predict(
      MPI_COMM_SELF, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  const ThermophysicalPredictorFailure& failure = call.diagnostics.failure;
  const double donor_divergence = mass_divergence(
      accepted_flux, {0, 0, 0});
  const double neighbor_divergence = mass_divergence(
      accepted_flux, {1, 0, 0});
  const double expected_rho_low =
      kSharedFaceRho - kSharedFaceDt * kSharedFaceDivergence;
  const std::uint32_t required_scalars =
      ThermophysicalPredictorFailure::scalar_density_current |
      ThermophysicalPredictorFailure::scalar_density_next |
      ThermophysicalPredictorFailure::scalar_divergence |
      ThermophysicalPredictorFailure::scalar_observed_value |
      ThermophysicalPredictorFailure::scalar_allowed_lower;
  const bool current_failure_witness =
      status.code == StatusCode::numerical_failure &&
      static_cast<unsigned>(status.code) == 5U && status.detail == 985U &&
      failure.valid &&
      failure.reason == ThermophysicalPredictorFailureReason::low_density_endpoint &&
      failure.field == ThermophysicalPredictorFailureField::density &&
      failure.global_index.x == 0 && failure.global_index.y == 0 &&
      failure.global_index.z == 0 &&
      (failure.scalar_mask & required_scalars) == required_scalars &&
      close(failure.density_current, kSharedFaceRho) &&
      close(failure.density_next, expected_rho_low) &&
      close(failure.observed_value, expected_rho_low) &&
      close(failure.divergence, kSharedFaceDivergence) &&
      close(donor_divergence, kSharedFaceDivergence) &&
      close(neighbor_divergence, -kSharedFaceDivergence) &&
      close(accepted_flux.x.unchecked({1, 0, 0}),
            kSharedFaceDivergence * (1.0 / 8.0));
  std::cerr << "shared_face_red_baseline status="
            << static_cast<unsigned>(status.code) << '/' << status.detail
            << " reason=" << static_cast<unsigned>(failure.reason)
            << " cell=" << failure.global_index.x << ','
            << failure.global_index.y << ',' << failure.global_index.z
            << " rho_n=" << failure.density_current
            << " rho_low=" << failure.density_next
            << " divergence=" << failure.divergence
            << " witness=" << (current_failure_witness ? 1 : 0) << '\n';

  bool all_density_positive = true;
  for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < fixture.patch.cells.x; ++x) {
        const double rho =
            call.output.density_workspace.unchecked({x, y, z}, 0U);
        all_density_positive = all_density_positive && std::isfinite(rho) &&
                               rho > 0.0;
      }
    }
  }
  const double input_pair_sum =
      data.rho.view.unchecked({0, 0, 0}, 0U) +
      data.rho.view.unchecked({1, 0, 0}, 0U);
  const double output_pair_sum =
      call.output.density_workspace.unchecked({0, 0, 0}, 0U) +
      call.output.density_workspace.unchecked({1, 0, 0}, 0U);
  const bool fixed_dt_expectation =
      static_cast<bool>(status) && !failure.valid &&
      call.input.dt == kSharedFaceDt && all_density_positive &&
      close(output_pair_sum, input_pair_sum);
  bool passed = expect(current_failure_witness || fixed_dt_expectation,
                       "shared-face trace reproduces baseline or fixed path");
  passed &= expect(fixed_dt_expectation,
                   "fixed dt keeps both shared-face densities positive and conserved");
  passed &= expect(static_cast<bool>(status) && !failure.valid,
                   "fixed dt shared-face predictor has no numerical failure");
  return passed;
}

bool test_bdf_accepted_rate_homotopy(Fixture& fixture,
                                      std::uintptr_t exchange_plan) {
  constexpr double target_divergence = 300.0;
  FluxHistory history;
  if (!expect(make_flux_history(fixture.patch.cells, history),
              "BDF-homotopy fixture allocates committed flux history")) {
    return false;
  }
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  if (!expect(commit_shared_face_transfer(
                  fixture.equations.kernels(), fixture.patch.cells, history,
                  previous_flux, target_divergence) &&
                  commit_shared_face_transfer(
                      fixture.equations.kernels(), fixture.patch.cells, history,
                      accepted_flux, target_divergence),
              "BDF-homotopy fixture publishes two transfer histories")) {
    return false;
  }

  // The accepted-rate BDF endpoint has a negative donor density at this
  // deliberately large CFL.  A single accepted-anchor coefficient must act
  // on the complete conserved tuple; the paired constant scalar therefore
  // remains exactly constant without damping mass flux.
  PredictorData data(fixture.patch.cells, 300000.0, 300000.0, 0.0,
                     4000000.0);
  fill(data.rho, kSharedFaceRho);
  fill(data.rho_previous, kSharedFaceRho);
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                 exchange_plan, 909U);
  call.input.dt = kSharedFaceDt;
  call.input.bdf = kSharedFaceBdf;
  const Status status = fixture.equations.thermophysical_predictor().predict(
      MPI_COMM_SELF, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  bool passed = expect(static_cast<bool>(status),
                       "BDF homotopy succeeds after endpoint rejection");
  passed &= expect(!call.diagnostics.failure.valid,
                   "BDF homotopy has no numerical failure provenance");
  passed &= expect(call.diagnostics.low_state ==
                           ThermophysicalLowStateKind::
                               bdf_accepted_rate_homotopy &&
                       call.diagnostics.mass_flux_scale == 1.0 &&
                       call.diagnostics.bdf_endpoint_alpha >= 0.0 &&
                       call.diagnostics.bdf_endpoint_alpha < 1.0,
                   "BDF homotopy publishes its state, alpha, and unit mass scale");
  passed &= expect(call.diagnostics.low_order_substeps == 1U &&
                       call.diagnostics.low_order_halo_exchanges == 0U,
                   "BDF homotopy reuses one transport pass and no halo");

  const double donor_divergence = mass_divergence(accepted_flux, {0, 0, 0});
  const double neighbor_divergence = mass_divergence(accepted_flux, {1, 0, 0});
  const double donor_low =
      call.output.low_order_density_workspace.unchecked({0, 0, 0}, 0U);
  const double neighbor_low =
      call.output.low_order_density_workspace.unchecked({1, 0, 0}, 0U);
  passed &= expect(close(donor_divergence, target_divergence) &&
                       close(neighbor_divergence, -target_divergence),
                   "BDF homotopy uses reciprocal shared-face mass divergence");

  const double bdf_donor =
      (-kSharedFaceBdf.a1 * kSharedFaceRho -
       kSharedFaceBdf.a2 * kSharedFaceRho - donor_divergence) /
      kSharedFaceBdf.a0;
  const double bdf_neighbor =
      (-kSharedFaceBdf.a1 * kSharedFaceRho -
       kSharedFaceBdf.a2 * kSharedFaceRho - neighbor_divergence) /
      kSharedFaceBdf.a0;
  constexpr double homotopy_safety =
      1.0 - 64.0 * std::numeric_limits<double>::epsilon();
  constexpr double accepted_species_fraction = 0.2;
  const double accepted_species_density =
      accepted_species_fraction * kSharedFaceRho;
  const double bdf_donor_species_density =
      accepted_species_fraction * bdf_donor;
  const double representation_scale =
      std::max({1.0, std::abs(kSharedFaceRho), std::abs(bdf_donor),
                std::abs(accepted_species_density),
                std::abs(bdf_donor_species_density)});
  const double target_species_margin = std::min(
      accepted_species_density,
      512.0 * std::numeric_limits<double>::epsilon() *
          representation_scale);
  const double expected_bdf_alpha = std::nextafter(
      ((accepted_species_density - target_species_margin) /
       (accepted_species_density - bdf_donor_species_density)) *
          homotopy_safety,
      0.0);
  const double expected_homotopy_donor =
      kSharedFaceRho +
      expected_bdf_alpha * (bdf_donor - kSharedFaceRho);
  const double expected_homotopy_neighbor =
      kSharedFaceRho +
      expected_bdf_alpha * (bdf_neighbor - kSharedFaceRho);
  passed &= expect(bdf_donor < 0.0 && bdf_neighbor > 0.0 &&
                       close(call.diagnostics.bdf_endpoint_alpha,
                             expected_bdf_alpha, 2.0e-13) &&
                       close(donor_low, expected_homotopy_donor, 2.0e-13) &&
                       close(neighbor_low, expected_homotopy_neighbor,
                             2.0e-13) &&
                       donor_low > 0.0,
                   "BDF guarded bundle selects the certified common alpha");

  const double input_pair_sum =
      data.rho.view.unchecked({0, 0, 0}, 0U) +
      data.rho.view.unchecked({1, 0, 0}, 0U);
  const double low_pair_sum = donor_low + neighbor_low;
  const double final_donor =
      call.output.density_workspace.unchecked({0, 0, 0}, 0U);
  const double final_neighbor =
      call.output.density_workspace.unchecked({1, 0, 0}, 0U);
  const auto high_density = [&](double rho_n, double rho_previous,
                                double divergence,
                                double previous_divergence) {
    return (-kSharedFaceBdf.a1 * rho_n -
            kSharedFaceBdf.a2 * rho_previous - 2.0 * divergence +
            previous_divergence) /
           kSharedFaceBdf.a0;
  };
  const double high_donor =
      high_density(kSharedFaceRho, kSharedFaceRho, target_divergence,
                   target_divergence);
  const double high_neighbor =
      high_density(kSharedFaceRho, kSharedFaceRho, -target_divergence,
                   -target_divergence);
  const double expected_final_donor =
      donor_low + call.diagnostics.theta * (high_donor - donor_low);
  const double expected_final_neighbor =
      neighbor_low + call.diagnostics.theta * (high_neighbor - neighbor_low);
  const double final_pair_sum = final_donor + final_neighbor;
  passed &= expect(close(low_pair_sum, input_pair_sum, 2.0e-12) &&
                       close(final_pair_sum, input_pair_sum, 2.0e-12),
                   "BDF homotopy and common theta remain pair-conservative");
  passed &= expect(close(final_donor, expected_final_donor, 2.0e-12) &&
                       close(final_neighbor, expected_final_neighbor, 2.0e-12),
                   "common theta blends homotopy mass with the high BDF mass");
  for (const Int3 cell : {Int3{0, 0, 0}, Int3{1, 0, 0}}) {
    const double rho =
        call.output.density_workspace.unchecked(cell, 0U);
    const double scalar = call.output.independent_species.data[0U]
                              .unchecked(cell, 0U);
    const double species_density = rho * scalar;
    passed &= expect(
        std::isfinite(scalar) && scalar >= 0.0 && scalar <= 1.0 &&
            std::abs(species_density - 0.2 * rho) <=
                256.0 * std::numeric_limits<double>::epsilon() *
                    kSharedFaceRho,
        "complete-tuple homotopy preserves paired conserved species mass");
  }
  std::cerr << "bdf_homotopy_oracle alpha="
            << call.diagnostics.bdf_endpoint_alpha
            << " donor_low=" << donor_low << " neighbor_low=" << neighbor_low
            << " theta=" << call.diagnostics.theta << '\n';
  return passed;
}

bool test_fast_path_bit_identity(Fixture& fixture,
                                 ConstFaceFluxView accepted_flux,
                                 ConstFaceFluxView previous_flux,
                                 std::uintptr_t exchange_plan) {
  PredictorData data(fixture.patch.cells, 300000.0, 299000.0, 0.0, 0.0);
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                  exchange_plan, 901U);
  const Status status = fixture.equations.thermophysical_predictor().predict(
      MPI_COMM_SELF, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  bool passed = expect(static_cast<bool>(status),
                       "normal BDF2/EX2 predictor succeeds");
  passed &= expect(call.certificate.valid(),
                   "normal predictor publishes certificate");
  passed &= expect(call.diagnostics.theta == 1.0 &&
                       call.diagnostics.mass_flux_scale == 1.0 &&
                       call.diagnostics.bdf_endpoint_alpha == 1.0 &&
                       call.diagnostics.source_endpoint_alpha == 1.0 &&
                       call.diagnostics.low_state ==
                           ThermophysicalLowStateKind::none &&
                       !call.diagnostics.limited &&
                       !call.diagnostics.failure.valid &&
                       call.diagnostics.enthalpy_endpoint_alpha == 1.0 &&
                       call.diagnostics.low_order_substeps == 0U &&
                       call.diagnostics.low_order_halo_exchanges == 0U &&
                       call.diagnostics.low_order_transport_passes == 0U &&
                       call.diagnostics.enthalpy_solve_calls == 0U &&
                       call.diagnostics.blocking_collectives == 1U,
                   "normal path has theta=1, no endpoint work, and one collective");

  const double expected_rho = high_rho(1.0, 1.0, 0.0, 0.0);
  const double expected_rhoh = high_rho_quantity(
      1.0, 300000.0, 1.0, 299000.0, 0.0, 0.0);
  const double expected_h = expected_rhoh / expected_rho;
  for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < fixture.patch.cells.x; ++x) {
        const Int3 cell{x, y, z};
        passed &= expect(same_bits(
                            call.output.density_workspace.unchecked(cell, 0U),
                            expected_rho),
                        "normal density follows the old BDF2 formula bitwise");
        passed &= expect(same_bits(call.output.enthalpy.unchecked(cell, 0U),
                                   expected_h),
                         "normal enthalpy follows the old BDF2/EX2 formula bitwise");
        passed &= expect(call.output.independent_species.data[0U]
                                     .unchecked(cell, 0U) == 0.2,
                         "normal species remains unchanged bitwise");
        passed &= expect(call.output.passive_scalars.data[0U]
                                     .unchecked(cell, 0U) == 0.4,
                         "normal passive scalar remains unchanged bitwise");
      }
    }
  }
  const ConstFaceFluxView paired = as_const(call.output.paired_mass_flux);
  for (const auto axis : {CartesianAxis::x, CartesianAxis::y,
                          CartesianAxis::z}) {
    const ConstFaceFieldView actual =
        axis == CartesianAxis::x
            ? paired.x
            : (axis == CartesianAxis::y ? paired.y : paired.z);
    const ConstFaceFieldView accepted =
        axis == CartesianAxis::x
            ? accepted_flux.x
            : (axis == CartesianAxis::y ? accepted_flux.y : accepted_flux.z);
    const ConstFaceFieldView previous =
        axis == CartesianAxis::x
            ? previous_flux.x
            : (axis == CartesianAxis::y ? previous_flux.y : previous_flux.z);
    for (std::int32_t z = 0; z < actual.extents.z; ++z)
      for (std::int32_t y = 0; y < actual.extents.y; ++y)
        for (std::int32_t x = 0; x < actual.extents.x; ++x) {
          const Int3 face{x, y, z};
          passed &= expect(
              same_bits(actual.unchecked(face),
                        2.0 * accepted.unchecked(face) -
                            previous.unchecked(face)),
              "normal paired face flux follows EX2 bitwise");
        }
  }
  return passed;
}

bool test_small_dt_paired_continuity_roundoff(
    Fixture& fixture, ConstFaceFluxView accepted_flux,
    ConstFaceFluxView previous_flux, std::uintptr_t exchange_plan) {
  PredictorData data(fixture.patch.cells, 300000.0, 300000.0, 0.0, 0.0);
  fill(data.rho, 1.0);
  fill(data.rho_previous, std::nextafter(1.0, 2.0));
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                  exchange_plan, 910U);
  constexpr double dt = 1.0e-12;
  call.input.dt = dt;
  call.input.bdf = {1.5 / dt, -2.0 / dt, 0.5 / dt, 2U};

  const Status status = fixture.equations.thermophysical_predictor().predict(
      MPI_COMM_SELF, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  if (!status) {
    std::cerr << "small-dt paired continuity status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " failure=" << call.diagnostics.failure.valid << '/'
              << static_cast<unsigned>(call.diagnostics.failure.reason)
              << " residual=" << call.diagnostics.failure.observed_value
              << " divergence=" << call.diagnostics.failure.divergence
              << '\n';
  }
  bool passed = expect(
      static_cast<bool>(status),
      "small-dt BDF2 paired continuity accepts roundoff from large history terms");
  passed &= expect(call.certificate.valid() &&
                       !call.diagnostics.failure.valid &&
                       call.diagnostics.theta == 1.0 &&
                       call.diagnostics.mass_flux_scale == 1.0,
                   "small-dt algebra audit preserves the unlimited predictor");
  for (std::int32_t z = 0; z < fixture.patch.cells.z && passed; ++z) {
    for (std::int32_t y = 0; y < fixture.patch.cells.y && passed; ++y) {
      for (std::int32_t x = 0; x < fixture.patch.cells.x; ++x) {
        const double rho =
            call.output.density_workspace.unchecked({x, y, z}, 0U);
        passed &= expect(std::isfinite(rho) && rho > 0.0,
                         "small-dt paired density remains admissible");
      }
    }
  }
  return passed;
}

bool test_local_donor_paired_flux(Fixture& fixture,
                                  std::uintptr_t exchange_plan,
                                  MPI_Comm communicator) {
  constexpr double target_divergence = 300.0;
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return false;
  }
  FluxHistory history;
  if (!expect(make_flux_history(fixture.patch.cells, history),
              "local-donor fixture allocates flux history"))
    return false;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  const Int3 global_cells = fixture.geometry.global_cells();
  const double volume =
      1.0 / (static_cast<double>(global_cells.x) * global_cells.y *
             global_cells.z);
  const double shared_flux = target_divergence * volume;
  const auto commit = [&](ConstFaceFluxView& flux) {
    if (size == 1) {
      return commit_shared_face_transfer(
          fixture.equations.kernels(), fixture.patch.cells, history, flux,
          target_divergence);
    }
    if (size != 2) return false;
    return commit_x_boundary_transfer(
        fixture, history, rank == 1 ? shared_flux : 0.0,
        rank == 0 ? shared_flux : 0.0, flux);
  };
  if (!expect(commit(previous_flux) && commit(accepted_flux),
              "local-donor fixture publishes paired histories"))
    return false;

  PredictorData data(fixture.patch.cells, 300000.0, 300000.0, 0.0, 0.0);
  fill(data.rho, kSharedFaceRho);
  fill(data.rho_previous, kSharedFaceRho);
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                 exchange_plan, 990U);
  call.input.dt = kSharedFaceDt;
  call.input.bdf = kSharedFaceBdf;
  const Status status = fixture.equations.thermophysical_predictor().predict(
      communicator, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  bool passed = expect(status && call.certificate.valid(),
                       "local donor route publishes a valid certificate");
  passed &= expect(
      call.diagnostics.low_state ==
              ThermophysicalLowStateKind::bdf_local_donor_flux &&
          call.diagnostics.low_order_halo_exchanges == 1U &&
          call.diagnostics.bdf_endpoint_alpha == 1.0 &&
          call.diagnostics.mass_flux_scale >= 0.0 &&
          call.diagnostics.mass_flux_scale < 1.0,
      "local donor route reports one factor halo and its minimum factor");
  const ConstFaceFluxView paired = as_const(call.output.paired_mass_flux);
  if (size == 2) {
    const Int3 shared_face =
        rank == 0 ? Int3{fixture.patch.cells.x, 0, 0} : Int3{0, 0, 0};
    const double actual_shared = paired.x.unchecked(shared_face);
    const double accepted_shared = accepted_flux.x.unchecked(shared_face);
    std::uint64_t shared_bits = 0U;
    std::memcpy(&shared_bits, &actual_shared, sizeof(shared_bits));
    std::uint64_t minimum_bits = 0U;
    std::uint64_t maximum_bits = 0U;
    const bool identical =
        MPI_Allreduce(&shared_bits, &minimum_bits, 1, MPI_UINT64_T, MPI_MIN,
                      communicator) == MPI_SUCCESS &&
        MPI_Allreduce(&shared_bits, &maximum_bits, 1, MPI_UINT64_T, MPI_MAX,
                      communicator) == MPI_SUCCESS &&
        minimum_bits == maximum_bits;
    const double expected_shared =
        accepted_shared *
        (call.diagnostics.mass_flux_scale +
         call.diagnostics.theta *
             (1.0 - call.diagnostics.mass_flux_scale));
    passed &= expect(
        identical && accepted_shared > 0.0 && actual_shared > 0.0 &&
            actual_shared < accepted_shared &&
            close(actual_shared, expected_shared, 2.0e-11),
        "MPI-periodic shared face uses the upstream donor factor bitwise");
  }
  for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < fixture.patch.cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double integral =
            paired.x.unchecked({x + 1, y, z}) - paired.x.unchecked(cell) +
            paired.y.unchecked({x, y + 1, z}) - paired.y.unchecked(cell) +
            paired.z.unchecked({x, y, z + 1}) - paired.z.unchecked(cell);
        const double rho_star =
            call.output.density_workspace.unchecked(cell, 0U);
        const double residual =
            volume * (kSharedFaceBdf.a0 * rho_star +
                      kSharedFaceBdf.a1 * kSharedFaceRho +
                      kSharedFaceBdf.a2 * kSharedFaceRho) +
            integral;
        const double species = call.output.independent_species.data[0U]
                                   .unchecked(cell, 0U);
        double lower = 0.0;
        double upper = 0.0;
        passed &= expect(
            std::abs(residual) <=
                    4096.0 * std::numeric_limits<double>::epsilon() *
                        std::max(1.0, std::abs(integral)) &&
                std::isfinite(rho_star) && rho_star > 0.0 &&
                std::isfinite(species) && species >= 0.0 &&
                species <= 1.0 &&
                bound_for_state(fixture.thermodynamics, rho_star,
                                rho_star * species, lower, upper) &&
                rho_star * call.output.enthalpy.unchecked(cell, 0U) >= lower &&
                rho_star * call.output.enthalpy.unchecked(cell, 0U) <= upper,
            "local donor state and paired flux satisfy one exact full-tuple identity");
      }
    }
  }
  return passed;
}

bool test_physical_inflow_outflow_factors(Fixture& fixture,
                                          std::uintptr_t exchange_plan) {
  constexpr double target_outflow_divergence = 300.0;
  const Int3 global = fixture.geometry.global_cells();
  const double volume =
      1.0 / (static_cast<double>(global.x) * global.y * global.z);
  const double outflow = target_outflow_divergence * volume;
  const double inflow = 0.25 * outflow;
  FluxHistory history;
  if (!expect(make_flux_history(fixture.patch.cells, history),
              "physical-factor fixture allocates flux history")) {
    return false;
  }
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  if (!expect(commit_x_boundary_transfer(fixture, history, inflow, outflow,
                                         previous_flux) &&
                  commit_x_boundary_transfer(fixture, history, inflow,
                                             outflow, accepted_flux),
              "physical-factor fixture publishes boundary flux histories")) {
    return false;
  }

  PredictorData data(fixture.patch.cells, 300000.0, 300000.0, 0.0, 0.0);
  fill(data.rho, kSharedFaceRho);
  fill(data.rho_previous, kSharedFaceRho);
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                 exchange_plan, 992U);
  call.input.dt = kSharedFaceDt;
  call.input.bdf = kSharedFaceBdf;
  const Status status = fixture.equations.thermophysical_predictor().predict(
      MPI_COMM_SELF, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  if (!expect(status && call.certificate.valid() &&
                  call.diagnostics.low_state ==
                      ThermophysicalLowStateKind::bdf_local_donor_flux &&
                  call.diagnostics.mass_flux_scale >= 0.0 &&
                  call.diagnostics.mass_flux_scale < 1.0 &&
                  call.diagnostics.theta >= 0.0 &&
                  call.diagnostics.theta < 1.0,
              "physical-factor fixture activates the local donor route")) {
    return false;
  }
  const ConstFaceFluxView paired = as_const(call.output.paired_mass_flux);
  const double accepted_inflow = accepted_flux.x.unchecked({0, 0, 0});
  const double paired_inflow = paired.x.unchecked({0, 0, 0});
  const Int3 high_face{fixture.patch.cells.x, 0, 0};
  const double accepted_outflow = accepted_flux.x.unchecked(high_face);
  const double paired_outflow = paired.x.unchecked(high_face);
  const double expected_outflow =
      accepted_outflow *
      (call.diagnostics.mass_flux_scale +
       call.diagnostics.theta *
           (1.0 - call.diagnostics.mass_flux_scale));
  return expect(accepted_inflow > 0.0 &&
                    same_bits(paired_inflow, accepted_inflow),
                "physical inflow keeps factor one") &&
         expect(accepted_outflow > 0.0 && paired_outflow > 0.0 &&
                    paired_outflow < accepted_outflow &&
                    close(paired_outflow, expected_outflow, 2.0e-11),
                "physical outflow uses its owner donor factor");
}

bool test_invalid_zero_outgoing_base(Fixture& fixture,
                                     ConstFaceFluxView accepted_flux,
                                     ConstFaceFluxView previous_flux,
                                     std::uintptr_t exchange_plan,
                                     MPI_Comm communicator) {
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return false;
  }
  constexpr double accepted_enthalpy = 300000.0;
  const double previous_enthalpy =
      rank == size - 1 ? 1000000.0 : accepted_enthalpy;
  PredictorData data(fixture.patch.cells, accepted_enthalpy,
                     previous_enthalpy, 0.0, 0.0);
  double lower = 0.0;
  double upper = 0.0;
  const bool input_states_admissible =
      bound_for_state(fixture.thermodynamics, 1.0, 0.2, lower, upper) &&
      accepted_enthalpy >= lower && accepted_enthalpy <= upper &&
      previous_enthalpy >= lower && previous_enthalpy <= upper;
  const double source_free_bdf =
      (-kBdf.a1 * accepted_enthalpy -
       kBdf.a2 * previous_enthalpy) /
      kBdf.a0;
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                 exchange_plan, 991U);
  const Status status = fixture.equations.thermophysical_predictor().predict(
      communicator, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  return expect(input_states_admissible &&
                    (rank == size - 1 ? source_free_bdf < lower
                                      : source_free_bdf >= lower &&
                                            source_free_bdf <= upper),
                "reason-17 fixture has admissible histories and only its highest rank has an inadmissible source-free BDF2 anchor") &&
         expect(status.code == StatusCode::numerical_failure &&
                    call.diagnostics.failure.valid &&
                    call.diagnostics.failure.rank == size - 1 &&
                    call.diagnostics.failure.reason ==
                        ThermophysicalPredictorFailureReason::
                            low_bdf_source_base_admissibility,
                "collective source-free BDF2 anchor failure preserves exact reason 17 and its owner");
}

bool test_complete_source_bundle_globalization(
    Fixture& fixture, std::uintptr_t exchange_plan, MPI_Comm communicator,
    int rank, int size) {
  FluxHistory history;
  if (!expect(make_flux_history(fixture.patch.cells, history),
              "source-bundle fixture allocates zero-flux history")) {
    return false;
  }
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  if (!expect(commit_zero_flux(fixture.equations.kernels(),
                               fixture.patch.cells, history, previous_flux) &&
                  commit_zero_flux(fixture.equations.kernels(),
                                   fixture.patch.cells, history,
                                   accepted_flux),
              "source-bundle fixture publishes zero mass-flux histories")) {
    return false;
  }

  constexpr double accepted_enthalpy = 300000.0;
  constexpr double accepted_species = 0.2;
  constexpr double accepted_passive = 0.4;
  constexpr double species_rhs = 0.5;
  constexpr double passive_rhs = -0.75;
  const double enthalpy_rhs = rank == size - 1 ? -2000000.0 : -1400000.0;
  PredictorData data(fixture.patch.cells, accepted_enthalpy,
                     accepted_enthalpy, enthalpy_rhs, 0.0);
  fill(data.species, accepted_species);
  fill(data.species_previous, accepted_species);
  fill(data.species_rhs, species_rhs);
  fill(data.passive, accepted_passive);
  fill(data.passive_previous, accepted_passive);
  fill(data.passive_rhs, passive_rhs);

  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                 exchange_plan, 993U);
  call.input.bdf = {1.0 / kDt, -1.0 / kDt, 0.0, 1U};
  call.input.density_previous = {};
  call.input.enthalpy_previous = {};
  call.input.mass_flux_previous = {};
  call.input.enthalpy_ghosts.previous = {};
  call.species_previous[0U] = {};
  call.passive_previous[0U] = {};
  call.species_ghosts[0U].previous = {};
  call.passive_ghosts[0U].previous = {};
  call.species_rates[0U].previous = {};
  call.passive_rates[0U].previous = {};
  call.input.enthalpy_nonadvective_rhs.previous = {};
  const Status status = fixture.equations.thermophysical_predictor().predict(
      communicator, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  if (!status) {
    std::cerr << "rank " << rank << " source_bundle_status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " failure_reason="
              << static_cast<unsigned>(call.diagnostics.failure.reason)
              << " failure_rank=" << call.diagnostics.failure.rank << '\n';
  }

  const double a0 = call.input.bdf.a0;
  const double base_rhoh = accepted_enthalpy;
  const double full_rhoh = base_rhoh + enthalpy_rhs / a0;
  const double base_species = accepted_species;
  const double full_species = base_species + species_rhs / a0;
  const double base_passive = accepted_passive;
  const double full_passive = base_passive + passive_rhs / a0;
  double base_lower = 0.0;
  double base_upper = 0.0;
  double full_lower = 0.0;
  double full_upper = 0.0;
  bool oracle_valid =
      bound_for_state(fixture.thermodynamics, 1.0, base_species, base_lower,
                      base_upper) &&
      bound_for_state(fixture.thermodynamics, 1.0, full_species, full_lower,
                      full_upper);
  constexpr double roundoff_safety =
      1.0 - 64.0 * std::numeric_limits<double>::epsilon();
  const auto oracle_factor = [&](double anchor_margin, double trial_margin,
                                 bool strict) {
    const double scale =
        std::max({1.0, std::abs(anchor_margin), std::abs(trial_margin)});
    const double guard =
        512.0 * std::numeric_limits<double>::epsilon() * scale;
    const double target =
        std::min(anchor_margin, std::max(guard, 0.1 * anchor_margin));
    if (!std::isfinite(anchor_margin) || !std::isfinite(trial_margin) ||
        (strict ? !(anchor_margin > 0.0) : anchor_margin < 0.0)) {
      oracle_valid = false;
      return -1.0;
    }
    if (trial_margin >= target) return 1.0;
    const double denominator = anchor_margin - trial_margin;
    if (!std::isfinite(denominator) || !(denominator > 0.0)) {
      oracle_valid = false;
      return -1.0;
    }
    return std::nextafter(
        std::clamp((anchor_margin - target) / denominator, 0.0, 1.0) *
            roundoff_safety,
        0.0);
  };
  double expected_local_alpha = 1.0;
  expected_local_alpha =
      std::min(expected_local_alpha, oracle_factor(1.0, 1.0, true));
  expected_local_alpha = std::min(
      expected_local_alpha,
      oracle_factor(base_species, full_species, false));
  expected_local_alpha = std::min(
      expected_local_alpha,
      oracle_factor(1.0 - base_species, 1.0 - full_species, false));
  expected_local_alpha = std::min(
      expected_local_alpha,
      oracle_factor(base_rhoh - base_lower, full_rhoh - full_lower, false));
  expected_local_alpha = std::min(
      expected_local_alpha,
      oracle_factor(base_upper - base_rhoh, full_upper - full_rhoh, false));
  double expected_global_alpha = -1.0;
  oracle_valid =
      oracle_valid &&
      MPI_Allreduce(&expected_local_alpha, &expected_global_alpha, 1,
                    MPI_DOUBLE, MPI_MIN, communicator) == MPI_SUCCESS;

  std::uint64_t alpha_bits = 0U;
  std::memcpy(&alpha_bits, &call.diagnostics.source_endpoint_alpha,
              sizeof(alpha_bits));
  std::uint64_t alpha_min_bits = 0U;
  std::uint64_t alpha_max_bits = 0U;
  const bool alpha_rank_identical =
      MPI_Allreduce(&alpha_bits, &alpha_min_bits, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) == MPI_SUCCESS &&
      MPI_Allreduce(&alpha_bits, &alpha_max_bits, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) == MPI_SUCCESS &&
      alpha_min_bits == alpha_max_bits;

  bool passed = expect(
      status && call.certificate.valid() && !call.diagnostics.failure.valid,
      "finite BE source crossing publishes a valid predictor certificate");
  passed &= expect(
      oracle_valid && alpha_rank_identical &&
          call.diagnostics.low_state ==
              ThermophysicalLowStateKind::
                  bdf_local_donor_flux_source_limited &&
          call.diagnostics.source_endpoint_alpha >= 0.0 &&
          call.diagnostics.source_endpoint_alpha < 1.0 &&
          close(call.diagnostics.source_endpoint_alpha,
                expected_global_alpha, 2.0e-13) &&
          call.diagnostics.bdf_endpoint_alpha == 1.0 &&
          call.diagnostics.enthalpy_endpoint_alpha == 1.0 &&
          call.diagnostics.enthalpy_solve_calls == 0U &&
          call.diagnostics.low_order_halo_exchanges == 1U,
      "source crossing selects one rank-common source alpha without BDF or implicit-endpoint relabeling");

  const double alpha = call.diagnostics.source_endpoint_alpha;
  bool bundle_reconstructed = true;
  for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < fixture.patch.cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho = call.output.low_order_density_workspace.unchecked(
            cell, 0U);
        const double low_h =
            call.output.low_order_enthalpy_workspace.unchecked(cell, 0U);
        const double low_species =
            call.output.low_order_independent_species.data[0U].unchecked(
                cell, 0U);
        const double low_passive =
            call.output.low_order_passive_scalars.data[0U].unchecked(cell,
                                                                     0U);
        bundle_reconstructed =
            bundle_reconstructed && same_bits(rho, 1.0) &&
            close(rho * low_h,
                  base_rhoh + alpha * (full_rhoh - base_rhoh)) &&
            close(rho * low_species,
                  base_species + alpha * (full_species - base_species)) &&
            close(rho * low_passive,
                  base_passive + alpha * (full_passive - base_passive));
      }
    }
  }
  passed &= expect(
      bundle_reconstructed,
      "enthalpy, species, and passive source increments use the same affine alpha");

  const ConstFaceFluxView paired = as_const(call.output.paired_mass_flux);
  bool mass_pair_unchanged = true;
  for (const auto axis : {CartesianAxis::x, CartesianAxis::y,
                          CartesianAxis::z}) {
    const ConstFaceFieldView actual =
        axis == CartesianAxis::x
            ? paired.x
            : (axis == CartesianAxis::y ? paired.y : paired.z);
    const ConstFaceFieldView accepted =
        axis == CartesianAxis::x
            ? accepted_flux.x
            : (axis == CartesianAxis::y ? accepted_flux.y : accepted_flux.z);
    for (std::int32_t z = 0; z < actual.extents.z; ++z)
      for (std::int32_t y = 0; y < actual.extents.y; ++y)
        for (std::int32_t x = 0; x < actual.extents.x; ++x) {
          const Int3 face{x, y, z};
          mass_pair_unchanged =
              mass_pair_unchanged &&
              same_bits(actual.unchecked(face), accepted.unchecked(face));
        }
  }
  passed &= expect(
      mass_pair_unchanged && call.diagnostics.mass_flux_scale == 1.0,
      "source alpha does not scale density or the paired face mass flux");

  PredictorData zero_data(fixture.patch.cells, accepted_enthalpy,
                          accepted_enthalpy, -1000000.0, 0.0);
  fill(zero_data.species, 0.0);
  fill(zero_data.species_previous, 0.0);
  fill(zero_data.species_rhs, rank == size - 1 ? -0.5 : 0.5);
  fill(zero_data.passive, accepted_passive);
  fill(zero_data.passive_previous, accepted_passive);
  fill(zero_data.passive_rhs, passive_rhs);
  PredictorCall zero_call = make_call(
      fixture, zero_data, accepted_flux, previous_flux, exchange_plan, 994U);
  zero_call.input.bdf = {1.0 / kDt, -1.0 / kDt, 0.0, 1U};
  zero_call.input.density_previous = {};
  zero_call.input.enthalpy_previous = {};
  zero_call.input.mass_flux_previous = {};
  zero_call.input.enthalpy_ghosts.previous = {};
  zero_call.species_previous[0U] = {};
  zero_call.passive_previous[0U] = {};
  zero_call.species_ghosts[0U].previous = {};
  zero_call.passive_ghosts[0U].previous = {};
  zero_call.species_rates[0U].previous = {};
  zero_call.passive_rates[0U].previous = {};
  zero_call.input.enthalpy_nonadvective_rhs.previous = {};
  const Status zero_status =
      fixture.equations.thermophysical_predictor().predict(
          communicator, Status{}, zero_call.input, zero_call.output,
          zero_call.slow_path, zero_call.diagnostics, zero_call.certificate);
  bool zero_anchor_preserved = zero_status && zero_call.certificate.valid();
  for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < fixture.patch.cells.x; ++x) {
        const Int3 cell{x, y, z};
        zero_anchor_preserved =
            zero_anchor_preserved &&
            same_bits(zero_call.output.low_order_density_workspace.unchecked(
                          cell, 0U),
                      1.0) &&
            same_bits(zero_call.output.low_order_enthalpy_workspace.unchecked(
                          cell, 0U),
                      accepted_enthalpy) &&
            same_bits(zero_call.output.low_order_independent_species.data[0U]
                          .unchecked(cell, 0U),
                      0.0) &&
            same_bits(zero_call.output.low_order_passive_scalars.data[0U]
                          .unchecked(cell, 0U),
                      accepted_passive);
      }
    }
  }
  passed &= expect(
      zero_anchor_preserved && !zero_call.diagnostics.failure.valid &&
          zero_call.diagnostics.low_state ==
              ThermophysicalLowStateKind::
                  bdf_local_donor_flux_source_limited &&
          zero_call.diagnostics.source_endpoint_alpha == 0.0 &&
          zero_call.diagnostics.bdf_endpoint_alpha == 1.0,
      "zero is a valid global source alpha and preserves the source-free anchor bit-for-bit");
  if (rank == 0) {
    std::cerr << "source_bundle_oracle alpha=" << alpha
              << " expected=" << expected_global_alpha
              << " collectives="
              << call.diagnostics.blocking_collectives << '\n';
  }
  return passed;
}

bool test_implicit_upwind_endpoint_red(Fixture& fixture,
                                       std::uintptr_t exchange_plan,
                                       MPI_Comm communicator) {
  constexpr double divergence = 60.0;
  constexpr double donor_h = 210000.0;
  constexpr double receiver_h = 360000.0;
  // Pinned independent oracle from the exact existing BDF2/EX2 endpoint and
  // the frozen two-row implicit-upwind M-matrix.  The high donor is below the
  // 200 K mixture bound; the implicit endpoint is strictly admissible.
  constexpr double high_donor = 0x1.6be5e50d79436p+17;
  constexpr double implicit_donor = 0x1.77p+17;
  constexpr double implicit_receiver = 0x1.4dfcp+18;

  FluxHistory history;
  if (!expect(make_flux_history(fixture.patch.cells, history),
              "implicit-upwind RED allocates committed flux history")) {
    return false;
  }
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  int communicator_size = 1;
  if (MPI_Comm_size(communicator, &communicator_size) != MPI_SUCCESS ||
      communicator_size <= 0) {
    return false;
  }
  const double fixture_divergence =
      divergence / static_cast<double>(communicator_size);
  if (!expect(commit_shared_face_transfer(
                  fixture.equations.kernels(), fixture.patch.cells, history,
                  previous_flux, fixture_divergence) &&
                  commit_shared_face_transfer(
                      fixture.equations.kernels(), fixture.patch.cells,
                      history, accepted_flux, fixture_divergence),
              "implicit-upwind RED publishes frozen shared-face histories")) {
    return false;
  }

  PredictorData data(fixture.patch.cells, 300000.0, 300000.0, 0.0, 0.0);
  fill(data.rho, 1.0);
  fill(data.rho_previous, 1.0);
  fill(data.enthalpy, 300000.0);
  fill(data.enthalpy_previous, 300000.0);
  data.enthalpy.view.unchecked({0, 0, 0}, 0U) = donor_h;
  data.enthalpy_previous.view.unchecked({0, 0, 0}, 0U) = donor_h;
  data.enthalpy.view.unchecked({1, 0, 0}, 0U) = receiver_h;
  data.enthalpy_previous.view.unchecked({1, 0, 0}, 0U) = receiver_h;
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                 exchange_plan, 919U);
  call.input.dt = kSharedFaceDt;
  call.input.bdf = kSharedFaceBdf;
  const Status status = fixture.equations.thermophysical_predictor().predict(
      communicator, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  double lower = 0.0;
  double upper = 0.0;
  bool passed = expect(static_cast<bool>(status),
                       "implicit-upwind controlled predictor succeeds");
  passed &= expect(bound_for_state(fixture.thermodynamics, 0.76, 0.152,
                                   lower, upper),
                   "implicit-upwind oracle obtains the donor bounds");
  passed &= expect(0.76 * high_donor < lower &&
                       0.76 * implicit_donor > lower &&
                       implicit_donor == 192000.0 &&
                       implicit_receiver == 342000.0,
                   "explicit high violates while independent implicit M-matrix endpoint is admissible");
  passed &= expect(call.diagnostics.constraint ==
                       ThermophysicalAdmissibilityConstraint::enthalpy_lower,
                   "controlled route is selected only by the enthalpy lower bound");
  passed &= expect(call.diagnostics.low_state ==
                       ThermophysicalLowStateKind::implicit_upwind,
                   "enthalpy-only high failure selects implicit_upwind endpoint");
  const LinearSolveResult& solve = call.diagnostics.enthalpy_endpoint;
  passed &= expect(call.diagnostics.enthalpy_solve_calls == 1U &&
                       solve.status &&
                       solve.termination == LinearTermination::converged &&
                       solve.iterations == 2U && solve.reduction_calls == 7U &&
                       solve.operator_applies == 4U &&
                       solve.preconditioner_applies == 2U &&
                       solve.final_true_residual == 0.0,
                   "implicit endpoint publishes the exact converged FGMRES work");
  passed &= expect(
      close(call.output.low_order_enthalpy_workspace.unchecked({0, 0, 0}, 0U),
            implicit_donor) &&
          close(call.output.low_order_enthalpy_workspace.unchecked(
                    {1, 0, 0}, 0U),
                implicit_receiver) &&
          call.diagnostics.mass_flux_scale == 1.0 &&
          call.diagnostics.enthalpy_endpoint_alpha == 1.0 &&
          call.diagnostics.blocking_collectives == 13U &&
          call.diagnostics.low_order_transport_passes == 1U,
      "implicit endpoint matches the independent solution without mass-flux damping");
  return passed;
}

bool test_implicit_source_limited_endpoint_red(
    Fixture& fixture, std::uintptr_t exchange_plan,
    MPI_Comm communicator) {
  constexpr double divergence = 60.0;
  constexpr double donor_h = 210000.0;
  constexpr double receiver_h = 360000.0;
  constexpr double accepted_source = -2500000.0;
  constexpr double implicit_donor = 172000.0;

  FluxHistory history;
  if (!expect(make_flux_history(fixture.patch.cells, history),
              "source-limited RED allocates committed flux history")) {
    return false;
  }
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  int communicator_size = 1;
  if (MPI_Comm_size(communicator, &communicator_size) != MPI_SUCCESS ||
      communicator_size <= 0) {
    return false;
  }
  const double fixture_divergence =
      divergence / static_cast<double>(communicator_size);
  if (!expect(commit_shared_face_transfer(
                  fixture.equations.kernels(), fixture.patch.cells, history,
                  previous_flux, fixture_divergence) &&
                  commit_shared_face_transfer(
                      fixture.equations.kernels(), fixture.patch.cells,
                      history, accepted_flux, fixture_divergence),
              "source-limited RED publishes frozen shared-face histories")) {
    return false;
  }

  PredictorData data(fixture.patch.cells, 300000.0, 300000.0,
                     accepted_source, 0.0);
  fill(data.rho, 1.0);
  fill(data.rho_previous, 1.0);
  fill(data.enthalpy, 300000.0);
  fill(data.enthalpy_previous, 300000.0);
  data.enthalpy.view.unchecked({0, 0, 0}, 0U) = donor_h;
  data.enthalpy_previous.view.unchecked({0, 0, 0}, 0U) = donor_h;
  data.enthalpy.view.unchecked({1, 0, 0}, 0U) = receiver_h;
  data.enthalpy_previous.view.unchecked({1, 0, 0}, 0U) = receiver_h;
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                 exchange_plan, 929U);
  call.input.dt = kSharedFaceDt;
  call.input.bdf = kSharedFaceBdf;
  const Status status = fixture.equations.thermophysical_predictor().predict(
      communicator, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);

  const double rho = call.output.low_order_density_workspace.unchecked(
      {0, 0, 0}, 0U);
  double lower = 0.0;
  double upper = 0.0;
  bool passed = expect(bound_for_state(fixture.thermodynamics, rho,
                                       0.2 * rho, lower, upper),
                       "source-limited RED obtains the donor bounds");
  passed &= expect(call.diagnostics.enthalpy_solve_calls == 1U &&
                       call.diagnostics.enthalpy_endpoint.status &&
                       call.diagnostics.enthalpy_endpoint.termination ==
                           LinearTermination::converged &&
                       rho * donor_h >= lower && rho * donor_h <= upper &&
                       rho * implicit_donor < lower,
                   "source-limited RED reproduces a converged endpoint below the nonzero bound");
  constexpr double kEndpointSafety =
      1.0 - 64.0 * std::numeric_limits<double>::epsilon();
  const double maximum_alpha =
      (rho * donor_h - lower) / (rho * (donor_h - implicit_donor));
  const double expected_alpha =
      std::nextafter(maximum_alpha * kEndpointSafety, 0.0);
  const double expected_limited =
      donor_h + expected_alpha * (implicit_donor - donor_h);
  double minimum_alpha = call.diagnostics.enthalpy_endpoint_alpha;
  double maximum_alpha_observed = call.diagnostics.enthalpy_endpoint_alpha;
  if (MPI_Allreduce(MPI_IN_PLACE, &minimum_alpha, 1, MPI_DOUBLE, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, &maximum_alpha_observed, 1, MPI_DOUBLE,
                    MPI_MAX, communicator) != MPI_SUCCESS) {
    return false;
  }
  passed &= expect(
      status && call.diagnostics.low_state ==
                    ThermophysicalLowStateKind::implicit_upwind_source_limited &&
          std::isfinite(call.diagnostics.enthalpy_endpoint_alpha) &&
          call.diagnostics.enthalpy_endpoint_alpha >= 0.0 &&
          call.diagnostics.enthalpy_endpoint_alpha < 1.0 &&
          close(call.diagnostics.enthalpy_endpoint_alpha, expected_alpha) &&
          same_bits(minimum_alpha, maximum_alpha_observed) &&
          same_bits(call.diagnostics.enthalpy_endpoint_alpha,
                    minimum_alpha) &&
          call.diagnostics.mass_flux_scale == 1.0 &&
          call.diagnostics.blocking_collectives == 13U &&
          close(call.output.low_order_enthalpy_workspace.unchecked(
                    {0, 0, 0}, 0U),
                expected_limited) &&
          rho * call.output.low_order_enthalpy_workspace.unchecked(
                    {0, 0, 0}, 0U) >= lower,
      "inadmissible implicit endpoint is certified by one global source homotopy");
  return passed;
}

bool test_implicit_endpoint_inactive_face(Fixture& fixture,
                                          MPI_Comm communicator) {
  int communicator_size = 1;
  if (MPI_Comm_size(communicator, &communicator_size) != MPI_SUCCESS ||
      communicator_size <= 0) {
    return false;
  }
  FluxHistory history;
  if (!expect(make_flux_history(fixture.patch.cells, history),
              "inactive-face endpoint allocates committed flux history")) {
    return false;
  }
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  const double fixture_divergence =
      60.0 / static_cast<double>(communicator_size);
  if (!expect(commit_shared_face_transfer(
                  fixture.equations.kernels(), fixture.patch.cells, history,
                  previous_flux, fixture_divergence) &&
                  commit_shared_face_transfer(
                      fixture.equations.kernels(), fixture.patch.cells,
                      history, accepted_flux, fixture_divergence),
              "inactive-face endpoint publishes poisoned shared-face flux")) {
    return false;
  }

  PredictorData data(fixture.patch.cells, 300000.0, 300000.0, 0.0, 0.0);
  fill(data.density_workspace, 1.0);
  fill(data.enthalpy, 300000.0);
  fill(data.high_enthalpy, 300000.0);
  data.high_enthalpy.view.unchecked({0, 0, 0}, 0U) = 290000.0;
  ResourceCounters resources;
  const ConservativeEnthalpyEndpointInput input{
      kSharedFaceBdf,
      931U,
      as_const(data.density_workspace.view),
      as_const(data.enthalpy.view),
      as_const(data.high_enthalpy.view),
      accepted_flux,
      data.accepted_advection_workspace.view,
      data.previous_advection_workspace.view,
      data.low_enthalpy_workspace.view,
      &resources};
  LinearSolveResult solve;
  const Status status = fixture.endpoint.solve(input, solve);
  bool passed = expect(static_cast<bool>(status) && solve.status &&
                           (solve.termination == LinearTermination::converged ||
                            solve.termination == LinearTermination::zero_rhs),
                       "inactive-face implicit endpoint converges collectively");
  passed &= expect(
      close(data.low_enthalpy_workspace.view.unchecked({0, 0, 0}, 0U),
            290000.0) &&
          close(data.low_enthalpy_workspace.view.unchecked({1, 0, 0}, 0U),
                300000.0),
      "inactive IBM face contributes no diagonal or neighbour coupling");
  passed &= expect(resources.allocations == 0U &&
                       solve.operator_applies > 0U &&
                       solve.preconditioner_applies > 0U &&
                       solve.reduction_calls > 0U,
                   "inactive-face solve reports work with zero allocation");
  return passed;
}

bool test_limited_conserved_bundle(Fixture& fixture,
                                   ConstFaceFluxView accepted_flux,
                                   ConstFaceFluxView previous_flux,
                                   std::uintptr_t exchange_plan) {
  constexpr double accepted_h = 300000.0;
  constexpr double accepted_rhs = -900000.0;
  PredictorData data(fixture.patch.cells, accepted_h, accepted_h,
                     accepted_rhs, 0.0);
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                  exchange_plan, 902U);
  const Status status = fixture.equations.thermophysical_predictor().predict(
      MPI_COMM_SELF, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  bool passed = expect(static_cast<bool>(status),
                       "stage-15-style lower-bound predictor succeeds");
  passed &= expect(call.diagnostics.limited && call.diagnostics.theta < 1.0 &&
                       call.diagnostics.theta > 0.0,
                   "lower-bound violation activates a common theta");
  passed &= expect(!call.diagnostics.failure.valid,
                   "limited predictor has no failure provenance");
  passed &= expect(call.diagnostics.constraint ==
                       ThermophysicalAdmissibilityConstraint::enthalpy_lower,
                   "enthalpy lower bound identifies the limiting constraint");
  passed &= expect(call.diagnostics.low_state ==
                       ThermophysicalLowStateKind::bdf_accepted_rate &&
                       call.diagnostics.mass_flux_scale == 1.0 &&
                       call.diagnostics.bdf_endpoint_alpha == 1.0,
                   "limited path selects the accepted-rate BDF endpoint");
  passed &= expect(call.diagnostics.low_order_substeps == 1U &&
                       call.diagnostics.low_order_halo_exchanges == 0U,
                   "limited path records one fixed-step low-order update");
  passed &= expect(call.diagnostics.blocking_collectives > 1U &&
                       call.diagnostics.low_order_transport_passes > 0U,
                   "limited path records extra collectives and low-order transport work");
  if (!passed) {
    return false;
  }

  double lower = 0.0;
  double upper = 0.0;
  passed &= expect(bound_for_state(fixture.thermodynamics, 1.0, 0.2, lower,
                                   upper),
                   "NASA endpoint bounds are independently available");
  const double high_rho = high_rho_quantity(1.0, 1.0, 1.0, 1.0, 0.0, 0.0);
  const double high_rhoh = high_rho_quantity(
      1.0, accepted_h, 1.0, accepted_h, -accepted_rhs, 0.0);
  const double low_rho = bdf_accepted_rate_density(1.0, 1.0, 0.0);
  const double low_rhoh = bdf_accepted_rate_quantity(
      1.0, accepted_h, 1.0, accepted_h, 0.0, accepted_rhs);
  passed &= expect(high_rhoh < lower && low_rhoh >= lower,
                   "high BDF2/EX2 violates while accepted-rate low endpoint is admissible");
  const double actual_low_rho =
      call.output.low_order_density_workspace.unchecked({0, 0, 0}, 0U);
  const double actual_low_rhoh =
      actual_low_rho *
      call.output.low_order_enthalpy_workspace.unchecked({0, 0, 0}, 0U);
  passed &= expect(close(actual_low_rho, low_rho) &&
                       close(actual_low_rhoh, low_rhoh),
                   "Tier B endpoint matches the independent accepted-rate formula");
  const double low_enthalpy_margin = low_rhoh - lower;
  const double high_enthalpy_margin = high_rhoh - lower;
  const double theta_scale =
      std::max({1.0, std::abs(low_rhoh), std::abs(high_rhoh),
                std::abs(lower)});
  const double theta_target_margin =
      std::min(low_enthalpy_margin,
               512.0 * std::numeric_limits<double>::epsilon() * theta_scale);
  constexpr double theta_roundoff_safety =
      1.0 - 64.0 * std::numeric_limits<double>::epsilon();
  const double expected_theta = std::nextafter(
      ((low_enthalpy_margin - theta_target_margin) /
       (low_enthalpy_margin - high_enthalpy_margin)) *
          theta_roundoff_safety,
      0.0);
  passed &= expect(close(call.diagnostics.theta, expected_theta, 2.0e-11),
                   "global theta matches the independent conserved oracle");

  for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < fixture.patch.cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho =
            call.output.density_workspace.unchecked(cell, 0U);
        const double rhoh = rho * call.output.enthalpy.unchecked(cell, 0U);
        const double expected_rhoh_blend =
            low_rhoh + call.diagnostics.theta * (high_rhoh - low_rhoh);
        passed &= expect(close(rhoh, expected_rhoh_blend, 2.0e-11),
                         "final rho*h equals the conservative common blend");
        passed &= expect(rhoh >= lower -
                                  2.0e-11 * std::max(1.0, std::abs(lower)) &&
                             rhoh <= upper +
                                         2.0e-11 * std::max(1.0, std::abs(upper)),
                         "final conserved enthalpy is inside the 200 K bound");
        passed &= expect(close(rho, low_rho +
                                     call.diagnostics.theta *
                                         (high_rho - low_rho),
                                 2.0e-11),
                         "density uses the same global blend coefficient");
      }
    }
  }
  return passed;
}

bool test_serial_source_aware_low_endpoint(
    Fixture& fixture, ConstFaceFluxView accepted_flux,
    ConstFaceFluxView previous_flux, std::uintptr_t exchange_plan) {
  constexpr double accepted_h = 300000.0;
  constexpr double accepted_rhs = -4000000.0;
  PredictorData data(fixture.patch.cells, accepted_h, accepted_h,
                     accepted_rhs, 0.0);
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                 exchange_plan, 905U);
  const Status status = fixture.equations.thermophysical_predictor().predict(
      MPI_COMM_SELF, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  const ThermophysicalPredictorFailure& failure = call.diagnostics.failure;
  bool passed = expect(static_cast<bool>(status),
                       "source-aware low endpoint succeeds");
  if (!passed) return false;
  passed &= expect(!failure.valid,
                   "source-aware low endpoint has no failure provenance");
  passed &= expect(call.diagnostics.low_state ==
                           ThermophysicalLowStateKind::
                               bdf_accepted_rate_homotopy &&
                       call.diagnostics.mass_flux_scale == 1.0,
                   "source-aware low endpoint publishes BDF homotopy state");

  constexpr double accepted_density = 1.0;
  constexpr double accepted_species_density = 0.2;
  const double accepted_conserved_enthalpy = accepted_density * accepted_h;
  const double bdf_density =
      (-kBdf.a1 * accepted_density - kBdf.a2 * accepted_density) /
      kBdf.a0;
  const double bdf_conserved_enthalpy =
      (-kBdf.a1 * accepted_conserved_enthalpy -
       kBdf.a2 * accepted_conserved_enthalpy - (0.0 - accepted_rhs)) /
      kBdf.a0;
  const double bdf_species_density = accepted_species_density;
  double lower = 0.0;
  double upper = 0.0;
  passed &= expect(bound_for_state(fixture.thermodynamics, accepted_density,
                                   accepted_species_density, lower,
                                   upper),
                   "source-aware oracle obtains independent thermodynamic bounds");
  passed &= expect(accepted_conserved_enthalpy > lower &&
                       bdf_conserved_enthalpy < lower,
                   "accepted enthalpy is admissible while BDF crosses the lower bound");
  const double accepted_margin = accepted_conserved_enthalpy - lower;
  const double trial_margin = bdf_conserved_enthalpy - lower;
  const double raw_bdf_alpha = accepted_margin /
                               (accepted_margin - trial_margin);
  const double bdf_enthalpy_scale =
      std::max({1.0, std::abs(accepted_conserved_enthalpy),
                std::abs(bdf_conserved_enthalpy), std::abs(lower)});
  const double bdf_target_margin =
      std::min(accepted_margin,
               512.0 * std::numeric_limits<double>::epsilon() *
                   bdf_enthalpy_scale);
  const double guarded_bdf_alpha =
      (accepted_margin - bdf_target_margin) /
      (accepted_margin - trial_margin);
  constexpr double bdf_roundoff_safety =
      1.0 - 64.0 * std::numeric_limits<double>::epsilon();
  double expected_bdf_alpha = std::clamp(guarded_bdf_alpha, 0.0, 1.0);
  if (expected_bdf_alpha < 1.0) {
    expected_bdf_alpha =
        std::nextafter(expected_bdf_alpha * bdf_roundoff_safety, 0.0);
  }
  const double expected_low_density =
      accepted_density +
      expected_bdf_alpha * (bdf_density - accepted_density);
  const double expected_low_species_density =
      accepted_species_density +
      expected_bdf_alpha *
          (bdf_species_density - accepted_species_density);
  const double expected_low_conserved_enthalpy =
      accepted_conserved_enthalpy +
      expected_bdf_alpha *
          (bdf_conserved_enthalpy - accepted_conserved_enthalpy);
  const double expected_low_enthalpy =
      expected_low_conserved_enthalpy / expected_low_density;
  passed &= expect(std::isfinite(raw_bdf_alpha) && raw_bdf_alpha > 0.0 &&
                       raw_bdf_alpha < 1.0 && expected_bdf_alpha > 0.0 &&
                       expected_bdf_alpha < 1.0 &&
                       close(call.diagnostics.bdf_endpoint_alpha,
                             expected_bdf_alpha, 2.0e-13) &&
                       call.diagnostics.mass_flux_scale == 1.0,
                   "source-aware oracle matches the complete BDF alpha");
  const double actual_low_density =
      call.output.low_order_density_workspace.unchecked({0, 0, 0}, 0U);
  const double actual_low_enthalpy =
      call.output.low_order_enthalpy_workspace.unchecked({0, 0, 0}, 0U);
  const double actual_low_species =
      call.output.low_order_independent_species.data[0U]
          .unchecked({0, 0, 0}, 0U);
  passed &= expect(close(actual_low_density, expected_low_density) &&
                       close(actual_low_enthalpy, expected_low_enthalpy) &&
                       close(actual_low_species,
                             expected_low_species_density /
                                 expected_low_density),
                   "source-aware low endpoint equals the affine conserved blend");
  const double actual_low_species_density =
      actual_low_density * actual_low_species;
  double low_lower = 0.0;
  double low_upper = 0.0;
  passed &= expect(
      bound_for_state(fixture.thermodynamics, actual_low_density,
                      actual_low_species_density, low_lower, low_upper) &&
          actual_low_density > 0.0 && actual_low_species_density > 0.0 &&
          actual_low_species_density < actual_low_density &&
          actual_low_density * actual_low_enthalpy >=
              low_lower - 2.0e-11 * std::max(1.0, std::abs(low_lower)) &&
          actual_low_density * actual_low_enthalpy <=
              low_upper + 2.0e-11 * std::max(1.0, std::abs(low_upper)),
      "source-aware low endpoint remains inside thermodynamic bounds");
  const double represented_low_conserved_enthalpy =
      actual_low_density * actual_low_enthalpy;
  const double representation_scale =
      std::max({1.0, std::abs(represented_low_conserved_enthalpy),
                std::abs(low_lower), std::abs(bdf_conserved_enthalpy),
                std::abs(accepted_conserved_enthalpy)});
  const double representation_guard =
      512.0 * std::numeric_limits<double>::epsilon() *
      representation_scale;
  passed &= expect(
      represented_low_conserved_enthalpy - low_lower >=
          representation_guard,
      "BDF homotopy reserves the represented enthalpy margin");
  const double final_density =
      call.output.density_workspace.unchecked({0, 0, 0}, 0U);
  const double final_enthalpy =
      call.output.enthalpy.unchecked({0, 0, 0}, 0U);
  const double final_species =
      call.output.independent_species.data[0U].unchecked({0, 0, 0}, 0U);
  const double final_species_density = final_density * final_species;
  double final_lower = 0.0;
  double final_upper = 0.0;
  passed &= expect(
      bound_for_state(fixture.thermodynamics, final_density,
                      final_species_density, final_lower, final_upper) &&
          std::isfinite(final_density) && final_density > 0.0 &&
          std::isfinite(final_species) && final_species > 0.0 &&
          final_species_density < final_density &&
          std::isfinite(final_enthalpy) &&
          final_density * final_enthalpy >=
              final_lower - 2.0e-11 * std::max(1.0, std::abs(final_lower)) &&
          final_density * final_enthalpy <=
              final_upper + 2.0e-11 * std::max(1.0, std::abs(final_upper)),
      "source-aware final tuple is admissible");
  std::cerr << "source_aware_idp_oracle accepted_rhoh="
            << accepted_conserved_enthalpy
            << " bdf_rhoh=" << bdf_conserved_enthalpy
            << " lower=" << lower << " raw_alpha=" << raw_bdf_alpha
            << " alpha=" << expected_bdf_alpha << '\n';
  return passed;
}

bool test_mpi_zero_alpha_homotopy(
    Fixture& fixture, ConstFaceFluxView accepted_flux,
    ConstFaceFluxView previous_flux, std::uintptr_t exchange_plan, int rank,
    int size) {
  if (size != 2) return true;
  double independent_min = 0.0;
  double independent_max = 0.0;
  double dependent_min = 0.0;
  double dependent_max = 0.0;
  if (!expect(
          static_cast<bool>(fixture.thermodynamics
                                .independent_species_enthalpy_bounds(
                                    0U, independent_min, independent_max)) &&
              static_cast<bool>(fixture.thermodynamics.species_enthalpy_bounds(
                  fixture.thermodynamics.dependent_species_index(),
                  dependent_min, dependent_max)),
          "MPI failure fixture obtains composition enthalpy bounds")) {
    return false;
  }
  // Make the limiting accepted tuple exactly meet the affine lower bound.
  // A finite negative source then selects the accepted anchor itself.  Zero
  // is a valid global homotopy coefficient and must be bit-identical on every
  // rank rather than being misclassified as a candidate failure.
  const double lower_bound_enthalpy =
      dependent_min + 0.2 * (independent_min - dependent_min);
  PredictorData one_rank_failure_data(
      fixture.patch.cells,
      rank == 1 ? lower_bound_enthalpy : 300000.0,
      rank == 1 ? lower_bound_enthalpy : 300000.0,
      rank == 1 ? -4000000.0 : -900000.0, 0.0);
  PredictorCall one_rank_call = make_call(
      fixture, one_rank_failure_data, accepted_flux, previous_flux,
      exchange_plan,
      906U);
  const Status one_rank_status =
      fixture.equations.thermophysical_predictor().predict(
          MPI_COMM_WORLD, Status{}, one_rank_call.input, one_rank_call.output,
          one_rank_call.slow_path, one_rank_call.diagnostics,
          one_rank_call.certificate);
  bool passed = expect(
      one_rank_status && one_rank_call.certificate.valid() &&
          !one_rank_call.diagnostics.failure.valid &&
          one_rank_call.diagnostics.low_state ==
              ThermophysicalLowStateKind::bdf_accepted_rate_homotopy &&
          one_rank_call.diagnostics.bdf_endpoint_alpha == 0.0 &&
          one_rank_call.diagnostics.mass_flux_scale == 1.0 &&
          one_rank_call.diagnostics.theta == 0.0 &&
          one_rank_call.diagnostics.limiting_rank == 1,
      "two-rank one-sided boundary state selects the common accepted anchor");
  bool bitwise_anchor_preserved = true;
  for (std::int32_t z = 0; z < fixture.patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < fixture.patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < fixture.patch.cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double accepted_density =
            one_rank_failure_data.rho.view.unchecked(cell, 0U);
        const double accepted_enthalpy =
            one_rank_failure_data.enthalpy.view.unchecked(cell, 0U);
        const double accepted_species =
            one_rank_failure_data.species.view.unchecked(cell, 0U);
        bitwise_anchor_preserved =
            bitwise_anchor_preserved &&
            same_bits(one_rank_call.output.low_order_density_workspace
                          .unchecked(cell, 0U),
                      accepted_density) &&
            same_bits(one_rank_call.output.low_order_enthalpy_workspace
                          .unchecked(cell, 0U),
                      accepted_enthalpy) &&
            same_bits(one_rank_call.output.low_order_independent_species
                          .data[0U]
                          .unchecked(cell, 0U),
                      accepted_species) &&
            same_bits(one_rank_call.output.density_workspace.unchecked(cell,
                                                                       0U),
                      accepted_density) &&
            same_bits(one_rank_call.output.enthalpy.unchecked(cell, 0U),
                      accepted_enthalpy) &&
            same_bits(one_rank_call.output.independent_species.data[0U]
                          .unchecked(cell, 0U),
                      accepted_species);
      }
    }
  }
  passed &= expect(
      bitwise_anchor_preserved,
      "zero alpha and theta preserve the accepted tuple bit-for-bit");
  std::array<int, 1U> local_pass{{passed ? 1 : 0}};
  int global_pass = 0;
  MPI_Allreduce(local_pass.data(), &global_pass, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
  if (global_pass == 0) passed = false;
  double one_rank_min = 1.0;
  double one_rank_max = -1.0;
  MPI_Allreduce(&one_rank_call.diagnostics.bdf_endpoint_alpha, &one_rank_min,
                1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&one_rank_call.diagnostics.bdf_endpoint_alpha, &one_rank_max,
                1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  passed &= expect(one_rank_min == 0.0 && one_rank_max == 0.0,
                   "two-rank one-sided homotopy alpha is rank-identical");

  PredictorData invalid_anchor_data(fixture.patch.cells, 300000.0, 300000.0,
                                    -900000.0, 0.0);
  fill(invalid_anchor_data.species, -0.01);
  fill(invalid_anchor_data.species_previous, -0.01);
  PredictorCall invalid_anchor_call = make_call(
      fixture, invalid_anchor_data, accepted_flux, previous_flux,
      exchange_plan, 907U);
  const Status invalid_anchor_status =
      fixture.equations.thermophysical_predictor().predict(
          MPI_COMM_WORLD, Status{}, invalid_anchor_call.input,
          invalid_anchor_call.output, invalid_anchor_call.slow_path,
          invalid_anchor_call.diagnostics, invalid_anchor_call.certificate);
  passed &= expect(
      invalid_anchor_status.code == StatusCode::numerical_failure &&
          invalid_anchor_status.detail == 985U &&
          invalid_anchor_call.diagnostics.failure.valid &&
          invalid_anchor_call.diagnostics.failure.rank == 0 &&
          invalid_anchor_call.diagnostics.failure.reason ==
              ThermophysicalPredictorFailureReason::low_tuple_admissibility,
      "invalid accepted species anchor is a collective numerical failure");
  return passed;
}

bool test_zero_high_density_uses_conserved_bundle(
    Fixture& fixture, ConstFaceFluxView accepted_flux,
    ConstFaceFluxView previous_flux, std::uintptr_t exchange_plan) {
  PredictorData data(fixture.patch.cells, 300000.0, 300000.0, 0.0, 0.0);
  // With the frozen constant-step BDF2 coefficients, rho_n=1 and
  // rho_nm1=4 produce an exactly zero high density.  The accepted-state low
  // endpoint remains rho=1 and must therefore recover the complete bundle.
  fill(data.rho_previous, 4.0);
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                  exchange_plan, 904U);
  const Status status = fixture.equations.thermophysical_predictor().predict(
      MPI_COMM_SELF, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  if (!status) {
    std::cerr << "zero-density status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " rho="
              << call.output.density_workspace.unchecked({0, 0, 0}, 0U)
              << " h=" << call.output.enthalpy.unchecked({0, 0, 0}, 0U)
              << " Y="
              << call.output.independent_species.data[0U].unchecked(
                     {0, 0, 0}, 0U)
              << " theta=" << call.diagnostics.theta << '\n';
  }
  bool passed = expect(static_cast<bool>(status),
                       "zero high density is recovered by the IDP bundle");
  passed &= expect(call.diagnostics.limited &&
                       call.diagnostics.theta > 0.0 &&
                       call.diagnostics.theta < 1.0 &&
                       call.diagnostics.low_state ==
                           ThermophysicalLowStateKind::
                               bdf_accepted_rate_homotopy &&
                       call.diagnostics.bdf_endpoint_alpha > 0.0 &&
                       call.diagnostics.bdf_endpoint_alpha < 1.0 &&
                       call.diagnostics.mass_flux_scale == 1.0 &&
                       call.diagnostics.constraint ==
                           ThermophysicalAdmissibilityConstraint::
                               independent_species,
                   "guarded zero-density bundle selects an inward common theta");
  passed &= expect(!call.diagnostics.failure.valid,
                   "zero-density recovery has no failure provenance");
  for (std::int32_t z = 0; z < fixture.patch.cells.z && passed; ++z) {
    for (std::int32_t y = 0; y < fixture.patch.cells.y && passed; ++y) {
      for (std::int32_t x = 0; x < fixture.patch.cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho =
            call.output.density_workspace.unchecked(cell, 0U);
        const double h = call.output.enthalpy.unchecked(cell, 0U);
        const double species =
            call.output.independent_species.data[0U].unchecked(cell, 0U);
        passed &= expect(std::isfinite(rho) && rho > 0.0 &&
                             std::isfinite(h) && std::isfinite(species) &&
                             species >= 0.0 && species <= 1.0,
                         "recovered zero-density bundle is finite and admissible");
      }
    }
  }
  return passed;
}

bool test_stale_previous_ghost_is_atomic(
    Fixture& fixture, ConstFaceFluxView accepted_flux,
    ConstFaceFluxView previous_flux, std::uintptr_t exchange_plan,
    MPI_Comm communicator, int rank, int size) {
  PredictorData data(fixture.patch.cells, 300000.0, 299000.0, 0.0, 0.0);
  PredictorCall call = make_call(fixture, data, accepted_flux, previous_flux,
                                  exchange_plan, 903U);
  const Status first = fixture.equations.thermophysical_predictor().predict(
      communicator, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  if (!expect(static_cast<bool>(first),
              "authority fixture first predictor publishes")) {
    return false;
  }
  const ThermophysicalPredictorCertificate marker = call.certificate;
  const std::vector<double> enthalpy_snapshot = data.high_enthalpy.storage;
  if (size == 1 || rank == size - 1) {
    call.input.enthalpy_ghosts.previous.state += 1U;
  }
  const Status stale = fixture.equations.thermophysical_predictor().predict(
      communicator, Status{}, call.input, call.output, call.slow_path,
      call.diagnostics, call.certificate);
  return expect(stale.code == StatusCode::invalid_plan,
                "stale previous ghost authority is rejected") &&
         expect(same_certificate(call.certificate, marker),
                "stale authority cannot publish a replacement certificate") &&
         expect(data.high_enthalpy.storage == enthalpy_snapshot,
                "stale authority rejects before trial output writes");
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 1;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Fixture fixture;
  bool passed = expect(make_fixture(fixture), "IDP predictor fixture compiles");
  if (passed) {
    FluxHistory history;
    passed &= expect(make_flux_history(fixture.patch.cells, history),
                     "IDP committed-flux history allocates");
    ConstFaceFluxView previous_flux;
    ConstFaceFluxView accepted_flux;
    if (passed) {
      passed &= expect(commit_zero_flux(fixture.equations.kernels(),
                                        fixture.patch.cells, history,
                                        previous_flux) &&
                           commit_zero_flux(fixture.equations.kernels(),
                                            fixture.patch.cells, history,
                                            accepted_flux),
                       "two committed zero-flux histories publish");
    }
    if (passed) {
      const std::uintptr_t exchange_plan =
          static_cast<std::uintptr_t>(0x1d7002U);
      passed &= test_fast_path_bit_identity(fixture, accepted_flux,
                                            previous_flux, exchange_plan);
      passed &= test_small_dt_paired_continuity_roundoff(
          fixture, accepted_flux, previous_flux, exchange_plan);
      passed &= test_local_donor_paired_flux(fixture, exchange_plan,
                                             MPI_COMM_SELF);
      passed &= test_invalid_zero_outgoing_base(
          fixture, accepted_flux, previous_flux, exchange_plan,
          MPI_COMM_SELF);
      passed &= test_stale_previous_ghost_is_atomic(
          fixture, accepted_flux, previous_flux, exchange_plan,
          MPI_COMM_WORLD, rank, size);
    }
  }

  Fixture collective_fixture;
  passed &= expect(
      make_fixture(collective_fixture, MPI_COMM_WORLD, {2 * size, 2, 2}) &&
          collective_fixture.patch.cells.x == 2 &&
          collective_fixture.patch.cells.y == 2 &&
          collective_fixture.patch.cells.z == 2,
      "collective implicit-upwind fixture compiles with two local x cells");
  if (passed)
    passed &= test_local_donor_paired_flux(
        collective_fixture, static_cast<std::uintptr_t>(0x1d7012U),
        MPI_COMM_WORLD);
  if (passed)
    passed &= test_complete_source_bundle_globalization(
        collective_fixture, static_cast<std::uintptr_t>(0x1d7013U),
        MPI_COMM_WORLD, rank, size);
  if (passed) {
    FluxHistory collective_zero_history;
    ConstFaceFluxView collective_previous_flux;
    ConstFaceFluxView collective_accepted_flux;
    passed &= expect(
        make_flux_history(collective_fixture.patch.cells,
                          collective_zero_history) &&
            commit_zero_flux(collective_fixture.equations.kernels(),
                             collective_fixture.patch.cells,
                             collective_zero_history,
                             collective_previous_flux) &&
            commit_zero_flux(collective_fixture.equations.kernels(),
                             collective_fixture.patch.cells,
                             collective_zero_history,
                             collective_accepted_flux),
        "collective reason-17 fixture publishes zero flux histories");
    if (passed)
      passed &= test_invalid_zero_outgoing_base(
          collective_fixture, collective_accepted_flux,
          collective_previous_flux,
          static_cast<std::uintptr_t>(0x1d7014U), MPI_COMM_WORLD);
  }

  Fixture physical_fixture;
  passed &= expect(make_fixture(physical_fixture, MPI_COMM_SELF, {2, 2, 2},
                                false, true),
                   "physical donor-factor fixture compiles");
  if (passed)
    passed &= test_physical_inflow_outflow_factors(
        physical_fixture, static_cast<std::uintptr_t>(0x1d7022U));

  int local = passed ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (rank == 0) {
    std::cout << (global != 0 ? "PASS" : "FAIL")
              << " v04_solver_thermophysical_idp\n";
  }
  MPI_Finalize();
  return global == 0 ? 1 : 0;
}
