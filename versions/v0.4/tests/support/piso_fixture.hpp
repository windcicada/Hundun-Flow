// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_flow.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

namespace hundun::v04::test {

inline PlanFingerprint exact_composition_identity_for_test(
    PlanFingerprint thermodynamics,
    Span<const ConstFieldView> independent_species, Int3 cells) noexcept {
  constexpr std::uint64_t offset = UINT64_C(1469598103934665603);
  constexpr std::uint64_t prime = UINT64_C(1099511628211);
  const auto mix = [=](std::uint64_t hash, std::uint64_t value) noexcept {
    hash ^= value;
    hash *= prime;
    return hash;
  };
  std::uint64_t hash = mix(offset, UINT64_C(0x706563636f6d706f));
  hash = mix(hash, thermodynamics);
  hash = mix(hash, static_cast<std::uint64_t>(independent_species.size));
  for (std::size_t species = 0U; species < independent_species.size;
       ++species) {
    const ConstFieldView field = independent_species.data[species];
    hash = mix(hash, field.components);
    for (std::uint8_t component = 0U; component < field.components;
         ++component)
      for (std::int32_t z = 0; z < cells.z; ++z)
        for (std::int32_t y = 0; y < cells.y; ++y)
          for (std::int32_t x = 0; x < cells.x; ++x) {
            std::uint64_t bits = 0U;
            const double value = field.unchecked({x, y, z}, component);
            std::memcpy(&bits, &value, sizeof(bits));
            hash = mix(hash, bits);
          }
  }
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

struct OwnedFaceField {
  std::vector<double> storage;
  FaceFieldView view{};
};

inline OwnedField make_field(FieldId field, Int3 cells,
                             std::uint8_t components,
                             std::uint8_t ghosts,
                             RevisionToken revision,
                             StorageIdentity storage_identity) {
  OwnedField field_storage;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  field_storage.storage.assign(nx * ny * nz * components, 0.0);
  field_storage.view.base = field_storage.storage.data() + ghosts + ghosts * nx +
                            ghosts * nx * ny;
  field_storage.view.interior = cells;
  field_storage.view.ghosts = {ghosts, ghosts, ghosts};
  field_storage.view.components = components;
  field_storage.view.stride_y = nx;
  field_storage.view.stride_z = nx * ny;
  field_storage.view.component_stride = nx * ny * nz;
  field_storage.view.field = field;
  field_storage.view.revision = revision;
  field_storage.view.storage_identity = storage_identity;
  field_storage.view.revision_domain = 15001U;
  return field_storage;
}

inline OwnedFaceField make_face_field(CartesianAxis axis, Int3 cells,
                                      StorageIdentity storage_identity) {
  OwnedFaceField field_storage;
  Int3 extents = cells;
  if (axis == CartesianAxis::x) {
    ++extents.x;
  } else if (axis == CartesianAxis::y) {
    ++extents.y;
  } else {
    ++extents.z;
  }
  field_storage.storage.assign(static_cast<std::size_t>(extents.x) * extents.y *
                                   extents.z,
                               0.0);
  field_storage.view = {
      field_storage.storage.data(),
      extents,
      static_cast<std::size_t>(extents.x),
      static_cast<std::size_t>(extents.x) * extents.y,
      axis,
      storage_identity,
      15002U};
  return field_storage;
}

inline void fill(OwnedField& field, double value) {
  std::fill(field.storage.begin(), field.storage.end(), value);
}

inline SpeciesThermophysicalSpec test_air() {
  SpeciesThermophysicalSpec species;
  species.stable_name = "air";
  species.molecular_weight = 28.96546;
  species.temperature_switch = 1000.0;
  species.nasa7_low[0U] = 3.5;
  species.nasa7_high[0U] = 3.5;
  species.viscosity_reference = 1.8e-5;
  species.conductivity = 0.026;
  return species;
}

inline ThermophysicalSpec test_thermophysical_spec() {
  ThermophysicalSpec spec;
  spec.data_file = "analytic.d";
  spec.minimum_temperature = 200.0;
  spec.maximum_temperature = 2000.0;
  spec.temperature_relative_tolerance = 1.0e-12;
  spec.maximum_temperature_iterations = 64U;
  spec.closed_mass_relative_tolerance = 1.0e-12;
  spec.maximum_closed_mass_iterations = 32U;
  spec.maximum_closed_mass_relative_step = 0.2;
  spec.species.push_back(test_air());
  return spec;
}

inline PisoPlanSpec test_piso_spec() {
  PisoPlanSpec spec;
  spec.pressure_correctors = 2U;
  spec.pressure_stage = 21U;
  spec.final_flux_slot = 0U;
  spec.pressure_solve = {1.0e-15, 1.0e-13, 400U, 4U, 12U};
  spec.eos_tolerance = 1.0e-10;
  spec.continuity_tolerance = 1.0e-10;
  spec.closed_mass_tolerance = 1.0e-10;
  spec.gauge_tolerance = 1.0e-12;
  return spec;
}

class PeriodicPisoFixture {
 public:
  bool initialize(std::int32_t cell_count,
                  MPI_Comm communicator = MPI_COMM_SELF,
                  bool multispecies = false) {
    if (cell_count < 2) {
      return false;
    }
    CartesianMeshSpec mesh;
    mesh.kind = GeometryKind::uniform;
    mesh.lower = {0.0, 0.0, 0.0};
    mesh.upper = {1.0, 1.0, 1.0};
    mesh.has_exact_cells = true;
    mesh.exact_cells = {cell_count, cell_count, cell_count};
    const double spacing = 1.0 / static_cast<double>(cell_count);
    mesh.minimum_spacing = {spacing, spacing, spacing};
    mesh.max_growth_ratio = 1.0;
    const std::uint64_t global_cells =
        static_cast<std::uint64_t>(cell_count) * cell_count * cell_count;
    mesh.limits.max_global_cells = global_cells;
    mesh.limits.max_memory_bytes_per_rank = 1U << 29U;

    ValidatedModel model;
    model.mesh = mesh;
    model.fingerprint = 0x15010001U +
                        static_cast<PlanFingerprint>(cell_count);
    model.pressure_reference = PressureReferenceKind::closed_mass;
    for (BoundaryFaceSpec& face : model.boundaries) {
      face.flow_kind = BoundaryKind::periodic;
      face.thermal_kind = BoundaryKind::none;
    }
    if (multispecies) {
      model.transported_scalars.push_back(
          {"sp0", TransportedScalarRole::species, 1.0, 1.0});
      model.transported_scalars.push_back(
          {"sp1", TransportedScalarRole::species, 1.0, 1.0});
    }
    model.schemes.momentum = ConvectionScheme::central2;
    model.schemes.enthalpy = ConvectionScheme::central2;
    model.schemes.species = ConvectionScheme::central2;
    model.schemes.passive_scalar = ConvectionScheme::central2;

    FieldRegistry registry;
    FieldId id = 0U;
    if (!registry.require_field("rho", 1U, 2U, id) || id != 0U ||
        !registry.require_field("U", 3U, 2U, id) || id != 1U ||
        !registry.require_field("pi", 1U, 2U, id) || id != 2U ||
        !registry.require_field("h", 1U, 2U, id) || id != 3U ||
        !registry.require_field("T", 1U, 2U, id) || id != 4U ||
        !registry.require_field("mu_eff", 1U, 2U, id) || id != 5U ||
        !registry.require_field("drho_dp", 1U, 2U, id) || id != 6U ||
        !registry.require_field("gradU", 9U, 2U, id) || id != 7U ||
        !CartesianGeometryCompiler::compile(communicator, mesh, {}, geometry,
                                            patch) ||
        !BoundaryCompiler::compile(communicator, model, geometry, patch,
                                   registry, boundary, schemes, time)) {
      return false;
    }
    if (multispecies &&
        (!registry.require_field("sp0", 1U, 2U, id) || id != 8U ||
         !registry.require_field("sp1", 1U, 2U, id) || id != 9U)) {
      return false;
    }
    ThermophysicalSpec thermophysics = test_thermophysical_spec();
    if (multispecies) {
      thermophysics.species.clear();
      for (std::size_t species_index = 0U; species_index < 3U;
           ++species_index) {
        SpeciesThermophysicalSpec species = test_air();
        species.stable_name = "sp" + std::to_string(species_index);
        species.molecular_weight += static_cast<double>(species_index) * 2.0;
        thermophysics.species.push_back(std::move(species));
      }
    }
    const Span<const TransportedScalarSpec> transported =
        multispecies
            ? Span<const TransportedScalarSpec>{
                  model.transported_scalars.data(),
                  model.transported_scalars.size()}
            : Span<const TransportedScalarSpec>{};
    if (!ThermodynamicsPlan::compile(thermophysics, transported,
                                     thermodynamics) ||
        !TransportPlan::compile(thermophysics, thermodynamics, transport)) {
      return false;
    }
    std::vector<FieldId> declared{0U, 1U, 2U, 3U, 4U, 5U, 6U, 7U};
    if (multispecies) {
      declared.push_back(8U);
      declared.push_back(9U);
    }
    if (!contributions.configure({declared.data(), declared.size()}) ||
        !contributions.freeze()) {
      return false;
    }
    const std::array<ScalarEquationSpec, 2U> scalar_equations{{
        {8U, TransportedScalarRole::species, 1.0, 1.0},
        {9U, TransportedScalarRole::species, 1.0, 1.0},
    }};
    EquationPlanSpec equation_spec;
    equation_spec.density = 0U;
    equation_spec.velocity = 1U;
    equation_spec.pressure_perturbation = 2U;
    equation_spec.enthalpy = 3U;
    equation_spec.temperature = 4U;
    equation_spec.effective_viscosity = 5U;
    equation_spec.pressure_compressibility = 6U;
    equation_spec.velocity_gradient = 7U;
    equation_spec.pressure_reference = PressureReferenceKind::closed_mass;
    if (multispecies) {
      equation_spec.scalars =
          {scalar_equations.data(), scalar_equations.size()};
    }
    equation_spec.closed_mass_service_stage = 1U;
    equation_spec.maximum_cells_per_rank =
        static_cast<std::size_t>(global_cells);
    if (!EquationPlanSet::compile(communicator, schemes, geometry, patch,
                                  boundary, contributions, thermodynamics,
                                  transport, equation_spec, equations) ||
        !PisoPlan::compile(communicator, equations, test_piso_spec(), piso)) {
      return false;
    }

    const Int3 cells = patch.cells;
    const std::uint8_t reach = equations.kernels().reach();
    density = make_field(0U, cells, 1U, reach, 1601U, 2601U);
    velocity = make_field(1U, cells, 3U, 0U, 1602U, 2602U);
    momentum_diagonal = make_field(30U, cells, 3U, 0U, 1603U, 2603U);
    momentum_rhs = make_field(31U, cells, 3U, 0U, 1604U, 2604U);
    r_au = make_field(40U, cells, 3U, 1U, 1605U, 2605U);
    h_by_a = make_field(41U, cells, 3U, reach, 1606U, 2606U);
    pressure_gradient = make_field(42U, cells, 3U, 0U, 1607U, 2610U);
    x_pressure_coefficient =
        make_face_field(CartesianAxis::x, cells, 2607U);
    y_pressure_coefficient =
        make_face_field(CartesianAxis::y, cells, 2608U);
    z_pressure_coefficient =
        make_face_field(CartesianAxis::z, cells, 2609U);
    if (!FaceFluxStorage::allocate_workspace(cells, 2U, phi_storage) ||
        !phi_storage.workspace_view(0U, 1610U, phi_h_by_a) ||
        !phi_storage.workspace_view(1U, 1702U, trial_flux)) {
      return false;
    }
    fill(density, 1.0);
    fill(velocity, 0.0);
    fill(momentum_diagonal, 1.0);
    fill(momentum_rhs, 0.0);

    const std::array<HaloFieldSpec, 3U> halo_fields{{
        {density.view.field, 1U, 1U},
        {r_au.view.field, 1U, 3U},
        {h_by_a.view.field, reach, 3U}}};
    if (!halo.reserve(communicator, patch,
                      {halo_fields.data(), halo_fields.size()},
                      boundary.halo_topology())) {
      return false;
    }
    constexpr FieldId correction_field = 90U;
    const std::array<HaloFieldSpec, 1U> correction_halo_fields{{
        {correction_field, 1U, 1U}}};
    if (!correction_halo.reserve(
            communicator, patch,
            {correction_halo_fields.data(), correction_halo_fields.size()},
            boundary.halo_topology())) {
      return false;
    }
    const PisoCouplerWorkspace workspace{
        r_au.view,
        h_by_a.view,
        pressure_gradient.view,
        x_pressure_coefficient.view,
        y_pressure_coefficient.view,
        z_pressure_coefficient.view,
        phi_h_by_a};
    const PisoCouplerServices services{
        communicator, &geometry, patch, &boundary, &thermodynamics, &halo,
        160U,
        density.view.field, &correction_halo, 162U, correction_field};
    return static_cast<bool>(PressureVelocityCoupler::bind(
               piso, equations, services, workspace, coupler)) &&
           initialize_committed_flux_history(communicator);
  }

  bool initialize_committed_flux_history(MPI_Comm communicator) {
    FieldRegistry registry;
    FieldSchema schema;
    if (!registry.declare_field("piso.fixture_flux_history", 1U, 0U,
                                flux_history_dependency) ||
        !registry.freeze(schema)) {
      return false;
    }
    const std::array requests{ArenaFieldRequest{
        flux_history_dependency, {1, 1, 1}, {0U},
        FieldLifetime::state_layer}};
    ArenaLayout layout;
    if (!ArenaLayout::compile(schema,
                              {requests.data(), requests.size()}, layout) ||
        !StateLayers::allocate(layout, flux_history_layers) ||
        !AttemptTransaction::create(
            flux_history_layers.field_count(), 1U,
            flux_history_layers.field_count(), flux_history_transaction) ||
        !FaceFluxStorage::allocate_final(patch.cells, committed_flux_storage) ||
        !committed_flux_authority.claim(171U, 0U, flux_history_transaction,
                                        committed_flux_writer)) {
      return false;
    }
    constexpr std::array<double, 3U> zero_velocity{};
    return commit_uniform_flux_history(communicator, zero_velocity,
                                       zero_velocity, 1711U, 1721U);
  }

  bool commit_uniform_flux_history(
      MPI_Comm communicator,
      const std::array<double, 3U>& previous_velocity,
      const std::array<double, 3U>& accepted_velocity,
      RevisionToken previous_seed, RevisionToken accepted_seed) {
    if (!commit_uniform_flux(communicator, previous_seed,
                             previous_velocity) ||
        !commit_uniform_flux(communicator, accepted_seed,
                             accepted_velocity)) {
      return false;
    }
    ConstFaceFluxView accepted;
    ConstFaceFluxView previous;
    if (!committed_flux_writer.committed(committed_flux_storage, accepted) ||
        !committed_flux_writer.committed_previous(committed_flux_storage,
                                                  previous)) {
      return false;
    }
    committed_flux_accepted = accepted;
    committed_flux_previous = previous;
    return true;
  }

  bool commit_uniform_flux(MPI_Comm communicator, RevisionToken seed,
                           const std::array<double, 3U>& velocity) {
    if (!flux_history_transaction.begin(flux_history_layers) ||
        !flux_history_transaction.revise_trial(flux_history_dependency)) {
      return false;
    }
    const RevisionDependency dependency{
        AttemptTransaction::field_revision_source(flux_history_dependency),
        flux_history_transaction.trial_revision(flux_history_dependency)};
    PendingFaceFluxView pending;
    if (!committed_flux_writer.begin_pending(
            flux_history_transaction, committed_flux_storage, pending)) {
      return false;
    }
    OwnedField history_density =
        make_field(91U, patch.cells, 1U, 2U, seed, seed + 1000U);
    OwnedField history_velocity =
        make_field(92U, patch.cells, 3U, 2U, seed + 1U, seed + 1001U);
    fill(history_density, 1.0);
    const std::size_t component_size =
        history_velocity.view.component_stride;
    for (std::size_t component = 0U; component < velocity.size();
         ++component) {
      const auto begin = history_velocity.storage.begin() +
                         static_cast<std::ptrdiff_t>(component *
                                                     component_size);
      std::fill_n(begin, component_size, velocity[component]);
    }
    const std::array<ConstFieldView, 2U> reads{
        as_const(history_density.view), as_const(history_velocity.view)};
    const KernelInvocation invocation{
        {reads.data(), reads.size()}, {}, {{0, 0, 0}, patch.cells},
        0U, 0U, 1U, 0U, nullptr};
    const std::array dependencies{dependency};
    return static_cast<bool>(reconstruct_mass_flux(
               equations.kernels(), invocation, pending)) &&
           static_cast<bool>(committed_flux_writer.publish_pending(
               {dependencies.data(), dependencies.size()}, pending)) &&
           static_cast<bool>(flux_history_transaction.collective_finish(
               communicator, Status{}));
  }

  PisoIntermediateInput intermediate_input(BdfCoefficients bdf,
                                           RevisionToken time_revision,
                                           std::uint8_t corrector = 1U,
                                           RevisionToken prior = 0U) {
    PisoIntermediateInput input;
    input.momentum = {equations.momentum().fingerprint(),
                      EquationAssemblyScope::momentum_predictor,
                      time_revision,
                      geometry.topology_revision(),
                      1702U,
                      1703U};
    input.predictor.plan =
        equations.thermophysical_predictor().fingerprint();
    input.predictor.time = time_revision;
    input.predictor.geometry = geometry.topology_revision();
    input.predictor.accepted_face_flux = committed_flux_accepted.revision;
    input.predictor.previous_face_flux =
        bdf.order == 2U ? committed_flux_previous.revision : RevisionToken{};
    input.predictor.committed_face_flux_authority =
        committed_flux_accepted.certificate.authority();
    input.predictor.committed_face_flux_storage =
        committed_flux_accepted.certificate.storage();
    input.predictor.committed_face_flux_revision_domain =
        committed_flux_accepted.certificate.revision_domain();
    input.predictor.predicted_density = density.view.revision;
    input.predictor.predicted_density_storage =
        density.view.storage_identity;
    input.predictor.predicted_density_revision_domain =
        density.view.revision_domain;
    input.predictor.paired_face_flux = trial_flux.revision;
    input.predictor.paired_face_flux_storage =
        trial_flux.x.storage_identity;
    input.predictor.paired_face_flux_revision_domain =
        trial_flux.x.revision_domain;
    input.predictor.state = 1705U;
    input.predictor.order = bdf.order;
    input.pressure_reference = {
        equations.pressure_reference().fingerprint(),
        equations.thermophysical_predictor().fingerprint(),
        thermodynamics.fingerprint(),
        1706U,
        time_revision,
        1707U,
        PressureReferenceKind::closed_mass};
    input.density = density.view;
    input.trial_velocity = as_const(velocity.view);
    input.trial_flux = as_const(trial_flux);
    input.temporal_reference =
        corrector == 1U ? as_const(phi_h_by_a) : ConstFaceFluxView{};
    input.committed_face_history.accepted =
        corrector == 1U ? committed_flux_accepted : ConstFaceFluxView{};
    input.committed_face_history.previous =
        corrector == 1U && bdf.order == 2U
            ? committed_flux_previous
            : ConstFaceFluxView{};
    input.momentum_system = {momentum_diagonal.view, momentum_rhs.view};
    input.bdf = bdf;
    input.numeric_boundary = boundary.revision();
    input.prior_corrector = prior;
    input.corrector = corrector;
    return input;
  }

  Status prepare_closed_gauge(
      const PressureCorrectionCertificate& pressure,
      const PressureReferenceCertificate& predecessor,
      double pressure_reference, ConstFieldView pressure_perturbation,
      ConstFieldView raw_pressure_correction,
      ConstFieldView candidate_pressure_compressibility,
      PlanFingerprint exact_eos_closure, ReductionEngine& reductions,
      ClosedGaugeCorrectionCertificate& certificate) const noexcept {
    const PressureEnergyCellActivity activity{
        {}, {}, {}};
    const ClosedGaugeCorrectionPrepareInput input{
        predecessor,
        pressure_reference,
        pressure.corrector,
        pressure.time,
        pressure.geometry,
        pressure.state,
        exact_eos_closure,
        pressure_perturbation,
        raw_pressure_correction,
        candidate_pressure_compressibility,
        activity};
    return equations.pressure_reference().prepare_closed_gauge_correction(
        input, reductions, certificate);
  }

  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  ContributionRegistry contributions;
  EquationPlanSet equations;
  PisoPlan piso;
  HaloEngine halo;
  HaloEngine correction_halo;
  PressureVelocityCoupler coupler;
  OwnedField density;
  OwnedField velocity;
  OwnedField momentum_diagonal;
  OwnedField momentum_rhs;
  OwnedField r_au;
  OwnedField h_by_a;
  OwnedField pressure_gradient;
  OwnedFaceField x_pressure_coefficient;
  OwnedFaceField y_pressure_coefficient;
  OwnedFaceField z_pressure_coefficient;
  FaceFluxStorage phi_storage;
  FaceFluxView phi_h_by_a{};
  FaceFluxView trial_flux{};
  FieldId flux_history_dependency{};
  StateLayers flux_history_layers;
  AttemptTransaction flux_history_transaction;
  FaceFluxStorage committed_flux_storage;
  FinalFaceFluxAuthority committed_flux_authority;
  FinalFaceFluxWriter committed_flux_writer;
  ConstFaceFluxView committed_flux_accepted{};
  ConstFaceFluxView committed_flux_previous{};
};

}  // namespace hundun::v04::test
