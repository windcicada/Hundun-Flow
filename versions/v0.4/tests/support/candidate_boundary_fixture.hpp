// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "piso_fixture.hpp"

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace hundun::v04::test {

enum class CandidateBoundaryInlet : std::uint8_t {
  velocity,
  mass_flow,
};

struct CandidateBoundaryFixtureSpec {
  CandidateBoundaryInlet inlet{CandidateBoundaryInlet::velocity};
  double inlet_velocity{1.0};
  double mass_flow_rate{0.2};
  double inlet_temperature{300.0};
  double outlet_pressure{101325.0};
  bool allow_backflow{};
  bool immersed{};
  bool immersed_touches_x_min{};
  bool immersed_touches_x_max{};
  bool reverse_open_boundaries{};
  bool multispecies{};
  double backflow_velocity{-0.25};
  double backflow_temperature{350.0};
  std::int32_t cells_per_axis{6};
};

struct CandidateBoundaryScratch {
  OwnedField raw_pressure;
  OwnedField scaled_pressure;
  OwnedField raw_enthalpy;
  OwnedField scaled_enthalpy;
  OwnedField pressure;
  OwnedField enthalpy;
  OwnedField density;
  OwnedField temperature;
  OwnedField velocity;
  OwnedField pressure_compressibility;
  std::vector<OwnedField> independent_species;
  std::vector<ConstFieldView> independent_species_views;
  std::vector<ConstFieldView> semantic_independent_species_views;
  std::vector<ConstFieldView> thermophysical_species_aliases;
  std::array<OwnedField, 6U> material;
  HaloEngine correction_halo;
  HaloEngine state_halo;
  FaceFluxStorage mechanical_storage;
  FaceFluxView mechanical_flux{};
  FaceFluxStorage final_storage;
  FaceFluxView final_flux{};
  PisoFrozenMomentumPressureStageCertificate pressure_stage{};
  PisoFrozenMomentumVelocityStageCertificate velocity_stage{};
  PisoFrozenMomentumFluxStageCertificate flux_stage{};
  BoundaryThermophysicalGhostCertificate boundary_thermophysics{};
  ConstFieldView thermophysical_pressure_alias{};
  ConstFieldView thermophysical_enthalpy_alias{};
  FinalBoundaryFluxCertificate final_boundary{};
  PisoExactThermodynamicCandidateView thermodynamic{};
};

class CandidateBoundaryFixture {
 public:
  bool initialize(MPI_Comm communicator,
                  CandidateBoundaryFixtureSpec requested = {}) {
    diagnostic_step = 1U;
    communicator_ = communicator;
    spec = requested;
    if (communicator == MPI_COMM_NULL || spec.cells_per_axis < 2 ||
        (spec.immersed_touches_x_min && spec.immersed_touches_x_max) ||
        ((spec.immersed_touches_x_min || spec.immersed_touches_x_max) &&
         !spec.immersed) ||
        !std::isfinite(spec.inlet_temperature) ||
        !(spec.inlet_temperature > 0.0) ||
        !std::isfinite(spec.outlet_pressure) ||
        !(spec.outlet_pressure > 0.0)) {
      return false;
    }
    int size = 0;
    if (MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0)
      return false;

    CartesianMeshSpec mesh;
    mesh.kind = GeometryKind::uniform;
    mesh.lower = {0.0, 0.0, 0.0};
    mesh.upper = {1.0, 1.0, 1.0};
    mesh.has_exact_cells = true;
    const std::int32_t nx =
        std::max(spec.cells_per_axis, 2 * size);
    mesh.exact_cells = {nx, spec.cells_per_axis, spec.cells_per_axis};
    mesh.minimum_spacing = {
        1.0 / static_cast<double>(mesh.exact_cells.x),
        1.0 / static_cast<double>(mesh.exact_cells.y),
        1.0 / static_cast<double>(mesh.exact_cells.z)};
    mesh.max_growth_ratio = 1.0;
    mesh.limits.max_global_cells =
        static_cast<std::uint64_t>(mesh.exact_cells.x) *
        static_cast<std::uint64_t>(mesh.exact_cells.y) *
        static_cast<std::uint64_t>(mesh.exact_cells.z);
    mesh.limits.max_memory_bytes_per_rank = 1U << 29U;

    ValidatedModel model;
    model.mesh = mesh;
    model.fingerprint = UINT64_C(0x63616e64626e6479) +
                        static_cast<std::uint8_t>(spec.inlet) +
                        (spec.allow_backflow ? 17U : 0U) +
                        (spec.multispecies ? 31U : 0U) +
                        (spec.reverse_open_boundaries ? 47U : 0U);
    model.pressure_reference = PressureReferenceKind::boundary_absolute;
    for (BoundaryFaceSpec& face : model.boundaries) {
      face.flow_kind = BoundaryKind::symmetry;
      face.thermal_kind = BoundaryKind::none;
    }
    const std::size_t inlet_face = spec.reverse_open_boundaries ? 1U : 0U;
    const std::size_t outlet_face = spec.reverse_open_boundaries ? 0U : 1U;
    BoundaryFaceSpec& inlet = model.boundaries[inlet_face];
    inlet.flow_kind =
        spec.inlet == CandidateBoundaryInlet::velocity
            ? BoundaryKind::velocity_inlet
            : BoundaryKind::mass_flow_inlet;
    const double inlet_sign = spec.reverse_open_boundaries ? -1.0 : 1.0;
    inlet.velocity = {inlet_sign * spec.inlet_velocity, 0.0, 0.0};
    inlet.direction = {inlet_sign, 0.0, 0.0};
    inlet.mass_flow_rate = spec.mass_flow_rate;
    inlet.temperature = spec.inlet_temperature;
    BoundaryFaceSpec& outlet = model.boundaries[outlet_face];
    outlet.flow_kind = BoundaryKind::pressure_outlet;
    outlet.pressure = spec.outlet_pressure;
    outlet.allow_backflow = spec.allow_backflow;
    outlet.backflow_velocity = {
        spec.reverse_open_boundaries ? -spec.backflow_velocity
                                     : spec.backflow_velocity,
        0.0, 0.0};
    outlet.backflow_temperature = spec.backflow_temperature;
    if (spec.multispecies) {
      model.transported_scalars.push_back(
          {"sp0", TransportedScalarRole::species, 1.0, 1.0});
      model.transported_scalars.push_back(
          {"sp1", TransportedScalarRole::species, 1.0, 1.0});
      for (std::size_t face_index = 0U;
           face_index < model.boundaries.size(); ++face_index) {
        BoundaryFaceSpec& face = model.boundaries[face_index];
        const bool is_inlet = face_index == inlet_face;
        const bool is_outlet = face_index == outlet_face;
        face.scalars.push_back(
            {"sp0",
             is_inlet ? ScalarBoundaryKind::dirichlet
                      : ScalarBoundaryKind::zero_gradient,
             is_inlet ? 0.2 : 0.0,
             is_outlet ? ScalarBoundaryKind::dirichlet
                       : ScalarBoundaryKind::zero_gradient,
             is_outlet ? 0.25 : 0.0});
        face.scalars.push_back(
            {"sp1",
             is_inlet ? ScalarBoundaryKind::dirichlet
                      : ScalarBoundaryKind::zero_gradient,
             is_inlet ? 0.3 : 0.0,
             is_outlet ? ScalarBoundaryKind::dirichlet
                       : ScalarBoundaryKind::zero_gradient,
             is_outlet ? 0.35 : 0.0});
      }
    }
    model.schemes.momentum = ConvectionScheme::central2;
    model.schemes.enthalpy = ConvectionScheme::central2;
    model.schemes.species = ConvectionScheme::central2;
    model.schemes.passive_scalar = ConvectionScheme::central2;

    diagnostic_step = 2U;
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
                                   registry, boundary, schemes, time) ||
        !CartesianKernelPlan::compile(schemes, geometry, patch, boundary,
                                      kernels)) {
      return false;
    }
    if (spec.multispecies &&
        (!registry.require_field("sp0", 1U, 2U, id) || id != 8U ||
         !registry.require_field("sp1", 1U, 2U, id) || id != 9U))
      return false;
    if (spec.immersed) {
      const std::array<TriangleInput, 12U> triangles =
          immersed_cube(spec.immersed_touches_x_min,
                        spec.immersed_touches_x_max);
      ImmersedPlanLimits limits;
      ImmersedDomainBoundaryPolicy boundary_policy;
      boundary_policy.allow_one_sided_quadratic[0U] =
          spec.immersed_touches_x_min;
      boundary_policy.allow_one_sided_quadratic[1U] =
          spec.immersed_touches_x_max;
      if (!StlScanCompiler::compile_triangles(
              geometry, patch, {triangles.data(), triangles.size()},
              CartesianAxis::y, immersed_scan_budget, immersed_scan) ||
          !ImmersedSurfaceCompiler::compile(immersed_scan,
                                             immersed_surface) ||
          !EBTopologyCompiler::compile(
              communicator, geometry, patch, immersed_scan,
              immersed_surface, ImmersedFluidSide::outside, limits,
              immersed_topology) ||
          !BoundaryStencilCompiler::compile(
              communicator, geometry, patch, immersed_surface,
              immersed_topology, boundary_policy, limits,
              immersed_boundary) ||
          !make_continuity_activity()) {
        return false;
      }
    }
    diagnostic_step = 3U;
    ThermophysicalSpec thermophysics = test_thermophysical_spec();
    if (spec.multispecies) {
      thermophysics.species.clear();
      for (std::size_t species_index = 0U; species_index < 3U;
           ++species_index) {
        SpeciesThermophysicalSpec species = test_air();
        species.stable_name = "sp" + std::to_string(species_index);
        species.molecular_weight += static_cast<double>(species_index) * 2.0;
        thermophysics.species.push_back(std::move(species));
      }
    }
    const Span<const TransportedScalarSpec> transported_scalars =
        spec.multispecies
            ? Span<const TransportedScalarSpec>{
                  model.transported_scalars.data(),
                  model.transported_scalars.size()}
            : Span<const TransportedScalarSpec>{};
    if (!ThermodynamicsPlan::compile(thermophysics, transported_scalars,
                                     thermodynamics) ||
        !TransportPlan::compile(thermophysics, thermodynamics, transport))
      return false;
    diagnostic_step = 4U;
    std::vector<FieldId> declared{0U, 1U, 2U, 3U, 4U, 5U, 6U, 7U};
    if (spec.multispecies) {
      declared.push_back(8U);
      declared.push_back(9U);
    }
    if (!contributions.configure({declared.data(), declared.size()}) ||
        !contributions.freeze())
      return false;
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
    equation_spec.pressure_reference =
        PressureReferenceKind::boundary_absolute;
    if (spec.multispecies)
      equation_spec.scalars =
          {scalar_equations.data(), scalar_equations.size()};
    equation_spec.closed_mass_service_stage = 0U;
    equation_spec.maximum_cells_per_rank =
        static_cast<std::size_t>(mesh.limits.max_global_cells);
    PisoPlanSpec piso_spec = test_piso_spec();
    piso_spec.energy_tolerance = 1.0e-10;
    if (!EquationPlanSet::compile(communicator, schemes, geometry, patch,
                                  boundary, contributions, thermodynamics,
                                  transport, equation_spec, equations) ||
        !PisoPlan::compile(communicator, equations, piso_spec, piso))
      return false;
    if (spec.immersed) {
      if (!IbmEquationInterfacePlan::compile(
              equations.kernels(), immersed_topology, immersed_boundary,
              immersed_topology.interface_metric(), immersed_interface) ||
          !IbmPhysicalBoundaryFluxAuthority::compile(
              communicator, geometry, patch, immersed_topology,
              immersed_interface, immersed_physical_boundary_flux))
        return false;
    }

    diagnostic_step = 5U;
    const Int3 cells = patch.cells;
    const std::uint8_t reach = equations.kernels().reach();
    const std::uint8_t state_reach =
        spec.immersed
            ? std::max(reach, immersed_boundary.maximum_halo_reach())
            : reach;
    density = make_field(0U, cells, 1U, state_reach, 1801U, 2801U);
    velocity = make_field(1U, cells, 3U, 0U, 1802U, 2802U);
    pressure = make_field(2U, cells, 1U, state_reach, 1803U, 2803U);
    enthalpy = make_field(3U, cells, 1U, state_reach, 1804U, 2804U);
    temperature = make_field(4U, cells, 1U, state_reach, 1805U, 2805U);
    momentum_diagonal = make_field(30U, cells, 3U, 0U, 1806U, 2806U);
    momentum_rhs = make_field(31U, cells, 3U, 0U, 1807U, 2807U);
    r_au = make_field(40U, cells, 3U, 1U, 1808U, 2808U);
    h_by_a = make_field(41U, cells, 3U, reach, 1809U, 2809U);
    pressure_gradient = make_field(42U, cells, 3U, 1U, 1810U, 2810U);
    x_pressure_coefficient =
        make_face_field(CartesianAxis::x, cells, 2811U);
    y_pressure_coefficient =
        make_face_field(CartesianAxis::y, cells, 2812U);
    z_pressure_coefficient =
        make_face_field(CartesianAxis::z, cells, 2813U);
    accepted_density = make_field(50U, cells, 1U, 0U, 1811U, 2814U);
    previous_density = make_field(51U, cells, 1U, 0U, 1812U, 2815U);
    pressure_compressibility =
        make_field(6U, cells, 1U, 0U, 1813U, 2816U);
    pressure_diagonal = make_field(53U, cells, 1U, 0U, 1814U, 2817U);
    pressure_rhs = make_field(54U, cells, 1U, 0U, 1815U, 2818U);
    for (std::size_t index = 0U; index < base_material.size(); ++index)
      base_material[index] = make_field(
          static_cast<FieldId>(60U + index), cells, 1U, state_reach,
          static_cast<RevisionToken>(1820U + index),
          static_cast<StorageIdentity>(2820U + index));
    if (spec.multispecies) {
      independent_species.resize(2U);
      independent_species[0U] =
          make_field(8U, cells, 1U, state_reach, 1827U, 2827U);
      independent_species[1U] =
          make_field(9U, cells, 1U, state_reach, 1828U, 2828U);
    }
    if (!FaceFluxStorage::allocate_workspace(cells, 2U, phi_storage) ||
        !phi_storage.workspace_view(0U, 1830U, phi_h_by_a) ||
        !phi_storage.workspace_view(1U, 1831U, trial_flux))
      return false;

    diagnostic_step = 6U;
    double base_h = 0.0;
    double base_cp = 0.0;
    double base_gas = 0.0;
    const std::array<double, 2U> base_composition{{0.2, 0.3}};
    const Span<const double> base_composition_view =
        spec.multispecies
            ? Span<const double>{base_composition.data(),
                                 base_composition.size()}
            : Span<const double>{};
    if (!thermodynamics.mixture_enthalpy(
            spec.inlet_temperature, base_composition_view, base_h, base_cp,
            base_gas))
      return false;
    boundary_scalar_values.assign(boundary.resolved_scalar_count(), 0.0);
    boundary_vector_values.assign(boundary.resolved_vector_count(), {});
    boundary_normal_gradient_values.assign(
        boundary.resolved_normal_gradient_count(), 0.0);
    const Span<const BoundaryIndexSpan> boundary_spans = boundary.spans();
    for (std::size_t span_index = 0U; span_index < boundary_spans.size;
         ++span_index) {
      const BoundaryIndexSpan& span = boundary_spans.data[span_index];
      const std::size_t begin = span.resolved_begin;
      const std::size_t end = begin + span.resolved_stride;
      if (span.value_source == BoundaryValueSource::resolved_scalar) {
        double value = 0.0;
        if (span.stage == BoundaryStage::pressure) {
          diagnostic_status = boundary.pressure_perturbation_target(
              span.face, absolute_pressure_reference, value);
          if (!diagnostic_status) return false;
        } else if (span.stage == BoundaryStage::enthalpy) {
          value = base_h;
        } else if (span.stage == BoundaryStage::scalar &&
                   span.parameter < boundary.parameter_count() &&
                   spec.multispecies &&
                   (span.field == 8U || span.field == 9U)) {
          // Conditional outlet scalar spans consume resolved values.  The
          // initial same-layer thermophysical closure is the outflow branch,
          // so its resolved value is the owned base composition rather than
          // the compiled ordinary zero-gradient placeholder.  The candidate
          // finalizer independently selects configured backflow Y from the
          // final/provisional mass-flux sign.
          value = base_composition[span.field == 8U ? 0U : 1U];
        } else {
          return false;
        }
        if (begin > boundary_scalar_values.size() ||
            end > boundary_scalar_values.size())
          return false;
        std::fill(boundary_scalar_values.begin() + begin,
                  boundary_scalar_values.begin() + end, value);
      } else if (span.value_source ==
                 BoundaryValueSource::resolved_vector) {
        if (span.parameter >= boundary.parameter_count() ||
            begin > boundary_vector_values.size() ||
            end > boundary_vector_values.size())
          return false;
        const Real3 value{
            boundary.velocity_x().data[span.parameter],
            boundary.velocity_y().data[span.parameter],
            boundary.velocity_z().data[span.parameter]};
        std::fill(boundary_vector_values.begin() + begin,
                  boundary_vector_values.begin() + end, value);
      } else if (span.value_source ==
                 BoundaryValueSource::resolved_normal_gradient) {
        if (begin > boundary_normal_gradient_values.size() ||
            end > boundary_normal_gradient_values.size())
          return false;
        std::fill(boundary_normal_gradient_values.begin() + begin,
                  boundary_normal_gradient_values.begin() + end, 0.0);
      }
    }
    fill(velocity, 0.0);
    fill(pressure, 0.0);
    fill(enthalpy, base_h);
    fill(temperature, spec.inlet_temperature);
    fill(density, 1.0);
    for (std::size_t species_index = 0U;
         species_index < independent_species.size(); ++species_index)
      fill(independent_species[species_index],
           base_composition[species_index]);
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          ThermoState state;
          if (!thermodynamics.evaluate(
                  absolute_pressure_reference,
                  enthalpy.view.unchecked(cell, 0U), base_composition_view,
                  {}, state))
            return false;
          density.view.unchecked(cell, 0U) = state.rho;
          temperature.view.unchecked(cell, 0U) = state.temperature;
          pressure_compressibility.view.unchecked(cell, 0U) =
              state.drho_dp_hY;
        }
    if (!independent_species.empty()) {
      std::vector<FieldView> species_views;
      species_views.reserve(independent_species.size());
      for (OwnedField& species : independent_species)
        species_views.push_back(species.view);
      if (!apply_boundary_ghosts(
              BoundaryStage::scalar, boundary,
              {species_views.data(), species_views.size()},
              resolved_boundary_values()))
        return false;
      for (std::size_t species_index = 0U;
           species_index < independent_species.size(); ++species_index)
        independent_species[species_index].view = species_views[species_index];
    }
    fill(accepted_density, density.view.unchecked({0, 0, 0}, 0U));
    fill(previous_density, density.view.unchecked({0, 0, 0}, 0U));
    fill(momentum_rhs, 0.0);
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const double volume = cell_volume(cell);
          for (std::uint8_t component = 0U; component < 3U; ++component)
            momentum_diagonal.view.unchecked(cell, component) = volume;
        }
    fill_face_flux(phi_h_by_a, 0.0);
    fill_face_flux(trial_flux, 0.0);

    diagnostic_step = 7U;
    const std::array<HaloFieldSpec, 4U> halo_fields{{
        {density.view.field, 1U, 1U},
        {r_au.view.field, 1U, 3U},
        {h_by_a.view.field, reach, 3U},
        {pressure_gradient.view.field, 1U, 3U}}};
    const std::array<HaloFieldSpec, 1U> correction_fields{{
        {correction_field, 1U, 1U}}};
    if (!halo.reserve(communicator, patch,
                      {halo_fields.data(), halo_fields.size()},
                      boundary.halo_topology()) ||
        !correction_halo.reserve(
            communicator, patch,
            {correction_fields.data(), correction_fields.size()},
            boundary.halo_topology()))
      return false;
    if (spec.immersed) {
      const std::array<RemoteDonorFieldSpec, 1U> live_fields{{
          {correction_field, 1U}}};
      const std::array<RemoteDonorFieldSpec, 1U> candidate_fields{{
          {candidate_pressure_correction_field, 1U}}};
      if (!RemoteDonorExchangePlan::compile(
              communicator, geometry.global_cells(), patch,
              immersed_boundary.reconstruction(),
              {live_fields.data(), live_fields.size()},
              live_pressure_donor_stage, live_pressure_donors) ||
          !RemoteDonorExchangePlan::compile(
              communicator, geometry.global_cells(), patch,
              immersed_boundary.reconstruction(),
              {candidate_fields.data(), candidate_fields.size()},
              candidate_pressure_donor_stage,
              candidate_pressure_donors))
        return false;
      candidate_pressure_donor_reach = candidate_pressure_donors.reach();
      if (candidate_pressure_donor_reach == 0U ||
          live_pressure_donors.reach() != candidate_pressure_donor_reach)
        return false;
    }
    diagnostic_step = 8U;
    const PisoCouplerWorkspace workspace{
        r_au.view, h_by_a.view, pressure_gradient.view,
        x_pressure_coefficient.view, y_pressure_coefficient.view,
        z_pressure_coefficient.view, phi_h_by_a};
    PisoCouplerServices services{
        communicator, &geometry, patch, &boundary, &thermodynamics, &halo,
        1840U,
        density.view.field, &correction_halo, 1841U, correction_field};
    if (spec.immersed) {
      services.continuity_activity = continuity_activity();
      services.pressure_correction_donors = &live_pressure_donors;
      services.pressure_correction_donor_stage =
          live_pressure_donor_stage;
      services.candidate_pressure_correction_donors =
          &candidate_pressure_donors;
      services.candidate_pressure_correction_donor_stage =
          candidate_pressure_donor_stage;
      services.candidate_pressure_correction_field =
          candidate_pressure_correction_field;
      services.candidate_pressure_correction_donor_reach =
          candidate_pressure_donor_reach;
      services.candidate_pressure_correction_donor_fingerprint =
          candidate_pressure_donors.fingerprint();
      services.immersed_interface = &immersed_interface;
    }
    if (!PressureVelocityCoupler::bind(
            piso, equations, services, workspace, coupler) ||
        !initialize_flux_history())
      return false;

    diagnostic_step = 9U;
    std::vector<ConstFieldView> base_species_views;
    base_species_views.reserve(independent_species.size());
    for (const OwnedField& species : independent_species)
      base_species_views.push_back(as_const(species.view));
    std::vector<BoundaryGhostFieldAuthority> base_authorities;
    base_authorities.reserve(2U + base_species_views.size());
    base_authorities.push_back(
        make_boundary_ghost_field_authority(as_const(pressure.view)));
    base_authorities.push_back(
        make_boundary_ghost_field_authority(as_const(enthalpy.view)));
    for (ConstFieldView species : base_species_views)
      base_authorities.push_back(
          make_boundary_ghost_field_authority(species));
    const BoundaryThermophysicalGhostAuthority base_authority{
        1842U, boundary.revision(), boundary.local_layout_fingerprint(),
        cells, {state_reach, state_reach, state_reach}, state_reach,
        {base_authorities.data(), base_authorities.size()}};
    const BoundaryThermophysicalGhostContext base_context{
        target_time, geometry.fingerprint(), pressure_reference_authority,
        boundary.revision(),
        BoundaryThermophysicalGhostPhase::corrector_one};
    diagnostic_status = BoundaryThermophysicalFaceClosure::close(
        boundary, thermodynamics, transport,
        {absolute_pressure_reference, as_const(pressure.view),
         as_const(enthalpy.view),
         {base_species_views.data(), base_species_views.size()},
         base_authority},
        {density.view, temperature.view, base_material[0U].view,
         base_material[1U].view, base_material[2U].view,
         base_material[3U].view, base_material[4U].view,
         base_material[5U].view},
        base_context, base_boundary_thermophysics);
    if (!diagnostic_status)
      return false;

    diagnostic_step = 10U;
    PisoIntermediateInput input;
    input.momentum = {equations.momentum().fingerprint(),
                      EquationAssemblyScope::momentum_predictor,
                      target_time, geometry.topology_revision(),
                      trial_flux.revision, 1843U,
                      fixture_time_step_for_bdf(bdf)};
    input.predictor.plan =
        equations.thermophysical_predictor().fingerprint();
    input.predictor.time = target_time;
    input.predictor.geometry = geometry.topology_revision();
    input.predictor.accepted_face_flux = committed_accepted.revision;
    input.predictor.committed_face_flux_authority =
        committed_accepted.certificate.authority();
    input.predictor.committed_face_flux_storage =
        committed_accepted.certificate.storage();
    input.predictor.committed_face_flux_revision_domain =
        committed_accepted.certificate.revision_domain();
    // The predictor density precedes the physical-face EOS refresh.  The
    // refresh input below carries the same storage lineage at its newer
    // thermophysical boundary revision.
    input.predictor.predicted_density = density.view.revision + 1U;
    input.predictor.predicted_density_storage = density.view.storage_identity;
    input.predictor.predicted_density_revision_domain =
        density.view.revision_domain;
    input.predictor.paired_face_flux = trial_flux.revision;
    input.predictor.paired_face_flux_storage =
        trial_flux.x.storage_identity;
    input.predictor.paired_face_flux_revision_domain =
        trial_flux.x.revision_domain;
    input.predictor.state = 1844U;
    input.predictor.order = bdf.order;
    input.pressure_reference = {
        equations.pressure_reference().fingerprint(),
        equations.thermophysical_predictor().fingerprint(),
        thermodynamics.fingerprint(), pressure_reference_revision,
        target_time, pressure_reference_authority,
        PressureReferenceKind::boundary_absolute};
    pressure_reference = input.pressure_reference;
    input.density = density.view;
    input.trial_velocity = as_const(velocity.view);
    input.trial_flux = as_const(trial_flux);
    input.temporal_reference = as_const(phi_h_by_a);
    input.committed_face_history.accepted = committed_accepted;
    input.momentum_system = {momentum_diagonal.view, momentum_rhs.view};
    input.bdf = bdf;
    input.numeric_boundary = boundary.revision();
    input.boundary_values = resolved_boundary_values();
    input.corrector = 1U;
    input.thermophysical_boundary = {
        base_boundary_thermophysics,
        {absolute_pressure_reference, as_const(pressure.view),
         as_const(enthalpy.view),
         {base_species_views.data(), base_species_views.size()},
         as_const(density.view)}};
    input.immersed_interface = spec.immersed ? &immersed_interface : nullptr;
    diagnostic_status = coupler.refresh(input, intermediate);
    if (!diagnostic_status)
      return false;

    diagnostic_step = 11U;
    const PressureCorrectionInput pressure_input{
        intermediate, pressure_reference, as_const(density.view),
        as_const(accepted_density.view), {},
        as_const(pressure_compressibility.view), bdf, target_time,
        geometry.topology_revision(), boundary.revision()};
    diagnostic_status = coupler.assemble_pressure_system(
        pressure_input, {pressure_diagonal.view, pressure_rhs.view},
        pressure_certificate);
    if (!diagnostic_status)
      return false;
    diagnostic_status = coupler.make_frozen_momentum_stage_authority(
        intermediate, pressure_certificate, authority);
    if (!diagnostic_status)
      return false;

    diagnostic_step = 12U;
    const PressureEnergyCandidateBoundaryFinalizerBinding finalizer_binding =
        this->finalizer_binding();
    if (!PressureEnergyCandidateBoundaryFinalizer::bind(
            finalizer_binding, finalizer))
      return false;
    diagnostic_step = 13U;
    const bool compiled_reductions = static_cast<bool>(ReductionEngine::compile(
        communicator, ReductionMode::mpi_allreduce, 8U, reductions));
    if (compiled_reductions) diagnostic_step = 0U;
    return compiled_reductions;
  }

  bool stage(double alpha, double raw_pressure_value,
             double raw_enthalpy_value, RevisionToken seed,
             CandidateBoundaryScratch& scratch) {
    diagnostic_step = 100U;
    diagnostic_status = {};
    if (!authority.valid() || !finalizer.ready() || seed == 0U)
      return false;
    const Int3 cells = patch.cells;
    const std::uint8_t reach = equations.kernels().reach();
    scratch.raw_pressure = make_field(
        correction_field, cells, 1U, 1U, seed + 1U, seed + 101U);
    scratch.scaled_pressure = make_field(
        candidate_pressure_correction_field, cells, 1U,
        spec.immersed ? candidate_pressure_donor_reach : 1U,
        seed + 2U, seed + 102U);
    scratch.raw_enthalpy =
        make_field(93U, cells, 1U, 0U, seed + 3U, seed + 103U);
    scratch.scaled_enthalpy =
        make_field(207U, cells, 1U, 0U, seed + 4U, seed + 104U);
    scratch.pressure =
        make_field(208U, cells, 1U, reach, seed + 5U, seed + 105U);
    scratch.enthalpy =
        make_field(203U, cells, 1U, reach, seed + 6U, seed + 106U);
    scratch.density =
        make_field(202U, cells, 1U, reach, seed + 7U, seed + 107U);
    scratch.temperature =
        make_field(204U, cells, 1U, reach, seed + 8U, seed + 108U);
    scratch.velocity =
        make_field(201U, cells, 3U, 1U, seed + 9U, seed + 109U);
    scratch.pressure_compressibility =
        make_field(206U, cells, 1U, 0U, seed + 10U, seed + 110U);
    for (std::size_t index = 0U; index < scratch.material.size(); ++index)
      scratch.material[index] = make_field(
          static_cast<FieldId>(220U + index), cells, 1U, reach,
          static_cast<RevisionToken>(seed + 20U + index),
          static_cast<StorageIdentity>(seed + 120U + index));
    if (spec.multispecies) {
      scratch.independent_species.resize(2U);
      scratch.thermophysical_species_aliases.resize(2U);
      scratch.independent_species[0U] =
          make_field(230U, cells, 1U, reach, seed + 26U, seed + 126U);
      scratch.independent_species[1U] =
          make_field(231U, cells, 1U, reach, seed + 27U, seed + 127U);
    }
    fill(scratch.raw_pressure, raw_pressure_value);
    fill(scratch.raw_enthalpy, raw_enthalpy_value);
    fill(scratch.scaled_pressure, -777.0);
    fill(scratch.scaled_enthalpy, 0.0);
    fill(scratch.pressure, pressure.view.unchecked({0, 0, 0}, 0U));
    fill(scratch.enthalpy, enthalpy.view.unchecked({0, 0, 0}, 0U));
    fill(scratch.density, 1.0);
    fill(scratch.temperature, spec.inlet_temperature);
    fill(scratch.velocity, 0.0);
    for (std::size_t species_index = 0U;
         species_index < scratch.independent_species.size(); ++species_index)
      fill(scratch.independent_species[species_index],
           species_index == 0U ? 0.2 : 0.3);

    diagnostic_step = 101U;
    const std::array<HaloFieldSpec, 1U> correction_contract{{
        {scratch.scaled_pressure.view.field, 1U, 1U}}};
    diagnostic_status = scratch.correction_halo.reserve(
        communicator_, patch,
        {correction_contract.data(), correction_contract.size()},
        boundary.halo_topology());
    if (!diagnostic_status)
      return false;
    diagnostic_status = coupler.form_frozen_momentum_scaled_pressure(
        authority, as_const(scratch.raw_pressure.view),
        scratch.correction_halo, alpha, scratch.scaled_pressure.view,
        scratch.pressure_stage);
    if (!diagnostic_status) return false;
    diagnostic_step = 102U;
    std::array<FieldView, 1U> correction_views{
        scratch.scaled_pressure.view};
    HaloTicket correction_ticket;
    Status status = scratch.correction_halo.begin(
        seed + 30U, {correction_views.data(), correction_views.size()},
        correction_ticket);
    if (status)
      status = scratch.correction_halo.finish(
          correction_ticket,
          {correction_views.data(), correction_views.size()});
    scratch.scaled_pressure.view = correction_views[0U];
    diagnostic_status = status;
    if (!diagnostic_status)
      return false;
    diagnostic_step = 103U;
    diagnostic_status = coupler.stage_frozen_momentum_velocity(
        authority, scratch.pressure_stage, scratch.correction_halo,
        scratch.scaled_pressure.view, scratch.velocity.view,
        scratch.velocity_stage);
    if (!diagnostic_status)
      return false;

    diagnostic_step = 104U;
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const double dp =
              scratch.scaled_pressure.view.unchecked(cell, 0U);
          const double dh = alpha == 0.0 ? 0.0 :
              alpha * scratch.raw_enthalpy.view.unchecked(cell, 0U);
          scratch.scaled_enthalpy.view.unchecked(cell, 0U) = dh;
          scratch.pressure.view.unchecked(cell, 0U) =
              pressure.view.unchecked(cell, 0U) + dp;
          scratch.enthalpy.view.unchecked(cell, 0U) =
              enthalpy.view.unchecked(cell, 0U) + dh;
          const Real3 candidate_velocity{
              scratch.velocity.view.unchecked(cell, 0U),
              scratch.velocity.view.unchecked(cell, 1U),
              scratch.velocity.view.unchecked(cell, 2U)};
          ThermoState state;
          const std::array<double, 2U> composition{{
              scratch.independent_species.empty()
                  ? 0.0
                  : scratch.independent_species[0U].view.unchecked(cell, 0U),
              scratch.independent_species.empty()
                  ? 0.0
                  : scratch.independent_species[1U].view.unchecked(cell,
                                                                    0U)}};
          const Span<const double> composition_view =
              scratch.independent_species.empty()
                  ? Span<const double>{}
                  : Span<const double>{composition.data(),
                                       composition.size()};
          diagnostic_status = thermodynamics.evaluate(
              absolute_pressure_reference +
                  scratch.pressure.view.unchecked(cell, 0U),
              scratch.enthalpy.view.unchecked(cell, 0U), composition_view,
              candidate_velocity, state);
          if (!diagnostic_status)
            return false;
          scratch.density.view.unchecked(cell, 0U) = state.rho;
          scratch.temperature.view.unchecked(cell, 0U) = state.temperature;
          scratch.pressure_compressibility.view.unchecked(cell, 0U) =
              state.drho_dp_hY;
        }

    FieldView semantic_pressure = scratch.pressure.view;
    semantic_pressure.field = boundary.pressure_field();
    FieldView semantic_enthalpy = scratch.enthalpy.view;
    semantic_enthalpy.field = boundary.enthalpy_field();
    if (spec.multispecies && spec.allow_backflow) {
      const bool configured_backflow =
          alpha > 0.0 && raw_pressure_value < 0.0;
      const std::array<double, 2U> outlet_composition{{
          configured_backflow ? 0.25 : 0.2,
          configured_backflow ? 0.35 : 0.3}};
      double outlet_enthalpy = enthalpy.view.unchecked({0, 0, 0}, 0U);
      if (configured_backflow) {
        double outlet_cp = 0.0;
        double outlet_gas = 0.0;
        diagnostic_status = thermodynamics.mixture_enthalpy(
            spec.backflow_temperature,
            {outlet_composition.data(), outlet_composition.size()},
            outlet_enthalpy, outlet_cp, outlet_gas);
        if (!diagnostic_status) return false;
      }
      for (std::size_t span_index = 0U;
           span_index < boundary.spans().size; ++span_index) {
        const BoundaryIndexSpan& span = boundary.spans().data[span_index];
        if (span.face != CartesianFace::x_max ||
            span.value_source != BoundaryValueSource::resolved_scalar)
          continue;
        const std::size_t begin = span.resolved_begin;
        const std::size_t end = begin + span.resolved_stride;
        if (begin > boundary_scalar_values.size() ||
            end > boundary_scalar_values.size())
          return false;
        if (span.stage == BoundaryStage::enthalpy) {
          std::fill(boundary_scalar_values.begin() + begin,
                    boundary_scalar_values.begin() + end,
                    outlet_enthalpy);
        } else if (span.stage == BoundaryStage::scalar &&
                   (span.field == 8U || span.field == 9U)) {
          std::fill(boundary_scalar_values.begin() + begin,
                    boundary_scalar_values.begin() + end,
                    outlet_composition[span.field == 8U ? 0U : 1U]);
        }
      }
    }
    diagnostic_step = 1041U;
    diagnostic_status = apply_boundary_ghosts(
        BoundaryStage::pressure, boundary, {&semantic_pressure, 1U},
        resolved_boundary_values());
    if (diagnostic_status) {
      diagnostic_step = 1042U;
      diagnostic_status = apply_boundary_ghosts(
          BoundaryStage::enthalpy, boundary, {&semantic_enthalpy, 1U},
          resolved_boundary_values());
    }
    std::vector<FieldView> semantic_species;
    if (diagnostic_status && !scratch.independent_species.empty()) {
      semantic_species.reserve(scratch.independent_species.size());
      for (std::size_t species_index = 0U;
           species_index < scratch.independent_species.size();
           ++species_index) {
        FieldView semantic = scratch.independent_species[species_index].view;
        semantic.field = independent_species[species_index].view.field;
        semantic_species.push_back(semantic);
      }
      diagnostic_status = apply_boundary_ghosts(
          BoundaryStage::scalar, boundary,
          {semantic_species.data(), semantic_species.size()},
          resolved_boundary_values());
    }
    if (!diagnostic_status) return false;
    scratch.thermophysical_pressure_alias = as_const(semantic_pressure);
    scratch.thermophysical_enthalpy_alias = as_const(semantic_enthalpy);
    for (std::size_t species_index = 0U;
         species_index < semantic_species.size(); ++species_index)
      scratch.thermophysical_species_aliases[species_index] =
          as_const(semantic_species[species_index]);
    std::vector<BoundaryGhostFieldAuthority> authorities;
    authorities.reserve(2U +
                        scratch.thermophysical_species_aliases.size());
    authorities.push_back(make_boundary_ghost_field_authority(
        scratch.thermophysical_pressure_alias));
    authorities.push_back(make_boundary_ghost_field_authority(
        scratch.thermophysical_enthalpy_alias));
    for (ConstFieldView species : scratch.thermophysical_species_aliases)
      authorities.push_back(make_boundary_ghost_field_authority(species));
    const BoundaryThermophysicalGhostAuthority boundary_authority{
        static_cast<std::uintptr_t>(seed + 31U), boundary.revision(),
        boundary.local_layout_fingerprint(), cells, {reach, reach, reach},
        reach, {authorities.data(), authorities.size()}};
    const BoundaryThermophysicalGhostContext boundary_context{
        target_time, geometry.fingerprint(), pressure_reference_authority,
        boundary.revision(),
        BoundaryThermophysicalGhostPhase::corrector_one};
    diagnostic_step = 105U;
    diagnostic_status = BoundaryThermophysicalFaceClosure::close(
        boundary, thermodynamics, transport,
        {absolute_pressure_reference, scratch.thermophysical_pressure_alias,
         scratch.thermophysical_enthalpy_alias,
         {scratch.thermophysical_species_aliases.data(),
          scratch.thermophysical_species_aliases.size()},
         boundary_authority},
        {scratch.density.view, scratch.temperature.view,
         scratch.material[0U].view, scratch.material[1U].view,
         scratch.material[2U].view, scratch.material[3U].view,
         scratch.material[4U].view, scratch.material[5U].view},
        boundary_context, scratch.boundary_thermophysics);
    if (!diagnostic_status)
      return false;

    diagnostic_step = 106U;
    std::vector<HaloFieldSpec> state_contract{
        {scratch.pressure.view.field, 1U, 1U},
        {scratch.enthalpy.view.field, 1U, 1U},
        {scratch.density.view.field, 1U, 1U},
        {scratch.temperature.view.field, 1U, 1U},
        {scratch.velocity.view.field, 1U, 3U}};
    for (const OwnedField& species : scratch.independent_species)
      state_contract.push_back({species.view.field, 1U, 1U});
    diagnostic_status = scratch.state_halo.reserve(
        communicator_, patch,
        {state_contract.data(), state_contract.size()},
        boundary.halo_topology());
    if (!diagnostic_status)
      return false;
    std::vector<FieldView> state_views{
        scratch.pressure.view, scratch.enthalpy.view,
        scratch.density.view, scratch.temperature.view,
        scratch.velocity.view};
    for (OwnedField& species : scratch.independent_species)
      state_views.push_back(species.view);
    HaloTicket state_ticket;
    status = scratch.state_halo.begin(
        seed + 32U, {state_views.data(), state_views.size()}, state_ticket);
    if (status)
      status = scratch.state_halo.finish(
          state_ticket, {state_views.data(), state_views.size()});
    scratch.pressure.view = state_views[0U];
    scratch.enthalpy.view = state_views[1U];
    scratch.density.view = state_views[2U];
    scratch.temperature.view = state_views[3U];
    scratch.velocity.view = state_views[4U];
    for (std::size_t species_index = 0U;
         species_index < scratch.independent_species.size(); ++species_index)
      scratch.independent_species[species_index].view =
          state_views[species_index + 5U];
    scratch.independent_species_views.clear();
    scratch.independent_species_views.reserve(
        scratch.independent_species.size());
    for (const OwnedField& species : scratch.independent_species)
      scratch.independent_species_views.push_back(as_const(species.view));
    scratch.semantic_independent_species_views.clear();
    scratch.semantic_independent_species_views.reserve(
        independent_species.size());
    for (const OwnedField& species : independent_species)
      scratch.semantic_independent_species_views.push_back(
          as_const(species.view));
    const PlanFingerprint composition_identity =
        exact_composition_identity_for_test(
            thermodynamics.fingerprint(),
            {scratch.independent_species_views.data(),
             scratch.independent_species_views.size()},
            cells);
    diagnostic_status = status;
    if (!diagnostic_status) return false;
    diagnostic_step = 107U;
    diagnostic_status = FaceFluxStorage::allocate_workspace(
        cells, 1U, scratch.mechanical_storage);
    if (diagnostic_status)
      diagnostic_status = scratch.mechanical_storage.workspace_view(
          0U, seed + 33U, scratch.mechanical_flux);
    if (diagnostic_status)
      diagnostic_status = coupler.stage_frozen_momentum_flux(
          authority, scratch.velocity_stage, as_const(scratch.density.view),
          scratch.mechanical_flux, scratch.flux_stage);
    if (diagnostic_status)
      diagnostic_status = FaceFluxStorage::allocate_workspace(
          cells, 1U, scratch.final_storage);
    if (diagnostic_status)
      diagnostic_status = scratch.final_storage.workspace_view(
          0U, seed + 34U, scratch.final_flux);
    if (!diagnostic_status)
      return false;
    fill_face_flux(scratch.final_flux, -991.0);
    const PressureEnergyCandidateBoundaryFinalizeInput finalizer_input{
        authority,
        scratch.pressure_stage,
        scratch.velocity_stage,
        scratch.flux_stage,
        pressure_reference,
        absolute_pressure_reference,
        as_const(scratch.pressure.view),
        as_const(scratch.enthalpy.view),
        as_const(scratch.density.view),
        as_const(scratch.temperature.view),
        as_const(scratch.velocity.view),
        {scratch.independent_species_views.data(),
         scratch.independent_species_views.size()},
        composition_identity,
        {scratch.boundary_thermophysics,
         {absolute_pressure_reference,
          scratch.thermophysical_pressure_alias,
          scratch.thermophysical_enthalpy_alias,
          {scratch.thermophysical_species_aliases.data(),
           scratch.thermophysical_species_aliases.size()},
          as_const(scratch.density.view)}},
        &scratch.state_halo,
        as_const(scratch.mechanical_flux),
        scratch.final_flux};
    diagnostic_step = 108U;
    diagnostic_status = finalizer.finalize(
        finalizer_input, reductions, scratch.final_boundary);
    if (!diagnostic_status) return false;

    PisoExactEosClosureIdentity closure;
    closure.thermodynamics = thermodynamics.fingerprint();
    closure.pressure_reference = pressure_reference_authority;
    closure.composition = composition_identity;
    closure.pressure_state =
        make_piso_field_revision_identity(as_const(pressure.view));
    closure.pressure_correction = make_piso_field_revision_identity(
        as_const(scratch.scaled_pressure.view));
    closure.enthalpy_state =
        make_piso_field_revision_identity(as_const(enthalpy.view));
    closure.enthalpy_correction = make_piso_field_revision_identity(
        as_const(scratch.scaled_enthalpy.view));
    closure.candidate_enthalpy = make_piso_field_revision_identity(
        as_const(scratch.enthalpy.view));
    closure.candidate_density = make_piso_field_revision_identity(
        as_const(scratch.density.view));
    closure.candidate_temperature = make_piso_field_revision_identity(
        as_const(scratch.temperature.view));
    closure.closure = seed + 35U;
    scratch.thermodynamic.enthalpy = as_const(scratch.enthalpy.view);
    scratch.thermodynamic.density = as_const(scratch.density.view);
    scratch.thermodynamic.temperature = as_const(scratch.temperature.view);
    scratch.thermodynamic.closure = closure;
    // Boundary-absolute candidates carry no closed-gauge derivative
    // authority.  Their EOS derivative remains a separate material scratch;
    // placing it in the exact view would make the open/closed contract
    // ambiguous.
    scratch.thermodynamic.pressure_compressibility = {};
    scratch.thermodynamic.closed_gauge = {};
    scratch.thermodynamic.independent_species = {
        scratch.semantic_independent_species_views.data(),
        scratch.semantic_independent_species_views.size()};
    diagnostic_step = 109U;
    if (!scratch.final_boundary.valid()) return false;
    diagnostic_step = 0U;
    return true;
  }

  PisoFrozenMomentumExactCandidateInput exact_input(
      const CandidateBoundaryScratch& scratch) const noexcept {
    return {as_const(scratch.raw_enthalpy.view),
            as_const(scratch.scaled_pressure.view),
            as_const(scratch.scaled_enthalpy.view),
            {velocity.view, pressure.view, enthalpy.view, density.view,
             temperature.view},
            as_const(scratch.pressure.view),
            scratch.thermodynamic,
            as_const(scratch.velocity.view),
            as_const(scratch.final_flux),
            scratch.final_boundary,
            {scratch.independent_species_views.data(),
             scratch.independent_species_views.size()}};
  }

  PressureEnergyCandidateBoundaryFinalizeInput finalizer_input(
      CandidateBoundaryScratch& scratch) const noexcept {
    const PlanFingerprint composition_identity =
        exact_composition_identity_for_test(
            thermodynamics.fingerprint(),
            {scratch.independent_species_views.data(),
             scratch.independent_species_views.size()},
            patch.cells);
    return {authority,
            scratch.pressure_stage,
            scratch.velocity_stage,
            scratch.flux_stage,
            pressure_reference,
            absolute_pressure_reference,
            as_const(scratch.pressure.view),
            as_const(scratch.enthalpy.view),
            as_const(scratch.density.view),
            as_const(scratch.temperature.view),
            as_const(scratch.velocity.view),
            {scratch.independent_species_views.data(),
             scratch.independent_species_views.size()},
            composition_identity,
            {scratch.boundary_thermophysics,
             {absolute_pressure_reference,
              scratch.thermophysical_pressure_alias,
              scratch.thermophysical_enthalpy_alias,
              {scratch.thermophysical_species_aliases.data(),
               scratch.thermophysical_species_aliases.size()},
              as_const(scratch.density.view)}},
            &scratch.state_halo,
            as_const(scratch.mechanical_flux),
            scratch.final_flux};
  }

  PressureEnergyCandidateBoundaryFinalizerBinding finalizer_binding()
      noexcept {
    PressureEnergyCandidateBoundaryFinalizerBinding binding{
        communicator_, &geometry, patch, &boundary, &equations.kernels(),
        &thermodynamics, &transport, &coupler};
    if (spec.immersed) {
      binding.immersed_interface = &immersed_interface;
      binding.immersed_physical_boundary_flux =
          &immersed_physical_boundary_flux;
      binding.candidate_pressure_correction_donors =
          &candidate_pressure_donors;
      binding.candidate_pressure_correction_donor_stage =
          candidate_pressure_donor_stage;
      binding.candidate_pressure_correction_field =
          candidate_pressure_correction_field;
      binding.candidate_pressure_correction_donor_reach =
          candidate_pressure_donor_reach;
    }
    return binding;
  }

  std::size_t immersed_link_count() const noexcept {
    return spec.immersed ? immersed_topology.links().size : 0U;
  }

  Status validate_zero_interface_flux(ConstFaceFluxView flux) const noexcept {
    return spec.immersed
               ? immersed_interface.validate_interface_flux(flux, 0.0)
               : Status{StatusCode::invalid_plan, 0U};
  }

  double cell_volume(Int3 local) const noexcept {
    const Int3 global{patch.begin.x + local.x, patch.begin.y + local.y,
                      patch.begin.z + local.z};
    return geometry.x().widths().data[global.x] *
           geometry.y().widths().data[global.y] *
           geometry.z().widths().data[global.z];
  }

  double face_area(CartesianAxis axis, Int3 local_face) const noexcept {
    const Int3 global{patch.begin.x + local_face.x,
                      patch.begin.y + local_face.y,
                      patch.begin.z + local_face.z};
    if (axis == CartesianAxis::x)
      return geometry.y().widths().data[global.y] *
             geometry.z().widths().data[global.z];
    if (axis == CartesianAxis::y)
      return geometry.x().widths().data[global.x] *
             geometry.z().widths().data[global.z];
    return geometry.x().widths().data[global.x] *
           geometry.y().widths().data[global.y];
  }

  bool local_face_owner(CartesianFace face) const noexcept {
    const BoundaryFacePlan* plan = nullptr;
    return boundary.face(face, plan) && plan != nullptr && plan->local_owner;
  }

  CandidateBoundaryFixtureSpec spec{};
  std::uint32_t diagnostic_step{};
  Status diagnostic_status{};
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  CartesianKernelPlan kernels;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  ContributionRegistry contributions;
  EquationPlanSet equations;
  PisoPlan piso;
  HaloEngine halo;
  HaloEngine correction_halo;
  PressureVelocityCoupler coupler;
  PressureEnergyCandidateBoundaryFinalizer finalizer;
  ReductionEngine reductions;
  StlScanPlan immersed_scan;
  ImmersedSurfacePlan immersed_surface;
  EBTopology immersed_topology;
  BoundaryStencilPlan immersed_boundary;
  IbmEquationInterfacePlan immersed_interface;
  IbmPhysicalBoundaryFluxAuthority immersed_physical_boundary_flux;
  RemoteDonorExchangePlan live_pressure_donors;
  RemoteDonorExchangePlan candidate_pressure_donors;
  OwnedField density;
  OwnedField velocity;
  OwnedField pressure;
  OwnedField enthalpy;
  OwnedField temperature;
  OwnedField momentum_diagonal;
  OwnedField momentum_rhs;
  OwnedField r_au;
  OwnedField h_by_a;
  OwnedField pressure_gradient;
  OwnedFaceField x_pressure_coefficient;
  OwnedFaceField y_pressure_coefficient;
  OwnedFaceField z_pressure_coefficient;
  OwnedField accepted_density;
  OwnedField previous_density;
  OwnedField pressure_compressibility;
  OwnedField pressure_diagonal;
  OwnedField pressure_rhs;
  std::array<OwnedField, 6U> base_material;
  std::vector<OwnedField> independent_species;
  std::vector<double> boundary_scalar_values;
  std::vector<Real3> boundary_vector_values;
  std::vector<double> boundary_normal_gradient_values;
  std::vector<std::uint8_t> continuity_cell_activity;
  std::vector<std::uint8_t> continuity_x_activity;
  std::vector<std::uint8_t> continuity_y_activity;
  std::vector<std::uint8_t> continuity_z_activity;
  PlanFingerprint continuity_local_fingerprint{};
  PlanFingerprint continuity_collective_fingerprint{};
  std::uint8_t candidate_pressure_donor_reach{};
  FaceFluxStorage phi_storage;
  FaceFluxView phi_h_by_a{};
  FaceFluxView trial_flux{};
  PisoIntermediateCertificate intermediate{};
  PressureCorrectionCertificate pressure_certificate{};
  PisoFrozenMomentumStageAuthority authority{};
  BoundaryThermophysicalGhostCertificate base_boundary_thermophysics{};
  PressureReferenceCertificate pressure_reference{};

  static constexpr double absolute_pressure_reference = 101325.0;
  static constexpr RevisionToken target_time = 1901U;
  static constexpr RevisionToken pressure_reference_revision = 1902U;
  static constexpr RevisionToken pressure_reference_authority = 1845U;
  static constexpr FieldId correction_field = 90U;
  static constexpr FieldId candidate_pressure_correction_field = 200U;
  static constexpr StageId live_pressure_donor_stage = 1911U;
  static constexpr StageId candidate_pressure_donor_stage = 1912U;
  static constexpr BdfCoefficients bdf{1.0, -1.0, 0.0, 1U};
  static constexpr StlScanBudget immersed_scan_budget{
      UINT64_C(268435456), UINT64_C(536870912), UINT64_C(4000000),
      UINT64_C(10000), 1U};

 private:
  static std::array<TriangleInput, 12U> immersed_cube(
      bool touches_x_min, bool touches_x_max) noexcept {
    const double x_min = touches_x_min ? -0.25 : 0.25;
    const double x_max = touches_x_max ? 1.25 : 0.75;
    const double inner_x_min = touches_x_max ? 0.55 : x_min;
    const double inner_x_max = touches_x_min ? 0.45 : x_max;
    const Real3 a{inner_x_min, 0.25, 0.25};
    const Real3 b{inner_x_min, 0.25, 0.75};
    const Real3 c{inner_x_min, 0.75, 0.25};
    const Real3 d{inner_x_min, 0.75, 0.75};
    const Real3 e{inner_x_max, 0.25, 0.25};
    const Real3 f{inner_x_max, 0.25, 0.75};
    const Real3 g{inner_x_max, 0.75, 0.25};
    const Real3 h{inner_x_max, 0.75, 0.75};
    return {TriangleInput{a, d, b}, TriangleInput{a, c, d},
            TriangleInput{e, f, h}, TriangleInput{e, h, g},
            TriangleInput{a, b, f}, TriangleInput{a, f, e},
            TriangleInput{c, g, h}, TriangleInput{c, h, d},
            TriangleInput{a, e, g}, TriangleInput{a, g, c},
            TriangleInput{b, d, h}, TriangleInput{b, h, f}};
  }

  static std::uint64_t mix_activity(std::uint64_t hash,
                                    std::uint64_t value) noexcept {
    hash ^= value;
    hash *= UINT64_C(1099511628211);
    return hash;
  }

  bool make_continuity_activity() {
    const Int3 cells = patch.cells;
    const Int3 x_faces{cells.x + 1, cells.y, cells.z};
    const Int3 y_faces{cells.x, cells.y + 1, cells.z};
    const Int3 z_faces{cells.x, cells.y, cells.z + 1};
    const auto count = [](Int3 shape) noexcept {
      return static_cast<std::size_t>(shape.x) *
             static_cast<std::size_t>(shape.y) *
             static_cast<std::size_t>(shape.z);
    };
    const auto offset = [](Int3 shape, Int3 value) noexcept {
      return static_cast<std::size_t>(value.x) +
             static_cast<std::size_t>(shape.x) *
                 (static_cast<std::size_t>(value.y) +
                  static_cast<std::size_t>(shape.y) *
                      static_cast<std::size_t>(value.z));
    };
    const Span<const std::uint8_t> region = immersed_topology.region();
    if (region.data == nullptr || region.size != count(cells)) return false;
    continuity_cell_activity.assign(region.data, region.data + region.size);
    continuity_x_activity.assign(count(x_faces), 1U);
    continuity_y_activity.assign(count(y_faces), 1U);
    continuity_z_activity.assign(count(z_faces), 1U);
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (continuity_cell_activity[offset(cells, cell)] != 0U) continue;
          continuity_x_activity[offset(x_faces, cell)] = 0U;
          continuity_x_activity[offset(x_faces, {x + 1, y, z})] = 0U;
          continuity_y_activity[offset(y_faces, cell)] = 0U;
          continuity_y_activity[offset(y_faces, {x, y + 1, z})] = 0U;
          continuity_z_activity[offset(z_faces, cell)] = 0U;
          continuity_z_activity[offset(z_faces, {x, y, z + 1})] = 0U;
        }
    const Span<const ImmersedLink> links = immersed_topology.links();
    for (std::size_t index = 0U; index < links.size; ++index) {
      const ImmersedLink& link = links.data[index];
      const Int3 cell = link.fluid_local_index;
      switch (link.direction) {
        case ImmersedFaceDirection::x_negative:
          continuity_x_activity[offset(x_faces, cell)] = 0U;
          break;
        case ImmersedFaceDirection::x_positive:
          continuity_x_activity[offset(
              x_faces, {cell.x + 1, cell.y, cell.z})] = 0U;
          break;
        case ImmersedFaceDirection::y_negative:
          continuity_y_activity[offset(y_faces, cell)] = 0U;
          break;
        case ImmersedFaceDirection::y_positive:
          continuity_y_activity[offset(
              y_faces, {cell.x, cell.y + 1, cell.z})] = 0U;
          break;
        case ImmersedFaceDirection::z_negative:
          continuity_z_activity[offset(z_faces, cell)] = 0U;
          break;
        case ImmersedFaceDirection::z_positive:
          continuity_z_activity[offset(
              z_faces, {cell.x, cell.y, cell.z + 1})] = 0U;
          break;
      }
    }
    std::uint64_t local = UINT64_C(1469598103934665603);
    const auto hash = [&](const std::vector<std::uint8_t>& values) {
      local = mix_activity(local, values.size());
      for (std::uint8_t value : values)
        local = mix_activity(local, value);
    };
    hash(continuity_cell_activity);
    hash(continuity_x_activity);
    hash(continuity_y_activity);
    hash(continuity_z_activity);
    local = mix_activity(local, immersed_topology.fingerprint());
    continuity_local_fingerprint = local == 0U ? 1U : local;
    std::uint64_t collective = UINT64_C(1469598103934665603);
    collective = mix_activity(collective,
                              immersed_topology.geometry_fingerprint());
    collective = mix_activity(collective,
                              immersed_topology.surface_fingerprint());
    collective = mix_activity(collective,
                              immersed_topology.geometry_revision());
    continuity_collective_fingerprint =
        collective == 0U ? 1U : collective;
    return true;
  }

  PressureContinuityActivityView continuity_activity() const noexcept {
    return {{continuity_cell_activity.data(),
             continuity_cell_activity.size()},
            {continuity_x_activity.data(), continuity_x_activity.size()},
            {continuity_y_activity.data(), continuity_y_activity.size()},
            {continuity_z_activity.data(), continuity_z_activity.size()},
            continuity_local_fingerprint,
            continuity_collective_fingerprint};
  }

  BoundaryResolvedValues resolved_boundary_values() const noexcept {
    return {{boundary_scalar_values.data(), boundary_scalar_values.size()},
            {boundary_vector_values.data(), boundary_vector_values.size()},
            {boundary_normal_gradient_values.data(),
             boundary_normal_gradient_values.size()}};
  }

  bool initialize_flux_history() {
    FieldRegistry registry;
    FieldSchema schema;
    if (!registry.declare_field("candidate.boundary.history", 1U, 0U,
                                history_dependency) ||
        !registry.freeze(schema))
      return false;
    const std::array requests{ArenaFieldRequest{
        history_dependency, {1, 1, 1}, {0U},
        FieldLifetime::state_layer}};
    ArenaLayout layout;
    if (!ArenaLayout::compile(schema,
                              {requests.data(), requests.size()}, layout) ||
        !StateLayers::allocate(layout, history_layers) ||
        !AttemptTransaction::create(
            history_layers.field_count(), 1U,
            history_layers.field_count(), history_transaction) ||
        !FaceFluxStorage::allocate_final(patch.cells,
                                         committed_storage) ||
        !committed_authority.claim(1910U, 0U, history_transaction,
                                   committed_writer))
      return false;
    if (!history_transaction.begin(history_layers) ||
        !history_transaction.revise_trial(history_dependency))
      return false;
    const RevisionDependency dependency{
        AttemptTransaction::field_revision_source(history_dependency),
        history_transaction.trial_revision(history_dependency)};
    PendingFaceFluxView pending;
    if (!committed_writer.begin_pending(history_transaction,
                                        committed_storage, pending))
      return false;
    OwnedField history_density =
        make_field(240U, patch.cells, 1U, 2U, 1911U, 2911U);
    OwnedField history_velocity =
        make_field(241U, patch.cells, 3U, 2U, 1912U, 2912U);
    fill(history_density, 1.0);
    fill(history_velocity, 0.0);
    const std::array<ConstFieldView, 2U> reads{
        as_const(history_density.view), as_const(history_velocity.view)};
    const KernelInvocation invocation{
        {reads.data(), reads.size()}, {}, {{0, 0, 0}, patch.cells},
        0U, 0U, 1U, 0U, nullptr};
    const std::array dependencies{dependency};
    if (!reconstruct_mass_flux(equations.kernels(), invocation, pending) ||
        !committed_writer.publish_pending(
            {dependencies.data(), dependencies.size()}, pending) ||
        !history_transaction.collective_finish(communicator_, Status{}))
      return false;
    return static_cast<bool>(committed_writer.committed(
        committed_storage, committed_accepted));
  }

  static void fill_face_flux(FaceFluxView flux, double value) {
    for (FaceFieldView face : {flux.x, flux.y, flux.z})
      for (std::int32_t z = 0; z < face.extents.z; ++z)
        for (std::int32_t y = 0; y < face.extents.y; ++y)
          for (std::int32_t x = 0; x < face.extents.x; ++x)
            face.unchecked({x, y, z}) = value;
  }

  MPI_Comm communicator_{MPI_COMM_NULL};
  FieldId history_dependency{};
  StateLayers history_layers;
  AttemptTransaction history_transaction;
  FaceFluxStorage committed_storage;
  FinalFaceFluxAuthority committed_authority;
  FinalFaceFluxWriter committed_writer;
  ConstFaceFluxView committed_accepted{};
};

}  // namespace hundun::v04::test
