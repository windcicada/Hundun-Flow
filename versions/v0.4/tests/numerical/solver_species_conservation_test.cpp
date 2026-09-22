// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_flow.hpp"
#include "../../src/solver_species_guess_detail.hpp"
#include "../../src/solver_ibm_scalar_transport_detail.hpp"
#include "../../src/solver_cold.hpp"
#include "../../src/solver_mixture_rows_detail.hpp"
#include "../../src/solver_mixture_bound_detail.hpp"
#include "../../src/solver_mixture_step_detail.hpp"
#include "../../src/solver_mixture_transport_detail.hpp"
#include "../../src/solver_statistical_detail.hpp"
#include "../../src/models_esf_detail.hpp"
#include "../../src/core_esf_energy_detail.hpp"
#include "../support/ibm_force_fixture.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <= 1.0e-13 * std::max(1.0, std::abs(expected));
}

struct OwnedField {
  std::vector<double> bytes;
  FieldView view{};
};

struct OwnedFaceField {
  std::vector<double> bytes;
  FaceFieldView view{};
};

bool output_is(const OwnedField&, const OwnedField&, const OwnedField&, double);
bool faces_are(const OwnedFaceField&, const OwnedFaceField&, const OwnedFaceField&, double);
bool same_certificate(const EquationAssemblyCertificate&, const EquationAssemblyCertificate&);

OwnedField make_field(FieldId id, Int3 cells, std::uint8_t components,
                      std::uint8_t ghosts, RevisionToken revision) {
  OwnedField field;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  field.bytes.assign(nx * ny * nz * components, 0.0);
  field.view.base = field.bytes.data() + ghosts + ghosts * nx + ghosts * nx * ny;
  field.view.interior = cells;
  field.view.ghosts = {ghosts, ghosts, ghosts};
  field.view.components = components;
  field.view.stride_y = nx;
  field.view.stride_z = nx * ny;
  field.view.component_stride = nx * ny * nz;
  field.view.field = id;
  field.view.revision = revision;
  field.view.storage_identity = 2000U + id;
  field.view.revision_domain = 9201U;
  return field;
}

OwnedFaceField make_face_field(CartesianAxis axis, Int3 cells,
                               StorageIdentity identity) {
  OwnedFaceField field;
  Int3 extents = cells;
  if (axis == CartesianAxis::x) {
    ++extents.x;
  } else if (axis == CartesianAxis::y) {
    ++extents.y;
  } else {
    ++extents.z;
  }
  field.bytes.assign(static_cast<std::size_t>(extents.x) * extents.y *
                         extents.z,
                     0.0);
  field.view = {field.bytes.data(),
                extents,
                static_cast<std::size_t>(extents.x),
                static_cast<std::size_t>(extents.x) * extents.y,
                axis,
                identity,
                9301U};
  return field;
}

void fill_field(OwnedField& field, double value) {
  std::fill(field.bytes.begin(), field.bytes.end(), value);
}

CartesianMeshSpec mesh_spec(std::int32_t n) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {n, n, n};
  mesh.minimum_spacing = {1.0 / n, 1.0 / n, 1.0 / n};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = static_cast<std::uint64_t>(n) * n * n;
  mesh.limits.max_memory_bytes_per_rank = 1U << 28U;
  return mesh;
}

SpeciesThermophysicalSpec thermophysical_species(std::string_view name,
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

ThermophysicalSpec thermo_spec() {
  ThermophysicalSpec value;
  value.data_file = "analytic.d";
  value.minimum_temperature = 200.0;
  value.maximum_temperature = 2000.0;
  value.temperature_relative_tolerance = 1.0e-12;
  value.maximum_temperature_iterations = 64U;
  value.closed_mass_relative_tolerance = 1.0e-12;
  value.maximum_closed_mass_iterations = 32U;
  value.maximum_closed_mass_relative_step = 0.2;
  value.species.push_back(thermophysical_species("species_a", 28.0));
  value.species.push_back(thermophysical_species("species_b", 32.0));
  value.species.back().nasa7_low[0U] = 4.25;
  value.species.back().nasa7_high[0U] = 4.25;
  return value;
}

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

struct ProductionFixture {
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

bool make_production_fixture(std::int32_t n, ProductionFixture& out,
                             bool mixture = false, int inlet_face = -1,
                             bool physical_material = true,
                             double mixing_cell_volume = 0,
                             MPI_Comm communicator = MPI_COMM_SELF,
                             bool interval_sources = false) {
  const auto checked = [](Status status, const char* stage) {
    if (!status) std::cerr << stage << " status=" << unsigned(status.code) << "/" << status.detail << '\n';
    return bool(status);
  };
  CartesianMeshSpec mesh = mesh_spec(n);
  if (mixing_cell_volume > 0) {
    const double dx=std::cbrt(mixing_cell_volume);
    mesh.upper={n*dx,n*dx,n*dx};
    mesh.minimum_spacing={dx,dx,dx};
  }
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x14a5ca1aU;
  model.pressure_reference = PressureReferenceKind::closed_mass;
  model.transported_scalars = {
      {"species_a", TransportedScalarRole::species, 0.5, 2.0},
      {"tracer", TransportedScalarRole::passive_scalar, 2.0, 4.0}};
  if (mixture)
    model.transported_scalars[1] =
        {"species_b",TransportedScalarRole::species,.5,2.0};
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::periodic;
    face.thermal_kind = BoundaryKind::none;
  }
  if (inlet_face >= 0) {
    model.pressure_reference = PressureReferenceKind::boundary_absolute;
    auto& inlet = model.boundaries[static_cast<std::size_t>(inlet_face)];
    auto& outlet = model.boundaries[static_cast<std::size_t>(inlet_face ^ 1)];
    inlet.flow_kind = BoundaryKind::velocity_inlet;
    inlet.temperature = 300.0;
    (inlet_face/2==0 ? inlet.velocity.x : inlet_face/2==1 ? inlet.velocity.y : inlet.velocity.z) =
        inlet_face%2 ? -1.0 : 1.0;
    outlet.flow_kind = BoundaryKind::pressure_outlet;
    outlet.pressure = 101325.0;
    for (auto& face : model.boundaries) {
      if (face.flow_kind == BoundaryKind::periodic) continue;
      for (const auto& scalar : model.transported_scalars)
        face.scalars.push_back({scalar.stable_name,
            &face == &inlet ? ScalarBoundaryKind::dirichlet : ScalarBoundaryKind::zero_gradient,
            0.25 + (inlet_face%2 ? 0.1 : 0.0), ScalarBoundaryKind::zero_gradient, 0.0});
    }
  }
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::central2;
  model.schemes.species = ConvectionScheme::central2;
  model.schemes.passive_scalar = ConvectionScheme::central2;
  if (mixture) model.schemes.species = ConvectionScheme::tvd2;

  FieldRegistry registry;
  FieldId density = 0U;
  FieldId velocity = 0U;
  FieldId pressure = 0U;
  FieldId enthalpy = 0U;
  FieldId temperature = 0U;
  if (!registry.require_field("rho", 1U, 2U, density) ||
      density != kDensity ||
      !registry.require_field("U", 3U, 2U, velocity) ||
      velocity != kVelocity ||
      !registry.require_field("pi", 1U, 2U, pressure) ||
      pressure != kPressure ||
      !registry.require_field("h", 1U, 2U, enthalpy) ||
      enthalpy != kEnthalpy ||
      !registry.require_field("T", 1U, 2U, temperature) ||
      temperature != kTemperature ||
      !checked(CartesianGeometryCompiler::compile(communicator, mesh, {}, out.geometry,
                                          out.patch),"scalar fixture geometry") ||
      !checked(BoundaryCompiler::compile(communicator, model, out.geometry, out.patch,
                                 registry, out.boundary, out.schemes,
                                 out.time), "scalar fixture boundary")) {
    return false;
  }

  ThermophysicalSpec thermophysics = thermo_spec();
  if (mixture)
    thermophysics.species.push_back(thermophysical_species("species_c",30.0));
  if (!ThermodynamicsPlan::compile(
          thermophysics,
          {model.transported_scalars.data(), model.transported_scalars.size()},
          out.thermodynamics) ||
      !TransportPlan::compile(thermophysics, out.thermodynamics,
                              out.transport)) {
    return false;
  }
  const std::array<FieldId, 14U> declared{
      kDensity,          kVelocity,       kPressure,
      kEnthalpy,         kTemperature,    kSpecies,
      kPassive,          kEffectiveViscosity,
      kCompressibility,  kVelocityGradient, 10U, 11U, 12U, 13U};
  if (!out.contributions.configure({declared.data(),
          mixing_cell_volume>0 || interval_sources ? declared.size() : declared.size()-4})) return false;
  if (mixing_cell_volume > 0) {
    ContributionSpec mixing;
    mixing.conserved_quantity=kSpecies;
    mixing.units.si_exponents={1,-3,-1,0,0,0,0};
    mixing.stage=9;
    const std::array<FieldId,2> reads{kSpecies,kDensity};
    mixing.reads={reads.data(),reads.size()};
    mixing.explicit_source=10U;
    mixing.implicit_diagonal=11U;
    mixing.supplies_implicit_diagonal=true;
    if (!checked(out.contributions.register_contribution(mixing),"IEM registration")) return false;
    mixing.conserved_quantity=kEnthalpy;
    mixing.units.si_exponents={1,-1,-3,0,0,0,0};
    const std::array<FieldId,2> heat_reads{kEnthalpy,kDensity};
    mixing.reads={heat_reads.data(),heat_reads.size()};
    mixing.explicit_source=12U; mixing.implicit_diagonal=13U;
    if (!checked(out.contributions.register_contribution(mixing),"enthalpy IEM registration")) return false;
  }
  if (interval_sources) {
    ContributionSpec source;
    source.conserved_quantity=kSpecies;
    source.units.si_exponents={1,-3,-1,0,0,0,0};
    const std::array<FieldId,2> reads{kSpecies,kDensity};
    source.reads={reads.data(),reads.size()};
    for (unsigned i=0;i<2;++i) {
      source.stage=i ? 47U : 31U;
      source.explicit_source=10U+i;
      if (!checked(out.contributions.register_contribution(source),"interval source registration")) return false;
    }
  }
  if (!checked(out.contributions.freeze(),"scalar contribution freeze")) return false;

  const std::array<ScalarEquationSpec, 2U> scalar_specs{{
      {kSpecies, TransportedScalarRole::species, 0.5, 2.0},
      {kPassive, mixture ? TransportedScalarRole::species : TransportedScalarRole::passive_scalar,
       mixture ? .5 : 2.0, mixture ? 2.0 : 4.0},
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
  spec.pressure_reference = model.pressure_reference;
  spec.physical_inlet_material = inlet_face >= 0 && physical_material;
  spec.scalars = {scalar_specs.data(), scalar_specs.size()};
  spec.closed_mass_service_stage = inlet_face >= 0 ? 0U : 1U;
  spec.maximum_cells_per_rank = static_cast<std::size_t>(n) * n * n;
  spec.unity_lewis_total_enthalpy = mixture;
  return checked(EquationPlanSet::compile(
      communicator, out.schemes, out.geometry, out.patch, out.boundary,
      out.contributions, out.thermodynamics, out.transport, spec,
      out.equations), "scalar fixture equations");
}

struct FinalFluxFixture {
  FieldId dependency{};
  StateLayers layers;
  AttemptTransaction transaction;
  FaceFluxStorage storage;
  FinalFaceFluxWriter writer;
};

bool make_linear_final_flux(const CartesianKernelPlan& kernels, Int3 cells,
                            FinalFluxFixture& fixture,
                            ConstFaceFluxView& committed, bool stationary = false,
                            bool uniform = false,
                            double prescribed_speed = std::numeric_limits<double>::quiet_NaN(),
                            double periodic_amplitude = 0) {
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
  if (!ArenaLayout::compile(schema, {requests.data(), requests.size()},
                            layout) ||
      !StateLayers::allocate(layout, fixture.layers) ||
      !AttemptTransaction::create(fixture.layers.field_count(), 1U,
                                  fixture.layers.field_count(),
                                  fixture.transaction) ||
      !FaceFluxStorage::allocate_final(cells, fixture.storage)) {
    return false;
  }
  FinalFaceFluxAuthority authority;
  if (!authority.claim(41U, 0U, fixture.transaction, fixture.writer) ||
      !fixture.transaction.begin(fixture.layers) ||
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

  OwnedField rho = make_field(80U, cells, 1U, 2U, 501U);
  OwnedField velocity = make_field(81U, cells, 3U, 2U, 502U);
  const double spacing = 1.0 / static_cast<double>(cells.x);
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const Int3 cell{i, j, k};
        rho.view.unchecked(cell, 0U) = 1.0;
        velocity.view.unchecked(cell, 0U) =
            std::isfinite(prescribed_speed) ? prescribed_speed :
            stationary ? 0.0 : uniform ? 1.0/(spacing*spacing)
                                     : (static_cast<double>(i) + 0.5) * spacing;
        if (periodic_amplitude!=0)
          velocity.view.unchecked(cell,0)+=periodic_amplitude*std::sin(
              2*std::acos(-1.)*detail::centre_coordinate(kernels,CartesianAxis::x,i));
        velocity.view.unchecked(cell, 1U) = 0.0;
        velocity.view.unchecked(cell, 2U) = 0.0;
      }
    }
  }
  const std::array<ConstFieldView, 2U> reads{as_const(rho.view),
                                             as_const(velocity.view)};
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

// Compare the complete correction matrix with an independent public residual
// difference, keeping history, material and limiter coefficients fixed.
bool check_scalar_correction_rows(const ProductionFixture& fixture,
    ConstFieldView q, EquationStateView state, EquationMaterialView material,
    Span<const EquationContributionView> contributions,
    EquationAssemblyContext context, EquationSystemView system, bool passive=false,
    bool specific=false, bool heat=false, bool statistical=false) {
  const auto cells=fixture.patch.cells;
  const auto boundary_stage=heat ? BoundaryStage::enthalpy : BoundaryStage::scalar;
  const auto count=std::size_t(cells.x)*cells.y*cells.z;
  auto shifted=make_field(q.field,cells,1,2,q.revision+1);
  auto direction=make_field(95,cells,1,2,1201);
  auto shifted_mid=make_field(q.field,cells,1,2,1211);
  ScalarMidpointView midpoint;
  auto variation=make_field(96,cells,1,1,1202);
  auto mass_divergence=make_field(97,cells,1,0,1203);
  std::vector<detail::ColdPressureRow> rows(count);
  std::vector<double> original(count);
  const auto assemble=[&](EquationAssemblyCertificate& certificate) {
    if (heat) return detail::StatisticalEnthalpy::assemble(fixture.equations.enthalpy(),state,
        material.enthalpy_diffusivity,contributions,context,system,certificate);
    if (statistical) return detail::StatisticalSpecies::assemble(fixture.equations.species(),0,
        state,material,contributions,context,system,certificate);
    return passive
        ? assemble_scalar(fixture.equations.scalars(),0,state,material,contributions,context,system,certificate)
        : assemble_species(fixture.equations.species(),0,state,material,contributions,context,system,certificate);
  };
  EquationAssemblyCertificate certificate;
  auto status=assemble(certificate);
  if (status && specific) status=detail::close_frozen_density_scalar_rows(fixture.equations.kernels(),
      fixture.boundary,boundary_stage,q.field,
      context.mixture_transport ? ConvectionScheme::tvd2 : ConvectionScheme::central2,
      state.velocity.trial,context,state.density,q,system,certificate,variation.view,{rows.data(),rows.size()});
  if (status && !specific) status=detail::close_mixture_scalar_rows(fixture.equations.kernels(),
      fixture.boundary,boundary_stage,q.field,
      context.mixture_transport ? ConvectionScheme::tvd2 : ConvectionScheme::central2,
      state.velocity.trial,context,system,certificate,variation.view,{rows.data(),rows.size()});
  if (!status) {std::cerr<<"correction rows status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
  if (specific) {
    KernelInvocation call{{},{&mass_divergence.view,1},{{0,0,0},cells},0,0,1,context.face_flux};
    status=cartesian_face_divergence(fixture.equations.kernels(),context.mass_flux,call);
    if (!status) return false;
    double correction{};
    for (const auto& row:rows) correction=std::max(correction,std::abs(row.rhs));
    std::cout<<"frozen density uniform scalar correction_rhs="<<correction<<'\n';
    const auto& row=rows[3+cells.x*(3+cells.y*3)];
    std::cout<<"frozen density reference row "<<std::setprecision(17)<<row.diagonal;
    for (double value:row.neighbour) std::cout<<' '<<value;
    std::cout<<' '<<row.rhs<<'\n';
    if (!expect(correction<1e-11,"frozen-density scalar preserves a uniform state with divergent mass flux")) return false;
    auto other_density=make_field(kDensity,cells,1,0,1204);
    fill_field(other_density,1.01);
    auto changing_density=state.density;
    changing_density.accepted=as_const(other_density.view);
    const auto before=rows[0].diagonal;
    const auto rejected=detail::close_frozen_density_scalar_rows(fixture.equations.kernels(),
        fixture.boundary,boundary_stage,q.field,ConvectionScheme::central2,
        state.velocity.trial,context,changing_density,q,system,certificate,variation.view,{rows.data(),rows.size()});
    if (!expect(rejected.code==StatusCode::invalid_plan && rows[0].diagonal==before,
        "frozen-density closure checks its time-layer contract before modifying rows")) return false;
    // Compare the complete frozen-density equation with its conservative
    // representation on the paired mass, including the matrix derivative.
    auto basis=make_field(kDensity,cells,1,0,1205);
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
      const Int3 c{x,y,z};
      basis.view.unchecked(c,0)=state.density.accepted.unchecked(c,0)-
          mass_divergence.view.unchecked(c,0)/context.bdf.a0;
    }
    const auto original_density=state.density;
    state.density.trial=as_const(basis.view);
    status=assemble(certificate);
    std::vector<detail::ColdPressureRow> conservative(count);
    if(status)status=detail::close_mixture_scalar_rows(fixture.equations.kernels(),
        fixture.boundary,boundary_stage,q.field,
        context.mixture_transport ? ConvectionScheme::tvd2 : ConvectionScheme::central2,
        state.velocity.trial,context,system,certificate,variation.view,
        {conservative.data(),conservative.size()});
    if(!status)return false;
    double equivalence{};
    for(std::size_t i=0;i<count;++i) {
      equivalence=std::max(equivalence,std::abs(conservative[i].diagonal-rows[i].diagonal)/std::max(1.,std::abs(rows[i].diagonal)));
      equivalence=std::max(equivalence,std::abs(conservative[i].rhs-rows[i].rhs)/std::max(1.,std::abs(rows[i].diagonal)));
      for(unsigned f=0;f<6;++f)
        equivalence=std::max(equivalence,std::abs(conservative[i].neighbour[f]-rows[i].neighbour[f])/std::max(1.,std::abs(rows[i].diagonal)));
    }
    std::cout<<"frozen density conservative basis full-row difference="<<equivalence<<'\n';
    if(!expect(equivalence<1e-11,"frozen-density full matrix and residual equal conservative paired-mass equation"))return false;
    state.density=original_density;
    status=assemble(certificate);
    if(!status)return false;
  }
  const double pi=std::acos(-1.);
  for (int z=-2;z<cells.z+2;++z) for (int y=-2;y<cells.y+2;++y) for (int x=-2;x<cells.x+2;++x)
    direction.view.unchecked({x,y,z},0)=-.02*q.unchecked({x,y,z},0)*
        (1+.2*std::sin(2*pi*(x+.5)/cells.x)+
         .1*std::cos(2*pi*(y+.5)/cells.y)-.1*std::sin(2*pi*(z+.5)/cells.z));
  const auto activity=context.immersed_interface
      ? context.immersed_interface->cell_activity() : Span<const std::uint8_t>{};
  std::size_t index{};
  for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y) for (int x=0;x<cells.x;++x,++index) {
    original[index]=system.residual.unchecked({x,y,z},0);
    if (specific) original[index]-=q.unchecked({x,y,z},0)*
        mass_divergence.view.unchecked({x,y,z},0)*detail::cell_volume(fixture.equations.kernels(),{x,y,z});
    if (activity.size && activity.data[index]==0) direction.view.unchecked({x,y,z},0)=0;
  }
  status=apply_homogeneous_scalar_boundary_ghosts(boundary_stage,
      fixture.boundary,q.field,direction.view,1,state.velocity.trial);
  if (!status) {std::cerr<<"correction boundary status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
  for (int z=-2;z<cells.z+2;++z) for (int y=-2;y<cells.y+2;++y) for (int x=-2;x<cells.x+2;++x)
    shifted.view.unchecked({x,y,z},0)=q.unchecked({x,y,z},0)+direction.view.unchecked({x,y,z},0);
  if(context.scalar_midpoint) {
    for(int z=-2;z<cells.z+2;++z)for(int y=-2;y<cells.y+2;++y)for(int x=-2;x<cells.x+2;++x) {
      const Int3 c{x,y,z};
      shifted_mid.view.unchecked(c,0)=context.scalar_midpoint->value.unchecked(c,0)+.5*direction.view.unchecked(c,0);
    }
    midpoint.value=as_const(shifted_mid.view);context.scalar_midpoint=&midpoint;
  }
  const auto histories=heat ? Span<const PrimitiveHistory>{&state.enthalpy,1}
      : passive ? state.passive_scalars : state.independent_species;
  std::vector<PrimitiveHistory> shifted_histories(histories.data,histories.data+histories.size);
  shifted_histories[0].trial=as_const(shifted.view);
  if (heat) state.enthalpy=shifted_histories[0];
  else if (passive) state.passive_scalars={shifted_histories.data(),shifted_histories.size()};
  else state.independent_species={shifted_histories.data(),shifted_histories.size()};
  status=assemble(certificate);
  if (!status) {std::cerr<<"shifted scalar status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
  double error{},scale{1.}; index=0;
  for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y) for (int x=0;x<cells.x;++x,++index) {
    const Int3 cell{x,y,z};
    const Int3 neighbors[]={{x-1,y,z},{x+1,y,z},{x,y-1,z},{x,y+1,z},{x,y,z-1},{x,y,z+1}};
    double action=rows[index].diagonal*direction.view.unchecked(cell,0);
    for (unsigned f=0;f<6;++f) action-=rows[index].neighbour[f]*direction.view.unchecked(neighbors[f],0);
    const double difference=(system.residual.unchecked(cell,0)-original[index])/
        detail::cell_volume(fixture.equations.kernels(),cell)-
        (specific ? shifted.view.unchecked(cell,0)*mass_divergence.view.unchecked(cell,0) : 0);
    error=std::max(error,std::abs(action-difference));
    scale=std::max(scale,std::abs(difference));
  }
  std::cout<<"scalar correction ibm="<<bool(context.immersed_interface)<<" passive="<<passive<<" enthalpy="<<heat
      <<" mixture="<<bool(context.mixture_transport)<<" relative_error="<<error/scale<<'\n';
  return expect(error/scale<1e-11,"complete scalar matrix differentiates the public residual with frozen history");
}

bool test_frozen_scalar_rows() {
  ProductionFixture fixture;
  if (!make_production_fixture(8,fixture,false,0)) return false;
  const auto cells=fixture.patch.cells;
  auto rho=make_field(kDensity,cells,1,2,1301);
  auto q=make_field(kSpecies,cells,1,2,1302);
  auto gamma=make_field(25,cells,1,2,1303);
  auto diagonal=make_field(30,cells,1,0,1304);
  auto rhs=make_field(31,cells,1,0,1305);
  auto residual=make_field(32,cells,1,0,1306);
  auto ax=make_face_field(CartesianAxis::x,cells,1307);
  auto ay=make_face_field(CartesianAxis::y,cells,1308);
  auto az=make_face_field(CartesianAxis::z,cells,1309);
  fill_field(rho,1.); fill_field(q,.25); fill_field(gamma,.1);
  const PrimitiveHistory history{as_const(q.view),as_const(q.view),as_const(q.view)};
  EquationStateView state;
  state.density={as_const(rho.view),as_const(rho.view),as_const(rho.view)};
  state.independent_species={&history,1};
  const auto diffusivity=as_const(gamma.view);
  EquationMaterialView material;material.scalar_mass_diffusivity={&diffusivity,1};
  FinalFluxFixture owner;
  ConstFaceFluxView flux;
  if (!make_linear_final_flux(fixture.equations.kernels(),cells,owner,flux)) return false;
  EquationAssemblyContext context;
  context.dt=.1;context.bdf={10,-10,0,1};context.time=1310;
  context.geometry=fixture.geometry.topology_revision();context.boundary=fixture.boundary.revision();
  context.thermo=fixture.thermodynamics.fingerprint();context.transport=fixture.transport.fingerprint();
  context.contribution_stage=1;context.scope=EquationAssemblyScope::final_conservative;
  context.mass_flux=flux;context.face_flux=flux.revision;
  context.face_flux_authority=flux.certificate.authority();context.face_flux_storage=flux.certificate.storage();
  context.face_flux_revision_domain=flux.certificate.revision_domain();
  return check_scalar_correction_rows(fixture,as_const(q.view),state,material,{},context,
      {diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view},false,true);
}

bool test_passive_midpoint() {
  ProductionFixture fixture;
  if(!make_production_fixture(8,fixture,false))return false;
  const auto cells=fixture.patch.cells;
  auto rho=make_field(kDensity,cells,1,2,1401),oldrho=make_field(kDensity,cells,1,2,1402);
  auto q=make_field(kPassive,cells,1,2,1403),old=make_field(kPassive,cells,1,2,1404);
  auto mid=make_field(kPassive,cells,1,2,1405),gamma=make_field(25,cells,1,2,1406);
  auto diagonal=make_field(30,cells,1,0,1407),rhs=make_field(31,cells,1,0,1408);
  auto residual=make_field(32,cells,1,0,1409);
  auto ax=make_face_field(CartesianAxis::x,cells,1410);
  auto ay=make_face_field(CartesianAxis::y,cells,1411);
  auto az=make_face_field(CartesianAxis::z,cells,1412);
  fill_field(rho,1.3);fill_field(oldrho,1.1);fill_field(gamma,.2);
  for(int z=-2;z<cells.z+2;++z)for(int y=-2;y<cells.y+2;++y)for(int x=-2;x<cells.x+2;++x) {
    const Int3 c{x,y,z};const double X=(x+.5)/8.,Y=(y+.5)/8.;
    q.view.unchecked(c,0)=.3+X*X+.2*Y;
    old.view.unchecked(c,0)=.2+.5*X*X+.1*Y;
    mid.view.unchecked(c,0)=.5*(q.view.unchecked(c,0)+old.view.unchecked(c,0));
  }
  PrimitiveHistory tracer{as_const(q.view),as_const(old.view),as_const(old.view)};
  EquationStateView state;
  state.density={as_const(rho.view),as_const(oldrho.view),as_const(oldrho.view)};
  state.passive_scalars={&tracer,1};
  const auto diffusivity=as_const(gamma.view);
  EquationMaterialView material;material.scalar_mass_diffusivity={&diffusivity,1};
  FinalFluxFixture owner;ConstFaceFluxView flux;
  if(!make_linear_final_flux(fixture.equations.kernels(),cells,owner,flux))return false;
  EquationAssemblyContext context;
  context.dt=.1;context.bdf={10,-10,0,1};context.time=1413;
  context.geometry=fixture.geometry.topology_revision();context.boundary=fixture.boundary.revision();
  context.thermo=fixture.thermodynamics.fingerprint();context.transport=fixture.transport.fingerprint();
  context.contribution_stage=1;context.scope=EquationAssemblyScope::final_conservative;
  context.mass_flux=flux;context.face_flux=flux.revision;
  context.face_flux_authority=flux.certificate.authority();context.face_flux_storage=flux.certificate.storage();
  context.face_flux_revision_domain=flux.certificate.revision_domain();
  ScalarMidpointView midpoint{as_const(mid.view)};context.scalar_midpoint=&midpoint;
  EquationSystemView system{diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view};
  EquationAssemblyCertificate certificate;
  auto status=assemble_scalar(fixture.equations.scalars(),0,state,material,{},context,system,certificate);
  if(!expect(bool(status),"passive CN assembly"))return false;
  bool passed=true;double error{};
  const double dx=1./8,V=dx*dx*dx,D=.2*dx;
  for(int z=1;z<7;++z)for(int y=1;y<7;++y)for(int x=1;x<7;++x) {
    const Int3 c{x,y,z},lo{x-1,y,z},hi{x+1,y,z};
    const double expected=V*(10*(1.3*q.view.unchecked(c,0)-1.1*old.view.unchecked(c,0))-.2*1.5)+
        flux.x.unchecked(hi)*.5*(mid.view.unchecked(c,0)+mid.view.unchecked(hi,0))-
        flux.x.unchecked(c)*.5*(mid.view.unchecked(c,0)+mid.view.unchecked(lo,0));
    error=std::max(error,std::abs(expected-residual.view.unchecked(c,0)));
    passed &= close(diagonal.view.unchecked(c,0),13*V+3*D);
    passed &= close(ax.view.unchecked(c),.5*D);
  }
  passed &= expect(error<1e-13,"CN conservative storage, analytic midpoint convection/diffusion and half matrix");
  passed &= check_scalar_correction_rows(fixture,as_const(q.view),state,material,{},context,system,true);
  // Solve the complete diffusion equation in its invariant Fourier
  // subspace using two independently assembled residuals for each step.
  // Compare with the exact semidiscrete exponential to isolate time order.
  FinalFluxFixture quiet_owner;ConstFaceFluxView quiet;
  if(!make_linear_final_flux(fixture.equations.kernels(),cells,quiet_owner,quiet,true))return false;
  context.mass_flux=quiet;context.face_flux=quiet.revision;
  context.face_flux_authority=quiet.certificate.authority();
  context.face_flux_storage=quiet.certificate.storage();
  context.face_flux_revision_domain=quiet.certificate.revision_domain();
  fill_field(rho,1.);fill_field(oldrho,1.);
  const double pi=std::acos(-1.),lambda=4*.2*64*std::pow(std::sin(pi/8),2);
  std::array<double,3> errors{};
  for(unsigned refinement=0;refinement<3;++refinement) {
    const unsigned steps=4u<<refinement;context.dt=.08/steps;
    context.bdf={1/context.dt,-1/context.dt,0,1};
    double amplitude=.1;
    const auto trial=[&](double value) {
      for(int z=-2;z<cells.z+2;++z)for(int y=-2;y<cells.y+2;++y)for(int x=-2;x<cells.x+2;++x) {
        const Int3 c{x,y,z};const double mode=std::sin(2*pi*(x+.5)/8);
        q.view.unchecked(c,0)=.5+value*mode;
        mid.view.unchecked(c,0)=.5*(q.view.unchecked(c,0)+old.view.unchecked(c,0));
      }
    };
    for(unsigned step=0;step<steps;++step) {
      for(int z=-2;z<cells.z+2;++z)for(int y=-2;y<cells.y+2;++y)for(int x=-2;x<cells.x+2;++x)
        old.view.unchecked({x,y,z},0)=.5+amplitude*std::sin(2*pi*(x+.5)/8);
      trial(0.);
      status=assemble_scalar(fixture.equations.scalars(),0,state,material,{},context,system,certificate);
      if(!status)return false;
      const double r0=residual.view.unchecked({2,2,2},0);
      trial(1.);
      status=assemble_scalar(fixture.equations.scalars(),0,state,material,{},context,system,certificate);
      if(!status)return false;
      amplitude=-r0/(residual.view.unchecked({2,2,2},0)-r0);
      trial(amplitude);
      status=assemble_scalar(fixture.equations.scalars(),0,state,material,{},context,system,certificate);
      if(!status)return false;
      double worst{};
      for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
        worst=std::max(worst,std::abs(residual.view.unchecked({x,y,z},0))/V);
      passed &= expect(worst<1e-12,"CN Fourier solve satisfies every complete scalar row");
    }
    errors[refinement]=std::abs(amplitude-.1*std::exp(-lambda*.08));
  }
  const double order1=std::log2(errors[0]/errors[1]),order2=std::log2(errors[1]/errors[2]);
  passed &= expect(order1>1.9 && order2>1.9,"passive CN has second order time convergence");
  std::cout<<"passive CN temporal orders="<<order1<<','<<order2<<'\n';
  midpoint.value=as_const(residual.view);
  const auto rejected=assemble_scalar(fixture.equations.scalars(),0,state,material,{},context,system,certificate);
  passed &= expect(rejected.code==StatusCode::invalid_plan,"CN rejects invalid midpoint storage");
  std::cout<<"passive CN analytic row error="<<error<<" passed="<<passed<<'\n';
  return passed;
}

bool mixture_face_probe(bool implicit=false) {
  ProductionFixture fixture;
  if (!make_production_fixture(8,fixture)) return false;
  const auto cells=fixture.patch.cells;
  double phi{},physical{};
  unsigned count{},thermal{},allow{};
  while (std::cin>>phi>>physical>>allow>>count>>thermal) {
    if (count>64 || thermal>1 || allow>1) return false;
    std::vector<OwnedField> fields;
    std::vector<ConstFieldView> composition;
    fields.reserve(count+thermal); composition.reserve(count);
    for (unsigned s=0; s<count+thermal; ++s) {
      std::array<double,4> samples{};
      for (auto& value:samples) if (!(std::cin>>value)) return false;
      fields.push_back(make_field(80U+s,cells,1U,2U,501U+s));
      auto& field=fields.back();
      for (int z=-2; z<cells.z+2; ++z) for (int y=-2; y<cells.y+2; ++y)
        for (int x=-2; x<cells.x+2; ++x)
          field.view.unchecked({x,y,z},0U)=samples[std::clamp(x-1,0,3)];
      if (s<count) composition.push_back(as_const(field.view));
    }
    MixtureFaceTransport face;
    const auto status=prepare_cartesian_mixture_face(fixture.equations.kernels(),
        {composition.data(),composition.size()},
        thermal ? as_const(fields.back().view) : ConstFieldView{},
        CartesianAxis::x,{3,3,3},phi,physical,face,allow!=0,
        implicit ? MixtureFlatStencilPolicy::upwind_constraint : MixtureFlatStencilPolicy::ignore_roundoff);
    if (!status) return false;
    std::cout<<std::setprecision(17)<<face.diffusion<<'\n';
  }
  return std::cin.eof();
}

bool test_mixture_bound_repair(MPI_Comm comm=MPI_COMM_SELF) {
  ProductionFixture fixture;
  if(!make_production_fixture(8,fixture,false,-1,true,0,comm))return false;
  const auto cells=fixture.patch.cells;
  // Two-dimensional periodic trace front. Frozen VLS permits the updated
  // transverse neighbour to drain a nearly empty cell through a central face.
  const std::array<double,16> initial{.2,.2,0,.1,0,1e-10,.01,1e-12,
      1e-10,1e-10,1e-8,0,0,1e-8,1e-13,1e-10};
  auto q=make_field(kSpecies,cells,1,2,8901);
  auto gamma=make_field(90,cells,1,2,8902);
  auto mask=make_field(91,cells,1,2,8903);
  fill_field(gamma,8e-6);fill_field(mask,0);
  const auto tile=[](int x,int y){return std::size_t((x+8)%4+4*((y+8)%4));};
  for(int z=-2;z<cells.z+2;++z)for(int y=-2;y<cells.y+2;++y)for(int x=-2;x<cells.x+2;++x)
    q.view.unchecked({x,y,z},0)=initial[tile(x+fixture.patch.begin.x,y+fixture.patch.begin.y)];
  FaceFluxStorage mass_storage,extra_storage;FaceFluxView mass,extra;
  auto status=FaceFluxStorage::allocate_workspace(cells,1,mass_storage);
  if(status)status=FaceFluxStorage::allocate_workspace(cells,1,extra_storage);
  if(status)status=mass_storage.workspace_view(0,8904,mass);
  if(status)status=extra_storage.workspace_view(0,8905,extra);
  for(auto face:{mass.x,mass.y,mass.z})
    for(int z=0;z<face.extents.z;++z)for(int y=0;y<face.extents.y;++y)for(int x=0;x<face.extents.x;++x)
      face.unchecked({x,y,z})=face.axis==CartesianAxis::z ? 0. : 1.;
  const auto species=as_const(q.view);MixtureTransportFaces mixture;
  if(status)status=prepare_cartesian_mixture_transport(fixture.equations.kernels(),
      {&species,1},{},as_const(gamma.view),as_const(mass),extra,8905,mixture,nullptr,
      MixtureFlatStencilPolicy::upwind_constraint);
  if(!expect(bool(status),"periodic trace-front common faces prepare"))return false;
  std::array<double,64> conductance{};
  const auto gather=[&]() {
    std::array<double,64> local{};
    for(int y=0;y<4;++y)for(int x=0;x<4;++x) {
      const Int3 c{x-fixture.patch.begin.x,y-fixture.patch.begin.y,2-fixture.patch.begin.z};
      if(c.x<0 || c.y<0 || c.z<0 || c.x>=cells.x || c.y>=cells.y || c.z>=cells.z)continue;
      for(unsigned axis=0;axis<2;++axis) {
        auto face=c;const auto values=axis==0 ? extra.x : extra.y;
        local[4*tile(x,y)+2*axis]=1e-6+values.unchecked(face);
        (axis==0 ? face.x : face.y)++;
        local[4*tile(x,y)+2*axis+1]=1e-6+values.unchecked(face);
      }
    }
    return MPI_Allreduce(local.data(),conductance.data(),64,MPI_DOUBLE,MPI_SUM,comm)==MPI_SUCCESS;
  };
  if(!gather())return false;
  const auto original_conductance=conductance;
  const auto solve=[&](double shift,double scale) {
    std::array<std::array<long double,17>,16> a{};
    for(int y=0;y<4;++y)for(int x=0;x<4;++x) {
      const auto i=tile(x,y);a[i][i]=10;a[i][16]=10*(shift+scale*initial[i]);
      for(unsigned axis=0;axis<2;++axis) {
        const auto left=tile(x-(axis==0),y-(axis==1));
        const auto right=tile(x+(axis==0),y+(axis==1));
        const double dl=conductance[4*i+2*axis];
        const double dr=conductance[4*i+2*axis+1];
        a[i][i]+=dl+dr;a[i][left]-=.5+dl;a[i][right]+=.5-dr;
      }
    }
    for(unsigned i=0;i<16;++i) {
      unsigned pivot=i;
      for(unsigned j=i+1;j<16;++j)if(std::abs(a[j][i])>std::abs(a[pivot][i]))pivot=j;
      std::swap(a[i],a[pivot]);
      for(unsigned j=i+1;j<16;++j) {
        const auto r=a[j][i]/a[i][i];
        for(unsigned k=i;k<=16;++k)a[j][k]-=r*a[i][k];
      }
    }
    std::array<long double,16> result{};
    for(int i=15;i>=0;--i) {
      auto rhs=a[i][16];for(int j=i+1;j<16;++j)rhs-=a[i][j]*result[j];
      result[i]=rhs/a[i][i];
    }
    return result;
  };
  const auto raw=solve(0,1);
  bool passed=expect(*std::min_element(raw.begin(),raw.end()) < -6e-4,
      "frozen VLS reproduces the GTMC trace-species loss of positivity");
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
    mask.view.unchecked({x,y,z},0)=raw[tile(x+fixture.patch.begin.x,y+fixture.patch.begin.y)]<0 ? 1 : 0;
  HaloEngine halo;const HaloFieldSpec mask_spec{mask.view.field,1,1};HaloTicket ticket;
  status=halo.reserve(comm,fixture.patch,{&mask_spec,1},fixture.boundary.halo_topology());
  if(status)status=halo.begin(140,{&mask.view,1},ticket);
  if(status)status=halo.finish(ticket,{&mask.view,1});
  if(!expect(bool(status),"periodic bound decision mask exchanges"))return false;
  status=detail::constrain_mixture_bounds(fixture.equations.kernels(),as_const(mask.view),
      as_const(gamma.view),as_const(mass),extra,nullptr);
  if(!gather())return false;
  const auto repaired=solve(0,1),dependent=solve(1,-1),enthalpy=solve(300,10);
  long double mass_error{};double affine_error{};
  for(unsigned i=0;i<16;++i) {
    mass_error+=repaired[i]-initial[i];
    affine_error=std::max(affine_error,double(std::abs(repaired[i]+dependent[i]-1)));
    affine_error=std::max(affine_error,double(std::abs(enthalpy[i]-(300+10*repaired[i]))));
  }
  passed &= expect(bool(status) && *std::min_element(repaired.begin(),repaired.end())>0 &&
      *std::max_element(repaired.begin(),repaired.end())<1 && std::abs(mass_error)<1e-14 &&
      affine_error<1e-12 && conductance[4*tile(1,1)]==original_conductance[4*tile(1,1)],
      "local common-face repair preserves bounds, periodic mass, composition and affine enthalpy");
  std::cout<<"mixture_bound initial_min="<<double(*std::min_element(raw.begin(),raw.end()))
      <<" repaired_min="<<double(*std::min_element(repaired.begin(),repaired.end()))
      <<" mass_error="<<double(mass_error)<<" affine_error="<<affine_error<<'\n';
  return passed;
}

bool test_mixture_face_closure() {
  ProductionFixture fixture;
  if (!expect(make_production_fixture(8, fixture),
              "mixture face production geometry compiles")) return false;
  const auto cells = fixture.patch.cells;
  auto first = make_field(80U, cells, 1U, 2U, 501U);
  auto second = make_field(81U, cells, 1U, 2U, 502U);
  bool passed = true;
  // REFERENCE gam_tvd/vls on this normalized three-species stencil selects
  // the common upwind state (.2,.7,.1) in either flow direction.
  for (unsigned axis = 0; axis < 3; ++axis)
    for (double direction : {1.0, -1.0}) {
      for (int z = -2; z < cells.z + 2; ++z)
        for (int y = -2; y < cells.y + 2; ++y)
          for (int x = -2; x < cells.x + 2; ++x) {
            const int coordinate = axis == 0 ? x : axis == 1 ? y : z;
            const int i = direction > 0 ? coordinate : 5 - coordinate;
            first.view.unchecked({x,y,z},0U) = i <= 1 ? 0.0 : i == 2 ? .2 : .5;
            second.view.unchecked({x,y,z},0U) = i <= 1 ? .6 : i == 2 ? .7 : .5;
          }
      std::array<double, 3> face{};
      const std::array<ConstFieldView,2> composition{
          as_const(first.view),as_const(second.view)};
      MixtureFaceTransport transport;
      auto status = prepare_cartesian_mixture_face(fixture.equations.kernels(),
          {composition.data(),composition.size()}, {},
          static_cast<CartesianAxis>(axis), {3,3,3},direction,0.0,transport);
      Int3 lower{3,3,3};
      (axis == 0 ? lower.x : axis == 1 ? lower.y : lower.z)--;
      for (unsigned s = 0; s < 2; ++s)
        face[s] = transport.lower_weight*composition[s].unchecked(lower,0U) +
            (1.0-transport.lower_weight)*composition[s].unchecked({3,3,3},0U);
      face[2] = 1.0 - face[0] - face[1];
      std::cout << "mixture_face axis=" << axis << " direction=" << direction
                << " Y=" << face[0] << ',' << face[1] << ',' << face[2] << '\n';
      passed &= expect(bool(status) && face[0] >= 0 && face[1] >= 0 &&
                           face[2] >= 0 && close(face[0],.2) &&
                           close(face[1],.7) && close(face[2],.1) &&
                           close(transport.diffusion,.5),
                       "transported face composition includes a bounded dependent species");
    }
  auto thermal = make_field(82U,cells,1U,2U,503U);
  for (int z=-2; z<cells.z+2; ++z) for (int y=-2; y<cells.y+2; ++y)
    for (int x=-2; x<cells.x+2; ++x) {
      first.view.unchecked({x,y,z},0U)=.1+.01*x;
      second.view.unchecked({x,y,z},0U)=.2;
      thermal.view.unchecked({x,y,z},0U)=300.0+(x%2 ? 1e-12 : -1e-12);
    }
  const std::array<ConstFieldView,2> linear{
      as_const(first.view),as_const(second.view)};
  MixtureFaceTransport smooth;
  auto smooth_status=prepare_cartesian_mixture_face(fixture.equations.kernels(),
      {linear.data(),linear.size()},as_const(thermal.view),CartesianAxis::x,
      {3,3,3},1.0,1e-6,smooth);
  passed &= expect(bool(smooth_status) && close(smooth.extra_diffusion,0.0),
      "constant composition and thermal roundoff preserve linear scalar transport");
  smooth_status=prepare_cartesian_mixture_face(fixture.equations.kernels(),
      {linear.data(),linear.size()},as_const(thermal.view),CartesianAxis::x,
      {3,3,3},1.0,1e-6,smooth,true,MixtureFlatStencilPolicy::upwind_constraint);
  passed &= expect(bool(smooth_status) && close(smooth.diffusion,.5),
      "implicit constant coordinate retains the REFERENCE VLS upwind constraint");
  smooth_status=prepare_cartesian_mixture_face(fixture.equations.kernels(),
      {linear.data(),linear.size()},as_const(thermal.view),CartesianAxis::x,
      {3,3,3},1.0,1e-6,smooth,true,static_cast<MixtureFlatStencilPolicy>(255));
  passed &= expect(!smooth_status && close(smooth.diffusion,.5),
      "invalid flat-stencil policy rejects atomically");
  // The coupled species solver resolves composition against the natural
  // unit mass-fraction scale. A trace below that FP64 resolution must have
  // the same neutral limiter contribution as its constant complement.
  for(double amplitude : {1e-15,1e-8}) {
    for(int z=-2;z<cells.z+2;++z) for(int y=-2;y<cells.y+2;++y)
      for(int x=-2;x<cells.x+2;++x)
        second.view.unchecked({x,y,z},0U)=amplitude*(x%2 ? 1.0 : 2.0);
    smooth_status=prepare_cartesian_mixture_face(fixture.equations.kernels(),
        {linear.data(),linear.size()},as_const(thermal.view),CartesianAxis::x,
        {3,3,3},1.0,1e-6,smooth);
    passed &= expect(bool(smooth_status) && close(smooth.extra_diffusion,
        amplitude<1e-14 ? 0.0 : .5-1e-6),
        "resolved species extrema limit flux; composition roundoff is neutral");
  }
  // The native mean-reactor plateau: an upstream gradient is resolved,
  // while the two face-adjacent compositions differ by at most one ULP.
  // All signs must give the zero-middle-gradient upwind limit, so the
  // common coefficient cannot inject finite changes into the thermal flux.
  const double plateau=.24900199739488843;
  for (unsigned axis=0;axis<3;++axis) for(double direction : {1.,-1.})
    for(double middle : {std::nextafter(plateau,0.),plateau,std::nextafter(plateau,1.)}) {
      for(int z=-2;z<cells.z+2;++z)for(int y=-2;y<cells.y+2;++y)
        for(int x=-2;x<cells.x+2;++x) {
          const int coordinate=axis==0 ? x : axis==1 ? y : z;
          const int i=direction>0 ? coordinate : 5-coordinate;
          first.view.unchecked({x,y,z},0)=i<=1 ? plateau+1.3747336602421001e-13 :
              i==3 ? middle : plateau;
          thermal.view.unchecked({x,y,z},0)=1.+.1*coordinate;
        }
      const auto q=as_const(first.view);
      const auto status=prepare_cartesian_mixture_face(fixture.equations.kernels(),
          {&q,1},as_const(thermal.view),static_cast<CartesianAxis>(axis),
          {3,3,3},direction,0.,smooth);
      passed &= expect(bool(status) && close(smooth.diffusion,.5),
          "resolved upstream to FP64 plateau has a sign-stable common face coefficient");
    }
  return passed;
}

bool test_species_storage_increment() {
  ProductionFixture fixture;
  if (!make_production_fixture(8,fixture,false,-1,true,0,MPI_COMM_SELF,true)) return false;
  const auto cells=fixture.patch.cells;
  auto rho=make_field(kDensity,cells,1U,2U,501U);
  auto q=make_field(kSpecies,cells,1U,2U,502U);
  auto old=make_field(kSpecies,cells,1U,2U,503U);
  auto gamma=make_field(25U,cells,1U,2U,504U);
  auto diagonal=make_field(30U,cells,1U,0U,505U);
  auto rhs=make_field(31U,cells,1U,0U,506U);
  auto residual=make_field(32U,cells,1U,0U,507U);
  const double density=.213389, before=.224;
  const double after=std::nextafter(before,1.0);
  fill_field(rho,density); fill_field(q,after); fill_field(old,before);
  fill_field(gamma,1e-6);
  FinalFluxFixture owner;
  ConstFaceFluxView flux;
  if (!make_linear_final_flux(fixture.equations.kernels(),cells,owner,flux,true))
    return false;
  const PrimitiveHistory history{as_const(q.view),as_const(old.view),as_const(old.view)};
  EquationStateView state;
  state.density={as_const(rho.view),as_const(rho.view),as_const(rho.view)};
  state.independent_species={&history,1U};
  const auto diffusion=as_const(gamma.view);
  EquationMaterialView material; material.scalar_mass_diffusivity={&diffusion,1U};
  EquationAssemblyContext context;
  context.dt=1e-5; context.bdf={1/context.dt,-1/context.dt,0,1U}; context.time=701U;
  context.geometry=fixture.geometry.topology_revision(); context.boundary=fixture.boundary.revision();
  context.thermo=fixture.thermodynamics.fingerprint(); context.transport=fixture.transport.fingerprint();
  context.contribution_stage=1U; context.scope=EquationAssemblyScope::final_conservative;
  context.mass_flux=flux; context.face_flux=flux.revision;
  context.face_flux_authority=flux.certificate.authority();
  context.face_flux_storage=flux.certificate.storage();
  context.face_flux_revision_domain=flux.certificate.revision_domain();
  EquationAssemblyCertificate certificate;
  auto ax=make_face_field(CartesianAxis::x,cells,9401U);
  auto ay=make_face_field(CartesianAxis::y,cells,9402U);
  auto az=make_face_field(CartesianAxis::z,cells,9403U);
  const auto status=assemble_species(fixture.equations.species(),0U,state,material,{},context,
      {diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view},certificate);
  if (!status) std::cerr << "storage_status=" << unsigned(status.code) << '/' << status.detail << '\n';
  if (!expect(bool(status),"one ULP species storage assembles")) return false;
  const Int3 c{3,3,3};
  const double volume=hundun::v04::detail::cell_volume(fixture.equations.kernels(),c);
  const double expected=static_cast<double>(static_cast<long double>(density)*
      (static_cast<long double>(after)-before)*context.bdf.a0*volume);
  const double actual=residual.view.unchecked(c,0);
  std::cout << "species_storage actual=" << std::setprecision(17) << actual
            << " expected=" << expected << '\n';
  bool passed=expect(std::abs(actual-expected)<=1e-13*std::abs(expected),
                "BE storage resolves a represented one ULP composition increment");
  auto endpoint=make_field(kSpecies,cells,1U,0U,601U);
  fill_field(endpoint,std::nextafter(after,1.0));
  context.reaction_endpoint=as_const(endpoint.view);
  const auto endpoint_status=assemble_species(fixture.equations.species(),0U,state,material,{},context,
      {diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view},certificate);
  const double split_expected=static_cast<double>(static_cast<long double>(density)*
      (static_cast<long double>(endpoint.view.unchecked(c,0))-before)*context.bdf.a0*volume);
  passed &= expect(bool(endpoint_status) &&
      std::abs(residual.view.unchecked(c,0)-split_expected)<=1e-13*std::abs(split_expected),
      "split species endpoint contributes conservative storage through the common assembler");
  fill_field(endpoint,-.01);
  fill_field(residual,91.0);
  const auto rejected=assemble_species(fixture.equations.species(),0U,state,material,{},context,
      {diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view},certificate);
  passed &= expect(!rejected && residual.view.unchecked(c,0)==91.0,
      "split endpoint bounds are checked before equation output writes");
  // A frozen chemical endpoint can cancel a large interval source while
  // a small parcel source still drives the physical species residual.
  auto first_source=make_field(10U,cells,1U,0U,610U);
  auto second_source=make_field(11U,cells,1U,0U,611U);
  std::array<EquationContributionView,2> sources{};
  for (unsigned i=0;i<2;++i) {
    sources[i].explicit_source_density=as_const(i ? second_source.view : first_source.view);
    sources[i].conserved_quantity=kSpecies;
    sources[i].units.si_exponents={1,-3,-1,0,0,0,0};
    sources[i].stage=i ? 47U : 31U;
    sources[i].explicit_source_field=10U+i;
  }
  fill_field(rho,1.); fill_field(q,.5); fill_field(old,.5);
  context.dt=std::ldexp(1.,-55);
  context.bdf={1/context.dt,-1/context.dt,0,1U};
  context.contribution_stage=31U; context.additional_contribution_stage=47U;
  for (double direction : {-1.,1.}) for (double parcel : {-.25,.25})
    for (bool chemical_first : {false,true}) {
      fill_field(endpoint,.5+direction*.25);
      const double chemical=direction*std::ldexp(1.,53);
      fill_field(first_source,chemical_first ? chemical : parcel);
      fill_field(second_source,chemical_first ? parcel : chemical);
      const auto assembled=assemble_species(fixture.equations.species(),0U,state,material,
          {sources.data(),sources.size()},context,
          {diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view},certificate);
      const double actual=residual.view.unchecked(c,0)/volume;
      std::cout << "species_interval chemical=" << chemical << " parcel=" << parcel
                << " chemical_first=" << chemical_first << " residual=" << actual << '\n';
      passed &= expect(bool(assembled) && actual == -parcel,
          "species residual preserves a small parcel source after chemical cancellation");
    }
  return passed;
}

bool test_mixture_material_floor() {
  ProductionFixture fixture;
  if (!expect(make_production_fixture(8,fixture,true),
      "different-diffusivity face fixture compiles")) return false;
  const auto cells=fixture.patch.cells;
  auto first=make_field(kSpecies,cells,1,2,501);
  auto second=make_field(kPassive,cells,1,2,502);
  auto thermal=make_field(23,cells,1,2,503);
  auto mu=make_field(24,cells,1,2,504);
  auto effective=make_field(25,cells,1,2,505);
  auto minimum=make_field(26,cells,1,2,506);
  auto physical=make_field(27,cells,1,2,507);
  FaceFluxStorage storage;
  FaceFluxView flow,work,transported;
  auto status=FaceFluxStorage::allocate_workspace(cells,3,storage);
  if(status)status=storage.workspace_view(0,801,flow);
  if(status)status=storage.workspace_view(1,802,work);
  if(status)status=storage.workspace_view(2,803,transported);
  if(!expect(bool(status),"different-diffusivity workspaces allocate"))return false;
  const std::array<ConstFieldView,2> composition{as_const(first.view),as_const(second.view)};
  const std::array<FaceFieldView,3> faces{flow.x,flow.y,flow.z};
  const std::array<FaceFieldView,3> output{transported.x,transported.y,transported.z};
  bool passed=true;
  for(double turbulent : {0.0,.0625}) for(double thermal_gamma : {2.0,.03125})
    for(unsigned axis=0;axis<3;++axis) for(double direction : {1.0,-1.0}) {
      fill_field(mu,.03125);fill_field(effective,.03125+turbulent);
      fill_field(thermal,thermal_gamma);
      // Fixture Sc=.5 and Sc_t=2; physical species diffusion stays unchanged.
      const double species_gamma=.0625+.5*turbulent;
      fill_field(physical,species_gamma);
      for(int z=-2;z<cells.z+2;++z)for(int y=-2;y<cells.y+2;++y)
        for(int x=-2;x<cells.x+2;++x){
          const int coordinate=axis==0?x:axis==1?y:z;
          const bool downstream=direction>0 ? coordinate>=3 : coordinate<3;
          first.view.unchecked({x,y,z},0)=downstream ? 5.6e-14 : 0;
          second.view.unchecked({x,y,z},0)=downstream ? .2 : .1;
        }
      for(unsigned a=0;a<3;++a){const auto f=faces[a];
        for(int z=0;z<f.extents.z;++z)for(int y=0;y<f.extents.y;++y)
          for(int x=0;x<f.extents.x;++x)f.unchecked({x,y,z})=a==axis?direction:0;
      }
      status=detail::prepare_mixture_diffusivity(fixture.equations.species(),
          as_const(thermal.view),as_const(mu.view),as_const(effective.view),minimum.view);
      MixtureTransportFaces mixture;
      if(status)status=prepare_cartesian_mixture_transport(fixture.equations.kernels(),
          {composition.data(),composition.size()},{},as_const(minimum.view),as_const(flow),
          work,901,mixture);
      if(status)status=form_cartesian_mixture_transport_flux(fixture.equations.kernels(),
          mixture,as_const(physical.view),as_const(flow),as_const(first.view),output);
      const double outward=direction*output[axis].unchecked({3,3,3});
      const double expected=-.125*std::max(0.0,species_gamma-thermal_gamma)*5.6e-14;
      passed &= expect(bool(status) && outward<=0 &&
          std::abs(outward-expected)<1e-28,
          "shared stabilization cannot remove an absent species with weaker physical diffusion");
    }
  return passed;
}

bool test_mixture_equation_faces() {
  ProductionFixture fixture;
  if (!expect(make_production_fixture(8,fixture,true),
              "three-species production equations compile")) return false;
  const auto cells = fixture.patch.cells;
  auto rho = make_field(kDensity,cells,1U,2U,501U);
  auto first = make_field(kSpecies,cells,1U,2U,502U);
  auto second = make_field(kPassive,cells,1U,2U,503U);
  auto gamma = make_field(25U,cells,1U,2U,504U);
  auto diagonal = make_field(30U,cells,1U,0U,505U);
  auto rhs = make_field(31U,cells,1U,0U,506U);
  auto residual = make_field(32U,cells,1U,0U,507U);
  fill_field(rho,1.0); fill_field(gamma,1e-6);
  for (int z=-2; z<cells.z+2; ++z) for (int y=-2; y<cells.y+2; ++y)
    for (int x=-2; x<cells.x+2; ++x) {
      first.view.unchecked({x,y,z},0U)=x<=1 ? 0.0 : x==2 ? .2 : .5;
      second.view.unchecked({x,y,z},0U)=x<=1 ? .6 : x==2 ? .7 : .5;
    }
  FinalFluxFixture owner;
  ConstFaceFluxView flux;
  if (!expect(make_linear_final_flux(fixture.equations.kernels(),cells,
                                     owner,flux,false,true),
              "unit-throughput final flux publishes")) return false;
  FaceFluxStorage storage;
  FaceFluxView workspace;
  auto status=FaceFluxStorage::allocate_workspace(cells,1U,storage);
  if (status) status=storage.workspace_view(0U,801U,workspace);
  const std::array<ConstFieldView,2> composition{as_const(first.view),as_const(second.view)};
  MixtureTransportFaces mixture;
  if (status) status=prepare_cartesian_mixture_transport(fixture.equations.kernels(),
      {composition.data(),composition.size()},{},as_const(gamma.view),flux,
      workspace,901U,mixture);
  if (!expect(bool(status),"common mixture face coefficients prepare")) return false;
  auto dependent=make_field(82U,cells,1U,2U,503U);
  for (int z=-2;z<cells.z+2;++z) for (int y=-2;y<cells.y+2;++y)
    for (int x=-2;x<cells.x+2;++x)
      dependent.view.unchecked({x,y,z},0U)=
          x<=1 ? .4 : x==2 ? .1 : 0.0;
  const std::array<ConstFieldView,3> full_composition{
      composition[0],composition[1],as_const(dependent.view)};
  std::vector<double> closure(std::size_t(cells.x)*cells.y*cells.z,0.0);
  const std::array<double,3> oracle{{.2,.1,-.3}};
  for (unsigned s=0;s<3;++s) {
    const auto q=full_composition[s];
    KernelInvocation call{{&q,1},{&residual.view,1},{{0,0,0},cells},
                          0,0,1,flux.revision};
    status=cartesian_mixture_transport(fixture.equations.kernels(),mixture,
        as_const(gamma.view),flux,call);
    // Complete REFERENCE VLS selects total conductance .5. The three flux
    // differences at x=2 are (.2,.1,-.3), independent of their split into
    // physical and artificial diffusion. At x=1 the upstream plateau holds.
    if (!expect(bool(status) &&
          close(residual.view.unchecked({2,3,3},0)/512.0,oracle[s]) &&
          std::abs(residual.view.unchecked({1,3,3},0))<1e-12,
          "complete mixture flux preserves upwind increments and zero species")) return false;
    for(int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x)
      closure[std::size_t((z*cells.y+y)*cells.x+x)]+=residual.view.unchecked({x,y,z},0);
  }
  for(double value:closure)
    if(!expect(std::abs(value)<2e-12,"complete species transport closes on every cell")) return false;
  const std::array<PrimitiveHistory,2> history{{
      {composition[0],composition[0],composition[0]},
      {composition[1],composition[1],composition[1]}}};
  EquationStateView state;
  state.density={as_const(rho.view),as_const(rho.view),as_const(rho.view)};
  state.independent_species={history.data(),history.size()};
  const std::array<ConstFieldView,2> diffusion{as_const(gamma.view),as_const(gamma.view)};
  EquationMaterialView material;
  material.scalar_mass_diffusivity={diffusion.data(),diffusion.size()};
  material.enthalpy_diffusivity=as_const(gamma.view);
  EquationAssemblyContext context;
  context.dt=.1; context.bdf={10.0,-10.0,0.0,1U}; context.time=701U;
  context.geometry=fixture.geometry.topology_revision();
  context.boundary=fixture.boundary.revision();
  context.thermo=fixture.thermodynamics.fingerprint();
  context.transport=fixture.transport.fingerprint();
  context.contribution_stage=1U;
  context.scope=EquationAssemblyScope::final_conservative;
  context.mass_flux=flux; context.face_flux=flux.revision;
  context.face_flux_authority=flux.certificate.authority();
  context.face_flux_storage=flux.certificate.storage();
  context.face_flux_revision_domain=flux.certificate.revision_domain();
  context.mixture_transport=&mixture;
  auto ax=make_face_field(CartesianAxis::x,cells,9401U);
  auto ay=make_face_field(CartesianAxis::y,cells,9402U);
  auto az=make_face_field(CartesianAxis::z,cells,9403U);
  bool passed=true;
  for (std::size_t s=0; s<2; ++s) {
    EquationAssemblyCertificate certificate;
    status=assemble_species(fixture.equations.species(),s,state,material,{},context,
        {diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view},certificate);
    const double r=residual.view.unchecked({2,3,3},0U);
    const double D=ax.view.unchecked({3,3,3});
    std::cout << "mixture_equation species=" << s << " residual=" << r
              << " conductance=" << D << '\n';
    passed &= expect(bool(status) && close(r,s==0 ? .2 : .1) && close(D,.5),
                     "production residual and matrix use the common REFERENCE face conductance");
  }
  auto h=make_field(kEnthalpy,cells,1U,2U,520U);
  auto correction_context=context;
  correction_context.time=mixture.linearization;
  passed &= check_scalar_correction_rows(fixture,composition[0],state,material,{},
      correction_context,{diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view});
  auto temperature=make_field(kTemperature,cells,1U,2U,521U);
  auto pressure=make_field(kPressure,cells,1U,2U,522U);
  auto velocity=make_field(kVelocity,cells,3U,2U,523U);
  auto gradient=make_field(kVelocityGradient,cells,9U,2U,524U);
  fill_field(temperature,300.0); fill_field(pressure,0.0);
  fill_field(velocity,0.0); fill_field(gradient,0.0);
  for (int z=-2; z<cells.z+2; ++z) for (int y=-2; y<cells.y+2; ++y)
    for (int x=-2; x<cells.x+2; ++x) {
      const Int3 c{x,y,z};
      h.view.unchecked(c,0U)=100.0+10.0*first.view.unchecked(c,0U)+
          20.0*second.view.unchecked(c,0U);
    }
  state.enthalpy={as_const(h.view),as_const(h.view),as_const(h.view)};
  state.temperature={as_const(temperature.view),as_const(temperature.view),as_const(temperature.view)};
  state.pressure_perturbation={as_const(pressure.view),as_const(pressure.view),as_const(pressure.view)};
  state.velocity={as_const(velocity.view),as_const(velocity.view),as_const(velocity.view)};
  state.pressure_reference=state.accepted_pressure_reference=state.previous_pressure_reference=101325.0;
  material.thermal_conductivity=material.molecular_viscosity=
      material.effective_viscosity=material.enthalpy_diffusivity=as_const(gamma.view);
  EquationAssemblyCertificate certificate;
  status=assemble_enthalpy(fixture.equations.enthalpy(),state,material,
      as_const(gradient.view),{},context,
      {diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view},certificate);
  const double h_residual=residual.view.unchecked({2,3,3},0U);
  std::cout << "mixture_equation enthalpy_residual=" << h_residual
            << " status=" << unsigned(status.code) << '/' << status.detail << '\n';
  passed &= expect(bool(status) && close(h_residual,4.0) &&
                       close(ax.view.unchecked({3,3,3}),.5),
                   "enthalpy references follow the same species face matrix and residual");
  auto old_pressure=make_field(kPressure,cells,1U,2U,525U);
  fill_field(old_pressure,0.0); fill_field(pressure,10.0);
  state.pressure_perturbation.accepted=state.pressure_perturbation.previous=as_const(old_pressure.view);
  status=detail::StatisticalEnthalpy::assemble(fixture.equations.enthalpy(),state,
      as_const(gamma.view),{},context,
      {diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view},certificate);
  std::cout << "statistical enthalpy transport residual=" << residual.view.unchecked({2,3,3},0) << '\n';
  passed &= expect(bool(status) && close(residual.view.unchecked({2,3,3},0),4.0),
      "statistical enthalpy matrix receives pressure work through its registered sources");
  passed &= check_scalar_correction_rows(fixture,as_const(h.view),state,material,{},
      correction_context,{diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view},false,false,true);
  // A sequential solve of A can transiently produce A+B>1 while B still
  // owns its previous iterate. Its B equation retains the frozen common
  // coefficients, and therefore matches the valid pre-update B equation.
  const EquationSystemView system{diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view};
  status=assemble_species(fixture.equations.species(),1,state,material,{},context,system,certificate);
  if(!expect(bool(status),"statistical component reference assembles")) return false;
  const auto reference_diagonal=diagonal.bytes,reference_rhs=rhs.bytes,reference_residual=residual.bytes;
  const auto reference_x=ax.bytes,reference_y=ay.bytes,reference_z=az.bytes;
  const auto public_certificate=certificate;
  // All histories use the same composition admitted by the thermodynamic
  // authority, including a trace smaller than one ulp of the major species.
  {
    auto major=make_field(kSpecies,cells,1U,2U,590U);
    auto trace=make_field(second.view.field,cells,1U,2U,591U);
    fill_field(major,1.0);
    fill_field(trace,1.0e-17);
    const std::array<double,2U> fractions{1.0,1.0e-17};
    ThermoState thermo;
    passed &= expect(bool(fixture.thermodynamics.evaluate(101325.0,300000.0,
        {fractions.data(),fractions.size()},{},thermo)),
        "EOS admits the represented major-plus-trace assembly fixture");
    const std::array<PrimitiveHistory,2U> trace_history{{
        {as_const(major.view),as_const(major.view),as_const(major.view)},
        {as_const(trace.view),as_const(trace.view),as_const(trace.view)}}};
    state.independent_species={trace_history.data(),trace_history.size()};
    for (std::size_t s=0; s<2; ++s) {
      status=assemble_species(fixture.equations.species(),s,state,material,{},context,system,certificate);
      passed &= expect(status && certificate.valid(),
          "public species assembly admits EOS-representable trace in every history");
    }
    passed &= expect(major.view.unchecked({0,0,0},0U)==1.0 &&
        trace.view.unchecked({0,0,0},0U)==1.0e-17,
        "species assembly does not normalize the input composition");
    state.independent_species={history.data(),history.size()};
  }
  auto trial_a=make_field(kSpecies,cells,1,2,527U);
  std::copy(first.bytes.begin(),first.bytes.end(),trial_a.bytes.begin());
  trial_a.view.unchecked({3,3,3},0)=.625;
  auto trial_history=history;
  trial_history[0].trial=as_const(trial_a.view);
  state.independent_species={trial_history.data(),trial_history.size()};
  const auto reset=[&] {
    fill_field(diagonal,91);fill_field(rhs,91);fill_field(residual,91);
    for(auto* face:{&ax,&ay,&az})std::fill(face->bytes.begin(),face->bytes.end(),91.);
    certificate=public_certificate;
  };
  const auto untouched=[&] {
    return output_is(diagonal,rhs,residual,91) && faces_are(ax,ay,az,91) &&
        same_certificate(certificate,public_certificate);
  };
  reset();
  status=assemble_species(fixture.equations.species(),1,state,material,{},context,system,certificate);
  passed &= expect(!status && untouched(),"public assembly keeps complete composition audit atomic");
  reset();
  status=detail::StatisticalSpecies::assemble(fixture.equations.species(),1,state,material,{},context,system,certificate);
  passed &= expect(bool(status) && diagonal.bytes==reference_diagonal && rhs.bytes==reference_rhs &&
      residual.bytes==reference_residual && ax.bytes==reference_x && ay.bytes==reference_y && az.bytes==reference_z,
      "statistical component assembly permits a finite off-simplex trial with fixed common faces");
  passed &= expect(bool(status) && certificate.plan!=public_certificate.plan,
      "statistical component certificate has its own plan identity");
  detail::StatisticalFaceBasis basis;
  context.statistical_face_basis=&basis;
  for(unsigned repeat=0;repeat<2;++repeat) {
    status=detail::StatisticalSpecies::assemble(fixture.equations.species(),1,state,material,{},context,system,certificate);
    passed &= expect(bool(status) && diagonal.bytes==reference_diagonal && rhs.bytes==reference_rhs &&
        residual.bytes==reference_residual && ax.bytes==reference_x && ay.bytes==reference_y && az.bytes==reference_z,
        "prepared common faces preserve the complete matrix, RHS and true residual exactly");
  }
  passed &= expect(basis.reuses==1,"second component assembly reuses frozen conductances");
  // Changing the material revision rebuilds even when values remain equal.
  auto revised_diffusion=diffusion;
  for(auto& value:revised_diffusion)++value.revision;
  material.scalar_mass_diffusivity={revised_diffusion.data(),revised_diffusion.size()};
  ++material.enthalpy_diffusivity.revision;
  status=detail::StatisticalSpecies::assemble(fixture.equations.species(),1,state,material,{},context,system,certificate);
  passed &= expect(bool(status) && basis.reuses==1 && ax.bytes==reference_x,
      "material revision invalidates prepared conductances");
  material.scalar_mass_diffusivity={diffusion.data(),diffusion.size()};
  material.enthalpy_diffusivity=as_const(gamma.view);
  context.statistical_face_basis=nullptr;

  std::cout<<"statistical species partial_sum="<<1.125<<" status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';
  passed &= check_scalar_correction_rows(fixture,as_const(trial_a.view),state,material,{},
      correction_context,system,false,false,false,true);
  const auto rejected=[&](const char* description) {
    reset();
    const auto result=detail::StatisticalSpecies::assemble(fixture.equations.species(),1,state,material,{},context,system,certificate);
    return expect(!result && untouched(),description);
  };
  trial_history[0].accepted=as_const(trial_a.view);
  passed &= rejected("statistical species rejects off-simplex accepted history atomically");
  trial_history=history;trial_history[0].previous=as_const(trial_a.view);
  passed &= rejected("statistical species rejects off-simplex previous history atomically");
  trial_history=history;trial_history[0].trial=as_const(trial_a.view);
  trial_a.view.unchecked({3,3,3},0)=std::numeric_limits<double>::quiet_NaN();
  passed &= rejected("statistical species rejects nonfinite component iterate atomically");
  trial_a.view.unchecked({3,3,3},0)=.625;
  auto stale_context=context;
  ++context.face_flux;
  passed &= rejected("statistical component assembly checks committed flux authority");
  context=stale_context;
  context.reaction_endpoint=history[1].accepted;
  passed &= rejected("statistical spatial assembly keeps the chemical endpoint separate");
  context.reaction_endpoint={};
  return passed;
}

bool test_inlet_scalar_material() {
  constexpr int n = 4;
  bool passed = true;
  for (bool physical : {false,true}) for (int inlet = 0; inlet < 6; ++inlet) {
    ProductionFixture fixture;
    if (!expect(make_production_fixture(n,fixture,false,inlet,physical),
                "physical inlet scalar fixture compiles")) return false;
    const auto cells = fixture.patch.cells;
    const auto axis = static_cast<std::size_t>(inlet / 2);
    const auto normal = [axis](Int3 c) { return axis==0 ? c.x : axis==1 ? c.y : c.z; };
    auto rho = make_field(kDensity,cells,1U,2U,601U);
    auto q = make_field(kSpecies,cells,1U,2U,602U);
    auto passive = make_field(kPassive,cells,1U,2U,603U);
    auto gamma = make_field(25U,cells,1U,2U,604U);
    fill_field(rho,1.0); fill_field(gamma,1.0);
    for (int z=-2; z<n+2; ++z) for (int y=-2; y<n+2; ++y)
      for (int x=-2; x<n+2; ++x) {
        const Int3 c{x,y,z};
        q.view.unchecked(c,0U)=passive.view.unchecked(c,0U)=
            0.25+0.1*(normal(c)+0.5)/n;
        if (normal(c)==(inlet%2 ? n : -1)) gamma.view.unchecked(c,0U)=9.0;
      }
    FinalFluxFixture owner;
    ConstFaceFluxView flux;
    if (!make_linear_final_flux(fixture.equations.kernels(),cells,owner,flux,true)) return false;
    const PrimitiveHistory species{as_const(q.view),as_const(q.view),as_const(q.view)};
    const PrimitiveHistory tracer{as_const(passive.view),as_const(passive.view),as_const(passive.view)};
    EquationStateView state;
    state.density={as_const(rho.view),as_const(rho.view),as_const(rho.view)};
    state.independent_species={&species,1U}; state.passive_scalars={&tracer,1U};
    const std::array<ConstFieldView,2> diffusivities{as_const(gamma.view),as_const(gamma.view)};
    EquationMaterialView material;
    material.scalar_mass_diffusivity={diffusivities.data(),diffusivities.size()};
    EquationAssemblyContext context;
    context.dt=.1; context.bdf={10.0,-10.0,0.0,1U}; context.time=701U;
    context.geometry=fixture.geometry.topology_revision(); context.boundary=fixture.boundary.revision();
    context.thermo=fixture.thermodynamics.fingerprint(); context.transport=fixture.transport.fingerprint();
    context.contribution_stage=1U; context.scope=EquationAssemblyScope::final_conservative;
    context.mass_flux=flux; context.face_flux=flux.revision;
    context.face_flux_authority=flux.certificate.authority();
    context.face_flux_storage=flux.certificate.storage();
    context.face_flux_revision_domain=flux.certificate.revision_domain();
    auto diagonal=make_field(30U,cells,1U,0U,610U);
    auto rhs=make_field(31U,cells,1U,0U,611U);
    auto residual=make_field(32U,cells,1U,0U,612U);
    auto ax=make_face_field(CartesianAxis::x,cells,613U);
    auto ay=make_face_field(CartesianAxis::y,cells,614U);
    auto az=make_face_field(CartesianAxis::z,cells,615U);
    EquationSystemView system{diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view};
    Int3 cell{1,1,1},face{1,1,1};
    (axis==0 ? cell.x : axis==1 ? cell.y : cell.z)=inlet%2 ? n-1 : 0;
    (axis==0 ? face.x : axis==1 ? face.y : face.z)=inlet%2 ? n : 0;
    const std::array<FaceFieldView,3> coefficients{ax.view,ay.view,az.view};
    for (int role=0; role<2; ++role) {
      EquationAssemblyCertificate certificate;
      const auto status=role==0
          ? assemble_species(fixture.equations.species(),0U,state,material,{},context,system,certificate)
          : assemble_scalar(fixture.equations.scalars(),0U,state,material,{},context,system,certificate);
      // Face-state storage uses Gamma=9 directly. Exterior-cell storage
      // uses two equal half-cell resistances, giving Gamma=1.8.
      const double face_gamma=physical ? 9.0 : 1.8;
      const double expected_diagonal=10.0/(n*n*n)+(5.0+face_gamma)/n;
      const double expected_residual=(inlet%2 ? -1.0 : 1.0)*(face_gamma-1.0)*0.1/(n*n);
      std::cout << "inlet_scalar face=" << inlet << " role=" << role << " physical=" << physical
                << " coefficient=" << coefficients[axis].unchecked(face)
                << " residual=" << residual.view.unchecked(cell,0U) << '\n';
      passed &= expect(bool(status) && close(coefficients[axis].unchecked(face),face_gamma/n) &&
          close(diagonal.view.unchecked(cell,0U),expected_diagonal) &&
          close(residual.view.unchecked(cell,0U),expected_residual) &&
          close(rhs.view.unchecked(cell,0U),expected_diagonal*q.view.unchecked(cell,0U)-expected_residual),
          "scalar matrix and complete residual match the configured inlet material storage");
      passed &= check_scalar_correction_rows(fixture,role==0 ? as_const(q.view) : as_const(passive.view),
          state,material,{},context,system,role!=0);
    }
  }
  return passed;
}

bool test_ibm_species_matrix(bool passive = false, bool periodic = false) {
  constexpr int n=16;
  ProductionFixture fixture;
  if (!expect(make_production_fixture(n,fixture,false,-1,true,passive ? 0 : 1./(n*n*n)),
              "IBM species equation fixture compiles")) return false;
  auto triangles=test::force_cube();
  const auto move=[&](Real3 p) {return Real3{(periodic ? 0.72 : 0.5)+0.4*p.x,0.5+0.4*p.y,0.5+0.4*p.z};};
  for(auto& t:triangles) t={move(t.a),move(t.b),move(t.c)};
  StlScanPlan scan; ImmersedSurfacePlan surface; EBTopology topology;
  BoundaryStencilPlan stencils; ImmersedPlanLimits limits;
  ImmersedDomainBoundaryPolicy policy;
  policy.allow_periodic_images.fill(periodic);
  limits.stencil.policy=IbmReconstructionPolicy::adaptive_order;
  IbmEquationInterfacePlan ibm;
  auto status=StlScanCompiler::compile_triangles(fixture.geometry,fixture.patch,
      {triangles.data(),triangles.size()},CartesianAxis::y,test::kForceScanBudget,scan);
  if(status) status=ImmersedSurfaceCompiler::compile(scan,surface);
  if(status) status=EBTopologyCompiler::compile(MPI_COMM_SELF,fixture.geometry,fixture.patch,
      scan,surface,ImmersedFluidSide::outside,limits,topology,policy);
  if(status) status=BoundaryStencilCompiler::compile(MPI_COMM_SELF,fixture.geometry,
      fixture.patch,surface,topology,policy,limits,stencils);
  if(status) status=IbmEquationInterfacePlan::compile(fixture.equations.kernels(),topology,
      stencils,topology.interface_metric(),ibm);
  if(!expect(bool(status),"IBM scalar topology compiles")) return false;
  const auto cells=fixture.patch.cells;
  auto rho=make_field(kDensity,cells,1U,2U,501U);
  auto q=make_field(passive ? kPassive : kSpecies,cells,1U,2U,502U);
  auto gamma=make_field(25U,cells,1U,2U,503U);
  auto diagonal=make_field(30U,cells,1U,0U,504U);
  auto rhs=make_field(31U,cells,1U,0U,505U);
  auto residual=make_field(32U,cells,1U,0U,506U);
  fill_field(rho,1.0); fill_field(q,0.25); fill_field(gamma,1.0);
  for(int z=0;z<n;++z) for(int y=0;y<n;++y) for(int x=0;x<n;++x)
    if(topology.region().data[x+n*(y+n*z)]==0U) q.view.unchecked({x,y,z},0U)=0.75;
  FaceFluxStorage flux_storage; FaceFluxView flux;
  status=FaceFluxStorage::allocate_workspace(cells,1U,flux_storage);
  if(status) status=flux_storage.workspace_view(0U,11U,flux);
  if(!expect(bool(status),"IBM scalar zero flux allocates")) return false;
  for(auto face : {flux.x,flux.y,flux.z})
    for(int z=0;z<face.extents.z;++z) for(int y=0;y<face.extents.y;++y)
      for(int x=0;x<face.extents.x;++x) face.unchecked({x,y,z})=0.0;
  const PrimitiveHistory q_history{as_const(q.view),as_const(q.view),as_const(q.view)};
  const auto d=as_const(gamma.view);
  EquationStateView state;
  state.density={as_const(rho.view),as_const(rho.view),as_const(rho.view)};
  state.independent_species={&q_history,1U};
  EquationMaterialView material; material.scalar_mass_diffusivity={&d,1U};
  EquationAssemblyContext context;
  context.dt=0.1; context.bdf={10.0,-10.0,0.0,1U}; context.time=701U;
  context.geometry=fixture.geometry.topology_revision(); context.boundary=fixture.boundary.revision();
  context.thermo=fixture.thermodynamics.fingerprint(); context.transport=fixture.transport.fingerprint();
  context.contribution_stage=1U; context.scope=EquationAssemblyScope::momentum_predictor;
  context.mass_flux=as_const(flux); context.face_flux=flux.revision; context.provisional_mass_flux=true;
  context.immersed_interface=&ibm;
  auto ax=make_face_field(CartesianAxis::x,cells,601U);
  auto ay=make_face_field(CartesianAxis::y,cells,602U);
  auto az=make_face_field(CartesianAxis::z,cells,603U);
  EquationSystemView system{diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view};
  FinalFluxFixture final_flux;
  ConstFaceFluxView committed;
  const std::array diffusivities{d,d};
  EquationAssemblyCertificate certificate;
  if(passive) {
    if(!make_linear_final_flux(fixture.equations.kernels(),cells,final_flux,committed,true)) return false;
    context.scope=EquationAssemblyScope::final_conservative;
    context.mass_flux=committed; context.face_flux=committed.revision;
    context.face_flux_authority=committed.certificate.authority();
    context.face_flux_storage=committed.certificate.storage();
    context.face_flux_revision_domain=committed.certificate.revision_domain();
    context.provisional_mass_flux=false;
    state.passive_scalars={&q_history,1U};
    material.scalar_mass_diffusivity={diffusivities.data(),diffusivities.size()};
    status=assemble_scalar(fixture.equations.scalars(),0U,state,material,{},context,system,certificate);
  } else {
    status=detail::assemble_species_coupling_rows(fixture.equations.species(),0U,state,material,context,system);
  }
  if(!status) std::cerr<<"IBM rows status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';
  if(!expect(bool(status),"IBM production species rows assemble")) return false;
  const auto region=topology.region();
  const auto fluid=[&](Int3 c) {return c.x<0 || c.y<0 || c.z<0 || c.x>=n || c.y>=n || c.z>=n ||
      region.data[c.x+n*(c.y+n*c.z)]!=0U;};
  double error=0.0, solid_residual=0.0, cut_coefficient=0.0;
  for(int z=0;z<n;++z) for(int y=0;y<n;++y) for(int x=0;x<n;++x) {
    const Int3 c{x,y,z}; unsigned neighbors=0U;
    for(int a=0;a<3;++a) for(int sign : {-1,1}) {
      Int3 neighbor=c; (a==0 ? neighbor.x : a==1 ? neighbor.y : neighbor.z)+=sign;
      neighbors+=fluid(neighbor);
      if(!fluid(c) || !fluid(neighbor)) {
        const auto face=a==0 ? ax.view : a==1 ? ay.view : az.view;
        cut_coefficient=std::max(cut_coefficient,std::abs(face.unchecked(sign>0 ? neighbor : c)));
      }
    }
    const double expected=fluid(c) ? (10.0+neighbors*n*n)*(passive ? 1.0/(n*n*n) : 1.0) : 1.0;
    error=std::max(error,std::abs(diagonal.view.unchecked(c,0U)-expected));
    solid_residual=std::max(solid_residual,std::abs(residual.view.unchecked(c,0U)));
  }
  std::cerr<<"IBM "<<(passive ? "passive" : "species")<<" matrix diagonal_error="<<error<<" cut_coefficient="<<cut_coefficient
           <<" solid_residual="<<solid_residual<<'\n';
  bool passed=expect(error<1e-11 && cut_coefficient==0.0 && solid_residual==0.0,
      "species matrix excludes solid edges and constrains solid corrections");
  if (passive) {
    ScalarMidpointView midpoint{as_const(q.view)};
    context.scalar_midpoint=&midpoint;
    passed &= check_scalar_correction_rows(fixture,as_const(q.view),state,material,
        {},context,system,true);
    return passed;
  }
  {
    auto source=make_field(10,cells,1,0,850);
    auto sink=make_field(11,cells,1,0,851);
    esf::detail::IemSource terms;
    if (esf::detail::iem_source(1./(n*n*n),1.8e-5,1e-4,2,1,.5,terms)!=
        portable::Status::success) return false;
    fill_field(source,terms.explicit_source_density);
    fill_field(sink,terms.implicit_sink_density);
    EquationContributionView mixing;
    mixing.explicit_source_density=as_const(source.view);
    mixing.implicit_sink_density=as_const(sink.view);
    mixing.has_implicit_sink=true; mixing.conserved_quantity=kSpecies;
    mixing.units.si_exponents={1,-3,-1,0,0,0,0}; mixing.stage=9;
    mixing.explicit_source_field=10; mixing.implicit_sink_field=11;
    auto mixed_context=context; mixed_context.contribution_stage=9;
    status=detail::assemble_species_coupling_rows(fixture.equations.species(),0,
        state,material,mixed_context,system,{&mixing,1});
    double mixed_error{},mixed_residual{};
    for(int z=0;z<n;++z) for(int y=0;y<n;++y) for(int x=0;x<n;++x) {
      const Int3 cell{x,y,z}; unsigned neighbors{};
      for (unsigned a=0;a<3;++a) for(int sign:{-1,1}) {
        auto next=cell; (a==0 ? next.x : a==1 ? next.y : next.z)+=sign;
        neighbors+=fluid(next);
      }
      const double expected=fluid(cell) ? 10+neighbors*n*n+terms.implicit_sink_density : 1;
      mixed_error=std::max(mixed_error,std::abs(diagonal.view.unchecked(cell,0)-expected));
      mixed_residual=std::max(mixed_residual,std::abs(residual.view.unchecked(cell,0)-
          (fluid(cell) ? -.25*terms.implicit_sink_density : 0)));
    }
    std::cout<<"IBM IEM diagonal_error="<<mixed_error<<" residual_error="<<mixed_residual<<'\n';
    passed &= expect(bool(status) && mixed_error<1e-11 && mixed_residual<1e-12,
        "registered IEM acts on fluid rows and preserves constrained solid rows");
    if (!make_linear_final_flux(fixture.equations.kernels(),cells,final_flux,committed,true)) return false;
    mixed_context.scope=EquationAssemblyScope::final_conservative;
    mixed_context.mass_flux=committed; mixed_context.face_flux=committed.revision;
    mixed_context.face_flux_authority=committed.certificate.authority();
    mixed_context.face_flux_storage=committed.certificate.storage();
    mixed_context.face_flux_revision_domain=committed.certificate.revision_domain();
    mixed_context.provisional_mass_flux=false;
    passed &= check_scalar_correction_rows(fixture,as_const(q.view),state,material,
        {&mixing,1},mixed_context,system);
  }
  Int3 face{-1,n/2,n/2};
  for (int x=2; x<n-2; ++x)
    if (fluid({x-1,n/2,n/2}) && fluid({x,n/2,n/2}) && !fluid({x+1,n/2,n/2})) {
      face.x=x; break;
    }
  if (!expect(face.x>=2,"IBM downstream stencil meets the immersed wall")) return false;
  fill_field(gamma,1e-6);
  for (int z=-2; z<n+2; ++z) for (int y=-2; y<n+2; ++y)
    for (int x=-2; x<n+2; ++x)
      q.view.unchecked({x,y,z},0U)=x<face.x-1 ? 0.0 : x==face.x-1 ? .2 : .5;
  for (int z=0; z<n; ++z) for (int y=0; y<n; ++y)
    for (int x=0; x<=n; ++x) flux.x.unchecked({x,y,z})=1.0;
  status=ibm.constrain_interface_flux(flux);
  FaceFluxStorage limiter_storage;
  FaceFluxView limiter_workspace;
  if (status) status=FaceFluxStorage::allocate_workspace(cells,1U,limiter_storage);
  if (status) status=limiter_storage.workspace_view(0U,901U,limiter_workspace);
  MixtureTransportFaces mixture;
  const auto fraction=as_const(q.view);
  if (status) status=prepare_cartesian_mixture_transport(fixture.equations.kernels(),
      {&fraction,1U},{},as_const(gamma.view),as_const(flux),limiter_workspace,
      902U,mixture,&ibm);
  const double total=status ? mixture.x.unchecked(face)+1e-6/n : -1.0;
  std::cout << "IBM mixture face=" << face.x << " conductance=" << total << '\n';
  passed &= expect(bool(status) && close(total,1e-6/n),
      "IBM mixture limiter uses the available fluid upstream gradient");
  const auto check_wall_transport=[&]() {
    KernelInvocation call{{&fraction,1},{&residual.view,1},{{0,0,0},cells},
                          0,0,1,flux.revision};
    auto transport_status=cartesian_mixture_transport(fixture.equations.kernels(),
        mixture,as_const(gamma.view),as_const(flux),call);
    if (transport_status) transport_status=detail::IbmScalarTransport::transport(
        ibm,0U,fraction,as_const(gamma.view),as_const(flux),
        {{0,0,0},cells},residual.view,mixture);
    // The open left face carries .35 - gamma*area/distance*(.5-.2).
    // The stationary solid right face is impermeable. Divide by 1/volume
    // to compare the physical net flux against this independent oracle.
    const double net=residual.view.unchecked(face,0U)/(n*n*n);
    std::cout << "IBM mixture wall net_flux=" << std::setprecision(17) << net << '\n';
    return expect(bool(transport_status) && close(net,-.35+.3e-6/n),
        "IBM complete transport retains fluid diffusion and wall impermeability");
  };
  passed &= check_wall_transport();
  // The REFERENCE selector's downstream fallback reaches the solid placeholder.
  // Fluid-face coefficients must stay invariant under changes to that state.
  for(int z=0;z<n;++z) for(int y=0;y<n;++y) for(int x=0;x<n;++x)
    if(!fluid({x,y,z})) q.view.unchecked({x,y,z},0U)=.8;
  status=prepare_cartesian_mixture_transport(fixture.equations.kernels(),
      {&fraction,1U},{},as_const(gamma.view),as_const(flux),limiter_workspace,
      903U,mixture,&ibm);
  passed &= expect(bool(status) && close(mixture.x.unchecked(face)+1e-6/n,total),
      "solid placeholder changes preserve the fluid face coefficient");
  passed &= check_wall_transport();
  for(int z=0;z<n;++z) for(int y=0;y<n;++y) for(int x=0;x<=n;++x)
    flux.x.unchecked({x,y,z})=-1.0;
  if(status) status=ibm.constrain_interface_flux(flux);
  if(status) status=prepare_cartesian_mixture_transport(fixture.equations.kernels(),
      {&fraction,1U},{},as_const(gamma.view),as_const(flux),limiter_workspace,
      904U,mixture,&ibm);
  passed &= expect(bool(status) && close(mixture.x.unchecked(face)+1e-6/n,.5),
      "an unavailable upstream donor gives a bounded one-sided fluid face");
  if (periodic) {
    // The x=14 solid is the second upstream donor of the periodic x=0 face.
    // Supplying a smooth placeholder gradient still requires a one-sided face.
    for (int z=-2; z<n+2; ++z) for (int y=-2; y<n+2; ++y)
      for (int x=-2; x<n+2; ++x) q.view.unchecked({x,y,z},0U)=.2+.01*x;
    for (int z=0; z<n; ++z) for (int y=0; y<n; ++y)
      for (int x=0; x<=n; ++x) flux.x.unchecked({x,y,z})=1.0;
    status=ibm.constrain_interface_flux(flux);
    if (status) status=prepare_cartesian_mixture_transport(fixture.equations.kernels(),
        {&fraction,1U},{},as_const(gamma.view),as_const(flux),limiter_workspace,
        905U,mixture,&ibm);
    passed &= expect(bool(status) && close(mixture.x.unchecked({0,n/2,n/2})+1e-6/n,.5),
        "periodic species face excludes the wrapped solid donor");
  }
  // The pressure matrix uses the same fluid-side stencil authority as
  // transported scalars. A solid placeholder carries no physical pressure.
  auto pressure=make_field(kPressure,cells,1U,2U,620U);
  auto velocity=make_field(kVelocity,cells,3U,2U,621U);
  auto mass_source=make_field(40U,cells,1U,0U,622U);
  for(int z=0;z<n;++z) for(int y=0;y<n;++y) for(int x=0;x<n;++x)
    mass_source.view.unchecked({x,y,z},0)=topology.is_fluid_global({x,y,z})
        ? 1.25 : std::numeric_limits<double>::quiet_NaN();
  for (unsigned axis=0; axis<3; ++axis) for (double sign : {-1.0,1.0}) {
    fill_field(velocity,0.0);
    for (int z=-2;z<n+2;++z) for(int y=-2;y<n+2;++y)
      for(int x=-2;x<n+2;++x)
        velocity.view.unchecked({x,y,z},axis)=sign;
    std::vector<detail::ColdPressureRow> baseline;
    for (double placeholder : {-1000.0,1000.0}) {
      for (int z=-2;z<n+2;++z) for(int y=-2;y<n+2;++y)
        for(int x=-2;x<n+2;++x) {
          const Int3 c{x,y,z};
          pressure.view.unchecked(c,0)=topology.is_fluid_stencil(c)
              ? 10.0*(x+y+z) : placeholder;
        }
      std::vector<detail::ColdPressureRow> rows;
      detail::ColdGridReport report;
      const bool assembled=detail::assemble_midpoint_cold_grid(
          fixture.equations.kernels(),fixture.patch,fixture.geometry.global_cells(),
          &topology,fixture.boundary,as_const(rho.view),as_const(rho.view),
          as_const(velocity.view),as_const(velocity.view),as_const(pressure.view),
          100000.0,0.01,as_const(flux),rows,report);
      passed &= expect(assembled,"IBM pressure matrix prepares for both flow directions");
      if (!assembled) return false;
      std::vector<detail::ColdPressureRow> sourced;
      detail::ColdGridReport sourced_report;
      const bool source_ok=detail::assemble_midpoint_cold_grid(
          fixture.equations.kernels(),fixture.patch,fixture.geometry.global_cells(),
          &topology,fixture.boundary,as_const(rho.view),as_const(rho.view),
          as_const(velocity.view),as_const(velocity.view),as_const(pressure.view),
          100000.0,0.01,as_const(flux),sourced,sourced_report,nullptr,
          {as_const(mass_source.view),71U,623U},623U);
      passed &= expect(source_ok && sourced.size()==rows.size(),
          "IBM mass-source pressure grid prepares");
      if(!source_ok || sourced.size()!=rows.size())return false;
      for(int z=0;z<n;++z) for(int y=0;y<n;++y) for(int x=0;x<n;++x) {
        const auto i=std::size_t(x+n*(y+n*z));
        const bool active=topology.is_fluid_global({x,y,z});
        passed &= expect(sourced[i].diagonal==rows[i].diagonal &&
            sourced[i].neighbour==rows[i].neighbour &&
            close(sourced[i].rhs,rows[i].rhs+(active ? 1.25 : 0.)),
            "IBM mass exchange affects fluid RHS while solid rows retain their identity");
      }
      if (baseline.empty()) baseline=rows;
      else {
        bool same=true;
        for(std::size_t i=0;i<rows.size();++i)
          same &= rows[i].diagonal==baseline[i].diagonal &&
                  rows[i].rhs==baseline[i].rhs && rows[i].neighbour==baseline[i].neighbour;
        passed &= expect(same,"solid pressure changes preserve fluid pressure matrix and RHS");
      }
    }
  }
  const double volume=1.0/(n*n*n),old_h=100.0,inlet_h=200.0;
  double stochastic_flux_error=0;
  auto h=make_field(kEnthalpy,cells,1U,2U,610U);
  auto old_rho=make_field(kDensity,cells,1U,2U,611U);
  fill_field(h,old_h); fill_field(old_rho,1.0);
  state.enthalpy={as_const(h.view),as_const(h.view),as_const(h.view)};
  state.density={as_const(rho.view),as_const(old_rho.view),as_const(old_rho.view)};
  fill_field(q,.5);
  for(int z=0;z<n;++z) for(int y=0;y<n;++y) for(int x=0;x<n;++x)
    if(!fluid({x,y,z})) q.view.unchecked({x,y,z},0U)=.8;
  for (unsigned direction=0; direction<6; ++direction) {
    const auto links=topology.links();
    const auto found=std::find_if(links.data,links.data+links.size,
        [&](const auto& item) {return static_cast<unsigned>(item.direction)==direction;});
    if (!expect(found!=links.data+links.size,"IBM cube provides every source direction")) return false;
    const auto& link=*found;
    const auto axis=direction/2;
    const double outward=direction%2 ? 1.0 : -1.0;
    const Int3 source_cell=link.fluid_local_index;
    const double mass_flux=-outward*1e-5;
    const double source_Y=.25;
    Real3 source_velocity{};
    (axis==0 ? source_velocity.x : axis==1 ? source_velocity.y : source_velocity.z)=
        mass_flux>0 ? 1.0 : -1.0;
    const IbmInterfaceInletState inlet{link.global_link,mass_flux,source_velocity,inlet_h,{&source_Y,1U}};
    IbmEquationInterfacePlan source_ibm;
    status=IbmEquationInterfacePlan::compile(fixture.equations.kernels(),topology,
        stencils,topology.interface_metric(),{&inlet,1U},1U,source_ibm);
    for (auto field : {flux.x,flux.y,flux.z,ax.view,ay.view,az.view})
      for (int z=0; z<field.extents.z; ++z) for (int y=0; y<field.extents.y; ++y)
        for (int x=0; x<field.extents.x; ++x) field.unchecked({x,y,z})=0.0;
    if (status) status=source_ibm.constrain_interface_flux(flux);
    if (status) status=prepare_cartesian_mixture_transport(fixture.equations.kernels(),
        {&fraction,1U},{},as_const(gamma.view),as_const(flux),limiter_workspace,
        906U+direction,mixture,&source_ibm);
    KernelInvocation source_call{{&fraction,1},{&residual.view,1},{{0,0,0},cells},
                                 0,0,1,flux.revision};
    if (status) status=cartesian_mixture_transport(fixture.equations.kernels(),
        mixture,as_const(gamma.view),as_const(flux),source_call);
    if (status) status=detail::IbmScalarTransport::transport(source_ibm,0U,
        fraction,as_const(gamma.view),as_const(flux),{{0,0,0},cells},
        residual.view,mixture);
    passed &= expect(bool(status) &&
        close(residual.view.unchecked(source_cell,0U)*volume,outward*mass_flux*source_Y),
        "IBM complete species transport uses the prescribed inlet composition on all axes");
    for (const auto quantity : {detail::IbmScalarTransport::Quantity::enthalpy,
                               detail::IbmScalarTransport::Quantity::dependent_species})
      for (const double placeholder : {-1000.,1000.})
        for (const bool common : {false,true}) {
          const bool heat=quantity==detail::IbmScalarTransport::Quantity::enthalpy;
          fill_field(h,heat ? old_h : .5);
          for (int z=-2;z<n+2;++z) for (int y=-2;y<n+2;++y) for (int x=-2;x<n+2;++x)
            if (!topology.is_fluid_stencil({x,y,z})) h.view.unchecked({x,y,z},0)=placeholder;
          const auto scalar=as_const(h.view);
          source_call.reads={&scalar,1};
          if (status) status=common
              ? cartesian_mixture_transport(fixture.equations.kernels(),mixture,
                    as_const(gamma.view),as_const(flux),source_call)
              : cartesian_provisional_convection(fixture.equations.kernels(),ConvectionScheme::central2,
                    as_const(flux),source_call);
          if (status) status=common
              ? detail::IbmScalarTransport::transport(source_ibm,{quantity,0},scalar,
                    as_const(gamma.view),as_const(flux),{{0,0,0},cells},residual.view,mixture)
              : detail::IbmScalarTransport::convection(source_ibm,{quantity,0},ConvectionScheme::central2,
                    scalar,as_const(flux),{{0,0,0},cells},residual.view);
          const double prescribed=heat ? inlet_h : 1-source_Y;
          if (!status) std::cerr << "IBM stochastic transport common=" << common
              << " heat=" << heat << " status=" << unsigned(status.code) << '/' << status.detail << '\n';
          const double exact=outward*mass_flux*prescribed;
          if(common) {
            auto ledger_x=make_face_field(CartesianAxis::x,cells,9901);
            auto ledger_y=make_face_field(CartesianAxis::y,cells,9902);
            auto ledger_z=make_face_field(CartesianAxis::z,cells,9903);
            const std::array<FaceFieldView,3> ledger{ledger_x.view,ledger_y.view,ledger_z.view};
            if(status)status=form_cartesian_mixture_transport_flux(fixture.equations.kernels(),
                mixture,as_const(gamma.view),as_const(flux),scalar,ledger);
            if(status)status=detail::IbmScalarTransport::constrain_flux(source_ibm,{quantity,0},as_const(flux),ledger);
            const auto cut_links=topology.links();
            for(std::size_t k=0;k<cut_links.size && status;++k) {
              const auto& cut=cut_links.data[k];const unsigned a=static_cast<unsigned>(cut.direction)/2;
              Int3 face=cut.fluid_local_index;
              if(static_cast<unsigned>(cut.direction)%2)++(a==0 ? face.x : a==1 ? face.y : face.z);
              const double expected=cut.global_link==link.global_link ? mass_flux*prescribed : 0.;
              passed &= expect(close(ledger[a].unchecked(face),expected),
                  "statistical IBM face ledger preserves prescribed inlet flux and sealed cut faces");
            }
            long double face_sum{};
            for(unsigned a=0;a<3;++a) {
              Int3 upper=source_cell;++(a==0 ? upper.x : a==1 ? upper.y : upper.z);
              face_sum+=static_cast<long double>(ledger[a].unchecked(upper))-ledger[a].unchecked(source_cell);
            }
            passed &= expect(status && close(static_cast<double>(face_sum),exact),
                "statistical IBM face ledger has the independent fluid-control-volume balance");
          }

          const double error=std::abs(residual.view.unchecked(source_cell,0U)*volume-exact)/std::abs(exact);
          stochastic_flux_error=std::max(stochastic_flux_error,error);
          passed &= expect(bool(status) && error<1e-11,
              "IBM stochastic enthalpy and dependent species use their inlet states on all axes");
        }
    fill_field(h,old_h);
    fill_field(rho,1.0); fill_field(diagonal,volume/context.dt); fill_field(residual,0.0);
    const double F=outward*mass_flux;
    const double new_rho=1.0-context.dt*F/volume;
    rho.view.unchecked(source_cell,0U)=new_rho;
    diagonal.view.unchecked(source_cell,0U)=new_rho*volume/context.dt;
    residual.view.unchecked(source_cell,0U)=F*(inlet_h-old_h);
    context.mixture_transport=&mixture;
    context.immersed_interface=&source_ibm;
    std::vector<detail::ColdPressureRow> rows;
    if (status) status=detail::close_cold_enthalpy_rows(fixture.equations.kernels(),
        fixture.patch,fixture.geometry.global_cells(),&topology,fixture.boundary,
        state,context,system,rows);
    const auto offset=std::size_t(source_cell.x+n*(source_cell.y+n*source_cell.z));
    if (!status) std::cerr<<"IBM source status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';
    if (!expect(bool(status),"IBM prescribed-source enthalpy row prepares")) return false;
    const auto& row=rows[offset];
    const double expected=(old_h-context.dt*F*inlet_h/volume)/new_rho;
    const double computed=old_h+row.rhs/row.diagonal;
    std::cout<<"IBM source axis="<<axis<<" flux="<<mass_flux
             <<" diagonal="<<row.diagonal<<" exact="<<new_rho/context.dt<<'\n';
    passed &= expect(close(row.diagonal,new_rho/context.dt) && close(computed,expected),
        "prescribed IBM enthalpy flux has zero derivative with respect to the fluid state");
    // A lagged mass iterate separates the advective search row from the
    // full conservative row requested by the terminal energy audit.
    const double perturbed_rho=new_rho+0.001;
    rho.view.unchecked(source_cell,0U)=perturbed_rho;
    diagonal.view.unchecked(source_cell,0U)=perturbed_rho*volume/context.dt;
    residual.view.unchecked(source_cell,0U)=
        (perturbed_rho-1.0)*old_h*volume/context.dt+F*inlet_h;
    status=detail::close_cold_enthalpy_rows(fixture.equations.kernels(),
        fixture.patch,fixture.geometry.global_cells(),&topology,fixture.boundary,
        state,context,system,rows,false);
    passed &= expect(bool(status) && close(rows[offset].diagonal,perturbed_rho/context.dt) &&
        close(rows[offset].rhs,-residual.view.unchecked(source_cell,0U)/volume),
        "full IBM enthalpy correction keeps conservative storage and unshifted residual");
  }
  std::cout << "IBM stochastic flux periodic=" << periodic
            << " cases=48 relative_error=" << stochastic_flux_error << '\n';
  return passed;
}

bool test_independent_species_closure() {
  const std::array<double, 3U> independent{0.15, 0.20, 0.25};
  const std::array<double, 3U> fluxes{1.5e-4, -4.0e-5, 7.0e-5};
  SpeciesClosure closure;
  bool passed = expect(static_cast<bool>(close_independent_species(
                           {independent.data(), independent.size()},
                           {fluxes.data(), fluxes.size()}, closure)),
                       "valid N-1 species close without clipping");
  passed &= expect(close(closure.dependent_mass_fraction, 0.40),
                   "dependent species is exactly one minus the independent sum");
  passed &= expect(close(closure.dependent_diffusive_flux, -1.8e-4),
                   "dependent diffusion flux is exactly minus the independent sum");
  for (const std::array<double,3U> fractions :
       {std::array<double,3U>{1.0,1.0e-17,3.0e-72},
        std::array<double,3U>{3.0e-72,1.0e-17,1.0}}) {
    passed &= expect(close_independent_species(
        {fractions.data(),fractions.size()}, {}, closure) &&
        closure.dependent_mass_fraction==0.0,
        "composition closure uses the same represented sum as EOS in either species order");
  }

  std::array<double, 2U> invalid_negative{-1.0e-12, 0.2};
  std::array<double, 2U> invalid_sum{0.8, 0.3};
  const SpeciesClosure sentinel{9.0, 11.0};
  closure = sentinel;
  passed &= expect(close_independent_species(
                       {invalid_negative.data(), invalid_negative.size()}, {},
                       closure).code == StatusCode::numerical_failure &&
                       close(closure.dependent_mass_fraction,
                             sentinel.dependent_mass_fraction),
                   "negative independent species rejects atomically, without clipping");
  closure = sentinel;
  passed &= expect(close_independent_species(
                       {invalid_sum.data(), invalid_sum.size()}, {}, closure)
                           .code == StatusCode::numerical_failure &&
                       close(closure.dependent_mass_fraction,
                             sentinel.dependent_mass_fraction),
                   "independent sum above one rejects atomically, without renormalizing");
  return passed;
}

bool test_composition_dependent_production_eos() {
  ProductionFixture fixture;
  bool passed = expect(make_production_fixture(2, fixture),
                       "composition-dependent production thermo plan compiles");
  if (!passed) {
    return false;
  }

  constexpr double p_ref = 101325.0;
  constexpr double pi = 1750.0;
  constexpr double p_abs = p_ref + pi;
  constexpr double h = 650000.0;
  constexpr double mw_a = 28.0;
  constexpr double mw_b = 32.0;
  constexpr double cp_over_r_a = 3.5;
  constexpr double cp_over_r_b = 4.25;
  const std::array<double, 2U> independent_a_mass_fractions{0.15, 0.75};
  std::array<ThermoState, 2U> states{};
  std::array<double, 2U> oracle_psi{};

  for (std::size_t case_index = 0U;
       case_index < independent_a_mass_fractions.size(); ++case_index) {
    const double y_a = independent_a_mass_fractions[case_index];
    const double y_b = 1.0 - y_a;
    const std::array<double, 1U> independent{y_a};
    const double gas_constant =
        kUniversalGasConstant * (y_a / mw_a + y_b / mw_b);
    const double cp = kUniversalGasConstant *
                      (y_a * cp_over_r_a / mw_a +
                       y_b * cp_over_r_b / mw_b);
    const double temperature = h / cp;
    oracle_psi[case_index] = 1.0 / (gas_constant * temperature);
    const double density = p_abs * oracle_psi[case_index];

    passed &= expect(
        static_cast<bool>(fixture.thermodynamics.evaluate_from_reference_pressure(
            p_ref, pi, h, {independent.data(), independent.size()}, {},
            states[case_index])) &&
            close(states[case_index].rho, density) &&
            close(states[case_index].drho_dp_hY, oracle_psi[case_index]),
        "production EOS matches the independent ideal-gas mixture oracle");
  }

  passed &= expect(
      !close(states[0U].rho, states[1U].rho) &&
          !close(states[0U].drho_dp_hY, states[1U].drho_dp_hY),
      "different valid N-1 compositions change rho and drho/dp through the dependent species closure");
  const double mutated_fixed_r_density = p_abs * oracle_psi[0U];
  passed &= expect(
      !close(mutated_fixed_r_density, states[1U].rho) &&
          !close(oracle_psi[0U], states[1U].drho_dp_hY),
      "a composition-blind fixed-gas-constant mutation is killed numerically");
  return passed;
}

bool test_scalar_catalog_contract() {
  std::array<ScalarEquationSpec, 2U> scalars{{
      {7U, TransportedScalarRole::species, 0.7, 0.9},
      {8U, TransportedScalarRole::passive_scalar, 1.1, 0.8},
  }};
  bool passed = expect(scalars[0].field != scalars[1].field &&
                           scalars[0].role == TransportedScalarRole::species &&
                           scalars[1].role ==
                               TransportedScalarRole::passive_scalar,
                       "species and passive scalar share one typed catalog");
  passed &= expect(scalars[0].molecular_schmidt > 0.0 &&
                       scalars[0].turbulent_schmidt > 0.0,
                   "molecular and turbulent Schmidt numbers are explicit");
  return passed;
}

bool test_scalar_mass_diffusivity_oracle_and_atomicity() {
  const Int3 cells{4, 3, 2};
  OwnedField molecular = make_field(20U, cells, 1U, 0U, 51U);
  OwnedField turbulent = make_field(21U, cells, 1U, 0U, 52U);
  OwnedField diffusivity = make_field(22U, cells, 1U, 0U, 53U);
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        molecular.view.unchecked({i, j, k}, 0U) = 2.0 + i;
        turbulent.view.unchecked({i, j, k}, 0U) = 3.0 + j;
      }
    }
  }
  const ScalarEquationSpec spec{7U, TransportedScalarRole::species, 0.5, 2.0};
  bool passed = expect(static_cast<bool>(form_scalar_mass_diffusivity(
                           spec, as_const(molecular.view),
                           as_const(turbulent.view), {{0, 0, 0}, cells},
                           diffusivity.view)),
                       "scalar diffusivity is formed by the production hot path");
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const double expected = (2.0 + i) / 0.5 + (3.0 + j) / 2.0;
        passed &= expect(close(diffusivity.view.unchecked({i, j, k}, 0U),
                               expected),
                         "rhoD=mu/Sc+mu_t/Sc_t exactly");
      }
    }
  }
  std::fill(diffusivity.bytes.begin(), diffusivity.bytes.end(), 91.0);
  turbulent.view.unchecked({2, 1, 1}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  passed &= expect(form_scalar_mass_diffusivity(
                       spec, as_const(molecular.view),
                       as_const(turbulent.view), {{0, 0, 0}, cells},
                       diffusivity.view).code == StatusCode::numerical_failure,
                   "non-finite turbulent viscosity rejects diffusivity formation");
  passed &= expect(std::all_of(diffusivity.bytes.begin(), diffusivity.bytes.end(),
                               [](double value) { return value == 91.0; }),
                   "diffusivity formation failure leaves all output bytes untouched");
  return passed;
}

bool test_outward_flux_sign_contract() {
  double heat_gradient = 0.0;
  double scalar_gradient = 0.0;
  bool passed = expect(static_cast<bool>(resolve_heat_flux_normal_gradient(
                           12.0, 3.0, heat_gradient)) &&
                           close(heat_gradient, -4.0),
                       "positive outward heat flux gives negative outward T gradient");
  passed &= expect(static_cast<bool>(resolve_scalar_flux_normal_gradient(
                       10.0, 2.0, scalar_gradient)) &&
                       close(scalar_gradient, -5.0),
                   "positive outward scalar flux gives negative outward q gradient");
  return passed;
}

bool test_conservative_global_balance_oracle() {
  // A one-cell finite-volume balance with supplied final outward mass fluxes.
  const std::array<double, 6U> outward_mdot{-0.5, 0.8, -0.1,
                                            0.2, -0.4, 0.3};
  constexpr double y = 0.25;
  long double boundary = 0.0L;
  for (const double flux : outward_mdot) {
    boundary += flux * y;
  }
  constexpr double rho_volume_old = 1.2;
  constexpr double dt = 0.1;
  const double rho_y_volume_new = rho_volume_old * y -
                                  dt * static_cast<double>(boundary);
  const double residual =
      (rho_y_volume_new - rho_volume_old * y) / dt +
      static_cast<double>(boundary);
  return expect(close(residual, 0.0),
                "N-1 transport global change equals supplied boundary flux");
}

bool output_is(const OwnedField& diagonal, const OwnedField& rhs,
               const OwnedField& residual, double value) {
  const auto matches = [value](const std::vector<double>& bytes) {
    return std::all_of(bytes.begin(), bytes.end(),
                       [value](double actual) { return actual == value; });
  };
  return matches(diagonal.bytes) && matches(rhs.bytes) &&
         matches(residual.bytes);
}

bool faces_are(const OwnedFaceField& x, const OwnedFaceField& y,
               const OwnedFaceField& z, double value) {
  const auto matches = [value](const std::vector<double>& bytes) {
    return std::all_of(bytes.begin(), bytes.end(),
                       [value](double actual) { return actual == value; });
  };
  return matches(x.bytes) && matches(y.bytes) && matches(z.bytes);
}

bool same_certificate(const EquationAssemblyCertificate& left,
                      const EquationAssemblyCertificate& right) {
  return left.plan == right.plan && left.scope == right.scope &&
         left.time == right.time && left.geometry == right.geometry &&
         left.face_flux == right.face_flux && left.state == right.state &&
         left.dt == right.dt;
}

bool test_production_species_and_passive_assembly() {
  constexpr std::int32_t n = 4;
  ProductionFixture fixture;
  bool passed = expect(make_production_fixture(n, fixture),
                       "species/passive production EquationPlanSet compiles");
  if (!passed) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  const double dx = 1.0 / static_cast<double>(n);
  const double volume = dx * dx * dx;

  OwnedField rho = make_field(kDensity, cells, 1U, 2U, 601U);
  OwnedField species = make_field(kSpecies, cells, 1U, 2U, 602U);
  OwnedField passive = make_field(kPassive, cells, 1U, 2U, 603U);
  OwnedField molecular = make_field(20U, cells, 1U, 1U, 604U);
  OwnedField turbulent = make_field(21U, cells, 1U, 1U, 605U);
  OwnedField species_diffusivity =
      make_field(22U, cells, 1U, 1U, 606U);
  OwnedField passive_diffusivity =
      make_field(23U, cells, 1U, 1U, 607U);
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const Int3 cell{i, j, k};
        const double x = (static_cast<double>(i) + 0.5) * dx;
        rho.view.unchecked(cell, 0U) = 1.0;
        species.view.unchecked(cell, 0U) = 0.20 + 0.10 * x;
        passive.view.unchecked(cell, 0U) = 0.30 - 0.04 * x;
        if (i >= -1 && i < cells.x + 1 && j >= -1 &&
            j < cells.y + 1 && k >= -1 && k < cells.z + 1) {
          molecular.view.unchecked(cell, 0U) = 2.0;
          turbulent.view.unchecked(cell, 0U) = 4.0;
          species_diffusivity.view.unchecked(cell, 0U) = 6.0;
          passive_diffusivity.view.unchecked(cell, 0U) = 2.0;
        }
      }
    }
  }
  fill_field(species_diffusivity, 6.0);
  fill_field(passive_diffusivity, 2.0);
  const ScalarEquationSpec* species_spec = fixture.equations.species().spec(0U);
  const ScalarEquationSpec* passive_spec = fixture.equations.scalars().spec(0U);
  passed &= expect(species_spec != nullptr && passive_spec != nullptr &&
                       static_cast<bool>(form_scalar_mass_diffusivity(
                           *species_spec, as_const(molecular.view),
                           as_const(turbulent.view), {{0, 0, 0}, cells},
                           species_diffusivity.view)) &&
                       static_cast<bool>(form_scalar_mass_diffusivity(
                           *passive_spec, as_const(molecular.view),
                           as_const(turbulent.view), {{0, 0, 0}, cells},
                           passive_diffusivity.view)),
                   "production Schmidt laws form both scalar diffusivities");
  passed &= expect(close(species_diffusivity.view.unchecked({1, 1, 1}, 0U),
                         2.0 / 0.5 + 4.0 / 2.0) &&
                       close(passive_diffusivity.view.unchecked({1, 1, 1}, 0U),
                             2.0 / 2.0 + 4.0 / 4.0),
                   "species and passive operators use their own Sc and Sc_t");

  FinalFluxFixture final_flux;
  ConstFaceFluxView committed;
  passed &= expect(make_linear_final_flux(fixture.equations.kernels(), cells,
                                          final_flux, committed),
                   "FinalFaceFluxWriter commits the production mass flux");
  if (!passed) {
    return false;
  }

  const PrimitiveHistory density_history{as_const(rho.view),
                                          as_const(rho.view),
                                          as_const(rho.view)};
  const PrimitiveHistory species_history{as_const(species.view),
                                          as_const(species.view),
                                          as_const(species.view)};
  const PrimitiveHistory passive_history{as_const(passive.view),
                                          as_const(passive.view),
                                          as_const(passive.view)};
  const std::array independent{species_history};
  const std::array passives{passive_history};
  EquationStateView state;
  state.density = density_history;
  state.independent_species = {independent.data(), independent.size()};
  state.passive_scalars = {passives.data(), passives.size()};
  const std::array diffusivities{as_const(species_diffusivity.view),
                                  as_const(passive_diffusivity.view)};
  EquationMaterialView material;
  material.scalar_mass_diffusivity = {diffusivities.data(),
                                      diffusivities.size()};
  EquationAssemblyContext context;
  context.dt = 0.1;
  context.bdf = {10.0, -10.0, 0.0, 1U};
  context.time = 701U;
  context.geometry = fixture.geometry.topology_revision();
  context.boundary = fixture.boundary.revision();
  context.thermo = fixture.thermodynamics.fingerprint();
  context.transport = fixture.transport.fingerprint();
  context.contribution_stage = 1U;
  context.face_flux = committed.revision;
  context.face_flux_authority = committed.certificate.authority();
  context.face_flux_storage = committed.certificate.storage();
  context.face_flux_revision_domain = committed.certificate.revision_domain();
  context.scope = EquationAssemblyScope::final_conservative;
  context.mass_flux = committed;
  context.provisional_mass_flux = false;

  OwnedField diagonal = make_field(30U, cells, 1U, 0U, 610U);
  OwnedField rhs = make_field(31U, cells, 1U, 0U, 611U);
  OwnedField residual = make_field(32U, cells, 1U, 0U, 612U);
  OwnedFaceField ax = make_face_field(CartesianAxis::x, cells, 613U);
  OwnedFaceField ay = make_face_field(CartesianAxis::y, cells, 614U);
  OwnedFaceField az = make_face_field(CartesianAxis::z, cells, 615U);
  EquationSystemView system{diagonal.view, rhs.view, residual.view,
                            ax.view, ay.view, az.view};

  EquationAssemblyCertificate species_certificate;
  passed &= expect(static_cast<bool>(assemble_species(
                       fixture.equations.species(), 0U, state, material, {},
                       context, system, species_certificate)),
                   "production N-1 species assembly succeeds");
  long double species_sum = 0.0L;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        const double x = (static_cast<double>(i) + 0.5) * dx;
        const double q = 0.20 + 0.10 * x;
        const double expected_residual = (0.20 + 0.20 * x) * volume;
        const double expected_diagonal = 10.0 * volume + 6.0 * 6.0 * dx;
        passed &= expect(close(residual.view.unchecked(cell, 0U),
                               expected_residual) &&
                             close(diagonal.view.unchecked(cell, 0U),
                                   expected_diagonal) &&
                             close(rhs.view.unchecked(cell, 0U),
                                   expected_diagonal * q - expected_residual),
                         "N-1 species cell-integral residual/operator matches oracle");
        species_sum += residual.view.unchecked(cell, 0U);
      }
    }
  }
  long double boundary_species = 0.0L;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      const double q_min = 0.5 *
          (species.view.unchecked({-1, j, k}, 0U) +
           species.view.unchecked({0, j, k}, 0U));
      const double q_max = 0.5 *
          (species.view.unchecked({cells.x - 1, j, k}, 0U) +
           species.view.unchecked({cells.x, j, k}, 0U));
      boundary_species += committed.x.unchecked({cells.x, j, k}) * q_max -
                          committed.x.unchecked({0, j, k}) * q_min;
    }
  }
  passed &= expect(close(static_cast<double>(species_sum),
                         static_cast<double>(boundary_species)),
                   "global N-1 residual equals the committed boundary flux");
  passed &= expect(species_certificate.valid() &&
                       close(ax.view.unchecked({2, 1, 1}), 6.0 * dx),
                   "N-1 assembly publishes certificate and Schmidt coefficient");

  EquationAssemblyCertificate passive_certificate;
  passed &= expect(static_cast<bool>(assemble_scalar(
                       fixture.equations.scalars(), 0U, state, material, {},
                       context, system, passive_certificate)),
                   "production passive-scalar assembly succeeds");
  const Int3 probe{2, 1, 1};
  const double probe_x = (static_cast<double>(probe.x) + 0.5) * dx;
  const double probe_q = 0.30 - 0.04 * probe_x;
  const double passive_residual = (0.30 - 0.08 * probe_x) * volume;
  const double passive_diagonal = 10.0 * volume + 6.0 * 2.0 * dx;
  passed &= expect(close(residual.view.unchecked(probe, 0U),
                         passive_residual) &&
                       close(diagonal.view.unchecked(probe, 0U),
                             passive_diagonal) &&
                       close(rhs.view.unchecked(probe, 0U),
                             passive_diagonal * probe_q - passive_residual) &&
                       close(ax.view.unchecked({2, 1, 1}), 2.0 * dx) &&
                       passive_certificate.valid(),
                   "passive residual/operator uses its independent Schmidt law");

  const auto reset_outputs = [&]() {
    fill_field(diagonal, 91.0);
    fill_field(rhs, 91.0);
    fill_field(residual, 91.0);
    std::fill(ax.bytes.begin(), ax.bytes.end(), 91.0);
    std::fill(ay.bytes.begin(), ay.bytes.end(), 91.0);
    std::fill(az.bytes.begin(), az.bytes.end(), 91.0);
  };

  reset_outputs();
  EquationAssemblyCertificate rejected = passive_certificate;
  PrimitiveHistory aliased_passive = passive_history;
  aliased_passive.accepted = as_const(residual.view);
  aliased_passive.accepted.field = kPassive;
  aliased_passive.accepted.revision = passive_history.accepted.revision;
  const std::array aliased_passives{aliased_passive};
  EquationStateView aliased_state = state;
  aliased_state.passive_scalars = {aliased_passives.data(),
                                   aliased_passives.size()};
  passed &= expect(
      assemble_scalar(fixture.equations.scalars(), 0U, aliased_state,
                      material, {}, context, system, rejected)
                  .code == StatusCode::invalid_plan &&
          output_is(diagonal, rhs, residual, 91.0) &&
          faces_are(ax, ay, az, 91.0) &&
          same_certificate(rejected, passive_certificate),
      "zero-contribution scalar history/output alias rejects atomically");

  reset_outputs();
  rejected = passive_certificate;
  context.box = {{1, 0, 0}, {cells.x - 1, cells.y, cells.z}};
  passed &= expect(
      assemble_scalar(fixture.equations.scalars(), 0U, state, material, {},
                      context, system, rejected)
                  .code == StatusCode::invalid_plan &&
          output_is(diagonal, rhs, residual, 91.0) &&
          faces_are(ax, ay, az, 91.0) &&
          same_certificate(rejected, passive_certificate),
      "zero-contribution partial-box scalar assembly rejects atomically");
  context.box = {};

  reset_outputs();
  rejected = passive_certificate;
  ConstFaceFluxView stale = committed;
  ++stale.revision;
  context.face_flux = stale.revision;
  context.mass_flux = stale;
  passed &= expect(assemble_scalar(fixture.equations.scalars(), 0U, state,
                                   material, {}, context, system, rejected)
                           .code == StatusCode::invalid_plan &&
                       output_is(diagonal, rhs, residual, 91.0) &&
                       faces_are(ax, ay, az, 91.0) &&
                       rejected.plan == passive_certificate.plan,
                   "stale final authority rejects atomically");

  FaceFluxStorage workspace;
  FaceFluxView provisional;
  reset_outputs();
  rejected = passive_certificate;
  passed &= expect(static_cast<bool>(FaceFluxStorage::allocate_workspace(
                           cells, 1U, workspace)) &&
                       static_cast<bool>(workspace.workspace_view(
                           0U, 777U, provisional)),
                   "provisional flux mutation allocates");
  context.face_flux = 777U;
  context.face_flux_authority = 0U;
  context.face_flux_storage = provisional.x.storage_identity;
  context.face_flux_revision_domain = provisional.x.revision_domain;
  context.mass_flux = as_const(provisional);
  context.provisional_mass_flux = true;
  passed &= expect(assemble_scalar(fixture.equations.scalars(), 0U, state,
                                   material, {}, context, system, rejected)
                           .code == StatusCode::invalid_plan &&
                       output_is(diagonal, rhs, residual, 91.0) &&
                       faces_are(ax, ay, az, 91.0) &&
                       rejected.plan == passive_certificate.plan,
                   "provisional final flux rejects atomically");

  context.scope=EquationAssemblyScope::momentum_predictor;
  rejected=species_certificate;
  passed &= expect(assemble_species(fixture.equations.species(),0U,state,
      material,{},context,system,rejected).code==StatusCode::invalid_plan &&
      output_is(diagonal,rhs,residual,91.0) && faces_are(ax,ay,az,91.0) &&
      rejected.plan==species_certificate.plan,
      "public species assembly does not promote a provisional guess");
  for(const auto face:{provisional.x,provisional.y,provisional.z})
    for(int z=0;z<face.extents.z;++z) for(int y=0;y<face.extents.y;++y)
      for(int x=0;x<face.extents.x;++x) face.unchecked({x,y,z})=0.0;
  EquationAssemblyCertificate guess_certificate;
  passed &= expect(static_cast<bool>(detail::assemble_species_guess(
      fixture.equations.species(),0U,state,material,context,system,guess_certificate)) &&
      guess_certificate.valid() && guess_certificate.scope==EquationAssemblyScope::momentum_predictor,
      "private species seed assembles without final conservation authority");

  reset_outputs();
  rejected = species_certificate;
  context.face_flux = committed.revision;
  context.face_flux_authority = committed.certificate.authority();
  context.face_flux_storage = committed.certificate.storage();
  context.face_flux_revision_domain = committed.certificate.revision_domain();
  context.mass_flux = committed;
  context.provisional_mass_flux = false;
  context.scope=EquationAssemblyScope::final_conservative;
  passed &= expect(detail::assemble_species_guess(fixture.equations.species(),
      0U,state,material,context,system,rejected).code==StatusCode::invalid_plan &&
      output_is(diagonal,rhs,residual,91.0) && faces_are(ax,ay,az,91.0) &&
      rejected.plan==species_certificate.plan,
      "private species seed cannot impersonate the final assembler");
  passed &= expect(static_cast<bool>(detail::assemble_species_coupling_rows(
      fixture.equations.species(),0U,state,material,context,system)),
      "private density rows require and accept real final-flux authority");
  for(int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x) {
    const Int3 cell{x,y,z};
    const double coordinate=(x+0.5)*dx;
    passed &= expect(close(residual.view.unchecked(cell,0U),0.20+0.20*coordinate) &&
        close(diagonal.view.unchecked(cell,0U),10.0+6.0*6.0*dx/volume),
        "private density row has the same physical equation and diagonal before volume scaling");
  }
  // The nonlinear iterate changes while the prepared material and density
  // remain fixed. Compare the retained-diagonal path against full assembly.
  auto varied = make_field(species.view.field, cells, 1U, 2U, 904U);
  for (std::size_t i = 0; i < species.bytes.size(); ++i)
    varied.bytes[i] = 0.9 * species.bytes[i];
  auto varied_history = species_history;
  varied_history.trial = as_const(varied.view);
  const std::array varied_species{varied_history};
  auto varied_state = state;
  varied_state.independent_species = {varied_species.data(), varied_species.size()};
  const auto saved_diagonal = diagonal.bytes;
  const auto original_residual = residual.bytes;
  std::fill(ax.bytes.begin(), ax.bytes.end(), 91.0);
  std::fill(ay.bytes.begin(), ay.bytes.end(), 91.0);
  std::fill(az.bytes.begin(), az.bytes.end(), 91.0);
  passed &= expect(static_cast<bool>(detail::assemble_species_coupling_residual(
      fixture.equations.species(), 0U, varied_state, material, context, system)) &&
      diagonal.bytes == saved_diagonal && faces_are(ax, ay, az, 91.0) &&
      residual.bytes != original_residual,
      "retained diagonal evaluates the new iterate and preserves coefficient storage");
  const auto cached_residual = residual.bytes;
  const auto cached_rhs = rhs.bytes;
  passed &= expect(static_cast<bool>(detail::assemble_species_coupling_rows(
      fixture.equations.species(), 0U, varied_state, material, context, system)) &&
      residual.bytes == cached_residual && rhs.bytes == cached_rhs &&
      diagonal.bytes == saved_diagonal,
      "retained diagonal matches the full physical species rows exactly");
  std::fill(residual.bytes.begin(), residual.bytes.end(), 91.0);
  std::fill(rhs.bytes.begin(), rhs.bytes.end(), 91.0);
  diagonal.view.unchecked({1, 1, 1}, 0U) = -1.0;
  passed &= expect(detail::assemble_species_coupling_residual(
      fixture.equations.species(), 0U, varied_state, material, context, system).code ==
      StatusCode::numerical_failure &&
      std::all_of(residual.bytes.begin(), residual.bytes.end(), [](double v) { return v == 91.0; }) &&
      std::all_of(rhs.bytes.begin(), rhs.bytes.end(), [](double v) { return v == 91.0; }),
      "invalid retained diagonal rejects before residual and rhs writes");

  reset_outputs();
  const auto saved_scope=context.scope;
  context.scope=EquationAssemblyScope::target_coupled;
  passed &= expect(detail::assemble_species_coupling_rows(fixture.equations.species(),
      0U,state,material,context,system).code==StatusCode::invalid_plan &&
      output_is(diagonal,rhs,residual,91.0),
      "private density rows cannot bypass final authority through target scope");
  context.scope=saved_scope;
  species.view.unchecked({1, 1, 1}, 0U) = 1.01;
  passed &= expect(assemble_species(fixture.equations.species(), 0U, state,
                                    material, {}, context, system, rejected)
                           .code == StatusCode::numerical_failure &&
                       output_is(diagonal, rhs, residual, 91.0) &&
                       faces_are(ax, ay, az, 91.0) &&
                       rejected.plan == species_certificate.plan,
                   "invalid N-1 composition rejects before all writes");
  return passed;
}

bool test_reacting_registration_stays_out_of_scope() {
  const std::array<FieldId, 3U> declared{0U, 1U, 2U};
  ContributionRegistry registry;
  bool passed = expect(static_cast<bool>(
                           registry.configure({declared.data(), declared.size()})),
                       "inert contribution registry configures");
  ContributionSpec chemistry;
  chemistry.conserved_quantity = 0U;
  chemistry.stage = 1U;
  chemistry.explicit_source = 1U;
  chemistry.capability = ContributionCapability::chemistry;
  passed &= expect(registry.register_contribution(chemistry).code ==
                       StatusCode::invalid_plan,
                   "v0.4 species equations reject chemistry registration");
  return passed;
}

bool test_statistical_thermal_operator() {
  ProductionFixture fixture;if(!make_production_fixture(8,fixture))return false;
  const auto cells=fixture.patch.cells;const auto& kernels=fixture.equations.kernels();
  auto total=make_field(94,cells,1,2,511);
  const Int3 face{3,3,3};double rates[2]{};
  for(unsigned side=0;side<2;++side) {
    for(int z=-2;z<cells.z+2;++z)for(int y=-2;y<cells.y+2;++y)for(int x=-2;x<cells.x+2;++x)
      total.view.unchecked({x,y,z},0)=x<2 ? 0. : x==2 ? 1. : x==3 ? 1.+(side ? 1e-10 : -1e-10) : 2.;
    MixtureFaceTransport chosen;
    const auto status=prepare_cartesian_mixture_face(kernels,{},as_const(total.view),CartesianAxis::x,face,.04,.0025,chosen,true);
    if(!status)return false;
    // The pressure-work increment has a finite gradient while total h
    // crosses its neighbour. Its source anchor carries the other gradient.
    rates[side]=-(.0025+chosen.extra_diffusion)*50.;
  }
  const double jump=std::abs(rates[1]-rates[0]);
  std::cout<<"statistical_thermal_total_gradient_crossing jump_W="<<jump<<'\n';
  bool passed=expect(jump>.8,"total-h limiter reproduces the finite statistical thermal flux jump");
  auto gamma=make_field(95,cells,1,1,512);fill_field(gamma,.02);
  std::array<OwnedFaceField,3> old_mass,new_mass,extra,workspace;
  for(unsigned axis=0;axis<3;++axis) {
    old_mass[axis]=make_face_field(static_cast<CartesianAxis>(axis),cells,8100);
    new_mass[axis]=make_face_field(static_cast<CartesianAxis>(axis),cells,8101);
    extra[axis]=make_face_field(static_cast<CartesianAxis>(axis),cells,8102+4*axis);
    workspace[axis]=make_face_field(static_cast<CartesianAxis>(axis),cells,8103);
    std::fill(old_mass[axis].bytes.begin(),old_mass[axis].bytes.end(),.04);
    std::fill(new_mass[axis].bytes.begin(),new_mass[axis].bytes.end(),.04);
  }
  ConstFaceFluxView old_flux{as_const(old_mass[0].view),as_const(old_mass[1].view),as_const(old_mass[2].view),711};
  ConstFaceFluxView new_flux{as_const(new_mass[0].view),as_const(new_mass[1].view),as_const(new_mass[2].view),712};
  FaceFluxView output{workspace[0].view,workspace[1].view,workspace[2].view,713};
  MixtureTransportFaces mixture{as_const(extra[0].view),as_const(extra[1].view),as_const(extra[2].view),711,714};
  detail::StatisticalEnergyLedger ledger;
  auto status=ledger.initialize(cells,2,715,fixture.equations.enthalpy());
  for(double coefficient:{.01,.03}) {
    for(auto& f:extra)std::fill(f.bytes.begin(),f.bytes.end(),coefficient);
    if(status)status=ledger.freeze(kernels,as_const(total.view),as_const(gamma.view),old_flux,&mixture,nullptr);
  }
  if(!expect(bool(status),"statistical thermal operator freezes both field coefficients"))return false;
  MixtureTransportFaces selected;
  double fixed_rates[2]{};
  for(unsigned side=0;side<2;++side) {
    total.view.unchecked(face,0)=1.+(side ? 1e-10 : -1e-10);
    status=ledger.correction_mixture(kernels,fixture.boundary,old_flux,new_flux,output,716+side,selected);
    if(!status)return false;
    fixed_rates[side]=-(.0025+selected.x.unchecked(face))*50.;
  }
  passed &= expect(fixed_rates[0]==fixed_rates[1] && close(selected.x.unchecked(face),.02),
      "frozen field-mean coefficients keep the thermal increment flux continuous");
  for(double increment:{-.02,0.,.02}) {
    new_mass[0].view.unchecked(face)=.04+increment;
    status=ledger.correction_mixture(kernels,fixture.boundary,old_flux,new_flux,output,719,selected);
    if(!status)return false;
    const double left=4.,right=54.;
    const double actual=(.04+increment)*.5*(left+right)-(.0025+selected.x.unchecked(face))*(right-left);
    const double independent=.04*.5*(left+right)-(.0025+.02)*(right-left)+increment*(increment>=0 ? left : right);
    passed &= expect(close(actual,independent),"common heat increment equals mean frozen field flux plus actual pressure upwind flux");
  }
  auto invalid=output;invalid.x=new_mass[0].view;
  passed &= expect(!ledger.correction_mixture(kernels,fixture.boundary,old_flux,new_flux,invalid,720,selected),
      "thermal operator rejects flux input alias");
  const auto bytes=ledger.owned_payload_bytes();
  status=ledger.initialize(cells,2,721,fixture.equations.enthalpy());
  passed &= expect(bool(status) && ledger.owned_payload_bytes()==bytes &&
      !ledger.correction_mixture(kernels,fixture.boundary,old_flux,new_flux,output,722,selected),
      "next attempt reuses allocation and awaits its own complete frozen fields");
  for(double coefficient:{.04,.08}) {
    for(auto& f:extra)std::fill(f.bytes.begin(),f.bytes.end(),coefficient);
    if(status)status=ledger.freeze(kernels,as_const(total.view),as_const(gamma.view),old_flux,&mixture,nullptr);
  }
  new_mass[0].view.unchecked(face)=.04;
  if(status)status=ledger.correction_mixture(kernels,fixture.boundary,old_flux,new_flux,output,723,selected);
  passed &= expect(bool(status) && close(selected.x.unchecked(face),.06),
      "reused ledger averages only this attempt's fields");
  ledger.discard_attempt();
  passed &= expect(!ledger.correction_mixture(kernels,fixture.boundary,old_flux,new_flux,output,724,selected),
      "discarded attempt exposes no thermal correction operator");
  std::cout<<"statistical_thermal_frozen jump_W="<<std::abs(fixed_rates[1]-fixed_rates[0])<<'\n';
  return passed;
}

bool test_mixture_flux_ledger() {
  if(!test_statistical_thermal_operator())return false;
  ProductionFixture fixture;
  if(!make_production_fixture(8,fixture))return false;
  const auto cells=fixture.patch.cells;
  const auto& kernels=fixture.equations.kernels();
  auto q=make_field(90,cells,1,1,501);
  auto gamma=make_field(91,cells,1,1,502);
  auto divergence=make_field(92,cells,1,0,503);
  fill_field(gamma,.02);
  std::array<OwnedFaceField,3> mass,extra,rate;
  for(unsigned a=0;a<3;++a) {
    mass[a]=make_face_field(static_cast<CartesianAxis>(a),cells,8000);
    extra[a]=make_face_field(static_cast<CartesianAxis>(a),cells,8001+3*a);
    rate[a]=make_face_field(static_cast<CartesianAxis>(a),cells,8002+3*a);
    auto end=cells;++(a==0 ? end.x : a==1 ? end.y : end.z);
    for(int z=0;z<end.z;++z)for(int y=0;y<end.y;++y)for(int x=0;x<end.x;++x) {
      const int normal=a==0 ? x : a==1 ? y : z;
      mass[a].view.unchecked({x,y,z})=.03*(normal-4);
      extra[a].view.unchecked({x,y,z})=.002+.001*((x+y+z)%3);
    }
  }
  ConstFaceFluxView flux{as_const(mass[0].view),as_const(mass[1].view),as_const(mass[2].view),601};
  MixtureTransportFaces mixture{as_const(extra[0].view),as_const(extra[1].view),as_const(extra[2].view),601,701};
  std::array<FaceFieldView,3> output{rate[0].view,rate[1].view,rate[2].view};
  bool passed=true;double analytic_error{},parity_error{},balance_error{};
  for(double formation:{0.,-4e6,4e6}) {
    const auto exact=[&](double x,double y,double z) {return 301400.+formation*(.3+.02*x-.03*y+.01*z);};
    for(int z=-1;z<=cells.z;++z)for(int y=-1;y<=cells.y;++y)for(int x=-1;x<=cells.x;++x)
      q.view.unchecked({x,y,z},0)=exact((x+.5)/8.,(y+.5)/8.,(z+.5)/8.);
    auto status=form_cartesian_mixture_transport_flux(kernels,mixture,as_const(gamma.view),flux,as_const(q.view),output);
    if(!expect(bool(status),"statistical face ledger exports common transport"))return false;
    for(unsigned a=0;a<3;++a) {
      auto end=cells;++(a==0 ? end.x : a==1 ? end.y : end.z);
      for(int z=0;z<end.z;++z)for(int y=0;y<end.y;++y)for(int x=0;x<end.x;++x) {
        const Int3 face{x,y,z};double r[3]{(x+.5)/8.,(y+.5)/8.,(z+.5)/8.};r[a]-=.5/8.;
        const double slope=a==0 ? .02 : a==1 ? -.03 : .01;
        const double expected=mass[a].view.unchecked(face)*exact(r[0],r[1],r[2])-
            (.02/8.+extra[a].view.unchecked(face))*formation*slope/8.;
        analytic_error=std::max(analytic_error,std::abs(output[a].unchecked(face)-expected)/std::max(1.,std::abs(expected)));
      }
    }
    const auto read=as_const(q.view);
    const KernelInvocation invocation{{&read,1},{&divergence.view,1},{{0,0,0},cells},0,0,1,601,nullptr};
    status=cartesian_mixture_transport(kernels,mixture,as_const(gamma.view),flux,invocation);
    if(!expect(bool(status),"common divergence evaluates the same state"))return false;
    long double sum{},boundary{},scale{};
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
      const Int3 c{x,y,z};long double integral{};
      for(unsigned a=0;a<3;++a) {
        Int3 upper=c;++(a==0 ? upper.x : a==1 ? upper.y : upper.z);
        integral+=static_cast<long double>(output[a].unchecked(upper))-output[a].unchecked(c);
        const auto normal=a==0 ? x : a==1 ? y : z;
        if(normal==0)boundary-=output[a].unchecked(c);
        if(normal==7)boundary+=output[a].unchecked(upper);
      }
      const double replay=static_cast<double>(integral*512.);
      parity_error=std::max(parity_error,std::abs(replay-divergence.view.unchecked(c,0)));
      sum+=divergence.view.unchecked(c,0)/512.;scale+=std::abs(integral);
    }
    balance_error=std::max(balance_error,double(std::abs(sum-boundary)/std::max(1.L,scale)));
  }
  passed &= expect(analytic_error<1e-12 && parity_error==0 && balance_error<1e-13,
      "face ledger matches affine analytic flux, bitwise divergence and independent boundary balance");
  for(auto& f:rate)std::fill(f.bytes.begin(),f.bytes.end(),123.);
  auto invalid=mixture;++invalid.face_flux;
  passed &= expect(!form_cartesian_mixture_transport_flux(kernels,invalid,as_const(gamma.view),flux,as_const(q.view),output),
      "stale face generation rejects");
  auto alias=output;alias[0]=mass[0].view;
  passed &= expect(!form_cartesian_mixture_transport_flux(kernels,mixture,as_const(gamma.view),flux,as_const(q.view),alias),
      "mass-flux input alias rejects");
  alias=output;alias[0]=extra[0].view;
  passed &= expect(!form_cartesian_mixture_transport_flux(kernels,mixture,as_const(gamma.view),flux,as_const(q.view),alias),
      "frozen coefficient input alias rejects");
  passed &= expect(faces_are(rate[0],rate[1],rate[2],123.),"invalid ledger requests preserve output");
  std::cout<<"mixture_flux_ledger analytic="<<analytic_error<<" divergence="<<parity_error<<" balance="<<balance_error<<'\n';
  return passed;
}

}  // namespace

bool test_implicit_mixture_solve(MPI_Comm comm, bool frozen=false, bool heat=false) {
  constexpr int global_n=9;
  constexpr double dt=.1,rho_value=1.,mu=.1;
  const double mean=heat ? 1e5 : .5,amplitude=heat ? 2e4 : .1;
  const double scale=heat ? 1e5 : 1.;
  const FieldId scalar_field=heat ? kEnthalpy : kSpecies;
  const FieldId source_field=heat ? 12U : 10U,sink_field=heat ? 13U : 11U;
  const auto boundary_stage=heat ? BoundaryStage::enthalpy : BoundaryStage::scalar;
  constexpr double dx=1./global_n,volume=dx*dx*dx;
  ProductionFixture fixture;
  if (!make_production_fixture(global_n,fixture,false,-1,true,volume,comm)) return false;
  const auto cells=fixture.patch.cells;
  const auto count=std::size_t(cells.x)*cells.y*cells.z;
  auto rho=make_field(kDensity,cells,1,2,1001);
  auto q=make_field(scalar_field,cells,1,2,1002);
  auto old=make_field(scalar_field,cells,1,2,1003);
  auto velocity=make_field(kVelocity,cells,3,2,1004);
  auto gamma=make_field(25,cells,1,2,1005);
  auto source=make_field(source_field,cells,1,0,1006);
  auto sink=make_field(sink_field,cells,1,0,1007);
  auto diagonal=make_field(30,cells,1,0,1008);
  auto rhs=make_field(31,cells,1,0,1009);
  auto residual=make_field(32,cells,1,0,1010);
  auto ax=make_face_field(CartesianAxis::x,cells,1011);
  auto ay=make_face_field(CartesianAxis::y,cells,1012);
  auto az=make_face_field(CartesianAxis::z,cells,1013);
  auto variation=make_field(90,cells,1,1,1014);
  auto correction=make_field(91,cells,1,1,1015);
  auto linear_rhs=make_field(94,cells,1,0,1016);
  auto mass_divergence=make_field(97,cells,1,0,1017);
  auto backup=make_field(98,cells,1,2,1018);
  fill_field(rho,rho_value);fill_field(gamma,mu/.7);
  esf::detail::IemSource terms;
  if (esf::detail::iem_source(volume,mu,0,2,1,mean,terms)!=portable::Status::success) return false;
  fill_field(source,terms.explicit_source_density);fill_field(sink,terms.implicit_sink_density);
  LinearWorkspaceRequirements requirements;
  auto status=make_linear_workspace_requirements(LinearAlgorithm::fgmres,cells,1,24,
      ReductionMode::mpi_allreduce,1020,requirements);
  ReductionEngine reductions;
  if (status) status=ReductionEngine::compile(comm,ReductionMode::mpi_allreduce,
      requirements.reduction_capacity,reductions);
  if (!status) {std::cerr<<"IEM solve requirements="<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
  auto vectors=make_field(92,cells,requirements.vector_slots,1,1021);
  auto scalars=make_field(93,{int(requirements.scalar_doubles),1,1},1,0,1022);
  scalars.view.storage_identity=vectors.view.storage_identity;
  SolverWorkspace workspace;
  status=SolverWorkspace::bind(requirements,vectors.view,scalars.view,workspace);
  HaloEngine halo, state_halo;
  const std::array<HaloFieldSpec,1> halo_fields{{{92,1,1}}};
  if (status) status=halo.reserve(comm,fixture.patch,{halo_fields.data(),halo_fields.size()},
      fixture.boundary.halo_topology());
  const HaloFieldSpec state_halo_field{scalar_field,2,1};
  if (status) status=state_halo.reserve(comm,fixture.patch,{&state_halo_field,1},
      fixture.boundary.halo_topology());
  if (!status) {std::cerr<<"IEM solve bind="<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
  std::vector<detail::ColdPressureRow> rows(count);
  const double kdx=2*std::acos(-1.)/global_n;
  double errors[3]{};
  unsigned iterations{};
  for (double speed : {-.2,.2}) {
    FinalFluxFixture flux_owner;
    ConstFaceFluxView flux;
    if (!make_linear_final_flux(fixture.equations.kernels(),cells,flux_owner,flux,false,false,speed,
                               frozen ? speed : 0)) {
      std::cerr<<"IEM solve final flux\n";return false;
    }
    if (frozen) {
      KernelInvocation call{{},{&mass_divergence.view,1},{{0,0,0},cells},0,0,1,flux.revision};
      status=cartesian_face_divergence(fixture.equations.kernels(),flux,call);
      if (!status) return false;
    }
    fill_field(velocity,0.);
    for (int z=-2;z<cells.z+2;++z) for (int y=-2;y<cells.y+2;++y) for (int x=-2;x<cells.x+2;++x)
      velocity.view.unchecked({x,y,z},0)=speed;
    for (double sign : {-1.,1.}) {
      for (int z=-2;z<cells.z+2;++z) for (int y=-2;y<cells.y+2;++y) for (int x=-2;x<cells.x+2;++x) {
        const double value=mean+sign*amplitude*(frozen ? 1. : std::cos(kdx*(x+fixture.patch.begin.x+.5)));
        q.view.unchecked({x,y,z},0)=old.view.unchecked({x,y,z},0)=value;
      }
      PrimitiveHistory species{as_const(q.view),as_const(old.view),as_const(old.view)};
      EquationStateView state;
      state.density={as_const(rho.view),as_const(rho.view),as_const(rho.view)};
      state.velocity={as_const(velocity.view),as_const(velocity.view),as_const(velocity.view)};
      state.independent_species={&species,1};
      const auto diffusivity=as_const(gamma.view);
      EquationMaterialView material;material.scalar_mass_diffusivity={&diffusivity,1};
      EquationContributionView mixing;
      mixing.explicit_source_density=as_const(source.view);mixing.implicit_sink_density=as_const(sink.view);
      mixing.has_implicit_sink=true;mixing.conserved_quantity=scalar_field;
      mixing.units.si_exponents=heat ? std::array<std::int8_t,7>{1,-1,-3,0,0,0,0} : std::array<std::int8_t,7>{1,-3,-1,0,0,0,0};mixing.stage=9;
      mixing.explicit_source_field=source_field;mixing.implicit_sink_field=sink_field;
      EquationAssemblyContext context;
      context.dt=dt;context.bdf={1/dt,-1/dt,0,1};context.time=1030;
      context.geometry=fixture.geometry.topology_revision();context.boundary=fixture.boundary.revision();
      context.thermo=fixture.thermodynamics.fingerprint();context.transport=fixture.transport.fingerprint();
      context.contribution_stage=9;context.scope=EquationAssemblyScope::final_conservative;
      context.mass_flux=flux;context.face_flux=flux.revision;
      context.face_flux_authority=flux.certificate.authority();context.face_flux_storage=flux.certificate.storage();
      context.face_flux_revision_domain=flux.certificate.revision_domain();
      EquationSystemView system{diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view};
      const auto assemble_equation=[&](EquationAssemblyCertificate& certificate) noexcept {
        if (heat) {
          state.enthalpy=species;
          return detail::StatisticalEnthalpy::assemble(fixture.equations.enthalpy(),state,
              diffusivity,{&mixing,1},context,system,certificate);
        }
        if (frozen) return detail::StatisticalSpecies::assemble(fixture.equations.species(),0,
            state,material,{&mixing,1},context,system,certificate);
        return assemble_species(fixture.equations.species(),0,state,material,{&mixing,1},context,system,certificate);
      };
      const LinearIdentity identity{1040,context.time,fixture.geometry.fingerprint(),workspace.fingerprint(),1041};
      detail::ColdPressureOperator op(rows,cells,halo,identity);
      detail::ColdPressureDilu pc(rows,cells,identity,true);
      for (int sweep=0;sweep<2;++sweep) {
        if(frozen) {
          detail::FrozenScalarProblem problem{fixture.equations.kernels(),fixture.boundary,
              boundary_stage,ConvectionScheme::central2,as_const(velocity.view),
              context,state.density,species,q.view};
          detail::ScalarCorrectionStorage storage{system,variation.view,linear_rhs.view,
              correction.view,backup.view};
          detail::ScalarCorrectionRuntime runtime{rows,pc,halo,workspace,reductions,identity,
              {1e-12,1e-12,400,1,24}};
          int rank{},ranks{};MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
          int fail_phase=-1;
          const auto assemble=[&](EquationAssemblyCertificate& certificate) noexcept {
            if(fail_phase==0 && rank==ranks-1) return Status{StatusCode::numerical_failure,19101};
            return assemble_equation(certificate);
          };
          const auto refresh=[&](FieldView& view) noexcept {
            HaloTicket ticket;
            auto refreshed=state_halo.begin(141,{&view,1},ticket);
            if(refreshed) refreshed=state_halo.finish(ticket,{&view,1});
            if(refreshed && fail_phase==1 && rank==ranks-1) {
              view.unchecked({-2,0,0},0)=73.;
              return Status{StatusCode::numerical_failure,19102};
            }
            return refreshed;
          };
          if(sweep==0) {
            const auto saved=q.bytes,accepted=old.bytes;
            const auto revision=q.view.revision;
            if (heat) {
              const auto units=mixing.units;
              mixing.units.si_exponents={1,-3,-1,0,0,0,0};
              const auto rejected=detail::correct_frozen_scalar(problem,storage,runtime,assemble,refresh);
              mixing.units=units;
              int valid=!rejected.status && rejected.status.code==StatusCode::invalid_plan &&
                  q.bytes==saved && old.bytes==accepted && q.view.revision==revision;
              MPI_Allreduce(MPI_IN_PLACE,&valid,1,MPI_INT,MPI_MIN,comm);
              if (!valid) {std::cerr<<"statistical enthalpy source units admission\n";return false;}
            }
            for(int phase:{0,1}) {
              fail_phase=phase;
              const auto failed=detail::correct_frozen_scalar(problem,storage,runtime,assemble,refresh);
              int okay=!failed.status && q.bytes==saved && old.bytes==accepted &&
                  q.view.revision==revision && species.trial.revision==revision &&
                  failed.lowest_failing_rank==ranks-1 && failed.status.detail==19101U+phase;
              MPI_Allreduce(MPI_IN_PLACE,&okay,1,MPI_INT,MPI_MIN,comm);
              if(!okay) {std::cerr<<"frozen correction rollback phase="<<phase<<'\n';return false;}
            }
            fail_phase=-1;
            runtime.control.maximum_iterations=1;
            const auto exhausted=detail::correct_frozen_scalar(problem,storage,runtime,assemble,refresh);
            runtime.control.maximum_iterations=400;
            int rolled_back=!exhausted.status &&
                exhausted.termination==LinearTermination::maximum_iterations &&
                q.bytes==saved && old.bytes==accepted && q.view.revision==revision &&
                species.trial.revision==revision;
            MPI_Allreduce(MPI_IN_PLACE,&rolled_back,1,MPI_INT,MPI_MIN,comm);
            if(!rolled_back) {std::cerr<<"frozen correction Krylov rollback termination="
                <<unsigned(exhausted.termination)<<'\n';return false;}
            // Aliasing the accepted state is rejected before assembly.
            const auto history=species.accepted;species.accepted=as_const(q.view);
            fail_phase=-1;
            const auto aliased=detail::correct_frozen_scalar(problem,storage,runtime,assemble,refresh);
            species.accepted=history;
            int okay=!aliased.status && q.bytes==saved && old.bytes==accepted;
            MPI_Allreduce(MPI_IN_PLACE,&okay,1,MPI_INT,MPI_MIN,comm);
            if(!okay) return false;
          }
          fail_phase=-1;
          const auto solved=detail::correct_frozen_scalar(problem,storage,runtime,assemble,refresh);
          if(!solved.status) {std::cerr<<"frozen correction status="<<unsigned(solved.status.code)
              <<'/'<<solved.status.detail<<'\n';return false;}
          iterations+=solved.iterations;
          // RHS changes retain the factor; an actual operator change rebuilds it.
          if (!pc.prepare() || !pc.reused_last_prepare()) return false;
          const auto saved_row=rows.front();
          rows.front().rhs+=1.;
          if (!pc.prepare() || !pc.reused_last_prepare()) return false;
          rows.front().diagonal+=1.;
          if (!pc.prepare() || pc.reused_last_prepare()) return false;
          rows.front()=saved_row;
          if (!pc.prepare() || pc.reused_last_prepare()) return false;
        } else {
        EquationAssemblyCertificate certificate;
        status=assemble_equation(certificate);
        if (status) status=detail::close_mixture_scalar_rows(fixture.equations.kernels(),fixture.boundary,
            boundary_stage,scalar_field,ConvectionScheme::central2,as_const(velocity.view),context,
            system,certificate,variation.view,{rows.data(),rows.size()});
        status=reductions.consensus(status);
        if (status) status=pc.prepare();
        status=reductions.consensus(status);
        if (!status) {std::cerr<<"IEM matrix status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
        std::size_t i{};
        for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y) for (int x=0;x<cells.x;++x,++i)
          linear_rhs.view.unchecked({x,y,z},0)=rows[i].rhs;
        fill_field(correction,0.);
        const LinearSolveControl control{1e-12,1e-12,400,1,24};
        const auto solved=solve_fgmres(op,pc,{as_const(linear_rhs.view),correction.view,identity,control},
            workspace,reductions);
        if (!solved.status) {std::cerr<<"IEM solve status="<<unsigned(solved.status.code)<<'/'<<solved.status.detail<<'\n';return false;}
        iterations+=solved.iterations;
        for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y) for (int x=0;x<cells.x;++x)
          q.view.unchecked({x,y,z},0)+=correction.view.unchecked({x,y,z},0);
        ++q.view.revision;species.trial=as_const(q.view);
        HaloTicket ticket;
        status=state_halo.begin(141,{&q.view,1},ticket);
        if (status) status=state_halo.finish(ticket,{&q.view,1});
        if (!status) return false;
        }
        const double a=1+dt*(terms.implicit_sink_density+4*(mu/.7)*std::pow(std::sin(kdx/2),2)/(dx*dx))/rho_value;
        const double b=dt*speed*std::sin(kdx)/dx;
        for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y) for (int x=0;x<cells.x;++x) {
          const double theta=kdx*(x+fixture.patch.begin.x+.5);
          const double expected=mean+sign*amplitude*(frozen
              ? 1/(1+dt*terms.implicit_sink_density/rho_value)
              : (a*std::cos(theta)+b*std::sin(theta))/(a*a+b*b));
          errors[0]=std::max(errors[0],std::abs(q.view.unchecked({x,y,z},0)-expected)/scale);
        }
      }
      EquationAssemblyCertificate certificate;
      status=assemble_equation(certificate);
      status=reductions.consensus(status);if (!status) return false;
      double sum{};
      for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y) for (int x=0;x<cells.x;++x) {
        const Int3 cell{x,y,z};
        errors[1]=std::max(errors[1],std::abs(residual.view.unchecked(cell,0)/detail::cell_volume(fixture.equations.kernels(),cell)-
            (frozen ? q.view.unchecked(cell,0)*mass_divergence.view.unchecked(cell,0) : 0))/scale);
        sum+=q.view.unchecked(cell,0)-mean;
      }
      MPI_Allreduce(MPI_IN_PLACE,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
      const double expected_mean_offset=frozen ? sign*amplitude/(1+dt*terms.implicit_sink_density/rho_value) : 0;
      errors[2]=std::max(errors[2],std::abs(sum/(global_n*global_n*global_n)-expected_mean_offset)/scale);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,errors,3,MPI_DOUBLE,MPI_MAX,comm);
  int rank{};MPI_Comm_rank(comm,&rank);
  if (rank==0)std::cout<<"IEM implicit transport frozen="<<frozen<<" enthalpy="<<heat<<" solution="<<errors[0]<<" residual="<<errors[1]
      <<" mean="<<errors[2]<<" iterations="<<iterations<<'\n';
  return errors[0]<1e-11 && errors[1]<1e-10 && errors[2]<1e-12;
}

bool iem_row_probe(std::istream& input=std::cin,std::ostream& output=std::cout) {
  double volume{},mu{},mut{},cd{},density{},dt{},mean{},old{},previous_volume{};
  std::unique_ptr<ProductionFixture> fixture;
  output << std::setprecision(17);
  while (input >> volume >> mu >> mut >> cd >> density >> dt >> mean >> old) {
    if (!fixture || volume!=previous_volume) {
      fixture=std::make_unique<ProductionFixture>();
      if (!make_production_fixture(4,*fixture,false,-1,true,volume)) {
        std::cerr<<"IEM fixture volume="<<volume<<'\n'; return false;
      }
      previous_volume=volume;
    }
    const auto cells=fixture->patch.cells;
    auto rho=make_field(kDensity,cells,1,2,801);
    auto q=make_field(kSpecies,cells,1,2,802);
    auto accepted=make_field(kSpecies,cells,1,2,803);
    auto gamma=make_field(25,cells,1,2,804);
    auto source=make_field(10,cells,1,0,805);
    auto sink=make_field(11,cells,1,0,806);
    auto diagonal=make_field(30,cells,1,0,807);
    auto rhs=make_field(31,cells,1,0,808);
    auto residual=make_field(32,cells,1,0,809);
    auto ax=make_face_field(CartesianAxis::x,cells,810);
    auto ay=make_face_field(CartesianAxis::y,cells,811);
    auto az=make_face_field(CartesianAxis::z,cells,812);
    fill_field(rho,density); fill_field(q,old); fill_field(accepted,old);
    fill_field(gamma,(mu+mut)/.7);
    esf::detail::IemSource terms;
    if (esf::detail::iem_source(volume,mu,mut,cd,1.,mean,terms)!=portable::Status::success)
      return false;
    fill_field(source,terms.explicit_source_density);
    fill_field(sink,terms.implicit_sink_density);
    FinalFluxFixture owner;
    ConstFaceFluxView flux;
    if (!make_linear_final_flux(fixture->equations.kernels(),cells,owner,flux,true)) {
      std::cerr<<"IEM zero final flux\n"; return false;
    }
    PrimitiveHistory species{as_const(q.view),as_const(accepted.view),as_const(accepted.view)};
    EquationStateView state;
    state.density={as_const(rho.view),as_const(rho.view),as_const(rho.view)};
    state.independent_species={&species,1};
    const auto diffusivity=as_const(gamma.view);
    EquationMaterialView material; material.scalar_mass_diffusivity={&diffusivity,1};
    EquationContributionView mixing;
    mixing.explicit_source_density=as_const(source.view);
    mixing.implicit_sink_density=as_const(sink.view);
    mixing.has_implicit_sink=true; mixing.conserved_quantity=kSpecies;
    mixing.units.si_exponents={1,-3,-1,0,0,0,0}; mixing.stage=9;
    mixing.explicit_source_field=10; mixing.implicit_sink_field=11;
    EquationAssemblyContext context;
    context.dt=dt; context.bdf={1/dt,-1/dt,0,1}; context.time=901;
    context.geometry=fixture->geometry.topology_revision(); context.boundary=fixture->boundary.revision();
    context.thermo=fixture->thermodynamics.fingerprint(); context.transport=fixture->transport.fingerprint();
    context.contribution_stage=9; context.scope=EquationAssemblyScope::final_conservative;
    context.mass_flux=flux; context.face_flux=flux.revision;
    context.face_flux_authority=flux.certificate.authority();
    context.face_flux_storage=flux.certificate.storage();
    context.face_flux_revision_domain=flux.certificate.revision_domain();
    EquationSystemView system{diagonal.view,rhs.view,residual.view,ax.view,ay.view,az.view};
    const Int3 cell{1,1,1};
    double first_a{},first_b{},solution{},final_residual{};
    for (int sweep=0;sweep<2;++sweep) {
      EquationAssemblyCertificate certificate;
      const auto status=assemble_species(fixture->equations.species(),0,state,material,
          {&mixing,1},context,system,certificate);
      if (!status) {
        std::cerr<<"IEM assembly status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';
        return false;
      }
      double neighbor{};
      const std::array<ConstFaceFieldView,3> faces{as_const(ax.view),as_const(ay.view),as_const(az.view)};
      for (unsigned a=0;a<3;++a) {
        auto upper=cell; (a==0 ? upper.x : a==1 ? upper.y : upper.z)++;
        neighbor+=faces[a].unchecked(cell)+faces[a].unchecked(upper);
      }
      const double physical_volume=detail::cell_volume(fixture->equations.kernels(),cell);
      const double a=(diagonal.view.unchecked(cell,0)-neighbor)/physical_volume;
      const double b=(rhs.view.unchecked(cell,0)-neighbor*q.view.unchecked(cell,0))/physical_volume;
      if (sweep==0) {
        first_a=a; first_b=b; solution=b/a;
        // The full original residual includes the implicit sink at the trial
        // state, while the RHS retains the same accepted history on jstep 2.
        if (std::abs(residual.view.unchecked(cell,0)/physical_volume-(a*old-b))>
            1e-11*std::max({1.,std::abs(a),std::abs(b)})) {
          std::cerr<<"IEM row residual="<<residual.view.unchecked(cell,0)/physical_volume
                   <<" a="<<a<<" b="<<b<<" volume="<<physical_volume<<'\n'; return false;
        }
        fill_field(q,solution);
        ++q.view.revision;
        species.trial=as_const(q.view);
      } else {
        if (std::abs(a-first_a)>1e-11*std::max(1.,std::abs(first_a)) ||
            std::abs(b-first_b)>1e-11*std::max(1.,std::abs(first_b))) {
          std::cerr<<"IEM repeated row a="<<a<<" first_a="<<first_a<<" b="<<b<<" first_b="<<first_b<<'\n'; return false;
        }
        final_residual=residual.view.unchecked(cell,0)/physical_volume;
      }
    }
    output<<first_a<<' '<<first_b<<' '<<solution<<' '<<final_residual<<'\n';
  }
  return input.eof();
}

bool test_implicit_mixing_rows() {
  std::istringstream input("1e-12 1e-5 1e-5 2 .1 1e-4 .2 1\n"
                          "1e-9 1.8e-5 1e-4 2 1 .001 .5 .2\n"
                          "1e-6 4e-5 0 0 5 1e-7 .8 0\n");
  std::ostringstream output;
  if (!iem_row_probe(input,output)) return false;
  std::istringstream rows(output.str());
  for (const auto& expected:std::array<std::array<double,2>,3>{{
      {{3000,1400}},{{1118,259}},{{5e7,0}}}}) {
    double a{},b{},q{},residual{};
    if (!(rows>>a>>b>>q>>residual) || !close(a,expected[0]) || !close(b,expected[1]) ||
        !close(q,expected[1]/expected[0]) ||
        std::abs(residual)>1e-11*std::max(1.,std::abs(expected[1]))) return false;
  }
  return expect(true,"registered IEM source, implicit diagonal and accepted history match analytic rows");
}

bool cold_scalar_row_probe() {
  double volume{}, density{}, old_value{}, dt{};
  while (std::cin >> volume >> density >> old_value >> dt) {
    std::array<hundun::v04::detail::ColdTransportFace,6> faces{};
    for (auto& face:faces)
      if (!(std::cin >> face.diffusion >> face.lower_weight >> face.mass_flux))
        return false;
    hundun::v04::detail::ColdPressureRow spatial,row;
    if (!hundun::v04::detail::assemble_cold_transport_row(volume,faces,0.0,spatial) ||
        !hundun::v04::detail::add_cold_backward_euler_storage(
            spatial,old_value,density,dt,row)) return false;
    std::cout << std::setprecision(17);
    for (unsigned axis=0;axis<3;++axis)
      std::cout << row.neighbour[2*axis+1] << ' ' << row.neighbour[2*axis] << ' ';
    std::cout << row.diagonal << ' ' << row.rhs << '\n';
  }
  return std::cin.eof();
}

bool test_density_pressure_rounding() {
  constexpr double rho=.21339443457677126, psi=1.2214966298693773e-6;
  constexpr double pressure=174699.15950532828;
  const double evaluated=rho/psi;
  bool pass=expect(pressure!=evaluated &&
      detail::cold_density_pressure_candidate(rho,psi,pressure,evaluated)==pressure,
      "rounding-compatible pressure iterate survives density EOS refresh");
  const double pressure_gap=std::abs(pressure-evaluated);
  pass &= expect(detail::cold_density_pressure_candidate(rho,psi,pressure,evaluated,
      .5*pressure_gap)==evaluated &&
      detail::cold_density_pressure_candidate(rho,psi,pressure,evaluated,
      2*pressure_gap)==pressure,
      "energy pressure-work budget resolves an EOS-compatible pressure offset");
  // A pressure perturbation update and the absolute-pressure rebase can
  // move the forward EOS by two density ULPs while preserving the existing
  // cold pressure update's 16-epsilon contract.
  constexpr double split_rho=.21339443457689178;
  constexpr double split_psi=1.1166874397185174e-6;
  constexpr double split_pressure=191095.93874422193;
  pass &= expect(detail::cold_density_pressure_candidate(split_rho,split_psi,
      split_pressure,split_rho/split_psi)==split_pressure,
      "pressure refresh retains a split-reference state inside the existing EOS gate");
  pass &= expect(detail::cold_density_pressure_candidate(rho,psi,pressure+1e-8,evaluated)==evaluated &&
      detail::cold_density_pressure_candidate(rho,psi,-1.,evaluated)==evaluated &&
      detail::cold_density_pressure_candidate(rho,psi,std::numeric_limits<double>::quiet_NaN(),evaluated)==evaluated,
      "a physical EOS change or invalid pressure uses the evaluated state");
  return pass;
}

bool test_cold_row_residual_precision() {
  // (1+2^-27)*(1-2^-27)-1 = -2^-54 exactly. Rounding the
  // matrix action before subtracting its RHS instead produces zero.
  const double small=std::ldexp(1.,-27), expected=-std::ldexp(1.,-54);
  detail::ColdPressureRow row;
  row.diagonal=1+small; row.rhs=1;
  std::array<double,6> neighbours{};
  bool pass=expect(detail::cold_row_residual(row,1-small,neighbours)==expected,
      "complete row retains the low product bits in its residual");
  // Cancellation may also occur between spatial neighbour products.
  row.diagonal=1; row.rhs=0; row.neighbour[0]=1+small;
  neighbours[0]=1-small;
  pass &= expect(detail::cold_row_residual(row,1,neighbours)==-expected,
      "neighbour product cancellation retains the exact residual");
  row={};row.diagonal=1;
  pass &= expect(detail::cold_row_residual(row,0,neighbours)==0,
      "isolated solid unit row keeps its exact zero residual");
  return pass;
}

bool test_species_search_quantization() {
  bool factor_passed=true;
  for (const bool pure_major : {false,true}) for (const bool reverse : {false,true}) {
    std::array<double,2> anchor{pure_major ? 1.0 : std::nextafter(1.0,0.0),
                                pure_major ? 5.417e-20 : 1.0e-16};
    std::array<double,2> proposal{anchor[0],pure_major ? 8.862e-20 : 1.5e-16};
    if (reverse) {std::swap(anchor[0],anchor[1]);std::swap(proposal[0],proposal[1]);}
    const double theta=detail::species_search_factor(anchor.size(),
        [&](std::size_t s){return anchor[s];},[&](std::size_t s){return proposal[s];});
    factor_passed &= expect(theta==1.0,
        "nonlinear search retains a full EOS-admissible trace update");
  }
  for (const std::array<double,2> proposal :
       {std::array<double,2>{.9,.3},std::array<double,2>{.8,-.1}}) {
    const std::array<double,2> anchor{.7,.2};
    const double theta=detail::species_search_factor(anchor.size(),
        [&](std::size_t s){return anchor[s];},[&](std::size_t s){return proposal[s];});
    const double a=anchor[0]+theta*(proposal[0]-anchor[0]);
    const double b=anchor[1]+theta*(proposal[1]-anchor[1]);
    factor_passed &= expect(theta>=0 && theta<1 && a>=0 && b>=0 && a+b<=1,
        "nonlinear search still bounds negative and overfull proposals together");
  }
  constexpr double midpoint_q=8.2308775753956881e-9;
  constexpr double midpoint_correction=8.2728821280182545e-25;
  const bool midpoint_pass=expect(
      detail::species_search_update(midpoint_q,midpoint_correction)==std::nextafter(midpoint_q,0.),
      "single FP64 rounding preserves a correction just beyond half an ULP");
  // Captured O2 row at dt=2.5e-6: the nearest represented composition
  // already minimizes this diagonal search, while the stricter inner
  // budget is 8.867724281788395e-13 in residual-density units.
  constexpr double q=.22399976475038841, diagonal=85357.781855068519;
  constexpr double residual=1.1425655268273257e-12;
  const double proposal=static_cast<double>(static_cast<long double>(q)-residual/diagonal);
  bool pass=midpoint_pass;
  pass &= expect(proposal==q && detail::species_search_quantized(q,proposal,diagonal,residual),
      "represented O2 search detects its half-ULP floor");
  pass &= expect(!detail::species_search_quantized(q,q,diagonal,1.3e-12),
      "a stalled search above half a diagonal ULP keeps iterating");
  pass &= expect(!detail::species_search_quantized(q,std::nextafter(q,0.),diagonal,residual),
      "a represented improving direction keeps iterating");
  pass &= expect(!detail::species_search_quantized(0,0,1,1e-320) &&
      !detail::species_search_quantized(1,1,1,-1e-16) &&
      !detail::species_search_quantized(q,q,0,residual) &&
      !detail::species_search_quantized(q,q,diagonal,std::numeric_limits<double>::quiet_NaN()),
      "composition bounds and finite positive search metric remain authoritative");
  // The H2O2 row at dt=1.25e-6 separates search quantization from
  // conservation: R_Y-Y*R_mass is at its floor while R_Y can still improve.
  constexpr double trace=1.6916092542060323e-11;
  constexpr double trace_diagonal=170715.54766318493;
  constexpr double full_residual=-3.9709886867014247e-22;
  constexpr double mass_residual=-9.4153959425404769e-12;
  const double reduced=full_residual-trace*mass_residual;
  const double full_proposal=static_cast<double>(static_cast<long double>(trace)-
      static_cast<long double>(full_residual)/trace_diagonal);
  const double reduced_proposal=static_cast<double>(static_cast<long double>(trace)-
      static_cast<long double>(reduced)/trace_diagonal);
  pass &= expect(detail::species_search_quantized(trace,reduced_proposal,trace_diagonal,reduced) &&
      !detail::species_search_quantized(trace,full_proposal,trace_diagonal,full_residual) &&
      std::abs(full_residual+trace_diagonal*(full_proposal-trace))<std::abs(full_residual),
      "full conservative trace row retains its represented improving direction");
  return pass && factor_passed;
}

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  if (argc==2 && std::string_view(argv[1])=="--bound") {
    const bool passed=test_mixture_bound_repair(MPI_COMM_WORLD);MPI_Finalize();return passed ? 0 : 1;
  }
  if (argc==2 && std::string_view(argv[1])=="--mixture-material") {
    const bool passed=test_mixture_material_floor();MPI_Finalize();return passed ? 0 : 1;
  }
  if (argc==2 && (std::string_view(argv[1])=="--face" || std::string_view(argv[1])=="--flat-face")) {
    const bool passed=mixture_face_probe(std::string_view(argv[1])=="--flat-face");
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc==2 && std::string_view(argv[1])=="--iem-solve") {
    const bool passed=test_implicit_mixture_solve(MPI_COMM_WORLD);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc==2 && std::string_view(argv[1])=="--esf-h-solve") {
    const bool passed=test_implicit_mixture_solve(MPI_COMM_WORLD,true,true);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc==2 && std::string_view(argv[1])=="--frozen-solve") {
    const bool passed=test_implicit_mixture_solve(MPI_COMM_WORLD,true);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc==2 && std::string_view(argv[1])=="--iem") {
    const bool passed=iem_row_probe();
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  if (argc==2 && std::string_view(argv[1])=="--row") {
    const bool passed=cold_scalar_row_probe();
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  bool passed = test_mixture_bound_repair();
  passed &= test_independent_species_closure();
  passed &= test_implicit_mixing_rows();
  passed &= test_species_search_quantization();
  passed &= test_cold_row_residual_precision();
  passed &= test_density_pressure_rounding();
  passed &= test_mixture_face_closure();
  passed &= test_mixture_flux_ledger();
  passed &= test_species_storage_increment();
  passed &= test_mixture_equation_faces();
  passed &= test_mixture_material_floor();
  passed &= test_frozen_scalar_rows();
  passed &= test_passive_midpoint();
  passed &= test_inlet_scalar_material();
  passed &= test_ibm_species_matrix();
  passed &= test_ibm_species_matrix(true);
  passed &= test_ibm_species_matrix(true,true);
  passed &= test_ibm_species_matrix(false,true);
  passed &= test_composition_dependent_production_eos();
  passed &= test_scalar_catalog_contract();
  passed &= test_scalar_mass_diffusivity_oracle_and_atomicity();
  passed &= test_outward_flux_sign_contract();
  passed &= test_conservative_global_balance_oracle();
  passed &= test_production_species_and_passive_assembly();
  passed &= test_reacting_registration_stays_out_of_scope();
  MPI_Finalize();
  return passed ? 0 : 1;
}
