// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/ibm_force_fixture.hpp"
#include "../support/product_fixture.hpp"
#include "../support/turbulence_fixture.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

bool expect(Status status, std::string_view description) {
  if (!status)
    std::cerr << "FAIL: " << description << " status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << '\n';
  return static_cast<bool>(status);
}

struct OwnedFace {
  std::vector<double> values;
  FaceFieldView view{};
};

OwnedFace make_face(CartesianAxis axis, Int3 cells,
                    StorageIdentity identity) {
  Int3 extents = cells;
  if (axis == CartesianAxis::x) ++extents.x;
  if (axis == CartesianAxis::y) ++extents.y;
  if (axis == CartesianAxis::z) ++extents.z;
  OwnedFace result;
  const std::size_t stride_y = static_cast<std::size_t>(extents.x);
  const std::size_t stride_z = stride_y * extents.y;
  result.values.assign(stride_z * extents.z, 7.0);
  result.view = {result.values.data(), extents, stride_y, stride_z, axis,
                 identity, identity + 1000U};
  return result;
}

std::size_t flat(Int3 cells, Int3 cell, std::uint8_t component = 0U) {
  const std::size_t count =
      static_cast<std::size_t>(cells.x) * cells.y * cells.z;
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(cells.y) * cell.z) +
         component * count;
}

FaceFieldView select(FaceFluxView flux, ImmersedFaceDirection direction) {
  if (direction == ImmersedFaceDirection::x_negative ||
      direction == ImmersedFaceDirection::x_positive)
    return flux.x;
  if (direction == ImmersedFaceDirection::y_negative ||
      direction == ImmersedFaceDirection::y_positive)
    return flux.y;
  return flux.z;
}

Int3 face_index(const ImmersedLink& link) {
  Int3 face = link.fluid_local_index;
  if (link.direction == ImmersedFaceDirection::x_positive) ++face.x;
  if (link.direction == ImmersedFaceDirection::y_positive) ++face.y;
  if (link.direction == ImmersedFaceDirection::z_positive) ++face.z;
  return face;
}

std::uint8_t face_axis(ImmersedFaceDirection direction) {
  if (direction == ImmersedFaceDirection::x_negative ||
      direction == ImmersedFaceDirection::x_positive)
    return 0U;
  if (direction == ImmersedFaceDirection::y_negative ||
      direction == ImmersedFaceDirection::y_positive)
    return 1U;
  return 2U;
}

bool positive_face(ImmersedFaceDirection direction) {
  return direction == ImmersedFaceDirection::x_positive ||
         direction == ImmersedFaceDirection::y_positive ||
         direction == ImmersedFaceDirection::z_positive;
}

bool test_prescribed_interface_mass_flux() {
  constexpr std::int32_t n = 16;
  IbmForceFixture fixture;
  bool passed = expect(fixture.initialize(MPI_COMM_SELF, n),
                       "prescribed-flux IBM fixture compiles");
  if (!passed) return false;

  ValidatedModel model = product_model({n, n, n});
  model.mesh = force_mesh(n);
  model.fingerprint = 88002U;
  FieldRegistry registry;
  BoundaryPlan physical_boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  passed &= expect(BoundaryCompiler::compile(
                       MPI_COMM_SELF, model, fixture.geometry, fixture.patch,
                       registry, physical_boundary, schemes, time),
                   "prescribed-flux physical boundary compiles");
  CartesianKernelPlan kernels;
  passed &= expect(CartesianKernelPlan::compile(
                       schemes, fixture.geometry, fixture.patch,
                       physical_boundary, kernels),
                   "prescribed-flux Cartesian kernel compiles");
  if (!passed) return false;

  const Span<const ImmersedLink> links = fixture.topology.links();
  std::size_t source_index = links.size;
  for (std::size_t index = 0U; index < links.size; ++index)
    if (links.data[index].direction == ImmersedFaceDirection::y_negative) {
      source_index = index;
      break;
    }
  passed &= expect(source_index < links.size && links.size > 1U,
                   "fixture exposes a GTMC-oriented source and a sealed link");
  if (!passed) return false;
  const std::size_t sealed_index = source_index == 0U ? 1U : 0U;

  constexpr double prescribed_phi = 2.7866666668e-6;
  const std::array<double, 2U> prescribed_species{{0.875, 0.125}};
  const std::array<IbmInterfaceInletState, 1U> sources{{
      {links.data[source_index].global_link, prescribed_phi,
       {1.25, 2.5, 3.75}, 450000.0,
       {prescribed_species.data(), prescribed_species.size()}},
  }};
  IbmEquationInterfacePlan interface;
  passed &= expect(IbmEquationInterfacePlan::compile(
                       kernels, fixture.topology, fixture.boundary,
                       fixture.topology.interface_metric(),
                       {sources.data(), sources.size()},
                       prescribed_species.size(), interface),
                   "prescribed interface mass-flux source compiles");
  if (!passed) return false;

  const Int3 cells = fixture.patch.cells;
  OwnedFace x_flux = make_face(CartesianAxis::x, cells, 211U);
  OwnedFace y_flux = make_face(CartesianAxis::y, cells, 211U);
  OwnedFace z_flux = make_face(CartesianAxis::z, cells, 211U);
  FaceFluxView flux{x_flux.view, y_flux.view, z_flux.view, 22U, {}};
  passed &= expect(interface.constrain_interface_flux(flux),
                   "prescribed and sealed IBM flux authority applies");
  for (std::size_t index = 0U; index < links.size; ++index) {
    const ImmersedLink& link = links.data[index];
    const double value =
        select(flux, link.direction).unchecked(face_index(link));
    if (index == source_index)
      passed &= expect(value == prescribed_phi,
                       "source link carries the exact prescribed phi");
    else
      passed &= expect(value == 0.0 && !std::signbit(value),
                       "sealed link carries canonical positive zero phi");
  }
  passed &= expect(interface.validate_interface_flux(as_const(flux)),
                   "prescribed and sealed interface flux validates");
  double queried_phi = -1.0;
  const ImmersedLink& source_link = links.data[source_index];
  const CartesianAxis source_axis =
      static_cast<CartesianAxis>(face_axis(source_link.direction));
  ForceOwnedField inlet_velocity = make_force_field(
      41U, cells, 3U, fixture.boundary.maximum_halo_reach(), 521U, 621U);
  ForceOwnedField inlet_gradient =
      make_force_field(42U, cells, 9U, 0U, 522U, 622U);
  const int reach = fixture.boundary.maximum_halo_reach();
  for (int z = -reach; z < cells.z + reach; ++z)
    for (int y = -reach; y < cells.y + reach; ++y)
      for (int x = -reach; x < cells.x + reach; ++x) {
        inlet_velocity.view.unchecked({x, y, z}, 0U) = 1.25;
        inlet_velocity.view.unchecked({x, y, z}, 1U) = 2.5;
        inlet_velocity.view.unchecked({x, y, z}, 2U) = 3.75;
      }
  passed &= expect(interface.correct_velocity_gradient(
                       as_const(inlet_velocity.view), inlet_gradient.view),
                   "prescribed inlet velocity gradient evaluates");
  passed &= expect(std::abs(inlet_gradient.view.unchecked(
                       source_link.fluid_local_index, 4U)) < 1.0e-10,
                   "matching constant inlet velocity is not a no-slip wall");
  ForceOwnedField inlet_mu = make_force_field(
      43U, cells, 1U, fixture.boundary.maximum_halo_reach(), 523U, 623U);
  constexpr double mu_inlet = 0.017;
  std::fill(inlet_mu.storage.begin(), inlet_mu.storage.end(), mu_inlet);
  double inlet_work = 123.0;
  passed &= expect(interface.inlet_viscous_work_input(
                       as_const(inlet_velocity.view), as_const(inlet_mu.view),
                       inlet_work) && std::abs(inlet_work) < 1.0e-10,
                   "constant matching inlet U produces zero viscous work");
  constexpr double inlet_slope = 2.0;
  for (int z = -reach; z < cells.z + reach; ++z)
    for (int y = -reach; y < cells.y + reach; ++y)
      for (int x = -reach; x < cells.x + reach; ++x) {
        const double centre = fixture.extrapolated_centre(
            fixture.geometry.y(), fixture.patch.begin.y + y);
        inlet_velocity.view.unchecked({x, y, z}, 1U) =
            2.5 + inlet_slope * (centre - source_link.wall_point.y);
      }
  const double expected_work = -2.5 * mu_inlet * (4.0 / 3.0) * inlet_slope *
      fixture.topology.interface_metric().links().data[source_index]
          .physical_quadrature_area;
  passed &= expect(interface.inlet_viscous_work_input(
                       as_const(inlet_velocity.view), as_const(inlet_mu.view),
                       inlet_work) && std::abs(inlet_work - expected_work) < 1.0e-10,
                   "linear normal inlet stress has matching signed energy work");
  ForceOwnedField inlet_temperature = make_force_field(
      44U, cells, 1U, fixture.boundary.maximum_halo_reach(), 524U, 624U);
  ForceOwnedField inlet_heat_rate =
      make_force_field(45U, cells, 1U, 0U, 525U, 625U);
  std::fill(inlet_temperature.storage.begin(), inlet_temperature.storage.end(), 300.0);
  const auto regions = fixture.topology.region();
  for (int z = 0; z < cells.z; ++z)
    for (int y = 0; y < cells.y; ++y)
      for (int x = 0; x < cells.x; ++x)
        if (regions.data[flat(cells, {x, y, z})] ==
            static_cast<std::uint8_t>(RegionFlag::solid))
          inlet_temperature.view.unchecked({x, y, z}, 0U) = 900.0;
  const std::array<ConstFieldView, 1U> thermal_reads{as_const(inlet_temperature.view)};
  const std::array<FieldView, 1U> thermal_writes{inlet_heat_rate.view};
  const auto reset_thermal = [&]() {
    return cartesian_diffusion(kernels, as_const(inlet_mu.view),
        {{thermal_reads.data(), thermal_reads.size()},
         {thermal_writes.data(), thermal_writes.size()},
         {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U, nullptr});
  };
  passed &= expect(reset_thermal() && interface.correct_zero_normal_diffusion(
                       as_const(inlet_temperature.view), as_const(inlet_mu.view),
                       inlet_heat_rate.view) &&
                       std::abs(inlet_heat_rate.view.unchecked(
                           source_link.fluid_local_index, 0U)) < 1.0e-9,
                   "inlet thermal flux excludes arbitrary hot solid placeholders");
  passed &= expect(reset_thermal() &&
                       interface.correct_positive_bounded_zero_normal_diffusion(
                           as_const(inlet_temperature.view), as_const(inlet_mu.view),
                           inlet_heat_rate.view) &&
                       std::abs(inlet_heat_rate.view.unchecked(
                           source_link.fluid_local_index, 0U)) < 1.0e-9,
                   "bounded inlet thermal closure uses the same zero diffusive flux");
  passed &= expect(interface.prescribed_face_flux(
                       source_axis, face_index(source_link), queried_phi) &&
                       queried_phi == prescribed_phi,
                   "source-face lookup returns the exact prescribed phi");
  double sealed_phi = -1.0;
  passed &= expect(!interface.prescribed_face_flux(
                       static_cast<CartesianAxis>(
                           face_axis(links.data[sealed_index].direction)),
                       face_index(links.data[sealed_index]), sealed_phi) &&
                       sealed_phi == -1.0,
                   "sealed face is absent from prescribed-flux authority");

  select(flux, links.data[source_index].direction)
      .unchecked(face_index(links.data[source_index])) = 0.0;
  passed &= expect(interface.validate_interface_flux(as_const(flux)).code ==
                       StatusCode::numerical_failure,
                   "missing prescribed source phi is rejected");
  passed &= expect(interface.constrain_interface_flux(flux),
                   "interface flux authority restores prescribed phi");
  select(flux, links.data[sealed_index].direction)
      .unchecked(face_index(links.data[sealed_index])) = 1.0e-7;
  passed &= expect(interface.validate_interface_flux(as_const(flux)).code ==
                       StatusCode::numerical_failure,
                   "nonzero sealed-link phi is rejected");

  OwnedFace x_values = make_face(CartesianAxis::x, cells, 212U);
  OwnedFace y_values = make_face(CartesianAxis::y, cells, 212U);
  OwnedFace z_values = make_face(CartesianAxis::z, cells, 212U);
  const FrozenConvectionFaceOutput face_values{
      x_values.view, y_values.view, z_values.view};
  const auto source_value = [&]() {
    const ImmersedLink& link = links.data[source_index];
    return select(FaceFluxView{face_values.x, face_values.y, face_values.z,
                               23U, {}},
                  link.direction)
        .unchecked(face_index(link));
  };
  const auto sealed_value = [&]() {
    const ImmersedLink& link = links.data[sealed_index];
    return select(FaceFluxView{face_values.x, face_values.y, face_values.z,
                               23U, {}},
                  link.direction)
        .unchecked(face_index(link));
  };
  passed &= expect(interface.override_source_face_values(
                       {IbmInterfaceInletFieldKind::velocity, 1U},
                       IbmInterfaceInletEvaluation::value, face_values) &&
                       source_value() == 2.5 && sealed_value() == 7.0,
                   "one inlet state overrides velocity on source faces only");
  passed &= expect(interface.override_source_face_values(
                       {IbmInterfaceInletFieldKind::enthalpy, 0U},
                       IbmInterfaceInletEvaluation::value, face_values) &&
                       source_value() == 450000.0 && sealed_value() == 7.0,
                   "the same inlet state overrides enthalpy");
  passed &= expect(interface.override_source_face_values(
                       {IbmInterfaceInletFieldKind::independent_species, 0U},
                       IbmInterfaceInletEvaluation::value, face_values) &&
                       source_value() == prescribed_species[0U] &&
                       sealed_value() == 7.0,
                   "the same inlet state overrides independent species");
  passed &= expect(interface.override_source_face_values(
                       {IbmInterfaceInletFieldKind::enthalpy, 0U},
                       IbmInterfaceInletEvaluation::fixed_state_variation,
                       face_values) &&
                       source_value() == 0.0 &&
                       !std::signbit(source_value()) &&
                       sealed_value() == 7.0,
                   "fixed inlet state has canonical zero face variation");

  const double prescribed_kinetic =
      0.5 * (1.25 * 1.25 + 2.5 * 2.5 + 3.75 * 3.75);
  passed &= expect(interface.override_source_face_values(
                       {IbmInterfaceInletFieldKind::kinetic_energy, 0U},
                       IbmInterfaceInletEvaluation::value, face_values) &&
                       source_value() == prescribed_kinetic,
                   "kinetic-energy inlet value is derived from the same U");

  ForceOwnedField transported =
      make_force_field(31U, cells, 1U, 2U, 511U, 611U);
  ForceOwnedField correction =
      make_force_field(32U, cells, 1U, 0U, 512U, 612U);
  for (std::int32_t z = -2; z < cells.z + 2; ++z)
    for (std::int32_t y = -2; y < cells.y + 2; ++y)
      for (std::int32_t x = -2; x < cells.x + 2; ++x)
        transported.view.unchecked({x, y, z}, 0U) =
            10.0 + 0.25 * x + 0.5 * y + 0.75 * z;
  const double volume = std::pow(3.0 / static_cast<double>(n), 3.0);
  const double divergence_sign = positive_face(source_link.direction)
                                     ? 1.0
                                     : -1.0;
  const std::array schemes_to_test{
      ConvectionScheme::central2, ConvectionScheme::limited_central2,
      ConvectionScheme::tvd2};
  for (ConvectionScheme scheme : schemes_to_test) {
    std::fill(correction.storage.begin(), correction.storage.end(), 7.0);
    double ordinary_face = 0.0;
    passed &= expect(reconstruct_cartesian_convection_face(
                         kernels, scheme, as_const(transported.view), 0U,
                         source_axis, face_index(source_link), prescribed_phi,
                         ordinary_face),
                     "single-face reconstruction uses production arithmetic");
    passed &= expect(interface.add_source_convection_correction(
                         {IbmInterfaceInletFieldKind::enthalpy, 0U}, scheme,
                         as_const(transported.view), 1.0, correction.view),
                     "prescribed inlet corrects high-order convection");
    const double expected =
        7.0 + divergence_sign * prescribed_phi *
                  (450000.0 - ordinary_face) / volume;
    passed &= expect(correction.view.unchecked(source_link.fluid_local_index,
                                               0U) == expected,
                     "high-order correction replaces exactly one face value");
  }
  std::fill(correction.storage.begin(), correction.storage.end(), 7.0);
  const double ordinary_upwind =
      transported.view.unchecked(source_link.solid_local_index, 0U);
  passed &= expect(interface.add_source_first_order_upwind_correction(
                       {IbmInterfaceInletFieldKind::enthalpy, 0U},
                       as_const(transported.view), 1.0, correction.view),
                   "prescribed inlet corrects first-order upwind convection");
  const double expected_upwind =
      7.0 + divergence_sign * prescribed_phi *
                (450000.0 - ordinary_upwind) / volume;
  passed &= expect(correction.view.unchecked(source_link.fluid_local_index,
                                             0U) == expected_upwind,
                   "upwind correction replaces the solid placeholder donor");

  // The nonlinear enthalpy residual uses the prescribed h_in on the source
  // face, while every E_h action must see the exact derivative +0 there.  The
  // same fixed-face authority must survive the generic and compiled paths.
  passed &= expect(interface.constrain_interface_flux(flux),
                   "source-aware frozen face starts from authoritative phi");
  OwnedFace frozen_x = make_face(CartesianAxis::x, cells, 213U);
  OwnedFace frozen_y = make_face(CartesianAxis::y, cells, 213U);
  OwnedFace frozen_z = make_face(CartesianAxis::z, cells, 213U);
  const FrozenConvectionFaceOutput frozen_output{
      frozen_x.view, frozen_y.view, frozen_z.view};
  OwnedFace direction_x = make_face(CartesianAxis::x, cells, 214U);
  OwnedFace direction_y = make_face(CartesianAxis::y, cells, 214U);
  OwnedFace direction_z = make_face(CartesianAxis::z, cells, 214U);
  const FrozenConvectionFaceOutput directional_output{
      direction_x.view, direction_y.view, direction_z.view};
  const FrozenConvectionContext frozen_context{91001U, 91002U};
  FrozenConvectionFaceField frozen;
  FrozenConvectionFaceDirectionalDerivative ordinary_derivative;
  passed &= expect(freeze_cartesian_target_convection_faces(
                       kernels, ConvectionScheme::limited_central2,
                       as_const(flux), as_const(transported.view), 0U,
                       frozen_context, frozen_output, frozen) &&
                       interface.override_source_face_values(
                           {IbmInterfaceInletFieldKind::enthalpy, 0U},
                           IbmInterfaceInletEvaluation::value, frozen_output),
                   "legacy freeze can be followed by an unsafe source overwrite");
  passed &= expect(differentiate_frozen_cartesian_target_convection_faces(
                       kernels, ConvectionScheme::limited_central2,
                       as_const(flux), as_const(transported.view), 0U,
                       frozen_context,
                       FrozenConvectionLinearizationPolicy::
                           semismooth_generalized_zero_slope,
                       frozen, as_const(transported.view), 0U,
                       directional_output, ordinary_derivative)
                           .code == StatusCode::invalid_plan,
                   "ordinary certificate rejects overwritten source h as stale");
  const auto frozen_source_value = [&]() {
    return select(FaceFluxView{frozen_output.x, frozen_output.y,
                               frozen_output.z, 24U, {}},
                  source_link.direction)
        .unchecked(face_index(source_link));
  };
  const auto directional_source_value = [&]() {
    return select(FaceFluxView{directional_output.x, directional_output.y,
                               directional_output.z, 25U, {}},
                  source_link.direction)
        .unchecked(face_index(source_link));
  };
  for (ConvectionScheme scheme : schemes_to_test) {
    FrozenConvectionFaceDirectionalDerivative derivative;
    passed &= expect(interface.freeze_source_convection_faces(
                         {IbmInterfaceInletFieldKind::enthalpy, 0U}, scheme,
                         as_const(flux), as_const(transported.view), 0U,
                         frozen_context, frozen_output, frozen) &&
                         frozen.fixed_faces.size == 1U &&
                         frozen.fixed_face_authority ==
                             interface.fingerprint() &&
                         frozen_source_value() == 450000.0 &&
                         interface.validate_frozen_source_face_values(
                             {IbmInterfaceInletFieldKind::enthalpy, 0U},
                             frozen) &&
                         differentiate_frozen_cartesian_target_convection_faces(
                             kernels, scheme, as_const(flux),
                             as_const(transported.view), 0U, frozen_context,
                             FrozenConvectionLinearizationPolicy::
                                 semismooth_generalized_zero_slope,
                             frozen, as_const(transported.view), 0U,
                             directional_output, derivative) &&
                         directional_source_value() == 0.0 &&
                         !std::signbit(directional_source_value()),
                     "all generic E_h schemes freeze h_in and return fixed +0");
  }
  passed &= expect(interface.freeze_source_convection_faces(
                       {IbmInterfaceInletFieldKind::enthalpy, 0U},
                       ConvectionScheme::limited_central2, as_const(flux),
                       as_const(transported.view), 0U, frozen_context,
                       frozen_output, frozen),
                   "limited E_h source face refreezes for compiled replay");
  const std::size_t branch_count = frozen_x.values.size() +
                                   frozen_y.values.size() +
                                   frozen_z.values.size();
  std::vector<std::uint16_t> branch_storage(branch_count);
  FrozenConvectionBranchPlan branches;
  passed &= expect(compile_frozen_limited_central2_branches(
                       kernels, as_const(flux), as_const(transported.view), 0U,
                       frozen_context,
                       FrozenConvectionLinearizationPolicy::
                           semismooth_generalized_zero_slope,
                       frozen,
                       {{branch_storage.data(), branch_storage.size()}},
                       branches) &&
                       apply_frozen_limited_central2_branches(
                           kernels, branches, as_const(transported.view), 0U,
                           directional_output) &&
                       directional_source_value() == 0.0 &&
                       !std::signbit(directional_source_value()),
                   "compiled E_h replays the same canonical fixed-source +0");
  const double certified_source_value = frozen_source_value();
  select(FaceFluxView{frozen_output.x, frozen_output.y, frozen_output.z, 24U,
                      {}},
         source_link.direction)
      .unchecked(face_index(source_link)) = certified_source_value + 1.0;
  passed &= expect(interface.validate_frozen_source_face_values(
                       {IbmInterfaceInletFieldKind::enthalpy, 0U}, frozen)
                           .code == StatusCode::invalid_plan,
                   "tampered frozen h_in is rejected by source authority");
  select(FaceFluxView{frozen_output.x, frozen_output.y, frozen_output.z, 24U,
                      {}},
         source_link.direction)
      .unchecked(face_index(source_link)) = certified_source_value;

  auto reversed_sources = sources;
  reversed_sources[0].velocity.y = -2.5;
  const auto prior_fingerprint = interface.fingerprint();
  passed &= expect(IbmEquationInterfacePlan::compile(
                       kernels, fixture.topology, fixture.boundary,
                       fixture.topology.interface_metric(),
                       {reversed_sources.data(), reversed_sources.size()},
                       prescribed_species.size(), interface).code ==
                       StatusCode::invalid_plan &&
                       interface.fingerprint() == prior_fingerprint,
                   "inward mass flux with outward inlet U is rejected atomically");
  {
    ThermodynamicsPlan thermo;
    TransportPlan transport;
    ContributionRegistry contributions;
    EquationPlanSet equations;
    EquationPlanSpec spec;
    spec.density = 60U;
    spec.velocity = physical_boundary.velocity_field();
    spec.pressure_perturbation = physical_boundary.pressure_field();
    spec.enthalpy = physical_boundary.enthalpy_field();
    spec.temperature = 59U;
    spec.effective_viscosity = 61U;
    spec.pressure_compressibility = 62U;
    spec.velocity_gradient = 63U;
    spec.pressure_reference = model.pressure_reference;
    spec.closed_mass_service_stage = 1U;
    spec.maximum_cells_per_rank = static_cast<std::size_t>(cells.x) * cells.y * cells.z;
    const std::array<FieldId, 8U> fields{spec.density, spec.velocity,
        spec.pressure_perturbation, spec.enthalpy, spec.temperature,
        spec.effective_viscosity, spec.pressure_compressibility, spec.velocity_gradient};
    passed &= expect(ThermodynamicsPlan::compile(model.thermophysics, {}, thermo) &&
        TransportPlan::compile(model.thermophysics, thermo, transport) &&
        contributions.configure({fields.data(), fields.size()}) && contributions.freeze() &&
        EquationPlanSet::compile(MPI_COMM_SELF, schemes, fixture.geometry, fixture.patch,
            physical_boundary, contributions, thermo, transport, spec, equations),
        "source-bearing momentum equation plans compile");
    if (!passed) return false;
    auto u = make_force_field(spec.velocity, cells, 3U, reach, 601U, 701U);
    auto rho = make_force_field(spec.density, cells, 1U, reach, 602U, 702U);
    auto p = make_force_field(spec.pressure_perturbation, cells, 1U, reach, 603U, 703U);
    auto mu = make_force_field(spec.effective_viscosity, cells, 1U, reach, 604U, 704U);
    auto grad = make_force_field(spec.velocity_gradient, cells, 9U, 1U, 605U, 705U);
    auto diagonal = make_force_field(70U, cells, 3U, 0U, 606U, 706U);
    auto rhs = make_force_field(71U, cells, 3U, 0U, 607U, 707U);
    auto residual = make_force_field(72U, cells, 3U, 0U, 608U, 708U);
    auto delta = make_force_field(73U, cells, 3U, 0U, 609U, 709U);
    std::fill(rho.storage.begin(), rho.storage.end(), 1.0);
    std::fill(mu.storage.begin(), mu.storage.end(), 0.017);
    std::fill(x_flux.values.begin(), x_flux.values.end(), 0.0);
    std::fill(y_flux.values.begin(), y_flux.values.end(), 0.0);
    std::fill(z_flux.values.begin(), z_flux.values.end(), 0.0);
    passed &= expect(interface.constrain_interface_flux(flux), "source-only momentum phi binds");
    auto ax = make_face(CartesianAxis::x, cells, 801U);
    auto ay = make_face(CartesianAxis::y, cells, 802U);
    auto az = make_face(CartesianAxis::z, cells, 803U);
    EquationStateView state;
    state.velocity = {as_const(u.view), as_const(u.view), as_const(u.view)};
    state.density = {as_const(rho.view), as_const(rho.view), as_const(rho.view)};
    state.pressure_perturbation = {as_const(p.view), as_const(p.view), as_const(p.view)};
    EquationMaterialView material;
    material.effective_viscosity = material.molecular_viscosity = as_const(mu.view);
    EquationAssemblyContext context;
    context.dt = 0.001;
    context.bdf = {1000.0, -1000.0, 0.0, 1U};
    context.time = 901U;
    context.geometry = fixture.geometry.topology_revision();
    context.boundary = physical_boundary.revision();
    context.transport = transport.fingerprint();
    context.face_flux = flux.revision;
    context.contribution_stage = 1U;
    context.scope = EquationAssemblyScope::momentum_predictor;
    context.mass_flux = as_const(flux);
    context.provisional_mass_flux = true;
    context.immersed_interface = &interface;
    EquationSystemView system;
    system.diagonal = diagonal.view;
    system.rhs = rhs.view;
    system.residual = residual.view;
    system.x_coefficient = ax.view;
    system.y_coefficient = ay.view;
    system.z_coefficient = az.view;
    EquationAssemblyCertificate certificate;
    select(flux, source_link.direction).unchecked(face_index(source_link)) = 0.0;
    passed &= expect(assemble_momentum_predictor(equations.momentum(), state, material,
        as_const(grad.view), {}, context, system, delta.view, certificate).code ==
        StatusCode::numerical_failure,
        "source-bearing momentum rejects a missing prescribed mass flux");
    passed &= expect(interface.constrain_interface_flux(flux), "restore momentum source phi");
    passed &= expect(assemble_momentum_predictor(equations.momentum(), state, material,
        as_const(grad.view), {}, context, system, delta.view, certificate),
        "full source-bearing momentum predictor assembles");
    passed &= expect(certificate.inlet_sources == interface.fingerprint(),
                     "momentum certificate binds the source-state authority");
    for (unsigned c = 0U; c < 3U; ++c)
      passed &= expect(std::abs(delta.view.unchecked(source_link.fluid_local_index, c)) < 1e-18,
                       "fixed inlet has no high-minus-low momentum antidiffusion");
    std::array<std::vector<std::uint8_t>, 3U> active_faces;
    std::array<OwnedFace, 12U> limiter_faces;
    for (unsigned group = 0U; group < 4U; ++group)
      for (unsigned axis = 0U; axis < 3U; ++axis)
        limiter_faces[3U * group + axis] =
            make_face(static_cast<CartesianAxis>(axis), cells, 1001U + group);
    const auto is_fluid = [&](Int3 q) {
      return q.x >= 0 && q.y >= 0 && q.z >= 0 &&
          q.x < cells.x && q.y < cells.y && q.z < cells.z &&
          regions.data[flat(cells, q)] == static_cast<std::uint8_t>(RegionFlag::fluid);
    };
    for (unsigned axis = 0U; axis < 3U; ++axis) {
      const auto shape = limiter_faces[axis].view.extents;
      auto& active = active_faces[axis];
      active.resize(static_cast<std::size_t>(shape.x) * shape.y * shape.z);
      for (int z = 0; z < shape.z; ++z)
        for (int y = 0; y < shape.y; ++y)
          for (int x = 0; x < shape.x; ++x) {
            const Int3 face{x, y, z};
            Int3 left = face;
            if (axis == 0U) --left.x;
            else if (axis == 1U) --left.y;
            else --left.z;
            active[flat(shape, face)] = is_fluid(left) && is_fluid(face);
          }
    }
    const MgDomainActivityView activity{regions,
        {active_faces[0].data(), active_faces[0].size()},
        {active_faces[1].data(), active_faces[1].size()},
        {active_faces[2].data(), active_faces[2].size()},
        fixture.topology.fingerprint(), fixture.topology.fingerprint()};
    auto ratios = make_force_field(74U, cells, 3U, 1U, 610U, 710U);
    MomentumPredictorLimiterWorkspace workspace;
    workspace.cell_ratios = ratios.view;
    for (unsigned c = 0U; c < 3U; ++c)
      workspace.high_order_faces[c] = {limiter_faces[3U*c].view,
          limiter_faces[3U*c+1U].view, limiter_faces[3U*c+2U].view, 1200U+c, {}};
    workspace.common_face_alpha = {limiter_faces[9].view, limiter_faces[10].view,
                                   limiter_faces[11].view, 1203U, {}};
    const std::array<HaloFieldSpec, 1U> halo_fields{{{ratios.view.field, 1U, 3U}}};
    HaloEngine halo;
    ReductionEngine reductions;
    passed &= expect(halo.reserve(MPI_COMM_SELF, fixture.patch,
        {halo_fields.data(), halo_fields.size()}, physical_boundary.halo_topology()) &&
        ReductionEngine::compile(MPI_COMM_SELF, ReductionMode::mpi_allreduce, 8U, reductions),
        "source-bearing limiter communication resources compile");
    MomentumPredictorLimiterReport report;
    const auto limit = [&]() {
      return limit_momentum_predictor_correction(MPI_COMM_SELF, equations.momentum(),
          physical_boundary, fixture.patch, certificate, as_const(u.view), as_const(rho.view),
          context.dt, 0.3, as_const(flux), activity, system, workspace, halo, reductions, report);
    };
    passed &= expect(limit().code == StatusCode::invalid_plan,
                     "source-bearing limiter rejects omitted inlet authority");
    workspace.immersed_interface = &interface;
    passed &= expect(limit(), "source-bearing limiter accepts fixed cut-face inflow");
    passed &= expect(report.advective_cfl.valid() &&
        std::abs(report.advective_cfl.absolute_max - context.dt * prescribed_phi / (2.0 * volume)) < 1e-15,
        "source inflow contributes exactly to absolute CFL without opening MG face");
  }
  return passed;
}

bool run() {
  constexpr std::int32_t n = 16;
  IbmForceFixture fixture;
  bool passed = expect(fixture.initialize(MPI_COMM_SELF, n),
                       "equation-interface IBM fixture compiles");
  if (!passed) return false;

  ValidatedModel model = product_model({n, n, n});
  model.mesh = force_mesh(n);
  model.fingerprint = 88001U;
  FieldRegistry registry;
  BoundaryPlan physical_boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  passed &= expect(BoundaryCompiler::compile(
                       MPI_COMM_SELF, model, fixture.geometry, fixture.patch,
                       registry, physical_boundary, schemes, time),
                   "physical boundary and schemes compile");
  CartesianKernelPlan kernels;
  passed &= expect(CartesianKernelPlan::compile(
                       schemes, fixture.geometry, fixture.patch,
                       physical_boundary, kernels),
                   "Cartesian metric kernel compiles");
  IbmEquationInterfacePlan interface;
  passed &= expect(IbmEquationInterfacePlan::compile(
                       kernels, fixture.topology, fixture.boundary,
                       fixture.topology.interface_metric(), interface),
                   "immutable IBM equation replacement compiles");
  const PlanFingerprint bound_interface = interface.fingerprint();
  IbmInterfaceMetricPlan unbound_metric;
  const Status unbound = IbmEquationInterfacePlan::compile(
      kernels, fixture.topology, fixture.boundary, unbound_metric, interface);
  passed &= expect(unbound.code == StatusCode::invalid_plan &&
                       interface.fingerprint() == bound_interface,
                   "equation compile requires metric binding atomically");
  if (!passed) return false;

  const Int3 cells = fixture.patch.cells;
  const auto interface_cells = fixture.topology.interface_cells();
  std::vector<std::uint8_t> interface_membership(
      static_cast<std::size_t>(cells.x) * cells.y * cells.z, 0U);
  for (std::size_t i = 0; i < interface_cells.size; ++i) {
    const auto cell = interface_cells.data[i];
    if (!expect(cell < interface_membership.size(), "interface cell index is in range")) return false;
    passed &= expect(interface_membership[cell]++ == 0U,
                     "interface cell list has no duplicate owners");
  }
  const auto all_links = fixture.topology.links();
  for (std::size_t i = 0; i < all_links.size; ++i)
    passed &= expect(interface_membership[flat(cells, all_links.data[i].fluid_local_index)] == 1U,
                     "interface list covers every cut-face owner");
  const auto interface_only = [&](const ForceOwnedField& value) {
    for (std::size_t i = 0; i < interface_membership.size(); ++i)
      if (interface_membership[i] == 0U && value.storage[i] != 0.0) return false;
    return true;
  };
  // Unbound public objects must reject before touching borrowed plans or data.
  // Also exercise a fresh destination after a failed compile (not a previously
  // valid plan, which the transactional compile deliberately preserves).
  IbmEquationInterfacePlan empty, failed_compile;
  passed &= expect(IbmEquationInterfacePlan::compile(
      kernels, fixture.topology, fixture.boundary, unbound_metric,
      failed_compile).code == StatusCode::invalid_plan,
      "fresh compile failure leaves an unbound object");
  auto sentinel = make_force_field(90U, cells, 3U, 0U, 1U, 901U);
  std::fill(sentinel.storage.begin(), sentinel.storage.end(), 7.25);
  const auto untouched = sentinel.storage;
  for (const auto* invalid : {&empty, &failed_compile}) {
    const auto rejects = [&](Status s) {
      return expect(s.code == StatusCode::invalid_plan &&
                        sentinel.storage == untouched,
                    "unbound IBM public call rejects without writing output");
    };
    passed &= rejects(invalid->zero_interface_flux({}));
    passed &= rejects(invalid->validate_interface_flux({}, 0.0));
    passed &= rejects(invalid->constrain_pressure_predictor(sentinel.view, {}));
    passed &= rejects(invalid->constrain_corrected_state(sentinel.view, {}));
    passed &= rejects(invalid->constrain_momentum({}, {}, {}, {}, {}, {}, nullptr,
                                                {sentinel.view, sentinel.view, sentinel.view}));
    passed &= rejects(invalid->correct_pressure_gradient({}, sentinel.view));
    passed &= rejects(invalid->correct_pressure_work({}, {}, sentinel.view));
    passed &= rejects(invalid->correct_velocity_gradient({}, sentinel.view, {}));
    passed &= rejects(invalid->correct_zero_normal_diffusion({}, {}, sentinel.view));
    passed &= rejects(invalid->correct_impermeable_scalar_diffusion({}, {}, sentinel.view));
    passed &= rejects(invalid->correct_positive_bounded_zero_normal_diffusion({}, {}, sentinel.view));
  }
  // Different base addresses can still overlap through strides/ghost storage.
  // Both transported/rate and diffusivity/rate overlap must reject pre-write.
  auto shared = make_force_field(91U, cells, 1U,
      fixture.boundary.maximum_halo_reach(), 1U, 902U);
  auto other = make_force_field(92U, cells, 1U,
      fixture.boundary.maximum_halo_reach(), 1U, 903U);
  std::fill(shared.storage.begin(), shared.storage.end(), 2.0);
  std::fill(other.storage.begin(), other.storage.end(), 3.0);
  const auto original_shared = shared.storage;
  FieldView offset = shared.view;
  ++offset.base;
  offset.ghosts = {0U, 0U, 0U};
  using Diffusion = Status (IbmEquationInterfacePlan::*)(
      ConstFieldView, ConstFieldView, FieldView) const noexcept;
  for (Diffusion operation : {&IbmEquationInterfacePlan::correct_zero_normal_diffusion,
                             &IbmEquationInterfacePlan::correct_impermeable_scalar_diffusion,
                             &IbmEquationInterfacePlan::correct_positive_bounded_zero_normal_diffusion}) {
    for (bool source_overlap : {false, true}) {
      const auto s = (interface.*operation)(
          as_const(source_overlap ? shared.view : other.view),
          as_const(source_overlap ? other.view : shared.view), offset);
      passed &= expect(s.code == StatusCode::invalid_plan &&
                           shared.storage == original_shared,
                       "partial diffusion input/output overlap rejects without writes");
    }
  }
  if (!passed) return false;
  const std::size_t local_cell_count =
      static_cast<std::size_t>(cells.x) * cells.y * cells.z;
  const std::uint8_t ghosts = fixture.boundary.maximum_halo_reach();
  const Span<const std::uint8_t> region = fixture.topology.region();

  // A non-constant fluid pressure makes the quadratic zero-normal ghost
  // distinguishable from the owning-fluid replacement. Solid values remain
  // deliberately arbitrary so the Cartesian stencil contains the deleted
  // solid-neighbour contribution that IBM must replace.
  ForceOwnedField& pressure = fixture.pressure;
  for (std::int32_t z = -ghosts; z < cells.z + ghosts; ++z)
    for (std::int32_t y = -ghosts; y < cells.y + ghosts; ++y)
      for (std::int32_t x = -ghosts; x < cells.x + ghosts; ++x) {
        const Int3 global{fixture.patch.begin.x + x,
                          fixture.patch.begin.y + y,
                          fixture.patch.begin.z + z};
        const double px = fixture.extrapolated_centre(fixture.geometry.x(),
                                                      global.x);
        const double py = fixture.extrapolated_centre(fixture.geometry.y(),
                                                      global.y);
        const double pz = fixture.extrapolated_centre(fixture.geometry.z(),
                                                      global.z);
        pressure.view.unchecked({x, y, z}, 0U) =
            37.25 + 0.5 * px + 0.25 * py * py - 0.1 * pz;
      }
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (region.data[flat(cells, cell)] ==
            static_cast<std::uint8_t>(RegionFlag::solid))
          pressure.view.unchecked(cell, 0U) =
              100.0 + 0.75 * static_cast<double>(x) +
              1.5 * static_cast<double>(y) +
              2.25 * static_cast<double>(z);
      }
  ForceOwnedField pressure_gradient =
      make_force_field(24U, cells, 3U, 0U, 25U, 115U);
  const auto gradient_from_pressure = [&](ConstFieldView source,
                                          FieldView target) {
    const std::array<ConstFieldView, 1U> reads{source};
    const std::array<FieldView, 1U> writes{target};
    const KernelInvocation invocation{
        {reads.data(), reads.size()},
        {writes.data(), writes.size()},
        {{0, 0, 0}, cells},
        0U,
        0U,
        1U,
        0U,
        nullptr};
    return cartesian_gradient(kernels, invocation);
  };
  passed &= expect(
      gradient_from_pressure(as_const(pressure.view), pressure_gradient.view),
      "ordinary Cartesian pressure gradient computes before IBM correction");
  const std::vector<double> pressure_gradient_before =
      pressure_gradient.storage;
  std::vector<double> expected_pressure_gradient = pressure_gradient_before;
  std::vector<std::size_t> owner_link_counts(local_cell_count, 0U);
  const Span<const ImmersedLink> links = fixture.topology.links();
  const Span<const IbmInterfaceLinkMetric> physical_links =
      fixture.topology.interface_metric().links();
  const Span<const BoundaryStencilLink> rows = fixture.boundary.links();
  const double inverse_width = 1.0 / fixture.geometry.x().uniform_width();
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const ImmersedLink& link = links.data[rows.data[index].topology_link];
    const std::uint8_t axis = face_axis(link.direction);
    const std::size_t owner = flat(cells, link.fluid_local_index);
    ++owner_link_counts[owner];
    const double derivative_weight =
        (positive_face(link.direction) ? 0.5 : -0.5) * inverse_width;
    double ghost = 0.0;
    passed &= expect(
        evaluate_quadratic_row(
            fixture.boundary.reconstruction(),
            rows.data[index].zero_normal_value_row, as_const(pressure.view),
            0U, 0.0, 0.0, ghost),
        "quadratic zero-normal pressure ghost evaluates");
    const double correction = derivative_weight *
                              (ghost - pressure.view.unchecked(
                                           link.solid_local_index, 0U));
    expected_pressure_gradient[flat(cells, link.fluid_local_index, axis)] +=
        correction;
  }
  const std::size_t maximum_owner_links =
      *std::max_element(owner_link_counts.begin(), owner_link_counts.end());
  passed &= expect(maximum_owner_links > 0U,
                   "IBM fixture contains immersed-link owners for the oracle");
  passed &= expect(
      interface.correct_pressure_gradient(as_const(pressure.view),
                                          pressure_gradient.view),
      "IBM pressure-gradient correction applies quadratic zero-normal ghost");
  double maximum_gradient_error = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (region.data[flat(cells, cell)] !=
            static_cast<std::uint8_t>(RegionFlag::fluid))
          continue;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const std::size_t entry = flat(cells, cell, component);
          maximum_gradient_error = std::max(
              maximum_gradient_error,
              std::abs(pressure_gradient.storage[entry] -
                       expected_pressure_gradient[entry]));
        }
      }
  passed &= expect(maximum_gradient_error < 5.0e-12,
                   "IBM pressure correction matches the quadratic-row oracle");

  ForceOwnedField dirichlet_probe =
      make_force_field(9U, cells, 1U, ghosts, 10U, 100U);
  std::fill(dirichlet_probe.storage.begin(), dirichlet_probe.storage.end(),
            1.0);
  const Span<const BoundaryStencilLink> boundary_rows =
      fixture.boundary.links();
  double minimum_dirichlet_derivative =
      std::numeric_limits<double>::infinity();
  for (std::size_t index = 0U; index < boundary_rows.size; ++index) {
    double derivative = 0.0;
    passed &= expect(
        evaluate_quadratic_row(
            fixture.boundary.reconstruction(),
            boundary_rows.data[index].wall_normal_gradient_row,
            as_const(dirichlet_probe.view), 0U, 0.0, 0.0, derivative),
        "Dirichlet wall-gradient row evaluates for each immersed link");
    minimum_dirichlet_derivative =
        std::min(minimum_dirichlet_derivative, derivative);
  }
  passed &= expect(
      std::isfinite(minimum_dirichlet_derivative) &&
          minimum_dirichlet_derivative > 1.0e-8,
      "zero wall value gives positive solid-to-fluid derivatives for a "
      "constant positive fluid field");
  ForceOwnedField velocity =
      make_force_field(10U, cells, 3U, ghosts, 11U, 101U);
  ForceOwnedField viscosity =
      make_force_field(11U, cells, 1U, ghosts, 12U, 102U);
  ForceOwnedField density =
      make_force_field(18U, cells, 1U, 0U, 19U, 109U);
  ForceOwnedField molecular =
      make_force_field(19U, cells, 1U, 0U, 20U, 110U);
  ForceOwnedField gradient =
      make_force_field(15U, cells, 9U, ghosts, 16U, 106U);
  ForceOwnedField diagonal =
      make_force_field(12U, cells, 3U, 0U, 13U, 103U);
  ForceOwnedField rhs =
      make_force_field(13U, cells, 3U, 0U, 14U, 104U);
  ForceOwnedField residual =
      make_force_field(14U, cells, 3U, 0U, 15U, 105U);
  const double width = fixture.geometry.x().uniform_width();
  constexpr double mu = 2.0;
  for (std::int32_t z = -ghosts; z < cells.z + ghosts; ++z)
    for (std::int32_t y = -ghosts; y < cells.y + ghosts; ++y)
      for (std::int32_t x = -ghosts; x < cells.x + ghosts; ++x) {
        const Int3 global{fixture.patch.begin.x + x,
                          fixture.patch.begin.y + y,
                          fixture.patch.begin.z + z};
        const double px = fixture.extrapolated_centre(fixture.geometry.x(),
                                                      global.x);
        const double py = fixture.extrapolated_centre(fixture.geometry.y(),
                                                      global.y);
        const double pz = fixture.extrapolated_centre(fixture.geometry.z(),
                                                      global.z);
        viscosity.view.unchecked({x, y, z}, 0U) = mu;
        const std::array<double, 9U> exact_gradient{
            1.0, 0.4 * py, 0.0, 0.0, 1.0, 0.2 * pz,
            0.3 * px, 0.0, 1.0};
        for (std::uint8_t component = 0U; component < 9U; ++component)
          gradient.view.unchecked({x, y, z}, component) =
              exact_gradient[component];
        velocity.view.unchecked({x, y, z}, 0U) =
            1.0 + px + 0.2 * py * py;
        velocity.view.unchecked({x, y, z}, 1U) =
            -0.4 + py + 0.1 * pz * pz;
        velocity.view.unchecked({x, y, z}, 2U) =
            0.3 + pz + 0.15 * px * px;
      }

  ForceOwnedField pressure_work_rate =
      make_force_field(27U, cells, 1U, 0U, 28U, 118U);
  std::fill(pressure_work_rate.storage.begin(),
            pressure_work_rate.storage.end(), 3.0);
  const std::vector<double> initial_pressure_work_rate =
      pressure_work_rate.storage;
  std::vector<double> expected_pressure_work_rate =
      initial_pressure_work_rate;
  double total_pressure_work_delta = 0.0;
  double maximum_pressure_work_delta = 0.0;
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const BoundaryStencilLink& row = rows.data[index];
    const ImmersedLink& link = links.data[row.topology_link];
    const std::uint8_t axis = face_axis(link.direction);
    const double derivative_weight =
        (positive_face(link.direction) ? 0.5 : -0.5) * inverse_width;
    double pressure_ghost = 0.0;
    passed &= expect(
        evaluate_quadratic_row(fixture.boundary.reconstruction(),
                               row.zero_normal_value_row,
                               as_const(pressure.view), 0U, 0.0, 0.0,
                               pressure_ghost),
        "pressure-work oracle evaluates the quadratic zero-normal ghost");
    const double pressure_gradient_correction =
        derivative_weight *
        (pressure_ghost - pressure.view.unchecked(link.solid_local_index, 0U));
    const double rate_correction =
        velocity.view.unchecked(link.fluid_local_index, axis) *
        pressure_gradient_correction;
    const std::size_t owner = flat(cells, link.fluid_local_index);
    expected_pressure_work_rate[owner] += rate_correction;
    total_pressure_work_delta += std::abs(rate_correction);
    maximum_pressure_work_delta = std::max(
        maximum_pressure_work_delta,
        std::abs(expected_pressure_work_rate[owner] -
                 initial_pressure_work_rate[owner]));
  }
  passed &= expect(total_pressure_work_delta > 1.0e-8 &&
                       maximum_pressure_work_delta > 1.0e-8,
                   "IBM pressure-work oracle has a nonzero per-link correction");
  passed &= expect(
      interface.correct_pressure_work(as_const(pressure.view),
                                      as_const(velocity.view),
                                      pressure_work_rate.view),
      "IBM pressure-work correction applies after parameter validation");
  double maximum_pressure_work_error = 0.0;
  for (std::size_t index = 0U; index < expected_pressure_work_rate.size();
       ++index)
    maximum_pressure_work_error = std::max(
        maximum_pressure_work_error,
        std::abs(pressure_work_rate.storage[index] -
                 expected_pressure_work_rate[index]));
  passed &= expect(maximum_pressure_work_error < 5.0e-12,
                   "IBM pressure-work correction matches the link-level oracle");
  std::fill(pressure_work_rate.storage.begin(), pressure_work_rate.storage.end(), 0.0);
  passed &= expect(interface.correct_pressure_work(as_const(pressure.view),
                       as_const(velocity.view), pressure_work_rate.view) &&
                       interface_only(pressure_work_rate),
                   "zero-based pressure work is supported only on unique interface cells");

  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const Int3 cell{x, y, z};
          density.view.unchecked(cell, 0U) = 1.2;
          molecular.view.unchecked(cell, 0U) = 1.8e-5;
          diagonal.view.unchecked(cell, component) = 9.0;
          residual.view.unchecked(cell, component) = 2.0 + component;
          rhs.view.unchecked(cell, component) =
              9.0 * velocity.view.unchecked(cell, component) -
              residual.view.unchecked(cell, component);
        }

  ForceOwnedField velocity_gradient_probe =
      make_force_field(26U, cells, 9U, 0U, 27U, 117U);
  {
    const std::array<ConstFieldView, 1U> reads{as_const(velocity.view)};
    const std::array<FieldView, 1U> writes{velocity_gradient_probe.view};
    passed &= expect(
        cartesian_gradient(
            kernels,
            {{reads.data(), reads.size()}, {writes.data(), writes.size()},
             {{0, 0, 0}, cells}, 0U, 0U, 3U, 0U, nullptr}),
        "ordinary Cartesian velocity gradient computes before IBM correction");
  }
  std::vector<double> expected_velocity_gradient =
      velocity_gradient_probe.storage;
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const BoundaryStencilLink& row = rows.data[index];
    const ImmersedLink& link = links.data[row.topology_link];
    const std::uint8_t axis = face_axis(link.direction);
    const double derivative_weight =
        (positive_face(link.direction) ? 0.5 : -0.5) * inverse_width;
    for (std::uint8_t component = 0U; component < 3U; ++component) {
      double ghost = 0.0;
      passed &= expect(
          evaluate_quadratic_row(fixture.boundary.reconstruction(),
                                 row.dirichlet_value_row,
                                 as_const(velocity.view), component, 0.0, 0.0,
                                 ghost),
          "quadratic no-slip velocity ghost evaluates");
      const std::uint8_t gradient_component =
          static_cast<std::uint8_t>(3U * component + axis);
      expected_velocity_gradient[flat(cells, link.fluid_local_index,
                                      gradient_component)] +=
          derivative_weight *
          (ghost - velocity.view.unchecked(link.solid_local_index,
                                           component));
    }
  }
  passed &= expect(
      interface.correct_velocity_gradient(as_const(velocity.view),
                                          velocity_gradient_probe.view),
      "IBM velocity-gradient correction applies quadratic no-slip ghosts");
  double maximum_velocity_gradient_error = 0.0;
  for (std::size_t index = 0U; index < expected_velocity_gradient.size();
       ++index)
    maximum_velocity_gradient_error =
        std::max(maximum_velocity_gradient_error,
                 std::abs(velocity_gradient_probe.storage[index] -
                          expected_velocity_gradient[index]));
  passed &= expect(maximum_velocity_gradient_error < 5.0e-12,
                   "IBM velocity gradient matches the quadratic-row oracle");

  const std::vector<double> initial_diagonal = diagonal.storage;
  const std::vector<double> initial_residual = residual.storage;
  const std::vector<double> initial_rhs = rhs.storage;
  std::vector<double> expected_diagonal = initial_diagonal;
  std::vector<double> expected_residual = initial_residual;
  std::vector<double> expected_rhs = initial_rhs;
  const double transmissibility = mu * width;
  double cartesian_traction_mutation_gap = 0.0;
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const ImmersedLink& link = links.data[rows.data[index].topology_link];
    const IbmInterfaceLinkMetric& physical =
        physical_links.data[rows.data[index].topology_link];
    const Int3 face = face_index(link);
    std::uint8_t axis = 2U;
    if (link.direction == ImmersedFaceDirection::x_negative ||
        link.direction == ImmersedFaceDirection::x_positive)
      axis = 0U;
    else if (link.direction == ImmersedFaceDirection::y_negative ||
             link.direction == ImmersedFaceDirection::y_positive)
      axis = 1U;
    Int3 left = face;
    if (axis == 0U)
      --left.x;
    else if (axis == 1U)
      --left.y;
    else
      --left.z;
    double normal_derivative[3]{};
    for (std::uint8_t component = 0U; component < 3U; ++component) {
      passed &= expect(evaluate_quadratic_row(
                           fixture.boundary.reconstruction(),
                           rows.data[index].wall_normal_gradient_row,
                           as_const(velocity.view), component, 0.0, 0.0,
                           normal_derivative[component]),
                       "wall-gradient row reconstructs each link independently");
    }
    const double normal[3]{link.solid_to_fluid_normal.x,
                           link.solid_to_fluid_normal.y,
                           link.solid_to_fluid_normal.z};
    const double normal_l1 = std::abs(normal[0U]) +
                             std::abs(normal[1U]) +
                             std::abs(normal[2U]);
    const double normal_projection = normal_derivative[0U] * normal[0U] +
                                     normal_derivative[1U] * normal[1U] +
                                     normal_derivative[2U] * normal[2U];
    const QuadraticStencilPlan& reconstruction =
        fixture.boundary.reconstruction();
    const QuadraticAffineRow& derivative_row =
        reconstruction.rows().data[rows.data[index].wall_normal_gradient_row];
    const QuadraticStencilGroup& derivative_group =
        reconstruction.groups().data[derivative_row.group];
    const double fluid_x = fixture.extrapolated_centre(
        fixture.geometry.x(), link.fluid_global_index.x);
    const double fluid_y = fixture.extrapolated_centre(
        fixture.geometry.y(), link.fluid_global_index.y);
    const double fluid_z = fixture.extrapolated_centre(
        fixture.geometry.z(), link.fluid_global_index.z);
    const double wall_dx = fluid_x - link.wall_point.x;
    const double wall_dy = fluid_y - link.wall_point.y;
    const double wall_dz = fluid_z - link.wall_point.z;
    const double inverse_wall_distance =
        1.0 / std::sqrt(wall_dx * wall_dx + wall_dy * wall_dy +
                        wall_dz * wall_dz);
    double correction_l1 = 0.0;
    std::size_t owner_count = 0U;
    for (std::size_t donor = 0U;
         donor < derivative_group.quality.donor_count; ++donor) {
      const bool owner =
          reconstruction.donor_global_cells()
              .data[derivative_group.donor_begin + donor] == link.fluid_cell;
      correction_l1 += std::abs(
          reconstruction.weights().data[derivative_row.weight_begin + donor] -
          (owner ? inverse_wall_distance : 0.0));
      owner_count += owner ? 1U : 0U;
    }
    if (owner_count == 0U) correction_l1 += inverse_wall_distance;
    const double gradient_majorant = inverse_wall_distance + correction_l1;
    passed &= expect(owner_count <= 1U && std::isfinite(gradient_majorant) &&
                         gradient_majorant > 0.0,
                     "wall-gradient low-order-plus-correction majorant is "
                     "positive and finite");
    const double wall_transmissibility =
        mu * physical.physical_quadrature_area * gradient_majorant;
    double pressure_ghost = 0.0;
    passed &= expect(
        evaluate_quadratic_row(fixture.boundary.reconstruction(),
                               rows.data[index].zero_normal_value_row,
                               as_const(pressure.view), 0U, 0.0, 0.0,
                               pressure_ghost),
        "momentum pressure ghost uses the quadratic zero-normal row");
    for (std::uint8_t component = 0U; component < 3U; ++component) {
      const double divergence =
          0.5 * (gradient.view.unchecked(left, 0U) +
                 gradient.view.unchecked(face, 0U)) +
          0.5 * (gradient.view.unchecked(left, 4U) +
                 gradient.view.unchecked(face, 4U)) +
          0.5 * (gradient.view.unchecked(left, 8U) +
                 gradient.view.unchecked(face, 8U));
      const std::uint8_t transpose =
          static_cast<std::uint8_t>(3U * axis + component);
      const double cross =
          mu * (0.5 * (gradient.view.unchecked(left, transpose) +
                       gradient.view.unchecked(face, transpose)) -
                (axis == component ? (2.0 / 3.0) * divergence : 0.0));
      const bool positive =
          link.direction == ImmersedFaceDirection::x_positive ||
          link.direction == ImmersedFaceDirection::y_positive ||
          link.direction == ImmersedFaceDirection::z_positive;
      const double pressure_derivative_weight =
          (positive ? 0.5 : -0.5) / width;
      const double pressure_gradient_correction =
          pressure_derivative_weight *
          (pressure_ghost -
           pressure.view.unchecked(link.solid_local_index, 0U));
      const double pressure_force_correction =
          width * width * width * pressure_gradient_correction;
      const double regular =
          transmissibility *
              (velocity.view.unchecked(link.fluid_local_index, component) -
               velocity.view.unchecked(link.solid_local_index, component)) +
          (positive ? -1.0 : 1.0) * cross *
              link.cartesian_control_face_area;
      double normal_stress = 0.0;
      double normal_row_l1 = 0.0;
      for (std::uint8_t derivative = 0U; derivative < 3U; ++derivative) {
        const double moment =
            physical.normal_second_moment[3U * component + derivative];
        normal_stress += moment * normal_derivative[derivative];
        normal_row_l1 += std::abs(moment);
      }
      normal_row_l1 /= physical.physical_quadrature_area;
      const double desired =
          mu * (physical.physical_quadrature_area *
                    normal_derivative[component] +
                (1.0 / 3.0) * normal_stress);
      const double cartesian_desired =
          mu * link.cartesian_control_face_area *
          (normal_derivative[component] +
           (1.0 / 3.0) * normal[component] * normal_projection);
      const double cartesian_diagonal =
          mu * link.cartesian_control_face_area * gradient_majorant *
          (1.0 + (1.0 / 3.0) * std::abs(normal[component]) * normal_l1);
      const double physical_diagonal =
          wall_transmissibility * (1.0 + (1.0 / 3.0) * normal_row_l1);
      cartesian_traction_mutation_gap = std::max(
          cartesian_traction_mutation_gap,
          std::max(std::abs(cartesian_desired - desired),
                   std::abs(cartesian_diagonal - physical_diagonal)));
      const double correction =
          desired - regular +
          (component == axis ? pressure_force_correction : 0.0);
      const std::size_t equation =
          flat(cells, link.fluid_local_index, component);
      expected_diagonal[equation] +=
          wall_transmissibility *
              (1.0 + (1.0 / 3.0) * normal_row_l1) -
          transmissibility;
      expected_residual[equation] += correction;
      expected_rhs[equation] =
          expected_diagonal[equation] *
              velocity.view.unchecked(link.fluid_local_index, component) -
          expected_residual[equation];
    }
  }
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (region.data[flat(cells, cell)] ==
            static_cast<std::uint8_t>(RegionFlag::solid))
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            expected_diagonal[flat(cells, cell, component)] = 1.0;
            expected_residual[flat(cells, cell, component)] =
                velocity.view.unchecked(cell, component);
            expected_rhs[flat(cells, cell, component)] = 0.0;
          }
      }
  passed &= expect(interface.constrain_momentum(
                       as_const(velocity.view), as_const(gradient.view),
                       as_const(pressure.view),
                       as_const(density.view), as_const(molecular.view),
                       as_const(viscosity.view), nullptr,
                       {diagonal.view, rhs.view, residual.view}),
                   "momentum interface replacement applies");
  double maximum_error = 0.0;
  for (std::size_t index = 0U; index < expected_rhs.size(); ++index) {
    maximum_error = std::max(
        maximum_error, std::abs(rhs.storage[index] - expected_rhs[index]));
    maximum_error = std::max(
        maximum_error,
        std::abs(residual.storage[index] - expected_residual[index]));
    maximum_error = std::max(
        maximum_error,
        std::abs(diagonal.storage[index] - expected_diagonal[index]));
  }
  passed &= expect(maximum_error < 5.0e-12,
                   "per-link momentum residual and positive deferred diagonal "
                   "match the independent oracle");
  passed &= expect(
      cartesian_traction_mutation_gap > 1.0e-6,
      "mutation feeding Cartesian control area to physical traction fails");

  TurbulenceFixture wall_fixture;
  passed &= expect(wall_fixture.initialize(TurbulencePlanSpec{}),
                   "Vreman wall-function authority compiles");
  ForceOwnedField wall_diagonal =
      make_force_field(20U, cells, 3U, 0U, 21U, 111U);
  ForceOwnedField wall_rhs =
      make_force_field(21U, cells, 3U, 0U, 22U, 112U);
  ForceOwnedField wall_residual =
      make_force_field(22U, cells, 3U, 0U, 23U, 113U);
  std::copy(initial_diagonal.begin(), initial_diagonal.end(),
            wall_diagonal.storage.begin());
  std::copy(initial_rhs.begin(), initial_rhs.end(), wall_rhs.storage.begin());
  std::copy(initial_residual.begin(), initial_residual.end(),
            wall_residual.storage.begin());
  passed &= expect(interface.constrain_momentum(
                       as_const(velocity.view), as_const(gradient.view),
                       as_const(pressure.view),
                       as_const(density.view), as_const(molecular.view),
                       as_const(viscosity.view), &wall_fixture.plan,
                       {wall_diagonal.view, wall_rhs.view,
                        wall_residual.view}),
                   "Vreman wall law enters the same per-link momentum operator");
  double wall_change = 0.0;
  double wall_diagonal_change = 0.0;
  for (std::size_t index = 0U; index < wall_rhs.storage.size(); ++index) {
    wall_change = std::max(
        wall_change, std::abs(wall_rhs.storage[index] - rhs.storage[index]));
    if (region.data[index % local_cell_count] ==
        static_cast<std::uint8_t>(RegionFlag::fluid))
      wall_diagonal_change =
          std::max(wall_diagonal_change,
                   std::abs(wall_diagonal.storage[index] -
                            initial_diagonal[index]));
    passed &= expect(std::isfinite(wall_diagonal.storage[index]) &&
                         wall_diagonal.storage[index] > 0.0,
                     "wall-law replacement preserves a positive diagonal");
  }
  passed &= expect(wall_change > 1.0e-8,
                   "wall-law mutation changes immersed tangential traction");
  passed &= expect(wall_diagonal_change > 1.0e-8,
                   "wall-law replacement removes the Cartesian solid-face "
                   "diagonal");

  ForceOwnedField transported =
      make_force_field(16U, cells, 1U, ghosts, 17U, 107U);
  ForceOwnedField diffusion_rate =
      make_force_field(17U, cells, 1U, 0U, 18U, 108U);
  for (std::int32_t z = -ghosts; z < cells.z + ghosts; ++z)
    for (std::int32_t y = -ghosts; y < cells.y + ghosts; ++y)
      for (std::int32_t x = -ghosts; x < cells.x + ghosts; ++x)
        transported.view.unchecked({x, y, z}, 0U) =
            velocity.view.unchecked({x, y, z}, 0U);
  std::fill(diffusion_rate.storage.begin(), diffusion_rate.storage.end(),
            5.0);
  std::vector<double> expected_rate = diffusion_rate.storage;
  const double volume = width * width * width;
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const ImmersedLink& link = links.data[rows.data[index].topology_link];
    double ghost = 0.0;
    passed &= expect(evaluate_quadratic_row(
                         fixture.boundary.reconstruction(),
                         rows.data[index].zero_normal_value_row,
                         as_const(transported.view), 0U, 0.0, 0.0, ghost),
                     "zero-normal row reconstructs each scalar link");
    expected_rate[flat(cells, link.fluid_local_index)] +=
        transmissibility *
        (ghost - transported.view.unchecked(link.solid_local_index, 0U)) /
        volume;
  }
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        if (region.data[flat(cells, {x, y, z})] ==
            static_cast<std::uint8_t>(RegionFlag::solid))
          expected_rate[flat(cells, {x, y, z})] = 0.0;
  passed &= expect(interface.correct_zero_normal_diffusion(
                       as_const(transported.view), as_const(viscosity.view),
                       diffusion_rate.view),
                   "adiabatic/zero-flux diffusion replacement applies");
  maximum_error = 0.0;
  for (std::size_t index = 0U; index < expected_rate.size(); ++index)
    maximum_error = std::max(
        maximum_error,
        std::abs(diffusion_rate.storage[index] - expected_rate[index]));
  passed &= expect(maximum_error < 5.0e-12,
                   "zero-normal diffusion replacement has correct sign and units");
  std::fill(diffusion_rate.storage.begin(), diffusion_rate.storage.end(), 0.0);
  passed &= expect(interface.correct_zero_normal_diffusion(as_const(transported.view),
                       as_const(viscosity.view), diffusion_rate.view) && interface_only(diffusion_rate),
                   "zero-based thermal correction is supported only on unique interface cells");

  // The thermal route keeps every donor strictly positive but deliberately
  // makes the quadratic data non-smooth.  A negative quadratic weight then
  // makes the unconstrained ghost leave the positive donor envelope, while
  // the bounded evaluator projects that same row back to the envelope.
  ForceOwnedField bounded_transported =
      make_force_field(29U, cells, 1U, ghosts, 30U, 120U);
  std::fill(bounded_transported.storage.begin(),
            bounded_transported.storage.end(), 1.0);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        if (region.data[flat(cells, {x, y, z})] ==
            static_cast<std::uint8_t>(RegionFlag::solid))
          bounded_transported.view.unchecked({x, y, z}, 0U) = 2.0;
  const QuadraticStencilPlan& reconstruction =
      fixture.boundary.reconstruction();
  const Span<const QuadraticAffineRow> reconstruction_rows =
      reconstruction.rows();
  const Span<const QuadraticStencilGroup> reconstruction_groups =
      reconstruction.groups();
  const Span<const Int3> reconstruction_donors =
      reconstruction.donor_local_indices();
  const Span<const double> reconstruction_weights = reconstruction.weights();
  std::uint32_t selected_boundary_row = kInvalidIbmIndex;
  std::uint32_t selected_donor = kInvalidIbmIndex;
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const std::uint32_t value_row = rows.data[index].zero_normal_value_row;
    if (value_row >= reconstruction_rows.size)
      continue;
    const QuadraticAffineRow& affine = reconstruction_rows.data[value_row];
    if (affine.group >= reconstruction_groups.size)
      continue;
    const QuadraticStencilGroup& group =
        reconstruction_groups.data[affine.group];
    for (std::size_t donor = 0U; donor < group.quality.donor_count; ++donor) {
      const std::size_t offset =
          static_cast<std::size_t>(affine.weight_begin) + donor;
      if (offset >= reconstruction_weights.size ||
          reconstruction_weights.data[offset] >= -1.0e-12)
        continue;
      selected_boundary_row = static_cast<std::uint32_t>(index);
      selected_donor = group.donor_begin + static_cast<std::uint32_t>(donor);
      break;
    }
    if (selected_boundary_row != kInvalidIbmIndex)
      break;
  }
  passed &= expect(selected_boundary_row != kInvalidIbmIndex,
                   "bounded thermal fixture finds a negative donor weight");
  if (selected_boundary_row == kInvalidIbmIndex)
    return false;
  const Int3 non_smooth_donor = reconstruction_donors.data[selected_donor];
  bounded_transported.view.unchecked(non_smooth_donor, 0U) = 1.0e-12;
  const QuadraticAffineRow& selected_affine = reconstruction_rows.data[
      rows.data[selected_boundary_row].zero_normal_value_row];
  const QuadraticStencilGroup& selected_group =
      reconstruction_groups.data[selected_affine.group];
  double donor_maximum = 0.0;
  for (std::size_t donor = 0U; donor < selected_group.quality.donor_count;
       ++donor)
    donor_maximum = std::max(
        donor_maximum,
        bounded_transported.view.unchecked(
            reconstruction_donors.data[selected_group.donor_begin + donor],
            0U));
  double raw_thermal_ghost = 0.0;
  passed &= expect(
      evaluate_quadratic_row(
          reconstruction, rows.data[selected_boundary_row].zero_normal_value_row,
          as_const(bounded_transported.view), 0U, 0.0, 0.0,
          raw_thermal_ghost),
      "ordinary zero-normal thermal ghost evaluates the donor oracle");
  passed &= expect(raw_thermal_ghost > donor_maximum + 1.0e-12,
                   "ordinary zero-normal thermal ghost overshoots positive "
                   "donor maximum");
  double bounded_thermal_ghost = 0.0;
  passed &= expect(
      evaluate_positive_bounded_quadratic_row(
          reconstruction, rows.data[selected_boundary_row].zero_normal_value_row,
          as_const(bounded_transported.view), 0U, bounded_thermal_ghost),
      "positive-bounded thermal ghost evaluates the donor envelope");
  passed &= expect(std::abs(bounded_thermal_ghost - donor_maximum) < 5.0e-12,
                   "positive-bounded thermal ghost uses donor maximum");

  ForceOwnedField bounded_rate =
      make_force_field(31U, cells, 1U, 0U, 32U, 122U);
  std::fill(bounded_rate.storage.begin(), bounded_rate.storage.end(), 5.0);
  std::vector<double> expected_bounded_rate = bounded_rate.storage;
  double total_bounded_correction = 0.0;
  double maximum_bounded_correction = 0.0;
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const BoundaryStencilLink& row = rows.data[index];
    const ImmersedLink& link = links.data[row.topology_link];
    double ghost = 0.0;
    passed &= expect(
        evaluate_positive_bounded_quadratic_row(
            reconstruction, row.zero_normal_value_row,
            as_const(bounded_transported.view), 0U, ghost),
        "positive-bounded row evaluates each scalar link");
    const double correction =
        transmissibility *
        (ghost - bounded_transported.view.unchecked(link.solid_local_index, 0U)) /
        volume;
    const std::size_t owner = flat(cells, link.fluid_local_index);
    expected_bounded_rate[owner] += correction;
    total_bounded_correction += std::abs(correction);
    maximum_bounded_correction = std::max(
        maximum_bounded_correction,
        std::abs(expected_bounded_rate[owner] - 5.0));
  }
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        if (region.data[flat(cells, {x, y, z})] ==
            static_cast<std::uint8_t>(RegionFlag::solid))
          expected_bounded_rate[flat(cells, {x, y, z})] = 0.0;
  passed &= expect(total_bounded_correction > 1.0e-8 &&
                       maximum_bounded_correction > 1.0e-8,
                   "bounded thermal diffusion oracle has a nonzero correction");
  passed &= expect(
      interface.correct_positive_bounded_zero_normal_diffusion(
          as_const(bounded_transported.view), as_const(viscosity.view),
          bounded_rate.view),
      "bounded thermal diffusion validates and applies");
  double maximum_bounded_error = 0.0;
  for (std::size_t index = 0U; index < expected_bounded_rate.size(); ++index)
    maximum_bounded_error = std::max(
        maximum_bounded_error,
        std::abs(bounded_rate.storage[index] - expected_bounded_rate[index]));
  passed &= expect(maximum_bounded_error < 5.0e-12,
                   "bounded thermal diffusion matches exact per-link oracle");

  // Accepted spatial history must use the same IBM thermal operator as the
  // target energy equation, not the positivity-only donor-envelope closure.
  // The preceding fixture proves that these two reconstructions differ.
  {
    ThermodynamicsPlan thermodynamics;
    TransportPlan transport;
    ContributionRegistry contributions;
    EquationPlanSet equations;
    EquationPlanSpec spec;
    spec.density = 60U;
    spec.velocity = physical_boundary.velocity_field();
    spec.pressure_perturbation = physical_boundary.pressure_field();
    spec.enthalpy = physical_boundary.enthalpy_field();
    spec.temperature = 59U;
    spec.effective_viscosity = 61U;
    spec.pressure_compressibility = 62U;
    spec.velocity_gradient = 63U;
    spec.pressure_reference = model.pressure_reference;
    spec.closed_mass_service_stage = 1U;
    spec.maximum_cells_per_rank = local_cell_count;
    const std::array<FieldId, 8U> declared{
        spec.density, spec.velocity, spec.pressure_perturbation, spec.enthalpy,
        spec.temperature, spec.effective_viscosity,
        spec.pressure_compressibility, spec.velocity_gradient};
    const bool compiled = ThermodynamicsPlan::compile(model.thermophysics, {}, thermodynamics) &&
        TransportPlan::compile(model.thermophysics, thermodynamics, transport) &&
        contributions.configure({declared.data(), declared.size()}) &&
        contributions.freeze() && EquationPlanSet::compile(
            MPI_COMM_SELF, schemes, fixture.geometry, fixture.patch,
            physical_boundary, contributions, thermodynamics, transport, spec, equations);
    if (!expect(compiled, "accepted IBM rate equation plans compile")) return false;
    auto t = make_force_field(spec.temperature, cells, 1U, ghosts, 2011U, 3011U);
    auto h = make_force_field(spec.enthalpy, cells, 1U, ghosts, 2012U, 3012U);
    auto rho = make_force_field(spec.density, cells, 1U, ghosts, 2013U, 3013U);
    auto u = make_force_field(spec.velocity, cells, 3U, ghosts, 2014U, 3014U);
    auto p = make_force_field(spec.pressure_perturbation, cells, 1U, ghosts, 2015U, 3015U);
    auto grad = make_force_field(spec.velocity_gradient, cells, 9U, 0U, 2016U, 3016U);
    auto mu = make_force_field(spec.effective_viscosity, cells, 1U, 1U, 2017U, 3017U);
    auto lambda = make_force_field(64U, cells, 1U, 1U, 2018U, 3018U);
    auto rhs = make_force_field(65U, cells, 1U, 0U, 2019U, 3019U);
    auto scratch = make_force_field(66U, cells, 1U, 0U, 2020U, 3020U);
    auto scalar_d = make_force_field(67U, cells, 1U, 1U, 2021U, 3021U);
    auto target_rate = make_force_field(68U, cells, 1U, 0U, 2022U, 3022U);
    for (auto* zero : {&u, &p, &grad})
      std::fill(zero->storage.begin(), zero->storage.end(), 0.0);
    std::fill(mu.storage.begin(), mu.storage.end(), 1.8e-5);
    std::fill(lambda.storage.begin(), lambda.storage.end(), 1.0);
    const double gas_constant = kUniversalGasConstant / 28.96546;
    for (std::int32_t z = -ghosts; z < cells.z + ghosts; ++z)
      for (std::int32_t y = -ghosts; y < cells.y + ghosts; ++y)
        for (std::int32_t x = -ghosts; x < cells.x + ghosts; ++x) {
          const Int3 cell{x, y, z};
          const double temperature = 300.0 +
              100.0 * bounded_transported.view.unchecked(cell, 0U);
          t.view.unchecked(cell, 0U) = temperature;
          h.view.unchecked(cell, 0U) = 3.5 * gas_constant * temperature;
          rho.view.unchecked(cell, 0U) = 101325.0 / (gas_constant * temperature);
        }
    const auto history = [](FieldView view) {
      return PrimitiveHistory{as_const(view), as_const(view), as_const(view)};
    };
    EquationStateView state;
    state.density = history(rho.view);
    state.velocity = history(u.view);
    state.pressure_perturbation = history(p.view);
    state.enthalpy = history(h.view);
    state.temperature = history(t.view);
    state.pressure_reference = state.accepted_pressure_reference =
        state.previous_pressure_reference = 101325.0;
    EquationMaterialView material;
    material.molecular_viscosity = material.effective_viscosity = as_const(mu.view);
    material.thermal_conductivity = as_const(lambda.view);
    ThermophysicalRateInput input{state, material, as_const(grad.view),
                                  {1.0, -1.0, 0.0, 1U}, 1U, {}, &interface};
    ThermophysicalRateOutput output{rhs.view, {}, {}, scratch.view, scalar_d.view};
    ThermophysicalRateCertificate certificate;
    passed &= expect(evaluate_thermophysical_rates(equations, input, output, certificate),
                     "accepted IBM spatial energy rate evaluates");
    const std::array<ConstFieldView, 1U> reads{as_const(t.view)};
    const std::array<FieldView, 1U> writes{target_rate.view};
    passed &= expect(cartesian_diffusion(kernels, as_const(lambda.view),
        {{reads.data(), reads.size()}, {writes.data(), writes.size()},
         {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U, nullptr}) &&
        interface.correct_zero_normal_diffusion(
            as_const(t.view), as_const(lambda.view), target_rate.view),
        "target energy equation's independently checked IBM thermal closure evaluates");
    double mismatch = 0.0;
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (region.data[flat(cells, cell)] == 0U) continue;
          mismatch = std::max(mismatch, std::abs(rhs.view.unchecked(cell, 0U) -
                                                target_rate.view.unchecked(cell, 0U)));
        }
    std::cerr << "accepted-IBM-thermal-rate-mismatch=" << mismatch << '\n';
    passed &= expect(mismatch < 1.0e-10,
                     "persisted energy history and target residual use one IBM thermal operator");
    // Frozen-lambda donor Jv: deltaT=T-350 contains both signs, and linearity
    // plus D(constant)=0 gives the independent identity D(deltaT)=D(T).
    // This checks the target donor response, not a claim that masked Schur Eh
    // already includes that response.
    const auto base_temperature = t.storage;
    constexpr double epsilon = 1.0e-3;
    const auto perturbed_diffusion = [&](double factor, FieldView destination) {
      for (std::size_t i = 0U; i < t.storage.size(); ++i)
        t.storage[i] = base_temperature[i] + factor * (base_temperature[i] - 350.0);
      const std::array<FieldView, 1U> destinations{destination};
      return cartesian_diffusion(kernels, as_const(lambda.view),
          {{reads.data(), reads.size()}, {destinations.data(), destinations.size()},
           {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U, nullptr}) &&
          interface.correct_zero_normal_diffusion(
              as_const(t.view), as_const(lambda.view), destination);
    };
    passed &= expect(perturbed_diffusion(epsilon, rhs.view) &&
                     perturbed_diffusion(-epsilon, scratch.view),
                     "positive target temperatures accept signed donor variations");
    double derivative_error = 0.0;
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (region.data[flat(cells, cell)] == 0U) continue;
          const double fd = (rhs.view.unchecked(cell, 0U) -
                             scratch.view.unchecked(cell, 0U)) / (2.0 * epsilon);
          derivative_error = std::max(derivative_error,
              std::abs(fd - target_rate.view.unchecked(cell, 0U)));
        }
    std::cerr << "IBM-donor-thermal-Jv-error=" << derivative_error << '\n';
    passed &= expect(derivative_error < 1.0e-6,
                     "frozen-material IBM donor thermal derivative matches finite differences");
    // A non-conjugate solid placeholder is not a material on the fluid side
    // of an adiabatic wall. Keep all fluid states and coefficients fixed.
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          if (region.data[flat(cells, {x, y, z})] == 0U)
            lambda.view.unchecked({x, y, z}, 0U) = 2.0;
    passed &= expect(perturbed_diffusion(0.0, rhs.view),
                     "legal solid placeholder conductivity perturbation evaluates");
    double placeholder_response = 0.0;
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (region.data[flat(cells, cell)] == 0U) continue;
          placeholder_response = std::max(placeholder_response,
              std::abs(rhs.view.unchecked(cell, 0U) -
                       target_rate.view.unchecked(cell, 0U)));
        }
    std::cerr << "IBM-solid-material-fluid-response=" << placeholder_response << '\n';
    passed &= expect(placeholder_response < 1.0e-10,
                     "adiabatic fluid equation is independent of solid placeholder material");
  }

  OwnedFace x_flux = make_face(CartesianAxis::x, cells, 201U);
  OwnedFace y_flux = make_face(CartesianAxis::y, cells, 201U);
  OwnedFace z_flux = make_face(CartesianAxis::z, cells, 201U);
  FaceFluxView flux{x_flux.view, y_flux.view, z_flux.view, 21U, {}};
  passed &= expect(interface.validate_interface_flux(as_const(flux)).code ==
                       StatusCode::numerical_failure,
                   "nonzero Restart/interface flux mutation is rejected");
  const Status constrained =
      interface.constrain_corrected_state(velocity.view, flux);
  passed &= expect(constrained,
                   "corrected state masks solids and interface flux");
  if (constrained)
    for (std::size_t index = 0U; index < links.size; ++index) {
      const ImmersedLink& link = links.data[index];
      passed &= expect(select(flux, link.direction)
                               .unchecked(face_index(link)) == 0.0,
                       "every fluid-solid face has exactly zero mass flux");
    }
  if (constrained)
    passed &= expect(interface.validate_interface_flux(as_const(flux)),
                     "zero interface flux validates for Restart");
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        if (region.data[flat(cells, {x, y, z})] ==
            static_cast<std::uint8_t>(RegionFlag::solid))
          for (std::uint8_t component = 0U; component < 3U; ++component)
            passed &= expect(velocity.view.unchecked({x, y, z}, component) ==
                                 0.0,
                             "solid velocity row is stationary");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  const bool passed = test_prescribed_interface_mass_flux() && run();
  MPI_Finalize();
  return passed ? 0 : 1;
}
