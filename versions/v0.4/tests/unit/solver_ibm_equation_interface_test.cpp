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
  const bool passed = run();
  MPI_Finalize();
  return passed ? 0 : 1;
}
