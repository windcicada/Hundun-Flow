// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/flow_immersed_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
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
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "src/ib_quadratic_reconstruction_detail.hpp"
#include "tests/support/stage3_decomposition_equality.hpp"
#include "tests/support/stage3_mms.hpp"
#include "tests/support/stage3_test_contracts.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

namespace {

using namespace hundun;

constexpr double kPi = 3.141592653589793238462643383279502884;

enum class EngineeringBody : std::uint8_t { cylinder, sphere };

struct EngineeringCase final {
  EngineeringBody body{};
  runtime::Int3 cells{};
  runtime::Int3 process_grid{};
  runtime::Real3 origin_m{};
  runtime::Real3 length_m{};
  double reynolds{};
  double final_time_s{};
  double statistics_window_s{};
  bool collect_cylinder_diagnostics{true};
  double fast_force_bound_N{std::numeric_limits<double>::infinity()};
};

struct EngineeringResult final {
  double total_drag_mean{};
  double pressure_drag_mean{};
  double viscous_drag_mean{};
  double lateral_y_mean{};
  double lateral_z_mean{};
  double total_drag_cv{};
  double pressure_drag_cv{};
  double viscous_drag_cv{};
  double recirculation_length_mean_m{};
  double recirculation_length_cv{};
  double upper_separation_angle_deg{};
  double lower_separation_angle_deg{};
  std::uint64_t classification_fingerprint{};
  std::uint64_t ghost_plan_fingerprint{};
  std::uint64_t wall_plan_fingerprint{};
  std::uint64_t accepted_steps{};
  double minimum_accepted_dt_s{};
  double maximum_accepted_dt_s{};
  std::uint64_t state_fingerprint{};
  std::vector<mesh::GlobalCellId> global_active_cell_ids;
  std::vector<double> global_velocity;
  immersed::ForceComponents final_force;
};

runtime::Int3 default_process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("unsupported engineering rank count");
}

void check_mpi(int code, const char *operation) {
  if (code != MPI_SUCCESS)
    throw runtime::Error(std::string(operation) + " failed");
}

void sum_in_place(const runtime::MpiContext &mpi, double *values,
                  std::size_t count) {
  HUNDUN_CHECK(count <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()));
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, values, static_cast<int>(count),
                          MPI_DOUBLE, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(engineering sum)");
}

double mean(const std::vector<double> &values) {
  HUNDUN_CHECK(!values.empty());
  for (const double value : values)
    HUNDUN_CHECK(std::isfinite(value));
  return std::accumulate(values.begin(), values.end(), 0.0) /
         static_cast<double>(values.size());
}

double coefficient_of_variation(const std::vector<double> &values) {
  const double average = mean(values);
  HUNDUN_CHECK(std::abs(average) > 0.0);
  double variance = 0.0;
  for (const double value : values) {
    const double delta = value - average;
    variance += delta * delta;
  }
  variance /= static_cast<double>(values.size());
  return std::sqrt(variance) / std::abs(average);
}

std::uint64_t fp64_bits(double value) noexcept {
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  std::uint64_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

void hash_u64(std::uint64_t &fingerprint, std::uint64_t value) noexcept {
  constexpr std::uint64_t prime = 1099511628211ULL;
  for (unsigned byte = 0U; byte < 8U; ++byte) {
    fingerprint ^= (value >> (byte * 8U)) & 0xffU;
    fingerprint *= prime;
  }
}

std::uint64_t
engineering_state_fingerprint(const std::vector<mesh::GlobalCellId> &active_ids,
                              const std::vector<double> &global_velocity) {
  std::uint64_t fingerprint = 14695981039346656037ULL;
  hash_u64(fingerprint, static_cast<std::uint64_t>(active_ids.size()));
  for (const auto global : active_ids) {
    const auto offset = static_cast<std::size_t>(global) * 3U;
    HUNDUN_CHECK(offset + 2U < global_velocity.size());
    hash_u64(fingerprint, global);
    for (std::size_t component = 0U; component < 3U; ++component)
      hash_u64(fingerprint, fp64_bits(global_velocity[offset + component]));
  }
  return fingerprint;
}

std::string fingerprint_hex(std::uint64_t fingerprint) {
  constexpr char digits[] = "0123456789abcdef";
  std::string result(16U, '0');
  for (std::size_t index = 0U; index < result.size(); ++index) {
    result[result.size() - index - 1U] = digits[fingerprint & 0xfU];
    fingerprint >>= 4U;
  }
  return result;
}

double committed_convective_rate(const runtime::MpiContext &mpi,
                                 const mesh::MeshTopology &topology,
                                 const mesh::MeshGeometry &geometry,
                                 const immersed::ImmersedDomain &domain,
                                 const flow::FlowLayerValues &committed) {
  HUNDUN_CHECK(committed.density.size() == topology.owned_cell_count());
  HUNDUN_CHECK(committed.face_mass_flux.size() == topology.local_face_count());
  std::vector<double> flux_sum(topology.owned_cell_count(), 0.0);
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const double magnitude = std::abs(committed.face_mass_flux[face]);
    HUNDUN_CHECK(std::isfinite(magnitude));
    const auto owner = topology.owner(face);
    if (topology.cell_ownership(owner) == mesh::EntityOwnership::owned &&
        domain.region(owner) == immersed::CellRegion::fluid)
      flux_sum[owner] += magnitude;
    const auto neighbour = topology.neighbour(face);
    if (!neighbour.has_value() || topology.periodic_pair(face).has_value() ||
        topology.cell_ownership(*neighbour) != mesh::EntityOwnership::owned ||
        domain.region(*neighbour) != immersed::CellRegion::fluid)
      continue;
    flux_sum[*neighbour] += magnitude;
  }
  double maximum_rate = 0.0;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const double density = committed.density[cell];
    const double volume = geometry.cell_volume_m3(cell);
    HUNDUN_CHECK(density > 0.0 && std::isfinite(density));
    HUNDUN_CHECK(volume > 0.0 && std::isfinite(volume));
    const double rate = flux_sum[cell] / (2.0 * density * volume);
    HUNDUN_CHECK(rate >= 0.0 && std::isfinite(rate));
    maximum_rate = std::max(maximum_rate, rate);
  }
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &maximum_rate, 1, MPI_DOUBLE, MPI_MAX,
                          mpi.comm()),
            "MPI_Allreduce(engineering convection rate)");
  HUNDUN_CHECK(maximum_rate >= 0.0 && std::isfinite(maximum_rate));
  return maximum_rate;
}

struct CommittedObservation final {
  double velocity_max{};
  double pressure_rms{};
  double pressure_max{};
  double pressure_parity{};
};

CommittedObservation
observe_committed_state(const runtime::MpiContext &mpi,
                        const mesh::MeshTopology &topology,
                        const immersed::ImmersedDomain &domain,
                        const flow::FlowLayerValues &committed) {
  HUNDUN_CHECK(committed.velocity.size() == topology.owned_cell_count() * 3U);
  HUNDUN_CHECK(committed.mechanical_pressure.size() ==
               topology.owned_cell_count());
  double velocity_max = 0.0;
  double pressure_sum_squared = 0.0;
  double pressure_max = 0.0;
  double parity_sum = 0.0;
  double active_count = 0.0;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    double speed_squared = 0.0;
    for (std::size_t component = 0U; component < 3U; ++component) {
      const double value = committed.velocity[cell * 3U + component];
      HUNDUN_CHECK(std::isfinite(value));
      speed_squared += value * value;
    }
    const double pressure = committed.mechanical_pressure[cell];
    HUNDUN_CHECK(std::isfinite(pressure));
    const auto logical = topology.global_cell(cell);
    const double parity =
        ((logical.x + logical.y + logical.z) & 1) == 0 ? 1.0 : -1.0;
    velocity_max = std::max(velocity_max, std::sqrt(speed_squared));
    pressure_sum_squared += pressure * pressure;
    pressure_max = std::max(pressure_max, std::abs(pressure));
    parity_sum += parity * pressure;
    active_count += 1.0;
  }
  std::array<double, 3> sums{pressure_sum_squared, parity_sum, active_count};
  sum_in_place(mpi, sums.data(), sums.size());
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &velocity_max, 1, MPI_DOUBLE, MPI_MAX,
                          mpi.comm()),
            "MPI_Allreduce(engineering velocity maximum)");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &pressure_max, 1, MPI_DOUBLE, MPI_MAX,
                          mpi.comm()),
            "MPI_Allreduce(engineering pressure maximum)");
  HUNDUN_CHECK(sums[0] >= 0.0 && std::isfinite(sums[0]));
  HUNDUN_CHECK(sums[2] > 0.0 && std::isfinite(sums[2]));
  const double pressure_rms = std::sqrt(sums[0] / sums[2]);
  const double parity_denominator = std::sqrt(sums[0] * sums[2]);
  const double pressure_parity =
      parity_denominator > 0.0 ? std::abs(sums[1]) / parity_denominator : 0.0;
  HUNDUN_CHECK(std::isfinite(velocity_max));
  HUNDUN_CHECK(std::isfinite(pressure_rms));
  HUNDUN_CHECK(std::isfinite(pressure_max));
  HUNDUN_CHECK(pressure_parity >= 0.0 && std::isfinite(pressure_parity));
  return {velocity_max, pressure_rms, pressure_max, pressure_parity};
}

test::stage3::BodySpec engineering_body(EngineeringBody kind) {
  if (kind == EngineeringBody::cylinder) {
    return {test::stage3::BodyKind::finite_cylinder,
            {4.0, 0.0, 0.0},
            0.5,
            4.0,
            {0.0, 0.0, 1.0},
            {}};
  }
  return {test::stage3::BodyKind::sphere, {3.0, 0.0, 0.0}, 0.5, 0.0, {}, {}};
}

config::FlowCaseConfig flow_config(const EngineeringCase &definition,
                                   int ranks) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = "stage3-task11-laminar-engineering";
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::constant;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = definition.process_grid;
  result.mesh.cells = definition.cells;
  result.mesh.origin_m = definition.origin_m;
  result.mesh.length_m = definition.length_m;
  result.physics.rho_ref_kg_per_m3 = 1.0;
  result.physics.dynamic_viscosity_pa_s = 1.0 / definition.reynolds;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t patch = 0U; patch < names.size(); ++patch) {
    result.boundaries[patch].patch = names[patch];
    result.boundaries[patch].type = config::BoundaryType::symmetry;
  }
  auto &inlet = result.boundaries[0];
  inlet.type = config::BoundaryType::velocity_inlet;
  inlet.velocity_m_per_s = runtime::Real3{1.0, 0.0, 0.0};
  inlet.thermal_authority = config::InletThermalAuthority::enthalpy;
  inlet.enthalpy_J_per_kg = 1.0;
  inlet.scalar_values = std::vector<config::InletScalarValue>{};
  auto &outlet = result.boundaries[1];
  outlet.type = config::BoundaryType::pressure_outlet;
  outlet.pressure_perturbation_pa = 0.0;
  return result;
}

runtime::FieldDescriptor cell_descriptor(const char *name,
                                         std::uint32_t components) {
  return {name,
          "1",
          "stage3_task11_engineering",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          4,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_descriptor(const char *name,
                                         std::uint32_t components) {
  return {name,
          "1",
          "stage3_task11_engineering",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

std::string collective_surface_path(
    const runtime::MpiContext &mpi,
    std::optional<test::Stage3TemporaryDirectory> &root_directory,
    const test::stage3::ManufacturedSurface &surface) {
  std::string path_text;
  if (mpi.rank() == 0) {
    root_directory.emplace("task11-laminar-engineering");
    const auto path = root_directory->path() / "body.stl";
    test::write_text(path, test::ascii_stl(surface.triangles, "body"));
    path_text = path.string();
  }
  std::uint64_t size = path_text.size();
  check_mpi(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()),
            "MPI_Bcast(engineering path size)");
  HUNDUN_CHECK(size <=
               static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
  path_text.resize(static_cast<std::size_t>(size));
  check_mpi(MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()),
                      MPI_BYTE, 0, mpi.comm()),
            "MPI_Bcast(engineering path)");
  return path_text;
}

double recirculation_length(const runtime::MpiContext &mpi,
                            const mesh::MeshTopology &topology,
                            const immersed::ImmersedDomain &domain,
                            const flow::FlowLayerValues &values,
                            const EngineeringCase &definition) {
  HUNDUN_CHECK(definition.body == EngineeringBody::cylinder);
  const int nx = definition.cells.x;
  const int ny = definition.cells.y;
  const int nz = definition.cells.z;
  const double dx = definition.length_m.x / static_cast<double>(nx);
  const double dy = definition.length_m.y / static_cast<double>(ny);
  const double dz = definition.length_m.z / static_cast<double>(nz);
  const auto lower_index = [](double origin, double spacing, int cells) {
    const double coordinate = (0.0 - origin) / spacing - 0.5;
    return std::clamp(static_cast<int>(std::floor(coordinate)), 0, cells - 2);
  };
  const int y0 = lower_index(definition.origin_m.y, dy, ny);
  const int z0 = lower_index(definition.origin_m.z, dz, nz);
  std::vector<double> velocity(static_cast<std::size_t>(nx), 0.0);
  std::vector<double> count(static_cast<std::size_t>(nx), 0.0);
  const auto owned = topology.owned_global_box();
  for (int i = 0; i < nx; ++i)
    for (const int j : {y0, y0 + 1})
      for (const int k : {z0, z0 + 1}) {
        const runtime::Int3 logical{i, j, k};
        if (logical.x < owned.begin.x || logical.x >= owned.end.x ||
            logical.y < owned.begin.y || logical.y >= owned.end.y ||
            logical.z < owned.begin.z || logical.z >= owned.end.z)
          continue;
        const auto local =
            topology.find_local_cell(topology.global_cell_id(logical));
        HUNDUN_CHECK(local.has_value());
        if (domain.region(*local) != immersed::CellRegion::fluid)
          continue;
        velocity[static_cast<std::size_t>(i)] += values.velocity[*local * 3U];
        count[static_cast<std::size_t>(i)] += 1.0;
      }
  sum_in_place(mpi, velocity.data(), velocity.size());
  sum_in_place(mpi, count.data(), count.size());
  for (int i = 0; i < nx; ++i)
    if (count[static_cast<std::size_t>(i)] > 0.0)
      velocity[static_cast<std::size_t>(i)] /=
          count[static_cast<std::size_t>(i)];

  const double downstream_surface = 4.5;
  bool negative_interval = false;
  double zero_position = 0.0;
  for (int i = 0; i + 1 < nx; ++i) {
    const double x =
        definition.origin_m.x + (static_cast<double>(i) + 0.5) * dx;
    const double next_x = x + dx;
    if (next_x <= downstream_surface ||
        count[static_cast<std::size_t>(i)] == 0.0 ||
        count[static_cast<std::size_t>(i + 1)] == 0.0)
      continue;
    const double left = velocity[static_cast<std::size_t>(i)];
    const double right = velocity[static_cast<std::size_t>(i + 1)];
    if (left < 0.0)
      negative_interval = true;
    if (negative_interval && left <= 0.0 && right > 0.0) {
      const double fraction = -left / (right - left);
      zero_position = x + fraction * dx;
      break;
    }
  }
  HUNDUN_CHECK(negative_interval);
  HUNDUN_CHECK(zero_position > downstream_surface);
  return zero_position - downstream_surface;
}

double unique_separation_crossing(const std::vector<double> &averaged) {
  HUNDUN_CHECK(averaged.size() == 180U);
  double maximum = 0.0;
  for (const double value : averaged)
    if (std::isfinite(value))
      maximum = std::max(maximum, std::abs(value));
  HUNDUN_CHECK(maximum > 0.0 && std::isfinite(maximum));
  const double noise =
      4096.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, maximum);
  int reference_sign = 0;
  std::size_t previous_bin = 0U;
  double previous_value = 0.0;
  std::vector<double> crossings;
  for (std::size_t bin = 5U; bin < averaged.size() - 5U; ++bin) {
    const double value = averaged[bin];
    if (!std::isfinite(value) || std::abs(value) <= noise)
      continue;
    const int sign = value > 0.0 ? 1 : -1;
    if (reference_sign == 0) {
      reference_sign = sign;
      previous_bin = bin;
      previous_value = value;
      continue;
    }
    const int previous_sign = previous_value > 0.0 ? 1 : -1;
    if (sign != previous_sign) {
      const double fraction = std::abs(previous_value) /
                              (std::abs(previous_value) + std::abs(value));
      crossings.push_back(static_cast<double>(previous_bin) + 0.5 +
                          fraction * static_cast<double>(bin - previous_bin));
    }
    previous_bin = bin;
    previous_value = value;
  }
  HUNDUN_CHECK(reference_sign != 0);
  HUNDUN_CHECK(crossings.size() == 1U);
  return crossings.front();
}

bool separation_oracle_is_mutation_sensitive() {
  std::vector<double> one_crossing(180U);
  for (std::size_t bin = 0U; bin < one_crossing.size(); ++bin)
    one_crossing[bin] = 90.0 - (static_cast<double>(bin) + 0.5);
  if (std::abs(unique_separation_crossing(one_crossing) - 90.0) > 1.0e-12)
    return false;
  const auto rejects = [](const std::vector<double> &values) {
    try {
      static_cast<void>(unique_separation_crossing(values));
    } catch (const std::exception &) {
      return true;
    }
    return false;
  };
  std::vector<double> no_crossing(180U, 1.0);
  auto two_crossings = one_crossing;
  for (std::size_t bin = 120U; bin < two_crossings.size(); ++bin)
    two_crossings[bin] = std::abs(two_crossings[bin]);
  return rejects(no_crossing) && rejects(two_crossings);
}

std::array<double, 2>
separation_angles(const runtime::MpiContext &mpi,
                  const immersed::WallQuadraturePlan &wall_plan,
                  const runtime::FieldView<const double> &velocity,
                  const EngineeringCase &definition, double mu) {
  HUNDUN_CHECK(definition.body == EngineeringBody::cylinder);
  constexpr std::size_t bins = 180U;
  std::array<std::vector<double>, 4> reduced{
      std::vector<double>(bins, 0.0), std::vector<double>(bins, 0.0),
      std::vector<double>(bins, 0.0), std::vector<double>(bins, 0.0)};
  const double dz =
      definition.length_m.z / static_cast<double>(definition.cells.z);
  const auto component = [](runtime::Real3 value, std::size_t direction) {
    return direction == 0U ? value.x : direction == 1U ? value.y : value.z;
  };
  for (const auto &point : wall_plan.local_points()) {
    if (std::abs(point.position_m.z) > 1.5 * dz)
      continue;
    const double dx = point.position_m.x - 4.0;
    const double dy = point.position_m.y;
    const double signed_angle = std::atan2(dy, -dx);
    const double absolute_angle = std::abs(signed_angle);
    if (!(absolute_angle > 0.0) || !(absolute_angle < kPi))
      continue;
    const auto bin = std::min<std::size_t>(
        bins - 1U, static_cast<std::size_t>(absolute_angle / kPi *
                                            static_cast<double>(bins)));
    std::array<runtime::Real3, 3> gradient{};
    for (std::size_t row = 0U; row < 3U; ++row)
      gradient[row] = immersed::detail::gradient_with_origin_constraint(
          point.reconstruction, point.position_m, velocity, row, 0.0);
    const double divergence = gradient[0].x + gradient[1].y + gradient[2].z;
    std::array<double, 3> traction{};
    for (std::size_t row = 0U; row < 3U; ++row)
      for (std::size_t direction = 0U; direction < 3U; ++direction) {
        double stress = mu * (component(gradient[row], direction) +
                              component(gradient[direction], row));
        if (row == direction)
          stress -= mu * (2.0 / 3.0) * divergence;
        stress *= component(point.solid_to_fluid_normal, direction);
        traction[row] += stress;
      }
    const runtime::Real3 tangent{std::sin(signed_angle), std::cos(signed_angle),
                                 0.0};
    const double shear = traction[0] * tangent.x + traction[1] * tangent.y;
    HUNDUN_CHECK(std::isfinite(shear));
    const std::size_t side = signed_angle > 0.0 ? 0U : 1U;
    reduced[side * 2U][bin] += shear * point.weight_m2;
    reduced[side * 2U + 1U][bin] += point.weight_m2;
  }
  for (auto &values : reduced)
    sum_in_place(mpi, values.data(), values.size());

  std::array<double, 2> result{};
  for (std::size_t side = 0U; side < 2U; ++side) {
    std::vector<double> averaged(bins,
                                 std::numeric_limits<double>::quiet_NaN());
    for (std::size_t bin = 0U; bin < bins; ++bin) {
      const double weight = reduced[side * 2U + 1U][bin];
      if (weight > 0.0)
        averaged[bin] = reduced[side * 2U][bin] / weight;
    }
    result[side] = unique_separation_crossing(averaged);
  }
  return result;
}

EngineeringResult run_engineering(const runtime::MpiContext &mpi,
                                  const EngineeringCase &definition) {
  HUNDUN_CHECK(definition.process_grid.x * definition.process_grid.y *
                   definition.process_grid.z ==
               mpi.size());
  HUNDUN_CHECK(definition.reynolds > 0.0);
  HUNDUN_CHECK(definition.final_time_s > 0.0);
  HUNDUN_CHECK(definition.statistics_window_s > 0.0 &&
               definition.statistics_window_s <= definition.final_time_s);
  const auto body = engineering_body(definition.body);
  const double h_max = std::max(
      {definition.length_m.x / static_cast<double>(definition.cells.x),
       definition.length_m.y / static_cast<double>(definition.cells.y),
       definition.length_m.z / static_cast<double>(definition.cells.z)});
  const auto surface_fixture =
      test::stage3::make_manufactured_surface(body, h_max);
  std::optional<test::Stage3TemporaryDirectory> root_directory;
  const auto surface_path =
      collective_surface_path(mpi, root_directory, surface_fixture);

  const auto decomposition = runtime::StructuredDecomposition::create(
      mpi, definition.cells, {false, false, false},
      runtime::DecompositionOptions{definition.process_grid});
  const mesh::MeshTopology topology(decomposition);
  const mesh::MeshGeometry geometry(
      topology,
      mesh::UniformBoxMapping(definition.origin_m, definition.length_m));
  const auto boundaries = boundary::BoundaryRegistry::create(
      flow_config(definition, mpi.size()), topology);
  const auto surface = immersed::ImmersedSurface::load_collective(
      std::filesystem::path(surface_path), 1.0, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  const auto ghost_plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  std::uint64_t local_link_count =
      static_cast<std::uint64_t>(domain.links().size());
  std::uint64_t global_link_count = 0U;
  check_mpi(MPI_Allreduce(&local_link_count, &global_link_count, 1,
                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(engineering link count)");
  std::size_t maximum_ghost_donors = 0U;
  for (immersed::ImmersedLinkId link = 0U; link < global_link_count; ++link)
    maximum_ghost_donors =
        std::max(maximum_ghost_donors,
                 ghost_plan.velocity_constraint(link, 0U).donors.size());
  if (mpi.rank() == 0)
    std::cerr << "laminar_engineering ghost_links=" << global_link_count
              << " maximum_ghost_donors=" << maximum_ghost_donors << '\n';
  const auto wall_plan = immersed::WallQuadraturePlan::create(
      surface, query, domain, topology, geometry, mpi);
  HUNDUN_CHECK(ghost_plan.maximum_halo_reach() <= 4U);

  double local_h_min = std::numeric_limits<double>::infinity();
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell)
    if (domain.region(cell) == immersed::CellRegion::fluid)
      local_h_min =
          std::min(local_h_min, std::cbrt(geometry.cell_volume_m3(cell)));
  double inverse_h_min = std::isfinite(local_h_min) ? 1.0 / local_h_min : 0.0;
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &inverse_h_min, 1, MPI_DOUBLE, MPI_MAX,
                          mpi.comm()),
            "MPI_Allreduce(engineering h-min)");
  HUNDUN_CHECK(inverse_h_min > 0.0 && std::isfinite(inverse_h_min));
  const double h_min = 1.0 / inverse_h_min;
  const double stable_dt = 0.25 * h_min;
  const auto nominal_step_count = static_cast<std::uint64_t>(
      std::ceil(definition.final_time_s / stable_dt));
  HUNDUN_CHECK(nominal_step_count > 0U);
  const double maximum_dt =
      definition.final_time_s / static_cast<double>(nominal_step_count);
  HUNDUN_CHECK(maximum_dt / h_min <=
               0.25 * (1.0 + 64.0 * std::numeric_limits<double>::epsilon()));

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_descriptor("rho", 1U));
  fields.velocity = registry.declare_field(cell_descriptor("velocity", 3U));
  fields.mechanical_pressure =
      registry.declare_field(cell_descriptor("pi", 1U));
  fields.face_velocity =
      registry.declare_field(face_descriptor("face_velocity", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  registry.freeze();

  flow::FlowLayerValues initial;
  initial.density.resize(topology.owned_cell_count(), 0.0);
  initial.velocity.resize(topology.owned_cell_count() * 3U, 0.0);
  initial.mechanical_pressure.resize(topology.owned_cell_count(), 0.0);
  initial.face_velocity.resize(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.resize(topology.local_face_count(), 0.0);
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    initial.density[cell] = 1.0;
    initial.velocity[cell * 3U] = 1.0;
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const bool owner_active =
        domain.region(topology.owner(face)) == immersed::CellRegion::fluid;
    const auto neighbour = topology.neighbour(face);
    const bool neighbour_active =
        neighbour.has_value() &&
        domain.region(*neighbour) == immersed::CellRegion::fluid;
    if (!owner_active || (neighbour.has_value() && !neighbour_active))
      continue;
    initial.face_velocity[face * 3U] = 1.0;
    const auto area = geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
    initial.face_mass_flux[face] = area.x;
  }
  {
    const auto baseline =
        observe_committed_state(mpi, topology, domain, initial);
    HUNDUN_CHECK(baseline.velocity_max == 1.0);
    HUNDUN_CHECK(baseline.pressure_rms == 0.0);
    HUNDUN_CHECK(baseline.pressure_max == 0.0);
    HUNDUN_CHECK(baseline.pressure_parity == 0.0);
    auto mutation = initial;
    const auto origin_cell = topology.find_local_cell(0U);
    if (origin_cell.has_value() &&
        topology.cell_ownership(*origin_cell) == mesh::EntityOwnership::owned)
      mutation.mechanical_pressure[*origin_cell] = 1.0;
    const auto mutated =
        observe_committed_state(mpi, topology, domain, mutation);
    HUNDUN_CHECK(mutated.velocity_max == baseline.velocity_max);
    HUNDUN_CHECK(mutated.pressure_rms > baseline.pressure_rms);
    HUNDUN_CHECK(mutated.pressure_max == 1.0);
    HUNDUN_CHECK(mutated.pressure_parity > baseline.pressure_parity);
  }

  auto state = flow::FlowState::create(
      registry, {decomposition.local_extent(), topology.local_face_count()},
      fields,
      {0U, 0.0, maximum_dt, 0.0, flow::MomentumTimeOrder::backward_euler});
  state.seed_accepted_layers(initial, initial);
  {
    const auto baseline = state.snapshot(flow::FlowLayer::committed);
    auto mutation = baseline;
    for (double &flux : mutation.face_mass_flux)
      flux *= 2.0;
    const double baseline_rate =
        committed_convective_rate(mpi, topology, geometry, domain, baseline);
    const double mutation_rate =
        committed_convective_rate(mpi, topology, geometry, domain, mutation);
    HUNDUN_CHECK(baseline_rate > 0.0);
    HUNDUN_CHECK_NEAR(mutation_rate, 2.0 * baseline_rate,
                      64.0 * std::numeric_limits<double>::epsilon() *
                          std::max(1.0, mutation_rate));
  }
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
  immersed::LocalFlowPatternTransform transform;
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, &domain, &ghost_plan,
      &wall_plan, &transform, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);
  const auto pressure_algebra =
      flow::test::ImmersedFlowTestAccess::active_pressure_algebra_probe(immersed_flow);
  HUNDUN_CHECK(std::isfinite(pressure_algebra[0]));
  HUNDUN_CHECK(std::isfinite(pressure_algebra[1]));
  HUNDUN_CHECK(std::isfinite(pressure_algebra[2]));
  HUNDUN_CHECK(pressure_algebra[0] <=
               512.0 * std::numeric_limits<double>::epsilon());
  HUNDUN_CHECK(pressure_algebra[1] > 0.0);
  HUNDUN_CHECK(pressure_algebra[2] > 0.0);
  if (mpi.rank() == 0)
    std::cerr << "laminar_engineering pressure_symmetry=" << pressure_algebra[0]
              << " rayleigh=" << pressure_algebra[1] << ','
              << pressure_algebra[2] << '\n';
  const double mu = 1.0 / definition.reynolds;
  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                    1.0,
                                    mu,
                                    std::nullopt,
                                    std::nullopt,
                                    std::nullopt};
  linear::SolveControl control;
  control.max_iterations = 5000U;
  control.residual_recompute_interval = 20U;

  std::vector<double> total_drag;
  std::vector<double> pressure_drag;
  std::vector<double> viscous_drag;
  std::vector<double> lateral_y;
  std::vector<double> lateral_z;
  std::vector<double> recirculation;
  const double statistics_begin =
      definition.final_time_s - definition.statistics_window_s;
  flow::ForceAttemptReport final_force_report{};
  double next_dt = maximum_dt;
  std::uint64_t accepted_steps = 0U;
  double minimum_accepted_dt = std::numeric_limits<double>::infinity();
  double maximum_accepted_dt = 0.0;
  while (state.metadata().time_s < definition.final_time_s) {
    const auto metadata_before = state.metadata();
    const double remaining = definition.final_time_s - metadata_before.time_s;
    const double finish_tolerance = 64.0 *
                                    std::numeric_limits<double>::epsilon() *
                                    std::max(1.0, definition.final_time_s);
    if (remaining <= finish_tolerance)
      break;
    const auto committed_before = state.snapshot(flow::FlowLayer::committed);
    const double convective_rate = committed_convective_rate(
        mpi, topology, geometry, domain, committed_before);
    const double cfl_dt = convective_rate > 0.0
                              ? 0.25 / convective_rate
                              : std::numeric_limits<double>::infinity();
    double attempt_dt = std::min({next_dt, remaining, cfl_dt});
    HUNDUN_CHECK(attempt_dt >= maximum_dt / 256.0);
    flow::ImmersedFlowStepAttemptReport attempt{};
    flow::MomentumTimeStencil accepted_stencil{};
    std::uint32_t retry = 0U;
    for (;;) {
      HUNDUN_CHECK(attempt_dt > 0.0 && std::isfinite(attempt_dt));
      HUNDUN_CHECK(attempt_dt * convective_rate <=
                   0.25 *
                       (1.0 + 64.0 * std::numeric_limits<double>::epsilon()));
      const auto metadata = state.metadata();
      const double ratio =
          metadata.step == 0U ? 0.0 : attempt_dt / metadata.dt_s;
      const auto order = metadata.step != 0U && ratio >= 0.5 && ratio <= 2.0
                             ? flow::MomentumTimeOrder::bdf2
                             : flow::MomentumTimeOrder::backward_euler;
      const auto stencil =
          flow::make_momentum_time_stencil(order, attempt_dt, metadata.dt_s);
      attempt = immersed_flow.attempt(state, physics, stencil, control, control);
      const auto &trial = std::get<flow::StepAttemptReport>(attempt.base);
      if (trial.disposition == flow::StepAttemptDisposition::committed) {
        accepted_stencil = stencil;
        break;
      }
      if (mpi.rank() == 0)
        std::cerr << "laminar_engineering retry step=" << metadata.step
                  << " attempt=" << retry
                  << " disposition=" << static_cast<int>(trial.disposition)
                  << " reason=" << static_cast<int>(trial.reason)
                  << " dt=" << attempt_dt
                  << " order=" << static_cast<int>(order);
      if (mpi.rank() == 0 && trial.final_backflow_evidence.has_value()) {
        const auto &evidence = *trial.final_backflow_evidence;
        std::cerr << " backflow_minimum="
                  << evidence.minimum_outward_mass_flux_kg_per_s
                  << " backflow_face=" << evidence.global_face_id;
        const auto local_face =
            topology.find_local_face(evidence.global_face_id);
        if (local_face.has_value()) {
          const auto center = geometry.face_center_m(*local_face);
          std::cerr << " backflow_center=(" << center.x << ',' << center.y
                    << ',' << center.z << ')';
        }
      }
      if (mpi.rank() == 0)
        std::cerr << '\n';
      HUNDUN_CHECK(trial.disposition ==
                   flow::StepAttemptDisposition::recoverable_failure);
      HUNDUN_CHECK(retry < 8U);
      attempt_dt *= 0.5;
      HUNDUN_CHECK(attempt_dt >= maximum_dt / 256.0);
      ++retry;
    }
    const auto &base = std::get<flow::StepAttemptReport>(attempt.base);
    if (base.disposition != flow::StepAttemptDisposition::committed &&
        mpi.rank() == 0) {
      std::cerr << "laminar_engineering failed step=" << accepted_steps
                << " body=" << static_cast<int>(definition.body)
                << " disposition=" << static_cast<int>(base.disposition)
                << " reason=" << static_cast<int>(base.reason)
                << " lowest_failing_rank=" << base.lowest_failing_rank
                << " correctors=" << base.pressure_corrector_count
                << " continuity=" << base.final_continuity_normalized_l2
                << " pressure_residual=" << base.final_pressure_residual_l2
                << " momentum_residuals="
                << base.final_momentum_normalized_l2[0] << ','
                << base.final_momentum_normalized_l2[1] << ','
                << base.final_momentum_normalized_l2[2];
      for (std::size_t component = 0U; component < 3U; ++component) {
        const auto &solve = base.momentum.components[component];
        std::cerr << " momentum[" << component
                  << "]={reason=" << static_cast<int>(solve.reason)
                  << ",iterations=" << solve.iterations
                  << ",initial=" << solve.initial_residual
                  << ",recursive=" << solve.recursive_residual
                  << ",final=" << solve.final_residual << '}';
      }
      for (std::size_t correction = 0U; correction < 2U; ++correction) {
        const auto &solve = base.pressure[correction];
        std::cerr << " pressure[" << correction
                  << "]={reason=" << static_cast<int>(solve.reason)
                  << ",iterations=" << solve.iterations
                  << ",initial=" << solve.initial_residual
                  << ",recursive=" << solve.recursive_residual
                  << ",final=" << solve.final_residual << '}';
      }
      if (base.final_backflow_evidence.has_value()) {
        const auto &evidence = *base.final_backflow_evidence;
        std::cerr << " backflow={patch=" << evidence.patch_id
                  << ",step=" << evidence.step << ",time=" << evidence.time_s
                  << ",minimum=" << evidence.minimum_outward_mass_flux_kg_per_s
                  << ",face=" << evidence.global_face_id
                  << ",rank=" << evidence.lowest_failing_rank << '}';
        const auto local_face =
            topology.find_local_face(evidence.global_face_id);
        if (local_face.has_value()) {
          const auto center = geometry.face_center_m(*local_face);
          std::cerr << " backflow_center=(" << center.x << ',' << center.y
                    << ',' << center.z << ')';
        }
      }
      std::cerr << '\n';
    }
    HUNDUN_CHECK(base.disposition == flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(base.reason == flow::StepFailureReason::none);
    HUNDUN_CHECK(base.pressure_corrector_count == 2U);
    HUNDUN_CHECK(attempt.force.has_value());
    final_force_report = *attempt.force;
    const auto &accepted_force = attempt.force->surface_traction;
    for (const auto &part :
         {accepted_force.pressure_N, accepted_force.viscous_N,
          accepted_force.total_N}) {
      HUNDUN_CHECK(std::isfinite(part.x));
      HUNDUN_CHECK(std::isfinite(part.y));
      HUNDUN_CHECK(std::isfinite(part.z));
    }
    HUNDUN_CHECK(std::abs(accepted_force.total_N.x) <=
                 definition.fast_force_bound_N);
    ++accepted_steps;
    minimum_accepted_dt = std::min(minimum_accepted_dt, accepted_stencil.dt_s);
    maximum_accepted_dt = std::max(maximum_accepted_dt, accepted_stencil.dt_s);
    const bool inexpensive_solves =
        std::all_of(base.momentum.components.begin(),
                    base.momentum.components.end(),
                    [&](const auto &solve) {
                      return solve.iterations <= control.max_iterations / 2U;
                    }) &&
        std::all_of(base.pressure.begin(), base.pressure.end(),
                    [&](const auto &solve) {
                      return solve.iterations <= control.max_iterations / 2U;
                    });
    next_dt = inexpensive_solves
                  ? std::min(maximum_dt, 1.25 * accepted_stencil.dt_s)
                  : accepted_stencil.dt_s;
    const double time = state.metadata().time_s;
    const auto reporting_interval =
        std::max<std::uint64_t>(1U, nominal_step_count / 10U);
    const bool report_step = accepted_steps % reporting_interval == 0U ||
                             definition.final_time_s - time <= finish_tolerance;
    if (report_step) {
      const auto accepted_state = state.snapshot(flow::FlowLayer::committed);
      const auto observation =
          observe_committed_state(mpi, topology, domain, accepted_state);
      const double accepted_rate = committed_convective_rate(
          mpi, topology, geometry, domain, accepted_state);
      if (mpi.rank() == 0)
        std::cerr << "laminar_engineering step=" << accepted_steps
                  << " time=" << time << " dt=" << accepted_stencil.dt_s
                  << " force_x=" << accepted_force.total_N.x
                  << " pressure_force_x=" << accepted_force.pressure_N.x
                  << " viscous_force_x=" << accepted_force.viscous_N.x
                  << " velocity_max=" << observation.velocity_max
                  << " convection_rate=" << accepted_rate
                  << " pressure_rms=" << observation.pressure_rms
                  << " pressure_max=" << observation.pressure_max
                  << " pressure_parity=" << observation.pressure_parity
                  << " continuity=" << base.final_continuity_normalized_l2
                  << " pressure_residual=" << base.final_pressure_residual_l2
                  << '\n';
    }
    if (time + 64.0 * std::numeric_limits<double>::epsilon() < statistics_begin)
      continue;
    const double denominator = definition.body == EngineeringBody::cylinder
                                   ? 0.5 * 1.0 * 1.0 * 1.0 * 4.0
                                   : 0.5 * 1.0 * 1.0 * kPi * 0.25;
    const auto &force = attempt.force->surface_traction;
    total_drag.push_back(force.total_N.x / denominator);
    pressure_drag.push_back(force.pressure_N.x / denominator);
    viscous_drag.push_back(force.viscous_N.x / denominator);
    lateral_y.push_back(force.total_N.y / denominator);
    lateral_z.push_back(force.total_N.z / denominator);
    if (definition.body == EngineeringBody::cylinder &&
        definition.collect_cylinder_diagnostics) {
      const auto committed = state.snapshot(flow::FlowLayer::committed);
      recirculation.push_back(
          recirculation_length(mpi, topology, domain, committed, definition));
    }
  }
  HUNDUN_CHECK(!total_drag.empty());
  HUNDUN_CHECK(total_drag.size() == pressure_drag.size());
  HUNDUN_CHECK(total_drag.size() == viscous_drag.size());

  EngineeringResult result;
  result.total_drag_mean = mean(total_drag);
  result.pressure_drag_mean = mean(pressure_drag);
  result.viscous_drag_mean = mean(viscous_drag);
  result.lateral_y_mean = mean(lateral_y);
  result.lateral_z_mean = mean(lateral_z);
  result.total_drag_cv = coefficient_of_variation(total_drag);
  result.pressure_drag_cv = coefficient_of_variation(pressure_drag);
  result.viscous_drag_cv = coefficient_of_variation(viscous_drag);
  if (!recirculation.empty()) {
    result.recirculation_length_mean_m = mean(recirculation);
    result.recirculation_length_cv = coefficient_of_variation(recirculation);
  }
  result.classification_fingerprint = domain.classification_fingerprint();
  result.ghost_plan_fingerprint = ghost_plan.fingerprint();
  result.wall_plan_fingerprint = wall_plan.fingerprint();
  result.accepted_steps = accepted_steps;
  result.minimum_accepted_dt_s = minimum_accepted_dt;
  result.maximum_accepted_dt_s = maximum_accepted_dt;
  HUNDUN_CHECK(result.accepted_steps > 0U);
  HUNDUN_CHECK(result.minimum_accepted_dt_s > 0.0 &&
               std::isfinite(result.minimum_accepted_dt_s));
  HUNDUN_CHECK(result.maximum_accepted_dt_s >= result.minimum_accepted_dt_s &&
               std::isfinite(result.maximum_accepted_dt_s));
  result.final_force = final_force_report.surface_traction;
  if (definition.body == EngineeringBody::cylinder &&
      definition.collect_cylinder_diagnostics) {
    const auto velocity =
        state.layer(flow::FlowLayer::committed).view<double>(fields.velocity);
    const auto angles =
        separation_angles(mpi, wall_plan, velocity, definition, mu);
    result.upper_separation_angle_deg = angles[0];
    result.lower_separation_angle_deg = angles[1];
  }

  const std::size_t global_cells =
      static_cast<std::size_t>(definition.cells.x) *
      static_cast<std::size_t>(definition.cells.y) *
      static_cast<std::size_t>(definition.cells.z);
  result.global_velocity.assign(global_cells * 3U, 0.0);
  std::vector<double> active_mask(global_cells, 0.0);
  const auto committed = state.snapshot(flow::FlowLayer::committed);
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const auto global = static_cast<std::size_t>(topology.global_cell_id(cell));
    active_mask[global] = 1.0;
    for (std::size_t component = 0U; component < 3U; ++component)
      result.global_velocity[global * 3U + component] =
          committed.velocity[cell * 3U + component];
  }
  sum_in_place(mpi, active_mask.data(), active_mask.size());
  sum_in_place(mpi, result.global_velocity.data(),
               result.global_velocity.size());
  for (std::size_t global = 0U; global < global_cells; ++global) {
    HUNDUN_CHECK(active_mask[global] == 0.0 || active_mask[global] == 1.0);
    if (active_mask[global] == 1.0)
      result.global_active_cell_ids.push_back(
          static_cast<mesh::GlobalCellId>(global));
  }
  result.state_fingerprint = engineering_state_fingerprint(
      result.global_active_cell_ids, result.global_velocity);
  const auto exact_copy = result.global_velocity;
  HUNDUN_CHECK(engineering_state_fingerprint(result.global_active_cell_ids,
                                             exact_copy) ==
               result.state_fingerprint);
  auto mutation = exact_copy;
  HUNDUN_CHECK(!result.global_active_cell_ids.empty());
  const auto mutation_offset =
      static_cast<std::size_t>(result.global_active_cell_ids.front()) * 3U;
  mutation[mutation_offset] = std::nextafter(
      mutation[mutation_offset], std::numeric_limits<double>::infinity());
  HUNDUN_CHECK(
      engineering_state_fingerprint(result.global_active_cell_ids, mutation) !=
      result.state_fingerprint);
  if (mpi.rank() == 0)
    std::cerr << "laminar_engineering final accepted_steps="
              << result.accepted_steps
              << " min_dt=" << result.minimum_accepted_dt_s
              << " max_dt=" << result.maximum_accepted_dt_s
              << " state_fingerprint="
              << fingerprint_hex(result.state_fingerprint) << '\n';
  return result;
}

EngineeringCase sphere_case(runtime::Int3 cells, runtime::Int3 process_grid,
                            double final_time_s, double statistics_window_s) {
  return {EngineeringBody::sphere, cells, process_grid, {0.0, -3.0, -3.0},
          {12.0, 6.0, 6.0},        1.0,   final_time_s, statistics_window_s};
}

EngineeringCase cylinder_case(double reynolds, runtime::Int3 process_grid) {
  return {EngineeringBody::cylinder,
          {192, 144, 72},
          process_grid,
          {0.0, -6.0, -3.0},
          {16.0, 12.0, 6.0},
          reynolds,
          40.0,
          10.0};
}

void run_smoke(const runtime::MpiContext &mpi) {
  const auto basic = run_engineering(
      mpi, sphere_case({24, 12, 12}, default_process_grid(mpi.size()), 0.125,
                       0.125));
  HUNDUN_CHECK(std::isfinite(basic.total_drag_mean));
  HUNDUN_CHECK(std::isfinite(basic.pressure_drag_mean));
  HUNDUN_CHECK(std::isfinite(basic.viscous_drag_mean));
  HUNDUN_CHECK(!basic.global_active_cell_ids.empty());
  HUNDUN_CHECK(!basic.global_velocity.empty());
  const auto stability = run_engineering(
      mpi,
      sphere_case({48, 24, 24}, default_process_grid(mpi.size()), 0.5, 0.5));
  HUNDUN_CHECK(std::isfinite(stability.total_drag_mean));
  HUNDUN_CHECK(std::abs(stability.total_drag_mean) < 100.0);
}

void run_cylinder_smoke(const runtime::MpiContext &mpi) {
  const EngineeringCase definition{EngineeringBody::cylinder,
                                   {96, 72, 36},
                                   default_process_grid(mpi.size()),
                                   {0.0, -6.0, -3.0},
                                   {16.0, 12.0, 6.0},
                                   20.0,
                                   8.0,
                                   2.0,
                                   false,
                                   100.0};
  const auto result = run_engineering(mpi, definition);
  HUNDUN_CHECK(std::isfinite(result.total_drag_mean));
  HUNDUN_CHECK(!result.global_active_cell_ids.empty());
  HUNDUN_CHECK(!result.global_velocity.empty());
}

void compare_decomposition_case(const runtime::MpiContext &world,
                                const EngineeringCase &distributed_case,
                                const EngineeringCase &reference_case) {
  const auto distributed = run_engineering(world, distributed_case);
  std::optional<EngineeringResult> reference;
  if (world.rank() == 0) {
    const auto self = runtime::MpiContext::duplicate(MPI_COMM_SELF);
    reference.emplace(run_engineering(self, reference_case));
  }
  int comparison_ok = 1;
  if (world.rank() == 0) {
    try {
      HUNDUN_CHECK(reference.has_value());
      HUNDUN_CHECK(reference->global_active_cell_ids ==
                   distributed.global_active_cell_ids);
      HUNDUN_CHECK(reference->classification_fingerprint ==
                   distributed.classification_fingerprint);
      HUNDUN_CHECK(reference->ghost_plan_fingerprint ==
                   distributed.ghost_plan_fingerprint);
      HUNDUN_CHECK(reference->wall_plan_fingerprint ==
                   distributed.wall_plan_fingerprint);
      test::stage3::require_decomposition_field("engineering velocity",
                                                reference->global_velocity,
                                                distributed.global_velocity);
      test::stage3::require_decomposition_force(
          "engineering force", reference->final_force, distributed.final_force);
    } catch (...) {
      comparison_ok = 0;
    }
  }
  check_mpi(MPI_Bcast(&comparison_ok, 1, MPI_INT, 0, world.comm()),
            "MPI_Bcast(engineering comparison)");
  HUNDUN_CHECK(comparison_ok == 1);
}

void run_decomposition_fast(const runtime::MpiContext &world) {
  const auto world_grid = default_process_grid(world.size());
  compare_decomposition_case(world,
                             {EngineeringBody::cylinder,
                              {24, 18, 24},
                              world_grid,
                              {0.0, -6.0, -3.0},
                              {16.0, 12.0, 6.0},
                              20.0,
                              0.5,
                              0.5,
                              false,
                              100.0},
                             {EngineeringBody::cylinder,
                              {24, 18, 24},
                              {1, 1, 1},
                              {0.0, -6.0, -3.0},
                              {16.0, 12.0, 6.0},
                              20.0,
                              0.5,
                              0.5,
                              false,
                              100.0});
  compare_decomposition_case(world,
                             sphere_case({24, 12, 12}, world_grid, 0.5, 0.5),
                             sphere_case({24, 12, 12}, {1, 1, 1}, 0.5, 0.5));
}

void run_cylinder_acceptance(const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(mpi.size() == 4);
  const auto grid = default_process_grid(mpi.size());
  const auto reynolds_20 = run_engineering(mpi, cylinder_case(20.0, grid));
  const auto reynolds_40 = run_engineering(mpi, cylinder_case(40.0, grid));
  for (const auto *result : {&reynolds_20, &reynolds_40}) {
    HUNDUN_CHECK(result->total_drag_mean > 0.0);
    HUNDUN_CHECK(result->total_drag_cv <= 1.0e-2);
    HUNDUN_CHECK(result->recirculation_length_cv <= 1.0e-2);
    HUNDUN_CHECK(std::abs(result->lateral_y_mean) <=
                 0.02 * result->total_drag_mean);
    HUNDUN_CHECK(std::abs(result->upper_separation_angle_deg -
                          result->lower_separation_angle_deg) <= 5.0);
  }
  HUNDUN_CHECK(reynolds_20.total_drag_mean > reynolds_40.total_drag_mean);
  HUNDUN_CHECK(reynolds_40.recirculation_length_mean_m >=
               reynolds_20.recirculation_length_mean_m + 0.2);
}

void run_sphere_acceptance(const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(mpi.size() == 4);
  const auto result = run_engineering(
      mpi,
      sphere_case({192, 96, 96}, default_process_grid(mpi.size()), 24.0, 6.0));
  HUNDUN_CHECK(result.total_drag_mean > 0.0);
  constexpr double stokes_drag = 24.0;
  HUNDUN_CHECK(result.total_drag_mean / stokes_drag >= 0.90);
  HUNDUN_CHECK(result.total_drag_mean / stokes_drag <= 1.60);
  const double pressure_fraction =
      result.pressure_drag_mean / result.total_drag_mean;
  const double viscous_fraction =
      result.viscous_drag_mean / result.total_drag_mean;
  HUNDUN_CHECK(pressure_fraction >= 0.20 && pressure_fraction <= 0.48);
  HUNDUN_CHECK(viscous_fraction >= 0.52 && viscous_fraction <= 0.80);
  HUNDUN_CHECK(std::abs(result.pressure_drag_mean + result.viscous_drag_mean -
                        result.total_drag_mean) /
                   result.total_drag_mean <=
               5.0e-11);
  HUNDUN_CHECK(std::abs(result.lateral_y_mean) <=
               0.02 * result.total_drag_mean);
  HUNDUN_CHECK(std::abs(result.lateral_z_mean) <=
               0.02 * result.total_drag_mean);
  HUNDUN_CHECK(result.total_drag_cv <= 1.0e-2);
  HUNDUN_CHECK(result.pressure_drag_cv <= 1.0e-2);
  HUNDUN_CHECK(result.viscous_drag_cv <= 1.0e-2);
}

void run(int argc, char **argv, const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(argc == 2);
  const std::string selector = argv[1];
  HUNDUN_CHECK(selector == "smoke" || selector == "cylinder_smoke" ||
               selector == "decomposition_fast" ||
               selector == "cylinder_acceptance" ||
               selector == "sphere_acceptance");
  HUNDUN_CHECK(
      test::stage3::decomposition_equality_oracle_is_mutation_sensitive());
  HUNDUN_CHECK(separation_oracle_is_mutation_sensitive());
  if (selector == "smoke")
    run_smoke(mpi);
  else if (selector == "cylinder_smoke")
    run_cylinder_smoke(mpi);
  else if (selector == "decomposition_fast")
    run_decomposition_fast(mpi);
  else if (selector == "cylinder_acceptance")
    run_cylinder_acceptance(mpi);
  else
    run_sphere_acceptance(mpi);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  const auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] { run(argc, argv, mpi); });
}
