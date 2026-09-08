// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "core_spray_events_detail.hpp"
#include "core_spray_gas_detail.hpp"
#include "models_spray_mechanics_detail.hpp"

#include <mpi.h>

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
using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool near(double actual, double expected, double tolerance = 1.0e-12) {
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::abs(actual - expected) <= tolerance;
}

bool near(Vector3 actual, Vector3 expected, double tolerance = 1.0e-12) {
  return near(actual[0], expected[0], tolerance) &&
         near(actual[1], expected[1], tolerance) &&
         near(actual[2], expected[2], tolerance);
}

CartesianMeshSpec uniform_mesh(Int3 cells = {4, 4, 4}) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {4.0, 4.0, 4.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = cells;
  mesh.minimum_spacing = {1.0e-9, 1.0e-9, 1.0e-9};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {1000000U, UINT64_C(1073741824)};
  return mesh;
}

CartesianMeshSpec tensor_mesh() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::coast_runtime_axes_v1;
  mesh.axes_file = "in-memory-axes.dat";
  mesh.coast_runtime_faces = {
      std::vector<double>{0.0, 0.4, 1.2, 2.5, 4.0},
      std::vector<double>{-2.0, -0.5, 0.25, 2.0},
      std::vector<double>{1.0, 1.1, 1.8, 3.0}};
  mesh.lower = {0.0, -2.0, 1.0};
  mesh.upper = {4.0, 2.0, 3.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {4, 3, 3};
  mesh.minimum_spacing = {1.0e-9, 1.0e-9, 1.0e-9};
  mesh.max_growth_ratio = 10.0;
  mesh.limits = {1000000U, UINT64_C(1073741824)};
  return mesh;
}

bool compile_geometry(const CartesianMeshSpec& spec,
                      CartesianGeometryPlan& geometry,
                      MeshPatch& patch) {
  return static_cast<bool>(CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, spec, GeometryBudget{}, geometry, patch));
}

SprayParcelState parcel(Vector3 position, Vector3 velocity) {
  SprayParcelState value;
  value.id = {11U, 29U};
  value.position_m = position;
  value.velocity_m_per_s = velocity;
  value.droplet_mass_kg = 2.0e-9;
  value.droplet_diameter_m = 1.0e-4;
  value.multiplicity = 5.0;
  value.temperature_k = 300.0;
  value.liquid_material_fingerprint = 7U;
  value.owner_global_cell = 0U;
  value.age_s = 0.25;
  return value;
}

ParcelAccelerationSample zero_acceleration(
    const ParcelKinematicSample&, const void*) noexcept {
  return {true, {0.0, 0.0, 0.0}};
}

ParcelAccelerationSample constant_acceleration(
    const ParcelKinematicSample&, const void* context) noexcept {
  if (context == nullptr) {
    return {};
  }
  return {true, *static_cast<const Vector3*>(context)};
}

ParcelAccelerationSample unavailable_acceleration(
    const ParcelKinematicSample&, const void*) noexcept {
  return {};
}

struct StokesRelaxationContext {
  Vector3 carrier_velocity_m_per_s{};
  double inverse_relaxation_time_per_s{};
};

ParcelAccelerationSample stokes_relaxation(
    const ParcelKinematicSample& sample, const void* context) noexcept {
  if (context == nullptr) return {};
  const auto& model = *static_cast<const StokesRelaxationContext*>(context);
  Vector3 acceleration{};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    acceleration[axis] = model.inverse_relaxation_time_per_s *
                         (model.carrier_velocity_m_per_s[axis] -
                          sample.velocity_m_per_s[axis]);
  }
  return {true, acceleration};
}

bool test_native_gas_sampler() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  if (!compile_geometry(tensor_mesh(), geometry, patch))
    return false;
  const auto cells = patch.cells;
  const std::size_t count = std::size_t(cells.x) * cells.y * cells.z;
  std::array<std::vector<double>, 5> buffers;
  std::array<ConstFieldView, 5> views{};
  for (std::size_t f = 0; f < views.size(); ++f) {
    const unsigned components = f == 2 ? 3 : 1;
    buffers[f].resize(count * components);
    auto &v = views[f];
    v.base = buffers[f].data();
    v.interior = cells;
    v.components = components;
    v.stride_y = cells.x;
    v.stride_z = cells.x * cells.y;
    v.component_stride = count;
    v.field = f;
    v.revision = 9;
    v.storage_identity = 100 + f;
    v.revision_domain = 1000;
    for (int z = 0; z < cells.z; ++z)
      for (int y = 0; y < cells.y; ++y)
        for (int x = 0; x < cells.x; ++x) {
          const double cx = geometry.x().centres().data[x];
          const double cy = geometry.y().centres().data[y];
          const auto i = std::size_t(x) + cells.x * (y + cells.y * z);
          buffers[f][i] = f == 0   ? 10 * cx
                          : f == 1 ? 300000 + 100 * cx
                          : f == 2 ? 2 * cx
                          : f == 3 ? .2 + .01 * cx
                                   : .3 + .01 * cy;
          if (f == 2) {
            buffers[f][count + i] = cy;
            buffers[f][2 * count + i] = 0;
          }
        }
  }
  const std::array<std::size_t, 2> mapping{2, 0};
  hundun::v04::detail::ProductParcelGas sampler;
  const portable::Revision revision{11, 71, 1};
  if (!sampler.configure(geometry, patch, {true, false, false}, 100,
                         {mapping.data(), mapping.size()}, 1) ||
      !sampler.bind(revision, .1, 101325, views[0], views[1], views[2],
                    {views.data() + 3, 2}))
    return false;
  const auto p = parcel({.75, -.2, 1.5}, {0, 0, 0});
  std::array<double, 3> ys{-1, -1, -1};
  const auto gas = sampler.sample(p, .05, ParcelPass::corrector, revision,
                                  ys.data(), ys.size());
  bool passed = expect(
      gas.status == portable::Status::success && gas.revision == revision &&
          gas.composition_fingerprint == 100 && gas.species_count == 3 &&
          near(gas.pressure_pa, 101332.5) &&
          near(gas.enthalpy_j_per_kg, 300075) &&
          near(gas.velocity_m_per_s, {1.5, -.2, 0}) && near(ys[0], .298) &&
          near(ys[2], .2075) && near(ys[1], .4945),
      "native nonuniform PH/U/Y fields feed the full-species parcel query by "
      "name mapping");
  const auto before = ys;
  auto stale = revision;
  ++stale.input_revision;
  const auto rejected = sampler.sample(p, .05, ParcelPass::predictor, stale,
                                       ys.data(), ys.size());
  passed &= expect(
      rejected.status == portable::Status::stale_revision && ys == before,
      "stale accepted field revision cannot publish parcel gas values");
  spray::detail::ParcelLocation location;
  passed &= expect(sampler.locate(p.position_m, revision, location) &&
                       location.global_cell == 17 && location.owner_rank == 0,
                   "actual stretched faces locate parcel ownership");
  auto wrapped = p.position_m;
  wrapped[0] += 4;
  passed &= expect(sampler.locate(wrapped, revision, location) &&
                       location.global_cell == 17,
                   "periodic image resolves to the same physical owner");
  auto probe = p;
  const auto faces = geometry.y().faces();
  probe.position_m[1] = faces.data[0] - .25 * (faces.data[1] - faces.data[0]);
  passed &= expect(sampler.sample(probe, .05, ParcelPass::predictor, revision,
                                  ys.data(), ys.size())
                           .status != portable::Status::success,
                   "strict gas sampler rejects an outside-domain probe");
  passed &=
      bool(sampler.configure(geometry, patch, {true, false, false}, 100,
                             {mapping.data(), mapping.size()}, 1, true)) &&
      bool(sampler.bind(revision, .1, 101325, views[0], views[1], views[2],
                        {views.data() + 3, 2}));
  const auto extended = sampler.sample(probe, .05, ParcelPass::predictor,
                                       revision, ys.data(), ys.size());
  passed &= expect(extended.status == portable::Status::success &&
                       near(ys[0], .3 + .01 * geometry.y().centres().data[0]) &&
                       !sampler.stencil(probe.position_m).succeeded() &&
                       !sampler.locate(probe.position_m, revision, location),
                   "bounded trial continuation is not an accepted owner or "
                   "deposition permission");
  probe.position_m[1] = faces.data[0] - 2 * (faces.data[1] - faces.data[0]);
  passed &= expect(sampler.sample(probe, .05, ParcelPass::corrector, revision,
                                  ys.data(), ys.size())
                           .status != portable::Status::success,
                   "trial continuation cannot mask an unbounded domain escape");
  ys = before;
  for (double &v : buffers[3])
    v = std::numeric_limits<double>::quiet_NaN();
  passed &= expect(sampler.sample(p, .05, ParcelPass::corrector, revision,
                                  ys.data(), ys.size())
                               .status != portable::Status::success &&
                       ys == before,
                   "invalid native composition does not partially overwrite "
                   "the sample buffer");
  return passed;
}

bool test_periodic_stencil() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  if (!compile_geometry(tensor_mesh(), geometry, patch))
    return false;
  const double first = geometry.x().centres().data[0];
  const double last_image = geometry.x().centres().data[3] - 4;
  const double y0 = geometry.y().centres().data[0];
  const double z0 = geometry.z().centres().data[0];
  const auto stencil = build_parcel_grid_coupling_stencil(
      geometry, patch, {0, y0, z0}, {true, false, false});
  bool passed = expect(stencil.succeeded() && stencil.entry_count == 2 &&
                           !stencil.boundary_clamped,
                       "periodic seam interpolates across both end cells");
  double mass = 0, sampled = 0;
  for (std::size_t i = 0; i < stencil.entry_count; ++i) {
    const auto &e = stencil.entries[i];
    const double expected =
        (e.global_index.x == 0 ? -last_image : first) / (first - last_image);
    passed &= expect((e.global_index.x == 0 || e.global_index.x == 3) &&
                         near(e.weight, expected),
                     "stretched periodic weights use end-cell centre spacing");
    // Tensor-linear continuation across this seam, independently manufactured.
    const double x = e.global_index.x == 0 ? first : last_image;
    sampled += e.weight * (7 + 3 * x);
    mass += e.weight * 9.0;
  }
  passed &=
      expect(near(sampled, 7) && near(mass, 9),
             "same periodic stencil interpolates and deposits conservatively");
  const auto wrapped = build_parcel_grid_coupling_stencil(
      geometry, patch, {8, y0, z0}, {true, false, false});
  passed &=
      expect(wrapped.succeeded() && wrapped.entry_count == stencil.entry_count,
             "unwrapped periodic parcel position resolves without clamping");
  for (std::size_t i = 0; i < wrapped.entry_count && i < stencil.entry_count;
       ++i)
    passed &= expect(wrapped.entries[i].global_cell ==
                             stencil.entries[i].global_cell &&
                         wrapped.entries[i].weight == stencil.entries[i].weight,
                     "whole-period translations preserve stencil identity");
  passed &= expect(!build_parcel_grid_coupling_stencil(
                        geometry, patch, {0, -3, z0}, {true, false, false})
                        .available,
                   "periodicity does not mask a physical-domain exit");
  return passed;
}

bool test_shared_stencil() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(compile_geometry(uniform_mesh(), geometry, patch),
                       "uniform stencil geometry compiles");
  const Vector3 position{1.25, 2.0, 2.75};
  const auto stencil =
      build_parcel_grid_coupling_stencil(geometry, patch, position);
  passed &= expect(stencil.succeeded() && stencil.entry_count == 8U,
                   "interior parcel receives one shared trilinear stencil");
  double weight_sum = 0.0;
  double interpolated = 0.0;
  std::array<double, kMaximumStencilEntries> deposited{};
  constexpr double source = 7.5;
  for (std::size_t entry = 0U; entry < stencil.entry_count; ++entry) {
    const auto& node = stencil.entries[entry];
    const double x = geometry.x().centres().data[node.global_index.x];
    const double y = geometry.y().centres().data[node.global_index.y];
    const double z = geometry.z().centres().data[node.global_index.z];
    const double linear = 2.0 * x - 3.0 * y + 0.5 * z + 9.0;
    weight_sum += node.weight;
    interpolated += node.weight * linear;
    deposited[entry] = node.weight * source;
    passed &= expect(std::isfinite(node.weight) && node.weight >= 0.0 &&
                         node.owner_rank == 0,
                     "stencil nodes have finite non-negative global ownership");
  }
  double deposited_sum = 0.0;
  for (const double value : deposited) deposited_sum += value;
  const double exact =
      2.0 * position[0] - 3.0 * position[1] + 0.5 * position[2] + 9.0;
  passed &= expect(near(weight_sum, 1.0) && near(interpolated, exact),
                   "one stencil preserves constants and tensor-linear fields");
  passed &= expect(near(deposited_sum, source),
                   "the same stencil deposits a conservative scalar source");

  MeshPatch decomposed = patch;
  decomposed.begin = {0, 0, 0};
  decomposed.cells = {4, 4, 2};
  // The current decomposition's lexicographic tie break selects {1,1,2}
  // for two ranks on this cubic mesh.
  decomposed.process_grid = {1, 1, 2};
  decomposed.process_coord = {0, 0, 0};
  const auto tie = build_parcel_grid_coupling_stencil(
      geometry, decomposed, {2.0, 2.0, 2.0});
  passed &= expect(tie.succeeded() && tie.entry_count == 8U,
                   "rank-boundary face/edge/corner tie remains a global stencil");
  bool saw_rank_zero = false;
  bool saw_rank_one = false;
  std::uint64_t previous = 0U;
  for (std::size_t entry = 0U; entry < tie.entry_count; ++entry) {
    saw_rank_zero |= tie.entries[entry].owner_rank == 0;
    saw_rank_one |= tie.entries[entry].owner_rank == 1;
    passed &= expect(entry == 0U || tie.entries[entry].global_cell > previous,
                     "tie nodes are ordered only by global cell id");
    previous = tie.entries[entry].global_cell;
  }
  passed &= expect(saw_rank_zero && saw_rank_one,
                   "rank-boundary stencil never falls back to rank-local cells");

  MeshPatch rank_one = decomposed;
  rank_one.begin = {0, 0, 2};
  rank_one.cells = {4, 4, 2};
  rank_one.process_coord = {0, 0, 1};
  const auto rank_one_tie = build_parcel_grid_coupling_stencil(
      geometry, rank_one, {2.0, 2.0, 2.0});
  bool identical = rank_one_tie.succeeded() &&
                   rank_one_tie.entry_count == tie.entry_count;
  for (std::size_t entry = 0U;
       identical && entry < tie.entry_count; ++entry) {
    identical &= rank_one_tie.entries[entry].global_cell ==
                     tie.entries[entry].global_cell &&
                 rank_one_tie.entries[entry].owner_rank ==
                     tie.entries[entry].owner_rank &&
                 rank_one_tie.entries[entry].weight ==
                     tie.entries[entry].weight;
  }
  passed &= expect(identical,
                   "stencil tie identity is independent of observing rank");

  MeshPatch non_authoritative = decomposed;
  non_authoritative.begin = {0, 0, 0};
  non_authoritative.cells = {2, 4, 4};
  non_authoritative.process_grid = {2, 1, 1};
  non_authoritative.process_coord = {0, 0, 0};
  const auto rejected_partition = build_parcel_grid_coupling_stencil(
      geometry, non_authoritative, {2.0, 2.0, 2.0});
  passed &= expect(!rejected_partition.available &&
                       rejected_partition.entry_count == 0U,
                   "non-authoritative Cartesian process grid is rejected");

  CartesianGeometryPlan tensor_geometry;
  MeshPatch tensor_patch;
  passed &= expect(compile_geometry(tensor_mesh(), tensor_geometry,
                                    tensor_patch),
                   "non-uniform tensor geometry compiles");
  const Vector3 tensor_position{
      0.75 * tensor_geometry.x().centres().data[1U] +
          0.25 * tensor_geometry.x().centres().data[2U],
      0.4 * tensor_geometry.y().centres().data[0U] +
          0.6 * tensor_geometry.y().centres().data[1U],
      0.2 * tensor_geometry.z().centres().data[1U] +
          0.8 * tensor_geometry.z().centres().data[2U]};
  const auto tensor = build_parcel_grid_coupling_stencil(
      tensor_geometry, tensor_patch, tensor_position);
  double tensor_linear = 0.0;
  for (std::size_t entry = 0U; entry < tensor.entry_count; ++entry) {
    const auto index = tensor.entries[entry].global_index;
    const double x = tensor_geometry.x().centres().data[index.x];
    const double y = tensor_geometry.y().centres().data[index.y];
    const double z = tensor_geometry.z().centres().data[index.z];
    tensor_linear += tensor.entries[entry].weight *
                     (-1.5 * x + 0.25 * y + 4.0 * z - 2.0);
  }
  const double tensor_exact = -1.5 * tensor_position[0U] +
                              0.25 * tensor_position[1U] +
                              4.0 * tensor_position[2U] - 2.0;
  passed &= expect(tensor.succeeded() && tensor.entry_count == 8U &&
                       near(tensor_linear, tensor_exact),
                   "non-uniform tensor-axis stencil is linear exact");

  const auto boundary = build_parcel_grid_coupling_stencil(
      geometry, patch, {0.0, 2.0, 2.0});
  passed &= expect(boundary.succeeded() && boundary.boundary_clamped,
                   "cell-centred boundary closure is reported explicitly");

  const double nan = std::numeric_limits<double>::quiet_NaN();
  passed &= expect(!build_parcel_grid_coupling_stencil(
                        geometry, patch, {nan, 1.0, 1.0})
                         .available,
                   "non-finite stencil input produces no publishable result");
  const auto invalid_stencil = build_parcel_grid_coupling_stencil(
      geometry, patch, {nan, 1.0, 1.0});
  passed &= expect(invalid_stencil.entry_count == 0U,
                   "failed stencil contains no partial node list");
  return passed;
}

bool test_trajectory_candidates() {
  ParcelTrajectoryInput input;
  input.parcel = parcel({1.0, 2.0, 3.0}, {2.0, -1.0, 0.5});
  input.kind = TrajectoryCandidateKind::predictor;
  input.duration_s = 2.0;
  input.characteristic_length_m = 0.4;
  input.maximum_particle_cfl = 0.5;
  input.acceleration = zero_acceleration;
  const auto ballistic = make_parcel_trajectory_candidate(input);
  bool passed = expect(ballistic.succeeded() && ballistic.segment_count > 1U,
                       "particle CFL creates bounded ordered substeps");
  passed &= expect(near(ballistic.parcel.position_m, {5.0, 0.0, 4.0}) &&
                       near(ballistic.parcel.velocity_m_per_s,
                            input.parcel.velocity_m_per_s) &&
                       near(ballistic.parcel.age_s, 2.25),
                   "ballistic trajectory is analytical");
  for (std::size_t segment = 0U; segment < ballistic.segment_count; ++segment) {
    const auto& item = ballistic.segments[segment];
    const Vector3 displacement{
        item.end_position_m[0] - item.begin_position_m[0],
        item.end_position_m[1] - item.begin_position_m[1],
        item.end_position_m[2] - item.begin_position_m[2]};
    const double distance = std::sqrt(displacement[0] * displacement[0] +
                                      displacement[1] * displacement[1] +
                                      displacement[2] * displacement[2]);
    passed &= expect(item.ordinal == segment &&
                         item.end_time_s > item.begin_time_s &&
                         distance <= 0.2 + 1.0e-12,
                     "each trajectory segment is ordered and CFL bounded");
  }

  const Vector3 acceleration{1.0, -2.0, 0.5};
  input.kind = TrajectoryCandidateKind::corrector;
  input.duration_s = 1.5;
  input.characteristic_length_m = 10.0;
  input.maximum_particle_cfl = 1.0;
  input.acceleration = constant_acceleration;
  input.immutable_acceleration_context = &acceleration;
  const auto accelerated = make_parcel_trajectory_candidate(input);
  const Vector3 expected_position{
      1.0 + 2.0 * 1.5 + 0.5 * acceleration[0] * 1.5 * 1.5,
      2.0 - 1.0 * 1.5 + 0.5 * acceleration[1] * 1.5 * 1.5,
      3.0 + 0.5 * 1.5 + 0.5 * acceleration[2] * 1.5 * 1.5};
  const Vector3 expected_velocity{2.0 + acceleration[0] * 1.5,
                                  -1.0 + acceleration[1] * 1.5,
                                  0.5 + acceleration[2] * 1.5};
  passed &= expect(accelerated.succeeded() &&
                       accelerated.kind == TrajectoryCandidateKind::corrector &&
                       near(accelerated.parcel.position_m, expected_position) &&
                       near(accelerated.parcel.velocity_m_per_s,
                            expected_velocity),
                   "midpoint trajectory is exact for constant acceleration");

  ParcelTrajectoryInput shifted = input;
  const Vector3 frame_velocity{10.0, -4.0, 3.0};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    shifted.parcel.velocity_m_per_s[axis] += frame_velocity[axis];
  }
  const auto galilean = make_parcel_trajectory_candidate(shifted);
  Vector3 shifted_expected = expected_position;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    shifted_expected[axis] += frame_velocity[axis] * input.duration_s;
  }
  passed &= expect(galilean.succeeded() &&
                       near(galilean.parcel.position_m, shifted_expected),
                   "constant-acceleration path is Galilean covariant");

  StokesRelaxationContext stokes{{1.0, 0.25, -0.5}, 2.0};
  input.parcel = parcel({1.0, 2.0, 3.0}, {2.0, -1.0, 0.5});
  input.duration_s = 0.1;
  input.characteristic_length_m = 100.0;
  input.maximum_particle_cfl = 1.0;
  input.event_stop_time_s = -1.0;
  input.acceleration = stokes_relaxation;
  input.immutable_acceleration_context = &stokes;
  const auto stokes_base = make_parcel_trajectory_candidate(input);
  StokesRelaxationContext shifted_stokes = stokes;
  ParcelTrajectoryInput stokes_shifted = input;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    shifted_stokes.carrier_velocity_m_per_s[axis] += frame_velocity[axis];
    stokes_shifted.parcel.velocity_m_per_s[axis] += frame_velocity[axis];
  }
  stokes_shifted.immutable_acceleration_context = &shifted_stokes;
  const auto stokes_frame = make_parcel_trajectory_candidate(stokes_shifted);
  bool stokes_covariant = stokes_base.succeeded() && stokes_frame.succeeded();
  for (std::size_t axis = 0U; stokes_covariant && axis < 3U; ++axis) {
    stokes_covariant &= near(
        stokes_frame.parcel.position_m[axis] -
            stokes_base.parcel.position_m[axis],
        frame_velocity[axis] * input.duration_s);
    stokes_covariant &= near(
        stokes_frame.parcel.velocity_m_per_s[axis] -
            stokes_base.parcel.velocity_m_per_s[axis],
        frame_velocity[axis]);
  }
  passed &= expect(stokes_covariant,
                   "Stokes-relaxation callback remains Galilean covariant");

  input.duration_s = 3.0;
  input.event_stop_time_s = 0.75;
  const auto stopped = make_parcel_trajectory_candidate(input);
  passed &= expect(stopped.succeeded() && near(stopped.advanced_duration_s, 0.75) &&
                       near(stopped.parcel.age_s, 1.0),
                   "known event time stops the candidate exactly");

  input.duration_s = 0.0;
  input.event_stop_time_s = -1.0;
  input.characteristic_length_m = 0.0;
  input.maximum_particle_cfl = 0.0;
  input.maximum_substeps = 0U;
  input.acceleration = nullptr;
  const auto zero_dt = make_parcel_trajectory_candidate(input);
  passed &= expect(zero_dt.succeeded() && zero_dt.segment_count == 0U &&
                       near(zero_dt.parcel.position_m,
                            input.parcel.position_m),
                   "zero dt is an exact no-op candidate");

  input.duration_s = 1.0;
  input.characteristic_length_m = 1.0;
  input.maximum_particle_cfl = 1.0;
  input.maximum_substeps =
      static_cast<std::uint32_t>(kMaximumTrajectorySegments);
  input.acceleration = unavailable_acceleration;
  const auto failed = make_parcel_trajectory_candidate(input);
  passed &= expect(!failed.available &&
                       failed.status ==
                           ParcelTrajectoryStatus::acceleration_unavailable &&
                       failed.segment_count == 0U &&
                       failed.parcel.id == ParcelId{},
                   "acceleration failure cannot publish a partial path");
  return passed;
}

void append_box(std::vector<TriangleInput>& triangles, Real3 lower,
                Real3 upper) {
  const Real3 p000{lower.x, lower.y, lower.z};
  const Real3 p001{lower.x, lower.y, upper.z};
  const Real3 p010{lower.x, upper.y, lower.z};
  const Real3 p011{lower.x, upper.y, upper.z};
  const Real3 p100{upper.x, lower.y, lower.z};
  const Real3 p101{upper.x, lower.y, upper.z};
  const Real3 p110{upper.x, upper.y, lower.z};
  const Real3 p111{upper.x, upper.y, upper.z};
  const std::array<TriangleInput, 12U> box{
      TriangleInput{p000, p001, p011}, TriangleInput{p000, p011, p010},
      TriangleInput{p100, p110, p111}, TriangleInput{p100, p111, p101},
      TriangleInput{p000, p100, p101}, TriangleInput{p000, p101, p001},
      TriangleInput{p010, p011, p111}, TriangleInput{p010, p111, p110},
      TriangleInput{p000, p010, p110}, TriangleInput{p000, p110, p100},
      TriangleInput{p001, p101, p111}, TriangleInput{p001, p111, p011}};
  triangles.insert(triangles.end(), box.begin(), box.end());
}

bool compile_cube_surface(const CartesianGeometryPlan& geometry,
                          const MeshPatch& patch,
                          ImmersedSurfacePlan& surface) {
  std::vector<TriangleInput> triangles;
  append_box(triangles, {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0});
  StlScanPlan scan;
  const StlScanBudget budget{UINT64_C(1073741824), UINT64_C(1073741824),
                             UINT64_C(1000000), UINT64_C(1000000), 1U};
  return static_cast<bool>(StlScanCompiler::compile_triangles(
             geometry, patch, {triangles.data(), triangles.size()},
             CartesianAxis::y, budget, scan)) &&
         static_cast<bool>(ImmersedSurfaceCompiler::compile(scan, surface));
}

bool test_native_ibm_gas_sampler() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  ImmersedSurfacePlan surface;
  if (!compile_geometry(uniform_mesh(), geometry, patch) ||
      !compile_cube_surface(geometry, patch, surface))
    return false;
  const auto n = patch.cells;
  const std::size_t count = std::size_t(n.x) * n.y * n.z;
  std::array<std::vector<double>, 5> data;
  std::array<ConstFieldView, 5> views{};
  for (unsigned f = 0; f < 5; ++f) {
    const unsigned components = f == 2 ? 3 : 1;
    data[f].assign(count * components, f == 1   ? 101000.
                                       : f == 3 ? .01
                                       : f == 4 ? 1.
                                                : 0.);
    auto &v = views[f];
    v.base = data[f].data();
    v.interior = n;
    v.components = components;
    v.stride_y = n.x;
    v.stride_z = std::size_t(n.x) * n.y;
    v.component_stride = count;
    v.field = f + 1;
    v.revision = 1;
    v.storage_identity = f + 101;
    v.revision_domain = 10;
  }
  // The cube contains the cell centred at (1.5, 1.5, 1.5). No gas state is
  // admissible there; masked sampling must never even read its NaNs.
  const auto solid = std::size_t(1 + n.x * (1 + n.y));
  data[4][solid] = 0;
  for (unsigned f = 0; f < 4; ++f)
    for (unsigned c = 0; c < views[f].components; ++c)
      data[f][solid + count * c] = std::numeric_limits<double>::quiet_NaN();
  hundun::v04::detail::ProductParcelGas sampler;
  const std::size_t species = 0;
  const portable::Revision revision{0, 1, 1};
  bool ok = bool(sampler.configure(geometry, patch, {}, 100, {&species, 1}, 1,
                                   true)) &&
            bool(sampler.bind(revision, .1, 100000, views[0], views[1],
                              views[2], {&views[3], 1}));
  ok &= bool(sampler.bind_immersed(surface, views[4], geometry.fingerprint()));
  auto p = parcel({.9, 1.5, 1.5}, {1, 0, 0});
  double y[2]{};
  const auto value =
      sampler.sample(p, .05, ParcelPass::predictor, revision, y, 2);
  auto weights = sampler.stencil(p.position_m);
  ok &= expect(
      value.status == portable::Status::success && near(y[0], .01) &&
          weights.succeeded() && weights.entry_count == 1 &&
          weights.entries[0].global_cell != solid &&
          weights.entries[0].weight == 1,
      "IBM gas sampling and deposition share a positive fluid-only stencil");
  p.position_m[0] = .5;
  ok &= expect(sampler.stencil(p.position_m).succeeded(),
               "fluid cell centre needs no zero-length ray query");
  p.position_m[0] = 1.1;
  spray::detail::ParcelLocation location;
  ok &= expect(
      sampler.sample(p, .05, ParcelPass::corrector, revision, y, 2).status ==
              portable::Status::success &&
          !sampler.stencil(p.position_m).succeeded() &&
          !sampler.locate(p.position_m, revision, location),
      "IBM trial gas extends from the visible fluid surface without admitting "
      "solid positions");
  return ok;
}

ParcelTrajectoryCandidate straight_candidate(Vector3 position,
                                              Vector3 velocity,
                                              double duration) {
  ParcelTrajectoryInput input;
  input.parcel = parcel(position, velocity);
  input.duration_s = duration;
  input.characteristic_length_m = 100.0;
  input.maximum_particle_cfl = 1.0;
  input.acceleration = zero_acceleration;
  return make_parcel_trajectory_candidate(input);
}

struct BallisticInterval final : ParcelIntervalProvider {
  ParcelIntervalReport
  advance(const SprayParcelState &begin, double, double dt, ParcelPass,
          portable::Revision revision) const noexcept override {
    ParcelIntervalReport out;
    out.available = true;
    out.revision = revision;
    out.parcel = begin;
    for (unsigned d = 0; d < 3; ++d)
      out.parcel.position_m[d] += dt * begin.velocity_m_per_s[d];
    out.parcel.age_s += dt;
    out.elapsed_duration_s = dt;
    out.initial_liquid_absolute_enthalpy_j_per_kg = 100;
    out.liquid_absolute_enthalpy_j_per_kg = 100;
    out.exchange.available = true;
    return out;
  }
};
bool test_native_event_geometry() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  ImmersedSurfacePlan surface;
  if (!compile_geometry(uniform_mesh(), geometry, patch) ||
      !compile_cube_surface(geometry, patch, surface))
    return false;
  hundun::v04::detail::ProductParcelGeometry events;
  const portable::Revision revision{11, 37, 1};
  bool ok = bool(events.configure(geometry, {true, true, true})) &&
            bool(events.bind_revision(revision));
  BallisticInterval physics;
  ParcelEventsInput input;
  input.accepted_parcel = parcel({1, .5, .5}, {-1, 0, 0});
  input.accepted_parcel.owner_global_cell = 1;
  input.revision = revision;
  input.interval = &physics;
  input.geometry = &events;
  input.duration_s = .25;
  input.initial_substep_s = .25;
  auto result = integrate_parcel_events(input);
  ok &= expect(
      result.available && result.parcel.owner_global_cell == 0 &&
          result.segment_count == 1 && result.segments[0].global_cell == 0 &&
          near(result.parcel.position_m, {.75, .5, .5}),
      "native exact-face negative crossing changes owner at zero time");
  input.accepted_parcel = parcel({3.75, .5, .5}, {1, 0, 0});
  input.accepted_parcel.owner_global_cell = 3;
  input.duration_s = .5;
  input.initial_substep_s = .5;
  result = integrate_parcel_events(input);
  ok &= expect(result.available && result.segment_count == 2 &&
                   result.segments[0].global_cell == 3 &&
                   result.segments[1].global_cell == 0 &&
                   near(result.parcel.position_m, {4.25, .5, .5}),
               "native periodic crossing retains unwrapped trajectory and "
               "correct owner");
  // Two coincident faces need two topological updates without creating a
  // zero-length exchange segment or losing any physical elapsed time.
  input.accepted_parcel = parcel({1, 1, .5}, {-1, -1, 0});
  input.accepted_parcel.owner_global_cell = 5;
  result = integrate_parcel_events(input);
  ok &= expect(
      result.available && result.parcel.owner_global_cell == 0 &&
          result.segment_count == 1 &&
          near(result.parcel.position_m, {.5, .5, .5}),
      "coincident native cell crossings preserve duration and inventory");
  ok &= bool(events.configure(geometry, {false, false, false}));
  input.accepted_parcel = parcel({3.75, .5, .5}, {1, 0, 0});
  input.accepted_parcel.owner_global_cell = 3;
  result = integrate_parcel_events(input);
  ok &= expect(
      result.available && result.physical_outlet && result.parent_removed &&
          near(result.advanced_duration_s, .25) &&
          near(result.outlet_inventory.mass_kg, 1e-8, 1e-22) &&
          result.exchange.gas_mass_delta_kg == 0,
      "native physical outlet exports inventory once at its actual face");
  input.accepted_parcel = parcel({4, .5, .5}, {1, 0, 0});
  input.accepted_parcel.owner_global_cell = 3;
  result = integrate_parcel_events(input);
  ok &= expect(
      result.available && result.physical_outlet && result.segment_count == 0 &&
          result.advanced_duration_s == 0 &&
          near(result.outlet_inventory.mass_kg, 1e-8, 1e-22),
      "initial outward face event exports without a fictitious interval");
  std::array<bool, 6> walls{};
  walls[1] = true;
  ok &= bool(events.configure(geometry, {false, false, false}, walls));
  input.accepted_parcel = parcel({4, .5, .5}, {1, 0, 0});
  input.accepted_parcel.owner_global_cell = 3;
  result = integrate_parcel_events(input);
  ok &= expect(
      result.available && !result.parent_removed &&
          near(result.parcel.position_m, {3.5, .5, .5}) &&
          near(result.parcel.velocity_m_per_s, {-1, 0, 0}) &&
          near(result.wall_exchange.momentum_kg_m_per_s[0], 2e-8, 1e-22),
      "initial incoming native wall hit rebounds with separate impulse ledger");
  ok &= bool(events.configure(geometry, {false, false, false}, {}, &surface));
  input.accepted_parcel = parcel({.75, 1.5, 1.5}, {1, 0, 0});
  input.accepted_parcel.owner_global_cell = 20;
  result = integrate_parcel_events(input);
  ok &=
      expect(result.available && !result.parent_removed &&
                 near(result.parcel.position_m, {.75, 1.5, 1.5}) &&
                 near(result.parcel.velocity_m_per_s, {-1, 0, 0}) &&
                 near(result.wall_exchange.momentum_kg_m_per_s[0], 2e-8, 1e-22),
             "native IBM surface collision wins over coincident cell crossing");
  input.duration_s = input.initial_substep_s = .25 + 5e-13;
  result = integrate_parcel_events(input);
  ok &= expect(result.available && result.parcel.position_m[0] == 1 &&
                   result.parcel.velocity_m_per_s[0] == -1,
               "native wall contact within event-time tolerance lands on the "
               "actual surface");
  input.duration_s = input.initial_substep_s = .5;
  auto stale = revision;
  ++stale.input_revision;
  auto end = input.accepted_parcel;
  end.position_m[0] = 1.25;
  const auto failed = events.query(input.accepted_parcel, end, 0, .5,
                                   ParcelPass::predictor, stale);
  ok &= expect(!failed.available,
               "native event geometry rejects stale revisions");
  CartesianGeometryPlan stretched;
  MeshPatch stretched_patch;
  ok &= compile_geometry(tensor_mesh(), stretched, stretched_patch);
  ok &= bool(events.configure(stretched, {false, false, false}));
  auto a = parcel({.1, -1, 1.05}, {2, 0, 0});
  auto b = a;
  b.position_m[0] = 2.1;
  const auto crossing =
      events.query(a, b, 0, 1, ParcelPass::corrector, revision);
  const double xface = stretched.x().faces().data[1];
  ok &= expect(crossing.available && crossing.count == 1 &&
                   near(crossing.events[0].elapsed_time_s, (xface - .1) / 2) &&
                   crossing.events[0].next_global_cell == 1,
               "native event locator uses actual nonuniform face coordinates");
  CartesianGeometryPlan one;
  MeshPatch one_patch;
  ok &= compile_geometry(uniform_mesh({1, 1, 1}), one, one_patch);
  ok &= bool(events.configure(one, {true, true, true}));
  input.accepted_parcel = parcel({3.75, .5, .5}, {1, 0, 0});
  result = integrate_parcel_events(input);
  ok &= expect(result.available && result.event_count == 0 &&
                   result.parcel.owner_global_cell == 0 &&
                   near(result.parcel.position_m, {4.25, .5, .5}),
               "one-cell periodic crossing does not invent an owner change");
  return ok;
}
bool test_static_ibm_rebound() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  ImmersedSurfacePlan surface;
  bool passed = expect(compile_geometry(uniform_mesh(), geometry, patch) &&
                           compile_cube_surface(geometry, patch, surface),
                       "closed test surface compiles");

  const auto incoming =
      straight_candidate({3.0, 1.5, 1.5}, {-2.0, 0.5, 0.0}, 1.0);
  StaticIbmReboundConfig config;
  config.fluid_side = ImmersedFluidSide::outside;
  config.normal_restitution = 0.5;
  config.tangential_restitution = 0.25;
  config.maximum_events = 2U;
  const auto bounced = resolve_static_ibm_rebounds(incoming, surface, config);
  passed &= expect(bounced.succeeded() && bounced.event_count == 1U &&
                       bounced.events[0].triangle != kInvalidSurfaceTriangle &&
                       near(bounced.events[0].position_m,
                            {2.0, 1.75, 1.5}) &&
                       near(bounced.events[0].fluid_side_normal,
                            {1.0, 0.0, 0.0}) &&
                       near(bounced.events[0].outgoing_velocity_m_per_s,
                            {1.0, 0.125, 0.0}) &&
                       near(bounced.parcel.position_m,
                            {2.5, 1.8125, 1.5}),
                   "earliest plane hit uses fluid-side normal/tangent restitution");

  ParcelTrajectoryInput segmented_input;
  segmented_input.parcel = parcel({3.0, 1.5, 1.5}, {-2.0, 0.5, 0.0});
  segmented_input.duration_s = 1.0;
  segmented_input.characteristic_length_m = 0.2;
  segmented_input.maximum_particle_cfl = 1.0;
  segmented_input.acceleration = zero_acceleration;
  const auto segmented_path =
      make_parcel_trajectory_candidate(segmented_input);
  const auto segmented_bounce =
      resolve_static_ibm_rebounds(segmented_path, surface, config);
  passed &= expect(segmented_path.segment_count > 1U &&
                       segmented_bounce.succeeded() &&
                       segmented_bounce.event_count == 1U &&
                       near(segmented_bounce.parcel.position_m,
                            bounced.parcel.position_m) &&
                       near(segmented_bounce.parcel.velocity_m_per_s,
                            bounced.parcel.velocity_m_per_s),
                   "rebound transform continues across ordered CFL segments");

  const auto miss_path =
      straight_candidate({3.0, 3.0, 1.5}, {-2.0, 0.0, 0.0}, 1.0);
  const auto missed = resolve_static_ibm_rebounds(miss_path, surface, config);
  passed &= expect(missed.succeeded() && missed.event_count == 0U &&
                       near(missed.parcel.position_m, {1.0, 3.0, 1.5}),
                   "surface miss preserves the trajectory candidate");

  config.maximum_events = 0U;
  const auto limited = resolve_static_ibm_rebounds(incoming, surface, config);
  passed &= expect(!limited.available &&
                       limited.status == StaticIbmReboundStatus::event_limit,
                   "maximum-event bound rejects an incomplete rebound path");

  const auto repeated_hit =
      straight_candidate({1.5, 1.5, 1.5}, {1.0, 0.0, 0.0}, 2.0);
  config.fluid_side = ImmersedFluidSide::inside;
  config.normal_restitution = 1.0;
  config.tangential_restitution = 1.0;
  config.maximum_events = 1U;
  const auto partially_limited =
      resolve_static_ibm_rebounds(repeated_hit, surface, config);
  passed &= expect(!partially_limited.available &&
                       partially_limited.status ==
                           StaticIbmReboundStatus::event_limit &&
                       partially_limited.event_count == 0U &&
                       partially_limited.segment_count == 0U &&
                       partially_limited.parcel.id == ParcelId{},
                   "late event-limit failure discards every partial result");

  const auto wall_departure =
      straight_candidate({2.0, 1.5, 1.5}, {-1.0, 0.0, 0.0}, 1.5);
  config.maximum_events = 2U;
  const auto departure =
      resolve_static_ibm_rebounds(wall_departure, surface, config);
  passed &= expect(departure.succeeded() && departure.event_count == 1U &&
                       near(departure.events[0].event_time_s, 1.0) &&
                       near(departure.events[0].position_m,
                            {1.0, 1.5, 1.5}) &&
                       near(departure.parcel.position_m,
                            {1.5, 1.5, 1.5}),
                   "fraction-zero departure is skipped before the next wall");

  ImmersedSurfacePlan unavailable_surface;
  const auto query_failure =
      resolve_static_ibm_rebounds(incoming, unavailable_surface, {});
  passed &= expect(!query_failure.available && query_failure.event_count == 0U,
                   "unavailable IBM query never publishes a partial path");
  return passed;
}

bool test_outlet_ledger() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(compile_geometry(uniform_mesh(), geometry, patch),
                       "outlet geometry compiles");
  const auto outgoing =
      straight_candidate({3.5, 2.0, 2.0}, {2.0, 0.0, 0.0}, 1.0);
  const auto outlet = resolve_physical_domain_outlet(outgoing, geometry);
  passed &= expect(outlet.succeeded() && !outlet.parcel_retained &&
                       outlet.ledger.recorded &&
                       outlet.ledger.requires_parcel_removal &&
                       !outlet.ledger.permits_gas_source_deposition &&
                       outlet.ledger.face == PhysicalDomainFace::x_max &&
                       near(outlet.ledger.event_time_s, 0.25) &&
                       near(outlet.ledger.position_m, {4.0, 2.0, 2.0}) &&
                       near(outlet.parcel.age_s, 0.5) &&
                       near(outlet.ledger.represented_liquid_mass_kg, 1.0e-8) &&
                       near(outlet.ledger.represented_momentum_kg_m_per_s,
                            {2.0e-8, 0.0, 0.0}),
                   "physical outlet creates a boundary ledger without deposition");

  const auto two_face_path =
      straight_candidate({3.0, 3.0, 2.0}, {2.0, 4.0, 0.0}, 1.0);
  const auto earliest =
      resolve_physical_domain_outlet(two_face_path, geometry);
  passed &= expect(earliest.succeeded() && earliest.ledger.recorded &&
                       earliest.ledger.face == PhysicalDomainFace::y_max &&
                       near(earliest.ledger.event_time_s, 0.25) &&
                       near(earliest.ledger.position_m, {3.5, 4.0, 2.0}),
                   "outlet resolver stops at the earliest crossed face");

  const auto retained_path =
      straight_candidate({1.0, 1.0, 1.0}, {0.25, 0.0, 0.0}, 1.0);
  const auto retained =
      resolve_physical_domain_outlet(retained_path, geometry);
  passed &= expect(retained.succeeded() && retained.parcel_retained &&
                       !retained.ledger.recorded &&
                       near(retained.parcel.position_m, {1.25, 1.0, 1.0}),
                   "in-domain path produces no outlet record");

  ParcelTrajectoryCandidate failed;
  const auto unavailable = resolve_physical_domain_outlet(failed, geometry);
  passed &= expect(!unavailable.available && !unavailable.ledger.recorded,
                   "failed path cannot silently remove or deposit a parcel");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    std::cerr << "FAIL: MPI_Init\n";
    return 1;
  }
  bool passed = test_shared_stencil();
  passed &= test_periodic_stencil();
  passed &= test_native_gas_sampler();
  passed &= test_native_ibm_gas_sampler();
  passed &= test_trajectory_candidates();
  passed &= test_static_ibm_rebound();
  passed &= test_native_event_geometry();
  passed &= test_outlet_ledger();
  if (MPI_Finalize() != MPI_SUCCESS) {
    std::cerr << "FAIL: MPI_Finalize\n";
    return 1;
  }
  if (!passed) return 1;
  std::cout << "PASS: spray mechanics detail\n";
  return 0;
}
