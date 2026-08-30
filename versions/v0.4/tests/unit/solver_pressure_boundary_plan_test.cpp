// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_flow.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {4.0, 3.0, 2.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {4, 3, 2};
  mesh.minimum_spacing = {1.0, 1.0, 1.0};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {24U, 1U << 20U};
  return mesh;
}

BoundaryFaceSpec wall() {
  BoundaryFaceSpec face;
  face.flow_kind = BoundaryKind::no_slip_wall;
  face.thermal_kind = BoundaryKind::adiabatic_wall;
  face.mach_limit = 0.95;
  return face;
}

ValidatedModel model_spec(BoundaryKind outlet_kind) {
  ValidatedModel model;
  model.fingerprint = outlet_kind == BoundaryKind::nscbc_outlet
                          ? UINT64_C(0x5043424e53434243)
                          : UINT64_C(0x5043425052455353);
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : model.boundaries) face = wall();

  BoundaryFaceSpec& outlet = model.boundaries[1U];
  outlet.flow_kind = outlet_kind;
  outlet.thermal_kind = BoundaryKind::none;
  outlet.pressure = 101325.0;
  outlet.temperature = 300.0;
  outlet.backflow_velocity = {-0.1, 0.0, 0.0};
  outlet.backflow_temperature = 301.0;
  outlet.relaxation = 0.1;
  outlet.mach_limit = 0.9;
  outlet.allow_backflow = true;

  model.boundaries[4U].flow_kind = BoundaryKind::periodic;
  model.boundaries[4U].thermal_kind = BoundaryKind::none;
  model.boundaries[5U].flow_kind = BoundaryKind::periodic;
  model.boundaries[5U].thermal_kind = BoundaryKind::none;
  model.schemes = SchemeSpec{};
  model.time = TimeControlSpec{};
  return model;
}

bool compile_boundary(MPI_Comm communicator, const ValidatedModel& model,
                      const CartesianGeometryPlan& geometry,
                      const MeshPatch& patch, BoundaryPlan& boundary) {
  FieldRegistry registry;
  SchemePlan schemes;
  TimeSchemePlan time;
  return static_cast<bool>(BoundaryCompiler::compile(
      communicator, model, geometry, patch, registry, boundary, schemes,
      time));
}

struct OwnedField {
  std::vector<double> values;
  FieldView view{};
};

OwnedField make_field(Int3 cells) {
  constexpr std::int32_t ghosts = 1;
  OwnedField field;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  field.values.assign(nx * ny * nz, -777.0);
  field.view = {field.values.data() + 1U + nx + nx * ny,
                cells,
                {ghosts, ghosts, ghosts},
                1U,
                nx,
                nx * ny,
                nx * ny * nz,
                0U,
                101U,
                201U,
                301U,
                401U};
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        field.view.unchecked({x, y, z}, 0U) =
            10.0 + x + 2.0 * y + 3.0 * z;
      }
    }
  }
  return field;
}

bool close(double actual, double expected) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <=
             1.0e-13 * std::max(1.0, std::abs(expected));
}

bool run_test() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_SELF, mesh_spec(), GeometryBudget{0U, 1U}, geometry,
          patch)),
      "geometry compiles");
  BoundaryPlan boundary;
  passed &= expect(compile_boundary(MPI_COMM_SELF,
                                    model_spec(BoundaryKind::nscbc_outlet),
                                    geometry, patch, boundary),
                   "NSCBC boundary compiles");

  PressureCorrectionBoundaryPlan pressure_boundary;
  passed &= expect(static_cast<bool>(PressureCorrectionBoundaryPlan::compile(
                       geometry, patch, boundary, pressure_boundary)),
                   "pressure-correction boundary authority compiles");
  const PressureCorrectionBoundaryCertificate certificate =
      pressure_boundary.certificate();
  passed &= expect(certificate.valid() &&
                       certificate.geometry == geometry.topology_revision() &&
                       certificate.source_revision == boundary.revision(),
                   "authority certifies semantic and source revisions");
  passed &= expect(
      pressure_boundary.kind(CartesianFace::x_max) ==
          PressureCorrectionFaceKind::homogeneous_dirichlet,
      "NSCBC outlet is homogeneous pressure-correction Dirichlet");
  passed &= expect(
      pressure_boundary.kind(CartesianFace::x_min) ==
          PressureCorrectionFaceKind::homogeneous_neumann,
      "wall remains homogeneous pressure-correction Neumann");
  passed &= expect(pressure_boundary.kind(CartesianFace::z_min) ==
                           PressureCorrectionFaceKind::periodic &&
                       pressure_boundary.kind(CartesianFace::z_max) ==
                           PressureCorrectionFaceKind::periodic,
                   "periodic pair remains exchanged periodic correction");

  MgBoundarySet mg{};
  passed &= expect(static_cast<bool>(pressure_boundary.mg_boundaries(mg)) &&
                       mg.x_min == MgBoundaryKind::neumann &&
                       mg.x_max == MgBoundaryKind::dirichlet &&
                       mg.y_min == MgBoundaryKind::neumann &&
                       mg.y_max == MgBoundaryKind::neumann &&
                       mg.z_min == MgBoundaryKind::periodic &&
                       mg.z_max == MgBoundaryKind::periodic,
                   "MG consumes the same six-face authority");

  OwnedField correction = make_field(patch.cells);
  correction.view.unchecked({1, 1, -1}, 0U) = -4.0;
  passed &= expect(
      static_cast<bool>(pressure_boundary.fill_ghosts(correction.view)),
      "authority fills physical correction ghosts");
  const double x_min_owner = correction.view.unchecked({0, 1, 1}, 0U);
  const double x_max_owner =
      correction.view.unchecked({patch.cells.x - 1, 1, 1}, 0U);
  passed &= expect(
      close(correction.view.unchecked({-1, 1, 1}, 0U), x_min_owner) &&
          close(correction.view.unchecked({patch.cells.x, 1, 1}, 0U),
                -x_max_owner) &&
          close(correction.view.unchecked({1, 1, -1}, 0U), -4.0),
      "Neumann copies, Dirichlet reflects, periodic ghost is untouched");
  passed &= expect(
      close(pressure_boundary.neighbor_value(
                as_const(correction.view), {0, 1, 1}, CartesianAxis::x, -1),
            x_min_owner) &&
          close(pressure_boundary.neighbor_value(
                    as_const(correction.view), {patch.cells.x - 1, 1, 1},
                    CartesianAxis::x, 1),
                0.0) &&
          close(pressure_boundary.neighbor_value(
                    as_const(correction.view), {1, 1, 0}, CartesianAxis::z,
                    -1),
                -4.0),
      "operator neighbours use Neumann owner, Dirichlet zero, and periodic "
      "halo values");

  const double dirichlet_jump = pressure_boundary.jump(
      as_const(correction.view), CartesianAxis::x,
      {patch.cells.x, 1, 1});
  const double periodic_jump = pressure_boundary.jump(
      as_const(correction.view), CartesianAxis::z, {1, 1, 0});
  passed &= expect(close(dirichlet_jump, -x_max_owner) &&
                       close(pressure_boundary.jump(as_const(correction.view),
                                                    CartesianAxis::x,
                                                    {0, 1, 1}),
                             0.0) &&
                       close(periodic_jump,
                             correction.view.unchecked({1, 1, 0}, 0U) + 4.0),
                   "face jump follows Dirichlet, Neumann, and halo semantics");
  passed &= expect(
      close(pressure_boundary.mass_flux_response(
                as_const(correction.view), CartesianAxis::x,
                {patch.cells.x, 1, 1}, 3.0),
            -3.0 * dirichlet_jump),
      "mass-flux response is derived from the same jump");

  BoundaryPlan replacement;
  passed &= expect(compile_boundary(MPI_COMM_SELF,
                                    model_spec(BoundaryKind::pressure_outlet),
                                    geometry, patch, replacement),
                   "replacement boundary compiles");
  boundary = std::move(replacement);
  const std::vector<double> before = correction.values;
  MgBoundarySet untouched{MgBoundaryKind::periodic, MgBoundaryKind::periodic,
                          MgBoundaryKind::periodic, MgBoundaryKind::periodic,
                          MgBoundaryKind::periodic, MgBoundaryKind::periodic};
  passed &= expect(!pressure_boundary.current() &&
                       !pressure_boundary.mg_boundaries(untouched) &&
                       !pressure_boundary.fill_ghosts(correction.view) &&
                       correction.values == before &&
                       std::isnan(pressure_boundary.jump(
                           as_const(correction.view), CartesianAxis::x,
                           {patch.cells.x, 1, 1})),
                   "stale authority fails closed without partial ghost writes");
  return passed;
}

bool run_decomposition_test(int rank, int size) {
  CartesianMeshSpec mesh = mesh_spec();
  mesh.upper.x = static_cast<double>(8 * size);
  mesh.exact_cells.x = 8 * size;
  mesh.limits.max_global_cells =
      static_cast<std::uint64_t>(mesh.exact_cells.x) *
      static_cast<std::uint64_t>(mesh.exact_cells.y) *
      static_cast<std::uint64_t>(mesh.exact_cells.z);
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool local = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_WORLD, mesh, GeometryBudget{0U, 1U}, geometry, patch)),
      "decomposed geometry compiles");
  BoundaryPlan boundary;
  local &= expect(
      compile_boundary(MPI_COMM_WORLD,
                       model_spec(BoundaryKind::nscbc_outlet), geometry,
                       patch, boundary),
      "decomposed NSCBC boundary compiles");
  PressureCorrectionBoundaryPlan authority;
  local &= expect(static_cast<bool>(PressureCorrectionBoundaryPlan::compile(
                      geometry, patch, boundary, authority)),
                  "decomposed pressure boundary authority compiles");
  const int setup_value = local ? 1 : 0;
  int global_setup = 0;
  MPI_Allreduce(&setup_value, &global_setup, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
  if (global_setup == 0) return false;

  const PressureCorrectionBoundaryCertificate certificate =
      authority.certificate();
  std::uint64_t semantic_min = certificate.semantic;
  std::uint64_t semantic_max = certificate.semantic;
  MPI_Allreduce(MPI_IN_PLACE, &semantic_min, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &semantic_max, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  local &= expect(semantic_min == semantic_max && semantic_min != 0U,
                  "semantic authority is decomposition independent");

  std::vector<std::uint64_t> layouts(static_cast<std::size_t>(size));
  MPI_Allgather(&certificate.rank_local_layout, 1, MPI_UINT64_T,
                layouts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
  std::sort(layouts.begin(), layouts.end());
  local &= expect(size == 1 ||
                      std::adjacent_find(layouts.begin(), layouts.end(),
                                         std::not_equal_to<>{}) !=
                          layouts.end(),
                  "rank-local layout authority differs across decomposition");

  PressureCorrectionFaceRule x_low;
  PressureCorrectionFaceRule x_high;
  local &= expect(authority.face_rule(CartesianAxis::x, {0, 0, 0}, x_low) &&
                      authority.face_rule(CartesianAxis::x,
                                          {patch.cells.x, 0, 0}, x_high) &&
                      x_low.physical == (patch.begin.x == 0) &&
                      x_high.physical ==
                          (patch.begin.x + patch.cells.x ==
                           geometry.global_cells().x) &&
                      (!x_low.physical ||
                       x_low.kind ==
                           PressureCorrectionFaceKind::homogeneous_neumann) &&
                      (!x_high.physical ||
                       x_high.kind ==
                           PressureCorrectionFaceKind::homogeneous_dirichlet),
                  "local edge rule distinguishes MPI, wall, and NSCBC faces");

  if (size > 1) {
    CartesianAxis interface_axis = CartesianAxis::x;
    bool high = false;
    bool found = false;
    for (std::uint8_t axis_index = 0U; axis_index < 3U && !found;
         ++axis_index) {
      const auto axis = static_cast<CartesianAxis>(axis_index);
      const std::int32_t begin = axis == CartesianAxis::x
                                     ? patch.begin.x
                                     : (axis == CartesianAxis::y
                                            ? patch.begin.y
                                            : patch.begin.z);
      const std::int32_t cells = axis == CartesianAxis::x
                                     ? patch.cells.x
                                     : (axis == CartesianAxis::y
                                            ? patch.cells.y
                                            : patch.cells.z);
      const std::int32_t global = axis == CartesianAxis::x
                                      ? geometry.global_cells().x
                                      : (axis == CartesianAxis::y
                                             ? geometry.global_cells().y
                                             : geometry.global_cells().z);
      if (begin > 0) {
        interface_axis = axis;
        high = false;
        found = true;
      } else if (begin + cells < global) {
        interface_axis = axis;
        high = true;
        found = true;
      }
    }
    local &= expect(found, "multi-rank patch owns an MPI interface");
    if (found) {
      OwnedField correction = make_field(patch.cells);
      Int3 face{0, 0, 0};
      if (interface_axis == CartesianAxis::x) face.x = high ? patch.cells.x : 0;
      if (interface_axis == CartesianAxis::y) face.y = high ? patch.cells.y : 0;
      if (interface_axis == CartesianAxis::z) face.z = high ? patch.cells.z : 0;
      PressureCorrectionFaceRule rule;
      local &= expect(authority.face_rule(interface_axis, face, rule) &&
                          !rule.physical &&
                          rule.kind == PressureCorrectionFaceKind::exchanged,
                      "MPI interface is exchanged rather than physical");
      Int3 left = face;
      if (interface_axis == CartesianAxis::x) --left.x;
      if (interface_axis == CartesianAxis::y) --left.y;
      if (interface_axis == CartesianAxis::z) --left.z;
      const Int3 ghost = high ? face : left;
      const Int3 owner = high ? left : face;
      const double owner_value = correction.view.unchecked(owner, 0U);
      correction.view.unchecked(ghost, 0U) = -100.0 - rank;
      const double ghost_value = correction.view.unchecked(ghost, 0U);
      local &= expect(authority.fill_ghosts(correction.view) &&
                          close(correction.view.unchecked(ghost, 0U),
                                ghost_value) &&
                          close(authority.neighbor_value(
                                    as_const(correction.view), owner,
                                    interface_axis, high ? 1 : -1),
                                ghost_value) &&
                          close(authority.jump(as_const(correction.view),
                                               interface_axis, face),
                                high ? ghost_value - owner_value
                                     : owner_value - ghost_value),
                      "MPI ghost supplies the operator neighbour and face "
                      "jump");
    }

    CartesianGeometryPlan foreign_geometry;
    MeshPatch foreign_patch;
    BoundaryPlan foreign_boundary;
    PressureCorrectionBoundaryPlan rejected;
    CartesianMeshSpec foreign_mesh = mesh;
    foreign_mesh.upper = {static_cast<double>(patch.cells.x),
                          static_cast<double>(patch.cells.y),
                          static_cast<double>(patch.cells.z)};
    foreign_mesh.exact_cells = patch.cells;
    foreign_mesh.limits.max_global_cells =
        static_cast<std::uint64_t>(patch.cells.x) *
        static_cast<std::uint64_t>(patch.cells.y) *
        static_cast<std::uint64_t>(patch.cells.z);
    local &= expect(
        CartesianGeometryCompiler::compile(MPI_COMM_SELF, foreign_mesh,
                                           GeometryBudget{0U, 1U},
                                           foreign_geometry, foreign_patch) &&
            compile_boundary(MPI_COMM_SELF,
                             model_spec(BoundaryKind::nscbc_outlet),
                             foreign_geometry, foreign_patch,
                             foreign_boundary) &&
            foreign_boundary.local_cells().x == patch.cells.x &&
            foreign_boundary.local_cells().y == patch.cells.y &&
            foreign_boundary.local_cells().z == patch.cells.z &&
            !PressureCorrectionBoundaryPlan::compile(
                geometry, patch, foreign_boundary, rejected),
        "same-shape foreign rank-local ownership fails closed");
  }

  const int local_value = local ? 1 : 0;
  int global_value = 0;
  MPI_Allreduce(&local_value, &global_value, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
  return global_value != 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  bool passed = run_test();
  passed = run_decomposition_test(rank, size) && passed;
  MPI_Finalize();
  return passed ? 0 : 1;
}
