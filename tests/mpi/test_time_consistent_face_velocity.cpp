// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/momentum_predictor.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace {

using hundun::boundary::BoundaryRegistry;
using hundun::config::BoundaryType;
using hundun::config::DensityModel;
using hundun::config::FlowBoundaryConfig;
using hundun::config::FlowCaseConfig;
using hundun::config::FlowScalarConfig;
using hundun::config::InletScalarValue;
using hundun::config::InletThermalAuthority;
using hundun::config::PatchName;
using hundun::config::SimulationType;
using hundun::finite_volume::declare_face_mass_flux;
using hundun::flow::make_momentum_time_stencil;
using hundun::flow::MomentumFaceHistory;
using hundun::flow::MomentumTimeOrder;
using hundun::flow::TimeConsistentFaceVelocity;
using hundun::mesh::AnalyticWarpedBoxMapping;
using hundun::mesh::EntityOwnership;
using hundun::mesh::FaceAxis;
using hundun::mesh::FaceSide;
using hundun::mesh::GlobalFaceId;
using hundun::mesh::LocalCellId;
using hundun::mesh::LocalFaceId;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::mesh::UniformBoxMapping;
using hundun::runtime::AccessMode;
using hundun::runtime::ActorId;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::ExchangePlan;
using hundun::runtime::FieldAccessPlan;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldLayoutSet;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::HaloExchange;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::OutputPolicy;
using hundun::runtime::PhaseId;
using hundun::runtime::Real3;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::StructuredDecomposition;
namespace allocation_probe = hundun::test::allocation_probe;

constexpr PhaseId kAssemblePhase = 17U;
constexpr ActorId kAssembleActor = 170U;
constexpr PhaseId kReadPhase = 18U;
constexpr ActorId kReadActor = 180U;
constexpr PhaseId kInitPhase = 19U;
constexpr ActorId kInitActor = 190U;

template <class Function> void expect_error(Function &&function) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const hundun::runtime::Error &) {
    threw = true;
  }
  HUNDUN_CHECK(threw);
}

Int3 process_grid_for(int ranks, Int3 extent) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return extent.x >= 2 ? Int3{2, 1, 1} : Int3{1, 2, 1};
  if (ranks == 4)
    return extent.x >= 2 ? Int3{2, 2, 1} : Int3{1, 2, 2};
  throw hundun::runtime::Error("unsupported Task 17 test rank count");
}

FlowCaseConfig periodic_case(double rho_ref) {
  FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type = SimulationType::variable_density_flow;
  config.density_model = DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = rho_ref;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back(FlowScalarConfig{"alpha", 0.125});
  constexpr std::array<PatchName, 6> names{PatchName::x_min, PatchName::x_max,
                                           PatchName::y_min, PatchName::y_max,
                                           PatchName::z_min, PatchName::z_max};
  for (std::size_t index = 0; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type = BoundaryType::periodic;
  }
  return config;
}

FlowCaseConfig mixed_case(double rho_ref) {
  auto config = periodic_case(rho_ref);
  constexpr std::array<BoundaryType, 6> kinds{
      BoundaryType::velocity_inlet, BoundaryType::pressure_outlet,
      BoundaryType::no_slip_wall,   BoundaryType::symmetry,
      BoundaryType::periodic,       BoundaryType::periodic};
  for (std::size_t index = 0; index < kinds.size(); ++index) {
    config.boundaries[index].type = kinds[index];
  }
  config.boundaries[0].velocity_m_per_s = Real3{1.0, 0.2, -0.1};
  config.boundaries[0].thermal_authority = InletThermalAuthority::enthalpy;
  config.boundaries[0].enthalpy_J_per_kg = 12.0;
  config.boundaries[0].scalar_values =
      std::vector<InletScalarValue>{InletScalarValue{"alpha", 0.3}};
  config.boundaries[1].pressure_perturbation_pa = 17.0;
  return config;
}

FieldDescriptor cell_descriptor(std::string name, std::uint32_t components) {
  return {std::move(name),
          "1",
          "task17",
          FunctionSpace::cell_average,
          ScalarType::float64,
          components,
          1,
          false,
          RestartPolicy::transient,
          OutputPolicy::never};
}

FieldDescriptor face_descriptor(std::string name, std::uint32_t components) {
  return {std::move(name),
          "1",
          "task17",
          FunctionSpace::face_value,
          ScalarType::float64,
          components,
          0,
          false,
          RestartPolicy::transient,
          OutputPolicy::never};
}

struct Fields final {
  FieldId predictor{};
  FieldId pressure{};
  FieldId pressure_gradient{};
  FieldId diagonal{};
  FieldId velocity_n{};
  FieldId velocity_nm1{};
  FieldId face_n{};
  FieldId face_nm1{};
  FieldId trial_face{};
  FieldId mass_flux{};
};

Fields declare_fields(FieldRegistry &registry) {
  Fields ids{};
  ids.predictor = registry.declare_field(cell_descriptor("u_predictor", 3U));
  ids.pressure = registry.declare_field(cell_descriptor("pi", 1U));
  ids.pressure_gradient =
      registry.declare_field(cell_descriptor("grad_pi", 3U));
  ids.diagonal = registry.declare_field(cell_descriptor("momentum_aP", 3U));
  ids.velocity_n = registry.declare_field(cell_descriptor("u_n", 3U));
  ids.velocity_nm1 = registry.declare_field(cell_descriptor("u_nm1", 3U));
  ids.face_n = registry.declare_field(face_descriptor("u_face_n", 3U));
  ids.face_nm1 = registry.declare_field(face_descriptor("u_face_nm1", 3U));
  ids.trial_face = registry.declare_field(face_descriptor("u_face_trial", 3U));
  ids.mass_flux = declare_face_mass_flux(registry);
  registry.freeze();
  return ids;
}

FieldAccessPlan make_access_plan(const FieldRegistry &registry,
                                 const Fields &ids) {
  FieldAccessPlan plan(registry);
  plan.declare_access(kInitPhase, kInitActor, ids.face_n,
                      AccessMode::read_write);
  plan.declare_access(kInitPhase, kInitActor, ids.face_nm1,
                      AccessMode::read_write);
  plan.declare_access(kAssemblePhase, kAssembleActor, ids.trial_face,
                      AccessMode::write);
  plan.declare_access(kAssemblePhase, kAssembleActor, ids.mass_flux,
                      AccessMode::write);
  plan.declare_access(kReadPhase, kReadActor, ids.mass_flux, AccessMode::read);
  plan.freeze();
  return plan;
}

std::uint64_t bits(double value) {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

Real3 subtract(Real3 left, Real3 right) {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

Real3 multiply(double scale, Real3 value) {
  return {scale * value.x, scale * value.y, scale * value.z};
}

Real3 add(Real3 left, Real3 right) {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

double dot(Real3 left, Real3 right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

Real3 cross(Real3 left, Real3 right) {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

double norm(Real3 value) { return std::sqrt(dot(value, value)); }

std::array<Int3, 4> face_vertices(hundun::mesh::LogicalFace face) {
  const Int3 c = face.coordinate;
  if (face.axis == FaceAxis::x) {
    return {Int3{c.x, c.y, c.z}, Int3{c.x, c.y + 1, c.z},
            Int3{c.x, c.y + 1, c.z + 1}, Int3{c.x, c.y, c.z + 1}};
  }
  if (face.axis == FaceAxis::y) {
    return {Int3{c.x, c.y, c.z}, Int3{c.x, c.y, c.z + 1},
            Int3{c.x + 1, c.y, c.z + 1}, Int3{c.x + 1, c.y, c.z}};
  }
  return {Int3{c.x, c.y, c.z}, Int3{c.x + 1, c.y, c.z},
          Int3{c.x + 1, c.y + 1, c.z}, Int3{c.x, c.y + 1, c.z}};
}

std::array<Real3, 4>
canonical_periodic_vertices(const MeshGeometry &geometry,
                            hundun::mesh::LogicalFace logical) {
  if (logical.axis == FaceAxis::x)
    logical.coordinate.x = 0;
  if (logical.axis == FaceAxis::y)
    logical.coordinate.y = 0;
  if (logical.axis == FaceAxis::z)
    logical.coordinate.z = 0;
  const auto coordinates = face_vertices(logical);
  std::array<Real3, 4> result{};
  for (std::size_t index = 0; index < result.size(); ++index) {
    result[index] = geometry.vertex_position_m(coordinates[index]);
  }
  return result;
}

Real3 canonical_center(const std::array<Real3, 4> &vertices, bool uniform) {
  const Real3 anchor = vertices[0];
  if (uniform)
    return add(anchor, multiply(0.5, subtract(vertices[2], anchor)));
  Real3 center = anchor;
  for (std::size_t index = 1; index < vertices.size(); ++index) {
    center = add(center, multiply(0.25, subtract(vertices[index], anchor)));
  }
  return center;
}

Real3 positive_area(const std::array<Real3, 4> &vertices) {
  return add(multiply(0.5, cross(subtract(vertices[1], vertices[0]),
                                 subtract(vertices[2], vertices[0]))),
             multiply(0.5, cross(subtract(vertices[2], vertices[0]),
                                 subtract(vertices[3], vertices[0]))));
}

struct StructuredIndex final {
  int i{};
  int j{};
  int k{};
};

StructuredIndex map_cell(Int3 global, hundun::runtime::Box3 box, Int3 extent) {
  const Int3 local{box.end.x - box.begin.x, box.end.y - box.begin.y,
                   box.end.z - box.begin.z};
  const auto axis = [](int coordinate, int begin, int end, int global_n,
                       int local_n) {
    if (coordinate >= begin && coordinate < end)
      return coordinate - begin;
    if (coordinate == begin - 1 || (begin == 0 && coordinate == global_n - 1))
      return -1;
    if (coordinate == end || (end == global_n && coordinate == 0))
      return local_n;
    throw hundun::runtime::Error("test cell mapping failed");
  };
  return {axis(global.x, box.begin.x, box.end.x, extent.x, local.x),
          axis(global.y, box.begin.y, box.end.y, extent.y, local.y),
          axis(global.z, box.begin.z, box.end.z, extent.z, local.z)};
}

StructuredIndex shift(StructuredIndex value, FaceAxis axis, int amount) {
  if (axis == FaceAxis::x)
    value.i += amount;
  if (axis == FaceAxis::y)
    value.j += amount;
  if (axis == FaceAxis::z)
    value.k += amount;
  return value;
}

template <class T>
double value(const hundun::runtime::FieldView<T> &view, StructuredIndex index,
             int component) {
  return view(index.i, index.j, index.k, component);
}

struct OracleFace final {
  StructuredIndex p{};
  StructuredIndex n{};
  bool reversed{};
  Real3 area{};
  Real3 displacement{};
  Real3 unit_normal{};
  double normal_distance{};
  double weight_p{};
  double weight_n{};
  double volume_p{};
  double volume_n{};
};

OracleFace oracle_face(const MeshTopology &topology,
                       const MeshGeometry &geometry, LocalFaceId face) {
  const auto neighbour = topology.neighbour(face);
  HUNDUN_CHECK(neighbour.has_value());
  const auto box = topology.owned_global_box();
  StructuredIndex owner = map_cell(topology.global_cell(topology.owner(face)),
                                   box, topology.global_extent());
  StructuredIndex other =
      map_cell(topology.global_cell(*neighbour), box, topology.global_extent());
  if (topology.periodic_pair(face).has_value()) {
    const auto logical = topology.logical_face(face);
    const int plane = logical.axis == FaceAxis::x   ? logical.coordinate.x
                      : logical.axis == FaceAxis::y ? logical.coordinate.y
                                                    : logical.coordinate.z;
    other = shift(owner, logical.axis, plane == 0 ? -1 : 1);
  }
  const bool reversed =
      topology.periodic_pair(face).has_value() &&
      topology.global_face_id(face) > *topology.periodic_pair(face);
  const LocalCellId p_cell = reversed ? *neighbour : topology.owner(face);
  const LocalCellId n_cell = reversed ? topology.owner(face) : *neighbour;
  OracleFace result{};
  result.p = reversed ? other : owner;
  result.n = reversed ? owner : other;
  result.reversed = reversed;
  Real3 p_to_face{};
  if (topology.periodic_pair(face).has_value()) {
    const auto vertices =
        canonical_periodic_vertices(geometry, topology.logical_face(face));
    const Real3 center =
        canonical_center(vertices, geometry.mapping_kind() ==
                                       hundun::mesh::MappingKind::uniform_box);
    result.area = multiply(-1.0, positive_area(vertices));
    const Real3 p_center = geometry.cell_center_m(p_cell);
    result.displacement = subtract(geometry.cell_center_m(n_cell), p_center);
    const Real3 length = geometry.length_m();
    const FaceAxis axis = topology.logical_face(face).axis;
    if (axis == FaceAxis::x)
      result.displacement.x -= length.x;
    if (axis == FaceAxis::y)
      result.displacement.y -= length.y;
    if (axis == FaceAxis::z)
      result.displacement.z -= length.z;
    p_to_face = subtract(center, p_center);
  } else {
    result.area = geometry.face_area_vector_m2(face, FaceSide::owner);
    result.displacement = geometry.face_displacement_m(face);
    p_to_face = subtract(geometry.face_center_m(face),
                         geometry.cell_center_m(topology.owner(face)));
  }
  result.unit_normal = multiply(1.0 / norm(result.area), result.area);
  result.normal_distance = dot(result.displacement, result.unit_normal);
  const double a = norm(p_to_face);
  const double b = norm(subtract(result.displacement, p_to_face));
  result.weight_p = b / (a + b);
  result.weight_n = a / (a + b);
  result.volume_p = geometry.cell_volume_m3(p_cell);
  result.volume_n = geometry.cell_volume_m3(n_cell);
  return result;
}

double interpolate(double wp, double p, double wn, double n) {
  return p == n ? p : wp * p + wn * n;
}

void exchange_cells(const StructuredDecomposition &decomposition,
                    FieldStorage &storage, const Fields &ids) {
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  for (FieldId field : {ids.predictor, ids.pressure, ids.pressure_gradient,
                        ids.diagonal, ids.velocity_n, ids.velocity_nm1}) {
    halo.exchange(storage, field);
  }
}

void initialize_periodic_fields(const StructuredDecomposition &decomposition,
                                const MeshTopology &topology,
                                const MeshGeometry &geometry,
                                FieldStorage &storage, const Fields &ids,
                                const FieldAccessPlan &access, double rho_ref,
                                double dt, double alpha0) {
  auto predictor = storage.view<double>(ids.predictor);
  auto pressure = storage.view<double>(ids.pressure);
  auto gradient = storage.view<double>(ids.pressure_gradient);
  auto diagonal = storage.view<double>(ids.diagonal);
  auto velocity_n = storage.view<double>(ids.velocity_n);
  auto velocity_nm1 = storage.view<double>(ids.velocity_nm1);
  const auto box = topology.owned_global_box();
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        const Int3 global{box.begin.x + i, box.begin.y + j, box.begin.z + k};
        const LocalCellId local =
            (static_cast<std::size_t>(k) *
                 static_cast<std::size_t>(decomposition.local_extent().y) +
             static_cast<std::size_t>(j)) *
                static_cast<std::size_t>(decomposition.local_extent().x) +
            static_cast<std::size_t>(i);
        const double volume = geometry.cell_volume_m3(local);
        pressure(i, j, k, 0) =
            ((global.x + global.y + global.z) % 2 == 0) ? 1.0 : -1.0;
        for (int component = 0; component < 3; ++component) {
          predictor(i, j, k, component) =
              0.125 * static_cast<double>(component + 1) +
              0.01 * static_cast<double>(global.x + 2 * global.y + global.z);
          gradient(i, j, k, component) = 0.0;
          velocity_n(i, j, k, component) = 0.0;
          velocity_nm1(i, j, k, component) = 0.0;
          const double factor =
              alpha0 * rho_ref / dt + 0.25 * static_cast<double>(component);
          diagonal(i, j, k, component) = volume * factor;
        }
      }
    }
  }
  auto face_n = storage.acquire_face_write<double>(access, kInitPhase,
                                                   kInitActor, ids.face_n);
  auto face_nm1 = storage.acquire_face_write<double>(access, kInitPhase,
                                                     kInitActor, ids.face_nm1);
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    constexpr std::array<double, 3> current{0.75, -0.5, 0.25};
    constexpr std::array<double, 3> older{0.25, 0.125, -0.375};
    for (int component = 0; component < 3; ++component) {
      face_n(face, component) = current[static_cast<std::size_t>(component)];
      face_nm1(face, component) = older[static_cast<std::size_t>(component)];
    }
  }
  exchange_cells(decomposition, storage, ids);
}

struct FaceRecord final {
  std::uint64_t key{};
  std::uint64_t global{};
  std::uint64_t periodic_pair{std::numeric_limits<std::uint64_t>::max()};
  std::array<std::uint64_t, 3> velocity{};
  std::uint64_t flux{};
};

void verify_canonical_records(
    const MpiContext &mpi, const MeshTopology &topology,
    const hundun::runtime::FaceFieldView<double> &trial,
    const hundun::runtime::FaceFieldView<const double> &flux) {
  std::vector<FaceRecord> local;
  local.reserve(topology.local_face_count());
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    if (!topology.neighbour(face).has_value())
      continue;
    const auto global = topology.global_face_id(face);
    const auto pair = topology.periodic_pair(face);
    const auto key = pair.has_value() ? std::min(global, *pair) : global;
    FaceRecord record{};
    record.key = key;
    record.global = global;
    if (pair.has_value())
      record.periodic_pair = *pair;
    for (int component = 0; component < 3; ++component) {
      record.velocity[static_cast<std::size_t>(component)] =
          bits(trial(face, component));
    }
    record.flux = bits(flux(face, 0));
    local.push_back(record);
  }
  const int local_bytes = static_cast<int>(local.size() * sizeof(FaceRecord));
  std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
  HUNDUN_CHECK(MPI_Allgather(&local_bytes, 1, MPI_INT, counts.data(), 1,
                             MPI_INT, mpi.comm()) == MPI_SUCCESS);
  std::vector<int> offsets(counts.size());
  int total = 0;
  for (std::size_t index = 0; index < counts.size(); ++index) {
    offsets[index] = total;
    total += counts[index];
  }
  std::vector<std::byte> gathered(static_cast<std::size_t>(total));
  HUNDUN_CHECK(MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE,
                              gathered.data(), counts.data(), offsets.data(),
                              MPI_BYTE, mpi.comm()) == MPI_SUCCESS);
  const auto *records = reinterpret_cast<const FaceRecord *>(gathered.data());
  const std::size_t record_count = gathered.size() / sizeof(FaceRecord);
  std::vector<FaceRecord> sorted(records, records + record_count);
  std::sort(sorted.begin(), sorted.end(),
            [](const auto &left, const auto &right) {
              if (left.key != right.key)
                return left.key < right.key;
              return left.global < right.global;
            });
  for (std::size_t begin = 0; begin < sorted.size();) {
    std::size_t end = begin + 1U;
    while (end < sorted.size() && sorted[end].key == sorted[begin].key)
      ++end;
    for (std::size_t index = begin + 1U; index < end; ++index) {
      HUNDUN_CHECK(sorted[index].velocity == sorted[begin].velocity);
      if (sorted[index].global == sorted[begin].global) {
        HUNDUN_CHECK(sorted[index].flux == sorted[begin].flux);
      } else {
        const double first_flux = [&] {
          double v;
          std::memcpy(&v, &sorted[begin].flux, sizeof(v));
          return v;
        }();
        const double other_flux = [&] {
          double v;
          std::memcpy(&v, &sorted[index].flux, sizeof(v));
          return v;
        }();
        HUNDUN_CHECK(bits(other_flux) == bits(-first_flux));
      }
    }
    begin = end;
  }
}

void run_periodic_case(const MpiContext &mpi, int ranks, Int3 extent,
                       bool warped) {
  constexpr double rho_ref = 1.25;
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{true, true, true},
      DecompositionOptions{process_grid_for(ranks, extent)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry =
      warped ? MeshGeometry(topology,
                            AnalyticWarpedBoxMapping(Real3{0.0, 0.0, 0.0},
                                                     Real3{1.0, 1.0, 1.0},
                                                     Real3{0.02, -0.015, 0.01}))
             : MeshGeometry(topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0},
                                                        Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(periodic_case(rho_ref), topology);
  auto interpolation = TimeConsistentFaceVelocity::create(topology, geometry);
  FieldRegistry registry;
  const Fields ids = declare_fields(registry);
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  const FieldAccessPlan access = make_access_plan(registry, ids);
  const auto stencil =
      make_momentum_time_stencil(MomentumTimeOrder::bdf2, 0.5, 1.0);
  initialize_periodic_fields(decomposition, topology, geometry, storage, ids,
                             access, rho_ref, stencil.dt_s, stencil.alpha0);
  const FieldStorage &read_storage = storage;
  const auto predictor = read_storage.view<double>(ids.predictor);
  const auto pressure = read_storage.view<double>(ids.pressure);
  const auto gradient = read_storage.view<double>(ids.pressure_gradient);
  const auto diagonal = read_storage.view<double>(ids.diagonal);
  const auto velocity_n = read_storage.view<double>(ids.velocity_n);
  const auto velocity_nm1 = read_storage.view<double>(ids.velocity_nm1);
  const auto face_n = storage.acquire_face_read<double>(access, kInitPhase,
                                                        kInitActor, ids.face_n);
  const auto face_nm1 = storage.acquire_face_read<double>(
      access, kInitPhase, kInitActor, ids.face_nm1);
  auto trial = storage.acquire_face_write<double>(
      access, kAssemblePhase, kAssembleActor, ids.trial_face);
  const MomentumFaceHistory history{velocity_n, face_n, &velocity_nm1,
                                    &face_nm1};
  interpolation.assemble_constant_density(
      boundaries, rho_ref, stencil, predictor, pressure, gradient, diagonal,
      history, trial, registry, storage, access, kAssemblePhase, kAssembleActor,
      ids.mass_flux);
  const auto flux = storage.acquire_face_read<double>(
      access, kReadPhase, kReadActor, ids.mass_flux);

  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * 16.0;
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    HUNDUN_CHECK(topology.neighbour(face).has_value());
    const OracleFace oracle = oracle_face(topology, geometry, face);
    Real3 grad{};
    for (int component = 0; component < 3; ++component) {
      const double p = value(gradient, oracle.p, component);
      const double n = value(gradient, oracle.n, component);
      const double g = interpolate(oracle.weight_p, p, oracle.weight_n, n);
      if (component == 0)
        grad.x = g;
      if (component == 1)
        grad.y = g;
      if (component == 2)
        grad.z = g;
    }
    const double defect =
        (value(pressure, oracle.n, 0) - value(pressure, oracle.p, 0) -
         dot(grad, oracle.displacement)) /
        oracle.normal_distance;
    Real3 expected{};
    for (int component = 0; component < 3; ++component) {
      const double lambda =
          oracle.weight_p *
              (value(diagonal, oracle.p, component) / oracle.volume_p) +
          oracle.weight_n *
              (value(diagonal, oracle.n, component) / oracle.volume_n);
      const double mobility = 1.0 / lambda;
      const double u =
          interpolate(oracle.weight_p, value(predictor, oracle.p, component),
                      oracle.weight_n, value(predictor, oracle.n, component));
      const double discrepancy_n =
          face_n(face, component) -
          interpolate(oracle.weight_p, value(velocity_n, oracle.p, component),
                      oracle.weight_n, value(velocity_n, oracle.n, component));
      const double discrepancy_nm1 =
          face_nm1(face, component) -
          interpolate(oracle.weight_p, value(velocity_nm1, oracle.p, component),
                      oracle.weight_n,
                      value(velocity_nm1, oracle.n, component));
      const double normal = component == 0   ? oracle.unit_normal.x
                            : component == 1 ? oracle.unit_normal.y
                                             : oracle.unit_normal.z;
      const double answer = u - mobility * defect * normal +
                            mobility * (rho_ref / stencil.dt_s) *
                                ((-stencil.alpha1) * discrepancy_n +
                                 (-stencil.alpha2) * discrepancy_nm1);
      HUNDUN_CHECK_NEAR(trial(face, component), answer, tolerance);
      if (component == 0)
        expected.x = answer;
      if (component == 1)
        expected.y = answer;
      if (component == 2)
        expected.z = answer;
    }
    const double canonical_flux = rho_ref * dot(expected, oracle.area);
    HUNDUN_CHECK_NEAR(flux(face, 0),
                      oracle.reversed ? -canonical_flux : canonical_flux,
                      tolerance);
  }
  verify_canonical_records(mpi, topology, trial, flux);

  std::size_t allocations = 0U;
  {
    allocation_probe::AllocationAttemptGuard guard;
    interpolation.assemble_constant_density(
        boundaries, rho_ref, stencil, predictor, pressure, gradient, diagonal,
        history, trial, registry, storage, access, kAssemblePhase,
        kAssembleActor, ids.mass_flux);
    allocations = guard.attempts();
  }
  HUNDUN_CHECK(allocations == 0U);

  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    for (int component = 0; component < 3; ++component) {
      trial(face, component) = -91.0;
    }
  }
  auto mass_write = storage.acquire_face_write<double>(
      access, kAssemblePhase, kAssembleActor, ids.mass_flux);
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    mass_write(face, 0) = -92.0;
  }
  expect_error([&] {
    interpolation.assemble_constant_density(
        boundaries, 0.0, stencil, predictor, pressure, gradient, diagonal,
        history, trial, registry, storage, access, kAssemblePhase,
        kAssembleActor, ids.mass_flux);
  });
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    HUNDUN_CHECK(mass_write(face, 0) == -92.0);
    for (int component = 0; component < 3; ++component) {
      HUNDUN_CHECK(trial(face, component) == -91.0);
    }
  }
  const MomentumFaceHistory missing{velocity_n, face_n, nullptr, nullptr};
  expect_error([&] {
    interpolation.assemble_constant_density(
        boundaries, rho_ref, stencil, predictor, pressure, gradient, diagonal,
        missing, trial, registry, storage, access, kAssemblePhase,
        kAssembleActor, ids.mass_flux);
  });
  auto diagonal_write = storage.view<double>(ids.diagonal);
  diagonal_write(0, 0, 0, 0) = 0.0;
  expect_error([&] {
    interpolation.assemble_constant_density(
        boundaries, rho_ref, stencil, predictor, pressure, gradient, diagonal,
        history, trial, registry, storage, access, kAssemblePhase,
        kAssembleActor, ids.mass_flux);
  });
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    HUNDUN_CHECK(mass_write(face, 0) == -92.0);
    for (int component = 0; component < 3; ++component) {
      HUNDUN_CHECK(trial(face, component) == -91.0);
    }
  }
  diagonal_write(0, 0, 0, 0) = 1.0;
  expect_error([&] {
    interpolation.assemble_constant_density(
        boundaries, rho_ref, stencil, predictor, pressure, gradient, diagonal,
        history, trial, registry, storage, access, kAssemblePhase,
        kAssembleActor, ids.trial_face);
  });
  FieldAccessPlan denied(registry);
  denied.freeze();
  expect_error([&] {
    interpolation.assemble_constant_density(
        boundaries, rho_ref, stencil, predictor, pressure, gradient, diagonal,
        history, trial, registry, storage, denied, kAssemblePhase,
        kAssembleActor, ids.mass_flux);
  });
}

void run_timestep_shrink_case(const MpiContext &mpi, int ranks) {
  constexpr Int3 extent{4, 4, 4};
  constexpr double rho_ref = 1.0;
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{true, true, true},
      DecompositionOptions{process_grid_for(ranks, extent)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(
      topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(periodic_case(rho_ref), topology);
  auto interpolation = TimeConsistentFaceVelocity::create(topology, geometry);
  FieldRegistry registry;
  const Fields ids = declare_fields(registry);
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  const FieldAccessPlan access = make_access_plan(registry, ids);
  auto pressure = storage.view<double>(ids.pressure);
  auto predictor = storage.view<double>(ids.predictor);
  auto gradient = storage.view<double>(ids.pressure_gradient);
  auto diagonal = storage.view<double>(ids.diagonal);
  auto velocity_n = storage.view<double>(ids.velocity_n);
  auto velocity_nm1 = storage.view<double>(ids.velocity_nm1);
  const auto box = topology.owned_global_box();
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        const LocalCellId local =
            (static_cast<std::size_t>(k) *
                 static_cast<std::size_t>(decomposition.local_extent().y) +
             static_cast<std::size_t>(j)) *
                static_cast<std::size_t>(decomposition.local_extent().x) +
            static_cast<std::size_t>(i);
        const double volume = geometry.cell_volume_m3(local);
        pressure(i, j, k, 0) = 0.0;
        for (int component = 0; component < 3; ++component) {
          predictor(i, j, k, component) = 0.0;
          gradient(i, j, k, component) = 0.0;
          velocity_n(i, j, k, component) = 0.0;
          velocity_nm1(i, j, k, component) = 0.0;
          diagonal(i, j, k, component) = volume;
        }
      }
    }
  }
  static_cast<void>(box);
  auto face_n = storage.acquire_face_write<double>(access, kInitPhase,
                                                   kInitActor, ids.face_n);
  auto face_nm1 = storage.acquire_face_write<double>(access, kInitPhase,
                                                     kInitActor, ids.face_nm1);
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    for (int component = 0; component < 3; ++component) {
      face_n(face, component) = component == 0 ? 1.0 : 0.0;
      face_nm1(face, component) = component == 0 ? 0.25 : 0.0;
    }
  }
  exchange_cells(decomposition, storage, ids);
  const FieldStorage &read = storage;
  const auto predictor_read = read.view<double>(ids.predictor);
  const auto pressure_read = read.view<double>(ids.pressure);
  const auto gradient_read = read.view<double>(ids.pressure_gradient);
  const auto velocity_n_read = read.view<double>(ids.velocity_n);
  const auto velocity_nm1_read = read.view<double>(ids.velocity_nm1);
  const auto face_n_read = storage.acquire_face_read<double>(
      access, kInitPhase, kInitActor, ids.face_n);
  const auto face_nm1_read = storage.acquire_face_read<double>(
      access, kInitPhase, kInitActor, ids.face_nm1);
  auto trial = storage.acquire_face_write<double>(
      access, kAssemblePhase, kAssembleActor, ids.trial_face);
  const MomentumFaceHistory history{velocity_n_read, face_n_read,
                                    &velocity_nm1_read, &face_nm1_read};

  const auto run = [&](double dt, double previous_dt, MomentumTimeOrder order) {
    const auto stencil = make_momentum_time_stencil(order, dt, previous_dt);
    auto diagonal_write = storage.view<double>(ids.diagonal);
    for (int k = 0; k < decomposition.local_extent().z; ++k) {
      for (int j = 0; j < decomposition.local_extent().y; ++j) {
        for (int i = 0; i < decomposition.local_extent().x; ++i) {
          const LocalCellId local =
              (static_cast<std::size_t>(k) *
                   static_cast<std::size_t>(decomposition.local_extent().y) +
               static_cast<std::size_t>(j)) *
                  static_cast<std::size_t>(decomposition.local_extent().x) +
              static_cast<std::size_t>(i);
          const double value =
              geometry.cell_volume_m3(local) * stencil.alpha0 * rho_ref / dt;
          for (int component = 0; component < 3; ++component) {
            diagonal_write(i, j, k, component) = value;
          }
        }
      }
    }
    auto halo = HaloExchange::create(
        decomposition,
        ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
    halo.exchange(storage, ids.diagonal);
    const auto diagonal_read = read.view<double>(ids.diagonal);
    const MomentumFaceHistory selected =
        order == MomentumTimeOrder::backward_euler
            ? MomentumFaceHistory{velocity_n_read, face_n_read, nullptr,
                                  nullptr}
            : history;
    interpolation.assemble_constant_density(
        boundaries, rho_ref, stencil, predictor_read, pressure_read,
        gradient_read, diagonal_read, selected, trial, registry, storage,
        access, kAssemblePhase, kAssembleActor, ids.mass_flux);
    double maximum = 0.0;
    for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
      maximum = std::max(maximum, std::abs(trial(face, 0)));
    }
    return maximum;
  };
  const double be_full = run(1.0, 0.0, MomentumTimeOrder::backward_euler);
  const double be_half = run(0.5, 0.0, MomentumTimeOrder::backward_euler);
  HUNDUN_CHECK_NEAR(be_full, 1.0,
                    32.0 * std::numeric_limits<double>::epsilon());
  HUNDUN_CHECK_NEAR(be_half, be_full,
                    32.0 * std::numeric_limits<double>::epsilon());
  const double bdf_full = run(1.0, 1.0, MomentumTimeOrder::bdf2);
  const double bdf_half = run(0.5, 1.0, MomentumTimeOrder::bdf2);
  HUNDUN_CHECK(bdf_full > 1.0);
  HUNDUN_CHECK(bdf_half > 0.8 * bdf_full);
}

void run_physical_case(const MpiContext &mpi, int ranks) {
  constexpr Int3 extent{4, 4, 4};
  constexpr double rho_ref = 1.5;
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{false, false, true},
      DecompositionOptions{process_grid_for(ranks, extent)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(
      topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(mixed_case(rho_ref), topology);
  auto interpolation = TimeConsistentFaceVelocity::create(topology, geometry);
  FieldRegistry registry;
  const Fields ids = declare_fields(registry);
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  const FieldAccessPlan access = make_access_plan(registry, ids);
  const auto stencil =
      make_momentum_time_stencil(MomentumTimeOrder::backward_euler, 1.0, 0.0);
  initialize_periodic_fields(decomposition, topology, geometry, storage, ids,
                             access, rho_ref, stencil.dt_s, stencil.alpha0);
  auto predictor_write = storage.view<double>(ids.predictor);
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        predictor_write(i, j, k, 0) = -1.0;
        predictor_write(i, j, k, 1) = 0.5;
        predictor_write(i, j, k, 2) = 0.25;
      }
    }
  }
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  halo.exchange(storage, ids.predictor);
  const FieldStorage &read = storage;
  const auto predictor = read.view<double>(ids.predictor);
  const auto pressure = read.view<double>(ids.pressure);
  const auto gradient = read.view<double>(ids.pressure_gradient);
  const auto diagonal = read.view<double>(ids.diagonal);
  const auto velocity_n = read.view<double>(ids.velocity_n);
  const auto face_n = storage.acquire_face_read<double>(access, kInitPhase,
                                                        kInitActor, ids.face_n);
  auto trial = storage.acquire_face_write<double>(
      access, kAssemblePhase, kAssembleActor, ids.trial_face);
  const MomentumFaceHistory history{velocity_n, face_n, nullptr, nullptr};
  interpolation.assemble_constant_density(
      boundaries, rho_ref, stencil, predictor, pressure, gradient, diagonal,
      history, trial, registry, storage, access, kAssemblePhase, kAssembleActor,
      ids.mass_flux);
  const auto flux = storage.acquire_face_read<double>(
      access, kReadPhase, kReadActor, ids.mass_flux);
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    if (topology.neighbour(face).has_value() ||
        topology.periodic_pair(face).has_value()) {
      continue;
    }
    const auto patch = topology.patch_id(face);
    HUNDUN_CHECK(patch.has_value());
    const auto owner = topology.owner(face);
    const StructuredIndex owner_index =
        map_cell(topology.global_cell(owner), topology.owned_global_box(),
                 topology.global_extent());
    const Real3 interior{value(predictor, owner_index, 0),
                         value(predictor, owner_index, 1),
                         value(predictor, owner_index, 2)};
    const Real3 area = geometry.face_area_vector_m2(face, FaceSide::owner);
    const auto expected = boundaries.evaluate_velocity(*patch, interior, area);
    HUNDUN_CHECK(bits(trial(face, 0)) == bits(expected.face.x));
    HUNDUN_CHECK(bits(trial(face, 1)) == bits(expected.face.y));
    HUNDUN_CHECK(bits(trial(face, 2)) == bits(expected.face.z));
    const auto rule = boundaries.patch(*patch).mass_flux_rule();
    if (rule == hundun::boundary::MassFluxRule::identically_zero) {
      HUNDUN_CHECK(bits(flux(face, 0)) == bits(0.0));
    } else {
      HUNDUN_CHECK(bits(flux(face, 0)) ==
                   bits(rho_ref * dot(expected.face, area)));
      if (rule == hundun::boundary::MassFluxRule::outflow_only) {
        HUNDUN_CHECK(flux(face, 0) < 0.0);
      }
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  MpiEnvironment environment(argc, argv);
  MpiContext mpi = MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    const int ranks = mpi.size();
    run_periodic_case(mpi, ranks, Int3{4, 4, 4}, true);
    run_periodic_case(mpi, ranks, Int3{1, 4, 4}, false);
    run_timestep_shrink_case(mpi, ranks);
    run_physical_case(mpi, ranks);
  });
}
