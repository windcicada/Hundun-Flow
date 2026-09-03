// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/flow_material_density_piso_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_material_density_piso.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

hundun::runtime::Int3 grid(int ranks) {
  return ranks == 1 ? hundun::runtime::Int3{1, 1, 1}
                    : ranks == 2 ? hundun::runtime::Int3{2, 1, 1}
                                 : hundun::runtime::Int3{4, 1, 1};
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<hundun::config::PatchName, 6> names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t i = 0; i < names.size(); ++i) {
    config.boundaries[i].patch = names[i];
    config.boundaries[i].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, const char *unit,
                                      std::uint32_t components,
                                      bool conservative) {
  return {name, unit, "task20-vortex",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64, components, 2, conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face(const char *name,
                                      std::uint32_t components) {
  return {name, "m/s", "task20-vortex",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64, components, 0, false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::Real3 velocity(double x, double y) {
  return {std::sin(x) * std::cos(y), -std::cos(x) * std::sin(y), 0.0};
}

double density(double x, double y) {
  return 1.0 + 0.1 * std::sin(x) * std::sin(y);
}

double density_cell_average(double x, double y, double spacing) {
  const double factor = std::sin(0.5 * spacing) / (0.5 * spacing);
  return 1.0 + 0.1 * factor * factor * std::sin(x) * std::sin(y);
}

std::uint64_t bits(double value) {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

struct Errors final {
  double rho{};
  double rho_u{};
  double rho_v{};
  double rho_w{};
  double pressure{};
  double pressure_mean{};
  double parity_pressure{};
  double parity_pressure_mean{};
};

struct VortexReference final {
  std::uint32_t resolution{};
  std::vector<double> cell_fields;
  std::vector<double> face_mass_flux;
};

struct CaseResult final {
  Errors errors;
  VortexReference reference;
};

constexpr std::array<unsigned char, 8> reference_magic{'H', 'V', 'T', '2',
                                                       '0', 'R', '1', 0};
constexpr std::uint32_t reference_version = 1U;
constexpr std::uint64_t fnv_offset = UINT64_C(14695981039346656037);
constexpr std::uint64_t fnv_prime = UINT64_C(1099511628211);

std::uint64_t fnv1a(const unsigned char *data, std::size_t size,
                    std::uint64_t hash = fnv_offset) noexcept {
  for (std::size_t index = 0; index < size; ++index) {
    hash ^= data[index];
    hash *= fnv_prime;
  }
  return hash;
}

std::uint64_t executable_fingerprint(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input)
    throw std::runtime_error("cannot open vortex test executable");
  std::array<char, 8192> buffer{};
  std::uint64_t hash = fnv_offset;
  while (input) {
    input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const auto count = input.gcount();
    if (count > 0)
      hash = fnv1a(reinterpret_cast<const unsigned char *>(buffer.data()),
                   static_cast<std::size_t>(count), hash);
  }
  if (!input.eof())
    throw std::runtime_error("cannot read vortex test executable");
  return hash;
}

void append_u32(std::vector<unsigned char> &bytes, std::uint32_t value) {
  for (unsigned shift = 0; shift < 32U; shift += 8U)
    bytes.push_back(static_cast<unsigned char>((value >> shift) & 0xffU));
}

void append_u64(std::vector<unsigned char> &bytes, std::uint64_t value) {
  for (unsigned shift = 0; shift < 64U; shift += 8U)
    bytes.push_back(static_cast<unsigned char>((value >> shift) & 0xffU));
}

void append_double(std::vector<unsigned char> &bytes, double value) {
  append_u64(bytes, bits(value));
}

class ReferenceReader final {
public:
  explicit ReferenceReader(const std::vector<unsigned char> &bytes)
      : bytes_(&bytes) {}

  unsigned char byte() {
    require(1U);
    return (*bytes_)[offset_++];
  }
  std::uint32_t u32() {
    std::uint32_t result{};
    for (unsigned shift = 0; shift < 32U; shift += 8U)
      result |= static_cast<std::uint32_t>(byte()) << shift;
    return result;
  }
  std::uint64_t u64() {
    std::uint64_t result{};
    for (unsigned shift = 0; shift < 64U; shift += 8U)
      result |= static_cast<std::uint64_t>(byte()) << shift;
    return result;
  }
  double fp64() {
    const auto value = u64();
    double result{};
    std::memcpy(&result, &value, sizeof(result));
    return result;
  }
  std::size_t offset() const noexcept { return offset_; }

private:
  void require(std::size_t count) const {
    if (count > bytes_->size() - std::min(bytes_->size(), offset_))
      throw std::runtime_error("vortex reference artifact is truncated");
  }
  const std::vector<unsigned char> *bytes_;
  std::size_t offset_{};
};

void write_reference_artifact(const std::filesystem::path &path,
                              std::uint64_t executable,
                              const std::array<CaseResult, 3> &cases) {
  std::vector<unsigned char> bytes;
  std::size_t reserve = 64U;
  for (const auto &result : cases)
    reserve += result.reference.cell_fields.size() * sizeof(double) +
               result.reference.face_mass_flux.size() *
                   (sizeof(std::uint64_t) + sizeof(double));
  bytes.reserve(reserve);
  bytes.insert(bytes.end(), reference_magic.begin(), reference_magic.end());
  append_u32(bytes, reference_version);
  append_u64(bytes, executable);
  append_u32(bytes, static_cast<std::uint32_t>(cases.size()));
  for (const auto &result : cases) {
    const auto &record = result.reference;
    if (record.cell_fields.size() % 5U != 0U)
      throw std::runtime_error("vortex reference cell layout is invalid");
    const std::uint64_t cells = record.cell_fields.size() / 5U;
    const std::uint64_t faces = record.face_mass_flux.size();
    append_u32(bytes, record.resolution);
    append_u64(bytes, cells);
    append_u64(bytes, faces);
    for (std::uint64_t cell = 0; cell < cells; ++cell) {
      append_u64(bytes, cell);
      for (std::size_t component = 0; component < 5U; ++component)
        append_double(bytes, record.cell_fields[cell * 5U + component]);
    }
    for (std::uint64_t face = 0; face < faces; ++face) {
      append_u64(bytes, face);
      append_double(bytes, record.face_mass_flux[face]);
    }
  }
  append_u64(bytes, fnv1a(bytes.data(), bytes.size()));
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output)
    throw std::runtime_error("cannot create vortex reference artifact");
  output.write(reinterpret_cast<const char *>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
  output.close();
  if (!output)
    throw std::runtime_error("cannot publish vortex reference artifact");
}

std::array<VortexReference, 3>
read_reference_artifact(const std::filesystem::path &path,
                        std::uint64_t executable,
                        const std::array<int, 3> &resolutions) {
  std::ifstream input(path, std::ios::binary);
  if (!input)
    throw std::runtime_error("vortex reference artifact is missing");
  const std::vector<unsigned char> bytes{std::istreambuf_iterator<char>(input),
                                         std::istreambuf_iterator<char>()};
  if (input.bad())
    throw std::runtime_error("cannot read vortex reference artifact");
  if (bytes.size() < reference_magic.size() + 4U + 8U + 4U + 8U)
    throw std::runtime_error("vortex reference artifact is truncated");
  std::uint64_t stored_checksum{};
  for (unsigned shift = 0; shift < 64U; shift += 8U)
    stored_checksum |=
        static_cast<std::uint64_t>(bytes[bytes.size() - 8U + shift / 8U])
        << shift;
  if (stored_checksum != fnv1a(bytes.data(), bytes.size() - 8U))
    throw std::runtime_error("vortex reference artifact checksum differs");

  ReferenceReader reader(bytes);
  for (const auto expected : reference_magic)
    if (reader.byte() != expected)
      throw std::runtime_error("vortex reference artifact magic differs");
  if (reader.u32() != reference_version)
    throw std::runtime_error("vortex reference artifact version differs");
  if (reader.u64() != executable)
    throw std::runtime_error("vortex reference artifact executable differs");
  if (reader.u32() != resolutions.size())
    throw std::runtime_error("vortex reference artifact record count differs");
  std::array<VortexReference, 3> result{};
  for (std::size_t record_index = 0; record_index < result.size();
       ++record_index) {
    auto &record = result[record_index];
    record.resolution = reader.u32();
    if (record.resolution !=
        static_cast<std::uint32_t>(resolutions[record_index]))
      throw std::runtime_error("vortex reference artifact resolution differs");
    const auto cells = reader.u64();
    const auto faces = reader.u64();
    const std::uint64_t expected_cells =
        static_cast<std::uint64_t>(record.resolution) * record.resolution * 4U;
    const std::uint64_t n = record.resolution;
    const std::uint64_t expected_faces =
        2U * (n + 1U) * n * 4U + n * n * 5U;
    if (cells != expected_cells || faces != expected_faces ||
        cells > SIZE_MAX / 5U || faces > SIZE_MAX)
      throw std::runtime_error("vortex reference artifact counts differ");
    record.cell_fields.resize(static_cast<std::size_t>(cells) * 5U);
    record.face_mass_flux.resize(static_cast<std::size_t>(faces));
    for (std::uint64_t cell = 0; cell < cells; ++cell) {
      if (reader.u64() != cell)
        throw std::runtime_error(
            "vortex reference cell IDs are duplicate, missing or unordered");
      for (std::size_t component = 0; component < 5U; ++component) {
        const double value = reader.fp64();
        if (!std::isfinite(value))
          throw std::runtime_error("vortex reference cell value is non-finite");
        record.cell_fields[cell * 5U + component] = value;
      }
    }
    for (std::uint64_t face = 0; face < faces; ++face) {
      if (reader.u64() != face)
        throw std::runtime_error(
            "vortex reference face IDs are duplicate, missing or unordered");
      const double value = reader.fp64();
      if (!std::isfinite(value))
        throw std::runtime_error("vortex reference face value is non-finite");
      record.face_mass_flux[face] = value;
    }
  }
  if (reader.offset() != bytes.size() - 8U)
    throw std::runtime_error("vortex reference artifact has trailing bytes");
  return result;
}

std::vector<unsigned char> artifact_bytes(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input)
    throw std::runtime_error("cannot open vortex artifact mutation input");
  std::vector<unsigned char> result{std::istreambuf_iterator<char>(input),
                                    std::istreambuf_iterator<char>()};
  if (input.bad())
    throw std::runtime_error("cannot read vortex artifact mutation input");
  return result;
}

void write_artifact_bytes(const std::filesystem::path &path,
                          const std::vector<unsigned char> &bytes) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output)
    throw std::runtime_error("cannot create vortex artifact mutation");
  output.write(reinterpret_cast<const char *>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
  output.close();
  if (!output)
    throw std::runtime_error("cannot write vortex artifact mutation");
}

void replace_artifact_checksum(std::vector<unsigned char> &bytes) {
  if (bytes.size() < 8U)
    throw std::runtime_error("vortex artifact mutation is too short");
  bytes.resize(bytes.size() - 8U);
  append_u64(bytes, fnv1a(bytes.data(), bytes.size()));
}

bool artifact_rejected(const std::filesystem::path &path,
                       std::uint64_t executable,
                       const std::array<int, 3> &resolutions) noexcept {
  try {
    static_cast<void>(read_reference_artifact(path, executable, resolutions));
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

bool decomposition_equal(const VortexReference &reference,
                         const VortexReference &candidate) noexcept {
  if (reference.resolution != candidate.resolution ||
      reference.cell_fields.size() != candidate.cell_fields.size() ||
      reference.face_mass_flux.size() != candidate.face_mass_flux.size() ||
      reference.cell_fields.size() % 5U != 0U)
    return false;
  for (std::size_t component = 0; component < 5U; ++component) {
    double maximum_difference{};
    double infinity_norm{};
    for (std::size_t cell = 0; cell < reference.cell_fields.size() / 5U;
         ++cell) {
      const double expected = reference.cell_fields[cell * 5U + component];
      const double observed = candidate.cell_fields[cell * 5U + component];
      if (!std::isfinite(expected) || !std::isfinite(observed))
        return false;
      maximum_difference =
          std::max(maximum_difference, std::abs(observed - expected));
      infinity_norm = std::max(infinity_norm, std::abs(expected));
    }
    if (maximum_difference > 5.0e-12 * std::max(1.0, infinity_norm))
      return false;
  }
  double maximum_difference{};
  double infinity_norm{};
  for (std::size_t face = 0; face < reference.face_mass_flux.size(); ++face) {
    const double expected = reference.face_mass_flux[face];
    const double observed = candidate.face_mass_flux[face];
    if (!std::isfinite(expected) || !std::isfinite(observed))
      return false;
    maximum_difference =
        std::max(maximum_difference, std::abs(observed - expected));
    infinity_norm = std::max(infinity_norm, std::abs(expected));
  }
  return maximum_difference <= 5.0e-12 * std::max(1.0, infinity_norm);
}
} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);

    const auto direct =
        hundun::flow::test::MaterialDensityPisoTestAccess::vortex_source(
            0.25 * pi, 0.5 * pi, 0.01);
    HUNDUN_CHECK_NEAR(
        direct.x,
        density(0.25 * pi, 0.5 * pi) * std::sin(0.25 * pi) *
                std::cos(0.25 * pi) +
            0.02 * std::sin(0.25 * pi) * std::cos(0.5 * pi),
        64.0 * std::numeric_limits<double>::epsilon());
    HUNDUN_CHECK(direct.z == 0.0);

    {
      constexpr double predictor_velocity_interpolation = 0.4;
      constexpr double predictor_momentum_interpolation = 0.72;
      constexpr double face_density = 1.2;
      constexpr double cell_momentum_n = 0.66;
      constexpr double face_momentum_n = 0.9;
      constexpr double cell_momentum_old = -0.15;
      constexpr double face_momentum_old = 0.05;
      constexpr double mobility = 0.08;
      constexpr double dt_s = 0.2;
      constexpr double alpha1 = -2.0;
      constexpr double alpha2 = 0.5;
      constexpr double pressure_correction = -0.03;
      const double value =
          hundun::flow::test::MaterialDensityPisoTestAccess::
              material_face_value(
                  predictor_velocity_interpolation,
                  predictor_momentum_interpolation, face_density,
                  cell_momentum_n, face_momentum_n, cell_momentum_old,
                  face_momentum_old, mobility, dt_s, alpha1, alpha2,
                  pressure_correction);
      const double delta_n = face_momentum_n - cell_momentum_n;
      const double delta_old = face_momentum_old - cell_momentum_old;
      const double time_term =
          (mobility / dt_s) *
          ((-alpha1) * delta_n + (-alpha2) * delta_old);
      const double expected = predictor_momentum_interpolation / face_density +
                              pressure_correction + time_term;
      HUNDUN_CHECK_NEAR(value, expected,
                        32.0 * std::numeric_limits<double>::epsilon());
      const double omitted = predictor_momentum_interpolation / face_density +
                             pressure_correction;
      const double plain_velocity_base = predictor_velocity_interpolation +
                                         pressure_correction + time_term;
      const double wrong_delta_n =
          face_density *
          (face_momentum_n / face_density -
           predictor_velocity_interpolation);
      const double wrong_delta_old =
          face_density *
          (face_momentum_old / face_density - (-0.1));
      const double wrong_discrepancy =
          predictor_momentum_interpolation / face_density +
          pressure_correction +
          (mobility / dt_s) *
              ((-alpha1) * wrong_delta_n + (-alpha2) * wrong_delta_old);
      HUNDUN_CHECK(std::abs(value - omitted) > 1.0e-3);
      HUNDUN_CHECK(std::abs(value - plain_velocity_base) > 1.0e-3);
      HUNDUN_CHECK(std::abs(value - wrong_discrepancy) > 1.0e-3);
    }

    const auto run_case = [&](int n) {
    const hundun::runtime::Int3 extent{n, n, 4};
    auto decomposition = hundun::runtime::StructuredDecomposition::create(
        mpi, extent, {true, true, true},
        hundun::runtime::DecompositionOptions{grid(mpi.size())});
    hundun::mesh::MeshTopology topology(decomposition);
    hundun::mesh::MeshGeometry geometry(
        topology, hundun::mesh::UniformBoxMapping(
                      {0.0, 0.0, 0.0}, {2.0 * pi, 2.0 * pi, 1.0}));
    auto boundaries = hundun::boundary::BoundaryRegistry::create(
        periodic_case(), topology);
    hundun::runtime::FieldRegistry registry;
    hundun::flow::FlowFieldIds fields;
    fields.density = registry.declare_field(cell("rho", "kg/m3", 1U, true));
    fields.velocity =
        registry.declare_field(cell("velocity", "m/s", 3U, false));
    fields.mechanical_pressure =
        registry.declare_field(cell("pi", "Pa", 1U, false));
    fields.face_velocity = registry.declare_field(face("face_velocity", 3U));
    fields.face_mass_flux =
        hundun::finite_volume::declare_face_mass_flux(registry);
    const auto rho_h =
        registry.declare_field(cell("rho_h", "J/m3", 1U, true));
    const auto initial_face_density =
        registry.declare_field(face("initial_face_density", 1U));
    fields.transported_cell_fields = {rho_h};
    registry.freeze();

    const auto box = decomposition.owned_box();
    const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                      box.end.y - box.begin.y,
                                      box.end.z - box.begin.z};
    auto halo = hundun::runtime::HaloExchange::create(
        decomposition,
        hundun::runtime::ExchangePlan::create(decomposition, local, 2));
    auto state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 1.0e-4, 1.0e-4,
         hundun::flow::MomentumTimeOrder::bdf2});
    hundun::flow::FlowLayerValues exact;
    exact.density.resize(topology.owned_cell_count());
    exact.velocity.resize(topology.owned_cell_count() * 3U);
    exact.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
    exact.transported_cell_fields.resize(1U);
    exact.transported_cell_fields[0].resize(topology.owned_cell_count());
    for (hundun::mesh::LocalCellId cell_id = 0;
         cell_id < topology.owned_cell_count(); ++cell_id) {
      const auto centre = geometry.cell_center_m(cell_id);
      const auto u = velocity(centre.x, centre.y);
      const double rho =
          n == 8 ? density(centre.x, centre.y)
                 : density_cell_average(
                       centre.x, centre.y,
                       2.0 * pi / static_cast<double>(n));
      exact.density[cell_id] = rho;
      exact.velocity[cell_id * 3U] = u.x;
      exact.velocity[cell_id * 3U + 1U] = u.y;
      exact.velocity[cell_id * 3U + 2U] = 0.0;
      exact.transported_cell_fields[0][cell_id] = rho;
    }
    exact.face_velocity.resize(topology.local_face_count() * 3U);
    exact.face_mass_flux.resize(topology.local_face_count());
    for (hundun::mesh::LocalFaceId face_id = 0;
         face_id < topology.local_face_count(); ++face_id) {
      const auto centre = geometry.face_center_m(face_id);
      const auto u = velocity(centre.x, centre.y);
      exact.face_velocity[face_id * 3U] = u.x;
      exact.face_velocity[face_id * 3U + 1U] = u.y;
      exact.face_velocity[face_id * 3U + 2U] = 0.0;
      const auto area =
          geometry.face_area_vector_m2(face_id, hundun::mesh::FaceSide::owner);
      exact.face_mass_flux[face_id] =
          u.x * area.x + u.y * area.y + u.z * area.z;
    }
    constexpr hundun::runtime::PhaseId initial_phase = 2020U;
    constexpr hundun::runtime::ActorId initial_actor = 2020U;
    hundun::runtime::FieldAccessPlan initial_access(registry);
    initial_access.declare_access(initial_phase, initial_actor, fields.density,
                                  hundun::runtime::AccessMode::read_write);
    initial_access.declare_access(initial_phase, initial_actor,
                                  fields.face_mass_flux,
                                  hundun::runtime::AccessMode::read_write);
    initial_access.declare_access(initial_phase, initial_actor,
                                  initial_face_density,
                                  hundun::runtime::AccessMode::read_write);
    initial_access.freeze();
    hundun::runtime::FieldStorage initial_storage(
        registry, {local, topology.local_face_count()});
    {
      auto rho = initial_storage.acquire_write<double>(
          initial_access, initial_phase, initial_actor, fields.density);
      for (hundun::mesh::LocalCellId cell_id = 0;
           cell_id < topology.owned_cell_count(); ++cell_id) {
        const auto global = topology.global_cell(cell_id);
        rho(global.x - box.begin.x, global.y - box.begin.y,
            global.z - box.begin.z, 0) = exact.density[cell_id];
      }
      auto direction = initial_storage.acquire_face_write<double>(
          initial_access, initial_phase, initial_actor, fields.face_mass_flux);
      for (hundun::mesh::LocalFaceId face_id = 0;
           face_id < topology.local_face_count(); ++face_id)
        direction(face_id, 0) = exact.face_mass_flux[face_id];
    }
    halo.exchange(initial_storage, fields.density);
    {
      const auto direction = hundun::finite_volume::FaceMassFlux::acquire(
          registry, initial_storage, initial_access, initial_phase,
          initial_actor, fields.face_mass_flux, topology);
      const auto rho = initial_storage.acquire_read<double>(
          initial_access, initial_phase, initial_actor, fields.density);
      auto rho_face = initial_storage.acquire_face_write<double>(
          initial_access, initial_phase, initial_actor, initial_face_density);
      auto fvm = hundun::finite_volume::CellCenteredFvmOperators::create(
          topology, geometry);
      fvm.reconstruct_transport_faces(
          hundun::finite_volume::FiniteVolumeQuantity::density(), boundaries,
          direction, rho, rho_face);
      for (hundun::mesh::LocalFaceId face_id = 0;
           face_id < topology.local_face_count(); ++face_id) {
        exact.face_mass_flux[face_id] *= rho_face(face_id, 0);
      }
    }
    state.seed_accepted_layers(exact, exact);

    hundun::execution::CpuReferenceContext execution;
    hundun::linear::ConjugateGradientSolver momentum_solver(execution, mpi);
    hundun::linear::ConjugateGradientSolver pressure_solver(execution, mpi);
    hundun::linear::JacobiPreconditioner mx(execution), my(execution),
        mz(execution), pressure_preconditioner(execution);
    hundun::flow::MaterialDensityTransportSpec specification;
    specification.enthalpy_density = rho_h;
    specification.enthalpy_diffusivity_kg_per_m_s = 0.0;
    auto flow = hundun::flow::FixedStepMaterialDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_preconditioner, registry, fields, specification);
    hundun::flow::test::MaterialDensityPisoTestAccess::enable_vortex_source(
        flow, true);
    const auto stencil = hundun::flow::make_momentum_time_stencil(
        hundun::flow::MomentumTimeOrder::bdf2, 1.0e-4, 1.0e-4);
    const auto report = flow.attempt(state, 0.01, stencil, {}, {});
    if (mpi.rank() == 0 &&
        report.flow().disposition !=
            hundun::flow::StepAttemptDisposition::committed) {
      std::cerr << "Task20 vortex failure disposition="
                << static_cast<int>(report.flow().disposition)
                << " reason=" << static_cast<int>(report.flow().reason)
                << " correctors=" << report.flow().pressure_corrector_count
                << " pressure-solve="
                << static_cast<int>(report.flow().pressure[0].reason)
                << " pressure-final="
                << report.flow().pressure[0].final_residual
                << " pressure-initial="
                << report.flow().pressure[0].initial_residual
                << " pressure-iterations="
                << report.flow().pressure[0].iterations
                << " material="
                << static_cast<int>(report.material_failure_reason())
                << " continuity="
                << report.flow().final_continuity_normalized_l2
                << " pressure=" << report.final_pressure_normalized_residual()
                << " momentum=" << report.flow().final_momentum_normalized_l2[0]
                << ',' << report.flow().final_momentum_normalized_l2[1] << ','
                << report.flow().final_momentum_normalized_l2[2]
                << " mass-defect="
                << report.flow().final_mass_relative_conservation_defect
                << " material-mass="
                << (report.material_report_available()
                        ? report.material_report()
                              .mass_relative_conservation_defect()
                        : -1.0)
                << " momentum-defect="
                << report.flow().final_momentum_relative_conservation_defect[0]
                << ','
                << report.flow().final_momentum_relative_conservation_defect[1]
                << ','
                << report.flow().final_momentum_relative_conservation_defect[2]
                << '\n';
    }
    HUNDUN_CHECK(report.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(report.flow().pressure_corrector_count == 2U);
    HUNDUN_CHECK(report.final_continuity_residual_available());
    HUNDUN_CHECK(report.flow().final_continuity_normalized_l2 <= 1.0e-10);
    HUNDUN_CHECK(report.final_pressure_residual_available());
    HUNDUN_CHECK(report.final_pressure_normalized_residual() <= 1.0);
    for (double residual : report.flow().final_momentum_normalized_l2)
      HUNDUN_CHECK(residual <= 1.0e-9);
    HUNDUN_CHECK(report.material_report().disposition() ==
                 hundun::flow::MaterialTransportDisposition::finalized);
    HUNDUN_CHECK(report.flux_provenance() ==
                 hundun::flow::MaterialFluxProvenance::final_corrected);
    HUNDUN_CHECK(report.shared_face_mass_flux_field() ==
                 fields.face_mass_flux);
    HUNDUN_CHECK(report.material_report().density_normalized_l2() <= 1.0e-10);
    for (const auto residual :
         report.material_report().transport_normalized_l2())
      HUNDUN_CHECK(residual <= 1.0e-9);
    HUNDUN_CHECK(report.flow().final_mass_relative_conservation_defect <=
                 5.0e-12);
    HUNDUN_CHECK(
        report.material_report().mass_relative_conservation_defect() <=
        5.0e-12);
    for (const auto defect :
         report.flow().final_momentum_relative_conservation_defect)
      HUNDUN_CHECK(defect <= 5.0e-11);
    HUNDUN_CHECK(report.material_report().minimum_density_kg_per_m3() > 0.0);
    const auto committed = state.snapshot(hundun::flow::FlowLayer::committed);
    const auto &finalizer_input =
        hundun::flow::test::MaterialDensityPisoTestAccess::
            finalizer_flux_evidence(flow);
    HUNDUN_CHECK(finalizer_input.size() == committed.face_mass_flux.size());
    for (std::size_t face_index = 0; face_index < finalizer_input.size();
         ++face_index)
      HUNDUN_CHECK(bits(finalizer_input[face_index]) ==
                   bits(committed.face_mass_flux[face_index]));

    const std::size_t global_cells = static_cast<std::size_t>(n) *
                                     static_cast<std::size_t>(n) * 4U;
    std::vector<double> cell_cover(global_cells, 0.0);
    std::vector<double> global_fields(global_cells * 5U, 0.0);
    for (hundun::mesh::LocalCellId cell_id = 0;
         cell_id < topology.owned_cell_count(); ++cell_id) {
      const auto global =
          static_cast<std::size_t>(topology.global_cell_id(cell_id));
      HUNDUN_CHECK(global < global_cells);
      cell_cover[global] = 1.0;
      global_fields[global * 5U] = committed.density[cell_id];
      global_fields[global * 5U + 1U] =
          committed.density[cell_id] * committed.velocity[cell_id * 3U];
      global_fields[global * 5U + 2U] =
          committed.density[cell_id] * committed.velocity[cell_id * 3U + 1U];
      global_fields[global * 5U + 3U] =
          committed.density[cell_id] * committed.velocity[cell_id * 3U + 2U];
      global_fields[global * 5U + 4U] =
          committed.mechanical_pressure[cell_id];
    }
    mpi.allreduce_fp64_in_place(cell_cover.data(), cell_cover.size(),
                                hundun::runtime::Fp64ReductionOperation::sum);
    mpi.allreduce_fp64_in_place(global_fields.data(), global_fields.size(),
                                hundun::runtime::Fp64ReductionOperation::sum);
    HUNDUN_CHECK(std::all_of(cell_cover.begin(), cell_cover.end(),
                            [](double value) { return value == 1.0; }));
    std::vector<double> face_cover(topology.global_face_count(), 0.0);
    std::vector<double> global_face_flux(topology.global_face_count(), 0.0);
    for (hundun::mesh::LocalFaceId face_id = 0;
         face_id < topology.local_face_count(); ++face_id) {
      if (topology.cell_ownership(topology.owner(face_id)) !=
          hundun::mesh::EntityOwnership::owned)
        continue;
      const auto global =
          static_cast<std::size_t>(topology.global_face_id(face_id));
      HUNDUN_CHECK(global < face_cover.size());
      face_cover[global] = 1.0;
      global_face_flux[global] = committed.face_mass_flux[face_id];
    }
    mpi.allreduce_fp64_in_place(face_cover.data(), face_cover.size(),
                                hundun::runtime::Fp64ReductionOperation::sum);
    mpi.allreduce_fp64_in_place(
        global_face_flux.data(), global_face_flux.size(),
        hundun::runtime::Fp64ReductionOperation::sum);
    HUNDUN_CHECK(std::all_of(face_cover.begin(), face_cover.end(),
                            [](double value) { return value == 1.0; }));
    std::array<double, 12> sums{};
    for (hundun::mesh::LocalCellId cell_id = 0;
         cell_id < topology.owned_cell_count(); ++cell_id) {
      const auto centre = geometry.cell_center_m(cell_id);
      const auto exact_u = velocity(centre.x, centre.y);
      const double exact_rho = density(centre.x, centre.y);
      const double volume = geometry.cell_volume_m3(cell_id);
      const double numeric_rho = committed.density[cell_id];
      const double numeric_u = committed.velocity[cell_id * 3U];
      const double numeric_v = committed.velocity[cell_id * 3U + 1U];
      const double numeric_w = committed.velocity[cell_id * 3U + 2U];
      const double rho_error = numeric_rho - exact_rho;
      const double rho_u_error = numeric_rho * numeric_u - exact_rho * exact_u.x;
      const double rho_v_error = numeric_rho * numeric_v - exact_rho * exact_u.y;
      const double rho_w_error = numeric_rho * numeric_w;
      const double pressure_value = committed.mechanical_pressure[cell_id];
      const auto global = topology.global_cell(cell_id);
      const double parity = ((global.x + global.y + global.z) & 1) == 0
                                ? 1.0
                                : -1.0;
      const double mutated_pressure = pressure_value + parity;
      sums[0] += volume * rho_error * rho_error;
      sums[1] += volume * exact_rho * exact_rho;
      sums[2] += volume * rho_u_error * rho_u_error;
      sums[3] += volume * (exact_rho * exact_u.x) * (exact_rho * exact_u.x);
      sums[4] += volume * rho_v_error * rho_v_error;
      sums[5] += volume * (exact_rho * exact_u.y) * (exact_rho * exact_u.y);
      sums[6] += volume * rho_w_error * rho_w_error;
      sums[7] += volume * pressure_value * pressure_value;
      sums[8] += volume;
      sums[9] += volume * pressure_value;
      sums[10] += volume * mutated_pressure * mutated_pressure;
      sums[11] += volume * mutated_pressure;
    }
    mpi.allreduce_fp64_in_place(sums.data(), sums.size(),
                                hundun::runtime::Fp64ReductionOperation::sum);
    return CaseResult{
        Errors{std::sqrt(sums[0] / sums[1]), std::sqrt(sums[2] / sums[3]),
               std::sqrt(sums[4] / sums[5]), std::sqrt(sums[6] / sums[8]),
               std::sqrt(sums[7] / sums[8]), sums[9] / sums[8],
               std::sqrt(sums[10] / sums[8]), sums[11] / sums[8]},
        VortexReference{static_cast<std::uint32_t>(n),
                        std::move(global_fields),
                        std::move(global_face_flux)}};
    };

    const bool full = argc > 1 && std::string_view(argv[1]) == "--full";
    if (!full) {
      static_cast<void>(run_case(8));
      return;
    }
    const std::array<int, 3> resolutions{16, 32, 64};
    std::array<CaseResult, 3> cases{};
    for (std::size_t index = 0; index < resolutions.size(); ++index)
      cases[index] = run_case(resolutions[index]);
    const auto require_order = [&](auto select) {
      for (std::size_t index = 0; index + 1U < cases.size(); ++index) {
        const double coarse = select(cases[index].errors);
        const double fine = select(cases[index + 1U].errors);
        HUNDUN_CHECK(std::isfinite(coarse) && coarse > 0.0);
        HUNDUN_CHECK(std::isfinite(fine) && fine > 0.0 && fine < coarse);
        HUNDUN_CHECK(std::log(coarse / fine) / std::log(2.0) >= 1.8);
      }
    };
    require_order([](const Errors &value) { return value.rho; });
    require_order([](const Errors &value) { return value.rho_u; });
    require_order([](const Errors &value) { return value.rho_v; });
    require_order([](const Errors &value) { return value.pressure; });
    for (std::size_t index = 0; index < cases.size(); ++index) {
      const double h = 2.0 * pi / static_cast<double>(resolutions[index]);
      HUNDUN_CHECK(cases[index].errors.rho_w <= 1.0e-12);
      HUNDUN_CHECK(cases[index].errors.pressure <= 2.0 * 1.1 * h * h);
      HUNDUN_CHECK(std::abs(cases[index].errors.pressure_mean) <= 1.0e-12);
      HUNDUN_CHECK(cases[index].errors.parity_pressure >
                   2.0 * 1.1 * h * h);
      HUNDUN_CHECK(
          std::abs(cases[index].errors.parity_pressure_mean) <= 1.0e-12);
    }
    HUNDUN_CHECK(argc == 4);
    HUNDUN_CHECK(
        (mpi.size() == 1 && std::string_view(argv[2]) == "--reference-write") ||
        (mpi.size() > 1 && std::string_view(argv[2]) == "--reference-read"));
    const auto executable = executable_fingerprint(argv[0]);
    const std::filesystem::path artifact(argv[3]);
    if (mpi.size() == 1) {
      write_reference_artifact(artifact, executable, cases);
      const auto round_trip =
          read_reference_artifact(artifact, executable, resolutions);
      for (std::size_t index = 0; index < cases.size(); ++index)
        HUNDUN_CHECK(
            decomposition_equal(round_trip[index], cases[index].reference));
      HUNDUN_CHECK(
          artifact_rejected(artifact, executable ^ UINT64_C(1), resolutions));
      const auto missing = artifact.string() + ".missing";
      std::filesystem::remove(missing);
      HUNDUN_CHECK(artifact_rejected(missing, executable, resolutions));

      const auto original = artifact_bytes(artifact);
      const auto reject_mutation = [&](std::string_view name,
                                       std::vector<unsigned char> bytes) {
        const std::filesystem::path candidate =
            artifact.string() + ".reject-" + std::string(name);
        write_artifact_bytes(candidate, bytes);
        const bool rejected =
            artifact_rejected(candidate, executable, resolutions);
        std::filesystem::remove(candidate);
        HUNDUN_CHECK(rejected);
      };
      auto truncated = original;
      truncated.pop_back();
      reject_mutation("truncated", std::move(truncated));
      auto corrupt = original;
      corrupt[50U] ^= 0x1U;
      reject_mutation("checksum", std::move(corrupt));
      auto trailing = original;
      trailing.resize(trailing.size() - 8U);
      trailing.push_back(0x5aU);
      append_u64(trailing, fnv1a(trailing.data(), trailing.size()));
      reject_mutation("trailing", std::move(trailing));
      auto mismatched_resolution = original;
      mismatched_resolution[24U] ^= 0x1U;
      replace_artifact_checksum(mismatched_resolution);
      reject_mutation("resolution", std::move(mismatched_resolution));
      auto mismatched_count = original;
      mismatched_count[28U] ^= 0x1U;
      replace_artifact_checksum(mismatched_count);
      reject_mutation("cell-count", std::move(mismatched_count));
      auto mismatched_face_count = original;
      mismatched_face_count[36U] ^= 0x1U;
      replace_artifact_checksum(mismatched_face_count);
      reject_mutation("face-count", std::move(mismatched_face_count));
      auto duplicate_id = original;
      duplicate_id[44U] = 1U;
      replace_artifact_checksum(duplicate_id);
      reject_mutation("duplicate-id", std::move(duplicate_id));
      auto non_finite = original;
      constexpr std::uint64_t quiet_nan = UINT64_C(0x7ff8000000000000);
      for (unsigned shift = 0; shift < 64U; shift += 8U)
        non_finite[52U + shift / 8U] =
            static_cast<unsigned char>((quiet_nan >> shift) & 0xffU);
      replace_artifact_checksum(non_finite);
      reject_mutation("non-finite", std::move(non_finite));
    } else {
      const auto references =
          read_reference_artifact(artifact, executable, resolutions);
      for (std::size_t index = 0; index < cases.size(); ++index)
        HUNDUN_CHECK(
            decomposition_equal(references[index], cases[index].reference));
    }

    HUNDUN_CHECK(
        decomposition_equal(cases.front().reference, cases.front().reference));
    auto changed_cell = cases.front().reference;
    double cell_scale = 1.0;
    for (double value : changed_cell.cell_fields)
      cell_scale = std::max(cell_scale, std::abs(value));
    changed_cell.cell_fields.front() += 1.0e-9 * cell_scale;
    HUNDUN_CHECK(!decomposition_equal(cases.front().reference, changed_cell));
    auto changed_face = cases.front().reference;
    double face_scale = 1.0;
    for (double value : changed_face.face_mass_flux)
      face_scale = std::max(face_scale, std::abs(value));
    changed_face.face_mass_flux.front() += 1.0e-9 * face_scale;
    HUNDUN_CHECK(!decomposition_equal(cases.front().reference, changed_face));
  });
}
