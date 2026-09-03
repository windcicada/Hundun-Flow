// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/flow_ideal_gas_closure_test_access.hpp"
#include "tests/support/flow_material_density_piso_test_access.hpp"
#include "tests/support/flow_material_density_transport_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_ideal_gas_piso.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_mpi_operation_error.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace task21_allocation_oracle {

std::atomic<const hundun::flow::FixedStepIdealGasFlow *> observed_flow{};
std::atomic<std::uint64_t> allocation_calls{};

void record() noexcept {
  const auto *flow = observed_flow.load(std::memory_order_relaxed);
  if (flow != nullptr &&
      hundun::flow::test::IdealGasClosureTestAccess::
          allocation_observation_active(*flow))
    allocation_calls.fetch_add(1U, std::memory_order_relaxed);
}

void *allocate(std::size_t size) {
  record();
  if (void *result = std::malloc(size == 0U ? 1U : size))
    return result;
  throw std::bad_alloc();
}

void *allocate_aligned(std::size_t size, std::size_t alignment) {
  record();
  void *result = nullptr;
  if (posix_memalign(&result, alignment, size == 0U ? 1U : size) == 0)
    return result;
  throw std::bad_alloc();
}

class Observation final {
public:
  explicit Observation(const hundun::flow::FixedStepIdealGasFlow &flow) {
    allocation_calls.store(0U, std::memory_order_relaxed);
    observed_flow.store(&flow, std::memory_order_relaxed);
  }
  ~Observation() noexcept {
    observed_flow.store(nullptr, std::memory_order_relaxed);
  }
  Observation(const Observation &) = delete;
  Observation &operator=(const Observation &) = delete;
  std::uint64_t count() const noexcept {
    return allocation_calls.load(std::memory_order_relaxed);
  }
};

} // namespace task21_allocation_oracle

void *operator new(std::size_t size) {
  return task21_allocation_oracle::allocate(size);
}
void *operator new[](std::size_t size) {
  return task21_allocation_oracle::allocate(size);
}
void operator delete(void *pointer) noexcept { std::free(pointer); }
void operator delete[](void *pointer) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete[](void *pointer, std::size_t) noexcept {
  std::free(pointer);
}
void *operator new(std::size_t size, std::align_val_t alignment) {
  return task21_allocation_oracle::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void *operator new[](std::size_t size, std::align_val_t alignment) {
  return task21_allocation_oracle::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void operator delete(void *pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void *pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void *pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void *pointer, std::size_t,
                       std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

#ifndef HUNDUN_TASK21_DIAGNOSTIC_MATRIX
#define HUNDUN_TASK21_DIAGNOSTIC_MATRIX(source, mpi, state, report, flow,     \
                                         stencil)                             \
  static_cast<void>(source)
#endif
#ifndef HUNDUN_TASK21_DIAGNOSTIC_EARLY_RETURN
#define HUNDUN_TASK21_DIAGNOSTIC_EARLY_RETURN() false
#endif
#ifndef HUNDUN_TASK21_DIAGNOSTIC_STATUS
#define HUNDUN_TASK21_DIAGNOSTIC_STATUS(source, mpi, status, classification,  \
                                        code, rank)                           \
  static_cast<void>(source)
#endif

hundun::runtime::Int3 grid(int ranks) { return {ranks, 1, 1}; }

void require_complete_final_gates(
    const hundun::flow::IdealGasStepAttemptReport &report) {
  const auto &material = report.flow();
  const auto &flow = material.flow();
  const auto &closure = report.closure_report();
  HUNDUN_CHECK(flow.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(flow.pressure_corrector_count == 2U);
  HUNDUN_CHECK(material.final_continuity_residual_available());
  HUNDUN_CHECK(flow.final_continuity_normalized_l2 <= 1.0e-10);
  HUNDUN_CHECK(material.final_pressure_residual_available());
  HUNDUN_CHECK(material.final_pressure_normalized_residual() <= 1.0);
  HUNDUN_CHECK(std::all_of(
      material.final_momentum_residual_availability().begin(),
      material.final_momentum_residual_availability().end(),
      [](std::uint8_t available) { return available == 1U; }));
  for (const double residual : flow.final_momentum_normalized_l2)
    HUNDUN_CHECK(residual <= 1.0e-9);
  for (const double residual : flow.final_transport_normalized_l2)
    HUNDUN_CHECK(residual <= 1.0e-9);
  HUNDUN_CHECK(material.mass_conservation_available());
  HUNDUN_CHECK(flow.final_mass_relative_conservation_defect <= 5.0e-11);
  HUNDUN_CHECK(std::all_of(
      material.momentum_conservation_availability().begin(),
      material.momentum_conservation_availability().end(),
      [](std::uint8_t available) { return available == 1U; }));
  for (const double defect : flow.final_momentum_relative_conservation_defect)
    HUNDUN_CHECK(defect <= 5.0e-11);
  for (const double defect :
       flow.final_transport_relative_conservation_defect)
    HUNDUN_CHECK(defect <= 5.0e-11);
  HUNDUN_CHECK(closure.final_metrics_available());
  HUNDUN_CHECK(closure.rho_remap_normalized_l2() <= 1.0e-10);
  HUNDUN_CHECK(closure.rho_h_remap_normalized_l2() <= 1.0e-9);
  HUNDUN_CHECK(closure.rho_remap_relative_conservation_defect() <= 5.0e-11);
  HUNDUN_CHECK(closure.rho_h_remap_relative_conservation_defect() <=
               5.0e-11);
  HUNDUN_CHECK(closure.enthalpy_temperature_max_relative_error() <= 1.0e-12);
  HUNDUN_CHECK(closure.eos_max_relative_error() <= 1.0e-12);
}

void require_halo_trace(
    const hundun::flow::FixedStepIdealGasFlow &flow,
    const hundun::flow::FlowFieldIds &fields,
    std::initializer_list<hundun::flow::IdealGasClosureStage> stages) {
  const auto trace =
      hundun::flow::test::IdealGasClosureTestAccess::halo_trace(flow);
  HUNDUN_CHECK(trace.size() == stages.size());
  std::size_t index = 0U;
  for (const auto stage : stages) {
    HUNDUN_CHECK(trace[index].stage == stage);
    HUNDUN_CHECK(trace[index].density == fields.density);
    HUNDUN_CHECK(trace[index].enthalpy_density ==
                 fields.transported_cell_fields.front());
    ++index;
  }
}

double global_extensive_total(const std::vector<double> &values,
                              const hundun::mesh::MeshGeometry &geometry,
                              const hundun::runtime::MpiContext &mpi) {
  double local = 0.0;
  for (std::size_t cell = 0; cell < values.size(); ++cell)
    local += values[cell] * geometry.cell_volume_m3(cell);
  double total = 0.0;
  HUNDUN_CHECK(MPI_Allreduce(&local, &total, 1, MPI_DOUBLE, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  return total;
}

struct DecompositionReference final {
  std::uint64_t components{};
  std::vector<double> cell_fields;
  std::vector<double> face_mass_flux;
};

DecompositionReference make_reference(
    const hundun::mesh::MeshTopology &topology,
    const hundun::flow::FlowLayerValues &values,
    const hundun::runtime::MpiContext &mpi) {
  const std::size_t cells =
      static_cast<std::size_t>(topology.global_cell_count());
  const std::size_t faces =
      static_cast<std::size_t>(topology.global_face_count());
  const std::size_t components = 5U + values.transported_cell_fields.size();
  DecompositionReference result{components,
                                std::vector<double>(cells * components),
                                std::vector<double>(faces)};
  std::vector<double> cell_cover(cells);
  std::vector<double> face_cover(faces);
  for (hundun::mesh::LocalCellId cell_id = 0;
       cell_id < topology.owned_cell_count(); ++cell_id) {
    const auto global =
        static_cast<std::size_t>(topology.global_cell_id(cell_id));
    cell_cover[global] = 1.0;
    result.cell_fields[global * components] = values.density[cell_id];
    for (std::size_t direction = 0; direction < 3U; ++direction)
      result.cell_fields[global * components + 1U + direction] =
          values.velocity[cell_id * 3U + direction];
    result.cell_fields[global * components + 4U] =
        values.mechanical_pressure[cell_id];
    for (std::size_t field = 0;
         field < values.transported_cell_fields.size(); ++field)
      result.cell_fields[global * components + 5U + field] =
          values.transported_cell_fields[field][cell_id];
  }
  for (hundun::mesh::LocalFaceId face_id = 0;
       face_id < topology.local_face_count(); ++face_id) {
    if (topology.face_ownership(face_id) !=
        hundun::mesh::EntityOwnership::owned)
      continue;
    const auto global =
        static_cast<std::size_t>(topology.global_face_id(face_id));
    face_cover[global] = 1.0;
    result.face_mass_flux[global] = values.face_mass_flux[face_id];
  }
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, cell_cover.data(),
                             static_cast<int>(cell_cover.size()), MPI_DOUBLE,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, result.cell_fields.data(),
                             static_cast<int>(result.cell_fields.size()),
                             MPI_DOUBLE, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, face_cover.data(),
                             static_cast<int>(face_cover.size()), MPI_DOUBLE,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, result.face_mass_flux.data(),
                             static_cast<int>(result.face_mass_flux.size()),
                             MPI_DOUBLE, MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(std::all_of(cell_cover.begin(), cell_cover.end(),
                           [](double value) { return value == 1.0; }));
  HUNDUN_CHECK(std::all_of(face_cover.begin(), face_cover.end(),
                           [](double value) { return value == 1.0; }));
  return result;
}

void write_reference(const std::filesystem::path &path,
                     const DecompositionReference &reference) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output)
    throw std::runtime_error("cannot create ideal-gas reference artifact");
  constexpr std::uint64_t magic = UINT64_C(0x48554e4449474153);
  const std::uint64_t cells = reference.cell_fields.size();
  const std::uint64_t faces = reference.face_mass_flux.size();
  output.write(reinterpret_cast<const char *>(&magic), sizeof(magic));
  output.write(reinterpret_cast<const char *>(&reference.components),
               sizeof(reference.components));
  output.write(reinterpret_cast<const char *>(&cells), sizeof(cells));
  output.write(reinterpret_cast<const char *>(&faces), sizeof(faces));
  output.write(reinterpret_cast<const char *>(reference.cell_fields.data()),
               static_cast<std::streamsize>(cells * sizeof(double)));
  output.write(
      reinterpret_cast<const char *>(reference.face_mass_flux.data()),
      static_cast<std::streamsize>(faces * sizeof(double)));
  if (!output)
    throw std::runtime_error("cannot write ideal-gas reference artifact");
}

DecompositionReference read_reference(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input)
    throw std::runtime_error("ideal-gas reference artifact is missing");
  std::uint64_t magic{}, components{}, cells{}, faces{};
  input.read(reinterpret_cast<char *>(&magic), sizeof(magic));
  input.read(reinterpret_cast<char *>(&components), sizeof(components));
  input.read(reinterpret_cast<char *>(&cells), sizeof(cells));
  input.read(reinterpret_cast<char *>(&faces), sizeof(faces));
  if (magic != UINT64_C(0x48554e4449474153) || components == 0U)
    throw std::runtime_error("ideal-gas reference header differs");
  DecompositionReference result{components, std::vector<double>(cells),
                                std::vector<double>(faces)};
  input.read(reinterpret_cast<char *>(result.cell_fields.data()),
             static_cast<std::streamsize>(cells * sizeof(double)));
  input.read(reinterpret_cast<char *>(result.face_mass_flux.data()),
             static_cast<std::streamsize>(faces * sizeof(double)));
  if (!input || input.peek() != std::ifstream::traits_type::eof())
    throw std::runtime_error("ideal-gas reference artifact is malformed");
  return result;
}

void compare_reference(const DecompositionReference &reference,
                       const DecompositionReference &candidate) {
  HUNDUN_CHECK(reference.components == candidate.components);
  HUNDUN_CHECK(reference.cell_fields.size() == candidate.cell_fields.size());
  HUNDUN_CHECK(reference.face_mass_flux.size() ==
               candidate.face_mass_flux.size());
  const std::size_t cells =
      reference.cell_fields.size() / reference.components;
  for (std::size_t component = 0; component < reference.components;
       ++component) {
    double difference{};
    double infinity_norm{};
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const double expected =
          reference.cell_fields[cell * reference.components + component];
      const double observed =
          candidate.cell_fields[cell * reference.components + component];
      difference = std::max(difference, std::abs(observed - expected));
      infinity_norm = std::max(infinity_norm, std::abs(expected));
    }
    HUNDUN_CHECK(difference <= 5.0e-12 * std::max(1.0, infinity_norm));
  }
  double difference{};
  double infinity_norm{};
  for (std::size_t face = 0; face < reference.face_mass_flux.size(); ++face) {
    difference = std::max(
        difference,
        std::abs(candidate.face_mass_flux[face] -
                 reference.face_mass_flux[face]));
    infinity_norm =
        std::max(infinity_norm, std::abs(reference.face_mass_flux[face]));
  }
  HUNDUN_CHECK(difference <= 5.0e-12 * std::max(1.0, infinity_norm));
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::ideal_gas;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.physics.cp_J_per_kg_K = 1000.0;
  config.physics.gas_constant_J_per_kg_K = 287.05;
  config.physics.thermodynamic_pressure_pa = 101325.0;
  constexpr std::array names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t patch = 0; patch < names.size(); ++patch) {
    config.boundaries[patch].patch = names[patch];
    config.boundaries[patch].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::config::FlowCaseConfig open_case(bool with_scalar) {
  auto config = periodic_case();
  if (with_scalar)
    config.scalars.push_back({"phi", 0.0});
  config.boundaries[0].type = hundun::config::BoundaryType::velocity_inlet;
  config.boundaries[0].velocity_m_per_s = {1.0, 0.0, 0.0};
  config.boundaries[0].thermal_authority =
      hundun::config::InletThermalAuthority::temperature;
  config.boundaries[0].temperature_K = 300.0;
  config.boundaries[0].enthalpy_J_per_kg = 300000.0;
  config.boundaries[0].density_kg_per_m3 = 101325.0 / (287.05 * 300.0);
  config.boundaries[0].scalar_values =
      with_scalar
          ? std::vector<hundun::config::InletScalarValue>{{"phi", 0.25}}
          : std::vector<hundun::config::InletScalarValue>{};
  config.boundaries[1].type = hundun::config::BoundaryType::pressure_outlet;
  config.boundaries[1].pressure_perturbation_pa = 0.0;
  for (std::size_t patch = 2U; patch < config.boundaries.size(); ++patch)
    config.boundaries[patch].type = hundun::config::BoundaryType::symmetry;
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, const char *unit,
                                      std::uint32_t components,
                                      bool conservative) {
  return {name,
          unit,
          "task21",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          components,
          2,
          conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face(const char *name, const char *unit,
                                      std::uint32_t components) {
  return {name,
          unit,
          "task21",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

void run(const hundun::runtime::MpiContext &mpi, bool open_domain,
         bool exhaust_generation, bool full_case,
         std::string_view reference_mode = {},
         const std::filesystem::path &reference_path = {}) {
  const hundun::runtime::Int3 extent{
      open_domain && full_case ? 16 : 8, 4, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {!open_domain, !open_domain, !open_domain},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology topology(decomposition);
  hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping(
          {0.0, 0.0, 0.0},
          open_domain ? hundun::runtime::Real3{1.0, 0.25, 0.25}
                      : hundun::runtime::Real3{1.0, 1.0, 1.0}));
  auto boundaries = hundun::boundary::BoundaryRegistry::create(
      open_domain ? open_case(full_case) : periodic_case(), topology);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("rho", "kg/m3", 1U, true));
  fields.velocity = registry.declare_field(cell("velocity", "m/s", 3U, false));
  fields.mechanical_pressure =
      registry.declare_field(cell("pi", "Pa", 1U, false));
  fields.face_velocity =
      registry.declare_field(face("face_velocity", "m/s", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  const auto rho_h = registry.declare_field(cell("rho_h", "J/m3", 1U, true));
  fields.transported_cell_fields = {rho_h};
  hundun::runtime::FieldId rho_phi{};
  if (open_domain && full_case) {
    rho_phi = registry.declare_field(cell("rho_phi", "kg/m3", 1U, true));
    fields.transported_cell_fields.push_back(rho_phi);
  }
  registry.freeze();

  const double dt = open_domain && full_case ? 6.25e-3 : 1.0e-3;
  constexpr double cp = 1000.0;
  constexpr double gas_constant = 287.05;
  constexpr double initial_temperature = 300.0;
  constexpr double pressure = 101325.0;
  const double density = pressure / (gas_constant * initial_temperature);
  auto state = hundun::flow::FlowState::create(
      registry, {decomposition.local_extent(), topology.local_face_count()},
      fields,
      {0U, 0.0, dt, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(topology.owned_cell_count(), density);
  initial.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  if (open_domain)
    for (std::size_t cell_index = 0; cell_index < topology.owned_cell_count();
         ++cell_index)
      initial.velocity[cell_index * 3U] = 1.0;
  initial.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
  if (open_domain) {
    for (hundun::mesh::LocalFaceId face = 0;
         face < topology.local_face_count(); ++face) {
      const auto area = geometry.face_area_vector_m2(
          face, hundun::mesh::FaceSide::owner);
      initial.face_velocity[face * 3U] = 1.0;
      initial.face_mass_flux[face] = density * area.x;
    }
  }
  initial.transported_cell_fields = {std::vector<double>(
      topology.owned_cell_count(), density * cp * initial_temperature)};
  if (open_domain && full_case)
    initial.transported_cell_fields.push_back(std::vector<double>(
        topology.owned_cell_count(), density * 0.25));
  state.seed_accepted_layers(initial, initial);
  const double initial_enthalpy_integral = global_extensive_total(
      initial.transported_cell_fields.front(), geometry, mpi);

  auto closure = hundun::flow::IdealGasClosure::create(
      topology, geometry, boundaries, mpi, registry, fields, state,
      {rho_h, cp, gas_constant, pressure});

  // The closure evaluates without Halo exchange here, so this directly
  // verifies the post-gate publication boundary: only owned interior rho and
  // rho_h may change. Ghosts and unrelated fields must retain their exact
  // bit patterns.
  state.begin_attempt();
  {
    auto &trial = state.trial_layer();
    auto rho = trial.view<double>(fields.density);
    auto q = trial.view<double>(rho_h);
    auto pi = trial.view<double>(fields.mechanical_pressure);
    const auto local = decomposition.local_extent();
    for (int k = -2; k < local.z + 2; ++k)
      for (int j = -2; j < local.y + 2; ++j)
        for (int i = -2; i < local.x + 2; ++i) {
          const bool interior = i >= 0 && i < local.x && j >= 0 &&
                                j < local.y && k >= 0 && k < local.z;
          if (!interior) {
            rho(i, j, k, 0) = -101.0;
            q(i, j, k, 0) = -202.0;
          }
          pi(i, j, k, 0) = interior ? 303.0 : -303.0;
        }
    q(0, 0, 0, 0) *= 1.01;
    const auto rho_before = state.snapshot(hundun::flow::FlowLayer::trial);
    hundun::flow::test::IdealGasClosureTestAccess::begin_attempt(
        closure, state, 9001U);
    const auto publication =
        hundun::flow::test::IdealGasClosureTestAccess::evaluate(
            closure, state, hundun::flow::IdealGasClosureStage::predictor);
    HUNDUN_CHECK(publication.disposition() ==
                 hundun::flow::IdealGasClosureDisposition::closed);
    const double candidate_pressure = publication.candidate_pressure_pa();
    bool rho_interior_changed = false;
    bool rho_h_interior_changed = false;
    for (int k = 0; k < local.z; ++k)
      for (int j = 0; j < local.y; ++j)
        for (int i = 0; i < local.x; ++i) {
          const auto offset =
              (static_cast<std::size_t>(k) *
                   static_cast<std::size_t>(local.y) +
               static_cast<std::size_t>(j)) *
                  static_cast<std::size_t>(local.x) +
              static_cast<std::size_t>(i);
          const double rho_tilde = rho_before.density[offset];
          const double q_tilde =
              rho_before.transported_cell_fields.front()[offset];
          const double enthalpy = q_tilde / rho_tilde;
          const double temperature = enthalpy / cp;
          const double expected_rho = static_cast<double>(
              static_cast<long double>(candidate_pressure) /
              (static_cast<long double>(gas_constant) * temperature));
          const double expected_rho_h = static_cast<double>(
              static_cast<long double>(expected_rho) * enthalpy);
          HUNDUN_CHECK(hundun::test::fp64_bits(rho(i, j, k, 0)) ==
                       hundun::test::fp64_bits(expected_rho));
          HUNDUN_CHECK(hundun::test::fp64_bits(q(i, j, k, 0)) ==
                       hundun::test::fp64_bits(expected_rho_h));
          rho_interior_changed =
              rho_interior_changed ||
              hundun::test::fp64_bits(rho(i, j, k, 0)) !=
                  hundun::test::fp64_bits(rho_tilde);
          rho_h_interior_changed =
              rho_h_interior_changed ||
              hundun::test::fp64_bits(q(i, j, k, 0)) !=
                  hundun::test::fp64_bits(q_tilde);
          HUNDUN_CHECK(
              hundun::test::fp64_bits(pi(i, j, k, 0)) ==
              hundun::test::fp64_bits(
                  rho_before.mechanical_pressure[offset]));
        }
    HUNDUN_CHECK(rho_interior_changed);
    HUNDUN_CHECK(rho_h_interior_changed);
    for (int k = -2; k < local.z + 2; ++k)
      for (int j = -2; j < local.y + 2; ++j)
        for (int i = -2; i < local.x + 2; ++i) {
          const bool interior = i >= 0 && i < local.x && j >= 0 &&
                                j < local.y && k >= 0 && k < local.z;
          if (interior)
            continue;
          HUNDUN_CHECK(hundun::test::fp64_bits(rho(i, j, k, 0)) ==
                       hundun::test::fp64_bits(-101.0));
          HUNDUN_CHECK(hundun::test::fp64_bits(q(i, j, k, 0)) ==
                       hundun::test::fp64_bits(-202.0));
          HUNDUN_CHECK(hundun::test::fp64_bits(pi(i, j, k, 0)) ==
                       hundun::test::fp64_bits(-303.0));
        }
    hundun::flow::test::IdealGasClosureTestAccess::rollback(closure);
  }
  state.rollback_attempt();
  auto halo = hundun::runtime::HaloExchange::create(
      decomposition, hundun::runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 2));
  hundun::execution::CpuReferenceContext execution;
  hundun::linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  hundun::linear::ConjugateGradientSolver pressure_solver(execution, mpi);
  hundun::linear::JacobiPreconditioner mx(execution), my(execution),
      mz(execution), pressure_preconditioner(execution);
  hundun::flow::MaterialDensityTransportSpec material_spec;
  material_spec.enthalpy_density = rho_h;
  material_spec.enthalpy_diffusivity_kg_per_m_s = 0.0;
  if (open_domain && full_case) {
    material_spec.scalar_densities = {rho_phi};
    material_spec.scalar_diffusivities_kg_per_m_s = {0.0};
  }
  hundun::mesh::MeshTopology alternate_topology(decomposition);
  hundun::mesh::MeshGeometry alternate_geometry(
      alternate_topology,
      hundun::mesh::UniformBoxMapping(
          {0.0, 0.0, 0.0},
          open_domain ? hundun::runtime::Real3{1.0, 0.25, 0.25}
                      : hundun::runtime::Real3{1.0, 1.0, 1.0}));
  auto alternate_boundaries = hundun::boundary::BoundaryRegistry::create(
      open_domain ? open_case(full_case) : periodic_case(), topology);
  const auto expect_collaborator_rejected =
      [&](const hundun::mesh::MeshTopology &candidate_topology,
          const hundun::mesh::MeshGeometry &candidate_geometry,
          const hundun::boundary::BoundaryRegistry &candidate_boundaries) {
        bool rejected = false;
        try {
          static_cast<void>(hundun::flow::FixedStepIdealGasFlow::create(
              decomposition, candidate_topology, candidate_geometry,
              candidate_boundaries, mpi, execution, halo, momentum_solver,
              {&mx, &my, &mz}, pressure_solver, pressure_preconditioner,
              registry, fields, material_spec, std::move(closure)));
        } catch (const hundun::runtime::Error &) {
          rejected = true;
        }
        HUNDUN_CHECK(rejected);
      };
  expect_collaborator_rejected(topology, alternate_geometry, boundaries);
  expect_collaborator_rejected(alternate_topology, geometry, boundaries);
  expect_collaborator_rejected(topology, geometry, alternate_boundaries);
  for (const int fault_rank : {0, mpi.size() - 1}) {
    for (const bool nested_material_factory : {true, false}) {
      auto fault_closure = hundun::flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, registry, fields, state,
          {rho_h, cp, gas_constant, pressure});
      if (nested_material_factory)
        hundun::flow::test::IdealGasClosureTestAccess::
            set_material_factory_create_fault(fault_closure, fault_rank);
      else
        hundun::flow::test::IdealGasClosureTestAccess::set_facade_create_fault(
            fault_closure, fault_rank);
      bool rejected = false;
      try {
        static_cast<void>(hundun::flow::FixedStepIdealGasFlow::create(
            decomposition, topology, geometry, boundaries, mpi, execution,
            halo, momentum_solver, {&mx, &my, &mz}, pressure_solver,
            pressure_preconditioner, registry, fields, material_spec,
            std::move(fault_closure)));
      } catch (const hundun::runtime::Error &error) {
        rejected =
            std::string(error.what()) ==
            "ideal-gas flow facade construction failed on rank " +
                std::to_string(fault_rank);
      }
      HUNDUN_CHECK(rejected);

      auto successor_state = hundun::flow::FlowState::create(
          registry,
          {decomposition.local_extent(), topology.local_face_count()}, fields,
          {0U, 0.0, dt, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      successor_state.seed_accepted_layers(initial, initial);
      auto successor_closure = hundun::flow::IdealGasClosure::create(
          topology, geometry, boundaries, mpi, registry, fields,
          successor_state,
          {rho_h, cp, gas_constant, pressure});
      auto successor = hundun::flow::FixedStepIdealGasFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&mx, &my, &mz}, pressure_solver,
          pressure_preconditioner, registry, fields, material_spec,
          std::move(successor_closure));
      const auto successor_report = successor.attempt(
          successor_state, 0.0,
          hundun::flow::make_momentum_time_stencil(
              hundun::flow::MomentumTimeOrder::backward_euler, dt, 0.0),
          {}, {});
      HUNDUN_CHECK(successor_report.flow().flow().disposition ==
                   hundun::flow::StepAttemptDisposition::committed);
    }
    if (mpi.size() == 1)
      break;
  }
  auto flow = hundun::flow::FixedStepIdealGasFlow::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner, registry, fields, material_spec,
      std::move(closure));
  HUNDUN_CHECK(
      hundun::test::ideal_gas_state_equality_oracle_is_mutation_sensitive());
  const auto stencil = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, dt, 0.0);
  for (const int allocation_rank : {0, mpi.size() - 1}) {
    auto allocation_state = hundun::flow::FlowState::create(
        registry, {decomposition.local_extent(), topology.local_face_count()},
        fields,
        {0U, 0.0, dt, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    allocation_state.seed_accepted_layers(initial, initial);
    auto allocation_closure = hundun::flow::IdealGasClosure::create(
        topology, geometry, boundaries, mpi, registry, fields,
        allocation_state, {rho_h, cp, gas_constant, pressure});
    auto allocation_flow = hundun::flow::FixedStepIdealGasFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_preconditioner, registry, fields, material_spec,
        std::move(allocation_closure));
    const hundun::test::IdealGasStateSnapshot before_allocation_probe{
        allocation_state.snapshot(hundun::flow::FlowLayer::history),
        allocation_state.snapshot(hundun::flow::FlowLayer::committed),
        allocation_state.snapshot(hundun::flow::FlowLayer::trial),
        allocation_state.metadata(), allocation_flow.closure_state()};
    hundun::flow::test::IdealGasClosureTestAccess::set_controlled_allocation(
        allocation_flow, allocation_rank);
    task21_allocation_oracle::Observation mutation_observation(allocation_flow);
    const auto allocation_probe =
        allocation_flow.attempt(allocation_state, 0.0, stencil, {}, {});
    HUNDUN_CHECK((mutation_observation.count() > 0U) ==
                 (mpi.rank() == allocation_rank));
    HUNDUN_CHECK(allocation_probe.flow().flow().disposition !=
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(allocation_probe.flow().flow().reason ==
                 hundun::flow::StepFailureReason::density_closure_failure);
    HUNDUN_CHECK(allocation_probe.flow().flow().lowest_failing_rank ==
                 allocation_rank);
    HUNDUN_CHECK(allocation_probe.closure_report_available());
    HUNDUN_CHECK(allocation_probe.closure_report().reason() ==
                 hundun::flow::IdealGasClosureFailureReason::invalid_input);
    const hundun::test::IdealGasStateSnapshot after_allocation_probe{
        allocation_state.snapshot(hundun::flow::FlowLayer::history),
        allocation_state.snapshot(hundun::flow::FlowLayer::committed),
        allocation_state.snapshot(hundun::flow::FlowLayer::trial),
        allocation_state.metadata(), allocation_flow.closure_state()};
    HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
        before_allocation_probe, after_allocation_probe));
    task21_allocation_oracle::Observation successor_observation(
        allocation_flow);
    const auto allocation_successor = allocation_flow.attempt(
        allocation_state, 0.0, stencil, {}, {});
    HUNDUN_CHECK(successor_observation.count() == 0U);
    HUNDUN_CHECK(allocation_successor.flow().flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    if (mpi.size() == 1)
      break;
  }
  hundun::flow::test::IdealGasClosureTestAccess::set_uniform_enthalpy_rate(
      flow, open_domain ? 0.0 : cp * 5.0 / dt);
  task21_allocation_oracle::Observation clean_observation(flow);
  const auto report = flow.attempt(state, 0.0, stencil, {}, {});
  HUNDUN_CHECK(clean_observation.count() == 0U);
  HUNDUN_CHECK(report.flow().flow().disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(report.flow().flow().reason ==
               hundun::flow::StepFailureReason::none);
  HUNDUN_CHECK(report.flow().flow().pressure_corrector_count == 2U);
  HUNDUN_CHECK(report.closure_report_available());
  const auto &closure_report = report.closure_report();
  HUNDUN_CHECK(closure_report.disposition() ==
               hundun::flow::IdealGasClosureDisposition::closed);
  HUNDUN_CHECK(closure_report.stage() ==
               hundun::flow::IdealGasClosureStage::final);
  HUNDUN_CHECK(closure_report.evaluation_count() == 3U);
  HUNDUN_CHECK(closure_report.collective_count() == 14U);
  HUNDUN_CHECK(closure_report.candidate_pressure_available());
  if (open_domain) {
    const double candidate = closure_report.candidate_pressure_pa();
    HUNDUN_CHECK(std::memcmp(&candidate, &pressure, sizeof(double)) == 0);
  } else
    HUNDUN_CHECK_NEAR(closure_report.candidate_pressure_pa(),
                      pressure * 305.0 / 300.0, 1.0e-12 * pressure);
  HUNDUN_CHECK(closure_report.final_metrics_available());
  require_halo_trace(flow, fields,
                     {hundun::flow::IdealGasClosureStage::predictor,
                      hundun::flow::IdealGasClosureStage::provisional,
                      hundun::flow::IdealGasClosureStage::final});
  HUNDUN_CHECK(closure_report.rho_remap_normalized_l2() <= 1.0e-10);
  HUNDUN_CHECK(closure_report.rho_h_remap_normalized_l2() <= 1.0e-9);
  HUNDUN_CHECK(closure_report.eos_max_relative_error() <= 1.0e-12);
  require_complete_final_gates(report);
  HUNDUN_CHECK(
      hundun::flow::test::IdealGasClosureTestAccess::report_authenticated(
          report));
  HUNDUN_CHECK(hundun::flow::test::IdealGasClosureTestAccess::
                   post_eos_evidence_authenticated(report.flow()));
  for (std::uint8_t mutation = 0U;
       mutation < static_cast<std::uint8_t>(
                      hundun::flow::test::IdealGasPostEvidenceMutation::count);
       ++mutation) {
    if (!hundun::flow::test::IdealGasClosureTestAccess::
            post_evidence_mutation_rejected(
                report,
                static_cast<hundun::flow::test::IdealGasPostEvidenceMutation>(
                    mutation)))
      throw std::runtime_error("post-EOS success mutation was accepted: " +
                               std::to_string(mutation));
  }
  HUNDUN_CHECK(state.metadata().step == 1U);
  auto closure_source = flow.closure_diagnostic_source(state, report);
  HUNDUN_CHECK(closure_source.fingerprint_field_count() == 3U);
  HUNDUN_CHECK(closure_source.fingerprint_field_id(0U) ==
               std::string_view("p0"));
  HUNDUN_CHECK(closure_source.fingerprint_field_id(1U) ==
               std::string_view("rho"));
  HUNDUN_CHECK(closure_source.fingerprint_field_id(2U) ==
               std::string_view("rho_h"));
  HUNDUN_CHECK(closure_source.fingerprint_field_item_count(0U) ==
               (mpi.rank() == 0 ? 1U : 0U));
  HUNDUN_CHECK(closure_source.sample_field_count() == 5U);
  HUNDUN_CHECK(closure_source.sample_field_id(4U) ==
               std::string_view("temperature"));
  const double first_temperature = open_domain ? 300.0 : 305.0;
  HUNDUN_CHECK_NEAR(closure_source.sample_field_value(4U, 0U),
                    first_temperature, 5.0e-12 * first_temperature);
  HUNDUN_TASK21_DIAGNOSTIC_MATRIX(closure_source, mpi, state, report, flow,
                                  stencil);
  if (HUNDUN_TASK21_DIAGNOSTIC_EARLY_RETURN())
    return;
  auto material_source = flow.flow_diagnostic_source(state, report);
  HUNDUN_CHECK(material_source.report().attempt_identity() ==
               report.attempt_identity());
  const auto closure_state = flow.closure_state();
  HUNDUN_CHECK(closure_state.revision == 1U);
  HUNDUN_CHECK_NEAR(closure_state.thermodynamic_pressure_pa,
                    open_domain ? pressure : pressure * 305.0 / 300.0,
                    1.0e-12 * pressure);
  HUNDUN_CHECK(closure_state.target_mass_kg.has_value() != open_domain);
  const auto committed = state.snapshot(hundun::flow::FlowLayer::committed);
  HUNDUN_CHECK(std::all_of(
      committed.density.begin(), committed.density.end(), [&](double value) {
        return std::abs(value - density) <= 5.0e-12 * density;
      }));
  for (std::size_t cell_index = 0; cell_index < committed.density.size();
       ++cell_index) {
    const double h = committed.transported_cell_fields.front()[cell_index] /
                     committed.density[cell_index];
    HUNDUN_CHECK_NEAR(h / cp, first_temperature, 5.0e-12 * first_temperature);
  }
  if (!open_domain) {
    const double expected =
        *flow.closure_state().target_mass_kg * cp * 5.0;
    const double observed = global_extensive_total(
                                committed.transported_cell_fields.front(),
                                geometry, mpi) -
                            initial_enthalpy_integral;
    HUNDUN_CHECK(std::abs(observed - expected) /
                     std::max(std::abs(expected),
                              std::numeric_limits<double>::min()) <=
                 5.0e-11);
    // The mutation-sensitive source oracle must reject a balance that omits
    // the configured source altogether.
    HUNDUN_CHECK(std::abs(observed) /
                     std::max(std::abs(expected),
                              std::numeric_limits<double>::min()) >
                 5.0e-11);
  }

  const auto bdf2 = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::bdf2, dt, dt);
  const auto second = flow.attempt(state, 0.0, bdf2, {}, {});
  HUNDUN_CHECK(second.flow().flow().disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  require_complete_final_gates(second);
  HUNDUN_CHECK(second.closure_report().evaluation_count() == 3U);
  HUNDUN_CHECK(second.closure_report().collective_count() == 14U);
  require_halo_trace(flow, fields,
                     {hundun::flow::IdealGasClosureStage::predictor,
                      hundun::flow::IdealGasClosureStage::provisional,
                      hundun::flow::IdealGasClosureStage::final});
  HUNDUN_CHECK(second.closure_report().candidate_pressure_available());
  if (open_domain) {
    const double candidate = second.closure_report().candidate_pressure_pa();
    HUNDUN_CHECK(std::memcmp(&candidate, &pressure, sizeof(double)) == 0);
  } else {
    HUNDUN_CHECK_NEAR(second.closure_report().candidate_pressure_pa(),
                      pressure * 310.0 / 300.0, 1.0e-12 * pressure);
  }
  HUNDUN_CHECK(flow.closure_state().revision == 2U);
  HUNDUN_CHECK_NEAR(flow.closure_state().thermodynamic_pressure_pa,
                    open_domain ? pressure : pressure * 310.0 / 300.0,
                    1.0e-12 * pressure);
  if (!open_domain) {
    const auto second_committed =
        state.snapshot(hundun::flow::FlowLayer::committed);
    const double expected =
        *flow.closure_state().target_mass_kg * cp * 10.0;
    const double observed = global_extensive_total(
                                second_committed
                                    .transported_cell_fields.front(),
                                geometry, mpi) -
                            initial_enthalpy_integral;
    HUNDUN_CHECK(std::abs(observed - expected) /
                     std::max(std::abs(expected),
                              std::numeric_limits<double>::min()) <=
                 5.0e-11);
    HUNDUN_CHECK(std::abs(observed) /
                     std::max(std::abs(expected),
                              std::numeric_limits<double>::min()) >
                 5.0e-11);
  }
  if (open_domain && full_case) {
    for (int step = 2; step < 4; ++step) {
      const auto continued = flow.attempt(state, 0.0, bdf2, {}, {});
      HUNDUN_CHECK(continued.flow().flow().disposition ==
                   hundun::flow::StepAttemptDisposition::committed);
      HUNDUN_CHECK(continued.flow().flow().pressure_corrector_count == 2U);
      require_complete_final_gates(continued);
      HUNDUN_CHECK(continued.closure_report().evaluation_count() == 3U);
      HUNDUN_CHECK(continued.closure_report().collective_count() == 14U);
      HUNDUN_CHECK(continued.flow().shared_face_mass_flux_field() ==
                   fields.face_mass_flux);
      HUNDUN_CHECK(continued.flow().flux_provenance() ==
                   hundun::flow::MaterialFluxProvenance::final_corrected);
    }
    HUNDUN_CHECK(state.metadata().step == 4U);
    const auto fixed_closure = flow.closure_state();
    HUNDUN_CHECK(std::memcmp(&fixed_closure.thermodynamic_pressure_pa,
                             &pressure, sizeof(double)) == 0);
    const auto final = state.snapshot(hundun::flow::FlowLayer::committed);
    HUNDUN_CHECK(final.transported_cell_fields.size() == 2U);
    for (std::size_t cell_index = 0; cell_index < final.density.size();
         ++cell_index) {
      HUNDUN_CHECK_NEAR(final.density[cell_index], density,
                        5.0e-12 * density);
      HUNDUN_CHECK_NEAR(
          final.transported_cell_fields[1][cell_index] /
              final.density[cell_index],
          0.25, 5.0e-12);
    }
    // Integrate the accepted shared face-mass-flux field independently of
    // the transport residual implementation.  Each extensive quantity uses
    // the same accepted flux and the uniform analytic face state.
    const double area = 0.25 * 0.25;
    const std::array multipliers{1.0, cp * initial_temperature, 1.0, 0.25};
    std::array<double, 16> local_boundary{};
    for (std::size_t patch = 0; patch < 2U; ++patch) {
      for (const auto face_id :
           topology.patch(static_cast<std::uint32_t>(patch)).local_faces()) {
        if (topology.face_ownership(face_id) !=
            hundun::mesh::EntityOwnership::owned)
          continue;
        const double flux = final.face_mass_flux[face_id];
        for (std::size_t quantity = 0; quantity < multipliers.size();
             ++quantity) {
          const double contribution = flux * multipliers[quantity];
          local_boundary[quantity * 2U] += contribution;
          local_boundary[quantity * 2U + 1U] += std::abs(contribution);
          local_boundary[8U + patch * 4U + quantity] += contribution;
        }
      }
    }
    std::array<double, 16> boundary{};
    HUNDUN_CHECK(MPI_Allreduce(local_boundary.data(), boundary.data(),
                               static_cast<int>(boundary.size()), MPI_DOUBLE,
                               MPI_SUM, mpi.comm()) == MPI_SUCCESS);
    const double mass_abs = 2.0 * density * area;
    const std::array expected_absolute{
        mass_abs, mass_abs * cp * initial_temperature, mass_abs,
        mass_abs * 0.25};
    for (std::size_t quantity = 0; quantity < expected_absolute.size();
         ++quantity) {
      const double tolerance =
          128.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, expected_absolute[quantity]);
      HUNDUN_CHECK(std::abs(boundary[quantity * 2U]) <= tolerance);
      HUNDUN_CHECK_NEAR(boundary[quantity * 2U + 1U],
                        expected_absolute[quantity], tolerance);
      // Omitting either inlet or outlet leaves the other nonzero boundary
      // contribution and therefore falsifies the signed balance oracle.
      HUNDUN_CHECK(std::abs(boundary[8U + quantity]) > tolerance);
      HUNDUN_CHECK(std::abs(boundary[12U + quantity]) > tolerance);
      HUNDUN_CHECK(std::abs(boundary[8U + quantity] +
                            boundary[12U + quantity]) <= tolerance);
    }

    using OpenTransportAccess =
        hundun::flow::test::MaterialDensityTransportTestAccess;
    for (const int fault_rank : {0, mpi.size() - 1}) {
      for (const bool conservation : {false, true}) {
        const hundun::test::IdealGasStateSnapshot before_scalar_gate{
            state.snapshot(hundun::flow::FlowLayer::history),
            state.snapshot(hundun::flow::FlowLayer::committed),
            state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
            flow.closure_state()};
        OpenTransportAccess::reset();
        if (conservation)
          OpenTransportAccess::set_transport_conservation_defect(
              1U, 1.0e-9, fault_rank);
        else
          OpenTransportAccess::set_transport_residual(1U, 2.0e-9,
                                                      fault_rank);
        const auto scalar_gate_failed =
            flow.attempt(state, 0.0, bdf2, {}, {});
        OpenTransportAccess::reset();
        HUNDUN_CHECK(scalar_gate_failed.flow().flow().disposition ==
                     hundun::flow::StepAttemptDisposition::recoverable_failure);
        HUNDUN_CHECK(
            scalar_gate_failed.flow().flow().reason ==
            (conservation
                 ? hundun::flow::StepFailureReason::final_conservation_defect
                 : hundun::flow::StepFailureReason::final_transport_residual));
        HUNDUN_CHECK(scalar_gate_failed.flow().flow().lowest_failing_rank ==
                     fault_rank);
        HUNDUN_CHECK(scalar_gate_failed.closure_report_available());
        HUNDUN_CHECK(scalar_gate_failed.closure_report().stage() ==
                     hundun::flow::IdealGasClosureStage::provisional);
        HUNDUN_CHECK(
            hundun::flow::test::IdealGasClosureTestAccess::report_authenticated(
                scalar_gate_failed));
        require_halo_trace(
            flow, fields,
            {hundun::flow::IdealGasClosureStage::predictor,
             hundun::flow::IdealGasClosureStage::provisional});
        const hundun::test::IdealGasStateSnapshot after_scalar_gate{
            state.snapshot(hundun::flow::FlowLayer::history),
            state.snapshot(hundun::flow::FlowLayer::committed),
            state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
            flow.closure_state()};
        HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
            before_scalar_gate, after_scalar_gate));
      }
      if (mpi.size() == 1)
        break;
    }

    const hundun::test::IdealGasStateSnapshot before_backflow{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    hundun::flow::test::IdealGasClosureTestAccess::set_outlet_backflow_fault(
        flow);
    const auto backflow = flow.attempt(state, 0.0, bdf2, {}, {});
    HUNDUN_CHECK(backflow.flow().flow().disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(backflow.flow().flow().reason ==
                 hundun::flow::StepFailureReason::boundary_backflow);
    HUNDUN_CHECK(backflow.flow().flow().final_backflow_evidence.has_value());
    HUNDUN_CHECK(backflow.closure_report_available());
    HUNDUN_CHECK(backflow.closure_report().stage() ==
                 hundun::flow::IdealGasClosureStage::final);
    require_halo_trace(flow, fields,
                       {hundun::flow::IdealGasClosureStage::predictor,
                        hundun::flow::IdealGasClosureStage::provisional,
                        hundun::flow::IdealGasClosureStage::final});
    HUNDUN_TASK21_DIAGNOSTIC_STATUS(
        flow.closure_diagnostic_source(state, backflow), mpi,
        hundun::diagnostics::DiagnosticStatus::warning,
        hundun::diagnostics::DiagnosticFailureClass::none, "none", -1);
    const hundun::test::IdealGasStateSnapshot after_backflow{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(before_backflow,
                                                             after_backflow));
  }
  if (full_case && !reference_mode.empty()) {
    const auto candidate = make_reference(
        topology, state.snapshot(hundun::flow::FlowLayer::committed), mpi);
    if (reference_mode == "--reference-write") {
      HUNDUN_CHECK(mpi.size() == 1);
      if (mpi.rank() == 0)
        write_reference(reference_path, candidate);
    } else {
      HUNDUN_CHECK(reference_mode == "--reference-read");
      compare_reference(read_reference(reference_path), candidate);
    }
  }
  bool stale_rejected = false;
  try {
    static_cast<void>(closure_source.committed_step());
  } catch (const hundun::runtime::Error &error) {
    stale_rejected = std::string_view(error.what()) ==
                     "ideal-gas closure diagnostic source is stale";
  }
  HUNDUN_CHECK(stale_rejected);
  if (exhaust_generation) {
    const hundun::test::IdealGasStateSnapshot before_precedence{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    hundun::flow::test::IdealGasClosureTestAccess::
        set_candidate_precedence_fault(flow, mpi.size() - 1);
    const auto precedence_failed = flow.attempt(state, 0.0, bdf2, {}, {});
    HUNDUN_CHECK(precedence_failed.flow().flow().reason ==
                 hundun::flow::StepFailureReason::density_closure_failure);
    HUNDUN_CHECK(precedence_failed.closure_report().reason() ==
                 hundun::flow::IdealGasClosureFailureReason::
                     non_finite_enthalpy);
    HUNDUN_CHECK(precedence_failed.closure_report().lowest_failing_rank() ==
                 mpi.size() - 1);
    HUNDUN_CHECK(precedence_failed.closure_report().stage() ==
                 hundun::flow::IdealGasClosureStage::predictor);
    HUNDUN_CHECK(precedence_failed.closure_report().collective_count() == 3U);
    HUNDUN_TASK21_DIAGNOSTIC_STATUS(
        flow.closure_diagnostic_source(state, precedence_failed), mpi,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::non_finite_state,
        "closure.enthalpy-non-finite", mpi.size() - 1);
    const hundun::test::IdealGasStateSnapshot after_precedence{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
        before_precedence, after_precedence));

#ifdef HUNDUN_TASK21_DIAGNOSTICS_EXHAUSTIVE
    struct ClosureStatusCase final {
      hundun::flow::IdealGasClosureFailureReason reason;
      hundun::diagnostics::DiagnosticFailureClass classification;
      std::string_view code;
    };
    const std::array closure_status_cases{
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::invalid_input,
            hundun::diagnostics::DiagnosticFailureClass::invalid_input,
            "closure.invalid-input"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::non_finite_enthalpy,
            hundun::diagnostics::DiagnosticFailureClass::non_finite_state,
            "closure.enthalpy-non-finite"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::non_positive_enthalpy,
            hundun::diagnostics::DiagnosticFailureClass::non_positive_state,
            "closure.enthalpy-non-positive"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::
                non_finite_temperature,
            hundun::diagnostics::DiagnosticFailureClass::non_finite_state,
            "closure.temperature-non-finite"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::
                non_positive_temperature,
            hundun::diagnostics::DiagnosticFailureClass::non_positive_state,
            "closure.temperature-non-positive"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::denominator_breakdown,
            hundun::diagnostics::DiagnosticFailureClass::numerical_breakdown,
            "closure.denominator-breakdown"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::non_finite_pressure,
            hundun::diagnostics::DiagnosticFailureClass::non_finite_state,
            "closure.pressure-non-finite"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::non_positive_pressure,
            hundun::diagnostics::DiagnosticFailureClass::non_positive_state,
            "closure.pressure-non-positive"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::non_finite_density,
            hundun::diagnostics::DiagnosticFailureClass::non_finite_state,
            "closure.density-non-finite"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::non_positive_density,
            hundun::diagnostics::DiagnosticFailureClass::non_positive_state,
            "closure.density-non-positive"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::eos_residual,
            hundun::diagnostics::DiagnosticFailureClass::non_convergence,
            "closure.eos-residual"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::remap_residual,
            hundun::diagnostics::DiagnosticFailureClass::non_convergence,
            "closure.remap-residual"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::mass_conservation,
            hundun::diagnostics::DiagnosticFailureClass::conservation,
            "closure.mass-conservation"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::enthalpy_conservation,
            hundun::diagnostics::DiagnosticFailureClass::conservation,
            "closure.enthalpy-conservation"},
        ClosureStatusCase{
            hundun::flow::IdealGasClosureFailureReason::collective_operation,
            hundun::diagnostics::DiagnosticFailureClass::collective_operation,
            "closure.collective-operation"}};
    for (const int target : {0, mpi.size() - 1}) {
      for (const auto &status_case : closure_status_cases) {
        const hundun::test::IdealGasStateSnapshot before_status{
            state.snapshot(hundun::flow::FlowLayer::history),
            state.snapshot(hundun::flow::FlowLayer::committed),
            state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
            flow.closure_state()};
        hundun::flow::test::IdealGasClosureTestAccess::set_stage_failure(
            flow, hundun::flow::IdealGasClosureStage::predictor,
            status_case.reason, target);
        const auto status_report = flow.attempt(state, 0.0, bdf2, {}, {});
        HUNDUN_CHECK(status_report.closure_report_available());
        HUNDUN_CHECK(status_report.closure_report().reason() ==
                     status_case.reason);
        HUNDUN_CHECK(status_report.closure_report().lowest_failing_rank() ==
                     target);
        HUNDUN_TASK21_DIAGNOSTIC_STATUS(
            flow.closure_diagnostic_source(state, status_report), mpi,
            hundun::diagnostics::DiagnosticStatus::failed,
            status_case.classification, status_case.code, target);
        const hundun::test::IdealGasStateSnapshot after_status{
            state.snapshot(hundun::flow::FlowLayer::history),
            state.snapshot(hundun::flow::FlowLayer::committed),
            state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
            flow.closure_state()};
        HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
            before_status, after_status));
      }
      if (mpi.size() == 1)
        break;
    }
#endif

    {
      const hundun::test::IdealGasStateSnapshot before{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      const auto before_predictor = flow.attempt(
          state, std::numeric_limits<double>::quiet_NaN(), bdf2, {}, {});
      HUNDUN_CHECK(before_predictor.flow().flow().disposition ==
                   hundun::flow::StepAttemptDisposition::non_retryable_failure);
      HUNDUN_CHECK(before_predictor.flow().flow().reason ==
                   hundun::flow::StepFailureReason::invalid_input);
      HUNDUN_CHECK(!before_predictor.closure_report_available());
      HUNDUN_TASK21_DIAGNOSTIC_STATUS(
          flow.closure_diagnostic_source(state, before_predictor), mpi,
          hundun::diagnostics::DiagnosticStatus::unavailable,
          hundun::diagnostics::DiagnosticFailureClass::unavailable,
          "closure.not-evaluated", 0);
      const hundun::test::IdealGasStateSnapshot after{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(before, after));
    }

    const std::array outer_failure_cases{
        std::tuple{hundun::flow::test::IdealGasOuterFailurePoint::
                       momentum_after_predictor,
                   hundun::flow::IdealGasClosureStage::predictor, 1U, 5U},
        std::tuple{hundun::flow::test::IdealGasOuterFailurePoint::
                       pressure_after_first_corrector,
                   hundun::flow::IdealGasClosureStage::provisional, 2U, 9U}};
    for (const auto &[point, stage, evaluations, collectives] :
         outer_failure_cases) {
      const hundun::test::IdealGasStateSnapshot before{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      hundun::flow::test::IdealGasClosureTestAccess::set_outer_failure(
          flow, point, mpi.size() - 1);
      const auto outer_failed = flow.attempt(state, 0.0, bdf2, {}, {});
      HUNDUN_CHECK(outer_failed.flow().flow().disposition ==
                   hundun::flow::StepAttemptDisposition::non_retryable_failure);
      HUNDUN_CHECK(outer_failed.flow().flow().reason ==
                   hundun::flow::StepFailureReason::invalid_input);
      HUNDUN_CHECK(outer_failed.flow().flow().lowest_failing_rank ==
                   mpi.size() - 1);
      HUNDUN_CHECK(outer_failed.closure_report_available());
      HUNDUN_CHECK(outer_failed.closure_report().stage() == stage);
      HUNDUN_CHECK(outer_failed.closure_report().disposition() ==
                   hundun::flow::IdealGasClosureDisposition::closed);
      HUNDUN_CHECK(outer_failed.closure_report().evaluation_count() ==
                   evaluations);
      HUNDUN_CHECK(outer_failed.closure_report().collective_count() ==
                   collectives);
      HUNDUN_CHECK(
          outer_failed.closure_report().candidate_pressure_available());
      HUNDUN_CHECK_NEAR(
          outer_failed.closure_report().candidate_pressure_pa(),
          pressure * 315.0 / 300.0, 1.0e-12 * pressure);
      if (stage == hundun::flow::IdealGasClosureStage::predictor)
        require_halo_trace(
            flow, fields,
            {hundun::flow::IdealGasClosureStage::predictor});
      else
        require_halo_trace(
            flow, fields,
            {hundun::flow::IdealGasClosureStage::predictor,
             hundun::flow::IdealGasClosureStage::provisional});
      HUNDUN_TASK21_DIAGNOSTIC_STATUS(
          flow.closure_diagnostic_source(state, outer_failed), mpi,
          hundun::diagnostics::DiagnosticStatus::warning,
          hundun::diagnostics::DiagnosticFailureClass::none, "none", -1);
      const hundun::test::IdealGasStateSnapshot after{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(before, after));
    }

    const auto expect_stage_failure =
        [&](hundun::flow::IdealGasClosureStage stage,
            std::uint32_t expected_evaluations,
            std::uint32_t expected_collectives) {
          const hundun::test::IdealGasStateSnapshot before{
              state.snapshot(hundun::flow::FlowLayer::history),
              state.snapshot(hundun::flow::FlowLayer::committed),
              state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
              flow.closure_state()};
          hundun::flow::test::IdealGasClosureTestAccess::set_stage_failure(
              flow, stage,
              hundun::flow::IdealGasClosureFailureReason::
                  non_finite_temperature,
              mpi.size() - 1);
          const auto failed_stage = flow.attempt(state, 0.0, bdf2, {}, {});
          HUNDUN_CHECK(failed_stage.flow().flow().disposition ==
                       hundun::flow::StepAttemptDisposition::
                           recoverable_failure);
          HUNDUN_CHECK(failed_stage.flow().flow().reason ==
                       hundun::flow::StepFailureReason::
                           density_closure_failure);
          HUNDUN_CHECK(failed_stage.closure_report_available());
          HUNDUN_CHECK(failed_stage.closure_report().stage() == stage);
          HUNDUN_CHECK(failed_stage.closure_report().reason() ==
                       hundun::flow::IdealGasClosureFailureReason::
                           non_finite_temperature);
          HUNDUN_CHECK(
              failed_stage.closure_report().lowest_failing_rank() ==
              mpi.size() - 1);
          HUNDUN_CHECK(failed_stage.closure_report().evaluation_count() ==
                       expected_evaluations);
          HUNDUN_CHECK(failed_stage.closure_report().collective_count() ==
                       expected_collectives);
          if (stage == hundun::flow::IdealGasClosureStage::predictor) {
            HUNDUN_CHECK(!failed_stage.closure_report()
                              .candidate_pressure_available());
            require_halo_trace(flow, fields, {});
          } else {
            HUNDUN_CHECK(failed_stage.closure_report()
                             .candidate_pressure_available());
            HUNDUN_CHECK_NEAR(
                failed_stage.closure_report().candidate_pressure_pa(),
                pressure * 315.0 / 300.0, 1.0e-12 * pressure);
            if (stage == hundun::flow::IdealGasClosureStage::provisional)
              require_halo_trace(
                  flow, fields,
                  {hundun::flow::IdealGasClosureStage::predictor});
            else
              require_halo_trace(
                  flow, fields,
                  {hundun::flow::IdealGasClosureStage::predictor,
                   hundun::flow::IdealGasClosureStage::provisional});
          }
          HUNDUN_TASK21_DIAGNOSTIC_STATUS(
              flow.closure_diagnostic_source(state, failed_stage), mpi,
              hundun::diagnostics::DiagnosticStatus::failed,
              hundun::diagnostics::DiagnosticFailureClass::non_finite_state,
              "closure.temperature-non-finite", mpi.size() - 1);
          const hundun::test::IdealGasStateSnapshot after{
              state.snapshot(hundun::flow::FlowLayer::history),
              state.snapshot(hundun::flow::FlowLayer::committed),
              state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
              flow.closure_state()};
          HUNDUN_CHECK(
              hundun::test::ideal_gas_state_bitwise_equal(before, after));
        };
    expect_stage_failure(hundun::flow::IdealGasClosureStage::predictor, 1U,
                         3U);
    expect_stage_failure(hundun::flow::IdealGasClosureStage::provisional, 2U,
                         7U);
    expect_stage_failure(hundun::flow::IdealGasClosureStage::final, 3U, 11U);

    for (const bool enthalpy_density : {false, true}) {
      const hundun::test::IdealGasStateSnapshot before_post_store{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      hundun::flow::test::IdealGasClosureTestAccess::set_post_store_corruption(
          flow, mpi.size() - 1, enthalpy_density);
      const auto post_store_failed = flow.attempt(state, 0.0, bdf2, {}, {});
      HUNDUN_CHECK(post_store_failed.flow().flow().disposition ==
                   hundun::flow::StepAttemptDisposition::recoverable_failure);
      HUNDUN_CHECK(post_store_failed.flow().flow().reason ==
                   hundun::flow::StepFailureReason::density_closure_failure);
      HUNDUN_CHECK(post_store_failed.closure_report_available());
      HUNDUN_CHECK(post_store_failed.closure_report().stage() ==
                   hundun::flow::IdealGasClosureStage::final);
      HUNDUN_CHECK(post_store_failed.closure_report().evaluation_count() == 3U);
      HUNDUN_CHECK(post_store_failed.closure_report().collective_count() ==
                   14U);
      HUNDUN_CHECK(post_store_failed.closure_report().reason() ==
                   hundun::flow::IdealGasClosureFailureReason::eos_residual);
      HUNDUN_CHECK(post_store_failed.closure_report().lowest_failing_rank() ==
                   mpi.size() - 1);
      require_halo_trace(
          flow, fields,
          {hundun::flow::IdealGasClosureStage::predictor,
           hundun::flow::IdealGasClosureStage::provisional});
      HUNDUN_TASK21_DIAGNOSTIC_STATUS(
          flow.closure_diagnostic_source(state, post_store_failed), mpi,
          hundun::diagnostics::DiagnosticStatus::failed,
          hundun::diagnostics::DiagnosticFailureClass::non_convergence,
          "closure.eos-residual", mpi.size() - 1);
      const hundun::test::IdealGasStateSnapshot after_post_store{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
          before_post_store, after_post_store));
    }

    using PostAssessmentFault = hundun::flow::test::IdealGasPostAssessmentFault;
    for (const int fault_rank : {0, mpi.size() - 1}) {
      for (const auto &[fault, expected_reason] :
           {std::pair{PostAssessmentFault::non_finite_state,
                      hundun::flow::StepFailureReason::transport_failure},
            std::pair{PostAssessmentFault::non_positive_density,
                      hundun::flow::StepFailureReason::transport_failure},
            std::pair{
                PostAssessmentFault::final_transport_residual,
                hundun::flow::StepFailureReason::final_transport_residual},
            std::pair{
                PostAssessmentFault::final_conservation_defect,
                hundun::flow::StepFailureReason::final_conservation_defect}}) {
        const hundun::test::IdealGasStateSnapshot before_post_assessment{
            state.snapshot(hundun::flow::FlowLayer::history),
            state.snapshot(hundun::flow::FlowLayer::committed),
            state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
            flow.closure_state()};
        hundun::flow::test::IdealGasClosureTestAccess::
            set_post_assessment_fault(flow, fault, fault_rank);
        const auto post_assessment_failed =
            flow.attempt(state, 0.0, bdf2, {}, {});
        HUNDUN_CHECK(post_assessment_failed.flow().flow().disposition ==
                     hundun::flow::StepAttemptDisposition::recoverable_failure);
        HUNDUN_CHECK(post_assessment_failed.flow().flow().reason ==
                     expected_reason);
        HUNDUN_CHECK(post_assessment_failed.flow().flow().lowest_failing_rank ==
                     fault_rank);
        HUNDUN_CHECK(post_assessment_failed.flow().material_failure_reason() ==
                     hundun::flow::MaterialTransportFailureReason::none);
        HUNDUN_CHECK(post_assessment_failed.closure_report_available());
        HUNDUN_CHECK(post_assessment_failed.closure_report().disposition() ==
                     hundun::flow::IdealGasClosureDisposition::closed);
        HUNDUN_CHECK(
            hundun::flow::test::IdealGasClosureTestAccess::report_authenticated(
                post_assessment_failed));
        HUNDUN_CHECK(
            hundun::flow::test::IdealGasClosureTestAccess::
                post_eos_evidence_authenticated(post_assessment_failed.flow()));
        require_halo_trace(flow, fields,
                           {hundun::flow::IdealGasClosureStage::predictor,
                            hundun::flow::IdealGasClosureStage::provisional,
                            hundun::flow::IdealGasClosureStage::final});
        for (std::uint8_t mutation = 0U;
             mutation <
             static_cast<std::uint8_t>(
                 hundun::flow::test::IdealGasPostEvidenceMutation::count);
             ++mutation)
          if (!hundun::flow::test::IdealGasClosureTestAccess::
                  post_evidence_mutation_rejected(
                      post_assessment_failed,
                      static_cast<
                          hundun::flow::test::IdealGasPostEvidenceMutation>(
                          mutation)))
            throw std::runtime_error(
                "post-EOS failure mutation was accepted: " +
                std::to_string(mutation));
        const hundun::test::IdealGasStateSnapshot after_post_assessment{
            state.snapshot(hundun::flow::FlowLayer::history),
            state.snapshot(hundun::flow::FlowLayer::committed),
            state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
            flow.closure_state()};
        HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
            before_post_assessment, after_post_assessment));
      }
      if (mpi.size() == 1)
        break;
    }

    using TerminalMode =
        hundun::flow::test::MaterialTerminalModeForTest;
    using TerminalPoint =
        hundun::flow::test::MaterialTerminalPointForTest;
    using MaterialAccess =
        hundun::flow::test::MaterialDensityPisoTestAccess;
    constexpr std::array final_terminal_points{
        TerminalPoint::final_continuity_reduction,
        TerminalPoint::final_continuity_status,
        TerminalPoint::final_pressure_entry,
        TerminalPoint::final_pressure_gamma_sum,
        TerminalPoint::final_pressure_gamma_count,
        TerminalPoint::final_pressure_residual_reduction,
        TerminalPoint::final_momentum_residual_reduction,
        TerminalPoint::final_momentum_conservation_reduction,
        TerminalPoint::final_momentum_status,
        TerminalPoint::final_conservation_status};
    for (const auto point : final_terminal_points) {
      const hundun::test::IdealGasStateSnapshot before_terminal{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      MaterialAccess::reset_terminal_fault();
      MaterialAccess::set_terminal_fault(point,
                                         TerminalMode::returned_rankless);
      bool typed = false;
      try {
        static_cast<void>(flow.attempt(state, 0.0, bdf2, {}, {}));
      } catch (const hundun::runtime::MpiOperationError &error) {
        typed = std::string_view(error.what()) ==
                "material-density flow collective failure has no reliable "
                "rank";
      }
      const auto terminal_calls = MaterialAccess::terminal_point_calls(point);
      MaterialAccess::reset_terminal_fault();
      HUNDUN_CHECK(typed);
      require_halo_trace(
          flow, fields,
          {hundun::flow::IdealGasClosureStage::predictor,
           hundun::flow::IdealGasClosureStage::provisional,
           hundun::flow::IdealGasClosureStage::final});
      const hundun::test::IdealGasStateSnapshot after_terminal{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
          before_terminal, after_terminal));
      HUNDUN_CHECK(terminal_calls == 1U);
    }

    using TransportAccess =
        hundun::flow::test::MaterialDensityTransportTestAccess;
    enum class FinalTransportFault { enthalpy_residual, mass, enthalpy_balance };
    for (const int fault_rank : {0, mpi.size() - 1}) {
      for (const auto fault : {FinalTransportFault::enthalpy_residual,
                               FinalTransportFault::mass,
                               FinalTransportFault::enthalpy_balance}) {
        const hundun::test::IdealGasStateSnapshot before_gate{
            state.snapshot(hundun::flow::FlowLayer::history),
            state.snapshot(hundun::flow::FlowLayer::committed),
            state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
            flow.closure_state()};
        TransportAccess::reset();
        if (fault == FinalTransportFault::enthalpy_residual)
          TransportAccess::set_transport_residual(0U, 2.0e-9, fault_rank);
        else if (fault == FinalTransportFault::mass)
          TransportAccess::set_mass_conservation_defect(1.0e-9, fault_rank);
        else
          TransportAccess::set_transport_conservation_defect(0U, 1.0e-9,
                                                              fault_rank);
        const auto gate_failed = flow.attempt(state, 0.0, bdf2, {}, {});
        TransportAccess::reset();
        HUNDUN_CHECK(gate_failed.flow().flow().disposition ==
                     hundun::flow::StepAttemptDisposition::recoverable_failure);
        HUNDUN_CHECK(gate_failed.flow().flow().reason ==
                     (fault == FinalTransportFault::enthalpy_residual
                          ? hundun::flow::StepFailureReason::
                                final_transport_residual
                          : hundun::flow::StepFailureReason::
                                final_conservation_defect));
        HUNDUN_CHECK(gate_failed.flow().flow().lowest_failing_rank ==
                     fault_rank);
        HUNDUN_CHECK(gate_failed.closure_report_available());
        HUNDUN_CHECK(gate_failed.closure_report().stage() ==
                     hundun::flow::IdealGasClosureStage::provisional);
        HUNDUN_CHECK(gate_failed.closure_report().disposition() ==
                     hundun::flow::IdealGasClosureDisposition::closed);
        HUNDUN_CHECK(
            hundun::flow::test::IdealGasClosureTestAccess::report_authenticated(
                gate_failed));
        require_halo_trace(
            flow, fields,
            {hundun::flow::IdealGasClosureStage::predictor,
             hundun::flow::IdealGasClosureStage::provisional});
        const hundun::test::IdealGasStateSnapshot after_gate{
            state.snapshot(hundun::flow::FlowLayer::history),
            state.snapshot(hundun::flow::FlowLayer::committed),
            state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
            flow.closure_state()};
        HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(before_gate,
                                                                 after_gate));
      }
      if (mpi.size() == 1)
        break;
    }

    for (const auto prepare_fault :
         {hundun::flow::test::IdealGasPrepareFault::state_prepare,
          hundun::flow::test::IdealGasPrepareFault::closure_prepare}) {
      const hundun::test::IdealGasStateSnapshot before_prepare{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      hundun::flow::test::IdealGasClosureTestAccess::set_prepare_fault(
          flow, prepare_fault, mpi.size() - 1);
      const auto prepare_failed = flow.attempt(state, 0.0, bdf2, {}, {});
      HUNDUN_CHECK(prepare_failed.flow().flow().disposition ==
                   hundun::flow::StepAttemptDisposition::
                       non_retryable_failure);
      HUNDUN_CHECK(prepare_failed.flow().flow().reason ==
                   hundun::flow::StepFailureReason::invalid_input);
      HUNDUN_CHECK(prepare_failed.flow().flow().lowest_failing_rank ==
                   mpi.size() - 1);
      const hundun::test::IdealGasStateSnapshot after_prepare{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
          before_prepare, after_prepare));
    }

    for (const int fault_rank : {0, mpi.size() - 1}) {
      const hundun::test::IdealGasStateSnapshot before_layout{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      hundun::flow::test::IdealGasClosureTestAccess::set_attempt_layout_fault(
          flow, fault_rank);
      const auto layout_failed = flow.attempt(state, 0.0, bdf2, {}, {});
      HUNDUN_CHECK(layout_failed.flow().flow().disposition ==
                   hundun::flow::StepAttemptDisposition::non_retryable_failure);
      HUNDUN_CHECK(layout_failed.flow().flow().reason ==
                   hundun::flow::StepFailureReason::invalid_input);
      HUNDUN_CHECK(layout_failed.flow().flow().lowest_failing_rank ==
                   fault_rank);
      HUNDUN_CHECK(!layout_failed.closure_report_available());
      const hundun::test::IdealGasStateSnapshot after_layout{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(before_layout,
                                                               after_layout));
    }

    for (const int fault_rank : {0, mpi.size() - 1}) {
      const hundun::test::IdealGasStateSnapshot before_mpi_failure{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      hundun::flow::test::IdealGasClosureTestAccess::set_post_store_mpi_fault(
          flow, fault_rank);
      bool typed_mpi_failure = false;
      try {
        static_cast<void>(flow.attempt(state, 0.0, bdf2, {}, {}));
      } catch (const hundun::runtime::MpiOperationError &error) {
        typed_mpi_failure =
            std::string_view(error.what()).find("post-store reduction") !=
            std::string_view::npos;
      }
      HUNDUN_CHECK(typed_mpi_failure);
      const hundun::test::IdealGasStateSnapshot after_mpi_failure{
          state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          flow.closure_state()};
      HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
          before_mpi_failure, after_mpi_failure));
    }

    using PreparationFault =
        hundun::flow::test::IdealGasAttemptPreparationFault;
    const std::array preparation_faults{
        PreparationFault::predictor_local,
        PreparationFault::predictor_write_capability,
        PreparationFault::provisional_local,
        PreparationFault::provisional_write_capability,
        PreparationFault::final_local,
        PreparationFault::final_write_capability,
        PreparationFault::final_readback,
        PreparationFault::post_density_views,
        PreparationFault::post_geometry,
        PreparationFault::post_transport_views,
        PreparationFault::post_next_integrals,
        PreparationFault::post_report_finalization,
        PreparationFault::pre_authority_preparation,
        PreparationFault::post_authority_preparation};
    for (const int fault_rank : {0, mpi.size() - 1}) {
      for (const auto fault : preparation_faults) {
        const hundun::test::IdealGasStateSnapshot before_preparation_fault{
            state.snapshot(hundun::flow::FlowLayer::history),
            state.snapshot(hundun::flow::FlowLayer::committed),
            state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
            flow.closure_state()};
        hundun::flow::test::IdealGasClosureTestAccess::
            set_attempt_preparation_fault(flow, fault, fault_rank);
        const auto preparation_failed =
            flow.attempt(state, 0.0, bdf2, {}, {});
        HUNDUN_CHECK(preparation_failed.flow().flow().disposition !=
                     hundun::flow::StepAttemptDisposition::committed);
        HUNDUN_CHECK(
            preparation_failed.flow().flow().lowest_failing_rank ==
            fault_rank);
        const hundun::test::IdealGasStateSnapshot after_preparation_fault{
            state.snapshot(hundun::flow::FlowLayer::history),
            state.snapshot(hundun::flow::FlowLayer::committed),
            state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
            flow.closure_state()};
        HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
            before_preparation_fault, after_preparation_fault));
      }
      const auto clean_successor =
          flow.attempt(state, 0.0, bdf2, {}, {});
      HUNDUN_CHECK(clean_successor.flow().flow().disposition ==
                   hundun::flow::StepAttemptDisposition::committed);
      if (mpi.size() == 1)
        break;
    }

    const hundun::test::IdealGasStateSnapshot before_failure{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    hundun::flow::test::IdealGasClosureTestAccess::set_uniform_enthalpy_rate(
        flow, -cp * 1000.0 / dt);
    const auto failed = flow.attempt(state, 0.0, bdf2, {}, {});
    HUNDUN_CHECK(failed.flow().flow().disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(failed.flow().flow().reason ==
                 hundun::flow::StepFailureReason::density_closure_failure);
    HUNDUN_CHECK(failed.closure_report_available());
    HUNDUN_CHECK(failed.closure_report().stage() ==
                 hundun::flow::IdealGasClosureStage::predictor);
    const hundun::test::IdealGasStateSnapshot after_failure{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(before_failure,
                                                             after_failure));
    hundun::flow::test::IdealGasClosureTestAccess::set_uniform_enthalpy_rate(
        flow, cp * 5.0 / dt);

    const hundun::test::IdealGasStateSnapshot before_finalization_wrap{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    hundun::flow::test::IdealGasClosureTestAccess::
        force_finalization_identity_wrap(flow);
    const auto finalization_wrap = flow.attempt(state, 0.0, bdf2, {}, {});
    HUNDUN_CHECK(finalization_wrap.flow().flow().disposition ==
                 hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(finalization_wrap.flow().flow().reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(finalization_wrap.flow().flow().lowest_failing_rank == 0);
    HUNDUN_CHECK(finalization_wrap.closure_report_available());
    HUNDUN_CHECK(finalization_wrap.closure_report().stage() ==
                 hundun::flow::IdealGasClosureStage::provisional);
    const hundun::test::IdealGasStateSnapshot after_finalization_wrap{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
        before_finalization_wrap, after_finalization_wrap));

    const hundun::test::IdealGasStateSnapshot before_attempt_wrap{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    hundun::flow::test::MaterialDensityTransportTestAccess::
        force_attempt_identity_wrap(state);
    const auto attempt_wrap = flow.attempt(state, 0.0, bdf2, {}, {});
    HUNDUN_CHECK(attempt_wrap.flow().flow().disposition ==
                 hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(attempt_wrap.flow().flow().reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(attempt_wrap.flow().flow().lowest_failing_rank == 0);
    HUNDUN_CHECK(!attempt_wrap.closure_report_available());
    const hundun::test::IdealGasStateSnapshot after_attempt_wrap{
        state.snapshot(hundun::flow::FlowLayer::history),
        state.snapshot(hundun::flow::FlowLayer::committed),
        state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
        flow.closure_state()};
    HUNDUN_CHECK(hundun::test::ideal_gas_state_bitwise_equal(
        before_attempt_wrap, after_attempt_wrap));

    hundun::flow::test::IdealGasClosureTestAccess::exhaust_source_generation(
        flow);
    bool exhausted_rejected = false;
    try {
      static_cast<void>(flow.attempt(state, 0.0, bdf2, {}, {}));
    } catch (const hundun::runtime::Error &) {
      exhausted_rejected = true;
    }
    HUNDUN_CHECK(exhausted_rejected);
  }
}

} // namespace

#ifndef HUNDUN_TASK21_IDEAL_GAS_FIXTURE_ONLY
int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    const bool open_domain =
        argc > 1 && std::string_view(argv[1]) == "--open-plug";
    const bool full_case = argc > 1;
    HUNDUN_CHECK(!full_case || argc == 4);
    run(mpi, open_domain, argc == 1, full_case,
        full_case ? std::string_view(argv[2]) : std::string_view{},
        full_case ? std::filesystem::path(argv[3]) : std::filesystem::path{});
  });
}
#endif
