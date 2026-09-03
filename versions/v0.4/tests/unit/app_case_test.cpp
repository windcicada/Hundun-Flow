// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_case.hpp"
#include "hundun/v04_field.hpp"

#include "app_case_detail.hpp"

#include <mpi.h>

#include <unistd.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace {

namespace fs = std::filesystem;
using hundun::v04::CaseCompiler;
using hundun::v04::CouplingKind;
using hundun::v04::FieldId;
using hundun::v04::FieldRegistry;
using hundun::v04::FieldSchema;
using hundun::v04::GeometryKind;
using hundun::v04::IbmReconstructionPolicy;
using hundun::v04::LinearAlgorithm;
using hundun::v04::MgCorrectionScaling;
using hundun::v04::PressureReferenceKind;
using hundun::v04::Status;
using hundun::v04::StatusCode;
using hundun::v04::TimeControlKind;
using hundun::v04::TimeScheme;
using hundun::v04::TransportLaw;
using hundun::v04::TransportedScalarRole;
using hundun::v04::TurbulenceKind;
using hundun::v04::ValidatedModel;

constexpr std::string_view kUniformMesh = R"json({
  "kind":"uniform",
  "domain":{"lower":[0,0,0],"upper":[1,2,3]},
  "exact_cells":[8,8,6],
  "base_spacing":null,
  "minimum_spacing":[0.1,0.2,0.4],
  "max_growth_ratio":1.0,
  "focus_regions":[],
  "limits":{"max_global_cells":4096,
            "max_memory_bytes_per_rank":67108864},
  "data_files":[],"immersed_boundary":null
})json";

constexpr std::string_view kTensorMesh = R"json({
  "kind":"tensor_stretched",
  "domain":{"lower":[0,0,0],"upper":[1,1,1]},
  "exact_cells":null,
  "base_spacing":[0.25,0.25,0.25],
  "minimum_spacing":[0.05,0.05,0.05],
  "max_growth_ratio":1.2,
  "focus_regions":[
    {"lower":[0.6,0.1,0.1],"upper":[1.2,0.5,0.5],
     "target_spacing":[0.2,0.15,0.12]},
    {"lower":[-1,0.2,0.2],"upper":[0.4,0.6,0.6],
     "target_spacing":[0.1,0.1,0.1]},
    {"lower":[-2,0.2,0.2],"upper":[0.4,0.6,0.6],
     "target_spacing":[0.1,0.1,0.1]}
  ],
  "limits":{"max_global_cells":100000,
            "max_memory_bytes_per_rank":134217728},
  "data_files":[],"immersed_boundary":null
})json";

constexpr std::string_view kCoastRuntimeAxesMesh = R"json({
  "kind":"coast_runtime_axes_v1",
  "axes_file":"axes.dat",
  "domain":{"lower":[0,0,0],"upper":[1,1,1]},
  "exact_cells":[2,2,1],
  "base_spacing":[0.5,0.5,1.0],
  "minimum_spacing":[0.01,0.01,0.01],
  "max_growth_ratio":1.2,
  "focus_regions":[],
  "limits":{"max_global_cells":16,
            "max_memory_bytes_per_rank":1048576},
  "data_files":[],"immersed_boundary":null
})json";

constexpr std::string_view kCoastRuntimeAxes = R"data(
COAST_RUNTIME_AXES 1
grid 2 2 1
x 3
0 0.1000000001 1
y 3
0 0.4 1
z 2
0 1
)data";

constexpr std::string_view kPlaceholderThermophysics = R"data(
HUNDUN_THERMOPHYSICS_V1
temperature_bounds 200 3000
temperature_inversion 1e-12 32
closed_mass_newton 1e-12 24 0.25
species_count 1
species air
molecular_weight 28.96546
temperature_switch 1000
nasa7_low 3.5 0 0 0 0 0 0
nasa7_high 3.5 0 0 0 0 0 0
transport_sutherland 1.716e-5 273.15 110.4 0.71
end_species
end
)data";

constexpr std::string_view kCoastNativeAirThermophysics = R"data(
HUNDUN_THERMOPHYSICS_V1
temperature_bounds 273.15 6000
temperature_inversion 1e-12 64
closed_mass_newton 1e-12 32 0.2
species_count 1
species air
molecular_weight 28.850334
temperature_switch 1000
nasa7_low 3.5838100068 -7.2700635412e-4 1.67056387003e-6 -1.091801341e-10 -4.317787988e-13 -1050.5394088 3.1124135035
nasa7_high 3.1013370688 1.24138813631e-3 -4.1882038804e-7 6.641656204e-11 -3.9127843272e-15 -985.27467132 5.3560174057
transport_coast_native_air
end_species
end
)data";

std::string case_json(std::string_view mesh,
                      std::string_view flow = R"json({
  "model":"single_phase_low_mach_compressible",
  "pressure_reference":"boundary_absolute",
  "reacting":false
})json",
                      std::string_view solver = R"json({
  "coupling":"PISO","pressure_correctors":2
})json",
                      std::string_view turbulence =
                          R"json({"model":"none"})json",
                      std::string_view time =
                          R"json({"control":"fixed","scheme":"backward_euler","initial_dt":0.001,"minimum_dt":1e-8,"maximum_dt":0.1,"convective_cfl":0.8,"viscous_cfl":0.5,"thermal_cfl":0.5,"species_cfl":0.5,"acoustic_cfl":0.8,"maximum_growth":1.25,"retry_factor":0.5,"maximum_retries":8,"minimum_bdf_ratio":0.2,"maximum_bdf_ratio":5.0})json") {
  constexpr std::string_view boundaries = R"json({
    "x_min":{"flow_kind":"velocity_inlet","thermal_kind":"none","velocity":[1,0,0],"direction":[1,0,0],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]},
    "x_max":{"flow_kind":"pressure_outlet","thermal_kind":"none","velocity":[0,0,0],"direction":[1,0,0],"backflow_velocity":[-1,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":true,"scalars":[]},
    "y_min":{"flow_kind":"symmetry","thermal_kind":"none","velocity":[0,0,0],"direction":[0,1,0],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]},
    "y_max":{"flow_kind":"symmetry","thermal_kind":"none","velocity":[0,0,0],"direction":[0,1,0],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]},
    "z_min":{"flow_kind":"symmetry","thermal_kind":"none","velocity":[0,0,0],"direction":[0,0,1],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]},
    "z_max":{"flow_kind":"symmetry","thermal_kind":"none","velocity":[0,0,0],"direction":[0,0,1],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[]}
  })json";
  constexpr std::string_view schemes = R"json({"momentum":"limited_central2","enthalpy":"limited_central2","species":"tvd2","passive_scalar":"tvd2","diffusion":"central2","limiter":1.0})json";
  return std::string{"{\"schema_version\":1,\"units\":\"SI\",\"mesh\":"} +
         std::string(mesh) + ",\"flow\":" + std::string(flow) +
         ",\"solver\":" + std::string(solver) +
         (turbulence.empty()
              ? std::string{}
              : ",\"turbulence\":" + std::string(turbulence)) +
         ",\"thermophysics\":{\"data_file\":\"thermophysics.d\"}" +
         ",\"transported_scalars\":[]" +
         ",\"boundaries\":" + std::string(boundaries) +
         ",\"schemes\":" + std::string(schemes) +
         ",\"time\":" + std::string(time) + "}";
}

class ScratchCase {
 public:
  explicit ScratchCase(std::string_view label) {
    static std::uint64_t sequence = 0;
    root_ = fs::temp_directory_path() /
            ("hundun-v04-case-" + std::to_string(::getpid()) + "-" +
             std::string(label) + "-" +
             std::to_string(++sequence));
    std::error_code error;
    fs::remove_all(root_, error);
    fs::create_directories(root_);
  }
  ~ScratchCase() {
    std::error_code error;
    fs::remove_all(root_, error);
  }
  const fs::path& root() const noexcept { return root_; }
  void write(std::string_view name, std::string_view bytes) const {
    std::ofstream output(root_ / name, std::ios::binary);
    output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  }

 private:
  fs::path root_;
};

fs::path g_swap_target;
fs::path g_swap_source;
int g_swap_open_count = 0;

void swap_reference_before_open(int) {
  ++g_swap_open_count;
  if (g_swap_open_count != 2) {
    return;
  }
  std::error_code error;
  fs::remove(g_swap_target, error);
  error.clear();
  fs::create_symlink(g_swap_source, g_swap_target, error);
}

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

Status compile(const fs::path& root, ValidatedModel& model) {
  return CaseCompiler::load_and_compile(MPI_COMM_SELF, root, model);
}

bool replace_once(std::string& text, std::string_view from,
                  std::string_view to) {
  const std::size_t position = text.find(from);
  if (position == std::string::npos) {
    return false;
  }
  text.replace(position, from.size(), to);
  return true;
}

void append_wire_real(std::vector<std::uint8_t>& bytes, double value) {
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    bytes.push_back(static_cast<std::uint8_t>((bits >> shift) & 0xffU));
  }
}

bool replace_unique_wire_real(std::vector<std::uint8_t>& payload, double from,
                              double to) {
  std::vector<std::uint8_t> source;
  std::vector<std::uint8_t> replacement;
  append_wire_real(source, from);
  append_wire_real(replacement, to);
  const auto found =
      std::search(payload.begin(), payload.end(), source.begin(), source.end());
  if (found == payload.end() ||
      std::search(std::next(found), payload.end(), source.begin(),
                  source.end()) != payload.end()) {
    return false;
  }
  std::copy(replacement.begin(), replacement.end(), found);
  return true;
}

bool forge_legacy_native_air_wire(std::vector<std::uint8_t>& payload,
                                  double canonical_high_a6) {
  std::vector<std::uint8_t> transport;
  transport.push_back(static_cast<std::uint8_t>(TransportLaw::sutherland));
  append_wire_real(transport, 1.23456789e-5);
  append_wire_real(transport, 321.25);
  append_wire_real(transport, 111.75);
  append_wire_real(transport, 0.73);
  append_wire_real(transport, 0.0);
  const auto found = std::search(payload.begin(), payload.end(),
                                 transport.begin(), transport.end());
  if (found == payload.end() ||
      std::search(std::next(found), payload.end(), transport.begin(),
                  transport.end()) != payload.end()) {
    return false;
  }
  *found = static_cast<std::uint8_t>(TransportLaw::coast_native_air);
  std::fill(std::next(found), found + transport.size(), 0U);
  return replace_unique_wire_real(payload, canonical_high_a6, -985.27467132);
}

bool compile_json_fingerprint(std::string_view label, std::string_view json,
                              std::uint64_t& out) {
  ScratchCase scratch(label);
  scratch.write("case.json", json);
  scratch.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel model;
  const Status status = compile(scratch.root(), model);
  if (!expect(static_cast<bool>(status), label)) {
    return false;
  }
  out = model.fingerprint;
  return true;
}

bool rejects(std::string_view label, std::string_view json,
             StatusCode expected = StatusCode::invalid_case) {
  ScratchCase scratch(label);
  scratch.write("case.json", json);
  scratch.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel model;
  model.fingerprint = 987654321U;
  const Status status = compile(scratch.root(), model);
  return expect(status.code == expected && model.fingerprint == 987654321U,
                label);
}

bool same_mesh(const ValidatedModel& left, const ValidatedModel& right) {
  if (left.mesh.kind != right.mesh.kind ||
      left.mesh.lower.x != right.mesh.lower.x ||
      left.mesh.lower.y != right.mesh.lower.y ||
      left.mesh.lower.z != right.mesh.lower.z ||
      left.mesh.upper.x != right.mesh.upper.x ||
      left.mesh.upper.y != right.mesh.upper.y ||
      left.mesh.upper.z != right.mesh.upper.z ||
      left.mesh.has_exact_cells != right.mesh.has_exact_cells ||
      left.mesh.exact_cells.x != right.mesh.exact_cells.x ||
      left.mesh.exact_cells.y != right.mesh.exact_cells.y ||
      left.mesh.exact_cells.z != right.mesh.exact_cells.z ||
      left.mesh.has_base_spacing != right.mesh.has_base_spacing ||
      left.mesh.base_spacing.x != right.mesh.base_spacing.x ||
      left.mesh.base_spacing.y != right.mesh.base_spacing.y ||
      left.mesh.base_spacing.z != right.mesh.base_spacing.z ||
      left.mesh.minimum_spacing.x != right.mesh.minimum_spacing.x ||
      left.mesh.minimum_spacing.y != right.mesh.minimum_spacing.y ||
      left.mesh.minimum_spacing.z != right.mesh.minimum_spacing.z ||
      left.mesh.max_growth_ratio != right.mesh.max_growth_ratio ||
      left.mesh.limits.max_global_cells !=
          right.mesh.limits.max_global_cells ||
      left.mesh.limits.max_memory_bytes_per_rank !=
          right.mesh.limits.max_memory_bytes_per_rank ||
      left.mesh.focus_regions.size() != right.mesh.focus_regions.size()) {
    return false;
  }
  for (std::size_t i = 0; i < left.mesh.focus_regions.size(); ++i) {
    const auto& a = left.mesh.focus_regions[i];
    const auto& b = right.mesh.focus_regions[i];
    if (a.lower.x != b.lower.x || a.lower.y != b.lower.y ||
        a.lower.z != b.lower.z || a.upper.x != b.upper.x ||
        a.upper.y != b.upper.y || a.upper.z != b.upper.z ||
        a.target_spacing.x != b.target_spacing.x ||
        a.target_spacing.y != b.target_spacing.y ||
        a.target_spacing.z != b.target_spacing.z) {
      return false;
    }
  }
  return true;
}

bool compile_fingerprint(std::string_view label, std::string_view mesh,
                         std::uint64_t& out) {
  ScratchCase scratch(label);
  scratch.write("case.json", case_json(mesh));
  scratch.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel model;
  const Status status = compile(scratch.root(), model);
  if (!expect(static_cast<bool>(status), label)) {
    return false;
  }
  out = model.fingerprint;
  return true;
}

bool test_valid_fixture() {
  const char* data_root = std::getenv("HUNDUN_V04_TEST_DATA");
  if (data_root == nullptr) {
    return expect(false, "test data path is configured");
  }
  std::ifstream input(fs::path(data_root) / "case_minimal_valid.json",
                      std::ios::binary);
  const std::string json{std::istreambuf_iterator<char>(input),
                         std::istreambuf_iterator<char>()};
  ScratchCase scratch("valid-fixture");
  scratch.write("case.json", json);
  std::ifstream thermo_input(fs::path(data_root) / "thermophysics.d",
                             std::ios::binary);
  const std::string thermo{std::istreambuf_iterator<char>(thermo_input),
                           std::istreambuf_iterator<char>()};
  scratch.write("thermophysics.d", thermo);
  ValidatedModel model;
  const Status status = compile(scratch.root(), model);
  bool passed = expect(static_cast<bool>(status), "minimal fixture compiles");
  passed &= expect(model.mesh.kind == GeometryKind::uniform,
                   "fixture selects uniform geometry");
  passed &= expect(model.mesh.has_exact_cells &&
                       model.mesh.exact_cells.x == 8 &&
                       model.mesh.exact_cells.y == 8 &&
                       model.mesh.exact_cells.z == 8,
                   "fixture publishes exact cells");
  passed &= expect(!model.mesh.has_base_spacing &&
                       model.mesh.minimum_spacing.x == 0.125 &&
                       model.mesh.max_growth_ratio == 1.0,
                   "fixture publishes uniform spacing controls");
  passed &= expect(model.thermophysics.data_file ==
                           fs::path("thermophysics.d") &&
                       model.thermophysics.minimum_temperature == 200.0 &&
                       model.thermophysics.maximum_temperature == 3000.0 &&
                       model.thermophysics.species.size() == 1U &&
                       model.thermophysics.species[0].stable_name == "air" &&
                       model.thermophysics.species[0].molecular_weight ==
                           28.96546 &&
                       model.thermophysics.species[0].transport_law ==
                           hundun::v04::TransportLaw::sutherland,
                   "fixture publishes compact typed thermophysical data");
  passed &= expect(model.mesh.limits.max_global_cells == 4096 &&
                       model.mesh.limits.max_memory_bytes_per_rank == 67108864,
                   "fixture publishes hard limits");
  passed &= expect(model.transported_scalars.empty(),
                   "fixture publishes an explicit empty scalar catalog");
  passed &= expect(model.solver.pressure.absolute_tolerance == 1.0e-13 &&
                       model.solver.pressure.relative_tolerance == 1.0e-13 &&
                       model.solver.pressure.maximum_iterations == 400U &&
                       model.solver.pressure.true_residual_interval == 4U &&
                       model.solver.pressure.krylov_restart == 12U &&
                       model.solver.terminal.eos == 1.0e-10 &&
                       model.solver.terminal.continuity == 1.0e-10 &&
                       model.solver.terminal.closed_mass == 1.0e-10 &&
                       model.solver.terminal.gauge == 1.0e-10,
                   "fixture publishes explicit solver work and terminal gates");
  passed &= expect(model.fingerprint != 0U, "fixture has fingerprint");
  return passed;
}

bool test_tensor_normalization_and_fingerprint() {
  ScratchCase first_case("tensor-first");
  first_case.write("case.json", case_json(kTensorMesh));
  first_case.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel first;
  const Status first_status = compile(first_case.root(), first);
  bool passed = expect(static_cast<bool>(first_status),
                       "tensor mesh with optional exact cells omitted compiles");
  passed &= expect(first.mesh.kind == GeometryKind::tensor_stretched &&
                       !first.mesh.has_exact_cells &&
                       first.mesh.has_base_spacing,
                   "tensor sizing flags are typed");
  passed &= expect(first.mesh.focus_regions.size() == 2U,
                   "focus regions are clipped and deduplicated");
  passed &= expect(first.mesh.focus_regions[0].lower.x == 0.0 &&
                       first.mesh.focus_regions[0].upper.x == 0.4 &&
                       first.mesh.focus_regions[1].lower.x == 0.6 &&
                       first.mesh.focus_regions[1].upper.x == 1.0,
                   "focus regions are clipped and sorted");

  constexpr std::string_view reordered = R"json({
    "kind":"tensor_stretched",
    "domain":{"lower":[-0.0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":[0.25,0.25,0.25],
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.2,
    "focus_regions":[
      {"lower":[-2,0.2,0.2],"upper":[0.4,0.6,0.6],"target_spacing":[0.1,0.1,0.1]},
      {"lower":[0.6,0.1,0.1],"upper":[1.2,0.5,0.5],"target_spacing":[0.2,0.15,0.12]}
    ],
    "limits":{"max_global_cells":100000,"max_memory_bytes_per_rank":134217728},
    "data_files":[],"immersed_boundary":null
  })json";
  ScratchCase second_case("tensor-second");
  second_case.write("case.json", case_json(reordered));
  second_case.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel second;
  const Status second_status = compile(second_case.root(), second);
  passed &= expect(static_cast<bool>(second_status) && same_mesh(first, second),
                   "equivalent focus inputs produce one canonical typed mesh");
  passed &= expect(first.fingerprint == second.fingerprint,
                   "focus order duplicates clipping and negative zero do not change fingerprint");

  constexpr std::string_view exact_tensor = R"json({
    "kind":"tensor_stretched",
    "domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,
    "focus_regions":[],
    "limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},
    "data_files":[],"immersed_boundary":null
  })json";
  ScratchCase exact_case("tensor-exact");
  exact_case.write("case.json", case_json(exact_tensor));
  exact_case.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel exact;
  const Status exact_status = compile(exact_case.root(), exact);
  passed &= expect(static_cast<bool>(exact_status) && exact.mesh.has_exact_cells &&
                       exact.mesh.has_base_spacing &&
                       exact.mesh.exact_cells.x == 12,
                   "tensor mesh accepts exact cells together with base spacing");
  passed &= expect(exact.fingerprint != first.fingerprint,
                   "typed mesh changes affect fingerprint");
  return passed;
}

bool test_coast_runtime_axes_case() {
  ScratchCase scratch("coast-runtime-axes");
  scratch.write("case.json", case_json(kCoastRuntimeAxesMesh));
  scratch.write("thermophysics.d", kPlaceholderThermophysics);
  scratch.write("axes.dat", kCoastRuntimeAxes);
  ValidatedModel model;
  const Status status = compile(scratch.root(), model);
  bool passed = expect(static_cast<bool>(status),
                       "COAST runtime axes case compiles");
  passed &= expect(
      status && model.mesh.kind == GeometryKind::coast_runtime_axes_v1 &&
          model.mesh.axes_file == fs::path("axes.dat") &&
          model.mesh.coast_runtime_faces[0U].size() == 3U,
      "COAST runtime axes source identity survives case wire");
  passed &= expect(
      status && model.mesh.coast_runtime_faces[0U][1U] ==
                    static_cast<double>(static_cast<float>(0.1000000001)),
      "COAST runtime axes wire carries float32-effective faces");

  constexpr std::string_view projected_mesh = R"json({
    "kind":"coast_runtime_axes_v1","axes_file":"axes.dat",
    "domain":{"lower":[0,0,0],"upper":[1,1,0.06]},
    "exact_cells":[2,2,1],"base_spacing":[0.5,0.5,0.06],
    "minimum_spacing":[0.01,0.01,0.01],"max_growth_ratio":1.2,
    "focus_regions":[{"lower":[0.1,0.1,0],"upper":[0.5,0.5,0.06],
                      "target_spacing":[0.1,0.1,0.03]}],
    "limits":{"max_global_cells":16,
              "max_memory_bytes_per_rank":1048576},
    "data_files":[],"immersed_boundary":null
  })json";
  constexpr std::string_view projected_axes = R"data(
COAST_RUNTIME_AXES 1
grid 2 2 1
x 3
0 0.1 1
y 3
0 0.4 1
z 2
0 0.06
)data";
  ScratchCase projected("coast-runtime-axes-projected-boundary");
  projected.write("case.json", case_json(projected_mesh));
  projected.write("thermophysics.d", kPlaceholderThermophysics);
  projected.write("axes.dat", projected_axes);
  ValidatedModel projected_model;
  const Status projected_status = compile(projected.root(), projected_model);
  passed &= expect(
      projected_status && projected_model.mesh.focus_regions.size() == 1U &&
          projected_model.mesh.focus_regions[0U].upper.z ==
              static_cast<double>(static_cast<float>(0.06)),
      "COAST float32 domain projection clips boundary-touching focus regions");
  return passed;
}

bool test_coast_native_air_requires_coast_axes_wire() {
  ScratchCase axes("coast-native-air-axes");
  axes.write("case.json", case_json(kCoastRuntimeAxesMesh));
  axes.write("thermophysics.d", kCoastNativeAirThermophysics);
  axes.write("axes.dat", kCoastRuntimeAxes);
  ValidatedModel native_model;
  bool passed = expect(
      compile(axes.root(), native_model) &&
          native_model.mesh.kind == GeometryKind::coast_runtime_axes_v1 &&
          native_model.thermophysics.species.size() == 1U &&
          native_model.thermophysics.species.front().transport_law ==
              TransportLaw::coast_native_air,
      "COAST-native air compiles through the COAST axes wire family");
  if (!passed) return false;

  ScratchCase non_axes("coast-native-air-non-axes");
  non_axes.write("case.json", case_json(kUniformMesh));
  non_axes.write("thermophysics.d", kCoastNativeAirThermophysics);
  ValidatedModel rejected;
  rejected.fingerprint = UINT64_C(987654321);
  passed &= expect(
      compile(non_axes.root(), rejected).code == StatusCode::invalid_case &&
          rejected.fingerprint == UINT64_C(987654321),
      "non-axes cases reject COAST-native air transactionally");

  ScratchCase ordinary("coast-native-air-serialize-gate");
  ordinary.write("case.json", case_json(kUniformMesh));
  ordinary.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel ordinary_model;
  passed &= expect(static_cast<bool>(compile(ordinary.root(), ordinary_model)),
                   "native-air serialization gate seed compiles");
  ordinary_model.thermophysics = native_model.thermophysics;
  std::vector<std::uint8_t> rejected_payload{0xa5U};
  passed &= expect(
      hundun::v04::detail::serialize_model_for_test(ordinary_model,
                                                    rejected_payload)
                  .code == StatusCode::invalid_case &&
          rejected_payload == std::vector<std::uint8_t>{0xa5U},
      "non-axes serialization rejects COAST-native air transactionally");

  ValidatedModel oversized = native_model;
  constexpr std::int32_t oversized_x_cells = 32768;
  oversized.mesh.exact_cells.x = oversized_x_cells;
  oversized.mesh.base_spacing.x =
      1.0 / static_cast<double>(oversized_x_cells);
  oversized.mesh.minimum_spacing.x = oversized.mesh.base_spacing.x;
  oversized.mesh.limits.max_global_cells =
      static_cast<std::uint64_t>(oversized_x_cells) * 2U;
  oversized.mesh.limits.max_memory_bytes_per_rank = UINT64_C(1073741824);
  oversized.mesh.coast_runtime_faces[0U].resize(
      static_cast<std::size_t>(oversized_x_cells) + 1U);
  for (std::int32_t index = 0; index <= oversized_x_cells; ++index) {
    oversized.mesh.coast_runtime_faces[0U][static_cast<std::size_t>(index)] =
        static_cast<double>(static_cast<float>(index) /
                            static_cast<float>(oversized_x_cells));
  }
  std::vector<std::uint8_t> oversized_payload{0xa5U};
  passed &= expect(
      hundun::v04::detail::serialize_model_for_test(oversized,
                                                    oversized_payload)
                  .code == StatusCode::invalid_case &&
          oversized_payload == std::vector<std::uint8_t>{0xa5U},
      "oversized COAST axes wire rejection preserves the caller payload");

  std::string sutherland{kCoastNativeAirThermophysics};
  passed &= expect(
      replace_once(sutherland, "transport_coast_native_air",
                   "transport_sutherland 1.23456789e-5 321.25 111.75 0.73"),
      "legacy wire seed selects a valid distinct transport law");
  ScratchCase seed("coast-native-air-legacy-wire-seed");
  seed.write("case.json", case_json(kUniformMesh));
  seed.write("thermophysics.d", sutherland);
  ValidatedModel seed_model;
  passed &= expect(static_cast<bool>(compile(seed.root(), seed_model)),
                   "legacy wire seed case compiles");
  if (!passed) return false;
  const double canonical_high_a6 =
      seed_model.thermophysics.species.front().nasa7_high[5U];
  const auto rejects_wire = [&](const std::vector<std::uint8_t>& payload,
                                std::string_view label) {
    ValidatedModel decoded;
    decoded.fingerprint = UINT64_C(987654321);
    const Status status =
        hundun::v04::detail::deserialize_model_for_test(payload, decoded);
    return expect(status.code == StatusCode::invalid_case &&
                      decoded.fingerprint == UINT64_C(987654321),
                  label);
  };
  struct WireVariant {
    CouplingKind coupling;
    LinearAlgorithm algorithm;
    std::uint8_t version;
  };
  constexpr std::array<WireVariant, 4U> variants{{
      {CouplingKind::piso, LinearAlgorithm::fgmres, 11U},
      {CouplingKind::piso, LinearAlgorithm::bicgstab, 12U},
      {CouplingKind::simple, LinearAlgorithm::fgmres, 13U},
      {CouplingKind::simple, LinearAlgorithm::bicgstab, 14U},
  }};
  for (const WireVariant& variant : variants) {
    ValidatedModel wire_model = seed_model;
    wire_model.solver.coupling = variant.coupling;
    wire_model.solver.pressure.algorithm = variant.algorithm;
    if (variant.algorithm == LinearAlgorithm::bicgstab) {
      wire_model.solver.pressure.mg_correction_scaling =
          MgCorrectionScaling::unit_linear;
      wire_model.solver.pressure.krylov_restart = 0U;
    }
    std::vector<std::uint8_t> payload;
    const bool forged =
        hundun::v04::detail::serialize_model_for_test(wire_model, payload) &&
        !payload.empty() && payload.front() == variant.version &&
        forge_legacy_native_air_wire(payload, canonical_high_a6);
    passed &= expect(forged, "legacy native-air wire fixture is well formed");
    if (!forged) continue;
    passed &= rejects_wire(payload,
                           "wire v11-v14 reject COAST-native air");
    if (variant.version <= 12U) {
      std::vector<std::uint8_t> legacy = payload;
      if (legacy.size() <= 9U) {
        passed &= expect(false, "legacy wire fixture has a complete tail");
        continue;
      }
      legacy.erase(legacy.end() - 9);
      legacy.front() = static_cast<std::uint8_t>(variant.version - 2U);
      passed &= rejects_wire(legacy, "wire v9-v10 reject COAST-native air");
    }
  }
  return passed;
}

bool coast_axes_rejects(std::string_view label, std::string_view axes,
                        std::string_view mesh = kCoastRuntimeAxesMesh) {
  ScratchCase scratch(label);
  scratch.write("case.json", case_json(mesh));
  scratch.write("thermophysics.d", kPlaceholderThermophysics);
  scratch.write("axes.dat", axes);
  ValidatedModel model;
  model.fingerprint = UINT64_C(987654321);
  const Status status = compile(scratch.root(), model);
  return expect(status.code == StatusCode::invalid_case &&
                    model.fingerprint == UINT64_C(987654321),
                label);
}

bool test_coast_runtime_axes_strictness_and_fingerprint() {
  bool passed = true;
  passed &= coast_axes_rejects(
      "COAST axes require exact header",
      "NOT_COAST_RUNTIME_AXES 1 grid 2 2 1 x 3 0 .1 1 y 3 0 .4 1 z 2 0 1");
  passed &= coast_axes_rejects(
      "COAST axes require version one",
      "COAST_RUNTIME_AXES 2 grid 2 2 1 x 3 0 .1 1 y 3 0 .4 1 z 2 0 1");
  passed &= coast_axes_rejects(
      "COAST grid count matches declared topology",
      "COAST_RUNTIME_AXES 1 grid 3 2 1 x 4 0 .1 .2 1 y 3 0 .4 1 z 2 0 1");
  passed &= coast_axes_rejects(
      "COAST axis count matches grid count",
      "COAST_RUNTIME_AXES 1 grid 2 2 1 x 4 0 .1 .2 1 y 3 0 .4 1 z 2 0 1");
  passed &= coast_axes_rejects(
      "COAST axes reject non-finite coordinates",
      "COAST_RUNTIME_AXES 1 grid 2 2 1 x 3 0 nan 1 y 3 0 .4 1 z 2 0 1");
  passed &= coast_axes_rejects(
      "COAST axes require strict monotonicity",
      "COAST_RUNTIME_AXES 1 grid 2 2 1 x 3 0 .4 .3 y 3 0 .4 1 z 2 0 1");
  passed &= coast_axes_rejects(
      "COAST axes reject float32 projection collapse",
      "COAST_RUNTIME_AXES 1 grid 2 2 1 x 3 0 1e-50 1 y 3 0 .4 1 z 2 0 1");
  passed &= coast_axes_rejects(
      "COAST axes reject trailing tokens",
      "COAST_RUNTIME_AXES 1 grid 2 2 1 x 3 0 .1 1 y 3 0 .4 1 z 2 0 1 extra");

  std::string memory_limited{kCoastRuntimeAxesMesh};
  passed &= expect(replace_once(memory_limited,
                                "\"max_memory_bytes_per_rank\":1048576",
                                "\"max_memory_bytes_per_rank\":63"),
                   "COAST axes memory-limit fixture mutation");
  passed &= coast_axes_rejects("COAST axes obey mesh memory limit",
                               kCoastRuntimeAxes, memory_limited);

  std::string nested_path{kCoastRuntimeAxesMesh};
  passed &= expect(replace_once(nested_path, "\"axes.dat\"",
                                "\"runtime/axes.dat\""),
                   "COAST nested axes fixture mutation");
  passed &= coast_axes_rejects("COAST axes path is a direct leaf",
                               kCoastRuntimeAxes, nested_path);

  std::string missing_key{kCoastRuntimeAxesMesh};
  passed &= expect(replace_once(missing_key, "\"axes_file\":\"axes.dat\",",
                                ""),
                   "COAST missing axes key fixture mutation");
  passed &= coast_axes_rejects("COAST axes key is mandatory",
                               kCoastRuntimeAxes, missing_key);

  ScratchCase first_case("coast-axes-fingerprint-first");
  first_case.write("case.json", case_json(kCoastRuntimeAxesMesh));
  first_case.write("thermophysics.d", kPlaceholderThermophysics);
  first_case.write("axes.dat", kCoastRuntimeAxes);
  ValidatedModel first;
  passed &= expect(static_cast<bool>(compile(first_case.root(), first)),
                   "first COAST axes fingerprint fixture compiles");

  std::string changed_faces{kCoastRuntimeAxes};
  passed &= expect(replace_once(changed_faces, "0.1000000001", "0.2"),
                   "COAST effective-face fingerprint fixture mutation");
  ScratchCase second_case("coast-axes-fingerprint-second");
  second_case.write("case.json", case_json(kCoastRuntimeAxesMesh));
  second_case.write("thermophysics.d", kPlaceholderThermophysics);
  second_case.write("axes.dat", changed_faces);
  ValidatedModel second;
  passed &= expect(static_cast<bool>(compile(second_case.root(), second)) &&
                       first.fingerprint != second.fingerprint,
                   "effective COAST faces affect case fingerprint");

  std::string reformatted{kCoastRuntimeAxes};
  reformatted.append("\n");
  ScratchCase third_case("coast-axes-fingerprint-third");
  third_case.write("case.json", case_json(kCoastRuntimeAxesMesh));
  third_case.write("thermophysics.d", kPlaceholderThermophysics);
  third_case.write("axes.dat", reformatted);
  ValidatedModel third;
  passed &= expect(static_cast<bool>(compile(third_case.root(), third)) &&
                       first.mesh.coast_runtime_faces ==
                           third.mesh.coast_runtime_faces &&
                       first.fingerprint != third.fingerprint,
                   "COAST source file bytes affect case fingerprint");
  return passed;
}

bool test_every_typed_mesh_field_affects_fingerprint() {
  constexpr std::string_view baseline = R"json({
    "kind":"tensor_stretched",
    "domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,
    "focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],
                      "target_spacing":[0.1,0.1,0.1]}],
    "limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},
    "data_files":[],"immersed_boundary":null
  })json";
  constexpr std::string_view variants[] = {
      R"json({"kind":"tensor_stretched","domain":{"lower":[-0.1,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1.1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":null,"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,11,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.21,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.04,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.15,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.11,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.41,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.11,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10001,"max_memory_bytes_per_rank":67108864},"data_files":[],"immersed_boundary":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108865},"data_files":[],"immersed_boundary":null})json",
  };

  std::uint64_t baseline_fingerprint = 0U;
  bool passed = compile_fingerprint("fingerprint-baseline", baseline,
                                    baseline_fingerprint);
  for (std::size_t index = 0U; index < std::size(variants); ++index) {
    std::uint64_t variant_fingerprint = 0U;
    const std::string label = "fingerprint-variant-" + std::to_string(index);
    passed &= compile_fingerprint(label, variants[index], variant_fingerprint);
    passed &= expect(variant_fingerprint != baseline_fingerprint,
                     "each typed mesh field contributes to fingerprint");
  }
  return passed;
}

bool test_typed_case_fields_and_fingerprint() {
  const std::string baseline = case_json(kUniformMesh);
  std::uint64_t baseline_fingerprint = 0U;
  bool passed = compile_json_fingerprint(
      "typed-fingerprint-baseline", baseline, baseline_fingerprint);
  const auto fingerprint_variant = [&](std::string_view label,
                                       std::string_view from,
                                       std::string_view to) {
    std::string variant = baseline;
    bool local = expect(replace_once(variant, from, to), label);
    std::uint64_t fingerprint = 0U;
    local &= compile_json_fingerprint(label, variant, fingerprint);
    local &= expect(fingerprint != baseline_fingerprint, label);
    return local;
  };
  passed &= fingerprint_variant(
      "transported scalar catalog affects fingerprint",
      "\"transported_scalars\":[]",
      "\"transported_scalars\":[{\"stable_name\":\"mixture_fraction\",\"role\":\"passive_scalar\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9}]");
  passed &= fingerprint_variant(
      "pressure reference affects fingerprint",
      "\"pressure_reference\":\"boundary_absolute\"",
      "\"pressure_reference\":\"closed_mass\"");
  passed &= fingerprint_variant(
      "boundary vector affects fingerprint", "\"velocity\":[1,0,0]",
      "\"velocity\":[2,0,0]");
  passed &= fingerprint_variant(
      "boundary enum affects fingerprint", "\"flow_kind\":\"velocity_inlet\"",
      "\"flow_kind\":\"static_state_inlet\"");
  passed &= fingerprint_variant(
      "boundary scalar data affects fingerprint", "\"scalars\":[]",
      "\"scalars\":[{\"stable_name\":\"mixture_fraction\",\"kind\":\"dirichlet\",\"value\":1,\"backflow_kind\":\"zero_gradient\",\"backflow_value\":0}]");
  passed &= fingerprint_variant(
      "boundary backflow flag affects fingerprint", "\"allow_backflow\":false",
      "\"allow_backflow\":true");
  passed &= fingerprint_variant(
      "scheme affects fingerprint", "\"enthalpy\":\"limited_central2\"",
      "\"enthalpy\":\"central2\"");
  passed &= fingerprint_variant(
      "scheme limiter affects fingerprint", "\"limiter\":1.0",
      "\"limiter\":0.5");
  passed &= fingerprint_variant(
      "time scheme affects fingerprint", "\"scheme\":\"backward_euler\"",
      "\"scheme\":\"variable_bdf2\"");
  passed &= fingerprint_variant(
      "time retry count affects fingerprint", "\"maximum_retries\":8",
      "\"maximum_retries\":9");

  ScratchCase typed("typed-publication");
  std::string typed_json = baseline;
  passed &= expect(replace_once(
                       typed_json, "\"scalars\":[]",
                       "\"scalars\":[{\"stable_name\":\"mixture_fraction\",\"kind\":\"normal_flux\",\"value\":0.125,\"backflow_kind\":\"zero_gradient\",\"backflow_value\":0}]"),
                   "typed scalar fixture mutation");
  typed.write("case.json", typed_json);
  typed.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel model;
  passed &= expect(static_cast<bool>(compile(typed.root(), model)) &&
                       model.boundaries[0].velocity.x == 1.0 &&
                       model.boundaries[0].scalars.size() == 1U &&
                       model.boundaries[0].scalars[0].stable_name ==
                           "mixture_fraction" &&
                       model.boundaries[1].allow_backflow &&
                       model.schemes.momentum ==
                           hundun::v04::ConvectionScheme::limited_central2 &&
                       model.time.initial_dt == 0.001 &&
                       model.time.maximum_retries == 8U,
                   "all typed case sections are published");
  return passed;
}

bool test_transported_scalar_catalog() {
  const std::string baseline = case_json(kUniformMesh);
  std::string valid = baseline;
  bool passed = expect(
      replace_once(
          valid, "\"transported_scalars\":[]",
          "\"transported_scalars\":[{\"stable_name\":\"O2\",\"role\":\"species\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9},{\"stable_name\":\"mixture_fraction\",\"role\":\"passive_scalar\",\"molecular_schmidt\":0.8,\"turbulent_schmidt\":0.6}]"),
      "catalog valid fixture mutation");
  ScratchCase valid_case("transported-scalars-valid");
  valid_case.write("case.json", valid);
  valid_case.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel model;
  passed &= expect(
      static_cast<bool>(compile(valid_case.root(), model)) &&
          model.transported_scalars.size() == 2U &&
          model.transported_scalars[0].stable_name == "O2" &&
          model.transported_scalars[0].role == TransportedScalarRole::species &&
          model.transported_scalars[0].molecular_schmidt == 0.7 &&
          model.transported_scalars[0].turbulent_schmidt == 0.9 &&
          model.transported_scalars[1].stable_name == "mixture_fraction" &&
          model.transported_scalars[1].role ==
              TransportedScalarRole::passive_scalar &&
          model.transported_scalars[1].molecular_schmidt == 0.8 &&
          model.transported_scalars[1].turbulent_schmidt == 0.6,
      "catalog is retained in declared stable order");

  const auto reject_catalog = [&](std::string_view label,
                                  std::string_view catalog) {
    std::string json = baseline;
    return expect(replace_once(json, "\"transported_scalars\":[]", catalog),
                  label) &&
           rejects(label, json);
  };
  std::string missing = baseline;
  passed &= expect(replace_once(missing, "\"transported_scalars\":[],", ""),
                   "missing transported scalar catalog mutation") &&
            rejects("missing transported scalar catalog", missing);
  passed &= reject_catalog("transported scalar catalog must be an array",
                            "\"transported_scalars\":{}");
  passed &= reject_catalog(
      "duplicate transported scalar names",
      "\"transported_scalars\":[{\"stable_name\":\"z\",\"role\":\"species\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9},{\"stable_name\":\"z\",\"role\":\"passive_scalar\",\"molecular_schmidt\":0.8,\"turbulent_schmidt\":0.6}]");
  passed &= reject_catalog(
      "reserved transported scalar U",
      "\"transported_scalars\":[{\"stable_name\":\"U\",\"role\":\"species\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9}]");
  passed &= reject_catalog(
      "reserved transported scalar pi",
      "\"transported_scalars\":[{\"stable_name\":\"pi\",\"role\":\"species\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9}]");
  passed &= reject_catalog(
      "reserved transported scalar h",
      "\"transported_scalars\":[{\"stable_name\":\"h\",\"role\":\"passive_scalar\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9}]");
  passed &= reject_catalog(
      "unknown transported scalar role",
      "\"transported_scalars\":[{\"stable_name\":\"z\",\"role\":\"temperature\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9}]");
  passed &= reject_catalog(
      "unknown transported scalar key",
      "\"transported_scalars\":[{\"stable_name\":\"z\",\"role\":\"species\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9,\"units\":\"1\"}]");
  passed &= reject_catalog(
      "invalid transported scalar name",
      "\"transported_scalars\":[{\"stable_name\":\"bad/name\",\"role\":\"species\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9}]");

  std::string oversized = "\"transported_scalars\":[";
  for (std::size_t index = 0U; index < 65U; ++index) {
    if (index != 0U) {
      oversized += ',';
    }
    oversized += "{\"stable_name\":\"s_" + std::to_string(index) +
                 "\",\"role\":\"species\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9}";
  }
  oversized += ']';
  passed &= reject_catalog("transported scalar catalog limit", oversized);

  std::string species = baseline;
  std::string passive = baseline;
  passed &= expect(
      replace_once(species, "\"transported_scalars\":[]",
                   "\"transported_scalars\":[{\"stable_name\":\"z\",\"role\":\"species\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9}]"),
      "species fingerprint fixture mutation");
  passed &= expect(
      replace_once(passive, "\"transported_scalars\":[]",
                   "\"transported_scalars\":[{\"stable_name\":\"z\",\"role\":\"passive_scalar\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9}]"),
      "passive scalar fingerprint fixture mutation");
  std::uint64_t species_fingerprint = 0U;
  std::uint64_t passive_fingerprint = 0U;
  passed &= compile_json_fingerprint("species catalog fingerprint", species,
                                     species_fingerprint);
  passed &= compile_json_fingerprint("passive catalog fingerprint", passive,
                                     passive_fingerprint);
  passed &= expect(species_fingerprint != passive_fingerprint,
                   "transported scalar role affects fingerprint");

  std::string transport_closed = baseline;
  passed &= expect(
      replace_once(
          transport_closed, "\"transported_scalars\":[]",
          "\"transported_scalars\":[{\"stable_name\":\"z\",\"role\":\"passive_scalar\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0.9}]"),
      "transport closure fixture mutation");
  ScratchCase transport_case("transported-scalar-closure");
  transport_case.write("case.json", transport_closed);
  transport_case.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel transport_model;
  passed &= expect(
      static_cast<bool>(compile(transport_case.root(), transport_model)) &&
          transport_model.transported_scalars.size() == 1U &&
          transport_model.transported_scalars[0].molecular_schmidt == 0.7 &&
          transport_model.transported_scalars[0].turbulent_schmidt == 0.9,
      "transported scalar publishes both molecular and turbulent closures");
  passed &= reject_catalog(
      "nonpositive molecular Schmidt number",
      "\"transported_scalars\":[{\"stable_name\":\"z\",\"role\":\"passive_scalar\",\"molecular_schmidt\":0,\"turbulent_schmidt\":0.9}]");
  passed &= reject_catalog(
      "nonpositive turbulent Schmidt number",
      "\"transported_scalars\":[{\"stable_name\":\"z\",\"role\":\"passive_scalar\",\"molecular_schmidt\":0.7,\"turbulent_schmidt\":0}]");
  return passed;
}

bool test_typed_case_rejections() {
  const std::string baseline = case_json(kUniformMesh);
  const auto reject_variant = [&](std::string_view label,
                                  std::string_view from,
                                  std::string_view to) {
    std::string variant = baseline;
    return expect(replace_once(variant, from, to), label) &&
           rejects(label, variant);
  };
  bool passed = true;
  passed &= reject_variant(
      "missing boundary field", "\"flow_kind\":\"velocity_inlet\",", "");
  passed &= reject_variant(
      "unknown boundary field", "\"flow_kind\":\"velocity_inlet\"",
      "\"flow_kind\":\"velocity_inlet\",\"unknown\":0");
  passed &= reject_variant(
      "unknown boundary enum", "\"flow_kind\":\"velocity_inlet\"",
      "\"flow_kind\":\"farfield\"");
  passed &= reject_variant(
      "invalid boundary Mach limit", "\"mach_limit\":0.95",
      "\"mach_limit\":1.0");
  passed &= reject_variant(
      "negative boundary relaxation", "\"relaxation\":1",
      "\"relaxation\":-0.1");
  passed &= reject_variant(
      "duplicate scalar names", "\"scalars\":[]",
      "\"scalars\":[{\"stable_name\":\"z\",\"kind\":\"dirichlet\",\"value\":0,\"backflow_kind\":\"zero_gradient\",\"backflow_value\":0},{\"stable_name\":\"z\",\"kind\":\"zero_gradient\",\"value\":1,\"backflow_kind\":\"zero_gradient\",\"backflow_value\":0}]");
  passed &= reject_variant(
      "invalid scalar stable name", "\"scalars\":[]",
      "\"scalars\":[{\"stable_name\":\"bad/name\",\"kind\":\"dirichlet\",\"value\":0,\"backflow_kind\":\"zero_gradient\",\"backflow_value\":0}]");
  passed &= reject_variant(
      "unknown scalar enum", "\"scalars\":[]",
      "\"scalars\":[{\"stable_name\":\"z\",\"kind\":\"fixed_value\",\"value\":0,\"backflow_kind\":\"zero_gradient\",\"backflow_value\":0}]");
  passed &= reject_variant(
      "missing scheme field", "\"diffusion\":\"central2\",", "");
  passed &= reject_variant(
      "unknown convection scheme", "\"momentum\":\"limited_central2\"",
      "\"momentum\":\"upwind1\"");
  passed &= reject_variant(
      "scheme limiter out of range", "\"limiter\":1.0",
      "\"limiter\":1.1");
  passed &= reject_variant(
      "missing time field", "\"maximum_retries\":8,", "");
  passed &= reject_variant(
      "unknown time scheme", "\"scheme\":\"backward_euler\"",
      "\"scheme\":\"crank_nicolson\"");
  passed &= reject_variant(
      "time step outside limits", "\"initial_dt\":0.001",
      "\"initial_dt\":1.0");
  passed &= reject_variant(
      "zero retry count", "\"maximum_retries\":8",
      "\"maximum_retries\":0");
  passed &= reject_variant(
      "invalid BDF ratio interval", "\"minimum_bdf_ratio\":0.2",
      "\"minimum_bdf_ratio\":6.0");
  return passed;
}

bool test_thermophysics_schema() {
  const std::string baseline = case_json(kUniformMesh);
  bool passed = true;

  std::string missing_root = baseline;
  passed &= expect(
      replace_once(
          missing_root,
          ",\"thermophysics\":{\"data_file\":\"thermophysics.d\"}", ""),
      "missing thermophysics root fixture mutation") &&
            rejects("missing thermophysics root key", missing_root);

  std::string missing_data_file = baseline;
  passed &= expect(
      replace_once(
          missing_data_file,
          "\"thermophysics\":{\"data_file\":\"thermophysics.d\"}",
          "\"thermophysics\":{}"),
      "missing thermophysics data_file fixture mutation") &&
            rejects("missing thermophysics data_file key", missing_data_file);

  std::string unknown = baseline;
  passed &= expect(
      replace_once(
          unknown,
          "\"thermophysics\":{\"data_file\":\"thermophysics.d\"}",
          "\"thermophysics\":{\"data_file\":\"thermophysics.d\",\"units\":\"SI\"}"),
      "unknown thermophysics key fixture mutation") &&
            rejects("unknown thermophysics key", unknown);

  std::string wrong_suffix = baseline;
  passed &= expect(replace_once(wrong_suffix, "\"thermophysics.d\"",
                                "\"thermophysics.txt\""),
                   "wrong thermophysics suffix fixture mutation") &&
            rejects("wrong thermophysics suffix", wrong_suffix);
  return passed;
}

bool test_mesh_rejections() {
  bool passed = true;
  passed &= rejects("unknown mesh key", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null,"extra":0})json"));
  passed &= rejects("unordered domain", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,1],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("non-finite mesh value", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1e999,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("non-finite derived domain span", case_json(R"json({
    "kind":"uniform","domain":{"lower":[-1.7e308,0,0],"upper":[1.7e308,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("uniform requires exact cells", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("uniform rejects base spacing", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":[0.25,0.25,0.25],
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("uniform minimum exceeds width", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.3,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("uniform requires unit growth", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1.1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("exact cell limit", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":63,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("tensor requires base spacing", case_json(R"json({
    "kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":null,
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.2,
    "focus_regions":[],"limits":{"max_global_cells":100,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("focus target below minimum", case_json(R"json({
    "kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":[0.2,0.2,0.2],
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.2,
    "focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.04,0.1,0.1]}],
    "limits":{"max_global_cells":1000,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("focus target exceeds base", case_json(R"json({
    "kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":[0.2,0.2,0.2],
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.2,
    "focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.21,0.1,0.1]}],
    "limits":{"max_global_cells":1000,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  passed &= rejects("focus wholly outside domain", case_json(R"json({
    "kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":[0.2,0.2,0.2],
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.2,
    "focus_regions":[{"lower":[2,2,2],"upper":[3,3,3],"target_spacing":[0.1,0.1,0.1]}],
    "limits":{"max_global_cells":1000,"max_memory_bytes_per_rank":1},
    "data_files":[],"immersed_boundary":null})json"));
  return passed;
}

bool test_case_and_reference_security() {
  bool passed = true;
  {
    ScratchCase missing("missing-json");
    ValidatedModel model;
    passed &= expect(compile(missing.root(), model).code ==
                         StatusCode::invalid_case,
                     "missing case.json is rejected");
  }
  passed &= rejects("recursive duplicate keys",
                    case_json(kUniformMesh).replace(
                        case_json(kUniformMesh).find("\"model\":"),
                        std::string{"\"model\":"}.size(),
                        "\"model\":\"single_phase_low_mach_compressible\",\"model\":"));
  std::string unknown_units = case_json(kUniformMesh);
  unknown_units.replace(unknown_units.find("\"SI\""), 4U, "\"CGS\"");
  passed &= rejects("unknown units", unknown_units);
  std::string unknown_top = case_json(kUniformMesh);
  unknown_top.insert(unknown_top.find("\"mesh\""), "\"extra\":0,");
  passed &= rejects("unknown top-level key",
                    unknown_top);
  passed &= rejects("constant density flow", case_json(
      kUniformMesh,
      R"json({"model":"constant_density","pressure_reference":"boundary_absolute","reacting":false})json"));
  passed &= rejects("missing pressure reference", case_json(
      kUniformMesh,
      R"json({"model":"single_phase_low_mach_compressible","reacting":false})json"));
  passed &= rejects("reacting flow", case_json(
      kUniformMesh,
      R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"boundary_absolute","reacting":true})json"));
  passed &= rejects("PIMPLE coupling", case_json(
      kUniformMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"boundary_absolute","reacting":false})json",
      R"json({"coupling":"PIMPLE","pressure_correctors":2})json"));
  passed &= rejects("wrong corrector count", case_json(
      kUniformMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"boundary_absolute","reacting":false})json",
      R"json({"coupling":"PISO","pressure_correctors":3})json"));
  passed &= rejects("zero pressure tolerance", case_json(
      kUniformMesh,
      R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"boundary_absolute","reacting":false})json",
      R"json({"coupling":"PISO","pressure_correctors":2,"pressure_linear":{"absolute_tolerance":0,"relative_tolerance":1e-6,"maximum_iterations":500,"true_residual_interval":4,"krylov_restart":12},"terminal_tolerances":{"eos":1e-6,"continuity":1e-6,"closed_mass":1e-6,"gauge":1e-6}})json"));
  passed &= rejects("pressure tolerance looser than terminal gate", case_json(
      kUniformMesh,
      R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"boundary_absolute","reacting":false})json",
      R"json({"coupling":"PISO","pressure_correctors":2,"pressure_linear":{"absolute_tolerance":1e-8,"relative_tolerance":1e-4,"maximum_iterations":500,"true_residual_interval":4,"krylov_restart":12},"terminal_tolerances":{"eos":1e-6,"continuity":1e-6,"closed_mass":1e-6,"gauge":1e-6}})json"));
  passed &= rejects("invalid Krylov restart", case_json(
      kUniformMesh,
      R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"boundary_absolute","reacting":false})json",
      R"json({"coupling":"PISO","pressure_correctors":2,"pressure_linear":{"absolute_tolerance":1e-8,"relative_tolerance":1e-6,"maximum_iterations":500,"true_residual_interval":4,"krylov_restart":65},"terminal_tolerances":{"eos":1e-6,"continuity":1e-6,"closed_mass":1e-6,"gauge":1e-6}})json"));
  passed &= rejects("unsupported turbulence", case_json(
      kUniformMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"boundary_absolute","reacting":false})json",
      R"json({"coupling":"PISO","pressure_correctors":2})json",
      R"json({"model":"smagorinsky"})json"));

  const auto reference_mesh = [](std::string_view data,
                                 std::string_view immersed) {
    return std::string{R"json({"kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[4,4,4],"base_spacing":null,"minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,"focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},"data_files":)json"} +
           std::string(data) + ",\"immersed_boundary\":" +
           std::string(immersed) + "}";
  };
  passed &= rejects("parent path", case_json(reference_mesh("[\"../table.d\"]", "null")));
  passed &= rejects("absolute data path", case_json(reference_mesh("[\"/tmp/table.d\"]", "null")));
  passed &= rejects("nested data path", case_json(reference_mesh("[\"data/table.d\"]", "null")));
  passed &= rejects(
      "nested STL",
      case_json(reference_mesh(
          "[]",
          R"json({"stl_file":"geometry/body.stl","fluid_side":"outside"})json")));
  passed &= rejects(
      "invalid immersed fluid side",
      case_json(reference_mesh(
          "[]",
          R"json({"stl_file":"body.stl","fluid_side":"middle"})json")));
  passed &= rejects(
      "immersed boundary requires fluid side",
      case_json(reference_mesh("[]", R"json({"stl_file":"body.stl"})json")));
  passed &= rejects(
      "immersed boundary rejects unknown key",
      case_json(reference_mesh(
          "[]",
          R"json({"stl_file":"body.stl","fluid_side":"outside","extra":0})json")));
  passed &= rejects("wrong data suffix", case_json(reference_mesh("[\"table.txt\"]", "null")));
  passed &= rejects("missing data file", case_json(reference_mesh("[\"missing.d\"]", "null")));

  ScratchCase outside("outside-symlink");
  ScratchCase inside("inside-symlink");
  outside.write("profile.d", "outside bytes\n");
  std::error_code error;
  fs::create_symlink(outside.root() / "profile.d",
                     inside.root() / "profile.d", error);
  inside.write("case.json", case_json(reference_mesh("[\"profile.d\"]", "null")));
  inside.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel model;
  passed &= expect(!error && compile(inside.root(), model).code ==
                                StatusCode::invalid_case,
                   "referenced-file symlink escape is rejected");

  ScratchCase outside_json("outside-json");
  ScratchCase inside_json("inside-json");
  outside_json.write("case.json", case_json(kUniformMesh));
  outside_json.write("thermophysics.d", kPlaceholderThermophysics);
  error.clear();
  fs::create_symlink(outside_json.root() / "case.json",
                     inside_json.root() / "case.json", error);
  passed &= expect(!error && compile(inside_json.root(), model).code ==
                                StatusCode::invalid_case,
                   "case.json symlink escape is rejected");

  ScratchCase swap_outside("swap-outside");
  ScratchCase swap_inside("swap-inside");
  swap_outside.write("escaped.d", "outside bytes\n");
  swap_inside.write("profile.d", "inside bytes\n");
  swap_inside.write("case.json", case_json(reference_mesh(
                                   "[\"profile.d\"]", "null")));
  swap_inside.write("thermophysics.d", kPlaceholderThermophysics);
  g_swap_target = swap_inside.root() / "profile.d";
  g_swap_source = swap_outside.root() / "escaped.d";
  g_swap_open_count = 0;
  hundun::v04::detail::set_file_open_observer_for_test(
      swap_reference_before_open);
  const Status swap_status = compile(swap_inside.root(), model);
  hundun::v04::detail::set_file_open_observer_for_test(nullptr);
  passed &= expect(swap_status.code == StatusCode::invalid_case,
                   "reference swap after validation is rejected");

  ScratchCase hardlinks("hardlinks");
  hardlinks.write("first.d", "same bytes\n");
  error.clear();
  fs::create_hard_link(hardlinks.root() / "first.d",
                       hardlinks.root() / "second.d", error);
  hardlinks.write("case.json", case_json(reference_mesh(
                                  "[\"first.d\",\"second.d\"]", "null")));
  hardlinks.write("thermophysics.d", kPlaceholderThermophysics);
  passed &= expect(!error && compile(hardlinks.root(), model).code ==
                                StatusCode::invalid_case,
                   "hard-link aliases are rejected");

  ScratchCase thermophysics_alias("thermophysics-alias");
  thermophysics_alias.write("thermophysics.d", kPlaceholderThermophysics);
  error.clear();
  fs::create_hard_link(thermophysics_alias.root() / "thermophysics.d",
                       thermophysics_alias.root() / "profile.d", error);
  thermophysics_alias.write(
      "case.json", case_json(reference_mesh("[\"profile.d\"]", "null")));
  passed &= expect(
      !error && compile(thermophysics_alias.root(), model).code ==
                    StatusCode::invalid_case,
      "thermophysics and generic-data hard-link aliases are rejected");
  return passed;
}

bool test_wire_rejects_duplicate_data_paths() {
  ScratchCase scratch("wire-duplicate-data-path");
  scratch.write("case.json", case_json(kUniformMesh));
  scratch.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel model;
  if (!expect(static_cast<bool>(compile(scratch.root(), model)),
              "wire mutation fixture compiles")) {
    return false;
  }

  // The encoder must reject duplicate references too, so construct the
  // malformed payload by serializing two distinct same-length names and
  // mutating the second name in-place. This reaches the independent decoder
  // invariant without relying on filesystem parsing.
  model.data_files = {"first.d", "other.d"};
  std::vector<std::uint8_t> payload;
  if (!expect(static_cast<bool>(
                  hundun::v04::detail::serialize_model_for_test(model,
                                                               payload)),
              "wire mutation seed serializes")) {
    return false;
  }
  constexpr std::string_view first{"first.d"};
  constexpr std::string_view other{"other.d"};
  const auto occurrence = [&](std::string_view text) {
    return std::search(payload.begin(), payload.end(), text.begin(),
                       text.end());
  };
  const auto first_path = occurrence(first);
  const auto second_path = occurrence(other);
  if (!expect(first_path != payload.end() && second_path != payload.end(),
              "wire mutation locates both typed path records")) {
    return false;
  }
  std::copy(first.begin(), first.end(), second_path);

  ValidatedModel decoded;
  decoded.fingerprint = UINT64_C(987654321);
  const Status status =
      hundun::v04::detail::deserialize_model_for_test(payload, decoded);
  return expect(status.code == StatusCode::invalid_case &&
                    decoded.fingerprint == UINT64_C(987654321),
                "wire decoder rejects duplicate data paths transactionally");
}

bool test_defaults_and_enums() {
  ScratchCase defaults("defaults");
  defaults.write("case.json", case_json(kUniformMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"boundary_absolute","reacting":false})json",
                                         R"json({"coupling":"PISO","pressure_correctors":2})json", ""));
  defaults.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel default_model;
  bool passed = expect(static_cast<bool>(compile(defaults.root(), default_model)) &&
                           default_model.solver.coupling ==
                               CouplingKind::piso &&
                           default_model.turbulence ==
                               TurbulenceKind::vreman_wall_function,
                       "omitted turbulence selects the production default");
  ScratchCase simple("simple-coupling");
  simple.write(
      "case.json",
      case_json(
          kUniformMesh,
          R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"boundary_absolute","reacting":false})json",
          R"json({"coupling":"SIMPLE","pressure_correctors":2})json"));
  simple.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel simple_model;
  passed &= expect(static_cast<bool>(compile(simple.root(), simple_model)) &&
                       simple_model.solver.coupling == CouplingKind::simple &&
                       simple_model.fingerprint != default_model.fingerprint,
                   "SIMPLE coupling compiles as a distinct typed model");
  ScratchCase enums("enums");
  enums.write("case.json", case_json(kTensorMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"closed_mass","reacting":false})json",
                                      R"json({"coupling":"PISO","pressure_correctors":2})json",
                                      R"json({"model":"wale"})json",
                                      R"json({"control":"adaptive_acoustic","scheme":"variable_bdf2","initial_dt":0.001,"minimum_dt":1e-8,"maximum_dt":0.1,"convective_cfl":0.8,"viscous_cfl":0.5,"thermal_cfl":0.5,"species_cfl":0.5,"acoustic_cfl":0.8,"maximum_growth":1.25,"retry_factor":0.5,"maximum_retries":8,"minimum_bdf_ratio":0.2,"maximum_bdf_ratio":5.0})json"));
  enums.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel enum_model;
  passed &= expect(static_cast<bool>(compile(enums.root(), enum_model)) &&
                       enum_model.turbulence == TurbulenceKind::wale &&
                       enum_model.pressure_reference ==
                           PressureReferenceKind::closed_mass &&
                       enum_model.time.control ==
                           TimeControlKind::adaptive_acoustic &&
                       enum_model.time.scheme == TimeScheme::variable_bdf2,
                   "supported turbulence and time enums compile");
  ScratchCase controls("solver-controls");
  controls.write(
      "case.json",
      case_json(
          kUniformMesh,
          R"json({"model":"single_phase_low_mach_compressible","pressure_reference":"boundary_absolute","reacting":false})json",
          R"json({"coupling":"PISO","pressure_correctors":2,"pressure_linear":{"absolute_tolerance":1e-8,"relative_tolerance":1e-6,"maximum_iterations":500,"true_residual_interval":5,"krylov_restart":16},"terminal_tolerances":{"eos":2e-6,"continuity":3e-6,"closed_mass":4e-6,"gauge":5e-6}})json"));
  controls.write("thermophysics.d", kPlaceholderThermophysics);
  ValidatedModel controlled_model;
  passed &= expect(
      static_cast<bool>(compile(controls.root(), controlled_model)) &&
          controlled_model.solver.pressure.absolute_tolerance == 1.0e-8 &&
          controlled_model.solver.pressure.relative_tolerance == 1.0e-6 &&
          controlled_model.solver.pressure.maximum_iterations == 500U &&
          controlled_model.solver.pressure.true_residual_interval == 5U &&
          controlled_model.solver.pressure.krylov_restart == 16U &&
          controlled_model.solver.terminal.eos == 2.0e-6 &&
          controlled_model.solver.terminal.continuity == 3.0e-6 &&
          controlled_model.solver.terminal.closed_mass == 4.0e-6 &&
          controlled_model.solver.terminal.gauge == 5.0e-6,
      "typed solver controls compile without product hard-coding");
  return passed;
}

bool test_immersed_reconstruction_policy_is_typed_and_hashed() {
  const auto immersed_mesh = [](std::string_view policy) {
    std::string mesh{kUniformMesh};
    const std::string immersed =
        std::string{R"json({"stl_file":"body.stl","fluid_side":"outside","reconstruction_policy":")json"} +
        std::string(policy) + R"json("})json";
    return replace_once(mesh, "\"immersed_boundary\":null",
                        "\"immersed_boundary\":" + immersed)
               ? mesh
               : std::string{};
  };

  ScratchCase strict_case("ibm-strict-policy");
  strict_case.write("case.json",
                    case_json(immersed_mesh("strict_quadratic")));
  strict_case.write("thermophysics.d", kPlaceholderThermophysics);
  strict_case.write("body.stl", "strict policy identity\n");
  ValidatedModel strict_model;
  bool passed = expect(
      static_cast<bool>(compile(strict_case.root(), strict_model)) &&
          strict_model.immersed_boundary.has_value() &&
          strict_model.immersed_boundary->reconstruction_policy ==
              IbmReconstructionPolicy::strict_quadratic,
      "strict_quadratic is a typed immersed reconstruction policy");

  ScratchCase adaptive_case("ibm-adaptive-policy");
  adaptive_case.write("case.json",
                      case_json(immersed_mesh("adaptive_order")));
  adaptive_case.write("thermophysics.d", kPlaceholderThermophysics);
  adaptive_case.write("body.stl", "strict policy identity\n");
  ValidatedModel adaptive_model;
  passed &= expect(
      static_cast<bool>(compile(adaptive_case.root(), adaptive_model)) &&
          adaptive_model.immersed_boundary.has_value() &&
          adaptive_model.immersed_boundary->reconstruction_policy ==
              IbmReconstructionPolicy::adaptive_order &&
          adaptive_model.fingerprint != strict_model.fingerprint,
      "adaptive_order is typed and changes immutable case identity");

  ScratchCase invalid_case("ibm-invalid-policy");
  invalid_case.write("case.json", case_json(immersed_mesh("nearest_copy")));
  invalid_case.write("thermophysics.d", kPlaceholderThermophysics);
  invalid_case.write("body.stl", "strict policy identity\n");
  ValidatedModel invalid_model;
  passed &= expect(compile(invalid_case.root(), invalid_model).code ==
                       StatusCode::invalid_case,
                   "unapproved nearest-copy fallback is rejected");
  return passed;
}

bool test_field_registry() {
  FieldRegistry registry;
  FieldId first = 99;
  FieldId second = 99;
  bool passed = expect(static_cast<bool>(registry.declare_field("alpha", 1, 0, first)) &&
                           static_cast<bool>(registry.declare_field("beta", 3, 2, second)) &&
                           first == 0 && second == 1,
                       "field IDs follow declaration order");
  FieldId ignored = 0;
  passed &= expect(registry.declare_field("alpha", 1, 0, ignored).code ==
                       StatusCode::invalid_plan,
                   "duplicate field is rejected");
  passed &= expect(registry.declare_field("", 1, 0, ignored).code ==
                       StatusCode::invalid_plan,
                   "empty field is rejected");
  passed &= expect(registry.declare_field("bad", 0, 0, ignored).code ==
                       StatusCode::invalid_plan,
                   "zero components are rejected");
  passed &= expect(registry.declare_field("wide", 1, 255, ignored).code ==
                       StatusCode::invalid_plan,
                   "unsupported ghost width is rejected");
  FieldSchema schema;
  passed &= expect(static_cast<bool>(registry.freeze(schema)) &&
                       schema.size() == 2U && schema[0].id == 0U &&
                       schema[0].stable_name == "alpha" &&
                       schema[0].components == 1U &&
                       schema[0].ghost_width == 0U && schema[1].id == 1U &&
                       schema[1].stable_name == "beta" &&
                       schema[1].components == 3U &&
                       schema[1].ghost_width == 2U,
                   "freeze creates an ordered immutable snapshot");
  passed &= expect(registry.declare_field("late", 1, 0, ignored).code ==
                       StatusCode::invalid_plan,
                   "post-freeze mutation is rejected");
  return passed;
}

bool test_field_id_overflow() {
  FieldRegistry registry;
  FieldId id = 0;
  for (std::uint32_t value = 0;
       value <= std::numeric_limits<FieldId>::max(); ++value) {
    const Status status =
        registry.declare_field("field_" + std::to_string(value), 1, 0, id);
    if (!status || id != static_cast<FieldId>(value)) {
      return expect(false, "all representable FieldIds are assigned");
    }
  }
  return expect(registry.declare_field("overflow", 1, 0, id).code ==
                    StatusCode::invalid_plan,
                "FieldId overflow is rejected");
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = true;
  passed &= test_valid_fixture();
  passed &= test_tensor_normalization_and_fingerprint();
  passed &= test_coast_runtime_axes_case();
  passed &= test_coast_native_air_requires_coast_axes_wire();
  passed &= test_coast_runtime_axes_strictness_and_fingerprint();
  passed &= test_every_typed_mesh_field_affects_fingerprint();
  passed &= test_typed_case_fields_and_fingerprint();
  passed &= test_transported_scalar_catalog();
  passed &= test_typed_case_rejections();
  passed &= test_thermophysics_schema();
  passed &= test_mesh_rejections();
  passed &= test_case_and_reference_security();
  passed &= test_wire_rejects_duplicate_data_paths();
  passed &= test_defaults_and_enums();
  passed &= test_immersed_reconstruction_policy_is_typed_and_hashed();
  passed &= test_field_registry();
  passed &= test_field_id_overflow();
  const int finalize_status = MPI_Finalize();
  return passed && finalize_status == MPI_SUCCESS ? 0 : 1;
}
