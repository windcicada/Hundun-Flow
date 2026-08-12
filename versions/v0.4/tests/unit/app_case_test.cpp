// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_case.hpp"
#include "hundun/v04_field.hpp"

#include "app_case_detail.hpp"

#include <mpi.h>

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>

namespace {

namespace fs = std::filesystem;
using hundun::v04::CaseCompiler;
using hundun::v04::FieldId;
using hundun::v04::FieldRegistry;
using hundun::v04::FieldSchema;
using hundun::v04::GeometryKind;
using hundun::v04::Status;
using hundun::v04::StatusCode;
using hundun::v04::TimeControlKind;
using hundun::v04::TurbulenceKind;
using hundun::v04::ValidatedModel;

constexpr std::string_view kValidCase = R"json({
  "schema_version": 1,
  "units": "SI",
  "mesh": {"kind": "uniform", "data_files": [], "stl_file": null},
  "flow": {
    "model": "single_phase_low_mach_compressible",
    "pressure_closure": "local_absolute_pressure_drho_dp",
    "reacting": false
  },
  "solver": {"coupling": "PISO", "pressure_correctors": 2},
  "turbulence": {"model": "none"},
  "time": {"control": "fixed"}
})json";

constexpr std::uint64_t kFingerprintOffset = 14695981039346656037ULL;
constexpr std::uint64_t kFingerprintPrime = 1099511628211ULL;
constexpr std::string_view kSemanticContract =
    "HUNDUN-FLOW-v0.4-case-schema-v1|units=SI|"
    "flow=single_phase_low_mach_compressible|"
    "pressure_closure=local_absolute_pressure_drho_dp|reacting=false|"
    "coupling=PISO|pressure_correctors=2";

class FingerprintOracle {
 public:
  void bytes(const void* data, std::size_t size) {
    const auto* input = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; ++index) {
      value_ ^= static_cast<std::uint64_t>(input[index]);
      value_ *= kFingerprintPrime;
    }
  }

  template <class Integer>
  void integer(Integer value) {
    using Unsigned = std::make_unsigned_t<Integer>;
    const Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t byte = 0; byte < sizeof(Unsigned); ++byte) {
      const auto part = static_cast<unsigned char>(
          (bits >> (byte * 8U)) & static_cast<Unsigned>(0xffU));
      bytes(&part, 1U);
    }
  }

  void text(std::string_view value) {
    integer(static_cast<std::uint64_t>(value.size()));
    bytes(value.data(), value.size());
  }

  std::uint64_t finish() const { return value_ == 0U ? 1U : value_; }

 private:
  std::uint64_t value_{kFingerprintOffset};
};

std::uint64_t empty_model_fingerprint(GeometryKind geometry,
                                      TurbulenceKind turbulence,
                                      TimeControlKind time_control) {
  FingerprintOracle hash;
  hash.text(kSemanticContract);
  hash.integer(static_cast<std::uint8_t>(geometry));
  hash.integer(static_cast<std::uint8_t>(turbulence));
  hash.integer(static_cast<std::uint8_t>(time_control));
  hash.integer(std::uint16_t{0});
  hash.text("no-stl");
  return hash.finish();
}

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

class ScratchCase {
 public:
  explicit ScratchCase(std::string_view label) {
    static std::uint64_t sequence = 0;
    root_ = fs::temp_directory_path() /
            ("hundun-v04-case-" + std::string(label) + "-" +
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

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

Status compile(const fs::path& root, ValidatedModel& model) {
  return CaseCompiler::load_and_compile(MPI_COMM_SELF, root, model);
}

bool rejects(std::string_view label, std::string_view json,
             StatusCode expected = StatusCode::invalid_case) {
  ScratchCase scratch(label);
  scratch.write("case.json", json);
  ValidatedModel model{};
  const Status status = compile(scratch.root(), model);
  return expect(status.code == expected, label);
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
  ValidatedModel model{};
  const Status status = compile(scratch.root(), model);
  bool passed = expect(static_cast<bool>(status), "minimal fixture compiles");
  passed &= expect(model.geometry == GeometryKind::uniform,
                   "fixture selects uniform geometry");
  passed &= expect(model.turbulence == TurbulenceKind::none,
                   "fixture selects no turbulence");
  passed &= expect(model.time_control == TimeControlKind::fixed,
                   "fixture selects fixed time control");
  passed &= expect(model.data_files.empty(), "fixture has no data files");
  passed &= expect(!model.stl_file.has_value(), "fixture has no STL file");
  passed &= expect(model.fingerprint != 0, "fixture has a fingerprint");
  passed &= expect(
      model.fingerprint ==
          empty_model_fingerprint(GeometryKind::uniform,
                                  TurbulenceKind::none,
                                  TimeControlKind::fixed),
      "fingerprint includes the complete fixed numerical contract");
  return passed;
}

bool test_case_json_symlink_escape() {
  ScratchCase outside("outside-json");
  ScratchCase inside("inside-json");
  outside.write("case.json", kValidCase);

  std::error_code error;
  fs::create_symlink(outside.root() / "case.json",
                     inside.root() / "case.json", error);
  if (error) {
    return expect(false, "case.json symlink fixture is created");
  }

  ValidatedModel model{};
  return expect(compile(inside.root(), model).code == StatusCode::invalid_case,
                "case.json symlink escape is rejected");
}

bool test_reference_symlink_escape() {
  ScratchCase outside("outside-reference");
  ScratchCase inside("inside-reference");
  outside.write("profile.d", "0 1\n");

  std::error_code error;
  fs::create_symlink(outside.root() / "profile.d",
                     inside.root() / "profile.d", error);
  if (error) {
    return expect(false, "reference symlink fixture is created");
  }
  inside.write("case.json", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":["profile.d"],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");

  ValidatedModel model{};
  return expect(compile(inside.root(), model).code == StatusCode::invalid_case,
                "referenced-file symlink escape is rejected");
}

bool test_reference_swap_after_validation() {
  ScratchCase outside("outside-swap");
  ScratchCase inside("inside-swap");
  outside.write("escaped.d", "outside bytes\n");
  inside.write("profile.d", "inside bytes\n");
  inside.write("case.json", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":["profile.d"],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");

  g_swap_target = inside.root() / "profile.d";
  g_swap_source = outside.root() / "escaped.d";
  g_swap_open_count = 0;
  hundun::v04::detail::set_file_open_observer_for_test(
      swap_reference_before_open);
  ValidatedModel model{};
  const Status status = compile(inside.root(), model);
  hundun::v04::detail::set_file_open_observer_for_test(nullptr);
  return expect(status.code == StatusCode::invalid_case,
                "reference swap after validation is rejected");
}

bool test_hard_link_alias_rejected() {
  ScratchCase scratch("hard-link-alias");
  scratch.write("first.d", "shared bytes\n");
  std::error_code error;
  fs::create_hard_link(scratch.root() / "first.d",
                       scratch.root() / "second.d", error);
  if (error) {
    return expect(false, "hard-link alias fixture is created");
  }
  scratch.write("case.json", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":["first.d","second.d"],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");

  ValidatedModel model{};
  return expect(compile(scratch.root(), model).code == StatusCode::invalid_case,
                "hard-link aliases are rejected by file identity");
}

bool test_supported_enum_values() {
  ScratchCase scratch("supported-enums");
  scratch.write("case.json", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"tensor_stretched","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"wale"},"time":{"control":"adaptive_acoustic"}
  })json");

  ValidatedModel model{};
  const Status status = compile(scratch.root(), model);
  bool passed = expect(static_cast<bool>(status),
                       "supported non-default enums compile");
  passed &= expect(model.geometry == GeometryKind::tensor_stretched,
                   "tensor-stretched enum is published");
  passed &= expect(model.turbulence == TurbulenceKind::wale,
                   "WALE enum is published");
  passed &= expect(model.time_control == TimeControlKind::adaptive_acoustic,
                   "acoustic time-control enum is published");
  return passed;
}

bool test_default_turbulence() {
  ScratchCase scratch("default-turbulence");
  scratch.write("case.json", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "time":{"control":"fixed"}
  })json");

  ValidatedModel model{};
  const Status status = compile(scratch.root(), model);
  return expect(static_cast<bool>(status) &&
                    model.turbulence == TurbulenceKind::vreman_wall_function &&
                    model.fingerprint == empty_model_fingerprint(
                        GeometryKind::uniform,
                        TurbulenceKind::vreman_wall_function,
                        TimeControlKind::fixed),
                "omitted turbulence selects the Vreman wall-function default");
}

bool test_invalid_cases() {
  bool passed = true;
  {
    ScratchCase scratch("missing-json");
    ValidatedModel model{};
    passed &= expect(compile(scratch.root(), model).code ==
                         StatusCode::invalid_case,
                     "missing case.json is rejected");
  }

  passed &= rejects("recursive duplicate keys", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","kind":"tensor_stretched","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("unknown top-level key", R"json({
    "schema_version":1,"units":"SI","extra":0,
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("unknown nested key", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null,"extra":0},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("unknown units", R"json({
    "schema_version":1,"units":"CGS",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("parent path", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":["../table.d"],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("absolute data path", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":["/tmp/table.d"],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("wrong data suffix", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":["table.txt"],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("missing data file", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":["missing.d"],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("nested data file", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":["data/table.d"],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("nested STL file", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":"geometry/shape.stl"},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("body fitted geometry", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"body_fitted","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("AMR geometry", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"amr","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("constant density flow", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null},
    "flow":{"model":"constant_density","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("SIMPLE coupling", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"SIMPLE","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("PIMPLE coupling", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PIMPLE","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("wrong PISO corrector count", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":3},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("reacting flow", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":true},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("unsupported turbulence", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"smagorinsky"},"time":{"control":"fixed"}
  })json");
  passed &= rejects("missing pressure closure", R"json({
    "schema_version":1,"units":"SI",
    "mesh":{"kind":"uniform","data_files":[],"stl_file":null},
    "flow":{"model":"single_phase_low_mach_compressible","reacting":false},
    "solver":{"coupling":"PISO","pressure_correctors":2},
    "turbulence":{"model":"none"},"time":{"control":"fixed"}
  })json");
  return passed;
}

bool test_field_registry() {
  FieldRegistry registry;
  FieldId first = 99;
  FieldId second = 99;
  bool passed = expect(static_cast<bool>(registry.declare_field("alpha", 1, 0, first)),
                       "first field declaration succeeds");
  passed &= expect(static_cast<bool>(registry.declare_field("beta", 3, 2, second)),
                   "second field declaration succeeds");
  passed &= expect(first == 0 && second == 1,
                   "field IDs follow registration order");

  FieldId ignored = 0;
  passed &= expect(registry.declare_field("alpha", 1, 0, ignored).code ==
                       StatusCode::invalid_plan,
                   "duplicate field name is rejected");
  passed &= expect(registry.declare_field("", 1, 0, ignored).code ==
                       StatusCode::invalid_plan,
                   "empty field name is rejected");
  passed &= expect(registry.declare_field("zero_components", 0, 0, ignored).code ==
                       StatusCode::invalid_plan,
                   "zero components are rejected");
  passed &= expect(registry.declare_field("bad_ghost", 1, 255, ignored).code ==
                       StatusCode::invalid_plan,
                   "unsupported ghost width is rejected");

  FieldSchema schema;
  passed &= expect(static_cast<bool>(registry.freeze_for_test(schema)),
                   "synthetic schema freezes");
  passed &= expect(schema.size() == 2, "frozen schema is an immutable snapshot");
  passed &= expect(schema[0].id == 0 && schema[0].stable_name == "alpha" &&
                       schema[0].components == 1 && schema[0].ghost_width == 0,
                   "snapshot preserves first declaration");
  passed &= expect(schema[1].id == 1 && schema[1].stable_name == "beta" &&
                       schema[1].components == 3 && schema[1].ghost_width == 2,
                   "snapshot preserves second declaration");
  passed &= expect(registry.declare_field("late", 1, 0, ignored).code ==
                       StatusCode::invalid_plan,
                   "mutation after freeze is rejected");
  return passed;
}

bool test_field_id_overflow() {
  FieldRegistry registry;
  FieldId id = 0;
  for (std::uint32_t value = 0;
       value <= std::numeric_limits<FieldId>::max(); ++value) {
    const std::string name = "field_" + std::to_string(value);
    const Status status = registry.declare_field(name, 1, 0, id);
    if (!status || id != static_cast<FieldId>(value)) {
      return expect(false, "all representable FieldIds can be assigned");
    }
  }
  return expect(registry.declare_field("overflow", 1, 0, id).code ==
                    StatusCode::invalid_plan,
                "FieldId overflow is rejected");
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    std::cerr << "FAIL: MPI_Init\n";
    return 1;
  }

  bool passed = true;
  passed &= test_valid_fixture();
  passed &= test_case_json_symlink_escape();
  passed &= test_reference_symlink_escape();
  passed &= test_reference_swap_after_validation();
  passed &= test_hard_link_alias_rejected();
  passed &= test_supported_enum_values();
  passed &= test_default_turbulence();
  passed &= test_invalid_cases();
  passed &= test_field_registry();
  passed &= test_field_id_overflow();

  if (MPI_Finalize() != MPI_SUCCESS) {
    std::cerr << "FAIL: MPI_Finalize\n";
    return 1;
  }
  return passed ? 0 : 1;
}
