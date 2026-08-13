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
  "data_files":[],"stl_file":null
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
  "data_files":[],"stl_file":null
})json";

std::string case_json(std::string_view mesh,
                      std::string_view flow = R"json({
  "model":"single_phase_low_mach_compressible",
  "pressure_closure":"local_absolute_pressure_drho_dp",
  "reacting":false
})json",
                      std::string_view solver = R"json({
  "coupling":"PISO","pressure_correctors":2
})json",
                      std::string_view turbulence =
                          R"json({"model":"none"})json",
                      std::string_view time =
                          R"json({"control":"fixed"})json") {
  return std::string{"{\"schema_version\":1,\"units\":\"SI\",\"mesh\":"} +
         std::string(mesh) + ",\"flow\":" + std::string(flow) +
         ",\"solver\":" + std::string(solver) +
         (turbulence.empty()
              ? std::string{}
              : ",\"turbulence\":" + std::string(turbulence)) +
         ",\"time\":" + std::string(time) + "}";
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

bool rejects(std::string_view label, std::string_view json,
             StatusCode expected = StatusCode::invalid_case) {
  ScratchCase scratch(label);
  scratch.write("case.json", json);
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
  passed &= expect(model.mesh.limits.max_global_cells == 4096 &&
                       model.mesh.limits.max_memory_bytes_per_rank == 67108864,
                   "fixture publishes hard limits");
  passed &= expect(model.fingerprint != 0U, "fixture has fingerprint");
  return passed;
}

bool test_tensor_normalization_and_fingerprint() {
  ScratchCase first_case("tensor-first");
  first_case.write("case.json", case_json(kTensorMesh));
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
    "data_files":[],"stl_file":null
  })json";
  ScratchCase second_case("tensor-second");
  second_case.write("case.json", case_json(reordered));
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
    "data_files":[],"stl_file":null
  })json";
  ScratchCase exact_case("tensor-exact");
  exact_case.write("case.json", case_json(exact_tensor));
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

bool test_every_typed_mesh_field_affects_fingerprint() {
  constexpr std::string_view baseline = R"json({
    "kind":"tensor_stretched",
    "domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,
    "focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],
                      "target_spacing":[0.1,0.1,0.1]}],
    "limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},
    "data_files":[],"stl_file":null
  })json";
  constexpr std::string_view variants[] = {
      R"json({"kind":"tensor_stretched","domain":{"lower":[-0.1,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1.1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":null,"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,11,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.21,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.04,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.15,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.11,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.41,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.11,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10001,"max_memory_bytes_per_rank":67108864},"data_files":[],"stl_file":null})json",
      R"json({"kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[12,10,8],"base_spacing":[0.2,0.2,0.2],"minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.1,"focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.1,0.1,0.1]}],"limits":{"max_global_cells":10000,"max_memory_bytes_per_rank":67108865},"data_files":[],"stl_file":null})json",
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

bool test_mesh_rejections() {
  bool passed = true;
  passed &= rejects("unknown mesh key", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null,"extra":0})json"));
  passed &= rejects("unordered domain", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,1],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("non-finite mesh value", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1e999,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("non-finite derived domain span", case_json(R"json({
    "kind":"uniform","domain":{"lower":[-1.7e308,0,0],"upper":[1.7e308,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("uniform requires exact cells", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("uniform rejects base spacing", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":[0.25,0.25,0.25],
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("uniform minimum exceeds width", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.3,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("uniform requires unit growth", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1.1,
    "focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("exact cell limit", case_json(R"json({
    "kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":[4,4,4],"base_spacing":null,
    "minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,
    "focus_regions":[],"limits":{"max_global_cells":63,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("tensor requires base spacing", case_json(R"json({
    "kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":null,
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.2,
    "focus_regions":[],"limits":{"max_global_cells":100,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("focus target below minimum", case_json(R"json({
    "kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":[0.2,0.2,0.2],
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.2,
    "focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.04,0.1,0.1]}],
    "limits":{"max_global_cells":1000,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("focus target exceeds base", case_json(R"json({
    "kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":[0.2,0.2,0.2],
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.2,
    "focus_regions":[{"lower":[0.1,0.1,0.1],"upper":[0.4,0.4,0.4],"target_spacing":[0.21,0.1,0.1]}],
    "limits":{"max_global_cells":1000,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
  passed &= rejects("focus wholly outside domain", case_json(R"json({
    "kind":"tensor_stretched","domain":{"lower":[0,0,0],"upper":[1,1,1]},
    "exact_cells":null,"base_spacing":[0.2,0.2,0.2],
    "minimum_spacing":[0.05,0.05,0.05],"max_growth_ratio":1.2,
    "focus_regions":[{"lower":[2,2,2],"upper":[3,3,3],"target_spacing":[0.1,0.1,0.1]}],
    "limits":{"max_global_cells":1000,"max_memory_bytes_per_rank":1},
    "data_files":[],"stl_file":null})json"));
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
                    std::string{"{\"schema_version\":1,\"units\":\"SI\",\"mesh\":"} +
                        std::string(kUniformMesh) +
                        ",\"flow\":{\"model\":\"single_phase_low_mach_compressible\",\"model\":\"constant_density\",\"pressure_closure\":\"local_absolute_pressure_drho_dp\",\"reacting\":false},\"solver\":{\"coupling\":\"PISO\",\"pressure_correctors\":2},\"turbulence\":{\"model\":\"none\"},\"time\":{\"control\":\"fixed\"}}");
  passed &= rejects("unknown units", std::string{"{\"schema_version\":1,\"units\":\"CGS\",\"mesh\":"} +
                        std::string(kUniformMesh) +
                        ",\"flow\":{\"model\":\"single_phase_low_mach_compressible\",\"pressure_closure\":\"local_absolute_pressure_drho_dp\",\"reacting\":false},\"solver\":{\"coupling\":\"PISO\",\"pressure_correctors\":2},\"turbulence\":{\"model\":\"none\"},\"time\":{\"control\":\"fixed\"}}");
  passed &= rejects("unknown top-level key",
                    std::string{"{\"schema_version\":1,\"units\":\"SI\",\"extra\":0,\"mesh\":"} +
                        std::string(kUniformMesh) +
                        ",\"flow\":{\"model\":\"single_phase_low_mach_compressible\",\"pressure_closure\":\"local_absolute_pressure_drho_dp\",\"reacting\":false},\"solver\":{\"coupling\":\"PISO\",\"pressure_correctors\":2},\"turbulence\":{\"model\":\"none\"},\"time\":{\"control\":\"fixed\"}}");
  passed &= rejects("constant density flow", case_json(
      kUniformMesh,
      R"json({"model":"constant_density","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false})json"));
  passed &= rejects("missing pressure closure", case_json(
      kUniformMesh,
      R"json({"model":"single_phase_low_mach_compressible","reacting":false})json"));
  passed &= rejects("reacting flow", case_json(
      kUniformMesh,
      R"json({"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":true})json"));
  passed &= rejects("SIMPLE coupling", case_json(
      kUniformMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false})json",
      R"json({"coupling":"SIMPLE","pressure_correctors":2})json"));
  passed &= rejects("PIMPLE coupling", case_json(
      kUniformMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false})json",
      R"json({"coupling":"PIMPLE","pressure_correctors":2})json"));
  passed &= rejects("wrong corrector count", case_json(
      kUniformMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false})json",
      R"json({"coupling":"PISO","pressure_correctors":3})json"));
  passed &= rejects("unsupported turbulence", case_json(
      kUniformMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false})json",
      R"json({"coupling":"PISO","pressure_correctors":2})json",
      R"json({"model":"smagorinsky"})json"));

  const auto reference_mesh = [](std::string_view data, std::string_view stl) {
    return std::string{R"json({"kind":"uniform","domain":{"lower":[0,0,0],"upper":[1,1,1]},"exact_cells":[4,4,4],"base_spacing":null,"minimum_spacing":[0.1,0.1,0.1],"max_growth_ratio":1,"focus_regions":[],"limits":{"max_global_cells":64,"max_memory_bytes_per_rank":1},"data_files":)json"} +
           std::string(data) + ",\"stl_file\":" + std::string(stl) + "}";
  };
  passed &= rejects("parent path", case_json(reference_mesh("[\"../table.d\"]", "null")));
  passed &= rejects("absolute data path", case_json(reference_mesh("[\"/tmp/table.d\"]", "null")));
  passed &= rejects("nested data path", case_json(reference_mesh("[\"data/table.d\"]", "null")));
  passed &= rejects("nested STL", case_json(reference_mesh("[]", "\"geometry/body.stl\"")));
  passed &= rejects("wrong data suffix", case_json(reference_mesh("[\"table.txt\"]", "null")));
  passed &= rejects("missing data file", case_json(reference_mesh("[\"missing.d\"]", "null")));

  ScratchCase outside("outside-symlink");
  ScratchCase inside("inside-symlink");
  outside.write("profile.d", "outside bytes\n");
  std::error_code error;
  fs::create_symlink(outside.root() / "profile.d",
                     inside.root() / "profile.d", error);
  inside.write("case.json", case_json(reference_mesh("[\"profile.d\"]", "null")));
  ValidatedModel model;
  passed &= expect(!error && compile(inside.root(), model).code ==
                                StatusCode::invalid_case,
                   "referenced-file symlink escape is rejected");

  ScratchCase outside_json("outside-json");
  ScratchCase inside_json("inside-json");
  outside_json.write("case.json", case_json(kUniformMesh));
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
  passed &= expect(!error && compile(hardlinks.root(), model).code ==
                                StatusCode::invalid_case,
                   "hard-link aliases are rejected");
  return passed;
}

bool test_defaults_and_enums() {
  ScratchCase defaults("defaults");
  defaults.write("case.json", case_json(kUniformMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false})json",
                                         R"json({"coupling":"PISO","pressure_correctors":2})json", ""));
  ValidatedModel default_model;
  bool passed = expect(static_cast<bool>(compile(defaults.root(), default_model)) &&
                           default_model.turbulence ==
                               TurbulenceKind::vreman_wall_function,
                       "omitted turbulence selects the production default");
  ScratchCase enums("enums");
  enums.write("case.json", case_json(kTensorMesh, R"json({"model":"single_phase_low_mach_compressible","pressure_closure":"local_absolute_pressure_drho_dp","reacting":false})json",
                                      R"json({"coupling":"PISO","pressure_correctors":2})json",
                                      R"json({"model":"wale"})json",
                                      R"json({"control":"adaptive_acoustic"})json"));
  ValidatedModel enum_model;
  passed &= expect(static_cast<bool>(compile(enums.root(), enum_model)) &&
                       enum_model.turbulence == TurbulenceKind::wale &&
                       enum_model.time_control == TimeControlKind::adaptive_acoustic,
                   "supported turbulence and time enums compile");
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
  passed &= expect(static_cast<bool>(registry.freeze_for_test(schema)) &&
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
  passed &= test_every_typed_mesh_field_affects_fingerprint();
  passed &= test_mesh_rejections();
  passed &= test_case_and_reference_security();
  passed &= test_defaults_and_enums();
  passed &= test_field_registry();
  passed &= test_field_id_overflow();
  const int finalize_status = MPI_Finalize();
  return passed && finalize_status == MPI_SUCCESS ? 0 : 1;
}
