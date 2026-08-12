# HUNDUN-FLOW Stage 0 and Stage 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (- [ ]) syntax for tracking.

**Goal:** Create an isolated, copyright-clean HUNDUN-FLOW repository and deliver its first independently testable C++17 runtime: strict JSON configuration, MPI lifecycle, three-dimensional structured decomposition, typed fields, arbitrary-width halo exchange, basic mesh/I/O, a stable plugin ABI, and a conservative second-order passive-scalar demonstration.

**Architecture:** Stage 0 establishes legal provenance and rejects legacy source contamination before solver code exists. Stage 1 builds only the downward runtime path applications -> sdk -> mesh/solver -> runtime. The first executable advances a periodic passive scalar, proving the runtime without importing reacting-flow, IBM, chemistry, TPDF, or spray implementation.

**Tech Stack:** C++17, CMake 3.21+, MPI-3, CTest, POSIX/Linux, dlopen, vendored yyjson 0.12.0 at tag commit 7871d321ff4cd8068c1f777c97975dc2fb640ab3, Apache-2.0 project license, MIT third-party license.

## Global Constraints

- New public repository: /home/wyf/code_dev/hundun-flow.
- Private audit directory: /home/wyf/code_dev/hundun-flow-private-audit. It is never nested in or committed to the public repository.
- Fixed read-only legal baseline: /home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray.
- Do not copy, translate, mechanically refactor, or imitate BOFFIN, COAST, COAST-2, their Fortran control flow, legacy ABI, array layout, comments, messages, input formats, Decomp, or Restart adapters.
- Stage 1 reuses zero COAST source files. Later selective C++ reuse requires a separate file-by-file audit and user-approved plan.
- Use current COAST only as a black-box source of functional requirements.
- Public name is HUNDUN-FLOW; repository is hundun-flow; executable and C++ root namespace are hundun.
- Minimums: C++17, CMake 3.21, MPI-3, Linux CPU, GCC 9+ or Clang 12+.
- No Python, Conda, PETSc, HYPRE, Boost, Kokkos, OpenFOAM runtime, online dependency fetch, or GPU requirement in configure, build, test, or runtime paths.
- A private audit script may use host Python standard library, but no Python file enters the public repository or installed runtime.
- All public builds work offline after yyjson is vendored.
- Do not implement pressure coupling, LES, IBM, chemistry, TPDF, TCR, particles, spray, WENO, or DG in this plan.
- Public source files use SPDX-License-Identifier: Apache-2.0. Vendored yyjson retains its MIT notice.
- NOTICE contains only WANG YUDONG and wangyudong@buaa.edu.cn.
- Use one fresh worker per task, sequentially. The coordinator reviews and closes it before dispatching the next. At most five workers may exist; only one implementation worker is active.
- Workers do not contact the user for approval. Questions return to the coordinator.

---

## Locked File Structure

~~~text
hundun-flow/
  .github/workflows/build.yml
  .gitignore
  CMakeLists.txt
  CMakePresets.json
  CONTRIBUTING.md
  DCO.md
  LICENSE
  NOTICE
  README.md
  THIRD_PARTY.md
  VERSION
  LICENSES/yyjson-MIT.txt
  cmake/HundunProvenanceGuard.cmake
  cmake/HundunSanitizers.cmake
  cmake/HundunWarnings.cmake
  third_party/yyjson/UPSTREAM.json
  third_party/yyjson/yyjson.c
  third_party/yyjson/yyjson.h
  applications/hundun/cli_options.hpp
  applications/hundun/cli_options.cpp
  applications/hundun/case_config_broadcast.hpp
  applications/hundun/case_config_broadcast.cpp
  applications/hundun/main.cpp
  config/include/hundun/config/case_config.hpp
  config/include/hundun/config/case_config_loader.hpp
  config/src/case_config_loader.cpp
  runtime/include/hundun/runtime/error.hpp
  runtime/include/hundun/runtime/types.hpp
  runtime/include/hundun/runtime/mpi_environment.hpp
  runtime/include/hundun/runtime/collective_status.hpp
  runtime/include/hundun/runtime/structured_decomposition.hpp
  runtime/include/hundun/runtime/field_descriptor.hpp
  runtime/include/hundun/runtime/field_registry.hpp
  runtime/include/hundun/runtime/field_storage.hpp
  runtime/include/hundun/runtime/field_view.hpp
  runtime/include/hundun/runtime/exchange_plan.hpp
  runtime/include/hundun/runtime/halo_exchange.hpp
  runtime/include/hundun/runtime/restart_binary.hpp
  runtime/include/hundun/runtime/vtk_legacy.hpp
  runtime/src/mpi_environment.cpp
  runtime/src/collective_status.cpp
  runtime/src/structured_decomposition.cpp
  runtime/src/field_registry.cpp
  runtime/src/field_storage.cpp
  runtime/src/exchange_plan.cpp
  runtime/src/halo_exchange.cpp
  runtime/src/restart_binary.cpp
  runtime/src/vtk_legacy.cpp
  mesh/include/hundun/mesh/uniform_structured_mesh.hpp
  mesh/src/uniform_structured_mesh.cpp
  sdk/include/hundun/sdk/plugin_api.h
  sdk/include/hundun/sdk/plugin_loader.hpp
  sdk/src/plugin_loader.cpp
  solver/include/hundun/solver/passive_scalar.hpp
  solver/src/passive_scalar.cpp
  cases/passive_scalar/case.json
  cases/passive_scalar/case_restart.json
  cases/passive_scalar/README.md
  tests/support/test_main.hpp
  tests/unit/test_case_config.cpp
  tests/unit/test_cli_options.cpp
  tests/unit/test_field_storage.cpp
  tests/unit/test_mesh.cpp
  tests/unit/test_plugin_loader.cpp
  tests/unit/test_restart.cpp
  tests/unit/test_vtk.cpp
  tests/unit/test_passive_scalar.cpp
  tests/mpi/test_collective_status.cpp
  tests/mpi/test_case_config_broadcast.cpp
  tests/mpi/test_decomposition.cpp
  tests/mpi/test_halo_exchange.cpp
  tests/mpi/test_passive_scalar_mpi.cpp
  tests/plugins/mock_plugin.c
  tests/cmake/provenance_fixture.cmake
  tests/acceptance/compare_scalar_vtk.cpp
  tests/acceptance/stage1_acceptance.sh
  docs/architecture/runtime.md
  docs/development/source-policy.md
  docs/references/scientific-sources.md
~~~

Ownership:

- config owns JSON syntax, strict key checks, cross-field validation, and resolved output.
- runtime owns MPI, indexing, fields, communication, Restart, and primitive output.
- mesh owns physical coordinates.
- sdk owns the external binary interface and loader; physical model callbacks are deferred.
- solver/passive_scalar owns only the verification equation.
- applications/hundun owns orchestration and no numerical kernels.

---

### Task 1: Isolated Repository And Legal Evidence

**Files:**
- Preserve: /home/wyf/code_dev/hundun-flow/AGENTS.md
- Preserve: /home/wyf/code_dev/hundun-flow/docs/design/2026-07-16-hundun-flow-clean-cpp-solver-design.md
- Preserve: /home/wyf/code_dev/hundun-flow/docs/plans/2026-07-16-hundun-flow-stage0-stage1.md
- Preserve: /home/wyf/code_dev/hundun-flow/docs/handoff/2026-07-16-hundun-flow-agent-handoff.md
- Create: /home/wyf/code_dev/hundun-flow/.gitignore
- Create: /home/wyf/code_dev/hundun-flow/LICENSE
- Create: /home/wyf/code_dev/hundun-flow/NOTICE
- Create: /home/wyf/code_dev/hundun-flow/DCO.md
- Create: /home/wyf/code_dev/hundun-flow/CONTRIBUTING.md
- Create: /home/wyf/code_dev/hundun-flow/README.md
- Create: /home/wyf/code_dev/hundun-flow/VERSION
- Create privately: /home/wyf/code_dev/hundun-flow-private-audit/BASELINE.txt
- Create privately: /home/wyf/code_dev/hundun-flow-private-audit/baseline_sha256.txt
- Create privately: /home/wyf/code_dev/hundun-flow-private-audit/BOFFIN_INDEPENDENCE_AUDIT.md

**Interfaces:**
- Consumes: approved HUNDUN-FLOW design.
- Produces: a standalone Git repository whose first commit belongs to HUNDUN-FLOW, plus an external immutable baseline manifest.

- [ ] **Step 1: Verify target safety**

~~~bash
test -d /home/wyf/code_dev/hundun-flow
test ! -e /home/wyf/code_dev/hundun-flow/.git
test -d /home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray
find /home/wyf/code_dev/hundun-flow -mindepth 1 -type f -printf '%P\n' | sort
git -C /home/wyf/code_dev/Coast_software status --short
~~~

Expected: the target contains only `AGENTS.md` and the three approved
design/plan/handoff documents listed above; the first three test commands exit
0. Stop if source, build artifacts, a Git repository, or any other file is
present. Record Coast_software status without cleaning it.

- [ ] **Step 2: Initialize the documentation-only public skeleton and private audit directory**

~~~bash
git -C /home/wyf/code_dev/hundun-flow init -b main
mkdir -p /home/wyf/code_dev/hundun-flow-private-audit
~~~

Expected: a main branch with no commits and only the approved handoff documents
in the working tree. Do not copy a source template.

- [ ] **Step 3: Freeze the baseline outside the repository**

BASELINE.txt:

~~~text
Legal comparison baseline:
/home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray

Purpose:
Private independence audit only. This directory is not an implementation reference.

Policy:
Stage 1 reuses zero files from BOFFIN, COAST, or COAST-2.
~~~

Generate:

~~~bash
find /home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray -type f -print0 | sort -z | xargs -0 sha256sum > /home/wyf/code_dev/hundun-flow-private-audit/baseline_sha256.txt
sha256sum /home/wyf/code_dev/hundun-flow-private-audit/baseline_sha256.txt
~~~

Expected: non-empty manifest and a printed digest. Record the digest in the audit report.

- [ ] **Step 4: Write legal files**

NOTICE is exactly:

~~~text
Copyright (c) 2026 WANG YUDONG
Contact: wangyudong@buaa.edu.cn
~~~

VERSION is:

~~~text
0.0.0-stage1
~~~

LICENSE is the unmodified Apache License 2.0 text. DCO.md contains DCO 1.1. CONTRIBUTING.md requires Signed-off-by and rejects uncertain-origin code.

README begins:

~~~markdown
# HUNDUN-FLOW

HUNDUN-FLOW is an extensible C++ framework for immersed reacting-flow simulation.

The name draws on the archaic image of Hundun associated with Dijiang in the
Classic of Mountains and Seas: a pouch-like form, red as cinnabar fire. The
project reinterprets that image as a bounded chamber containing intense
combustion, turbulent transport, and interacting physical processes.

This engineering interpretation is not presented as a literal translation.
The primary text reference is:
https://ctext.org/text.pl?if=gb&node=83583&show=parallel

Stage 1 currently provides the independent MPI runtime and passive-scalar
verification kernel. Reacting flow, IBM, TPDF-TCR, and spray are roadmap items,
not current capabilities.
~~~

.gitignore:

~~~gitignore
/build/
/build-*/
/output/
/.cache/
/.idea/
/.vscode/
*.o
*.a
*.so
*.dylib
*.mod
compile_commands.json
CTestTestfile.cmake
Testing/
~~~

- [ ] **Step 5: Create the private audit header**

BOFFIN_INDEPENDENCE_AUDIT.md records baseline path, manifest digest, public initial commit, zero Stage 1 reuse, and this rule:

~~~text
Any exact normalized five-line source match outside third_party is a manual-review failure until explained. No explanation may describe COAST as a HUNDUN-FLOW upstream.
~~~

- [ ] **Step 6: Scan the empty repository**

~~~bash
rg -n -i 'boffin|coast_legacy|domxch|coalesced_legacy_block|input\.d|ibm_mesh\.d' /home/wyf/code_dev/hundun-flow
find /home/wyf/code_dev/hundun-flow -type f \( -iname '*.f' -o -iname '*.for' -o -iname '*.f90' -o -iname '*.inc' \)
~~~

Expected: no matches.

- [ ] **Step 7: Commit**

~~~bash
git -C /home/wyf/code_dev/hundun-flow add .gitignore LICENSE NOTICE DCO.md CONTRIBUTING.md README.md VERSION AGENTS.md docs
git -C /home/wyf/code_dev/hundun-flow commit -s -m 'chore: establish HUNDUN-FLOW legal foundation'
git -C /home/wyf/code_dev/hundun-flow rev-parse HEAD
~~~

Expected: one root commit with Signed-off-by. Add its hash to the private report.

---

### Task 2: Offline CMake, Test Harness, And yyjson

**Files:**
- Create: CMakeLists.txt
- Create: CMakePresets.json
- Create: cmake/HundunWarnings.cmake
- Create: cmake/HundunSanitizers.cmake
- Create: cmake/HundunProvenanceGuard.cmake
- Create: tests/support/test_main.hpp
- Create: tests/cmake/provenance_fixture.cmake
- Create: third_party/yyjson/yyjson.c
- Create: third_party/yyjson/yyjson.h
- Create: third_party/yyjson/UPSTREAM.json
- Create: LICENSES/yyjson-MIT.txt
- Create: THIRD_PARTY.md

**Interfaces:**
- Produces CMake targets yyjson, hundun_options, hundun_warnings and helpers
  hundun_add_test and hundun_add_mpi_test.

- [ ] **Step 1: Write a failing provenance fixture**

The fixture creates one clean C++ tree and one dirty tree containing legacy.F90. It includes HundunProvenanceGuard.cmake, calls hundun_assert_clean_tree on both, and is registered twice: clean must pass; dirty has CTest property WILL_FAIL TRUE.

Run:

~~~bash
cmake -P tests/cmake/provenance_fixture.cmake
~~~

Expected: FAIL because the guard does not exist.

- [ ] **Step 2: Implement the provenance guard**

Required CMake function:

~~~cmake
function(hundun_assert_clean_tree root)
  file(GLOB_RECURSE forbidden_fortran
    "${root}/*.f" "${root}/*.F" "${root}/*.for"
    "${root}/*.f90" "${root}/*.F90" "${root}/*.inc")
  if(forbidden_fortran)
    message(FATAL_ERROR "Forbidden legacy-language files")
  endif()
endfunction()
~~~

Extend it to scan non-third-party C/C++/CMake text for case-insensitive tokens
boffin, coast_legacy, domxch, and coalesced_legacy_block. Exclude
`cmake/HundunProvenanceGuard.cmake` itself because it owns that policy list;
exclude no production source directory. Error messages include the offending
path and token. Add a fixture proving that a forbidden token in an ordinary
CMakeLists.txt is rejected while the guard does not reject itself.

- [ ] **Step 3: Vendor yyjson 0.12.0**

~~~bash
git clone --depth 1 --branch 0.12.0 https://github.com/ibireme/yyjson.git /tmp/hundun-yyjson-0.12.0
git -C /tmp/hundun-yyjson-0.12.0 rev-parse HEAD
install -m 0644 /tmp/hundun-yyjson-0.12.0/src/yyjson.c third_party/yyjson/yyjson.c
install -m 0644 /tmp/hundun-yyjson-0.12.0/src/yyjson.h third_party/yyjson/yyjson.h
install -m 0644 /tmp/hundun-yyjson-0.12.0/LICENSE LICENSES/yyjson-MIT.txt
~~~

Expected commit: 7871d321ff4cd8068c1f777c97975dc2fb640ab3.

UPSTREAM.json:

~~~json
{
  "name": "yyjson",
  "version": "0.12.0",
  "repository": "https://github.com/ibireme/yyjson",
  "commit": "7871d321ff4cd8068c1f777c97975dc2fb640ab3",
  "license": "MIT",
  "vendored_files": ["yyjson.c", "yyjson.h"]
}
~~~

THIRD_PARTY.md records the same provenance and LICENSES path.

- [ ] **Step 4: Implement deterministic build settings**

CMake contract:

- project languages C and CXX;
- C++17 required, compiler extensions off;
- options HUNDUN_BUILD_TESTS, HUNDUN_ENABLE_ASAN, HUNDUN_ENABLE_UBSAN;
- find_package MPI 3 REQUIRED COMPONENTS CXX;
- yyjson is a local STATIC library with position-independent code;
- GNU/Clang C++ warnings: -Wall -Wextra -Wpedantic -Wconversion -Wshadow;
- hundun_add_test creates an executable, links options/warnings, registers the
  test, and applies label unit;
- hundun_add_mpi_test accepts a rank count, registers through MPIEXEC, and
  applies label mpi;
- no FetchContent or ExternalProject.

CMakePresets defines debug, release, and asan with build/debug, build/release, and build/asan.

- [ ] **Step 5: Add the dependency-free test harness**

tests/support/test_main.hpp:

~~~cpp
#pragma once

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace hundun::test {
inline void check(bool condition, const char* expression,
                  const char* file, int line) {
  if (!condition) {
    throw std::runtime_error(std::string(file) + ":" + std::to_string(line) +
                             " check failed: " + expression);
  }
}
inline void check_near(double actual, double expected, double tolerance,
                       const char* file, int line) {
  if (std::abs(actual - expected) > tolerance) {
    throw std::runtime_error(std::string(file) + ":" + std::to_string(line) +
                             " numerical check failed");
  }
}
template <class Function>
int run(Function&& function) {
  try {
    function();
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
}
#define HUNDUN_CHECK(expr) \
  ::hundun::test::check(static_cast<bool>(expr), #expr, __FILE__, __LINE__)
#define HUNDUN_CHECK_NEAR(actual, expected, tolerance) \
  ::hundun::test::check_near((actual), (expected), (tolerance), __FILE__, __LINE__)
~~~

- [ ] **Step 6: Verify offline configure**

~~~bash
cmake --preset debug
cmake --build --preset debug -j 8
ctest --preset debug --output-on-failure
rg -n 'FetchContent|ExternalProject|pip|python|conda' CMakeLists.txt cmake third_party
~~~

Expected: build/tests pass; dependency scan has no matches.

- [ ] **Step 7: Commit**

~~~bash
git add CMakeLists.txt CMakePresets.json cmake tests/support tests/cmake third_party LICENSES THIRD_PARTY.md
git commit -s -m 'build: add offline CMake and vendored JSON parser'
~~~

---

### Task 3: Strict Typed Case Configuration

**Files:**
- Create: runtime/include/hundun/runtime/error.hpp
- Create: runtime/include/hundun/runtime/types.hpp
- Create: config/include/hundun/config/case_config.hpp
- Create: config/include/hundun/config/case_config_loader.hpp
- Create: config/src/case_config_loader.cpp
- Create: tests/unit/test_case_config.cpp

**Interfaces:**
- Produces:
  - CaseConfig load_case_config(const std::filesystem::path&)
  - void validate_case_config(const CaseConfig&)
  - std::string to_resolved_json(const CaseConfig&)

- [ ] **Step 1: Define shared types**

~~~cpp
namespace hundun::runtime {
struct Int3 { int x{}; int y{}; int z{}; };
struct Real3 { double x{}; double y{}; double z{}; };
struct Box3 { Int3 begin{}; Int3 end{}; };
inline std::int64_t volume(Int3 e) {
  return static_cast<std::int64_t>(e.x) * e.y * e.z;
}
class Error : public std::runtime_error {
 public:
  explicit Error(const std::string& message) : std::runtime_error(message) {}
};
class ConfigError final : public Error {
 public:
  ConfigError(std::string pointer, std::string message);
  const std::string& pointer() const noexcept;
 private:
  std::string pointer_;
};
}
~~~

CaseConfig:

~~~cpp
struct MeshConfig {
  runtime::Int3 cells;
  runtime::Real3 origin_m;
  runtime::Real3 length_m;
  std::array<bool, 3> periodic;
};
struct TimeConfig { double dt_s; int steps; };
struct TransportConfig {
  runtime::Real3 velocity_m_per_s;
  double diffusivity_m2_per_s;
};
struct OutputConfig {
  std::filesystem::path directory;
  int write_interval;
  int restart_interval;
};
struct CaseConfig {
  int schema_version;
  std::string case_name;
  std::optional<int> expected_ranks;
  MeshConfig mesh;
  TimeConfig time;
  TransportConfig transport;
  std::string initial_condition;
  OutputConfig output;
};
~~~

- [ ] **Step 2: Write failing tests**

Valid JSON contains schema_version, case.name, optional resources.expected_ranks, mesh cells/origin_m/length_m/periodic, time dt_s/steps, transport velocity/diffusivity, initial_condition.type, and output directory/intervals.

Tests require exact JSON Pointers for:

- malformed JSON;
- unknown root and nested keys;
- wrong scalar/array type;
- three-vector wrong length;
- cell count < 1;
- non-positive length or dt;
- steps < 0;
- write/restart interval < 1;
- absolute output path;
- output path containing ..;
- initial condition other than sine_x.

Representative assertion:

~~~cpp
bool rejected = false;
try {
  static_cast<void>(load_case_config(path_with_unknown_mesh_key));
} catch (const hundun::runtime::ConfigError& error) {
  rejected = error.pointer() == "/mesh/unexpected";
}
HUNDUN_CHECK(rejected);
~~~

Run and expect compile failure before implementation.

- [ ] **Step 3: Implement strict yyjson traversal**

The loader:

1. reads the complete file;
2. uses strict yyjson_read flags;
3. checks each object against an explicit allowed-key set;
4. preserves integer type checks;
5. rejects path escape after lexical normalization;
6. throws ConfigError with exact pointer;
7. emits deterministic resolved JSON key order.

Every relative input/output path is resolved against the directory containing
the case JSON, never against the caller's current working directory. The
resolved path must remain below that case directory. `--print-resolved` emits
normalized case-relative paths so its output is independent of the launch
directory.

Required internal helpers:

~~~cpp
void reject_unknown_keys(yyjson_val*, std::initializer_list<std::string_view>,
                         std::string_view pointer);
double require_number(yyjson_val*, const char*, std::string_view pointer);
int require_integer(yyjson_val*, const char*, std::string_view pointer);
std::string require_string(yyjson_val*, const char*, std::string_view pointer);
~~~

- [ ] **Step 4: Verify**

~~~bash
cmake --build --preset debug -j 8
ctest --test-dir build/debug -R test_case_config --output-on-failure
cmake --preset asan
cmake --build --preset asan -j 8
ctest --test-dir build/asan -R test_case_config --output-on-failure
~~~

Expected: all cases pass without sanitizer output.

- [ ] **Step 5: Commit**

~~~bash
git add runtime/include config tests/unit/test_case_config.cpp CMakeLists.txt
git commit -s -m 'feat: add strict typed case configuration'
~~~

---

### Task 4: CLI, MPI Lifecycle, And Collective Failure

**Files:**
- Create: applications/hundun/cli_options.hpp
- Create: applications/hundun/cli_options.cpp
- Create: runtime/include/hundun/runtime/mpi_environment.hpp
- Create: runtime/include/hundun/runtime/collective_status.hpp
- Create: runtime/src/mpi_environment.cpp
- Create: runtime/src/collective_status.cpp
- Create: tests/unit/test_cli_options.cpp
- Create: tests/mpi/test_collective_status.cpp

**Interfaces:**
- Produces parse_cli, MpiEnvironment, collective_status, require_expected_ranks.

- [ ] **Step 1: Write CLI tests**

~~~cpp
struct CliOptions {
  std::filesystem::path case_path;
  bool validate_only{false};
  bool print_resolved{false};
  bool show_version{false};
};
~~~

Declare `CliOptions`, `parse_cli(int argc, char** argv)`, and the usage
constant in `cli_options.hpp`; keep parsing implementation in
`cli_options.cpp`. `main.cpp` includes the header and does not redeclare the
type.

Require:

- hundun case.json;
- hundun case.json --validate;
- hundun case.json --print-resolved;
- hundun --version;
- reject unknown option, two paths, no arguments, and both modes.

Exact usage:

~~~text
usage: hundun <case.json> [--validate|--print-resolved] | hundun --version
~~~

- [ ] **Step 2: Implement CLI parsing**

Parsing happens before MPI. --version needs no case. All invalid forms throw runtime::Error with the exact usage line.

- [ ] **Step 3: Write MPI collective test**

~~~cpp
int main(int argc, char** argv) {
  hundun::runtime::MpiEnvironment mpi(argc, argv);
  return hundun::test::run([&] {
    const bool local_ok = mpi.rank() != 1;
    const auto status = hundun::runtime::collective_status(
        mpi.comm(), local_ok, local_ok ? "" : "rank-one failure");
    HUNDUN_CHECK(!status.ok);
    HUNDUN_CHECK(status.failing_rank == 1);
    HUNDUN_CHECK(status.message == "rank-one failure");
  });
}
~~~

Register with mpiexec -n 2. Expect link failure initially.

- [ ] **Step 4: Implement MPI contracts**

~~~cpp
class MpiEnvironment final {
 public:
  MpiEnvironment(int& argc, char**& argv);
  ~MpiEnvironment();
  MpiEnvironment(const MpiEnvironment&) = delete;
  MPI_Comm comm() const noexcept { return MPI_COMM_WORLD; }
  int rank() const noexcept { return rank_; }
  int size() const noexcept { return size_; }
  void barrier() const;
 private:
  bool owns_mpi_{false};
  int rank_{0};
  int size_{1};
};

struct CollectiveStatus {
  bool ok;
  int failing_rank;
  std::string message;
};
~~~

MpiEnvironment requests MPI_THREAD_FUNNELED, validates provided thread level, and finalizes only if it initialized MPI.

collective_status chooses the lowest failing rank with MPI_Allreduce(MPI_MIN), then broadcasts its message length and bytes. All-success returns {true,-1,""}. require_expected_ranks returns the same error on every rank without hanging.

- [ ] **Step 5: Verify**

~~~bash
cmake --build --preset debug -j 8
ctest --test-dir build/debug -R 'test_cli_options|test_collective_status' --output-on-failure
mpiexec -n 2 build/debug/test_collective_status
~~~

Expected: exit 0, no MPI abort warning, no residual process.

- [ ] **Step 6: Commit**

~~~bash
git add applications runtime tests CMakeLists.txt
git commit -s -m 'feat: add CLI and MPI lifecycle'
~~~

---

### Task 5: Three-Dimensional Structured Decomposition

**Files:**
- Create: runtime/include/hundun/runtime/structured_decomposition.hpp
- Create: runtime/src/structured_decomposition.cpp
- Create: tests/mpi/test_decomposition.cpp

**Interfaces:**
- Produces Cartesian process grid, half-open owned boxes, 26-neighbor lookup, and partition-independent global_cell_id.

- [ ] **Step 1: Write decomposition invariants**

Run with 1, 2, and 4 ranks on global extent 17x11x7. Gather owned boxes and prove:

- every local extent is positive;
- summed volume is 17*11*7;
- boxes do not overlap;
- every global cell appears once;
- global_cell_id values are unique and in range.

Representative:

~~~cpp
const auto decomposition = StructuredDecomposition::create(
    mpi.comm(), Int3{17, 11, 7}, std::array<bool, 3>{true, false, true});
std::int64_t local = volume(decomposition.local_extent());
std::int64_t global = 0;
MPI_Allreduce(&local, &global, 1, MPI_INT64_T, MPI_SUM, mpi.comm());
HUNDUN_CHECK(global == 17LL * 11LL * 7LL);
~~~

- [ ] **Step 2: Implement exact interface**

~~~cpp
class StructuredDecomposition final {
 public:
  static StructuredDecomposition create(
      MPI_Comm, Int3 global_extent, std::array<bool, 3> periodic);
  ~StructuredDecomposition();
  StructuredDecomposition(StructuredDecomposition&&) noexcept;
  StructuredDecomposition(const StructuredDecomposition&) = delete;
  MPI_Comm comm() const noexcept;
  Int3 global_extent() const noexcept;
  Int3 process_grid() const noexcept;
  Int3 process_coordinates() const noexcept;
  Box3 owned_box() const noexcept;
  Int3 local_extent() const noexcept;
  int neighbor_rank(Int3 offset) const;
  std::uint64_t global_cell_id(Int3 local_cell) const;
  Int3 global_cell(Int3 local_cell) const;
};
~~~

Use MPI_Dims_create and MPI_Cart_create with reorder=0. Axis split is q=n/p, r=n%p; coordinate c owns q+1 if c<r, with begin=c*q+min(c,r). Reject process-grid dimensions larger than global cells.

neighbor_rank accepts offsets in {-1,0,1}^3 except zero, wraps periodic axes, and returns MPI_PROC_NULL at physical boundaries.

global_cell_id is ((k*global_y)+j)*global_x+i with checked uint64 arithmetic.

- [ ] **Step 3: Verify**

~~~bash
cmake --build --preset debug -j 8
mpiexec -n 1 build/debug/test_decomposition
mpiexec -n 2 build/debug/test_decomposition
mpiexec -n 4 build/debug/test_decomposition
~~~

Expected: all pass.

- [ ] **Step 4: Commit**

~~~bash
git add runtime tests/mpi/test_decomposition.cpp CMakeLists.txt
git commit -s -m 'feat: add structured MPI decomposition'
~~~

---

### Task 6: Typed Field Registry, Storage, And Views

**Files:**
- Create: runtime/include/hundun/runtime/field_descriptor.hpp
- Create: runtime/include/hundun/runtime/field_registry.hpp
- Create: runtime/include/hundun/runtime/field_storage.hpp
- Create: runtime/include/hundun/runtime/field_view.hpp
- Create: runtime/src/field_registry.cpp
- Create: runtime/src/field_storage.cpp
- Create: tests/unit/test_field_storage.cpp

**Interfaces:**
- Produces declare_field, freeze, FieldStorage and typed FieldView.

- [ ] **Step 1: Freeze descriptor types**

~~~cpp
enum class FunctionSpace { cell_average, face_value, vertex_value,
                           element_dof, quadrature_point, particle };
enum class ScalarType { float64, int32, uint8 };
enum class RestartPolicy { persistent, transient };
enum class OutputPolicy { never, selected, always };
using FieldId = std::uint32_t;
struct FieldDescriptor {
  std::string name;
  std::string unit;
  std::string owner;
  FunctionSpace space;
  ScalarType scalar_type;
  std::uint32_t components;
  int ghost_width;
  bool conservative;
  RestartPolicy restart;
  OutputPolicy output;
};
~~~

Stage 1 allocates cell_average only. Other spaces may be declared but allocation rejects them explicitly.

- [ ] **Step 2: Write failing tests**

Cover duplicate/empty names, missing unit/owner, zero components, negative ghost, declaration after freeze, stable IDs, float64/int32/uint8, negative and positive ghost indexing, wrong template type, out-of-bounds access, and size overflow.

~~~cpp
FieldRegistry registry;
const auto id = registry.declare_field(FieldDescriptor{
  "passive_scalar", "1", "passive_scalar_solver",
  FunctionSpace::cell_average, ScalarType::float64, 1, 2, true,
  RestartPolicy::persistent, OutputPolicy::selected});
registry.freeze();
FieldStorage storage(registry, Int3{4,3,2});
auto q = storage.view<double>(id);
q(-2,-2,-2,0) = -1.0;
q(3,2,1,0) = 7.0;
HUNDUN_CHECK_NEAR(q(-2,-2,-2,0), -1.0, 0.0);
HUNDUN_CHECK_NEAR(q(3,2,1,0), 7.0, 0.0);
~~~

- [ ] **Step 3: Implement memory/index contract**

FieldRegistry owns descriptors and name lookup. FieldStorage owns one byte vector per field. FieldView contains pointer, interior extent, ghost width, components, and strides.

Linear index:

~~~cpp
auto ii = i + ghost_width;
auto jj = j + ghost_width;
auto kk = k + ghost_width;
auto linear = (((kk * padded_y + jj) * padded_x + ii) * components)
              + component;
~~~

Provide explicit ScalarType mappings for double, int32_t, and uint8_t.

- [ ] **Step 4: Verify**

~~~bash
cmake --build --preset debug -j 8
ctest --test-dir build/debug -R test_field_storage --output-on-failure
cmake --build --preset asan -j 8
ctest --test-dir build/asan -R test_field_storage --output-on-failure
~~~

Expected: pass without sanitizer report.

- [ ] **Step 5: Commit**

~~~bash
git add runtime tests/unit/test_field_storage.cpp CMakeLists.txt
git commit -s -m 'feat: add typed field storage'
~~~

---

### Task 7: Arbitrary-Width 26-Neighbor Halo Exchange

**Files:**
- Create: runtime/include/hundun/runtime/exchange_plan.hpp
- Create: runtime/include/hundun/runtime/halo_exchange.hpp
- Create: runtime/src/exchange_plan.cpp
- Create: runtime/src/halo_exchange.cpp
- Create: tests/mpi/test_halo_exchange.cpp

**Interfaces:**
- Produces ExchangePlan::create and HaloExchange::exchange.

- [ ] **Step 1: Define region interface**

~~~cpp
struct ExchangeRegion {
  Int3 offset;
  int neighbor_rank;
  Box3 send_box;
  Box3 receive_box;
};
class ExchangePlan {
 public:
  static ExchangePlan create(const StructuredDecomposition&,
                             Int3 local_extent, int ghost_width);
  const std::vector<ExchangeRegion>& regions() const noexcept;
  int ghost_width() const noexcept;
};
~~~

Field coordinates use interior [0,n), lower ghosts [-g,0), upper ghosts [n,n+g).

- [ ] **Step 2: Write failing MPI tests**

On 12x10x8, periodic all axes, 1/2/4 ranks, ghost widths 1 and 2:

~~~cpp
q(i,j,k,0) = static_cast<double>(
    decomposition.global_cell_id(Int3{i,j,k}));
~~~

Set all ghosts -1, exchange, then verify every face/edge/corner ghost equals the wrapped global_cell_id.

- [ ] **Step 3: Implement exchange**

Create 26 regions for offsets in {-1,0,1}^3 except zero. MPI_PROC_NULL regions remain visible but inactive. send_box takes owned boundary thickness g on nonzero axes and full span on zero axes; receive_box takes matching ghosts. Reject a local dimension smaller than g.

HaloExchange packs contiguous byte buffers, posts all Irecv, posts all Isend, waits all, then unpacks. Support every Stage 1 ScalarType and component count.

Tag:

~~~cpp
int offset_code(Int3 o) {
  return (o.z + 1) * 9 + (o.y + 1) * 3 + (o.x + 1);
}
~~~

Send uses offset_code(offset); receive uses offset_code(-offset). Check every MPI return code.

- [ ] **Step 4: Verify**

~~~bash
cmake --build --preset debug -j 8
mpiexec -n 1 build/debug/test_halo_exchange
mpiexec -n 2 build/debug/test_halo_exchange
mpiexec -n 4 build/debug/test_halo_exchange
~~~

Expected: all widths, faces, edges, corners pass.

- [ ] **Step 5: Commit**

~~~bash
git add runtime tests/mpi/test_halo_exchange.cpp CMakeLists.txt
git commit -s -m 'feat: add arbitrary-width halo exchange'
~~~

---

### Task 8: Uniform Structured Mesh Baseline

**Files:**
- Create: mesh/include/hundun/mesh/uniform_structured_mesh.hpp
- Create: mesh/src/uniform_structured_mesh.cpp
- Create: tests/unit/test_mesh.cpp

**Interfaces:**
- Produces spacing, cell center, cell volume, local extent, owned box.

- [ ] **Step 1: Write geometry tests**

~~~cpp
class UniformStructuredMesh final {
 public:
  UniformStructuredMesh(Int3 global_cells, Real3 origin_m, Real3 length_m,
                        const StructuredDecomposition&);
  Real3 spacing_m() const noexcept;
  Real3 cell_center(Int3 local_cell) const;
  double cell_volume_m3() const noexcept;
  Int3 local_extent() const noexcept;
  Box3 owned_global_box() const noexcept;
};
~~~

For 4x2x2, origin (-1,0,2), length (2,1,4), require spacing (0.5,0.5,2), first center (-0.75,0.25,3), last global center (0.75,0.75,5), volume 0.5.

- [ ] **Step 2: Implement**

Cell center:

~~~cpp
x = origin.x + (global_i + 0.5) * spacing.x;
y = origin.y + (global_j + 0.5) * spacing.y;
z = origin.z + (global_k + 0.5) * spacing.z;
~~~

Reject non-positive dimensions/lengths. Do not add curvilinear metrics, STL, or IBM.

- [ ] **Step 3: Verify and commit**

~~~bash
cmake --build --preset debug -j 8
ctest --test-dir build/debug -R test_mesh --output-on-failure
git add mesh tests/unit/test_mesh.cpp CMakeLists.txt
git commit -s -m 'feat: add uniform structured mesh baseline'
~~~

Expected: test passes before commit.

---

### Task 9: Stable C Plugin ABI

**Files:**
- Create: sdk/include/hundun/sdk/plugin_api.h
- Create: sdk/include/hundun/sdk/plugin_loader.hpp
- Create: sdk/src/plugin_loader.cpp
- Create: tests/plugins/mock_plugin.c
- Create: tests/unit/test_plugin_loader.cpp

**Interfaces:**
- Produces symbol hundun_plugin_entry_v1 and RAII PluginLoader.

- [ ] **Step 1: Define ABI**

~~~c
#ifndef HUNDUN_SDK_PLUGIN_API_H
#define HUNDUN_SDK_PLUGIN_API_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
#define HUNDUN_PLUGIN_ABI_VERSION 1u
typedef struct HundunPluginDescriptorV1 {
  uint32_t abi_version;
  const char* name;
  const char* version;
  const char* description;
} HundunPluginDescriptorV1;
typedef const HundunPluginDescriptorV1* (*HundunPluginEntryV1)(void);
const HundunPluginDescriptorV1* hundun_plugin_entry_v1(void);
#ifdef __cplusplus
}
#endif
#endif
~~~

No C++ type, MPI handle, FieldView, allocator, or model callback crosses Stage 1 ABI.

- [ ] **Step 2: Write loader tests**

Cover valid mock plugin, missing file, missing entry, ABI version 2, empty name, and independent load/unload.

- [ ] **Step 3: Implement**

~~~cpp
class PluginLoader final {
 public:
  explicit PluginLoader(const std::filesystem::path&);
  ~PluginLoader();
  PluginLoader(PluginLoader&&) noexcept;
  PluginLoader(const PluginLoader&) = delete;
  const HundunPluginDescriptorV1& descriptor() const noexcept;
};
~~~

Use dlopen with RTLD_NOW|RTLD_LOCAL, clear dlerror before dlsym, validate descriptor before storing, dlclose once.

- [ ] **Step 4: Verify and commit**

~~~bash
cmake --build --preset debug -j 8
ctest --test-dir build/debug -R test_plugin_loader --output-on-failure
nm -D build/debug/tests/mock_plugin.so | rg 'hundun_plugin_entry_v1'
git add sdk tests/plugins tests/unit/test_plugin_loader.cpp CMakeLists.txt
git commit -s -m 'feat: establish plugin ABI version one'
~~~

Expected: one entry symbol and passing tests.

---

### Task 10: Versioned Restart And Primitive VTK

**Files:**
- Create: runtime/include/hundun/runtime/restart_binary.hpp
- Create: runtime/include/hundun/runtime/vtk_legacy.hpp
- Create: runtime/src/restart_binary.cpp
- Create: runtime/src/vtk_legacy.cpp
- Create: tests/unit/test_restart.cpp
- Create: tests/unit/test_vtk.cpp

**Interfaces:**
- Produces write_restart_rank/read_restart_rank and write_vtk_rank.

- [ ] **Step 1: Freeze Restart v1**

Write fields individually, never dump a C++ struct:

~~~text
magic[8] HUNDUNR1
uint32 version 1
uint32 endian 0x01020304
int32 rank, ranks
int32 global_nx, global_ny, global_nz
int32 begin_x, begin_y, begin_z
int32 local_nx, local_ny, local_nz
int64 step
float64 time_s
uint32 field_count
~~~

Per field:

~~~text
uint32 name_bytes
UTF-8 name
uint32 scalar_type
uint32 components
uint64 value_count
owned values in k-j-i-component order
~~~

Ghosts and transient fields are not persisted.

- [ ] **Step 2: Write failing Restart tests**

Round-trip Float64 and Int32 persistent fields while a transient field and all ghosts retain sentinels. Check exact step/time. Corrupt magic and truncate file; both must fail. Reader also rejects version/endian/rank/box/type/component mismatch and trailing bytes.

- [ ] **Step 3: Implement checked atomic I/O**

Use write_exact/read_exact. Write final-name.tmp, flush/close, rename atomically. Rank filenames:

~~~text
Restart/step00000010/restart.rank000000.bin
Restart/step00000010/restart.rank000001.bin
~~~

The `step%08d` directory is part of Restart v1. A later write never overwrites
an earlier checkpoint. The reader accepts one exact step directory and requires
one rank file for the calling rank.

- [ ] **Step 4: Write failing VTK tests**

Legacy ASCII STRUCTURED_POINTS:

~~~text
# vtk DataFile Version 3.0
HUNDUN-FLOW Stage 1
ASCII
DATASET STRUCTURED_POINTS
DIMENSIONS nx ny nz
ORIGIN x y z
SPACING dx dy dz
CELL_DATA count
SCALARS passive_scalar double 1
LOOKUP_TABLE default
~~~

For cell data, dimensions are local cells + 1. Test dimensions, lower owned vertex, spacing, count, field name, and exact value count.

- [ ] **Step 5: Implement and verify**

Writer supports selected Float64 cell-average fields with one component, owned values only, setprecision(17). Unsupported Stage 1 fields fail explicitly.

VTK filenames are
`scalar.step%08d.rank%06d.vtk` below `output.directory`, so outputs from
different steps and ranks never overwrite one another.

~~~bash
cmake --build --preset debug -j 8
ctest --test-dir build/debug -R 'test_restart|test_vtk' --output-on-failure
~~~

Expected: pass.

- [ ] **Step 6: Commit**

~~~bash
git add runtime tests/unit/test_restart.cpp tests/unit/test_vtk.cpp CMakeLists.txt
git commit -s -m 'feat: add versioned restart and VTK output'
~~~

---

### Task 11: Conservative Second-Order Passive Scalar

**Files:**
- Create: solver/include/hundun/solver/passive_scalar.hpp
- Create: solver/src/passive_scalar.cpp
- Create: tests/unit/test_passive_scalar.cpp
- Create: tests/mpi/test_passive_scalar_mpi.cpp

**Interfaces:**
- Produces mc_limiter, PassiveScalarSolver::advance_ssprk2, global_mass, global_l1_error.

- [ ] **Step 1: Write limiter tests**

~~~cpp
HUNDUN_CHECK_NEAR(mc_limiter(1.0,1.0), 1.0, 0.0);
HUNDUN_CHECK_NEAR(mc_limiter(-1.0,-1.0), -1.0, 0.0);
HUNDUN_CHECK_NEAR(mc_limiter(1.0,-1.0), 0.0, 0.0);
HUNDUN_CHECK_NEAR(mc_limiter(1.0,3.0), 2.0, 0.0);
HUNDUN_CHECK_NEAR(mc_limiter(3.0,1.0), 2.0, 0.0);
~~~

Formula:

~~~cpp
if (left * right <= 0.0) return 0.0;
double sign = left > 0.0 ? 1.0 : -1.0;
return sign * std::min({2.0*std::abs(left), 2.0*std::abs(right),
                        0.5*std::abs(left+right)});
~~~

- [ ] **Step 2: Write periodic convergence test**

q(x,0)=1+0.2 sin(2*pi*x/L), velocity=(1,0,0), diffusion=0, all periodic. Run one period at CFL 0.35 with nx=32,64,128 and ny=nz=4.

Require volume-weighted L1 orders:

~~~cpp
double order_a = std::log(error32/error64) / std::log(2.0);
double order_b = std::log(error64/error128) / std::log(2.0);
HUNDUN_CHECK(order_a > 1.70);
HUNDUN_CHECK(order_b > 1.70);
HUNDUN_CHECK(relative_mass_error < 1.0e-12);
~~~

- [ ] **Step 3: Implement MUSCL and SSPRK2**

At each face:

~~~cpp
double slope_i = mc_limiter(q_i-q_im1, q_ip1-q_i);
double slope_ip1 = mc_limiter(q_ip1-q_i, q_ip2-q_ip1);
double q_left = q_i + 0.5*slope_i;
double q_right = q_ip1 - 0.5*slope_ip1;
double flux = velocity >= 0.0 ? velocity*q_left : velocity*q_right;
~~~

Apply in x/y/z and compute negative flux divergence. Reject nonzero diffusivity in Stage 1.

SSPRK2:

1. exchange q;
2. q_stage=q+dt*L(q);
3. exchange q_stage;
4. q=0.5*q_old+0.5*(q_stage+dt*L(q_stage)).

- [ ] **Step 4: Write MPI invariance test**

Run 64x8x8 for 40 steps on 1/2/4 ranks. Gather sorted global_cell_id,value pairs and require:

~~~text
max_abs(q1-q2) < 1e-13
max_abs(q1-q4) < 1e-13
relative mass difference < 1e-13
~~~

- [ ] **Step 5: Verify and commit**

~~~bash
cmake --build --preset release -j 8
ctest --test-dir build/release -R test_passive_scalar --output-on-failure
mpiexec -n 1 build/release/test_passive_scalar_mpi
mpiexec -n 2 build/release/test_passive_scalar_mpi
mpiexec -n 4 build/release/test_passive_scalar_mpi
git add solver tests/unit/test_passive_scalar.cpp tests/mpi/test_passive_scalar_mpi.cpp CMakeLists.txt
git commit -s -m 'feat: add conservative passive scalar solver'
~~~

Expected: both orders >1.70, mass and decomposition tolerances pass.

---

### Task 12: Integrated Executable And Stage 1 Gate

**Files:**
- Create: applications/hundun/case_config_broadcast.hpp
- Create: applications/hundun/case_config_broadcast.cpp
- Create: applications/hundun/main.cpp
- Create: cases/passive_scalar/case.json
- Create: cases/passive_scalar/case_restart.json
- Create: cases/passive_scalar/README.md
- Create: tests/mpi/test_case_config_broadcast.cpp
- Create: tests/acceptance/compare_scalar_vtk.cpp
- Create: tests/acceptance/stage1_acceptance.sh
- Create: docs/architecture/runtime.md
- Create: docs/development/source-policy.md
- Create: docs/references/scientific-sources.md
- Create: .github/workflows/build.yml
- Modify: README.md
- Modify: config types/loader for Restart input

**Interfaces:**
- Produces: mpiexec -n N ./hundun cases/passive_scalar/case.json.

- [ ] **Step 1: Create canonical case**

~~~json
{
  "schema_version": 1,
  "case": {"name": "periodic_passive_scalar"},
  "resources": {"expected_ranks": 2},
  "mesh": {
    "cells": [64, 8, 8],
    "origin_m": [0.0, 0.0, 0.0],
    "length_m": [1.0, 0.125, 0.125],
    "periodic": [true, true, true]
  },
  "time": {"dt_s": 0.003125, "steps": 20},
  "transport": {
    "velocity_m_per_s": [1.0, 0.0, 0.0],
    "diffusivity_m2_per_s": 0.0
  },
  "initial_condition": {"type": "sine_x"},
  "restart": {"read": false, "write_directory": "Restart"},
  "output": {
    "directory": "output",
    "write_interval": 10,
    "restart_interval": 10
  }
}
~~~

Extend CaseConfig with:

~~~cpp
struct RestartConfig {
  bool read;
  std::optional<std::filesystem::path> read_directory;
  std::filesystem::path write_directory;
};
~~~

Apply the same path safety rules.
Add `RestartConfig restart;` to `CaseConfig`, include every Restart member in
resolved JSON and typed MPI broadcast, and add exact-pointer tests for missing
`read_directory`, forbidden `read_directory` when `read=false`, unknown Restart
keys, absolute paths, and `..` escapes.

`time.steps` is the absolute target step. With a Restart at step 10 and
`time.steps=20`, the executable advances steps 11 through 20. Writes go below
`restart.write_directory/step%08d`; when `read=true`,
`restart.read_directory` is required and names one exact step directory. When
`read=false`, `read_directory` must be absent.

`case_restart.json` is identical to `case.json` except:

~~~json
"restart": {
  "read": true,
  "read_directory": "Restart/step00000010",
  "write_directory": "Restart.resumed"
},
"output": {
  "directory": "output.resumed",
  "write_interval": 10,
  "restart_interval": 10
}
~~~

- [ ] **Step 2: Write failing acceptance script**

The Bash script uses set -euo pipefail and:

1. builds Release from clean;
2. checks --version, --validate, deterministic --print-resolved;
3. runs two ranks;
4. requires STEP 20 and FINISHED;
5. requires two-rank Restart and VTK outputs;
6. restarts from step 10 and compares final values with uninterrupted step 20 at 1e-13;
7. runs all CTest tests;
8. checks ldd for Python, missing libraries, COAST, BOFFIN;
9. runs provenance guard.

Required checks:

~~~bash
test "$(find output -name 'scalar.*.vtk' | wc -l)" -ge 2
test "$(find Restart/step00000020 -name 'restart.rank*.bin' | wc -l)" -eq 2
rg -n '^HUNDUN-FLOW 0\.0\.0-stage1$' output/run.log
rg -n '^FINISHED step=20 ' output/run.log
if ldd build/release/hundun | rg -i 'python|not found|coast|boffin'; then
  exit 1
fi
for rank in 000000 000001; do
  build/release/compare_scalar_vtk \
    "output/scalar.step00000020.rank${rank}.vtk" \
    "output.resumed/scalar.step00000020.rank${rank}.vtk" \
    1e-13
done
~~~

Expected: fail before main integration.

The acceptance script works in a fresh temporary case directory, copies the
two canonical JSON files into it, and removes the temporary directory through
a shell `trap`. It saves uninterrupted step-20 VTK fields before the restart
branch. The restart branch points to `Restart/step00000010`, writes to
`Restart.resumed`, and writes VTK to `output.resumed`.
`compare_scalar_vtk` reads the numeric scalar payload after
`LOOKUP_TABLE default`, checks identical value counts, and requires every
absolute difference to be at most its command-line tolerance. The script runs
it on each uninterrupted/resumed rank-file pair with tolerance `1e-13`. No
Python or JSON text rewriting is used.

- [ ] **Step 3: Implement rank-zero parse and typed broadcast**

Main flow:

1. parse CLI before MPI for --version;
2. initialize MpiEnvironment;
3. rank 0 loads CaseConfig or records ConfigError;
4. collective_status synchronizes success/failure;
5. broadcast typed CaseConfig members, strings as length+bytes;
6. validate expected_ranks collectively;
7. --validate prints VALID and exits before field allocation;
8. --print-resolved prints only rank 0;
9. construct decomposition, mesh, registry, storage, exchange, solver;
10. initialize sine field or read Restart;
11. advance and write diagnostics/VTK/Restart;
12. final barrier, then rank 0 prints FINISHED.

Do not broadcast JSON text for reparsing.

The broadcast interface is declared once in
`applications/hundun/case_config_broadcast.hpp`:

~~~cpp
namespace hundun::application {
config::CaseConfig broadcast_case_config(
    MPI_Comm comm, int root, const config::CaseConfig* root_config);
}
~~~

`case_config_broadcast.cpp` broadcasts every typed scalar and array member in
schema order. Strings and relative paths use checked `uint64_t` byte lengths
followed by bytes. Tests cover empty optional rank, populated optional rank,
empty and non-empty strings, and a non-root null `root_config`.

Stdout must match these regular expressions:

~~~text
^HUNDUN-FLOW 0\.0\.0-stage1$
^CASE name=periodic_passive_scalar ranks=2 cells=64x8x8$
^STEP 10 time_s=[0-9.eE+-]+ mass=[0-9.eE+-]+ relative_mass_error=[0-9.eE+-]+$
^STEP 20 time_s=[0-9.eE+-]+ mass=[0-9.eE+-]+ relative_mass_error=[0-9.eE+-]+$
^FINISHED step=20 time_s=[0-9.eE+-]+$
~~~

No screen file.

- [ ] **Step 4: Document architecture and provenance**

runtime.md records dependency arrows, interfaces, MPI thread level, field coordinates, halo tags, Restart v1, and non-goals.

source-policy.md states:

~~~text
HUNDUN-FLOW is independently implemented.
BOFFIN is used only as a private legal comparison baseline.
COAST and COAST-2 are black-box scientific references, not source ancestors.
Stage 1 reuses zero files from those programs.
No Fortran source or legacy adapter is accepted.
~~~

scientific-sources.md records these sources and the decision each supports; it
does not cite COAST code:

- R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems* (2002),
  DOI `10.1017/CBO9780511791253`: conservative finite-volume formulation;
- B. van Leer, *Towards the Ultimate Conservative Difference Scheme. V. A
  Second-Order Sequel to Godunov's Method* (1979),
  DOI `10.1016/0021-9991(79)90145-1`: MUSCL reconstruction;
- B. van Leer, *Towards the Ultimate Conservative Difference Scheme. IV. A
  New Approach to Numerical Convection* (1977),
  DOI `10.1016/0021-9991(77)90095-X`: monotonized slope limiting;
- S. Gottlieb and C.-W. Shu, *Total Variation Diminishing Runge-Kutta
  Schemes* (1998), DOI `10.1090/S0025-5718-98-00913-2`: SSPRK2;
- MPI Forum, *MPI: A Message-Passing Interface Standard, Version 3.1*
  (2015), `https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report.pdf`:
  Cartesian topology and nonblocking communication contracts;
- C. M. Rhie and W. L. Chow, *Numerical Study of the Turbulent Flow Past an
  Airfoil with Trailing Edge Separation* (1983), DOI `10.2514/3.8284`:
  future collocated pressure-velocity coupling;
- R. I. Issa, *Solution of the Implicitly Discretised Fluid Flow Equations by
  Operator-Splitting* (1986), DOI `10.1016/0021-9991(86)90099-9`: future PISO
  pressure correction.

- [ ] **Step 5: Add CI**

`.github/workflows/build.yml` is complete and contains no dependency download
other than Ubuntu packages:

~~~yaml
name: stage1

on:
  push:
  pull_request:

jobs:
  build:
    runs-on: ubuntu-22.04
    strategy:
      fail-fast: false
      matrix:
        include:
          - compiler: gcc
            cc: gcc
            cxx: g++
            preset: debug
          - compiler: gcc
            cc: gcc
            cxx: g++
            preset: release
          - compiler: clang
            cc: clang
            cxx: clang++
            preset: debug
          - compiler: clang
            cc: clang
            cxx: clang++
            preset: release
    env:
      CC: ${{ matrix.cc }}
      CXX: ${{ matrix.cxx }}
    steps:
      - uses: actions/checkout@v4
      - name: Install build dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y cmake ninja-build openmpi-bin libopenmpi-dev
      - name: Configure
        run: cmake --preset ${{ matrix.preset }}
      - name: Build
        run: cmake --build --preset ${{ matrix.preset }} -j 2
      - name: Test
        run: ctest --preset ${{ matrix.preset }} --output-on-failure
      - name: Acceptance
        if: matrix.preset == 'release'
        run: bash tests/acceptance/stage1_acceptance.sh
      - name: Configure ASan
        if: matrix.compiler == 'gcc' && matrix.preset == 'debug'
        run: cmake --preset asan
      - name: Build ASan
        if: matrix.compiler == 'gcc' && matrix.preset == 'debug'
        run: cmake --build --preset asan -j 2
      - name: Test ASan
        if: matrix.compiler == 'gcc' && matrix.preset == 'debug'
        run: ctest --preset asan -L unit --output-on-failure
~~~

The vendored yyjson tree means CI does not use FetchContent, git clone, pip, or
another source dependency fetch.

- [ ] **Step 6: Run complete Stage 1 gate**

~~~bash
bash tests/acceptance/stage1_acceptance.sh
cmake --preset asan
cmake --build --preset asan -j 8
ctest --test-dir build/asan -L unit --output-on-failure
git status --short
~~~

Expected:

- all tests pass;
- scalar orders >1.70;
- mass error <1e-12;
- decomposition and Restart differences <1e-13;
- no Python linkage;
- no forbidden source/token;
- only intended Task 12 files are dirty before commit.

- [ ] **Step 7: Run private similarity audit**

Create private audit_similarity.py outside the public repository. It reads
Fortran, C, and C++ files from the fixed baseline and C/C++ files from the
candidate, excluding `third_party`; removes comments with language-aware
lexing; lowercases; collapses whitespace; hashes normalized five-line windows;
writes path/line matches; and exits nonzero for unexplained matches. It must not
compare only like-language file extensions, because the audit is intended to
detect translated legacy blocks as well as literal reuse.

~~~bash
python3 /home/wyf/code_dev/hundun-flow-private-audit/audit_similarity.py \
  --baseline /home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray \
  --candidate /home/wyf/code_dev/hundun-flow \
  --exclude third_party \
  --report /home/wyf/code_dev/hundun-flow-private-audit/stage1_similarity.json
~~~

Expected: zero unexplained matches. Record script hash, candidate commit, count, and human conclusion privately.

- [ ] **Step 8: Commit and tag**

~~~bash
git add applications cases config runtime mesh sdk solver tests docs .github README.md CMakeLists.txt
git commit -s -m 'feat: complete independent stage one runtime'
git tag -a stage1-runtime -m 'HUNDUN-FLOW independent Stage 1 runtime'
~~~

Do not publish until coordinator review.

---

## Coordinator Review Gates

After every task:

1. confirm only task files changed;
2. inspect for legacy names, Fortran-shaped arrays, copied comments, and hidden dependencies;
3. run exact tests instead of trusting worker output;
4. use a requirements reviewer, then a code-quality reviewer for runtime changes;
5. resolve findings before closing the worker;
6. close the worker immediately after acceptance or rejection.

Final commands:

~~~bash
git -C /home/wyf/code_dev/hundun-flow log --oneline --decorate
for commit in $(git -C /home/wyf/code_dev/hundun-flow rev-list stage1-runtime); do
  git -C /home/wyf/code_dev/hundun-flow show -s --format=%B "$commit" |
    rg '^Signed-off-by: .+ <[^>]+>$' >/dev/null || exit 1
done
git -C /home/wyf/code_dev/hundun-flow status --short
bash /home/wyf/code_dev/hundun-flow/tests/acceptance/stage1_acceptance.sh
cmake --build /home/wyf/code_dev/hundun-flow/build/asan -j 8
ctest --test-dir /home/wyf/code_dev/hundun-flow/build/asan -L unit --output-on-failure
~~~

Acceptance:

- clean public worktree;
- all commits have Signed-off-by;
- Stage 1 tag points to accepted commit;
- zero source reuse from BOFFIN/COAST/COAST-2;
- private evidence exists outside public repo;
- offline configure/build;
- all unit, MPI, numerical, Restart, plugin, I/O, and provenance tests pass;
- no residual MPI process;
- no claim that reacting flow, IBM, TPDF-TCR, or spray exists.

## Deferred Work

Separate approved plans are required for:

- Stage 2 variable-density pressure-based flow;
- Stage 3 static single-STL IBM and LES;
- Stage 4 reacting flow and ChemistryBackend;
- Stage 5 tightly bound TPDF-TCR;
- Stage 6 spray and evaporation;
- selective reuse of audited user-owned C++ assets;
- WENO, DG, multi-part IBM, or moving walls.
