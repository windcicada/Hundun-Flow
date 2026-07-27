# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

foreach(required_variable IN ITEMS
    HUNDUN_SOURCE_ROOT
    HUNDUN_BINARY_ROOT
    HUNDUN_STAGE2_ACCEPTANCE_SCRIPT
    HUNDUN_STAGE2_CAPABILITY_LEDGER
    HUNDUN_CTEST_COMMAND)
  if(NOT DEFINED ${required_variable} OR "${${required_variable}}" STREQUAL "")
    message(FATAL_ERROR "Task 26 contract: missing ${required_variable}")
  endif()
endforeach()

if(NOT EXISTS "${HUNDUN_STAGE2_ACCEPTANCE_SCRIPT}")
  message(FATAL_ERROR
    "Task 26 contract: Stage 2 acceptance script does not exist")
endif()
if(NOT EXISTS "${HUNDUN_STAGE2_CAPABILITY_LEDGER}")
  message(FATAL_ERROR
    "Task 26 contract: Stage 2 capability ledger does not exist")
endif()

find_program(HUNDUN_BASH_COMMAND bash REQUIRED)
execute_process(
  COMMAND "${HUNDUN_BASH_COMMAND}" -n
          "${HUNDUN_STAGE2_ACCEPTANCE_SCRIPT}"
  RESULT_VARIABLE bash_syntax_result
  OUTPUT_VARIABLE bash_syntax_stdout
  ERROR_VARIABLE bash_syntax_stderr)
if(NOT bash_syntax_result EQUAL 0)
  message(FATAL_ERROR
    "Task 26 contract: acceptance script does not parse:\n"
    "${bash_syntax_stdout}${bash_syntax_stderr}")
endif()

# This final gate is deliberately a frozen artifact, not a shell language
# classifier. The exact-content guard covers command indirection, redirection
# and unenumerated external operations that a token denylist cannot classify
# soundly. The inventory below independently documents every accepted
# external/read-only command site plus the shell entry/control/output
# boundaries; each signature must occur exactly once in the accepted bytes.
# Reviewed shell-language sites are assignments and parameter expansion,
# functions, if/test/then/else/fi, case/esac, for/do/done, cd, command -v,
# printf and exit. External sites are dirname, pwd, sed, awk, sort, the one
# read-only Git query and the CTest registry/run operations named below.
set(HUNDUN_TASK26_APPROVED_GATE_SHA256
    "f45b4861a4ac59e9dd9a9d1bbb52efc3d5c5c23e58e218fc4f7394e87bbe9657")
set(HUNDUN_TASK26_ALLOWED_COMMAND_SITES
  [=[bash-shebang|#!/usr/bin/env bash]=]
  [=[shell-strict-mode|set -euo pipefail]=]
  [=[shell-control-boundary|fail() {]=]
  [=[bounded-failure-output|printf 'stage2 acceptance: %s\n' "$*" >&2]=]
  [=[source-path-query|CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P]=]
  [=[build-path-query|build_root=$(CDPATH= cd -- "${requested_build}" && pwd -P)]=]
  [=[cache-source-read|cache_source=$(sed -n 's/^HundunFlow_SOURCE_DIR:STATIC=//p' "${cache}")]=]
  [=[cache-type-read|cache_type=$(sed -n 's/^CMAKE_BUILD_TYPE:STRING=//p' "${cache}")]=]
  [=[cache-tests-read|cache_tests=$(sed -n 's/^HUNDUN_BUILD_TESTS:BOOL=//p' "${cache}")]=]
  [=[build-identity-read|awk '/performance_source_commit[ \t]*=/{getline]=]
  [=[source-head-read|source_commit=$(git -C "${source_root}" rev-parse HEAD]=]
  [=[ctest-path-query|command -v "${ctest_command}" >/dev/null 2>&1]=]
  [=[ctest-registry-query|"${ctest_command}" --test-dir "${build_root}" -N -R "${selection_regex}" 2>&1]=]
  [=[registry-name-read|  printf '%s\n' "${inventory_output}" |
    sed -n]=]
  [=[registry-count-read|  printf '%s\n' "${actual_tests}" |
    awk]=]
  [=[expected-name-sort|printf '%s\n' "${required_tests[@]}" | LC_ALL=C sort]=]
  [=[actual-name-sort|printf '%s\n' "${actual_tests}" | LC_ALL=C sort]=]
  [=[bounded-inventory-output|printf 'STAGE2_ACCEPTANCE_INVENTORY cardinality=%s\n' "${required_count}"]=]
  [=[bounded-inventory-names|  printf '%s\n' "${required_tests[@]}"
  exit 0]=]
  [=[bounded-run-output|printf 'STAGE2_ACCEPTANCE_RUN cardinality=%s build=%s\n']=]
  [=[ctest-serial-run|"${ctest_command}" --test-dir "${build_root}" --output-on-failure]=]
  [=[ctest-serial-limit|-j1 -R "${selection_regex}"]=])

function(hundun_task26_classify_gate_policy script_path allowed_sites
         output_class)
  if(NOT EXISTS "${script_path}")
    set(${output_class} "missing-file" PARENT_SCOPE)
    return()
  endif()

  file(SHA256 "${script_path}" candidate_sha256)
  if(NOT candidate_sha256 STREQUAL HUNDUN_TASK26_APPROVED_GATE_SHA256)
    set(${output_class} "unapproved-gate-content" PARENT_SCOPE)
    return()
  endif()

  set(candidate_allowed_sites "${allowed_sites}")
  list(JOIN candidate_allowed_sites "\n" candidate_inventory)
  list(JOIN HUNDUN_TASK26_ALLOWED_COMMAND_SITES "\n" approved_inventory)
  if(NOT candidate_inventory STREQUAL approved_inventory)
    set(${output_class} "allowed-command-inventory-mismatch" PARENT_SCOPE)
    return()
  endif()

  file(READ "${script_path}" candidate_script)
  foreach(command_site IN LISTS HUNDUN_TASK26_ALLOWED_COMMAND_SITES)
    string(FIND "${command_site}" "|" delimiter_offset)
    if(delimiter_offset LESS 1)
      set(${output_class} "allowed-command-inventory-mismatch" PARENT_SCOPE)
      return()
    endif()
    math(EXPR signature_offset "${delimiter_offset} + 1")
    string(SUBSTRING "${command_site}" "${signature_offset}" -1 signature)
    string(FIND "${candidate_script}" "${signature}" first_match)
    if(first_match LESS 0)
      set(${output_class} "allowed-command-inventory-mismatch" PARENT_SCOPE)
      return()
    endif()
    string(LENGTH "${signature}" signature_length)
    math(EXPR remainder_offset "${first_match} + ${signature_length}")
    string(SUBSTRING "${candidate_script}" "${remainder_offset}" -1 remainder)
    string(FIND "${remainder}" "${signature}" second_match)
    if(NOT second_match LESS 0)
      set(${output_class} "allowed-command-inventory-mismatch" PARENT_SCOPE)
      return()
    endif()
  endforeach()

  set(${output_class} "approved" PARENT_SCOPE)
endfunction()

function(hundun_task26_require_unapproved_fixture fixture_name mutation_text)
  set(fixture_path
      "${HUNDUN_BINARY_ROOT}/task26-contract-mutation-${fixture_name}.sh")
  file(READ "${HUNDUN_STAGE2_ACCEPTANCE_SCRIPT}" fixture_contents)
  string(REPLACE "#!/usr/bin/env bash\n"
                 "#!/usr/bin/env bash\n${mutation_text}\n"
                 fixture_contents "${fixture_contents}")
  file(WRITE "${fixture_path}" "${fixture_contents}")
  hundun_task26_classify_gate_policy(
    "${fixture_path}" "${HUNDUN_TASK26_ALLOWED_COMMAND_SITES}" observed_class)
  file(REMOVE "${fixture_path}")
  if(NOT observed_class STREQUAL "unapproved-gate-content")
    message(FATAL_ERROR
      "Task 26 contract: ${fixture_name} mutation classification mismatch: "
      "expected unapproved-gate-content, observed ${observed_class}")
  endif()
endfunction()

hundun_task26_classify_gate_policy(
  "${HUNDUN_STAGE2_ACCEPTANCE_SCRIPT}"
  "${HUNDUN_TASK26_ALLOWED_COMMAND_SITES}" acceptance_script_class)
if(NOT acceptance_script_class STREQUAL "approved")
  message(FATAL_ERROR
    "Task 26 contract: acceptance gate policy rejected the script: "
    "${acceptance_script_class}")
endif()

set(missing_inventory "${HUNDUN_TASK26_ALLOWED_COMMAND_SITES}")
list(REMOVE_ITEM missing_inventory
  [=[source-head-read|source_commit=$(git -C "${source_root}" rev-parse HEAD]=])
hundun_task26_classify_gate_policy(
  "${HUNDUN_STAGE2_ACCEPTANCE_SCRIPT}" "${missing_inventory}"
  missing_inventory_class)
if(NOT missing_inventory_class STREQUAL
   "allowed-command-inventory-mismatch")
  message(FATAL_ERROR
    "Task 26 contract: allowed-command inventory mutation was not rejected")
endif()

string(CONCAT script_runtime "py" "thon")
string(CONCAT download_tool "cu" "rl")
set(dynamic_runtime_mutation
  "if false\nthen\n  forbidden=${script_runtime}3\n  \"\${forbidden}\" --version\nfi")
set(git_fetch_mutation [=[if false
then
  /usr/bin/git -C . fetch
fi]=])
set(source_redirect_mutation [=[if false
then
  printf x > "${source_root}/CMakeLists.txt"
fi]=])
set(harmless_argument_mutation
  "if false\nthen\n  printf '%s\\n' ${script_runtime}3\nfi")
set(path_runtime_mutation
  "if false\nthen\n  /usr/bin/${script_runtime}3 --version\nfi")
set(env_download_mutation
  "if false\nthen\n  /usr/bin/env ${download_tool} --version\nfi")
set(path_file_mutation [=[if false
then
  /usr/bin/rm never-created
fi]=])
set(env_git_mutation [=[if false
then
  /usr/bin/env git -C . clean -fdx
fi]=])

hundun_task26_require_unapproved_fixture(
  dynamic-runtime "${dynamic_runtime_mutation}")
hundun_task26_require_unapproved_fixture(
  git-fetch "${git_fetch_mutation}")
hundun_task26_require_unapproved_fixture(
  source-redirection "${source_redirect_mutation}")
hundun_task26_require_unapproved_fixture(
  harmless-runtime-argument "${harmless_argument_mutation}")
hundun_task26_require_unapproved_fixture(
  path-runtime "${path_runtime_mutation}")
hundun_task26_require_unapproved_fixture(
  env-download "${env_download_mutation}")
hundun_task26_require_unapproved_fixture(
  path-file-operation "${path_file_mutation}")
hundun_task26_require_unapproved_fixture(
  env-git-operation "${env_git_mutation}")

file(READ "${HUNDUN_STAGE2_ACCEPTANCE_SCRIPT}" acceptance_script)
if(acceptance_script MATCHES "stage2_acceptance\\.sh")
  message(FATAL_ERROR
    "Task 26 contract: script recursively names or invokes itself")
endif()
if(acceptance_script MATCHES "acceptance_complete_gate")
  message(FATAL_ERROR
    "Task 26 contract: script selects the legacy acceptance gate")
endif()

set(missing_build
    "${HUNDUN_BINARY_ROOT}/task26-contract-deliberately-missing")
file(REMOVE_RECURSE "${missing_build}")
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env
          "HUNDUN_STAGE2_BUILD_DIR=${missing_build}"
          "CTEST_COMMAND=${HUNDUN_CTEST_COMMAND}"
          "${HUNDUN_BASH_COMMAND}" "${HUNDUN_STAGE2_ACCEPTANCE_SCRIPT}"
  RESULT_VARIABLE missing_build_result
  OUTPUT_VARIABLE missing_build_stdout
  ERROR_VARIABLE missing_build_stderr)
if(missing_build_result EQUAL 0)
  message(FATAL_ERROR
    "Task 26 contract: normal path accepted a missing build tree")
endif()
set(expected_missing_message
    "stage2 acceptance: build directory not found: ${missing_build}\n")
if(NOT "${missing_build_stdout}" STREQUAL ""
    OR NOT "${missing_build_stderr}" STREQUAL "${expected_missing_message}")
  message(FATAL_ERROR
    "Task 26 contract: missing-build rejection is not bounded and deterministic")
endif()

set(empty_build "${HUNDUN_BINARY_ROOT}/task26-contract-empty-build")
file(REMOVE_RECURSE "${empty_build}")
file(MAKE_DIRECTORY "${empty_build}")
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env
          "HUNDUN_STAGE2_BUILD_DIR=${empty_build}"
          "CTEST_COMMAND=${HUNDUN_CTEST_COMMAND}"
          "${HUNDUN_BASH_COMMAND}" "${HUNDUN_STAGE2_ACCEPTANCE_SCRIPT}"
  RESULT_VARIABLE empty_build_result
  OUTPUT_VARIABLE empty_build_stdout
  ERROR_VARIABLE empty_build_stderr)
file(REMOVE_RECURSE "${empty_build}")
if(empty_build_result EQUAL 0)
  message(FATAL_ERROR
    "Task 26 contract: normal path accepted an empty build tree")
endif()
set(expected_empty_message
    "stage2 acceptance: CMake cache not found: ${empty_build}/CMakeCache.txt\n")
if(NOT "${empty_build_stdout}" STREQUAL ""
    OR NOT "${empty_build_stderr}" STREQUAL "${expected_empty_message}")
  message(FATAL_ERROR
    "Task 26 contract: empty-build rejection is not bounded and deterministic")
endif()

set(required_tests
  test_stage2_task2_dispatch
  test_mesh_geometry_1_rank
  test_mesh_geometry_2_rank
  test_mesh_geometry_4_rank
  test_conjugate_gradient_1_rank
  test_conjugate_gradient_2_rank
  test_conjugate_gradient_4_rank
  test_bicgstab_1_rank
  test_bicgstab_2_rank
  test_bicgstab_4_rank
  test_matrix_free_poisson_1_rank
  test_matrix_free_poisson_2_rank
  test_matrix_free_poisson_4_rank
  test_basic_boundary_1_rank
  test_basic_boundary_2_rank
  test_basic_boundary_4_rank
  test_cell_centered_fvm_1_rank
  test_cell_centered_fvm_2_rank
  test_cell_centered_fvm_4_rank
  test_fixed_step_piso_1_rank
  test_fixed_step_piso_2_rank
  test_fixed_step_piso_4_rank
  test_taylor_green_piso_1_rank
  test_taylor_green_piso_2_rank
  test_taylor_green_piso_4_rank
  test_material_density_wave_1_rank
  test_material_density_wave_2_rank
  test_material_density_wave_4_rank
  test_variable_density_vortex_full_1_rank
  test_variable_density_vortex_full_2_rank
  test_variable_density_vortex_full_4_rank
  test_ideal_gas_closed_heating_full_1_rank
  test_ideal_gas_closed_heating_full_2_rank
  test_ideal_gas_closed_heating_full_4_rank
  test_ideal_gas_open_plug_full_1_rank
  test_ideal_gas_open_plug_full_2_rank
  test_ideal_gas_open_plug_full_4_rank
  test_adaptive_time_control_1_rank_acceptance
  test_adaptive_time_control_2_rank_acceptance
  test_adaptive_time_control_4_rank_acceptance
  test_checkpoint_v2_acceptance_1_rank
  test_checkpoint_v2_acceptance_2_rank
  test_checkpoint_v2_acceptance_4_rank
  test_checkpoint_v2_diagnostics_acceptance_1_rank
  test_checkpoint_v2_diagnostics_acceptance_2_rank
  test_checkpoint_v2_diagnostics_acceptance_4_rank
  test_task24_flow_models_1rank
  test_task24_flow_models_2rank
  test_task24_flow_models_4rank
  test_task24_restart_1rank
  test_task24_restart_2rank
  test_task24_restart_4rank
  test_diagnostic_session_mpi_1rank
  test_diagnostic_session_mpi_2rank
  test_diagnostic_session_mpi_4rank
  test_stage2_module_diagnostics
  test_task25_performance_1rank
  test_task25_performance_2rank
  test_task25_performance_4rank)
list(LENGTH required_tests required_test_count)
if(NOT required_test_count EQUAL 59)
  message(FATAL_ERROR
    "Task 26 contract: internal required-test cardinality is not 59")
endif()
list(JOIN required_tests "\n" expected_test_lines)
set(expected_inventory
    "STAGE2_ACCEPTANCE_INVENTORY cardinality=59\n${expected_test_lines}\n")

execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env
          "HUNDUN_STAGE2_LIST_ONLY=1"
          "HUNDUN_STAGE2_BUILD_DIR=${HUNDUN_BINARY_ROOT}"
          "CTEST_COMMAND=${HUNDUN_CTEST_COMMAND}"
          "${HUNDUN_BASH_COMMAND}" "${HUNDUN_STAGE2_ACCEPTANCE_SCRIPT}"
  RESULT_VARIABLE inventory_result
  OUTPUT_VARIABLE inventory_stdout
  ERROR_VARIABLE inventory_stderr)
if(NOT inventory_result EQUAL 0)
  message(FATAL_ERROR
    "Task 26 contract: inventory-only validation failed:\n"
    "${inventory_stdout}${inventory_stderr}")
endif()
if(NOT "${inventory_stderr}" STREQUAL "")
  message(FATAL_ERROR
    "Task 26 contract: inventory-only validation wrote stderr")
endif()
if(NOT "${inventory_stdout}" STREQUAL "${expected_inventory}")
  message(FATAL_ERROR
    "Task 26 contract: inventory is not the exact registered set\n"
    "expected:\n${expected_inventory}\nactual:\n${inventory_stdout}")
endif()

file(READ "${HUNDUN_STAGE2_CAPABILITY_LEDGER}" capability_ledger)

function(hundun_task26_classify_field_diagnostic_ledger
         ledger_text output_class)
  if(ledger_text MATCHES "deterministic layout/generation summaries")
    set(${output_class} "field-generation-overclaim" PARENT_SCOPE)
    return()
  endif()

  string(REPLACE "\n" ";" ledger_lines "${ledger_text}")
  set(provider_claim_present FALSE)
  set(deferred_disposition_present FALSE)
  set(deferred_boundary_present FALSE)
  foreach(ledger_line IN LISTS ledger_lines)
    if(ledger_line MATCHES "^\\| CAP-DIAG-FIELD \\|"
       AND ledger_line MATCHES
           "deterministic layout/descriptor/role summaries")
      set(provider_claim_present TRUE)
    endif()
    if(ledger_line MATCHES "^\\| CAP-DEFER-FIELD-GENERATION \\|"
       AND ledger_line MATCHES
           "field-storage/access/checked-view generation diagnostics"
       AND ledger_line MATCHES "\\| explicitly-deferred \\|")
      set(deferred_disposition_present TRUE)
      if(ledger_line MATCHES
          "registry/layout/roles-only")
        set(deferred_boundary_present TRUE)
      endif()
    endif()
  endforeach()

  if(NOT provider_claim_present)
    set(${output_class} "field-provider-claim-missing" PARENT_SCOPE)
    return()
  endif()
  if(NOT deferred_disposition_present)
    set(${output_class} "field-generation-disposition-missing"
        PARENT_SCOPE)
    return()
  endif()
  if(NOT deferred_boundary_present
     OR NOT ledger_text MATCHES
        "no `FieldStorage` generation query or test-only seam")
    set(${output_class} "field-generation-boundary-missing" PARENT_SCOPE)
    return()
  endif()
  set(${output_class} "" PARENT_SCOPE)
endfunction()

hundun_task26_classify_field_diagnostic_ledger(
  "${capability_ledger}" field_diagnostic_ledger_class)
if(NOT field_diagnostic_ledger_class STREQUAL "")
  message(FATAL_ERROR
    "Task 26 contract: field diagnostic ledger is invalid: "
    "${field_diagnostic_ledger_class}")
endif()

string(REPLACE "deterministic layout/descriptor/role summaries"
               "deterministic layout/generation summaries"
               overclaim_ledger "${capability_ledger}")
hundun_task26_classify_field_diagnostic_ledger(
  "${overclaim_ledger}" overclaim_ledger_class)
if(NOT overclaim_ledger_class STREQUAL "field-generation-overclaim")
  message(FATAL_ERROR
    "Task 26 contract: field-generation overclaim mutation was not rejected")
endif()

string(REPLACE "\n" ";" ledger_lines "${capability_ledger}")
set(missing_deferred_ledger "")
foreach(ledger_line IN LISTS ledger_lines)
  if(NOT ledger_line MATCHES "^\\| CAP-DEFER-FIELD-GENERATION \\|")
    string(APPEND missing_deferred_ledger "${ledger_line}\n")
  endif()
endforeach()
hundun_task26_classify_field_diagnostic_ledger(
  "${missing_deferred_ledger}" missing_deferred_ledger_class)
if(NOT missing_deferred_ledger_class STREQUAL
   "field-generation-disposition-missing")
  message(FATAL_ERROR
    "Task 26 contract: missing field-generation disposition was not rejected")
endif()

foreach(required_heading IN ITEMS
    "Approved requirement"
    "Disposition"
    "Implementation location"
    "Positive test"
    "Failure / rollback test"
    "MPI / numerical acceptance")
  if(NOT capability_ledger MATCHES "${required_heading}")
    message(FATAL_ERROR
      "Task 26 contract: ledger lacks ${required_heading}")
  endif()
endforeach()

string(REPLACE "\n" ";" ledger_lines "${capability_ledger}")
set(capability_row_count 0)
foreach(ledger_line IN LISTS ledger_lines)
  if(ledger_line MATCHES "^\\| CAP-[A-Z0-9-]+ \\|")
    math(EXPR capability_row_count "${capability_row_count} + 1")
    set(disposition_count 0)
    foreach(disposition IN ITEMS
        implemented-and-accepted
        explicitly-deferred
        out-of-scope)
      if(ledger_line MATCHES "\\| ${disposition} \\|")
        math(EXPR disposition_count "${disposition_count} + 1")
      endif()
    endforeach()
    if(NOT disposition_count EQUAL 1)
      message(FATAL_ERROR
        "Task 26 contract: every capability row needs exactly one disposition")
    endif()
    if(ledger_line MATCHES "\\|[ \t]*(TBD|none|N/A|—)[ \t]*\\|")
      message(FATAL_ERROR
        "Task 26 contract: capability row has a missing evidence field")
    endif()
  endif()
endforeach()
if(capability_row_count LESS 20)
  message(FATAL_ERROR
    "Task 26 contract: ledger does not cover the required capability set")
endif()

foreach(required_ledger_token IN ITEMS
    "schema-v2 same-executable dispatch"
    "frozen schema-v1 path"
    "field epoch/capability/kernel-view"
    "topology/geometry/warped mapping"
    "Buffer Halo"
    "VectorOps"
    "CG and BiCGStab"
    "five boundary kinds"
    "shared finite-volume mass flux"
    "constant-density PISO"
    "material-density transport and vortex"
    "ideal-gas closure"
    "adaptive BDF2 collective retry"
    "Checkpoint v2"
    "flow-driver"
    "MeshDiag v2"
    "warped numerical and performance evidence"
    "shared generic provider"
    "full model callback ABI"
    "production device/GPU execution"
    "vendor solver backend"
    "rank-changing checkpoint"
    "complex thermochemical boundary"
    "later-stage physics")
  if(NOT capability_ledger MATCHES "${required_ledger_token}")
    message(FATAL_ERROR
      "Task 26 contract: ledger lacks ${required_ledger_token}")
  endif()
endforeach()

foreach(provider IN ITEMS
    runtime/MPI
    mesh
    field
    execution
    linear
    boundary
    finite-volume
    PISO
    material-transport/flow
    ideal-gas-closure
    time-control
    checkpoint
    driver-integration)
  if(NOT capability_ledger MATCHES "provider: ${provider}")
    message(FATAL_ERROR
      "Task 26 contract: diagnostic coverage lacks provider: ${provider}")
  endif()
endforeach()
