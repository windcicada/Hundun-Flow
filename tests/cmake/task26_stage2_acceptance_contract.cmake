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

file(READ "${HUNDUN_STAGE2_ACCEPTANCE_SCRIPT}" acceptance_script)
string(CONCAT script_runtime "py" "thon")
if(acceptance_script MATCHES
    "(^|[\n;|&])[ \t]*(${script_runtime}([0-9.]*)?|pip([0-9.]*)?)([ \t]|$)")
  message(FATAL_ERROR
    "Task 26 contract: script invokes a forbidden scripting runtime")
endif()
if(acceptance_script MATCHES
    "(^|[\n;|&])[ \t]*(curl|wget)([ \t]|$)")
  message(FATAL_ERROR
    "Task 26 contract: script invokes a network download command")
endif()
if(acceptance_script MATCHES
    "(^|[\n;|&])[ \t]*(rm|mv|cp|touch|install)([ \t]|$)"
    OR acceptance_script MATCHES
       "(^|[\n;|&])[ \t]*git[ \t]+(add|commit|checkout|reset|clean|restore|switch)"
    OR acceptance_script MATCHES
       "(^|[\n;|&])[ \t]*(sed|perl)[ \t]+-[A-Za-z]*i")
  message(FATAL_ERROR
    "Task 26 contract: script contains a source-modification command")
endif()
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
