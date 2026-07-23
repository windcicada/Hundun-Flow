# SPDX-License-Identifier: Apache-2.0
if(NOT DEFINED HUNDUN_SOURCE_ROOT OR
   NOT DEFINED HUNDUN_CXX_COMPILER OR
   NOT DEFINED HUNDUN_CXX_FLAGS OR
   NOT DEFINED HUNDUN_MPI_INCLUDE_DIRS)
  message(FATAL_ERROR "Task22 tests-off contract inputs are incomplete")
endif()

separate_arguments(cxx_flags NATIVE_COMMAND "${HUNDUN_CXX_FLAGS}")
string(REPLACE "|" ";" mpi_include_dirs "${HUNDUN_MPI_INCLUDE_DIRS}")
set(command "${HUNDUN_CXX_COMPILER}" ${cxx_flags} -std=c++17)
foreach(include_dir IN ITEMS
    "${HUNDUN_SOURCE_ROOT}/diagnostics/include"
    "${HUNDUN_SOURCE_ROOT}/config/include"
    "${HUNDUN_SOURCE_ROOT}/flow/include"
    "${HUNDUN_SOURCE_ROOT}/boundary/include"
    "${HUNDUN_SOURCE_ROOT}/execution/include"
    "${HUNDUN_SOURCE_ROOT}/finite_volume/include"
    "${HUNDUN_SOURCE_ROOT}/linear/include"
    "${HUNDUN_SOURCE_ROOT}/mesh/include"
    "${HUNDUN_SOURCE_ROOT}/runtime/include")
  list(APPEND command "-I${include_dir}")
endforeach()
foreach(include_dir IN LISTS mpi_include_dirs)
  list(APPEND command "-I${include_dir}")
endforeach()

set(consumer "${CMAKE_CURRENT_BINARY_DIR}/task22_tests_off_consumer.cpp")
file(WRITE "${consumer}"
"#include \"hundun/flow/adaptive_time_control.hpp\"\n"
"#include \"hundun/diagnostics/time_control_diagnostics.hpp\"\n"
"int main(){ return 0; }\n")
execute_process(
  COMMAND ${command} -E -P "${consumer}"
  RESULT_VARIABLE preprocess_result
  OUTPUT_VARIABLE preprocessed
  ERROR_VARIABLE preprocess_error)
if(NOT preprocess_result EQUAL 0)
  message(FATAL_ERROR "Task22 tests-off preprocessing failed: ${preprocess_error}")
endif()
foreach(test_name IN ITEMS
    AdaptiveTimeControlTestAccess
    TimeControlDiagnosticsTestAccess
    TimeControlDiagnosticFault
    TimeControlWireMutation
    TimeControlLocalMutation
    set_recoverable_failure_reason
    set_outcome
    set_raw_fault
    set_fault)
  string(FIND "${preprocessed}" "${test_name}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task22 public headers retain ${test_name}")
  endif()
endforeach()
execute_process(
  COMMAND ${command} -fsyntax-only "${consumer}"
  RESULT_VARIABLE compile_result
  ERROR_VARIABLE compile_error)
if(NOT compile_result EQUAL 0)
  message(FATAL_ERROR "Task22 public consumer failed: ${compile_error}")
endif()

foreach(archive IN ITEMS
    "${HUNDUN_TESTS_OFF_FLOW}"
    "${HUNDUN_TESTS_OFF_DIAGNOSTICS}")
  if(NOT archive OR NOT EXISTS "${archive}")
    message(FATAL_ERROR
      "required Task22 tests-off archive is missing: ${archive}")
  endif()
  execute_process(
    COMMAND nm -C "${archive}"
    RESULT_VARIABLE nm_result
    OUTPUT_VARIABLE symbols
    ERROR_VARIABLE nm_error)
  if(NOT nm_result EQUAL 0)
    message(FATAL_ERROR "nm failed for ${archive}: ${nm_error}")
  endif()
  if(symbols MATCHES
     "AdaptiveTimeControlTestAccess|TimeControlDiagnosticsTestAccess|diagnostic_fault|set_raw_fault|raw_fault_observation|set_fault|set_recoverable_failure_reason|wire_mutation|set_local_mutation|set_outcome|set_attempt_observer|exercise_trusted_tail_attempt_observer|facade_cache_snapshot|material_facade_cache_values_for_ideal")
    message(FATAL_ERROR "Task22 test-only symbol leaked into ${archive}")
  endif()
endforeach()
