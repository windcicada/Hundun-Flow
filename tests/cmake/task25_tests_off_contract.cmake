# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT OR
   NOT DEFINED HUNDUN_CXX_COMPILER OR
   NOT DEFINED HUNDUN_CXX_FLAGS OR
   NOT DEFINED HUNDUN_TESTS_OFF_EXECUTION OR
   NOT DEFINED HUNDUN_TESTS_OFF_RUNTIME OR
   NOT DEFINED HUNDUN_TESTS_OFF_DIAGNOSTICS OR
   NOT DEFINED HUNDUN_TESTS_OFF_LINEAR OR
   NOT DEFINED HUNDUN_TESTS_OFF_FLOW OR
   NOT DEFINED HUNDUN_TESTS_OFF_HUNDUN)
  message(FATAL_ERROR "Task 25 tests-off contract inputs are incomplete")
endif()

separate_arguments(cxx_flags NATIVE_COMMAND "${HUNDUN_CXX_FLAGS}")
set(command "${HUNDUN_CXX_COMPILER}" ${cxx_flags} -std=c++17
  "-I${HUNDUN_SOURCE_ROOT}/diagnostics/include"
  "-I${HUNDUN_SOURCE_ROOT}/runtime/include")
set(consumer
  "${CMAKE_CURRENT_BINARY_DIR}/task25-tests-off-public-consumer.cpp")
file(WRITE "${consumer}"
  "#include \"hundun/runtime/halo_performance_counters.hpp\"\n"
  "#include \"hundun/diagnostics/performance_correctness.hpp\"\n"
  "#include <type_traits>\n"
  "static_assert(std::is_trivially_copyable_v<hundun::runtime::HaloPerformanceCounters>);\n"
  "int main(){ hundun::diagnostics::PerformanceWorkRecord w; return w.phase == 'W' ? 0 : 1; }\n")
execute_process(
  COMMAND ${command} -fsyntax-only "${consumer}"
  RESULT_VARIABLE compile_result
  ERROR_VARIABLE compile_error)
if(NOT compile_result EQUAL 0)
  message(FATAL_ERROR
    "Task 25 tests-off public consumer failed: ${compile_error}")
endif()

foreach(archive IN ITEMS
    "${HUNDUN_TESTS_OFF_EXECUTION}"
    "${HUNDUN_TESTS_OFF_RUNTIME}"
    "${HUNDUN_TESTS_OFF_DIAGNOSTICS}"
    "${HUNDUN_TESTS_OFF_LINEAR}"
    "${HUNDUN_TESTS_OFF_FLOW}")
  if(NOT EXISTS "${archive}")
    message(FATAL_ERROR
      "required Task 25 tests-off archive is missing: ${archive}")
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
     "set_allocation_counters_for_test|set_next_halo_performance_counters|material_pressure_halo_counters|combine_pressure_halo_counters|PerformanceArtifactTestAccess")
    message(FATAL_ERROR
      "Task 25 test-only symbol leaked into ${archive}")
  endif()
endforeach()

if(NOT EXISTS "${HUNDUN_TESTS_OFF_HUNDUN}")
  message(FATAL_ERROR
    "required Task 25 tests-off executable is missing: "
    "${HUNDUN_TESTS_OFF_HUNDUN}")
endif()
execute_process(
  COMMAND nm -C "${HUNDUN_TESTS_OFF_HUNDUN}"
  RESULT_VARIABLE executable_nm_result
  OUTPUT_VARIABLE executable_symbols
  ERROR_VARIABLE executable_nm_error)
if(NOT executable_nm_result EQUAL 0)
  message(FATAL_ERROR
    "nm failed for ${HUNDUN_TESTS_OFF_HUNDUN}: ${executable_nm_error}")
endif()
if(executable_symbols MATCHES
   "Stage2DiagnosticObserver|capture_diagnostic_observer_snapshot|performance_failure_injection")
  message(FATAL_ERROR
    "Task 25 observer or injection symbol leaked into tests-off hundun")
endif()
file(STRINGS "${HUNDUN_TESTS_OFF_HUNDUN}" executable_strings)
string(JOIN "\n" executable_text ${executable_strings})
if(executable_text MATCHES
   "HUNDUN_TASK25_DIAGNOSTIC_OBSERVER|HUNDUN_TASK25_PERFORMANCE_FAILURE")
  message(FATAL_ERROR
    "Task 25 observer or injection selector leaked into tests-off hundun")
endif()
