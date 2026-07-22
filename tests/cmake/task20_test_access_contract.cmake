# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT OR
   NOT DEFINED HUNDUN_CXX_COMPILER OR
   NOT DEFINED HUNDUN_CXX_FLAGS OR
   NOT DEFINED HUNDUN_MPI_INCLUDE_DIRS)
  message(FATAL_ERROR "Task 20 tests-off contract inputs are incomplete")
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

execute_process(
  COMMAND ${command} -E -P
          "${HUNDUN_SOURCE_ROOT}/tests/cmake/task20_header_preprocess.cpp"
  RESULT_VARIABLE preprocess_result
  OUTPUT_VARIABLE preprocessed
  ERROR_VARIABLE preprocess_error)
if(NOT preprocess_result EQUAL 0)
  message(FATAL_ERROR "Task 20 tests-off preprocessing failed: ${preprocess_error}")
endif()
foreach(test_name IN ITEMS
    MaterialDensityPisoTestAccess
    MaterialDensityPisoDiagnosticsTestAccess
    vortex_source)
  string(FIND "${preprocessed}" "${test_name}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task 20 tests-off header retains ${test_name}")
  endif()
endforeach()

execute_process(
  COMMAND ${command} -fsyntax-only
          "${HUNDUN_SOURCE_ROOT}/tests/cmake/task20_header_preprocess.cpp"
  RESULT_VARIABLE compile_result
  ERROR_VARIABLE compile_error)
if(NOT compile_result EQUAL 0)
  message(FATAL_ERROR "Task 20 public consumer failed: ${compile_error}")
endif()
