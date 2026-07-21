# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT OR
   NOT DEFINED HUNDUN_CXX_COMPILER OR
   NOT DEFINED HUNDUN_CXX_FLAGS OR
   NOT DEFINED HUNDUN_MPI_INCLUDE_DIRS)
  message(FATAL_ERROR "Task 18 test-access contract inputs are incomplete")
endif()

separate_arguments(cxx_flags NATIVE_COMMAND "${HUNDUN_CXX_FLAGS}")
string(REPLACE "|" ";" mpi_include_dirs "${HUNDUN_MPI_INCLUDE_DIRS}")

set(preprocess_command
    "${HUNDUN_CXX_COMPILER}"
    ${cxx_flags}
    -std=c++17
    -E
    -P)
foreach(include_dir IN ITEMS
    "${HUNDUN_SOURCE_ROOT}/flow/include"
    "${HUNDUN_SOURCE_ROOT}/boundary/include"
    "${HUNDUN_SOURCE_ROOT}/config/include"
    "${HUNDUN_SOURCE_ROOT}/execution/include"
    "${HUNDUN_SOURCE_ROOT}/finite_volume/include"
    "${HUNDUN_SOURCE_ROOT}/linear/include"
    "${HUNDUN_SOURCE_ROOT}/mesh/include"
    "${HUNDUN_SOURCE_ROOT}/runtime/include")
  list(APPEND preprocess_command "-I${include_dir}")
endforeach()
foreach(include_dir IN LISTS mpi_include_dirs)
  list(APPEND preprocess_command "-I${include_dir}")
endforeach()
list(APPEND preprocess_command
     "${HUNDUN_SOURCE_ROOT}/tests/cmake/task18_header_preprocess.cpp")

execute_process(
  COMMAND ${preprocess_command}
  RESULT_VARIABLE preprocess_result
  OUTPUT_VARIABLE preprocessed_header
  ERROR_VARIABLE preprocess_error)
if(NOT preprocess_result EQUAL 0)
  message(FATAL_ERROR
    "Task 18 tests-off header preprocessing failed: ${preprocess_error}")
endif()

string(FIND "${preprocessed_header}" "ConstantDensityPisoTestAccess"
       test_access_position)
if(NOT test_access_position EQUAL -1)
  message(FATAL_ERROR
    "Task 18 tests-off header retains ConstantDensityPisoTestAccess")
endif()
