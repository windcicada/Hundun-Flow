# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "HUNDUN_SOURCE_ROOT is required")
endif()
if(NOT DEFINED HUNDUN_TEST_BINARY_ROOT)
  message(FATAL_ERROR "HUNDUN_TEST_BINARY_ROOT is required")
endif()

get_filename_component(HUNDUN_SOURCE_ROOT
  "${HUNDUN_SOURCE_ROOT}" ABSOLUTE)
get_filename_component(HUNDUN_TEST_BINARY_ROOT
  "${HUNDUN_TEST_BINARY_ROOT}" ABSOLUTE)
set(HUNDUN_MUTATION_ROOT
  "${HUNDUN_TEST_BINARY_ROOT}/source-layout-mutation-fixture")

file(REMOVE_RECURSE "${HUNDUN_MUTATION_ROOT}")
file(MAKE_DIRECTORY
  "${HUNDUN_MUTATION_ROOT}/include/hundun"
  "${HUNDUN_MUTATION_ROOT}/src")
file(WRITE "${HUNDUN_MUTATION_ROOT}/include/hundun/app_probe.hpp"
  "#pragma once\n")
file(WRITE "${HUNDUN_MUTATION_ROOT}/src/app_probe.cpp"
  "#define HUNDUN_TEST_HEADER \"tests/support/app_probe_test.hpp\"\n"
  "#include HUNDUN_TEST_HEADER\n")

execute_process(
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_LAYOUT_ROOT=${HUNDUN_MUTATION_ROOT}"
    -P "${HUNDUN_SOURCE_ROOT}/tests/cmake/source_layout_fixture.cmake"
  RESULT_VARIABLE HUNDUN_MUTATION_RESULT
  OUTPUT_VARIABLE HUNDUN_MUTATION_STDOUT
  ERROR_VARIABLE HUNDUN_MUTATION_STDERR)
set(HUNDUN_MUTATION_OUTPUT
  "${HUNDUN_MUTATION_STDOUT}${HUNDUN_MUTATION_STDERR}")

if(HUNDUN_MUTATION_RESULT EQUAL 0)
  message(FATAL_ERROR
    "source layout fixture accepted a macro-hidden tests/support path")
endif()
if(NOT HUNDUN_MUTATION_OUTPUT MATCHES
    "product source references a tests/ path")
  message(FATAL_ERROR
    "source layout fixture rejected the mutation for the wrong reason:\n"
    "${HUNDUN_MUTATION_OUTPUT}")
endif()

file(REMOVE_RECURSE "${HUNDUN_MUTATION_ROOT}")
