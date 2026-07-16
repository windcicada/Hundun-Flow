# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

project(HundunFlow LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(HUNDUN_BUILD_TESTS "Build HUNDUN-FLOW tests" ON)
option(HUNDUN_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(HUNDUN_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(HundunProvenanceGuard)
include(HundunSanitizers)
include(HundunWarnings)

hundun_assert_clean_tree("${PROJECT_SOURCE_DIR}")

find_package(MPI 3 REQUIRED COMPONENTS CXX)

add_library(hundun_options INTERFACE)
target_compile_features(hundun_options INTERFACE cxx_std_17)
hundun_configure_sanitizers(hundun_options)

add_library(hundun_warnings INTERFACE)
hundun_configure_warnings(hundun_warnings)

add_library(yyjson STATIC third_party/yyjson/yyjson.c)
set_target_properties(yyjson PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(yyjson PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third_party/yyjson>")

function(hundun_add_test target)
  add_executable("${target}" ${ARGN})
  target_link_libraries("${target}" PRIVATE hundun_options hundun_warnings)
  add_test(NAME "${target}" COMMAND "${target}")
  set_tests_properties("${target}" PROPERTIES LABELS unit)
endfunction()

function(hundun_add_mpi_test target ranks)
  add_executable("${target}" ${ARGN})
  target_link_libraries("${target}" PRIVATE
    hundun_options hundun_warnings MPI::MPI_CXX)
  add_test(NAME "${target}"
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:${target}>" ${MPIEXEC_POSTFLAGS})
  set_tests_properties("${target}" PROPERTIES LABELS mpi)
endfunction()

if(HUNDUN_BUILD_TESTS)
  enable_testing()

  add_test(NAME provenance_clean
    COMMAND "${CMAKE_COMMAND}" -DHUNDUN_PROVENANCE_CASE=clean
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  add_test(NAME provenance_forbidden_fortran
    COMMAND "${CMAKE_COMMAND}" -DHUNDUN_PROVENANCE_CASE=fortran
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  add_test(NAME provenance_forbidden_token
    COMMAND "${CMAKE_COMMAND}" -DHUNDUN_PROVENANCE_CASE=token
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  add_test(NAME provenance_forbidden_uppercase_suffix
    COMMAND "${CMAKE_COMMAND}" -DHUNDUN_PROVENANCE_CASE=uppercase
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  add_test(NAME provenance_forbidden_fortran_diagnostic
    COMMAND "${CMAKE_COMMAND}" -DHUNDUN_PROVENANCE_CASE=verify_fortran
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  add_test(NAME provenance_forbidden_token_diagnostic
    COMMAND "${CMAKE_COMMAND}" -DHUNDUN_PROVENANCE_CASE=verify_token
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  add_test(NAME provenance_forbidden_uppercase_suffix_diagnostic
    COMMAND "${CMAKE_COMMAND}" -DHUNDUN_PROVENANCE_CASE=verify_uppercase
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")

  set_tests_properties(
    provenance_clean
    provenance_forbidden_fortran
    provenance_forbidden_token
    provenance_forbidden_uppercase_suffix
    provenance_forbidden_fortran_diagnostic
    provenance_forbidden_token_diagnostic
    provenance_forbidden_uppercase_suffix_diagnostic
    PROPERTIES LABELS provenance)
  set_tests_properties(
    provenance_forbidden_fortran
    provenance_forbidden_token
    provenance_forbidden_uppercase_suffix
    PROPERTIES WILL_FAIL TRUE)
endif()
