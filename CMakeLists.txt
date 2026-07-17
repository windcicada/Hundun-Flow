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

add_library(hundun_config STATIC config/src/case_config_loader.cpp)
target_include_directories(hundun_config PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/config/include>"
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/runtime/include>")
target_link_libraries(hundun_config
  PRIVATE yyjson hundun_warnings
  PUBLIC hundun_options)

add_library(hundun_cli STATIC applications/hundun/cli_options.cpp)
target_include_directories(hundun_cli PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>")
target_link_libraries(hundun_cli
  PRIVATE hundun_warnings
  PUBLIC hundun_options)

add_library(hundun_runtime STATIC
  runtime/src/mpi_error.cpp
  runtime/src/mpi_environment.cpp
  runtime/src/mpi_context.cpp
  runtime/src/collective_status.cpp
  runtime/src/process_grid.cpp
  runtime/src/structured_decomposition.cpp
  runtime/src/field_registry.cpp
  runtime/src/field_storage.cpp
  runtime/src/halo_detail.cpp
  runtime/src/exchange_plan.cpp
  runtime/src/halo_exchange.cpp)
target_include_directories(hundun_runtime PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/runtime/include>")
target_link_libraries(hundun_runtime
  PRIVATE hundun_warnings
  PUBLIC hundun_options MPI::MPI_CXX)

add_library(hundun_mesh STATIC mesh/src/uniform_structured_mesh.cpp)
target_include_directories(hundun_mesh PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/mesh/include>")
target_link_libraries(hundun_mesh
  PRIVATE hundun_warnings
  PUBLIC hundun_options hundun_runtime)

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

  hundun_add_test(test_case_config tests/unit/test_case_config.cpp)
  target_include_directories(test_case_config PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_case_config PRIVATE hundun_config)

  hundun_add_test(test_process_grid
    tests/unit/test_process_grid.cpp
    runtime/src/process_grid.cpp)
  target_include_directories(test_process_grid PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/runtime/include")

  hundun_add_test(test_cli_options tests/unit/test_cli_options.cpp)
  target_include_directories(test_cli_options PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_cli_options PRIVATE hundun_cli)

  hundun_add_test(test_field_storage tests/unit/test_field_storage.cpp)
  target_include_directories(test_field_storage PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_field_storage PRIVATE hundun_runtime)

  hundun_add_test(test_halo_helpers tests/unit/test_halo_helpers.cpp)
  target_include_directories(test_halo_helpers PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/runtime/include")
  target_link_libraries(test_halo_helpers PRIVATE hundun_runtime)
  set_tests_properties(test_halo_helpers PROPERTIES TIMEOUT 20)

  hundun_add_mpi_test(test_mesh 2 tests/unit/test_mesh.cpp)
  target_include_directories(test_mesh PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_mesh PRIVATE hundun_mesh)
  set_tests_properties(test_mesh PROPERTIES
    TIMEOUT 30
    PROCESSORS 2)

  add_executable(test_halo_exchange tests/mpi/test_halo_exchange.cpp)
  target_include_directories(test_halo_exchange PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_halo_exchange PRIVATE
    hundun_options hundun_warnings hundun_runtime MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(halo_test "test_halo_exchange_${ranks}_rank")
    add_test(NAME "${halo_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_halo_exchange>"
              full ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${halo_test}" PROPERTIES
      LABELS mpi
      TIMEOUT 60
      PROCESSORS "${ranks}")
  endforeach()

  foreach(mode IN ITEMS state lifetime mismatch small_chunks failure)
    set(halo_mode_test "test_halo_exchange_${mode}")
    add_test(NAME "${halo_mode_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_halo_exchange>"
              "${mode}" ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${halo_mode_test}" PROPERTIES
      LABELS mpi
      TIMEOUT 60
      PROCESSORS 4)
  endforeach()

  add_test(NAME test_halo_exchange_finalized_idle
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_halo_exchange>"
            finalized_idle ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_halo_exchange_finalized_idle PROPERTIES
    LABELS mpi
    TIMEOUT 30
    PROCESSORS 1)

  add_test(NAME test_halo_exchange_unknown_mode
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_halo_exchange>"
            unknown_mode ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_halo_exchange_unknown_mode PROPERTIES
    LABELS mpi
    TIMEOUT 20
    PROCESSORS 1
    WILL_FAIL TRUE)

  hundun_add_mpi_test(test_collective_status 2
    tests/mpi/test_collective_status.cpp)
  target_include_directories(test_collective_status PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_collective_status PRIVATE hundun_runtime)
  set_tests_properties(test_collective_status PROPERTIES TIMEOUT 20)

  hundun_add_mpi_test(test_mpi_context 4
    tests/mpi/test_mpi_context.cpp)
  target_include_directories(test_mpi_context PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_mpi_context PRIVATE hundun_runtime)
  set_tests_properties(test_mpi_context PROPERTIES TIMEOUT 30 PROCESSORS 4)

  hundun_add_mpi_test(test_mpi_external_lifecycle 1
    tests/mpi/test_mpi_external_lifecycle.cpp)
  target_include_directories(test_mpi_external_lifecycle PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_mpi_external_lifecycle PRIVATE hundun_runtime)
  set_tests_properties(test_mpi_external_lifecycle PROPERTIES TIMEOUT 20)

  add_executable(test_decomposition tests/mpi/test_decomposition.cpp)
  target_include_directories(test_decomposition PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_decomposition PRIVATE
    hundun_options hundun_warnings hundun_runtime MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4 12)
    set(decomposition_test "test_decomposition_${ranks}_rank")
    add_test(NAME "${decomposition_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_decomposition>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${decomposition_test}" PROPERTIES
      LABELS mpi
      TIMEOUT 45
      PROCESSORS "${ranks}")
  endforeach()

  add_test(NAME test_decomposition_rank_mismatch
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_decomposition>"
            mismatch ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_decomposition_rank_mismatch PROPERTIES
    LABELS mpi
    TIMEOUT 30
    PROCESSORS 4)

  add_test(NAME test_decomposition_12_rank_shape
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "12"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_decomposition>"
            shape ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_decomposition_12_rank_shape PROPERTIES
    LABELS mpi
    TIMEOUT 30
    PROCESSORS 12)

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
    PROPERTIES
      LABELS provenance
      RESOURCE_LOCK hundun_provenance_fixture_tree)
  set_tests_properties(
    provenance_forbidden_fortran
    provenance_forbidden_token
    provenance_forbidden_uppercase_suffix
    PROPERTIES WILL_FAIL TRUE)
endif()
