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
hundun_assert_public_dependency_policy("${PROJECT_SOURCE_DIR}")

find_package(MPI 3 REQUIRED COMPONENTS CXX)

add_library(hundun_options INTERFACE)
target_compile_features(hundun_options INTERFACE cxx_std_17)
hundun_configure_sanitizers(hundun_options)

add_library(hundun_warnings INTERFACE)
hundun_configure_warnings(hundun_warnings)

add_library(hundun_diagnostics STATIC
  diagnostics/src/performance_artifact.cpp)
target_include_directories(hundun_diagnostics PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/diagnostics/include>")
target_link_libraries(hundun_diagnostics
  PRIVATE hundun_warnings
  PUBLIC hundun_options)

add_library(yyjson STATIC third_party/yyjson/yyjson.c)
set_target_properties(yyjson PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(yyjson PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third_party/yyjson>")

add_library(hundun_config STATIC
  config/src/case_config_loader.cpp
  config/src/resolved_case_loader.cpp)
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

add_library(hundun_application STATIC
  applications/hundun/case_config_broadcast.cpp)
target_include_directories(hundun_application
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
  PRIVATE "${PROJECT_SOURCE_DIR}/runtime/src")
target_link_libraries(hundun_application
  PRIVATE hundun_runtime hundun_warnings
  PUBLIC hundun_config hundun_options MPI::MPI_CXX)

add_library(hundun_mesh STATIC mesh/src/uniform_structured_mesh.cpp)
target_include_directories(hundun_mesh PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/mesh/include>")
target_link_libraries(hundun_mesh
  PRIVATE hundun_warnings
  PUBLIC hundun_options hundun_runtime)

add_library(hundun_io STATIC
  runtime/src/restart_binary.cpp
  runtime/src/vtk_legacy.cpp)
target_include_directories(hundun_io PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/runtime/include>")
target_link_libraries(hundun_io
  PRIVATE hundun_warnings
  PUBLIC hundun_options hundun_runtime hundun_mesh)

add_library(hundun_solver STATIC solver/src/passive_scalar.cpp)
target_include_directories(hundun_solver
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/solver/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/runtime/src")
target_link_libraries(hundun_solver
  PRIVATE hundun_options hundun_warnings
  PUBLIC hundun_mesh hundun_runtime)

add_executable(hundun
  applications/hundun/main.cpp
  applications/hundun/stage1_driver.cpp)
target_include_directories(hundun PRIVATE "${PROJECT_SOURCE_DIR}")
target_link_libraries(hundun PRIVATE
  hundun_application
  hundun_cli
  hundun_config
  hundun_io
  hundun_solver
  hundun_options
  hundun_warnings
  MPI::MPI_CXX)

add_library(hundun_plugin_sdk STATIC sdk/src/plugin_loader.cpp)
target_include_directories(hundun_plugin_sdk
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/sdk/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/runtime/include")
target_link_libraries(hundun_plugin_sdk
  PRIVATE hundun_warnings
  PUBLIC hundun_options ${CMAKE_DL_LIBS})

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

  find_package(Threads REQUIRED)

  function(hundun_add_mock_plugin target)
    add_library("${target}" MODULE tests/plugins/mock_plugin.c)
    target_include_directories("${target}" PRIVATE
      "${PROJECT_SOURCE_DIR}/sdk/include")
    target_link_libraries("${target}" PRIVATE hundun_options)
    set_target_properties("${target}" PROPERTIES
      C_STANDARD 11
      C_STANDARD_REQUIRED ON
      C_EXTENSIONS OFF
      PREFIX ""
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
      target_compile_options("${target}" PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    endif()
  endfunction()

  hundun_add_mock_plugin(mock_plugin)
  hundun_add_mock_plugin(mock_plugin_larger)
  target_compile_definitions(mock_plugin_larger PRIVATE HUNDUN_MOCK_LARGER)
  hundun_add_mock_plugin(mock_plugin_overlap)
  target_compile_definitions(mock_plugin_overlap PRIVATE HUNDUN_MOCK_OVERLAP)
  hundun_add_mock_plugin(mock_plugin_short)
  target_compile_definitions(mock_plugin_short PRIVATE HUNDUN_MOCK_SHORT)
  hundun_add_mock_plugin(mock_plugin_range_only)
  target_compile_definitions(mock_plugin_range_only PRIVATE
    HUNDUN_MOCK_RANGE_ONLY)
  hundun_add_mock_plugin(mock_plugin_no_overlap)
  target_compile_definitions(mock_plugin_no_overlap PRIVATE
    HUNDUN_MOCK_NO_OVERLAP)
  hundun_add_mock_plugin(mock_plugin_reversed)
  target_compile_definitions(mock_plugin_reversed PRIVATE
    HUNDUN_MOCK_REVERSED)
  hundun_add_mock_plugin(mock_plugin_capability)
  target_compile_definitions(mock_plugin_capability PRIVATE
    HUNDUN_MOCK_CAPABILITY)
  hundun_add_mock_plugin(mock_plugin_null_entry)
  target_compile_definitions(mock_plugin_null_entry PRIVATE
    HUNDUN_MOCK_NULL_ENTRY)
  hundun_add_mock_plugin(mock_plugin_null_name)
  target_compile_definitions(mock_plugin_null_name PRIVATE
    HUNDUN_MOCK_NULL_NAME)
  hundun_add_mock_plugin(mock_plugin_empty_name)
  target_compile_definitions(mock_plugin_empty_name PRIVATE
    HUNDUN_MOCK_EMPTY_NAME)
  hundun_add_mock_plugin(mock_plugin_null_version)
  target_compile_definitions(mock_plugin_null_version PRIVATE
    HUNDUN_MOCK_NULL_VERSION)
  hundun_add_mock_plugin(mock_plugin_empty_version)
  target_compile_definitions(mock_plugin_empty_version PRIVATE
    HUNDUN_MOCK_EMPTY_VERSION)
  hundun_add_mock_plugin(mock_plugin_null_description)
  target_compile_definitions(mock_plugin_null_description PRIVATE
    HUNDUN_MOCK_NULL_DESCRIPTION)
  hundun_add_mock_plugin(mock_plugin_empty_description)
  target_compile_definitions(mock_plugin_empty_description PRIVATE
    HUNDUN_MOCK_EMPTY_DESCRIPTION)
  hundun_add_mock_plugin(mock_plugin_missing_entry)
  target_compile_definitions(mock_plugin_missing_entry PRIVATE
    HUNDUN_MOCK_MISSING_ENTRY)

  add_executable(test_plugin_loader tests/unit/test_plugin_loader.cpp)
  set(hundun_dlclose_wrap_supported FALSE)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
       CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 9)
      set(hundun_dlclose_wrap_supported TRUE)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND
           CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
      set(hundun_dlclose_wrap_supported TRUE)
    endif()
  endif()
  if(NOT hundun_dlclose_wrap_supported)
    message(FATAL_ERROR
      "test_plugin_loader requires Linux and GNU 9+ or Clang 12+ linker wrapping")
  endif()
  target_include_directories(test_plugin_loader PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/sdk/include")
  target_link_libraries(test_plugin_loader PRIVATE
    hundun_plugin_sdk hundun_options hundun_warnings Threads::Threads)
  target_link_options(test_plugin_loader PRIVATE
    "-Wl,--wrap=dlclose"
    "-Wl,--wrap=dlerror"
    "-Wl,--wrap=dlsym")
  if(HUNDUN_ENABLE_ASAN AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    execute_process(
      COMMAND "${CMAKE_CXX_COMPILER}"
              --print-file-name=libclang_rt.asan.so
      OUTPUT_VARIABLE hundun_clang_asan_runtime
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT EXISTS "${hundun_clang_asan_runtime}")
      message(FATAL_ERROR "Clang shared AddressSanitizer runtime was not found")
    endif()
    get_filename_component(hundun_clang_asan_runtime_dir
      "${hundun_clang_asan_runtime}" DIRECTORY)
    target_link_options(test_plugin_loader PRIVATE "-shared-libasan")
    set_property(TARGET test_plugin_loader APPEND PROPERTY BUILD_RPATH
      "${hundun_clang_asan_runtime_dir}")
  endif()
  add_dependencies(test_plugin_loader
    mock_plugin
    mock_plugin_larger
    mock_plugin_overlap
    mock_plugin_short
    mock_plugin_range_only
    mock_plugin_no_overlap
    mock_plugin_reversed
    mock_plugin_capability
    mock_plugin_null_entry
    mock_plugin_null_name
    mock_plugin_empty_name
    mock_plugin_null_version
    mock_plugin_empty_version
    mock_plugin_null_description
    mock_plugin_empty_description
    mock_plugin_missing_entry)
  add_test(NAME test_plugin_loader COMMAND test_plugin_loader
    "$<TARGET_FILE:mock_plugin>"
    "$<TARGET_FILE:mock_plugin_larger>"
    "$<TARGET_FILE:mock_plugin_overlap>"
    "$<TARGET_FILE:mock_plugin_short>"
    "$<TARGET_FILE:mock_plugin_range_only>"
    "$<TARGET_FILE:mock_plugin_no_overlap>"
    "$<TARGET_FILE:mock_plugin_reversed>"
    "$<TARGET_FILE:mock_plugin_capability>"
    "$<TARGET_FILE:mock_plugin_null_entry>"
    "$<TARGET_FILE:mock_plugin_null_name>"
    "$<TARGET_FILE:mock_plugin_empty_name>"
    "$<TARGET_FILE:mock_plugin_null_version>"
    "$<TARGET_FILE:mock_plugin_empty_version>"
    "$<TARGET_FILE:mock_plugin_null_description>"
    "$<TARGET_FILE:mock_plugin_empty_description>"
    "$<TARGET_FILE:mock_plugin_missing_entry>")
  set_tests_properties(test_plugin_loader PROPERTIES
    LABELS unit
    TIMEOUT 30)

  hundun_add_test(test_case_config tests/unit/test_case_config.cpp)
  target_include_directories(test_case_config PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_case_config PRIVATE hundun_config)

  hundun_add_test(test_resolved_case tests/unit/test_resolved_case.cpp)
  target_include_directories(test_resolved_case PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_resolved_case PRIVATE hundun_config)

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

  hundun_add_test(test_performance_artifact
    tests/unit/test_performance_artifact.cpp)
  target_include_directories(test_performance_artifact PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/diagnostics/include")
  target_link_libraries(test_performance_artifact PRIVATE
    hundun_diagnostics yyjson)

  add_executable(test_passive_scalar tests/unit/test_passive_scalar.cpp)
  target_include_directories(test_passive_scalar PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/solver/include")
  target_link_libraries(test_passive_scalar PRIVATE
    hundun_solver hundun_options hundun_warnings MPI::MPI_CXX)
  add_test(NAME test_passive_scalar_limiter
    COMMAND test_passive_scalar limiter)
  set_tests_properties(test_passive_scalar_limiter PROPERTIES
    LABELS unit TIMEOUT 20)
  add_test(NAME test_passive_scalar_serial
    COMMAND test_passive_scalar serial)
  set_tests_properties(test_passive_scalar_serial PROPERTIES
    LABELS numerical TIMEOUT 300 PROCESSORS 1)
  add_test(NAME test_passive_scalar_oracle
    COMMAND test_passive_scalar oracle)
  set_tests_properties(test_passive_scalar_oracle PROPERTIES
    LABELS numerical TIMEOUT 60 PROCESSORS 1)

  add_executable(test_passive_scalar_mpi
    tests/mpi/test_passive_scalar_mpi.cpp)
  target_include_directories(test_passive_scalar_mpi PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/runtime/src"
    "${PROJECT_SOURCE_DIR}/solver/include")
  target_link_libraries(test_passive_scalar_mpi PRIVATE
    hundun_solver hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(passive_scalar_test "test_passive_scalar_mpi_${ranks}_rank")
    add_test(NAME "${passive_scalar_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_passive_scalar_mpi>"
              full ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${passive_scalar_test}" PROPERTIES
      LABELS "mpi;numerical" TIMEOUT 300 PROCESSORS "${ranks}")
  endforeach()
  add_test(NAME test_passive_scalar_mpi_finalized
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_passive_scalar_mpi>"
            finalized ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_passive_scalar_mpi_finalized PROPERTIES
    LABELS mpi TIMEOUT 30 PROCESSORS 1)
  foreach(mode IN ITEMS halo_compatibility halo_diagnostics)
    set(passive_scalar_mode_test "test_passive_scalar_mpi_${mode}")
    add_test(NAME "${passive_scalar_mode_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_passive_scalar_mpi>"
              "${mode}" ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${passive_scalar_mode_test}" PROPERTIES
      LABELS "mpi;numerical" TIMEOUT 60 PROCESSORS 4)
  endforeach()

  hundun_add_test(test_restart tests/unit/test_restart.cpp)
  target_include_directories(test_restart PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_restart PRIVATE hundun_io)
  set_tests_properties(test_restart PROPERTIES TIMEOUT 30)

  hundun_add_mpi_test(test_vtk 1 tests/unit/test_vtk.cpp)
  target_include_directories(test_vtk PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_vtk PRIVATE hundun_io)
  set_tests_properties(test_vtk PROPERTIES TIMEOUT 30 PROCESSORS 1)

  add_executable(test_restart_checkpoint
    tests/mpi/test_restart_checkpoint.cpp)
  target_include_directories(test_restart_checkpoint PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_restart_checkpoint PRIVATE
    hundun_io hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2)
    set(restart_test "test_restart_checkpoint_${ranks}_rank")
    add_test(NAME "${restart_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_restart_checkpoint>"
              full ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${restart_test}" PROPERTIES
      LABELS mpi TIMEOUT 90 PROCESSORS "${ranks}")
  endforeach()
  add_test(NAME test_restart_checkpoint_failure
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_restart_checkpoint>"
            failure ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_restart_checkpoint_failure PROPERTIES
    LABELS mpi TIMEOUT 90 PROCESSORS 2)
  add_test(NAME test_restart_checkpoint_finalized
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_restart_checkpoint>"
            finalized ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_restart_checkpoint_finalized PROPERTIES
    LABELS mpi TIMEOUT 30 PROCESSORS 1)

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

  add_executable(test_case_config_broadcast
    tests/mpi/test_case_config_broadcast.cpp)
  target_include_directories(test_case_config_broadcast PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_case_config_broadcast PRIVATE
    hundun_application hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(mode IN ITEMS full mismatch)
    set(case_broadcast_test "test_case_config_broadcast_${mode}")
    add_test(NAME "${case_broadcast_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_case_config_broadcast>"
              "${mode}" ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${case_broadcast_test}" PROPERTIES
      LABELS mpi TIMEOUT 45 PROCESSORS 4)
  endforeach()
  add_test(NAME test_case_config_broadcast_path_failures
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_case_config_broadcast>"
            path_failures ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_case_config_broadcast_path_failures PROPERTIES
    LABELS mpi TIMEOUT 15 PROCESSORS 4)
  add_test(NAME test_case_config_broadcast_finalized
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_case_config_broadcast>"
            finalized ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_case_config_broadcast_finalized PROPERTIES
    LABELS mpi TIMEOUT 20 PROCESSORS 1)

  add_executable(test_resolved_case_broadcast
    tests/mpi/test_resolved_case_broadcast.cpp)
  target_include_directories(test_resolved_case_broadcast PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_resolved_case_broadcast PRIVATE
    hundun_application hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(resolved_broadcast_test
      "test_resolved_case_broadcast_${ranks}_rank")
    add_test(NAME "${resolved_broadcast_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_resolved_case_broadcast>"
              full ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${resolved_broadcast_test}" PROPERTIES
      LABELS mpi TIMEOUT 45 PROCESSORS "${ranks}")
  endforeach()
  add_test(NAME test_resolved_case_broadcast_errors
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_resolved_case_broadcast>"
            errors ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_resolved_case_broadcast_errors PROPERTIES
    LABELS mpi TIMEOUT 45 PROCESSORS 4)
  add_test(NAME test_resolved_case_broadcast_finalized
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_resolved_case_broadcast>"
            finalized ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_resolved_case_broadcast_finalized PROPERTIES
    LABELS mpi TIMEOUT 20 PROCESSORS 1)

  add_test(NAME test_stage2_task2_dispatch
    COMMAND bash
            "${PROJECT_SOURCE_DIR}/tests/acceptance/stage2_task2_dispatch.sh"
            "$<TARGET_FILE:hundun>")
  set_tests_properties(test_stage2_task2_dispatch PROPERTIES
    LABELS acceptance TIMEOUT 30 PROCESSORS 1)

  add_executable(compare_scalar_vtk
    tests/acceptance/compare_scalar_vtk.cpp)
  target_link_libraries(compare_scalar_vtk PRIVATE
    hundun_options hundun_warnings)

  add_test(NAME acceptance_complete_gate
    COMMAND "${CMAKE_COMMAND}"
            -DHUNDUN_PROVENANCE_CASE=verify_acceptance_complete_gate
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  set_tests_properties(acceptance_complete_gate PROPERTIES
    LABELS acceptance
    TIMEOUT 10)

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

  add_test(NAME source_policy_clean
    COMMAND "${CMAKE_COMMAND}" -DHUNDUN_PROVENANCE_CASE=source_clean
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  string(CONCAT source_policy_script_runtime "py" "thon")
  set(source_policy_cases
    "${source_policy_script_runtime}_file"
    "${source_policy_script_runtime}_header"
    "${source_policy_script_runtime}_package"
    "${source_policy_script_runtime}_codegen"
    fetch external download missing tool_a tool_a3 tool_b vendor_control
    make_control make_at_control make_plus_control make_combined_control
    vendor_configure vendor_conditional vendor_path_retrieval)
  list(APPEND source_policy_cases
    vendor_quoted_path_retrieval
    vendor_quoted_action_retrieval
    vendor_quoted_bare_retrieval
    vendor_bracket_argument_retrieval
    yaml_inline_run
    yaml_block_run
    yaml_folded_retrieval)
  add_test(NAME source_policy_cmake_message_prose
    COMMAND "${CMAKE_COMMAND}"
            -DHUNDUN_PROVENANCE_CASE=source_cmake_message_prose
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  add_test(NAME source_policy_cmake_bracket_comment_prose
    COMMAND "${CMAKE_COMMAND}"
            -DHUNDUN_PROVENANCE_CASE=source_cmake_bracket_comment_prose
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  add_test(NAME source_policy_cmake_bracket_message_prose
    COMMAND "${CMAKE_COMMAND}"
            -DHUNDUN_PROVENANCE_CASE=source_cmake_bracket_message_prose
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  add_test(NAME source_policy_yaml_notes_prose
    COMMAND "${CMAKE_COMMAND}"
            -DHUNDUN_PROVENANCE_CASE=source_yaml_notes_prose
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  add_test(NAME source_policy_yaml_literal_split_prose
    COMMAND "${CMAKE_COMMAND}"
            -DHUNDUN_PROVENANCE_CASE=source_yaml_literal_split_prose
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  set(source_policy_tests
    source_policy_clean
    source_policy_cmake_message_prose
    source_policy_cmake_bracket_comment_prose
    source_policy_cmake_bracket_message_prose
    source_policy_yaml_notes_prose
    source_policy_yaml_literal_split_prose)
  foreach(source_policy_case IN LISTS source_policy_cases)
    add_test(NAME "source_policy_${source_policy_case}"
      COMMAND "${CMAKE_COMMAND}"
              "-DHUNDUN_PROVENANCE_CASE=verify_source_${source_policy_case}"
              -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
    list(APPEND source_policy_tests "source_policy_${source_policy_case}")
  endforeach()

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
  set_tests_properties(${source_policy_tests}
    PROPERTIES
      LABELS "provenance;source_policy"
      RESOURCE_LOCK hundun_provenance_fixture_tree)
  set_tests_properties(
    provenance_forbidden_fortran
    provenance_forbidden_token
    provenance_forbidden_uppercase_suffix
    PROPERTIES WILL_FAIL TRUE)
endif()
