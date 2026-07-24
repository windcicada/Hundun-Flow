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
find_package(Threads REQUIRED)

add_library(hundun_options INTERFACE)
target_compile_features(hundun_options INTERFACE cxx_std_17)
hundun_configure_sanitizers(hundun_options)

add_library(hundun_warnings INTERFACE)
hundun_configure_warnings(hundun_warnings)

add_library(hundun_execution STATIC
  execution/src/execution.cpp)
target_include_directories(hundun_execution PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/execution/include>"
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/runtime/include>")
target_link_libraries(hundun_execution
  PRIVATE hundun_warnings
  PUBLIC hundun_options Threads::Threads)

add_library(hundun_diagnostics STATIC
  diagnostics/src/performance_artifact.cpp)
target_include_directories(hundun_diagnostics PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/diagnostics/include>")
target_link_libraries(hundun_diagnostics
  PRIVATE hundun_warnings
  PUBLIC hundun_options)

add_library(hundun_diagnostics_core STATIC
  diagnostics/src/structured_diagnostics.cpp)
target_include_directories(hundun_diagnostics_core PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/diagnostics/include>")
target_link_libraries(hundun_diagnostics_core
  PRIVATE hundun_warnings
  PUBLIC hundun_options)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_diagnostics_core PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
endif()

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
  runtime/src/checkpoint_v2_protocol.cpp
  runtime/src/mpi_error.cpp
  runtime/src/mpi_environment.cpp
  runtime/src/mpi_context.cpp
  runtime/src/collective_status.cpp
  runtime/src/process_grid.cpp
  runtime/src/structured_decomposition.cpp
  runtime/src/field_registry.cpp
  runtime/src/field_access_plan.cpp
  runtime/src/field_storage.cpp
  runtime/src/halo_detail.cpp
  runtime/src/exchange_plan.cpp
  runtime/src/halo_exchange.cpp)
target_include_directories(hundun_runtime PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/runtime/include>")
target_link_libraries(hundun_runtime
  PRIVATE hundun_warnings
  PUBLIC hundun_options MPI::MPI_CXX)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_runtime PRIVATE
    HUNDUN_RUNTIME_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_application STATIC
  applications/hundun/case_config_broadcast.cpp)
target_include_directories(hundun_application
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
  PRIVATE "${PROJECT_SOURCE_DIR}/runtime/src")
target_link_libraries(hundun_application
  PRIVATE hundun_runtime hundun_warnings
  PUBLIC hundun_config hundun_options MPI::MPI_CXX)

add_library(hundun_mesh STATIC
  mesh/src/mesh_geometry.cpp
  mesh/src/mesh_topology.cpp
  mesh/src/uniform_structured_mesh.cpp)
target_include_directories(hundun_mesh PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/mesh/include>")
target_link_libraries(hundun_mesh
  PRIVATE hundun_warnings
  PUBLIC hundun_options hundun_runtime)

add_library(hundun_boundary STATIC
  boundary/src/basic_boundary.cpp)
target_include_directories(hundun_boundary PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/boundary/include>")
target_link_libraries(hundun_boundary
  PRIVATE hundun_warnings
  PUBLIC hundun_config hundun_mesh hundun_runtime hundun_options)

add_library(hundun_linear STATIC
  linear/src/bicgstab.cpp
  linear/src/conjugate_gradient.cpp
  linear/src/ghosted_vector.cpp
  linear/src/ghosted_vector_halo.cpp
  linear/src/ghosted_vector_halo_detail.cpp
  linear/src/preconditioners.cpp
  linear/src/solve_collective_detail.cpp
  linear/src/vector_ops.cpp)
target_include_directories(hundun_linear
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/linear/include>"
  PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/runtime/src")
target_link_libraries(hundun_linear
  PRIVATE hundun_warnings
  PUBLIC
    hundun_execution
    hundun_mesh
    hundun_runtime
    hundun_options
    MPI::MPI_CXX)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_linear PUBLIC
    HUNDUN_LINEAR_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_finite_volume STATIC
  finite_volume/src/cell_centered_fvm.cpp
  finite_volume/src/matrix_free_poisson.cpp
  finite_volume/src/poisson_boundary_adapter.cpp)
target_include_directories(hundun_finite_volume PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/finite_volume/include>"
  PRIVATE "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/runtime/src>")
target_link_libraries(hundun_finite_volume
  PRIVATE hundun_warnings
  PUBLIC
    hundun_boundary
    hundun_linear
    hundun_mesh
    hundun_runtime
    hundun_execution
    hundun_options)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_finite_volume PUBLIC
    HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_flow STATIC
  flow/src/checkpoint_v2.cpp
  flow/src/adaptive_time_control.cpp
  flow/src/constant_density_piso.cpp
  flow/src/flow_state.cpp
  flow/src/ideal_gas_closure.cpp
  flow/src/ideal_gas_piso.cpp
  flow/src/material_density_piso.cpp
  flow/src/material_density_transport.cpp
  flow/src/momentum_predictor.cpp)
target_include_directories(hundun_flow PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/flow/include>")
target_include_directories(hundun_flow PRIVATE
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/runtime/src>")
target_link_libraries(hundun_flow
  PRIVATE hundun_warnings
  PUBLIC
    hundun_finite_volume
    hundun_boundary
    hundun_linear
    hundun_mesh
    hundun_runtime
    hundun_execution
    hundun_options)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_flow PUBLIC
    HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_material_diagnostics STATIC
  diagnostics/src/ideal_gas_closure_diagnostics.cpp
  diagnostics/src/material_density_piso_diagnostics.cpp
  diagnostics/src/material_density_transport_diagnostics.cpp
  diagnostics/src/time_control_diagnostics.cpp)
target_include_directories(hundun_material_diagnostics PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/diagnostics/include>")
target_include_directories(hundun_material_diagnostics PRIVATE
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/flow/src>")
target_link_libraries(hundun_material_diagnostics
  PRIVATE hundun_warnings
  PUBLIC
    hundun_diagnostics_core
    hundun_flow
    hundun_runtime
    hundun_options
    MPI::MPI_CXX)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_material_diagnostics PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_checkpoint_diagnostics STATIC
  diagnostics/src/checkpoint_v2_diagnostics.cpp)
target_include_directories(hundun_checkpoint_diagnostics PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/diagnostics/include>")
target_include_directories(hundun_checkpoint_diagnostics PRIVATE
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/flow/src>"
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/runtime/src>")
target_link_libraries(hundun_checkpoint_diagnostics
  PRIVATE hundun_warnings
  PUBLIC
    hundun_diagnostics_core
    hundun_flow
    hundun_runtime
    hundun_options
    MPI::MPI_CXX)

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

  hundun_add_test(test_stage1_dispatch_contract
    tests/unit/test_stage1_dispatch_contract.cpp)
  target_include_directories(test_stage1_dispatch_contract PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_stage1_dispatch_contract PRIVATE
    hundun_application)

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

  hundun_add_test(test_checkpoint_v2_protocol
    tests/unit/test_checkpoint_v2_protocol.cpp)
  target_include_directories(test_checkpoint_v2_protocol PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/runtime/src")
  target_link_libraries(test_checkpoint_v2_protocol PRIVATE hundun_runtime)
  target_compile_definitions(test_checkpoint_v2_protocol PRIVATE
    HUNDUN_RUNTIME_ENABLE_TEST_ACCESS=1)
  set_tests_properties(test_checkpoint_v2_protocol PROPERTIES
    LABELS "unit;task23;task23-fast" TIMEOUT 30)
  add_test(NAME test_checkpoint_v2_protocol_acceptance
    COMMAND "$<TARGET_FILE:test_checkpoint_v2_protocol>")
  set_tests_properties(test_checkpoint_v2_protocol_acceptance PROPERTIES
    LABELS "unit;task23;task23-acceptance;task23-cluster-a" TIMEOUT 60)

  hundun_add_test(test_checkpoint_v2_epoch
    tests/unit/test_checkpoint_v2_epoch.cpp)
  target_include_directories(test_checkpoint_v2_epoch PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_checkpoint_v2_epoch PRIVATE hundun_runtime)
  set_tests_properties(test_checkpoint_v2_epoch PROPERTIES
    LABELS "unit;task23;task23-fast" TIMEOUT 30)

  hundun_add_test(test_checkpoint_v2_header_contract
    tests/unit/test_checkpoint_v2_header_contract.cpp)
  target_include_directories(test_checkpoint_v2_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/flow/src")
  target_link_libraries(test_checkpoint_v2_header_contract PRIVATE
    hundun_checkpoint_diagnostics hundun_flow)
  set_tests_properties(test_checkpoint_v2_header_contract PROPERTIES
    LABELS "unit;contract;task23;task23-fast" TIMEOUT 30)
  add_test(NAME test_task23_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task23_source_policy.cmake")
  set_tests_properties(test_task23_source_policy PROPERTIES
    LABELS "unit;source_policy;task23;task23-fast" TIMEOUT 30)
  add_test(NAME test_task23_test_access_contract
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task23_test_access_contract.cmake")
  set_tests_properties(test_task23_test_access_contract PROPERTIES
    LABELS "unit;source_policy;task23;task23-fast" TIMEOUT 30)

  hundun_add_test(test_field_access tests/unit/test_field_access.cpp)
  target_include_directories(test_field_access PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_field_access PRIVATE hundun_runtime)

  hundun_add_test(test_field_kernel_view
    tests/unit/test_field_kernel_view.cpp)
  target_include_directories(test_field_kernel_view PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_field_kernel_view PRIVATE hundun_runtime)

  hundun_add_test(test_field_view_header_contract
    tests/unit/test_field_view_header_contract.cpp)
  target_include_directories(test_field_view_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/runtime/include")

  hundun_add_test(test_performance_artifact
    tests/unit/test_performance_artifact.cpp)
  target_include_directories(test_performance_artifact PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/diagnostics/include")
  target_link_libraries(test_performance_artifact PRIVATE
    hundun_diagnostics yyjson)

  hundun_add_test(test_structured_diagnostics_header_contract
    tests/unit/test_structured_diagnostics_header_contract.cpp)
  target_include_directories(test_structured_diagnostics_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/diagnostics/include")
  set_tests_properties(test_structured_diagnostics_header_contract PROPERTIES
    LABELS "unit;task19" TIMEOUT 30)

  hundun_add_test(test_structured_diagnostics
    tests/unit/test_structured_diagnostics.cpp)
  target_include_directories(test_structured_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/diagnostics/include")
  target_link_libraries(test_structured_diagnostics PRIVATE
    hundun_diagnostics_core)
  target_compile_definitions(test_structured_diagnostics PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
  set_tests_properties(test_structured_diagnostics PROPERTIES
    LABELS "unit;task19" TIMEOUT 30)

  hundun_add_test(test_material_density_transport_header_contract
    tests/unit/test_material_density_transport_header_contract.cpp)
  target_include_directories(test_material_density_transport_header_contract
    PRIVATE "${PROJECT_SOURCE_DIR}/diagnostics/include")
  target_link_libraries(test_material_density_transport_header_contract PRIVATE
    hundun_flow)
  set_tests_properties(test_material_density_transport_header_contract PROPERTIES
    LABELS "unit;task19" TIMEOUT 30)

  hundun_add_test(test_material_density_piso_header_contract
    tests/unit/test_material_density_piso_header_contract.cpp)
  target_include_directories(test_material_density_piso_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/diagnostics/include")
  target_link_libraries(test_material_density_piso_header_contract PRIVATE
    hundun_material_diagnostics hundun_flow)
  set_tests_properties(test_material_density_piso_header_contract PROPERTIES
    LABELS "unit;task20;task20-fast" TIMEOUT 30)

  hundun_add_test(test_ideal_gas_header_contract
    tests/unit/test_ideal_gas_header_contract.cpp)
  target_include_directories(test_ideal_gas_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/diagnostics/include")
  target_link_libraries(test_ideal_gas_header_contract PRIVATE
    hundun_material_diagnostics hundun_flow)
  set_tests_properties(test_ideal_gas_header_contract PROPERTIES
    LABELS "unit;task21;task21-fast" TIMEOUT 30)

  hundun_add_test(test_mpi_operation_error
    tests/unit/test_mpi_operation_error.cpp)
  target_include_directories(test_mpi_operation_error PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_mpi_operation_error PRIVATE
    hundun_runtime MPI::MPI_CXX)
  set_tests_properties(test_mpi_operation_error PROPERTIES
    LABELS "unit;task19" TIMEOUT 30)

  string(REPLACE ";" "|" task19_mpi_include_dirs
         "${MPI_CXX_INCLUDE_DIRS}")
  string(REPLACE ";" "|" task19_mpi_libraries "${MPI_CXX_LIBRARIES}")
  set(task19_sanitizer_flags)
  if(HUNDUN_ENABLE_ASAN)
    list(APPEND task19_sanitizer_flags -fsanitize=address
         -fno-omit-frame-pointer)
  endif()
  if(HUNDUN_ENABLE_UBSAN)
    list(APPEND task19_sanitizer_flags -fsanitize=undefined
         -fno-omit-frame-pointer)
  endif()
  string(REPLACE ";" "|" task19_sanitizer_flags_argument
         "${task19_sanitizer_flags}")
  add_test(NAME test_task19_tests_off_header_contract
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            "-DHUNDUN_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
            "-DHUNDUN_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
            "-DHUNDUN_MPI_INCLUDE_DIRS=${task19_mpi_include_dirs}"
            "-DHUNDUN_MPI_LIBRARIES=${task19_mpi_libraries}"
            "-DHUNDUN_MPI_LINK_FLAGS=${MPI_CXX_LINK_FLAGS}"
            "-DHUNDUN_SANITIZER_FLAGS=${task19_sanitizer_flags_argument}"
            "-DHUNDUN_RUNTIME_LIBRARY=$<TARGET_FILE:hundun_runtime>"
            "-DHUNDUN_CONTRACT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task19_test_access_contract.cmake")
  set_tests_properties(test_task19_tests_off_header_contract PROPERTIES
    LABELS "unit;source_policy;task19" TIMEOUT 30)
  add_test(NAME test_task19_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task19_source_policy.cmake")
  set_tests_properties(test_task19_source_policy PROPERTIES
    LABELS "unit;source_policy;task19" TIMEOUT 30)

  add_test(NAME test_task20_tests_off_header_contract
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            "-DHUNDUN_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
            "-DHUNDUN_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
            "-DHUNDUN_MPI_INCLUDE_DIRS=${task19_mpi_include_dirs}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task20_test_access_contract.cmake")
  set_tests_properties(test_task20_tests_off_header_contract PROPERTIES
    LABELS "unit;source_policy;task20;task20-fast" TIMEOUT 30)
  add_test(NAME test_task20_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task20_source_policy.cmake")
  set_tests_properties(test_task20_source_policy PROPERTIES
    LABELS "unit;source_policy;task20;task20-fast" TIMEOUT 30)

  add_test(NAME test_task21_tests_off_header_contract
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            "-DHUNDUN_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
            "-DHUNDUN_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
            "-DHUNDUN_MPI_INCLUDE_DIRS=${task19_mpi_include_dirs}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task21_test_access_contract.cmake")
  set_tests_properties(test_task21_tests_off_header_contract PROPERTIES
    LABELS "unit;source_policy;task21;task21-fast" TIMEOUT 30)
  add_test(NAME test_task21_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task21_source_policy.cmake")
  set_tests_properties(test_task21_source_policy PROPERTIES
    LABELS "unit;source_policy;task21;task21-fast" TIMEOUT 30)

  hundun_add_test(test_execution tests/unit/test_execution.cpp)
  target_include_directories(test_execution PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_execution PRIVATE hundun_execution)
  set_tests_properties(test_execution PROPERTIES TIMEOUT 30)

  hundun_add_test(test_ghosted_vector
    tests/unit/test_ghosted_vector.cpp)
  target_include_directories(test_ghosted_vector PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_ghosted_vector PRIVATE hundun_linear)
  set_tests_properties(test_ghosted_vector PROPERTIES TIMEOUT 30)

  hundun_add_test(test_ghosted_vector_header_contract
    tests/unit/test_ghosted_vector_header_contract.cpp)
  target_include_directories(test_ghosted_vector_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/linear/include")
  target_link_libraries(test_ghosted_vector_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_ghosted_vector_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_vector_ops_header_contract
    tests/unit/test_vector_ops_header_contract.cpp)
  target_link_libraries(test_vector_ops_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_vector_ops_header_contract PROPERTIES TIMEOUT 30)

  hundun_add_test(test_vector_ops
    tests/unit/test_vector_ops.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_vector_ops PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_vector_ops PRIVATE hundun_linear)
  set_tests_properties(test_vector_ops PROPERTIES TIMEOUT 30)

  hundun_add_test(test_linear_system_header_contract
    tests/unit/test_linear_system_header_contract.cpp)
  target_link_libraries(test_linear_system_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_linear_system_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_preconditioners
    tests/unit/test_preconditioners.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_preconditioners PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_preconditioners PRIVATE hundun_linear)
  set_tests_properties(test_preconditioners PROPERTIES TIMEOUT 30)

  hundun_add_test(test_conjugate_gradient_header_contract
    tests/unit/test_conjugate_gradient_header_contract.cpp)
  target_link_libraries(test_conjugate_gradient_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_conjugate_gradient_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_bicgstab_header_contract
    tests/unit/test_bicgstab_header_contract.cpp)
  target_link_libraries(test_bicgstab_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_bicgstab_header_contract PROPERTIES TIMEOUT 30)

  hundun_add_test(test_matrix_free_poisson_header_contract
    tests/unit/test_matrix_free_poisson_header_contract.cpp)
  target_link_libraries(test_matrix_free_poisson_header_contract PRIVATE
    hundun_finite_volume)
  set_tests_properties(test_matrix_free_poisson_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_basic_boundary_header_contract
    tests/unit/test_basic_boundary_header_contract.cpp)
  target_include_directories(test_basic_boundary_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_compile_definitions(test_basic_boundary_header_contract PRIVATE
    HUNDUN_SOURCE_ROOT="${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_basic_boundary_header_contract PRIVATE
    hundun_config hundun_finite_volume hundun_mesh hundun_runtime)
  set_tests_properties(test_basic_boundary_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_cell_centered_fvm_header_contract
    tests/unit/test_cell_centered_fvm_header_contract.cpp)
  target_link_libraries(test_cell_centered_fvm_header_contract PRIVATE
    hundun_finite_volume)
  set_tests_properties(test_cell_centered_fvm_header_contract PROPERTIES
    TIMEOUT 30)
  add_test(NAME test_cell_centered_fvm_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task16_source_policy.cmake")
  set_tests_properties(test_cell_centered_fvm_source_policy PROPERTIES
    LABELS "unit;source_policy"
    TIMEOUT 30)

  hundun_add_test(test_momentum_predictor_header_contract
    tests/unit/test_momentum_predictor_header_contract.cpp)
  target_link_libraries(test_momentum_predictor_header_contract PRIVATE
    hundun_flow)
  set_tests_properties(test_momentum_predictor_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_task18_header_contract
    tests/unit/test_task18_header_contract.cpp)
  target_link_libraries(test_task18_header_contract PRIVATE hundun_flow)
  set_tests_properties(test_task18_header_contract PROPERTIES TIMEOUT 30)
  string(REPLACE ";" "|" task18_mpi_include_dirs
         "${MPI_CXX_INCLUDE_DIRS}")
  add_test(NAME test_task18_tests_off_header_contract
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            "-DHUNDUN_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
            "-DHUNDUN_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
            "-DHUNDUN_MPI_INCLUDE_DIRS=${task18_mpi_include_dirs}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task18_test_access_contract.cmake")
  set_tests_properties(test_task18_tests_off_header_contract PROPERTIES
    LABELS "unit;source_policy;task18"
    TIMEOUT 30)
  add_test(NAME test_task18_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task18_source_policy.cmake")
  set_tests_properties(test_task18_source_policy PROPERTIES
    LABELS "unit;source_policy;task18"
    TIMEOUT 30)

  hundun_add_test(test_flow_state
    tests/unit/test_flow_state.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_flow_state PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_flow_state PRIVATE hundun_flow)
  set_tests_properties(test_flow_state PROPERTIES TIMEOUT 30)

  hundun_add_test(test_adaptive_time_control_header_contract
    tests/unit/test_adaptive_time_control_header_contract.cpp)
  target_include_directories(test_adaptive_time_control_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_adaptive_time_control_header_contract PRIVATE
    hundun_flow)
  set_tests_properties(test_adaptive_time_control_header_contract PROPERTIES
    LABELS "unit;task22;task22-fast" TIMEOUT 30)

  hundun_add_test(test_adaptive_time_control
    tests/unit/test_adaptive_time_control.cpp)
  target_include_directories(test_adaptive_time_control PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_adaptive_time_control PRIVATE hundun_flow)
  set_tests_properties(test_adaptive_time_control PROPERTIES
    LABELS "unit;task22;task22-fast" TIMEOUT 30)

  add_executable(test_adaptive_time_control_mpi
    tests/mpi/test_adaptive_time_control.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_adaptive_time_control_mpi PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_compile_definitions(test_adaptive_time_control_mpi PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS)
  target_link_libraries(test_adaptive_time_control_mpi PRIVATE
    hundun_material_diagnostics hundun_flow hundun_options hundun_warnings
    MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    foreach(mode IN ITEMS fast acceptance)
      set(task22_controller_test
          "test_adaptive_time_control_${ranks}_rank_${mode}")
      add_test(NAME "${task22_controller_test}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_adaptive_time_control_mpi>" "${mode}"
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties("${task22_controller_test}" PROPERTIES
        LABELS "mpi;numerical;task22;task22-${mode}"
        TIMEOUT 240 PROCESSORS "${ranks}" RUN_SERIAL TRUE)
    endforeach()
  endforeach()
  foreach(selector_case IN ITEMS missing unknown extra)
    if(selector_case STREQUAL "missing")
      set(selector_args)
    elseif(selector_case STREQUAL "unknown")
      set(selector_args "not-a-task22-mode")
    else()
      set(selector_args "fast" "unexpected-extra-argument")
    endif()
    add_test(NAME "test_adaptive_time_control_selector_${selector_case}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_adaptive_time_control_mpi>"
              ${selector_args} ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_adaptive_time_control_selector_${selector_case}" PROPERTIES
      LABELS "mpi;contract;task22;task22-fast"
      TIMEOUT 30 PROCESSORS 1 RUN_SERIAL TRUE WILL_FAIL TRUE)
  endforeach()
  add_executable(test_adaptive_time_control_diagnostics_mpi
    tests/mpi/test_adaptive_time_control_diagnostics.cpp)
  target_include_directories(test_adaptive_time_control_diagnostics_mpi PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_compile_definitions(test_adaptive_time_control_diagnostics_mpi PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS)
  target_link_libraries(test_adaptive_time_control_diagnostics_mpi PRIVATE
    hundun_material_diagnostics hundun_flow hundun_options hundun_warnings
    MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    foreach(mode IN ITEMS fast acceptance)
      set(task22_diagnostics_test
          "test_adaptive_time_control_diagnostics_${ranks}_rank_${mode}")
      add_test(NAME "${task22_diagnostics_test}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_adaptive_time_control_diagnostics_mpi>"
                "${mode}" ${MPIEXEC_POSTFLAGS})
      set_tests_properties("${task22_diagnostics_test}" PROPERTIES
        LABELS "mpi;diagnostics;task22;task22-${mode}"
        TIMEOUT 240 PROCESSORS "${ranks}" RUN_SERIAL TRUE)
    endforeach()
  endforeach()
  foreach(selector_case IN ITEMS missing unknown extra)
    if(selector_case STREQUAL "missing")
      set(selector_args)
    elseif(selector_case STREQUAL "unknown")
      set(selector_args "not-a-task22-mode")
    else()
      set(selector_args "acceptance" "unexpected-extra-argument")
    endif()
    add_test(
      NAME "test_adaptive_time_control_diagnostics_selector_${selector_case}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_adaptive_time_control_diagnostics_mpi>"
              ${selector_args} ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_adaptive_time_control_diagnostics_selector_${selector_case}"
      PROPERTIES LABELS "mpi;contract;diagnostics;task22;task22-fast"
      TIMEOUT 30 PROCESSORS 1 RUN_SERIAL TRUE WILL_FAIL TRUE)
  endforeach()
  add_test(NAME test_task22_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task22_source_policy.cmake")
  set_tests_properties(test_task22_source_policy PROPERTIES
    LABELS "unit;source_policy;task22" TIMEOUT 30)
  add_test(NAME test_task22_tests_off_path_contract
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SETUP_SCRIPT=${PROJECT_SOURCE_DIR}/tests/cmake/task22_tests_off_setup.cmake"
            "-DHUNDUN_TEST_BINARY_ROOT=${CMAKE_BINARY_DIR}"
            "-DHUNDUN_GENERATOR=${CMAKE_GENERATOR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task22_tests_off_path_contract.cmake")
  set_tests_properties(test_task22_tests_off_path_contract PROPERTIES
    LABELS "unit;contract;source_policy;task22"
    RUN_SERIAL TRUE
    TIMEOUT 60)
  string(REPLACE ";" "|" task22_mpi_include_dirs
         "${MPI_CXX_INCLUDE_DIRS}")
  set(task22_tests_off_binary_dir
      "${CMAKE_BINARY_DIR}/task22-tests-off")
  set(task22_tests_off_archive_dir
      "${task22_tests_off_binary_dir}/archives")
  add_test(NAME test_task22_tests_off_setup
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            "-DHUNDUN_OUTER_BINARY_DIR=${CMAKE_BINARY_DIR}"
            "-DHUNDUN_NESTED_BINARY_DIR=${task22_tests_off_binary_dir}"
            "-DHUNDUN_ARCHIVE_DIR=${task22_tests_off_archive_dir}"
            "-DHUNDUN_NESTED_CONFIG=Release"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task22_tests_off_setup.cmake")
  set_tests_properties(test_task22_tests_off_setup PROPERTIES
    FIXTURES_SETUP task22_tests_off_archives
    LABELS "unit;source_policy;task22"
    RUN_SERIAL TRUE
    TIMEOUT 600)
  add_test(NAME test_task22_tests_off_header_contract
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            "-DHUNDUN_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
            "-DHUNDUN_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
            "-DHUNDUN_MPI_INCLUDE_DIRS=${task22_mpi_include_dirs}"
            "-DHUNDUN_TESTS_OFF_FLOW=${task22_tests_off_archive_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}hundun_flow${CMAKE_STATIC_LIBRARY_SUFFIX}"
            "-DHUNDUN_TESTS_OFF_LINEAR=${task22_tests_off_archive_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}hundun_linear${CMAKE_STATIC_LIBRARY_SUFFIX}"
            "-DHUNDUN_TESTS_OFF_DIAGNOSTICS=${task22_tests_off_archive_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}hundun_material_diagnostics${CMAKE_STATIC_LIBRARY_SUFFIX}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task22_test_access_contract.cmake")
  set_tests_properties(test_task22_tests_off_header_contract PROPERTIES
    FIXTURES_REQUIRED task22_tests_off_archives
    LABELS "unit;source_policy;task22"
    RUN_SERIAL TRUE
    TIMEOUT 60)

  hundun_add_test(test_momentum_predictor
    tests/unit/test_momentum_predictor.cpp)
  target_include_directories(test_momentum_predictor PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_momentum_predictor PRIVATE hundun_flow)
  set_tests_properties(test_momentum_predictor PROPERTIES TIMEOUT 30)
  add_test(NAME test_momentum_predictor_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task17_source_policy.cmake")
  set_tests_properties(test_momentum_predictor_source_policy PROPERTIES
    LABELS "unit;source_policy"
    TIMEOUT 30)

  add_executable(test_time_consistent_face_velocity
    tests/mpi/test_time_consistent_face_velocity.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_time_consistent_face_velocity PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_time_consistent_face_velocity PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(time_consistent_face_test
        "test_time_consistent_face_velocity_${ranks}_rank")
    add_test(NAME "${time_consistent_face_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_time_consistent_face_velocity>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${time_consistent_face_test}" PROPERTIES
      LABELS "mpi;numerical"
      TIMEOUT 120
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_momentum_predictor_collective
    tests/mpi/test_momentum_predictor_collective.cpp)
  target_include_directories(test_momentum_predictor_collective PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_momentum_predictor_collective PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 2 4)
    set(momentum_predictor_collective_test
        "test_momentum_predictor_collective_${ranks}_rank")
    add_test(NAME "${momentum_predictor_collective_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_momentum_predictor_collective>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${momentum_predictor_collective_test}" PROPERTIES
      LABELS "mpi"
      TIMEOUT 60
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_fixed_step_piso
    tests/mpi/test_fixed_step_piso.cpp)
  target_include_directories(test_fixed_step_piso PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_fixed_step_piso PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(fixed_step_piso_test "test_fixed_step_piso_${ranks}_rank")
    add_test(NAME "${fixed_step_piso_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_fixed_step_piso>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${fixed_step_piso_test}" PROPERTIES
      LABELS "mpi;numerical;task18"
      TIMEOUT 120
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_material_density_transport
    tests/mpi/test_material_density_transport.cpp)
  target_include_directories(test_material_density_transport PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_material_density_transport PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(material_transport_test
        "test_material_density_transport_${ranks}_rank")
    add_test(NAME "${material_transport_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_material_density_transport>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${material_transport_test}" PROPERTIES
      LABELS "mpi;numerical;task19" TIMEOUT 120 PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_material_density_piso
    tests/mpi/test_material_density_piso.cpp)
  target_include_directories(test_material_density_piso PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_material_density_piso PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  add_executable(test_variable_density_vortex
    tests/mpi/test_variable_density_vortex.cpp)
  target_include_directories(test_variable_density_vortex PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_variable_density_vortex PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  add_executable(test_material_density_piso_diagnostics
    tests/mpi/test_material_density_piso_diagnostics.cpp)
  target_include_directories(test_material_density_piso_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_material_density_piso_diagnostics PRIVATE
    hundun_material_diagnostics hundun_flow hundun_options hundun_warnings
    MPI::MPI_CXX)
  target_compile_definitions(test_material_density_piso_diagnostics PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
  add_executable(test_ideal_gas_closure
    tests/mpi/test_ideal_gas_closure.cpp)
  target_include_directories(test_ideal_gas_closure PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_ideal_gas_closure PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  add_executable(test_ideal_gas_piso
    tests/mpi/test_ideal_gas_piso.cpp)
  target_include_directories(test_ideal_gas_piso PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_ideal_gas_piso PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  add_executable(test_ideal_gas_diagnostics
    tests/mpi/test_ideal_gas_diagnostics.cpp)
  target_include_directories(test_ideal_gas_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_ideal_gas_diagnostics PRIVATE
    hundun_material_diagnostics hundun_flow hundun_options hundun_warnings
    MPI::MPI_CXX)
  target_compile_definitions(test_ideal_gas_diagnostics PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
  add_executable(test_checkpoint_v2
    tests/mpi/test_checkpoint_v2.cpp)
  target_include_directories(test_checkpoint_v2 PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/flow/src"
    "${PROJECT_SOURCE_DIR}/runtime/src")
  target_link_libraries(test_checkpoint_v2 PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  target_compile_definitions(test_checkpoint_v2 PRIVATE
    HUNDUN_RUNTIME_ENABLE_TEST_ACCESS=1)
  foreach(ranks IN ITEMS 1 2 4)
    set(task23_checkpoint_test "test_checkpoint_v2_${ranks}_rank")
    add_test(NAME "${task23_checkpoint_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_checkpoint_v2>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${task23_checkpoint_test}" PROPERTIES
      LABELS "mpi;numerical;task23;task23-fast"
      TIMEOUT 180 PROCESSORS "${ranks}" RUN_SERIAL TRUE)
    set(task23_checkpoint_acceptance
        "test_checkpoint_v2_acceptance_${ranks}_rank")
    add_test(NAME "${task23_checkpoint_acceptance}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_checkpoint_v2>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${task23_checkpoint_acceptance}" PROPERTIES
      LABELS
        "mpi;numerical;task23;task23-acceptance;task23-cluster-b;task23-cluster-c;task23-cluster-d"
      TIMEOUT 240 PROCESSORS "${ranks}" RUN_SERIAL TRUE)
  endforeach()
  add_executable(test_checkpoint_v2_diagnostics
    tests/mpi/test_checkpoint_v2_diagnostics.cpp)
  target_include_directories(test_checkpoint_v2_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_checkpoint_v2_diagnostics PRIVATE
    hundun_checkpoint_diagnostics hundun_flow hundun_options hundun_warnings
    MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(task23_diagnostics_test
        "test_checkpoint_v2_diagnostics_${ranks}_rank")
    add_test(NAME "${task23_diagnostics_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_checkpoint_v2_diagnostics>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${task23_diagnostics_test}" PROPERTIES
      LABELS "mpi;diagnostics;task23;task23-fast"
      TIMEOUT 120 PROCESSORS "${ranks}" RUN_SERIAL TRUE)
    set(task23_diagnostics_acceptance
        "test_checkpoint_v2_diagnostics_acceptance_${ranks}_rank")
    add_test(NAME "${task23_diagnostics_acceptance}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_checkpoint_v2_diagnostics>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${task23_diagnostics_acceptance}" PROPERTIES
      LABELS "mpi;diagnostics;task23;task23-acceptance;task23-cluster-e"
      TIMEOUT 180 PROCESSORS "${ranks}" RUN_SERIAL TRUE)
  endforeach()
  set(task20_vortex_reference
      "${CMAKE_CURRENT_BINARY_DIR}/task20-variable-density-vortex-reference.v1.bin")
  set(task21_closed_reference
      "${CMAKE_CURRENT_BINARY_DIR}/task21-ideal-gas-closed-reference.v1.bin")
  set(task21_open_reference
      "${CMAKE_CURRENT_BINARY_DIR}/task21-ideal-gas-open-reference.v1.bin")
  set(task21_diagnostic_reference
      "${CMAKE_CURRENT_BINARY_DIR}/task21-ideal-gas-diagnostic-reference.v1.bin")
  foreach(ranks IN ITEMS 1 2 4)
    foreach(target IN ITEMS material_density_piso
                            variable_density_vortex
                            material_density_piso_diagnostics)
      set(task20_test "test_${target}_${ranks}_rank")
      add_test(NAME "${task20_test}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_${target}>"
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties("${task20_test}" PROPERTIES
        LABELS "mpi;numerical;task20;task20-fast"
        TIMEOUT 240 PROCESSORS "${ranks}")
    endforeach()
    foreach(target IN ITEMS ideal_gas_closure ideal_gas_piso
                            ideal_gas_diagnostics)
      set(task21_test "test_${target}_${ranks}_rank")
      add_test(NAME "${task21_test}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_${target}>"
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties("${task21_test}" PROPERTIES
        LABELS "mpi;numerical;task21;task21-fast"
        TIMEOUT 300 PROCESSORS "${ranks}")
    endforeach()
    foreach(case_name IN ITEMS closed_heating open_plug)
      if(case_name STREQUAL "open_plug")
        set(task21_case_argument --open-plug)
      else()
        set(task21_case_argument --closed-heating)
      endif()
      set(task21_full "test_ideal_gas_${case_name}_full_${ranks}_rank")
      if(case_name STREQUAL "open_plug")
        set(task21_reference "${task21_open_reference}")
        set(task21_reference_fixture task21_ideal_gas_open_reference)
      else()
        set(task21_reference "${task21_closed_reference}")
        set(task21_reference_fixture task21_ideal_gas_closed_reference)
      endif()
      if(ranks EQUAL 1)
        set(task21_reference_arguments
            --reference-write "${task21_reference}")
      else()
        set(task21_reference_arguments
            --reference-read "${task21_reference}")
      endif()
      add_test(NAME "${task21_full}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_ideal_gas_piso>"
                "${task21_case_argument}" ${task21_reference_arguments}
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties("${task21_full}" PROPERTIES
        LABELS "mpi;numerical;acceptance;task21"
        TIMEOUT 1800 PROCESSORS "${ranks}")
      if(ranks EQUAL 1)
        set_tests_properties("${task21_full}" PROPERTIES
          FIXTURES_SETUP "${task21_reference_fixture}")
      else()
        set_tests_properties("${task21_full}" PROPERTIES
          FIXTURES_REQUIRED "${task21_reference_fixture}")
      endif()
    endforeach()
    set(task21_diagnostic_decomposition
        "test_ideal_gas_diagnostics_decomposition_${ranks}_rank")
    if(ranks EQUAL 1)
      set(task21_diagnostic_reference_arguments
          --reference-write "${task21_diagnostic_reference}")
    else()
      set(task21_diagnostic_reference_arguments
          --reference-read "${task21_diagnostic_reference}")
    endif()
    add_test(NAME "${task21_diagnostic_decomposition}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_ideal_gas_diagnostics>"
              ${task21_diagnostic_reference_arguments}
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${task21_diagnostic_decomposition}" PROPERTIES
      LABELS "mpi;numerical;acceptance;task21"
      TIMEOUT 300 PROCESSORS "${ranks}")
    if(ranks EQUAL 1)
      set_tests_properties("${task21_diagnostic_decomposition}" PROPERTIES
        FIXTURES_SETUP task21_ideal_gas_diagnostic_reference)
    else()
      set_tests_properties("${task21_diagnostic_decomposition}" PROPERTIES
        FIXTURES_REQUIRED task21_ideal_gas_diagnostic_reference)
    endif()
    set(task20_vortex_full "test_variable_density_vortex_full_${ranks}_rank")
    if(ranks EQUAL 1)
      set(task20_vortex_reference_arguments
          --reference-write "${task20_vortex_reference}")
    else()
      set(task20_vortex_reference_arguments
          --reference-read "${task20_vortex_reference}")
    endif()
    add_test(NAME "${task20_vortex_full}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_variable_density_vortex>"
              --full ${task20_vortex_reference_arguments}
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${task20_vortex_full}" PROPERTIES
      LABELS "mpi;numerical;task20"
      TIMEOUT 1800 PROCESSORS "${ranks}")
    if(ranks EQUAL 1)
      set_tests_properties("${task20_vortex_full}" PROPERTIES
        FIXTURES_SETUP task20_vortex_reference)
    else()
      set_tests_properties("${task20_vortex_full}" PROPERTIES
        FIXTURES_REQUIRED task20_vortex_reference)
    endif()
  endforeach()

  add_executable(test_material_density_wave
    tests/mpi/test_material_density_wave.cpp)
  target_include_directories(test_material_density_wave PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_material_density_wave PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(material_wave_test "test_material_density_wave_${ranks}_rank")
    add_test(NAME "${material_wave_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_material_density_wave>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${material_wave_test}" PROPERTIES
      LABELS "mpi;numerical;task19" TIMEOUT 240 PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_material_density_diagnostics
    tests/mpi/test_material_density_diagnostics.cpp)
  target_include_directories(test_material_density_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_material_density_diagnostics PRIVATE
    hundun_material_diagnostics hundun_flow hundun_options hundun_warnings
    MPI::MPI_CXX)
  target_compile_definitions(test_material_density_diagnostics PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1
    HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)
  foreach(ranks IN ITEMS 1 2 4)
    set(material_diagnostics_test
        "test_material_density_diagnostics_${ranks}_rank")
    add_test(NAME "${material_diagnostics_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_material_density_diagnostics>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${material_diagnostics_test}" PROPERTIES
      LABELS "mpi;diagnostics;task19" TIMEOUT 120 PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_piso_decomposition
    tests/mpi/test_piso_decomposition.cpp)
  target_include_directories(test_piso_decomposition PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_piso_decomposition PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  add_test(NAME test_piso_decomposition_1_2_4
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_piso_decomposition>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_piso_decomposition_1_2_4 PROPERTIES
    LABELS "mpi;numerical;task18"
    TIMEOUT 30
    PROCESSORS 4)

  add_executable(test_taylor_green_piso
    tests/mpi/test_taylor_green_piso.cpp)
  target_include_directories(test_taylor_green_piso PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_taylor_green_piso PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(taylor_green_piso_test "test_taylor_green_piso_${ranks}_rank")
    add_test(NAME "${taylor_green_piso_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_taylor_green_piso>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${taylor_green_piso_test}" PROPERTIES
      LABELS "mpi;numerical;task18"
      TIMEOUT 300
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_basic_boundary
    tests/mpi/test_basic_boundary.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_basic_boundary PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_basic_boundary PRIVATE
    hundun_config hundun_finite_volume hundun_mesh hundun_runtime
    hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(basic_boundary_test "test_basic_boundary_${ranks}_rank")
    add_test(NAME "${basic_boundary_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_basic_boundary>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${basic_boundary_test}" PROPERTIES
      LABELS "mpi;numerical"
      TIMEOUT 90
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_cell_centered_fvm
    tests/mpi/test_cell_centered_fvm.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_cell_centered_fvm PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_cell_centered_fvm PRIVATE
    hundun_finite_volume hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(cell_centered_fvm_test "test_cell_centered_fvm_${ranks}_rank")
    add_test(NAME "${cell_centered_fvm_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_cell_centered_fvm>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${cell_centered_fvm_test}" PROPERTIES
      LABELS "mpi;numerical"
      TIMEOUT 120
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_matrix_free_poisson
    tests/mpi/test_matrix_free_poisson.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_matrix_free_poisson PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_matrix_free_poisson PRIVATE
    hundun_finite_volume hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(matrix_free_poisson_test
        "test_matrix_free_poisson_${ranks}_rank")
    add_test(NAME "${matrix_free_poisson_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_matrix_free_poisson>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${matrix_free_poisson_test}" PROPERTIES
      LABELS "mpi;numerical"
      TIMEOUT 90
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_poisson_sine
    tests/numerical/test_poisson_sine.cpp)
  target_include_directories(test_poisson_sine PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_poisson_sine PRIVATE
    hundun_finite_volume hundun_options hundun_warnings MPI::MPI_CXX)
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    foreach(ranks IN ITEMS 1 2 4)
      set(poisson_sine_test "test_poisson_sine_${ranks}_rank")
      add_test(NAME "${poisson_sine_test}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_poisson_sine>"
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties("${poisson_sine_test}" PROPERTIES
        LABELS "mpi;numerical;acceptance"
        TIMEOUT 240
        PROCESSORS "${ranks}")
    endforeach()
  endif()

  add_executable(test_conjugate_gradient
    tests/mpi/test_conjugate_gradient.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_conjugate_gradient PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_conjugate_gradient PRIVATE
    hundun_linear hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(conjugate_gradient_test "test_conjugate_gradient_${ranks}_rank")
    add_test(NAME "${conjugate_gradient_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_conjugate_gradient>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${conjugate_gradient_test}" PROPERTIES
      LABELS mpi TIMEOUT 60 PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_bicgstab
    tests/mpi/test_bicgstab.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_bicgstab PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_bicgstab PRIVATE
    hundun_linear hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(bicgstab_test "test_bicgstab_${ranks}_rank")
    add_test(NAME "${bicgstab_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_bicgstab>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${bicgstab_test}" PROPERTIES
      LABELS mpi TIMEOUT 60 PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_vector_ops_mpi
    tests/mpi/test_vector_ops_mpi.cpp
    tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_vector_ops_mpi PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_vector_ops_mpi PRIVATE
    hundun_linear hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(vector_ops_mpi_test "test_vector_ops_mpi_${ranks}_rank")
    add_test(NAME "${vector_ops_mpi_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_vector_ops_mpi>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${vector_ops_mpi_test}" PROPERTIES
      LABELS mpi TIMEOUT 60 PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_ghosted_vector_halo
    tests/mpi/test_ghosted_vector_halo.cpp)
  target_include_directories(test_ghosted_vector_halo PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_ghosted_vector_halo PRIVATE
    hundun_linear hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(vector_halo_test "test_ghosted_vector_halo_${ranks}_rank")
    add_test(NAME "${vector_halo_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_ghosted_vector_halo>"
              full ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${vector_halo_test}" PROPERTIES
      LABELS mpi TIMEOUT 60 PROCESSORS "${ranks}")
  endforeach()
  foreach(mode IN ITEMS mismatch failure)
    set(vector_halo_mode_test "test_ghosted_vector_halo_${mode}")
    add_test(NAME "${vector_halo_mode_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_ghosted_vector_halo>"
              "${mode}" ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${vector_halo_mode_test}" PROPERTIES
      LABELS mpi TIMEOUT 60 PROCESSORS 4)
  endforeach()
  add_test(NAME test_ghosted_vector_halo_finalized_idle
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_ghosted_vector_halo>"
            finalized_idle ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_ghosted_vector_halo_finalized_idle PROPERTIES
    LABELS mpi TIMEOUT 30 PROCESSORS 1)

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

  add_executable(test_mesh_topology tests/mpi/test_mesh_topology.cpp)
  target_include_directories(test_mesh_topology PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_mesh_topology PRIVATE
    hundun_mesh hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(mesh_topology_test "test_mesh_topology_${ranks}_rank")
    add_test(NAME "${mesh_topology_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_mesh_topology>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${mesh_topology_test}" PROPERTIES
      LABELS mpi
      TIMEOUT 60
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_mesh_geometry tests/mpi/test_mesh_geometry.cpp)
  target_include_directories(test_mesh_geometry PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_mesh_geometry PRIVATE
    hundun_mesh hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(mesh_geometry_test "test_mesh_geometry_${ranks}_rank")
    add_test(NAME "${mesh_geometry_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_mesh_geometry>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${mesh_geometry_test}" PROPERTIES
      LABELS "mpi;numerical"
      TIMEOUT 90
      PROCESSORS "${ranks}")
  endforeach()
  add_test(NAME test_mesh_geometry_injected_sample_failure
    COMMAND "${CMAKE_COMMAND}"
            "-DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}"
            "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
            "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
            "-DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}"
            "-DTEST_EXECUTABLE=$<TARGET_FILE:test_mesh_geometry>"
            -P "${PROJECT_SOURCE_DIR}/tests/support/expected_geometry_sample_failure.cmake")
  set_tests_properties(test_mesh_geometry_injected_sample_failure PROPERTIES
    LABELS "mpi;numerical"
    TIMEOUT 10
    PROCESSORS 2
    WILL_FAIL TRUE)

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

    set(resolved_v1_mismatch_test
      "test_resolved_case_broadcast_v1_mismatch_${ranks}_rank")
    add_test(NAME "${resolved_v1_mismatch_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_resolved_case_broadcast>"
              v1_mismatch ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${resolved_v1_mismatch_test}" PROPERTIES
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
            "$<TARGET_FILE:hundun>" "${MPIEXEC_EXECUTABLE}")
  set_tests_properties(test_stage2_task2_dispatch PROPERTIES
    LABELS acceptance TIMEOUT 45 PROCESSORS 2)

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
