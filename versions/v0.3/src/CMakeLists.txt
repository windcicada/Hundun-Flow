# SPDX-License-Identifier: Apache-2.0

include(GNUInstallDirs)

target_compile_features(hundun_options INTERFACE cxx_std_17)
target_include_directories(hundun_options INTERFACE
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>")
hundun_configure_sanitizers(hundun_options)

add_library(hundun_warnings INTERFACE)
hundun_configure_warnings(hundun_warnings)

add_library(hundun_chemistry INTERFACE)
target_link_libraries(hundun_chemistry INTERFACE hundun_options)
if(HUNDUN_ENABLE_CANTERA)
  add_library(hundun_cantera_backend STATIC
    chem_cantera_backend.cpp)
  target_include_directories(hundun_cantera_backend
    PRIVATE "${PROJECT_SOURCE_DIR}/src"
    PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>")
  target_link_libraries(hundun_cantera_backend
    PRIVATE hundun_warnings hundun_third_party_cantera
    PUBLIC hundun_chemistry hundun_config hundun_options)
endif()

add_library(hundun_execution STATIC
  exec_execution.cpp)
target_include_directories(hundun_execution
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_execution
  PRIVATE hundun_warnings
  PUBLIC hundun_options Threads::Threads)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_execution PRIVATE
    HUNDUN_EXECUTION_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_diagnostics STATIC
  diag_performance_artifact.cpp
  diag_stage3_evidence_manifest.cpp
  diag_stage3_performance.cpp
  diag_performance_correctness.cpp)
target_include_directories(hundun_diagnostics
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_diagnostics
  PRIVATE hundun_warnings
  PUBLIC hundun_options)

add_library(hundun_diagnostics_core STATIC
  diag_structured.cpp
  diag_provider_registry.cpp
  diag_reacting.cpp)
target_include_directories(hundun_diagnostics_core
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_diagnostics_core
  PRIVATE hundun_warnings
  PUBLIC hundun_options)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_diagnostics_core PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
endif()

add_library(yyjson STATIC ${PROJECT_SOURCE_DIR}/third_party/yyjson/yyjson.c)
set_target_properties(yyjson PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(yyjson PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third_party/yyjson>")

add_library(hundun_config STATIC
  cfg_case_config_loader.cpp
  cfg_resolved_case_loader.cpp
  cfg_resolved_case_v3_loader.cpp
  cfg_resolved_case_v4_loader.cpp)
target_include_directories(hundun_config
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_config
  PRIVATE yyjson hundun_warnings
  PUBLIC hundun_options)

add_library(hundun_cli STATIC app_cli_options.cpp)
target_include_directories(hundun_cli
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_cli
  PRIVATE hundun_warnings
  PUBLIC hundun_options)

add_library(hundun_runtime STATIC
  rt_checkpoint_v2_protocol.cpp
  rt_mpi_error.cpp
  rt_mpi_environment.cpp
  rt_mpi_context.cpp
  rt_collective_status.cpp
  rt_process_grid.cpp
  rt_structured_decomposition.cpp
  rt_field_registry.cpp
  rt_field_access_plan.cpp
  rt_field_storage.cpp
  rt_halo_detail.cpp
  rt_exchange_plan.cpp
  rt_halo_exchange.cpp)
target_include_directories(hundun_runtime
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_runtime
  PRIVATE hundun_warnings
  PUBLIC hundun_options MPI::MPI_CXX)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_runtime PRIVATE
    HUNDUN_RUNTIME_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_application STATIC
  app_case_config_broadcast.cpp
  app_resolved_case_v3_broadcast.cpp
  app_resolved_case_v4_broadcast.cpp)
target_include_directories(hundun_application
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_application
  PRIVATE hundun_runtime hundun_warnings
  PUBLIC hundun_config hundun_options MPI::MPI_CXX)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_application PRIVATE
    HUNDUN_APPLICATION_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_mesh STATIC
  mesh_geometry.cpp
  mesh_topology.cpp
  mesh_uniform_structured.cpp)
target_include_directories(hundun_mesh
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_mesh
  PRIVATE hundun_warnings
  PUBLIC hundun_options hundun_runtime)

add_library(hundun_boundary STATIC
  bc_basic_boundary.cpp)
target_include_directories(hundun_boundary
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_boundary
  PRIVATE hundun_warnings
  PUBLIC hundun_config hundun_mesh hundun_runtime hundun_options)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_boundary PRIVATE
    HUNDUN_BOUNDARY_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_immersed STATIC
  ib_ghost_stencil_plan.cpp
  ib_surface.cpp
  ib_domain.cpp
  ib_local_flow_pattern.cpp
  ib_quadratic_reconstruction.cpp
  ib_surface_query.cpp
  ib_wall_force.cpp)
target_include_directories(hundun_immersed
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_immersed
  PRIVATE
    hundun_warnings
    MPI::MPI_CXX
    Threads::Threads
  PUBLIC
    hundun_boundary
    hundun_config
    hundun_mesh
    hundun_runtime
    hundun_options)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_immersed PRIVATE
    HUNDUN_IMMERSED_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_linear STATIC
  lin_bicgstab.cpp
  lin_conjugate_gradient.cpp
  lin_restarted_gmres.cpp
  lin_ghosted_vector.cpp
  lin_ghosted_vector_halo.cpp
  lin_ghosted_vector_halo_detail.cpp
  lin_preconditioners.cpp
  lin_solve_collective_detail.cpp
  lin_vector_ops.cpp)
target_include_directories(hundun_linear
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
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

add_library(hundun_fvm STATIC
  fvm_cell_centered.cpp
  fvm_immersed_operator.cpp
  fvm_immersed_reconstruction.cpp
  fvm_matrix_free_poisson.cpp
  fvm_poisson_boundary_adapter.cpp)
target_include_directories(hundun_fvm PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>")
target_link_libraries(hundun_fvm
  PRIVATE hundun_warnings
  PUBLIC
    hundun_boundary
    hundun_immersed
    hundun_linear
    hundun_mesh
    hundun_runtime
    hundun_execution
    hundun_options)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_fvm PUBLIC
    HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_flow STATIC
  flow_checkpoint_v2.cpp
  flow_checkpoint_v3.cpp
  flow_checkpoint_v4.cpp
  flow_adaptive_time_control.cpp
  flow_constant_density_piso.cpp
  flow_state.cpp
  flow_ideal_gas_closure.cpp
  flow_ideal_gas_piso.cpp
  flow_material_density_piso.cpp
  flow_material_density_transport.cpp
  flow_momentum_predictor.cpp
  flow_reacting_boundary.cpp
  flow_reacting_coupling.cpp
  flow_reacting_immersed.cpp
  flow_reacting_pressure_constraint.cpp
  flow_reacting_state.cpp
  flow_reacting_transport.cpp
  flow_reacting_transaction.cpp
  flow_immersed.cpp)
target_include_directories(hundun_flow
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_flow
  PRIVATE hundun_warnings
  PUBLIC
    hundun_fvm
    hundun_boundary
    hundun_linear
    hundun_mesh
    hundun_runtime
    hundun_chemistry
    hundun_execution
    hundun_options)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_flow PUBLIC
    HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_les STATIC
  les_wale.cpp)
target_include_directories(hundun_les
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_les
  PRIVATE hundun_warnings
  PUBLIC hundun_execution hundun_mesh hundun_runtime hundun_options)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_les PRIVATE
    HUNDUN_LES_ENABLE_TEST_ACCESS=1)
endif()
target_link_libraries(hundun_flow PUBLIC hundun_les)

add_library(hundun_material_diagnostics STATIC
  diag_ideal_gas_closure.cpp
  diag_material_density_piso.cpp
  diag_material_density_transport.cpp
  diag_time_control.cpp)
target_include_directories(hundun_material_diagnostics
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
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
  diag_checkpoint_v2.cpp
  diag_checkpoint_v3.cpp)
target_include_directories(hundun_checkpoint_diagnostics
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_checkpoint_diagnostics
  PRIVATE hundun_warnings
  PUBLIC
    hundun_diagnostics_core
    hundun_flow
    hundun_runtime
    hundun_options
    MPI::MPI_CXX)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_checkpoint_diagnostics PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
endif()

add_library(hundun_session_diagnostics STATIC
  diag_session.cpp
  diag_mesh_v2.cpp
  diag_module.cpp
  diag_immersed_module.cpp
  diag_immersed_static.cpp
  diag_les_wale.cpp)
target_include_directories(hundun_session_diagnostics
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_session_diagnostics
  PRIVATE hundun_warnings
  PUBLIC
    hundun_diagnostics_core
    hundun_boundary
    hundun_fvm
    hundun_flow
    hundun_checkpoint_diagnostics
    hundun_linear
    hundun_mesh
    hundun_runtime
    hundun_execution
    hundun_options
    MPI::MPI_CXX)

add_library(hundun_io STATIC
  rt_restart_binary.cpp
  rt_vtk_legacy.cpp)
target_include_directories(hundun_io
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_io
  PRIVATE hundun_warnings
  PUBLIC hundun_options hundun_runtime hundun_mesh)

add_library(hundun_solver STATIC flow_passive_scalar.cpp)
target_include_directories(hundun_solver
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_solver
  PRIVATE hundun_options hundun_warnings
  PUBLIC hundun_mesh hundun_runtime)

set(HUNDUN_SOURCE_COMMIT "unavailable")
set(HUNDUN_SOURCE_CLEAN "false")
set(HUNDUN_SOURCE_DIRTY_SUMMARY "unavailable")
include(HundunPerformanceSourceMetadata)
find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --verify HEAD
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE HUNDUN_GIT_COMMIT_RESULT
    OUTPUT_VARIABLE HUNDUN_GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if(HUNDUN_GIT_COMMIT_RESULT EQUAL 0 AND
     HUNDUN_GIT_COMMIT MATCHES "^[0-9a-fA-F]+$")
    string(TOLOWER "${HUNDUN_GIT_COMMIT}" HUNDUN_SOURCE_COMMIT)
  endif()
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" status --porcelain=v1 --untracked-files=normal
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE HUNDUN_GIT_STATUS_RESULT
    OUTPUT_VARIABLE HUNDUN_GIT_STATUS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  hundun_source_status_metadata(
    "${HUNDUN_GIT_STATUS_RESULT}" "${HUNDUN_GIT_STATUS}"
    HUNDUN_SOURCE_CLEAN HUNDUN_SOURCE_DIRTY_SUMMARY)
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" HUNDUN_BUILD_TYPE_UPPER)
set(HUNDUN_PERFORMANCE_COMPILER_FLAGS
  "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${HUNDUN_BUILD_TYPE_UPPER}}")
if(MPI_CXX_COMPILE_OPTIONS)
  string(JOIN " " HUNDUN_MPI_CXX_COMPILE_OPTIONS
    ${MPI_CXX_COMPILE_OPTIONS})
  string(APPEND HUNDUN_PERFORMANCE_COMPILER_FLAGS
    " ${HUNDUN_MPI_CXX_COMPILE_OPTIONS}")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  string(APPEND HUNDUN_PERFORMANCE_COMPILER_FLAGS
    " -Wall -Wextra -Wpedantic -Wconversion -Wshadow")
endif()
string(APPEND HUNDUN_PERFORMANCE_COMPILER_FLAGS " -std=c++17")
string(STRIP "${HUNDUN_PERFORMANCE_COMPILER_FLAGS}"
  HUNDUN_PERFORMANCE_COMPILER_FLAGS)
set(HUNDUN_PERFORMANCE_LINK_FLAGS
  "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${HUNDUN_BUILD_TYPE_UPPER}}")
if(MPI_CXX_LINK_FLAGS)
  string(APPEND HUNDUN_PERFORMANCE_LINK_FLAGS " ${MPI_CXX_LINK_FLAGS}")
endif()
string(STRIP "${HUNDUN_PERFORMANCE_LINK_FLAGS}"
  HUNDUN_PERFORMANCE_LINK_FLAGS)
foreach(HUNDUN_METADATA_FIELD IN ITEMS
    HUNDUN_SOURCE_COMMIT
    HUNDUN_SOURCE_DIRTY_SUMMARY
    HUNDUN_PERFORMANCE_COMPILER_FLAGS
    HUNDUN_PERFORMANCE_LINK_FLAGS)
  hundun_escape_cpp_string(
    "${${HUNDUN_METADATA_FIELD}}" HUNDUN_METADATA_ESCAPED)
  set("${HUNDUN_METADATA_FIELD}_ESCAPED" "${HUNDUN_METADATA_ESCAPED}")
endforeach()
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
configure_file(
  app_version.hpp.in
  generated/app_version.hpp
  @ONLY)
configure_file(
  app_performance_build.hpp.in
  generated/app_performance_build.hpp
  @ONLY)

add_executable(hundun
  app_main.cpp
  app_passive_scalar_driver.cpp
  app_flow_driver.cpp
  app_immersed_flow_driver.cpp
  app_reacting_flow_driver.cpp)
target_include_directories(hundun PRIVATE
  "${PROJECT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_BINARY_DIR}/generated")
target_link_libraries(hundun PRIVATE
  hundun_application
  hundun_cli
  hundun_config
  hundun_checkpoint_diagnostics
  hundun_diagnostics
  hundun_material_diagnostics
  hundun_session_diagnostics
  hundun_flow
  hundun_io
  hundun_solver
  hundun_options
  hundun_warnings
  MPI::MPI_CXX)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun PRIVATE
    HUNDUN_APPLICATION_ENABLE_TEST_ACCESS=1)
endif()
if(HUNDUN_ENABLE_CANTERA)
  hundun_apply_relocatable_target(hundun)
endif()

add_library(hundun_sdk STATIC sdk_plugin_loader.cpp)
target_include_directories(hundun_sdk
  PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(hundun_sdk
  PRIVATE hundun_warnings
  PUBLIC hundun_options ${CMAKE_DL_LIBS})

set(HUNDUN_PRODUCT_STATIC_LIBRARIES
  hundun_chemistry
  hundun_execution
  hundun_diagnostics
  hundun_diagnostics_core
  hundun_config
  hundun_cli
  hundun_runtime
  hundun_application
  hundun_mesh
  hundun_boundary
  hundun_immersed
  hundun_linear
  hundun_fvm
  hundun_flow
  hundun_les
  hundun_material_diagnostics
  hundun_checkpoint_diagnostics
  hundun_session_diagnostics
  hundun_io
  hundun_solver
  hundun_sdk)
if(HUNDUN_ENABLE_CANTERA)
  list(APPEND HUNDUN_PRODUCT_STATIC_LIBRARIES hundun_cantera_backend)
endif()

install(TARGETS hundun yyjson ${HUNDUN_PRODUCT_STATIC_LIBRARIES}
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/hundun/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/hundun"
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp")
