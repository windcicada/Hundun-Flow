# SPDX-License-Identifier: Apache-2.0

include_directories(
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_SOURCE_DIR}/include"
  "${PROJECT_SOURCE_DIR}/src")

include("${PROJECT_SOURCE_DIR}/tests/cmake/cmake_include_authority_fixture.cmake")
hundun_validate_configured_product_include_authority()

add_test(NAME source_layout_fixture
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_LAYOUT_ROOT=${PROJECT_SOURCE_DIR}"
    -P "${PROJECT_SOURCE_DIR}/tests/cmake/source_layout_fixture.cmake")
set_tests_properties(source_layout_fixture PROPERTIES
  LABELS "unit;governance;layout")

add_test(NAME source_layout_mutation_fixture
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
    "-DHUNDUN_TEST_BINARY_ROOT=${CMAKE_BINARY_DIR}"
    -P "${PROJECT_SOURCE_DIR}/tests/cmake/source_layout_mutation_fixture.cmake")
set_tests_properties(source_layout_mutation_fixture PROPERTIES
  LABELS "unit;governance;layout;mutation"
  RUN_SERIAL TRUE
  TIMEOUT 10)

add_test(NAME cmake_include_authority_fixture
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_LAYOUT_ROOT=${PROJECT_SOURCE_DIR}"
    -P "${PROJECT_SOURCE_DIR}/tests/cmake/cmake_include_authority_fixture.cmake")
set_tests_properties(cmake_include_authority_fixture PROPERTIES
  LABELS "unit;governance;layout;cmake_authority"
  TIMEOUT 10)

add_test(NAME cmake_include_authority_mutation_fixture
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
    "-DHUNDUN_TEST_BINARY_ROOT=${CMAKE_BINARY_DIR}"
    -P "${PROJECT_SOURCE_DIR}/tests/cmake/cmake_include_authority_mutation_fixture.cmake")
set_tests_properties(cmake_include_authority_mutation_fixture PROPERTIES
  LABELS "unit;governance;layout;cmake_authority;mutation"
  RUN_SERIAL TRUE
  TIMEOUT 10)

add_test(NAME product_projection_fixture
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
    "-DHUNDUN_TEST_BINARY_ROOT=${CMAKE_BINARY_DIR}"
    "-DHUNDUN_GENERATOR=${CMAKE_GENERATOR}"
    "-DHUNDUN_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}"
    "-DHUNDUN_C_COMPILER=${CMAKE_C_COMPILER}"
    "-DHUNDUN_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    "-DHUNDUN_C_COMPILER_ID=${CMAKE_C_COMPILER_ID}"
    "-DHUNDUN_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}"
    "-DHUNDUN_C_FLAGS=${CMAKE_C_FLAGS}"
    "-DHUNDUN_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
    "-DHUNDUN_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}"
    -P "${PROJECT_SOURCE_DIR}/tests/cmake/product_projection_fixture.cmake")
set_tests_properties(product_projection_fixture PROPERTIES
  LABELS "unit;governance;layout;product_projection"
  RUN_SERIAL TRUE
  TIMEOUT 120)

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

  add_test(NAME test_stage3_registration_contract
    COMMAND "${CMAKE_COMMAND}"
      "-DHUNDUN_STAGE3_REGISTRATION_ROOT=${PROJECT_SOURCE_DIR}"
      -P "${PROJECT_SOURCE_DIR}/tests/cmake/stage3_registration_contract.cmake")
  set_tests_properties(test_stage3_registration_contract PROPERTIES
    LABELS "unit;red;stage3;s3-p0;fast"
    TIMEOUT 60)

  add_test(NAME test_stage3_registration_contract_mutation
    COMMAND "${CMAKE_COMMAND}"
      "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
      "-DHUNDUN_TEST_BINARY_ROOT=${CMAKE_BINARY_DIR}"
      -P "${PROJECT_SOURCE_DIR}/tests/cmake/stage3_registration_contract_fixture.cmake")
  set_tests_properties(test_stage3_registration_contract_mutation PROPERTIES
    LABELS "unit;red;stage3;s3-p0;fast;mutation"
    RUN_SERIAL TRUE
    TIMEOUT 60)

  include("${PROJECT_SOURCE_DIR}/tests/cmake/stage3_science_registration.cmake")
  include("${PROJECT_SOURCE_DIR}/tests/cmake/stage3_checkpoint_registration.cmake")
  include("${PROJECT_SOURCE_DIR}/tests/cmake/stage3_diagnostics_registration.cmake")
  include("${PROJECT_SOURCE_DIR}/tests/cmake/stage3_framework_registration.cmake")
  include("${PROJECT_SOURCE_DIR}/tests/cmake/stage3_acceptance_registration.cmake")

  if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND
     NOT HUNDUN_ENABLE_ASAN AND NOT HUNDUN_ENABLE_UBSAN)
    add_test(NAME test_task26_stage2_acceptance_contract
      COMMAND "${CMAKE_COMMAND}"
        "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
        "-DHUNDUN_BINARY_ROOT=${CMAKE_BINARY_DIR}"
        "-DHUNDUN_STAGE2_ACCEPTANCE_SCRIPT=${PROJECT_SOURCE_DIR}/tests/acceptance/stage2_acceptance.sh"
        "-DHUNDUN_STAGE2_CAPABILITY_LEDGER=${PROJECT_SOURCE_DIR}/docs/numerics/stage2-capability-ledger.md"
        "-DHUNDUN_CTEST_COMMAND=${CMAKE_CTEST_COMMAND}"
        -P "${PROJECT_SOURCE_DIR}/tests/cmake/task26_stage2_acceptance_contract.cmake")
    set_tests_properties(test_task26_stage2_acceptance_contract PROPERTIES
      LABELS "unit;contract;acceptance;task26"
      RUN_SERIAL TRUE
      TIMEOUT 30)
  endif()

  hundun_add_test(test_task24_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_task24_header_contract.cpp)
  target_include_directories(test_task24_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_task24_header_contract PRIVATE
    hundun_application hundun_session_diagnostics)
  set_tests_properties(test_task24_header_contract PROPERTIES
    LABELS "unit;task24;header")

  hundun_add_test(test_mesh_diagnostic_v2
    ${PROJECT_SOURCE_DIR}/tests/unit/test_mesh_diagnostic_v2.cpp)
  target_include_directories(test_mesh_diagnostic_v2 PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_mesh_diagnostic_v2 PRIVATE
    hundun_session_diagnostics)
  set_tests_properties(test_mesh_diagnostic_v2 PROPERTIES
    LABELS "unit;task24;io")

  hundun_add_test(test_diagnostic_session
    ${PROJECT_SOURCE_DIR}/tests/unit/test_diagnostic_session.cpp)
  target_include_directories(test_diagnostic_session PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_diagnostic_session PRIVATE
    hundun_session_diagnostics)
  set_tests_properties(test_diagnostic_session PROPERTIES
    LABELS "unit;task24;diagnostics")

  hundun_add_test(test_stage2_module_diagnostics
    ${PROJECT_SOURCE_DIR}/tests/unit/test_stage2_module_diagnostics.cpp)
  target_include_directories(test_stage2_module_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_stage2_module_diagnostics PRIVATE
    hundun_session_diagnostics)
  set_tests_properties(test_stage2_module_diagnostics PROPERTIES
    LABELS "unit;task24;diagnostics")

  add_executable(test_diagnostic_session_mpi
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_diagnostic_session_mpi.cpp)
  target_include_directories(test_diagnostic_session_mpi PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_diagnostic_session_mpi PRIVATE
    hundun_session_diagnostics MPI::MPI_CXX)
  foreach(task24_session_ranks IN ITEMS 1 2 4)
    add_test(NAME
      "test_diagnostic_session_mpi_${task24_session_ranks}rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${task24_session_ranks}"
        "$<TARGET_FILE:test_diagnostic_session_mpi>")
    set_tests_properties(
      "test_diagnostic_session_mpi_${task24_session_ranks}rank"
      PROPERTIES LABELS "io;mpi;task24;diagnostics"
      TIMEOUT 60 PROCESSORS "${task24_session_ranks}")
  endforeach()

  foreach(task24_ranks IN ITEMS 1 2 4)
    add_test(NAME "test_task24_flow_models_${task24_ranks}rank"
      COMMAND bash
        "${PROJECT_SOURCE_DIR}/tests/acceptance/task24_flow_models.sh"
        "$<TARGET_FILE:hundun>" "${MPIEXEC_EXECUTABLE}" "${task24_ranks}")
    set_tests_properties("test_task24_flow_models_${task24_ranks}rank"
      PROPERTIES LABELS "application;mpi;task24;task25"
      TIMEOUT 120 PROCESSORS "${task24_ranks}")
    add_test(NAME "test_task24_restart_${task24_ranks}rank"
      COMMAND bash
        "${PROJECT_SOURCE_DIR}/tests/acceptance/task24_restart_continuation.sh"
        "$<TARGET_FILE:hundun>" "${MPIEXEC_EXECUTABLE}" "${task24_ranks}")
    set_tests_properties("test_task24_restart_${task24_ranks}rank"
      PROPERTIES LABELS "application;io;mpi;task24"
      TIMEOUT 180 PROCESSORS "${task24_ranks}")
  endforeach()

  function(hundun_add_mock_plugin target)
    add_library("${target}" MODULE ${PROJECT_SOURCE_DIR}/tests/plugins/mock_plugin.c)
    target_include_directories("${target}" PRIVATE
      "${PROJECT_SOURCE_DIR}/include")
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

  add_executable(test_plugin_loader ${PROJECT_SOURCE_DIR}/tests/unit/test_plugin_loader.cpp)
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
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_plugin_loader PRIVATE
    hundun_sdk hundun_options hundun_warnings Threads::Threads)
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

  hundun_add_test(test_case_config ${PROJECT_SOURCE_DIR}/tests/unit/test_case_config.cpp)
  target_include_directories(test_case_config PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_case_config PRIVATE hundun_config)

  hundun_add_test(test_resolved_case ${PROJECT_SOURCE_DIR}/tests/unit/test_resolved_case.cpp)
  target_include_directories(test_resolved_case PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_resolved_case PRIVATE hundun_config)

  hundun_add_test(test_resolved_case_v3
    ${PROJECT_SOURCE_DIR}/tests/unit/test_resolved_case_v3.cpp)
  target_include_directories(test_resolved_case_v3 PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_resolved_case_v3 PRIVATE hundun_config)

  hundun_add_test(test_resolved_case_v4
    ${PROJECT_SOURCE_DIR}/tests/unit/test_resolved_case_v4.cpp)
  target_include_directories(test_resolved_case_v4 PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_resolved_case_v4 PRIVATE hundun_config)
  set_tests_properties(test_resolved_case_v4 PROPERTIES
    LABELS "unit;stage4;task4f4" TIMEOUT 30)

  hundun_add_test(test_checkpoint_v4
    ${PROJECT_SOURCE_DIR}/tests/unit/test_checkpoint_v4.cpp)
  target_include_directories(test_checkpoint_v4 PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_checkpoint_v4 PRIVATE hundun_flow hundun_runtime)
  set_tests_properties(test_checkpoint_v4 PROPERTIES
    LABELS "unit;header;stage4;task4f5" TIMEOUT 30)

  hundun_add_mpi_test(test_checkpoint_v4_mpi_1_rank 1
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_checkpoint_v4_mpi.cpp)
  target_link_libraries(test_checkpoint_v4_mpi_1_rank PRIVATE
    hundun_flow hundun_runtime)
  set_tests_properties(test_checkpoint_v4_mpi_1_rank PROPERTIES
    LABELS "mpi;stage4;4A-1" TIMEOUT 30)

  hundun_add_mpi_test(test_checkpoint_v4_mpi_2_rank 2
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_checkpoint_v4_mpi.cpp)
  target_link_libraries(test_checkpoint_v4_mpi_2_rank PRIVATE
    hundun_flow hundun_runtime)
  set_tests_properties(test_checkpoint_v4_mpi_2_rank PROPERTIES
    LABELS "mpi;stage4;4A-1" TIMEOUT 30)

  hundun_add_test(test_reacting_coupling
    ${PROJECT_SOURCE_DIR}/tests/unit/test_reacting_coupling.cpp)
  target_include_directories(test_reacting_coupling PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_reacting_coupling PRIVATE
    hundun_flow hundun_runtime MPI::MPI_CXX)
  set_tests_properties(test_reacting_coupling PROPERTIES
    LABELS "unit;stage4;4R-0;4R-3" TIMEOUT 30)
  add_test(NAME test_reacting_coupling_2_rank
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_reacting_coupling>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_reacting_coupling_2_rank PROPERTIES
    LABELS "mpi;stage4;4R-3" TIMEOUT 30)

  hundun_add_test(test_reacting_mms
    ${PROJECT_SOURCE_DIR}/tests/numerical/test_reacting_mms.cpp)
  target_include_directories(test_reacting_mms PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_reacting_mms PRIVATE hundun_flow hundun_runtime)
  set_tests_properties(test_reacting_mms PROPERTIES
    LABELS "numerical;stage4;4R-0;4R-3" TIMEOUT 30)

  add_executable(test_reacting_state
    ${PROJECT_SOURCE_DIR}/tests/unit/test_reacting_state.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/chem_analytic_backend.cpp)
  target_include_directories(test_reacting_state PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_reacting_state PRIVATE
    hundun_flow hundun_warnings hundun_options)
  add_test(NAME test_reacting_state COMMAND test_reacting_state)
  set_tests_properties(test_reacting_state PROPERTIES
    LABELS "unit;header;stage4;4R-1" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_transport 1
    ${PROJECT_SOURCE_DIR}/tests/unit/test_reacting_transport.cpp)
  target_link_libraries(test_reacting_transport PRIVATE hundun_flow)
  set_tests_properties(test_reacting_transport PROPERTIES
    LABELS "unit;stage4;4R-2" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_decomposition_1_rank 1
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_reacting_decomposition.cpp)
  target_link_libraries(test_reacting_decomposition_1_rank PRIVATE hundun_flow)
  set_tests_properties(test_reacting_decomposition_1_rank PROPERTIES
    LABELS "mpi;stage4;4R-2" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_decomposition_2_rank 2
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_reacting_decomposition.cpp)
  target_link_libraries(test_reacting_decomposition_2_rank PRIVATE hundun_flow)
  set_tests_properties(test_reacting_decomposition_2_rank PROPERTIES
    LABELS "mpi;stage4;4R-2" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_boundary 1
    ${PROJECT_SOURCE_DIR}/tests/unit/test_reacting_boundary.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/chem_analytic_backend.cpp)
  target_link_libraries(test_reacting_boundary PRIVATE
    hundun_flow hundun_chemistry hundun_mesh hundun_runtime)
  set_tests_properties(test_reacting_boundary PROPERTIES
    LABELS "unit;stage4;4R-4" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_immersed_1_rank 1
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_reacting_immersed.cpp)
  target_link_libraries(test_reacting_immersed_1_rank PRIVATE
    hundun_flow hundun_les)
  set_tests_properties(test_reacting_immersed_1_rank PROPERTIES
    LABELS "mpi;stage4;4R-5" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_immersed_2_rank 2
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_reacting_immersed.cpp)
  target_link_libraries(test_reacting_immersed_2_rank PRIVATE
    hundun_flow hundun_les)
  set_tests_properties(test_reacting_immersed_2_rank PROPERTIES
    LABELS "mpi;stage4;4R-5" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_closed_pressure_1_rank 1
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_reacting_closed_pressure.cpp)
  target_link_libraries(test_reacting_closed_pressure_1_rank PRIVATE
    hundun_flow hundun_runtime)
  set_tests_properties(test_reacting_closed_pressure_1_rank PROPERTIES
    LABELS "mpi;stage4;4R-6" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_closed_pressure_2_rank 2
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_reacting_closed_pressure.cpp)
  target_link_libraries(test_reacting_closed_pressure_2_rank PRIVATE
    hundun_flow hundun_runtime)
  set_tests_properties(test_reacting_closed_pressure_2_rank PROPERTIES
    LABELS "mpi;stage4;4R-6" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_closed_pressure_4_rank 4
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_reacting_closed_pressure.cpp)
  target_link_libraries(test_reacting_closed_pressure_4_rank PRIVATE
    hundun_flow hundun_runtime)
  set_tests_properties(test_reacting_closed_pressure_4_rank PROPERTIES
    LABELS "mpi;stage4;4R-6" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_driver 1
    ${PROJECT_SOURCE_DIR}/tests/unit/test_reacting_driver.cpp
    ${PROJECT_SOURCE_DIR}/src/app_reacting_flow_driver.cpp)
  target_link_libraries(test_reacting_driver PRIVATE
    hundun_config hundun_flow hundun_runtime)
  set_tests_properties(test_reacting_driver PROPERTIES
    LABELS "unit;stage4;4R-7" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_diagnostics_1_rank 1
    ${PROJECT_SOURCE_DIR}/tests/unit/test_reacting_diagnostics.cpp)
  target_link_libraries(test_reacting_diagnostics_1_rank PRIVATE
    hundun_diagnostics_core hundun_runtime)
  set_tests_properties(test_reacting_diagnostics_1_rank PROPERTIES
    LABELS "unit;stage4;4A-2" TIMEOUT 30)

  hundun_add_mpi_test(test_reacting_diagnostics_2_rank 2
    ${PROJECT_SOURCE_DIR}/tests/unit/test_reacting_diagnostics.cpp)
  target_link_libraries(test_reacting_diagnostics_2_rank PRIVATE
    hundun_diagnostics_core hundun_runtime)
  set_tests_properties(test_reacting_diagnostics_2_rank PROPERTIES
    LABELS "mpi;stage4;4A-2" TIMEOUT 30)

  add_test(NAME stage4_acceptance_selector_contract
    COMMAND bash "${PROJECT_SOURCE_DIR}/tests/acceptance/stage4_acceptance.sh"
            --list)
  set_tests_properties(stage4_acceptance_selector_contract PROPERTIES
    LABELS "unit;stage4;4A-3" TIMEOUT 30)

  hundun_add_test(test_diag_provider_registry
    ${PROJECT_SOURCE_DIR}/tests/unit/test_diag_provider_registry.cpp)
  target_include_directories(test_diag_provider_registry PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_diag_provider_registry PRIVATE
    hundun_diagnostics_core)
  set_tests_properties(test_diag_provider_registry PROPERTIES
    LABELS "unit;stage4;task4f5" TIMEOUT 30)

  hundun_add_test(test_resolved_case_v3_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_resolved_case_v3_header_contract.cpp)
  target_link_libraries(test_resolved_case_v3_header_contract PRIVATE
    hundun_application)

  hundun_add_test(test_immersed_surface_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_immersed_surface_header_contract.cpp)
  target_link_libraries(test_immersed_surface_header_contract PRIVATE
    hundun_immersed)
  set_tests_properties(test_immersed_surface_header_contract PROPERTIES
    LABELS "unit;header;stage3;task2")

  hundun_add_test(test_immersed_domain_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_immersed_domain_header_contract.cpp)
  target_link_libraries(test_immersed_domain_header_contract PRIVATE
    hundun_immersed)
  set_tests_properties(test_immersed_domain_header_contract PROPERTIES
    LABELS "unit;header;stage3;task3"
    TIMEOUT 30)

  hundun_add_test(test_quadratic_reconstruction_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_quadratic_reconstruction_header_contract.cpp)
  target_link_libraries(test_quadratic_reconstruction_header_contract PRIVATE
    hundun_immersed)
  set_tests_properties(
    test_quadratic_reconstruction_header_contract PROPERTIES
    LABELS "unit;header;stage3;task4"
    TIMEOUT 30)

  hundun_add_test(test_ghost_stencil_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_ghost_stencil_header_contract.cpp)
  target_link_libraries(test_ghost_stencil_header_contract PRIVATE
    hundun_immersed)
  set_tests_properties(test_ghost_stencil_header_contract PROPERTIES
    LABELS "unit;header;stage3;task5"
    TIMEOUT 30)

  hundun_add_test(test_local_flow_pattern_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_local_flow_pattern_header_contract.cpp)
  target_link_libraries(test_local_flow_pattern_header_contract PRIVATE
    hundun_immersed)
  set_tests_properties(test_local_flow_pattern_header_contract PROPERTIES
    LABELS "unit;header;stage3;task6"
    TIMEOUT 30)

  hundun_add_test(test_local_flow_pattern
    ${PROJECT_SOURCE_DIR}/tests/unit/test_local_flow_pattern.cpp)
  target_include_directories(test_local_flow_pattern PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_local_flow_pattern PRIVATE hundun_immersed)
  target_compile_definitions(test_local_flow_pattern PRIVATE
    HUNDUN_IMMERSED_ENABLE_TEST_ACCESS=1)
  set_tests_properties(test_local_flow_pattern PROPERTIES
    LABELS "unit;stage3;task6"
    TIMEOUT 60)

  add_library(stage3_stl_fixture STATIC
    ${PROJECT_SOURCE_DIR}/tests/support/stage3_stl_fixture.cpp)
  target_include_directories(stage3_stl_fixture PUBLIC
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(stage3_stl_fixture
    PRIVATE hundun_warnings
    PUBLIC hundun_options hundun_runtime)

  add_library(stage3_mms STATIC
    ${PROJECT_SOURCE_DIR}/tests/support/stage3_mms.cpp)
  target_include_directories(stage3_mms PUBLIC
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(stage3_mms
    PRIVATE hundun_warnings
    PUBLIC hundun_options hundun_flow hundun_mesh hundun_runtime)

  hundun_add_test(test_stage3_decomposition_equality
    ${PROJECT_SOURCE_DIR}/tests/unit/test_stage3_decomposition_equality.cpp)
  target_include_directories(test_stage3_decomposition_equality PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_stage3_decomposition_equality PRIVATE
    hundun_immersed)
  set_tests_properties(test_stage3_decomposition_equality PROPERTIES
    LABELS "unit;stage3;task11"
    TIMEOUT 30)

  hundun_add_test(test_immersed_surface
    ${PROJECT_SOURCE_DIR}/tests/unit/test_immersed_surface.cpp)
  target_include_directories(test_immersed_surface PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_immersed_surface PRIVATE
    hundun_immersed
    stage3_stl_fixture)
  target_compile_definitions(test_immersed_surface PRIVATE
    HUNDUN_IMMERSED_ENABLE_TEST_ACCESS=1)
  set_tests_properties(test_immersed_surface PROPERTIES
    LABELS "unit;stage3;task2"
    TIMEOUT 60)

  hundun_add_test(test_surface_query
    ${PROJECT_SOURCE_DIR}/tests/unit/test_surface_query.cpp)
  target_include_directories(test_surface_query PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_surface_query PRIVATE
    hundun_immersed
    stage3_stl_fixture)
  target_compile_definitions(test_surface_query PRIVATE
    HUNDUN_IMMERSED_ENABLE_TEST_ACCESS=1)
  set_tests_properties(test_surface_query PROPERTIES
    LABELS "unit;stage3;task2"
    TIMEOUT 60)

  add_executable(test_immersed_surface_mpi
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_immersed_surface_mpi.cpp)
  target_include_directories(test_immersed_surface_mpi PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_immersed_surface_mpi PRIVATE
    hundun_immersed
    stage3_stl_fixture
    hundun_options
    hundun_warnings
    MPI::MPI_CXX)
  target_compile_definitions(test_immersed_surface_mpi PRIVATE
    HUNDUN_IMMERSED_ENABLE_TEST_ACCESS=1)
  foreach(surface_ranks IN ITEMS 1 2 4)
    add_test(NAME "test_immersed_surface_mpi_${surface_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${surface_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_immersed_surface_mpi>"
              success ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_immersed_surface_mpi_${surface_ranks}_rank" PROPERTIES
      LABELS "mpi;stage3;task2"
      TIMEOUT 60
      PROCESSORS "${surface_ranks}")
  endforeach()
  add_test(NAME test_immersed_surface_mpi_failures
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_surface_mpi>"
            failures ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_immersed_surface_mpi_failures PROPERTIES
    LABELS "mpi;stage3;task2"
    TIMEOUT 60
    PROCESSORS 4)

  add_executable(test_immersed_domain
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_immersed_domain.cpp)
  target_include_directories(test_immersed_domain PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_immersed_domain PRIVATE
    hundun_immersed
    stage3_stl_fixture
    hundun_options
    hundun_warnings
    MPI::MPI_CXX)
  foreach(domain_ranks IN ITEMS 1 2)
    add_test(NAME "test_immersed_domain_${domain_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${domain_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_immersed_domain>"
              success default ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_immersed_domain_${domain_ranks}_rank" PROPERTIES
      LABELS "mpi;stage3;task3"
      TIMEOUT 180
      PROCESSORS "${domain_ranks}")
  endforeach()
  add_test(NAME test_immersed_domain_process_grid_2x2x1
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_domain>"
            success default ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_immersed_domain_process_grid_2x2x1 PROPERTIES
    LABELS "mpi;stage3;task3"
    TIMEOUT 180
    PROCESSORS 4)
  add_test(NAME test_immersed_domain_process_grid_4x1x1
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_domain>"
            success alternate ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_immersed_domain_process_grid_4x1x1 PROPERTIES
    LABELS "mpi;stage3;task3"
    TIMEOUT 180
    PROCESSORS 4)
  foreach(domain_failure_ranks IN ITEMS 2 4)
    add_test(NAME "test_immersed_domain_failures_${domain_failure_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${domain_failure_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_immersed_domain>"
              failures ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_immersed_domain_failures_${domain_failure_ranks}_rank" PROPERTIES
      LABELS "mpi;stage3;task3"
      TIMEOUT 120
      PROCESSORS "${domain_failure_ranks}")
  endforeach()

  add_executable(test_active_boundary_layout
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_active_boundary_layout.cpp)
  target_include_directories(test_active_boundary_layout PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_active_boundary_layout PRIVATE
    hundun_immersed
    stage3_stl_fixture
    hundun_options
    hundun_warnings
    MPI::MPI_CXX)
  foreach(active_boundary_ranks IN ITEMS 1 2)
    add_test(NAME "test_active_boundary_layout_${active_boundary_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${active_boundary_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_active_boundary_layout>"
              success default ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_active_boundary_layout_${active_boundary_ranks}_rank" PROPERTIES
      LABELS "mpi;stage3;task3"
      TIMEOUT 180
      PROCESSORS "${active_boundary_ranks}")
  endforeach()
  add_test(NAME test_active_boundary_layout_process_grid_2x2x1
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:test_active_boundary_layout>"
            success default ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    test_active_boundary_layout_process_grid_2x2x1 PROPERTIES
    LABELS "mpi;stage3;task3"
    TIMEOUT 180
    PROCESSORS 4)
  add_test(NAME test_active_boundary_layout_process_grid_4x1x1
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:test_active_boundary_layout>"
            success alternate ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    test_active_boundary_layout_process_grid_4x1x1 PROPERTIES
    LABELS "mpi;stage3;task3"
    TIMEOUT 180
    PROCESSORS 4)
  foreach(active_boundary_failure_ranks IN ITEMS 2 4)
    add_test(
      NAME
        "test_active_boundary_layout_failures_${active_boundary_failure_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${active_boundary_failure_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_active_boundary_layout>"
              failures ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_active_boundary_layout_failures_${active_boundary_failure_ranks}_rank"
      PROPERTIES
      LABELS "mpi;stage3;task3"
      TIMEOUT 120
      PROCESSORS "${active_boundary_failure_ranks}")
  endforeach()

  hundun_add_test(test_quadratic_reconstruction
    ${PROJECT_SOURCE_DIR}/tests/unit/test_quadratic_reconstruction.cpp)
  target_include_directories(test_quadratic_reconstruction PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_quadratic_reconstruction PRIVATE
    hundun_immersed)
  set_tests_properties(test_quadratic_reconstruction PROPERTIES
    LABELS "unit;stage3;task4"
    TIMEOUT 120)

  hundun_add_test(test_periodic_surface_window
    ${PROJECT_SOURCE_DIR}/benchmarks/cylinder_re3900/tests/test_periodic_surface_window.cpp)
  target_include_directories(test_periodic_surface_window PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src")
  target_compile_definitions(test_periodic_surface_window PRIVATE
    HUNDUN_CYLINDER_CASE_STL="${PROJECT_SOURCE_DIR}/benchmarks/cylinder_re3900/cases/hundun/geometry/cylinder_re3900.stl")
  target_link_libraries(test_periodic_surface_window PRIVATE
    hundun_immersed MPI::MPI_CXX)
  set_tests_properties(test_periodic_surface_window PROPERTIES
    LABELS "unit;stage3;periodic_ibm;red;fast"
    TIMEOUT 30)

  hundun_add_mpi_test(test_periodic_reconstruction_donors 1
    ${PROJECT_SOURCE_DIR}/benchmarks/cylinder_re3900/tests/test_periodic_reconstruction_donors.cpp)
  target_include_directories(test_periodic_reconstruction_donors PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_periodic_reconstruction_donors PRIVATE
    hundun_immersed hundun_mesh hundun_runtime)
  set_tests_properties(test_periodic_reconstruction_donors PROPERTIES
    LABELS "mpi;stage3;periodic_ibm;red;fast"
    TIMEOUT 30
    PROCESSORS 1)

  add_executable(test_quadratic_reconstruction_mpi
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_quadratic_reconstruction_mpi.cpp)
  target_include_directories(test_quadratic_reconstruction_mpi PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_quadratic_reconstruction_mpi PRIVATE
    hundun_immersed
    hundun_options
    hundun_warnings
    MPI::MPI_CXX)
  foreach(reconstruction_ranks IN ITEMS 1 2 4)
    add_test(NAME "test_quadratic_reconstruction_mpi_${reconstruction_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${reconstruction_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_quadratic_reconstruction_mpi>"
              success default ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_quadratic_reconstruction_mpi_${reconstruction_ranks}_rank"
      PROPERTIES
      LABELS "mpi;stage3;task4"
      TIMEOUT 120
      PROCESSORS "${reconstruction_ranks}")
  endforeach()
  foreach(reconstruction_grid IN ITEMS 2x2x1 4x1x1)
    add_test(
      NAME "test_quadratic_reconstruction_mpi_process_grid_${reconstruction_grid}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_quadratic_reconstruction_mpi>"
              success "${reconstruction_grid}" ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_quadratic_reconstruction_mpi_process_grid_${reconstruction_grid}"
      PROPERTIES
      LABELS "mpi;stage3;task4"
      TIMEOUT 120
      PROCESSORS 4)
  endforeach()
  foreach(reconstruction_failure_ranks IN ITEMS 1 2 4)
    add_test(
      NAME
        "test_quadratic_reconstruction_mpi_failures_${reconstruction_failure_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${reconstruction_failure_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_quadratic_reconstruction_mpi>"
              failures default ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_quadratic_reconstruction_mpi_failures_${reconstruction_failure_ranks}_rank"
      PROPERTIES
      LABELS "mpi;stage3;task4"
      TIMEOUT 120
      PROCESSORS "${reconstruction_failure_ranks}")
  endforeach()

  add_executable(test_ghost_stencil_plan
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_ghost_stencil_plan.cpp)
  target_include_directories(test_ghost_stencil_plan PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_ghost_stencil_plan PRIVATE
    hundun_immersed
    stage3_stl_fixture
    hundun_options
    hundun_warnings
    MPI::MPI_CXX)

  add_executable(test_wall_quadrature_plan
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_wall_quadrature_plan.cpp)
  target_include_directories(test_wall_quadrature_plan PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_wall_quadrature_plan PRIVATE
    hundun_immersed
    stage3_stl_fixture
    hundun_options
    hundun_warnings
    MPI::MPI_CXX)

  foreach(task5_ranks IN ITEMS 1 2 4)
    foreach(task5_target IN ITEMS ghost_stencil_plan wall_quadrature_plan)
      add_test(NAME "test_${task5_target}_${task5_ranks}_rank"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
                "${task5_ranks}" ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_${task5_target}>"
                success default ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        "test_${task5_target}_${task5_ranks}_rank" PROPERTIES
        LABELS "mpi;stage3;task5"
        TIMEOUT 180
        PROCESSORS "${task5_ranks}")
    endforeach()
  endforeach()
  foreach(task5_grid IN ITEMS 2x2x1 4x1x1)
    foreach(task5_target IN ITEMS ghost_stencil_plan wall_quadrature_plan)
      add_test(NAME "test_${task5_target}_process_grid_${task5_grid}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_${task5_target}>"
                success "${task5_grid}" ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        "test_${task5_target}_process_grid_${task5_grid}" PROPERTIES
        LABELS "mpi;stage3;task5"
        TIMEOUT 180
        PROCESSORS 4)
    endforeach()
  endforeach()
  foreach(task5_failure_ranks IN ITEMS 1 2 4)
    add_test(NAME
      "test_ghost_stencil_plan_failures_${task5_failure_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${task5_failure_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_ghost_stencil_plan>"
              failures default ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_ghost_stencil_plan_failures_${task5_failure_ranks}_rank"
      PROPERTIES
      LABELS "mpi;stage3;task5"
      TIMEOUT 180
      PROCESSORS "${task5_failure_ranks}")
  endforeach()
  foreach(task5_periodic_ranks IN ITEMS 1 2 4)
    add_test(NAME
      "test_ghost_stencil_plan_periodic_boundary_${task5_periodic_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${task5_periodic_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_ghost_stencil_plan>"
              periodic_boundary default ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_ghost_stencil_plan_periodic_boundary_${task5_periodic_ranks}_rank"
      PROPERTIES
      LABELS "mpi;stage3;task11"
      TIMEOUT 180
      PROCESSORS "${task5_periodic_ranks}")
  endforeach()
  foreach(task5_warped_ranks IN ITEMS 1 2 4)
    foreach(task5_target IN ITEMS ghost_stencil_plan wall_quadrature_plan)
      add_test(NAME "test_${task5_target}_warped_${task5_warped_ranks}_rank"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
                "${task5_warped_ranks}" ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_${task5_target}>"
                warped default ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        "test_${task5_target}_warped_${task5_warped_ranks}_rank"
        PROPERTIES
        LABELS "mpi;stage3;task5"
        TIMEOUT 180
        PROCESSORS "${task5_warped_ranks}")
    endforeach()
  endforeach()
  foreach(task5_disagreement_ranks IN ITEMS 2 4)
    add_test(NAME
      "test_ghost_stencil_plan_rank_disagreement_${task5_disagreement_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${task5_disagreement_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_ghost_stencil_plan>"
              rank_disagreement default ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_ghost_stencil_plan_rank_disagreement_${task5_disagreement_ranks}_rank"
      PROPERTIES
      LABELS "mpi;stage3;task5"
      TIMEOUT 180
      PROCESSORS "${task5_disagreement_ranks}")
  endforeach()

  add_executable(test_local_flow_pattern_mpi
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_local_flow_pattern_mpi.cpp)
  target_include_directories(test_local_flow_pattern_mpi PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_local_flow_pattern_mpi PRIVATE
    hundun_immersed
    hundun_options
    hundun_warnings
    MPI::MPI_CXX)
  foreach(task6_ranks IN ITEMS 1 2 4)
    add_test(NAME "test_local_flow_pattern_mpi_${task6_ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
              "${task6_ranks}" ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_local_flow_pattern_mpi>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_local_flow_pattern_mpi_${task6_ranks}_rank" PROPERTIES
      LABELS "mpi;stage3;task6"
      TIMEOUT 60
      PROCESSORS "${task6_ranks}")
  endforeach()

  add_test(NAME test_stage3_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/stage3_source_policy.cmake")
  set_tests_properties(test_stage3_source_policy PROPERTIES
    LABELS "unit;source_policy;stage3;task1" TIMEOUT 30)
  foreach(stage3_policy_fixture IN ITEMS
      script_dependency vendor_header stage4_key private_path
      missing_contract)
    set(stage3_policy_test
      "test_stage3_source_policy_rejects_${stage3_policy_fixture}")
    add_test(NAME "${stage3_policy_test}"
      COMMAND "${CMAKE_COMMAND}"
              "-DHUNDUN_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
              "-DHUNDUN_STAGE3_POLICY_VERIFY_FIXTURE=${stage3_policy_fixture}"
              -P
              "${PROJECT_SOURCE_DIR}/tests/cmake/stage3_source_policy.cmake")
    set_tests_properties("${stage3_policy_test}" PROPERTIES
      LABELS "unit;source_policy;stage3;task1"
      TIMEOUT 30)
  endforeach()

  add_test(NAME stage4_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/stage4_source_policy.cmake")
  set_tests_properties(stage4_source_policy PROPERTIES
    LABELS "unit;source_policy;stage4;4F-0" TIMEOUT 30)
  add_test(NAME stage4_cantera_builder_contract
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
            -P
            "${PROJECT_SOURCE_DIR}/tests/cmake/stage4_cantera_builder_contract.cmake")
  set_tests_properties(stage4_cantera_builder_contract PROPERTIES
    LABELS "unit;cmake;stage4;4P-2" TIMEOUT 60)
  add_test(NAME stage4_package_contract
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/stage4_package_contract.cmake")
  set_tests_properties(stage4_package_contract PROPERTIES
    LABELS "unit;cmake;stage4;4P-4" TIMEOUT 30)

  if(HUNDUN_ENABLE_CANTERA)
    if(NOT EXISTS "${HUNDUN_STAGE4_CANTERA_TEST_MECHANISM}")
      message(FATAL_ERROR
        "Stage 4 Cantera tests require HUNDUN_STAGE4_CANTERA_TEST_MECHANISM")
    endif()
    file(SHA256 "${HUNDUN_STAGE4_CANTERA_TEST_MECHANISM}"
      stage4_test_mechanism_sha256)
    if(NOT stage4_test_mechanism_sha256 STREQUAL
       "c518a07cada5f1bddcdb308f0a2f695d92cc6373e173ffd87e96312530b52aee")
      message(FATAL_ERROR "Stage 4 synthetic mechanism SHA-256 mismatch")
    endif()

    add_executable(cantera_cpp_link_spike
      ${PROJECT_SOURCE_DIR}/tests/integration/cantera_cpp_link_spike.cpp)
    target_link_libraries(cantera_cpp_link_spike PRIVATE
      hundun_third_party_cantera hundun_warnings hundun_options)
    add_test(NAME cantera_cpp_link_spike
      COMMAND cantera_cpp_link_spike
              "${HUNDUN_STAGE4_CANTERA_TEST_MECHANISM}")
    set_tests_properties(cantera_cpp_link_spike PROPERTIES
      LABELS "integration;stage4;4P-3" TIMEOUT 30)

    add_executable(test_cantera_backend
      ${PROJECT_SOURCE_DIR}/tests/unit/test_cantera_backend.cpp)
    target_include_directories(test_cantera_backend PRIVATE
      "${PROJECT_SOURCE_DIR}")
    target_link_libraries(test_cantera_backend PRIVATE
      hundun_cantera_backend hundun_warnings hundun_options)
    add_test(NAME test_cantera_backend_runtime
      COMMAND test_cantera_backend
              "${HUNDUN_STAGE4_CANTERA_TEST_MECHANISM}")
    set_tests_properties(test_cantera_backend_runtime PROPERTIES
      LABELS "unit;stage4;4C-1" TIMEOUT 30)

    add_executable(test_cantera_thermodynamics
      ${PROJECT_SOURCE_DIR}/tests/unit/test_cantera_thermodynamics.cpp)
    target_include_directories(test_cantera_thermodynamics PRIVATE
      "${PROJECT_SOURCE_DIR}")
    target_link_libraries(test_cantera_thermodynamics PRIVATE
      hundun_cantera_backend hundun_third_party_cantera
      hundun_warnings hundun_options)
    add_test(NAME test_cantera_thermodynamics
      COMMAND test_cantera_thermodynamics
              "${HUNDUN_STAGE4_CANTERA_TEST_MECHANISM}")
    set_tests_properties(test_cantera_thermodynamics PROPERTIES
      LABELS "unit;stage4;4C-2" TIMEOUT 30)

    add_executable(test_cantera_transport
      ${PROJECT_SOURCE_DIR}/tests/unit/test_cantera_transport.cpp)
    target_include_directories(test_cantera_transport PRIVATE
      "${PROJECT_SOURCE_DIR}")
    target_link_libraries(test_cantera_transport PRIVATE
      hundun_cantera_backend hundun_third_party_cantera
      hundun_warnings hundun_options)
    add_test(NAME test_cantera_transport
      COMMAND test_cantera_transport
              "${HUNDUN_STAGE4_CANTERA_TEST_MECHANISM}")
    set_tests_properties(test_cantera_transport PROPERTIES
      LABELS "unit;stage4;4C-3" TIMEOUT 30)

    add_executable(test_cantera_chemistry_interval
      ${PROJECT_SOURCE_DIR}/tests/unit/test_cantera_chemistry_interval.cpp)
    target_include_directories(test_cantera_chemistry_interval PRIVATE
      "${PROJECT_SOURCE_DIR}")
    target_link_libraries(test_cantera_chemistry_interval PRIVATE
      hundun_cantera_backend hundun_third_party_cantera
      hundun_warnings hundun_options)
    add_test(NAME test_cantera_chemistry_interval
      COMMAND test_cantera_chemistry_interval
              "${HUNDUN_STAGE4_CANTERA_TEST_MECHANISM}")
    set_tests_properties(test_cantera_chemistry_interval PROPERTIES
      LABELS "unit;stage4;4C-4" TIMEOUT 30)

    foreach(stage4_backend_conformance IN ITEMS
        reacting_zero_dimensional reacting_psr)
      add_executable(test_${stage4_backend_conformance}
        "${PROJECT_SOURCE_DIR}/tests/numerical/test_${stage4_backend_conformance}.cpp")
      target_include_directories(test_${stage4_backend_conformance} PRIVATE
        "${PROJECT_SOURCE_DIR}")
      target_link_libraries(test_${stage4_backend_conformance} PRIVATE
        hundun_cantera_backend hundun_third_party_cantera
        hundun_warnings hundun_options)
      add_test(NAME test_${stage4_backend_conformance}
        COMMAND test_${stage4_backend_conformance}
                "${HUNDUN_STAGE4_CANTERA_TEST_MECHANISM}")
      set_tests_properties(test_${stage4_backend_conformance} PROPERTIES
        LABELS "numerical;stage4;4C-5" TIMEOUT 45)
    endforeach()

    add_executable(test_cantera_workspace_isolation
      ${PROJECT_SOURCE_DIR}/tests/mpi/test_cantera_workspace_isolation.cpp)
    target_link_libraries(test_cantera_workspace_isolation PRIVATE
      hundun_third_party_cantera hundun_warnings hundun_options
      MPI::MPI_CXX Threads::Threads)
    foreach(ranks IN ITEMS 1 2)
      add_test(NAME "cantera_workspace_isolation_${ranks}_rank"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_cantera_workspace_isolation>"
                "${HUNDUN_STAGE4_CANTERA_TEST_MECHANISM}"
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties("cantera_workspace_isolation_${ranks}_rank"
        PROPERTIES LABELS "mpi;stage4;4P-3" TIMEOUT 45
                   PROCESSORS "${ranks}")
    endforeach()
  endif()
  foreach(stage4_policy_fixture IN ITEMS
      missing_allowlist
      public_cantera_type
      public_sundials_type
      public_python_header
      public_pybind_type
      second_plugin_abi
      private_path
      unregistered_prefix
      cantera_archive_hash
      cantera_missing_license
      cantera_product_source)
    set(stage4_policy_test
      "stage4_source_policy_rejects_${stage4_policy_fixture}")
    add_test(NAME "${stage4_policy_test}"
      COMMAND "${CMAKE_COMMAND}"
              "-DHUNDUN_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
              "-DHUNDUN_STAGE4_POLICY_VERIFY_FIXTURE=${stage4_policy_fixture}"
              -P
              "${PROJECT_SOURCE_DIR}/tests/cmake/stage4_source_policy.cmake")
    set_tests_properties("${stage4_policy_test}" PROPERTIES
      LABELS "unit;source_policy;stage4;4F-0"
      TIMEOUT 30)
  endforeach()

  hundun_add_test(test_chem_composition
    ${PROJECT_SOURCE_DIR}/tests/unit/test_chem_composition.cpp)
  target_include_directories(test_chem_composition PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_chem_composition PRIVATE hundun_chemistry)
  set_tests_properties(test_chem_composition PROPERTIES
    LABELS "unit;header;stage4;4F-1")

  hundun_add_test(test_chem_composition_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_chem_composition_header_contract.cpp)
  target_include_directories(test_chem_composition_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_chem_composition_header_contract PRIVATE
    hundun_chemistry)
  set_tests_properties(test_chem_composition_header_contract PROPERTIES
    LABELS "unit;header;stage4;4F-1")

  hundun_add_test(test_chem_services
    ${PROJECT_SOURCE_DIR}/tests/unit/test_chem_services.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/chem_analytic_backend.cpp)
  target_include_directories(test_chem_services PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_chem_services PRIVATE hundun_chemistry)
  set_tests_properties(test_chem_services PROPERTIES
    LABELS "unit;stage4;4F-2")

  hundun_add_test(test_reacting_transaction
    ${PROJECT_SOURCE_DIR}/tests/unit/test_reacting_transaction.cpp)
  target_include_directories(test_reacting_transaction PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_reacting_transaction PRIVATE
    hundun_flow hundun_chemistry hundun_runtime)
  set_tests_properties(test_reacting_transaction PROPERTIES
    LABELS "unit;stage4;4F-3")

  hundun_add_mpi_test(test_reacting_rollback 2
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_reacting_rollback.cpp)
  target_include_directories(test_reacting_rollback PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_reacting_rollback PRIVATE
    hundun_flow hundun_chemistry hundun_runtime)
  set_tests_properties(test_reacting_rollback PROPERTIES
    LABELS "mpi;stage4;4F-3" TIMEOUT 60 PROCESSORS 2)

  hundun_add_test(test_stage1_dispatch_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_stage1_dispatch_contract.cpp)
  target_include_directories(test_stage1_dispatch_contract PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_stage1_dispatch_contract PRIVATE
    hundun_application)

  hundun_add_test(test_process_grid
    ${PROJECT_SOURCE_DIR}/tests/unit/test_process_grid.cpp
    ${PROJECT_SOURCE_DIR}/src/rt_process_grid.cpp)
  target_include_directories(test_process_grid PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")

  hundun_add_test(test_cli_options ${PROJECT_SOURCE_DIR}/tests/unit/test_cli_options.cpp)
  target_include_directories(test_cli_options PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_cli_options PRIVATE hundun_cli)

  hundun_add_test(test_field_storage ${PROJECT_SOURCE_DIR}/tests/unit/test_field_storage.cpp)
  target_include_directories(test_field_storage PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_field_storage PRIVATE hundun_runtime)

  add_executable(test_checkpoint_v2_protocol
    ${PROJECT_SOURCE_DIR}/tests/unit/test_checkpoint_v2_protocol.cpp)
  target_include_directories(test_checkpoint_v2_protocol PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_checkpoint_v2_protocol PRIVATE
    hundun_flow hundun_runtime)
  target_compile_definitions(test_checkpoint_v2_protocol PRIVATE
    HUNDUN_RUNTIME_ENABLE_TEST_ACCESS=1)
  add_test(NAME test_checkpoint_v2_protocol
    COMMAND "$<TARGET_FILE:test_checkpoint_v2_protocol>" fast)
  set_tests_properties(test_checkpoint_v2_protocol PROPERTIES
    LABELS "unit;task23;task23-fast" TIMEOUT 30)
  add_test(NAME test_checkpoint_v2_protocol_acceptance
    COMMAND "$<TARGET_FILE:test_checkpoint_v2_protocol>" acceptance)
  set_tests_properties(test_checkpoint_v2_protocol_acceptance PROPERTIES
    LABELS "unit;task23;task23-acceptance;task23-cluster-a" TIMEOUT 60)
  foreach(selector_case IN ITEMS missing unknown extra)
    if(selector_case STREQUAL "missing")
      set(selector_args)
    elseif(selector_case STREQUAL "unknown")
      set(selector_args "not-a-task23-mode")
    else()
      set(selector_args "acceptance" "unexpected-extra-argument")
    endif()
    add_test(NAME "test_checkpoint_v2_protocol_selector_${selector_case}"
      COMMAND "$<TARGET_FILE:test_checkpoint_v2_protocol>" ${selector_args})
    set_tests_properties(
      "test_checkpoint_v2_protocol_selector_${selector_case}" PROPERTIES
      LABELS "unit;contract;task23;task23-fast" TIMEOUT 30 WILL_FAIL TRUE)
  endforeach()

  hundun_add_test(test_checkpoint_v2_epoch
    ${PROJECT_SOURCE_DIR}/tests/unit/test_checkpoint_v2_epoch.cpp)
  target_include_directories(test_checkpoint_v2_epoch PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_checkpoint_v2_epoch PRIVATE hundun_runtime)
  set_tests_properties(test_checkpoint_v2_epoch PROPERTIES
    LABELS "unit;task23;task23-fast" TIMEOUT 30)

  hundun_add_test(test_checkpoint_v2_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_checkpoint_v2_header_contract.cpp)
  target_include_directories(test_checkpoint_v2_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_checkpoint_v2_header_contract PRIVATE
    hundun_checkpoint_diagnostics hundun_flow)
  set_tests_properties(test_checkpoint_v2_header_contract PROPERTIES
    LABELS "unit;contract;task23;task23-fast" TIMEOUT 30)

  hundun_add_test(test_checkpoint_v3_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_checkpoint_v3_header_contract.cpp)
  target_include_directories(test_checkpoint_v3_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_checkpoint_v3_header_contract PRIVATE
    hundun_flow)
  set_tests_properties(test_checkpoint_v3_header_contract PROPERTIES
    LABELS "unit;header;stage3;task17a" TIMEOUT 30)

  hundun_add_test(test_checkpoint_v3_diagnostics_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_checkpoint_v3_diagnostics_header_contract.cpp)
  target_link_libraries(test_checkpoint_v3_diagnostics_header_contract PRIVATE
    hundun_checkpoint_diagnostics)
  set_tests_properties(test_checkpoint_v3_diagnostics_header_contract PROPERTIES
    LABELS "unit;header;stage3;task18a" TIMEOUT 30)

  hundun_add_test(test_checkpoint_v3_diagnostics
    ${PROJECT_SOURCE_DIR}/tests/unit/test_checkpoint_v3_diagnostics.cpp)
  target_include_directories(test_checkpoint_v3_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_checkpoint_v3_diagnostics PRIVATE
    hundun_checkpoint_diagnostics)
  set_tests_properties(test_checkpoint_v3_diagnostics PROPERTIES
    LABELS "unit;red;stage3;task18a" TIMEOUT 30)

  hundun_add_test(test_immersed_diagnostics_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_immersed_diagnostics_header_contract.cpp)
  target_link_libraries(test_immersed_diagnostics_header_contract PRIVATE
    hundun_session_diagnostics)
  set_tests_properties(test_immersed_diagnostics_header_contract PROPERTIES
    LABELS "unit;header;stage3;task18a" TIMEOUT 30)

  hundun_add_test(test_wale_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_wale_header_contract.cpp)
  target_link_libraries(test_wale_header_contract PRIVATE hundun_les)
  set_tests_properties(test_wale_header_contract PROPERTIES
    LABELS "unit;header;stage3;task12" TIMEOUT 30)

  hundun_add_test(test_wale
    ${PROJECT_SOURCE_DIR}/tests/unit/test_wale.cpp)
  target_include_directories(test_wale PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_wale PRIVATE hundun_les)
  set_tests_properties(test_wale PROPERTIES
    LABELS "unit;red;stage3;task12" TIMEOUT 30)

  add_executable(test_wale_mpi
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_wale_mpi.cpp)
  target_include_directories(test_wale_mpi PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_wale_mpi PRIVATE
    hundun_les hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2)
    add_test(NAME "test_wale_mpi_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_wale_mpi>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_wale_mpi_${ranks}_rank" PROPERTIES
      LABELS "mpi;red;stage3;task12" TIMEOUT 60 PROCESSORS "${ranks}")
  endforeach()
  hundun_add_test(test_checkpoint_v3_codec
    ${PROJECT_SOURCE_DIR}/tests/unit/test_checkpoint_v3_codec.cpp)
  target_include_directories(test_checkpoint_v3_codec PRIVATE
    "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_checkpoint_v3_codec PRIVATE
    hundun_flow)
  set_tests_properties(test_checkpoint_v3_codec PROPERTIES
    LABELS "unit;red;stage3;task17a" TIMEOUT 30)

  add_test(NAME test_task23_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task23_source_policy.cmake")
  set_tests_properties(test_task23_source_policy PROPERTIES
    LABELS "unit;source_policy;task23;task23-fast" TIMEOUT 30)
  foreach(task23_r4_section IN ITEMS report fingerprint diagnostics)
    add_test(NAME "test_task23_r4_${task23_r4_section}_contract"
      COMMAND "${CMAKE_COMMAND}"
              "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
              "-DHUNDUN_R4_SECTION=${task23_r4_section}"
              -P
              "${PROJECT_SOURCE_DIR}/tests/cmake/task23_r4_oracle_contract.cmake")
    set_tests_properties(
      "test_task23_r4_${task23_r4_section}_contract" PROPERTIES
      LABELS "unit;source_policy;task23;task23-fast" TIMEOUT 30)
  endforeach()
  add_test(NAME test_task23_test_access_contract
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task23_test_access_contract.cmake")
  set_tests_properties(test_task23_test_access_contract PROPERTIES
    LABELS "unit;source_policy;task23;task23-fast" TIMEOUT 30)

  hundun_add_test(test_field_access ${PROJECT_SOURCE_DIR}/tests/unit/test_field_access.cpp)
  target_include_directories(test_field_access PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_field_access PRIVATE hundun_runtime)

  hundun_add_test(test_field_kernel_view
    ${PROJECT_SOURCE_DIR}/tests/unit/test_field_kernel_view.cpp)
  target_include_directories(test_field_kernel_view PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_field_kernel_view PRIVATE hundun_runtime)

  hundun_add_test(test_field_view_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_field_view_header_contract.cpp)
  target_include_directories(test_field_view_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/include")

  hundun_add_test(test_performance_artifact
    ${PROJECT_SOURCE_DIR}/tests/unit/test_performance_artifact.cpp)
  target_include_directories(test_performance_artifact PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_performance_artifact PRIVATE
    hundun_diagnostics yyjson)
  set_tests_properties(test_performance_artifact PROPERTIES
    LABELS "unit;task25;performance")

  hundun_add_test(test_task25_performance_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_task25_performance_contract.cpp)
  target_include_directories(test_task25_performance_contract PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_task25_performance_contract PRIVATE
    hundun_diagnostics hundun_execution)
  target_compile_definitions(test_task25_performance_contract PRIVATE
    HUNDUN_EXECUTION_ENABLE_TEST_ACCESS=1)
  set_tests_properties(test_task25_performance_contract PROPERTIES
    LABELS "unit;task25;performance" TIMEOUT 30)

  hundun_add_test(test_task25_performance_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_task25_performance_header_contract.cpp
    ${PROJECT_SOURCE_DIR}/tests/unit/test_task25_correctness_header_contract.cpp)
  target_include_directories(test_task25_performance_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/include")
  set_tests_properties(test_task25_performance_header_contract PROPERTIES
    LABELS "unit;task25;performance" TIMEOUT 30)

  add_test(NAME test_task25_tests_off_contract
    COMMAND "${CMAKE_COMMAND}"
      "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
      "-DHUNDUN_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
      "-DHUNDUN_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
      "-DHUNDUN_TESTS_OFF_EXECUTION=${CMAKE_BINARY_DIR}/task22-tests-off/archives/${CMAKE_STATIC_LIBRARY_PREFIX}hundun_execution${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "-DHUNDUN_TESTS_OFF_RUNTIME=${CMAKE_BINARY_DIR}/task22-tests-off/archives/${CMAKE_STATIC_LIBRARY_PREFIX}hundun_runtime${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "-DHUNDUN_TESTS_OFF_DIAGNOSTICS=${CMAKE_BINARY_DIR}/task22-tests-off/archives/${CMAKE_STATIC_LIBRARY_PREFIX}hundun_diagnostics${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "-DHUNDUN_TESTS_OFF_LINEAR=${CMAKE_BINARY_DIR}/task22-tests-off/archives/${CMAKE_STATIC_LIBRARY_PREFIX}hundun_linear${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "-DHUNDUN_TESTS_OFF_FLOW=${CMAKE_BINARY_DIR}/task22-tests-off/archives/${CMAKE_STATIC_LIBRARY_PREFIX}hundun_flow${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "-DHUNDUN_TESTS_OFF_HUNDUN=${CMAKE_BINARY_DIR}/task22-tests-off/archives/hundun${CMAKE_EXECUTABLE_SUFFIX}"
      -P "${PROJECT_SOURCE_DIR}/tests/cmake/task25_tests_off_contract.cmake")
  set_tests_properties(test_task25_tests_off_contract PROPERTIES
    FIXTURES_REQUIRED task22_tests_off_archives
    LABELS "unit;source_policy;task25"
    RUN_SERIAL TRUE
    TIMEOUT 60)

  add_test(NAME test_task25_source_metadata_contract
    COMMAND "${CMAKE_COMMAND}"
      "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
      -P "${PROJECT_SOURCE_DIR}/tests/cmake/task25_source_metadata_contract.cmake")
  set_tests_properties(test_task25_source_metadata_contract PROPERTIES
    LABELS "unit;source_policy;task25;performance" TIMEOUT 30)

  hundun_add_test(test_task25_state_fingerprint
    ${PROJECT_SOURCE_DIR}/tests/unit/test_task25_state_fingerprint.cpp)
  target_include_directories(test_task25_state_fingerprint PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_task25_state_fingerprint PRIVATE
    hundun_flow hundun_options hundun_warnings)
  set_tests_properties(test_task25_state_fingerprint PROPERTIES
    LABELS "unit;task25;performance" TIMEOUT 30)

  add_executable(task25_performance_evidence
    ${PROJECT_SOURCE_DIR}/tests/support/task25_performance_evidence.cpp)
  target_include_directories(task25_performance_evidence PRIVATE
    "${PROJECT_SOURCE_DIR}/third_party/yyjson")
  target_link_libraries(task25_performance_evidence PRIVATE
    hundun_diagnostics hundun_options hundun_warnings yyjson)
  add_test(NAME test_task25_performance_evidence_self
    COMMAND "$<TARGET_FILE:task25_performance_evidence>" --self-test)
  set_tests_properties(test_task25_performance_evidence_self PROPERTIES
    LABELS "unit;task25;performance" TIMEOUT 30)

  add_executable(task25_performance_oracle
    ${PROJECT_SOURCE_DIR}/tests/mpi/task25_performance_oracle.cpp)
  target_include_directories(task25_performance_oracle PRIVATE
    "${PROJECT_SOURCE_DIR}/third_party/yyjson")
  target_link_libraries(task25_performance_oracle PRIVATE
    hundun_diagnostics hundun_mesh hundun_runtime yyjson
    hundun_options hundun_warnings)
  add_test(NAME test_task25_performance_oracle_self
    COMMAND "$<TARGET_FILE:task25_performance_oracle>" --self-test)
  set_tests_properties(test_task25_performance_oracle_self PROPERTIES
    LABELS "unit;task25;performance" TIMEOUT 30)

  hundun_add_test(test_task25_diagnostic_branch
    ${PROJECT_SOURCE_DIR}/tests/unit/test_task25_diagnostic_branch.cpp)
  target_include_directories(test_task25_diagnostic_branch PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_task25_diagnostic_branch PRIVATE
    hundun_execution)
  set_tests_properties(test_task25_diagnostic_branch PROPERTIES
    LABELS "unit;task25;diagnostics" TIMEOUT 30)
  foreach(task25_diagnostic_ranks IN ITEMS 1 2 4)
    add_test(NAME
        "test_task25_diagnostic_branch_${task25_diagnostic_ranks}rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${task25_diagnostic_ranks}" ${MPIEXEC_PREFLAGS}
        "$<TARGET_FILE:test_task25_diagnostic_branch>" ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_task25_diagnostic_branch_${task25_diagnostic_ranks}rank"
      PROPERTIES LABELS "mpi;task25;diagnostics"
      TIMEOUT 30 PROCESSORS "${task25_diagnostic_ranks}")
  endforeach()

  foreach(task25_ranks IN ITEMS 1 2 4)
    add_test(NAME "test_task25_performance_${task25_ranks}rank"
      COMMAND bash
        "${PROJECT_SOURCE_DIR}/tests/acceptance/task25_performance.sh"
        "$<TARGET_FILE:hundun>" "${MPIEXEC_EXECUTABLE}" "${task25_ranks}")
    set_tests_properties("test_task25_performance_${task25_ranks}rank"
      PROPERTIES LABELS "acceptance;mpi;task25;performance"
      TIMEOUT 300 PROCESSORS "${task25_ranks}" RUN_SERIAL TRUE)
  endforeach()

  hundun_add_test(test_structured_diagnostics_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_structured_diagnostics_header_contract.cpp)
  target_include_directories(test_structured_diagnostics_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/include")
  set_tests_properties(test_structured_diagnostics_header_contract PROPERTIES
    LABELS "unit;task19" TIMEOUT 30)

  hundun_add_test(test_structured_diagnostics
    ${PROJECT_SOURCE_DIR}/tests/unit/test_structured_diagnostics.cpp)
  target_include_directories(test_structured_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_structured_diagnostics PRIVATE
    hundun_diagnostics_core)
  target_compile_definitions(test_structured_diagnostics PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
  set_tests_properties(test_structured_diagnostics PROPERTIES
    LABELS "unit;task19" TIMEOUT 30)

  hundun_add_test(test_material_density_transport_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_material_density_transport_header_contract.cpp)
  target_include_directories(test_material_density_transport_header_contract
    PRIVATE "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_material_density_transport_header_contract PRIVATE
    hundun_flow)
  set_tests_properties(test_material_density_transport_header_contract PROPERTIES
    LABELS "unit;task19" TIMEOUT 30)

  hundun_add_test(test_material_density_piso_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_material_density_piso_header_contract.cpp)
  target_include_directories(test_material_density_piso_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_material_density_piso_header_contract PRIVATE
    hundun_material_diagnostics hundun_flow)
  set_tests_properties(test_material_density_piso_header_contract PROPERTIES
    LABELS "unit;task20;task20-fast" TIMEOUT 30)

  hundun_add_test(test_ideal_gas_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_ideal_gas_header_contract.cpp)
  target_include_directories(test_ideal_gas_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_ideal_gas_header_contract PRIVATE
    hundun_material_diagnostics hundun_flow)
  set_tests_properties(test_ideal_gas_header_contract PROPERTIES
    LABELS "unit;task21;task21-fast" TIMEOUT 30)

  hundun_add_test(test_mpi_operation_error
    ${PROJECT_SOURCE_DIR}/tests/unit/test_mpi_operation_error.cpp)
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

  hundun_add_test(test_execution ${PROJECT_SOURCE_DIR}/tests/unit/test_execution.cpp)
  target_include_directories(test_execution PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_execution PRIVATE hundun_execution)
  set_tests_properties(test_execution PROPERTIES TIMEOUT 30)

  hundun_add_test(test_ghosted_vector
    ${PROJECT_SOURCE_DIR}/tests/unit/test_ghosted_vector.cpp)
  target_include_directories(test_ghosted_vector PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_ghosted_vector PRIVATE hundun_linear)
  set_tests_properties(test_ghosted_vector PROPERTIES TIMEOUT 30)

  hundun_add_test(test_ghosted_vector_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_ghosted_vector_header_contract.cpp)
  target_include_directories(test_ghosted_vector_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_ghosted_vector_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_ghosted_vector_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_vector_ops_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_vector_ops_header_contract.cpp)
  target_link_libraries(test_vector_ops_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_vector_ops_header_contract PROPERTIES TIMEOUT 30)

  hundun_add_test(test_vector_ops
    ${PROJECT_SOURCE_DIR}/tests/unit/test_vector_ops.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_vector_ops PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_vector_ops PRIVATE hundun_linear)
  set_tests_properties(test_vector_ops PROPERTIES TIMEOUT 30)

  hundun_add_test(test_linear_system_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_linear_system_header_contract.cpp)
  target_link_libraries(test_linear_system_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_linear_system_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_preconditioners
    ${PROJECT_SOURCE_DIR}/tests/unit/test_preconditioners.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_preconditioners PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_preconditioners PRIVATE hundun_linear)
  set_tests_properties(test_preconditioners PROPERTIES TIMEOUT 30)

  hundun_add_test(test_conjugate_gradient_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_conjugate_gradient_header_contract.cpp)
  target_link_libraries(test_conjugate_gradient_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_conjugate_gradient_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_bicgstab_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_bicgstab_header_contract.cpp)
  target_link_libraries(test_bicgstab_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_bicgstab_header_contract PROPERTIES TIMEOUT 30)

  hundun_add_test(test_restarted_gmres_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_restarted_gmres_header_contract.cpp)
  target_link_libraries(test_restarted_gmres_header_contract PRIVATE
    hundun_linear)
  set_tests_properties(test_restarted_gmres_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_matrix_free_poisson_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_matrix_free_poisson_header_contract.cpp)
  target_link_libraries(test_matrix_free_poisson_header_contract PRIVATE
    hundun_fvm)
  set_tests_properties(test_matrix_free_poisson_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_basic_boundary_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_basic_boundary_header_contract.cpp)
  target_include_directories(test_basic_boundary_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_compile_definitions(test_basic_boundary_header_contract PRIVATE
    HUNDUN_SOURCE_ROOT="${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_basic_boundary_header_contract PRIVATE
    hundun_config hundun_fvm hundun_mesh hundun_runtime)
  set_tests_properties(test_basic_boundary_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_cell_centered_fvm_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_cell_centered_fvm_header_contract.cpp)
  target_link_libraries(test_cell_centered_fvm_header_contract PRIVATE
    hundun_fvm)
  set_tests_properties(test_cell_centered_fvm_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_immersed_reconstruction_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_immersed_reconstruction_header_contract.cpp)
  target_link_libraries(test_immersed_reconstruction_header_contract PRIVATE
    hundun_fvm)
  set_tests_properties(test_immersed_reconstruction_header_contract PROPERTIES
    LABELS "unit;task7"
    TIMEOUT 30)

  hundun_add_test(test_immersed_operator_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_immersed_operator_header_contract.cpp)
  target_link_libraries(test_immersed_operator_header_contract PRIVATE
    hundun_fvm)
  set_tests_properties(test_immersed_operator_header_contract PROPERTIES
    LABELS "unit;task8"
    TIMEOUT 30)
  add_test(NAME test_cell_centered_fvm_source_policy
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/task16_source_policy.cmake")
  set_tests_properties(test_cell_centered_fvm_source_policy PROPERTIES
    LABELS "unit;source_policy"
    TIMEOUT 30)

  hundun_add_test(test_momentum_predictor_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_momentum_predictor_header_contract.cpp)
  target_link_libraries(test_momentum_predictor_header_contract PRIVATE
    hundun_flow)
  set_tests_properties(test_momentum_predictor_header_contract PROPERTIES
    TIMEOUT 30)

  hundun_add_test(test_task18_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_task18_header_contract.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/unit/test_flow_state.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_flow_state PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_flow_state PRIVATE hundun_flow)
  set_tests_properties(test_flow_state PROPERTIES TIMEOUT 30)

  hundun_add_test(test_adaptive_time_control_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_adaptive_time_control_header_contract.cpp)
  target_include_directories(test_adaptive_time_control_header_contract PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_adaptive_time_control_header_contract PRIVATE
    hundun_flow)
  set_tests_properties(test_adaptive_time_control_header_contract PROPERTIES
    LABELS "unit;task22;task22-fast" TIMEOUT 30)

  hundun_add_test(test_adaptive_time_control
    ${PROJECT_SOURCE_DIR}/tests/unit/test_adaptive_time_control.cpp)
  target_include_directories(test_adaptive_time_control PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_adaptive_time_control PRIVATE hundun_flow)
  set_tests_properties(test_adaptive_time_control PROPERTIES
    LABELS "unit;task22;task22-fast" TIMEOUT 30)

  add_executable(test_adaptive_time_control_mpi
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_adaptive_time_control.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
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

  add_executable(test_restarted_gmres
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_restarted_gmres.cpp)
  target_include_directories(test_restarted_gmres PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_restarted_gmres PRIVATE
    hundun_linear hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(gmres_test "test_restarted_gmres_${ranks}_rank")
    add_test(NAME "${gmres_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_restarted_gmres>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${gmres_test}" PROPERTIES
      LABELS "mpi;linear;nonsymmetric;red;fast" TIMEOUT 60
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_wale_body_fitted
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_wale_body_fitted.cpp)
  target_include_directories(test_wale_body_fitted PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_wale_body_fitted PRIVATE
    hundun_flow hundun_les hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2)
    add_test(NAME "test_wale_body_fitted_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_wale_body_fitted>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_wale_body_fitted_${ranks}_rank" PROPERTIES
      LABELS "mpi;red;stage3;task13;fast" TIMEOUT 120 PROCESSORS "${ranks}")
  endforeach()
  foreach(cells IN ITEMS 12 24)
    add_test(NAME "test_wale_taylor_green_${cells}_screen"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_wale_body_fitted>"
              "${cells}" ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_wale_taylor_green_${cells}_screen" PROPERTIES
      LABELS "numerical;screen;stage3;task13;wale-milestone"
      TIMEOUT 300 PROCESSORS 1)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_adaptive_time_control_diagnostics.cpp)
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
  foreach(ranks IN ITEMS 1 2 4)
    add_test(NAME "test_checkpoint_v3_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_transaction>"
              checkpoint_v3 ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_checkpoint_v3_${ranks}_rank" PROPERTIES
      LABELS "mpi;numerical;red;stage3;task17a" TIMEOUT 900
      PROCESSORS "${ranks}")
  endforeach()
  add_test(NAME test_checkpoint_v3_collective_preflight_2_rank
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_transaction>"
            checkpoint_v3_collective_preflight ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_checkpoint_v3_collective_preflight_2_rank
    PROPERTIES LABELS "mpi;red;stage3;task17a" TIMEOUT 30 PROCESSORS 2)
  add_test(NAME test_checkpoint_v3_invalid_state_1_rank
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_transaction>"
            checkpoint_v3_invalid_state ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_checkpoint_v3_invalid_state_1_rank
    PROPERTIES LABELS "mpi;red;stage3;task17a" TIMEOUT 30 PROCESSORS 1)
  add_test(NAME test_checkpoint_v3_path_agreement_2_rank
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_transaction>"
            checkpoint_v3_path_agreement ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_checkpoint_v3_path_agreement_2_rank
    PROPERTIES LABELS "mpi;red;stage3;task17a" TIMEOUT 30 PROCESSORS 2)
  add_test(NAME test_checkpoint_v3_active_attempt_1_rank
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_transaction>"
            checkpoint_v3_active_attempt ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_checkpoint_v3_active_attempt_1_rank
    PROPERTIES LABELS "mpi;red;stage3;task17a" TIMEOUT 30 PROCESSORS 1)
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
    ${PROJECT_SOURCE_DIR}/tests/unit/test_momentum_predictor.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_time_consistent_face_velocity.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_momentum_predictor_collective.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_fixed_step_piso.cpp)
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
      LABELS "mpi;numerical;task18;task25"
      TIMEOUT 120
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_material_density_transport
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_material_density_transport.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_material_density_piso.cpp)
  target_include_directories(test_material_density_piso PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_material_density_piso PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  add_executable(test_variable_density_vortex
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_variable_density_vortex.cpp)
  target_include_directories(test_variable_density_vortex PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_variable_density_vortex PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  add_executable(test_material_density_piso_diagnostics
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_material_density_piso_diagnostics.cpp)
  target_include_directories(test_material_density_piso_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_material_density_piso_diagnostics PRIVATE
    hundun_material_diagnostics hundun_flow hundun_options hundun_warnings
    MPI::MPI_CXX)
  target_compile_definitions(test_material_density_piso_diagnostics PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
  add_executable(test_ideal_gas_closure
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_ideal_gas_closure.cpp)
  target_include_directories(test_ideal_gas_closure PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_ideal_gas_closure PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  add_executable(test_ideal_gas_piso
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_ideal_gas_piso.cpp)
  target_include_directories(test_ideal_gas_piso PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_ideal_gas_piso PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  add_executable(test_ideal_gas_diagnostics
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_ideal_gas_diagnostics.cpp)
  target_include_directories(test_ideal_gas_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_ideal_gas_diagnostics PRIVATE
    hundun_material_diagnostics hundun_flow hundun_options hundun_warnings
    MPI::MPI_CXX)
  target_compile_definitions(test_ideal_gas_diagnostics PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
  add_executable(test_checkpoint_v2
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_checkpoint_v2.cpp)
  target_include_directories(test_checkpoint_v2 PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/src")
  target_link_libraries(test_checkpoint_v2 PRIVATE
    hundun_checkpoint_diagnostics hundun_flow hundun_options hundun_warnings
    MPI::MPI_CXX)
  target_compile_definitions(test_checkpoint_v2 PRIVATE
    HUNDUN_RUNTIME_ENABLE_TEST_ACCESS=1)
  foreach(ranks IN ITEMS 1 2 4)
    set(task23_checkpoint_test "test_checkpoint_v2_${ranks}_rank")
    add_test(NAME "${task23_checkpoint_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_checkpoint_v2>"
              fast
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${task23_checkpoint_test}" PROPERTIES
      LABELS "mpi;numerical;task23;task23-fast"
      TIMEOUT 180 PROCESSORS "${ranks}" RUN_SERIAL TRUE)
    set(task23_checkpoint_acceptance
        "test_checkpoint_v2_acceptance_${ranks}_rank")
    add_test(NAME "${task23_checkpoint_acceptance}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_checkpoint_v2>"
              acceptance
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${task23_checkpoint_acceptance}" PROPERTIES
      LABELS
        "mpi;numerical;task23;task23-acceptance;task23-cluster-b;task23-cluster-c;task23-cluster-d"
      TIMEOUT 1200 PROCESSORS "${ranks}" RUN_SERIAL TRUE)
  endforeach()
  foreach(selector_case IN ITEMS missing unknown extra)
    if(selector_case STREQUAL "missing")
      set(selector_args)
    elseif(selector_case STREQUAL "unknown")
      set(selector_args "not-a-task23-mode")
    else()
      set(selector_args "acceptance" "unexpected-extra-argument")
    endif()
    add_test(NAME "test_checkpoint_v2_selector_${selector_case}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_checkpoint_v2>"
              ${selector_args} ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_checkpoint_v2_selector_${selector_case}"
      PROPERTIES LABELS "mpi;contract;task23;task23-fast"
      TIMEOUT 30 PROCESSORS 1 RUN_SERIAL TRUE WILL_FAIL TRUE)
  endforeach()
  add_executable(test_checkpoint_v2_diagnostics
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_checkpoint_v2_diagnostics.cpp)
  target_include_directories(test_checkpoint_v2_diagnostics PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_checkpoint_v2_diagnostics PRIVATE
    hundun_checkpoint_diagnostics hundun_flow hundun_options hundun_warnings
    MPI::MPI_CXX)
  target_compile_definitions(test_checkpoint_v2_diagnostics PRIVATE
    HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS=1)
  foreach(ranks IN ITEMS 1 2 4)
    set(task23_diagnostics_test
        "test_checkpoint_v2_diagnostics_${ranks}_rank")
    add_test(NAME "${task23_diagnostics_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_checkpoint_v2_diagnostics>"
              fast
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
              acceptance
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${task23_diagnostics_acceptance}" PROPERTIES
      LABELS "mpi;diagnostics;task23;task23-acceptance;task23-cluster-e"
      TIMEOUT 300 PROCESSORS "${ranks}" RUN_SERIAL TRUE)
  endforeach()
  foreach(selector_case IN ITEMS missing unknown extra)
    if(selector_case STREQUAL "missing")
      set(selector_args)
    elseif(selector_case STREQUAL "unknown")
      set(selector_args "not-a-task23-mode")
    else()
      set(selector_args "acceptance" "unexpected-extra-argument")
    endif()
    add_test(NAME "test_checkpoint_v2_diagnostics_selector_${selector_case}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_checkpoint_v2_diagnostics>"
              ${selector_args} ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_checkpoint_v2_diagnostics_selector_${selector_case}" PROPERTIES
      LABELS "mpi;diagnostics;contract;task23;task23-fast"
      TIMEOUT 30 PROCESSORS 1 RUN_SERIAL TRUE WILL_FAIL TRUE)
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
      if(target STREQUAL "material_density_piso")
        set(task20_labels "mpi;numerical;task20;task20-fast;task25")
      else()
        set(task20_labels "mpi;numerical;task20;task20-fast")
      endif()
      set_tests_properties("${task20_test}" PROPERTIES
        LABELS "${task20_labels}" TIMEOUT 240 PROCESSORS "${ranks}")
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_material_density_wave.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_material_density_diagnostics.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_piso_decomposition.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_taylor_green_piso.cpp)
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
      LABELS "mpi;numerical;task18;task25"
      TIMEOUT 300
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_basic_boundary
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_basic_boundary.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_basic_boundary PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_basic_boundary PRIVATE
    hundun_config hundun_fvm hundun_mesh hundun_runtime
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_cell_centered_fvm.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_cell_centered_fvm PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_cell_centered_fvm PRIVATE
    hundun_fvm hundun_options hundun_warnings MPI::MPI_CXX)
  set(cell_centered_fvm_timeout 300)
  if(HUNDUN_ENABLE_ASAN OR HUNDUN_ENABLE_UBSAN)
    set(cell_centered_fvm_timeout 600)
  endif()
  foreach(ranks IN ITEMS 1 2 4)
    set(cell_centered_fvm_test "test_cell_centered_fvm_${ranks}_rank")
    add_test(NAME "${cell_centered_fvm_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_cell_centered_fvm>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${cell_centered_fvm_test}" PROPERTIES
      LABELS "mpi;numerical;task25"
      TIMEOUT "${cell_centered_fvm_timeout}"
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_immersed_reconstruction
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_immersed_reconstruction.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_immersed_reconstruction PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_immersed_reconstruction PRIVATE
    hundun_fvm stage3_stl_fixture hundun_options hundun_warnings
    MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    foreach(mode IN ITEMS uniform warped periodic symmetry open)
      set(immersed_reconstruction_test
          "test_immersed_reconstruction_${mode}_${ranks}_rank")
      add_test(NAME "${immersed_reconstruction_test}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_immersed_reconstruction>" "${mode}"
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties("${immersed_reconstruction_test}" PROPERTIES
        LABELS "mpi;numerical;task7"
        TIMEOUT 300
        PROCESSORS "${ranks}")
    endforeach()
  endforeach()

  add_executable(test_immersed_operator
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_immersed_operator.cpp)
  target_include_directories(test_immersed_operator PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_immersed_operator PRIVATE
    hundun_fvm stage3_stl_fixture hundun_options hundun_warnings
    MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    foreach(mode IN ITEMS uniform warped)
      set(immersed_operator_test
          "test_immersed_operator_${mode}_${ranks}_rank")
      add_test(NAME "${immersed_operator_test}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_immersed_operator>" "${mode}"
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties("${immersed_operator_test}" PROPERTIES
        LABELS "mpi;numerical;task8"
        TIMEOUT 300
        PROCESSORS "${ranks}")
    endforeach()
  endforeach()

  hundun_add_test(test_stage3_flow_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_stage3_flow_header_contract.cpp)
  target_link_libraries(test_stage3_flow_header_contract PRIVATE hundun_flow)
  set_tests_properties(test_stage3_flow_header_contract PROPERTIES
    LABELS "unit;stage3;task9")

  hundun_add_test(test_wall_force_header_contract
    ${PROJECT_SOURCE_DIR}/tests/unit/test_wall_force_header_contract.cpp)
  target_link_libraries(test_wall_force_header_contract PRIVATE
    hundun_flow hundun_immersed)
  set_tests_properties(test_wall_force_header_contract PROPERTIES
    LABELS "unit;stage3;task10")

  add_executable(test_wall_force ${PROJECT_SOURCE_DIR}/tests/mpi/test_wall_force.cpp)
  target_include_directories(test_wall_force PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_wall_force PRIVATE
    hundun_immersed stage3_stl_fixture hundun_options hundun_warnings
    MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    add_test(NAME "test_wall_force_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_wall_force>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_wall_force_${ranks}_rank" PROPERTIES
      LABELS "mpi;numerical;stage3;task10"
      TIMEOUT 300
      PROCESSORS "${ranks}")
  endforeach()
  foreach(ranks IN ITEMS 2 4)
    add_test(NAME "test_wall_force_${ranks}_rank_alternate"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_wall_force>" alternate
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_wall_force_${ranks}_rank_alternate" PROPERTIES
      LABELS "mpi;numerical;stage3;task10"
      TIMEOUT 300
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_immersed_pressure_operator
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_immersed_pressure_operator.cpp)
  target_include_directories(test_immersed_pressure_operator PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_immersed_pressure_operator PRIVATE
    hundun_flow hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    add_test(NAME "test_immersed_pressure_operator_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_immersed_pressure_operator>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_immersed_pressure_operator_${ranks}_rank" PROPERTIES
      LABELS "mpi;numerical;stage3;task9"
      TIMEOUT 120
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_immersed_piso
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_immersed_piso.cpp)
  target_include_directories(test_immersed_piso PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_immersed_piso PRIVATE
    hundun_flow stage3_stl_fixture hundun_options hundun_warnings
    MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    add_test(NAME "test_immersed_piso_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_piso>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_immersed_piso_${ranks}_rank" PROPERTIES
      LABELS "mpi;numerical;stage3;task9"
      TIMEOUT 1800
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_immersed_transaction
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_immersed_transaction.cpp)
  target_include_directories(test_immersed_transaction PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_immersed_transaction PRIVATE
    hundun_checkpoint_diagnostics hundun_session_diagnostics hundun_flow
    stage3_stl_fixture hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    add_test(NAME "test_immersed_transaction_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_transaction>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_immersed_transaction_${ranks}_rank" PROPERTIES
      LABELS "mpi;numerical;stage3;task11"
      TIMEOUT 1800
      PROCESSORS "${ranks}")
  endforeach()

  foreach(ranks IN ITEMS 1 2)
    add_test(NAME "test_immersed_diagnostics_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_transaction>"
              diagnostics ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_immersed_diagnostics_${ranks}_rank" PROPERTIES
      LABELS "mpi;red;stage3;task18a" TIMEOUT 900 PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_task11_signed_force_authority
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_task11_signed_force_authority.cpp)
  target_include_directories(test_task11_signed_force_authority PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_task11_signed_force_authority PRIVATE
    hundun_flow stage3_stl_fixture hundun_options hundun_warnings
    MPI::MPI_CXX)
  target_compile_definitions(test_task11_signed_force_authority PRIVATE
    HUNDUN_IMMERSED_ENABLE_TEST_ACCESS=1)
  foreach(ranks IN ITEMS 1 2 4)
    add_test(NAME "fixture_preflight_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_task11_signed_force_authority>"
              fixture_preflight ${MPIEXEC_POSTFLAGS})
    set_tests_properties("fixture_preflight_${ranks}_rank" PROPERTIES
      LABELS "mpi;numerical;stage3;task11;preflight"
      TIMEOUT 300
      PROCESSORS "${ranks}")
    add_test(NAME "test_task11_signed_force_authority_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_task11_signed_force_authority>"
              behavioral_red ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      "test_task11_signed_force_authority_${ranks}_rank" PROPERTIES
      LABELS "mpi;numerical;stage3;task11;red"
      TIMEOUT 300
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_laminar_ibm_order
    ${PROJECT_SOURCE_DIR}/tests/numerical/test_laminar_ibm_order.cpp)
  target_include_directories(test_laminar_ibm_order PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_laminar_ibm_order PRIVATE
    hundun_flow stage3_mms stage3_stl_fixture hundun_options hundun_warnings
    MPI::MPI_CXX)
  add_test(NAME test_laminar_ibm_order_oracle
    COMMAND "$<TARGET_FILE:test_laminar_ibm_order>" oracle)
  set_tests_properties(test_laminar_ibm_order_oracle PROPERTIES
    LABELS "unit;numerical;stage3;task11"
    TIMEOUT 30)
  add_test(NAME test_laminar_ibm_signed_force_reference_preflight
    COMMAND "$<TARGET_FILE:test_laminar_ibm_order>"
            signed_force_reference_preflight)
  set_tests_properties(test_laminar_ibm_signed_force_reference_preflight
    PROPERTIES
      LABELS "unit;numerical;stage3;task11;fast"
      TIMEOUT 30)
  add_test(NAME test_laminar_ibm_exact_force_consistency_fast
    COMMAND "$<TARGET_FILE:test_laminar_ibm_order>"
            exact_force_consistency_fast)
  set_tests_properties(test_laminar_ibm_exact_force_consistency_fast
    PROPERTIES
      LABELS "numerical;stage3;task11;fast"
      TIMEOUT 1800)
  add_test(NAME test_laminar_ibm_surface_acceptance
    COMMAND "$<TARGET_FILE:test_laminar_ibm_order>" surface_acceptance)
  set_tests_properties(test_laminar_ibm_surface_acceptance PROPERTIES
    LABELS "unit;numerical;stage3;task11"
    TIMEOUT 120)
  add_test(NAME test_laminar_ibm_order_sphere_fast
    COMMAND "$<TARGET_FILE:test_laminar_ibm_order>" sphere_fast)
  set_tests_properties(test_laminar_ibm_order_sphere_fast PROPERTIES
    LABELS "numerical;stage3;task11;fast"
    TIMEOUT 900)
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_test(NAME test_laminar_ibm_functional_selection_fast
      COMMAND "$<TARGET_FILE:test_laminar_ibm_order>"
              functional_selection_fast)
    set_tests_properties(test_laminar_ibm_functional_selection_fast PROPERTIES
      LABELS "numerical;stage3;task11;fast"
      TIMEOUT 1800
      RUN_SERIAL TRUE)
    add_test(NAME test_laminar_ibm_functional_closure_fast
      COMMAND "$<TARGET_FILE:test_laminar_ibm_order>"
              functional_closure_fast)
    set_tests_properties(test_laminar_ibm_functional_closure_fast PROPERTIES
      LABELS "numerical;stage3;task11;fast"
      TIMEOUT 7200
      RUN_SERIAL TRUE)
    add_test(NAME test_laminar_ibm_near_wall_pressure_diagnostic
      COMMAND "$<TARGET_FILE:test_laminar_ibm_order>"
              near_wall_pressure_diagnostic)
    set_tests_properties(test_laminar_ibm_near_wall_pressure_diagnostic
      PROPERTIES
      LABELS "numerical;stage3;task11;diagnostic"
      TIMEOUT 7200
      RUN_SERIAL TRUE)
    foreach(selector IN ITEMS
        sphere_uniform_acceptance
        sphere_warped_acceptance
        prism_warped_acceptance
        )
      add_test(NAME "test_laminar_ibm_order_${selector}"
        COMMAND "$<TARGET_FILE:test_laminar_ibm_order>" "${selector}")
      set_tests_properties("test_laminar_ibm_order_${selector}" PROPERTIES
        LABELS "numerical;acceptance;stage3;task11"
        TIMEOUT 21600
        RUN_SERIAL TRUE)
    endforeach()
  endif()

  add_executable(test_laminar_ibm_decomposition
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_laminar_ibm_decomposition.cpp)
  target_include_directories(test_laminar_ibm_decomposition PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_laminar_ibm_decomposition PRIVATE
    stage3_mms stage3_stl_fixture hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(decomposition_case IN ITEMS
      "1:1:1:1"
      "2:2:1:1"
      "2:1:2:1"
      "4:4:1:1"
      "4:2:2:1")
    string(REPLACE ":" ";" decomposition_parts "${decomposition_case}")
    list(GET decomposition_parts 0 ranks)
    list(GET decomposition_parts 1 px)
    list(GET decomposition_parts 2 py)
    list(GET decomposition_parts 3 pz)
    set(decomposition_test
      "test_laminar_ibm_decomposition_${ranks}_rank_${px}x${py}x${pz}_fast")
    add_test(NAME "${decomposition_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_laminar_ibm_decomposition>"
              fast "${px}" "${py}" "${pz}"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${decomposition_test}" PROPERTIES
      LABELS "mpi;numerical;stage3;task11;fast"
      TIMEOUT 1800
      PROCESSORS "${ranks}")
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
      set(decomposition_acceptance
        "test_laminar_ibm_decomposition_${ranks}_rank_${px}x${py}x${pz}_acceptance")
      add_test(NAME "${decomposition_acceptance}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_laminar_ibm_decomposition>"
                acceptance "${px}" "${py}" "${pz}"
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties("${decomposition_acceptance}" PROPERTIES
        LABELS "mpi;numerical;acceptance;stage3;task11"
        TIMEOUT 3600
        PROCESSORS "${ranks}"
        RUN_SERIAL TRUE)
    endif()
  endforeach()

  add_executable(test_laminar_ibm_engineering
    ${PROJECT_SOURCE_DIR}/tests/numerical/test_laminar_ibm_engineering.cpp)
  target_include_directories(test_laminar_ibm_engineering PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_laminar_ibm_engineering PRIVATE
    hundun_flow stage3_mms stage3_stl_fixture hundun_options hundun_warnings
    MPI::MPI_CXX)
  add_test(NAME test_laminar_ibm_engineering_fast_smoke
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:test_laminar_ibm_engineering>" smoke
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_laminar_ibm_engineering_fast_smoke PROPERTIES
    LABELS "mpi;numerical;stage3;task11;fast"
    TIMEOUT 1800
    PROCESSORS 1)
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    foreach(ranks IN ITEMS 1 2 4)
      add_test(NAME
        "test_laminar_ibm_engineering_${ranks}_rank_decomposition_fast"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_laminar_ibm_engineering>"
                decomposition_fast
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        "test_laminar_ibm_engineering_${ranks}_rank_decomposition_fast"
        PROPERTIES
        LABELS "mpi;numerical;stage3;task11;fast"
        TIMEOUT 21600
        PROCESSORS "${ranks}"
        RUN_SERIAL TRUE)
    endforeach()
    foreach(selector IN ITEMS cylinder_acceptance sphere_acceptance)
      add_test(NAME "test_laminar_ibm_engineering_${selector}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_laminar_ibm_engineering>" "${selector}"
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        "test_laminar_ibm_engineering_${selector}" PROPERTIES
        LABELS "mpi;numerical;acceptance;stage3;task11"
        TIMEOUT 43200
        PROCESSORS 4
        RUN_SERIAL TRUE)
    endforeach()
  endif()

  add_executable(test_immersed_closed_transient
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_immersed_closed_transient.cpp)
  target_include_directories(test_immersed_closed_transient PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_immersed_closed_transient PRIVATE
    hundun_flow stage3_mms stage3_stl_fixture hundun_options hundun_warnings
    MPI::MPI_CXX)
  add_test(NAME test_immersed_closed_transient_smoke
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:test_immersed_closed_transient>" smoke
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_immersed_closed_transient_smoke PROPERTIES
    LABELS "mpi;numerical;stage3;task11;fast"
    TIMEOUT 900
    PROCESSORS 1)
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    foreach(ranks IN ITEMS 1 2 4)
      add_test(NAME "test_immersed_closed_transient_${ranks}_rank_acceptance"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_immersed_closed_transient>" acceptance
                ${MPIEXEC_POSTFLAGS})
      set_tests_properties(
        "test_immersed_closed_transient_${ranks}_rank_acceptance" PROPERTIES
        LABELS "mpi;numerical;acceptance;stage3;task11"
        TIMEOUT 21600
        PROCESSORS "${ranks}"
        RUN_SERIAL TRUE)
    endforeach()
  endif()

  add_executable(test_matrix_free_poisson
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_matrix_free_poisson.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
  target_include_directories(test_matrix_free_poisson PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_matrix_free_poisson PRIVATE
    hundun_fvm hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(matrix_free_poisson_test
        "test_matrix_free_poisson_${ranks}_rank")
    add_test(NAME "${matrix_free_poisson_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_matrix_free_poisson>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${matrix_free_poisson_test}" PROPERTIES
      LABELS "mpi;numerical;task25"
      TIMEOUT 90
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(test_poisson_sine
    ${PROJECT_SOURCE_DIR}/tests/numerical/test_poisson_sine.cpp)
  target_include_directories(test_poisson_sine PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_poisson_sine PRIVATE
    hundun_fvm hundun_options hundun_warnings MPI::MPI_CXX)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_conjugate_gradient.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_bicgstab.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_vector_ops_mpi.cpp
    ${PROJECT_SOURCE_DIR}/tests/support/allocation_attempt_guard.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_ghosted_vector_halo.cpp)
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
      LABELS "mpi;task25" TIMEOUT 60 PROCESSORS "${ranks}")

    set(sparse_vector_halo_test
      "test_ghosted_vector_halo_sparse_${ranks}_rank")
    add_test(NAME "${sparse_vector_halo_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_ghosted_vector_halo>"
              sparse ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${sparse_vector_halo_test}" PROPERTIES
      LABELS "mpi;task25;periodic_ibm;red;fast" TIMEOUT 60
      PROCESSORS "${ranks}")
  endforeach()
  add_test(NAME test_ghosted_vector_halo_sparse_failures
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_ghosted_vector_halo>"
            sparse_failures ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_ghosted_vector_halo_sparse_failures PROPERTIES
    LABELS "mpi;task25;periodic_ibm;red;fast;mutation" TIMEOUT 60
    PROCESSORS 4)
  foreach(mode IN ITEMS mismatch failure)
    set(vector_halo_mode_test "test_ghosted_vector_halo_${mode}")
    add_test(NAME "${vector_halo_mode_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_ghosted_vector_halo>"
              "${mode}" ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${vector_halo_mode_test}" PROPERTIES
      LABELS "mpi;task25" TIMEOUT 60 PROCESSORS 4)
  endforeach()
  add_test(NAME test_ghosted_vector_halo_finalized_idle
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_ghosted_vector_halo>"
            finalized_idle ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_ghosted_vector_halo_finalized_idle PROPERTIES
    LABELS "mpi;task25" TIMEOUT 30 PROCESSORS 1)

  add_executable(test_passive_scalar ${PROJECT_SOURCE_DIR}/tests/unit/test_passive_scalar.cpp)
  target_include_directories(test_passive_scalar PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_passive_scalar_mpi.cpp)
  target_include_directories(test_passive_scalar_mpi PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_passive_scalar_mpi PRIVATE
    hundun_solver hundun_options hundun_warnings MPI::MPI_CXX)
  target_compile_definitions(test_passive_scalar_mpi PRIVATE
    HUNDUN_RUNTIME_ENABLE_TEST_ACCESS=1)
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

  hundun_add_test(test_restart ${PROJECT_SOURCE_DIR}/tests/unit/test_restart.cpp)
  target_include_directories(test_restart PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_restart PRIVATE hundun_io)
  set_tests_properties(test_restart PROPERTIES TIMEOUT 30)

  hundun_add_mpi_test(test_vtk 1 ${PROJECT_SOURCE_DIR}/tests/unit/test_vtk.cpp)
  target_include_directories(test_vtk PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_vtk PRIVATE hundun_io)
  set_tests_properties(test_vtk PROPERTIES TIMEOUT 30 PROCESSORS 1)

  add_executable(test_restart_checkpoint
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_restart_checkpoint.cpp)
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

  hundun_add_test(test_halo_helpers ${PROJECT_SOURCE_DIR}/tests/unit/test_halo_helpers.cpp)
  target_include_directories(test_halo_helpers PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_halo_helpers PRIVATE hundun_runtime)
  set_tests_properties(test_halo_helpers PROPERTIES TIMEOUT 20)

  hundun_add_mpi_test(test_mesh 2 ${PROJECT_SOURCE_DIR}/tests/unit/test_mesh.cpp)
  target_include_directories(test_mesh PRIVATE "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_mesh PRIVATE hundun_mesh)
  set_tests_properties(test_mesh PROPERTIES
    TIMEOUT 30
    PROCESSORS 2)

  add_executable(test_mesh_topology ${PROJECT_SOURCE_DIR}/tests/mpi/test_mesh_topology.cpp)
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

  add_executable(test_mesh_geometry ${PROJECT_SOURCE_DIR}/tests/mpi/test_mesh_geometry.cpp)
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
      LABELS "mpi;numerical;task25"
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
    LABELS "mpi;numerical;task25"
    TIMEOUT 10
    PROCESSORS 2
    WILL_FAIL TRUE)

  add_executable(test_halo_exchange ${PROJECT_SOURCE_DIR}/tests/mpi/test_halo_exchange.cpp)
  target_include_directories(test_halo_exchange PRIVATE
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(test_halo_exchange PRIVATE
    hundun_options hundun_warnings hundun_runtime hundun_execution
    MPI::MPI_CXX)
  target_compile_definitions(test_halo_exchange PRIVATE
    HUNDUN_RUNTIME_ENABLE_TEST_ACCESS=1)
  foreach(ranks IN ITEMS 1 2 4)
    set(halo_test "test_halo_exchange_${ranks}_rank")
    add_test(NAME "${halo_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_halo_exchange>"
              full ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${halo_test}" PROPERTIES
      LABELS "mpi;task25"
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
      LABELS "mpi;task25"
      TIMEOUT 60
      PROCESSORS 4)
  endforeach()

  add_test(NAME test_halo_exchange_finalized_idle
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "1"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_halo_exchange>"
            finalized_idle ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_halo_exchange_finalized_idle PROPERTIES
    LABELS "mpi;task25"
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_collective_status.cpp)
  target_include_directories(test_collective_status PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_collective_status PRIVATE hundun_runtime)
  set_tests_properties(test_collective_status PROPERTIES TIMEOUT 20)

  hundun_add_mpi_test(test_mpi_context 4
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_mpi_context.cpp)
  target_include_directories(test_mpi_context PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_mpi_context PRIVATE hundun_runtime)
  set_tests_properties(test_mpi_context PROPERTIES TIMEOUT 30 PROCESSORS 4)

  hundun_add_mpi_test(test_mpi_external_lifecycle 1
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_mpi_external_lifecycle.cpp)
  target_include_directories(test_mpi_external_lifecycle PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(test_mpi_external_lifecycle PRIVATE hundun_runtime)
  set_tests_properties(test_mpi_external_lifecycle PROPERTIES TIMEOUT 20)

  add_executable(test_case_config_broadcast
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_case_config_broadcast.cpp)
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
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_resolved_case_broadcast.cpp)
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

  add_executable(test_resolved_case_v3_broadcast
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_resolved_case_v3_broadcast.cpp)
  target_include_directories(test_resolved_case_v3_broadcast PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_compile_definitions(test_resolved_case_v3_broadcast PRIVATE
    HUNDUN_APPLICATION_ENABLE_TEST_ACCESS=1)
  target_link_libraries(test_resolved_case_v3_broadcast PRIVATE
    hundun_application hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2 4)
    set(resolved_v3_broadcast_test
      "test_resolved_case_v3_broadcast_${ranks}_rank")
    add_test(NAME "${resolved_v3_broadcast_test}"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_resolved_case_v3_broadcast>"
              success ${MPIEXEC_POSTFLAGS})
    set_tests_properties("${resolved_v3_broadcast_test}" PROPERTIES
      LABELS "mpi;stage3;task1" TIMEOUT 45 PROCESSORS "${ranks}")
  endforeach()
  add_test(NAME test_resolved_case_v3_broadcast_failures
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:test_resolved_case_v3_broadcast>"
            failures ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_resolved_case_v3_broadcast_failures PROPERTIES
    LABELS "mpi;stage3;task1" TIMEOUT 45 PROCESSORS 4)

  add_executable(test_resolved_case_v4_broadcast
    ${PROJECT_SOURCE_DIR}/tests/mpi/test_resolved_case_v4_broadcast.cpp)
  target_include_directories(test_resolved_case_v4_broadcast PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_compile_definitions(test_resolved_case_v4_broadcast PRIVATE
    HUNDUN_APPLICATION_ENABLE_TEST_ACCESS=1)
  target_link_libraries(test_resolved_case_v4_broadcast PRIVATE
    hundun_application hundun_options hundun_warnings MPI::MPI_CXX)
  foreach(ranks IN ITEMS 1 2)
    foreach(mode IN ITEMS success failures)
      set(resolved_v4_broadcast_test
        "test_resolved_case_v4_broadcast_${mode}_${ranks}_rank")
      add_test(NAME "${resolved_v4_broadcast_test}"
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:test_resolved_case_v4_broadcast>"
                "${mode}" ${MPIEXEC_POSTFLAGS})
      set_tests_properties("${resolved_v4_broadcast_test}" PROPERTIES
        LABELS "mpi;stage4;task4f4" TIMEOUT 45 PROCESSORS "${ranks}")
    endforeach()
  endforeach()

  add_test(NAME test_stage2_task2_dispatch
    COMMAND bash
            "${PROJECT_SOURCE_DIR}/tests/acceptance/stage2_task2_dispatch.sh"
            "$<TARGET_FILE:hundun>" "${MPIEXEC_EXECUTABLE}")
  set_tests_properties(test_stage2_task2_dispatch PROPERTIES
    LABELS acceptance TIMEOUT 45 PROCESSORS 2)

  add_executable(write_task19a_stl
    ${PROJECT_SOURCE_DIR}/tests/support/write_task19a_stl.cpp)
  target_include_directories(write_task19a_stl PRIVATE
    "${PROJECT_SOURCE_DIR}")
  target_link_libraries(write_task19a_stl PRIVATE
    stage3_stl_fixture hundun_options hundun_warnings)

  foreach(ranks IN ITEMS 1 2 4)
    set(task19a_dispatch_test
      "test_task19a_immersed_flow_dispatch_${ranks}_rank")
    add_test(NAME "${task19a_dispatch_test}"
      COMMAND bash
              "${PROJECT_SOURCE_DIR}/tests/acceptance/task19a_immersed_flow_dispatch.sh"
              "$<TARGET_FILE:hundun>" "$<TARGET_FILE:write_task19a_stl>"
              "${MPIEXEC_EXECUTABLE}" "${ranks}")
    set_tests_properties("${task19a_dispatch_test}" PROPERTIES
      LABELS "acceptance;mpi;stage3;task19a;task19b;task17a" TIMEOUT 120
      PROCESSORS "${ranks}")
  endforeach()

  add_executable(compare_scalar_vtk
    ${PROJECT_SOURCE_DIR}/tests/acceptance/compare_scalar_vtk.cpp)
  target_link_libraries(compare_scalar_vtk PRIVATE
    hundun_options hundun_warnings)

  add_test(NAME acceptance_complete_gate
    COMMAND "${CMAKE_COMMAND}"
            -DHUNDUN_PROVENANCE_CASE=verify_acceptance_complete_gate
            -P "${PROJECT_SOURCE_DIR}/tests/cmake/provenance_fixture.cmake")
  set_tests_properties(acceptance_complete_gate PROPERTIES
    LABELS acceptance
    TIMEOUT 10)

  add_executable(test_decomposition ${PROJECT_SOURCE_DIR}/tests/mpi/test_decomposition.cpp)
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
