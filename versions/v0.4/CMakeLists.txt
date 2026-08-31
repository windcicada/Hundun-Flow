# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

project(HundunFlow VERSION 0.4.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(HUNDUN_BUILD_TESTS "Build HUNDUN-FLOW tests" OFF)
option(HUNDUN_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(HUNDUN_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(HUNDUN_ENABLE_HYPRE "Enable the isolated HYPRE Struct adapter" OFF)

find_package(MPI 3 REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)

if(HUNDUN_ENABLE_HYPRE)
  find_package(HYPRE REQUIRED)
endif()

get_filename_component(hundun_repository_root
  "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(hundun_v04_yyjson_root "${hundun_repository_root}/third_party/yyjson")
foreach(hundun_v04_yyjson_file IN ITEMS yyjson.c yyjson.h UPSTREAM.json)
  if(NOT EXISTS "${hundun_v04_yyjson_root}/${hundun_v04_yyjson_file}")
    message(FATAL_ERROR
      "v0.4 requires bundled yyjson at the fixed repository path: "
      "${hundun_v04_yyjson_root}/${hundun_v04_yyjson_file}")
  endif()
endforeach()

# Runtime Evidence V6 binds an exact executable digest to canonical build
# material.  This generated header contains no path or timestamp: those are
# locations/observations, not candidate identity.  The executable SHA-256 is
# computed over /proc/self/exe by the product at run start.
set(HUNDUN_SOURCE_COMMIT "unavailable")
set(HUNDUN_SOURCE_TREE "unavailable")
set(HUNDUN_SOURCE_CLEAN "false")
set(HUNDUN_SOURCE_DIRTY_SUMMARY "unavailable")
include("${hundun_repository_root}/cmake/HundunPerformanceSourceMetadata.cmake")
find_package(Git QUIET)
if(Git_FOUND)
  # Candidate identity is configure-time material, so an ordinary incremental
  # build must reconfigure when either the checked-out Git authority or any
  # v0.4 executable input changes.  In particular, linked worktrees keep HEAD
  # symbolic while the referenced branch file moves on commit.
  foreach(HUNDUN_GIT_AUTHORITY IN ITEMS HEAD index packed-refs)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" rev-parse --git-path
              "${HUNDUN_GIT_AUTHORITY}"
      WORKING_DIRECTORY "${hundun_repository_root}"
      RESULT_VARIABLE HUNDUN_GIT_PATH_RESULT
      OUTPUT_VARIABLE HUNDUN_GIT_PATH
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
    if(HUNDUN_GIT_PATH_RESULT EQUAL 0 AND EXISTS "${HUNDUN_GIT_PATH}")
      set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
        "${HUNDUN_GIT_PATH}")
    endif()
  endforeach()
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" symbolic-ref -q HEAD
    WORKING_DIRECTORY "${hundun_repository_root}"
    RESULT_VARIABLE HUNDUN_GIT_SYMBOLIC_REF_RESULT
    OUTPUT_VARIABLE HUNDUN_GIT_SYMBOLIC_REF
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if(HUNDUN_GIT_SYMBOLIC_REF_RESULT EQUAL 0)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" rev-parse --git-path
              "${HUNDUN_GIT_SYMBOLIC_REF}"
      WORKING_DIRECTORY "${hundun_repository_root}"
      RESULT_VARIABLE HUNDUN_GIT_REF_PATH_RESULT
      OUTPUT_VARIABLE HUNDUN_GIT_REF_PATH
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
    if(HUNDUN_GIT_REF_PATH_RESULT EQUAL 0 AND EXISTS "${HUNDUN_GIT_REF_PATH}")
      set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
        "${HUNDUN_GIT_REF_PATH}")
    endif()
  endif()
  file(GLOB_RECURSE HUNDUN_V04_IDENTITY_INPUTS CONFIGURE_DEPENDS
    "${hundun_repository_root}/CMakeLists.txt"
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*"
    "${hundun_repository_root}/cmake/*"
    "${hundun_v04_yyjson_root}/*")
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    ${HUNDUN_V04_IDENTITY_INPUTS})
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --verify HEAD
    WORKING_DIRECTORY "${hundun_repository_root}"
    RESULT_VARIABLE HUNDUN_GIT_COMMIT_RESULT
    OUTPUT_VARIABLE HUNDUN_GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if(HUNDUN_GIT_COMMIT_RESULT EQUAL 0 AND
     HUNDUN_GIT_COMMIT MATCHES "^[0-9a-fA-F]+$")
    string(TOLOWER "${HUNDUN_GIT_COMMIT}" HUNDUN_SOURCE_COMMIT)
  endif()
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --verify HEAD^{tree}
    WORKING_DIRECTORY "${hundun_repository_root}"
    RESULT_VARIABLE HUNDUN_GIT_TREE_RESULT
    OUTPUT_VARIABLE HUNDUN_GIT_TREE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if(HUNDUN_GIT_TREE_RESULT EQUAL 0 AND
     HUNDUN_GIT_TREE MATCHES "^[0-9a-fA-F]+$")
    string(TOLOWER "${HUNDUN_GIT_TREE}" HUNDUN_SOURCE_TREE)
  endif()
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" status --porcelain=v1 --untracked-files=normal
            -- CMakeLists.txt cmake versions/v0.4/CMakeLists.txt
            versions/v0.4/include versions/v0.4/src third_party/yyjson
    WORKING_DIRECTORY "${hundun_repository_root}"
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
string(STRIP "${HUNDUN_PERFORMANCE_COMPILER_FLAGS}"
  HUNDUN_PERFORMANCE_COMPILER_FLAGS)
set(HUNDUN_PERFORMANCE_LINK_FLAGS
  "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${HUNDUN_BUILD_TYPE_UPPER}}")
string(STRIP "${HUNDUN_PERFORMANCE_LINK_FLAGS}"
  HUNDUN_PERFORMANCE_LINK_FLAGS)
string(SHA256 HUNDUN_SOURCE_HEAD_SHA256
  "hundun-git-head-v1:${HUNDUN_SOURCE_COMMIT}")
string(SHA256 HUNDUN_SOURCE_TREE_SHA256
  "hundun-git-tree-v1:${HUNDUN_SOURCE_TREE}")
set(HUNDUN_BUILD_MANIFEST_CANONICAL
  "schema=HUNDUN_V04_BUILD_MANIFEST_V1\n")
string(APPEND HUNDUN_BUILD_MANIFEST_CANONICAL
  "head=${HUNDUN_SOURCE_HEAD_SHA256}\n"
  "tree=${HUNDUN_SOURCE_TREE_SHA256}\n"
  "source_clean=${HUNDUN_SOURCE_CLEAN}\n"
  "dirty_summary=${HUNDUN_SOURCE_DIRTY_SUMMARY}\n"
  "compiler_id=${CMAKE_CXX_COMPILER_ID}\n"
  "compiler_version=${CMAKE_CXX_COMPILER_VERSION}\n"
  "compiler_flags=${HUNDUN_PERFORMANCE_COMPILER_FLAGS}\n"
  "link_flags=${HUNDUN_PERFORMANCE_LINK_FLAGS}\n"
  "build_type=${CMAKE_BUILD_TYPE}\n"
  "source_version=v0.4\n"
  "tests=${HUNDUN_BUILD_TESTS}\n"
  "asan=${HUNDUN_ENABLE_ASAN}\n"
  "ubsan=${HUNDUN_ENABLE_UBSAN}\n"
  "hypre=${HUNDUN_ENABLE_HYPRE}\n")
string(SHA256 HUNDUN_BUILD_MANIFEST_SHA256
  "${HUNDUN_BUILD_MANIFEST_CANONICAL}")
foreach(HUNDUN_METADATA_FIELD IN ITEMS
    HUNDUN_SOURCE_COMMIT
    HUNDUN_SOURCE_TREE
    HUNDUN_SOURCE_HEAD_SHA256
    HUNDUN_SOURCE_TREE_SHA256
    HUNDUN_SOURCE_DIRTY_SUMMARY
    HUNDUN_PERFORMANCE_COMPILER_FLAGS
    HUNDUN_PERFORMANCE_LINK_FLAGS
    HUNDUN_BUILD_MANIFEST_SHA256)
  hundun_escape_cpp_string(
    "${${HUNDUN_METADATA_FIELD}}" HUNDUN_METADATA_ESCAPED)
  set("${HUNDUN_METADATA_FIELD}_ESCAPED" "${HUNDUN_METADATA_ESCAPED}")
endforeach()
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/app_identity_build.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/app_identity_build.hpp"
  @ONLY)

add_library(hundun_v04_build_options INTERFACE)
target_compile_definitions(hundun_v04_build_options INTERFACE
  OMPI_SKIP_MPICXX=1
  MPICH_SKIP_MPICXX=1)
if(HUNDUN_ENABLE_ASAN OR HUNDUN_ENABLE_UBSAN)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    message(FATAL_ERROR
      "v0.4 sanitizers require a GNU or Clang C++ compiler; found "
      "${CMAKE_CXX_COMPILER_ID}")
  endif()

  set(hundun_v04_sanitizer_flags "")
  if(HUNDUN_ENABLE_ASAN)
    list(APPEND hundun_v04_sanitizer_flags -fsanitize=address)
  endif()
  if(HUNDUN_ENABLE_UBSAN)
    list(APPEND hundun_v04_sanitizer_flags -fsanitize=undefined)
  endif()
  target_compile_options(hundun_v04_build_options INTERFACE
    ${hundun_v04_sanitizer_flags})
  target_link_options(hundun_v04_build_options INTERFACE
    ${hundun_v04_sanitizer_flags})
endif()

set(hundun_v04_core_sources
  src/common_terminal_audit.c
  src/app_case.cpp
  src/app_driver.cpp
  src/app_identity.cpp
  src/bc_apply.cpp
  src/bc_compile.cpp
  src/bc_nscbc.cpp
  src/bc_thermo.cpp
  src/bc_time.cpp
  src/core_arena.cpp
  src/core_execution_graph.cpp
  src/core_field_registry.cpp
  src/core_product_freeze.cpp
  src/core_resource_contract.cpp
  src/core_revision.cpp
  src/core_status.cpp
  src/core_transaction.cpp
  src/mesh_cartesian.cpp
  src/mesh_decomposition.cpp
  src/mesh_focus.cpp
  src/mesh_ibm_plan.cpp
  src/mesh_ibm_reconstruction.cpp
  src/mesh_ibm_surface.cpp
  src/mesh_ibm_topology.cpp
  src/mesh_stl_scan.cpp
  src/parallel_cpu.cpp
  src/parallel_halo.cpp
  src/parallel_ibm_donor.cpp
  src/physics_contribution.cpp
  src/physics_derived.cpp
  src/physics_input.cpp
  src/physics_thermo.cpp
  src/physics_transport.cpp
  src/physics_vreman.cpp
  src/physics_wale.cpp
  src/physics_wall_function.cpp
  src/io_service_plan.cpp
  src/io_restart.cpp
  src/io_visit.cpp
  src/io_screen.cpp
  src/io_monitor.cpp
  src/io_evidence.cpp
  src/solver_cartesian.cpp
  src/solver_candidate_boundary.cpp
  src/solver_diffusion.cpp
  src/solver_enthalpy.cpp
  src/solver_equations.cpp
  src/solver_face_flux.cpp
  src/solver_hypre.cpp
  src/solver_ibm_force.cpp
  src/solver_ibm_equations.cpp
  src/solver_ibm_pressure.cpp
  src/solver_initialization_projection.cpp
  src/solver_krylov.cpp
  src/solver_linear_state.cpp
  src/solver_mg.cpp
  src/solver_mg_workspace.cpp
  src/solver_momentum.cpp
  src/solver_piso.cpp
  src/solver_pressure.cpp
  src/solver_pressure_energy.cpp
  src/solver_reduction.cpp
  src/solver_species.cpp
  src/solver_thermophysical_predictor.cpp
  src/solver_thermophysical_rates.cpp
  "${hundun_v04_yyjson_root}/yyjson.c")

function(hundun_v04_configure_core target_name)
  target_include_directories(${target_name} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include")
  target_include_directories(${target_name} PRIVATE
    "${hundun_v04_yyjson_root}"
    "${CMAKE_CURRENT_BINARY_DIR}/generated")
  target_link_libraries(${target_name} PUBLIC
    hundun_v04_build_options
    MPI::MPI_C)
  target_link_libraries(${target_name} PRIVATE
    Threads::Threads)
  if(HUNDUN_ENABLE_HYPRE)
    target_compile_definitions(${target_name} PRIVATE
      HUNDUN_V04_HAVE_HYPRE=1)
    target_link_libraries(${target_name} PRIVATE HYPRE::HYPRE)
  endif()
endfunction()

# The production core is independent of whether tests are configured.
add_library(hundun_v04_core STATIC ${hundun_v04_core_sources})
hundun_v04_configure_core(hundun_v04_core)

add_executable(hundun
  src/app_main.cpp)
target_link_libraries(hundun PRIVATE hundun_v04_core)
target_compile_definitions(hundun PRIVATE
  HUNDUN_V04_SOURCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}")

# Benchmark observer for segmented Re=3900 thin-domain runs.  It links the
# production core and has no test-access hooks, so adding the observer cannot
# change the numerical product executable or its compile definitions.
add_executable(v04_thin_domain_runner
  "${hundun_repository_root}/tools/v04_thin_domain_runner.cpp")
target_link_libraries(v04_thin_domain_runner PRIVATE hundun_v04_core)

if(HUNDUN_BUILD_TESTS)
  # Test-only hooks live in a separate core variant so they cannot alter the
  # production executable or its compile/link inputs.
  add_library(hundun_v04_test_core STATIC ${hundun_v04_core_sources})
  hundun_v04_configure_core(hundun_v04_test_core)
  target_compile_definitions(hundun_v04_test_core PRIVATE
    HUNDUN_V04_ENABLE_TEST_ACCESS=1)
  enable_testing()
  add_subdirectory(tests)
endif()
