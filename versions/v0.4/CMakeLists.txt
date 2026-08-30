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
    "${hundun_v04_yyjson_root}")
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
