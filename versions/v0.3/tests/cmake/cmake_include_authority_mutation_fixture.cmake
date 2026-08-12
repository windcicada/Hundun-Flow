# SPDX-License-Identifier: Apache-2.0

foreach(HUNDUN_REQUIRED IN ITEMS HUNDUN_SOURCE_ROOT HUNDUN_TEST_BINARY_ROOT)
  if(NOT DEFINED ${HUNDUN_REQUIRED} OR "${${HUNDUN_REQUIRED}}" STREQUAL "")
    message(FATAL_ERROR "${HUNDUN_REQUIRED} is required")
  endif()
endforeach()

get_filename_component(HUNDUN_SOURCE_ROOT "${HUNDUN_SOURCE_ROOT}" ABSOLUTE)
get_filename_component(HUNDUN_TEST_BINARY_ROOT
  "${HUNDUN_TEST_BINARY_ROOT}" ABSOLUTE)
set(HUNDUN_MUTATION_ROOT
  "${HUNDUN_TEST_BINARY_ROOT}/cmake-include-authority-mutation-fixture")
cmake_path(GET HUNDUN_MUTATION_ROOT PARENT_PATH HUNDUN_MUTATION_PARENT)
cmake_path(GET HUNDUN_MUTATION_ROOT FILENAME HUNDUN_MUTATION_NAME)
if(NOT HUNDUN_MUTATION_PARENT STREQUAL HUNDUN_TEST_BINARY_ROOT OR
   NOT HUNDUN_MUTATION_NAME STREQUAL
     "cmake-include-authority-mutation-fixture")
  message(FATAL_ERROR "mutation fixture escaped its binary root")
endif()
if(IS_SYMLINK "${HUNDUN_MUTATION_ROOT}")
  message(FATAL_ERROR "mutation fixture root cannot be a symlink")
endif()
file(REMOVE_RECURSE "${HUNDUN_MUTATION_ROOT}")
file(MAKE_DIRECTORY "${HUNDUN_MUTATION_ROOT}")

set(HUNDUN_AUTHORITY_FIXTURE
  "${HUNDUN_SOURCE_ROOT}/tests/cmake/cmake_include_authority_fixture.cmake")
if(NOT EXISTS "${HUNDUN_AUTHORITY_FIXTURE}")
  message(FATAL_ERROR "include-authority fixture is missing")
endif()

set(HUNDUN_PROJECT_TEMPLATE [=[
cmake_minimum_required(VERSION 3.25)
project(hundun_authority_mutation LANGUAGES CXX)

@HUNDUN_MUTATION_BEFORE@

add_library(hundun_options INTERFACE)
target_include_directories(hundun_options INTERFACE
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>")

set(HUNDUN_PRODUCT_LIBRARY_TARGETS
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
  hundun_material_diagnostics
  hundun_checkpoint_diagnostics
  hundun_session_diagnostics
  hundun_io
  hundun_solver
  hundun_sdk)
foreach(HUNDUN_TARGET IN LISTS HUNDUN_PRODUCT_LIBRARY_TARGETS)
  add_library("${HUNDUN_TARGET}" STATIC dummy.cpp)
  target_include_directories("${HUNDUN_TARGET}"
    PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    PRIVATE "${PROJECT_SOURCE_DIR}/src")
endforeach()

add_executable(hundun dummy.cpp)
target_include_directories(hundun PRIVATE
  "${PROJECT_SOURCE_DIR}/src"
  "${PROJECT_BINARY_DIR}/src/generated")

@HUNDUN_MUTATION_AFTER@

include("@HUNDUN_AUTHORITY_FIXTURE@")
hundun_validate_configured_product_include_authority()
]=])

foreach(HUNDUN_CASE IN ITEMS
    control
    repository_root
    tests_root
    directory_root
    target_property_tests
    included_module_root)
  set(HUNDUN_CASE_ROOT "${HUNDUN_MUTATION_ROOT}/${HUNDUN_CASE}")
  set(HUNDUN_CASE_SOURCE "${HUNDUN_CASE_ROOT}/source")
  set(HUNDUN_CASE_BUILD "${HUNDUN_CASE_ROOT}/build")
  file(MAKE_DIRECTORY
    "${HUNDUN_CASE_SOURCE}/include"
    "${HUNDUN_CASE_SOURCE}/src"
    "${HUNDUN_CASE_SOURCE}/cmake")
  file(WRITE "${HUNDUN_CASE_SOURCE}/dummy.cpp" "int probe() { return 0; }\n")

  set(HUNDUN_MUTATION_BEFORE "")
  set(HUNDUN_MUTATION_AFTER "")
  set(HUNDUN_EXPECTED "")
  if(HUNDUN_CASE STREQUAL "repository_root")
    set(HUNDUN_MUTATION_AFTER [=[
target_include_directories(hundun_linear PRIVATE "${PROJECT_SOURCE_DIR}")
]=])
    set(HUNDUN_EXPECTED "repository-root include authority")
  elseif(HUNDUN_CASE STREQUAL "tests_root")
    set(HUNDUN_MUTATION_AFTER [=[
target_include_directories(hundun_flow PRIVATE
  "${PROJECT_SOURCE_DIR}/tests/support")
]=])
    set(HUNDUN_EXPECTED "tests include authority")
  elseif(HUNDUN_CASE STREQUAL "directory_root")
    set(HUNDUN_MUTATION_BEFORE [=[
include_directories("${PROJECT_SOURCE_DIR}")
]=])
    set(HUNDUN_EXPECTED "repository-root include authority")
  elseif(HUNDUN_CASE STREQUAL "target_property_tests")
    set(HUNDUN_MUTATION_AFTER [=[
set_property(TARGET hundun_flow APPEND PROPERTY INCLUDE_DIRECTORIES
  "${PROJECT_SOURCE_DIR}/tests/support")
]=])
    set(HUNDUN_EXPECTED "tests include authority")
  elseif(HUNDUN_CASE STREQUAL "included_module_root")
    file(WRITE "${HUNDUN_CASE_SOURCE}/cmake/injected.cmake" [=[
add_library(hundun_injected_helper INTERFACE)
target_include_directories(hundun_injected_helper INTERFACE
  "${PROJECT_SOURCE_DIR}")
target_link_libraries(hundun_flow PUBLIC hundun_injected_helper)
]=])
    set(HUNDUN_MUTATION_AFTER [=[
include("${PROJECT_SOURCE_DIR}/cmake/injected.cmake")
]=])
    set(HUNDUN_EXPECTED "repository-root include authority")
  endif()

  string(CONFIGURE "${HUNDUN_PROJECT_TEMPLATE}" HUNDUN_PROJECT_CMAKE @ONLY)
  file(WRITE "${HUNDUN_CASE_SOURCE}/CMakeLists.txt"
    "${HUNDUN_PROJECT_CMAKE}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}"
      -S "${HUNDUN_CASE_SOURCE}"
      -B "${HUNDUN_CASE_BUILD}"
    RESULT_VARIABLE HUNDUN_RESULT
    OUTPUT_VARIABLE HUNDUN_STDOUT
    ERROR_VARIABLE HUNDUN_STDERR)
  set(HUNDUN_OUTPUT "${HUNDUN_STDOUT}${HUNDUN_STDERR}")

  if(HUNDUN_CASE STREQUAL "control")
    if(NOT HUNDUN_RESULT EQUAL 0)
      message(FATAL_ERROR
        "include-authority control project failed:\n${HUNDUN_OUTPUT}")
    endif()
  else()
    if(HUNDUN_RESULT EQUAL 0)
      message(FATAL_ERROR
        "include-authority fixture accepted '${HUNDUN_CASE}' mutation")
    endif()
    if(NOT HUNDUN_OUTPUT MATCHES "${HUNDUN_EXPECTED}")
      message(FATAL_ERROR
        "include-authority fixture rejected '${HUNDUN_CASE}' for the wrong reason:\n${HUNDUN_OUTPUT}")
    endif()
  endif()
endforeach()

file(REMOVE_RECURSE "${HUNDUN_MUTATION_ROOT}")
