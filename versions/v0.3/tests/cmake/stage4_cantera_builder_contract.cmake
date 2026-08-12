# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

if(NOT DEFINED HUNDUN_SOURCE_DIR)
  message(FATAL_ERROR "HUNDUN_SOURCE_DIR is required")
endif()
get_filename_component(HUNDUN_SOURCE_DIR "${HUNDUN_SOURCE_DIR}" ABSOLUTE)
set(cantera_module "${HUNDUN_SOURCE_DIR}/cmake/HundunCanteraPackage.cmake")
if(NOT EXISTS "${cantera_module}")
  message(FATAL_ERROR "Stage 4 Cantera package module is missing")
endif()
include("${cantera_module}")

if(DEFINED HUNDUN_CANTERA_CONSUMER_MUTATION)
  if(HUNDUN_CANTERA_CONSUMER_MUTATION STREQUAL "fetchcontent")
    hundun_validate_cantera_consumer_text(
      "fetchcontent" "FetchContent_Declare(cantera GIT_REPOSITORY local)")
  elseif(HUNDUN_CANTERA_CONSUMER_MUTATION STREQUAL "network_url")
    hundun_validate_cantera_consumer_text(
      "network_url" "set(source https://example.invalid/cantera.tar.gz)")
  elseif(HUNDUN_CANTERA_CONSUMER_MUTATION STREQUAL "builder_invocation")
    hundun_validate_cantera_consumer_text(
      "builder_invocation" "execute_process(COMMAND build_bundled_cantera_linux_cpu.sh)")
  else()
    message(FATAL_ERROR "unknown Cantera consumer mutation")
  endif()
  message(FATAL_ERROR "Cantera consumer mutation was unexpectedly accepted")
endif()

foreach(mutation IN ITEMS fetchcontent network_url builder_invocation)
  execute_process(
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_DIR=${HUNDUN_SOURCE_DIR}"
            "-DHUNDUN_CANTERA_CONSUMER_MUTATION=${mutation}"
            -P "${CMAKE_CURRENT_LIST_FILE}"
    RESULT_VARIABLE mutation_result
    OUTPUT_VARIABLE mutation_stdout
    ERROR_VARIABLE mutation_stderr)
  if(mutation_result EQUAL 0 OR
     NOT "${mutation_stdout}${mutation_stderr}" MATCHES
         "Cantera consumer policy rejected ${mutation}")
    message(FATAL_ERROR "Cantera ${mutation} mutation was not rejected")
  endif()
endforeach()

set(contract_root "${HUNDUN_SOURCE_DIR}/build/stage4-cantera-contract")
set(empty_root "${contract_root}/empty-package")
file(REMOVE_RECURSE "${empty_root}")
file(MAKE_DIRECTORY "${empty_root}")
set(fixture_source "${contract_root}/fixture-source")
set(fixture_build "${contract_root}/fixture-build")
file(REMOVE_RECURSE "${fixture_source}" "${fixture_build}")
file(MAKE_DIRECTORY "${fixture_source}")
file(WRITE "${fixture_source}/CMakeLists.txt"
  "cmake_minimum_required(VERSION 3.21)\n"
  "project(stage4_cantera_fixture LANGUAGES CXX)\n"
  "include(\"${cantera_module}\")\n"
  "hundun_configure_cantera_package(\"${empty_root}\")\n")
execute_process(
  COMMAND "${CMAKE_COMMAND}" -S "${fixture_source}" -B "${fixture_build}"
  RESULT_VARIABLE empty_result
  OUTPUT_VARIABLE empty_stdout
  ERROR_VARIABLE empty_stderr)
if(empty_result EQUAL 0 OR
   NOT "${empty_stdout}${empty_stderr}" MATCHES
       "verified Cantera package is incomplete")
  message(FATAL_ERROR "empty Cantera package root did not produce precise error")
endif()

message(STATUS "Stage 4 Cantera builder/consumer contract passed")
