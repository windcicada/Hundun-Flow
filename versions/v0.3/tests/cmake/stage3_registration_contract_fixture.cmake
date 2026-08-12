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
  "${HUNDUN_TEST_BINARY_ROOT}/stage3-registration-contract-fixture")
cmake_path(GET HUNDUN_MUTATION_ROOT PARENT_PATH HUNDUN_MUTATION_PARENT)
cmake_path(GET HUNDUN_MUTATION_ROOT FILENAME HUNDUN_MUTATION_NAME)
if(NOT HUNDUN_MUTATION_PARENT STREQUAL HUNDUN_TEST_BINARY_ROOT OR
   NOT HUNDUN_MUTATION_NAME STREQUAL
     "stage3-registration-contract-fixture")
  message(FATAL_ERROR "stage3 registration fixture escaped its binary root")
endif()
if(IS_SYMLINK "${HUNDUN_MUTATION_ROOT}")
  message(FATAL_ERROR "stage3 registration fixture root cannot be a symlink")
endif()

file(REMOVE_RECURSE "${HUNDUN_MUTATION_ROOT}")
file(MAKE_DIRECTORY "${HUNDUN_MUTATION_ROOT}/tests/cmake")
file(COPY "${HUNDUN_SOURCE_ROOT}/tests/CMakeLists.txt"
  DESTINATION "${HUNDUN_MUTATION_ROOT}/tests")

set(HUNDUN_STAGE3_REGISTRATION_FRAGMENTS
  stage3_science_registration.cmake
  stage3_checkpoint_registration.cmake
  stage3_diagnostics_registration.cmake
  stage3_framework_registration.cmake
  stage3_acceptance_registration.cmake)
foreach(HUNDUN_STAGE3_FRAGMENT IN LISTS
    HUNDUN_STAGE3_REGISTRATION_FRAGMENTS)
  set(HUNDUN_STAGE3_SOURCE_FRAGMENT
    "${HUNDUN_SOURCE_ROOT}/tests/cmake/${HUNDUN_STAGE3_FRAGMENT}")
  if(NOT EXISTS "${HUNDUN_STAGE3_SOURCE_FRAGMENT}")
    message(FATAL_ERROR
      "stage3 registration fixture source is missing: ${HUNDUN_STAGE3_FRAGMENT}")
  endif()
  file(COPY "${HUNDUN_STAGE3_SOURCE_FRAGMENT}"
    DESTINATION "${HUNDUN_MUTATION_ROOT}/tests/cmake")
endforeach()

set(HUNDUN_MUTATED_CMAKE
  "${HUNDUN_MUTATION_ROOT}/tests/CMakeLists.txt")
file(READ "${HUNDUN_MUTATED_CMAKE}" HUNDUN_MUTATED_CMAKE_TEXT)
set(HUNDUN_REMOVED_INCLUDE
  "include(\"\${PROJECT_SOURCE_DIR}/tests/cmake/stage3_acceptance_registration.cmake\")")
string(FIND "${HUNDUN_MUTATED_CMAKE_TEXT}"
  "${HUNDUN_REMOVED_INCLUDE}" HUNDUN_REMOVED_INCLUDE_POSITION)
if(HUNDUN_REMOVED_INCLUDE_POSITION EQUAL -1)
  message(FATAL_ERROR
    "stage3 registration fixture could not find its mutation target")
endif()
string(REPLACE "${HUNDUN_REMOVED_INCLUDE}" ""
  HUNDUN_MUTATED_CMAKE_TEXT "${HUNDUN_MUTATED_CMAKE_TEXT}")
file(WRITE "${HUNDUN_MUTATED_CMAKE}" "${HUNDUN_MUTATED_CMAKE_TEXT}")

execute_process(
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_STAGE3_REGISTRATION_ROOT=${HUNDUN_MUTATION_ROOT}"
    -P "${HUNDUN_SOURCE_ROOT}/tests/cmake/stage3_registration_contract.cmake"
  RESULT_VARIABLE HUNDUN_MUTATION_RESULT
  OUTPUT_VARIABLE HUNDUN_MUTATION_STDOUT
  ERROR_VARIABLE HUNDUN_MUTATION_STDERR)
set(HUNDUN_MUTATION_OUTPUT
  "${HUNDUN_MUTATION_STDOUT}${HUNDUN_MUTATION_STDERR}")

if(HUNDUN_MUTATION_RESULT EQUAL 0)
  message(FATAL_ERROR
    "stage3 registration contract accepted a missing include mutation")
endif()
if(NOT HUNDUN_MUTATION_OUTPUT MATCHES
    "missing stage3 registration include: stage3_acceptance_registration\\.cmake")
  message(FATAL_ERROR
    "stage3 registration contract rejected the mutation for the wrong reason:\n${HUNDUN_MUTATION_OUTPUT}")
endif()

file(REMOVE_RECURSE "${HUNDUN_MUTATION_ROOT}")
