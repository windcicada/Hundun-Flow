# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_STAGE3_REGISTRATION_ROOT OR
   HUNDUN_STAGE3_REGISTRATION_ROOT STREQUAL "")
  message(FATAL_ERROR "HUNDUN_STAGE3_REGISTRATION_ROOT is required")
endif()

get_filename_component(HUNDUN_STAGE3_REGISTRATION_ROOT
  "${HUNDUN_STAGE3_REGISTRATION_ROOT}" ABSOLUTE)
set(HUNDUN_STAGE3_TEST_CMAKE
  "${HUNDUN_STAGE3_REGISTRATION_ROOT}/tests/CMakeLists.txt")

if(NOT EXISTS "${HUNDUN_STAGE3_TEST_CMAKE}")
  message(FATAL_ERROR
    "stage3 registration root has no tests/CMakeLists.txt")
endif()

file(READ "${HUNDUN_STAGE3_TEST_CMAKE}" HUNDUN_STAGE3_TEST_CMAKE_TEXT)

function(_hundun_stage3_count_exact haystack needle output)
  set(remaining "${haystack}")
  string(LENGTH "${needle}" needle_length)
  set(count 0)
  set(searching 1)
  while(searching)
    string(FIND "${remaining}" "${needle}" position)
    if(position EQUAL -1)
      set(searching 0)
      continue()
    endif()
    math(EXPR next "${position} + ${needle_length}")
    string(SUBSTRING "${remaining}" "${next}" -1 remaining)
    math(EXPR count "${count} + 1")
  endwhile()
  set("${output}" "${count}" PARENT_SCOPE)
endfunction()

set(HUNDUN_STAGE3_REGISTRATION_FRAGMENTS
  stage3_science_registration.cmake
  stage3_checkpoint_registration.cmake
  stage3_diagnostics_registration.cmake
  stage3_framework_registration.cmake
  stage3_acceptance_registration.cmake)

foreach(HUNDUN_STAGE3_FRAGMENT IN LISTS
    HUNDUN_STAGE3_REGISTRATION_FRAGMENTS)
  set(HUNDUN_STAGE3_INCLUDE
    "include(\"\${PROJECT_SOURCE_DIR}/tests/cmake/${HUNDUN_STAGE3_FRAGMENT}\")")
  _hundun_stage3_count_exact(
    "${HUNDUN_STAGE3_TEST_CMAKE_TEXT}"
    "${HUNDUN_STAGE3_INCLUDE}"
    HUNDUN_STAGE3_INCLUDE_COUNT)
  if(HUNDUN_STAGE3_INCLUDE_COUNT EQUAL 0)
    message(FATAL_ERROR
      "missing stage3 registration include: ${HUNDUN_STAGE3_FRAGMENT}")
  endif()
  if(NOT HUNDUN_STAGE3_INCLUDE_COUNT EQUAL 1)
    message(FATAL_ERROR
      "duplicate stage3 registration include: ${HUNDUN_STAGE3_FRAGMENT}")
  endif()

  set(HUNDUN_STAGE3_FRAGMENT_PATH
    "${HUNDUN_STAGE3_REGISTRATION_ROOT}/tests/cmake/${HUNDUN_STAGE3_FRAGMENT}")
  if(NOT EXISTS "${HUNDUN_STAGE3_FRAGMENT_PATH}")
    message(FATAL_ERROR
      "missing stage3 registration fragment: ${HUNDUN_STAGE3_FRAGMENT}")
  endif()

  file(READ "${HUNDUN_STAGE3_FRAGMENT_PATH}" HUNDUN_STAGE3_FRAGMENT_TEXT)
  if(NOT HUNDUN_STAGE3_FRAGMENT_TEXT MATCHES
      "^# SPDX-License-Identifier: Apache-2\\.0([\r\n])")
    message(FATAL_ERROR
      "stage3 registration fragment lacks SPDX: ${HUNDUN_STAGE3_FRAGMENT}")
  endif()
  _hundun_stage3_count_exact(
    "${HUNDUN_STAGE3_FRAGMENT_TEXT}"
    "include_guard(GLOBAL)"
    HUNDUN_STAGE3_GUARD_COUNT)
  if(NOT HUNDUN_STAGE3_GUARD_COUNT EQUAL 1)
    message(FATAL_ERROR
      "stage3 registration fragment lacks one include guard: ${HUNDUN_STAGE3_FRAGMENT}")
  endif()
endforeach()
