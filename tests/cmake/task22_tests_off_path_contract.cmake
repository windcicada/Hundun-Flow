# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SETUP_SCRIPT OR
   NOT DEFINED HUNDUN_TEST_BINARY_ROOT OR
   NOT DEFINED HUNDUN_GENERATOR)
  message(FATAL_ERROR "Task22 path-contract inputs are incomplete")
endif()

set(contract_root
    "${HUNDUN_TEST_BINARY_ROOT}/task22-tests-off-path-contract")
file(MAKE_DIRECTORY "${contract_root}")

function(task22_write_outer_cache outer)
  file(MAKE_DIRECTORY "${outer}")
  file(WRITE "${outer}/CMakeCache.txt"
       "CMAKE_GENERATOR:INTERNAL=${HUNDUN_GENERATOR}\n")
endfunction()

function(task22_require_unsafe_rejection name source outer nested archive)
  set(case_root "${contract_root}/${name}")
  file(REMOVE_RECURSE "${case_root}")
  file(MAKE_DIRECTORY "${source}")
  task22_write_outer_cache("${outer}")
  file(MAKE_DIRECTORY "${nested}")
  set(sentinel "${nested}/sentinel")
  file(WRITE "${sentinel}" "Task22 path sentinel: ${name}\n")
  file(SHA256 "${sentinel}" expected_sha)

  execute_process(
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_ROOT=${source}"
            "-DHUNDUN_OUTER_BINARY_DIR=${outer}"
            "-DHUNDUN_NESTED_BINARY_DIR=${nested}"
            "-DHUNDUN_ARCHIVE_DIR=${archive}"
            "-DHUNDUN_NESTED_CONFIG=Release"
            -P "${HUNDUN_SETUP_SCRIPT}"
    RESULT_VARIABLE setup_result
    OUTPUT_QUIET
    ERROR_QUIET)
  if(setup_result EQUAL 0)
    message(FATAL_ERROR
      "Task22 unsafe path case '${name}' unexpectedly succeeded")
  endif()
  if(NOT EXISTS "${sentinel}")
    message(FATAL_ERROR
      "Task22 unsafe path case '${name}' removed its sentinel")
  endif()
  file(SHA256 "${sentinel}" observed_sha)
  if(NOT observed_sha STREQUAL expected_sha)
    message(FATAL_ERROR
      "Task22 unsafe path case '${name}' changed its sentinel")
  endif()
endfunction()

set(sibling_root "${contract_root}/nested-sibling")
task22_require_unsafe_rejection(
  nested-sibling
  "${sibling_root}/source"
  "${sibling_root}/outer"
  "${sibling_root}/victim"
  "${sibling_root}/victim/archives")

set(equal_outer_root "${contract_root}/nested-equals-outer")
task22_require_unsafe_rejection(
  nested-equals-outer
  "${equal_outer_root}/source"
  "${equal_outer_root}/outer"
  "${equal_outer_root}/outer"
  "${equal_outer_root}/outer/archives")

set(equal_source_root "${contract_root}/nested-equals-source")
task22_require_unsafe_rejection(
  nested-equals-source
  "${equal_source_root}/source"
  "${equal_source_root}/outer"
  "${equal_source_root}/source"
  "${equal_source_root}/source/archives")

set(wrong_suffix_root "${contract_root}/wrong-nested-suffix")
task22_require_unsafe_rejection(
  wrong-nested-suffix
  "${wrong_suffix_root}/source"
  "${wrong_suffix_root}/outer"
  "${wrong_suffix_root}/outer/not-task22-tests-off"
  "${wrong_suffix_root}/outer/not-task22-tests-off/archives")

set(archive_outside_root "${contract_root}/archive-outside")
task22_require_unsafe_rejection(
  archive-outside
  "${archive_outside_root}/source"
  "${archive_outside_root}/outer"
  "${archive_outside_root}/outer/task22-tests-off"
  "${archive_outside_root}/archive-victim")

set(archive_equal_root "${contract_root}/archive-equals-nested")
task22_require_unsafe_rejection(
  archive-equals-nested
  "${archive_equal_root}/source"
  "${archive_equal_root}/outer"
  "${archive_equal_root}/outer/task22-tests-off"
  "${archive_equal_root}/outer/task22-tests-off")

set(symlink_root "${contract_root}/nested-symlink")
file(REMOVE_RECURSE "${symlink_root}")
file(MAKE_DIRECTORY "${symlink_root}/source")
task22_write_outer_cache("${symlink_root}/outer")
file(MAKE_DIRECTORY "${symlink_root}/outside-victim")
set(symlink_sentinel "${symlink_root}/outside-victim/sentinel")
file(WRITE "${symlink_sentinel}" "Task22 symlink sentinel\n")
file(SHA256 "${symlink_sentinel}" symlink_expected_sha)
file(CREATE_LINK
     "${symlink_root}/outside-victim"
     "${symlink_root}/outer/task22-tests-off"
     SYMBOLIC
     RESULT symlink_result)
if(NOT symlink_result STREQUAL "0")
  message(FATAL_ERROR "Task22 path contract could not create its symlink")
endif()
execute_process(
  COMMAND "${CMAKE_COMMAND}"
          "-DHUNDUN_SOURCE_ROOT=${symlink_root}/source"
          "-DHUNDUN_OUTER_BINARY_DIR=${symlink_root}/outer"
          "-DHUNDUN_NESTED_BINARY_DIR=${symlink_root}/outer/task22-tests-off"
          "-DHUNDUN_ARCHIVE_DIR=${symlink_root}/outer/task22-tests-off/archives"
          "-DHUNDUN_NESTED_CONFIG=Release"
          -P "${HUNDUN_SETUP_SCRIPT}"
  RESULT_VARIABLE symlink_setup_result
  OUTPUT_QUIET
  ERROR_QUIET)
if(symlink_setup_result EQUAL 0 OR
   NOT IS_SYMLINK "${symlink_root}/outer/task22-tests-off" OR
   NOT EXISTS "${symlink_sentinel}")
  message(FATAL_ERROR
    "Task22 symlink path case was not rejected without mutation")
endif()
file(SHA256 "${symlink_sentinel}" symlink_observed_sha)
if(NOT symlink_observed_sha STREQUAL symlink_expected_sha)
  message(FATAL_ERROR "Task22 symlink path case changed its sentinel")
endif()

set(root_source "${contract_root}/filesystem-root/source")
file(MAKE_DIRECTORY "${root_source}")
execute_process(
  COMMAND "${CMAKE_COMMAND}"
          "-DHUNDUN_SOURCE_ROOT=${root_source}"
          "-DHUNDUN_OUTER_BINARY_DIR=/"
          "-DHUNDUN_NESTED_BINARY_DIR=/task22-tests-off"
          "-DHUNDUN_ARCHIVE_DIR=/task22-tests-off/archives"
          "-DHUNDUN_NESTED_CONFIG=Release"
          "-DHUNDUN_TASK22_PATH_CONTRACT_VALIDATE_ONLY=ON"
          -P "${HUNDUN_SETUP_SCRIPT}"
  RESULT_VARIABLE root_setup_result
  OUTPUT_QUIET
  ERROR_QUIET)
if(root_setup_result EQUAL 0)
  message(FATAL_ERROR "Task22 filesystem-root outer path was not rejected")
endif()
