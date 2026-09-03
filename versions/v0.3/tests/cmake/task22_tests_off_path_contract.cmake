# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

if(NOT DEFINED HUNDUN_SETUP_SCRIPT OR
   NOT DEFINED HUNDUN_TEST_BINARY_ROOT OR
   NOT DEFINED HUNDUN_GENERATOR)
  message(FATAL_ERROR "Task22 path-contract inputs are incomplete")
endif()
if("${HUNDUN_TEST_BINARY_ROOT}" STREQUAL "")
  message(FATAL_ERROR "Task22 path-contract binary root is empty")
endif()

set(test_binary_root_input "${HUNDUN_TEST_BINARY_ROOT}")
cmake_path(ABSOLUTE_PATH test_binary_root_input NORMALIZE
           OUTPUT_VARIABLE test_binary_root_absolute)
if(NOT IS_DIRECTORY "${test_binary_root_absolute}")
  message(FATAL_ERROR
    "Task22 path-contract binary root is not an existing directory")
endif()
file(REAL_PATH "${test_binary_root_absolute}" test_binary_root)
cmake_path(GET test_binary_root ROOT_PATH test_binary_filesystem_root)
if(test_binary_root STREQUAL test_binary_filesystem_root)
  message(FATAL_ERROR
    "Task22 path-contract binary root cannot be a filesystem root")
endif()

set(contract_root "")
set(contract_nonce "")
foreach(collision_attempt RANGE 1 16)
  string(RANDOM LENGTH 32 ALPHABET 0123456789abcdef candidate_nonce)
  set(candidate_root
      "${test_binary_root}/task22-tests-off-path-contract-${candidate_nonce}")
  if(NOT EXISTS "${candidate_root}" AND NOT IS_SYMLINK "${candidate_root}")
    set(contract_nonce "${candidate_nonce}")
    set(contract_root "${candidate_root}")
    break()
  endif()
endforeach()
if(contract_root STREQUAL "")
  message(FATAL_ERROR
    "Task22 path-contract could not select a fresh owned root")
endif()

set(contract_root_name
    "task22-tests-off-path-contract-${contract_nonce}")
set(contract_marker "${contract_root}/.task22-path-contract-owner")
file(MAKE_DIRECTORY "${contract_root}")
file(WRITE "${contract_marker}" "${contract_nonce}\n")

function(task22_require_owned_root)
  if(NOT IS_DIRECTORY "${contract_root}" OR IS_SYMLINK "${contract_root}")
    message(FATAL_ERROR
      "Task22 path-contract owned root is missing or is a symlink")
  endif()
  file(REAL_PATH "${contract_root}" observed_contract_root)
  cmake_path(GET observed_contract_root PARENT_PATH observed_contract_parent)
  cmake_path(GET observed_contract_root FILENAME observed_contract_name)
  if(NOT observed_contract_parent STREQUAL test_binary_root OR
     NOT observed_contract_name STREQUAL contract_root_name)
    message(FATAL_ERROR
      "Task22 path-contract owned root escaped its validated binary root")
  endif()
  string(FIND "${observed_contract_name}"
         "task22-tests-off-path-contract-" contract_prefix_position)
  if(NOT contract_prefix_position EQUAL 0)
    message(FATAL_ERROR
      "Task22 path-contract owned root has an unexpected name")
  endif()
  if(NOT EXISTS "${contract_marker}" OR IS_SYMLINK "${contract_marker}")
    message(FATAL_ERROR
      "Task22 path-contract ownership marker is missing or is a symlink")
  endif()
  file(READ "${contract_marker}" observed_contract_marker)
  if(NOT observed_contract_marker STREQUAL "${contract_nonce}\n")
    message(FATAL_ERROR
      "Task22 path-contract ownership marker does not match this invocation")
  endif()
endfunction()

function(task22_prepare_case output_variable name)
  task22_require_owned_root()
  set(case_root "${contract_root}/${name}")
  if(EXISTS "${case_root}" OR IS_SYMLINK "${case_root}")
    message(FATAL_ERROR
      "Task22 path-contract case '${name}' was not freshly created")
  endif()
  file(MAKE_DIRECTORY "${case_root}")
  set("${output_variable}" "${case_root}" PARENT_SCOPE)
endfunction()

function(task22_write_outer_cache outer)
  file(MAKE_DIRECTORY "${outer}")
  file(WRITE "${outer}/CMakeCache.txt"
       "CMAKE_GENERATOR:INTERNAL=${HUNDUN_GENERATOR}\n")
endfunction()

function(task22_rejection_matches
         result_variable stdout_variable stderr_variable expected output_variable)
  set(combined_output
      "${${stdout_variable}}\n${${stderr_variable}}")
  string(REGEX REPLACE "[ \t\r\n]+" " " normalized_output
         "${combined_output}")
  string(REGEX REPLACE "[ \t\r\n]+" " " normalized_expected
         "${expected}")
  string(FIND "${normalized_output}" "${normalized_expected}"
         expected_position)
  if(NOT "${${result_variable}}" STREQUAL "0" AND
     NOT expected_position EQUAL -1)
    set("${output_variable}" TRUE PARENT_SCOPE)
  else()
    set("${output_variable}" FALSE PARENT_SCOPE)
  endif()
endfunction()

function(task22_require_unsafe_rejection
         name source outer nested archive expected_message)
  task22_prepare_case(case_root "${name}")
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
    OUTPUT_VARIABLE setup_stdout
    ERROR_VARIABLE setup_stderr)
  task22_rejection_matches(
    setup_result setup_stdout setup_stderr "${expected_message}"
    rejection_matches)
  if(NOT rejection_matches)
    message(FATAL_ERROR
      "Task22 unsafe path case '${name}' did not reach its expected guard "
      "'${expected_message}':\n${setup_stdout}${setup_stderr}")
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

task22_require_unsafe_rejection(
  nested-sibling
  "${contract_root}/nested-sibling/source"
  "${contract_root}/nested-sibling/outer"
  "${contract_root}/nested-sibling/victim"
  "${contract_root}/nested-sibling/victim/archives"
  "Task22 tests-off nested binary directory must be the dedicated outer/task22-tests-off path")

task22_require_unsafe_rejection(
  nested-equals-outer
  "${contract_root}/nested-equals-outer/source"
  "${contract_root}/nested-equals-outer/outer"
  "${contract_root}/nested-equals-outer/outer"
  "${contract_root}/nested-equals-outer/outer/archives"
  "Task22 tests-off nested binary directory must be the dedicated outer/task22-tests-off path")

task22_require_unsafe_rejection(
  nested-equals-source
  "${contract_root}/nested-equals-source/source"
  "${contract_root}/nested-equals-source/outer"
  "${contract_root}/nested-equals-source/source"
  "${contract_root}/nested-equals-source/source/archives"
  "Task22 tests-off nested binary directory must be the dedicated outer/task22-tests-off path")

task22_require_unsafe_rejection(
  wrong-nested-suffix
  "${contract_root}/wrong-nested-suffix/source"
  "${contract_root}/wrong-nested-suffix/outer"
  "${contract_root}/wrong-nested-suffix/outer/not-task22-tests-off"
  "${contract_root}/wrong-nested-suffix/outer/not-task22-tests-off/archives"
  "Task22 tests-off nested binary directory must be the dedicated outer/task22-tests-off path")

task22_require_unsafe_rejection(
  archive-outside
  "${contract_root}/archive-outside/source"
  "${contract_root}/archive-outside/outer"
  "${contract_root}/archive-outside/outer/task22-tests-off"
  "${contract_root}/archive-outside/archive-victim"
  "Task22 tests-off archive directory must be nested/archives")

task22_require_unsafe_rejection(
  archive-equals-nested
  "${contract_root}/archive-equals-nested/source"
  "${contract_root}/archive-equals-nested/outer"
  "${contract_root}/archive-equals-nested/outer/task22-tests-off"
  "${contract_root}/archive-equals-nested/outer/task22-tests-off"
  "Task22 tests-off archive directory must be nested/archives")

task22_prepare_case(symlink_root nested-symlink)
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
  OUTPUT_VARIABLE symlink_setup_stdout
  ERROR_VARIABLE symlink_setup_stderr)
task22_rejection_matches(
  symlink_setup_result symlink_setup_stdout symlink_setup_stderr
  "Task22 tests-off nested binary directory cannot be a symlink"
  symlink_rejection_matches)
if(NOT symlink_rejection_matches OR
   NOT IS_SYMLINK "${symlink_root}/outer/task22-tests-off" OR
   NOT EXISTS "${symlink_sentinel}")
  message(FATAL_ERROR
    "Task22 symlink path case was not rejected by its guard without mutation:"
    "\n${symlink_setup_stdout}${symlink_setup_stderr}")
endif()
file(SHA256 "${symlink_sentinel}" symlink_observed_sha)
if(NOT symlink_observed_sha STREQUAL symlink_expected_sha)
  message(FATAL_ERROR "Task22 symlink path case changed its sentinel")
endif()

task22_prepare_case(filesystem_root_case filesystem-root)
set(root_source "${filesystem_root_case}/source")
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
  OUTPUT_VARIABLE root_setup_stdout
  ERROR_VARIABLE root_setup_stderr)
task22_rejection_matches(
  root_setup_result root_setup_stdout root_setup_stderr
  "Task22 tests-off outer binary directory cannot be a filesystem root"
  root_rejection_matches)
if(NOT root_rejection_matches)
  message(FATAL_ERROR
    "Task22 filesystem-root case did not reach its exact guard:"
    "\n${root_setup_stdout}${root_setup_stderr}")
endif()

task22_prepare_case(later_failure_root unrelated-later-failure)
file(MAKE_DIRECTORY "${later_failure_root}/source")
file(MAKE_DIRECTORY "${later_failure_root}/outer")
file(WRITE "${later_failure_root}/outer/CMakeCache.txt"
     "CMAKE_GENERATOR:INTERNAL=Task22 Invalid Generator\n")
execute_process(
  COMMAND "${CMAKE_COMMAND}"
          "-DHUNDUN_SOURCE_ROOT=${later_failure_root}/source"
          "-DHUNDUN_OUTER_BINARY_DIR=${later_failure_root}/outer"
          "-DHUNDUN_NESTED_BINARY_DIR=${later_failure_root}/outer/task22-tests-off"
          "-DHUNDUN_ARCHIVE_DIR=${later_failure_root}/outer/task22-tests-off/archives"
          "-DHUNDUN_NESTED_CONFIG=Release"
          -P "${HUNDUN_SETUP_SCRIPT}"
  RESULT_VARIABLE later_setup_result
  OUTPUT_VARIABLE later_setup_stdout
  ERROR_VARIABLE later_setup_stderr)
task22_rejection_matches(
  later_setup_result later_setup_stdout later_setup_stderr
  "Task22 tests-off nested binary directory must be the dedicated outer/task22-tests-off path"
  later_failure_matches_guard)
if(later_setup_result EQUAL 0 OR later_failure_matches_guard)
  message(FATAL_ERROR
    "Task22 unrelated later failure unexpectedly satisfied a path guard:"
    "\n${later_setup_stdout}${later_setup_stderr}")
endif()
string(FIND "${later_setup_stdout}\n${later_setup_stderr}"
       "Task22 tests-off configure failed:" later_failure_position)
if(later_failure_position EQUAL -1)
  message(FATAL_ERROR
    "Task22 unrelated later-failure probe did not reach configure failure:"
    "\n${later_setup_stdout}${later_setup_stderr}")
endif()

task22_require_owned_root()
file(REMOVE_RECURSE "${contract_root}")
if(EXISTS "${contract_root}" OR IS_SYMLINK "${contract_root}")
  message(FATAL_ERROR
    "Task22 path-contract failed to remove its invocation-owned root")
endif()
