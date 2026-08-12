# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

get_filename_component(HUNDUN_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(TEST_ROOT "${HUNDUN_ROOT}/build/version-dispatch-test")
file(REMOVE_RECURSE "${TEST_ROOT}")

set_property(GLOBAL PROPERTY HUNDUN_DISPATCH_FAILURES "")

function(record_failure message_text)
  set_property(GLOBAL APPEND_STRING PROPERTY HUNDUN_DISPATCH_FAILURES
    "\n${message_text}")
endfunction()

function(check_successful_selection case_name expected_version)
  set(build_dir "${TEST_ROOT}/${case_name}")
  set(configure_arguments
    -S "${HUNDUN_ROOT}"
    -B "${build_dir}"
    -DHUNDUN_BUILD_TESTS=OFF)
  if(ARGC GREATER 2)
    list(APPEND configure_arguments
      "-DHUNDUN_SOURCE_VERSION=${ARGV2}")
  endif()

  execute_process(
    COMMAND "${CMAKE_COMMAND}" ${configure_arguments}
    RESULT_VARIABLE configure_result
    OUTPUT_VARIABLE configure_stdout
    ERROR_VARIABLE configure_stderr)

  if(NOT configure_result EQUAL 0)
    record_failure(
      "${case_name}: configure failed (${configure_result})\n"
      "stdout:\n${configure_stdout}\n"
      "stderr:\n${configure_stderr}")
    return()
  endif()

  file(READ "${build_dir}/CMakeCache.txt" cache_text)
  if(NOT cache_text MATCHES
     "(^|\n)HUNDUN_SOURCE_VERSION:STRING=${expected_version}(\n|$)")
    record_failure(
      "${case_name}: root does not recognize HUNDUN_SOURCE_VERSION as a "
      "STRING cache selection for ${expected_version}")
  endif()

  if(NOT EXISTS "${build_dir}/versions/${expected_version}/CMakeFiles")
    record_failure(
      "${case_name}: HUNDUN_SOURCE_VERSION=${expected_version} did not "
      "select versions/${expected_version}")
  endif()
endfunction()

function(check_invalid_selection)
  set(build_dir "${TEST_ROOT}/invalid-v9")
  execute_process(
    COMMAND "${CMAKE_COMMAND}"
      -S "${HUNDUN_ROOT}"
      -B "${build_dir}"
      -DHUNDUN_BUILD_TESTS=OFF
      -DHUNDUN_SOURCE_VERSION=v9
    RESULT_VARIABLE configure_result
    OUTPUT_VARIABLE configure_stdout
    ERROR_VARIABLE configure_stderr)
  set(configure_output "${configure_stdout}\n${configure_stderr}")

  if(configure_result EQUAL 0)
    record_failure(
      "invalid-v9: configure succeeded; unsupported "
      "HUNDUN_SOURCE_VERSION was not rejected")
  elseif(NOT configure_output MATCHES
         "unsupported HUNDUN_SOURCE_VERSION")
    record_failure(
      "invalid-v9: configure failed without the required diagnostic\n"
      "${configure_output}")
  endif()
endfunction()

check_successful_selection(default v0.4)
check_successful_selection(explicit-v0.4 v0.4 v0.4)
check_successful_selection(explicit-v0.3 v0.3 v0.3)
check_invalid_selection()

get_property(failures GLOBAL PROPERTY HUNDUN_DISPATCH_FAILURES)
if(NOT failures STREQUAL "")
  message(FATAL_ERROR "Version dispatch assertions failed:${failures}")
endif()

file(REMOVE_RECURSE "${TEST_ROOT}")
message(STATUS "Version dispatch configure checks passed")
