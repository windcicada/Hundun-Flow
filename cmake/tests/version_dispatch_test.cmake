# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

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

function(check_sanitizer_configuration case_name option_name sanitizer_flag)
  set(build_dir "${TEST_ROOT}/${case_name}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}"
      -G "Unix Makefiles"
      -S "${HUNDUN_ROOT}"
      -B "${build_dir}"
      -DCMAKE_BUILD_TYPE=Debug
      -DHUNDUN_BUILD_TESTS=OFF
      -DHUNDUN_SOURCE_VERSION=v0.4
      "-D${option_name}=ON"
    RESULT_VARIABLE configure_result
    OUTPUT_VARIABLE configure_stdout
    ERROR_VARIABLE configure_stderr)
  set(configure_output "${configure_stdout}\n${configure_stderr}")

  if(NOT configure_result EQUAL 0)
    record_failure(
      "${case_name}: configure failed (${configure_result})\n"
      "${configure_output}")
    return()
  endif()

  if(configure_output MATCHES
     "Manually-specified variables were not used")
    record_failure(
      "${case_name}: configure reported manually-specified unused variables\n"
      "${configure_output}")
  endif()

  set(core_flags_file
    "${build_dir}/versions/v0.4/CMakeFiles/hundun_v04_core.dir/flags.make")
  if(NOT EXISTS "${core_flags_file}")
    record_failure(
      "${case_name}: core compile flags were not generated at "
      "${core_flags_file}")
  else()
    file(READ "${core_flags_file}" core_flags)
    if(NOT core_flags MATCHES "(^|[ =])${sanitizer_flag}([ \n]|$)")
      record_failure(
        "${case_name}: core compile flags do not contain "
        "${sanitizer_flag}\n${core_flags}")
    endif()
  endif()

  set(hundun_link_file
    "${build_dir}/versions/v0.4/CMakeFiles/hundun.dir/link.txt")
  if(NOT EXISTS "${hundun_link_file}")
    record_failure(
      "${case_name}: hundun link flags were not generated at "
      "${hundun_link_file}")
  else()
    file(READ "${hundun_link_file}" hundun_link)
    if(NOT hundun_link MATCHES "(^|[ ])${sanitizer_flag}([ \n]|$)")
      record_failure(
        "${case_name}: hundun link flags do not contain "
        "${sanitizer_flag}\n${hundun_link}")
    endif()
  endif()
endfunction()

check_successful_selection(default v0.4)
check_successful_selection(explicit-v0.4 v0.4 v0.4)
check_successful_selection(explicit-v0.3 v0.3 v0.3)
check_invalid_selection()
check_sanitizer_configuration(
  v0.4-asan HUNDUN_ENABLE_ASAN -fsanitize=address)
check_sanitizer_configuration(
  v0.4-ubsan HUNDUN_ENABLE_UBSAN -fsanitize=undefined)

get_property(failures GLOBAL PROPERTY HUNDUN_DISPATCH_FAILURES)
if(NOT failures STREQUAL "")
  message(FATAL_ERROR "Version dispatch assertions failed:${failures}")
endif()

file(REMOVE_RECURSE "${TEST_ROOT}")
message(STATUS "Version dispatch configure checks passed")
