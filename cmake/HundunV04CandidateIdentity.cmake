# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

function(hundun_v04_write_build_identity)
  cmake_parse_arguments(ARG "" "OUTPUT;SOURCE_ROOT;BINARY;TESTS_ON_BINARY;TESTS_OFF_BINARY" "INPUTS" ${ARGN})
  foreach(required OUTPUT SOURCE_ROOT BINARY TESTS_ON_BINARY TESTS_OFF_BINARY)
    if(NOT ARG_${required})
      message(FATAL_ERROR "hundun_v04_write_build_identity requires ${required}")
    endif()
  endforeach()
  if(EXISTS "${ARG_OUTPUT}")
    message(FATAL_ERROR "candidate identity is immutable: ${ARG_OUTPUT} already exists")
  endif()
  execute_process(
    COMMAND git -C "${ARG_SOURCE_ROOT}" status --porcelain=v1
    OUTPUT_VARIABLE dirty OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE git_status)
  if(NOT git_status EQUAL 0 OR NOT dirty STREQUAL "")
    message(FATAL_ERROR "candidate identity requires a clean worktree")
  endif()
  execute_process(COMMAND git -C "${ARG_SOURCE_ROOT}" rev-parse HEAD
    OUTPUT_VARIABLE head OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY)
  execute_process(COMMAND git -C "${ARG_SOURCE_ROOT}" rev-parse HEAD^{tree}
    OUTPUT_VARIABLE tree OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY)
  foreach(binary BINARY TESTS_ON_BINARY TESTS_OFF_BINARY)
    if(NOT EXISTS "${ARG_${binary}}")
      message(FATAL_ERROR "missing candidate binary: ${ARG_${binary}}")
    endif()
    file(SHA256 "${ARG_${binary}}" ${binary}_sha256)
    file(REAL_PATH "${ARG_${binary}}" ${binary}_path)
    string(REPLACE "\\" "\\\\" ${binary}_path_escaped "${${binary}_path}")
    string(REPLACE "\"" "\\\"" ${binary}_path_escaped "${${binary}_path_escaped}")
  endforeach()
  set(input_json "")
  foreach(input IN LISTS ARG_INPUTS)
    if(NOT EXISTS "${input}")
      message(FATAL_ERROR "missing candidate input: ${input}")
    endif()
    file(SHA256 "${input}" input_sha256)
    file(REAL_PATH "${input}" input_path)
    string(REPLACE "\\" "\\\\" escaped "${input_path}")
    string(REPLACE "\"" "\\\"" escaped "${escaped}")
    if(NOT input_json STREQUAL "")
      string(APPEND input_json ",")
    endif()
    string(APPEND input_json "\n    {\"path\":\"${escaped}\",\"sha256\":\"${input_sha256}\"}")
  endforeach()
  set(compiler "${CMAKE_CXX_COMPILER}")
  set(linker "${CMAKE_LINKER}")
  string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type_upper)
  set(flags "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${build_type_upper}} ${CMAKE_EXE_LINKER_FLAGS}")
  foreach(value compiler linker flags)
    string(REPLACE "\\" "\\\\" ${value}_escaped "${${value}}")
    string(REPLACE "\"" "\\\"" ${value}_escaped "${${value}_escaped}")
  endforeach()
  file(WRITE "${ARG_OUTPUT}.pending"
    "{\n"
    "  \"schema\":\"HUNDUN_V04_BUILD_IDENTITY_V1\",\n"
    "  \"head\":\"${head}\",\n"
    "  \"tree\":\"${tree}\",\n"
    "  \"compiler\":\"${compiler_escaped}\",\n"
    "  \"linker\":\"${linker_escaped}\",\n"
    "  \"flags\":\"${flags_escaped}\",\n"
    "  \"compiler_id\":\"${CMAKE_CXX_COMPILER_ID}\",\n"
    "  \"compiler_version\":\"${CMAKE_CXX_COMPILER_VERSION}\",\n"
    "  \"build_type\":\"${CMAKE_BUILD_TYPE}\",\n"
    "  \"binary_path\":\"${BINARY_path_escaped}\",\n"
    "  \"binary_sha256\":\"${BINARY_sha256}\",\n"
    "  \"tests_on_path\":\"${TESTS_ON_BINARY_path_escaped}\",\n"
    "  \"tests_on_sha256\":\"${TESTS_ON_BINARY_sha256}\",\n"
    "  \"tests_off_path\":\"${TESTS_OFF_BINARY_path_escaped}\",\n"
    "  \"tests_off_sha256\":\"${TESTS_OFF_BINARY_sha256}\",\n"
    "  \"inputs\":[${input_json}\n  ]\n"
    "}\n")
  file(RENAME "${ARG_OUTPUT}.pending" "${ARG_OUTPUT}" NO_REPLACE)
endfunction()
