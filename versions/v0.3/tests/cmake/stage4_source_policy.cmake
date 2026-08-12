# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

if(NOT DEFINED HUNDUN_SOURCE_DIR)
  message(FATAL_ERROR "HUNDUN_SOURCE_DIR is required")
endif()

get_filename_component(HUNDUN_SOURCE_DIR "${HUNDUN_SOURCE_DIR}" ABSOLUTE)

function(stage4_reject label reason)
  message(FATAL_ERROR "Stage 4 source policy rejected ${label}: ${reason}")
endfunction()

function(stage4_parse_prefix_allowlist label contents output_variable)
  string(REGEX MATCH
    "stage4_product_prefix_allowlist=([a-z0-9_,]+)"
    allowlist_match "${contents}")
  if(allowlist_match STREQUAL "")
    stage4_reject("${label}" "missing Stage 4 product prefix allowlist")
  endif()
  set(allowlist_csv "${CMAKE_MATCH_1}")
  string(REPLACE "," ";" allowlist "${allowlist_csv}")
  list(REMOVE_DUPLICATES allowlist)
  set(${output_variable} "${allowlist}" PARENT_SCOPE)
endfunction()

function(stage4_check_product_name label name allowlist)
  if(NOT name MATCHES "^([a-z0-9]+)_[a-z0-9_]+\\.(h|hpp|cpp)$")
    stage4_reject("${label}" "unregistered product filename ${name}")
  endif()
  set(prefix "${CMAKE_MATCH_1}")
  if(NOT prefix IN_LIST allowlist)
    stage4_reject("${label}" "unregistered product prefix ${prefix}_")
  endif()
endfunction()

function(stage4_check_public_header label contents)
  set(forbidden_public_patterns
    "Python\\.h"
    "pybind(11)?"
    "Cantera::"
    "[#]include[ \t]*[<\"][^>\"]*cantera/"
    "N_Vector"
    "SUN(Context|Matrix|LinearSolver|NonlinearSolver)"
    "SUNDIALS_[A-Za-z0-9_]+"
    "[#]include[ \t]*[<\"][^>\"]*sundials/")
  foreach(pattern IN LISTS forbidden_public_patterns)
    if(contents MATCHES "${pattern}")
      stage4_reject("${label}" "public header exposes forbidden third-party type")
    endif()
  endforeach()
endfunction()

function(stage4_check_private_source label contents)
  string(TOLOWER "${contents}" lower_contents)
  foreach(pattern IN ITEMS
      "boffin" "coast_software" "compress_boffin" "/coast/" "\\coast\\")
    string(FIND "${lower_contents}" "${pattern}" found)
    if(NOT found EQUAL -1)
      stage4_reject("${label}" "private source or path reference")
    endif()
  endforeach()
endfunction()

function(stage4_check_plugin_abi label contents)
  string(REGEX MATCHALL
    "HUNDUN_PLUGIN_[A-Z0-9_]*ABI_V[0-9]+"
    plugin_abi_symbols "${contents}")
  foreach(symbol IN LISTS plugin_abi_symbols)
    if(NOT symbol STREQUAL "HUNDUN_PLUGIN_METADATA_ABI_V1")
      stage4_reject("${label}" "second plugin ABI ${symbol}")
    endif()
  endforeach()
  if(contents MATCHES "hundun_plugin_entry_v[2-9][0-9]*")
    stage4_reject("${label}" "second plugin entry ABI")
  endif()
endfunction()

function(stage4_check_cantera_provenance label contents)
  string(JSON provenance_type ERROR_VARIABLE provenance_error
    TYPE "${contents}")
  if(NOT provenance_error STREQUAL "NOTFOUND" OR
     NOT provenance_type STREQUAL "OBJECT")
    stage4_reject("${label}" "UPSTREAM.json is not a valid JSON object")
  endif()
  foreach(required IN ITEMS
      "\"schema\": \"hundun.cantera.upstream.v1\""
      "\"tag\": \"v3.2.0\""
      "\"commit\": \"4a8358eb80cfeb50474386b5f9ec0b3a83519889\""
      "\"archive_sha256\": \"a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b\""
      "\"artifact_sha256\": \"093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760\""
      "\"license_file\": \"LICENSES/cantera-BSD-3-Clause.txt\""
      "\"mechanisms_bundled\": false")
    string(FIND "${contents}" "${required}" found)
    if(found EQUAL -1)
      stage4_reject("${label}" "missing or changed provenance field ${required}")
    endif()
  endforeach()
  foreach(dependency IN ITEMS fmt yaml-cpp sundials eigen boost)
    string(FIND "${contents}" "\"name\": \"${dependency}\"" found)
    if(found EQUAL -1)
      stage4_reject("${label}" "missing dependency ${dependency}")
    endif()
  endforeach()
  string(REGEX MATCHALL "\"license\": \"[^\"]+\"" licenses "${contents}")
  list(LENGTH licenses license_count)
  if(license_count LESS 5)
    stage4_reject("${label}" "dependency license closure is incomplete")
  endif()
endfunction()

if(DEFINED HUNDUN_STAGE4_POLICY_VERIFY_FIXTURE)
  execute_process(
    COMMAND "${CMAKE_COMMAND}"
            "-DHUNDUN_SOURCE_DIR=${HUNDUN_SOURCE_DIR}"
            "-DHUNDUN_STAGE4_POLICY_FIXTURE=${HUNDUN_STAGE4_POLICY_VERIFY_FIXTURE}"
            -P "${CMAKE_CURRENT_LIST_FILE}"
    RESULT_VARIABLE fixture_result
    OUTPUT_VARIABLE fixture_stdout
    ERROR_VARIABLE fixture_stderr)
  if(fixture_result EQUAL 0)
    message(FATAL_ERROR
      "Stage 4 source-policy fixture was unexpectedly accepted: "
      "${HUNDUN_STAGE4_POLICY_VERIFY_FIXTURE}")
  endif()
  set(fixture_output "${fixture_stdout}${fixture_stderr}")
  string(FIND "${fixture_output}"
    "Stage 4 source policy rejected ${HUNDUN_STAGE4_POLICY_VERIFY_FIXTURE}:"
    fixture_marker)
  if(fixture_marker EQUAL -1)
    message(FATAL_ERROR
      "Stage 4 source-policy fixture failed for an unexpected reason: "
      "${fixture_output}")
  endif()
  return()
endif()

set(stage4_ledger_path
  "${HUNDUN_SOURCE_DIR}/docs/numerics/stage4-capability-ledger.md")
if(NOT EXISTS "${stage4_ledger_path}")
  stage4_reject("capability ledger" "missing required file")
endif()
file(READ "${stage4_ledger_path}" stage4_ledger_contents)
stage4_parse_prefix_allowlist(
  "capability ledger" "${stage4_ledger_contents}" stage4_prefix_allowlist)
string(CONCAT stage4_python_header "#include <Py" "thon.h>")

if(DEFINED HUNDUN_STAGE4_POLICY_FIXTURE)
  if(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "missing_allowlist")
    stage4_parse_prefix_allowlist(
      "missing_allowlist" "# empty fixture ledger" ignored_allowlist)
  elseif(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "public_cantera_type")
    stage4_check_public_header(
      "public_cantera_type" "struct Leak { Cantera::Solution* value; };")
  elseif(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "public_sundials_type")
    stage4_check_public_header(
      "public_sundials_type" "struct Leak { N_Vector value; };")
  elseif(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "public_python_header")
    stage4_check_public_header(
      "public_python_header" "${stage4_python_header}")
  elseif(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "public_pybind_type")
    stage4_check_public_header(
      "public_pybind_type" "pybind11::object value;")
  elseif(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "second_plugin_abi")
    stage4_check_plugin_abi(
      "second_plugin_abi" "#define HUNDUN_PLUGIN_MODEL_ABI_V1 1u")
  elseif(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "private_path")
    stage4_check_private_source(
      "private_path" "/home/example/Coast_software/private/source.f90")
  elseif(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "unregistered_prefix")
    stage4_check_product_name(
      "unregistered_prefix" "react_stage4_backend.cpp"
      "${stage4_prefix_allowlist}")
  elseif(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "cantera_archive_hash")
    file(READ "${HUNDUN_SOURCE_DIR}/third_party/cantera/UPSTREAM.json"
      cantera_fixture)
    string(REPLACE
      "a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b"
      "b94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b"
      cantera_fixture "${cantera_fixture}")
    stage4_check_cantera_provenance("cantera_archive_hash" "${cantera_fixture}")
  elseif(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "cantera_missing_license")
    file(READ "${HUNDUN_SOURCE_DIR}/third_party/cantera/UPSTREAM.json"
      cantera_fixture)
    string(REGEX REPLACE
      "\"license\": \"[^\"]+\"" "\"omitted\": true"
      cantera_fixture "${cantera_fixture}")
    stage4_check_cantera_provenance(
      "cantera_missing_license" "${cantera_fixture}")
  elseif(HUNDUN_STAGE4_POLICY_FIXTURE STREQUAL "cantera_product_source")
    stage4_reject("cantera_product_source"
      "Cantera source is forbidden under src")
  else()
    message(FATAL_ERROR
      "unknown Stage 4 source-policy fixture: ${HUNDUN_STAGE4_POLICY_FIXTURE}")
  endif()
  return()
endif()

file(GLOB stage4_public_headers
  "${HUNDUN_SOURCE_DIR}/include/hundun/*.h"
  "${HUNDUN_SOURCE_DIR}/include/hundun/*.hpp")
file(GLOB stage4_product_sources
  "${HUNDUN_SOURCE_DIR}/src/*.cpp"
  "${HUNDUN_SOURCE_DIR}/src/*.hpp")

foreach(path IN LISTS stage4_public_headers stage4_product_sources)
  get_filename_component(name "${path}" NAME)
  stage4_check_product_name("${name}" "${name}" "${stage4_prefix_allowlist}")
  file(READ "${path}" contents)
  stage4_check_private_source("${name}" "${contents}")
  stage4_check_plugin_abi("${name}" "${contents}")
endforeach()

set(stage4_cantera_upstream
  "${HUNDUN_SOURCE_DIR}/third_party/cantera/UPSTREAM.json")
set(stage4_cantera_license
  "${HUNDUN_SOURCE_DIR}/LICENSES/cantera-BSD-3-Clause.txt")
if(NOT EXISTS "${stage4_cantera_upstream}")
  stage4_reject("Cantera provenance" "missing UPSTREAM.json")
endif()
if(NOT EXISTS "${stage4_cantera_license}")
  stage4_reject("Cantera provenance" "missing Cantera license")
endif()
file(SHA256 "${stage4_cantera_license}" stage4_cantera_license_sha256)
if(NOT stage4_cantera_license_sha256 STREQUAL
   "e92980b9712ce20e73898a97b0116889e84e07f548d6be8591e87dcad79c41bb")
  stage4_reject("Cantera provenance" "Cantera license identity changed")
endif()
file(READ "${stage4_cantera_upstream}" stage4_cantera_contents)
stage4_check_cantera_provenance(
  "Cantera provenance" "${stage4_cantera_contents}")

file(GLOB stage4_forbidden_cantera_sources
  "${HUNDUN_SOURCE_DIR}/src/*cantera*"
  "${HUNDUN_SOURCE_DIR}/include/hundun/*cantera*")
foreach(path IN LISTS stage4_forbidden_cantera_sources)
  get_filename_component(name "${path}" NAME)
  if(NOT name STREQUAL "chem_cantera_backend.cpp")
    stage4_reject("Cantera product source"
      "unapproved Cantera-named product file ${name}")
  endif()
endforeach()

foreach(path IN LISTS stage4_public_headers)
  get_filename_component(name "${path}" NAME)
  file(READ "${path}" contents)
  stage4_check_public_header("${name}" "${contents}")
endforeach()
