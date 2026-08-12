# SPDX-License-Identifier: Apache-2.0
cmake_minimum_required(VERSION 3.21)
function(classify text output)
  set(result accepted)
  string(REPLACE "\n" ";" lines "${text}")
  foreach(line IN LISTS lines)
    if(line STREQUAL "")
      continue()
    endif()
    string(REPLACE "\t" ";" columns "${line}")
    list(LENGTH columns width)
    if(NOT width EQUAL 4)
      set(result invalid-width)
      break()
    endif()
    list(GET columns 0 path)
    list(GET columns 3 relation)
    if(path STREQUAL "path")
      continue()
    endif()
    if("${path}" MATCHES "(^|/)tests(/|$)|(^|/)\\.superpowers(/|$)|private|token" OR NOT "${relation}" MATCHES "^(identical|product_tests_off_preset|product_only_override)$")
      set(result illegal-projection)
      break()
    endif()
  endforeach()
  set(${output} "${result}" PARENT_SCOPE)
endfunction()
file(READ "${HUNDUN_STAGE3_PRODUCT_PROJECTION_MANIFEST}" manifest)
classify("${manifest}" result)
if(NOT result STREQUAL accepted)
  message(FATAL_ERROR "Stage 3 product projection manifest rejected: ${result}")
endif()
classify("${manifest}\ntests/private-token.cpp\ta\tb\tidentical" mutation)
if(NOT mutation STREQUAL illegal-projection)
  message(FATAL_ERROR "Stage 3 projection illegal-path mutation survived")
endif()
file(READ "${HUNDUN_STAGE3_PRODUCT_PROJECTOR}" projector)
string(FIND "${projector}" "product-projection-manifest" projector_reference)
if(projector_reference LESS 0)
  message(FATAL_ERROR "Stage 3 final projector does not bind the manifest")
endif()

set(fixture_root
    "${CMAKE_CURRENT_BINARY_DIR}/stage3-product-projection-contract-fixture")
set(owner_marker "${fixture_root}/.hundun-stage3-product-projection-contract")
set(expected_owner "HUNDUN Stage 3 product projection contract\n")
if(EXISTS "${fixture_root}")
  if(NOT EXISTS "${owner_marker}" OR IS_SYMLINK "${owner_marker}")
    message(FATAL_ERROR "refusing to replace an unowned projection fixture")
  endif()
  file(READ "${owner_marker}" observed_owner)
  if(NOT observed_owner STREQUAL expected_owner)
    message(FATAL_ERROR "projection fixture owner mismatch")
  endif()
  file(REMOVE_RECURSE "${fixture_root}")
endif()
file(MAKE_DIRECTORY "${fixture_root}")
file(WRITE "${owner_marker}" "${expected_owner}")

function(run_checked)
  execute_process(
    COMMAND ${ARGN}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "fixture command failed: ${ARGN}\n${output}${error}")
  endif()
endfunction()

function(write_fixture repository kind)
  file(MAKE_DIRECTORY "${repository}/include/hundun")
  if(kind STREQUAL "governance")
    file(MAKE_DIRECTORY "${repository}/.superpowers")
    file(MAKE_DIRECTORY "${repository}/docs/development")
    file(WRITE "${repository}/CMakeLists.txt" "governance candidate\n")
    file(WRITE "${repository}/CMakePresets.json" "governance tests on\n")
    file(WRITE "${repository}/VERSION" "0.1.0\n")
    file(WRITE "${repository}/include/hundun/new.hpp" "stage3 product header\n")
    file(WRITE "${repository}/docs/development/source-policy.md"
      "governance-only source policy\n")
    file(WRITE "${repository}/.superpowers/governance-only.md" "do not copy\n")
  elseif(kind STREQUAL "product")
    file(WRITE "${repository}/CMakeLists.txt" "baseline product\n")
    file(WRITE "${repository}/CMakePresets.json" "product tests off\n")
    file(WRITE "${repository}/VERSION" "0.1.0\n")
  else()
    message(FATAL_ERROR "unknown fixture kind: ${kind}")
  endif()
  run_checked(git -C "${repository}" init -q)
  run_checked(git -C "${repository}" config user.name "Stage 3 Fixture")
  run_checked(git -C "${repository}" config user.email "fixture@example.invalid")
  run_checked(git -C "${repository}" add .)
  run_checked(git -C "${repository}" commit -q -m "fixture")
endfunction()

function(blob_id repository path output)
  execute_process(
    COMMAND git -C "${repository}" rev-parse "HEAD:${path}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE blob
    ERROR_VARIABLE error
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "cannot read fixture blob ${path}: ${error}")
  endif()
  set(${output} "${blob}" PARENT_SCOPE)
endfunction()

set(governance_root "${fixture_root}/governance")
set(product_root "${fixture_root}/product")
write_fixture("${governance_root}" governance)
write_fixture("${product_root}" product)
blob_id("${product_root}" "CMakeLists.txt" cmake_blob)
blob_id("${product_root}" "CMakePresets.json" preset_blob)
blob_id("${product_root}" "VERSION" version_blob)
string(CONCAT fixture_manifest
  "path\tproduct_blob\tgovernance_blob\trelation\n"
  "CMakeLists.txt\t${cmake_blob}\t${cmake_blob}\tidentical\n"
  "CMakePresets.json\t${preset_blob}\t${preset_blob}\tproduct_tests_off_preset\n"
  "VERSION\t${version_blob}\t${version_blob}\tidentical\n")
file(WRITE "${governance_root}/.superpowers/base.tsv" "${fixture_manifest}")
run_checked(git -C "${governance_root}" add .superpowers/base.tsv)
run_checked(git -C "${governance_root}" commit -q -m "add base manifest")

execute_process(
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_GOVERNANCE_ROOT=${governance_root}"
    "-DHUNDUN_PRODUCT_ROOT=${product_root}"
    "-DHUNDUN_BASE_MANIFEST=.superpowers/base.tsv"
    "-DHUNDUN_OUTPUT_MANIFEST=.superpowers/final.tsv"
    -P "${HUNDUN_STAGE3_PRODUCT_PROJECTOR}"
  RESULT_VARIABLE projection_result
  OUTPUT_VARIABLE projection_output
  ERROR_VARIABLE projection_error)
if(NOT projection_result EQUAL 0)
  message(FATAL_ERROR
    "valid product projection failed:\n${projection_output}${projection_error}")
endif()
foreach(expected IN ITEMS
    "${product_root}/include/hundun/new.hpp"
    "${governance_root}/.superpowers/final.tsv")
  if(NOT EXISTS "${expected}")
    message(FATAL_ERROR "product projection omitted ${expected}")
  endif()
endforeach()
file(READ "${product_root}/CMakeLists.txt" projected_cmake)
file(READ "${product_root}/CMakePresets.json" projected_preset)
file(READ "${product_root}/VERSION" projected_version)
file(READ "${governance_root}/.superpowers/final.tsv" final_manifest)
if(NOT projected_cmake STREQUAL "governance candidate\n")
  message(FATAL_ERROR "product projection did not copy governance content")
endif()
if(NOT projected_preset STREQUAL "product tests off\n")
  message(FATAL_ERROR "product projection did not preserve tests-off preset")
endif()
if(NOT projected_version STREQUAL "0.2.0\n")
  message(FATAL_ERROR "product projection did not set product VERSION=0.2.0")
endif()
if(EXISTS "${product_root}/.superpowers" OR
   EXISTS "${product_root}/tests")
  message(FATAL_ERROR "governance path leaked into product projection")
endif()
if(EXISTS "${product_root}/docs/development/source-policy.md")
  message(FATAL_ERROR "governance source-policy document leaked into product")
endif()
string(FIND "${final_manifest}"
  "docs/development/source-policy.md" source_policy_index)
if(NOT source_policy_index EQUAL -1)
  message(FATAL_ERROR "final manifest contains governance source-policy document")
endif()
foreach(fragment IN ITEMS
    "CMakeLists.txt\t"
    "CMakePresets.json\t"
    "VERSION\t"
    "include/hundun/new.hpp\t"
    "product_tests_off_preset"
    "product_only_override")
  string(FIND "${final_manifest}" "${fragment}" fragment_index)
  if(fragment_index LESS 0)
    message(FATAL_ERROR "final projection manifest omitted '${fragment}'")
  endif()
endforeach()

# A tracked product path not authenticated by the base manifest must fail before
# changing any baseline file.
run_checked(git -C "${governance_root}" clean -fdq)
run_checked(git -C "${product_root}" reset -q --hard HEAD)
run_checked(git -C "${product_root}" clean -fdq)
file(WRITE "${product_root}/unexpected.txt" "unexpected\n")
run_checked(git -C "${product_root}" add unexpected.txt)
run_checked(git -C "${product_root}" commit -q -m "unexpected path mutation")
execute_process(
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_GOVERNANCE_ROOT=${governance_root}"
    "-DHUNDUN_PRODUCT_ROOT=${product_root}"
    "-DHUNDUN_BASE_MANIFEST=.superpowers/base.tsv"
    "-DHUNDUN_OUTPUT_MANIFEST=.superpowers/rejected.tsv"
    -P "${HUNDUN_STAGE3_PRODUCT_PROJECTOR}"
  RESULT_VARIABLE rejected_result
  OUTPUT_VARIABLE rejected_output
  ERROR_VARIABLE rejected_error)
if(rejected_result EQUAL 0)
  message(FATAL_ERROR "projector accepted an unexpected tracked product path")
endif()
string(CONCAT rejected_diagnostic "${rejected_output}" "${rejected_error}")
if(NOT rejected_diagnostic MATCHES "unexpected tracked product path")
  message(FATAL_ERROR
    "unexpected-path mutation failed for the wrong reason:\n${rejected_diagnostic}")
endif()
file(READ "${product_root}/CMakeLists.txt" rejected_sentinel)
if(NOT rejected_sentinel STREQUAL "baseline product\n")
  message(FATAL_ERROR "rejected projection partially modified the product")
endif()

# A product blob that no longer matches the authenticated baseline must also
# fail without copying governance content.
run_checked(git -C "${product_root}" reset -q --hard HEAD~1)
string(REPLACE
  "CMakeLists.txt\t${cmake_blob}"
  "CMakeLists.txt\t0000000000000000000000000000000000000000"
  bad_manifest "${fixture_manifest}")
file(WRITE "${governance_root}/.superpowers/bad-base.tsv" "${bad_manifest}")
run_checked(git -C "${governance_root}" add .superpowers/bad-base.tsv)
run_checked(git -C "${governance_root}" commit -q -m "bad blob mutation")
execute_process(
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_GOVERNANCE_ROOT=${governance_root}"
    "-DHUNDUN_PRODUCT_ROOT=${product_root}"
    "-DHUNDUN_BASE_MANIFEST=.superpowers/bad-base.tsv"
    "-DHUNDUN_OUTPUT_MANIFEST=.superpowers/rejected-blob.tsv"
    -P "${HUNDUN_STAGE3_PRODUCT_PROJECTOR}"
  RESULT_VARIABLE bad_blob_result
  OUTPUT_VARIABLE bad_blob_output
  ERROR_VARIABLE bad_blob_error)
if(bad_blob_result EQUAL 0)
  message(FATAL_ERROR "projector accepted a mismatched baseline product blob")
endif()
string(CONCAT bad_blob_diagnostic "${bad_blob_output}" "${bad_blob_error}")
if(NOT bad_blob_diagnostic MATCHES "product baseline blob mismatch")
  message(FATAL_ERROR
    "bad-blob mutation failed for the wrong reason:\n${bad_blob_diagnostic}")
endif()
file(READ "${product_root}/CMakeLists.txt" bad_blob_sentinel)
if(NOT bad_blob_sentinel STREQUAL "baseline product\n")
  message(FATAL_ERROR "bad-blob rejection partially modified the product")
endif()

run_checked(git -C "${governance_root}" reset -q --hard HEAD~1)
string(APPEND illegal_manifest
  "${fixture_manifest}"
  "tests/private-token.cpp\t${cmake_blob}\t${cmake_blob}\tidentical\n")
file(WRITE "${governance_root}/.superpowers/illegal-base.tsv"
  "${illegal_manifest}")
run_checked(git -C "${governance_root}" add .superpowers/illegal-base.tsv)
run_checked(git -C "${governance_root}" commit -q -m "illegal path mutation")
execute_process(
  COMMAND "${CMAKE_COMMAND}"
    "-DHUNDUN_GOVERNANCE_ROOT=${governance_root}"
    "-DHUNDUN_PRODUCT_ROOT=${product_root}"
    "-DHUNDUN_BASE_MANIFEST=.superpowers/illegal-base.tsv"
    "-DHUNDUN_OUTPUT_MANIFEST=.superpowers/rejected-path.tsv"
    -P "${HUNDUN_STAGE3_PRODUCT_PROJECTOR}"
  RESULT_VARIABLE illegal_path_result
  OUTPUT_VARIABLE illegal_path_output
  ERROR_VARIABLE illegal_path_error)
if(illegal_path_result EQUAL 0)
  message(FATAL_ERROR "projector accepted an illegal manifest path")
endif()
string(CONCAT illegal_path_diagnostic
  "${illegal_path_output}" "${illegal_path_error}")
if(NOT illegal_path_diagnostic MATCHES "illegal product projection path")
  message(FATAL_ERROR
    "illegal-path mutation failed for the wrong reason:\n"
    "${illegal_path_diagnostic}")
endif()
