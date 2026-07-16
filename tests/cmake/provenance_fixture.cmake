# SPDX-License-Identifier: Apache-2.0

include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/HundunProvenanceGuard.cmake")

if(NOT DEFINED HUNDUN_PROVENANCE_CASE)
  set(HUNDUN_PROVENANCE_CASE clean)
endif()

set(fixture_base "$ENV{TMPDIR}")
if(fixture_base STREQUAL "")
  set(fixture_base "/tmp")
endif()
set(fixture_root
  "${fixture_base}/hundun-flow-provenance-fixture-${HUNDUN_PROVENANCE_CASE}")
set(fixture_prefix "${fixture_base}/hundun-flow-provenance-fixture")
set(clean_root "${fixture_root}/clean")
set(fortran_root "${fixture_root}/fortran")
set(token_root "${fixture_root}/token")
set(uppercase_root "${fixture_root}/uppercase")

file(REMOVE_RECURSE "${fixture_root}")
file(MAKE_DIRECTORY
  "${clean_root}/cmake" "${fortran_root}" "${token_root}" "${uppercase_root}")
file(WRITE "${clean_root}/main.cpp"
  "// SPDX-License-Identifier: Apache-2.0\nint main() { return 0; }\n")
file(WRITE "${clean_root}/CMakeLists.txt"
  "# SPDX-License-Identifier: Apache-2.0\n"
  "cmake_minimum_required(VERSION 3.21)\n"
  "project(HundunFixture LANGUAGES CXX)\n")
file(COPY "${CMAKE_CURRENT_LIST_DIR}/../../cmake/HundunProvenanceGuard.cmake"
  DESTINATION "${clean_root}/cmake")
file(WRITE "${fortran_root}/legacy.F90" "program legacy\nend program legacy\n")
string(CONCAT forbidden_fixture_token "CoAsT" "_LeGaCy")
file(WRITE "${token_root}/CMakeLists.txt"
  "# SPDX-License-Identifier: Apache-2.0\n"
  "set(reference_name \"${forbidden_fixture_token}\")\n")
string(CONCAT uppercase_fixture_token "DoM" "xCh")
file(WRITE "${uppercase_root}/legacy.CpP"
  "const char* reference_name = \"${uppercase_fixture_token}\";\n")

function(hundun_verify_negative child_case expected_relative_path
         expected_message expected_token)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" "-DHUNDUN_PROVENANCE_CASE=${child_case}"
            -P "${CMAKE_CURRENT_LIST_FILE}"
    RESULT_VARIABLE child_result
    OUTPUT_VARIABLE child_stdout
    ERROR_VARIABLE child_stderr)

  if(child_result STREQUAL "0")
    message(FATAL_ERROR
      "Provenance child '${child_case}' unexpectedly succeeded")
  endif()

  string(CONCAT child_output "${child_stdout}" "${child_stderr}")
  set(expected_path
    "${fixture_prefix}-${child_case}/${expected_relative_path}")
  foreach(expected_fragment IN ITEMS
      "${expected_message}" "${expected_path}" "${expected_token}")
    if(NOT expected_fragment STREQUAL "")
      string(FIND "${child_output}" "${expected_fragment}" fragment_position)
      if(fragment_position EQUAL -1)
        message(FATAL_ERROR
          "Provenance child '${child_case}' did not emit expected diagnostic "
          "fragment '${expected_fragment}'. Output:\n${child_output}")
      endif()
    endif()
  endforeach()
endfunction()

if(HUNDUN_PROVENANCE_CASE STREQUAL "clean")
  hundun_assert_clean_tree("${clean_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "fortran")
  hundun_assert_clean_tree("${fortran_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "token")
  hundun_assert_clean_tree("${token_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "uppercase")
  hundun_assert_clean_tree("${uppercase_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_fortran")
  hundun_verify_negative(
    fortran "fortran/legacy.F90" "Forbidden legacy-language files" "")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_token")
  string(CONCAT expected_token "co" "ast" "_le" "gacy")
  hundun_verify_negative(
    token "token/CMakeLists.txt" "Forbidden provenance token" "${expected_token}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_uppercase")
  string(CONCAT expected_token "do" "mx" "ch")
  hundun_verify_negative(
    uppercase "uppercase/legacy.CpP" "Forbidden provenance token"
    "${expected_token}")
else()
  message(FATAL_ERROR
    "Unknown HUNDUN_PROVENANCE_CASE: ${HUNDUN_PROVENANCE_CASE}")
endif()
