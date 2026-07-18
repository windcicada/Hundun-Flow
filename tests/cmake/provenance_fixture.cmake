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
set(source_clean_root "${fixture_root}/source-clean")
set(source_python_file_root "${fixture_root}/source-python-file")
set(source_python_header_root "${fixture_root}/source-python-header")
set(source_python_package_root "${fixture_root}/source-python-package")
set(source_python_codegen_root "${fixture_root}/source-python-codegen")
set(source_fetch_root "${fixture_root}/source-fetch")
set(source_external_root "${fixture_root}/source-external")
set(source_download_root "${fixture_root}/source-download")
set(source_missing_root "${fixture_root}/source-missing")

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

file(MAKE_DIRECTORY
  "${source_clean_root}/src"
  "${source_python_file_root}/tools"
  "${source_python_header_root}/src"
  "${source_python_package_root}/cmake"
  "${source_python_codegen_root}/cmake"
  "${source_fetch_root}/cmake"
  "${source_external_root}/cmake"
  "${source_download_root}/cmake")
file(WRITE "${source_clean_root}/CMakeLists.txt"
  "cmake_minimum_required(VERSION 3.21)\n"
  "project(CleanPublicTree LANGUAGES CXX)\n")
file(WRITE "${source_clean_root}/src/main.cpp"
  "int main() { return 0; }\n")
file(WRITE "${source_python_file_root}/tools/generate.py"
  "# fixture-only public script\n")
string(CONCAT python_header_name "Py" "ThOn" ".h")
file(WRITE "${source_python_header_root}/src/main.cpp"
  "#" " include <${python_header_name}>\nint main() { return 0; }\n")
string(CONCAT python_package_command "find" "_package")
string(CONCAT python_package_name "Py" "ThOn3")
file(WRITE "${source_python_package_root}/cmake/dependency.cmake"
  "${python_package_command} ( ${python_package_name} REQUIRED COMPONENTS Interpreter )\n")
string(CONCAT python_interpreter_command "Py" "ThOn3")
file(WRITE "${source_python_codegen_root}/cmake/generate.cmake"
  "add_custom_command(OUTPUT generated.cpp COMMAND "
  "${python_interpreter_command} generate.py)\n")
string(CONCAT fetch_module "Fetch" "Content")
string(CONCAT fetch_declare "${fetch_module}" "_" "Declare")
file(WRITE "${source_fetch_root}/cmake/dependency.cmake"
  "include ( ${fetch_module} )\n"
  "${fetch_declare} ( remote SOURCE_DIR remote )\n")
string(CONCAT external_module "External" "Project")
string(CONCAT external_add "${external_module}" "_" "Add")
file(WRITE "${source_external_root}/cmake/dependency.cmake"
  "include ( ${external_module} )\n"
  "${external_add} ( remote URL https://invalid.example/source.tar )\n")
string(CONCAT download_operation "DOWN" "LOAD")
file(WRITE "${source_download_root}/cmake/dependency.cmake"
  "file ( ${download_operation} https://invalid.example/source.tar downloaded.tar )\n")

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

function(hundun_verify_source_policy_negative child_case
         expected_relative_path expected_message)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" "-DHUNDUN_PROVENANCE_CASE=${child_case}"
            -P "${CMAKE_CURRENT_LIST_FILE}"
    RESULT_VARIABLE child_result
    OUTPUT_VARIABLE child_stdout
    ERROR_VARIABLE child_stderr)

  if(child_result STREQUAL "0")
    message(FATAL_ERROR
      "Source-policy child '${child_case}' unexpectedly succeeded")
  endif()

  string(CONCAT child_output "${child_stdout}" "${child_stderr}")
  set(expected_path
    "${fixture_prefix}-${child_case}/${expected_relative_path}")
  foreach(expected_fragment IN ITEMS "${expected_message}" "${expected_path}")
    string(FIND "${child_output}" "${expected_fragment}" fragment_position)
    if(fragment_position EQUAL -1)
      message(FATAL_ERROR
        "Source-policy child '${child_case}' did not emit expected fragment "
        "'${expected_fragment}'. Output:\n${child_output}")
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
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_clean")
  hundun_assert_public_dependency_policy("${source_clean_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_python_file")
  hundun_assert_public_dependency_policy("${source_python_file_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_python_header")
  hundun_assert_public_dependency_policy("${source_python_header_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_python_package")
  hundun_assert_public_dependency_policy("${source_python_package_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_python_codegen")
  hundun_assert_public_dependency_policy("${source_python_codegen_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_fetch")
  hundun_assert_public_dependency_policy("${source_fetch_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_external")
  hundun_assert_public_dependency_policy("${source_external_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_download")
  hundun_assert_public_dependency_policy("${source_download_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_missing")
  hundun_assert_public_dependency_policy("${source_missing_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_python_file")
  hundun_verify_source_policy_negative(
    source_python_file "source-python-file/tools/generate.py"
    "Forbidden public Python source file")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_python_header")
  hundun_verify_source_policy_negative(
    source_python_header "source-python-header/src/main.cpp"
    "Forbidden Python header dependency")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_python_package")
  hundun_verify_source_policy_negative(
    source_python_package "source-python-package/cmake/dependency.cmake"
    "Forbidden Python package discovery")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_python_codegen")
  hundun_verify_source_policy_negative(
    source_python_codegen "source-python-codegen/cmake/generate.cmake"
    "Forbidden Python interpreter dependency")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_fetch")
  hundun_verify_source_policy_negative(
    source_fetch "source-fetch/cmake/dependency.cmake"
    "Forbidden FetchContent source retrieval")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_external")
  hundun_verify_source_policy_negative(
    source_external "source-external/cmake/dependency.cmake"
    "Forbidden ExternalProject source retrieval")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_download")
  hundun_verify_source_policy_negative(
    source_download "source-download/cmake/dependency.cmake"
    "Forbidden file download source retrieval")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_missing")
  hundun_verify_source_policy_negative(
    source_missing "source-missing"
    "Public dependency scan root does not exist")
else()
  message(FATAL_ERROR
    "Unknown HUNDUN_PROVENANCE_CASE: ${HUNDUN_PROVENANCE_CASE}")
endif()
