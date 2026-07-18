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
set(source_tool_a_root "${fixture_root}/source-tool-a")
set(source_tool_a3_root "${fixture_root}/source-tool-a3")
set(source_tool_b_root "${fixture_root}/source-tool-b")
set(source_vendor_control_root "${fixture_root}/source-vendor-control")
set(source_make_control_root "${fixture_root}/source-make-control")
set(source_vendor_configure_root "${fixture_root}/source-vendor-configure")
set(source_vendor_conditional_root "${fixture_root}/source-vendor-conditional")
set(source_vendor_path_retrieval_root
  "${fixture_root}/source-vendor-path-retrieval")
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
  "${source_download_root}/cmake"
  "${source_tool_a_root}/cmake"
  "${source_tool_a3_root}/cmake"
  "${source_tool_b_root}/cmake"
  "${source_vendor_control_root}/third_party/remote"
  "${source_make_control_root}/control"
  "${source_vendor_configure_root}/third_party/remote"
  "${source_vendor_conditional_root}/third_party/remote"
  "${source_vendor_path_retrieval_root}/third_party/remote")
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
string(CONCAT package_tool_a "p" "ip")
string(CONCAT package_tool_a3 "p" "ip3")
string(CONCAT package_tool_b "con" "da")
file(WRITE "${source_tool_a_root}/cmake/dependency.cmake"
  "execute_process ( COMMAND ${package_tool_a} install local-package )\n")
file(WRITE "${source_tool_a3_root}/cmake/dependency.cmake"
  "add_custom_command ( OUTPUT generated COMMAND /opt/tools/${package_tool_a3} install local-package )\n")
file(WRITE "${source_tool_b_root}/cmake/dependency.cmake"
  "execute_process ( COMMAND /usr/bin/env ${package_tool_b} install local-package )\n")
file(WRITE
  "${source_vendor_control_root}/third_party/remote/dependency.cmake"
  "include ( ${fetch_module} )\n"
  "${fetch_declare} ( remote SOURCE_DIR remote )\n")
file(WRITE "${source_make_control_root}/control/Makefile"
  "all:\n\t${package_tool_a} install local-package\n")
file(WRITE "${source_vendor_configure_root}/third_party/remote/configure"
  "#!/bin/sh\n${package_tool_a} install local-package\n")
file(WRITE "${source_vendor_conditional_root}/third_party/remote/build.sh"
  "#!/bin/sh\n"
  "if ${package_tool_a} install local-package; then\n  :\nfi\n")
string(CONCAT command_line_retrieval_tool "g" "it")
string(CONCAT command_line_retrieval_action "cl" "one")
file(WRITE
  "${source_vendor_path_retrieval_root}/third_party/remote/dependency.cmake"
  "execute_process ( COMMAND /usr/bin/${command_line_retrieval_tool} "
  "${command_line_retrieval_action} https://invalid.example/source.git )\n")

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
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_tool_a")
  hundun_assert_public_dependency_policy("${source_tool_a_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_tool_a3")
  hundun_assert_public_dependency_policy("${source_tool_a3_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_tool_b")
  hundun_assert_public_dependency_policy("${source_tool_b_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_vendor_control")
  hundun_assert_public_dependency_policy("${source_vendor_control_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_make_control")
  hundun_assert_public_dependency_policy("${source_make_control_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_vendor_configure")
  hundun_assert_public_dependency_policy("${source_vendor_configure_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_vendor_conditional")
  hundun_assert_public_dependency_policy("${source_vendor_conditional_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "source_vendor_path_retrieval")
  hundun_assert_public_dependency_policy("${source_vendor_path_retrieval_root}")
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
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_tool_a")
  string(CONCAT expected_tool "p" "ip")
  hundun_verify_source_policy_negative(
    source_tool_a "source-tool-a/cmake/dependency.cmake"
    "Forbidden direct package-manager command '${expected_tool}'")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_tool_a3")
  string(CONCAT expected_tool "p" "ip3")
  hundun_verify_source_policy_negative(
    source_tool_a3 "source-tool-a3/cmake/dependency.cmake"
    "Forbidden direct package-manager command '${expected_tool}'")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_tool_b")
  string(CONCAT expected_tool "con" "da")
  hundun_verify_source_policy_negative(
    source_tool_b "source-tool-b/cmake/dependency.cmake"
    "Forbidden direct package-manager command '${expected_tool}'")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_vendor_control")
  string(CONCAT expected_retrieval "Fetch" "Content")
  hundun_verify_source_policy_negative(
    source_vendor_control
    "source-vendor-control/third_party/remote/dependency.cmake"
    "Forbidden ${expected_retrieval} source retrieval")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_make_control")
  string(CONCAT expected_tool "p" "ip")
  hundun_verify_source_policy_negative(
    source_make_control "source-make-control/control/Makefile"
    "Forbidden direct package-manager command '${expected_tool}'")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_vendor_configure")
  string(CONCAT expected_tool "p" "ip")
  hundun_verify_source_policy_negative(
    source_vendor_configure
    "source-vendor-configure/third_party/remote/configure"
    "Forbidden direct package-manager command '${expected_tool}'")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_vendor_conditional")
  string(CONCAT expected_tool "p" "ip")
  hundun_verify_source_policy_negative(
    source_vendor_conditional
    "source-vendor-conditional/third_party/remote/build.sh"
    "Forbidden direct package-manager command '${expected_tool}'")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_vendor_path_retrieval")
  hundun_verify_source_policy_negative(
    source_vendor_path_retrieval
    "source-vendor-path-retrieval/third_party/remote/dependency.cmake"
    "Forbidden command-line source retrieval")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_source_missing")
  hundun_verify_source_policy_negative(
    source_missing "source-missing"
    "Public dependency scan root does not exist")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "verify_acceptance_complete_gate")
  set(configure_marker "${fixture_root}/configure-attempted")
  set(failing_command "${fixture_root}/known-failing-command.sh")
  file(WRITE "${failing_command}"
    "#!/bin/sh\n"
    ": >\"${configure_marker}\"\n"
    "exit 1\n")
  file(CHMOD "${failing_command}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env
            "HUNDUN_ACCEPTANCE_FAST_FIXTURES_ONLY=1"
            "CMAKE_COMMAND=${failing_command}"
            "CTEST_COMMAND=${failing_command}"
            "MPIEXEC_COMMAND=${failing_command}"
            "LDD_COMMAND=${failing_command}"
            bash
            "${CMAKE_CURRENT_LIST_DIR}/../acceptance/stage1_acceptance.sh"
    RESULT_VARIABLE acceptance_result
    OUTPUT_VARIABLE acceptance_stdout
    ERROR_VARIABLE acceptance_stderr)
  if(NOT acceptance_result STREQUAL "1")
    message(FATAL_ERROR
      "Acceptance complete-gate fixture returned ${acceptance_result}. "
      "stdout:\n${acceptance_stdout}\nstderr:\n${acceptance_stderr}")
  endif()
  if(NOT EXISTS "${configure_marker}")
    message(FATAL_ERROR
      "Acceptance complete-gate fixture did not attempt configure")
  endif()
else()
  message(FATAL_ERROR
    "Unknown HUNDUN_PROVENANCE_CASE: ${HUNDUN_PROVENANCE_CASE}")
endif()
