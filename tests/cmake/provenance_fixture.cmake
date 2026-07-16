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
set(clean_root "${fixture_root}/clean")
set(fortran_root "${fixture_root}/fortran")
set(token_root "${fixture_root}/token")

file(REMOVE_RECURSE "${fixture_root}")
file(MAKE_DIRECTORY "${clean_root}/cmake" "${fortran_root}" "${token_root}")
file(WRITE "${clean_root}/main.cpp"
  "// SPDX-License-Identifier: Apache-2.0\nint main() { return 0; }\n")
file(COPY "${CMAKE_CURRENT_LIST_DIR}/../../cmake/HundunProvenanceGuard.cmake"
  DESTINATION "${clean_root}/cmake")
file(WRITE "${fortran_root}/legacy.F90" "program legacy\nend program legacy\n")
string(CONCAT forbidden_fixture_token "CoAsT" "_LeGaCy")
file(WRITE "${token_root}/CMakeLists.txt"
  "# SPDX-License-Identifier: Apache-2.0\n"
  "set(reference_name \"${forbidden_fixture_token}\")\n")

if(HUNDUN_PROVENANCE_CASE STREQUAL "clean")
  hundun_assert_clean_tree("${clean_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "fortran")
  hundun_assert_clean_tree("${fortran_root}")
elseif(HUNDUN_PROVENANCE_CASE STREQUAL "token")
  hundun_assert_clean_tree("${token_root}")
else()
  message(FATAL_ERROR
    "Unknown HUNDUN_PROVENANCE_CASE: ${HUNDUN_PROVENANCE_CASE}")
endif()
