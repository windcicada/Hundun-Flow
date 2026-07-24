# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "HUNDUN_SOURCE_ROOT is required")
endif()

file(READ "${HUNDUN_SOURCE_ROOT}/tests/mpi/test_checkpoint_v2.cpp"
  checkpoint_test)
file(READ
  "${HUNDUN_SOURCE_ROOT}/tests/mpi/test_checkpoint_v2_diagnostics.cpp"
  diagnostics_test)

if(NOT DEFINED HUNDUN_R4_SECTION OR HUNDUN_R4_SECTION STREQUAL "report")
  foreach(required IN ITEMS "require_exact_product_report")
    string(FIND "${checkpoint_test}" "${required}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R4 checkpoint product oracle is missing '${required}'")
    endif()
  endforeach()
endif()

if(NOT DEFINED HUNDUN_R4_SECTION OR HUNDUN_R4_SECTION STREQUAL "fingerprint")
  foreach(required IN ITEMS
      "require_constructible_fingerprint_failure"
      "face_mass_flux"
      "checkpoint_v2_authenticate_global_payload_for_test"
      "checkpoint_v2_authenticate_rank_payload_for_test")
    string(FIND "${checkpoint_test}" "${required}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R4 checkpoint product oracle is missing '${required}'")
    endif()
  endforeach()
endif()

if(NOT DEFINED HUNDUN_R4_SECTION OR HUNDUN_R4_SECTION STREQUAL "diagnostics")
  foreach(required IN ITEMS
      "require_exact_product_diagnostic_record"
      "product_failure_reports"
      "diagnostic_oracle_is_mutation_sensitive")
    string(FIND "${diagnostics_test}" "${required}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R4 diagnostic product oracle is missing '${required}'")
    endif()
  endforeach()
endif()
