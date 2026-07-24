# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "HUNDUN_SOURCE_ROOT is required")
endif()

file(READ "${HUNDUN_SOURCE_ROOT}/tests/mpi/test_checkpoint_v2.cpp"
  checkpoint_test)
file(READ
  "${HUNDUN_SOURCE_ROOT}/tests/mpi/test_checkpoint_v2_diagnostics.cpp"
  diagnostics_test)
file(READ
  "${HUNDUN_SOURCE_ROOT}/tests/unit/test_checkpoint_v2_protocol.cpp"
  protocol_test)
file(READ
  "${HUNDUN_SOURCE_ROOT}/flow/src/checkpoint_v2.cpp"
  checkpoint_product)

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
      "checkpoint_v2_authenticate_rank_payload_for_test"
      "r5_authoritative_thermal_numeric"
      "r5_density_optional_presence"
      "r5_density_numeric_value"
      "r5_scalar_list_membership"
      "r5_scalar_name"
      "r5_inlet_scalar_value")
    string(FIND "${checkpoint_test}" "${required}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R4 checkpoint product oracle is missing '${required}'")
    endif()
  endforeach()
  foreach(required IN ITEMS
      "checkpoint_v2_authenticate_rank_wrapper_for_test"
      "checkpoint_v2_authenticate_manifest_for_test"
      "checkpoint_v2_authenticate_completed_marker_for_test")
    string(FIND "${protocol_test}" "${required}" test_position)
    string(FIND "${checkpoint_product}" "${required}" product_position)
    if(test_position EQUAL -1 OR product_position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R5 authenticated product call site is missing '${required}'")
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
