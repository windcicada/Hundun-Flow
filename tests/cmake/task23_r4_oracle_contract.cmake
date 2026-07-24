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
  if(NOT DEFINED HUNDUN_R6_SECTION OR HUNDUN_R6_SECTION STREQUAL "boundary")
    string(FIND "${checkpoint_test}"
      "if (inlet_patch == 0U && authority_case == 0U)"
      scalar_patch_zero_gate)
    if(NOT scalar_patch_zero_gate EQUAL -1)
      message(FATAL_ERROR
        "Task 23 R6 scalar schema mutations remain guarded to patch 0")
    endif()
    foreach(required IN ITEMS
        "for (std::size_t material_inlet_patch"
        "r6_material_density_all_patches"
        "boundaries[material_inlet_patch]"
        "r5_scalar_list_membership"
        "r5_scalar_name")
      string(FIND "${checkpoint_test}" "${required}" position)
      if(position EQUAL -1)
        message(FATAL_ERROR
          "Task 23 R6 six-patch boundary matrix is missing '${required}'")
      endif()
    endforeach()
  endif()
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
  if(NOT DEFINED HUNDUN_R6_SECTION OR HUNDUN_R6_SECTION STREQUAL "limits")
    string(REGEX MATCHALL "expected_manifest_actual_size\\("
      manifest_formula_calls "${checkpoint_product}")
    list(LENGTH manifest_formula_calls manifest_formula_call_count)
    if(manifest_formula_call_count LESS 4)
      message(FATAL_ERROR
        "Task 23 R6 manifest authentication lacks the shared checked formula")
    endif()
    foreach(required IN ITEMS
        "test_authenticated_rank_wrapper_limits"
        "rank_wrapper_expected_actual_size_mismatch"
        "rank_wrapper_expected_rank_mismatch"
        "rank_wrapper_expected_rank_count_mismatch"
        "rank_wrapper_declared_rank_count"
        "rank_wrapper_declared_payload_size"
        "rank_wrapper_checked_sum_overflow"
        "rank_wrapper_platform_size_rejection"
        "test_authenticated_manifest_limits"
        "checkpoint_v2_authenticate_manifest_limits_for_test("
        "manifest_expected_actual_size_mismatch"
        "manifest_expected_rank_count_mismatch"
        "manifest_expected_global_payload_size_mismatch"
        "manifest_declared_rank_count"
        "manifest_declared_global_payload_size"
        "manifest_declared_record_count"
        "manifest_rank_product_limit"
        "manifest_checked_sum_overflow"
        "manifest_platform_size_rejection")
      string(FIND "${protocol_test}" "${required}" position)
      if(position EQUAL -1)
        message(FATAL_ERROR
          "Task 23 R6 authenticated limit call site is missing '${required}'")
      endif()
    endforeach()
    string(REGEX MATCHALL
      "checkpoint_v2_authenticate_rank_wrapper_for_test\\("
      rank_auth_calls "${protocol_test}")
    list(LENGTH rank_auth_calls rank_auth_call_count)
    string(REGEX MATCHALL
      "checkpoint_v2_authenticate_manifest_for_test\\("
      manifest_auth_calls "${protocol_test}")
    list(LENGTH manifest_auth_calls manifest_auth_call_count)
    string(REGEX MATCHALL
      "checkpoint_v2_authenticate_manifest_limits_for_test\\("
      manifest_limit_auth_calls "${protocol_test}")
    list(LENGTH manifest_limit_auth_calls manifest_limit_auth_call_count)
    if(rank_auth_call_count LESS 2 OR manifest_auth_call_count LESS 1 OR
       manifest_limit_auth_call_count LESS 1)
      message(FATAL_ERROR
        "Task 23 R6 limit sections do not call both product authenticators")
    endif()
  endif()
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
