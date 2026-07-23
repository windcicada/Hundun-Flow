# SPDX-License-Identifier: Apache-2.0
if(NOT DEFINED HUNDUN_SOURCE_ROOT OR NOT DEFINED HUNDUN_TESTS_OFF_FLOW)
  message(FATAL_ERROR "Task22 test-access inputs are required")
endif()
execute_process(
  COMMAND nm -C "${HUNDUN_TESTS_OFF_FLOW}"
  RESULT_VARIABLE nm_result OUTPUT_VARIABLE symbols ERROR_VARIABLE nm_error)
if(NOT nm_result EQUAL 0)
  message(FATAL_ERROR "nm failed: ${nm_error}")
endif()
if(symbols MATCHES "AdaptiveTimeControlTestAccess|reseal|failure.schedule")
  message(FATAL_ERROR "Task22 test-only symbol leaked into tests-off product")
endif()
