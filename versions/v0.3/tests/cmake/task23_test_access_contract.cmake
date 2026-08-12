# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "HUNDUN_SOURCE_ROOT is required")
endif()

file(READ
  "${HUNDUN_SOURCE_ROOT}/tests/support/flow_checkpoint_v2_test_access.hpp"
  test_access)
string(FIND "${test_access}" "#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS" guard)
if(guard EQUAL -1)
  message(FATAL_ERROR "Checkpoint v2 test access lacks tests-off guard")
endif()
foreach(required IN ITEMS
    "CheckpointV2DeepSnapshot"
    "checkpoint_v2_deep_snapshot_equal"
    "checkpoint_v2_failed_read_preserved_values"
    "CheckpointV2PathCodeObservation"
    "checkpoint_v2_path_code_observation_for_test"
    "set_committed_density_ghost"
    "set_rollback_density_ghost"
    "force_generation")
  string(FIND "${test_access}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR
      "Checkpoint v2 deep transaction test access misses ${required}")
  endif()
endforeach()

file(READ
  "${HUNDUN_SOURCE_ROOT}/tests/support/flow_ideal_gas_closure_test_access.hpp"
  ideal_gas_test_access)
string(FIND "${ideal_gas_test_access}"
  "#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS" ideal_gas_guard)
if(ideal_gas_guard EQUAL -1)
  message(FATAL_ERROR "Ideal-gas closure test access lacks tests-off guard")
endif()
foreach(required IN ITEMS
    "set_restore_preparation_fault"
    "set_restore_snapshot_shape_fault"
    "preflight_failure_rank"
    "create_validation_failure_reason")
  string(FIND "${ideal_gas_test_access}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR
      "Ideal-gas restore test access misses ${required}")
  endif()
endforeach()

foreach(public_header IN ITEMS
    "${HUNDUN_SOURCE_ROOT}/include/hundun/flow_checkpoint_v2.hpp"
    "${HUNDUN_SOURCE_ROOT}/include/hundun/flow_ideal_gas_closure.hpp"
    "${HUNDUN_SOURCE_ROOT}/include/hundun/diag_checkpoint_v2.hpp")
  file(READ "${public_header}" content)
  foreach(private_header IN ITEMS
      "checkpoint_v2_test_access"
      "ideal_gas_closure_test_access")
    string(FIND "${content}" "${private_header}" leak)
    if(NOT leak EQUAL -1)
      message(FATAL_ERROR "Task 23 test access leaked into ${public_header}")
    endif()
  endforeach()
endforeach()

foreach(product_source IN ITEMS
    "${HUNDUN_SOURCE_ROOT}/src/flow_checkpoint_v2.cpp"
    "${HUNDUN_SOURCE_ROOT}/src/flow_state.cpp")
  file(READ "${product_source}" content)
  string(FIND "${content}" "CheckpointV2TestAccess" leak)
  if(NOT leak EQUAL -1)
    message(FATAL_ERROR
      "Checkpoint v2 test-only access leaked into product source ${product_source}")
  endif()
endforeach()

file(READ "${HUNDUN_SOURCE_ROOT}/src/CMakeLists.txt" cmake_text)
string(FIND "${cmake_text}" "hundun_checkpoint_diagnostics" target)
if(target EQUAL -1)
  message(FATAL_ERROR "Task 23 diagnostics target is missing")
endif()
