# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "HUNDUN_SOURCE_ROOT is required")
endif()

file(READ
  "${HUNDUN_SOURCE_ROOT}/flow/src/checkpoint_v2_test_access.hpp"
  test_access)
string(FIND "${test_access}" "#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS" guard)
if(guard EQUAL -1)
  message(FATAL_ERROR "Checkpoint v2 test access lacks tests-off guard")
endif()
foreach(required IN ITEMS
    "CheckpointV2DeepSnapshot"
    "checkpoint_v2_deep_snapshot_equal"
    "checkpoint_v2_failed_read_preserved_values"
    "set_committed_density_ghost"
    "set_rollback_density_ghost"
    "force_generation")
  string(FIND "${test_access}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR
      "Checkpoint v2 deep transaction test access misses ${required}")
  endif()
endforeach()

foreach(public_header IN ITEMS
    "${HUNDUN_SOURCE_ROOT}/flow/include/hundun/flow/checkpoint_v2.hpp"
    "${HUNDUN_SOURCE_ROOT}/diagnostics/include/hundun/diagnostics/checkpoint_v2_diagnostics.hpp")
  file(READ "${public_header}" content)
  string(FIND "${content}" "checkpoint_v2_test_access" leak)
  if(NOT leak EQUAL -1)
    message(FATAL_ERROR "Task 23 test access leaked into ${public_header}")
  endif()
endforeach()

foreach(product_source IN ITEMS
    "${HUNDUN_SOURCE_ROOT}/flow/src/checkpoint_v2.cpp"
    "${HUNDUN_SOURCE_ROOT}/flow/src/flow_state.cpp")
  file(READ "${product_source}" content)
  string(FIND "${content}" "CheckpointV2TestAccess" leak)
  if(NOT leak EQUAL -1)
    message(FATAL_ERROR
      "Checkpoint v2 test-only access leaked into product source ${product_source}")
  endif()
endforeach()

file(READ "${HUNDUN_SOURCE_ROOT}/CMakeLists.txt" cmake_text)
string(FIND "${cmake_text}" "hundun_checkpoint_diagnostics" target)
if(target EQUAL -1)
  message(FATAL_ERROR "Task 23 diagnostics target is missing")
endif()
