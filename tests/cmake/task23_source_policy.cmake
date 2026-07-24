# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "HUNDUN_SOURCE_ROOT is required")
endif()

set(task23_product
  "${HUNDUN_SOURCE_ROOT}/runtime/src/checkpoint_v2_protocol.hpp"
  "${HUNDUN_SOURCE_ROOT}/runtime/src/checkpoint_v2_protocol.cpp"
  "${HUNDUN_SOURCE_ROOT}/flow/include/hundun/flow/checkpoint_v2.hpp"
  "${HUNDUN_SOURCE_ROOT}/flow/include/hundun/flow/flow_state.hpp"
  "${HUNDUN_SOURCE_ROOT}/flow/include/hundun/flow/ideal_gas_closure.hpp"
  "${HUNDUN_SOURCE_ROOT}/flow/src/checkpoint_v2.cpp"
  "${HUNDUN_SOURCE_ROOT}/flow/src/checkpoint_v2_detail.hpp"
  "${HUNDUN_SOURCE_ROOT}/flow/src/flow_state.cpp"
  "${HUNDUN_SOURCE_ROOT}/flow/src/ideal_gas_closure.cpp"
  "${HUNDUN_SOURCE_ROOT}/runtime/include/hundun/runtime/field_storage.hpp"
  "${HUNDUN_SOURCE_ROOT}/runtime/src/field_storage.cpp"
  "${HUNDUN_SOURCE_ROOT}/diagnostics/include/hundun/diagnostics/checkpoint_v2_diagnostics.hpp"
  "${HUNDUN_SOURCE_ROOT}/diagnostics/src/checkpoint_v2_diagnostics.cpp")

foreach(path IN LISTS task23_product)
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "missing Task 23 product file: ${path}")
  endif()
  file(READ "${path}" content)
  foreach(forbidden IN ITEMS
      "Python.h" "pybind" "cuda" "hip/" "HDF5" "MPI_File"
      "restart_binary")
    string(FIND "${content}" "${forbidden}" position)
    if(NOT position EQUAL -1)
      message(FATAL_ERROR
        "Task 23 product contains forbidden dependency/reference '${forbidden}': ${path}")
    endif()
  endforeach()
endforeach()

file(READ
  "${HUNDUN_SOURCE_ROOT}/runtime/src/checkpoint_v2_protocol.cpp"
  protocol)
foreach(required IN ITEMS
    "0x42F0E1EBA9EA3693" "kRankMagic" "kManifestMagic"
    "kCompletedMagic")
  string(FIND "${protocol}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR "Checkpoint v2 protocol misses ${required}")
  endif()
endforeach()
