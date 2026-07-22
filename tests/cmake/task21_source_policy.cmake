# SPDX-License-Identifier: Apache-2.0

set(paths
  "${HUNDUN_SOURCE_ROOT}/flow/include/hundun/flow/ideal_gas_closure.hpp"
  "${HUNDUN_SOURCE_ROOT}/flow/include/hundun/flow/ideal_gas_piso.hpp"
  "${HUNDUN_SOURCE_ROOT}/flow/src/ideal_gas_closure.cpp"
  "${HUNDUN_SOURCE_ROOT}/flow/src/ideal_gas_piso.cpp"
  "${HUNDUN_SOURCE_ROOT}/diagnostics/include/hundun/diagnostics/ideal_gas_closure_diagnostics.hpp"
  "${HUNDUN_SOURCE_ROOT}/diagnostics/src/ideal_gas_closure_diagnostics.cpp")
set(text "")
foreach(path IN LISTS paths)
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "Task 21 product file is missing: ${path}")
  endif()
  file(READ "${path}" part)
  string(APPEND text "${part}")
endforeach()
foreach(forbidden IN ITEMS
    "Python.h" "Py_" "CUDA" "HIP_" "PETSc" "HYPRE"
    "Checkpoint" "DiagnosticSession" "SourceCallback")
  string(FIND "${text}" "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task 21 product source contains ${forbidden}")
  endif()
endforeach()
foreach(required IN ITEMS
    "IdealGasClosure" "FixedStepIdealGasFlow"
    "IdealGasClosureDiagnosticSource" "flow.ideal-gas-closure")
  string(FIND "${text}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR "Task 21 product path misses ${required}")
  endif()
endforeach()
