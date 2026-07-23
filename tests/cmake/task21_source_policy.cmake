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

# Repair R2 authority guards.  These source-level checks are deliberately
# narrow: the MPI/product tests carry the behavioural evidence, while this
# fast RED catches regression to the exact rejected implementation patterns.
set(closure_source "${HUNDUN_SOURCE_ROOT}/flow/src/ideal_gas_closure.cpp")
set(piso_source "${HUNDUN_SOURCE_ROOT}/flow/src/ideal_gas_piso.cpp")
set(diagnostic_source
    "${HUNDUN_SOURCE_ROOT}/diagnostics/src/ideal_gas_closure_diagnostics.cpp")
file(READ "${closure_source}" closure_text)
file(READ "${piso_source}" piso_text)
file(READ "${diagnostic_source}" diagnostic_text)
foreach(forbidden IN ITEMS
    "const double stored_h = q_eos[cell] / rho_eos[cell]"
    "impl_->topology->local_face_count()};")
  string(FIND "${closure_text}" "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task 21 Repair R2 closure authority guard found ${forbidden}")
  endif()
endforeach()
string(FIND "${piso_text}" "const_cast<FlowState &>(state)" position)
if(NOT position EQUAL -1)
  message(FATAL_ERROR "Task 21 diagnostics must use const FlowState acquisition")
endif()
foreach(forbidden IN ITEMS
    "closure.diagnostics.communicator"
    "closure.diagnostics.ownership"
    "closure.diagnostics.collective-operation"
    "std::vector<double> rho;"
    "std::vector<double> rho_h;")
  string(FIND "${diagnostic_text}" "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task 21 Repair R2 diagnostics guard found ${forbidden}")
  endif()
endforeach()
