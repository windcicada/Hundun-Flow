# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

set(paths
  "${HUNDUN_SOURCE_ROOT}/include/hundun/flow_ideal_gas_closure.hpp"
  "${HUNDUN_SOURCE_ROOT}/include/hundun/flow_ideal_gas_piso.hpp"
  "${HUNDUN_SOURCE_ROOT}/src/flow_ideal_gas_closure.cpp"
  "${HUNDUN_SOURCE_ROOT}/src/flow_ideal_gas_piso.cpp"
  "${HUNDUN_SOURCE_ROOT}/include/hundun/diag_ideal_gas_closure.hpp"
  "${HUNDUN_SOURCE_ROOT}/src/diag_ideal_gas_closure.cpp")
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
set(closure_source "${HUNDUN_SOURCE_ROOT}/src/flow_ideal_gas_closure.cpp")
set(piso_source "${HUNDUN_SOURCE_ROOT}/src/flow_ideal_gas_piso.cpp")
set(diagnostic_source
    "${HUNDUN_SOURCE_ROOT}/src/diag_ideal_gas_closure.cpp")
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
foreach(forbidden IN ITEMS
    "MPI_Send(" "MPI_Recv(" "MPI_TAG_UB"
    "MPI_Comm_get_attr(ideal-gas closure preflight tag)"
    "constexpr int tag = 31741")
  string(FIND "${closure_text}" "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR
      "Task 21 bounded preflight exchange guard found ${forbidden}")
  endif()
endforeach()
foreach(required IN ITEMS
    "MPI_Allgather(" "preflight_workspace"
    "preflight_wire_exchange_count")
  string(FIND "${closure_text}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR
      "Task 21 bounded preflight exchange guard misses ${required}")
  endif()
endforeach()
string(REGEX MATCHALL "MPI_Allgather\\(&local" preflight_exchanges
  "${closure_text}")
list(LENGTH preflight_exchanges preflight_exchange_count)
if(NOT preflight_exchange_count EQUAL 1)
  message(FATAL_ERROR
    "Task 21 preflight must have exactly one bounded wire-exchange implementation")
endif()
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
