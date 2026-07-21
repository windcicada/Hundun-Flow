# SPDX-License-Identifier: Apache-2.0

set(task18_files
    "${HUNDUN_SOURCE_ROOT}/flow/include/hundun/flow/constant_density_piso.hpp"
    "${HUNDUN_SOURCE_ROOT}/flow/include/hundun/flow/flow_state.hpp"
    "${HUNDUN_SOURCE_ROOT}/flow/src/constant_density_piso.cpp"
    "${HUNDUN_SOURCE_ROOT}/flow/src/flow_state.cpp")

set(task18_text "")
foreach(task18_file IN LISTS task18_files)
  file(READ "${task18_file}" contents)
  string(APPEND task18_text "${contents}")
endforeach()

foreach(forbidden IN ITEMS
    "#include <mpi.h>" "MPI_Comm" "MPI_Request" "MPI_Allreduce"
    "CUDA" "cuda/" "HIP_" "hip/" "SYCL" "sycl/" "PETSc" "HYPRE"
    "Py_" "Python.h" "density closure"
    "retry_count" "max_retries" "third_corrector")
  string(FIND "${task18_text}" "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task 18 public/product source contains ${forbidden}")
  endif()
endforeach()

string(REGEX MATCHALL "coupler\\.correct\\(" corrector_calls
       "${task18_text}")
list(LENGTH corrector_calls corrector_call_count)
if(NOT corrector_call_count EQUAL 2)
  message(FATAL_ERROR
    "Task 18 fixed-step transaction must contain exactly two corrector calls")
endif()

foreach(legacy_view IN ITEMS ".view(" "->view(")
  string(FIND "${task18_text}" "FieldStorage${legacy_view}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR
      "Task 18 public/product source contains legacy FieldStorage view")
  endif()
endforeach()
