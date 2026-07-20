# SPDX-License-Identifier: Apache-2.0

file(READ "${HUNDUN_SOURCE_ROOT}/flow/include/hundun/flow/momentum_predictor.hpp"
     task17_header)
file(READ "${HUNDUN_SOURCE_ROOT}/flow/src/momentum_predictor.cpp"
     task17_source)

foreach(forbidden IN ITEMS
    "#include <mpi.h>" "MPI_Comm" "MPI_Request" "CUDA" "cuda/" "HIP_"
    "hip/" "SYCL" "sycl/" "PETSc" "HYPRE" "FlowState" "PisoCoupler")
  string(FIND "${task17_header}${task17_source}" "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task 17 public/product source contains ${forbidden}")
  endif()
endforeach()

foreach(dead_copy IN ITEMS "struct TopologySignature" "make_signature(topology)"
                           "TopologySignature signature")
  string(FIND "${task17_source}" "${dead_copy}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR
      "Task 17 face object retains dead topology copy token: ${dead_copy}")
  endif()
endforeach()
