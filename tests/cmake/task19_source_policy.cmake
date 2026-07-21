# SPDX-License-Identifier: Apache-2.0

set(core_header
    "${HUNDUN_SOURCE_ROOT}/diagnostics/include/hundun/diagnostics/structured_diagnostics.hpp")
set(core_source
    "${HUNDUN_SOURCE_ROOT}/diagnostics/src/structured_diagnostics.cpp")
set(adapter_source
    "${HUNDUN_SOURCE_ROOT}/diagnostics/src/material_density_transport_diagnostics.cpp")
set(flow_source
    "${HUNDUN_SOURCE_ROOT}/flow/src/material_density_transport.cpp")
file(READ "${core_header}" core_header_text)
file(READ "${core_source}" core_source_text)
file(READ "${adapter_source}" adapter_text)
file(READ "${flow_source}" flow_text)

foreach(forbidden IN ITEMS
    "#include <mpi.h>" "hundun/runtime/" "hundun/flow/"
    "yyjson" "Python.h" "Py_" "CUDA" "HIP_" "PETSc" "HYPRE")
  string(FIND "${core_header_text}${core_source_text}" "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Structured diagnostics core contains ${forbidden}")
  endif()
endforeach()

foreach(ignored_result IN ITEMS
    "static_cast<void>(MPI_" "(void)MPI_" "== MPI_SUCCESS"
    "!= MPI_SUCCESS")
  string(FIND "${adapter_text}${flow_text}" "${ignored_result}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task 19 source bypasses typed result checking")
  endif()
endforeach()

string(FIND "${core_source_text}"
            "std::string to_canonical_json(const DiagnosticRecord &record)"
            serializer_start)
if(serializer_start EQUAL -1)
  message(FATAL_ERROR "Task 19 serializer body was not found")
endif()
string(SUBSTRING "${core_source_text}" ${serializer_start} -1 serializer_text)
foreach(hidden_copy IN ITEMS
    "DiagnosticRecord copy" "DiagnosticRecord normalized" "std::sort("
    "std::stable_sort(" "normalize_record(")
  string(FIND "${serializer_text}" "${hidden_copy}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task 19 serializer hides record normalization")
  endif()
endforeach()

foreach(forbidden IN ITEMS "Task20" "Stage3" "Python.h" "Py_")
  string(FIND "${adapter_text}${flow_text}" "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Task 19 product source contains ${forbidden}")
  endif()
endforeach()
