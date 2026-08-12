# SPDX-License-Identifier: Apache-2.0

set(core_header
    "${HUNDUN_SOURCE_ROOT}/include/hundun/diag_structured.hpp")
set(core_source
    "${HUNDUN_SOURCE_ROOT}/src/diag_structured.cpp")
set(adapter_source
    "${HUNDUN_SOURCE_ROOT}/src/diag_material_density_transport.cpp")
set(adapter_test_access
    "${HUNDUN_SOURCE_ROOT}/tests/support/diag_material_density_transport_test_access.hpp")
set(flow_source
    "${HUNDUN_SOURCE_ROOT}/src/flow_material_density_transport.cpp")
file(READ "${core_header}" core_header_text)
file(READ "${core_source}" core_source_text)
file(READ "${adapter_source}" adapter_text)
file(READ "${adapter_test_access}" adapter_test_access_text)
file(READ "${flow_source}" flow_text)

foreach(forbidden IN ITEMS
    "#include <mpi.h>" "hundun/rt_" "hundun/flow_"
    "yyjson" "Python.h" "Py_" "CUDA" "HIP_" "PETSc" "HYPRE")
  string(FIND "${core_header_text}${core_source_text}" "${forbidden}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR "Structured diagnostics core contains ${forbidden}")
  endif()
endforeach()

foreach(arbitrary_count_seam IN ITEMS
    "override_reported_transport_total_count"
    "override_reported_sample_wire_bytes"
    "reported_transport_total_count"
    "reported_sample_wire_bytes"
    "reported_sample_wire_byte_count")
  string(FIND "${adapter_text}${adapter_test_access_text}"
              "${arbitrary_count_seam}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR
            "Task 19 retains arbitrary count seam ${arbitrary_count_seam}")
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

foreach(native_wire IN ITEMS
    "struct SampleWire" "sizeof(SampleWire)" "std::memcpy"
    "reinterpret_cast<const unsigned char" "global_max_id"
    "material diagnostic ownership IDs")
  string(FIND "${adapter_text}" "${native_wire}" position)
  if(NOT position EQUAL -1)
    message(FATAL_ERROR
            "Task 19 adapter retains forbidden native/full-ID path ${native_wire}")
  endif()
endforeach()

foreach(required_portable_path IN ITEMS
    "kRequestWireSchemaV1" "kSampleWireSchemaV1"
    "decode_request_wire" "decode_sample_wire"
    "MPI_Allgather(material diagnostic ownership boxes)"
    "material.diagnostics.transport-total-count")
  string(FIND "${adapter_text}" "${required_portable_path}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR
            "Task 19 adapter misses portable path ${required_portable_path}")
  endif()
endforeach()
