# SPDX-License-Identifier: Apache-2.0
if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "HUNDUN_SOURCE_ROOT is required")
endif()

file(READ "${HUNDUN_SOURCE_ROOT}/CMakeLists.txt" cmake_text)
file(READ "${HUNDUN_SOURCE_ROOT}/flow/src/adaptive_time_control.cpp" flow_text)
file(READ "${HUNDUN_SOURCE_ROOT}/flow/src/adaptive_time_control_detail.hpp"
     detail_text)
file(READ "${HUNDUN_SOURCE_ROOT}/diagnostics/src/time_control_diagnostics.cpp"
     diagnostics_text)

if(cmake_text MATCHES
   "target_link_libraries\\(hundun_flow[^\\)]*hundun_diagnostics")
  message(FATAL_ERROR "hundun_flow must not link diagnostics")
endif()
if(flow_text MATCHES "hundun/diagnostics|time_control_diagnostics")
  message(FATAL_ERROR "flow time control must not depend on diagnostics")
endif()
if(diagnostics_text MATCHES "acquire_write|trial_layer|mutable.*FieldView")
  message(FATAL_ERROR "time-control diagnostics uses mutable field access")
endif()
if(diagnostics_text MATCHES
   "#include \"time_control_diagnostics_test_access.hpp\"" AND
   NOT diagnostics_text MATCHES
   "#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS[\n\r]+#include \"time_control_diagnostics_test_access.hpp\"")
  message(FATAL_ERROR "diagnostics test access is not compile-gated")
endif()
if(cmake_text MATCHES
   "add_executable\\(test_adaptive_time_control_mpi[^\\)]*test_adaptive_time_control_diagnostics.cpp")
  message(FATAL_ERROR "controller and diagnostics MPI roles are combined")
endif()
foreach(required IN ITEMS
    "add_executable(test_adaptive_time_control_mpi"
    "add_executable(test_adaptive_time_control_diagnostics_mpi"
    "test_task22_tests_off_header_contract"
    "HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS")
  string(FIND "${cmake_text}" "${required}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR "Task22 CMake authority missing: ${required}")
  endif()
endforeach()
foreach(required IN ITEMS
    "struct TimeControlStateCodec final"
    "struct detail::AdaptiveTimeControlEngine final")
  string(FIND "${detail_text}${flow_text}" "${required}" first)
  if(first EQUAL -1)
    message(FATAL_ERROR "Task22 authority missing: ${required}")
  endif()
  math(EXPR next "${first} + 1")
  string(SUBSTRING "${detail_text}${flow_text}" ${next} -1 tail)
  string(FIND "${tail}" "${required}" duplicate)
  if(NOT duplicate EQUAL -1)
    message(FATAL_ERROR "Task22 authority duplicated: ${required}")
  endif()
endforeach()
foreach(required IN ITEMS
    "request_projection("
    "provider_projection("
    "std::array<Tuple, 9>"
    "render_owned_cells(")
  string(FIND "${flow_text}${detail_text}${diagnostics_text}" "${required}"
       position)
  if(position EQUAL -1)
    message(FATAL_ERROR "Task22 product authority missing: ${required}")
  endif()
endforeach()
foreach(forbidden IN ITEMS
    "Python.h" "pybind" "manifest.v2" "COMPLETED"
    "adaptive_time_control_test_support")
  if(flow_text MATCHES "${forbidden}" OR
     diagnostics_text MATCHES "${forbidden}" OR
     detail_text MATCHES "${forbidden}")
    message(FATAL_ERROR "Task22 product contains forbidden term ${forbidden}")
  endif()
endforeach()
