# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

cmake_minimum_required(VERSION 3.21)
get_filename_component(HUNDUN_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(CMAKE_CURRENT_SOURCE_DIR "${HUNDUN_ROOT}")
include("${HUNDUN_ROOT}/cmake/HundunVersionDispatch.cmake")

# Exercise the real selector without configuring compilers or starting MPI.
# Full configure/sanitizer coverage remains in version_dispatch_test.cmake.
function(add_subdirectory source binary)
  if(NOT source STREQUAL "${HUNDUN_ROOT}/versions/v0.4" OR
     NOT binary STREQUAL "versions/v0.4")
    message(FATAL_ERROR "selector changed the current source/build location")
  endif()
  set_property(GLOBAL PROPERTY HUNDUN_CURRENT_SELECTED TRUE)
endfunction()

if(DEFINED REJECT_VERSION)
  hundun_add_selected_source_line("${REJECT_VERSION}")
  message(FATAL_ERROR "retired/unknown source was accepted")
endif()

hundun_add_selected_source_line(v0.4)
get_property(selected GLOBAL PROPERTY HUNDUN_CURRENT_SELECTED)
if(NOT selected)
  message(FATAL_ERROR "current source was not selected")
endif()
foreach(version IN ITEMS v0.3 v0.2 v9)
  execute_process(COMMAND "${CMAKE_COMMAND}" "-DREJECT_VERSION=${version}"
    -P "${CMAKE_CURRENT_LIST_FILE}"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(result EQUAL 0 OR NOT stderr MATCHES "unsupported HUNDUN_SOURCE_VERSION")
    message(FATAL_ERROR "${version}: missing explicit rejection: ${stdout}${stderr}")
  endif()
endforeach()

foreach(retired IN ITEMS versions/v0.3 src include)
  if(EXISTS "${HUNDUN_ROOT}/${retired}")
    message(FATAL_ERROR "retired implementation remains: ${retired}")
  endif()
endforeach()
foreach(current IN ITEMS versions/v0.4/CMakeLists.txt
    versions/v0.4/src/core_product_freeze.cpp
    versions/v0.4/include/hundun/v04_product.hpp
    versions/v0.4/tests/CMakeLists.txt tools/v04_thin_domain_runner.cpp
    third_party/yyjson/yyjson.c third_party/yyjson/yyjson.h)
  if(NOT EXISTS "${HUNDUN_ROOT}/${current}")
    message(FATAL_ERROR "current input missing: ${current}")
  endif()
endforeach()
foreach(name IN ITEMS case.json thermophysics.d)
  file(SHA256 "${HUNDUN_ROOT}/examples/minimal/${name}" example)
  if(name STREQUAL "case.json")
    set(fixture "case_minimal_valid.json")
  else()
    set(fixture "${name}")
  endif()
  file(SHA256 "${HUNDUN_ROOT}/versions/v0.4/tests/data/${fixture}" tested)
  if(NOT example STREQUAL tested)
    message(FATAL_ERROR "minimal example differs from the tested fixture: ${name}")
  endif()
endforeach()
message(STATUS "Current-source routing, retired-source rejection and fixture checks passed")
