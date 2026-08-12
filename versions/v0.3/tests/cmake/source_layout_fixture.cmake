# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_LAYOUT_ROOT)
  message(FATAL_ERROR "HUNDUN_LAYOUT_ROOT is required")
endif()

get_filename_component(HUNDUN_LAYOUT_ROOT
  "${HUNDUN_LAYOUT_ROOT}" ABSOLUTE)
set(HUNDUN_PUBLIC_ROOT "${HUNDUN_LAYOUT_ROOT}/include/hundun")
set(HUNDUN_SOURCE_ROOT "${HUNDUN_LAYOUT_ROOT}/src")

if(NOT IS_DIRECTORY "${HUNDUN_PUBLIC_ROOT}")
  message(FATAL_ERROR "flat public include directory is missing")
endif()
if(NOT IS_DIRECTORY "${HUNDUN_SOURCE_ROOT}")
  message(FATAL_ERROR "flat product source directory is missing")
endif()

file(GLOB HUNDUN_PUBLIC_ENTRIES LIST_DIRECTORIES true
  "${HUNDUN_PUBLIC_ROOT}/*")
foreach(HUNDUN_ENTRY IN LISTS HUNDUN_PUBLIC_ENTRIES)
  if(IS_DIRECTORY "${HUNDUN_ENTRY}")
    message(FATAL_ERROR "public include directory is not flat: ${HUNDUN_ENTRY}")
  endif()
endforeach()

file(GLOB HUNDUN_SOURCE_ENTRIES LIST_DIRECTORIES true
  "${HUNDUN_SOURCE_ROOT}/*")
foreach(HUNDUN_ENTRY IN LISTS HUNDUN_SOURCE_ENTRIES)
  if(IS_DIRECTORY "${HUNDUN_ENTRY}")
    message(FATAL_ERROR "product source directory is not flat: ${HUNDUN_ENTRY}")
  endif()
endforeach()

set(HUNDUN_PREFIX_PATTERN
  "^(app|cfg|rt|exec|mesh|bc|ib|lin|fvm|flow|les|diag|sdk|chem)_[a-z0-9_]+\\.(h|hpp|cpp)$")
file(GLOB HUNDUN_PUBLIC_FILES
  "${HUNDUN_PUBLIC_ROOT}/*.h" "${HUNDUN_PUBLIC_ROOT}/*.hpp")
file(GLOB HUNDUN_PRODUCT_FILES
  "${HUNDUN_SOURCE_ROOT}/*.cpp" "${HUNDUN_SOURCE_ROOT}/*.hpp")
foreach(HUNDUN_FILE IN LISTS HUNDUN_PUBLIC_FILES HUNDUN_PRODUCT_FILES)
  get_filename_component(HUNDUN_NAME "${HUNDUN_FILE}" NAME)
  if(NOT HUNDUN_NAME MATCHES "${HUNDUN_PREFIX_PATTERN}")
    message(FATAL_ERROR "unregistered product file prefix: ${HUNDUN_NAME}")
  endif()
endforeach()

foreach(HUNDUN_FILE IN LISTS HUNDUN_PRODUCT_FILES)
  file(READ "${HUNDUN_FILE}" HUNDUN_PRODUCT_TEXT)
  if(HUNDUN_PRODUCT_TEXT MATCHES "tests/")
    message(FATAL_ERROR
      "product source references a tests/ path: ${HUNDUN_FILE}")
  endif()
endforeach()

set(HUNDUN_SOURCE_CMAKE "${HUNDUN_SOURCE_ROOT}/CMakeLists.txt")
if(EXISTS "${HUNDUN_SOURCE_CMAKE}")
  file(READ "${HUNDUN_SOURCE_CMAKE}" HUNDUN_SOURCE_CMAKE_TEXT)
  if(HUNDUN_SOURCE_CMAKE_TEXT MATCHES "tests/")
    message(FATAL_ERROR
      "product build references a tests/ path: ${HUNDUN_SOURCE_CMAKE}")
  endif()
endif()

foreach(HUNDUN_HEADER IN LISTS HUNDUN_PUBLIC_FILES)
  file(READ "${HUNDUN_HEADER}" HUNDUN_HEADER_TEXT)
  if(HUNDUN_HEADER_TEXT MATCHES
      "#include[ \\t]+[<\"][^>\"]*(_detail|_test_access|_test_seam)\\.hpp[>\"]")
    message(FATAL_ERROR "public header includes a private header: ${HUNDUN_HEADER}")
  endif()
endforeach()

file(GLOB_RECURSE HUNDUN_CPP_TEXT_FILES
  "${HUNDUN_LAYOUT_ROOT}/include/*.h"
  "${HUNDUN_LAYOUT_ROOT}/include/*.hpp"
  "${HUNDUN_LAYOUT_ROOT}/src/*.cpp"
  "${HUNDUN_LAYOUT_ROOT}/src/*.hpp"
  "${HUNDUN_LAYOUT_ROOT}/tests/*.cpp"
  "${HUNDUN_LAYOUT_ROOT}/tests/*.hpp")
foreach(HUNDUN_FILE IN LISTS HUNDUN_CPP_TEXT_FILES)
  file(READ "${HUNDUN_FILE}" HUNDUN_TEXT)
  if(HUNDUN_TEXT MATCHES
      "hundun/(application|sdk|config|runtime|execution|mesh|boundary|immersed|linear|finite_volume|flow|diagnostics|solver)/")
    message(FATAL_ERROR "old module include path remains: ${HUNDUN_FILE}")
  endif()
endforeach()
