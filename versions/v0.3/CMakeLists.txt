# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

cmake_minimum_required(VERSION 3.21)

set(HUNDUN_VERSION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")
if(NOT EXISTS "${HUNDUN_VERSION_FILE}")
  message(FATAL_ERROR
    "HUNDUN-FLOW VERSION file is missing: ${HUNDUN_VERSION_FILE}")
endif()
file(READ "${HUNDUN_VERSION_FILE}" HUNDUN_PROJECT_VERSION)
string(STRIP "${HUNDUN_PROJECT_VERSION}" HUNDUN_PROJECT_VERSION)
if(HUNDUN_PROJECT_VERSION STREQUAL "")
  message(FATAL_ERROR
    "HUNDUN-FLOW VERSION file is empty: ${HUNDUN_VERSION_FILE}")
endif()

project(HundunFlow VERSION "${HUNDUN_PROJECT_VERSION}" LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(HUNDUN_BUILD_TESTS "Build HUNDUN-FLOW tests" OFF)
option(HUNDUN_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(HUNDUN_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(HUNDUN_ENABLE_CANTERA "Enable the verified Cantera backend" OFF)
set(HUNDUN_CANTERA_PACKAGE_ROOT "" CACHE PATH
  "Root of the verified prebuilt Cantera package")
set(HUNDUN_STAGE4_CANTERA_TEST_MECHANISM "" CACHE FILEPATH
  "External hash-verified synthetic mechanism used by Stage 4 tests")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(HundunProvenanceGuard)
include(HundunSanitizers)
include(HundunWarnings)

if(HUNDUN_ENABLE_CANTERA)
  include(HundunCanteraPackage)
  hundun_configure_cantera_package("${HUNDUN_CANTERA_PACKAGE_ROOT}")
  include(HundunRelocatablePackage)
  hundun_configure_relocatable_package("${HUNDUN_CANTERA_PACKAGE_ROOT}")
endif()

hundun_assert_clean_tree("${PROJECT_SOURCE_DIR}")
hundun_assert_public_dependency_policy("${PROJECT_SOURCE_DIR}")

find_package(MPI 3 REQUIRED COMPONENTS CXX)
find_package(Threads REQUIRED)

add_library(hundun_options INTERFACE)

add_subdirectory(src)

if(HUNDUN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
