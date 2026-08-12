# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

project(HundunFlowVersions LANGUAGES NONE)

set(HUNDUN_SOURCE_VERSION "v0.4" CACHE STRING
  "HUNDUN-FLOW source line to configure")
set_property(CACHE HUNDUN_SOURCE_VERSION PROPERTY STRINGS v0.4 v0.3)
option(HUNDUN_BUILD_TESTS "Build HUNDUN-FLOW tests" OFF)

if(HUNDUN_BUILD_TESTS)
  enable_testing()
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(HundunVersionDispatch)
hundun_add_selected_source_line("${HUNDUN_SOURCE_VERSION}")
