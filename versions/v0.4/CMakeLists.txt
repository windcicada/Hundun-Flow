# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

project(HundunFlow VERSION 0.4.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(HUNDUN_BUILD_TESTS "Build HUNDUN-FLOW tests" OFF)

add_library(hundun_v04_core STATIC
  src/core_status.cpp)
target_include_directories(hundun_v04_core PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_executable(hundun
  src/app_main.cpp)
target_link_libraries(hundun PRIVATE hundun_v04_core)

if(HUNDUN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
