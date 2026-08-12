# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

project(HundunFlow VERSION 0.4.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(HUNDUN_BUILD_TESTS "Build HUNDUN-FLOW tests" OFF)
option(HUNDUN_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(HUNDUN_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)

add_library(hundun_v04_build_options INTERFACE)
if(HUNDUN_ENABLE_ASAN OR HUNDUN_ENABLE_UBSAN)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    message(FATAL_ERROR
      "v0.4 sanitizers require a GNU or Clang C++ compiler; found "
      "${CMAKE_CXX_COMPILER_ID}")
  endif()

  set(hundun_v04_sanitizer_flags "")
  if(HUNDUN_ENABLE_ASAN)
    list(APPEND hundun_v04_sanitizer_flags -fsanitize=address)
  endif()
  if(HUNDUN_ENABLE_UBSAN)
    list(APPEND hundun_v04_sanitizer_flags -fsanitize=undefined)
  endif()
  target_compile_options(hundun_v04_build_options INTERFACE
    ${hundun_v04_sanitizer_flags})
  target_link_options(hundun_v04_build_options INTERFACE
    ${hundun_v04_sanitizer_flags})
endif()

add_library(hundun_v04_core STATIC
  src/core_status.cpp)
target_include_directories(hundun_v04_core PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(hundun_v04_core PUBLIC hundun_v04_build_options)

add_executable(hundun
  src/app_main.cpp)
target_link_libraries(hundun PRIVATE hundun_v04_core)

if(HUNDUN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
