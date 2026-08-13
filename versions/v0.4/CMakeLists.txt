# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.21)

project(HundunFlow VERSION 0.4.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(HUNDUN_BUILD_TESTS "Build HUNDUN-FLOW tests" OFF)
option(HUNDUN_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(HUNDUN_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)

find_package(MPI 3 REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)

get_filename_component(hundun_repository_root
  "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(hundun_v04_yyjson_root "${hundun_repository_root}/third_party/yyjson")
foreach(hundun_v04_yyjson_file IN ITEMS yyjson.c yyjson.h UPSTREAM.json)
  if(NOT EXISTS "${hundun_v04_yyjson_root}/${hundun_v04_yyjson_file}")
    message(FATAL_ERROR
      "v0.4 requires bundled yyjson at the fixed repository path: "
      "${hundun_v04_yyjson_root}/${hundun_v04_yyjson_file}")
  endif()
endforeach()

add_library(hundun_v04_build_options INTERFACE)
target_compile_definitions(hundun_v04_build_options INTERFACE
  OMPI_SKIP_MPICXX=1
  MPICH_SKIP_MPICXX=1)
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
  src/app_case.cpp
  src/core_arena.cpp
  src/core_field_registry.cpp
  src/core_revision.cpp
  src/core_status.cpp
  src/core_transaction.cpp
  src/mesh_cartesian.cpp
  src/mesh_decomposition.cpp
  src/mesh_focus.cpp
  src/mesh_stl_scan.cpp
  "${hundun_v04_yyjson_root}/yyjson.c")
target_include_directories(hundun_v04_core PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(hundun_v04_core PRIVATE
  "${hundun_v04_yyjson_root}")
target_link_libraries(hundun_v04_core PUBLIC
  hundun_v04_build_options
  MPI::MPI_C)
target_link_libraries(hundun_v04_core PRIVATE
  Threads::Threads)
if(HUNDUN_BUILD_TESTS)
  target_compile_definitions(hundun_v04_core PRIVATE
    HUNDUN_V04_ENABLE_TEST_ACCESS=1)
endif()

add_executable(hundun
  src/app_main.cpp)
target_link_libraries(hundun PRIVATE hundun_v04_core)

if(HUNDUN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
