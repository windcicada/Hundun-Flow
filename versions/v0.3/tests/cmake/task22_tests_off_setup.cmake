# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

if(NOT DEFINED HUNDUN_SOURCE_ROOT OR
   NOT DEFINED HUNDUN_OUTER_BINARY_DIR OR
   NOT DEFINED HUNDUN_NESTED_BINARY_DIR OR
   NOT DEFINED HUNDUN_ARCHIVE_DIR OR
   NOT DEFINED HUNDUN_NESTED_CONFIG)
  message(FATAL_ERROR "Task22 tests-off setup inputs are incomplete")
endif()

foreach(path_variable IN ITEMS
    HUNDUN_SOURCE_ROOT
    HUNDUN_OUTER_BINARY_DIR
    HUNDUN_NESTED_BINARY_DIR
    HUNDUN_ARCHIVE_DIR)
  if("${${path_variable}}" STREQUAL "")
    message(FATAL_ERROR
      "Task22 tests-off setup path '${path_variable}' is empty")
  endif()
  set(path_value "${${path_variable}}")
  cmake_path(ABSOLUTE_PATH path_value NORMALIZE
             OUTPUT_VARIABLE "${path_variable}_ABSOLUTE")
endforeach()

if(NOT IS_DIRECTORY "${HUNDUN_SOURCE_ROOT_ABSOLUTE}")
  message(FATAL_ERROR "Task22 tests-off source root is not a directory")
endif()
if(NOT IS_DIRECTORY "${HUNDUN_OUTER_BINARY_DIR_ABSOLUTE}")
  message(FATAL_ERROR "Task22 tests-off outer binary root is not a directory")
endif()

file(REAL_PATH "${HUNDUN_SOURCE_ROOT_ABSOLUTE}" task22_source_root)
file(REAL_PATH "${HUNDUN_OUTER_BINARY_DIR_ABSOLUTE}" task22_outer_binary_dir)
cmake_path(GET task22_outer_binary_dir ROOT_PATH task22_filesystem_root)
if(task22_outer_binary_dir STREQUAL task22_filesystem_root)
  message(FATAL_ERROR
    "Task22 tests-off outer binary directory cannot be a filesystem root")
endif()
if(task22_outer_binary_dir STREQUAL task22_source_root)
  message(FATAL_ERROR
    "Task22 tests-off outer binary directory cannot equal the source root")
endif()

if(IS_SYMLINK "${HUNDUN_NESTED_BINARY_DIR_ABSOLUTE}")
  message(FATAL_ERROR
    "Task22 tests-off nested binary directory cannot be a symlink")
endif()
if(EXISTS "${HUNDUN_NESTED_BINARY_DIR_ABSOLUTE}")
  file(REAL_PATH "${HUNDUN_NESTED_BINARY_DIR_ABSOLUTE}"
       task22_nested_binary_dir)
else()
  set(task22_nested_binary_dir "${HUNDUN_NESTED_BINARY_DIR_ABSOLUTE}")
endif()
set(task22_expected_nested "${task22_outer_binary_dir}/task22-tests-off")
cmake_path(NORMAL_PATH task22_expected_nested)
if(NOT task22_nested_binary_dir STREQUAL task22_expected_nested)
  message(FATAL_ERROR
    "Task22 tests-off nested binary directory must be the dedicated "
    "outer/task22-tests-off path")
endif()
if(task22_nested_binary_dir STREQUAL task22_source_root OR
   task22_nested_binary_dir STREQUAL task22_outer_binary_dir OR
   task22_nested_binary_dir STREQUAL task22_filesystem_root)
  message(FATAL_ERROR
    "Task22 tests-off nested binary directory has an unsafe authority")
endif()

if(IS_SYMLINK "${HUNDUN_ARCHIVE_DIR_ABSOLUTE}")
  message(FATAL_ERROR
    "Task22 tests-off archive directory cannot be a symlink")
endif()
if(EXISTS "${HUNDUN_ARCHIVE_DIR_ABSOLUTE}")
  file(REAL_PATH "${HUNDUN_ARCHIVE_DIR_ABSOLUTE}" task22_archive_dir)
else()
  set(task22_archive_dir "${HUNDUN_ARCHIVE_DIR_ABSOLUTE}")
endif()
set(task22_expected_archive "${task22_nested_binary_dir}/archives")
cmake_path(NORMAL_PATH task22_expected_archive)
if(NOT task22_archive_dir STREQUAL task22_expected_archive OR
   task22_archive_dir STREQUAL task22_nested_binary_dir)
  message(FATAL_ERROR
    "Task22 tests-off archive directory must be nested/archives")
endif()

set(HUNDUN_SOURCE_ROOT "${task22_source_root}")
set(HUNDUN_OUTER_BINARY_DIR "${task22_outer_binary_dir}")
set(HUNDUN_NESTED_BINARY_DIR "${task22_nested_binary_dir}")
set(HUNDUN_ARCHIVE_DIR "${task22_archive_dir}")

if(DEFINED HUNDUN_TASK22_PATH_CONTRACT_VALIDATE_ONLY AND
   HUNDUN_TASK22_PATH_CONTRACT_VALIDATE_ONLY)
  return()
endif()

load_cache(
  "${HUNDUN_OUTER_BINARY_DIR}"
  READ_WITH_PREFIX outer_
  CMAKE_GENERATOR
  CMAKE_GENERATOR_INSTANCE
  CMAKE_GENERATOR_PLATFORM
  CMAKE_GENERATOR_TOOLSET
  CMAKE_MAKE_PROGRAM
  CMAKE_BUILD_TYPE
  CMAKE_TOOLCHAIN_FILE
  CMAKE_SYSROOT
  CMAKE_C_COMPILER
  CMAKE_C_COMPILER_TARGET
  CMAKE_C_FLAGS
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_RELWITHDEBINFO
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_CXX_COMPILER
  CMAKE_CXX_COMPILER_TARGET
  CMAKE_CXX_FLAGS
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_CXX_FLAGS_RELWITHDEBINFO
  CMAKE_CXX_FLAGS_MINSIZEREL
  CMAKE_EXE_LINKER_FLAGS
  CMAKE_EXE_LINKER_FLAGS_DEBUG
  CMAKE_EXE_LINKER_FLAGS_RELEASE
  CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO
  CMAKE_EXE_LINKER_FLAGS_MINSIZEREL
  CMAKE_MODULE_LINKER_FLAGS
  CMAKE_MODULE_LINKER_FLAGS_DEBUG
  CMAKE_MODULE_LINKER_FLAGS_RELEASE
  CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO
  CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL
  CMAKE_SHARED_LINKER_FLAGS
  CMAKE_SHARED_LINKER_FLAGS_DEBUG
  CMAKE_SHARED_LINKER_FLAGS_RELEASE
  CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO
  CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL
  CMAKE_STATIC_LINKER_FLAGS
  CMAKE_STATIC_LINKER_FLAGS_DEBUG
  CMAKE_STATIC_LINKER_FLAGS_RELEASE
  CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO
  CMAKE_STATIC_LINKER_FLAGS_MINSIZEREL
  MPI_CXX_COMPILER
  MPI_CXX_COMPILER_FLAGS
  MPI_CXX_SKIP_MPICXX
  MPIEXEC_EXECUTABLE
  HUNDUN_ENABLE_ASAN
  HUNDUN_ENABLE_UBSAN)

if(NOT outer_CMAKE_GENERATOR)
  message(FATAL_ERROR "Task22 tests-off setup cannot recover the outer generator")
endif()

set(configure_command
    "${CMAKE_COMMAND}"
    -S "${HUNDUN_SOURCE_ROOT}"
    -B "${HUNDUN_NESTED_BINARY_DIR}"
    -G "${outer_CMAKE_GENERATOR}"
    "-DHUNDUN_BUILD_TESTS:BOOL=OFF"
    "-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${HUNDUN_ARCHIVE_DIR}"
    "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${HUNDUN_ARCHIVE_DIR}")

if(outer_CMAKE_GENERATOR_PLATFORM)
  list(APPEND configure_command -A "${outer_CMAKE_GENERATOR_PLATFORM}")
endif()
if(outer_CMAKE_GENERATOR_TOOLSET)
  list(APPEND configure_command -T "${outer_CMAKE_GENERATOR_TOOLSET}")
endif()

function(task22_forward_cache variable type)
  if(DEFINED outer_${variable} AND NOT "${outer_${variable}}" STREQUAL "")
    list(APPEND configure_command
         "-D${variable}:${type}=${outer_${variable}}")
    set(configure_command "${configure_command}" PARENT_SCOPE)
  endif()
endfunction()

foreach(variable IN ITEMS
    CMAKE_C_FLAGS
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_C_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_EXE_LINKER_FLAGS
    CMAKE_EXE_LINKER_FLAGS_DEBUG
    CMAKE_EXE_LINKER_FLAGS_RELEASE
    CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO
    CMAKE_EXE_LINKER_FLAGS_MINSIZEREL
    CMAKE_MODULE_LINKER_FLAGS
    CMAKE_MODULE_LINKER_FLAGS_DEBUG
    CMAKE_MODULE_LINKER_FLAGS_RELEASE
    CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO
    CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL
    CMAKE_SHARED_LINKER_FLAGS
    CMAKE_SHARED_LINKER_FLAGS_DEBUG
    CMAKE_SHARED_LINKER_FLAGS_RELEASE
    CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO
    CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL
    CMAKE_STATIC_LINKER_FLAGS
    CMAKE_STATIC_LINKER_FLAGS_DEBUG
    CMAKE_STATIC_LINKER_FLAGS_RELEASE
    CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO
    CMAKE_STATIC_LINKER_FLAGS_MINSIZEREL
    MPI_CXX_COMPILER_FLAGS)
  task22_forward_cache("${variable}" STRING)
endforeach()

foreach(variable IN ITEMS
    CMAKE_C_COMPILER
    CMAKE_CXX_COMPILER
    CMAKE_MAKE_PROGRAM
    CMAKE_TOOLCHAIN_FILE
    MPI_CXX_COMPILER
    MPIEXEC_EXECUTABLE)
  task22_forward_cache("${variable}" FILEPATH)
endforeach()

foreach(variable IN ITEMS
    CMAKE_C_COMPILER_TARGET
    CMAKE_CXX_COMPILER_TARGET
    CMAKE_GENERATOR_INSTANCE
    CMAKE_SYSROOT)
  task22_forward_cache("${variable}" STRING)
endforeach()

foreach(variable IN ITEMS
    HUNDUN_ENABLE_ASAN
    HUNDUN_ENABLE_UBSAN
    MPI_CXX_SKIP_MPICXX)
  task22_forward_cache("${variable}" BOOL)
endforeach()

if(HUNDUN_NESTED_CONFIG)
  list(APPEND configure_command
       "-DCMAKE_BUILD_TYPE:STRING=${HUNDUN_NESTED_CONFIG}")
  string(TOUPPER "${HUNDUN_NESTED_CONFIG}" config_upper)
  list(APPEND configure_command
       "-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_${config_upper}:PATH=${HUNDUN_ARCHIVE_DIR}"
       "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_${config_upper}:PATH=${HUNDUN_ARCHIVE_DIR}")
elseif(outer_CMAKE_BUILD_TYPE)
  list(APPEND configure_command
       "-DCMAKE_BUILD_TYPE:STRING=${outer_CMAKE_BUILD_TYPE}")
endif()

file(REMOVE_RECURSE "${HUNDUN_NESTED_BINARY_DIR}")
file(MAKE_DIRECTORY "${HUNDUN_ARCHIVE_DIR}")

execute_process(
  COMMAND ${configure_command}
  RESULT_VARIABLE configure_result
  OUTPUT_VARIABLE configure_output
  ERROR_VARIABLE configure_error)
if(NOT configure_result EQUAL 0)
  message(FATAL_ERROR
    "Task22 tests-off configure failed:\n${configure_output}${configure_error}")
endif()

set(build_command
    "${CMAKE_COMMAND}"
    --build "${HUNDUN_NESTED_BINARY_DIR}"
    --target hundun hundun_flow hundun_linear hundun_material_diagnostics
    --parallel 2)
if(HUNDUN_NESTED_CONFIG)
  list(APPEND build_command --config "${HUNDUN_NESTED_CONFIG}")
endif()
execute_process(
  COMMAND ${build_command}
  RESULT_VARIABLE build_result
  OUTPUT_VARIABLE build_output
  ERROR_VARIABLE build_error)
if(NOT build_result EQUAL 0)
  message(FATAL_ERROR
    "Task22 tests-off build failed:\n${build_output}${build_error}")
endif()
