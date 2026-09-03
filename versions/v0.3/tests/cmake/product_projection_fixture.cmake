# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

foreach(required_variable IN ITEMS
    HUNDUN_SOURCE_ROOT
    HUNDUN_TEST_BINARY_ROOT
    HUNDUN_GENERATOR
    HUNDUN_C_COMPILER
    HUNDUN_CXX_COMPILER
    HUNDUN_C_COMPILER_ID
    HUNDUN_CXX_COMPILER_ID)
  if(NOT DEFINED ${required_variable} OR
     "${${required_variable}}" STREQUAL "")
    message(FATAL_ERROR
      "product projection input '${required_variable}' is required")
  endif()
endforeach()

foreach(path_variable IN ITEMS HUNDUN_SOURCE_ROOT HUNDUN_TEST_BINARY_ROOT)
  set(path_value "${${path_variable}}")
  cmake_path(ABSOLUTE_PATH path_value NORMALIZE
             OUTPUT_VARIABLE "${path_variable}_ABSOLUTE")
endforeach()

if(NOT IS_DIRECTORY "${HUNDUN_SOURCE_ROOT_ABSOLUTE}")
  message(FATAL_ERROR "product projection source root is not a directory")
endif()
if(NOT IS_DIRECTORY "${HUNDUN_TEST_BINARY_ROOT_ABSOLUTE}")
  message(FATAL_ERROR "product projection binary root is not a directory")
endif()

file(REAL_PATH "${HUNDUN_SOURCE_ROOT_ABSOLUTE}" product_source_root)
file(REAL_PATH "${HUNDUN_TEST_BINARY_ROOT_ABSOLUTE}" test_binary_root)
cmake_path(GET test_binary_root ROOT_PATH filesystem_root)
if(test_binary_root STREQUAL filesystem_root OR
   test_binary_root STREQUAL product_source_root)
  message(FATAL_ERROR "product projection binary root has unsafe authority")
endif()

if(NOT HUNDUN_C_COMPILER_ID STREQUAL "Clang" OR
   NOT HUNDUN_CXX_COMPILER_ID STREQUAL "Clang")
  message(FATAL_ERROR "product projection requires the Clang toolchain")
endif()
if(NOT DEFINED HUNDUN_CXX_FLAGS OR
   NOT HUNDUN_CXX_FLAGS MATCHES "(^| )-stdlib=libc\\+\\+($| )")
  message(FATAL_ERROR "product projection requires the libc++ C++ flags")
endif()

set(fixture_name "product-projection-fixture")
set(fixture_root "${test_binary_root}/${fixture_name}")
set(owner_marker "${fixture_root}/.hundun-product-projection-owner")
set(expected_marker "HUNDUN product projection fixture\n")
cmake_path(GET fixture_root PARENT_PATH fixture_parent)
cmake_path(GET fixture_root FILENAME observed_fixture_name)
if(NOT fixture_parent STREQUAL test_binary_root OR
   NOT observed_fixture_name STREQUAL fixture_name)
  message(FATAL_ERROR "product projection fixture escaped its binary root")
endif()
if(IS_SYMLINK "${fixture_root}")
  message(FATAL_ERROR "product projection fixture root cannot be a symlink")
endif()
if(EXISTS "${fixture_root}")
  if(NOT EXISTS "${owner_marker}" OR IS_SYMLINK "${owner_marker}")
    message(FATAL_ERROR
      "refusing to replace an unowned product projection fixture root")
  endif()
  file(READ "${owner_marker}" observed_marker)
  if(NOT observed_marker STREQUAL expected_marker)
    message(FATAL_ERROR
      "refusing to replace a product projection with another owner")
  endif()
  file(REMOVE_RECURSE "${fixture_root}")
endif()

set(projection_source "${fixture_root}/source")
set(projection_build "${fixture_root}/build")
file(MAKE_DIRECTORY "${projection_source}")
file(WRITE "${owner_marker}" "${expected_marker}")

foreach(required_file IN ITEMS
    CMakeLists.txt
    LICENSE
    NOTICE
    THIRD_PARTY.md)
  if(NOT EXISTS "${product_source_root}/${required_file}")
    message(FATAL_ERROR
      "required product projection file is missing: ${required_file}")
  endif()
  file(COPY "${product_source_root}/${required_file}"
       DESTINATION "${projection_source}")
endforeach()

foreach(required_directory IN ITEMS cmake include src third_party LICENSES)
  if(NOT IS_DIRECTORY "${product_source_root}/${required_directory}")
    message(FATAL_ERROR
      "required product projection directory is missing: ${required_directory}")
  endif()
  file(COPY "${product_source_root}/${required_directory}"
       DESTINATION "${projection_source}")
endforeach()

set(public_manual "docs/development/naming-and-style.md")
if(NOT EXISTS "${product_source_root}/${public_manual}")
  message(FATAL_ERROR "public development manual is missing")
endif()
file(MAKE_DIRECTORY "${projection_source}/docs/development")
file(COPY "${product_source_root}/${public_manual}"
     DESTINATION "${projection_source}/docs/development")

foreach(excluded_path IN ITEMS
    tests
    AGENTS.md
    CONTRIBUTING.md
    .github
    .superpowers)
  if(EXISTS "${projection_source}/${excluded_path}" OR
     IS_SYMLINK "${projection_source}/${excluded_path}")
    message(FATAL_ERROR
      "governance path leaked into product projection: ${excluded_path}")
  endif()
endforeach()

set(configure_command
    "${CMAKE_COMMAND}"
    -S "${projection_source}"
    -B "${projection_build}"
    -G "${HUNDUN_GENERATOR}"
    "-DCMAKE_BUILD_TYPE:STRING=Release"
    "-DCMAKE_C_COMPILER:FILEPATH=${HUNDUN_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER:FILEPATH=${HUNDUN_CXX_COMPILER}"
    "-DCMAKE_C_FLAGS:STRING=${HUNDUN_C_FLAGS}"
    "-DCMAKE_CXX_FLAGS:STRING=${HUNDUN_CXX_FLAGS}"
    "-DCMAKE_EXE_LINKER_FLAGS:STRING=${HUNDUN_EXE_LINKER_FLAGS}")
if(DEFINED HUNDUN_MAKE_PROGRAM AND NOT HUNDUN_MAKE_PROGRAM STREQUAL "")
  list(APPEND configure_command
    "-DCMAKE_MAKE_PROGRAM:FILEPATH=${HUNDUN_MAKE_PROGRAM}")
endif()

execute_process(
  COMMAND ${configure_command}
  RESULT_VARIABLE configure_result
  OUTPUT_VARIABLE configure_output
  ERROR_VARIABLE configure_error)
if(NOT configure_result EQUAL 0)
  message(FATAL_ERROR
    "tests-free product projection default configure failed:\n"
    "${configure_output}${configure_error}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" --build "${projection_build}"
          --target hundun --parallel 2
  RESULT_VARIABLE build_result
  OUTPUT_VARIABLE build_output
  ERROR_VARIABLE build_error)
if(NOT build_result EQUAL 0)
  message(FATAL_ERROR
    "tests-free product projection focused build failed:\n"
    "${build_output}${build_error}")
endif()
