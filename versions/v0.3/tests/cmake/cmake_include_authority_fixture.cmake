# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

cmake_policy(SET CMP0057 NEW)

function(_hundun_authority_product_lists libraries_output targets_output)
  set(product_libraries
    hundun_execution
    hundun_diagnostics
    hundun_diagnostics_core
    hundun_config
    hundun_cli
    hundun_runtime
    hundun_application
    hundun_mesh
    hundun_boundary
    hundun_immersed
    hundun_linear
    hundun_fvm
    hundun_flow
    hundun_material_diagnostics
    hundun_checkpoint_diagnostics
    hundun_session_diagnostics
    hundun_io
    hundun_solver
    hundun_sdk)
  set(product_targets hundun_options ${product_libraries} hundun)
  set("${libraries_output}" "${product_libraries}" PARENT_SCOPE)
  set("${targets_output}" "${product_targets}" PARENT_SCOPE)
endfunction()

function(_hundun_authority_unwrap_build_interface value output)
  set(normalized "${value}")
  if(normalized MATCHES "^\\$<BUILD_INTERFACE:(.*)>$")
    set(normalized "${CMAKE_MATCH_1}")
  elseif(normalized MATCHES "^\\$<INSTALL_INTERFACE:.*>$")
    set(normalized "")
  endif()
  if(NOT normalized STREQUAL "" AND IS_ABSOLUTE "${normalized}")
    cmake_path(NORMAL_PATH normalized OUTPUT_VARIABLE normalized)
  endif()
  set("${output}" "${normalized}" PARENT_SCOPE)
endfunction()

function(_hundun_authority_assert_graph_path owner value)
  _hundun_authority_unwrap_build_interface("${value}" normalized)
  if(normalized STREQUAL "")
    return()
  endif()

  cmake_path(NORMAL_PATH PROJECT_SOURCE_DIR OUTPUT_VARIABLE source_root)
  cmake_path(NORMAL_PATH PROJECT_BINARY_DIR OUTPUT_VARIABLE binary_root)

  string(FIND "${normalized}" "${binary_root}/" binary_prefix)
  if(binary_prefix EQUAL 0)
    if(NOT normalized MATCHES "(^|/)generated(/|$)")
      message(FATAL_ERROR
        "target '${owner}' has unauthorized generated include authority: ${normalized}")
    endif()
    return()
  endif()

  if(normalized STREQUAL source_root)
    message(FATAL_ERROR
      "target '${owner}' has repository-root include authority")
  endif()

  string(FIND "${normalized}" "${source_root}/tests" tests_prefix)
  if(tests_prefix EQUAL 0)
    string(LENGTH "${source_root}/tests" tests_root_length)
    string(LENGTH "${normalized}" normalized_length)
    if(normalized_length EQUAL tests_root_length)
      message(FATAL_ERROR "target '${owner}' has tests include authority")
    endif()
    string(SUBSTRING "${normalized}" ${tests_root_length} 1 tests_suffix)
    if(tests_suffix STREQUAL "/")
      message(FATAL_ERROR "target '${owner}' has tests include authority")
    endif()
  endif()

  string(FIND "${normalized}" "${source_root}/" source_prefix)
  if(source_prefix EQUAL 0 AND
     NOT normalized STREQUAL "${source_root}/include" AND
     NOT normalized STREQUAL "${source_root}/src" AND
     NOT normalized MATCHES "^${source_root}/third_party(/|$)")
    message(FATAL_ERROR
      "target '${owner}' has unauthorized project include authority: ${normalized}")
  endif()

endfunction()

function(_hundun_authority_get_property target property output)
  get_target_property(value "${target}" "${property}")
  if(value MATCHES "-NOTFOUND$")
    set(value "")
  endif()
  set("${output}" "${value}" PARENT_SCOPE)
endfunction()

function(_hundun_authority_validate_product_target target is_library)
  _hundun_authority_get_property("${target}" INCLUDE_DIRECTORIES build_entries)
  _hundun_authority_get_property(
    "${target}" INTERFACE_INCLUDE_DIRECTORIES interface_entries)
  _hundun_authority_get_property(
    "${target}" INTERFACE_SYSTEM_INCLUDE_DIRECTORIES system_entries)

  set(expected_public "${PROJECT_SOURCE_DIR}/include")
  set(expected_private "${PROJECT_SOURCE_DIR}/src")
  set(expected_generated "${PROJECT_BINARY_DIR}/src/generated")
  cmake_path(NORMAL_PATH expected_public OUTPUT_VARIABLE expected_public)
  cmake_path(NORMAL_PATH expected_private OUTPUT_VARIABLE expected_private)
  cmake_path(NORMAL_PATH expected_generated OUTPUT_VARIABLE expected_generated)

  set(public_count 0)
  set(private_count 0)
  set(generated_count 0)
  foreach(entry IN LISTS build_entries)
    _hundun_authority_assert_graph_path("${target}" "${entry}")
    _hundun_authority_unwrap_build_interface("${entry}" normalized)
    if(normalized STREQUAL expected_public)
      math(EXPR public_count "${public_count} + 1")
    elseif(normalized STREQUAL expected_private)
      math(EXPR private_count "${private_count} + 1")
    elseif(target STREQUAL "hundun" AND
           normalized STREQUAL expected_generated)
      math(EXPR generated_count "${generated_count} + 1")
    else()
      message(FATAL_ERROR
        "product target '${target}' has unauthorized direct include authority: ${normalized}")
    endif()
  endforeach()

  if(is_library)
    if(NOT public_count EQUAL 1 OR NOT private_count EQUAL 1 OR
       NOT generated_count EQUAL 0)
      message(FATAL_ERROR
        "product library '${target}' must have one public include and one private src authority")
    endif()
  elseif(target STREQUAL "hundun_options")
    if(NOT public_count EQUAL 0 OR NOT private_count EQUAL 0 OR
       NOT generated_count EQUAL 0)
      message(FATAL_ERROR
        "hundun_options may not have direct include authority")
    endif()
  elseif(target STREQUAL "hundun")
    if(NOT public_count EQUAL 0 OR NOT private_count EQUAL 1 OR
       NOT generated_count EQUAL 1)
      message(FATAL_ERROR
        "hundun executable must have only private src and generated-header authority")
    endif()
  endif()

  set(interface_count 0)
  foreach(entry IN LISTS interface_entries)
    _hundun_authority_assert_graph_path("${target}" "${entry}")
    _hundun_authority_unwrap_build_interface("${entry}" normalized)
    if(normalized STREQUAL expected_public)
      math(EXPR interface_count "${interface_count} + 1")
    elseif(NOT normalized STREQUAL "")
      message(FATAL_ERROR
        "product target '${target}' exposes a non-public include directory: ${normalized}")
    endif()
  endforeach()
  if((is_library OR target STREQUAL "hundun_options") AND
     NOT interface_count EQUAL 1)
    message(FATAL_ERROR
      "product target '${target}' must expose exactly one public include authority")
  elseif(target STREQUAL "hundun" AND NOT interface_count EQUAL 0)
    message(FATAL_ERROR "hundun executable may not expose include authority")
  endif()

  foreach(entry IN LISTS system_entries)
    _hundun_authority_assert_graph_path("${target}" "${entry}")
    message(FATAL_ERROR
      "product target '${target}' may not directly expose a system include directory")
  endforeach()
endfunction()

function(_hundun_authority_validate_reachable_target target)
  get_property(visited GLOBAL PROPERTY HUNDUN_AUTHORITY_VISITED_TARGETS)
  if(target IN_LIST visited)
    return()
  endif()
  set_property(GLOBAL APPEND PROPERTY
    HUNDUN_AUTHORITY_VISITED_TARGETS "${target}")

  foreach(property IN ITEMS
      INCLUDE_DIRECTORIES
      INTERFACE_INCLUDE_DIRECTORIES
      INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    _hundun_authority_get_property("${target}" "${property}" entries)
    foreach(entry IN LISTS entries)
      _hundun_authority_assert_graph_path("${target}" "${entry}")
    endforeach()
  endforeach()

  foreach(property IN ITEMS LINK_LIBRARIES INTERFACE_LINK_LIBRARIES)
    _hundun_authority_get_property("${target}" "${property}" link_entries)
    foreach(link_entry IN LISTS link_entries)
      string(REGEX MATCHALL
        "[A-Za-z_][A-Za-z0-9_]*(::[A-Za-z_][A-Za-z0-9_]*)?"
        link_candidates "${link_entry}")
      foreach(candidate IN LISTS link_candidates)
        if(TARGET "${candidate}")
          _hundun_authority_validate_reachable_target("${candidate}")
        endif()
      endforeach()
    endforeach()
  endforeach()
endfunction()

function(hundun_validate_configured_product_include_authority)
  _hundun_authority_product_lists(product_libraries product_targets)
  foreach(target IN LISTS product_targets)
    if(NOT TARGET "${target}")
      message(FATAL_ERROR "required product target '${target}' is missing")
    endif()
    if(target IN_LIST product_libraries)
      _hundun_authority_validate_product_target("${target}" TRUE)
    else()
      _hundun_authority_validate_product_target("${target}" FALSE)
    endif()
  endforeach()

  set_property(GLOBAL PROPERTY HUNDUN_AUTHORITY_VISITED_TARGETS "")
  foreach(target IN LISTS product_targets)
    _hundun_authority_validate_reachable_target("${target}")
  endforeach()
endfunction()

if(NOT CMAKE_SCRIPT_MODE_FILE)
  return()
endif()

if(NOT DEFINED HUNDUN_LAYOUT_ROOT)
  message(FATAL_ERROR "HUNDUN_LAYOUT_ROOT is required")
endif()

get_filename_component(HUNDUN_LAYOUT_ROOT "${HUNDUN_LAYOUT_ROOT}" ABSOLUTE)
set(HUNDUN_SOURCE_CMAKE "${HUNDUN_LAYOUT_ROOT}/src/CMakeLists.txt")
if(NOT EXISTS "${HUNDUN_SOURCE_CMAKE}")
  message(FATAL_ERROR "product source CMakeLists.txt is missing")
endif()

file(READ "${HUNDUN_SOURCE_CMAKE}" HUNDUN_SOURCE_CMAKE_TEXT)
string(REGEX MATCHALL
  "target_include_directories\\([^\\)]*\\)"
  HUNDUN_INCLUDE_BLOCKS
  "${HUNDUN_SOURCE_CMAKE_TEXT}")
_hundun_authority_product_lists(
  HUNDUN_PRODUCT_LIBRARY_TARGETS HUNDUN_PRODUCT_TARGETS)

set(HUNDUN_EXPECTED_PUBLIC
  "$<BUILD_INTERFACE:\${PROJECT_SOURCE_DIR}/include>")
set(HUNDUN_EXPECTED_PRIVATE "\${PROJECT_SOURCE_DIR}/src")
set(HUNDUN_EXPECTED_GENERATED "\${CMAKE_CURRENT_BINARY_DIR}/generated")

foreach(HUNDUN_BLOCK IN LISTS HUNDUN_INCLUDE_BLOCKS)
  string(REGEX REPLACE
    "^target_include_directories\\([ \t\r\n]*" ""
    HUNDUN_BLOCK_BODY "${HUNDUN_BLOCK}")
  string(REGEX REPLACE "\\)[ \t\r\n]*$" ""
    HUNDUN_BLOCK_BODY "${HUNDUN_BLOCK_BODY}")
  string(REGEX REPLACE "[ \t\r\n]+" ";"
    HUNDUN_BLOCK_TOKENS "${HUNDUN_BLOCK_BODY}")
  list(POP_FRONT HUNDUN_BLOCK_TOKENS HUNDUN_TARGET)
  if(NOT HUNDUN_TARGET IN_LIST HUNDUN_PRODUCT_TARGETS)
    continue()
  endif()

  set(HUNDUN_SCOPE "")
  set(HUNDUN_PUBLIC_COUNT 0)
  set(HUNDUN_PRIVATE_COUNT 0)
  set(HUNDUN_GENERATED_COUNT 0)
  foreach(HUNDUN_TOKEN IN LISTS HUNDUN_BLOCK_TOKENS)
    string(REGEX REPLACE "^\"|\"$" "" HUNDUN_TOKEN "${HUNDUN_TOKEN}")
    if(HUNDUN_TOKEN MATCHES "^(PUBLIC|PRIVATE|INTERFACE)$")
      set(HUNDUN_SCOPE "${HUNDUN_TOKEN}")
      continue()
    endif()
    if(HUNDUN_TOKEN STREQUAL "\${PROJECT_SOURCE_DIR}" OR
       HUNDUN_TOKEN STREQUAL "$<BUILD_INTERFACE:\${PROJECT_SOURCE_DIR}>")
      message(FATAL_ERROR
        "product target '${HUNDUN_TARGET}' has repository-root include authority")
    endif()
    if(HUNDUN_TOKEN MATCHES
        "(^|[:/])\\$\\{PROJECT_SOURCE_DIR\\}/tests([/>]|$)")
      message(FATAL_ERROR
        "product target '${HUNDUN_TARGET}' has tests include authority")
    endif()

    if(HUNDUN_TARGET STREQUAL "hundun_options")
      if(NOT HUNDUN_SCOPE STREQUAL "INTERFACE" OR
         NOT HUNDUN_TOKEN STREQUAL "${HUNDUN_EXPECTED_PUBLIC}")
        message(FATAL_ERROR
          "hundun_options may expose only the public include directory")
      endif()
      continue()
    endif()
    if(HUNDUN_SCOPE STREQUAL "PUBLIC")
      if(NOT HUNDUN_TOKEN STREQUAL "${HUNDUN_EXPECTED_PUBLIC}")
        message(FATAL_ERROR
          "product target '${HUNDUN_TARGET}' exposes a non-public include directory")
      endif()
      math(EXPR HUNDUN_PUBLIC_COUNT "${HUNDUN_PUBLIC_COUNT} + 1")
    elseif(HUNDUN_SCOPE STREQUAL "PRIVATE")
      if(HUNDUN_TOKEN STREQUAL "${HUNDUN_EXPECTED_PRIVATE}" OR
         HUNDUN_TOKEN STREQUAL
           "$<BUILD_INTERFACE:\${PROJECT_SOURCE_DIR}/src>")
        math(EXPR HUNDUN_PRIVATE_COUNT "${HUNDUN_PRIVATE_COUNT} + 1")
      elseif(HUNDUN_TARGET STREQUAL "hundun" AND
             HUNDUN_TOKEN STREQUAL "${HUNDUN_EXPECTED_GENERATED}")
        math(EXPR HUNDUN_GENERATED_COUNT "${HUNDUN_GENERATED_COUNT} + 1")
      else()
        message(FATAL_ERROR
          "product target '${HUNDUN_TARGET}' has an unauthorized private include directory")
      endif()
    else()
      message(FATAL_ERROR
        "product target '${HUNDUN_TARGET}' has an include without explicit authority")
    endif()
  endforeach()

  if(HUNDUN_TARGET IN_LIST HUNDUN_PRODUCT_LIBRARY_TARGETS)
    if(NOT HUNDUN_PUBLIC_COUNT EQUAL 1 OR
       NOT HUNDUN_PRIVATE_COUNT EQUAL 1)
      message(FATAL_ERROR
        "product library '${HUNDUN_TARGET}' must have one public include and one private src authority")
    endif()
  elseif(HUNDUN_TARGET STREQUAL "hundun")
    if(NOT HUNDUN_PUBLIC_COUNT EQUAL 0 OR
       NOT HUNDUN_PRIVATE_COUNT EQUAL 1 OR
       NOT HUNDUN_GENERATED_COUNT EQUAL 1)
      message(FATAL_ERROR
        "hundun executable must have only private src and generated-header authority")
    endif()
  endif()
endforeach()

foreach(HUNDUN_TARGET IN LISTS HUNDUN_PRODUCT_TARGETS)
  set(HUNDUN_TARGET_SEEN FALSE)
  foreach(HUNDUN_BLOCK IN LISTS HUNDUN_INCLUDE_BLOCKS)
    if(HUNDUN_BLOCK MATCHES
        "^target_include_directories\\([ \t\r\n]*${HUNDUN_TARGET}([ \t\r\n]|\\))")
      set(HUNDUN_TARGET_SEEN TRUE)
      break()
    endif()
  endforeach()
  if(NOT HUNDUN_TARGET_SEEN)
    message(FATAL_ERROR
      "product target '${HUNDUN_TARGET}' has no explicit include authority")
  endif()
endforeach()
