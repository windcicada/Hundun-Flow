# SPDX-License-Identifier: Apache-2.0

function(hundun_add_selected_source_line source_version)
  set(supported_versions v0.4 v0.3)
  if(NOT source_version IN_LIST supported_versions)
    message(FATAL_ERROR
      "unsupported HUNDUN_SOURCE_VERSION: ${source_version}; "
      "expected one of: v0.4, v0.3")
  endif()

  set(source_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/versions/${source_version}")
  if(NOT IS_DIRECTORY "${source_directory}")
    message(FATAL_ERROR
      "selected HUNDUN-FLOW source directory is missing: ${source_directory}")
  endif()

  add_subdirectory("${source_directory}" "versions/${source_version}")
endfunction()
