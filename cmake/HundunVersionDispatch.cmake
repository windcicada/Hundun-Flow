# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

function(hundun_add_selected_source_line source_version)
  if(NOT source_version STREQUAL "v0.4")
    message(FATAL_ERROR
      "unsupported HUNDUN_SOURCE_VERSION: ${source_version}; "
      "this checkout contains only the current implementation (v0.4); "
      "use Git history for retired implementations")
  endif()

  set(source_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/versions/${source_version}")
  if(NOT IS_DIRECTORY "${source_directory}")
    message(FATAL_ERROR
      "selected HUNDUN-FLOW source directory is missing: ${source_directory}")
  endif()

  add_subdirectory("${source_directory}" "versions/${source_version}")
endfunction()
