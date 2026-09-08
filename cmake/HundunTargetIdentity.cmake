# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

# Recurse only through source directories. A recursive CMakeLists.txt glob
# from the repository root also matches nested build/compiler-probe files.
function(hundun_v04_identity_inputs root output)
  file(GLOB_RECURSE inputs CONFIGURE_DEPENDS
    "${root}/versions/v0.4/include/*"
    "${root}/versions/v0.4/src/*"
    "${root}/cmake/*"
    "${root}/third_party/yyjson/*")
  list(APPEND inputs "${root}/CMakeLists.txt" "${root}/versions/v0.4/CMakeLists.txt")
  set(${output} "${inputs}" PARENT_SCOPE)
endfunction()

# Sorted relative names and bytes, not timestamps or checkout locations.
function(hundun_source_content_digest root output)
  set(inputs ${ARGN})
  list(SORT inputs)
  set(material "HUNDUN_SOURCE_CONTENT_V1\n")
  foreach(path IN LISTS inputs)
    if(NOT IS_DIRECTORY "${path}")
      file(RELATIVE_PATH name "${root}" "${path}")
      file(SHA256 "${path}" digest)
      string(APPEND material "${name}=${digest}\n")
    endif()
  endforeach()
  string(SHA256 digest "${material}")
  set(${output} "${digest}" PARENT_SCOPE)
endfunction()

function(hundun_target_manifest target entry output)
  set(material "${HUNDUN_BUILD_MANIFEST_CANONICAL}")
  string(APPEND material "target=${target}\nentry_sha256=${entry}\n")
  # Core clean status alone says nothing about a runner-only dirty edit.
  set(target_clean "${HUNDUN_SOURCE_CLEAN}")
  if(ARGC GREATER 3 AND Git_FOUND)
    execute_process(COMMAND "${GIT_EXECUTABLE}" status --porcelain=v1
      --untracked-files=normal -- "${ARGV3}"
      WORKING_DIRECTORY "${hundun_repository_root}"
      RESULT_VARIABLE status_result OUTPUT_VARIABLE entry_status ERROR_QUIET)
    if(NOT status_result EQUAL 0 OR NOT entry_status STREQUAL "")
      set(target_clean false)
    endif()
  endif()
  string(APPEND material "target_source_clean=${target_clean}\n")
  string(SHA256 digest "${material}")
  file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/${target}.build-manifest.txt"
       CONTENT "${material}")
  set(${output} "${digest}" PARENT_SCOPE)
endfunction()
