# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

function(hundun_source_status_metadata status_result status
         output_clean output_summary)
  if(NOT "${status_result}" STREQUAL "0")
    set(clean "false")
    set(summary "unavailable")
  elseif("${status}" STREQUAL "")
    set(clean "true")
    set(summary "")
  else()
    string(SHA256 digest "${status}")
    set(clean "false")
    set(summary "git-status-sha256:${digest}")
  endif()
  set("${output_clean}" "${clean}" PARENT_SCOPE)
  set("${output_summary}" "${summary}" PARENT_SCOPE)
endfunction()

function(hundun_escape_cpp_string input output)
  set(escaped "${input}")
  string(REPLACE "\\" "\\\\" escaped "${escaped}")
  string(REPLACE "\"" "\\\"" escaped "${escaped}")
  string(REPLACE "\n" "\\n" escaped "${escaped}")
  string(REPLACE "\r" "\\r" escaped "${escaped}")
  string(REPLACE "\t" "\\t" escaped "${escaped}")
  set("${output}" "${escaped}" PARENT_SCOPE)
endfunction()
