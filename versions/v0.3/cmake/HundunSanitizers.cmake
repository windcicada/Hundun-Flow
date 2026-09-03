# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

function(hundun_configure_sanitizers target)
  set(enabled_sanitizers)
  if(HUNDUN_ENABLE_ASAN)
    list(APPEND enabled_sanitizers address)
  endif()
  if(HUNDUN_ENABLE_UBSAN)
    list(APPEND enabled_sanitizers undefined)
  endif()

  if(enabled_sanitizers)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      message(FATAL_ERROR "Sanitizers require a GNU or Clang C++ compiler")
    endif()

    list(JOIN enabled_sanitizers "," sanitizer_list)
    target_compile_options("${target}" INTERFACE
      "-fsanitize=${sanitizer_list}" -fno-omit-frame-pointer)
    target_link_options("${target}" INTERFACE
      "-fsanitize=${sanitizer_list}")
  endif()
endfunction()
