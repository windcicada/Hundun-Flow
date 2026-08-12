# SPDX-License-Identifier: Apache-2.0

function(hundun_configure_warnings target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options("${target}" INTERFACE
      "$<$<COMPILE_LANGUAGE:CXX>:-Wall>"
      "$<$<COMPILE_LANGUAGE:CXX>:-Wextra>"
      "$<$<COMPILE_LANGUAGE:CXX>:-Wpedantic>"
      "$<$<COMPILE_LANGUAGE:CXX>:-Wconversion>"
      "$<$<COMPILE_LANGUAGE:CXX>:-Wshadow>")
  endif()
endfunction()
