# SPDX-License-Identifier: Apache-2.0

function(hundun_assert_clean_tree root)
  file(GLOB_RECURSE forbidden_fortran
    "${root}/*.f" "${root}/*.F" "${root}/*.for"
    "${root}/*.f90" "${root}/*.F90" "${root}/*.inc")
  if(forbidden_fortran)
    list(GET forbidden_fortran 0 offending_path)
    message(FATAL_ERROR "Forbidden legacy-language files: ${offending_path}")
  endif()

  file(GLOB_RECURSE candidate_files LIST_DIRECTORIES false "${root}/*")
  set(candidate_suffixes c h cc hh cpp hpp cxx hxx cmake)

  set(forbidden_tokens
    "boffin"
    "coast_legacy"
    "domxch"
    "coalesced_legacy_block")

  foreach(candidate IN LISTS candidate_files)
    file(RELATIVE_PATH relative_path "${root}" "${candidate}")
    if(relative_path STREQUAL "cmake/HundunProvenanceGuard.cmake"
       OR relative_path MATCHES "(^|/)third_party(/|$)")
      continue()
    endif()

    get_filename_component(candidate_name "${candidate}" NAME)
    get_filename_component(candidate_suffix "${candidate}" LAST_EXT)
    string(REGEX REPLACE "^\\." "" candidate_suffix "${candidate_suffix}")
    string(TOLOWER "${candidate_suffix}" candidate_suffix)
    list(FIND candidate_suffixes "${candidate_suffix}" candidate_suffix_index)
    if(NOT candidate_name STREQUAL "CMakeLists.txt"
       AND candidate_suffix_index EQUAL -1)
      continue()
    endif()

    file(READ "${candidate}" candidate_text)
    string(TOLOWER "${candidate_text}" candidate_text_lower)
    foreach(forbidden_token IN LISTS forbidden_tokens)
      string(FIND "${candidate_text_lower}" "${forbidden_token}" token_position)
      if(NOT token_position EQUAL -1)
        message(FATAL_ERROR
          "Forbidden provenance token '${forbidden_token}' in ${candidate}")
      endif()
    endforeach()
  endforeach()
endfunction()

function(hundun_assert_public_dependency_policy root)
  if(NOT EXISTS "${root}")
    message(FATAL_ERROR
      "Public dependency scan root does not exist: ${root}")
  endif()
  if(NOT IS_DIRECTORY "${root}")
    message(FATAL_ERROR
      "Public dependency scan root is not a directory: ${root}")
  endif()

  file(GLOB_RECURSE candidate_files LIST_DIRECTORIES false
    "${root}/*" "${root}/.github/*")
  set(control_suffixes
    c h cc hh cpp hpp cxx hxx cmake sh yml yaml json)
  set(vendored_control_suffixes cmake sh yml yaml)
  set(control_basenames
    CMakeLists.txt Makefile makefile GNUmakefile configure)
  string(CONCAT scripting_runtime "py" "thon")
  string(CONCAT package_tool_a "p" "ip")
  string(CONCAT package_tool_a3 "p" "ip3")
  string(CONCAT package_tool_b "con" "da")
  string(CONCAT retrieval_a_display "Fetch" "Content")
  string(CONCAT retrieval_b_display "External" "Project")
  string(TOLOWER "${retrieval_a_display}" retrieval_a_token)
  string(TOLOWER "${retrieval_b_display}" retrieval_b_token)

  foreach(candidate IN LISTS candidate_files)
    file(RELATIVE_PATH relative_path "${root}" "${candidate}")
    if(relative_path MATCHES "(^|/)\\.git(/|$)"
       OR relative_path MATCHES "(^|/)\\.codegraphf(/|$)"
       OR relative_path MATCHES "(^|/)\\.agents(/|$)"
       OR relative_path MATCHES "(^|/)\\.codex(/|$)"
       OR relative_path MATCHES "(^|/)\\.superpowers(/|$)"
       OR relative_path MATCHES "(^|/)\\.cache(/|$)"
       OR relative_path MATCHES "(^|/)build([-/]|$)")
      continue()
    endif()

    get_filename_component(candidate_name "${candidate}" NAME)
    get_filename_component(candidate_suffix "${candidate}" LAST_EXT)
    string(REGEX REPLACE "^\\." "" candidate_suffix "${candidate_suffix}")
    string(TOLOWER "${candidate_suffix}" candidate_suffix)
    list(FIND control_basenames "${candidate_name}" basename_index)
    if(candidate_suffix STREQUAL "py")
      message(FATAL_ERROR
        "Forbidden public Python source file: ${candidate}")
    endif()

    if(relative_path STREQUAL "cmake/HundunProvenanceGuard.cmake"
       OR relative_path MATCHES "(^|/)docs(/|$)")
      continue()
    endif()

    if(relative_path MATCHES "(^|/)third_party(/|$)")
      list(FIND vendored_control_suffixes
        "${candidate_suffix}" suffix_index)
      if(basename_index EQUAL -1 AND suffix_index EQUAL -1)
        continue()
      endif()
    else()
      list(FIND control_suffixes "${candidate_suffix}" suffix_index)
      if(basename_index EQUAL -1 AND suffix_index EQUAL -1)
        continue()
      endif()
    endif()

    file(READ "${candidate}" candidate_text)
    string(TOLOWER "${candidate_text}" candidate_text_lower)
    string(REPLACE "\n" ";" command_text "${candidate_text_lower}")
    set(cmake_command_text "${candidate_text_lower}")
    foreach(command_separator IN ITEMS "\n" "\r" "\t" "(" ")" "\"" "'")
      string(REPLACE "${command_separator}" " " cmake_command_text
        "${cmake_command_text}")
    endforeach()
    string(REGEX REPLACE " +" " " cmake_command_text
      "${cmake_command_text}")
    set(shell_command_text "${command_text}")
    foreach(command_quote IN ITEMS "\"" "'")
      string(REPLACE "${command_quote}" "" shell_command_text
        "${shell_command_text}")
    endforeach()
    string(REPLACE "\t" " " shell_command_text "${shell_command_text}")
    set(shell_command_prefix
      "(^|[:;|&]) *-* *((if|elif|while|until) +)?")
    set(shell_executable_prefix
      "(/usr/bin/env +|[^ ;:|&()]+[/])?")
    set(cmake_executable_prefix
      "(/usr/bin/env[ \\t\\r\\n]+|[-+./_a-z0-9]+[/])?")

    if(candidate_text_lower MATCHES
       "#[ \\t\\r\\n]*include[ \\t\\r\\n]*[<\"][ \\t\\r\\n]*(${scripting_runtime}[0-9.]*[/])?${scripting_runtime}[.]h[ \\t\\r\\n]*[>\"]")
      message(FATAL_ERROR
        "Forbidden Python header dependency in ${candidate}")
    endif()
    if(candidate_text_lower MATCHES
       "find_package[ \\t\\r\\n]*\\([ \\t\\r\\n]*(${scripting_runtime}([0-9]+)?|${scripting_runtime}interp|${scripting_runtime}libs(new)?|pybind11|nanobind)([ \\t\\r\\n]|\\))")
      message(FATAL_ERROR
        "Forbidden Python package discovery in ${candidate}")
    endif()
    if(candidate_text_lower MATCHES
          "${scripting_runtime}([0-9]+)?::interpreter"
       OR candidate_text_lower MATCHES
          "${scripting_runtime}([0-9]+)?_executable"
       OR candidate_text_lower MATCHES
          "find_program[ \\t\\r\\n]*\\([^)]*${scripting_runtime}"
       OR candidate_text_lower MATCHES
          "${scripting_runtime}([0-9]+)?_add_(library|executable|module)[ \\t\\r\\n]*\\("
       OR candidate_text_lower MATCHES
          "command[ \\t\\r\\n]+${scripting_runtime}([0-9]+([.][0-9]+)*)?([ \\t\\r\\n]|\\))"
       OR command_text MATCHES
          "(^|[:;|&])[ \\t-]*(/usr/bin/env[ \\t]+)?${scripting_runtime}([0-9]+([.][0-9]+)*)?([ \\t]|$)")
      message(FATAL_ERROR
        "Forbidden Python interpreter dependency in ${candidate}")
    endif()

    foreach(package_tool IN ITEMS
        "${package_tool_a}" "${package_tool_a3}" "${package_tool_b}")
      if(cmake_command_text MATCHES
           "(^| )command +${package_tool}( |$)"
         OR cmake_command_text MATCHES
           "(^| )command +[^ ]+[/]${package_tool}( |$)"
         OR cmake_command_text MATCHES
           "(^| )command +/usr/bin/env +${package_tool}( |$)"
         OR shell_command_text MATCHES
           "${shell_command_prefix}${shell_executable_prefix}${package_tool}( |[:;|&]|$)")
        message(FATAL_ERROR
          "Forbidden direct package-manager command '${package_tool}' in ${candidate}")
      endif()
    endforeach()

    if(candidate_text_lower MATCHES
       "include[ \\t\\r\\n]*\\([ \\t\\r\\n]*${retrieval_a_token}([ \\t\\r\\n]|\\))"
       OR candidate_text_lower MATCHES
          "${retrieval_a_token}_(declare|makeavailable|populate)[ \\t\\r\\n]*\\(")
      message(FATAL_ERROR
        "Forbidden ${retrieval_a_display} source retrieval in ${candidate}")
    endif()
    if(candidate_text_lower MATCHES
       "include[ \\t\\r\\n]*\\([ \\t\\r\\n]*${retrieval_b_token}([ \\t\\r\\n]|\\))"
       OR candidate_text_lower MATCHES
          "${retrieval_b_token}_(add|add_step)[ \\t\\r\\n]*\\(")
      message(FATAL_ERROR
        "Forbidden ${retrieval_b_display} source retrieval in ${candidate}")
    endif()
    if(candidate_text_lower MATCHES
       "file[ \\t\\r\\n]*\\([ \\t\\r\\n]*download([ \\t\\r\\n]|\\))")
      message(FATAL_ERROR
        "Forbidden file download source retrieval in ${candidate}")
    endif()
    if(candidate_text_lower MATCHES
       "command[ \\t\\r\\n]+${cmake_executable_prefix}git[ \\t\\r\\n]+clone([ \\t\\r\\n]|\\))"
       OR candidate_text_lower MATCHES
          "command[ \\t\\r\\n]+${cmake_executable_prefix}(curl|wget)[ \\t\\r\\n]+[^)]*https?://"
       OR shell_command_text MATCHES
          "${shell_command_prefix}${shell_executable_prefix}git +clone( |$)"
       OR shell_command_text MATCHES
          "${shell_command_prefix}${shell_executable_prefix}(curl|wget) +[^;]*https?://")
      message(FATAL_ERROR
        "Forbidden command-line source retrieval in ${candidate}")
    endif()
  endforeach()
endfunction()
