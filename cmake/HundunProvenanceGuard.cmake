# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

function(_hundun_collect_cmake_command_views candidate_text
         syntax_output execution_output)
  set(syntax_text "")
  set(execution_text "")
  set(execution_commands
    execute_process add_custom_command add_custom_target)
  set(pending_command "")
  set(current_command "")
  set(current_tokens "")
  set(current_token "")
  set(parenthesis_depth 0)
  set(in_comment FALSE)
  set(in_quote FALSE)
  set(quote_escape FALSE)
  set(in_bracket_comment FALSE)
  set(in_bracket_argument FALSE)
  set(bracket_closer "")

  string(LENGTH "${candidate_text}" candidate_length)
  set(candidate_index 0)
  while(candidate_index LESS candidate_length)
    string(SUBSTRING "${candidate_text}" ${candidate_index} 1 character)
    math(EXPR candidate_index "${candidate_index} + 1")

    if(in_bracket_comment OR in_bracket_argument)
      set(closes_bracket FALSE)
      if(character STREQUAL "]")
        string(LENGTH "${bracket_closer}" bracket_closer_length)
        math(EXPR bracket_candidate_index "${candidate_index} - 1")
        math(EXPR bracket_remaining_length
          "${candidate_length} - ${bracket_candidate_index}")
        if(bracket_remaining_length GREATER_EQUAL bracket_closer_length)
          string(SUBSTRING "${candidate_text}" ${bracket_candidate_index}
            ${bracket_closer_length} bracket_candidate)
          if(bracket_candidate STREQUAL bracket_closer)
            set(closes_bracket TRUE)
          endif()
        endif()
      endif()

      if(closes_bracket)
        math(EXPR candidate_index
          "${candidate_index} + ${bracket_closer_length} - 1")
        if(in_bracket_argument)
          list(APPEND current_tokens "${current_token}")
          set(current_token "")
        endif()
        set(in_bracket_comment FALSE)
        set(in_bracket_argument FALSE)
        set(bracket_closer "")
      elseif(in_bracket_argument)
        if(character MATCHES "[ \t\r\n;]")
          string(APPEND current_token "~")
        else()
          string(APPEND current_token "${character}")
        endif()
      endif()
      continue()
    endif()

    if(in_comment)
      if(character STREQUAL "\n")
        set(in_comment FALSE)
      endif()
      continue()
    endif()

    if(parenthesis_depth EQUAL 0)
      if(character STREQUAL "#")
        set(pending_command "")
        string(SUBSTRING "${candidate_text}" ${candidate_index} -1
          bracket_remainder)
        if(bracket_remainder MATCHES "^\\[(=*)\\[")
          set(bracket_equals "${CMAKE_MATCH_1}")
          string(LENGTH "${bracket_equals}" bracket_equals_length)
          set(bracket_closer "]${bracket_equals}]")
          set(in_bracket_comment TRUE)
          math(EXPR candidate_index
            "${candidate_index} + ${bracket_equals_length} + 2")
        else()
          set(in_comment TRUE)
        endif()
      elseif(character MATCHES "[a-z0-9_]")
        string(APPEND pending_command "${character}")
      elseif(character MATCHES "[ \t\r\n]")
        # Whitespace between a command name and its opening parenthesis is valid.
      elseif(character MATCHES "^[(]$")
        if(NOT pending_command STREQUAL "")
          set(current_command "${pending_command}")
          set(pending_command "")
          set(current_tokens "")
          set(current_token "")
          set(parenthesis_depth 1)
        endif()
      else()
        set(pending_command "")
      endif()
      continue()
    endif()

    if(in_quote)
      if(quote_escape)
        if(character MATCHES "[ \t\r\n;]")
          string(APPEND current_token "~")
        else()
          string(APPEND current_token "${character}")
        endif()
        set(quote_escape FALSE)
      elseif(character STREQUAL "\\")
        set(quote_escape TRUE)
      elseif(character STREQUAL "\"")
        set(in_quote FALSE)
      elseif(character MATCHES "[ \t\r\n;]")
        string(APPEND current_token "~")
      else()
        string(APPEND current_token "${character}")
      endif()
      continue()
    endif()

    if(character STREQUAL "#")
      if(NOT current_token STREQUAL "")
        list(APPEND current_tokens "${current_token}")
        set(current_token "")
      endif()
      string(SUBSTRING "${candidate_text}" ${candidate_index} -1
        bracket_remainder)
      if(bracket_remainder MATCHES "^\\[(=*)\\[")
        set(bracket_equals "${CMAKE_MATCH_1}")
        string(LENGTH "${bracket_equals}" bracket_equals_length)
        set(bracket_closer "]${bracket_equals}]")
        set(in_bracket_comment TRUE)
        math(EXPR candidate_index
          "${candidate_index} + ${bracket_equals_length} + 2")
      else()
        set(in_comment TRUE)
      endif()
    elseif(character STREQUAL "\"")
      set(in_quote TRUE)
    elseif(character STREQUAL "[")
      math(EXPR bracket_candidate_index "${candidate_index} - 1")
      string(SUBSTRING "${candidate_text}" ${bracket_candidate_index} -1
        bracket_remainder)
      if(bracket_remainder MATCHES "^\\[(=*)\\[")
        if(NOT current_token STREQUAL "")
          list(APPEND current_tokens "${current_token}")
          set(current_token "")
        endif()
        set(bracket_equals "${CMAKE_MATCH_1}")
        string(LENGTH "${bracket_equals}" bracket_equals_length)
        set(bracket_closer "]${bracket_equals}]")
        set(in_bracket_argument TRUE)
        math(EXPR candidate_index
          "${candidate_index} + ${bracket_equals_length} + 1")
      else()
        string(APPEND current_token "${character}")
      endif()
    elseif(character MATCHES "[ \t\r\n]")
      if(NOT current_token STREQUAL "")
        list(APPEND current_tokens "${current_token}")
        set(current_token "")
      endif()
    elseif(character MATCHES "^[(]$")
      if(NOT current_token STREQUAL "")
        list(APPEND current_tokens "${current_token}")
        set(current_token "")
      endif()
      math(EXPR parenthesis_depth "${parenthesis_depth} + 1")
    elseif(character MATCHES "^[)]$")
      if(NOT current_token STREQUAL "")
        list(APPEND current_tokens "${current_token}")
        set(current_token "")
      endif()
      math(EXPR parenthesis_depth "${parenthesis_depth} - 1")
      if(parenthesis_depth EQUAL 0)
        string(JOIN " " command_arguments ${current_tokens})
        set(command_record "${current_command}")
        if(NOT command_arguments STREQUAL "")
          string(APPEND command_record " ${command_arguments}")
        endif()
        string(APPEND syntax_text ";${command_record}")
        list(FIND execution_commands "${current_command}"
          execution_command_index)
        if(NOT execution_command_index EQUAL -1)
          string(APPEND execution_text ";${command_record}")
        endif()
        set(current_command "")
        set(current_tokens "")
      endif()
    elseif(character STREQUAL ";")
      string(APPEND current_token "~")
    else()
      string(APPEND current_token "${character}")
    endif()
  endwhile()

  set(${syntax_output} "${syntax_text}" PARENT_SCOPE)
  set(${execution_output} "${execution_text}" PARENT_SCOPE)
endfunction()

function(_hundun_collect_yaml_run_view candidate_text output_variable)
  string(REPLACE "\r\n" "\n" yaml_text "${candidate_text}")
  string(REPLACE "\r" "\n" yaml_text "${yaml_text}")
  set(yaml_line_marker "__hundun_yaml_line__")
  string(REPLACE ";" "\\;" yaml_lines "${yaml_text}")
  string(REPLACE "\n" ";${yaml_line_marker}" yaml_lines "${yaml_lines}")
  set(yaml_lines "${yaml_line_marker}${yaml_lines}")
  set(run_text "")
  set(in_run_block FALSE)
  set(run_block_indent 0)
  set(run_block_style "")
  set(folded_record_open FALSE)

  string(LENGTH "${yaml_line_marker}" yaml_line_marker_length)
  foreach(yaml_line_record IN LISTS yaml_lines)
    string(SUBSTRING "${yaml_line_record}" ${yaml_line_marker_length} -1
      yaml_line)
    string(LENGTH "${yaml_line}" yaml_line_length)
    string(REGEX REPLACE "^[ \t]+" "" yaml_line_without_indent
      "${yaml_line}")
    string(LENGTH "${yaml_line_without_indent}" yaml_content_length)
    math(EXPR yaml_indent_length
      "${yaml_line_length} - ${yaml_content_length}")
    string(STRIP "${yaml_line}" yaml_line_stripped)
    set(check_run_key TRUE)

    if(in_run_block)
      if(yaml_line_stripped STREQUAL "")
        if(run_block_style STREQUAL ">")
          set(folded_record_open FALSE)
        endif()
        set(check_run_key FALSE)
      elseif(yaml_indent_length GREATER run_block_indent)
        string(REGEX REPLACE "^[ \t]+" "" run_line "${yaml_line}")
        if(run_block_style STREQUAL ">")
          if(folded_record_open)
            string(APPEND run_text " ${run_line}")
          else()
            string(APPEND run_text ";${run_line}")
            set(folded_record_open TRUE)
          endif()
        else()
          string(APPEND run_text ";${run_line}")
        endif()
        set(check_run_key FALSE)
      else()
        set(in_run_block FALSE)
        set(run_block_style "")
        set(folded_record_open FALSE)
      endif()
    endif()

    if(check_run_key AND yaml_line MATCHES
       "^([ \t]*)(-[ \t]+)?run:[ \t]*(.*)$")
      set(run_key_prefix "${CMAKE_MATCH_1}${CMAKE_MATCH_2}")
      string(LENGTH "${run_key_prefix}" run_block_indent)
      string(STRIP "${CMAKE_MATCH_3}" run_value)
      if(run_value MATCHES "^([|>])[-+]?$")
        set(in_run_block TRUE)
        set(run_block_style "${CMAKE_MATCH_1}")
        set(folded_record_open FALSE)
      elseif(NOT run_value STREQUAL "")
        string(APPEND run_text ";${run_value}")
      endif()
    endif()
  endforeach()

  set(${output_variable} "${run_text}" PARENT_SCOPE)
endfunction()

function(hundun_assert_clean_tree root)
  file(GLOB_RECURSE forbidden_fortran
    "${root}/*.f" "${root}/*.F" "${root}/*.for"
    "${root}/*.f90" "${root}/*.F90" "${root}/*.inc")
  if(forbidden_fortran)
    list(GET forbidden_fortran 0 offending_path)
    message(FATAL_ERROR "Forbidden legacy-language files: ${offending_path}")
  endif()
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
    if(NOT relative_path MATCHES "^\\.github(/|$)"
       AND relative_path MATCHES "(^|/)\\.[^/]+(/|$)")
      continue()
    endif()
    if(relative_path MATCHES "(^|/)build([-/]|$)")
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

    set(is_cmake_control FALSE)
    set(is_shell_control FALSE)
    set(is_make_control FALSE)
    set(is_yaml_control FALSE)
    if(candidate_name STREQUAL "CMakeLists.txt"
       OR candidate_suffix STREQUAL "cmake")
      set(is_cmake_control TRUE)
    elseif(candidate_suffix STREQUAL "sh"
           OR candidate_name STREQUAL "configure")
      set(is_shell_control TRUE)
    elseif(candidate_name STREQUAL "Makefile"
           OR candidate_name STREQUAL "makefile"
           OR candidate_name STREQUAL "GNUmakefile")
      set(is_shell_control TRUE)
      set(is_make_control TRUE)
    elseif(candidate_suffix STREQUAL "yml"
           OR candidate_suffix STREQUAL "yaml")
      set(is_yaml_control TRUE)
    endif()

    file(READ "${candidate}" candidate_text)
    string(TOLOWER "${candidate_text}" candidate_text_lower)
    set(cmake_syntax_text "")
    set(cmake_execution_text "")
    set(shell_command_text "")
    if(is_cmake_control)
      _hundun_collect_cmake_command_views("${candidate_text_lower}"
        cmake_syntax_text cmake_execution_text)
    elseif(is_shell_control)
      set(shell_command_text "${candidate_text_lower}")
      if(is_make_control)
        string(REGEX REPLACE "(^|\n)\t[ \t]*[@+-]+" "\\1\t"
          shell_command_text "${shell_command_text}")
      endif()
      string(REPLACE "\n" ";" shell_command_text "${shell_command_text}")
    elseif(is_yaml_control)
      _hundun_collect_yaml_run_view("${candidate_text_lower}"
        shell_command_text)
    endif()
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
    if(is_cmake_control AND cmake_syntax_text MATCHES
       "(^|;) *find_package +(${scripting_runtime}([0-9]+)?|${scripting_runtime}interp|${scripting_runtime}libs(new)?|pybind11|nanobind)( |;|$)")
      message(FATAL_ERROR
        "Forbidden Python package discovery in ${candidate}")
    endif()
    set(has_interpreter_dependency FALSE)
    if(is_cmake_control AND
       (cmake_syntax_text MATCHES
          "${scripting_runtime}([0-9]+)?::interpreter"
        OR cmake_syntax_text MATCHES
          "${scripting_runtime}([0-9]+)?_executable"
        OR cmake_syntax_text MATCHES
          "(^|;) *find_program +[^;]*${scripting_runtime}"
        OR cmake_syntax_text MATCHES
          "(^|;) *${scripting_runtime}([0-9]+)?_add_(library|executable|module)( |;|$)"
        OR cmake_execution_text MATCHES
          "(^| )command +${scripting_runtime}([0-9]+([.][0-9]+)*)?( |;|$)"))
      set(has_interpreter_dependency TRUE)
    endif()
    if(shell_command_text MATCHES
       "${shell_command_prefix}${shell_executable_prefix}${scripting_runtime}([0-9]+([.][0-9]+)*)?( |[:;|&]|$)")
      set(has_interpreter_dependency TRUE)
    endif()
    if(has_interpreter_dependency)
      message(FATAL_ERROR
        "Forbidden Python interpreter dependency in ${candidate}")
    endif()

    foreach(package_tool IN ITEMS
        "${package_tool_a}" "${package_tool_a3}" "${package_tool_b}")
      set(has_package_command FALSE)
      if(is_cmake_control AND
         (cmake_execution_text MATCHES
           "(^| )command +${package_tool}( |$)"
          OR cmake_execution_text MATCHES
           "(^| )command +[^ ]+[/]${package_tool}( |$)"
          OR cmake_execution_text MATCHES
           "(^| )command +/usr/bin/env +${package_tool}( |$)"))
        set(has_package_command TRUE)
      endif()
      if(shell_command_text MATCHES
         "${shell_command_prefix}${shell_executable_prefix}${package_tool}( |[:;|&]|$)")
        set(has_package_command TRUE)
      endif()
      if(has_package_command)
        message(FATAL_ERROR
          "Forbidden direct package-manager command '${package_tool}' in ${candidate}")
      endif()
    endforeach()

    if(is_cmake_control AND
       (cmake_syntax_text MATCHES
          "(^|;) *include +${retrieval_a_token}( |;|$)"
        OR cmake_syntax_text MATCHES
          "(^|;) *${retrieval_a_token}_(declare|makeavailable|populate)( |;|$)"))
      message(FATAL_ERROR
        "Forbidden ${retrieval_a_display} source retrieval in ${candidate}")
    endif()
    if(is_cmake_control AND
       (cmake_syntax_text MATCHES
          "(^|;) *include +${retrieval_b_token}( |;|$)"
        OR cmake_syntax_text MATCHES
          "(^|;) *${retrieval_b_token}_(add|add_step)( |;|$)"))
      message(FATAL_ERROR
        "Forbidden ${retrieval_b_display} source retrieval in ${candidate}")
    endif()
    if(is_cmake_control AND cmake_syntax_text MATCHES
       "(^|;) *file +download( |;|$)")
      message(FATAL_ERROR
        "Forbidden file download source retrieval in ${candidate}")
    endif()
    if((is_cmake_control AND
        (cmake_execution_text MATCHES
           "(^| )command +${cmake_executable_prefix}git +clone( |;|$)"
         OR cmake_execution_text MATCHES
           "(^| )command +${cmake_executable_prefix}(curl|wget) +[^;]*https?://"))
       OR shell_command_text MATCHES
          "${shell_command_prefix}${shell_executable_prefix}git +clone( |$)"
       OR shell_command_text MATCHES
          "${shell_command_prefix}${shell_executable_prefix}(curl|wget) +[^;]*https?://")
      message(FATAL_ERROR
        "Forbidden command-line source retrieval in ${candidate}")
    endif()
  endforeach()
endfunction()
