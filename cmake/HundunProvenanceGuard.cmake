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
    if(candidate_suffix STREQUAL "py")
      message(FATAL_ERROR
        "Forbidden public Python source file: ${candidate}")
    endif()

    if(relative_path STREQUAL "cmake/HundunProvenanceGuard.cmake"
       OR relative_path MATCHES "(^|/)third_party(/|$)"
       OR relative_path MATCHES "(^|/)docs(/|$)")
      continue()
    endif()

    list(FIND control_suffixes "${candidate_suffix}" suffix_index)
    if(NOT candidate_name STREQUAL "CMakeLists.txt"
       AND suffix_index EQUAL -1)
      continue()
    endif()

    file(READ "${candidate}" candidate_text)
    string(TOLOWER "${candidate_text}" candidate_text_lower)
    string(REPLACE "\n" ";" command_text "${candidate_text_lower}")

    if(candidate_text_lower MATCHES
       "#[ \\t\\r\\n]*include[ \\t\\r\\n]*[<\"][ \\t\\r\\n]*(python[0-9.]*[/])?python[.]h[ \\t\\r\\n]*[>\"]")
      message(FATAL_ERROR
        "Forbidden Python header dependency in ${candidate}")
    endif()
    if(candidate_text_lower MATCHES
       "find_package[ \\t\\r\\n]*\\([ \\t\\r\\n]*(python([0-9]+)?|pythoninterp|pythonlibs(new)?|pybind11|nanobind)([ \\t\\r\\n]|\\))")
      message(FATAL_ERROR
        "Forbidden Python package discovery in ${candidate}")
    endif()
    if(candidate_text_lower MATCHES "python([0-9]+)?::interpreter"
       OR candidate_text_lower MATCHES "python([0-9]+)?_executable"
       OR candidate_text_lower MATCHES
          "find_program[ \\t\\r\\n]*\\([^)]*python"
       OR candidate_text_lower MATCHES
          "python([0-9]+)?_add_(library|executable|module)[ \\t\\r\\n]*\\("
       OR candidate_text_lower MATCHES
          "command[ \\t\\r\\n]+python([0-9]+([.][0-9]+)*)?([ \\t\\r\\n]|\\))"
       OR command_text MATCHES
          "(^|[:;|&])[ \\t-]*(/usr/bin/env[ \\t]+)?python([0-9]+([.][0-9]+)*)?([ \\t]|$)")
      message(FATAL_ERROR
        "Forbidden Python interpreter dependency in ${candidate}")
    endif()

    if(candidate_text_lower MATCHES
       "include[ \\t\\r\\n]*\\([ \\t\\r\\n]*fetchcontent([ \\t\\r\\n]|\\))"
       OR candidate_text_lower MATCHES
          "fetchcontent_(declare|makeavailable|populate)[ \\t\\r\\n]*\\(")
      message(FATAL_ERROR
        "Forbidden FetchContent source retrieval in ${candidate}")
    endif()
    if(candidate_text_lower MATCHES
       "include[ \\t\\r\\n]*\\([ \\t\\r\\n]*externalproject([ \\t\\r\\n]|\\))"
       OR candidate_text_lower MATCHES
          "externalproject_(add|add_step)[ \\t\\r\\n]*\\(")
      message(FATAL_ERROR
        "Forbidden ExternalProject source retrieval in ${candidate}")
    endif()
    if(candidate_text_lower MATCHES
       "file[ \\t\\r\\n]*\\([ \\t\\r\\n]*download([ \\t\\r\\n]|\\))")
      message(FATAL_ERROR
        "Forbidden file download source retrieval in ${candidate}")
    endif()
    if(candidate_text_lower MATCHES
       "command[ \\t\\r\\n]+git[ \\t\\r\\n]+clone([ \\t\\r\\n]|\\))"
       OR command_text MATCHES
          "(^|[:;|&])[ \\t-]*(/usr/bin/env[ \\t]+)?git[ \\t]+clone([ \\t]|$)"
       OR command_text MATCHES
          "(^|[:;|&])[ \\t-]*(/usr/bin/env[ \\t]+)?(curl|wget)[ \\t]+[^;]*https?://")
      message(FATAL_ERROR
        "Forbidden command-line source retrieval in ${candidate}")
    endif()
  endforeach()
endfunction()
