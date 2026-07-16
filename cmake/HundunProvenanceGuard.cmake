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
