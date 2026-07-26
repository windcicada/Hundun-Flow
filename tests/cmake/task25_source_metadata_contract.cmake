# SPDX-License-Identifier: Apache-2.0

if(NOT DEFINED HUNDUN_SOURCE_ROOT)
  message(FATAL_ERROR "Task 25 source metadata contract input is missing")
endif()

include("${HUNDUN_SOURCE_ROOT}/applications/hundun/stage2_performance_source_metadata.cmake")

hundun_source_status_metadata(0 "" clean summary)
if(NOT clean OR NOT summary STREQUAL "")
  message(FATAL_ERROR "clean source metadata classification differs")
endif()

set(dirty_status
  " M path with quote-\" and backslash-\\ and semicolon-;\n?? new file")
hundun_source_status_metadata(0 "${dirty_status}" dirty_a summary_a)
hundun_source_status_metadata(0 "${dirty_status}" dirty_b summary_b)
string(LENGTH "${summary_a}" summary_a_length)
if(dirty_a OR dirty_b OR
   NOT summary_a STREQUAL summary_b OR
   NOT summary_a MATCHES "^git-status-sha256:[0-9a-f]+$" OR
   NOT summary_a_length EQUAL 82)
  message(FATAL_ERROR "stable dirty source metadata classification differs")
endif()

string(APPEND dirty_status "\n M another-file")
hundun_source_status_metadata(0 "${dirty_status}" dirty_c summary_c)
string(LENGTH "${summary_c}" summary_c_length)
if(dirty_c OR summary_c STREQUAL summary_a OR
   NOT summary_c MATCHES "^git-status-sha256:[0-9a-f]+$" OR
   NOT summary_c_length EQUAL 82)
  message(FATAL_ERROR "different dirty source metadata classification differs")
endif()

hundun_source_status_metadata(1 "ignored status" unavailable_clean
                              unavailable_summary)
if(unavailable_clean OR NOT unavailable_summary STREQUAL "unavailable")
  message(FATAL_ERROR "unavailable source metadata classification differs")
endif()

set(unescaped "quote-\" backslash-\\ newline-\n carriage-\r tab-\t")
hundun_escape_cpp_string("${unescaped}" escaped)
if(NOT escaped STREQUAL
   "quote-\\\" backslash-\\\\ newline-\\n carriage-\\r tab-\\t")
  message(FATAL_ERROR "C++ metadata escaping differs: ${escaped}")
endif()
