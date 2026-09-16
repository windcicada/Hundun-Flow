# SPDX-License-Identifier: Apache-2.0
cmake_minimum_required(VERSION 3.21)
include("${CMAKE_CURRENT_LIST_DIR}/../HundunPerformanceSourceMetadata.cmake")
if(NOT DEFINED WORK)
  message(FATAL_ERROR "WORK fixture directory required")
endif()
file(MAKE_DIRECTORY "${WORK}")
string(REPEAT "a" 40 head)
string(REPEAT "B" 40 tree)
file(WRITE "${WORK}/.git_archival.txt"
  "HUNDUN_SOURCE_ARCHIVE_V1\ncommit=${head}\ntree=${tree}\n")
hundun_archive_source_metadata("${WORK}" actual_head actual_tree available)
string(TOLOWER "${tree}" tree_lower)
if(NOT available OR NOT actual_head STREQUAL head OR NOT actual_tree STREQUAL tree_lower)
  message(FATAL_ERROR "expanded archive identity")
endif()
foreach(text IN ITEMS
    "HUNDUN_SOURCE_ARCHIVE_V1\ncommit=$Format:%H$\ntree=$Format:%T$\n"
    "HUNDUN_SOURCE_ARCHIVE_V1\ncommit=ab\ntree=${tree}\n"
    "HUNDUN_SOURCE_ARCHIVE_V1\ncommit=${head}\ntree=zzzz\n"
    "HUNDUN_SOURCE_ARCHIVE_V2\ncommit=${head}\ntree=${tree}\n")
  file(WRITE "${WORK}/.git_archival.txt" "${text}")
  hundun_archive_source_metadata("${WORK}" actual_head actual_tree available)
  if(available OR NOT actual_head STREQUAL "unavailable" OR NOT actual_tree STREQUAL "unavailable")
    message(FATAL_ERROR "invalid archive identity accepted")
  endif()
endforeach()
file(REMOVE "${WORK}/.git_archival.txt")
hundun_archive_source_metadata("${WORK}" actual_head actual_tree available)
if(available)
  message(FATAL_ERROR "missing archive identity accepted")
endif()
message(STATUS "Archive provenance valid identity and five invalid/missing cases passed")
