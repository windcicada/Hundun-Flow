# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

cmake_minimum_required(VERSION 3.21)

if(NOT DEFINED HUNDUN_SOURCE_DIR)
  message(FATAL_ERROR "HUNDUN_SOURCE_DIR is required")
endif()
get_filename_component(HUNDUN_SOURCE_DIR "${HUNDUN_SOURCE_DIR}" ABSOLUTE)
set(module "${HUNDUN_SOURCE_DIR}/cmake/HundunRelocatablePackage.cmake")
if(NOT EXISTS "${module}")
  message(FATAL_ERROR "Stage 4 relocatable package module is missing")
endif()
file(READ "${module}" module_text)
foreach(required IN ITEMS
    [[$ORIGIN/../lib]]
    [[CanteraRuntime]]
    [[PREBUILT-LINUX-X86_64.json]]
    [[THIRD_PARTY_NOTICES]])
  string(FIND "${module_text}" "${required}" found)
  if(found EQUAL -1)
    message(FATAL_ERROR "relocatable package contract missing ${required}")
  endif()
endforeach()
foreach(forbidden IN ITEMS
    [[libc.so]]
    [[libstdc++.so]]
    [[libgcc_s.so]]
    [[*.a]]
    [[mechanism]])
  string(FIND "${module_text}" "${forbidden}" found)
  if(NOT found EQUAL -1)
    message(FATAL_ERROR
      "relocatable package contract contains forbidden payload ${forbidden}")
  endif()
endforeach()

if(DEFINED HUNDUN_PACKAGE_ROOT)
  get_filename_component(HUNDUN_PACKAGE_ROOT "${HUNDUN_PACKAGE_ROOT}" ABSOLUTE)
  if(HUNDUN_PACKAGE_MUTATE_LIBRARY)
    file(APPEND
      "${HUNDUN_PACKAGE_ROOT}/lib/libcantera_shared.so.3.2.0"
      "single-library-upgrade-mutation")
  endif()
  foreach(required IN ITEMS
      lib/libcantera_shared.so.3.2.0
      lib/libcantera_shared.so
      lib/libcantera_shared.so.3
      share/cantera/data/README.md
      share/hundun/hundun-cantera-abi.json
      share/hundun/PREBUILT-LINUX-X86_64.json
      share/hundun/THIRD_PARTY_NOTICES)
    if(NOT EXISTS "${HUNDUN_PACKAGE_ROOT}/${required}")
      message(FATAL_ERROR "relocated package missing ${required}")
    endif()
  endforeach()
  file(SHA256
    "${HUNDUN_PACKAGE_ROOT}/lib/libcantera_shared.so.3.2.0"
    installed_library_sha256)
  if(NOT installed_library_sha256 STREQUAL
     "093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760")
    message(FATAL_ERROR
      "relocated package library does not match ABI manifest")
  endif()
  foreach(link IN ITEMS libcantera_shared.so libcantera_shared.so.3)
    if(NOT IS_SYMLINK "${HUNDUN_PACKAGE_ROOT}/lib/${link}")
      message(FATAL_ERROR "relocated package link ${link} is not a symlink")
    endif()
    file(READ_SYMLINK "${HUNDUN_PACKAGE_ROOT}/lib/${link}" target)
    if(NOT target STREQUAL "libcantera_shared.so.3.2.0")
      message(FATAL_ERROR "relocated package link ${link} changed generation")
    endif()
  endforeach()
  file(GLOB_RECURSE forbidden_payload LIST_DIRECTORIES false
    "${HUNDUN_PACKAGE_ROOT}/*.a"
    "${HUNDUN_PACKAGE_ROOT}/*.py"
    "${HUNDUN_PACKAGE_ROOT}/*.yaml"
    "${HUNDUN_PACKAGE_ROOT}/*.yml")
  if(forbidden_payload)
    message(FATAL_ERROR "relocated package contains forbidden payload")
  endif()
endif()
message(STATUS "Stage 4 relocatable package contract passed")
