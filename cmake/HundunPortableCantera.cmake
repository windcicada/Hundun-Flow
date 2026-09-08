# SPDX-License-Identifier: Apache-2.0

include_guard(GLOBAL)

function(hundun_validate_cantera_consumer_text label contents)
  string(TOLOWER "${contents}" lower_contents)
  if(lower_contents MATCHES "fetchcontent")
    message(FATAL_ERROR
      "Cantera consumer policy rejected ${label}: FetchContent is forbidden")
  endif()
  if(lower_contents MATCHES "https?://")
    message(FATAL_ERROR
      "Cantera consumer policy rejected ${label}: network URLs are forbidden")
  endif()
  if(lower_contents MATCHES "build_bundled_cantera_linux_cpu[.]sh")
    message(FATAL_ERROR
      "Cantera consumer policy rejected ${label}: normal configure cannot invoke the maintainer builder")
  endif()
endfunction()

function(_hundun_require_cantera_path root relative)
  if(NOT EXISTS "${root}/${relative}")
    message(FATAL_ERROR
      "verified Cantera package is incomplete: missing ${relative}")
  endif()
endfunction()

function(hundun_configure_cantera_package package_root)
  string(CONCAT abi_flags " " "${CMAKE_CXX_FLAGS}" " "
    "${CMAKE_CXX_FLAGS_DEBUG}" " " "${CMAKE_CXX_FLAGS_RELEASE}")
  if(abi_flags MATCHES "-stdlib=libc[+][+]")
    message(FATAL_ERROR
      "verified Cantera package rejects Clang/libc++ mixing")
  endif()
  if(abi_flags MATCHES "_GLIBCXX_USE_CXX11_ABI(=| )0")
    message(FATAL_ERROR "verified Cantera package rejects ABI=0")
  endif()
  if(abi_flags MATCHES "_GLIBCXX_DEBUG|-fno-exceptions|-fno-rtti")
    message(FATAL_ERROR
      "verified Cantera package rejects debug-ABI or exceptions/RTTI mismatch")
  endif()
  if(package_root STREQUAL "")
    message(FATAL_ERROR
      "verified Cantera package is incomplete: HUNDUN_CANTERA_PACKAGE_ROOT is empty")
  endif()
  get_filename_component(package_root "${package_root}" ABSOLUTE)
  if(NOT IS_DIRECTORY "${package_root}")
    message(FATAL_ERROR
      "verified Cantera package is incomplete: root does not exist")
  endif()

  set(required_paths
    include/cantera/core.h
    include/cantera/base/Solution.h
    lib/libcantera_shared.so.3.2.0
    share/cantera/data/README.md
    licenses/Cantera/License.txt
    licenses/UPSTREAM-COMBINED-LICENSES.txt
    licenses/Boost/LICENSE_1_0.txt
    licenses/Eigen/COPYING.MPL2
    licenses/SUNDIALS/LICENSE
    licenses/SUNDIALS/NOTICE
    licenses/fmt/LICENSE.rst
    licenses/yaml-cpp/LICENSE)
  foreach(path IN LISTS required_paths)
    _hundun_require_cantera_path("${package_root}" "${path}")
  endforeach()

  file(SHA256 "${package_root}/lib/libcantera_shared.so.3.2.0"
    library_sha256)
  if(NOT library_sha256 STREQUAL
     "093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760")
    message(FATAL_ERROR
      "verified Cantera package is incomplete: shared-library SHA-256 mismatch")
  endif()
  file(SHA256 "${package_root}/licenses/Cantera/License.txt"
    license_sha256)
  if(NOT license_sha256 STREQUAL
     "e92980b9712ce20e73898a97b0116889e84e07f548d6be8591e87dcad79c41bb")
    message(FATAL_ERROR
      "verified Cantera package is incomplete: license SHA-256 mismatch")
  endif()

  foreach(link IN ITEMS libcantera_shared.so libcantera_shared.so.3)
    set(link_path "${package_root}/lib/${link}")
    if(NOT IS_SYMLINK "${link_path}")
      message(FATAL_ERROR
        "verified Cantera package is incomplete: ${link} is not a symlink")
    endif()
    file(READ_SYMLINK "${link_path}" target)
    if(NOT target STREQUAL "libcantera_shared.so.3.2.0")
      message(FATAL_ERROR
        "verified Cantera package is incomplete: ${link} escapes pinned library identity")
    endif()
  endforeach()

  if(NOT TARGET hundun_third_party_fmt)
    add_library(hundun_third_party_fmt INTERFACE IMPORTED GLOBAL)
    add_library(hundun_third_party_yamlcpp INTERFACE IMPORTED GLOBAL)
    add_library(hundun_third_party_sundials INTERFACE IMPORTED GLOBAL)
    add_library(hundun_third_party_eigen INTERFACE IMPORTED GLOBAL)
  endif()
  if(NOT TARGET hundun_third_party_cantera)
    add_library(hundun_third_party_cantera SHARED IMPORTED GLOBAL)
    set_target_properties(hundun_third_party_cantera PROPERTIES
      IMPORTED_LOCATION "${package_root}/lib/libcantera_shared.so.3.2.0"
      IMPORTED_SONAME "libcantera_shared.so.3"
      INTERFACE_INCLUDE_DIRECTORIES "${package_root}/include"
      INTERFACE_LINK_LIBRARIES
        "hundun_third_party_fmt;hundun_third_party_yamlcpp;hundun_third_party_sundials;hundun_third_party_eigen")
  endif()
  set(HUNDUN_CANTERA_DATA_DIR "${package_root}/share/cantera/data"
    CACHE PATH "Verified Cantera data directory" FORCE)
  set(HUNDUN_CANTERA_LICENSE_DIR "${package_root}/licenses"
    CACHE PATH "Verified Cantera license directory" FORCE)
endfunction()
