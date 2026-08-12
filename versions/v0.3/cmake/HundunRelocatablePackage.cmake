# SPDX-License-Identifier: Apache-2.0

include_guard(GLOBAL)
include(GNUInstallDirs)

function(hundun_configure_relocatable_package package_root)
  get_filename_component(package_root "${package_root}" ABSOLUTE)
  set(CMAKE_INSTALL_RPATH "\$ORIGIN/../lib"
    PARENT_SCOPE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE PARENT_SCOPE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE PARENT_SCOPE)

  set(abi_manifest "${CMAKE_CURRENT_BINARY_DIR}/hundun-cantera-abi.json")
  file(WRITE "${abi_manifest}"
    "{\n"
    "  \"schema\": \"hundun.cantera.bundle.v1\",\n"
    "  \"bundle_generation\": \"cantera-3.2.0-gcc11-release-v4\",\n"
    "  \"library_sha256\": \"093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760\",\n"
    "  \"soname\": \"libcantera_shared.so.3\",\n"
    "  \"glibcxx_cxx11_abi\": 1,\n"
    "  \"relative_rpath\": \"$ORIGIN/../lib\"\n"
    "}\n")

  install(DIRECTORY "${package_root}/lib/"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    COMPONENT CanteraRuntime
    FILES_MATCHING PATTERN "libcantera_shared.so*")
  install(DIRECTORY "${package_root}/share/cantera/data/"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/cantera/data"
    COMPONENT CanteraRuntime)
  install(DIRECTORY "${package_root}/licenses/"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/cantera"
    COMPONENT CanteraRuntime)
  install(FILES
    "${PROJECT_SOURCE_DIR}/third_party/cantera/PREBUILT-LINUX-X86_64.json"
    "${abi_manifest}"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/hundun"
    COMPONENT CanteraRuntime)
  install(FILES "${PROJECT_SOURCE_DIR}/THIRD_PARTY_NOTICES"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/hundun"
    COMPONENT CanteraRuntime)
endfunction()

function(hundun_apply_relocatable_target target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "unknown relocatable target: ${target}")
  endif()
  set_target_properties("${target}" PROPERTIES
    INSTALL_RPATH "\$ORIGIN/../lib"
    BUILD_WITH_INSTALL_RPATH FALSE
    INSTALL_RPATH_USE_LINK_PATH FALSE)
endfunction()
