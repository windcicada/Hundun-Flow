# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
cmake_minimum_required(VERSION 3.21)

foreach(required_variable IN ITEMS MANIFEST EXTERNAL_ROOT EXPECTED_CANTERA_SHA)
  if(NOT DEFINED ${required_variable} OR "${${required_variable}}" STREQUAL "")
    message(FATAL_ERROR "${required_variable} is required")
  endif()
endforeach()

if(NOT EXISTS "${MANIFEST}")
  message(FATAL_ERROR "Manifest does not exist: ${MANIFEST}")
endif()
if(NOT IS_DIRECTORY "${EXTERNAL_ROOT}")
  message(FATAL_ERROR "External root does not exist: ${EXTERNAL_ROOT}")
endif()
file(REAL_PATH "${EXTERNAL_ROOT}" canonical_external_root)
if(canonical_external_root STREQUAL "/")
  message(FATAL_ERROR "EXTERNAL_ROOT must not resolve to the filesystem root")
endif()
string(REGEX REPLACE "/+$" "" canonical_external_root
  "${canonical_external_root}")
set(canonical_external_root_prefix "${canonical_external_root}/")

file(READ "${MANIFEST}" manifest_json)

function(manifest_get output)
  string(JSON value ERROR_VARIABLE error GET "${manifest_json}" ${ARGN})
  if(NOT error STREQUAL "NOTFOUND")
    string(JOIN "." json_path ${ARGN})
    message(FATAL_ERROR "Missing or invalid manifest field ${json_path}: ${error}")
  endif()
  set(${output} "${value}" PARENT_SCOPE)
endfunction()

function(require_equal label actual expected)
  if(NOT "${actual}" STREQUAL "${expected}")
    message(FATAL_ERROR "${label}: expected '${expected}', got '${actual}'")
  endif()
endfunction()

function(verify_external_file relative_path expected_sha label)
  if("${relative_path}" STREQUAL "")
    message(FATAL_ERROR "${label} path is empty")
  endif()
  if(IS_ABSOLUTE "${relative_path}" OR
     "${relative_path}" MATCHES "(^|/)\\.\\.(/|$)" OR
     "${relative_path}" MATCHES "^file://")
    message(FATAL_ERROR "${label} path is not external-root relative: ${relative_path}")
  endif()
  set(absolute_path "${EXTERNAL_ROOT}/${relative_path}")
  if(NOT EXISTS "${absolute_path}" OR IS_DIRECTORY "${absolute_path}")
    message(FATAL_ERROR "${label} file is missing: ${absolute_path}")
  endif()
  file(REAL_PATH "${absolute_path}" resolved_path)
  string(FIND "${resolved_path}" "${canonical_external_root_prefix}"
    root_prefix_position)
  if(NOT root_prefix_position EQUAL 0)
    message(FATAL_ERROR
      "${label} resolves outside EXTERNAL_ROOT: ${relative_path} -> ${resolved_path}")
  endif()
  file(SHA256 "${resolved_path}" actual_sha)
  if(NOT actual_sha STREQUAL "${expected_sha}")
    message(FATAL_ERROR
      "${label} disk SHA mismatch: expected ${expected_sha}, got ${actual_sha}")
  endif()
endfunction()

function(verify_external_file_bytes relative_path expected_sha expected_bytes label)
  verify_external_file("${relative_path}" "${expected_sha}" "${label}")
  file(SIZE "${EXTERNAL_ROOT}/${relative_path}" actual_bytes)
  require_equal("${label} disk size mismatch" "${actual_bytes}"
    "${expected_bytes}")
endfunction()

function(verify_ruamel_wheel_tree tree_relative manifest_relative manifest_sha
    expected_file_count label)
  verify_external_file("${manifest_relative}" "${manifest_sha}"
    "${label} wheel tree manifest")

  set(tree_path "${EXTERNAL_ROOT}/${tree_relative}")
  if(NOT IS_DIRECTORY "${tree_path}")
    message(FATAL_ERROR "${label} wheel tree is missing: ${tree_path}")
  endif()
  file(REAL_PATH "${tree_path}" resolved_tree_path)
  string(FIND "${resolved_tree_path}" "${canonical_external_root_prefix}"
    tree_root_prefix_position)
  if(NOT tree_root_prefix_position EQUAL 0)
    message(FATAL_ERROR
      "${label} wheel tree resolves outside EXTERNAL_ROOT: ${tree_relative}")
  endif()
  string(REGEX REPLACE "/+$" "" resolved_tree_path "${resolved_tree_path}")
  set(tree_prefix "${resolved_tree_path}/")

  execute_process(
    COMMAND find -P "${resolved_tree_path}" -type l -print
    OUTPUT_VARIABLE symlink_output
    RESULT_VARIABLE symlink_result
    ERROR_VARIABLE symlink_error)
  if(NOT symlink_result EQUAL 0)
    message(FATAL_ERROR
      "${label} wheel tree symlink scan failed: ${symlink_error}")
  endif()
  string(STRIP "${symlink_output}" symlink_output)
  if(NOT symlink_output STREQUAL "")
    message(FATAL_ERROR
      "${label} wheel tree contains symlink(s): ${symlink_output}")
  endif()

  execute_process(
    COMMAND find -P "${resolved_tree_path}" -type f -printf "%P\\n"
    OUTPUT_VARIABLE actual_file_output
    RESULT_VARIABLE actual_file_result
    ERROR_VARIABLE actual_file_error)
  if(NOT actual_file_result EQUAL 0)
    message(FATAL_ERROR
      "${label} wheel tree file scan failed: ${actual_file_error}")
  endif()
  string(STRIP "${actual_file_output}" actual_file_output)
  if(actual_file_output STREQUAL "")
    set(actual_files)
  else()
    string(REPLACE "\n" ";" actual_files "${actual_file_output}")
    list(SORT actual_files)
  endif()

  set(manifest_path "${EXTERNAL_ROOT}/${manifest_relative}")
  file(STRINGS "${manifest_path}" manifest_lines)
  set(expected_files)
  foreach(manifest_line IN LISTS manifest_lines)
    if(NOT manifest_line MATCHES "^([0-9a-fA-F]+)[ \t]+(.+)$")
      message(FATAL_ERROR
        "${label} wheel tree manifest line is malformed: ${manifest_line}")
    endif()
    set(expected_sha "${CMAKE_MATCH_1}")
    set(expected_relative "${CMAKE_MATCH_2}")
    string(LENGTH "${expected_sha}" expected_sha_length)
    if(NOT expected_sha_length EQUAL 64)
      message(FATAL_ERROR
        "${label} wheel tree manifest SHA must be 64 hex characters:"
        " ${expected_sha}")
    endif()
    if(NOT expected_relative MATCHES "^\\./" OR
       expected_relative MATCHES "(^|/)\\.\\.(/|$)" OR
       expected_relative MATCHES "^/" OR
       expected_relative MATCHES "^file://")
      message(FATAL_ERROR
        "${label} wheel tree manifest path is invalid: ${expected_relative}")
    endif()
    string(REGEX REPLACE "^\\./" "" expected_relative
      "${expected_relative}")
    if(expected_relative STREQUAL "")
      message(FATAL_ERROR "${label} wheel tree manifest path is empty")
    endif()
    list(FIND expected_files "${expected_relative}" duplicate_index)
    if(NOT duplicate_index EQUAL -1)
      message(FATAL_ERROR
        "${label} wheel tree manifest contains duplicate path: ${expected_relative}")
    endif()
    list(APPEND expected_files "${expected_relative}")

    set(expected_path "${resolved_tree_path}/${expected_relative}")
    if(NOT EXISTS "${expected_path}" OR IS_DIRECTORY "${expected_path}")
      message(FATAL_ERROR
        "${label} wheel tree manifest file is missing: ${expected_relative}")
    endif()
    file(REAL_PATH "${expected_path}" resolved_expected_path)
    string(FIND "${resolved_expected_path}" "${tree_prefix}"
      expected_path_prefix_position)
    if(NOT expected_path_prefix_position EQUAL 0)
      message(FATAL_ERROR
        "${label} wheel tree manifest file escapes tree: ${expected_relative}")
    endif()
    file(SHA256 "${resolved_expected_path}" actual_sha)
    if(NOT actual_sha STREQUAL "${expected_sha}")
      message(FATAL_ERROR
        "${label} wheel tree file SHA mismatch for ${expected_relative}:"
        " expected ${expected_sha}, got ${actual_sha}")
    endif()
  endforeach()

  list(LENGTH expected_files expected_count)
  require_equal("${label} wheel tree manifest file count mismatch"
    "${expected_count}" "${expected_file_count}")
  list(LENGTH actual_files actual_count)
  require_equal("${label} wheel tree disk file count mismatch"
    "${actual_count}" "${expected_file_count}")
  foreach(expected_relative IN LISTS expected_files)
    list(FIND actual_files "${expected_relative}" actual_index)
    if(actual_index EQUAL -1)
      message(FATAL_ERROR
        "${label} wheel tree manifest path is not a regular disk file:"
        " ${expected_relative}")
    endif()
  endforeach()
  foreach(actual_relative IN LISTS actual_files)
    list(FIND expected_files "${actual_relative}" expected_index)
    if(expected_index EQUAL -1)
      message(FATAL_ERROR
        "${label} wheel tree has an unmanifested file: ${actual_relative}")
    endif()
  endforeach()
endfunction()

function(verify_metadata_identity relative_path expected_name expected_version
    expected_requires_python label)
  set(metadata_path "${EXTERNAL_ROOT}/${relative_path}")
  file(READ "${metadata_path}" metadata_text)
  string(FIND "${metadata_text}" "Name: ${expected_name}\n" name_position)
  if(name_position EQUAL -1)
    message(FATAL_ERROR "${label} METADATA Name field mismatch")
  endif()
  string(FIND "${metadata_text}" "Version: ${expected_version}\n"
    version_position)
  if(version_position EQUAL -1)
    message(FATAL_ERROR "${label} METADATA Version field mismatch")
  endif()
  string(FIND "${metadata_text}" "Requires-Python: ${expected_requires_python}\n"
    python_position)
  if(python_position EQUAL -1)
    message(FATAL_ERROR "${label} METADATA Requires-Python field mismatch")
  endif()
endfunction()

function(verify_pypi_metadata relative_path expected_name expected_version
    expected_license label)
  set(metadata_path "${EXTERNAL_ROOT}/${relative_path}")
  file(READ "${metadata_path}" pypi_json)
  string(JSON pypi_name ERROR_VARIABLE pypi_name_error
    GET "${pypi_json}" info name)
  if(NOT pypi_name_error STREQUAL "NOTFOUND")
    message(FATAL_ERROR "${label} PyPI metadata JSON name is invalid:"
      " ${pypi_name_error}")
  endif()
  string(JSON pypi_version ERROR_VARIABLE pypi_version_error
    GET "${pypi_json}" info version)
  if(NOT pypi_version_error STREQUAL "NOTFOUND")
    message(FATAL_ERROR "${label} PyPI metadata JSON version is invalid:"
      " ${pypi_version_error}")
  endif()
  string(JSON pypi_license ERROR_VARIABLE pypi_license_error
    GET "${pypi_json}" info license)
  if(NOT pypi_license_error STREQUAL "NOTFOUND")
    message(FATAL_ERROR "${label} PyPI metadata JSON license is invalid:"
      " ${pypi_license_error}")
  endif()
  require_equal("${label} PyPI metadata name mismatch" "${pypi_name}"
    "${expected_name}")
  require_equal("${label} PyPI metadata version mismatch" "${pypi_version}"
    "${expected_version}")
  require_equal("${label} PyPI metadata license mismatch" "${pypi_license}"
    "${expected_license}")
endfunction()

manifest_get(schema schema)
require_equal("Manifest schema mismatch" "${schema}" "hundun.stage4_p0.inputs.v1")

manifest_get(result result)
if(NOT result STREQUAL "PREFLIGHT_PARTIAL" AND
   NOT result STREQUAL "PREFLIGHT_PASS")
  message(FATAL_ERROR "Manifest result is invalid: ${result}")
endif()
manifest_get(stage4_product_accepted stage4_product_accepted)
if(stage4_product_accepted)
  message(FATAL_ERROR "stage4_product_accepted must remain false")
endif()
manifest_get(product_changes product_changes)
require_equal("Product-change boundary mismatch" "${product_changes}" "none")

manifest_get(target_os target os)
manifest_get(target_glibc target glibc_floor)
manifest_get(target_arch target arch)
manifest_get(target_compiler target compiler)
manifest_get(target_stdlib target stdlib)
manifest_get(target_cxx target cxx_standard)
manifest_get(target_abi target glibcxx_cxx11_abi)
manifest_get(target_isa target isa)
require_equal("Target OS mismatch" "${target_os}" "Ubuntu 22.04")
require_equal("Target glibc mismatch" "${target_glibc}" "2.35")
require_equal("Target architecture mismatch" "${target_arch}" "x86_64")
require_equal("Target compiler mismatch" "${target_compiler}" "GCC 11")
require_equal("Target standard library mismatch" "${target_stdlib}" "libstdc++")
require_equal("Target C++ standard mismatch" "${target_cxx}" "17")
require_equal("Target libstdc++ ABI mismatch" "${target_abi}" "1")
require_equal("Target ISA mismatch" "${target_isa}" "x86-64")

manifest_get(rootfs_status ubuntu_rootfs signature_status)
manifest_get(rootfs_fingerprint ubuntu_rootfs signing_key_fingerprint)
require_equal("Rootfs signature status mismatch" "${rootfs_status}" "verified_gpgv")
require_equal("Rootfs signing-key mismatch" "${rootfs_fingerprint}"
  "D2EB44626FDDC30B513D5BB71A5D6C4C7DB87C81")

foreach(rootfs_field IN ITEMS image signed_sums signature)
  if(rootfs_field STREQUAL "image")
    manifest_get(rootfs_path ubuntu_rootfs image_path)
    manifest_get(rootfs_sha ubuntu_rootfs image_sha256)
  elseif(rootfs_field STREQUAL "signed_sums")
    manifest_get(rootfs_path ubuntu_rootfs signed_sums_path)
    manifest_get(rootfs_sha ubuntu_rootfs signed_sums_sha256)
  else()
    manifest_get(rootfs_path ubuntu_rootfs signature_path)
    manifest_get(rootfs_sha ubuntu_rootfs signature_sha256)
  endif()
  verify_external_file("${rootfs_path}" "${rootfs_sha}" "Ubuntu ${rootfs_field}")
endforeach()

manifest_get(keyring_path ubuntu_rootfs keyring_path)
manifest_get(keyring_sha ubuntu_rootfs keyring_sha256)
set(expected_keyring_path "/usr/share/keyrings/ubuntu-cloudimage-keyring.gpg")
set(expected_keyring_sha
  "2ddbc33fdd3acfa0715914e3970b6a033faade0de25985eb17995b3aa85f455e")
require_equal("Ubuntu keyring path mismatch" "${keyring_path}"
  "${expected_keyring_path}")
require_equal("Ubuntu keyring manifest SHA mismatch" "${keyring_sha}"
  "${expected_keyring_sha}")
if(NOT EXISTS "${keyring_path}" OR IS_DIRECTORY "${keyring_path}")
  message(FATAL_ERROR "Ubuntu keyring file is missing: ${keyring_path}")
endif()
file(REAL_PATH "${keyring_path}" resolved_keyring_path)
file(REAL_PATH "${expected_keyring_path}" expected_resolved_keyring_path)
require_equal("Ubuntu keyring resolved path mismatch" "${resolved_keyring_path}"
  "${expected_resolved_keyring_path}")
file(SHA256 "${resolved_keyring_path}" actual_keyring_sha)
if(NOT actual_keyring_sha STREQUAL keyring_sha)
  message(FATAL_ERROR
    "Ubuntu keyring SHA mismatch: expected ${keyring_sha}, got ${actual_keyring_sha}")
endif()

manifest_get(cantera_tag cantera tag)
manifest_get(cantera_commit cantera commit)
manifest_get(cantera_manifest_sha cantera source_sha256)
require_equal("Cantera tag mismatch" "${cantera_tag}" "v3.2.0")
require_equal("Cantera commit mismatch" "${cantera_commit}"
  "4a8358eb80cfeb50474386b5f9ec0b3a83519889")
if(NOT cantera_manifest_sha STREQUAL EXPECTED_CANTERA_SHA)
  message(FATAL_ERROR "Cantera archive SHA mismatch")
endif()
manifest_get(cantera_path cantera source_path)
verify_external_file("${cantera_path}" "${cantera_manifest_sha}" "Cantera archive")

manifest_get(cantera_license_status cantera license status)
require_equal("Cantera license status mismatch" "${cantera_license_status}"
  "verified_redistributable")
manifest_get(cantera_license_path cantera license path)
manifest_get(cantera_license_sha cantera license sha256)
verify_external_file("${cantera_license_path}" "${cantera_license_sha}"
  "Cantera license")

string(JSON dependency_count ERROR_VARIABLE dependency_error
  LENGTH "${manifest_json}" dependencies)
if(NOT dependency_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR "Invalid dependency array: ${dependency_error}")
endif()
require_equal("Dependency count mismatch" "${dependency_count}" "4")

set(expected_names fmt yaml-cpp sundials eigen)
set(expected_commits
  a33701196adfad74917046096bf5a2aa0ab0bb50
  0579ae3d976091d7d664aa9d2527e0d0cff25763
  887af4374af2271db9310d31eaa9b5aeff49e829
  3147391d946bb4b6c68edd901f2add6ac1f31f8c)
set(expected_license_statuses
  verified_redistributable
  verified_redistributable
  verified_redistributable
  candidate_consumed_file_scope_audit_required)

math(EXPR dependency_last "${dependency_count} - 1")
foreach(index RANGE 0 ${dependency_last})
  list(GET expected_names ${index} expected_name)
  list(GET expected_commits ${index} expected_commit)
  list(GET expected_license_statuses ${index} expected_license_status)
  manifest_get(name dependencies ${index} name)
  manifest_get(commit dependencies ${index} commit)
  require_equal("Dependency ${index} name mismatch" "${name}" "${expected_name}")
  require_equal("${name} commit mismatch" "${commit}" "${expected_commit}")

  manifest_get(archive_path dependencies ${index} archive_path)
  manifest_get(archive_sha dependencies ${index} archive_sha256)
  verify_external_file("${archive_path}" "${archive_sha}" "${name} archive")

  manifest_get(license_status dependencies ${index} license_status)
  require_equal("${name} license status mismatch" "${license_status}"
    "${expected_license_status}")
  string(JSON license_count ERROR_VARIABLE license_error
    LENGTH "${manifest_json}" dependencies ${index} licenses)
  if(NOT license_error STREQUAL "NOTFOUND" OR license_count LESS 1)
    message(FATAL_ERROR "${name} license records are missing")
  endif()
  math(EXPR license_last "${license_count} - 1")
  foreach(license_index RANGE 0 ${license_last})
    manifest_get(license_spdx dependencies ${index} licenses ${license_index} spdx)
    manifest_get(license_path dependencies ${index} licenses ${license_index} path)
    manifest_get(license_sha dependencies ${index} licenses ${license_index} sha256)
    manifest_get(license_record_status
      dependencies ${index} licenses ${license_index} status)
    if(license_spdx STREQUAL "" OR license_record_status STREQUAL "")
      message(FATAL_ERROR "${name} license SPDX/status is missing")
    endif()
    verify_external_file("${license_path}" "${license_sha}"
      "${name} license ${license_index}")
  endforeach()

  if(name STREQUAL "eigen")
    manifest_get(eigen_consumption_status dependencies ${index}
      consumption_status)
    manifest_get(eigen_bundle_status dependencies ${index} bundle_status)
    require_equal("Eigen consumption status mismatch"
      "${eigen_consumption_status}" "required_source_candidate_for_p0_2")
    require_equal("Eigen bundle status mismatch" "${eigen_bundle_status}"
      "candidate_pending_p0_2_consumed_file_audit")

    set(required_eigen_obligation
      "compile P0-2 with EIGEN_MPL2_ONLY and treat any NonMPL2.h failure as a design blocker")
    string(JSON eigen_obligation_count ERROR_VARIABLE eigen_obligation_error
      LENGTH "${manifest_json}" dependencies ${index} obligations)
    if(NOT eigen_obligation_error STREQUAL "NOTFOUND" OR
       eigen_obligation_count LESS 1)
      message(FATAL_ERROR "Eigen obligations are missing")
    endif()
    set(has_required_eigen_obligation FALSE)
    math(EXPR eigen_obligation_last "${eigen_obligation_count} - 1")
    foreach(obligation_index RANGE 0 ${eigen_obligation_last})
      manifest_get(eigen_obligation dependencies ${index} obligations
        ${obligation_index})
      if(eigen_obligation STREQUAL required_eigen_obligation)
        set(has_required_eigen_obligation TRUE)
      endif()
    endforeach()
    if(NOT has_required_eigen_obligation)
      message(FATAL_ERROR "Eigen EIGEN_MPL2_ONLY obligation is missing")
    endif()
  endif()
endforeach()

string(JSON builder_header_dependency_count
  ERROR_VARIABLE builder_header_dependency_error
  LENGTH "${manifest_json}" builder_header_dependencies)
if(NOT builder_header_dependency_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR
    "Invalid builder-header dependency array: ${builder_header_dependency_error}")
endif()
require_equal("Builder-header dependency count mismatch"
  "${builder_header_dependency_count}" "1")

manifest_get(boost_name builder_header_dependencies 0 name)
manifest_get(boost_role builder_header_dependencies 0 role)
manifest_get(boost_package builder_header_dependencies 0 package)
manifest_get(boost_source_package builder_header_dependencies 0 source_package)
manifest_get(boost_version builder_header_dependencies 0 version)
manifest_get(boost_architecture builder_header_dependencies 0 architecture)
manifest_get(boost_repository builder_header_dependencies 0 repository)
manifest_get(boost_archive_url builder_header_dependencies 0 archive_url)
manifest_get(boost_archive_path builder_header_dependencies 0 archive_path)
manifest_get(boost_archive_sha builder_header_dependencies 0 archive_sha256)
manifest_get(boost_archive_bytes builder_header_dependencies 0 archive_bytes)
manifest_get(boost_copyright_holder
  builder_header_dependencies 0 copyright_holder)
manifest_get(boost_license_status
  builder_header_dependencies 0 license_status)
manifest_get(boost_consumption_status
  builder_header_dependencies 0 consumption_status)
manifest_get(boost_consumed_header_audit_status
  builder_header_dependencies 0 consumed_header_audit_status)
manifest_get(boost_bundle_status builder_header_dependencies 0 bundle_status)
manifest_get(boost_runtime_status builder_header_dependencies 0 runtime_status)

require_equal("Boost builder dependency name mismatch" "${boost_name}" "boost")
require_equal("Boost builder dependency role mismatch" "${boost_role}"
  "builder_header_only")
require_equal("Boost builder package mismatch" "${boost_package}"
  "libboost1.74-dev")
require_equal("Boost source package mismatch" "${boost_source_package}"
  "boost1.74")
require_equal("Boost version mismatch" "${boost_version}"
  "1.74.0-14ubuntu3")
require_equal("Boost architecture mismatch" "${boost_architecture}" "amd64")
require_equal("Boost repository mismatch" "${boost_repository}"
  "Canonical Ubuntu Jammy main")
require_equal("Boost archive URL mismatch" "${boost_archive_url}"
  "http://archive.ubuntu.com/ubuntu/pool/main/b/boost1.74/libboost1.74-dev_1.74.0-14ubuntu3_amd64.deb")
require_equal("Boost archive manifest SHA mismatch" "${boost_archive_sha}"
  "4d9c90e43f0d25db6280d1ee326771cbb76462f73b9430f06bac1de8d05b7a78")
require_equal("Boost archive manifest size mismatch" "${boost_archive_bytes}"
  "9608510")
require_equal("Boost copyright-holder mismatch" "${boost_copyright_holder}"
  "Boost project contributors and file-specific copyright holders")
require_equal("Boost license status mismatch" "${boost_license_status}"
  "candidate_consumed_header_scope_audit_required")
require_equal("Boost consumption status mismatch" "${boost_consumption_status}"
  "required_builder_header_only_for_p0_2")
require_equal("Boost consumed-header audit status mismatch"
  "${boost_consumed_header_audit_status}" "pending_p0_2")
require_equal("Boost bundle status mismatch" "${boost_bundle_status}"
  "forbidden_from_artifact_install")
require_equal("Boost runtime status mismatch" "${boost_runtime_status}"
  "forbidden_runtime_dependency")
verify_external_file("${boost_archive_path}" "${boost_archive_sha}"
  "Boost builder package")
file(SIZE "${EXTERNAL_ROOT}/${boost_archive_path}" boost_archive_actual_bytes)
require_equal("Boost archive disk size mismatch" "${boost_archive_actual_bytes}"
  "${boost_archive_bytes}")

string(JSON boost_license_count ERROR_VARIABLE boost_license_error
  LENGTH "${manifest_json}" builder_header_dependencies 0 licenses)
if(NOT boost_license_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR "Boost license array is invalid: ${boost_license_error}")
endif()
require_equal("Boost license count mismatch" "${boost_license_count}" "2")
manifest_get(boost_inventory_spdx
  builder_header_dependencies 0 licenses 0 spdx)
manifest_get(boost_inventory_url
  builder_header_dependencies 0 licenses 0 url)
manifest_get(boost_inventory_source_member
  builder_header_dependencies 0 licenses 0 source_member)
manifest_get(boost_inventory_path
  builder_header_dependencies 0 licenses 0 path)
manifest_get(boost_inventory_sha
  builder_header_dependencies 0 licenses 0 sha256)
manifest_get(boost_inventory_status
  builder_header_dependencies 0 licenses 0 status)
require_equal("Boost package inventory SPDX mismatch" "${boost_inventory_spdx}"
  "NOASSERTION")
require_equal("Boost package inventory URL mismatch" "${boost_inventory_url}"
  "http://archive.ubuntu.com/ubuntu/pool/main/b/boost1.74/libboost1.74-dev_1.74.0-14ubuntu3_amd64.deb")
require_equal("Boost package inventory member mismatch"
  "${boost_inventory_source_member}"
  "./usr/share/doc/libboost1.74-dev/copyright")
require_equal("Boost package inventory manifest SHA mismatch"
  "${boost_inventory_sha}"
  "17369eeac3938acb31085c8a1d4f1a40dc88a518a9389958a12c0592bb2d5766")
require_equal("Boost package inventory status mismatch"
  "${boost_inventory_status}" "verified_package_license_inventory")
verify_external_file("${boost_inventory_path}" "${boost_inventory_sha}"
  "Boost package license inventory")

manifest_get(boost_primary_spdx
  builder_header_dependencies 0 licenses 1 spdx)
manifest_get(boost_primary_url
  builder_header_dependencies 0 licenses 1 url)
manifest_get(boost_primary_path
  builder_header_dependencies 0 licenses 1 path)
manifest_get(boost_primary_sha
  builder_header_dependencies 0 licenses 1 sha256)
manifest_get(boost_primary_status
  builder_header_dependencies 0 licenses 1 status)
require_equal("Boost primary license SPDX mismatch" "${boost_primary_spdx}"
  "BSL-1.0")
require_equal("Boost primary license URL mismatch" "${boost_primary_url}"
  "https://www.boost.org/LICENSE_1_0.txt")
require_equal("Boost primary license manifest SHA mismatch"
  "${boost_primary_sha}"
  "c9bff75738922193e67fa726fa225535870d2aa1059f91452c411736284ad566")
require_equal("Boost primary license status mismatch" "${boost_primary_status}"
  "verified_primary_license")
verify_external_file("${boost_primary_path}" "${boost_primary_sha}"
  "Boost primary license")

set(expected_boost_obligations
  "preserve the exact Ubuntu package copyright inventory and official Boost Software License 1.0 evidence"
  "audit the actual Boost headers included or compiled before accepting the P0 artifact"
  "treat Boost as a builder-only header input and reject Boost headers from the P0 artifact install"
  "reject a Boost runtime dependency in the P0 artifact")
string(JSON boost_obligation_count ERROR_VARIABLE boost_obligation_error
  LENGTH "${manifest_json}" builder_header_dependencies 0 obligations)
if(NOT boost_obligation_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR "Boost obligation array is invalid: ${boost_obligation_error}")
endif()
require_equal("Boost obligation count mismatch" "${boost_obligation_count}" "4")
foreach(boost_obligation_index RANGE 0 3)
  list(GET expected_boost_obligations ${boost_obligation_index}
    expected_boost_obligation)
  manifest_get(boost_obligation builder_header_dependencies 0 obligations
    ${boost_obligation_index})
  require_equal("Boost obligation ${boost_obligation_index} mismatch"
    "${boost_obligation}" "${expected_boost_obligation}")
endforeach()

string(JSON boost_patch_count ERROR_VARIABLE boost_patch_error
  LENGTH "${manifest_json}" builder_header_dependencies 0 local_patches)
if(NOT boost_patch_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR "Boost local-patch array is invalid: ${boost_patch_error}")
endif()
require_equal("Boost local-patch count mismatch" "${boost_patch_count}" "0")

string(JSON builder_tool_dependency_count
  ERROR_VARIABLE builder_tool_dependency_error
  LENGTH "${manifest_json}" builder_tool_dependencies)
if(NOT builder_tool_dependency_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR
    "Invalid builder-tool dependency array: ${builder_tool_dependency_error}")
endif()
require_equal("Builder-tool dependency count mismatch"
  "${builder_tool_dependency_count}" "6")

set(expected_builder_tool_names
  doxygen
  libclang-cpp14
  libclang1-14
  libllvm14
  libxapian30
  libxml2)
set(expected_builder_tool_source_packages
  doxygen
  llvm-toolchain-14
  llvm-toolchain-14
  llvm-toolchain-14
  xapian-core
  libxml2)
set(expected_builder_tool_versions
  1.9.1-2ubuntu2
  1:14.0.0-1ubuntu1.1
  1:14.0.0-1ubuntu1.1
  1:14.0.0-1ubuntu1.1
  1.4.18-4
  2.9.13+dfsg-1ubuntu0.12)
set(expected_builder_tool_repositories
  "Canonical Ubuntu Jammy universe"
  "Canonical Ubuntu Jammy updates universe"
  "Canonical Ubuntu Jammy updates universe"
  "Canonical Ubuntu Jammy updates main"
  "Canonical Ubuntu Jammy universe"
  "Canonical Ubuntu Jammy updates main")
set(expected_builder_tool_urls
  "http://archive.ubuntu.com/ubuntu/pool/universe/d/doxygen/doxygen_1.9.1-2ubuntu2_amd64.deb"
  "http://archive.ubuntu.com/ubuntu/pool/universe/l/llvm-toolchain-14/libclang-cpp14_14.0.0-1ubuntu1.1_amd64.deb"
  "http://archive.ubuntu.com/ubuntu/pool/universe/l/llvm-toolchain-14/libclang1-14_14.0.0-1ubuntu1.1_amd64.deb"
  "http://archive.ubuntu.com/ubuntu/pool/main/l/llvm-toolchain-14/libllvm14_14.0.0-1ubuntu1.1_amd64.deb"
  "http://archive.ubuntu.com/ubuntu/pool/universe/x/xapian-core/libxapian30_1.4.18-4_amd64.deb"
  "http://archive.ubuntu.com/ubuntu/pool/main/libx/libxml2/libxml2_2.9.13%2bdfsg-1ubuntu0.12_amd64.deb")
set(expected_builder_tool_archive_paths
  "inputs/dependencies/doxygen-builder-v1/doxygen_1.9.1-2ubuntu2_amd64.deb"
  "inputs/dependencies/doxygen-builder-v1/libclang-cpp14_1%3a14.0.0-1ubuntu1.1_amd64.deb"
  "inputs/dependencies/doxygen-builder-v1/libclang1-14_1%3a14.0.0-1ubuntu1.1_amd64.deb"
  "inputs/dependencies/doxygen-builder-v1/libllvm14_1%3a14.0.0-1ubuntu1.1_amd64.deb"
  "inputs/dependencies/doxygen-builder-v1/libxapian30_1.4.18-4_amd64.deb"
  "inputs/dependencies/doxygen-builder-v1/libxml2_2.9.13+dfsg-1ubuntu0.12_amd64.deb")
set(expected_builder_tool_archive_shas
  d3fe7f77f505262db0872ba810b691feb4a607887ceb368e17a1da66ee0387f7
  462546f6149fbef99b47f3e1e01e0743d256ae279337c7f89a80a976dfd02175
  b75b743f5d5effaab97790c1379fb1855d1a20bd5432a2387bf4dc82d86d45e3
  9044b614a6c7fb6262e7cbeb13dc731fc0c92bed96281c1a3920dd706442ee8e
  5cfe52f4ca570e85efa828efda6d6831ceb0f667a32faf5438887cfaf528b7c2
  b3678e6e4b166bc0e4226fb118d489ab51802c772914421397c9dcb2dd0e0d2b)
set(expected_builder_tool_archive_bytes
  4620198
  12053266
  6792182
  23967046
  700928
  764660)
set(expected_builder_tool_holders
  "Doxygen contributors and package file-specific copyright holders"
  "LLVM contributors and package file-specific copyright holders"
  "LLVM contributors and package file-specific copyright holders"
  "LLVM contributors and package file-specific copyright holders"
  "Xapian contributors and package file-specific copyright holders"
  "libxml2 contributors and package file-specific copyright holders")
set(expected_builder_tool_copyright_members
  "./usr/share/doc/doxygen/copyright"
  "./usr/share/doc/libclang-cpp14/copyright"
  "./usr/share/doc/libclang1-14/copyright"
  "./usr/share/doc/libllvm14/copyright"
  "./usr/share/doc/libxapian30/copyright"
  "./usr/share/doc/libxml2/copyright")
set(expected_builder_tool_copyright_paths
  "logs/doxygen-builder-v1/extract/doxygen/usr/share/doc/doxygen/copyright"
  "logs/doxygen-builder-v1/extract/libclang-cpp14/usr/share/doc/libclang-cpp14/copyright"
  "logs/doxygen-builder-v1/extract/libclang1-14/usr/share/doc/libclang1-14/copyright"
  "logs/doxygen-builder-v1/extract/libllvm14/usr/share/doc/libllvm14/copyright"
  "logs/doxygen-builder-v1/extract/libxapian30/usr/share/doc/libxapian30/copyright"
  "logs/doxygen-builder-v1/extract/libxml2/usr/share/doc/libxml2/copyright")
set(expected_builder_tool_copyright_shas
  286b6aadd5b010bb3df2e226dc67fa8cf5819af8870b3ac093442370214c50b9
  dc2a171c75ca72818f7f4312602f27cb0e4bf6d010a86e7073a3797cb228ddff
  dc2a171c75ca72818f7f4312602f27cb0e4bf6d010a86e7073a3797cb228ddff
  dc2a171c75ca72818f7f4312602f27cb0e4bf6d010a86e7073a3797cb228ddff
  d530c1aa427b7e55f170f48b43f3657f847255438496847796eb8a0496bdd5d6
  ee746b96cfa5be73c3ea3e4cfb1285e9b315d4c9267f99b2ee9c5d911d9fe3f4)
set(expected_builder_tool_runtime_statuses
  forbidden_runtime_dependency
  forbidden_runtime_dependency
  forbidden_runtime_dependency
  forbidden_runtime_dependency
  forbidden_runtime_dependency
  audit_pending_p0_2)

foreach(builder_tool_index RANGE 0 5)
  list(GET expected_builder_tool_names ${builder_tool_index}
    expected_builder_tool_name)
  list(GET expected_builder_tool_source_packages ${builder_tool_index}
    expected_builder_tool_source_package)
  list(GET expected_builder_tool_versions ${builder_tool_index}
    expected_builder_tool_version)
  list(GET expected_builder_tool_repositories ${builder_tool_index}
    expected_builder_tool_repository)
  list(GET expected_builder_tool_urls ${builder_tool_index}
    expected_builder_tool_url)
  list(GET expected_builder_tool_archive_paths ${builder_tool_index}
    expected_builder_tool_archive_path)
  list(GET expected_builder_tool_archive_shas ${builder_tool_index}
    expected_builder_tool_archive_sha)
  list(GET expected_builder_tool_archive_bytes ${builder_tool_index}
    expected_builder_tool_archive_size)
  list(GET expected_builder_tool_holders ${builder_tool_index}
    expected_builder_tool_holder)
  list(GET expected_builder_tool_copyright_members ${builder_tool_index}
    expected_builder_tool_copyright_member)
  list(GET expected_builder_tool_copyright_paths ${builder_tool_index}
    expected_builder_tool_copyright_path)
  list(GET expected_builder_tool_copyright_shas ${builder_tool_index}
    expected_builder_tool_copyright_sha)
  list(GET expected_builder_tool_runtime_statuses ${builder_tool_index}
    expected_builder_tool_runtime_status)

  manifest_get(builder_tool_name
    builder_tool_dependencies ${builder_tool_index} name)
  manifest_get(builder_tool_role
    builder_tool_dependencies ${builder_tool_index} role)
  manifest_get(builder_tool_package
    builder_tool_dependencies ${builder_tool_index} package)
  manifest_get(builder_tool_source_package
    builder_tool_dependencies ${builder_tool_index} source_package)
  manifest_get(builder_tool_version
    builder_tool_dependencies ${builder_tool_index} version)
  manifest_get(builder_tool_architecture
    builder_tool_dependencies ${builder_tool_index} architecture)
  manifest_get(builder_tool_repository
    builder_tool_dependencies ${builder_tool_index} repository)
  manifest_get(builder_tool_archive_url
    builder_tool_dependencies ${builder_tool_index} archive_url)
  manifest_get(builder_tool_archive_path
    builder_tool_dependencies ${builder_tool_index} archive_path)
  manifest_get(builder_tool_archive_sha
    builder_tool_dependencies ${builder_tool_index} archive_sha256)
  manifest_get(builder_tool_archive_bytes
    builder_tool_dependencies ${builder_tool_index} archive_bytes)
  manifest_get(builder_tool_copyright_holder
    builder_tool_dependencies ${builder_tool_index} copyright_holder)
  manifest_get(builder_tool_license_status
    builder_tool_dependencies ${builder_tool_index} license_status)
  manifest_get(builder_tool_inventory_spdx
    builder_tool_dependencies ${builder_tool_index} copyright_inventory spdx)
  manifest_get(builder_tool_inventory_member
    builder_tool_dependencies ${builder_tool_index} copyright_inventory source_member)
  manifest_get(builder_tool_inventory_path
    builder_tool_dependencies ${builder_tool_index} copyright_inventory path)
  manifest_get(builder_tool_inventory_sha
    builder_tool_dependencies ${builder_tool_index} copyright_inventory sha256)
  manifest_get(builder_tool_inventory_status
    builder_tool_dependencies ${builder_tool_index} copyright_inventory status)
  manifest_get(builder_tool_consumption_status
    builder_tool_dependencies ${builder_tool_index} consumption_status)
  manifest_get(builder_tool_artifact_audit_status
    builder_tool_dependencies ${builder_tool_index} artifact_audit_status)
  manifest_get(builder_tool_bundle_status
    builder_tool_dependencies ${builder_tool_index} bundle_status)
  manifest_get(builder_tool_runtime_status
    builder_tool_dependencies ${builder_tool_index} runtime_status)

  require_equal("Builder tool ${builder_tool_index} name mismatch"
    "${builder_tool_name}" "${expected_builder_tool_name}")
  require_equal("${builder_tool_name} builder dependency role mismatch"
    "${builder_tool_role}" "builder_tool_only")
  require_equal("${builder_tool_name} builder package mismatch"
    "${builder_tool_package}" "${expected_builder_tool_name}")
  require_equal("${builder_tool_name} source package mismatch"
    "${builder_tool_source_package}" "${expected_builder_tool_source_package}")
  require_equal("${builder_tool_name} version mismatch"
    "${builder_tool_version}" "${expected_builder_tool_version}")
  require_equal("${builder_tool_name} architecture mismatch"
    "${builder_tool_architecture}" "amd64")
  require_equal("${builder_tool_name} repository mismatch"
    "${builder_tool_repository}" "${expected_builder_tool_repository}")
  require_equal("${builder_tool_name} archive URL mismatch"
    "${builder_tool_archive_url}" "${expected_builder_tool_url}")
  require_equal("${builder_tool_name} archive path mismatch"
    "${builder_tool_archive_path}" "${expected_builder_tool_archive_path}")
  require_equal("${builder_tool_name} archive manifest SHA mismatch"
    "${builder_tool_archive_sha}" "${expected_builder_tool_archive_sha}")
  require_equal("${builder_tool_name} archive manifest size mismatch"
    "${builder_tool_archive_bytes}" "${expected_builder_tool_archive_size}")
  require_equal("${builder_tool_name} copyright-holder mismatch"
    "${builder_tool_copyright_holder}" "${expected_builder_tool_holder}")
  require_equal("${builder_tool_name} license status mismatch"
    "${builder_tool_license_status}"
    "verified_package_multi_license_inventory")
  require_equal("${builder_tool_name} copyright inventory SPDX mismatch"
    "${builder_tool_inventory_spdx}" "NOASSERTION")
  require_equal("${builder_tool_name} copyright inventory member mismatch"
    "${builder_tool_inventory_member}" "${expected_builder_tool_copyright_member}")
  require_equal("${builder_tool_name} copyright inventory path mismatch"
    "${builder_tool_inventory_path}" "${expected_builder_tool_copyright_path}")
  require_equal("${builder_tool_name} copyright inventory manifest SHA mismatch"
    "${builder_tool_inventory_sha}" "${expected_builder_tool_copyright_sha}")
  require_equal("${builder_tool_name} copyright inventory status mismatch"
    "${builder_tool_inventory_status}" "verified_package_license_inventory")
  require_equal("${builder_tool_name} consumption status mismatch"
    "${builder_tool_consumption_status}"
    "required_builder_tool_only_for_p0_2")
  require_equal("${builder_tool_name} artifact audit status mismatch"
    "${builder_tool_artifact_audit_status}" "pending_p0_2")
  require_equal("${builder_tool_name} bundle status mismatch"
    "${builder_tool_bundle_status}" "forbidden_from_artifact_install")
  require_equal("${builder_tool_name} runtime status mismatch"
    "${builder_tool_runtime_status}" "${expected_builder_tool_runtime_status}")

  verify_external_file("${builder_tool_archive_path}"
    "${builder_tool_archive_sha}" "${builder_tool_name} builder package")
  file(SIZE "${EXTERNAL_ROOT}/${builder_tool_archive_path}"
    builder_tool_archive_actual_bytes)
  require_equal("${builder_tool_name} archive disk size mismatch"
    "${builder_tool_archive_actual_bytes}" "${builder_tool_archive_bytes}")
  verify_external_file("${builder_tool_inventory_path}"
    "${builder_tool_inventory_sha}"
    "${builder_tool_name} package copyright inventory")

  string(JSON builder_tool_patch_count ERROR_VARIABLE builder_tool_patch_error
    LENGTH "${manifest_json}" builder_tool_dependencies ${builder_tool_index}
    local_patches)
  if(NOT builder_tool_patch_error STREQUAL "NOTFOUND")
    message(FATAL_ERROR
      "${builder_tool_name} local-patch array is invalid: ${builder_tool_patch_error}")
  endif()
  require_equal("${builder_tool_name} local-patch count mismatch"
    "${builder_tool_patch_count}" "0")
endforeach()

string(JSON builder_python_dependency_count
  ERROR_VARIABLE builder_python_dependency_error
  LENGTH "${manifest_json}" builder_python_dependencies)
if(NOT builder_python_dependency_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR
    "Invalid builder-Python dependency array: ${builder_python_dependency_error}")
endif()
require_equal("Builder-Python dependency count mismatch"
  "${builder_python_dependency_count}" "3")

manifest_get(typing_name builder_python_dependencies 0 name)
manifest_get(typing_distribution builder_python_dependencies 0 distribution)
manifest_get(typing_role builder_python_dependencies 0 role)
manifest_get(typing_version builder_python_dependencies 0 version)
manifest_get(typing_repository builder_python_dependencies 0 repository)
manifest_get(typing_project_url builder_python_dependencies 0 project_url)
manifest_get(typing_metadata_url builder_python_dependencies 0 metadata_url)
manifest_get(typing_metadata_path builder_python_dependencies 0 metadata_path)
manifest_get(typing_metadata_sha builder_python_dependencies 0 metadata_sha256)
manifest_get(typing_wheel_url builder_python_dependencies 0 wheel_url)
manifest_get(typing_wheel_path builder_python_dependencies 0 wheel_path)
manifest_get(typing_wheel_sha builder_python_dependencies 0 wheel_sha256)
manifest_get(typing_wheel_bytes builder_python_dependencies 0 wheel_bytes)
manifest_get(typing_sdist_url builder_python_dependencies 0 sdist_url)
manifest_get(typing_sdist_path builder_python_dependencies 0 sdist_path)
manifest_get(typing_sdist_sha builder_python_dependencies 0 sdist_sha256)
manifest_get(typing_sdist_bytes builder_python_dependencies 0 sdist_bytes)
manifest_get(typing_python_requires builder_python_dependencies 0 python_requires)
manifest_get(typing_copyright_holder
  builder_python_dependencies 0 copyright_holder)
manifest_get(typing_api_evidence_path
  builder_python_dependencies 0 api_evidence path)
manifest_get(typing_api_evidence_sha
  builder_python_dependencies 0 api_evidence sha256)
manifest_get(typing_api_evidence_status
  builder_python_dependencies 0 api_evidence status)
manifest_get(typing_license_spdx builder_python_dependencies 0 license spdx)
manifest_get(typing_wheel_license_path
  builder_python_dependencies 0 license wheel_path)
manifest_get(typing_wheel_license_sha
  builder_python_dependencies 0 license wheel_sha256)
manifest_get(typing_sdist_license_path
  builder_python_dependencies 0 license sdist_path)
manifest_get(typing_sdist_license_sha
  builder_python_dependencies 0 license sdist_sha256)
manifest_get(typing_license_status builder_python_dependencies 0 license status)
manifest_get(typing_consumption_status
  builder_python_dependencies 0 consumption_status)
manifest_get(typing_injection builder_python_dependencies 0 injection)
manifest_get(typing_rootfs_install_status
  builder_python_dependencies 0 rootfs_install_status)
manifest_get(typing_artifact_audit_status
  builder_python_dependencies 0 artifact_audit_status)
manifest_get(typing_bundle_status builder_python_dependencies 0 bundle_status)
manifest_get(typing_runtime_status builder_python_dependencies 0 runtime_status)

require_equal("typing_extensions name mismatch" "${typing_name}"
  "typing_extensions")
require_equal("typing_extensions distribution mismatch" "${typing_distribution}"
  "typing-extensions")
require_equal("typing_extensions role mismatch" "${typing_role}"
  "builder_pythonpath_only")
require_equal("typing_extensions version mismatch" "${typing_version}" "4.15.0")
require_equal("typing_extensions repository mismatch" "${typing_repository}"
  "Python Package Index")
require_equal("typing_extensions project URL mismatch" "${typing_project_url}"
  "https://pypi.org/project/typing-extensions/4.15.0/")
require_equal("typing_extensions metadata URL mismatch" "${typing_metadata_url}"
  "https://pypi.org/pypi/typing-extensions/4.15.0/json")
require_equal("typing_extensions metadata SHA mismatch" "${typing_metadata_sha}"
  "e97e0b1087254aa1c7e8b2074c3796124dfd7d26e0f54ffcdc3a975b53047938")
require_equal("typing_extensions wheel URL mismatch" "${typing_wheel_url}"
  "https://files.pythonhosted.org/packages/18/67/36e9267722cc04a6b9f15c7f3441c2363321a3ea07da7ae0c0707beb2a9c/typing_extensions-4.15.0-py3-none-any.whl")
require_equal("typing_extensions wheel manifest SHA mismatch" "${typing_wheel_sha}"
  "f0fa19c6845758ab08074a0cfa8b7aecb71c999ca73d62883bc25cc018c4e548")
require_equal("typing_extensions wheel manifest size mismatch"
  "${typing_wheel_bytes}" "44614")
require_equal("typing_extensions sdist URL mismatch" "${typing_sdist_url}"
  "https://files.pythonhosted.org/packages/72/94/1a15dd82efb362ac84269196e94cf00f187f7ed21c242792a923cdb1c61f/typing_extensions-4.15.0.tar.gz")
require_equal("typing_extensions sdist manifest SHA mismatch" "${typing_sdist_sha}"
  "0cea48d173cc12fa28ecabc3b837ea3cf6f38c6d1136f85cbaaf598984861466")
require_equal("typing_extensions sdist manifest size mismatch"
  "${typing_sdist_bytes}" "109391")
require_equal("typing_extensions Python requirement mismatch"
  "${typing_python_requires}" ">=3.9")
require_equal("typing_extensions copyright-holder mismatch"
  "${typing_copyright_holder}"
  "Python Software Foundation and historical licensors named in LICENSE")
require_equal("typing_extensions API-evidence SHA mismatch"
  "${typing_api_evidence_sha}"
  "f61a6f3540b43f1545d42866718b3dc19b674afa4745d873d7e2827152409c8b")
require_equal("typing_extensions API-evidence status mismatch"
  "${typing_api_evidence_status}" "self_runtime_support_added_in_4_0_0")
require_equal("typing_extensions license SPDX mismatch" "${typing_license_spdx}"
  "PSF-2.0")
require_equal("typing_extensions wheel license SHA mismatch"
  "${typing_wheel_license_sha}"
  "3b2f81fe21d181c499c59a256c8e1968455d6689d269aa85373bfb6af41da3bf")
require_equal("typing_extensions sdist license SHA mismatch"
  "${typing_sdist_license_sha}"
  "3b2f81fe21d181c499c59a256c8e1968455d6689d269aa85373bfb6af41da3bf")
require_equal("typing_extensions license status mismatch" "${typing_license_status}"
  "verified_matching_wheel_and_sdist")
require_equal("typing_extensions consumption status mismatch"
  "${typing_consumption_status}"
  "required_cantera_sourcegen_python310_for_p0_2")
require_equal("typing_extensions injection mismatch" "${typing_injection}"
  "exact_wheel_via_pythonpath")
require_equal("typing_extensions rootfs install status mismatch"
  "${typing_rootfs_install_status}" "forbidden")
require_equal("typing_extensions artifact audit status mismatch"
  "${typing_artifact_audit_status}" "pending_p0_2")
require_equal("typing_extensions bundle status mismatch" "${typing_bundle_status}"
  "forbidden_from_artifact_install")
require_equal("typing_extensions runtime status mismatch" "${typing_runtime_status}"
  "forbidden_runtime_dependency")

verify_external_file("${typing_metadata_path}" "${typing_metadata_sha}"
  "typing_extensions PyPI metadata")
verify_external_file("${typing_wheel_path}" "${typing_wheel_sha}"
  "typing_extensions wheel")
file(SIZE "${EXTERNAL_ROOT}/${typing_wheel_path}" typing_wheel_actual_bytes)
require_equal("typing_extensions wheel disk size mismatch"
  "${typing_wheel_actual_bytes}" "${typing_wheel_bytes}")
verify_external_file("${typing_sdist_path}" "${typing_sdist_sha}"
  "typing_extensions sdist")
file(SIZE "${EXTERNAL_ROOT}/${typing_sdist_path}" typing_sdist_actual_bytes)
require_equal("typing_extensions sdist disk size mismatch"
  "${typing_sdist_actual_bytes}" "${typing_sdist_bytes}")
verify_external_file("${typing_api_evidence_path}" "${typing_api_evidence_sha}"
  "typing_extensions API evidence")
verify_external_file("${typing_wheel_license_path}" "${typing_wheel_license_sha}"
  "typing_extensions wheel license")
verify_external_file("${typing_sdist_license_path}" "${typing_sdist_license_sha}"
  "typing_extensions sdist license")

string(JSON typing_patch_count ERROR_VARIABLE typing_patch_error
  LENGTH "${manifest_json}" builder_python_dependencies 0 local_patches)
if(NOT typing_patch_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR
    "typing_extensions local-patch array is invalid: ${typing_patch_error}")
endif()
require_equal("typing_extensions local-patch count mismatch"
  "${typing_patch_count}" "0")

manifest_get(ruamel_yaml_name builder_python_dependencies 1 name)
manifest_get(ruamel_yaml_distribution builder_python_dependencies 1 distribution)
manifest_get(ruamel_yaml_role builder_python_dependencies 1 role)
manifest_get(ruamel_yaml_version builder_python_dependencies 1 version)
manifest_get(ruamel_yaml_repository builder_python_dependencies 1 repository)
manifest_get(ruamel_yaml_project_url builder_python_dependencies 1 project_url)
manifest_get(ruamel_yaml_metadata_url builder_python_dependencies 1 metadata_url)
manifest_get(ruamel_yaml_metadata_path builder_python_dependencies 1 metadata_path)
manifest_get(ruamel_yaml_metadata_sha builder_python_dependencies 1 metadata_sha256)
manifest_get(ruamel_yaml_metadata_bytes builder_python_dependencies 1 metadata_bytes)
manifest_get(ruamel_yaml_wheel_url builder_python_dependencies 1 wheel_url)
manifest_get(ruamel_yaml_wheel_path builder_python_dependencies 1 wheel_path)
manifest_get(ruamel_yaml_wheel_sha builder_python_dependencies 1 wheel_sha256)
manifest_get(ruamel_yaml_wheel_bytes builder_python_dependencies 1 wheel_bytes)
manifest_get(ruamel_yaml_wheel_tag builder_python_dependencies 1 wheel_tag)
manifest_get(ruamel_yaml_wheel_metadata_path builder_python_dependencies 1
  wheel_metadata_path)
manifest_get(ruamel_yaml_wheel_metadata_sha builder_python_dependencies 1
  wheel_metadata_sha256)
manifest_get(ruamel_yaml_wheel_tree_path builder_python_dependencies 1
  wheel_tree_path)
manifest_get(ruamel_yaml_wheel_tree_manifest_path builder_python_dependencies 1
  wheel_tree_manifest_path)
manifest_get(ruamel_yaml_wheel_tree_manifest_sha builder_python_dependencies 1
  wheel_tree_manifest_sha256)
manifest_get(ruamel_yaml_wheel_tree_file_count builder_python_dependencies 1
  wheel_tree_file_count)
manifest_get(ruamel_yaml_sdist_url builder_python_dependencies 1 sdist_url)
manifest_get(ruamel_yaml_sdist_path builder_python_dependencies 1 sdist_path)
manifest_get(ruamel_yaml_sdist_sha builder_python_dependencies 1 sdist_sha256)
manifest_get(ruamel_yaml_sdist_bytes builder_python_dependencies 1 sdist_bytes)
manifest_get(ruamel_yaml_python_requires builder_python_dependencies 1
  python_requires)
manifest_get(ruamel_yaml_requires_dist builder_python_dependencies 1
  requires_dist)
manifest_get(ruamel_yaml_copyright_holder builder_python_dependencies 1
  copyright_holder)
manifest_get(ruamel_yaml_license_spdx builder_python_dependencies 1 license spdx)
manifest_get(ruamel_yaml_wheel_license_path builder_python_dependencies 1
  license wheel_path)
manifest_get(ruamel_yaml_wheel_license_sha builder_python_dependencies 1
  license wheel_sha256)
manifest_get(ruamel_yaml_sdist_license_path builder_python_dependencies 1
  license sdist_path)
manifest_get(ruamel_yaml_sdist_license_sha builder_python_dependencies 1
  license sdist_sha256)
manifest_get(ruamel_yaml_license_status builder_python_dependencies 1
  license status)
manifest_get(ruamel_yaml_consumption_status builder_python_dependencies 1
  consumption_status)
manifest_get(ruamel_yaml_injection builder_python_dependencies 1 injection)
manifest_get(ruamel_yaml_rootfs_status builder_python_dependencies 1
  rootfs_install_status)
manifest_get(ruamel_yaml_artifact_status builder_python_dependencies 1
  artifact_audit_status)
manifest_get(ruamel_yaml_bundle_status builder_python_dependencies 1 bundle_status)
manifest_get(ruamel_yaml_runtime_status builder_python_dependencies 1
  runtime_status)

require_equal("ruamel.yaml name mismatch" "${ruamel_yaml_name}" "ruamel.yaml")
require_equal("ruamel.yaml distribution mismatch" "${ruamel_yaml_distribution}"
  "ruamel.yaml")
require_equal("ruamel.yaml role mismatch" "${ruamel_yaml_role}"
  "builder_pythonpath_only")
require_equal("ruamel.yaml version mismatch" "${ruamel_yaml_version}" "0.17.21")
require_equal("ruamel.yaml repository mismatch" "${ruamel_yaml_repository}"
  "Python Package Index")
require_equal("ruamel.yaml project URL mismatch" "${ruamel_yaml_project_url}"
  "https://pypi.org/project/ruamel.yaml/0.17.21/")
require_equal("ruamel.yaml metadata URL mismatch" "${ruamel_yaml_metadata_url}"
  "https://pypi.org/pypi/ruamel.yaml/0.17.21/json")
require_equal("ruamel.yaml metadata path mismatch" "${ruamel_yaml_metadata_path}"
  "inputs/dependencies/ruamel-yaml-0.17.21/pypi-metadata.json")
require_equal("ruamel.yaml metadata SHA mismatch" "${ruamel_yaml_metadata_sha}"
  "e140839f8da85e5ba08ae2ce1d875bca7b4797a780dbf85c4319633ad33540c3")
require_equal("ruamel.yaml metadata size mismatch" "${ruamel_yaml_metadata_bytes}"
  "15193")
require_equal("ruamel.yaml wheel URL mismatch" "${ruamel_yaml_wheel_url}"
  "https://files.pythonhosted.org/packages/9e/cb/938214ac358fbef7058343b3765c79a1b7ed0c366f7f992ce7ff38335652/ruamel.yaml-0.17.21-py3-none-any.whl")
require_equal("ruamel.yaml wheel path mismatch" "${ruamel_yaml_wheel_path}"
  "inputs/dependencies/ruamel-yaml-0.17.21/ruamel.yaml-0.17.21-py3-none-any.whl")
require_equal("ruamel.yaml wheel SHA mismatch" "${ruamel_yaml_wheel_sha}"
  "742b35d3d665023981bd6d16b3d24248ce5df75fdb4e2924e93a05c1f8b61ca7")
require_equal("ruamel.yaml wheel size mismatch" "${ruamel_yaml_wheel_bytes}"
  "109478")
require_equal("ruamel.yaml wheel tag mismatch" "${ruamel_yaml_wheel_tag}"
  "py3-none-any")
require_equal("ruamel.yaml wheel METADATA SHA mismatch"
  "${ruamel_yaml_wheel_metadata_sha}"
  "10efd31d77a2b726f57cfc35ac1e2c85db47f49b5adb9f2f290da9f55d99f0ad")
require_equal("ruamel.yaml wheel METADATA path mismatch"
  "${ruamel_yaml_wheel_metadata_path}"
  "inputs/dependencies/ruamel-yaml-0.17.21/wheel-extract/ruamel.yaml-0.17.21.dist-info/METADATA")
require_equal("ruamel.yaml wheel tree path mismatch"
  "${ruamel_yaml_wheel_tree_path}"
  "inputs/dependencies/ruamel-yaml-0.17.21/wheel-extract")
require_equal("ruamel.yaml wheel tree manifest path mismatch"
  "${ruamel_yaml_wheel_tree_manifest_path}"
  "manifests/ruamel-yaml-0.17.21-wheel-extract.sha256")
require_equal("ruamel.yaml wheel tree manifest SHA mismatch"
  "${ruamel_yaml_wheel_tree_manifest_sha}"
  "ce5fef34f3fb24dd65ee21902b0d5c2b33123333c3ed95b6bd6241644702460f")
require_equal("ruamel.yaml wheel tree file count mismatch"
  "${ruamel_yaml_wheel_tree_file_count}" "36")
require_equal("ruamel.yaml sdist URL mismatch" "${ruamel_yaml_sdist_url}"
  "https://files.pythonhosted.org/packages/46/a9/6ed24832095b692a8cecc323230ce2ec3480015fbfa4b79941bd41b23a3c/ruamel.yaml-0.17.21.tar.gz")
require_equal("ruamel.yaml sdist path mismatch" "${ruamel_yaml_sdist_path}"
  "inputs/dependencies/ruamel-yaml-0.17.21/ruamel.yaml-0.17.21.tar.gz")
require_equal("ruamel.yaml sdist SHA mismatch" "${ruamel_yaml_sdist_sha}"
  "8b7ce697a2f212752a35c1ac414471dc16c424c9573be4926b56ff3f5d23b7af")
require_equal("ruamel.yaml sdist size mismatch" "${ruamel_yaml_sdist_bytes}"
  "128123")
require_equal("ruamel.yaml Python requirement mismatch"
  "${ruamel_yaml_python_requires}" ">=3")
require_equal("ruamel.yaml declared dependency mismatch"
  "${ruamel_yaml_requires_dist}"
  "ruamel.yaml.clib (>=0.2.6) ; platform_python_implementation==\"CPython\" and python_version<\"3.11\"")
require_equal("ruamel.yaml copyright-holder mismatch"
  "${ruamel_yaml_copyright_holder}" "Anthon van der Neut, Ruamel bvba")
require_equal("ruamel.yaml license SPDX mismatch" "${ruamel_yaml_license_spdx}"
  "MIT")
require_equal("ruamel.yaml wheel license SHA mismatch"
  "${ruamel_yaml_wheel_license_sha}"
  "ab837b032c5aae84503fc0c733a116a26fd272e90dc4402fa68d3c9e51aed3b0")
require_equal("ruamel.yaml wheel license path mismatch"
  "${ruamel_yaml_wheel_license_path}"
  "inputs/dependencies/ruamel-yaml-0.17.21/wheel-extract/ruamel.yaml-0.17.21.dist-info/LICENSE")
require_equal("ruamel.yaml sdist license SHA mismatch"
  "${ruamel_yaml_sdist_license_sha}"
  "ab837b032c5aae84503fc0c733a116a26fd272e90dc4402fa68d3c9e51aed3b0")
require_equal("ruamel.yaml sdist license path mismatch"
  "${ruamel_yaml_sdist_license_path}"
  "inputs/dependencies/ruamel-yaml-0.17.21/sdist-extract/ruamel.yaml-0.17.21/LICENSE")
require_equal("ruamel.yaml license status mismatch" "${ruamel_yaml_license_status}"
  "verified_matching_wheel_and_sdist")
require_equal("ruamel.yaml consumption status mismatch"
  "${ruamel_yaml_consumption_status}"
  "required_cantera_clib_sourcegen_safe_yaml_for_p0_2")
require_equal("ruamel.yaml injection mismatch" "${ruamel_yaml_injection}"
  "exact_wheel_extract_via_pythonpath")
require_equal("ruamel.yaml rootfs install status mismatch"
  "${ruamel_yaml_rootfs_status}" "forbidden")
require_equal("ruamel.yaml artifact audit status mismatch"
  "${ruamel_yaml_artifact_status}" "pending_p0_2")
require_equal("ruamel.yaml bundle status mismatch" "${ruamel_yaml_bundle_status}"
  "forbidden_from_artifact_install")
require_equal("ruamel.yaml runtime status mismatch"
  "${ruamel_yaml_runtime_status}" "forbidden_runtime_dependency")

verify_external_file_bytes("${ruamel_yaml_metadata_path}"
  "${ruamel_yaml_metadata_sha}" "${ruamel_yaml_metadata_bytes}"
  "ruamel.yaml PyPI metadata")
verify_pypi_metadata("${ruamel_yaml_metadata_path}" "ruamel.yaml" "0.17.21"
  "MIT license" "ruamel.yaml")
verify_external_file_bytes("${ruamel_yaml_wheel_path}"
  "${ruamel_yaml_wheel_sha}" "${ruamel_yaml_wheel_bytes}"
  "ruamel.yaml wheel")
verify_external_file("${ruamel_yaml_wheel_metadata_path}"
  "${ruamel_yaml_wheel_metadata_sha}" "ruamel.yaml wheel METADATA")
verify_metadata_identity("${ruamel_yaml_wheel_metadata_path}" "ruamel.yaml"
  "0.17.21" ">=3" "ruamel.yaml")
verify_ruamel_wheel_tree("${ruamel_yaml_wheel_tree_path}"
  "${ruamel_yaml_wheel_tree_manifest_path}"
  "${ruamel_yaml_wheel_tree_manifest_sha}"
  "${ruamel_yaml_wheel_tree_file_count}" "ruamel.yaml")
verify_external_file_bytes("${ruamel_yaml_sdist_path}"
  "${ruamel_yaml_sdist_sha}" "${ruamel_yaml_sdist_bytes}"
  "ruamel.yaml sdist")
verify_external_file("${ruamel_yaml_wheel_license_path}"
  "${ruamel_yaml_wheel_license_sha}" "ruamel.yaml wheel license")
verify_external_file("${ruamel_yaml_sdist_license_path}"
  "${ruamel_yaml_sdist_license_sha}" "ruamel.yaml sdist license")

string(JSON ruamel_yaml_patch_count ERROR_VARIABLE ruamel_yaml_patch_error
  LENGTH "${manifest_json}" builder_python_dependencies 1 local_patches)
if(NOT ruamel_yaml_patch_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR
    "ruamel.yaml local-patch array is invalid: ${ruamel_yaml_patch_error}")
endif()
require_equal("ruamel.yaml local-patch count mismatch"
  "${ruamel_yaml_patch_count}" "0")

manifest_get(ruamel_clib_name builder_python_dependencies 2 name)
manifest_get(ruamel_clib_distribution builder_python_dependencies 2 distribution)
manifest_get(ruamel_clib_role builder_python_dependencies 2 role)
manifest_get(ruamel_clib_version builder_python_dependencies 2 version)
manifest_get(ruamel_clib_repository builder_python_dependencies 2 repository)
manifest_get(ruamel_clib_project_url builder_python_dependencies 2 project_url)
manifest_get(ruamel_clib_metadata_url builder_python_dependencies 2 metadata_url)
manifest_get(ruamel_clib_metadata_path builder_python_dependencies 2 metadata_path)
manifest_get(ruamel_clib_metadata_sha builder_python_dependencies 2 metadata_sha256)
manifest_get(ruamel_clib_metadata_bytes builder_python_dependencies 2 metadata_bytes)
manifest_get(ruamel_clib_wheel_url builder_python_dependencies 2 wheel_url)
manifest_get(ruamel_clib_wheel_path builder_python_dependencies 2 wheel_path)
manifest_get(ruamel_clib_wheel_sha builder_python_dependencies 2 wheel_sha256)
manifest_get(ruamel_clib_wheel_bytes builder_python_dependencies 2 wheel_bytes)
manifest_get(ruamel_clib_wheel_tag builder_python_dependencies 2 wheel_tag)
manifest_get(ruamel_clib_wheel_metadata_path builder_python_dependencies 2
  wheel_metadata_path)
manifest_get(ruamel_clib_wheel_metadata_sha builder_python_dependencies 2
  wheel_metadata_sha256)
manifest_get(ruamel_clib_wheel_tree_path builder_python_dependencies 2
  wheel_tree_path)
manifest_get(ruamel_clib_wheel_tree_manifest_path builder_python_dependencies 2
  wheel_tree_manifest_path)
manifest_get(ruamel_clib_wheel_tree_manifest_sha builder_python_dependencies 2
  wheel_tree_manifest_sha256)
manifest_get(ruamel_clib_wheel_tree_file_count builder_python_dependencies 2
  wheel_tree_file_count)
manifest_get(ruamel_clib_sdist_url builder_python_dependencies 2 sdist_url)
manifest_get(ruamel_clib_sdist_path builder_python_dependencies 2 sdist_path)
manifest_get(ruamel_clib_sdist_sha builder_python_dependencies 2 sdist_sha256)
manifest_get(ruamel_clib_sdist_bytes builder_python_dependencies 2 sdist_bytes)
manifest_get(ruamel_clib_python_requires builder_python_dependencies 2
  python_requires)
manifest_get(ruamel_clib_requires_dist builder_python_dependencies 2
  requires_dist)
manifest_get(ruamel_clib_copyright_holder builder_python_dependencies 2
  copyright_holder)
manifest_get(ruamel_clib_license_spdx builder_python_dependencies 2 license spdx)
manifest_get(ruamel_clib_wheel_license_path builder_python_dependencies 2
  license wheel_path)
manifest_get(ruamel_clib_wheel_license_sha builder_python_dependencies 2
  license wheel_sha256)
manifest_get(ruamel_clib_sdist_license_path builder_python_dependencies 2
  license sdist_path)
manifest_get(ruamel_clib_sdist_license_sha builder_python_dependencies 2
  license sdist_sha256)
manifest_get(ruamel_clib_license_status builder_python_dependencies 2
  license status)
manifest_get(ruamel_clib_abi_python builder_python_dependencies 2 abi python)
manifest_get(ruamel_clib_abi_arch builder_python_dependencies 2 abi architecture)
manifest_get(ruamel_clib_abi_manylinux builder_python_dependencies 2 abi manylinux_floor)
manifest_get(ruamel_clib_abi_needed builder_python_dependencies 2 abi needed)
manifest_get(ruamel_clib_abi_glibc builder_python_dependencies 2 abi maximum_glibc_symbol)
manifest_get(ruamel_clib_consumption_status builder_python_dependencies 2
  consumption_status)
manifest_get(ruamel_clib_injection builder_python_dependencies 2 injection)
manifest_get(ruamel_clib_rootfs_status builder_python_dependencies 2
  rootfs_install_status)
manifest_get(ruamel_clib_artifact_status builder_python_dependencies 2
  artifact_audit_status)
manifest_get(ruamel_clib_bundle_status builder_python_dependencies 2 bundle_status)
manifest_get(ruamel_clib_runtime_status builder_python_dependencies 2
  runtime_status)

require_equal("ruamel.yaml.clib name mismatch" "${ruamel_clib_name}"
  "ruamel.yaml.clib")
require_equal("ruamel.yaml.clib distribution mismatch"
  "${ruamel_clib_distribution}" "ruamel.yaml.clib")
require_equal("ruamel.yaml.clib role mismatch" "${ruamel_clib_role}"
  "builder_pythonpath_only")
require_equal("ruamel.yaml.clib version mismatch" "${ruamel_clib_version}"
  "0.2.6")
require_equal("ruamel.yaml.clib repository mismatch" "${ruamel_clib_repository}"
  "Python Package Index")
require_equal("ruamel.yaml.clib project URL mismatch"
  "${ruamel_clib_project_url}"
  "https://pypi.org/project/ruamel.yaml.clib/0.2.6/")
require_equal("ruamel.yaml.clib metadata URL mismatch"
  "${ruamel_clib_metadata_url}"
  "https://pypi.org/pypi/ruamel.yaml.clib/0.2.6/json")
require_equal("ruamel.yaml.clib metadata path mismatch"
  "${ruamel_clib_metadata_path}"
  "inputs/dependencies/ruamel-yaml-clib-0.2.6/pypi-metadata.json")
require_equal("ruamel.yaml.clib metadata SHA mismatch"
  "${ruamel_clib_metadata_sha}"
  "677b9e08e0143864232cccc1122275d0bf33f79af9aaf464ba19fc6789953c88")
require_equal("ruamel.yaml.clib metadata size mismatch"
  "${ruamel_clib_metadata_bytes}" "28289")
require_equal("ruamel.yaml.clib wheel URL mismatch" "${ruamel_clib_wheel_url}"
  "https://files.pythonhosted.org/packages/e1/c6/a8ed2b252c9a1018ea1758bbfa6bcd1b4965009e4f9040e1d0456417d7ef/ruamel.yaml.clib-0.2.6-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.manylinux_2_24_x86_64.whl")
require_equal("ruamel.yaml.clib wheel path mismatch"
  "${ruamel_clib_wheel_path}"
  "inputs/dependencies/ruamel-yaml-clib-0.2.6/ruamel.yaml.clib-0.2.6-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.manylinux_2_24_x86_64.whl")
require_equal("ruamel.yaml.clib wheel SHA mismatch" "${ruamel_clib_wheel_sha}"
  "221eca6f35076c6ae472a531afa1c223b9c29377e62936f61bc8e6e8bdc5f9e7")
require_equal("ruamel.yaml.clib wheel size mismatch"
  "${ruamel_clib_wheel_bytes}" "519289")
require_equal("ruamel.yaml.clib wheel tag mismatch" "${ruamel_clib_wheel_tag}"
  "cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.manylinux_2_24_x86_64")
require_equal("ruamel.yaml.clib wheel METADATA SHA mismatch"
  "${ruamel_clib_wheel_metadata_sha}"
  "ee7c146559d808531a6b8729460ba7bc62140e8dc7fe9279ef6001ea32ede75c")
require_equal("ruamel.yaml.clib wheel METADATA path mismatch"
  "${ruamel_clib_wheel_metadata_path}"
  "inputs/dependencies/ruamel-yaml-clib-0.2.6/wheel-extract/ruamel.yaml.clib-0.2.6.dist-info/METADATA")
require_equal("ruamel.yaml.clib wheel tree path mismatch"
  "${ruamel_clib_wheel_tree_path}"
  "inputs/dependencies/ruamel-yaml-clib-0.2.6/wheel-extract")
require_equal("ruamel.yaml.clib wheel tree manifest path mismatch"
  "${ruamel_clib_wheel_tree_manifest_path}"
  "manifests/ruamel-yaml-clib-0.2.6-wheel-extract.sha256")
require_equal("ruamel.yaml.clib wheel tree manifest SHA mismatch"
  "${ruamel_clib_wheel_tree_manifest_sha}"
  "e1d31c44b7d6e12418793fcf1340f39b1a020efd7c095a381b9c77202415ff8a")
require_equal("ruamel.yaml.clib wheel tree file count mismatch"
  "${ruamel_clib_wheel_tree_file_count}" "7")
require_equal("ruamel.yaml.clib sdist URL mismatch" "${ruamel_clib_sdist_url}"
  "https://files.pythonhosted.org/packages/8b/25/08e5ad2431a028d0723ca5540b3af6a32f58f25e83c6dda4d0fcef7288a3/ruamel.yaml.clib-0.2.6.tar.gz")
require_equal("ruamel.yaml.clib sdist path mismatch"
  "${ruamel_clib_sdist_path}"
  "inputs/dependencies/ruamel-yaml-clib-0.2.6/ruamel.yaml.clib-0.2.6.tar.gz")
require_equal("ruamel.yaml.clib sdist SHA mismatch" "${ruamel_clib_sdist_sha}"
  "4ff604ce439abb20794f05613c374759ce10e3595d1867764dd1ae675b85acbd")
require_equal("ruamel.yaml.clib sdist size mismatch"
  "${ruamel_clib_sdist_bytes}" "180695")
require_equal("ruamel.yaml.clib Python requirement mismatch"
  "${ruamel_clib_python_requires}" ">=3.5")
require_equal("ruamel.yaml.clib declared dependency mismatch"
  "${ruamel_clib_requires_dist}" "none")
require_equal("ruamel.yaml.clib copyright-holder mismatch"
  "${ruamel_clib_copyright_holder}" "Anthon van der Neut, Ruamel bvba")
require_equal("ruamel.yaml.clib license SPDX mismatch"
  "${ruamel_clib_license_spdx}" "MIT")
require_equal("ruamel.yaml.clib wheel license SHA mismatch"
  "${ruamel_clib_wheel_license_sha}"
  "16174d2cf8c2ee4b900bf8573106bc61c1ce0092da4fa781cc0ee81047a46539")
require_equal("ruamel.yaml.clib wheel license path mismatch"
  "${ruamel_clib_wheel_license_path}"
  "inputs/dependencies/ruamel-yaml-clib-0.2.6/wheel-extract/ruamel.yaml.clib-0.2.6.dist-info/LICENSE")
require_equal("ruamel.yaml.clib sdist license SHA mismatch"
  "${ruamel_clib_sdist_license_sha}"
  "16174d2cf8c2ee4b900bf8573106bc61c1ce0092da4fa781cc0ee81047a46539")
require_equal("ruamel.yaml.clib sdist license path mismatch"
  "${ruamel_clib_sdist_license_path}"
  "inputs/dependencies/ruamel-yaml-clib-0.2.6/sdist-extract/ruamel.yaml.clib-0.2.6/LICENSE")
require_equal("ruamel.yaml.clib license status mismatch"
  "${ruamel_clib_license_status}" "verified_matching_wheel_and_sdist")
require_equal("ruamel.yaml.clib ABI Python mismatch" "${ruamel_clib_abi_python}"
  "CPython 3.10")
require_equal("ruamel.yaml.clib ABI architecture mismatch"
  "${ruamel_clib_abi_arch}" "x86_64")
require_equal("ruamel.yaml.clib ABI manylinux floor mismatch"
  "${ruamel_clib_abi_manylinux}" "glibc 2.17")
require_equal("ruamel.yaml.clib ABI needed mismatch" "${ruamel_clib_abi_needed}"
  "libpthread.so.0,libc.so.6")
require_equal("ruamel.yaml.clib ABI glibc symbol mismatch"
  "${ruamel_clib_abi_glibc}" "GLIBC_2.14")
require_equal("ruamel.yaml.clib consumption status mismatch"
  "${ruamel_clib_consumption_status}"
  "required_ruamel_yaml_declared_closure_for_python310_p0_2")
require_equal("ruamel.yaml.clib injection mismatch" "${ruamel_clib_injection}"
  "exact_wheel_extract_via_pythonpath")
require_equal("ruamel.yaml.clib rootfs install status mismatch"
  "${ruamel_clib_rootfs_status}" "forbidden")
require_equal("ruamel.yaml.clib artifact audit status mismatch"
  "${ruamel_clib_artifact_status}" "pending_p0_2")
require_equal("ruamel.yaml.clib bundle status mismatch"
  "${ruamel_clib_bundle_status}" "forbidden_from_artifact_install")
require_equal("ruamel.yaml.clib runtime status mismatch"
  "${ruamel_clib_runtime_status}" "forbidden_runtime_dependency")

verify_external_file_bytes("${ruamel_clib_metadata_path}"
  "${ruamel_clib_metadata_sha}" "${ruamel_clib_metadata_bytes}"
  "ruamel.yaml.clib PyPI metadata")
verify_pypi_metadata("${ruamel_clib_metadata_path}" "ruamel.yaml.clib" "0.2.6"
  "MIT" "ruamel.yaml.clib")
verify_external_file_bytes("${ruamel_clib_wheel_path}"
  "${ruamel_clib_wheel_sha}" "${ruamel_clib_wheel_bytes}"
  "ruamel.yaml.clib wheel")
verify_external_file("${ruamel_clib_wheel_metadata_path}"
  "${ruamel_clib_wheel_metadata_sha}" "ruamel.yaml.clib wheel METADATA")
verify_metadata_identity("${ruamel_clib_wheel_metadata_path}"
  "ruamel.yaml.clib" "0.2.6" ">=3.5" "ruamel.yaml.clib")
verify_ruamel_wheel_tree("${ruamel_clib_wheel_tree_path}"
  "${ruamel_clib_wheel_tree_manifest_path}"
  "${ruamel_clib_wheel_tree_manifest_sha}"
  "${ruamel_clib_wheel_tree_file_count}" "ruamel.yaml.clib")
verify_external_file_bytes("${ruamel_clib_sdist_path}"
  "${ruamel_clib_sdist_sha}" "${ruamel_clib_sdist_bytes}"
  "ruamel.yaml.clib sdist")
verify_external_file("${ruamel_clib_wheel_license_path}"
  "${ruamel_clib_wheel_license_sha}" "ruamel.yaml.clib wheel license")
verify_external_file("${ruamel_clib_sdist_license_path}"
  "${ruamel_clib_sdist_license_sha}" "ruamel.yaml.clib sdist license")

string(JSON ruamel_clib_patch_count ERROR_VARIABLE ruamel_clib_patch_error
  LENGTH "${manifest_json}" builder_python_dependencies 2 local_patches)
if(NOT ruamel_clib_patch_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR
    "ruamel.yaml.clib local-patch array is invalid: ${ruamel_clib_patch_error}")
endif()
require_equal("ruamel.yaml.clib local-patch count mismatch"
  "${ruamel_clib_patch_count}" "0")

string(JSON mechanism_count ERROR_VARIABLE mechanism_error
  LENGTH "${manifest_json}" mechanisms)
if(NOT mechanism_error STREQUAL "NOTFOUND")
  message(FATAL_ERROR "Mechanism separation record is missing")
endif()
require_equal("Mechanism count mismatch" "${mechanism_count}" "1")
manifest_get(mechanism_name mechanisms 0 name)
manifest_get(mechanism_phase mechanisms 0 phase)
manifest_get(mechanism_path mechanisms 0 path)
manifest_get(mechanism_sha mechanisms 0 sha256)
manifest_get(mechanism_status mechanisms 0 license_status)
manifest_get(mechanism_spdx mechanisms 0 spdx)
manifest_get(mechanism_consumption_status mechanisms 0 consumption_status)
manifest_get(mechanism_license_path mechanisms 0 license_path)
manifest_get(mechanism_license_sha mechanisms 0 license_sha256)
manifest_get(terms_path mechanisms 0 terms_evidence_path)
manifest_get(terms_sha mechanisms 0 terms_evidence_sha256)
require_equal("Mechanism name mismatch" "${mechanism_name}" "h2o2.yaml")
require_equal("Mechanism phase mismatch" "${mechanism_phase}" "ohmech")
require_equal("Mechanism license status mismatch" "${mechanism_status}"
  "candidate_user_supplied")
require_equal("Mechanism SPDX mismatch" "${mechanism_spdx}" "NOASSERTION")
require_equal("Mechanism consumption status mismatch"
  "${mechanism_consumption_status}"
  "blocked_until_asset_specific_redistribution_permission")
require_equal("Mechanism license path mismatch" "${mechanism_license_path}" "")
require_equal("Mechanism license SHA mismatch" "${mechanism_license_sha}" "")
verify_external_file("${mechanism_path}" "${mechanism_sha}" "Mechanism candidate")
verify_external_file("${terms_path}" "${terms_sha}" "Mechanism terms evidence")

message(STATUS "Stage 4 P0 input manifest validated")
