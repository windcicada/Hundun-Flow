# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

include_guard(GLOBAL)

# Exclusive registration owners: S3-G1, S3-DOC contract, S3-V0, and S3-V1.

hundun_add_test(test_stage3_evidence_manifest
  ${PROJECT_SOURCE_DIR}/tests/unit/test_stage3_evidence_manifest.cpp)
target_include_directories(test_stage3_evidence_manifest PRIVATE
  "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_stage3_evidence_manifest PRIVATE hundun_diagnostics)
set_tests_properties(test_stage3_evidence_manifest PROPERTIES
  LABELS "unit;governance;stage3;s3-g1;mutation" TIMEOUT 30)

add_test(NAME test_stage3_acceptance_contract COMMAND "${CMAKE_COMMAND}"
  "-DHUNDUN_STAGE3_ACCEPTANCE_SCRIPT=${PROJECT_SOURCE_DIR}/tests/acceptance/stage3_acceptance.sh"
  "-DHUNDUN_STAGE3_ACCEPTANCE_INVENTORY=${PROJECT_SOURCE_DIR}/tests/acceptance/stage3_acceptance_inventory.tsv"
  -P "${PROJECT_SOURCE_DIR}/tests/cmake/stage3_acceptance_contract.cmake")
add_test(NAME test_stage3_capability_ledger COMMAND "${CMAKE_COMMAND}"
  "-DHUNDUN_STAGE3_CAPABILITY_LEDGER=${PROJECT_SOURCE_DIR}/docs/numerics/stage3-capability-ledger.md"
  -P "${PROJECT_SOURCE_DIR}/tests/cmake/stage3_capability_ledger_contract.cmake")
add_test(NAME test_stage3_product_projection_contract COMMAND "${CMAKE_COMMAND}"
  "-DHUNDUN_STAGE3_PRODUCT_PROJECTION_MANIFEST=${PROJECT_SOURCE_DIR}/.superpowers/product-projection-manifest-2026-08-09.tsv"
  "-DHUNDUN_STAGE3_PRODUCT_PROJECTOR=${PROJECT_SOURCE_DIR}/tests/cmake/stage3_product_projection.cmake"
  -P "${PROJECT_SOURCE_DIR}/tests/cmake/stage3_product_projection_contract.cmake")
foreach(test IN ITEMS test_stage3_acceptance_contract test_stage3_capability_ledger test_stage3_product_projection_contract)
  set_tests_properties(${test} PROPERTIES LABELS "unit;contract;governance;stage3;s3-g1" TIMEOUT 30)
endforeach()

add_test(NAME test_stage3_documentation_contract COMMAND "${CMAKE_COMMAND}"
  "-DHUNDUN_SOURCE_ROOT=${PROJECT_SOURCE_DIR}"
  -P "${PROJECT_SOURCE_DIR}/tests/cmake/stage3_documentation_contract.cmake")
set_tests_properties(test_stage3_documentation_contract PROPERTIES
  LABELS "unit;contract;documentation;stage3;s3-doc" TIMEOUT 30)
