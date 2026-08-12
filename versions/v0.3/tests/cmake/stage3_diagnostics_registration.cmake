# SPDX-License-Identifier: Apache-2.0

include_guard(GLOBAL)

# Exclusive registration owners: S3-O1 and S3-O2.

hundun_add_test(test_wale_diagnostics_header_contract
  ${PROJECT_SOURCE_DIR}/tests/unit/test_wale_diagnostics_header_contract.cpp)
target_include_directories(test_wale_diagnostics_header_contract PRIVATE
  "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_wale_diagnostics_header_contract PRIVATE
  hundun_session_diagnostics)
set_tests_properties(test_wale_diagnostics_header_contract PROPERTIES
  LABELS "unit;header;stage3;s3-o1;fast"
  TIMEOUT 30)

hundun_add_test(test_wale_diagnostics
  ${PROJECT_SOURCE_DIR}/tests/unit/test_wale_diagnostics.cpp)
target_include_directories(test_wale_diagnostics PRIVATE
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_wale_diagnostics PRIVATE
  hundun_session_diagnostics)
set_tests_properties(test_wale_diagnostics PROPERTIES
  LABELS "unit;diagnostics;stage3;s3-o1;fast"
  TIMEOUT 30)

hundun_add_test(test_immersed_static_diagnostics_header_contract
  ${PROJECT_SOURCE_DIR}/tests/unit/test_immersed_static_diagnostics_header_contract.cpp)
target_include_directories(test_immersed_static_diagnostics_header_contract
  PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_immersed_static_diagnostics_header_contract PRIVATE
  hundun_session_diagnostics)
set_tests_properties(test_immersed_static_diagnostics_header_contract PROPERTIES
  LABELS "unit;header;stage3;s3-o2;fast"
  TIMEOUT 30)

hundun_add_test(test_stage3_provider_inventory
  ${PROJECT_SOURCE_DIR}/tests/unit/test_stage3_provider_inventory.cpp)
target_include_directories(test_stage3_provider_inventory PRIVATE
  "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_stage3_provider_inventory PRIVATE
  hundun_session_diagnostics hundun_checkpoint_diagnostics)
set_tests_properties(test_stage3_provider_inventory PROPERTIES
  LABELS "unit;diagnostics;stage3;s3-o2;fast"
  TIMEOUT 30)

add_executable(test_stage3_diagnostics_mpi
  ${PROJECT_SOURCE_DIR}/tests/mpi/test_stage3_diagnostics_mpi.cpp)
target_include_directories(test_stage3_diagnostics_mpi PRIVATE
  "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_stage3_diagnostics_mpi PRIVATE
  hundun_session_diagnostics stage3_stl_fixture
  hundun_options hundun_warnings MPI::MPI_CXX)
foreach(ranks IN ITEMS 1 2)
  add_test(NAME "test_stage3_diagnostics_mpi_${ranks}_rank"
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_stage3_diagnostics_mpi>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties("test_stage3_diagnostics_mpi_${ranks}_rank" PROPERTIES
    LABELS "mpi;diagnostics;stage3;s3-o2;fast"
    TIMEOUT 180
    PROCESSORS "${ranks}"
    RESOURCE_LOCK hundun_stage3_mpi_m)
endforeach()
