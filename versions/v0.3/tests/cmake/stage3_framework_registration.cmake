# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

include_guard(GLOBAL)

# Exclusive registration owners: S3-A1 and S3-E1.

add_executable(test_stage3_dispatch_inventory
  ${PROJECT_SOURCE_DIR}/tests/unit/test_stage3_dispatch_inventory.cpp
  ${PROJECT_SOURCE_DIR}/src/app_immersed_flow_driver.cpp)
target_include_directories(test_stage3_dispatch_inventory PRIVATE
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_SOURCE_DIR}/include"
  "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(test_stage3_dispatch_inventory PRIVATE
  hundun_application
  hundun_cli
  hundun_config
  hundun_checkpoint_diagnostics
  hundun_diagnostics
  hundun_material_diagnostics
  hundun_session_diagnostics
  hundun_flow
  hundun_io
  hundun_solver
  hundun_options
  hundun_warnings
  MPI::MPI_CXX)
add_test(NAME test_stage3_dispatch_inventory
  COMMAND "$<TARGET_FILE:test_stage3_dispatch_inventory>")
set_tests_properties(test_stage3_dispatch_inventory PROPERTIES
  LABELS "unit;dispatch;stage3;s3-a1;fast;mutation"
  TIMEOUT 30)

hundun_add_test(test_stage3_performance_header_contract
  ${PROJECT_SOURCE_DIR}/tests/unit/test_stage3_performance_header_contract.cpp)
target_include_directories(test_stage3_performance_header_contract PRIVATE
  "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_stage3_performance_header_contract PRIVATE
  hundun_diagnostics)
set_tests_properties(test_stage3_performance_header_contract PROPERTIES
  LABELS "unit;header;performance;stage3;s3-e1;fast"
  TIMEOUT 30)

hundun_add_test(test_stage3_exact_counters
  ${PROJECT_SOURCE_DIR}/tests/unit/test_stage3_exact_counters.cpp)
target_include_directories(test_stage3_exact_counters PRIVATE
  "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(test_stage3_exact_counters PRIVATE hundun_diagnostics)
set_tests_properties(test_stage3_exact_counters PROPERTIES
  LABELS "unit;performance;stage3;s3-e1;fast;mutation"
  TIMEOUT 30)

add_executable(test_stage3_performance
  ${PROJECT_SOURCE_DIR}/tests/mpi/test_stage3_performance.cpp
  ${PROJECT_SOURCE_DIR}/tests/support/stage3_case_generator.cpp
  ${PROJECT_SOURCE_DIR}/tests/support/stage3_performance_evidence.cpp)
target_include_directories(test_stage3_performance PRIVATE
  "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include"
  "${PROJECT_SOURCE_DIR}/src" "${PROJECT_BINARY_DIR}/src/generated")
target_link_libraries(test_stage3_performance PRIVATE
  hundun_diagnostics stage3_stl_fixture hundun_flow
  hundun_options hundun_warnings MPI::MPI_CXX)
target_compile_definitions(test_stage3_performance PRIVATE
  HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)

foreach(ranks IN ITEMS 1 2)
  add_test(NAME "test_stage3_performance_${ranks}_rank_fast"
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_stage3_performance>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties("test_stage3_performance_${ranks}_rank_fast" PROPERTIES
    LABELS "mpi;performance;stage3;s3-e1;fast"
    TIMEOUT 300
    PROCESSORS "${ranks}"
    RESOURCE_LOCK hundun_stage3_mpi_m)
endforeach()

foreach(ranks IN ITEMS 1 2 4)
  add_test(NAME "test_stage3_performance_${ranks}_rank_formal"
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_stage3_performance>"
            formal 24 2 3 1 ${MPIEXEC_POSTFLAGS})
  set_tests_properties("test_stage3_performance_${ranks}_rank_formal"
    PROPERTIES
    LABELS "mpi;performance;stage3;s3-e1;formal"
    TIMEOUT 1800
    PROCESSORS "${ranks}"
    RESOURCE_LOCK hundun_stage3_performance_m)
endforeach()
