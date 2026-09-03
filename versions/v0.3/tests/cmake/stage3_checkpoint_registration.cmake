# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

include_guard(GLOBAL)

# Exclusive registration owners: S3-R1 and S3-R2.

add_executable(write_stage3_restart_stl
  ${PROJECT_SOURCE_DIR}/tests/support/write_task19a_stl.cpp)
target_include_directories(write_stage3_restart_stl PRIVATE
  "${PROJECT_SOURCE_DIR}")
target_link_libraries(write_stage3_restart_stl PRIVATE
  stage3_stl_fixture hundun_options hundun_warnings)

add_executable(test_checkpoint_v3_wale
  ${PROJECT_SOURCE_DIR}/tests/mpi/test_checkpoint_v3_wale.cpp)
target_include_directories(test_checkpoint_v3_wale PRIVATE
  "${PROJECT_SOURCE_DIR}")
target_compile_definitions(test_checkpoint_v3_wale PRIVATE
  HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)
target_link_libraries(test_checkpoint_v3_wale PRIVATE
  hundun_flow hundun_les stage3_mms stage3_stl_fixture
  hundun_options hundun_warnings MPI::MPI_CXX)

foreach(ranks IN ITEMS 1 2)
  add_test(NAME "test_checkpoint_v3_wale_${ranks}_rank"
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_checkpoint_v3_wale>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties("test_checkpoint_v3_wale_${ranks}_rank" PROPERTIES
    LABELS "mpi;numerical;stage3;s3-r1;fast"
    TIMEOUT 120
    PROCESSORS "${ranks}"
    RESOURCE_LOCK hundun_stage3_mpi_m)
endforeach()

foreach(ranks IN ITEMS 1 2)
  add_test(NAME "test_stage3_restart_fast_${ranks}_rank"
    COMMAND bash
            "${PROJECT_SOURCE_DIR}/tests/acceptance/stage3_restart_fast.sh"
            "$<TARGET_FILE:hundun>" "$<TARGET_FILE:write_stage3_restart_stl>"
            "${MPIEXEC_EXECUTABLE}" "${ranks}")
  set_tests_properties("test_stage3_restart_fast_${ranks}_rank" PROPERTIES
    LABELS "acceptance;mpi;numerical;restart;stage3;s3-r2;fast"
    TIMEOUT 180
    PROCESSORS "${ranks}"
    RESOURCE_LOCK hundun_stage3_mpi_m)
endforeach()

add_executable(test_checkpoint_v3_density_profiles
  ${PROJECT_SOURCE_DIR}/tests/mpi/test_checkpoint_v3_density_profiles.cpp)
target_include_directories(test_checkpoint_v3_density_profiles PRIVATE
  "${PROJECT_SOURCE_DIR}")
target_compile_definitions(test_checkpoint_v3_density_profiles PRIVATE
  HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)
target_link_libraries(test_checkpoint_v3_density_profiles PRIVATE
  hundun_flow hundun_les stage3_stl_fixture
  hundun_options hundun_warnings MPI::MPI_CXX)
add_dependencies(test_checkpoint_v3_density_profiles write_stage3_restart_stl)

foreach(ranks IN ITEMS 1 2)
  add_test(NAME "test_checkpoint_v3_density_profiles_${ranks}_rank"
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:test_checkpoint_v3_density_profiles>" fast
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    "test_checkpoint_v3_density_profiles_${ranks}_rank" PROPERTIES
    LABELS "mpi;numerical;restart;stage3;s3-r2;fast"
    TIMEOUT 180
    PROCESSORS "${ranks}"
    RESOURCE_LOCK hundun_stage3_mpi_m)
endforeach()

foreach(ranks IN ITEMS 1 2 4)
  add_test(NAME "checkpoint-continuation-n12-r${ranks}"
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:test_checkpoint_v3_density_profiles>" formal 12
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties("checkpoint-continuation-n12-r${ranks}" PROPERTIES
    LABELS "formal;scientific;restart;stage3"
    TIMEOUT 1800
    PROCESSORS "${ranks}"
    RESOURCE_LOCK hundun_stage3_mpi_m)
endforeach()
