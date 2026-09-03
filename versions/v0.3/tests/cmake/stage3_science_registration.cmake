# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

include_guard(GLOBAL)

# Exclusive registration owners: S3-C1, S3-D1, S3-C2, S3-D2, S3-C3, S3-S1.

add_library(stage3_scientific_row STATIC
  ${PROJECT_SOURCE_DIR}/tests/support/stage3_scientific_row.cpp)
target_include_directories(stage3_scientific_row PUBLIC
  "${PROJECT_SOURCE_DIR}")
target_link_libraries(stage3_scientific_row PUBLIC
  hundun_runtime hundun_options hundun_warnings)

add_executable(test_stage3_scientific_row_contract
  ${PROJECT_SOURCE_DIR}/tests/unit/test_stage3_scientific_row.cpp)
target_include_directories(test_stage3_scientific_row_contract PRIVATE
  "${PROJECT_SOURCE_DIR}")
target_link_libraries(test_stage3_scientific_row_contract PRIVATE
  stage3_scientific_row hundun_options hundun_warnings)
add_test(NAME test_stage3_scientific_row_contract
  COMMAND "$<TARGET_FILE:test_stage3_scientific_row_contract>")
set_tests_properties(test_stage3_scientific_row_contract PROPERTIES
  LABELS "unit;contract;mutation;stage3;s3-s1;fast"
  TIMEOUT 30)

add_executable(test_wale_taylor_green
  ${PROJECT_SOURCE_DIR}/tests/numerical/test_wale_taylor_green.cpp)
target_include_directories(test_wale_taylor_green PRIVATE
  "${PROJECT_SOURCE_DIR}")
target_compile_definitions(test_wale_taylor_green PRIVATE
  HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)
target_link_libraries(test_wale_taylor_green PRIVATE
  hundun_flow hundun_les stage3_scientific_row
  hundun_options hundun_warnings MPI::MPI_CXX)
add_test(NAME test_wale_taylor_green_12_smoke_1_rank
  COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
          ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_wale_taylor_green>"
          smoke 12 ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_wale_taylor_green_12_smoke_1_rank PROPERTIES
  LABELS "mpi;numerical;smoke;stage3;s3-s1;fast"
  TIMEOUT 300
  PROCESSORS 1
  RESOURCE_LOCK hundun_stage3_mpi_m)

function(hundun_add_stage3_scientific_mpi_test
    name target ranks timeout resource_lock)
  add_test(NAME "${name}"
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:${target}>" ${ARGN}
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties("${name}" PROPERTIES
    LABELS "formal;scientific;mpi;numerical;stage3;s3-s1"
    TIMEOUT "${timeout}"
    PROCESSORS "${ranks}"
    RESOURCE_LOCK "${resource_lock}")
endfunction()

hundun_add_stage3_scientific_mpi_test(
  test_wale_taylor_green_convergence_1_rank_formal
  test_wale_taylor_green 1 43200 hundun_stage3_mpi_h
  formal convergence)
foreach(ranks IN ITEMS 1 2 4)
  hundun_add_stage3_scientific_mpi_test(
    "test_wale_taylor_green_n24_${ranks}_rank_formal"
    test_wale_taylor_green "${ranks}" 7200 hundun_stage3_mpi_m
    formal 24)
endforeach()

# This target is declared later in tests/CMakeLists.txt.  Resolve its S1-only
# row-support dependency after all tests in this directory have been declared.
cmake_language(DEFER CALL target_link_libraries
  test_wale_body_fitted PRIVATE stage3_scientific_row)
hundun_add_stage3_scientific_mpi_test(
  test_wale_channel_n48_1_rank_formal
  test_wale_body_fitted 1 43200 hundun_stage3_mpi_h
  formal channel)

add_executable(test_immersed_wale_constant
  ${PROJECT_SOURCE_DIR}/tests/mpi/test_immersed_wale_constant.cpp)
target_include_directories(test_immersed_wale_constant PRIVATE
  "${PROJECT_SOURCE_DIR}")
target_compile_definitions(test_immersed_wale_constant PRIVATE
  HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)
target_link_libraries(test_immersed_wale_constant PRIVATE
  hundun_flow hundun_les stage3_mms stage3_stl_fixture
  stage3_scientific_row
  hundun_options hundun_warnings MPI::MPI_CXX)

hundun_add_stage3_scientific_mpi_test(
  test_constant_ibm_wale_n48_1_rank_formal
  test_immersed_wale_constant 1 43200 hundun_stage3_mpi_h
  formal 48)
foreach(ranks IN ITEMS 1 2 4)
  hundun_add_stage3_scientific_mpi_test(
    "test_constant_ibm_wale_n24_${ranks}_rank_formal"
    test_immersed_wale_constant "${ranks}" 7200 hundun_stage3_mpi_m
    formal 24)
endforeach()

add_executable(test_material_wale_composition
  ${PROJECT_SOURCE_DIR}/tests/mpi/test_material_wale_composition.cpp)
target_include_directories(test_material_wale_composition PRIVATE
  "${PROJECT_SOURCE_DIR}")
target_compile_definitions(test_material_wale_composition PRIVATE
  HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)
target_link_libraries(test_material_wale_composition PRIVATE
  hundun_flow hundun_les stage3_mms stage3_stl_fixture
  stage3_scientific_row
  hundun_options hundun_warnings MPI::MPI_CXX)

foreach(target IN ITEMS
    test_ideal_gas_wale_composition
    test_immersed_combined_retry)
  add_executable(${target}
    ${PROJECT_SOURCE_DIR}/tests/mpi/${target}.cpp)
  target_include_directories(${target} PRIVATE "${PROJECT_SOURCE_DIR}")
  target_compile_definitions(${target} PRIVATE
    HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)
  target_link_libraries(${target} PRIVATE
    hundun_flow hundun_les stage3_mms stage3_stl_fixture
    stage3_scientific_row
    hundun_options hundun_warnings MPI::MPI_CXX)
endforeach()

foreach(model IN ITEMS material ideal_gas)
  if(model STREQUAL "material")
    set(target test_material_wale_composition)
    set(row_prefix material_ibm_wale)
  else()
    set(target test_ideal_gas_wale_composition)
    set(row_prefix ideal_ibm_wale)
  endif()
  foreach(cells IN ITEMS 12 24)
    if(cells EQUAL 12)
      set(timeout 1800)
    else()
      set(timeout 7200)
    endif()
    foreach(ranks IN ITEMS 1 2 4)
      hundun_add_stage3_scientific_mpi_test(
        "test_${row_prefix}_n${cells}_${ranks}_rank_formal"
        "${target}" "${ranks}" "${timeout}" hundun_stage3_mpi_m
        formal "${cells}")
    endforeach()
  endforeach()
endforeach()

foreach(target IN ITEMS
    test_immersed_material_density
    test_immersed_material_transaction
    test_immersed_ideal_gas
    test_immersed_ideal_gas_transaction)
  add_executable(${target}
    ${PROJECT_SOURCE_DIR}/tests/mpi/${target}.cpp)
  target_include_directories(${target} PRIVATE "${PROJECT_SOURCE_DIR}")
  target_compile_definitions(${target} PRIVATE
    HUNDUN_FLOW_ENABLE_TEST_ACCESS=1)
  target_link_libraries(${target} PRIVATE
    hundun_flow hundun_les stage3_mms stage3_stl_fixture
    hundun_options hundun_warnings MPI::MPI_CXX)
endforeach()

foreach(kind IN ITEMS ideal_gas ideal_gas_transaction)
  foreach(ranks IN ITEMS 1 2)
    add_test(NAME "test_immersed_${kind}_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_immersed_${kind}>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_immersed_${kind}_${ranks}_rank"
      PROPERTIES
      LABELS "mpi;numerical;stage3;s3-d2;fast"
      TIMEOUT 180
      PROCESSORS "${ranks}"
      RESOURCE_LOCK hundun_stage3_mpi_m)
  endforeach()
endforeach()

add_test(NAME test_immersed_wale_constant_1_rank
  COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
          ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_wale_constant>"
          ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_immersed_wale_constant_1_rank PROPERTIES
  LABELS "mpi;numerical;red;stage3;s3-c1;fast"
  TIMEOUT 120
  PROCESSORS 1
  RESOURCE_LOCK hundun_stage3_mpi_m)

foreach(kind IN ITEMS density transaction)
  foreach(ranks IN ITEMS 1 2)
    add_test(NAME "test_immersed_material_${kind}_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS}
              "$<TARGET_FILE:test_immersed_material_${kind}>"
              ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_immersed_material_${kind}_${ranks}_rank"
      PROPERTIES
      LABELS "mpi;numerical;stage3;s3-d1;fast"
      TIMEOUT 120
      PROCESSORS "${ranks}"
      RESOURCE_LOCK hundun_stage3_mpi_m)
  endforeach()
endforeach()

foreach(ranks IN ITEMS 1 2)
  add_test(NAME "test_immersed_wale_constant_${ranks}_rank_fast"
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_immersed_wale_constant>"
            fast ${MPIEXEC_POSTFLAGS})
  set_tests_properties("test_immersed_wale_constant_${ranks}_rank_fast"
    PROPERTIES
    LABELS "mpi;numerical;stage3;s3-c1;fast"
    TIMEOUT 300
    PROCESSORS "${ranks}"
    RESOURCE_LOCK hundun_stage3_mpi_m)
endforeach()

foreach(ranks IN ITEMS 1 2)
  add_test(NAME "test_material_wale_composition_${ranks}_rank"
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
            ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_material_wale_composition>"
            fast ${MPIEXEC_POSTFLAGS})
  set_tests_properties("test_material_wale_composition_${ranks}_rank"
    PROPERTIES
    LABELS "mpi;numerical;stage3;s3-c2;fast"
    TIMEOUT 300
    PROCESSORS "${ranks}"
    RESOURCE_LOCK hundun_stage3_mpi_m)
endforeach()

foreach(kind IN ITEMS ideal_gas_wale_composition immersed_combined_retry)
  foreach(ranks IN ITEMS 1 2)
    add_test(NAME "test_${kind}_${ranks}_rank"
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ranks}"
              ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:test_${kind}>"
              fast ${MPIEXEC_POSTFLAGS})
    set_tests_properties("test_${kind}_${ranks}_rank"
      PROPERTIES
      LABELS "mpi;numerical;stage3;s3-c3;fast"
      TIMEOUT 300
      PROCESSORS "${ranks}"
      RESOURCE_LOCK hundun_stage3_mpi_m)
  endforeach()
endforeach()

foreach(ranks IN ITEMS 1 2)
  add_test(NAME "test_stage3_flow_models_fast_${ranks}_rank"
    COMMAND bash
            "${PROJECT_SOURCE_DIR}/tests/acceptance/stage3_flow_models_fast.sh"
            "$<TARGET_FILE:hundun>" "$<TARGET_FILE:write_task19a_stl>"
            "${MPIEXEC_EXECUTABLE}" "${ranks}")
  set_tests_properties("test_stage3_flow_models_fast_${ranks}_rank"
    PROPERTIES
    LABELS "acceptance;mpi;numerical;stage3;s3-c1;s3-d1;s3-c2;s3-d2;s3-c3;fast"
    TIMEOUT 300
    PROCESSORS "${ranks}"
    RESOURCE_LOCK hundun_stage3_mpi_m)
endforeach()
