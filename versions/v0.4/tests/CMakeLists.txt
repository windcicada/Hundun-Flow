# SPDX-License-Identifier: Apache-2.0

add_executable(v04_core_status_test
  unit/core_status_test.cpp)
target_link_libraries(v04_core_status_test PRIVATE hundun_v04_core)

add_test(NAME v04_core_status COMMAND v04_core_status_test)

add_executable(v04_core_arena_test
  unit/core_arena_test.cpp)
target_link_libraries(v04_core_arena_test PRIVATE hundun_v04_core)

add_test(NAME v04_core_arena COMMAND v04_core_arena_test)

add_executable(v04_core_transaction_test
  unit/core_transaction_test.cpp)
target_link_libraries(v04_core_transaction_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)

add_test(NAME v04_core_transaction COMMAND v04_core_transaction_test)

add_executable(v04_core_transaction_mpi_test
  mpi/core_transaction_mpi_test.cpp)
target_link_libraries(v04_core_transaction_mpi_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)

foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_core_transaction_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_transaction_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_core_transaction_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}")
endforeach()

add_executable(v04_app_case_test
  unit/app_case_test.cpp)
target_link_libraries(v04_app_case_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)
target_include_directories(v04_app_case_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_app_case_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)

add_test(NAME v04_app_case COMMAND v04_app_case_test)
set_tests_properties(v04_app_case PROPERTIES
  ENVIRONMENT "HUNDUN_V04_TEST_DATA=${CMAKE_CURRENT_SOURCE_DIR}/data")

add_executable(v04_app_case_broadcast_test
  mpi/app_case_broadcast_test.cpp)
target_include_directories(v04_app_case_broadcast_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_app_case_broadcast_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_app_case_broadcast_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)

foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_app_case_broadcast_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_app_case_broadcast_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_app_case_broadcast_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}")
endforeach()
