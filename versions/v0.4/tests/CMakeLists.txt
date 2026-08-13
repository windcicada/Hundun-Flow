# SPDX-License-Identifier: Apache-2.0

add_executable(v04_core_status_test
  unit/core_status_test.cpp)
target_link_libraries(v04_core_status_test PRIVATE hundun_v04_core)

add_test(NAME v04_core_status COMMAND v04_core_status_test)

add_executable(v04_physics_thermo_test
  unit/physics_thermo_test.cpp)
target_link_libraries(v04_physics_thermo_test PRIVATE hundun_v04_core)
add_test(NAME v04_physics_thermo COMMAND v04_physics_thermo_test)

add_executable(v04_physics_thermophysical_compiler_test
  mpi/physics_thermophysical_compiler_test.cpp)
target_include_directories(v04_physics_thermophysical_compiler_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_link_libraries(v04_physics_thermophysical_compiler_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_physics_thermophysical_compiler_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_physics_thermophysical_compiler_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_physics_thermophysical_compiler_${rank_count}
    PROPERTIES
      PROCESSORS "${rank_count}"
      TIMEOUT 30)
endforeach()

add_executable(v04_physics_transport_test
  unit/physics_transport_test.cpp)
target_link_libraries(v04_physics_transport_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)
add_test(NAME v04_physics_transport COMMAND v04_physics_transport_test)
set_tests_properties(v04_physics_transport PROPERTIES TIMEOUT 30)

add_executable(v04_physics_contribution_test
  unit/physics_contribution_test.cpp)
target_link_libraries(v04_physics_contribution_test PRIVATE hundun_v04_core)
add_test(NAME v04_physics_contribution COMMAND v04_physics_contribution_test)

add_executable(v04_physics_closed_mass_test
  numerical/physics_closed_mass_test.cpp)
target_link_libraries(v04_physics_closed_mass_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_physics_closed_mass_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_physics_closed_mass_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_physics_closed_mass_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 30)
endforeach()

add_executable(v04_bc_compile_test
  unit/bc_compile_test.cpp)
target_link_libraries(v04_bc_compile_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)
add_test(NAME v04_bc_compile COMMAND v04_bc_compile_test)

add_executable(v04_bc_apply_test
  unit/bc_apply_test.cpp)
target_link_libraries(v04_bc_apply_test PRIVATE hundun_v04_core)
add_test(NAME v04_bc_apply COMMAND v04_bc_apply_test)

add_executable(v04_bc_nscbc_test
  unit/bc_nscbc_test.cpp)
target_link_libraries(v04_bc_nscbc_test PRIVATE hundun_v04_core)
add_test(NAME v04_bc_nscbc COMMAND v04_bc_nscbc_test)

add_executable(v04_bc_time_test
  unit/bc_time_test.cpp)
target_link_libraries(v04_bc_time_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_bc_time_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_bc_time_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_bc_time_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}")
endforeach()

add_executable(v04_bc_decomposition_mpi_test
  mpi/bc_decomposition_mpi_test.cpp)
target_link_libraries(v04_bc_decomposition_mpi_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)

foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_bc_decomposition_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_bc_decomposition_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_bc_decomposition_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}")
endforeach()

add_executable(v04_parallel_cpu_test
  unit/parallel_cpu_test.cpp)
target_link_libraries(v04_parallel_cpu_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)

foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_parallel_cpu_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_parallel_cpu_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_parallel_cpu_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}")
endforeach()

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

add_executable(v04_mesh_cartesian_test
  unit/mesh_cartesian_test.cpp)
target_link_libraries(v04_mesh_cartesian_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)
add_test(NAME v04_mesh_cartesian COMMAND v04_mesh_cartesian_test)

add_executable(v04_mesh_focus_test
  unit/mesh_focus_test.cpp)
target_link_libraries(v04_mesh_focus_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)
add_test(NAME v04_mesh_focus COMMAND v04_mesh_focus_test)

add_executable(v04_mesh_stl_scan_test
  unit/mesh_stl_scan_test.cpp)
target_include_directories(v04_mesh_stl_scan_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_mesh_stl_scan_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_mesh_stl_scan_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)
add_test(NAME v04_mesh_stl_scan COMMAND v04_mesh_stl_scan_test)
set_tests_properties(v04_mesh_stl_scan PROPERTIES
  ENVIRONMENT "HUNDUN_V04_TEST_DATA=${CMAKE_CURRENT_SOURCE_DIR}/data")

add_executable(v04_mesh_stl_scan_mpi_test
  mpi/mesh_stl_scan_mpi_test.cpp)
target_include_directories(v04_mesh_stl_scan_mpi_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_mesh_stl_scan_mpi_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_mesh_stl_scan_mpi_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)

foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_mesh_stl_scan_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_mesh_stl_scan_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_mesh_stl_scan_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    ENVIRONMENT "HUNDUN_V04_TEST_DATA=${CMAKE_CURRENT_SOURCE_DIR}/data")
endforeach()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_executable(v04_mesh_stl_scan_release_test
    performance/mesh_stl_scan_release_test.cpp)
  target_include_directories(v04_mesh_stl_scan_release_test PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../src")
  target_compile_definitions(v04_mesh_stl_scan_release_test PRIVATE
    HUNDUN_V04_ENABLE_TEST_ACCESS=1)
  target_link_libraries(v04_mesh_stl_scan_release_test PRIVATE
    hundun_v04_core
    MPI::MPI_C)
  add_test(NAME v04_mesh_stl_scan_release
    COMMAND v04_mesh_stl_scan_release_test)
  set_tests_properties(v04_mesh_stl_scan_release PROPERTIES
    ENVIRONMENT "HUNDUN_V04_TEST_DATA=${CMAKE_CURRENT_SOURCE_DIR}/data"
    LABELS "release;performance"
    TIMEOUT 60)
endif()

add_executable(v04_mesh_decomposition_mpi_test
  mpi/mesh_decomposition_mpi_test.cpp)
target_link_libraries(v04_mesh_decomposition_mpi_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)

foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_mesh_decomposition_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_mesh_decomposition_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_mesh_decomposition_mpi_${rank_count} PROPERTIES
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

add_executable(v04_parallel_halo_mpi_test
  mpi/parallel_halo_mpi_test.cpp)
target_include_directories(v04_parallel_halo_mpi_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_parallel_halo_mpi_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_parallel_halo_mpi_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)

add_executable(v04_parallel_halo_failure_test
  mpi/parallel_halo_failure_test.cpp)
target_include_directories(v04_parallel_halo_failure_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_parallel_halo_failure_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_parallel_halo_failure_test PRIVATE
  hundun_v04_core
  MPI::MPI_C)

foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_parallel_halo_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_parallel_halo_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_parallel_halo_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}")

  add_test(NAME v04_parallel_halo_failure_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_parallel_halo_failure_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_parallel_halo_failure_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}")
endforeach()
