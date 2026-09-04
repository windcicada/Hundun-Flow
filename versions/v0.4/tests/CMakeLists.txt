# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

find_program(HUNDUN_V04_PYTHON_EXECUTABLE
  NAMES python3.6 python3
  REQUIRED)

add_test(NAME v04_thin_domain_runner_self_test
  COMMAND $<TARGET_FILE:v04_thin_domain_runner> --self-test)
set_tests_properties(v04_thin_domain_runner_self_test PROPERTIES TIMEOUT 30)

add_test(NAME v04_evidence_workflow
  COMMAND "${HUNDUN_V04_PYTHON_EXECUTABLE}"
          "${hundun_repository_root}/tools/v04_evidence_validate.py"
          self-test)
add_test(NAME v04_evidence_v4_product_validator
  COMMAND "${CMAKE_COMMAND}"
          "-DPRODUCT=$<TARGET_FILE:hundun>"
          "-DPYTHON=${HUNDUN_V04_PYTHON_EXECUTABLE}"
          "-DVALIDATOR=${hundun_repository_root}/tools/v04_evidence_validate.py"
          "-DPROBE_ROOT=${CMAKE_CURRENT_BINARY_DIR}/v04-evidence-v4-product-validator"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/integration/evidence_v3_product_validator_test.cmake")
set_tests_properties(v04_evidence_v4_product_validator PROPERTIES TIMEOUT 120)
add_test(NAME v04_app_cli_continuity_witness
  COMMAND "${CMAKE_COMMAND}"
          "-DPRODUCT=$<TARGET_FILE:hundun>"
          "-DCASE_JSON=${CMAKE_CURRENT_SOURCE_DIR}/data/case_minimal_valid.json"
          "-DTHERMOPHYSICS=${CMAKE_CURRENT_SOURCE_DIR}/data/thermophysics.d"
          "-DSOURCE_ROOT=${hundun_repository_root}"
          "-DPROBE_ROOT=${CMAKE_CURRENT_BINARY_DIR}/v04-app-cli-continuity-witness"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/integration/app_cli_continuity_witness_test.cmake")
set_tests_properties(v04_app_cli_continuity_witness PROPERTIES TIMEOUT 120)
add_test(NAME v04_evidence_performance_policy
  COMMAND "${HUNDUN_V04_PYTHON_EXECUTABLE}"
          "${hundun_repository_root}/tools/v04_evidence_validate.py"
          policy
          "${hundun_repository_root}/docs/verification/v0.4-re3900-performance-policy.json")
add_test(NAME v04_literature_workflow
  COMMAND "${HUNDUN_V04_PYTHON_EXECUTABLE}"
          "${hundun_repository_root}/tools/v04_literature_extract.py"
          self-test)
add_test(NAME v04_literature_partial_receipt
  COMMAND "${HUNDUN_V04_PYTHON_EXECUTABLE}"
          "${hundun_repository_root}/tools/v04_literature_extract.py"
          receipt-validate
          "${hundun_repository_root}/docs/verification/v0.4-literature-data-receipt-r4-partial.json")

add_executable(v04_core_status_test
  unit/core_status_test.cpp)
target_link_libraries(v04_core_status_test PRIVATE hundun_v04_test_core)

add_test(NAME v04_core_status COMMAND v04_core_status_test)

add_executable(v04_app_mpi_runtime_test
  unit/app_mpi_runtime_test.cpp)
target_link_libraries(v04_app_mpi_runtime_test PRIVATE hundun_v04_test_core)
add_test(NAME v04_app_mpi_runtime COMMAND v04_app_mpi_runtime_test)

add_executable(v04_public_headers_test
  unit/public_headers_test.cpp)
target_link_libraries(v04_public_headers_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_public_headers COMMAND v04_public_headers_test)

add_executable(v04_mesh_ibm_quadratic_test
  unit/mesh_ibm_quadratic_test.cpp)
target_link_libraries(v04_mesh_ibm_quadratic_test PRIVATE hundun_v04_test_core)
add_test(NAME v04_mesh_ibm_quadratic COMMAND v04_mesh_ibm_quadratic_test)

add_executable(v04_mesh_ibm_surface_test
  unit/mesh_ibm_surface_test.cpp)
target_link_libraries(v04_mesh_ibm_surface_test PRIVATE hundun_v04_test_core)
add_test(NAME v04_mesh_ibm_surface COMMAND v04_mesh_ibm_surface_test)

add_executable(v04_mesh_ibm_order_test
  numerical/mesh_ibm_order_test.cpp)
target_link_libraries(v04_mesh_ibm_order_test PRIVATE hundun_v04_test_core)
add_test(NAME v04_mesh_ibm_order COMMAND v04_mesh_ibm_order_test)

add_executable(v04_solver_low_mach_mms_test
  numerical/solver_low_mach_mms_test.cpp)
target_link_libraries(v04_solver_low_mach_mms_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_low_mach_mms
  COMMAND v04_solver_low_mach_mms_test)
set_tests_properties(v04_solver_low_mach_mms PROPERTIES TIMEOUT 120)

add_executable(v04_solver_enthalpy_terms_test
  numerical/solver_enthalpy_terms_test.cpp)
target_link_libraries(v04_solver_enthalpy_terms_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_enthalpy_terms
  COMMAND v04_solver_enthalpy_terms_test)
set_tests_properties(v04_solver_enthalpy_terms PROPERTIES TIMEOUT 120)

add_executable(v04_solver_source_stage_test
  numerical/solver_source_stage_test.cpp)
target_link_libraries(v04_solver_source_stage_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_source_stage
  COMMAND v04_solver_source_stage_test)
set_tests_properties(v04_solver_source_stage PROPERTIES TIMEOUT 60)

add_executable(v04_solver_species_conservation_test
  numerical/solver_species_conservation_test.cpp)
target_link_libraries(v04_solver_species_conservation_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_species_conservation
  COMMAND v04_solver_species_conservation_test)
set_tests_properties(v04_solver_species_conservation PROPERTIES TIMEOUT 120)

add_executable(v04_solver_pressure_storage_test
  numerical/solver_pressure_storage_test.cpp)
target_link_libraries(v04_solver_pressure_storage_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_pressure_storage
  COMMAND v04_solver_pressure_storage_test)
set_tests_properties(v04_solver_pressure_storage PROPERTIES TIMEOUT 120)

add_executable(v04_solver_closed_gauge_prepare_mpi_test
  mpi/solver_closed_gauge_prepare_mpi_test.cpp)
target_compile_definitions(v04_solver_closed_gauge_prepare_mpi_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_solver_closed_gauge_prepare_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_closed_gauge_prepare_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_closed_gauge_prepare_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_closed_gauge_prepare_mpi_${rank_count}
    PROPERTIES PROCESSORS "${rank_count}" TIMEOUT 120)
endforeach()

add_executable(v04_solver_thermophysical_coupling_test
  numerical/solver_thermophysical_coupling_test.cpp)
target_link_libraries(v04_solver_thermophysical_coupling_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_thermophysical_coupling
  COMMAND v04_solver_thermophysical_coupling_test)
set_tests_properties(v04_solver_thermophysical_coupling PROPERTIES TIMEOUT 120)

add_executable(v04_solver_pressure_energy_schur_test
  numerical/solver_pressure_energy_schur_test.cpp)
target_link_libraries(v04_solver_pressure_energy_schur_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_pressure_energy_schur
  COMMAND v04_solver_pressure_energy_schur_test)
set_tests_properties(v04_solver_pressure_energy_schur PROPERTIES TIMEOUT 120)

add_executable(v04_solver_pressure_energy_globalization_test
  numerical/solver_pressure_energy_globalization_test.cpp)
target_link_libraries(v04_solver_pressure_energy_globalization_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_pressure_energy_globalization
  COMMAND v04_solver_pressure_energy_globalization_test)
set_tests_properties(v04_solver_pressure_energy_globalization PROPERTIES
  TIMEOUT 60)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_pressure_energy_globalization_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_pressure_energy_globalization_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    v04_solver_pressure_energy_globalization_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}" TIMEOUT 60)
endforeach()

add_executable(v04_solver_thermophysical_idp_test
  numerical/solver_thermophysical_idp_test.cpp)
target_link_libraries(v04_solver_thermophysical_idp_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_thermophysical_idp
  COMMAND v04_solver_thermophysical_idp_test)
set_tests_properties(v04_solver_thermophysical_idp PROPERTIES TIMEOUT 120)

add_test(NAME v04_solver_thermophysical_idp_mpi_2
  COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
          ${MPIEXEC_PREFLAGS}
          "$<TARGET_FILE:v04_solver_thermophysical_idp_test>"
          ${MPIEXEC_POSTFLAGS})
set_tests_properties(v04_solver_thermophysical_idp_mpi_2 PROPERTIES
  PROCESSORS 2
  TIMEOUT 120)

add_executable(v04_solver_thermophysical_ghost_authority_mpi_test
  mpi/solver_thermophysical_ghost_authority_mpi_test.cpp)
target_link_libraries(v04_solver_thermophysical_ghost_authority_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_thermophysical_ghost_authority_mpi_2
  COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
          ${MPIEXEC_PREFLAGS}
          "$<TARGET_FILE:v04_solver_thermophysical_ghost_authority_mpi_test>"
          ${MPIEXEC_POSTFLAGS})
set_tests_properties(v04_solver_thermophysical_ghost_authority_mpi_2 PROPERTIES
  PROCESSORS 2
  TIMEOUT 120)

add_executable(v04_solver_piso_authority_test
  unit/solver_piso_authority_test.cpp)
target_link_libraries(v04_solver_piso_authority_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_piso_authority
  COMMAND v04_solver_piso_authority_test)
set_tests_properties(v04_solver_piso_authority PROPERTIES TIMEOUT 60)

add_executable(v04_solver_pressure_boundary_plan_test
  unit/solver_pressure_boundary_plan_test.cpp)
target_link_libraries(v04_solver_pressure_boundary_plan_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_pressure_boundary_plan
  COMMAND v04_solver_pressure_boundary_plan_test)
set_tests_properties(v04_solver_pressure_boundary_plan PROPERTIES TIMEOUT 60)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_pressure_boundary_plan_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_pressure_boundary_plan_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_pressure_boundary_plan_mpi_${rank_count}
    PROPERTIES PROCESSORS "${rank_count}" TIMEOUT 120)
endforeach()

add_executable(v04_solver_initialization_projection_test
  unit/solver_initialization_projection_test.cpp)
target_link_libraries(v04_solver_initialization_projection_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_initialization_projection
  COMMAND v04_solver_initialization_projection_test)
set_tests_properties(v04_solver_initialization_projection PROPERTIES
  TIMEOUT 120)

add_executable(v04_solver_initialization_projection_mpi_test
  mpi/solver_initialization_projection_mpi_test.cpp)
target_link_libraries(v04_solver_initialization_projection_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_initialization_projection_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_initialization_projection_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_initialization_projection_mpi_${rank_count}
    PROPERTIES PROCESSORS "${rank_count}" TIMEOUT 120)
endforeach()

add_executable(v04_solver_piso_mutation_test
  unit/solver_piso_mutation_test.cpp)
target_link_libraries(v04_solver_piso_mutation_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_piso_mutation
  COMMAND v04_solver_piso_mutation_test)
set_tests_properties(v04_solver_piso_mutation PROPERTIES TIMEOUT 60)

add_executable(v04_solver_piso_checkerboard_test
  numerical/solver_piso_checkerboard_test.cpp)
target_link_libraries(v04_solver_piso_checkerboard_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_piso_checkerboard
  COMMAND v04_solver_piso_checkerboard_test)
set_tests_properties(v04_solver_piso_checkerboard PROPERTIES TIMEOUT 60)
foreach(rank_count IN ITEMS 2 4)
  add_test(NAME v04_solver_piso_checkerboard_mpi_${rank_count}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${rank_count}
            ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:v04_solver_piso_checkerboard_test>
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_piso_checkerboard_mpi_${rank_count}
    PROPERTIES PROCESSORS "${rank_count}" TIMEOUT 60)
endforeach()

add_executable(v04_solver_piso_temporal_order_test
  numerical/solver_piso_temporal_order_test.cpp)
target_link_libraries(v04_solver_piso_temporal_order_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_piso_temporal_order
  COMMAND v04_solver_piso_temporal_order_test)
set_tests_properties(v04_solver_piso_temporal_order PROPERTIES TIMEOUT 60)

add_executable(v04_solver_ibm_pressure_operator_test
  unit/solver_ibm_pressure_operator_test.cpp)
target_link_libraries(v04_solver_ibm_pressure_operator_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_ibm_pressure_operator
  COMMAND v04_solver_ibm_pressure_operator_test)
set_tests_properties(v04_solver_ibm_pressure_operator PROPERTIES TIMEOUT 120)

add_executable(v04_solver_ibm_equation_interface_test
  unit/solver_ibm_equation_interface_test.cpp)
target_link_libraries(v04_solver_ibm_equation_interface_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_ibm_equation_interface
  COMMAND v04_solver_ibm_equation_interface_test)
set_tests_properties(v04_solver_ibm_equation_interface PROPERTIES TIMEOUT 120)

add_executable(v04_solver_ibm_force_oracle_test
  numerical/solver_ibm_force_oracle_test.cpp)
target_link_libraries(v04_solver_ibm_force_oracle_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_ibm_force_oracle
  COMMAND v04_solver_ibm_force_oracle_test)
set_tests_properties(v04_solver_ibm_force_oracle PROPERTIES TIMEOUT 120)

add_executable(v04_solver_ibm_force_mutation_test
  numerical/solver_ibm_force_mutation_test.cpp)
target_link_libraries(v04_solver_ibm_force_mutation_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_ibm_force_mutation
  COMMAND v04_solver_ibm_force_mutation_test)
set_tests_properties(v04_solver_ibm_force_mutation PROPERTIES TIMEOUT 120)

add_executable(v04_solver_equations_mpi_test
  mpi/solver_equations_mpi_test.cpp)
target_link_libraries(v04_solver_equations_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_equations_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_equations_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_equations_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
endforeach()

add_executable(v04_solver_piso_mpi_test
  mpi/solver_piso_mpi_test.cpp)
target_include_directories(v04_solver_piso_mpi_test PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_compile_definitions(v04_solver_piso_mpi_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_solver_piso_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_piso_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_piso_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_piso_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 180)
endforeach()

add_executable(v04_solver_ibm_physical_boundary_flux_mpi_test
  mpi/solver_ibm_physical_boundary_flux_mpi_test.cpp)
target_include_directories(v04_solver_ibm_physical_boundary_flux_mpi_test PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_compile_definitions(v04_solver_ibm_physical_boundary_flux_mpi_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_solver_ibm_physical_boundary_flux_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_ibm_physical_boundary_flux_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_ibm_physical_boundary_flux_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    v04_solver_ibm_physical_boundary_flux_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 180)
endforeach()

add_executable(v04_solver_ibm_force_mpi_test
  mpi/solver_ibm_force_mpi_test.cpp)
target_link_libraries(v04_solver_ibm_force_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_ibm_force_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_ibm_force_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_ibm_force_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 180)
endforeach()

add_executable(v04_mesh_ibm_stencil_mpi_test
  mpi/mesh_ibm_stencil_mpi_test.cpp)
target_link_libraries(v04_mesh_ibm_stencil_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_mesh_ibm_stencil_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_mesh_ibm_stencil_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_mesh_ibm_stencil_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
endforeach()

add_executable(v04_mesh_ibm_plan_test
  integration/mesh_ibm_plan_test.cpp)
target_link_libraries(v04_mesh_ibm_plan_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_mesh_ibm_plan_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_mesh_ibm_plan_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_mesh_ibm_plan_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
endforeach()

add_executable(v04_core_execution_graph_test
  unit/core_execution_graph_test.cpp)
target_link_libraries(v04_core_execution_graph_test PRIVATE hundun_v04_test_core)
target_compile_definitions(v04_core_execution_graph_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
add_test(NAME v04_core_execution_graph COMMAND v04_core_execution_graph_test)

add_executable(v04_core_product_freeze_test
  integration/core_product_freeze_test.cpp)
target_link_libraries(v04_core_product_freeze_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_core_product_freeze_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_core_product_freeze_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
add_test(NAME v04_core_product_freeze COMMAND v04_core_product_freeze_test)
set_tests_properties(v04_core_product_freeze PROPERTIES TIMEOUT 120)

add_executable(v04_core_hot_resource_test
  integration/core_hot_resource_test.cpp)
target_link_libraries(v04_core_hot_resource_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_compile_definitions(v04_core_hot_resource_test PRIVATE
  HUNDUN_V04_SOURCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/..")
add_test(NAME v04_core_hot_resource COMMAND v04_core_hot_resource_test)
set_tests_properties(v04_core_hot_resource PROPERTIES TIMEOUT 120)

add_executable(v04_core_product_freeze_mpi_test
  mpi/core_product_freeze_mpi_test.cpp)
target_link_libraries(v04_core_product_freeze_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_core_product_freeze_mpi_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_core_product_freeze_mpi_test PRIVATE
  HUNDUN_V04_SOURCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/.."
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_core_product_freeze_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_product_freeze_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_core_product_freeze_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 180)
  add_test(NAME v04_core_temporal_fallback_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_product_freeze_mpi_test>"
            "--temporal-fallback-only"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_core_temporal_fallback_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
  add_test(NAME v04_core_candidate_globalization_red_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_product_freeze_mpi_test>"
            "--candidate-globalization-red-only"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    v04_core_candidate_globalization_red_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
  add_test(NAME v04_core_multispecies_open_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_product_freeze_mpi_test>"
            "--multispecies-open-only"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    v04_core_multispecies_open_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 180)
  add_test(NAME v04_core_fresh_open_eos_flux_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_product_freeze_mpi_test>"
            "--fresh-open-eos-flux-only"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    v04_core_fresh_open_eos_flux_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
  add_test(NAME v04_core_fresh_initialize_contract_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_product_freeze_mpi_test>"
            "--fresh-initialize-contract-only"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    v04_core_fresh_initialize_contract_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
  add_test(NAME v04_core_restart_allocation_contract_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_product_freeze_mpi_test>"
            "--restart-allocation-contract-only"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    v04_core_restart_allocation_contract_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
  add_test(NAME v04_core_unsupported_boundary_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_product_freeze_mpi_test>"
            "--unsupported-boundary-only"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    v04_core_unsupported_boundary_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
endforeach()

foreach(rank_count IN ITEMS 2 4)
  add_test(NAME v04_core_product_localized_immersed_compile_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_product_freeze_mpi_test>"
            "--localized-immersed-compile-only"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    v04_core_product_localized_immersed_compile_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 60)
  add_test(
    NAME v04_core_product_localized_immersed_bind_fail_close_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_core_product_freeze_mpi_test>"
            "--localized-immersed-bind-fail-close-only"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    v04_core_product_localized_immersed_bind_fail_close_mpi_${rank_count}
    PROPERTIES PROCESSORS "${rank_count}" TIMEOUT 60)
endforeach()

add_executable(v04_io_restart_mpi_test
  mpi/io_restart_mpi_test.cpp)
target_link_libraries(v04_io_restart_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_io_restart_mpi_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_io_restart_mpi_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
add_test(NAME v04_io_restart_mpi_4
  COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "4"
          ${MPIEXEC_PREFLAGS}
          "$<TARGET_FILE:v04_io_restart_mpi_test>"
          ${MPIEXEC_POSTFLAGS})
set_tests_properties(v04_io_restart_mpi_4 PROPERTIES
  PROCESSORS 4
  TIMEOUT 180)

add_executable(v04_io_product_path_test
  integration/io_product_path_test.cpp)
target_link_libraries(v04_io_product_path_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_io_product_path COMMAND v04_io_product_path_test)
set_tests_properties(v04_io_product_path PROPERTIES TIMEOUT 120)

add_executable(v04_app_init_case_test
  integration/app_init_case_test.cpp)
target_link_libraries(v04_app_init_case_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_compile_definitions(v04_app_init_case_test PRIVATE
  HUNDUN_V04_SOURCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/..")
add_test(NAME v04_app_init_case COMMAND v04_app_init_case_test)
set_tests_properties(v04_app_init_case PROPERTIES TIMEOUT 120)

add_executable(v04_app_driver_test
  integration/app_driver_test.cpp)
target_link_libraries(v04_app_driver_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_app_driver_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_app_driver_test PRIVATE
  HUNDUN_V04_SOURCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/.."
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
add_test(NAME v04_app_driver COMMAND v04_app_driver_test)
set_tests_properties(v04_app_driver PROPERTIES TIMEOUT 120)

add_executable(v04_full_half_restart_compare
  integration/full_half_restart_compare.cpp)
target_link_libraries(v04_full_half_restart_compare PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_compile_definitions(v04_full_half_restart_compare PRIVATE
  HUNDUN_V04_SOURCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/..")
add_test(NAME v04_full_half_restart_compare_selftest
  COMMAND v04_full_half_restart_compare --self-test)
set_tests_properties(v04_full_half_restart_compare_selftest
  PROPERTIES TIMEOUT 180)

add_test(NAME v04_restart_rank_change_continue
  COMMAND "${CMAKE_COMMAND}"
          "-DPRODUCT=$<TARGET_FILE:hundun>"
          "-DCOMPARATOR=$<TARGET_FILE:v04_full_half_restart_compare>"
          "-DPYTHON=${HUNDUN_V04_PYTHON_EXECUTABLE}"
          "-DVALIDATOR=${hundun_repository_root}/tools/v04_evidence_validate.py"
          "-DPROBE_ROOT=${CMAKE_CURRENT_BINARY_DIR}/v04-restart-rank-change"
          "-DSOURCE_ROOT=${CMAKE_CURRENT_SOURCE_DIR}/.."
          "-DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}"
          "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
          "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
          "-DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/integration/restart_rank_change_continue_test.cmake")
set_tests_properties(v04_restart_rank_change_continue PROPERTIES
  PROCESSORS 4
  TIMEOUT 300)

add_executable(v04_product_fresh_projection_test
  integration/product_fresh_projection_test.cpp)
target_link_libraries(v04_product_fresh_projection_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_product_fresh_projection_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_product_fresh_projection_test PRIVATE
  HUNDUN_V04_SOURCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/.."
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
add_test(NAME v04_product_fresh_projection
  COMMAND v04_product_fresh_projection_test)
set_tests_properties(v04_product_fresh_projection PROPERTIES TIMEOUT 180)

add_executable(v04_product_fresh_projection_mpi_test
  mpi/product_fresh_projection_mpi_test.cpp)
target_link_libraries(v04_product_fresh_projection_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_product_fresh_projection_mpi_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_product_fresh_projection_mpi_test PRIVATE
  HUNDUN_V04_SOURCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/.."
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_product_fresh_projection_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_product_fresh_projection_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_product_fresh_projection_mpi_${rank_count}
    PROPERTIES PROCESSORS "${rank_count}" TIMEOUT 180)
endforeach()

add_executable(v04_product_pressure_energy_retry_mpi_test
  mpi/product_pressure_energy_retry_mpi_test.cpp)
target_link_libraries(v04_product_pressure_energy_retry_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_product_pressure_energy_retry_mpi_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_product_pressure_energy_retry_mpi_test PRIVATE
  HUNDUN_V04_SOURCE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/.."
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_product_pressure_energy_retry_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_product_pressure_energy_retry_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_product_pressure_energy_retry_mpi_${rank_count}
    PROPERTIES PROCESSORS "${rank_count}" TIMEOUT 180)
endforeach()

add_executable(v04_product_pressure_energy_temporal_convergence_test
  integration/product_pressure_energy_temporal_convergence_test.cpp)
target_link_libraries(v04_product_pressure_energy_temporal_convergence_test
  PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_product_pressure_energy_temporal_convergence
  COMMAND v04_product_pressure_energy_temporal_convergence_test)
set_tests_properties(v04_product_pressure_energy_temporal_convergence
  PROPERTIES TIMEOUT 180)

add_executable(v04_core_resource_contract_test
  unit/core_resource_contract_test.cpp)
target_link_libraries(v04_core_resource_contract_test PRIVATE hundun_v04_test_core)
add_test(NAME v04_core_resource_contract COMMAND v04_core_resource_contract_test)

add_executable(v04_solver_linear_lifecycle_test
  unit/solver_linear_lifecycle_test.cpp)
target_link_libraries(v04_solver_linear_lifecycle_test PRIVATE
  hundun_v04_test_core)
add_test(NAME v04_solver_linear_lifecycle
  COMMAND v04_solver_linear_lifecycle_test)

add_executable(v04_solver_mg_reuse_test
  unit/solver_mg_reuse_test.cpp)
target_link_libraries(v04_solver_mg_reuse_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_solver_mg_reuse_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_solver_mg_reuse_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
add_test(NAME v04_solver_mg_reuse COMMAND v04_solver_mg_reuse_test)

add_executable(v04_solver_mg_chebyshev_test
  unit/solver_mg_chebyshev_test.cpp)
target_link_libraries(v04_solver_mg_chebyshev_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_solver_mg_chebyshev_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_solver_mg_chebyshev_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
add_test(NAME v04_solver_mg_chebyshev COMMAND v04_solver_mg_chebyshev_test)

add_executable(v04_solver_mg_convergence_test
  numerical/solver_mg_convergence_test.cpp)
target_link_libraries(v04_solver_mg_convergence_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_mg_convergence
  COMMAND v04_solver_mg_convergence_test)
set_tests_properties(v04_solver_mg_convergence PROPERTIES TIMEOUT 120)

add_executable(v04_solver_mg_coarse_use_test
  numerical/solver_mg_coarse_use_test.cpp)
target_link_libraries(v04_solver_mg_coarse_use_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_mg_coarse_use
  COMMAND v04_solver_mg_coarse_use_test)
set_tests_properties(v04_solver_mg_coarse_use PROPERTIES TIMEOUT 120)

add_executable(v04_solver_mg_line_relaxation_test
  numerical/solver_mg_line_relaxation_test.cpp)
target_link_libraries(v04_solver_mg_line_relaxation_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_solver_mg_line_relaxation_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_solver_mg_line_relaxation_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
add_test(NAME v04_solver_mg_line_relaxation
  COMMAND v04_solver_mg_line_relaxation_test)
set_tests_properties(v04_solver_mg_line_relaxation PROPERTIES TIMEOUT 120)

add_executable(v04_solver_mg_krylov_isolation_test
  integration/solver_mg_krylov_isolation_test.cpp)
target_link_libraries(v04_solver_mg_krylov_isolation_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_mg_krylov_isolation
  COMMAND v04_solver_mg_krylov_isolation_test)

add_executable(v04_solver_krylov_basis_update_test
  unit/solver_krylov_basis_update_test.cpp)
target_link_libraries(v04_solver_krylov_basis_update_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_solver_krylov_basis_update_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_solver_krylov_basis_update_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
add_test(NAME v04_solver_krylov_basis_update
  COMMAND v04_solver_krylov_basis_update_test)

add_executable(v04_solver_hypre_isolation_test
  integration/solver_hypre_isolation_test.cpp)
target_link_libraries(v04_solver_hypre_isolation_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_solver_hypre_isolation_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_solver_hypre_isolation_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
if(HUNDUN_ENABLE_HYPRE)
  target_compile_definitions(v04_solver_hypre_isolation_test PRIVATE
    HUNDUN_V04_TEST_HAVE_HYPRE=1)
  target_link_libraries(v04_solver_hypre_isolation_test PRIVATE HYPRE::HYPRE)
endif()
add_test(NAME v04_solver_hypre_isolation
  COMMAND v04_solver_hypre_isolation_test)
foreach(rank_count IN ITEMS 2 4)
  add_test(NAME v04_solver_hypre_communicator_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_hypre_isolation_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_hypre_communicator_mpi_${rank_count}
    PROPERTIES
      PROCESSORS "${rank_count}"
      TIMEOUT 60)
endforeach()

add_executable(v04_solver_mg_mpi_test
  mpi/solver_mg_mpi_test.cpp)
target_link_libraries(v04_solver_mg_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_solver_mg_mpi_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_solver_mg_mpi_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_mg_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_mg_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_mg_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
endforeach()

add_executable(v04_solver_mg_update_contract_mpi_test
  mpi/solver_mg_update_contract_mpi_test.cpp)
target_link_libraries(v04_solver_mg_update_contract_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 2 4)
  add_test(NAME v04_solver_mg_update_contract_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_mg_update_contract_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_mg_update_contract_mpi_${rank_count}
    PROPERTIES
      PROCESSORS "${rank_count}"
      TIMEOUT 60)
endforeach()

add_executable(v04_solver_mg_odd_partition_mpi_test
  mpi/solver_mg_odd_partition_mpi_test.cpp)
target_link_libraries(v04_solver_mg_odd_partition_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 2 4)
  add_test(NAME v04_solver_mg_odd_partition_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_mg_odd_partition_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_mg_odd_partition_mpi_${rank_count}
    PROPERTIES
      PROCESSORS "${rank_count}"
      TIMEOUT 120)
endforeach()

add_executable(v04_solver_krylov_mpi_test
  mpi/solver_krylov_mpi_test.cpp)
target_link_libraries(v04_solver_krylov_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
target_include_directories(v04_solver_krylov_mpi_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_solver_krylov_mpi_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)

add_executable(v04_solver_reduction_mpi_test
  mpi/solver_reduction_mpi_test.cpp)
target_link_libraries(v04_solver_reduction_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)

foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_krylov_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_krylov_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_krylov_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 90)

  add_test(NAME v04_solver_reduction_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_reduction_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_reduction_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 60)
endforeach()

add_executable(v04_physics_thermo_test
  unit/physics_thermo_test.cpp)
target_link_libraries(v04_physics_thermo_test PRIVATE hundun_v04_test_core)
add_test(NAME v04_physics_thermo COMMAND v04_physics_thermo_test)

add_executable(v04_physics_thermophysical_compiler_test
  mpi/physics_thermophysical_compiler_test.cpp)
target_include_directories(v04_physics_thermophysical_compiler_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_link_libraries(v04_physics_thermophysical_compiler_test PRIVATE
  hundun_v04_test_core
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
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_physics_transport COMMAND v04_physics_transport_test)
set_tests_properties(v04_physics_transport PROPERTIES TIMEOUT 30)

add_executable(v04_physics_contribution_test
  unit/physics_contribution_test.cpp)
target_link_libraries(v04_physics_contribution_test PRIVATE hundun_v04_test_core)
add_test(NAME v04_physics_contribution COMMAND v04_physics_contribution_test)

add_executable(v04_physics_wale_test
  unit/physics_wale_test.cpp)
target_link_libraries(v04_physics_wale_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_physics_wale COMMAND v04_physics_wale_test)

add_executable(v04_physics_vreman_test
  unit/physics_vreman_test.cpp)
target_link_libraries(v04_physics_vreman_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_physics_vreman COMMAND v04_physics_vreman_test)

add_executable(v04_physics_wall_function_test
  unit/physics_wall_function_test.cpp)
target_link_libraries(v04_physics_wall_function_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_physics_wall_function COMMAND v04_physics_wall_function_test)

add_executable(v04_physics_turbulence_mpi_test
  mpi/physics_turbulence_mpi_test.cpp)
target_link_libraries(v04_physics_turbulence_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_physics_turbulence_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_physics_turbulence_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_physics_turbulence_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 120)
endforeach()

add_executable(v04_solver_manufactured_cartesian_test
  numerical/solver_manufactured_cartesian_test.cpp)
target_link_libraries(v04_solver_manufactured_cartesian_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_manufactured_cartesian
  COMMAND v04_solver_manufactured_cartesian_test)
set_tests_properties(v04_solver_manufactured_cartesian PROPERTIES TIMEOUT 60)

add_executable(v04_solver_conservation_test
  numerical/solver_conservation_test.cpp)
target_link_libraries(v04_solver_conservation_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_solver_conservation COMMAND v04_solver_conservation_test)
set_tests_properties(v04_solver_conservation PROPERTIES TIMEOUT 60)

add_executable(v04_solver_decomposition_mpi_test
  mpi/solver_decomposition_mpi_test.cpp)
target_link_libraries(v04_solver_decomposition_mpi_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
foreach(rank_count IN ITEMS 1 2 4)
  add_test(NAME v04_solver_decomposition_mpi_${rank_count}
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}"
            ${MPIEXEC_PREFLAGS}
            "$<TARGET_FILE:v04_solver_decomposition_mpi_test>"
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(v04_solver_decomposition_mpi_${rank_count} PROPERTIES
    PROCESSORS "${rank_count}"
    TIMEOUT 60)
endforeach()

add_executable(v04_physics_closed_mass_test
  numerical/physics_closed_mass_test.cpp)
target_link_libraries(v04_physics_closed_mass_test PRIVATE
  hundun_v04_test_core
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
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_bc_compile COMMAND v04_bc_compile_test)

add_executable(v04_bc_apply_test
  unit/bc_apply_test.cpp)
target_link_libraries(v04_bc_apply_test PRIVATE hundun_v04_test_core)
add_test(NAME v04_bc_apply COMMAND v04_bc_apply_test)

add_executable(v04_bc_thermo_test
  unit/bc_thermo_test.cpp)
target_link_libraries(v04_bc_thermo_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_bc_thermo COMMAND v04_bc_thermo_test)
add_test(NAME v04_bc_thermo_mpi_2
  COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2"
          ${MPIEXEC_PREFLAGS}
          "$<TARGET_FILE:v04_bc_thermo_test>"
          ${MPIEXEC_POSTFLAGS})
set_tests_properties(v04_bc_thermo_mpi_2 PROPERTIES
  PROCESSORS 2
  TIMEOUT 60)

add_executable(v04_bc_nscbc_test
  unit/bc_nscbc_test.cpp)
target_link_libraries(v04_bc_nscbc_test PRIVATE hundun_v04_test_core)
add_test(NAME v04_bc_nscbc COMMAND v04_bc_nscbc_test)

add_executable(v04_bc_time_test
  unit/bc_time_test.cpp)
target_link_libraries(v04_bc_time_test PRIVATE
  hundun_v04_test_core
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
  hundun_v04_test_core
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
  hundun_v04_test_core
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
target_link_libraries(v04_core_arena_test PRIVATE hundun_v04_test_core)

add_test(NAME v04_core_arena COMMAND v04_core_arena_test)

add_executable(v04_core_transaction_test
  unit/core_transaction_test.cpp)
target_compile_definitions(v04_core_transaction_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_core_transaction_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)

add_test(NAME v04_core_transaction COMMAND v04_core_transaction_test)

add_executable(v04_core_transaction_mpi_test
  mpi/core_transaction_mpi_test.cpp)
target_compile_definitions(v04_core_transaction_mpi_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_core_transaction_mpi_test PRIVATE
  hundun_v04_test_core
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
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_mesh_cartesian COMMAND v04_mesh_cartesian_test)

add_executable(v04_mesh_focus_test
  unit/mesh_focus_test.cpp)
target_link_libraries(v04_mesh_focus_test PRIVATE
  hundun_v04_test_core
  MPI::MPI_C)
add_test(NAME v04_mesh_focus COMMAND v04_mesh_focus_test)

add_executable(v04_mesh_stl_scan_test
  unit/mesh_stl_scan_test.cpp)
target_include_directories(v04_mesh_stl_scan_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_mesh_stl_scan_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_mesh_stl_scan_test PRIVATE
  hundun_v04_test_core
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
  hundun_v04_test_core
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
    hundun_v04_test_core
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
  hundun_v04_test_core
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
  hundun_v04_test_core
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
  hundun_v04_test_core
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
  hundun_v04_test_core
  MPI::MPI_C)

add_executable(v04_parallel_halo_failure_test
  mpi/parallel_halo_failure_test.cpp)
target_include_directories(v04_parallel_halo_failure_test PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_compile_definitions(v04_parallel_halo_failure_test PRIVATE
  HUNDUN_V04_ENABLE_TEST_ACCESS=1)
target_link_libraries(v04_parallel_halo_failure_test PRIVATE
  hundun_v04_test_core
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

# The release selector is data, not a regex assembled by the runner.  Validate
# every entry at configure time and attach one stable label for execution.
include("${CMAKE_CURRENT_LIST_DIR}/focused_manifest.cmake")
foreach(focused_test IN LISTS HUNDUN_V04_FOCUSED_TESTS)
  if(NOT TEST "${focused_test}")
    message(FATAL_ERROR
      "v0.4 focused manifest names an unregistered test: ${focused_test}")
  endif()
  set_property(TEST "${focused_test}" APPEND PROPERTY LABELS v04_focused)
endforeach()
foreach(sanitizer_test IN LISTS HUNDUN_V04_FOCUSED_SANITIZER_TESTS)
  if(NOT sanitizer_test IN_LIST HUNDUN_V04_FOCUSED_TESTS)
    message(FATAL_ERROR
      "v0.4 sanitizer test is absent from focused manifest: ${sanitizer_test}")
  endif()
endforeach()
