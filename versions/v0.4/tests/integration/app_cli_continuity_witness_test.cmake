# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

foreach(required IN ITEMS PRODUCT CASE_JSON THERMOPHYSICS PROBE_ROOT SOURCE_ROOT)
  if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
    message(FATAL_ERROR "missing ${required}")
  endif()
endforeach()

get_filename_component(probe_root_absolute "${PROBE_ROOT}" ABSOLUTE)
get_filename_component(source_root_absolute "${SOURCE_ROOT}" ABSOLUTE)
if(probe_root_absolute STREQUAL "/" OR
   probe_root_absolute STREQUAL source_root_absolute OR
   NOT probe_root_absolute MATCHES "v04-app-cli-continuity-witness")
  message(FATAL_ERROR
    "refusing unsafe CLI witness probe root: ${probe_root_absolute}")
endif()
set(PROBE_ROOT "${probe_root_absolute}")

file(REMOVE_RECURSE "${PROBE_ROOT}")
file(MAKE_DIRECTORY "${PROBE_ROOT}/case")
configure_file("${CASE_JSON}" "${PROBE_ROOT}/case/case.json" COPYONLY)
configure_file("${THERMOPHYSICS}"
               "${PROBE_ROOT}/case/thermophysics.d" COPYONLY)

# Make the public CLI take a deterministic non-continuity failure path: the
# isothermal wall gives the first pressure solve a nontrivial RHS, while one
# allowed Krylov iteration is deliberately insufficient.  No terminal
# physical audit has run, so terminal residuals are unavailable and the
# report's continuity witness is invalid by contract.
file(READ "${PROBE_ROOT}/case/case.json" case_json)
set(original_case_json "${case_json}")
string(REPLACE
  [=["y_min": {
      "flow_kind": "symmetry",
      "thermal_kind": "none"]=]
  [=["y_min": {
      "flow_kind": "no_slip_wall",
      "thermal_kind": "isothermal_wall"]=]
  case_json "${case_json}")
string(REGEX REPLACE
  [=[("y_min"[^}]*"temperature": )0\.0]=]
  [=[\1400.0]=]
  case_json "${case_json}")
string(REPLACE [=["maximum_iterations": 400]=]
               [=["maximum_iterations": 1]=]
               case_json "${case_json}")
string(REPLACE [=["true_residual_interval": 4]=]
               [=["true_residual_interval": 1]=]
               case_json "${case_json}")
if(case_json STREQUAL original_case_json OR
   NOT case_json MATCHES [=["maximum_iterations": 1]=] OR
   NOT case_json MATCHES [=["thermal_kind": "isothermal_wall"]=] OR
   NOT case_json MATCHES [=["temperature": 400\.0]=])
  message(FATAL_ERROR "failed to construct the CLI failure fixture")
endif()
file(WRITE "${PROBE_ROOT}/case/case.json" "${case_json}")

execute_process(
  COMMAND "${PRODUCT}" run "${PROBE_ROOT}/case"
          --output "${PROBE_ROOT}/run"
          --steps 1 --output-interval 0 --restart-interval 0
  RESULT_VARIABLE run_status
  OUTPUT_VARIABLE run_output
  ERROR_VARIABLE run_error)

if(run_status EQUAL 0)
  message(FATAL_ERROR
    "CLI failure fixture unexpectedly completed: ${run_output}${run_error}")
endif()
if(NOT run_error MATCHES "failed_stage=44" OR
   NOT run_error MATCHES "terminal_audit=unavailable" OR
   NOT run_error MATCHES "pressure_calls=1" OR
   NOT run_error MATCHES "p1_iterations=1" OR
   NOT run_error MATCHES "rejected step detail=604")
  message(FATAL_ERROR
    "CLI did not reach the intended pressure-solve failure: ${run_error}")
endif()
if(run_error MATCHES " eos=" OR
   run_error MATCHES " continuity=" OR
   run_error MATCHES " energy=" OR
   run_error MATCHES " closed_mass=" OR
   run_error MATCHES " gauge=" OR
   run_error MATCHES "continuity_witness=")
  message(FATAL_ERROR
    "CLI exposed unavailable terminal audit data: ${run_error}")
endif()

# A post-commit file failure must not look like an unaccepted numerical step.
file(WRITE "${PROBE_ROOT}/case/case.json" "${original_case_json}")
file(MAKE_DIRECTORY "${PROBE_ROOT}/output-failure/evidence.jsonl")
execute_process(
  COMMAND "${PRODUCT}" run "${PROBE_ROOT}/case"
          --output "${PROBE_ROOT}/output-failure"
          --steps 1 --output-interval 0 --restart-interval 0
  RESULT_VARIABLE output_status
  OUTPUT_VARIABLE output_text
  ERROR_VARIABLE output_error)
if(output_status EQUAL 0 OR
   NOT output_error MATCHES "termination_phase=Evidence committed_step=1" OR
   output_error MATCHES "committed_time=0([ \\n]|$)" OR
   output_error MATCHES "failed_stage=")
  message(FATAL_ERROR
    "CLI lost the accepted step or mislabeled an output failure: ${output_error}")
endif()

# The actual CLI must compute the flow-based bound before its first attempt.
# On this dx=0.125 grid with U=10, Co=0.8 gives dt=0.01, not the input 0.1.
set(adaptive_json "${original_case_json}")
string(REPLACE "[0.1, 0.0, 0.0]" "[10.0, 0.0, 0.0]"
               adaptive_json "${adaptive_json}")
string(REPLACE "\"initial_dt\": 0.0001" "\"initial_dt\": 0.1"
               adaptive_json "${adaptive_json}")
string(REPLACE "\"maximum_dt\": 0.01" "\"maximum_dt\": 0.1"
               adaptive_json "${adaptive_json}")
string(REPLACE "\"maximum_retries\": 8" "\"maximum_retries\": 1"
               adaptive_json "${adaptive_json}")
file(WRITE "${PROBE_ROOT}/case/case.json" "${adaptive_json}")
execute_process(
  COMMAND "${PRODUCT}" run "${PROBE_ROOT}/case"
          --output "${PROBE_ROOT}/adaptive-flow"
          --steps 1 --output-interval 1 --restart-interval 0
  RESULT_VARIABLE adaptive_status
  OUTPUT_VARIABLE adaptive_output
  ERROR_VARIABLE adaptive_error)
if(NOT adaptive_status EQUAL 0)
  message(FATAL_ERROR "CLI flow-based dt failed: ${adaptive_output}${adaptive_error}")
endif()
file(READ "${PROBE_ROOT}/adaptive-flow/screen.log" adaptive_screen)
if(NOT adaptive_screen MATCHES "attempts=1" OR
   NOT adaptive_output MATCHES " time=0\\.01([ \\n]|$)")
  message(FATAL_ERROR "CLI did not use one dt=0.01 attempt: ${adaptive_output}${adaptive_screen}")
endif()

# Reject the same flow-based dt before any solve, and explain the actual
# proposal versus its configured lower bound instead of reporting only 454.
string(REGEX REPLACE [=["minimum_dt": [^,}]+]=]
                     [=["minimum_dt": 0.02]=]
                     minimum_json "${adaptive_json}")
file(WRITE "${PROBE_ROOT}/case/case.json" "${minimum_json}")
execute_process(
  COMMAND "${PRODUCT}" run "${PROBE_ROOT}/case"
          --output "${PROBE_ROOT}/minimum-dt"
          --steps 1 --output-interval 0 --restart-interval 0
  RESULT_VARIABLE minimum_status
  OUTPUT_VARIABLE minimum_output
  ERROR_VARIABLE minimum_error)
if(minimum_status EQUAL 0 OR
   NOT minimum_error MATCHES "rejected step detail=454" OR
   NOT minimum_error MATCHES "termination_phase=advance committed_step=0" OR
   NOT minimum_error MATCHES "attempts=0" OR
   NOT minimum_error MATCHES "initial_proposed_dt=0\\.01" OR
   NOT minimum_error MATCHES "minimum_dt=0\\.02" OR
   minimum_error MATCHES "failed_stage=")
  message(FATAL_ERROR "CLI lost pre-solve dt rejection diagnostics: ${minimum_error}")
endif()
