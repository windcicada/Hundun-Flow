# SPDX-License-Identifier: Apache-2.0

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
# continuity audit has run, so the report's witness is invalid by contract.
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
   NOT run_error MATCHES "continuity=0" OR
   NOT run_error MATCHES "pressure_calls=1" OR
   NOT run_error MATCHES "p1_iterations=1" OR
   NOT run_error MATCHES "rejected step detail=604")
  message(FATAL_ERROR
    "CLI did not reach the intended pressure-solve failure: ${run_error}")
endif()
if(run_error MATCHES "continuity_witness=")
  message(FATAL_ERROR
    "CLI exposed an invalid continuity witness: ${run_error}")
endif()
