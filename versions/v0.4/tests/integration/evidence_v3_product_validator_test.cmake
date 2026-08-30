# SPDX-License-Identifier: Apache-2.0

foreach(required IN ITEMS PRODUCT PYTHON VALIDATOR PROBE_ROOT)
  if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
    message(FATAL_ERROR "missing ${required}")
  endif()
endforeach()

set(v5_product_row_check [=[
import json
import math
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
scenario = sys.argv[2]
expected = {
    "fresh": ((1, 1, 1, True, False), (2, 2, 2, False, False)),
    "restart": ((3, 1, 1, False, True), (4, 2, 2, False, False)),
}
if scenario not in expected:
    raise SystemExit("unknown V5 producer scenario: " + scenario)
try:
    rows = [json.loads(line) for line in path.read_text(
        encoding="utf-8").splitlines() if line.strip()]
except (OSError, ValueError) as error:
    raise SystemExit("cannot read V5 producer evidence: " + str(error))
if len(rows) != len(expected[scenario]):
    raise SystemExit("{} evidence must contain exactly {} rows; found {}".format(
        scenario, len(expected[scenario]), len(rows)))
for index, (row, wanted) in enumerate(zip(rows, expected[scenario]), 1):
    step, requested, executed, startup, restart_recovery = wanted
    if not isinstance(row, dict):
        raise SystemExit("{} row {} is not an object".format(scenario, index))
    if row.get("schema") != "HUNDUN_V04_EVIDENCE_V5":
        raise SystemExit("{} row {} is not V5".format(scenario, index))
    observed = (row.get("step"), row.get("requested_bdf_order"),
                row.get("bdf_order"), row.get("startup"),
                row.get("restart_recovery"))
    if observed != wanted:
        raise SystemExit("{} row {} lifecycle {} != {}".format(
            scenario, index, observed, wanted))
    if row.get("temporal_method_fallback") is not False:
        raise SystemExit("{} step {} unexpectedly used temporal fallback".format(
            scenario, step))
    if row.get("pressure_solve_contract") != "continuity_energy_coupled":
        raise SystemExit("{} step {} lost coupled pressure contract".format(
            scenario, step))
    if (row.get("pressure_energy_refinement_termination") !=
            "component_residuals_converged"):
        raise SystemExit("{} step {} did not close refinement".format(
            scenario, step))
    refinement_count = row.get("pressure_energy_refinement_solve_calls")
    refinements = row.get("pressure_energy_refinement")
    if (not isinstance(refinement_count, int) or
            isinstance(refinement_count, bool) or refinement_count < 0 or
            not isinstance(refinements, list) or
            len(refinements) != refinement_count):
        raise SystemExit("{} step {} has invalid refinement prefix".format(
            scenario, step))
    terminal = row.get("terminal_physical_audit")
    if (not isinstance(terminal, dict) or terminal.get("present") is not True or
            not isinstance(terminal.get("final_flux_revision"), int) or
            isinstance(terminal.get("final_flux_revision"), bool) or
            terminal["final_flux_revision"] <= 0):
        raise SystemExit("{} step {} lacks terminal physical audit".format(
            scenario, step))
    cfl = terminal.get("committed_convective_cfl")
    if (not isinstance(cfl, dict) or
            any(isinstance(cfl.get(field), bool) or
                not isinstance(cfl.get(field), (int, float)) or
                not math.isfinite(cfl.get(field))
                for field in ("out_max", "abs_max", "limit")) or
            cfl["out_max"] < 0 or cfl["abs_max"] < 0 or cfl["limit"] <= 0 or
            cfl["out_max"] > cfl["limit"] * (1 + 64 * sys.float_info.epsilon)):
        raise SystemExit("{} step {} has invalid committed CFL".format(
            scenario, step))
    limiter = row.get("momentum_predictor_limiter")
    if (not isinstance(limiter, dict) or
            limiter.get("scheme") != "common_face_afc_v2" or
            "theta" in limiter or "activations" in limiter):
        raise SystemExit("{} step {} lost V5 AFC provenance".format(
            scenario, step))
    advective = limiter.get("advective_cfl")
    if not isinstance(advective, dict) or advective.get("present") is not True:
        raise SystemExit("{} step {} lacks advective CFL evidence".format(
            scenario, step))
    for field in ("plan", "time_revision", "density_revision",
                  "face_flux_revision"):
        value = advective.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise SystemExit("{} step {} has invalid advective {}".format(
                scenario, step, field))
    activity = advective.get("activity_collective")
    if (isinstance(activity, bool) or not isinstance(activity, int) or
            activity < 0):
        raise SystemExit("{} step {} has invalid advective activity".format(
            scenario, step))
    if advective["face_flux_revision"] == terminal["final_flux_revision"]:
        raise SystemExit("{} step {} reused the final CFL revision".format(
            scenario, step))
    for field in ("dt", "out_max", "abs_max", "limit"):
        value = advective.get(field)
        if (isinstance(value, bool) or
                not isinstance(value, (int, float)) or
                not math.isfinite(value)):
            raise SystemExit("{} step {} has invalid advective {}".format(
                scenario, step, field))
    if (advective["dt"] <= 0 or advective["out_max"] < 0 or
            advective["abs_max"] < 0 or advective["limit"] <= 0 or
            advective["out_max"] >
                advective["limit"] * (1 + 64 * sys.float_info.epsilon)):
        raise SystemExit("{} step {} has inadmissible advective CFL".format(
            scenario, step))
]=])

function(assert_v5_product_rows evidence_path scenario)
  execute_process(
    COMMAND "${PYTHON}" -c "${v5_product_row_check}"
            "${evidence_path}" "${scenario}"
    RESULT_VARIABLE row_status
    OUTPUT_VARIABLE row_output
    ERROR_VARIABLE row_error)
  if(NOT row_status EQUAL 0)
    message(FATAL_ERROR
      "V5 ${scenario} row semantics failed (${row_status}): "
      "${row_output}${row_error}")
  endif()
endfunction()

file(REMOVE_RECURSE "${PROBE_ROOT}")

execute_process(
  COMMAND "${PRODUCT}" init-case --output "${PROBE_ROOT}/case"
  RESULT_VARIABLE init_status
  OUTPUT_VARIABLE init_output
  ERROR_VARIABLE init_error)
if(NOT init_status EQUAL 0)
  message(FATAL_ERROR
    "V5 producer init failed (${init_status}): ${init_output}${init_error}")
endif()

execute_process(
  COMMAND "${PRODUCT}" run "${PROBE_ROOT}/case"
          --output "${PROBE_ROOT}/run"
          --steps 2 --output-interval 0 --restart-interval 1
  RESULT_VARIABLE run_status
  OUTPUT_VARIABLE run_output
  ERROR_VARIABLE run_error)
if(NOT run_status EQUAL 0)
  message(FATAL_ERROR
    "V5 producer run failed (${run_status}): ${run_output}${run_error}")
endif()

assert_v5_product_rows("${PROBE_ROOT}/run/evidence.jsonl" fresh)

execute_process(
  COMMAND "${PYTHON}" "${VALIDATOR}" runtime
          "${PROBE_ROOT}/run/evidence.jsonl"
  RESULT_VARIABLE validation_status
  OUTPUT_VARIABLE validation_output
  ERROR_VARIABLE validation_error)
if(NOT validation_status EQUAL 0)
  message(FATAL_ERROR
    "V5 producer evidence rejected (${validation_status}): "
    "${validation_output}${validation_error}")
endif()

execute_process(
  COMMAND "${PRODUCT}" run "${PROBE_ROOT}/case"
          --output "${PROBE_ROOT}/run-resumed"
          --restart "${PROBE_ROOT}/run/Restart"
          --steps 2 --output-interval 0 --restart-interval 0
  RESULT_VARIABLE resumed_status
  OUTPUT_VARIABLE resumed_output
  ERROR_VARIABLE resumed_error)
if(NOT resumed_status EQUAL 0)
  message(FATAL_ERROR
    "V5 resumed producer run failed (${resumed_status}): "
    "${resumed_output}${resumed_error}")
endif()

assert_v5_product_rows("${PROBE_ROOT}/run-resumed/evidence.jsonl" restart)

execute_process(
  COMMAND "${PYTHON}" "${VALIDATOR}" runtime
          "${PROBE_ROOT}/run-resumed/evidence.jsonl"
  RESULT_VARIABLE resumed_validation_status
  OUTPUT_VARIABLE resumed_validation_output
  ERROR_VARIABLE resumed_validation_error)
if(NOT resumed_validation_status EQUAL 0)
  message(FATAL_ERROR
    "V5 resumed evidence rejected (${resumed_validation_status}): "
    "${resumed_validation_output}${resumed_validation_error}")
endif()
