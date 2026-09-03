# SPDX-License-Identifier: Apache-2.0

foreach(required IN ITEMS PRODUCT PYTHON VALIDATOR PROBE_ROOT)
  if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
    message(FATAL_ERROR "missing ${required}")
  endif()
endforeach()

set(v8_product_row_check [=[
import hashlib
import json
import math
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
scenario = sys.argv[2]
expected = {
    "fresh": ((1, 1, 1, True, False), (2, 2, 2, False, False)),
    "restart": ((3, 2, 2, False, False), (4, 2, 2, False, False)),
}
if scenario not in expected:
    raise SystemExit("unknown V8 producer scenario: " + scenario)
try:
    rows = [json.loads(line) for line in path.read_text(
        encoding="utf-8").splitlines() if line.strip()]
except (OSError, ValueError) as error:
    raise SystemExit("cannot read V8 producer evidence: " + str(error))
if len(rows) != len(expected[scenario]):
    raise SystemExit("{} evidence must contain exactly {} rows; found {}".format(
        scenario, len(expected[scenario]), len(rows)))
for index, (row, wanted) in enumerate(zip(rows, expected[scenario]), 1):
    step, requested, executed, startup, restart_recovery = wanted
    if not isinstance(row, dict):
        raise SystemExit("{} row {} is not an object".format(scenario, index))
    if row.get("schema") != "HUNDUN_V04_EVIDENCE_V8":
        raise SystemExit("{} row {} is not V8".format(scenario, index))
    observed = (row.get("step"), row.get("requested_bdf_order"),
                row.get("bdf_order"), row.get("startup"),
                row.get("restart_recovery"))
    if observed != wanted:
        raise SystemExit("{} row {} lifecycle {} != {}".format(
            scenario, index, observed, wanted))
    anchor = row.get("run_start")
    if not isinstance(anchor, dict):
        raise SystemExit("{} row {} lacks run-start anchor".format(
            scenario, index))
    if scenario == "fresh":
        expected_anchor = {
            "kind": "fresh", "previous_step": 0,
            "previous_time": 0, "restart_manifest_sha256": None,
        }
        if anchor != expected_anchor:
            raise SystemExit("fresh row {} has false anchor {}".format(
                index, anchor))
    else:
        restart_root = path.parent.parent / "run" / "Restart"
        current = (restart_root / "current").read_text(
            encoding="utf-8").strip()
        manifest = restart_root / current / "manifest.bin"
        expected_manifest = hashlib.sha256(manifest.read_bytes()).hexdigest()
        if (anchor.get("kind") != "restart" or
                anchor.get("previous_step") != 2 or
                anchor.get("previous_time") != rows[0].get(
                    "previous_committed_time") or
                anchor.get("restart_manifest_sha256") != expected_manifest):
            raise SystemExit("restart row {} has false snapshot anchor {}".format(
                index, anchor))
    if index > 1 and anchor != rows[index - 2].get("run_start"):
        raise SystemExit("{} row {} changes run-start anchor".format(
            scenario, index))
    if row.get("temporal_method_fallback") is not False:
        raise SystemExit("{} step {} unexpectedly used temporal fallback".format(
            scenario, step))
    if row.get("pressure_solve_contract") != "continuity_energy_coupled":
        raise SystemExit("{} step {} lost coupled pressure contract".format(
            scenario, step))
    if (row.get("coupling") != "PISO" or
            row.get("momentum_predictor_passes") != 1):
        raise SystemExit("{} step {} has invalid PISO coupling evidence".format(
            scenario, step))
    identity = row.get("candidate_identity")
    if (not isinstance(identity, dict) or
            identity.get("schema") !=
                "HUNDUN_V04_RUNTIME_CANDIDATE_IDENTITY_V2" or
            any(not isinstance(identity.get(field), str)
                for field in ("head", "tree", "build_manifest_sha256",
                              "executable_sha256", "identity_sha256"))):
        raise SystemExit("{} step {} lacks candidate identity".format(
            scenario, step))
    payload = (
        "schema=HUNDUN_V04_RUNTIME_CANDIDATE_IDENTITY_V2\n"
        "evidence_schema=HUNDUN_V04_EVIDENCE_V8\n"
        "head={}\ntree={}\nbuild_manifest_sha256={}\n"
        "executable_sha256={}\n").format(
            identity["head"], identity["tree"],
            identity["build_manifest_sha256"],
            identity["executable_sha256"])
    if (hashlib.sha256(payload.encode()).hexdigest() !=
            identity["identity_sha256"] or
            row.get("build") != int(identity["build_manifest_sha256"][:16], 16) or
            row.get("binary") != int(identity["executable_sha256"][:16], 16)):
        raise SystemExit("{} step {} has invalid candidate identity".format(
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
    if (not isinstance(cfl, dict) or cfl.get("valid") is not True or
            cfl.get("final_flux_revision") != terminal["final_flux_revision"] or
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
            limiter.get("scheme") != "common_face_afc_v3_owner" or
            "theta" in limiter or "activations" in limiter):
        raise SystemExit("{} step {} lost V6 AFC provenance".format(
            scenario, step))
    applicability = limiter.get("correction_metrics_applicability")
    active = limiter.get("active_correction_faces")
    limited_faces = limiter.get("limited_faces")
    if applicability == "not_applicable":
        if (active != 0 or limited_faces != 0 or
                limiter.get("retained_correction_l1_ratio") is not None or
                limiter.get("minimum_face_alpha") is not None or
                limiter.get("limited_face_fraction") is not None):
            raise SystemExit("{} step {} fabricates AFC metrics".format(
                scenario, step))
    elif applicability == "applicable":
        if (not isinstance(active, int) or active <= 0 or
                not isinstance(limited_faces, int) or
                limited_faces < 0 or limited_faces > active):
            raise SystemExit("{} step {} has invalid AFC counts".format(
                scenario, step))
    else:
        raise SystemExit("{} step {} lacks AFC applicability".format(
            scenario, step))
    advective = limiter.get("advective_cfl")
    if not isinstance(advective, dict) or advective.get("present") is not True:
        raise SystemExit("{} step {} lacks advective CFL evidence".format(
            scenario, step))
    for field in ("plan", "time_revision_collective",
                  "density_view_collective", "face_flux_view_collective"):
        value = advective.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise SystemExit("{} step {} has invalid advective {}".format(
                scenario, step, field))
    activity = advective.get("activity_collective")
    if (isinstance(activity, bool) or not isinstance(activity, int) or
            activity < 0):
        raise SystemExit("{} step {} has invalid advective activity".format(
            scenario, step))
    face_view = cfl.get("face_flux_view")
    if (not isinstance(face_view, dict) or
            not isinstance(face_view.get("collective"), int) or
            face_view["collective"] <= 0 or
            advective["face_flux_view_collective"] ==
                face_view["collective"]):
        raise SystemExit("{} step {} reused the final CFL view".format(
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
    previous = row.get("previous_committed_time")
    current = row.get("time")
    if (isinstance(previous, bool) or
            not isinstance(previous, (int, float)) or
            not math.isfinite(previous) or
            abs((current - previous) - advective["dt"]) >
                128 * sys.float_info.epsilon * max(
                    1.0, abs(current), abs(previous), abs(advective["dt"]))):
        raise SystemExit("{} step {} has false first/previous time".format(
            scenario, step))
    if index > 1 and previous != rows[index - 2]["time"]:
        raise SystemExit("{} step {} loses adjacent time authority".format(
            scenario, step))
    if index > 1 and identity != rows[index - 2]["candidate_identity"]:
        raise SystemExit("{} step {} changes immutable identity".format(
            scenario, step))
]=])

function(assert_v8_product_rows evidence_path scenario)
  execute_process(
    COMMAND "${PYTHON}" -c "${v8_product_row_check}"
            "${evidence_path}" "${scenario}"
    RESULT_VARIABLE row_status
    OUTPUT_VARIABLE row_output
    ERROR_VARIABLE row_error)
  if(NOT row_status EQUAL 0)
    message(FATAL_ERROR
      "V8 ${scenario} row semantics failed (${row_status}): "
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
    "V8 producer init failed (${init_status}): ${init_output}${init_error}")
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
    "V8 producer run failed (${run_status}): ${run_output}${run_error}")
endif()

assert_v8_product_rows("${PROBE_ROOT}/run/evidence.jsonl" fresh)

execute_process(
  COMMAND "${PYTHON}" "${VALIDATOR}" runtime
          "${PROBE_ROOT}/run/evidence.jsonl"
  RESULT_VARIABLE validation_status
  OUTPUT_VARIABLE validation_output
  ERROR_VARIABLE validation_error)
if(NOT validation_status EQUAL 0)
  message(FATAL_ERROR
    "V8 producer evidence rejected (${validation_status}): "
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
    "V8 resumed producer run failed (${resumed_status}): "
    "${resumed_output}${resumed_error}")
endif()

assert_v8_product_rows("${PROBE_ROOT}/run-resumed/evidence.jsonl" restart)

file(READ "${PROBE_ROOT}/run/Restart/current" restart_generation)
string(STRIP "${restart_generation}" restart_generation)
set(restart_manifest
  "${PROBE_ROOT}/run/Restart/${restart_generation}/manifest.bin")
execute_process(
  COMMAND "${PYTHON}" "${VALIDATOR}" runtime
          "${PROBE_ROOT}/run-resumed/evidence.jsonl"
          --run-start-manifest "${restart_manifest}"
  RESULT_VARIABLE resumed_validation_status
  OUTPUT_VARIABLE resumed_validation_output
  ERROR_VARIABLE resumed_validation_error)
if(NOT resumed_validation_status EQUAL 0)
  message(FATAL_ERROR
    "V8 resumed evidence rejected (${resumed_validation_status}): "
    "${resumed_validation_output}${resumed_validation_error}")
endif()
