# SPDX-License-Identifier: Apache-2.0

foreach(required IN ITEMS
    PRODUCT COMPARATOR PYTHON VALIDATOR PROBE_ROOT SOURCE_ROOT
    MPIEXEC_EXECUTABLE MPIEXEC_NUMPROC_FLAG)
  if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
    message(FATAL_ERROR "missing ${required}")
  endif()
endforeach()

foreach(required_file IN ITEMS PRODUCT COMPARATOR PYTHON VALIDATOR)
  if(NOT EXISTS "${${required_file}}")
    message(FATAL_ERROR "${required_file} does not exist: ${${required_file}}")
  endif()
endforeach()

get_filename_component(probe_root_absolute "${PROBE_ROOT}" ABSOLUTE)
get_filename_component(source_root_absolute "${SOURCE_ROOT}" ABSOLUTE)
if(probe_root_absolute STREQUAL "/" OR
   probe_root_absolute STREQUAL source_root_absolute OR
   NOT probe_root_absolute MATCHES "restart-rank-change")
  message(FATAL_ERROR
    "refusing unsafe rank-change probe root: ${probe_root_absolute}")
endif()
set(PROBE_ROOT "${probe_root_absolute}")

function(run_checked description)
  execute_process(
    COMMAND ${ARGN}
    RESULT_VARIABLE command_status
    OUTPUT_VARIABLE command_output
    ERROR_VARIABLE command_error)
  if(NOT command_status EQUAL 0)
    message(FATAL_ERROR
      "${description} failed (${command_status}):\n"
      "stdout:\n${command_output}\n"
      "stderr:\n${command_error}")
  endif()
endfunction()

function(run_mpi description rank_count)
  set(command
    "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${rank_count}")
  if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
    list(APPEND command ${MPIEXEC_PREFLAGS})
  endif()
  list(APPEND command ${ARGN})
  if(DEFINED MPIEXEC_POSTFLAGS AND NOT "${MPIEXEC_POSTFLAGS}" STREQUAL "")
    list(APPEND command ${MPIEXEC_POSTFLAGS})
  endif()
  run_checked("${description}" ${command})
endfunction()

set(case_root "${PROBE_ROOT}/case")
set(seed_root "${PROBE_ROOT}/seed-mpi2")
set(control_root "${PROBE_ROOT}/control-mpi2")
set(rank_change_root "${PROBE_ROOT}/rank-change-mpi4")
set(comparison_json "${PROBE_ROOT}/equivalent-endpoints.json")
set(fixed_dt "0.000125")
set(maximum_relative_rms "1e-11")
set(maximum_relative_linf "1e-10")

file(REMOVE_RECURSE "${PROBE_ROOT}")
run_checked("rank-change fixture init"
  "${PRODUCT}" init-case --output "${case_root}")
run_checked("rank-change IBM surface copy"
  "${CMAKE_COMMAND}" -E copy_if_different
  "${SOURCE_ROOT}/tests/data/cylinder_ascii.stl"
  "${case_root}/cylinder_ascii.stl")

set(configure_case [=[
import json
import pathlib
import sys

case_path = pathlib.Path(sys.argv[1])
fixed_dt = float(sys.argv[2])
document = json.loads(case_path.read_text(encoding="utf-8"))

mesh = document.get("mesh")
flow = document.get("flow")
time = document.get("time")
boundaries = document.get("boundaries")
if not isinstance(mesh, dict) or not isinstance(flow, dict) or \
        not isinstance(time, dict) or not isinstance(boundaries, dict):
    raise SystemExit("init-case JSON lacks the frozen mesh/flow/time/boundary objects")
if (mesh.get("kind") != "uniform" or
        mesh.get("domain") != {"lower": [0, 0, 0], "upper": [1, 1, 1]} or
        mesh.get("exact_cells") != [8, 8, 8] or
        mesh.get("minimum_spacing") != [0.125, 0.125, 0.125] or
        mesh.get("data_files") != [] or
        mesh.get("immersed_boundary") is not None):
    raise SystemExit("init-case mesh template changed, refusing an unreviewed mutation")
if flow.get("pressure_reference") != "closed_mass":
    raise SystemExit("init-case pressure-reference template changed")
if set(boundaries) != {"x_min", "x_max", "y_min", "y_max", "z_min", "z_max"}:
    raise SystemExit("init-case boundary catalog changed")
if any(face.get("flow_kind") != "periodic" for face in boundaries.values()):
    raise SystemExit("init-case boundary template changed")
if (time.get("scheme") != "variable_bdf2" or
        time.get("initial_dt") != 0.001 or
        time.get("minimum_dt") != 1e-8 or
        time.get("maximum_dt") != 0.1):
    raise SystemExit("init-case time template changed, refusing an unreviewed mutation")

mesh["domain"] = {"lower": [-2.0, -2.0, -2.0],
                  "upper": [2.0, 2.0, 2.0]}
mesh["exact_cells"] = [16, 16, 16]
mesh["minimum_spacing"] = [0.25, 0.25, 0.25]
mesh["limits"]["max_global_cells"] = 4096
mesh["limits"]["max_memory_bytes_per_rank"] = 1073741824
mesh["immersed_boundary"] = {
    "stl_file": "cylinder_ascii.stl",
    "fluid_side": "outside",
}
flow["pressure_reference"] = "boundary_absolute"
document["turbulence"]["model"] = "none"

directions = {
    "x_min": [1.0, 0.0, 0.0],
    "x_max": [1.0, 0.0, 0.0],
    "y_min": [0.0, -1.0, 0.0],
    "y_max": [0.0, 1.0, 0.0],
    "z_min": [0.0, 0.0, -1.0],
    "z_max": [0.0, 0.0, 1.0],
}
for name, face in boundaries.items():
    face["flow_kind"] = "symmetry"
    face["thermal_kind"] = "none"
    face["velocity"] = [0.0, 0.0, 0.0]
    face["direction"] = directions[name]
    face["backflow_velocity"] = [0.0, 0.0, 0.0]
    face["mass_flow_rate"] = 0.0
    face["pressure"] = 101325.0
    face["temperature"] = 300.0
    face["backflow_temperature"] = 300.0
    face["allow_backflow"] = False

boundaries["x_min"]["flow_kind"] = "velocity_inlet"
boundaries["x_min"]["velocity"] = [0.1, 0.0, 0.0]
boundaries["x_max"]["flow_kind"] = "pressure_outlet"

time["initial_dt"] = fixed_dt
time["minimum_dt"] = fixed_dt
time["maximum_dt"] = fixed_dt

temporary = case_path.with_suffix(".json.pending")
temporary.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
temporary.replace(case_path)
]=])
run_checked("rank-change case specialization"
  "${PYTHON}" -c "${configure_case}" "${case_root}/case.json" "${fixed_dt}")

# The seed is one real two-rank Product run.  Both continuations below read
# this exact Restart generation; no uninterrupted or different-time endpoint
# is used as an oracle.
run_mpi("rank-change seed MPI2" 2
  "${PRODUCT}" run "${case_root}"
  --output "${seed_root}"
  --steps 2 --output-interval 0 --restart-interval 2)

run_mpi("same-rank Restart control MPI2" 2
  "${PRODUCT}" run "${case_root}"
  --output "${control_root}"
  --restart "${seed_root}/Restart"
  --steps 2 --output-interval 0 --restart-interval 2)

run_mpi("rank-change Restart continuation MPI4" 4
  "${PRODUCT}" run "${case_root}"
  --output "${rank_change_root}"
  --restart "${seed_root}/Restart"
  --steps 2 --output-interval 0 --restart-interval 2)

file(READ "${seed_root}/Restart/current" seed_generation)
string(STRIP "${seed_generation}" seed_generation)
set(seed_manifest
  "${seed_root}/Restart/${seed_generation}/manifest.bin")
foreach(branch IN ITEMS control-mpi2 rank-change-mpi4)
  run_checked("${branch} V6 runtime evidence validation"
    "${PYTHON}" "${VALIDATOR}" runtime
    "${PROBE_ROOT}/${branch}/evidence.jsonl"
    --run-start-manifest "${seed_manifest}")
endforeach()

set(assert_evidence [=[
import json
import math
import pathlib
import sys

fixed_dt = float(sys.argv[1])
paths = [pathlib.Path(value) for value in sys.argv[2:]]
expected = (
    (3, 1, 1, True),
    (4, 2, 2, False),
)
all_rows = []
for path in paths:
    rows = [json.loads(line) for line in path.read_text(
        encoding="utf-8").splitlines() if line.strip()]
    if len(rows) != 2:
        raise SystemExit(f"{path}: expected exactly two committed rows, found {len(rows)}")
    revisions = []
    for index, (row, wanted) in enumerate(zip(rows, expected), 1):
        step, requested, executed, restart_recovery = wanted
        observed = (
            row.get("step"), row.get("requested_bdf_order"),
            row.get("bdf_order"), row.get("restart_recovery"),
        )
        if row.get("schema") != "HUNDUN_V04_EVIDENCE_V6" or observed != wanted:
            raise SystemExit(f"{path}: row {index} lifecycle {observed} != {wanted}")
        if (row.get("startup") is not False or row.get("retry") is not False or
                row.get("temporal_method_fallback") is not False or
                row.get("thermophysical_predictor_calls") != 1):
            raise SystemExit(f"{path}: step {step} has a fallback/retry/startup lifecycle")
        if row.get("pressure_solve_contract") != "continuity_energy_coupled":
            raise SystemExit(f"{path}: step {step} lost the coupled pressure-energy contract")
        refinement_count = row.get("pressure_energy_refinement_solve_calls")
        refinements = row.get("pressure_energy_refinement")
        if (isinstance(refinement_count, bool) or
                not isinstance(refinement_count, int) or
                refinement_count < 0 or refinement_count > 6 or
                not isinstance(refinements, list) or
                len(refinements) != refinement_count or
                row.get("pressure_energy_refinement_termination") !=
                    "component_residuals_converged"):
            raise SystemExit(f"{path}: step {step} has invalid refinement closure")
        for ordinal, refinement in enumerate(refinements, 1):
            if (refinement.get("ordinal") != ordinal or
                    not isinstance(refinement.get("target_generation"), int) or
                    refinement["target_generation"] <= 0 or
                    not isinstance(refinement.get("collective_lineage"), int) or
                    refinement["collective_lineage"] <= 0):
                raise SystemExit(
                    f"{path}: step {step} refinement {ordinal} lost typed provenance")
        accepted_time = row.get("time")
        expected_time = step * fixed_dt
        time_tolerance = 128.0 * sys.float_info.epsilon * max(1.0, abs(expected_time))
        if (isinstance(accepted_time, bool) or
                not isinstance(accepted_time, (int, float)) or
                not math.isfinite(accepted_time) or
                abs(accepted_time - expected_time) > time_tolerance):
            raise SystemExit(
                f"{path}: step {step} time {accepted_time} is not the fixed-dt target {expected_time}")
        previous_time = row.get("previous_committed_time")
        if (not isinstance(previous_time, (int, float)) or
                isinstance(previous_time, bool) or
                abs((accepted_time - previous_time) - fixed_dt) > time_tolerance):
            raise SystemExit(
                f"{path}: step {step} loses previous committed time authority")
        terminal = row.get("terminal_physical_audit")
        if not isinstance(terminal, dict) or terminal.get("present") is not True:
            raise SystemExit(f"{path}: step {step} lacks the terminal physical audit")
        revision = terminal.get("final_flux_revision")
        if isinstance(revision, bool) or not isinstance(revision, int) or revision <= 0:
            raise SystemExit(f"{path}: step {step} has no final-flux authority")
        revisions.append(revision)
        limiter = row.get("momentum_predictor_limiter")
        if (not isinstance(limiter, dict) or
                limiter.get("scheme") != "common_face_afc_v3_owner" or
                "theta" in limiter or "activations" in limiter):
            raise SystemExit(f"{path}: step {step} lost V6 AFC provenance")
        advective = limiter.get("advective_cfl")
        if not isinstance(advective, dict) or advective.get("present") is not True:
            raise SystemExit(f"{path}: step {step} lacks advective CFL evidence")
        for field in ("plan", "time_revision_collective",
                      "density_view_collective",
                      "face_flux_view_collective"):
            value = advective.get(field)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise SystemExit(
                    f"{path}: step {step} advective CFL {field} is invalid")
        activity = advective.get("activity_collective")
        if (isinstance(activity, bool) or not isinstance(activity, int) or
                activity < 0):
            raise SystemExit(
                f"{path}: step {step} advective CFL activity is invalid")
        cfl = terminal.get("committed_convective_cfl")
        face_view = cfl.get("face_flux_view") if isinstance(cfl, dict) else None
        if (not isinstance(face_view, dict) or
                not isinstance(face_view.get("collective"), int) or
                face_view["collective"] <= 0 or
                advective["face_flux_view_collective"] ==
                    face_view["collective"]):
            raise SystemExit(
                f"{path}: step {step} advective and final CFL views coincide")
        for field in ("dt", "out_max", "abs_max", "limit"):
            value = advective.get(field)
            if (isinstance(value, bool) or
                    not isinstance(value, (int, float)) or
                    not math.isfinite(value)):
                raise SystemExit(
                    f"{path}: step {step} advective CFL {field} is invalid")
        if (advective["dt"] <= 0.0 or advective["out_max"] < 0.0 or
                advective["abs_max"] < 0.0 or advective["limit"] <= 0.0 or
                advective["out_max"] >
                    advective["limit"] *
                        (1.0 + 64.0 * sys.float_info.epsilon)):
            raise SystemExit(f"{path}: step {step} advective CFL failed")
        if (not isinstance(cfl, dict) or cfl.get("valid") is not True or
                cfl.get("final_flux_revision") != revision):
            raise SystemExit(f"{path}: step {step} lacks committed CFL evidence")
        for field in ("out_max", "abs_max", "limit"):
            value = cfl.get(field)
            if (isinstance(value, bool) or
                    not isinstance(value, (int, float)) or
                    not math.isfinite(value)):
                raise SystemExit(
                    f"{path}: step {step} committed CFL {field} is invalid")
        if (cfl["out_max"] < 0.0 or cfl["abs_max"] < 0.0 or
                cfl["limit"] <= 0.0 or
                cfl["out_max"] >
                    cfl["limit"] * (1.0 + 64.0 * sys.float_info.epsilon)):
            raise SystemExit(f"{path}: step {step} committed CFL failed")
        tolerances = {}
        for metric in ("eos", "continuity", "energy", "closed_mass", "gauge"):
            residual = terminal.get(f"{metric}_residual")
            tolerance = terminal.get(f"{metric}_tolerance")
            if (isinstance(residual, bool) or isinstance(tolerance, bool) or
                    not isinstance(residual, (int, float)) or
                    not isinstance(tolerance, (int, float)) or
                    not math.isfinite(residual) or not math.isfinite(tolerance) or
                    residual < 0.0 or tolerance <= 0.0 or residual > tolerance):
                raise SystemExit(
                    f"{path}: step {step} terminal {metric}={residual}/{tolerance} is inadmissible")
            tolerances[metric] = tolerance
        if tolerances["energy"] != tolerances["continuity"]:
            raise SystemExit(f"{path}: step {step} energy gate diverges from continuity")
    if revisions[1] != revisions[0] + 1:
        raise SystemExit(f"{path}: final-flux revisions are not one publication per accepted step")
    all_rows.append(rows)

for index, wanted in enumerate(expected):
    control = all_rows[0][index]
    changed = all_rows[1][index]
    if control["time"] != changed["time"]:
        raise SystemExit(
            f"rank-count branches reached different physical times at step {wanted[0]}")
    control_revision = control["terminal_physical_audit"]["final_flux_revision"]
    changed_revision = changed["terminal_physical_audit"]["final_flux_revision"]
    if control_revision != changed_revision:
        raise SystemExit(
            f"rank-count branches published different final-flux revisions at step {wanted[0]}")
    control_advective = control["momentum_predictor_limiter"]["advective_cfl"]
    changed_advective = changed["momentum_predictor_limiter"]["advective_cfl"]
    for field in ("plan", "activity_collective"):
        if control_advective[field] != changed_advective[field]:
            raise SystemExit(
                f"rank-count branches published different advective CFL {field} "
                f"at step {wanted[0]}")
]=])
run_checked("rank-change coupled lifecycle certificate"
  "${PYTHON}" -c "${assert_evidence}" "${fixed_dt}"
  "${control_root}/evidence.jsonl"
  "${rank_change_root}/evidence.jsonl")

run_mpi("rank-neutral equivalent endpoint comparison MPI4" 4
  "${COMPARATOR}"
  "${case_root}" "${control_root}/Restart"
  "${case_root}" "${rank_change_root}/Restart"
  --json "${comparison_json}"
  --max-relative-rms "${maximum_relative_rms}"
  --max-relative-linf "${maximum_relative_linf}"
  --equivalent-endpoints)

set(assert_comparison [=[
import json
import math
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
fixed_dt = float(sys.argv[2])
maximum_relative_rms = float(sys.argv[3])
maximum_relative_linf = float(sys.argv[4])
report = json.loads(path.read_text(encoding="utf-8"))
if (report.get("schema") != "HUNDUN_V04_FULL_HALF_RESTART_COMPARE_V1" or
        report.get("comparison_mode") != "equivalent_endpoints" or
        report.get("passed") is not True or report.get("failures") != []):
    raise SystemExit(f"{path}: equivalent endpoint certificate did not pass")
thresholds = report.get("comparison_thresholds", {})
if (thresholds.get("maximum_relative_rms") != maximum_relative_rms or
        thresholds.get("maximum_relative_linf") != maximum_relative_linf):
    raise SystemExit(f"{path}: comparison thresholds were not recorded exactly")
authority = report.get("authority", {})
zero_fingerprint = "0x0000000000000000"
for name in ("source_stl", "surface", "topology", "interface_metric"):
    if authority.get(name) in (None, zero_fingerprint):
        raise SystemExit(f"{path}: IBM {name} authority is absent")

required_metrics = {"U", "pi", "p_abs", "h", "rho", "T", "q", "rhoU",
                    "phi", "phi_over_A"}
metrics = report.get("metrics")
if not isinstance(metrics, dict) or set(metrics) != required_metrics:
    raise SystemExit(f"{path}: incomplete physical endpoint metric catalog")
for name, metric in metrics.items():
    weight = metric.get("weight")
    values = (metric.get("rms"), metric.get("linf"),
              metric.get("relative_rms"), metric.get("relative_linf"))
    if (isinstance(weight, bool) or not isinstance(weight, (int, float)) or
            not math.isfinite(weight) or weight <= 0.0 or
            any(isinstance(value, bool) or not isinstance(value, (int, float)) or
                not math.isfinite(value) or value < 0.0 for value in values)):
        raise SystemExit(f"{path}: {name} comparison is non-finite")
    if values[2] > maximum_relative_rms or values[3] > maximum_relative_linf:
        raise SystemExit(f"{path}: {name} comparison exceeds its explicit threshold")

for endpoint_name in ("full", "half"):
    endpoint = report.get(endpoint_name)
    if not isinstance(endpoint, dict):
        raise SystemExit(f"{path}: missing {endpoint_name} endpoint")
    expected_time = 4.0 * fixed_dt
    time_tolerance = 128.0 * sys.float_info.epsilon * max(1.0, abs(expected_time))
    if (endpoint.get("step") != 4 or endpoint.get("dt") != fixed_dt or
            not isinstance(endpoint.get("time"), (int, float)) or
            abs(endpoint["time"] - expected_time) > time_tolerance):
        raise SystemExit(f"{path}: {endpoint_name} is not the common step-4 target")
    extrema = endpoint.get("extrema", {})
    for name in ("p_abs", "rho", "T"):
        bounds = extrema.get(name)
        if (not isinstance(bounds, list) or len(bounds) != 2 or
                any(isinstance(value, bool) or not isinstance(value, (int, float)) or
                    not math.isfinite(value) for value in bounds) or bounds[0] <= 0.0):
            raise SystemExit(f"{path}: {endpoint_name} has inadmissible {name} extrema")
    violations = endpoint.get("violations")
    if not isinstance(violations, dict) or any(value != 0 for value in violations.values()):
        raise SystemExit(f"{path}: {endpoint_name} has positivity/IBM/flux violations")
    boundaries = endpoint.get("boundaries")
    if not isinstance(boundaries, list) or len(boundaries) != 6:
        raise SystemExit(f"{path}: {endpoint_name} boundary certificate is incomplete")
    for boundary in boundaries:
        for field in ("inlet_reversal_faces", "velocity_inlet_flux_mismatch_faces",
                      "outlet_backflow_faces", "prohibited_backflow_faces",
                      "symmetry_nonpositive_zero", "symmetry_negative_zero"):
            if boundary.get(field) != 0:
                raise SystemExit(
                    f"{path}: {endpoint_name} {boundary.get('face')} fails {field}")
]=])
run_checked("rank-change endpoint physics certificate"
  "${PYTHON}" -c "${assert_comparison}" "${comparison_json}" "${fixed_dt}"
  "${maximum_relative_rms}" "${maximum_relative_linf}")
