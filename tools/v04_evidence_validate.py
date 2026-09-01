#!/usr/bin/env python3
"""Validate and advance immutable HUNDUN-FLOW v0.4 release evidence."""

import argparse
import hashlib
import json
import math
import os
import random
import re
import statistics
import struct
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Set


SCHEMA = "HUNDUN_V04_CANDIDATE_V1"
RUNTIME_SCHEMA = "HUNDUN_V04_EVIDENCE_V7"
V6_RUNTIME_SCHEMA = "HUNDUN_V04_EVIDENCE_V6"
V5_RUNTIME_SCHEMA = "HUNDUN_V04_EVIDENCE_V5"
V4_RUNTIME_SCHEMA = "HUNDUN_V04_EVIDENCE_V4"
V3_RUNTIME_SCHEMA = "HUNDUN_V04_EVIDENCE_V3"
V2_RUNTIME_SCHEMA = "HUNDUN_V04_EVIDENCE_V2"
LEGACY_RUNTIME_SCHEMA = "HUNDUN_V04_EVIDENCE_V1"
EQUIVALENCE_SCHEMA = "HUNDUN_V04_COAST_EQUIVALENCE_V2"
GATES = ("focused", "full2", "frozen", "full20", "literature", "final")
REQUIRED_CANDIDATE = (
    "head",
    "tree",
    "compiler",
    "linker",
    "flags",
    "binary_path",
    "binary_sha256",
    "tests_on_path",
    "tests_on_sha256",
    "tests_off_path",
    "tests_off_sha256",
    "inputs",
    "compiled_plan",
    "cpu_plan",
    "mpi_version",
    "mpi_thread_level",
    "ranks",
    "process_grid",
    "cpu",
    "numa",
    "affinity",
    "environment",
    "command",
    "output_state",
    "checkpoint_state",
    "start_timestamp",
    "performance_policy",
)

REQUIRED_GATE_CHECKS = {
    "focused": (
        "release_focused_passed",
        "asan_passed",
        "ubsan_passed",
        "force_mutation_red_green_passed",
        "positive_normal_donor_passed",
    ),
    "full2": (
        "equivalence_receipt_complete",
        "both_products_completed",
        "resource_contracts_passed",
        "scientific_work_matched",
        "directional_screen_passed",
    ),
    "frozen": (
        "clean_worktree",
        "candidate_files_revalidated",
        "candidate_identity_sealed",
    ),
    "full20": (
        "paired_policy_matched",
        "paired_statistics_accept",
        "scientific_work_matched",
    ),
    "literature": (
        "literature_authority_complete",
        "physical_accuracy_accept",
    ),
    "final": (
        "numerical_correctness_accept",
        "robustness_accept",
        "coast_short_performance_accept",
        "literature_physical_accuracy_accept",
        "provenance_accept",
    ),
}

FROZEN_STAGE_IDS = (10, 12, 15, 20, 30, 40, 45, 50, 60, 70)
UINT8_MAX = (1 << 8) - 1
UINT16_MAX = (1 << 16) - 1
UINT32_MAX = (1 << 32) - 1
UINT64_MAX = (1 << 64) - 1
INT32_MIN = -(1 << 31)
INT32_MAX = (1 << 31) - 1
MAX_FINITE_DOUBLE = float.fromhex("0x1.fffffffffffffp+1023")
PRESSURE_EXTENSION_FIELDS = (
    "reduction_calls", "operator_applies", "preconditioner_applies",
    "norm_breakdown_restarts", "recycle_offered_directions",
    "recycle_retained_directions", "recycle_operator_applies",
    "recycle_reduction_calls", "recycle_projection_attempted",
    "recycle_projection_accepted", "recycle_projected_true_residual",
    "recycle_cycle_corrections", "recycle_capture_vector_passes",
    "recycle_capture_cycle_attempts", "recycle_capture_reduction_calls",
    "recycle_capture_blocking_operations",
)
V3_REQUIRED_FIELDS = (
    "schema", "build", "binary", "case", "stl", "product", "cpu_plan", "step",
    "time", "requested_bdf_order", "bdf_order",
    "temporal_method_fallback", "thermophysical_predictor_calls",
    "launcher_ns", "max_rank_step_ns", "max_rank_rss_bytes",
    "max_node_rss_bytes", "structured_messages", "structured_bytes",
    "ibm_messages", "ibm_bytes", "blocking_collectives",
    "nonblocking_collectives", "reduction_ns", "linear_iterations",
    "exact_numeric_refills", "coarse_numeric_refills",
    "preconditioner_setups", "preconditioner_reuses", "heap_allocations",
    "pressure_solve_calls", "pressure_solve_contract",
    "terminal_physical_audit", "momentum_predictor_solve_calls",
    "momentum_predictor_limiter", "thermophysical_predictor",
    "thermophysical_enthalpy_endpoint", "pressure", "momentum_predictor",
    "startup", "retry", "restart_recovery", "statistics_eligible", "stages",
)
V3_TERMINAL_FIELDS = (
    "present", "final_flux_revision", "eos_residual", "eos_tolerance",
    "continuity_residual", "continuity_tolerance", "energy_residual",
    "energy_tolerance", "closed_mass_residual", "closed_mass_tolerance",
    "gauge_residual", "gauge_tolerance",
)
V3_MOMENTUM_LIMITER_FIELDS = ("limited", "theta", "activations")
V5_MOMENTUM_LIMITER_FIELDS = (
    "scheme", "limited", "retained_correction_l1_ratio", "limited_faces",
    "advective_cfl",
)
V5_COMMITTED_CONVECTIVE_CFL_FIELDS = ("out_max", "abs_max", "limit")
V5_ADVECTIVE_CONVECTIVE_CFL_FIELDS = (
    "present", "plan", "time_revision", "density_revision",
    "face_flux_revision", "activity_collective", "dt", "out_max",
    "abs_max", "limit",
)
V6_CANDIDATE_IDENTITY_FIELDS = (
    "schema", "head", "tree", "build_manifest_sha256",
    "executable_sha256", "identity_sha256",
)
V6_RUN_START_FIELDS = (
    "kind", "previous_step", "previous_time", "restart_manifest_sha256",
)
V6_MOMENTUM_LIMITER_FIELDS = (
    "scheme", "limited", "correction_metrics_applicability",
    "retained_correction_l1_ratio", "minimum_face_alpha",
    "active_correction_faces", "limited_faces", "limited_face_fraction",
    "advective_cfl",
)
V6_ADVECTIVE_CONVECTIVE_CFL_FIELDS = (
    "present", "plan", "time_revision_collective",
    "density_view_collective", "face_flux_view_collective",
    "activity_collective", "dt", "out_max", "abs_max", "limit",
)
V6_COMMITTED_CONVECTIVE_CFL_FIELDS = (
    "valid", "density_revision",
    "final_flux_revision", "density_view", "face_flux_view",
    "activity_collective", "dt", "out_max", "abs_max", "limit",
    "out_winner", "abs_winner",
)
V6_CFL_WINNER_FIELDS = (
    "valid", "global_cell", "rank", "out", "abs", "density_volume",
    "outgoing_mass_flow", "absolute_mass_flow",
)
V3_THERMOPHYSICAL_PREDICTOR_FIELDS = (
    "limited", "theta", "low_state", "mass_flux_scale", "constraint",
    "limiting_rank", "limiting_global_cell", "low_margin", "high_margin",
    "low_order_substeps", "low_order_transport_passes",
    "low_order_halo_exchanges", "blocking_collectives",
    "enthalpy_endpoint_alpha", "bdf_endpoint_alpha", "source_endpoint_alpha",
    "enthalpy_solve_calls",
)
V3_LINEAR_SOLVE_FIELDS = (
    "status_code", "termination", "iterations", "initial_true_residual",
    "final_true_residual", "recursive_residual", "reduction_calls",
    "operator_applies", "preconditioner_applies", "norm_breakdown_restarts",
)
V3_PRESSURE_SOLVE_FIELDS = (
    "corrector", "status_code", "termination", "iterations",
    "initial_true_residual", "final_true_residual", "recursive_residual",
    "final_convergence_metric", "convergence_limit", "convergence_audits",
    "convergence_rejections", "reduction_calls", "operator_applies",
    "preconditioner_applies", "norm_breakdown_restarts",
) + PRESSURE_EXTENSION_FIELDS
V3_STAGE_FIELDS = ("id", "min_ns", "mean_ns", "max_ns")
V4_REFINEMENT_FIELDS = (
    "pressure_energy_refinement_solve_calls",
    "pressure_energy_refinement_termination",
    "pressure_energy_refinement",
)
V4_REFINEMENT_SOLVE_FIELDS = (
    "ordinal", "target_generation", "collective_lineage", "status_code",
    "termination", "iterations", "initial_true_residual",
    "final_true_residual", "recursive_residual",
    "final_convergence_metric", "convergence_limit", "convergence_audits",
    "convergence_rejections", "reduction_calls", "operator_applies",
    "preconditioner_applies", "norm_breakdown_restarts",
) + PRESSURE_EXTENSION_FIELDS
PRESSURE_ENERGY_REFINEMENT_CAPACITY = 12
HISTORICAL_PRESSURE_ENERGY_REFINEMENT_CAPACITY = 6
LINEAR_TERMINATIONS = (
    "converged", "zero_rhs", "maximum_iterations", "breakdown",
    "non_finite", "operator_failure", "preconditioner_failure",
    "collective_failure", "invalid_plan", "convergence_audit_failure",
)
PREDICTOR_LOW_STATES = (
    "none", "bdf_accepted_rate", "scaled_euler", "implicit_upwind",
    "implicit_upwind_source_limited", "bdf_accepted_rate_homotopy",
    "bdf_local_donor_flux", "bdf_local_donor_flux_source_limited",
)
PREDICTOR_CONSTRAINTS = (
    "none", "density", "independent_species", "dependent_species",
    "enthalpy_lower", "enthalpy_upper",
)


class EvidenceError(ValueError):
    pass


def require_object_fields(value: Any, fields, name: str) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise EvidenceError(f"{name} must be an object")
    missing = [field for field in fields if field not in value]
    if missing:
        raise EvidenceError(f"{name} missing required field {missing[0]}")
    return value


def require_integer(value: Any, name: str, minimum: int = 0,
                    maximum: int = UINT64_MAX) -> int:
    if (not isinstance(value, int) or isinstance(value, bool) or
            value < minimum or value > maximum):
        raise EvidenceError(
            f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


def require_boolean(value: Any, name: str) -> bool:
    if not isinstance(value, bool):
        raise EvidenceError(f"{name} must be a boolean")
    return value


def require_finite_number(value: Any, name: str) -> Any:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise EvidenceError(f"{name} must be a finite number")
    if isinstance(value, float) and not math.isfinite(value):
        raise EvidenceError(f"{name} must be a finite number")
    if value < -MAX_FINITE_DOUBLE or value > MAX_FINITE_DOUBLE:
        raise EvidenceError(f"{name} must be a finite number")
    return value


def require_nonnegative_finite_number(value: Any, name: str) -> Any:
    value = require_finite_number(value, name)
    if value < 0:
        raise EvidenceError(f"{name} must be nonnegative")
    return value


def validate_v3_required_shape(record: Any, line_number: int) -> None:
    prefix = f"line {line_number}: V3 runtime"
    record = require_object_fields(record, V3_REQUIRED_FIELDS, prefix)
    require_object_fields(record["terminal_physical_audit"],
                          V3_TERMINAL_FIELDS,
                          f"{prefix}.terminal_physical_audit")
    require_object_fields(record["momentum_predictor_limiter"],
                          V3_MOMENTUM_LIMITER_FIELDS,
                          f"{prefix}.momentum_predictor_limiter")
    require_object_fields(record["thermophysical_predictor"],
                          V3_THERMOPHYSICAL_PREDICTOR_FIELDS,
                          f"{prefix}.thermophysical_predictor")
    require_object_fields(record["thermophysical_enthalpy_endpoint"],
                          V3_LINEAR_SOLVE_FIELDS,
                          f"{prefix}.thermophysical_enthalpy_endpoint")
    momentum = record["momentum_predictor"]
    if not isinstance(momentum, list) or len(momentum) != 3:
        raise EvidenceError(f"{prefix}.momentum_predictor must have 3 solves")
    for component, solve in enumerate(momentum):
        require_object_fields(
            solve, ("component",) + V3_LINEAR_SOLVE_FIELDS,
            f"{prefix}.momentum_predictor[{component}]")
    pressure = record["pressure"]
    if not isinstance(pressure, list) or len(pressure) != 2:
        raise EvidenceError(f"{prefix}.pressure must have 2 solves")
    for corrector, solve in enumerate(pressure):
        require_object_fields(
            solve, V3_PRESSURE_SOLVE_FIELDS,
            f"{prefix}.pressure[{corrector}]")
    stages = record["stages"]
    if not isinstance(stages, list) or not stages:
        raise EvidenceError(f"{prefix}.stages must be a nonempty array")
    for index, stage in enumerate(stages):
        require_object_fields(stage, V3_STAGE_FIELDS,
                              f"{prefix}.stages[{index}]")


def validate_v3_top_level_types(record: Dict[str, Any],
                                line_number: int) -> None:
    prefix = f"line {line_number}: V3 runtime"
    for field in ("build", "binary", "case", "product", "cpu_plan"):
        require_integer(record[field], f"{prefix}.{field}", 1)
    require_integer(record["stl"], f"{prefix}.stl")
    require_integer(record["step"], f"{prefix}.step", 1)
    require_finite_number(record["time"], f"{prefix}.time")
    for field in ("requested_bdf_order", "bdf_order",
                  "thermophysical_predictor_calls", "pressure_solve_calls",
                  "momentum_predictor_solve_calls"):
        require_integer(record[field], f"{prefix}.{field}", 0, UINT8_MAX)
    for field in (
            "launcher_ns", "max_rank_step_ns", "max_rank_rss_bytes",
            "max_node_rss_bytes", "structured_messages", "structured_bytes",
            "ibm_messages", "ibm_bytes", "blocking_collectives",
            "nonblocking_collectives", "reduction_ns", "linear_iterations",
            "exact_numeric_refills", "coarse_numeric_refills",
            "preconditioner_setups", "preconditioner_reuses",
            "heap_allocations"):
        require_integer(record[field], f"{prefix}.{field}")
    for field in ("temporal_method_fallback", "startup", "retry",
                  "restart_recovery", "statistics_eligible"):
        require_boolean(record[field], f"{prefix}.{field}")


def validate_v3_serialized_solve_types(solve: Dict[str, Any],
                                       prefix: str) -> None:
    require_integer(solve["status_code"], f"{prefix}.status_code", 0,
                    UINT16_MAX)
    termination = solve["termination"]
    if not isinstance(termination, str) or termination not in LINEAR_TERMINATIONS:
        raise EvidenceError(f"{prefix}.termination is invalid")
    require_integer(solve["iterations"], f"{prefix}.iterations", 0,
                    UINT32_MAX)
    for field in ("reduction_calls", "operator_applies",
                  "preconditioner_applies", "norm_breakdown_restarts"):
        require_integer(solve[field], f"{prefix}.{field}")
    for field in ("initial_true_residual", "final_true_residual",
                  "recursive_residual"):
        require_nonnegative_finite_number(solve[field], f"{prefix}.{field}")


def validate_v3_nested_types(record: Dict[str, Any], line_number: int) -> None:
    prefix = f"line {line_number}: V3 runtime"
    limiter = record["momentum_predictor_limiter"]
    require_boolean(limiter["limited"], f"{prefix}.momentum limiter.limited")
    require_finite_number(limiter["theta"],
                          f"{prefix}.momentum limiter.theta")
    require_integer(limiter["activations"],
                    f"{prefix}.momentum limiter.activations", 0, UINT32_MAX)

    predictor = record["thermophysical_predictor"]
    require_boolean(predictor["limited"], f"{prefix}.predictor.limited")
    for field in ("theta", "mass_flux_scale", "low_margin", "high_margin",
                  "enthalpy_endpoint_alpha", "bdf_endpoint_alpha",
                  "source_endpoint_alpha"):
        require_finite_number(predictor[field], f"{prefix}.predictor.{field}")
    for field in ("low_state", "constraint"):
        if not isinstance(predictor[field], str):
            raise EvidenceError(f"{prefix}.predictor.{field} must be a string")
    require_integer(predictor["limiting_rank"],
                    f"{prefix}.predictor.limiting_rank", INT32_MIN, INT32_MAX)
    limiting_cell = predictor["limiting_global_cell"]
    if not isinstance(limiting_cell, list) or len(limiting_cell) != 3:
        raise EvidenceError(
            f"{prefix}.predictor.limiting_global_cell must have 3 entries")
    for coordinate, value in enumerate(limiting_cell):
        require_integer(value,
                        f"{prefix}.predictor.limiting_global_cell[{coordinate}]",
                        INT32_MIN, INT32_MAX)
    require_integer(predictor["low_order_transport_passes"],
                    f"{prefix}.predictor.low_order_transport_passes")
    for field in ("low_order_substeps", "low_order_halo_exchanges",
                  "blocking_collectives"):
        require_integer(predictor[field], f"{prefix}.predictor.{field}",
                        0, UINT32_MAX)
    require_integer(predictor["enthalpy_solve_calls"],
                    f"{prefix}.predictor.enthalpy_solve_calls", 0, UINT8_MAX)

    validate_v3_serialized_solve_types(
        record["thermophysical_enthalpy_endpoint"],
        f"{prefix}.thermophysical_enthalpy_endpoint")
    for component, solve in enumerate(record["momentum_predictor"]):
        require_integer(solve["component"],
                        f"{prefix}.momentum_predictor[{component}].component",
                        0, UINT8_MAX)
        validate_v3_serialized_solve_types(
            solve, f"{prefix}.momentum_predictor[{component}]")

    for corrector, solve in enumerate(record["pressure"]):
        solve_prefix = f"{prefix}.pressure[{corrector}]"
        require_integer(solve["corrector"], f"{solve_prefix}.corrector",
                        0, UINT8_MAX)
        require_integer(solve["status_code"], f"{solve_prefix}.status_code",
                        0, UINT16_MAX)
        termination = solve["termination"]
        if (not isinstance(termination, str) or
                termination not in LINEAR_TERMINATIONS):
            raise EvidenceError(f"{solve_prefix}.termination is invalid")
        for field in (
                "convergence_audits", "convergence_rejections",
                "reduction_calls", "operator_applies",
                "preconditioner_applies", "norm_breakdown_restarts",
                "recycle_offered_directions", "recycle_retained_directions",
                "recycle_operator_applies", "recycle_reduction_calls",
                "recycle_cycle_corrections",
                "recycle_capture_vector_passes",
                "recycle_capture_cycle_attempts",
                "recycle_capture_reduction_calls",
                "recycle_capture_blocking_operations"):
            require_integer(solve[field], f"{solve_prefix}.{field}")
        require_integer(solve["iterations"], f"{solve_prefix}.iterations", 0,
                        UINT32_MAX)
        for field in ("recycle_projection_attempted",
                      "recycle_projection_accepted"):
            require_boolean(solve[field], f"{solve_prefix}.{field}")
        for field in (
                "initial_true_residual", "final_true_residual",
                "recursive_residual", "final_convergence_metric",
                "convergence_limit", "recycle_projected_true_residual"):
            require_nonnegative_finite_number(solve[field],
                                              f"{solve_prefix}.{field}")

    for index, stage in enumerate(record["stages"]):
        stage_prefix = f"{prefix}.stages[{index}]"
        require_integer(stage["id"], f"{stage_prefix}.id", 1, UINT16_MAX)
        for field in ("min_ns", "mean_ns", "max_ns"):
            require_integer(stage[field], f"{stage_prefix}.{field}")


def validate_v3_temporal_invariants(record: Dict[str, Any],
                                    line_number: int) -> None:
    prefix = f"line {line_number}: V3 runtime"
    requested = record["requested_bdf_order"]
    executed = record["bdf_order"]
    predictor_calls = record["thermophysical_predictor_calls"]
    fallback = record["temporal_method_fallback"]
    predictor_collectives = record["thermophysical_predictor"][
        "blocking_collectives"]
    valid_method = (
        requested in (1, 2) and executed in (1, 2) and
        ((fallback and requested == 2 and executed == 1 and
          predictor_calls == 2) or
         (not fallback and requested == executed and predictor_calls == 1)))
    if (not valid_method or
            record["blocking_collectives"] < predictor_collectives):
        raise EvidenceError(f"{prefix} has invalid temporal-method evidence")


def validate_v3_momentum_invariants(record: Dict[str, Any],
                                    line_number: int) -> None:
    prefix = f"line {line_number}: V3 runtime"
    limiter = record["momentum_predictor_limiter"]
    theta = limiter["theta"]
    limited = limiter["limited"]
    if (record["momentum_predictor_solve_calls"] != 3 or theta < 0 or
            theta > 1 or limiter["activations"] != (1 if limited else 0) or
            limited != (theta < 1)):
        raise EvidenceError(f"{prefix} has invalid momentum limiter evidence")
    for component, solve in enumerate(record["momentum_predictor"]):
        if (solve["component"] != component or solve["status_code"] != 0 or
                solve["termination"] not in ("converged", "zero_rhs") or
                solve["final_true_residual"] >
                solve["initial_true_residual"]):
            raise EvidenceError(
                f"{prefix}.momentum_predictor[{component}] was not accepted")


def validate_v3_predictor_invariants(record: Dict[str, Any],
                                     line_number: int) -> None:
    prefix = f"line {line_number}: V3 runtime"
    predictor = record["thermophysical_predictor"]
    endpoint = record["thermophysical_enthalpy_endpoint"]
    limited = predictor["limited"]
    theta = predictor["theta"]
    mass_scale = predictor["mass_flux_scale"]
    low_state = predictor["low_state"]
    constraint = predictor["constraint"]
    substeps = predictor["low_order_substeps"]
    halo_exchanges = predictor["low_order_halo_exchanges"]
    if (theta < 0 or theta > 1 or mass_scale < 0 or mass_scale > 1 or
            low_state not in PREDICTOR_LOW_STATES or
            constraint not in PREDICTOR_CONSTRAINTS):
        raise EvidenceError(f"{prefix} has invalid predictor values")

    if not limited:
        valid_fast_collectives = (
            predictor["blocking_collectives"] > 1
            if record["temporal_method_fallback"]
            else predictor["blocking_collectives"] == 1)
        valid_predictor_branch = (
            theta == 1 and mass_scale == 1 and
            predictor["low_margin"] == 0 and
            predictor["high_margin"] == 0 and substeps == 0 and
            predictor["low_order_halo_exchanges"] == 0 and
            predictor["low_order_transport_passes"] == 0 and
            predictor["limiting_rank"] == -1 and constraint == "none" and
            low_state == "none" and valid_fast_collectives)
    else:
        donor_transport = low_state in (
            "bdf_local_donor_flux",
            "bdf_local_donor_flux_source_limited")
        valid_halo_count = (
            halo_exchanges == 1 if donor_transport
            else halo_exchanges == substeps - 1)
        valid_predictor_branch = (
            theta < 1 and predictor["low_margin"] >= 0 and
            predictor["high_margin"] < 0 and substeps > 0 and
            valid_halo_count and predictor["blocking_collectives"] > 1 and
            predictor["low_order_transport_passes"] > 0 and
            all(coordinate >= 0 for coordinate in
                predictor["limiting_global_cell"]) and
            predictor["limiting_rank"] >= 0 and constraint != "none" and
            low_state != "none" and
            (low_state not in (
                "bdf_accepted_rate", "bdf_accepted_rate_homotopy",
                "implicit_upwind", "implicit_upwind_source_limited") or
             mass_scale == 1))
    if not valid_predictor_branch:
        raise EvidenceError(f"{prefix} has invalid predictor branch evidence")

    uses_enthalpy_endpoint = low_state in (
        "implicit_upwind", "implicit_upwind_source_limited")
    enthalpy_alpha = predictor["enthalpy_endpoint_alpha"]
    if uses_enthalpy_endpoint:
        valid_enthalpy_endpoint = (
            predictor["enthalpy_solve_calls"] == 1 and
            endpoint["status_code"] == 0 and
            endpoint["termination"] in ("converged", "zero_rhs") and
            endpoint["operator_applies"] > 0 and
            endpoint["reduction_calls"] > 0 and 0 <= enthalpy_alpha <= 1 and
            ((low_state == "implicit_upwind" and enthalpy_alpha == 1) or
             (low_state == "implicit_upwind_source_limited" and
              enthalpy_alpha < 1)))
    else:
        valid_enthalpy_endpoint = (
            predictor["enthalpy_solve_calls"] == 0 and enthalpy_alpha == 1)

    bdf_alpha = predictor["bdf_endpoint_alpha"]
    valid_bdf_endpoint = (
        0 <= bdf_alpha <= 1 and
        ((low_state == "bdf_accepted_rate_homotopy" and bdf_alpha < 1) or
         (low_state != "bdf_accepted_rate_homotopy" and bdf_alpha == 1)))
    source_alpha = predictor["source_endpoint_alpha"]
    valid_source_endpoint = (
        0 <= source_alpha <= 1 and
        ((low_state == "bdf_local_donor_flux_source_limited" and
          source_alpha < 1) or
         (low_state != "bdf_local_donor_flux_source_limited" and
          source_alpha == 1)))
    if (not valid_enthalpy_endpoint or not valid_bdf_endpoint or
            not valid_source_endpoint):
        raise EvidenceError(f"{prefix} has invalid predictor endpoint evidence")


def require_nonnegative_number(value: Any, name: str) -> Any:
    return require_nonnegative_finite_number(value, name)


def validate_v3_terminal_audit(value: Any, contract: str,
                               line_number: int) -> None:
    prefix = f"line {line_number}: terminal_physical_audit"
    if not isinstance(value, dict) or value.get("present") is not True:
        raise EvidenceError(f"{prefix} is missing")
    require_integer(value.get("final_flux_revision"),
                    f"{prefix}.final_flux_revision", 1, UINT64_MAX)
    tolerances = {}
    for metric in ("eos", "continuity", "energy", "closed_mass", "gauge"):
        residual = require_nonnegative_number(
            value.get(f"{metric}_residual"), f"{prefix}.{metric}_residual")
        tolerance = require_nonnegative_number(
            value.get(f"{metric}_tolerance"), f"{prefix}.{metric}_tolerance")
        tolerances[metric] = tolerance
        energy_evidence_only = (metric == "energy" and
                                contract == "pressure_continuity" and
                                tolerance == 0.0)
        if not energy_evidence_only:
            if tolerance == 0.0:
                raise EvidenceError(
                    f"{prefix}.{metric}_tolerance must be positive")
            if residual > tolerance:
                raise EvidenceError(f"{prefix}.{metric} exceeds tolerance")
    if (contract == "continuity_energy_coupled" and
            tolerances["energy"] != tolerances["continuity"]):
        raise EvidenceError(
            f"{prefix}.energy_tolerance must equal continuity_tolerance")


def validate_v3_runtime_record(record: Dict[str, Any],
                               line_number: int) -> str:
    validate_v3_required_shape(record, line_number)
    validate_v3_top_level_types(record, line_number)
    validate_v3_nested_types(record, line_number)
    validate_v3_temporal_invariants(record, line_number)
    validate_v3_momentum_invariants(record, line_number)
    validate_v3_predictor_invariants(record, line_number)
    contract = record["pressure_solve_contract"]
    if contract not in ("pressure_continuity", "continuity_energy_coupled"):
        raise EvidenceError(
            f"line {line_number}: invalid pressure solve contract")
    validate_v3_terminal_audit(record["terminal_physical_audit"], contract,
                               line_number)
    return contract


def reject_v4_fields_in_v3(record: Dict[str, Any], line_number: int) -> None:
    for field in V4_REFINEMENT_FIELDS:
        if field in record:
            raise EvidenceError(
                f"line {line_number}: V3 runtime carries V4 field {field}")


def validate_v4_refinement(record: Dict[str, Any], line_number: int,
                           capacity: int) -> None:
    prefix = f"line {line_number}: V4 runtime"
    require_object_fields(record, V4_REFINEMENT_FIELDS, prefix)
    count = require_integer(
        record["pressure_energy_refinement_solve_calls"],
        f"{prefix}.pressure_energy_refinement_solve_calls", 0,
        capacity)
    termination = record["pressure_energy_refinement_termination"]
    if termination != "component_residuals_converged":
        raise EvidenceError(
            f"{prefix}.pressure_energy_refinement_termination is not accepted")
    refinements = record["pressure_energy_refinement"]
    if not isinstance(refinements, list) or len(refinements) != count:
        raise EvidenceError(
            f"{prefix}.pressure_energy_refinement must match solve count")

    target_generation = None
    collective_lineages = set()  # type: Set[int]
    refinement_iterations = 0
    for index, solve in enumerate(refinements):
        solve_prefix = f"{prefix}.pressure_energy_refinement[{index}]"
        solve = require_object_fields(
            solve, V4_REFINEMENT_SOLVE_FIELDS, solve_prefix)
        ordinal = require_integer(solve["ordinal"], f"{solve_prefix}.ordinal",
                                  1, capacity)
        target = require_integer(solve["target_generation"],
                                 f"{solve_prefix}.target_generation", 1)
        lineage = require_integer(solve["collective_lineage"],
                                  f"{solve_prefix}.collective_lineage", 1)
        if ordinal != index + 1:
            raise EvidenceError(f"{solve_prefix}.ordinal is out of order")
        if target_generation is None:
            target_generation = target
        elif target != target_generation:
            raise EvidenceError(
                f"{solve_prefix}.target_generation changed within one target")
        if lineage in collective_lineages:
            raise EvidenceError(
                f"{solve_prefix}.collective_lineage is not unique")
        collective_lineages.add(lineage)

        require_integer(solve["status_code"], f"{solve_prefix}.status_code",
                        0, UINT16_MAX)
        if solve["status_code"] != 0 or solve["termination"] not in (
                "converged", "zero_rhs"):
            raise EvidenceError(f"{solve_prefix} was not accepted")
        iterations = require_integer(solve["iterations"],
                                     f"{solve_prefix}.iterations", 0,
                                     UINT32_MAX)
        refinement_iterations += iterations
        for field in (
                "reduction_calls", "operator_applies",
                "preconditioner_applies", "norm_breakdown_restarts",
                "convergence_audits", "convergence_rejections",
                "recycle_offered_directions", "recycle_retained_directions",
                "recycle_operator_applies", "recycle_reduction_calls",
                "recycle_cycle_corrections",
                "recycle_capture_vector_passes",
                "recycle_capture_cycle_attempts",
                "recycle_capture_reduction_calls",
                "recycle_capture_blocking_operations"):
            require_integer(solve[field], f"{solve_prefix}.{field}")
        for field in ("recycle_projection_attempted",
                      "recycle_projection_accepted"):
            require_boolean(solve[field], f"{solve_prefix}.{field}")
        for field in (
                "initial_true_residual", "final_true_residual",
                "recursive_residual", "final_convergence_metric",
                "convergence_limit", "recycle_projected_true_residual"):
            require_nonnegative_finite_number(solve[field],
                                              f"{solve_prefix}.{field}")

        capture_attempts = solve["recycle_capture_cycle_attempts"]
        valid_capture = (
            capture_attempts >= solve["recycle_cycle_corrections"] and
            capture_attempts <= UINT64_MAX // 2 and
            solve["recycle_capture_vector_passes"] == 2 * capture_attempts and
            solve["recycle_capture_reduction_calls"] == capture_attempts and
            solve["recycle_capture_blocking_operations"] ==
            2 * capture_attempts)
        no_audit_or_projection = (
            solve["convergence_audits"] == 0 and
            solve["convergence_rejections"] == 0 and
            solve["final_convergence_metric"] == 0 and
            solve["convergence_limit"] == 0 and
            solve["recycle_offered_directions"] == 0 and
            solve["recycle_retained_directions"] == 0 and
            solve["recycle_operator_applies"] == 0 and
            solve["recycle_reduction_calls"] == 0 and
            solve["recycle_projection_attempted"] is False and
            solve["recycle_projection_accepted"] is False and
            solve["recycle_projected_true_residual"] == 0)
        if not valid_capture or not no_audit_or_projection:
            raise EvidenceError(
                f"{solve_prefix} violates refinement capture provenance")

    if record["linear_iterations"] < refinement_iterations:
        raise EvidenceError(
            f"{prefix}.linear_iterations omits refinement solve work")


def validate_v5_runtime_record(record: Dict[str, Any],
                               line_number: int) -> str:
    """Validate V5 additions, then reuse the frozen V3/V4 base contracts."""
    prefix = f"line {line_number}: V5 runtime"
    record = require_object_fields(record, V3_REQUIRED_FIELDS, prefix)
    terminal = require_object_fields(
        record["terminal_physical_audit"],
        V3_TERMINAL_FIELDS + ("committed_convective_cfl",),
        f"{prefix}.terminal_physical_audit")
    cfl = require_object_fields(
        terminal["committed_convective_cfl"],
        V5_COMMITTED_CONVECTIVE_CFL_FIELDS,
        f"{prefix}.terminal_physical_audit.committed_convective_cfl")
    out_max = require_nonnegative_finite_number(
        cfl["out_max"], f"{prefix}.committed CFL.out_max")
    require_nonnegative_finite_number(
        cfl["abs_max"], f"{prefix}.committed CFL.abs_max")
    limit = require_nonnegative_finite_number(
        cfl["limit"], f"{prefix}.committed CFL.limit")
    if limit == 0.0:
        raise EvidenceError(f"{prefix}.committed CFL.limit must be positive")
    if out_max > limit * (1.0 + 64.0 * float.fromhex("0x1.0p-52")):
        raise EvidenceError(f"{prefix}.committed CFL exceeds its limit")

    limiter = require_object_fields(
        record["momentum_predictor_limiter"],
        V5_MOMENTUM_LIMITER_FIELDS,
        f"{prefix}.momentum_predictor_limiter")
    if limiter["scheme"] != "common_face_afc_v2":
        raise EvidenceError(f"{prefix}.momentum limiter scheme is invalid")
    limited = require_boolean(
        limiter["limited"], f"{prefix}.momentum limiter.limited")
    ratio = require_finite_number(
        limiter["retained_correction_l1_ratio"],
        f"{prefix}.momentum limiter.retained_correction_l1_ratio")
    faces = require_integer(
        limiter["limited_faces"], f"{prefix}.momentum limiter.limited_faces",
        0, UINT32_MAX)
    if (ratio <= 0.0 or ratio > 1.0 or
            (limited and (faces == 0 or ratio >= 1.0)) or
            (not limited and (faces != 0 or ratio != 1.0))):
        raise EvidenceError(f"{prefix} has invalid common-face AFC evidence")
    advective = require_object_fields(
        limiter["advective_cfl"], V5_ADVECTIVE_CONVECTIVE_CFL_FIELDS,
        f"{prefix}.momentum_predictor_limiter.advective_cfl")
    if not require_boolean(
            advective["present"], f"{prefix}.advective CFL.present"):
        raise EvidenceError(f"{prefix}.advective CFL must be present")
    for field in ("plan", "time_revision", "density_revision",
                  "face_flux_revision"):
        require_integer(advective[field], f"{prefix}.advective CFL.{field}",
                        1, UINT64_MAX)
    require_integer(advective["activity_collective"],
                    f"{prefix}.advective CFL.activity_collective",
                    0, UINT64_MAX)
    if ((record["stl"] == 0) !=
            (advective["activity_collective"] == 0)):
        raise EvidenceError(
            f"{prefix}.advective CFL IBM activity disagrees with stl")
    if advective["face_flux_revision"] == terminal["final_flux_revision"]:
        raise EvidenceError(
            f"{prefix}.advective and committed CFL revisions are identical")
    dt = require_nonnegative_finite_number(
        advective["dt"], f"{prefix}.advective CFL.dt")
    advective_out = require_nonnegative_finite_number(
        advective["out_max"], f"{prefix}.advective CFL.out_max")
    require_nonnegative_finite_number(
        advective["abs_max"], f"{prefix}.advective CFL.abs_max")
    advective_limit = require_nonnegative_finite_number(
        advective["limit"], f"{prefix}.advective CFL.limit")
    if dt == 0.0 or advective_limit == 0.0:
        raise EvidenceError(
            f"{prefix}.advective CFL dt and limit must be positive")
    if advective_limit != terminal["committed_convective_cfl"]["limit"]:
        raise EvidenceError(
            f"{prefix}.advective and committed CFL limits differ")
    if advective_out > advective_limit * (
            1.0 + 64.0 * float.fromhex("0x1.0p-52")):
        raise EvidenceError(f"{prefix}.advective CFL exceeds its limit")
    if "theta" in limiter or "activations" in limiter:
        raise EvidenceError(f"{prefix} carries legacy V4 limiter fields")

    # The pressure, thermophysical predictor, solver, timing and terminal
    # metric contracts are unchanged.  Project only the renamed limiter for
    # the frozen V3 validator; validate the actual face count above.
    legacy_projection = dict(record)
    legacy_projection["terminal_physical_audit"] = dict(terminal)
    legacy_projection["terminal_physical_audit"].pop(
        "committed_convective_cfl")
    legacy_projection["momentum_predictor_limiter"] = {
        "limited": limited,
        "theta": ratio,
        "activations": 1 if limited else 0,
    }
    return validate_v3_runtime_record(legacy_projection, line_number)


def _v6_close(left: float, right: float) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= 128.0 * float.fromhex("0x1.0p-52") * scale


def validate_v6_cfl_winner(value: Any, dt: float, maximum: float,
                           maximum_field: str, name: str) -> None:
    winner = require_object_fields(value, V6_CFL_WINNER_FIELDS, name)
    if not require_boolean(winner["valid"], f"{name}.valid"):
        raise EvidenceError(f"{name} must be valid")
    cell = winner["global_cell"]
    if (not isinstance(cell, list) or len(cell) != 3 or
            any(not isinstance(coordinate, int) or isinstance(coordinate, bool)
                or coordinate < 0 or coordinate > INT32_MAX
                for coordinate in cell)):
        raise EvidenceError(f"{name}.global_cell is invalid")
    require_integer(winner["rank"], f"{name}.rank", 0, INT32_MAX)
    out = require_nonnegative_finite_number(winner["out"], f"{name}.out")
    absolute = require_nonnegative_finite_number(
        winner["abs"], f"{name}.abs")
    density_volume = require_nonnegative_finite_number(
        winner["density_volume"], f"{name}.density_volume")
    outgoing = require_nonnegative_finite_number(
        winner["outgoing_mass_flow"], f"{name}.outgoing_mass_flow")
    absolute_flow = require_nonnegative_finite_number(
        winner["absolute_mass_flow"], f"{name}.absolute_mass_flow")
    if density_volume == 0.0:
        raise EvidenceError(f"{name}.density_volume must be positive")
    predicted_out = dt * outgoing / density_volume
    predicted_absolute = dt * absolute_flow / (2.0 * density_volume)
    if (not _v6_close(out, predicted_out) or
            not _v6_close(absolute, predicted_absolute)):
        raise EvidenceError(f"{name} violates the cell-CFL formula")
    observed = out if maximum_field == "out" else absolute
    if not _v6_close(observed, maximum):
        raise EvidenceError(f"{name} does not bind {maximum_field}_max")


def validate_v6_v7_candidate_identity(record: Dict[str, Any],
                                      line_number: int) -> None:
    prefix = f"line {line_number}: V6/V7 candidate_identity"
    identity = require_object_fields(
        record.get("candidate_identity"), V6_CANDIDATE_IDENTITY_FIELDS,
        prefix)
    runtime_schema = record.get("schema")
    if runtime_schema == RUNTIME_SCHEMA:
        identity_schema = "HUNDUN_V04_RUNTIME_CANDIDATE_IDENTITY_V2"
    elif runtime_schema == V6_RUNTIME_SCHEMA:
        identity_schema = "HUNDUN_V04_RUNTIME_CANDIDATE_IDENTITY_V1"
    else:
        raise EvidenceError(f"{prefix} has invalid runtime schema")
    if identity["schema"] != identity_schema:
        raise EvidenceError(f"{prefix}.schema is invalid")
    require_git_object(identity["head"], f"{prefix}.head")
    require_git_object(identity["tree"], f"{prefix}.tree")
    for field in ("build_manifest_sha256", "executable_sha256",
                  "identity_sha256"):
        require_sha256(identity[field], f"{prefix}.{field}")
        if re.fullmatch(r"[0-9a-f]{64}", identity[field]) is None:
            raise EvidenceError(f"{prefix}.{field} is not lowercase hex")
    payload = (
        f"schema={identity_schema}\n" +
        (f"evidence_schema={RUNTIME_SCHEMA}\n"
         if runtime_schema == RUNTIME_SCHEMA else "") +
        f"head={identity['head']}\n"
        f"tree={identity['tree']}\n"
        f"build_manifest_sha256={identity['build_manifest_sha256']}\n"
        f"executable_sha256={identity['executable_sha256']}\n")
    if hashlib.sha256(payload.encode("utf-8")).hexdigest() != \
            identity["identity_sha256"]:
        raise EvidenceError(f"{prefix}.identity_sha256 is invalid")
    if record.get("build") != int(identity["build_manifest_sha256"][:16], 16):
        raise EvidenceError(f"{prefix} disagrees with build")
    if record.get("binary") != int(identity["executable_sha256"][:16], 16):
        raise EvidenceError(f"{prefix} disagrees with binary")


def validate_v6_run_start(record: Dict[str, Any],
                          line_number: int) -> Dict[str, Any]:
    prefix = f"line {line_number}: V6 run_start"
    anchor = require_object_fields(record.get("run_start"),
                                   V6_RUN_START_FIELDS, prefix)
    previous_step = require_integer(anchor["previous_step"],
                                    f"{prefix}.previous_step", 0)
    previous_time = require_nonnegative_finite_number(
        anchor["previous_time"], f"{prefix}.previous_time")
    step = require_integer(record.get("step"),
                           f"line {line_number}: V6 step", 1)
    if step <= previous_step:
        raise EvidenceError(f"{prefix} does not precede the current step")
    first = step == previous_step + 1
    kind = anchor["kind"]
    if kind == "fresh":
        if (previous_step != 0 or previous_time != 0.0 or
                anchor["restart_manifest_sha256"] is not None or
                record.get("restart_recovery") is not False or
                record.get("startup") is not (step == 1)):
            raise EvidenceError(f"{prefix} has an invalid fresh anchor")
    elif kind == "restart":
        require_sha256(anchor["restart_manifest_sha256"],
                       f"{prefix}.restart_manifest_sha256")
        if (re.fullmatch(r"[0-9a-f]{64}",
                         anchor["restart_manifest_sha256"]) is None or
                int(anchor["restart_manifest_sha256"], 16) == 0 or
                previous_step == 0 or record.get("startup") is not False or
                (record.get("restart_recovery") is True and not first)):
            raise EvidenceError(f"{prefix} has an invalid restart anchor")
        if (first and record.get("restart_recovery") is True and
                (record.get("requested_bdf_order") != 1 or
                 record.get("bdf_order") != 1)):
            raise EvidenceError(
                f"{prefix} lifecycle disagrees with legacy Restart recovery")
    else:
        raise EvidenceError(f"{prefix}.kind is invalid")
    if first and not _v6_close(record["previous_committed_time"],
                               previous_time):
        raise EvidenceError(
            f"{prefix} does not bind the first previous committed time")
    return anchor


def load_v04_restart_manifest(path: Path) -> Dict[str, Any]:
    """Read the frozen V0.4 manifest authority needed by a restart run."""
    data = path.read_bytes()
    header_format = "<8sIIiiiQQQdddQQI"
    header_size = struct.calcsize(header_format)
    if len(data) < header_size + 8:
        raise EvidenceError(f"restart manifest {path} is truncated")
    fields = struct.unpack_from(header_format, data)
    (magic, version, rank_count, cells_x, cells_y, cells_z,
     plan, schema, geometry, time, dt, pressure_reference,
     step, controller_state, field_count) = fields
    cursor = header_size + 4 * field_count
    exact_history = version == 2
    valid_exact = True
    if exact_history and len(data) >= cursor + struct.calcsize("<ddQQI") + 8:
        (previous_pressure_reference, closed_mass_target,
         final_mass_flux_revision, previous_mass_flux_revision,
         rate_count) = struct.unpack_from("<ddQQI", data, cursor)
        cursor += struct.calcsize("<ddQQI")
        valid_exact = (
            controller_state != 0 and
            math.isfinite(previous_pressure_reference) and
            previous_pressure_reference > 0.0 and
            math.isfinite(closed_mass_target) and closed_mass_target > 0.0 and
            previous_mass_flux_revision > 0 and
            previous_mass_flux_revision < final_mass_flux_revision and
            final_mass_flux_revision < UINT64_MAX and
            1 <= rate_count <= 64)
        cursor += 4 * rate_count
    elif exact_history:
        valid_exact = False
    expected_size = cursor + 40 * rank_count + 8
    if (magic != b"H4MANI01" or version not in (1, 2) or
            rank_count == 0 or
            min(cells_x, cells_y, cells_z) <= 0 or
            min(plan, schema, geometry) == 0 or
            not math.isfinite(time) or time < 0.0 or
            not math.isfinite(dt) or dt <= 0.0 or
            not math.isfinite(pressure_reference) or step == 0 or
            field_count < 3 or field_count > 64 or
            not valid_exact or
            len(data) != expected_size):
        raise EvidenceError(f"restart manifest {path} has an invalid header")
    offset = 1469598103934665603
    prime = 1099511628211
    checksum = offset
    for value in data[:-8]:
        checksum = ((checksum ^ value) * prime) & UINT64_MAX
    stored_checksum = struct.unpack_from("<Q", data, len(data) - 8)[0]
    if checksum != stored_checksum:
        raise EvidenceError(f"restart manifest {path} failed integrity")
    return {
        "sha256": hashlib.sha256(data).hexdigest(),
        "step": step,
        "time": time,
        "backward_euler_recovery": version == 1,
    }


def validate_v6_v7_runtime_record(record: Dict[str, Any], line_number: int,
                                  refinement_capacity: int) -> str:
    """Validate immutable identity, first-step time, AFC and CFL V6/V7 data."""
    prefix = f"line {line_number}: V6/V7 runtime"
    for field in ("candidate_identity", "run_start",
                  "previous_committed_time"):
        if field not in record:
            raise EvidenceError(f"{prefix} missing required field {field}")
    validate_v6_v7_candidate_identity(record, line_number)
    validate_v6_run_start(record, line_number)

    terminal = require_object_fields(
        record.get("terminal_physical_audit"),
        V3_TERMINAL_FIELDS + ("committed_convective_cfl",),
        f"{prefix}.terminal_physical_audit")
    cfl = require_object_fields(
        terminal["committed_convective_cfl"],
        V6_COMMITTED_CONVECTIVE_CFL_FIELDS,
        f"{prefix}.committed_convective_cfl")
    if not require_boolean(cfl["valid"], f"{prefix}.committed CFL.valid"):
        raise EvidenceError(f"{prefix}.committed CFL must be valid")
    for field in ("density_revision",
                  "final_flux_revision"):
        require_integer(cfl[field], f"{prefix}.committed CFL.{field}", 1)
    if cfl["final_flux_revision"] != terminal["final_flux_revision"]:
        raise EvidenceError(f"{prefix}.committed CFL final flux drifted")
    density_view = require_object_fields(
        cfl["density_view"], ("field", "collective"),
        f"{prefix}.committed CFL.density_view")
    require_integer(density_view["field"],
                    f"{prefix}.committed CFL.density_view.field", 0,
                    UINT32_MAX)
    require_integer(density_view["collective"],
                    f"{prefix}.committed CFL.density_view.collective", 1)
    flux_view = require_object_fields(
        cfl["face_flux_view"], ("collective",),
        f"{prefix}.committed CFL.face_flux_view")
    require_integer(flux_view["collective"],
                    f"{prefix}.committed CFL.face_flux_view.collective", 1)
    activity = require_integer(
        cfl["activity_collective"],
        f"{prefix}.committed CFL.activity_collective", 0)
    if ((record["stl"] == 0) != (activity == 0)):
        raise EvidenceError(f"{prefix}.committed CFL IBM authority drifted")
    dt = require_nonnegative_finite_number(cfl["dt"],
                                           f"{prefix}.committed CFL.dt")
    out_max = require_nonnegative_finite_number(
        cfl["out_max"], f"{prefix}.committed CFL.out_max")
    abs_max = require_nonnegative_finite_number(
        cfl["abs_max"], f"{prefix}.committed CFL.abs_max")
    limit = require_nonnegative_finite_number(
        cfl["limit"], f"{prefix}.committed CFL.limit")
    if dt == 0.0 or limit == 0.0:
        raise EvidenceError(f"{prefix}.committed CFL dt/limit must be positive")
    if out_max > limit * (1.0 + 64.0 * float.fromhex("0x1.0p-52")):
        raise EvidenceError(f"{prefix}.committed CFL exceeds its limit")
    validate_v6_cfl_winner(cfl["out_winner"], dt, out_max, "out",
                           f"{prefix}.committed CFL.out_winner")
    validate_v6_cfl_winner(cfl["abs_winner"], dt, abs_max, "abs",
                           f"{prefix}.committed CFL.abs_winner")

    limiter = require_object_fields(
        record.get("momentum_predictor_limiter"),
        V6_MOMENTUM_LIMITER_FIELDS,
        f"{prefix}.momentum_predictor_limiter")
    if limiter["scheme"] != "common_face_afc_v3_owner":
        raise EvidenceError(f"{prefix}.momentum limiter scheme is invalid")
    limited = require_boolean(limiter["limited"],
                              f"{prefix}.momentum limiter.limited")
    applicability = limiter["correction_metrics_applicability"]
    active = require_integer(
        limiter["active_correction_faces"],
        f"{prefix}.momentum limiter.active_correction_faces", 0, UINT32_MAX)
    faces = require_integer(
        limiter["limited_faces"],
        f"{prefix}.momentum limiter.limited_faces", 0, UINT32_MAX)
    if applicability == "not_applicable":
        if (active != 0 or faces != 0 or limited or
                limiter["retained_correction_l1_ratio"] is not None or
                limiter["minimum_face_alpha"] is not None or
                limiter["limited_face_fraction"] is not None):
            raise EvidenceError(f"{prefix} fabricates inactive AFC metrics")
        ratio = 1.0
    elif applicability == "applicable":
        ratio = require_finite_number(
            limiter["retained_correction_l1_ratio"],
            f"{prefix}.momentum limiter.retained_correction_l1_ratio")
        minimum_alpha = require_finite_number(
            limiter["minimum_face_alpha"],
            f"{prefix}.momentum limiter.minimum_face_alpha")
        fraction = require_finite_number(
            limiter["limited_face_fraction"],
            f"{prefix}.momentum limiter.limited_face_fraction")
        expected_fraction = faces / active if active else math.nan
        if (active == 0 or faces > active or ratio <= 0.0 or ratio > 1.0 or
                minimum_alpha < 0.0 or minimum_alpha > 1.0 or
                not _v6_close(fraction, expected_fraction) or
                limited != (faces > 0) or
                (limited and (ratio >= 1.0 or minimum_alpha >= 1.0)) or
                (not limited and (ratio != 1.0 or minimum_alpha != 1.0))):
            raise EvidenceError(f"{prefix} has invalid owner-face AFC metrics")
    else:
        raise EvidenceError(f"{prefix}.AFC applicability is invalid")

    previous = require_nonnegative_finite_number(
        record["previous_committed_time"],
        f"{prefix}.previous_committed_time")
    current = require_nonnegative_finite_number(record.get("time"),
                                                f"{prefix}.time")
    if current <= previous or not _v6_close(current - previous, dt):
        raise EvidenceError(
            f"{prefix} first/previous committed time disagrees with CFL dt")

    advective = require_object_fields(
        limiter["advective_cfl"], V6_ADVECTIVE_CONVECTIVE_CFL_FIELDS,
        f"{prefix}.momentum_predictor_limiter.advective_cfl")
    if not require_boolean(advective["present"],
                           f"{prefix}.advective CFL.present"):
        raise EvidenceError(f"{prefix}.advective CFL must be present")
    for field in ("plan", "time_revision_collective",
                  "density_view_collective", "face_flux_view_collective"):
        require_integer(advective[field], f"{prefix}.advective CFL.{field}",
                        1, UINT64_MAX)
    provisional_activity = require_integer(
        advective["activity_collective"],
        f"{prefix}.advective CFL.activity_collective", 0, UINT64_MAX)
    if ((record["stl"] == 0) != (provisional_activity == 0)):
        raise EvidenceError(
            f"{prefix}.advective CFL IBM activity disagrees with stl")
    if advective["face_flux_view_collective"] == \
            flux_view["collective"]:
        raise EvidenceError(
            f"{prefix}.provisional and committed CFL views are identical")
    provisional_dt = require_nonnegative_finite_number(
        advective["dt"], f"{prefix}.advective CFL.dt")
    provisional_out = require_nonnegative_finite_number(
        advective["out_max"], f"{prefix}.advective CFL.out_max")
    require_nonnegative_finite_number(
        advective["abs_max"], f"{prefix}.advective CFL.abs_max")
    provisional_limit = require_nonnegative_finite_number(
        advective["limit"], f"{prefix}.advective CFL.limit")
    if (provisional_dt == 0.0 or provisional_limit == 0.0 or
            provisional_out > provisional_limit *
            (1.0 + 64.0 * float.fromhex("0x1.0p-52"))):
        raise EvidenceError(f"{prefix}.advective CFL failed")

    # Reuse the frozen V5/V4/V3 contracts through a projection.  This does
    # not regenerate a historical oracle; it supplies only the legacy names
    # whose semantics V6 intentionally retains.
    projection = json.loads(json.dumps(record))
    projection["schema"] = V5_RUNTIME_SCHEMA
    projection.pop("candidate_identity")
    projection.pop("run_start")
    projection.pop("previous_committed_time")
    projection["terminal_physical_audit"]["committed_convective_cfl"] = {
        "out_max": out_max, "abs_max": abs_max, "limit": limit,
    }
    legacy_advective = {
        "present": advective["present"],
        "plan": advective["plan"],
        "time_revision": advective["time_revision_collective"],
        "density_revision": advective["density_view_collective"],
        "face_flux_revision": advective["face_flux_view_collective"],
        "activity_collective": advective["activity_collective"],
        "dt": advective["dt"],
        "out_max": advective["out_max"],
        "abs_max": advective["abs_max"],
        "limit": advective["limit"],
    }
    projection["momentum_predictor_limiter"] = {
        "scheme": "common_face_afc_v2",
        "limited": limited,
        "retained_correction_l1_ratio": ratio,
        "limited_faces": faces,
        "advective_cfl": legacy_advective,
    }
    contract = validate_v5_runtime_record(projection, line_number)
    validate_v4_refinement(projection, line_number, refinement_capacity)
    if (not _v6_close(advective["dt"], dt) or
            not _v6_close(advective["limit"], limit) or
            advective["activity_collective"] != activity):
        raise EvidenceError(f"{prefix}.provisional/committed CFL drifted")
    return contract


def reject_v6_fields_in_pre_v6(record: Dict[str, Any],
                               line_number: int) -> None:
    schema = record.get("schema")
    for field in ("candidate_identity", "run_start",
                  "previous_committed_time"):
        if field in record:
            raise EvidenceError(
                f"line {line_number}: {schema} runtime carries V6 {field}")
    limiter = record.get("momentum_predictor_limiter")
    if isinstance(limiter, dict):
        for field in ("correction_metrics_applicability",
                      "minimum_face_alpha", "active_correction_faces",
                      "limited_face_fraction"):
            if field in limiter:
                raise EvidenceError(
                    f"line {line_number}: {schema} runtime carries V6 "
                    f"limiter field {field}")


def reject_v5_fields_in_pre_v5(record: Dict[str, Any],
                               line_number: int) -> None:
    schema = record.get("schema")
    terminal = record.get("terminal_physical_audit")
    if isinstance(terminal, dict) and "committed_convective_cfl" in terminal:
        raise EvidenceError(
            f"line {line_number}: {schema} runtime carries V5 committed CFL")
    limiter = record.get("momentum_predictor_limiter")
    if isinstance(limiter, dict):
        for field in ("scheme", "retained_correction_l1_ratio",
                      "limited_faces", "advective_cfl"):
            if field in limiter:
                raise EvidenceError(
                    f"line {line_number}: {schema} runtime carries V5 limiter "
                    f"field {field}")


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def digest_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def require_sha256(value: Any, name: str) -> None:
    if not isinstance(value, str) or len(value) != 64:
        raise EvidenceError(f"{name} is not a SHA-256 hex digest")
    try:
        int(value, 16)
    except ValueError as error:
        raise EvidenceError(f"{name} is not a SHA-256 hex digest") from error


def require_git_object(value: Any, name: str) -> None:
    if (not isinstance(value, str) or len(value) not in (40, 64) or
            re.fullmatch(r"[0-9a-f]+", value) is None):
        raise EvidenceError(f"{name} is not a Git object id")


def require_text(value: Any, name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise EvidenceError(f"{name} must be nonempty text")


def require_file_hash(path_value: Any, digest_value: Any, name: str) -> None:
    require_text(path_value, f"{name}_path")
    require_sha256(digest_value, f"{name}_sha256")
    path = Path(path_value)
    if not path.is_absolute() or not path.is_file():
        raise EvidenceError(f"{name}_path must name an existing absolute file")
    if digest_file(path) != digest_value:
        raise EvidenceError(f"{name}_sha256 does not match {path}")


def validate_performance_policy(policy: Any) -> str:
    if (not isinstance(policy, dict) or
            policy.get("schema") != "HUNDUN_V04_PERFORMANCE_POLICY_V1"):
        raise EvidenceError("performance_policy has the wrong schema")
    if policy.get("grid") != [480, 480, 48] or policy.get("ranks") != 64:
        raise EvidenceError("performance policy must freeze 480x480x48 on 64 ranks")
    if policy.get("full2_steps") != [1, 2]:
        raise EvidenceError("performance policy must freeze full2 steps [1,2]")
    warmup = policy.get("full20_warmup_steps")
    measured = policy.get("full20_measured_steps")
    if (not isinstance(warmup, list) or not warmup or
            not isinstance(measured, list) or not measured or
            any(not isinstance(value, int) or value <= 0
                for value in warmup + measured) or
            len(set(warmup + measured)) != len(warmup) + len(measured) or
            sorted(warmup + measured) != list(range(1, 21))):
        raise EvidenceError("performance policy must partition steps 1..20 into warmup/measured")
    minimum_pairs = policy.get("minimum_pairs")
    maximum_pairs = policy.get("maximum_pairs")
    if (not isinstance(minimum_pairs, int) or minimum_pairs < 5 or
            not isinstance(maximum_pairs, int) or maximum_pairs < minimum_pairs):
        raise EvidenceError("performance policy requires 5 <= minimum_pairs <= maximum_pairs")
    if policy.get("first_pair_order") not in ("HC", "CH"):
        raise EvidenceError("performance policy first_pair_order must be HC or CH")
    if policy.get("timing_authority") not in (
            "max_rank_hot_step_sum", "launcher_hot_region"):
        raise EvidenceError("unsupported performance timing authority")
    if policy.get("serialized_output") is not False:
        raise EvidenceError("performance policy must disable serialized output")
    directional = policy.get("full2_max_directional_ratio")
    if (not isinstance(directional, (int, float)) or
            not math.isfinite(float(directional)) or directional <= 0.0):
        raise EvidenceError("full2 directional ratio must be finite and positive")
    return digest_bytes(canonical(policy))


def null_paths(value: Any, prefix: str = "") -> List[str]:
    paths = []  # type: List[str]
    if value is None:
        paths.append(prefix or "<root>")
    elif isinstance(value, dict):
        for key, child in value.items():
            child_prefix = f"{prefix}.{key}" if prefix else key
            paths.extend(null_paths(child, child_prefix))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            paths.extend(null_paths(child, f"{prefix}[{index}]"))
    return paths


def validate_named_sha256(value: Any, prefix: str = "") -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            child_prefix = f"{prefix}.{key}" if prefix else key
            if key.endswith("sha256") and child is not None:
                require_sha256(child, child_prefix)
            else:
                validate_named_sha256(child, child_prefix)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            validate_named_sha256(child, f"{prefix}[{index}]")


def validate_equivalence(path: Path, require_sealed: bool) -> str:
    receipt = load_json(path)
    if (not isinstance(receipt, dict) or
            receipt.get("schema") != EQUIVALENCE_SCHEMA or
            receipt.get("rule_status") != "FROZEN"):
        raise EvidenceError("COAST equivalence rule/schema is not frozen V2")
    status = receipt.get("receipt_status")
    if status not in ("UNSEALED", "SEALED"):
        raise EvidenceError("COAST equivalence receipt status is invalid")
    if require_sealed and status != "SEALED":
        raise EvidenceError("COAST equivalence receipt is not SEALED")

    decision = receipt.get("decision")
    algorithms = receipt.get("native_algorithms")
    common = receipt.get("common_problem")
    linear = receipt.get("linear_work")
    terminal = receipt.get("terminal_work")
    timing = receipt.get("timing")
    output = receipt.get("output_state")
    inventory = receipt.get("inventory")
    if any(not isinstance(value, dict) for value in
           (decision, inventory, algorithms, common, linear, terminal, timing,
            output)):
        raise EvidenceError("COAST equivalence receipt is missing a rule section")
    if (receipt.get("comparison_unit") !=
            "one_complete_accepted_physical_time_step" or
            decision.get("equal_pressure_coupling_iteration_counts_required")
            is not False or
            decision.get("native_algorithm_work_is_timed") is not True or
            decision.get("failed_or_retried_steps_are_statistics_eligible")
            is not False):
        raise EvidenceError("COAST scientific-work comparison decision drifted")

    rule_path_value = decision.get("rule_document")
    if not isinstance(rule_path_value, str) or not rule_path_value:
        raise EvidenceError("COAST equivalence rule document is missing")
    repository = path.resolve().parents[2]
    rule_path = repository / rule_path_value
    require_file_hash(str(rule_path), decision.get("rule_document_sha256"),
                      "equivalence_rule_document")

    inventory_path_value = inventory.get("document")
    if not isinstance(inventory_path_value, str) or not inventory_path_value:
        raise EvidenceError("COAST capability inventory document is missing")
    inventory_path = repository / inventory_path_value
    require_file_hash(str(inventory_path), inventory.get("document_sha256"),
                      "equivalence_capability_inventory")
    require_git_object(inventory.get("coast_head"), "inventory.coast_head")
    require_git_object(inventory.get("coast_index_tree"),
                       "inventory.coast_index_tree")
    require_sha256(inventory.get("coast_status_sha256"),
                   "inventory.coast_status_sha256")
    if status == "UNSEALED":
        if (inventory.get("status") != "REVIEWED_FAIL_CLOSED" or
                inventory.get("formal_pairing_authorized") is not False):
            raise EvidenceError("UNSEALED COAST inventory must remain fail closed")
    elif (inventory.get("status") != "VERIFIED_COMPLETE" or
          inventory.get("formal_pairing_authorized") is not True):
        raise EvidenceError("SEALED COAST inventory is not verified complete")

    common_receipt_meta = receipt.get("common_input_method_receipt")
    if not isinstance(common_receipt_meta, dict):
        raise EvidenceError("COAST common input/method receipt is missing")
    common_receipt_path_value = common_receipt_meta.get("document")
    if (not isinstance(common_receipt_path_value, str) or
            not common_receipt_path_value):
        raise EvidenceError("COAST common input/method receipt path is missing")
    common_receipt_path = repository / common_receipt_path_value
    require_file_hash(str(common_receipt_path),
                      common_receipt_meta.get("document_sha256"),
                      "equivalence_common_input_method_receipt")
    common_receipt = load_json(common_receipt_path)
    if (not isinstance(common_receipt, dict) or
            common_receipt.get("schema") !=
            "HUNDUN_V04_RE3900_COMMON_INPUT_METHOD_V1" or
            common_receipt.get("status") != "COMPLETE_PRE_IDENTITY" or
            common_receipt_meta.get("status") != "COMPLETE_PRE_IDENTITY"):
        raise EvidenceError("COAST common input/method receipt is incomplete")
    receipt_file_hash = digest_file(common_receipt_path)
    if terminal.get("audit_inputs_sha256") != receipt_file_hash:
        raise EvidenceError("common terminal audit input receipt hash drifted")
    common_terminal = common_receipt.get("common_terminal_evaluator")
    if (not isinstance(common_terminal, dict) or
            terminal.get("audit_tool_sha256") !=
            common_terminal.get("source_sha256")):
        raise EvidenceError("common terminal evaluator hash drifted")

    canonical_subhashes = (
        (common_receipt.get("boundaries"),
         common.get("boundaries", {}).get("resolved_values_sha256"),
         "common boundary values"),
        (common_receipt.get("turbulence_and_wall"),
         common.get("turbulence_and_wall", {}).get(
             "constants_and_clipping_sha256"),
         "common turbulence/wall constants"),
        (common_receipt.get("spatial_method_receipts", {}).get("hundun"),
         common.get("spatial_discretization", {}).get(
             "hundun_method_receipt", {}).get("canonical_object_sha256"),
         "HUNDUN method receipt"),
        (common_receipt.get("spatial_method_receipts", {}).get("coast"),
         common.get("spatial_discretization", {}).get(
             "coast_method_receipt", {}).get("canonical_object_sha256"),
         "COAST method receipt"))
    for value, expected, label in canonical_subhashes:
        require_sha256(expected, label + " sha256")
        if not isinstance(value, dict) or digest_bytes(canonical(value)) != expected:
            raise EvidenceError(label + " hash drifted")

    mesh = common.get("mesh")
    reference = common.get("reference_state")
    equation_scope = common.get("equation_scope")
    time = common.get("time")
    boundaries = common.get("boundaries")
    if any(not isinstance(value, dict) for value in
           (mesh, reference, equation_scope, time, boundaries)):
        raise EvidenceError("COAST common-problem section is incomplete")
    if (mesh.get("domain") != [[-5.0, -10.0, 0.0],
                                [15.0, 10.0, math.pi]] or
            mesh.get("cells") != [480, 480, 48] or
            mesh.get("diameter") != 1.0 or
            reference.get("reynolds_number") != 3900 or
            time.get("time_step") != 0.006 or
            time.get("startup_scheme") != "backward_euler" or
            time.get("subsequent_scheme") != "bdf2"):
        raise EvidenceError("COAST common mesh/reference/time problem drifted")
    required_scope = (
        "transient_continuity", "three_momentum_components",
        "sensible_enthalpy", "inert_single_species_eos_closure",
        "single_phase_low_mach_compressible", "local_absolute_pressure_eos",
        "pressure_storage_drho_dp_at_fixed_hY", "full_bdf_density_defect")
    if any(equation_scope.get(key) is not True for key in required_scope):
        raise EvidenceError("COAST common equation scope was weakened")
    if equation_scope.get("reacting") is not False:
        raise EvidenceError("COAST performance baseline must remain nonreacting")
    if boundaries != {
            "x_min": "velocity_inlet",
            "x_max": "absolute_pressure_outlet",
            "y_min": "symmetry", "y_max": "symmetry",
            "z_min": "symmetry", "z_max": "symmetry",
            "resolved_values_sha256": boundaries.get("resolved_values_sha256")
            }:
        raise EvidenceError("COAST common boundary kinds drifted")

    hundun_algorithm = algorithms.get("hundun")
    matching = algorithms.get("matching_rule")
    if (not isinstance(hundun_algorithm, dict) or
            not isinstance(matching, dict) or
            hundun_algorithm.get("pressure_velocity_coupling") != "PISO" or
            hundun_algorithm.get("pressure_corrections_per_attempt") != 2 or
            matching.get(
                "pressure_correction_or_outer_iteration_counts_are_not_matched")
            is not True or
            matching.get("all_registered_native_iterations_are_inside_the_step_timer")
            is not True or
            matching.get("both_products_must_pass_the_common_terminal_audit")
            is not True):
        raise EvidenceError("COAST/HUNDUN native-algorithm rule drifted")

    true_residual = linear.get("canonical_true_residual")
    pressure = linear.get("pressure")
    if (not isinstance(true_residual, dict) or
            not isinstance(pressure, dict) or
            true_residual.get("arithmetic") != "FP64" or
            true_residual.get("definition") !=
            "global_l2_norm_of_b_minus_Ax_on_active_algebraic_rows" or
            true_residual.get("termination_limit") !=
            "max(absolute_tolerance,relative_tolerance*global_l2_norm_of_b)" or
            true_residual.get("recursive_or_preconditioned_residual_alone_may_terminate")
            is not False or
            true_residual.get("true_residual_required_at_termination") is not True or
            pressure.get("absolute_tolerance") != 1e-8 or
            pressure.get("relative_tolerance") != 1e-6 or
            pressure.get("maximum_iterations_per_linear_solve") != 500):
        raise EvidenceError("COAST canonical true-residual rule drifted")

    tolerances = terminal.get("tolerances")
    if (terminal.get("audit_authority") !=
            "independent_common_post_step_evaluator_over_candidate_fields_before_commit" or
            terminal.get("all_four_gates_required") is not True or
            not isinstance(tolerances, dict) or
            any(tolerances.get(key) != 1e-6
                for key in ("eos", "continuity", "closed_mass", "gauge"))):
        raise EvidenceError("COAST common terminal audit rule drifted")
    if (timing.get("authority") != "max_rank_hot_step_sum" or
            timing.get("rank_reduction") != "MPI_MAX_after_local_timer_stop" or
            output != {"serialized_restart": False,
                       "serialized_visit": False,
                       "serialized_screen": False}):
        raise EvidenceError("COAST timer/output rule drifted")

    validate_named_sha256(receipt)
    if require_sealed:
        unresolved = null_paths(receipt)
        if unresolved:
            preview = ", ".join(unresolved[:8])
            if len(unresolved) > 8:
                preview += f", ... ({len(unresolved)} total)"
            raise EvidenceError("SEALED COAST receipt has unresolved fields: " +
                                preview)
        coast_algorithm = algorithms.get("coast")
        identities = receipt.get("immutable_identity")
        if (not isinstance(coast_algorithm, dict) or
                not isinstance(coast_algorithm.get(
                    "outer_iterations_per_attempt"), int) or
                coast_algorithm["outer_iterations_per_attempt"] <= 0):
            raise EvidenceError("SEALED COAST receipt needs a positive SIMPLE niter")
        require_text(coast_algorithm.get("outer_iteration_selection_authority"),
                     "COAST SIMPLE niter authority")
        require_text(timing.get("coast_max_rank_timer_authority"),
                     "COAST max-rank timer authority")
        if not isinstance(identities, dict):
            raise EvidenceError("SEALED COAST receipt lacks immutable identity")
        for key in ("hundun_head", "hundun_tree", "coast_head", "coast_tree"):
            require_git_object(identities.get(key), key)
    return digest_bytes(canonical(receipt))


def validate_candidate(candidate: Any) -> str:
    if not isinstance(candidate, dict):
        raise EvidenceError("candidate must be an object")
    missing = [key for key in REQUIRED_CANDIDATE if key not in candidate]
    if missing:
        raise EvidenceError("candidate missing: " + ", ".join(missing))
    require_git_object(candidate["head"], "head")
    require_git_object(candidate["tree"], "tree")
    for key in ("compiler", "linker", "flags", "mpi_version",
                "mpi_thread_level", "start_timestamp"):
        require_text(candidate[key], key)
    if re.match(r"^\d{4}-\d{2}-\d{2}T", candidate["start_timestamp"]) is None:
        raise EvidenceError("start_timestamp must use an ISO-8601 date-time")
    require_file_hash(candidate["binary_path"], candidate["binary_sha256"],
                      "binary")
    require_file_hash(candidate["tests_on_path"], candidate["tests_on_sha256"],
                      "tests_on")
    require_file_hash(candidate["tests_off_path"], candidate["tests_off_sha256"],
                      "tests_off")
    if not isinstance(candidate["inputs"], list) or not candidate["inputs"]:
        raise EvidenceError("candidate inputs must be a nonempty list")
    seen = set()  # type: Set[str]
    for item in candidate["inputs"]:
        if not isinstance(item, dict) or not isinstance(item.get("path"), str):
            raise EvidenceError("every input requires a path and sha256")
        require_sha256(item.get("sha256"), f"input {item.get('path')}")
        path = Path(item["path"])
        if not path.is_absolute() or not path.is_file():
            raise EvidenceError("candidate input must be an existing absolute file")
        if digest_file(path) != item["sha256"]:
            raise EvidenceError(f"candidate input hash mismatch: {path}")
        canonical_path = str(path.resolve())
        if canonical_path in seen:
            raise EvidenceError(f"duplicate input path: {item['path']}")
        seen.add(canonical_path)
    if not isinstance(candidate["ranks"], int) or candidate["ranks"] <= 0:
        raise EvidenceError("ranks must be positive")
    if (not isinstance(candidate["process_grid"], list) or
            len(candidate["process_grid"]) != 3 or
            any(not isinstance(value, int) or value <= 0
                for value in candidate["process_grid"]) or
            (candidate["process_grid"][0] * candidate["process_grid"][1] *
             candidate["process_grid"][2]) != candidate["ranks"]):
        raise EvidenceError("process_grid product must equal ranks")
    for key in ("cpu", "numa", "affinity", "output_state", "checkpoint_state"):
        if not isinstance(candidate[key], dict) or not candidate[key]:
            raise EvidenceError(f"{key} must be an explicit nonempty object")
    if (not isinstance(candidate["environment"], dict) or
            any(not isinstance(key, str) or not isinstance(value, str)
                for key, value in candidate["environment"].items())):
        raise EvidenceError("environment must be an explicit allowlist object")
    if (not isinstance(candidate["command"], list) or not candidate["command"] or
            any(not isinstance(value, str) or not value
                for value in candidate["command"])):
        raise EvidenceError("command must be a nonempty argv array")
    for key in ("compiled_plan", "cpu_plan"):
        if not isinstance(candidate[key], int) or candidate[key] <= 0:
            raise EvidenceError(f"{key} must be a positive fingerprint")
    validate_performance_policy(candidate["performance_policy"])
    return digest_bytes(canonical(candidate))


def validate_gate_receipt(receipt: Any, gate_name: str,
                          candidate_sha256: str) -> None:
    if (not isinstance(receipt, dict) or
            receipt.get("schema") != "HUNDUN_V04_GATE_RECEIPT_V1" or
            receipt.get("gate") != gate_name or
            receipt.get("accepted") is not True or
            receipt.get("candidate_sha256") != candidate_sha256):
        raise EvidenceError("receipt schema/gate/candidate/acceptance mismatch")
    checks = receipt.get("checks")
    if not isinstance(checks, dict):
        raise EvidenceError("receipt checks must be an object")
    for check in REQUIRED_GATE_CHECKS[gate_name]:
        if checks.get(check) is not True:
            raise EvidenceError(f"gate {gate_name} requires check {check}=true")
    if any(value is not True for value in checks.values()):
        raise EvidenceError("an accepted receipt cannot contain a failed check")
    require_text(receipt.get("start_timestamp"), "receipt start_timestamp")
    require_text(receipt.get("end_timestamp"), "receipt end_timestamp")


def validate_manifest(document: Any) -> None:
    if not isinstance(document, dict) or document.get("schema") != SCHEMA:
        raise EvidenceError(f"manifest schema must be {SCHEMA}")
    candidate_sha = validate_candidate(document.get("candidate"))
    if document.get("candidate_sha256") != candidate_sha:
        raise EvidenceError("candidate_sha256 does not match candidate")
    gates = document.get("gates")
    if not isinstance(gates, list) or len(gates) > len(GATES):
        raise EvidenceError("gates must be an ordered list")
    evidence_paths = set()  # type: Set[str]
    for index, gate in enumerate(gates):
        if not isinstance(gate, dict) or gate.get("name") != GATES[index]:
            raise EvidenceError("gate order is not focused->full2->frozen->full20->literature->final")
        if gate.get("candidate_sha256") != candidate_sha or gate.get("accepted") is not True:
            raise EvidenceError(f"gate {GATES[index]} has wrong candidate or is not accepted")
        require_sha256(gate.get("receipt_sha256"), "receipt_sha256")
        receipt_path = gate.get("receipt_path")
        if not isinstance(receipt_path, str):
            raise EvidenceError("gate receipt_path must be absolute")
        receipt = Path(receipt_path)
        if (not receipt.is_absolute() or not receipt.is_file() or
                digest_file(receipt) != gate["receipt_sha256"] or
                str(receipt.resolve()) in evidence_paths):
            raise EvidenceError("gate receipt is missing, changed or reused")
        evidence_paths.add(str(receipt.resolve()))
        for item in gate.get("evidence", []):
            path = item.get("path") if isinstance(item, dict) else None
            artifact = Path(path) if isinstance(path, str) else None
            if (artifact is None or not artifact.is_absolute() or
                    not artifact.is_file() or
                    str(artifact.resolve()) in evidence_paths):
                raise EvidenceError("evidence paths must be nonempty and immutable")
            require_sha256(item.get("sha256"), f"evidence {path}")
            if digest_file(artifact) != item["sha256"]:
                raise EvidenceError(f"evidence artifact changed: {artifact}")
            evidence_paths.add(str(artifact.resolve()))


def write_new(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(path, flags, 0o644)
    try:
        payload = json.dumps(document, indent=2, sort_keys=True,
                             ensure_ascii=False).encode("utf-8") + b"\n"
        cursor = 0
        while cursor < len(payload):
            cursor += os.write(descriptor, payload[cursor:])
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def replace_atomic(path: Path, document: Any) -> None:
    pending = path.with_name(path.name + f".pending-{os.getpid()}")
    write_new(pending, document)
    os.replace(pending, path)
    descriptor = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def initialize(candidate_path: Path, output: Path) -> None:
    candidate = load_json(candidate_path)
    candidate_sha = validate_candidate(candidate)
    document = {"schema": SCHEMA, "candidate": candidate,
                "candidate_sha256": candidate_sha, "gates": []}
    validate_manifest(document)
    write_new(output, document)


def advance(manifest_path: Path, gate_name: str, receipt_path: Path) -> None:
    document = load_json(manifest_path)
    validate_manifest(document)
    next_index = len(document["gates"])
    if next_index >= len(GATES) or gate_name != GATES[next_index]:
        expected = GATES[next_index] if next_index < len(GATES) else "none"
        raise EvidenceError(f"next gate is {expected}, not {gate_name}")
    receipt = load_json(receipt_path)
    validate_gate_receipt(receipt, gate_name, document["candidate_sha256"])
    evidence = receipt.get("evidence")
    if not isinstance(evidence, list) or not evidence:
        raise EvidenceError("receipt must bind at least one evidence artifact")
    already = {old["receipt_path"] for old in document["gates"]}
    already.update(item["path"] for old in document["gates"]
                   for item in old.get("evidence", []))
    receipt_absolute = str(receipt_path.resolve())
    if receipt_absolute in already:
        raise EvidenceError(f"receipt path already registered: {receipt_absolute}")
    checked = []
    for item in evidence:
        if not isinstance(item, dict) or not isinstance(item.get("path"), str):
            raise EvidenceError("invalid receipt evidence entry")
        path = Path(item["path"]).resolve()
        if str(path) in already or str(path) == receipt_absolute:
            raise EvidenceError(f"evidence path already registered: {path}")
        if not path.is_file():
            raise EvidenceError(f"missing evidence artifact: {path}")
        actual = digest_file(path)
        if item.get("sha256") != actual:
            raise EvidenceError(f"evidence hash mismatch: {path}")
        checked.append({"path": str(path), "sha256": actual})
    gate = {"name": gate_name, "accepted": True,
            "candidate_sha256": document["candidate_sha256"],
            "receipt_path": receipt_absolute,
            "receipt_sha256": digest_file(receipt_path),
            "evidence": checked}
    candidate = json.loads(json.dumps(document))
    candidate["gates"].append(gate)
    validate_manifest(candidate)
    replace_atomic(manifest_path, candidate)


def validate_runtime(path: Path, run_start_manifest: Path = None) -> None:
    prior_step = 0
    run_schema = None
    prior_modern_record = None
    prior_v5_record = None
    restart_authority = (load_v04_restart_manifest(run_start_manifest)
                         if run_start_manifest is not None else None)
    restart_authority_used = False
    seen = 0
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            record = json.loads(line)
            seen += 1
            if not isinstance(record, dict):
                raise EvidenceError(
                    f"line {line_number}: runtime record must be an object")
            schema = record.get("schema")
            if schema not in (RUNTIME_SCHEMA, V6_RUNTIME_SCHEMA,
                              V5_RUNTIME_SCHEMA,
                              V4_RUNTIME_SCHEMA,
                              V3_RUNTIME_SCHEMA, V2_RUNTIME_SCHEMA,
                              LEGACY_RUNTIME_SCHEMA):
                raise EvidenceError(f"line {line_number}: wrong runtime schema")
            if run_schema is None:
                run_schema = schema
            elif schema != run_schema:
                raise EvidenceError(
                    f"line {line_number}: runtime schema changed within one file")
            is_v7 = schema == RUNTIME_SCHEMA
            is_v6 = schema == V6_RUNTIME_SCHEMA
            is_modern = is_v7 or is_v6
            is_v5 = schema == V5_RUNTIME_SCHEMA
            is_v4 = schema == V4_RUNTIME_SCHEMA
            is_v3 = schema == V3_RUNTIME_SCHEMA
            is_v2 = schema == V2_RUNTIME_SCHEMA
            if not is_modern:
                reject_v6_fields_in_pre_v6(record, line_number)
            if not (is_modern or is_v5):
                reject_v5_fields_in_pre_v5(record, line_number)
            requires_pressure = is_modern or is_v5 or is_v4 or is_v3 or is_v2
            if is_modern:
                pressure_contract = validate_v6_v7_runtime_record(
                    record, line_number,
                    PRESSURE_ENERGY_REFINEMENT_CAPACITY if is_v7 else
                    HISTORICAL_PRESSURE_ENERGY_REFINEMENT_CAPACITY)
            elif is_v5:
                pressure_contract = validate_v5_runtime_record(
                    record, line_number)
            elif is_v4 or is_v3:
                pressure_contract = validate_v3_runtime_record(
                    record, line_number)
            else:
                pressure_contract = None
            if is_modern:
                validate_v4_refinement(
                    record, line_number,
                    PRESSURE_ENERGY_REFINEMENT_CAPACITY if is_v7 else
                    HISTORICAL_PRESSURE_ENERGY_REFINEMENT_CAPACITY)
            elif is_v5:
                validate_v4_refinement(
                    record, line_number,
                    HISTORICAL_PRESSURE_ENERGY_REFINEMENT_CAPACITY)
            elif is_v4:
                validate_v4_refinement(
                    record, line_number,
                    HISTORICAL_PRESSURE_ENERGY_REFINEMENT_CAPACITY)
            elif is_v3:
                reject_v4_fields_in_v3(record, line_number)
            for identity in ("build", "binary", "case", "product", "cpu_plan"):
                if not isinstance(record.get(identity), int) or record[identity] <= 0:
                    raise EvidenceError(f"line {line_number}: zero {identity}")
            step = record.get("step")
            if not isinstance(step, int) or step <= prior_step:
                raise EvidenceError(f"line {line_number}: steps not strictly increasing")
            prior_step = step
            if is_modern:
                if prior_modern_record is None:
                    anchor = record["run_start"]
                    if step != anchor["previous_step"] + 1:
                        raise EvidenceError(
                            f"line {line_number}: V6/V7 first row is not the "
                            "successor of its run-start anchor")
                    if not _v6_close(record["previous_committed_time"],
                                     anchor["previous_time"]):
                        raise EvidenceError(
                            f"line {line_number}: V6/V7 first row disagrees "
                            "with its run-start time")
                    if anchor["kind"] == "restart":
                        if restart_authority is None:
                            raise EvidenceError(
                                f"line {line_number}: V6/V7 restart evidence "
                                "requires --run-start-manifest")
                        if (anchor["restart_manifest_sha256"] !=
                                restart_authority["sha256"] or
                                anchor["previous_step"] !=
                                restart_authority["step"] or
                                not _v6_close(anchor["previous_time"],
                                              restart_authority["time"]) or
                                record["restart_recovery"] is not
                                restart_authority[
                                    "backward_euler_recovery"]):
                            raise EvidenceError(
                                f"line {line_number}: V6/V7 run-start anchor "
                                "disagrees with the frozen restart manifest")
                        restart_authority_used = True
                    elif restart_authority is not None:
                        raise EvidenceError(
                            f"line {line_number}: fresh V6/V7 evidence was "
                            "given a restart manifest")
                else:
                    for field in ("build", "binary", "candidate_identity",
                                  "run_start", "case", "stl", "product",
                                  "cpu_plan"):
                        if record[field] != prior_modern_record[field]:
                            raise EvidenceError(
                                f"line {line_number}: V6/V7 immutable run "
                                f"identity field {field} changed")
                    if step != prior_modern_record["step"] + 1:
                        raise EvidenceError(
                            f"line {line_number}: V6/V7 steps are not contiguous")
                    if not _v6_close(record["previous_committed_time"],
                                     prior_modern_record["time"]):
                        raise EvidenceError(
                            f"line {line_number}: V6/V7 previous committed "
                            "time does not bind the adjacent row")
                    current_advective = record[
                        "momentum_predictor_limiter"]["advective_cfl"]
                    prior_advective = prior_modern_record[
                        "momentum_predictor_limiter"]["advective_cfl"]
                    for field in ("plan", "activity_collective", "limit"):
                        if current_advective[field] != prior_advective[field]:
                            raise EvidenceError(
                                f"line {line_number}: V6/V7 static advective "
                                f"authority field {field} changed")
                prior_modern_record = record
                prior_v5_record = None
            elif is_v5:
                if prior_v5_record is not None:
                    for field in ("build", "binary", "case", "stl",
                                  "product", "cpu_plan"):
                        if record[field] != prior_v5_record[field]:
                            raise EvidenceError(
                                f"line {line_number}: V5 run identity "
                                f"field {field} changed")
                    if step != prior_v5_record["step"] + 1:
                        raise EvidenceError(
                            f"line {line_number}: V5 steps are not contiguous")
                    current_advective = record[
                        "momentum_predictor_limiter"]["advective_cfl"]
                    prior_advective = prior_v5_record[
                        "momentum_predictor_limiter"]["advective_cfl"]
                    for field in ("plan", "activity_collective", "limit"):
                        if current_advective[field] != prior_advective[field]:
                            raise EvidenceError(
                                f"line {line_number}: V5 static advective "
                                f"authority field {field} changed")
                    observed_dt = record["time"] - prior_v5_record["time"]
                    certified_dt = record["momentum_predictor_limiter"][
                        "advective_cfl"]["dt"]
                    scale = max(1.0, abs(observed_dt), abs(certified_dt))
                    if (not math.isfinite(observed_dt) or
                            observed_dt <= 0.0 or
                            abs(observed_dt - certified_dt) >
                            128.0 * float.fromhex("0x1.0p-52") * scale):
                        raise EvidenceError(
                            f"line {line_number}: V5 advective CFL dt "
                            "disagrees with the adjacent committed time delta")
                prior_v5_record = record
                prior_modern_record = None
            else:
                prior_modern_record = None
                prior_v5_record = None
            if ((record.get("startup") or record.get("retry") or
                 record.get("restart_recovery")) and
                    record.get("statistics_eligible")):
                raise EvidenceError(f"line {line_number}: transient record marked eligible")
            stages = record.get("stages")
            if not isinstance(stages, list) or not stages:
                raise EvidenceError(f"line {line_number}: missing stage timings")
            ids = set()
            for stage in stages:
                values = (stage.get("min_ns"), stage.get("mean_ns"),
                          stage.get("max_ns"))
                if (stage.get("id") in ids or any(not isinstance(value, int)
                                                   for value in values) or
                        not values[0] <= values[1] <= values[2]):
                    raise EvidenceError(f"line {line_number}: invalid stage timing")
                ids.add(stage["id"])
            if tuple(sorted(ids)) != FROZEN_STAGE_IDS:
                raise EvidenceError(f"line {line_number}: stage set is not frozen")
            for counter in (
                    "launcher_ns", "max_rank_step_ns", "max_rank_rss_bytes",
                    "max_node_rss_bytes", "structured_messages",
                    "structured_bytes", "ibm_messages", "ibm_bytes",
                    "blocking_collectives", "nonblocking_collectives",
                    "reduction_ns", "linear_iterations", "exact_numeric_refills",
                    "coarse_numeric_refills", "preconditioner_setups",
                    "preconditioner_reuses", "heap_allocations"):
                value = record.get(counter)
                if not isinstance(value, int) or value < 0:
                    raise EvidenceError(
                        f"line {line_number}: invalid nonnegative counter {counter}")
            pressure = record.get("pressure")
            pressure_calls = record.get("pressure_solve_calls")
            if requires_pressure and pressure is None and \
                    pressure_calls is None:
                raise EvidenceError(
                    f"line {line_number}: {schema} lacks pressure evidence")
            if pressure is not None or pressure_calls is not None:
                if pressure_calls != 2 or not isinstance(pressure, list) or \
                        len(pressure) != 2:
                    raise EvidenceError(
                        f"line {line_number}: invalid two-PISO pressure evidence")
                extension_enabled = False
                if is_v2 or is_v3 or is_v4 or is_v5 or is_modern:
                    extension_fields = set(PRESSURE_EXTENSION_FIELDS)
                    extension_presence = []
                    for solve in pressure:
                        if not isinstance(solve, dict):
                            raise EvidenceError(
                                f"line {line_number}: invalid pressure solve")
                        present = extension_fields.intersection(solve)
                        if present and present != extension_fields:
                            raise EvidenceError(
                                f"line {line_number}: partial pressure extension")
                        extension_presence.append(bool(present))
                    if extension_presence[0] != extension_presence[1]:
                        raise EvidenceError(
                            f"line {line_number}: inconsistent pressure extension presence")
                    extension_enabled = extension_presence[0]
                    if (is_v3 or is_v4 or is_v5 or is_modern) and not extension_enabled:
                        raise EvidenceError(
                            f"line {line_number}: V3 lacks pressure extension")
                for corrector, solve in enumerate(pressure, 1):
                    if not isinstance(solve, dict) or \
                            not isinstance(solve.get("corrector"), int) or \
                            isinstance(solve.get("corrector"), bool) or \
                            solve["corrector"] != corrector or \
                            not isinstance(solve.get("status_code"), int) or \
                            isinstance(solve.get("status_code"), bool) or \
                            solve["status_code"] != 0 or \
                            solve.get("termination") not in (
                                "converged", "zero_rhs"):
                        raise EvidenceError(
                            f"line {line_number}: invalid pressure solve {corrector}")
                    counters = ["iterations", "convergence_audits",
                                "convergence_rejections"]
                    if extension_enabled:
                        counters.extend((
                            "reduction_calls", "operator_applies",
                            "preconditioner_applies",
                            "norm_breakdown_restarts",
                            "recycle_offered_directions",
                            "recycle_retained_directions",
                            "recycle_operator_applies",
                            "recycle_reduction_calls",
                            "recycle_cycle_corrections",
                            "recycle_capture_vector_passes",
                            "recycle_capture_cycle_attempts",
                            "recycle_capture_reduction_calls",
                            "recycle_capture_blocking_operations"))
                    for counter in counters:
                        value = solve.get(counter)
                        if requires_pressure and \
                                (not isinstance(value, int) or
                                 isinstance(value, bool) or value < 0):
                            raise EvidenceError(
                                f"line {line_number}: invalid pressure {counter}")
                    if extension_enabled:
                        for flag in ("recycle_projection_attempted",
                                     "recycle_projection_accepted"):
                            if not isinstance(solve.get(flag), bool):
                                raise EvidenceError(
                                    f"line {line_number}: invalid pressure {flag}")
                        projected = solve.get("recycle_projected_true_residual")
                        if (not isinstance(projected, (int, float)) or
                                isinstance(projected, bool) or
                                not math.isfinite(projected) or projected < 0):
                            raise EvidenceError(
                                f"line {line_number}: invalid pressure projected residual")
                        offered = solve["recycle_offered_directions"]
                        retained = solve["recycle_retained_directions"]
                        if (offered > 4 or retained > offered or
                                solve["recycle_capture_cycle_attempts"] <
                                solve["recycle_cycle_corrections"] or
                                solve["recycle_capture_reduction_calls"] !=
                                solve["recycle_capture_cycle_attempts"] or
                                solve["recycle_capture_vector_passes"] !=
                                2 * solve["recycle_capture_cycle_attempts"] or
                                solve["recycle_capture_blocking_operations"] !=
                                2 * solve["recycle_capture_cycle_attempts"]):
                            raise EvidenceError(
                                f"line {line_number}: invalid recycle count relation")
                        attempted = solve["recycle_projection_attempted"]
                        accepted = solve["recycle_projection_accepted"]
                        if corrector == 1:
                            if (offered != 0 or retained != 0 or
                                    solve["recycle_operator_applies"] != 0 or
                                    solve["recycle_reduction_calls"] != 0 or
                                    attempted or accepted or projected != 0):
                                raise EvidenceError(
                                    f"line {line_number}: PISO1 has projection work")
                        elif (solve["recycle_cycle_corrections"] != 0 or
                              solve["recycle_capture_cycle_attempts"] != 0 or
                              solve["recycle_capture_vector_passes"] != 0 or
                              solve["recycle_capture_reduction_calls"] != 0 or
                              solve["recycle_capture_blocking_operations"] != 0):
                            raise EvidenceError(
                                f"line {line_number}: PISO2 has capture work")
                        if not attempted:
                            if (retained != 0 or
                                    solve["recycle_operator_applies"] != 0 or
                                    solve["recycle_reduction_calls"] != 0 or
                                    accepted or projected != 0):
                                raise EvidenceError(
                                    f"line {line_number}: skipped projection has work")
                        elif (corrector != 2 or offered == 0 or
                              solve["recycle_reduction_calls"] == 0 or
                              solve["recycle_operator_applies"] !=
                              offered + (1 if retained else 0)):
                            raise EvidenceError(
                                f"line {line_number}: invalid projection work accounting")
                        if attempted and retained == 0:
                            if (accepted or projected != 0 or
                                    solve["recycle_reduction_calls"] != offered):
                                raise EvidenceError(
                                    f"line {line_number}: invalid all-deflated projection")
                        elif attempted:
                            improved = projected < solve["initial_true_residual"]
                            if accepted != improved:
                                raise EvidenceError(
                                    f"line {line_number}: invalid projection admission")
                            if (accepted and solve["iterations"] == 0 and
                                    solve["final_true_residual"] != projected):
                                raise EvidenceError(
                                    f"line {line_number}: invalid zero-iteration projection residual")
                    for scalar in ("initial_true_residual",
                                   "final_true_residual",
                                   "recursive_residual",
                                   "final_convergence_metric",
                                   "convergence_limit"):
                        value = solve.get(scalar)
                        if not isinstance(value, (int, float)) or \
                                isinstance(value, bool) or \
                                not math.isfinite(value) or value < 0:
                            raise EvidenceError(
                                f"line {line_number}: invalid pressure {scalar}")
                    audits = solve["convergence_audits"]
                    rejections = solve["convergence_rejections"]
                    coupled_contract = (
                        (is_v3 or is_v4 or is_v5 or is_modern) and
                        pressure_contract == "continuity_energy_coupled")
                    if corrector == 1 or coupled_contract:
                        if audits != 0 or rejections != 0:
                            raise EvidenceError(
                                f"line {line_number}: contract-incompatible supplemental audit")
                        if (solve["final_convergence_metric"] != 0 or
                                solve["convergence_limit"] != 0):
                            raise EvidenceError(
                                f"line {line_number}: unaudited solve has audit metric")
                    elif audits <= 0 or rejections > audits or \
                            solve["convergence_limit"] <= 0 or \
                            solve["final_convergence_metric"] > \
                            solve["convergence_limit"]:
                        raise EvidenceError(
                            f"line {line_number}: PISO2 continuity audit failed")
                if extension_enabled and \
                        pressure[1]["recycle_offered_directions"] != \
                        min(pressure[0]["recycle_cycle_corrections"], 4):
                    raise EvidenceError(
                        f"line {line_number}: projection lacks matching capture")
            if record.get("statistics_eligible"):
                if (record["launcher_ns"] <= 0 and
                        record["max_rank_step_ns"] <= 0):
                    raise EvidenceError(
                        f"line {line_number}: eligible record lacks timing authority")
    if seen == 0:
        raise EvidenceError("runtime evidence is empty")
    if restart_authority is not None and not restart_authority_used:
        raise EvidenceError(
            "--run-start-manifest does not authorize this evidence")


def percentile(values, probability):
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def paired_statistics(path):
    document = load_json(path)
    if (not isinstance(document, dict) or
            document.get("schema") != "HUNDUN_V04_PAIRED_INPUT_V1"):
        raise EvidenceError("paired input has the wrong schema")
    policy = document.get("policy")
    policy_sha256 = validate_performance_policy(policy)
    if document.get("policy_sha256") != policy_sha256:
        raise EvidenceError("paired input policy_sha256 mismatch")
    pairs = document.get("pairs") if isinstance(document, dict) else None
    if (not isinstance(pairs, list) or
            len(pairs) < policy["minimum_pairs"] or
            len(pairs) > policy["maximum_pairs"]):
        raise EvidenceError("paired statistics violate the frozen pair-count bounds")
    ratios = []  # type: List[float]
    prior_order = None
    hundun_steps = []  # type: List[float]
    coast_steps = []  # type: List[float]
    expected_samples = len(policy["full20_measured_steps"])
    for index, pair in enumerate(pairs):
        if not isinstance(pair, dict) or pair.get("order") not in ("HC", "CH"):
            raise EvidenceError(f"pair {index}: order must be HC or CH")
        expected_order = policy["first_pair_order"]
        if index % 2 == 1:
            expected_order = "CH" if expected_order == "HC" else "HC"
        if pair["order"] != expected_order or prior_order == pair["order"]:
            raise EvidenceError("launch order must alternate")
        prior_order = pair["order"]
        hundun = float(pair.get("hundun_hot_ns", 0.0))
        coast = float(pair.get("coast_hot_ns", 0.0))
        if not math.isfinite(hundun) or not math.isfinite(coast) or min(hundun, coast) <= 0.0:
            raise EvidenceError(f"pair {index}: timings must be finite and positive")
        ratios.append(hundun / coast)
        steps = pair.get("hundun_step_ns")
        coast_pair_steps = pair.get("coast_step_ns")
        if (not isinstance(steps, list) or len(steps) != expected_samples or
                not isinstance(coast_pair_steps, list) or
                len(coast_pair_steps) != expected_samples):
            raise EvidenceError(f"pair {index}: wrong measured-step sample count")
        converted_hundun = [float(value) for value in steps]
        converted_coast = [float(value) for value in coast_pair_steps]
        if (any(not math.isfinite(value) or value <= 0
                for value in converted_hundun + converted_coast) or
                abs(sum(converted_hundun) - hundun) > max(1.0, hundun * 1.0e-9) or
                abs(sum(converted_coast) - coast) > max(1.0, coast * 1.0e-9)):
            raise EvidenceError(f"pair {index}: hot total/sample mismatch")
        hundun_steps.extend(converted_hundun)
        coast_steps.extend(converted_coast)
    generator = random.Random(0x48464C4F57)
    medians = []
    for _ in range(20000):
        sample = [ratios[generator.randrange(len(ratios))]
                  for _ in range(len(ratios))]
        medians.append(statistics.median(sample))
    median = statistics.median(ratios)
    lower = percentile(medians, 0.025)
    upper = percentile(medians, 0.975)
    decision = "ACCEPT" if upper <= 1.0 else ("NEAR" if median <= 1.10 else "REJECT")
    return {"schema": "HUNDUN_V04_PAIRED_STATS_V1", "pairs": len(ratios),
            "policy_sha256": policy_sha256,
            "ratios": ratios, "median_ratio": median,
            "bootstrap_95": [lower, upper],
            "hundun_p90_step_ns": percentile(hundun_steps, 0.90),
            "coast_p90_step_ns": percentile(coast_steps, 0.90),
            "decision": decision}


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="hundun-v04-evidence-") as root_text:
        root = Path(root_text)
        binary = root / "hundun"
        tests_on = root / "tests-on"
        tests_off = root / "tests-off"
        case = root / "case.json"
        for path, payload in ((binary, b"binary"), (tests_on, b"tests-on"),
                              (tests_off, b"tests-off"), (case, b"{}\n")):
            path.write_bytes(payload)
        policy = {
            "schema": "HUNDUN_V04_PERFORMANCE_POLICY_V1",
            "grid": [480, 480, 48],
            "ranks": 64,
            "full2_steps": [1, 2],
            "full20_warmup_steps": list(range(1, 6)),
            "full20_measured_steps": list(range(6, 21)),
            "minimum_pairs": 5,
            "maximum_pairs": 9,
            "first_pair_order": "HC",
            "timing_authority": "max_rank_hot_step_sum",
            "serialized_output": False,
            "full2_max_directional_ratio": 1.25,
        }
        candidate = {key: "value" for key in REQUIRED_CANDIDATE}
        candidate.update({
            "head": "0" * 40,
            "tree": "1" * 40,
            "binary_path": str(binary.resolve()),
            "binary_sha256": digest_file(binary),
            "tests_on_path": str(tests_on.resolve()),
            "tests_on_sha256": digest_file(tests_on),
            "tests_off_path": str(tests_off.resolve()),
            "tests_off_sha256": digest_file(tests_off),
            "inputs": [{"path": str(case.resolve()),
                        "sha256": digest_file(case)}],
            "compiled_plan": 11,
            "cpu_plan": 12,
            "ranks": 64,
            "process_grid": [8, 8, 1],
            "cpu": {"model": "synthetic"},
            "numa": {"nodes": 1},
            "affinity": {"policy": "one-rank-per-core"},
            "environment": {},
            "command": [str(binary.resolve()), "run", str(case.resolve())],
            "output_state": {"restart": False, "visit": False},
            "checkpoint_state": {"kind": "fresh"},
            "start_timestamp": "2026-08-21T00:00:00Z",
            "performance_policy": policy,
        })
        candidate_path = root / "candidate.json"
        candidate_path.write_text(json.dumps(candidate), encoding="utf-8")
        manifest = root / "manifest.json"
        initialize(candidate_path, manifest)
        artifact = root / "focused.txt"
        artifact.write_text("accepted\n", encoding="utf-8")
        receipt = root / "receipt.json"
        focused_checks = {key: True for key in REQUIRED_GATE_CHECKS["focused"]}
        receipt.write_text(json.dumps({
            "schema": "HUNDUN_V04_GATE_RECEIPT_V1",
            "gate": "focused", "accepted": True,
            "candidate_sha256": load_json(manifest)["candidate_sha256"],
            "checks": focused_checks,
            "start_timestamp": "2026-08-21T00:00:01Z",
            "end_timestamp": "2026-08-21T00:00:02Z",
            "evidence": [{"path": str(artifact.resolve()),
                          "sha256": digest_file(artifact)}]}), encoding="utf-8")
        before = manifest.read_bytes()
        try:
            advance(manifest, "full2", receipt)
            raise AssertionError("out-of-order gate accepted")
        except EvidenceError:
            assert manifest.read_bytes() == before
        advance(manifest, "focused", receipt)
        validate_manifest(load_json(manifest))
        accepted = manifest.read_bytes()
        mutated = load_json(manifest)
        mutated["candidate"]["tree"] = "changed"
        try:
            validate_manifest(mutated)
            raise AssertionError("candidate mutation accepted")
        except EvidenceError:
            assert manifest.read_bytes() == accepted
        reused = root / "reused.json"
        reused.write_text(json.dumps({
            "schema": "HUNDUN_V04_GATE_RECEIPT_V1",
            "gate": "full2", "accepted": True,
            "candidate_sha256": load_json(manifest)["candidate_sha256"],
            "checks": {key: True for key in REQUIRED_GATE_CHECKS["full2"]},
            "start_timestamp": "2026-08-21T00:00:03Z",
            "end_timestamp": "2026-08-21T00:00:04Z",
            "evidence": [{"path": str(artifact.resolve()),
                          "sha256": digest_file(artifact)}]}), encoding="utf-8")
        try:
            advance(manifest, "full2", reused)
            raise AssertionError("reused evidence path accepted")
        except EvidenceError:
            assert manifest.read_bytes() == accepted
        pairs = root / "pairs.json"
        pair_rows = []
        sample_count = len(policy["full20_measured_steps"])
        for index in range(5):
            hundun_sample = 90 + index
            coast_sample = 100 + index
            pair_rows.append({
                "order": "HC" if index % 2 == 0 else "CH",
                "hundun_hot_ns": hundun_sample * sample_count,
                "coast_hot_ns": coast_sample * sample_count,
                "hundun_step_ns": [hundun_sample] * sample_count,
                "coast_step_ns": [coast_sample] * sample_count,
            })
        pairs.write_text(json.dumps({
            "schema": "HUNDUN_V04_PAIRED_INPUT_V1",
            "policy": policy,
            "policy_sha256": validate_performance_policy(policy),
            "pairs": pair_rows,
        }), encoding="utf-8")
        assert paired_statistics(pairs)["decision"] == "ACCEPT"

        pressure_base = {
            "status_code": 0, "termination": "converged", "iterations": 1,
            "initial_true_residual": 1.0, "final_true_residual": 1.0e-7,
            "recursive_residual": 1.0e-7,
            "final_convergence_metric": 0.0, "convergence_limit": 0.0,
            "convergence_audits": 0, "convergence_rejections": 0,
            "reduction_calls": 0, "operator_applies": 0,
            "preconditioner_applies": 0, "norm_breakdown_restarts": 0,
            "recycle_offered_directions": 0,
            "recycle_retained_directions": 0,
            "recycle_operator_applies": 0,
            "recycle_reduction_calls": 0,
            "recycle_projection_attempted": False,
            "recycle_projection_accepted": False,
            "recycle_projected_true_residual": 0.0,
            "recycle_cycle_corrections": 0,
            "recycle_capture_vector_passes": 0,
            "recycle_capture_cycle_attempts": 0,
            "recycle_capture_reduction_calls": 0,
            "recycle_capture_blocking_operations": 0,
        }
        linear_solve_base = {
            "status_code": 0, "termination": "converged", "iterations": 1,
            "initial_true_residual": 1.0, "final_true_residual": 1.0e-7,
            "recursive_residual": 1.0e-7, "reduction_calls": 0,
            "operator_applies": 0, "preconditioner_applies": 0,
            "norm_breakdown_restarts": 0,
        }
        runtime = {
            # Frozen historical V3 fixture.  This remains a literal validator
            # oracle and is not regenerated through the current V5 writer.
            "schema": V3_RUNTIME_SCHEMA,
            "build": 1, "binary": 2, "case": 3, "product": 4,
            "stl": 0, "cpu_plan": 5, "step": 1, "time": 0.1,
            "requested_bdf_order": 1, "bdf_order": 1,
            "temporal_method_fallback": False,
            "thermophysical_predictor_calls": 1,
            "startup": False, "retry": False,
            "restart_recovery": False, "statistics_eligible": False,
            "stages": [{"id": stage, "min_ns": 1, "mean_ns": 1,
                        "max_ns": 1} for stage in FROZEN_STAGE_IDS],
            "pressure_solve_contract": "continuity_energy_coupled",
            "terminal_physical_audit": {
                "present": True, "final_flux_revision": 17,
                "eos_residual": 1.0e-9, "eos_tolerance": 1.0e-6,
                "continuity_residual": 2.0e-9,
                "continuity_tolerance": 1.0e-6,
                "energy_residual": 3.0e-9,
                "energy_tolerance": 1.0e-6,
                "closed_mass_residual": 4.0e-9,
                "closed_mass_tolerance": 1.0e-6,
                "gauge_residual": 5.0e-9,
                "gauge_tolerance": 1.0e-6,
            },
            "momentum_predictor_solve_calls": 3,
            "momentum_predictor_limiter": {
                "limited": False, "theta": 1.0, "activations": 0,
            },
            "thermophysical_predictor": {
                "limited": False, "theta": 1.0,
                "low_state": "none", "mass_flux_scale": 1.0,
                "constraint": "none", "limiting_rank": -1,
                "limiting_global_cell": [0, 0, 0],
                "low_margin": 0.0, "high_margin": 0.0,
                "low_order_substeps": 0,
                "low_order_transport_passes": 0,
                "low_order_halo_exchanges": 0,
                "blocking_collectives": 1,
                "enthalpy_endpoint_alpha": 1.0,
                "bdf_endpoint_alpha": 1.0,
                "source_endpoint_alpha": 1.0,
                "enthalpy_solve_calls": 0,
            },
            "thermophysical_enthalpy_endpoint": dict(linear_solve_base),
            "momentum_predictor": [
                dict(linear_solve_base, component=component)
                for component in range(3)
            ],
        }
        for counter in (
                "launcher_ns", "max_rank_step_ns", "max_rank_rss_bytes",
                "max_node_rss_bytes", "structured_messages",
                "structured_bytes", "ibm_messages", "ibm_bytes",
                "blocking_collectives", "nonblocking_collectives",
                "reduction_ns", "linear_iterations", "exact_numeric_refills",
                "coarse_numeric_refills", "preconditioner_setups",
                "preconditioner_reuses", "heap_allocations"):
            runtime[counter] = 0
        runtime["blocking_collectives"] = 1
        runtime["pressure_solve_calls"] = 2
        runtime["pressure"] = [dict(pressure_base, corrector=1),
                               dict(pressure_base, corrector=2)]
        runtime_path = root / "runtime.jsonl"

        def reject_runtime(mutated, message, run_start_manifest=None):
            runtime_path.write_text(json.dumps(mutated) + "\n",
                                    encoding="utf-8")
            try:
                validate_runtime(runtime_path, run_start_manifest)
                raise AssertionError(message)
            except EvidenceError:
                pass

        def reject_runtime_rows(rows, message):
            runtime_path.write_text(
                "".join(json.dumps(row) + "\n" for row in rows),
                encoding="utf-8")
            try:
                validate_runtime(runtime_path)
                raise AssertionError(message)
            except EvidenceError:
                pass

        def reject_v5_contamination(fixture, schema_name):
            contaminated_cfl = json.loads(json.dumps(fixture))
            contaminated_cfl.setdefault("terminal_physical_audit", {})[
                "committed_convective_cfl"] = {
                    "out_max": 0.5, "abs_max": 0.375, "limit": 0.5,
                }
            reject_runtime(
                contaminated_cfl,
                f"{schema_name} record carrying V5 committed CFL accepted")
            for field, value in (
                    ("scheme", "common_face_afc_v2"),
                    ("retained_correction_l1_ratio", 1.0),
                    ("limited_faces", 0),
                    ("advective_cfl", {
                        "present": True, "plan": 31, "time_revision": 32,
                        "density_revision": 33, "face_flux_revision": 34,
                        "activity_collective": 0,
                        "dt": 0.006, "out_max": 0.5, "abs_max": 0.375,
                        "limit": 0.8,
                    })):
                contaminated_limiter = json.loads(json.dumps(fixture))
                contaminated_limiter.setdefault(
                    "momentum_predictor_limiter", {})[field] = value
                reject_runtime(
                    contaminated_limiter,
                    f"{schema_name} record carrying V5 limiter field "
                    f"{field} accepted")

        runtime_path.write_text(json.dumps(runtime) + "\n", encoding="utf-8")
        validate_runtime(runtime_path)
        reject_runtime([], "non-object runtime row accepted")

        runtime_v4 = json.loads(json.dumps(runtime))
        runtime_v4.update({
            "schema": V4_RUNTIME_SCHEMA,
            "pressure_energy_refinement_solve_calls": 0,
            "pressure_energy_refinement_termination":
                "component_residuals_converged",
            "pressure_energy_refinement": [],
        })
        runtime_path.write_text(json.dumps(runtime_v4) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        # Preserve the historical V4 global-theta oracle, including its
        # formerly legal full-collapse value.  V5 deliberately rejects that
        # state under the common-face AFC contract.
        frozen_v4_theta_zero = json.loads(json.dumps(runtime_v4))
        frozen_v4_theta_zero["momentum_predictor_limiter"].update({
            "limited": True, "theta": 0.0, "activations": 1,
        })
        runtime_path.write_text(json.dumps(frozen_v4_theta_zero) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        runtime_v5 = json.loads(json.dumps(runtime_v4))
        runtime_v5["schema"] = V5_RUNTIME_SCHEMA
        runtime_v5["terminal_physical_audit"][
            "committed_convective_cfl"] = {
                "out_max": 0.5, "abs_max": 0.375, "limit": 0.5,
            }
        runtime_v5["momentum_predictor_limiter"] = {
            "scheme": "common_face_afc_v2",
            "limited": False,
            "retained_correction_l1_ratio": 1.0,
            "limited_faces": 0,
            "advective_cfl": {
                "present": True,
                "plan": 31,
                "time_revision": 32,
                "density_revision": 33,
                "face_flux_revision": 34,
                "activity_collective": 0,
                "dt": 0.006,
                "out_max": 0.5,
                "abs_max": 0.375,
                "limit": 0.5,
            },
        }
        runtime_path.write_text(json.dumps(runtime_v5) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        cfl_limit = runtime_v5["terminal_physical_audit"][
            "committed_convective_cfl"]["limit"]
        cfl_edge = cfl_limit * (
            1.0 + 64.0 * float.fromhex("0x1.0p-52"))
        roundoff_edge_v5 = json.loads(json.dumps(runtime_v5))
        roundoff_edge_v5["terminal_physical_audit"][
            "committed_convective_cfl"]["out_max"] = cfl_edge
        runtime_path.write_text(json.dumps(roundoff_edge_v5) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)
        beyond_roundoff_v5 = json.loads(json.dumps(roundoff_edge_v5))
        beyond_roundoff_v5["terminal_physical_audit"][
            "committed_convective_cfl"]["out_max"] = (
                cfl_edge + float.fromhex("0x1.0p-52"))
        reject_runtime(beyond_roundoff_v5,
                       "V5 committed CFL beyond 64-epsilon slack accepted")

        valid_limited_v5 = json.loads(json.dumps(runtime_v5))
        valid_limited_v5["momentum_predictor_limiter"].update({
            "limited": True,
            "retained_correction_l1_ratio": 0.5,
            "limited_faces": 7,
        })
        runtime_path.write_text(json.dumps(valid_limited_v5) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        for field in V5_COMMITTED_CONVECTIVE_CFL_FIELDS:
            missing_v5_cfl = json.loads(json.dumps(runtime_v5))
            missing_v5_cfl["terminal_physical_audit"][
                "committed_convective_cfl"].pop(field)
            reject_runtime(missing_v5_cfl,
                           f"V5 committed CFL missing field {field}")
        for field, value in (
                ("out_max", -1.0), ("out_max", math.nan),
                ("abs_max", -1.0), ("abs_max", math.inf),
                ("limit", 0.0), ("limit", math.nan)):
            invalid_v5_cfl = json.loads(json.dumps(runtime_v5))
            invalid_v5_cfl["terminal_physical_audit"][
                "committed_convective_cfl"][field] = value
            reject_runtime(invalid_v5_cfl,
                           f"V5 invalid committed CFL {field} accepted")
        over_limit_v5_cfl = json.loads(json.dumps(runtime_v5))
        over_limit_v5_cfl["terminal_physical_audit"][
            "committed_convective_cfl"]["out_max"] = 0.5000000001
        reject_runtime(over_limit_v5_cfl,
                       "V5 over-limit committed CFL accepted")

        for field in V5_ADVECTIVE_CONVECTIVE_CFL_FIELDS:
            missing_advective = json.loads(json.dumps(runtime_v5))
            missing_advective["momentum_predictor_limiter"][
                "advective_cfl"].pop(field)
            reject_runtime(missing_advective,
                           f"V5 advective CFL missing field {field}")
        for field, value in (
                ("present", False), ("plan", 0), ("time_revision", 0),
                ("density_revision", 0), ("face_flux_revision", 0),
                ("activity_collective", -1),
                ("dt", 0.0), ("out_max", -1.0), ("abs_max", math.inf),
                ("limit", 0.0)):
            invalid_advective = json.loads(json.dumps(runtime_v5))
            invalid_advective["momentum_predictor_limiter"][
                "advective_cfl"][field] = value
            reject_runtime(invalid_advective,
                           f"V5 invalid advective CFL {field} accepted")
        same_revision_advective = json.loads(json.dumps(runtime_v5))
        same_revision_advective["momentum_predictor_limiter"][
            "advective_cfl"]["face_flux_revision"] = runtime_v5[
                "terminal_physical_audit"]["final_flux_revision"]
        reject_runtime(same_revision_advective,
                       "V5 identical advective/final revision accepted")
        mismatched_limit_advective = json.loads(json.dumps(runtime_v5))
        mismatched_limit_advective["momentum_predictor_limiter"][
            "advective_cfl"]["limit"] = 0.75
        reject_runtime(mismatched_limit_advective,
                       "V5 mismatched provisional/final CFL limits accepted")
        missing_ibm_activity = json.loads(json.dumps(runtime_v5))
        missing_ibm_activity["stl"] = 1
        reject_runtime(missing_ibm_activity,
                       "V5 IBM row without an activity collective accepted")
        fabricated_ibm_activity = json.loads(json.dumps(runtime_v5))
        fabricated_ibm_activity["momentum_predictor_limiter"][
            "advective_cfl"]["activity_collective"] = 1
        reject_runtime(fabricated_ibm_activity,
                       "V5 no-IBM row with an activity collective accepted")
        over_limit_advective = json.loads(json.dumps(runtime_v5))
        over_limit_advective["momentum_predictor_limiter"][
            "advective_cfl"]["out_max"] = 0.5000000001
        reject_runtime(over_limit_advective,
                       "V5 over-limit advective CFL accepted")

        second_v5 = json.loads(json.dumps(runtime_v5))
        second_v5["step"] = runtime_v5["step"] + 1
        second_v5["time"] = runtime_v5["time"] + 0.006
        second_v5["startup"] = False
        second_v5["momentum_predictor_limiter"]["advective_cfl"].update({
            "time_revision": 132,
            "density_revision": 133,
            "face_flux_revision": 134,
        })
        second_v5["terminal_physical_audit"]["final_flux_revision"] = 135
        runtime_path.write_text(
            json.dumps(runtime_v5) + "\n" + json.dumps(second_v5) + "\n",
            encoding="utf-8")
        validate_runtime(runtime_path)
        wrong_adjacent_dt = json.loads(json.dumps(second_v5))
        wrong_adjacent_dt["momentum_predictor_limiter"]["advective_cfl"][
            "dt"] = 0.003
        reject_runtime_rows(
            [runtime_v5, wrong_adjacent_dt],
            "V5 adjacent time delta/advective CFL dt mismatch accepted")
        changed_identity_v5 = json.loads(json.dumps(second_v5))
        changed_identity_v5["binary"] += 1
        reject_runtime_rows(
            [runtime_v5, changed_identity_v5],
            "V5 cross-binary row concatenation accepted")
        changed_cfl_plan_v5 = json.loads(json.dumps(second_v5))
        changed_cfl_plan_v5["momentum_predictor_limiter"]["advective_cfl"][
            "plan"] += 1
        reject_runtime_rows(
            [runtime_v5, changed_cfl_plan_v5],
            "V5 momentum CFL plan changed within one run")
        changed_cfl_limit_v5 = json.loads(json.dumps(second_v5))
        changed_cfl_limit_v5["momentum_predictor_limiter"]["advective_cfl"][
            "limit"] = 0.75
        changed_cfl_limit_v5["terminal_physical_audit"][
            "committed_convective_cfl"]["limit"] = 0.75
        reject_runtime_rows(
            [runtime_v5, changed_cfl_limit_v5],
            "V5 configured CFL ceiling changed within one run")
        ibm_first_v5 = json.loads(json.dumps(runtime_v5))
        ibm_first_v5["stl"] = 1
        ibm_first_v5["momentum_predictor_limiter"]["advective_cfl"][
            "activity_collective"] = 41
        ibm_second_v5 = json.loads(json.dumps(second_v5))
        ibm_second_v5["stl"] = 1
        ibm_second_v5["momentum_predictor_limiter"]["advective_cfl"][
            "activity_collective"] = 41
        runtime_path.write_text(
            json.dumps(ibm_first_v5) + "\n" +
            json.dumps(ibm_second_v5) + "\n", encoding="utf-8")
        validate_runtime(runtime_path)
        changed_activity_v5 = json.loads(json.dumps(ibm_second_v5))
        changed_activity_v5["momentum_predictor_limiter"]["advective_cfl"][
            "activity_collective"] = 42
        reject_runtime_rows(
            [ibm_first_v5, changed_activity_v5],
            "V5 IBM activity authority changed within one run")
        nonmonotone_time_v5 = json.loads(json.dumps(second_v5))
        nonmonotone_time_v5["time"] = runtime_v5["time"]
        reject_runtime_rows(
            [runtime_v5, nonmonotone_time_v5],
            "V5 non-increasing committed time accepted")
        skipped_step_v5 = json.loads(json.dumps(second_v5))
        skipped_step_v5["step"] += 1
        reject_runtime_rows(
            [runtime_v5, skipped_step_v5],
            "V5 non-contiguous committed steps accepted")
        interleaved_v4 = json.loads(json.dumps(runtime_v4))
        interleaved_v4["step"] = runtime_v5["step"] + 1
        interleaved_v4["time"] = runtime_v5["time"] + 0.006
        third_v5 = json.loads(json.dumps(second_v5))
        third_v5["step"] = runtime_v5["step"] + 2
        third_v5["time"] = runtime_v5["time"] + 0.012
        third_v5["binary"] += 1
        reject_runtime_rows(
            [runtime_v5, interleaved_v4, third_v5],
            "mixed-schema rows bypassed the V5 run identity contract")

        for field in V5_MOMENTUM_LIMITER_FIELDS:
            missing_v5_limiter = json.loads(json.dumps(runtime_v5))
            missing_v5_limiter["momentum_predictor_limiter"].pop(field)
            reject_runtime(missing_v5_limiter,
                           f"V5 momentum limiter missing field {field}")
        for limiter_update in (
                {"scheme": "global_theta_v1"},
                {"retained_correction_l1_ratio": 0.0},
                {"retained_correction_l1_ratio": 0.5},
                {"limited": True},
                {"limited_faces": 1}):
            invalid_v5_limiter = json.loads(json.dumps(runtime_v5))
            invalid_v5_limiter["momentum_predictor_limiter"].update(
                limiter_update)
            reject_runtime(invalid_v5_limiter,
                           f"V5 invalid common-face limiter accepted: "
                           f"{limiter_update}")

        v5_with_legacy_limiter = json.loads(json.dumps(runtime_v5))
        v5_with_legacy_limiter["momentum_predictor_limiter"]["theta"] = 1.0
        reject_runtime(v5_with_legacy_limiter,
                       "V5 record carrying V4 limiter field accepted")
        reject_v5_contamination(runtime_v4, "V4")
        reject_v5_contamination(runtime, "V3")

        def make_candidate_identity(executable, runtime_schema,
                                    build_manifest="3" * 64):
            if runtime_schema == RUNTIME_SCHEMA:
                identity_schema = \
                    "HUNDUN_V04_RUNTIME_CANDIDATE_IDENTITY_V2"
                evidence_binding = f"evidence_schema={RUNTIME_SCHEMA}\n"
            elif runtime_schema == V6_RUNTIME_SCHEMA:
                identity_schema = \
                    "HUNDUN_V04_RUNTIME_CANDIDATE_IDENTITY_V1"
                evidence_binding = ""
            else:
                raise AssertionError("unsupported candidate identity schema")
            identity = {
                "schema": identity_schema,
                "head": "1" * 40,
                "tree": "2" * 40,
                "build_manifest_sha256": build_manifest,
                "executable_sha256": executable,
            }
            payload = (
                f"schema={identity_schema}\n"
                f"{evidence_binding}"
                f"head={identity['head']}\n"
                f"tree={identity['tree']}\n"
                "build_manifest_sha256="
                f"{identity['build_manifest_sha256']}\n"
                f"executable_sha256={identity['executable_sha256']}\n")
            identity["identity_sha256"] = hashlib.sha256(
                payload.encode("utf-8")).hexdigest()
            return identity

        runtime_v7 = json.loads(json.dumps(runtime_v5))
        runtime_v7.update({
            "schema": RUNTIME_SCHEMA,
            "candidate_identity": make_candidate_identity(
                "4" * 64, RUNTIME_SCHEMA),
            "build": int(("3" * 64)[:16], 16),
            "binary": int(("4" * 64)[:16], 16),
            "run_start": {
                "kind": "fresh", "previous_step": 0,
                "previous_time": 0.0,
                "restart_manifest_sha256": None,
            },
            "step": 1,
            "previous_committed_time": 0.0,
            "time": 0.5,
            "startup": True,
            "restart_recovery": False,
        })
        runtime_v7["terminal_physical_audit"][
            "committed_convective_cfl"] = {
                "valid": True,
                "density_revision": 33,
                "final_flux_revision": 17,
                "density_view": {
                    "field": 0, "collective": 41,
                },
                "face_flux_view": {
                    "collective": 42,
                },
                "activity_collective": 0,
                "dt": 0.5,
                "out_max": 0.5,
                "abs_max": 0.375,
                "limit": 0.5,
                "out_winner": {
                    "valid": True, "global_cell": [0, 0, 0], "rank": 0,
                    "out": 0.5, "abs": 0.375, "density_volume": 1.0,
                    "outgoing_mass_flow": 1.0,
                    "absolute_mass_flow": 1.5,
                },
                "abs_winner": {
                    "valid": True, "global_cell": [0, 0, 0], "rank": 0,
                    "out": 0.5, "abs": 0.375, "density_volume": 1.0,
                    "outgoing_mass_flow": 1.0,
                    "absolute_mass_flow": 1.5,
                },
            }
        runtime_v7["momentum_predictor_limiter"] = {
            "scheme": "common_face_afc_v3_owner",
            "limited": False,
            "correction_metrics_applicability": "not_applicable",
            "retained_correction_l1_ratio": None,
            "minimum_face_alpha": None,
            "active_correction_faces": 0,
            "limited_faces": 0,
            "limited_face_fraction": None,
            "advective_cfl": {
                "present": True, "plan": 31,
                "time_revision_collective": 32,
                "density_view_collective": 33,
                "face_flux_view_collective": 34,
                "activity_collective": 0,
                "dt": 0.5, "out_max": 0.5, "abs_max": 0.375,
                "limit": 0.5,
            },
        }
        runtime_path.write_text(json.dumps(runtime_v7) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        twelve_refinements_v7 = json.loads(json.dumps(runtime_v7))
        twelve_refinements_v7["pressure_energy_refinement_solve_calls"] = 12
        twelve_refinements_v7["linear_iterations"] = 12
        twelve_refinements_v7["pressure_energy_refinement"] = [
            dict(pressure_base, ordinal=ordinal, target_generation=71,
                 collective_lineage=80 + ordinal)
            for ordinal in range(1, 13)
        ]
        runtime_path.write_text(json.dumps(twelve_refinements_v7) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        thirteen_refinements_v7 = json.loads(json.dumps(
            twelve_refinements_v7))
        thirteen_refinements_v7[
            "pressure_energy_refinement_solve_calls"] = 13
        thirteen_refinements_v7["pressure_energy_refinement"].append(dict(
            pressure_base, ordinal=13, target_generation=71,
            collective_lineage=93))
        reject_runtime(
            thirteen_refinements_v7,
            "V7 refinement prefix above the twelve-solve bound accepted")

        frozen_six_refinements_v6 = json.loads(json.dumps(
            twelve_refinements_v7))
        frozen_six_refinements_v6["schema"] = V6_RUNTIME_SCHEMA
        frozen_six_refinements_v6["candidate_identity"] = \
            make_candidate_identity("4" * 64, V6_RUNTIME_SCHEMA)
        frozen_six_refinements_v6[
            "pressure_energy_refinement_solve_calls"] = 6
        frozen_six_refinements_v6["pressure_energy_refinement"] = \
            frozen_six_refinements_v6[
                "pressure_energy_refinement"][:6]
        runtime_path.write_text(json.dumps(frozen_six_refinements_v6) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)
        seventh_refinement_v6 = json.loads(json.dumps(
            frozen_six_refinements_v6))
        seventh_refinement_v6[
            "pressure_energy_refinement_solve_calls"] = 7
        seventh_refinement_v6["pressure_energy_refinement"].append(dict(
            pressure_base, ordinal=7, target_generation=71,
            collective_lineage=87))
        reject_runtime(
            seventh_refinement_v6,
            "frozen V6 refinement prefix above six accepted")

        relabelled_v6_as_v7 = json.loads(json.dumps(
            frozen_six_refinements_v6))
        relabelled_v6_as_v7["schema"] = RUNTIME_SCHEMA
        reject_runtime(
            relabelled_v6_as_v7,
            "frozen V6 row relabelled as V7 without a V2 identity accepted")

        first_anchor_mismatch = json.loads(json.dumps(runtime_v7))
        first_anchor_mismatch["time"] = 0.006
        first_anchor_mismatch["momentum_predictor_limiter"][
            "advective_cfl"]["dt"] = 0.003
        first_cfl = first_anchor_mismatch["terminal_physical_audit"][
            "committed_convective_cfl"]
        first_cfl.update({"dt": 0.003, "out_max": 0.003,
                          "abs_max": 0.00225})
        for winner in (first_cfl["out_winner"], first_cfl["abs_winner"]):
            winner.update({"out": 0.003, "abs": 0.00225})
        reject_runtime(
            first_anchor_mismatch,
            "V7 first row time=0.006 with dt=0.003 accepted")

        invalid_identity_hash = json.loads(json.dumps(runtime_v7))
        invalid_identity_hash["candidate_identity"][
            "identity_sha256"] = "0" * 64
        reject_runtime(invalid_identity_hash,
                       "V7 invalid candidate identity hash accepted")
        invalid_winner_formula = json.loads(json.dumps(runtime_v7))
        invalid_winner_formula["terminal_physical_audit"][
            "committed_convective_cfl"]["out_winner"][
                "outgoing_mass_flow"] = 2.0
        reject_runtime(invalid_winner_formula,
                       "V7 invalid committed CFL winner formula accepted")
        aliased_cfl_view = json.loads(json.dumps(runtime_v7))
        aliased_cfl_view["momentum_predictor_limiter"]["advective_cfl"][
            "face_flux_view_collective"] = 42
        reject_runtime(aliased_cfl_view,
                       "V7 aliased provisional/committed CFL view accepted")
        fabricated_inactive_metric = json.loads(json.dumps(runtime_v7))
        fabricated_inactive_metric["momentum_predictor_limiter"][
            "minimum_face_alpha"] = 0.0
        reject_runtime(fabricated_inactive_metric,
                       "V7 inactive AFC zero metric accepted")

        applicable_v7 = json.loads(json.dumps(runtime_v7))
        applicable_v7["momentum_predictor_limiter"].update({
            "limited": True,
            "correction_metrics_applicability": "applicable",
            "retained_correction_l1_ratio": 0.5,
            "minimum_face_alpha": 0.25,
            "active_correction_faces": 14,
            "limited_faces": 7,
            "limited_face_fraction": 0.5,
        })
        runtime_path.write_text(json.dumps(applicable_v7) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        second_v7 = json.loads(json.dumps(runtime_v7))
        second_v7.update({"step": 2, "previous_committed_time": 0.5,
                          "time": 1.0, "startup": False})
        second_v7["momentum_predictor_limiter"]["advective_cfl"].update({
            "time_revision_collective": 132,
            "density_view_collective": 133,
            "face_flux_view_collective": 134,
        })
        second_v7["terminal_physical_audit"]["final_flux_revision"] = 135
        second_v7["terminal_physical_audit"]["committed_convective_cfl"].update({
            "density_revision": 133, "final_flux_revision": 135,
        })
        runtime_path.write_text(
            json.dumps(runtime_v7) + "\n" + json.dumps(second_v7) + "\n",
            encoding="utf-8")
        validate_runtime(runtime_path)

        mixed_candidate_v7 = json.loads(json.dumps(second_v7))
        mixed_candidate_v7["candidate_identity"] = make_candidate_identity(
            "5" * 64, RUNTIME_SCHEMA, "6" * 64)
        mixed_candidate_v7["build"] = int(("6" * 64)[:16], 16)
        mixed_candidate_v7["binary"] = int(("5" * 64)[:16], 16)
        reject_runtime_rows(
            [runtime_v7, mixed_candidate_v7],
            "V7 mixed executable/build candidate identities accepted")

        restart_v7 = json.loads(json.dumps(runtime_v7))
        restart_v7.update({"step": 10, "previous_committed_time": 4.5,
                           "time": 5.0, "startup": False,
                           "restart_recovery": True})
        restart_manifest_path = root / "restart-manifest.bin"
        restart_manifest = bytearray(struct.pack(
            "<8sIIiiiQQQdddQQI", b"H4MANI01", 1, 1, 8, 8, 8,
            11, 12, 13, 4.5, 0.5, 101325.0, 9, 14, 3))
        restart_manifest.extend(struct.pack("<BHB", 0, 0, 3))
        restart_manifest.extend(struct.pack("<BHB", 1, 1, 1))
        restart_manifest.extend(struct.pack("<BHB", 3, 2, 1))
        restart_manifest.extend(struct.pack(
            "<iiiiiiQQ", 0, 0, 0, 8, 8, 8, 16, 1))
        checksum = 1469598103934665603
        for value in restart_manifest:
            checksum = ((checksum ^ value) * 1099511628211) & UINT64_MAX
        restart_manifest.extend(struct.pack("<Q", checksum))
        restart_manifest_path.write_bytes(restart_manifest)
        restart_manifest_sha256 = hashlib.sha256(restart_manifest).hexdigest()
        restart_v7["run_start"] = {
            "kind": "restart", "previous_step": 9,
            "previous_time": 4.5,
            "restart_manifest_sha256": restart_manifest_sha256,
        }
        continuation_v7 = json.loads(json.dumps(second_v7))
        continuation_v7.update({"step": 11, "previous_committed_time": 5.0,
                                "time": 5.5, "startup": False,
                                "restart_recovery": False})
        continuation_v7["run_start"] = json.loads(json.dumps(
            restart_v7["run_start"]))
        runtime_path.write_text(
            json.dumps(restart_v7) + "\n" +
            json.dumps(continuation_v7) + "\n", encoding="utf-8")
        validate_runtime(runtime_path, restart_manifest_path)

        false_exact_history = json.loads(json.dumps(restart_v7))
        false_exact_history.update({
            "requested_bdf_order": 2,
            "bdf_order": 2,
            "restart_recovery": False,
        })
        reject_runtime(
            false_exact_history,
            "legacy Restart manifest accepted as exact BDF2 history",
            restart_manifest_path)

        reject_runtime(
            restart_v7,
            "V7 restart evidence accepted without its frozen manifest")

        shifted_restart_anchor = json.loads(json.dumps(restart_v7))
        shifted_restart_anchor["previous_committed_time"] += 100.0
        shifted_restart_anchor["time"] += 100.0
        reject_runtime(
            shifted_restart_anchor,
            "V7 restart row time and previous time shifted away from the "
            "snapshot anchor", restart_manifest_path)

        shifted_entire_anchor = json.loads(json.dumps(restart_v7))
        shifted_entire_anchor["run_start"]["previous_time"] += 100.0
        shifted_entire_anchor["previous_committed_time"] += 100.0
        shifted_entire_anchor["time"] += 100.0
        reject_runtime(
            shifted_entire_anchor,
            "V7 restart anchor and row times shifted away from the same "
            "manifest", restart_manifest_path)

        v5_with_v6_identity = json.loads(json.dumps(runtime_v5))
        v5_with_v6_identity["candidate_identity"] = make_candidate_identity(
            "4" * 64, V6_RUNTIME_SCHEMA)
        reject_runtime(v5_with_v6_identity,
                       "frozen V5 accepted a V6 candidate identity")

        one_refinement = json.loads(json.dumps(runtime_v4))
        one_refinement["linear_iterations"] = 1
        one_refinement["pressure_energy_refinement_solve_calls"] = 1
        one_refinement["pressure_energy_refinement"] = [dict(
            pressure_base, ordinal=1, target_generation=71,
            collective_lineage=81)]
        runtime_path.write_text(json.dumps(one_refinement) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        two_refinements = json.loads(json.dumps(one_refinement))
        two_refinements["linear_iterations"] = 2
        two_refinements["pressure_energy_refinement_solve_calls"] = 2
        two_refinements["pressure_energy_refinement"].append(dict(
            pressure_base, ordinal=2, target_generation=71,
            collective_lineage=82))

        wrong_refinement_order = json.loads(json.dumps(two_refinements))
        wrong_refinement_order["pressure_energy_refinement"][0]["ordinal"] = 2
        wrong_refinement_order["pressure_energy_refinement"][1]["ordinal"] = 1
        reject_runtime(wrong_refinement_order,
                       "V4 out-of-order refinement accepted")

        wrong_refinement_target = json.loads(json.dumps(two_refinements))
        wrong_refinement_target["pressure_energy_refinement"][1][
            "target_generation"] = 72
        reject_runtime(wrong_refinement_target,
                       "V4 mixed refinement targets accepted")

        wrong_refinement_lineage = json.loads(json.dumps(two_refinements))
        wrong_refinement_lineage["pressure_energy_refinement"][1][
            "collective_lineage"] = 81
        reject_runtime(wrong_refinement_lineage,
                       "V4 duplicate refinement lineage accepted")

        wrong_refinement_count = json.loads(json.dumps(one_refinement))
        wrong_refinement_count[
            "pressure_energy_refinement_solve_calls"] = 2
        reject_runtime(wrong_refinement_count,
                       "V4 refinement count/array mismatch accepted")

        wrong_refinement_termination = json.loads(json.dumps(one_refinement))
        wrong_refinement_termination[
            "pressure_energy_refinement_termination"] = \
            "iteration_capacity_exhausted"
        reject_runtime(wrong_refinement_termination,
                       "V4 non-converged refinement termination accepted")

        nonfinite_refinement = json.loads(json.dumps(one_refinement))
        nonfinite_refinement["pressure_energy_refinement"][0][
            "final_true_residual"] = math.nan
        reject_runtime(nonfinite_refinement,
                       "V4 non-finite refinement residual accepted")

        missing_refinement_resource = json.loads(json.dumps(one_refinement))
        missing_refinement_resource["linear_iterations"] = 0
        reject_runtime(missing_refinement_resource,
                       "V4 refinement omitted from linear resources")

        for v4_field in V4_REFINEMENT_FIELDS:
            contaminated_v3 = json.loads(json.dumps(runtime))
            contaminated_v3[v4_field] = runtime_v4[v4_field]
            reject_runtime(contaminated_v3,
                           f"V3 carried V4 field {v4_field}")

        for v4_field in V4_REFINEMENT_FIELDS:
            missing_v4_field = json.loads(json.dumps(runtime_v4))
            missing_v4_field.pop(v4_field)
            reject_runtime(missing_v4_field,
                           f"V4 missing refinement field {v4_field}")

        for required_field in (
                "schema", "build", "binary", "case", "stl", "product", "cpu_plan",
                "step", "time", "requested_bdf_order", "bdf_order",
                "temporal_method_fallback", "thermophysical_predictor_calls",
                "launcher_ns", "max_rank_step_ns", "max_rank_rss_bytes",
                "max_node_rss_bytes", "structured_messages",
                "structured_bytes", "ibm_messages", "ibm_bytes",
                "blocking_collectives", "nonblocking_collectives",
                "reduction_ns", "linear_iterations", "exact_numeric_refills",
                "coarse_numeric_refills", "preconditioner_setups",
                "preconditioner_reuses", "heap_allocations",
                "pressure_solve_calls", "pressure_solve_contract",
                "terminal_physical_audit", "momentum_predictor_solve_calls",
                "momentum_predictor_limiter", "thermophysical_predictor",
                "thermophysical_enthalpy_endpoint", "pressure",
                "momentum_predictor", "startup", "retry",
                "restart_recovery", "statistics_eligible", "stages"):
            missing_field = json.loads(json.dumps(runtime))
            missing_field.pop(required_field)
            reject_runtime(missing_field,
                           f"V3 missing required field accepted: {required_field}")

        for object_field, required_fields in (
                ("terminal_physical_audit", (
                    "present", "final_flux_revision", "eos_residual",
                    "eos_tolerance", "continuity_residual",
                    "continuity_tolerance", "energy_residual",
                    "energy_tolerance", "closed_mass_residual",
                    "closed_mass_tolerance", "gauge_residual",
                    "gauge_tolerance")),
                ("momentum_predictor_limiter",
                 ("limited", "theta", "activations")),
                ("thermophysical_predictor", (
                    "limited", "theta", "low_state", "mass_flux_scale",
                    "constraint", "limiting_rank", "limiting_global_cell",
                    "low_margin", "high_margin", "low_order_substeps",
                    "low_order_transport_passes", "low_order_halo_exchanges",
                    "blocking_collectives", "enthalpy_endpoint_alpha",
                    "bdf_endpoint_alpha", "source_endpoint_alpha",
                    "enthalpy_solve_calls")),
                ("thermophysical_enthalpy_endpoint", (
                    "status_code", "termination", "iterations",
                    "initial_true_residual", "final_true_residual",
                    "recursive_residual", "reduction_calls",
                    "operator_applies", "preconditioner_applies",
                    "norm_breakdown_restarts"))):
            for required_field in required_fields:
                missing_field = json.loads(json.dumps(runtime))
                missing_field[object_field].pop(required_field)
                reject_runtime(
                    missing_field,
                    f"V3 {object_field} missing field accepted: {required_field}")

        for required_field in (
                "component", "status_code", "termination", "iterations",
                "initial_true_residual", "final_true_residual",
                "recursive_residual", "reduction_calls", "operator_applies",
                "preconditioner_applies", "norm_breakdown_restarts"):
            missing_field = json.loads(json.dumps(runtime))
            missing_field["momentum_predictor"][0].pop(required_field)
            reject_runtime(
                missing_field,
                f"V3 momentum solve missing field accepted: {required_field}")

        for required_field in ("id", "min_ns", "mean_ns", "max_ns"):
            missing_field = json.loads(json.dumps(runtime))
            missing_field["stages"][0].pop(required_field)
            reject_runtime(
                missing_field,
                f"V3 stage timing missing field accepted: {required_field}")

        for required_field in (
                "corrector", "status_code", "termination", "iterations",
                "initial_true_residual", "final_true_residual",
                "recursive_residual", "final_convergence_metric",
                "convergence_limit", "convergence_audits",
                "convergence_rejections", "reduction_calls",
                "operator_applies", "preconditioner_applies",
                "norm_breakdown_restarts", "recycle_offered_directions",
                "recycle_retained_directions", "recycle_operator_applies",
                "recycle_reduction_calls", "recycle_projection_attempted",
                "recycle_projection_accepted",
                "recycle_projected_true_residual",
                "recycle_cycle_corrections",
                "recycle_capture_vector_passes",
                "recycle_capture_cycle_attempts",
                "recycle_capture_reduction_calls",
                "recycle_capture_blocking_operations"):
            missing_field = json.loads(json.dumps(runtime))
            missing_field["pressure"][0].pop(required_field)
            reject_runtime(
                missing_field,
                f"V3 pressure solve missing field accepted: {required_field}")

        for integer_field in (
                "build", "binary", "case", "stl", "product", "cpu_plan",
                "step", "requested_bdf_order", "bdf_order",
                "thermophysical_predictor_calls", "launcher_ns",
                "max_rank_step_ns", "max_rank_rss_bytes",
                "max_node_rss_bytes", "structured_messages",
                "structured_bytes", "ibm_messages", "ibm_bytes",
                "blocking_collectives", "nonblocking_collectives",
                "reduction_ns", "linear_iterations", "exact_numeric_refills",
                "coarse_numeric_refills", "preconditioner_setups",
                "preconditioner_reuses", "heap_allocations",
                "pressure_solve_calls", "momentum_predictor_solve_calls"):
            wrong_type = json.loads(json.dumps(runtime))
            wrong_type[integer_field] = True
            reject_runtime(
                wrong_type,
                f"V3 boolean integer accepted: {integer_field}")

        nonfinite_time = json.loads(json.dumps(runtime))
        nonfinite_time["time"] = math.nan
        reject_runtime(nonfinite_time, "V3 non-finite time accepted")

        for boolean_field in (
                "temporal_method_fallback", "startup", "retry",
                "restart_recovery", "statistics_eligible"):
            wrong_type = json.loads(json.dumps(runtime))
            wrong_type[boolean_field] = 0
            reject_runtime(
                wrong_type,
                f"V3 integer boolean accepted: {boolean_field}")

        nested_type_cases = (
            ("momentum_predictor_limiter", "limited", 0),
            ("momentum_predictor_limiter", "theta", True),
            ("momentum_predictor_limiter", "activations", True),
            ("thermophysical_predictor", "limited", 0),
            ("thermophysical_predictor", "theta", True),
            ("thermophysical_predictor", "low_state", 0),
            ("thermophysical_predictor", "mass_flux_scale", True),
            ("thermophysical_predictor", "constraint", 0),
            ("thermophysical_predictor", "limiting_rank", True),
            ("thermophysical_predictor", "low_margin", True),
            ("thermophysical_predictor", "high_margin", True),
            ("thermophysical_predictor", "low_order_substeps", True),
            ("thermophysical_predictor", "low_order_transport_passes", True),
            ("thermophysical_predictor", "low_order_halo_exchanges", True),
            ("thermophysical_predictor", "blocking_collectives", True),
            ("thermophysical_predictor", "enthalpy_endpoint_alpha", True),
            ("thermophysical_predictor", "bdf_endpoint_alpha", True),
            ("thermophysical_predictor", "source_endpoint_alpha", True),
            ("thermophysical_predictor", "enthalpy_solve_calls", True),
        )
        for object_field, field, invalid_value in nested_type_cases:
            wrong_type = json.loads(json.dumps(runtime))
            wrong_type[object_field][field] = invalid_value
            reject_runtime(
                wrong_type,
                f"V3 invalid nested type accepted: {object_field}.{field}")

        for invalid_cell in ([0, 0], [0, True, 0], "0,0,0"):
            wrong_type = json.loads(json.dumps(runtime))
            wrong_type["thermophysical_predictor"][
                "limiting_global_cell"] = invalid_cell
            reject_runtime(wrong_type,
                           "V3 invalid limiting-global-cell accepted")

        for solve_field, invalid_value in (
                ("status_code", True), ("termination", 0),
                ("iterations", True), ("initial_true_residual", True),
                ("final_true_residual", True),
                ("recursive_residual", True), ("reduction_calls", True),
                ("operator_applies", True),
                ("preconditioner_applies", True),
                ("norm_breakdown_restarts", True)):
            wrong_endpoint = json.loads(json.dumps(runtime))
            wrong_endpoint["thermophysical_enthalpy_endpoint"][
                solve_field] = invalid_value
            reject_runtime(
                wrong_endpoint,
                f"V3 invalid enthalpy endpoint type accepted: {solve_field}")
            wrong_momentum = json.loads(json.dumps(runtime))
            wrong_momentum["momentum_predictor"][0][
                solve_field] = invalid_value
            reject_runtime(
                wrong_momentum,
                f"V3 invalid momentum solve type accepted: {solve_field}")

        wrong_component = json.loads(json.dumps(runtime))
        wrong_component["momentum_predictor"][0]["component"] = True
        reject_runtime(wrong_component,
                       "V3 boolean momentum component accepted")

        wrong_corrector = json.loads(json.dumps(runtime))
        wrong_corrector["pressure"][0]["corrector"] = True
        reject_runtime(wrong_corrector, "V3 boolean pressure corrector accepted")
        wrong_pressure_count = json.loads(json.dumps(runtime))
        wrong_pressure_count["pressure"][0]["iterations"] = True
        reject_runtime(wrong_pressure_count,
                       "V3 boolean pressure count accepted")
        wrong_stage_timing = json.loads(json.dumps(runtime))
        wrong_stage_timing["stages"][0]["min_ns"] = True
        reject_runtime(wrong_stage_timing, "V3 boolean stage timing accepted")

        oversized_endpoint_status = json.loads(json.dumps(runtime))
        oversized_endpoint_status["thermophysical_enthalpy_endpoint"][
            "status_code"] = 1 << 16
        reject_runtime(oversized_endpoint_status,
                       "V3 status code wider than uint16 accepted")
        oversized_momentum_iterations = json.loads(json.dumps(runtime))
        oversized_momentum_iterations["momentum_predictor"][0][
            "iterations"] = 1 << 32
        reject_runtime(oversized_momentum_iterations,
                       "V3 momentum iterations wider than uint32 accepted")
        oversized_pressure_iterations = json.loads(json.dumps(runtime))
        oversized_pressure_iterations["pressure"][0][
            "iterations"] = 1 << 32
        reject_runtime(oversized_pressure_iterations,
                       "V3 pressure iterations wider than uint32 accepted")

        for temporal_updates in (
                {"requested_bdf_order": 3},
                {"bdf_order": 3},
                {"requested_bdf_order": 2, "bdf_order": 1},
                {"thermophysical_predictor_calls": 2},
                {"temporal_method_fallback": True},
                {"blocking_collectives": 0}):
            invalid_temporal = json.loads(json.dumps(runtime))
            invalid_temporal.update(temporal_updates)
            reject_runtime(invalid_temporal,
                           f"V3 invalid temporal relation accepted: {temporal_updates}")

        valid_temporal_fallback = json.loads(json.dumps(runtime))
        valid_temporal_fallback.update({
            "requested_bdf_order": 2,
            "bdf_order": 1,
            "thermophysical_predictor_calls": 2,
            "temporal_method_fallback": True,
            "blocking_collectives": 2,
        })
        valid_temporal_fallback["thermophysical_predictor"][
            "blocking_collectives"] = 2
        runtime_path.write_text(json.dumps(valid_temporal_fallback) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        invalid_momentum_calls = json.loads(json.dumps(runtime))
        invalid_momentum_calls["momentum_predictor_solve_calls"] = 2
        reject_runtime(invalid_momentum_calls,
                       "V3 non-three momentum solve count accepted")
        for limiter_updates in (
                {"theta": -0.1}, {"theta": 1.1}, {"theta": 0.5},
                {"limited": True}, {"activations": 1}):
            invalid_limiter = json.loads(json.dumps(runtime))
            invalid_limiter["momentum_predictor_limiter"].update(
                limiter_updates)
            reject_runtime(invalid_limiter,
                           f"V3 invalid momentum limiter accepted: {limiter_updates}")

        for solve_updates in (
                {"component": 1}, {"status_code": 1},
                {"termination": "maximum_iterations"},
                {"initial_true_residual": 1.0e-8,
                 "final_true_residual": 1.0e-7}):
            invalid_momentum = json.loads(json.dumps(runtime))
            invalid_momentum["momentum_predictor"][0].update(solve_updates)
            reject_runtime(invalid_momentum,
                           f"V3 invalid momentum solve accepted: {solve_updates}")

        valid_limited_momentum = json.loads(json.dumps(runtime))
        valid_limited_momentum["momentum_predictor_limiter"].update({
            "limited": True, "theta": 0.5, "activations": 1,
        })
        runtime_path.write_text(json.dumps(valid_limited_momentum) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        for predictor_updates in (
                {"theta": 0.5}, {"mass_flux_scale": 0.5},
                {"low_margin": 1.0}, {"high_margin": -1.0},
                {"low_order_substeps": 1},
                {"low_order_transport_passes": 1},
                {"low_order_halo_exchanges": 1},
                {"limiting_rank": 0}, {"constraint": "density"},
                {"low_state": "scaled_euler"},
                {"enthalpy_endpoint_alpha": 0.5},
                {"bdf_endpoint_alpha": 0.5},
                {"source_endpoint_alpha": 0.5},
                {"enthalpy_solve_calls": 1},
                {"constraint": "invalid"}, {"low_state": "invalid"}):
            invalid_predictor = json.loads(json.dumps(runtime))
            invalid_predictor["thermophysical_predictor"].update(
                predictor_updates)
            reject_runtime(invalid_predictor,
                           f"V3 invalid fast predictor accepted: {predictor_updates}")

        limited_predictor = json.loads(json.dumps(runtime))
        limited_predictor["blocking_collectives"] = 3
        limited_predictor["thermophysical_predictor"].update({
            "limited": True, "theta": 0.5,
            "low_state": "scaled_euler", "mass_flux_scale": 0.25,
            "constraint": "density", "limiting_rank": 0,
            "limiting_global_cell": [1, 2, 3],
            "low_margin": 1.0, "high_margin": -1.0,
            "low_order_substeps": 1,
            "low_order_transport_passes": 1,
            "low_order_halo_exchanges": 0,
            "blocking_collectives": 3,
        })
        runtime_path.write_text(json.dumps(limited_predictor) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        for predictor_updates in (
                {"theta": 1.0}, {"theta": -0.1}, {"theta": 1.1},
                {"mass_flux_scale": -0.1}, {"mass_flux_scale": 1.1},
                {"low_margin": -0.1}, {"high_margin": 0.0},
                {"low_order_substeps": 0},
                {"low_order_transport_passes": 0},
                {"low_order_halo_exchanges": 1},
                {"blocking_collectives": 1}, {"limiting_rank": -1},
                {"limiting_global_cell": [-1, 2, 3]},
                {"constraint": "none"}, {"low_state": "none"}):
            invalid_predictor = json.loads(json.dumps(limited_predictor))
            invalid_predictor["thermophysical_predictor"].update(
                predictor_updates)
            reject_runtime(
                invalid_predictor,
                f"V3 invalid limited predictor accepted: {predictor_updates}")

        bdf_homotopy = json.loads(json.dumps(limited_predictor))
        bdf_homotopy["thermophysical_predictor"].update({
            "low_state": "bdf_accepted_rate_homotopy",
            "mass_flux_scale": 1.0,
            "bdf_endpoint_alpha": 0.375,
        })
        runtime_path.write_text(json.dumps(bdf_homotopy) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        implicit_upwind = json.loads(json.dumps(limited_predictor))
        implicit_upwind["thermophysical_predictor"].update({
            "low_state": "implicit_upwind", "mass_flux_scale": 1.0,
            "enthalpy_solve_calls": 1,
        })
        implicit_upwind["thermophysical_enthalpy_endpoint"].update({
            "reduction_calls": 7, "operator_applies": 4,
        })
        runtime_path.write_text(json.dumps(implicit_upwind) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        implicit_source_limited = json.loads(json.dumps(implicit_upwind))
        implicit_source_limited["thermophysical_predictor"].update({
            "low_state": "implicit_upwind_source_limited",
            "enthalpy_endpoint_alpha": 0.625,
        })
        runtime_path.write_text(json.dumps(implicit_source_limited) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        donor_source_limited = json.loads(json.dumps(limited_predictor))
        donor_source_limited["thermophysical_predictor"].update({
            "low_state": "bdf_local_donor_flux_source_limited",
            "low_order_halo_exchanges": 1,
            "source_endpoint_alpha": 0.375,
        })
        runtime_path.write_text(json.dumps(donor_source_limited) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        invalid_implicit_endpoint = json.loads(json.dumps(implicit_upwind))
        invalid_implicit_endpoint["thermophysical_enthalpy_endpoint"][
            "operator_applies"] = 0
        reject_runtime(invalid_implicit_endpoint,
                       "V3 implicit predictor without endpoint apply accepted")
        invalid_source_alpha = json.loads(json.dumps(donor_source_limited))
        invalid_source_alpha["thermophysical_predictor"][
            "source_endpoint_alpha"] = 1.0
        reject_runtime(invalid_source_alpha,
                       "V3 source-limited predictor with unit alpha accepted")

        maximum_revision = json.loads(json.dumps(runtime))
        maximum_revision["terminal_physical_audit"][
            "final_flux_revision"] = UINT64_MAX
        runtime_path.write_text(json.dumps(maximum_revision) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)
        oversized_revision = json.loads(json.dumps(runtime))
        oversized_revision["terminal_physical_audit"][
            "final_flux_revision"] = UINT64_MAX + 1
        reject_runtime(oversized_revision,
                       "V3 oversized final-flux revision accepted")

        lossy_metric_comparison = json.loads(json.dumps(runtime))
        lossy_metric_comparison["terminal_physical_audit"].update({
            "eos_residual": (1 << 53) + 1,
            "eos_tolerance": 1 << 53,
        })
        reject_runtime(lossy_metric_comparison,
                       "V3 residual above tolerance hidden by FP64 cast")

        widened_coupled_energy_gate = json.loads(json.dumps(runtime))
        widened_coupled_energy_gate["terminal_physical_audit"][
            "energy_tolerance"] = 2.0e-6
        reject_runtime(widened_coupled_energy_gate,
                       "V3 coupled energy tolerance widened independently")

        historical_v2 = json.loads(json.dumps(runtime))
        historical_v2["schema"] = V2_RUNTIME_SCHEMA
        historical_v2.pop("pressure_solve_contract")
        historical_v2.pop("terminal_physical_audit")
        historical_v2["pressure"][1].update({
            "convergence_audits": 1,
            "final_convergence_metric": 5.0e-7,
            "convergence_limit": 1.0e-6,
        })
        for solve in historical_v2["pressure"]:
            for field in PRESSURE_EXTENSION_FIELDS:
                solve.pop(field)
        runtime_path.write_text(json.dumps(historical_v2) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)
        reject_v5_contamination(historical_v2, "V2")

        historical_v2_extended = json.loads(json.dumps(runtime))
        historical_v2_extended["schema"] = V2_RUNTIME_SCHEMA
        historical_v2_extended.pop("pressure_solve_contract")
        historical_v2_extended.pop("terminal_physical_audit")
        historical_v2_extended["pressure"][1].update({
            "convergence_audits": 1,
            "final_convergence_metric": 5.0e-7,
            "convergence_limit": 1.0e-6,
        })
        runtime_path.write_text(json.dumps(historical_v2_extended) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        missing_contract = json.loads(json.dumps(runtime))
        missing_contract.pop("pressure_solve_contract")
        reject_runtime(missing_contract, "V3 missing pressure contract accepted")

        unknown_contract = json.loads(json.dumps(runtime))
        unknown_contract["pressure_solve_contract"] = "pressure_only"
        reject_runtime(unknown_contract, "V3 unknown pressure contract accepted")

        missing_terminal = json.loads(json.dumps(runtime))
        missing_terminal.pop("terminal_physical_audit")
        reject_runtime(missing_terminal, "V3 missing terminal audit accepted")

        absent_terminal = json.loads(json.dumps(runtime))
        absent_terminal["terminal_physical_audit"]["present"] = False
        reject_runtime(absent_terminal, "V3 absent terminal audit accepted")

        for invalid_revision in (0, True, -1):
            bad_revision = json.loads(json.dumps(runtime))
            bad_revision["terminal_physical_audit"][
                "final_flux_revision"] = invalid_revision
            reject_runtime(bad_revision,
                           "V3 invalid final-flux revision accepted")

        for metric in ("eos", "continuity", "energy", "closed_mass",
                       "gauge"):
            over_limit = json.loads(json.dumps(runtime))
            over_limit["terminal_physical_audit"][
                f"{metric}_residual"] = 2.0e-6
            reject_runtime(over_limit,
                           f"V3 over-limit {metric} residual accepted")
            for suffix in ("residual", "tolerance"):
                for invalid_number in (-1.0, math.nan, math.inf, True):
                    invalid_metric = json.loads(json.dumps(runtime))
                    invalid_metric["terminal_physical_audit"][
                        f"{metric}_{suffix}"] = invalid_number
                    reject_runtime(
                        invalid_metric,
                        f"V3 invalid {metric} {suffix} accepted")

        disabled_coupled_energy = json.loads(json.dumps(runtime))
        disabled_coupled_energy["terminal_physical_audit"][
            "energy_tolerance"] = 0.0
        reject_runtime(disabled_coupled_energy,
                       "V3 coupled energy gate disabled")

        fabricated_coupled_audit = json.loads(json.dumps(runtime))
        fabricated_coupled_audit["pressure"][1].update({
            "convergence_audits": 1,
            "final_convergence_metric": 5.0e-7,
            "convergence_limit": 1.0e-6,
        })
        reject_runtime(fabricated_coupled_audit,
                       "V3 coupled supplemental audit accepted")

        pressure_continuity_v3 = json.loads(json.dumps(runtime))
        pressure_continuity_v3["pressure_solve_contract"] = \
            "pressure_continuity"
        pressure_continuity_v3["pressure"][1].update({
            "convergence_audits": 1,
            "final_convergence_metric": 5.0e-7,
            "convergence_limit": 1.0e-6,
        })
        pressure_continuity_v3["terminal_physical_audit"].update({
            "energy_residual": 2.0,
            "energy_tolerance": 0.0,
        })
        runtime_path.write_text(json.dumps(pressure_continuity_v3) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        pressure_continuity_energy_enabled = json.loads(
            json.dumps(pressure_continuity_v3))
        pressure_continuity_energy_enabled["terminal_physical_audit"].update({
            "energy_residual": 5.0e-7,
            "energy_tolerance": 1.0e-6,
        })
        runtime_path.write_text(
            json.dumps(pressure_continuity_energy_enabled) + "\n",
            encoding="utf-8")
        validate_runtime(runtime_path)

        missing_pressure_continuity_audit = json.loads(
            json.dumps(pressure_continuity_v3))
        missing_pressure_continuity_audit["pressure"][1].update({
            "convergence_audits": 0,
            "final_convergence_metric": 0.0,
            "convergence_limit": 0.0,
        })
        reject_runtime(missing_pressure_continuity_audit,
                       "V3 pressure-continuity C2 audit omission accepted")

        for field in PRESSURE_EXTENSION_FIELDS:
            partial_both = json.loads(json.dumps(runtime))
            partial_both["pressure"][0].pop(field)
            partial_both["pressure"][1].pop(field)
            reject_runtime(partial_both,
                           f"partial pressure extension accepted: {field}")

        inconsistent_presence = json.loads(json.dumps(runtime))
        for field in PRESSURE_EXTENSION_FIELDS:
            inconsistent_presence["pressure"][1].pop(field)
        reject_runtime(inconsistent_presence,
                       "C1/C2 extension presence mismatch accepted")

        legal_capture = json.loads(json.dumps(runtime))
        legal_capture["pressure"][0].update({
            "recycle_cycle_corrections": 1,
            "recycle_capture_vector_passes": 4,
            "recycle_capture_cycle_attempts": 2,
            "recycle_capture_reduction_calls": 2,
            "recycle_capture_blocking_operations": 4,
        })
        legal_capture["pressure"][1]["recycle_offered_directions"] = 1
        runtime_path.write_text(json.dumps(legal_capture) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        legal_accepted = json.loads(json.dumps(runtime))
        legal_accepted["pressure"][0].update({
            "recycle_cycle_corrections": 2,
            "recycle_capture_vector_passes": 4,
            "recycle_capture_cycle_attempts": 2,
            "recycle_capture_reduction_calls": 2,
            "recycle_capture_blocking_operations": 4,
        })
        legal_accepted["pressure"][1].update({
            "iterations": 0,
            "initial_true_residual": 1.0,
            "final_true_residual": 0.25,
            "recursive_residual": 0.25,
            "recycle_offered_directions": 2,
            "recycle_retained_directions": 1,
            "recycle_operator_applies": 3,
            "recycle_reduction_calls": 7,
            "recycle_projection_attempted": True,
            "recycle_projection_accepted": True,
            "recycle_projected_true_residual": 0.25,
        })
        runtime_path.write_text(json.dumps(legal_accepted) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        # Exact frozen identity-operator fixture: r0=e0+e1, U=span(e0),
        # projected residual=e1, then one ordinary audited FGMRES iteration.
        legal_accepted_ordinary = json.loads(json.dumps(runtime))
        legal_accepted_ordinary["pressure"][0].update({
            "recycle_cycle_corrections": 1,
            "recycle_capture_vector_passes": 2,
            "recycle_capture_cycle_attempts": 1,
            "recycle_capture_reduction_calls": 1,
            "recycle_capture_blocking_operations": 2,
        })
        legal_accepted_ordinary["pressure"][1].update({
            "iterations": 1,
            "initial_true_residual": 1.4142135623730951,
            "final_true_residual": 0.0,
            "recursive_residual": 0.0,
            "operator_applies": 3,
            "preconditioner_applies": 1,
            "reduction_calls": 7,
            "recycle_offered_directions": 1,
            "recycle_retained_directions": 1,
            "recycle_operator_applies": 2,
            "recycle_reduction_calls": 5,
            "recycle_projection_attempted": True,
            "recycle_projection_accepted": True,
            "recycle_projected_true_residual": 1.0,
        })
        runtime_path.write_text(json.dumps(legal_accepted_ordinary) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        legal_all_deflated = json.loads(json.dumps(runtime))
        legal_all_deflated["pressure"][0].update({
            "recycle_cycle_corrections": 2,
            "recycle_capture_vector_passes": 4,
            "recycle_capture_cycle_attempts": 2,
            "recycle_capture_reduction_calls": 2,
            "recycle_capture_blocking_operations": 4,
        })
        legal_all_deflated["pressure"][1].update({
            "initial_true_residual": 1.0,
            "final_true_residual": 1.0,
            "recursive_residual": 1.0,
            "recycle_offered_directions": 2,
            "recycle_operator_applies": 2,
            "recycle_reduction_calls": 2,
            "recycle_projection_attempted": True,
            "recycle_projected_true_residual": 0.0,
        })
        runtime_path.write_text(json.dumps(legal_all_deflated) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        legal_non_improving = json.loads(json.dumps(runtime))
        legal_non_improving["pressure"][0].update({
            "recycle_cycle_corrections": 1,
            "recycle_capture_vector_passes": 2,
            "recycle_capture_cycle_attempts": 1,
            "recycle_capture_reduction_calls": 1,
            "recycle_capture_blocking_operations": 2,
        })
        legal_non_improving["pressure"][1].update({
            "initial_true_residual": 1.0,
            "final_true_residual": 1.0,
            "recursive_residual": 1.0,
            "recycle_offered_directions": 1,
            "recycle_retained_directions": 1,
            "recycle_operator_applies": 2,
            "recycle_reduction_calls": 5,
            "recycle_projection_attempted": True,
            "recycle_projected_true_residual": 1.0,
        })
        runtime_path.write_text(json.dumps(legal_non_improving) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        missing_pressure = dict(runtime)
        missing_pressure.pop("pressure")
        missing_pressure.pop("pressure_solve_calls")
        runtime_path.write_text(json.dumps(missing_pressure) + "\n",
                                encoding="utf-8")
        try:
            validate_runtime(runtime_path)
            raise AssertionError("V2 runtime without pressure evidence accepted")
        except EvidenceError:
            pass

        legacy = dict(missing_pressure, schema=LEGACY_RUNTIME_SCHEMA)
        runtime_path.write_text(json.dumps(legacy) + "\n", encoding="utf-8")
        validate_runtime(runtime_path)
        reject_v5_contamination(legacy, "V1")

        invalid_pressure = json.loads(json.dumps(runtime))
        invalid_pressure["pressure"][1]["final_convergence_metric"] = 2.0e-6
        runtime_path.write_text(json.dumps(invalid_pressure) + "\n",
                                encoding="utf-8")
        try:
            validate_runtime(runtime_path)
            raise AssertionError("rejected PISO2 continuity evidence accepted")
        except EvidenceError:
            pass

        invalid_role = json.loads(json.dumps(runtime))
        invalid_role["pressure"][0]["recycle_projection_attempted"] = True
        invalid_role["pressure"][0]["recycle_offered_directions"] = 1
        reject_runtime(invalid_role, "PISO1 projection work accepted")

        invalid_skip = json.loads(json.dumps(runtime))
        invalid_skip["pressure"][1]["recycle_operator_applies"] = 1
        reject_runtime(invalid_skip, "skipped projection work accepted")

        valid_skip = json.loads(json.dumps(runtime))
        valid_skip["pressure"][0].update({
            "recycle_cycle_corrections": 2,
            "recycle_capture_vector_passes": 4,
            "recycle_capture_cycle_attempts": 2,
            "recycle_capture_reduction_calls": 2,
            "recycle_capture_blocking_operations": 4,
        })
        valid_skip["pressure"][1]["recycle_offered_directions"] = 2
        runtime_path.write_text(json.dumps(valid_skip) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        legal_capped = json.loads(json.dumps(runtime))
        legal_capped["pressure"][0].update({
            "recycle_cycle_corrections": 5,
            "recycle_capture_vector_passes": 10,
            "recycle_capture_cycle_attempts": 5,
            "recycle_capture_reduction_calls": 5,
            "recycle_capture_blocking_operations": 10,
        })
        legal_capped["pressure"][1]["recycle_offered_directions"] = 4
        runtime_path.write_text(json.dumps(legal_capped) + "\n",
                                encoding="utf-8")
        validate_runtime(runtime_path)

        projection_without_capture = json.loads(json.dumps(runtime))
        projection_without_capture["pressure"][1][
            "recycle_offered_directions"] = 1
        reject_runtime(projection_without_capture,
                       "projection without C1 capture accepted")

        invalid_work = json.loads(json.dumps(runtime))
        invalid_work["pressure"][1].update({
            "recycle_offered_directions": 1,
            "recycle_projection_attempted": True,
            "recycle_reduction_calls": 1,
            "recycle_operator_applies": 2,
            "recycle_projected_true_residual": 1.0,
        })
        reject_runtime(invalid_work, "all-deflated projection work mismatch accepted")

        invalid_operator_formula = json.loads(json.dumps(legal_accepted))
        invalid_operator_formula["pressure"][1]["recycle_operator_applies"] = 2
        reject_runtime(invalid_operator_formula,
                       "projection operator formula mismatch accepted")

        invalid_admission = json.loads(json.dumps(runtime))
        invalid_admission["pressure"][1].update({
            "recycle_offered_directions": 1,
            "recycle_retained_directions": 1,
            "recycle_projection_attempted": True,
            "recycle_projection_accepted": True,
            "recycle_reduction_calls": 1,
            "recycle_operator_applies": 2,
            "recycle_projected_true_residual": 1.0,
            "final_true_residual": 1.0,
        })
        reject_runtime(invalid_admission,
                       "non-improving projection admission accepted")

        improved_but_rejected = json.loads(json.dumps(legal_accepted))
        improved_but_rejected["pressure"][1][
            "recycle_projection_accepted"] = False
        reject_runtime(improved_but_rejected,
                       "improved but rejected projection accepted")

        zero_iteration_mismatch = json.loads(json.dumps(legal_accepted))
        zero_iteration_mismatch["pressure"][1]["final_true_residual"] = 0.1
        reject_runtime(zero_iteration_mismatch,
                       "zero-iteration final/projected mismatch accepted")

        cap_exceeded = json.loads(json.dumps(legal_capped))
        cap_exceeded["pressure"][1]["recycle_offered_directions"] = 5
        reject_runtime(cap_exceeded, "projection cap above four accepted")

        offered_too_few = json.loads(json.dumps(valid_skip))
        offered_too_few["pressure"][1]["recycle_offered_directions"] = 1
        reject_runtime(offered_too_few,
                       "offered below min(published,4) accepted")

        offered_too_many = json.loads(json.dumps(valid_skip))
        offered_too_many["pressure"][1]["recycle_offered_directions"] = 3
        reject_runtime(offered_too_many,
                       "offered above min(published,4) accepted")

        invalid_capture = json.loads(json.dumps(runtime))
        invalid_capture["pressure"][0].update({
            "recycle_cycle_corrections": 1,
            "recycle_capture_vector_passes": 0,
            "recycle_capture_cycle_attempts": 1,
            "recycle_capture_reduction_calls": 2,
        })
        reject_runtime(invalid_capture, "partial capture evidence accepted")

        invalid_capture_blocking = json.loads(json.dumps(runtime))
        invalid_capture_blocking["pressure"][1][
            "recycle_capture_blocking_operations"] = 1
        reject_runtime(invalid_capture_blocking,
                       "PISO2 capture blocking work accepted")


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    init_parser = subparsers.add_parser("init")
    init_parser.add_argument("--candidate", type=Path, required=True)
    init_parser.add_argument("--output", type=Path, required=True)
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("manifest", type=Path)
    gate_parser = subparsers.add_parser("advance")
    gate_parser.add_argument("manifest", type=Path)
    gate_parser.add_argument("--gate", choices=GATES, required=True)
    gate_parser.add_argument("--receipt", type=Path, required=True)
    runtime_parser = subparsers.add_parser("runtime")
    runtime_parser.add_argument("evidence", type=Path)
    runtime_parser.add_argument("--run-start-manifest", type=Path)
    paired_parser = subparsers.add_parser("paired-stats")
    paired_parser.add_argument("pairs", type=Path)
    policy_parser = subparsers.add_parser("policy")
    policy_parser.add_argument("path", type=Path)
    equivalence_parser = subparsers.add_parser("equivalence")
    equivalence_parser.add_argument("path", type=Path)
    equivalence_parser.add_argument("--require-sealed", action="store_true")
    subparsers.add_parser("self-test")
    arguments = parser.parse_args()
    if arguments.command is None:
        parser.error("a command is required")
    try:
        if arguments.command == "init":
            initialize(arguments.candidate, arguments.output)
        elif arguments.command == "validate":
            validate_manifest(load_json(arguments.manifest))
        elif arguments.command == "advance":
            advance(arguments.manifest, arguments.gate, arguments.receipt)
        elif arguments.command == "runtime":
            validate_runtime(arguments.evidence,
                             arguments.run_start_manifest)
        elif arguments.command == "paired-stats":
            print(json.dumps(paired_statistics(arguments.pairs), indent=2,
                             sort_keys=True))
        elif arguments.command == "policy":
            policy = load_json(arguments.path)
            print(validate_performance_policy(policy))
        elif arguments.command == "equivalence":
            print(validate_equivalence(arguments.path,
                                       arguments.require_sealed))
        else:
            self_test()
    except (EvidenceError, OSError, json.JSONDecodeError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
