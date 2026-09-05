#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Summarize one diagnostic window, without creating a benchmark acceptance."""
import argparse
import bisect
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess

STEP_PHASES = ("advance_ns", "observables_ns", "visit_ns", "evidence_resources_ns",
               "csv_ns", "checkpoint_ns")
CANDIDATE_PHASES = ("state_copy_ns", "velocity_halo_ns", "thermo_ns",
                    "boundary_derived_ns", "flux_ns", "certificate_ns",
                    "residual_ns", "equivalence_hash_ns")


def step_summary(rows):
    if not rows or len({row["rank"] for row in rows}) != len(rows):
        raise ValueError("empty step or duplicate rank")
    for row in rows:
        if (any(value < 0 for value in row.values()) or
                sum(row[key] for key in STEP_PHASES) != row["full_step_ns"] or
                sum(row[key] for key in CANDIDATE_PHASES) > row["candidate_ns"]):
            raise ValueError("phase accounting mismatch")
    critical = max(rows, key=lambda row: row["full_step_ns"])
    return {"step": critical["step"], "ranks": len(rows),
            "maximum_full_step_ns": critical["full_step_ns"],
            "maximum_advance_ns": max(row["advance_ns"] for row in rows),
            "mean_advance_ns": sum(row["advance_ns"] for row in rows) / len(rows),
            "full_step_critical_rank": critical["rank"],
            "critical_rank_phases": {key: critical[key] for key in STEP_PHASES},
            "rank_mean_candidate_phases": {
                key: sum(row[key] for row in rows) / len(rows)
                for key in CANDIDATE_PHASES + ("candidate_ns", "schur_prepare_ns",
                                              "krylov_ns", "recovery_ns")},
            "work_counts": {key: sorted({row[key] for row in rows}) for key in
                ("baseline_evaluations", "extrapolation_evaluations", "ladder_evaluations",
                 "incomplete_evaluations", "rejected_extrapolations")}}


def symbols(binary):
    text = subprocess.check_output(["nm", "-n", "-S", "--defined-only", "-C", str(binary)],
                                   universal_newlines=True)
    entries = []
    for line in text.splitlines():
        match = re.match(r"^([0-9a-fA-F]+)\s+([0-9a-fA-F]+)\s+[TtWw]\s+(.+)$", line)
        if match:
            entries.append((int(match[1], 16), int(match[2], 16), match[3]))
    return sorted(entries)


def summarize(arguments):
    with arguments.performance.open() as stream:
        rows = [{key: int(value) for key, value in row.items()}
                for row in csv.DictReader(stream)]
    steps = sorted({row["step"] for row in rows})
    if not steps:
        raise ValueError("no performance samples")
    result = {"schema": "HUNDUN_V04_DIAGNOSTIC_OBSERVATION_V1",
              "formal_performance_comparison": False,
              "steps": [step_summary([row for row in rows if row["step"] == step])
                        for step in steps],
              "binary_sha256": hashlib.sha256(arguments.binary.read_bytes()).hexdigest(),
              "limits": ["One instrumented window, no repeated-round median",
                         "Full-step timer excludes its own observation gather/write",
                         "PMPI elapsed excludes observer bookkeeping; advance includes it",
                         "Nine blocking MPI primitives, observed only inside driver.advance",
                         "Candidate/solve details describe the latest numerical attempt",
                         "Per-rank phase maxima must not be added as whole-step wall time"]}
    if arguments.mpi_directory:
        entries = symbols(arguments.binary)
        starts = [item[0] for item in entries]
        observed = {}
        owners = {}
        for path in sorted(arguments.mpi_directory.glob("rank-*.tsv")):
            rank = int(path.stem.split("-")[1])
            lines = path.read_text().splitlines()
            if not lines or lines[-1] != "# dropped_calls=0":
                raise ValueError("incomplete/dropped MPI records: {}".format(path))
            for record in csv.DictReader(lines[:-1], delimiter="\t"):
                key = (int(record["step"]), rank)
                calls, elapsed = int(record["calls"]), int(record["ns"])
                observed[key] = observed.get(key, 0) + calls
                site = int(record["site"], 16) - 1  # A call's return address.
                index = bisect.bisect_right(starts, site) - 1
                owner = "unresolved:" + record["module"] + ":" + record["site"]
                if (Path(record["module"]).name == arguments.binary.name and index >= 0 and
                        site < entries[index][0] + entries[index][1]):
                    owner = entries[index][2]
                bucket = owners.setdefault((owner, record["kind"]), {"calls": 0, "ns": 0})
                bucket["calls"] += calls
                bucket["ns"] += elapsed
        if set(observed) != {(row["step"], row["rank"]) for row in rows}:
            raise ValueError("MPI and performance step/rank coverage mismatch")
        result["mpi_counts"] = [{"step": row["step"], "rank": row["rank"],
            "observed": observed[(row["step"], row["rank"])],
            "reported_subtotal": row["reported_collective_subtotal"],
            "gap": observed[(row["step"], row["rank"])] - row["reported_collective_subtotal"]}
            for row in rows]
        result["mpi_callers_all_rank_sum"] = [dict(owner=owner, kind=kind, **value)
            for (owner, kind), value in sorted(owners.items(), key=lambda item: -item[1]["ns"])]
    payload = json.dumps(result, indent=2) + "\n"
    if arguments.output:
        with arguments.output.open("x") as stream:
            stream.write(payload)
    else:
        print(payload, end="")


def self_test():
    row = dict.fromkeys(STEP_PHASES + CANDIDATE_PHASES + (
        "candidate_ns", "schur_prepare_ns", "krylov_ns", "recovery_ns",
        "baseline_evaluations", "extrapolation_evaluations", "ladder_evaluations",
        "incomplete_evaluations", "rejected_extrapolations"), 0)
    row.update(step=1, rank=0, advance_ns=10, visit_ns=20, full_step_ns=30)
    other = dict(row, rank=1, advance_ns=20, visit_ns=0, full_step_ns=20)
    summary = step_summary([row, other])
    assert summary["maximum_full_step_ns"] == 30
    assert summary["maximum_advance_ns"] == 20
    assert sum(summary["critical_rank_phases"].values()) == 30
    for invalid in ([row, row], [dict(row, full_step_ns=29)]):
        try:
            step_summary(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid accounting accepted")
    print("performance accounting self-test PASS")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--performance", type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--mpi-directory", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    if arguments.self_test:
        self_test()
    elif not arguments.performance or not arguments.binary:
        parser.error("--performance and --binary are required")
    else:
        summarize(arguments)


if __name__ == "__main__":
    main()
