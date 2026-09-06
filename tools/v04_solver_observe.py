#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Attribute one instrumented run by step/attempt/C1-C2/refinement; not a speed gate."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def ratio(final, initial):
    return float(final) / float(initial) if float(initial) != 0 else None


def summarize(root):
    paths = sorted(root.glob("solver-rank-*.csv"))
    if not paths:
        raise ValueError("no per-loop observations; old aggregate costs cannot supply C1/C2 attribution")
    grouped = defaultdict(list)
    by_step_rank = defaultdict(list)
    for path in paths:
        for row in rows(path):
            if int(row["dropped_loops"]):
                raise ValueError("diagnostic loop capacity exceeded; attribution is incomplete")
            key = tuple(int(row[k]) for k in ("step", "attempt", "corrector", "refinement", "kind"))
            grouped[key].append(row)
            by_step_rank[(int(row["step"]), int(row["rank"]))].append(row)
    performance = rows(root / "performance.csv")
    for row in performance:
        if int(row["dropped_loops"]):
            raise ValueError("diagnostic loop capacity exceeded")
        loops = by_step_rank[(int(row["step"]), int(row["rank"]))]
        if sum(int(v["invoked"]) for v in loops) != sum(int(row[k + "_calls"]) for k in ("pressure", "diagonal", "spatial")):
            raise ValueError("solve-count attribution mismatch")
        for detailed, total in (("A_ns", "A_apply_ns"), ("M_ns", "M_apply_ns")):
            if sum(int(v[detailed]) for v in loops) != int(row[total]):
                raise ValueError("linear timing attribution mismatch: " + detailed)
    summaries = []
    for key, group in sorted(grouped.items()):
        if len({row["rank"] for row in group}) != len(group):
            raise ValueError("duplicate rank at the same loop")
        timings = [k for k in group[0] if k.endswith("_ns")]
        summaries.append({"step": key[0], "attempt": key[1], "corrector": key[2],
            "refinement": key[3], "kind": ("pressure", "diagonal", "spatial")[key[4]],
            "ranks": len(group),
            "counts_by_rank": {k: sorted({int(row[k]) for row in group}) for k in
                ("invoked", "iterations", "A_calls", "M_calls", "baseline_candidates",
                 "extrapolated_candidates", "ladder_candidates", "incomplete_candidates")},
            "rank_mean_ns": {k: sum(int(row[k]) for row in group) / len(group) for k in timings},
            "rank_max_ns": {k: max(int(row[k]) for row in group) for k in timings},
            "linear_contraction": [ratio(row["linear_final"], row["linear_initial"]) for row in group],
            "physical_contraction": [{"rank": int(row["rank"]),
                "continuity": ratio(row["selected_continuity"], row["baseline_continuity"]),
                "energy": ratio(row["selected_energy"], row["baseline_energy"])}
                for row in group if int(row["globalization_valid"])]})
    added = ("final_momentum_ns", "terminal_metrics_ns", "boundary_ledger_ns",
             "structured_wait_ns", "structured_control_ns")
    return {"schema": "HUNDUN_LOOP_OBSERVATION_V1", "formal_speed_comparison": False,
        "files_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in paths + [root / "performance.csv", root / "conservation.csv"]},
        "loops": summaries,
        "rank_step_mean_ns": {k: sum(int(row[k]) for row in performance) / len(performance) for k in added},
        "scope": ["M/A include their communication; do not add them again to solve time",
                  "MG copy is a subset of refill, itself in solve preparation",
                  "Structured wait/control excludes IBM donor and MG-specific MPI; not total MPI time",
                  "Candidate assembly is outside the associated linear solve",
                  "No independent phase maxima are summed into full-step wall time",
                  "64-loop capacity per advance; nonzero drop count is rejected",
                  "Additional diagnostics do not change physical acceptance thresholds"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(summarize(args.run_root), indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
