#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Validate frozen coverage, then attribute a diagnostic window (not a speed gate)."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path

COUNTS = ("invoked", "iterations", "A_calls", "M_calls", "baseline_candidates",
          "extrapolated_candidates", "ladder_candidates", "incomplete_candidates")
ADDED = ("final_momentum_ns", "terminal_metrics_ns", "boundary_ledger_ns",
         "structured_wait_ns", "structured_control_ns", "scalar_remap_ns")
FLOATS = ("dt", "linear_initial", "linear_final", "baseline_continuity",
          "baseline_energy", "selected_continuity", "selected_energy", "selected_alpha")
SOURCE = "source_meta_sha256"


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metadata(root, legacy_ranks):
    path = root / "RUN.meta"
    if path.stat().st_size > 65536:
        raise ValueError("oversized RUN.meta")
    lines = path.read_text().splitlines()
    if not lines or lines[0] != "HUNDUN_V04_THIN_DOMAIN_RUN_V1" or lines[-1] != "end":
        raise ValueError("invalid RUN.meta envelope")
    values = {}
    for line in lines[1:-1]:
        key, value = line.split(" ", 1)
        if key in values:
            raise ValueError("duplicate metadata key: " + key)
        values[key] = value
    if values.get("observation_schema") not in (None, "3"):
        raise ValueError("unsupported observation schema")
    bound = values.get("observation_schema") == "3"
    if bound:
        ranks = int(values["expected_ranks"])
        if legacy_ranks is not None and legacy_ranks != ranks:
            raise ValueError("CLI rank count conflicts with frozen metadata")
    elif legacy_ranks is not None:
        ranks = legacy_ranks
    else:
        raise ValueError("legacy observations need --allow-partial --expected-ranks; coverage is not frozen")
    first = int(values["starting_step"]) + 1
    count = int(values["requested_steps"])
    if ranks <= 0 or ranks > 2**31 - 1 or first <= 0 or count <= 0 or first + count - 1 > 2**64 - 1:
        raise ValueError("invalid frozen rank/step range")
    return ranks, first, first + count - 1, bound, file_hash(path)


def checked_rows(stream, source, loop, rank=None):
    def lines():
        for line in stream:
            if not line.endswith("\n"):
                raise ValueError("truncated CSV tail: " + stream.name)
            yield line
    reader = csv.DictReader(lines())
    required = {"step", "rank", "dropped_loops"}
    required.update(("attempt", "corrector", "refinement", "kind", "A_ns", "M_ns",
                     "linear_initial", "linear_final", "globalization_valid") + COUNTS if loop else
                    ("pressure_calls", "diagonal_calls", "spatial_calls", "A_apply_ns", "M_apply_ns"))
    header = reader.fieldnames
    if not header or len(set(header)) != len(header) or not required.issubset(header):
        raise ValueError("missing/duplicate CSV columns: " + stream.name)
    if source is not None and SOURCE not in header:
        raise ValueError("missing source identity: " + stream.name)
    for raw in reader:
        if None in raw or any(v is None or v == "" for v in raw.values()):
            raise ValueError("malformed CSV row: " + stream.name)
        row = {}
        for key, value in raw.items():
            if key == SOURCE:
                if source is not None and value != source:
                    raise ValueError("source identity mismatch: " + stream.name)
                continue
            number = float(value) if key in FLOATS else int(value)
            if not math.isfinite(number) or number < 0 or (key not in FLOATS and number > 2**64 - 1):
                raise ValueError("invalid count/time/residual: " + key)
            row[key] = number
        if row["step"] == 0 or (rank is not None and row["rank"] != rank):
            raise ValueError("invalid step or rank/file mismatch")
        if row["dropped_loops"]:
            raise ValueError("diagnostic loop capacity exceeded; attribution is incomplete")
        if loop:
            if (row["kind"] not in (0, 1, 2) or row["invoked"] not in (0, 1) or
                    row["globalization_valid"] not in (0, 1) or row["attempt"] == 0 or
                    row.get("scalar_coupling_sweep", 1) == 0 or row.get("dt", 1) <= 0):
                raise ValueError("invalid loop discriminator")
            if row["globalization_valid"] and not set(FLOATS[3:7]).issubset(row):
                raise ValueError("missing physical residuals")
        yield row


def step_groups(rows, maximum):
    step, group = None, []
    for row in rows:
        current = row["step"]
        if step is not None and current != step:
            if current < step:
                raise ValueError("nonmonotonic step order")
            yield step, group
            group = []
        step = current
        group.append(row)
        if len(group) > maximum:
            raise ValueError("step exceeds frozen row capacity")
    if group:
        yield step, group


def ratio(final, initial):
    return final / initial if initial else None


def summarize_step(step, loop_rows, performance, expected_ranks):
    grouped, by_rank = defaultdict(list), defaultdict(list)
    for row in loop_rows:
        key = (row["attempt"], row.get("scalar_coupling_sweep", 1),
               row["corrector"], row["refinement"], row["kind"])
        grouped[key].append(row)
        by_rank[row["rank"]].append(row)
    if len(performance) != len(expected_ranks) or {r["rank"] for r in performance} != expected_ranks:
        raise ValueError("incomplete or duplicate performance rank coverage")
    for row in performance:
        phases = ("advance_ns", "observables_ns", "visit_ns",
                  "evidence_resources_ns", "csv_ns", "checkpoint_ns")
        if "full_step_ns" in row and (not set(phases).issubset(row) or
                sum(row[k] for k in phases) != row["full_step_ns"]):
            raise ValueError("full-step timing attribution mismatch")
        loops = by_rank[row["rank"]]
        for kind, name in enumerate(("pressure", "diagonal", "spatial")):
            if sum(v["invoked"] for v in loops if v["kind"] == kind) != row[name + "_calls"]:
                raise ValueError("solve-count attribution mismatch: " + name)
        for detailed, total in (("A_ns", "A_apply_ns"), ("M_ns", "M_apply_ns")):
            if sum(v[detailed] for v in loops) != row[total]:
                raise ValueError("linear timing attribution mismatch: " + detailed)
    summaries = []
    for key, group in sorted(grouped.items()):
        if len(group) != len(expected_ranks) or {r["rank"] for r in group} != expected_ranks:
            raise ValueError("incomplete or duplicate rank coverage at a logical loop")
        timings = [k for k in group[0] if k.endswith("_ns")]
        if any(set(row) != set(group[0]) for row in group):
            raise ValueError("inconsistent per-rank loop schema")
        summaries.append({"step": step, "attempt": key[0], "scalar_coupling_sweep": key[1],
            "corrector": key[2], "refinement": key[3], "kind": ("pressure", "diagonal", "spatial")[key[4]],
            "ranks": len(group), "counts_by_rank": {k: sorted({r[k] for r in group}) for k in COUNTS},
            "rank_mean_ns": {k: sum(r[k] for r in group) / len(group) for k in timings},
            "rank_max_ns": {k: max(r[k] for r in group) for k in timings},
            "linear_contraction": [ratio(r["linear_final"], r["linear_initial"]) for r in group],
            "physical_contraction": [{"rank": r["rank"],
                "continuity": ratio(r["selected_continuity"], r["baseline_continuity"]),
                "energy": ratio(r["selected_energy"], r["baseline_energy"])}
                for r in group if r["globalization_valid"]]})
    return summaries


def summarize(root, allow_partial=False, expected_ranks=None, details=None):
    result = {"schema": "HUNDUN_LOOP_OBSERVATION_V3", "complete": False,
              "formal_speed_comparison": False, "loops": [], "issues": [],
              "validated_steps": 0, "loop_count": 0}
    totals, samples = dict.fromkeys(ADDED, 0), 0
    paths = sorted(root.glob("solver-rank-*.csv"))
    try:
        ranks, first, last, bound, source = metadata(root, expected_ranks)
        result.update(expected_rank_count=ranks, expected_step_range=[first, last])
        if len(paths) != ranks:
            raise ValueError("per-rank files do not cover the frozen rank set")
        result["expected_ranks"] = list(range(ranks))
        if not bound:
            if not allow_partial:
                raise ValueError("legacy source binding cannot be called complete")
            result["issues"].append("legacy metadata/CSV lacks frozen coverage and source binding")
        expected = set(range(ranks))
        if {p.name for p in paths} != {"solver-rank-{}.csv".format(r) for r in expected}:
            raise ValueError("per-rank files do not cover the frozen rank set")
        paths += [root / "performance.csv", root / "conservation.csv", root / "RUN.meta"]
        stamps = {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in paths}
        with ExitStack() as stack:
            streams = [stack.enter_context((root / "solver-rank-{}.csv".format(r)).open()) for r in range(ranks)]
            streams.append(stack.enter_context((root / "performance.csv").open()))
            groups = [step_groups(checked_rows(s, source if bound else None, True, r), 64)
                      for r, s in enumerate(streams[:-1])]
            groups.append(step_groups(checked_rows(streams[-1], source if bound else None, False), ranks))
            for step in range(first, last + 1):
                batches = [next(g, None) for g in groups]
                if any(b is None or b[0] != step for b in batches):
                    raise ValueError("missing/truncated expected step {}".format(step))
                loops = [row for b in batches[:-1] for row in b[1]]
                performance = batches[-1][1]
                summaries = summarize_step(step, loops, performance, expected)
                for summary in summaries:
                    if details is None:
                        result["loops"].append(summary)
                    else:
                        details.write(json.dumps(summary, allow_nan=False) + "\n")
                result["loop_count"] += len(summaries)
                result["validated_steps"] += 1
                for row in performance:
                    for k in ADDED:
                        totals[k] += row.get(k, 0)
                samples += len(performance)
            if any(next(g, None) is not None for g in groups):
                raise ValueError("rows beyond frozen requested step range")
        result["files_sha256"] = {p.name: file_hash(p) for p in paths}
        if any((p.stat().st_size, p.stat().st_mtime_ns) != stamps[p] for p in paths):
            raise ValueError("source changed during inspection; freeze a diagnostic copy")
        result["complete"] = bound
    except (ValueError, KeyError, OSError, OverflowError, csv.Error) as error:
        if not allow_partial:
            raise ValueError(str(error)) from error
        result["issues"].append(str(error))
    result["rank_step_mean_ns"] = {k: totals[k] / samples for k in ADDED} if samples else {}
    result["scope"] = [
        "complete means the frozen requested diagnostic window, not scientific acceptance",
        "partial results contain only validated whole steps, never a complete attribution",
        "M/A include communication; MG copy is within refill and solve preparation",
        "Structured wait/control excludes IBM donor and MG-specific MPI",
        "Composition sweeps and dt retries are separate; old CSV defaults to sweep 1",
        "64-loop per-advance overflow remains a rejection; no discarded loops are hidden",
        "Input memory is bounded by one step per rank; hashes use 1 MiB blocks",
        "Use --details-output for bounded output memory; default loops JSON grows with the window",
        "Independent maxima are never summed into step time; no physical thresholds change"]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--expected-ranks", type=int, help="explicit legacy diagnostic expectation; never upgrades to complete")
    parser.add_argument("--details-output", type=Path, help="exclusive JSONL sidecar; bounded-memory loop details")
    args = parser.parse_args()
    created = []
    try:
        with ExitStack() as stack:
            def exclusive(path):
                if path is None:
                    return None
                stream = stack.enter_context(path.open("x"))
                created.append(path)
                return stream
            output, details = exclusive(args.output), exclusive(args.details_output)
            result = summarize(args.run_root, args.allow_partial, args.expected_ranks, details)
            if details:
                details.flush()
                result["details"] = {"path": str(args.details_output), "sha256": file_hash(args.details_output),
                                     "complete": result["complete"]}
            payload = json.dumps(result, indent=2, allow_nan=False) + "\n"
            if output:
                output.write(payload)
            else:
                print(payload, end="")
    except (ValueError, OSError) as error:
        # Remove only files exclusively created by this invocation. A failed
        # strict validation must not leave an apparently publishable sidecar.
        for path in created:
            try:
                path.unlink()
            except OSError:
                pass
        parser.exit(2, "observation rejected: {}\n".format(error))


if __name__ == "__main__":
    main()
