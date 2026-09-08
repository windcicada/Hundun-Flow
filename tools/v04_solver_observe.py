#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Validate frozen coverage, then attribute a diagnostic window (not a speed gate)."""
import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path

COUNTS = ("invoked", "iterations", "A_calls", "M_calls", "baseline_candidates",
          "extrapolated_candidates", "ladder_candidates", "incomplete_candidates")
ADDED = ("final_momentum_ns", "terminal_metrics_ns", "boundary_ledger_ns",
         "structured_wait_ns", "structured_control_ns", "scalar_remap_ns")
FLOATS = ("dt", "linear_initial", "linear_final", "baseline_continuity",
          "baseline_energy", "selected_continuity", "selected_energy", "selected_alpha")
CRITERION_FLOATS = ("linear_rhs_norm", "linear_atol", "linear_rtol", "linear_residual_limit")
CRITERION = ("linear_criterion_valid",) + CRITERION_FLOATS
SOURCE = "source_meta_sha256"
MG_MAX_LEVELS = 32
MG_PHASES = ("pre_smooth", "residual", "restriction", "prolongation", "post_smooth", "terminal", "direct_mpi")
MG_TOTALS = ("attempts", "successes", "failures", "apply_ns", "reduction_ns")
MG_LEVELS = tuple(p + s for p in MG_PHASES for s in ("_calls", "_ns")) + (
    "halo_wait_ns", "halo_control_ns", "halo_control_calls")
MG_LOOP = ("mg_enabled", "mg_complete") + tuple("mg_" + k for k in MG_TOTALS + MG_LEVELS) + (
    "mg_finest_pre_smooth_ns", "mg_finest_post_smooth_ns")


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
    schema = values.get("observation_schema")
    if schema not in (None, "3", "4", "5"):
        raise ValueError("unsupported observation schema")
    bound = schema in ("3", "4", "5")
    if schema == "5" and (values.get("observe_performance") != "1" or values.get("observe_mg_cost") != "1"):
        raise ValueError("MG observation flags missing from frozen metadata")
    if schema != "5" and values.get("observe_mg_cost", "0") != "0":
        raise ValueError("MG observation needs schema 5")
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
    return ranks, first, first + count - 1, bound, file_hash(path), schema, values.get("mg_layout_sha256")


def mg_layout(path, ranks, digest):
    # Bound the cold description before reading. This scales with the frozen
    # rank/level count, never with the number of observed steps or field cells.
    if path.stat().st_size > 128 + ranks * (64 + MG_MAX_LEVELS * 256):
        raise ValueError("oversized MG layout")
    lines = path.read_text().splitlines()
    if (lines[:3] != ["HUNDUN_MG_LAYOUT_V1", "expected_ranks {}".format(ranks),
                      "maximum_levels {}".format(MG_MAX_LEVELS)] or not lines or lines[-1] != "end"):
        raise ValueError("invalid MG layout envelope")
    if file_hash(path) != digest:
        raise ValueError("MG layout source identity mismatch")
    counts, geometry, cursor = {}, {}, 3
    for rank in range(ranks):
        if cursor >= len(lines) - 1:
            raise ValueError("missing MG layout rank")
        header = lines[cursor].split()
        if len(header) != 3 or header[:2] != ["rank", str(rank)]:
            raise ValueError("invalid MG layout rank order")
        count = int(header[2])
        if count <= 0 or count > MG_MAX_LEVELS:
            raise ValueError("invalid MG layout level count")
        counts[rank] = count
        cursor += 1
        for level in range(count):
            if cursor >= len(lines) - 1:
                raise ValueError("missing MG layout level")
            words = lines[cursor].split()
            if len(words) != 11 or words[:3] != ["level", str(rank), str(level)]:
                raise ValueError("invalid MG layout level order")
            shape = list(map(int, words[3:]))
            if (any(v <= 0 or v > 2**31 - 1 for v in shape[:6]) or
                    any(shape[c + 3] > shape[c] for c in range(3)) or
                    shape[6] not in range(7) or shape[7] not in range(8)):
                raise ValueError("invalid MG layout geometry")
            geometry[rank, level] = dict(rank=rank, level=level, global_shape=shape[:3],
                local_shape=shape[3:6], coarsening=shape[6], line_axis_mask=shape[7])
            cursor += 1
    if cursor != len(lines) - 1:
        raise ValueError("extra MG layout records")
    return counts, geometry


def checked_mg_rows(stream, source, rank):
    def lines():
        for line in stream:
            if not line.endswith("\n"):
                raise ValueError("truncated MG CSV tail")
            yield line
    reader = csv.DictReader(lines())
    expected = {"step", "rank", "level", "initialized", "complete", "visits", SOURCE}
    expected.update(MG_TOTALS + MG_LEVELS)
    if not reader.fieldnames or len(reader.fieldnames) != len(expected) or set(reader.fieldnames) != expected:
        raise ValueError("missing/duplicate MG CSV columns")
    for raw in reader:
        if None in raw or any(v is None or v == "" for v in raw.values()):
            raise ValueError("malformed MG CSV row")
        if raw.pop(SOURCE) != source:
            raise ValueError("MG source identity mismatch")
        row = {k: int(v) for k, v in raw.items()}
        if any(v < (-1 if k == "level" else 0) or v >= 2**64 - 1 for k, v in row.items()):
            raise ValueError("invalid/saturated MG counter")
        if row["step"] == 0 or row["rank"] != rank:
            raise ValueError("MG step or rank/file mismatch")
        if row["initialized"] != 1 or row["complete"] != 1:
            raise ValueError("MG observation uninitialized/incomplete")
        yield row


def validate_mg_step(loops, batches, counts):
    by_rank = defaultdict(list)
    for row in loops:
        by_rank[row["rank"]].append(row)
    for rank, (_, rows) in enumerate(batches):
        if len(rows) != counts[rank] + 1 or {r["level"] for r in rows} != set(range(-1, counts[rank])):
            raise ValueError("missing/duplicate MG level coverage")
        total = next(r for r in rows if r["level"] == -1)
        levels = [r for r in rows if r["level"] != -1]
        if total["attempts"] != total["successes"] + total["failures"]:
            raise ValueError("MG apply outcome mismatch")
        if total["visits"] or any(total[k] for k in MG_LEVELS) or any(r[k] for r in levels for k in MG_TOTALS):
            raise ValueError("MG total and level columns overlap")
        for key in MG_TOTALS:
            if sum(r["mg_" + key] for r in by_rank[rank]) != total[key]:
                raise ValueError("MG loop/step attribution mismatch: " + key)
        for key in MG_LEVELS:
            if sum(r["mg_" + key] for r in by_rank[rank]) != sum(r[key] for r in levels):
                raise ValueError("MG loop/level attribution mismatch: " + key)
        finest = next(r for r in levels if r["level"] == 0)
        for phase in ("pre_smooth", "post_smooth"):
            if sum(r["mg_finest_" + phase + "_ns"] for r in by_rank[rank]) != finest[phase + "_ns"]:
                raise ValueError("MG finest attribution mismatch")


def checked_rows(stream, source, loop, rank=None, require_criterion=False, require_mg=False):
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
    if loop and (require_criterion or set(CRITERION).intersection(header)) and not set(CRITERION).issubset(header):
        raise ValueError("missing/partial linear criterion columns: " + stream.name)
    if require_mg and not set(MG_LOOP).issubset(header):
        raise ValueError("missing MG loop columns")
    if not require_mg and set(MG_LOOP).intersection(header):
        raise ValueError("MG loop columns need frozen schema 5")
    for raw in reader:
        if None in raw or any(v is None or v == "" for v in raw.values()):
            raise ValueError("malformed CSV row: " + stream.name)
        row = {}
        for key, value in raw.items():
            if key == SOURCE:
                if source is not None and value != source:
                    raise ValueError("source identity mismatch: " + stream.name)
                continue
            floating = key in FLOATS or key in CRITERION_FLOATS
            number = float(value) if floating else int(value)
            if not math.isfinite(number) or number < 0 or (not floating and number > 2**64 - 1):
                raise ValueError("invalid count/time/residual: " + key)
            row[key] = number
        if row["step"] == 0 or (rank is not None and row["rank"] != rank):
            raise ValueError("invalid step or rank/file mismatch")
        if row["dropped_loops"]:
            raise ValueError("diagnostic loop capacity exceeded; attribution is incomplete")
        if loop:
            if require_mg:
                if row["mg_enabled"] != 1 or row["mg_complete"] != 1:
                    raise ValueError("MG loop observation disabled/incomplete")
                if any(row[k] == 2**64 - 1 for k in MG_LOOP):
                    raise ValueError("saturated MG loop counter")
                if (row["mg_attempts"] != row["M_calls"] or
                        row["mg_attempts"] != row["mg_successes"] + row["mg_failures"]):
                    raise ValueError("MG loop apply count mismatch")
                phase_ns = sum(row["mg_" + p + "_ns"] for p in MG_PHASES[:-1])
                if phase_ns > row["mg_apply_ns"] or row["mg_apply_ns"] > row["M_ns"]:
                    raise ValueError("MG phase/apply timing overlap")
            if (row["kind"] not in (0, 1, 2) or row["invoked"] not in (0, 1) or
                    row["globalization_valid"] not in (0, 1) or row["attempt"] == 0 or
                    row.get("scalar_coupling_sweep", 1) == 0 or row.get("dt", 1) <= 0):
                raise ValueError("invalid loop discriminator")
            if row["globalization_valid"] and not set(FLOATS[3:7]).issubset(row):
                raise ValueError("missing physical residuals")
            if "linear_criterion_valid" in row:
                if row["linear_criterion_valid"] not in (0, 1):
                    raise ValueError("invalid linear criterion flag")
                if not row["linear_criterion_valid"]:
                    if any(row[k] != 0 for k in CRITERION_FLOATS):
                        raise ValueError("unavailable linear criterion must have zero fields")
                else:
                    if row["linear_atol"] == 0 and row["linear_rtol"] == 0:
                        raise ValueError("linear criterion needs at least one positive tolerance")
                    expected = max(row["linear_atol"], row["linear_rtol"] * row["linear_rhs_norm"])
                    # Permit binary64 multiplication/decimal round-trip rounding,
                    # including subnormal results, while preserving exact zero;
                    # this is not a solver tolerance.
                    limit = row["linear_residual_limit"]
                    if ((limit == 0) != (expected == 0) or not math.isclose(limit, expected,
                            rel_tol=8 * sys.float_info.epsilon,
                            abs_tol=8 * sys.float_info.min * sys.float_info.epsilon)):
                        raise ValueError("linear criterion limit does not match max(atol, rtol * rhs_norm)")
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
        criterion = None
        if "linear_criterion_valid" in group[0]:
            if any(any(row[k] != group[0][k] for k in CRITERION) for row in group):
                raise ValueError("inconsistent per-rank linear criterion")
            if group[0]["linear_criterion_valid"]:
                criterion = {k[len("linear_"):]: group[0][k] for k in CRITERION_FLOATS}
        summaries.append({"step": step, "attempt": key[0], "scalar_coupling_sweep": key[1],
            "corrector": key[2], "refinement": key[3], "kind": ("pressure", "diagonal", "spatial")[key[4]],
            "ranks": len(group), "counts_by_rank": {k: sorted({r[k] for r in group}) for k in COUNTS},
            "rank_mean_ns": {k: sum(r[k] for r in group) / len(group) for k in timings},
            "rank_max_ns": {k: max(r[k] for r in group) for k in timings},
            "linear_criterion_available": criterion is not None,
            "linear_criterion": criterion,
            "linear_contraction": [ratio(r["linear_final"], r["linear_initial"]) for r in group],
            "physical_contraction": [{"rank": r["rank"],
                "continuity": ratio(r["selected_continuity"], r["baseline_continuity"]),
                "energy": ratio(r["selected_energy"], r["baseline_energy"])}
                for r in group if r["globalization_valid"]]})
        if "mg_enabled" in group[0]:
            summaries[-1]["mg_counts_by_rank"] = {
                k: sorted({r[k] for r in group}) for k in MG_LOOP
                if k.endswith("_calls") or k in ("mg_attempts", "mg_successes", "mg_failures")}
    return summaries


def summarize(root, allow_partial=False, expected_ranks=None, details=None):
    result = {"schema": "HUNDUN_LOOP_OBSERVATION_V3", "complete": False,
              "formal_speed_comparison": False, "loops": [], "issues": [],
              "validated_steps": 0, "loop_count": 0}
    totals, samples = dict.fromkeys(ADDED, 0), 0
    schema = None  # Metadata itself may be invalid in an explicit partial read.
    paths = sorted(root.glob("solver-rank-*.csv"))
    try:
        ranks, first, last, bound, source, schema, mg_digest = metadata(root, expected_ranks)
        if schema in ("4", "5"):
            result["schema"] = "HUNDUN_LOOP_OBSERVATION_V" + schema
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
        mg_paths = []
        if schema == "5":
            mg_paths = [root / "mg-rank-{}.csv".format(r) for r in range(ranks)]
            if set(root.glob("mg-rank-*.csv")) != set(mg_paths):
                raise ValueError("MG files do not cover the frozen rank set")
            paths += mg_paths + [root / "mg-layout.meta"]
        stamps = {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in paths}
        if schema == "5":
            mg_counts, geometry = mg_layout(root / "mg-layout.meta", ranks, mg_digest)
            mg_levels = {key: dict(value, visits=0, **dict.fromkeys(MG_LEVELS, 0))
                         for key, value in geometry.items()}
            mg_totals = [dict(rank=r, **dict.fromkeys(MG_TOTALS, 0)) for r in range(ranks)]
            result["mg"] = {"layout_sha256": mg_digest, "totals_by_rank": mg_totals,
                            "levels": [mg_levels[key] for key in sorted(mg_levels)]}
        with ExitStack() as stack:
            streams = [stack.enter_context((root / "solver-rank-{}.csv".format(r)).open()) for r in range(ranks)]
            streams.append(stack.enter_context((root / "performance.csv").open()))
            groups = [step_groups(checked_rows(s, source if bound else None, True, r,
                                              schema in ("4", "5"), schema == "5"), 64)
                      for r, s in enumerate(streams[:-1])]
            groups.append(step_groups(checked_rows(streams[-1], source if bound else None, False), ranks))
            mg_groups = [step_groups(checked_mg_rows(stack.enter_context(p.open()), source, r),
                                     mg_counts[r] + 1) for r, p in enumerate(mg_paths)]
            for step in range(first, last + 1):
                batches = [next(g, None) for g in groups]
                if any(b is None or b[0] != step for b in batches):
                    raise ValueError("missing/truncated expected step {}".format(step))
                loops = [row for b in batches[:-1] for row in b[1]]
                performance = batches[-1][1]
                summaries = summarize_step(step, loops, performance, expected)
                mg_batches = [next(g, None) for g in mg_groups]
                if schema == "5":
                    if any(b is None or b[0] != step for b in mg_batches):
                        raise ValueError("missing/truncated MG step {}".format(step))
                    validate_mg_step(loops, mg_batches, mg_counts)
                    # Commit diagnostics only after the whole step reconciles.
                    for r, (_, rows) in enumerate(mg_batches):
                        for row in rows:
                            if row["level"] == -1:
                                for k in MG_TOTALS:
                                    mg_totals[r][k] += row[k]
                            else:
                                for k in MG_LEVELS + ("visits",):
                                    mg_levels[r, row["level"]][k] += row[k]
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
            if any(next(g, None) is not None for g in groups + mg_groups):
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
        "linear_criterion reports max(atol, rtol * rhs_norm); availability and contraction are not acceptance gates",
        "partial results contain only validated whole steps, never a complete attribution",
        "M/A include communication; MG copy is within refill and solve preparation",
        "Structured wait/control excludes IBM donor and MG-specific MPI",
        "Composition sweeps and dt retries are separate; old CSV defaults to sweep 1",
        "64-loop per-advance overflow remains a rejection; no discarded loops are hidden",
        "Input memory is bounded by one step per rank; hashes use 1 MiB blocks",
        "Use --details-output for bounded output memory; default loops JSON grows with the window",
        "Independent maxima are never summed into step time; no physical thresholds change"]
    if schema == "5":
        result["scope"] += [
            "MG layout is cold ProductDriver authority; step level rows reconcile to compact loop totals",
            "MG six phases are disjoint subsets of apply; halo, reductions and direct MPI are nested",
            "MG totals and levels sum only validated steps; per-loop full level distribution is not retained",
            "MG observations exclude Fresh projection and refill/copy; no solver status is inferred from completeness"]
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
