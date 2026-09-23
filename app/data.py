"""Bounded, read-only views of native run files."""

import csv
import hashlib
import io
import json
import math
import os
import stat
from functools import lru_cache
from pathlib import Path

LIMIT = 4 * 1024 * 1024
SKIP = {"src", "build", "Restart", "Visit", ".git", "node_modules", "lib", "bin"}


def ident(path):
    return hashlib.sha256(str(Path(path).resolve()).encode()).hexdigest()[:20]


def read_json(path, default=None):
    try:
        info = path.lstat()
        if not stat.S_ISREG(info.st_mode) or info.st_size > LIMIT:
            return default
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def tail(path, limit=LIMIT):
    try:
        if not stat.S_ISREG(path.lstat().st_mode):
            return "", False
        with path.open("rb") as stream:
            stream.seek(0, 2)
            size = stream.tell()
            stream.seek(max(0, size - limit))
            data = stream.read(limit)
        if size > limit:
            data = data.partition(b"\n")[2]
        return data.decode("utf-8", "replace"), size > limit
    except OSError:
        return "", False


def rows(path, limit=LIMIT):
    try:
        stat = path.stat()
        reader = _rows if limit <= 65536 else _rows.__wrapped__
        return reader(str(path), stat.st_mtime_ns, stat.st_size, stat.st_ino, limit)
    except OSError:
        return [], False


@lru_cache(maxsize=192)
def _rows(raw_path, mtime, size, inode, limit):
    path = Path(raw_path)
    text, truncated = tail(path, limit)
    result = []
    # Only complete lines are published; an in-progress tail is not a record.
    for line in text.splitlines(keepends=True):
        if not line.endswith("\n"):
            continue
        try:
            value = json.loads(line)
            if isinstance(value, dict):
                result.append(value)
        except ValueError:
            pass
    return result, truncated


def safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [safe(v) for v in value]
    return value


def is_run(path):
    return any(
        (path / name).is_file() and not (path / name).is_symlink()
        for name in (
            "status.json",
            "monitor.jsonl",
            "evidence.jsonl",
            "diagnostics.jsonl",
            "force.csv",
        )
    ) or has_fields(path)


def has_fields(path):
    """Recognize published field metadata without reading any mesh payload."""
    root = (
        path / "Visit"
        if (path / "Visit").is_dir() and not (path / "Visit").is_symlink()
        else path
    )
    try:
        return _has_fields(str(root), root.stat().st_mtime_ns)
    except OSError:
        return False


@lru_cache(maxsize=192)
def _has_fields(raw_root, changed):
    root = Path(raw_root)
    try:
        for entry in root.iterdir():
            if (
                entry.suffix.lower() not in (".visit", ".pvtu", ".pvd", ".vtm")
                or entry.is_symlink()
                or not entry.is_file()
            ):
                continue
            if entry.stat().st_size > 65536:
                continue
            text = entry.read_text(errors="replace", encoding="utf-8")
            if entry.suffix.lower() == ".visit":
                lines = [line.strip() for line in text.splitlines() if line.strip()]
                if not any(line.startswith("!NBLOCKS ") for line in lines):
                    continue
                for line in lines:
                    if line.startswith("!"):
                        continue
                    target = (root / line).resolve()
                    if (
                        root.resolve() in target.parents
                        and target.suffix.lower() in (".vtk", ".vtu")
                        and target.is_file()
                    ):
                        return True
            elif "<VTKFile" in text and ("Source=" in text or "file=" in text):
                return True
    except OSError:
        pass
    return False


def walk(roots, depth=3):
    seen = set()
    for root in roots:
        root = root.resolve()
        if not root.is_dir():
            continue
        for current, dirs, _ in os.walk(str(root), followlinks=False):
            path = Path(current)
            level = len(path.relative_to(root).parts)
            dirs[:] = (
                sorted(
                    d
                    for d in dirs
                    if d not in SKIP
                    and not d.startswith(".")
                    and not (path / d).is_symlink()
                )
                if level < depth
                else []
            )
            if path not in seen:
                seen.add(path)
                yield path


def history(path, limit=LIMIT):
    monitor, truncated = rows(path / "monitor.jsonl", limit)
    evidence, ev_truncated = rows(path / "evidence.jsonl", limit)
    indexed = {}
    for row in evidence:
        audit = row.get("terminal_physical_audit", {})
        cold = row.get("cold", {})
        cfl = audit.get("committed_convective_cfl", {})
        indexed[row.get("step")] = {
            "step": row.get("step"),
            "time": row.get("time"),
            "seconds": row["max_rank_step_ns"] / 1e9
            if isinstance(row.get("max_rank_step_ns"), (int, float))
            else None,
            "dt": cfl.get("dt"),
            "continuity": audit.get("continuity_residual"),
            "energy": audit.get("energy_residual"),
            "eos": audit.get("eos_residual"),
            "cfl": cfl.get("directional_max", cfl.get("out_max")),
            "pressure_iterations": cold.get("pressure_iterations"),
            "outer": cold.get("outer_iterations"),
        }
    for row in monitor:
        p = row.get("payload", row)
        item = {
            "step": row.get("step"),
            "time": row.get("time"),
            "seconds": p.get("seconds"),
            "dt": p.get("dt"),
            "continuity": p.get("continuity"),
            "energy": p.get("energy"),
            "eos": p.get("eos"),
            "cfl": p.get(
                "committed_convective_cfl_directional",
                p.get("committed_convective_cfl_out"),
            ),
            "pressure_iterations": p.get("pressure_iterations"),
            "outer": p.get("outer_iterations"),
        }
        indexed.setdefault(row.get("step"), item).update(
            {k: v for k, v in item.items() if v is not None}
        )
    values = [
        indexed[k] for k in sorted(k for k in indexed if isinstance(k, (int, float)))
    ]
    try:
        if not stat.S_ISREG((path / "force.csv").lstat().st_mode):
            raise OSError("not a regular force table")
        with (path / "force.csv").open() as stream:
            header = stream.readline(16384)
        text, cut = tail(path / "force.csv", limit)
        if cut:
            text = header + text
        force = {
            int(r["step"]): r
            for r in csv.DictReader(io.StringIO(text))
            if r.get("step", "").isdigit()
            and r.get("included", "1") not in ("0", "false")
        }
        for row in values:
            for key in ("cd", "cl"):
                if row["step"] in force and force[row["step"]].get(key):
                    row[key] = float(force[row["step"]][key])
    except (OSError, ValueError, KeyError):
        pass
    count = len(values)
    stride = max(1, math.ceil(count / 600))
    selected = values[::stride]
    if values and (not selected or selected[-1] != values[-1]):
        selected.append(values[-1])
    return (
        safe(selected),
        {
            "sampled": stride > 1 or truncated or ev_truncated,
            "available_rows": count,
            "returned_rows": len(selected),
            "tail_limited": truncated or ev_truncated,
            "timing_scope": "max_rank_advance; excludes output",
            "first_step": values[0]["step"] if values else None,
            "last_step": values[-1]["step"] if values else None,
        },
        values,
        evidence,
    )


def detail(path, summary, case):
    hist, meta, _, evidence = history(path)
    diagnostic, _ = rows(path / "diagnostics.jsonl")
    controls, _ = rows(path / "control.jsonl")
    ev = evidence[-1] if evidence else {}
    resources = {
        k: ev[k]
        for k in (
            "max_rank_rss_bytes",
            "max_node_rss_bytes",
            "structured_bytes",
            "ibm_bytes",
            "linear_iterations",
            "heap_allocations",
        )
        if k in ev
    }
    if resources:
        resources["memory_scope"] = "peak_rss"
    logs = []
    for name in ("stderr.log", "stdout.log", "run.log", "screen.log"):
        text, _ = tail(path / name, 32000)
        logs.extend(text.splitlines()[-100:])
    files = []
    try:
        for child in sorted(path.iterdir())[:100]:
            if child.is_file() and not child.is_symlink():
                files.append({"name": child.name, "size": child.stat().st_size})
    except OSError:
        pass
    return safe(
        {
            "run": summary,
            "history": hist,
            "history_meta": meta,
            "budgets": diagnostic[-1].get("payload", {}) if diagnostic else {},
            "resources": resources or None,
            "logs": logs[-200:],
            "controls": controls[-50:],
            "files": files,
            "case": case,
        }
    )
