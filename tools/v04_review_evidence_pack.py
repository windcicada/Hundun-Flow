#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Package an existing review's logs and identities; never run or certify tests."""

import argparse
import gzip
import hashlib
import json
import pathlib
import re
import subprocess
import tempfile


def digest(data):
    return hashlib.sha256(data).hexdigest()


def run(repo, *args):
    return subprocess.check_output(args, cwd=str(repo))


def ctest_summary(text):
    summaries = re.findall(r"(\d+)% tests passed, (\d+) tests failed out of (\d+)", text)
    if not summaries:
        return {"complete_summary_found": False}
    _, failed, total = map(int, summaries[-1])
    return {
        "complete_summary_found": True,
        "total": total,
        "passed": total - failed,
        "failed": failed,
        "failed_tests": re.findall(r"^\s*\d+ - (\S+) \(([^)]+)\)", text, re.M),
    }


def trace_records(text):
    result = {"memory": [], "capacity": [], "collectives": []}
    for prefix, key in (("MEMORY", "memory"), ("MEMORY_CAPACITY", "capacity"),
                        ("COLLECTIVES", "collectives")):
        for match in re.finditer(r"\b" + prefix + r" (\{[^\n]+\})", text):
            result[key].append(json.loads(match.group(1)))
    return result


def write_payload(root, relative, data, compressed=True):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as output:
        if compressed:
            with gzip.GzipFile(filename="", fileobj=output, mode="wb", mtime=0) as stream:
                stream.write(data)
        else:
            output.write(data)
    return {"path": relative, "raw_bytes": len(data), "raw_sha256": digest(data),
            "stored_bytes": path.stat().st_size, "stored_sha256": digest(path.read_bytes())}


def historical_receipt_audit(repo):
    receipt = repo / "docs/verification/v0.4-literature-data-receipt-r4-partial.json"
    result = []
    for item in json.loads(receipt.read_text())["artifacts"]:
        original = pathlib.Path(item["path"])
        candidate = original
        prefix = "/home/wyf/code_dev/hundun-flow/"
        if not original.is_file() and str(original).startswith(prefix):
            candidate = repo / str(original)[len(prefix):]
        found = candidate.is_file()
        actual = digest(candidate.read_bytes()) if found else None
        result.append({"original_path": str(original), "checked_path": str(candidate),
                       "relocated_candidate": candidate != original, "exists": found,
                       "expected_sha256": item["sha256"], "actual_sha256": actual,
                       "hash_matches": actual == item["sha256"]})
    return result


def package(args):
    repo = args.repo.resolve()
    root = args.output.resolve()
    logs = sorted(args.logs.glob(args.log_prefix + "*.log"))
    if not logs:
        raise ValueError("No review logs found; refusing to create empty evidence")
    root.mkdir(parents=True, exist_ok=False)
    manifest = []
    summaries = {}
    for path in logs:
        data = path.read_bytes()
        entry = write_payload(root, "logs/" + path.name + ".gz", data)
        entry["original_path"] = str(path)
        manifest.append(entry)
        parsed = ctest_summary(data.decode("utf-8", errors="replace"))
        if parsed["complete_summary_found"]:
            summaries[path.name] = parsed
    source_paths = ["versions/v0.4", "tools/v04_thin_domain_runner.cpp",
                    "tools/v04_review_evidence_pack.py", "README.md"]
    source_paths.extend(args.source_path)
    patch = run(repo, "git", "diff", "--binary", args.baseline, "--", *source_paths)
    manifest.append(write_payload(root, "implementation.patch.gz", patch))
    builds = {}
    binaries = {}
    for name in ("build-re3900-fix-test", "build-audit-sanitize"):
        build = repo / name
        cache = build / "CMakeCache.txt"
        if cache.is_file():
            manifest.append(write_payload(root, "builds/" + name + "-CMakeCache.txt.gz", cache.read_bytes()))
            selected = {}
            for line in cache.read_text().splitlines():
                if re.match(r"(?:CMAKE_(?:BUILD_TYPE|CXX_COMPILER|C_COMPILER|CXX_FLAGS|C_FLAGS|EXE_LINKER_FLAGS)|HUNDUN_|MPI_)[^:]*:", line):
                    key, _, value = line.partition("=")
                    selected[key] = value
            builds[name] = selected
        for relative in ("hundun", "v04_thin_domain_runner",
                         "tests/libv04_mpi_observer.so", "tests/v04_mpi_observer_probe",
                         "tests/v04_product_restart_storage_mpi_test",
                         "tests/v04_app_output_failure_mpi_test",
                         "tests/v04_io_product_path_test",
                         "tests/v04_product_allocation_failure_mpi_test",
                         "tests/v04_product_memory_profile_mpi_test",
                         "tests/v04_product_mg_capacity_mpi_test",
                         "tests/v04_product_pressure_energy_retry_mpi_test",
                         "tests/v04_parallel_halo_mpi_test"):
            path = build / "versions/v0.4" / relative
            if path.is_file():
                binaries[str(path.relative_to(repo))] = {"bytes": path.stat().st_size,
                                                         "sha256": digest(path.read_bytes())}
    trace = args.logs / (args.log_prefix + "final-release-full-details.log")
    records = trace_records(trace.read_text(errors="replace")) if trace.is_file() else {}
    command_statuses = args.logs / (args.log_prefix + "command-statuses.log")
    if args.run_root:
        for name in ("RUN.meta", "evidence.jsonl", "health.csv", "force.csv", "probe.csv",
                     "performance.csv", "Restart/current"):
            path = args.run_root / name
            if path.is_file():
                manifest.append(write_payload(root, "diagnostic/" + name + ".gz", path.read_bytes()))
        for path in sorted(args.run_root.glob("step-*.complete")):
            manifest.append(write_payload(root, "diagnostic/" + path.name + ".gz", path.read_bytes()))
    if args.mpi_directory:
        for path in sorted(args.mpi_directory.glob("rank-*.tsv")):
            manifest.append(write_payload(root, "diagnostic/mpi/" + path.name + ".gz", path.read_bytes()))
    observation = json.loads(args.observation.read_text()) if args.observation else None
    results = {
        "schema": "HUNDUN_MAIN_REVIEW_EVIDENCE_V1",
        "baseline_sha": args.baseline,
        "head_when_packed_not_a_claim_of_build_head": run(repo, "git", "rev-parse", "HEAD").decode().strip(),
        "branch": run(repo, "git", "branch", "--show-current").decode().strip(),
        "implementation_patch_sha256": digest(patch),
        "source_patch_paths": source_paths,
        "source_identity_note": "Baseline plus patch records packaged source, including postprocessing updates; untracked source must be staged before packing. Use binary hashes and RUN.meta for the diagnostic build identity. A later documentation/evidence commit does not relabel existing binaries.",
        "worktree_status_when_packed": run(repo, "git", "status", "--short").decode().splitlines(),
        "builds": builds, "binaries": binaries, "ctest_summaries": summaries,
        "recorded_command_statuses": (json.loads(command_statuses.read_text())
                                      if command_statuses.is_file() else None),
        "final_runtime_records": records,
        "historical_receipt_read_only_audit": historical_receipt_audit(repo),
        "historical_receipt_audit_scope": "Current-path comparison only, not receipt-validate --historical; archived bytes and current tools intentionally differ",
        "diagnostic_observation": observation,
        "formal_performance_comparison": {"performed": False, "paths": ["Hundun PISO", "Hundun SIMPLE", "COAST ordinary compressible"]},
        "limits": ["C++ new tracker excludes MPI/libc malloc, stacks and allocator overhead",
                   "ru_maxrss is a historical process high-water mark, not synchronized node peak",
                   "PMPI subtotal gaps are explicit and not assigned to unrelated counters",
                   "Migration and short-run outcomes must be read from their bound logs, not inferred by the packer",
                   "No long Re3900 run, push, tag or release"],
    }
    (root / "manifest.json").write_text(json.dumps({"schema": "HUNDUN_REVIEW_LOG_MANIFEST_V1", "artifacts": manifest}, indent=2) + "\n")
    (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    checksums = []
    for path in sorted(root.rglob("*")):
        if path.is_file():
            checksums.append(digest(path.read_bytes()) + "  " + str(path.relative_to(root)))
    (root / "SHA256SUMS").write_text("\n".join(checksums) + "\n")
    print(json.dumps({"output": str(root), "logs": len(logs), "artifacts": len(checksums),
                      "ctest_summaries": summaries}, indent=2))


def self_test():
    summary = ctest_summary("99% tests passed, 1 tests failed out of 247\n 7 - old_receipt (Failed)\n")
    assert summary["passed"] == 246 and summary["failed_tests"] == [("old_receipt", "Failed")]
    assert not ctest_summary("still running")["complete_summary_found"]
    assert trace_records('2: MEMORY {"rank":0}\n')["memory"] == [{"rank": 0}]
    with tempfile.TemporaryDirectory() as directory:
        root = pathlib.Path(directory)
        record = write_payload(root, "test.log.gz", b"original\n")
        assert gzip.decompress((root / record["path"]).read_bytes()) == b"original\n"
        assert record["raw_sha256"] == digest(b"original\n")
        try:
            write_payload(root, "test.log.gz", b"overwrite")
        except FileExistsError:
            pass
        else:
            raise AssertionError("must not overwrite evidence")
    print("review evidence pack self-test PASS")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path(__file__).resolve().parents[1])
    parser.add_argument("--logs", type=pathlib.Path, default=pathlib.Path("/tmp"))
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--baseline", default="69d8eee715c9565ccc033ea7fb95d2bf3b95aae6")
    parser.add_argument("--log-prefix", default="hundun-review-")
    parser.add_argument("--source-path", action="append", default=[])
    parser.add_argument("--run-root", type=pathlib.Path)
    parser.add_argument("--mpi-directory", type=pathlib.Path)
    parser.add_argument("--observation", type=pathlib.Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    elif args.output is None:
        parser.error("--output is required for packaging")
    else:
        package(args)


if __name__ == "__main__":
    main()
