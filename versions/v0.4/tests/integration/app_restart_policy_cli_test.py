#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Ordinary CLI: explicit recovery policy, exact continuation and read-only source."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--mpi", required=True)
    parser.add_argument("--ranks", type=int, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="hundun-app-history-") as directory:
        root = Path(directory)
        case = root / "case"
        binary = str(args.binary.resolve())
        result = subprocess.run([binary, "init-case", "--output", str(case)],
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                universal_newlines=True, timeout=20)
        assert result.returncode == 0, result.stdout

        def run(name, source=None, method=False, extra=()):
            path = root / name
            command = [args.mpi, "-n", str(args.ranks), binary, "run", str(case),
                       "--output", str(path), "--steps", "2", "--output-interval", "0"]
            if source:
                command += ["--restart", str(source / "Restart")]
            if method:
                command += ["--restart-method-recovery"]
            result = subprocess.run(command + list(extra), stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, universal_newlines=True, timeout=35)
            return path, result

        def accepted(name, source=None, method=False):
            path, result = run(name, source, method)
            assert result.returncode == 0, result.stdout
            rows = [json.loads(line) for line in (path / "evidence.jsonl").read_text().splitlines()]
            assert len(rows) == 2
            if source:
                assert [r["bdf_order"] for r in rows] == ([1, 2] if method else [2, 2])
                assert [r["restart_recovery"] for r in rows] == [method, False]
                assert all(r["run_start"]["history"]["policy"] ==
                           ("rebuild_method_history" if method else "require_compatible") for r in rows)
            return path

        source = accepted("source")
        exact = accepted("exact", source)
        hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in exact.rglob("*") if p.is_file()}
        recovered = accepted("recovered", exact, True)
        accepted("exact-again", recovered)
        assert all(hashlib.sha256(p.read_bytes()).hexdigest() == digest for p, digest in hashes.items())
        initial = ("--initial-state", "101325,300,0.1,0,0")
        path, result = run("explicit-initial", extra=initial)
        assert result.returncode == 0, result.stdout
        for name, origin, extra in (
                ("initial-with-restart", source, initial),
                ("initial-duplicate", None, initial + initial),
                ("initial-short", None, ("--initial-state", "101325,300,0,0")),
                ("initial-nan", None, ("--initial-state", "101325,nan,0,0,0")),
                ("initial-trailing", None, ("--initial-state", "101325,300,0,0,0,")),
                ("initial-extra-scalar", None, ("--initial-state", "101325,300,0,0,0,0.2"))):
            path, result = run(name, origin, extra=extra)
            assert result.returncode != 0 and not path.exists(), result.stdout
        for name, origin, extra in (("missing", None, ()),
                                     ("duplicate", source, ("--restart-method-recovery",)),
                                     ("incompatible-policies", source,
                                      ("--restart-storage-compatibility", "mg-bundle-ghost-v1"))):
            path, result = run(name, origin, True, extra)
            assert result.returncode != 0 and not path.exists(), result.stdout
        print("PASS ordinary CLI history chain ranks={}".format(args.ranks))


if __name__ == "__main__":
    main()
