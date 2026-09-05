#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile

with tempfile.TemporaryDirectory(prefix="hundun-pmpi-") as root:
    environment = os.environ.copy()
    environment["LD_PRELOAD"] = sys.argv[2]
    environment["HUNDUN_MPI_OBSERVER_DIR"] = root
    subprocess.run([sys.argv[1], "-n", "2", sys.argv[3]], env=environment,
                   check=True, timeout=30)
    files = sorted(Path(root).glob("rank-*.tsv"))
    assert len(files) == 2
    for path in files:
        lines = path.read_text().splitlines()
        assert lines[-1] == "# dropped_calls=0"
        records = list(csv.DictReader(lines[:-1], delimiter="\t"))
        counts = {}
        for row in records:
            assert int(row["step"]) == 7
            assert int(row["site"], 16) > 0
            assert int(row["ns"]) >= 0 and row["module"]
            counts[row["kind"]] = counts.get(row["kind"], 0) + int(row["calls"])
        assert counts == {"Allreduce": 1, "Bcast": 2, "Barrier": 1}, counts
    print("scoped PMPI observer 2 ranks: exactly 4 in-scope calls per rank")
