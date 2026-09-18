#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Deterministic parser and conservative-accumulator gate for the AECSC bridge."""

import json
import subprocess
import sys


executable = sys.argv[1]
result = subprocess.run(
    [executable, "--self-test"], text=True, stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT, check=False)
assert result.returncode == 0, result.stdout
lines = [line for line in result.stdout.splitlines() if line.strip()]
assert len(lines) == 1, result.stdout
report = json.loads(lines[0])
assert report == {
    "status": "pass",
    "fortran_records": 3,
    "source_samples": 3,
    "target_fluid_cells": 3,
    "filled_target_cells": 1,
    "source_mass": 10.0,
    "target_mass": 10.0,
    "momentum_x": 26.0,
    "species_sum_error": 0.0,
}

bad = subprocess.run(
    [executable, "--unknown"], text=True, stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT, check=False)
assert bad.returncode != 0 and "usage:" in bad.stdout, bad.stdout
print("AECSC PDF transfer parser/conservation self-test passed")
