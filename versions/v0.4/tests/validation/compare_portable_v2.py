#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Compare already-run 1/2/4-rank configurations; never launches simulations.

FNV fingerprints are reproducibility observations, not scientific reference
answers. Physical correctness is checked separately against analytic budgets.
"""
import re
import sys
import xml.etree.ElementTree as ET


def compare(path):
    groups = {}
    for test in ET.parse(path).iter("testcase"):
        name = test.get("name", "")
        match = re.fullmatch(
            r"v04_portable_v2_(spray_(?:alpha|beta)|closed_analytic)_n([24])_r([124])",
            name,
        )
        if not match:
            continue
        if test.find("failure") is not None or test.find("skipped") is not None:
            raise ValueError("non-passing configuration: " + name)
        records = dict(re.findall(
            r"^(partition_state_accepted[12])=(\d+)$",
            test.findtext("system-out", ""), re.MULTILINE))
        required = 1 if match.group(1) == "closed_analytic" else 2
        if len(records) != required:
            raise ValueError("missing state observations: " + name)
        group = groups.setdefault(match.group(1, 2), {})
        if match.group(3) in group:
            raise ValueError("duplicate configuration: " + name)
        group[match.group(3)] = records
    if len(groups) != 6:
        raise ValueError("expected six independent pack/field groups")
    for group, ranks in sorted(groups.items()):
        if set(ranks) != {"1", "2", "4"}:
            raise ValueError("missing rank configuration: " + repr(group))
        if not ranks["1"] == ranks["2"] == ranks["4"]:
            raise ValueError("partition-dependent accepted state: " + repr(group))
    print("PASS: six groups have identical accepted-state fingerprints on 1/2/4 ranks")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: compare_portable_v2.py ctest-junit.xml")
    try:
        compare(sys.argv[1])
    except (ValueError, OSError, ET.ParseError) as error:
        raise SystemExit(str(error))
