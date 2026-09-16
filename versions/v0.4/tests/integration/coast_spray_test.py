#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exercise legacy slot ownership, SI units, REAL daughter cells and failures."""
import importlib.util
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

tool = Path(sys.argv[1]).resolve()
spec = importlib.util.spec_from_file_location("spray_inventory", str(tool))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture():
    def part(d, n):
        return [1., 2., 3., d, n, 400., d, 400., .01, 0., 2., 3., 4., .02, 6., 7., 8., 9., 10., 11., 12., 13.]
    points = [1] + [2]*13 + [4]*15
    rows = [bytes(27*4)]*5 + [struct.pack("<2i", 4, 25500),
            struct.pack("<4i", 99, 99, 1, 2), struct.pack("<f", 1e-10),
            struct.pack("<10f", .1, .2, .3, 1., 2., 3., 10., 2., 400., 14.),
            struct.pack("<10f", .1, .2, .3, 1., 2., 3., 0., 2., 400., 14.),
            struct.pack("<29i", *points),
            struct.pack("<88f", *(part(1000., 1.)+part(20., 2.)+part(0., 2.)+part(30., 1.))),
            struct.pack("<4i", 4, 2, 3, 1), struct.pack("<4i", -1, 14, 14, 1),
            struct.pack("<12f", *([.1, .2, .3]*4)), struct.pack("<i", 1)]
    return rows


def pack(rows):
    return b"".join(struct.pack("<i", len(row))+row+struct.pack("<i", len(row)) for row in rows)


rows = fixture()
report, parcels = module.decode(pack(rows), (3, 3, 3), 0)
assert report["counts"] == dict(slots=4, retired=1, indexed=3, positive=2,
                               zero_inventory=1, interior_positive=1, halo_positive=1,
                               dummy_positive=0, pending_children=2, positive_children=1,
                               zero_inventory_children=1)
assert [(p["kind"], p["source_ordinal"]) for p in parcels] == [("retained", 4), ("retained", 2), ("daughter", 1)]
assert all(p["coast_cell"] == 14 for p in parcels[1:])
assert math.isclose(parcels[0]["droplet_diameter_m"], 30e-6, rel_tol=1e-15)
assert parcels[1]["multiplicity"] == 2 and parcels[2]["multiplicity"] == 2
assert parcels[1]["legacy_state"] == list(struct.unpack("<22f", rows[11][88:176]))
assert parcels[2]["legacy_state"] == list(struct.unpack("<10f", rows[8]))
expected_volume = math.pi/6*((30e-6)**3+2*(20e-6)**3+2*(10e-6)**3)
assert math.isclose(sum(math.pi/6*p["droplet_diameter_m"]**3*p["multiplicity"] for p in parcels), expected_volume, rel_tol=1e-15)

failures = [pack(rows)[:-1]]
for index, replacement in ((12, struct.pack("<4i", 4, 1, 3, 2)),
                           (12, struct.pack("<4i", 4, 2, 2, 1)),
                           (13, struct.pack("<4i", -1, 13, 14, 1)),
                           (8, struct.pack("<9fi", .1, .2, .3, 1., 2., 3., 10., 2., 400., 14)),
                           (10, struct.pack("<29i", *([1]+[2]*13+[5]*15))),
                           (6, struct.pack("<4i", 99, 99, 5, 2))):
    changed = list(rows)
    changed[index] = replacement
    failures.append(pack(changed))
for data in failures:
    try:
        module.decode(data, (3, 3, 3), 0)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid input admitted")

with tempfile.TemporaryDirectory(prefix="hf-spray-") as temp:
    root = Path(temp)
    source = root/"in"
    source.mkdir()
    (source/"Restart_spray_PDF.000").write_bytes(pack(rows))
    def run(out, ranks=1):
        return subprocess.run([sys.executable, str(tool), str(source), "--shape", "3", "3", "3",
                               "--ranks", str(ranks), "--output", str(out)],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert run(root/"out").returncode == 0
    inventory = json.loads((root/"out/inventory.json").read_text())
    assert inventory["counts"] == report["counts"]
    assert math.isclose(inventory["liquid_volume_m3"], expected_volume, rel_tol=1e-15)
    saved = (root/"out/parcels.jsonl").read_bytes()
    assert len(saved.splitlines()) == 3
    assert run(root/"out").returncode != 0
    assert (root/"out/parcels.jsonl").read_bytes() == saved
    assert run(root/"missing", 2).returncode != 0
    (source/"Restart_spray_PDF.001").write_bytes(pack(rows)[:-1])
    assert run(root/"bad", 2).returncode != 0
    assert not (root/"bad/parcels.jsonl").exists() and not (root/"bad/parcels.tmp").exists()
    busy = root/"busy"
    busy.mkdir()
    (busy/"parcels.tmp").write_text("existing work")
    assert run(busy).returncode != 0
    assert (busy/"parcels.tmp").read_text() == "existing work"
print("coast_spray slots=partitioned units=SI daughters=REAL histories=exact corruptions=7 publication=checked passed=1")
