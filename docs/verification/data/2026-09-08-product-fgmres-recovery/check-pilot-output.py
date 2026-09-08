#!/usr/bin/env python3
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Read-only, post-completion checks; does not replace runtime/MPI validation."""
import csv
import hashlib
import json
import math
from pathlib import Path
import struct
import sys
import xml.etree.ElementTree as ET

root = Path(sys.argv[1])
output = Path(sys.argv[2])
assert not output.exists(), "do not overwrite evidence"
expected_steps = list(range(9501, 9511))
meta = dict(line.split(" ", 1) for line in (root / "RUN.meta").read_text().splitlines() if " " in line)
assert meta["starting_step"] == "9500" and meta["requested_steps"] == "10"
assert meta["expected_ranks"] == "128"
assert meta["restart_method_recovery"] == "0"
assert meta["statistics_epoch_start_step"] == "7000"
assert meta["statistics_sampling_start_step"] == "17001"
assert meta["statistics_reset_reason"] == "method_recovery"
assert meta["starting_sample_steps"] == "0"
assert "COMPLETED" in Path(str(root) + ".log").read_text()
assert meta["restart_requires_recovery"] == "0"
assert meta["observation_schema"] == "6"
assert meta["observe_fgmres_recovery"] == "1" and meta["pressure_linear_algorithm"] == "fgmres"
completion = (root / "step-00000000000000009510.complete").read_text().splitlines()
assert completion[0] == "HUNDUN_V04_THIN_DOMAIN_CHECKPOINT_V1" and completion[-1] == "end"
done = dict(line.split(" ", 1) for line in completion[1:-1])
assert done["step"] == "9510"
assert done["restart_generation"] == (root / "Restart" / "current").read_text().strip()
generation = root / "Restart" / done["restart_generation"]
rank_files = sorted(generation.glob("rank-*.bin"))
assert len(rank_files) == 128 and all(p.stat().st_size > 0 for p in rank_files)
assert (generation / "manifest.bin").stat().st_size > 0
assert done["statistics"] == "step-00000000000000009510.statistics.json"
assert done["accumulator"] == "step-00000000000000009510.accumulator"
stats = json.loads((root / done["statistics"]).read_text())
assert stats["snapshot_step"] == 9510 and stats["sample_steps"] == 0
assert stats["statistics_epoch"] == {"start_step": 7000, "development_steps": 10000,
    "sampling_start_step": 17001, "reset_reason": "method_recovery"}
assert (root / done["accumulator"]).stat().st_size > 0

def read_csv(path):
    with path.open(newline="") as source:
        for row in csv.DictReader(source):
            assert None not in row and all(v is not None for v in row.values()), path
            yield row

health = list(read_csv(root / "health.csv"))
assert [int(r["step"]) for r in health] == expected_steps
metrics = {k: 0.0 for k in ("eos", "continuity", "energy", "closed_mass", "gauge")}
fluid = {k: [math.inf, -math.inf] for k in ("T", "p", "rho")}
solid_reference = None
for row in health:
    step = int(row["step"])
    assert int(row["attempts"]) == 1 and int(row["retry"]) == 0
    assert int(row["bdf_order"]) == 2
    assert int(row["restart_recovery"]) == 0
    assert int(row["terminal_audit_present"]) == 1
    assert int(row["fluid_count"]) == 5994352
    assert int(row["solid_placeholder_count"]) == 75920
    assert int(row["included"]) == 0
    assert int(row["temporal_fallback"]) == 0
    for key in metrics:
        value, tolerance = float(row[key]), float(row[key + "_tolerance"])
        assert math.isfinite(value) and 0.0 <= value <= tolerance
        metrics[key] = max(metrics[key], value)
    assert float(row["committed_cfl_out"]) <= float(row["committed_cfl_limit"])
    for key in fluid:
        lo, hi = float(row["fluid_" + key + "_min"]), float(row["fluid_" + key + "_max"])
        assert math.isfinite(lo) and math.isfinite(hi) and lo <= hi
        fluid[key][0] = min(fluid[key][0], lo)
        fluid[key][1] = max(fluid[key][1], hi)
    for region in ("fluid", "solid_placeholder"):
        for key in ("T", "p", "rho"):
            for end in ("min", "max"):
                assert 0 <= int(row[region + "_" + key + "_" + end + "_cell"]) < 6070272
    solid = {k: row[k] for k in row if k.startswith("solid_placeholder_")}
    if solid_reference is None:
        solid_reference = solid
    assert solid == solid_reference, "solid regional extrema drift"

conservation = list(read_csv(root / "conservation.csv"))
assert [int(r["step"]) for r in conservation] == expected_steps
for row in conservation:
    # Conservation integration is a run-local diagnostic ledger; statistics
    # retain their independently persisted epoch=7000 across exact restart.
    assert int(row["epoch_start_step"]) == 9500
    assert int(row["normalization_valid"]) == 1
    assert int(row["ibm_adjacent_cells"]) + int(row["interior_cells"]) == 5994352
    assert all(math.isfinite(float(v)) for v in row.values())

trace_summary = []
assert not list(root.glob("cell-trace-rank-*.csv")), "trace must remain disabled"

with (root / "evidence.jsonl").open() as stream:
    seen = []
    for line in stream:
        row = json.loads(line)
        seen.append(row["step"])
        assert row["candidate_identity"]["head"] == "7c03a54a9f1061d1de0e0922ed38a3bcc8f2412f"
        assert row["heap_allocations"] == 0 and not row["retry"]
        assert row["run_start"]["history"]["policy"] == "require_compatible"
    assert seen == expected_steps

# Check actual VTK binary blocks and partition extents one file at a time.
vtk_files = sorted((root / "Visit").rglob("*.vtr"))
assert len(vtk_files) == 128, len(vtk_files)
visit_lists = list((root / "Visit").glob("*.visit"))
assert len(visit_lists) == 1
references = visit_lists[0].read_text().splitlines()
assert references[0] == "!NBLOCKS 128" and len(references) == 129
assert {(visit_lists[0].parent / value).resolve() for value in references[1:]} == {p.resolve() for p in vtk_files}
boxes, total_cells, fields = [], 0, None
visit_digests = {}
for path in vtk_files:
    data = path.read_bytes()
    marker = b'<AppendedData encoding="raw">_'
    pos = data.index(marker)
    header = ET.fromstring(data[:pos] + b'</VTKFile>')
    assert header.attrib["byte_order"] == "LittleEndian" and header.attrib["header_type"] == "UInt64"
    grid = header.find("RectilinearGrid")
    assert grid.attrib["WholeExtent"] == "0 456 0 256 0 52"
    piece = grid.find("Piece")
    box = tuple(map(int, piece.attrib["Extent"].split()))
    n = [box[2 * a + 1] - box[2 * a] for a in range(3)]
    assert all(v > 0 for v in n)
    for a, size in enumerate((456, 256, 52)):
        assert 0 <= box[2 * a] < box[2 * a + 1] <= size
    boxes.append(box)
    cells = n[0] * n[1] * n[2]
    total_cells += cells
    names = [(e.attrib["Name"], int(e.attrib["NumberOfComponents"])) for e in piece.find("CellData")]
    if fields is None:
        fields = names
    assert names == fields
    payload = pos + len(marker)
    end = payload
    arrays = list(piece.find("CellData")) + list(piece.find("Coordinates"))
    for i, entry in enumerate(arrays):
        offset = int(entry.attrib["offset"])
        assert payload + offset == end, "gap/overlap in appended payload"
        length = struct.unpack_from("<Q", data, payload + offset)[0]
        count = cells * int(entry.attrib["NumberOfComponents"]) if i < len(names) else n[i - len(names)] + 1
        assert length == 8 * count
        end += 8 + length
        assert end <= len(data)
        previous = None
        for (value,) in struct.iter_unpack("<d", memoryview(data)[end - length:end]):
            assert math.isfinite(value)
            if i >= len(names):
                assert previous is None or value > previous
                previous = value
    assert data[end:] == b'</AppendedData>\n</VTKFile>\n'
    visit_digests[str(path.relative_to(root))] = hashlib.sha256(data).hexdigest()
assert total_cells == 6070272
for i, a in enumerate(boxes):
    for b in boxes[:i]:
        assert not all(max(a[2*k], b[2*k]) < min(a[2*k+1], b[2*k+1]) for k in range(3)), "partition overlap"

result = dict(schema="HUNDUN_RE3900_FGMRES_PILOT_AUDIT_V1", complete=True,
              step_range=[9501, 9510], steps=10, ranks=128, attempts_per_step=1,
              bdf1_steps=0, bdf2_steps=10, terminal_maxima=metrics, fluid_ranges=fluid,
              solid_extrema=solid_reference, traces=trace_summary, visit_files=128,
              visit_cells=total_cells, visit_fields=fields, visit_sha256=visit_digests,
              checkpoint_rank_files=128, completion_generation=done["restart_generation"],
              inherited_statistics_epoch=stats["statistics_epoch"],
              checkpoint_all_solid_check="not repeated for observation-only change; regional extrema checked",
              statistical_validation=False)
with output.open("x") as stream:
    json.dump(result, stream, indent=2, allow_nan=False)
    stream.write("\n")
print(json.dumps({k: v for k, v in result.items() if k != "visit_sha256"}, indent=2))
