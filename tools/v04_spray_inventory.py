#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Decode the audited COAST FP32 spray restart into an SI transfer inventory.

NPOINT/INDEX identify retained slots; INDEX_OUT identifies retired slots.
PART diameter and daughter diameter are micrometres. Daughter cell labels
are stored as REAL, despite referring to integer cells. This is a decoded
transfer asset: target IBM ownership, model histories and native Restart
construction follow the selected migration policy.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct


def require(condition, message):
    if not condition:
        raise ValueError(message)


def records(data):
    result, offset = [], 0
    while offset < len(data):
        require(offset+8 <= len(data), "truncated record marker")
        size, = struct.unpack_from("<i", data, offset)
        require(size >= 0 and offset+size+8 <= len(data), "invalid record length")
        require(struct.unpack_from("<i", data, offset+size+4)[0] == size,
                "record marker mismatch")
        result.append(memoryview(data)[offset+4:offset+size+4])
        offset += size+8
    require(len(result) >= 14, "missing spray records")
    return result


def decode(data, shape, rank):
    require(len(shape) == 3 and all(n >= 3 for n in shape), "shape includes two ghost layers")
    ncell = shape[0]*shape[1]*shape[2]
    r = records(data)
    require(all(len(r[i]) == 4*ncell for i in range(5)), "source field shape")
    require(len(r[5]) == 8 and len(r[6]) == 16 and len(r[7]) == 4, "spray header size")
    npar, step = struct.unpack("<2i", r[5])
    nsum, added, npout, children = struct.unpack("<4i", r[6])
    excess, = struct.unpack("<f", r[7])
    require(min(npar, step, nsum, added, npout, children) >= 0 and npout <= npar,
            "negative counter or retired count")
    require(math.isfinite(excess), "nonfinite injection excess")
    require(len(r) == 14+children, "record count")
    require(all(len(r[8+i]) == 40 for i in range(children)), "daughter record size")
    lengths = ((ncell+2)*4, npar*22*4, npar*4, npar*4, npar*3*4, npout*4)
    require(all(len(r[8+children+i]) == length for i, length in enumerate(lengths)),
            "particle record layout")
    point = struct.unpack("<{}i".format(ncell+2), r[8+children])
    index = struct.unpack("<{}i".format(npar), r[10+children])
    cell_index = struct.unpack("<{}i".format(npar), r[11+children])
    retired = struct.unpack("<{}i".format(npout), r[13+children])
    retained = npar-npout
    require(point[0] == 1 and point[-1] == retained+1 and
            all(a <= b for a, b in zip(point, point[1:])), "NPOINT range")
    seen = bytearray(npar)
    for slot in retired:
        require(1 <= slot <= npar and not seen[slot-1], "duplicate/invalid retired slot")
        seen[slot-1] = 1
    for slot in index[:retained]:
        require(1 <= slot <= npar and not seen[slot-1], "duplicate/retired retained slot")
        seen[slot-1] = 1
    require(all(seen), "slot partition coverage")
    digest = hashlib.sha256(data).hexdigest()
    counts = dict(slots=npar, retired=npout, indexed=retained, positive=0,
                  zero_inventory=0, interior_positive=0, halo_positive=0,
                  dummy_positive=0, pending_children=children, positive_children=0,
                  zero_inventory_children=0)
    parcels = []

    def append(kind, ordinal, xyz, uvw, diameter, number, temperature, cell, history, bin_kind):
        require(all(math.isfinite(v) for v in (*xyz, *uvw, diameter, number, temperature, *history)),
                "nonfinite parcel state")
        require(diameter >= 0 and number >= 0, "negative parcel diameter or weight")
        if diameter == 0 or number == 0:
            counts["zero_inventory_children" if kind == "daughter" else "zero_inventory"] += 1
            return
        require(temperature > 0, "nonpositive parcel temperature")
        counts["positive_children" if kind == "daughter" else "positive"] += 1
        if kind == "retained":
            counts[bin_kind+"_positive"] += 1
        parcels.append(dict(source_sha256=digest, source_rank=rank, source_ordinal=ordinal,
                            kind=kind, coast_cell=cell, source_bin=bin_kind,
                            position_m=xyz, velocity_m_per_s=uvw,
                            droplet_diameter_m=diameter*1e-6, multiplicity=number,
                            temperature_k=temperature, legacy_state=list(history)))

    for cell in range(ncell+1) if retained else ():
        i, j, k = cell % shape[0], (cell//shape[0]) % shape[1], cell//(shape[0]*shape[1])
        bin_kind = "dummy" if cell == ncell else "interior" if (
            0 < i < shape[0]-1 and 0 < j < shape[1]-1 and 0 < k < shape[2]-1) else "halo"
        for entry in range(point[cell]-1, point[cell+1]-1):
            slot = index[entry]
            require(cell_index[slot-1] == cell+1, "NPOINT/INDEX_CELL disagreement")
            part = struct.unpack_from("<22f", r[9+children], (slot-1)*88)
            xyz = struct.unpack_from("<3f", r[12+children], (slot-1)*12)
            append("retained", slot, xyz, part[:3], part[3], part[4], part[5], cell+1, part, bin_kind)
    for i in range(children):
        child = struct.unpack("<10f", r[8+i])
        require(math.isfinite(child[9]) and child[9].is_integer() and 1 <= child[9] <= ncell,
                "daughter REAL cell label")
        append("daughter", i+1, child[:3], child[3:6], child[6], child[7], child[8],
               int(child[9]), child, "pending")
    report = dict(rank=rank, sha256=digest, bytes=len(data), spray_step=step,
                  cumulative_injected=nsum, cumulative_added=added, excess_mass_kg=excess,
                  counts=counts)
    return report, parcels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("restart", type=Path)
    parser.add_argument("--shape", type=int, nargs=3, required=True, metavar=("NX", "NY", "NZ"))
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(args.ranks > 0, "positive rank count")
    paths = sorted(args.restart.glob("Restart_spray_PDF.[0-9][0-9][0-9]"))
    require([int(p.suffix[1:]) for p in paths] == list(range(args.ranks)), "complete rank set")
    args.output.mkdir(parents=True, exist_ok=True)
    destination = args.output / "parcels.jsonl"
    require(not destination.exists() and not (args.output/"inventory.json").exists(), "fresh transfer output")
    temporary = args.output / "parcels.tmp"
    rows, volumes, number = [], [], []
    temporary_created = False
    try:
        with temporary.open("x") as stream:
            temporary_created = True
            for rank, path in enumerate(paths):
                row, parcels = decode(path.read_bytes(), args.shape, rank)
                rows.append(row)
                for parcel in parcels:
                    stream.write(json.dumps(parcel, separators=(",", ":"), allow_nan=False)+"\n")
                    volumes.append(math.pi/6*parcel["droplet_diameter_m"]**3*parcel["multiplicity"])
                    number.append(parcel["multiplicity"])
        require(len({row["spray_step"] for row in rows}) == 1, "consistent spray step")
        counts = {key: sum(row["counts"][key] for row in rows) for key in rows[0]["counts"]}
        report = dict(schema="hundun_coast_spray_transfer_v1", shape_with_ghosts=args.shape,
                      spray_step=rows[0]["spray_step"], counts=counts,
                      represented_droplet_number=math.fsum(number), liquid_volume_m3=math.fsum(volumes),
                      fields="SI primitive state; retained PART[1:22] and pending daughter REAL[1:10]",
                      scope="decoded transfer inventory; target ownership and native model restoration follow",
                      transfer_sha256=hashlib.sha256(temporary.read_bytes()).hexdigest(), ranks=rows)
        temporary.rename(destination)
        (args.output/"inventory.json").write_text(json.dumps(report, indent=2)+"\n")
    except BaseException:
        if temporary_created and temporary.exists():
            temporary.unlink()
        raise
    print(json.dumps({key: report[key] for key in ("schema", "spray_step", "counts", "liquid_volume_m3")}, indent=2))


if __name__ == "__main__":
    main()
