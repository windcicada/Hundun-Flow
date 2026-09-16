#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Native THICK_EX with real kerosene thermo and 1/2/4-rank restart."""
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET



def validate_esf_records(records):
    for label, rows in records.items():
        for i, row in enumerate(rows):
            p = row["payload"]
            budget = p["composition_balance"]
            assert budget["after_parcel_exchange"]
            assert budget["density"] == "field0"
            assert budget["composition"] == "physical_ensemble_mean"
            assert len(budget["species"]) == 7 and len(budget["elements"]) == 4
            for group in ("species", "elements"):
                for entry in budget[group]:
                    assert all(math.isfinite(v) for k, v in entry.items() if k != "name")
                    assert entry["relative_defect"] < 1e-6, (label, row["step"], entry)
            previous = (rows[i - 1] if i else records["1"][-1] if label in ("2", "r")
                        else records["2"][-1] if label == "4" else None)
            if previous:
                inventory = {e["name"]: e["current_inventory"] for e in
                             previous["payload"]["composition_balance"]["species"]}
                # THICK_EX exchanges pure C12H23. The starting gas-stage
                # inventory includes this same parcel source exactly once.
                for entry in budget["species"]:
                    delta = p["phase_mass_input_kg_s"] * p["dt"] if entry["name"] == "C12H23" else 0.
                    error = entry["accepted_inventory"] - inventory[entry["name"]] - delta
                    assert abs(error) < 1e-12 * max(1., abs(entry["accepted_inventory"]))
    for group in ("species", "elements"):
        a, b = (records[label][-1]["payload"]["composition_balance"][group] for label in ("2", "r"))
        for x, y in zip(a, b):
            assert x["name"] == y["name"]
            for key in ("accepted_inventory", "current_inventory"):
                assert abs(x[key] - y[key]) < 1e-11 * max(1., abs(x[key]), abs(y[key]))


def validate_sgs_output(work, case):
    """Read actual parallel output and compare the same-step derived fields."""
    names = ("nu_sgs", "k_sgs", "eps_sgs_volume", "eps_sgs_specific")
    shape = case["mesh"]["exact_cells"]
    lower, upper = (case["mesh"]["domain"][key] for key in ("lower", "upper"))
    spacing = [(hi-lo)/n for lo, hi, n in zip(lower, upper, shape)]
    vertices = [tuple(map(float, line.split()[1:]))
                for line in (work / "cube.stl").read_text().splitlines()
                if line.strip().startswith("vertex ")]
    assert len(set(vertices)) == 8
    box = [(min(v[a] for v in vertices), max(v[a] for v in vertices)) for a in range(3)]
    frames, files = {}, {}
    solid_samples, active_samples = 0, 0
    maxima = dict.fromkeys(names, 0.)
    for label, count in (("1", 3), ("2", 2), ("r", 1), ("4", 4)):
        paths = sorted((work.with_name(work.name + label) / "Visit").glob("*.vti"))
        assert len(paths) == count, (label, len(paths))
        for path in paths:
            step = int(path.name.split("-")[1])
            cells = frames.setdefault((label, step), {})
            data = path.read_bytes()
            split = data.index(b"<AppendedData")
            tree = ET.fromstring(data[:split] + b"</VTKFile>")
            piece = tree.find("./ImageData/Piece")
            extent = list(map(int, piece.attrib["Extent"].split()))
            begin = extent[::2]
            counts = [extent[a+1]-extent[a] for a in (0, 2, 4)]
            n = counts[0]*counts[1]*counts[2]
            base = data.index(b"_", data.index(b">", split)) + 1
            fields = {}
            for node in piece.find("CellData"):
                name = node.attrib["Name"]
                if name not in names:
                    continue
                offset = base + int(node.attrib["offset"])
                length, = struct.unpack_from("<Q", data, offset)
                assert node.attrib["type"] == "Float64"
                assert int(node.attrib["NumberOfComponents"]) == 1 and length == n*8
                fields[name] = struct.unpack_from("<{}d".format(n), data, offset+8)
            assert set(fields) == set(names)
            for cell in range(n):
                local = (cell % counts[0], (cell // counts[0]) % counts[1], cell // (counts[0]*counts[1]))
                global_cell = tuple(b+i for b, i in zip(begin, local))
                assert all(0 <= i < limit for i, limit in zip(global_cell, shape))
                assert global_cell not in cells
                center = [lo+(i+.5)*dx for lo, i, dx in zip(lower, global_cell, spacing)]
                solid = all(lo < x < hi for x, (lo, hi) in zip(center, box))
                values = tuple(fields[name][cell] for name in names)
                assert all(math.isfinite(v) and v >= 0 for v in values)
                if solid:
                    assert values == (0., 0., 0., 0.)
                    solid_samples += 1
                else:
                    active_samples += values[0] > 0
                for name, value in zip(names, values):
                    maxima[name] = max(maxima[name], value)
                cells[global_cell] = values
            files[str(path)] = hashlib.sha256(data).hexdigest()
    assert set(frames) == {("1", 1), ("1", 2), ("1", 3), ("2", 4), ("r", 4), ("4", 5)}
    assert all(len(cells) == shape[0]*shape[1]*shape[2] for cells in frames.values())
    assert solid_samples > 0 and active_samples > 0 and all(v > 0 for v in maxima.values())
    a, b = frames[("2", 4)], frames[("r", 4)]
    errors = {name: max(abs(a[cell][i]-b[cell][i]) for cell in a) /
              max(1., max(max(abs(a[cell][i]), abs(b[cell][i])) for cell in a))
              for i, name in enumerate(names)}
    assert max(errors.values()) < 1e-11, errors
    return dict(model="vreman", frames=len(frames), files=len(files),
                active_samples=active_samples, solid_samples=solid_samples,
                maximum=maxima, same_step_normalized_difference=errors, sha256=files)


def main():
    binary, fixture, mpi, work, validator = map(Path, sys.argv[1:6])
    mode = sys.argv[6] if len(sys.argv) > 6 else "backward_euler"
    assert mode in ("backward_euler", "cn_be", "cn_be_esf", "cn_be_esf_ibm", "cn_be_esf_ibm_sgs")
    sgs = mode == "cn_be_esf_ibm_sgs"
    ibm = mode in ("cn_be_esf_ibm", "cn_be_esf_ibm_sgs")
    esf = mode == "cn_be_esf" or ibm
    scheme = "cn_be" if esf else mode
    method = "CN/BE+ESF" if esf else "CN/BE" if scheme == "cn_be" else "BE/PISO"
    if ibm:
        method += "+IBM"
    if sgs:
        method += "+Vreman"
    work.mkdir(parents=True, exist_ok=True)
    case = json.loads((fixture / "case.json").read_text())
    if sgs:
        case["turbulence"] = dict(model="vreman")
    for name in ("gas.yaml", "thermophysics.d", "viscosity.dat"):
        shutil.copyfile(str(fixture / name), str(work / name))
    thermo = (fixture / "thermophysics.d").read_text()
    species = re.findall(r"^species (\S+)$", thermo, re.M)
    fuel = thermo.split("species C12H23\n")[1].split("end_species")[0]
    mw = float(re.search(r"molecular_weight (\S+)", fuel).group(1))
    cp = [float(x) for x in re.search(r"nasa7_low (.*)", fuel).group(1).split()]
    temperature = 298.15
    vapor_h = 8314.46261815324 / mw * (
        sum(cp[i] * temperature ** (i + 1) / (i + 1) for i in range(5)) + cp[5])
    latent = 250183 * ((684.26 - temperature) / (684.26 - 483.15)) ** .38
    asset = """HUNDUN_LIQUID_ASSET_V1
units SI
source Rachner-kerosene-phase-transport-fixture
gas_sha {sha}
gas_phase {phase}
gas_species {count} {species}
enthalpy_reference absolute-standard-formation-298.15K-v1
vapor C12H23
molecular_weight {mw:.17g}
temperature_range 200 680
liquid_reference 298.15 {href:.17g}
density kerosene_density_v1 298.15 0 0 0 0
cp kerosene_cp_v1 298.15 0 0 0 0
latent kerosene_latent_v1 298.15 0 0 0 0
surface_tension constant 298.15 0.025 0 0 0
viscosity constant 298.15 0.0008 0 0 0
saturation kerosene_v1
end
""".format(sha=case["reaction"]["mechanism_sha256"],
           phase=case["reaction"]["phase"], count=len(species),
           species=" ".join(species), mw=mw, href=vapor_h - latent)
    (work / "liquid.asset").write_text(asset)
    fingerprint = 14695981039346656037
    for byte in asset.encode():
        fingerprint = ((fingerprint ^ byte) * 1099511628211) & ((1 << 64) - 1)
    case["spray"] = dict(
        evaporation="thick_exchange", liquid_file="liquid.asset",
        liquid_fingerprint=fingerprint, seed=42, maximum_local_parcels=32,
        maximum_local_segments=8192, maximum_substep_s=1e-9,
        minimum_substep_s=1e-14, relative_tolerance=1e-6, tab_breakup=False,
        injectors=[dict(id=1, origin_m=[.49, .49, .49], axis=[1, 0, 0],
                        cone_half_angle_rad=0., speed_m_per_s=2.,
                        mass_flow_rate_kg_per_s=1., represented_mass_per_parcel_kg=1e-9,
                        droplet_diameter_m=1e-4, temperature_k=350.)])
    if scheme == "cn_be":
        # A resolved synthetic exchange makes omission of gas kinetic energy
        # visible to the total-inventory balance, while retaining one parcel
        # per step and the same thermophysical assets.
        case["spray"]["injectors"][0].update(
            speed_m_per_s=200., mass_flow_rate_kg_per_s=1e5,
            represented_mass_per_parcel_kg=1e-4)
    if esf:
        independent = len(species) - 1
        offsets = [0.] * (2 * independent)
        offsets[0], offsets[independent] = .0001, -.0001
        case["reaction"]["model"] = "esf_tpdf"
        case["reaction"]["ensemble"] = dict(
            fields=2, seed=1234, initial_species_offsets=offsets, tcr=dict(mode="off"))
        case["reaction"]["mixing"] = dict(c_z=.25, turbulent_schmidt=.7)
    if ibm:
        wall_case = fixture.with_name("reacting-spray-esf-ibm")
        shutil.copyfile(str(wall_case / "cube.stl"), str(work / "cube.stl"))
        wall_model = json.loads((wall_case / "case.json").read_text())
        case["mesh"]["immersed_boundary"] = wall_model["mesh"]["immersed_boundary"]
        case["mesh"]["limits"]["max_memory_bytes_per_rank"] = wall_model["mesh"]["limits"]["max_memory_bytes_per_rank"]
        case["mesh"]["exact_cells"] = wall_model["mesh"]["exact_cells"]
        case["mesh"]["minimum_spacing"] = wall_model["mesh"]["minimum_spacing"]
        # The first interval reaches the cube wall. Deposition must select
        # the fluid side while the common source and field transport retain
        # their IBM matrix/flux constraints.
        case["spray"]["injectors"][0]["origin_m"] = [.375 - 1e-7, .5, .5]
    case["solver"]["coupling"] = "outer_corrected" if scheme == "cn_be" else "PISO"
    case["solver"].pop("cold_stopping", None)
    case["time"]["scheme"] = scheme
    for key in ("initial_dt", "minimum_dt", "maximum_dt"):
        case["time"][key] = 1e-9
    (work / "case.json").write_text(json.dumps(case, indent=2) + "\n")

    def command(args, path):
        with path.open("w") as log:
            result = subprocess.run([str(x) for x in args], stdout=log, stderr=log)
        if result.returncode:
            raise RuntimeError(path.read_text())

    command([binary, "check", work], work / "check.log")
    check = (work / "check.log").read_text()
    assert "evaporation=thick_exchange" in check
    reaction_model = "kerosene_4_step_v1/frozen_material"
    assert "reaction_model=" + reaction_model in check
    if sgs:
        assert "sgs=vreman" in check
    records = {}

    def run(label, ranks, count, restart=None):
        output = work.with_name(work.name + label)
        # The CLI appends evidence for real continuation. Each independent
        # test invocation owns fresh output; restart inputs are other labels.
        if output.exists():
            shutil.rmtree(str(output))
        args = [mpi, "--oversubscribe", "--bind-to", "none", "-n", str(ranks),
                binary, "run", work, "--output", output, "--steps", str(count),
                "--max-dt", "1e-9", "--output-interval", "1" if sgs else "0",
                "--restart-interval", "1", "--diagnostics-interval", "1"]
        if restart:
            args += ["--restart", restart / "Restart"]
        else:
            args += ["--initial-state", "100000,900,0,0,0,0.001,0.05,0.005,0.01,0.2,0.02"]
        command(args, work / (label + ".log"))
        validation = [sys.executable, validator, "runtime", output / "evidence.jsonl"]
        if restart:
            generation = (restart / "Restart/current").read_text().strip()
            validation += ["--run-start-manifest", restart / "Restart" / generation / "manifest.bin"]
        command(validation, work / (label + "-validate.log"))
        evidence = [json.loads(line) for line in (output / "evidence.jsonl").read_text().splitlines()]
        rows = [json.loads(line) for line in (output / "diagnostics.jsonl").read_text().splitlines()]
        assert len(rows) == count and all(not row["retry"] for row in evidence)
        for row in rows:
            p = row["payload"]
            assert p["phase_mass_input_kg_s"] > 0
            assert abs(p["mass_balance_defect_kg_s"] * p["dt"]) / p["mass_kg"] < 1e-12
            scale = max(1., abs(p["internal_energy_J"]) + p["kinetic_energy_J"])
            assert abs(p["total_energy_balance_defect_W"] * p["dt"]) / scale < 1e-12
            if scheme == "cn_be":
                if not ibm:
                    assert p["kinetic_energy_J"] > 1e-11 * scale
                assert abs(p["total_energy_balance_defect_W"]) / max(
                    1., abs(p["phase_energy_input_W"])) < 1e-6
        records[label] = rows
        return output

    initial = run("1", 1, 3)
    two = run("2", 2, 1, initial)
    run("r", 1, 1, initial)
    run("4", 4, 1, two)
    for key in ("mass_kg", "internal_energy_J", "kinetic_energy_J",
                "phase_mass_input_kg_s", "phase_energy_input_W"):
        a, b = (records[label][-1]["payload"][key] for label in ("2", "r"))
        assert abs(a-b) < 1e-10 * max(1., abs(a), abs(b)), (key, a, b)
    if esf:
        validate_esf_records(records)
        if len(sys.argv) > 7:
            command([mpi, "--oversubscribe", "--bind-to", "none", "-n", "1", sys.argv[7],
                     work, work.with_name(work.name + "r") / "Restart", two / "Restart"]
                    + (["--wall"] if ibm else []), work / "compare.log")
    report = dict(scope=method + " THICK_EX spray with native kerosene four-step chemistry",
                  reaction_model=reaction_model,
                  binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  ranks=[1, 2, 4], restart_steps=[4, 5], passed=True,
                  records=records)
    if sgs:
        report["sgs"] = validate_sgs_output(work, case)
    (work / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print("native_thick thermo=kerosene time=" + method + " phase_balance=pass restart=1/2/4 pass")


if __name__ == "__main__":
    main()
