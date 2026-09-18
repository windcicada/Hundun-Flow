#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""AECSC-style mechanical pressure stays separate from fixed EOS pressure."""

import argparse
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile


P0 = 790216.58
STEP = 25000
TIME_S = 0.008621884509921074
DT_S = 2.1487912249540386e-7
CELLS = (8, 8, 8)
FIELDS = 2
SPECIES = ("A", "B")
RU = 8314.46261815324
MW = 28.014
CP_R = 3.3693097541669275
H_A_R = -667.62872778817666
H_B_R = -1004.5597032048694


def write_f64(path: Path, values) -> None:
    with path.open("wb") as stream:
        for value in values:
            stream.write(struct.pack("<d", value))


def enthalpy(temperature: float, y_a: float) -> float:
    return RU / MW * (
        CP_R * temperature + y_a * H_A_R + (1.0 - y_a) * H_B_R
    )


def make_case(fixture: Path, root: Path) -> None:
    shutil.copytree(fixture, root)
    model = json.loads((root / "case.json").read_text(encoding="utf-8"))
    outlet = json.loads(
        (fixture.parent / "reacting-spray-esf-outlet" / "case.json").read_text(
            encoding="utf-8"
        )
    )
    model["flow"]["pressure_reference"] = "boundary_absolute"
    model["flow"]["thermodynamic_pressure_pa"] = P0
    model["boundaries"]["x_min"] = outlet["boundaries"]["x_min"]
    model["boundaries"]["x_max"] = outlet["boundaries"]["x_max"]
    model["boundaries"]["x_max"]["pressure"] = P0
    model["solver"]["coupling"] = "outer_corrected"
    model["time"]["scheme"] = "cn_be"
    model["reaction"]["ensemble"]["fields"] = FIELDS
    model["reaction"]["ensemble"]["initial_species_offsets"] = [0.01, -0.01]
    # Keep the full PDF order A,B while making A the dependent species.  This
    # exercises name-based import instead of the invalid "dependent is last"
    # assumption.
    model["transported_scalars"][0]["stable_name"] = "B"
    for boundary in model["boundaries"].values():
        for scalar in boundary["scalars"]:
            if scalar["stable_name"] == "A":
                scalar["stable_name"] = "B"
                scalar["value"] = 1.0 - scalar["value"]
                scalar["backflow_value"] = 1.0 - scalar["backflow_value"]
    for key in ("initial_dt", "minimum_dt", "maximum_dt"):
        model["time"][key] = DT_S
    (root / "case.json").write_text(
        json.dumps(model, indent=2) + "\n", encoding="utf-8"
    )


def make_spray_case(fixture: Path, root: Path) -> None:
    shutil.copytree(fixture.parent / "reacting-spray-esf-outlet", root)
    model = json.loads((root / "case.json").read_text(encoding="utf-8"))
    model["reaction"]["ensemble"]["fields"] = FIELDS
    model["reaction"]["ensemble"]["initial_species_offsets"] = [0.01, -0.01]
    model["reaction"]["ensemble"]["tcr"] = {"mode": "off"}
    model["turbulence"] = {"model": "vreman"}
    model["spray"].pop("tab_breakup")
    model["spray"]["breakup"] = "stochastic_sgs"
    model["spray"]["relative_tolerance"] = 1.0e-3
    model["spray"]["injectors"][0]["origin_m"] = [0.5625, 0.5625, 0.5625]
    for key in ("initial_dt", "minimum_dt", "maximum_dt"):
        model["time"][key] = DT_S
    (root / "case.json").write_text(
        json.dumps(model, indent=2) + "\n", encoding="utf-8"
    )


def make_spray_transfer(root: Path) -> None:
    root.mkdir()
    diameter = 1.0e-4
    mass = 800.0 * math.pi / 6.0 * diameter ** 3
    parcel = {
        "droplet_mass_kg": mass,
        "source": {
            "kind": "retained",
            "multiplicity": 2.0,
            "position_m": [0.1875, 0.1875, 0.1875],
            "source_ordinal": 17,
            "source_rank": 3,
            "source_sha256": "0123456789abcdef" * 4,
            "velocity_m_per_s": [1.0, -2.0, 3.0],
        },
        "target_thermodynamics": {
            "droplet_diameter_m": diameter,
            "liquid_material_fingerprint": 6004043157121730787,
            "temperature_k": 298.15,
        },
    }
    (root / "parcels.jsonl").write_text(
        json.dumps(parcel, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "format": "hundun_spray_restart_import_v1",
        "parcels": 1,
        "parcels_file": "parcels.jsonl",
        "history_policy": {
            "age": "reset_zero",
            "breakup_ordinal": "reset_zero",
            "sgs": "reset_zero",
            "tab": "reset_zero",
        },
        "injectors": [
            {
                "id": 9007199254740997,
                "next_ordinal": 1234567,
                "residual_mass_kg": 2.5e-11,
            }
        ],
    }
    (root / "restart.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def make_transfer(
    root: Path,
    version: int,
    pressure: float = P0,
    positive_absolute: bool = False,
    base_y_a: float = 0.23,
) -> tuple[float, float]:
    root.mkdir()
    count = math.prod(CELLS)
    flow = []
    rho = []
    pdf = [[] for _ in range(FIELDS)]
    pressure_min = math.inf
    pressure_max = -math.inf
    for z in range(CELLS[2]):
        for y in range(CELLS[1]):
            for x in range(CELLS[0]):
                # The magnitude deliberately matches the real 624CF mechanical
                # field. It is a finite gauge variable, not an EOS pressure.
                mechanical = (
                    pressure + 11.0 * x - 7.0 * y + 3.0 * z
                    if positive_absolute
                    else -17_800_000.0 + 11.0 * x - 7.0 * y + 3.0 * z
                )
                pressure_min = min(pressure_min, mechanical)
                pressure_max = max(pressure_max, mechanical)
                flow.extend((20.0 + x, -2.0 + 0.1 * y, 0.05 * z, mechanical))
                rho.append(pressure * MW / (RU * 800.0))
                for field in range(FIELDS):
                    # A valid pure-species restart at the first local cell
                    # must not be perturbed again by fresh-start ESF offsets.
                    y_a = (
                        1.0
                        if x == 0 and y == 0 and z == 0
                        else base_y_a + 0.01 * field + 0.001 * (x + y + z)
                    )
                    temperature = 760.0 + 20.0 * field + x + 2.0 * y + z
                    pdf[field].extend((y_a, 1.0 - y_a,
                                       enthalpy(temperature, y_a)))
    header = (
        "HUNDUN_PDF_TRANSFER {} {} {} {} {} {:.17g} {:.17g} {:.17g} {} {}\n"
        "{}\n"
    ).format(
        version, *CELLS, STEP, TIME_S, DT_S, pressure, FIELDS, len(SPECIES),
        " ".join(SPECIES)
    )
    (root / "state.txt").write_text(header, encoding="ascii")
    write_f64(root / "flow.f64", flow)
    write_f64(root / "rho_ref.f64", rho)
    (root / "fluid.u8").write_bytes(bytes([1]) * count)
    for field, values in enumerate(pdf):
        write_f64(root / f"pdf{field}.f64", values)
    return pressure_min, pressure_max


def invoke(command, success: bool) -> subprocess.CompletedProcess:
    result = subprocess.run(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        timeout=120, check=False
    )
    if success != (result.returncode == 0):
        raise AssertionError(
            "unexpected import result {}\ncommand: {}\nstdout:\n{}\nstderr:\n{}".format(
                result.returncode, " ".join(map(str, command)),
                result.stdout, result.stderr
            )
        )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--product", required=True, type=Path)
    parser.add_argument("--mpi", required=True)
    parser.add_argument("--ranks", required=True, type=int)
    parser.add_argument("--fixture", required=True, type=Path)
    options = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="hundun-pdf-p0-") as temporary:
        root = Path(temporary)
        case = root / "case"
        spray_case = root / "spray-case"
        transfer = root / "transfer-v2"
        legacy = root / "transfer-v1"
        spray_gas = root / "spray-gas-v1"
        spray_transfer = root / "spray-transfer"
        make_case(options.fixture, case)
        make_spray_case(options.fixture, spray_case)
        pressure_min, pressure_max = make_transfer(transfer, 2)
        make_transfer(legacy, 1)
        make_transfer(
            spray_gas,
            1,
            pressure=101325.0,
            positive_absolute=True,
            base_y_a=0.01,
        )
        make_spray_transfer(spray_transfer)
        prefix = [options.mpi, "-n", str(options.ranks), str(options.binary)]
        result = invoke(prefix + [str(case), str(transfer), str(root / "restart")], True)
        if "pressure=mechanical_perturbation eos=fixed_thermodynamic" not in result.stdout:
            raise AssertionError("missing pressure-separation receipt:\n" + result.stdout)
        audit = json.loads((transfer / "native.json").read_text(encoding="utf-8"))
        assert audit["transfer_format_version"] == 2
        assert audit["pressure_semantics"] == "mechanical_perturbation"
        assert audit["eos_pressure_pa"] == P0
        assert audit["mechanical_pressure_range_pa"] == [pressure_min, pressure_max]
        assert audit["fluid_cells"] == math.prod(CELLS)
        rejected = invoke(
            prefix + [str(case), str(legacy), str(root / "legacy-restart")], False
        )
        if "reconstruct status=1/24106" not in rejected.stdout:
            raise AssertionError(
                "legacy absolute-pressure contract did not fail closed:\n"
                + rejected.stdout + rejected.stderr
            )
        spray_result = invoke(
            prefix
            + [
                str(spray_case),
                str(spray_gas),
                str(root / "spray-restart"),
                str(spray_transfer),
            ],
            True,
        )
        if "history=V5_migration spray_parcels=1" not in spray_result.stdout:
            raise AssertionError("missing spray migration receipt:\n" + spray_result.stdout)
        # Import/readback alone accepts the canonical empty (version-zero) SGS
        # record.  The first stochastic-SGS advance requires an active,
        # version-one record, so exercise that runtime boundary explicitly.
        invoke(
            [
                options.mpi,
                "-n",
                str(options.ranks),
                str(options.product),
                "run",
                str(spray_case),
                "--output",
                str(root / "spray-run"),
                "--steps",
                "1",
                "--restart",
                str(root / "spray-restart"),
                "--output-interval",
                "0",
                "--restart-interval",
                "1",
                "--diagnostics-interval",
                "1",
            ],
            True,
        )
    print(f"pdf fixed-pressure import ranks={options.ranks} ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
