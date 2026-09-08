#!/usr/bin/env python3
"""Focused native V1/BE-recovery test for v04_coast_restart_import."""

import argparse
import json
import math
import pathlib
import shutil
import struct
import subprocess
import tempfile


AIR_O2 = 0.23291751145757963
COAST_RU = 8314.3
MOLECULAR_WEIGHTS = (16.04308, 31.9988, 28.0134)


def write_f32(path: pathlib.Path, values) -> None:
    with path.open("wb") as stream:
        for value in values:
            stream.write(struct.pack("<f", value))


def scalar_boundary(name: str, value: float, inlet: bool) -> dict:
    return {
        "stable_name": name,
        "kind": "dirichlet" if inlet else "zero_gradient",
        "value": value,
        "backflow_kind": "dirichlet",
        "backflow_value": value,
    }


def make_case(template: pathlib.Path, thermo: pathlib.Path,
              case_root: pathlib.Path) -> None:
    model = json.loads(template.read_text(encoding="utf-8"))
    shutil.copyfile(thermo, case_root / "thermophysics.d")
    model["transported_scalars"] = [
        {
            "stable_name": "CH4",
            "role": "species",
            "molecular_schmidt": 0.7,
            "turbulent_schmidt": 0.7,
        },
        {
            "stable_name": "O2",
            "role": "species",
            "molecular_schmidt": 0.7,
            "turbulent_schmidt": 0.7,
        },
    ]
    for face, boundary in model["boundaries"].items():
        inlet = face == "x_min"
        boundary["scalars"] = [
            scalar_boundary("CH4", 0.0, inlet),
            scalar_boundary("O2", AIR_O2, inlet),
        ]
    (case_root / "case.json").write_text(
        json.dumps(model, indent=2) + "\n", encoding="utf-8"
    )


def make_transfer(root: pathlib.Path, invalid_y: bool = False) -> None:
    shape = (8, 8, 8)
    fields = {name: [] for name in (
        "u", "v", "w", "p_abs", "temperature", "Y_CH4", "rho", "flow_h"
    )}
    for z in range(shape[2]):
        for y in range(shape[1]):
            for x in range(shape[0]):
                methane = (x + 2.0 * y + 3.0 * z) / 42.0
                if invalid_y and x == y == z == 0:
                    methane = 1.25
                pressure = 100000.0 + 3.0 * x - 2.0 * y + z
                temperature = 295.0 + 0.2 * x + 0.1 * y + 0.05 * z
                oxygen = (1.0 - methane) * AIR_O2
                nitrogen = 1.0 - methane - oxygen
                gas_constant = COAST_RU * (
                    methane / MOLECULAR_WEIGHTS[0]
                    + oxygen / MOLECULAR_WEIGHTS[1]
                    + nitrogen / MOLECULAR_WEIGHTS[2]
                )
                fields["u"].append(0.1 + 0.01 * x)
                fields["v"].append(-0.02 + 0.002 * y)
                fields["w"].append(0.003 * z)
                fields["p_abs"].append(pressure)
                fields["temperature"].append(temperature)
                fields["Y_CH4"].append(methane)
                fields["rho"].append(pressure / (gas_constant * temperature))
                fields["flow_h"].append(-3185.0 - 100000.0 * methane)
    for name, values in fields.items():
        write_f32(root / f"{name}.f32", values)
    (root / "fluid_mask.u8").write_bytes(bytes([1]) * math.prod(shape))


def run_import(arguments, expect_success: bool) -> subprocess.CompletedProcess:
    result = subprocess.run(
        arguments, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False, timeout=120
    )
    if expect_success and result.returncode != 0:
        raise AssertionError(
            f"import failed ({result.returncode})\nstdout:\n{result.stdout}"
            f"\nstderr:\n{result.stderr}"
        )
    if not expect_success and result.returncode == 0:
        raise AssertionError("invalid species input unexpectedly succeeded")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=pathlib.Path)
    parser.add_argument("--mpi", required=True)
    parser.add_argument("--ranks", required=True, type=int)
    parser.add_argument("--case-template", required=True, type=pathlib.Path)
    parser.add_argument("--thermo", required=True, type=pathlib.Path)
    options = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="hundun-coast-import-") as temporary:
        root = pathlib.Path(temporary)
        case = root / "case"
        transfer = root / "transfer"
        restart = root / "restart"
        invalid = root / "invalid-transfer"
        invalid_restart = root / "invalid-restart"
        case.mkdir()
        transfer.mkdir()
        invalid.mkdir()
        make_case(options.case_template, options.thermo, case)
        make_transfer(transfer)
        make_transfer(invalid, invalid_y=True)
        prefix = [options.mpi, "-n", str(options.ranks), str(options.binary)]
        result = run_import(
            prefix + ["--case", str(case), "--transfer", str(transfer),
                      "--output", str(restart)], True
        )
        if "source_format=1 backward_euler_recovery=true" not in result.stdout:
            raise AssertionError(f"missing V1 recovery receipt: {result.stdout}")
        receipt = json.loads(
            (restart / "mean-field-transfer.json").read_text(encoding="utf-8")
        )
        assert receipt["status"] == "ok"
        assert receipt["restart_format_version"] == 1
        assert receipt["backward_euler_recovery_verified"] is True
        assert receipt["cell_counts"] == {"fluid": 512, "solid_placeholder": 0}
        assert receipt["native_readback_max_absolute"] == {
            "fields": 0, "mass_flux": 0
        }
        assert receipt["native_flux_adjustment"] == {
            "changed_owned_faces": 0,
            "sum_absolute_change_kg_s": 0,
            "internal_zero_to_nonzero_source_faces": 0,
            "internal_source_absolute_mass_flow_kg_s": 0,
            "configured_immersed_source_mass_flow_kg_s": 0,
        }
        memory = receipt["memory_estimate_max_per_rank_bytes"]
        assert memory["bridge_payload"] > 0
        assert memory["sealed_arena"] > 0
        assert memory["single_product_sum_upper_estimate"] >= (
            memory["bridge_payload"] + memory["sealed_arena"]
        )
        assert receipt["gas_constants_J_per_kmol_K"]["COAST"] == COAST_RU
        invalid_result = run_import(
            prefix + ["--case", str(case), "--transfer", str(invalid),
                      "--output", str(invalid_restart)], False
        )
        if "transfer_status=5/24004" not in invalid_result.stderr:
            raise AssertionError(
                f"unexpected fail-closed diagnostic: {invalid_result.stderr}"
            )
        if (invalid_restart / "current").exists():
            raise AssertionError("invalid input published a restart")
    print(f"coast_restart_import_cli ranks={options.ranks} ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
