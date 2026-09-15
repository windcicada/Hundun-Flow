# SPDX-License-Identifier: Apache-2.0
"""Exercise mechanism admission and equivalent concentration units via the backend."""
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile

program, fixture = sys.argv[1:]
source = Path(fixture).read_text()
changes = [
    ("model: concentration_c1_v1", "model: unknown"),
    ("negative-floor: 0.1 kmol/m^3", "negative-floor: 0 kmol/m^3"),
    ("negative-floor: 0.1 kmol/m^3", "negative-floor: -1 kmol/m^3"),
    ("fractional-floor: 0.01 kmol/m^3", "fractional-floor: .nan"),
    ("fractional-floor: 0.01 kmol/m^3", "fractional-floor: 1 s"),
    ("negative-floor: 0.1 kmol/m^3", "extra-floor: 0.1 kmol/m^3"),
    ("model: concentration_c1_v1", "model: concentration_c1_v1\n    unused: 1"),
    ("A => B", "A <=> B"),
    ("orders: {A: 0.25, B: -0.75}", "orders: {A: 1.0, B: 0.0}"),
]
with tempfile.TemporaryDirectory(prefix="orders-") as directory:
    path = Path(directory) / "gas.yaml"
    for old, new in changes:
        assert source.count(old) == 1
        text = source.replace(old, new)
        path.write_text(text)
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        subprocess.run([program, "--reject", str(path), sha], check=True)
    text = source.replace("0.1 kmol/m^3", "1e-4 mol/cm^3").replace(
        "0.01 kmol/m^3", "1e-5 mol/cm^3")
    path.write_text(text)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    subprocess.run([program, str(path), sha], check=True)
print("9 invalid inputs rejected; equivalent units passed query and integration")
