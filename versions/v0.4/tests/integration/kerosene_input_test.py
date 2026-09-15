# SPDX-License-Identifier: Apache-2.0
"""Admit the hashed four-step model and its exact chemical identities."""
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile

program, fixture = sys.argv[1:]
source = Path(fixture).read_text()
changes = [
    ("model: kerosene_4_step_v1", "model: other"),
    ("interval: frozen_material", "interval: other"),
    ("interval: frozen_material", "interval: frozen_material, extra: 1"),
    ("model: kerosene_4_step_v1, ", ""),
    (", interval: frozen_material", ""),
    ("composition: {C: 12, H: 23}", "composition: {C: 12, H: 22}"),
    ("reactions: none", "reactions: all"),
]
with tempfile.TemporaryDirectory(prefix="ker-") as directory:
    path = Path(directory) / "gas.yaml"
    for old, new in changes:
        assert source.count(old) == 1, old
        candidate = source.replace(old, new)
        if old == "reactions: none":
            candidate += "\nreactions:\n- equation: 2 H2 + O2 => 2 H2O\n  rate-constant: {A: 1, b: 0, Ea: 0}\n"
        path.write_text(candidate)
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        subprocess.run([program, "--reject", str(path), sha], check=True)
print("seven model and chemical identity admission checks passed")
