#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Pure-species zero-diffusivity limit through real Cantera and native Restart."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

with tempfile.TemporaryDirectory(prefix='hundun-pasr-pure-') as directory:
    root = Path(directory) / 'case'
    fixture = Path(__file__).resolve().parents[1] / 'fixtures' / 'reacting-cantera-isomer'
    shutil.copytree(str(fixture), str(root))
    model = json.loads((root/'case.json').read_text())
    model['reaction']['model'] = 'auto'
    model['reaction']['ensemble'] = {'fields': 1}
    model['reaction']['mixing'] = {'c_z': 1., 'turbulent_schmidt': .7}
    (root/'case.json').write_text(json.dumps(model))
    run = Path(directory)/'run'
    command = [sys.argv[1], 'run', str(root), '--output', str(run), '--steps', '1',
               '--initial-state', '101325,300,0,0,0,1', '--restart-interval', '1',
               '--output-interval', '0']
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
    assert result.returncode == 0, result.stdout
    result = subprocess.run([sys.argv[1], 'run', str(root), '--output', str(Path(directory)/'next'),
                             '--steps', '1', '--restart', str(run/'Restart'),
                             '--output-interval', '0', '--restart-interval', '0'],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
    assert result.returncode == 0, result.stdout
print('PASS: real Cantera pure species initializes, advances and restores')
