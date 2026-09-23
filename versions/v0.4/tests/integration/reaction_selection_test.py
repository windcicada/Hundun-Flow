#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exercise canonical reaction selection through the production case compiler."""
import copy
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

binary = sys.argv[1]
fixtures = Path(__file__).resolve().parents[1] / 'fixtures'


def validate(root, spec, accepted=True):
    (root / 'case.json').write_text(json.dumps(spec))
    result = subprocess.run([binary, 'validate', str(root), '--dry-plan'],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
    assert (result.returncode == 0) == accepted, result.stdout
    if not accepted:
        return None
    match = re.search(r'VALID case=(\d+) product=(\d+)', result.stdout)
    assert match, result.stdout
    return match.groups(), result.stdout


with tempfile.TemporaryDirectory(prefix='hundun-reaction-selection-') as directory:
    root = Path(directory) / 'case'
    shutil.copytree(str(fixtures / 'reacting-pasr'), str(root))
    base = json.loads((root / 'case.json').read_text())
    expected, _ = validate(root, base)
    count = 1
    for model in ('auto', 'esf_tpdf', None):
        for ensemble in (None, {}, {'fields': 1},
                         {'fields': 1, 'seed': 42, 'initial_species_offsets': [],
                          'tcr': {'mode': 'off'}}):
            spec = copy.deepcopy(base)
            if model is None:
                del spec['reaction']['model']
            else:
                spec['reaction']['model'] = model
            if ensemble is not None:
                spec['reaction']['ensemble'] = ensemble
            identity, output = validate(root, spec)
            assert identity == expected and 'chemistry=frozen_pasr ' in output, output
            count += 1
    defaults = copy.deepcopy(base)
    defaults['reaction']['mixing'] = {'c_z': 1., 'turbulent_schmidt': .7}
    expected_defaults, _ = validate(root, defaults)
    defaults['reaction'].pop('model')
    defaults['reaction'].pop('mixing')
    identity, output = validate(root, defaults)
    assert identity == expected_defaults, output
    count += 2
    auto = copy.deepcopy(base)
    auto['reaction']['model'] = 'auto'
    auto['reaction']['ensemble'] = {'fields': 1}
    # The normalized plan must also support native continuation across spellings.
    (root / 'case.json').write_text(json.dumps(auto))
    first = Path(directory) / 'first'
    command = [binary, 'run', str(root), '--output', str(first), '--steps', '1',
               '--initial-state', '101325,300,0,0,0,0.2', '--output-interval', '0',
               '--restart-interval', '1']
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
    assert result.returncode == 0, result.stdout
    (root / 'case.json').write_text(json.dumps(base))
    result = subprocess.run([binary, 'run', str(root), '--output', str(Path(directory)/'next'),
                             '--steps', '1', '--restart', str(first/'Restart'),
                             '--output-interval', '0', '--restart-interval', '0'],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
    assert result.returncode == 0, result.stdout
    for fields in (0, 3, 65, 1.5, -1):
        spec = copy.deepcopy(auto)
        spec['reaction']['ensemble']['fields'] = fields
        validate(root, spec, False)
        count += 1
    for ensemble in ({'fields': 1, 'initial_species_offsets': [0.]},
                     {'fields': 1, 'tcr': {'mode': 'experimental'}},
                     {'fields': 1, 'typo': 2}):
        spec = copy.deepcopy(auto)
        spec['reaction']['ensemble'] = ensemble
        validate(root, spec, False)
        count += 1
    # Even-field dispatch preserves the established ESF model identity.
    for fields in (2, 4):
        spec = copy.deepcopy(base)
        spec['reaction']['model'] = 'esf_tpdf'
        spec['reaction']['ensemble'] = {'fields': fields, 'seed': 0,
                                       'initial_species_offsets': [], 'tcr': {'mode':'off'}}
        expected_esf, _ = validate(root, spec)
        spec['reaction']['model'] = 'auto'
        spec['reaction']['ensemble'] = {'fields': fields}
        identity, output = validate(root, spec)
        assert identity == expected_esf, output
        count += 2
    print('PASS: {} selection cases and PaSR native continuation'.format(count))
